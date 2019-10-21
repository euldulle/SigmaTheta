/*   sigma_theta.c                                F. Vernotte - 2010/10/25  */
/*                                    Modified by FV for MDEV - 2014/10/01  */
/*       Modified for selecting the deviation by argument, FV - 2015/06/25  */
/*                                        */
/*   Computation of the (modified) Allan Deviation measurement of a         */
/*   frequency deviation measurement serie                                  */
/*                                                                          */
/*                                                   - SIGMA-THETA Project  */
/*                                                                          */
/* Copyright or © or Copr. Université de Franche-Comté, Besançon, France    */
/* Contributor: François Vernotte, UTINAM/OSU THETA (2012/07/17)            */
/* Contact: francois.vernotte@obs-besancon.fr                               */
/*                                                                          */
/* This software, SigmaTheta, is a collection of computer programs for      */
/* time and frequency metrology.                                            */
/*                                                                          */
/* This software is governed by the CeCILL license under French law and     */
/* abiding by the rules of distribution of free software.  You can  use,    */
/* modify and/ or redistribute the software under the terms of the CeCILL   */
/* license as circulated by CEA, CNRS and INRIA at the following URL        */
/* "http://www.cecill.info".                                                */
/*                                                                          */
/* As a counterpart to the access to the source code and  rights to copy,   */
/* modify and redistribute granted by the license, users are provided only  */
/* with a limited warranty  and the software's author,  the holder of the   */
/* economic rights,  and the successive licensors  have only  limited       */
/* liability.                                                               */
/*                                                                          */
/* In this respect, the user's attention is drawn to the risks associated   */
/* with loading,  using,  modifying and/or developing or reproducing the    */
/* software by the user in light of its specific status of free software,   */
/* that may mean  that it is complicated to manipulate,  and  that  also    */
/* therefore means  that it is reserved for developers  and  experienced    */
/* professionals having in-depth computer knowledge. Users are therefore    */
/* encouraged to load and test the software's suitability as regards their  */
/* requirements in conditions enabling the security of their systems and/or */
/* data to be ensured and,  more generally, to use and operate it in the    */
/* same conditions as regards security.                                     */
/*                                                                          */
/* The fact that you are presently reading this means that you have had     */
/* knowledge of the CeCILL license and that you accept its terms.           */
/*                                                                          */
/*                                                                          */
                                    
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include "sigma_theta.h"

#define db(x) ((double)(x))
#define sisig(x) ( (x) == 0 ) ? (db(0)) : (  ( (x) > 0 ) ? (db(1)) : (db(-1))  )

void usage(void)
/* Help message */
    {
    printf("##############################################################################################################\n\n");
    printf(" SigmaTheta : a tool from the SigmaTheta suite\n\n");
    printf("    Usage : SigmaTheta [-amHph] [-x xscalingfactor] [SOURCE [TARGET]]\n\n");
    printf("     Computes the (modified) Allan Deviation of a sequence of normalized frequency deviation measurements.\n\n");
    printf("     If TARGET file is not specified, builds target name out of SOURCE name\n");
    printf("     If SOURCE file is not specified, reads input from stdin and builds up a generic name for TARGET \n");
    printf("     output goes to both stdout and the TARGET file\n\n");
    printf("           -a use the standard Allan variance (default). \n");
    printf("           -m use the modified Allan variance. \n");
    printf("           -H use the Hadamard variance. \n");
    printf("           -p use the Parabolic variance. \n");
    printf("           -x xscalingfactor\n");
    printf("                Units are SI units (s) by default ; should the input data be in other units (MJD, ns, ...)\n");
    printf("                x option allows to properly normalize output : \n");
    printf("                scaling factor is one of : \n");
    printf("                    d : days  \n");
    printf("                    H : hours  \n");
    printf("                    M : minutes  \n");
    printf("                    m : millisecond \n");
    printf("                    u : microsecond \n");
    printf("                    n : nanosecond \n");
    printf("                    p : picosecond \n");
    printf("                  A file containing data as MJD.XXXXX vs freq_dev can be processed with : \n");
    printf("                  SigmaTheta -x d datafile\n\n");
    printf("           -h : this message\n\n");
    printf("    The input file SOURCE contains a 2-column table with time values (dates) in the first column\n");
    printf("         and normalized frequency deviation measurements in the second column.\n\n");
    printf("    A 7-column table is sent to the standard output and/or the output file with:\n");
    printf("            1st column: tau values\n");
    printf("            2nd column: (modified) Allan deviation estimate\n");
    printf("            3rd column: unbiased estimate\n");
    printf("            4th column: 2.5 %% bound\n");
    printf("            5th column: 16 %% bound\n");
    printf("            6th column: 84 %% bound\n");
    printf("            7th column: 97.5 %% bound.\n\n");
    printf("    The file TARGET.gnu is generated and used as an input to gnuplot.\n");
    printf("    The configuration file \".SigmaTheta.conf\" is taken into account \n");
    printf("        (e.g. selection of the modified Allan varance).\n\n");
    printf("    The file TARGET.pdf is the pdf file of the gnuplot graph\n");
    printf("        if the PDF option has been chosen in the configuration file.\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    printf("##############################################################################################################\n\n");
    }

int main(int argc, char *argv[]) {
    int i,j,nto,N,err,alpha[256];
    double tsas,asympt,tau[256], dev[256], avar[256], edf[256], bmin[256], bmax[256],bi1s[256],bx1s[256],adc[256],w[256];
    char scalex=0, scaley=0;   
    char fsw=0, fv=AVAR, command[32], source[MAXCHAR], outfile[MAXCHAR];
    char varchar='a', c;
    struct conf_int rayl;
    FILE *ofd;
    int index;
    size_t lensrc, stridx;

    opterr = 0;

    while ((c = getopt (argc, argv, "amHphx:")) != -1)
        switch (c) {
            case 'p':
                fsw=1;
                fv=PVAR;
                varchar=c;
                break;
            case 'H':
                fsw=1;
                fv=HVAR;
                varchar=c;
                break;
            case 'm':
                fsw=1;
                fv=MVAR;
                varchar=c;
                break;
            case 'a':
                fsw=1;
                fv=AVAR;
                varchar=c;
                break;
            case 'h':
                usage();
                break;
            case 'x':
                scalex=optarg[0];
                break;
            case '?':
                if (optopt == 'x'){
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                    usage();
                }
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return 1;
            default:
                abort ();
        }

    index=optind;
    if (index<argc){
        lensrc=strlen(strncpy(source,argv[index], MAXCHAR));
        printf ("# Input file %s\n", argv[index]);

        index++;
        if (index<argc){
            strncpy(outfile,argv[index], MAXCHAR);
            printf ("# Output file %s\n", outfile);
        }
        else {
            //
            // SOURCE is specified but not TARGET ; 
            // lets build a reasonable output file name :
            //
            stridx=lensrc-3;
            if (source[stridx++]=='y' && source[stridx++]=='k' && source[stridx++]=='t' ){
                strncpy(outfile, source, MAXCHAR);
                snprintf(outfile+lensrc-3, 5, "%cdev", varchar);
            }
            else{
                snprintf(outfile, lensrc+6, "%s.%cdev", source, varchar);
            }
            printf ("# Output file %s\n", outfile);
        }

        index++;
        if (index<argc)
        {
            fprintf(stderr,"Extra parameter.\n");
            usage();
        }
    }
    else{ // neither source nor target specified 
        // build a generic filename, needed for gener_gplt
        snprintf(outfile, MAXCHAR , "st_generic_target.%cdev", varchar);
    }

    if (strlen(source)==0){
        fprintf(stderr,"#\n#\n# No input file given , expecting data on stdin...\n#\n", argv[0]);
    }

    err=0;
    nto=load_ykt(source,scalex,scaley);

    if (nto==-1)
        printf("# File not found\n");
    else
    {
        if (nto<2)
        {
            printf("# Unrecognized file\n");
            if (nto==1)
                printf("# Use 1col2col command\n\n");
            usage();
        }
        else
        {
            err=init_flag();
            if (err==-1) printf("# ~/.SigmaTheta.conf not found, default values selected\n");
            if (err==-2) printf("# ~/.SigmaTheta.conf improper, default values selected\n");
            if (fsw)
            { 
                flag_variance=fv;
                //        printf("flag_variance=%d\n",flag_variance);
            }
            if (flag_variance&1)
            {
                flag_slopes[0]=1;
                //        printf("flag_slopes[0]=%d\n",flag_slopes[0]);
            }
            else
            {
                flag_slopes[0]=0;
                //        printf("flag_slopes[0]=%d\n",flag_slopes[0]);
            }
            if (flag_variance==3)
                flag_conf=0;
            /*        for (i=0;i<6;++i) printf("%d ",flag_slopes[i]);
                    printf("\n");*/
            N=serie_dev(nto, tau, dev);
            for(i=0;i<N;++i) avar[i]=dev[i]*dev[i];
            err=relatfit(N,tau,avar,tau,6);
            /*        printf("# Asymptote coefficients:\n");
                    printf("# tau^-1   \t tau^-1/2 \t tau^0    \t tau^1/2  \t tau^1\n");
                    for(i=0;i<5;++i) printf("%12.6e \t ",coeff[i]);
                    printf("\n");*/
            for(i=0;i<N;++i)
            {
                asympt=0;
                if (flag_variance)
                {
                    for(j=0;j<5;++j) /* The drift is not concerned */
                    {
                        tsas=coeff[j]*interpo(tau[i],j);
                        if (tsas>asympt)
                        {
                            asympt=tsas;
                            alpha[i]=2-j;
                        }
                    }
                }
                else
                {
                    for(j=1;j<5;++j) /* Neither the tau^-3/2 asymptote nor the drift are concerned */
                    {
                        tsas=coeff[j]*interpo(tau[i],j);
                        if (tsas>asympt)
                        {
                            asympt=tsas;
                            if (j==1) alpha[i]=+2;
                            else alpha[i]=2-j;
                        }
                    }
                }
            }
            avardof(N, tau, alpha, edf);
            if (flag_fit&0x8)
            {
                for(i=0;i<N;++i) w[i]=((double)1)/edf[i];
                err=relatfit(N,tau,avar,w,6);
                for(i=0;i<N;++i)
                {
                    asympt=0;
                    if (flag_variance)
                    {
                        for(j=0;j<5;++j) /* The drift is not concerned */
                        {
                            tsas=coeff[j]*interpo(tau[i],j);
                            if (tsas>asympt)
                            {
                                asympt=tsas;
                                alpha[i]=2-j;
                            }
                        }
                    }
                    else
                    {
                        for(j=1;j<5;++j) /* Neither the tau^-3/2 asymptote nor the drift are concerned */
                        {
                            tsas=coeff[j]*interpo(tau[i],j);
                            if (tsas>asympt)
                            {
                                asympt=tsas;
                                if (j==1) alpha[i]=+2;
                                else alpha[i]=2-j;
                            }
                        }
                    }
                }
                /*            for(i=0;i<N;++i) printf("%d \t ",alpha[i]);
                            printf("\n");*/
                avardof(N, tau, alpha, edf);
                /*        printf("# Asymptote coefficients:\n");
                        printf("# tau^-1   \t tau^-1/2 \t tau^0    \t tau^1/2  \t tau^1\n");
                        for(i=0;i<5;++i) printf("%12.6e \t ",coeff[i]);
                        printf("\n");*/
            }
            ofd=fopen(outfile, "w");

            switch(flag_variance)
            {
                case 1 : 
                    printf("# Tau       \t Mdev       \t Mdev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
                    fprintf(ofd,"# Tau       \t Mdev       \t Mdev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
                    break;

                case 2 :
                    printf("# Tau       \t Hdev       \t Hdev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
                    fprintf(ofd,"# Tau       \t Hdev       \t Hdev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
                    break;

                case 3 :
                    printf("# Tau       \t Pdev       \t Pdev unbiased\n");
                    fprintf(ofd,"# Tau       \t Pdev       \t Pdev unbiased\n");
                    break;

                default :
                    printf("# Tau       \t Adev       \t Adev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
                    fprintf(ofd,"# Tau       \t Adev       \t Adev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
            }
            for(i=0;i<N;++i) {
                rayl=raylconfint(edf[i]);
                bmin[i]=dev[i]*sqrt(edf[i])/rayl.sup_bound;
                bmax[i]=dev[i]*sqrt(edf[i])/rayl.inf_bound;
                rayl=raylconfint1s(edf[i]);
                bi1s[i]=dev[i]*sqrt(edf[i])/rayl.sup_bound;
                bx1s[i]=dev[i]*sqrt(edf[i])/rayl.inf_bound;
                adc[i]=dev[i]/rayl.mean;
                if (flag_variance!=3) {
                    printf("%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],dev[i],adc[i],bmin[i],bi1s[i],bx1s[i],bmax[i]);
                    fprintf(ofd,"%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],dev[i],adc[i],bmin[i],bi1s[i],bx1s[i],bmax[i]);
                }
                else {
                    printf("%12.6e \t %12.6e \t %12.6e\n",tau[i],dev[i],adc[i]);
                    fprintf(ofd,"%12.6e \t %12.6e \t %12.6e\n",tau[i],dev[i],adc[i]);
                }
            }
            fclose(ofd);

            if (flag_bias){
                for(i=0;i<N;++i) {
                    avar[i]=adc[i]*adc[i];
                }
            }

            for(i=0;i<N;++i) {
                w[i]=((double)1)/edf[i];
            }
            err=relatfit(N,tau,avar,w,6);
            printf("# Asymptote coefficients:\n");
            printf("# tau^-3/2 \t tau^-1   \t tau^-1/2 \t tau^0    \t tau^1/2  \t tau^1\n");
            for(i=0;i<6;++i) {
                printf("%12.6e \t ",coeff[i]);
            }    
            printf("\n");
            if (flag_graph) {
                /* Use of gnuplot for generating the graph as a ps file */
                err=gener_gplt(outfile,N,tau,dev,bmax,"unbiased",1);
                if (err) printf("# Error %d: ps file not created\n",err);
            }
        }
    }
}

