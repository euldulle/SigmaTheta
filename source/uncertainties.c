/*   uncertainties.c                              F. Vernotte - 2010/10/24  */
/*   Computation of the uncertainties of a adev or mdev measurement serie   */
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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "sigma_theta.h"

#define db(x) ((double)(x))
#define sisig(x) ( (x) == 0 ) ? (db(0)) : (  ( (x) > 0 ) ? (db(1)) : (db(-1))  )

void usage(void)
/* Help message */
    {
    printf("##############################################################################################################\n\n");
    printf(" uncertainties : a tool from the SigmaTheta suite\n\n");
    printf("     Usage: uncertainties [-mnch] [-o outfile] [SOURCE [TARGET]]\n\n");
    printf("     Computes the 95 %% confidence intervals of a sequence of \n");
    printf("     (modified) Allan Deviations, the asymptotes and plot a graph.\n\n");
    printf("     The input file SOURCE contains a 2-column table with tau values (integration time)\n");
    printf("     in the first column and (modified) Allan deviation measurement in the second column.\n");
    printf("     If the input file is missing data is read out of stdin.\n\n");

    printf("     The first tau value is assumed to be equal to the sampling step.\n");
    printf("     The last tau value is assumed to be equal to the half.\n");
    printf("        of the whole time sequence duration.\n\n");
    printf("     A 7-column table is sent to the standard output with:\n");
    printf("         1st column: tau values\n");
    printf("         2nd column: (modified) Allan deviation estimate\n");
    printf("         3rd column: unbiased estimate\n");
    printf("         4th column: 2.5 %% bound\n");
    printf("         5th column: 16 %% bound\n");
    printf("         6th column: 84 %% bound\n");
    printf("         7th column: 97.5 %% bound.\n\n");
    printf("    The file TARGET.gnu is generated for invoking gnuplot. \n");
    printf("    The configuration file \".SigmaTheta.conf\" is taken into account.\n");
    printf("    The file TARGET.ps is the postscript file of the gnuplot graph.\n\n");
    printf("    By default, the variance is assumed to be the one\n"); 
    printf("    selected in the configuration file \".SigmaTheta.conf\".\n\n"); 
    printf("    Options :\n");
    printf("        -m : use the modified Allan variance.\n");
    printf("        -c : use the classical Allan variance.\n");
    printf("        -L : specific ltfb output format (numline, tau, avar, 2.5%, 97.5%)\n");
    printf("             consequently, -L inhibits plot generation by sigmatheta \n");
    printf("        -n : the plot is not built (gnuplot file(s) might still be generated).\n");
    printf("        -N : dont write asymptote coefficients to stdout.\n");
    printf("        -o output : when reading from stdin, TARGET cannot be passed as arg : \n");
    printf("                    -o is the way user can specify the TARGET output \n");
    printf("                    If not specified, a random filename will be generated as TARGET.\n");
    printf("                    (Specifying TARGET as arg along with SOURCE on the command line\n");
    printf("                     is still possible for backward compatibility)\n");
    printf("       Options for gnuplot output :\n");
    printf("        -p : insert png build in gnuplot file.\n\n");
    printf("        -d : insert pdf build in gnuplot file.\n\n");
    printf("        -X : insert x11 build in gnuplot file.\n\n");
    printf("        -h : this message.\n\n");
    printf("           SigmaTheta %s %s \n",st_version,st_date);
    printf("           FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n");
    printf("##############################################################################################################\n\n");
    exit(-1);
    }

int main(int argc, char *argv[])
    /* Compute the coefficients of the tau^-3/2, tau^-1, tau^-1/2, tau^0, tau^1/2 and tau of a {tau,adev} serie */
    /* Input : tau and adev value file   */
    /* Output: 6 asymptotes coefficients */
{
    int i,j,N,err,alpha[32];
    double tsas,asympt,tau[32], adev[32], avar[32], edf[32], bmin[32], bmax[32],bi1s[32],bx1s[32],adc[32],w[32];
    char pre_flag_v, source[MAXCHAR]="", tmpoutfile[16]="st_unc_XXXXXX", outfile[MAXCHAR]="", command[32], varname[32];
    struct conf_int rayl;
    FILE *ofd;
    int8_t c, stdo;
    uint8_t doplot=GPLDOPLOT; // by default, 1 (GPDDOPLOT flag on) will build gnuplot file and run it
    //
    // gadgets for gnss certificates at ltfb (might still prove useful elsewhere though)
    uint8_t asymout=1;        // by default, output asymptotes coefficients on stdout ; specify -N option for no output
                              //
    uint8_t numlineout=0;     // by default, dont output numline to stdout ; specify -N option for no output

    char flagchar='a';
    int index, stridx;
    size_t lensrc;

    opterr = 0;

    while ((c = getopt (argc, argv, ":cmnLNho:pdX")) != -1)
        switch (c)
        {
            case 'c': // classical Allan Dev
                pre_flag_v=2;
                flagchar='a';
                break;

            case 'L': // output line number as column 1 in the lines of dev and uncertainties values for each tau
                      // WARNING : this disables gnuplot file and plot generation
                numlineout=1;
                break;

            case 'm': // modified Allan Dev
                pre_flag_v=1;
                flagchar='m';
                break;

            case 'n': // clear DOPLT flag = no actual plot (gnuplot file might still be generated but not run)
                doplot&=~GPLDOPLOT;
                break;

            case 'N': // prevents the output of asymptote coefficients to stdout 
                asymout=0;
                break;

            case 'p': // term png : insert png generation in gnuplot file
                doplot=doplot|GPLPNG;
                break;

            case 'd': // term pdf : insert pdf generation in gnuplot file
                doplot=doplot|GPLPDF;
                break;

            case 'X': // term pdf : insert x11 generation in gnuplot file 
                doplot=doplot|GPLX11;
                break;

            case 'o':
                strncpy(outfile,optarg, MAXCHAR);
                break;

            case 'h':
                usage();
                break;

            case ':':
                fprintf (stderr, "option -%c needs a parameter.\n", optopt);
                usage();
                break;

            case '?':
                if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                return 1;
            default:
                abort ();
        }

    index=optind;
    if (index<argc){
        lensrc=strlen(strncpy(source,argv[index], MAXCHAR));
        printf ("Input file %s\n", argv[index]);

        index++;
        if (index<argc){
            printf ("Output file %s\n", argv[index]);
            strncpy(outfile,argv[index], MAXCHAR);
        }
        else {
            if (stdo != 1){
                stridx=lensrc-3;
                if (source[stridx++]=='d' && source[stridx++]=='e' && source[stridx++]=='v' ){
                    snprintf(outfile, lensrc+1, "%su", source);
                }
                else{
                    snprintf(outfile, lensrc+5, "%s.%cdevu", source, flagchar);
                }
            }
        }

        index++;
        if (index<argc)
        {
            fprintf(stderr,"Extra parameter.\n");
            usage();
        }
    }

    if (strlen(source)==0){
        fprintf(stderr,"#\n#\n# No input file given , expecting data on stdin...\n#  (%s -h to show usage)\n", argv[0]);
        mktemp(tmpoutfile);
        if (outfile[0]==0){ // if -o file was not given, outfile[0] is 0 
                            // so we need to build an output filename
            snprintf(outfile,MAXCHAR,"%s.%cdevu", tmpoutfile, flagchar);
            }
        stdo=1;
    }

    printf("# Output file : %s\n", outfile);

    N=load_adev(source,tau,adev);
    if (N==-1)
        printf("# File not found\n");
    else
    {
        if (N<2)
        {
            printf("# Unrecognized file\n\n");
            usage();
        }
        else
        {
            err=init_flag();

            if (err==-1) 
                printf("# ~/.SigmaTheta.conf not found, default values selected\n");

            if (err==-2) 
                printf("# ~/.SigmaTheta.conf improper, default values selected\n");

            if (pre_flag_v==1) 
                flag_variance=1;
            else if (pre_flag_v==2) 
                flag_variance=0;

            for(i=0;i<N;++i) avar[i]=adev[i]*adev[i];
            err=relatfit(N,tau,avar,tau,6);
            for(i=0;i<N;++i)
            {
                asympt=0;
                /*              for(j=0;j<5;++j)
                                {
                                tsas=coeff[j]*interpo(tau[i],j);
                                if (tsas>asympt)
                                {
                                asympt=tsas;
                                if (j==0) alpha[i]=+2;
                                else alpha[i]=1-j;
                                }
                                }*/
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
                    /*                    for(j=0;j<4;++j)
                                          {
                                          tsas=coeff[j]*interpo(tau[i],j);
                                          if (tsas>asympt)
                                          {
                                          asympt=tsas;
                                          if (j==0) alpha[i]=+2;
                                          else alpha[i]=1-j;
                                          }
                                          }*/
                }
                avardof(N, tau, alpha, edf);
            }
            ofd=fopen(outfile, "w");

            if (ofd==NULL){
                fprintf(stderr,"Cannot open file %s for writing\nRedirecting to stderr\n");
                ofd=stderr;
                }

            strncpy(varname, "Adev", 5);
            if (flag_variance) {
                strncpy(varname, "Mdev", 5);
            }
            if (numlineout==0){ // regular output header with all bounds
                printf("# Tau       \t %4s       \t Adev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n", varname);
                fprintf(ofd,"# Tau       \t %4s       \t Adev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n", varname);
                }
            else{               // short output header only with 2 sigma bounds (ltfb)
                printf("# num \tTau       \t %4s      \t 2.5 %% bound \t 97.5 %% bound\n", varname);
                fprintf(ofd,"# num \tTau       \t %4s     \t 2.5 %% bound \t 97.5 %% bound\n", varname);
                }

            double plot_upperbound=1e-99, plot_lowerbound=1e99;
            for(i=0;i<N;++i) {
                rayl=raylconfint(edf[i]);
                bmin[i]=adev[i]*sqrt(edf[i])/rayl.sup_bound;
                bmax[i]=adev[i]*sqrt(edf[i])/rayl.inf_bound;

                if (bmin[i]<plot_lowerbound)
                    plot_lowerbound=bmin[i];

                if (bmax[i]>plot_upperbound)
                    plot_upperbound=bmax[i];

                rayl=raylconfint1s(edf[i]);
                bi1s[i]=adev[i]*sqrt(edf[i])/rayl.sup_bound;
                bx1s[i]=adev[i]*sqrt(edf[i])/rayl.inf_bound;
                adc[i]=adev[i]/rayl.mean;

                if (numlineout!=0){
                    //
                    // specific output format to be used in ltfb certificate processing
                    // numline, tau, adev, 2.5%, 97.5% bounds
                    //
                    printf("%d \t%8.2f \t %6.2e \t %6.2e \t %6.2e\n",i+1,tau[i],adev[i],bmin[i],bmax[i]);
                    fprintf(ofd, "%d \t%8.2f \t %6.2e \t %6.2e \t %6.2e\n",i+1,tau[i],adev[i],bmin[i],bmax[i]);
                }
                else{
                    // regular output, tau, adev, unbiased, 2.5%, 16%, 84%, 97.5% bounds
                printf("%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],adev[i],adc[i],bmin[i],bi1s[i],bx1s[i],bmax[i]);
                fprintf(ofd,"%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],adev[i],adc[i],bmin[i],bi1s[i],bx1s[i],bmax[i]);
                }
            }
            //
            // output (commented out) gnuplot command in outfile for use by external process
            //
            printf("# set yrange [%6.1e:%6.1e]\n", plot_lowerbound*.9, plot_upperbound*1.1);
            fprintf(ofd, "# set yrange [%6.1e:%6.1e]\n", plot_lowerbound*.9, plot_upperbound*1.1);
            if (ofd!=stderr)
                fclose(ofd);

            if (flag_bias)
                for(i=0;i<N;++i) 
                    avar[i]=adc[i]*adc[i];

            for(i=0;i<N;++i) 
                w[i]=((double)1)/edf[i];

            err=relatfit(N,tau,avar,w,6);
            if (asymout!=0){
                printf("# Asymptote coefficients:\n");
                printf("# tau^-3/2 \t tau^-1   \t tau^-1/2 \t tau^0    \t tau^1/2  \t tau^1\n");

                for(i=0;i<6;++i) 
                    printf("%12.6e \t ",coeff[i]);

                printf("\n");
                }
            
            /* Use of gnuplot for generating the graph as a ps file */

            if (numlineout==0){
                err=gener_gplt(outfile,N,tau,adev,bmax,"unbiased",doplot);
                if (err) 
                    printf("# Error %d: ps file not created\n",err);
                /*	    printf("tau: ");
                        for(i=0;i<N;++i) printf("%12.6e \t ",tau[i]);
                        printf("\n");
                        printf("edf: ");
                        for(i=0;i<N;++i) printf("%12.6e \t ",edf[i]);
                        printf("\n");*/
                }
            else{
                fprintf(stderr,"# Not generating gplot file since line numbers were requested on data lines with -L\n");
            }
        }
    }
}

