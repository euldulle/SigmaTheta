/*   adev.c                                        F. Vernotte - 2010/10/20 */
/*   Computation of the Allan Deviation (ADev)                              */
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
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include "sigma_theta.h"

//	void usage(void)
//	/* Help message */
//		{
//		printf("Usage: ADev SOURCE\n\n");
//		printf("Computes the Allan Deviations of a sequence of normalized frequency deviation measurements.\n\n");
//		printf("The input file SOURCE contains a 2-column table with time values (dates) in the first column and normalized frequency deviation measurements in the second column.\n\n");
//		printf("A 2-column table containing tau values (integration time) in the first column and Allan deviation measurement in the second column is sent to the standard output.\n\n");
//		printf("A redirection should be used for saving the results in a TARGET file: ADeV SOURCE > TARGET\n\n");
//		printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
//		}

void usage(void)
/* Help message */
    {                                                                                                                      
	printf("##############################################################################################################\n\n");
    printf(" ADev : a tool from the SigmaTheta suite\n\n");
    printf("     Usage: ADev [-ch] [-x xscalingfactor] [SOURCE [TARGET]]\n\n");
    printf("        Computes the Allan Deviations of a sequence of \n");
    printf("        normalized frequency deviation measurements.\n\n");
    printf("      Default behaviour (no file specified) is a filter, taking stdin as input and stdout as output.\n\n");
	printf("      If SOURCE and TARGET are specified, output goes to TARGET unless -c option (output to stdout) is given \n");
	printf("      If only SOURCE is specified, output goes to SOURCE.adev unless -c option (output to stdout) is given \n\n");
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
    printf("                  A file containing data as MJD.XXXXX vs freq_dev can be processed with \n");
    printf("                  ADev datafile -x d  \n\n");
    printf("           -c : output to stdout only, TARGET file is ignored even if specified ;\n");
    printf("                this is the default if SOURCE is unspecified (stdin)\n\n");
    printf("           -h : this message\n\n");
    printf("     Input consists in an N line / 2 column table with time values (dates) in the first column\n");
    printf("                                   and normalized frequency deviation samples in the second column.\n\n");
    printf("     Ouput is a 2 column table with averaging times tau in the first column \n");
    printf("                                   and Allan deviation measurement in column 2.\n\n");
    printf("           SigmaTheta %s %s \n",st_version,st_date);
    printf("           FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n");                                   
	printf("##############################################################################################################\n\n");
    exit(-1);
    }


int main(int argc, char *argv[])
    /* Compute ADEV serie from tau=tau0 (tau0=sampling step) to tau=N*tau0/2 (N=number of samples) by octave (log step : tau_n+1=2*tau_n) */
    /* Input : file name of the normalized frequency deviation samples*/
    /* Output : tau \t ADev(tau) */
    /*          (for tau=tau0 to tau=N*tau0/2 by octave) */
{
    int err,i,nbv,N,nto,tomax;
    long int dtmx;
    char gm[100], scalex=0;
    FILE *ofd;
    double v1,v2,smpt,rslt,tau[256],dev[256];
    uint8_t stdo=0, index;
    int8_t c;
    int16_t stridx, lensrc;
    uint8_t source[MAXCHAR]="", outfile[MAXCHAR]="";

    opterr = 0;

    while ((c = getopt (argc, argv, "hcx:")) != -1)
        switch (c)
        {
            case 'c':
                stdo = 1;
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

        index++;
        if (index<argc){
            strncpy(outfile,argv[index], MAXCHAR);
            printf ("# Output file %s\n", outfile);
        }
        else {
            if (stdo != 1){
                stridx=lensrc-3;
                if (source[stridx++]=='y' && source[stridx++]=='k' && source[stridx++]=='t' ){
                    strncpy(outfile, source, MAXCHAR);
                    snprintf(outfile+lensrc-3, 5, "adev");
                }
                else{
                    snprintf(outfile, lensrc+6, "%s.adev", source);
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
        fprintf(stderr,"#\n#\n# No input file given , expecting data on stdin...\n# %s -h to show usage )\n#\n", argv[0]);
        stdo=1;
        }

    err=0;

    N=load_ykt(source,scalex,0);

    if (N==-1)
        printf("# File %s not found\n", source);
    else
    {
        if (N<2)
        {
            printf("# not enough columns found in file %s \n\n", source);
            if (N==1)
                printf("# only one column found (tip : use 1col2col command)\n\n");
            usage();
        }
        else {
            printf ("# Input file %s = %d lines\n", source, N);
            if (stdo==1){
                ofd=stdout;
                printf("# Output to stdout: \n#\n");
                }
            else{
                ofd=fopen(outfile, "w");
                }

            err=init_flag();
            if (err==-1) printf("# ~/.SigmaTheta.conf not found, default values selected\n");
            if (err==-2) printf("# ~/.SigmaTheta.conf improper, default values selected\n");
            if (err) {
                flag_log_inc=1;
                log_inc=(double)2;
            }
            flag_variance=AVAR;

            nto=serie_dev(N, tau, dev);

            for(i=0;i<nto;++i)
                fprintf(ofd,"%24.16e \t %24.16e\n",tau[i],dev[i]);

            printf ("# Output file %s = %d lines\n", outfile, nto);

			if (ofd!=stdout)
				fclose(ofd);
        }
    }
}

