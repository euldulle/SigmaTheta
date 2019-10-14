/*   mdev.c                                        F. Vernotte - 2014/09/24 */
/*   Computation of the Modified Allan Deviation (MDev)                     */
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
#include <math.h>
#include <unistd.h>
#include <stdint.h>
#include "sigma_theta.h"

#define DATAMAX 16384
#define GRANMAX 67108864

void usage(uint8_t *message)
/* Help message */
    {
    printf("\n   Sigma-Theta %s %s - UTINAM/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE)\n\n",st_version,st_date);
    printf("Usage: %s [-t] SOURCE\n\n", message);
    printf(" Computes the Modified Allan Deviations of a sequence of \n");
    printf(" normalized frequency deviation measurements.\n\n");
    printf("  options  \n");
    printf("     -t : outputs TDEV instead of MDEV\n\n");

    printf("  input file SOURCE is a 2-column table : \n");
    printf("     . time values (timestamps, in seconds) in the first column \n");
    printf("     . normalized frequency deviation measurements in the second column.\n\n");
    printf("  output consists of a 2-column table containing : \n");
    printf("     . tau values (integration time) in the first column  \n");
    printf("     . Modified Allan deviation (or TDEV if -t) measurement in the second column \n");
    printf("    which is sent to the standard output.\n\n");
    printf(" A redirection should be used to save the results in a TARGET file: MDeV SOURCE > TARGET\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
	/* Compute MDEV serie from tau=tau0 (tau0=sampling step) to tau=N*tau0/3 (N=number of samples) by octave (log step : tau_n+1=2*tau_n) */
	/* Input : file name of the normalized frequency deviation samples*/
	/* Output : tau \t MDev(tau) */
	/*          (for tau=tau0 to tau=N*tau0/3 by octave) */
{
	int i,nbv,N,nto,tomax;
	long int dtmx;
	char source[256], gm[100];
    char scalex=0, scaley=0;
	FILE *ofd;
	double v1,v2,smpt,rslt,tau[256],dev[256];
	uint8_t stddev=0, tdev=0;
    int8_t c;
	extern char *optarg;
	extern int opterr;

	while(1) {
		c=getopt(argc, argv, "t");
		if (c=='?') {
			usage(argv[0]);
		}

		if (c==-1) break;
		switch(c) {
			case 't':
				tdev=1;
				break;

			default:
				usage(argv[0]);
				break;
		}
	}

	if (optind>=argc)
	{
		usage(argv[0]);
		exit(-1);
	}
	else
		strcpy(source,argv[optind]);

	N=load_ykt(source,scalex,scaley);
	if (N==-1)
		printf("# File not found\n");
	else
	{
		if (N<2)
		{
			printf("# Unrecognized file\n\n");
			if (N==1)
				printf("# Use 1col2col command\n\n");
			usage(argv[0]);
		}
		else
		{
			flag_log_inc=1;
			flag_variance=1;
			log_inc=(double)2;
			nto=serie_dev(N, tau, dev);
			for(i=0;i<nto;++i){
                if (tdev==1)
                    printf("%24.16e \t %24.16e\n",tau[i],tau[i]*dev[i]/sqrt(3));
                else 
                    printf("%24.16e \t %24.16e\n",tau[i],dev[i]);
			}
		}
	}
}

