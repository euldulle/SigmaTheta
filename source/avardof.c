/*   avardof.c                                     F. Vernotte - 2010/10/24 */
/*   Computation of the degrees of freedom of an Allan variance sequence    */
/*   following the paper "Uncertainty of Stability Variances Based on       */
/*   Finite Differences" by C. Greenhall and W. Riley (35th PTTI, San       */
/*   Diego, 2003, pp. 267-278)                                              */
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
#include "sigma_theta.h"

#define db(x) ((double)(x))
#define sisig(x) ( (x) == 0 ) ? (db(0)) : (  ( (x) > 0 ) ? (db(1)) : (db(-1))  )

void usage(void)
/* Help message */
    {
    printf("Usage: AVarDOF [-a|m|h|p] SOURCE\n\n");
    printf("Compute the degrees of freedom of the variance computations over a sequence of tau values.\n\n");
    printf("The input file SOURCE contains a 2-column table with tau values (integration time) in the first column and the exponent of the power law of the dominating noise in the second column.\n\n");
    printf("The first tau value is assumed to be equal to the sampling step.\n");
    printf("The last tau value is assumed to be equal to the half of the whole time sequence duration.\n\n");
    printf("A 2-column table containing tau values (integration time) in the first column and the equivalent degrees of freedom in the second column is sent to the standard output.\n\n");
    printf("If the option '-a' is selected, the variance is assumed to be the classical Allan variance (default).\n"); 
    printf("If the option '-m' is selected, the variance is assumed to be the modified Allan variance.\n"); 
    printf("If the option '-h' is selected, the variance is assumed to be the Hadamard variance.\n"); 
    printf("If the option '-p' is selected, the variance is assumed to be the parabolic variance.\n\n"); 
    printf("A redirection should be used for saving the results in a TARGET file: AVarDOF SOURCE > TARGET\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
/* Computation of the degrees of freedom of the Allan variance estimates    */
/* (based on "Uncertainty of stability variances based on finite            */
/* differences" by C.A.Greenhall and W.J.Riley, 35th PTTI)                  */
/* input:  2-column table containing the tau values and the dominant FM     */
/*         power law exponent for each tau value                            */
/* output: 2-column table containing the tau values and the equivalent      */
/*         degrees of freedom of each Allan variance measurement            */
    {
    int i, N, alphint[32];
    double tau[32], edf[32], alpha[32];
    char gm[256], source[256], command[32];
    FILE *ofd;

    flag_variance=0;
    if ((argc<2)||(argc>3))
        {
	usage();
	exit(-1);
        }
    else
	{
	if (argc==2)
		{
		strcpy(source,*++argv);
		flag_variance=AVAR;
		}
	else
		{
		strcpy(command,*++argv);
		if (command[0]=='-')
			{
			if ((!strcmp(command,"-a"))||(!strcmp(command,"-m"))||(!strcmp(command,"-h"))||(!strcmp(command,"-p"))) 
				{
				switch(command[1])
					{
					case 'p':
						flag_variance=PVAR;
						break;
					case 'h':
						flag_variance=HVAR;
						break;
					case 'm':
						flag_variance=MVAR;
						break;
					case 'a':
					default:
						flag_variance=AVAR;
					}
				strcpy(source,*++argv);
				}
			else
				{
				printf("Unknown option '%s'\n",command);
				usage();
				exit(-1);
				}
			}
		else
			{
			usage();
			exit(-1);
			}
		}
	}
/*if (source[0]=='-')*/
/*			{*/
/*			if (!strcmp(source,"-m"))*/
/*				{*/
/*				flag_variance=1;*/
/*				strcpy(source,*++argv);*/
/*				}*/
/*			else*/
/*				{*/
/*				printf("Unknown option '%s'\n",source);*/
/*				usage();*/
/*				exit(-1);*/
/*				}*/
/*			}*/
/*		else*/
/*			{*/
/*			strcpy(command,*++argv);*/
/*			if (command[0]=='-')*/
/*				{*/
/*				if (!strcmp(command,"-m")) flag_variance=1;*/
/*				else*/
/*					{*/
/*					printf("Unknown option '%s'\n",command);*/
/*					usage();*/
/*					exit(-1);*/
/*					}*/
/*				}*/
/*			else*/
/*				{*/
/*				usage();*/
/*				exit(-1);*/
/*				}*/
/*			}*/
/*		}*/
/*	}*/

    N=load_adev(source,tau,alpha);
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
	    for(i=0;i<N;++i) alphint[i]=(int)alpha[i];
	    if (flag_variance!=PVAR)
	    	avardof(N, tau, alphint, edf);
	    else
		pvardof(N, tau, alphint, edf);
            for(i=0;i<N;++i)
	        printf("%24.16e \t %24.16e\n",tau[i],edf[i]);
            }
	}
    }
		  
