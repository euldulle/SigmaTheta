/*   asym2alpha.c                                  F. Vernotte - 2010/12/30 */
/*   Finds the dominating alpha power law noise versus tau.                 */
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
    printf("Usage: Asym2Alpha [-a|m|h|p] DevFILE FitFILE\n\n");
    printf("Finds the dominating power law noise (alpha) versus tau.\n\n");
    printf("The input file DevFILE contains a 2-column table with tau values (integration time) in the first column and Allan, modified Allan, Hadamard or Parabolic deviation measurements in the second column.\n\n");
    printf("The input file FitFILE contains the 6 asymptote coefficients (from tau^-3/2 to tau^+1) in a 1-line 6-column table.\n\n");
    printf("A a 2-column table with tau values in the first column and the dominating power law noise (from alpha=-2 to +2) in the second column is sent to the standard output.\n\n");
    printf("A redirection should be used for saving the results in a TARGET file: Asym2Alpha DevFILE FitFILE > TARGET\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
    {
    char source[256], fitfile[256], command[32];
    int flt,N,nto,i,j,nc,alpha[32];
    double tau[32], avar[32], truc[32], tsas, asympt;
    FILE *ofd;

    if ((argc<3)||(argc>4))
        {
        usage();
	exit(-1);
	}
    else
	if (argc==3)
        	{
		flag_variance=AVAR;
		strcpy(source,*++argv);
		strcpy(fitfile,*++argv);
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
				strcpy(fitfile,*++argv);
				}
			else
				{
				printf("Unknown option '%s'\n",command);
				usage();
				exit(-1);
				}
			}
		}

    N=load_3col(source,tau,avar,truc);
    if (N==-1)
      printf("# File %s not found\n",source);
    else
        {
        if (N<2)
	    {
	    printf("# %s: unrecognized file\n\n", source);
	    usage();
	    }
        else
	    {
            nc=load_coef(fitfile);
	    if (nc==-1)
	      {
              printf("# File %s not found\n",fitfile);
	      exit(-1);
	      }
	    else
	      if (nc!=6)
		{
	        printf("# %s: unrecognized file\n\n", fitfile);
	        usage();
		exit(-1);
		}
            for(i=0;i<N;++i)
                {
                asympt=0;
		if ((flag_variance==MVAR)||(flag_variance==PVAR))
			{
			flt=1;
			for(j=0;j<6;++j)
                    		{
                    		tsas=coeff[j]*interpo(tau[i],j);
                    		if (tsas>asympt)
                        		{
                        		asympt=tsas;
                        		alpha[i]=2-j;
					if (j==5) alpha[i]=-2;
                         		}
                    		}
			}
		else
			{
			flt=2;
	                for(j=1;j<6;++j) 
                    		{
                    		tsas=coeff[j]*interpo(tau[i],j);
                    		if (tsas>asympt)
                        		{
                        		asympt=tsas;
                        		if (j==1) alpha[i]=+2;
                        		else alpha[i]=2-j;
					if (j==5) alpha[i]=-2;
                        		}
                    		}
	                }
                }
//	    fprintf(stdout,"flt: %d\n",flt);
	    for(i=0;i<N;++i) printf("%24.16e \t %d\n",tau[i],alpha[i]);
	    }
	}
    }
		  
