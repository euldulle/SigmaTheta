/*   aduncert.c                                    F. Vernotte - 2010/12/31 */
/*   Computation of the uncertainties of a adev measurement serie           */
/*   without fitting nor plotting a graph                                   */
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
    printf("Usage: ADUncert [-m] DevFILE EdfFILE\n\n");
    printf("Computes the 95 %% (2 sigma) and the 68%% (1 sigma) confidence intervals of a sequence of (modified) Allan Deviations.\n\n");
    printf("The input file DevFILE contains a 2-column table with tau values (integration time) in the first column and (modified) Allan deviation measurement in the second column.\n\n");
    printf("The input file EdfFILE contains a 2-column table with tau values in the first column and the equivalent degrees of freedom (edf) of the deviation measurement in the second column.\n\n");
    printf("A 7-column table is sent to the standard output with:\n");
    printf("\t1st column: tau values\n");
    printf("\t2nd column: (modified) Allan deviation estimate\n");
    printf("\t3rd column: unbiased estimate\n");
    printf("\t4th column: 2.5 %% bound\n");
    printf("\t5th column: 16 %% bound\n");
    printf("\t6th column: 84 %% bound\n");
    printf("\t7th column: 97.5 %% bound.\n\n");
    printf("If the option '-m' is selected, the variance is assumed to be the modified Allan variance. Otherwise, the variance is assumed to be the classical Allan variance.\n\n"); 
    printf("A redirection should be used for saving the results in a TARGET file: ADUncert AdevFILE EdfFILE > TARGET\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
/* Compute the coefficients of the tau^-1, tau^-1/2, tau^0, tau^1/2 and tau of a {tau,adev} serie */
/* Input : tau and adev value file   */
/* Output: 5 asymptotes coefficients */
    {
      int i, j, N, N2, err, ineq;
      double tau[32], adev[32], truc[32], to2[32], edf[32], bmin[32], bmax[32],bi1s[32],bx1s[32],adc[32],w[32];
    char source[256], edffile[256], command[32];
    struct conf_int rayl;

    flag_variance=0;
    if ((argc<3)||(argc>4))
        {
        usage();
	exit(-1);
	}
    else
	if (argc==3)
        	{
		strcpy(source,*++argv);
		strcpy(edffile,*++argv);
		}
	else
		{
		strcpy(command,*++argv);
		if (command[0]=='-')
			{
			if (!strcmp(command,"-m")) 
				{
				flag_variance=1;
				strcpy(source,*++argv);
				strcpy(edffile,*++argv);
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

    N=load_3col(source,tau,adev,truc);
    if (N==-1)
      printf("# File %s not found\n",source);
    else
        {
        if (N<2)
	    {
	    printf("# Unrecognized file %s\n\n",source);
	    usage();
	    }
        else
	    {
	    N2=load_3col(edffile,to2,edf,truc);
	    if (N2>0)
	        {
	        ineq=0;
	        for (i=0;i<N2;++i) if ((tau[i]-to2[i])/tau[i]>1e-3) ++ineq;
	        }
	    else
	        {
		printf("# File %s not found\n",edffile);
		exit(-1);
		}
	    if ((N2!=N)||(N2<=0)||(ineq>0))
	        {
	        printf("# %s incompatible with %s\n\n",edffile,source);
	        usage();
	        }
	    if (flag_variance)
	    	printf("# Tau       \t Mdev       \t Mdev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
	    else
	    	printf("# Tau       \t Adev       \t Adev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
	    for(i=0;i<N;++i)
	        {
		rayl=raylconfint(edf[i]);
		bmin[i]=adev[i]*sqrt(edf[i])/rayl.sup_bound;
		bmax[i]=adev[i]*sqrt(edf[i])/rayl.inf_bound;
		rayl=raylconfint1s(edf[i]);
		bi1s[i]=adev[i]*sqrt(edf[i])/rayl.sup_bound;
		bx1s[i]=adev[i]*sqrt(edf[i])/rayl.inf_bound;
		adc[i]=adev[i]/rayl.mean;
		printf("%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],adev[i],adc[i],bmin[i],bi1s[i],bx1s[i],bmax[i]);
		}
	    }
	}
    }
		  
