/*   asymptote.c                                   F. Vernotte - 2010/10/22 */
/*                                              modified by FV - 2010/12/30 */
/*   Computation of the asymptotes of a adev measurement serie              */
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
    printf("Usage: Asymptote SOURCE [edfFILE]\n\n");
    printf("Computes the tau^-3/2, tau^-1, tau^-1/2, tau^0, tau^1/2 and tau asymptotes of a sequence of (modified) Allan Deviations.\n\n");
    printf("The input file SOURCE contains a 2-column table with tau values (integration time) in the first column and (modified) Allan deviation measurements in the second column.\n\n");
    printf("The optional input file edfFILE contains a 2-column table with tau values in the first column and the equivalent degrees of freedom (edf) of the Allan deviation measurement in the second column.\n\n");
    printf("The Allan deviation estimates are weigted by the tau-values in the default case and by the inverse of the \"edf\" if the optional file edfFILE has been used.\n\n");
    printf("If the configuration file \".SigmaTheta.conf\" contains the line \"Unbiased estimates: ON\", the fit is performed over the unbiased estimates.\n\n");
    printf("The 6 asymptote coefficients are sent to the standard output separated by a tabulation.\n\n");
    printf("A redirection should be used for saving the results in a TARGET file: Asymptote SOURCE [edfFILE] > TARGET\n\n");
    printf("Sigma-Theta %s %s - UTINAM/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
/* Compute the coefficients of the tau^-1, tau^-1/2, tau^0, tau^1/2 and tau of a {tau,adev} serie */
/* Input : tau and adev value file   */
/* Output: 5 asymptotes coefficients */
    {
    int i,nbv,N,N2,err,ineq;
    double tau[32], avar[32], ubad[32], to2[32], edf[32], wght[32];
    char gm[256], source[256], weight[256];
    FILE *ofd;

    if (argc<2)
        {
	usage();
	exit(-1);
        }
    else
        {
	strcpy(source,*++argv);
	if (argc>2) 
	  strcpy(weight,*++argv);
	}
    
    if (init_flag()<0) printf("# ~/.SigmaTheta.conf not found or improper, default values selected\n");
    N=load_3col(source,tau,avar,ubad);
    if (N==-1)
        printf("# File not found\n");
    else
        {
        if (N<2)
	    {
	    printf("# %s: Unrecognized file\n\n",source);
	    usage();
	    }
        else
	    {
	    if (argc>2)
		{
	        N2=load_adev(weight,to2,edf);
	        if (N2>0)
	            {
	            ineq=0;
	            for (i=0;i<N2;++i) if ((tau[i]-to2[i])/tau[i]>1e-3) ++ineq;
	            }
	        if ((N2!=N)||(N2<=0)||(ineq>0))
	            {
		    printf("# %s incompatible with %s: weighted by 1/tau\n",weight,source);
		    for(i=0;i<N;++i) wght[i]=tau[i];
	            }
	        else
		    for(i=0;i<N;++i) wght[i]=1/edf[i];
		}
	    else
		for(i=0;i<N;++i) wght[i]=tau[i];
            if ((flag_bias)&&(ubad[0]!=-1)) 
                for(i=0;i<N;++i)avar[i]=ubad[i];
	    for(i=0;i<N;++i) avar[i]*=avar[i];
	    err=relatfit(N,tau,avar,wght,6);
	    printf("# Asymptote coefficients:\n");
	    printf("# tau^-3/2 \t tau^-1   \t tau^-1/2 \t tau^0    \t tau^1/2  \t tau^1\n");
	    for(i=0;i<6;++i) printf("%24.16e \t ",coeff[i]);
	    printf("\n");
	    }
	}
    }
		  
