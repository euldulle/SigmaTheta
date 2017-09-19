/*   adgraph.c                                     F. Vernotte - 2010/12/31 */
/*   Plot the graph of the Allan Deviation estimates versus tau             */
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
    printf("Usage: ADGraph [-m] DevFILE FitFILE\n\n");
    printf("Plots the graph of the (modified) Allan Deviation estimates versus tau.\n\n");
    printf("The input file DevFILE contains a 2-column table with tau values (integration time) in the first column and deviation measurement in the second column.\n\n");
    printf("The input file FitFILE contains the 5 asymptote coefficients (from tau^-1 to tau^+1) in a 1-line 5-column table.\n\n");
    printf("The file DevFILE.gnu is generated for invoking gnuplot.\n");
    printf("The file DevFILE.ps is the postscript file of the gnuplot graph.\n\n");
    printf("If the option '-m' is selected, the variance is assumed to be the modified Allan variance. Otherwise, the variance is assumed to be the classical Allan variance.\n\n"); 
    printf("Sigma-Theta %s %s - UTINAM/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
    {
    char source[256], fitfile[256], command[32];
    int N,err,i,j,nc;
    double tau[32], adev[32], truc[32], bmax[32];
    char fv;

    fv=0;
    if ((argc<3)||(argc>4))
        {
        usage();
	exit(-1);
	}
    else
	if (argc==3)
        	{
		strcpy(source,*++argv);
		strcpy(fitfile,*++argv);
		}
	else
		{
		strcpy(command,*++argv);
		if (command[0]=='-')
			{
			if (!strcmp(command,"-m")) 
				{
				fv=1;
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
		else
			{
			usage();
			exit(-1);
			}
		}

    err=init_flag();
    if (err==-1) printf("# ~/.SigmaTheta.conf not found, default values selected\n");
    if (err==-2) printf("# ~/.SigmaTheta.conf improper, default values selected\n");
    flag_variance=fv;
    N=load_7col(source,tau,adev,truc,truc,truc,truc,bmax);
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
/* Use of gnuplot for generating the graph as a ps file */
	    err=gener_gplt(source,N,tau,adev,bmax);
	    if (err) printf("# Error %d: ps file not created\n",err);
	    }
	}
    }
		  
