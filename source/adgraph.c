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
    printf("The input file DevFILE contains a 7-column table with tau values (integration time) in the first column, deviation measurements in the second column, unbiased estimates in the 3rd, 68 and 95 %% confidence interval bounds in the following columns.\n\n");
    printf("The input file FitFILE contains the 6 asymptote coefficients (from tau^-3/2 to tau^+1) in a 1-line 6-column table.\n\n");
    printf("The file DevFILE.gnu is generated for invoking gnuplot.\n");
    printf("The file DevFILE.pdf is the pdf file of the gnuplot graph (if the PDF option has been chosen in the configuration file).\n\n");
    printf("If the option '-m' is selected, the variance is assumed to be the modified Allan variance. Otherwise, the variance is assumed to be the classical Allan variance.\n\n"); 
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
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
//    N=load_adev(source,tau,adev);
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
	    err=gener_gplt(source,N,tau,adev,bmax,"unbiased");
	    if (err) printf("# Error %d: ps file not created\n",err);
	    }
	}
    }
		  
