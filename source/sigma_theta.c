/*   sigma_theta.c                                F. Vernotte - 2010/10/25  */
/*                                    Modified by FV for MDEV - 2014/10/01  */
/*       Modified for selecting the deviation by argument, FV - 2015/06/25  */
/*									    */
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
#include <math.h>
#include "sigma_theta.h"

#define db(x) ((double)(x))
#define sisig(x) ( (x) == 0 ) ? (db(0)) : (  ( (x) > 0 ) ? (db(1)) : (db(-1))  )

void usage(void)
/* Help message */
    {
    printf("Usage: SigmaTheta [-m|h|p] SOURCE TARGET\n\n");
    printf("Computes the Allan, Modified Allan, Hadamard or Parabolic Deviation of a sequence of normalized frequency deviation measurements.\n\n");
    printf("The input file SOURCE contains a 2-column table with time values (dates) in the first column and normalized frequency deviation measurements in the second column.\n\n");
    printf("A 7-column table is sent to the standard output with:\n");
    printf("\t1st column: tau values\n");
    printf("\t2nd column: deviation estimate\n");
    printf("\t3rd column: unbiased estimate\n");
    printf("\t4th column: 2.5 %% bound\n");
    printf("\t5th column: 16 %% bound\n");
    printf("\t6th column: 84 %% bound\n");
    printf("\t7th column: 97.5 %% bound.\n\n");
    printf("The file TARGET.gnu is generated for invoking gnuplot. The configuration file \".SigmaTheta.conf\" is taken into account (e.g. selection of the modified Allan varance).\n");
    printf("The file TARGET.pdf is the pdf file of the gnuplot graph (if the PDF option has been chosen in the configuration file).\n");
    printf("If the option '-m' is selected, the modified Allan variance is invoked. \n");
    printf("If the option '-h' is selected, the Hadamard variance is invoked. \n");
    printf("If the option '-p' is selected, the Parabolic variance is invoked. \n");
    printf("Otherwise, the classical Allan variance is invoked.\n\n");    
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
    {
    int i,j,nto,N,err,alpha[256];
    double tsas,asympt,tau[256], dev[256], avar[256], edf[256], bmin[256], bmax[256],bi1s[256],bx1s[256],adc[256],w[256];
    char fsw, fv, command[32], source[256], outfile[256];
    struct conf_int rayl;
    FILE *ofd;

    fsw=0;
    if ((argc<3)||(argc>4))
        {
        usage();
	exit(-1);
	}
    else
	if (argc==3)
        	{
		strcpy(source,*++argv);
		strcpy(outfile,*++argv);
		}
	else
		{
		strcpy(command,*++argv);
		if (command[0]=='-')
			{
			if ((!strcmp(command,"-a"))||(!strcmp(command,"-m"))||(!strcmp(command,"-h"))||(!strcmp(command,"-p"))) 
				{
				fsw=1;
				switch(command[1])
					{
					case 'p':
						fv=PVAR;
						break;
					case 'h':
						fv=HVAR;
						break;
					case 'm':
						fv=MVAR;
						break;
					case 'a':
					default:
						fv=AVAR;
					}
				strcpy(source,*++argv);
				strcpy(outfile,*++argv);
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
    nto=load_ykt(source);
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
//		printf("flag_variance=%d\n",flag_variance);
		}
	    if (flag_variance&1)
		{
		flag_slopes[0]=1;
//		printf("flag_slopes[0]=%d\n",flag_slopes[0]);
		}
	    else
		{
		flag_slopes[0]=0;
//		printf("flag_slopes[0]=%d\n",flag_slopes[0]);
		}
/*	    if (flag_variance==3)
		flag_conf=0;*/
/*	    for (i=0;i<6;++i) printf("%d ",flag_slopes[i]);
	    printf("\n");*/
	    N=serie_dev(nto, tau, dev);
	    for(i=0;i<N;++i) avar[i]=dev[i]*dev[i];
	    err=relatfit(N,tau,avar,tau,6);
/*	    printf("# Asymptote coefficients:\n");
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
	    if (flag_variance==PVAR)
	    	pvardof(N, tau, alpha, edf);
	    else
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
/*	        for(i=0;i<N;++i) printf("%d \t ",alpha[i]);
                printf("\n");*/
	    	if (flag_variance==PVAR)
	    		pvardof(N, tau, alpha, edf);
	    	else
	    		avardof(N, tau, alpha, edf);
/*		printf("# Asymptote coefficients:\n");
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
	    		printf("# Tau       \t Pdev       \t Pdev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
	    		fprintf(ofd,"# Tau       \t Pdev       \t Pdev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
			break;

		default :
	    		printf("# Tau       \t Adev       \t Adev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
	    		fprintf(ofd,"# Tau       \t Adev       \t Adev unbiased\t 2.5 %% bound \t 16 %% bound \t 84 %% bound \t 97.5 %% bound\n");
		}
	    for(i=0;i<N;++i)
	        {
		rayl=raylconfint(edf[i]);
		bmin[i]=dev[i]*sqrt(edf[i])/rayl.sup_bound;
		bmax[i]=dev[i]*sqrt(edf[i])/rayl.inf_bound;
		rayl=raylconfint1s(edf[i]);
		bi1s[i]=dev[i]*sqrt(edf[i])/rayl.sup_bound;
		bx1s[i]=dev[i]*sqrt(edf[i])/rayl.inf_bound;
		adc[i]=dev[i]/rayl.mean;
/*		if (flag_variance!=3)
			{*/
			printf("%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],dev[i],adc[i],bmin[i],bi1s[i],bx1s[i],bmax[i]);
			fprintf(ofd,"%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],dev[i],adc[i],bmin[i],bi1s[i],bx1s[i],bmax[i]);
/*			}
		else
			{
			printf("%12.6e \t %12.6e \t %12.6e\n",tau[i],dev[i],adc[i]);
			fprintf(ofd,"%12.6e \t %12.6e \t %12.6e\n",tau[i],dev[i],adc[i]);
			}*/
		}
	    fclose(ofd);
	    if (flag_bias)
	      for(i=0;i<N;++i) avar[i]=adc[i]*adc[i];
	    for(i=0;i<N;++i) w[i]=((double)1)/edf[i];
	    err=relatfit(N,tau,avar,w,6);
	    printf("# Variance asymptote coefficients:\n");
	    printf("# tau^-3 \t tau^-2   \t tau^-1 \t tau^0    \t tau^1  \t tau^2\n");
	    for(i=0;i<6;++i) printf("%12.6e \t ",coeff[i]);
	    printf("\n");
	    printf("# Deviation asymptote coefficients:\n");
	    printf("# tau^-3/2 \t tau^-1   \t tau^-1/2 \t tau^0    \t tau^1/2  \t tau^1\n");
	    for(i=0;i<6;++i) printf("%12.6e \t ",sqrt(coeff[i]));
	    printf("\n");
	    if (flag_graph)
	        {
/* Use of gnuplot for generating the graph as a ps file */
		err=gener_gplt(outfile,N,tau,dev,bmax,"unbiased");
		if (err) printf("# Error %d: ps file not created\n",err);
		}
	    }
	}
    }
		  
