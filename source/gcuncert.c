/*   gcuncert.c                                F. Vernotte - 2019/07/28     */
/*   Computation of confidence intervals from Groslambert covariance        */
/*   estimates by using the KLTG and KLTS methods			    */
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
#include "filtre.h"
#include "3ch_sbr.h"
#ifdef DATAMAX
#undef DATAMAX
#define DATAMAX 10000000   /* Maximum sequence size 			    */
#endif
#define NUMAX 1024    /* Maximum EDF numer */

long GR, cpt_nan;
double *montecarl[3], *pmont[3], *asort[3], *cd, *expal;
long *indices[3], *indix;
double *x;
double hm3,hm2,hm1,h0,hp1,hp2,C1,C0;
double tau0;

void usage(void)
/* Help message */
    {
    printf("Usage: GCUncert SOURCE1 SOURCE2 SOURCE3 NOISE EDF TARGET1 TARGET2 TARGET3\n\n");
    printf("Computes the 95 %% confidence interval from the GCov measurements contained in the 2 column files SOURCE1-3, from the measurement noise adev contained in the 2 column file NOISE and from the EDF numbers contained in the 2 column file EDF. \n\n");
    printf("The results are stored in the 7-column files TARGET1-3. The 7 columns contain: tau, GCov or classical 3-conered hat, mean or median estimate, 68 and 95 %% confidence interval bounds.\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
    {
    int i, j, k, nbv, N, N2, N3, N4, N5, indeq, err, klt;
    long nmont;
    char source1[256], source2[256], source3[256], noisefile[256], edffile[256];
    char target1[256], target2[256], target3[256];
    double tau[256], tau2[256], tau3[256], tau4[256], tau5[256];
    double gcod1[256], gcod2[256], gcod3[256], mesti[3], varnoise[256], edf[256];
    double est[3][3], ci[3][5];
    FILE *tarf1, *tarf2, *tarf3;

    if (argc<9)
        {
	usage();
	exit(-1);
        }
    else
        {
	strcpy(source1,*++argv);
	strcpy(source2,*++argv);
	strcpy(source3,*++argv);
	strcpy(noisefile,*++argv);
	strcpy(edffile,*++argv);
	strcpy(target1,*++argv);
	strcpy(target2,*++argv);
	strcpy(target3,*++argv);
	}
    err=0;

// File SOURCE1
    N=load_adev(source1,tau,gcod1);
    if (N==-1)
	{
        printf("# File %s not found\n",source1);
	exit(-1);
	}

// File SOURCE2
    N2=load_adev(source2,tau2,gcod2);
    if (N2==-1)
	{
        printf("# File %s not found\n",source2);
	exit(-2);
	}
    if (N!=N2)
	{
	printf("# Different number of lines in %s and %s\n",source1,source2);
	exit(-3);
	}
    indeq=0;
    for(i=0;i<N;++i)
	if (tau[i]!=tau2[i])
	    {
	    indeq=1;
	    break;
	    }
    if (indeq)
	{
	printf("# Incompatible 1st column in %s and %s\n",source1,source2);
	exit(-4);
	}

// File SOURCE3
    N3=load_adev(source3,tau3,gcod3);
    if (N3==-1)
	{
        printf("# File %s not found\n",source3);
	exit(-5);
	}
    if (N!=N3)
	{
	printf("# Different number of lines in %s and %s\n",source1,source3);
	exit(-6);
	}
    indeq=0;
    for(i=0;i<N;++i)
	if (tau[i]!=tau3[i])
	    {
	    indeq=1;
	    break;
	    }
    if (indeq)
	{
	printf("# Incompatible 1st column in %s and %s\n",source1,source3);
	exit(-7);
	}

// File NOISE
    N4=load_adev(noisefile,tau4,varnoise);
    if (N4==-1)
	{
        printf("# File %s not found\n",noisefile);
	exit(-8);
	}
    if (N!=N4)
	{
	printf("# Different number of lines in %s and %s\n",source1,noisefile);
	exit(-9);
	}
    indeq=0;
    for(i=0;i<N;++i)
	if (tau[i]!=tau4[i])
	    {
	    indeq=1;
	    break;
	    }
    if (indeq)
	{
	printf("# Incompatible 1st column in %s and %s\n",source1,noisefile);
	exit(-10);
	}
    for(i=0;i<N4;++i) varnoise[i]*=varnoise[i];

// File EDF
    N5=load_adev(edffile,tau5,edf);
    if (N5==-1)
	{
        printf("# File %s not found\n",edffile);
	exit(-11);
	}
    if (N!=N5)
	{
	printf("# Different number of lines in %s and %s\n",source1,edffile);
	exit(-12);
	}
    indeq=0;
    for(i=0;i<N;++i)
	if (tau[i]!=tau5[i])
	    {
	    indeq=1;
	    break;
	    }
    if (indeq)
	{
	printf("# Incompatible 1st column in %s and %s\n",source1,edffile);
	exit(-13);
	}

// Files TARGET1, 2 and 3
    tarf1=fopen(target1,"w");
    tarf2=fopen(target2,"w");
    tarf3=fopen(target3,"w");

// Computation of the GCod confidence intervals
    nmont=(long)DATAMAX;
    for(k=0;k<3;++k)
	{
	montecarl[k] = (double *)malloc(nmont * sizeof(double));
	pmont[k] = (double *)malloc(nmont * sizeof(double));
	asort[k] = (double *)malloc(nmont * sizeof(double));
	indices[k] = (long *)malloc(nmont * sizeof(long));
	}
    cd = (double *)malloc(nmont * sizeof(double));
    indix = (long *)malloc(nmont * sizeof(long));
    expal = (double *)malloc(nmont * sizeof(double));
    for (i=0;i<N;++i)
	{
	mesti[0]=-gcod1[i];
	mesti[1]=-gcod2[i];
	mesti[2]=-gcod3[i];
	if (edf[i]<100)
		klt=1;
	else
		klt=0;

	for(j=0;j<3;++j)
		if (mesti[j]<0)
			mesti[j]*=-mesti[j];
		else
			mesti[j]*=mesti[j];
	err=ci_kltgs(mesti, varnoise[i], edf[i], klt, est, ci);
	for(j=0;j<3;++j)
		{
		for(k=0;k<3;++k)
			if (est[j][k]<0)
			    est[j][k]=-sqrt(-est[j][k]);
			else
			    est[j][k]=sqrt(est[j][k]);
//		printf("est[%d]=%e %e %e\n",j,est[j][0],est[j][1],est[j][2]);
		for(k=0;k<5;++k)
			ci[j][k]=sqrt(ci[j][k]);
//		printf("ci[%d]=%e %e %e %e %e\n",j,ci[j][0],ci[j][1],ci[j][2],ci[j][3],ci[j][3]);
		}

	printf("%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %.3e \t %d\n",tau[i],est[0][0],est[0][1],ci[0][0],ci[0][1],ci[0][2],ci[0][4],edf[i],klt);
	fflush(stdout);
	fprintf(tarf1,"%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],est[0][0],est[0][1],ci[0][0],ci[0][1],ci[0][2],ci[0][4]);
	fflush(tarf1);
	fprintf(tarf2,"%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],est[1][0],est[1][1],ci[1][0],ci[1][1],ci[1][2],ci[1][4]);
	fflush(tarf2);
	fprintf(tarf3,"%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],est[2][0],est[2][1],ci[2][0],ci[2][1],ci[2][2],ci[2][4]);
	fflush(tarf3);
	}
    fclose(tarf1);
    fclose(tarf2);
    fclose(tarf3);
    exit(err);		
    }
    
