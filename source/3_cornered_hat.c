/*   3_cornered_hat.c                          F. Vernotte - 2019/07/28     */
/*   Computation of the Groslambert covariances of a 3-cornered hat 	    */
/*   system. The confidence intervals are estimated by using the KLTG	    */
/*   and KLTS methods.							    */
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
#include <unistd.h>
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
    printf("Usage: 3CorneredHat [-g|c] [-a|f] SOURCE12 SOURCE23 SOURCE32 TARGET\n\n");
    printf("Computes the Allan deviation of each clock of a 3-cornered hat system and estimates the confidence intervals.\n\n");
    printf("The input files SOURCEij contains a 2-column table with time values (dates) in the first column and normalized frequency deviation measurements of clock 'j' compared to clock 'i' in the second column. All input files must have the same dates.\n\n");
    printf("The root TARGET is used to build the 3 output files for clocks 1, 2 and 3: TARGET1.gcod, TARGET2.gcod and TARGET3.gcod.\n\n"); 
    printf("Each of these files is a 7-column table with:\n");
    printf("\t1st column: tau values\n");
    printf("\t2nd column: Groslambert codeviation estimate\n");
    printf("\t3rd column: mean or 50 %% estimate\n");
    printf("\t4th column: 2.5 %% bound\n");
    printf("\t5th column: 16 %% bound\n");
    printf("\t6th column: 84 %% bound\n");
    printf("\t7th column: 97.5 %% bound.\n");
    printf("The configuration file \".SigmaTheta.conf\" is taken into account (e.g. choice of the increment of the tau values).\n\n");
    printf("The file TARGET.gnu is generated for invoking gnuplot. \n");
    printf("The file TARGET.pdf is the pdf file of the gnuplot graph (if the PDF option has been chosen in the configuration file).\n\n");
    printf("If the option '-g' is selected, the Groslambert covariance method is invoked (default). \n");
    printf("If the option '-c' is selected, the classical 3-cornered hat method is invoked. \n");
    printf("If the option '-m' is selected, the mean estimate is displayed in the 3rd column of the output files (default).\n");
    printf("If the option '-f' is selected, the 50 %% estimate (median) is displayed in the 3rd column of the output files.\n");
    printf("If the option '-d' is selected, a drift is removed from all SOURCEij data.\n");
    printf("If the option '-i' is selected, 3 individual graphs are displayed (one for each clock).\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
	{
	char source1[256], source2[256], source3[256];
	char output[256], target[4][256], chestim[256];
	int err, rep, gflag, cflag, aflag, fflag;
	int dflag, iflag, ind3ch, ind50;
	int i, j, k, N, nto, alpha[256], klt;
	long nmont;
	double tau[256], measdev[256], varnoise[256];
	double adev12[256], adev23[256], adev31[256];
	double avar[256], edfm[256], tsas, asympt;
	double edf12[256], edf23[256], edf31[256];
	double mesti[3], est[3][3], ci[3][5];
	double gadev[4][256], bmin[3][256], bmax[3][256];
	FILE *tarf1, *tarf2, *tarf3;

	err=gflag=cflag=aflag=fflag=dflag=iflag=ind3ch=ind50=0;
	strcpy(chestim,"mean  ");
	while ((rep = getopt (argc, argv, "gcmfdi")) != -1)
		switch(rep)
			{
			case 'g':
				gflag = 1;
				break;
			case 'c':
				cflag = ind3ch = 1;
				break;
			case 'm':
				aflag = 1;
				break;
			case 'f':
				fflag = ind50 = 1;
				strcpy(chestim,"median");
				break;
			case 'd':
				dflag = 1;
				break;
			case 'i':
				iflag = 1;
				break;
			case '?':
				printf("# Unknown option '-%c'\n",optopt);
				usage();
			default:
				exit(-1);
			}
	if (gflag&cflag)
		{
		printf("# Incompatible options '-g' and '-c'\n");
		usage();
		exit(-1);
		}
	if (aflag&fflag)
		{
		printf("# Incompatible options '-m' and '-f'\n");
		usage();
		exit(-1);
		}
	if (argc-optind<4)
		{
		printf("# Missing arguments\n");
		usage();
		exit(-1);
		}
	if (argc-optind>4)
		{
		printf("# Too many arguments\n");
		usage();
		exit(-1);
		}

	strcpy(source1,argv[optind]);
	strcpy(source2,argv[optind+1]);
	strcpy(source3,argv[optind+2]);
	strcpy(output,argv[optind+3]);

	for (i=0;i<4;++i) strcpy(target[i],output);
	if (ind3ch)
		{
		strcat(target[0],"_1.3ch");
		strcat(target[1],"_2.3ch");
		strcat(target[2],"_3.3ch");
		}
	else
		{
		strcat(target[0],"_1.gcod");
		strcat(target[1],"_2.gcod");
		strcat(target[2],"_3.gcod");
		}
	strcat(target[3],"_cls.adev");

	err=init_flag();
	if (err==-1) printf("# ~/.SigmaTheta.conf not found, default values selected\n");
	if (err==-2) printf("# ~/.SigmaTheta.conf improper, default values selected\n");
	if (err)
		{
	        flag_log_inc=1;
	        log_inc=(double)2;
		}

	N=load_3yk(source1,source2,source3);
	if (N<1)
		{
		switch(N)
			{
			case 0:
				printf("# Empty file\n");
				break;
			case -1:
				printf("# File not found\n");
				break;
			case -2:
				printf("# Different file lengths\n");
			}
		exit(N);
		}
	Y=(double *)malloc(N*sizeof(double)); 
// Optional drift removal
	if (dflag)
		{
		printf("# Drift coefficients (y=a.t+b)\n");
		// Y12
		for(i=0;i<N;++i) Y[i]=Y12[i];
		tchebyfit((long)N);
		for(i=0;i<N;++i) Y12[i]=Y[i];
		printf("# %s: a=%24.16e \t b=%24.16e\n",source1,coeff[1],coeff[0]);
		// Y23
		for(i=0;i<N;++i) Y[i]=Y23[i];
		tchebyfit((long)N);
		for(i=0;i<N;++i) Y23[i]=Y[i];
		printf("# %s: a=%24.16e \t b=%24.16e\n",source2,coeff[1],coeff[0]);
		// Y31
		for(i=0;i<N;++i) Y[i]=Y31[i];
		tchebyfit((long)N);
		for(i=0;i<N;++i) Y31[i]=Y[i];
		printf("# %s: a=%24.16e \t b=%24.16e\n\n",source3,coeff[1],coeff[0]);
		}

// ADev of the measurement noise
	for (i=0;i<N;++i)
		Y[i]=(Y12[i]+Y23[i]+Y31[i])/sqrt((double)3);
	flag_variance=0;
	nto=serie_dev(N, tau, gadev[3]);
	tarf1=fopen(target[3],"w");
	for (i=0;i<nto;++i) 
		{
		fprintf(tarf1,"%12.6e %12.6e\n",tau[i],gadev[3][i]);
		varnoise[i]=gadev[3][i]*gadev[3][i];
		}
	fclose(tarf1);

// ADev of each file
	for (i=0;i<N;++i)
		Y[i]=Y12[i];
	nto=serie_dev(N, tau, adev12);
	for (i=0;i<N;++i)
		Y[i]=Y23[i];
	nto=serie_dev(N, tau, adev23);
	for (i=0;i<N;++i)
		Y[i]=Y31[i];
	nto=serie_dev(N, tau, adev31);

// GCoDev of each file pair or 3-cornered hat
	if (ind3ch)
	// Classical 3-cornered hat method
		for(i=0;i<nto;++i)
			{
			/* gcod are multiplied by -1 for the sake of compatibility 
			with the Groslambert covariance*/
			gadev[0][i]=-(+adev12[i]-adev23[i]+adev31[i])/((double)2);
			gadev[1][i]=-(+adev12[i]+adev23[i]-adev31[i])/((double)2);
			gadev[2][i]=-(-adev12[i]+adev23[i]+adev31[i])/((double)2);
			}
	else
	// Groslambert covariance method
		{
		Y1=(double *)malloc(N*sizeof(double)); 
		Y2=(double *)malloc(N*sizeof(double)); 
		for (i=0;i<N;++i)
			{
			Y1[i]=Y12[i];
			Y2[i]=Y23[i];
			}
		flag_variance=4;
		nto=serie_dev(N, tau, gadev[1]);
		for (i=0;i<N;++i)
			Y1[i]=Y31[i];
		nto=serie_dev(N, tau, gadev[2]);
		for (i=0;i<N;++i)
			Y2[i]=Y12[i];
		nto=serie_dev(N, tau, gadev[0]);
		}
	// Release of the memory allocated to the arrays Y, Y1, Y2, Y12, Y23 and Y31
	free(Y);
	free(Y1);
	free(Y2);
	free(Y12);
	free(Y23);
	free(Y31);

// Estimation of the EDF
	flag_variance=0;
    // File 12
	for(i=0;i<nto;++i) avar[i]=adev12[i]*adev12[i];
	err=relatfit(nto,tau,avar,tau,6);
	for(i=0;i<nto;++i)
		{
		asympt=0;
		for(j=1;j<5;++j)
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
	avardof(nto, tau, alpha, edf12);
    // File 23
	for(i=0;i<nto;++i) avar[i]=adev23[i]*adev23[i];
	err=relatfit(nto,tau,avar,tau,6);
	for(i=0;i<nto;++i)
		{
		asympt=0;
		for(j=1;j<5;++j)
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
	avardof(nto, tau, alpha, edf23);
    // File 31
	for(i=0;i<nto;++i) avar[i]=adev31[i]*adev31[i];
	err=relatfit(nto,tau,avar,tau,6);
	for(i=0;i<nto;++i)
		{
		asympt=0;
		for(j=1;j<5;++j)
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
	avardof(nto, tau, alpha, edf31);
    // EDF mean
	for(i=0;i<nto;++i)
		edfm[i]=(edf12[i]+edf23[i]+edf31[i])/((double)3);

// Computation of the GCod confidence intervals
    // Initialization
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
    // Files TARGET1, 2 and 3
	tarf1=fopen(target[0],"w");
	tarf2=fopen(target[1],"w");
	tarf3=fopen(target[2],"w");

	printf("#  Tau \t B 2.5%%   < Clock 1  < B97.5%%   \t B 2.5%%   < Clock 2  < B97.5%%   \t B 2.5%%   < Clock 3  < B97.5%%\n");
	fprintf(tarf1,"# Tau        \t GCod         \t %s      \t B 2.5 %%     \t B 16 %%       \t B 84 %%       \t B 97.5 %%\n",chestim);
	fprintf(tarf2,"# Tau        \t GCod         \t %s      \t B 2.5 %%     \t B 16 %%       \t B 84 %%       \t B 97.5 %%\n",chestim);
	fprintf(tarf3,"# Tau        \t GCod         \t %s      \t B 2.5 %%     \t B 16 %%       \t B 84 %%       \t B 97.5 %%\n",chestim);
	for (i=0;i<nto;++i)
		{
		mesti[0]=-gadev[0][i];
		mesti[1]=-gadev[1][i];
		mesti[2]=-gadev[2][i];
		for(j=0;j<3;++j)
			if (mesti[j]<0)
				mesti[j]*=-mesti[j];
			else
				mesti[j]*=mesti[j];
		if (edfm[i]<100)
			klt=1;
		else
			klt=0;

		err=ci_kltgs(mesti, varnoise[i], edfm[i], klt, est, ci);
		for(j=0;j<3;++j)
			{
			for(k=0;k<3;++k)
				if (est[j][k]<0)
				    est[j][k]=-sqrt(-est[j][k]);
				else
				    est[j][k]=sqrt(est[j][k]);
			for(k=0;k<5;++k)
				ci[j][k]=sqrt(ci[j][k]);
			bmin[j][i]=est[j][0];
			if (est[j][1+ind50]<bmin[j][i]) bmin[j][i]=est[j][1+ind50];
			bmax[j][i]=ci[j][4];
			}

		printf("%6.0f \t %8.2e < %8.2e < %8.2e \t %8.2e < %8.2e < %8.2e \t %8.2e < %8.2e < %8.2e\n",tau[i],ci[0][0],est[0][1+ind50],ci[0][4],ci[1][0],est[1][1+ind50],ci[1][4],ci[2][0],est[2][1+ind50],ci[2][4]);
		fflush(stdout);
		fprintf(tarf1,"%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],est[0][0],est[0][1+ind50],ci[0][0],ci[0][1],ci[0][2],ci[0][4]);
		fflush(tarf1);
		fprintf(tarf2,"%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],est[1][0],est[1][1+ind50],ci[1][0],ci[1][1],ci[1][2],ci[1][4]);
		fflush(tarf2);
		fprintf(tarf3,"%12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n",tau[i],est[2][0],est[2][1+ind50],ci[2][0],ci[2][1],ci[2][2],ci[2][4]);
		fflush(tarf3);
		}
	fclose(tarf1);
	fclose(tarf2);
	fclose(tarf3);

	// Plot the graph(s)
	if (iflag)
		{
		// Clock 1
		for(i=0;i<nto;++i) avar[i]=gadev[0][i]*gadev[0][i];
		err=relatfit(nto,tau,avar,tau,6);
		err=gener_gplt(target[0],nto,tau,bmin[0],bmax[0],chestim,1);
		if (err) printf("# Error %d: pdf file not created\n",err);
		// Clock 2
		for(i=0;i<nto;++i) avar[i]=gadev[1][i]*gadev[1][i];
		err=relatfit(nto,tau,avar,tau,6);
		err=gener_gplt(target[1],nto,tau,bmin[1],bmax[1],chestim,1);
		if (err) printf("# Error %d: pdf file not created\n",err);
		// Clock 3
		for(i=0;i<nto;++i) avar[i]=gadev[2][i]*gadev[2][i];
		err=relatfit(nto,tau,avar,tau,6);
		err=gener_gplt(target[2],nto,tau,bmin[2],bmax[2],chestim,1);
		if (err) printf("# Error %d: pdf file not created\n",err);
		}
	else
		{
		err=gen_3chplt(target, output, nto, tau, bmin, bmax, ind50);
		if (err) printf("# Error %d: pdf file not created\n",err);
		}

	// Release of the memory
	for(k=0;k<3;++k)
		{
		free(montecarl[k]);
		free(pmont[k]);
		free(asort[k]);
		free(indices[k]);
		}
	free(cd);
	free(indix);
	free(expal);
	exit(err);		
	}

