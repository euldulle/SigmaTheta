/*   3ch_sbr.c                      F. Vernotte - First release: 2019/07/28 */
/* 			    						    */
/*   Subroutines for 3cornered_hat (KLTG & KLTS)			    */
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "xorshift1024star.h"
#include "ziggurat.h"
#include "filtre.h"

#define DIRECT 1
#define INVERSE 0
#define DATAMAX 10000000   /* Maximum sequence size 			    */

extern long GR, cpt_nan, *indices[3], *indix;
extern double *montecarl[3], *pmont[3], *asort[3], *cd, *expal;
double *table;
const double Pi=3.14159265358979323846;

int cmpfunc (const void *a, const void *b)
	{
    	int ia = *(int *)a;
    	int ib = *(int *)b;
    	return (table[ia] > table[ib]) - (table[ia] < table[ib]);
	}

long *matlab_sort(double *tab, long nbr, double *tsort, long *isort)
	{
	long ii;

	for(ii=0;ii<nbr;++ii) 
		indix[ii]=ii;
	table=tab;
	qsort(indix,nbr,sizeof(long),cmpfunc);
	for(ii=0;ii<nbr;++ii)
		{
		isort[ii]=indix[ii];
		tsort[ii]=tab[isort[ii]];
		}
	return(isort);
	}

double *eig(double *cov, int n, double *egval, double *egvec)
	{
	int ii, jj;
	gsl_eigen_symmv_workspace *w;
	gsl_matrix_view m;
	gsl_vector *eval;
	gsl_matrix *evec;

	w=gsl_eigen_symmv_alloc(n);
	m = gsl_matrix_view_array (cov, n, n);
	eval = gsl_vector_alloc (n);
	evec = gsl_matrix_alloc (n,n);
	gsl_eigen_symmv (&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
	
	for(ii=0;ii<n;++ii)
		{
		egval[ii]=gsl_vector_get(eval, ii);
		for(jj=0;jj<n;++jj)
			egvec[ii*n+jj]=gsl_matrix_get(evec, ii, jj);
		}
	gsl_vector_free (eval);
	gsl_matrix_free (evec);
	return(egvec);
	}

double normpdf(double x, double m, double s)
	{
	double res;

	res=((double)1)/(s*sqrt(((double)2)*Pi))*exp(-(x-m)*(x-m)/(((double)2)*s*s));
	return(res);
	}

// KLTG method
double pgaussdirectfin(double *triv, double *x0, double nmesures, double varnoise)
	{
	int ii, jj, k, corner;
	double covariances[3][3]={{0,0,0},{0,0,0},{0,0,0}};
	double vardifftheo[3], v[3][3], d[3], x[3], theo[3], p[3];
	double vara, varb, varc, probadirect;

	vara=triv[0];varb=triv[1];varc=triv[2];
	corner=1;
	if (corner)
		{
		covariances[0][0]=(vara*(((double)2)*vara+varb+varc)+varb*varc+varnoise*(((double)2)*vara+varnoise+varb+varc))/nmesures;
		covariances[1][1]=(varb*(((double)2)*varb+varc+vara)+varc*vara+varnoise*(((double)2)*varb+varnoise+vara+varc))/nmesures;
		covariances[2][2]=(varc*(((double)2)*varc+vara+varb)+vara*varb+varnoise*(((double)2)*varc+varnoise+vara+varb))/nmesures;
		covariances[0][1]=(vara*varb-varc*(vara+varb+varnoise))/nmesures;
		covariances[1][0]=covariances[0][1];
		covariances[1][2]=(varb*varc-vara*(varb+varc+varnoise))/nmesures;
		covariances[2][1]=covariances[1][2];
		covariances[2][0]=(varc*vara-varb*(varc+vara+varnoise))/nmesures;
		covariances[0][2]=covariances[2][0];
		}
	else
		{
		triv[0]=vardifftheo[0]=vara+varb;
 		triv[1]=vardifftheo[1]=varc+vara;
 		triv[2]=vardifftheo[2]=varb+varc;
		covariances[0][0]=((double)2)*(vara+varb)*(vara+varb)/nmesures;
		covariances[1][1]=((double)2)*(varc+vara)*(varc+vara)/nmesures;
		covariances[2][2]=((double)2)*(varb+varc)*(varb+varc)/nmesures;

		covariances[0][1]=((double)2)*vara*vara/nmesures;
		covariances[1][0]=covariances[0][1];
		covariances[1][2]=((double)2)*varc*varc/nmesures;
		covariances[2][1]=covariances[1][2];
		covariances[2][0]=((double)2)*varb*varb/nmesures;
		covariances[0][2]=covariances[2][0];
		}
	eig(&covariances[0][0], 3, d, &v[0][0]); 
	for(ii=0;ii<3;++ii) // Produits matriciels x=v'*x0 et theo=v'*triv'
		{
		x[ii]=theo[ii]=(double)0;
		for(jj=0;jj<3;++jj)
			{
			x[ii]+=v[jj][ii]*x0[jj];
			theo[ii]+=v[jj][ii]*triv[jj];
			}
		}

	for(k=0;k<3;++k)
		p[k]=normpdf(x[k],theo[k],sqrt(d[k]));

	probadirect=p[0]*p[1]*p[2];
	if (isnan(probadirect))
		{
		++cpt_nan;
		probadirect=0;
		}
	return(probadirect);
	}

// KLTS method
double pgaussdirectavecestimbr(double *triv, double estima, double estimb, double estimc, double var_br, int nmesures, double ddl)
	{
	int ii, jj, k, l, nk;
	double covariances[3][3]={{0,0,0},{0,0,0},{0,0,0}}, p[3]={0,0,0};
	double cov2[2][2]={{0,0},{0,0}};
	double v2[2][2],d2[2],coef2[2][2],xcar2[2];
	double v[3][3], d[3], x[3][256], coeff[3][3], xcarre[3];
	double vara, varb, varc, probadirect;
	double varab, varbc, varac;

	if (var_br==(double)0)
		nk=2;
	else
	    	nk=3;

	varab=estima+estimb+var_br;
	varbc=estimb+estimc+var_br;
	varac=estimc+estima+var_br;

	vara=triv[0];varb=triv[1];varc=triv[2];

	covariances[0][0]=vara+varb+var_br;
	cov2[0][0]=covariances[0][0];
	covariances[1][1]=vara+varc+var_br;
	cov2[1][1]=covariances[1][1];
	covariances[2][2]=varc+varb+var_br;
	covariances[0][1]=-vara;	// FV's choice (alter)
	covariances[1][0]=covariances[0][1];
	cov2[0][1]=covariances[0][1];
	cov2[1][0]=covariances[0][1];
	covariances[2][0]=-varb;
	covariances[0][2]=covariances[2][0];
	covariances[1][2]=-varc;	// FV's choice (alter)
	covariances[2][1]=covariances[1][2];

	if (var_br==(double)0)
		{
		eig(&cov2[0][0], 2, d2, &v2[0][0]); 
		for(k=0;k<2;++k)
			{
			for(ii=0;ii<2;++ii)
				for(jj=0;jj<2;++jj)
					coef2[ii][jj]=v2[ii][k]*v2[jj][k];
			xcar2[k]=coef2[0][0]*varab+coef2[1][1]*varac-((double)2)*coef2[0][1]*estima;
			p[k]=pow(d2[k],-ddl/((double)2))*exp(-xcar2[k]*ddl/((double)(2*nmesures))/d2[k]);
			}
		}
	else
		{
		eig(&covariances[0][0], 3, d, &v[0][0]); 
		for(ii=0;ii<3;++ii)
			{
			for(k=0;k<3;++k)
				for(l=0;l<3;++l)
					coeff[k][l]=v[k][ii]*v[l][ii];
			xcarre[ii]=coeff[0][0]*varab+coeff[1][1]*varac+coeff[2][2]*varbc-((double)2)*coeff[0][1]*estima-((double)2)*coeff[0][2]*estimb-((double)2)*coeff[1][2]*estimc;
			p[ii]=pow(d[ii],ddl/((double)-2))*exp(-xcarre[ii]*ddl/((double)(2*nmesures))/d[ii]);
			}
		}

	probadirect=(double)1;
	for(k=0;k<nk;++k)
		{
		probadirect*=p[k];
		}
	if (isnan(probadirect))
		{
		++cpt_nan;
		probadirect=0;
		}
	return(probadirect);
	}

double xinvy(double *x, double *y, double yexp, long nbr)
	{
	double xexp, f;
	long ii, k;

	f=(double)0;
	k=0;
	for(ii=0;ii<nbr;++ii)
		if ((y[ii]>f)&&(y[ii]<yexp))
			{
			f=y[ii];
			k=ii;
			}
	return(x[k]);
	}

int ci_kltgs(double mest[3], double varnoise, double ddl, int klt, double est[3][3], double ci[3][5])
	{
	int k, graine;
	long nmont, ii, jj, l, nb_vr, indutile;
	double mlogest, minest, maxest, thetac, B2;
	double exp_min, exp_pas, exp_max, expos;
	double mc[3], vraies[10000];

//	printf("mest: %e %e %e, varnoise: %e, ddl: %e, klt: %d\n",mest[0],mest[1],mest[2],varnoise,ddl,klt);
	mlogest=(double)0;
	for(k=0;k<3;++k)
		mlogest+=log(fabs(mest[k]));
	mlogest/=(double)3;
	mlogest=exp(mlogest);
	for(k=0;k<3;++k)
		mest[k]/=mlogest;
	varnoise/=mlogest;
	minest=maxest=fabs(mest[0]);
	for (k=1;k<3;++k)
		{
		if (fabs(mest[k])<minest)
			minest=fabs(mest[k]);
		if (maxest<fabs(mest[k]))
			maxest=fabs(mest[k]);
		}

// Calcul des post en supposant les trois lois de même variance
	if ((varnoise<minest/((double)10))&&(ddl>300)&&(minest>maxest/((double)10)))
		{
		thetac=minest/((double)10);
		B2=maxest*((double)10);
		}
	else
		{
		thetac=minest*((double)1e-5);
		B2=maxest*((double)1000);
		}
	exp_min=log(thetac)/log((double)10);
	exp_max=log(B2)/log((double)10);
	exp_pas=(exp_max-exp_min)/((double)1e4);
	ii=0;
	for(expos=exp_min;expos<exp_max;expos+=exp_pas)
		{
	  	vraies[ii]=pow((double)10,expos);
		++ii;
		}
	nb_vr=ii;
	
// Montecarlo avec des prior rand
	GR=init_rnd(graine);
	nmont=(long)DATAMAX; // nombre de tirages de Monte-Carlo
	for(ii=0;ii<nmont;++ii)
		for(jj=0;jj<3;++jj)
			montecarl[jj][ii]=indices[jj][ii]=asort[jj][ii]=(double)0;
	for(k=0;k<3;++k)
		{
		for(jj=0;jj<nmont;++jj)
			{
			indutile=(rand()*nb_vr)/RAND_MAX;
			montecarl[k][jj]=vraies[indutile];
			}
		matlab_sort(montecarl[k],nmont,asort[k],indices[k]);
		}
	if (!klt)
		{	// KLTG method
		for(ii=0;ii<nmont;++ii)
			for(jj=0;jj<3;++jj)
				pmont[jj][ii]=(double)0;

		for(l=1;l<nmont;++l)
			{
			for(k=0;k<3;++k) mc[k]=montecarl[k][l];
			pmont[0][l]=pmont[1][l]=pmont[2][l]=pgaussdirectfin(mc,mest,ddl,varnoise);
			}
		}
	else
		{	// KLTS method
		for(ii=0;ii<nmont;++ii)
			for(jj=0;jj<3;++jj)
				pmont[jj][ii]=(double)0;

		for(l=0;l<nmont;++l)
			{
			for(k=0;k<3;++k) mc[k]=montecarl[k][l];
			pmont[0][l]=pmont[1][l]=pmont[2][l]=pgaussdirectavecestimbr(mc, mest[0], mest[1], mest[2], varnoise, 1, ddl);
			}
		}

	for (k=0;k<3;++k)
		{
		est[k][0]=mest[k]*mlogest;
		est[k][1]=(double)0;
		for(ii=1;ii<nmont;++ii)
			est[k][1]+=pmont[k][indices[k][ii]]*montecarl[k][indices[k][ii]];
		cd[0]=pmont[k][indices[k][0]]; // cumsum
		for(ii=1;ii<nmont;++ii) cd[ii]=cd[ii-1]+pmont[k][indices[k][ii]];
		est[k][1]/=cd[nmont-1];
		est[k][1]*=mlogest;
		cd[0]=pmont[k][indices[k][0]]; // cumsum
		for(ii=1;ii<nmont;++ii) cd[ii]=cd[ii-1]+pmont[k][indices[k][ii]];
		for(ii=0;ii<nmont;++ii) cd[ii]/=cd[nmont-1];
		ci[k][0]=xinvy(asort[k],cd,(double)0.023,nmont)*mlogest;
		ci[k][1]=xinvy(asort[k],cd,(double)0.159,nmont)*mlogest;
		est[k][2]=xinvy(asort[k],cd,(double)0.5,nmont)*mlogest;
		ci[k][2]=xinvy(asort[k],cd,(double)0.841,nmont)*mlogest;
		ci[k][3]=xinvy(asort[k],cd,(double)0.955,nmont)*mlogest;
		ci[k][4]=xinvy(asort[k],cd,(double)0.977,nmont)*mlogest;
		}
	return(1);
	}


