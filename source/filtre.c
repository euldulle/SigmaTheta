/*   filtre.c                       F. Vernotte - First release: 1994/02/02 */
/*					    Sigma-Theta version: 2015/10/26 */
/*   Filtering subroutines of bruiteur 					    */
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
#include "xorshift1024star.h"
#include "ziggurat.h"

#define DIRECT 1
#define INVERSE 0

extern double *x;
fftw_complex *freq_x;
extern double hm3,hm2,hm1,h0,hp1,hp2,C1,C0;
extern double tau0;
extern long GR;

/*	Subroutines for generating random noises.			     */
/*						FV	1989/02/21	     */

/* Initialization of the random number generator using /dev/urandom	     */
long init_rnd(long graine)
	{
	long x0,x1;
	FILE *rifich;
	char *homepath;
	char filepath[255];

	strncpy(filepath, "/dev/urandom",sizeof(filepath));
	
	rifich=fopen(filepath,"r");
	if (rifich==NULL)
		{
		printf("File %s not found.\n",filepath);
		exit(-1);
		}
	fread(&x0, sizeof(x0), 1, rifich);
	fclose(rifich);
	x1=x0+1;

	xorshift1024_init64(x0);
	return(x0);
        }



/*  Generation of a sequence of 'nbr_dat' Gaussian random numbers, centered */
/*  and with unity rms.                                                     */
double gausseq(long nbr_dat,int graine)
	{
	long int i,ordre;
	double xm,x2m,var,sig;

	GR=init_rnd(graine);
	for(i=0;i<nbr_dat;++i)
	    {
	    x[i]=ojr_next_normal();
	    }
/*  We check that the RMS is one:                         		    */
	xm=x2m=(double)0;
	for(i=0;i<nbr_dat;++i)
	    {
	    xm+=(double)x[i];
	    x2m+=((double)x[i])*((double)x[i]);
	    }
	xm/=(double)nbr_dat;
	x2m/=(double)nbr_dat;
	var=x2m-xm*xm;
	sig=sqrt(((double)2)*var);
/* Theoretically, 'sig' should be a random variable centered around 1.      */
/* It is the classical estimator of the RMS of the computed sequence.       */
	return(sig);
	}

/*  Subroutine filtering the white noise sequence of 'nbr_dat' Gaussian     */
/*  data sampled with a step 'tau0' in order to obtain a sequence with      */
/*  a linear frequency drift and a Power Spectral Density following power   */
/*  laws (from f^{-2} to f^{+2}) according to the entered coefficients.     */
/*									    */
/*							    1994/02/02	    */
/*									    */
/*						    Francois Vernotte	    */
/*									    */

void filtreur(long nbr_dat, double tau0)
	{
	long i,limite,saut,finsaut;
	double x2m,Ri,R2i,Rx,cor,ksx,rho,phi, one_over_nbr;
	fftw_plan p_forward, p_backward;

	limite=nbr_dat/2;
	one_over_nbr = 1/((double)nbr_dat);

	if ((hm3)||(hm2)||(hm1)||(hp1)||(hp2)) /* If only white noise is    */
/* needed, we don't have to filter the sequence!                            */
		{
		freq_x = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * limite);
		p_forward = fftw_plan_dft_r2c_1d(nbr_dat, x, freq_x, FFTW_ESTIMATE);
		p_backward = fftw_plan_dft_c2r_1d(nbr_dat, freq_x, x, FFTW_ESTIMATE);
		fftw_execute(p_forward);

		cor=sqrt(h0/tau0);
		freq_x[0] *= cor * one_over_nbr; /* Filtering of the null frequency amplitude:              */
		/* only h0 is used for avoiding to multiply or divide this */
/* amplitude by 0 with respectively positive and negative exponents of the power      */
/* laws. This amplitude corresponds to the mean of the sequence and should not be set */
/* to 0.									      */
		for(i=1;i<limite;++i)
			{
			cor=((double)i) * one_over_nbr;
			Ri=cor/tau0;
			R2i=Ri*Ri;
			Rx=sqrt(fabs((double)(hm3/R2i/Ri + hm2/R2i + hm1/Ri + h0 + hp1*Ri + hp2*R2i))/tau0) * one_over_nbr;
			freq_x[i]*=Rx;           /* "positive" frequencies */
			}
		Ri=((double).5)/tau0;
		R2i=Ri*Ri;
		Rx = sqrt(fabs((double)(hm3/R2i/Ri + hm2/R2i + hm1/Ri + h0 + hp1*Ri + hp2*R2i))/tau0) * one_over_nbr;
		freq_x[limite] *= Rx;

		fftw_execute(p_backward); /* After filtering, the inverse FFT is computed.*/

		fftw_destroy_plan(p_forward);
		fftw_destroy_plan(p_backward);
		fftw_free(freq_x);
		}
	else
		if ((h0!=1)||(tau0!=1))
			{
			cor=sqrt(h0/tau0); /* Noise level of the pure white noise.    */
			for(i=0;i<nbr_dat;++i) x[i]*=cor;
			}
	if ((C0)||(C1)) /* Then, the drift is added in the time domain.		      */
		{
		for(i=0;i<nbr_dat;++i)
		    {
		    Ri=((double)i)*tau0;
		    x[i]+= ((double)C0) + ((double)C1)*Ri;
		    }
		}
	}

/* Conversion of a frequency deviation Yk sequence into a time error x(t) sequence.   */
/* The arbitrarily convention x(t0) = Y0 is used.                                     */
/* At the beginning the YK are in x[], at the end the x(t) are also in x[].           */
void yk_xt(long nbr_dat, double x0, double tau0)
    {
    long i;
    double xint,yint;

    yint=x[0];
    x[0]=x0;

    for(i=0;i<nbr_dat-1;++i)
	{
	xint=((double)x[i])+yint*tau0;
	yint=x[i+1];
	x[i+1]=xint;
	}
    }


/*  Adding the suffix to the file name.					    */
/*							    FV	1992/04/17  */

char *mod_nom(nom_brut,suffixe)
char *nom_brut,*suffixe;
	{
	int i,point;
	static char nomfich[64];

	for (i=0;i<64;++i)         /* Seeking the '.' in the file name. */
		if (nom_brut[i]==0)
			{
			point=i;
			break;
			}
	for(i=0;i<=point;++i) nomfich[i]=nom_brut[i];
	i=0;
	do
		++i;
	while((nomfich[i]!=0)&&(nomfich[i]!='.'));
	if (nomfich[i]==0)
		nomfich[i]='.';
	strcpy(&nomfich[i+1],suffixe); /* Adding the suffix. */

	return(nomfich);
	}








