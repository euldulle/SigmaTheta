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

#define DIRECT 1
#define INVERSE 0

extern double *x, *y;
extern double hm3,hm2,hm1,h0,hp1,hp2,C1,C0;
extern double tau0, vm_gauss, sigma_gauss;
extern long GR;

/*	   This subroutine performs the Fourier transform of   		*/
/*	   the double tables x[] and y[] containing 			*/
/*	   respectively the real and the imaginary parts of 		*/
/*	   the input signal. The result is sent to x[] (real		*/
/*	   part) and y[] (imaginary part).				*/
/*	                                                     		*/
/*	   The data number n and the type (DIRECT or INVERSE)		*/
/*	   are given by the input parameters. The output is 		*/
/*	   the number of computed frequencies.				*/
/*	                                                     		*/
/*	F.VERNOTTE                                         1987/12/04	*/
/*	                                   modified by FV, 1989/02/21	*/

long fft(long n,int type)
	{
	double tempr,tempi,wr,wi,pi,theta,signe;
	unsigned long nn,m,i,fn,j,mmax,istep;

/*	FFT directe ou FFT inverse					*/
	
	if (type)
	    signe=(double)-1;	/* direct FFT	*/
	else
	    signe=(double)1;	/* inverse FFT	*/

/*	   Is n a power of 2 ?		                  		*/

	tempr = log((double)n)/log((double)2);
	j = (long) (tempr + .5);
	if (tempr != j)	
		n =(long)(pow((double)2,(double)j)+.1);
	pi = ((double)4)*atan((double)1);
	fn = n;

/*	   Binary inversion of the elements of table X     		*/

	j = 0;	
	for (i=0; i<n; ++i)
		{
		if (i < j)
			{
			tempr = x[j] ;
			tempi = y[j] ;
			x[j]  = x[i] ;
			y[j]  = y[i] ;
			x[i]  = tempr;
			y[i]  = tempi;
			}
		m = n/2.;
		while (j >= m)
			{
			j = j-m;
			m = (m+1)/2.;
			}
		j = j+m;
		}

/*	   COOLEY - TUKEY algorithm                     		*/

	mmax = 1;
	while (mmax < n)
		{
		istep = 2*mmax;
		for (m=0; m<mmax; ++m)
			{
			theta = signe*pi*((double)m)/((double)mmax);
			wr = cos(theta);
			wi = sin(theta);
			for (i=m; i<n; i+=istep)
                        	{
                        	j = i+mmax;
        	               	tempr = wr*x[j]-wi*y[j];
                        	tempi = wi*x[j]+wr*y[j];
                        	x[j]  = x[i]-tempr;
                        	y[j]  = y[i]-tempi;
                        	x[i]  = x[i]+tempr;
                        	y[i]  = y[i]+tempi;
                        	}
			}
        	mmax = istep;
        	}

/*	   Normalization                                       		*/

	for (i=0; i<fn; ++i)
	       {
	       x[i]/=sqrt((double)fn);
	       y[i]/=sqrt((double)fn);
	       }
	return (fn);
	}

/*	Subroutines for generating random noises.			     */
/*						FV	1989/02/21	     */

/* Initialization of the random number generator:			     */
/* The file "$HOME/.randinit2" contains the last used 'seed'.		     */
/* It is incremented and stored then the generator is initialized with the   */ 
/* new seed (=former+1).						     */
long init_rnd(long graine)
	{
	long x0,x1;
	FILE *rifich;
	char *homepath;
	char filepath[255];

	homepath=getenv("HOME");
	if (homepath!=NULL)
	{
		strncpy(filepath, homepath, sizeof(filepath));
		strncat(filepath, "/.randinit2",sizeof(filepath));
	}
	else
		strncpy(filepath,".randinit2",sizeof(filepath)); /* If the environment variable $HOME   */
/* is not set, the .randinit2 file is searched in the current directory.	      */
	if (graine) x0=graine; /* If 'graine' not equal to 0, it will be used for     */
/* initializing the generator rather than using the value contained in the file	      */
/* "$HOME/exe/.randinit2".                                                            */
	else
		{
		rifich=fopen(filepath,"r");
		if (rifich==NULL)
			{
			printf("File %s not found.\n",homepath);
			exit(-1);
			}
		fscanf(rifich,"%ld",&x0);
		fclose(rifich);
		x1=x0+1;
		if (x1>RAND_MAX-2) x1=2;
		rifich=fopen(filepath,"w");
		fprintf(rifich,"%ld\n",x1);
		fclose(rifich);
		}
	srand(x0);
	return(x0);
        }


/*	Generation of a Gaussian random number by averaging N uniformly     */
/* distributed random numbers (application of the central limit theorem).   */
double gauss(N)
int N;
	{
	int i;
	double x;


	x=(double)0;
	for (i=0;i<N;++i)
		x+=(double)rand();
	x=(x-vm_gauss)/sigma_gauss;
	return(x);
	}		
		

/*  Generation of a sequence of 'nbr_dat' Gaussian random numbers, centered */
/*  and with unity rms.                                                     */
double gausseq(long nbr_dat,int graine)
	{
	long int i,ordre;
	double xm,x2m,var,sig;

	GR=init_rnd(graine);
	ordre=16;
	vm_gauss=((double)ordre)*((double)RAND_MAX)/((double)2); /* Average of 'ordre' numbers between 0 and 'RAND_MAX' */
	sigma_gauss=sqrt( ((double)ordre) * ((double)RAND_MAX) * ((double)(RAND_MAX-2)) / ((double)6) ); /* RMS of 'ordre' numbers between 0 and 'RAND_MAX' */
	for(i=0;i<nbr_dat;++i)
	    {
	    x[i]=(double)gauss(ordre);
	    y[i]=(double)0;
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
	double x2m,Ri,R2i,Rx,pi,cor,ksx,rho,phi;

	limite=nbr_dat/2;
	pi=((double)4)*atan((double)1); /* Computation of pi.		    */
	if ((hm3)||(hm2)||(hm1)||(hp1)||(hp2)) /* If only white noise is    */
/* needed, we don't have to filter the sequence!                            */
		{
		fft(nbr_dat,DIRECT);
		cor=sqrt(h0/tau0);
		x[0]*=cor; /* Filtering of the null frequency amplitude:              */
		y[0]*=cor; /* only h0 is used for avoiding to multiply or divide this */
/* amplitude by 0 with respectively positive and negative exponents of the power      */
/* laws. This amplitude corresponds to the mean of the sequence and should not be set */
/* to 0.									      */
		for(i=1;i<limite;++i)
			{
			cor=((double)i)/((double)nbr_dat);
			Ri=cor/tau0;
			R2i=Ri*Ri;
			Rx=sqrt(fabs((double)(hm3/R2i/Ri + hm2/R2i + hm1/Ri + h0 + hp1*Ri + hp2*R2i))/tau0);
			x[i]*=Rx;           /* "positive" frequencies */
			y[i]*=Rx;
			x[nbr_dat-i]*=Rx;   /* "negative" frequencies */
			y[nbr_dat-i]*=Rx;
			}
		Ri=((double).5)/tau0;
		R2i=Ri*Ri;
		Rx=sqrt(fabs((double)(hm3/R2i/Ri + hm2/R2i + hm1/Ri + h0 + hp1*Ri + hp2*R2i))/tau0);
		x[limite]*=Rx;
		y[limite]*=Rx;
		fft(nbr_dat,INVERSE); /* After filtering, the inverse FFT is computed.*/
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








