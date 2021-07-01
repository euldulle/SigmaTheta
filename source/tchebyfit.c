/*   tchebyfit.c                                  F. Vernotte - 2000/06/08  */
/*   Linear fit subroutines using the first 2 Tchebytchev polynomials       */
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

#define lg(x) ((long)(x))
#define db(x) ((double)(x))

extern double *T, *Y;
extern double coeff[10];

double (* phi[2]) (long,long);

double tcheby0(long k, long N)
    {
    double result;
    result=db(1)/sqrt(db(N));
    return(result);
    }

double tcheby1(long k, long N)
    {
    double result;

    result=db(2)*db(k) - (db(N) - db(1));
    result*=sqrt( db(3)/((db(N)-db(1))*db(N)*(db(N)+db(1))) );
    return(result);
    }

void initcheb(void)
    {
    phi[0]=tcheby0;
    phi[1]=tcheby1;
    }

double extrapol(long k, long N)
    {
    int i;
    double result;

    result=db(0);
    for(i=0;i<2;++i)
        result+=coeff[i]*phi[i](k,N);
    return(result);
    }

double tchebyfit(long int N)
    {
    int i, j;
    double moy,var,coco,dx,dy,a,b;

    initcheb();
    for(i=0;i<2;++i) /* Computation of the most probable coefficients */
        {            /* by the least square method */
	coeff[i]=db(0);
	for(j=0;j<N;++j)
	    coeff[i]+=phi[i]((long int)j,(long int)N)*Y[j];
        }
    var=moy=db(0);
    for(j=0;j<N;++j) /* Computation of the residuals */
        moy+=Y[j]=Y[j]-extrapol(j,N);
    moy/=db(N);
    for(j=0;j<N;++j) /* Computation of the variance of the residuals   */
        {
        coco=Y[j]-moy;
	coco*=coco;
	var+=coco;
        }
    var/=db(N-lg(1));
    dx=T[N-1]-T[0];
    dy=extrapol(N-1,N)-extrapol(0,N);
    a=dy/dx;
    b=extrapol(0,N)-a*T[0];
    coeff[0]=b;
    coeff[1]=a;
    return(var);
    }

