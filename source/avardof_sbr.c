/*   avardof_sbr.c                                 F. Vernotte - 2010/10/24 */
/*   Subroutines for the computation of the degrees of freedom of an Allan  */
/*   variance sequence following the paper "Uncertainty of Stability        */
/*   Variances Based on Finite Differences" by C. Greenhall and W. Riley    */
/*   (35th PTTI, San Diego, 2003, pp. 267-278)                              */
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

#define db(x) ((double)(x))
#define sisig(x) ( (x) == 0 ) ? (db(0)) : (  ( (x) > 0 ) ? (db(1)) : (db(-1))  )

extern char flag_variance;

double sw(double t,int alpha)
/* Sw(t,alpha) function defined in 35th PTTI p.270 */
    {
    double rslt;

    switch(alpha)
        {
        case 2 :
            rslt=-fabs(t);
	    break;
        case 1 :
            if (t!=0)
                rslt=t*t*log(fabs(t));
            else
	        rslt=db(0);
	    break;
        case 0 :
	    rslt=fabs(pow(t,db(3)));
	    break;
        case -1 :
            if (t!=0)
	        rslt=-pow(t,db(4))*log(fabs(t));
            else
                rslt=0;
	    break;
        case -2 :
	    rslt=-fabs(pow(t,db(5)));
	    break;
	default :
	    rslt=0;
	}
    return(rslt);
    }

double sx(double t, int Flt, int alpha)
/* Sx(t,F,alpha) function defined in 35th PTTI p.270 */
    {
    double rslt;

    rslt=pow(db(Flt),db(2))*(db(2)*sw(t,alpha)-sw(t-db(1)/db(Flt),alpha)-sw(t+db(1)/db(Flt),alpha));
    return(rslt);
    }

double sz(double t, int Flt, int alpha, int d)
/* Sz(t,F,alpha,d) function defined in 35th PTTI p.270*/
    {
    double rslt;

    switch(d)
	{
	case 1 :
		rslt=db(2)*sx(t,Flt,alpha)-sx(t-db(1),Flt,alpha)-sx(t+db(1),Flt,alpha);
		break;

	case 2 :
    		rslt=db(6)*sx(t,Flt,alpha)-db(4)*sx(t-db(1),Flt,alpha)-db(4)*sx(t+db(1),Flt,alpha)+sx(t-db(2),Flt,alpha)+sx(t+db(2),Flt,alpha);
		break;

	case 3 :
		rslt=db(20)*sx(t,Flt,alpha)-db(15)*sx(t-db(1),Flt,alpha)-db(15)*sx(t+db(1),Flt,alpha)+db(6)*sx(t-db(2),Flt,alpha)+db(6)*sx(t+db(2),Flt,alpha)-sx(t-db(3),Flt,alpha)-sx(t+db(3),Flt,alpha);
		break;

	default :
    		rslt=db(6)*sx(t,Flt,alpha)-db(4)*sx(t-db(1),Flt,alpha)-db(4)*sx(t+db(1),Flt,alpha)+sx(t-db(2),Flt,alpha)+sx(t+db(2),Flt,alpha);
	}
    return(rslt);
    }

double BasicSum(int J, int M, int S, int F, int alpha, int d)
/* BasicSum(J,M,S,F,alpha,d) function defined in 35th PTTI p.270 */
    {
    int jj;
    double ssj, rslt ;

    ssj=0;
    for(jj=1;jj<=(J-1);++jj)
        ssj+=(db(1)-db(jj)/db(M))*pow(sz(db(jj)/db(S),F,alpha,d),db(2));
    rslt=db(2)*ssj+pow(sz(db(0),F,alpha,d),db(2))+(db(1)-db(J)/db(M))*pow(sz(db(J)/db(S),F,alpha,d),db(2));
    return(rslt);
    }

void avardof(int la, double tau[32], int alpha[32], double edf[32])
    {
/* Computation of the degrees of freedom of the Allan variance estimates    */
/* (based on "Uncertainty of stability variances based on finite            */
/* differences" by C.A.Greenhall and W.J.Riley, 35th PTTI)                  */
    int i, N, m[32], d, F[32], S[32], L[32], M[32], J[32];

    for(i=0;i<la;++i)
	{
	m[i]=(int)(tau[i]/tau[0]);
	if (flag_variance==1) F[i]=1;
	else F[i]=m[i];
	S[i]=m[i];
	}
    N=2*m[la-1]+1;
    if (flag_variance) N=3*m[la-1]+1;
    else N=2*m[la-1]+1;
    if (flag_variance==2) d=3;
    else d=2;

/* Initial steps */
    for(i=0;i<la;++i)
	{
        L[i]=m[i]/F[i]+m[i]*d;
        M[i]=1+N-L[i];
	if (M[i]<((d+1)*S[i]))
	    J[i]=M[i];
	else
	    J[i]=(d+1)*S[i];
	}

/* Main procedure, simplified version */
    for(i=0;i<la;++i)
	{
	edf[i]=M[i]*pow(sz(db(0),F[i],alpha[i],d),db(2))/BasicSum(J[i],M[i],S[i],F[i],alpha[i],d);
	}
    }
