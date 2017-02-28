/*   rayleigh.c                                   F. Vernotte - 2010/10/19  */
/*   Rayleigh distribution subroutines                                      */
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
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

struct conf_int
    {
    double inf_bound;
    double mean;
    double sup_bound;
    } ;

double cdf_rayleigh(double x, double nu)
/* Cumulative distribution function of a Rayleigh distribution with 'nu' degrees of freedom for the value 'x'*/
    {
    double rslt;
    rslt=gsl_sf_gamma_inc_P(nu/2,x*x/2);
    return(rslt);
    }

double dcdfr(double x, double nu)
/* Derivative of the cumulative distribution function of a Rayleigh distribution with 'nu' degrees of freedom for the value 'x'*/
    {
    double rslt;
    rslt=pow(2.,1.-nu/2.)*x*pow(x*x,-1.+nu/2.)/(exp(x*x/2.)*gsl_sf_gamma(nu/2.));
    return(rslt);
    }

struct conf_int raylconfint(double nu)
/* 95 % confidence interval bounds for a Rayleigh distribution with 'nu' degrees of freedom
   input : nu 
   output : inf bound, mean, sup bound */
    {
    int i;
    double mu,ecty,b1,b2,x0,x1,tstop;
    struct conf_int rslt;

    x0=1.;
    if (nu<250)
      mu=sqrt(2.)*gsl_sf_gamma((nu+1.)/2.)/gsl_sf_gamma(nu/2.);
    else
      mu=sqrt(nu-1.);
    ecty=sqrt(nu-mu*mu);
    b1=mu-2*ecty;
    b2=mu+2*ecty;
    if (nu<275)
        {
        if (b1<=1)
            x1=sqrt(nu/2);
        else
            x1=b1;
        do
            {
            x0=x1;
             x1=x0-(cdf_rayleigh(x0,nu)-((double)0.025))/dcdfr(x0,nu);
            x1=fabs(x1);
	    tstop=fabs((x1-x0)/x1);
	    }
        while(tstop>1e-9); 
        b1=fabs(x1);
        x1=b2;
        do
            {
            x0=x1;
	    x1=x0-(cdf_rayleigh(x0,nu)-0.975)/dcdfr(x0,nu);
	    x1=fabs(x1);
	    tstop=fabs((x1-x0)/x1);
	    }
        while(tstop>1e-9); 
        b2=fabs(x1);
	}
    rslt.inf_bound=b1;
    rslt.mean=sqrt(((double)2))*exp(gsl_sf_psi(nu/((double)2))/((double)2))/sqrt(nu);
    rslt.sup_bound=b2;
    return(rslt);
    }

struct conf_int raylconfint1s(double nu)
/* 95 % confidence interval bounds for a Rayleigh distribution with 'nu' degrees of freedom
   input : nu 
   output : inf bound, mean, sup bound */
    {
    int i;
    double mu,ecty,b1,b2,x0,x1,tstop;
    struct conf_int rslt;

    x0=1.;
    if (nu<250)
      mu=sqrt(2.)*gsl_sf_gamma((nu+1.)/2.)/gsl_sf_gamma(nu/2.);
    else
      mu=sqrt(nu-1.);
    ecty=sqrt(nu-mu*mu);
    b1=mu-ecty;
    b2=mu+ecty;
    if (nu<275)
        {
        if (b1<=1)
            x1=sqrt(nu/2);
        else
            x1=b1;
        do
            {
            x0=x1;
             x1=x0-(cdf_rayleigh(x0,nu)-((double)0.16))/dcdfr(x0,nu);
            x1=fabs(x1);
	    tstop=fabs((x1-x0)/x1);
	    }
        while(tstop>1e-9); 
        b1=fabs(x1);
        x1=b2;
        do
            {
            x0=x1;
	    x1=x0-(cdf_rayleigh(x0,nu)-0.84)/dcdfr(x0,nu);
	    x1=fabs(x1);
	    tstop=fabs((x1-x0)/x1);
	    }
        while(tstop>1e-9); 
        b2=fabs(x1);
	}
    rslt.inf_bound=b1;
    rslt.mean=sqrt(((double)2))*exp(gsl_sf_psi(nu/((double)2))/((double)2))/sqrt(nu);
    rslt.sup_bound=b2;
    return(rslt);
    }

