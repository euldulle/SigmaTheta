/*   rayleigh.c                                   F. Vernotte - 2010/10/19  */
/*   Computation of a 95 % confidence interval for a Rayleigh               */
/*   distribution with 'nu' degrees of freedom:                             */
/*     input : nu                                                           */
/*     output : inf bound, log mean, sup bound                              */
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
     #include "sigma_theta.h"
  
void usage(void)
/* Help message */
    {
    printf("Usage: RaylConfInt 'value'\n\n");
    printf("Computes the mean and the 95 %% confidence interval of a Chi distribution with 'value' degrees of freedom, normalized by the square root of 'value'.\n\n");
    printf("The input 'value' is a floating point number.\n\n");
    printf("The inferior bound, the logarithmic mean and the superior bound are sent to the standard output separated by a tabulation.\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
/* 95 % confidence interval bounds for a Rayleigh distribution with 'nu' degrees of freedom
   input : nu 
   output : inf bound, mean, sup bound */
    {
    int i;
    double nu,a,b;
    struct conf_int rayl;

    if (argc<2)
	usage();
    else
        {
	nu=atof (*++argv);
        rayl=raylconfint(nu);
        printf("%g \t %g \t %g\n",sqrt(nu)/rayl.sup_bound,((double)1)/rayl.mean,sqrt(nu)/rayl.inf_bound);
	/*        rayl=raylconfint1s(nu);
		  printf("%g \t %g \t %g\n",sqrt(nu)/rayl.sup_bound,((double)1)/rayl.mean,sqrt(nu)/rayl.inf_bound);*/
        }
    }
