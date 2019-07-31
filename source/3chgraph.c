/*   3CHGraph.c                                   F. Vernotte - 2019/07/31  */
/*   Plot 3-cornered hat results with confidence intervals 		    */
/*   in a single graph 							    */
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

#define db(x) ((double)(x))
#define sisig(x) ( (x) == 0 ) ? (db(0)) : (  ( (x) > 0 ) ? (db(1)) : (db(-1))  )

void usage(void)
/* Help message */
    {
    printf("Usage: 3CHGraph [-m|f] OutFILE InFILE1 InFILE2 InFILE3 [InFILE4]\n\n");
    printf("Plots 3-cornered hat results (computed by Groslambert Covariance or classical 3-cornered hat) with confidence intervals in a single graph.\n\n");
    printf("The 3 input files InFILE1-3 have 7 columns. Each of these files should contain the results of each of the 3 clocks of the 3-cornered hat system: tau, GCov or classical 3-conered hat, mean or median estimate, 68 and 95 %% confidence interval bounds.\n\n");
    printf("The optional 4th input file InFILE4 has 2 columns. It should contain the ADev of the measurement noise.\n\n");
    printf("If the option '-m' is selected, the mean estimate is meant to be displayed in the 3rd column of the input files.\n");
    printf("If the option '-f' is selected, the 50 %% estimate (median) is meant to be displayed in the 3rd column of the input files.\n\n");
    printf("The OutFILE is given to store:\n");
    printf("- the file OutFILE.gnu for invoking gnuplot,\n");
    printf("- the file OutFILE.pdf which is the pdf file of the gnuplot graph (if the PDF option has been chosen in the configuration file).\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
    {
    int i, j, rep, nbf, err, N, Nref, aflag, fflag, ind50;
    char outfile[256], infiles[4][256];
    double tau[256], gcod[4][256], bmin[3][256], bmax[3][256], truc[256];

    aflag=fflag=0;
    ind50=2;
    while ((rep = getopt (argc, argv, "mf")) != -1)
	switch(rep)
		{
		case 'm':
			aflag = 1;
			ind50=0;
			break;
		case 'f':
			fflag = ind50 = 1;
			break;
		case '?':
			printf("# Unknown option '-%c'\n",optopt);
			usage();
		default:
			exit(-1);
		}
    if (aflag&fflag)
		{
		printf("# Incompatible options '-m' and '-f'\n");
		usage();
		exit(-1);
		}
    nbf=argc-optind;
    if (nbf<4)
        {
	printf("# 3CHGraph requires at least 4 file names as arguments\n\n");
        usage();
	exit(-1);
	}
    if (nbf>5)
	{
	printf("# Too many file names\n\n");
	usage();
	exit(-1);
	}
    strcpy(outfile,argv[optind]);
    for (i=0;i<3;++i)
	{
	strcpy(infiles[i],argv[optind+i+1]);
	N=load_7col(infiles[i],tau,gcod[i],bmin[i],truc,truc,truc,bmax[i]);
    	if (N==-1)
		{
      		printf("# File %s not found\n", infiles[i]);
		exit(-1);
		}
        if (N<2)
		{
		printf("# %s: unrecognized file\n\n", infiles[i]);
		usage();
		exit(-1);
		}
	for (j=0;j<N;++j)
		if (gcod[i][j]<bmin[i][j]) bmin[i][j]=gcod[i][j];
	if (i)
		{
		if (N!=Nref)
			{
			printf("# Different number of data in %s and %s\n",infiles[0],infiles[i]);
			exit(-1);
			}
		}
	else
		Nref=N;
	}
    if (nbf==5)
	{
	strcpy(infiles[3],argv[optind+4]);
	N=load_adev(infiles[3],tau,gcod[3]);
    	if (N==-1)
		{
      		printf("# File %s not found\n", infiles[3]);
		exit(-1);
		}
        if (N<2)
		{
		printf("# %s: unrecognized file\n\n", infiles[3]);
		usage();
		exit(-1);
		}
	}
    else
	infiles[3][0]=0;

    err=init_flag();
    if (err==-1) printf("# ~/.SigmaTheta.conf not found, default values selected\n");
    if (err==-2) printf("# ~/.SigmaTheta.conf improper, default values selected\n");


/* Use of gnuplot for generating the graph as a ps file */
    err=gen_3chplt(infiles, outfile, N, tau, bmin, bmax, ind50);
    if (err) printf("# Error %d: pdf file not created\n",err);
    }
		  
