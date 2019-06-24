/*   devgraph.c                                   F. Vernotte - 2019/06/21  */
/*   Plot one or several ADEV files                                         */
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

#define db(x) ((double)(x))
#define sisig(x) ( (x) == 0 ) ? (db(0)) : (  ( (x) > 0 ) ? (db(1)) : (db(-1))  )

void usage(void)
/* Help message */
    {
    printf("Usage: DevGraph OutFILE DevFILE1 [Dev FILE2] ... [DevFILE8]\n\n");
    printf("Plots the graph of 1 to 8 2-columns Dev files (suitable to plot the 3 ADEV computed by GCoDev on a three-cornered hat system and the ADEV of the measurement noise computed by Aver).\n\n");
    printf("If a file name is preceded by a dash (-DevFILEn), the contain of the second column of this file is multiplied by -1.\n\n");
    printf("The OutFILE is given to store:\n");
    printf("- the file OutFILE.gnu for invoking gnuplot,\n");
    printf("- the file OutFILE.pdf which is the pdf file of the gnuplot graph (if the PDF option has been chosen in the configuration file).\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
    {
    int i, j, nbf, N, err, ind_gcod;
    char outfile[256], infiles[8][256], intername[256];
    double tau[256], gcod[8][256];

    if (argc<2)
        {
        usage();
	exit(-1);
	}
    else
	{
	if (argc>10)
		{
		printf("# Too many file names\n\n");
		usage();
		exit(-1);
		}
	else
		{
		nbf=argc-2;
		strcpy(outfile,*++argv);
		for(i=0;i<argc-2;++i)
			strcpy(infiles[i],*++argv);
		}
	}

    err=init_flag();
    if (err==-1) printf("# ~/.SigmaTheta.conf not found, default values selected\n");
    if (err==-2) printf("# ~/.SigmaTheta.conf improper, default values selected\n");

    for(i=0;i<nbf;++i)
	{
	if (infiles[i][0]!='-')
		{
		intername[0]='+';
		intername[1]=(char)0;
		strcat(intername,infiles[i]);
		strcpy(infiles[i],intername);
		}
	N=load_adev(&infiles[i][1],tau,gcod[i]);
    	if (N==-1)
		{
		printf("# File %s not found\n",&infiles[i][1]);
		exit(-1);
		}
        if (N<2)
	    	{
		printf("# %s: unrecognized file\n\n", &infiles[i][1]);
	    	usage();
		exit(-1);
		}
	}
	if (nbf>4) ind_gcod=0;
	else ind_gcod=1;

/* Use of gnuplot for generating the graph as a ps file */
    err=gen_gcodplt(outfile, infiles, N, nbf, tau, gcod, ind_gcod);
    if (err) printf("# Error %d: ps file not created\n",err);
    }
		  
