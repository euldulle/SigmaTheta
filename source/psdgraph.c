/*   psdgraph.c                                   F. Vernotte - 2015/06/24  */
/*   Plot the Power Spectral Density of a .ykt file (normalized frequency   */
/*   deviation) versus the frequency                                        */
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
#include <gsl/gsl_fft_real.h>
#include "sigma_theta.h"

#define db(x) ((double)(x))
#define sisig(x) ( (x) == 0 ) ? (db(0)) : (  ( (x) > 0 ) ? (db(1)) : (db(-1))  )

void usage(void)
/* Help message */
    {
    printf("Usage: PSDGraph yktFILE\n\n");
    printf("Plots the graph of the Power Spectrum Density (PSD) of normalized frequency deviation versus the frequency.\n\n");
    printf("The input file yktFILE contains a 2-column table with dates (in s) in the first column and frequency deviation samples (yk) in the second column.\n\n");
    printf("The PSD versus the frequency are stored in the output file yktFILE.psd.\n");
    printf("The file yktFILE.psd.gnu is generated for invoking gnuplot.\n");
    printf("The file yktFILE.psd.pdf is the pdf file of the gnuplot graph (if the PDF option has been chosen in the configuration file)\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
    {
    char source[256], psdfile[256], command[32];
    int N,M,err,i,j,nc;
    double stride,l2n,tot_dur,ksy;
    char fv;
    FILE *ofd;
    char scalex=0, scaley=0;

    fv=0;
    if (argc!=2)
        {
        usage();
	exit(-1);
	}
    else
	strcpy(source,*++argv);

    err=init_flag();
    if (err==-1) printf("# ~/.SigmaTheta.conf not found, default values selected\n");
    if (err==-2) printf("# ~/.SigmaTheta.conf improper, default values selected\n");
    flag_variance=fv;
    strcpy(psdfile,source);
    strcat(psdfile,".psd");
    N=load_ykt(source,scalex,scaley);
    if (N==-1)
      printf("# File %s not found\n",source);
    else
        {
        if (N<2)
	    {
	    printf("# %s: unrecognized file\n\n", source);
	    usage();
	    }
        else
	    {
	    stride=T[1]-T[0];
	    l2n=floor(log((double)N)/log((double)2));
	    N=(int)pow((double)2,l2n);
	    gsl_fft_real_radix2_transform (Y, 1, N);
	    M=N/2;
	    for (i=1;i<N/2;++i) Y[i-1]=Y[i]*Y[i]+Y[N-i]*Y[N-i];
	    Y[N/2-1]=Y[N/2]*Y[N/2];
	    tot_dur=((double)N)*stride;
	    ksy=((double)N)/(((double)2)*stride);
    	    ofd=fopen(psdfile, "w");
	    for (i=0;i<M;++i)
		{
		T[i]=((double)i+1)/tot_dur;
		Y[i]/=ksy;
		fprintf(ofd,"%e %e\n",T[i],Y[i]);
		}
	    fclose(ofd);
/* Use of gnuplot for generating the graph as a ps file */
	    err=gen_psdplt(psdfile,M,T,Y);
	    if (err) printf("# Error %d: ps file not created\n",err);
	    }
	}
    }
		  
