/*   drirem.c                                     F. Vernotte - 2010/10/20  */
/*   Removal of a linear frequency drift                                    */
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

#define DATAMAX 16384
#define GRANMAX 67108864

void usage(void)
/* Help message */
    {
    printf("Usage: DriRem SOURCE TARGET\n\n");
    printf("Removes the linear drift of a sequence of normalized frequency deviation measurements.\n\n");
    printf("The input file SOURCE contains a 2-column table with time values (dates) in the first column and normalized frequency deviation measurements in the second column.\n\n");
    printf("The output file TARGET contains a 2-column table with the same time values in the first column and the normalized frequency deviation measurements from which the linear drift was removed in the second column.\n\n");
    printf("The coefficients a and b of the drift (equation y=a.x+b) are given in the header of the TARGET file.\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
/* Remove the linear drift of a normalized frequency deviation sequence */
    {
    int i,nbv,N;
    long int dtmx;
    char source[256], outfile[256], gm[100];
    char scalex=0, scaley=0;
    FILE *ofd;

    if (argc<3)
        {
	usage();
	exit(-1);
        }
    else
        {
	strcpy(source,*++argv);
	strcpy(outfile,*++argv);
	}
    N=load_ykt(source,scalex,scaley,0);
    if (N==-1)
        printf("# File not found\n");
    else
        {
        if (N<2)
	    {
            printf("# Unrecognized file\n");
	    if (N==1)
	      printf("# Use 1col2col command\n\n");
	    usage();
	    }
        else
            {
	    tchebyfit((long)N);
	    ofd=fopen(outfile, "w");
            if (ofd==NULL)
                printf("# Incorrect file name\n");
	    else
		{
		fprintf(ofd,"# File generated by DriRem from the file %s\n",source);
		fprintf(ofd,"# Drift coefficients (y=a.t+b): a=%24.16e \t b=%24.16e\n",coeff[1],coeff[0]);
		printf("# Drift coefficients (y=a.t+b): \n");
		printf("# a=%24.16e \t b=%24.16e\n",coeff[1],coeff[0]);
		for(i=0;i<N;++i)
		    fprintf(ofd,"%24.16e \t %24.16e\n",T[i],Y[i]);
                fclose(ofd);
		}
	    }
        }
    }
    
