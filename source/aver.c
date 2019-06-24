/*   aver.c                                    F. Vernotte - 2018/10/19     */
/*   Average rounded to integer of 2nd columns of file 1 and file 2         */
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
    printf("Usage: AVer FILE1 FILE2 FILE3\n\n");
    printf("Computes the average of the 2nd columns of FILE1, FILE2 and FILE3 (suitable to compute the measurement noise of a three-cornered system thanks to the closure relationship).\n\n");
    printf("The input files FILE1, FILE2 and FILE3 contain a N line / 2 column table, the 1st columns are assumed identical.\n\n");
    printf("A 2-column table containing the 1st column of FILE1 as the first column of the output table and the average of the 2nd columns of FILE1, FILE2, FILE3 as the the second column of the output table is sent to the standard output.\n\n");
    printf("A redirection should be used for saving the results in a TARGET file: AVer FILE1 FILE2 > TARGET.\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
    {
    int indeq, err;
    int i, nbv, N, N2, N3;
    char source1[256], source2[256], source3[256], gm[100];

    if (argc<4)
        {
	usage();
	exit(-1);
        }
    else
        {
	strcpy(source1,*++argv);
	strcpy(source2,*++argv);
	strcpy(source3,*++argv);
	}
    err=0;
    N=load_3yk(source1,source2,source3);
    if (N<1)
	{
	switch(N)
		{
	        case 0:
			printf("# Empty file\n");
			break;
		case -1:
			printf("# File not found\n");
			break;
		case -2:
			printf("# Different file lengths\n");
		}
	exit(N);
	}
    for (i=0;i<N-1;++i)
	printf("%24.16e \t %24.16e\n",T[i],(Y1[i]+Y2[i]+Y[i])/((double)3));
    exit(err);		
    }
    
