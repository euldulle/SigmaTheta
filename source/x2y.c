/*   x2y.c                                        F. Vernotte - 2010/10/27  */
/*   Transformation of a time error sequence {x(t)} into a normalized       */
/*   frequency deviation sequence {Yk}                                      */
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
                                    
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "sigma_theta.h"

#define DATAMAX 16384
#define GRANMAX 67108864

#define MAXCHAR 512

void usage(void)
/* Help message */
    {
    printf("##############################################################################################################\n\n");
    printf(" X2Y : a tool from the SigmaTheta suite\n\n");
    printf("     Usage: X2Y [-ch] [-x xscalingfactor] [-y yscalingfactor] [SOURCE [TARGET]]\n\n");
    printf("     Transforms a time error sequence {x(t)} into a normalized frequency deviation sequence {Yk}.\n\n");
    printf("      Default behaviour (no file specified) is a filter, taking stdin as input and stdout as output.\n\n");
    printf("      If only SOURCE is specified, output goes to SOURCE.ykt unless -c option (output to stdout) is given \n\n");
    printf("      Options :\n");
    printf("           -x xscalingfactor\n");
    printf("           -y xscalingfactor\n");
    printf("            	Units are SI units (s) by default ; should the input data be in other units (MJD, ns, ...)\n");
    printf("            	x and y options allow to properly normalize output : \n");
    printf("            	scaling factor is one of : \n");
    printf("           	  	  	d : days  \n");
    printf("           	  	  	H : hours  \n");
    printf("           	  	  	M : minutes  \n");
    printf("           	  	  	m : millisecond \n");
    printf("           	  	  	u : microsecond \n");
    printf("           	  	  	n : nanosecond \n");
    printf("           	  	  	p : picosecond \n");
	printf("            	  A file containing data as MJD.XXXXX vs time in ns can be processed with \n");
	printf("            	  X2Y datafile -x d -y n \n\n");
    printf("           -c : stdout output ; this is the default if SOURCE is stdin\n");
    printf("                if SOURCE is specified, both stdout and SOURCE.ykt are fed with the results.)\n");
    printf("                if SOURCE and TARGET are specified, both stdout and TARGET are fed with the results.)\n\n");
    printf("           -h : this message\n\n");
                
	            
    printf("     Input consists in an N line / 2 column table with time values (dates) in the first column\n");
    printf("                                                   and time error samples in the second column.\n\n");
    printf("     Ouput is a N-1 line / 2 column table with time values (dates) in the first column \n");
    printf("                                           and normalized frequency samples in the second column.\n\n");
    printf("           SigmaTheta %s %s \n",st_version,st_date);
    printf("           FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n");
    printf("##############################################################################################################\n\n");
    exit(-1);
    }

int main(int argc, char *argv[])
    {
    int i, nbv, N;
    long int dtmx;
    char source[MAXCHAR]="", outfile[MAXCHAR]="", gm[100];
    FILE *ofd;
    double tau;
    char scalex=0, scaley=0;
    uint8_t stdo=0;
    int8_t c, index;
    int16_t stridx, lensrc;

	opterr = 0;

    while ((c = getopt (argc, argv, "hcx:y:")) != -1)
        switch (c)
        {
            case 'c':
                stdo = 1;
                break;
            case 'h':
                usage();
                break;
            case 'x':
                scalex=optarg[0];
                break;
            case 'y':
                scaley=optarg[0];
                break;
            case '?':
                if (optopt == 'x'){
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                    usage();
                }
                else if (optopt == 'y'){
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                    usage();
                }
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return 1;
            default:
                abort ();
        }

    index=optind;
	if (index<argc){
		lensrc=strlen(strncpy(source,argv[index], MAXCHAR));
        printf ("Input file %s\n", argv[index]);

        index++;
        if (index<argc){
            printf ("Output file %s\n", argv[index]);
            strncpy(outfile,argv[index], MAXCHAR);
            }
        else {
            if (stdo != 1){
                stridx=lensrc-3;
                if (source[stridx++]=='x' && source[stridx++]=='k' && source[stridx++]=='t' ){
                    snprintf(outfile, lensrc+1, "%s", source);
                    outfile[lensrc-3]='y';
                    }
                else{
                    snprintf(outfile, lensrc+5, "%s.ykt", source);
                }
            }
        }

        index++;
        if (index<argc)
        {
            fprintf(stderr,"Extra parameter.\n");
            usage();
        }
    }

    if (strlen(source)==0){
        fprintf(stderr,"#\n#\n# No input file given , expecting data on stdin...\n#   (X2Y -h to show usage)\n");
        stdo=1;
        }

    //fprintf(stderr,"# scalex %d %c / scaley %d %c", scalex, scalex, scaley, scaley);
    N=load_ykt(source, scalex , scaley);
	if (N==-1)
		printf("# File not found\n");
	else
	{
		if (N<2)
		{
			printf("# Unrecognized file\n\n");
			if (N==1)
				printf("# Use 1col2col command\n\n");
			usage();
		}
		else
		{
			tau=T[1]-T[0];
            if (stdo==1)
                ofd=stdout;
            else
                ofd=fopen(outfile, "w");

			if (ofd==NULL)
				fprintf(stderr, "# Could not open file %s\n", source);
			else
			{
				fprintf(ofd,"# File generated by X2Y from the file %s\n",source);
				for(i=0;i<N-1;++i)
					fprintf(ofd,"%24.16e \t %24.16e\n",(T[i]-T[0]),(Y[i+1]-Y[i])/tau);
                if (ofd!=stdout)
                    fclose(ofd);
			}
		}
	}
	}

