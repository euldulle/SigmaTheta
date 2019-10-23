/*   1col2col.c                                    F. Vernotte - 2010/10/31 */
/*   Transformation of a 1 column file into a 2 column file                 */
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
#include <unistd.h>
#include <ctype.h>
#include <stdint.h>
#include <math.h>
#include "sigma_theta.h"

#define DATAMAX 16384
#define GRANMAX 67108864

void usage(void)
/* Help message */
    {                                                                                                                      
	printf("##############################################################################################################\n\n");
    printf(" 1col2col : a tool from the SigmaTheta suite\n\n");
    printf("     Usage: 1col2col [-ch] [-n normfreq] [-o normfreq] [-t tau] [SOURCE [TARGET]]\n\n");
    printf("     Transforms a 1 column file into a 2 column file.\n\n");
    printf("      Default behaviour (no file specified) is a filter, taking stdin as input and stdout as output.\n\n");
	printf("      If SOURCE and TARGET are specified, output goes to TARGET\n");
	printf("      If only SOURCE is specified, output goes to SOURCE.2col unless -c option (output to stdout) is given \n\n");
    printf("           -n normfreq\n");
	printf("                If input data are absolute frequency deviations, -n normalizes to the frequency given as arg\n");
    printf("           -o normfreq\n");
	printf("                If input data are absolute frequency samples, -o offsets and normalizes them to the frequency given as arg\n");
    printf("           -c : output to stdout only, TARGET file is ignored even if specified ;\n");
    printf("                this is the default if SOURCE is unspecified (stdin)\n\n");
    printf("           -t tau\n");
	printf("                sampling time to be used to generate timing coordinates (defaults to 1s) \n");
    printf("           -h : this message\n\n");
    printf("     Input consists in an N line / 2 column table with time values (dates) in the first column\n");
    printf("                                   and normalized frequency deviation samples in the second column.\n\n");
    printf("     Ouput is a 2 column table with averaging times tau in the first column \n");
    printf("                                   and Allan deviation measurement in column 2.\n\n");
    printf("           SigmaTheta %s %s \n",st_version,st_date);
    printf("           FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n");                                   
	printf("##############################################################################################################\n\n");
    exit(-1);
    }

int main(int argc, char *argv[]) {
    int i, nbv, N;
    long int dtmx;
    char source[MAXCHAR]="", outfile[MAXCHAR]="";
    FILE *ofd;
    double tau=(double)1;
    char offset=0, c, stdo=0;
    size_t index, lensrc, stridx;
    double nominalfreq=(double) 1;

    opterr = 0;

    while ((c = getopt (argc, argv, "hcn:o:t:")) != -1)
        switch (c)
        {
            case 'c':
                stdo = 1;
                break;
            case 'h':
                usage();
                break;
            case 'n':
                nominalfreq=atof(optarg);
                break;
            case 'o':
                offset=1;
                nominalfreq=atof(optarg);
                break;
            case 't':
                tau=atof(optarg);
                break;
            case '?':
                if (optopt == 't'){
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                    usage();
                }
                else if (optopt == 'o'){
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                    usage();
                }
                else if (optopt == 'n'){
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

        index++;
        if (index<argc){
            strncpy(outfile,argv[index], MAXCHAR);
            fprintf (stderr, "# Output file %s\n", outfile);
        }
        else {
            if (stdo != 1){
                stridx=lensrc-3;
                if (source[stridx++]=='y' && source[stridx++]=='k' && source[stridx++]=='t' ){
                    strncpy(outfile, source, MAXCHAR);
                    snprintf(outfile+lensrc-3, 5, "adev");
                }
                else{
                    snprintf(outfile, lensrc+6, "%s.adev", source);
                }
            }
        }

        index++;
        if (index<argc) { // for compatibility reasons we process the 3rd arg as back in the days
            // though it should be carried via -t option
            tau=atof(argv[index]);
        }

        index++;
        if (index<argc) {
            fprintf(stderr,"Extra parameter.\n");
        }
    }

    if (strlen(source)==0){
        fprintf(stderr,"#\n#\n# No input file given , expecting data on stdin...\n# %s -h to show usage )\n#\n", argv[0]);
        stdo=1;
    }

    fprintf(stderr,"in @%s@ out @%s@ normf %le off %d tau %le \n", source, outfile, nominalfreq, offset, tau);

    N=load_1col(source);
    if (N==-1)
        printf("# File not found\n");
    else
    {
        if (N<2)
        {
            if (N==-2) printf("# %s is already a 2-column table\n",source);
            else
            {
                printf("# Unrecognized file\n\n");
                usage();
            }
        }
        else
        {
            if (stdo == 0){
                ofd=fopen(outfile, "w");
            }
            else{
                ofd=stdout;
            }

            if (ofd==NULL)
                printf("# Incorrect file name\n");
            else
            {
                fprintf(ofd,"# File generated by 1col2col from the file %s with tau=%g s\n",source,tau);
                for(i=0;i<N;++i)
                    fprintf(ofd,"%24.16e \t %24.16e\n",tau*((double)i),(Y[i]-offset*nominalfreq)/nominalfreq);
                fclose(ofd);
            }
        }
    }
}

