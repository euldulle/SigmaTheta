/*   find_tau_xdev.c                               F. Meyer    - 2021/04/22 */
/*   Interpolates a given tau from a 2col deviation table                   */
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
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include "sigma_theta.h"
#define NMAX 256
/*

find_tau_xdev.c

Filter aimed at interpolating a deviation value for a given tau
that is not present in the original data

data are read from either stdin (default) or from the optionally specified file

column 1 : tau (in s by default)
column 2 : deviation data (Hz/Hz)

tau are supposed in ascending order

*/

void usage(void)
/* Help message */
    {
	printf("##############################################################################################################\n\n");
    printf(" find_tau_allan : filtre pour interpoler une valeur donnée de tau d'un tableau de variance d'Allan\n\n");
    printf("     reads a 2 col table with tau in seconds and adev in Hz/Hz\n");
    printf("     Usage: find_tau_allan [-f tau ] [SOURCE]\n\n");
    printf("      Default behaviour (no file specified) is a filter, taking stdin as input and stdout as output.\n\n");
	printf("      If SOURCE is as an input file, specified output still goes to stdout \n");
    printf("           -x xscalingfactor\n");
	printf("                Units are SI units (s) by default ; should the input tau be in other units (MJD, ns, ...)\n");
    printf("                x option allows to properly rescale tau\n");
    printf("                  allowed values are the same as for the -f option and are listed below\n");
    printf("           -f tau to extract, under the form numvalue+suffix eg -f 86400s or -f 1d\n");
	printf("                Units are SI units (s) by default ; should the input data be in other units (MJD, ns, ...)\n");
    printf("                suffix is one of : \n");
    printf("                    d : days  \n");
    printf("                    H : hours  \n");
    printf("                    M : minutes  \n");
    printf("                    m : millisecond \n");
    printf("                    u : microsecond \n");
    printf("                    n : nanosecond \n");
    printf("                    p : picosecond \n");
    printf("                  A file containing data as tau in days vs sigma_y(tau) can be processed with \n");
    printf("                   find_tau_allan datafile -x s -f 1d \n\n");
    printf("                   to interpolate the adev value for tau = 1d \n\n");
    printf("           -X output scalingfactor for tau\n");
	printf("                output tau are are in SI units (s) by default \n");
    printf("                should the output data be wanted in other units (MJD, ns, ...)\n");
    printf("                X option allows to properly normalize output taus \n");
    printf("                see f option above for valid scaling factor\n");
    printf("           -h : this message\n\n");
    printf("     Input consists in an N line / 2 column table with tau values in second in the first column\n");
    printf("                                   and adev values (Hz/Hz) in the second column.\n\n");
    printf("     Ouput is a 2 column line with requested time tau in the first column \n");
    printf("                                   and Allan deviation measurement in column 2.\n\n");
//    printf("           SigmaTheta %s %s \n",st_version,st_date);
//    printf("           FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n");
	printf("##############################################################################################################\n\n");
    exit(-1);
    }

double suffix_process(uint8_t *suffix){
    uint8_t *last;
    uint16_t length;
    double scale_unit;
    double scale_factor;

    if (suffix == NULL)
        return (double)-1;

    last=suffix+strlen(suffix)-1;

    if (last <0)
        return (double)-1;

    scale_unit=(double)1; // by defaults tau are supposed in seconds
    if (isalpha(*last)){      // if there is a scale specifier in the suffix then
                             // get its value from scale()
        scale_unit=scale(*last);
        *last=0;
        }

    if (last==suffix && *last==0)   // only one char was given and it was a letter
        scale_factor=scale_unit;
    else
        scale_factor=atoi(suffix)*scale_unit;

    return scale_factor;
    }

int main(int argc, char **argv) {
    double input_scale=1, outscaletau=1, first_tau, last_tau, target_tau, cur_tau, prev_tau, target_dev;
    uint8_t stdo=1, index, rescale_input=0, rescale_tau=0; 
    char c;
    uint16_t lensrc;
    int nbv;
    uint8_t source[MAXCHAR],outfile[MAXCHAR];
    double tau[NMAX], dev[NMAX];
	char line[MAXCHAR];

	long i, j, nb;
	long nn,n=0;
	extern char *optarg;
	extern int opterr;

  	while(1) {
		c=getopt(argc, argv, "f:x:X:");
		if (c=='?') {
			usage();
			exit(0);
			}

		if (c==-1)
			break;

		switch(c) {
				case 'f':
                    target_tau=suffix_process(optarg);
					break;
					
				case 'x':
                    input_scale=suffix_process(optarg);
                    rescale_input=1;
					break;
					
				case 'X':
                    outscaletau=suffix_process(optarg);
                    rescale_tau=1;
					break;
					
				default:
					usage();
					exit(0);
				}
			}

    index=optind;
    if (index<argc){
        lensrc=strlen(strncpy(source,argv[index], MAXCHAR));

        index++;
        if (index<argc){
            strncpy(outfile,argv[index], MAXCHAR);
            printf ("# Output file %s\n", outfile);
            }
        if (index<argc)
        {
            fprintf(stderr,"Extra parameter.\n");
            usage();
        }
    }

    if (strlen(source)==0){
        fprintf(stderr,"#\n#\n# No input file given , expecting data on stdin...\n# %s -h to show usage )\n#\n", argv[0]);
        stdo=1;
        }
    //
    // target tau in seconds :
    //
    nbv=load_adev(source, tau, dev);
    if (rescale_input!=0){
        for (i=0;i<nbv; ++i){
            tau[i]*=input_scale;
            }
        }
    first_tau=tau[0];
    last_tau=tau[nbv-1];
    if (first_tau > target_tau){
        fprintf(stderr, "# \%s : Required value for tau %le is lower than lowest data tau (%le)\n", argv[0], target_tau, first_tau);
        printf("%le %le\n", target_tau, -dev[0]);
        return -(double) dev[0];
        }
    if (last_tau < target_tau){
        fprintf(stderr, "# \%s : Required value for tau %le is higher than highest data tau (%le)\n", argv[0], target_tau, last_tau);
            printf("%le %le\n", target_tau, -dev[nbv-1]);
        return -(double) dev[nbv-1];
        }
    cur_tau=first_tau;
    for (i=1;i<nbv; ++i){
        prev_tau=cur_tau;
        cur_tau=tau[i];
        if (cur_tau> target_tau){ // linear interpolation
            target_dev=dev[i-1] + (dev[i]-dev[i-1]) * (target_tau-tau[i-1])/(tau[i]-tau[i-1]) ;
            printf("%le %le\n", target_tau/outscaletau, target_dev);
            return target_dev;
            }
        }
    return -(double)1;
	}
