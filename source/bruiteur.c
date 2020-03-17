/*   bruiteur.c                     F. Vernotte - First release: 1994/02/02 */
/*					    Sigma-Theta version: 2015/10/26 */
/*                                    Modified by Attila Kinali: 2017/06/11 */
/*   Simulation of time error x(t) and/or frequency deviation yk samples    */
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/sysinfo.h>
#include "filtre.h"
#include "sigma_theta.h"


double *x;
double hm3,hm2,hm1,h0,hp1,hp2,C1,C0;
double tau0;
long GR;

struct sysinfo sinfo;

void usage(void)
/* Help message */
    {
    printf("Usage: bruiteur TARGET\n\n");
    printf("Computes a sequence of time error x(t) and/or normalized frequency deviation yk samples.\n\n");
    printf("A 2-column table containing dates in the first column and time errors in the second column is saved in the file TARGET.xtt and/or a 2-column table containing dates in the first column and frequency deviations in the second column is saved in the file TARGET.ykt.\n\n");
    printf("A redirection should be used for loading the input parameters from a SOURCE file: bruiteur TARGET < SOURCE\n\n");
    printf("SigmaTheta %s %s - FEMTO-ST/OSU THETA/Universite de Franche-Comte/CNRS - FRANCE\n",st_version,st_date);
    }

int main(int argc, char *argv[])
/* Simulation of time error x(t) samples and/or frequency deviation yk      */
/* samples.								    */
/* The samples are assumed to be composed of 5 power law noises from 	    */
/* f^{-2} FM to f^{+2} FM, i.e. from f^{-4} PM to white PM, and of a 	    */
/* deterministic linear frequency drift.				    */ 
	{
	FILE *fiptxt, *fiptyk;
	static char suffxt[4]= "xtt\0";
	static char suffyk[4]= "ykt\0";
	char nxtfich[64],nykfich[64];
	char fifich[64],rep;
	long i,ordre,nbr_dat,nbr_stat,debut,size_mem;
	int nb_seq,is;
	double mx,m2x,sigx,sigy,di,dx,dN,SX,SY,SXY,SX2,den,num,a0,a1,kdb,cpu_time;
	double Sm,Sm2,Ssig,Ssig2,Sa0,Sa02,Sa1,Sa12,adjsig,Sadj,Sadj2,un,deux,six;
	clock_t start, end;

	if (argc!=2)
        	{
		usage();
		exit(-1);
        	}
	else
		{
		strcpy(fifich,*++argv);
		}
	// get sysinfo to determine memory size
	if(sysinfo(&sinfo))
	{
		printf("sysinfo failed, something is wrong, seriously wrong!\n");
		exit(1);
	}
	strcpy(nxtfich,mod_nom(fifich,suffxt));
	strcpy(nykfich,mod_nom(fifich,suffyk));

	hm3=hm2=hm1=h0=hp1=hp2=C1=C0=(double)0;
	printf("Inverse of the low cut-off frequency (power of 2, < %ld): ", sinfo.totalram/2/(32+8));
/* For generating a sequence of N samples with a low cut-off frequency much lower than 1/(N tau0), we generate a sequence of M samples; we keep only a subsequence of N consecutive data (the beginning is randomly chosen). The M-N other data are sent to the trash!                                                   */
	scanf("%ld",&nbr_dat);
	if (nbr_dat<4) nbr_dat=4;
	// ensure that the allocated data is smaller than half of total memory
	// in order to avoid heavy swapping
	if (nbr_dat * (32+8) * 2 >sinfo.totalram)
	{
		printf("\n\nError: Requested number of samples would require to store data more than\nhalf the total RAM size\n");
		printf("Aborting to prevent the system from comming to a grinding halt\n");
		exit(1);
	}
	nbr_dat=(long)ceil(log((double)nbr_dat)/log((double)2));
	nbr_dat=(long)pow((double)2,(double)nbr_dat);
	printf("Low cut-off frequency: 1/%ld.\n",nbr_dat);
	if (nbr_dat>65536) size_mem=nbr_dat;
	else size_mem=(long)65536;
	printf("Memory size: %ld\n",size_mem);
	x=(double *)malloc(size_mem*sizeof(double)); 
	if ((x==NULL))
		{
		printf("Not enough memory for %ld data\n",nbr_dat);
		exit(-1);
		}
	printf("Number of samples (< or =  %ld): ",nbr_dat);
	scanf("%ld",&nbr_stat);
	if (nbr_stat<4) nbr_stat=4;
	if (nbr_stat>nbr_dat) nbr_stat=nbr_dat;
	printf("%ld samples selected\n",nbr_stat);
	dN=(double)nbr_stat;

	printf("Sampling step: ");
	scanf("%lf",&tau0);

	printf("\n\nEnter the values of the filtering coefficients: \n");
/*	printf("\th-3  ( coef. en 1/f3 ) : "); (simulation of millisecond pulsars)
	scanf("%lf",&hm3);*/
	hm3=0;
	printf("\th-2  (1/f2 coef.):     ");
	scanf("%lf",&hm2);
	printf("\th-1  (1/f coef.):      ");
	scanf("%lf",&hm1);
	printf("\th0   (constant coef.): ");
	scanf("%lf",&h0);
	printf("\th+1  (f coef.):        ");
	scanf("%lf",&hp1);
	printf("\th+2  (f2 coef.):       ");
	scanf("%lf",&hp2);
/* The high cut-off frequency is not arbitrarily set at the Nyquist frequency 	*/
/*	hp2*=tau0;								*/
	printf("\tC1 (linear frequency drift coef.):   ");
	scanf("%lf",&C1);
	printf("\tC0 (constant frequency drift coef.): ");
	scanf("%lf",&C0);

//	start=clock();
	kdb=((double)(nbr_dat-nbr_stat))/((double)RAND_MAX);/*Maximum value of the beginning index.*/
	adjsig=gausseq(nbr_dat,0); /* Generation of a unity white noise sequence. */
       	filtreur(nbr_dat,tau0); /* Filtering of the sequence.                           */
       	if (nbr_dat==nbr_stat) debut = 0;
       	else debut = (long int) (rand()*kdb); /* Setting of the beginning of the subsequence.  */
/*	end=clock();
	cpu_time=((double)(end-start))/CLOCKS_PER_SEC;
	printf("CPU time: %.6f s\n",cpu_time);*/
       	printf("Beginning at sample %ld\n",debut);

	printf("\nFor storing:\n");
	printf("\t\t- the time errors          \t\t->\t\t0\n");
	printf("\t\t- the frequency deviations \t\t->\t\t1\n");
	printf("\t\t- both                     \t\t->\t\t2\n");
	do
	    scanf("%c",&rep);
	while((rep<48)||(rep>50));
	if (rep!='0')
	    {
	    printf("\nCreation of the file %s \n",nykfich);
	    fiptyk = fopen(nykfich,"w");
	    fprintf(fiptyk,"%% %s\n",nykfich);
	    fprintf(fiptyk,"%% Sequence of simulated frequency deviation samples Yk.\n");
	    fprintf(fiptyk,"%% Inverse of the low cut-off frequency:...... %24ld\n",nbr_dat);
	    fprintf(fiptyk,"%% Data number:............................... %24ld\n",nbr_stat);
	    fprintf(fiptyk,"%% Sampling step:............................. %24.16le\n",tau0);
	    fprintf(fiptyk,"%% h-3  (1/f3 coef.):......................... %24.16le\n",hm3);
	    fprintf(fiptyk,"%% h-2  (1/f2 coef.):......................... %24.16le\n",hm2);
	    fprintf(fiptyk,"%% h-1  (1/f coef.):.......................... %24.16le\n",hm1);
	    fprintf(fiptyk,"%% h0   (constant coef.):..................... %24.16le\n",h0);
	    fprintf(fiptyk,"%% h+1  (f coef.):............................ %24.16le\n",hp1);
	    fprintf(fiptyk,"%% h+2  (f2 coef.):........................... %24.16le\n",hp2);
	    fprintf(fiptyk,"%% C1 (linear frequency drift coef.):......... %24.16le\n",C1);
	    fprintf(fiptyk,"%% C0 (constant frequency drift coef.):....... %24.16le\n",C0);
	    fprintf(fiptyk,"%% Initialization seed of the random sequence: %24ld\n",GR);
	    fprintf(fiptyk,"%% Beginning at sample:....................... %24ld\n",debut);
	    fprintf(fiptyk,"%% ____________________________________________________________________\n");
	    for(i=0;i<nbr_stat;++i)
	        fprintf(fiptyk,"%24.16le \t %24.16le\n",((double)(debut+i))*tau0,x[debut+i]);
	    fclose(fiptyk);
            }
	if (rep!='1')
	    {
	    yk_xt(nbr_dat,0,tau0);
	    printf("\nCreation of the file %s \n",nxtfich);
	    fiptxt = fopen(nxtfich,"w");
	    fprintf(fiptxt,"%% %s\n",nxtfich);
	    fprintf(fiptxt,"%% Sequence of simulated time error samples x(t).\n");
	    fprintf(fiptxt,"%% Inverse of the low cut-off frequency:...... %24ld\n",nbr_dat);
	    fprintf(fiptxt,"%% Data number:............................... %24ld\n",nbr_stat);
	    fprintf(fiptxt,"%% Sampling step:............................. %24.16le\n",tau0);
	    fprintf(fiptxt,"%% h-3  (1/f3 coef.):......................... %24.16le\n",hm3);
	    fprintf(fiptxt,"%% h-2  (1/f2 coef.):......................... %24.16le\n",hm2);
	    fprintf(fiptxt,"%% h-1  (1/f coef.):.......................... %24.16le\n",hm1);
	    fprintf(fiptxt,"%% h0   (constant coef.):..................... %24.16le\n",h0);
	    fprintf(fiptxt,"%% h+1  (f coef.):............................ %24.16le\n",hp1);
	    fprintf(fiptxt,"%% h+2  (f2 coef.):........................... %24.16le\n",hp2);
	    fprintf(fiptxt,"%% C1 (linear frequency drift coef.):......... %24.16le\n",C1);
	    fprintf(fiptxt,"%% C0 (constant frequency drift coef.):....... %24.16le\n",C0);
	    fprintf(fiptxt,"%% Initialization seed of the random sequence: %24ld\n",GR);
	    fprintf(fiptxt,"%% Beginning at sample:....................... %24ld\n",debut);
	    fprintf(fiptxt,"%% ____________________________________________________________________\n");
	    for(i=0;i<nbr_stat;++i)
	        fprintf(fiptxt,"%24.16le \t %24.16le\n",((double)(debut+i))*tau0,x[debut+i]);
	    fclose(fiptxt);
	    }
	free(x);

	return 0;
	}










