/*   stio_sbr.c                                   F. Vernotte - 2010/10/31  */
/*   Input/output subroutines                                               */
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

#define DATAMAX 16384
#define GRANMAX 67108864

extern double *T, *Y, coeff[], ortau[], log_inc;
extern char st_version[];
extern char flag_graph, flag_conf, flag_bias, flag_title, flag_fit, flag_asymptote, flag_slopes[], flag_variance, flag_log_inc;
extern int ntau;

#define db(x) ((double)(x))
#define sisig(x) ( (x) == 0 ) ? (db(0)) : (  ( (x) > 0 ) ? (db(1)) : (db(-1))  )

int init_flag()
/* Initialize flags by reading the configuration file $HOME/.SigmaTheta.conf */
    {
    char *homedir, filepath[512], file[]="/.SigmaTheta.conf", filwin[]=".SigmaTheta.conf", gm[100], *fg, tg[100];
    int tfs, i, indpb;
    long deb_file;
    FILE *ofd;

    homedir=getenv("HOME");
    if (homedir != NULL)
        {
        strcpy(filepath,homedir);
        strcat(filepath,file);
        }
    else
      strcpy(filepath,filwin);
    flag_conf=1;       /* default value: 95%                   */
    flag_bias=1;       /* default value: ON                    */
    flag_graph=1;      /* default value: ON                    */
    flag_title=1;      /* default value: ON                    */
    flag_fit=9;        /* default value: double fit            */
    flag_asymptote=1;  /* default value: ON                    */
    for(i=0;i<6;++i) 
      flag_slopes[i]=1;/* default value: fit all slopes        */
    flag_variance=0;   /* default value: Allan variance        */
    flag_slopes[0]=0;  /* No tau^-3/2 slope for Allan variance */
    ofd=fopen(filepath, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
	do
	    {
	    fg=fgets(gm,100,ofd);
	    if (fg==NULL) return(-2);
	    }
       	while(strcmp(fg,"# Contact: francois.vernotte@obs-besancon.fr\n"));
	deb_file=ftell(ofd);
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[12]=0;
	    }
       	while(strcmp(tg,"95% bounds: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[12]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_conf=0;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[12]=0;
	    }
       	while(strcmp(tg,"68% bounds: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[12]);
	    if (!strcmp(tg,"ON\n"))
	        flag_conf+=8;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[20]=0;
	    }
       	while(strcmp(tg,"Unbiased estimates: ")&&(fg)); 
	if (fg)
	    {
	    strcpy(tg,&gm[20]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_bias=0;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[7]=0;
	    }
       	while(strcmp(tg,"Graph: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[7]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_graph=0;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[7]=0;
	    }
       	while(strcmp(tg,"Title: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[7]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_title=0;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[5]=0;
	    }
       	while(strcmp(tg,"Fit: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[5]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_fit-=1;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[12]=0;
	    }
       	while(strcmp(tg,"Double Fit: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[12]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_fit-=8;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	flag_asymptote=1; /* default value : with asymptotes */
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[12]=0;
	    }
       	while(strcmp(tg,"Asymptotes: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[12]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_asymptote=0;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[16]=0;
	    }
       	while(strcmp(tg,"Tau^-3/2 slope: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[16]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_slopes[0]=0;
	    else
	        flag_slopes[0]=1; /* set to 0 at initialization */
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[14]=0;
	    }
       	while(strcmp(tg,"Tau^-1 slope: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[14]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_slopes[1]=0;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[16]=0;
	    }
       	while(strcmp(tg,"Tau^-1/2 slope: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[16]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_slopes[2]=0;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[13]=0;
	    }
       	while(strcmp(tg,"Tau^0 slope: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[13]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_slopes[3]=0;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[15]=0;
	    }
       	while(strcmp(tg,"Tau^1/2 slope: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[15]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_slopes[4]=0;
	    }
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[13]=0;
	    }
       	while(strcmp(tg,"Tau^1 slope: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[13]);
	    if (!strcmp(tg,"OFF\n"))
	        flag_slopes[5]=0;
	    }
	indpb=0;
	for(i=0;i<6;++i) indpb+=flag_slopes[i];
	if (!indpb) for(i=0;i<6;++i) flag_slopes[i]=1;
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[25]=0;
	    }
       	while(strcmp(tg,"Modified Allan variance: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[25]);
	    if (!strcmp(tg,"ON\n"))
		{
	        flag_variance=1;
		printf("MVAR ON  -> flag_variance=%d flag_slopes[0]=%d\n",flag_variance,flag_slopes[0]);
		}
	    else
		{
		flag_slopes[0]=0;
		printf("MVAR OFF -> flag_variance=%d flag_slopes[0]=%d\n",flag_variance,flag_slopes[0]);
		}
	    }
	flag_log_inc=1;
	log_inc=(double)2; /* default value: tau increment by octave */
	tfs=fseek(ofd,deb_file,SEEK_SET); 
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[5]=0;
	    }
       	while(strcmp(tg,"Tau: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[5]);
            ntau=sscanf(tg,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&ortau[0],&ortau[1],&ortau[2],&ortau[3],&ortau[4],&ortau[5],&ortau[6],&ortau[6],&ortau[8],&ortau[9]);
	    if (ntau>0)
		{
		if (ortau[0]>1) log_inc=ortau[0];
		else if ((ortau[0]==(double)1)&&(ntau>1)) flag_log_inc=0;
		}
	    }
	return(0);
        }
    }

int load_ykt(char *source)
/* Load the file pointed by 'source' and transfer its contain into the global 'T' and 'Y' tables. */
/* Output values: -1 file not found       */
/*                 0 unrecognized file    */
/*                 N length of the tables */
    {
    int i, nbv, N;
    long int dtmx;
    double tst;
    char gm[256];    
    FILE *ofd;

    dtmx=DATAMAX;
    T=(double *)malloc(dtmx*sizeof(double)); 
    Y=(double *)malloc(dtmx*sizeof(double)); 
    ofd=fopen(source, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    fgets(gm,100,ofd);
	while((gm[0]=='#')||(gm[0]=='%'));
        i=0;
        nbv=sscanf(gm,"%lf %lf %lf",&T[i],&Y[i],&tst);
	if (nbv!=2)
	    {
	    if (nbv!=1) nbv=-nbv;
	    return(nbv);
	    }
	else
	    {
            do
                {
                i++;
	        if (i>=dtmx)
	            {
	            dtmx+=DATAMAX;
	            if (dtmx>GRANMAX)
	                {
		        printf("# File trucated to %ld elements\n",dtmx);
                        break;
	                }
	            T=(double *)realloc(T,dtmx*sizeof(double));
	            Y=(double *)realloc(Y,dtmx*sizeof(double));
	            }
	        }
            while(fscanf(ofd,"%lf %lf",&T[i],&Y[i])==2);
            fclose(ofd);
            N=i;
	    }
        }
    return(N);
    }

int load_1col(char *source)
/* Load the file pointed by 'source' and transfer its contain into the global 'Y' table. */
/* Output values: -2 file with 2 columns  */
/*                -1 file not found       */
/*                 0 unrecognized file    */
/*                 N length of the tables */
    {
    int i, nbv, N;
    double tst;
    long int dtmx;
    char gm[256];    
    FILE *ofd;

    dtmx=DATAMAX;
    Y=(double *)malloc(dtmx*sizeof(double)); 
    ofd=fopen(source, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    fgets(gm,100,ofd);
	while((gm[0]=='#')||(gm[0]=='%'));
        i=0;
        nbv=sscanf(gm,"%lf %lf %lf",&Y[i],&tst,&tst);
	if (nbv!=1)
	    {
	    return(-nbv);
	    }
	else
	    {
            do
                {
                i++;
	        if (i>=dtmx)
	            {
	            dtmx+=DATAMAX;
	            if (dtmx>GRANMAX)
	                {
		        printf("# File trucated to %ld elements\n",dtmx);
                        break;
	                }
	            Y=(double *)realloc(Y,dtmx*sizeof(double));
	            }
	        }
            while(fscanf(ofd,"%lf",&Y[i])==1);
            fclose(ofd);
            N=i;
	    }
        }
    return(N);
    }

int load_adev(char *source, double tau[], double adev[])
/* Load the file pointed by 'source' and transfer its contain into the 'tau' and 'adev' tables. */
/* Output values: -1 file not found       */
/*                 0 unrecognized file    */
/*                 N length of the tables */
    {
    int i, nbv, N;
    double tst;
    char gm[256];    
    FILE *ofd;

    ofd=fopen(source, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    fgets(gm,100,ofd);
	while((gm[0]=='#')||(gm[0]=='%'));
        i=0;
        nbv=sscanf(gm,"%lf %lf %lf",&tau[i],&adev[i],&tst);
	if (nbv!=2)
	    {
	    if (nbv!=1) nbv=-nbv;
	    return(nbv);
	    }
	else
	    {
            do
		  i++;
            while(fscanf(ofd,"%lf %lf",&tau[i],&adev[i])==2);
            fclose(ofd);
            N=i;
	    }
        }
    return(N);
    }

int load_coef(char *source)
/* Load the file pointed by 'source' and transfer its contain into the 'coeff' global table. */
/* Output values: -1 file not found       */
/*                 0 unrecognized file    */
/*                 N length of the tables */
    {
    int nbv, N;
    char gm[512],tst[256];    
    FILE *ofd;

    nbv=0;
    ofd=fopen(source, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    fgets(gm,256,ofd);
	while((gm[0]=='#')||(gm[0]=='%'));
        nbv=sscanf(gm,"%lf %lf %lf %lf %lf %lf %s",&coeff[0],&coeff[1],&coeff[2],&coeff[3],&coeff[4],&coeff[5],tst);
	fclose(ofd);
	if (nbv!=6)
	    if (nbv!=1) nbv=-nbv;
        }
    return(nbv);
    }

int load_3col(char *source, double tau[], double adev[], double ubad[])
/* Load the file pointed by 'source' and transfer its contain into the 'tau', 'adev' and 'ubad' tables. */
/* Output values: -1 file not found       */
/*                 0 unrecognized file    */
/*                 N length of the tables */
    {
    int i, nbv, N;
    char tst[512];
    char gm[256];    
    FILE *ofd;

    ofd=fopen(source, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    fgets(gm,256,ofd);
	while((gm[0]=='#')||(gm[0]=='%'));
        i=0;
        nbv=sscanf(gm,"%lf %lf %lf %s",&tau[i],&adev[i],&ubad[i],tst);
	if (nbv<2)
	    {
	    if (nbv!=1) nbv=-nbv;
	    return(nbv);
	    }
	else
	    {
            do
	        {
		if (nbv==2) ubad[i]=-1;
                nbv=sscanf(gm,"%lf %lf %lf %s",&tau[i],&adev[i],&ubad[i],tst);
		i++;
	        }
            while(fgets(gm,256,ofd)!=NULL);
            fclose(ofd);
            N=i;
	    }
        }
    return(N);
    }

int load_7col(char *source, double tau[], double adev[], double ubad[], double b1[], double b2[], double b3[], double b4[])
/* Load the file pointed by 'source' and transfer its contain into the 'tau', 'adev', 'ubad' and bounds (b1, b2, b3, b4) tables. */
/* Output values: -1 file not found       */
/*                 0 unrecognized file    */
/*                 N length of the tables */
    {
    int i, nbv, N;
    char tst[512];
    char gm[256];    
    FILE *ofd;

    ofd=fopen(source, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    fgets(gm,256,ofd);
	while((gm[0]=='#')||(gm[0]=='%'));
        i=0;
        nbv=sscanf(gm,"%lf %lf %lf %lf %lf %lf %lf",&tau[i],&adev[i],&ubad[i],&b1[i],&b2[i],&b3[i],&b4[i]);
	if (nbv<7)
	    {
	    if (nbv!=1) nbv=-nbv;
	    return(nbv);
	    }
	else
	    {
            do
	        {
		nbv=sscanf(gm,"%lf %lf %lf %lf %lf %lf %lf",&tau[i],&adev[i],&ubad[i],&b1[i],&b2[i],&b3[i],&b4[i]);
		i++;
	        }
            while(fgets(gm,256,ofd)!=NULL);
            fclose(ofd);
            N=i;
	    }
        }
    return(N);
    }

int gener_gplt(char *outfile, int N, double tau[], double adev[], double bmax[])
/* Generate a gnuplot file (.gnu) and invoke gnuplot for creating a postscript file */
    {
    int i,mii,mxi,err;
    double minx, maxx, miny, maxy, lmix, lmax, lmiy, lmay, ltx, lty, lmx, lmy, rtmx, rtmy;
    char gptfile[256], psfile[256], gpt_cmd[65536], sys_cmd[256];
    FILE *ofd;

    strcpy(gptfile,outfile);
    strcat(gptfile,".gnu");
    strcpy(psfile,outfile);
    strcat(psfile,".ps");
    ofd=fopen(gptfile, "w");
    if (ofd==NULL) return(-1);
    fprintf(ofd,"set terminal postscript landscape enhanced color solid \"Helvetica\" 18\n");
    fprintf(ofd,"set output \"%s\"\n",psfile);
    fprintf(ofd,"set logscale xy\n");
    fprintf(ofd,"set format xy \"%%.0e\"\n");
    fprintf(ofd,"set grid\n");
    minx=miny=1e99;
    maxx=maxy=db(0);
    for(i=0;i<N;++i)
	{
	if (tau[i]<minx) minx=tau[i];
	if (tau[i]>maxx) maxx=tau[i];
	if (adev[i]<miny)
	    {
	    mii=i;
	    miny=adev[i];
	    }
	if (bmax[i]>maxy)
	    {
	    mxi=i;
	    maxy=bmax[i];
	    }
	}
    if (mxi<mii)
	fprintf(ofd,"set key right\n");
    else
	fprintf(ofd,"set key left\n");
    lmix=log(minx);
    lmax=log(maxx);
    lmiy=log(miny);
    lmay=log(maxy);
    ltx=db(0.05)*(lmax-lmix);
    lmix=lmix-ltx;
    lmax=lmax+ltx;
    lty=db(0.1)*(lmay-lmiy);
    lmiy=lmiy-lty;
    lmay=lmay+lty;
    lmx=lmax-db(0.15)*ltx;
    lmy=lmiy+db(0.25)*lty;
    rtmx=exp(lmx);
    rtmy=exp(lmy);
    minx=exp(lmix);
    maxx=exp(lmax);
    miny=exp(lmiy);
    maxy=exp(lmay);
    fprintf(ofd,"set xrange[%9.3e:%9.3e]\n",minx,maxx);
    fprintf(ofd,"set yrange[%9.3e:%9.3e]\n",miny,maxy);
    if (flag_title) fprintf(ofd,"set title \"%s\"\n",outfile);
    fprintf(ofd,"set xlabel \"Integration time {/Symbol t} [s]\"\n");
    fprintf(ofd,"set ylabel \"");
    switch(flag_variance)
	{
	case 1 : 
		fprintf(ofd,"M");
		break;

	case 2 :
		fprintf(ofd,"H");
		break;

	case 3 :
		fprintf(ofd,"P");
		break;

	default :
    		fprintf(ofd,"A");
	}
    fprintf(ofd,"DEV {/Symbol s}_");
    switch(flag_variance)
	{
	case 1 : 
		fprintf(ofd,"M");
		break;

	case 2 :
		fprintf(ofd,"H");
		break;

	case 3 :
		fprintf(ofd,"P");
		break;

	default :
    		fprintf(ofd,"A");
	}
    fprintf(ofd,"({/Symbol t})\"\n");
    fprintf(ofd,"set style line 1 pt 2 lc 7 lw 3\n");
    fprintf(ofd,"set style line 2 pt 2 lc 7 lw 2\n");
    fprintf(ofd,"set style line 3 pt 6 lc 2 lw 3\n");
    fprintf(ofd,"set label \"Sigma Theta %s\" at %9.3e,%9.3e right font \"Helvetica,10\"\n",st_version,rtmx,rtmy);
    fprintf(ofd,"plot ");
    switch(flag_conf)
        {
	case 0 :
            fprintf(ofd,"\"%s\" using 1:2 notitle with points ls 2 ",outfile);
	    break;
	case 1 :
            fprintf(ofd,"\"%s\" using 1:2:4:7 title \"95 %% confidence interval\" with yerrorbars ls 2 ",outfile);
	    break;
	case 8 : 
            fprintf(ofd,"\"%s\" using 1:2:5:6 title \"68 %% confidence interval\" with yerrorbars ls 1 ",outfile);
	    break;
	case 9 :
            fprintf(ofd,"\"%s\" using 1:2:4:7 title \"95 %% confidence interval\" with yerrorbars ls 2 ",outfile);
            fprintf(ofd,", \"%s\" using 1:2:5:6 title \"68 %% confidence interval\" with yerrorbars ls 1 ",outfile);
	    break;
	default :
            fprintf(ofd,"\"%s\" using 1:2:4:7 title \"95 %% confidence interval\" with yerrorbars ls 2 ",outfile);
	}
    if (flag_bias)
        fprintf(ofd,", \"%s\" using 1:3 title \"unbiased estimates\" with points ls 3",outfile);
    if (flag_fit)
        {
        fprintf(ofd,", sqrt(%12.6e/x**3+",coeff[0]);
        fprintf(ofd,"%12.6e/x**2+",coeff[1]);
        fprintf(ofd,"%12.6e/x+",coeff[2]);
        fprintf(ofd,"%12.6e+",coeff[3]);
        fprintf(ofd,"%12.6e*x+",coeff[4]);
        fprintf(ofd,"%12.6e*x**2) notitle with line lt 1 lw 2",coeff[5]);
        }
    if (flag_asymptote)
        {
        if (coeff[0]!=0)
            fprintf(ofd,", sqrt(%12.6e/x**3) title \"%7.1e {/Symbol t}^{-3/2} \" with line lt 7",coeff[0],sqrt(coeff[0]));
        if (coeff[1]!=0)
            fprintf(ofd,", sqrt(%12.6e/x**2) title \"%7.1e {/Symbol t}^{-1} \" with line lt 3",coeff[1],sqrt(coeff[1]));
        if (coeff[2]!=0)
            fprintf(ofd,", sqrt(%12.6e/x) title \"%7.1e {/Symbol t}^{-1/2}\" with line lt 4",coeff[2],sqrt(coeff[2]));
        if (coeff[3]!=0)
	    fprintf(ofd,", sqrt(%12.6e) title \"%7.1e      \" with line lt 5",coeff[3],sqrt(coeff[3]));
        if (coeff[4]!=0)
	    fprintf(ofd,", sqrt(%12.6e*x) title \"%7.1e {/Symbol t}^{1/2} \" with line lt 6",coeff[4],sqrt(coeff[4]));
        if (coeff[5]!=0)
	    fprintf(ofd,", sqrt(%12.6e*x**2) title \"%7.1e {/Symbol t}    \" with line lt 2",coeff[5],sqrt(coeff[5]));
	}
    fprintf(ofd,"\n");
    fprintf(ofd,"exit\n");
    fclose(ofd);
    strcpy(sys_cmd,"gnuplot ");
    strcat(sys_cmd,gptfile);
    err=system(sys_cmd);
    return(err);
    }

int gen_psdplt(char *outfile, int N, double freq[], double syf[])
/* Generate a gnuplot file (.gnu) and invoke gnuplot for creating a postscript file */
    {
    int i,mii,mxi,err;
    double minx, maxx, miny, maxy, lmix, lmax, lmiy, lmay, ltx, lty, lmx, lmy, rtmx, rtmy;
    char gptfile[256], psfile[256], gpt_cmd[65536], sys_cmd[256];
    FILE *ofd;

    strcpy(gptfile,outfile);
    strcat(gptfile,".gnu");
    strcpy(psfile,outfile);
    strcat(psfile,".ps");
    ofd=fopen(gptfile, "w");
    if (ofd==NULL) return(-1);
    fprintf(ofd,"set terminal postscript landscape enhanced color solid \"Helvetica\" 18\n");
    fprintf(ofd,"set output \"%s\"\n",psfile);
    fprintf(ofd,"set logscale xy\n");
    fprintf(ofd,"set format xy \"%%.0e\"\n");
    fprintf(ofd,"set grid\n");
    minx=miny=1e99;
    maxx=maxy=db(0);
    for(i=0;i<N;++i)
	{
	if (freq[i]<minx) minx=freq[i];
	if (freq[i]>maxx) maxx=freq[i];
	if (syf[i]<miny)
	    {
	    mii=i;
	    miny=syf[i];
	    }
	if (syf[i]>maxy)
	    {
	    mxi=i;
	    maxy=syf[i];
	    }
	}
    if (mxi<mii)
	fprintf(ofd,"set key right\n");
    else
	fprintf(ofd,"set key left\n");
    lmix=log(minx);
    lmax=log(maxx);
    lmiy=log(miny);
    lmay=log(maxy);
    ltx=db(0.05)*(lmax-lmix);
    lmix=lmix-ltx;
    lmax=lmax+ltx;
    lty=db(0.1)*(lmay-lmiy);
    lmiy=lmiy-lty;
    lmay=lmay+lty;
    lmx=lmax-db(0.15)*ltx;
    lmy=lmiy+db(0.25)*lty;
    rtmx=exp(lmx);
    rtmy=exp(lmy);
    minx=exp(lmix);
    maxx=exp(lmax);
    miny=exp(lmiy);
    maxy=exp(lmay);
    fprintf(ofd,"set xrange[%9.3e:%9.3e]\n",minx,maxx);
    fprintf(ofd,"set yrange[%9.3e:%9.3e]\n",miny,maxy);
    if (flag_title) fprintf(ofd,"set title \"%s\"\n",outfile);
    fprintf(ofd,"set xlabel \"Frequency f [Hz]\"\n");
    fprintf(ofd,"set ylabel \"PSD S_y(f)\"\n");
    fprintf(ofd,"set style line 1 pt 6 lc rgb \"#308010\" lw 3\n");
    fprintf(ofd,"set label \"Sigma Theta %s\" at %9.3e,%9.3e right font \"Helvetica,10\"\n",st_version,rtmx,rtmy);
    fprintf(ofd,"plot ");
    fprintf(ofd,"\"%s\" using 1:2 notitle with lines ls 1\n",outfile);
    fprintf(ofd,"exit\n");
    fclose(ofd);
    strcpy(sys_cmd,"gnuplot ");
    strcat(sys_cmd,gptfile);
    err=system(sys_cmd);
    return(err);
    }

int gen_linplt(char *outfile, int N, double tt[], double xy[], int xory)
/* Generate a gnuplot file (.gnu) and invoke gnuplot for creating a postscript file */
    {
    int i,mii,mxi,err;
    double minx, maxx, miny, maxy, lmix, lmax, lmiy, lmay, ltx, lty, lmx, lmy, rtmx, rtmy;
    char gptfile[256], psfile[256], gpt_cmd[65536], sys_cmd[256];
    FILE *ofd;

    strcpy(gptfile,outfile);
    strcat(gptfile,".gnu");
    strcpy(psfile,outfile);
    strcat(psfile,".ps");
    ofd=fopen(gptfile, "w");
    if (ofd==NULL) return(-1);
    fprintf(ofd,"set terminal postscript landscape enhanced color solid \"Helvetica\" 18\n");
    fprintf(ofd,"set output \"%s\"\n",psfile);
    fprintf(ofd,"set format xy \"%%g\"\n");
    fprintf(ofd,"set grid\n");
    minx=maxx=tt[0];
    miny=maxy=xy[0];
    for(i=0;i<N;++i)
	{
	if (tt[i]<minx) minx=tt[i];
	if (tt[i]>maxx) maxx=tt[i];
	if (xy[i]<miny)
	    {
	    mii=i;
	    miny=xy[i];
	    }
	if (xy[i]>maxy)
	    {
	    mxi=i;
	    maxy=xy[i];
	    }
	}
    if (mxi<mii)
	fprintf(ofd,"set key right\n");
    else
	fprintf(ofd,"set key left\n");
    ltx=db(0.05)*(maxx-minx);
    lmix=minx-ltx;
    lmax=maxx+ltx;
    lty=db(0.1)*(maxy-miny);
    lmiy=miny-lty;
    lmay=maxy+lty;
    minx=lmix;
    maxx=lmax;
    miny=lmiy;
    maxy=lmay;
    rtmx=lmax-db(0.15)*ltx;
    rtmy=lmiy+db(0.25)*lty;
    fprintf(ofd,"set xrange[%9.3e:%9.3e]\n",minx,maxx);
    fprintf(ofd,"set yrange[%9.3e:%9.3e]\n",miny,maxy);
    if (flag_title) fprintf(ofd,"set title \"%s\"\n",outfile);
    fprintf(ofd,"set xlabel \"Time t [s]\"\n");
    fprintf(ofd,"set ylabel ");
    if (xory) fprintf(ofd,"\"Frequency deviation Y_k\"\n");
    else fprintf(ofd,"\"Time error x(t) [s]\"\n");
    fprintf(ofd,"set style line 1 pt 6 lc rgb ");
    if (xory) fprintf(ofd,"\"#D01000\"");
    else fprintf(ofd,"\"#0010D0\"");
    fprintf(ofd," lw 3\n");
    fprintf(ofd,"set label \"Sigma Theta %s\" at %9.3e,%9.3e right font \"Helvetica,10\"\n",st_version,rtmx,rtmy);
    fprintf(ofd,"plot ");
    fprintf(ofd,"\"%s\" using 1:2 notitle with lines ls 1\n",outfile);
    fprintf(ofd,"exit\n");
    fclose(ofd);
    strcpy(sys_cmd,"gnuplot ");
    strcat(sys_cmd,gptfile);
    err=system(sys_cmd);
    return(err);
    }

