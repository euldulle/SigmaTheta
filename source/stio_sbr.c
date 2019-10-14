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
#include <stdint.h>
#include <math.h>

#define DATAMAX 16384
#define GRANMAX 67108864
#define MAXLINELENGTH 256
#define MAXCHAR 512


extern double *T, *Y, *Y1, *Y2, *Y12, *Y23, *Y31, coeff[], ortau[], log_inc;
extern char st_version[];
extern char flag_graph, flag_conf, flag_bias, flag_title, flag_fit, flag_asymptote, flag_slopes[], flag_variance, flag_log_inc, flag_display;
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
//		printf("MVAR ON  -> flag_variance=%d flag_slopes[0]=%d\n",flag_variance,flag_slopes[0]);
		}
	    else
		{
		flag_slopes[0]=0;
//		printf("MVAR OFF -> flag_variance=%d flag_slopes[0]=%d\n",flag_variance,flag_slopes[0]);
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
	flag_display=1;
	tfs=fseek(ofd,deb_file,SEEK_SET);
	do
	    {
	    fg=fgets(gm,100,ofd);
	    strcpy(tg,gm);
	    tg[8]=0;
	    }
       	while(strcmp(tg,"Output: ")&&(fg));
	if (fg)
	    {
	    strcpy(tg,&gm[8]);
	    if (!strncmp(tg,"PDF",3)) flag_display=0;
	    }
	return(0);
        }
    }

double scale(char code){
    //
    // returns a scaling factor as a double defined according to
    // the following conventions :
    //
    // code char :
    //  
    //
    //
    
    double factor;


    switch(code){
        case 'd': 
            //
            // d for days, 86400 s ; useful for timestamps in MJDs
            //
            factor=86400.;
            break;
        case 'H':
            //
            // H for hours, 3600 s ; useful for timestamps in hours
            //
            factor=3600.;
            break;
        case 'M':
            //
            // M for minutes, 60 s ; useful for timestamps in minutes
            //
            factor=60.;
            break;
        case 'm':
            //
            // m for milliseconds ; useful for data in milliseconds
            //
            factor=1.e-3;
            break;
        case 'u':
            //
            // u for microseconds ; useful for data in microseconds
            //
            factor=1.e-6;
            break;
        case 'n':
            //
            // n for nanoseconds ; useful for data in nanoseconds
            //
            factor=1.e-9;
            break;
            //
        case 'p':
            //
            // p for picoseconds ; useful for data in picoseconds
            factor=1.e-12;
            break;
        case 'f':
            //
            // f for femtoseconds ; useful for data in femtoseconds
            //
            factor=1.e-15;
            break;
        case 'a':
            //
            // a for attoseconds ; useful for data in attoseconds. Why not ?..
            //
            factor=1.e-18;
            break;
        default:
            // default: no scaling. 
            //    No scaling is normally dealt with code = 0 
            //    so no multiplication by 1 of a full dataset 
            //    should ever occur under normal use.
            factor=1;
        }
    // fprintf(stderr,"# scale: factor %le input @%c@%d@\n", factor, code, code );
    return factor;
    }


int load_ykt(char *source, char scalex, char scaley)
/* Load the file pointed by 'source' and transfer its contain into the global 'T' and 'Y' tables. */
/* Output values: -1 file not found       */
/*                 0 unrecognized file    */
/*                 N length of the tables */
//
// file reading switched from fscanf to readline and sscanf 
// reading only checks that there are at least 2 readable columns
// the existence of additional columns is not checked, and hence 
// no longer considered a fatal error
//
// Scaling can be done using 2 chars scalex scaley
// to deal with timestamps available in MJD, hours or minutes
// and data available in ms, us, ns, ps, fs
//
// if those are 0 (decimal value 0, not char '0'), no scaling 
// applies ; see double scale(char) for details.
//
{
    int i, nbv, N, linecount=0, valcount=0;
    char *rep;
    long int dtmx;
    double tst, xscale=1, yscale=1;
    size_t n;
    //char line[MAXLINELENGTH+1];    
    char *line=NULL;
    FILE *ofd;

    dtmx=DATAMAX;
    T=(double *)malloc(dtmx*sizeof(double)); 
    Y=(double *)malloc(dtmx*sizeof(double)); 
    if (scalex!=0) xscale=scale(scalex);
    if (scaley!=0) yscale=scale(scaley);
    // fprintf(stderr,"# stio_sbr: xscale %le yscale %le\n", xscale, yscale);

    if (strlen(source)==0)
        ofd=stdin;
    else
        ofd=fopen(source, "r");

    if (ofd==NULL){
        fprintf(stderr,"Could not open file %s", source);
        return(-1);
    }
    else {
        i=0;
        while(getline(&line, &n, ofd)!=-1) {
            //
            // read an entire line
            //
            linecount++;
            if (line[0] == '#' || line[0] == '%')
                //
                // ignore comment lines
                //
                continue;
            if ((valcount=sscanf(line,"%lf %lf",&T[i],&Y[i])==2)){
                // reads 2 values out of line

                //
                // proceed to rescaling if needed (scalex != 0 or scaley != 0)
                //
                if (scalex) T[i]*=xscale;
                if (scaley) Y[i]*=yscale;
                i++;
                if (i>=dtmx)
                {
                    dtmx+=DATAMAX;
                    if (dtmx>GRANMAX)
                    {
                        printf("# File truncated to %ld elements\n",dtmx);
                        break;
                    }
                    T=(double *)realloc(T,dtmx*sizeof(double));
                    Y=(double *)realloc(Y,dtmx*sizeof(double));
                }
            }
            else{
                fprintf(stderr," # ignoring line %d : only %d value(s) read #%s# \n", linecount, valcount, line);
            }
        }

        if (ofd!=stdin)
            fclose(ofd);
        N=i;
    }
    return(N);
}

int load_2yk(char *source1, char *source2)
/* Load the file pointed by 'source' and transfer its content into the global 'T' and 'Y' tables. */
/* Output values: -1 file not found       */
/*                 0 unrecognized file    */
/*                 N length of the tables */
    {
    int i, nbv, n1, n2, N;
    char *rep;
    long int dtmx;
    double tst;
    char gm[256];    
    FILE *ofd, *of2;

    dtmx=DATAMAX;
    T=(double *)malloc(dtmx*sizeof(double)); 
    Y1=(double *)malloc(dtmx*sizeof(double)); 
    Y2=(double *)malloc(dtmx*sizeof(double));
/* First file */ 
    ofd=fopen(source1, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    rep=fgets(gm,100,ofd);
	while((gm[0]=='#')||(gm[0]=='%'));
        i=0;
        nbv=sscanf(gm,"%lf %lf %lf",&T[i],&Y1[i],&tst);
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
	            Y1=(double *)realloc(Y1,dtmx*sizeof(double));
	            }
	        }
            while(fscanf(ofd,"%lf %lf",&T[i],&Y1[i])==2);
            fclose(ofd);
            n1=i;
	    }
        }
/* Second file */ 
    dtmx=DATAMAX;
    of2=fopen(source2, "r");
    if (of2==NULL)
        return(-1);
    else
        {
       	do
	    rep=fgets(gm,100,of2);
	while((gm[0]=='#')||(gm[0]=='%'));
        i=0;
        nbv=sscanf(gm,"%lf %lf %lf",&T[i],&Y2[i],&tst);
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
//	            T=(double *)realloc(T,dtmx*sizeof(double));
	            Y2=(double *)realloc(Y2,dtmx*sizeof(double));
	            }
	        }
            while(fscanf(of2,"%lf %lf",&T[i],&Y2[i])==2);
            fclose(of2);
            n2=i;
	    }
        }
    if (n1==n2) N=n1;
    else N=-2;
    return(N);
    }

int load_3yk(char *source1, char *source2, char *source3)
/* Load the file pointed by 'source' and transfer its contain into the global 'T' and 'Y' tables. */
/* Output values: -1 file not found       */
/*                 0 unrecognized file    */
/*                 N length of the tables */
    {
    int i, nbv, n1, n2, n3, N;
    char *rep;
    long int dtmx;
    double tst;
    char gm[256];    
    FILE *ofd, *of2, *of3;

    dtmx=DATAMAX;
    T=(double *)malloc(dtmx*sizeof(double)); 
    Y12=(double *)malloc(dtmx*sizeof(double)); 
    Y23=(double *)malloc(dtmx*sizeof(double)); 
    Y31=(double *)malloc(dtmx*sizeof(double));
/* First file */ 
    ofd=fopen(source1, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    rep=fgets(gm,100,ofd);
	while((gm[0]=='#')||(gm[0]=='%'));
        i=0;
        nbv=sscanf(gm,"%lf %lf %lf",&T[i],&Y12[i],&tst);
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
	            Y12=(double *)realloc(Y12,dtmx*sizeof(double));
	            }
	        }
            while(fscanf(ofd,"%lf %lf",&T[i],&Y12[i])==2);
            fclose(ofd);
            n1=i;
	    }
        }
/* Second file */ 
    dtmx=DATAMAX;
    of2=fopen(source2, "r");
    if (of2==NULL)
        return(-1);
    else
        {
       	do
	    rep=fgets(gm,100,of2);
	while((gm[0]=='#')||(gm[0]=='%'));
        i=0;
        nbv=sscanf(gm,"%lf %lf %lf",&T[i],&Y23[i],&tst);
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
//	            T=(double *)realloc(T,dtmx*sizeof(double));
	            Y23=(double *)realloc(Y23,dtmx*sizeof(double));
	            }
	        }
            while(fscanf(of2,"%lf %lf",&T[i],&Y23[i])==2);
            fclose(of2);
            n2=i;
	    }
        }
    if (n1==n2) N=n1;
    else
	{ 
	N=-2;
	return(N);
	}
/* Third file */ 
    dtmx=DATAMAX;
    of3=fopen(source3, "r");
    if (of2==NULL)
        return(-1);
    else
        {
       	do
	    rep=fgets(gm,100,of2);
	while((gm[0]=='#')||(gm[0]=='%'));
        i=0;
        nbv=sscanf(gm,"%lf %lf %lf",&T[i],&Y31[i],&tst);
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
//	            T=(double *)realloc(T,dtmx*sizeof(double));
	            Y31=(double *)realloc(Y31,dtmx*sizeof(double));
	            }
	        }
            while(fscanf(of2,"%lf %lf",&T[i],&Y31[i])==2);
            fclose(of2);
            n3=i;
	    }
        }
    if (n3!=N) N=-2;
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
    char *rep;
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
	    rep=fgets(gm,100,ofd);
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
    char *rep;
    double tst;
    char gm[256];    
    FILE *ofd;

    if (strlen(source)==0)
        ofd=stdin;
    else
        ofd=fopen(source, "r");
    
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    rep=fgets(gm,100,ofd);
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
    char gm[512],tst[256], *rep;    
    FILE *ofd;

    nbv=0;
    ofd=fopen(source, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    rep=fgets(gm,256,ofd);
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
    char tst[512], *rep;
    char gm[256];    
    FILE *ofd;

    ofd=fopen(source, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    rep=fgets(gm,256,ofd);
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
    char tst[512], *rep;
    char gm[256];    
    FILE *ofd;

    ofd=fopen(source, "r");
    if (ofd==NULL)
        return(-1);
    else
        {
       	do
	    rep=fgets(gm,256,ofd);
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

int gener_gplt(char *outfile, int N, double tau[], double adev[], double bmax[], char *est_typ)
/* Generate a gnuplot file (.gnu) and invoke gnuplot for creating a postscript file */
    {
    int i,mii,mxi,err;
    double minx, maxx, miny, maxy, lmix, lmax, lmiy, lmay, ltx, lty, lmx, lmy, rtmx, rtmy;
    char gptfile[256], psfile[256], gpt_cmd[65536], sys_cmd[256];
    FILE *ofd;

    strcpy(gptfile,outfile);
    strcat(gptfile,".gnu");
    strcpy(psfile,outfile);
    strcat(psfile,".pdf");
    ofd=fopen(gptfile, "w");
    if (ofd==NULL) return(-1);
//    fprintf(ofd,"set terminal postscript landscape enhanced color solid \"Helvetica\" 18\n");
    if (flag_display)
	    fprintf(ofd,"set terminal wxt size 1024,768 enhanced font \"Helvetica\" fontscale 1.5 persist\n");
    else
	{
	fprintf(ofd,"set terminal pdfcairo size 172,128 enhanced color font \"Helvetica\" fontscale 18\n");
	fprintf(ofd,"set output \"%s\"\n",psfile);
	}
    fprintf(ofd,"set logscale xy\n");
    fprintf(ofd,"set format xy \"10^{%%+T}\"\n");
    fprintf(ofd,"set grid\n");
    minx=miny=1e99;
    maxx=maxy=db(0);
    for(i=0;i<N;++i)
	{
	if (tau[i]<minx) minx=tau[i];
	if (tau[i]>maxx) maxx=tau[i];
	if (fabs(adev[i])<miny)
	    {
	    mii=i;
	    miny=fabs(adev[i]);
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
    fprintf(ofd,"set mxtics 10\n");
    fprintf(ofd,"set mytics 10\n");
    if (flag_title) fprintf(ofd,"set title \"%s\" noenhanced\n",outfile);
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
    if (flag_display)
	{
	fprintf(ofd,"set style line 1 pt 2 ps 1 lc 7 lw 3\n");
	fprintf(ofd,"set style line 2 pt 2 ps 1 lc 7 lw 2\n");
	fprintf(ofd,"set style line 3 pt 6 ps 2 lc rgb \"#30D015\" lw 4\n");
	fprintf(ofd,"set style line 4 lc rgb \"#D01000\" lw 5\n");
	fprintf(ofd,"set style line 5 lc rgb \"#00A0A0\" lw 3\n");
	fprintf(ofd,"set style line 6 lc rgb \"#FFE000\" lw 3\n");
	fprintf(ofd,"set style line 7 lc rgb \"#109010\" lw 3\n");
	fprintf(ofd,"set style line 8 lc rgb \"#A000A0\" lw 3\n");
	fprintf(ofd,"set style line 9 lc rgb \"#0010D0\" lw 3\n");
	fprintf(ofd,"set style line 10 lc rgb \"#FF8000\" lw 3\n");
	fprintf(ofd,"set style line 11 pt 0 ps 1 lc 7 lw 3\n");
	fprintf(ofd,"set style line 12 pt 0 ps 1 lc rgb \"#A0A0A0\" lw 2\n");
	fprintf(ofd,"set label \"SigmaTheta %s\" at %9.3e,%9.3e right font \"Helvetica,6\"\n",st_version,rtmx,rtmy);
	}
    else
	{
	fprintf(ofd,"set style line 1 pt 2 ps 20 lc 7 lw 60\n");
	fprintf(ofd,"set style line 2 pt 2 ps 20 lc 7 lw 40\n");
	fprintf(ofd,"set style line 3 pt 6 ps 20 lc rgb \"#30D015\" lw 80\n");
	fprintf(ofd,"set style line 4 lc rgb \"#D01000\" lw 80\n");
	fprintf(ofd,"set style line 5 lc rgb \"#00C0FF\" lw 60\n");
	fprintf(ofd,"set style line 6 lc rgb \"#FFE000\" lw 60\n");
	fprintf(ofd,"set style line 7 lc rgb \"#109010\" lw 60\n");
	fprintf(ofd,"set style line 8 lc rgb \"#A000A0\" lw 60\n");
	fprintf(ofd,"set style line 9 lc rgb \"#0010D0\" lw 60\n");
	fprintf(ofd,"set style line 10 lc rgb \"#FF8000\" lw 60\n");
	fprintf(ofd,"set style line 11 pt 0 ps 20 lc 7 lw 60\n");
	fprintf(ofd,"set style line 12 pt 0 ps 20 lc rgb \"#A0A0A0\" lw 40\n");
	fprintf(ofd,"set label \"SigmaTheta %s\" at %9.3e,%9.3e right font \"Helvetica,8\"\n",st_version,rtmx,rtmy);
	}
    fprintf(ofd,"plot ");
    switch(flag_conf)
        {
	case 0 :
            fprintf(ofd,"\"%s\" using 1:2 notitle with points ls 2 ",outfile);
	    break;
	case 1 :
            fprintf(ofd,"\"%s\" using 1:2 notitle with points ls 2, ",outfile);
            fprintf(ofd,"\"%s\" using 1:3:4:7 title \"95 %% confidence interval\" with yerrorbars ls 12 ",outfile);
	    break;
	case 8 : 
            fprintf(ofd,"\"%s\" using 1:2 notitle with points ls 2, ",outfile);
            fprintf(ofd,"\"%s\" using 1:3:5:6 title \"68 %% confidence interval\" with yerrorbars ls 11 ",outfile);
	    break;
	case 9 :
            fprintf(ofd,"\"%s\" using 1:2 notitle with points ls 2, ",outfile);
            fprintf(ofd,"\"%s\" using 1:3:4:7 title \"95 %% confidence interval\" with yerrorbars ls 12 ",outfile);
            fprintf(ofd,", \"%s\" using 1:3:5:6 title \"68 %% confidence interval\" with yerrorbars ls 11 ",outfile);
	    break;
	default :
            fprintf(ofd,"\"%s\" using 1:2 notitle with points ls 2, ",outfile);
            fprintf(ofd,"\"%s\" using 1:3:4:7 title \"95 %% confidence interval\" with yerrorbars ls 12 ",outfile);
	}
    if (flag_bias)
        fprintf(ofd,", \"%s\" using 1:3 title \"%s estimates\" with points ls 3",outfile, est_typ);
    if (flag_fit)
        {
        fprintf(ofd,", sqrt(%12.6e/x**3+",coeff[0]);
        fprintf(ofd,"%12.6e/x**2+",coeff[1]);
        fprintf(ofd,"%12.6e/x+",coeff[2]);
        fprintf(ofd,"%12.6e+",coeff[3]);
        fprintf(ofd,"%12.6e*x+",coeff[4]);
        fprintf(ofd,"%12.6e*x**2) notitle with line ls 4",coeff[5]);
        }
    if (flag_asymptote)
        {
        if (coeff[0]!=0)
            fprintf(ofd,", sqrt(%12.6e/x**3) title \"%7.1e {/Symbol t}^{-3/2} \" with line ls 5",coeff[0],sqrt(coeff[0]));
        if (coeff[1]!=0)
            fprintf(ofd,", sqrt(%12.6e/x**2) title \"%7.1e {/Symbol t}^{-1} \" with line ls 6",coeff[1],sqrt(coeff[1]));
        if (coeff[2]!=0)
            fprintf(ofd,", sqrt(%12.6e/x) title \"%7.1e {/Symbol t}^{-1/2}\" with line ls 7",coeff[2],sqrt(coeff[2]));
        if (coeff[3]!=0)
	    fprintf(ofd,", sqrt(%12.6e) title \"%7.1e      \" with line ls 8",coeff[3],sqrt(coeff[3]));
        if (coeff[4]!=0)
	    fprintf(ofd,", sqrt(%12.6e*x) title \"%7.1e {/Symbol t}^{1/2} \" with line ls 9",coeff[4],sqrt(coeff[4]));
        if (coeff[5]!=0)
	    fprintf(ofd,", sqrt(%12.6e*x**2) title \"%7.1e {/Symbol t}    \" with line ls 10",coeff[5],sqrt(coeff[5]));
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
    strcat(psfile,".pdf");
    ofd=fopen(gptfile, "w");
    if (ofd==NULL) return(-1);
//    fprintf(ofd,"set terminal postscript landscape enhanced color solid \"Helvetica\" 18\n");
    if (flag_display)
	fprintf(ofd,"set terminal wxt size 1024,768 enhanced font \"Helvetica\" fontscale 1.5 persist\n");
    else
	{
	fprintf(ofd,"set terminal pdfcairo size 172,128 enhanced color font \"Helvetica\" fontscale 18\n");
	fprintf(ofd,"set output \"%s\"\n",psfile);
	}
    fprintf(ofd,"set logscale xy\n");
    fprintf(ofd,"set format xy \"10^{%%+T}\"\n");
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
    fprintf(ofd,"set mxtics 10\n");
    fprintf(ofd,"set mytics 10\n");
    if (flag_title) fprintf(ofd,"set title \"%s\" noenhanced\n",outfile);
    fprintf(ofd,"set xlabel \"Frequency f [Hz]\"\n");
    fprintf(ofd,"set ylabel \"PSD S_y(f)\"\n");
    fprintf(ofd,"set style line 1 pt 6 lc rgb \"#308015\" lw 3\n");
    fprintf(ofd,"set label \"SigmaTheta %s\" at %9.3e,%9.3e right font \"Helvetica,8\"\n",st_version,rtmx,rtmy);
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
    strcat(psfile,".pdf");
    ofd=fopen(gptfile, "w");
    if (ofd==NULL) return(-1);
//    fprintf(ofd,"set terminal postscript landscape enhanced color solid \"Helvetica\" 18\n");
    if (flag_display)
	fprintf(ofd,"set terminal wxt size 1024,768 enhanced font \"Helvetica\" fontscale 1.5 persist\n");
    else
	{
	fprintf(ofd,"set terminal pdfcairo size 172,128 enhanced color font \"Helvetica\" fontscale 18\n");
	fprintf(ofd,"set output \"%s\"\n",psfile);
	}
    fprintf(ofd,"set format xy \"%%g\"\n");
//    fprintf(ofd,"set format xy \"10^{%%+T}\"\n");
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
    if (flag_title) fprintf(ofd,"set title \"%s\" noenhanced\n",outfile);
    fprintf(ofd,"set xlabel \"Time t [s]\"\n");
    fprintf(ofd,"set ylabel ");
    if (xory) fprintf(ofd,"\"Frequency deviation Y_k\"\n");
    else fprintf(ofd,"\"Time error x(t) [s]\"\n");
    fprintf(ofd,"set style line 1 pt 6 lc rgb ");
    if (xory) fprintf(ofd,"\"#D01000\"");
    else fprintf(ofd,"\"#0010D0\"");
    fprintf(ofd," lw 3\n");
    fprintf(ofd,"set label \"SigmaTheta %s\" at %9.3e,%9.3e right font \"Helvetica,8\"\n",st_version,rtmx,rtmy);
    fprintf(ofd,"plot ");
    fprintf(ofd,"\"%s\" using 1:2 notitle with lines ls 1\n",outfile);
    fprintf(ofd,"exit\n");
    fclose(ofd);
    strcpy(sys_cmd,"gnuplot ");
    strcat(sys_cmd,gptfile);
    err=system(sys_cmd);
    return(err);
    }

int gen_gcodplt(char *outfile, char names[][256], int N, int nbf, double tau[], double gcod[][256], int ind_gcod)
/* Generate a gnuplot file (.gnu) and invoke gnuplot for creating a postscript file */
    {
    int i,j,k,mii,mxi,err;
    double minx, maxx, miny, maxy, lmix, lmax, lmiy, lmay, ltx, lty, lmx, lmy, rtmx, rtmy;
    char gptfile[256], psfile[256], gpt_cmd[65536], sys_cmd[256];
    FILE *ofd;

    strcpy(gptfile,outfile);
    strcat(gptfile,".gnu");
    strcpy(psfile,outfile);
    strcat(psfile,".pdf");
    ofd=fopen(gptfile, "w");
    if (ofd==NULL) return(-1);

    if (flag_display)
	fprintf(ofd,"set terminal wxt size 1024,768 enhanced font \"Helvetica\" fontscale 1.5 persist\n");
    else
	{
	fprintf(ofd,"set terminal pdfcairo size 172,128 enhanced color font \"Helvetica\" fontscale 18\n");
	fprintf(ofd,"set output \"%s\"\n",psfile);
	}
    fprintf(ofd,"set logscale xy\n");
    fprintf(ofd,"set format xy \"10^{%%+T}\"\n");
    fprintf(ofd,"set grid\n");
    minx=miny=1e99;
    maxx=maxy=db(0);
    for(i=0;i<N;++i)
	{
	if (tau[i]<minx) minx=tau[i];
	if (tau[i]>maxx) maxx=tau[i];
	}
    for(j=0;j<nbf;++j)
	for(i=0;i<N;++i)
		{
		if (fabs(gcod[j][i])<miny)
	    		{
	    		mii=i;
	    		miny=fabs(gcod[j][i]);
	    		}
		if (fabs(gcod[j][i])>maxy)
	    		{
	    		mxi=i;
	    		maxy=fabs(gcod[j][i]);
	    		}
		}
/*    if (mxi<mii)
	fprintf(ofd,"set key right\n");
    else*/
	fprintf(ofd,"set key left bottom reverse\n");
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
    fprintf(ofd,"set mxtics 10\n");
    fprintf(ofd,"set mytics 10\n");
    if (flag_title) fprintf(ofd,"set title \"%s\" noenhanced\n",outfile);
    fprintf(ofd,"set xlabel \"Integration time {/Symbol t} [s]\"\n");
    fprintf(ofd,"set ylabel \"ADEV {/Symbol s}_A({/Symbol t})\"\n");
    if (flag_display)
	{
	fprintf(ofd,"set style line 1 lt 1 pt 4 ps 1.1 lc rgb \"#D01000\" lw 2\n");
	fprintf(ofd,"set style line 2 lt 1 pt 12 ps 1.6 lc rgb \"#308015\" lw 2\n");
	fprintf(ofd,"set style line 3 lt 1 pt 6 ps 1.3 lc rgb \"#0010D0\" lw 2\n");
	fprintf(ofd,"set style line 4 lt 1 pt 13 ps 1 lc rgb \"#806020\" lw 2\n");
	fprintf(ofd,"set style line 5 lt 2 dt 2 pt 3 ps 1.2 lc rgb \"#D01000\" lw 2\n");
	fprintf(ofd,"set style line 6 lt 2 dt 2 pt 1 ps 1.4 lc rgb \"#308015\" lw 2\n");
	fprintf(ofd,"set style line 7 lt 2 dt 2 pt 2 ps 1.1 lc rgb \"#0010D0\" lw 2\n");
	fprintf(ofd,"set style line 8 lt 2 dt 2 pt 5 ps 1 lc rgb \"#806020\" lw 2\n");
	fprintf(ofd,"set label \"SigmaTheta %s\" at %9.3e,%9.3e right font \"Helvetica,6\"\n",st_version,rtmx,rtmy);
	}
    else
	{
	fprintf(ofd,"set style line 1 lt 1 pt 4 ps 22 lc rgb \"#D01000\" lw 60\n");
	fprintf(ofd,"set style line 2 lt 1 pt 12 ps 32 lc rgb \"#308015\" lw 60\n");
	fprintf(ofd,"set style line 3 lt 1 pt 6 ps 26 lc rgb \"#0010D0\" lw 60\n");
	fprintf(ofd,"set style line 4 lt 1 pt 13 ps 20 lc rgb \"#806020\" lw 60\n");
	fprintf(ofd,"set style line 5 lt 2 dt 2 pt 3 ps 24 lc rgb \"#D01000\" lw 60\n");
	fprintf(ofd,"set style line 6 lt 2 dt 2 pt 1 ps 28 lc rgb \"#308015\" lw 60\n");
	fprintf(ofd,"set style line 7 lt 2 dt 2 pt 2 ps 22 lc rgb \"#0010D0\" lw 60\n");
	fprintf(ofd,"set style line 8 lt 2 dt 2 pt 5 ps 20 lc rgb \"#806020\" lw 60\n");
	fprintf(ofd,"set label \"SigmaTheta %s\" at %9.3e,%9.3e right font \"Helvetica,8\"\n",st_version,rtmx,rtmy);
	}
    fprintf(ofd,"plot ");
    for(i=0;i<nbf;++i)
	{
	fprintf(ofd,"\"%s\" using 1:",&names[i][1]);
	if (names[i][0]=='-') fprintf(ofd,"(-$2)");
	else fprintf(ofd,"2");
	fprintf(ofd," title \"%s\" with linespoints ls %d",names[i], i+1);
	if (i<nbf-1) fprintf(ofd,",");
	}
    if (ind_gcod)
	{
	if (nbf==4) --nbf;
	fprintf(ofd,", ");
    	for(i=0;i<nbf;++i)
		{
		if (names[i][0]=='-') names[i][0]='+';
		else names[i][0]='-';
		fprintf(ofd,"\"%s\" using 1:",&names[i][1]);
		if (names[i][0]=='-') fprintf(ofd,"(-$2)");
		else fprintf(ofd,"2");
		fprintf(ofd," title \"%s\" with linespoints ls %d",names[i], i+5);
		if (i<nbf-1) fprintf(ofd,",");
		}
	}
    fprintf(ofd,"\n");
    fprintf(ofd,"exit\n");
    fclose(ofd);
    strcpy(sys_cmd,"gnuplot ");
    strcat(sys_cmd,gptfile);
    err=system(sys_cmd);
    return(err);
    }

int gen_3chplt(char input[][256], char *outfile, int N, double tau[], double gcod[][256], double bmax[][256], int flagest)
/* Generate a gnuplot file (.gnu) and invoke gnuplot for creating a pdf file */
    {
    int i,j,k,mii,mxi,err, nbf=4, no_cls=0;
    double minx, maxx, miny, maxy, lmix, lmax, lmiy, lmay, ltx, lty, lmx, lmy, rtmx, rtmy;
    char nomest[256], gptfile[256], psfile[256], gpt_cmd[65536], sys_cmd[256];
    FILE *ofd;

    if (!input[3][0])
	no_cls=1;
    strcpy(gptfile,outfile);
    strcat(gptfile,".gnu");
    strcpy(psfile,outfile);
    strcat(psfile,".pdf");
    ofd=fopen(gptfile, "w");
    if (ofd==NULL) 
	{
	printf("# File %s note created\n",gptfile);
	return(-1);
	}
    if (flag_display)
	fprintf(ofd,"set terminal wxt size 1024,768 enhanced font \"Helvetica\" fontscale 1.5 persist\n");
    else
	{
	fprintf(ofd,"set terminal pdfcairo size 172,128 enhanced color font \"Helvetica\" fontscale 18\n");
	fprintf(ofd,"set output \"%s\"\n",psfile);
	}
    fprintf(ofd,"set logscale xy\n");
    fprintf(ofd,"set format xy \"10^{%%+T}\"\n");
    fprintf(ofd,"set grid\n");
    minx=miny=1e99;
    maxx=maxy=db(0);
    for(i=0;i<N;++i)
	{
	if (tau[i]<minx) minx=tau[i];
	if (tau[i]>maxx) maxx=tau[i];
	}
    for(j=0;j<nbf-1;++j)
	for(i=0;i<N;++i)
		{
		if (fabs(gcod[j][i])<miny)
	    		{
	    		mii=i;
	    		miny=fabs(gcod[j][i]);
	    		}
		if (fabs(bmax[j][i])>maxy)
	    		{
	    		mxi=i;
	    		maxy=fabs(bmax[j][i]);
	    		}
		}
/*    if (mxi<mii)
	fprintf(ofd,"set key right\n");
    else*/
    fprintf(ofd,"set key left bottom Left reverse\n");
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
    fprintf(ofd,"set mxtics 10\n");
    fprintf(ofd,"set mytics 10\n");
    if (flag_title) fprintf(ofd,"set title \"%s\" noenhanced\n",outfile);
    fprintf(ofd,"set xlabel \"Integration time {/Symbol t} [s]\"\n");
    fprintf(ofd,"set ylabel \"ADEV {/Symbol s}_A({/Symbol t})\"\n");
    fprintf(ofd,"set style fill transparent solid 0.25 border\n");
    if (flag_display)
	{
	fprintf(ofd,"set style line 1 lt 1 pt 7 ps 1.5 lc rgb \"#FF0000\" lw 1\n");
	fprintf(ofd,"set style line 2 lt 1 pt 7 ps 1.5 lc rgb \"#00FF00\" lw 1\n");
	fprintf(ofd,"set style line 3 lt 1 pt 7 ps 1.5 lc rgb \"#0000FF\" lw 1\n");
	fprintf(ofd,"set style line 11 lt 1 pt 7 ps 1.5 lc rgb \"#FFA0A0\" lw 1\n");
	fprintf(ofd,"set style line 12 lt 1 pt 7 ps 1.5 lc rgb \"#A0FFA0\" lw 1\n");
	fprintf(ofd,"set style line 13 lt 1 pt 7 ps 1.5 lc rgb \"#A0A0FF\" lw 1\n");
	fprintf(ofd,"set style line 21 lt 1 pt 4 ps 1.5 lc rgb \"#D01000\" lw 2\n");
	fprintf(ofd,"set style line 22 lt 1 pt 12 ps 1.5 lc rgb \"#308015\" lw 2\n");
	fprintf(ofd,"set style line 23 lt 1 pt 6 ps 1.5 lc rgb \"#0010D0\" lw 2\n");
	fprintf(ofd,"set style line 24 lt 1 pt 13 ps 1.5 lc rgb \"#806020\" lw 2\n");
	fprintf(ofd,"set style line 31 lt 2 dt 2 pt 3 ps 1.5 lc rgb \"#D01000\" lw 2\n");
	fprintf(ofd,"set style line 32 lt 2 dt 2 pt 1 ps 1.5 lc rgb \"#308015\" lw 2\n");
	fprintf(ofd,"set style line 33 lt 2 dt 2 pt 2 ps 1.5 lc rgb \"#0010D0\" lw 2\n");
	fprintf(ofd,"set label \"SigmaTheta %s\" at %9.3e,%9.3e right font \"Helvetica,6\"\n",st_version,rtmx,rtmy);
	}
    else
	{
	fprintf(ofd,"set style line 1 lt 1 pt 7 ps 22 lc rgb \"#FF0000\" lw 20\n");
	fprintf(ofd,"set style line 2 lt 1 pt 7 ps 22 lc rgb \"#00FF00\" lw 20\n");
	fprintf(ofd,"set style line 3 lt 1 pt 7 ps 22 lc rgb \"#0000FF\" lw 20\n");
	fprintf(ofd,"set style line 11 lt 1 pt 7 ps 22 lc rgb \"#FFA0A0\" lw 20\n");
	fprintf(ofd,"set style line 12 lt 1 pt 7 ps 22 lc rgb \"#A0FFA0\" lw 20\n");
	fprintf(ofd,"set style line 13 lt 1 pt 7 ps 22 lc rgb \"#A0A0FF\" lw 20\n");
	fprintf(ofd,"set style line 21 lt 1 pt 4 ps 22 lc rgb \"#D01000\" lw 60\n");
	fprintf(ofd,"set style line 22 lt 1 pt 12 ps 32 lc rgb \"#308015\" lw 60\n");
	fprintf(ofd,"set style line 23 lt 1 pt 6 ps 26 lc rgb \"#0010D0\" lw 60\n");
	fprintf(ofd,"set style line 24 lt 1 pt 13 ps 20 lc rgb \"#806020\" lw 60\n");
	fprintf(ofd,"set style line 31 lt 2 dt 2 pt 3 ps 24 lc rgb \"#D01000\" lw 60\n");
	fprintf(ofd,"set style line 32 lt 2 dt 2 pt 1 ps 28 lc rgb \"#308015\" lw 60\n");
	fprintf(ofd,"set style line 33 lt 2 dt 2 pt 2 ps 22 lc rgb \"#0010D0\" lw 60\n");
	fprintf(ofd,"set label \"SigmaTheta %s\" at %9.3e,%9.3e right font \"Helvetica,8\"\n",st_version,rtmx,rtmy);
	}
    fprintf(ofd,"plot ");
// 68% confidence interval
    for(i=0;i<nbf-1;++i)
	{
	fprintf(ofd,"\"%s\" using 1:5:6",input[i]);
	fprintf(ofd," title \"68%% conf. int. %d\" with filledcurves ls %d, ", i+1, i+1);
	}
// 95% confidence interval
    for(i=0;i<nbf-1;++i)
	{
	fprintf(ofd,"\"%s\" using 1:4:7",input[i]);
	fprintf(ofd," title \"95%% c. i. %d\" with filledcurves ls %d, ", i+1, i+11);
	}
// Direct estimates
    for(i=0;i<nbf-1;++i)
	{
	fprintf(ofd,"\"%s\" using 1:2",input[i]);
	fprintf(ofd," title \"direct estimates %d\" with linespoints ls %d, ", i+1, i+21);
	}
// Mean or median estimates
    switch (flagest)
	{
	case 0:
		strcpy(nomest,"mean");
		break;
	case 1:
		strcpy(nomest,"median");
		break;
	default:
		strcpy(nomest,"KLT");
	}
    for(i=0;i<nbf-2;++i)
	{
	fprintf(ofd,"\"%s\" using 1:3",input[i]);
	fprintf(ofd," title \"%s estimates %d\" with linespoints ls %d, ",nomest, i+1, i+31);
	}
    fprintf(ofd,"\"%s\" using 1:3",input[nbf-2]);
    fprintf(ofd," title \"%s estimates %d\" with linespoints ls %d",nomest, nbf-1, nbf+29);
// Closure
    if (no_cls)
	fprintf(ofd,"\n");
    else
	{
	fprintf(ofd," ,\"%s\" using 1:2",input[3]);
	fprintf(ofd," title \"measurement noise\" with linespoints ls %d\n", 24);
	}
    fprintf(ofd,"exit\n");
    fclose(ofd);
    strcpy(sys_cmd,"gnuplot ");
    strcat(sys_cmd,gptfile);
    err=system(sys_cmd);
    return(err);
    }
