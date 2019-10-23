/*   dev_sbr.c                                     F. Vernotte - 2014/10/02 */
/*                               From adev_sbr.c, created by FV, 2010/10/20 */
/*                               From mdev_sbr.c, created by FV, 2014/09/24 */
/*		     Adding of Parabolic deviation (PDev) by FV, 2015/02/06 */
/*       Subroutines for the computation of several Deviations (ADEV, MDev) */
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
/* requirements in conditions enabling the security of their systoms and/or */
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

#ifndef AVAR
// Classical Allan variance
#define AVAR 0 
// Modified Allan variance
#define MVAR 1
// Hadamard variance
#define HVAR 2
// Parabolic variance
#define PVAR 3
// Groslambert covariance
#define GVAR 4
#endif

extern double *T, *Y, *Y1, *Y2, ortau[], log_inc;
extern char flag_log_inc, flag_variance;
extern int ntau;

double adev_y(int tau, int ny)
/* adev_y(tau,ny) : compute the Allan deviation of the 'ny' frequency deviation elements of the vector 'Y' with an integration time 'tau'.*/
    {
    int i,j,nt;
    double Myd,Myf,al;

    nt=ny-2*tau+1;
    Myd=Myf=(double)0;
    for (i=0;i<tau;++i)
	{
	Myd+=Y[i];
	Myf+=Y[i+tau];
	}
    if (nt==1)
      al=sqrt((double)2)*fabs(Myd-Myf)/((double)ny);
    else
        {
	al=pow(Myf-Myd,((double)2));
        for(i=1;i<nt;++i)
	    {
            Myd+=Y[tau+i-1]-Y[i-1];
            Myf+=Y[2*tau+i-1]-Y[tau+i-1];
	    al+=pow(Myf-Myd,((double)2));
	    }
        al=sqrt(al)/(((double)tau)*sqrt((double)(2*nt)));
	}
    return(al);
    }
  
double gcodev_y(int tau, int ny)
/* gcodev_y(tau,ny) : compute the Groslambert codeviation (see arXiv:1904.05849) of the 'ny' frequency deviation elements of the vectors 'Y1' and 'Y2' with an integration time 'tau'.*/
    {
    int i,j,nt;
    double Myd1,Myf1,Myd2,Myf2,al,bl;

    nt=ny-2*tau+1;
    Myd1=Myf1=Myd2=Myf2=(double)0;
    for (i=0;i<tau;++i)
	{
	Myd1+=Y1[i];
	Myf1+=Y1[i+tau];
	Myd2+=Y2[i];
	Myf2+=Y2[i+tau];
	}
    if (nt==1)
	{
        al=(Myd1-Myf1)*(Myd2-Myf2);
	if (al>0)
		bl=sqrt( ((double)2) * al ) / ((double)ny);
	else
		bl=-sqrt( ((double)-2) * al ) / ((double)ny);
	}
    else
        {
	al=(Myf1-Myd1)*(Myf2-Myd2);
        for(i=1;i<nt;++i)
	    {
            Myd1+=Y1[tau+i-1]-Y1[i-1];
            Myf1+=Y1[2*tau+i-1]-Y1[tau+i-1];
            Myd2+=Y2[tau+i-1]-Y2[i-1];
            Myf2+=Y2[2*tau+i-1]-Y2[tau+i-1];
	    al+=(Myf1-Myd1)*(Myf2-Myd2);
	    }
    	if (al>0)
        	bl=sqrt(al)/(((double)tau)*sqrt((double)(2*nt)));
    	else
        	bl=-sqrt(-al)/(((double)tau)*sqrt((double)(2*nt)));
	}
//    if (isnan(bl)) printf("Myd1=%12.6le, Myf1=%12.6le, Myd2=%12.6le, Myf2=%12.6le, al=%12.6le, ny=%d, nt=%d\n",Myd1,Myf1,Myd2,Myf2,al,ny,nt);
    return(bl);
    }

double mdev_y(int tau, int ny) /* Dudullized modified Allan deviation (Francois Meyer, 199?) */
/* mdev_y(tau,ny) : compute the modified Allan deviation of the 'ny' frequency deviation elements of the vector 'Y' with an integration time 'tau'.*/
	{
	int i,j,k,l,m,to2, to3m1, to2m1, tom1;
	double result,moypon,coco, qto_ajto;

	moypon=0;
	to2=2*tau;

	l=3*tau-2;
	for(i=0;i<tau;++i)
		moypon+=(i+1)*(Y[i]-Y[l-i]);

	k=tau-2;
	l=2*tau-1;
	for(i=tau;i<l;++i)
		{
		moypon+=k*Y[i];
		k-=2;
		}

	result=moypon*moypon;

	qto_ajto=0;

	for(j=0;j<tau;++j)
		qto_ajto+=-Y[j+to2]-Y[j];

	for(j=tau;j<to2;++j)
		qto_ajto+=2*Y[j];

	l=ny-tau*3;
	to3m1=3*tau-1;
	to2m1=2*tau-1;
	tom1=tau-1;

	for(i=1;i<l;++i)
		{
		moypon+=qto_ajto;
		result+=moypon*moypon;
		qto_ajto+=Y[i-1]+3*(Y[i+to2m1]-Y[i+tom1])-Y[i+to3m1];
		}

	coco=(double)tau;
	coco*=coco;
	coco*=coco;
	result/=((double)2)*coco*((double)(ny-3*tau+1));
	return(sqrt(result));
	}

double hadamard_y(int tau, int ny)
/* mdev_y(tau,ny) : compute the Hadamard deviation of the 'ny' frequency deviation elements of the vector 'Y' with an integration time 'tau'.*/
	{
	int i,j,ifin,cpi,te;
	double y1,y2,y3,dy,dify,tau2;

	tau2=((double)tau)*((double)tau);
	if ((tau<1)||(tau>(ny/3))) dify=0;
	else
		{
		y1=y2=y3=(double)0;
		for (i=0;i<tau;++i)
			{
			y1+=(double)Y[i];
			y2+=(double)Y[i+tau];
			y3+=(double)Y[i+2*tau];
			}
		ifin=ny-2*tau;
		dify=0;
		cpi=1;
		if (tau<(ny/3))
			{
			for(i=tau;i<ifin;++i)
				{
				++cpi;
				dy=((double)2)*y2-y1-y3;
				dify+=dy*dy/tau2;
				y1-=(double)Y[i-tau];
				y1+=(double)Y[i];
				y2-=(double)Y[i];
				y2+=(double)Y[i+tau];
				y3-=(double)Y[i+tau];
				y3+=(double)Y[i+2*tau];
				}
			dy=((double)2)*y2-y1-y3;
			dify+=dy*dy/tau2;
			dify/=(double)(cpi);
			}
		else
			dify=(((double)2)*y2-y1-y3)*(((double)2)*y2-y1-y3)/tau2;
		dify/=(double)9;
		}
	return(sqrt(dify));
	}
  
double pdev_y(int tau, int ny) /* Parabolic deviation (Vernotte et al. 2015, arXiv:1506.00687) */
/* pdev_y(tau,ny) : compute the Parabolic deviation of the 'ny' frequency deviation elements of the vector 'Y' with an integration time 'tau'.*/
	{
	double *H0,OV,Ovi;
	int dtmx, Tm, Tp, *Tc, bb, N, fin, i, ia, ja;

	dtmx=2*tau-2;
	H0=(double *)malloc(dtmx*sizeof(double)); 
	Tc=(int *)malloc(dtmx*sizeof(int)); 
	if (tau==1)
		{
		Tc[0]=0;
		Tc[1]=1;
		H0[0]=(double)1;
		H0[1]=(double)-1;
		bb=2;
		N=ny-1;
		fin=0;
		}
	else
		{
		for(i=0;i<tau-1;++i)
			{
			Tm=-tau+1+i;
			Tp=i+1;
			H0[i]=((double)(-6*Tm*(Tm+tau)))/pow((double)tau,3);
			H0[i+tau-1]=((double)(6*Tp*(Tp-tau)))/pow((double)tau,3);
			Tc[i]=Tm+tau-1;
			Tc[i+tau-1]=Tp+tau-1;
			bb=2*tau-2;
			N=ny-2*tau+1;
			fin=0;
			}
		}
	OV=0;
	for(ia=0;ia<N;++ia)
		{
		Ovi=0;
		for(ja=0;ja<bb;++ja)
			Ovi+=H0[ja]*Y[Tc[ja]+ia+fin];
		OV+=Ovi*Ovi;
		}
	OV/=(double)(2*N);
	return(sqrt(OV));
	}

int serie_dev(int N, double *tau, double *dev)
/* Compute DEV (ADEV or MDEV) serie from tau=tau0 (tau0=sampling step) to tau=N*tau0/2 or N*tau0/3 (N=number of samples) according to the variance type (flag_variance) by octave (log step : tau_n+1=2*tau_n, default setting) or other if specified in the configuration file '.SigmaTheta.conf' */
/* Input  : number of frequency deviation  measurements */
/* Output : number of Dev measurements */
/*          tau \t Dev(tau) */
    {
    int i,j,toi[256],nto,tomax,ndec,d0,psup,indt,toto;
    double tn0,smpt;

    smpt=T[1]-T[0];
    toi[0]=1;
    if ((flag_variance==1)||(flag_variance==2)) tomax=(int)floor((double)(N/3));
    else tomax=(int)floor((double)(N/2));
    if (flag_log_inc)
	{
    	i=0;
    	do
        	{
		i++;
		toi[i]=(int)pow(log_inc,((double)i))*((double)toi[0]);
        	}
    	while(toi[i]<tomax);
    	if (toi[i-1]<tomax)
        	{
        	toi[i]=tomax;
        	nto=i+1;
		}
    	else
        	nto=i;
	}
    else
	{
	ndec=(int)(log((double)(tomax-toi[0]))/log((double)10));
	d0=(int)floor(log(((double)toi[0])*smpt)/log(10));
	tn0=((double)toi[0])*smpt/pow((double)10,(double)d0);
	psup=0;
	for (i=1;i<=ntau;++i) 
		if (ortau[ntau-i]>tn0) psup=ntau-i;
	if (!psup) ++d0;
	indt=1;
	for (i=psup;i<ntau;++i)
		{
		toi[indt]=(int)(ortau[i]*pow((double)10,(double)d0)/smpt);
		if (toi[indt]==toi[indt-1]) ++toi[indt];
		++indt;
		}
	for (j=d0+1;j<d0+ndec+1;++j)
		{
		for(i=0;i<ntau;++i)
			{
			toi[indt]=(int)(ortau[i]*pow((double)10,(double)j)/smpt);
			++indt;
			if (toi[indt-1]>=tomax)
				{
				toi[indt-1]=tomax;
				break;
				}
			}
		}
	nto=indt;
    }
    for(i=0;i<nto;++i) {
        tau[i]=((double)toi[i])*smpt;
        switch(flag_variance)
        {
            case MVAR :
                dev[i]=mdev_y(toi[i],N);
                break;
            case HVAR :
                dev[i]=hadamard_y(toi[i],N);
                break;
            case PVAR :
                dev[i]=pdev_y(toi[i],N);
                break;
            case GVAR :
                dev[i]=gcodev_y(toi[i],N);
                break;
            case AVAR :
            default :
                dev[i]=adev_y(toi[i],N);
        }
    }
    return(nto);
    }  
