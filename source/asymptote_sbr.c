/*   asymptote_sbr.c                               F. Vernotte - 1989/02/21 */
/*                            modified for fit selection by FV - 2011/01/19 */
/*   Subroutines for the computation of the asymptotes of a adev            */
/*   measurement serie                                                      */
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

#define db(x) ((double)(x))
#define sisig(x) ( (x) == 0 ) ? (db(0)) : (  ( (x) > 0 ) ? (db(1)) : (db(-1))  )

extern double coeff[];
extern char flag_slopes[];

int resol(np,mat,r)
int np;
double mat[32][32],r[32];
/* Resolution of a linear system with N equations and N unknowns (N<=32)     */
	{
	int k,i,m,erreur,j,kp1,ieie,ip1;
	double rr,tx;

	erreur=0;
	for (i=0;i<=np;++i) coeff[i]=db(0);
	for (k=0;k<np;++k)
		{
		i=k;m=k;
		do
			{
			i+=1;
			if(fabs(mat[i][k])>fabs(mat[m][k])) m=i;
			}
		while(i<np);
		if (mat[m][k]==0)
			{
			erreur=1;
			break;
			}
		if(m!=k)
			{
			for(j=k;j<=np;++j)
				{
				rr=mat[k][j];
				mat[k][j]=mat[m][j];
				mat[m][j]=rr;
				}
			rr=r[k];r[k]=r[m];r[m]=rr;
			}
		kp1=k+1;
		for (i=kp1;i<=np;++i)
			{
			rr=mat[i][k]/mat[k][k];
			mat[i][k]=0.;
			for(j=kp1;j<=np;++j)
				mat[i][j]-=rr*mat[k][j];
			r[i]-=rr*r[k];
			}
		}
	if (erreur==0)
		{
		coeff[np]=r[np]/mat[np][np];
		for(ieie=0;ieie<np;++ieie)
			{		
			i=np-ieie-1;
			tx=0.;
			ip1=i+1;
			for (j=ip1;j<=np;++j)
				tx-=mat[i][j]*coeff[j];
			if (mat[i][i]==0) erreur=1;
			else coeff[i]=(r[i]+tx)/mat[i][i];
			}
		}
	return(erreur);
	}

double interpo(double t, int intyp)
    {
    double rslt;
    switch(intyp)
        {
        case 0 :
	    rslt=db(1)/t/t/t;
	    break;
        case 1 :
	    rslt=db(1)/t/t;
	    break;
        case 2 :
	    rslt=db(1)/t;
	    break;
        case 3 :
	    rslt=db(1);
	    break;
        case 4 :
	    rslt=t;
	    break;
        case 5 :
	    rslt=t*t;
	    break;
        default :
            rslt=0;
        }
    return(rslt);
    }
	

/*	Recherche des coefficients du polynome passant par les		*/
/*	points M(donnee[0][i],donnee[1][i],donnee[2][i]) par		*/
/*	une methode derivee des moindres carres (somme des		*/
/*	incertitudes relatives minimum).				*/

int relatfit(nbm,to,av,war,ord)
int nbm,ord;
double to[32],av[32],war[32];
	{
	double mat[32][32],r[32],tx,rr,kcr,rat,max;
	int i,j,k,m,i1,jj,jjj,nm,error,kp1,ieie,ip1,drex;
	int indneg,indpos,cn[32];

	--nbm;
	--ord;
	for(i=0;i<32;++i)
		cn[i]=0;
	indneg=indpos=0;
	for(i=0;i<=ord;++i)
	  if (!flag_slopes[i])
	      {
	      ++indneg;
	      cn[indneg]=i;
	      }

/*    1) Matrix element computation				*/
	do
		{
		for (j=0;j<=ord;++j)
		    for(m=0;m<=ord;++m)
			{
			mat[j][m]=(double)0;
			for(k=0;k<=nbm;++k)
				{
				kcr=1/war[k];
				mat[j][m]+=kcr*interpo(to[k],j)*interpo(to[k],m)/av[k]/av[k];
				}
			}
		for(m=0;m<=ord;++m)
		    {
		    r[m]=(double)0;
		    for(k=0;k<=nbm;++k)
			    {
			    r[m]+=interpo(to[k],m)/av[k]/war[k];
			    }
		    }
	if (indneg)
		{
		for(i=1;i<=indneg;++i)
			{
			for(j=cn[i]-i+1;j<ord;++j)
				{
				for(k=0;k<=ord;++k) mat[k][j]=mat[k][j+1];
				for(k=0;k<=ord;++k) mat[j][k]=mat[j+1][k];
				r[j]=r[j+1];
				}
			--ord;
			}
		}

/*    2) System resolution 					*/
	error = resol(ord,mat,r);
	if (error==0)
		{
		if (indneg)
			{
			for(i=1;i<=indneg;++i)
				{
				++ord;
				for(j=ord;j>=cn[i];--j)
					coeff[j+1]=coeff[j];
				coeff[cn[i]]=0;
				}
			indneg=0;
			}
		for(i=0;i<=ord;++i)
		    {
		    max=db(0);
		    for(j=0;j<=nbm;++j)
		        {
			rat=(coeff[i]*interpo(to[j],i))/av[j];
		        if (max<rat) max=rat;
			}
		    /*max=1.;*/
		    if ((coeff[i]<=0)||(max<0.15))
				{
				++indneg;
				cn[indneg]=i;
				}
		    }
		}
	indpos=0;
	for(i=0;i<=ord;++i)
	    if(coeff[i]>=0) ++indpos;
	if (indpos==ord+1) indneg=0;
	if (indpos==0)
	    {
	    indneg=0;
	    error=1;
	    }
	}
	while(indneg);
	return(error);		
	}
