/*   3ch_sbr.h                      F. Vernotte - First release: 2019/07/28 */
/* 			    						    */
/*   Subroutines for 3cornered_hat (KLTG & KLTS)			    */
/*									    */
/*                                                   - SIGMA-THETA Project  */
/*                                                                          */
/* Copyright or ? or Copr. Universit? de Franche-Comt?, Besan?on, France    */
/* Contributor: Fran?ois Vernotte, UTINAM/OSU THETA (2012/07/17)            */
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

long *matlab_sort(double *, long, double *, long *);
double pgaussdirectfin(double *, double *, int, double);
double pgaussdirectavecestimbr(double *, double, double, double, double, int, double);
double *eig(double *, int, double *, double *);
double xinvy(double *, double *, double, long);
int ci_kltgs(double *, double, double, int, double (*)[3], double (*)[5]);

