/*   sigma_theta.h                                F. Vernotte - 2010/10/25  */
/*   Declaration of global variables and functions for sigma_theta.c        */
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

// #define ST_VERSION "4.0"
//#define ST_DATE "2019/06/22"

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

#define GPLDOPLOT (0x01) // whether we actually run the gnuplot file or just build it
#define GPLPNG (0x02)    // whether we include png as an output format to the gnuplot file 
#define GPLPDF (0x04)    // whether we include pdf as an output format to the gnuplot file 
#define GPLX11 (0x08)    // whether we include X11 as an output (opens a graphic window whan running the gnuplot file

#define DATAMAX 16384
#define GRANMAX 67108864

#define MAXCHAR 512
#include <stdint.h>
struct conf_int
    {
    double inf_bound;
    double mean;
    double sup_bound;
    };

double coeff[10], ortau[10], log_inc;
double *T, *Y, *Y1, *Y2, *Y12, *Y23, *Y31, *U;
char st_version[]=ST_VERSION;
char st_date[]=ST_DATE;
char flag_graph, flag_conf, flag_bias, flag_title, flag_fit, flag_asymptote, flag_slopes[6], flag_variance, flag_log_inc, flag_display;
int ntau;

int init_flag();
double scale(char);
double tchebyfit(long);
double adev_y(int, int);
double gcodev_y(int, int);
double mdev_y(int, int);
double hadamard_y(int, int);
double pdev_y(int, int);
int serie_dev(int, double *, double *);
int resol(int, double *, double*);
double interpo(double, int);
int relatfit(int, double *, double *, double *, int);
double sw(double, int);
double sx(double, int, int);
double sz(double, int, int);
double BasicSum(int, int, int, int, int);
void avardof(int , double *, int *, double *);
void pvardof(int , double *, int *, double *);
double cdf_rayleigh(double , double);
double dcdfr(double x, double nu);
struct conf_int raylconfint(double);
struct conf_int raylconfint1s(double);
int load_ykt(char *, char, char, uint8_t);
int load_2yk(char *, char *);
int load_3yk(char *, char *, char *);
int load_adev(char *, double *, double *);
int load_coef(char *);
int load_1col(char *);
int load_3col(char *, double *, double *, double *);
int load_7col(char *, double *, double *, double *, double *, double *, double *, double *);
int gener_gplt(char *, int, double *, double *, double *, char *, uint8_t);
int gen_psdplt(char *, int, double *, double *);
int gen_linplt(char *, int, double *, double *, int);
int gen_gcodplt(char *, char (*)[], int, int, double *, double (*)[], int);
int gen_3chplt(char (*)[], char *, int, double *, double (*)[], double (*)[], int);

