#ifndef _GRID_H_
#define _GRID_H_ 1

#include "global.h"

#define MAXNKL         50
#define MAXKL          512
#define MAXNUSR        30
#define MAXNE          20
#define MAXNTE         6

#define TE_MIN_MAX     0.2

int AddPW(int *nkl0, double *kl, double *logkl, 
	  int maxkl, int n, int step);
int SetPWGrid(int *nkl0, double *kl, double *logkl, 
	      int maxkl, int *ns, int *n, int *step);
int SetTEGridDetail(double *te, double *logte, int n, double *x);
int SetTEGrid(double *te, double *logte, int n, double emin, double emax);
int SetEGridDetail(double *e, double *log_e, int n, double *xg);
int SetEGrid(double *e, double *log_e, 
	     int n, double emin, double emax, double eth);

#endif
