#ifndef _RATES_H_
#define _RATES_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"

#define MAX_DIST  8

#define VelocityFromE(e) (sqrt((e))*(5.93096895E-3))

typedef struct _DISTRIBUTION_ {
  int nparams;
  double *params;
  double (*dist)(double, double *);
} DISTRIBUTION;


int SetEleDist(int i, int np, double *p);
int SetPhoDist(int i, int np, double *p);
DISTRIBUTION *GetEleDist(int *i);
DISTRIBUTION *GetPhoDist(int *i);
int SetRateAccuracy(double epsrel, double epsabs);
double IntegrateRate(int idist, double eth, double bound, 
		     int np, int nshells, void *params,
		     int i0, int f0, int type,
		     double (*Rate1E)(double, double, int, int, void *));
double IntegrateRate2(int idist, double e, int np, 
		      int nshells, void *params,
		      int i0, int f0, int type,
		      double (*Rate2E)(double, double, 
				       double, int, int, void *));

double CERate1E(double e, double eth, int np, int ns, void *p);
double DERate1E(double e, double eth, int np, int ns, void *p);
int CERate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, double *params, int i0, int f0);

int TRRate(double *dir, double *inv, int iinv,
	   int j1, int j2, double e, float strength);

double CIRate1E(double e, double eth, int np, int ns, void *p);
double R3BRate1E(double e1, double e2, double eth, int np, int ns, void *p);
int CIRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, int nshells, float *params, int i0, int f0);

double RRRate1E(double e, double eth, int np, int ns, void *p);
double PIRate1E(double e, double eth, int np, int ns, void *p);
int RRRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, int nshells, float *params, int i0, int f0);
int AIRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e, float rate);

int InitRates(void);	   

#endif
