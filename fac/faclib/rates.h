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
		     int np, void *params, int i0, int f0, int type,
		     double (*Rate1E)(double, double, int, void *));
double IntegrateRate2(int idist, double e, int np, 
		      void *params, int i0, int f0, int type,
		      double (*Rate2E)(double, double, double, int, void *));

double CERate1E(double e, double eth, int np, void *p);
double DERate1E(double e, double eth, int np, void *p);
int CERate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, double *params, int i0, int f0);

int TRRate(double *dir, double *inv, int iinv,
	   int j1, int j2, double e, float strength);

double CIRate1E(double e, double eth, int np, void *p);
double R3BRate1E(double e1, double e2, double eth, int np, void *p);
int CIRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, float *params, int i0, int f0);
double RRRateHydrogenic(double t, double z, int n, double *top);
double RRRate1E(double e, double eth, int np, void *p);
double PIRate1E(double e, double eth, int np, void *p);
int RRRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, double *params, int i0, int f0);
int AIRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e, float rate);
double RRFit(int z, int nele, double t);
double DRFit(int z, int nele, double t);
double PhFit2(int z, int nele, int is, double e);
double CFit(int z, int nele, double t, double *a, double *dir);
double ColFit(int z, int nele, int is, double t, double *a, double *dir);
double Ionis(int z, int nele, double t, double *a, double *dir, int m);
double Recomb(int z, int nele, double t, double *rr, double *dr, int m);
int FracAbund(int z, double t, double *a, int im, int rm);
double MaxAbund(int z, int nele, double *a, double eps, int im, int rm);
double TwoPhotonRate(double z, int t);
int InitRates(void);	   

#endif
