#ifndef _RATES_H_
#define _RATES_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "dbase.h"

#define MAX_DIST  32

#define RT_CE 1
#define RT_CI 2
#define RT_RR 3

typedef struct _DISTRIBUTION_ {
  int xlog;
  int nparams;
  double *params;
  double (*dist)(double, double *);
} DISTRIBUTION;

double VelocityFromE(double e, double bms);
int EleDist(char *fn, int n);
int PhoDist(char *fn, int n);
int SetEleDist(int i, int np, double *p);
int SetPhoDist(int i, int np, double *p);
int DistFromFile(char *fn, double **p);

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
int CIRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, float *params, int i0, int f0);
double RRRateHydrogenic(double t, double z, int n, double *top);
double RRRate1E(double e, double eth, int np, void *p);
double PIRate1E(double e, double eth, int np, void *p);
double PIRateKramers(double e, double eth, int np, void *p);
int RRRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, double *params, int i0, int f0);
int AIRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e, float rate);
double RRFit(int z, int nele, double t);
double DRFit(int z, int nele, double t);
double NRRFit(int z, int nele, double t);
double NDRFit(int z, int nele, double t);
double BremssNR(int z, double te, double e);
double PhFit2(int z, int nele, int is, double e);
double CBeli(int z, int nele, double ene, 
	     double *a, double *dir, double *err);
double EBeli(int z, int nele);
double EColFit(int z, int nele, int is);
double EPhFit2(int z, int nele, int is);
double RBeli(int z, int nele, double t, double *a, double *dir);
double CFit(int z, int nele, double t, double *a, double *dir);
double ColFit(int z, int nele, int is, double t, double *a, double *dir);
double CColFit(int z, int nele, int is, double t, double *a, double *dir);
double Ionis(int z, int nele, double t, double *a, double *dir, int m);
double Recomb(int z, int nele, double t, double *rr, double *dr, int m);
int FracAbund(int z, double t, double *a, int im, int rm);
double MaxAbund(int z, int nele, double *a, double eps, int im, int rm);
double TwoPhotonRate(double z, int t);
int InitRates(void);	   

#endif
