#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "global.h"
#include "dbase.h"

void spline_work(double *x, double *y, int n, 
		 double yp1, double ypn, double *y2, double *work);
void spline(double *x, double *y, int n, 
	    double yp1, double ypn, double *y2);
int splint(double *xa, double *ya, double *y2a, 
	   int n, double x, double *y);
void splie2(double *x1a, double *x2a, double **ya, 
	    int m, int n, double **y2a);
void splin2(double *x1a, double *x2a, double **ya, double **y2a,
	   int m, int n, double x1, double x2, double *y);
void PolyBasis(int n, double *c, double x, double logx);
void PolyFit(int n, double *c, int nd, double *x, double *y);
void SVDFit(int np, double *coeff, double *chisq, double tol,
	    int nd, double *x, double *logx, double *y, double *sig,
	    void Basis(int, double *, double, double));
int NLSQFit(int np, double *p, double tol, int *ipvt,
	    double *fvec, double *fjac, int ldfjac, double *wa, int lwa,
	    int n, double *x, double *logx, double *y, double *sig,
	    void func(int, double *, int , double *, double *, 
		      double *, double *, int, void *), 
	    void *extra);
double Simpson(double *y, int ia, int ib);
int NewtonCotes(double *r, double *x, int i0, int i1, int m, int id);

double RRCrossHn(double z, double e, int n);
int PrintRRTable(FILE *f1, FILE *f2, int v, int swp);
void PrepCECrossHeader(CE_HEADER *h, double *data);
void PrepCECrossRecord(int k, CE_RECORD *r, CE_HEADER *h,
		       double *data);
double InterpolateCECross(double e, CE_RECORD *r, CE_HEADER *h,
			  double *data, double *ratio);
int CECross(char *ifn, char *ofn, int i0, int i1, 
	    int negy, double *egy, int mp);
int CEMaxwell(char *ifn, char *ofn, int i0, int i1, 
	      int nt, double *temp);
void PrepCEFCrossHeader(CEF_HEADER *h, double *data);
void PrepCEFCrossRecord(CEF_RECORD *r, CEF_HEADER *h,
			double *data);
double InterpolateCEFCross(double e, CEF_RECORD *r, CEF_HEADER *h,
			   double *data);
int CEFCross(char *ifn, char *ofn, int i0, int i1, 
	     int negy, double *egy, int mp);
int CEFMaxwell(char *ifn, char *ofn, int i0, int i1, 
	       int nt, double *temp);
int TotalCICross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int imin, int imax);
int TotalPICross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int imin, int imax);
int TotalRRCross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int n0, int n1, 
		 int nmax, int imin, int imax);
double voigt(double a, double v);

#endif
