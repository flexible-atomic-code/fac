#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "global.h"

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
int NewtonCotes(double *r, double *x, int i0, int i1, int m, int maxrp);

#endif
