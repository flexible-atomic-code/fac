/*
 *   FAC - Flexible Atomic Code
 *   Copyright (C) 2001-2015 Ming Feng Gu
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 * 
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "global.h"
#include "array.h"
#include "dbase.h"

typedef struct _MOD_RECORD_ {
  int m;
  int i0;
  int i1;
  int op;
  double c;
  void *r, *h;
} MOD_RECORD;

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
double InterpolateCICross(double e, double eth, CI_RECORD *r, CI_HEADER *h);
double InterpolateCIMCross(double e, double eth, CIM_RECORD *r, CIM_HEADER *h, 
			   int q);
int TotalCICross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int imin, int imax);
int TotalPICross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int imin, int imax);
int TotalRRCross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int n0, int n1, 
		 int nmax, int imin, int imax);
int CICross(char *ifn, char *ofn, int i0, int i1, 
	    int negy, double *egy, int mp);
int CIMCross(char *ifn, char *ofn, int i0, int i1, 
	     int negy, double *egy, int mp);
int CIMaxwell(char *ifn, char *ofn, int i0, int i1, 
	      int nt, double *temp);
int RRCross(char *ifn, char *ofn, int i0, int i1, 
	    int negy, double *egy, int mp);
int RRMaxwell(char *ifn, char *ofn, int i0, int i1, 
	      int nt, double *temp);
int InterpCross(char *ifn, char *ofn, int i0, int i1, 
		int negy, double *egy, int mp);
int MaxwellRate(char *ifn, char *ofn, int i0, int i1, 
		int nt, double *temp);
double voigt(double a, double v);
void ModifyTable(char *fh, char *fn0, char *fn1, char *fnm);

#endif
