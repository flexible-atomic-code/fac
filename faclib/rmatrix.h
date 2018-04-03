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

#ifndef _RMATRIX_H_ 
#define _RMATRIX_H_

#include <float.h>
#include "global.h"
#include "structure.h"
#include "radial.h"

#define NBTERMS 5
#define NBFIT 7
typedef struct _RBASIS_ {
  int kmax, nbk, nkappa, nbuttle;
  int ib0, ib1;
  double rb0, rb1, bqp;
  int **basis, **bnode;
  double cp0[2], cp1[2], cp2[3];
  double **ebuttle, **cbuttle[NBTERMS];
  double **ek, **w0, **w1;
} RBASIS;

typedef struct _RMATRIX_ {
  int nts, nkappa, nchan, nchan0, mchan, ncs;
  int *chans, *ilev, *kappa, *ts, *jts, *cs, *jcs;
  int ndim, nop, nlam;
  int nsym, isym, p, j;
  double **aij;
  double et0, *et, *ec, *ek, **w0, **w1;
  double *rmatrix[3];
  double z, energy;
} RMATRIX;

typedef struct _DCFG_ {
  int *iwork, diag, nop;
  double *dwork, *rwork;
  double *fs, *fc, *gs, *gc;
  double *fs0, *fc0, *gs0, *gc0;
  double *a, *b, *c, *d, *e, *p, *p2, *rm;
  int nmultipoles, ngailitis, nlam, pdirection;
  double rgailitis, degenerate, accuracy;
  int lrw, liw;
  RMATRIX *rmx;
} DCFG;

int InitRMatrix(void);
void ClearRMatrixBasis(RBASIS *rbs);
void ReadRMtraixBasis(char *fn, RBASIS *rbs, int fmt);
void WriteRMatrixBasis(char *fn, int fmt);
void RMatrixBoundary(double r0, double r1, double b);
void ExtrapolateButtle(RBASIS *rbs, int t, int m, double *e,
		       double *r0, double *r1, double *r2);
int RMatrixBasis(char *fn, int kmax, int nb);
int IndexFromKappa(int k);
int KappaFromIndex(int i);
void RMatrixTargets(int nt, int *kt, int nc, int *kc);
void RMatrixNMultipoles(int n);
void ClearRMatrixSurface(RMATRIX *rmx);
int ReadRMatrixSurface(FILE *f, RMATRIX *rmx, int m, int fmt);
int WriteRMatrixSurface(FILE *f, double **wik0, double **wik1, int m, 
			int fmt, RMATRIX *rmx, HAMILTON *h);
int RMatrixSurface(char *fn);
int RMatrix(double e, RMATRIX *rmx, RBASIS *rbs, int m);
int RMatrixPropogate(double *r0, double *r1, RMATRIX *rmx1);
int RMatrixKMatrix(RMATRIX *rmx0, RBASIS *rbs, double *r0);
int SMatrix(RMATRIX *rmx0);
void RMatrixExpansion(int n, double d, double a, double r);
int GailitisExp(RMATRIX *rmx, double r);
int IntegrateExternal(RMATRIX *rmx, double r1, double r0);
void TransformQ(RMATRIX *rmx, double b, double r, int m);
void PropogateDirection(int m);
int PropogateExternal(RMATRIX *rmx, RBASIS *rbs);
void RMatrixNBatch(int n);
void RMatrixFMode(int m);
int RMatrixCE(char *fn, int np, char *bfn[], char *rfn[], 
	      double emin, double emax, double de, int m, int mb);
int RMatrixConvert(char *ifn, char *ofn, int m);
void TestRMatrix(double e, int m, char *fn1, char *fn2, char *fn3);

#endif
