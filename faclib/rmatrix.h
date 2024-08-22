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
  int ki, kmin, kmax, kapmin, nbi, nbk, nkappa, nbuttle;
  int ib0, ib1;
  double rb0, rb1, bqp, emin;
  int *kappa, **basis, **bnode;
  double cp0[2], cp1[2], cp2[3];
  double **ebuttle, **cbuttle[NBTERMS];
  double **ek, **w0, **w1;
} RBASIS;

typedef struct _RMATRIX_ {
  int nts, nkappa, nchan, nchan0, mchan, ncs;
  int *chans, *ilev, *kappa, *ts, *pts, *jts, *cs, *pcs, *jcs;
  int ndim, nop, nlam;
  int nsym, isym, p, j;
  double **aij;
  double et0, *et, *ec, *ek, **w0, **w1;
  double **rmatrix[3];
  double z;
} RMATRIX;

typedef struct _SMATRIX_ {
  int pp, jj, isym, *jmin, *jmax, *nj, *nk, *sp;
  double ***rp, ***ip;
} SMATRIX;

typedef struct _RMXCE_ {
  int nke, nsp, *isp, nes;
  double *e, *es, **s, **ecs, ****ap, ***sp, ***sp0, de;
  int *rcs;
  SMATRIX *smx;
  IDXARY its;
} RMXCE;
  
typedef struct _DCFG_ {
  int *iwork, diag, nop;
  double *dwork, *rwork;
  double *fs, *fc, *gs, *gc;
  double *fs0, *fc0, *gs0, *gc0;
  double *a, *b, *c, *d, *e, *p, *p2, *rm;
  double *afs, *afc, *ags, *agc;
  double *va, *vb;
  int *xch, *ich;
  int nmultipoles, ngailitis, nlam, pdirection;
  double rgailitis, degenerate, accuracy;
  int nch, nov, lrw, liw;
  RMATRIX *rmx;
  int nr, mr, ierr;  
  double energy;
  int nts, nka, nke, ntk, ike, n0;
  double *ek, *eo, eo0, eo1;
} DCFG;

int InitRMatrix(void);
void ClearRMatrixBasis(RBASIS *rbs);
void ReadRMtraixBasis(char *fn, RBASIS *rbs, int fmt);
void WriteRMatrixBasis(char *fn, int fmt);
void RMatrixBoundary(double r0, double r1, double b);
void ExtrapolateButtle(RBASIS *rbs, int t, int m, double *e,
		       double *r0, double *r1, double *r2);
int RMatrixBasis(char *fn, int kmax, int nb);
int IndexFromKappa(int k, int k0);
int KappaFromIndex(int i, int k0);
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
void IntegrateDiracCoulomb(double z, int ka, double r, double rt,
			   double e, double *t1, double *c1);
void PrepDiracCoulomb(RMATRIX *rmx, RBASIS *rbs, double r);
int GailitisExp(RMATRIX *rmx, RBASIS *rbs, double r);
int IntegrateExternal(RMATRIX *rmx, double r1, double r0);
void TransformQ(RMATRIX *rmx, double b, double r, int m);
void PropogateDirection(int m);
int PropogateExternal(RMATRIX *rmx, RBASIS *rbs);
void RMatrixFMode(int m);
void RMatrixRefine(int n, int m, double r);
int RefineRMatrixEGrid(double **er, int iter,
		       RMXCE *rs, RBASIS *rbs, RMATRIX *rmx);
void SaveRMatrixCE(RMXCE *rs, RBASIS *rbs, RMATRIX *rmx,
		   int iter, char *fn, double wt0);
int RMatrixCE(char *fn, int np, char *bfn[], char *rfn[],	      
	      double emin, double emax, int nst, double *sde,
	      int m, int mb);
int RMatrixCEW(int np, RBASIS *rbs, RMATRIX *rmx,
	       FILE **f, FILE **f1, char *fn, RMXCE *rs, 
	       int nke, double *e, int m, int mb, int idep);
int RMatrixConvert(char *ifn, char *ofn, int m);
void TestRMatrix(double e, int m, char *fn1, char *fn2, char *fn3);
void SetOptionRMatrix(char *s, char *sp, int ip, double dp);
double InterpLinear(double de, int *isp, int nke, double *ea,
		    double *ra, double e);
void SortGroupEnergy(RMXCE *rs, RBASIS *rbs, RMATRIX *rmx);
int GetChanL(int j, int j0, int k0);
int IgnoreChan(RMATRIX *rmx, RBASIS *rbs, int i);

#endif
