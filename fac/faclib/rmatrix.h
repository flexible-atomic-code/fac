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
  int **basis;
  double cp0[2], cp1[2], cp2[3];
  double **ebuttle, **cbuttle[NBTERMS];
  double **ek, **w0, **w1;
} RBASIS;

typedef struct _RMATRIX_ {
  int nts, nkappa, nchan, nchan0, mchan;
  int *chans, *ilev, *kappa, *ts, *jts;
  int ndim, nop, nlam;
  int nsym, isym, p, j;
  double **aij;
  double et0, *et, *ek, **w0, **w1;
  double *rmatrix[3];
  double z, energy;
} RMATRIX;

typedef struct _DCFG_ {
  int *iwork, diag;
  double *dwork;
  double *fs, *fc, *gs, *gc;
  double *fs0, *fc0, *gs0, *gc0;
  double *a, *b, *c, *d, *e, *p, *p2, *rm;
  int nmultipoles, ngailitis;
  double rgailitis, degenerate, accuracy;
} DCFG;

int InitRMatrix(void);
void ClearRMatrixBasis(RBASIS *rbs);
void ReadRMtraixBasis(char *fn, RBASIS *rbs);
void WriteRMatrixBasis(char *fn);
void RMatrixBoundary(double r0, double r1, double b);
void ExtrapolateButtle(RBASIS *rbs, int t, int m, double *e,
		       double *r0, double *r1, double *r2);
int RMatrixBasis(char *fn, int kmax, int nb);
int IndexFromKappa(int k);
int KappaFromIndex(int i);
void RMatrixTargets(int nt, int *kt, int nc, int *kc);
void RMatrixNMultipoles(int n);
void ClearRMatrixSurface(RMATRIX *rmx);
int ReadRMatrixSurface(FILE *f, RMATRIX *rmx, int m);
int WriteRMatrixSurface(FILE *f, double **wik0, double **wik1, int m);
int RMatrixSurface(char *fn);
int RMatrix(double e, RMATRIX *rmx, RBASIS *rbs, int m);
int RMatrixPropogate(double *r0, double *r1, RMATRIX *rmx1);
int RMatrixKMatrix(RMATRIX *rmx0, RBASIS *rbs, double *r0);
int SMatrix(RMATRIX *rmx0);
void RMatrixExpansion(int n, double d, double a, double r);
int GailitisExp(RMATRIX *rmx, double r);
int RMatrixCE(char *fn, int np, char *bfn[], char *rfn[], 
	      double emin, double emax, double de, int m, int mb);
void TestRMatrix(double e, int m, char *fn1, char *fn2, char *fn3);

#endif
