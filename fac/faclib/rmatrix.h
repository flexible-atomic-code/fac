#ifndef _RMATRIX_H_ 
#define _RMATRIX_H_

#include "structure.h"
#include "radial.h"

#define NBTERMS 5
#define NBFIT 6
typedef struct _RBASIS_ {
  int kmax, nbk, nkappa, nbuttle;
  int ib0, ib1;
  double rb0, rb1, bqp;
  int **basis;
  double **ebuttle, **cbuttle[NBTERMS];
  double **ek, **w0, **w1;
} RBASIS;

int InitRMatrix(void);
void ClearRMatrixBasis(void);
void WriteRMatrixBasis(char *fn);
void RMatrixBoundary(double r0, double r1, double b);
int RMatrixBasis(char *fn, int kmax, int nb);
int IndexFromKappa(int k);
int KappaFromIndex(int i);
void RMatrixTargets(int nt, int *kt, int nc, int *kc);
void WriteRMatrixSurface(FILE *f, double **wik0, double **wik1);
int RMatrixSurface(char *fn);

#endif
