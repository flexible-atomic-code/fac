#ifndef _RMATRIX_H_ 
#define _RMATRIX_H_

#include "structure.h"
#include "radial.h"

typedef struct _RBASIS_ {
  int kmax, nbk, nkappa, nbuttle;
  int **basis;
  double **ebuttle, **cbuttle;
} RBASIS;

int RMatrixBasis(int kmax, int nb);
int IndexKappa(int k);
int RMatrix(char *fn, int nt, int *kt, int nc, int *kc);

#endif
