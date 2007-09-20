#ifndef _MBPT_H_
#define _MBPT_H_

#include "dbase.h"
#include "structure.h"
#include "transition.h"

typedef struct _TR_OPT_{
  int mktr;
  int naw;
  char tfn[1024];
  int nlow, nup;
  int *low, *up;
} TR_OPT;

typedef struct _MBPT_TR_ {
  int m, nsym1;
  int isym0;
  int *isym1;
  SYMMETRY *sym0;
  SYMMETRY **sym1;
  double **tma;
  double **rma;
} MBPT_TR;

typedef struct _CORR_CONFIG_ {
  CONFIG *c;
  int inp, inq, np, nq;
  int ig, ncs, kp, kq;
} CORR_CONFIG;

typedef struct _MBPT_HAM_ {
  int n, n2, n3;
  int *ng, *ng2;
  int isym, dim;
  int ibra, iket;
  double a, b, c;
  double *hab1, *hba1;
  double *hab, *hba;
  MBPT_TR *mtr;
} MBPT_HAM;

typedef struct _MBPT_EFF_ {
  int n, n2;
  int nbasis, *basis;
  double *h0, *e0, *heff;
  /* effective hamilton elements for 1-virtual */
  double **hab1, **hba1;
  /* effective hamilton elements for 2-virtual */
  double **hab, **hba;
} MBPT_EFF;

void InitMBPT(void);
int StructureMBPT0(char *fn, double de, double ccut, int n, int *s0, int kmax, 
		   int n1, int *nm, int n2, int *nmp, int n3, int *n3g,
		   int n4, int *n4g, char *gn);
int StructureMBPT1(char *fn, char *fn1, int n, int *s0, int nk, int *nkm, 
		   int n1, int *nm, int n2, int *nmp, int n0);
int StructureReadMBPT(char *fn, char *fn2, int nf, char *fn1[], 
		      int nkg, int *kg, int nkg0);
void SetExtraMBPT(int m);
void SetOptMBPT(int i3rd, int n3, double c);
void SetSymMBPT(int nlev, int *ilev);
void TransitionMBPT(int mk, int naw);
void TRTableMBPT(char *fn, int nlow, int *low, int nup, int *up);

#endif
