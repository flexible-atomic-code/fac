#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "global.h"
#include "nucleus.h"
#include "angular.h"
#include "config.h"
#include "structure.h"
#include "excitation.h"
#include "transition.h"
#include "coulomb.h"

#define MAX_COMPLEX 500
typedef struct _REC_COMPLEX_ {
  int n;
  ARRAY *rg;
  int s0;
  int s1;
} REC_COMPLEX;

int InitRecombination();
int FreeRecPk();
int FreeRecAngZ();
int SetAICut(double c);
int SetRRTEGrid(int n, double emin, double emax);
int SetRRTEGridDetail(int n, double *x);
int SetPEGrid(int n, double emin, double emax, double eth);
int SetPEGridDetail(int n, double *x);
int SetPEGridLimits(double min, double max, int type);
int SetUsrPEGridType(int type);
int SetUsrPEGrid(int n, double emin, double emax, double eth);
int SetUsrPEGridDetail(int n, double *x);
int AddRecPW(int n, int step);
int SetRecQkMode(int m, double tol);
int SetRecPWOptions(int kl_interp, int max_kl);
int SetRecSpectator(int n_max, int n_frozen, int n_spec);
int ConstructRecGroupName(char *rgn, char *gn, int n);
int RecStates(int n, int k, int *kg, char *fn);
int RecStatesFrozen(int n, int k, int *kg, char *fn);
int RRRadialQkTable(double *qr, int k0, int k1, int m);
int RRRadialQk(double *rqc, double te, int k0, int k1, int m);
void RRRadialQkFromFit(int np, double *p, int n, double *x, double *logx, 
		       double *y, double *dy, int ndy, void *extra);
void RRRadialQkHydrogenicParams(int np, double *p, int n, int klb);
int BoundFreeOS(double *rqu, int *np, double *p, 
		double *eb, int rec, int f, int m);
int PrepRREGrids(double eth);
int SaveRecRR(int nlow, int *low, int nup, int *up, char *fn, int m);
int SaveAI(int nlow, int *low, int nup, int *up, char *fn, int channel);
int AIRadial1E(double *pk, int kb, int kappaf);
int AIRadialPk(double **pk, int k0, int k1, int kb, int kappaf, int k);
int AIRate(double *rate, double *e, int rec, int f);
int SaveDR(int nf, int *f, int na, int *a, int nb, int *b, int ng, int *g, 
	   char *fna, char *fnt, int channel);
int DROpen(int n, int *nlev, int **ops);
 
#endif
