#ifndef _EXCITATION_H_
#define _EXCITATION_H_

#include "structure.h"

#define MAX_USR_EGRID 20
#define MAX_EGRID  10
#define MAX_KL     512
#define MAX_TEGRID 5
#define MAX_NKL    25
#define MAX_K      12
#define MAX_NKAPPA (MAX_NKL+1)*(MAX_K+1)*4
#define MAX_MSUB   200

typedef struct _EXCIT_TIMING_ {
  clock_t rad_pk;
  clock_t rad_qk;
  clock_t set_kappa;
} EXCIT_TIMING;

int FreeExcitationPk(int ie);
int InitExcitation();
int SetTEGrid(int n, double emin, double emax);
int SetCEPWOptions(int qr, int max, int max_1, int max_f, double tolerence);
int AddCEPW(int n, int step);
int SetCEPWGrid(int ns, int *n, int *step);
int SetEGridDetail(int n, double *x);
int SetEGrid(int n, double emin, double emax, int type);
int SetUsrEGridDetail(int n, double *x, int type);
int SetUsrEGrid(int n, double emin, double emax, int type);
int CEContinuaKappas(int ie, int k, int *nkl, int *nkappa,
		     int *kappa0, int *kappa1);
int CERadialPk(int ie, int k0, int k1, int k, int nkappa,
	       int *kappa0, int *kappa1);
double CERadialQk(int ie, double te, int k0, int k1, int k2, int k3, int k);
int CERadialQkMSub(double *rq, int ie, double te, int k0, int k1,
		   int k2, int k3, int k, int kp, int nq, int *q);
double InterpolatePk(double te, int type, double *pk);

int CollisionStrength(double *s, double *e, int lower, int upper, int msub);
int SaveExcitation(int nlow, int *low, int nup, int *up, int msub, char *fn);
int CEQkTable(char *fn, int k, double te);
void spline(double *x, double *y, int n, 
	    double yp1, double ypn, double *y2);
int splint(double *xa, double *ya, double *y2a, 
	   int n, double x, double *y);
void splie2(double *x1a, double *x2a, double **ya, 
	    int m, int n, double **y2a);
void splin2(double *x1a, double *x2a, double **ya, double **y2a,
	   int m, int n, double x1, double x2, double *y);

#endif
