#ifndef _EXCITATION_H_
#define _EXCITATION_H_

#include "structure.h"

#define MAX_USR_EGRID 20
#define MAX_EGRID  15
#define MAX_KL     512
#define MAX_TEGRID 5
#define MAX_NKL    50
#define MAX_MSUB   200

typedef struct _EXCIT_TIMING_ {
  clock_t rad_pk;
  clock_t rad_qk;
  clock_t set_kappa;
} EXCIT_TIMING;

typedef struct _CEPW_SCRATCH_ {
  int qr;
  int max_kl;
  double eps_dipole;
  double eps_allowed;
  double eps_forbidden;
  int nkl0;
  int nkl;
  double kl[MAX_NKL+1];
  double log_kl[MAX_NKL];
  double qk[MAX_NKL];
} CEPW_SCRATCH;

CEPW_SCRATCH *GetCEPWScratch();
int FreeExcitationPk(int ie);
int InitExcitation();
int SetTEGrid(int n, double emin, double emax);
int SetTEGridDetail(int n, double *x);
int SetCEPWOptions(int qr, int max, double eps_dipole,
		   double eps_allowed, double eps_forbidden);
int AddCEPW(int n, int step);
int SetCEPWGrid(int ns, int *n, int *step);
int SetEGridDetail(int n, double *x);
int SetEGrid(int n, double emin, double emax, int type);
int SetUsrEGridDetail(int n, double *x, int type);
int SetUsrEGrid(int n, double emin, double emax, int type);

int CERadialPk(int *nkappa, int *nkl, double **pk, 
	       short **kappa0, short **kappa1, int ie,
	       int k0, int k1, int k, int mode);
double CERadialQk(int ie, double te, 
		  int k0, int k1, int k2, int k3, int k);
int CERadialQkMSub(double *rq, int ie, double te, int k0, int k1,
		   int k2, int k3, int k, int kp, int nq, int *q);
double CEMSubAngNR(int lf, int li1, int li2, int k, int q);
int CERadialQkMSubRatio(int k, double *rq, int lf, 
			double z, double e1, double te);
double InterpolatePk(double te, int type, double *pk);

int CollisionStrength(double *s, double *e, int lower, int upper, int msub);
int SaveExcitation(int nlow, int *low, int nup, int *up, int msub, char *fn);
int CEQkTable(char *fn, int k, double te);

#endif
