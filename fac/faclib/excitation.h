#ifndef _EXCITATION_H_
#define _EXCITATION_H_

#include "coulomb.h"
#include "structure.h"

typedef struct _CEPW_SCRATCH_ {
  int qr;
  int max_kl;
  int kl_cb;
  double tolerence;
  int nkl0;
  int nkl;
  int ns;
  double kl[MAXNKL+1];
  double log_kl[MAXNKL];
  double qk[MAXNKL];
  double y2[MAXNKL];
} CEPW_SCRATCH;

#ifdef PERFORM_STATISTICS
typedef struct _EXCIT_TIMING_ {
  clock_t rad_pk;
  clock_t rad_qk;
  clock_t set_kappa;
} EXCIT_TIMING;
#endif

CEPW_SCRATCH *GetCEPWScratch();
int FreeExcitationPk(int ie);
int InitExcitation();
int SetCETEGrid(int n, double emin, double emax);
int SetCETEGridDetail(int n, double *x);
int SetCEPWOptions(int qr, int max, int kl_cb, double tol);
int AddCEPW(int n, int step);
int SetCEPWGrid(int ns, int *n, int *step);
int SetCEEGridType(int utype, int etype, int ltype);
int SetCEEGridDetail(int n, double *x);
int SetCEEGrid(int n, double emin, double emax, double eth);
int SetUsrCEEGridDetail(int n, double *x);
int SetUsrCEEGrid(int n, double emin, double emax, double eth);

int CERadialPk(int *nkappa, int *nkl, double **pk, 
	       short **kappa0, short **kappa1, int ie,
	       int k0, int k1, int k);
int CERadialQk(double *r, int ie, double te, 
		  int k0, int k1, int k2, int k3, int k);
int CERadialQkMSub(double *rq, int ie, double te, int k0, int k1,
		   int k2, int k3, int k, int kp, int nq, int *q);
double InterpolatePk(double te, int type, double *pk);
int CollisionStrength(double *s, double *e, int lower, int upper, int msub);
int SaveExcitation(int nlow, int *low, int nup, int *up, int msub, char *fn);

#endif
