#ifndef _COULOMB_H_
#define _COULOMB_H_ 1

#include "global.h"
#include "angular.h"

#define MAXNKL         50
#define MAXKL          1000
#define MAXNUSR        30
#define MAXNE          15
#define MAXNTE         6
#define CBMULTIPOLES   2
#define MAXNCB         (CBMULTIPOLES*(CBMULTIPOLES+3)/2)

int AddPW(int *nkl0, double *kl, double *logkl, 
	  int maxkl, int n, int step);
int SetPWGrid(int *nkl0, double *kl, double *logkl, 
	      int maxkl, int *ns, int *n, int *step);
int SetTEGridDetail(double *te, double *logte, int n, double *x);
int SetTEGrid(double *te, double *logte, int n, double emin, double emax);
int SetEGridDetail(double *e, double *log_e, int n, double *xg);
int SetEGrid(double *e, double *log_e, 
	     int n, double emin, double emax, double eth);
double CoulombPhaseShift(double z, double e, int kappa);
double AngularMSub(int lf, int li1, int li2, int q);
double *GetCoulombBethe(int ie2, int ite, int ie1, int m, int k);
double GetCoulombBetheAsymptotic(double te, double e1);
int PrepCBIndex(int mode);
int CoulombBetheTail(int n, double *w, int nkl, double *kl, double *tcb);
int PrepCoulombBethe(int ne2, int nte, int ne1, double z,
		     double *e2, double *te, double *e1,
		     int nkl, double *kl, 
		     int etype, int ltype, int mode);
int TestCoulomb();
int InitCoulomb();

#endif
