#ifndef _COULOMB_H_
#define _COULOMB_H_ 1

#include "global.h"
#include "angular.h"
#include "grid.h"

#define CBMULTIPOLES   2
#define MAXNCB         (CBMULTIPOLES*(CBMULTIPOLES+3)/2)

void SetHydrogenicNL(int n, int kl);
void GetHydrogenicNL(int *n, int *kl);
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
int TestCoulomb(char *s);
int InitCoulomb(void);

#endif
