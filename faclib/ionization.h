#ifndef _IONIZATION_H_
#define _IONIZATION_H_

#include "structure.h"
#include "excitation.h"
#include "recombination.h"

int FreeIonizationQk(void);
int InitIonization(void);
int ReinitIonization(int m);
int SetIEGrid(int n, double emin, double emax);
int SetIEGridDetail(int n, double *x);
void SetCILQR(int m);
void SetCILMax(int m);
void SetCILMaxEject(int m);
void SetCILCB(int m);
void SetCITol(double t);
int SetCIPWOptions(int qr, int max, int max_eject, int kl_cb, double tol);
int SetCIPWGrid(int ns, int *n, int *step);
int SetCIFormat(int m);
int SetCIEGrid(int n, double emin, double emax, double eth);
int SetCIEGridDetail(int n, double *x);
int SetCIFormat(int m);
int SetCIEGridLimits(double min, double max, int type);
int SetUsrCIEGridType(int type);
int SetUsrCIEGrid(int n, double emin, double emax, double eth);
int SetUsrCIEGridDetail(int n, double *x);
int SetCIQkMode(int m, double tol);
int SetCIMaxK(int k);
int CIRadialQk(double *qk, int ie1, int ie2, int kb, int kbp, int k);
int CIRadialQkIntegrated(double *qku, double te, int kb, int kbp);
void CIRadialQkBasis(int npar, double *yb, double x, double logx);
void CIRadialQkFromFit(int np, double *p, int n, 
		       double *x, double *logx, double *y);
int CIRadialQkBED(double *dp, double *bethe, double *b0, int kl,
		  double *logxe, double *q, double *p, double te);
double *CIRadialQkIntegratedTable(int kb, int kbp);
double *CIRadialQkTable(int kb, int kbp);
int IonizeStrength(double *qku, double *p, double *e, int b, int f);
int SaveIonization(int nb, int *b, int nf, int *f, char *fn);

#endif
