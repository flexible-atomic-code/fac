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
void BEDFromFit(int np, double *p, int n, double *x, double *logx,
		double *y, double *dy, int ndy, void *extra);
int CIRadialQkBED(double *dp, double *b, int kl, double *p);
int SaveCIRadialQkIntegrated(int n, char *s);
int LoadCIRadialQkIntegrated(int n, char *s);
int PrepCIRadialQkIntegrated(int nz, double *z, int na, double *a,
			     int np, int *n, int nte, 
			     double emin, double emax, char *s);
double *CIRadialQkIntegratedTable(int kb, int kbp);
double IntegrateQk(double *qk);
double *CIRadialQkTable(int kb, int kbp);
int IonizeStrength(double *qku, int *n, double **p, double *e, int b, int f);
int SaveIonization(int nb, int *b, int nf, int *f, char *fn);

#endif
