#ifndef _IONIZATION_H_
#define _IONIZATION_H_

#include "structure.h"
#include "excitation.h"

#define MAX_USR_CIEGRID 50
#define MAX_CIEGRID 15
#define MAX_CIKL 500
#define MAX_CINKL 30
#define MAX_IEGRID 5
#define N_INTEGRATE 32

int FreeIonizationQk();
int InitIonization();
int SetIEGrid(int n, double emin, double emax);
int SetCIPWOptions(int qr, int max, int max_eject, int kl_cb, double tol);
int AddCIPW(int n, int step);
int SetCIPWGrid(int ns, int *n, int *step);
int SetCIEGrid(int n, double emin, double emax);
int SetCIEGridDetail(int n, double *x);
int SetUsrCIEGrid(int n, double emin, double emax, int type);
int SetUsrEGrid(int n, double emin, double emax, int type);
int CIRadialQk(int ie1, int ie2, int kb, int kbp, int k);
int IonizeStrength(double *s, double *e, int b, int f);
int SaveIonization(int nb, int *b, int nf, int *f, char *fn);

#endif
