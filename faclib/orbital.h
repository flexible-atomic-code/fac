#ifndef _ORBITAL_H_
#define _ORBITAL_H_

#include "global.h"
#include "nucleus.h"

#define MAX_POINTS 650
#define CUTOFF_POINTS 10

typedef struct _POTENTIAL_ {
  int flag;
  double max_points;
  double Z[MAX_POINTS]; /*effective atomic number*/
  double N; /*number of electrons*/
  double lambda, lambdap;
  double a, ap;
  double rad[MAX_POINTS], ar, br;
  double dr_drho[MAX_POINTS];
  double drho_dr2[MAX_POINTS];
  double Vc[MAX_POINTS];
  double dVc[MAX_POINTS];
  double dVc2[MAX_POINTS];
  double Vtail[MAX_POINTS];
  double dVtail[MAX_POINTS];
  double dVtail2[MAX_POINTS];
  double U[MAX_POINTS];
  double dU[MAX_POINTS];
  double dU2[MAX_POINTS];
  double W[MAX_POINTS];
  double dW[MAX_POINTS];
  double dW2[MAX_POINTS];
} POTENTIAL;

typedef struct _ORBITAL_ {
  int n;
  int kappa;
  double energy;
  double qr_norm;
  double phase;
  double *wfun;
  int ilast;
} ORBITAL;

#define Large(orb) ((orb)->wfun)
#define Small(orb) ((orb)->wfun + MAX_POINTS)

double *GetVEffective();
int RadialSolver(ORBITAL *orb,  POTENTIAL *pot, double tol);
int RadialBound(ORBITAL *orb, POTENTIAL *pot, double tol);
int RadialFree(ORBITAL *orb, POTENTIAL *pot, double tol);
double InnerProduct(int n, double *p1, double *p2, POTENTIAL *pot);
double Simpson(double *y, int ia, int ib);
int NewtonCotes(double *r, double *x, int i0, int i1, int m);
int SetOrbitalRGrid(POTENTIAL *pot, double rmin, double rmax);
double GetRFromRho(double rho, double a, double b, double r0);
int SetPotentialZ(POTENTIAL *pot, double c);
int SetPotentialVc(POTENTIAL *pot);
int SetPotentialVTail(POTENTIAL *pot);
int SetPotentialU(POTENTIAL *pot, int n, double *u);
int SetPotentialW (POTENTIAL *pot, double e, int kappa);

#endif









