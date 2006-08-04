#ifndef _ORBITAL_H_
#define _ORBITAL_H_

#include "global.h"
#include "nucleus.h"
#include "coulomb.h"
#include "interpolation.h"

typedef struct _POTENTIAL_ {
  int flag;
  int r_core;
  int maxrp;
  double ratio, asymp, rmin;
  double Z[MAXRP]; /*effective atomic number*/
  double N; /*number of electrons*/
  double lambda, a; /* parameter for the Vc */
  double ar, br; /* parameter for the transformation */
  int ib, nb, ib1; 
  double bqp; /* boundary condition */
  double rad[MAXRP];
  double dr_drho[MAXRP];
  double dr_drho2[MAXRP];
  double Vc[MAXRP];
  double dVc[MAXRP];
  double dVc2[MAXRP];
  double U[MAXRP];
  double dU[MAXRP];
  double dU2[MAXRP];
  double W[MAXRP];
  double dW[MAXRP];
  double dW2[MAXRP];
  double uehling[MAXRP];
} POTENTIAL;

typedef struct _ORBITAL_ {
  int n;
  int kappa;
  double energy;
  double qr_norm;
  double *phase;
  double *wfun;
  int ilast;
} ORBITAL;

int GetNMax(void);
double *GetVEffective(void);
double RadialDiracCoulomb(int npts, double *p, double *q, double *r,
			  double z, int n, int kappa);
int RadialSolver(ORBITAL *orb,  POTENTIAL *pot);
int RadialBasis(ORBITAL *orb, POTENTIAL *pot);
int RadialBasis(ORBITAL *orb, POTENTIAL *pot);
int RadialRydberg(ORBITAL *orb, POTENTIAL *pot);
int RadialBound(ORBITAL *orb, POTENTIAL *pot);
int RadialFreeInner(ORBITAL *orb, POTENTIAL *pot);
int RadialFree(ORBITAL *orb, POTENTIAL *pot);
double InnerProduct(int i1, int n, 
		    double *p1, double *p2, POTENTIAL *pot);
void Differential(double *p, double *dp, int i1, int i2);
int SetOrbitalRGrid(POTENTIAL *pot);
double GetRFromRho(double rho, double a, double b, double r0);
int SetPotentialZ(POTENTIAL *pot, double c);
int SetPotentialUehling(POTENTIAL *pot, int vp);
int SetPotentialVc(POTENTIAL *pot);
int SetPotentialU(POTENTIAL *pot, int n, double *u);
int SetPotentialW (POTENTIAL *pot, double e, int kappa);

#endif


