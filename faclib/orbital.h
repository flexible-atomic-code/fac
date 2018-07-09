/*
 *   FAC - Flexible Atomic Code
 *   Copyright (C) 2001-2015 Ming Feng Gu
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 * 
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _ORBITAL_H_
#define _ORBITAL_H_

#include "global.h"
#include "nucleus.h"
#include "config.h"
#include "interpolation.h"

typedef struct _POTENTIAL_ {
  int nfrozen;
  int npseudo;
  int mpseudo;
  int dpseudo;
  int mode;
  int flag;
  int r_core;
  int maxrp;
  int nmax;
  int nse, mse, pse, mvp, pvp, hpvs, hlike;
  double hx0, hx1, chx;
  double hxs, ahx, ihx, rhx, dhx, ratio, asymp, rmin;
  double Z[MAXRP]; /*effective atomic number*/
  double dZ[MAXRP], dZ2[MAXRP];
  double N; /*number of electrons*/
  double lambda; /* parameter for the Vc */
  double a; /* previously used to parameterize the Vc, but now always 0.0*/
  double ar, br, qr; /* parameter for the transformation */
  int ib, nb, ib1;
  double bqp, rb; /* boundary condition */
  double rad[MAXRP];
  double mqrho[MAXRP];
  double dr_drho[MAXRP];
  double dr_drho2[MAXRP];
  double vtr[MAXRP];
  double Vc[MAXRP];
  double dVc[MAXRP];
  double dVc2[MAXRP];
  double qdist[MAXRP];
  double U[MAXRP];
  double dU[MAXRP];
  double dU2[MAXRP];
  double W[MAXRP];
  double dW[MAXRP];
  double dW2[MAXRP];
  double ZVP[MAXRP];
  double dZVP[MAXRP];
  double dZVP2[MAXRP];
  double rfn[NKSEP];
  int nfn[NKSEP];
  double ZSE[NKSEP][MAXRP];
  double dZSE[NKSEP][MAXRP];
  double dZSE2[NKSEP][MAXRP];
  double VT[NKSEP+1][MAXRP];
  double dVT[NKSEP+1][MAXRP];
  double dVT2[NKSEP+1][MAXRP];
  NUCLEUS *atom;
} POTENTIAL;

typedef struct _ORBITAL_ {
  int n;
  int kappa, kv;
  double energy, se, ose, qed;
  double qr_norm;
  double *phase;
  double *wfun;
  double bqp0, bqp1, pdx;
  int ilast, idx, isol;
  double rfn;
  struct _ORBITAL_ *horb;
  struct _ORBITAL_ *rorb;
} ORBITAL;

void InitOrbitalData(void *p, int n);
double *GetVEffective(void);
double RadialDiracCoulomb(int npts, double *p, double *q, double *r,
			  double z, int n, int kappa);
int RadialSolver(ORBITAL *orb,  POTENTIAL *pot);
int RadialBasis(ORBITAL *orb, POTENTIAL *pot);
int RadialBasisOuter(ORBITAL *orb, POTENTIAL *pot);
int RadialRydberg(ORBITAL *orb, POTENTIAL *pot);
int RadialBound(ORBITAL *orb, POTENTIAL *pot);
int RadialFreeInner(ORBITAL *orb, POTENTIAL *pot);
int RadialFree(ORBITAL *orb, POTENTIAL *pot);
double InnerProduct(int i1, int n, 
		    double *p1, double *p2, POTENTIAL *pot);
void Differential(double *p, double *dp, int i1, int i2, double *drdrho);
void DrLargeSmall(ORBITAL *orb, POTENTIAL *pot, double *pr, double *qr);
int SetOrbitalRGrid(POTENTIAL *pot);
double GetRFromRho(double rho, double a, double b, double q, double r0);
int SetPotentialExtraZ(POTENTIAL *pot, int iep);
int SetPotentialZ(POTENTIAL *pot);
int SetPotentialVP(POTENTIAL *pot);
int SetPotentialSE(POTENTIAL *pot);
int SetPotentialVc(POTENTIAL *pot);
int SetPotentialVT(POTENTIAL *pot);
int SetPotentialU(POTENTIAL *pot, int n, double *u);
int SetPotentialW (POTENTIAL *pot, double e, int kappa, int kv);
int IdxVT(int kappa);
int DiracSmall(ORBITAL *orb, POTENTIAL *pot, int i2, int kv);
double EneTol(double e);
  
#endif


