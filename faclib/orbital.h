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
  int maxrp;
  int nfrozen;
  int npseudo;
  int mpseudo;
  int dpseudo;
  int mode;
  int flag;
  int r_core;
  int nmax;
  int nse, mse, pse, mvp, pvp, hpvs, hlike;
  int ifm;
  double cfm, xfm;
  double hx0, hx1, chx;
  double hxs, ahx, ihx, rhx, dhx, bhx, ratio, asymp, rmin;
  double N; /*number of electrons*/
  double N1;
  double NC;
  double lambda; /* parameter for the Vc */
  double a; /* previously used to parameterize the Vc, but now always 0.0*/
  double ar, br, qr; /* parameter for the transformation */
  int ib, nb, ib1, ib0;
  double bqp, rb; /* boundary condition */
  double rfn[NKSEP];
  int nfn[NKSEP];
  double sturm_ene[NKSEP1];
  double sturm_idx;
  int nws;  
  double *dws;
  double zps, nps, tps, rps, dps, aps, bps, nbt, nbs, nqf;
  double fps, ups, xps, jps, qps, gps, sf0, sf1, efm, eth, ewd;
  int mps, kps, vxf, ips, sps, iqf, miter, vxm;
  double *Z, *dZ, *dZ2, *rad, *rho, *mqrho, *dr_drho, *dr_drho2, *vtr;
  double *a1r, *a2r, *a3r, *a1dr, *a2dr, *a3dr;
  double *Vc, *dVc, *dVc2, *qdist, *U, *dU, *dU2, *W, *dW, *dW2;
  double *ZVP, *dZVP, *dZVP2;
  double *NPS, *VPS, *EPS, *VXF, *ICF, *ZPS, *dZPS, *dZPS2;
  double *ZSE[NKSEP], *dZSE[NKSEP], *dZSE2[NKSEP];
  double *VT[NKSEP1], *dVT[NKSEP1], *dVT2[NKSEP1];
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
  double rfn, dn;
  struct _ORBITAL_ *horb;
  struct _ORBITAL_ *rorb;
} ORBITAL;

void InitOrbitalData(void *p, int n);
double *GetVEffective(void);
double RadialDiracCoulomb(int npts, double *p, double *q, double *r,
			  double z, int n, int kappa);
int RadialSolver(ORBITAL *orb,  POTENTIAL *pot);
int RadialBasis(ORBITAL *orb, POTENTIAL *pot);
int RadialSturm(ORBITAL *orb, POTENTIAL *pot);
int RadialBasisOuter(ORBITAL *orb, POTENTIAL *pot);
int RadialRydberg(ORBITAL *orb, POTENTIAL *pot);
int RadialBound(ORBITAL *orb, POTENTIAL *pot);
int RadialFreeInner(ORBITAL *orb, POTENTIAL *pot);
int RadialFree(ORBITAL *orb, POTENTIAL *pot);
int RadialFreeZMC(ORBITAL *orb, POTENTIAL *pot);
int CountNodes(double *p, POTENTIAL *pot, int i1, int i2);
double InnerProduct(int i1, int n, 
		    double *p1, double *p2, POTENTIAL *pot);
void Differential(double *p, double *dp, int i1, int i2, POTENTIAL *pot);
void DrLargeSmall(ORBITAL *orb, POTENTIAL *pot, double *pr, double *qr);
void InitializePS(POTENTIAL *pot);
int SetOrbitalRGrid(POTENTIAL *pot);
double GetRFromRho(double rho, double a, double b, double q, double r0);
int SetPotentialExtraZ(POTENTIAL *pot, int iep);
int SetPotentialZ(POTENTIAL *pot);
int SetPotentialVP(POTENTIAL *pot);
int SetPotentialSE(POTENTIAL *pot);
int SetPotentialPS(POTENTIAL *pot, int iter);
void FreeIonDensity(POTENTIAL *pot, double *de, double *di);
double FreeElectronDensity(POTENTIAL *pot, double e0, double u,
			   double zn, int md, int iter);
double FreeElectronIntegral(POTENTIAL *pot, int i0, int i1, int i2,
			    double *eps, double *icf,
			    double eth, double tps, double eref, double ewd,
			    double u, double zn, int md, int iter,
			    double ups, double fa, double rps,
			    double dps, double *aps);
double StewartPyattIntegrand(double x, double a, double fa, double y,
			     double y0, double g, double z, double nb);
double FermiDegeneracy(double ne, double te, double *yi);
double FermiIntegral(double x, double y, double g);
double InterpFermiRM1(double y, int m);
int SetPotentialVc(POTENTIAL *pot);
int SetPotentialVT(POTENTIAL *pot);
int SetPotentialU(POTENTIAL *pot, int n, double *u);
int SetPotentialW(POTENTIAL *pot, double e, int kappa, int kv);
int SetPotentialWZMC (POTENTIAL *pot, double e, int kappa, int kv);
int SetPotentialSturm(POTENTIAL *pot);
int AddPotentialSturm(POTENTIAL *pot, int k, double e);
int IdxVT(int kappa);
void MaxRGrid(POTENTIAL *pot, double gasymp, double gratio, double z,
	      double rmin, int maxrp, double *ap, double *cp, double *rmp);
int DiracSmall(ORBITAL *orb, POTENTIAL *pot, int i2, int kv);
int DiracSmallZMC(ORBITAL *orb, POTENTIAL *pot, int i2, int kv);
double EneTol(double e);
void SetOrbitalWorkSpace(double *p, int n);
double EnergyH(double z, double n, int ka);
double QuantumDefect(double z, int n, int ka, double e);
void SetOptionOrbital(char *s, char *sp, int ip, double dp);
void PrepFermiRM1(double a, double fa, double t);
double SetSPZW(int n, double *zw);
void IonCF(POTENTIAL *pot, double *y, double rmax, int n, double *w);
void SetPotentialIPS(POTENTIAL *pot, int iter);
int SPMode();
int OnErrorOrb();
void SetOnErrorOrb(int e);
double ZColl();
double MColl();
void SetZColl(double z);
void SetMColl(double z);
int DensityToSZ(POTENTIAL *pot, double *d,
		double *z, double *zx, double *jps, int md);
double XCPotential(double t, double n, int md);
double SCRSF(void);
double SCRBF(void);
double SCBQP(void);
double SCEWF(void);
int SCEWM(void);
int OrbPWA(void);
double SPZU(void);
int IonSphBMode(void);
void PrepFermiNR();
double FreeEta0(double ne, double te);
double FreeEta1(double ne, double ke);
double FreeEta2(double ne, double ve);
double FreeKe(double u, double te);
double FreeVe(double u, double te);
double FreeTe(double u, double ne);
double InterpFermiNR(int m, double u);
double BoundFactor(double e, double eth, double de, int mps);
double DrRho(POTENTIAL *p, int i, double x);
double RadRho(POTENTIAL *p, int i, double x);
double ExpM1(double x);
int FirstMaximum(double *p, int i1, int i2, POTENTIAL *pot);
int LastMaximum(double *p, int i1, int i2, POTENTIAL *pot);
#endif


