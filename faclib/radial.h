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

#ifndef _RADIAL_H_
#define _RADIAL_H_

#include <time.h>
#include <math.h>

#include "global.h"
#include "nucleus.h"
#include "interpolation.h"
#include "grid.h"
#include "coulomb.h"
#include "orbital.h"
#include "config.h"
#include "angular.h"
#include "recouple.h"

typedef struct _SLATER_YK_ {
  short npts;
  float *yk;
  float coeff[2];
} SLATER_YK;

#ifdef PERFORM_STATISTICS
typedef struct _RAD_TIMING_ {
  double radial_1e;
  double radial_2e;
  double dirac;
  double radial_slater;
} RAD_TIMING;

int GetRadTiming(RAD_TIMING *t);
#endif

double *WLarge(ORBITAL *orb);
double *WSmall(ORBITAL *orb);
int GetBoundary(double *rb, double *b, int *nmax, double *dr);
int SetBoundary(int nmax, double p, double bqp);
void PrintQED();
int RadialOverlaps(char *fn, int kappa);
void SetSlaterCut(int k0, int k1);
void SetPotentialMode(int m, double h);
void SetSE(int n, int m);
void SetVP(int n);
void SetBreit(int n, int m);
void SetMS(int nms, int sms);
int SetAWGrid(int n, double min, double max);
int GetAWGrid(double **a);
int SetRadialGrid(int maxrp, double ratio, double asymp, double rmin);
double SetPotential(AVERAGE_CONFIG *acfg, int iter);
POTENTIAL *RadialPotential(void);
int GetPotential(char *s);
double GetResidualZ(void);
double GetRMax(void);

/* solve the dirac equation for the given orbital */
int SolveDirac(ORBITAL *orb);
int WaveFuncTable(char *s, int n, int kappa, double e);

/* get the index of the given orbital in the table */
int OrbitalIndex(int n, int kappa, double energy);
int OrbitalExists(int n, int kappa, double energy);
int AddOrbital(ORBITAL *orb);
ORBITAL *GetOrbital(int k);
ORBITAL *GetOrbitalSolved(int k);
ORBITAL *GetNewOrbital(void);
int GetNumBounds(void);
int GetNumOrbitals(void);
int GetNumContinua(void);

double CoulombEnergyShell(CONFIG *cfg, int i);
void ShiftOrbitalEnergy(CONFIG *cfg);
double GetPhaseShift(int k);

/* radial optimization */
int SetAverageConfig(int nshells, int *n, int *kappa, double *nq);
void SetOptimizeMaxIter(int m);
void SetOptimizeStabilizer(double m);
void SetOptimizeTolerance(double c);
void SetOptimizePrint(int m);
void SetOptimizeControl(double tolerence, double stablizer, 
			int maxiter, int iprint);
void SetScreening(int n_screen, int *screened_n, 
		  double screened_harge, int kl);
int OptimizeRadial(int ng, int *kg, double *weight);
int RefineRadial(int maxfun, int msglvl);
double ConfigEnergyShiftCI(int nrs0, int nrs1);
double ConfigEnergyShift(int ns, SHELL *bra, int ia, int ib, int m2);
double ConfigEnergyVariance(int ns, SHELL *bra, int ia, int ib, int m2);
int ConfigEnergy(int m, int mr, int ng, int *kg);
double TotalEnergyGroup(int kg);
double ZerothEnergyConfig(CONFIG *cfg);
double ZerothResidualConfig(CONFIG *cfg);
double AverageEnergyConfig(CONFIG *cfg);
double AverageEnergyAvgConfig(AVERAGE_CONFIG *cfg);
void DiExAvgConfig(AVERAGE_CONFIG *cfg, double *d0, double *d1);

/* routines for radial integral calculations */
int GetYk(int k, double *yk, ORBITAL *orb1, ORBITAL *orb2, 
	  int k1, int k2, int type);
int Integrate(double *f, ORBITAL *orb1, ORBITAL *orb2, int type, double *r, int id);
int IntegrateSubRegion(int i0, int i1, 
		       double *f, ORBITAL *orb1, ORBITAL *orb2,
		       int t, double *r, int m, double *ext);
int IntegrateSinCos(int j, double *x, double *y, 
		    double *phase, double *dphase, 
		    int i0, double *r, int t, double *ext);
int SlaterTotal(double *sd, double *se, int *js, int *ks, int k, int mode);
double Vinti(int k0, int k1);
double QED1E(int k0, int k1);
double SelfEnergyRatio(ORBITAL *orb);
int Slater(double *s, int k0, int k1, int k2, int k3, int k, int mode);
void BreitX(ORBITAL *orb0, ORBITAL *orb1, int k, int m, double e, double *r);
double BreitC(int n, int m, int k, int k0, int k1, int k2, int k3);
double BreitS(int k0, int k1, int k2, int k3, int k);
double BreitI(int n, int k0, int k1, int k2, int k3, int m);
double Breit(int k0, int k1, int k2, int k3, int k,
	     int kp0, int kp1, int kp2, int kp3,
	     int kl0, int kl1, int kl2, int kl3);
void SortSlaterKey(int *kd);
void PrepSlater(int ib0, int iu0, int ib1, int iu1,
		int ib2, int iu2, int ib3, int iu3);
int ResidualPotential(double *s, int k0, int k1);
double MeanPotential(int k0, int k1);
int FreeResidualArray(void);
int FreeMultipoleArray(void);
int FreeSlaterArray(void);
int FreeSimpleArray(MULTI *ma);
int FreeMomentsArray(void);
int FreeGOSArray(void);

double RadialMoments(int m, int k1, int k2);
double MultipoleRadialNR(int m, int k1, int k2, int guage);
int MultipoleRadialFRGrid(double **p, int m, int k1, int k2, int guage);
double MultipoleRadialFR(double aw, int m, int k1, int k2, int guage);
double InterpolateMultipole(double aw2, int n, double *x, double *y);
double *GeneralizedMoments(int k0, int k1, int m);
void PrintGeneralizedMoments(char *fn, int m, int n0, int k0, int n1, int k1, 
			     double e1);
int SaveOrbital(int i);
int RestoreOrbital(int i); 
int FreeOrbital(int i);
int SaveAllContinua(int mode); 
int SaveContinua(double e, int mode);
int FreeAllContinua(void);
int FreeContinua(double e);
int ClearOrbitalTable(int m);
void LimitArrayRadial(int m, double n);
int InitRadial(void);
int ReinitRadial(int m);
int TestIntegrate(void);

#endif

