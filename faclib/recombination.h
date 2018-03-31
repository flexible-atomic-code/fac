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

#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "global.h"
#include "nucleus.h"
#include "angular.h"
#include "config.h"
#include "structure.h"
#include "excitation.h"
#include "transition.h"
#include "coulomb.h"

#define MAX_COMPLEX 512
typedef struct _REC_COMPLEX_ {
  int n;
  ARRAY *rg;
  int s0;
  int s1;
} REC_COMPLEX;

typedef struct _AICACHE_ {
  int nc;
  int *nz, *nzf;
  ANGULAR_ZxZMIX **az;
  ANGULAR_ZFB **azf;
  int *low;
  int *up;
  double *e;
} AICACHE;

int InitRecombination(void);
int ReinitRecombination(int m);
void SetMaxAICache(int n);
void AllocAICache(void);
void FreeAICache(int m);
int FreeRecPk(void);
int FreeRecQk(void);
int FreeRecAngZ(void);
int SetAICut(double c);
int SetRRTEGrid(int n, double emin, double emax);
int SetRRTEGridDetail(int n, double *x);
int SetPEGrid(int n, double emin, double emax, double eth);
int SetPEGridDetail(int n, double *x);
int SetPEGridLimits(double min, double max, int type);
int SetUsrPEGridType(int type);
int SetUsrPEGrid(int n, double emin, double emax, double eth);
int SetUsrPEGridDetail(int n, double *x);
int AddRecPW(int n, int step);
int SetRecQkMode(int m, double tol);
int SetRecPWOptions(int kl_interp, int max_kl);
int SetRecPWLimits(int m1, int m2);
int SetRecSpectator(int n_max, int n_frozen, int n_spec);
int ConstructRecGroupName(char *rgn, char *gn, int n);
int RecStates(int n, int k, int *kg, char *fn);
int RecStatesFrozen(int n, int k, int *kg, char *fn);
int RRRadialMultipoleTable(double *qr, int k0, int k1, int m);
int RRRadialQkTable(double *qr, int k0, int k1, int m);
int RRRadialMultipole(double *rqc, double te, int k0, int k1, int m);
int RRRadialQk(double *rqc, double te, int k0, int k1, int m);
void RRRadialQkFromFit(int np, double *p, int n, double *x, double *logx, 
		       double *y, double *dy, int ndy, void *extra);
void RRRadialQkHydrogenicParams(int np, double *p, double z, int n, int klb);
double PICrossH(double z, int n0, int kl0, double e, int os);
double RRCrossH(double z, int n0, int kl0, double e);
int BoundFreeMultipole(FILE *fp, int rec, int f, int m);
int BoundFreeOS(double *rqu, double *p, 
		double *eb, int rec, int f, int m);
int BoundFreeOSUTA(double *rqu, double *rqc, double *eb, 
		   int rec, int f, int m);
int PrepRREGrids(double eth, double emax0);
int SaveRRMultipole(int nlow, int *low, int nup, int *up, char *fn, int m);
int SaveRecRR(int nlow, int *low, int nup, int *up, char *fn, int m);
int SaveAI(int nlow, int *low, int nup, int *up, char *fn, 
	   double eref, int msub);
int AsymmetryPI(int k0, double e, int mx, int m, double *b);
int SaveAsymmetry(char *fn, char *s, int mx);
int AIRadial1E(double *pk, int kb, int kappaf);
int AIRadialPk(double **pk, int k0, int k1, int kb, int kappaf,
	       int k, int trylock);
int AutoionizeRateUTA(double *rate, double *e, int rec, int f);
int AutoionizeRate(double *rate, double *e, int rec, int f, int msub);
void ProcessAICache(int msub, int iuta, TFILE *f);
int DROpen(int n, int *nlev, int **ops);
 
#endif
