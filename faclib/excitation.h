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

#ifndef _EXCITATION_H_
#define _EXCITATION_H_

#include "coulomb.h"
#include "structure.h"
#include "transition.h"

typedef struct _CEPW_SCRATCH_ {
  int qr;
  int max_kl;
  int kl_cb;
  double tolerance;
  int nkl0;
  int nkl;
  int ns;
  double kl[MAXNKL+1];
  double log_kl[MAXNKL];
} CEPW_SCRATCH;

typedef struct _CEPK_ {
  short nkl;
  short nkappa;
  short *kappa0;
  short *kappa1;
  double *pkd;
  double *pke;
} CEPK;

typedef struct _CEPKK_ {
  short kmin, kmax;
  CEPK **pk;
} CEPKK;

typedef struct _CEQKK_ {
  short kmin, kmax;
  short kminp, kmaxp;
  double **qk;
} CEQKK;

typedef struct _CECACHE_ {
  int m;
  int nc;
  int *low;
  int *up;
  int *nz;
  ANGULAR_ZMIX **az;
  int *nmk;
  double **mbk;
  IDXARY *ks, *kd0, *kd1;
  int msub;
  CEPKK *pk;
  CEQKK *qk;
} CECACHE;

#ifdef PERFORM_STATISTICS
typedef struct _EXCIT_TIMING_ {
  double rad_pk;
  double rad_qk;
  double set_kappa;
} EXCIT_TIMING;
#endif

CEPW_SCRATCH *GetCEPWScratch(void);
void SetMaxCECache(int n);
void AllocCECache(int msub);
void FreeCECache(int m);
void FreeCEQKK(CEQKK *qk, int m);
void FreeCEPKK(CEPKK *pk);
int FreeExcitationQk(void);
int InitExcitation(void);
int ReinitExcitation(int m);
int SetCETEGrid(int n, double emin, double emax);
int SetCETEGridDetail(int n, double *x);
int SetAngleGrid(int m, int n, double xmin, double xmax);
int SetAngleGridDetail(int m, int n, double *xg);
int SetCEPWFile(char *fn);
int SetCEBorn(double e, double x, double x1, double x0);
void SetCELQR(int m);
void SetCELMax(int m);
void SetCELCB(int m);
void SetCETol(double t);
int SetCEPWOptions(int qr, int max, int kl_cb, double tol);
int AddCEPW(int n, int step);
int SetCEFormat(int m);
int SetCEQkMode(int m, double tol);
int SetCEPWGrid(int ns, int *n, int *step);
int SetCEEGridLimits(double min, double max, int type);
int SetCEEGridType(int type);
int SetUsrCEEGridType(int type);
int SetCEPWGridType(int type);
int SetCEEGridDetail(int n, double *x);
int SetCEEGrid(int n, double emin, double emax, double eth);
int SetUsrCEEGridDetail(int n, double *x);
int SetUsrCEEGrid(int n, double emin, double emax, double eth);

int CERadialPk(CEPK **pk, int ie, int k0, int k1, int k, int trylock);
int CERadialQkBorn(int k0, int k1, int k2, int k3, int k, 
		   double te, double e1, double *qk, int m);
int CERadialQkBornMSub(int k0, int k1, int k2, int k3, int k, int kp,
		       double te, double e1, 
		       int nq, int *q, double *qk, int m);
double *CERadialQkTable(int k0, int k1, int k2, int k3, int k, int trylock);
double *CERadialQkMSubTable(int k0, int k1, int k2, int k3, int k, int kp,
			    int trylock);
int CERadialQk(double *r, double te, 
	       int k0, int k1, int k2, int k3, int k, int trylock);
int CERadialQkMSub(double *rq, double te, int k0, int k1,
		   int k2, int k3, int k, int kp, int trylock);
void CERadialQkFromFit(int np, double *p, int n, double *x, double *logx,
		       double *y, double *dy, int ndy, void *extra);
int CollisionStrengthUTA(double *qkt, double *params, double *e, double *bethe,
			 int lower, int upper);
int CollisionStrength(double *s, double *p, double *e, double *bethe,
		      int lower, int upper, int msub);
void ProcessCECache(int msub, int iuta, TFILE *f);
int SaveExcitation(int nlow, int *low, int nup, int *up, int msub, char *fn);
int CollisionStrengthEB(double *s, double *e, double *bethe, int lower, int upper);
int CollisionStrengthEBD(double *s, double *e, double *bethe, double *born,
			 int lower, int upper);
int SaveExcitationEB(int nlow, int *low, int nup, int *up, char *fn);
int SaveExcitationEBD(int nlow, int *low, int nup, int *up, char *fn);

#endif
