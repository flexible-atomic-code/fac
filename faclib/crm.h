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

#ifndef _CRM_H_
#define _CRM_H_ 1

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "array.h"
#include "dbase.h"
#include "rates.h"
#include "nucleus.h"
#include "interpolation.h"
#include "coulomb.h"

#define RATES_BLOCK   1024
#define ION_BLOCK     4
#define LBLOCK_BLOCK  1024

#define MAXNREC 128
typedef struct _RECOMBINED_ {
  int bmin, bmax;
  int n, n_ext;
  int imin[MAXNREC];
  int imax[MAXNREC];
  int nrec[MAXNREC];
} RECOMBINED;

typedef struct _NCOMPLEX_ {
  short n;
  short nq;
} NCOMPLEX;

#define MAXNCOMPLEX 8
typedef struct _LBLOCK_ {
  int ib;
  int iion;
  int irec;
  int ionized;
  RECOMBINED *rec;
  int imin;
  int nlevels;
  double nb;
  double *n, *n0;
  double *r;
  double *total_rate;
  NCOMPLEX ncomplex[MAXNCOMPLEX];
} LBLOCK;

typedef struct _BLK_RATE_ {
  LBLOCK *iblock;
  LBLOCK *fblock;
  ARRAY *rates;
} BLK_RATE;

typedef struct _ION_ {
  int iground; /* ionized ground state of this ion */
  int nlevels;
  LBLOCK **iblock;
  int *ilev;
  int *j;
  short *vnl;
  short *ibase;
  double *energy;
  ARRAY *ce_rates;
  ARRAY *tr_rates;
  ARRAY *tr_sdev;
  ARRAY *tr2_rates;
  ARRAY *ci_rates;
  ARRAY *rr_rates;
  ARRAY *ai_rates;
  ARRAY *recombined;
  int nele;
  char *dbfiles[NDB];
  double n, nt, n0;
  NCOMPLEX ce_max;
  int KLN_min, KLN_max;
  int KLN_bmin, KLN_bmax;
  int KLN_amin, KLN_amax;
  double *KLN_ai;
  int *KLN_nai;
  double ace, atr, aci, arr, aai;
} ION;

typedef struct _IONIZED_ {
  int nele;
  char symbol[4];
  double atom;
  char *dbfiles[NDB];
  int nionized;
  int *ionized_map[2];
  int imin[2], imax[2];
  double *energy;
  double n, nt, n0;
  double ace, atr, aci, arr, aai;
} IONIZED;

typedef struct _RATE_ {
  int i, f;
  double dir;
  double inv;
} RATE;

int SetNumSingleBlocks(int n);
int SetEleDensity(double ele);
int SetPhoDensity(double pho);
int SetCascade(int c, double a);
int SetIteration(double acc, double s, int max);
int InitCRM(void);
int ReinitCRM(int m);
int AddIon(int nele, double n, char *pref);
int IonIndex(ION *ion, int i, int k);
int IonizedIndex(int i, int m);
int FindLevelBlock(int n, EN_RECORD *r0, EN_RECORD *r1, 
		   int nele, char *ifn);
void GetRecombined(int *b, int *nrec, char *name);
int CopyNComplex(NCOMPLEX *dest, NCOMPLEX *src);
int GetNComplex(NCOMPLEX *c, char *s);
int CompareNComplex(NCOMPLEX *c1, NCOMPLEX *c2);
int StrNComplex(char *s, NCOMPLEX *c);
int TransitionType(NCOMPLEX *ic, NCOMPLEX *fc);
void ExtrapolateEN(int i, ION *ion);
int FindFinalTR(ION *ion, int f, int n1, int n0);
void ExtrapolateTR(ION *ion, int inv, int **irb);
void ExtrapolateRR(ION *ion, int inv, int **irb);
void ExtrapolateAI(ION *ion, int inv, int **irb);
int SetBlocks(double ni, char *ifn);
void SetRateMultiplier(int nele, int t, double a);
int SetAbund(int nele, double abund);
int InitBlocks(void);
int AddRate(ION *ion, ARRAY *rts, RATE *r, int m, int **irb);
int SetCERates(int inv);
int SetTRRates(int inv);
int SetCIRates(int inv);
int SetRRRates(int inv);
int SetAIRates(int inv);
int SetAIRatesInner(char *fn);
int RateTable(char *fn, int nc, char *sc[], int md);
int BlockMatrix(void);
int BlockPopulation(int n);
double BlockRelaxation(int iter);
int LevelPopulation(void);
int Cascade(void);
int SpecTable(char *fn, int rrc, double smin);
int SelectLines(char *ifn, char *ofn, int nele, int type, 
		double emin, double emax, double fmin);
int PlotSpec(char *ifn, char *ofn, int nele, int type,
	     double emin, double emax, double de, double smin);
int DRBranch(void);
int DRStrength(char *fn, int nele, int mode, int ilev0);
int DumpRates(char *fn, int k, int m, int imax, int a);
int ModifyRates(char *fn);
int SetInnerAuger(int i);
int SetExtrapolate(int e);
void TabNLTE(char *fn1, char *fn2, char *fn3, char *fn,
	     double xmin, double xmax, double dx);
int SetEMinAI(double e);
int DRSuppression(char *fn, double z, int nmax);
int RydBranch(char *fn, char *ofn, int n0, int n1);
int NormalizeMode(int m);
void FixNorm(int m);

ARRAY* _GetIons();  // Add an access to ions for testing purpose
#endif

