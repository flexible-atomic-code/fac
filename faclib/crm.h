#ifndef _CRM_H_
#define _CRM_H_ 1

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "array.h"
#include "dbase.h"
#include "rates.h"

#define RATES_BLOCK   1024
#define ION_BLOCK     16
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
  RECOMBINED *rec;
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
  int nlevels;
  LBLOCK **iblock;
  int *ilev;
  short *j;
  double *energy;
  ARRAY *ce_rates;
  ARRAY *tr_rates;
  ARRAY *ci_rates;
  ARRAY *rr_rates;
  ARRAY *ai_rates;
  ARRAY *recombined;
  int nele;
  char *dbfiles[NDB];
  double n;
  NCOMPLEX ce_max;
} ION;

typedef struct _IONIZED_ {
  int nele;
  char symbol[4];
  int atom;
  char *dbfiles[NDB];
  int nionized;
  int *ionized_map[2];
  double *energy;
  double n;
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
void ExtrapolateTR(ION *ion, int inv);
void ExtrapolateRR(ION *ion, int inv);
void ExtrapolateAI(ION *ion, int inv);
int SetBlocks(double ni, char *ifn);
int SetAbund(int nele, double abund);
int InitBlocks(void);
void AddRate(ION *ion, ARRAY *rts, RATE *r, int m);
int SetCERates(int inv);
int SetTRRates(int inv);
int SetCIRates(int inv);
int SetRRRates(int inv);
int SetAIRates(int inv);
int RateTable(char *fn);
int BlockMatrix(void);
int BlockPopulation(void);
double BlockRelaxation(int iter);
int LevelPopulation(void);
int Cascade(void);
int SpecTable(char *fn, double smin);
int SelectLines(char *ifn, char *ofn, int nele, int type, 
		double emin, double emax, double fmin);
int PlotSpec(char *ifn, char *ofn, int nele, int type,
	     double emin, double emax, double de, double smin);

#endif

