#ifndef _CRM_H_
#define _CRM_H_ 1

#include <stdlib.h>
#include <stdio.h>
#include "array.h"
#include "dbase.h"
#include "rates.h"

#define RATES_BLOCK   1024
#define ION_BLOCK     16
#define LBLOCK_BLOCK  1024

typedef struct _LBLOCK_ {
  int iion;
  int nlevels;
  double nb;
  double *n, *n0;
  double *r;
  double *total_rate;
} LBLOCK;

typedef struct _ION_ {
  int nlevels;
  int *iblock;
  int *ilev;
  short *j;
  double *energy;
  ARRAY *ce_rates;
  ARRAY *tr_rates;
  ARRAY *ci_rates;
  ARRAY *rr_rates;
  ARRAY *ai_rates;
  int nele;
  char *dbfiles[NDB];
  double n;
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
int SetIteration(double acc, double s, int max);
int InitCRM(void);
int ReinitCRM(int m);
int AddIon(int nele, double n, char *pref);
int IonizedIndex(int i, int m);
int FindLevelBlock(int n, EN_RECORD *r0, EN_RECORD *r1, 
		   int nele, char *ifn);
int SetBlocks(double ni, char *ifn);
int SetAbund(int nele, double abund);
int InitBlocks(void);
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
int SpecTable(char *fn, double smin);
int SelectLines(char *ifn, char *ofn, int type, 
		double emin, double emax);
int PlotSpec(char *ifn, char *ofn, int type,
	     double emin, double emax, double de, double smin);

#endif

