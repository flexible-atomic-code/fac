#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include "global.h"
#include "config.h"
#include "nucleus.h"
#include "radial.h"
#include "angular.h"
#include "rcfp.h"
#include "recouple.h"
#include <time.h>

#define MAX_ENERGY_CORRECTION 100
#define LEVELS_BLOCK    5000

typedef struct _HAMILTON_ {
  int pj;
  int dim;
  int n_basis;
  int hsize;
  int msize;
  int dim0;
  int n_basis0;
  int hsize0;
  int msize0;
  int *basis;
  double *hamilton;
  double *mixing;
  double *work;
  int *iwork;
} HAMILTON;

typedef struct _LEVEL_ {
  int pj;
  int n_basis;
  int *basis;
  double *mixing;
  int major_component;
  double energy;
} LEVEL;

#define ANGZ_BLOCK 4096
typedef struct _ANGZ_DATUM_ {
  void *angz;
  short nz;
} ANGZ_DATUM;

typedef struct _ANGULAR_ZMIX_ {
  double coeff;
  short k;
  short k0;
  short k1;
} ANGULAR_ZMIX;

typedef struct _ANGULAR_ZFB_ {
  double coeff;
  short kb;
} ANGULAR_ZFB;

typedef struct _ANGULAR_ZxZMIX_ {
  double coeff;
  short k;
  short k0;
  short k1;
  short k2;
  short k3;
} ANGULAR_ZxZMIX;

#ifdef PERFORM_STATISTICS
typedef struct _STRUCT_TIMING_ {
  clock_t angz_mix;
  clock_t angzxz_mix;
  clock_t angz_fb;
  clock_t angzxz_fb;
  clock_t angz_states;
  clock_t angzfb_states;
  clock_t angzxzfb_states;
  clock_t add_angz;
  clock_t add_angzxz;
  clock_t diag_ham;
  clock_t set_ham;
} STRUCT_TIMING;
int GetStructTiming(STRUCT_TIMING *t);
#endif

int ConstructHamilton(int isym, int k, int *kg, int kp, int *kgp);
int ValidBasis(STATE *s, int k, int *kg, int n);
int ConstructHamiltonFrozen(int isym, int k, int *kg, int n);
double HamiltonElement(int isym, int isi, int isj);
double HamiltonElementFrozen(int isym, int isi, int isj);
double Hamilton2E(int n_shells, SHELL_STATE *sbra, 
		  SHELL_STATE *sket,INTERACT_SHELL *s);
double Hamilton1E(int n_shells, SHELL_STATE *sbra, 
		  SHELL_STATE *sket,INTERACT_SHELL *s);

int DiagnolizeHamilton();
int AddToLevels();
int CorrectEnergy(int n, int *k, double *e);
LEVEL *GetLevel(int k);
int LevelTotalJ(int k);
int GetNumLevels();
int GetPrincipleBasis(double *mix, int d);
int SortLevels(int start, int n);
int GetBaseJ(STATE *s);

int AddToAngularZMix(int *n, int *nz, ANGULAR_ZMIX **ang,
		     int k, int k0, int k1, double coeff);
int AngularZMix(ANGULAR_ZMIX **ang, int lower, int upper, int mink, int maxk);
int CompareAngularZMix(const void *c1, const void *c2);
int AngularZFreeBound(ANGULAR_ZFB **ang, int lower, int upper);
int AngularZMixStates(ANGULAR_ZMIX **ang, STATE *s1, STATE *s2);
int AngZSwapBraKet(int nz, ANGULAR_ZMIX *ang, int p);
int AngularZFreeBoundStates(ANGULAR_ZFB **ang, STATE *slow, STATE *sup);
int AngularZxZMixStates(ANGULAR_ZxZMIX **ang, STATE *slow, STATE *sup);
int AngularZxZFreeBoundStates(ANGULAR_ZxZMIX **ang, STATE *slow, STATE *sup);
int AddToAngularZxZ(int *n, int *nz, ANGULAR_ZxZMIX **ang, 
		    int n_shells, int phase, SHELL_STATE *sbra, 
		    SHELL_STATE *sket, INTERACT_SHELL *s, int m);
int AddToAngularZxZMix(int *n, int *nz, ANGULAR_ZxZMIX **ang, 
		       int k, int k0, int k1, int k2, int k3, double r);
int AngularZxZFreeBound(ANGULAR_ZxZMIX **ang, int lower, int upper);

int GetBasisTable(char *fn);
int ConstructLevelName(char *name, char *sname, STATE *basis);
int SaveLevelsToAscii(char *fn, int m, int n);
int SetAngZOptions(int n, double mc, double c);
int FreeAngZ(int g, int which_array);
int InitStructure();
int ClearLevelTable();

#endif



