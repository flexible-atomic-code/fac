#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include "global.h"
#include "config.h"
#include "nucleus.h"
#include "radial.h"
#include "angular.h"
#include "dbase.h"
#include "rcfp.h"
#include "recouple.h"
#include <time.h>

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
  int pb;
  int kpb[NPRINCIPLE];
  int ibase;
  int *basis;
  double *mixing;
  double energy;
  int ngp, *igp;
} LEVEL;

typedef struct _LEVEL_ION_ {
  int imin;
  int imax;
} LEVEL_ION;

typedef struct _ANGZ_DATUM_ {
  int ns;
  int *ic;
  int *nz;
  void **angz;
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

typedef struct _ANGULAR_FROZEN_ {
  int nts, ncs;
  int *ts, *cs;
  int *nz, *nzfb, *nzxzfb;
  ANGULAR_ZMIX **z;
  ANGULAR_ZFB **zfb;
  ANGULAR_ZxZMIX **zxzfb;
} ANGULAR_FROZEN;

typedef struct _ECORRECTION_ {
  int iref;
  int ilev;
  double e;
  int nmin;
  STATE *s;
} ECORRECTION;

typedef struct _CORR_CONFIG_ {
  CONFIG *c;
  int inp, inq, np, nq;
  int ig, ncs, kp, kq;
} CORR_CONFIG;

#ifdef PERFORM_STATISTICS
typedef struct _STRUCT_TIMING_ {
  double angz_mix;
  double angzxz_mix;
  double angz_fb;
  double angzxz_fb;
  double angz_states;
  double angz_states_load;
  long n_angz_states;
  long n_angz_states_load;
  double angzfb_states;
  double angzxzfb_states;
  double add_angz;
  double add_angzxz;
  double diag_ham;
  double set_ham;
} STRUCT_TIMING;
int GetStructTiming(STRUCT_TIMING *t);
#endif

double MBPT0(int isym, SYMMETRY *sym, int q1, int q2, int kg, int ic);
double MBPT1(LEVEL *lev, SYMMETRY *sym, int kg, int m);
void MBPT2(int ilev, LEVEL *lev, int kg, int n0, int n1, double *de, int nk);
int MBPT(char *fn, int n, int *s, int k, int *kg, int *n0, int n1, 
	 int kmax, int nt, int m);
int MBPTS(char *fn, char *fn1, int n, int *s, int k, int *kg,
	  int *n0, int nmax, int kmax, int nt);
int StructureMBPT(char *fn, char *fn1, int n, int *s0, int k, int *kg,
		  int *n0, int *ni, int nmax, int kmax, int nt, 
		  int n2, int nt2, char *gn0, double eps, double eps1);
int IBisect(int k, int n, int *a);
int ConstructHamilton(int isym, int k0, int k, int *kg, int kp, int *kgp);
int ConstructHamiltonDiagonal(int isym, int k, int *kg);
int ValidBasis(STATE *s, int k, int *kg, int n);
int ConstructHamiltonFrozen(int isym, int k, int *kg, int n, int nc, int *kc);
double HamiltonElement(int isym, int isi, int isj);
double HamiltonElementFrozen(int isym, int isi, int isj);
double HamiltonElementFB(int isym, int isi, int isj);
double Hamilton2E2(int n_shells, SHELL_STATE *sbra, 
		   SHELL_STATE *sket,INTERACT_SHELL *s);
double Hamilton2E(int n_shells, SHELL_STATE *sbra, 
		  SHELL_STATE *sket,INTERACT_SHELL *s);
double Hamilton1E(int n_shells, SHELL_STATE *sbra, 
		  SHELL_STATE *sket,INTERACT_SHELL *s);
HAMILTON *GetHamilton(void);
int DiagnolizeHamilton(void);
int AddToLevels(int ng, int *kg);
int AddECorrection(int kref, int k, double e, int nmin);
LEVEL *GetLevel(int k);
int LevelTotalJ(int k);
int GetNumLevels(void);
int GetNumElectrons(int k);
int SortMixing(int start, int n, int *basis, double *mix, SYMMETRY *sym);
int GetPrincipleBasis(double *mix, int d, int *kpb);
int CompareLevels(LEVEL *lev1, LEVEL *lev2);
int SortLevels(int start, int n);
int GetBaseJ(STATE *s);
void AngularFrozen(int nts, int *ts, int ncs, int *cs);
void ClearAngularFrozeb(void);
int AngularZMix(ANGULAR_ZMIX **ang, int lower, int upper, int mink, int maxk);
int CompareAngularZMix(const void *c1, const void *c2);
int CompareAngularZxZMix(const void *c1, const void *c2);
int CompareAngularZFB(const void *c1, const void *c2);
int PackAngularZxZMix(int *n, ANGULAR_ZxZMIX **ang, int nz);
int PackAngularZMix(int *n, ANGULAR_ZMIX **ang, int nz);
int PackAngularZFB(int *n, ANGULAR_ZFB **ang, int nz);
int AngularZFreeBound(ANGULAR_ZFB **ang, int lower, int upper);
int AngularZMixStates(ANGZ_DATUM **ad, 
		      int kg1, int kg2, 
		      int kp1, int kp2);
int AngZSwapBraKet(int nz, ANGULAR_ZMIX *ang, int p);
int AngularZFreeBoundStates(ANGZ_DATUM **ad, 
			    int kg1, int kg2, int kc1, int kc2);
int AngularZxZMixStates(ANGZ_DATUM **ad, 
			int kg1, int kg2, 
			int kp1, int kp2);
int AngularZxZFreeBoundStates(ANGZ_DATUM **ad, 
			      int kg1, int kg2, 
			      int kp1, int kp2);
int AddToAngularZxZ(int *n, int *nz, ANGULAR_ZxZMIX **ang, 
		    int n_shells, int phase, SHELL_STATE *sbra, 
		    SHELL_STATE *sket, INTERACT_SHELL *s, int m);
int AddToAngularZxZMix(int *n, int *nz, ANGULAR_ZxZMIX **ang, 
		       int k, int k0, int k1, int k2, int k3, double coeff);
int AddToAngularZMix(int *n, int *nz, ANGULAR_ZMIX **ang,
		     int k, int k0, int k1, double coeff);
int AddToAngularZFB(int *n, int *nz, ANGULAR_ZFB **ang,
		    int kb, double coeff);
int AngularZxZFreeBound(ANGULAR_ZxZMIX **ang, int lower, int upper);

int GetBasisTable(char *fn);
int ConstructLevelName(char *name, char *sname, char *nc, 
		       int *vnl, STATE *basis);
int SaveLevels(char *fn, int m, int n);
int SetAngZOptions(int n, double mc, double c);
int SetAngZCut(double c);
int SetMixCut(double c);
int FreeAngZ(int g, int which_array);
int ClearLevelTable(void);
int InitStructure(void);
int ReinitStructure(int m);
int TestHamilton(void);

#endif



