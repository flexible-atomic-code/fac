#ifndef _RCFP_H_
#define _RCFP_H_

#include "global.h"
#include "angular.h"

/* this specifies a subshell. the subshellMQ value determines
   its angular momentum */
typedef struct _RCFP_STATE_ {
  int state; /* identifier */
  int n; /* principle quantum number */
  int nq; /* occupation number */
  int subshellMQ; /*projection of quasi-spin */
} RCFP_STATE;

typedef struct _RCFP_TERM_ {
  int j; /* angular momentum of the subshell */
  int Q; /* quasi-spin */
  int nu; /* seniority */
  int subshellJ; /* total angular momentum */
  int Nr; /* additional quantum numbers */
} RCFP_TERM;

/* reduced coeff. of fractional parentage. the value is 
   (-1)^phase * nom/denom */
typedef struct _REDUCED_COEFF_ {
  int phase;
  int nom;
  int denom;
} REDUCED_COEFF;

/* data structure for a general coupled operator */
typedef struct _RCFP_OPERATOR_ {
  int rank;
  int nops;
  int qm;
  struct _RCFP_OPERATOR_ *left;
  struct _RCFP_OPERATOR_ *right;
} RCFP_OPERATOR;


double ReducedCFP(int no_bra, int no_ket);

double CompleteReducedW(int no_bra, int no_ket, int k_q, int k_j);
double CompleteReducedWFromTable(int no_bra, int no_ket, int k1, int k2);
double CompleteReducedWAll(REDUCED_COEFF *w1, REDUCED_COEFF *w3,
			   REDUCED_COEFF *w5, REDUCED_COEFF *w7,
			   int no_bra, int no_ket);
double CompleteReducedW3(REDUCED_COEFF *w3_o, REDUCED_COEFF *w5_e,
			 int no_bra, int no_ket);
double CompleteReducedW5(REDUCED_COEFF *w5_o, REDUCED_COEFF *w5_e,
			 int no_bra, int no_ket);
double CompleteReducedW7(REDUCED_COEFF *w7_o, REDUCED_COEFF *w7_e,
			 int no_bra, int no_ket);

double ReducedW(RCFP_STATE *bra, RCFP_STATE *ket, 
		int k_j, int q_m1, int q_m2);
double ReducedWxW0(RCFP_STATE *bra, RCFP_STATE *ket,
		   int k_j, int q_m1, 
		   int q_m2, int q_m3, int q_m4);
double ReducedAxW(RCFP_STATE *bra, RCFP_STATE *ket,
		  int k_j1, int kk_j2, int q_m1, 
		  int q_m2, int q_m3);
double ReducedWxA(RCFP_STATE *bra, RCFP_STATE *ket,
		   int k_j1, int kk_j2, int q_m1,
		  int q_m2, int q_m3);
double ReducedA(RCFP_STATE *bra, RCFP_STATE *ket, int q_m);

int QSpaceDelta(RCFP_STATE *sub_state);
int RCFPTermIndex(int j, int nu, int Nr, int subshellJ);
double ClebschGordanQuasispin(int ja, int ma, int jb, int mb, 
			      int jab, int mab);
void UnpackRCFPState(int s, int *j, int *nu, int *J);
int PackRCFPState(int j, int nu, int J);

void CoupleOperators(RCFP_OPERATOR *op1, RCFP_OPERATOR *op2, 
		     RCFP_OPERATOR *op, int rank);
double ReducedOperator(RCFP_STATE *bra, RCFP_STATE *ket, RCFP_OPERATOR *op);

#endif

