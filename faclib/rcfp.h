#ifndef _RCFP_H_
#define _RCFP_H_ 1

/*************************************************************
  Header for module "rcfp".
  This module calculates the reduced coefficients of fractional
  parentage, the reduced matrix elements of creation, and
  annihilation operators between the single subshell states.

  It is mainly translated from the F90 package of 
  Gaigalas et al. CPC 139 (2001) 263.

  Author: M. F. Gu, mfgu@space.mit.edu
**************************************************************/

/* 
<** The following format is used for documenting the source **>
*/

/* documenting a struct */
/*
** STRUCT:      
** PURPOSE:     
** FIELDS:      
** NOTE:        
*/

/* documenting a function */
/* 
** FUNCTION:    
** PURPOSE:     
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/

/* documenting a macro function */
/* 
** MACRO:       
** PURPOSE:     
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/

/* documenting a global, static varialbe or a macro constant */
/*
** VARIABLE:    
** TYPE:        
** PURPOSE:     
** NOTE:        
*/

#include "global.h"
#include "angular.h"


/*
** STRUCT:      RCFP_STATE
** PURPOSE:     specify a jj coupled subshell state.
** FIELDS:      {int state},
**              index identifies a particular term in the 
**              RCFP_TERMS table which this state belongs to. 
**              {int n},
**              principal quantum number of the subshell.
**              {int nq},
**              occupation number of the subshell.
**              {int subshellMQ},
**              projection of the quasi-spin of the subshell.
** NOTE:        the angular momentum of the subshell j is 
**              related to subshellMQ and nq through:
**              2*j+1 = 2*(nq - 2*MQ)
*/
typedef struct _RCFP_STATE_ {
  int state; 
  int n; 
  int nq; 
  int subshellMQ; 
} RCFP_STATE;


/*
** STRUCT:      RCFP_TERM
** PURPOSE:     a particular term of a coupled subshell state
** FIELDS:      {int j},
**              angular momentum of the subshell.
**              {int Q},
**              quasi-spin of the subshell.
**              {int nu},
**              seneority of the state.
**              {int subshellJ},
**              total angular momentum of the coupled state.
**              {int Nr},
**              additional quantum numbers needed.
** NOTE:        
*/
typedef struct _RCFP_TERM_ {
  int j;
  int Q;
  int nu;
  int subshellJ;
  int Nr;
} RCFP_TERM;

/*
** STRUCT:      REDUCED_COEFF
** PURPOSE:     phase, numerator, and denominator of the coeff.
** FIELDS:      {int phase},
**              phase.
**              {int nom},
**              numerator.
**              {int denom},
**              denominator.
** NOTE:        coeff. = phase * nom/denom.
*/
typedef struct _REDUCED_COEFF_ {
  int phase;
  int nom;
  int denom;
} REDUCED_COEFF;

/* data structure for a general coupled operator */
/*
** STRUCT:      RCFP_OPERATOR
** PURPOSE:     a general coupled operator.
** FIELDS:      {int rank},
**              the shperical rank of the operator.
**              {int nops},
**              Num. of creation and annihilation operators
**              in the coupling.
**              {int qm},
**              the quasi-spin projection of the coupled operator.
**              {RCFP_OPERATOR *left},
**              pointer to the left component.
**              {RCFP_OPERATOR *right},
**              pointer to the right component.
** NOTE:        Such a general definition is not needed in actual
**              application. one only needs the operators of type
**              A, W=AxA, AxW, WxA, WxW, where A is a creation 
**              or annihilation operator. These cases are handled
**              by more specialized routines. This general 
**              definition exists just for the sake of completeness.
*/
typedef struct _RCFP_OPERATOR_ {
  int rank;
  int nops;
  int qm;
  struct _RCFP_OPERATOR_ *left;
  struct _RCFP_OPERATOR_ *right;
} RCFP_OPERATOR;


/* 
** public functions provided by "rcfp"
*/
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
int    QSpaceDelta(RCFP_STATE *sub_state);
int    RCFPTermIndex(int j, int nu, int Nr, int subshellJ);
double ClebschGordanQuasispin(int ja, int ma, int jb, int mb, 
			      int jab, int mab);
void   UnpackRCFPState(int s, int *j, int *nu, int *J);
int    PackRCFPState(int j, int nu, int J);
void   CoupleOperators(RCFP_OPERATOR *op1, RCFP_OPERATOR *op2, 
		       RCFP_OPERATOR *op, int rank);
double ReducedOperator(RCFP_STATE *bra, RCFP_STATE *ket, RCFP_OPERATOR *op);

#endif

