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

#ifndef _RECOUPLE_H_
#define _RECOUPLE_H_ 1

/*************************************************************
  Header for module "recouple".
  This module calculates the recoupling coefficients. 

  The main task is to determine which electrons are the 
  interacting ones, and calculate the reduced matrix elements
  of the operator Z and ZxZ0, 

  Author: M. F. Gu, mfgu@stanford.edu
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

#include "config.h"
#include "angular.h"
#include "rcfp.h"
#include <time.h>

/*
** STRUCT:      INTERACT_SHELL
** PURPOSE:     stores information of an itneracting shell.
** FIELDS:      {int index},
**              the index of the shell within the SHELL_STATE.
**              {int n, j, kl, kappa},
**              the principle quantum number, the angular momentum,
**              the orbital angular momentum, and the relativistic 
**              angular quantum number of the shell.
**              {int nq_bra, nq_ket}, 
**              the occupation numbers int the bra and the ket states.
** NOTE:        
*/
typedef struct _INTERACT_SHELL_ {
  int index; /* the index of the shell within the SHELL_STATE */
  int n; /* principle quantum number */
  int j; /* the angular momentum of the shell, double of its actual value */
  int kl; /* orbital angular momentum */
  int kappa; /* relativistic angular quantum number kappa */
  int nq_bra; /* the occupation number in the bra state */
  int nq_ket; /* occupation number in the ket state */
} INTERACT_SHELL;

/*
** STRUCT:      INTERACT_DATUM
** PURPOSE:     the information about interacting shells to be 
**              saved for later use.
** FIELDS:      {SHELL *bra},
**              the shell structure of the bra state.
**              {short s[4]},
**              the indexes of all interacting shells,
**              s[0] and s[2] are shells from the bra state,
**              s[1] and s[3] are shells from the ket state.
**              {short n_shells}, 
**              number of shells in the bra state.
**              {short phase},
**              the phase resulting from the decoupling that depends
**              on the shell structure of the states.
** NOTE:        
*/
typedef struct _INTERACT_DATUM_ {
  SHELL *bra;
  INTERACT_SHELL s[4];
  short n_shells;
  short phase;
} INTERACT_DATUM;

#define MAXJ 80
#define MAXNJGD 2000
typedef int TRIADS[MAXJ][4];
typedef struct _FORMULA_ {
  int njgdata[MAXNJGD];
  int ifree[MAXJ], js[MAXJ], ns, ninter, j1, j2;
  int order[MAXJ], inter[MAXJ], interp[MAXJ], irank[MAXJ];
  TRIADS tr1, tr2;  
  double coeff;
  int phase;
} FORMULA;

#ifdef PERFORM_STATISTICS
/*
** STRUCT:      RECOUPLE_TIMING
** PURPOSE:     the timing statistics of the module "recouple".
** FIELDS:      {clock_t angz, angzxz, decouple, interact},
**              time used in the functions AngZ, AngZxZ0, 
**              DecoupleShell, and GetInteract respectively.
** NOTE:        it is only compiled in when the macro 
**              PERFORM_STATISTICS is defined in "global.h".
*/
typedef struct _RECOUPLE_TIMING_ {
  double angz;
  double angzxz;
  double decouple;
  double interact;
} RECOUPLE_TIMING;

int     GetRecoupleTiming(RECOUPLE_TIMING *t);
#endif

/* the recoupling matrix going from the coupled operators to uncoupled ones */
double  DecoupleShell(int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
		      int n_interact, int *interact, int *rank);
double  DecoupleShellRecursive(int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
			       int n_interact, int *interact, int *rank);
int IsPresent(int i, int n, int *m);

/* check if the recoupling matrix is non-zero */
int     IsShellInteracting(int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
			   int n_interact, int *interact, int *rank);

/* the coeff of type (Z^k dot Z^k) */
int     AngularZxZ0(double **coeff, int **kk, int nk, 
		    int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
		    INTERACT_SHELL *s);

/* the coeff of type Z^k */
int     AngularZ(double **coeff, int **kk, int nk,
		 int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
		 INTERACT_SHELL *s1, INTERACT_SHELL *s2);

void    SumCoeff(double *coeff,  int *kk,  int nk,  int p, 
		 double *coeff1, int *kk1, int nk1, int p1, 
		 int phase, int j1, int j2, int j3, int j4);

int     SortShell(int ns, INTERACT_SHELL *s, int *order);

/* analyze the structure of the configuration of bra and ket to 
   determine if they can interact, get the interacting shells that
   must interact, and determine the phase factor of the recoupling 
   coeff. which does not include the phase resulting from the reordering
   of operators. that is calculated in the AngularZxZ0 and AngularZ0 */
int InteractingShells(INTERACT_DATUM **idatum,
		      SHELL_STATE **sbra, 
		      SHELL_STATE **sket, 
		      CONFIG *ci, CONFIG *cj,
		      SHELL_STATE *csf_i, SHELL_STATE *csf_j);
int GetInteract(INTERACT_DATUM **idatum,
		SHELL_STATE **sbra, 
		SHELL_STATE **sket, 
		int kgi, int kgj,
		int kci, int kcj, 
		int ki, int kj, int bf);
void CompactInteractShell(char c[4], INTERACT_SHELL *s, int m);

/* only compile these test routines if the debug flag is on */
void    TestAngular(void);
void    CheckAngularConsistency(int n_shells, SHELL *bra, 
			     SHELL_STATE *sbra, SHELL_STATE *sket,
			     INTERACT_SHELL *s, int phase);
int     SetMaxRank(int k);
int     GetMaxRank(void);
int     InitRecouple(void);
int     ReinitRecouple(int m);
double EvaluateFormula(FORMULA *fm);
int GenerateFormula(FORMULA *fm);
int FixJsZ(INTERACT_SHELL *s, FORMULA *fm);
int TriadsZ(int n1, int n2, FORMULA *fm);
int CoupleSuccessive(int n, int *ik, int itr, TRIADS tr, int *i0);
int RecoupleTensor(int ns, INTERACT_SHELL *s, FORMULA *fm);
void EvaluateTensor(int nshells, SHELL_STATE *bra, SHELL_STATE *ket,
		    INTERACT_SHELL *s, int itr, FORMULA *fm);
#endif




