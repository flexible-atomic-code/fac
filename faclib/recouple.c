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

#include "recouple.h"
#include "structure.h"
#include "cf77.h"

static char *rcsid="$Id: recouple.c,v 1.28 2006/08/28 23:44:17 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

/*************************************************************
  Implementation of the module "recouple".
  This module calculates the recoupling coefficients. 

  The main task is to determine which electrons are the 
  interacting ones, and calculate the reduced matrix elements
  of the operator Z and ZxZ0, 

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

/*
** VARIABLE:    max_rank
** TYPE:        static int
** PURPOSE:     the maximum rank of the operators allowed.
** NOTE:        the ranks are represented by an integer that is
**              twice of its actuall value. The default value of 
**              16 means a maximum rank of 8.
*/
static int max_rank = MAXRANK;

/*
** VARIABLE:    interact_shells
** TYPE:        static MULTI *
** PURPOSE:     the multi-dimensional array stores the interacting 
**              shells information.
** NOTE:        
*/
static MULTI *interact_shells;

#ifdef PERFORM_STATISTICS
static RECOUPLE_TIMING timing = {0, 0, 0, 0};
/* 
** FUNCTION:    GetRecoupleTiming
** PURPOSE:     retrieve the timing information of the module "recouple".
** INPUT:       {RECOUPLE_TIMING *t},
**              pointer of the structure which will hold the result.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        
*/
int GetRecoupleTiming(RECOUPLE_TIMING *t) {
  memcpy(t, &timing, sizeof(timing));
  return 0;
}
#endif

static void InitInteractDatum(void *p, int n) {
  INTERACT_DATUM *d;
  int i;

  d = (INTERACT_DATUM *) p;
  for (i = 0; i < n; i++) {
    d[i].n_shells = 0;
    d[i].bra = NULL;
  }
}

/* 
** FUNCTION:    FreeInteractDatum
** PURPOSE:     free memory of an INTERACT_DATUM struct.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
static void FreeInteractDatum(void *p) {
  INTERACT_DATUM *d;
  
  if (!p) return;
  d = (INTERACT_DATUM *) p;
  if (d->n_shells > 0) {
    free(d->bra);
    d->bra = NULL;
    d->n_shells = -1;
  }
}
  
/* 
** FUNCTION:    SetMaxRank
** PURPOSE:     set the maximum rank of the operators.
** INPUT:       {int k},
**              the maximum rank. it should be twice the 
**              actual value.
** RETURN:      {int}
**              always 0.
** SIDE EFFECT: the static max_rank is set to k.
** NOTE:        
*/
int SetMaxRank(int k) {
  max_rank = k;
  return 0;
}

/* 
** FUNCTION:    GetMaxRank
** PURPOSE:     retrieve the maximum rank.
** INPUT:       
** RETURN:      {int},
**              the maximum rank.
** SIDE EFFECT: 
** NOTE:        
*/
int GetMaxRank(void) {
  return max_rank;
}

/* 
** MACRO:       IsOrder
** PURPOSE:     check the ordering of the shells.
** INPUT:       {int order[4]},
**              the shell indexes to be checked.
**              {int a, b, c, d},
**              the ordering required.
** RETURN:      {int},
**              0: if the order is not of the required type.
**              1: otherwise.
** SIDE EFFECT: 
** NOTE:        
*/
#define IsOrder(order, a, b, c, d) (((order)[0] == (a)) && \
				    ((order)[1] == (b)) && \
				    ((order)[2] == (c)) && \
				    ((order)[3] == (d)))

int IsPresent(int i, int n, int *m) {
  int k, j;

  k = 0;
  for (j = 0; j < n; j++) {
    if (i == m[j]) k++;
  }
  return k;
}

/* 
** FUNCTION:    DecoupleShell
** PURPOSE:     decouple the operators, so that their reduced 
**              matrix elements are expressed in terms of those 
**              involve individual shells.
** INPUT:       {int n_shells},
**              number of shells in the bra and ket states.
**              {SHELL_STATE *bra, *ket},
**              the shell states. if the bra and ket have different 
**              shell structures, they must be padded with empty shells
**              to make them identical. this is done in GetInteract.
**              {int n_interact},
**              number of interacting shells.
**              {int *interact},
**              the indexes of the interacting shells.
**              {int *rank},
**              the ranks of the operators.
** RETURN:      {double},
**              the decoupling coefficient.
** SIDE EFFECT: 
** NOTE:        
*/
double DecoupleShell(int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
		     int n_interact, int *interact, int *rank) {
  int i, j, k;
  double coeff;

  /* check the delta function for non-interacting shells first */
  i = 0;
  for (j = 0; j < n_interact; j++) {
    k = n_shells - interact[j] - 1;
    for (; i < k; i++) {
      if (bra[i].shellJ != ket[i].shellJ ||
	  bra[i].nu != ket[i].nu ||
	  bra[i].Nr != ket[i].Nr) return 0.0;
    }
    i++;
  }
  
  for (; i < n_shells; i++) {
    if (bra[i].shellJ != ket[i].shellJ ||
	bra[i].totalJ != ket[i].totalJ ||
	bra[i].nu != ket[i].nu ||
	bra[i].Nr != ket[i].Nr) return 0.0;
  }

  coeff = DecoupleShellRecursive(n_shells, bra, ket, n_interact, interact, rank);

  return coeff;
}
  
double DecoupleShellRecursive(int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
			      int n_interact, int *interact, int *rank) {
  double coeff, a;
  int Jbra, Jket, j1bra, j1ket, j2bra, j2ket, i;
  int k1, k2, k;

  coeff = 0.0;
  if (n_interact > 0) {
    if (n_shells >= n_interact && interact[0] < n_shells) {
      if (n_shells == 1) return 1.0; /* coeff. is 1, if one shell left */
      Jbra = (bra[0]).totalJ;
      Jket = (ket[0]).totalJ;
      j1bra = (bra[1]).totalJ;
      j1ket = (ket[1]).totalJ;
      j2bra = (bra[0]).shellJ;
      j2ket = (ket[0]).shellJ;      

      k = rank[0];
      if (interact[0] < n_shells-1) {
	k2 = 0;
	k1 = rank[0];
      } else {
	k2 = rank[1];
	k1 = (n_interact == 1) ? 0: rank[2];
      }
      
      /* it is a 9j symbol, although in most cases it reduces to a 
	 6j symbol or even a number. shall distinguish them if maximum
	 efficiency is needed here */
      coeff = (sqrt((Jbra+1.0)*(Jket+1.0)*(rank[0]+1.0)) *
	       W9j(j1bra, j2bra, Jbra, 
		   j1ket, j2ket, Jket,
		   k1, k2, k));
      if (fabs(coeff) < EPS30) return 0.0;
      if (interact[0] < n_shells - 1) {
	/* if the current shell is not an interacting one, the two shells
	   in bra and ket state should be identical. and the reduced matrix
	   element of the identity operator should be included */
	coeff *= sqrt(j2bra + 1);	
      } else {
	/* otherwise proceed to the next interacting one */
	n_interact--;
	interact++;
	rank += 2;
      }
      /* strip the outmost shell and call DecoupleShell recursively */
      n_shells--;
      bra++;
      ket++;
      a = DecoupleShellRecursive(n_shells, bra, ket, 
				 n_interact, interact, rank);
      coeff *= a;
    }
  } else {
    coeff = sqrt(bra[0].totalJ + 1.0);
  }
  return coeff;
}

/* 
** FUNCTION:    IsShellInteracting
** PURPOSE:     check if the states lead to non-zero matrix elements.
** INPUT:       {int n_shells},
**              number of shells in the bra and ket states.
**              {SHELL_STATE *bra, *ket},
**              the shell states. if the bra and ket have different 
**              shell structures, they must be padded with empty shells
**              to make them identical. this is done in GetInteract.
**              {int n_interact},
**              number of interacting shells.
**              {int *interact},
**              the indexes of the interacting shells.
**              {int *rank},
**              the ranks of the operators.
** RETURN:      {int},
**              0: the states alwas results in zero matrix elements.
**              1: otherwise.
** SIDE EFFECT: 
** NOTE:        
*/
int IsShellInteracting(int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
		       int n_interact, int *interact, int *rank) {
  int i, j, j0;

  j0 = 0;
  for (i = 0; i < n_interact; i++) {
    for (j = j0; j < n_shells - interact[i] -1; j++) {
      if ((bra[i]).totalJ != (ket[i]).totalJ ||
	  (bra[i]).shellJ != (ket[i]).shellJ ||
	  (bra[i]).nu != (ket[1]).nu)
	return 0;
    }
    j0 = n_shells - interact[i];
  }

  return 1;
}

/* 
** FUNCTION:    AngularZ
** PURPOSE:     calculate the reduced matrix element of Z operator.
** INPUT:       {double **coeff},
**              a pointer to a double array, which holds the 
**              results on exit.
**              {int **k},
**              a pointer to an int array, which holds all possible 
**              ranks of the operator.
**              {int n_shells},
**              number of the shells in the bra and ket states.
**              {SHELL_STATE *bra, *ket},
**              the bra and ket states.
**              {INTERACT_SHELL *s1, *s2},
**              two interacting shells involved in the operator.
** RETURN:      {int},
**              0: error occured.
**              1: succeeded.
** SIDE EFFECT: 
** NOTE:        For the definition of the operator see 
**              Bar-Shalom et al. Phys. Rev A. 38, 1773. 

**              if the rank of the operator is passed in, the total 
**              number of ranks and ranks should be set in nk, and **kk, 
**              and the storage for the coeff be provided. 
**              otherwise all possible ranks are determined and 
**              storage for the coeff. and kk allocated in this routine. 

**              Note that the reduced matrix element calculated here 
**              does not include the overall phase factor that arises 
**              from interchanging of electron shells, which depends 
**              on the occupation numbers of the config.
**              However, the phase factor arise from the interchanging 
**              of operators are included. This is also the same for 
**              the next routine, AngularZxZ0
*/   
int AngularZ(double **coeff, int **kk, int nk,
	     int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
	     INTERACT_SHELL *s1, INTERACT_SHELL *s2){
  SHELL_STATE st1, st2;
  int rank[4];
  int interact[2];
  int n_interact;
  int kmin, kmax, k, m;
  double coeff1, coeff2;
  RCFP_STATE rcfp_bra, rcfp_ket;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  if (nk == 0) {
    k = (bra[0]).totalJ;
    m = (ket[0]).totalJ;
    kmin = Max(abs(s1->j - s2->j), abs(k - m));
    kmax = Min(s1->j + s2->j, k + m);
    kmax = Min(kmax, max_rank);
    if (IsOdd(kmin)) kmin++; 
    if (kmax < kmin) return -1;
    nk = (kmax - kmin)/2 + 1;

    (*kk) = malloc(sizeof(int) * nk);

    if (!(*kk)) {
      return -1;
    }
    
    m = 0;
    for (k = kmin; k <= kmax; k += 2) {
      (*kk)[m++] = k;
    }

    (*coeff) = malloc(sizeof(double) * nk);
    if (!(*coeff)) return -1;  
  }

  if (s1->index == s2->index) { /* interacting shells are the same */
    n_interact = 1;
    interact[0] = s1->index;
    for (m = 0; m < nk; m++) {
      rank[0] = (*kk)[m];
      rank[1] = (*kk)[m];
      (*coeff)[m] = DecoupleShell(n_shells, bra, ket, 
				  n_interact, interact, rank);
#if (FAC_DEBUG >= DEBUG_RECOUPLE)
    fprintf(debug_log, "AngularZ: Decouple %d %lf\n", m, (*coeff)[m]);
#endif
     
    }
    st1 = bra[n_shells - s1->index -1];
    st2 = ket[n_shells - s1->index -1];
    rcfp_bra.state = RCFPTermIndex(s1->j, (st1).nu, 
				   (st1).Nr, (st1).shellJ);
    rcfp_bra.nq = s1->nq_bra;
    rcfp_bra.subshellMQ = rcfp_bra.nq - (s1->j + 1)/2;

    rcfp_ket.state = RCFPTermIndex(s1->j, (st2).nu, 
				   (st2).Nr, (st2).shellJ);
    rcfp_ket.nq = s1->nq_ket;
    rcfp_ket.subshellMQ = rcfp_ket.nq - (s1->j + 1)/2;
    for (k = 0; k < nk; k++) {
      (*coeff)[k] *= ReducedW(&rcfp_bra, &rcfp_ket, (*kk)[k]/2, 1, -1);
#if (FAC_DEBUG >= DEBUG_RECOUPLE)
    fprintf(debug_log, "AngularZ: %d %lf\n", k, (*coeff)[k]);
#endif
    }

  } else { /* two different interacting shells */
    if (s1->index < s2->index) {
      n_interact = 2;
      interact[0] = s2->index;
      interact[1] = s1->index;
      rank[1] = s2->j;
      rank[2] = s1->j;
      rank[3] = s1->j;
    } else {
      n_interact = 2;
      interact[0] = s1->index;
      interact[1] = s2->index;
      rank[1] = s1->j;
      rank[2] = s2->j;
      rank[3] = s2->j;
    }

    for (m = 0; m < nk; m++) {
      rank[0] = (*kk)[m];
      (*coeff)[m] = DecoupleShell(n_shells, bra, ket, 
				  n_interact, interact, rank);  
    }
 
    st1 = bra[n_shells - s1->index -1];
    st2 = ket[n_shells - s1->index -1];
    rcfp_bra.state = RCFPTermIndex(s1->j, (st1).nu, 
				   (st1).Nr, (st1).shellJ);
    rcfp_bra.nq = s1->nq_bra;
    rcfp_bra.subshellMQ = rcfp_bra.nq - (s1->j + 1)/2;

    rcfp_ket.state = RCFPTermIndex(s1->j, (st2).nu, 
				   (st2).Nr, (st2).shellJ);
    rcfp_ket.nq = s1->nq_ket;
    rcfp_ket.subshellMQ = rcfp_ket.nq - (s1->j + 1)/2;  
    coeff1 = ReducedA(&rcfp_bra, &rcfp_ket, 1);

    st1 = bra[n_shells - s2->index -1];
    st2 = ket[n_shells - s2->index -1];
    rcfp_bra.state = RCFPTermIndex(s2->j, (st1).nu, 
				   (st1).Nr, (st1).shellJ);
    rcfp_bra.nq = s2->nq_bra;
    rcfp_bra.subshellMQ = rcfp_bra.nq - (s2->j + 1)/2;
    rcfp_ket.state = RCFPTermIndex(s2->j, (st2).nu, 
				   (st2).Nr, (st2).shellJ);
    rcfp_ket.nq = s2->nq_ket;
    rcfp_ket.subshellMQ = rcfp_ket.nq - (s2->j + 1)/2;  
    coeff2 = ReducedA(&rcfp_bra, &rcfp_ket, -1);
    for (m = 0; m < nk; m++) {
      if (s1->index > s2->index && IsEven((s1->j + s2->j - (*kk)[m])/2)) {
	(*coeff)[m] = -(*coeff)[m];
      }
      (*coeff)[m] *= coeff1*coeff2;
    }
  }
 
  for (k = 0; k < nk; k++) {
    (*coeff)[k] /= -sqrt((*kk)[k] + 1);
#if FAC_DEBUG    
    fprintf(debug_log, "AngularZ: %d %lf\n", (*kk)[k], (*coeff)[k]);
#endif
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.angz += stop-start;
#endif
  return nk;
}

/* 
** FUNCTION:    AngularZxZ0
** PURPOSE:     calculate the reduced matrix element of 
**              (Z \dot Z) operator.
** INPUT:       {double **coeff},
**              a pointer to a double array, which holds the 
**              results on exit.
**              {int **k},
**              a pointer to an int array, which holds all possible 
**              ranks of the operator.
**              {int n_shells},
**              number of the shells in the bra and ket states.
**              {SHELL_STATE *bra, *ket},
**              the bra and ket states.
**              {INTERACT_SHELL s[4]},
**              4 interacting shells involved in the operator.
** RETURN:      {int},
**              0: error occured.
**              1: succeeded.
** SIDE EFFECT: 
** NOTE:        see the notes of the routine AngularZ. 
**              Note that the order of the interacting shells is 
**              different from that in the radial part. 
**              e.g. if the order in the slater integral is R(ab, cd), 
**              then the order passing into this routine is a, c, b, d 
*/
int AngularZxZ0(double **coeff, int **kk, int nk,
		int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
		INTERACT_SHELL *s) {
  
  SHELL_STATE st1, st2;
  int rank[8];
  int interact[4];
  int nops[4]={0,0,0,0};
  int n_interact;
  int order[4] = {0, 1, 2, 3};
  int kmin, kmax, k, m;
  int qm[4] = {1, -1, 1, -1};
  int recouple_operators = 0;
  int phase;
  double *coeff1;
  int *kk1, nk1;
  int i, j;
  double r, a, b;  
  RCFP_STATE rcfp_bra[4], rcfp_ket[4];

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  /* if nk is non positive, allocate the memory for the coeff. and kk */
  if (nk <= 0) {
    kmin = Max(abs(s[0].j-s[1].j), abs(s[2].j-s[3].j));
    kmax = Min(s[0].j+s[1].j, s[2].j+s[3].j);
    kmax = Min(kmax, max_rank);
    if (IsOdd(kmin)) kmin++; 
    if (kmax < kmin) return -1;
    nk = (kmax - kmin)/2 + 1;
    (*kk) = malloc(sizeof(int)*nk);
    if (!(*kk)) {
      return -1;
    }
    m = 0;
    for (k = kmin; k <= kmax; k += 2) {
      (*kk)[m++] = k;
    }

    (*coeff) = malloc(sizeof(double)*nk);
    if (!(*coeff)) return -1;
  }

  coeff1 = *coeff;
  kk1 = *kk;
  nk1 = nk;

  /* sort the order of interacting shells */
  phase = SortShell(4, s, order);
  n_interact = 0;
  j = -1;
  for (i = 3; i >= 0; i--) {
    if (i == 3 || s[order[i]].index != s[order[i+1]].index) {
      n_interact++;
      j++;
      interact[j] = s[order[i]].index;
      nops[j] = 1;
      st1 = bra[n_shells - interact[j] - 1];
      st2 = ket[n_shells - interact[j] - 1];
      rcfp_bra[j].state = RCFPTermIndex(s[order[i]].j, (st1).nu,
					(st1).Nr, (st1).shellJ);
      rcfp_bra[j].nq = s[order[i]].nq_bra;
      rcfp_bra[j].subshellMQ = rcfp_bra[j].nq - (s[order[i]].j + 1)/2;
      rcfp_ket[j].state = RCFPTermIndex(s[order[i]].j, (st2).nu,
					(st2).Nr, (st2).shellJ);
      rcfp_ket[j].nq = s[order[i]].nq_ket;
      rcfp_ket[j].subshellMQ = rcfp_ket[j].nq - (s[order[i]].j + 1)/2;
    } else {
      nops[j]++;
    }
  }
  

  rank[0] = 0;
  
  if (nops[0] == 4) { /* 4 identical interacting shells */
    for (m = 0; m < nk1; m++){
      rank[1] = 0;      
      a = DecoupleShell(n_shells, bra, ket, 
			n_interact, interact, rank); 
      b = ReducedWxW0(rcfp_bra, rcfp_ket, kk1[m]/2, 1, -1, 1, -1);
      coeff1[m] = a*b;
    }
    recouple_operators = 0;
  } else if (nops[0] == 3 && nops[1] == 1) { /* 1 x 3 */
#if FAC_DEBUG    
    fprintf(debug_log, "%d %d  %d %d  %d %d  %d %d\n",
	    s[0].nq_bra, s[0].nq_ket, s[1].nq_bra, s[1].nq_ket, 
	    s[2].nq_bra, s[2].nq_ket, s[3].nq_bra, s[3].nq_ket);
#endif
    for (m = 0; m < nk1; m++){
      rank[1] = s[order[0]].j;
      rank[2] = rank[1];
      rank[3] = rank[1];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 
#if FAC_DEBUG
      fprintf(debug_log, "  1x3 rank %d %lf\n", kk1[m], coeff1[m]);
#endif
      if (order[0] == 0 || order[0] == 1) { 
	coeff1[m] *= ReducedAxW(rcfp_bra, rcfp_ket, kk1[m]/2, rank[1],
				qm[order[1]], qm[order[2]], qm[order[3]]);
#if FAC_DEBUG
	fprintf(debug_log, "1x3:1 rank %d %lf\n", kk1[m], coeff1[m]);
#endif
	coeff1[m] *= ReducedA(rcfp_bra+1, rcfp_ket+1, qm[order[0]]);
#if FAC_DEBUG
	fprintf(debug_log, "1x3:1 rank %d %lf\n", kk1[m], coeff1[m]);
#endif
      } else {     
	coeff1[m] *= ReducedWxA(rcfp_bra, rcfp_ket, kk1[m]/2, rank[1],
				qm[order[1]], qm[order[2]], qm[order[3]]);
#if FAC_DEBUG
	fprintf(debug_log, "1x3:2 rank %d %lf\n", kk1[m], coeff1[m]);
#endif
	coeff1[m] *= ReducedA(rcfp_bra+1, rcfp_ket+1, qm[order[0]]);
#if FAC_DEBUG
	fprintf(debug_log, "1x3:2 rank %d %lf\n", kk1[m], coeff1[m]);
#endif
      }
    }
    if (order[0] == 0 || order[0] == 1) recouple_operators = 1;
    else recouple_operators = 2;
  } else if (nops[0] == 1 && nops[1] == 3) { /* 3 x 1 */
    for (m = 0; m < nk1; m++){
      rank[1] = s[order[3]].j;
      rank[2] = rank[1];
      rank[3] = rank[1];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 
#if FAC_DEBUG
      fprintf(debug_log, "  3x1 rank %d %lf\n", kk1[m], coeff1[m]);
#endif
      if (order[3] == 0 || order[3] == 1) {
	coeff1[m] *= ReducedA(rcfp_bra, rcfp_ket, qm[order[3]]);
#if FAC_DEBUG
	fprintf(debug_log, "3x1:3 rank %d %lf\n", kk1[m], coeff1[m]);
#endif
	coeff1[m] *= ReducedAxW(rcfp_bra+1, rcfp_ket+1, kk1[m]/2, rank[1],
				qm[order[0]], qm[order[1]], qm[order[2]]);
#if FAC_DEBUG
	fprintf(debug_log, "3x1:3 rank %d %lf\n", kk1[m], coeff1[m]);
#endif
      } else {
	coeff1[m] *= ReducedA(rcfp_bra, rcfp_ket, qm[order[3]]);
#if FAC_DEBUG
	fprintf(debug_log, "3x1:4 rank %d %lf\n", kk1[m], coeff1[m]);
#endif
	coeff1[m] *= ReducedWxA(rcfp_bra+1, rcfp_ket+1, kk1[m]/2, rank[1],
				qm[order[0]], qm[order[1]], qm[order[2]]);
#if FAC_DEBUG
	fprintf(debug_log, "3x1:4 rank %d %lf\n", kk1[m], coeff1[m]);
#endif
      }
    }
    if (order[3] == 0 || order[3] == 1) recouple_operators = 3;
    else recouple_operators = 4;  
  } else if (nops[0] == 1 && nops[1] == 2) { /* 1 x 2 x 1 */
    if (!(order[1] == 0 && order[2] == 1) &&
	!(order[1] == 2 && order[2] == 3)) {
      m = 0;
      kmin = abs(s[order[0]].j - s[order[3]].j);
      kmax = Min(2*s[order[1]].j, s[order[0]].j+s[order[3]].j);
      if (IsOdd(kmin)) kmin++; 
      if (kmax < kmin) return -1;
      nk1 = (kmax - kmin)/2 + 1;
      coeff1 = malloc(sizeof(double)*nk1);
      kk1 = malloc(sizeof(int)*nk1);
      if (!kk1 || !coeff1) return -1;
      for (k = kmin; k <= kmax; k += 2) {
	kk1[m++] = k;
      }
    }      
    for (m = 0; m < nk1; m++){
      rank[1] = s[order[3]].j;
      rank[2] = s[order[3]].j;
      rank[3] = kk1[m];
      rank[4] = s[order[0]].j;
      rank[5] = rank[4];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 
      coeff1[m] *= ReducedA(rcfp_bra, rcfp_ket, qm[order[3]]);
      coeff1[m] *= ReducedW(rcfp_bra+1, rcfp_ket+1, kk1[m]/2,
			    qm[order[1]], qm[order[2]]);
      coeff1[m] *= ReducedA(rcfp_bra+2, rcfp_ket+2, qm[order[0]]);
    }
    recouple_operators = 5;
  } else if (nops[0] == 2 && nops[1] == 2) { /* 2 x 2 */
    if (!(order[0] == 0 && order[1] == 1) &&
	!(order[0] == 2 && order[1] == 3)) {
      m = 0;
      kmin = 0;
      kmax = 2 * Min(s[order[0]].j, s[order[2]].j);
      if (IsOdd(kmin)) kmin++; 
      if (kmax < kmin) return -1;
      nk1 = (kmax - kmin)/2 + 1;
      coeff1 = malloc(sizeof(double)*nk1);
      kk1 = malloc(sizeof(int)*nk1);
      if (!kk1 || !coeff1) return -1;
      for (k = kmin; k <= kmax; k += 2) {
	kk1[m++] = k;
      }
    }

#if (FAC_DEBUG >= DEBUG_RECOUPLE)
    fprintf(debug_log, "AngularZxZ0: Order, %d %d %d %d, nk1 = %d\n", 
	    order[0], order[1], order[2], order[3], nk1);
#endif

    for (m = 0; m < nk1; m++){
      rank[1] = kk1[m];
      rank[2] = rank[1];
      rank[3] = rank[1];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 

#if (FAC_DEBUG >= DEBUG_RECOUPLE)
      fprintf(debug_log, "2x2 rank: %d, %lf \n", kk1[m], coeff1[m]);
#endif

      r = ReducedW(rcfp_bra, rcfp_ket, kk1[m]/2, 
		   qm[order[2]], qm[order[3]]);
      coeff1[m] *= r;

#if (FAC_DEBUG >= DEBUG_RECOUPLE)
      fprintf(debug_log, "2x2 rank: %d, %lf \n", kk1[m], r);
      fprintf(debug_log, "%d %d %d %d %d %d %d %d \n", 
	      rcfp_bra[0].state, rcfp_bra[0].nq, rcfp_bra[0].subshellMQ,
	      rcfp_ket[0].state, rcfp_ket[0].nq, rcfp_ket[0].subshellMQ,
	      qm[order[2]], qm[order[3]]);
#endif

      r = ReducedW(rcfp_bra+1, rcfp_ket+1, kk1[m]/2,
		   qm[order[0]], qm[order[1]]);
      coeff1[m] *= r;

#if (FAC_DEBUG >= DEBUG_RECOUPLE)
      fprintf(debug_log, "2x2 rank: %d, %lf \n", kk1[m], r);
      fprintf(debug_log, "%d %d %d %d %d %d %d %d \n", 
	      rcfp_bra[1].state, rcfp_bra[1].nq, rcfp_bra[1].subshellMQ,
	      rcfp_ket[1].state, rcfp_ket[1].nq, rcfp_ket[1].subshellMQ,
	      qm[order[0]], qm[order[1]]);
#endif

    }
    recouple_operators = 6;
  } else if (nops[0] ==  2 && nops[1] == 1) { /* 1 x 1 x 2 */
    if (!(order[2] == 0 && order[3] == 1) &&
	!(order[2] == 2 && order[3] == 3)) {
      m = 0;
      kmin = abs(s[order[0]].j - s[order[1]].j);
      kmax = Min(2*s[order[3]].j, s[order[0]].j+s[order[1]].j);
      if (IsOdd(kmin)) kmin++; 
      if (kmax < kmin) return -1;
      nk1 = (kmax - kmin)/2 + 1;
      coeff1 = malloc(sizeof(double)*nk1);
      kk1 = malloc(sizeof(int)*nk1);
      if (!kk1 || !coeff1) return -1;
      for (k = kmin; k <= kmax; k += 2) {
	kk1[m++] = k;
      }
    }
    for (m = 0; m < nk1; m++) {
      rank[1] = kk1[m];
      rank[2] = kk1[m];
      rank[3] = s[order[1]].j;
      rank[4] = s[order[0]].j;
      rank[5] = rank[4];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 
      coeff1[m] *= ReducedW(rcfp_bra, rcfp_ket, kk1[m]/2,
			    qm[order[2]], qm[order[3]]);
      coeff1[m] *= ReducedA(rcfp_bra+1, rcfp_ket+1, qm[order[1]]);
      coeff1[m] *= ReducedA(rcfp_bra+2, rcfp_ket+2, qm[order[0]]);
    }
    recouple_operators = 7;
  } else if (nops[0] == 1 && nops[1] == 1 && nops[2] == 2) { /* 2 x 1 x 1 */
    if (!(order[0] == 0 && order[1] == 1) &&
	!(order[0] == 2 && order[1] == 3)) {
      m = 0;
      kmin = abs(s[order[2]].j - s[order[3]].j);
      kmax = Min(2*s[order[0]].j, s[order[2]].j+s[order[3]].j);
      if (IsOdd(kmin)) kmin++; 
      if (kmax < kmin) return -1;
      nk1 = (kmax - kmin)/2 + 1;
      coeff1 = malloc(sizeof(double)*nk1);
      kk1 = malloc(sizeof(int)*nk1);
      if (!kk1 || !coeff1) return -1;
      for (k = kmin; k <= kmax; k += 2) {
	kk1[m++] = k;
      }
    }
    for (m = 0; m < nk1; m++) {
      rank[1] = s[order[3]].j;
      rank[2] = s[order[3]].j;
      rank[3] = s[order[2]].j;
      rank[4] = kk1[m];
      rank[5] = rank[4];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 

#if (FAC_DEBUG > DEBUG_RECOUPLE)
      fprintf(debug_log, "2x1x1 rank, %d %lf \n", kk1[m], coeff1[m]);
#endif

      r = ReducedA(rcfp_bra, rcfp_ket, qm[order[3]]);
      coeff1[m] *= r;

#if (FAC_DEBUG > DEBUG_RECOUPLE)
      fprintf(debug_log, "2x1x1 rank, %d %lf \n", kk1[m], r);
#endif

      r = ReducedA(rcfp_bra+1, rcfp_ket+1, qm[order[2]]);
      coeff1[m] *= r;

#if (FAC_DEBUG > DEBUG_RECOUPLE)
      fprintf(debug_log, "2x1x1 rank, %d %lf \n", kk1[m], r);
#endif

      r = ReducedW(rcfp_bra+2, rcfp_ket+2, kk1[m]/2,
		   qm[order[0]], qm[order[1]]);
      coeff1[m] *= r;

#if (FAC_DEBUG > DEBUG_RECOUPLE)
      fprintf(debug_log, "2x1x1 rank, %d %lf \n", kk1[m], r);
#endif

    }
    recouple_operators = 8;
  } else { /* 1 x 1 x 1 x 1 */
    i = order[0] + order[1];
    if (i != 1 && i != 5) {
      m = 0;
      kmin = Max(abs(s[order[0]].j - s[order[1]].j), 
		 abs(s[order[2]].j - s[order[3]].j));
      kmax = Min(s[order[0]].j + s[order[1]].j,
		 s[order[2]].j + s[order[3]].j);
      if (IsOdd(kmin)) kmin++; 
      if (kmax < kmin) return -1;
      nk1 = (kmax - kmin)/2 + 1;
      coeff1 = malloc(sizeof(double)*nk1);
      kk1 = malloc(sizeof(int)*nk1);
      if (!kk1 || !coeff1) return -1;
      for (k = kmin; k <= kmax; k += 2) {
	kk1[m++] = k;
      }
    }
    for (m = 0; m < nk1; m++) {
      rank[1] = s[order[3]].j;
      rank[2] = rank[1];
      rank[3] = s[order[2]].j;
      rank[4] = kk1[m];
      rank[5] = s[order[1]].j;
      rank[6] = s[order[0]].j;
      rank[7] = rank[6];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank);
      coeff1[m] *= ReducedA(rcfp_bra, rcfp_ket, qm[order[3]]);
      coeff1[m] *= ReducedA(rcfp_bra+1, rcfp_ket+1, qm[order[2]]);
      coeff1[m] *= ReducedA(rcfp_bra+2, rcfp_ket+2, qm[order[1]]);
      coeff1[m] *= ReducedA(rcfp_bra+3, rcfp_ket+3, qm[order[0]]);
    }
    recouple_operators = 9;
  }

#if FAC_DEBUG
  fprintf(debug_log, "shell order: %d %d %d %d %d %d\n",
	  order[0], order[1], order[2], order[3], phase, recouple_operators);
#endif
  /* now adjust the phase factor, and do the recoupling if necessary */
  switch (recouple_operators) {
  case 0:
    break;
  case 1:
    if (IsOrder(order, 0, 1, 2, 3)) {
      if (IsOdd(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 1, 0, 2, 3)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[0].j+s[1].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    }
    break;
  case 2:
    if (IsOrder(order, 2, 0, 1, 3)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[3].j-s[2].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 3, 0, 1, 2)) {
      if (IsEven(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    }
    break;
  case 3:
    if (IsOrder(order, 1, 2, 3, 0)) {
      if (IsEven(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 0, 2, 3, 1)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[0].j-s[1].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    }
    break;
  case 4:
    if (IsOrder(order, 0, 1, 2, 3)) {
      if (IsOdd(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 0, 1, 3, 2)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[2].j+s[3].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    }
    break;
  case 5:
    if (IsOrder(order, 2, 0, 1, 3)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[2].j-s[3].j+kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 3, 0, 1, 2)) {     
      if (IsEven(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 0, 2, 3, 1)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[0].j-s[1].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 1, 2, 3, 0)) { 
      if (IsEven(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    } else {
      if (IsOrder(order, 1, 0, 2, 3)) {
	phase += (s[2].j + s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 3, 0, 2, 1)) {
	phase += (s[1].j-s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 1, 0, 3, 2)) {
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 2, 0, 3, 1)) {
	phase += (s[1].j+s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 0, 1, 2, 3)) {
	phase += (s[0].j+s[1].j+s[2].j+s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 3, 1, 2, 0)) {
	phase += (s[1].j+s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 0, 1, 3, 2)) {
	phase += (s[0].j+s[1].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 2, 1, 3, 0)) {
	phase += (s[1].j - s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      }
      free(coeff1);
      free(kk1);
    }
    break;
  case 6:
    if (IsOrder(order, 0, 1, 2, 3) ||
	IsOrder(order, 2, 3, 0, 1)) {
      if (IsOdd(phase)) {
	for (m = 0; m < nk1; m++) coeff1[m] = -coeff1[m];
      }
    } else {
      if (IsOrder(order, 0, 2, 1, 3) ||
	  IsOrder(order, 1, 3, 0, 2)) {
	phase += (s[1].j+s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 3, 1, 2) ||
		 IsOrder(order, 1, 2, 0, 3)) {
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);

#if (FAC_DEBUG >= DEBUG_RECOUPLE)
	fprintf(debug_log, "AngularZxZ0: recouple_operators = 6, %lf, %lf\n",
		(*coeff)[0], coeff1[0]);
#endif

      }
      free(coeff1);
      free(kk1);
    }
    break;      
  case 7:
    if (IsOrder(order, 0, 1, 2, 3) ||
	IsOrder(order, 2, 3, 0, 1)) {
      /* do nothing in this case */
    } else if (IsOrder(order, 1, 0, 2, 3)) {
      for (m = 0; m < nk1; m++) {
#if FAC_DEBUG
	fprintf(debug_log, "rank*: %d, phase: %d %d %lf \n", kk1[m], phase,
		(s[0].j + s[1].j - kk1[m])/2, coeff1[m]);
#endif
	if (IsOdd(phase + (s[0].j+s[1].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 3, 2, 0, 1)) {
      for (m = 0; m < nk1; m++) {
	if (IsOdd(phase + (s[2].j+s[3].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
      }
    } else {
      if (IsOrder(order, 0, 2, 1, 3) ||
	  IsOrder(order, 1, 3, 0, 2)) {
	phase += (s[1].j+s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 2, 0, 1, 3)) {
	phase += (s[0].j-s[1].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 3, 1, 0, 2)) {
	phase += (s[2].j - s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 3, 1, 2)) {
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 3, 0, 1, 2)) {
	phase += (s[0].j+s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase, 
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 1, 2, 0, 3)) {
	phase += (s[1].j - s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 2, 1, 0, 3)) {
	phase += 1;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      }
      free(coeff1);
      free(kk1);
    }
    break; 
  case 8:
    if (IsOrder(order, 0, 1, 2, 3) ||
	IsOrder(order, 2, 3, 0, 1)) {
      if (IsOdd(phase)) {
	for (m = 0; m < nk1; m++) coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 0, 1, 3, 2)) {

#if (FAC_DEBUG > DEBUG_RECOUPLE)
      fprintf(debug_log, "phase: %d %lf\n", phase, coeff1[0]);
#endif

      for (m = 0; m < nk1; m++) {
	if (IsOdd(phase + (s[2].j+s[3].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 2, 3, 1, 0)) {

#if (FAC_DEBUG > DEBUG_RECOUPLE)
      fprintf(debug_log, "phase: %d %lf\n", phase, coeff1[0]);
#endif

      for (m = 0; m < nk1; m++) {
	if (IsOdd(phase + (s[0].j+s[1].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
      }
    } else {
      if (IsOrder(order, 0, 2, 1, 3) ||
	  IsOrder(order, 1, 3, 0, 2)) {
	phase += (s[1].j+s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 2, 3, 1)) {
	phase += (s[2].j-s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 1, 3, 2, 0)) {
	phase += (s[0].j - s[1].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 3, 1, 2)) {
	phase += (s[1].j - s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 0, 3, 2, 1)) {
	phase += 1;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase, 
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 1, 2, 0, 3)) {
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 1, 2, 3, 0)) {
	phase += (s[0].j+s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      }
      free(coeff1);
      free(kk1);
    }
    break; 
  case 9:
    if (IsOrder(order, 0, 1, 2, 3) ||
	IsOrder(order, 2, 3, 0, 1)) {
      if (IsOdd(phase)) 
	for (m = 0; m < nk1; m++) coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 0, 1, 3, 2) ||
	       IsOrder(order, 3, 2, 0, 1)) {
      for (m = 0; m < nk1; m++) 
	if (IsOdd(phase + (s[2].j + s[3].j - kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 1, 0, 2, 3) ||
	       IsOrder(order, 2, 3, 1, 0)) {
      for (m = 0; m < nk1; m++) 
	if (IsOdd(phase + (s[0].j + s[1].j - kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 1, 0, 3, 2) ||
	       IsOrder(order, 3, 2, 1, 0)) {
      for (m = 0; m < nk1; m++) 
	if (IsOdd(phase+(s[0].j+s[1].j+s[2].j+s[3].j)/2)) 
	  coeff1[m] = -coeff1[m];
    } else {
      if (IsOrder(order, 0, 2, 1, 3) ||
	  IsOrder(order, 1, 3, 0, 2)) {
	phase += (s[1].j + s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 2, 3, 1) ||
		 IsOrder(order, 3, 1, 0, 2)) {
	phase += (s[2].j - s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 2, 0, 1, 3) ||
		 IsOrder(order, 1, 3, 2, 0)) {
	phase += (s[0].j - s[1].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 2, 0, 3, 1) ||
		 IsOrder(order, 3, 1, 2, 0)) {
	phase += (s[0].j + s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 3, 1, 2) ||
		 IsOrder(order, 1, 2, 0, 3)) {
	phase += (s[1].j - s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 0, 3, 2, 1) ||
		 IsOrder(order, 2, 1, 0, 3)) {
	phase += 1;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 3, 0, 1, 2) ||
		 IsOrder(order, 1, 2, 3, 0)) {
	phase += (s[0].j + s[3].j + s[1].j - s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 3, 0, 2, 1) ||
		 IsOrder(order, 2, 1, 3, 0)) {
	phase += (s[0].j - s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      }      
      free(coeff1);
      free(kk1);
    }
    break;
  }

  /* adjust the prefactor */
  for (m = 0; m < nk; m++) {
    /* this factor arise from the definition of Z^k and the 
       scalar product. */
    (*coeff)[m] /= sqrt((*kk)[m] + 1);
    if (IsOdd((*kk)[m]/2)) (*coeff)[m] = -(*coeff)[m];
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.angzxz += stop - start;
#endif

  return nk;
}

/* 
** FUNCTION:    SumCoeff
** PURPOSE:     perform the summation due to the exchange of 
**              two operators.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
void SumCoeff(double *coeff,  int *kk,  int nk,  int p, 
	      double *coeff1, int *kk1, int nk1, int p1, 
	      int phase, int j1, int j2, int j3, int j4) {
  int i, j;
  double x;
  
  for (i = 0; i < nk; i++) {
    coeff[i] = 0.0;
    for (j = 0; j < nk1; j++) {
      if (fabs(coeff1[j]) > 0.0) {
	x = W6j(j1, j2, kk[i], j3, j4, kk1[j]);
	if (fabs(x) > 0.0) {
	  x *= sqrt(kk1[j] + 1.0);
	  if (p1 && IsOdd(kk1[j]/2)) x = -x;
	  coeff[i] += x * coeff1[j];
	}
      }
    }
    if (fabs(coeff[i]) > 0.0) {
      coeff[i] *= sqrt(kk[i]+1.0);
      if (p) {
	if (IsOdd(phase + kk[i]/2)) coeff[i] = -coeff[i];
      } else {
	  if (IsOdd(phase)) coeff[i] = -coeff[i];
      }
    }
  }
}

/* 
** FUNCTION:    SortShell
** PURPOSE:     sort the interacting shells. 
**              calculate the phase of the interchange.
** INPUT:       {int ns},
**              number of shells.
**              {INTERACT_SHELL *s},
**              shell indexes to be sorted.
**              {int *order},
**              the order returned.
** RETURN:      {int},
**              the phase.
** SIDE EFFECT: 
** NOTE:        
*/
int SortShell(int ns, INTERACT_SHELL *s, int *order) {
  int i, j, k;
  int phase;

  phase = 0;
  for (j = 0; j < ns-1; j++){
    for (i = ns-1; i > j; i--) {
      if (s[order[i]].index < s[order[i-1]].index) {
	k = order[i];
	order[i] = order[i-1];
	order[i-1] = k;
	phase++;
      }
    }
  }
  return phase;
}

int InteractingShells(INTERACT_DATUM **idatum,
		      SHELL_STATE **sbra, 
		      SHELL_STATE **sket, 
		      CONFIG *ci, CONFIG *cj,
		      SHELL_STATE *csf_i, SHELL_STATE *csf_j) {
  int i, j, k, m, pb, pk;
  int n, kl, jj, nq, nq_plus, nq_minus, qd;
  SHELL *bra;
  INTERACT_SHELL *s;
  int interaction[4] = {-1, -1, -1, -1};
  int n_shells;

  s = (*idatum)->s;

  for (i = 0; i < 4; i++) s[i].index = -1;
  /* allocate memory. use the sum of two configs to avoid couting the
     exact size in advance */
  k = ci->n_shells + cj->n_shells;
  (*idatum)->bra = malloc(sizeof(SHELL)*k);
  if ((*idatum)->bra == NULL) {
    printf("error allocating idatam->bra, likely too many configurations.\n");
    exit(1);
  }
  bra = (*idatum)->bra;
  s = (*idatum)->s;
  i = 0; 
  j = 0;
  m = 0;
  pb = 0;
  pk = 1;
  nq_plus = 0;
  nq_minus = 0;
  while (1) {
    if (i >= ci->n_shells) {
      if (j < cj->n_shells) {
	k = -1;
      } else {
	break;      
      }
    } else {
      if (j >= cj->n_shells) {
	k = 1;
      } else {
	k = CompareShell(ci->shells+i, cj->shells+j);
      }
    }
    if (k > 0) { /* bra has a shell that does not exist in ket */
      if (nq_plus >= 2) break;
      if (ci->shells[i].nq > 0) {
	memcpy(bra+m, ci->shells+i, sizeof(SHELL));
	UnpackShell(bra+m, &n, &kl, &jj, &nq);
	nq_plus += nq;
	if (nq_plus > 2) break;
	s[pb].index = m;
	s[pb].n = n;
	s[pb].kappa = bra[m].kappa;
	s[pb].j = jj;
	s[pb].kl = kl;
	s[pb].nq_bra = nq;
	s[pb].nq_ket = 0;
	pb += 2;
	if (nq == 2) {
	  memcpy(s+pb, s, sizeof(INTERACT_SHELL));
	} else {
	  interaction[pb-2] = m;
	}
	m++;
      }
      i++;
    } else if (k < 0) { /* ket has a shell that does not exist in bra */
      if (nq_minus >= 2) break;
      if (cj->shells[j].nq > 0) {
	memcpy(bra+m, cj->shells+j, sizeof(SHELL));
	UnpackShell(bra+m, &n, &kl, &jj, &nq);
	nq_minus += nq;
	if (nq_minus > 2) break;
	s[pk].index = m;
	s[pk].n = n;
	s[pk].kappa = bra[m].kappa;
	s[pk].j = jj;
	s[pk].kl = kl;
	s[pk].nq_bra = 0;
	bra[m].nq = 0;
	s[pk].nq_ket = nq;
	pk += 2;
	if (nq == 2) {
	  memcpy(s+pk, s+1, sizeof(INTERACT_SHELL));
	} else {
	  interaction[pk-2] = m;
	}
	m++;
      }
      j++;
    } else { /* both bra and ket has the shell */
      memcpy(bra+m, ci->shells+i, sizeof(SHELL));
      qd = ci->shells[i].nq - cj->shells[j].nq;
      if (qd > 0) { /* bra has more electrons in the shell */
	if (nq_plus >= 2) break;
	UnpackShell(bra+m, &n, &kl, &jj, &nq);
	nq_plus += qd;
	if (nq_plus > 2) break;
	s[pb].index = m;
	s[pb].n = n;
	s[pb].kappa = bra[m].kappa;
	s[pb].j = jj;
	s[pb].kl = kl;
	s[pb].nq_bra = nq;
	s[pb].nq_ket = cj->shells[j].nq;
	pb += 2;
	if (qd == 2) {
	  memcpy(s+pb, s, sizeof(INTERACT_SHELL));
	} else {
	  interaction[pb-2] = m;
	}
      } else if (qd < 0) { /* ket has more electrons in the shell */
	if (nq_minus >= 2) break;
	UnpackShell(bra+m, &n, &kl, &jj, &nq);
	nq_minus += -qd;
	if (nq_minus > 2) break;
	s[pk].index = m;
	s[pk].n = n;
	s[pk].kappa = bra[m].kappa;
	s[pk].j = jj;
	s[pk].kl = kl;
	s[pk].nq_bra = nq;
	s[pk].nq_ket = cj->shells[j].nq;
	pk += 2;
	if (qd == -2) {
	  memcpy(s+pk, s+1, sizeof(INTERACT_SHELL));
	} else {
	  interaction[pk-2] = m;
	}
      }
      i++;
      j++;
      m++;
    }
  }
  n_shells = m;
  if (nq_plus != nq_minus ||
      nq_plus > 2 ||
      nq_minus > 2 ||
      i < ci->n_shells ||
      j < cj->n_shells) {
    free(bra);
    bra = NULL;
    n_shells = -1;
    goto END;
  }    

  if (csf_i == NULL) goto END;

#if (FAC_DEBUG > DEBUG_RECOUPLE)
  fprintf(debug_log, "before sort: %d %d %d %d\n", 
	  interaction[0], interaction[1], 
	  interaction[2], interaction[3]);
#endif

  /* determine the phase factor */
  for (j = 0; j < 3; j++){
    for (i = 3; i > j; i--) {
      if (interaction[i] < interaction[i-1]) {
	k = interaction[i];
	interaction[i] = interaction[i-1];
	interaction[i-1] = k;
      }
    }
  }

#if (FAC_DEBUG > DEBUG_RECOUPLE)
  fprintf(debug_log, "after sort: %d %d %d %d\n", 
	  interaction[0], interaction[1], 
	  interaction[2], interaction[3]);
#endif
  
  (*idatum)->phase = 0;
  if (interaction[0] >= 0) {
    for (j = interaction[0]+1; j <= interaction[1]; j++) {
      (*idatum)->phase += bra[j].nq;
    }
    for (j = interaction[2]+1; j <= interaction[3]; j++) {
      (*idatum)->phase += bra[j].nq;
    }
  } else if (interaction[2] >= 0) {
    for (j = interaction[2]+1; j <= interaction[3]; j++) {
      (*idatum)->phase += bra[j].nq;
    }
    (*idatum)->phase += 1;
  }

  if (n_shells > 0) {
    i = 0;
    j = 0;
      
    (*sbra) = calloc(n_shells, sizeof(SHELL_STATE));
    (*sket) = calloc(n_shells, sizeof(SHELL_STATE));
    for (m = 0; m < n_shells; m++) {
      if (i < ci->n_shells) {
	if (bra[m].n == ci->shells[i].n &&
	      bra[m].kappa == ci->shells[i].kappa) {
	  memcpy((*sbra)+m, csf_i+i, sizeof(SHELL_STATE));
	  i++;
	} else {
	  (*sbra)[m].totalJ = csf_i[i].totalJ;
	}
      }
      if (j < cj->n_shells) {
	if (bra[m].n == cj->shells[j].n &&
	    bra[m].kappa == cj->shells[j].kappa) {
	  memcpy((*sket)+m, csf_j+j, sizeof(SHELL_STATE));
	  j++;
	} else {
	  (*sket)[m].totalJ = csf_j[j].totalJ;
	}	
      }
    }
  }

 END:
  if (n_shells == 0) n_shells = -1;
  if (n_shells > 0) {
    (*idatum)->bra = (SHELL *) ReallocNew(bra, sizeof(SHELL)*n_shells);
    /* adjust the index so that it counts from inner shells */
    for (i = 0; i < 4; i++) {
      if (s[i].index >= 0) 
	s[i].index = n_shells - 1 - s[i].index;
    }     
  }

  return n_shells;
}

void CompactInteractShell(char c[4], INTERACT_SHELL *s, int m) {
  c[0] = (char) s->n;
  c[1] = (char) s->kl/2;
  if (s->j > s->kl) c[2] = '+';
  else c[2] = '-';
  if (m == 0) {
    c[3] = (char) s->nq_bra;
  } else {
    c[3] = (char) s->nq_ket;
  }
}

/* 
** FUNCTION:    GetInteract
** PURPOSE:     determing which shells can be interacting.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int GetInteract(INTERACT_DATUM **idatum,
		SHELL_STATE **sbra, 
		SHELL_STATE **sket, 
		int kgi, int kgj,
		int kci, int kcj, 
		int ki, int kj, int ifb) {
  int i, j, m;
  CONFIG *ci, *cj, cip;
  SHELL_STATE *csf_i, *csf_j, *csf_ip;
  SHELL *bra;
  INTERACT_SHELL *s;
  int n_shells;
  int index[4];
  LOCK *lock = NULL;
  int locked = 0;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;  
  start = clock();
#endif

  ci = GetConfigFromGroup(kgi, kci);
  cj = GetConfigFromGroup(kgj, kcj);
  if (ci->n_csfs > 0) {
    csf_i = ci->csfs + ki;
    csf_j = cj->csfs + kj;
  } else {
    csf_i = NULL;
    csf_j = NULL;
  }
  if (ci->n_shells <= 0 || cj->n_shells <= 0) return -1;
  if (abs(ci->n_shells+ifb - cj->n_shells) > 2) return -1;

  if (csf_i != NULL) {
    n_shells = -1;
    /* check if this is a repeated call,
     * if not, search in the array.
     */
    if (*idatum == NULL) {
      index[0] = kgi;
      index[1] = kgj;
      index[2] = kci;
      index[3] = kcj;
      (*idatum) = (INTERACT_DATUM *) MultiSet(interact_shells, index, 
					      NULL, &lock, InitInteractDatum, 
					      FreeInteractDatum);
    }
    if (lock && (*idatum)->n_shells == 0) {
      SetLock(lock);
      locked = 1;
    }
    if ((*idatum)->n_shells < 0) {
      if (locked) ReleaseLock(lock);
      return -1;
    }
  } else {
    (*idatum) = malloc(sizeof(INTERACT_DATUM));
    (*idatum)->n_shells = 0;
  }  
  n_shells = (*idatum)->n_shells;
  if (n_shells > 0) {
    bra = (*idatum)->bra;
    s = (*idatum)->s;
    i = 0;
    j = 0;
    (*sbra) = calloc(n_shells, sizeof(SHELL_STATE));
    (*sket) = calloc(n_shells, sizeof(SHELL_STATE));
    for (m = 0; m < n_shells; m++) {
      if (i < ci->n_shells) {
	if (bra[m].n == ci->shells[i].n &&
	    bra[m].kappa == ci->shells[i].kappa) {
	  memcpy((*sbra)+m, csf_i+i, sizeof(SHELL_STATE));
	  i++;
	} else {
	  (*sbra)[m].totalJ = csf_i[i].totalJ;
	}
      }
      if (j < cj->n_shells) {
	if (bra[m].n == cj->shells[j].n &&
	    bra[m].kappa == cj->shells[j].kappa) {
	  memcpy((*sket)+m, csf_j+j, sizeof(SHELL_STATE));
	  j++;
	} else {
	  (*sket)[m].totalJ = csf_j[j].totalJ;
	}	
      }
    }
    /* the quantum number for the free electron must be reset */
    if (ifb) {
      (*sbra)[0].shellJ = 1;
      (*sbra)[0].totalJ = 1;
      (*sbra)[0].nu = 1;
      (*sbra)[0].Nr = 0;
    }
  } else {
    if (ifb) {
      cip.n_shells = ci->n_shells + 1;
      cip.shells = malloc(sizeof(SHELL)*cip.n_shells);
      memcpy(cip.shells+1, ci->shells, sizeof(SHELL)*ci->n_shells);
      cip.shells[0].n = 9999;
      cip.shells[0].nq = 1;
      cip.shells[0].kappa = -1;
      if (csf_i) {
	cip.csfs = malloc(sizeof(SHELL_STATE)*cip.n_shells);
	cip.n_csfs = 1;
	csf_ip = cip.csfs;
	memcpy(csf_ip+1, csf_i, sizeof(SHELL_STATE)*ci->n_shells);
	csf_ip[0].shellJ = 1;
	csf_ip[0].totalJ = 1;
	csf_ip[0].nu = 1;
	csf_ip[0].Nr = 0;
      } else {
	cip.n_csfs = 0;
	csf_ip = NULL;
      }
      n_shells = InteractingShells(idatum, sbra, sket, 
				   &cip, cj, csf_ip, csf_j);
      free(cip.shells);
      if (csf_i) {
	free(cip.csfs);
      }
    } else {
      n_shells = InteractingShells(idatum, sbra, sket, 
				   ci, cj, csf_i, csf_j);
    }
  }

  if (n_shells < 0 && csf_i == NULL) {
    free((*idatum));
    idatum = NULL;
  } else {  
    (*idatum)->n_shells = n_shells;
  }
  if (locked) ReleaseLock(lock);
#pragma omp flush

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.interact += stop -start;
#endif
  return n_shells;
}
    
double EvaluateFormula(FORMULA *fm) {
  double r;

  if (fm->js[0]) {
    CPYDAT(MAXNJGD, fm->njgdata, 1);
  }
  NJSUM(fm->js+1, &r);
  return r;
}

int GenerateFormula(FORMULA *fm) {
  int ij2[MAXJ*3], ij3[MAXJ*3];
  int n, i, j, k, ns;
  
  ns = fm->ns;  
  n = 3*(ns-1);
  k = 0;
  for (i = 1; i < ns; i++) {
    for (j = 1; j <= 3; j++) {
      ij2[k] = fm->tr1[i][j];
      ij3[k] = fm->tr2[i][j];
      k++;
    }
    /*
    printf("%2d: %2d %2d %2d   %2d %2d %2d\n",
	   i, fm->tr1[i][1], fm->tr1[i][2], fm->tr1[i][3], 
	   fm->tr2[i][1], fm->tr2[i][2], fm->tr2[i][3]);
    */
  }
  NJFORM(n, ns, ij2, ij3, fm->ifree+1);
  if (fm->js[0]) {
    CPYDAT(MAXNJGD, fm->njgdata, 0);
  }
  return 0;
}
   
int FixJsZ(INTERACT_SHELL *s, FORMULA *fm) {
  int *js, i, ns;

  ns = fm->ns;
  js = fm->js;
  for (i = 1; i <= ns; i++) {
    js[i] = s[i-1].j;
    if (IsOdd(i)) {
      s[i-1].n = 1;
    } else {
      s[i-1].n = -1;
    }
  }

  return 0;
}

int TriadsZ(int n1, int n2, FORMULA *fm) {
  int ns, i0, i;

  ns = 2*n1 + 2*n2;
  fm->ns = ns;
  for (i = 1; i <= 3*(ns-1); i++) {
    fm->js[i] = 0;
    fm->ifree[i] = 1;
  }
  i0 = n1+n2;
  for (i = 0; i < i0; i++) {
    fm->tr1[i+1][1] = 2*i + 1;
    fm->tr1[i+1][2] = 2*i + 2;
    fm->tr1[i+1][3] = ns + i + 1;
  }

  i0 = ns;
  ns++;
  if (n1 == 1 && n2 == 1) {
    fm->tr1[3][1] = ns;
    fm->tr1[3][2] = ns+1;
    fm->tr1[3][3] = ns+2;
    fm->ifree[ns] = 0;
    fm->ifree[ns+1] = 0;
    fm->ifree[ns+2] = 0;
  } else if (n1 == 1 && n2 == 2) {
    fm->tr1[4][1] = ns+1;
    fm->tr1[4][2] = ns+2;
    fm->tr1[4][3] = ns+3;
    fm->tr1[5][1] = ns;
    fm->tr1[5][2] = ns+3;
    fm->tr1[5][3] = ns+4;
    fm->ifree[ns] = 0;
    fm->ifree[ns+3] = 0;
    fm->ifree[ns+4] = 0;
  } else if (n1 == 2 && n2 == 1) {
    fm->tr1[4][1] = ns;
    fm->tr1[4][2] = ns+1;
    fm->tr1[4][3] = ns+3;
    fm->tr1[5][1] = ns+3;
    fm->tr1[5][2] = ns+2;
    fm->tr1[5][3] = ns+4;
    fm->ifree[ns+2] = 0;
    fm->ifree[ns+3] = 0;
    fm->ifree[ns+4] = 0;
  } else if (n1 == 2 && n2 == 2) {
    fm->tr1[5][1] = ns;
    fm->tr1[5][2] = ns+1;
    fm->tr1[5][3] = ns+4;
    fm->tr1[6][1] = ns+2;
    fm->tr1[6][2] = ns+3;
    fm->tr1[6][3] = ns+5;
    fm->tr1[7][1] = ns+4;
    fm->tr1[7][2] = ns+5;
    fm->tr1[7][3] = ns+6;
    fm->ifree[ns+4] = 0;
    fm->ifree[ns+5] = 0;
    fm->ifree[ns+6] = 0;
  }

  return i0;
}

int CoupleSuccessive(int n, int *ik, int itr, TRIADS tr, int *i0) {
  int i;

  if (n <= 1) return itr;

  tr[itr][1] = ik[0];
  tr[itr][2] = ik[1];
  tr[itr][3] = (*i0);
  itr++;
  (*i0)++;
  for (i = 2; i < n; i++) {
    tr[itr][1] = tr[itr-1][3];
    tr[itr][2] = ik[i];
    tr[itr][3] = (*i0);
    itr++;
    (*i0)++;
  }

  return itr;
}

int RecoupleTensor(int ns, INTERACT_SHELL *s, FORMULA *fm) {
  int ninter, nop;
  int *inter, *interp, *irank, *order;
  int i, j, ij, itr, m;
  int ik[MAXJ], fk[MAXJ];

  order = fm->order;
  inter = fm->inter;
  interp = fm->interp;
  irank = fm->irank;
  
  for (i = 0; i < ns; i++) order[i] = i;
  fm->phase = SortShell(ns, s, order);  
  ninter = 1;
  inter[0] = s[order[0]].index;
  interp[0] = 0;

  i = 0;
  for (i = 1; i < ns; i++) {
    if (s[order[i]].index != s[order[i-1]].index) {
      inter[ninter] = s[order[i]].index;
      interp[ninter] = i;
      ninter++;
    }
  }
  interp[ninter] = ns;
  ij = ninter/2;
  if (IsOdd(ninter)) ij++;
  for (i = 0; i < ij; i++) {
    j = inter[i];
    inter[i] = inter[ninter-i-1];
    inter[ninter-i-1] = j;
  }

  itr = 1;
  ij = 2*ns;
  for (i = 0; i < ninter; i++) {
    nop = interp[i+1] - interp[i];
    if (nop > 1) {
      for (j = interp[i]; j < interp[i+1]; j++) {
	ik[j-interp[i]] = order[j]+1;
      }
      itr = CoupleSuccessive(nop, ik, itr, fm->tr2, &ij);
      fk[i] = fm->tr2[itr-1][3];
    } else {
      fk[i] = order[interp[i]]+1;
    }
  }
  if (ninter > 1) {
    itr = CoupleSuccessive(ninter, fk, itr, fm->tr2, &ij);
  }
  fm->tr2[itr-1][3] = 2*ns-1;

  j = 0;
  irank[j++] = fm->tr2[itr-1][3];
  for (i = 1; i < ninter; i++) {
    irank[j++] = fm->tr2[itr-i][2];
    irank[j++] = fm->tr2[itr-i][1];
  }
  irank[j] = irank[j-1];

  fm->ns = ns;
  fm->ninter = ninter;
  if (GenerateFormula(fm) < 0) return -1;

  return ninter;
}

void EvaluateTensor(int nshells, SHELL_STATE *bra, SHELL_STATE *ket,
		    INTERACT_SHELL *s, int itr, FORMULA *fm) {
  int i, j, k, i0, j0, j1, kmin, kmax, ncp;
  int *js, m, nop, k1, k2, ns, ninter;
  int *order, *inter, *interp, *irank;
  double a, b, a1, a2, *r;
  SHELL_STATE st1, st2;
  RCFP_STATE rcfp_bra, rcfp_ket;
  RCFP_OPERATOR ops[2*MAXJ];
  int rank[2*MAXJ];

  r = &(fm->coeff);
  order = fm->order;
  ninter = fm->ninter;
  inter = fm->inter;
  interp = fm->interp;
  irank = fm->irank;
  ns = fm->ns;
  if (itr == 1) {
    *r = 0.0;
  }
  js = fm->js;
  i0 = 2*ns;
  ncp = ns-1;
  j0 = js[fm->tr2[itr][1]];
  j1 = js[fm->tr2[itr][2]];
  kmin = abs(j0-j1);
  kmax = j0+j1;
  if (fm->tr2[itr][3] >= i0) {
    for (k = kmin; k <= kmax; k += 2) {
      js[fm->tr2[itr][3]] = k;
      for (i = itr+1; i <= ncp; i++) {
	if (fm->tr2[i][1] == fm->tr2[itr][3]) {
	  js[fm->tr2[i][1]] = k;
	  break;
	}	
	if (fm->tr2[i][2] == fm->tr2[itr][3]) {
	  js[fm->tr2[i][2]] = k;
	  break;
	}
      }
      EvaluateTensor(nshells, bra, ket, s, itr+1, fm);
    }
  } else {
    k = js[fm->tr2[itr][3]];
    if (!(kmin <= k && kmax >= k)) return;
    k = 2*ninter;
    for (i = 0; i < k; i++) {
      rank[i] = js[irank[i]];
    }
    a1 = DecoupleShell(nshells, bra, ket, ninter, inter, rank);
    if (fabs(a1) < EPS30) return;
    a2 = EvaluateFormula(fm);
    if (fabs(a2) < EPS30) return;
    i0 = 1;
    a = 1.0;
    for (i = 0; i < ninter; i++) {
      nop = interp[i+1] - interp[i];
      k = order[interp[i]];
      k1 = nshells - s[k].index -1;
      st1 = bra[k1];
      st2 = ket[k1];
      rcfp_bra.state = RCFPTermIndex(s[k].j, st1.nu, st1.Nr, st1.shellJ);
      rcfp_bra.nq = s[k].nq_bra;
      rcfp_bra.subshellMQ = rcfp_bra.nq - (s[k].j + 1)/2;
      rcfp_ket.state = RCFPTermIndex(s[k].j, st2.nu, st2.Nr, st2.shellJ);
      rcfp_ket.nq = s[k].nq_ket;
      rcfp_ket.subshellMQ = rcfp_ket.nq - (s[k].j + 1)/2;
      b = 0.0;
      switch (nop) {
      case 1:
	b = ReducedA(&rcfp_bra, &rcfp_ket, s[k].n);
	break;
      case 2:
	k1 = order[interp[i]+1];
	b = ReducedW(&rcfp_bra, &rcfp_ket, js[fm->tr2[i0][3]]/2, s[k].n, s[k1].n);
	break;
      case 3:
	k1 = order[interp[i]+1];
	k2 = order[interp[i]+2];
	b = ReducedWxA(&rcfp_bra, &rcfp_ket, js[fm->tr2[i0][3]]/2, 
		       js[fm->tr2[i0+1][3]], s[k].n, s[k1].n, s[k2].n);
	break;
      default:
	for (m = 0; m < nop; m++) {
	  ops[m].rank = s[k].j;
	  ops[m].left = NULL;
	  ops[m].right = NULL;
	  ops[m].nops = 1;
	  k1 = order[interp[i]+m];
	  ops[m].qm = s[k1].n;
	}
	CoupleOperators(&ops[0], &ops[1], &ops[nop], js[fm->tr2[i0][3]]);
	for (m = 2; m < nop; m++) {
	  CoupleOperators(&ops[nop+m-2], &ops[m], &ops[nop+m-1], 
			  js[fm->tr2[i0+m-1][3]]);
	}
	m = 2*nop-2;
	b = ReducedOperator(&rcfp_bra, &rcfp_ket, &ops[m]);
	break;
      }
      if (fabs(b) < EPS30) return;
      a *= b;
      i0 += nop - 1;
    }
    if (IsOdd(fm->phase)) a = -a;    
    (*r) += a*a1*a2;
  }

  return;
}

void TestAngular(void) {
  int i, j, k, ns, t;  
  SYMMETRY *sym;  
  CONFIG *c1, *c2;
  ARRAY *s1, *s2;
  STATE *st1, *st2;
  int n_shells;
  SHELL *bra;
  SHELL_STATE *sbra, *sket;
  INTERACT_DATUM *idatum;
  INTERACT_SHELL *s;
  INTERACT_SHELL es[4];
  char name1[LEVEL_NAME_LEN];
  char name2[LEVEL_NAME_LEN];
  int phase, q;

  ns = MAX_SYMMETRIES;
  
  for (i = 0; i < ns; i++) {
    sym = GetSymmetry(i);
    s1 = &(sym->states);
    s2 = s1;
    for (t = 0; t < sym->n_states; t++) {
      st1 = ArrayGet(s1, t);
      c1 = GetConfig(st1);      
      ConstructLevelName(name1, NULL, NULL, NULL, st1);
      for (q = 0; q < sym->n_states; q++) {
	st2 = ArrayGet(s2, q);
	c2 = GetConfig(st2);
	ConstructLevelName(name2, NULL, NULL, NULL, st2);
	printf("Bra: %s \nKet: %s\n", name1, name2);
	n_shells = GetInteract(&idatum, &sbra, &sket, 
			       st1->kgroup, st2->kgroup, 
			       st1->kcfg, st2->kcfg, 
			       st1->kstate, st2->kstate, 0);
	phase = idatum->phase;
	s = idatum->s;
	bra = idatum->bra;
	printf("n_shells: %d \n", n_shells);
	if (n_shells < 0) continue;
	if (s[0].index >= 0 && s[3].index >= 0) {
	  CheckAngularConsistency(n_shells, bra, sbra, sket, s, phase);
	  /*
	  if (s[0].index != s[2].index &&
	      s[1].index != s[3].index) {
	    memcpy(es, s, sizeof(INTERACT_SHELL));
	    memcpy(es+1, s+3, sizeof(INTERACT_SHELL));
	    memcpy(es+2, s+2, sizeof(INTERACT_SHELL));
	    memcpy(es+3, s+1, sizeof(INTERACT_SHELL));
	    CheckAngularConsistency(n_shells, bra, sbra, sket, es, phase);
	  }
	  */
	} else if (s[0].index >= 0) {
	  for (j = 0; j < n_shells; j++) {
	    s[2].index = n_shells - j - 1;
	    s[3].index = s[2].index;
	    s[2].n = bra[j].n;
	    s[3].n = s[2].n;
	    s[2].kappa = bra[j].kappa;
	    s[3].kappa = s[2].kappa;
	    s[2].j = GetJ(bra+j);
	    s[3].j = s[2].j;
	    s[2].kl = GetL(bra+j);
	    s[3].kl = s[2].kl;
	    s[2].nq_bra = GetNq(bra+j);
	    if (s[2].index == s[0].index) {
	      s[2].nq_ket = s[2].nq_bra - 1;
	    } else if (s[2].index == s[1].index) {
	      s[2].nq_ket = s[2].nq_bra + 1;
	    } else {
	      s[2].nq_ket = s[2].nq_bra;
	    }
	    s[3].nq_bra = s[2].nq_bra;
	    s[3].nq_ket = s[2].nq_ket;

	    printf("In Pfac: %d %d  %d %d  %d %d  %d %d\n",
		    s[0].nq_bra, s[0].nq_ket, s[1].nq_bra, s[1].nq_ket, 
		    s[2].nq_bra, s[2].nq_ket, s[3].nq_bra, s[3].nq_ket);
		    
	    
	    CheckAngularConsistency(n_shells, bra, sbra, sket, s, phase);
	    /*
	    memcpy(es, s+2, sizeof(INTERACT_SHELL));
	    memcpy(es+1, s+1, sizeof(INTERACT_SHELL));
	    memcpy(es+2, s, sizeof(INTERACT_SHELL));
	    memcpy(es+3, s+3, sizeof(INTERACT_SHELL));

	    printf("In Pfac: %d %d  %d %d  %d %d  %d %d\n",
		    s[0].nq_bra, s[0].nq_ket, s[1].nq_bra, s[1].nq_ket, 
		    s[2].nq_bra, s[2].nq_ket, s[3].nq_bra, s[3].nq_ket);
	    printf("In Pfac ex: %d %d  %d %d  %d %d  %d %d\n",
		    s[0].nq_bra, es[0].nq_ket, es[1].nq_bra, es[1].nq_ket, 
		    s[2].nq_bra, es[2].nq_ket, es[3].nq_bra, es[3].nq_ket);
	    CheckAngularConsistency(n_shells, bra, sbra, sket, es, phase);
	    */
	  }
	} else {
	  for (j = 0; j < n_shells; j++) {
	    s[0].index = n_shells - j - 1;
	    s[1].index = s[0].index;
	    s[0].n = bra[j].n;
	    s[1].n = s[0].n;
	    s[0].kappa = bra[j].kappa;
	    s[1].kappa = s[0].kappa;
	    s[0].j = GetJ(bra+j);
	    s[1].j = s[0].j;
	    s[0].kl = GetL(bra+j);
	    s[1].kl = s[0].kl;
	    s[0].nq_bra = GetNq(bra+j);
	    s[0].nq_ket = s[0].nq_bra;
	    s[1].nq_bra = s[0].nq_bra;
	    s[1].nq_ket = s[1].nq_bra;
	    for (k = 0; k < n_shells; k++) {
	      s[2].index = n_shells - k - 1;
	      s[3].index = s[2].index;
	      s[2].n = bra[k].n;
	      s[3].n = s[2].n;
	      s[2].kappa = bra[k].kappa;
	      s[3].kappa = s[2].kappa;
	      s[2].j = GetJ(bra+k);
	      s[3].j = s[2].j;
	      s[2].kl = GetL(bra+k);
	      s[3].kl = s[2].kl;
	      s[2].nq_bra = GetNq(bra+k);
	      s[2].nq_ket = s[2].nq_bra;
	      s[3].nq_bra = s[2].nq_bra;
	      s[3].nq_ket = s[3].nq_bra;
	      
	      CheckAngularConsistency(n_shells, bra, sbra, sket, s, phase);
	      /*
	      memcpy(es, s+2, sizeof(INTERACT_SHELL));
	      memcpy(es+1, s+1, sizeof(INTERACT_SHELL));
	      memcpy(es+2, s, sizeof(INTERACT_SHELL));
	      memcpy(es+3, s+3, sizeof(INTERACT_SHELL));
	      CheckAngularConsistency(n_shells, bra, sbra, sket, es, phase);
	      */
	    }
	  }
	}
	free(sbra);
	free(sket);
      }
    }
  }
}

void CheckAngularConsistency(int n_shells, SHELL *bra, 
			     SHELL_STATE *sbra, SHELL_STATE *sket,
			     INTERACT_SHELL *ss, int phase) {
  int nk1, *kk1, nk2, *kk2, nk3, *kk3, i, j, k, nk0, *k0;
  int k1, k2, k3, j1, j2, j3, imin, imax, p, p2, p3, err;
  double *ang1, *ang2, *ang3, *ang4, z0, *x, y;
  CONFIG icfg;
  INTERACT_SHELL s[4];
  SHELL *t;
  int maxq;
  FORMULA fm;

  fm.js[0] = 0;
  memcpy(s, ss, sizeof(INTERACT_SHELL)*4);
  nk1 = AngularZxZ0(&ang1, &kk1, 0, n_shells, sbra, sket, s);
  if (nk1 > 0) {
    TriadsZ(1, 1, &fm);
    RecoupleTensor(4, s, &fm);
    FixJsZ(s, &fm);
    kk2 = fm.js;
    printf("js = ");
    for (j = 1; j <= 9; j++) printf("%d ", kk2[j]);
    printf("\n");
    printf("%d %d %d %d\n", s[0].index, s[1].index, s[2].index, s[3].index);
    for (i = 0; i < nk1; i++) {
      k = kk1[i];
      kk2[5] = k;
      kk2[6] = k;
      kk2[7] = 0;      
      EvaluateTensor(n_shells, sbra, sket, s, 1, &fm);
      y = fm.coeff;
      y /= sqrt(k+1.0);
      if (IsOdd(k/2)) y = -y;
      printf("Ang: %2d %15.8E %15.8E\n", k, ang1[i], y);
    }
  }

#if FAC_DEBUG
  memcpy(s, ss, sizeof(INTERACT_SHELL)*4);

  nk1 = AngularZxZ0(&ang1, &kk1, 0, n_shells, sbra, sket, s);

  if (nk1 <= 0) {
    fprintf(debug_log, "no angular coeff.\n\n");
    return;
  }
  ang4 = malloc(sizeof(double)*nk1);
  for (i = 0; i < nk1; i++) ang4[i] = 0.0;

  z0 = 0;
  if (s[1].index == s[2].index) {
    nk0 = 1;
    k = 0;
    k0 = &k;
    x = &z0;
    nk0 = AngularZ(&x, &k0, nk0, n_shells, sbra, sket, s, s+3);
    z0 = z0/sqrt(s[0].j+1.0);
    if (IsOdd((s[0].j + s[2].j)/2)) z0 = -z0;
  }
  
  t = (SHELL *) malloc(sizeof(SHELL)*n_shells);
  if (t == NULL) {
    free(ang1);
    free(ang4);
    free(kk1);
    exit(1);
  }

  icfg.shells = t;
  icfg.n_shells = n_shells;
  memcpy(t, bra, sizeof(SHELL)*n_shells);
  icfg.shells[n_shells-1-s[0].index].nq -= 1;
  icfg.shells[n_shells-1-s[1].index].nq += 1;
  s[0].nq_ket = icfg.shells[n_shells-1-s[0].index].nq;
  s[1].nq_ket = icfg.shells[n_shells-1-s[1].index].nq;
  s[2].nq_bra = icfg.shells[n_shells-1-s[2].index].nq;
  s[3].nq_bra = icfg.shells[n_shells-1-s[3].index].nq;


  fprintf(debug_log, "Interacting Shells: %d %d %d %d\n", 
	  s[0].index, s[1].index, s[2].index, s[3].index);
  fprintf(debug_log, "Intermediate Config.: ");
  for (i = n_shells-1; i >= 0; i--) {
    fprintf(debug_log, "(%d %d %d), ", icfg.shells[i].n, 
	   icfg.shells[i].kappa, icfg.shells[i].nq);
  }
  fprintf(debug_log, "\n");
  
  for (i = 0; i < n_shells; i++) {
    maxq = 2*abs(icfg.shells[i].kappa);
    if (icfg.shells[i].nq < 0 || icfg.shells[i].nq > maxq) {
      fprintf(debug_log, "\n");
      for (i = 0; i < nk1; i++) {
	if (IsOdd(p)) ang4[i] = -ang4[i];
	fprintf(debug_log, "Rank: %d, ZxZ0: %12.8lf, %12.8lf; Z0: %12.8lf \n", 
		kk1[i]/2, ang1[i], ang4[i], z0);    
	if (fabs(ang1[i] - ang4[i]) > EPS10) {
	  fprintf(debug_log, "***********Error*********\n");
	}
      }
      fprintf(debug_log, "No CSFs for this intermediate Config.\n");
      fprintf(debug_log, "\n");
      free(ang1);
      free(ang4);
      free(kk1);
      free(t);
      return;
    }
  }

  err = Couple(&icfg);

  if (err < 0) {  
    fprintf(debug_log, "\n");
    for (i = 0; i < nk1; i++) {
      if (IsOdd(p)) ang4[i] = -ang4[i];
      fprintf(debug_log, "Rank: %d, ZxZ0: %12.8lf, %12.8lf; Z0: %12.8lf \n", 
	      kk1[i]/2, ang1[i], ang4[i], z0);    
      if (fabs(ang1[i] - ang4[i]) > EPS10) {
	fprintf(debug_log, "***********Error*********\n");
      }
    }
    fprintf(debug_log, "No CSFs for this intermediate Config.\n");
    fprintf(debug_log, "\n");
    free(ang1);
    free(ang4);
    free(kk1);
    free(t);
    return;
  }
  
  if (s[0].index < s[1].index) {
    imin = s[0].index;
    imax = s[1].index;
  } else {
    imin = s[1].index;
    imax = s[0].index;
  }

  p = 0;
  if (s[2].index < imax && s[2].index >= imin) p += 1;
  if (s[3].index < imax && s[3].index >= imin) p += 1;
  
  p2 = 0;
  for (i = imin; i < imax; i++) {
    p2 += icfg.shells[n_shells-1-i].nq;
  }

  if (s[2].index < s[3].index) {
    imin = s[2].index;
    imax = s[3].index;
  } else {
    imin = s[3].index;
    imax = s[2].index;
  }

  if (imin == imax) p3 = 0; 
  else {
    p3 = 1;
    for (i = imin; i < imax; i++) {
      p3 += icfg.shells[n_shells-1-i].nq;
    }
  }
  p = (IsOdd(p)?-1:1);
  p2 = (IsOdd(p2)?-1:1);
  p3 = (IsOdd(p3)?-1:1);
  phase = (IsOdd(phase)?-1:1);
  
  for (i = 0; i < n_shells*icfg.n_csfs; i += n_shells) {
    fprintf(debug_log, "bra state: ");
    for (k = n_shells-1; k >= 0; k--) {
      fprintf(debug_log, "(%d %d), ", sbra[k].shellJ, sbra[k].totalJ);
    }
    fprintf(debug_log, "\n");
    
    fprintf(debug_log, "ket state: ");
    for (k = n_shells-1; k >= 0; k--) {
      fprintf(debug_log, "(%d %d), ", sket[k].shellJ, sket[k].totalJ);
    }
    fprintf(debug_log, "\n");    
    
    fprintf(debug_log, "Intermediate State %d: ", i);
    for (k = n_shells-1; k >= 0; k--) {
      fprintf(debug_log, "(%d %d), ",  
	      icfg.csfs[i+k].shellJ, icfg.csfs[i+k].totalJ);
    }
    fprintf(debug_log, "\n");
    
    nk2 = AngularZ(&ang2, &kk2, 0, n_shells, sbra, icfg.csfs+i, 
		   s, s+1);
    if (nk2 <= 0) continue;

    nk3 = AngularZ(&ang3, &kk3, 0, n_shells, icfg.csfs+i, sket, 
		   s+2, s+3);
    if (nk3 <= 0) {
      free(ang2);
      free(kk2);
      continue;
    }

    if (kk2[0] > kk3[0]) {
      k2 = 0;
      k3 = (kk2[0] - kk3[0])/2;
    } else {
      k2 = (kk3[0] - kk2[0])/2;
      k3 = 0;
    } 

    if (k2 >= nk2 || k3 >= nk3) {
      free(ang2);
      free(kk2);
      free(ang3);
      free(kk3);
      continue;
    }

    if (kk2[k2] > kk1[0]) {
      k1 = (kk2[k2] - kk1[0]) / 2;
    } else {
      k1 = 0;
      k = (kk1[0] - kk2[k2])/2;
      k2 += k;
      k3 += k;
    } 
 
    if (k1 >= nk1 || k2 >= nk2 || k3 >= nk3) {
      free(ang2);
      free(ang3);
      free(kk2);
      free(kk3);
      continue;
    }

    fprintf(debug_log, 
	    "Sum Range: 1: %d %d,  2: %d %d, 3: %d %d .. %d %d %d\n", 
	    k1, nk1, k2, nk2, k3, nk3, kk1[k1], kk2[k2], kk3[k3]);

    for (j1 = k1, j2 = k2, j3 = k3; 
	 j1 < nk1 && j2 < nk2 && j3 < nk3;
	 j1++, j2++, j3++) {
      y = ang2[j2] * ang3[j3];
      if (IsOdd((sbra[0].totalJ - icfg.csfs[i].totalJ) / 2)) y = -y;
      fprintf(debug_log, "Sum Terms; %d %d %d %lf\n",
	      kk1[j1], kk2[j2], kk3[j3], y);
      ang4[j1] += y/sqrt(sbra[0].totalJ + 1.0);
    }
    free(ang2);
    free(ang3);
    free(kk2);
    free(kk3);
  }
  fprintf(debug_log, "\n");
  for (i = 0; i < nk1; i++) {
    fprintf(debug_log, 
	    "Rank: %d, ZxZ0: %+d %10.6lf;  %+d %+d %10.6lf; %+d \n", 
	    kk1[i]/2, phase, ang1[i], p2, p3, ang4[i], p);    
    if (fabs(phase*ang1[i] - p2*p3*ang4[i]) > EPS10) {
      fprintf(debug_log, "***********Error*********\n");
    }
    if (phase != p*p2*p3) {
      fprintf(debug_log, "Phase Error\n");
    }
  }
  fprintf(debug_log, "\n");

  free(ang4);
  free(t);
  free(icfg.csfs);
#endif /* FAC_DEBUG */

  free(kk1);
  free(ang1);
}

/* 
** FUNCTION:    InitRecouple
** PURPOSE:     Initialize the module "recouple"
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int InitRecouple(void) {
  int blocks[4] = {10, 10, 64, 64};
  int ndim = 4;
  
  FACTT();
  interact_shells = (MULTI *) malloc(sizeof(MULTI));
  return MultiInit(interact_shells, sizeof(INTERACT_DATUM),
		   ndim, blocks, "interact_shells");
}

/* 
** FUNCTION:    ReinitRecouple
** PURPOSE:     Reinitialize the module "recouple"
** INPUT:       {int m},
**              >=0: do a full reinitialization.
**               <0: do nothing.
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int ReinitRecouple(int m) {
  if (m < 0) return 0;
#pragma omp barrier
#pragma omp master
  MultiFreeData(interact_shells, FreeInteractDatum);  
#pragma omp barrier
  return 0;
}
  
