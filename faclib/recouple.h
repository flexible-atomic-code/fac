#ifndef _RECOUPLE_H_
#define _RECOUPLE_H_

#include "config.h"
#include "angular.h"
#include "rcfp.h"
#include <time.h>

/* qunatum # of an interacting shell */
typedef struct _INTERACT_SHELL_ {
  int index; /* the index of the shell within the SHELL_STATE */
  int n;
  int j; /* the angular momentum of the shell, double of its actual value */
  int kl;
  int kappa;
  int nq_bra; /* the occupation number in the bra state */
  int nq_ket;
} INTERACT_SHELL;

typedef struct _RECOUPLE_TIMING_ {
  clock_t angz;
  clock_t angzxz;
  clock_t decouple;
  clock_t interact;
} RECOUPLE_TIMING;

typedef struct _INTERACT_DATUM_ {
  SHELL *bra;
  short s[4];
  short n_shells;
  short phase;
} INTERACT_DATUM;

int GetRecoupleTiming(RECOUPLE_TIMING *t);

/* the recoupling matrix going from the coupled operators to uncoupled ones */
double DecoupleShell(int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
		     int n_interact, int *interact, int *rank);

/* check if the recoupling matrix is non-zero */
int IsShellInteracting(int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
		       int n_interact, int *interact, int *rank);

/* the coeff of type (Z^k dot Z^k) */
int AngularZxZ0(double **coeff, int **kk, int nk, 
		int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
		INTERACT_SHELL *s);

/* the coeff of type Z^k */
int AngularZ(double **coeff, int **kk, int nk,
	     int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
	     INTERACT_SHELL *s1, INTERACT_SHELL *s2);


void SumCoeff(double *coeff,  int *kk,  int nk,  int p, 
	      double *coeff1, int *kk1, int nk1, int p1, 
	      int phase, int j1, int j2, int j3, int j4);

int SortShell(INTERACT_SHELL *s, int *order);

/* analyze the structure of the configuration of bra and ket to 
   determine if they can interact, get the interacting shells that
   must interact, and determine the phase factor of the recoupling 
   coeff. which does not include the phase resulting from the reordering
   of operators. that is calculated in the AngularZxZ0 and AngularZ0 */
int GetInteract(int *phase, INTERACT_SHELL *s, SHELL **bra, 
		SHELL_STATE **sbra, SHELL_STATE **sket, CONFIG *ci, 
		int ki, CONFIG *cj, int kj, STATE *s1, STATE *s2);

void TestAngular();
void CheckAngularConsistency(int n_shells, SHELL *bra, 
			     SHELL_STATE *sbra, SHELL_STATE *sket,
			     INTERACT_SHELL *s, int phase);
int SetMaxRank(int k);
int GetMaxRank();
int InitRecouple();

#endif




