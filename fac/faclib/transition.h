#ifndef _TRANSITION_H_
#define _TRANSITION_H_

#include "global.h"
#include "dbase.h"
#include "nucleus.h"
#include "angular.h"
#include "config.h"
#include "structure.h"

int SetTransitionCut(double c);
double GetTransitionCut();
void SetTransitionOption(int gauge, int mode, 
			 int max_e, int max_m);
int TransitionGauge();
int OscillatorStrength(double *strength, double *energy,
			  int *multipole, int low, int up);
int SaveTransition(int nlow, int *low, int nup, int *up,
		   char *fn, int multipole);

/*
int SaveTransitionGroups(char *fn, int multipole, int glow, int gup); 
*/

int GetLowestMultipole(int p1, int j1, int p2, int j2);
double StrengthInPureCoupling(double aw, int m, STATE *low, STATE *up);
double MultipoleReduced(double aw, int m, int n_shells, 
			SHELL_STATE *sbra, SHELL_STATE *sket, 
			INTERACT_SHELL *s1, INTERACT_SHELL *s2);


#endif
