#ifndef _TRANSITION_H_
#define _TRANSITION_H_

#include "global.h"
#include "nucleus.h"
#include "angular.h"
#include "config.h"
#include "structure.h"

int SetTransitionCut(double c);
double GetTransitionCut();
void SetTransitionOption(int gauge, int mode, int max_e, int max_m);
int GetTransitionGauge();
int GetTransitionMode();
int OscillatorStrength(double *strength, double *energy, double aw0,
		       int m, int low, int up);
int SaveTransition(int nlow, int *low, int nup, int *up,
		   char *fn, int multipole);
int GetLowestMultipole(int p1, int j1, int p2, int j2);


#endif
