#ifndef _TRANSITION_H_
#define _TRANSITION_H_

#include "global.h"
#include "nucleus.h"
#include "angular.h"
#include "config.h"
#include "structure.h"

int SetTransitionCut(double c);
double GetTransitionCut(void);
void SetTransitionMode(int m);
void SetTransitionGauge(int m);
void SetTransitionMaxE(int m);
void SetTransitionMaxM(int m);
void SetTransitionOptions(int gauge, int mode, int max_e, int max_m);
int GetTransitionGauge(void);
int GetTransitionMode(void);
int TRMultipole(double *strength, double *energy,
		int m, int low, int up);
int OverlapLowUp(int nlow, int *low, int nup, int *up);
int SaveTransition(int nlow, int *low, int nup, int *up,
		   char *fn, int multipole);
int GetLowestMultipole(int p1, int j1, int p2, int j2);


#endif
