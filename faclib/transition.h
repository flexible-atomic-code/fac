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

#ifndef _TRANSITION_H_
#define _TRANSITION_H_

#include "global.h"
#include "nucleus.h"
#include "angular.h"
#include "config.h"
#include "structure.h"

int SetTransitionCut(double c0, double c);
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
int TRMultipoleEB(double *strength, double *energy,
		  int m, int lower, int upper);
int OverlapLowUp(int nlow, int *low, int nup, int *up);
int SaveTransition(int nlow, int *low, int nup, int *up,
		   char *fn, int multipole);
int SaveTransitionEB(int nlow, int *low, int nup, int *up,
		     char *fn, int multipole);
int GetLowUpEB(int *nlow, int **low, int *nup, int **up, 
	       int nlow0, int *low0, int nup0, int *up0);
int PolarizeCoeff(char *ifn, char *ofn, int i0, int i1);

#endif
