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

#ifndef _GRID_H_
#define _GRID_H_ 1

#include "global.h"

int AddPW(int *nkl0, double *kl, double *logkl, 
	  int maxkl, int n, int step);
int SetPWGrid(int *nkl0, double *kl, double *logkl, 
	      int maxkl, int *ns, int *n, int *step);
int SetTEGridDetail(double *te, double *logte, int n, double *x);
int SetTEGrid(double *te, double *logte, int n, double emin, double emax);
int SetEGridDetail(double *e, double *log_e, int n, double *xg);
int SetEGrid(double *e, double *log_e, 
	     int n, double emin, double emax, double eth);
int SetLinearGrid(double *x, int n, double xmin, double xmax);

#endif
