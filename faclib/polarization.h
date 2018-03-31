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

#ifndef _POLARIZATION_H_
#define _POLARIZATION_H_

#include "dbase.h"

typedef struct _MLEVEL_ {
  short nele;
  short j;
  short p;
  double energy;
  double dtotal;
  double *rtotal;
  double *pop;
  double *npop;
  int ic;
} MLEVEL;

typedef struct _MTR_ {
  int multipole;
  int lower;
  int upper;
  int n;
  double rtotal;
  double *rates;
} MTR;

typedef struct _MCE_ {
  int lower;
  int upper;
  int n;
  double *rates;
} MCE;

typedef struct _MAI_ {
  int f;
  int b;
  int n;
  double *rates;
} MAI;


int InitPolarization(void);
int SetMaxLevels(int m);
int SetMIteration(double a, int m);
int SetEnergy(double energy, double esigma);
int SetDensity(double eden);
int SetIDR(int idr, int ndr, double *pdr);
int SetMLevels(char *fn, char *tfn);
int SetMCERates(char *fn);
int SetMAIRates(char *fn);
int PopulationTable(char *fn);
int Orientation(char *fn, double e);
int PolarizationTable(char *fn, char *ifn, int n, char **sc);

#endif
