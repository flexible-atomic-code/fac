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

#ifndef _NUCLEUS_H_
#define _NUCLEUS_H_

#include "global.h"
#include "interpolation.h"

typedef struct _NUCLEUS_ {
  char symbol[5];
  double z0, atomic_number;
  double m0, mass;
  double rn, z1, rms, rms0;
  double a, b, c;
  int npr;
  double rfermi[5];
} NUCLEUS;

void PrintNucleus();
int InitNucleus();
double GraspRRMS(double z, double m);
int SetAtom(char *s, double z, double mass, double rn, double a, int npr);
char *GetAtomicSymbolTable(void);
double *GetAtomicMassTable(void);
double GetAtomicNumber(void);
double GetAtomicMass(void);
double GetAtomicR(void);
char *GetAtomicSymbol(void);
double GetAtomicEffectiveZ(double r);
double GetAtomicChargeDist(double r);
NUCLEUS *GetAtomicNucleus(void);

#endif

