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

#define NEP  5
#define NEPP 2
typedef struct _NUCLEUS_ {
  char symbol[5];
  double z0, atomic_number;
  double m0, mass;
  double rn, z1, rms, rms0, rmse;
  double a, b, c;
  double rfermi[5];
  int nep;
  int epm[NEP];
  double epp[NEP][NEPP];
  int nepr, cepr;
  double *epr, *epv;
} NUCLEUS;

typedef struct _CXTGT_ {
  char symbol[128];
  double z, m, a, b, e, x;
} CXTGT;
  
void PrintNucleus(int m, char *fn);
int InitNucleus(void);
double GraspRRMS(double z, double m);
int SetAtom(char *s, double z, double mass, double rn, double a, double rse);
void SetAtomicChargeDist(double a, double rse);
void ScaleAtomicChargeDist(double a);
char *GetAtomicSymbolTable(void);
double *GetAtomicMassTable(void);
double GetAtomicNumber(void);
double GetAtomicMass(void);
double GetAtomicR(void);
char *GetAtomicSymbol(void);
double GetAtomicEffectiveZ(double r);
double GetExtraZ(double r, int iep);
void SetExtraPotential(int m, int n, double *p, char *fn);
double GetAtomicChargeDist(double r);
NUCLEUS *GetAtomicNucleus(void);
int SetCXTarget(char *s, double a, double b, double e, double x,
		double z, double m);
CXTGT *GetCXTarget(void);
void PrintCXTarget(char *fn);
double NucleusRadius(double z, double m, int md);
void SetupFermi(void);
double FermiRMS(double c, double a);
double FermiParamC(double r, double a);
void SetOptionNucleus(char *s, char *sp, int ip, double dp);
int GetGround2J(int z, int k, int md);
int GetGroundParity(int z, int k, int md);
double GetGroundIP(int z, int k, int md);
char *GetGroundLev(int z, int k, int md);
char *GetGroundCfg(int z, int k, int md);

#endif

