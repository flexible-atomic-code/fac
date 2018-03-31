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

#ifndef _COULOMB_H_
#define _COULOMB_H_ 1

/*************************************************************
  Header for module "coulomb". 
  This module calculates quatities related to the H-like ions.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

/* 
<** The following format is used for documenting the source **>
*/

/* documenting a struct */
/*
** STRUCT:      
** PURPOSE:     
** FIELDS:      
** NOTE:        
*/

/* documenting a function */
/* 
** FUNCTION:    
** PURPOSE:     
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/

/* documenting a macro function */
/* 
** MACRO:       
** PURPOSE:     
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/

/* documenting a global, static varialbe or a macro constant */
/*
** VARIABLE:    
** TYPE:        
** PURPOSE:     
** NOTE:        
*/

#include "global.h"
#include "array.h"
#include "angular.h"
#include "grid.h"
#include "orbital.h"
#include "radial.h"

double NucleusRRMS(double z);
void    SetHydrogenicNL(int n, int kl, int nm, int klm);
void    GetHydrogenicNL(int *n, int *kl, int *nm, int *klm);
double  HydrogenicDipole(double z, int n0, int kl0, 
			int n1, int kl1);
double HydrogenicExpectation(double z, int m, int n, int kl);
double HydrogenicSelfEnergy(int md, int pse, double scl,
			    POTENTIAL *pot, ORBITAL *orb, ORBITAL *orbp);
double  TRRateHydrogenic(double z, int n0, int kl0,
			int n1, int kl1, int s);
double  CoulombPhaseShift(double z, double e, int kappa);
int CoulombMultip(char *fn, double z, double te, double e1,
		  int k, int q0, int q1, int m);
double *GetCoulombBethe(int ie2, int ite, int ie1, int t, int q);
double  GetCoulombBetheAsymptotic(double te, double e1);
void    PrepCBIndex(void);
int     CoulombBetheTail(int n, double *w, int nkl, double *kl, double *tcb);
int     PrepCoulombBethe(int ne2, int nte, int ne1, double z,
			 double *e2, double *te, double *e1,
			 int nkl, double *kl, int mode);
int     CoulombBethe(char *s, double z, double te, double e1);
int     InitCoulomb(void);

#endif
