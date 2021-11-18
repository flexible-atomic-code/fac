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

#ifndef _ANGULAR_H_
#define _ANGULAR_H_ 1

/*************************************************************
  Header for module "angular". 
  This module calculates Wigner 3j, 6j, 9j symbols, 
  and related vector coupling coefficients.
  
  This module is mainly translated from the F90 package of
  Gaigalas et al. CPC 139 (2001) 263. with slight modification 
  for the Wigner 3j symbol which avoids the overflow for 
  large angular momenta.

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

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "global.h"

/*
** VARIABLE:    MAX_FACTORIAL = 1000
** TYPE:        macro constant.
** PURPOSE:     the maximum integer whose ln and ln(factorial) 
**              is tabulated.
** NOTE:        
*/
#define MAX_FACTORIAL 1000

/*
** VARIABLE:    ln_factorial[MAX_FACTORIAL]
** TYPE:        global array.
** PURPOSE:     holds the ln(factorial) of integers.
** NOTE:        
*/
extern double ln_factorial[MAX_FACTORIAL];

/*
** VARIABLE:    ln_integer[MAX_FACTORIAL]
** TYPE:        global array.
** PURPOSE:     holds the ln of integers.
** NOTE:        
*/
extern double ln_integer[MAX_FACTORIAL];

/* 
** MACRO:       LnFactorial
** PURPOSE:     access the ln_factorial array.
** INPUT:       {int n},
**              argument of ln_factorial.
** RETURN:      {double},
**              result.
** SIDE EFFECT: 
** NOTE:        
*/
#define LnFactorial(n) ln_factorial[(n)]

/* 
** MACRO:       LnInteger
** PURPOSE:     access the ln_integer array.
** INPUT:       {int n},
**              argument of ln.
** RETURN:      {int},
**              result.
** SIDE EFFECT: 
** NOTE:        
*/
#define LnInteger(n) ln_integer[(n)]

#ifdef PERFORM_STATISTICS
/*
** STRUCT:      ANGULAR_TIMING
** PURPOSE:     tracking the time spent in the angular module.
** FIELDS:      {clock_t w3j},
**              time spent in W3j.
**              {clock_t w6j},
**              time spent in w6j.
**              {clock_t w9j},
**              time spent in w9j.
** NOTE:        this is used for profiling. 
**              it is only compiled in when the macro 
**              PERFORM_STATISTICS is defined in "global.h".
*/
typedef struct _ANGULAR_TIMING_ {
  clock_t w3j;
  clock_t w6j;
  clock_t w9j;
} ANGULAR_TIMING;

int    GetAngularTiming(ANGULAR_TIMING *t);
#endif

/*
** Public functions provided by *angular*
*/
int    InitAngular(void);
int    Triangle(int j1, int j2, int j3);
double W3j(int j1, int j2, int j3, int m1, int m2, int m3);
double W6j(int j1, int j2, int j3, int i1, int i2, int i3);
int    W6jTriangle(int j1, int j2, int j3, int i1, int i2, int i3);
double W9j(int j1, int j2, int j3,
	   int i1, int i2, int i3,
	   int k1, int k2, int k3);
int    W9jTriangle(int j1, int j2, int j3,
		   int i1, int i2, int i3,
		   int k1, int k2, int k3);
double WignerEckartFactor(int jf, int k, int ji,
			  int mf, int q, int mi);
double ClebschGordan(int j1, int m1, int j2, int m2, int jf, int mf);
double ReducedCL(int ja, int k, int jb);
double WignerDMatrix(double a, int j2, int m2, int n2);

#endif
