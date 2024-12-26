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

#ifndef _CFP_H_
#define _CFP_H_ 1

/*************************************************************
  Header for module "cfp"

  This module calculates the coefficients of fractional 
  parentage for relativistic subshell with j <= 9/2, by 
  looking up the table. The numerical values in the table 
  are deduced from the book, "Nuclear Shell Theory", 
  A. de-Shalit and I. Talmi. New York, Academic Press, 1963.

  the routines provided are actually no longer used in FAC.
  instead, the reduced coefficients of fractional parentage
  approach of Gaigalas et al. CPC 139 (2001) 263. is used.
  See the module "rcfp.c"

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

#include <stdio.h>
#include <math.h>
#include "global.h"

/*
** public functions provided by "cfp"
*/
int CFP(double *coeff, int j2, int q, int dj, 
	int dw, int pj, int pw);
int GetIndex(int j2, int q, int tj2, int w);

#endif

