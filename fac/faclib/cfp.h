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

  Author: M. F. Gu, mfgu@space.mit.edu
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

