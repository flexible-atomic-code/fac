#ifndef _CFP_H_
#define _CFP_H_

#include <stdio.h>
#include <math.h>
#include "global.h"

/* this file implement the table look up of CFPs. this is obsolete 
   since the rcfp package can be used insteak */
int CFP(double *coeff, int j2, int q, int dj, 
	int dw, int pj, int pw);
int GetIndex(int j2, int q, int tj2, int w);

#endif

