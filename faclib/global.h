#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* the following macros for Odd Even integer determination is Paltform
   dependent. this is good for UNIX systems */
#define IsOdd(x)  ((abs((x))&0x01)?1:0)
#define IsEven(x) ((abs((x))&0x01)?0:1)

#define Max(a, b) (((a)>(b))?(a):(b))
#define Min(a, b) (((a)<(b))?(a):(b))

#define EPS16 1E-16
#define EPS12 1E-12
#define EPS10 1E-10
#define EPS8  1E-08
#define EPS6  1E-06
#define EPS5  1E-05
#define EPS4  1E-04
#define EPS3  1E-03
#define EPS2  1E-02
#define EPS1  1E-01

#define PI         3.14159265359
#define TWO_PI     6.28318530718
#define SQRT2      1.414214

#define HARTREE_EV 27.2113962 
#define RATE_AU    4.13413733E16
#define RATE_AU10  4.13413733E06
#define RATE_AU12  4.13413733E04
#define AREA_AU20  2.80028560859E3
#define RBOHR      0.529177249
#define FINE_STRUCTURE_CONST  7.29735308E-3
#define FINE_STRUCTURE_CONST2 5.32513620E-5

#define G_COULOMB   1
#define G_BABUSHKIN 2

#define DEBUG_RECOUPLE  10
#define DEBUG_STRUCTURE 20
#define FAC_DEBUG 0
#if FAC_DEBUG
extern FILE *debug_log;
#endif

/*
#define PERFORM_STATISTICS 1
*/

#endif


