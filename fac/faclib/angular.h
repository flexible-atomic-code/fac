#ifndef _ANGULAR_H_
#define _ANGULAR_H_

#include <math.h>
#include <stdio.h>
#include "config.h"
#include <time.h>

/** The tabulated ln(n!) goes up to n = MAX_FACTORIAL - 1 **/
#define MAX_FACTORIAL 1000

double ln_factorial[MAX_FACTORIAL];
double ln_integer[MAX_FACTORIAL];

#define LnFactorial(n) ln_factorial[(n)]
#define LnInteger(n) ln_integer[(n)]

typedef struct _ANGULAR_TIMING_ {
  clock_t w3j;
  clock_t w6j;
  clock_t w9j;
} ANGULAR_TIMING;

int GetAngularTiming(ANGULAR_TIMING *t);

int InitAngular();
int Triangle(int j1, int j2, int j3);
double W3j(int j1, int j2, int j3, int m1, int m2, int m3);
double W6jDelta(int j1, int j2, int j3);
double W6j(int j1, int j2, int j3, int i1, int i2, int i3);
int  W6jTriangle(int j1, int j2, int j3, int i1, int i2, int i3);
double W9j(int j1, int j2, int j3,
	   int i1, int i2, int i3,
	   int k1, int k2, int k3);
int W9jTriangle(int j1, int j2, int j3,
		int i1, int i2, int i3,
		int k1, int k2, int k3);
double WignerEckartFactor(int jf, int k, int ji,
			  int mf, int q, int mi);
double ClebschGordan(int j1, int m1, int j2, int m2, int jf, int mf);
double ReducedCL(int ja, int k, int jb);

#endif
