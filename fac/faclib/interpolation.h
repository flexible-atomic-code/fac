#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void spline(double *x, double *y, int n, 
	    double yp1, double ypn, double *y2);
int splint(double *xa, double *ya, double *y2a, 
	   int n, double x, double *y);
void splie2(double *x1a, double *x2a, double **ya, 
	    int m, int n, double **y2a);
void splin2(double *x1a, double *x2a, double **ya, double **y2a,
	   int m, int n, double x1, double x2, double *y);

#endif
