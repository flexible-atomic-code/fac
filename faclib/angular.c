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

#include "angular.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

/************************************************************
  Implementation of module "angular".

  All angular momentum arguments are twice their actual
  values to represent half-integer using integers.

  Author: M. F. Gu, mfgu@stanford.edu
*************************************************************/

double ln_factorial[MAX_FACTORIAL];
double ln_integer[MAX_FACTORIAL];

#ifdef PERFORM_STATISTICS
static ANGULAR_TIMING timing = {0, 0, 0};

/* 
** FUNCTION:    GetAngularTiming.
** PURPOSE:     Get the profiling information for 
**              module *angular*.
**              
** INPUT:       {ANGULAR_TIMING *t},
**              pointer to the struct
**              which holds the result on output.
** RETURN:      {int}, 
**              always 0.
** SIDE EFFECT: 
** NOTE:        included only if the macro PERFOR_STATISTICS 
**              is defined in "global.h".
**              the input pointer must have the storage allocated.
*/
int GetAngularTiming(ANGULAR_TIMING *t) {
  memcpy(t, &timing, sizeof(timing));
  return 0;
}
#endif

/* 
** FUNCTION:    InitAngular.
** PURPOSE:     initialize the nature log of factorial 
**              and integer arrays.
** INPUT:       
** RETURN:      {int}, 
**              always 0.
** SIDE EFFECT: 
** NOTE:        
*/
int InitAngular(void) {
  int n;

  ln_factorial[0] = 0.0;
  ln_integer[0] = -100.0;
  for (n = 1; n < MAX_FACTORIAL; n++) {
    ln_integer[n] = log((double) n);
    ln_factorial[n] = ln_factorial[n-1] + ln_integer[n]; 
  }
  return 0;
}
 
/* 
** FUNCTION:    Triangle.
** PURPOSE:     check for triangular relation.
** INPUT:       {int j1},
**              angular momentum.
**              {int j2},
**              angular momentum.
**              {int j3},
**              angular momentum.
** RETURN:      {int},
**              0: triangular relation failed.
**              1: triangular relation holds.
** SIDE EFFECT: 
** NOTE:        
*/
int Triangle(int j1, int j2, int j3) {
  int i;

  j1++;
  j2++;
  j3++;

  i = j2 - j3;
  if (j1 >= abs(i) + 1 && j1 <= j2 + j3 - 1) 
    return 1;
  else 
    return 0;
}

/* 
** calculate the Wigner 3j symbols.
** maximum summation terms of 512 should allow 
** the angular momentum up to about 500.
*/
#define MAXTERM 512
static double _sumk[MAXTERM];
#pragma omp threadprivate(_sumk)

/* 
** FUNCTION:    W3j.
** PURPOSE:     calculate the Wigner 3j symbol.
** INPUT:       {int j1},
**              angular momentum.
**              {int j2},
**              angular momentum.
**              {int j3},
**              angular momentum.
**              {int m1},
**              projection of j1.
**              {int m2},
**              projection of j2.
**              {int m3},
**              projection of j3.
** RETURN:      {double},
**              3j coefficients.
** SIDE EFFECT: 
** NOTE:        the _sumk array is used to store all 
**              summation terms to avoid overflow. 
**              the predefined MAXTERM=512 allows the 
**              maximum angular momentum of about 500.
**              if this limit is exceeded, the routine
**              issues a warning.
*/
double W3j(int j1, int j2, int j3, int m1, int m2, int m3) {
  int i, k, kmin, kmax, ik[14];
  double delta, qsum, a, b;

#ifdef PERFORM_STATISTICS
  clock_t start, stop; 
  start = clock();
#endif

  if (m1 + m2 + m3) return 0.0;
  if (!Triangle(j1, j2, j3)) return 0.0;
  if (abs(m1) > j1) return 0.0;
  if (abs(m2) > j2) return 0.0;
  if (abs(m3) > j3) return 0.0;

  ik[0] = j1 + j2 - j3;
  ik[1] = j1 - j2 + j3;
  ik[2] = -j1 + j2 + j3;
  ik[3] = j1 + j2 + j3 + 2;
  ik[4] = j1 - m1;
  ik[5] = j1 + m1;
  ik[6] = j2 - m2;
  ik[7] = j2 + m2;
  ik[8] = j3 - m3;
  ik[9] = j3 + m3;
  ik[10] = j2 - j3 - m1;
  ik[11] = j1 - j3 + m2;
  ik[12] = j3 - j2 + m1;
  ik[13] = j1 - j2 - m3;
  
  for (i = 0; i < 14; i++) {
    ik[i] = ik[i] / 2;
  }

  delta = - LnFactorial(ik[3]);
  for (i = 0; i < 3; i++) delta += LnFactorial(ik[i]);
  for (i = 4; i < 10; i++) delta += LnFactorial(ik[i]);

  kmin = Max(0, ik[10]);
  kmin = Max(kmin, ik[11]);
  kmax = Min(ik[0], ik[4]);
  kmax = Min(kmax, ik[7]);
  
  qsum = 0.0;
  a = 1E30;
  for (k = kmin, i = 0; k <= kmax && i < MAXTERM; k++, i++) {
    _sumk[i] = LnFactorial(k) +
      LnFactorial(ik[0]-k) +
      LnFactorial(ik[4]-k) +
      LnFactorial(ik[7]-k) +
      LnFactorial(ik[12]+k) +
      LnFactorial(k-ik[11]);
    if (_sumk[i] < a) a = _sumk[i];
  }
  if (i == MAXTERM) {
    printf("Maximum terms in the 3j symbol sum reached\n");
    printf("Results may be inaccurate\n");
  }
  for (k = kmin, i = 0; k <= kmax; k++, i++) {
    b = exp(-(_sumk[i]-a));    
    if (IsOdd(k)) b = -b;    
    qsum += b;
  }

  if (IsOdd(ik[13])) qsum = -qsum;
  
  b = exp(0.5*delta-a);
  b *= qsum; 

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.w3j += stop -start;
#endif

  return b;
}

/* 
** FUNCTION:    W6j.
** PURPOSE:     calculate the 6j symbol.
** INPUT:       {int j1},
**              angular momentum.
**              {int j2},
**              angular momentum.
**              {int j3},
**              angular momentum.
**              {int i1},
**              angular momentum.
**              {int i2},
**              angular momentum.
**              {int i3},
**              angular momentum.
** RETURN:      {double},
**              6j symbol.
** SIDE EFFECT: 
** NOTE:        
*/
double W6j(int j1, int j2, int j3, int i1, int i2, int i3) {
  int n1, n2, n3, n4, n5, n6, n7, k, kmin, kmax, ic, ki;
  double r, a;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  if (!(Triangle(j1, j2, j3) &&
	Triangle(j1, i2, i3) &&
	Triangle(i1, j2, i3) &&
	Triangle(i1, i2, j3)))
    return 0.0;

  n1 = (j1 + j2 + j3) / 2;
  n2 = (i2 + i1 + j3) / 2;
  n3 = (j1 + i2 + i3) / 2;
  n4 = (j2 + i1 + i3) / 2;
  n5 = (j1 + j2 + i2 + i1) / 2;
  n6 = (j1 + i1 + j3 + i3) / 2;
  n7 = (j2 + i2 + j3 + i3) / 2;
  
  kmin = Max(n1, n2);
  kmin = Max(kmin, n3);
  kmin = Max(kmin, n4) + 1;
  kmax = Min(n5, n6);
  kmax = Min(kmax, n7) + 1;

  r = 1.0;
  ic = 0;
  for (k = kmin + 1; k <= kmax; k++) {
    ki = kmax -ic;
    r = 1.0 - (r * ki * (n5-ki+2.0) * (n6-ki+2.0) * (n7-ki+2.0))/
      ((ki-1.0-n1) * (ki-1.0-n2) * (ki-1.0-n3) * (ki - 1.0 - n4));
    ic++;
  }

  a = (LnFactorial(kmin) -
       LnFactorial(kmin-n1-1) -
       LnFactorial(kmin-n2-1) -
       LnFactorial(kmin-n3-1) -
       LnFactorial(kmin-n4-1) -
       LnFactorial(n5+1-kmin) -
       LnFactorial(n6+1-kmin) -
       LnFactorial(n7+1-kmin)) +
    ((LnFactorial(n1-j1) + LnFactorial(n1-j2) +
      LnFactorial(n1-j3) - LnFactorial(n1+1) +
      LnFactorial(n2-i2) + LnFactorial(n2-i1) +
      LnFactorial(n2-j3) - LnFactorial(n2+1) +
      LnFactorial(n3-j1) + LnFactorial(n3-i2) +
      LnFactorial(n3-i3) - LnFactorial(n3+1) +
      LnFactorial(n4-j2) + LnFactorial(n4-i1) +
      LnFactorial(n4-i3) - LnFactorial(n4+1))/2.0);

  r = r * exp(a);
  
  if (IsEven(n5+kmin)) r = -r;
  if (IsOdd(((j1+j2+i1+i2)/2))) r = -r;

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.w6j += stop - start;
#endif
  return r;

}

/* 
** FUNCTION:    W6jTriangle.
** PURPOSE:     determine if 6j symbol is permitted
**              by the triangular constraints.
** INPUT:       {int j1},
**              angular momentum.
**              {int j2},
**              angular momentum.
**              {int j3},
**              angular momentum.
**              {int i1},
**              angular momentum.
**              {int i2},
**              angular momentum.
**              {int i3},
**              angular momentum.
** RETURN:      {int},
**              0: 6j symbol forbidden.
**              1: 6j symbol allowed.
** SIDE EFFECT: 
** NOTE:        
*/
int W6jTriangle(int j1, int j2, int j3, int i1, int i2, int i3) {
  return (Triangle(j1, j2, j3) &&
	  Triangle(j1, i2, i3) &&
	  Triangle(i1, j2, i3) &&
	  Triangle(i1, i2, j3));
}
  
/* 
** FUNCTION:    W9j.
** PURPOSE:     calculate the 9j symbol.
** INPUT:       {int j1},
**              angular momentum.
**              {int j2},
**              angular momentum.
**              {int j3},
**              angular momentum.
**              {int i1},
**              angular momentum.
**              {int i2},
**              angular momentum.
**              {int i3},
**              angular momentum.
**              {int k1},
**              angular momentum.
**              {int k2},
**              angular momentum.
**              {int k3},
**              angular momentum.
** RETURN:      {double},
**              9j symbol.
** SIDE EFFECT: 
** NOTE:        
*/     
double W9j(int j1, int j2, int j3,
	   int i1, int i2, int i3,
	   int k1, int k2, int k3) {
  int j, jmin, jmax;
  double r;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;  
  start = clock();
#endif

  if (!Triangle(j1, j2, j3) ||
      !Triangle(i1, i2, i3) ||
      !Triangle(k1, k2, k3) ||
      !Triangle(j1, i1, k1) ||
      !Triangle(j2, i2, k2) ||
      !Triangle(j3, i3, k3))
    return 0.0;

  jmin = Max(abs(j1-k3), abs(j2-i3));
  jmin = Max(jmin, abs(k2-i1));
  jmax = Min(j1+k3, j2+i3);
  jmax = Min(jmax, k2+i1);

  r = 0.0;
  for (j = jmin; j <= jmax; j += 2) {
    r = r + ((j+1.0) * 
	     W6j(j1, i1, k1, k2, k3, j) *
	     W6j(j2, i2, k2, i1, j, i3) *
	     W6j(j3, i3, k3, j, j1, j2));
  }

  if (IsOdd(jmin)) r = -r;

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.w9j += stop - start;
#endif

  return r;
}

/* 
** FUNCTION:    W9jTriangle.
** PURPOSE:     determine if 9j symbol is allowed by
**              the triangular constraints.
** INPUT:       {int j1},
**              angular momentum.
**              {int j2},
**              angular momentum.
**              {int j3},
**              angular momentum.
**              {int i1},
**              angular momentum.
**              {int i2},
**              angular momentum.
**              {int i3},
**              angular momentum.
**              {int k1},
**              angular momentum.
**              {int k2},
**              angular momentum.
**              {int k3},
**              angular momentum.
** RETURN:      {int},
**              0: 9j symbol forbidden.
**              1: 9j symbol allowed.
** SIDE EFFECT: 
** NOTE:        
*/     
int W9jTriangle(int j1, int j2, int j3,
		int i1, int i2, int i3,
		int k1, int k2, int k3) {
  return (Triangle(j1, j2, j3) &&
	  Triangle(i1, i2, i3) &&
	  Triangle(k1, k2, k3) &&
	  Triangle(j1, i1, k1) &&
	  Triangle(j2, i2, k2) &&
	  Triangle(j3, i3, k3));
}

/* 
** FUNCTION:    WignerEckartFactor.
** PURPOSE:     calculate the geometric prefactor in 
**              Wigner Eckart theorem, 
**              (-1)^{jf-mf}sqrt(2*jf+1)W3j(jf, k, ji, -mf, q, mi)
** INPUT:       {int jf},
**              angular momentum.
**              {int k },
**              angular momentum.
**              {int ji},
**              angular momentum.
**              {int mf},
**              projection of jf.
**              {int q },
**              projection of k.
**              {int mi},
**              projection of mi.
** RETURN:      {double},
**              prefactor.
** SIDE EFFECT: 
** NOTE:        
*/
double WignerEckartFactor(int jf, int k, int ji, int mf, int q, int mi) {
  double r;

  if (!Triangle(jf, k, ji)) return 0.0;
  if (mi + q - mf) return 0.0;

  r = sqrt(jf + 1.0);
  if (IsOdd((jf-mf)/2)) r = -r;
  r *= W3j(jf, k, ji, -mf, q, mi);
  return r;
}

/* 
** FUNCTION:    ClebschGordan.
** PURPOSE:     claculate the Clebsch Gordan coeff.
** INPUT:       {int j1},
**              angular momentum.
**              {int m1},
**              projection of j1.
**              {int j2},
**              angular momentum.
**              {int m2},
**              projection of j2.
**              {int jf},
**              angular momentum, final result 
**              of the coupling of j1 and j2.
**              {int mf},
**              projection of jf.
** RETURN:      {double},
**              CG coefficients.
** SIDE EFFECT: 
** NOTE:        
*/
double ClebschGordan(int j1, int m1, int j2, int m2, int jf, int mf) {
  double r;
  r = sqrt(jf+1.0);
  r *= W3j(j1, j2, jf, m1, m2, -mf);

  if (IsOdd((j1-j2+mf)/2)) r = -r;
  return r;
}

/* 
** FUNCTION:    ReducedCL.
** PURPOSE:     calculate the reduced matrix element
**              of the normalized spherical harmonics <ja||C^L||jb>
** INPUT:       {int ja},
**              angular momentum.
**              {int k },
**              rank of the spherical harmonics.
**              {int jb},
**              angular momentum.
** RETURN:      {double},
**              reduced matrix element.
** SIDE EFFECT: 
** NOTE:        it does not check for the triangular delta 
**              involving the orbital angular momenta.
*/
double ReducedCL(int ja, int k, int jb) {
  double r;

  r = sqrt((ja+1.0)*(jb+1.0))*W3j(ja, k, jb, 1, 0, -1);
  if (IsOdd((ja+1)/2)) r = -r;
  return r;
}

/*
** Wigner d-matrix <jm|exp(-iJ_y*a)|jn>
** j2 = j*2, m2 = m*2, n2 = n*2
*/
double WignerDMatrix(double a, int j2, int m2, int n2) {
  double b, c, ca, sa, x;
  int k, kmin, kmax;

  a *= 0.5;
  kmin = Max(0, (m2+n2)/2);
  kmax = Min((j2+m2)/2, (j2+n2)/2);
  ca = cos(a);
  sa = sin(a);
  x = 0.0;
  for (k = kmin; k <= kmax; k++) {
    b = pow(ca, (2*k-(m2+n2)/2));
    b *= pow(sa, (j2+(m2+n2)/2-2*k));
    c = LnFactorial(k);    
    c += LnFactorial((j2+m2)/2-k);
    c += LnFactorial((j2+n2)/2-k);
    c += LnFactorial(k-(m2+n2)/2);
    b /= exp(c);
    if (IsOdd(k)) b = -b;
    x += b;
  }
  c = LnFactorial((j2+m2)/2);
  c += LnFactorial((j2-m2)/2);
  c += LnFactorial((j2+n2)/2);
  c += LnFactorial((j2-n2)/2);
  c = exp(0.5*c);
  if (IsOdd((j2+m2)/2)) c = -c;
  x *= c;

  return x;
}
