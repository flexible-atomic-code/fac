#include "angular.h"

static ANGULAR_TIMING timing = {0, 0, 0};
static double _sumk[500];

int GetAngularTiming(ANGULAR_TIMING *t) {
  memcpy(t, &timing, sizeof(timing));
  return 0;
}

int InitAngular() {
  int n;

  ln_factorial[0] = 0.0;
  ln_integer[0] = -100.0;
  for (n = 1; n < MAX_FACTORIAL; n++) {
    ln_integer[n] = log((double) n);
    ln_factorial[n] = ln_factorial[n-1] + ln_integer[n]; 
  }
  return 0;
}

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
  for (k = kmin, i = 0; k <= kmax; k++, i++) {
    _sumk[i] = LnFactorial(k) +
      LnFactorial(ik[0]-k) +
      LnFactorial(ik[4]-k) +
      LnFactorial(ik[7]-k) +
      LnFactorial(ik[12]+k) +
      LnFactorial(k-ik[11]);
    if (_sumk[i] < a) a = _sumk[i];
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

double W6jDelta(j1, j2, j3) {
  int i, ik[4];
  double delta;

  ik[0] = j1 + j2 - j3;
  ik[1] = j1 - j2 + j3;
  ik[2] = -j1 + j2 + j3;
  ik[3] = j1 + j2 + j3 + 2;

  for (i = 0; i < 4; i++) {
    ik[i] = ik[i] / 2;
  }
  
  delta = (LnFactorial(ik[0]) +
	   LnFactorial(ik[1]) +
	   LnFactorial(ik[2]) -
	   LnFactorial(ik[3]));

  delta = exp(0.5*delta);
  return delta;
}


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

int W6jTriangle(int j1, int j2, int j3, int i1, int i2, int i3) {
  return (Triangle(j1, j2, j3) &&
	  Triangle(j1, i2, i3) &&
	  Triangle(i1, j2, i3) &&
	  Triangle(i1, i2, j3));
}
       

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

double WignerEckartFactor(int jf, int k, int ji, int mf, int q, int mi) {
  double r;

  if (!Triangle(jf, k, ji)) return 0.0;
  if (mi + q - mf) return 0.0;

  r = sqrt(jf + 1.0);
  if (IsOdd((jf-mf)/2)) r = -r;
  r *= W3j(jf, k, ji, -mf, q, mi);
  return r;
}

double ClebschGordan(int j1, int m1, int j2, int m2, int jf, int mf) {
  double r;
  r = sqrt(jf+1.0);
  r *= W3j(j1, j2, jf, m1, m2, -mf);

  if (IsOdd((j1-j2+mf)/2)) r = -r;
  return r;
}

/** Reduced matrix element of C^L in jj coupled states.
    it does not include the delta factor involving the 
    orbital angular momenta **/
double ReducedCL(int ja, int k, int jb) {
  double r;

  r = sqrt((ja+1.0)*(jb+1.0))*W3j(ja, k, jb, 1, 0, -1);
  if (IsOdd((ja+1)/2)) r = -r;
  return r;
}

