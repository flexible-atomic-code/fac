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

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "angular.h"
#include "dbase.h"
#include "parser.h"
#include "polarization.h"
#include "interpolation.h"
#include "rates.h"
#include "cf77.h"

#define NEINT 240

static int nlevels=0;
static MLEVEL *levels=NULL;
static int ntr=0;
static MTR *tr_rates=NULL;
static int nce=0;
static MCE *ce_rates=NULL;
static int nai = 0;
static MAI *ai_rates=NULL;

static int max_iter = 128;
static double iter_accuracy = EPS6;
static int maxlevels = 0;
static int nmlevels=0;
static int nmlevels1;
static double *rmatrix=NULL;
static double *BL[MAXPOL+1];
static double PL[MAXPOL+1];
static double PL2[MAXPOL+1];

static struct {
  double energy;
  double esigma;
  double etrans;
  double density;
  int idr, ndr;
  double *pdr;
} params;

int InitPolarization(void) {
  int k;
  double pqa[MAXPOL*2+1];
  int ipqa[MAXPOL*2+1], ierr;
  double nu1, theta;
  int nudiff, mu1;

  InitDBase();
  InitAngular();
  InitRates();
  
  params.energy = 1E3;
  params.esigma = 0;
  params.etrans = 0;
  params.density = 1.0;

  params.idr = -1;
  params.ndr = 0;
  params.pdr = NULL;

  theta = acos(0.0);
  nu1 = 0;
  nudiff = MAXPOL*2;
  mu1 = 0;
  DXLEGF(nu1, nudiff, mu1, mu1, theta, 3, pqa, ipqa, &ierr);
  for (k = 0; k <= MAXPOL; k++) {
    PL[k] = pqa[k*2];
  }
  nu1 = 2;
  nudiff = MAXPOL*2-2;
  mu1 = 2;
  DXLEGF(nu1, nudiff, mu1, mu1, theta, 3, pqa, ipqa, &ierr);
  PL2[0] = 0.0;
  for (k = 1; k <= MAXPOL; k++) {
    PL2[k] = pqa[k*2-2];
  }

  return 0;
}

int SetMIteration(double a, int m) {
  iter_accuracy = a;
  if (m > 0) max_iter = m;
  
  return 0;
}

int SetMaxLevels(int m) {
  maxlevels = m;
  
  return 0;
}

int SetIDR(int idr, int ndr, double *pdr) {
  int n;

  if (params.ndr > 0) free(params.pdr);

  if (ndr > 0 && idr > 0) {
    if (idr < nlevels) {
      n = levels[idr].j/2 + 1;
      if (ndr != n) {
	printf("ndr=%d does not equal no. sublevels for idr=%d\n", ndr, idr);
	free(pdr);
	return -1;
      }
    } else {
      printf("idr=%d exceeds level table size\n", idr);
      free(pdr);
      return -1;
    }
  }
  params.idr = idr;
  params.ndr = ndr;
  params.pdr = pdr;
  
  return 0;
}

int SetEnergy(double energy, double esigma) {
  params.energy = energy;
  params.esigma = esigma;
  
  return 0;
}

int SetDensity(double eden) {
  params.density = eden;
  
  return 0;
}

int SetMLevels(char *fn, char *tfn) {
  F_HEADER fh;  
  EN_HEADER h;
  EN_RECORD r;
  TR_HEADER h1;
  TR_RECORD r1;
  TR_EXTRA r1x;
  TFILE *f;  
  int n, k, m, t, t0, p, k0;
  int m1, m2, j1, j2;
  int swp;
  double a, b, z, e;

  MemENTable(fn);

  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }

  if (fh.type != DB_EN) {
    printf("File type is not DB_EN\n");
    FCLOSE(f);
    return -1;
  }

  if (nlevels > 0) {
    for (k = 0; k <= MAXPOL; k++) {
      free(BL[k]);
    }
    for (t = 0; t < nlevels; t++) {
      free(levels[t].rtotal);
      free(levels[t].pop);
      free(levels[t].npop);
    }
    free(levels);
    nlevels = 0;
  }
  if (ntr > 0) {
    for (t = 0; t < ntr; t++) {
      free(tr_rates[t].rates);
    }
    free(tr_rates);
    ntr = 0;
  }

  while (1) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    nlevels += h.nlevels;
    FSEEK(f, h.length, SEEK_CUR);
  }
  
  if (nlevels == 0) {
    FCLOSE(f);
    return -1;
  }

  FSEEK(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);

  for (k = 0; k <= MAXPOL; k++) {
    BL[k] = (double *) malloc(sizeof(double)*nlevels);
  }
  levels = (MLEVEL *) malloc(sizeof(MLEVEL)*nlevels);
  t = 0;
  while (1) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    for (k = 0; k < h.nlevels; k++) {
      n = ReadENRecord(f, &r, swp);
      if (n == 0) break;
      levels[t].nele = h.nele;
      levels[t].j = r.j;
      levels[t].p = r.p;
      levels[t].energy = r.energy;
      m = r.j/2 + 1;
      levels[t].rtotal = malloc(sizeof(double)*m);
      levels[t].pop = malloc(sizeof(double)*m);
      levels[t].npop = malloc(sizeof(double)*m);
      levels[t].dtotal = 0.0;
      for (p = 0; p < m; p++) {
	levels[t].rtotal[p] = 0.0;
	levels[t].pop[p] = 0.0;
	levels[t].npop[p] = 0.0;
      }
      t++;
    }
  }
  FCLOSE(f);

  if (t != nlevels) {
    printf("Energy file %s corrupted\n", fn);
    return -1;
  }

  if (nmlevels > 0) {
    nmlevels = 0;
    free(rmatrix);
  }

  levels[0].ic = 0;
  for (t = 1; t < nlevels; t++) {
    levels[t].ic = levels[t-1].ic + levels[t-1].j/2 + 1;
  }
  if (maxlevels > 0 && maxlevels < nlevels) {
    nmlevels = levels[maxlevels-1].ic + levels[maxlevels-1].j/2+1 + 1;
    nmlevels1 = nmlevels-1;
  } else {
    nmlevels = levels[nlevels-1].ic + levels[nlevels-1].j/2 + 1;
    nmlevels1 = nmlevels;
  }
  rmatrix = (double *) malloc(sizeof(double)*nmlevels*(2+nmlevels));

  f = OpenFileRO(tfn, &fh, &swp);
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }

  if (fh.type != DB_TR) {
    printf("File type is not DB_TR\n");
    FCLOSE(f);
    return -1;
  }
  
  while (1) {
    n = ReadTRHeader(f, &h1, swp);
    if (n == 0) break;
    ntr += h1.ntransitions;
    FSEEK(f, h1.length, SEEK_CUR);
  }

  if (ntr == 0) {
    FCLOSE(f);
    return -1;
  }

  FSEEK(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);

  z = fh.atom;
  t0 = 0;
  if (h.nele == 1) {
    k = FindLevelByName(fn, 1, "1*1", "1s1", "1s+1(1)1");
    t = FindLevelByName(fn, 1, "2*1", "2s1", "2s+1(1)1");
    if (k >= 0 && t >= 0) {
      ntr += 1;
      tr_rates = (MTR *) malloc(sizeof(MTR)*ntr);
      tr_rates[0].lower = k;
      tr_rates[0].upper = t;
      tr_rates[0].multipole = 0;
      tr_rates[0].n = 1;
      tr_rates[0].rates = (double *) malloc(sizeof(double)*tr_rates[0].n);
      a = TwoPhotonRate(z, 0);
      tr_rates[0].rtotal = a;
      tr_rates[0].rates[0] = a;
      t0 = 1;
    } else {
      tr_rates = (MTR *) malloc(sizeof(MTR)*ntr);
    }
  } else if (h.nele == 2) {
    k = FindLevelByName(fn, 2, "1*2", "1s2", "1s+2(0)0");
    t = FindLevelByName(fn, 2, "1*1.2*1", "1s1.2s1", "1s+1(1)1.2s+1(1)0");
    if (k >= 0 && t >= 0) {
      ntr += 1;
      tr_rates = (MTR *) malloc(sizeof(MTR)*ntr);
      tr_rates[0].lower = k;
      tr_rates[0].upper = t;
      tr_rates[0].multipole = 0;
      tr_rates[0].n = 1;
      tr_rates[0].rates = (double *) malloc(sizeof(double)*tr_rates[0].n);
      a = TwoPhotonRate(z, 1);
      tr_rates[0].rtotal = a;
      tr_rates[0].rates[0] = a;
      t0 = 1;
    } else {
      tr_rates = (MTR *) malloc(sizeof(MTR)*ntr);
    }
  } else {
    tr_rates = (MTR *) malloc(sizeof(MTR)*ntr);
  }

  while (1) {
    n = ReadTRHeader(f, &h1, swp);
    if (n == 0) break;
    k0 = abs(h1.multipole)*2;
    for (t = 0; t < h1.ntransitions; t++) {
      n = ReadTRRecord(f, &r1, &r1x, swp);
      if (n == 0) break;
      j1 = levels[r1.lower].j;
      j2 = levels[r1.upper].j;
      if (k0 != 0) {
	k = k0;
	tr_rates[t0].multipole = h1.multipole;
      } else {
	k = abs(j1-j2);
	if (k == 0) k = 2;
	if (IsOdd(levels[r1.lower].p+levels[r1.upper].p+k/2)) {
	  tr_rates[t0].multipole = k/2;
	} else {
	  tr_rates[t0].multipole = -k/2;
	}
      }
      tr_rates[t0].lower = r1.lower;
      tr_rates[t0].upper = r1.upper;
      e = levels[r1.upper].energy-levels[r1.lower].energy;
      b = OscillatorStrength(h1.multipole, e, r1.strength, &a);
      a *= RATE_AU;
      tr_rates[t0].rtotal = a/(j2+1.0);
      tr_rates[t0].n = (j1/2+1)*(j2/2+1);
      tr_rates[t0].rates = (double *) malloc(sizeof(double)*tr_rates[t0].n);
      p = 0;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	for (m2 = -j2; m2 <= 0; m2 += 2) {
	  b = W3j(j1, k, j2, -m1, m1-m2, m2);
	  b = b*b;
	  tr_rates[t0].rates[p] = a*b;
	  if (m1 != 0) {
	    b = W3j(j1, k, j2, -m1, m1+m2, -m2);
	    b = b*b;
	    tr_rates[t0].rates[p] += a*b;
	  }
	  p++;
	}
      }
      t0++;
    }
  }
  if (t0 < ntr) ntr = t0;

  FCLOSE(f);

  return 0;
}  

int SetMCERates(char *fn) {
  F_HEADER fh;  
  CE_HEADER h;
  CE_RECORD r;
  TFILE *f;
  int n, k, m, t, p, i, q;
  int m1, m2, j1, j2, i0;
  int swp;
  double data[2+(1+MAXNUSR)*4];
  double cs1[128], cs2[128];
  double e1, e2, e, a, v, ratio;
  double esigma, energy;
  double egrid[NEINT], rint[NEINT], fint[NEINT];
  
  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }

  if (fh.type != DB_CE) {
    printf("File type is not DB_CE\n");
    FCLOSE(f);
    return -1;
  }

  energy = params.energy;
  esigma = params.esigma;
  
  e1 = energy;
  v = 0.0;
  if (esigma > 0) {
    a = 20.0*esigma/(NEINT-1.0);
    egrid[0] = energy - 10.0*esigma;
    for (i = 1; i < NEINT; i++) {
      egrid[i] = egrid[i-1] + a;
    }
    for (i = 0; i < NEINT; i++) {
      a = (egrid[i] - energy)/esigma;
      v = VelocityFromE(egrid[i], 0.0);
      fint[i] = (1.0/(sqrt(2.0*PI)*esigma))*exp(-0.5*a*a)*v;
    }
  } else {
    v = VelocityFromE(e1, 0.0);
  }

  if (nce > 0) {
    for (t = 0; t < nce; t++) {
      free(ce_rates[t].rates);
    }
    free(ce_rates);
    nce = 0;
  }

  while (1) {
    n = ReadCEHeader(f, &h, swp);
    if (n == 0) break;
    nce += h.ntransitions;
    FSEEK(f, h.length, SEEK_CUR);
  }
  FSEEK(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);

  ce_rates = (MCE *) malloc(sizeof(MCE)*nce);
  
  t = 0;
  while (1) {
    n = ReadCEHeader(f, &h, swp);
    if (n == 0) break;
    PrepCECrossHeader(&h, data);
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f, &r, swp, &h);
      e = levels[r.upper].energy - levels[r.lower].energy;
      e *= HARTREE_EV;
      for (k = 0; k < r.nsub; k++) {
	PrepCECrossRecord(k, &r, &h, data);
	if (esigma > 0) {
	  i0 = -1;
	  for (q = 0; q < NEINT; q++) {
	    e2 = egrid[q] - e;
	    if (e2 >= 0) {
	      rint[q] = InterpolateCECross(e2, &r, &h, data, &ratio);
	      a = egrid[q]/HARTREE_EV;
	      a *= (1.0 + 0.5*FINE_STRUCTURE_CONST2*a);
	      a = PI*AREA_AU20/(2.0*a);
	      rint[q] *= a*fint[q];
	      if (i0 < 0) i0 = q;
	    }
	  }
	  if (i0 >= 0) {
	    cs1[k] = Simpson(rint, i0, NEINT-1)*(egrid[1]-egrid[0]);
	    cs1[k] += (egrid[i0]-e)*rint[i0];
	  } else {
	    cs1[k] = 0.0;
	  }
	  for (q = 0; q < NEINT; q++) {
	    e2 = egrid[q];
	    rint[q] = InterpolateCECross(e2, &r, &h, data, &ratio);
	    a = e2/HARTREE_EV;
	    a *= (1.0 + 0.5*FINE_STRUCTURE_CONST2*a);
	    a = PI*AREA_AU20/(2.0*a);
	    rint[q] *= a*fint[q];
	  }
	  cs2[k] = Simpson(rint, 0, NEINT-1)*(egrid[1]-egrid[0]);
	} else {
	  e2 = e1 - e;
	  cs1[k] = InterpolateCECross(e2, &r, &h, data, &ratio);
	  a = e1/HARTREE_EV;
	  a *= (1.0 + 0.5*FINE_STRUCTURE_CONST2*a);
	  a = PI*AREA_AU20/(2.0*a);
	  cs1[k] *= a*v;
	  e2 = e1;
	  cs2[k] = InterpolateCECross(e2, &r, &h, data, &ratio);
	  a = e2/HARTREE_EV;
	  a *= (1.0 + 0.5*FINE_STRUCTURE_CONST2*a);
	  a = PI*AREA_AU20/(2.0*a);
	  cs2[k] *= a*v;
	}
      }
      ce_rates[t].lower = r.lower;
      ce_rates[t].upper = r.upper;
      j1 = levels[r.lower].j;
      j2 = levels[r.upper].j;
      ce_rates[t].n = 2*(j1/2+1)*(j2/2+1);
      ce_rates[t].rates = (double *) malloc(sizeof(double)*ce_rates[t].n);
      k = 0;
      p = 0;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	for (m2 = -j2; m2 <= j2; m2 += 2) {
	  if (m2 > 0) {
	    k++;
	    continue;
	  }
	  ce_rates[t].rates[p] = cs1[k];
	  if (m2 != 0) {
	    ce_rates[t].rates[p] += cs1[k-m2];
	  }
	  p++;
	  ce_rates[t].rates[p] = cs2[k];
	  if (m2 != 0 && m1 != 0) {
	    ce_rates[t].rates[p] += cs2[k-m2];
	  }
	  p++;
	  k++;
	}
      }
      t++;
      free(r.params);
      free(r.strength);
    }
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
  }
  FCLOSE(f);

  if (t < nce) nce = t;
  
  return 0;
}

int SetMAIRates(char *fn) {
  F_HEADER fh;
  AIM_HEADER h;
  AIM_RECORD r;
  TFILE *f;
  int n, k, m, t, p, i;
  int m1, m2, j1, j2;
  int swp;
  double e, v, x, pf=1.0;
  double esigma, energy;
    
  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }

  if (fh.type != DB_AIM) {
    printf("File type is not DB_AIM\n");
    FCLOSE(f);
    return -1;
  }
  
  energy = params.energy;
  esigma = params.esigma;
  if (esigma > 0) {
    pf = 1.0/(sqrt(2.0*PI)*esigma);
  }

  if (nai > 0) {
    for (t = 0; t < nai; t++) {
      free(ai_rates[t].rates);
    }
    free(ai_rates);
    nai = 0;
  }

  while (1) {
    n = ReadAIMHeader(f, &h, swp);
    if (n == 0) break;
    nai += h.ntransitions;
    FSEEK(f, h.length, SEEK_CUR);
  }
  FSEEK(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);
  
  ai_rates = (MAI *) malloc(sizeof(MAI)*nai);
  
  t = 0;
  while (1) {
    n = ReadAIMHeader(f, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadAIMRecord(f, &r, swp);
      e = levels[r.b].energy - levels[r.f].energy;
      if (e < 0) e -= h.emin;
      e *= HARTREE_EV;
      v = VelocityFromE(e, 0.0)*AREA_AU20*HARTREE_EV;
      if (esigma > 0.0) {
	x = (e - energy)/esigma;
	x = 0.5*x*x;
	if (x > 50.0) v = 0.0;
	v *= pf*exp(-x);
      }
      ai_rates[t].f = r.f;
      ai_rates[t].b = r.b;
      j1 = levels[r.b].j;
      j2 = levels[r.f].j;
      ai_rates[t].n = 2*(j1/2+1)*(j2/2+1);
      ai_rates[t].rates = (double *) malloc(sizeof(double)*ai_rates[t].n);
      k = 0;
      p = 0;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	for (m2 = -j2; m2 <= j2; m2 += 2) {
	  if (m2 > 0) {
	    k += 2;
	    continue;
	  } 
	  ai_rates[t].rates[p] = (r.rate[k])*RATE_AU;
	  if (m2 != 0) {
	    ai_rates[t].rates[p] += (r.rate[k-m2*2])*RATE_AU;
	  }
	  p++;
	  k++;
	  ai_rates[t].rates[p] = (r.rate[k])*v;
	  if (m2 != 0 && m1 != 0) {
	    ai_rates[t].rates[p] += (r.rate[k-m2*2])*v;
	  }
	  p++;
	  k++;
	}
      }
      t++;
      free(r.rate);
    }
    free(h.egrid);
  }

  FCLOSE(f);
  
  if (t < nai) nai = t;
  
  return 0;
}

static double Population(int iter) {
  int i, p, i1, i2;
  int j1, j2, m1, m2;
  int q1, q2, t, idr;
  FILE *f;
  double *b, a, c, eden;
  int *ipiv, info;
  int nmax, nm;

  eden = params.density;

  p = nmlevels*nmlevels;
  b = rmatrix + p;
  ipiv = (int *) (b+nmlevels);
  for (i = 0; i < p; i++) {
    rmatrix[i] = 0.0;
  }
  
  if (maxlevels > 0) nmax = maxlevels;
  else nmax = nlevels;

  for (i = 0; i < ntr; i++) {
    i1 = tr_rates[i].lower;
    i2 = tr_rates[i].upper;
    j1 = levels[i1].j;
    j2 = levels[i2].j;
    t = 0;
    q1 = Min(levels[i1].ic, nmlevels1);
    for (m1 = -j1; m1 <= 0; m1 += 2) {
      q2 = Min(levels[i2].ic, nmlevels1);
      for (m2 = -j2; m2 <= 0; m2 += 2) {
	a = tr_rates[i].rates[t++];
	p = q2*nmlevels+q1;
	if (q2 < nmlevels1) {
	  rmatrix[p] += a;
	  q2++;
	} else {
	  rmatrix[p] += levels[i2].pop[(m2+j2)/2]*a;
	}
	if (iter == 0) {
	  levels[i2].rtotal[(m2+j2)/2] += a;
	}
      }
      if (q1 < nmlevels1) {
	q1++;
      }
    }
  }
  
  for (i = 0; i < nmlevels; i++) {
    ipiv[i] = 0;
  }
  idr = -1;
  for (i = 0; i < nai; i++) {
    i1 = ai_rates[i].b;
    i2 = ai_rates[i].f;
    j1 = levels[i1].j;
    j2 = levels[i2].j;
    if (params.idr >= 0) {
      q2 = levels[i2].ic;
      for (m2 = -j2; m2 <= 0; m2 += 2) {
	ipiv[q2++] = 1;
      }
      if (params.idr == i2) idr = params.idr;
    }
    t = 0;
    q1 = Min(levels[i1].ic, nmlevels1);
    for (m1 = -j1; m1 <= 0; m1 += 2) {
      q2 = Min(levels[i2].ic, nmlevels1);
      for (m2 = -j2; m2 <= 0; m2 += 2) {
	p = q1*nmlevels+q2;
	a = ai_rates[i].rates[t++];
	if (q1 < nmlevels1) {
	  rmatrix[p] += a;
	} else {
	  rmatrix[p] += levels[i1].pop[(m1+j1)/2]*a;
	}
	if (iter == 0) {
	  levels[i1].rtotal[(m1+j1)/2] += a;
	}
	p = q2*nmlevels+q1;
	a = eden*ai_rates[i].rates[t++];
	if (q2 < nmlevels1) {
	  rmatrix[p] += a;
	} else {
	  rmatrix[p] += levels[i2].pop[(m2+j2)/2]*a;
	}
	if (iter == 0) {
	  levels[i2].rtotal[(m2+j2)/2] += a;
	}
	if (q2 < nmlevels1) {
	  q2++;
	}
      }
      if (q1 < nmlevels1) {
	q1++;
      }
    }
  }
  
  if (params.idr >= 0 && idr < 0) {
    printf("idr=%d/%d is not a DR target level\n", params.idr, idr);
    return -1;
  }
  
  if (idr < 0) {
    for (i = 0; i < nce; i++) {
      i1 = ce_rates[i].lower;
      i2 = ce_rates[i].upper;
      j1 = levels[i1].j;
      j2 = levels[i2].j;
      t = 0;
      q1 = Min(levels[i1].ic, nmlevels1);
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	q2 = Min(levels[i2].ic, nmlevels1);
	for (m2 = -j2; m2 <= 0; m2 += 2) {
	  p = q1*nmlevels+q2;
	  a = eden*ce_rates[i].rates[t++];
	  if (q1 < nmlevels1) {
	    rmatrix[p] += a;
	  } else {
	    rmatrix[p] += levels[i1].pop[(m1+j1)/2]*a;
	  }
	  if (iter == 0) {
	    levels[i1].rtotal[(m1+j1)/2] += a;
	  }
	  p = q2*nmlevels+q1;
	  a = eden*ce_rates[i].rates[t++];
	  if (q2 < nmlevels1) {
	    rmatrix[p] += a;
	  } else {
	    rmatrix[p] += levels[i2].pop[(m2+j2)/2]*a;
	  }
	  if (iter == 0) {
	    levels[i2].rtotal[(m2+j2)/2] += a;
	  }
	  if (q2 < nmlevels1) {
	    q2++;
	  }
	}
	if (q1 < nmlevels1) {
	  q1++;
	}
      }
    }  
  }

  if (nlevels > nmax) {
    a = 0;
    for (t = nmax; t < nlevels; t++) {
      a += levels[t].dtotal;
    }
    if (a) {
      for (q2 = 0; q2 < nmlevels; q2++) {
	p = nmlevels1*nmlevels+q2;
	rmatrix[p] /= a;
      }
    }
  }

  for (q1 = 0; q1 < nmlevels; q1++) {
    p = q1*nmlevels + q1;
    rmatrix[p] = 0.0;
    for (q2 = 0; q2 < nmlevels; q2++) {
      if (q2 != q1) {
	t = q1*nmlevels + q2;
	rmatrix[p] -= rmatrix[t];
      }
    }
  }

  if (idr < 0) {
    i = -1;
    for (q1 = 0; q1 < nmlevels; q1++) {
      p = q1*nmlevels + q1;
      if (rmatrix[p]) {
	i = q1;
	break;
      }
    }

    for (q1 = 0; q1 < nmlevels; q1++) {
      p = q1*nmlevels+i;
      rmatrix[p] = 1.0;
      b[q1] = 0.0;
      t = q1*nmlevels+q1;
      if (!rmatrix[t]) {
	for (q2 = 0; q2 < nmlevels; q2++) {
	  p = q2*nmlevels+q1;
	  rmatrix[p] = 0.0;
	}
	rmatrix[t] = 1E50;
      }
    }
    b[i] = 1.0;
  } else {
    q1 = levels[idr].ic;
    q2 = levels[idr].ic + (levels[idr].j/2 + 1);
    for (t = 0; t < nmlevels; t++) {
      if (ipiv[t] && !(t >= q1 && t < q2)) {
	p = t*nmlevels+t;
	rmatrix[p] = 0.0;
      }
    }

    for (q1 = 0; q1 < nmlevels; q1++) {
      t = q1*nmlevels+q1;
      if (!rmatrix[t]) {
	for (q2 = 0; q2 < nmlevels; q2++) {
	  p = q2*nmlevels+q1;
	  rmatrix[p] = 0.0;
	}
	rmatrix[t] = 1E50;
      }
    }

    j2 = levels[idr].j;
    t = j2/2 + 1;
    a = 1.0/(j2 + 1.0);
    q1 = levels[idr].ic;
    i = 0;
    for (m2 = -j2; m2 <= 0; m2 += 2) {
      for (q2 = 0; q2 < nmlevels; q2++) {
	p = q2*nmlevels+q1;
	rmatrix[p] = 0.0;
      }
      rmatrix[q1*nmlevels+q1] = 1.0;
      if (params.ndr != t) {
	if (m2 != 0) {
	  b[q1] = a*2.0;
	} else {
	  b[q1] = a;
	}
      } else {
	if (m2 != 0) {
	  b[q1] = params.pdr[i]*2.0;
	} else {
	  b[q1] = params.pdr[i];
	}
      }
      q1++;
      i++;
    }
  }

  DGESV(nmlevels, 1, rmatrix, nmlevels, ipiv, b, nmlevels, &info);
  c = 0.0;
  if (nlevels > nmax) {
    nm = 0;
    for (i = 0; i < nmax; i++) {
      p = levels[i].ic;
      j1 = levels[i].j;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	t = (m1+j1)/2;
	if (b[p]) {
	  c += fabs((b[p]-levels[i].pop[t])/b[p]);
	  nm++;
	}
	levels[i].pop[t] = b[p];
	p++;
      }
    }
  
    for (i = nmax; i < nlevels; i++) {
      j1 = levels[i].j;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	t = (m1+j1)/2;
	levels[i].npop[t] = 0.0;
      }
    }

    for (i = ntr-1; i >= 0; i--) {
      i1 = tr_rates[i].lower;
      if (i1 < nmax) continue;
      i2 = tr_rates[i].upper;
      j1 = levels[i1].j;
      j2 = levels[i2].j;
      t = 0;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	for (m2 = -j2; m2 <= 0; m2 += 2) {
	  a = levels[i2].pop[(m2+j2)/2]*tr_rates[i].rates[t++];
	  levels[i1].npop[(m1+j1)/2] += a;
	}
      }
    }
    
    for (i = nai-1; i >= 0; i--) {
      i1 = ai_rates[i].b;
      i2 = ai_rates[i].f;
      if (i1 < nmax && i2 < nmax) continue;
      j1 = levels[i1].j;
      j2 = levels[i2].j;
      t = 0;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	for (m2 = -j2; m2 <= 0; m2 += 2) {
	  if (i2 >= nmax) {
	    a = levels[i1].pop[(m1+j1)/2]*ai_rates[i].rates[t];
	    levels[i2].npop[(m2+j2)/2] += a;
	  }
	  t++;
	  if (i1 >= nmax) {
	    a = levels[i2].pop[(m2+j2)/2]*eden*ai_rates[i].rates[t];
	    levels[i1].npop[(m1+j1)/2] += a;
	  }
	  t++;
	}
      }
    }

    for (i = nce-1; i >= 0; i--) {
      i1 = ce_rates[i].lower;
      i2 = ce_rates[i].upper;
      if (i1 < nmax && i2 < nmax) continue;
      j1 = levels[i1].j;
      j2 = levels[i2].j;
      t = 0;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	for (m2 = -j2; m2 <= 0; m2 += 2) {
	  if (i2 >= nmax) {
	    a = levels[i1].pop[(m1+j1)/2]*eden*ce_rates[i].rates[t];
	    levels[i2].npop[(m2+j2)/2] += a;
	  }
	  t++;
	  if (i1 >= nmax) {
	    a = levels[i2].pop[(m2+j2)/2]*eden*ce_rates[i].rates[t];
	    levels[i1].npop[(m1+j1)/2] += a;
	  }
	  t++;
	}
      }
    }

    for (i = nmax; i < nlevels; i++) {
      j1 = levels[i].j;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	t = (m1+j1)/2;
	if (levels[i].rtotal[t]) {
	  levels[i].npop[t] = levels[i].npop[t]/levels[i].rtotal[t];
	} else {
	  levels[i].npop[t] = 0.0;
	}
	if (levels[i].npop[t]) {
	  c += fabs((levels[i].npop[t]-levels[i].pop[t])/levels[i].npop[t]);
	  nm++;
	}
	levels[i].pop[t] = levels[i].npop[t];
      }
    }
    c /= nm;
  } else {
    for (i = 0; i < nlevels; i++) {
      p = levels[i].ic;
      j1 = levels[i].j;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	t = (m1+j1)/2;
	levels[i].pop[t] = b[p];
	p++;
      }
    }
  }

  for (i = 0; i < nlevels; i++) {
    p = levels[i].ic;
    j1 = levels[i].j;
    a = 0.0;
    for (m1 = -j1; m1 <= 0; m1 += 2) {
      t = (m1+j1)/2;
      a += levels[i].pop[t];
    }
    levels[i].dtotal = a;
  }
    
  return c;
}

int PopulationTable(char *fn) {
  int i, t, j1, m1, p;
  FILE *f;
  double c;

  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }

  for (i = 0; i < max_iter; i++) {
    c = Population(i);
    printf("%5d %11.4E\n", i, c);
    fflush(stdout);
    if (c < iter_accuracy) break;
  }

  fprintf(f, "# FAC %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
  fprintf(f, "# Energy  = %-12.5E\n", params.energy);
  fprintf(f, "# ESigma  = %-12.5E\n", params.esigma);
  fprintf(f, "# Density = %-12.5E\n", params.density);
  fprintf(f, "# IDR     = %-3d\n", params.idr);
  fprintf(f, "\n");
  for (i = 0; i < nlevels; i++) {
    j1 = levels[i].j;
    fprintf(f, "%7d\t%12.5E\n", i, levels[i].dtotal);
    p = levels[i].ic;
    for (m1 = -j1; m1 <= 0; m1 += 2) {
      t = (m1+j1)/2;
      if (levels[i].dtotal) {
	c = levels[i].pop[t]/levels[i].dtotal;
	if (m1 != 0) c = 0.5*c;
      } else {
	c = 0.0;
      }
      fprintf(f, "%3d %3d\t%12.5E\n", p, -m1, c);
      p++;
    }
    fprintf(f, "\n");
  }

  fclose(f);
  return 0;
}

int Orientation(char *fn, double etrans) {
  int k, i, k2, t, j1, m1;
  double a, b;
  double pqa[MAXPOL*2+1];
  int ipqa[MAXPOL*2+1], ierr;
  double nu1, x, theta;
  int nudiff, mu1;
  FILE *f;
  
  if (fn) {
    f = fopen(fn, "w");
    if (f == NULL) {
      printf("cannot open file %s\n", fn);
      return -1;
    }
  } else {
    f = NULL;
  }
  params.etrans = etrans;
  x = sqrt(etrans/params.energy); 
  if (x >= 1.0) {
    printf("ETrans (%11.4E) >= Energy (%11.4E)\n", etrans, params.energy);
    x = 1.0-EPS10;
  }
  if (x > EPS10) {
    theta = asin(x);  
    nu1 = 0;
    nudiff = MAXPOL*2;
    mu1 = 0;
    DXLEGF(nu1, nudiff, mu1, mu1, theta, 3, pqa, ipqa, &ierr);
  }

  for (i = 0; i < nlevels; i++) {
    j1 = levels[i].j;
    for (k = 0; k <= MAXPOL; k++) {
      k2 = k*4;
      BL[k][i] = 0.0;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	t = (m1+j1)/2;
	if (levels[i].dtotal) {
	  b = levels[i].pop[t]/levels[i].dtotal;
	} else {
	  b = 0.0;
	}
	a = W3j(j1, j1, k2, -m1, m1, 0)*b;
	a *= sqrt(k2+1.0);
	if (IsOdd((j1+m1)/2)) a = -a;
	BL[k][i] += a;
      }
      BL[k][i] *= sqrt(j1 + 1.0);
      if (x > EPS10) {
	BL[k][i] *= pqa[k*2];
      }
    }
  }    

  if (f) {
    fprintf(f, "# FAC %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
    fprintf(f, "# Energy  = %-12.5E\n", params.energy);
    fprintf(f, "# ESigma  = %-12.5E\n", params.esigma);
    fprintf(f, "# ETrans  = %-12.5E\n", params.etrans);
    fprintf(f, "# Density = %-12.5E\n", params.density);
    fprintf(f, "# IDR     = %-3d\n", params.idr);
    fprintf(f, "\n");
    for (i = 0; i < nlevels; i++) {
      fprintf(f, "%5d ", i);
      for (k = 0; k <= MAXPOL; k++) {
	fprintf(f, " %10.3E", BL[k][i]);
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }
  
  return 0;
}

static int InBetween(int k, int *t) {
  int p1, p2;

  p1 = 1;
  p2 = 1;
  if (t[0] >= 0 && k < t[0]) p1 = 0;
  if (t[1] >= 0 && k > t[1]) p2 = 0;

  return p1 && p2;
}

static int InTrans(int n, int *trans, int i1, int i2, int k) {
  int i, p1, p2, p3, p4, *t;

  t = trans;
  for (i = 0; i < n; i += 7) {
    p1 = InBetween(levels[i1].nele, t);
    t += 2;
    p2 = InBetween(i1, t);
    t += 2;
    p3 = InBetween(i2, t);
    t += 2;
    if (*t == 0 || *t == k) p4 = 1;
    else p4 = 0;
    t += 1;
    if (p1 && p2 && p3 && p4) return 1;
  }

  return 0;
}
    
static int GetTrans(int *trans, char *buf) {
  char *s, *p;
  int k, t, i;

  s = buf;
  while (*s) {
    if (*s == '\t') *s = ' ';
    s++;
  }

  k = StrSplit(buf, ' ');
  if (k != 4) return -1;

  if (k == 4) {
    s = buf;
    i = 0;
    for (t = 0; t < k; t++) {
      while (*s == ' ') s++;
      if (t == 3) {
	if (*s == '*') {
	  trans[i] = 0;
	} else {
	  trans[i] = atoi(s);
	} 
      } else {
	p = s;
	while (*p && *p != '-') p++;
	if (*p) {
	  *p = '\0';
	  if (*s == '*') trans[i] = -1;
	  else trans[i] = atoi(s);
	  if (*(p+1) == '*') trans[i+1] = -1;
	  else trans[i+1] = atoi(p+1);
	  *p = '-';
	} else {
	  if (*s == '*') trans[i] = -1;
	  else trans[i] = atoi(s);
	  trans[i+1] = trans[i];
	}
      }
      i += 2;
      while (*s) s++;
      s++;
    }
  }

  return 0;
}

int PolarizationTable(char *fn, char *ifn, int n, char **sc) {
  int i, k, k2, t, t2;
  int j1, j2, m1, i1, i2;
  FILE *f;
  double AL[MAXPOL+1];
  double FL[MAXPOL+1];
  double a, b, tem, e;
  char buf[128], *s;
  int *trans, *tp, ns;

  ns = 0;
  if (ifn) {
    n = 0;
    f = fopen(ifn, "r");
    if (f == NULL) {
      printf("cannot open file %s\n", ifn);
      return -1;
    }
    while (1) {
      if (fgets(buf, 128, f) == NULL) break;
      s = buf;
      while (*s) {
	if (*s == '\t') *s = ' ';
	s++;
      }
      k = StrSplit(buf, ' ');
      if (k == 4) {
	n++;
      }
    }
    if (n > 0) {
      trans = malloc(sizeof(int)*n*7);
      fseek(f, 0, SEEK_SET);
      tp = trans;
      while (1) {
	if (fgets(buf, 128, f) == NULL) break;
	k = GetTrans(tp, buf);
	if (k == 0) {
	  tp += 7;
	  ns++;
	}
      }
    }
    fclose(f);
  } else if (n > 0) {
    trans = malloc(sizeof(int)*n*7);
    tp = trans;
    for (i = 0; i < n; i++) {
      k = GetTrans(tp, sc[i]);
      if (k == 0) {
	tp += 7;
	ns++;
      }
    }
  }

  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  
  fprintf(f, "# FAC %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
  fprintf(f, "# Energy  = %-12.5E\n", params.energy);
  fprintf(f, "# ESigma  = %-12.5E\n", params.esigma);
  fprintf(f, "# ETrans  = %-12.5E\n", params.etrans);
  fprintf(f, "# Density = %-12.5E\n", params.density);
  fprintf(f, "# IDR     = %-3d\n", params.idr);
  fprintf(f, "\n");
  for (i = 0; i < ntr; i++) {
    i1 = tr_rates[i].upper;
    i2 = tr_rates[i].lower;
    k = tr_rates[i].multipole;
    if (k == 0) continue;
    if (ns > 0) {
      if (!(InTrans(ns, trans, i1, i2, k))) continue;
    }
    k2 = 2*abs(k);
    j1 = levels[i1].j;
    j2 = levels[i2].j;
    for (t = 0; t <= MAXPOL; t++) {
      t2 = 4*t;
      b = W3j(k2, k2, t2, 2, -2, 0);
      a = b*W6j(k2, k2, t2, j1, j1, j2);
      if (a) {
	a *= k2+1.0;
	a *= sqrt(j1+1.0);
	a *= sqrt(t2+1.0);
	if (IsEven((j1+j2)/2)) a = -a;
      }
      AL[t] = a;
      FL[t] = 0;
      if (t > 0) {
	if (b) {
	  a = W3j(k2, k2, t2, 2, 2, -4);
	  a = a/b;
	  a *= exp(0.5*(LnFactorial(2*t-2)-LnFactorial(2*t+2)));
	  if (k < 0) a = -a;
	  FL[t] = a;
	}
      }
    }
    a = 0.0;
    b = 0.0;
    for (t = 0; t <= MAXPOL; t++) {
      a += FL[t]*PL2[t]*AL[t]*BL[t][i1];
      b += PL[t]*AL[t]*BL[t][i1];
    }
    if (a) {
      a = a/b;
    }
    tem = levels[i1].dtotal*tr_rates[i].rtotal;
    e = (levels[i1].energy - levels[i2].energy)*HARTREE_EV;
    fprintf(f, "%3d %5d %5d %2d %12.5E %12.5E %10.3E %10.3E\n",
	    levels[i1].nele, i1, i2, k, e, tem, b, a);
  }

  if (n > 0) free(trans);
  fclose(f);
  return 0;
}
