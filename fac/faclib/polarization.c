static char *rcsid="$Id: polarization.c,v 1.6 2003/07/15 17:59:20 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "angular.h"
#include "dbase.h"
#include "polarization.h"
#include "rates.h"
#include "cf77.h"

static int nlevels=0;
static MLEVEL *levels=NULL;
static int ntr=0;
static MTR *tr_rates=NULL;
static int nce=0;
static MCE *ce_rates=NULL;

static int nmlevels=0;
static double *rmatrix=NULL;

static struct {
  double energy;
  double density;
} params;

int InitPolarization(void) {
  InitDBase();
  InitAngular();
  InitRates();

  return 0;
}

int SetMLevels(char *fn, char *tfn) {
  F_HEADER fh;  
  EN_HEADER h;
  EN_RECORD r;
  TR_HEADER h1;
  TR_RECORD r1;
  FILE *f;  
  int n, k, m, t, t0, p;
  int m1, m2, j1, j2;
  int swp;
  double a, b, z, e;

  MemENTable(fn);

  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  
  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) {
    fclose(f);
    return 0;
  }

  if (fh.type != DB_EN) {
    printf("File type is not DB_EN\n");
    fclose(f);
    return -1;
  }

  if (nlevels > 0) {
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
    fseek(f, h.length, SEEK_CUR);
  }
  
  if (nlevels == 0) {
    fclose(f);
    return -1;
  }

  fseek(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);

  levels = (MLEVEL *) malloc(sizeof(MLEVEL)*nlevels);
  t = 0;
  while (1) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    for (k = 0; k < h.nlevels; k++) {
      n = ReadENRecord(f, &r, swp);
      if (n == 0) break;
      levels[t].p = r.p;
      levels[t].j = r.j;
      levels[t].energy = r.energy;
      t++;
    }
  }
  fclose(f);

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
    levels[t].ic = levels[t-1].ic + levels[t-1].j+1;
  }
  nmlevels = levels[nlevels-1].ic + levels[nlevels-1].j+1;
  rmatrix = (double *) malloc(sizeof(double)*nmlevels*(2+nmlevels));

  f = fopen(tfn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  
  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) {
    fclose(f);
    return 0;
  }

  if (fh.type != DB_TR) {
    printf("File type is not DB_TR\n");
    fclose(f);
    return -1;
  }
  
  while (1) {
    n = ReadTRHeader(f, &h1, swp);
    if (n == 0) break;
    ntr += h1.ntransitions;
    fseek(f, h1.length, SEEK_CUR);
  }

  if (ntr == 0) {
    fclose(f);
    return -1;
  }

  fseek(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);

  z = fh.atom;
  t0 = 0;
  if (h.nele == 1) {
    k = FindLevelByName(fn, 1, "1*1", "1s1", "1s+1(1)1");
    t = FindLevelByName(fn, 1, "2*1", "2s1", "1s+1(1)1");
    if (k >= 0 && t >= 0) {
      ntr += 1;
      tr_rates = (MTR *) malloc(sizeof(MTR)*ntr);
      tr_rates[0].lower = k;
      tr_rates[0].upper = t;
      tr_rates[0].multipole = 0;
      tr_rates[0].n = (levels[k].j+1)*(levels[t].j+1);
      tr_rates[0].rates = (double *) malloc(sizeof(tr_rates[0].n));
      a = TwoPhotonRate(z, 0);
      tr_rates[0].rtotal = a;
      a /= levels[k].j+1.0;
      for (m = 0; m < tr_rates[0].n; m++) {
	tr_rates[0].rates[m] = a;
      }
      t0 = 1;
    } else {
      tr_rates = (MTR *) malloc(sizeof(MTR)*ntr);
    }
  } else if (h.nele == 2) {
    k = FindLevelByName(fn, 2, "1*2", "1s2", "1s+2(0)0");
    t = FindLevelByName(fn, 2, "1*1 2*1", "1s1 2s1", "1s+1(1)1 2s+1(1)0");
    if (k >= 0 && t >= 0) {
      ntr += 1;
      tr_rates = (MTR *) malloc(sizeof(MTR)*ntr);
      tr_rates[0].lower = k;
      tr_rates[0].upper = t;
      tr_rates[0].multipole = 0;
      tr_rates[0].n = (levels[k].j+1)*(levels[t].j+1);
      tr_rates[0].rates = (double *) malloc(sizeof(tr_rates[0].n));
      a = TwoPhotonRate(z, 1);
      tr_rates[0].rtotal = a;
      a /= (levels[k].j+1.0);
      for (m = 0; m < tr_rates[0].n; m++) {
	tr_rates[0].rates[m] = a;
      }
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
    k = abs(h1.multipole)*2;
    for (t = 0; t < h1.ntransitions; t++) {
      n = ReadTRRecord(f, &r1, swp);
      if (n == 0) break;
      tr_rates[t0].multipole = h1.multipole;
      tr_rates[t0].lower = r1.lower;
      tr_rates[t0].upper = r1.upper;
      j1 = levels[r1.lower].j;
      j2 = levels[r1.upper].j;
      e = levels[r1.upper].energy-levels[r1.lower].energy;
      a = 2.0*pow((FINE_STRUCTURE_CONST*e),2)*FINE_STRUCTURE_CONST;
      a *= r1.strength;
      a *= RATE_AU;
      tr_rates[t0].rtotal = a/(j2+1.0);
      tr_rates[t0].n = (j1+1)*(j2+1);
      tr_rates[t0].rates = (double *) malloc(sizeof(double)*tr_rates[t0].n);
      p = 0;
      for (m1 = -j1; m1 <= j1; m1 += 2) {
	for (m2 = -j2; m2 <= j2; m2 += 2) {
	  b = W3j(j1, k, j2, -m1, m1-m2, m2);
	  tr_rates[t0].rates[p] = a*b*b;
	  p++;
	}
      }
      t0++;
    }
  }
  if (t0 < ntr) ntr = t0;

  fclose(f);

  return 0;
}  

int SetMCERates(char *fn, double energy) {  
  F_HEADER fh;  
  CE_HEADER h;
  CE_RECORD r;
  FILE *f;
  int n, k, m, t, p, i;
  int m1, m2, j1, j2;
  int swp;
  double data[2+(1+MAXNUSR)*4];
  double cs1[128], cs2[128];
  double e1, e2, e, a, v, ratio;

  
  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  
  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) {
    fclose(f);
    return 0;
  }

  if (fh.type != DB_CE) {
    printf("File type is not DB_CE\n");
    fclose(f);
    return -1;
  }

  params.energy = energy;

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
    fseek(f, h.length, SEEK_CUR);
  }
  fseek(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);

  ce_rates = (MCE *) malloc(sizeof(MCE)*nce);
  
  t = 0;
  e1 = energy;
  v = VelocityFromE(e1);
  while (1) {
    n = ReadCEHeader(f, &h, swp);
    if (n == 0) break;
    PrepCECrossHeader(&h, data);
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f, &r, swp, &h);
      e = levels[r.upper].energy - levels[r.lower].energy;
      e *= HARTREE_EV;
      e2 = e1 + e;
      for (k = 0; k < r.nsub; k++) {
	PrepCECrossRecord(k, &r, &h, data);
	cs1[k] = InterpolateCECross(e1, &r, &h, data, &ratio);
	a = e1/HARTREE_EV;
	a = a*(1.0+0.5*FINE_STRUCTURE_CONST2*a);
	a = PI*AREA_AU20/(2.0*a);
	cs1[k] *= a*v;
	cs2[k] = InterpolateCECross(e2, &r, &h, data, &ratio);
	a = e2/HARTREE_EV;
	a = a*(1.0+0.5*FINE_STRUCTURE_CONST2*a);
	a = PI*AREA_AU20/(2.0*e2/HARTREE_EV);
	cs2[k] *= a*v;	
      }
      ce_rates[t].lower = r.lower;
      ce_rates[t].upper = r.upper;
      j1 = levels[r.lower].j;
      j2 = levels[r.upper].j;
      ce_rates[t].n = 2*(j1+1)*(j2+1);
      ce_rates[t].rates = (double *) malloc(sizeof(double)*ce_rates[t].n);
      k = 0;
      p = 0;
      for (m1 = -j1; m1 <= j1; m1 += 2) {
	for (m2 = -j2; m2 <= j2; m2 += 2) {
	  ce_rates[t].rates[p++] = cs1[k];
	  ce_rates[t].rates[p++] = cs2[k];
	  if (m1 <= 0) {
	    k++;
	  } else {
	    k--;
	  }
	}
	if (m1 == 0) {
	  k -= j2+2;
	} else if (m1 == 1) {
	  k--;
	}
      }
      t++;
    }
  }
  fclose(f);

  if (t < nce) nce = t;
  
  return 0;
}

int PopulationTable(char *fn, double eden) {
  int i, p, i1, i2;
  int j1, j2, m1, m2;
  int q1, q2, t;
  FILE *f;
  double *b, a, c;
  int *ipiv, info;

  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }

  params.density = eden;

  p = nmlevels*nmlevels;
  b = rmatrix + p;
  ipiv = (int *) (b+nmlevels);
  for (i = 0; i < p; i++) {
    rmatrix[i] = 0.0;
  }

  for (i = 0; i < ntr; i++) {
    i1 = tr_rates[i].lower;
    i2 = tr_rates[i].upper;
    j1 = levels[i1].j;
    j2 = levels[i2].j;
    t = 0;
    q1 = levels[i1].ic;
    for (m1 = -j1; m1 <= j1; m1 += 2) {
      q2 = levels[i2].ic;
      for (m2 = -j2; m2 <= j2; m2 += 2) {
	p = q2*nmlevels+q1;
	rmatrix[p] += tr_rates[i].rates[t++];
	q2++;
      }
      q1++;
    }
  }

  for (i = 0; i < nce; i++) {
    i1 = ce_rates[i].lower;
    i2 = ce_rates[i].upper;
    j1 = levels[i1].j;
    j2 = levels[i2].j;
    t = 0;
    q1 = levels[i1].ic;
    for (m1 = -j1; m1 <= j1; m1 += 2) {
      q2 = levels[i2].ic;
      for (m2 = -j2; m2 <= j2; m2 += 2) {
	p = q1*nmlevels+q2;
	a = eden*ce_rates[i].rates[t++];
	rmatrix[p] += a;
	p = q2*nmlevels+q1;
	a = eden*ce_rates[i].rates[t++];
	rmatrix[p] += a;
	q2++;
      }
      q1++;
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

  for (q1 = 0; q1 < nmlevels; q1++) {
    p = q1*nmlevels;
    rmatrix[p] = 1.0;
    b[q1] = 0.0;
  }
  b[0] = 1.0;

  DGESV(nmlevels, 1, rmatrix, nmlevels, ipiv, b, nmlevels, &info);

  fprintf(f, "Energy  = %12.5E\n", params.energy);
  fprintf(f, "Density = %12.5E\n", params.density);
  fprintf(f, "\n");
  for (i = 0; i < nlevels; i++) {
    j1 = levels[i].j;
    a = 0.0;
    p = levels[i].ic;
    for (m1 = -j1; m1 <= j1; m1 += 2) {
      a += b[p];
      p++;
    }
    levels[i].dtotal = a;
    p = levels[i].ic;
    fprintf(f, "%5d\t%12.5E\n", i, a);
    for (m1 = -j1; m1 <= j1; m1 += 2) {
      if (a) {
	c = b[p]/a;
      } else {
	c = 0.0;
      }
      fprintf(f, "%5d\t%12.5E\n", m1, c);
      p++;
    }
    fprintf(f, "\n");
  }

  fclose(f);
  return 0;
}

int PolarizationTable(char *fn) {
  int i, k, k2, t, t2;
  int j1, j2, m1, i1, i2;
  FILE *f;
  double *BL[MAXPOL+1];
  double AL[MAXPOL+1];
  double FL[MAXPOL+1];
  double PL[MAXPOL+1];
  double PL2[MAXPOL+1];
  double pqa[MAXPOL*2+1];
  int ipqa[MAXPOL*2+1], ierr;
  double nu1, theta;
  int nudiff, mu1;
  double a, b, tem, e, *x;

  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  
  x = rmatrix+nmlevels*nmlevels;

  for (k = 0; k <= MAXPOL; k++) {
    BL[k] = (double *) malloc(sizeof(double)*nlevels);
  }
  for (i = 0; i < nlevels; i++) {
    j1 = levels[i].j;
    for (k = 0; k <= MAXPOL; k++) {
      k2 = k*4;
      t = levels[i].ic;    
      BL[k][i] = 0.0;
      for (m1 = -j1; m1 <= j1; m1 += 2) {
	if (levels[i].dtotal) {
	  b = x[t]/levels[i].dtotal;
	} else {
	  b = 0.0;
	}
	a = W3j(j1, j1, k2, -m1, m1, 0)*b;
	a *= sqrt(k2+1.0);
	if (IsOdd((j1+m1)/2)) a = -a;
	BL[k][i] += a;
	t++;
      }
      BL[k][i] *= sqrt(j1 + 1.0);
    }
  }

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
  
  fprintf(f, "Energy  = %12.5E\n", params.energy);
  fprintf(f, "Density = %12.5E\n", params.density);
  fprintf(f, "\n");
  for (i = 0; i < ntr; i++) {
    k = tr_rates[i].multipole;
    if (k == 0) continue;
    k2 = 2*abs(k);
    i1 = tr_rates[i].upper;
    i2 = tr_rates[i].lower;
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
    tem = levels[i1].dtotal*tr_rates[i].rtotal/params.density;
    e = (levels[i1].energy - levels[i2].energy)*HARTREE_EV;
    fprintf(f, "%5d %5d %2d %12.5E %12.5E %10.3E %10.3E\n",
	    i1, i2, k, e, tem, b, a);
  }
  
  for (k = 0; k <= MAXPOL; k++) {
    free(BL[k]);
  }
  fclose(f);
  return 0;
}
