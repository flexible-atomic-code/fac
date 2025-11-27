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

#include "stoken.h"
#include "crm.h"

static int PPrint(int argc, char *argv[], int argt[], ARRAY *variables) {
  int i;
  for (i = 0; i < argc; i++) {
    switch (argt[i]) {
    case NUMBER:
      printf("%s", argv[i]);
      break;
    case STRING:
      printf("%s", argv[i]);
      break;
    case LIST:
      printf("[%s]", argv[i]);
      break;
    case TUPLE:
      printf("(%s)", argv[i]);
      break;
    case KEYWORD:
      printf("%s = ", argv[i]);
    }
    if (i != argc-1 && argt[i] != KEYWORD) {
      printf(" ");
    }
  }
  if (argc > 0) printf("\n");
  fflush(stdout);
  return 0;
}

static int PExit(int argc, char *argv[], int argt[], ARRAY *variables) {
  if (argc != 0) return -1;
  exit(0);
}

static int PCheckEndian(int argc, char *argv[], int argt[], ARRAY *variables) {
  TFILE *f;
  F_HEADER fh;
  int i, swp;

  if (argc == 0) {
    i = CheckEndian(NULL);
  } else {
    f = FOPEN(argv[0], "rb");
    if (f == NULL) {
      printf("Cannot open file %s\n", argv[0]);
      return -1;
    }
    ReadFHeader(f, &fh, &swp);
    i = CheckEndian(&fh);
    FCLOSE(f);
  }

  printf("Endian: %d\n", i);

  return 0;
}

static int IntFromList(char *argv, int tp, ARRAY *variables, int **k) {
  int i;
  char *v[MAXNARGS];
  int t[MAXNARGS], n;
  
  if (tp == LIST) {
    n = DecodeArgs(argv, v, t, variables);
    if (n > 0) {
      *k = malloc(sizeof(int)*n);
      for (i = 0; i < n; i++) {
	(*k)[i] = atoi(v[i]);
	free(v[i]);
      }
    }
  } else {
    n = 1;
    *k = malloc(sizeof(int));
    (*k)[0] = atoi(argv);
  }
  
  return n;
}

static int DoubleFromList(char *argv, int tp, ARRAY *variables, double **k) {
  int i;
  char *v[MAXNARGS];
  int t[MAXNARGS], n;
  
  if (tp == LIST) {
    n = DecodeArgs(argv, v, t, variables);
    if (n > 0) {
      *k = malloc(sizeof(double)*n);
      for (i = 0; i < n; i++) {
	(*k)[i] = atof(v[i]);
	free(v[i]);
      }
    }
  } else {
    n = 1;
    *k = malloc(sizeof(double));
    (*k)[0] = atof(argv);
  }
  
  return n;
}

static int PEleDist(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  int n;

  if (argc != 2) return -1;
  n = atoi(argv[1]);
  
  EleDist(argv[0], n);
  
  return 0;
}

static int PCxtDist(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  int n;

  if (argc != 2) return -1;
  n = atoi(argv[1]);
  
  CxtDist(argv[0], n);
  
  return 0;
}

static int PSetEleDist(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int k, i, np;
  double *p = NULL;

  if (argc < 1) return -1;

  i = atoi(argv[0]);
  if (i == -1 ) {
    np = DistFromFile(argv[1], &p);
    if (np <= 0) return -1;
  } else {
    np = argc - 1;
    if (np > 0) p = (double *) malloc(sizeof(double)*np);
    for (k = 1; k < argc; k++) {
      p[k-1] = atof(argv[k]);
    }
  }
  if (SetEleDist(i, np, p) < 0) return -1;
  if (np > 0) free(p);
  return 0;
}

static int PSetCxtDist(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int k, i, np;
  double *p = NULL;

  if (argc < 1) return -1;

  i = atoi(argv[0]);
  if (i == -1 ) {
    np = DistFromFile(argv[1], &p);
    if (np <= 0) return -1;
  } else {
    np = argc - 1;
    if (np > 0) p = (double *) malloc(sizeof(double)*np);
    for (k = 1; k < argc; k++) {
      p[k-1] = atof(argv[k]);
    }
  }
  if (SetCxtDist(i, np, p) < 0) return -1;
  if (np > 0) free(p);
  return 0;
}

static int PPhoDist(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  int n;

  if (argc != 2) return -1;
  n = atoi(argv[1]);
  
  PhoDist(argv[0], n);
  
  return 0;
}

static int PSetPhoDist(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int k, i, np;
  double *p = NULL;

  if (argc < 1) return -1;

  i = atoi(argv[0]);
  if (i == -1) {
    np = DistFromFile(argv[1], &p);
    if (np <= 0) return -1;
  } else {
    np = argc - 1;
    if (np > 0) p = (double *) malloc(sizeof(double)*np);
    for (k = 1; k < argc; k++) {
      p[k-1] = atof(argv[k]);
    }
  }
  if (SetPhoDist(i, np, p) < 0) return -1;
  if (np > 0) free(p);
  return 0;
}

static int PSetNumSingleBlocks(int argc, char *argv[], int argt[], 
			       ARRAY *variables) {
  int n;

  if (argc != 1) return -1;
  n = atoi(argv[0]);
  SetNumSingleBlocks(n);
  return 0;
}

static int PSetExtrapolate(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  int n;

  if (argc != 1) return -1;
  n = atoi(argv[0]);
  SetExtrapolate(n);
  return 0;
}

static int PSetEMinAI(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  double e;

  if (argc != 1) return -1;
  e = atof(argv[0]);
  SetEMinAI(e);
  return 0;
}

static int PSetInnerAuger(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  int n;

  if (argc != 1) return -1;
  n = atoi(argv[0]);
  SetInnerAuger(n);
  return 0;
}

static int PSetEleDensity(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  double den;

  if (argc != 1) return -1;
  den = atof(argv[0]);
  SetEleDensity(den);
  return 0;
}

static int PSetCxtDensity(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  double den;

  if (argc != 1) return -1;
  den = atof(argv[0]);
  SetCxtDensity(den);
  return 0;
}
 
static int PSetPhoDensity(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  double den;

  if (argc != 1) return -1;
  den = atof(argv[0]);
  SetPhoDensity(den);
  return 0;
} 
 
static int PSetRateAccuracy(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  double epsrel, epsabs;

  if (argc < 1 || argc > 2) return -1;
  epsrel = atof(argv[0]);
  if (argc > 1) epsabs = atof(argv[1]);
  else epsabs = -1.0;

  SetRateAccuracy(epsrel, epsabs);
  
  return 0;
}

static int PSetCascade(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int c;
  double a;

  if (argc == 0 || argc > 2) return -1;

  c = atoi(argv[0]);
  a = 0.0;
  if (argc > 1) a = atof(argv[1]);
  
  SetCascade(c, a);
  
  return 0;
}

static int PSetIteration(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int maxiter;
  double a, s;

  if (argc < 1 || argc > 3) return -1;

  a = atof(argv[0]);
  maxiter = -1;
  s = -1.0;
  if (argc > 1) {
    s = atof(argv[1]);
    if (argc > 2) {
      maxiter = atoi(argv[2]);
    }
  }

  SetIteration(a, s, maxiter);
  return 0;
} 

static int PSetBlocks(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  char *ifn;
  double n;
  
  n = 0.0;
  ifn = NULL;
  if (argc > 0) {
    n = atof(argv[0]);
    if (argc > 1) {
      ifn = argv[1];
    }
  }
    
  SetBlocks(n, ifn);
  return 0;
}

static int PInitBlocks(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  InitBlocks();
  return 0;
}

static int PLevelPopulation(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  LevelPopulation();
  return 0;
}

static int PCascade(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  Cascade();
  return 0;
}

static int PSpecTable(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  char *fn;
  double smin;
  int rrc;

  smin = EPS10;
  rrc = 0;
  if (argc == 1) {
   fn = argv[0];
  } else if (argc == 2) {
    fn = argv[0];
    rrc = atoi(argv[1]);
  } else if (argc == 3) {
    fn = argv[0];
    rrc = atoi(argv[1]);
    smin = atof(argv[2]);
  } else {
    return -1;
  }
  SpecTable(fn, rrc, smin);
  return 0;
}

static int PSelectLines(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  double emin, emax, fmin;
  int nele, type;
  
  if (argc != 6 && argc != 7) return -1;
  nele = atoi(argv[2]);
  type = atoi(argv[3]);
  emin = atof(argv[4]);
  emax = atof(argv[5]);
  fmin = EPS6;
  if (argc == 7) fmin = atof(argv[6]);
  SelectLines(argv[0], argv[1], nele, type, emin, emax, fmin);
  return 0;
}

static int PPrepInterpSpec(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  if (argc != 11) return -1;
  int nd, nt, ds, ts;
  double d0, d1, t0, t1, smin, maxmem;

  nd = atoi(argv[0]);
  d0 = atof(argv[1]);
  d1 = atof(argv[2]);
  ds = atoi(argv[3]);
  nt = atoi(argv[4]);
  t0 = atof(argv[5]);
  t1 = atof(argv[6]);
  ts = atoi(argv[7]);
  smin = atof(argv[8]);
  maxmem = atof(argv[9]);

  PrepInterpSpec(nd, d0, d1, ds, nt, t0, t1, ts, smin, maxmem, argv[10]);
  return 0;
}

static int PInterpSpec(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  if (argc != 12) return -1;

  int nele = atoi(argv[1]);
  int type = atoi(argv[2]);
  int nmin = atoi(argv[3]);
  int nmax = atoi(argv[4]);
  double c = atof(argv[5]);
  double d = atof(argv[6]);
  double t = atof(argv[7]);
  double s = atof(argv[8]);
  int n = atoi(argv[9]);
  double emin = atof(argv[10]);
  double emax = atof(argv[11]);

  InterpSpecWF(argv[0], nele, type, nmin, nmax, c, d, t, s, n, emin, emax);
  return 0;
}

static int PPlotSpec(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  double emin, emax, de, smin;
  int nele, type;

  if (argc != 7 && argc != 8) return -1;
  
  nele = atoi(argv[2]);
  type = atoi(argv[3]);
  emin = atof(argv[4]);
  emax = atof(argv[5]);
  de = atof(argv[6]);
  if (argc == 8) {
    smin = atof(argv[7]);
  } else {
    smin = EPS6;
  }

  PlotSpec(argv[0], argv[1], nele, type, emin, emax, de, smin);
  
  return 0;
}

static int PAddIon(int argc, char *argv[], int argt[], 
		   ARRAY *variables) {
  int i;
  double n;

  if (argc != 3) return -1;

  i = atoi(argv[0]);
  n = atof(argv[1]);
  
  AddIon(i, n, argv[2]);

  return 0;
}

static int PSetCERates(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int inv;
  
  if (argc != 1) return -1;
  inv = atoi(argv[0]);
  SetCERates(inv);
  return 0;
}

static int PSetCXRates(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int m;
  char *s;
  
  if (argc < 1) return -1;
  m = atoi(argv[0]);
  s = NULL;
  if (argc > 1) {
    s = argv[1];
  }
  SetCXRates(m, s);
  return 0;
}

static int PSetTRRates(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int inv;
  
  if (argc != 1) return -1;
  inv = atoi(argv[0]);
  SetTRRates(inv);
  return 0;
}

static int PSetCIRates(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int inv;
  
  if (argc != 1) return -1;
  inv = atoi(argv[0]);
  SetCIRates(inv);
  return 0;
}

static int PSetRRRates(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int inv;
  
  if (argc != 1) return -1;
  inv = atoi(argv[0]);
  SetRRRates(inv);
  return 0;
}

static int PSetAIRates(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int inv;
  
  if (argc != 1) return -1;
  inv = atoi(argv[0]);
  SetAIRates(inv);
  return 0;
}

static int PSetAIRatesInner(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {  

  if (argc != 1) return -1;

  SetAIRatesInner(argv[0]);
  return 0;
}
static int PPrintTable(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int v;

  if (argc != 2 && argc != 3) return -1;
  if (argt[0] != STRING || argt[1] != STRING) return -1;
  
  v = 1;
  if (argc == 3) {
    if (argt[2] != NUMBER) return -1;
    v = atoi(argv[2]);
  }

  PrintTable(argv[0], argv[1], v);
  return 0;
}

static int PReinitCRM(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int m;

  if (argc == 0) m = 0;
  else m = atoi(argv[0]);

  ReinitCRM(m);
  return 0;
}

static int PRateTable(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int i, nc, md;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];
  
  if (argc != 1 && argc != 2 && argc != 3) return -1;
  if (argt[0] != STRING) return -1;
  
  nc = 0;
  md = 0;
  if (argc >= 2) {
    if (argt[1] != LIST) return -1;
    nc = DecodeArgs(argv[1], vg, ig, variables);
    if (argc == 3) md = atoi(argv[2]);
  }
  
  RateTable(argv[0], nc, vg, md);
  
  for (i = 0; i < nc; i++) {
    free(vg[i]);
  }
  
  return 0;
}

static int PTabNLTE(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  char *fn3;
  double xmin, xmax, dx;

  if (argc < 6 || argc > 7) return -1;
  xmin = atof(argv[3]);
  xmax = atof(argv[4]);
  dx = atof(argv[5]);
  if (argc == 7) {
    fn3 = argv[6];
  } else {
    fn3 = NULL;
  }
  
  TabNLTE(argv[0], argv[1], fn3, argv[2], xmin, xmax, dx);
  
  return 0;
}

static int PSetAbund(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  int nele;
  double a;
  
  if (argc != 2) return -1;
  
  nele = atoi(argv[0]);
  a = atof(argv[1]);
  
  SetAbund(nele, a);
  
  return 0;
}

static int PSetRateMultiplier(int argc, char *argv[], int argt[], 
			      ARRAY *variables) {
  int nele, t;
  double a;
  
  if (argc != 3) return -1;
  
  nele = atoi(argv[0]);
  t = atoi(argv[1]);
  a = atof(argv[2]);
  
  SetRateMultiplier(nele, t, a);
  
  return 0;
}

static int PDRBranch(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  
  DRBranch();
  
  return 0;
}

static int PDRStrength(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int n, i, m;

  if (argc < 2) return -1;

  n = atoi(argv[1]);
  if (argc == 4) {
    m = atoi(argv[2]);
    i = atoi(argv[3]);
  } else if (argc == 3) {
    m = atoi(argv[2]);
    i = 0;
  } else {
    i = 0;
    m = 0;
  }

  DRStrength(argv[0], n, m, i);
  
  return 0;
}

static int PRydBranch(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int n0, n1;

  if (argc == 3) n1 = -1;
  else if (argc == 4) n1 = atoi(argv[3]);
  else return -1;

  n0 = atoi(argv[2]);
  
  RydBranch(argv[0], argv[1], n0, n1);
  
  return 0;
}

static int PDumpRates(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int k, m, imax, a;

  if (argc < 3 || argc > 5) return -1;
  k = atoi(argv[1]);
  m = atoi(argv[2]);
  imax = -1;
  a = 0;
  if (argc > 3) {
    imax = atoi(argv[3]);
    if (argc > 4) {
      a = atoi(argv[4]);
    }
  }
  
  DumpRates(argv[0], k, m, imax, a);
  
  return 0;
}
 
static int PModifyRates(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  if (argc != 1) return -1;
  ModifyRates(argv[0]);
  
  return 0;
}
 
static int PSetUTA(int argc, char *argv[], int argt[], 
		   ARRAY *variables) {
  int m, mci;

  if (argc == 1) {
    if (argt[0] != NUMBER) return -1;
    m = atoi(argv[0]);
    mci = 1;
  } else if (argc == 2) {
    if (argt[0] != NUMBER || argt[1] != NUMBER) return -1;
    m = atoi(argv[0]);
    mci = atoi(argv[1]);
  }

  SetUTA(m, mci);
  
  return 0;
}

static int PDRSuppression(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int nmax;
  double z;
  
  if (argc != 3) return -1;

  z = atof(argv[1]);
  nmax = atoi(argv[2]);

  DRSuppression(argv[0], z, nmax);

  return 0;
}

static int PNormalizeMode(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {

  if (argc != 1) return -1;
  NormalizeMode(atoi(argv[0]));
  
  return 0;
}

static int PSetBornFormFactor(int argc, char *argv[], int argt[],
			      ARRAY *variables) {
  double te;
  char *fn;

  if (argc < 1 || argc > 2) return -1;
  if (argt[0] != NUMBER) return -1;
  te = atof(argv[0]);
  if (argc == 2) fn = argv[1];
  else fn = NULL;
  
  SetBornFormFactor(te, fn);
  
  return 0;
}

static int PSetBornMass(int argc, char *argv[], int argt[],
			ARRAY *variables) {
  double m;

  if (argc != 1) return -1;
  if (argt[0] != NUMBER) return -1;
  m = atof(argv[0]);

  SetBornMass(m);
  
  return 0;
}

static int PSetGamma3B(int argc, char *argv[], int argt[],
		       ARRAY *variables) {
  double m;

  if (argc != 1) return -1;
  if (argt[0] != NUMBER) return -1;
  m = atof(argv[0]);

  SetGamma3B(m);
  
  return 0;
}

static int PWallTime(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  int m = 0;
  if (argc < 1) return -1;
  if (argc > 1) {
    m = atoi(argv[1]);
  }
  PrintWallTime(argv[0], m);
  return 0;
}

static int PInitializeMPI(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  int n = -1;
  if (argc > 0) {
    n = atoi(argv[0]);
  }
  InitializeMPI(n, 1);

  return 0;
}

static int PMPIRank(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  int n, k;
  k = MPIRank(&n);
  MPrintf(-1, "%d of %d\n", k, n);
  return 0;
}

static int PMemUsed(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  MPrintf(-1, "mem used %g\n", msize());
  return 0;
}

static int PFinalizeMPI(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
#if USE_MPI == 1
  FinalizeMPI();
#endif
  return 0;
}

static int PSystem(int argc, char *argv[], int argt[], ARRAY *variables) {
  if (argc != 1) return 0;
  int r = system(argv[0]);
  return 0;
}

static int PSetProcID(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int id;
  if (argc != 1) return -1;
  id = atoi(argv[0]);
  SetProcID(id);
  return 0;
}

static int PReadKronos(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  char *dn, *prj, *tgt, *cxm;
  int z, k, md, ilog;
  dn = argv[0];
  z = atoi(argv[1]);
  k = atoi(argv[2]);
  prj = argv[3];
  tgt = argv[4];
  cxm = NULL;
  md = 1;
  ilog = 3;
  if (argc > 5) {
    cxm = argv[5];
    if (argc > 6) {
      md = atoi(argv[6]);
      if (argc > 7) {
	ilog = atoi(argv[7]);
      }
    }
  }
  int r = ReadKronos(dn, z, k, prj, tgt, cxm, md, ilog);
  return r;
}

static int PSetStarkZMP(int argc, char *argv[], int argt[],
			ARRAY *variables) {
  if (argc == 0) {
    SetStarkZMP(0, NULL);
    return 0;
  }
  if (argc > 1) return -1;
  int np;
  double *wp;
  np = DoubleFromList(argv[0], argt[0], variables, &wp);
  
  if (np < 3) {
    double wp1[3];
    wp1[0] = 1.0;
    wp1[1] = 1.0;
    wp1[2] = 1.0;
    if (np > 0)  {
      wp1[1] = wp[0];
      if (np > 1) {
	wp1[2] = wp[1];
      }
    }
    SetStarkZMP(1, wp1);
  } else {
    if (np%3 != 0) {
      printf("SetStarkZMP args must be multiples of 3\n");
      return -1;
    }
    SetStarkZMP(np/3, wp);
  }
  if (np > 0) free(wp);
  return 0;
}

static int PSetOption(int argc, char *argv[], int argt[],
		      ARRAY *variables) {
  int ip;
  double dp;
  if (argc != 2) return -1;
  if (argt[0] != STRING) return -1;
  if (argt[1] != NUMBER && argt[1] != STRING) return -1;

  if (argt[1] == NUMBER) {
    ip = atoi(argv[1]);
    dp = atof(argv[1]);
  } else {
    ip = 0;
    dp = 0;
  }

  SetOption(argv[0], argv[1], ip, dp);

  return 0;
}

static int PSetCXLDist(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int n;
  double *w;
  n = DoubleFromList(argv[0], argt[0], variables, &w);
  SetCXLDist(n, w);
  if (n > 0) {
    free(w);
  }
  return 0;
}
	    
static int PMemENTable(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {

  if (argc != 1) return -1;
  if (argt[0] != STRING) return -1;

  MemENTable(argv[0]);

  return 0;
}

static int PRateCoefficients(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {  
  int k0, k1, nexc, ncap, nt, nd, md;
  double d0, d1, t0, t1;

  if (argc != 12) return -1;
  k0 = atoi(argv[1]);
  k1 = atoi(argv[2]);
  nexc = atoi(argv[3]);
  ncap = atoi(argv[4]);
  nt = atoi(argv[5]);
  t0 = atof(argv[6]);
  t1 = atof(argv[7]);
  nd = atoi(argv[8]);
  d0 = atof(argv[9]);
  d1 = atof(argv[10]);
  md = atof(argv[11]);
  RateCoefficients(argv[0], k0, k1, nexc, ncap, nt, t0, t1, nd, d0, d1, md);
  return 0;
}

static METHOD methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"SetUTA", PSetUTA, METH_VARARGS}, 
  {"Exit", PExit, METH_VARARGS},
  {"CheckEndian", PCheckEndian, METH_VARARGS},
  {"DRSuppression", PDRSuppression, METH_VARARGS},
  {"EleDist", PEleDist, METH_VARARGS},
  {"CxtDist", PCxtDist, METH_VARARGS},
  {"PhoDist", PPhoDist, METH_VARARGS},
  {"SetEleDist", PSetEleDist, METH_VARARGS},
  {"SetCxtDist", PSetCxtDist, METH_VARARGS},
  {"SetPhoDist", PSetPhoDist, METH_VARARGS},
  {"SetNumSingleBlocks", PSetNumSingleBlocks, METH_VARARGS},
  {"SetExtrapolate", PSetExtrapolate, METH_VARARGS},
  {"SetEMinAI", PSetEMinAI, METH_VARARGS},
  {"SetInnerAuger", PSetInnerAuger, METH_VARARGS},
  {"SetEleDensity", PSetEleDensity, METH_VARARGS},
  {"SetCxtDensity", PSetCxtDensity, METH_VARARGS},
  {"SetPhoDensity", PSetPhoDensity, METH_VARARGS},
  {"SetCascade", PSetCascade, METH_VARARGS},
  {"SetIteration", PSetIteration, METH_VARARGS},
  {"SetRateAccuracy", PSetRateAccuracy, METH_VARARGS},
  {"SetBlocks", PSetBlocks, METH_VARARGS},
  {"RateTable", PRateTable, METH_VARARGS},
  {"AddIon", PAddIon, METH_VARARGS},
  {"SetCERates", PSetCERates, METH_VARARGS},
  {"SetCXRates", PSetCXRates, METH_VARARGS},
  {"SetTRRates", PSetTRRates, METH_VARARGS},
  {"SetCIRates", PSetCIRates, METH_VARARGS},
  {"SetRRRates", PSetRRRates, METH_VARARGS},
  {"SetAIRates", PSetAIRates, METH_VARARGS},
  {"SetAIRatesInner", PSetAIRates, METH_VARARGS},
  {"SetAbund", PSetAbund, METH_VARARGS},
  {"SetRateMultiplier", PSetRateMultiplier, METH_VARARGS},
  {"InitBlocks", PInitBlocks, METH_VARARGS},
  {"LevelPopulation", PLevelPopulation, METH_VARARGS},
  {"Cascade", PCascade, METH_VARARGS},
  {"SpecTable", PSpecTable, METH_VARARGS},
  {"PlotSpec", PPlotSpec, METH_VARARGS},
  {"PrepInterpSpec", PPrepInterpSpec, METH_VARARGS},
  {"InterpSpec", PInterpSpec, METH_VARARGS},
  {"TabNLTE", PTabNLTE, METH_VARARGS},
  {"SelectLines", PSelectLines, METH_VARARGS},
  {"PrintTable", PPrintTable, METH_VARARGS},
  {"ReinitCRM", PReinitCRM, METH_VARARGS},
  {"DRBranch", PDRBranch, METH_VARARGS},
  {"DRStrength", PDRStrength, METH_VARARGS},
  {"DumpRates", PDumpRates, METH_VARARGS},
  {"ModifyRates", PModifyRates, METH_VARARGS},
  {"RydBranch", PRydBranch, METH_VARARGS},
  {"NormalizeMode", PNormalizeMode, METH_VARARGS},
  {"SetBornFormFactor", PSetBornFormFactor, METH_VARARGS},
  {"SetBornMass", PSetBornMass, METH_VARARGS},
  {"SetGamma3B", PSetGamma3B, METH_VARARGS},
  {"WallTime", PWallTime, METH_VARARGS},
  {"InitializeMPI", PInitializeMPI, METH_VARARGS},
  {"MPIRank", PMPIRank, METH_VARARGS},
  {"MemUsed", PMemUsed, METH_VARARGS},
  {"FinalizeMPI", PFinalizeMPI, METH_VARARGS},
  {"System", PSystem, METH_VARARGS},
  {"SetProcID", PSetProcID, METH_VARARGS},
  {"ReadKronos", PReadKronos, METH_VARARGS},
  {"SetStarkZMP", PSetStarkZMP, METH_VARARGS},
  {"SetOption", PSetOption, METH_VARARGS},
  {"SetCXLDist", PSetCXLDist, METH_VARARGS},
  {"MemENTable", PMemENTable, METH_VARARGS},
  {"RateCoefficients", PRateCoefficients, METH_VARARGS},
  {"", NULL, METH_VARARGS}
};


int main(int argc, char *argv[]) {
  int i;
  FILE *f;

#if PMALLOC_CHECK >= 2
  pmalloc_open();
#endif

  InitCRM();

  SetModName("crm");

  ParseArgs(argc, argv, methods);

#if PMALLOC_CHECK >= 2
  pmalloc_check();
#endif

  return 0;
}
