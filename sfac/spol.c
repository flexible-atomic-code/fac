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

static char *rcsid="$Id: spol.c,v 1.10 2005/04/06 03:34:25 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "stoken.h"
#include "polarization.h"

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

static int PSetMaxLevels(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int m;

  if (argc != 1) return -1;
  if (argt[0] != NUMBER) return -1;

  m = atoi(argv[0]);

  SetMaxLevels(m);
  
  return 0;
}

static int PSetMIteration(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int m;
  double a;

  if (argc != 1 && argc != 2) return -1;
  if (argt[0] != NUMBER) return -1;
  if (argc == 2 && argt[1] != NUMBER) return -1;

  a = atof(argv[0]);
  if (argc == 2) {
    m = atoi(argv[1]);
  } else {
    m = 0;
  }

  SetMIteration(a, m);
  
  return 0;
}

static int PSetIDR(int argc, char *argv[], int argt[], 
		   ARRAY *variables) {
  int idr, ndr, i;
  double *pdr;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];
  
  if (argc < 1) return -1;
  if (argc > 2) return -1;

  if (argt[0] != NUMBER) return -1;
  if (argc == 2 && argt[1] != LIST && argt[1] != TUPLE) return -1;

  idr = atoi(argv[0]);
  ndr = DecodeArgs(argv[1], vg, ig, variables);
  if (ndr > 0) {
    pdr = (double *) malloc(sizeof(double)*ndr);
    for (i = 0; i < ndr; i++) {
      pdr[i] = atof(vg[i]);
      free(vg[i]);
    }
  } else {
    ndr = 0;
    pdr = NULL;
  }
  
  return SetIDR(idr, ndr, pdr);
}

static int PSetEnergy(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int i;
  double e, es;

  if (argc < 1 || argc > 2) return -1;
  if (argt[0] != NUMBER) return -1;
  if (argc == 2 && argt[1] != NUMBER) return -1;
  
  e = atof(argv[0]);
  if (argc == 2) {
    es = atof(argv[1]);
  } else {
    es = 0.0;
  }

  i = SetEnergy(e, es);
  
  return i;
}

static int PSetDensity(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int i;
  double d;

  if (argc != 1) return -1;
  if (argt[0] != NUMBER) return -1;
  
  d = atof(argv[0]);
  
  i = SetDensity(d);
  
  return i;
}

static int PSetMLevels(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int i;
  char *fn, *tfn;
  
  if (argc > 2) return -1;
  fn = NULL;
  tfn = NULL;
  if (argc > 0) {
    fn = argv[0];
    if (argc > 1) {
      tfn = argv[1];
    }
  }

  i = SetMLevels(fn, tfn);

  return i;
}

static int PSetMCERates(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int i;

  if (argc != 1) return -1;
  if (argt[0] != STRING) return -1;

  i = SetMCERates(argv[0]);

  return i;
}

static int PSetMAIRates(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int i;

  if (argc != 1) return -1;
  if (argt[0] != STRING) return -1;

  i = SetMAIRates(argv[0]);

  return i;
}

static int PPolarizationTable(int argc, char *argv[], int argt[], 
			      ARRAY *variables) {
  int i, n, k;
  char *ifn;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];
  double emin, emax, sth;

  if (argc < 1 || argc > 5) return -1;
  if (argt[0] != STRING) return -1;

  ifn = NULL;
  n = 0;
  emin = 0;
  emax = 0;
  sth = 0;
  if (argc > 1) {
    if (argt[1] == STRING) {
      ifn = argv[1];
    } else if (argt[1] == LIST) {
      n = DecodeArgs(argv[1], vg, ig, variables);
    } else {
      return -1;
    }
    if (argc > 2) {
      emin = atof(argv[2]);
      if (argc > 3) {
	emax = atof(argv[3]);
	if (argc > 4) {
	  sth = atof(argv[4]);
	}
      }
    }
  }

  i = PolarizationTable(argv[0], ifn, n, vg, emin, emax, sth);
  for (k = 0; k < n; k++) {
    free(vg[k]);
  }
  
  return i;
}

static int POrientation(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int i;
  double e;
  char *fn;

  if (argc > 2) return -1;
  e = 0.0;
  fn = NULL;
  if (argc > 0) {
    if (argt[0] != NUMBER) return -1;
    e = atof(argv[0]);
    if (argc > 1) {
      if (argt[1] != STRING) return -1;
      fn = argv[1];
    }
  }
  i = Orientation(fn, e);

  return i;
}

static int PPopulationTable(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  int i;

  if (argc != 1) return -1;
  if (argt[0] != STRING) return -1;

  i = PopulationTable(argv[0]);

  return i;
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
#ifdef USE_MPI
  int n = -1;
  if (argc > 0) {
    n = atoi(argv[0]);
  }
  InitializeMPI(n, 1);
#endif
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

static METHOD methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"Exit", PExit, METH_VARARGS},
  {"Orientation", POrientation, METH_VARARGS},
  {"SetIDR", PSetIDR, METH_VARARGS},
  {"SetMaxLevels", PSetMaxLevels, METH_VARARGS},
  {"SetMIteration", PSetMIteration, METH_VARARGS},
  {"SetEnergy", PSetEnergy, METH_VARARGS},
  {"SetDensity", PSetDensity, METH_VARARGS},
  {"SetMLevels", PSetMLevels, METH_VARARGS},
  {"SetMCERates", PSetMCERates, METH_VARARGS},
  {"SetMAIRates", PSetMAIRates, METH_VARARGS},
  {"PopulationTable", PPopulationTable, METH_VARARGS}, 
  {"PolarizationTable", PPolarizationTable, METH_VARARGS},  
  {"WallTime", PWallTime, METH_VARARGS},
  {"InitializeMPI", PInitializeMPI, METH_VARARGS},
  {"MPIRank", PMPIRank, METH_VARARGS},
  {"MemUsed", PMemUsed, METH_VARARGS},
  {"FinalizeMPI", PFinalizeMPI, METH_VARARGS},
  {"System", PSystem, METH_VARARGS},
  {"", NULL, METH_VARARGS}
};


int main(int argc, char *argv[]) {
  int i;
  FILE *f;

#if PMALLOC_CHECK == 1
  pmalloc_open();
#endif

  InitPolarization();

  SetModName("spol");

  if (argc == 1) {
    EvalFile(stdin, 1, methods);
  } else {
    for (i = 1; i < argc; i++) {
      f = fopen(argv[i], "r");
      if (!f) {
	printf("Cannot open file %s, Skipping\n", argv[i]);
	continue;
      }
      EvalFile(f, 0, methods);
    }
  }

#if PMALLOC_CHECK == 1
  pmalloc_check();
#endif

  return 0;
}
