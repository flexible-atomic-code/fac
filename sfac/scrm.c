static char *rcsid="$Id: scrm.c,v 1.25 2005/04/01 00:17:48 mfgu Exp $";
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
      printf(", ");
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
  FILE *f;
  F_HEADER fh;
  int i, swp;

  if (argc == 0) {
    i = CheckEndian(NULL);
  } else {
    f = fopen(argv[0], "rb");
    if (f == NULL) {
      printf("Cannot open file %s\n", argv[0]);
      return -1;
    }
    ReadFHeader(f, &fh, &swp);
    i = CheckEndian(&fh);
    fclose(f);
  }

  printf("Endian: %d\n", i);

  return 0;
}

static int PEleDist(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  int n;

  if (argc != 2) return -1;
  n = atoi(argv[1]);
  
  EleDist(argv[0], n);
  
  return 0;
}

static int PSetEleDist(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int k, i, np;
  double *p = NULL;

  if (argc < 1) return -1;

  i = atoi(argv[0]);
  np = argc - 1;
  if (np > 0) p = (double *) malloc(sizeof(double)*np);
  for (k = 1; k < argc; k++) {
    p[k-1] = atof(argv[k]);
  }
  if (SetEleDist(i, np, p) < 0) return -1;
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
  np = argc - 1;
  if (np > 0) p = (double *) malloc(sizeof(double)*np);
  for (k = 1; k < argc; k++) {
    p[k-1] = atof(argv[k]);
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
  int i, nc;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];
  
  if (argc != 1 && argc != 2) return -1;
  if (argt[0] != STRING) return -1;
  
  nc = 0;
  if (argc == 2) {
    if (argt[1] != LIST) return -1;
    nc = DecodeArgs(argv[1], vg, ig, variables);
  }
  
  RateTable(argv[0], nc, vg);
  
  for (i = 0; i < nc; i++) {
    free(vg[i]);
  }
  
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

static METHOD methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"SetUTA", PSetUTA, METH_VARARGS}, 
  {"Exit", PExit, METH_VARARGS},
  {"CheckEndian", PCheckEndian, METH_VARARGS},
  {"EleDist", PEleDist, METH_VARARGS},
  {"PhoDist", PPhoDist, METH_VARARGS},
  {"SetEleDist", PSetEleDist, METH_VARARGS},
  {"SetPhoDist", PSetPhoDist, METH_VARARGS},
  {"SetNumSingleBlocks", PSetNumSingleBlocks, METH_VARARGS},
  {"SetExtrapolate", PSetExtrapolate, METH_VARARGS},
  {"SetInnerAuger", PSetInnerAuger, METH_VARARGS},
  {"SetEleDensity", PSetEleDensity, METH_VARARGS},
  {"SetPhoDensity", PSetPhoDensity, METH_VARARGS},
  {"SetCascade", PSetCascade, METH_VARARGS},
  {"SetIteration", PSetIteration, METH_VARARGS},
  {"SetRateAccuracy", PSetRateAccuracy, METH_VARARGS},
  {"SetBlocks", PSetBlocks, METH_VARARGS},
  {"RateTable", PRateTable, METH_VARARGS},
  {"AddIon", PAddIon, METH_VARARGS},
  {"SetCERates", PSetCERates, METH_VARARGS},
  {"SetTRRates", PSetTRRates, METH_VARARGS},
  {"SetCIRates", PSetCIRates, METH_VARARGS},
  {"SetRRRates", PSetRRRates, METH_VARARGS},
  {"SetAIRates", PSetAIRates, METH_VARARGS},
  {"SetAIRatesInner", PSetAIRates, METH_VARARGS},
  {"SetAbund", PSetAbund, METH_VARARGS},
  {"InitBlocks", PInitBlocks, METH_VARARGS},
  {"LevelPopulation", PLevelPopulation, METH_VARARGS},
  {"Cascade", PCascade, METH_VARARGS},
  {"SpecTable", PSpecTable, METH_VARARGS},
  {"PlotSpec", PPlotSpec, METH_VARARGS},
  {"SelectLines", PSelectLines, METH_VARARGS},
  {"PrintTable", PPrintTable, METH_VARARGS},
  {"ReinitCRM", PReinitCRM, METH_VARARGS},
  {"DRBranch", PDRBranch, METH_VARARGS},
  {"DRStrength", PDRStrength, METH_VARARGS},
  {"DumpRates", PDumpRates, METH_VARARGS},
  {"", NULL, METH_VARARGS}
};


int main(int argc, char *argv[]) {
  int i;
  FILE *f;

#ifdef PMALLOC_CHECK
  pmalloc_open();
#endif

  InitCRM();

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

#ifdef PMALLOC_CHECK
  pmalloc_check();
#endif

  return 0;
}
