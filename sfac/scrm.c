static char *rcsid="$Id: scrm.c,v 1.3 2002/01/21 18:33:51 mfgu Exp $";
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
  SetEleDist(i, np, p);
  if (np > 0) free(p);
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
  SetPhoDist(i, np, p);
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

static int PSpecTable(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  char *fn;
  double smin;

  smin = EPS10;
  if (argc == 1) {
   fn = argv[0];
  } else if (argc == 2) {
    fn = argv[0];
    smin = atof(argv[1]);
  } else {
    return -1;
  }
  SpecTable(fn, smin);
  return 0;
}

static int PPlotSpec(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  double emin, emax, de, smin;
  int type;

  if (argc != 6 && argc != 7) return -1;
  
  type = atoi(argv[2]);
  emin = atof(argv[3]);
  emax = atof(argv[4]);
  de = atof(argv[5]);
  if (argc == 7) {
    smin = atof(argv[6]);
  } else {
    smin = EPS10;
  }

  PlotSpec(argv[0], argv[1], type, emin, emax, de, smin);
  
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

static int PFreeMemENTable(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  
  if (argc != 0) return -1;
  FreeMemENTable();
  return 0;
}

static int PMemENTable(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  
  if (argc != 1 || argt[0] != STRING) return -1;
  MemENTable(argv[0]);
  return 0;
}

static int PPrintTable(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int v;

  if (argc != 2 && argc != 3) return -1;
  if (argt[0] != STRING || argt[1] != STRING) return -1;
  
  v = 0;
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

  if (argc != 1) return -1;
  RateTable(argv[0]);
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

static METHOD methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"Exit", PExit, METH_VARARGS},
  {"SetEleDist", PSetEleDist, METH_VARARGS},
  {"SetPhoDist", PSetPhoDist, METH_VARARGS},
  {"SetNumSingleBlocks", PSetNumSingleBlocks, METH_VARARGS},
  {"SetEleDensity", PSetEleDensity, METH_VARARGS},
  {"SetPhoDensity", PSetPhoDensity, METH_VARARGS},
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
  {"SetAbund", PSetAbund, METH_VARARGS},
  {"InitBlocks", PInitBlocks, METH_VARARGS},
  {"LevelPopulation", PLevelPopulation, METH_VARARGS},
  {"SpecTable", PSpecTable, METH_VARARGS},
  {"PlotSpec", PPlotSpec, METH_VARARGS},
  {"FreeMemENTable", PFreeMemENTable, METH_VARARGS},
  {"MemENTable", PMemENTable, METH_VARARGS},
  {"PrintTable", PPrintTable, METH_VARARGS},
  {"ReinitCRM", PReinitCRM, METH_VARARGS},
  {"", NULL, METH_VARARGS}
};


int main(int argc, char *argv[]) {
  int i;
  FILE *f;

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
