static char *rcsid="$Id: spol.c,v 1.5 2003/08/05 16:29:52 mfgu Exp $";
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

  if (argc != 2) return -1;
  if (argt[0] != STRING || argt[1] != STRING) return -1;

  i = SetMLevels(argv[0], argv[1]);

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
  int i;

  if (argc != 1) return -1;
  if (argt[0] != STRING) return -1;

  i = PolarizationTable(argv[0]);

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

static METHOD methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"Exit", PExit, METH_VARARGS},
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
  {"", NULL, METH_VARARGS}
};


int main(int argc, char *argv[]) {
  int i;
  FILE *f;

#ifdef PMALLOC_CHECK
  pmalloc_open();
#endif

  InitPolarization();

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
