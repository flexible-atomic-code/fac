static char *rcsid="$Id: spol.c,v 1.1 2003/07/14 16:27:35 mfgu Exp $";
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
  double e;

  if (argc != 2) return -1;
  if (argt[0] != STRING || argt[1] != NUMBER) return -1;

  e = atof(argv[1]);
  i = SetMCERates(argv[0], e);

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
  double d;

  if (argc != 2) return -1;
  if (argt[0] != STRING || argt[1] != NUMBER) return -1;

  d = atof(argv[1]);
  i = PopulationTable(argv[0], d);

  return i;
}

static METHOD methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"Exit", PExit, METH_VARARGS},
  {"SetMLevels", PSetMLevels, METH_VARARGS},
  {"SetMCERates", PSetMCERates, METH_VARARGS},
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
