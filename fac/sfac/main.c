static char *rcsid="$Id: main.c,v 1.2 2001/11/07 16:34:07 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "sfac.h"

int main(int argc, char *argv[]) {
  int i;
  FILE *f;

  if (InitFac() < 0) {
    printf("initialization failed\n");
    exit(1);
  }

  if (argc == 1) {
    EvalFile(stdin, 1);
  } else {
    for (i = 1; i < argc; i++) {
      f = fopen(argv[i], "r");
      if (!f) {
	printf("Cannot open file %s, Skipping\n", argv[i]);
	continue;
      }
      EvalFile(f, 0);
    }
  }

  return 0;
}

