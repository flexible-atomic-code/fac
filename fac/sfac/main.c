static char *rcsid="$Id: main.c,v 1.1 2001/11/06 23:55:22 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "sfac.h"

int main(int argc, char *argv[]) {
  int i;
  FILE *f;

#if FAC_DEBUG
  debug_log = fopen("debug.log", "w");
  if (!debug_log) {
    printf("can not open the debug log file\n");
    return -1;
  }
#endif
  
  if (InitConfig() < 0) {
    printf("initialize failed in InitConfig\n");
    return -1;
  }
  InitCoulomb();
  InitAngular();
  InitRecouple(); 
  if (InitRadial() < 0) {
    printf("initialize failed in InitRadial\n");
    return -1;
  }
  InitStructure();
  InitExcitation();
  InitRecombination();
  InitIonization();

  for (i = 1; i < argc; i++) {
    f = fopen(argv[i], "r");
    if (!f) {
      printf("Cannot open file %s, Skipping\n", argv[i]);
      continue;
    }
    EvalFile(f);
  }

  return 0;
}

