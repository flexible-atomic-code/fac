static char *rcsid="$Id: init.c,v 1.1 2001/11/07 16:34:05 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "init.h"

int Infor() {
  printf("========================================\n");
  printf("The Flexible Atomic Code (FAC)\n");
  printf("Version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
  printf("Bugs and suggestions, please contact:\n");
  printf("Ming Feng Gu, mfgu@space.mit.edu\n");
  printf("========================================\n");
  return 0;
}

int InitFac() {
  int ierr;

#if FAC_DEBUG
  debug_log = fopen("debug.log", "w");
#endif

  ierr = InitConfig();
  if ( ierr < 0) {
    printf("initialize failed in InitConfig\n");
    return ierr;
  }

  InitCoulomb();
  InitAngular();
  InitRecouple();

  ierr = InitRadial();
  if (ierr < 0) {
    printf("initialize failed in InitRadial\n");
    return ierr;
  }

  InitStructure();
  InitExcitation();
  InitRecombination();
  InitIonization();

  return 0;
}
