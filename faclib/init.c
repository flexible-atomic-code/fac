static char *rcsid="$Id: init.c,v 1.5 2003/09/26 19:42:52 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "init.h"

#if FAC_DEBUG
  FILE *debug_log = NULL;
#endif

#ifdef PERFORM_STATISTICS
  FILE *perform_log = NULL;
#endif

int Info(void) {
  printf("========================================\n");
  printf("The Flexible Atomic Code (FAC)\n");
  printf("Version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
  printf("Bugs and suggestions, please contact:\n");
  printf("Ming Feng Gu, mfgu@stanford.edu\n");
  printf("========================================\n");
  return 0;
}

int InitFac(void) {
  int ierr;

#if FAC_DEBUG
  debug_log = fopen("debug.log", "w");
#endif

#ifdef PERFORM_STATISTICS
  perform_log = fopen("perform.log", "w");
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
  
  InitDBase();
  InitStructure();
  InitExcitation();
  InitRecombination();
  InitIonization();

  return 0;
}

int ReinitFac(int m_config, int m_recouple, int m_radial,
	      int m_dbase, int m_structure, int m_excitation,
	      int m_recombination, int m_ionization) {
  ReinitConfig(m_config);
  ReinitRecouple(m_recouple);
  ReinitRadial(m_radial);
  ReinitDBase(m_dbase);
  ReinitStructure(m_structure);
  ReinitExcitation(m_excitation);
  ReinitRecombination(m_recombination);
  ReinitIonization(m_ionization);

  return 0;
}
