#include "dbase.h"

int SetDBase(int i, char *s, char mode) {
  return 0;
}

int ClearDBase(int i) {
  return 0;
}

int InitDBase() {
  int i;

#if FAC_DEBUG
  debug_log = fopen("debug.log", "w");
  if (!debug_log) return -2;
#endif

  return 0;  
}








