#ifndef _DBASE_H_
#define _DBASE_H_ 1

#include <time.h>

#include "global.h"
#include "config.h"
#include "radial.h"


#define DB_EN 1
#define DB_TR 2
#define DB_CE 3
#define DB_RR 4
#define DB_AI 5
#define DB_CI 6

#define LNCOMPLEX   20
#define LSNAME      20
#define LNAME       50

typedef struct _F_HEADER_ {
  time_t tsession;
  int version;
  int sversion;
  int ssversion;
  int type;
  int atom;
  int nblocks;
} F_HEADER;

typedef struct _EN_HEADER_ {
  off_t position;
  off_t length;
  int nele;
  int nlevels;
} EN_HEADER;

typedef struct _EN_RECORD_ {
  short p;
  short j;
  int ilev;
  float energy;
  char ncomplex[LNCOMPLEX];
  char sname[LSNAME];
  char name[LNAME];
} EN_RECORD;


int InitDBase();
FILE *InitFile(char *fn, int type, int nele);
int CloseFile(FILE *f, int type);
int PrintTable(char *ifn, char *ofn);
int WriteENRecord(FILE *f, EN_RECORD *r);
int PrintENTable(FILE *f1, FILE *f2);

#endif

