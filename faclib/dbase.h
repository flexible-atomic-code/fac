#ifndef _DBASE_H_
#define _DBASE_H_ 1

#include <time.h>

#include "global.h"

#define DB_EN 1
#define DB_TR 2
#define DB_CE 3
#define DB_RR 4
#define DB_AI 5
#define DB_CI 6
#define DB_SP 7
#define DB_RT 8
#define NDB   DB_RT

#define LNCOMPLEX   20
#define LSNAME      20
#define LNAME       50

typedef struct _F_HEADER_ {
  time_t tsession;
  int version;
  int sversion;
  int ssversion;
  int type;
  char symbol[4];
  int atom;
  int nblocks;
} F_HEADER;

typedef struct _EN_HEADER_ {
  long int position;
  long int length;
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

typedef struct _EN_SRECORD_ {
  short p;
  short j;
  float energy;
} EN_SRECORD;

typedef struct _TR_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int gauge;
  int mode;
  int multipole;
} TR_HEADER;

typedef struct _TR_RECORD_ {
  int lower;
  int upper;
  float strength;
} TR_RECORD;

typedef struct _CE_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int qk_mode;
  int n_tegrid;
  int n_egrid;
  int egrid_type;
  int n_usr;
  int usr_egrid_type;
  int nparams;
  int pw_type;
  int msub;
  double *tegrid;
  double *egrid;
  double *usr_egrid;
} CE_HEADER;

typedef struct _CE_RECORD_ {
  int lower;
  int upper;
  int nsub;
  float *params;
  float *strength;
} CE_RECORD;

typedef struct _RR_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int qk_mode;
  int multipole;
  int n_tegrid;
  int n_egrid;
  int egrid_type;
  int n_usr;
  int usr_egrid_type;
  int nparams;
  double *tegrid;
  double *egrid;
  double *usr_egrid;
} RR_HEADER;

typedef struct _RR_RECORD_ {
  int b;
  int f;
  int nshells;
  float *params;
  float *strength;
} RR_RECORD;

typedef struct _AI_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int channel;
  int n_egrid;
  double *egrid;
} AI_HEADER;

typedef struct _AI_RECORD_ {
  int b;
  int f;
  float rate;
} AI_RECORD;

typedef struct _CI_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int qk_mode;
  int n_tegrid;
  int n_egrid;
  int egrid_type;
  int n_usr;
  int usr_egrid_type;
  int nparams;
  int pw_type;
  double *tegrid;
  double *egrid;
  double *usr_egrid;
} CI_HEADER;

typedef struct _CI_RECORD_ {
  int b;
  int f;
  int nshells;
  float *params;
  float *strength;
} CI_RECORD;

typedef struct _SP_HEADER_ { 
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int type;
  int iedist;
  int np_edist;
  double *p_edist;
  double eden;
  int ipdist;
  int np_pdist;
  double *p_pdist;
  double pden;
} SP_HEADER;

typedef struct _SP_RECORD_ {
  int lower;
  int upper;
  float energy;
  float strength;
} SP_RECORD;

typedef struct _RT_HEADER_ { 
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int iedist;
  int np_edist;
  double *p_edist;
  int ipdist;
  int np_pdist;
  double *p_pdist;
  int iblock;
  double nb;
} RT_HEADER;

typedef struct _RT_RECORD_ {
  int iblock;
  float tr;
  float ce;
  float rr;
  float ai;
  float ci;
} RT_RECORD;

int InitDBase(void);
int ReinitDBase(int m);
FILE *InitFile(char *fn, F_HEADER *fhdr, void *rhdr);
int CloseFile(FILE *f, F_HEADER *fhdr);
int PrintTable(char *ifn, char *ofn, int v);
int WriteENRecord(FILE *f, EN_RECORD *r);
int FreeMemENTable(void);
int MemENTable(char *fn);
int PrintENTable(FILE *f1, FILE *f2, int v);
int WriteTRRecord(FILE *f, TR_RECORD *r);
int PrintTRTable(FILE *f1, FILE *f2, int v);
int WriteCERecord(FILE *f, CE_RECORD *r);
int PrintCETable(FILE *f1, FILE *f2, int v);
int WriteRRRecord(FILE *f, RR_RECORD *r);
int PrintRRTable(FILE *f1, FILE *f2, int v);
int WriteAIRecord(FILE *f, AI_RECORD *r);
int PrintAITable(FILE *f1, FILE *f2, int v);
int WriteCIRecord(FILE *f, CI_RECORD *r);
int PrintCITable(FILE *f1, FILE *f2, int v);
int WriteSPRecord(FILE *f, SP_RECORD *r);
int PrintSPTable(FILE *f1, FILE *f2, int v);
int WriteRTRecord(FILE *f, RT_RECORD *r);
int PrintRTTable(FILE *f1, FILE *f2, int v);

#endif

