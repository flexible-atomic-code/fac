#ifndef _DBASE_H_
#define _DBASE_H_ 1

#include <time.h>
#include "global.h"
#include "coulomb.h"

#define DB_EN 1
#define DB_TR 2
#define DB_CE 3
#define DB_RR 4
#define DB_AI 5
#define DB_CI 6
#define DB_SP 7
#define DB_RT 8
#define DB_DR 9
#define NDB   9

#define LNCOMPLEX   32
#define LSNAME      24
#define LNAME       56

typedef struct _F_HEADER_ {
  long int tsession;
  int version;
  int sversion;
  int ssversion;
  int type;
  float atom;
  char symbol[4];
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
  double energy;
  char ncomplex[LNCOMPLEX];
  char sname[LSNAME];
  char name[LNAME];
} EN_RECORD;

typedef struct _EN_SRECORD_ {
  short p;
  short j;
  double energy;
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
  float te0;
  double *tegrid;
  double *egrid;
  double *usr_egrid;
} CE_HEADER;

typedef struct _CE_RECORD_ {
  int lower;
  int upper;
  int nsub;
  float bethe;
  float born[2];
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
  int kl;
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
  int kl;
  float *params;
  float *strength;
} CI_RECORD;

typedef struct _SP_HEADER_ { 
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int iblock;
  int fblock;
  char icomplex[LNCOMPLEX];
  char fcomplex[LNCOMPLEX];
  int type;
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
  int iblock;
  int ilev;
  char icomplex[LNCOMPLEX];
  int iedist;
  int np_edist;
  double *p_edist;
  float eden;
  int ipdist;
  int np_pdist;
  double *p_pdist;
  float pden;
  float nb;
} RT_HEADER;

typedef struct _RT_RECORD_ {
  int iblock;
  float nb;
  float tr;
  float ce;
  float rr;
  float ai;
  float ci;
  char icomplex[LNCOMPLEX];
} RT_RECORD;

typedef struct _DR_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ilev;
  int ntransitions;
  int vn;
  int j;
  float energy;
} DR_HEADER;

typedef struct _DR_RECORD_ {
  int ilev;
  int ibase;
  int flev;
  short vl;
  short j;
  float energy;
  float br;
  float ai;
  float total_rate;
} DR_RECORD;  

int CheckEndian(F_HEADER *fh);
void SwapEndian(char *p, int size);
int SwapEndianFHeader(F_HEADER *h);
int InitDBase(void);
int ReinitDBase(int m);
FILE *OpenFile(char *fn, F_HEADER *fhdr);
int CloseFile(FILE *f, F_HEADER *fhdr);
int InitFile(FILE *f, F_HEADER *fhdr, void *rhdr);
int DeinitFile(FILE *f, F_HEADER *fhdr);
int PrintTable(char *ifn, char *ofn, int v);
int WriteENRecord(FILE *f, EN_RECORD *r);
int FreeMemENTable(void);
int MemENTable(char *fn);
int LevelInfor(char *fn, int ilev, EN_RECORD *r0);
int FindLevelByName(char *fn, int nele, char *nc, char *cnr, char *cr);
int PrintENTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianENHeader(EN_HEADER *h);
int SwapEndianENRecord(EN_RECORD *r);
int WriteTRRecord(FILE *f, TR_RECORD *r);
int PrintTRTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianTRHeader(TR_HEADER *h);
int SwapEndianTRRecord(TR_RECORD *r);
int WriteCERecord(FILE *f, CE_RECORD *r);
int PrintCETable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianCEHeader(CE_HEADER *h);
int SwapEndianCERecord(CE_RECORD *r);
int WriteRRRecord(FILE *f, RR_RECORD *r);
int PrintRRTable(FILE *f1, FILE *f2, int v, int swp);
int TotalRRCross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int n0, int n1, int nmax);
int SwapEndianRRHeader(RR_HEADER *h);
int SwapEndianRRRecord(RR_RECORD *r);
int WriteAIRecord(FILE *f, AI_RECORD *r);
int PrintAITable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianAIHeader(AI_HEADER *h);
int SwapEndianAIRecord(AI_RECORD *r);
int WriteCIRecord(FILE *f, CI_RECORD *r);
int PrintCITable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianCIHeader(CI_HEADER *h);
int SwapEndianCIRecord(CI_RECORD *r);
int WriteSPRecord(FILE *f, SP_RECORD *r);
int PrintSPTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianSPHeader(SP_HEADER *h);
int SwapEndianSPRecord(SP_RECORD *r);
int WriteRTRecord(FILE *f, RT_RECORD *r);
int PrintRTTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianRTHeader(RT_HEADER *h);
int SwapEndianRTRecord(RT_RECORD *r);
int WriteDRRecord(FILE *f, DR_RECORD *r);
int PrintDRTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianDRHeader(DR_HEADER *h);
int SwapEndianDRRecord(DR_RECORD *r);

#endif

