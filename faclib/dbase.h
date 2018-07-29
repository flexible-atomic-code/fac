/*
 *   FAC - Flexible Atomic Code
 *   Copyright (C) 2001-2015 Ming Feng Gu
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 * 
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _DBASE_H_
#define _DBASE_H_ 1

#include <time.h>
#include <stdio.h>
#include <math.h>
#include "consts.h"

#define DB_EN 1
#define DB_TR 2
#define DB_CE 3
#define DB_RR 4
#define DB_AI 5
#define DB_CI 6
#define DB_SP 7
#define DB_RT 8
#define DB_DR 9
#define DB_AIM 10
#define DB_CIM 11
#define DB_ENF 12
#define DB_TRF 13
#define DB_CEF 14
#define DB_CEMF 15
#define NDB   15

#define LNCOMPLEX   32
#define LSNAME0 24
#define LSNAME      48
#define LNAME0 56
#define LNAME       128

typedef struct _FORM_FACTOR_ {  
  double te;
  int nk;
  double *k, *logk, *fk;
} FORM_FACTOR;

typedef struct _IDX_RECORD_ {
  int i0;
  int i1;
  long int position;
} IDX_RECORD;
  
typedef struct _F_HEADER_ {
  long int tsession;
  int version;
  int sversion;
  int ssversion;
  int type;
  float atom;
  char symbol[4];
  int nblocks;
  int nthreads;
} F_HEADER;
#define SIZE_F_HEADER (sizeof(long int)+5*sizeof(int)+sizeof(float)+4)

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
  int ibase;
  double energy;
  char ncomplex[LNCOMPLEX];
  char sname[LSNAME];
  char name[LNAME];
} EN_RECORD;
#define SIZE_EN_RECORD \
  (sizeof(short)*2+sizeof(int)*2+sizeof(double)+LNCOMPLEX+LSNAME+LNAME)
#define SIZE_EN_RECORD0 \
  (sizeof(short)*2+sizeof(int)*2+sizeof(double)+LNCOMPLEX+LSNAME0+LNAME0)
  
typedef struct _EN_SRECORD_ {
  int p;
  int j;
  int ibase;
  double energy;
} EN_SRECORD;

typedef struct _ENF_HEADER_ {
  long int position;
  long int length;
  int nele;
  int nlevels;
  double efield;
  double bfield;
  double fangle;
} ENF_HEADER;

typedef struct _ENF_RECORD_ {
  int ilev;
  double energy;
  int pbasis;
} ENF_RECORD;
#define SIZE_ENF_RECORD \
  (sizeof(int)*2 + sizeof(double))

typedef struct _TR_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int gauge;
  int mode;
  int multipole;
} TR_HEADER;

typedef struct _TRF_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int gauge;
  int mode;
  int multipole;
  double efield;
  double bfield;
  double fangle;
} TRF_HEADER;

typedef struct _TR_RECORD_ {
  int lower;
  int upper;
  float strength;
} TR_RECORD;
#define SIZE_TR_RECORD (sizeof(int)+sizeof(int)+sizeof(float))

typedef struct _TR_EXTRA_ {
  float energy;
  float sdev;
  float sci;
} TR_EXTRA;

typedef struct _TRF_RECORD_ {
  int lower;
  int upper;
  float *strength;
} TRF_RECORD;

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

typedef struct _CEF_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int n_tegrid;
  int n_egrid;
  float te0;
  double efield;
  double bfield;
  double fangle;
  double *tegrid;
  double *egrid;
} CEF_HEADER;

typedef struct _CE_RECORD_ {
  int lower;
  int upper;
  int nsub;
  float bethe;
  float born[2];
  float *params;
  float *strength;
} CE_RECORD;

typedef struct _CEF_RECORD_ {
  int lower;
  int upper;
  float bethe;
  float born[2];
  float *strength;
} CEF_RECORD;

typedef struct _CEMF_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int n_tegrid;
  int n_egrid;
  int n_thetagrid;
  int n_phigrid;
  float te0;
  double efield;
  double bfield;
  double fangle;
  double *tegrid;
  double *egrid;
  double *thetagrid;
  double *phigrid;
} CEMF_HEADER;

typedef struct _CEMF_RECORD_ {
  int lower;
  int upper;
  float *bethe;
  float *born;
  float *strength;
} CEMF_RECORD;

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
  float emin;
  int n_egrid;
  double *egrid;
} AI_HEADER;

typedef struct _AI_RECORD_ {
  int b;
  int f;
  float rate;
} AI_RECORD;

typedef struct _AIM_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  float emin;
  int n_egrid;
  double *egrid;
} AIM_HEADER;

typedef struct _AIM_RECORD_ {
  int b;
  int f;
  int nsub;
  float *rate;
} AIM_RECORD;

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

typedef struct _CIM_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int n_egrid;
  int egrid_type;
  int n_usr;
  int usr_egrid_type;
  double *egrid;
  double *usr_egrid;
} CIM_HEADER;

typedef struct _CIM_RECORD_ {
  int b;
  int f;
  int nsub;
  float *strength;
} CIM_RECORD;

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
  float rrate;
  float trate;
} SP_RECORD;
#define SIZE_SP_RECORD (sizeof(int)+sizeof(int)+sizeof(float)*4)

typedef struct _SP_EXTRA_ {
  float sdev;
} SP_EXTRA;

typedef struct _RT_HEADER_ { 
  long int position;
  long int length;
  int ntransitions;
  int iedist;
  int np_edist;
  double *p_edist;
  float eden;
  int ipdist;
  int np_pdist;
  double *p_pdist;
  float pden;
} RT_HEADER;

typedef struct _RT_RECORD_ {
  int dir;
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
  int flev;
  int ibase;
  int fbase;
  int vl;
  int j;
  float energy;
  float etrans;
  float br;
  float ai;
  float total_rate;
} DR_RECORD;  

/* these read functions interface with the binary data files.
 * they can be used in custom c/c++ codes to read the binary 
 * files directly. to do so, copy consts.h, dbase.h, and dbase.c
 * into a working directory, and compile and link dbase.c against the
 * custom code using these functions.
 */
void         *ReallocNew(void *p, int s);
int CompIdxRecord(const void *r1, const void *r2);
int ReadFHeader(TFILE *f, F_HEADER *fh, int *swp);
int ReadENHeader(TFILE *f, EN_HEADER *h, int swp);
int ReadENRecord(TFILE *f, EN_RECORD *r, int swp);
int ReadENFHeader(TFILE *f, ENF_HEADER *h, int swp);
int ReadENFRecord(TFILE *f, ENF_RECORD *r, int swp);
int ReadTRHeader(TFILE *f, TR_HEADER *h, int swp);
int ReadTRRecord(TFILE *f, TR_RECORD *r, TR_EXTRA *rx, int swp);
int ReadTRFHeader(TFILE *f, TRF_HEADER *h, int swp);
int ReadTRFRecord(TFILE *f, TRF_RECORD *r, int swp, TRF_HEADER *h);
int ReadCEHeader(TFILE *f, CE_HEADER *h, int swp);
int ReadCERecord(TFILE *f, CE_RECORD *r, int swp, CE_HEADER *h);
int ReadCEFHeader(TFILE *f, CEF_HEADER *h, int swp);
int ReadCEFRecord(TFILE *f, CEF_RECORD *r, int swp, CEF_HEADER *h);
int ReadCEMFHeader(TFILE *f, CEMF_HEADER *h, int swp);
int ReadCEMFRecord(TFILE *f, CEMF_RECORD *r, int swp, CEMF_HEADER *h);
int ReadRRHeader(TFILE *f, RR_HEADER *h, int swp);
int ReadRRRecord(TFILE *f, RR_RECORD *r, int swp, RR_HEADER *h);
int ReadAIHeader(TFILE *f, AI_HEADER *h, int swp);
int ReadAIRecord(TFILE *f, AI_RECORD *r, int swp);
int ReadAIMHeader(TFILE *f, AIM_HEADER *h, int swp);
int ReadAIMRecord(TFILE *f, AIM_RECORD *r, int swp);
int ReadCIHeader(TFILE *f, CI_HEADER *h, int swp);
int ReadCIMHeader(TFILE *f, CIM_HEADER *h, int swp);
int ReadCIRecord(TFILE *f, CI_RECORD *r, int swp, CI_HEADER *h);
int ReadCIMRecord(TFILE *f, CIM_RECORD *r, int swp, CIM_HEADER *h);
int ReadSPHeader(TFILE *f, SP_HEADER *h, int swp);
int ReadSPRecord(TFILE *f, SP_RECORD *r, SP_EXTRA *rx, int swp);
int ReadRTHeader(TFILE *f, RT_HEADER *h, int swp);
int ReadRTRecord(TFILE *f, RT_RECORD *r, int swp);
int ReadDRHeader(TFILE *f, DR_HEADER *h, int swp);
int ReadDRRecord(TFILE *f, DR_RECORD *r, int swp);

void CEMF2CEFHeader(CEMF_HEADER *mh, CEF_HEADER *h);
void CEMF2CEFRecord(CEMF_RECORD *mr, CEF_RECORD *r, CEMF_HEADER *mh, 
		    int ith, int iph);
/* to accommadate for the possible larger statistical weight of UTA levels.
 * the eqivalent 2j value for them are stored in r.ibase of the EN_RECORD for 
 * UTA. The two functions here are wrappers to determine the 2j and ibase for 
 * the level, depending on whether r.j < 0.
 */
int JFromENRecord(EN_RECORD *r);
int IBaseFromENRecord(EN_RECORD *r);

/* these are the write functions, which shouldn't be of much interest.
 * unless one needs to format the external data into FAC binary format.
 */
int WriteFHeader(TFILE *f, F_HEADER *fh);
int WriteENHeader(TFILE *f, EN_HEADER *h);
int WriteENFHeader(TFILE *f, ENF_HEADER *h);
int WriteTRHeader(TFILE *f, TR_HEADER *h);
int WriteTRFHeader(TFILE *f, TRF_HEADER *h);
int WriteCEHeader(TFILE *f, CE_HEADER *h);
int WriteCEFHeader(TFILE *f, CEF_HEADER *h);
int WriteCEMFHeader(TFILE *f, CEMF_HEADER *h);
int WriteRRHeader(TFILE *f, RR_HEADER *h);
int WriteAIHeader(TFILE *f, AI_HEADER *h);
int WriteAIMHeader(TFILE *f, AIM_HEADER *h);
int WriteCIHeader(TFILE *f, CI_HEADER *h);
int WriteCIMHeader(TFILE *f, CIM_HEADER *h);
int WriteSPHeader(TFILE *f, SP_HEADER *h);
int WriteRTHeader(TFILE *f, RT_HEADER *h);
int WriteDRHeader(TFILE *f, DR_HEADER *h);

int CheckEndian(F_HEADER *fh);
void SwapEndian(char *p, int size);
int SwapEndianFHeader(F_HEADER *h);
int InitDBase(void);
int ReinitDBase(int m);
TFILE *OpenFile(char *fn, F_HEADER *fhdr);
TFILE *OpenFileRO(char *fn, F_HEADER *fhdr, int *swp);
int CloseFile(TFILE *f, F_HEADER *fhdr);
int InitFile(TFILE *f, F_HEADER *fhdr, void *rhdr);
int DeinitFile(TFILE *f, F_HEADER *fhdr);
int PrintTable(char *ifn, char *ofn, int v);
int FreeMemENTable(void);
int MemENTable(char *fn);
int MemENFTable(char *fn);
EN_SRECORD *GetMemENTable(int *s);
EN_SRECORD *GetMemENFTable(int *s);
int WriteENRecord(TFILE *f, EN_RECORD *r);
int WriteENFRecord(TFILE *f, ENF_RECORD *r);
int PrintENTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int PrintENFTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianENHeader(EN_HEADER *h);
int SwapEndianENRecord(EN_RECORD *r);
int SwapEndianENFHeader(ENF_HEADER *h);
int SwapEndianENFRecord(ENF_RECORD *r);
int WriteTRRecord(TFILE *f, TR_RECORD *r, TR_EXTRA *rx);
double OscillatorStrength(int m, double e, double s, double *ga);
int PrintTRTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianTRHeader(TR_HEADER *h);
int SwapEndianTRRecord(TR_RECORD *r, TR_EXTRA *rx);
int WriteTRFRecord(TFILE *f, TRF_RECORD *r);
int PrintTRFTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianTRFHeader(TRF_HEADER *h);
int SwapEndianTRFRecord(TRF_RECORD *r);
int WriteCERecord(TFILE *f, CE_RECORD *r);
int PrintCETable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianCEHeader(CE_HEADER *h);
int SwapEndianCERecord(CE_RECORD *r);
int WriteCEFRecord(TFILE *f, CEF_RECORD *r);
int PrintCEFTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianCEFHeader(CEF_HEADER *h);
int SwapEndianCEFRecord(CEF_RECORD *r);
int WriteCEMFRecord(TFILE *f, CEMF_RECORD *r);
int PrintCEMFTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianCEMFHeader(CEMF_HEADER *h);
int SwapEndianCEMFRecord(CEMF_RECORD *r);
int WriteRRRecord(TFILE *f, RR_RECORD *r);
int PrintRRTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianRRHeader(RR_HEADER *h);
int SwapEndianRRRecord(RR_RECORD *r);
int WriteAIRecord(TFILE *f, AI_RECORD *r);
int PrintAITable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianAIHeader(AI_HEADER *h);
int SwapEndianAIRecord(AI_RECORD *r);
int WriteAIMRecord(TFILE *f, AIM_RECORD *r);
int PrintAIMTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianAIMHeader(AIM_HEADER *h);
int SwapEndianAIMRecord(AIM_RECORD *r);
int WriteCIRecord(TFILE *f, CI_RECORD *r);
int PrintCITable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianCIHeader(CI_HEADER *h);
int SwapEndianCIRecord(CI_RECORD *r);
int WriteCIMRecord(TFILE *f, CIM_RECORD *r);
int PrintCIMTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianCIMHeader(CIM_HEADER *h);
int SwapEndianCIMRecord(CIM_RECORD *r);
int WriteSPRecord(TFILE *f, SP_RECORD *r, SP_EXTRA *rx);
int PrintSPTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianSPHeader(SP_HEADER *h);
int SwapEndianSPRecord(SP_RECORD *r, SP_EXTRA *rx);
int WriteRTRecord(TFILE *f, RT_RECORD *r);
int PrintRTTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianRTHeader(RT_HEADER *h);
int SwapEndianRTRecord(RT_RECORD *r);
int WriteDRRecord(TFILE *f, DR_RECORD *r);
int PrintDRTable(TFILE *f1, FILE *f2, int v, int vs, int swp);
int SwapEndianDRHeader(DR_HEADER *h);
int SwapEndianDRRecord(DR_RECORD *r);
double IonDensity(char *fn, int k);
double IonRadiation(char *fn, int k, int m);
void SetUTA(int m, int mci);
int IsUTA(void);
void SetTRF(int m);
int AppendTable(char *fn);
int JoinTable(char *fn1, char *fn2, char *fn);
int TRBranch(char *fn, int i, int j, double *te, double *pa, double *ta);
int AIBranch(char *fn, int i, int j, double *te, double *pa, double *ta);
int LevelInfor(char *fn, int ilev, EN_RECORD *r0);
int FindLevelByName(char *fn, int nele, char *nc, char *cnr, char *cr);
int AdjustEnergy(int nlevs, int *ilevs, double *e, 
		 char *efn0, char *efn1, char *afn0, char *afn1);
int ISearch(int i, int n, int *ia);
void SetBornFormFactor(double te, char *fn);
int BornFormFactorTE(double *te);
FORM_FACTOR *BornFormFactor(void);
void SetBornMass(double m);
double BornMass(void);
int CodeBasisEB(int s, int m);
void DecodeBasisEB(int k, int *s, int *m);
double ClockStart(void);
double ClockLast(void);
double ClockNow(int m);
void PrintWallTime(char *s, int m);
int IsNewV114(TFILE *f);

#endif

