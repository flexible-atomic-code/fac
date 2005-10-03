#include "dbase.h"

static char *rcsid="$Id: dbase.c,v 1.83 2005/10/03 20:42:33 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int version_read[NDB];
static F_HEADER fheader[NDB];
static EN_HEADER en_header;
static TR_HEADER tr_header;
static CE_HEADER ce_header;
static RR_HEADER rr_header;
static AI_HEADER ai_header;
static AIM_HEADER aim_header;
static CI_HEADER ci_header;
static CIM_HEADER cim_header;
static SP_HEADER sp_header;
static RT_HEADER rt_header;
static DR_HEADER dr_header;

static EN_SRECORD *mem_en_table = NULL;
static int mem_en_table_size = 0;
static int iground;
static int iuta = 0;
static int utaci = 1;
static int itrf = 0;

#define _WSF0(sv, f) {					\
    n = fwrite(&(sv), sizeof(sv), 1, f);		\
    if (n != 1) return 0;				\
    m += sizeof(sv);					\
  }while(0)
#define _RSF0(sv, f) {					\
    n = fread(&(sv), sizeof(sv), 1, f);			\
    if (n != 1) return 0;				\
    m += sizeof(sv);					\
  }while(0)
#define _WSF1(sv, s, k, f) {				\
    n = fwrite(sv, s, k, f);				\
    if ((n) != (k)) return 0;				\
    m += (s)*(k);					\
  }while(0)
#define _RSF1(sv, s, k, f) {				\
    n = fread(sv, s, k, f);				\
    if ((n) != (k)) return 0;				\
    m += (s)*(k);					\
  }while(0)
#define WSF0(sv) _WSF0(sv, f)
#define WSF1(sv, s, k) _WSF1(sv, s, k, f)
#define RSF0(sv) _RSF0(sv, f)
#define RSF1(sv, s, k) _RSF1(sv, s, k, f)

void SetVersionRead(int t, int v) {
  version_read[t-1] = v;
}

void SetTRF(int m) {
  itrf = m;
}

void SetUTA(int m, int mci) {
  iuta = m;
  utaci = mci;
}

int IsUTA(void) {
  return iuta;
}

int CheckEndian(F_HEADER *fh) {
  unsigned short t = 0x01;
  char *p;

  if (fh) {
    if (fh->version > 0 || fh->sversion >= 7) {
      return (int) (fh->symbol[3]);
    }
  }
       
  p = (char *) &t;
  p += sizeof(unsigned short)-1;
  if ((unsigned short) (*p) == 1) return 1;
  else return 0;
}
 
void SwapEndian(char *p, int size) {
  int t1, t2;
  char tmp;

  t1 = 0;
  t2 = size-1;
  while (t2 > t1) {
    tmp = p[t1];
    p[t1] = p[t2];
    p[t2] = tmp;
    t1++;
    t2--;
  }
}

int SwapEndianFHeader(F_HEADER *h) {
  SwapEndian((char *) &(h->tsession), sizeof(long int));
  SwapEndian((char *) &(h->version), sizeof(int));
  SwapEndian((char *) &(h->sversion), sizeof(int));
  SwapEndian((char *) &(h->ssversion), sizeof(int));
  SwapEndian((char *) &(h->type), sizeof(int));
  SwapEndian((char *) &(h->atom), sizeof(float));
  SwapEndian((char *) &(h->nblocks), sizeof(int));
  return 0;
}

int SwapEndianENHeader(EN_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->nlevels), sizeof(int));
  return 0;
}

int SwapEndianENRecord(EN_RECORD *r) {
  SwapEndian((char *) &(r->p), sizeof(short));
  SwapEndian((char *) &(r->j), sizeof(short));
  SwapEndian((char *) &(r->ilev), sizeof(int));
  SwapEndian((char *) &(r->ibase), sizeof(int));
  SwapEndian((char *) &(r->energy), sizeof(double));
  return 0;
}

int SwapEndianTRHeader(TR_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->gauge), sizeof(int));
  SwapEndian((char *) &(h->mode), sizeof(int));
  SwapEndian((char *) &(h->multipole), sizeof(int));
  return 0;
}

int SwapEndianTRRecord(TR_RECORD *r, TR_EXTRA *rx) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->strength), sizeof(float));
  if (iuta) {
    SwapEndian((char *) &(rx->energy), sizeof(float));
    SwapEndian((char *) &(rx->sdev), sizeof(float));
    SwapEndian((char *) &(rx->sci), sizeof(float));
  }
  return 0;
}

int SwapEndianCEHeader(CE_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->qk_mode), sizeof(int));
  SwapEndian((char *) &(h->n_tegrid), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->egrid_type), sizeof(int));
  SwapEndian((char *) &(h->n_usr), sizeof(int));
  SwapEndian((char *) &(h->usr_egrid_type), sizeof(int));
  SwapEndian((char *) &(h->nparams), sizeof(int));
  SwapEndian((char *) &(h->pw_type), sizeof(int));
  SwapEndian((char *) &(h->msub), sizeof(int));
  SwapEndian((char *) &(h->te0), sizeof(float));
  return 0;
}

int SwapEndianCERecord(CE_RECORD *r) {
  int m;

  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->nsub), sizeof(int));
  SwapEndian((char *) &(r->bethe), sizeof(float));
  for (m = 0; m < 2; m++) {
    SwapEndian((char *) &(r->born[m]), sizeof(float));
  }
  return 0;
}
 
int SwapEndianRRHeader(RR_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->qk_mode), sizeof(int));
  SwapEndian((char *) &(h->multipole), sizeof(int));
  SwapEndian((char *) &(h->n_tegrid), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->egrid_type), sizeof(int));
  SwapEndian((char *) &(h->n_usr), sizeof(int));
  SwapEndian((char *) &(h->usr_egrid_type), sizeof(int));
  SwapEndian((char *) &(h->nparams), sizeof(int));
  return 0;
}
  
int SwapEndianRRRecord(RR_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->kl), sizeof(int));
  return 0;
}

int SwapEndianAIHeader(AI_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->channel), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  return 0;
}

int SwapEndianAIMHeader(AIM_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->channel), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  return 0;
}

int SwapEndianAIRecord(AI_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->rate), sizeof(float));
  return 0;
}

int SwapEndianAIMRecord(AIM_RECORD *r) {
  int i;

  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->nsub), sizeof(int));
  return 0;
}

int SwapEndianCIHeader(CI_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->qk_mode), sizeof(int));
  SwapEndian((char *) &(h->n_tegrid), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->egrid_type), sizeof(int));
  SwapEndian((char *) &(h->n_usr), sizeof(int));
  SwapEndian((char *) &(h->usr_egrid_type), sizeof(int));
  SwapEndian((char *) &(h->nparams), sizeof(int));
  SwapEndian((char *) &(h->pw_type), sizeof(int));
  return 0;
}

int SwapEndianCIRecord(CI_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->kl), sizeof(int));
  return 0;
}

int SwapEndianCIMHeader(CIM_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->egrid_type), sizeof(int));
  SwapEndian((char *) &(h->n_usr), sizeof(int));
  SwapEndian((char *) &(h->usr_egrid_type), sizeof(int));
  return 0;
}

int SwapEndianCIMRecord(CIM_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->nsub), sizeof(int));
  return 0;
}

int SwapEndianSPHeader(SP_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->iblock), sizeof(int));
  SwapEndian((char *) &(h->fblock), sizeof(int));
  SwapEndian((char *) &(h->type), sizeof(int));
  return 0;
}

int SwapEndianSPRecord(SP_RECORD *r, SP_EXTRA *rx) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->energy), sizeof(float));
  SwapEndian((char *) &(r->strength), sizeof(float));
  SwapEndian((char *) &(r->rrate), sizeof(float));
  SwapEndian((char *) &(r->trate), sizeof(float));
  if (iuta) {
    SwapEndian((char *) &(rx->sdev), sizeof(float));
  }
  return 0;
}

int SwapEndianRTHeader(RT_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->iblock), sizeof(int));
  SwapEndian((char *) &(h->ilev), sizeof(int));
  SwapEndian((char *) &(h->iedist), sizeof(int));
  SwapEndian((char *) &(h->np_edist), sizeof(int));
  SwapEndian((char *) &(h->eden), sizeof(float));
  SwapEndian((char *) &(h->ipdist), sizeof(int));
  SwapEndian((char *) &(h->np_pdist), sizeof(int));
  SwapEndian((char *) &(h->pden), sizeof(float));
  SwapEndian((char *) &(h->nb), sizeof(float));
  return 0;
}

int SwapEndianRTRecord(RT_RECORD *r) {
  SwapEndian((char *) &(r->iblock), sizeof(int));
  SwapEndian((char *) &(r->nb), sizeof(float));
  SwapEndian((char *) &(r->tr), sizeof(float));
  SwapEndian((char *) &(r->ce), sizeof(float));
  SwapEndian((char *) &(r->rr), sizeof(float));
  SwapEndian((char *) &(r->ai), sizeof(float));
  SwapEndian((char *) &(r->ci), sizeof(float));
  return 0;
}

int SwapEndianDRHeader(DR_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ilev), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->vn), sizeof(int));
  SwapEndian((char *) &(h->j), sizeof(int));
  SwapEndian((char *) &(h->energy), sizeof(float));
  return 0;
}

int SwapEndianDRRecord(DR_RECORD *r) {
  SwapEndian((char *) &(r->ilev), sizeof(int));
  SwapEndian((char *) &(r->flev), sizeof(int));
  SwapEndian((char *) &(r->ibase), sizeof(int));
  SwapEndian((char *) &(r->fbase), sizeof(int));
  SwapEndian((char *) &(r->vl), sizeof(int));
  SwapEndian((char *) &(r->j), sizeof(int));
  SwapEndian((char *) &(r->energy), sizeof(float));
  SwapEndian((char *) &(r->etrans), sizeof(float));
  SwapEndian((char *) &(r->br), sizeof(float));
  SwapEndian((char *) &(r->ai), sizeof(float));
  SwapEndian((char *) &(r->total_rate), sizeof(float));
}

int InitDBase(void) {
  int i;
  for (i = 0; i < NDB; i++) {
    fheader[i].tsession = (long int) time(0);
    fheader[i].version = VERSION;
    fheader[i].sversion = SUBVERSION;
    fheader[i].ssversion = SUBSUBVERSION;
    fheader[i].symbol[2] = '\0';
    fheader[i].symbol[3] = (char) (CheckEndian(NULL));
    fheader[i].type = 0;
    fheader[i].atom = 0;
    fheader[i].nblocks = 0;
  }

  mem_en_table = NULL;
  mem_en_table_size = 0;
  iground = 0;
  itrf = 0;

  return 0;
}

int ReinitDBase(int m) {
  int i;

  if (m < 0) return 0;
  if (mem_en_table) {
    free(mem_en_table);
    mem_en_table = NULL;
    mem_en_table_size = 0;
  }
  if (m == 0) {
    return InitDBase();
  } else {
    iground = 0;
    itrf = 0;
    if (m > NDB) return -1;
    i = m-1;
    fheader[i].tsession = (long int) time(0);
    fheader[i].version = VERSION;
    fheader[i].sversion = SUBVERSION;
    fheader[i].ssversion = SUBSUBVERSION;
    fheader[i].type = 0;
    fheader[i].atom = 0;
    fheader[i].nblocks = 0;
    return 0;
  }
}


int ReadENHeaderOld(FILE *f, EN_HEADER *h, int swp) {
  int n;

  n = fread(h, sizeof(EN_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianENHeader(h);
  
  return sizeof(EN_HEADER);
}
  
int ReadENRecordOld(FILE *f, EN_RECORD *r, int swp) {
  int n;

  n = fread(r, sizeof(EN_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianENRecord(r);
  
  return sizeof(EN_RECORD);
}

int ReadTRHeaderOld(FILE *f, TR_HEADER *h, int swp) {
  int n;

  n = fread(h, sizeof(TR_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianTRHeader(h);
  if (h->length/h->ntransitions > sizeof(TR_RECORD)) iuta = 1;
  else iuta =0;
  return sizeof(TR_HEADER);
}

int ReadTRRecordOld(FILE *f, TR_RECORD *r, TR_EXTRA *rx, int swp) {
  int n;
  
  n = fread(r, sizeof(TR_RECORD), 1, f);
  if (n != 1) return 0;
  if (iuta) {
    n = fread(rx, sizeof(TR_EXTRA), 1, f);
    if (n != 1) return 0;
  }
  if (swp) SwapEndianTRRecord(r, rx);
  
  return sizeof(TR_RECORD);
}

int ReadCEHeaderOld(FILE *f, CE_HEADER *h, int swp) {
  int i, n, m;

  n = fread(h, sizeof(CE_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianCEHeader(h);
  m = sizeof(CE_HEADER);

  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  n = fread(h->tegrid, sizeof(double), h->n_tegrid, f);
  if (n != h->n_tegrid) {
    free(h->tegrid);
    return 0;
  }
  m += sizeof(double)*h->n_tegrid;
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = fread(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) {
    free(h->tegrid);
    free(h->egrid);
    return 0;
  }
  m += sizeof(double)*h->n_egrid;
  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  n = fread(h->usr_egrid, sizeof(double), h->n_usr, f);
  if (n != h->n_usr) {
    free(h->tegrid);
    free(h->egrid);
    free(h->usr_egrid);
    return 0;
  }
  m += sizeof(double)*h->n_usr;
  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCERecordOld(FILE *f, CE_RECORD *r, int swp, CE_HEADER *h) {
  int i, n, m, m0;

  n = fread(r, sizeof(CE_RECORD), 1, f);
  if (n != 1) return 0;
  m0 = sizeof(CE_RECORD);
  if (swp) SwapEndianCERecord(r);
  
  if (h->msub) {
    m = r->nsub;
    r->params = (float *) malloc(sizeof(float)*m);
    n = fread(r->params, sizeof(float), m, f);
    if (n != m) {
      free(r->params);
      return 0;
    }
    if (swp) {
      for (i = 0; i < m; i++) {
	SwapEndian((char *) &(r->params[i]), sizeof(float));
      }
    }
    m0 += sizeof(float)*m;
  } else if (h->qk_mode == QK_FIT) {
    m = h->nparams * r->nsub;
    r->params = (float *) malloc(sizeof(float)*m);
    n = fread(r->params, sizeof(float), m, f);
    if (n != m) {
      free(r->params);
      return 0;
    }
    if (swp) {
      for (i = 0; i < m; i++) {
	SwapEndian((char *) &(r->params[i]), sizeof(float));
      }
    }
    m0 += sizeof(float)*m;
  }
  
  m = h->n_usr * r->nsub;
  r->strength = (float *) malloc(sizeof(float)*m);
  n = fread(r->strength, sizeof(float), m, f);
  if (n != m) {
    if (h->qk_mode) free(r->params);
    free(r->strength);
    return 0;
  }
  if (swp) {
    for (i = 0; i < m; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }
  m0 += sizeof(float)*m;
  
  return m0;
}

int ReadRRHeaderOld(FILE *f, RR_HEADER *h, int swp) {
  int i, n, m;

  n = fread(h, sizeof(RR_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianRRHeader(h);
  m = sizeof(RR_HEADER);

  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  n = fread(h->tegrid, sizeof(double), h->n_tegrid, f);
  if (n != h->n_tegrid) {
    free(h->tegrid);
    return 0;
  }
  m += sizeof(double)*h->n_tegrid;
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = fread(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) {
    free(h->tegrid);
    free(h->egrid);
    return 0;
  }
  m += sizeof(double)*h->n_egrid;
  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  n = fread(h->usr_egrid, sizeof(double), h->n_usr, f);
  if (n != h->n_usr) {
    free(h->tegrid);
    free(h->egrid);
    free(h->usr_egrid);
    return 0;
  }
  m += sizeof(double)*h->n_usr;
  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadRRRecordOld(FILE *f, RR_RECORD *r, int swp, RR_HEADER *h) {
  int i, n, m, m0;
  
  n = fread(r, sizeof(RR_RECORD), 1, f);
  if (n != 1) return 0;
  m0 = sizeof(RR_RECORD);
  if (swp) SwapEndianRRRecord(r);
  
  if (h->qk_mode == QK_FIT) {
    m = h->nparams;
    r->params = (float *) malloc(sizeof(float)*m);
    n = fread(r->params, sizeof(float), m, f);
    if (n != m) {
      free(r->params);
      return 0;
    }
    if (swp) {
      for (i = 0; i < m; i++) {
	SwapEndian((char *) &(r->params[i]), sizeof(float));
      }
    }
    m0 += sizeof(float)*m;
  }
  m = h->n_usr;
  r->strength = (float *) malloc(sizeof(float)*m);
  n = fread(r->strength, sizeof(float), m, f);
  if (n != m) {
    if (h->qk_mode) free(r->params);
    free(r->strength);
    return 0;
  }
  if (swp) {
    for (i = 0; i < m; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }
  m0 += sizeof(float)*m;

  return m0;
}

int ReadAIHeaderOld(FILE *f, AI_HEADER *h, int swp) {
  int i, n, m;

  n = fread(h, sizeof(AI_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianAIHeader(h);
  m = sizeof(AI_HEADER);
  
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = fread(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) {
    free(h->egrid);
    return 0;
  }
  if (swp) {
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
  }
  m += sizeof(double)*h->n_egrid;
  
  return m;
}

int ReadAIMHeaderOld(FILE *f, AIM_HEADER *h, int swp) {
  int i, n, m;

  n = fread(h, sizeof(AIM_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianAIMHeader(h);
  m = sizeof(AIM_HEADER);
  
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = fread(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) {
    free(h->egrid);
    return 0;
  }
  if (swp) {
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
  }
  m += sizeof(double)*h->n_egrid;
  
  return m;
}

int ReadAIRecordOld(FILE *f, AI_RECORD *r, int swp) {
  int n;

  n = fread(r, sizeof(AI_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianAIRecord(r);
  
  return sizeof(AI_RECORD);
}

int ReadAIMRecordOld(FILE *f, AIM_RECORD *r, int swp) {
  int n, i;

  n = fread(r, sizeof(AIM_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) {
    SwapEndianAIMRecord(r);
  }
  r->rate = (float *) malloc(sizeof(float)*r->nsub);
  n = fread(r->rate, sizeof(float), r->nsub, f);
  if (n != r->nsub) return 0;
  if (swp) {
    for (i = 0; i < r->nsub; i++) {
      SwapEndian((char *) &(r->rate[i]), sizeof(float));
    }
  }
  return sizeof(AIM_RECORD);
}

int ReadCIHeaderOld(FILE *f, CI_HEADER *h, int swp) {
  int i, n, m;

  n = fread(h, sizeof(CI_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianCIHeader(h);
  m = sizeof(CI_HEADER);
  
  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  n = fread(h->tegrid, sizeof(double), h->n_tegrid, f);
  if (n != h->n_tegrid) {
    free(h->tegrid);
    return 0;
  }
  m += sizeof(double)*h->n_tegrid;
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = fread(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) {
    free(h->tegrid);
    free(h->egrid);
    return 0;
  }
  m += sizeof(double)*h->n_egrid;
  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  n = fread(h->usr_egrid, sizeof(double), h->n_usr, f);
  if (n != h->n_usr) {
    free(h->tegrid);
    free(h->egrid);
    free(h->usr_egrid);
    return 0;
  }
  m += sizeof(double)*h->n_usr;
  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCIRecordOld(FILE *f, CI_RECORD *r, int swp, CI_HEADER *h) {
  int i, n, m, m0;

  n = fread(r, sizeof(CI_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianCIRecord(r);
  m0 = sizeof(CI_RECORD);

  m = h->nparams;
  r->params = (float *) malloc(sizeof(float)*m);
  n = fread(r->params, sizeof(float), m, f);
  if (n != m) {
    free(r->params);
    return 0;
  }
  if (swp) {
    for (i = 0; i < m; i++) {
      SwapEndian((char *) &(r->params[i]), sizeof(float));
    }
  }
  m0 += sizeof(float)*m;

  m = h->n_usr;
  r->strength = (float *) malloc(sizeof(float)*m);
  n = fread(r->strength, sizeof(float), m, f);
  if (n != m) {
    free(r->params);
    free(r->strength);
    return 0;
  }
  if (swp) {
    for (i = 0; i < m; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }
  m0 += sizeof(float)*m;

  return m0;
}

int ReadCIMHeaderOld(FILE *f, CIM_HEADER *h, int swp) {
  int i, n, m;

  n = fread(h, sizeof(CIM_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianCIMHeader(h);
  m = sizeof(CIM_HEADER);
  
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = fread(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) {
    free(h->egrid);
    return 0;
  }
  m += sizeof(double)*h->n_egrid;
  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  n = fread(h->usr_egrid, sizeof(double), h->n_usr, f);
  if (n != h->n_usr) {
    free(h->egrid);
    free(h->usr_egrid);
    return 0;
  }
  m += sizeof(double)*h->n_usr;
  if (swp) {
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCIMRecordOld(FILE *f, CIM_RECORD *r, int swp, CIM_HEADER *h) {
  int i, n, m, m0;

  n = fread(r, sizeof(CIM_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianCIMRecord(r);
  m0 = sizeof(CIM_RECORD);

  m = h->n_usr*r->nsub;
  r->strength = (float *) malloc(sizeof(float)*m);
  n = fread(r->strength, sizeof(float), m, f);
  if (n != m) {
    free(r->strength);
    return 0;
  }
  if (swp) {
    for (i = 0; i < m; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }
  m0 += sizeof(float)*m;
 
  return m0;
} 

int ReadSPHeaderOld(FILE *f, SP_HEADER *h, int swp) {
  int n;

  n = fread(h, sizeof(SP_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianSPHeader(h);
  if (h->length/h->ntransitions > sizeof(SP_RECORD)) iuta = 1;
  else iuta = 0;
  return sizeof(SP_HEADER);
}

int ReadSPRecordOld(FILE *f, SP_RECORD *r, SP_EXTRA *rx, int swp) {
  int n;

  n = fread(r, sizeof(SP_RECORD)-2*sizeof(float), 1, f);
  if (n != 1) return 0;
  if (iuta) {
    n = fread(rx, sizeof(SP_EXTRA), 1, f);
    if (n != 1) return 0;
  }
  if (swp) SwapEndianSPRecord(r, rx);
  return sizeof(SP_RECORD);
}

int ReadRTHeaderOld(FILE *f, RT_HEADER *h, int swp) {
  int i, n, m;

  n = fread(h, sizeof(RT_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianRTHeader(h);
  m = sizeof(RT_HEADER);
  
  h->p_edist = (double *) malloc(sizeof(double)*h->np_edist);
  n = fread(h->p_edist, sizeof(double), h->np_edist, f);
  if (n != h->np_edist) {
    free(h->p_edist);
    return 0;
  }
  m += sizeof(double)*h->np_edist;

  h->p_pdist = (double *) malloc(sizeof(double)*h->np_pdist);
  n = fread(h->p_pdist, sizeof(double), h->np_pdist, f);
  if (n != h->np_pdist) {
    free(h->p_edist);
    free(h->p_pdist);
    return 0;
  }
  m += sizeof(double)*h->np_pdist;

  if (swp) {
    for (i = 0; i < h->np_edist; i++) {
      SwapEndian((char *) &(h->p_edist[i]), sizeof(double));
    }
    for (i = 0; i < h->np_pdist; i++) {
      SwapEndian((char *) &(h->p_pdist[i]), sizeof(double));
    }
  }
  
  return m;
}

int ReadRTRecordOld(FILE *f, RT_RECORD *r, int swp) {
  int n;

  n = fread(r, sizeof(RT_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianRTRecord(r);
  return sizeof(RT_RECORD);
}

int ReadDRHeaderOld(FILE *f, DR_HEADER *h, int swp) {
  int n;
  
  n = fread(h, sizeof(DR_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianDRHeader(h);
  
  return sizeof(DR_HEADER);
}

int ReadDRRecordOld(FILE *f, DR_RECORD *r, int swp) {
  int n;

  n = fread(r, sizeof(DR_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianDRRecord(r);
  
  return sizeof(DR_RECORD);
}   

int WriteFHeader(FILE *f, F_HEADER *fh) {
  int n, m = 0;

  WSF0(fh->tsession);
  WSF0(fh->version);
  WSF0(fh->sversion);
  WSF0(fh->ssversion);
  WSF0(fh->type);
  WSF0(fh->atom);
  WSF1(fh->symbol, sizeof(char), 4);
  WSF0(fh->nblocks);

  return m;
}

int ReadFHeader(FILE *f, F_HEADER *fh, int *swp) {
  int n, m = 0;

  RSF0(fh->tsession);
  RSF0(fh->version);
  RSF0(fh->sversion);
  RSF0(fh->ssversion);
  RSF0(fh->type);
  RSF0(fh->atom);
  RSF1(fh->symbol, sizeof(char), 4);
  RSF0(fh->nblocks);
  
  *swp = 0;

  if (CheckEndian(fh) != (int) (fheader[0].symbol[3])) {
    *swp = 1;
    SwapEndianFHeader(fh);
  }

  SetVersionRead(fh->type, fh->version*100+fh->sversion*10+fh->ssversion);
  if (fh->type == DB_TR && itrf >= 0) {
    if (VersionLE(fh, 1, 0, 6)) itrf = 1;
    else itrf = 0;
  }

  return m;
}

int WriteENHeader(FILE *f, EN_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->nlevels);

  return m;
}

int WriteTRHeader(FILE *f, TR_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->gauge);
  WSF0(h->mode);
  WSF0(h->multipole);
  
  return m;
}

int WriteCEHeader(FILE *f, CE_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->qk_mode);
  WSF0(h->n_tegrid);
  WSF0(h->n_egrid);
  WSF0(h->egrid_type);
  WSF0(h->n_usr);
  WSF0(h->usr_egrid_type);
  WSF0(h->nparams);
  WSF0(h->pw_type);
  WSF0(h->msub);
  WSF0(h->te0);
  WSF1(h->tegrid, sizeof(double), h->n_tegrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  WSF1(h->usr_egrid, sizeof(double), h->n_usr);
  
  return m;
}

int WriteRRHeader(FILE *f, RR_HEADER *h) {
  int n, m = 0;
  
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->qk_mode);
  WSF0(h->multipole);
  WSF0(h->n_tegrid);
  WSF0(h->n_egrid);
  WSF0(h->egrid_type);
  WSF0(h->n_usr);
  WSF0(h->usr_egrid_type);
  WSF0(h->nparams);
  WSF1(h->tegrid, sizeof(double), h->n_tegrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  WSF1(h->usr_egrid, sizeof(double), h->n_usr);

  return m;
}

int WriteAIHeader(FILE *f, AI_HEADER *h) {
  int n, m = 0;
 
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->channel);
  WSF0(h->n_egrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  
  return m;
}

int WriteAIMHeader(FILE *f, AIM_HEADER *h) {
  int n, m = 0;
    
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->channel);
  WSF0(h->n_egrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  
  return m;
}

int WriteCIHeader(FILE *f, CI_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->qk_mode);
  WSF0(h->n_tegrid);
  WSF0(h->n_egrid);
  WSF0(h->egrid_type);
  WSF0(h->n_usr);
  WSF0(h->usr_egrid_type);
  WSF0(h->nparams);
  WSF0(h->pw_type);
  WSF1(h->tegrid, sizeof(double), h->n_tegrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  WSF1(h->usr_egrid, sizeof(double), h->n_usr);
  
  return m;
}

int WriteCIMHeader(FILE *f, CIM_HEADER *h) {
  int n, m = 0;
    
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->n_egrid);
  WSF0(h->egrid_type);
  WSF0(h->n_usr);
  WSF0(h->usr_egrid_type);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  WSF1(h->usr_egrid, sizeof(double), h->n_usr);

  return m;
}

int WriteSPHeader(FILE *f, SP_HEADER *h) {
  int n, m = 0;
     
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->iblock);
  WSF0(h->fblock);
  WSF1(h->icomplex, sizeof(char), LNCOMPLEX);
  WSF1(h->fcomplex, sizeof(char), LNCOMPLEX);  
  WSF0(h->type);

  return m;
}

int WriteRTHeader(FILE *f, RT_HEADER *h) {
  int n, m = 0;
         
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->iblock);
  WSF0(h->ilev);
  WSF1(h->icomplex, sizeof(char), LNCOMPLEX);
  WSF0(h->iedist);
  WSF0(h->np_edist);
  WSF0(h->eden);
  WSF0(h->ipdist);
  WSF0(h->np_pdist);
  WSF0(h->pden);
  WSF0(h->nb);
  WSF1(h->p_edist, sizeof(double), h->np_edist);
  WSF1(h->p_pdist, sizeof(double), h->np_pdist);
  
  return m;
}

int WriteDRHeader(FILE *f, DR_HEADER *h) {
  int n, m = 0;
         
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ilev);
  WSF0(h->ntransitions);
  WSF0(h->vn);
  WSF0(h->j);
  WSF0(h->energy);
  
  return m;
}

int WriteENRecord(FILE *f, EN_RECORD *r) {
  int n, m = 0;

  if (en_header.length == 0) {
    fheader[DB_EN-1].nblocks++;
    n = WriteENHeader(f, &en_header);
  }
  
  WSF0(r->p);
  WSF0(r->j);
  WSF0(r->ilev);
  WSF0(r->ibase);
  WSF0(r->energy);
  WSF1(r->ncomplex, sizeof(char), LNCOMPLEX);
  WSF1(r->sname, sizeof(char), LSNAME);
  WSF1(r->name, sizeof(char), LNAME);

  en_header.nlevels += 1;
  en_header.length += m;

  return m;
}

int WriteTRRecord(FILE *f, TR_RECORD *r, TR_EXTRA *rx) {
  int n, m = 0;

  if (tr_header.length == 0) {
    fheader[DB_TR-1].nblocks++;
    n = WriteTRHeader(f, &tr_header);
  }
  
  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->strength);

  if (iuta) {
    WSF0(rx->energy);
    WSF0(rx->sdev);
    WSF0(rx->sci);
  }

  tr_header.ntransitions += 1;
  tr_header.length += m;

  return m;
}

int WriteCERecord(FILE *f, CE_RECORD *r) {
  int n;
  int m0, m = 0;

  if (ce_header.length == 0) {
    fheader[DB_CE-1].nblocks++;
    n = WriteCEHeader(f, &ce_header);
  }
  
  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->nsub);
  WSF0(r->bethe);
  WSF1(r->born, sizeof(float), 2);

  if (ce_header.msub) {
    m0 = r->nsub;
  } else if (ce_header.qk_mode == QK_FIT) {
    m0 = ce_header.nparams * r->nsub;
  } else m0 = 0;
  if (m0) {
    WSF1(r->params, sizeof(float), m0);
  }
  m0 = ce_header.n_usr * r->nsub;
  WSF1(r->strength, sizeof(float), m0);

  ce_header.ntransitions += 1;
  ce_header.length += m;

  return m;
}

int WriteRRRecord(FILE *f, RR_RECORD *r) {
  int n;
  int m = 0, m0;

  if (rr_header.length == 0) {
    fheader[DB_RR-1].nblocks++;
    n = WriteRRHeader(f, &rr_header);
  }

  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->kl);

  if (rr_header.qk_mode == QK_FIT) {
    m0 = rr_header.nparams;
    WSF1(r->params, sizeof(float), m0);
  }
  m0 = rr_header.n_usr;
  WSF1(r->strength, sizeof(float), m0);

  rr_header.ntransitions += 1;
  rr_header.length += m;

  return m;
}

int WriteAIRecord(FILE *f, AI_RECORD *r) {
  int n, m = 0;

  if (ai_header.length == 0) {
    fheader[DB_AI-1].nblocks++;
    WriteAIHeader(f, &ai_header);
  }
  
  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->rate);

  ai_header.ntransitions += 1;
  ai_header.length += m;

  return m;
}

int WriteAIMRecord(FILE *f, AIM_RECORD *r) {
  int n, m = 0;

  if (aim_header.length == 0) {
    fheader[DB_AIM-1].nblocks++;
    WriteAIMHeader(f, &aim_header);
  }

  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->nsub);
  WSF1(r->rate, sizeof(float), r->nsub);
  
  aim_header.ntransitions += 1;
  aim_header.length += m;

  return m;
}

int WriteCIRecord(FILE *f, CI_RECORD *r) {
  int n;
  int m = 0, m0;

  if (ci_header.length == 0) {
    fheader[DB_CI-1].nblocks++;
    WriteCIHeader(f, &ci_header);
  }

  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->kl);
  m0 = ci_header.nparams;
  WSF1(r->params, sizeof(float), m0);
  m0 = ci_header.n_usr;
  WSF1(r->strength, sizeof(float), m0);

  ci_header.ntransitions += 1;
  ci_header.length += m;

  return m;
}

int WriteCIMRecord(FILE *f, CIM_RECORD *r) {
  int n;
  int m = 0, m0;

  if (cim_header.length == 0) {
    fheader[DB_CIM-1].nblocks++;
    WriteCIMHeader(f, &cim_header);
  }

  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->nsub);
  m0 = r->nsub*cim_header.n_usr;
  WSF1(r->strength, sizeof(float), m0);
  
  cim_header.ntransitions += 1;
  cim_header.length += m;

  return m;
}

int WriteSPRecord(FILE *f, SP_RECORD *r, SP_EXTRA *rx) {
  int n, m = 0;

  if (sp_header.length == 0) {
    fheader[DB_SP-1].nblocks++;
    WriteSPHeader(f, &sp_header);
  }
  
  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->energy);
  WSF0(r->strength);
  WSF0(r->rrate);
  WSF0(r->trate);
  if (iuta) {
    WSF0(rx->sdev);
  }

  sp_header.ntransitions += 1;
  sp_header.length += m;
  return m;
}

int WriteRTRecord(FILE *f, RT_RECORD *r) {
  int n, m = 0;

  if (rt_header.length == 0) {
    fheader[DB_RT-1].nblocks++;
    WriteRTHeader(f, &rt_header);
  }

  WSF0(r->iblock);
  WSF0(r->nb);
  WSF0(r->tr);
  WSF0(r->ce);
  WSF0(r->rr);
  WSF0(r->ai);
  WSF0(r->ci);
  WSF1(r->icomplex, sizeof(char), LNCOMPLEX);

  rt_header.ntransitions += 1;
  rt_header.length += m;

  return m;
}

int WriteDRRecord(FILE *f, DR_RECORD *r) {
  int n, m = 0;

  if (dr_header.length == 0) {
    fheader[DB_DR-1].nblocks++;
    WriteDRHeader(f, &dr_header);
  }

  WSF0(r->ilev);
  WSF0(r->flev);
  WSF0(r->ibase);
  WSF0(r->fbase);
  WSF0(r->vl);
  WSF0(r->j);
  WSF0(r->energy);
  WSF0(r->etrans);
  WSF0(r->br);
  WSF0(r->ai);
  WSF0(r->total_rate);

  dr_header.ntransitions += 1;
  dr_header.length += m;
  
  return m;
}

int ReadENHeader(FILE *f, EN_HEADER *h, int swp) {
  int n, m = 0;
  
  if (version_read[DB_EN-1] < 109) return ReadENHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->nlevels);

  if (swp) SwapEndianENHeader(h);
  
  return m;
}
  
int ReadENRecord(FILE *f, EN_RECORD *r, int swp) {
  int n, m = 0;

  if (version_read[DB_EN-1] < 109) return ReadENRecordOld(f, r, swp);

  RSF0(r->p);
  RSF0(r->j);
  RSF0(r->ilev);
  RSF0(r->ibase);
  RSF0(r->energy);
  RSF1(r->ncomplex, sizeof(char), LNCOMPLEX);
  RSF1(r->sname, sizeof(char), LSNAME);
  RSF1(r->name, sizeof(char), LNAME);

  if (swp) SwapEndianENRecord(r);
  
  return m;
}

int ReadTRHeader(FILE *f, TR_HEADER *h, int swp) {
  int n, m = 0;
  
  if (version_read[DB_TR-1] < 109) return ReadTRHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->gauge);
  RSF0(h->mode);
  RSF0(h->multipole);

  if (swp) SwapEndianTRHeader(h);

  if (h->length/h->ntransitions > SIZE_TR_RECORD) iuta = 1;
  else iuta =0;
  
  return m;
}

int ReadTRRecord(FILE *f, TR_RECORD *r, TR_EXTRA *rx, int swp) {
  int n, m = 0;
    
  if (version_read[DB_TR-1] < 109) return ReadTRRecordOld(f, r, rx, swp);

  RSF0(r->lower);
  RSF0(r->upper);
  RSF0(r->strength);

  if (iuta) {
    RSF0(rx->energy);
    RSF0(rx->sdev);
    RSF0(rx->sci);
  }

  if (swp) SwapEndianTRRecord(r, rx);

  if (utaci == 0) {
    rx->sci = 1.0;
  }

  return m;
}

int ReadCEHeader(FILE *f, CE_HEADER *h, int swp) {
  int i, n, m = 0;

  if (version_read[DB_CE-1] < 109) return ReadCEHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->qk_mode);
  RSF0(h->n_tegrid);
  RSF0(h->n_egrid);
  RSF0(h->egrid_type);
  RSF0(h->n_usr);
  RSF0(h->usr_egrid_type);
  RSF0(h->nparams);
  RSF0(h->pw_type);
  RSF0(h->msub);
  RSF0(h->te0);
  
  if (swp) SwapEndianCEHeader(h);
  
  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  RSF1(h->tegrid, sizeof(double), h->n_tegrid);
 
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  RSF1(h->usr_egrid, sizeof(double), h->n_usr);

  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCERecord(FILE *f, CE_RECORD *r, int swp, CE_HEADER *h) {
  int i, n, m = 0, m0;
  
  if (version_read[DB_CE-1] < 109) return ReadCERecordOld(f, r, swp, h);

  RSF0(r->lower);
  RSF0(r->upper);
  RSF0(r->nsub);
  RSF0(r->bethe);
  RSF1(r->born, sizeof(float), 2);

  if (swp) SwapEndianCERecord(r);
  
  if (h->msub) {
    m0 = r->nsub;
  }  else if (h->qk_mode == QK_FIT) {
    m0 = h->nparams * r->nsub;
  } else m0 = 0;
  if (m0) {
    r->params = (float *) malloc(sizeof(float)*m0);
    RSF1(r->params, sizeof(float), m0);
    if (swp) {
      for (i = 0; i < m0; i++) {
	SwapEndian((char *) &(r->params[i]), sizeof(float));
      }
    }
  }
  
  m0 = h->n_usr * r->nsub;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }
  
  return m;
}

int ReadRRHeader(FILE *f, RR_HEADER *h, int swp) {
  int i, n, m = 0;
  
  if (version_read[DB_RR-1] < 109) return ReadRRHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->qk_mode);
  RSF0(h->multipole);
  RSF0(h->n_tegrid);
  RSF0(h->n_egrid);
  RSF0(h->egrid_type);
  RSF0(h->n_usr);
  RSF0(h->usr_egrid_type);
  RSF0(h->nparams);

  if (swp) SwapEndianRRHeader(h);

  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  RSF1(h->tegrid, sizeof(double), h->n_tegrid);

  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  RSF1(h->usr_egrid, sizeof(double), h->n_usr);

  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadRRRecord(FILE *f, RR_RECORD *r, int swp, RR_HEADER *h) {
  int i, n, m = 0, m0;
  
  if (version_read[DB_RR-1] < 109) return ReadRRRecordOld(f, r, swp, h);

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->kl);

  if (swp) SwapEndianRRRecord(r);

  if (h->qk_mode == QK_FIT) {
    m0 = h->nparams;
    r->params = (float *) malloc(sizeof(float)*m0);
    RSF1(r->params, sizeof(float), m0);
    if (swp) {
      for (i = 0; i < m0; i++) {
	SwapEndian((char *) &(r->params[i]), sizeof(float));
      }
    }
  }
  m0 = h->n_usr;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }

  return m;
}

int ReadAIHeader(FILE *f, AI_HEADER *h, int swp) {
  int i, n, m = 0;

  if (version_read[DB_AI-1] < 109) return ReadAIHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->channel);
  RSF0(h->n_egrid);

  if (swp) SwapEndianAIHeader(h);
  
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);
  
  if (swp) {
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
  }
  
  return m;
}

int ReadAIMHeader(FILE *f, AIM_HEADER *h, int swp) {
  int i, n, m = 0;

  if (version_read[DB_AIM-1] < 109) return ReadAIMHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->channel);
  RSF0(h->n_egrid);

  if (swp) SwapEndianAIMHeader(h);
  
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);
  
  if (swp) {
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
  }
  
  return m;
}

int ReadAIRecord(FILE *f, AI_RECORD *r, int swp) {
  int n, m = 0;

  if (version_read[DB_AI-1] < 109) return ReadAIRecordOld(f, r, swp);

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->rate);

  if (swp) SwapEndianAIRecord(r);
  
  return m;
}

int ReadAIMRecord(FILE *f, AIM_RECORD *r, int swp) {
  int n, i, m = 0;

  if (version_read[DB_AIM-1] < 109) return ReadAIMRecordOld(f, r, swp);

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->nsub);
  if (swp) {
    SwapEndianAIMRecord(r);
  }
  
  r->rate = (float *) malloc(sizeof(float)*r->nsub);
  RSF1(r->rate, sizeof(float), r->nsub);
  if (swp) {
    for (i = 0; i < r->nsub; i++) {
      SwapEndian((char *) &(r->rate[i]), sizeof(float));
    }
  }
  return m;
}

int ReadCIHeader(FILE *f, CI_HEADER *h, int swp) {
  int i, n, m = 0;

  if (version_read[DB_CI-1] < 109) return ReadCIHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->qk_mode);
  RSF0(h->n_tegrid);
  RSF0(h->n_egrid);
  RSF0(h->egrid_type);
  RSF0(h->n_usr);
  RSF0(h->usr_egrid_type);
  RSF0(h->nparams);
  RSF0(h->pw_type);
  
  if (swp) SwapEndianCIHeader(h);
  
  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  RSF1(h->tegrid, sizeof(double), h->n_tegrid);
 
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  RSF1(h->usr_egrid, sizeof(double), h->n_usr);
 
  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCIRecord(FILE *f, CI_RECORD *r, int swp, CI_HEADER *h) {
  int i, n, m = 0, m0;
  
  if (version_read[DB_CI-1] < 109) return ReadCIRecordOld(f, r, swp, h);

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->kl);

  if (swp) SwapEndianCIRecord(r);

  m0 = h->nparams;
  r->params = (float *) malloc(sizeof(float)*m0);
  RSF1(r->params, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->params[i]), sizeof(float));
    }
  }

  m0 = h->n_usr;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }

  return m;
}

int ReadCIMHeader(FILE *f, CIM_HEADER *h, int swp) {
  int i, n, m = 0;

  if (version_read[DB_CIM-1] < 109) return ReadCIMHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->n_egrid);
  RSF0(h->egrid_type);
  RSF0(h->n_usr);
  RSF0(h->usr_egrid_type);
  
  if (swp) SwapEndianCIMHeader(h);
  
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  RSF1(h->usr_egrid, sizeof(double), h->n_usr);
 
  if (swp) {
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCIMRecord(FILE *f, CIM_RECORD *r, int swp, CIM_HEADER *h) {
  int i, n, m = 0, m0;

  if (version_read[DB_CIM-1] < 109) return ReadCIMRecordOld(f, r, swp, h);

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->nsub);

  if (swp) SwapEndianCIMRecord(r);

  m0 = r->nsub*h->n_usr;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }
  
  return m;
} 

int ReadSPHeader(FILE *f, SP_HEADER *h, int swp) {
  int n, m = 0;
     
  if (version_read[DB_SP-1] < 109) return ReadSPHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->iblock);
  RSF0(h->fblock);
  RSF1(h->icomplex, sizeof(char), LNCOMPLEX);
  RSF1(h->fcomplex, sizeof(char), LNCOMPLEX);  
  RSF0(h->type);

  if (swp) SwapEndianSPHeader(h);

  if (h->length/h->ntransitions > SIZE_SP_RECORD) iuta = 1;
  else iuta = 0;
  
  return m;
}

int ReadSPRecord(FILE *f, SP_RECORD *r, SP_EXTRA *rx, int swp) {
  int n, m = 0;

  if (version_read[DB_SP-1] < 109) return ReadSPRecordOld(f, r, rx, swp);

  RSF0(r->lower);
  RSF0(r->upper);
  RSF0(r->energy);
  RSF0(r->strength);
  RSF0(r->rrate);
  RSF0(r->trate);
  if (iuta) {
    RSF0(rx->sdev);
  }

  if (swp) SwapEndianSPRecord(r, rx);

  return m;
}

int ReadRTHeader(FILE *f, RT_HEADER *h, int swp) {
  int i, n, m = 0;

  if (version_read[DB_RT-1] < 109) return ReadRTHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->iblock);
  RSF0(h->ilev);
  RSF1(h->icomplex, sizeof(char), LNCOMPLEX);
  RSF0(h->iedist);
  RSF0(h->np_edist);
  RSF0(h->eden);
  RSF0(h->ipdist);
  RSF0(h->np_pdist);
  RSF0(h->pden);
  RSF0(h->nb);
  
  if (swp) SwapEndianRTHeader(h);
  
  h->p_edist = (double *) malloc(sizeof(double)*h->np_edist);
  RSF1(h->p_edist, sizeof(double), h->np_edist);

  h->p_pdist = (double *) malloc(sizeof(double)*h->np_pdist);
  RSF1(h->p_pdist, sizeof(double), h->np_pdist);

  if (swp) {
    for (i = 0; i < h->np_edist; i++) {
      SwapEndian((char *) &(h->p_edist[i]), sizeof(double));
    }
    for (i = 0; i < h->np_pdist; i++) {
      SwapEndian((char *) &(h->p_pdist[i]), sizeof(double));
    }
  }
  
  return m;
}

int ReadRTRecord(FILE *f, RT_RECORD *r, int swp) {
  int n, m = 0;

  if (version_read[DB_RT-1] < 109) return ReadRTRecordOld(f, r, swp);

  RSF0(r->iblock);
  RSF0(r->nb);
  RSF0(r->tr);
  RSF0(r->ce);
  RSF0(r->rr);
  RSF0(r->ai);
  RSF0(r->ci);
  RSF1(r->icomplex, sizeof(char), LNCOMPLEX);

  if (swp) SwapEndianRTRecord(r);

  return m;
}

int ReadDRHeader(FILE *f, DR_HEADER *h, int swp) {
  int n, m = 0;
           
  if (version_read[DB_DR-1] < 109) return ReadDRHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ilev);
  RSF0(h->ntransitions);
  RSF0(h->vn);
  RSF0(h->j);
  RSF0(h->energy);

  if (swp) SwapEndianDRHeader(h);
  
  return m;
}

int ReadDRRecord(FILE *f, DR_RECORD *r, int swp) {
  int n, m = 0;

  if (version_read[DB_DR-1] < 109) return ReadDRRecordOld(f, r, swp);

  RSF0(r->ilev);
  RSF0(r->flev);
  RSF0(r->ibase);
  RSF0(r->fbase);
  RSF0(r->vl);
  RSF0(r->j);
  RSF0(r->energy);
  RSF0(r->etrans);
  RSF0(r->br);
  RSF0(r->ai);
  RSF0(r->total_rate);

  if (swp) SwapEndianDRRecord(r);
  
  return m;
} 
  
FILE *OpenFile(char *fn, F_HEADER *fhdr) {
  int ihdr;
  FILE *f;

  ihdr = fhdr->type - 1;

  f = fopen(fn, "r+");
  if (f == NULL) {
    if (fheader[ihdr].nblocks > 0) {
      printf("A single file for one DB type must be used in one session.\n");
      exit(1);
    }
    f = fopen(fn, "w");
    if (f == NULL) {
      printf("cannot open file %s\n", fn);
      exit(1);
    }
  } else {
    if (fheader[ihdr].nblocks == 0) {
      fclose(f);
      f = fopen(fn, "w");
    }
  }

  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    exit(1);
  }

  fheader[ihdr].type = fhdr->type;
  strncpy(fheader[ihdr].symbol, fhdr->symbol, 2);
  fheader[ihdr].atom = fhdr->atom;
  WriteFHeader(f, &(fheader[ihdr]));

  return f;
}

int CloseFile(FILE *f, F_HEADER *fhdr) {
  int ihdr;
 
  ihdr = fhdr->type-1;
  fseek(f, 0, SEEK_SET);
  fheader[ihdr].type = fhdr->type;
  WriteFHeader(f, &(fheader[ihdr]));
  
  fclose(f);
  return 0;
}

int InitFile(FILE *f, F_HEADER *fhdr, void *rhdr) {
  EN_HEADER *en_hdr;
  TR_HEADER *tr_hdr;
  CE_HEADER *ce_hdr;
  RR_HEADER *rr_hdr;
  AI_HEADER *ai_hdr;
  AIM_HEADER *aim_hdr;
  CI_HEADER *ci_hdr;
  CIM_HEADER *cim_hdr;
  SP_HEADER *sp_hdr;
  RT_HEADER *rt_hdr;
  DR_HEADER *dr_hdr;
  long int p;
  int ihdr;

  if (f == NULL) return 0;
  
  ihdr = fhdr->type - 1;
  fseek(f, 0, SEEK_END);
  p = ftell(f);

  switch (fhdr->type) {
  case DB_EN:
    en_hdr = (EN_HEADER *) rhdr;
    en_header.position = p;
    en_header.length = 0;
    en_header.nele = en_hdr->nele;
    en_header.nlevels = 0;
    break;
  case DB_TR:
    tr_hdr = (TR_HEADER *) rhdr;
    memcpy(&tr_header, tr_hdr, sizeof(TR_HEADER));
    tr_header.position = p;
    tr_header.length = 0;
    tr_header.ntransitions = 0;
    break;
  case DB_CE:
    ce_hdr = (CE_HEADER *) rhdr;
    memcpy(&ce_header, ce_hdr, sizeof(CE_HEADER));
    ce_header.position = p;
    ce_header.length = 0;
    ce_header.ntransitions = 0;
    break;
  case DB_RR:
    rr_hdr = (RR_HEADER *) rhdr;
    memcpy(&rr_header, rr_hdr, sizeof(RR_HEADER));
    rr_header.position = p;
    rr_header.length = 0;
    rr_header.ntransitions = 0;
    break;
  case DB_AI:
    ai_hdr = (AI_HEADER *) rhdr;
    memcpy(&ai_header, ai_hdr, sizeof(AI_HEADER));
    ai_header.position = p;
    ai_header.length = 0;
    ai_header.ntransitions = 0;
    break;
  case DB_CI:    
    ci_hdr = (CI_HEADER *) rhdr;
    memcpy(&ci_header, ci_hdr, sizeof(CI_HEADER));
    ci_header.position = p;
    ci_header.length = 0;
    ci_header.ntransitions = 0;
    break;
  case DB_SP:
    sp_hdr = (SP_HEADER *) rhdr;
    memcpy(&sp_header, sp_hdr, sizeof(SP_HEADER));
    sp_header.position = p;
    sp_header.length = 0;
    sp_header.ntransitions = 0;
    break;
  case DB_RT:
    rt_hdr = (RT_HEADER *) rhdr;
    memcpy(&rt_header, rt_hdr, sizeof(RT_HEADER));
    rt_header.position = p;
    rt_header.length = 0;
    rt_header.ntransitions = 0;
    break;
  case DB_DR:
    dr_hdr = (DR_HEADER *) rhdr;
    memcpy(&dr_header, dr_hdr, sizeof(DR_HEADER));
    dr_header.position = p;
    dr_header.length = 0;
    dr_header.ntransitions = 0;
    break;
  case DB_AIM:
    aim_hdr = (AIM_HEADER *) rhdr;
    memcpy(&aim_header, aim_hdr, sizeof(AIM_HEADER));
    aim_header.position = p;
    aim_header.length = 0;
    aim_header.ntransitions = 0;
    break;
  case DB_CIM:
    cim_hdr = (CIM_HEADER *) rhdr;
    memcpy(&cim_header, cim_hdr, sizeof(CIM_HEADER));
    cim_header.position = p;
    cim_header.length = 0;
    cim_header.ntransitions = 0;
    break;
  default:
    break;
  }

  return 0;
}

int DeinitFile(FILE *f, F_HEADER *fhdr) {
  int n;

  if (f == NULL || fhdr->type <= 0) return 0;

  switch (fhdr->type) {
  case DB_EN:
    fseek(f, en_header.position, SEEK_SET);
    if (en_header.length > 0) {
      n = WriteENHeader(f, &en_header);
    }
    break;
  case DB_TR:
    fseek(f, tr_header.position, SEEK_SET);
    if (tr_header.length > 0) {
      n = WriteTRHeader(f, &tr_header);
    }
    break;
  case DB_CE:
    fseek(f, ce_header.position, SEEK_SET);
    if (ce_header.length > 0) {
      n = WriteCEHeader(f, &ce_header);
    }
    break;
  case DB_RR:
    fseek(f, rr_header.position, SEEK_SET);
    if (rr_header.length > 0) {
      n = WriteRRHeader(f, &rr_header);
    }
    break;
  case DB_AI:
    fseek(f, ai_header.position, SEEK_SET);
    if (ai_header.length > 0) {
      n = WriteAIHeader(f, &ai_header);
    }
    break;
  case DB_CI:
    fseek(f, ci_header.position, SEEK_SET);
    if (ci_header.length > 0) {
      n = WriteCIHeader(f, &ci_header);
    }
    break;
  case DB_SP:
    fseek(f, sp_header.position, SEEK_SET);
    if (sp_header.length > 0) {
      n = WriteSPHeader(f, &sp_header);
    }
    break;
  case DB_RT:
    fseek(f, rt_header.position, SEEK_SET);
    if (rt_header.length > 0) {
      n = WriteRTHeader(f, &rt_header);
    }
    break;
  case DB_DR:
    fseek(f, dr_header.position, SEEK_SET);
    if (dr_header.length > 0) {
      n = WriteDRHeader(f, &dr_header);
    }
    break;
  case DB_AIM:
    fseek(f, aim_header.position, SEEK_SET);
    if (aim_header.length > 0) {
      n = WriteAIMHeader(f, &aim_header);
    }
    break;
  case DB_CIM:
    fseek(f, cim_header.position, SEEK_SET);
    if (cim_header.length > 0) {
      n = WriteCIMHeader(f, &cim_header);
    }
    break;
  default:
    break;
  }
  return 0;
}

int PrintTable(char *ifn, char *ofn, int v) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;

  f1 = fopen(ifn, "r");
  if (f1 == NULL) return -1;

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) return -1;

  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    goto DONE;
  }

  if (v && (fh.type < DB_SP || fh.type > DB_DR)) {
    if (mem_en_table == NULL) {
      printf("Energy table has not been built in memory.\n");
      goto DONE;
    }
  }

  fprintf(f2, "FAC %d.%d.%d\n", fh.version, fh.sversion, fh.ssversion);
  fprintf(f2, "Endian\t= %d\n", (int) CheckEndian(&fh));
  fprintf(f2, "TSess\t= %lu\n", fh.tsession);
  fprintf(f2, "Type\t= %d\n", fh.type);
  fprintf(f2, "Verbose\t= %d\n", v);
  fprintf(f2, "%s Z\t= %5.1f\n", fh.symbol, fh.atom);
  fprintf(f2, "NBlocks\t= %d\n", fh.nblocks);
  switch (fh.type) {
  case DB_EN:
    if (v) {
      fprintf(f2, "E0\t= %-d, %15.8E\n", 
	      iground, (mem_en_table[iground].energy * HARTREE_EV));
    }
    n = PrintENTable(f1, f2, v, swp);
    break;
  case DB_TR:
    n = PrintTRTable(f1, f2, v, swp);
    break;
  case DB_CE:
    n = PrintCETable(f1, f2, v, swp);
    break;
  case DB_RR:
    n = PrintRRTable(f1, f2, v, swp);
    break;
  case DB_AI:
    n = PrintAITable(f1, f2, v, swp);
    break;
  case DB_CI:
    n = PrintCITable(f1, f2, v, swp);
    break;
  case DB_SP:
    n = PrintSPTable(f1, f2, v, swp);
    break;
  case DB_RT:
    n = PrintRTTable(f1, f2, v, swp);
    break;
  case DB_DR:
    n = PrintDRTable(f1, f2, v, swp);
    break;
  case DB_AIM:
    n = PrintAIMTable(f1, f2, v, swp);
    break;
  case DB_CIM:
    n = PrintCIMTable(f1, f2, v, swp);
  default:
    break;
  }

 DONE:
  fclose(f1);
  if (f2 != stdout) fclose(f2);
  else fflush(f2);

  return n;
}

int FreeMemENTable(void) {
  if (mem_en_table) free(mem_en_table);
  mem_en_table = NULL;
  mem_en_table_size = 0;
  return 0;
}

static int StrTrimCmp(char *s1, char *s2) {
  int i, j;

  i = 0;
  while (s1[i] == ' ' || s1[i] == '\t') i++;
  j = 0;
  while (s2[j] == ' ' || s2[j] == '\t') j++;
  while (s1[i] && s2[j]) {
    if (s1[i] != s2[j]) {
      return 1;
    }
    i++;
    j++;
  }
  if (s1[i] == '\0') {
    while (s2[j]) {
      if (s2[j] != ' ' && s2[j] != '\t') {
	return 1;
      }
      j++;
    }
  }
  if (s2[j] == '\0') {
    while (s1[i]) {
      if (s1[i] != ' ' && s1[i] != '\t') {
	return 1;
      }
      i++;
    }
  }
  return 0;
}
   
int FindLevelByName(char *fn, int nele, char *nc, char *cnr, char *cr) {
  F_HEADER fh;  
  EN_HEADER h;
  EN_RECORD r;
  FILE *f;
  int n, k;
  int swp;
  
  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) {
    fclose(f);
    return 0;
  }

  if (fh.type != DB_EN) {
    printf("File type is not DB_EN\n");
    fclose(f);
    return -1;
  }

  while (1) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    if (h.nele != nele) {
      fseek(f, h.length, SEEK_CUR);
      continue;
    }
    for (k = 0; k < h.nlevels; k++) {
      n = ReadENRecord(f, &r, swp);
      if (StrTrimCmp(r.ncomplex, nc) == 0 &&
	  StrTrimCmp(r.sname, cnr) == 0 &&
	  StrTrimCmp(r.name, cr) == 0) {
	fclose(f);
	return r.ilev;
      }
    }
  }
  
  fclose(f);
  return -1;
}
      
int LevelInfor(char *fn, int ilev, EN_RECORD *r0) {
  F_HEADER fh;  
  EN_HEADER h;
  EN_RECORD r;
  FILE *f;
  int n, i, k, nlevels;
  int swp, sr;
  
  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) {
    fclose(f);
    return 0;
  }
  if (fh.type != DB_EN) {
    printf("File type is not DB_EN\n");
    fclose(f);
    return -1;
  }
  if (version_read[DB_EN-1] < 109) sr = sizeof(EN_RECORD);
  else sr = SIZE_EN_RECORD;

  if (ilev >= 0) {
    k = ilev;
    nlevels = 0;
    for (i = 0; i < fh.nblocks; i++) {
      n = ReadENHeader(f, &h, swp);
      if (n == 0) break;
      nlevels += h.nlevels;
      if (k < h.nlevels) {
	if (k > 0) fseek(f, sr*k, SEEK_CUR);
	n = ReadENRecord(f, &r, swp);
	if (n == 0) break;
	if (r.ilev != ilev) {
	  fclose(f);
	  return -1;
	}
	memcpy(r0, &r, sizeof(EN_RECORD));
	break;
      } else {
	k -= h.nlevels;
	fseek(f, h.length, SEEK_CUR);
      }
    }
    fclose(f);
    if (i == fh.nblocks) return -1;
    return 0;
  } else {
    k = -ilev;
    if (k >= 1000) k = 0;
    nlevels = 0;
    for (i = 0; i < fh.nblocks; i++) {
      n = ReadENHeader(f, &h, swp);
      if (n == 0) break;
      if (h.nele == k) break;
      nlevels += h.nlevels;
      fseek(f, h.length, SEEK_CUR);
    }
    fclose(f);
    if (i == fh.nblocks) return -1;
    return nlevels;
  }  
}

EN_SRECORD *GetMemENTable(int *s) {
  *s = mem_en_table_size;
  return mem_en_table;
}

int JFromENRecord(EN_RECORD *r) {
  if (r->j < 0) return r->ibase;
  else return r->j;
}

int IBaseFromENRecord(EN_RECORD *r) {
  if (r->j < 0) return -1;
  else return r->ibase;
}


int MemENTable(char *fn) {
  F_HEADER fh;  
  EN_HEADER h;
  EN_RECORD r;
  FILE *f;
  char *s;
  int n, i, nlevels;
  float e0;
  int swp, sr;

  f = fopen(fn, "r");
  if (f == NULL) return -1;

  n = ReadFHeader(f, &fh, &swp);  
  if (n == 0) return 0;
  if (fh.type != DB_EN) return -1;
  if (version_read[DB_EN-1] < 109) sr = sizeof(EN_RECORD);
  else sr = SIZE_EN_RECORD;

  if (mem_en_table) free(mem_en_table);

  nlevels = 0;
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    if (h.length > sr) {
      fseek(f, h.length-sr, SEEK_CUR);
    }
    n = ReadENRecord(f, &r, swp);
    if (r.ilev >= nlevels) nlevels = r.ilev+1;
  }

  mem_en_table = (EN_SRECORD *) malloc(sizeof(EN_SRECORD)*nlevels);
  mem_en_table_size = nlevels;

  e0 = 0.0;
  if (version_read[DB_EN-1] < 109) {
    fseek(f, sizeof(F_HEADER), SEEK_SET);
  } else {
    fseek(f, SIZE_F_HEADER, SEEK_SET);
  }
  while (1) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.nlevels; i++) {
      n = ReadENRecord(f, &r, swp);
      if (n == 0) break;
      if (r.energy < e0) {
	e0 = r.energy;
	iground = r.ilev;
      }
      mem_en_table[r.ilev].energy = r.energy;
      mem_en_table[r.ilev].p = r.p;
      mem_en_table[r.ilev].j = JFromENRecord(&r);
    }
  }

  if (nlevels > 0) {
    s = r.name;
    iuta = 1;
    while (*s) {
      if (*s == '(') {
	iuta = 0;
	break;
      }
      s++;
    }
  }
  return 0;
}    

int PrintENTable(FILE *f1, FILE *f2, int v, int swp) {
  EN_HEADER h;
  EN_RECORD r;
  int n, i;
  int nb;
  double e;
  int p, vnl, j, ibase;

  nb = 0;
  while (1) {
    n = ReadENHeader(f1, &h, swp);
    if (n == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NLEV\t= %d\n", h.nlevels);
    fprintf(f2, "  ILEV  IBASE    ENERGY       P   VNL   2J\n");
    for (i = 0; i < h.nlevels; i++) {
      n = ReadENRecord(f1, &r, swp);
      if (n == 0) break;
      e = r.energy;
      if (v) {
	e -= mem_en_table[iground].energy;
	e *= HARTREE_EV;
      }
      if (r.p < 0) {
	p = 1;
	vnl = -r.p;
      } else {
	p = 0;
	vnl = r.p;
      }
      j = JFromENRecord(&r);
      ibase = IBaseFromENRecord(&r);
      fprintf(f2, "%6d %6d %15.8E %1d %5d %4d %-20s %-20s %-s\n",
	      r.ilev, ibase, e, p, vnl, j, r.ncomplex, r.sname, r.name);
    }
    nb += 1;
  }
  
  return nb;
}

double OscillatorStrength(int m, double e, double s, double *ga) {
  int m2;
  double aw, x;

  aw = FINE_STRUCTURE_CONST * e;
  if (itrf == 0) {
    m2 = 2*abs(m);
    x = s*s/(m2+1.0);
    x *= e;
    m2 -= 2;
    if (m2) {
      x *= pow(aw, m2);
    }
  } else {
    x = s;
  }
  if (ga) {
    *ga = x*2.0*pow(aw,2)*FINE_STRUCTURE_CONST;
  }  

  return x;
}  

int PrintTRTable(FILE *f1, FILE *f2, int v, int swp) {
  TR_HEADER h;
  TR_RECORD r;
  TR_EXTRA rx;
  int n, i;
  int nb;
  double e, a, gf;
  char s0[10], s1[10];
  char slow[20], sup[20];
  int j0, j1, n0, n1, nq0, nq1;

  nb = 0;
  
  while (1) {
    n = ReadTRHeader(f1, &h, swp);
    if (n == 0) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "MULTIP\t= %d\n", (int)h.multipole);
    fprintf(f2, "GAUGE\t= %d\n", (int)h.gauge);
    fprintf(f2, "MODE\t= %d\n", (int)h.mode);

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadTRRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      if (iuta) {
	if (v) {
	  e = rx.energy;
	  gf = OscillatorStrength(h.multipole, e, r.strength, &a);
	  a /= (mem_en_table[r.upper].j + 1.0);
	  a *= RATE_AU;
	  fprintf(f2, "%5d %4d %5d %4d %13.6E %11.4E %13.6E %13.6E %13.6E %10.3E\n",
		  r.upper, mem_en_table[r.upper].j, 
		  r.lower, mem_en_table[r.lower].j,
		  (e*HARTREE_EV), 
		  (rx.sdev*HARTREE_EV), gf, a, r.strength, rx.sci);
	} else {
	  e = rx.energy;
	  fprintf(f2, "%5d %6s %5d %6s %13.6E %11.4E %13.6E %10.3E\n",
		  r.upper, sup, r.lower, slow, e, rx.sdev, r.strength, rx.sci);
	}
      } else {
	if (v) {
	  e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	  gf = OscillatorStrength(h.multipole, e, (double)r.strength, &a);
	  a /= (mem_en_table[r.upper].j + 1.0);
	  a *= RATE_AU;
	  fprintf(f2, "%6d %2d %6d %2d %13.6E %13.6E %13.6E %13.6E\n",
		  r.upper, mem_en_table[r.upper].j,
		  r.lower, mem_en_table[r.lower].j,
		  (e*HARTREE_EV), gf, a, r.strength);
	} else {
	  fprintf(f2, "%6d %6d %13.6E\n", 
		  r.upper, r.lower, r.strength);
	}
      }
    }
    nb += 1;
  }

  return nb;
}

int TRBranch(char *fn, int upper, int lower, 
	     double *te, double *pa, double *ta) {
  F_HEADER fh;
  TR_HEADER h;
  TR_RECORD r;
  TR_EXTRA rx;
  FILE *f;
  int n, i, k;
  double a, b, c, e;
  int swp;
 
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) {
    fclose(f);
    return 0;
  }
  if (fh.type != DB_TR) {
    printf("File type is not DB_TR\n");
    fclose(f);
    return -1;
  }
  
  a = 0.0;
  c = 0.0;
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadTRHeader(f, &h, swp);
    if (n == 0) break;
    for (k = 0; k < h.ntransitions; k++) {
      n = ReadTRRecord(f, &r, &rx, swp);
      if (n == 0) break;
      if (r.upper == upper) {
	e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	OscillatorStrength(h.multipole, e, r.strength, &b);
	b /= (mem_en_table[r.upper].j + 1.0);
	b *= RATE_AU;
	a += b;
	if (r.lower == lower) {
	  c += b;
	}
      }
    }
  }
  
  *pa = c;
  *ta = a;
  if (lower >= 0) {
    *te = mem_en_table[upper].energy - mem_en_table[lower].energy;
    *te *= HARTREE_EV;
  } else {
    *te = 0.0;
  }

  fclose(f);

  return 0;
}
  
int PrintCETable(FILE *f1, FILE *f2, int v, int swp) {
  CE_HEADER h;
  CE_RECORD r;
  int n, i, t;
  int nb;
  int m, k, p1, p2;
  float a, e;

  nb = 0;
 
  while (1) {
    n = ReadCEHeader(f1, &h, swp);
    if (n == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "QKMODE\t= %d\n", h.qk_mode);
    fprintf(f2, "NPARAMS\t= %d\n", h.nparams);
    fprintf(f2, "MSUB\t= %d\n", h.msub);
    fprintf(f2, "PWTYPE\t= %d\n", h.pw_type);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);

    for (i = 0; i < h.n_tegrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "TE0\t= %15.8E\n", h.te0 * HARTREE_EV);
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (v) {
	e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	fprintf(f2, "%6d %2d %6d %2d %11.4E %d\n",
		r.lower, mem_en_table[r.lower].j,
		r.upper, mem_en_table[r.upper].j,
		e*HARTREE_EV, r.nsub);
	fprintf(f2, "%11.4E %11.4E %11.4E\n", 
		r.bethe, r.born[0], r.born[1]*HARTREE_EV);
      } else {
	fprintf(f2, "%6d %6d %d\n", 
		r.lower, r.upper, r.nsub);
	fprintf(f2, "%11.4E %11.4E %11.4E\n", 
		r.bethe, r.born[0], r.born[1]);
      }
      
      p1 = 0;
      p2 = 0;
      for (k = 0; k < r.nsub; k++) {
	if (h.msub) {
	  fprintf(f2, "%11.4E\n", r.params[k]);
	} else if (h.qk_mode == QK_FIT) {
	  for (t = 0; t < h.nparams; t++) {
	    fprintf(f2, "%11.4E ", r.params[p1]);
	    p1++;
	  }
	  fprintf(f2, "\n");
	}
	for (t = 0; t < h.n_usr; t++) {
	  if (v) {
	    a = h.usr_egrid[t];
	    if (h.usr_egrid_type == 1) a += e;
	    a *= 1.0 + 0.5*FINE_STRUCTURE_CONST2 * a;
	    a = PI * AREA_AU20/(2.0*a);
	    if (!h.msub) a /= (mem_en_table[r.lower].j+1.0);
	    a *= r.strength[p2];
	    fprintf(f2, "%11.4E %11.4E %11.4E\n",
		    h.usr_egrid[t]*HARTREE_EV,
		    r.strength[p2], a);
	  } else {
	    fprintf(f2, "%11.4E %11.4E\n", h.usr_egrid[t], r.strength[p2]);
	  }
	  p2++;
	}
	if (k < r.nsub-1) {
	  fprintf(f2, "--------------------------------------------\n");
	}
      }      
      if (h.msub || h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
    nb += 1;
  }

  return nb;
}

int PrintRRTable(FILE *f1, FILE *f2, int v, int swp) {
  RR_HEADER h;
  RR_RECORD r;
  int n, i, t;
  int nb, k, m;
  float e, eph, ee, phi, rr;

  nb = 0;
  while (1) {
    n = ReadRRHeader(f1, &h, swp);
    if (n == 0) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "QKMODE\t= %d\n", h.qk_mode);
    fprintf(f2, "MULTIP\t= %d\n", h.multipole);
    fprintf(f2, "NPARAMS\t= %d\n", h.nparams);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);
    for (i = 0; i < h.n_tegrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }
    
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadRRRecord(f1, &r, swp, &h);
      if (n == 0) break;

      if (v) {
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "%6d %2d %6d %2d %11.4E %2d\n",
		r.b, mem_en_table[r.b].j, 
		r.f, mem_en_table[r.f].j,
		(e*HARTREE_EV), r.kl);
	
      } else {
	fprintf(f2, "%6d %6d %2d\n", r.b, r.f, r.kl);
      }
      
      if (h.qk_mode == QK_FIT) {
	for (t = 0; t < h.nparams; t++) {
	  if (v && t == h.nparams-1) {
	    fprintf(f2, "%11.4E ", r.params[t]*HARTREE_EV);
	  } else {
	    fprintf(f2, "%11.4E ", r.params[t]);
	  }
	}
	fprintf(f2, "\n");
      }
      
      for (t = 0; t < h.n_usr; t++) {
	if (v) {
	  if (h.usr_egrid_type == 0) {
	    eph = h.usr_egrid[t];
	    ee = eph - e;
	  } else {
	    ee = h.usr_egrid[t];
	    eph = ee + e;
	  }
	  phi = FINE_STRUCTURE_CONST2*ee;	  
	  phi = (1.0+phi)/(1.0+0.5*phi);
	  phi *= 2.0*PI*FINE_STRUCTURE_CONST*r.strength[t]*AREA_AU20;
	  rr = phi * pow(FINE_STRUCTURE_CONST*eph, 2) / (2.0*ee);
	  rr /= 1.0+0.5*FINE_STRUCTURE_CONST2*ee;
	  phi /= (mem_en_table[r.b].j + 1.0);
	  rr /= (mem_en_table[r.f].j + 1.0);
	  fprintf(f2, "%11.4E %11.4E %11.4E %11.4E\n",
		  h.usr_egrid[t]*HARTREE_EV, rr, phi, r.strength[t]);
	} else {
	  fprintf(f2, "%11.4E %11.4E\n", h.usr_egrid[t], r.strength[t]);
	}
      }
      if (h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }

    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);

    nb++;
  }

  return nb;
}

int PrintAITable(FILE *f1, FILE *f2, int v, int swp) {
  AI_HEADER h;
  AI_RECORD r;
  int n, i;
  int nb;
  float e, sdr;
  
  nb = 0;
  
  while (1) {
    n = ReadAIHeader(f1, &h, swp);
    if (n == 0) break;
 
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "CHANNE\t= %d\n", h.channel);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
       
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadAIRecord(f1, &r, swp);
      if (n == 0) break;
      if (v) {
	e = mem_en_table[r.b].energy - mem_en_table[r.f].energy;
	sdr = 0.5*(mem_en_table[r.b].j + 1.0);
	sdr *= PI*PI*r.rate/(e*(mem_en_table[r.f].j + 1.0));
	sdr *= AREA_AU20*HARTREE_EV;
	fprintf(f2, "%6d %2d %6d %2d %11.4E %11.4E %11.4E\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, (RATE_AU*r.rate), sdr);
      } else {
	fprintf(f2, "%6d %6d %15.8E\n", r.b, r.f, r.rate);
      }
    }
    
    free(h.egrid);
    nb++;
  }

  return nb;
}

int AIBranch(char *fn, int ib, int ia,
	     double *te, double *pa, double *ta) {
  F_HEADER fh;
  AI_HEADER h;
  AI_RECORD r;
  FILE *f;
  int n, i, k;
  double a, b, c, e;
  int swp;
    
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) {
    fclose(f);
    return 0;
  }
  if (fh.type != DB_AI) {
    printf("File type is not DB_AI\n");
    fclose(f);
    return -1;
  }
   
  a = 0.0;
  c = 0.0;
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadAIHeader(f, &h, swp);
    if (n == 0) break;
    for (k = 0; k < h.ntransitions; k++) {
      n = ReadAIRecord(f, &r, swp);
      if (n == 0) break;
      if (r.b == ib) {
	e = mem_en_table[r.b].energy - mem_en_table[r.f].energy;
	b = RATE_AU*r.rate;
	a += b;
	if (r.f == ia) {
	  c += b;
	}
      }
    }    
    free(h.egrid);
  }
  
  *pa = c;
  *ta = a;
  if (ia >= 0) {
    *te = mem_en_table[ib].energy - mem_en_table[ia].energy;
    *te *= HARTREE_EV;
  } else {
    *te = 0.0;
  }

  fclose(f);
  
  return 0;
}
  
int PrintAIMTable(FILE *f1, FILE *f2, int v, int swp) {
  AIM_HEADER h;
  AIM_RECORD r;
  int n, i, m;
  int nb;
  float e;
  double u = AREA_AU20*HARTREE_EV;
  
  nb = 0;
  while (1) {
    n = ReadAIMHeader(f1, &h, swp);
    if (n == 0) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "CHANNE\t= %d\n", h.channel);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadAIMRecord(f1, &r, swp);
      if (n == 0) break;
      if (v) {
	e = mem_en_table[r.b].energy - mem_en_table[r.f].energy;
	fprintf(f2, "%6d %2d %6d %2d %11.4E %2d\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, r.nsub);
	for (m = 0; m < r.nsub; m += 2) {
	  fprintf(f2, "%11.4E %11.4E\n", 
		  r.rate[m]*RATE_AU, r.rate[m+1]*u);
	}
      } else {
	fprintf(f2, "%6d %6d %2d\n", r.b, r.f, r.nsub);
	for (m = 0; m < r.nsub; m += 2) {
	  fprintf(f2, "%11.4E %11.4E\n", r.rate[m], r.rate[m+1]);
	}
      }
      free(r.rate);
    }

    free(h.egrid);
    nb++;
  }

  return nb;
}

int PrintCITable(FILE *f1, FILE *f2, int v, int swp) {
  CI_HEADER h;
  CI_RECORD r;
  int n, i, t;
  int nb, m;
  float e, a;

  nb = 0;
  while (1) {
    n = ReadCIHeader(f1, &h, swp);
    if (n == 0) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "QKMODE\t= %d\n", h.qk_mode);
    fprintf(f2, "NPARAMS\t= %d\n", h.nparams);
    fprintf(f2, "PWTYPE\t= %d\n", h.pw_type);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);
    for (i = 0; i < h.n_tegrid; i++) {      
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCIRecord(f1, &r, swp, &h);
      if (n == 0) break;
      
      if (v) {
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "%6d %2d %6d %2d %11.4E %2d\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, r.kl);
      } else {
	fprintf(f2, "%6d %6d %2d\n", r.b, r.f, r.kl);
      }
      
      for (t = 0; t < h.nparams; t++) {
	fprintf(f2, "%11.4E ", r.params[t]);
      }
      fprintf(f2, "\n");
      for (t = 0; t < h.n_usr; t++) {
	if (v) {
	  a = h.usr_egrid[t];
	  if (h.usr_egrid_type == 1) a += e;
	  a *= 1.0 + 0.5*FINE_STRUCTURE_CONST2*a;
	  a = AREA_AU20/(2.0*a*(mem_en_table[r.b].j + 1.0));
	  a *= r.strength[t];
	  fprintf(f2, "%11.4E %11.4E %11.4E\n",
		  h.usr_egrid[t]*HARTREE_EV, r.strength[t], a);
	} else {
	  fprintf(f2, "%11.4E %11.4E\n", h.usr_egrid[t], r.strength[t]);
	}
      }
      free(r.params); 
      free(r.strength);
    }
    
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
    
    nb++;
  }

  return nb;
}

int PrintCIMTable(FILE *f1, FILE *f2, int v, int swp) {
  CIM_HEADER h;
  CIM_RECORD r;
  int n, i, t, q;
  int nb, m, k;
  float e, a;

  nb = 0;
  while (1) {
    n = ReadCIMHeader(f1, &h, swp);
    if (n == 0) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCIMRecord(f1, &r, swp, &h);
      if (n == 0) break;
      
      if (v) {
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "%6d %2d %6d %2d %11.4E %2d\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, r.nsub);
      } else {
	fprintf(f2, "%6d %6d %2d\n", r.b, r.f, r.nsub);
      }

      if (v) {
	q = 0;
	for (k = 0; k < r.nsub; k ++) {
	  for (t = 0; t < h.n_usr; t++) {
	    a = h.usr_egrid[t];
	    if (h.usr_egrid_type == 1) a += e;
	    a *= 1.0 + 0.5*FINE_STRUCTURE_CONST2*a;
	    a = AREA_AU20/(2.0*a);
	    a *= r.strength[q];
	    fprintf(f2, "%11.4E %11.4E %11.4E\n",
		    h.usr_egrid[t]*HARTREE_EV, r.strength[q], a);
	    q++;
	  }
	  if (k < r.nsub-1) {
	    fprintf(f2, "--------------------------------------------\n");
	  }
	}
      } else {
	q = 0;
	for (k = 0; k < r.nsub; k++) {
	  for (t = 0; t < h.n_usr; t++) {
	    fprintf(f2, "%11.4E %11.4E\n", h.usr_egrid[t], r.strength[q]);
	    q++;
	  }
	}
      }
      free(r.strength);
    }
    
    free(h.egrid);
    free(h.usr_egrid);
    
    nb++;
  }

  return nb;
}

int PrintSPTable(FILE *f1, FILE *f2, int v, int swp) {
  SP_HEADER h;
  SP_RECORD r;
  SP_EXTRA rx;
  int n, i;
  int nb;
  float e, a;

  nb = 0;
  
  while (1) {
    n = ReadSPHeader(f1, &h, swp);
    if (n == 0) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "TYPE\t= %07d\n", h.type);
    fprintf(f2, "IBLK\t= %d\n", h.iblock);
    fprintf(f2, "ICOMP\t= %s\n", h.icomplex);
    fprintf(f2, "FBLK\t= %d\n", h.fblock);
    fprintf(f2, "FCOMP\t= %s\n", h.fcomplex);
    
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      e = r.energy;
      if (v) e *= HARTREE_EV;
      a = r.strength;
      if (iuta) {
	fprintf(f2, "%6d %6d %13.6E %11.4E %11.4E %11.4E %11.4E\n", 
		r.upper, r.lower, e, rx.sdev*HARTREE_EV, a, r.rrate, r.trate);
      } else {
	fprintf(f2, "%6d %6d %13.6E %11.4E\n", 
		r.upper, r.lower, e, a, r.rrate, r.trate);
      }
    }
    nb += 1;
  }
  
  return nb;
}

double IonDensity(char *fn, int k) {
  F_HEADER fh;
  SP_HEADER h;
  SP_RECORD r;
  SP_EXTRA rx;
  int i;
  FILE *f;
  int n, swp;
  double d;

  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1.0;
  }
  
  n = ReadFHeader(f, &fh, &swp);
  
  if (fh.type != DB_SP) {
    printf("file not of type DB_SP\n");
    fclose(f);
    return -1.0;
  }

  d = 0.0;
  while (1) {
    n = ReadSPHeader(f, &h, swp);
    if (n == 0) break;
    if (h.type != 0 || h.nele != k) {
      fseek(f, h.length, SEEK_CUR);
      continue;
    }
    
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f, &r, &rx, swp);
      if (n == 0) break;
      d += r.strength;
    }
  }
 
  fclose(f);

  return d;
}
		     
int PrintRTTable(FILE *f1, FILE *f2, int v, int swp) {
  RT_HEADER h;
  RT_RECORD r;
  int n, i;
  int nb, nele;
  double dc, re, ea, ci, rr, pi;

  nb = 0;
  nele = -1;
  while (1) {
    n = ReadRTHeader(f1, &h, swp);
    if (n == 0) break;
    if (h.nele != nele) {
      if (nele != -1) {
	fprintf(f2, "\n");
	fprintf(f2, " SUM  %10.4E %10.4E %10.4E %10.4E %10.4E %10.4E %3d\n",
		rr, dc, re, pi, ea, ci, nele);
	fprintf(f2, "\n");
      }
      nele = h.nele;
      dc = 0.0;
      re = 0.0;
      ea = 0.0;
      rr = 0.0;
      pi = 0.0;
      ci = 0.0;
    }
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "IBLK\t= %d\n", h.iblock);
    fprintf(f2, "ILEV\t= %d\n", h.ilev);
    fprintf(f2, "ICOMP\t= %s\n", h.icomplex);
    fprintf(f2, "EDEN\t= %15.8E\n", h.eden);
    fprintf(f2, "EDIST\t= %d\n", h.iedist);
    fprintf(f2, "NPEDIS\t= %d\n", h.np_edist);
    for (i = 0; i < h.np_edist; i++) {
      fprintf(f2, "\t %15.8E\n", h.p_edist[i]);
    }
    fprintf(f2, "PDEN\t= %15.8E\n", h.pden);
    fprintf(f2, "PDIST\t= %d\n", h.ipdist);
    fprintf(f2, "NPPDIS\t= %d\n", h.np_pdist);
    for (i = 0; i < h.np_pdist; i++) {
      fprintf(f2, "\t %15.8E\n", h.p_pdist[i]);
    }
    free(h.p_edist);
    free(h.p_pdist);

    fprintf(f2, "DENS\t= %15.8E\n", h.nb);
    fprintf(f2,"         NB         TR         CE");
    fprintf(f2, "         RR         AI         CI\n");
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadRTRecord(f1, &r, swp);
      if (n == 0) break;
      fprintf(f2, "%4d  %10.4E %10.4E %10.4E %10.4E %10.4E %10.4E %s\n",
	      r.iblock, r.nb, r.tr, r.ce, r.rr, r.ai, r.ci, r.icomplex);
      if (r.iblock == -1) {
	dc += r.ai;
	rr += r.rr;
      } else if (r.iblock == -2) {
	re += r.ai;
      } else if (r.iblock == -3) {
	pi += r.rr;
	ci += r.ci;
	ea += r.ai;
      }
    }
    nb += 1;
  }
  fprintf(f2, "\n SUM  %10.4E %10.4E %10.4E %10.4E %10.4E %10.4E %3d\n\n",
	  rr, dc, re, pi, ea, ci, nele);
  return nb;
}

int PrintDRTable(FILE *f1, FILE *f2, int v, int swp) {
  DR_HEADER h;
  DR_RECORD r;
  int n, i;
  int nb;
  double e, e1;
  
  nb = 0;
  while (1) {
    n = ReadDRHeader(f1, &h, swp);
    if (n == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "ILEV\t= %d\n", h.ilev);
    e = h.energy;
    if (v) e *= HARTREE_EV;
    fprintf(f2, "E\t= %15.8E\n", e);
    fprintf(f2, "JLEV\t= %d\n", h.j);
    fprintf(f2, "NREC\t= %d\n", h.vn);
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadDRRecord(f1, &r, swp);
      if (n == 0) break;
      e = r.energy;
      e1 = r.etrans;
      if (v) {
	e *= HARTREE_EV;
	e1 *= HARTREE_EV;
	fprintf(f2, "%6d %2d %4d %2d %3d %4d %3d %2d %2d %10.4E %10.4E %10.4E %10.4E %10.4E\n",
		r.ilev, r.j, h.ilev, h.j, r.ibase, r.flev, r.fbase, 
		h.vn, r.vl, e, e1, r.ai, r.total_rate, r.br);
      } else {
	fprintf(f2, "%6d %2d %3d %4d %3d %2d %10.4E %10.4E %10.4E %10.4E %10.4E\n",
		r.ilev, r.j, r.ibase, r.flev, r.fbase, r.vl, 
		e, e1, r.ai, r.total_rate, r.br);
      }
    }
    nb++;
  }

  return nb;
}

int AppendTable(char *fn) {
  F_HEADER fh;
  FILE *f;
  int n, swp;
    
  f = fopen(fn, "r");
  if (f == NULL) return -1;
  n = ReadFHeader(f, &fh, &swp);
  if (swp) {
    printf("File %s is in different byte-order\n");
    fclose(f);
    return -1;
  }
  memcpy(&(fheader[fh.type-1]), &fh, sizeof(F_HEADER));
  fclose(f);
  
  return 0;
}
  
int JoinTable(char *fn1, char *fn2, char *fn) {
  F_HEADER fh1, fh2;
  FILE *f1, *f2, *f;
  int n, swp1, swp2;
#define NBUF 8192
  char buf[NBUF];

  f1 = fopen(fn1, "r");
  if (f1 == NULL) return -1;
  f2 = fopen(fn2, "r");
  if (f2 == NULL) return -1;

  n = ReadFHeader(f1, &fh1, &swp1);
  if (n == 0) {
    fclose(f1);
    fclose(f2);
    return 0;
  }
  n = ReadFHeader(f2, &fh2, &swp2);
  if (n == 0) {
    fclose(f1);
    fclose(f2);
    return 0;
  }
  if (swp1 != swp2) {
    printf("Files %s and %s have different byte-order\n", fn1, fn2);
    return -1;
  }
  if (fh1.type != fh2.type) {
    printf("Files %s and %s are of different type\n", fn1, fn2);
    return -1;
  }
  if (fh1.atom != fh2.atom) {
    printf("Files %s and %s are for different element\n", fn1, fn2);
    return -1;
  }

  f = fopen(fn, "w");
  if (f == NULL) return -1;
  fh1.nblocks += fh2.nblocks;
  
  WriteFHeader(f, &fh1);
  while (1) {
    n = fread(buf, 1, NBUF, f1);
    if (n > 0) {
      if (n > fwrite(buf, 1, n, f)) {
	printf("write error\n");
	return -1;
      }
    }
    if (n < NBUF) break;
  }
  while (1) {
    n = fread(buf, 1, NBUF, f2);
    if (n > 0) {
      if (n > fwrite(buf, 1, n, f)) {
	printf("write error\n");
	return -1;
      }
    }
    if (n < NBUF) break;
  }

  fclose(f1);
  fclose(f2);
  fclose(f);
  
  return 0;
#undef NBUF
}
