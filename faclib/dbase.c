#include "dbase.h"
#include "cf77.h"

static char *rcsid="$Id: dbase.c,v 1.71 2005/01/15 01:23:35 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

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
static int itrf = 0;

void SetTRF(int m) {
  itrf = m;
}

void SetUTA(int m) {
  iuta = m;
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
  SwapEndian((char *) &(r->ibase), sizeof(short));
  SwapEndian((char *) &(r->fbase), sizeof(short));
  SwapEndian((char *) &(r->vl), sizeof(short));
  SwapEndian((char *) &(r->j), sizeof(short));
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

int WriteFHeader(FILE *f, F_HEADER *fh) {
  int n;

  n = fwrite(fh, sizeof(F_HEADER), 1, f);
  if (n != 1) return 0;
  return sizeof(F_HEADER);
}

int ReadFHeader(FILE *f, F_HEADER *fh, int *swp) {
  int n;

  n = fread(fh, sizeof(F_HEADER), 1, f);
  if (n != 1) return 0;
  *swp = 0;
  if (CheckEndian(fh) != (int) (fheader[0].symbol[3])) {
    *swp = 1;
    SwapEndianFHeader(fh);
  }
  if (fh->type == DB_TR && itrf >= 0) {
    if (VersionLE(fh, 1, 0, 6)) itrf = 1;
    else itrf = 0;
  }
  return sizeof(F_HEADER);
}

int WriteENHeader(FILE *f, EN_HEADER *h) {
  int n;
  
  n = fwrite(h, sizeof(EN_HEADER), 1, f);
  if (n != 1) return 0;

  return sizeof(EN_HEADER);
}

int WriteTRHeader(FILE *f, TR_HEADER *h) {
  int n;
  
  n = fwrite(h, sizeof(TR_HEADER), 1, f);
  if (n != 1) return 0;
  
  return sizeof(TR_HEADER);
}

int WriteCEHeader(FILE *f, CE_HEADER *h) {
  int n, m;

  n = fwrite(h, sizeof(CE_HEADER), 1, f);
  if (n != 1) return 0;
  m = sizeof(CE_HEADER);
  n = fwrite(h->tegrid, sizeof(double), h->n_tegrid, f);
  if (n != h->n_tegrid) return 0;
  m += sizeof(double)*h->n_tegrid;
  n = fwrite(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) return 0;
  m += sizeof(double)*h->n_egrid;
  n = fwrite(h->usr_egrid, sizeof(double), h->n_usr, f);
  if (n != h->n_usr) return 0;
  m += sizeof(double)*h->n_usr;
  
  return m;
}

int WriteRRHeader(FILE *f, RR_HEADER *h) {
  int n, m;

  n = fwrite(h, sizeof(RR_HEADER), 1, f);
  if (n != 1) return 0;
  m = sizeof(RR_HEADER);
  n = fwrite(h->tegrid, sizeof(double), h->n_tegrid, f);
  if (n != h->n_tegrid) return 0;
  m += sizeof(double)*h->n_tegrid;
  n = fwrite(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) return 0;
  m += sizeof(double)*h->n_egrid;
  n = fwrite(h->usr_egrid, sizeof(double), h->n_usr, f);
  if (n != h->n_usr) return 0;
  m += sizeof(double)*h->n_usr;
  
  return m;
}

int WriteAIHeader(FILE *f, AI_HEADER *h) {
  int n, m;
  
  n = fwrite(h, sizeof(AI_HEADER), 1, f);
  if (n != 1) return 0;
  m = sizeof(AI_HEADER);
  n = fwrite(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) return 0;
  m += sizeof(double)*h->n_egrid;
  
  return m;
}

int WriteAIMHeader(FILE *f, AIM_HEADER *h) {
  int n, m;
  
  n = fwrite(h, sizeof(AIM_HEADER), 1, f);
  if (n != 1) return 0;
  m = sizeof(AIM_HEADER);
  n = fwrite(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) return 0;
  m += sizeof(double)*h->n_egrid;
  
  return m;
}

int WriteCIHeader(FILE *f, CI_HEADER *h) {
  int n, m;

  n = fwrite(h, sizeof(CI_HEADER), 1, f);
  if (n != 1) return 0;
  m = sizeof(CI_HEADER);
  n = fwrite(h->tegrid, sizeof(double), h->n_tegrid, f);
  if (n != h->n_tegrid) return 0;
  m += sizeof(double)*h->n_tegrid;
  n = fwrite(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) return 0;
  m += sizeof(double)*h->n_egrid;
  n = fwrite(h->usr_egrid, sizeof(double), h->n_usr, f);
  if (n != h->n_usr) return 0;
  m += sizeof(double)*h->n_usr;
  
  return m;
}

int WriteCIMHeader(FILE *f, CIM_HEADER *h) {
  int n, m;
  
  n = fwrite(h, sizeof(CIM_HEADER), 1, f);
  if (n != 1) return 0;
  m = sizeof(CIM_HEADER);
  n = fwrite(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) return 0;
  m += sizeof(double)*h->n_egrid;
  n = fwrite(h->usr_egrid, sizeof(double), h->n_usr, f);
  if (n != h->n_usr) return 0;
  m += sizeof(double)*h->n_usr;
  return m;
}

int WriteSPHeader(FILE *f, SP_HEADER *h) {
  int n;
  
  n = fwrite(h, sizeof(SP_HEADER), 1, f);
  if (n != 1) return 0;
  
  return sizeof(SP_HEADER);
}

int WriteRTHeader(FILE *f, RT_HEADER *h) {
  int n, m;
  
  n = fwrite(h, sizeof(RT_HEADER), 1, f);
  if (n != 1) return 0;
  m = sizeof(RT_HEADER);
  n = fwrite(h->p_edist, sizeof(double), h->np_edist, f);
  if (n != h->np_edist) return 0;
  m += sizeof(double)*h->np_edist;
  n = fwrite(h->p_pdist, sizeof(double), h->np_pdist, f);
  if (n != h->np_pdist) return 0;
  m += sizeof(double)*h->np_pdist;
  
  return m;
}

int WriteDRHeader(FILE *f, DR_HEADER *h) {
  int n;

  n = fwrite(h, sizeof(DR_HEADER), 1, f);
  if (n != 1) return 0;
  
  return sizeof(DR_HEADER);
}

int WriteENRecord(FILE *f, EN_RECORD *r) {
  int n;

  if (en_header.length == 0) {
    fheader[DB_EN-1].nblocks++;
    n = WriteENHeader(f, &en_header);
  }
  n = fwrite(r, sizeof(EN_RECORD), 1, f);
  if (n != 1) return 0;
  en_header.nlevels += 1;
  en_header.length += sizeof(EN_RECORD);
  return sizeof(EN_RECORD);
}

int WriteTRRecord(FILE *f, TR_RECORD *r, TR_EXTRA *rx) {
  int n;

  if (tr_header.length == 0) {
    fheader[DB_TR-1].nblocks++;
    n = WriteTRHeader(f, &tr_header);
  }
  n = fwrite(r, sizeof(TR_RECORD), 1, f);  
  if (n != 1) return 0;
  tr_header.ntransitions += 1;
  tr_header.length += sizeof(TR_RECORD);
  if (iuta) {
    n = fwrite(rx, sizeof(TR_EXTRA), 1, f);
    if (n != 1) return 0;
    tr_header.length += sizeof(TR_EXTRA);
  }
  return sizeof(TR_RECORD);
}

int WriteCERecord(FILE *f, CE_RECORD *r) {
  int n;
  int m, m0;

  if (ce_header.length == 0) {
    fheader[DB_CE-1].nblocks++;
    n = WriteCEHeader(f, &ce_header);
  }
  m = sizeof(CE_RECORD);
  n = fwrite(r, m, 1, f);
  if (n != 1) return 0;
  ce_header.ntransitions += 1;
  ce_header.length += m;
  m0 = m;
  if (ce_header.msub) {
    m = r->nsub;
    ce_header.length += sizeof(float)*m;
    n = fwrite(r->params, sizeof(float), m, f);
    if (n != m) return 0;
    m0 += sizeof(float)*m;
  } else if (ce_header.qk_mode == QK_FIT) {
    m = ce_header.nparams * r->nsub;
    ce_header.length += sizeof(float)*m;
    n = fwrite(r->params, sizeof(float), m, f);
    if (n != m) return 0;
    m0 += sizeof(float)*m;
  }
  m = ce_header.n_usr * r->nsub;
  ce_header.length += sizeof(float)*m;
  n = fwrite(r->strength, sizeof(float), m, f);
  if (n != m) return 0;
  m0 += sizeof(float)*m;

  return m0;
}

int WriteRRRecord(FILE *f, RR_RECORD *r) {
  int n;
  int m, m0;

  if (rr_header.length == 0) {
    fheader[DB_RR-1].nblocks++;
    n = WriteRRHeader(f, &rr_header);
  }
  rr_header.ntransitions += 1;
  m = sizeof(RR_RECORD);
  rr_header.length += m;
  n = fwrite(r, m, 1, f);
  if (n != 1) return 0;
  m0 = m;
  if (rr_header.qk_mode == QK_FIT) {
    m = rr_header.nparams;
    rr_header.length += sizeof(float)*m;
    n = fwrite(r->params, sizeof(float), m, f);
    if (n != m) return 0;
    m0 += sizeof(float)*m;
  }
  m = rr_header.n_usr;
  rr_header.length += sizeof(float)*m;
  n = fwrite(r->strength, sizeof(float), m, f);
  if (n != m) return 0;
  m0 += sizeof(float)*m;

  return m0;
}

int WriteAIRecord(FILE *f, AI_RECORD *r) {
  int n;

  if (ai_header.length == 0) {
    fheader[DB_AI-1].nblocks++;
    WriteAIHeader(f, &ai_header);
  }
  ai_header.ntransitions += 1;
  ai_header.length += sizeof(AI_RECORD);
  n = fwrite(r, sizeof(AI_RECORD), 1, f);
  if (n != 1) return 0;

  return sizeof(AI_RECORD);
}

int WriteAIMRecord(FILE *f, AIM_RECORD *r) {
  int n;

  if (aim_header.length == 0) {
    fheader[DB_AIM-1].nblocks++;
    WriteAIMHeader(f, &aim_header);
  }
  aim_header.ntransitions += 1;
  aim_header.length += sizeof(AIM_RECORD);
  n = fwrite(r, sizeof(AIM_RECORD), 1, f);
  if (n != 1) return 0;
  aim_header.length += sizeof(float)*r->nsub;
  n = fwrite(r->rate, sizeof(float), r->nsub, f);
  if (n != r->nsub) return 0;

  return sizeof(AIM_RECORD);
}

int WriteCIRecord(FILE *f, CI_RECORD *r) {
  int n;
  int m, m0;

  if (ci_header.length == 0) {
    fheader[DB_CI-1].nblocks++;
    WriteCIHeader(f, &ci_header);
  }
  ci_header.ntransitions += 1;
  m = sizeof(CI_RECORD);
  ci_header.length += m;
  n = fwrite(r, m, 1, f);
  if (n != 1) return 0;
  m0 = m;
  m = ci_header.nparams;
  ci_header.length += sizeof(float)*m;
  n = fwrite(r->params, sizeof(float), m, f);
  if (n != m) return 0;
  m0 += sizeof(float)*m;
  m = ci_header.n_usr;
  ci_header.length += sizeof(float)*m;
  n = fwrite(r->strength, sizeof(float), m, f);
  if (n != m) return 0;
  m0 += sizeof(float)*m;

  return m0;
}

int WriteCIMRecord(FILE *f, CIM_RECORD *r) {
  int n;
  int m, m0;

  if (cim_header.length == 0) {
    fheader[DB_CIM-1].nblocks++;
    WriteCIMHeader(f, &cim_header);
  }
  cim_header.ntransitions += 1;
  m = sizeof(CIM_RECORD);
  cim_header.length += m;
  n = fwrite(r, m, 1, f);
  if (n != 1) return 0;
  m0 = m;
  m = r->nsub*cim_header.n_usr;
  cim_header.length += sizeof(float)*m;
  n = fwrite(r->strength, sizeof(float), m, f);
  if (n != m) return 0;
  m0 += sizeof(float)*m;
  return m0;
}

int WriteSPRecord(FILE *f, SP_RECORD *r, SP_EXTRA *rx) {
  int n;

  if (sp_header.length == 0) {
    fheader[DB_SP-1].nblocks++;
    WriteSPHeader(f, &sp_header);
  }
  sp_header.ntransitions += 1;
  sp_header.length += sizeof(SP_RECORD);
  n = fwrite(r, sizeof(SP_RECORD), 1, f);  
  if (n != 1) return 0;
  if (iuta) {
    sp_header.length += sizeof(SP_EXTRA);
    n = fwrite(rx, sizeof(SP_EXTRA), 1, f);
    if (n != 1) return 0;
  }

  return sizeof(SP_RECORD);
}

int WriteRTRecord(FILE *f, RT_RECORD *r) {
  int n;

  if (rt_header.length == 0) {
    fheader[DB_RT-1].nblocks++;
    WriteRTHeader(f, &rt_header);
  }
  rt_header.ntransitions += 1;
  rt_header.length += sizeof(RT_RECORD);
  n = fwrite(r, sizeof(RT_RECORD), 1, f);  
  if (n != 1) return 0;

  return sizeof(RT_RECORD);
}

int WriteDRRecord(FILE *f, DR_RECORD *r) {
  int n;

  if (dr_header.length == 0) {
    fheader[DB_DR-1].nblocks++;
    WriteDRHeader(f, &dr_header);
  }
  dr_header.ntransitions += 1;
  dr_header.length += sizeof(DR_RECORD);
  n = fwrite(r, sizeof(DR_RECORD), 1, f);
  if (n != 1) return 0;
  
  return sizeof(DR_RECORD);
}

int ReadENHeader(FILE *f, EN_HEADER *h, int swp) {
  int n;

  n = fread(h, sizeof(EN_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianENHeader(h);
  
  return sizeof(EN_HEADER);
}
  
int ReadENRecord(FILE *f, EN_RECORD *r, int swp) {
  int n;

  n = fread(r, sizeof(EN_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianENRecord(r);
  
  return sizeof(EN_RECORD);
}

int ReadTRHeader(FILE *f, TR_HEADER *h, int swp) {
  int n;

  n = fread(h, sizeof(TR_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianTRHeader(h);
  if (h->length/h->ntransitions > sizeof(TR_RECORD)) iuta = 1;
  else iuta =0;
  return sizeof(TR_HEADER);
}

int ReadTRRecord(FILE *f, TR_RECORD *r, TR_EXTRA *rx, int swp) {
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

int ReadCEHeader(FILE *f, CE_HEADER *h, int swp) {
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

int ReadCERecord(FILE *f, CE_RECORD *r, int swp, CE_HEADER *h) {
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

int ReadRRHeader(FILE *f, RR_HEADER *h, int swp) {
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

int ReadRRRecord(FILE *f, RR_RECORD *r, int swp, RR_HEADER *h) {
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

int ReadAIHeader(FILE *f, AI_HEADER *h, int swp) {
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

int ReadAIMHeader(FILE *f, AIM_HEADER *h, int swp) {
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

int ReadAIRecord(FILE *f, AI_RECORD *r, int swp) {
  int n;

  n = fread(r, sizeof(AI_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianAIRecord(r);
  
  return sizeof(AI_RECORD);
}

int ReadAIMRecord(FILE *f, AIM_RECORD *r, int swp) {
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

int ReadCIHeader(FILE *f, CI_HEADER *h, int swp) {
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

int ReadCIRecord(FILE *f, CI_RECORD *r, int swp, CI_HEADER *h) {
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

int ReadCIMHeader(FILE *f, CIM_HEADER *h, int swp) {
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

int ReadCIMRecord(FILE *f, CIM_RECORD *r, int swp, CIM_HEADER *h) {
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

int ReadSPHeader(FILE *f, SP_HEADER *h, int swp) {
  int n;

  n = fread(h, sizeof(SP_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianSPHeader(h);
  if (h->length/h->ntransitions > sizeof(SP_RECORD)) iuta = 1;
  else iuta = 0;
  return sizeof(SP_HEADER);
}

int ReadSPRecord(FILE *f, SP_RECORD *r, SP_EXTRA *rx, int swp) {
  int n;

  n = fread(r, sizeof(SP_RECORD), 1, f);
  if (n != 1) return 0;
  if (iuta) {
    n = fread(rx, sizeof(SP_EXTRA), 1, f);
    if (n != 1) return 0;
  }
  if (swp) SwapEndianSPRecord(r, rx);
  return sizeof(SP_RECORD);
}

int ReadRTHeader(FILE *f, RT_HEADER *h, int swp) {
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

int ReadRTRecord(FILE *f, RT_RECORD *r, int swp) {
  int n;

  n = fread(r, sizeof(RT_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianRTRecord(r);
  return sizeof(RT_RECORD);
}

int ReadDRHeader(FILE *f, DR_HEADER *h, int swp) {
  int n;
  
  n = fread(h, sizeof(DR_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianDRHeader(h);
  
  return sizeof(DR_HEADER);
}

int ReadDRRecord(FILE *f, DR_RECORD *r, int swp) {
  int n;

  n = fread(r, sizeof(DR_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianDRRecord(r);
  
  return sizeof(DR_RECORD);
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
    if (f == NULL) return NULL;
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
  size_t n;
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

void PrepCECrossHeader(CE_HEADER *h, double *data) {
  double *eusr, *x;
  int m, m1, j;

  eusr = h->usr_egrid;
  m = h->n_usr;
  m1 = m + 1;
  x = data+2+m1;
  data[0] = h->te0*HARTREE_EV;
  for (j = 0; j < m; j++) {
    x[j] = log((h->te0 + eusr[j])/h->te0);
  }
  x[m] = eusr[m-1]/(h->te0+eusr[m-1]);
}

void PrepCECrossRecord(int k, CE_RECORD *r, CE_HEADER *h, 
		       double *data) {
  double *eusr, *x, *y, *w, e;
  float *cs;
  int m, m1, j, t;
  int j1, j2, t1, t2;

  eusr = h->usr_egrid;
  m = h->n_usr;
  m1 = m + 1;
  y = data + 2;
  x = y + m1;
  w = x + m1;
  e = mem_en_table[r->upper].energy - mem_en_table[r->lower].energy;
  data[1] = r->bethe;
  
  cs = r->strength;
  if (k == 0) {
    if (h->msub) {
      j1 = mem_en_table[r->lower].j;
      j2 = mem_en_table[r->upper].j;
      for (j = 0; j < m; j++) {
	y[j] = 0.0;
      }
      t = 0;
      for (t1 = -j1; t1 <= 0; t1 += 2) {
	for (t2 = -j2; t2 <= j2; t2 += 2) {
	  for (j = 0; j < m; j++) {
	    y[j] += cs[t];
	    if (t1 != 0) y[j] += cs[t];
	    t++;
	  }
	}
      }
    } else {
      for (j = 0; j < m; j++) {	
	y[j] = cs[j];
      }
    }
    y[m] = r->born[0];
  }

  if (h->msub) {
    for (j = 0; j < m; j++) {
      if (y[j]) {
	w[j] = cs[k*m+j]/y[j];
      } else {
	w[j] = 1.0;
      }
    }
    if (r->bethe < 0) {
      w[j] = w[j-1];
    } else {
      w[j] = r->params[k];
    }
  }
}

double InterpolateCECross(double e, CE_RECORD *r, CE_HEADER *h, 
			  double *data, double *ratio) {
  double *x, *y, *w;
  int m, m1, n, one;
  double a, b, x0, y0, eth, e0, c, d;

  eth = mem_en_table[r->upper].energy - mem_en_table[r->lower].energy;
  eth = eth * HARTREE_EV;

  a = 0.0;
  *ratio = 1.0;

  if (e < 0.0) return a;

  m = h->n_usr;
  m1 = m + 1;
  x0 = log((data[0]+e)/data[0]);
  y = data + 2;
  x = y + m1;
  w = x + m1;
  
  if (x0 < x[m-1]) {
    n = 2;
    one = 1;
    UVIP3P(n, m, x, y, one, &x0, &a);
    if (h->msub) {
      UVIP3P(n, m, x, w, one, &x0, &b);
      if (b < 0.0) b = 0.0;
      a *= b;
      *ratio = b;
    }
  } else {
    x0 = e/(data[0] + e);
    y0 = y[m-1];
    if (data[1] > 0) {
      e0 = ((x[m]*data[0]/(1.0-x[m]))+eth)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth) - c/(1.0+c);
      y0 /= 1.0 + c;
      y0 -= data[1]*b;
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
      e0 = (e + eth)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth) - c/(1.0+c);
      a += data[1]*b;
      a *= 1.0 + c;
    } else if (data[1]+1.0 == 1.0) {
      e0 = ((x[m]*data[0]/(1.0-x[m]))+eth)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      y0 /= 1.0 + c;
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
      e0 = (e + eth)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      a *= 1.0 + c;
    } else {
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
    }
    if (h->msub) {
      b = w[m] + (x0-1.0)*(w[m-1]-w[m])/(x[m]-1.0);
      a *= b;
      *ratio = b;
    }
  }

  return a;
}

int CECross(char *ifn, char *ofn, int i0, int i1, 
	    int negy, double *egy, int mp) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;
  CE_HEADER h;
  CE_RECORD r;
  int i, t, m, k;
  double data[2+(1+MAXNUSR)*3], e, cs, a, ratio;
  double eth, a1, cs1, k2, rp, e1, e0;
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }
  
  f1 = fopen(ifn, "r");
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }
    
  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    fclose(f1);
    fclose(f2);
    return 0;  
  }
  
  if (fh.type != DB_CE || fh.nblocks == 0) {
    printf("File %s is not of DB_CE type\n", ifn);
    fclose(f1);
    fclose(f2);
    return 0;
  }
   
  while (1) {
    n = ReadCEHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f1, &r, swp, &h);
      if (r.lower == i0 && r.upper == i1) {
	PrepCECrossHeader(&h, data);
	eth = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	e = eth*HARTREE_EV;
	fprintf(f2, "#%5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\t%d\n",
		r.lower, mem_en_table[r.lower].j,
		r.upper, mem_en_table[r.upper].j,
		e, negy, r.nsub);
	for (k = 0; k < r.nsub; k++) {
	  PrepCECrossRecord(k, &r, &h, data);
	  for (t = 0; t < negy; t++) {
	    if (mp == 0) {
	      e0 = egy[t];
	      e1 = e0 - e;
	    } else {
	      e1 = egy[t];
	      e0 = e1 + e;
	    }
	    if (e1 > 0) {
	      cs = InterpolateCECross(e1, &r, &h, data, &ratio);
	      a = e0/HARTREE_EV;
	      a = a*(1.0+0.5*FINE_STRUCTURE_CONST2*a);
	      a = PI*AREA_AU20/(2.0*a);
	      if (!h.msub) a /= (mem_en_table[r.lower].j+1.0);
	      a *= cs;
	      if (data[1] > 0.0) {	      
		cs1 = data[1]*log(e0/e) + r.born[0];
		k2 = e0/HARTREE_EV;
		k2 = 2.0*k2*(1.0+0.5*FINE_STRUCTURE_CONST2*k2);
		a1 = FINE_STRUCTURE_CONST2*k2;
		a1 = a1/(1.0+a1);
		a1 = data[1]*(log(0.5*k2/eth) - a1);
		a1 += r.born[0];
		a1 *= 1.0 + FINE_STRUCTURE_CONST2*k2;
		k2 = cs1/a1;
		rp = k2*(1.0+0.5*FINE_STRUCTURE_CONST2*e0/HARTREE_EV);
		if (rp > 1.0) {
		  k2 /= rp;
		  rp = 1.0;
		}
		cs1 = cs*k2;
		a1 = a*rp;
	      } else {
		k2 = (egy[t]+e)/HARTREE_EV;
		k2 = 2.0*k2*(1.0+0.5*FINE_STRUCTURE_CONST2*k2);
		k2 = 1.0 + FINE_STRUCTURE_CONST2*k2;
		cs1 = cs/k2;
		a1 = a*(1.0+0.5*FINE_STRUCTURE_CONST2*e0/HARTREE_EV)/k2;
	      }
	    } else {
	      cs = 0.0;
	      a = 0.0;
	      cs1 = 0.0;
	      a1 = 0.0;
	      ratio = 0.0;
	    }
	    fprintf(f2, "%11.4E %11.4E %11.4E %11.4E %11.4E %11.4E\n",
		    e0, cs, a, cs1, a1, ratio);
	  }
	  fprintf(f2, "\n\n");
	}
	goto DONE;
      }
    }
  }

 DONE:
  fclose(f1);
  fclose(f2);
  return 0;
}

int CEMaxwell(char *ifn, char *ofn, int i0, int i1, 
	      int nt, double *temp) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;
  CE_HEADER h;
  CE_RECORD r;
  int i, t, m, k, p;
  double data[2+(1+MAXNUSR)*4], e, cs, a, b, c, ratio;
  double xg[15] = {.933078120172818E-01,
		   .492691740301883E+00,
		   .121559541207095E+01,
		   .226994952620374E+01,
		   .366762272175144E+01,
		   .542533662741355E+01,
		   .756591622661307E+01,
		   .101202285680191E+02,
		   .131302824821757E+02,
		   .166544077083300E+02,
		   .207764788994488E+02,
		   .256238942267288E+02,
		   .314075191697539E+02,
		   .385306833064860E+02,
		   .480260855726858E+02};
  double wg[15] = {.218234885940086E+00,
		   .342210177922884E+00,
		   .263027577941681E+00,
		   .126425818105931E+00,
		   .402068649210010E-01,
		   .856387780361184E-02,
		   .121243614721425E-02,
		   .111674392344251E-03,
		   .645992676202287E-05,
		   .222631690709627E-06,
		   .422743038497938E-08,
		   .392189726704110E-10,
		   .145651526407313E-12,
		   .148302705111330E-15,
		   .160059490621113E-19};
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }
  
  f1 = fopen(ifn, "r");
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }
    
  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    fclose(f1);
    fclose(f2);
    return 0;  
  }
  
  if (fh.type != DB_CE || fh.nblocks == 0) {
    printf("File %s is not of DB_CE type\n", ifn);
    fclose(f1);
    fclose(f2);
    return 0;
  }
   
  while (1) {
    n = ReadCEHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f1, &r, swp, &h);
      if (r.lower == i0 && r.upper == i1) {
	PrepCECrossHeader(&h, data);
	e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	e *= HARTREE_EV;
	fprintf(f2, "#%5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\t%d\n",
		r.lower, mem_en_table[r.lower].j,
		r.upper, mem_en_table[r.upper].j,
		e, nt, r.nsub);
	for (k = 0; k < r.nsub; k++) {
	  PrepCECrossRecord(k, &r, &h, data);
	  for (t = 0; t < nt; t++) {
	    cs = 0.0;
	    for (p = 0; p < 15; p++) {
	      a = temp[t]*xg[p];
	      b = a/HARTREE_EV;
	      a = InterpolateCECross(a, &r, &h, data, &ratio);
	      c = 2.0*b*(1.0+0.5*FINE_STRUCTURE_CONST2*b);
	      c = 1.0+FINE_STRUCTURE_CONST2*c;
	      c = c*(1.0+0.5*FINE_STRUCTURE_CONST2*b);
	      c = sqrt(c);
	      cs += wg[p]*a/c;
	    }
	    a = 217.16*sqrt(HARTREE_EV/(2.0*temp[t]));
	    a *= cs*exp(-e/temp[t]);
	    if (!h.msub) a /= (mem_en_table[r.lower].j+1.0);
	    fprintf(f2, "%11.4E\t%11.4E\t%11.4E\n", 
		    temp[t], cs, a);
	  }
	  fprintf(f2, "\n\n");
	}
	goto DONE;
      }
    }
  }

 DONE:
  fclose(f1);
  fclose(f2);
  return 0;
}

int TotalCICross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int imin, int imax) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;
  CI_HEADER h;
  CI_RECORD r;
  int i, t, nb, m;
  double *c, tc, a, b, x, e;
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  f1 = fopen(ifn, "r");
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }
  
  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    fclose(f1);
    fclose(f2);
    return 0;  
  }

  if (fh.type != DB_CI || fh.nblocks == 0) {
    printf("File %s is not of DB_CI type\n", ifn);
    fclose(f1);
    fclose(f2);
    return 0;
  }
  
  c = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    c[i] = 0.0;
  }

  nb = 0;
  
  if (imin < 0) imin = 0;
  if (imax < 0) imax = mem_en_table_size - 1;

  while (1) {
    n = ReadCIHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCIRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (r.b != ilev) continue;
      if (r.f < imin || r.f > imax) continue;
      e = mem_en_table[r.f].energy - mem_en_table[r.b].energy; 
      
      for (t = 0; t < negy; t++) {
	if (egy[t] < e) continue;
	x = egy[t]/e;
	a = 1.0/x;
	b = 1.0 - a;
	tc = r.params[0]*log(x) + r.params[1]*b*b;
	tc += r.params[2]*a*b + r.params[3]*a*a*b;
	a = egy[t]*(1.0 + FINE_STRUCTURE_CONST2*egy[t]);
	tc *= AREA_AU20/(2.0*a*(mem_en_table[r.b].j + 1.0));
	c[t] += tc;
      }
      free(r.params);
      free(r.strength);
    }

    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
    
    nb++;
  }

  fprintf(f2, "#Energy (eV)   CI Cross (10^-20 cm2)\n");
  for (t = 0; t < negy; t++) {
    fprintf(f2, " %11.4E    %15.8E\n", egy[t]*HARTREE_EV, c[t]);
  }

  free(c);  

  fclose(f1);
  fclose(f2);
  return nb;
} 

int TotalPICross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int imin, int imax) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;
  RR_HEADER h;
  RR_RECORD r;
  int i, t, nb, m;
  float e, eph, ee, phi;
  double *xusr, *dstrength, *c, tc, emax;
  double x, y, a;
  int np=3, one=1, nele;
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  f1 = fopen(ifn, "r");
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    fclose(f1);
    fclose(f2);
    return 0;  
  }

  if (fh.type != DB_RR || fh.nblocks == 0) {
    printf("File %s is not of DB_RR type\n", ifn);
    fclose(f1);
    fclose(f2);
    return 0;
  }
  
  c = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    c[i] = 0.0;
  }

  nb = 0;
  
  if (imin < 0) imin = 0;
  if (imax < 0) imax = mem_en_table_size - 1;

  while(1) {
    n = ReadRRHeader(f1, &h, swp);
    if (n == 0) break;
    nele = h.nele;
    xusr = (double *) malloc(sizeof(double)*h.n_usr); 
    dstrength = (double *) malloc(sizeof(double)*h.n_usr);
    emax = h.usr_egrid[h.n_usr-1];
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadRRRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (r.b != ilev) continue;
      if (r.f < imin || r.f > imax) continue;
      e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;

      for (t = 0; t < h.n_usr; t++) {
	dstrength[t] = log(r.strength[t]);
	xusr[t] = log(1.0 + h.usr_egrid[t]/e);
      }
      
      for (t = 0; t < negy; t++) {
	eph = egy[t];
	ee = eph - e;
	if (ee <= 0.0) continue;
	a = FINE_STRUCTURE_CONST2*ee;
	a = (1.0+a)/(1+0.5*a);
	if (h.qk_mode != QK_FIT || ee <= emax) {
	  x = log(eph/e);
	  UVIP3P(np, h.n_usr, xusr, dstrength, one, &x, &tc);
	  tc = exp(tc);
	} else {
	  x = (ee + r.params[3])/r.params[3];
	  y = (1 + r.params[2])/(sqrt(x) + r.params[2]);
	  tc = (-3.5 - r.kl + 0.5*r.params[1])*log(x) + r.params[1]*log(y);
	  if (r.params[0] > 0.0) {
	    tc = tc + log(r.params[0]*(eph/(ee+r.params[3])));
	    tc = exp(tc);
	  } else {
	    tc = 0.0;
	  }
	}
	phi = a*2.0*PI*FINE_STRUCTURE_CONST*tc*AREA_AU20;	
	c[t] += phi/(mem_en_table[r.b].j + 1.0);
      }
      if (h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }

    free(dstrength);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);      
    free(xusr);
    
    nb++;
  }

  fprintf(f2, "# Energy (eV)   PI Cross (10^-20 cm2)\n");
  for (t = 0; t < negy; t++) {
    fprintf(f2, " %12.5E    %15.8E\n", egy[t]*HARTREE_EV, c[t]);
  }

  free(c);  

  fclose(f1);
  fclose(f2);
  return nb;
}
  
int TotalRRCross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int n0, int n1, int nmax,
		 int imin, int imax) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;
  RR_HEADER h;
  RR_RECORD r;
  int i, t, nb, m;
  float e, eph, ee, phi, rr;
  double *xusr, *dstrength, *c, tc, emax;
  double x, y, a;
  int np=3, one=1, nele;

  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  f1 = fopen(ifn, "r");
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    fclose(f1);
    fclose(f2);
    return 0;  
  }

  if (fh.type != DB_RR || fh.nblocks == 0) {
    printf("File %s is not of DB_RR type\n", ifn);
    fclose(f1);
    fclose(f2);
    return 0;
  }
  
  c = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    c[i] = 0.0;
  }

  nb = 0;
  
  if (imin < 0) imin = 0;
  if (imax < 0) imax = mem_en_table_size - 1;

  while(1) {
    n = ReadRRHeader(f1, &h, swp);
    if (n == 0) break;
    nele = h.nele;
    xusr = (double *) malloc(sizeof(double)*h.n_usr); 
    dstrength = (double *) malloc(sizeof(double)*h.n_usr);
    emax = h.usr_egrid[h.n_usr-1];
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadRRRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (r.f != ilev) continue;
      if (r.b < imin || r.b > imax) continue;
      e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
      
      for (t = 0; t < h.n_usr; t++) {
	dstrength[t] = log(r.strength[t]);
	xusr[t] = log(1.0 + h.usr_egrid[t]/e);
      }
      
      for (t = 0; t < negy; t++) {
	ee = egy[t];
	a = FINE_STRUCTURE_CONST2*ee;
	a = (1.0+a)/(1+0.5*a);
	eph = ee + e;
	if (h.qk_mode != QK_FIT || ee <= emax) {
	  x = log(eph/e);
	  UVIP3P(np, h.n_usr, xusr, dstrength, one, &x, &tc);
	  tc = exp(tc);
	} else {
	  x = (ee + r.params[3])/r.params[3];
	  y = (1 + r.params[2])/(sqrt(x) + r.params[2]);
	  tc = (-3.5 - r.kl + 0.5*r.params[1])*log(x) + r.params[1]*log(y);
	  if (r.params[0] > 0.0) {
	    tc = tc + log(r.params[0]*(eph/(ee+r.params[3])));
	    tc = exp(tc);
	  } else {
	    tc = 0.0;
	  }
	}
	phi = 2.0*PI*FINE_STRUCTURE_CONST*tc*AREA_AU20;
	rr = a * phi * pow(FINE_STRUCTURE_CONST*eph, 2) / (2.0*ee);
	rr /= (mem_en_table[r.f].j + 1.0);
	c[t] += rr;
      }
      if (h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }  

    free(dstrength);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);      
    free(xusr);
    
    nb++;
  }
  
  x = fh.atom - nele + 1.0;
  for (n = n0+1; n < n1; n++) {
    for (t = 0; t < negy; t++) {
      c[t] += AREA_AU20*RRCrossHn(x, egy[t], n);
    }
  }
  for (n = n1+1; n <= nmax; n++) {
    for (t = 0; t < negy; t++) {
      c[t] += AREA_AU20*RRCrossHn(x, egy[t], n);
    }
  }
  
  fprintf(f2, "#Energy (eV)   RR Cross (10^-20 cm2)\n");
  for (t = 0; t < negy; t++) {
    fprintf(f2, " %11.4E    %15.8E\n", egy[t]*HARTREE_EV, c[t]);
  }

  free(c);

  fclose(f1);
  fclose(f2);
  return nb;
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
    fclose(f1);
    fclose(f2);
    return 0;  
  }

  if (v && (fh.type < DB_SP || fh.type > DB_DR)) {
    if (mem_en_table == NULL) {
      printf("Energy table has not been built in memory.\n");
      fclose(f1);
      fclose(f2);
      return -1;
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

  fclose(f1);
  if (f2 != stdout) fclose(f2);
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

  k = ilev;
  nlevels = 0;
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    nlevels += h.nlevels;
    if (k < h.nlevels) {
      if (k > 0) fseek(f, sizeof(EN_RECORD)*k, SEEK_CUR);
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
}

int MemENTable(char *fn) {
  F_HEADER fh;  
  EN_HEADER h;
  EN_RECORD r;
  FILE *f;
  char *s;
  int n, i, nlevels;
  float e0;
  int swp;

  f = fopen(fn, "r");
  if (f == NULL) return -1;

  n = ReadFHeader(f, &fh, &swp);  
  if (n == 0) return 0;
  if (fh.type != DB_EN) return -1;

  if (mem_en_table) free(mem_en_table);

  nlevels = 0;
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    n = sizeof(EN_RECORD);
    if (h.length > n) {
      fseek(f, h.length-n, SEEK_CUR);
    }
    n = ReadENRecord(f, &r, swp);
    if (r.ilev >= nlevels) nlevels = r.ilev+1;
  }

  mem_en_table = (EN_SRECORD *) malloc(sizeof(EN_SRECORD)*nlevels);
  mem_en_table_size = nlevels;

  e0 = 0.0;
  fseek(f, sizeof(F_HEADER), SEEK_SET);
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
      mem_en_table[r.ilev].j = r.j;
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
  int p, vnl;

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
      fprintf(f2, "%6d %6d %15.8E %1d %5d %4d %-20s %-20s %-s\n",
	      r.ilev, r.ibase, e, p, vnl, r.j, r.ncomplex, r.sname, r.name);
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
	fprintf(f2, "%6d %6d %13.6E %11.4E %11.4E\n", 
		r.upper, r.lower, e, rx.sdev*HARTREE_EV, a);
      } else {
	fprintf(f2, "%6d %6d %13.6E %11.4E\n", r.upper, r.lower, e, a);
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
