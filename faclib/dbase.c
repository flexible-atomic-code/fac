#include "dbase.h"

static char *rcsid="$Id: dbase.c,v 1.31 2002/11/13 22:28:26 mfgu Exp $";
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
static CI_HEADER ci_header;
static SP_HEADER sp_header;
static RT_HEADER rt_header;
static DR_HEADER dr_header;

static EN_SRECORD *mem_en_table = NULL;
static int iground;

void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);

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

int SwapEndianTRRecord(TR_RECORD *r) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->strength), sizeof(float));
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

int SwapEndianAIRecord(AI_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->rate), sizeof(float));
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

int SwapEndianSPRecord(SP_RECORD *r) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->energy), sizeof(float));
  SwapEndian((char *) &(r->strength), sizeof(float));
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
  SwapEndian((char *) &(r->ibase), sizeof(int));
  SwapEndian((char *) &(r->vl), sizeof(short));
  SwapEndian((char *) &(r->j), sizeof(short));
  SwapEndian((char *) &(r->energy), sizeof(float));
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
  iground = 0;

  return 0;
}

int ReinitDBase(int m) {
  int i;

  if (m < 0) return 0;
  if (mem_en_table) {
    free(mem_en_table);
    mem_en_table = NULL;
  }
  if (m == 0) {
    return InitDBase();
  } else {
    iground = 0;
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

FILE *OpenFile(char *fn, F_HEADER *fhdr) {
  int ihdr;
  FILE *f;

  ihdr = fhdr->type - 1;

  f = fopen(fn, "r+");
  if (f == NULL) {
    if (fheader[ihdr].nblocks != 0) {
      printf("File %s written in this session has been erased\n", fn);
      printf("Exitting");
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

  fheader[ihdr].type = fhdr->type;
  strncpy(fheader[ihdr].symbol, fhdr->symbol, 2);
  fheader[ihdr].atom = fhdr->atom;
  fwrite(&(fheader[ihdr]), sizeof(F_HEADER), 1, f);
  return f;
}

int CloseFile(FILE *f, F_HEADER *fhdr) {
  int ihdr;
 
  ihdr = fhdr->type-1;
  fseek(f, 0, SEEK_SET);
  fheader[ihdr].type = fhdr->type;
  fwrite(&(fheader[ihdr]), sizeof(F_HEADER), 1, f);
  
  fclose(f);
  return 0;
}

int InitFile(FILE *f, F_HEADER *fhdr, void *rhdr) {
  EN_HEADER *en_hdr;
  TR_HEADER *tr_hdr;
  CE_HEADER *ce_hdr;
  RR_HEADER *rr_hdr;
  AI_HEADER *ai_hdr;
  CI_HEADER *ci_hdr;
  SP_HEADER *sp_hdr;
  RT_HEADER *rt_hdr;
  DR_HEADER *dr_hdr;
  long int p;
  size_t n;
  int ihdr;

  if (f == NULL) return 0;

  ihdr = fhdr->type - 1;
  fheader[ihdr].nblocks += 1;
  fseek(f, 0, SEEK_END);
  p = ftell(f);

  switch (fhdr->type) {
  case DB_EN:
    en_hdr = (EN_HEADER *) rhdr;
    en_header.position = p;
    en_header.length = 0;
    en_header.nele = en_hdr->nele;
    en_header.nlevels = 0;
    n = fwrite(&en_header, sizeof(EN_HEADER), 1, f);
    break;
  case DB_TR:
    tr_hdr = (TR_HEADER *) rhdr;
    memcpy(&tr_header, tr_hdr, sizeof(TR_HEADER));
    tr_header.position = p;
    tr_header.length = 0;
    tr_header.ntransitions = 0;
    n = fwrite(&tr_header, sizeof(TR_HEADER), 1, f);
    break;
  case DB_CE:
    ce_hdr = (CE_HEADER *) rhdr;
    memcpy(&ce_header, ce_hdr, sizeof(CE_HEADER));
    ce_header.position = p;
    ce_header.length = 0;
    ce_header.ntransitions = 0;
    n = fwrite(&ce_header, sizeof(CE_HEADER), 1, f);
    n = fwrite(ce_header.tegrid, sizeof(double), ce_header.n_tegrid, f);
    n = fwrite(ce_header.egrid, sizeof(double), ce_header.n_egrid, f);
    n = fwrite(ce_header.usr_egrid, sizeof(double), ce_header.n_usr, f);
    break;
  case DB_RR:
    rr_hdr = (RR_HEADER *) rhdr;
    memcpy(&rr_header, rr_hdr, sizeof(RR_HEADER));
    rr_header.position = p;
    rr_header.length = 0;
    rr_header.ntransitions = 0;
    n = fwrite(&rr_header, sizeof(RR_HEADER), 1, f);
    n = fwrite(rr_header.tegrid, sizeof(double), rr_header.n_tegrid, f);
    n = fwrite(rr_header.egrid, sizeof(double), rr_header.n_egrid, f);
    n = fwrite(rr_header.usr_egrid, sizeof(double), rr_header.n_usr, f);
    break;
  case DB_AI:
    ai_hdr = (AI_HEADER *) rhdr;
    memcpy(&ai_header, ai_hdr, sizeof(AI_HEADER));
    ai_header.position = p;
    ai_header.length = 0;
    ai_header.ntransitions = 0;
    n = fwrite(&ai_header, sizeof(AI_HEADER), 1, f);
    n = fwrite(ai_header.egrid, sizeof(double), ai_header.n_egrid, f);
    break;
  case DB_CI:    
    ci_hdr = (CI_HEADER *) rhdr;
    memcpy(&ci_header, ci_hdr, sizeof(CI_HEADER));
    ci_header.position = p;
    ci_header.length = 0;
    ci_header.ntransitions = 0;
    n = fwrite(&ci_header, sizeof(CI_HEADER), 1, f);
    n = fwrite(ci_header.tegrid, sizeof(double), ci_header.n_tegrid, f);
    n = fwrite(ci_header.egrid, sizeof(double), ci_header.n_egrid, f);
    n = fwrite(ci_header.usr_egrid, sizeof(double), ci_header.n_usr, f);
    break;
  case DB_SP:
    sp_hdr = (SP_HEADER *) rhdr;
    memcpy(&sp_header, sp_hdr, sizeof(SP_HEADER));
    sp_header.position = p;
    sp_header.length = 0;
    sp_header.ntransitions = 0;
    n = fwrite(&sp_header, sizeof(SP_HEADER), 1, f);
    break;
  case DB_RT:
    rt_hdr = (RT_HEADER *) rhdr;
    memcpy(&rt_header, rt_hdr, sizeof(RT_HEADER));
    rt_header.position = p;
    rt_header.length = 0;
    rt_header.ntransitions = 0;
    n = fwrite(&rt_header, sizeof(RT_HEADER), 1, f);
    n = fwrite(rt_header.p_edist, sizeof(double), rt_header.np_edist, f);
    n = fwrite(rt_header.p_pdist, sizeof(double), rt_header.np_pdist, f);
    break;
  case DB_DR:
    dr_hdr = (DR_HEADER *) rhdr;
    memcpy(&dr_header, dr_hdr, sizeof(DR_HEADER));
    dr_header.position = p;
    dr_header.length = 0;
    dr_header.ntransitions = 0;
    n = fwrite(&dr_header, sizeof(DR_HEADER), 1, f);
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
    n = fwrite(&en_header, sizeof(EN_HEADER), 1, f);
    break;
  case DB_TR:
    fseek(f, tr_header.position, SEEK_SET);
    n = fwrite(&tr_header, sizeof(TR_HEADER), 1, f);
    break;
  case DB_CE:
    fseek(f, ce_header.position, SEEK_SET);
    n = fwrite(&ce_header, sizeof(CE_HEADER), 1, f);
    n = fwrite(ce_header.tegrid, sizeof(double), ce_header.n_tegrid, f);
    n = fwrite(ce_header.egrid, sizeof(double), ce_header.n_egrid, f);
    n = fwrite(ce_header.usr_egrid, sizeof(double), ce_header.n_usr, f);
    break;
  case DB_RR:
    fseek(f, rr_header.position, SEEK_SET);
    n = fwrite(&rr_header, sizeof(RR_HEADER), 1, f);
    n = fwrite(rr_header.tegrid, sizeof(double), rr_header.n_tegrid, f);
    n = fwrite(rr_header.egrid, sizeof(double), rr_header.n_egrid, f);
    n = fwrite(rr_header.usr_egrid, sizeof(double), rr_header.n_usr, f);
    break;
  case DB_AI:
    fseek(f, ai_header.position, SEEK_SET);
    n = fwrite(&ai_header, sizeof(AI_HEADER), 1, f);
    n = fwrite(ai_header.egrid, sizeof(double), ai_header.n_egrid, f);
    break;
  case DB_CI:
    fseek(f, ci_header.position, SEEK_SET);
    n = fwrite(&ci_header, sizeof(CI_HEADER), 1, f);
    n = fwrite(ci_header.tegrid, sizeof(double), ci_header.n_tegrid, f);
    n = fwrite(ci_header.egrid, sizeof(double), ci_header.n_egrid, f);
    n = fwrite(ci_header.usr_egrid, sizeof(double), ci_header.n_usr, f);
    break;
  case DB_SP:
    fseek(f, sp_header.position, SEEK_SET);
    n = fwrite(&sp_header, sizeof(SP_HEADER), 1, f);
    break;
  case DB_RT:
    fseek(f, rt_header.position, SEEK_SET);
    n = fwrite(&rt_header, sizeof(RT_HEADER), 1, f);
    break;
  case DB_DR:
    fseek(f, dr_header.position, SEEK_SET);
    n = fwrite(&dr_header, sizeof(DR_HEADER), 1, f);
    break;
  default:
    break;
  }
  return 0;
}

int TotalRRCross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;
  RR_HEADER h;
  RR_RECORD r;
  int i, t, nb, m, k;
  float *params, *strength;
  float e, eph, ee, phi, rr;
  double *dstrength, *c, tc, emax;
  double x, y;
  int np=3, one=1;

  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  f1 = fopen(ifn, "r");
  if (f1 == NULL) return -1;

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) return -1;

  n = fread(&fh, sizeof(F_HEADER), 1, f1);
  if (n != 1) {
    fclose(f1);
    fclose(f2);
    return 0;  
  }

  swp = 0;
  if (CheckEndian(&fh) != (int) (fheader[0].symbol[3])) {
    swp = 1;
    SwapEndianFHeader(&fh);
  }

  if (fh.type != DB_RR || fh.nblocks == 0) {
    fclose(f1);
    fclose(f2);
    return 0;
  }
  
  c = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    c[i] = 0.0;
  }

  nb = 0;
  
  while(1) {
    n = fread(&h, sizeof(RR_HEADER), 1, f1);
    if (n != 1) break;
    if (swp) SwapEndianRRHeader(&h);
    h.tegrid = (double *) malloc(sizeof(double)*h.n_tegrid);
    n = fread(h.tegrid, sizeof(double), h.n_tegrid, f1);
    h.egrid = (double *) malloc(sizeof(double)*h.n_egrid);
    n = fread(h.egrid, sizeof(double), h.n_egrid, f1);
    h.usr_egrid = (double *) malloc(sizeof(double)*h.n_usr);
    n = fread(h.usr_egrid, sizeof(double), h.n_usr, f1);
    if (swp) {
      for (i = 0; i < h.n_tegrid; i++) {
	SwapEndian((char *) &(h.tegrid[i]), sizeof(double));
      }
      for (i = 0; i < h.n_egrid; i++) {
	SwapEndian((char *) &(h.egrid[i]), sizeof(double));
      }
      for (i = 0; i < h.n_usr; i++) {
	SwapEndian((char *) &(h.usr_egrid[i]), sizeof(double));
      }
    }      
    if (h.qk_mode == QK_FIT) {
      m = h.nparams;
      params = (float *) malloc(sizeof(float)*m);
    }
    m = h.n_usr;
    strength = (float *) malloc(sizeof(float)*m);
    dstrength = (double *) malloc(sizeof(double)*m);
    emax = h.usr_egrid[h.n_usr-1];
    for (t = 0; t < h.n_usr; t++) {
      h.usr_egrid[t] = log(h.usr_egrid[t]);
    }

    for (i = 0; i < h.ntransitions; i++) {
      n = fread(&r, sizeof(RR_RECORD), 1, f1);
      if (swp) SwapEndianRRRecord(&r);
      if (h.qk_mode == QK_FIT) {
	m = h.nparams;
	n = fread(params, sizeof(float), m, f1);
      }
      m = h.n_usr;
      n = fread(strength, sizeof(float), m, f1);	
      if (r.f != ilev) continue;
      e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
      
      if (swp) {
	if (h.qk_mode == QK_FIT) {
	  for (t = 0; t < h.nparams; t++) {
	    SwapEndian((char *) &(params[t]), sizeof(float));
	  }
	}
	for (t = 0; t < h.n_usr; t++) {
	  SwapEndian((char *) &(strength[t]), sizeof(float));
	}
      }
      
      for (t = 0; t < h.n_usr; t++) {
	dstrength[t] = strength[t];
      }
      
      for (t = 0; t < negy; t++) {
	ee = egy[t];
	if (ee <= emax) {
	  x = log(ee);
	  uvip3p_(&np, &(h.n_usr), h.usr_egrid, dstrength, &one, &x, &tc);
	  eph = ee + e;
	} else {
	  eph = ee + e;
	  x = (ee + params[3])/params[3];
	  y = (1 + params[2])/(sqrt(x) + params[2]);
	  tc = (-3.5 - r.kl + 0.5*params[1])*log(x) + params[1]*log(y);
	  tc = params[0]*exp(tc)*(eph/(ee+params[3]));
	}
	phi = 2.0*PI*FINE_STRUCTURE_CONST*tc*AREA_AU20;
	rr = phi * pow(FINE_STRUCTURE_CONST*eph, 2) / (2.0*ee);
	rr /= (mem_en_table[r.f].j + 1.0);
	c[t] += rr;
      }
    }  

    if (h.qk_mode == QK_FIT) free(params);
    free(strength);
    free(dstrength);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);      
    
    nb++;
  }

  fprintf(f2, "Energy (eV)   RR Cross (10^-20 cm2)\n");
  for (t = 0; t < negy; t++) {
    fprintf(f2, "%11.4E    %15.8E\n", egy[t]*HARTREE_EV, c[t]);
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

  n = fread(&fh, sizeof(F_HEADER), 1, f1);
  if (n != 1) {
    fclose(f1);
    fclose(f2);
    return 0;  
  }

  swp = 0;
  if (CheckEndian(&fh) != (int) (fheader[0].symbol[3])) {
    swp = 1;
    SwapEndianFHeader(&fh);
  }

  if (fh.nblocks == 0) {
    fclose(f1);
    fclose(f2);
    return 0;
  }

  if (v && fh.type < DB_SP) {
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
  default:
    break;
  }

  fclose(f1);
  if (f2 != stdout) fclose(f2);
  return n;
}

int WriteENRecord(FILE *f, EN_RECORD *r) {
  int n;
  en_header.nlevels += 1;
  en_header.length += sizeof(EN_RECORD);
  n = fwrite(r, sizeof(EN_RECORD), 1, f);  
  return n;
}

int FreeMemENTable(void) {
  if (mem_en_table) free(mem_en_table);
  mem_en_table = NULL;
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
  int n, i, k;
  int swp;
  
  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  n = fread(&fh, sizeof(F_HEADER), 1, f);
  if (n != 1) {
    fclose(f);
    return 0;
  }
  if (CheckEndian(&fh) != (int) (fheader[0].symbol[3])) {
    swp = 1;
    SwapEndianFHeader(&fh);
  } else {
    swp = 0;
  }
  if (fh.type != DB_EN) {
    printf("File type is not DB_EN\n");
    fclose(f);
    return -1;
  }

  for (i = 0; i < fh.nblocks; i++) {
    n = fread(&h, sizeof(EN_HEADER), 1, f);
    if (swp) {
      SwapEndianENHeader(&h);
    }
    if (h.nele != nele) {
      fseek(f, h.length, SEEK_CUR);
      continue;
    }
    for (k = 0; k < h.nlevels; k++) {
      n = fread(&r, sizeof(EN_RECORD), 1, f);
      if (swp) SwapEndianENRecord(&r);
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
  n = fread(&fh, sizeof(F_HEADER), 1, f);
  if (n != 1) {
    fclose(f);
    return 0;
  }
  if (CheckEndian(&fh) != (int) (fheader[0].symbol[3])) {
    swp = 1;
    SwapEndianFHeader(&fh);
  } else {
    swp = 0;
  }
  if (fh.type != DB_EN) {
    printf("File type is not DB_EN\n");
    fclose(f);
    return -1;
  }
  k = ilev;
  nlevels = 0;
  for (i = 0; i < fh.nblocks; i++) {
    n = fread(&h, sizeof(EN_HEADER), 1, f);
    if (swp) {
      SwapEndianENHeader(&h);
    }
    nlevels += h.nlevels;
    if (k < h.nlevels) {
      if (k > 0) fseek(f, sizeof(EN_RECORD)*k, SEEK_CUR);
      n = fread(&r, sizeof(EN_RECORD), 1, f);
      if (swp) {
	SwapEndianENRecord(&r);
      }	
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
  if (i == fh.nblocks) return nlevels;

  return 0;
}

int MemENTable(char *fn) {
  F_HEADER fh;  
  EN_HEADER h;
  EN_RECORD r;
  FILE *f;
  int n, i, nlevels;
  float e0;
  int swp;

  f = fopen(fn, "r");
  if (f == NULL) return -1;
  
  n = fread(&fh, sizeof(F_HEADER), 1, f);
  if (n != 1) return 0;
  if (CheckEndian(&fh) != (int) (fheader[0].symbol[3])) {
    swp = 1;
    SwapEndianFHeader(&fh);
  } else {
    swp = 0;
  }
  if (fh.type != DB_EN) return -1;

  if (mem_en_table) free(mem_en_table);  
  nlevels = 0;
  for (i = 0; i < fh.nblocks; i++) {
    n = fread(&h, sizeof(EN_HEADER), 1, f);
    if (swp) {
      SwapEndianENHeader(&h);
    }
    nlevels += h.nlevels;
    fseek(f, h.length, SEEK_CUR);
  }
  mem_en_table = (EN_SRECORD *) malloc(sizeof(EN_SRECORD)*nlevels);
    
  e0 = 0.0;
  fseek(f, sizeof(F_HEADER), SEEK_SET);
  while (1) {
    n = fread(&h, sizeof(EN_HEADER), 1, f);
    if (swp) {
      SwapEndianENHeader(&h);
    }
    if (n != 1) break;
    for (i = 0; i < h.nlevels; i++) {
      n = fread(&r, sizeof(EN_RECORD), 1, f);
      if (swp) {
	SwapEndianENRecord(&r);
      }
      if (n != 1) break;
      if (r.energy < e0) {
	e0 = r.energy;
	iground = r.ilev;
      }
      mem_en_table[r.ilev].energy = r.energy;
      mem_en_table[r.ilev].p = r.p;
      mem_en_table[r.ilev].j = r.j;
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
    n = fread(&h, sizeof(EN_HEADER), 1, f1);
    if (n != 1) break;
    if (swp) SwapEndianENHeader(&h);
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NLEV\t= %d\n", h.nlevels);
    fprintf(f2, "         Energy       P   VNL 2J\n");
    for (i = 0; i < h.nlevels; i++) {
      n = fread(&r, sizeof(EN_RECORD), 1, f1);
      if (swp) SwapEndianENRecord(&r);
      if (n != 1) break;
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
      fprintf(f2, "%5d %15.8E %1d %5d %2d %-20s %-20s %-s\n",
	      r.ilev, e, p, vnl, r.j, r.ncomplex, r.sname, r.name);
    }
    nb += 1;
  }
  
  return nb;
}

int WriteTRRecord(FILE *f, TR_RECORD *r) {
  int n;
  tr_header.ntransitions += 1;
  tr_header.length += sizeof(TR_RECORD);
  n = fwrite(r, sizeof(TR_RECORD), 1, f);  
  return n;
}

int PrintTRTable(FILE *f1, FILE *f2, int v, int swp) {
  TR_HEADER h;
  TR_RECORD r;
  int n, i;
  int nb;
  float e, a;

  nb = 0;
  
  while (1) {
    n = fread(&h, sizeof(TR_HEADER), 1, f1);
    if (n != 1) break;
    if (swp) SwapEndianTRHeader(&h);
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "Multip\t= %d\n", (int)h.multipole);
    fprintf(f2, "Gauge\t= %d\n", (int)h.gauge);
    fprintf(f2, "Mode\t= %d\n", (int)h.mode);

    for (i = 0; i < h.ntransitions; i++) {
      n = fread(&r, sizeof(TR_RECORD), 1, f1);
      if (swp) SwapEndianTRRecord(&r);
      if (v) {
	e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	a = 2.0*pow((FINE_STRUCTURE_CONST*e),2)*FINE_STRUCTURE_CONST;
	a *= r.strength/(mem_en_table[r.upper].j + 1.0);
	a *= RATE_AU;
	fprintf(f2, "%5d\t%2d\t%5d\t%2d\t%11.4E %15.8E %15.8E\n",
		r.upper, mem_en_table[r.upper].j,
		r.lower, mem_en_table[r.lower].j,
		(e*HARTREE_EV), r.strength, a);
      } else {
	fprintf(f2, "%5d\t%5d\t%15.8E\n", 
		r.upper, r.lower, r.strength);
      }
    }
    nb += 1;
  }

  return nb;
}

int WriteCERecord(FILE *f, CE_RECORD *r) {
  int n;
  int m;

  ce_header.ntransitions += 1;
  m = sizeof(CE_RECORD);
  ce_header.length += m;
  n = fwrite(r, m, 1, f);
  if (ce_header.qk_mode == QK_FIT) {
    m = ce_header.nparams * r->nsub;
    ce_header.length += sizeof(float)*m;
    n = fwrite(r->params, sizeof(float), m, f);
  }
  m = ce_header.n_usr * r->nsub;
  ce_header.length += sizeof(float)*m;
  n = fwrite(r->strength, sizeof(float), m, f);

  return n;
}
  
int PrintCETable(FILE *f1, FILE *f2, int v, int swp) {
  CE_HEADER h;
  CE_RECORD r;
  int n, i, t;
  int nb, nsub;
  int m, k, p1, p2;
  float *params, *strength;
  float a, e;

  nb = 0;
 
  while (1) {
    n = fread(&h, sizeof(CE_HEADER), 1, f1);
    if (n != 1) break;
    if (swp) SwapEndianCEHeader(&h);

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "QKMODE\t= %d\n", h.qk_mode);
    fprintf(f2, "NPARAMS\t= %d\n", h.nparams);
    fprintf(f2, "MSUB\t= %d\n", h.msub);
    fprintf(f2, "PWTYPE\t= %d\n", h.pw_type);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);
    h.tegrid = (double *) malloc(sizeof(double)*h.n_tegrid);
    n = fread(h.tegrid, sizeof(double), h.n_tegrid, f1);
    for (i = 0; i < h.n_tegrid; i++) {
      if (swp) SwapEndian((char *) &(h.tegrid[i]), sizeof(double));
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "TE0\t= %15.8E\n", h.te0 * HARTREE_EV);
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    h.egrid = (double *) malloc(sizeof(double)*h.n_egrid);
    n = fread(h.egrid, sizeof(double), h.n_egrid, f1);
    for (i = 0; i < h.n_egrid; i++) {
      if (swp) SwapEndian((char *) &(h.egrid[i]), sizeof(double));
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    h.usr_egrid = (double *) malloc(sizeof(double)*h.n_usr);
    n = fread(h.usr_egrid, sizeof(double), h.n_usr, f1);
    for (i = 0; i < h.n_usr; i++) {
      if (swp) SwapEndian((char *) &(h.usr_egrid[i]), sizeof(double));
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }

    nsub = 1;
    if (h.qk_mode == QK_FIT) {
      m = h.nparams * nsub;
      params = (float *) malloc(sizeof(float)*m);
    }
    m = h.n_usr * nsub;
    strength = (float *) malloc(sizeof(float)*m);
    for (i = 0; i < h.ntransitions; i++) {
      n = fread(&r, sizeof(CE_RECORD), 1, f1);
      if (swp) SwapEndianCERecord(&r);
      if (r.nsub > nsub) {
	if (h.qk_mode == QK_FIT) {
	  m = h.nparams * r.nsub;
	  params = (float *) realloc(params, sizeof(float)*m);
	}
	m = h.n_usr * r.nsub;
	strength = (float *) realloc(strength, sizeof(float)*m);
	nsub = r.nsub;
      }

      if (h.qk_mode == QK_FIT) {
	m = h.nparams * r.nsub;
	n = fread(params, sizeof(float), m, f1);
      }
      m = h.n_usr * r.nsub;
      n = fread(strength, sizeof(float), m, f1);
      if (v) {
	e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	fprintf(f2, "%5d\t%2d\t%5d\t%2d\t%11.4E\t%d\n",
		r.lower, mem_en_table[r.lower].j,
		r.upper, mem_en_table[r.upper].j,
		e*HARTREE_EV, r.nsub);
	fprintf(f2, "%11.4E %11.4E %11.4E\n", 
		r.bethe, r.born[0], r.born[1]);
      } else {
	fprintf(f2, "%5d\t%5d\t%d\n", 
		r.lower, r.upper, r.nsub);
	fprintf(f2, "%11.4E %11.4E %11.4E\n", 
		r.bethe, r.born[0], r.born[1]);
      }
      
      p1 = 0;
      p2 = 0;
      for (k = 0; k < r.nsub; k++) {
	if (h.qk_mode == QK_FIT) {
	  for (t = 0; t < h.nparams; t++) {
	    if (swp) SwapEndian((char *) &(params[p1]), sizeof(float));
	    fprintf(f2, "%11.4E ", params[p1]);
	    p1++;
	  }
	  fprintf(f2, "\n");
	}
	for (t = 0; t < h.n_usr; t++) {
	  if (swp) SwapEndian((char *) &(strength[p2]), sizeof(float));
	  if (v) {
	    a = h.usr_egrid[t];
	    if (h.usr_egrid_type == 1) a += e;
	    a *= 1.0 + 0.5*FINE_STRUCTURE_CONST2 * a;
	    a = PI * AREA_AU20/(2.0*a*(mem_en_table[r.lower].j+1.0));
	    a *= strength[p2];
	    fprintf(f2, "%11.4E\t%11.4E\t%11.4E\n",
		    h.usr_egrid[t]*HARTREE_EV,
		    strength[p2], a);
	  } else {
	    fprintf(f2, "%11.4E\t%11.4E\n", h.usr_egrid[t], strength[p2]);
	  }
	  p2++;
	}
	if (k < r.nsub-1) {
	  fprintf(f2, "--------------------------------------------\n");
	}
      }      
    }
    if (h.qk_mode == QK_FIT) free(params);
    free(strength);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
    nb += 1;
  }

  return nb;
}

int WriteRRRecord(FILE *f, RR_RECORD *r) {
  int n;
  int m;

  rr_header.ntransitions += 1;
  m = sizeof(RR_RECORD);
  rr_header.length += m;
  n = fwrite(r, m, 1, f);
  if (rr_header.qk_mode == QK_FIT) {
    m = rr_header.nparams;
    rr_header.length += sizeof(float)*m;
    n = fwrite(r->params, sizeof(float), m, f);
  }
  
  m = rr_header.n_usr;
  rr_header.length += sizeof(float)*m;
  n = fwrite(r->strength, sizeof(float), m, f);
  
  return n;
}

int PrintRRTable(FILE *f1, FILE *f2, int v, int swp) {
  RR_HEADER h;
  RR_RECORD r;
  int n, i, t;
  int nb;
  int m, k;
  float *params, *strength;
  float e, eph, ee, phi, rr;

  nb = 0;

  while (1) {
    n = fread(&h, sizeof(RR_HEADER), 1, f1);
    if (n != 1) break;
    if (swp) SwapEndianRRHeader(&h);
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "QKMODE\t= %d\n", h.qk_mode);
    fprintf(f2, "MULTIP\t= %d\n", h.multipole);
    fprintf(f2, "NPARAMS\t= %d\n", h.nparams);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);
    h.tegrid = (double *) malloc(sizeof(double)*h.n_tegrid);
    n = fread(h.tegrid, sizeof(double), h.n_tegrid, f1);
    for (i = 0; i < h.n_tegrid; i++) {
      if (swp) SwapEndian((char *) &(h.tegrid[i]), sizeof(double));
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    h.egrid = (double *) malloc(sizeof(double)*h.n_egrid);
    n = fread(h.egrid, sizeof(double), h.n_egrid, f1);
    for (i = 0; i < h.n_egrid; i++) {
      if (swp) SwapEndian((char *) &(h.egrid[i]), sizeof(double));
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    h.usr_egrid = (double *) malloc(sizeof(double)*h.n_usr);
    n = fread(h.usr_egrid, sizeof(double), h.n_usr, f1);
    for (i = 0; i < h.n_usr; i++) {
      if (swp) SwapEndian((char *) &(h.usr_egrid[i]), sizeof(double));
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }
    
    if (h.qk_mode == QK_FIT) {
      m = h.nparams;
      params = (float *) malloc(sizeof(float)*m);
    }
    m = h.n_usr;
    strength = (float *) malloc(sizeof(float)*m);
    for (i = 0; i < h.ntransitions; i++) {
      n = fread(&r, sizeof(RR_RECORD), 1, f1);
      if (swp) SwapEndianRRRecord(&r);
      
      if (h.qk_mode == QK_FIT) {
	m = h.nparams;
	n = fread(params, sizeof(float), m, f1);
      }
      m = h.n_usr;
      n = fread(strength, sizeof(float), m, f1);
      if (v) {
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "%5d\t%2d\t%5d\t%2d\t%11.4E\t%2d\n",
		r.b, mem_en_table[r.b].j, 
		r.f, mem_en_table[r.f].j,
		(e*HARTREE_EV), r.kl);
	
      } else {
	fprintf(f2, "%5d\t%5d\t%2d\n", r.b, r.f, r.kl);
      }
      
      if (h.qk_mode == QK_FIT) {
	for (t = 0; t < h.nparams; t++) {
	  if (swp) SwapEndian((char *) &(params[t]), sizeof(float));
	  if (v && t == h.nparams-1) {
	    fprintf(f2, "%11.4E ", params[t]*HARTREE_EV);
	  } else {
	    fprintf(f2, "%11.4E ", params[t]);
	  }
	}
	fprintf(f2, "\n");
      }
      
      for (t = 0; t < h.n_usr; t++) {
	if (swp) SwapEndian((char *) &(strength[t]), sizeof(float));
	if (v) {
	  if (h.usr_egrid_type == 0) {
	    eph = h.usr_egrid[t];
	    ee = eph - e;
	  } else {
	    ee = h.usr_egrid[t];
	    eph = ee + e;
	  }
	  phi = 2.0*PI*FINE_STRUCTURE_CONST*strength[t]*AREA_AU20;
	  rr = phi * pow(FINE_STRUCTURE_CONST*eph, 2) / (2.0*ee);
	  phi /= (mem_en_table[r.b].j + 1.0);
	  rr /= (mem_en_table[r.f].j + 1.0);
	  fprintf(f2, "%11.4E\t%11.4E\t%11.4E\t%11.4E\n",
		  h.usr_egrid[t]*HARTREE_EV, rr, phi, strength[t]);
	} else {
	  fprintf(f2, "%11.4E\t%11.4E\n", h.usr_egrid[t], strength[t]);
	}
      }
    }

    if (h.qk_mode == QK_FIT) free(params);
    free(strength);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);

    nb++;
  }

  return nb;
}

int WriteAIRecord(FILE *f, AI_RECORD *r) {
  int n;
  ai_header.ntransitions += 1;
  ai_header.length += sizeof(AI_RECORD);
  n = fwrite(r, sizeof(AI_RECORD), 1, f);
  return n;
}

int PrintAITable(FILE *f1, FILE *f2, int v, int swp) {
  AI_HEADER h;
  AI_RECORD r;
  int n, i;
  int nb;
  float e, sdr;
  
  nb = 0;
  
  while (1) {
    n = fread(&h, sizeof(AI_HEADER), 1, f1);
    if (n != 1) break;
    if (swp) SwapEndianAIHeader(&h);
 
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "CHANNE\t= %d\n", h.channel);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    h.egrid = (double *) malloc(sizeof(double)*h.n_egrid);
    n = fread(h.egrid, sizeof(double), h.n_egrid, f1);
    for (i = 0; i < h.n_egrid; i++) {
      if (swp) SwapEndian((char *) &(h.egrid[i]), sizeof(double));
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
       
    for (i = 0; i < h.ntransitions; i++) {
      n = fread(&r, sizeof(AI_RECORD), 1, f1);
      if (swp) SwapEndianAIRecord(&r);
      if (v) {
	e = mem_en_table[r.b].energy - mem_en_table[r.f].energy;
	sdr = 0.5*(mem_en_table[r.b].j + 1.0);
	sdr *= PI*PI*r.rate/(e*(mem_en_table[r.f].j + 1.0));
	sdr *= AREA_AU20*HARTREE_EV;
	fprintf(f2, "%5d\t%2d%5d\t%2d\t%11.4E\t%11.4E\t%11.4E\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, (RATE_AU*r.rate), sdr);
      } else {
	fprintf(f2, "%5d\t%5d\t%15.8E\n", r.b, r.f, r.rate);
      }
    }
    
    free(h.egrid);
    nb++;
  }

  return nb;
}

int WriteCIRecord(FILE *f, CI_RECORD *r) {
  int n;
  int m;

  ci_header.ntransitions += 1;
  m = sizeof(CI_RECORD);
  ci_header.length += m;
  n = fwrite(r, m, 1, f);
  m = ci_header.nparams;
  ci_header.length += sizeof(float)*m;
  n = fwrite(r->params, sizeof(float), m, f);
  m = ci_header.n_usr;
  ci_header.length += sizeof(float)*m;
  n = fwrite(r->strength, sizeof(float), m, f);
  
  return n;
}

int PrintCITable(FILE *f1, FILE *f2, int v, int swp) {
  CI_HEADER h;
  CI_RECORD r;
  int n, i, t;
  int nb;
  int m, k;
  float *params, *strength;
  float e, a;

  nb = 0;

  while (1) {
    n = fread(&h, sizeof(CI_HEADER), 1, f1);
    if (n != 1) break;
    if (swp) SwapEndianCIHeader(&h);
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "QKMODE\t= %d\n", h.qk_mode);
    fprintf(f2, "NPARAMS\t= %d\n", h.nparams);
    fprintf(f2, "PWTYPE\t= %d\n", h.pw_type);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);
    h.tegrid = (double *) malloc(sizeof(double)*h.n_tegrid);
    n = fread(h.tegrid, sizeof(double), h.n_tegrid, f1);
    for (i = 0; i < h.n_tegrid; i++) {      
      if (swp) SwapEndian((char *) &(h.tegrid[i]), sizeof(double));
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    h.egrid = (double *) malloc(sizeof(double)*h.n_egrid);
    n = fread(h.egrid, sizeof(double), h.n_egrid, f1);
    for (i = 0; i < h.n_egrid; i++) {
      if (swp) SwapEndian((char *) &(h.egrid[i]), sizeof(double));
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    h.usr_egrid = (double *) malloc(sizeof(double)*h.n_usr);
    n = fread(h.usr_egrid, sizeof(double), h.n_usr, f1);
    for (i = 0; i < h.n_usr; i++) {
      if (swp) SwapEndian((char *) &(h.usr_egrid[i]), sizeof(double));
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }

    m = h.nparams;
    params = (float *) malloc(sizeof(float)*m);
    m = h.n_usr;
    strength = (float *) malloc(sizeof(float)*m);
    
    for (i = 0; i < h.ntransitions; i++) {
      n = fread(&r, sizeof(CI_RECORD), 1, f1);
      if (swp) SwapEndianCIRecord(&r);
      m = h.nparams;
      n = fread(params, sizeof(float), m, f1);
      m = h.n_usr;
      n = fread(strength, sizeof(float), m, f1);
      
      if (v) {
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "%5d\t%2d\t%5d\t%2d\t%11.4E\t%2d\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, r.kl);
      } else {
	fprintf(f2, "%5d\t%5d\t%2d\n", r.b, r.f, r.kl);
      }
      
      for (t = 0; t < h.nparams; t++) {
	if (swp) SwapEndian((char *) &(params[t]), sizeof(float));
	fprintf(f2, "%11.4E ", params[t]);
      }
      fprintf(f2, "\n");
      for (t = 0; t < h.n_usr; t++) {
	if (swp) SwapEndian((char *) &(strength[t]), sizeof(float));
	if (v) {
	  a = h.usr_egrid[t];
	  if (h.usr_egrid_type == 1) a += e;
	  a *= 1.0 + FINE_STRUCTURE_CONST2*a;
	  a = AREA_AU20/(2.0*a*(mem_en_table[r.b].j + 1.0));
	  a *= strength[t];
	  fprintf(f2, "%11.4E\t%11.4E\t%11.4E\n",
		  h.usr_egrid[t]*HARTREE_EV, strength[t], a);
	} else {
	  fprintf(f2, "%11.4E\t%11.4E\n", h.usr_egrid[t], strength[t]);
	}
      }
    }
    
    free(params); 
    free(strength);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
    
    nb++;
  }

  return nb;
}

int WriteSPRecord(FILE *f, SP_RECORD *r) {
  int n;
  sp_header.ntransitions += 1;
  sp_header.length += sizeof(SP_RECORD);
  n = fwrite(r, sizeof(SP_RECORD), 1, f);  
  return n;
}

int PrintSPTable(FILE *f1, FILE *f2, int v, int swp) {
  SP_HEADER h;
  SP_RECORD r;
  int n, i;
  int nb;
  float e, a;

  nb = 0;
  
  while (1) {
    n = fread(&h, sizeof(SP_HEADER), 1, f1);
    if (n != 1) break;
    if (swp) SwapEndianSPHeader(&h);
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "TYPE\t= %06d\n", h.type);
    fprintf(f2, "IBLK\t= %d\n", h.iblock);
    fprintf(f2, "ICOMP\t= %s\n", h.icomplex);
    fprintf(f2, "FBLK\t= %d\n", h.fblock);
    fprintf(f2, "FCOMP\t= %s\n", h.fcomplex);

    for (i = 0; i < h.ntransitions; i++) {
      n = fread(&r, sizeof(SP_RECORD), 1, f1);
      if (swp) SwapEndianSPRecord(&r);
      e = r.energy;
      if (v) e *= HARTREE_EV;
      a = r.strength;
      fprintf(f2, "%5d\t%5d\t%15.8E\t%15.8E\n", r.upper, r.lower, e, a);
    }
    nb += 1;
  }
  
  return nb;
}

int WriteRTRecord(FILE *f, RT_RECORD *r) {
  int n;
  rt_header.ntransitions += 1;
  rt_header.length += sizeof(RT_RECORD);
  n = fwrite(r, sizeof(RT_RECORD), 1, f);  
  return n;
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
    n = fread(&h, sizeof(RT_HEADER), 1, f1);
    if (n != 1) break;
    if (swp) SwapEndianRTHeader(&h);
    if (h.nele != nele) {
      if (nele != -1) {
	fprintf(f2, "\n");
	fprintf(f2, " Sum  %10.4E %10.4E %10.4E %10.4E %10.4E %10.4E %3d\n",
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
    h.p_edist = (double *) malloc(sizeof(double)*h.np_edist);
    n = fread(h.p_edist, sizeof(double), h.np_edist, f1);
    for (i = 0; i < h.np_edist; i++) {
      if (swp) SwapEndian((char *) &(h.p_edist[i]), sizeof(double));
      fprintf(f2, "\t %15.8E\n", h.p_edist[i]);
    }
    fprintf(f2, "PDEN\t= %15.8E\n", h.pden);
    fprintf(f2, "PDIST\t= %d\n", h.ipdist);
    fprintf(f2, "NPPDIS\t= %d\n", h.np_pdist);
    h.p_pdist = (double *) malloc(sizeof(double)*h.np_pdist);
    n = fread(h.p_pdist, sizeof(double), h.np_pdist, f1);
    for (i = 0; i < h.np_pdist; i++) {
      if (swp) SwapEndian((char *) &(h.p_pdist[i]), sizeof(double));
      fprintf(f2, "\t %15.8E\n", h.p_pdist[i]);
    }
    free(h.p_edist);
    free(h.p_pdist);

    fprintf(f2, "Dens\t= %15.8E\n", h.nb);
    fprintf(f2,"         NB         TR         CE");
    fprintf(f2, "         RR         AI         CI\n");
    for (i = 0; i < h.ntransitions; i++) {
      n = fread(&r, sizeof(RT_RECORD), 1, f1);
      if (swp) SwapEndianRTRecord(&r);
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
  fprintf(f2, "\n Sum  %10.4E %10.4E %10.4E %10.4E %10.4E %10.4E %3d\n\n",
	  rr, dc, re, pi, ea, ci, nele);
  return nb;
}

int WriteDRRecord(FILE *f, DR_RECORD *r) {
  int n;
  dr_header.ntransitions += 1;
  dr_header.length += sizeof(DR_RECORD);
  n = fwrite(r, sizeof(DR_RECORD), 1, f);
  return n;
}

int PrintDRTable(FILE *f1, FILE *f2, int v, int swp) {
  DR_HEADER h;
  DR_RECORD r;
  int n, i;
  int nb;
  double e;
  
  nb = 0;
  while (1) {
    n = fread(&h, sizeof(DR_HEADER), 1, f1);
    if (n != 1) break;
    if (swp) SwapEndianDRHeader(&h);
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
      n = fread(&r, sizeof(DR_RECORD), 1, f1);
      if (n != 1) break;
      if (swp) SwapEndianDRRecord(&r);
      e = r.energy;
      if (v) {
	e *= HARTREE_EV;
	fprintf(f2, "%5d\t%2d\t%5d\t%2d\t%5d\t%2d %2d\t%11.4E %11.4E %11.4E %11.4E\n",
		r.ilev, r.j, h.ilev, h.j, r.ibase, h.vn, r.vl, 
		e, r.ai, r.total_rate, r.br);
      } else {
	fprintf(f2, "%5d\t%2d\t%5d\t%2d\t%11.4E %11.4E %11.4E %11.4E\n",
		r.ilev, r.j, r.ibase, r.vl, e, r.ai, r.total_rate, r.br);
      }
    }
    nb++;
  }

  return nb;
}
