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

#include "dbase.h"
#include "parser.h"
#include "angular.h"
#include "nucleus.h"
#include "structure.h"
#include "cf77.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int version_read[NDB];
static F_HEADER fheader[NDB];
static EN_HEADER en_header;
static ENF_HEADER enf_header;
static TR_HEADER tr_header;
static TRF_HEADER trf_header;
static CE_HEADER ce_header;
static CEF_HEADER cef_header;
static CEMF_HEADER cemf_header;
static RR_HEADER rr_header;
static AI_HEADER ai_header;
static AIM_HEADER aim_header;
static CI_HEADER ci_header;
static CIM_HEADER cim_header;
static SP_HEADER sp_header;
static RT_HEADER rt_header;
static DR_HEADER dr_header;
static RO_HEADER ro_header;
static CX_HEADER cx_header;

static EN_SRECORD *mem_en_table = NULL;
static int mem_en_table_size = 0;
static EN_SRECORD *mem_enf_table = NULL;
static int mem_enf_table_size = 0;
static int iground;
static int iuta = 0;
static int utaci = 1;
static int itrf = 0;
static double clock_start=0, clock_last=0;

static double _tpbez[] = {0.0, 2.30258509, 2.89037176, 3.25809654,
			  3.68887945, 3.98898405,
			  4.36944785, 4.52178858};
static double _tpber[] = {-41.609, -20.88578477, -15.24262691,
			  -12.59173513, -10.05431044,
			  -8.29404964,  -5.62682143,  -4.50986001};

static double born_mass = 1.0;
static FORM_FACTOR bform = {0.0, -1, NULL, NULL, NULL};

static IDXMAP _idxmap = {0, 0, 0, NULL};

#define _WSF0(sv, f) do{				\
    n = FWRITE(&(sv), sizeof(sv), 1, f);		\
    m += sizeof(sv);					\
  }while(0)
#define _RSF0(sv, f) do{			       	\
    n = FREAD(&(sv), sizeof(sv), 1, f);			\
    if (n != 1) return 0;                               \
    m += sizeof(sv);					\
  }while(0)
#define _WSF1(sv, s, k, f) do{				\
    n = FWRITE(sv, s, k, f);				\
    m += (s)*(k);					\
  }while(0)
#define _RSF1(sv, s, k, f) do{				\
    n = FREAD(sv, s, k, f);				\
    if ((n) != (k)) return 0;				\
    m += (s)*(k);					\
  }while(0)
#define WSF0(sv) _WSF0(sv, f)
#define WSF1(sv, s, k) _WSF1(sv, s, k, f)
#define RSF0(sv) _RSF0(sv, f)
#define RSF1(sv, s, k) _RSF1(sv, s, k, f)

void *ReallocNew(void *p, int s) {
  void *q;

  q = malloc(s);
  memcpy(q, p, s);
  free(p);
  
  return q;
}

int CompIdxRecord(const void *r1, const void *r2) {
  IDX_RECORD *i1, *i2;
  i1 = (IDX_RECORD *) r1;
  i2 = (IDX_RECORD *) r2;
  if (i1->i0 < i2->i0) return -1;
  else if (i1->i0 > i2->i0) return 1;
  else {
    if (i1->i1 < i2->i1) return -1;
    else if (i1->i1 > i2->i1) return 1;
    return 0;
  }
}

void SetBornMass(double m) {
  if (m > 0) born_mass = m*AMU;
  else if (m < 0) born_mass = -m;
  else born_mass = 1.0;
}

double BornMass(void) {
  return born_mass;
}

void SetBornFormFactor(double te, char *fn) {
  int i, n;
  double dk, df;
  FILE *f;
  char buf[1024];

  if (bform.nk > 0) {
    free(bform.k);
    free(bform.logk);
    free(bform.fk);
    bform.k = NULL;
    bform.logk = NULL;
    bform.fk = NULL;
  }
  if (te < 0) {
    bform.nk = -1;
    bform.te = 0.0;
  } else {
    bform.te = te/HARTREE_EV;
    if (fn == NULL) {
      bform.nk = 0;
    } else {
      f = fopen(fn, "r");
      if (f == NULL) {
	printf("cannot open file %s\n", fn);
	exit(1);
      }
      i = 0;
      while (1) {
	if (NULL == fgets(buf, 1024, f)) break;
	n = sscanf(buf, "%lf %lf\n", &dk, &df);
	if (n != 2) continue;
	i++;
      }
      fseek(f, 0, SEEK_SET);
      bform.nk = i;
      bform.k = malloc(sizeof(double)*bform.nk);
      bform.logk = malloc(sizeof(double)*bform.nk);
      bform.fk = malloc(sizeof(double)*bform.nk);
      i = 0;
      while (1) {
	if (NULL == fgets(buf, 1024, f)) break;
	n = sscanf(buf, "%lf %lf\n", &dk, &df);
	if (n != 2) continue;
	bform.k[i] = dk;
	bform.fk[i] = df;
	bform.logk[i] = log(dk);
	i++;
      }
      fclose(f);
    }
  }
}

int BornFormFactorTE(double *bte) {
  if (bte) *bte = bform.te;
  return bform.nk;
}

FORM_FACTOR *BornFormFactor(void) {
  return &bform;
}

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

int SwapEndianENFHeader(ENF_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->nlevels), sizeof(int));
  SwapEndian((char *) &(h->efield), sizeof(double));
  SwapEndian((char *) &(h->bfield), sizeof(double));
  SwapEndian((char *) &(h->fangle), sizeof(double));
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

int SwapEndianENFRecord(ENF_RECORD *r) {
  SwapEndian((char *) &(r->ilev), sizeof(int));
  SwapEndian((char *) &(r->energy), sizeof(double));
  SwapEndian((char *) &(r->pbasis), sizeof(int));
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

int SwapEndianTRFHeader(TRF_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->gauge), sizeof(int));
  SwapEndian((char *) &(h->mode), sizeof(int));
  SwapEndian((char *) &(h->multipole), sizeof(int));
  SwapEndian((char *) &(h->efield), sizeof(double));
  SwapEndian((char *) &(h->bfield), sizeof(double));
  SwapEndian((char *) &(h->fangle), sizeof(double));
  return 0;
}

int SwapEndianTRFRecord(TRF_RECORD *r) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
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
 
int SwapEndianCEFHeader(CEF_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->n_tegrid), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->te0), sizeof(float));
  SwapEndian((char *) &(h->efield), sizeof(double));
  SwapEndian((char *) &(h->bfield), sizeof(double));
  SwapEndian((char *) &(h->fangle), sizeof(double));
  return 0;
}

int SwapEndianCEFRecord(CEF_RECORD *r) {
  int m;

  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->bethe), sizeof(float));
  for (m = 0; m < 2; m++) {
    SwapEndian((char *) &(r->born[m]), sizeof(float));
  }
  return 0;
}

int SwapEndianCEMFHeader(CEMF_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->n_tegrid), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->n_thetagrid), sizeof(int));
  SwapEndian((char *) &(h->n_phigrid), sizeof(int));  
  SwapEndian((char *) &(h->te0), sizeof(float));
  SwapEndian((char *) &(h->efield), sizeof(double));
  SwapEndian((char *) &(h->bfield), sizeof(double));
  SwapEndian((char *) &(h->fangle), sizeof(double));
  return 0;
}

int SwapEndianCEMFRecord(CEMF_RECORD *r) {
  int m;

  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
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
 
int SwapEndianROHeader(RO_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  return 0;
}
 
int SwapEndianCXHeader(CX_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->tgtz), sizeof(double));
  SwapEndian((char *) &(h->tgtm), sizeof(double));
  SwapEndian((char *) &(h->tgta), sizeof(double));
  SwapEndian((char *) &(h->tgtb), sizeof(double));
  SwapEndian((char *) &(h->tgte), sizeof(double));
  SwapEndian((char *) &(h->tgtx), sizeof(double));
  SwapEndian((char *) &(h->ldist), sizeof(int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ne0), sizeof(int));
  return 0;
}

int SwapEndianRORecord(RO_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->n), sizeof(int));
  return 0;
}

int SwapEndianCXRecord(CX_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->vnl), sizeof(int));
  return 0;
}

int SwapEndianAIHeader(AI_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->emin), sizeof(float));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  return 0;
}

int SwapEndianAIMHeader(AIM_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->emin), sizeof(float));
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
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->iedist), sizeof(int));
  SwapEndian((char *) &(h->np_edist), sizeof(int));
  SwapEndian((char *) &(h->eden), sizeof(float));
  SwapEndian((char *) &(h->ipdist), sizeof(int));
  SwapEndian((char *) &(h->np_pdist), sizeof(int));
  SwapEndian((char *) &(h->pden), sizeof(float));
  return 0;
}

int SwapEndianRTRecord(RT_RECORD *r) {
  SwapEndian((char *) &(r->dir), sizeof(int));
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
  return 0;
}

void CEMF2CEFHeader(CEMF_HEADER *mh, CEF_HEADER *h) {
  h->position = mh->position;
  h->length = mh->length;
  h->nele = mh->nele;
  h->ntransitions = mh->ntransitions;
  h->n_tegrid = mh->n_tegrid;
  h->n_egrid = mh->n_egrid;
  h->te0 = mh->te0;
  h->efield = mh->efield;
  h->bfield = mh->bfield;
  h->fangle = mh->fangle;
  h->tegrid = mh->tegrid;
  h->egrid = mh->egrid;
}

void CEMF2CEFRecord(CEMF_RECORD *mr, CEF_RECORD *r, CEMF_HEADER *mh, 
		    int ith, int iph) {
  int k;

  r->lower = mr->lower;
  r->upper = mr->upper;
  k = ith*mh->n_phigrid + iph;
  r->bethe = mr->bethe[k];
  r->born[0] = mr->born[k];
  r->born[1] = mr->born[mh->n_phigrid*mh->n_thetagrid];
  r->strength = &(mr->strength[k*mh->n_egrid]);
}

double ClockStart(void) {
  return clock_start;
}

double ClockLast(void) {
  return clock_last;
}

double ClockNow(int m) {
  double t = WallTime();
  double dt;
  if (m%2 == 0) {
    dt = t-clock_start;
  } else {
    dt = t-clock_last;
  }
  if (m < 2) clock_last = t;
  return dt;
}

void PrintWallTime(char *s, int m) {
  double t = ClockNow(m);
  double ttskip = 0, ttlock = 0, mtskip = 0, mtlock = 0;
  long long tnlock = 0, mnlock = 0;
#pragma omp parallel default(shared)
  {
    double tskip = TimeSkip();
    double tlock = TimeLock();
    long long nlock = NumLock();
#pragma omp critical
    {
      ttskip += tskip;
      ttlock += tlock;
      tnlock += nlock;
      if (tskip > mtskip) mtskip = tskip;
      if (tlock > mtlock) mtlock = tlock;
      if (nlock > mnlock) mnlock = nlock;
    }
  }
  printf("WallTime%d: %9.3E %9.3E %9.3E/%9.3E %9.3E/%9.3E %8lld %8lld ... %s\n",
	 m, t, TotalSize(), ttskip, mtskip, ttlock, mtlock, tnlock, mnlock, s);
  fflush(stdout);
}

int InitDBase(void) {
  int i;

  clock_start = WallTime();
  clock_last = clock_start;
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
  mem_enf_table = NULL;
  mem_enf_table_size = 0;
  iground = 0;
  itrf = 0;

  if (_idxmap.imap == NULL) {
    _idxmap.imap = malloc(sizeof(IDXMAP));
    ArrayInit(_idxmap.imap, sizeof(IDXDAT), 5000);
  }
  ClearIdxMap();
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
  if (mem_enf_table) {
    free(mem_enf_table);
    mem_enf_table = NULL;
    mem_enf_table_size = 0;
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
    ClearIdxMap();
    return 0;
  }
}


int ReadENHeaderOld(TFILE *f, EN_HEADER *h, int swp) {
  int n;

  n = FREAD(h, sizeof(EN_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianENHeader(h);
  
  return sizeof(EN_HEADER);
}
  
int ReadENRecordOld(TFILE *f, EN_RECORD *r, int swp) {
  int n;

  n = FREAD(r, sizeof(EN_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianENRecord(r);
  
  return sizeof(EN_RECORD);
}

int ReadTRHeaderOld(TFILE *f, TR_HEADER *h, int swp) {
  int n;

  n = FREAD(h, sizeof(TR_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianTRHeader(h);
  if (h->length/h->ntransitions > sizeof(TR_RECORD)) iuta = 1;
  else iuta =0;
  return sizeof(TR_HEADER);
}

int ReadTRRecordOld(TFILE *f, TR_RECORD *r, TR_EXTRA *rx, int swp) {
  int n;
  
  n = FREAD(r, sizeof(TR_RECORD), 1, f);
  if (n != 1) return 0;
  if (iuta) {
    n = FREAD(rx, sizeof(TR_EXTRA), 1, f);
    if (n != 1) return 0;
  }
  if (swp) SwapEndianTRRecord(r, rx);
  
  return sizeof(TR_RECORD);
}

int ReadCEHeaderOld(TFILE *f, CE_HEADER *h, int swp) {
  int i, n, m;

  n = FREAD(h, sizeof(CE_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianCEHeader(h);
  m = sizeof(CE_HEADER);

  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  n = FREAD(h->tegrid, sizeof(double), h->n_tegrid, f);
  if (n != h->n_tegrid) {
    free(h->tegrid);
    return 0;
  }
  m += sizeof(double)*h->n_tegrid;
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = FREAD(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) {
    free(h->tegrid);
    free(h->egrid);
    return 0;
  }
  m += sizeof(double)*h->n_egrid;
  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  n = FREAD(h->usr_egrid, sizeof(double), h->n_usr, f);
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

int ReadCERecordOld(TFILE *f, CE_RECORD *r, int swp, CE_HEADER *h) {
  int i, n, m, m0;

  n = FREAD(r, sizeof(CE_RECORD), 1, f);
  if (n != 1) return 0;
  m0 = sizeof(CE_RECORD);
  if (swp) SwapEndianCERecord(r);
  
  if (h->msub) {
    m = r->nsub;
    r->params = (float *) malloc(sizeof(float)*m);
    n = FREAD(r->params, sizeof(float), m, f);
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
    n = FREAD(r->params, sizeof(float), m, f);
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
  n = FREAD(r->strength, sizeof(float), m, f);
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

int ReadRRHeaderOld(TFILE *f, RR_HEADER *h, int swp) {
  int i, n, m;

  n = FREAD(h, sizeof(RR_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianRRHeader(h);
  m = sizeof(RR_HEADER);

  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  n = FREAD(h->tegrid, sizeof(double), h->n_tegrid, f);
  if (n != h->n_tegrid) {
    free(h->tegrid);
    return 0;
  }
  m += sizeof(double)*h->n_tegrid;
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = FREAD(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) {
    free(h->tegrid);
    free(h->egrid);
    return 0;
  }
  m += sizeof(double)*h->n_egrid;
  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  n = FREAD(h->usr_egrid, sizeof(double), h->n_usr, f);
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

int ReadRRRecordOld(TFILE *f, RR_RECORD *r, int swp, RR_HEADER *h) {
  int i, n, m, m0;
  
  n = FREAD(r, sizeof(RR_RECORD), 1, f);
  if (n != 1) return 0;
  m0 = sizeof(RR_RECORD);
  if (swp) SwapEndianRRRecord(r);
  
  if (h->qk_mode == QK_FIT) {
    m = h->nparams;
    r->params = (float *) malloc(sizeof(float)*m);
    n = FREAD(r->params, sizeof(float), m, f);
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
  n = FREAD(r->strength, sizeof(float), m, f);
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

int ReadAIHeaderOld(TFILE *f, AI_HEADER *h, int swp) {
  int i, n, m;

  n = FREAD(h, sizeof(AI_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianAIHeader(h);
  m = sizeof(AI_HEADER);
  
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = FREAD(h->egrid, sizeof(double), h->n_egrid, f);
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

int ReadAIMHeaderOld(TFILE *f, AIM_HEADER *h, int swp) {
  int i, n, m;

  n = FREAD(h, sizeof(AIM_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianAIMHeader(h);
  m = sizeof(AIM_HEADER);
  
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = FREAD(h->egrid, sizeof(double), h->n_egrid, f);
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

int ReadAIRecordOld(TFILE *f, AI_RECORD *r, int swp) {
  int n;

  n = FREAD(r, sizeof(AI_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianAIRecord(r);
  
  return sizeof(AI_RECORD);
}

int ReadAIMRecordOld(TFILE *f, AIM_RECORD *r, int swp) {
  int n, i;

  n = FREAD(r, sizeof(AIM_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) {
    SwapEndianAIMRecord(r);
  }
  r->rate = (float *) malloc(sizeof(float)*r->nsub);
  n = FREAD(r->rate, sizeof(float), r->nsub, f);
  if (n != r->nsub) return 0;
  if (swp) {
    for (i = 0; i < r->nsub; i++) {
      SwapEndian((char *) &(r->rate[i]), sizeof(float));
    }
  }
  return sizeof(AIM_RECORD);
}

int ReadCIHeaderOld(TFILE *f, CI_HEADER *h, int swp) {
  int i, n, m;

  n = FREAD(h, sizeof(CI_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianCIHeader(h);
  m = sizeof(CI_HEADER);
  
  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  n = FREAD(h->tegrid, sizeof(double), h->n_tegrid, f);
  if (n != h->n_tegrid) {
    free(h->tegrid);
    return 0;
  }
  m += sizeof(double)*h->n_tegrid;
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = FREAD(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) {
    free(h->tegrid);
    free(h->egrid);
    return 0;
  }
  m += sizeof(double)*h->n_egrid;
  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  n = FREAD(h->usr_egrid, sizeof(double), h->n_usr, f);
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

int ReadCIRecordOld(TFILE *f, CI_RECORD *r, int swp, CI_HEADER *h) {
  int i, n, m, m0;

  n = FREAD(r, sizeof(CI_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianCIRecord(r);
  m0 = sizeof(CI_RECORD);

  m = h->nparams;
  r->params = (float *) malloc(sizeof(float)*m);
  n = FREAD(r->params, sizeof(float), m, f);
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
  n = FREAD(r->strength, sizeof(float), m, f);
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

int ReadCIMHeaderOld(TFILE *f, CIM_HEADER *h, int swp) {
  int i, n, m;

  n = FREAD(h, sizeof(CIM_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianCIMHeader(h);
  m = sizeof(CIM_HEADER);
  
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  n = FREAD(h->egrid, sizeof(double), h->n_egrid, f);
  if (n != h->n_egrid) {
    free(h->egrid);
    return 0;
  }
  m += sizeof(double)*h->n_egrid;
  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  n = FREAD(h->usr_egrid, sizeof(double), h->n_usr, f);
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

int ReadCIMRecordOld(TFILE *f, CIM_RECORD *r, int swp, CIM_HEADER *h) {
  int i, n, m, m0;

  n = FREAD(r, sizeof(CIM_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianCIMRecord(r);
  m0 = sizeof(CIM_RECORD);

  m = h->n_usr*r->nsub;
  r->strength = (float *) malloc(sizeof(float)*m);
  n = FREAD(r->strength, sizeof(float), m, f);
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

int ReadSPHeaderOld(TFILE *f, SP_HEADER *h, int swp) {
  int n;

  n = FREAD(h, sizeof(SP_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianSPHeader(h);
  if (h->length/h->ntransitions > sizeof(SP_RECORD)) iuta = 1;
  else iuta = 0;
  return sizeof(SP_HEADER);
}

int ReadSPRecordOld(TFILE *f, SP_RECORD *r, SP_EXTRA *rx, int swp) {
  int n;

  n = FREAD(r, sizeof(SP_RECORD)-2*sizeof(float), 1, f);
  if (n != 1) return 0;
  if (iuta) {
    n = FREAD(rx, sizeof(SP_EXTRA), 1, f);
    if (n != 1) return 0;
  }
  if (swp) SwapEndianSPRecord(r, rx);
  return sizeof(SP_RECORD);
}

int ReadRTHeaderOld(TFILE *f, RT_HEADER *h, int swp) {
  int i, n, m;

  n = FREAD(h, sizeof(RT_HEADER)-8, 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianRTHeader(h);
  m = sizeof(RT_HEADER);
  
  h->p_edist = (double *) malloc(sizeof(double)*h->np_edist);
  n = FREAD(h->p_edist, sizeof(double), h->np_edist, f);
  if (n != h->np_edist) {
    free(h->p_edist);
    return 0;
  }
  m += sizeof(double)*h->np_edist;

  h->p_pdist = (double *) malloc(sizeof(double)*h->np_pdist);
  n = FREAD(h->p_pdist, sizeof(double), h->np_pdist, f);
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

int ReadRTRecordOld(TFILE *f, RT_RECORD *r, int swp) {
  int n;

  n = FREAD(r, sizeof(RT_RECORD)-4, 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianRTRecord(r);
  return sizeof(RT_RECORD);
}

int ReadDRHeaderOld(TFILE *f, DR_HEADER *h, int swp) {
  int n;
  
  n = FREAD(h, sizeof(DR_HEADER), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianDRHeader(h);
  
  return sizeof(DR_HEADER);
}

int ReadDRRecordOld(TFILE *f, DR_RECORD *r, int swp) {
  int n;

  n = FREAD(r, sizeof(DR_RECORD), 1, f);
  if (n != 1) return 0;
  if (swp) SwapEndianDRRecord(r);
  
  return sizeof(DR_RECORD);
}   

int WriteFHeader(TFILE *f, F_HEADER *fh) {
  int n, m = 0;
  int v = fh->version;

  if (fh->nthreads > 1) {
    v |= fh->nthreads<<16;
  }

  FFLUSH(f);
  WSF0(fh->tsession);
  WSF0(v);
  WSF0(fh->sversion);
  WSF0(fh->ssversion);
  WSF0(fh->type);
  WSF0(fh->atom);
  WSF1(fh->symbol, sizeof(char), 4);
  WSF0(fh->nblocks);
  FFLUSH(f);
  return m;
}

int ReadFHeader(TFILE *f, F_HEADER *fh, int *swp) {
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

  fh->nthreads = (fh->version&0xFFFF0000)>>16;
  fh->version &= 0xFFFF;
  SetVersionRead(fh->type, fh->version*100+fh->sversion*10+fh->ssversion);

  if (fh->type == DB_EN && version_read[DB_EN-1] == 114) {
    if (IsNewV114(f)) {
      SetVersionRead(DB_EN, 115);
    }
    FSEEK(f, m, SEEK_SET);    
  }
  if (fh->type == DB_TR && itrf >= 0) {
    if (VersionLE(fh, 1, 0, 6)) itrf = 1;
    else itrf = 0;
  }

  return m;
}

int WriteENHeader(TFILE *f, EN_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->nlevels);

  return m;
}

int WriteENFHeader(TFILE *f, ENF_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->nlevels);
  WSF0(h->efield);
  WSF0(h->bfield);
  WSF0(h->fangle);

  return m;
}

int WriteTRHeader(TFILE *f, TR_HEADER *h) {
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

int WriteTRFHeader(TFILE *f, TRF_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->gauge);
  WSF0(h->mode);
  WSF0(h->multipole);
  WSF0(h->efield);
  WSF0(h->bfield);
  WSF0(h->fangle);
  
  return m;
}

int WriteCEHeader(TFILE *f, CE_HEADER *h) {
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

int WriteCEFHeader(TFILE *f, CEF_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->n_tegrid);
  WSF0(h->n_egrid);
  WSF0(h->te0);
  WSF0(h->efield);
  WSF0(h->bfield);
  WSF0(h->fangle);
  WSF1(h->tegrid, sizeof(double), h->n_tegrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  
  return m;
}

int WriteCEMFHeader(TFILE *f, CEMF_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->n_tegrid);
  WSF0(h->n_egrid);
  WSF0(h->n_thetagrid);
  WSF0(h->n_phigrid);
  WSF0(h->te0);
  WSF0(h->efield);
  WSF0(h->bfield);
  WSF0(h->fangle);
  WSF1(h->tegrid, sizeof(double), h->n_tegrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  WSF1(h->thetagrid, sizeof(double), h->n_thetagrid);
  WSF1(h->phigrid, sizeof(double), h->n_phigrid);
  
  return m;
}

int WriteROHeader(TFILE *f, RO_HEADER *h) {
  int n, m = 0;
  
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);

  return m;
}

int WriteCXHeader(TFILE *f, CX_HEADER *h) {
  int n, m = 0;
  
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF1(h->tgts, sizeof(char), 128);
  WSF0(h->tgtz);
  WSF0(h->tgtm);
  WSF0(h->tgta);
  WSF0(h->tgtb);
  WSF0(h->tgte);
  WSF0(h->tgtx);
  WSF0(h->ldist);
  WSF0(h->te0);
  WSF0(h->ne0);
  WSF1(h->e0, sizeof(double), h->ne0);
  
  return m;
}
  
int WriteRRHeader(TFILE *f, RR_HEADER *h) {
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

int WriteAIHeader(TFILE *f, AI_HEADER *h) {
  int n, m = 0;
 
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->emin);
  WSF0(h->n_egrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  return m;
}

int WriteAIMHeader(TFILE *f, AIM_HEADER *h) {
  int n, m = 0;
    
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->emin);
  WSF0(h->n_egrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  
  return m;
}

int WriteCIHeader(TFILE *f, CI_HEADER *h) {
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

int WriteCIMHeader(TFILE *f, CIM_HEADER *h) {
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

int WriteSPHeader(TFILE *f, SP_HEADER *h) {
  int i, n, m = 0;
     
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->iblock);
  WSF0(h->fblock);
  for (i = 0; i < LNCOMPLEX; i++) {
    if (h->icomplex[i] == '\0') break;
  } 
  for (i++; i < LNCOMPLEX; i++) h->icomplex[i] = '\0';
  for (i = 0; i < LNCOMPLEX; i++) {
    if (h->fcomplex[i] == '\0') break;
  } 
  for (i++; i < LNCOMPLEX; i++) h->fcomplex[i] = '\0';      
  WSF1(h->icomplex, sizeof(char), LNCOMPLEX);
  WSF1(h->fcomplex, sizeof(char), LNCOMPLEX);
  WSF0(h->type);

  return m;
}

int WriteRTHeader(TFILE *f, RT_HEADER *h) {
  int n, m = 0;
         
  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->ntransitions);
  WSF0(h->iedist);
  WSF0(h->np_edist);
  WSF0(h->eden);
  WSF0(h->ipdist);
  WSF0(h->np_pdist);
  WSF0(h->pden);
  WSF1(h->p_edist, sizeof(double), h->np_edist);
  WSF1(h->p_pdist, sizeof(double), h->np_pdist);
  
  return m;
}

int WriteDRHeader(TFILE *f, DR_HEADER *h) {
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

int WriteENRecord(TFILE *f, EN_RECORD *r) {
  int i, n, m = 0;

  if (en_header.nlevels == 0) {
    SetLockMPI();
    if (en_header.nlevels == 0) {      
      fheader[DB_EN-1].nblocks++;
      FFLUSH(f);
      n = WriteENHeader(f, &en_header);
      FFLUSH(f);
    }
#pragma omp atomic
    en_header.nlevels += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    en_header.nlevels += 1;
  }
  
#ifdef USEBF
  BFileCheckBuf(f,		
		sizeof(r->p)+
		sizeof(r->j)+
		sizeof(r->ilev)+
		sizeof(r->energy)+
		sizeof(char)*(LNCOMPLEX+LSNAME+LNAME));
#endif
  WSF0(r->p);
  WSF0(r->j);
  WSF0(r->ilev);
  WSF0(r->ibase);
  WSF0(r->energy);
  /* make sure the strings zeroes out after NULL */
  for (i = 0; i < LNCOMPLEX; i++) {
    if (r->ncomplex[i] == '\0') break;
  }
  for (i++; i < LNCOMPLEX; i++) r->ncomplex[i] = '\0';

  for (i = 0; i < LSNAME; i++) {
    if (r->sname[i] == '\0') break;
  }
  for (i++; i < LSNAME; i++) r->sname[i] = '\0';

  for (i = 0; i < LNAME; i++) {
    if (r->name[i] == '\0') break;
  }
  for (i++; i < LNAME; i++) r->name[i] = '\0';
  
  WSF1(r->ncomplex, sizeof(char), LNCOMPLEX);
  WSF1(r->sname, sizeof(char), LSNAME);
  WSF1(r->name, sizeof(char), LNAME);

#pragma omp atomic
  en_header.length += m;

  return m;
}

int WriteENFRecord(TFILE *f, ENF_RECORD *r) {
  int n, m = 0;

  if (enf_header.nlevels == 0) {
    SetLockMPI();
    if (enf_header.nlevels == 0) {
      fheader[DB_ENF-1].nblocks++;
      FFLUSH(f);
      n = WriteENFHeader(f, &enf_header);
      FFLUSH(f);
    }
#pragma omp atomic
    enf_header.nlevels += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    enf_header.nlevels += 1;
  }
  
#ifdef USEBF
  BFileCheckBuf(f,		
		sizeof(r->ilev)+
		sizeof(r->energy)+
		sizeof(r->pbasis));
#endif
  WSF0(r->ilev);
  WSF0(r->energy);
  WSF0(r->pbasis);
  
#pragma omp atomic
  enf_header.length += m;

  return m;
}

int WriteTRRecord(TFILE *f, TR_RECORD *r, TR_EXTRA *rx) {
  int n, m = 0;

  if (tr_header.ntransitions == 0) {
    SetLockMPI();
    if (tr_header.ntransitions == 0) {
      fheader[DB_TR-1].nblocks++;
      FFLUSH(f);
      n = WriteTRHeader(f, &tr_header);
      FFLUSH(f);
    }
#pragma omp atomic
    tr_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    tr_header.ntransitions += 1;
  }
  
#ifdef USEBF
  if (!iuta) {
    BFileCheckBuf(f,		
		  sizeof(r->lower)+
		  sizeof(r->upper)+
		  sizeof(r->strength));
  } else {
    BFileCheckBuf(f,		
		  sizeof(r->lower)+
		  sizeof(r->upper)+
		  sizeof(r->strength)+
		  sizeof(rx->energy)+
		  sizeof(rx->sdev)+
		  sizeof(rx->sci));
  }
#endif
  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->strength);

  if (iuta) {
    WSF0(rx->energy);
    WSF0(rx->sdev);
    WSF0(rx->sci);
  }

#pragma omp atomic
  tr_header.length += m;

  return m;
}

int WriteTRFRecord(TFILE *f, TRF_RECORD *r) {
  int n, m = 0;

  if (trf_header.ntransitions == 0) {
    SetLockMPI();
    if (trf_header.ntransitions == 0) {
      fheader[DB_TRF-1].nblocks++;
      FFLUSH(f);
      n = WriteTRFHeader(f, &trf_header);
      FFLUSH(f);
    }
#pragma omp atomic
    trf_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    trf_header.ntransitions += 1;
  }
  
#ifdef USEBF
  BFileCheckBuf(f,		
		sizeof(r->lower)+
		sizeof(r->upper)+
		sizeof(float)*(2*abs(trf_header.multipole)+1));
#endif
  WSF0(r->lower);
  WSF0(r->upper);
  WSF1(r->strength, sizeof(float), 2*abs(trf_header.multipole)+1);

#pragma omp atomic
  trf_header.length += m;

  return m;
}

int WriteCERecord(TFILE *f, CE_RECORD *r) {
  int n;
  int m0, m = 0;
    
  if (ce_header.ntransitions == 0) {
    SetLockMPI();
    if (ce_header.ntransitions == 0) {
      fheader[DB_CE-1].nblocks++;
      FFLUSH(f);
      n = WriteCEHeader(f, &ce_header);
      FFLUSH(f);
    }
#pragma omp atomic
    ce_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    ce_header.ntransitions += 1;
  }
  
  if (ce_header.msub) {
    m0 = r->nsub;
  } else if (ce_header.qk_mode == QK_FIT) {
    m0 = ce_header.nparams * r->nsub;
  } else m0 = 0;
  
  int m1 = ce_header.n_usr*r->nsub;
  
#ifdef USEBF
  BFileCheckBuf(f,		
		sizeof(r->lower)+
		sizeof(r->upper)+
		sizeof(r->nsub)+
		sizeof(r->bethe)+
		sizeof(float)*(2+m0+m1));  
#endif
  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->nsub);
  WSF0(r->bethe);
  WSF1(r->born, sizeof(float), 2);

  if (m0) {
    WSF1(r->params, sizeof(float), m0);
  }
  WSF1(r->strength, sizeof(float), m1);

#pragma omp atomic
  ce_header.length += m;

  return m;
}

int WriteCEFRecord(TFILE *f, CEF_RECORD *r) {
  int n;
  int m0, m = 0;

  if (cef_header.ntransitions == 0) {
    SetLockMPI();
    if (cef_header.ntransitions == 0) {
      fheader[DB_CEF-1].nblocks++;
      FFLUSH(f);
      n = WriteCEFHeader(f, &cef_header);
      FFLUSH(f);
    }
#pragma omp atomic
    cef_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    cef_header.ntransitions += 1;
  }
  
  m0 = cef_header.n_egrid;
#ifdef USEBF
  BFileCheckBuf(f,		
		sizeof(r->lower)+
		sizeof(r->upper)+
		sizeof(r->bethe)+
		sizeof(float)*(2+m0));  
#endif
  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->bethe);
  WSF1(r->born, sizeof(float), 2);
  WSF1(r->strength, sizeof(float), m0);
  
#pragma omp atomic
  cef_header.length += m;

  return m;
}

int WriteCEMFRecord(TFILE *f, CEMF_RECORD *r) {
  int n;
  int m0, m = 0;

  if (cemf_header.ntransitions == 0) {
    SetLockMPI();
    if (cemf_header.ntransitions == 0) {
      fheader[DB_CEMF-1].nblocks++;
      FFLUSH(f);
      n = WriteCEMFHeader(f, &cemf_header);
      FFLUSH(f);
    }
#pragma omp atomic
    cemf_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    cemf_header.ntransitions += 1;
  }
  
  m0 = cemf_header.n_thetagrid * cemf_header.n_phigrid;
  int m1 = cemf_header.n_egrid * m0;
#ifdef USEBF
  BFileCheckBuf(f,		
		sizeof(r->lower)+
		sizeof(r->upper)+
		sizeof(float)*(2*m0+1+m1));  
#endif
  WSF0(r->lower);
  WSF0(r->upper);
  
  WSF1(r->bethe, sizeof(float), m0);
  WSF1(r->born, sizeof(float), m0+1);

  WSF1(r->strength, sizeof(float), m1);
  
#pragma omp atomic
  cemf_header.length += m;

  return m;
}

int WriteRORecord(TFILE *f, RO_RECORD *r) {
  int n;
  int m = 0;

  if (ro_header.ntransitions == 0) {
    SetLockMPI();
    if (ro_header.ntransitions == 0) {
      fheader[DB_RO-1].nblocks++;
      FFLUSH(f);
      n = WriteROHeader(f, &ro_header);
      FFLUSH(f);
    }
#pragma omp atomic
    ro_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    ro_header.ntransitions += 1;
  }
  
#ifdef USEBF
  BFileCheckBuf(f,		
		sizeof(r->b)+
		sizeof(r->f)+
		sizeof(r->n)+
		sizeof(int)*r->n+
		sizeof(double)*2*r->n);
#endif
  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->n);
  WSF1(r->nk, sizeof(int), r->n);
  WSF1(r->nq, sizeof(double), r->n);
  WSF1(r->dn, sizeof(double), r->n);
#pragma omp atomic
  ro_header.length += m;

  return m;
}

int WriteCXRecord(TFILE *f, CX_RECORD *r) {
  int n;
  int m = 0;

  if (cx_header.ntransitions == 0) {
    SetLockMPI();
    if (cx_header.ntransitions == 0) {
      fheader[DB_CX-1].nblocks++;
      FFLUSH(f);
      n = WriteCXHeader(f, &cx_header);
      FFLUSH(f);
    }
#pragma omp atomic
    cx_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    cx_header.ntransitions += 1;
  }
  
#ifdef USEBF
  BFileCheckBuf(f,		
		sizeof(r->b)+
		sizeof(r->f)+
		sizeof(r->vnl)+
		sizeof(double)*cx_header.ne0);
#endif
  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->vnl);
  WSF1(r->cx, sizeof(double), cx_header.ne0);
#pragma omp atomic
  cx_header.length += m;

  return m;
}

int WriteRRRecord(TFILE *f, RR_RECORD *r) {
  int n;
  int m = 0, m0;

  if (rr_header.ntransitions == 0) {
    SetLockMPI();
    if (rr_header.ntransitions == 0) {
      fheader[DB_RR-1].nblocks++;
      FFLUSH(f);
      n = WriteRRHeader(f, &rr_header);
      FFLUSH(f);
    }
#pragma omp atomic
    rr_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    rr_header.ntransitions += 1;
  }
  
#ifdef USEBF
  if (rr_header.qk_mode == QK_FIT) {
    BFileCheckBuf(f,		
		  sizeof(r->b)+
		  sizeof(r->f)+
		  sizeof(r->kl)+
		  sizeof(float)*(rr_header.nparams+rr_header.n_usr));
  } else {
    BFileCheckBuf(f,		
		  sizeof(r->b)+
		  sizeof(r->f)+
		  sizeof(r->kl)+
		  sizeof(float)*rr_header.n_usr);
  }
#endif
  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->kl);

  if (rr_header.qk_mode == QK_FIT) {
    m0 = rr_header.nparams;
    WSF1(r->params, sizeof(float), m0);
  }
  m0 = rr_header.n_usr;
  WSF1(r->strength, sizeof(float), m0);

#pragma omp atomic
  rr_header.length += m;

  return m;
}

int WriteAIRecord(TFILE *f, AI_RECORD *r) {
  int n, m = 0;

  if (ai_header.ntransitions == 0) {
    SetLockMPI();
    if (ai_header.ntransitions == 0) {
      fheader[DB_AI-1].nblocks++;
      FFLUSH(f);
      WriteAIHeader(f, &ai_header);
      FFLUSH(f);
    }
#pragma omp atomic
    ai_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    ai_header.ntransitions += 1;
  }
#ifdef USEBF
  BFileCheckBuf(f,		
		sizeof(r->b)+
		sizeof(r->f)+
		sizeof(r->rate));
#endif
  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->rate);

#pragma omp atomic
  ai_header.length += m;
  return m;
}

int WriteAIMRecord(TFILE *f, AIM_RECORD *r) {
  int n, m = 0;

  if (aim_header.ntransitions == 0) {
    SetLockMPI();
    if (aim_header.ntransitions == 0) {
      fheader[DB_AIM-1].nblocks++;
      FFLUSH(f);
      WriteAIMHeader(f, &aim_header);
      FFLUSH(f);
    }
#pragma omp atomic
    aim_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    aim_header.ntransitions += 1;
  }
#ifdef USEBF
  BFileCheckBuf(f,
	       sizeof(r->b)+
	       sizeof(r->f)+
	       sizeof(r->nsub)+
	       sizeof(float)*r->nsub);
#endif
  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->nsub);
  WSF1(r->rate, sizeof(float), r->nsub);
  
#pragma omp atomic
  aim_header.length += m;

  return m;
}

int WriteCIRecord(TFILE *f, CI_RECORD *r) {
  int n;
  int m = 0, m0;

  if (ci_header.ntransitions == 0) {
    SetLockMPI();
    if (ci_header.ntransitions == 0) {
      fheader[DB_CI-1].nblocks++;
      FFLUSH(f);
      WriteCIHeader(f, &ci_header);
      FFLUSH(f);
    }
#pragma omp atomic
    ci_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    ci_header.ntransitions += 1;
  }
#ifdef USEBF
  BFileCheckBuf(f,
		sizeof(r->b)+
		sizeof(r->f)+
		sizeof(r->kl)+
		sizeof(float)*(ci_header.nparams+ci_header.n_usr));
#endif
  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->kl);
  m0 = ci_header.nparams;
  WSF1(r->params, sizeof(float), m0);
  m0 = ci_header.n_usr;
  WSF1(r->strength, sizeof(float), m0);

#pragma omp atomic
  ci_header.length += m;

  return m;
}

int WriteCIMRecord(TFILE *f, CIM_RECORD *r) {
  int n;
  int m = 0, m0;

  if (cim_header.ntransitions == 0) {
    SetLockMPI();
    if (cim_header.ntransitions == 0) {
      fheader[DB_CIM-1].nblocks++;
      FFLUSH(f);
      WriteCIMHeader(f, &cim_header);
      FFLUSH(f);
    }
#pragma omp atomic
    cim_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    cim_header.ntransitions += 1;
  }
  
  m0 = r->nsub*cim_header.n_usr;
#ifdef USEBF
  BFileCheckBuf(f,
		sizeof(r->b)+
		sizeof(r->f)+
		sizeof(r->nsub)+
		sizeof(float)*m0);
#endif
  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->nsub);
  WSF1(r->strength, sizeof(float), m0);
  
#pragma omp atomic
  cim_header.length += m;

  return m;
}

int WriteSPRecord(TFILE *f, SP_RECORD *r, SP_EXTRA *rx) {
  int n, m = 0;

  if (sp_header.ntransitions == 0) {
    SetLockMPI();
    if (sp_header.ntransitions == 0) {
      fheader[DB_SP-1].nblocks++;
      FFLUSH(f);
      WriteSPHeader(f, &sp_header);
      FFLUSH(f);
    }
#pragma omp atomic
    sp_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    sp_header.ntransitions += 1;
  }
#ifdef USEBF
  if (iuta) {
    BFileCheckBuf(f,
		  sizeof(r->lower)+
		  sizeof(r->upper)+
		  sizeof(r->energy)+
		  sizeof(r->strength)+
		  sizeof(r->rrate)+
		  sizeof(r->trate)+
		  sizeof(rx->sdev));
  } else {
    BFileCheckBuf(f,
		  sizeof(r->lower)+
		  sizeof(r->upper)+
		  sizeof(r->energy)+
		  sizeof(r->strength)+
		  sizeof(r->rrate)+
		  sizeof(r->trate));
  }
#endif
  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->energy);
  WSF0(r->strength);
  WSF0(r->rrate);
  WSF0(r->trate);
  if (iuta) {
    WSF0(rx->sdev);
  }

#pragma omp atomic
  sp_header.length += m;

  return m;
}

int WriteRTRecord(TFILE *f, RT_RECORD *r) {
  int i, n, m = 0;

  if (rt_header.ntransitions == 0) {
    SetLockMPI();
    if (rt_header.ntransitions == 0) {
      fheader[DB_RT-1].nblocks++;
      FFLUSH(f);
      WriteRTHeader(f, &rt_header);
      FFLUSH(f);
    }
#pragma omp atomic
    rt_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
  rt_header.ntransitions += 1;
  }
#ifdef USEBF
  BFileCheckBuf(f,
		sizeof(r->dir)+
		sizeof(r->iblock)+
		sizeof(r->nb)+
		sizeof(r->tr)+
		sizeof(r->ce)+
		sizeof(r->rr)+
		sizeof(r->ai)+
		sizeof(r->ci)+
		sizeof(char)*LNCOMPLEX);
#endif
  WSF0(r->dir);
  WSF0(r->iblock);
  WSF0(r->nb);
  WSF0(r->tr);
  WSF0(r->ce);
  WSF0(r->rr);
  WSF0(r->ai);
  WSF0(r->ci);
  for (i = 0; i < LNCOMPLEX; i++) {
    if (r->icomplex[i] == '\0') break;
  } 
  for (i++; i < LNCOMPLEX; i++) r->icomplex[i] = '\0';
  WSF1(r->icomplex, sizeof(char), LNCOMPLEX);

#pragma omp atomic
  rt_header.length += m;

  return m;
}

int WriteDRRecord(TFILE *f, DR_RECORD *r) {
  int n, m = 0;

  if (dr_header.ntransitions == 0) {
    SetLockMPI();
    if (dr_header.ntransitions == 0) {
      fheader[DB_DR-1].nblocks++;
      FFLUSH(f);
      WriteDRHeader(f, &dr_header);
      FFLUSH(f);
    }
#pragma omp atomic
    dr_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    dr_header.ntransitions += 1;
  }
#ifdef USEBF
  BFileCheckBuf(f,
		sizeof(r->ilev)+
		sizeof(r->flev)+
		sizeof(r->ibase)+
		sizeof(r->fbase)+
		sizeof(r->vl)+
		sizeof(r->j)+
		sizeof(r->energy)+
		sizeof(r->etrans)+
		sizeof(r->br)+
		sizeof(r->ai)+
		sizeof(r->total_rate));
#endif
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

#pragma omp atomic
  dr_header.length += m;

  return m;
}

int ReadENHeader(TFILE *f, EN_HEADER *h, int swp) {
  int n, m = 0;
  
  if (version_read[DB_EN-1] < 109) return ReadENHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->nlevels);

  if (swp) SwapEndianENHeader(h);
  
  return m;
}
 
int ReadENFHeader(TFILE *f, ENF_HEADER *h, int swp) {
  int n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->nlevels);
  RSF0(h->efield);
  RSF0(h->bfield);
  RSF0(h->fangle);
  
  if (swp) SwapEndianENFHeader(h);

  return m;
}

int IsNewV114(TFILE *f) {
  EN_HEADER h;
  EN_RECORD r;
  int i, swp, n, m = 0;

  while (1) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.nlevels; i++) {
      n = ReadENRecord(f, &r, swp);
      if (n == 0) break;
      if (!isdigit(r.name[0])) {
	return 1;
      }
      if (!isalpha(r.name[1]) && r.name[1] != '[') {
	return 1;
      }
      if (r.name[1] == '[') {
	char *c = &r.name[2];
	while (*c && *c != ']') c++;
	if (c[0] == ']' && (c[1] != '-' && c[1] != '+')) {
	  return 1;
	}
	return 0;
      }
      if (r.name[2] != '-' && r.name[2] != '+') {
	return 1;
      }
      return 0;
    }
  }
  return 0;
}

int ReadENRecord(TFILE *f, EN_RECORD *r, int swp) {
  int n, m = 0;

  if (version_read[DB_EN-1] < 109) return ReadENRecordOld(f, r, swp);

  RSF0(r->p);
  RSF0(r->j);
  RSF0(r->ilev);
  RSF0(r->ibase);
  RSF0(r->energy);
  RSF1(r->ncomplex, sizeof(char), LNCOMPLEX);
  if (version_read[DB_EN-1] < 115) {
    RSF1(r->sname, sizeof(char), LSNAME0);
    RSF1(r->name, sizeof(char), LNAME0);
    StrReplace(LNCOMPLEX, r->ncomplex, ' ', '.', '.', '\0');
    StrReplace(LSNAME0, r->sname, ' ', '.', '.', '\0');
    StrReplace(LNAME0, r->name, ' ', '.', '.', '\0');
  } else {
    RSF1(r->sname, sizeof(char), LSNAME);
    RSF1(r->name, sizeof(char), LNAME);    
  }
  if (swp) SwapEndianENRecord(r);
  
  return m;
}

int ReadENFRecord(TFILE *f, ENF_RECORD *r, int swp) {
  int n, m = 0;

  RSF0(r->ilev);
  RSF0(r->energy);
  RSF0(r->pbasis);
  
  if (swp) SwapEndianENFRecord(r);

  return m;
}

int ReadTRHeader(TFILE *f, TR_HEADER *h, int swp) {
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

int ReadTRFHeader(TFILE *f, TRF_HEADER *h, int swp) {
  int n, m = 0;
  
  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->gauge);
  RSF0(h->mode);
  RSF0(h->multipole);
  RSF0(h->efield);
  RSF0(h->bfield);
  RSF0(h->fangle);
  if (swp) SwapEndianTRFHeader(h);

  iuta = 0;
  
  return m;
}

int ReadTRRecord(TFILE *f, TR_RECORD *r, TR_EXTRA *rx, int swp) {
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

int ReadTRFRecord(TFILE *f, TRF_RECORD *r, int swp, TRF_HEADER *h) {
  int n, m = 0, i, nq;
    
  RSF0(r->lower);
  RSF0(r->upper);
  nq = 2*abs(h->multipole) + 1;
  r->strength = (float *) malloc(sizeof(float)*nq);
  RSF1(r->strength, sizeof(float), nq);

  if (swp) {
    SwapEndianTRFRecord(r);
    for (i = 0; i < nq; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }

  return m;
}

int ReadCEHeader(TFILE *f, CE_HEADER *h, int swp) {
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

int ReadCEFHeader(TFILE *f, CEF_HEADER *h, int swp) {
  int i, n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->n_tegrid);
  RSF0(h->n_egrid);
  RSF0(h->te0);
  RSF0(h->efield);
  RSF0(h->bfield);
  RSF0(h->fangle);

  if (swp) SwapEndianCEFHeader(h);
  
  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  RSF1(h->tegrid, sizeof(double), h->n_tegrid);
 
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCEMFHeader(TFILE *f, CEMF_HEADER *h, int swp) {
  int i, n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->n_tegrid);
  RSF0(h->n_egrid);
  RSF0(h->n_thetagrid);
  RSF0(h->n_phigrid);
  RSF0(h->te0);
  RSF0(h->efield);
  RSF0(h->bfield);
  RSF0(h->fangle);

  if (swp) SwapEndianCEMFHeader(h);
  
  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  RSF1(h->tegrid, sizeof(double), h->n_tegrid);
 
  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  h->thetagrid = (double *) malloc(sizeof(double)*h->n_thetagrid);
  RSF1(h->thetagrid, sizeof(double), h->n_thetagrid);

  h->phigrid = (double *) malloc(sizeof(double)*h->n_phigrid);
  RSF1(h->phigrid, sizeof(double), h->n_phigrid);
  
  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_thetagrid; i++) {
      SwapEndian((char *) &(h->thetagrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_phigrid; i++) {
      SwapEndian((char *) &(h->phigrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCERecord(TFILE *f, CE_RECORD *r, int swp, CE_HEADER *h) {
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
  r->params = NULL;
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

int ReadCEFRecord(TFILE *f, CEF_RECORD *r, int swp, CEF_HEADER *h) {
  int i, n, m = 0, m0;
  
  RSF0(r->lower);
  RSF0(r->upper);
  RSF0(r->bethe);
  RSF1(r->born, sizeof(float), 2);

  if (swp) SwapEndianCEFRecord(r);  
  
  m0 = h->n_egrid;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }
  
  return m;
}

int ReadCEMFRecord(TFILE *f, CEMF_RECORD *r, int swp, CEMF_HEADER *h) {
  int i, n, m = 0, m0;
  
  RSF0(r->lower);
  RSF0(r->upper);
  if (swp) SwapEndianCEMFRecord(r);  

  m0 = h->n_thetagrid * h->n_phigrid;
  r->bethe = (float *) malloc(sizeof(float)*m0);
  RSF1(r->bethe, sizeof(float), m0);
  r->born = (float *) malloc(sizeof(float)*(m0+1));
  RSF1(r->born, sizeof(float), m0+1);

  m0 = h->n_egrid*m0;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
    m0 = h->n_thetagrid * h->n_phigrid;
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->bethe[i]), sizeof(float));
    }
    for (i = 0; i <= m0; i++) {
      SwapEndian((char *) &(r->born[i]), sizeof(float));
    }
  }

  return m;
}

int ReadROHeader(TFILE *f, RO_HEADER *h, int swp) {
  int i, n, m = 0;
  
  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);

  if (swp) SwapEndianROHeader(h);
  return m;
}

int ReadCXHeader(TFILE *f, CX_HEADER *h, int swp) {
  int i, n, m = 0;
  
  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF1(h->tgts, sizeof(char), 128);
  RSF0(h->tgtz);
  RSF0(h->tgtm);
  RSF0(h->tgta);
  RSF0(h->tgtb);
  RSF0(h->tgte);
  RSF0(h->tgtx);
  RSF0(h->ldist);
  RSF0(h->te0);
  RSF0(h->ne0);
  
  if (swp) SwapEndianCXHeader(h);
  h->e0 = malloc(sizeof(double)*h->ne0);
  RSF1(h->e0, sizeof(double), h->ne0);
  
  if (swp) {
    for (i = 0; i < h->ne0; i++) {
      SwapEndian((char *) &(h->e0[i]), sizeof(double));
    }
  }
  return m;
}

int ReadRRHeader(TFILE *f, RR_HEADER *h, int swp) {
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

int ReadRORecord(TFILE *f, RO_RECORD *r, int swp) {
  int n, m = 0;
  
  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->n);
  r->nk = malloc(sizeof(int)*r->n);
  r->nq = malloc(sizeof(double)*r->n);
  r->dn = malloc(sizeof(double)*r->n);
  if (swp) SwapEndianRORecord(r);
  RSF1(r->nk, sizeof(int), r->n);
  RSF1(r->nq, sizeof(double), r->n);
  RSF1(r->dn, sizeof(double), r->n);
  if (swp) {
    int i;
    for (i = 0; i < r->n; i++) {
      SwapEndian((char *) &(r->nk[i]), sizeof(int));
      SwapEndian((char *) &(r->nq[i]), sizeof(double));
      SwapEndian((char *) &(r->dn[i]), sizeof(double));
    }
  }
  return m;
}

int ReadCXRecord(TFILE *f, CX_RECORD *r, int swp, CX_HEADER *h) {
  int i, n, m = 0;
    
  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->vnl);
  if (swp) SwapEndianCXRecord(r);
  r->cx = malloc(sizeof(double)*h->ne0);
  RSF1(r->cx, sizeof(double), h->ne0);
  if (swp) {
    for (i = 0; i < h->ne0; i++) {
      SwapEndian((char *)&(r->cx[i]), sizeof(double));
    }
  }
  return m;
}

int ReadRRRecord(TFILE *f, RR_RECORD *r, int swp, RR_HEADER *h) {
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

int ReadAIHeader(TFILE *f, AI_HEADER *h, int swp) {
  int i, n, m = 0;

  if (version_read[DB_AI-1] < 109) return ReadAIHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->emin);
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

int ReadAIMHeader(TFILE *f, AIM_HEADER *h, int swp) {
  int i, n, m = 0;

  if (version_read[DB_AIM-1] < 109) return ReadAIMHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->emin);
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

int ReadAIRecord(TFILE *f, AI_RECORD *r, int swp) {
  int n, m = 0;

  if (version_read[DB_AI-1] < 109) return ReadAIRecordOld(f, r, swp);

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->rate);

  if (swp) SwapEndianAIRecord(r);
  
  return m;
}

int ReadAIMRecord(TFILE *f, AIM_RECORD *r, int swp) {
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

int ReadCIHeader(TFILE *f, CI_HEADER *h, int swp) {
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

int ReadCIRecord(TFILE *f, CI_RECORD *r, int swp, CI_HEADER *h) {
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

int ReadCIMHeader(TFILE *f, CIM_HEADER *h, int swp) {
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

int ReadCIMRecord(TFILE *f, CIM_RECORD *r, int swp, CIM_HEADER *h) {
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

int ReadSPHeader(TFILE *f, SP_HEADER *h, int swp) {
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

int ReadSPRecord(TFILE *f, SP_RECORD *r, SP_EXTRA *rx, int swp) {
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

int ReadRTHeader(TFILE *f, RT_HEADER *h, int swp) {
  int i, n, m = 0;

  if (version_read[DB_RT-1] < 109) return ReadRTHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->ntransitions);
  RSF0(h->iedist);
  RSF0(h->np_edist);
  RSF0(h->eden);
  RSF0(h->ipdist);
  RSF0(h->np_pdist);
  RSF0(h->pden);

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

int ReadRTRecord(TFILE *f, RT_RECORD *r, int swp) {
  int n, m = 0;

  if (version_read[DB_RT-1] < 109) return ReadRTRecordOld(f, r, swp);

  RSF0(r->dir);
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

int ReadDRHeader(TFILE *f, DR_HEADER *h, int swp) {
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

int ReadDRRecord(TFILE *f, DR_RECORD *r, int swp) {
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

TFILE *OpenFileRO(char *fn, F_HEADER *fhdr, int *swp) {
  TFILE *f;
  int n;
  
  f = FOPEN(fn, "r");
  if (f == NULL) return f;
  n = ReadFHeader(f, fhdr, swp);
  if (n == 0) {
    FCLOSE(f);
    return NULL;
  }
  return f;
}

TFILE *OpenFile(char *fn, F_HEADER *fhdr) {
  return OpenFileWTN(fn, fhdr, 0);
}

TFILE *OpenFileWTN(char *fn, F_HEADER *fhdr, int nth) {
  int ihdr;
  TFILE *f;

  ihdr = fhdr->type - 1;

  f = FOPEN(fn, "r+b");
  if (f == NULL) {
    if (fheader[ihdr].nblocks > 0) {
      printf("A single file for one DB type must be used in one session.\n");
      exit(1);
    }
    f = FOPEN(fn, "wb");
    if (f == NULL) {
      printf("cannot open file %s\n", fn);
      exit(1);
    }
  } else {
    if (fheader[ihdr].nblocks == 0) {
      FCLOSE(f);
      f = FOPEN(fn, "wb");
    }
  }

  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    exit(1);
  }

  fheader[ihdr].type = fhdr->type;
  strncpy(fheader[ihdr].symbol, fhdr->symbol, 2);
  fheader[ihdr].atom = fhdr->atom;
  int nr = 0, mr = 0;
  if (nth > 0) {
    nr = nth;
  } else {
#ifdef USE_MPI
    mr = MPIRank(&nr);
#endif
  }
  fhdr->nthreads = nr;
  fheader[ihdr].nthreads = nr;
  WriteFHeader(f, &(fheader[ihdr]));
  return f;
}

int CloseFile(TFILE *f, F_HEADER *fhdr) {
  int ihdr;
 
  ihdr = fhdr->type-1;
  FSEEK(f, 0, SEEK_SET); 
  fheader[ihdr].type = fhdr->type;
  WriteFHeader(f, &(fheader[ihdr]));
  FCLOSE(f);
  return 0;
}

int InitFile(TFILE *f, F_HEADER *fhdr, void *rhdr) {
  EN_HEADER *en_hdr;
  ENF_HEADER *enf_hdr;
  TR_HEADER *tr_hdr;
  TRF_HEADER *trf_hdr;
  CE_HEADER *ce_hdr;
  CEF_HEADER *cef_hdr;
  CEMF_HEADER *cemf_hdr;
  RR_HEADER *rr_hdr;
  AI_HEADER *ai_hdr;
  AIM_HEADER *aim_hdr;
  CI_HEADER *ci_hdr;
  CIM_HEADER *cim_hdr;
  SP_HEADER *sp_hdr;
  RT_HEADER *rt_hdr;
  DR_HEADER *dr_hdr;
  RO_HEADER *ro_hdr;
  CX_HEADER *cx_hdr;
  long int p;
  int ihdr;

  if (f == NULL) return 0;
  
  ihdr = fhdr->type - 1;
  FSEEK(f, 0, SEEK_END);
  p = FTELL(f);
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
  case DB_ENF:
    enf_hdr = (ENF_HEADER *) rhdr;
    memcpy(&enf_header, enf_hdr, sizeof(ENF_HEADER));
    enf_header.position = p;
    enf_header.length = 0;
    enf_header.nele = enf_hdr->nele;
    enf_header.nlevels = 0;
    break;
  case DB_TRF:
    trf_hdr = (TRF_HEADER *) rhdr;
    memcpy(&trf_header, trf_hdr, sizeof(TRF_HEADER));
    trf_header.position = p;
    trf_header.length = 0;
    trf_header.ntransitions = 0;
    break;
  case DB_CEF:
    cef_hdr = (CEF_HEADER *) rhdr;
    memcpy(&cef_header, cef_hdr, sizeof(CEF_HEADER));
    cef_header.position = p;
    cef_header.length = 0;
    cef_header.ntransitions = 0;
    break;
  case DB_CEMF:
    cemf_hdr = (CEMF_HEADER *) rhdr;
    memcpy(&cemf_header, cemf_hdr, sizeof(CEMF_HEADER));
    cemf_header.position = p;
    cemf_header.length = 0;
    cemf_header.ntransitions = 0;
    break;
  case DB_RO:
    ro_hdr = (RO_HEADER *) rhdr;
    memcpy(&ro_header, ro_hdr, sizeof(RO_HEADER));
    ro_header.position = p;
    ro_header.length = 0;
    ro_header.ntransitions = 0;
    break;
  case DB_CX:
    cx_hdr = (CX_HEADER *) rhdr;
    memcpy(&cx_header, cx_hdr, sizeof(CX_HEADER));
    cx_header.position = p;
    cx_header.length = 0;
    cx_header.ntransitions = 0;
    break;
  default:
    break;
  }

  return 0;
}

int DeinitFile(TFILE *f, F_HEADER *fhdr) {
  int n;

  if (f == NULL || fhdr->type <= 0) return 0;

  switch (fhdr->type) {
  case DB_EN:
    FSEEK(f, en_header.position, SEEK_SET);
    if (en_header.length > 0) {
      n = WriteENHeader(f, &en_header);
    }
    break;
  case DB_TR:
    FSEEK(f, tr_header.position, SEEK_SET);
    if (tr_header.length > 0) {
      n = WriteTRHeader(f, &tr_header);
    }
    break;
  case DB_CE:
    FSEEK(f, ce_header.position, SEEK_SET);
    if (ce_header.length > 0) {
      n = WriteCEHeader(f, &ce_header);
    }
    break;
  case DB_RR:
    FSEEK(f, rr_header.position, SEEK_SET);
    if (rr_header.length > 0) {
      n = WriteRRHeader(f, &rr_header);
    }
    break;
  case DB_AI:
    FSEEK(f, ai_header.position, SEEK_SET);
    if (ai_header.length > 0) {
      n = WriteAIHeader(f, &ai_header);
    }
    break;
  case DB_CI:
    FSEEK(f, ci_header.position, SEEK_SET);
    if (ci_header.length > 0) {
      n = WriteCIHeader(f, &ci_header);
    }
    break;
  case DB_SP:
    FSEEK(f, sp_header.position, SEEK_SET);
    if (sp_header.length > 0) {
      n = WriteSPHeader(f, &sp_header);
    }
    break;
  case DB_RT:
    FSEEK(f, rt_header.position, SEEK_SET);
    if (rt_header.length > 0) {
      n = WriteRTHeader(f, &rt_header);
    }
    break;
  case DB_DR:
    FSEEK(f, dr_header.position, SEEK_SET);
    if (dr_header.length > 0) {
      n = WriteDRHeader(f, &dr_header);
    }
    break;
  case DB_AIM:
    FSEEK(f, aim_header.position, SEEK_SET);
    if (aim_header.length > 0) {
      n = WriteAIMHeader(f, &aim_header);
    }
    break;
  case DB_CIM:
    FSEEK(f, cim_header.position, SEEK_SET);
    if (cim_header.length > 0) {
      n = WriteCIMHeader(f, &cim_header);
    }
    break;
  case DB_ENF:
    FSEEK(f, enf_header.position, SEEK_SET);
    if (enf_header.length > 0) {
      n = WriteENFHeader(f, &enf_header);
    }
    break;
  case DB_TRF:
    FSEEK(f, trf_header.position, SEEK_SET);
    if (trf_header.length > 0) {
      n = WriteTRFHeader(f, &trf_header);
    }
    break;
  case DB_CEF:
    FSEEK(f, cef_header.position, SEEK_SET);
    if (cef_header.length > 0) {
      n = WriteCEFHeader(f, &cef_header);
    }
    break;
  case DB_CEMF:
    FSEEK(f, cemf_header.position, SEEK_SET);
    if (cemf_header.length > 0) {
      n = WriteCEMFHeader(f, &cemf_header);
    }
    break;
  case DB_RO:
    FSEEK(f, ro_header.position, SEEK_SET);
    if (ro_header.length > 0) {
      n = WriteROHeader(f, &ro_header);
    }
    break;
  case DB_CX:
    FSEEK(f, cx_header.position, SEEK_SET);
    if (cx_header.length > 0) {
      n = WriteCXHeader(f, &cx_header);
    }
    break;
  default:
    break;
  }
  FFLUSH(f);
  return 0;
}

int PrintTable(char *ifn, char *ofn, int v0) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp, v, vs;

  f1 = FOPEN(ifn, "r");
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
  v = v0%2;
  v0 /= 2;
  if (v0 == 0) {
    vs = fh.nthreads>1;
  } else if (v0 == 1) {
    vs = 0;
  } else {
    vs = 1;
  }
  if (v && (fh.type < DB_SP ||
	    (fh.type > DB_DR && fh.type < DB_ENF)
	    || fh.type >= DB_RO)) {
    if (mem_en_table == NULL) {
      printf("Energy table has not been built in memory.\n");
      goto DONE;
    }
  }

  if (v && fh.type > DB_CIM && fh.type < DB_RO) {
    if (mem_enf_table == NULL) {
      printf("Field dependent energy table has not been built in memory.\n");
      goto DONE;
    }
  }

  if (fh.nthreads > 0) {
    fprintf(f2, "FAC %d.%d.%d[%d]\n",
	    fh.version, fh.sversion, fh.ssversion, fh.nthreads);
  } else {
    fprintf(f2, "FAC %d.%d.%d\n",
	    fh.version, fh.sversion, fh.ssversion);
  }   
  fprintf(f2, "Endian\t= %d\n", (int) CheckEndian(&fh));
  fprintf(f2, "TSess\t= %lu\n", fh.tsession);
  fprintf(f2, "Type\t= %d\n", fh.type);
  fprintf(f2, "Verbose\t= %d\n", v);
  fprintf(f2, "%s Z\t= %5.1f\n", fh.symbol, fh.atom);
  fprintf(f2, "NBlocks\t= %d\n", fh.nblocks);
  fflush(f2);
  switch (fh.type) {
  case DB_EN:
    if (v) {
      fprintf(f2, "E0\t= %-d, %15.8E\n", 
	      iground, (mem_en_table[iground].energy * HARTREE_EV));
    }
    n = PrintENTable(f1, f2, v, 0, swp);
    break;
  case DB_TR:
    n = PrintTRTable(f1, f2, v, vs, swp);
    break;
  case DB_CE:
    n = PrintCETable(f1, f2, v, vs, swp);
    break;
  case DB_RR:
    n = PrintRRTable(f1, f2, v, vs, swp);
    break;
  case DB_AI:
    n = PrintAITable(f1, f2, v, vs, swp);
    break;
  case DB_CI:
    n = PrintCITable(f1, f2, v, vs, swp);
    break;
  case DB_SP:
    n = PrintSPTable(f1, f2, v, vs, swp);
    break;
  case DB_RT:
    n = PrintRTTable(f1, f2, v, vs, swp);
    break;
  case DB_DR:
    n = PrintDRTable(f1, f2, v, vs, swp);
    break;
  case DB_AIM:
    n = PrintAIMTable(f1, f2, v, vs, swp);
    break;
  case DB_CIM:
    n = PrintCIMTable(f1, f2, v, vs, swp);
    break;
  case DB_ENF:
    n = PrintENFTable(f1, f2, v, vs, swp);
    break;
  case DB_TRF:
    n = PrintTRFTable(f1, f2, v, vs, swp);
    break;
  case DB_CEF:
    n = PrintCEFTable(f1, f2, v, vs, swp);
    break;
  case DB_CEMF:
    n = PrintCEMFTable(f1, f2, v, vs, swp);
    break;
  case DB_RO:
    n = PrintROTable(f1, f2, v, vs, swp);
    break;
  case DB_CX:
    n = PrintCXTable(f1, f2, v, vs, swp);
    break;
  default:
    break;
  }

 DONE:
  FCLOSE(f1);
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

int LevelMatchByName(EN_RECORD *r, char *nc, char*cnr, char *cr) {
  return (StrTrimCmp(r->ncomplex, nc) == 0 &&
	  StrTrimCmp(r->sname, cnr) == 0 &&
	  StrTrimCmp(r->name, cr) == 0);
}

void Match2PhotonLevels(int k, EN_RECORD *r, int *ilow2ph, int *iup2ph,
			double *elow2ph, double *eup2ph) {
  if (k == 1) {
    if (LevelMatchByName(r, "1*1", "1s1", "1s+1(1)1")) {
      ilow2ph[0] = r->ilev;
      elow2ph[0] = r->energy;
    } else if (LevelMatchByName(r, "2*1", "2s1", "2s+1(1)1")) {
      iup2ph[0] = r->ilev;
      eup2ph[0] = r->energy;
    }
  } else if (k == 2) {
    if (LevelMatchByName(r, "1*2", "1s2", "1s+2(0)0")) {
      ilow2ph[1] = r->ilev;
      elow2ph[1] = r->energy;
    } else if (LevelMatchByName(r, "1*1.2*1", "1s1.2s1", "1s+1(1)1.2s+1(1)0")) {
      iup2ph[1] = r->ilev;
      eup2ph[1] = r->energy;
    }
  } else if (k == 4) {
    if (LevelMatchByName(r, "1*2.2*2", "2s2", "2s+2(0)0") ||
	LevelMatchByName(r, "1*2.2*2", "1s2.2s2", "1s+2(0)0.2s+2(0)0")) {
      ilow2ph[2] = r->ilev;
      elow2ph[2] = r->energy;
    } else if (LevelMatchByName(r,
				"1*2.2*2", "2s1.2p1",
				"2s+1(1)1.2p-1(1)0") ||
	       LevelMatchByName(r,
				"1*2.2*2", "1s2.2s1.2p1",
				"1s+2(0)0.2s+1(1)1.2p-1(1)0")) {
      iup2ph[2] = r->ilev;
      eup2ph[2] = r->energy;
    }
  }
}

int FindLevelByName(char *fn, int nele, char *nc, char *cnr, char *cr) {
  F_HEADER fh;  
  EN_HEADER h;
  EN_RECORD r;
  TFILE *f;
  int n, k;
  int swp;
  
  f = FOPEN(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) {
    FCLOSE(f);
    return 0;
  }

  if (fh.type != DB_EN) {
    printf("File type is not DB_EN\n");
    FCLOSE(f);
    return -1;
  }

  while (1) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    if (h.nele != nele) {
      FSEEK(f, h.length, SEEK_CUR);
      continue;
    }
    for (k = 0; k < h.nlevels; k++) {
      n = ReadENRecord(f, &r, swp);
      if (StrTrimCmp(r.ncomplex, nc) == 0 &&
	  StrTrimCmp(r.sname, cnr) == 0 &&
	  StrTrimCmp(r.name, cr) == 0) {
	FCLOSE(f);
	return r.ilev;
      }
    }
  }
  
  FCLOSE(f);
  return -1;
}
      
int LevelInfor(char *fn, int ilev, EN_RECORD *r0) {
  F_HEADER fh;  
  EN_HEADER h;
  EN_RECORD r;
  TFILE *f;
  int n, i, k, nlevels;
  int swp, sr;
  
  f = FOPEN(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) {
    FCLOSE(f);
    return 0;
  }
  if (fh.type != DB_EN) {
    printf("File type is not DB_EN\n");
    FCLOSE(f);
    return -1;
  }
  if (version_read[DB_EN-1] < 109) sr = sizeof(EN_RECORD);
  else if (version_read[DB_EN-1] < 115) sr = SIZE_EN_RECORD0;
  else sr = SIZE_EN_RECORD;

  if (ilev >= 0) {
    k = ilev;
    nlevels = 0;
    for (i = 0; i < fh.nblocks; i++) {
      n = ReadENHeader(f, &h, swp);
      if (n == 0) break;
      nlevels += h.nlevels;
      if (k < h.nlevels) {
	if (k > 0) FSEEK(f, sr*k, SEEK_CUR);
	n = ReadENRecord(f, &r, swp);
	if (n == 0) break;
	if (r.ilev != ilev) {
	  FCLOSE(f);
	  return -1;
	}
	memcpy(r0, &r, sizeof(EN_RECORD));
	break;
      } else {
	k -= h.nlevels;
	FSEEK(f, h.length, SEEK_CUR);
      }
    }
    FCLOSE(f);
    if (i == fh.nblocks) return -1;
    return 0;
  } else {
    k = -ilev;
    if (k == 1000) k = 0;
    if (k < 1000) {
      nlevels = 0;
      for (i = 0; i < fh.nblocks; i++) {
	n = ReadENHeader(f, &h, swp);
	if (n == 0) break;
	if (h.nele == k) break;
	nlevels += h.nlevels;
	FSEEK(f, h.length, SEEK_CUR);
      }
      FCLOSE(f);
      if (i == fh.nblocks) return -1;
      return nlevels;
    } else {
      nlevels = 0;
      k -= 1000;
      if (k == 1) {
	for (i = 0; i < fh.nblocks; i++) {
	  n = ReadENHeader(f, &h, swp);
	  if (n == 0) break;
	  nlevels += h.nlevels;
	  FSEEK(f, h.length, SEEK_CUR);
	}
      } else if (k == 2) {
	nlevels = fh.nblocks;
      } else if (k >= 1000) {
	k -= 1000;
	for (i = 0; i < fh.nblocks; i++) {
	  if (i >= k) break;
	  n = ReadENHeader(f, &h, swp);
	  if (n == 0) break;
	  nlevels += h.nlevels;
	  FSEEK(f, h.length, SEEK_CUR);
	}
      }
      FCLOSE(f);
      return nlevels;
    }
  }  
}

EN_SRECORD *GetMemENTable(int *s) {
  *s = mem_en_table_size;
  return mem_en_table;
}

EN_SRECORD *GetMemENFTable(int *s) {
  *s = mem_enf_table_size;
  return mem_enf_table;
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
  TFILE *f;
  char *s;
  int n, i, nlevels;
  float e0;
  int swp, sr;

  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) return -1;

  if (fh.type == DB_ENF) {
    FCLOSE(f);
    return MemENFTable(fn);
  }
  if (fh.type != DB_EN) {
    printf("File type is not DB_EN: %d\n", fh.type);
    FCLOSE(f);
    return -1;
  }
  if (version_read[DB_EN-1] < 109) sr = sizeof(EN_RECORD);
  else if (version_read[DB_EN-1] < 115) sr = SIZE_EN_RECORD0;
  else sr = SIZE_EN_RECORD;

  if (mem_en_table) free(mem_en_table);

  nlevels = 0;
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    if (h.length > sr) {
      FSEEK(f, h.length-sr, SEEK_CUR);
    }
    n = ReadENRecord(f, &r, swp);
    if (r.ilev >= nlevels) nlevels = r.ilev+1;
  }

  mem_en_table = (EN_SRECORD *) malloc(sizeof(EN_SRECORD)*nlevels);
  mem_en_table_size = nlevels;

  e0 = 0.0;
  if (version_read[DB_EN-1] < 109) {
    FSEEK(f, sizeof(F_HEADER), SEEK_SET);
  } else {
    FSEEK(f, SIZE_F_HEADER, SEEK_SET);
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
      mem_en_table[r.ilev].ibase = r.ibase;
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

  FCLOSE(f);
  return 0;
}    

int MemENFTable(char *fn) {
  F_HEADER fh;  
  ENF_HEADER h;
  ENF_RECORD r;
  TFILE *f;
  char *s;
  int n, i, nlevels;
  float e0;
  int swp, sr;

  f = FOPEN(fn, "r");
  if (f == NULL) return -1;

  n = ReadFHeader(f, &fh, &swp);  
  if (n == 0) return 0;

  sr = SIZE_ENF_RECORD;

  if (mem_enf_table) free(mem_enf_table);

  nlevels = 0;
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadENFHeader(f, &h, swp);
    if (n == 0) break;
    if (h.length > sr) {
      FSEEK(f, h.length-sr, SEEK_CUR);
    }
    n = ReadENFRecord(f, &r, swp);
    if (r.ilev >= nlevels) nlevels = r.ilev+1;
  }
  
  mem_enf_table = (EN_SRECORD *) malloc(sizeof(EN_SRECORD)*nlevels);
  mem_enf_table_size = nlevels;

  FSEEK(f, SIZE_F_HEADER, SEEK_SET);
  while (1) {
    n = ReadENFHeader(f, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.nlevels; i++) {
      n = ReadENFRecord(f, &r, swp);
      if (n == 0) break;
      mem_enf_table[r.ilev].energy = r.energy;
      DecodeBasisEB(r.pbasis, &mem_enf_table[r.ilev].p, &mem_enf_table[r.ilev].j);
    }
  }

  FCLOSE(f);

  return 0;
}    

int PrintENTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
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
      fprintf(f2, "%6d %6d %15.8E %1d %5d %4d %-32s %-48s %-s\n",
	      r.ilev, ibase, e, p, vnl, j, r.ncomplex, r.sname, r.name);
    }
    nb += 1;
  }
  
  return nb;
}

int CodeBasisEB(int s, int m) {
  int k;
  
  k = s + MAXLEVEB*abs(m);
  if (m < 0) k = -k;

  return k;
}

void DecodeBasisEB(int k, int *s, int *m) {
  *m = abs(k)/MAXLEVEB;
  *s = abs(k)%MAXLEVEB;
  if (k < 0) *m = -(*m);
}

int PrintENFTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  ENF_HEADER h;
  ENF_RECORD r;
  int n, i;
  int nb, ilev, mlev;
  double e;
  
  nb = 0;
  while (1) {
    n = ReadENFHeader(f1, &h, swp);
    if (n == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NLEV\t= %d\n", h.nlevels);
    fprintf(f2, "EFIELD\t= %15.8E\n", h.efield);
    fprintf(f2, "BFIELD\t= %15.8E\n", h.bfield);
    fprintf(f2, "FANGLE\t= %15.8E\n", h.fangle);
    fprintf(f2, "%6s %22s %6s %4s\n", "ILEV", "ENERGY", "PBASIS", "M");
    for (i = 0; i < h.nlevels; i++) {
      n = ReadENFRecord(f1, &r, swp);
      if (n == 0) break;
      e = r.energy;
      if (v) {
	e -= mem_en_table[iground].energy;
	e *= HARTREE_EV;
      }
      DecodeBasisEB(r.pbasis, &ilev, &mlev);
      fprintf(f2, "%6d %22.15E %6d %4d\n", r.ilev, e, ilev, mlev);
    }
    nb++;
  }

  return nb;
}
    
double OscillatorStrength(int m, double e, double s, double *ga) {
  int m2;
  double aw, x;

  aw = FINE_STRUCTURE_CONST * e;
  if (itrf == 0 && m != 0) {
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

int PrintTRTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  TR_HEADER h;
  TR_RECORD r;
  TR_EXTRA rx;
  int n, i;
  int nb, nh;
  double e, a, gf;

  nb = 0;
  
  while (1) {
    nh = ReadTRHeader(f1, &h, swp);
    if (nh == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "MULTIP\t= %d\n", (int)h.multipole);
    fprintf(f2, "GAUGE\t= %d\n", (int)h.gauge);
    fprintf(f2, "MODE\t= %d\n", (int)h.mode);

    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadTRRecord(f1, &r, &rx, swp);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.upper;
	idx[i].i1 = r.lower;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
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
	  fprintf(f2, "%5d %5d %13.6E %11.4E %13.6E %10.3E\n",
		  r.upper, r.lower, e, rx.sdev, r.strength, rx.sci);
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
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb += 1;
  }

  return nb;
}

int PrintTRFTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  TRF_HEADER h;
  TRF_RECORD r;
  int n, nh, i, j;
  int nb, nq;
  double e, a, gf, ta;

  nb = 0;
  
  while (1) {
    nh = ReadTRFHeader(f1, &h, swp);
    if (nh == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "MULTIP\t= %d\n", (int)h.multipole);
    fprintf(f2, "GAUGE\t= %d\n", (int)h.gauge);
    fprintf(f2, "MODE\t= %d\n", (int)h.mode);
    fprintf(f2, "EFIELD\t= %15.8E\n", h.efield);
    fprintf(f2, "BFIELD\t= %15.8E\n", h.bfield);
    fprintf(f2, "FANGLE\t= %15.8E\n", h.fangle);
    nq = 2*abs(h.multipole)+1;
    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadTRFRecord(f1, &r, swp, &h);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.upper;
	idx[i].i1 = r.lower;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
      n = ReadTRFRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (v) {
	e = mem_enf_table[r.upper].energy - mem_enf_table[r.lower].energy;
	ta = 0.0;
	for (j = 0; j < nq; j++) {
	  gf = OscillatorStrength(h.multipole, e, (double)(r.strength[j]), &a);
	  a *= RATE_AU;
	  ta += a;
	  fprintf(f2, "%6d %6d %3d %6d %6d %3d %2d %13.6E %13.6E %13.6E %13.6E %13.6E\n",
		  r.upper, mem_enf_table[r.upper].p, mem_enf_table[r.upper].j,
		  r.lower, mem_enf_table[r.lower].p, mem_enf_table[r.lower].j,
		  j-abs(h.multipole),(e*HARTREE_EV), gf, a, r.strength[j], ta);
	}
      } else {
	for (j = 0; j < nq; j++) {
	  fprintf(f2, "%6d %6d %13.6E\n", 
		  r.upper, r.lower, r.strength[j]);
	}
      }
      free(r.strength);
    }
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
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
  TFILE *f;
  int n, i, k;
  double a, b, c, e;
  int swp;
 
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) {
    return -1;
  }
  if (fh.type != DB_TR) {
    printf("File type is not DB_TR\n");
    FCLOSE(f);
    return -1;
  }
  
  a = 0.0;
  c = 0.0;
  int nt = 0;
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
	nt++;
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

  FCLOSE(f);

  return nt;
}
  
int PrintCETable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  CE_HEADER h;
  CE_RECORD r;
  int n, i, t, nb, nh;
  int m, k, p1, p2;
  float a, e;
  double bte, bms, be;

  nb = 0;
  BornFormFactorTE(&bte);
  bms = BornMass();  
  while (1) {
    nh = ReadCEHeader(f1, &h, swp);
    if (nh == 0) break;

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
    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadCERecord(f1, &r, swp, &h);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.lower;
	idx[i].i1 = r.upper;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
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
      
      be = (e + bte)/bms;
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
	    if (h.usr_egrid_type == 1) a += be;
	    a *= 2.0*(1.0 + 0.5*FINE_STRUCTURE_CONST2 * a);
	    a = PI * AREA_AU20/a;
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
      fflush(f2);
      if (h.msub || h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb += 1;
  }

  return nb;
}
  
int PrintCEFTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  CEF_HEADER h;
  CEF_RECORD r;
  int n, i, t;
  int nb, nh;
  int m;
  float a, e;
  double bte, bms, be;

  nb = 0;
 
  BornFormFactorTE(&bte);
  bms = BornMass();  
  while (1) {
    nh = ReadCEFHeader(f1, &h, swp);
    if (nh == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);

    for (i = 0; i < h.n_tegrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "TE0\t= %15.8E\n", h.te0 * HARTREE_EV);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }

    fprintf(f2, "EFIELD\t= %15.8E\n", h.efield);
    fprintf(f2, "BFIELD\t= %15.8E\n", h.bfield);
    fprintf(f2, "FANGLE\t= %15.8E\n", h.fangle);

    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadCEFRecord(f1, &r, swp, &h);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.lower;
	idx[i].i1 = r.upper;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
      n = ReadCEFRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (v) {
	e = mem_enf_table[r.upper].energy - mem_enf_table[r.lower].energy;
	fprintf(f2, "%6d %6d %3d %6d %2d %3d %11.4E\n",
		r.lower, mem_enf_table[r.lower].p, mem_enf_table[r.lower].j,
		r.upper, mem_enf_table[r.upper].p, mem_enf_table[r.upper].j,
		e*HARTREE_EV);
	fprintf(f2, "%11.4E %11.4E %11.4E\n", 
		r.bethe, r.born[0], r.born[1]*HARTREE_EV);
      } else {
	fprintf(f2, "%6d %6d\n", 
		r.lower, r.upper);
	fprintf(f2, "%11.4E %11.4E %11.4E\n", 
		r.bethe, r.born[0], r.born[1]);
      }
      
      be = (e + bte)/bms;
      for (t = 0; t < h.n_egrid; t++) {
	if (v) {
	  a = h.egrid[t];
	  a += be;
	  a *= 2.0*(1.0 + 0.5*FINE_STRUCTURE_CONST2 * a);
	  a = PI * AREA_AU20/a;
	  a *= r.strength[t];
	    fprintf(f2, "%11.4E %11.4E %11.4E\n",
		    h.egrid[t]*HARTREE_EV,
		    r.strength[t], a);
	} else {
	  fprintf(f2, "%11.4E %11.4E\n", h.egrid[t], r.strength[t]);
	}
      }      
      free(r.strength);
    }
    free(h.tegrid);
    free(h.egrid);
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb += 1;
  }

  return nb;
}
  
int PrintCEMFTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  CEMF_HEADER h;
  CEMF_RECORD r;
  int n, i, t;
  int nb, nh;
  int m, k, na, ith, iph;
  float a, e;
  double bte, bms, be;

  nb = 0;
 
  BornFormFactorTE(&bte);
  bms = BornMass(); 
  while (1) {
    nh = ReadCEMFHeader(f1, &h, swp);
    if (nh == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);

    for (i = 0; i < h.n_tegrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "TE0\t= %15.8E\n", h.te0 * HARTREE_EV);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    
    fprintf(f2, "NTHETA\t= %d\n", h.n_thetagrid);
    for (i = 0; i < h.n_thetagrid; i++) {
      fprintf(f2, "\t %15.8E\n", h.thetagrid[i]);
    }
    fprintf(f2, "NPHI\t= %d\n", h.n_phigrid);
    for (i = 0; i < h.n_phigrid; i++) {
      fprintf(f2, "\t %15.8E\n", h.phigrid[i]);
    }
    
    fprintf(f2, "EFIELD\t= %15.8E\n", h.efield);
    fprintf(f2, "BFIELD\t= %15.8E\n", h.bfield);
    fprintf(f2, "FANGLE\t= %15.8E\n", h.fangle);

    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadCEMFRecord(f1, &r, swp, &h);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.lower;
	idx[i].i1 = r.upper;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    na = h.n_thetagrid * h.n_phigrid;
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
      n = ReadCEMFRecord(f1, &r, swp, &h);
      if (n == 0) break;
      k = 0;
      for (ith = 0; ith < h.n_thetagrid; ith++) {
	for (iph = 0; iph < h.n_phigrid; iph++) {	    
	  if (v) {
	    e = mem_enf_table[r.upper].energy - mem_enf_table[r.lower].energy;
	    fprintf(f2, "%6d %6d %3d %6d %2d %3d %11.4E\n",
		    r.lower, mem_enf_table[r.lower].p, mem_enf_table[r.lower].j,
		    r.upper, mem_enf_table[r.upper].p, mem_enf_table[r.upper].j,
		    e*HARTREE_EV);
	    fprintf(f2, "%11.4E %11.4E %11.4E %11.4E %11.4E\n", 
		    h.thetagrid[ith]*180.0/PI, h.phigrid[iph]*180.0/PI,
		    r.bethe[k], r.born[k], r.born[na]*HARTREE_EV);
	  } else {
	    fprintf(f2, "%6d %6d\n", 
		    r.lower, r.upper);
	    fprintf(f2, "%11.4E %11.4E %11.4E\n", 
		    r.bethe[k], r.born[k], r.born[na]);
	  }      
	  be = (e + bte)/bms;
	  for (t = 0; t < h.n_egrid; t++) {
	    if (v) {
	      a = h.egrid[t];
	      a += be;
	      a *= 2.0*(1.0 + 0.5*FINE_STRUCTURE_CONST2 * a);
	      a = PI * AREA_AU20/a;
	      a *= r.strength[t+h.n_egrid*k];
	      fprintf(f2, "%11.4E %11.4E %11.4E\n",
		      h.egrid[t]*HARTREE_EV,
		      r.strength[t+h.n_egrid*k], a);
	    } else {
	      fprintf(f2, "%11.4E %11.4E\n", h.egrid[t], r.strength[t+h.n_egrid*k]);
	    }
	  }  
	  k++;
	}
      }    
      free(r.strength);
      free(r.bethe);
      free(r.born);
    }
    free(h.tegrid);
    free(h.egrid);
    free(h.thetagrid);
    free(h.phigrid);
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb += 1;
  }

  return nb;
}

int PrintROTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  RO_HEADER h;
  RO_RECORD r;
  int n, i, nb, nh, t;
  double e;
  
  nb = 0;
  while (1) {
    nh = ReadROHeader(f1, &h, swp);
    if (nh == 0) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);

    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadRORecord(f1, &r, swp);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.f;
	idx[i].i1 = r.b;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
      n = ReadRORecord(f1, &r, swp);
      if (n == 0) break;

      if (v) {	
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	for (t = 0; t < r.n; t++) {
	  fprintf(f2, "%6d %6d %15.8E %6d %6d %6d %12.5E %12.5E\n",
		  r.b, r.f, e*HARTREE_EV, t, r.n, r.nk[t], r.nq[t], r.dn[t]);
	}
      } else {
	for (t = 0; t < r.n; t++) {
	  fprintf(f2, "%6d %6d %6d %12.5E %12.5E\n",
		  r.b, r.f, r.nk[t], r.nq[t], r.dn[t]);
	}
      }
      free(r.nk);
      free(r.nq);
    }
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb++;
  }
  return nb;
}

int PrintCXTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  CX_HEADER h;
  CX_RECORD r;
  int n, i, nb, nh, t;
  double e;
  
  nb = 0;
  while (1) {
    nh = ReadCXHeader(f1, &h, swp);
    if (nh == 0) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "TGTS\t= %s\n", h.tgts);
    fprintf(f2, "TGTZ\t= %g\n", h.tgtz);
    fprintf(f2, "TGTM\t= %g\n", h.tgtm);
    fprintf(f2, "TGTA\t= %g\n", h.tgta);
    fprintf(f2, "TGTB\t= %g\n", h.tgtb);
    fprintf(f2, "TGTE\t= %g\n", h.tgte);
    fprintf(f2, "TGTX\t= %g\n", h.tgtx);
    fprintf(f2, "LDIST\t= %d\n", h.ldist);
    fprintf(f2, "TE0\t= %d\n", h.te0);
    fprintf(f2, "NE0\t= %d\n", h.ne0);      
    for (i = 0; i < h.ne0; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.e0[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.e0[i]);
      }
    }
    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadCXRecord(f1, &r, swp, &h);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.f;
	idx[i].i1 = r.b;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
      n = ReadCXRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (v) {	
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "%6d %2d %6d %2d %4d %15.8E\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		r.vnl, e*HARTREE_EV);
	for (t = 0; t < h.ne0; t++) {
	  fprintf(f2, "%12.5E %12.5E\n",
		  h.e0[t]*HARTREE_EV, r.cx[t]*AREA_AU20);
	}
      } else {
	fprintf(f2, "%6d %6d %4d\n", r.b, r.f, r.vnl);
	for (t = 0; t < h.ne0; t++) {
	  fprintf(f2, "%12.5E %12.5E\n", h.e0[t], r.cx[t]);
	}
      }
      free(r.cx);
    }
    free(h.e0);
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb++;
  }
  return nb;
}

int PrintRRTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  RR_HEADER h;
  RR_RECORD r;
  int n, i, t;
  int nb, nh, k, m;
  float e, eph, ee, phi, rr;

  nb = 0;
  while (1) {
    nh = ReadRRHeader(f1, &h, swp);
    if (nh == 0) break;
    
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
    
    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadRRRecord(f1, &r, swp, &h);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.f;
	idx[i].i1 = r.b;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
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
	  phi = 2.0*PI*FINE_STRUCTURE_CONST*r.strength[t]*AREA_AU20;
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
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb++;
  }

  return nb;
}

int PrintAITable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  AI_HEADER h;
  AI_RECORD r;
  int n, i;
  int nb, nh;
  float e, sdr, er;
  
  nb = 0;
  
  while (1) {
    nh = ReadAIHeader(f1, &h, swp);
    if (nh == 0) break;
 
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "EMIN\t= %15.8E\n", h.emin*HARTREE_EV);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
       
    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadAIRecord(f1, &r, swp);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.b;
	idx[i].i1 = r.f;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
      n = ReadAIRecord(f1, &r, swp);
      if (n == 0) break;
      if (v) {
	e = mem_en_table[r.b].energy - mem_en_table[r.f].energy;
	if (e < 0) er = e - h.emin;
	else er = e;
	sdr = 0.5*(mem_en_table[r.b].j + 1.0);
	sdr *= PI*PI*r.rate/(er*(mem_en_table[r.f].j + 1.0));
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
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb++;
  }

  return nb;
}

int AIBranch(char *fn, int ib, int ia,
	     double *te, double *pa, double *ta) {
  F_HEADER fh;
  AI_HEADER h;
  AI_RECORD r;
  TFILE *f;
  int n, i, k, nt;
  double a, b, c, e;
  int swp;
    
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  if (fh.type != DB_AI) {
    printf("File type is not DB_AI\n");
    FCLOSE(f);
    return -1;
  }
   
  a = 0.0;
  c = 0.0;
  nt = 0;
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
	nt++;
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

  FCLOSE(f);
  
  return nt;
}
      
int PrintAIMTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  AIM_HEADER h;
  AIM_RECORD r;
  int n, i, m;
  int nb, nh;
  float e;
  double u = AREA_AU20*HARTREE_EV;
  
  nb = 0;
  while (1) {
    nh = ReadAIMHeader(f1, &h, swp);
    if (nh == 0) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "EMIN\t= %15.8E\n", h.emin*HARTREE_EV);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    
    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadAIMRecord(f1, &r, swp);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.b;
	idx[i].i1 = r.f;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
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
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb++;
  }

  return nb;
}

int PrintCITable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  CI_HEADER h;
  CI_RECORD r;
  int n, i, t;
  int nb, nh, m;
  float e, a;
  double bte, bms, be;

  nb = 0;
 
  BornFormFactorTE(&bte);
  bms = BornMass(); 
  while (1) {
    nh = ReadCIHeader(f1, &h, swp);
    if (nh == 0) break;
    
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

    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadCIRecord(f1, &r, swp, &h);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.b;
	idx[i].i1 = r.f;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
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
      be = (e + bte)/bms;
      for (t = 0; t < h.n_usr; t++) {
	if (v) {
	  a = h.usr_egrid[t];
	  if (h.usr_egrid_type == 1) a += be;
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
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb++;
  }

  return nb;
}

int PrintCIMTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  CIM_HEADER h;
  CIM_RECORD r;
  int n, i, t, q;
  int nb, nh, m, k;
  float e, a;
  double bte, bms, be;

  nb = 0;
 
  BornFormFactorTE(&bte);
  bms = BornMass(); 
  while (1) {
    nh = ReadCIMHeader(f1, &h, swp);
    if (nh == 0) break;
    
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

    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadCIMRecord(f1, &r, swp, &h);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.b;
	idx[i].i1 = r.f;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }

    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
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

      be = (e + bte)/bms;
      if (v) {
	q = 0;
	for (k = 0; k < r.nsub; k ++) {
	  for (t = 0; t < h.n_usr; t++) {
	    a = h.usr_egrid[t];
	    if (h.usr_egrid_type == 1) a += be;
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
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb++;
  }

  return nb;
}

int PrintSPTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  SP_HEADER h;
  SP_RECORD r;
  SP_EXTRA rx;
  int n, i;
  int nb, nh;
  float e, a;

  nb = 0;
  
  while (1) {
    nh = ReadSPHeader(f1, &h, swp);
    if (nh == 0) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "TYPE\t= %07d\n", h.type);
    fprintf(f2, "IBLK\t= %d\n", h.iblock);
    fprintf(f2, "ICOMP\t= %s\n", h.icomplex);
    fprintf(f2, "FBLK\t= %d\n", h.fblock);
    fprintf(f2, "FCOMP\t= %s\n", h.fcomplex);

    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadSPRecord(f1, &r, &rx, swp);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.upper;
	idx[i].i1 = r.lower;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
      n = ReadSPRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      e = r.energy;
      if (v) e *= HARTREE_EV;
      a = r.strength;
      if (iuta) {
	fprintf(f2, "%6d %6d %13.6E %11.4E %11.4E %11.4E %11.4E\n", 
		r.upper, r.lower, e, rx.sdev*HARTREE_EV, a, r.rrate, r.trate);
      } else {
	fprintf(f2, "%6d %6d %13.6E %11.4E %11.4E %11.4E\n", 
		r.upper, r.lower, e, a, r.rrate, r.trate);
      }
    }
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
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
  TFILE *f;
  int n, swp;
  double d;

  f = FOPEN(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1.0;
  }
  
  n = ReadFHeader(f, &fh, &swp);
  
  if (fh.type != DB_SP) {
    printf("file not of type DB_SP\n");
    FCLOSE(f);
    return -1.0;
  }

  d = 0.0;
  while (1) {
    n = ReadSPHeader(f, &h, swp);
    if (n == 0) break;
    if (h.type != 0 || h.nele != k) {
      FSEEK(f, h.length, SEEK_CUR);
      continue;
    }
    
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f, &r, &rx, swp);
      if (n == 0) break;
      d += r.strength;
    }
  }
 
  FCLOSE(f);

  return d;
}

double IonRadiation(char *fn, int k, int m) {
  F_HEADER fh;
  SP_HEADER h;
  SP_RECORD r;
  SP_EXTRA rx;
  int i;
  TFILE *f;
  int n, swp;
  double d;

  f = FOPEN(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1.0;
  }
  
  n = ReadFHeader(f, &fh, &swp);
  
  if (fh.type != DB_SP) {
    printf("file not of type DB_SP\n");
    FCLOSE(f);
    return -1.0;
  }

  d = 0.0;
  while (1) {
    n = ReadSPHeader(f, &h, swp);
    if (n == 0) break;
    if (h.type == 0 || h.nele != k) {
      FSEEK(f, h.length, SEEK_CUR);
      continue;
    }      
    if (m == 1 && h.type < 100) {
      FSEEK(f, h.length, SEEK_CUR);
      continue;
    }
    if (m == 2 && h.type > 100) {
      FSEEK(f, h.length, SEEK_CUR);
      continue;
    }
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f, &r, &rx, swp);
      if (n == 0) break;
      d += r.strength*r.energy;
    }
  }
 
  FCLOSE(f);

  d *= HARTREE_EV;
  return d;
}
		     
int PrintRTTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  RT_HEADER h;
  RT_RECORD r;
  int n, i;
  int nb, nh, nele;

  nb = 0;
  nele = -1;
  while (1) {
    nh = ReadRTHeader(f1, &h, swp);
    if (nh == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
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

    fprintf(f2,"              NB          TR          CE");
    fprintf(f2, "          RR          AI          CI\n");

    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadRTRecord(f1, &r, swp);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	if (r.dir >= 0) {
	  idx[i].i0 = 0;
	} else {
	  idx[i].i0 = 1;
	}
	idx[i].i1 = r.iblock;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
      n = ReadRTRecord(f1, &r, swp);
      if (n == 0) break;
      if (r.ci < 0) {
	r.ce *= HARTREE_EV;
      }
      fprintf(f2, "%6d %4d  %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %s\n",
	      r.dir, r.iblock, r.nb, r.tr, r.ce, r.rr, r.ai, r.ci, r.icomplex);
    }
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb += 1;
  }

  return nb;
}

int PrintDRTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  DR_HEADER h;
  DR_RECORD r;
  int n, i;
  int nb, nh;
  double e, e1;
  
  nb = 0;
  while (1) {
    nh = ReadDRHeader(f1, &h, swp);
    if (nh == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "ILEV\t= %d\n", h.ilev);
    e = h.energy;
    if (v) e *= HARTREE_EV;
    fprintf(f2, "E\t= %15.8E\n", e);
    fprintf(f2, "JLEV\t= %d\n", h.j);
    fprintf(f2, "NREC\t= %d\n", h.vn);

    IDX_RECORD *idx = NULL;
    int hp = 0;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadDRRecord(f1, &r, swp);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.ilev;
	idx[i].i1 = r.ibase;
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
      n = ReadDRRecord(f1, &r, swp);
      if (n == 0) break;
      e = r.energy;
      e1 = r.etrans;
      if (v) {
	e *= HARTREE_EV;
	e1 *= HARTREE_EV;
	fprintf(f2, "%6d %2d %4d %2d %3d %4d %3d %2d %2d %11.4E %10.4E %10.4E %10.4E %10.4E\n",
		r.ilev, r.j, h.ilev, h.j, r.ibase, r.flev, r.fbase, 
		h.vn, r.vl, e, e1, r.ai, r.total_rate, r.br);
      } else {
	fprintf(f2, "%6d %2d %3d %4d %3d %2d %11.4E %10.4E %10.4E %10.4E %10.4E\n",
		r.ilev, r.j, r.ibase, r.flev, r.fbase, r.vl, 
		e, e1, r.ai, r.total_rate, r.br);
      }
    }
    if (idx) {
      free(idx);
      FSEEK(f1, h.position+nh+h.length, SEEK_SET);
    }
    nb++;
  }

  return nb;
}

int AppendTable(char *fn) {
  F_HEADER fh;
  TFILE *f;
  int n, swp;
    
  f = FOPEN(fn, "r");
  if (f == NULL) return -1;
  n = ReadFHeader(f, &fh, &swp);
  if (swp) {
    printf("File %s is in different byte-order\n", fn);
    FCLOSE(f);
    return -1;
  }
  memcpy(&(fheader[fh.type-1]), &fh, sizeof(F_HEADER));
  FCLOSE(f);
  
  return 0;
}
  
int JoinTable(char *fn1, char *fn2, char *fn) {
  F_HEADER fh1, fh2;
  TFILE *f1, *f2, *f;
  int n, swp1, swp2;
#define NBUF 8192
  char buf[NBUF];

  f1 = FOPEN(fn1, "r");
  if (f1 == NULL) return -1;
  f2 = FOPEN(fn2, "r");
  if (f2 == NULL) return -1;

  n = ReadFHeader(f1, &fh1, &swp1);
  if (n == 0) {
    FCLOSE(f1);
    FCLOSE(f2);
    return 0;
  }
  n = ReadFHeader(f2, &fh2, &swp2);
  if (n == 0) {
    FCLOSE(f1);
    FCLOSE(f2);
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

  f = FOPEN(fn, "w");
  if (f == NULL) return -1;
  fh1.nblocks += fh2.nblocks;
  
  WriteFHeader(f, &fh1);
  while (1) {
    n = FREAD(buf, 1, NBUF, f1);
    if (n > 0) {
      if (n > FWRITE(buf, 1, n, f)) {
	printf("write error\n");
	return -1;
      }
    }
    if (n < NBUF) break;
  }
  while (1) {
    n = FREAD(buf, 1, NBUF, f2);
    if (n > 0) {
      if (n > FWRITE(buf, 1, n, f)) {
	printf("write error\n");
	return -1;
      }
    }
    if (n < NBUF) break;
  }

  FCLOSE(f1);
  FCLOSE(f2);
  FCLOSE(f);
  
  return 0;
#undef NBUF
}

int ISearch(int i, int n, int *ia) {
  int k;

  for (k = 0; k < n; k++) {
    if (ia[k] == i) {
      return k;
    }
  }

  return -1;
}

int AdjustEnergy(int nlevs, int *ilevs, double *e, 
		 char *efn0, char *efn1, char *afn0, char *afn1) {
  int i, k, k0, k1, n, ig, swp, nb;
  double ae0, ae1, e0, e1;
  TFILE *f0, *f1;
  F_HEADER efh, afh;
  AI_HEADER ah;
  AI_RECORD ar;
  EN_HEADER eh;
  EN_RECORD er;

  ig = -1;
  for (i = 0; i < nlevs; i++) {
    if (e[i] == 0) {
      ig = ilevs[i];
      break;
    }
  }

  MemENTable(efn0);
  if (ig >= 0) {
    for (i = 0; i < nlevs; i++) {
      e[i] = e[i] - (mem_en_table[ilevs[i]].energy - mem_en_table[ilevs[ig]].energy);
    }
  }
  f0 = FOPEN(efn0, "r");
  if (f0 == NULL) return -1;

  n = ReadFHeader(f0, &efh, &swp);  
  if (n == 0) {
    FCLOSE(f0);
    return -1;
  }
  
  f1 = OpenFile(efn1, &efh);
  
  while (1) {
    n = ReadENHeader(f0, &eh, swp);
    if (n == 0) break;
    InitFile(f1, &efh, &eh);
    for (i = 0; i < eh.nlevels; i++) {
      n = ReadENRecord(f0, &er, swp);
      if (n == 0) break;
      k = ISearch(er.ilev, nlevs, ilevs);
      if (k >= 0) {
	er.energy += e[k];
      }
      k = ISearch(er.ibase, nlevs, ilevs);
      if (k >= 0) {
	er.energy += e[k];
      }
      WriteENRecord(f1, &er);
    }
    DeinitFile(f1, &efh);
  }

  CloseFile(f1, &efh);
  FCLOSE(f0);

  f0 = FOPEN(afn0, "r");
  if (f0 == NULL) return -1;
  
  n = ReadFHeader(f0, &afh, &swp);
  if (n == 0) {
    FCLOSE(f0);
    return -1;
  }
  
  f1 = OpenFile(afn1, &afh);

  nb = 0;
  while (1) {
    n = ReadAIHeader(f0, &ah, swp);
    if (n == 0) break;
    /*
    if (nb == 0) printf("EMIN = %10.3E\n", ah.emin*HARTREE_EV);
    nb++;
    */
    InitFile(f1, &afh, &ah);
    for (i = 0; i < ah.ntransitions; i++) {
      n = ReadAIRecord(f0, &ar, swp);
      if (n == 0) break;
      ae0 = mem_en_table[ar.b].energy - mem_en_table[ar.f].energy;
      if (ae0 < 0) ae0 -= ah.emin;
      if (mem_en_table[ar.b].ibase >= 0) {
	e0 = 0.0;
	k0 = ISearch(ar.b, nlevs, ilevs);
	if (k0 >= 0) e0 += e[k0];
	k0 = ISearch(mem_en_table[ar.b].ibase, nlevs, ilevs);
	if (k0 >= 0) e0 += e[k0];
	e1 = 0.0;
	k1 = ISearch(ar.f, nlevs, ilevs);
	if (k1 >= 0) e1 += e[k1];
	k1 = ISearch(mem_en_table[ar.f].ibase, nlevs, ilevs);
	if (k1 >= 0) e1 += e[k1];      
	ae1 = ae0 + (e0 - e1);
      } else {
	ae1 = ae0;
      }
      if (ae1 < 0) ae1 -= ah.emin;	
      if (ae1 > 0) {
	ar.rate *= sqrt(ae0/ae1);
      } else {
	ar.rate = 0.0;
      }
      WriteAIRecord(f1, &ar);
    }
    DeinitFile(f1, &afh);
  }

  CloseFile(f1, &afh);
  FCLOSE(f0);

  ReinitDBase(0);

  return 0;
}

int ReadJJLSJ(char *fn, JJLSJ **lsj) {
  int i, p, tj, m, n, nmax, nf;
  int *ts, *tk;
  double *w;
  FILE *f;
  JJLSJ *pt;
  char buf[8192];

  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open LSJ file: %s\n", fn);
    return 0;
  }
  n = -1;
  nmax = 0;
  m = 0;
  while (1) {
    if (NULL == fgets(buf, 8192, f)) break;
    if (buf[0] == '#') {
      if (buf[1] == '#') continue;
      if (nmax < m) nmax = m;
      m = 0;
      n++;
    } else {
      if (strlen(buf) > 5) {
	m++;
      }
    }
  }
  if (nmax < m) nmax = m;
  n++;
  pt = malloc(sizeof(JJLSJ)*n);
  ts = malloc(sizeof(int)*nmax);
  tk = malloc(sizeof(int)*nmax);
  w = malloc(sizeof(double)*nmax);  
  fseek(f, 0, SEEK_SET);
  n = -1;
  m = 0;
  while (1) {
    if (NULL == fgets(buf, 8192, f)) break;
    if (buf[0] == '#') {
      if (buf[1] == '#') continue;
      if (m > 0) {
	pt[n].nks = m;
	pt[n].k = malloc(sizeof(int)*m);
	pt[n].s = malloc(sizeof(int)*m);
	pt[n].w = malloc(sizeof(double)*m);
	for (i = 0; i < m; i++) {
	  pt[n].k[i] = tk[i];
	  pt[n].s[i] = ts[i];
	  pt[n].w[i] = w[i];
	}
      }
      m = 0;
      n++;
    } else {
      nf = sscanf(buf, "%d %d %d %d %d %lg",
		  &i, &p, &tj, &ts[m], &tk[m], &w[m]);
      if (nf != 6) continue;
      if (m == 0) {
	pt[n].ilev = i;
	pt[n].j = tj;
      }
      m++;
    }
  }
  if (m > 0) {
    pt[n].nks = m;
    pt[n].k = malloc(sizeof(int)*m);
    pt[n].s = malloc(sizeof(int)*m);
    pt[n].w = malloc(sizeof(double)*m);
    for (i = 0; i < m; i++) {
      pt[n].k[i] = tk[i];
      pt[n].s[i] = ts[i];
      pt[n].w[i] = w[i];
    }  
    n++;
  }
  fclose(f);
  free(ts);
  free(tk);
  free(w);
  *lsj = pt;
  return n;
}

void RecoupleRORecord(RO_RECORD *r0, RO_RECORD *r1) {
  int i, k, kk, k2, j2, nmax, ssm, nn, n, tt;
  int ss, nk, t, jb, jf;
  int *nid;
  double a, **nq, **dn;
  
  r1->n = 0;
  jb = mem_en_table[r0->b].j;
  jf = mem_en_table[r0->f].j;
  nmax = 0;
  for (k = 0; k < r0->n; k++) {
    n = abs(r0->nk[k])/100;
    if (n > nmax) nmax = n;
  }
  ssm = 2;
  nq = malloc(sizeof(double *)*ssm);
  dn = malloc(sizeof(double *)*ssm);
  nk = nmax*(nmax+1)/2;
  for (k = 0; k < ssm; k++) {
    nq[k] = NULL;
    dn[k] = NULL;
  }
  nid = malloc(sizeof(int)*nk);
  k = 0;
  for (n = 1; n <= nmax; n++) {
    for (t = 0; t < n; t++) {
      nid[k++] = n;
    }
  }
  for (t = 0; t < r0->n; t++) {
    n = abs(r0->nk[t])/100;
    k2 = 2*(abs(r0->nk[t])%100);
    if (r0->nk[t] < 0) j2 = k2-1;
    else j2 = k2+1;
    for (ss = -1; ss <= 1; ss += 2) {
      k = jf + ss;
      if (k < 0) continue;
      a = W6j(jf, 1, k, k2, jb, j2);
      if (!a) continue;
      if (IsOdd((jf+1+k+jb)/2)) a = -a;
      i = (1+ss)/2;
      if (nq[i] == NULL) {
	nq[i] = malloc(sizeof(double)*nk);
	dn[i] = malloc(sizeof(double)*nk);
	for (tt = 0; tt < nk; tt++) {
	  nq[i][tt] = 0.0;
	  dn[i][tt] = 0.0;
	}
      }
      a *= sqrt((k+1.0)*(j2+1.0));
      //printf("lsj: %d %d %d %d %d %d %d %d %g %g %g\n", r0->b, r0->f, jb, jf, t, ss, j2, i, a, r0->nq[t], a*r0->nq[t]);
      tt = n*(n-1)/2 + k2/2;
      nq[i][tt] += a*r0->nq[t];
      dn[i][tt] = r0->dn[t];
    }
  }
  for (ss = 0; ss < ssm; ss++) {
    if (nq[ss] == NULL) continue;
    for (tt = 0; tt < nk; tt++) {
      if (1==1+nq[ss][tt]) continue;
      r1->n++;
    }
  }
  r1->nk = malloc(sizeof(int)*r1->n);
  r1->nq = malloc(sizeof(double)*r1->n);
  r1->dn = malloc(sizeof(double)*r1->n);
  k = 0;
  for (ss = 0; ss < ssm; ss++) {
    if (nq[ss] == NULL) continue;
    for (tt = 0; tt < nk; tt++) {
      if (1==1+nq[ss][tt]) continue;
      n = nid[tt];
      kk = tt - n*(n-1)/2;
      r1->nk[k] = (2*ss+1)*10000 + n*100 + kk;
      r1->nq[k] = nq[ss][tt];
      r1->dn[k] = dn[ss][tt];
      k++;
    }
    free(nq[ss]);
    free(dn[ss]);	 
  }
  free(nq);
  free(dn);
  free(nid);
  r1->b = r0->b;
  r1->f = r0->f;  
}

void RecoupleRO(char *ifn, char *ofn) {
  F_HEADER fh0;
  TFILE *f0, *f1;
  int n, nh, i, swp, nb;
  RO_HEADER h0;
  RO_RECORD r, rs;

  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return;
  }  
  
  f0 = FOPEN(ifn, "r");
  if (f0 == NULL) return;
  n = ReadFHeader(f0, &fh0, &swp);
  if (n == 0) {
    FCLOSE(f0);
    return;
  }

  f1 = OpenFile(ofn, &fh0);
  nb = 0;
  while (1) {
    nh = ReadROHeader(f0, &h0, swp);
    if (nh == 0) break;
    InitFile(f1, &fh0, &h0);
    for (i = 0; i < h0.ntransitions; i++) {
      n = ReadRORecord(f0, &r, swp);
      if (n == 0) break;
      RecoupleRORecord(&r, &rs);
      free(r.nk);
      free(r.nq);
      if (rs.n > 0) {
	WriteRORecord(f1, &rs);
	free(rs.nk);
	free(rs.nq);
      }
    }
    DeinitFile(f1, &fh0);
    nb++;
  }
  CloseFile(f1, &fh0);
  FCLOSE(f0);
}

int CompareENRecordEnergy(const void *p0, const void *p1) {
  EN_RECORD *r0, *r1;
  
  r0 = (EN_RECORD *) p0;
  r1 = (EN_RECORD *) p1;
  if (r0->energy < r1->energy) {
    return -1;
  } else if (r0->energy > r1->energy) {
    return 1;
  } else {
    return 0;
  }
}

int CompareENRecord(const void *p0, const void *p1) {
  EN_RECORD *r0, *r1;
  
  r0 = (EN_RECORD *) p0;
  r1 = (EN_RECORD *) p1;
  
  if (r0->j < r1->j) {
    return 1;
  } else if (r0->j > r1->j) {
    return -1;
  } else {
    if (r0->p < 0 && 0 < r1->p) {
      return -1;
    } else if (r0->p > 0 && 0 > r1->p) {
      return 1;
    } else {
      if (r0->energy < r1->energy) {
	return -1;
      } else if (r0->energy > r1->energy) {
	return 1;
      } else {
	return 0;
      }
    }
  }
}

int CompareENComplex(const void *c1, const void *c2) {
  EN_RECORD *r1, *r2;

  r1 = (EN_RECORD *) c1;
  r2 = (EN_RECORD *) c2;
  return strcmp(r1->ncomplex, r2->ncomplex);
}

int SortUniqNComplex(int n, EN_RECORD *a) {
  int i, j;
  EN_RECORD b;

  qsort(a, n, sizeof(EN_RECORD), CompareENComplex);
  j = 1;
  memcpy(&b, &a[0], sizeof(EN_RECORD));
  for (i = 1; i < n; i++) {
    if (CompareENComplex(&a[i], &b) != 0) {
      if (i != j) {
	memcpy(&a[j], &a[i], sizeof(EN_RECORD));
      }
      memcpy(&b, &a[i], sizeof(EN_RECORD));
      j++;
    }
  }

  return j;
}

int FindLevelBlock(int n0, EN_RECORD *r0, int n1, EN_RECORD *r1, 
		   int nele, char *ifn) {
  F_HEADER fh;
  EN_HEADER h;
  EN_RECORD g;
  TFILE *f;
  int i, k, j, nr, nb, nv;
  int swp, sfh;
  int mk0[1024], mk1[1024];

  f = OpenFileRO(ifn, &fh, &swp);
  if (f == NULL) {
    printf("File %s does not exist\n", ifn);
    return -1;
  }

  if (VersionLE((&fh), 1, 0, 8)) sfh = sizeof(F_HEADER);
  else sfh = SIZE_F_HEADER;

  for (i = 0; i < 1024; i++) {
    mk0[i] = -1;
    mk1[i] = -1;
  }
  for (i = 0; i < n0; i++) {
    nv = abs(r0[i].p);    
    j = nv/100;
    nv = nv%100;
    if (mk0[j] < nv) mk0[j] = nv;
  }

  EN_RECORD *r0c;
  r0c = malloc(sizeof(EN_RECORD)*n0);
  memcpy(r0c, r0, sizeof(EN_RECORD)*n0);
  int n0c = SortUniqNComplex(n0, r0c);
  
  k = 0;
  for (nb = 0; nb < fh.nblocks; nb++) {
    nr = ReadENHeader(f, &h, swp);
    if (h.nele != nele) {
      FSEEK(f, h.length, SEEK_CUR);
      continue;
    }
    for (i = 0; i < h.nlevels; i++) {
      nr = ReadENRecord(f, &r1[k], swp);
      for (j = 0; j < n0c; j++) {
	if (strcmp(r1[k].ncomplex, r0c[j].ncomplex) == 0) {
	  break;
	}
      }
      if (j < n0c) {
	nv = abs(r1[k].p);
	j = nv/100;
	nv = nv%100;
	if (mk0[j] >= nv) {
	  if (mk1[j] < nv) mk1[j] = nv;	
	  k++;
	  if (k == n1) break;
	}
      }
    }
    if (k == n1) break;
  }
  FCLOSE(f);
  free(r0c);
  n1 = k;
  for (i = 0; i < 1024; i++) {
    if (mk1[i] > mk0[i]) mk1[i] = mk0[i];
  }
  int nk0 = 0;
  for (i = 0; i < n0; i++) {
    nv = abs(r0[i].p);
    j = nv/100;
    nv = nv%100;
    if (nv > mk1[j]) {
      r0[i].j = -(r0[i].j+1);
    } else {      
      nk0++;
    }
  }
  int nk1 = 0;
  for (i = 0; i < n1; i++) {
    nv = abs(r1[i].p);
    j = nv/100;
    nv = nv%100;
    if (nv > mk1[j]) {
      r1[i].j = -(r1[i].j+1);
    } else {
      nk1++;
    }
  }
  qsort(r0, n0, sizeof(EN_RECORD), CompareENRecord);
  qsort(r1, n1, sizeof(EN_RECORD), CompareENRecord);

  double eb0=0, eb1=0, w0=0,w1=0;
  for (i = 0; i < nk0; i++) {
    eb0 += (r0[i].j+1)*r0[i].energy;
    w0 += r0[i].j+1;
  }
  for (i = 0; i < nk1; i++) {
    eb1 += (r1[i].j+1)*r1[i].energy;
    w1 += r1[i].j+1;
  }
  eb0 /= w0;
  eb1 /= w1;
  int na0, nb0, na1, nb1;
  na0 = 0;
  while (na0 < nk0) {
    for (na1 = na0; na1 < nk0; na1++) {
      if (r0[na1].j != r0[na0].j || r0[na1].p*r0[na0].p < 0) break;
    }
    for (nb0 = 0; nb0 < nk1; nb0++) {
      if (r1[nb0].j == r0[na0].j && r1[nb0].p*r0[na0].p > 0) break;
    }
    for (nb1 = nb0; nb1 < nk1; nb1++) {
      if (r1[nb1].j != r1[nb0].j || r1[nb1].p*r1[nb0].p < 0) break;
    }
    int ni0 = na1-na0;
    int ni1 = nb1-nb0;
    double de0 = eb1-eb0;
    double de;
    if (ni0 < ni1) {
      j = nb0;
      for (i = 0; i < ni0; i++) {
	for (; j < nb1; j++) {
	  de = de0+r0[i+na0].energy-r1[j].energy;
	  if (0 == strcmp(r0[i+na0].sname, r1[j].sname) &&
	      0 == strcmp(r0[i+na0].name, r1[j].name) &&
	      fabs(de) < 0.2) {
	    j++;
	    break;
	  } else {
	    r1[j].j = -(r1[j].j+1);
	  }
	}
      }
      for (; j < nb1; j++) {
	r1[j].j = -(r1[j].j+1);
      }
    } else if (ni0 > ni1) {
      j = na0;
      for (i = 0; i < ni1; i++) {
	for(; j < na1; j++) {
	  de = de0+r0[j].energy-r1[i+nb0].energy;
	  if (0 == strcmp(r0[j].sname, r1[i+nb0].sname) &&
	      0 == strcmp(r0[j].name, r1[i+nb0].name) &&
	      fabs(de) < 0.2) {
	    j++;
	    break;
	  } else {
	    r0[j].j = -(r0[j].j+1);
	  }
	}
      }
      for (; j < na1; j++) {
	r0[j].j = -(r0[j].j+1);
      }
    }
    na0 = na1;
  }  
  qsort(r0, n0, sizeof(EN_RECORD), CompareENRecord);
  qsort(r1, n1, sizeof(EN_RECORD), CompareENRecord);
  nk0 = 0;
  nk1 = 0;
  for (i = 0; i < n0; i++) {
    if (r0[i].j < 0) break;
  }
  nk0 = i;
  for (i = 0; i < n1; i++) {
    if (r1[i].j < 0) break;
  }
  nk1 = i;

  if (nk0 != nk1) return -1;
  
  int n = nk0;
  EN_RECORD *r2 = (EN_RECORD *) malloc(sizeof(EN_RECORD)*n*2);
  j = 0;
  for (i = 0; i < n; i++, j+=2) {
    r2[j] = r0[i];
    r2[j+1] = r1[i];
  }
  qsort(r2, n, 2*sizeof(EN_RECORD), CompareENRecordEnergy);
  j = 0;
  for (i = 0; i < n; i++, j+=2) {
    r0[i] = r2[j];
    r1[i] = r2[j+1];
  }
  free(r2);
  for (i = n; i < n0; i++) {
    r0[i].j = -(r0[i].j+1);
  }
  for (i = n; i < n1; i++) {
    r1[i].j = -(r1[i].j+1);
  }
  return n;
}

void CombineDBase(char *pref, int k0, int k1, int nexc, int ic) {
  int k, i, n, nb, nlevs, clevs, ni0, ni1, vn, z;
  char ifn[1024], ofn[1024], buf[1024];
  char a[8];
  F_HEADER fh, fh1[6];
  EN_HEADER h0;
  EN_RECORD r0, *ri0, *ri1;
  TR_HEADER h1;
  TR_RECORD r1;
  TR_EXTRA r1x;
  CE_HEADER h2;
  CE_RECORD r2;
  RR_HEADER h3;
  RR_RECORD r3;
  CI_HEADER h4;
  CI_RECORD r4;
  AI_HEADER h5;
  AI_RECORD r5;
  DR_HEADER h6;
  DR_RECORD r6;
  char ext[6][3] = {"en", "tr", "ce", "rr", "ci", "ai"};
  int types[6] = {DB_EN, DB_TR, DB_CE, DB_RR, DB_CI, DB_AI};
  TFILE *f0, *f1[6];
  int swp, *im, *imp, **ima, nim, nk, nilevs, nplevs, nth;
  double e0, e1, e0p, e1p, de;
  int ilow2ph[3], iup2ph[3];
  double elow2ph[3], eup2ph[3];
  int ncap, nt, nd;
  double t0, dt, d0, dd;

  ncap = 0;
  nt = 0;
  nd = 0;
  d0 = 0.0;
  t0 = 0.0;
  dd = 0.0;
  dt = 0.0;
  for (k = 0; k < 3; k++) {
    ilow2ph[k] = 0;
    iup2ph[k] = 0;
  }
  nk = k1-k0+1;
  ima = malloc(sizeof(int *)*nk);

  FILE *frp0, *frp1, *frc0, *frc1;
  
  sprintf(ofn, "%s%02d%02db.rc", pref, k0, k1);  
  frc1 = fopen(ofn, "w");
  if (frc1) {
    fprintf(frc1, "#MEXC: %d\n", nexc);
  }
    
  sprintf(ofn, "%s%02d%02db.rp", pref, k0, k1);
  frp1 = fopen(ofn, "w");
  if (frp1 != NULL) {
    for (k = k1; k >= k0; k--) {
      sprintf(ifn, "%s%02db.rp", pref, k);
      frp0 = fopen(ifn, "r");      
      if (frp0 != NULL) {
	while(NULL != fgets(buf, 1024, frp0)) {
	  fprintf(frp1, "%s", buf);
	}
	fclose(frp0);
      }      
    }
    fclose(frp1);
  }
  
  nth = 0;
  z = 0;
  for (k = k1; k >= k0; k--) {
    sprintf(ifn, "%s%02db.en", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (z == 0) {
      z = (int)(fh.atom);
      strcpy(a, fh.symbol);
    }
    if (fh.nthreads > nth) nth = fh.nthreads;
    FCLOSE(f0);
  }
  
  for (i = 0; i < 6; i++) {
    sprintf(ofn, "%s%02d%02db.%s", pref, k0, k1, ext[i]);
    fh1[i].atom = z;
    strcpy(fh1[i].symbol, a);
    fh1[i].type = types[i];
    f1[i] = OpenFileWTN(ofn, &fh1[i], nth);
  }
  clevs = 0;
  nlevs = 0;
  e0 = 0.0;
  e1 = 0.0;
  de = 0.0;
  for (k = k1; k >= k0; k--) {
    sprintf(ifn, "%s%02db.rc", pref, k);
    ncap = 0;
    frc0 = fopen(ifn, "r");
    if (frc0) {
      while (NULL != fgets(buf, 1024, frc0)) {
	if (!strstr(buf, ":")) break;
	if (strstr(buf, "NCAP:")) {
	  ncap = atoi(buf+7);
	}
	if (buf[0] == '#') {
	  fprintf(frc1, "%s", buf);
	}
      }	  
      fclose(frc0);
    }
    sprintf(ifn, "%s%02db.en", pref, k);
    printf("rebuild level indices %d %d %s %d %d\n", z, k, ifn, nexc, ncap);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 == NULL) {
      printf("cannot open file %s\n", ifn);
      continue;
    }
    if (fh.type != DB_EN) {
      printf("%s is not of type DB_EN\n");
      continue;
    }
    nplevs = nlevs;
    nlevs = 0;
    e0p = e0;
    e1p = e1;
    e0 = 0.0;
    e1 = 0.0;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadENHeader(f0, &h0, swp);
      if (n == 0) break;
      for (i = 0; i < h0.nlevels; i++) {
	n = ReadENRecord(f0, &r0, swp);
	if (h0.nele == k && r0.energy < e0) e0 = r0.energy;
	if (h0.nele == k-1 && r0.energy < e1) e1 = r0.energy;
      }
      nlevs += h0.nlevels;
    }
    if (k < k1) {
      de = e1p - e0;
    }
    FSEEK(f0, SIZE_F_HEADER, SEEK_SET);
    im = malloc(sizeof(int)*nlevs);
    ima[k-k0] = im;
    nilevs = 0;
    for (i = 0; i < nlevs; i++) im[i] = -1;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadENHeader(f0, &h0, swp);
      if (n == 0) break;
      if (h0.nele == k) {	
	InitFile(f1[0], &fh1[0], &h0);
	for (i = 0; i < h0.nlevels; i++) {
	  n = ReadENRecord(f0, &r0, swp);
	  vn = abs(r0.p)/100;
	  if ((nexc > 0 && vn > nexc) ||
	      (ncap > 0 && vn > ncap && r0.energy > e1)) {
	    continue;
	  }
	  im[r0.ilev] = clevs;
	  if (r0.ibase >= 0) {
	    r0.ibase = im[r0.ibase];
	  }
	  r0.ilev = im[r0.ilev];
	  r0.energy += de;
	  Match2PhotonLevels(k, &r0, ilow2ph, iup2ph, elow2ph, eup2ph);
	  WriteENRecord(f1[0], &r0);
	  clevs++;
	}
	DeinitFile(f1[0], &fh1[0]);
      } else {
	if (h0.nele == k-1 && k > k0) {
	  ni0 = h0.nlevels;
	  ni1 = ni0*5;
	  ri0 = malloc(sizeof(EN_RECORD)*ni0);
	  ri1 = malloc(sizeof(EN_RECORD)*ni1);
	  for (i = 0; i < h0.nlevels; i++) {
	    n = ReadENRecord(f0, &ri0[i], swp);
	    if (n == 0) break;
	  }
	  sprintf(ifn, "%s%02db.en", pref, k-1);
	  nim = FindLevelBlock(ni0, ri0, ni1, ri1, k-1, ifn);
	  for (i = 0; i < nim; i++) {
	    im[ri0[i].ilev] = -10-ri1[i].ilev;
	  }
	  free(ri0);
	  free(ri1);
	  nilevs += ni0;
	} else {
	  nilevs += h0.nlevels;
	  FSEEK(f0, h0.length, SEEK_CUR);
	}
      }
    }
    if (k < k1) {
      imp = ima[k+1-k0];
      for (i = 0; i < nplevs; i++) {
	if (imp[i] <= -10) {
	  imp[i] = im[-(10+imp[i])];
	}
      }
    }
    //clevs += nlevs-nilevs;
    if (k == k0) {
      FSEEK(f0, SIZE_F_HEADER, SEEK_SET);
      nilevs = 0;
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadENHeader(f0, &h0, swp);
	if (n == 0) break;
	if (h0.nele == k-1) {	
	  InitFile(f1[0], &fh1[0], &h0);
	  for (i = 0; i < h0.nlevels; i++) {
	    n = ReadENRecord(f0, &r0, swp);
	    vn = abs(r0.p)/100;
	    if ((nexc > 0 && vn > nexc) ||
		(ncap > 0 && vn > ncap && r0.energy > e1)) {
	      continue;
	    }
	    im[r0.ilev] = clevs;
	    if (r0.ibase >= 0) {
	      r0.ibase = im[r0.ibase];
	    }
	    r0.ilev = im[r0.ilev];
	    r0.energy += de;
	    WriteENRecord(f1[0], &r0);
	    clevs++;
	  }
	  DeinitFile(f1[0], &fh1[0]);
	} else {
	  nilevs += h0.nlevels;
	  FSEEK(f0, h0.length, SEEK_CUR);
	}
      }
    }
    FCLOSE(f0);
  }
  
  CloseFile(f1[0], &fh1[0]);
  sprintf(ifn, "%s%02d%02db.en", pref, k0, k1);
  MemENTable(ifn);

  if (frc1) {
    int nid, id0, id1, id2, id3, id4, id5, id6, ndr[128], nre[128];
    fprintf(frc1, "#STARTDATA\n");
    for (k = k1; k >= k0; k--) {
      im = ima[k-k0];
      sprintf(ifn, "%s%02db.rc", pref, k);
      frc0 = fopen(ifn, "r");
      ndr[k] = 0;
      if (frc0) {
	while(NULL != fgets(buf, 1024, frc0)) {
	  if (buf[0] == '#') continue;
	  nid = sscanf(buf, "%d %d %d %d %d %d %d",
		       &id0, &id1, &id2, &id3, &id4, &id5, &id6);
	  if (nid != 7) continue;
	  if (id0 != 6 && id0 != 4) continue;
	  id3 = im[id3];
	  id5 = im[id5];
	  if (id3 >= 0 && id5 >= 0) {
	    fprintf(frc1, "%d %3d %2d %8d %3d %8d %3d %s",
		    id0, id1, id2, id3, id4, id5, id6, buf+35);
	    if (id0 == 6) ndr[k]++;
	  }
	}
	fclose(frc0);
      }
    }
    for (k = k1; k >= k0; k--) {
      im = ima[k-k0];
      sprintf(ifn, "%s%02db.rc", pref, k);
      frc0 = fopen(ifn, "r");
      nre[k] = 0;
      if (frc0) {
	while(NULL != fgets(buf, 1024, frc0)) {
	  if (buf[0] == '#') continue;
	  nid = sscanf(buf, "%d %d %d %d %d %d %d",
		       &id0, &id1, &id2, &id3, &id4, &id5, &id6);
	  if (nid != 7) continue;
	  if (id0 != 5) continue;
	  id3 = im[id3];
	  id5 = im[id5];
	  if (id3 >= 0 && id5 >= 0) {
	    fprintf(frc1, "%d %3d %2d %8d %3d %8d %3d %s",
		    id0, id1, id2, id3, id4, id5, id6, buf+35);
	    nre[k]++;
	  }
	}
	fclose(frc0);
      }
    }
    fclose(frc1);
    sprintf(ofn, "%s%02d%02db.rc", pref, k0, k1);  
    frc1 = fopen(ofn, "r+");
    fseek(frc1, 0, SEEK_SET);
    k = k1;
    while(NULL != fgets(buf, 1024, frc1)) {
      if (strstr(buf, "#RDG:")) {
	printf("ndr/re: %d %d %d\n", k, ndr[k], nre[k]);
	fseek(frc1, 0, SEEK_CUR);
	sprintf(buf, "#NCE: %d %3d %8d\n", RC_CE, k, 0);
	fwrite(buf, 1, strlen(buf), frc1);
	sprintf(buf, "#NCI: %d %3d %8d\n", RC_CI, k, 0);
	fwrite(buf, 1, strlen(buf), frc1);
	sprintf(buf, "#NRR: %d %3d %8d\n", RC_RR, k, 0);
	fwrite(buf, 1, strlen(buf), frc1);
	sprintf(buf, "#NDR: %d %3d %8d\n", RC_DR, k, ndr[k]);
	fwrite(buf, 1, strlen(buf), frc1);
	sprintf(buf, "#NRE: %d %3d %8d\n", RC_RE, k, nre[k]);
	fwrite(buf, 1, strlen(buf), frc1);
	sprintf(buf, "#NEA: %d %3d %8d\n", RC_EA, k, ndr[k]);
	fwrite(buf, 1, strlen(buf), frc1);
	k--;
      } else if (strstr(buf, "#STARTDATA")) {
	break;
      }
    }
    fclose(frc1);
  }
  
  for (k = k1; k >= k0; k--) {
    printf("processing %d %d\n", z, k);
    im = ima[k-k0];
    sprintf(ifn, "%s%02db.tr", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_TR) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadTRHeader(f0, &h1, swp);
	if (n == 0) break;
	InitFile(f1[1], &fh1[1], &h1);
	for (i = 0; i < h1.ntransitions; i++) {
	  n = ReadTRRecord(f0, &r1, &r1x, swp);
	  if (n == 0) break;
	  if (im[r1.lower] >= 0 && im[r1.upper] >= 0) {
	    r1.lower = im[r1.lower];
	    r1.upper = im[r1.upper];
	    WriteTRRecord(f1[1], &r1, &r1x);
	  }
	}
	DeinitFile(f1[1], &fh1[1]);
      }
      double r2p = 0.0;
      i = -1;
      if (k == 1 && ilow2ph[0] && iup2ph[0]) {
	i = 0;
      } else if (k == 2 && ilow2ph[1] && iup2ph[1]) {
	i = 1;
      } else if (k == 4 && ilow2ph[2] && iup2ph[2]) {
	i = 2;
      }
      if (i >= 0) {
	r2p = TwoPhotonRate(z, i);
	if (i == 0) r2p *= 2;
	de = eup2ph[i]-elow2ph[i];	
	de = pow(de*FINE_STRUCTURE_CONST,2);
	r2p /= 2*de*FINE_STRUCTURE_CONST*RATE_AU;
	h1.ntransitions = 1;
	h1.gauge = G_TWOPHOTON;
	h1.multipole = 0;
	InitFile(f1[1], &fh1[1], &h1);
	r1.lower = ilow2ph[i];
	r1.upper = iup2ph[i];
	r1.strength = (float)r2p;
	WriteTRRecord(f1[1], &r1, &r1x);
      }
      DeinitFile(f1[1], &fh1[1]);
    }
    FCLOSE(f0);
    
    sprintf(ifn, "%s%02db.ce", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_CE) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadCEHeader(f0, &h2, swp);
	if (n == 0) break;
	InitFile(f1[2], &fh1[2], &h2);
	for (i = 0; i < h2.ntransitions; i++) {
	  n = ReadCERecord(f0, &r2, swp, &h2);
	  if (n == 0) break;
	  if (im[r2.lower] >= 0 && im[r2.upper] >= 0) {
	    r2.lower = im[r2.lower];
	    r2.upper = im[r2.upper];	    
	    WriteCERecord(f1[2], &r2);
	  }
	  if (h2.qk_mode == QK_FIT) free(r2.params);
	  free(r2.strength);
	}	
	DeinitFile(f1[2], &fh1[2]);
	free(h2.tegrid);
	free(h2.egrid);
	free(h2.usr_egrid);
      }
    }
    FCLOSE(f0);
    
    sprintf(ifn, "%s%02db.rr", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_RR) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadRRHeader(f0, &h3, swp);
	if (n == 0) break;
	InitFile(f1[3], &fh1[3], &h3);
	for (i = 0; i < h3.ntransitions; i++) {
	  n = ReadRRRecord(f0, &r3, swp, &h3);
	  if (n == 0) break;
	  if (im[r3.b] >= 0 && im[r3.f] >= 0) {
	    r3.b = im[r3.b];
	    r3.f = im[r3.f];
	    de = mem_en_table[r3.f].energy - mem_en_table[r3.b].energy;
	    if (de > 0) {
	      WriteRRRecord(f1[3], &r3);
	    }
	  }
	  free(r3.params);
	  free(r3.strength);
	}
	DeinitFile(f1[3], &fh1[3]);
	free(h3.tegrid);
	free(h3.egrid);
	free(h3.usr_egrid);
      }
    }
    FCLOSE(f0);
    
    sprintf(ifn, "%s%02db.ci", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_CI) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadCIHeader(f0, &h4, swp);
	if (n == 0) break;
	InitFile(f1[4], &fh1[4], &h4);
	for (i = 0; i < h4.ntransitions; i++) {
	  n = ReadCIRecord(f0, &r4, swp, &h4);
	  if (n == 0) break;
	  if (im[r4.b] >= 0 && im[r4.f] >= 0) {
	    r4.b = im[r4.b];
	    r4.f = im[r4.f];
	    de = mem_en_table[r4.f].energy-mem_en_table[r4.b].energy;
	    if (de > 0) {
	      WriteCIRecord(f1[4], &r4);
	    }
	  }
	  free(r4.params);
	  free(r4.strength);
	}
	DeinitFile(f1[4], &fh1[4]);
	free(h4.tegrid);
	free(h4.egrid);
	free(h4.usr_egrid);
      }
    }
    FCLOSE(f0);
      
    sprintf(ifn, "%s%02db.ai", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_AI) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadAIHeader(f0, &h5, swp);
	if (n == 0) break;
	InitFile(f1[5], &fh1[5], &h5);
	for (i = 0; i < h5.ntransitions; i++) {
	  n = ReadAIRecord(f0, &r5, swp);
	  if (n == 0) break;
	  if (im[r5.b] >= 0 && im[r5.f] >= 0) {
	    r5.b = im[r5.b];
	    r5.f = im[r5.f];
	    de = mem_en_table[r5.b].energy - mem_en_table[r5.f].energy;
	    if (de > 0) {
	      WriteAIRecord(f1[5], &r5);
	    }
	  }
	}
	DeinitFile(f1[5], &fh1[5]);
	free(h5.egrid);
      }
    }
    FCLOSE(f0);    
    free(im);
  }
  free(ima);
  for (i = 1; i < 6; i++) {
    CloseFile(f1[i], &fh1[i]);
  }
  
  if (ic >= 0) {
    for (i = 0; i < 6; i++) {
      sprintf(ifn, "%s%02d%02db.%s", pref, k0, k1, ext[i]);
      sprintf(ofn, "%s%02d%02da.%s", pref, k0, k1, ext[i]);
      PrintTable(ifn, ofn, 1);
    }
  }
}
  
/*
** two-photon rates are taken from G. W. F. Drake, PRA, 34, 2871, 1986.
*/
double TwoPhotonRate(double z, int t) {
  double a, a2, a4, z6;
  
  switch (t) {
  case 0: /* 2S_1/2 of H-like ion */
    z6 = z*z;
    z6 = z6*z6*z6;
    a = FINE_STRUCTURE_CONST*z;
    a2 = a*a;
    a4 = a2*a2;
    a = 8.22943*z6*(1.0 + 3.9448*a2 - 2.04*a4)/(1.0 + 4.6019*a2);
    break;
  case 1: /* 1s2s S_0 of He-like ion */
    a = (z - 0.806389);
    z6 = a*a;
    z6 = z6*z6*z6;
    a = FINE_STRUCTURE_CONST*a;
    a2 = a*a;
    a4 = (z+2.5);
    a4 = a4*a4;
    a = 16.458762*(z6*(1.0 + 1.539/a4) - 
		   z6*a2*(0.6571 + 2.04*a2)/(1.0 + 4.6019*a2));
    break;
  case 2: /* 2s2p J=0 of Be-like ion */
    a = log(z);
    UVIP3P(1, 8, _tpbez, _tpber, 1, &a, &a2);
    a = exp(a2);
    break;
  default:
    a = 0.0;
    break;
  }

  return a;
}

void ClearIdxMap(void) {
  ArrayFree(_idxmap.imap, NULL);
  if (_idxmap.nij > 0) free(_idxmap.mask);
  _idxmap.mask = NULL;
  _idxmap.ni = 0;
  _idxmap.nj = 0;
  _idxmap.i0 = 0;
  _idxmap.im0 = 0;
  _idxmap.j0 = 0;
  _idxmap.jm0 = 0;
  _idxmap.i1 = 0;
  _idxmap.im1 = 0;
  _idxmap.j1 = 0;
  _idxmap.jm1 = 0;
  int i;
  for (i = 0; i < 32; i++) {
    _idxmap.xm[i] = 1<<i;
  }
}

int PreloadEN(char *fn, int i0, int i1, int j0, int j1) {
  FILE *f;
  char buf[2048];
  int i, im, n;
  double e, em;
  IDXDAT d, *ip;
  
  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  ClearIdxMap();
  _idxmap.i0 = i0;
  _idxmap.j0 = j0;
  _idxmap.im0 = 1000000000;
  _idxmap.jm0 = 1000000000;
  _idxmap.i1 = 0;
  _idxmap.j1 = 0;
  _idxmap.im1 = 0;
  _idxmap.jm1 = 0;
  while (1) {
    if (NULL == fgets(buf, 2048, f)) break;
    n = sscanf(buf, "%d %d %lf %lf", &im, &i, &em, &e);
    if (n != 4) continue;
    if (e <= 0 && im > 0) continue;
    e /= HARTREE_EV;
    d.i = im;
    if (im == 0) {
      d.e = 0.0;
    } else {
      d.e = e;
    }
    if (i > 0 || im == 0) {
      ArraySet(_idxmap.imap, i, &d, NULL);
      if (_idxmap.i1 < i) _idxmap.i1 = i;
    }
    AddECorrection(0, im, e, 1);
  }
  fclose(f);
  _idxmap.j1 = _idxmap.i1;
  if (i1 > 0 && i1 < _idxmap.i1) _idxmap.i1 = i1;
  if (j1 > 0 && j1 < _idxmap.j1) _idxmap.j1 = j1;
  for (i = 0; i < _idxmap.imap->dim; i++) {
    ip = ArrayGet(_idxmap.imap, i);
    if (!ip) continue;
    if (i >= _idxmap.i0 && i <= _idxmap.i1) {
      if (_idxmap.im0 > i) _idxmap.im0 = i;
      if (_idxmap.im1 < i) _idxmap.im1 = i;
    }
    if (i >= _idxmap.j0 && i <= _idxmap.j1) {
      if (_idxmap.jm0 > i) _idxmap.jm0 = i;
      if (_idxmap.jm1 < i) _idxmap.jm1 = i;
    }
  }
  _idxmap.ni = 1 + (_idxmap.im1 - _idxmap.im0);
  _idxmap.nj = 1 + (_idxmap.jm1 - _idxmap.jm0);
  _idxmap.nij = _idxmap.ni*_idxmap.nj;
  if (_idxmap.nij > 0) {
    _idxmap.mask = malloc(sizeof(long)*_idxmap.nij);
    for (i = 0; i < _idxmap.nij; i++) _idxmap.mask[i] = 0;
  }
  return 0;
}

int SetPreloaded(int i, int j, int m) {
  if (_idxmap.nij == 0) return 0;
  if (i < _idxmap.im0) return 0;
  if (i > _idxmap.im1) return 0;
  if (j < _idxmap.jm0) return 0;
  if (j > _idxmap.jm1) return 0;
  i -= _idxmap.im0;
  j -= _idxmap.jm0;
  _idxmap.mask[j*_idxmap.ni+i] |= _idxmap.xm[m];
  return 1;
}

int SetPreloadedTR(int i, int j, int m) {
  if (m >= 8 || m <= -8) return 0;
  if (m < 0) {
    return SetPreloaded(i, j, (-m)-1);
  }
  if (m > 0) {
    return SetPreloaded(i, j, m+7);
  }
  return SetPreloaded(i, j, 16);
}

int SetPreloadedCE(int i, int j) {
  return SetPreloaded(i, j, 17);
}

int IsPreloaded(int i, int j, int m) {
  if (_idxmap.nij == 0) return 0;
  if (i < _idxmap.im0) return 0;
  if (i > _idxmap.im1) return 0;
  if (j < _idxmap.jm0) return 0;
  if (j > _idxmap.jm1) return 0;
  i -= _idxmap.im0;
  j -= _idxmap.jm0;
  return _idxmap.mask[j*_idxmap.ni+i]&_idxmap.xm[m];
}

int IsPreloadedTR(int i, int j, int m) {
  if (m >= 8 || m <= -8) return 0;
  if (m < 0) {
    return IsPreloaded(i, j, (-m)-1);
  }
  if (m > 0) {
    return IsPreloaded(i, j, m+7);
  }
  return IsPreloaded(i, j, 16);
}

int IsPreloadedCE(int i, int j) {
  return IsPreloaded(i, j, 17);
}

IDXDAT *IdxMap(int i) {
  return ArrayGet(_idxmap.imap, i);
}

int PreloadTable(char *tfn, char *sfn, int m) {
  TFILE *f0;
  F_HEADER fh;
  int swp;

  f0 = OpenFileRO(sfn, &fh, &swp);
  if (f0 == NULL) {
    printf("cannot open file %s\n", sfn);
    return -1;
  }
  FCLOSE(f0);
  switch(fh.type) {
  case DB_TR:
    PreloadTR(tfn, sfn, m);
    break;
  case DB_CE:
    PreloadCE(tfn, sfn);
    break;
  default:
    break;
  }
  return 0;
}

int PreloadTR(char *tfn, char *sfn, int m) {
  TFILE *f0, *f1;
  int ib, i, j, k, t, n, swp;
  float *rt;
  F_HEADER fh, fh1;
  TR_HEADER h;
  TR_RECORD r;
  TR_EXTRA rx;
  IDXDAT *di, *dj;

  if (_idxmap.nij == 0) {
    printf("index map not loaded\n");
    return -1;
  }
  f0 = OpenFileRO(sfn, &fh, &swp);
  if (f0 == NULL) {
    printf("cannot open file %s\n", sfn);
    return -1;
  }
  if (fh.type != DB_TR) {
    printf("file %s is not of type DB_TR\n", sfn);
    return -1;
  }
  fh1.atom = fh.atom;
  strcpy(fh1.symbol, fh.symbol);
  fh1.type = DB_TR;
  f1 = OpenFile(tfn, &fh1);
  if (f1 == NULL) {
    printf("cannot open file %s\n", tfn);
    FCLOSE(f0);
    return -1;
  }
  if (m == 0) {
    rt = malloc(sizeof(float)*_idxmap.nij);
    for (i = 0; i < _idxmap.nij; i++) {
      rt[i] = 0.0;
    }
  }
  
  for (ib = 0; ib < fh.nblocks; ib++) {
    n = ReadTRHeader(f0, &h, swp);
    if (n == 0) break;
    if (m) {
      InitFile(f1, &fh1, &h);
    }
    for (k = 0; k < h.ntransitions; k++) {
      n = ReadTRRecord(f0, &r, &rx, swp);
      di = IdxMap(r.lower);
      dj = IdxMap(r.upper);
      if (!di || !dj) continue;
      r.lower = di->i;
      r.upper = dj->i;
      if (SetPreloadedTR(r.lower, r.upper, h.multipole)) {
	if (m == 0) {
	  t = (r.upper-_idxmap.jm0)*_idxmap.ni + r.lower-_idxmap.im0;
	  rt[t] += OscillatorStrength(h.multipole,
				      dj->e-di->e, r.strength, NULL);
	  SetPreloadedTR(r.lower, r.upper, 0);
	} else {
	  WriteTRRecord(f1, &r, &rx);
	}
      }
    }
    if (m) {
      DeinitFile(f1, &fh1);
    }
  }
  FCLOSE(f0);
  if (m == 0) {
    h.multipole = 0;
    InitFile(f1, &fh1, &h);
    for (i = _idxmap.im0; i <= _idxmap.im1; i++) {
      for (j = _idxmap.jm0; j <= _idxmap.jm1; j++) {
	if (IsPreloadedTR(i, j, 0)) {
	  t = (j-_idxmap.jm0)*_idxmap.ni + i-_idxmap.im0;	
	  r.lower = i;
	  r.upper = j;
	  r.strength = rt[t];
	  WriteTRRecord(f1, &r, &rx);
	}
      }
    }
    DeinitFile(f1, &fh1);
    free(rt);
  }
  CloseFile(f1, &fh1);
  return 0;
}

int PreloadCE(char *tfn, char *sfn) {
  TFILE *f0, *f1;
  int ib, k, n, swp;
  F_HEADER fh, fh1;
  CE_HEADER h;
  CE_RECORD r;
  IDXDAT *di, *dj;
  
  if (_idxmap.nij == 0) {
    printf("index map not loaded\n");
    return -1;
  }
  f0 = OpenFileRO(sfn, &fh, &swp);
  if (f0 == NULL) {
    printf("cannot open file %s\n", sfn);
    return -1;
  }
  if (fh.type != DB_CE) {
    printf("file %s is not of type DB_CE\n", sfn);
    return -1;
  }
  fh1.atom = fh.atom;
  strcpy(fh1.symbol, fh.symbol);
  fh1.type = DB_CE;
  f1 = OpenFile(tfn, &fh1);
  if (f1 == NULL) {
    printf("cannot open file %s\n", tfn);
    FCLOSE(f0);
    return -1;
  }

  for (ib = 0; ib < fh.nblocks; ib++) {
    n = ReadCEHeader(f0, &h, swp);
    if (n == 0) break;
    InitFile(f1, &fh1, &h);
    for (k = 0; k < h.ntransitions; k++) {
      n = ReadCERecord(f0, &r, swp, &h);
      if (n == 0) break;
      di = IdxMap(r.lower);
      dj = IdxMap(r.upper);
      if (!di || !dj) continue;
      r.lower = di->i;
      r.upper = dj->i;
      if (SetPreloadedCE(r.lower, r.upper)) {
	WriteCERecord(f1, &r);
      }
      if (h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }
    DeinitFile(f1, &fh1);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
  }
  FCLOSE(f0);
  CloseFile(f1, &fh1);
  return 0;
}

void SetOptionDBase(char *s, char *sp, int ip, double dp) {
  
}

