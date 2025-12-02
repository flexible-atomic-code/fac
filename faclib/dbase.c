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
#include "excitation.h"
#include "recombination.h"
#include "ionization.h"
#include "cf77.h"
#include "global.h"

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
static RC_HEADER rc_header;

static EN_SRECORD *mem_en_table = NULL;
static int mem_en_table_size = 0;
static EN_SRECORD *mem_enf_table = NULL;
static int mem_enf_table_size = 0;
static int iground;
static double _eground[N_ELEMENTS1];
static int iuta = 0;
static int utaci = 1;
static int cuta = 0;
static int tuta = 0;
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
static double _cmpetol = 0.25;
static double _cmpnbm = -1;
static int _remove_closed = 0;
static int _ncombex[N_ELEMENTS1];
static char **_pcombex[N_ELEMENTS1];
static char _scombex[N_ELEMENTS1][2048];
static int _adj_ip = 1;
static int _exk2h = 0;

static int _ic_iai[100];
static int _ix_iai[100];
static int _nc_iai = 0;
static int _nilast = 1;
static double _grptol = 0.05;
static int _collapse_mask = 0xff;
static double _sfu_smin = -1.0;
static double _sfu_dmax = 1.0;
static double _sfu_minde = -1e30;
static double _sfu_maxde = 1e30;
static double _cwf = 1.0;
static char _transmod_file[1024] = "";
static double _ste = -1.0;
static char _groupmod_file[1024] = "";

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

  q = realloc(p, s);
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
  if (m >= 0) {
    iuta = Max(iuta, m);
    cuta = m;
  }
  if (mci >= 0) utaci = mci;
}

int TransUTA(void) {
  return tuta;
}

int IsUTA(void) {
  return iuta > 0;
}

int TrueUTA(int n) {
  if (cuta <= 0) return 0;
  if (n > 0 && n < cuta/10) return 0;
  return 1;
}

int CurrentUTA(int *iu, int *ici) {
  *iu = iuta;
  *ici = utaci;
  return cuta;
}

int IsUTARecord(EN_RECORD *r) {
  int i;

  for (i = 0; i < 10; i++) {
    if (r->name[i] == '.') return 1;
    if (r->name[i] == '\0') return 1;
    if (r->name[i] == '(') return 0;
  }
  
  return 1;
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

int SwapEndianTRRecord(TR_RECORD *r, TR_EXTRA *rx, int utr) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->strength), sizeof(float));
  if (utr && rx) {
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
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  return 0;
}

int SwapEndianRCHeader(RC_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->type), sizeof(int));
  SwapEndian((char *) &(h->nexc), sizeof(int));
  SwapEndian((char *) &(h->mexc), sizeof(int));
  SwapEndian((char *) &(h->ncap), sizeof(int));
  SwapEndian((char *) &(h->nte), sizeof(int));
  SwapEndian((char *) &(h->nde), sizeof(int));
  SwapEndian((char *) &(h->te0), sizeof(double));
  SwapEndian((char *) &(h->dte), sizeof(double));
  SwapEndian((char *) &(h->de0), sizeof(double));
  SwapEndian((char *) &(h->dde), sizeof(double));
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

int SwapEndianRCRecord(RC_RECORD *r) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
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

int SwapEndianSPRecord(SP_RECORD *r, SP_EXTRA *rx, int utr) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->energy), sizeof(float));
  SwapEndian((char *) &(r->strength), sizeof(float));
  SwapEndian((char *) &(r->rrate), sizeof(float));
  SwapEndian((char *) &(r->trate), sizeof(float));
  if (utr && rx) {
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
  printf("WallTime%d: %11.5E %15.9E %9.3E/%9.3E %9.3E/%9.3E %8lld %8lld ... %s\n",
	 m, t, TotalSize(), ttskip, mtskip, ttlock, mtlock, tnlock, mnlock, s);
  fflush(stdout);
}

int *InitTransReport(int *nproc) {
  int k, np;
  k = MPIRank(&np);
  int *n = malloc(sizeof(int)*(np));
  for (k = 0; k < np; k++) n[k] = 0;
  *nproc = np;
  return n;
}

void PrintTransReport(int nproc, double t0, int *ntrans,
		      char *sid, int ip) {
  int k;
  double nt;
  int n0, n1, md;
  if (nproc > 0) {
    n1 = n0 = ntrans[0];
    nt = 0;
    for (k = 0; k < nproc; k++) {
      nt += ntrans[k];
      if (n0 > ntrans[k]) n0 = ntrans[k];
      if (n1 < ntrans[k]) n1 = ntrans[k];
    }
    nt /= nproc;
    md = 0;
  } else {
    n0 = n1 = nt = ntrans[-nproc];
    md = -1;
  }
  double t1 = WallTime();
  double dt = t1 - t0;
  double mdt = 1e3*dt;
  double mdta=0.0, mdt0=0.0, mdt1=0.0;
  if (nt > 0) mdta = mdt/nt;
  if (n1 > 0) mdt1 = mdt/n1;
  if (n0 > 0) mdt0 = mdt/n0;
  MPrintf(md, "%s %06d: %09d(%09d,%09d)trans in %8.2Es, %8.2E(%8.2E,%8.2E)ms/tran @ %11.4Es m=%11.4E\n",
	  sid, ip, ((int)(nt+0.25)), n0, n1, dt, mdta, mdt0, mdt1, ClockNow(0), TotalSize());
  if (ip < 0) free(ntrans);
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
  for (i = 0; i < N_ELEMENTS1; i++) _eground[i] = 0.0;
  itrf = 0;

  if (_idxmap.imap == NULL) {
    _idxmap.imap = malloc(sizeof(IDXMAP));
    ArrayInit(_idxmap.imap, sizeof(IDXDAT), 5000);
  }
  ClearIdxMap();
  SetCombEx(NULL);
  SetInnerAI(NULL);
  
  return 0;
}

void SetInnerAI(char *s) {
  int i;
  char *p, buf[2048];

  while (s && (*s == ' ' || *s == '\t')) s++;
  if (s == NULL || strlen(s) == 0) {
    _nc_iai = 0;
    for (i = 0; i < 100; i++) {
      _ic_iai[i] = 0;
      _ix_iai[i] = -1;
    }
    return;
  }
  strncpy(buf, s, 2048);
  _nc_iai = StrSplit(buf, ',');
  p = buf;
  for (i = 0; i < _nc_iai; i++) {
    _ic_iai[i] = atoi(p);
    _ix_iai[_ic_iai[i]] = i;
    while (*p) p++;
    p++;
  }
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
    for (i = 0; i < N_ELEMENTS1; i++) _eground[i] = 0.0;
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
  //else iuta =0;
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
  if (swp) SwapEndianTRRecord(r, rx, 0);
  
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
  //else iuta = 0;
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
  if (swp) SwapEndianSPRecord(r, rx, 0);
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

  if (iuta) {
    v |= 1<<8;
    if (utaci) v |= 1<<9;
  }
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

  iuta = (fh->version&0x100)>>8;
  utaci = (fh->version&0x200)>>9;
  fh->nthreads = (fh->version&0xFFFF0000)>>16;
  fh->version &= 0xFF;
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

int WriteRCHeader(TFILE *f, RC_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->type);
  WSF0(h->nexc);
  WSF0(h->mexc);
  WSF0(h->ncap);
  WSF0(h->nte);
  WSF0(h->nde);
  WSF0(h->te0);
  WSF0(h->dte);
  WSF0(h->de0);
  WSF0(h->dde);

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
  if (isnan(r->energy) || isinf(r->energy)) {
    MPrintf(-1, "WriteENRecord invalid energy: %d %g\n", r->ilev, r->energy);
  }
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

  if (isnan(r->energy) || isinf(r->energy)) {
    MPrintf(-1, "WriteENFRecord invalid energy: %d %d %g\n",
	    r->ilev, r->pbasis, r->energy);
  }
  WSF0(r->ilev);
  WSF0(r->energy);
  WSF0(r->pbasis);
  
#pragma omp atomic
  enf_header.length += m;

  return m;
}

int WriteTRRecord(TFILE *f, TR_RECORD *r, TR_EXTRA *rx, int utr) {
  int n, m = 0;
  
  if (1.0 == 1.0 + (double)(r->strength)) {
    return 0;
  }
  
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
  if (!utr) {
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

  if (isnan(r->strength) || isinf(r->strength)) {
    MPrintf(-1, "WriteTRRecord invalid strength: %d %d %g\n",
	    r->lower, r->upper, r->strength);
  }

  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->strength);

  if (utr) {
    if (isnan(rx->sdev) || isinf(rx->sdev)) {
      MPrintf(-1, "WriteTRRecord invlid sdev: %d %d %g\n",
	      r->lower, r->upper, rx->sdev);
    }
    WSF0(rx->energy);
    WSF0(rx->sdev);
    WSF0(rx->sci);
  }

#pragma omp atomic
  tr_header.length += m;

  return m;
}

int WriteTRFRecord(TFILE *f, TRF_RECORD *r) {
  int n, i, m = 0;

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
  for (i = 0; i < 2*abs(trf_header.multipole)+1; i++) {
    if (isnan(r->strength[i]) || isinf(r->strength[i])) {
      MPrintf(-1, "WriteTRFRecord invalid strength: %d %d %d %g\n",
	      r->lower, r->upper, i, r->strength[i]);
    }
  }
  WSF0(r->lower);
  WSF0(r->upper);
  WSF1(r->strength, sizeof(float), 2*abs(trf_header.multipole)+1);

#pragma omp atomic
  trf_header.length += m;

  return m;
}

int WriteCERecord(TFILE *f, CE_RECORD *r) {
  int n, i;
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
  if (isnan(r->bethe) || isinf(r->bethe)) {
    MPrintf(-1, "WriteCERecord invalid bethe: %d %d %g\n",
	    r->lower, r->upper, r->bethe);
  }
  WSF0(r->bethe);
  if (isnan(r->born[0]) || isinf(r->born[0])) {
    MPrintf(-1, "WriteCERecord invalid born0: %d %d %g\n",
	    r->lower, r->upper, r->born[0]);
  }
  if (isnan(r->born[1]) || isinf(r->born[1])) {
    MPrintf(-1, "WriteCERecord invalid born1: %d %d %g\n",
	    r->lower, r->upper, r->born[1]);
  }    
  WSF1(r->born, sizeof(float), 2);

  if (m0) {
    for (i = 0; i < m0; i++) {
      if (isnan(r->params[i]) || isinf(r->params[i])) {
	MPrintf(-1, "WriteCERecord invalid param: %d %d %d %g\n",
		r->lower, r->upper, i, r->params[i]);
      }
    }    
    WSF1(r->params, sizeof(float), m0);
  }
  for (i = 0; i < m1; i++) {
    if (isnan(r->strength[i]) || isinf(r->strength[i])) {
      MPrintf(-1, "WriteCERecord invalid strength: %d %d %d %g\n",
	      r->lower, r->upper, i, r->strength[i]);
    }
  }
  WSF1(r->strength, sizeof(float), m1);

#pragma omp atomic
  ce_header.length += m;

  return m;
}

int WriteCEFRecord(TFILE *f, CEF_RECORD *r) {
  int n, i;
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
  if (isnan(r->bethe) || isinf(r->bethe)) {
    MPrintf(-1, "WriteCEFRecord invalid bethe: %d %d %g\n",
	    r->lower, r->upper, r->bethe);
  }
  WSF0(r->bethe);
  if (isnan(r->born[0]) || isinf(r->born[0])) {
    MPrintf(-1, "WriteCEFRecord invalid born0: %d %d %g\n",
	    r->lower, r->upper, r->born[0]);
  }
  if (isnan(r->born[1]) || isinf(r->born[1])) {
    MPrintf(-1, "WriteCEFRecord invalid born1: %d %d %g\n",
	    r->lower, r->upper, r->born[1]);
  }
  WSF1(r->born, sizeof(float), 2);
  for (i = 0; i < m0; i++) {
    if (isnan(r->strength[i]) || isinf(r->strength[i])) {
      MPrintf(-1, "WriteCEFRecord invalid strength: %d %d %d %g\n",
	      r->lower, r->upper, i, r->strength[i]);
    }
  }
  WSF1(r->strength, sizeof(float), m0);
  
#pragma omp atomic
  cef_header.length += m;

  return m;
}

int WriteCEMFRecord(TFILE *f, CEMF_RECORD *r) {
  int n, i;
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
  
  for (i = 0; i < m0; i++) {
    if (isnan(r->bethe[i]) || isinf(r->bethe[i])) {
      MPrintf(-1, "WriteCEMFRecord invalid bethe: %d %d %d %g\n",
	      r->lower, r->upper, i, r->bethe[i]);
    }
  }
  WSF1(r->bethe, sizeof(float), m0);
  for (i = 0; i <= m0; i++) {
    if (isnan(r->born[i]) || isinf(r->born[i])) {
      MPrintf(-1, "WriteCEMFRecord invalid strength: %d %d %d %g\n",
	      r->lower, r->upper, i, r->born[i]);
    }
  }
  WSF1(r->born, sizeof(float), m0+1);
  for (i = 0; i < m1; i++) {
    if (isnan(r->strength[i]) || isinf(r->strength[i])) {
      MPrintf(-1, "WriteCEFRecord invalid strength: %d %d %d %g\n",
	      r->lower, r->upper, i, r->strength[i]);
    }
  }
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
  int n, i;
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
  
  for (i = 0; i < cx_header.ne0; i++) {
    if (isnan(r->cx[i]) || isinf(r->cx[i])) {
      MPrintf(-1, "WriteCXRecord invalid cx: %d %d %d %g\n",
	      r->b, r->f, i, r->cx[i]);
    }
  }
  WSF1(r->cx, sizeof(double), cx_header.ne0);
#pragma omp atomic
  cx_header.length += m;

  return m;
}

int WriteRCRecord(TFILE *f, RC_RECORD *r) {
  int n, i, m0;
  int m = 0;

  if (rc_header.ntransitions == 0) {
    SetLockMPI();
    if (rc_header.ntransitions == 0) {
      fheader[DB_RC-1].nblocks++;
      FFLUSH(f);
      n = WriteRCHeader(f, &rc_header);
      FFLUSH(f);
    }
#pragma omp atomic
    rc_header.ntransitions += 1;
    ReleaseLockMPI();
  } else {
#pragma omp atomic
    rc_header.ntransitions += 1;
  }

  m0 = rc_header.nte*rc_header.nde;
#ifdef USEBF
  BFileCheckBuf(f,		
		sizeof(r->lower)+
		sizeof(r->upper)+
		sizeof(float)*m0);
#endif
  WSF0(r->lower);
  WSF0(r->upper);
  for (i = 0; i < m0; i++) {
    if (isnan(r->rc[i]) || isinf(r->rc[i])) {
      MPrintf(-1, "WriteRCRecord invalid strength: %d %d %d %d %d %d %g\n",
	      rc_header.type, rc_header.nte, rc_header.nde,
	      r->lower, r->upper, i, r->rc[i]);
    }
  }
  WSF1(r->rc, sizeof(float), m0);
  
#pragma omp atomic
  rc_header.length += m;

  return m;
}

int WriteRRRecord(TFILE *f, RR_RECORD *r) {
  int n, i;
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
    for (i = 0; i < m0; i++) {
      if (isnan(r->params[i]) || isinf(r->params[i])) {
	MPrintf(-1, "WriteRRRecord invalid param: %d %d %d %g\n",
		r->b, r->f, i, r->params[i]);
      }
    }
    WSF1(r->params, sizeof(float), m0);
  }
  m0 = rr_header.n_usr;
  for (i = 0; i < m0; i++) {
    if (isnan(r->strength[i]) || isinf(r->strength[i])) {
      MPrintf(-1, "WriteRRRecord invalid strength: %d %d %d %g\n",
	      r->b, r->f, i, r->strength[i]);
    }
  }
  WSF1(r->strength, sizeof(float), m0);

#pragma omp atomic
  rr_header.length += m;

  return m;
}

int WriteAIRecord(TFILE *f, AI_RECORD *r) {
  int n, m = 0;

  if (1.0 == 1.0 + (double)(r->rate)) {
    return 0;
  }
  
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
  if (isnan(r->rate) || isinf(r->rate)) {
    MPrintf(-1, "WriteAIRecord invalid rate: %d %d %g\n",
	    r->b, r->f, r->rate);
  }
  WSF0(r->rate);

#pragma omp atomic
  ai_header.length += m;
  return m;
}

int WriteAIMRecord(TFILE *f, AIM_RECORD *r) {
  int i, n, m = 0;

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
  
  for (i = 0; i < r->nsub; i++) {
    if (isnan(r->rate[i]) || isinf(r->rate[i])) {
      MPrintf(-1, "WriteAIMRecord invalid rate: %d %d %d %g\n",
	      r->b, r->f, i, r->rate[i]);
    }
  }
  WSF1(r->rate, sizeof(float), r->nsub);
  
#pragma omp atomic
  aim_header.length += m;

  return m;
}

int WriteCIRecord(TFILE *f, CI_RECORD *r) {
  int n, i;
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
  for (i = 0; i < m0; i++) {
    if (isnan(r->params[i]) || isinf(r->params[i])) {
      MPrintf(-1, "WriteCIRecord invalid param: %d %d %d %g\n",
	      r->b, r->f, i, r->params[i]);
    }
  }
  WSF1(r->params, sizeof(float), m0);
  m0 = ci_header.n_usr;
  for (i = 0; i < m0; i++) {
    if (isnan(r->strength[i]) || isinf(r->strength[i])) {
      MPrintf(-1, "WriteCIRecord invalid strength: %d %d %d %g\n",
	      r->b, r->f, i, r->strength[i]);
    }
  }
  WSF1(r->strength, sizeof(float), m0);

#pragma omp atomic
  ci_header.length += m;

  return m;
}

int WriteCIMRecord(TFILE *f, CIM_RECORD *r) {
  int n, i;
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
  for (i = 0; i < m0; i++) {
    if (isnan(r->strength[i]) || isinf(r->strength[i])) {
      MPrintf(-1, "WriteCIMRecord invalid strength: %d %d %d %g\n",
	      r->b, r->f, i, r->strength[i]);
    }
  }
  WSF1(r->strength, sizeof(float), m0);
  
#pragma omp atomic
  cim_header.length += m;

  return m;
}

int WriteSPRecord(TFILE *f, SP_RECORD *r, SP_EXTRA *rx, int utr) {
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
  if (utr) {
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
  if (isnan(r->strength) || isinf(r->strength)) {
    MPrintf(-1, "WriteSPRecord invalid strength: %d %d %g\n",
	    r->lower, r->upper, r->strength);
  }
  WSF0(r->strength);
  WSF0(r->rrate);
  WSF0(r->trate);
  if (utr) {
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

  if (r->j < 0) {
    if (r->ibase < 32767) {
      int i;
      for (i = 0; i < LNAME; i++) {
	if (r->name[i] == '\0') break;
	if (r->name[i] == '.') break;
	if (r->name[i] == '(') {
	  i = r->j;
	  r->j = r->ibase;
	  r->ibase = i;
	  break;
	}
      }
    }
    iuta = 1;
  }
  //else iuta = 0;
  
  if (_remove_closed) RemoveClosedShell(r);
  return m;
}

void RemoveClosedShell(EN_RECORD *r) {
  int ns, n, i, j, q, ic[16];
  char *s, b[LNAME], *c, *p, *d;
  
  p = r->ncomplex;
  strncpy(b, r->ncomplex, LNCOMPLEX);
  ns = StrSplit(b, '.');
  s = b;
  j = 0;
  n = 0;
  for (i = 0; i < 16; i++) ic[i] = 0;
  for (i = 0; i < ns; i++) {
    while (*s == ' ' || *s == '\t') s++;
    c = s;
    d = s;
    while (*s != '*') s++;
    s++;
    n = atoi(c);
    q = atoi(s);
    if (q < 2*n*n) {
      while (*c) {
	p[j++] = *(c++);
      }
      p[j++] = '.';
    } else {
      ic[n-1] = 1;
    }
    while (*s) s++;
    s++;
  }
  if (j > 0) p[j-1] = '\0';
  else {
    if (n == 0) p[j] = '\0';
    else {
      ic[n-1] = 0;
      while (*d) {
	p[j++] = *(d++);
      }
      p[j] = '\0';
    }
  }    
  
  p = r->sname;
  strncpy(b, r->sname, LSNAME);
  ns = StrSplit(b, '.');
  s = b;
  j = 0;
  for (i = 0; i < ns; i++) {
    while (*s == ' ' || *s == '\t') s++;
    n = atoi(s);
    if (ic[n-1] == 0) {
      while (*s) {
	p[j++] = *(s++);
      }
      p[j++] = '.';
    }
    while (*s) s++;
    s++;
  }
  if (j > 0) p[j-1] = '\0';
  else p[j] = '\0';
  
  p = r->name;
  strncpy(b, r->name, LNAME);
  ns = StrSplit(b, '.');
  s = b;
  j = 0;
  for (i = 0; i < ns; i++) {
    while (*s == ' ' || *s == '\t') s++;
    n = atoi(s);
    if (ic[n-1] == 0) {
      while (*s) {
	p[j++] = *(s++);
      }
      p[j++] = '.';
    }
    while (*s) s++;
    s++;
  }
  if (j > 0) p[j-1] = '\0';
  else p[j] = '\0';
}

int FindNRShells(int nele, EN_RECORD *r, int nm, short *nqc, short *nqs) {
  const int nmax=32;
  const int nmax2 = (nmax*(nmax+3))/2;
  int ncq[nmax], npq[nmax], nk[nmax];
  int nsq[nmax2], jmq[nmax2], jpq[nmax2];
  int nmj[nmax2], npj[nmax2], nmt[nmax2], npt[nmax2];
  int i, i0, n, ns, nq, k, kp, ik, nqt, nst, mk, kf, dq, wk, mi, ij;
  char c, *p0, *p1, nc0[LNCOMPLEX], sn0[LSNAME], nm0[LNAME], ss[16];

  if (nm > nmax) return -1;
  
  for (i = 0; i < nmax; i++) {
    ncq[i] = 0;
    npq[i] = 0;
  }
  if (nqs) {
    for (i = 0; i < nmax2; i++) {
      nsq[i] = 0;
      jmq[i] = 0;
      jpq[i] = 0;
      nmj[i] = 0;
      npj[i] = 0;
      nmt[i] = 0;
      npt[i] = 0;
    }
  }
  i0 = 0;
  n = 0;
  nqt = 0;
  p0 = r->ncomplex;
  p1 = nc0;
  for (i = 0; i < LNCOMPLEX; i++) {
    if (isspace(p0[i])) continue;
    if (p0[i] == '\0') break;
    *p1 = p0[i];
    p1++;
  }
  *p1 = '\0';

  ns = strlen(nc0);
  for (i = 0; i <= ns; i++) {
    c = nc0[i];
    if (c == '*') {
      n = atoi(&nc0[i0]);
      i0 = i+1;
    } else if (c == '.' || c == '\0') {
      nq = atoi(&nc0[i0]);
      ncq[n-1] = nq;
      nqt += nq;
      i0 = i+1;
      if (c == '\0') break;
    }
  }
  if (nqt > nele) return 1;
  if (nqt < nele) {
    mk = 8;
    mi = (1<<mk)-1;
    for (i = 0; i < mi; i++) {
      nst = 0;
      for (k = 0; k < mk; k++) {
	if ((i & (1<<k)) && ncq[k] == 0) {
	  nst += 2*(k+1)*(k+1);
	}
      }
      if (nst == nele-nqt) {
	for (k = 0; k < mk; k++) {
	  if ((i & (1<<k)) && ncq[k] == 0) {
	    ncq[k] = 2*(k+1)*(k+1);
	    nqt += ncq[k];
	  }
	}
	break;
      }
    }
    if (nqt != nele) return 2;
  }
  for (n = 1; n <= nm; n++) {
    nqc[n-1] = ncq[n-1];
  }
  
  if (nqs) {
    p0 = r->sname;
    p1 = sn0;
    for (i = 0; i < LSNAME; i++) {
      if (isspace(p0[i])) continue;
      if (p0[i] == '\0') break;
      *p1 = p0[i];
      p1++;
    }
    *p1 = '\0';
    
    i0 = 0;
    n = 0;
    k = -1;
    nst = 0;
    ns = strlen(sn0);
    for (i = 0; i <= ns; i++) {
      c = sn0[i];
      k = GetLFromSymbol(tolower(c));
      if (k >= 0) {
	n = atoi(&sn0[i0]);
	kp = k;
	ik = ((n+2)*(n-1))/2 + kp;
	i0 = i+1;
      } else if (c == 'a' || c == 'A') {
	n = atoi(&sn0[i0]);
	kp = n;
	ik = ((n+2)*(n-1))/2 + kp;
	i0 = i+1;
      } else if (c == '.' || c == '\0') {
	nq = atoi(&sn0[i0]);
	nsq[ik] = nq;
	npq[n-1] += nq;
	nst += nq;
	i0 = i+1;
	if (c == '\0') break;
      }
    }
  
    for (n = 1; n <= nmax; n++) {
      i0 = ((n+2)*(n-1))/2;
      if (ncq[n-1] > npq[n-1]) {
	dq = ncq[n-1] - npq[n-1];
	mk = -1;
	kf = -1;
	nqt = 0;
	for (k = 0; k < n; k++) {
	  if (nsq[k+i0] == 0) {
	    mk++;
	    nk[mk] = k;
	    wk = 2*(2*k+1);
	    if (wk > dq) {
	      nk[mk] = -1;
	      mk--;
	    }
	    nqt += wk;
	  }
	}
	mk++;
	if (nqt == dq) {
	  for (k = 0; k < mk; k++) {
	    wk = 2*(2*nk[k]+1);
	    nsq[nk[k]+i0] = wk;
	  }
	} else if (nqt < dq) {
	  return 3;
	} else {
	  mk = Min(mk,8);
	  mi = (1<<mk) - 1;
	  for (i = 0; i < mi; i++) {
	    nst = 0;
	    for (k = 0; k < mk; k++) {
	      if (i&(1<<k)) {
		nst += 2*(2*nk[k]+1);
	      }
	    }
	    if (nst == dq) {
	      for (k = 0; k < mk; k++) {
		if (i&(1<<k)) {
		  nsq[nk[k]+i0] = 2*(2*nk[k]+1);
		}
	      }
	      break;
	    }
	  }
	}
      }
    }
    if (nqs) {
      i0 = 0;
      for (n = 1; n <= nm; n++) {
	for (k = 0; k < n; k++) {
	  nqs[i0] = nsq[i0];
	  i0++;
	}
      }
    }
  }
  return 0;
}

int FillClosedShell(int nele, EN_RECORD *r, char *nc, char *sn, char *nm) {
  const int nmax=32;
  const int nmax2 = (nmax*(nmax+3))/2;
  int ncq[nmax], npq[nmax], nk[nmax];
  int nsq[nmax2], jmq[nmax2], jpq[nmax2];
  int nmj[nmax2], npj[nmax2], nmt[nmax2], npt[nmax2], nr[nmax2];
  int i, i0, n, ns, nq, k, kp, ik, nqt, nst, mk, kf, dq, wk, mi, ij;
  char c, c0, *p0, *p1, nc0[LNCOMPLEX], sn0[LSNAME], nm0[LNAME], ss[16];

  for (i = 0; i < nmax; i++) {
    ncq[i] = 0;
    npq[i] = 0;
  }
  for (i = 0; i < nmax2; i++) {
    nsq[i] = 0;
    jmq[i] = 0;
    jpq[i] = 0;
    nmj[i] = 0;
    npj[i] = 0;
    nmt[i] = 0;
    npt[i] = 0;
    nr[i] = 0;
  }
  i0 = 0;
  n = 0;
  nqt = 0;
  p0 = r->ncomplex;
  p1 = nc0;
  for (i = 0; i < LNCOMPLEX; i++) {
    if (isspace(p0[i])) continue;
    if (p0[i] == '\0') break;
    *p1 = p0[i];
    p1++;
  }
  *p1 = '\0';
  p0 = r->sname;
  p1 = sn0;
  for (i = 0; i < LSNAME; i++) {
    if (isspace(p0[i])) continue;
    if (p0[i] == '\0') break;
    *p1 = p0[i];
    p1++;
  }
  *p1 = '\0';
  p0 = r->name;
  p1 = nm0;
  for (i = 0; i < LNAME; i++) {
    if (isspace(p0[i])) continue;
    if (p0[i] == '\0') break;
    *p1 = p0[i];
    p1++;
  }
  *p1 = '\0';
  ns = strlen(nc0);
  for (i = 0; i <= ns; i++) {
    c = nc0[i];
    if (c == '*') {
      n = atoi(&nc0[i0]);
      i0 = i+1;
    } else if (c == '.' || c == '\0') {
      nq = atoi(&nc0[i0]);
      ncq[n-1] = nq;
      nqt += nq;
      i0 = i+1;
      if (c == '\0') break;
    }
  }
  if (nqt > nele) return 1;
  if (nqt < nele) {
    mk = 8;
    mi = (1<<mk)-1;
    for (i = 0; i < mi; i++) {
      nst = 0;
      for (k = 0; k < mk; k++) {
	if ((i & (1<<k)) && ncq[k] == 0) {
	  nst += 2*(k+1)*(k+1);
	}
      }
      if (nst == nele-nqt) {
	for (k = 0; k < mk; k++) {
	  if ((i & (1<<k)) && ncq[k] == 0) {
	    ncq[k] = 2*(k+1)*(k+1);
	    nqt += ncq[k];
	  }
	}
	break;
      }
    }
    if (nqt != nele) return 2;
  }
  i0 = 0;
  n = 0;
  k = -1;
  nst = 0;
  ns = strlen(sn0);
  for (i = 0; i <= ns; i++) {
    c = sn0[i];
    k = GetLFromSymbol(tolower(c));
    if (k >= 0) {
      n = atoi(&sn0[i0]);
      kp = k;
      ik = ((n+2)*(n-1))/2 + kp;
      i0 = i+1;
    } else if (c == 'a' || c == 'A') {
      n = atoi(&sn0[i0]);
      kp = n;
      ik = ((n+2)*(n-1))/2 + kp;
      i0 = i+1;
    } else if (c == '.' || c == '\0') {
      nq = atoi(&sn0[i0]);
      nsq[ik] = nq;
      npq[n-1] += nq;
      nst += nq;
      i0 = i+1;
      if (c == '\0') break;
    }
  }
  i0 = 0;
  n = 0;
  k = -1;
  ij = 0;
  ns = strlen(nm0);  
  for (i = 0; i <= ns; i++) {
    c = nm0[i];
    if (c == ')') {
      if (ij < 0) {
	nmj[ik] = atoi(&nm0[i0]);      
      } else if (ij > 0) {
	npj[ik] = atoi(&nm0[i0]);
      }
      i0 = i+1;
      continue;
    }
    if (c == '.' || c == '\0') {
      if (ij < 0) {
	nmt[ik] = atoi(&nm0[i0]);
      } else if (ij > 0) {
	npt[ik] = atoi(&nm0[i0]);
      } else {
	nq = atoi(&nm0[i0]);
	if (nq >= 0) {
	  jpq[ik] = nq;
	} else {
	  jmq[ik] = -nq;
	}
      }
      i0 = i+1;
      if (c == '\0') break;
      continue;
    }
    if (isupper(c)) {
      c0 = tolower(c);
    } else {
      c0 = c;
    }
    k = GetLFromSymbol(c0);
    if (k >= 0) {
      n = atoi(&nm0[i0]);
      kp = k;
      ik = ((n+2)*(n-1))/2 + kp;
      if (c != c0) {
	nr[ik] = 1;
      }
      i0 = i+1;
      continue;
    } else if (c == 'a' || c == 'A') {
      n = atoi(&nm0[i0]);
      kp = n;
      ik = ((n+2)*(n-1))/2 + kp;
      i0 = i+1;
      continue;
    }
    if (c == '(' || c == '\0') {
      nq = atoi(&nm0[i0]);      
      i0 = i+1;
      if (nq > 0) {
	ij = 1;
	jpq[ik] = nq;
      } else {
	ij = -1;
	jmq[ik] = -nq;
      }
      if (c == '\0') break;
      continue;
    }
  }
  
  for (n = 1; n <= nmax; n++) {
    i0 = ((n+2)*(n-1))/2;
    if (ncq[n-1] > npq[n-1]) {
      dq = ncq[n-1] - npq[n-1];
      mk = -1;
      kf = -1;
      nqt = 0;
      for (k = 0; k < n; k++) {
	if (nsq[k+i0] == 0) {
	  mk++;
	  nk[mk] = k;
	  wk = 2*(2*k+1);
	  if (wk > dq) {
	    nk[mk] = -1;
	    mk--;
	  }
	  nqt += wk;
	}
      }
      mk++;
      if (nqt == dq) {
	for (k = 0; k < mk; k++) {
	  wk = 2*(2*nk[k]+1);
	  nsq[nk[k]+i0] = wk;
	}
      } else if (nqt < dq) {
	return 3;
      } else {
	mk = Min(mk,8);
	mi = (1<<mk) - 1;
	for (i = 0; i < mi; i++) {
	  nst = 0;
	  for (k = 0; k < mk; k++) {
	    if (i&(1<<k)) {
	      nst += 2*(2*nk[k]+1);
	    }
	  }
	  if (nst == dq) {
	    for (k = 0; k < mk; k++) {
	      if (i&(1<<k)) {
		nsq[nk[k]+i0] = 2*(2*nk[k]+1);
	      }
	    }
	    break;
	  }
	}
      }
    }
    for (k = 0; k < n; k++) {
      if (nsq[i0+k] > jmq[i0+k]+jpq[i0+k]) {
	if (jmq[i0+k] == 0 && k > 0) jmq[i0+k] = (2*k-1)+1;
	if (jpq[i0+k] == 0) jpq[i0+k] = (2*k+1)+1;
	if (jmq[i0+k]+jpq[i0+k] > nsq[i0+k]) {
	  if (jmq[i0+k] == nsq[i0+k]) {
	    jpq[i0+k] = 0;
	  } else if (jpq[i0+k] == nsq[i0+k]) {
	    jmq[i0+k] = 0;
	  }
	}
      }
      if (jmq[i0+k] + jpq[i0+k] != nsq[i0+k]) return 4;	
    }
  }
  nc[0] = '\0';
  sn[0] = '\0';
  nm[0] = '\0';
  for (n = 1; n <= nmax; n++) {
    i0 = ((n+2)*(n-1))/2;
    if (ncq[n-1] > 0) {
      sprintf(nc, "%s%d*%d.", nc, n, ncq[n-1]);
    }
    for (k = 0; k < n; k++) {
      if (nsq[i0+k] > 0) {
	SpecSymbol(ss, k);
	if (nr[i0+k]) {
	  ss[0] = toupper(ss[0]);
	}
	sprintf(sn, "%s%d%s%d.", sn, n, ss, nsq[i0+k]);
	if (nr[i0+k] == 0) {
	  if (jmq[i0+k] > 0) {
	    if (ij) {
	      sprintf(nm, "%s%d%s-%d(%d)%d.",
		      nm, n, ss, jmq[i0+k], nmj[i0+k], nmt[i0+k]);
	    } else {
	      sprintf(nm, "%s%d%s-%d.", nm, n, ss, jmq[i0+k]);
	    }
	  }
	  if (jpq[i0+k] > 0) {
	    if (ij) {
	      sprintf(nm, "%s%d%s+%d(%d)%d.",
		      nm, n, ss, jpq[i0+k], npj[i0+k], npt[i0+k]);
	    } else {
	      sprintf(nm, "%s%d%s+%d.", nm, n, ss, jpq[i0+k]);
	    }
	  }
	} else {
	  sprintf(nm, "%s%d%s%d.", nm, n, ss, nsq[i0+k]);
	}
      }
    }
    if (nsq[i0+n] > 0) {
      sprintf(sn, "%s%da%d.", sn, n, nsq[i0+n]);
      sprintf(nm, "%s%da%d.", nm, n, nsq[i0+n]);
    }
  }
  n = strlen(nc);
  if (n > 0) nc[n-1] = '\0';
  n = strlen(sn);
  if (n > 0) sn[n-1] = '\0';
  n = strlen(nm);
  if (n > 0) nm[n-1] = '\0';
  
  return 0;
}

int ReadENFRecord(TFILE *f, ENF_RECORD *r, int swp) {
  int n, m = 0;

  RSF0(r->ilev);
  RSF0(r->energy);
  RSF0(r->pbasis);
  
  if (swp) SwapEndianENFRecord(r);

  return m;
}

int ReadTRHeader(TFILE *f, TR_HEADER *h, int swp, int *utr) {
  int n, m = 0;

  *utr = 0;
  if (version_read[DB_TR-1] < 109) return ReadTRHeaderOld(f, h, swp);

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->gauge);
  RSF0(h->mode);
  RSF0(h->multipole);

  if (swp) SwapEndianTRHeader(h);

  if (h->length/h->ntransitions > SIZE_TR_RECORD) *utr = 1;
  
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

int ReadTRRecord(TFILE *f, TR_RECORD *r, TR_EXTRA *rx, int swp, int utr) {
  int n, m = 0;
    
  if (version_read[DB_TR-1] < 109) return ReadTRRecordOld(f, r, rx, swp);

  RSF0(r->lower);
  RSF0(r->upper);
  RSF0(r->strength);

  if (rx) {
    rx->energy = 0.0;
    rx->sdev = 0.0;
    rx->sci = 1.0;
  }
  if (utr) {
    RSF0(rx->energy);
    RSF0(rx->sdev);
    RSF0(rx->sci);
  }

  if (swp) SwapEndianTRRecord(r, rx, utr);

  //if (utaci == 0 && rx) {
  //  rx->sci = 1.0;
  //}

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

int ReadRCHeader(TFILE *f, RC_HEADER *h, int swp) {
  int n, m = 0;
  
  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->type);
  RSF0(h->nexc);
  RSF0(h->mexc);
  RSF0(h->ncap);
  RSF0(h->nte);
  RSF0(h->nde);
  RSF0(h->te0);
  RSF0(h->dte);
  RSF0(h->de0);
  RSF0(h->dde);

  if (swp) SwapEndianRCHeader(h);
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

int ReadRCRecord(TFILE *f, RC_RECORD *r, int swp, RC_HEADER *h) {
  int m = 0;
  int n, i, m0;
  
  RSF0(r->lower);
  RSF0(r->upper);

  if (swp) SwapEndianRCRecord(r);

  m0 = h->nte*h->nde;
  r->rc = (float *) malloc(sizeof(float)*m0);
  RSF1(r->rc, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->rc[i]), sizeof(float));
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

int ReadSPHeader(TFILE *f, SP_HEADER *h, int swp, int *utr) {
  int n, m = 0;

  *utr = 0;
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

  if (h->length/h->ntransitions > SIZE_SP_RECORD) *utr = 1;
  
  return m;
}

int ReadSPRecord(TFILE *f, SP_RECORD *r, SP_EXTRA *rx, int swp, int utr) {
  int n, m = 0;

  if (version_read[DB_SP-1] < 109) return ReadSPRecordOld(f, r, rx, swp);

  RSF0(r->lower);
  RSF0(r->upper);
  RSF0(r->energy);
  RSF0(r->strength);
  RSF0(r->rrate);
  RSF0(r->trate);
  if (utr) {
    RSF0(rx->sdev);
  }

  if (swp) SwapEndianSPRecord(r, rx, utr);

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
  fheader[ihdr].atom = fhdr->atom;
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
  RC_HEADER *rc_hdr;
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
  case DB_RC:
    rc_hdr = (RC_HEADER *) rhdr;
    memcpy(&rc_header, rc_hdr, sizeof(RC_HEADER));
    rc_header.position = p;
    rc_header.length = 0;
    rc_header.ntransitions = 0;
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
  case DB_RC:
    FSEEK(f, rc_header.position, SEEK_SET);
    if (rc_header.length > 0 || rc_header.nde == 0) {
      n = WriteRCHeader(f, &rc_header);
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
  char tfn[1024];
  int n, swp, v, vs, i;

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
	    || (fh.type >= DB_RO && fh.type < DB_RC))) {
    if (mem_en_table == NULL) {
      if (fh.type == DB_EN) {
	FCLOSE(f1);
	MemENTable(ifn);
	f1 = FOPEN(ifn, "r");
	n = ReadFHeader(f1, &fh, &swp);
      } else {
	ConstructFileNameEN(ifn, tfn);
	MemENTable(tfn);
      }
      if (mem_en_table == NULL) {
	printf("Energy table has not been built in memory.\n");
	goto DONE;
      }
    }
  }

  if (v && fh.type > DB_CIM && fh.type < DB_RO) {
    if (mem_enf_table == NULL) {
      if (fh.type == DB_ENF) {
	FCLOSE(f1);
	MemENFTable(ifn);
	f1 = FOPEN(ifn, "r");
	n = ReadFHeader(f1, &fh, &swp);
      } else {
	ConstructFileNameEN(ifn, tfn);
	MemENFTable(tfn);
      }
      if (mem_enf_table == NULL) {
	printf("Field dependent energy table has not been built in memory.\n");
	goto DONE;
      }
    }
  }

  if (fh.nthreads > 0 || iuta) {
    fprintf(f2, "FAC %d.%d.%d[%d.%d.%d]\n",
	    fh.version, fh.sversion, fh.ssversion, fh.nthreads, iuta, utaci);
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
  case DB_RC:
    n = PrintRCTable(f1, f2, v, vs, swp);
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
    printf("%s File type is not DB_EN\n", fn);
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
    printf("%s File type is not DB_EN\n", fn);
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

void ConstructFileNameEN(char *ifn, char *tfn) {
  int i, n;

  n = strlen(ifn);
  n = Min(n, 1024);
  i = 0;
  while (i < n && ifn[i]) {
    tfn[i] = ifn[i];
    if (ifn[i] == '.') {
      tfn[i+1] = 'e';
      tfn[i+2] = 'n';
      tfn[i+3] = '\0';
      break;
    }
    i++;
  }
}

EN_SRECORD *GetOrLoadMemENTable(int *s, char *fn) {
  char efn[1024];
  if (mem_en_table && mem_en_table_size > 0) {
    *s = mem_en_table_size;
    return mem_en_table;
  }
  *s = 0;
  if (fn == NULL) return NULL;
  ConstructFileNameEN(fn, efn);
  MemENTable(efn);
  *s = mem_en_table_size;
  return mem_en_table;
}

EN_SRECORD *GetMemENFTable(int *s) {
  *s = mem_enf_table_size;
  return mem_enf_table;
}


EN_SRECORD *GetOrLoadMemENFTable(int *s, char *fn) {
  char efn[1024];
  if (mem_enf_table && mem_enf_table_size > 0) {
    *s = mem_enf_table_size;
    return mem_enf_table;
  }
  ConstructFileNameEN(fn, efn);
  MemENFTable(efn);
  *s = mem_enf_table_size;
  return mem_enf_table;
}

int JFromENRecord(EN_RECORD *r) {
  if (r->j < 0) return r->ibase;
  return r->j;
}

double WFromENRecord(EN_RECORD *r) {
  double w, k;
  if (r->j >= 0) w = 1.0 + r->j;
  else {
    k = (-r->j)/10;
    w = 1.0 + r->ibase;    
    if (k > 0) {
      w += k*0x80000000;
    }
  }
  return w;
}

int IBaseFromENRecord(EN_RECORD *r) {
  if (r->j < 0) return r->j;
  else return r->ibase;
}

int MemENTable(char *fn) {
  return MemENTableWC(fn, 0, NULL, NULL);
}

int MemENTableWC(char *fn, int k0, int *ifk, short ***nc) {
  F_HEADER fh;  
  EN_HEADER h;
  EN_RECORD r;
  NCOMPLEX ncomplex[MAXNCOMPLEX];
  TFILE *f;
  char *s;
  int n, i, nlevels, j, k;
  float e0;
  int swp, sr, ic, kk;

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

  for (i = 0; i < N_ELEMENTS1; i++) _eground[i] = 0.0;
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
      if (r.energy < _eground[h.nele]) _eground[h.nele] = r.energy;
      mem_en_table[r.ilev].energy = r.energy;
      mem_en_table[r.ilev].p = r.p;
      if (r.j < 0) {
	j = -((-r.j)%10);
	k = (-r.j)/10;
	mem_en_table[r.ilev].k = k;
	mem_en_table[r.ilev].j = r.ibase;
	if (j == -2) {
	  mem_en_table[r.ilev].ibase = -1;
	} else {
	  mem_en_table[r.ilev].ibase = -nlevels-1;
	}
      } else {
	mem_en_table[r.ilev].k = 0;
	mem_en_table[r.ilev].j = r.j;
	mem_en_table[r.ilev].ibase = r.ibase;
      }
      if (nc && h.nele >= k0) {
	kk = h.nele-k0;
	GetNComplex(ncomplex, r.ncomplex);
	for (ic = 0; ic < MAXNCOMPLEX; ic++) {
	  if (ncomplex[ic].n > 0 && ncomplex[ic].nq > 0) {
	    nc[kk][r.ilev-ifk[kk]][ncomplex[ic].n-1] = ncomplex[ic].nq;
	  }
	}
      }
    }
  }

  /*
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
  */
  FCLOSE(f);
  return 0;
}    

double GroundEnergy(int k) {
  return _eground[k];
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
    fprintf(f2, "  ILEV  IBASE    ENERGY       P   VNL         2J\n");
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
      fprintf(f2, "%6d %6d %15.8E %1d %5d %10d %-32s %-48s %-s\n",
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
  int nb, nh, utr;
  double e, a, gf;

  nb = 0;
  
  while (1) {
    nh = ReadTRHeader(f1, &h, swp, &utr);
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
	n = ReadTRRecord(f1, &r, &rx, swp, utr);
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
      n = ReadTRRecord(f1, &r, &rx, swp, utr);
      if (n == 0) break;
      if (utr) {
	if (v) {
	  if (rx.energy <= 0 || rx.sdev <= 0) {
	    e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	  } else {
	    e = rx.energy;
	  }
	  gf = OscillatorStrength(h.multipole, e, r.strength, &a);
	  a /= (mem_en_table[r.upper].j + 1.0);
	  a *= RATE_AU;
	  fprintf(f2, "%6d %10d %6d %10d %13.6E %11.4E %13.6E %13.6E %13.6E %10.3E\n",
		  r.upper, mem_en_table[r.upper].j, 
		  r.lower, mem_en_table[r.lower].j,
		  (e*HARTREE_EV), 
		  (rx.sdev*HARTREE_EV), gf, a, r.strength, rx.sci);
	} else {
	  e = rx.energy;
	  fprintf(f2, "%6d %6d %13.6E %11.4E %13.6E %10.3E\n",
		  r.upper, r.lower, e, rx.sdev, r.strength, rx.sci);
	}
      } else {
	if (v) {
	  e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	  gf = OscillatorStrength(h.multipole, e, (double)r.strength, &a);
	  a /= (mem_en_table[r.upper].j + 1.0);
	  a *= RATE_AU;
	  fprintf(f2, "%6d %6d %6d %6d %13.6E %13.6E %13.6E %13.6E\n",
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
	free(r.strength);
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
  int utr;
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadTRHeader(f, &h, swp, &utr);
    if (n == 0) break;
    for (k = 0; k < h.ntransitions; k++) {
      n = ReadTRRecord(f, &r, &rx, swp, utr);
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
  double bte, bms, be, eu;

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
      if (v && h.tegrid[0] >= 0) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "TE0\t= %15.8E\n", h.te0 * HARTREE_EV);
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v && h.tegrid[0] >= 0) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v && h.tegrid[0] >= 0) {
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
	if (h.msub || h.qk_mode == QK_FIT) free(r.params);
	free(r.strength);
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
	fprintf(f2, "%6d %10d %6d %10d %11.4E %d\n",
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
	    eu = h.usr_egrid[t];
	    if (h.tegrid[0] < 0) eu *= be;
	    a = eu;
	    if (h.usr_egrid_type == 1) a += be;
	    a *= 2.0*(1.0 + 0.5*FINE_STRUCTURE_CONST2 * a);
	    a = PI * AREA_AU20/a;
	    if (!h.msub) a /= (mem_en_table[r.lower].j+1.0);
	    a *= r.strength[p2];
	    fprintf(f2, "%11.4E %11.4E %11.4E\n",
		    eu*HARTREE_EV, r.strength[p2], a);
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
	free(r.strength);
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
	free(r.strength);
	free(r.bethe);
	free(r.born);
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
	free(r.nk);
	free(r.nq);
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
	free(r.cx);
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

int PrintRCTable(TFILE *f1, FILE *f2, int v, int vs, int swp) {
  RC_HEADER h;
  RC_RECORD r;
  int n, i, nb, nh, k, ite, ide;
  double e;

  nb = 0;
  while (1) {
    nh = ReadRCHeader(f1, &h, swp);
    if (nh == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "TYPE\t= %d\n", h.type);
    fprintf(f2, "NEXC\t= %d\n", h.nexc);
    fprintf(f2, "MEXC\t= %d\n", h.mexc);
    fprintf(f2, "NCAP\t= %d\n", h.ncap);
    fprintf(f2, "NTE\t= %d\n", h.nte);
    fprintf(f2, "TE0\t= %15.8E\n", h.te0);
    e = h.dte/100;
    h.dte = (e - (int)e)*100;
    e = ((int)(e))/1e4;
    fprintf(f2, "DTE\t= %15.8E\n", h.dte);
    fprintf(f2, "TZS\t= %15.8E\n", e);
    fprintf(f2, "NDE\t= %d\n", h.nde);
    fprintf(f2, "DE0\t= %15.8E\n", h.de0);
    fprintf(f2, "DDE\t= %15.8E\n", h.dde);
        
    IDX_RECORD *idx = NULL;
    if (vs && h.ntransitions > 1) {
      idx = malloc(sizeof(IDX_RECORD)*h.ntransitions);
      idx[0].position = h.position + nh;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadRCRecord(f1, &r, swp, &h);
	if (n == 0) break;
	if (i < h.ntransitions-1) {
	  idx[i+1].position = idx[i].position + n;
	}
	idx[i].i0 = r.lower;
	idx[i].i1 = r.upper;
	free(r.rc);
      }
      qsort(idx, h.ntransitions, sizeof(IDX_RECORD), CompIdxRecord);
      FSEEK(f1, h.position+nh, SEEK_SET);
    }
    for (i = 0; i < h.ntransitions; i++) {
      if (idx) {
	FSEEK(f1, idx[i].position, SEEK_SET);
      }
      n = ReadRCRecord(f1, &r, swp, &h);
      if (n == 0) break;
      k = 0;
      for (ide = 0; ide < h.nde; ide++) {
	fprintf(f2, "%d %4d %6d %6d", h.type, ide, r.lower, r.upper);
	for (ite = 0; ite < h.nte; ite++) {
	  fprintf(f2, " %12.5E", r.rc[k]);
	  k++;
	}
	fprintf(f2, "\n");
      }	
      free(r.rc);
    }
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
  float e, eph, ee, phi, rr, eu;

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
      if (v && h.tegrid[0] >= 0) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v && h.tegrid[0] >= 0) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v && h.tegrid[0] >= 0) {
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
	if (h.qk_mode == QK_FIT) free(r.params);
	free(r.strength);
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
	if (r.f >= 0) {
	  e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	  fprintf(f2, "%6d %10d %6d %10d %11.4E %6d\n",
		  r.b, mem_en_table[r.b].j, 
		  r.f, mem_en_table[r.f].j,
		  (e*HARTREE_EV), r.kl);
	} else {
	  e = (double)(*((float *) &r.kl));
	  r.kl = -r.f;
	  fprintf(f2, "%6d %10d %6d %10d %11.4E %6d\n",
		  r.b, mem_en_table[r.b].j, 
		  r.f, -1,
		  (e*HARTREE_EV), r.kl);
	}
      } else {
	fprintf(f2, "%6d %6d %6d\n", r.b, r.f, r.kl);
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
	    if (h.tegrid[0] < 0) eph *= e;
	    eu = eph;
	    ee = eph - e;
	  } else {
	    ee = h.usr_egrid[t];
	    if (h.tegrid[0] < 0) ee *= e;
	    eu = ee;
	    eph = ee + e;
	  }
	  phi = 2.0*PI*FINE_STRUCTURE_CONST*r.strength[t]*AREA_AU20;
	  rr = phi * pow(FINE_STRUCTURE_CONST*eph, 2) / (2.0*ee);
	  rr /= 1.0+0.5*FINE_STRUCTURE_CONST2*ee;
	  phi /= (mem_en_table[r.b].j + 1.0);
	  if (r.f >= 0) {
	    rr /= (mem_en_table[r.f].j + 1.0);
	  }
	  fprintf(f2, "%11.4E %11.4E %11.4E %11.4E\n",
		  eu*HARTREE_EV, rr, phi, r.strength[t]);
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
	if (r.f >= 0) {
	  e = mem_en_table[r.b].energy - mem_en_table[r.f].energy;
	} else {
	  e = h.egrid[0];
	}
	if (e < 0) er = e - h.emin;
	else er = e;
	if (r.f >= 0) {
	  sdr = 0.5*(mem_en_table[r.b].j + 1.0);
	  sdr *= PI*PI*r.rate/(er*(mem_en_table[r.f].j + 1.0));
	  sdr *= AREA_AU20*HARTREE_EV;
	} else {
	  sdr = 0.0;
	}
	fprintf(f2, "%6d %10d %6d %10d %11.4E %11.4E %11.4E\n",
		r.b, mem_en_table[r.b].j,
		r.f, r.f>=0?mem_en_table[r.f].j:-1,
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
	free(r.rate);
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
  float e, a, e1;
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
      if (v && h.tegrid[0] >= 0) {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v && h.tegrid[0] >= 0) {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
	fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v && h.tegrid[0] >= 0) {
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
	free(r.params); 
	free(r.strength);
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
	fprintf(f2, "%6d %10d %6d %10d %11.4E %2d\n",
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
	  e1 = h.usr_egrid[t];
	  if (h.tegrid[0] < 0) {
	    e1 *= be;
	  }
	  a = e1;
	  if (h.usr_egrid_type == 1) a += be;
	  a *= 1.0 + 0.5*FINE_STRUCTURE_CONST2*a;
	  a = AREA_AU20/(2.0*a*(mem_en_table[r.b].j + 1.0));
	  a *= r.strength[t];
	  fprintf(f2, "%11.4E %11.4E %11.4E\n",
		  e1*HARTREE_EV, r.strength[t], a);
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
	free(r.strength);
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
  int nb, nh, utr;
  float e, a;

  nb = 0;
  
  while (1) {
    nh = ReadSPHeader(f1, &h, swp, &utr);
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
	n = ReadSPRecord(f1, &r, &rx, swp, utr);
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
      n = ReadSPRecord(f1, &r, &rx, swp, utr);
      if (n == 0) break;
      e = r.energy;
      if (v) e *= HARTREE_EV;
      a = r.strength;
      if (utr) {
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
  int n, swp, utr;
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
    n = ReadSPHeader(f, &h, swp, &utr);
    if (n == 0) break;
    if (h.type != 0 || h.nele != k) {
      FSEEK(f, h.length, SEEK_CUR);
      continue;
    }
    
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f, &r, &rx, swp, utr);
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
  int n, swp, utr;
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
    n = ReadSPHeader(f, &h, swp, &utr);
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
      n = ReadSPRecord(f, &r, &rx, swp, utr);
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

int CompareENRecordEnergySU(const void *p0, const void *p1) {
  EN_RECORD *r0, *r1;

  r0 = (EN_RECORD *) p0;
  r1 = (EN_RECORD *) p1;

  if (r0->j != -1 && r1->j == -1) return -1;
  else if (r0->j == -1 && r1->j != -1) return 1;
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

  if (r0->j < 0 && r1->j >= 0) return 1;
  if (r0->j >= 0 && r1->j < 0) return -1;
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

int CompareENSName(const void *c1, const void *c2) {
  EN_RECORD *r1, *r2;

  r1 = (EN_RECORD *) c1;
  r2 = (EN_RECORD *) c2;
  return strcmp(r1->sname, r2->sname);
}

int CompareENName(const void *c1, const void *c2) {
  EN_RECORD *r1, *r2;

  r1 = (EN_RECORD *) c1;
  r2 = (EN_RECORD *) c2;
  return strcmp(r1->name, r2->name);
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

int SortUniqSName(int n, EN_RECORD *a) {
  int i, j;
  EN_RECORD b;

  qsort(a, n, sizeof(EN_RECORD), CompareENSName);
  j = 1;
  memcpy(&b, &a[0], sizeof(EN_RECORD));
  for (i = 1; i < n; i++) {
    if (CompareENSName(&a[i], &b) != 0) {
      if (i != j) {
	memcpy(&a[j], &a[i], sizeof(EN_RECORD));
      }
      memcpy(&b, &a[i], sizeof(EN_RECORD));
      j++;
    }
  }

  return j;
}

int MatchLevelsPJ(int n0, EN_RECORD *r0, int n2, EN_RECORD *r1) {
  int i0, i1, j, n, n1, im;
  double de, a, b;
  
  n1 = 0;
  for (i1 = 0; i1 < n2; i1++) {
    for (j = 0; j < n0; j++) {
      if (strcmp(r0[j].sname, r1[i1].sname) == 0) break;
    }
    if (j == n0) {
      r1[i1].j = -(r1[i1].j+1);
    } else {
      n1++;
    }
  }

  qsort(r1, n2, sizeof(EN_RECORD), CompareENRecord);

  i0 = 0;
  i1 = 0;
  n = 0;
  de = 0.0;
  while (i0 < n0 && i1 < n1) {
    if (n0-i0 == n1-i1) {
      for (; i0 < n0; i0++, i1++) {	
	de += r0[i0].energy - r1[i1].energy;
	n++;
      }
      break;
    }
    if (r0[i0].energy < r1[i1].energy - _cmpetol) {
      r0[i0].j = -(r0[i0].j+1);
      i0++;
      continue;
    }
    if (r0[i0].energy > r1[i1].energy + _cmpetol) {
      r1[i1].j = -(r1[i1].j+1);
      i1++;
      continue;
    }
    de += r0[i0].energy - r1[i1].energy;
    n++;
    i0++;
    i1++;    
  }
  for (; i0 < n0; i0++) r0[i0].j = -(r0[i0].j+1);
  for (; i1 < n1; i1++) r1[i1].j = -(r1[i1].j+1);
  if (n > 0) de /= n;
  for (i0 = 0; i0 < n0; i0++) {
    if (r0[i0].j < 0) {
      a = 1e30;
      im = -1;
      for (i1 = n1; i1 < n2; i1++) {
	if (r1[i1].j >= 0) continue;
	b = fabs(r1[i1].energy+de-r0[i0].energy);
	if (b < a) {
	  a = b;
	  im = i1;
	}
      }
      if (im >= n1) {
	r0[i0].j = -(r0[i0].j+1);
	r1[im].j = -(r1[im].j+1);
	n++;
      }      
    }
  }
  qsort(r0, n0, sizeof(EN_RECORD), CompareENRecord);
  qsort(r1, n2, sizeof(EN_RECORD), CompareENRecord);

  return n;
}
  
int FindLevelBlock(int n0, EN_RECORD *r0, EN_RECORD **r1p, 
		   int nele, char *ifn, int *nbs) {
  F_HEADER fh;
  EN_HEADER h;
  EN_RECORD g, *r1, *r0r, *r0c, *r1r;
  TFILE *f;
  int i, k, j, jp, nr, nb, nb0, nv, ni, nj, n1;
  int swp, sfh, utalev;
  int mk0[1024], mk1[1024];

  *r1p = NULL;
  if (n0 <= 0) return 0;
  
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
  utalev = 0;
  for (i = 0; i < n0; i++) {
    if (r0[i].j < 0) utalev = 1;
    nv = abs(r0[i].p);    
    j = nv/100;
    nv = nv%100;
    if (mk0[j] < nv) mk0[j] = nv;
  }

  n1 = 2*n0;
  r1 = malloc(sizeof(EN_RECORD)*n1);
  r1r = malloc(sizeof(EN_RECORD)*n1);
  r0c = malloc(sizeof(EN_RECORD)*n0);
  r0r = malloc(sizeof(EN_RECORD)*n0);
  memcpy(r0r, r0, sizeof(EN_RECORD)*n0);
  for (i = 0; i < n0; i++) {
    RemoveClosedShell(&r0r[i]);
    if (r0r[i].j < 0) r0r[i].j = 1;
  }
  memcpy(r0c, r0r, sizeof(EN_RECORD)*n0);
  int n0c = SortUniqNComplex(n0, r0c);
  
  k = 0;
  nb0 = -1;
  for (nb = 0; nb < fh.nblocks; nb++) {
    nr = ReadENHeader(f, &h, swp);
    if (h.nele != nele || nb < *nbs) {
      FSEEK(f, h.length, SEEK_CUR);
      continue;
    }
    for (i = 0; i < h.nlevels; i++) {
      nr = ReadENRecord(f, &r1[k], swp);
      memcpy(&r1r[k], &r1[k], sizeof(EN_RECORD));
      RemoveClosedShell(&r1r[k]);
      if (r1r[k].j < 0) r1r[k].j = 1;
      for (j = 0; j < n0c; j++) {
	if (strcmp(r1r[k].ncomplex, r0c[j].ncomplex) == 0) {
	  break;
	}
      }
      if (j < n0c) {
	nv = abs(r1[k].p);
	j = nv/100;
	nv = nv%100;
	if (mk0[j] >= nv) {
	  if (mk1[j] < nv) mk1[j] = nv;
	  if (k == 0) nb0 = nb;
	  k++;
	  if (k == n1) {
	    n1 += n0;
	    r1 = ReallocNew(r1, sizeof(EN_RECORD)*n1);
	    r1r = ReallocNew(r1r, sizeof(EN_RECORD)*n1);
	  }
	}
      }
    }
    if (_cmpnbm >= 0 && k > 0 && nb-nb0+1 >= _cmpnbm) break;
    if (k >= n0) break;
  }
  *nbs = nb0 + 1;
  FCLOSE(f);
  free(r0c);
  if (n0 <= 0 || k <= 0) {
    *r1p = r1;
    if (n0) free(r0r);
    if (n1) free(r1r);
    return 0;
  }
  n1 = k;
  for (i = 0; i < n0; i++) {
    r0r[i].ilev = i;
  }
  for (i = 0; i < n1; i++) {
    r1r[i].ilev = i;
  }
  //printf("flb: %d %d %d\n", utalev, n0, n1);
  if (utalev) {
    qsort(r0r, n0, sizeof(EN_RECORD), CompareENName);
    qsort(r1r, n1, sizeof(EN_RECORD), CompareENName);
    i = 0;
    j = 0;
    nr = 0;
    while (i < n0 && j < n1) {      
      k = CompareENName(&r0r[i], &r1r[j]);
      if (k > 0) {
	r1r[j].j = -(r1r[j].j+1);
	j++;	
      } else if (k < 0) {
	r0r[i].j = -(r0r[i].j+1);
	i++;
      } else {
	i++;
	j++;
	nr++;
      }
    }
  } else {  
    for (i = 0; i < 1024; i++) {
      if (mk1[i] > mk0[i]) mk1[i] = mk0[i];
    }
    int nk0 = 0;
    for (i = 0; i < n0; i++) {
      nv = abs(r0r[i].p);
      j = nv/100;
      nv = nv%100;
      if (nv > mk1[j]) {
	r0r[i].j = -(r0r[i].j+1);
      } else {      
	nk0++;
      }
    }
    int nk1 = 0;
    for (i = 0; i < n1; i++) {
      nv = abs(r1r[i].p);
      j = nv/100;
      nv = nv%100;
      if (nv > mk1[j]) {
	r1r[i].j = -(r1r[i].j+1);
      } else {
	nk1++;
      }
    }
    qsort(r0r, n0, sizeof(EN_RECORD), CompareENRecord);
    qsort(r1r, n1, sizeof(EN_RECORD), CompareENRecord);
    
    i = 0;
    j = 0;
    nr = 0;
    while (i < nk0 && j < nk1) {
      for (k = i+1; k < nk0; k++) {
	if (r0r[k].j != r0r[i].j ||
	    r0r[k].p*r0r[i].p < 0) break;
      }
      ni = k-i;
      k = 0;
      for (jp = j; jp < nk1; jp++) {
	if (r1r[jp].j < r0r[i].j) break;
	if (r1r[jp].j == r0r[i].j) {
	  if (r1r[jp].p < 0 && r0r[i].p > 0) break;
	  if (r1r[jp].p*r0r[i].p > 0) {
	    k = 1;
	    break;
	  }
	}
      }
      if (k == 0) {
	k = i+ni;
	for (; i < k; i++) {
	  r0r[i].j = -(r0r[i].j+1);
	}
	continue;
      }
      for (; j < jp; j++) {
	r1r[j].j = -(r1r[j].j+1);
      }
      for (k = j+1; k < nk1; k++) {
	if (r1r[k].j != r1r[j].j ||
	    r1r[k].p*r1r[j].p < 0) break;
      }
      nj = k-j;
      nr += MatchLevelsPJ(ni, &r0r[i], nj, &r1r[j]);
      i += ni;
      j += nj;
    }
    n0 = nk0;
    n1 = nk1;
  }
  for (; i < n0; i++) {
    r0r[i].j = -(r0r[i].j+1);
  }
  for (; j < n1; j++) {
    r1r[j].j = -(r1r[j].j+1);
  }  
  qsort(r0r, n0, sizeof(EN_RECORD), CompareENRecord);
  qsort(r1r, n1, sizeof(EN_RECORD), CompareENRecord);
  if (utalev) {
    for (i = 0; i < n1; i++) r1[i].j = -1;
  }
  if (nr > 0) {
    for (i = 0; i < nr; i++) {
      memcpy(&r0r[i], &r0[r0r[i].ilev], sizeof(EN_RECORD));
      memcpy(&r1r[i], &r1[r1r[i].ilev], sizeof(EN_RECORD));
    }
    memcpy(r0, r0r, sizeof(EN_RECORD)*nr);
    memcpy(r1, r1r, sizeof(EN_RECORD)*nr);
  }
  free(r0r);
  free(r1r);
  *r1p = r1;
  return nr;
}

int GetNComplex(NCOMPLEX *c, char *s) {
  int i, n, nq;
  char *p, buf[8192];
  n = strlen(s);
  for (i = 0; i <= n; i++) {
    if (s[i] == '.') buf[i] = ' ';
    else buf[i] = s[i];
  }
  s = buf;
  i = 0;
  while (1) {
    if (i == MAXNCOMPLEX-1) {
      printf("Num of NCOMPLEX shells exceeded the limit %d\n", MAXNCOMPLEX-1);
      exit(1);
    }
    n = strtol(s, &p, 10);
    if (n == 0) {
      int j;
      for (j = i; j < MAXNCOMPLEX; j++) {
	c[j].n = 0;
	c[j].nq = 0;
      }
      return i;
    }
    s = p+1;
    nq = strtol(s, &p, 10);
    c[i].n = n;
    c[i].nq = nq;
    s = p;
    i++;
  }
  return i;
}

int ChannelAI(int b, int nmb, short *ncb,
	      int f, int nmf, short *ncf,
	      int *cn0, int *cn1, int *cn2) {
  int n, n0, n1, n2, nm;

  nm = Min(nmb, nmf);
  n0 = -1;
  n1 = -1;
  n2 = -1;
  for (n = 0; n < nm; n++) {
    if (ncb[n] == ncf[n]) continue;
    if (ncb[n] < ncf[n] && n0 < 0) {
      n0 = n;
      continue;
    }
    if (ncb[n] > ncf[n]) {
      if (n1 < 0) n1 = n;
      else n2 = n;
      if (ncb[n] == ncf[n]+2 && n2 < 0) n2 = n;
    }
  }
  if (n2 < 0 && n < nmb) n2 = nmb;
  n0++;
  n1++;
  n2++;
  *cn0 = n0;
  *cn1 = n1;
  *cn2 = n2;

  if (n0 > 0 && n1 > 0 && n0 < 10 && n1 < 10 && n2 > 0) {
    n = _ix_iai[n0*10+n1];
    if (n >= 0) {
      return (n2-1)*_nc_iai+_ix_iai[n0*10+n1];
    }
  }
  return -1;
}

int JoinDBase(char *pref, int nk, int *ks, int ic) {
  char ifn[1024], ofn[1024], buf[1024], a[16];
  int k, i, j, n, nb, swp, k0, k1, z, nth, tlevs, nt;
  double wt0, wt1, tt0, tt1;
  F_HEADER fh, fh1[7];
  EN_HEADER h0;
  EN_RECORD r0;
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
  RC_HEADER h6;
  RC_RECORD r6;
  char ext[7][3] = {"en", "tr", "ce", "rr", "ci", "ai", "rc"};
  int types[7] = {DB_EN, DB_TR, DB_CE, DB_RR, DB_CI, DB_AI, DB_RC};
  TFILE *f0, *f1[7];
  int nlevs0, nlevs1, nlast, rdn, utr;
  double er0[N_ELEMENTS1], er1[N_ELEMENTS1], de[N_ELEMENTS1];

  tt0 = WallTime();
  z = 0;
  for (k = 0; k < nk-1; k++) {
    k0 = ks[k];
    k1 = ks[k+1];
    wt0 = WallTime();
    sprintf(ifn, "%s%02d%02db.en", pref, k0, k1);
    printf("check levels: %s ... ", ifn);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 == NULL) {
      printf("cannot open file: %s\n", ifn);
      return -1;
    }
    if (z == 0) {
      z = (int)(fh.atom);
      strcpy(a, fh.symbol);
    } else if (z != (int)(fh.atom)) {
      printf("atomic number does not match: %d %d\n", z, (int)(fh.atom));
      return -1;
    }
    if (fh.nthreads > nth) nth = fh.nthreads;
    tlevs = 0;
    nlevs0 = 0;
    nlevs1 = 0;
    j = -1;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadENHeader(f0, &h0, swp);
      if (n == 0) break;
      tlevs += h0.nlevels;
      nt = 0;
      if (h0.nele == k0) {
	nlevs0 += h0.nlevels;
	if (h0.nele != j) {
	  nt = ReadENRecord(f0, &r0, swp);
	  er0[k] = r0.energy;
	}
      } else if (h0.nele == k1) {
	nlevs1 += h0.nlevels;
	if (h0.nele != j) {
	  nt = ReadENRecord(f0, &r0, swp);
	  er1[k] = r0.energy;
	}
      }
      j = h0.nele;
      FSEEK(f0, h0.length-nt, SEEK_CUR);
    }
    if (k > 0) {
      if (nlevs0 != nlast) {
	printf("nlevels do not match: %d %d %d\n", k0, nlast, nlevs0);
	return -1;
      }
    }
    nlast = nlevs1;
    wt1 = WallTime();
    printf("%d %.3e\n", tlevs, wt1-wt0);
    FCLOSE(f0);
  }

  de[0] = 0.0;
  for (k = 1; k < nk-1; k++) {
    de[k] = de[k-1] + er1[k-1]-er0[k];
  }
  for (i = 0; i < 7; i++) {
    sprintf(ofn, "%s%02d%02db.%s", pref, ks[0], ks[nk-1], ext[i]);
    fh1[i].atom = z;
    strcpy(fh1[i].symbol, a);
    fh1[i].type = types[i];
    f1[i] = OpenFileWTN(ofn, &fh1[i], nth);
  }

  FILE *frp0, *frp1;
  sprintf(ofn, "%s%02d%02db.rp", pref, ks[0], ks[nk-1]);
  frp1 = fopen(ofn, "w");
  if (frp1) {
    for (k = nk-1; k > 0; k--) {
      k0 = ks[k-1];
      k1 = ks[k];
      sprintf(ifn, "%s%02d%02db.rp", pref, k0, k1);
      frp0 = fopen(ifn, "r");
      if (frp0) {
	while (NULL != fgets(buf, 1024, frp0)) {
	  n = sscanf(buf, "%d %d", &i, &j);
	  if (k > 1 && j == k0) break;
	  fprintf(frp1, "%s", buf);
	}
	fclose(frp0);
      }
    }
    fclose(frp1);
  }

  nlevs0 = 0;
  for (k = nk-1; k > 0; k--) {
    k0 = ks[k-1];
    k1 = ks[k];    
    sprintf(ifn, "%s%02d%0db.en", pref, k0, k1);
    wt0 = WallTime();
    printf("write EN: %s ... ", ifn);
    f0 = OpenFileRO(ifn, &fh, &swp);
    rdn = 0;
    tlevs = 0;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadENHeader(f0, &h0, swp);
      if (n == 0) break;
      if (k == 1 || h0.nele > k0) {
	InitFile(f1[0], &fh1[0], &h0);
      }
      for (i = 0; i < h0.nlevels; i++) {
	n = ReadENRecord(f0, &r0, swp);
        if (h0.nele == k0) {
	  if (k > 1) {
	    rdn = 1;
	    break;
	  }
	}
	r0.energy += de[k-1];
	tlevs++;
	r0.ilev += nlevs0;
	WriteENRecord(f1[0], &r0);
      }
      if (rdn) break;
      DeinitFile(f1[0], &fh1[0]);
    }
    FCLOSE(f0);
    wt1 = WallTime();
    printf("%d %.3e\n", tlevs, wt1-wt0);
    
    sprintf(ifn, "%s%02d%02db.tr", pref, k0, k1);
    wt0 = WallTime();
    printf("write TR: %s ... ", ifn);
    f0 = OpenFileRO(ifn, &fh, &swp);
    nt = 0;    
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadTRHeader(f0, &h1, swp, &utr);
      if (n == 0) break;
      if (k > 1 && h1.nele == k0) {
	break;
      }
      InitFile(f1[1], &fh1[1], &h1);
      for (i = 0; i < h1.ntransitions; i++) {
	n = ReadTRRecord(f0, &r1, &r1x, swp, utr);
	if (n == 0) break;
	r1.lower += nlevs0;
	r1.upper += nlevs0;
	WriteTRRecord(f1[1], &r1, &r1x, utr);
	nt++;
      }
      DeinitFile(f1[1], &fh1[1]);
    }
    FCLOSE(f0);
    wt1 = WallTime();
    printf("%d %.3e\n", nt, wt1-wt0);

    sprintf(ifn, "%s%02d%02db.ce", pref, k0, k1);
    wt0 = WallTime();
    printf("write CE: %s ... ", ifn);
    f0 = OpenFileRO(ifn, &fh, &swp);
    nt = 0;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadCEHeader(f0, &h2, swp);
      if (n == 0) break;
      if (k > 1 && h2.nele == k0) {
	break;
      }
      InitFile(f1[2], &fh1[2], &h2);
      for (i = 0; i < h2.ntransitions; i++) {
	n = ReadCERecord(f0, &r2, swp, &h2);
	if (n == 0) break;
	r2.lower += nlevs0;
	r2.upper += nlevs0;
	WriteCERecord(f1[2], &r2);
	nt++;
	if (h2.qk_mode == QK_FIT) free(r2.params);
	free(r2.strength);
      }
      DeinitFile(f1[2], &fh1[2]);
      free(h2.tegrid);
      free(h2.egrid);
      free(h2.usr_egrid);
    }
    FCLOSE(f0);
    wt1 = WallTime();
    printf("%d %.3e\n", nt, wt1-wt0);

    sprintf(ifn, "%s%02d%02db.rr", pref, k0, k1);
    f0 = OpenFileRO(ifn, &fh, &swp);
    wt0 = WallTime();
    printf("write RR: %s ... ", ifn);
    nt = 0;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadRRHeader(f0, &h3, swp);
      if (n == 0) break;
      if (k > 1 && h3.nele == k0) {
	break;
      }
      InitFile(f1[3], &fh1[3], &h3);
      for (i = 0; i < h3.ntransitions; i++) {
	n = ReadRRRecord(f0, &r3, swp, &h3);
	if (n == 0) break;
	r3.b += nlevs0;
	if (r3.f >= 0) {
	  r3.f += nlevs0;
	}
	WriteRRRecord(f1[3], &r3);
	nt++;
	free(r3.params);
	free(r3.strength);
      }
      DeinitFile(f1[3], &fh1[3]);
      free(h3.tegrid);
      free(h3.egrid);
      free(h3.usr_egrid);
    }
    FCLOSE(f0);
    wt1 = WallTime();
    printf("%d %.3e\n", nt, wt1-wt0);

    sprintf(ifn, "%s%02d%02db.ci", pref, k0, k1);
    f0 = OpenFileRO(ifn, &fh, &swp);
    wt0 = WallTime();
    printf("write CI: %s ... ", ifn);
    nt = 0;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadCIHeader(f0, &h4, swp);
      if (k > 1 && h4.nele == k0) {
	break;
      }
      InitFile(f1[4], &fh1[4], &h4);
      for (i = 0; i < h4.ntransitions; i++) {
	n = ReadCIRecord(f0, &r4, swp, &h4);
	if (n == 0) break;
	r4.b += nlevs0;
	r4.f += nlevs0;
	WriteCIRecord(f1[4], &r4);
	nt++;
	free(r4.params);
	free(r4.strength);
      }
      DeinitFile(f1[4], &fh1[4]);
      free(h4.tegrid);
      free(h4.egrid);
      free(h4.usr_egrid);
    }
    FCLOSE(f0);
    wt1 = WallTime();
    printf("%d %.3e\n", nt, wt1);

    sprintf(ifn, "%s%02d%02db.ai", pref, k0, k1);
    f0 = OpenFileRO(ifn, &fh, &swp);
    wt0 = WallTime();
    printf("write AI: %s ... ", ifn);
    nt = 0;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadAIHeader(f0, &h5, swp);
      if (k > 1 && h5.nele == k0) {
	break;
      }
      InitFile(f1[5], &fh1[5], &h5);
      for (i = 0; i < h5.ntransitions; i++) {
	n = ReadAIRecord(f0, &r5, swp);
	if (n == 0) break;
	r5.b += nlevs0;
	r5.f += nlevs0;
	WriteAIRecord(f1[5], &r5);
	nt++;
      }
      DeinitFile(f1[5], &fh1[5]);
      free(h5.egrid);
    }
    FCLOSE(f0);
    wt1 = WallTime();
    printf("%d %.3e\n", nt, wt1-wt0);

    sprintf(ifn, "%s%02d%02db.rc", pref, k0, k1);
    f0 = OpenFileRO(ifn, &fh, &swp);
    wt0 = WallTime();
    printf("write RC: %s ... ", ifn);
    nt = 0;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadRCHeader(f0, &h6, swp);
      if (k > 1 && h6.nele == k0) {
	break;
      }
      InitFile(f1[6], &fh1[6], &h6);
      for (i = 0; i < h6.ntransitions; i++) {
	n = ReadRCRecord(f0, &r6, swp, &h6);
	if (n == 0) break;
	r6.lower += nlevs0;
	r6.upper += nlevs0;
	WriteRCRecord(f1[6], &r6);
	nt++;
	free(r6.rc);
      }
      DeinitFile(f1[6], &fh1[6]);
    }
    FCLOSE(f0);
    wt1 = WallTime();
    printf("%d %.3e\n", nt, wt1-wt0);

    nlevs0 += tlevs;
  }

  for (i = 0; i < 7; i++) {
    CloseFile(f1[i], &fh1[i]);
  }
  
  if (ic > 0) {
    k0 = ks[0];
    k1 = ks[nk-1];
    n = ic;
    k = 0;    
    for (i = 0; i < 7; i++) {
      if (n < 10) {	
	k = i < n;
      } else {
	k = ic%10;
	ic = ic/10;
      }
      if (k) {
	wt0 = WallTime();	
	sprintf(ifn, "%s%02d%02db.%s", pref, k0, k1, ext[i]);
	sprintf(ofn, "%s%02d%02da.%s", pref, k0, k1, ext[i]);
	printf("print table: %s %s ...", ifn, ofn);
	fflush(stdout);
	PrintTable(ifn, ofn, 1);
	wt1 = WallTime();
	printf(" %.3e\n", wt1-wt0);
      }
    }
  }
  tt1 = WallTime();
  printf("total time: %.3e\n", tt1-tt0);

  return 0;
}

void CombineDBase(char *pref, int k0, int k1, int kic, int nexc0, int ic) {
  int k, i, j, t, n, nb, nlevs, clevs, ni0, ni1, vn, vni, vn0, z;
  char ifn[1024], ofn[1024], buf[1024];
  char a[8], nc0[10*LNCOMPLEX], sn0[10*LSNAME], nm0[10*LNAME];
  F_HEADER fh, fh1[7];
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
  RC_HEADER h6;
  RC_RECORD r6;
  char ext[7][3] = {"en", "tr", "ce", "rr", "ci", "ai", "rc"};
  int types[7] = {DB_EN, DB_TR, DB_CE, DB_RR, DB_CI, DB_AI, DB_RC};
  TFILE *f0, *fr, *f1[7], *f0u, *f1u;
  int swp, *im, *imp, **ima, nim, nk, nilevs, *nplevs, nth;
  double e0, e1, e0p, e1p, tde, *de, *ei, **eai;
  float fde;
  int ilow2ph[3], iup2ph[3];
  double elow2ph[3], eup2ph[3];
  int nexc, ncap, nt, nd, ix, ibx, ifx, ib, cn0, cn1, cn2;
  int ia, *nm, *nklevs, *igk, *ifk, kk, kk1;
  short ***nc, nqc[N_ELEMENTS1], *muta, iuta0;
  double zh, t0, dt, d0, dd, *egk, **eo, eip, wt0, wt1, tt0, tt1;
  float ***pai;
  int nt1, nt2, nt3, nt4, nt5, nt6, nt3a, nt3b;
  int nbk, kbk, ibk, mbk, sbk, nq, i0, i1, i2, i3, utr;
  unsigned char *rrq;

  tt0 = WallTime();
  ncap = 0;
  if (nexc0 <= 0) nexc0 = 99;
  nexc = nexc0;  
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
  if (nk <= 0) {
    printf("invalid k0 and k1: %d %d\n", k0, k1);
    return;
  }
  ima = malloc(sizeof(int *)*nk);
  for (i = 0; i < nk; i++) ima[i] = NULL;
  de = malloc(sizeof(double)*nk);
  ei = malloc(sizeof(double)*nk);
  nplevs = malloc(sizeof(int)*nk);
  FILE *frp0, *frp1;
  sprintf(ofn, "%s%02d%02db.rp", pref, k0, k1);
  frp1 = fopen(ofn, "w");
  mbk = 10;
  sbk = mbk*(mbk+1);
  eo = malloc(sizeof(double *)*nk);
  muta = malloc(sizeof(short)*nk);
  for (k = k0; k <= k1; k++) {
    eo[k-k0] = malloc(sizeof(double)*sbk);  
    for (i = 0; i < sbk; i++) eo[k-k0][i] = 0;
    muta[k-k0] = 0;
  }
  if (frp1 != NULL) {
    for (k = k1; k >= k0; k--) {
      sprintf(ifn, "%s%02db.rp", pref, k);
      frp0 = fopen(ifn, "r");      
      if (frp0 != NULL) {
	while(NULL != fgets(buf, 1024, frp0)) {
	  fprintf(frp1, "%s", buf);
	  ibk = sscanf(buf, "%d %d %d %d %d %s %lg",
		       &i0, &i1, &nbk, &kbk, &j, nm0, &tde);
	  ibk = 2*((nbk-1)*nbk/2 + kbk);
	  if (j < 0) ibk++;
	  if (ibk < sbk) {
	    eo[k-k0][ibk] = -tde;
	  }
	  //printf("eo: %d %d %d %d %s %d %d %g\n",
	  //	 k, nbk, kbk, j, nm0, ibk, sbk, tde);
	}
	fclose(frp0);
      }      
    }
    fclose(frp1);
  }
  sprintf(ofn, "%s%02d%02db.uf", pref, k0, k1);
  frp1 = fopen(ofn, "w");
  if (frp1 != NULL) {
    for (k = k1; k >= k0; k--) {
      sprintf(ifn, "%s%02db.uf", pref, k);
      frp0 = fopen(ifn, "r");
      if (frp0 != NULL) {
	while(NULL != fgets(buf, 1024, frp0)) {
	  fprintf(frp1, "%s", buf);
	}
	fclose(frp0);
	fprintf(frp1, "\n");
      }
    }
    fclose(frp1);
  }
  nth = 0;
  z = 0;
  for (k = k1; k >= k0; k--) {
    sprintf(ifn, "%s%02db.en", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 == NULL) continue;
    if (z == 0) {
      z = (int)(fh.atom);
      strcpy(a, fh.symbol);
    }
    if (fh.nthreads > nth) nth = fh.nthreads;
    muta[k-k0] = iuta;
    FCLOSE(f0);
  }
  iuta = 0;
  for (k = k1; k >= k0; k--) {
    if (muta[k-k0]) {
      iuta = muta[k-k0];
      break;
    }
  }
  iuta0 = iuta;
  for (i = 0; i < 7; i++) {
    sprintf(ofn, "%s%02d%02db.%s", pref, k0, k1, ext[i]);
    fh1[i].atom = z;
    strcpy(fh1[i].symbol, a);
    fh1[i].type = types[i];
    f1[i] = OpenFileWTN(ofn, &fh1[i], nth);
  }
  clevs = 0;
  e0 = 0.0;
  e1 = 0.0;
  nm = malloc(sizeof(int)*nk);
  nklevs = malloc(sizeof(int)*nk);
  igk = malloc(sizeof(int)*nk);
  ifk = malloc(sizeof(int)*nk);
  egk = malloc(sizeof(double)*nk);
  for (k = k0; k <= k1; k++) {
    de[k-k0] = 0.0;
    nm[k-k0] = 0;
    nklevs[k-k0] = 0;
    igk[k-k0] = -1;
    ifk[k-k0] = -1;
    egk[k-k0] = 0.0;
  }
  for (k = k1; k >= k0; k--) {
    wt0 = WallTime();
    sprintf(ifn, "%s%02db.en", pref, k);
    printf("check levels: %d %d %s ...", z, k, ifn);
    fflush(stdout);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 == NULL) {
      printf("cannot open file %s\n", ifn);
      continue;
    }
    if (fh.type != DB_EN) {
      printf("%s is not of type DB_EN\n", ifn);
      continue;
    }
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
    if (_adj_ip) {
      eip = GetGroundIP(z, k, 0)/HARTREE_EV;
      if (eip > 0) {
	e1 = e0 + eip;
      } else {
	eip = -1;
      }
    } else {
      eip = -1;
    }
    if (k < k1) {
      de[k-k0] = de[k+1-k0]+ (e1p - e0);
    }
    ei[k-k0] = e1;
    im = malloc(sizeof(int)*nlevs);
    ima[k-k0] = im;
    nplevs[k-k0] = nlevs;
    FCLOSE(f0);
    for (i = 0; i < nlevs; i++) im[i] = 0;
    sprintf(ifn, "%s%02db.rr", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_RR) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadRRHeader(f0, &h3, swp);
	if (n == 0) break;
	for (i = 0; i < h3.ntransitions; i++) {
	  n = ReadRRRecord(f0, &r3, swp, &h3);
	  if (n == 0) break;
	  im[r3.b]++;
	  free(r3.params);
	  free(r3.strength);
	}
      }
    }
    if (f0) FCLOSE(f0);
    sprintf(ifn, "%s%02db.ai", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_AI) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadAIHeader(f0, &h5, swp);
	if (n == 0) break;
	for (i = 0; i < h5.ntransitions; i++) {
	  n = ReadAIRecord(f0, &r5, swp);
	  if (n == 0) break;
	  im[r5.b]++;
	}
      }
    }
    if (f0) FCLOSE(f0);
    wt1 = WallTime();
    printf(" %d %.3e\n", nlevs, wt1-wt0);
  }
  if (kic < k0) kic = k0;
  if (kic > k1) kic = k1;
  tde = de[kic-k0];
  for (k = k1; k >= k0; k--) {
    if (k <= kic) de[k-k0] = 0;
    else de[k-k0] -= tde;
  }
  for (k = k1; k >= k0; k--) {
    sprintf(ifn, "%s%02db.rc", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 == NULL) {
      printf("cannot open file %s\n", ifn);
    } else {
      if (fh.type != DB_RC) {
	printf("%s is not of type DB_RC\n", ifn);
	f0 = NULL;
      }
    }
    if (f0) {
      n = ReadRCHeader(f0, &h6, swp);
      if (n > 0) {
	ncap = h6.ncap;
	nexc = Min(nexc0, h6.nexc);
      } else {
	ncap = 0;
	nexc = nexc0;
      }
      FCLOSE(f0);
    } else {
      ncap = 0;
      nexc = nexc0;
    }
    
    wt0 = WallTime();
    sprintf(ifn, "%s%02db.en", pref, k);
    printf("rebuild level indices %d %d %s %d %d ...",
	   z, k, ifn, nexc, ncap);
    fflush(stdout);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 == NULL) {
      printf("cannot open file %s\n", ifn);
      continue;
    }
    if (fh.type != DB_EN) {
      printf("%s is not of type DB_EN\n", ifn);
      continue;
    }
    e1 = ei[k-k0];
    nlevs = nplevs[k-k0];
    im = ima[k-k0];
    nilevs = 0;
    for (i = 0; i < nlevs; i++) {
      if (im[i] > 0) {
	im[i] = -1;
      } else {
	im[i] = -2;
      }
    }
    int nbs = 0;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadENHeader(f0, &h0, swp);
      if (n == 0) break;
      if (h0.nele == k) {	
	for (i = 0; i < h0.nlevels; i++) {
	  n = ReadENRecord(f0, &r0, swp);
	  vn = abs(r0.p)/100;
	  if (nb == 0 && i == 0) vn0 = vn;
	  vni = VNIFromSName(r0.sname)/100;	      
	  if ((nexc > 0 && vn > nexc) ||
	      (ncap > 0 && vn > ncap && vni > vn0 && r0.energy > e1)) {
	    continue;
	  }
	  if (im[r0.ilev] == -2) {
	    continue;
	  }
	  if (_ncombex[k] > 0) {
	    for (ix = 0; ix < _ncombex[k]; ix++) {
	      if (strstr(r0.sname, _pcombex[k][ix])) {
		break;
	      }
	    }
	    if (ix < _ncombex[k]) continue;
	  }
	  if ((k == 2 && (_exk2h&1)) ||
	      (k > 2 && k <= 10 && (_exk2h&2)) ||
	      (k > 10 && k <= 28 && (_exk2h&4)) ||
	      (k > 28 && (_exk2h&8))) {
	    if (0 == FindNRShells(k, &r0, vn+1, nqc, NULL)) {
	      if (nqc[0] == 0) continue;
	    }
	  }
	  im[r0.ilev] = clevs;
	  clevs++;
	}
      } else {
	if (h0.nele == k-1 && k > k0) {
	  ni0 = h0.nlevels;
	  ri0 = malloc(sizeof(EN_RECORD)*ni0);
	  for (i = 0; i < h0.nlevels; i++) {
	    n = ReadENRecord(f0, &ri0[i], swp);
	    if (n == 0) break;
	  }
	  sprintf(ifn, "%s%02db.en", pref, k-1);
	  nim = FindLevelBlock(ni0, ri0, &ri1, k-1, ifn, &nbs);
	  if (nim != ni0) {
	    printf("\nionized mismatch: %d %d %d %d %d ... ",
		   k, nb, ni0, nim, nbs);
	    Abort(1);
	  }
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
      for (i = 0; i < nplevs[k+1-k0]; i++) {
	if (imp[i] <= -10) {
	  imp[i] = im[-(10+imp[i])];
	}
      }
    }
    if (k == k0) {
      FSEEK(f0, SIZE_F_HEADER, SEEK_SET);
      nilevs = 0;
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadENHeader(f0, &h0, swp);
	if (n == 0) break;
	if (h0.nele == k-1) {	
	  for (i = 0; i < h0.nlevels; i++) {
	    n = ReadENRecord(f0, &r0, swp);
	    if (nilevs >= _nilast) continue;
	    vn = abs(r0.p)/100;
	    if (nb == 0 && i == 0) vn0 = vn;
	    vni = VNIFromSName(r0.sname)/100;	
	    t = k-1;
	    if (_ncombex[t] > 0) {
	      for (ix = 0; ix < _ncombex[t]; ix++) {
		if (strstr(r0.sname, _pcombex[t][ix])) {
		  break;
		}
	      }
	      if (ix < _ncombex[t]) continue;
	    }
	    if ((t == 2 && (_exk2h&1)) ||
		(t > 2 && t <= 10 && (_exk2h&2)) ||
		(t > 10 && t <= 28 && (_exk2h&4)) ||
		(t > 28 && (_exk2h&8))) {
	      if (0 == FindNRShells(t, &r0, vn+1, nqc, NULL)) {
		if (nqc[0] == 0) continue;
	      }
	    }	        
	    if ((nexc > 0 && vn > nexc) ||
		(ncap > 0 && vn > ncap && vni > vn0 && r0.energy > e1)) {
	      continue;
	    }  
	    im[r0.ilev] = clevs;
	    nilevs++;
	    clevs++;
	  }
	} else {
	  FSEEK(f0, h0.length, SEEK_CUR);
	}
      }
    }
    FCLOSE(f0);
    wt1 = WallTime();
    printf(" %d %.3e\n", clevs, wt1-wt0);
  }
  rrq = malloc(sbk*clevs);
  for (i = 0; i < sbk*clevs; i++) rrq[i] = 0;
  for (k = k1; k >= k0; k--) {
    sprintf(ifn, "%s%02db.en", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 == NULL) {
      continue;
    }
    if (fh.type != DB_EN) {
      continue;
    }
    im = ima[k-k0];
    wt0 = WallTime();
    printf("writing level file: %d %d %s ...", z, k, ifn);
    fflush(stdout);
    e0 = 1e30;
    for (nb = 0; nb < fh.nblocks; nb++) {
      iuta = muta[k-k0];
      n = ReadENHeader(f0, &h0, swp);
      if (n == 0) break;
      if (h0.nele == k) {
	iuta = iuta0;
	InitFile(f1[0], &fh1[0], &h0);	
	for (i = 0; i < h0.nlevels; i++) {
	  iuta = muta[k-k0];
	  n = ReadENRecord(f0, &r0, swp);
	  if (im[r0.ilev] < 0) continue;
	  vn = abs(r0.p)/100;
	  if (!iuta && r0.ibase >= 0) {
	    r0.ibase = im[r0.ibase];
	  }
	  r0.ilev = im[r0.ilev];
	  r0.energy += de[k-k0];
	  if (r0.energy < e0) e0 = r0.energy;
	  Match2PhotonLevels(k, &r0, ilow2ph, iup2ph, elow2ph, eup2ph);
	  iuta = iuta0;
	  WriteENRecord(f1[0], &r0);
	  if (nm[k-k0] < vn) nm[k-k0] = vn;
	  if (ifk[k-k0] < 0) {
	    ifk[k-k0] = r0.ilev;
	    igk[k-k0] = r0.ilev;
	    egk[k-k0] = r0.energy;
	  } else {
	    if (r0.energy < egk[k-k0]) {
	      igk[k-k0] = r0.ilev;
	      egk[k-k0] = r0.energy;
	    }
	  }
	  nklevs[k-k0]++;
	  if (k > 1) {
	    if (0 == FillClosedShell(k, &r0, nc0, sn0, nm0)) {
	      n = StrSplit(nm0, '.');
	      i0 = 0;
	      for (t = 0; t < n; t++) {
		i1 = i0;
		while(isdigit(nm0[i1])) i1++;
		i2 = i1;
		while(nm0[i2] != '-' && nm0[i2] != '+') i2++;
		i2++;
		i3 = i2;
		while(isdigit(nm0[i3])) i3++;
		char c3 = nm0[i3];
		nm0[i3] = '\0';
		nq = atoi(nm0+i2);
		nm0[i2] = '\0';		
		ibk = GetJLFromSymbol(nm0+i1, &j, &kbk);
		if (ibk >= 0) {
		  nm0[i1] = '\0';
		  nbk = atoi(nm0+i0);
		  ibk = 2*((nbk-1)*nbk/2+kbk)+(j>0);
		  if (ibk < sbk) {
		    rrq[sbk*r0.ilev+ibk] = (unsigned char) nq;
		  }
		  //printf("rq: %d %d %d, %d %d %d %d, %d %d %d %d\n", r0.ilev, t, ibk, nbk, kbk, j, nq, i0, i1, i2, i3);
		}
		if (c3) {
		  i3++;
		  while(nm0[i3]) i3++;
		}
		i0 = i3 + 1;
	      }
	    }
	  }
	}
	iuta = iuta0;
	DeinitFile(f1[0], &fh1[0]);
      } else {
	FSEEK(f0, h0.length, SEEK_CUR);
      }      
    }

    if (k == k0) {
      FSEEK(f0, SIZE_F_HEADER, SEEK_SET);
      for (nb = 0; nb < fh.nblocks; nb++) {
	iuta = muta[k-k0];
	n = ReadENHeader(f0, &h0, swp);
	if (n == 0) break;
	if (h0.nele == k-1) {
	  iuta = iuta0;
	  InitFile(f1[0], &fh1[0], &h0);
	  for (i = 0; i < h0.nlevels; i++) {
	    iuta = muta[k-k0];
	    n = ReadENRecord(f0, &r0, swp);
	    if (im[r0.ilev] < 0) continue;
	    if (!iuta && r0.ibase >= 0) {
	      r0.ibase = im[r0.ibase];
	    }
	    r0.ilev = im[r0.ilev];
	    if (eip > 0) {
	      r0.energy = e0 + eip;
	    }
	    iuta = iuta0;
	    WriteENRecord(f1[0], &r0);
	  }
	  iuta = iuta0;
	  DeinitFile(f1[0], &fh1[0]);
	} else {
	  FSEEK(f0, h0.length, SEEK_CUR);
	}
      }
    }
    FCLOSE(f0);
    wt1 = WallTime();
    printf(" %d %.3e\n", nklevs[k-k0], wt1-wt0);
  }

  nc = NULL;
  pai = NULL;
  eai = NULL;
  if (_nc_iai > 0) {
    nc = malloc(sizeof(short **)*nk);
    for (k = 0; k < nk; k++) {
      nc[k] = malloc(sizeof(short *)*nklevs[k]);
      for (i = 0; i < nklevs[k]; i++) {
	nc[k][i] = malloc(sizeof(short)*nm[k]);
	for (j = 0; j < nm[k]; j++) {
	  nc[k][i][j] = 0;
	}
      }
    }
    pai = malloc(sizeof(float **)*nk);
    for (k = 0; k < nk; k++) {
      pai[k] = malloc(sizeof(float *)*nklevs[k]);
      for (i = 0; i < nklevs[k]; i++) {
	pai[k][i] = malloc(sizeof(float)*nm[k]*_nc_iai);
	for (j = 0; j < _nc_iai*nm[k]; j++) pai[k][i][j] = 0.0;
      }
    }
    eai = malloc(sizeof(double *)*nk);
    for (k = 0; k < nk; k++) {
      eai[k] = malloc(sizeof(double)*nm[k]*_nc_iai);
      for (i = 0; i < _nc_iai*nm[k]; i++) {
	eai[k][i] = 0.0;
      }
    }
  }
  CloseFile(f1[0], &fh1[0]);
  wt0 = WallTime();
  sprintf(ifn, "%s%02d%02db.en", pref, k0, k1);
  printf("build EN mem: %d %d %s ...", k0, k1, ifn);
  fflush(stdout);
  MemENTableWC(ifn, k0, ifk, nc);
  int mesize = mem_en_table_size;
  EN_SRECORD *metable = malloc(sizeof(EN_SRECORD)*mesize);
  memcpy(metable, mem_en_table, sizeof(EN_SRECORD)*mesize);
  wt1 = WallTime();
  printf(" %d %.3e\n", mesize, wt1-wt0);
  for (k = k1; k >= k0; k--) {
    wt0 = WallTime();
    printf("processing %d %d ... ", z, k);
    fflush(stdout);
    im = ima[k-k0];
    iuta = muta[k-k0];
    sprintf(ifn, "%s%02db.en", pref, k);
    MemENTable(ifn);
    sprintf(ifn, "%s%02db.tr", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    nt1 = 0;
    nt2 = 0;
    nt3 = 0;
    nt4 = 0;
    nt5 = 0;
    nt6 = 0;
    nt3a = 0;
    nt3b = 0;
    if (f0 != NULL && fh.type == DB_TR) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	iuta = muta[k-k0];
	n = ReadTRHeader(f0, &h1, swp, &utr);
	if (n == 0) break;
	iuta = iuta0;
	InitFile(f1[1], &fh1[1], &h1);	
	for (i = 0; i < h1.ntransitions; i++) {
	  iuta = muta[k-k0];
	  n = ReadTRRecord(f0, &r1, &r1x, swp, utr);
	  if (n == 0) break;
	  if (im[r1.lower] >= 0 && im[r1.upper] >= 0) {
	    r1.lower = im[r1.lower];
	    r1.upper = im[r1.upper];
	    iuta = iuta0;
	    WriteTRRecord(f1[1], &r1, &r1x, utr);
	    nt1++;
	  }
	}
	iuta = iuta0;
	DeinitFile(f1[1], &fh1[1]);
      }
      double r2p = 0.0;
      i = -1;
      if (k == 1 && ilow2ph[0] >= 0 && iup2ph[0] >= 0) {
	i = 0;
      } else if (k == 2 && ilow2ph[1] >= 0 && iup2ph[1] >= 0) {
	i = 1;
      } else if (k == 4 && ilow2ph[2] >= 0 && iup2ph[2] >= 0) {
	i = 2;
      }
      if (i >= 0) {
	r2p = TwoPhotonRate(z, i);
	if (i == 0) r2p *= 2;
	tde = eup2ph[i]-elow2ph[i];	
	tde = pow(tde*FINE_STRUCTURE_CONST,2);
	r2p /= 2*tde*FINE_STRUCTURE_CONST*RATE_AU;
	h1.ntransitions = 1;
	h1.gauge = G_TWOPHOTON;
	h1.multipole = 0;
	iuta = iuta0;
	InitFile(f1[1], &fh1[1], &h1);
	r1.lower = ilow2ph[i];
	r1.upper = iup2ph[i];
	r1.strength = (float)r2p;
	WriteTRRecord(f1[1], &r1, &r1x, iuta);
	nt1++;
      }
      iuta = iuta0;
      DeinitFile(f1[1], &fh1[1]);
      FCLOSE(f0);
    }

    sprintf(ifn, "%s%02db.ce", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_CE) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	iuta = muta[k-k0];
	n = ReadCEHeader(f0, &h2, swp);
	if (n == 0) break;
	iuta = iuta0;
	InitFile(f1[2], &fh1[2], &h2);
	for (i = 0; i < h2.ntransitions; i++) {
	  iuta = muta[k-k0];
	  n = ReadCERecord(f0, &r2, swp, &h2);
	  if (n == 0) break;
	  if (im[r2.lower] >= 0 && im[r2.upper] >= 0) {
	    r2.lower = im[r2.lower];
	    r2.upper = im[r2.upper];
	    iuta = iuta0;
	    WriteCERecord(f1[2], &r2);
	    nt2++;
	  }
	  if (h2.qk_mode == QK_FIT) free(r2.params);
	  free(r2.strength);
	}
	iuta = iuta0;
	DeinitFile(f1[2], &fh1[2]);
	free(h2.tegrid);
	free(h2.egrid);
	free(h2.usr_egrid);
      }
      FCLOSE(f0);
    }
    
    sprintf(ifn, "%s%02db.rr", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_RR) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	iuta = muta[k-k0];
	n = ReadRRHeader(f0, &h3, swp);
	if (n == 0) break;
	iuta = iuta0;
	InitFile(f1[3], &fh1[3], &h3);
	for (i = 0; i < h3.ntransitions; i++) {
	  iuta = muta[k-k0];
	  n = ReadRRRecord(f0, &r3, swp, &h3);
	  if (n == 0) break;
	  ibx = r3.b;
	  ifx = r3.f;
	  if (im[r3.b] >= 0) {
	    if (im[r3.f] >= 0) {
	      r3.b = im[r3.b];
	      r3.f = im[r3.f];
	      tde = metable[r3.f].energy - metable[r3.b].energy;
	      if (tde > 0) {
		iuta = iuta0;
		WriteRRRecord(f1[3], &r3);
		nbk = r3.kl/1000;
		kbk = r3.kl%1000;
		ibk = 2*(nbk*(nbk-1)/2 + kbk);
		if (ibk < sbk) {
		  i0 = r3.b*sbk + ibk;
		  rrq[i0] = 0;
		  rrq[i0+1] = 0;		  
		}
		nt3++;
	      }
	    } else {
	      tde = mem_en_table[r3.f].energy - mem_en_table[r3.b].energy;
	      if (tde > 0) {
		r3.b = im[r3.b];
		nbk = r3.kl/1000;
		kbk = r3.kl%1000;
		ibk = 2*(nbk*(nbk-1)/2 + kbk);
		if (ibk < sbk) {
		  i0 = r3.b*sbk + ibk;
		  rrq[i0] = 0;
		  rrq[i0+1] = 0;
		}
		r3.f = -r3.kl;
		fde = tde;
		r3.kl = *((int *) &fde);
		iuta = iuta0;
		WriteRRRecord(f1[3], &r3);
		nt3a++;
	      }
	    }
	  }
	  free(r3.params);
	  free(r3.strength);
	}
	iuta = iuta0;
	DeinitFile(f1[3], &fh1[3]);
	free(h3.tegrid);
	free(h3.egrid);
	free(h3.usr_egrid);
      }
      FCLOSE(f0);
    }
    
    // add type 10 pi trans
    sprintf(ifn, "%s%02db.uen", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    int *ium = NULL;
    int *uim = NULL;
    if (f0) {
      ium = malloc(sizeof(int)*sbk);
      for (i = 0; i < sbk; i++) {
	ium[i] = -1;
      }
      nlevs = 0;
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadENHeader(f0, &h0, swp);
	nlevs += h0.nlevels;
	if (n == 0) break;
	if (h0.nele == 0) break;
	for (i = 0; i < h0.nlevels; i++) {
	  n = ReadENRecord(f0, &r0, swp);	  
	  if (h0.nele == 1) {
	    i0 = 0;
	    i1 = 0;
	    while(isdigit(r0.name[i1])) i1++;
	    i2 = i1;
	    while(r0.name[i2] != '-' && r0.name[i2] != '+') i2++;
	    i3 = i2;
	    while(isdigit(r0.name[i3])) i3++;
	    r0.name[i3] = '\0';
	    nq = atoi(r0.name+i2);
	    r0.name[i2] = '\0';
	    ibk = GetJLFromSymbol(r0.name+i1, &j, &kbk);
	    if (ibk >= 0) {
	      r0.name[i1] = '\0';
	      nbk = atoi(r0.name+i0);
	      ibk = 2*((nbk-1)*nbk/2 + kbk) + (j>0);
	      if (ibk < sbk) {
		ium[ibk] = r0.ilev;
		eo[k-k0][ibk] = -r0.energy;
	      }
	    }
	  }
	}	
      }
      FCLOSE(f0);
      uim = malloc(sizeof(int)*nlevs);
      for (i = 0; i < nlevs; i++) uim[i] = -1;
      for (ibk = 0; ibk < sbk; ibk++) {
	if (ium[ibk] >= 0) {
	  uim[ium[ibk]] = ibk;
	}
      }
      sprintf(ifn, "%s%02db.urr", pref, k);
      f0 = OpenFileRO(ifn, &fh, &swp);
      RR_RECORD *r3u = NULL;
      if (f0) {
	r3u = malloc(sizeof(RR_RECORD)*sbk);
	for (i = 0; i < sbk; i++) {
	  r3u[i].params = NULL;
	  r3u[i].strength = NULL;
	}
	float *tpar, *tstr;
	tpar = malloc(sizeof(float)*h3.nparams);
	tstr = malloc(sizeof(float)*h3.n_usr);
	int ibmin, ibmax;
	for (nb = 0; nb < fh.nblocks; nb++) {
	  n = ReadRRHeader(f0, &h3, swp);
	  if (n == 0) break;
	  ibmin = sbk;
	  ibmax = 0;
	  for (i = 0; i < h3.ntransitions; i++) {
	    n = ReadRRRecord(f0, &r3, swp, &h3);
	    ibk = uim[r3.b];
	    if (ibk >= 0) {
	      memcpy(&r3u[ibk], &r3, sizeof(RR_RECORD));
	      r3u[ibk].f = -r3.kl;
	      fde = eo[k-k0][ibk];
	      r3u[ibk].kl = *((int *) &fde);
	      r3u[ibk].params = malloc(sizeof(float)*h3.nparams);
	      r3u[ibk].strength = malloc(sizeof(float)*h3.n_usr);
	      memcpy(r3u[ibk].params, r3.params, sizeof(float)*h3.nparams);
	      memcpy(r3u[ibk].strength, r3.strength, sizeof(float)*h3.n_usr);
	      if (ibmin > ibk) ibmin = ibk;
	      if (ibmax < ibk) ibmax = ibk;
	    }
	    free(r3.params);
	    free(r3.strength);
	  }
	  h3.nele = k;
	  InitFile(f1[3], &fh1[3], &h3);
	  for (i = 0; i < nklevs[k-k0]; i++) {
	    i1 = i + ifk[k-k0];
	    for (nbk = 1; nbk <= mbk; nbk++) {
	      for (kbk = 0; kbk < nbk; kbk++) {
		ibk = 2*((nbk-1)*nbk/2 + kbk);
		if (ibk >= ibmin && ibk <= ibmax) {
		  if (kbk == 0) j = 1;
		  else j = 0;
		  for (; j <= 1; j++) {
		    i0 = i1*sbk + ibk + j;
		    if (rrq[i0] > 0) {		    
		      //e1 = (1.0/(kbk+j) - 0.75/nbk)*nbk;	      
		      //e0 = 2*eo[k-k0][ibk+j]*FINE_STRUCTURE_CONST2;
		      //zh = nbk*sqrt((sqrt(1.0+4*e1*e0) - 1.0)/(2*e1));
		      //zh /= FINE_STRUCTURE_CONST;
		      memcpy(&r3, &r3u[ibk], sizeof(RR_RECORD));
		      r3.b = i1;
		      r3.params = tpar;
		      r3.strength = tstr;
		      memcpy(r3.params, r3u[ibk].params, sizeof(float)*h3.nparams);
		      memcpy(r3.strength, r3u[ibk].strength, sizeof(float)*h3.n_usr);
		      tde = rrq[i0]*(metable[r3.b].j+1.0);
		      if (j == 0) tde /= (2*kbk-1.0);
		      else tde /= (2*kbk+1.0);
		      r3.params[0] *= tde;
		      for (t = 0; t < h3.n_usr; t++) {
			r3.strength[t] *= tde;
		      }
		      iuta = iuta0;
		      WriteRRRecord(f1[3], &r3);
		      nt3b++;
		    }
		  }	      
		}
	      }
	    }
	  }
	  iuta = iuta0;
	  DeinitFile(f1[3], &fh1[3]);
	  free(h3.tegrid);
	  free(h3.egrid);
	  free(h3.usr_egrid);
	}
	FCLOSE(f0);
	free(tpar);
	free(tstr);
	for (i = 0; i < sbk; i++) {
	  free(r3u[i].params);
	  free(r3u[i].strength);
	}
	free(r3u);
	free(ium);
	free(uim);
      }
    }
    
    sprintf(ifn, "%s%02db.ci", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_CI) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	iuta = muta[k-k0];
	n = ReadCIHeader(f0, &h4, swp);
	if (n == 0) break;
	iuta = iuta0;
	InitFile(f1[4], &fh1[4], &h4);
	for (i = 0; i < h4.ntransitions; i++) {
	  iuta = muta[k-k0];
	  n = ReadCIRecord(f0, &r4, swp, &h4);
	  if (n == 0) break;
	  if (im[r4.b] >= 0 && im[r4.f] >= 0) {
	    r4.b = im[r4.b];
	    r4.f = im[r4.f];
	    tde = metable[r4.f].energy-metable[r4.b].energy;
	    if (tde > 0) {
	      iuta = iuta0;
	      WriteCIRecord(f1[4], &r4);
	      nt4++;
	    }
	  }
	  free(r4.params);
	  free(r4.strength);
	}
	iuta = iuta0;
	DeinitFile(f1[4], &fh1[4]);
	free(h4.tegrid);
	free(h4.egrid);
	free(h4.usr_egrid);
      }
      FCLOSE(f0);
    }

    sprintf(ifn, "%s%02db.ai", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_AI) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	iuta = muta[k-k0];
	n = ReadAIHeader(f0, &h5, swp);
	if (n == 0) break;
	iuta = iuta0;
	InitFile(f1[5], &fh1[5], &h5);
	for (i = 0; i < h5.ntransitions; i++) {
	  iuta = muta[k-k0];
	  n = ReadAIRecord(f0, &r5, swp);
	  if (n == 0) break;
	  if (im[r5.b] >= 0 && im[r5.f] >= 0) {
	    r5.b = im[r5.b];
	    r5.f = im[r5.f];
	    tde = metable[r5.b].energy - metable[r5.f].energy;
	    if (tde > 0 && r5.rate > 0) {
	      iuta = iuta0;
	      WriteAIRecord(f1[5], &r5);
	      nt5++;
	      if (k > k0 && k > 2 && _nc_iai > 0) {
		ia = ChannelAI(r5.b, nm[k-k0], nc[k-k0][r5.b-ifk[k-k0]],
			       r5.f, nm[k-k0-1], nc[k-k0-1][r5.f-ifk[k-k0-1]],
			       &cn0, &cn1, &cn2);
		if (ia >= 0) {
		  pai[k-k0][r5.b-ifk[k-k0]][ia] += r5.rate;
		  eai[k-k0][ia] += r5.rate*tde;
		}
	      }
	    }
	  }
	}
	iuta = iuta0;
	DeinitFile(f1[5], &fh1[5]);
	free(h5.egrid);
      }
      FCLOSE(f0);
    }
    if (_nc_iai > 0) {
      for (ia = 0; ia < nm[k-k0]*_nc_iai; ia++) {
	dd = 0.0;
	for (i = 0; i < nklevs[k-k0]; i++) {
	  dd += pai[k-k0][i][ia];
	}
	if (dd > 0) eai[k-k0][ia] /= dd;
      }
    }
    sprintf(ifn, "%s%02db.rc", pref, k);
    f0 = OpenFileRO(ifn, &fh, &swp);
    if (f0 != NULL && fh.type == DB_RC) {
      for (nb = 0; nb < fh.nblocks; nb++) {
	iuta = muta[k-k0];
	n = ReadRCHeader(f0, &h6, swp);
	if (n == 0) break;
	if (h6.ntransitions == 0) continue;
	h6.mexc = nexc0;
	iuta = iuta0;
	InitFile(f1[6], &fh1[6], &h6);
	for (i = 0; i < h6.ntransitions; i++) {
	  iuta = muta[k-k0];
	  n = ReadRCRecord(f0, &r6, swp, &h6);
	  if (n == 0) break;
	  if (im[r6.lower] >= 0 && im[r6.upper] >= 0) {
	    r6.lower = im[r6.lower];
	    r6.upper = im[r6.upper];
	    tde = metable[r6.upper].energy - metable[r6.lower].energy;
	    if (tde > 0) {
	      iuta = iuta0;
	      WriteRCRecord(f1[6], &r6);
	      nt6++;
	    }
	  }
	  free(r6.rc);
	}
	iuta = iuta0;
	DeinitFile(f1[6], &fh1[6]);
      }
      FCLOSE(f0);
    }
    if (im) free(im);
    wt1 = WallTime();
    printf(" %d %d %d %d %d %d %d %d %.3e\n",
	   nt1, nt2, nt3, nt3a, nt3b, nt4, nt5, nt6, wt1-wt0);
  }

  if (_nc_iai > 0) {
    for (k = k0+1; k <= k1; k++) {
      if (k < 3) continue;
      kk1 = k-k0;
      kk = kk1-1;
      nb = Min(nm[kk], nm[kk1]);
      for (ia = 0; ia < _nc_iai*nb; ia++) {
	wt0 = WallTime();
	r5.f = -(_ic_iai[ia%_nc_iai]*1000 + (ia/_nc_iai) + 1);
	h5.nele = k;
	h5.n_egrid = 1;
	h5.egrid = &eai[kk][ia];
	iuta = iuta0;
	InitFile(f1[5], &fh1[5], &h5); 
	for (i = 0; i < nklevs[kk1]; i++) {
	  j = i + ifk[kk1];
	  ib = metable[j].ibase;
	  if (ib >= 0) {
	    d0 = pai[kk][ib-ifk[kk]][ia];
	    if (!pai[kk1][i][ia] && d0) {
	      pai[kk1][i][ia] = d0;
	      eai[kk1][ia] = eai[kk][ia];
	      r5.rate = d0;
	      r5.b = j;
	      iuta = iuta0;
	      WriteAIRecord(f1[5], &r5);
	    }
	  }
	}
	wt1 = WallTime();
	if (ai_header.ntransitions > 0) {
	  printf("inner shell ai: %3d %6d %5d %12.5E %.3e\n", k, -r5.f, ai_header.ntransitions, h5.egrid[0], wt1-wt0);
	}
	iuta = iuta0;
	DeinitFile(f1[5], &fh1[5]);      
      }
    }
  }
  free(ima);
  free(de);
  free(ei);
  free(nplevs);
  for (i = 1; i < 7; i++) {
    iuta = iuta0;
    CloseFile(f1[i], &fh1[i]);
  }
  if (_nc_iai > 0) {
    for (k = 0; k < nk; k++) {
      for (i = 0; i < nklevs[k]; i++) {
	free(nc[k][i]);
	free(pai[k][i]);
      }
      free(nc[k]);
      free(pai[k]);
      free(eai[k]);
    }
    free(nc);
    free(pai);
    free(eai);
  }
  free(nm);
  free(nklevs);
  free(igk);
  free(ifk);
  free(egk);
  free(metable);
  FreeMemENTable();
  for (k = k0; k <= k1; k++) {
    if (eo[k-k0]) free(eo[k-k0]);
  }
  if (eo) free(eo);
  free(muta);
  free(rrq);
  if (ic > 0) {
    n = ic;
    k = 0;    
    for (i = 0; i < 7; i++) {
      if (n < 10) {	
	k = i < n;
      } else {
	k = ic%10;
	ic = ic/10;
      }
      if (k) {
	wt0 = WallTime();	
	sprintf(ifn, "%s%02d%02db.%s", pref, k0, k1, ext[i]);
	sprintf(ofn, "%s%02d%02da.%s", pref, k0, k1, ext[i]);
	printf("print table: %s %s ...", ifn, ofn);
	fflush(stdout);
	PrintTable(ifn, ofn, 1);
	wt1 = WallTime();
	printf(" %.3e\n", wt1-wt0);
      }
    }
  }
  tt1 = WallTime();
  printf("total time: %.3e\n", tt1-tt0);
}

int GroupLevels(EN_RECORD *rs, int nr, double ei, double des,
		int minlev, LEVGRP *rg) {
  int i, ii, iu, ku, n, ng, j, k, j0;
  double e0, w0, w, wt;

  if (nr <= minlev) {
    ng = nr;
    if (rg) {
      for (i = 0; i < ng; i++) {
	memcpy(&rg[i].r, &rs[i], sizeof(EN_RECORD));
	if (rg[i].r.j >= 0) {
	  rg[i].r.ibase = rg[i].r.j;
	  rg[i].r.j = -2;
	}
	rg[i].nlev = 1;
	rg[i].ilev = malloc(sizeof(int));
	rg[i].ilev[0] = i;
      }
    }
    return ng;
  }

  ku = nr;
  for (i = 0; i < nr; i++) {
    if (rs[i].j == -1) {
      ku = i;
      break;
    }
  }
  
  ii = ku;
  for (i = 0; i < ku; i++) {
    if (rs[i].energy >= ei) {
      ii = i;
      break;
    }
  }

  iu = nr;
  for (i = ku; i < nr; i++) {
    if (rs[i].energy >= ei) {
      iu = i;
      break;
    }
  }

  if (ku < minlev ) ii = iu;
  ng = Min(minlev, ii);
  if (rg) {
    for (i = 0; i < ng; i++) {
      memcpy(&rg[i].r, &rs[i], sizeof(EN_RECORD));
      if (rg[i].r.j >= 0) {
	rg[i].r.ibase = rg[i].r.j;
	rg[i].r.j = -2;
      }
      rg[i].nlev = 1;
      rg[i].ilev = malloc(sizeof(int));
      rg[i].ilev[0] = i;
    }
  } else {
    i = ng;
  }
  
  e0 = rs[i++].energy;
  n = 1;
  for (; i <= nr; i++) {
    if (i == ii || i == ku || i == iu || i == nr || rs[i].energy-e0 > des) {
      if (rg && n > 0) {
	rg[ng].nlev = n;
	rg[ng].ilev = malloc(sizeof(int)*n);
	j0 = 0;
	w0 = 0.0;
	wt = 0.0;
	for (k = 0; k < n; k++) {
	  j = i-n+k;
	  rg[ng].ilev[k] = j;
	  w = WFromENRecord(&rs[j]);
	  if (w > w0) {
	    j0 = j;
	    w0 = w;
	  }
	  wt += w;
	}
	memcpy(&rg[ng].r, &(rs[j0]), sizeof(EN_RECORD));
	if (n == 1) {
	  if (rg[ng].r.j >= 0) {
	    rg[ng].r.j = -2;
	  }
	} else {
	  rg[ng].r.j = -1;
	  j = 0;
	  for (k = 0; k < LNAME; k++) {	    
	    if (rs[j0].name[k] == '(') {
	      k++;
	      while (rs[j0].name[k] != ')') k++;
	      while (rs[j0].name[k] != '.' && rs[j0].name[k] != '\0') k++;
	    }
	    rg[ng].r.name[j] = rs[j0].name[k];
	    j++;
	    if (rs[j0].name[k] == '\0') break;
	  }
	  for (; j < LNAME; j++) rg[ng].r.name[j] = '\0';
	}
	if (wt > 0x80000000) {
	  w = wt-1;
	  j = w/0x80000000;
	  w = j;
	  w *= 0x80000000;
	  w = wt-w;
	  if (j < 0 || j > 3276) j = 3276;
	  rg[ng].r.j = -(j*10 - rg[ng].r.j);
	  rg[ng].r.ibase = (int)(w+0.1);
	} else {
	  rg[ng].r.ibase = (int)(wt-0.5);
	}
	rg[ng].r.energy *= w0;
	for (k = 0; k < n; k++) {
	  j = i-n+k;
	  if (j == j0) continue;
	  w = WFromENRecord(&rs[j]);
	  rg[ng].r.energy += rs[j].energy*w;
	}
	rg[ng].r.energy /= wt;
      }
      n = 1;
      if (i < nr) e0 = rs[i].energy;
      ng++;
    } else {
      n++;
    }
  }

  return ng;
}

int LoadSFU(char *ipr, int ke, double **efu) {
  char buf[1024], ifn[1024];
  FILE *f;
  double d0, d1, d2, d3, d4;
  int nmx, i, j, k, t, n, i0, i1;
  int nm1, nm2, nm3, nm4, nm5;
  
  sprintf(ifn, "%sb.uf", ipr);
  f = fopen(ifn, "r");
  if (f == NULL) {
    printf("cannot open file: %s", ifn);
    return 0;
  }

  *efu = NULL;
  nmx = 0;
  while (1) {
    if (NULL == fgets(buf, 1024, f)) break;
    char *c = buf;
    while (c && isspace(*c)) c++;
    if (*c == '\0') continue;
    if (buf[0] == '#') {
      n = sscanf(buf+1, "%d %d", &k, &nmx);
      if (n != 2) {
	printf("invalid sfu file 0: %s\n", ifn);
	fclose(f);
	return 0;
      }
      if (k != ke) {
	for (i = 0; i < nmx; i++) {
	  for (j = i; j < nmx; j++) {
	    if (NULL == fgets(buf, 1024, f)) {
	      fclose(f);
	      return 0;
	    }
	  }
	}
	continue;
      }
      nm1 = nmx*nmx;
      nm2 = nm1*2;
      nm3 = nm1*3;
      nm4 = nm1*4;
      nm5 = nm1*5;
      *efu = malloc(sizeof(double)*6*nm1);
      for (i = 0; i < 6*nm1; i++) {
	(*efu)[i] = 0.0;
      }
      for (i = 0; i < nmx; i++) {
	for (j = 0; j < nmx; j++) {
	  char *r = fgets(buf, 1024, f);
	  if (NULL == r) {
	    free(*efu);
	    printf("invalid sfu file 1: %s\n", ifn);
	    return 0;
	  }
	  n=sscanf(buf, "%d %d %lg %lg %lg %lg %lg %d %d\n",
		   &t, &k, &d0, &d1, &d2, &d3, &d4, &i0, &i1);
	  k = j*nmx + i;
	  (*efu)[k] = d0/HARTREE_EV;
	  (*efu)[k+nm1] = d1/HARTREE_EV;
	  (*efu)[k+nm2] = d2;
	  (*efu)[k+nm3] = d3;
	  (*efu)[k+nm4] = i0;
	  (*efu)[k+nm5] = i1;
	}
      }
      break;
    }
  }
  fclose(f);
  if (*efu == NULL) nmx = 0;
  return nmx;
}

int LoadGroupMod(char *fn, int z0, GROUPMOD *gd) {
  FILE *f;
  char buf[1024], *p;

  printf("load groupmod data from %s for z=%d ... ", fn, z0);
  int i, j, n, zi, ki, k0, k1, n0, n1, z;
  double d0, d1;

  if (z0 > 0) {
    z = z0;
  } else {
    z = -z0;
    z0 = 0;
  }
  for (i = 0; i <= z; i++) {
    gd[i].minlevs = 0;
    gd[i].maxlevs = 0;
    gd[i].des0 = 0.0;
    gd[i].des1 = 0.0;
  }
  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file: %s\n", fn);
    return 0;
  }
  j = 0;
  while (1) {
    if (NULL == fgets(buf, 1024, f)) break;
    p = buf;
    while(isspace(*p)) p++;
    if (p[0] == '#') continue;
    n = sscanf(p, "%d %d %d %d %lg %lg",
	       &zi, &ki, &n0, &n1, &d0, &d1);
    if (n != 6) continue;
    if (z0 != zi) continue;
    if (ki <= 0) continue;
    if (ki > 1000) {
      k0 = ki/1000;
      k1 = ki%1000;
    } else {
      k0 = ki;
      k1 = ki;
    }
    for (ki = k0; ki <= k1; ki++) {
      if (ki < 0 || ki > z) continue;
      gd[ki].minlevs = n0;
      gd[ki].maxlevs = n1;
      gd[ki].des0 = d0/HARTREE_EV;
      gd[ki].des1 = d1/HARTREE_EV;
      j++;
    }
  }
  fclose(f);
  printf("%d\n", j);
  return j;
}

int LoadTransMod(char *fn, int z0, TRANSMOD *sd) {
  FILE *f;
  char buf[1024], *p;

  printf("load transmod data from %s for z=%d ... ", fn, z0);
  int i, j, n, ns, z, zi, ki, n0, n1, k0, k1;
  double e0, e1, d0, d1, f0, f1, w0;
  
  if (z0 > 0) {
    z = z0;
  } else {
    z = -z0;
    z0 = 0;
  }
  for (i = 0; i < z; i++) {
    sd[i].minlo = 1000000;
    sd[i].maxlo = 0;
    sd[i].minup = 1000000;
    sd[i].maxup = 0;
    sd[i].nde = 0;
    sd[i].emin = NULL;
    sd[i].emax = NULL;
    sd[i].fde = NULL;
    sd[i].ude = NULL;
    sd[i].fmf = NULL;
    sd[i].umf = NULL;
    sd[i].uwf = NULL;
    sd[i].fme = NULL;
    sd[i].ume = NULL;
    sd[i].fwe = NULL;
    sd[i].uwe = NULL;
    sd[i].fst = NULL;
    sd[i].ust = NULL;
    sd[i].usd = NULL;
  }

  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file: %s\n", fn);
    return 0;
  }
  
  while (1) {
    if (NULL == fgets(buf, 1024, f)) break;
    p = buf;
    while(isspace(*p)) p++;
    if (p[0] == '#') continue;
    n = sscanf(p, "%d %d %d %d %lg %lg %lg %lg %lg %lg %lg",
	       &zi, &ki, &n0, &n1, &e0, &e1, &d0, &d1, &f0, &f1, &w0);
    if (n != 11) continue;
    if (z0 != zi) continue;
    if (ki <= 0) continue;
    if (ki > 1000) {
      k0 = ki/1000;
      k1 = ki%1000;
    } else {
      k0 = ki;
      k1 = ki;
    }
    for (ki = k0; ki <= k1; ki++) {      
      if (ki < 1 || ki > z) continue;
      i = ki-1;
      if (n0 < sd[i].minlo) sd[i].minlo = n0;
      if (n0 > sd[i].maxlo) sd[i].maxlo = n0;
      if (n1 < sd[i].minup) sd[i].minup = n1;
      if (n1 > sd[i].maxup) sd[i].maxup = n1;
    }
  }
  ns = 0;
  for (i = 0; i < z; i++) {
    if (sd[i].maxlo >= sd[i].minlo &&
	sd[i].maxup >= sd[i].minup) {
      sd[i].nlo = 1 + sd[i].maxlo - sd[i].minlo;
      sd[i].nup = 1 + sd[i].maxup - sd[i].minup;
      sd[i].nde = sd[i].nlo * sd[i].nup;
      sd[i].emin = malloc(sizeof(double)*sd[i].nde);
      sd[i].emax = malloc(sizeof(double)*sd[i].nde);
      sd[i].fde = malloc(sizeof(double)*sd[i].nde);
      sd[i].ude = malloc(sizeof(double)*sd[i].nde);
      sd[i].fmf = malloc(sizeof(double)*sd[i].nde);
      sd[i].umf = malloc(sizeof(double)*sd[i].nde);
      sd[i].uwf = malloc(sizeof(double)*sd[i].nde);
      sd[i].fme = malloc(sizeof(double)*sd[i].nde);
      sd[i].fwe = malloc(sizeof(double)*sd[i].nde);
      sd[i].ume = malloc(sizeof(double)*sd[i].nde);
      sd[i].uwe = malloc(sizeof(double)*sd[i].nde);
      sd[i].fst = malloc(sizeof(double)*sd[i].nde);
      sd[i].ust = malloc(sizeof(double)*sd[i].nde);
      sd[i].usd = malloc(sizeof(double)*sd[i].nde);
      for (j = 0; j < sd[i].nde; j++) {
	sd[i].emin[j] = 0.0;
	sd[i].emax[j] = 0.0;
	sd[i].fde[j] = 0.0;
	sd[i].ude[j] = 0.0;
	sd[i].fmf[j] = 1.0;
	sd[i].umf[j] = 1.0;
	sd[i].uwf[j] = 1.0;
	sd[i].fme[j] = 0.0;
	sd[i].fwe[j] = 0.0;
	sd[i].ume[j] = 0.0;
	sd[i].uwe[j] = 0.0;
	sd[i].fst[j] = 0.0;
	sd[i].ust[j] = 0.0;
	sd[i].usd[j] = 0.0;
      }
      ns++;
    }
  }

  if (ns > 0) {
    fseek(f, 0, SEEK_SET);
    while (1) {
      if (NULL == fgets(buf, 1024, f)) break;
      p = buf;
      while(isspace(*p)) p++;
      if (p[0] == '#') continue;
      n = sscanf(buf, "%d %d %d %d %lg %lg %lg %lg %lg %lg %lg",
		 &zi, &ki, &n0, &n1, &e0, &e1, &d0, &d1, &f0, &f1, &w0);
      if (n != 11) continue;
      if (z0 != zi) continue;
      if (ki <= 0) continue;
      if (ki > 1000) {
	k0 = ki/1000;
	k1 = ki%1000;
      } else {
	k0 = ki;
	k1 = ki;
      }
      for (ki = k0; ki <= k1; ki++) {      
	if (ki < 1 || ki > z) continue;
	i = ki-1;
	j = (n1-sd[i].minup)*sd[i].nlo + n0-sd[i].minlo;
	sd[i].emin[j] = e0/HARTREE_EV;
	sd[i].emax[j] = e1/HARTREE_EV;
	sd[i].fde[j] = d0/HARTREE_EV;
	sd[i].ude[j] = d1/HARTREE_EV;
	sd[i].fmf[j] = f0;
	sd[i].umf[j] = f1;
	sd[i].uwf[j] = w0;
      }
    }
  }
  fclose(f);  
  printf("%d\n", ns);
  return ns;
}

		  
void CollapseDBase(char *ipr, char *opr, int k0, int k1,
		   double des00, double des01,
		   int minlevs0, int maxlevs0, int ic) {
  char ifn[1024], ofn[1024], buf[1024];
  F_HEADER fh, fh1[7];
  EN_HEADER h0;
  EN_RECORD r0;
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
  RC_HEADER h6;
  RC_RECORD r6;
  char a[8];
  char ext[7][3] = {"en", "tr", "ce", "rr", "ci", "ai", "rc"};
  int types[7] = {DB_EN, DB_TR, DB_CE, DB_RR, DB_CI, DB_AI, DB_RC};
  int k, i, j, n, nb, nlevs, vn, swp, z, nth;
  int ng0, ng1, ng, dm, ng2, ilo, iup, m, t, s;
  int klev[N_ELEMENTS1], ngrp[N_ELEMENTS1];
  int imin[N_ELEMENTS1], imax[N_ELEMENTS1], igrd[N_ELEMENTS1];
  LEVGRP *rg[N_ELEMENTS1];
  EN_RECORD *ra[N_ELEMENTS1];  
  TFILE *f0, *f1[7];
  double ei, des, te, e, cs, des0, des1, egrd[N_ELEMENTS1];
  float fte;
  int *im, tlevs, nt0, nt1, nt2, minlevs, maxlevs;
  double wt0, wt1, tt0, tt1;
  int nbk, kbk, n2, nm2;

  sprintf(ifn, "%sb.en", ipr);
  wt0 = WallTime();
  tt0 = wt0;
  printf("read levels: %s ...", ifn);
  fflush(stdout);
  f0 = OpenFileRO(ifn, &fh, &swp);
  if (f0 == NULL) {
    printf("cannot open file: %s\n", ifn);
    return;
  }
  z = (int)(fh.atom);
  strncpy(a, fh.symbol, 4);
  nth = fh.nthreads;

  if (k0 <= 0) k0 = 0;
  if (k1 <= 0) k1 = N_ELEMENTS;
  
  for (k = 0; k <= z; k++) {
    klev[k] = 0;
    imin[k] = 0;
    imax[k] = 0;
    igrd[k] = 0;
    egrd[k] = 1e30;
    ra[k] = NULL;
  }
  tlevs = 0;
  nt0 = z;
  nt1 = 0;
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadENHeader(f0, &h0, swp);
    if (n == 0) break;
    k = h0.nele;
    if (k < nt0) nt0 = k;
    if (k > nt1) nt1 = k;
    klev[k] += h0.nlevels;
    FSEEK(f0, h0.length, SEEK_CUR);
  }
  n = 0;
  for (k = 0; k <= z; k++) {
    if (klev[k] > 0) {
      ra[k] = malloc(sizeof(EN_RECORD)*klev[k]);
      n += klev[k];
      klev[k] = 0;
    }
  }
  tlevs = n;
  FCLOSE(f0);

  for (k = 0; k <= z; k++) {
    imin[k] = n;
    imax[k] = 0;
  }
  f0 = OpenFileRO(ifn, &fh, &swp);
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadENHeader(f0, &h0, swp);
    if (n == 0) break;
    k = h0.nele;
    for (i = 0; i < h0.nlevels; i++) {
      j = klev[k];
      n = ReadENRecord(f0, &(ra[k][j]), swp);
      if (i == 0) {
	if (ra[k][j].j == -1 && !IsUTARecord(&ra[k][j])) {
	  t = 1;
	} else {
	  t = 0;
	}
      }
      if (t) {
	ra[k][j].j = -2;
      }
      if (imin[k] > ra[k][j].ilev) imin[k] = ra[k][j].ilev;
      if (imax[k] < ra[k][j].ilev) imax[k] = ra[k][j].ilev;
      if (ra[k][j].energy < egrd[k]) {
	egrd[k] = ra[k][j].energy;
	igrd[k] = ra[k][j].ilev;
      }
      klev[k]++;
    }
  }
  FCLOSE(f0);
  wt1 = WallTime();
  printf(" %d %d %d %.3e\n", tlevs, nt0, nt1, wt1-wt0);

  GROUPMOD gd[N_ELEMENTS1];
  int ngm = 0;
  if (strlen(_groupmod_file) > 0) {
    ngm = LoadGroupMod(_groupmod_file, z, gd);
    if (ngm <= 0) {
      ngm = LoadGroupMod(_groupmod_file, -z, gd);
    }
  }
  for (k = z; k >= 0; k--) {
    if (gd[k].minlevs > 0 && gd[k].maxlevs > 0) {
      minlevs = gd[k].minlevs;
      maxlevs = gd[k].maxlevs;
      des0 = gd[k].des0;
      des1 = gd[k].des1;
    } else {
      minlevs = minlevs0;
      maxlevs = maxlevs0;
      des0 = des00/HARTREE_EV;
      des1 = des01/HARTREE_EV;
    }
    ngrp[k] = 0;
    if (klev[k] == 0) continue;
    wt0 = WallTime();
    printf("collapse EN: %d %d %d %.3e %.3e ...", k, minlevs, maxlevs,
	   des0*HARTREE_EV, des1*HARTREE_EV);
    fflush(stdout);
    if (klev[k] > 1) {
      j = igrd[k]-imin[k];
      if (ra[k][j].j == -1) {
	ra[k][j].j = -1000;
      }
      qsort(ra[k], klev[k], sizeof(EN_RECORD), CompareENRecordEnergySU);
      if (ra[k][0].j == -1000) {
	ra[k][0].j = -1;
      }
    }
    if (k > 0 && klev[k-1] > 0) {
      ei = ra[k-1][0].energy;
    } else {
      ei = ra[k][klev[k]-1].energy+des1;
    }
    ng0 = GroupLevels(ra[k], klev[k], ei, des0, minlevs, NULL);
    des = -1.0;
    if (ng0 > maxlevs) {
      ng1 = GroupLevels(ra[k], klev[k], ei, des1, minlevs, NULL);
    } else {
      des = des0;
      ng = ng0;
      ng1 = 0;
    }
    if (ng1 >= maxlevs) {
      des = des1;
      ng = ng1;
    }      
    if (des < 0) {
      dm = (int)(_grptol*maxlevs);
      if (dm < 1) dm = 1;
      while (des1-des0 > EPS6) {
	des = 0.5*(des0+des1);
	ng = GroupLevels(ra[k], klev[k], ei, des, minlevs, NULL);
	if (abs(ng-maxlevs) <= dm) break;
	if (ng > maxlevs) {
	  des0 = des;
	} else if (ng < maxlevs) {
	  des1 = des;
	}
      }
    }
    ngrp[k] = ng;
    rg[k] = malloc(sizeof(LEVGRP)*ng);
    ng = GroupLevels(ra[k], klev[k], ei, des, minlevs, rg[k]);
    wt1 = WallTime();
    printf(" %d %d %.2e %.3e\n", klev[k], ngrp[k], des*HARTREE_EV, wt1-wt0);
  }

  for (i = 0; i < 7; i++) {
    f1[i] = NULL;
    if (i > 0) {
      if (0 == (_collapse_mask & (1<<i))) {
	continue;
      }
      sprintf(ifn, "%sb.%s", ipr, ext[i]);
      f0 = OpenFileRO(ifn, &fh, &swp);
      if (f0 == NULL) {
	continue;
      }
      FCLOSE(f0);
    }
    sprintf(ofn, "%sb.%s", opr, ext[i]);
    fh1[i].atom = z;
    strcpy(fh1[i].symbol, a);
    fh1[i].type = types[i];
    f1[i] = OpenFileWTN(ofn, &fh1[i], nth);
  }

  iuta = 1;
  n = 0;
  for (k = z; k >= 0; k--) {
    if (ngrp[k] == 0) continue;
    h0.nele = k;
    InitFile(f1[0], &fh1[0], &h0);
    for (j = 0; j < ngrp[k]; j++) {	
      rg[k][j].r.ilev = n++;
      WriteENRecord(f1[0], &rg[k][j].r);
    }
    DeinitFile(f1[0], &fh1[0]);
  }
  CloseFile(f1[0], &fh1[0]);

  if (f1[1] || f1[2] || f1[3] || f1[4] || f1[5] || f1[6]) {
    sprintf(ifn, "%sb.en", ipr);
    MemENTable(ifn);
    im = malloc(sizeof(int)*tlevs);
    for (k = z; k >= 0; k--) {
      if (ngrp[k] == 0) continue;
      for (i = 0; i < ngrp[k]; i++) {
	for (j = 0; j < rg[k][i].nlev; j++) {
	  t = rg[k][i].ilev[j];
	  n = ra[k][t].ilev;
	  im[n] = i;
	}
      }
    }
    
    if (f1[1]) {
      TR_ALL **rt;
      int gauges[3] = {1, 2, 10};
      int i0, i1, ig, ib, nm, nm2, nmx, nmx1, nmx2, nmx3, nmx4, utr;
      double *efu;
      short **nq, *nqs, nq0, nq1;
      char c, nc0[LNCOMPLEX], *p0, *p1;
      TRANSMOD *sd;

      sd = NULL;
      nq = NULL;
      nqs = NULL;
      if (strlen(_transmod_file) > 0) {
	sd = malloc(sizeof(TRANSMOD)*z);
	nb = LoadTransMod(_transmod_file, z, sd);
	if (nb == 0) {
	  nb = LoadTransMod(_transmod_file, -z, sd);
	  if (nb == 0) {
	    free(sd);
	    sd = NULL;
	  }
	}
      }
      for (k = z; k >= 0; k--) {	
	if (ngrp[k] <= 1) continue;	
	nm = 0;
	nmx = 0;
	if (_sfu_smin >= 0) {
	  nmx = LoadSFU(ipr, k, &efu);
	}
	if (sd || nmx) {
	  nm = 0;
	  for (i = 0; i < klev[k]; i++) {
	    if (ra[k][i].p < 0) {
	      n = (-ra[k][i].p)/100;
	    } else {
	      n = (ra[k][i].p)/100;
	    }
	    if (n > nm) nm = n;
	  }
	  nm2 = (nm*(nm+1))/2;
	  n = imax[k]-imin[k]+1;
	  nq = malloc(sizeof(short *)*n);
	  for (i = 0; i < n; i++) nq[i] = NULL;
	  for (i = 0; i < klev[k]; i++) {
	    ig = ra[k][i].ilev - imin[k];
	    nq[ig] = malloc(sizeof(short)*(nm+nm2));
	    nqs = nq[ig] + nm;
	    for (j = 0; j < nm+nm2; j++) {
	      nq[ig][j] = 0;
	    }
	    j = FindNRShells(k, &ra[k][i], nm, nq[ig], nqs);
	  }
	}
	if (nmx > 0) {
	  nmx1 = nmx*nmx;
	  nmx2 = nmx1*2;
	  nmx3 = nmx1*3;
	  nmx4 = nmx1*4;
	}

	wt0 = WallTime();
	sprintf(ifn, "%sb.tr", ipr);
	for (ig = 0; ig < 3; ig++) {	  
	  f0 = OpenFileRO(ifn, &fh, &swp);
	  ib = 0;
	  for (nb = 0; nb < fh.nblocks; nb++) {
	    n = ReadTRHeader(f0, &h1, swp, &utr);
	    if (n == 0) break;
	    if (k != h1.nele) {
	      FSEEK(f0, h1.length, SEEK_CUR);
	      continue;
	    }
	    if (gauges[ig] != h1.gauge) {
	      FSEEK(f0, h1.length, SEEK_CUR);
	      continue;
	    }
	    ib = 1;
	    break;
	  }
	  FCLOSE(f0);
	  if (ib) {
	    printf("collapse TR: %d %d ...", k, gauges[ig]);
	    fflush(stdout);
	    ng2 = ngrp[k]*ngrp[k];
	    rt = malloc(sizeof(TR_ALL *)*ng2);
	    for (i = 0; i < ng2; i++) {
	      rt[i] = NULL;
	    }
	    f0 = OpenFileRO(ifn, &fh, &swp);
	    nt0 = 0;
	    for (nb = 0; nb < fh.nblocks; nb++) {
	      n = ReadTRHeader(f0, &h1, swp, &utr);
	      if (n == 0) break;
	      if (k != h1.nele) {
		FSEEK(f0, h1.length, SEEK_CUR);
		continue;
	      }
	      if (gauges[ig] != h1.gauge) {
		FSEEK(f0, h1.length, SEEK_CUR);
		continue;
	      }
	      for (i = 0; i < h1.ntransitions; i++) {
		n = ReadTRRecord(f0, &r1, &r1x, swp, utr);
		if (n == 0) break;
		ilo = im[r1.lower];
		iup = im[r1.upper];
		if (rg[k][iup].r.energy > rg[k][ilo].r.energy) {
		  j = iup*ngrp[k] + ilo;
		  if (rt[j] == NULL) {
		    rt[j] = malloc(sizeof(TR_ALL));
		    rt[j]->r.strength = 0.0;
		    rt[j]->x.energy = 0.0;
		    rt[j]->x.sdev = 0.0;
		    rt[j]->x.sci = 0.0;
		    rt[j]->r.lower = rg[k][ilo].r.ilev;
		    rt[j]->r.upper = rg[k][iup].r.ilev;
		  }
		  e = r1x.energy;
		  double gf = r1.strength;
		  double uw = 0.0;
		  if (e <= 0 || r1x.sdev <= 0) {
		    e = mem_en_table[r1.upper].energy -
		      mem_en_table[r1.lower].energy;
		  }
		  int isuta = 0;
		  if (mem_en_table[r1.lower].ibase == -mem_en_table_size-1 ||
		      mem_en_table[r1.upper].ibase == -mem_en_table_size-1) {
		    isuta = 1;
		    gf *= r1x.sci;
		    uw = r1x.sdev;
		  }
		  if (ig < 10 && (sd || (isuta && _sfu_smin >= 0 && nmx > 0))) {
		    s = r1.lower - imin[k];
		    t = r1.upper - imin[k];
		    int nlo = 0;
		    int nup = 0;
		    for (m = 0; m < nm; m++) {
		      if (nlo == 0) {
			if (nq[s][m] > nq[t][m]) {
			  nlo = m+1;
			}
		      }
		      if (nq[s][m] < nq[t][m]) {
			nup = m+1;
		      }
		    }
		    if (nlo == 0 && nup == 0) {
		      i0 = nm;
		      for (m = 1; m <= nm; m++) {
			for (i1 = 0; i1 < m; i1++) {
			  nq0 = nq[s][i0];
			  nq1 = nq[t][i0];
			  if (nq0 != nq1) {
			    nlo = m;
			    nup = nlo;
			    break;
			  }
			  i0++;
			}
			if (nlo > 0) break;
		      }
		    }		
		    if (nlo > 0 && nup > 0 && sd) {
		      int km = k-1;
		      if (sd[km].nde > 0) {
			double emin, emax, fde, ude;
			i0 = nlo - sd[km].minlo;
			i1 = nup - sd[km].minup;
			i0 = i1*sd[km].nlo + i0;
			emin = sd[km].emin[i0];
			emax = sd[km].emax[i0];
			fde = sd[km].fde[i0];
			ude = sd[km].ude[i0];
			if (e >= emin && e <= emax) {
			  if (_ste > 0) {
			    des = mem_en_table[r1.upper].energy - _eground[k];
			    cs = exp(-des/_ste)*gf;
			    if (isuta) {
			      sd[km].ume[i0] += cs*e;
			      sd[km].uwe[i0] += cs*e*e;
			      des = r1x.sdev;
			      sd[km].usd[i0] += cs*des*des;
			      sd[km].ust[i0] += cs;			    
			    } else {
			      sd[km].fme[i0] += cs*e;
			      sd[km].fwe[i0] += cs*e*e;
			      sd[km].fst[i0] += cs;	      
			    }
			  }
			  if (isuta) {
			    e += ude;
			    gf *= sd[km].umf[i0];
			    uw *= sd[km].uwf[i0];
			  } else {
			    e += fde;
			    gf *= sd[km].fmf[i0];
			  }			    
			}
		      }
		    }
		    if (isuta && _sfu_smin >= 0 && nmx > 0) {
		      if (nlo == 0 && nup == 0) {
			des = 1e30;
			for (m = 0; m < nmx; m++) {
			  nb = m*nmx + m;
			  if (efu[nb+nmx3] > 0) {
			    cs = fabs(e - efu[nb+nmx1]);
			    if (cs < des) {
			      des = cs;
			      nlo = m+1;
			    }
			  }
			}
			nup = nlo;
		      }
		      if (nlo > 0 && nup > 0 && nlo <= nmx && nup <= nmx) {
			nb = (nup-1)*nmx + nlo-1;
			if (efu[nb+nmx3] > _sfu_smin &&
			    efu[nb+nmx2] > _sfu_smin) {
			  des = fabs(e - efu[nb+nmx1]);
			  if (des < _sfu_dmax) {
			    cs = efu[nb]-efu[nb+nmx1];
			    cs = Max(cs, _sfu_minde);
			    cs = Min(cs, _sfu_maxde);
			    cs = 1 + cs/efu[nb+nmx1];
			    cs = Max(1e-2, cs);
			    cs = Min(1e2, cs);
			    e *= cs;
			  }
			}
		      }
		    }		    
		  }
		  if (gf > 0) {
		    rt[j]->x.sci += gf*e*e;
		    rt[j]->r.strength += gf;
		    rt[j]->x.energy += e*gf;
		    rt[j]->x.sdev += uw*uw*gf;
		  }
		}		
		nt0++;
	      }
	    }
	    FCLOSE(f0);
	    nt1 = 0;
	    h1.nele = k;
	    h1.gauge = gauges[ig];
	    InitFile(f1[1], &fh1[1], &h1);
	    for (i = 0; i < ng2; i++) {
	      if (rt[i]) {
		if (rt[i]->r.strength > 0) {
		  rt[i]->x.energy /= rt[i]->r.strength;
		  rt[i]->x.sdev /= rt[i]->r.strength;
		  rt[i]->x.sci /= rt[i]->r.strength;
		  e = rt[i]->x.sci - rt[i]->x.energy*rt[i]->x.energy;
		  if (e <= 0.0) e = 0.0;
		  else e *= _cwf;
		  e = sqrt(e + rt[i]->x.sdev);
		  rt[i]->x.sdev = e;
		  rt[i]->x.sci = 1.0;
		  WriteTRRecord(f1[1], &(rt[i]->r), &(rt[i]->x), iuta);
		  nt1++;
		}
		free(rt[i]);
	      }
	    }
	    DeinitFile(f1[1], &fh1[1]);
	    free(rt);
	    wt1 = WallTime();
	    printf(" %d %d %.3e\n", nt0, nt1, wt1-wt0);
	  }
	}
	if (nmx > 0) {
	  free(efu);
	}
	if (nq) {
	  n = imax[k]-imin[k]+1;
	  for (i = 0; i < n; i++) {
	    if (nq[i]) free(nq[i]);
	  }
	  free(nq);
	}		
      }
      FILE *f = NULL;
      if (sd) {
	sprintf(buf, "%s.out", _transmod_file);
	f = fopen(buf, "w");
      }
      for (i = 0; i < z; i++) {
	k = i+1;
	if (ngrp[k] <= 1) continue;
	if (sd && sd[i].nde > 0) {
	  if (f) {
	    for (ilo = sd[i].minlo; ilo <= sd[i].maxlo; ilo++) {
	      for (iup = sd[i].minup; iup <= sd[i].maxup; iup++) {
		m = (iup-sd[i].minup)*sd[i].nlo + (ilo-sd[i].minlo);
		if (sd[i].emax[m] > 0) {
		  if (sd[i].fst[m] > 0) {
		    sd[i].fme[m] /= sd[i].fst[m];
		    sd[i].fwe[m] /= sd[i].fst[m];
		    sd[i].fwe[m] = sd[i].fwe[m] - sd[i].fme[m]*sd[i].fme[m];
		    if (sd[i].fwe[m] < 0) sd[i].fwe[m] = 0.0;
		    sd[i].fwe[m] = sqrt(sd[i].fwe[m]);
		  }
		  if (sd[i].ust[m] > 0) {
		    sd[i].ume[m] /= sd[i].ust[m];
		    sd[i].uwe[m] /= sd[i].ust[m];
		    sd[i].usd[m] /= sd[i].ust[m];
		    sd[i].uwe[m] = sd[i].uwe[m] - sd[i].ume[m]*sd[i].ume[m];
		    if (sd[i].uwe[m] < 0) sd[i].uwe[m] = 0.0;
		    sd[i].uwe[m] = sqrt(sd[i].uwe[m] + sd[i].usd[m]);
		  }
		  fprintf(f, "%3d %3d %2d %2d %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E\n",
			  z, k, ilo, iup,
			  sd[i].emin[m]*HARTREE_EV,
			  sd[i].emax[m]*HARTREE_EV,
			  sd[i].fme[m]*HARTREE_EV,
			  sd[i].fwe[m]*HARTREE_EV,
			  sd[i].fst[m],
			  sd[i].ume[m]*HARTREE_EV,
			  sd[i].uwe[m]*HARTREE_EV,
			  sqrt(sd[i].usd[m])*HARTREE_EV,
			  sd[i].ust[m]);
		}
	      }
	    }
	  }
	  free(sd[i].emin);
	  free(sd[i].emax);
	  free(sd[i].fde);
	  free(sd[i].ude);
	  free(sd[i].fme);
	  free(sd[i].ume);
	  free(sd[i].fmf);
	  free(sd[i].umf);
	  free(sd[i].uwf);
	  free(sd[i].fst);
	  free(sd[i].ust);
	  free(sd[i].fwe);
	  free(sd[i].uwe);
	  free(sd[i].usd);
	}
      }
      if (sd) free(sd);
      if (f) fclose(f);
    }
    if (f1[2]) {
      CE_RECORD **rt;
      CE_HEADER ht;
      double data[2+(1+MAXNUSR)*3], ratio;
      int neg = 6;
      SetCEEGridType(1);
      SetUsrCEEGridType(1);
      SetCEPWGridType(0);
      SetCETEGrid(-1, 0.0, 0.0);
      SetCEEGrid(neg, 0.05, 25.0, 1.0);
      SetUsrCEEGrid(neg, 0.05, 25.0, 1.0);
      for (k = z; k >= 0; k--) {
	if (ngrp[k] == 0) continue;
	wt0 = WallTime();
	printf("collapse CE: %d ...", k);
	fflush(stdout);
	PrepCEHeader(&ht, k, 0);
	ht.nparams = 0;
	ng2 = ngrp[k]*ngrp[k];
	rt = malloc(sizeof(CE_RECORD *)*ng2);
	for (i = 0; i < ng2; i++) {
	  rt[i] = NULL;
	}
	nt0 = 0;
	sprintf(ifn, "%sb.ce", ipr);
	f0 = OpenFileRO(ifn, &fh, &swp);
	for (nb = 0; nb < fh.nblocks; nb++) {
	  n = ReadCEHeader(f0, &h2, swp);
	  if (n == 0) break;
	  if (k != h2.nele) {
	    FSEEK(f0, h2.length, SEEK_CUR);
	    free(h2.tegrid);
	    free(h2.egrid);
	    free(h2.usr_egrid);
	    continue;
	  }
	  PrepCECrossHeader(&h2, data);
	  for (i = 0; i < h2.ntransitions; i++) {
	    n = ReadCERecord(f0, &r2, swp, &h2);
	    if (n == 0) break;
	    ilo = im[r2.lower];
	    iup = im[r2.upper];
	    if (rg[k][iup].r.energy > rg[k][ilo].r.energy) {
	      j = iup*ngrp[k] + ilo;
	      if (rt[j] == NULL) {
		rt[j] = malloc(sizeof(CE_RECORD));
		rt[j]->strength = malloc(sizeof(float)*neg);
		for (t = 0; t < neg; t++) {
		  rt[j]->strength[t] = 0.0;
		}
		rt[j]->bethe = 0.0;
		rt[j]->born[0] = 0.0;
		rt[j]->born[1] = 0.0;
		rt[j]->nsub = 1;
		rt[j]->lower = rg[k][ilo].r.ilev;
		rt[j]->upper = rg[k][iup].r.ilev;
	      }
	      PrepCECrossRecord(0, &r2, &h2, data);
	      te = mem_en_table[r2.upper].energy - mem_en_table[r2.lower].energy;
	      for (t = 0; t < neg; t++) {
		e = te*ht.egrid[t]*HARTREE_EV;
		cs = InterpolateCECross(e, &r2, &h2, data, &ratio);
		rt[j]->strength[t] += cs;	      
	      }
	      if (r2.bethe >= 0.0) {
		rt[j]->bethe += r2.bethe;
		rt[j]->born[0] += r2.born[0];
		rt[j]->born[1] += r2.born[1]*r2.born[0];
	      }
	    }
	    free(r2.strength);
	    nt0++;
	  }
	  free(h2.tegrid);
	  free(h2.egrid);
	  free(h2.usr_egrid);
	}
	FCLOSE(f0);
	nt1 = 0;
	InitFile(f1[2], &fh1[2], &ht);
	for (i = 0; i < ng2; i++) {
	  if (rt[i] == NULL) continue;
	  t = 0;
	  for (j = 0; j < neg; j++) {
	    if (rt[i]->strength[j] > 0) {
	      t = 1;
	      break;
	    }
	  }
	  if (t) {
	    if (rt[i]->born[0] > 0) {
	      rt[i]->born[1] /= rt[i]->born[0];
	    } else {
	      rt[i]->bethe = -1.0;
	    }
	    WriteCERecord(f1[2], rt[i]);
	    nt1++;
	  }
	  free(rt[i]->strength);
	  free(rt[i]);
	}
	DeinitFile(f1[2], &fh1[2]);
	free(rt);
	wt1 = WallTime();
	printf(" %d %d %.3e\n", nt0, nt1, wt1-wt0);
      }
    }

    if (f1[3]) {
      RR_ALL **rt;
      RR_HEADER ht;
      int neg = 6, nqk, nt10;
      int iwe = 0;
      int iwk = neg-2;
      SetUsrPEGridType(1);
      SetRRTEGrid(-1, 0.0, 0.0);
      SetPEGrid(neg, 0.05, 10.0, 1.0);
      SetUsrPEGrid(neg, 0.05, 10.0, 1.0);
      SetRecQkMode(QK_DEFAULT, 0.1);
      for (k = z; k > 0; k--) {
	if (ngrp[k] == 0 || ngrp[k-1] == 0) continue;		
	short *nrmin, *nrmax;
	nrmin = malloc(sizeof(short)*ngrp[k]);
	nrmax = malloc(sizeof(short)*ngrp[k]);	
	wt0 = WallTime();
	printf("collapse RR: %d ...", k);
	fflush(stdout);
	nqk = PrepRRHeader(&ht, k, -1);
	ng2 = ngrp[k]*ngrp[k-1];
	rt = malloc(sizeof(RR_ALL *)*ng2);
	for (i = 0; i < ng2; i++) {
	  rt[i] = NULL;
	}
	for (i = 0; i < ngrp[k]; i++) {
	  nrmin[i] = 1000;
	  nrmax[i] = 0;
	}
	nt0 = 0;
	nt10 = 0;
	sprintf(ifn, "%sb.rr", ipr);
	f0 = OpenFileRO(ifn, &fh, &swp);
	for (nb = 0; nb < fh.nblocks; nb++) {
	  n = ReadRRHeader(f0, &h3, swp);
	  if (n == 0) break;
	  if (k != h3.nele) {
	    FSEEK(f0, h3.length, SEEK_CUR);
	    free(h3.tegrid);
	    free(h3.egrid);
	    free(h3.usr_egrid);
	    continue;
	  }
	  for (i = 0; i < h3.ntransitions; i++) {
	    n = ReadRRRecord(f0, &r3, swp, &h3);
	    if (n == 0) break;
	    if (r3.f >= 0) {
	      ilo = im[r3.b];
	      iup = im[r3.f];
	      if (rg[k-1][iup].r.energy > rg[k][ilo].r.energy) {
		j = iup*ngrp[k] + ilo;
		if (rt[j] == NULL) {
		  rt[j] = malloc(sizeof(RR_ALL));
		  rt[j]->dk = 0.0;
		  rt[j]->dn = 0.0;
		  rt[j]->r.strength = malloc(sizeof(float)*neg);
		  if (nqk > 0) rt[j]->r.params = malloc(sizeof(float)*nqk);
		  for (t = 0; t < neg; t++) {
		    rt[j]->r.strength[t] = 0.0;
		  }
		  for (t = 0; t < nqk; t++) {
		    rt[j]->r.params[t] = 0.0;
		  }
		  rt[j]->r.kl = 0;
		  rt[j]->r.b = rg[k][ilo].r.ilev;
		  rt[j]->r.f = rg[k-1][iup].r.ilev;
		}
		te = mem_en_table[r3.f].energy - mem_en_table[r3.b].energy;
		for (t = 0; t < neg; t++) {
		  e = te*ht.egrid[t];
		  cs = InterpolateRRCross(e, te, &r3, &h3);
		  rt[j]->r.strength[t] += cs;
		  if (t == iwk && cs > 0) {
		    rt[j]->dk += cs*(r3.kl%1000);		  
		    rt[j]->dn += cs*(r3.kl/1000);
		  }
		}
		if (nqk > 0) {
		  rt[j]->r.params[0] += r3.params[0];
		  rt[j]->r.params[1] += r3.params[1]*r3.params[0];
		  rt[j]->r.params[2] += r3.params[2]*r3.params[0];
		  rt[j]->r.params[3] += r3.params[3]*r3.params[0];
		}
	      }
	      nt0++;
	    } else {
	      n = (-r3.f)/1000;
	      t = im[r3.b];
	      if (n < nrmin[t]) nrmin[t] = n;
	      if (n > nrmax[t]) nrmax[t] = n;
	      nt10++;
	    }
	    free(r3.strength);
	    if (nqk > 0) free(r3.params);
	  }
	  free(h3.tegrid);
	  free(h3.egrid);
	  free(h3.usr_egrid);
	}
	FCLOSE(f0);
	nt1 = 0;
	InitFile(f1[3], &fh1[3], &ht);
	for (i = 0; i < ng2; i++) {
	  if (rt[i] == NULL) continue;
	  if (nqk > 0 && rt[i]->r.params[0] > 0) {
	    rt[i]->r.params[1] /= rt[i]->r.params[0];
	    rt[i]->r.params[2] /= rt[i]->r.params[0];
	    rt[i]->r.params[3] /= rt[i]->r.params[0];
	  }
	  if (rt[i]->r.strength[iwk] > 0) {
	    t = (int)(rt[i]->dk/rt[i]->r.strength[iwk]+0.25);
	    s = (int)(rt[i]->dn/rt[i]->r.strength[iwk]+0.25);
	    rt[i]->r.kl = s*1000 + t;	  
	    WriteRRRecord(f1[3], &(rt[i]->r));
	    nt1++;
	  }
	  free(rt[i]->r.strength);
	  if (nqk > 0) free(rt[i]->r.params);
	  free(rt[i]);
	}
	DeinitFile(f1[3], &fh1[3]);	
	free(rt);
	wt1 = WallTime();
	printf(" %d %d %.3e ...", nt0, nt1, wt1-wt0);
	fflush(stdout);
	wt0 = wt1;
	if (nt10 > 0) {
	  rt = malloc(sizeof(RR_ALL *)*ngrp[k]);
	  for (i = 0; i < ngrp[k]; i++) {
	    rt[i] = NULL;
	  }
	  nt2 = 0;
	  f0 = OpenFileRO(ifn, &fh, &swp);
	  for (nb = 0; nb < fh.nblocks; nb++) {
	    n = ReadRRHeader(f0, &h3, swp);
	    if (n == 0) break;
	    if (k != h3.nele) {
	      FSEEK(f0, h3.length, SEEK_CUR);
	      free(h3.tegrid);
	      free(h3.egrid);
	      free(h3.usr_egrid);
	      continue;
	    }
	    for (i = 0; i < h3.ntransitions; i++) {
	      n = ReadRRRecord(f0, &r3, swp, &h3);
	      if (n == 0) break;
	      if (r3.f < 0) {
		ilo = im[r3.b];
		int kl = -r3.f;
		fte = *((float *)(&(r3.kl)));
		s = ilo;
		nm2 = (nrmin[s]-1)*nrmin[s]/2;
		if (rt[s] == NULL) {
		  n = 1 + nrmax[s] - nrmin[s];
		  n2 = n*(nrmax[s]+nrmin[s])/2;
		  rt[s] = malloc(sizeof(RR_ALL)*n2);
		  for (nbk = nrmin[s]; nbk <= nrmax[s]; nbk++) {
		    for (kbk = 0; kbk < nbk; kbk++) {
		      t = (nbk-1)*nbk/2 + kbk - nm2;
		      rt[s][t].r.strength = malloc(sizeof(float)*neg);
		      if (nqk > 0) {
			rt[s][t].r.params = malloc(sizeof(float)*nqk);
		      }
		      for (j = 0; j < neg; j++) {
			rt[s][t].r.strength[j] = 0.0;
		      }
		      for (j = 0; j < nqk; j++) {
			rt[s][t].r.params[j] = 0.0;
		      }
		      rt[s][t].r.kl = 0;
		      rt[s][t].r.b = rg[k][ilo].r.ilev;
		      rt[s][t].r.f = -(nbk*1000 + kbk);
		      rt[s][t].dk = 0.0;
		      rt[s][t].dn = 0.0;
		    }
		  }
		}
		nbk = kl/1000;
		kbk = kl%1000;
		j = (nbk-1)*nbk/2 + kbk - nm2;
		te = fte;
		for (t = 0; t < neg; t++) {
		  e = te*ht.egrid[t];
		  cs = InterpolateRRCross(e, te, &r3, &h3);
		  rt[s][j].r.strength[t] += cs;
		  if (t == iwk && cs > 0) {
		    rt[s][j].dk += cs*(kl%1000);
		  }
		  if (t == iwe && cs > 0) {
		    rt[s][j].dn += cs*fte;
		    nt2++;
		  }
		}
		if (nqk > 0) {
		  rt[s][j].r.params[0] += r3.params[0];
		  rt[s][j].r.params[1] += r3.params[1]*r3.params[0];
		  rt[s][j].r.params[2] += r3.params[2]*r3.params[0];
		  rt[s][j].r.params[3] += r3.params[3]*r3.params[0];
		}
	      }
	      free(r3.strength);
	      if (nqk > 0) free(r3.params);
	    }
	    free(h3.tegrid);
	    free(h3.egrid);
	    free(h3.usr_egrid);
	  }
	  FCLOSE(f0);
	  nt2 = 0;
	  InitFile(f1[3], &fh1[3], &ht);
	  for (i = 0; i < ngrp[k]; i++) {
	    if (rt[i] == NULL) continue;
	    n = nrmax[i] - nrmin[i] + 1;
	    n2 = n*(nrmin[i]+nrmax[i])/2;
	    nm2 = (nrmin[i]-1)*nrmin[i]/2;
	    for (nbk = nrmin[i]; nbk <= nrmax[i]; nbk++) {
	      for (kbk = 0; kbk < nbk; kbk++) {
		j = (nbk-1)*nbk/2 + kbk - nm2;
		if (nqk > 0 && rt[i][j].r.params[0] > 0) {
		  rt[i][j].r.params[1] /= rt[i][j].r.params[0];
		  rt[i][j].r.params[2] /= rt[i][j].r.params[0];
		  rt[i][j].r.params[3] /= rt[i][j].r.params[0];
		}
		if (rt[i][j].r.strength[iwk] > 0 &&
		    rt[i][j].r.strength[iwe] > 0) {
		  fte = rt[i][j].dn/rt[i][j].r.strength[iwe];
		  rt[i][j].r.kl = *((int *) &fte);
		  WriteRRRecord(f1[3], &(rt[i][j].r));
		  nt2++;
		}
		free(rt[i][j].r.strength);
		if (nqk > 0) free(rt[i][j].r.params);
	      }
	    }
	    free(rt[i]);
	  }
	  DeinitFile(f1[3], &fh1[3]);
	  free(rt);
	} else {
	  nt2 = 0;
	}
	free(nrmin);
	free(nrmax);
	wt1 = WallTime();
	printf(" %d %d %.3e\n", nt10, nt2, wt1-wt0);
      }
    }

    if (f1[4]) {
      CI_ALL **rt;
      CI_HEADER ht;
      int neg = 6;
      int nqk;
      int iwk = neg-2;
      SetUsrCIEGridType(1);
      SetIEGrid(-1, 0.0, 0.0);
      SetCIEGrid(neg, 0.05, 10.0, 1.0);
      SetUsrCIEGrid(neg, 0.05, 10.0, 1.0);
      SetCIQkMode(QK_DEFAULT, 0.1);
      for (k = z; k > 0; k--) {
	if (ngrp[k] == 0 || ngrp[k-1] == 0) continue;
	wt0 = WallTime();
	printf("collapse CI: %d ...", k);
	fflush(stdout);
	nqk = PrepCIHeader(&ht, k);
	ng2 = ngrp[k]*ngrp[k-1];
	rt = malloc(sizeof(CI_ALL *)*ng2);
	for (i = 0; i < ng2; i++) {
	  rt[i] = NULL;
	}
	nt0 = 0;
	sprintf(ifn, "%sb.ci", ipr);
	f0 = OpenFileRO(ifn, &fh, &swp);
	for (nb = 0; nb < fh.nblocks; nb++) {
	  n = ReadCIHeader(f0, &h4, swp);
	  if (n == 0) break;
	  if (k != h4.nele) {
	    FSEEK(f0, h4.length, SEEK_CUR);
	    free(h4.tegrid);
	    free(h4.egrid);
	    free(h4.usr_egrid);
	    continue;
	  }
	  for (i = 0; i < h4.ntransitions; i++) {
	    n = ReadCIRecord(f0, &r4, swp, &h4);
	    if (n == 0) break;
	    ilo = im[r4.b];
	    iup = im[r4.f];
	    if (rg[k-1][iup].r.energy > rg[k][ilo].r.energy) {
	      j = iup*ngrp[k] + ilo;
	      te = mem_en_table[r4.f].energy - mem_en_table[r4.b].energy;	      
	      if (rt[j] == NULL) {
		rt[j] = malloc(sizeof(CI_ALL));
		rt[j]->r.strength = malloc(sizeof(float)*neg);
		rt[j]->r.params = malloc(sizeof(float)*nqk);
		for (t = 0; t < neg; t++) {
		  rt[j]->r.strength[t] = 0.0;
		}
		for (t = 0; t < nqk; t++) {
		  rt[j]->r.params[t] = 0.0;
		}
		rt[j]->dk = 0.0;
		rt[j]->r.b = rg[k][ilo].r.ilev;
		rt[j]->r.f = rg[k-1][iup].r.ilev;
	      }
	      for (t = 0; t < neg; t++) {
		e = te*ht.egrid[t];
		cs = InterpolateCICross(e, te, BornMass(), &r4, &h4);
		rt[j]->r.strength[t] += cs;
		if (t == iwk) {
		  rt[j]->dk += cs*r4.kl;
		}
	      }
	      for (t = 0; t < nqk; t++) {
		rt[j]->r.params[t] += r4.params[t];
	      }
	    }
	    free(r4.strength);
	    free(r4.params);
	    nt0++;
	  }
	  free(h4.tegrid);
	  free(h4.egrid);
	  free(h4.usr_egrid);
	}
	FCLOSE(f0);
	nt1 = 0;
	InitFile(f1[4], &fh1[4], &ht);
	for (i = 0; i < ng2; i++) {
	  if (rt[i] == NULL) continue;
	  if (rt[i]->r.strength[iwk] > 0) {
	    rt[i]->dk /= rt[i]->r.strength[iwk];
	    rt[i]->r.kl = (int)(rt[i]->dk+0.25);
	    WriteCIRecord(f1[4], &(rt[i]->r));
	    nt1++;
	  }
	  free(rt[i]->r.strength);
	  free(rt[i]->r.params);
	  free(rt[i]);
	}
	DeinitFile(f1[4], &fh1[4]);
	free(rt);
	wt1 = WallTime();
	printf(" %d %d %.3e\n", nt0, nt1, wt1-wt0);
      }      
    }

    if (f1[5]) {
      AI_RECORD **rt;
      AI_HEADER ht;
      SetPEGrid(-1, 0.0, 0.0, 0.0);
      for (k = z; k > 0; k--) {
	if (ngrp[k] == 0 || ngrp[k-1] == 0) continue;
	wt0 = WallTime();
	printf("collapse AI: %d ...", k);
	fflush(stdout);
	PrepAIHeader(&ht, k, 0.0);
	ng2 = ngrp[k]*ngrp[k-1];
	rt = malloc(sizeof(AI_RECORD *)*ng2);
	for (i = 0; i < ng2; i++) {
	  rt[i] = NULL;
	}

	nt0 = 0;
	sprintf(ifn, "%sb.ai", ipr);
	f0 = OpenFileRO(ifn, &fh, &swp);
	for (nb = 0; nb < fh.nblocks; nb++) {
	  n = ReadAIHeader(f0, &h5, swp);
	  if (n == 0) break;
	  if (k != h5.nele) {
	    FSEEK(f0, h5.length, SEEK_CUR);
	    continue;
	  }
	  for (i = 0; i < h5.ntransitions; i++) {
	    n = ReadAIRecord(f0, &r5, swp);
	    if (n == 0) break;
	    ilo = im[r5.b];
	    iup = im[r5.f];
	    if (rg[k-1][iup].r.energy < rg[k][ilo].r.energy) {
	      j = iup*ngrp[k] + ilo;
	      if (rt[j] == NULL) {
		rt[j] = malloc(sizeof(AI_RECORD));
		rt[j]->rate = 0.0;
		rt[j]->b = rg[k][ilo].r.ilev;
		rt[j]->f = rg[k-1][iup].r.ilev;
	      }	      
	      cs = r5.rate*(mem_en_table[r5.b].j+1.0);
	      rt[j]->rate += cs;
	    }
	    nt0++;
	  }
	}
	FCLOSE(f0);
	InitFile(f1[5], &fh1[5], &ht);
	nt1 = 0;
	for (i = 0; i < ng2; i++) {
	  if (rt[i] == NULL) continue;
	  if (rt[i]->rate > 0) {
	    ilo = i%ngrp[k];
	    rt[i]->rate /= rg[k][ilo].r.ibase+1.0;
	    WriteAIRecord(f1[5], rt[i]);
	    nt1++;
	  }
	  free(rt[i]);
	}
	DeinitFile(f1[5], &fh1[5]);
	free(rt);
	wt1 = WallTime();
	printf(" %d %d %.3e\n", nt0, nt1, wt1-wt0);
      }
    }

    if (f1[6]) {
      RC_RECORD **rt;
      int irc, ib, ntd, nte, nde;
      double wgt;
      sprintf(ifn, "%sb.rc", ipr);
      for (k = z; k > 0; k--) {
	if (ngrp[k] == 0) continue;
	for (irc = 1; irc <= RC_TT; irc++) {
	  f0 = OpenFileRO(ifn, &fh, &swp);
	  ib = 0;
	  for (nb = 0; nb < fh.nblocks; nb++) {
	    n = ReadRCHeader(f0, &h6, swp);
	    if (n == 0) break;
	    if (k != h6.nele) {
	      FSEEK(f0, h6.length, SEEK_CUR);
	      continue;
	    }
	    if (irc != h6.type) {
	      FSEEK(f0, h6.length, SEEK_CUR);
	      continue;
	    }
	    nte = h6.nte;
	    nde = h6.nde;
	    ib = 1;
	    break;
	  }
	  FCLOSE(f0);
	  if (ib) {
	    ntd = nte*nde;
	    wt0 = WallTime();
	    printf("collapse RC: %d %d %d %d ...", k, irc, nte, nde);
	    fflush(stdout);
	    if (irc == RC_CE) {
	      ng2 = ngrp[k]*ngrp[k];
	    } else if (irc == RC_RE) {
	      ng2 = ngrp[k-1]*ngrp[k-1];
	    } else {
	      ng2 = ngrp[k]*ngrp[k-1];
	    }
	    if (ng2 > 0) {
	      rt = malloc(sizeof(RC_RECORD *)*ng2);
	      for (i = 0; i < ng2; i++) {
		rt[i] = NULL;
	      }
	      nt0 = 0;
	      f0 = OpenFileRO(ifn, &fh, &swp);
	      for (nb = 0; nb < fh.nblocks; nb++) {
		n = ReadRCHeader(f0, &h6, swp);
		if (n == 0) break;
		if (k != h6.nele) {
		  FSEEK(f0, h6.length, SEEK_CUR);
		  continue;
		}
		if (irc != h6.type) {
		  FSEEK(f0, h6.length, SEEK_CUR);
		  continue;
		}
		for (i = 0; i < h6.ntransitions; i++) {
		  n = ReadRCRecord(f0, &r6, swp, &h6);
		  if (n == 0) break;
		  ilo = im[r6.lower];
		  iup = im[r6.upper];
		  if (irc == RC_CE) {
		    if (rg[k][iup].r.energy <= rg[k][ilo].r.energy) continue;
		  } else if (irc == RC_RE) {
		    if (rg[k-1][iup].r.energy <= rg[k-1][ilo].r.energy) continue;
		  } else {
		    if (rg[k-1][iup].r.energy <= rg[k][ilo].r.energy) continue;
		  }
		  if (irc == RC_RE) {
		    j = iup*ngrp[k-1] + ilo;
		  } else {
		    j = iup*ngrp[k] + ilo;
		  }
		  if (rt[j] == NULL) {
		    rt[j] = malloc(sizeof(RC_RECORD));
		    rt[j]->rc = malloc(sizeof(float)*ntd);
		    for (t = 0; t < ntd; t++) {
		      rt[j]->rc[t] = 0.0;
		    }
		    if (irc == RC_CE) {
		      rt[j]->lower = rg[k][ilo].r.ilev;
		      rt[j]->upper = rg[k][iup].r.ilev;
		    } else if (irc == RC_RE) {
		      rt[j]->lower = rg[k-1][ilo].r.ilev;
		      rt[j]->upper = rg[k-1][iup].r.ilev;
		    } else {
		      rt[j]->lower = rg[k][ilo].r.ilev;
		      rt[j]->upper = rg[k-1][iup].r.ilev;
		    }		      
		  }
		  wgt = mem_en_table[r6.upper].j+1.0;
		  for (t = 0; t < ntd; t++) {
		    rt[j]->rc[t] += r6.rc[t]*wgt;
		  }
		  nt0++;
		}
		free(r6.rc);
	      }
	    }
	    FCLOSE(f0);
	    nt1 = 0;
	    h6.type = irc;
	    h6.nde = nde;
	    h6.nte = nte;
	    InitFile(f1[6], &fh1[6], &h6);
	    for (i = 0; i < ng2; i++) {
	      if (rt[i]) {
		ilo = i%ngrp[k];
		iup = i/ngrp[k];
		if (irc == RC_CE || irc == RC_RE) {
		  wgt = rg[k][iup].r.ibase + 1.0;
		} else {
		  wgt = rg[k-1][iup].r.ibase + 1.0;
		}
		for (t = 0; t < ntd; t++) {
		  rt[i]->rc[t] /= wgt;
		}
		WriteRCRecord(f1[6], rt[i]);
		free(rt[i]->rc);
		free(rt[i]);
		nt1++;
	      }
	    }
	    DeinitFile(f1[6], &fh1[6]);
	    free(rt);
	    wt1 = WallTime();
	    printf(" %d %d %.3e\n", nt0, nt1, wt1-wt0);
	  }
	}
      }
    }
    
    free(im);
  }

  sprintf(ifn, "%sb.rp", ipr);
  FILE *f = fopen(ifn, "r");
  if (f != NULL) {
    printf("copy %sb.rp\n", ipr);
    sprintf(buf, "cp %sb.rp %sb.rp", ipr, opr);
    i = system(buf);
    fclose(f);
  }
  
  for (k = 0; k <= z; k++) {
    if (klev[k] > 0) free(ra[k]);
    if (ngrp[k] > 0) {
      for (i = 0; i < ngrp[k]; i++) {
	if (rg[k][i].nlev > 0) free(rg[k][i].ilev);
      }
      free(rg[k]);
    }
  }
  
  for (i = 1; i < 7; i++) {
    if (f1[i]) CloseFile(f1[i], &fh1[i]);
  }

  if (ic > 0) {
    wt0 = WallTime();
    FreeMemENTable();
    n = ic;
    k = 0;    
    for (i = 0; i < 7; i++) {
      if (n < 10) {	
	k = i < n;
      } else {
	k = ic%10;
	ic = ic/10;
      }
      if (k) {
	wt0 = WallTime();
	sprintf(ifn, "%sb.%s", opr, ext[i]);
	sprintf(ofn, "%sa.%s", opr, ext[i]);
	printf("print table %s %s ...", ifn, ofn);
	fflush(stdout);
	PrintTable(ifn, ofn, 1);
	wt1 = WallTime();
	printf(" %.3e\n", wt1-wt0);
      }
    }
  }
  tt1 = WallTime();
  printf("total time: %.3e\n", tt1-tt0);
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

void InitIdxDat(void *d, int nb) {
  IDXDAT *p;
  p = (IDXDAT *) d;
  int i;
  for (i = 0; i < nb; i++) {
    p[i].i = -1;
    p[i].e = 0;
  }
}

int PreloadEN(char *fn, int i0, int i1, int j0, int j1) {
  FILE *f;
  char buf[2048];
  int i, im, n, nb;
  double e, em, de;
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
    n = sscanf(buf, "%d %d %lf %lf %lf", &im, &i, &em, &e, &de);
    if (n != 5) continue;
    nb = strlen(buf);
    if (nb < 2) continue;
    if (buf[nb-1] != '\n') continue;
    if (1+e==1 && im > 0) continue;
    e /= HARTREE_EV;
    d.i = im;
    if (im == 0) {
      d.e = 0.0;
    } else {
      d.e = e;
    }
    if (i > 0 || im == 0) {
      ArraySet(_idxmap.imap, i, &d, InitIdxDat);
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

  if (_idxmap.ni > 0 && _idxmap.nj > 0 && _idxmap.nij > 0) {
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
  IDXDAT *d = ArrayGet(_idxmap.imap, i);
  if (d && d->i >= 0) return d;
  return NULL;
}

int PreloadTable(char *tfn, char *sfn, int m) {
  TFILE *f0;
  F_HEADER fh;
  int swp, iu0, iu1;

  iu0 = iuta;
  f0 = OpenFileRO(sfn, &fh, &swp);
  iu1 = iuta;
  if (f0 == NULL) {
    printf("cannot open file %s\n", sfn);
    return -1;
  }
  FCLOSE(f0);
  iuta = iu0;
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
  int ib, i, j, k, t, n, swp, iu0, iu1, utr;
  float *rt, *et;
  F_HEADER fh, fh1;
  TR_HEADER h;
  TR_RECORD r;
  TR_EXTRA rx;
  IDXDAT *di, *dj;

  if (_idxmap.nij == 0) {
    printf("index map not loaded\n");
    return -1;
  }
  iu0 = iuta;
  f0 = OpenFileRO(sfn, &fh, &swp);
  iu1 = iuta;
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
  iuta = iu0;
  f1 = OpenFile(tfn, &fh1);
  if (f1 == NULL) {
    printf("cannot open file %s\n", tfn);
    FCLOSE(f0);
    return -1;
  }
  if (m == 0) {
    rt = malloc(sizeof(float)*_idxmap.nij);
    et = malloc(sizeof(float)*_idxmap.nij);
    for (i = 0; i < _idxmap.nij; i++) {
      rt[i] = 0.0;
      et[i] = 0.0;
    }
  }
  rx.energy = 0.0;
  rx.sdev = 0.0;
  rx.sci = 1.0;
  for (ib = 0; ib < fh.nblocks; ib++) {
    iuta = iu1;
    n = ReadTRHeader(f0, &h, swp, &utr);
    if (n == 0) break;
    if (m) {
      iuta = iu0;
      InitFile(f1, &fh1, &h);
    }
    for (k = 0; k < h.ntransitions; k++) {
      iuta = iu1;
      n = ReadTRRecord(f0, &r, &rx, swp, utr);
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
	  et[t] = dj->e - di->e;
	  SetPreloadedTR(r.lower, r.upper, 0);
	} else {
	  iuta = iu0;
	  if (iu1 == 0) rx.energy = dj->e - di->e;
	  WriteTRRecord(f1, &r, &rx, iuta);
	}
      }
    }
    if (m) {
      iuta = iu0;
      DeinitFile(f1, &fh1);
    }
  }
  FCLOSE(f0);
  if (m == 0) {
    h.multipole = 0;
    iuta = iu0;
    InitFile(f1, &fh1, &h);
    for (i = _idxmap.im0; i <= _idxmap.im1; i++) {
      for (j = _idxmap.jm0; j <= _idxmap.jm1; j++) {
	if (IsPreloadedTR(i, j, 0)) {
	  t = (j-_idxmap.jm0)*_idxmap.ni + i-_idxmap.im0;	
	  r.lower = i;
	  r.upper = j;
	  r.strength = rt[t];
	  iuta = iu0;
	  if (iu1 == 0) rx.energy = et[t];
	  WriteTRRecord(f1, &r, &rx, iuta);
	}
      }
    }
    DeinitFile(f1, &fh1);
    free(rt);
    free(et);
  }
  CloseFile(f1, &fh1);
  iuta = iu0;
  return 0;
}

int PreloadCE(char *tfn, char *sfn) {
  TFILE *f0, *f1;
  int ib, k, n, swp, iu0, iu1;
  F_HEADER fh, fh1;
  CE_HEADER h;
  CE_RECORD r;
  IDXDAT *di, *dj;
  
  if (_idxmap.nij == 0) {
    printf("index map not loaded\n");
    return -1;
  }
  iu0 = iuta;
  f0 = OpenFileRO(sfn, &fh, &swp);
  iu1 = iuta;
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
  iuta = iu0;
  f1 = OpenFile(tfn, &fh1);
  if (f1 == NULL) {
    printf("cannot open file %s\n", tfn);
    FCLOSE(f0);
    return -1;
  }

  for (ib = 0; ib < fh.nblocks; ib++) {
    iuta = iu1;
    n = ReadCEHeader(f0, &h, swp);
    if (n == 0) break;
    iuta = iu0;
    InitFile(f1, &fh1, &h);
    for (k = 0; k < h.ntransitions; k++) {
      iuta = iu1;
      n = ReadCERecord(f0, &r, swp, &h);
      if (n == 0) break;
      di = IdxMap(r.lower);
      dj = IdxMap(r.upper);
      if (!di || !dj) continue;
      r.lower = di->i;
      r.upper = dj->i;
      if (SetPreloadedCE(r.lower, r.upper)) {
	iuta = iu0;
	WriteCERecord(f1, &r);
      }
      if (h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }
    iuta = iu0;
    DeinitFile(f1, &fh1);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
  }
  FCLOSE(f0);
  CloseFile(f1, &fh1);
  return 0;
}

void SetCombEx(char *s) {
  int k, n, i;
  char *p, buf[2048];
  
  if (s == NULL || strlen(s) == 0) {
    for (k = 0; k < N_ELEMENTS1; k++) {
      if (s && _ncombex[k] > 0) free(_pcombex[k]);
      _scombex[k][0] = '\0';
      _ncombex[k] = 0;
      _pcombex[k] = NULL;
    }
    return;
  }
  strncpy(buf, s, 2048); 
  n = StrSplit(buf, ',');
  k = atoi(buf);
  if (k <= 0 || k > N_ELEMENTS) {
    printf("invalid combex option: %d %s\n", k, s);
    return;
  }
  memcpy(_scombex[k], buf, 2048);
  if (_ncombex[k] > 0) free(_pcombex[k]);
  if (n <= 1) {
    _scombex[k][0] = '\0';
    _ncombex[k] = 0;
    _pcombex[k] = NULL;
    return;
  }
  _ncombex[k] = n-1;
  _pcombex[k] = malloc(sizeof(char *)*_ncombex[k]);
  p = _scombex[k];
  for (i = 0; i < n; i++) {
    while(*p) p++;
    p++;
    if (i < n-1) {
      _pcombex[k][i] = p;
    }
  }
}

void SetOptionDBase(char *s, char *sp, int ip, double dp) {
  if (0 == strcmp(s, "dbase:combex")) {
    SetCombEx(sp);
    return;
  }
  if (0 == strcmp(s, "dbase:exk2h")) {
    _exk2h = ip;
    return;
  }
  if (0 == strcmp(s, "dbase:cmpetol")) {
    _cmpetol = dp;
    return;
  }
  if (0 == strcmp(s, "dbase:cmpnbm")) {
    _cmpnbm = ip;
    return;
  }  
  if (0 == strcmp(s, "dbase:inner_ai")) {
    SetInnerAI(sp);
    return;
  }
  if (0 == strcmp(s, "dbase:remove_closed")) {
    _remove_closed = ip;
    return;
  }
  if (0 == strcmp(s, "dbase:adj_ip")) {
    _adj_ip = ip;
    return;
  }
  if (0 == strcmp(s, "dbase:tuta")) {
    tuta = ip;
    return;
  }
  if (0 == strcmp(s, "dbase:nilast")) {
    _nilast = ip;
    return;
  }
  if (0 == strcmp(s, "dbase:collapse_mask")) {
    _collapse_mask = ip;
    return;
  }
  if (0 == strcmp(s, "dbase:grptol")) {
    _grptol = dp;
    return;
  }
  if (0 == strcmp(s, "dbase:sfu_smin")) {
    _sfu_smin = dp;
    return;
  }
  if (0 == strcmp(s, "dbase:sfu_dmax")) {
    _sfu_dmax = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp(s, "dbase:sfu_minde")) {
    _sfu_minde = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp(s, "dbase:sfu_maxde")) {
    _sfu_maxde = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp(s, "dbase:cwf")) {
    _cwf = dp;
    return;
  }
  if (0 == strcmp(s, "dbase:transmod")) {
    strncpy(_transmod_file, sp, 1024);
    return;
  }
  if (0 == strcmp(s, "dbase:ste")) {
    _ste = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp(s, "dbase:groupmod")) {
    strncpy(_groupmod_file, sp, 1024);
    return;
  }
}

