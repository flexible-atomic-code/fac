#include "crm.h"
#include "grid.h"
#include "cf77.h"

static char *rcsid="$Id: crm.c,v 1.92 2005/10/07 01:01:50 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static IONIZED ion0;
static ARRAY *ions;
static ARRAY *blocks;
static double *bmatrix = NULL;

static int n_single_blocks = 64;

static int rec_cascade = 0;
static double cas_accuracy = EPS4;
static int max_iter = 256;
static double iter_accuracy = EPS4;
static double iter_stabilizer = 0.75;

/* electron density in 10^10 cm-3 */
static double electron_density = EPS3;
/* photon energy density in erg cm-3 */
/* if the distribution is blackbody, then the normalization constant
 * is the dilution factor, (r0/r)^2 */
static double photon_density = 0.0;
static int ai_extra_nmax = 400;
static int do_extrapolate = 100;
static int inner_auger = 0;

int SetInnerAuger(int i) {
  inner_auger = i;
  return 0;
}

int SetExtrapolate(int e) {
  do_extrapolate = e;
  return 0;
}

int SetNumSingleBlocks(int n) {
  n_single_blocks = n;
  return 0;
}

int SetEleDensity(double ele) {
  if (ele >= 0.0) electron_density = ele;
  return 0;
}

int SetPhoDensity(double pho) {
  if (pho >= 0.0) photon_density = pho;
  return 0;
}

int SetIteration(double acc, double s, int max) {
  if (max >= 0) max_iter = max;
  if (acc > 0) iter_accuracy = acc;
  if (s > 0.0 && s < 1.0) iter_stabilizer = s;
  return 0;
}

int SetCascade(int c, double a) {
  rec_cascade = c;
  if (a > 0.0) cas_accuracy = a;
  return 0;
}

int InitCRM(void) {
  int i;

  for (i = 0; i < NDB; i++) ion0.dbfiles[i] = NULL;
  ion0.nionized = 0;
  ion0.energy = NULL;
  ion0.atom = 0;

  ions = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ions, sizeof(ION), ION_BLOCK);
  blocks = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(blocks, sizeof(LBLOCK), LBLOCK_BLOCK);
  bmatrix = NULL;
  
  InitDBase();
  InitRates();

  return 0;
}

static void InitBlkRateData(void *p, int n) {
  BLK_RATE *r;
  int k;

  r = (BLK_RATE *) p;
  for (k = 0; k < n; k++, r++) {
    r->rates = NULL;
  }
}

static void FreeBlkRateData(void *p) {
  BLK_RATE *r;

  r = (BLK_RATE *) p;
  ArrayFree(r->rates, NULL);
  if (r->rates) free(r->rates);
  r->rates = NULL;
}

static void InitIonData(void *p, int n) {
  ION *ion;
  int i, k;
  
  ion = (ION *) p;
  for (k = 0; k < n; k++,ion++) {
    ion->nlevels = 0;
    ion->KLN_min = 0;
    ion->KLN_max = -1;
    ion->KLN_bmin = 0;
    ion->KLN_bmax = -1;
    ion->KLN_amin = 0;
    ion->KLN_amax = -1;
    ion->nlevels = 0;
    ion->iblock = NULL;
    ion->ilev = NULL;
    ion->j = NULL;
    ion->vnl = NULL;
    ion->ibase = NULL;
    ion->energy = NULL;
    for (i = 0; i < NDB; i++) {
      ion->dbfiles[i] = NULL;
    }
    ion->ce_rates = NULL;
    ion->tr_rates = NULL;
    ion->tr_sdev = NULL;
    ion->tr2_rates = NULL;
    ion->ci_rates = NULL;
    ion->rr_rates = NULL;
    ion->ai_rates = NULL;
    ion->recombined = NULL;
  }
}
    
static void FreeIonData(void *p) {
  ION *ion;
  int i;

  ion = (ION *) p;
  if (ion->nlevels > 0) {
    free(ion->iblock);
    free(ion->ilev);
    free(ion->j);
    free(ion->vnl);
    free(ion->ibase);
    free(ion->energy);
    if (ion->KLN_bmax >= ion->KLN_bmin) {
      free(ion->KLN_ai);
      free(ion->KLN_nai);
    }
    ion->KLN_min = 0;
    ion->KLN_max = -1;
    ion->KLN_bmin = 0;
    ion->KLN_bmax = -1;
    ion->KLN_amin = 0;
    ion->KLN_amax = -1;
    ion->nlevels = 0;
  }
  for (i = 0; i < NDB; i++) {
    if (ion->dbfiles[i]) free(ion->dbfiles[i]);
    ion->dbfiles[i] = NULL;
  }
  ArrayFree(ion->ce_rates, FreeBlkRateData);
  free(ion->ce_rates);
  ion->ce_rates = NULL;
  ArrayFree(ion->tr_rates, FreeBlkRateData);
  free(ion->tr_rates);
  ion->tr_rates = NULL;
  ArrayFree(ion->tr_sdev, FreeBlkRateData);
  free(ion->tr_sdev);
  ion->tr_sdev = NULL;
  free(ion->tr2_rates);
  ion->tr2_rates = NULL;
  ArrayFree(ion->ci_rates, FreeBlkRateData);
  free(ion->ci_rates);
  ion->ci_rates = NULL;
  ArrayFree(ion->rr_rates, FreeBlkRateData);
  free(ion->rr_rates);
  ion->rr_rates = NULL;
  ArrayFree(ion->ai_rates, FreeBlkRateData);
  free(ion->ai_rates);
  ion->ai_rates = NULL;
  ArrayFree(ion->recombined, NULL);
  free(ion->recombined);
  ion->recombined = NULL;
}

static void InitBlockData(void *p, int n) {
  LBLOCK *blk;
  int k;

  blk = (LBLOCK *) p;
  for (k = 0; k < n; k++, blk++) {
    blk->nlevels = 0;
  }
}

static void FreeBlockData(void *p) {
  LBLOCK *blk;

  blk = (LBLOCK *) p;
  if (blk->nlevels > 0) {
    free(blk->n);
    free(blk->n0);
    free(blk->r);
    free(blk->total_rate);
  }
  blk->nlevels = 0;
}

int ReinitCRM(int m) {
  ION *ion;
  int i, k;

  if (m < 0) return 0;

  ReinitDBase(0);
  if (m == 3) return 0;
  
  if (m == 1) {
    for (k = 0; k < ions->dim; k++) {
      ion = (ION *) ArrayGet(ions, k);
      ArrayFree(ion->ce_rates, FreeBlkRateData);
      ArrayFree(ion->tr_rates, FreeBlkRateData);
      ArrayFree(ion->tr_sdev, FreeBlkRateData);
      ArrayFree(ion->tr2_rates, FreeBlkRateData);
      ArrayFree(ion->ci_rates, FreeBlkRateData);
      ArrayFree(ion->rr_rates, FreeBlkRateData);
      ArrayFree(ion->ai_rates, FreeBlkRateData);
    }
    return 0;
  } else if (m == 2) {
    for (k = 0; k < ions->dim; k++) {
      ion = (ION *) ArrayGet(ions, k);
      ArrayFree(ion->ce_rates, FreeBlkRateData);
      ArrayFree(ion->ci_rates, FreeBlkRateData);
      ArrayFree(ion->rr_rates, FreeBlkRateData);
      ArrayFree(ion->ai_rates, FreeBlkRateData);
    }
    return 0;
  }

  for (i = 0; i < NDB; i++) {
    if (ion0.dbfiles[i]) free(ion0.dbfiles[i]);
    ion0.dbfiles[i] = NULL;
  }
  if (ion0.nionized > 0) {
    free(ion0.energy);
    free(ion0.ionized_map[0]);
    free(ion0.ionized_map[1]);
    ion0.nionized = 0;
  }
  ion0.atom = 0;
  ArrayFree(ions, FreeIonData);
  ArrayFree(blocks, FreeBlockData);
  if (bmatrix) free(bmatrix);
  bmatrix = NULL;
  
  return 0;
}

int AddIon(int nele, double n, char *pref) {
  ION ion;
  int i;
  int m;
  
  ion.nlevels = 0;
  ion.ce_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.ce_rates, sizeof(BLK_RATE), RATES_BLOCK);
  ion.tr_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.tr_rates, sizeof(BLK_RATE), RATES_BLOCK);
  ion.tr_sdev = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.tr_sdev, sizeof(BLK_RATE), RATES_BLOCK);
  ion.tr2_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.tr2_rates, sizeof(BLK_RATE), RATES_BLOCK);
  ion.rr_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.rr_rates, sizeof(BLK_RATE), RATES_BLOCK);
  ion.ci_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.ci_rates, sizeof(BLK_RATE), RATES_BLOCK);
  ion.ai_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.ai_rates, sizeof(BLK_RATE), RATES_BLOCK);

  ion.KLN_min = 0;
  ion.KLN_max = -1;
  ion.KLN_bmin = 0;
  ion.KLN_bmax = -1;
  ion.KLN_amin = 0;
  ion.KLN_amax = -1;
  ion.KLN_ai = NULL;
  ion.KLN_nai = NULL;

  ion.recombined = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.recombined, sizeof(RECOMBINED), 16);
  
  ion.nele = nele;
  m = strlen(pref);
  m = m+4;
  for (i = 0; i < NDB; i++) {
    ion.dbfiles[i] = malloc(m);
    switch (i+1) {
    case DB_EN:
      sprintf(ion.dbfiles[i], "%s.en", pref);
      break;
    case DB_TR:
      sprintf(ion.dbfiles[i], "%s.tr", pref);
      break;
    case DB_CE:
      sprintf(ion.dbfiles[i], "%s.ce", pref);
      break;      
    case DB_RR:
      sprintf(ion.dbfiles[i], "%s.rr", pref);
      break;
    case DB_AI:
      sprintf(ion.dbfiles[i], "%s.ai", pref);
      break;
    case DB_CI:
      sprintf(ion.dbfiles[i], "%s.ci", pref);
      break;
    default:
      break;
    }
  }

  ion.n = n;

  ArrayAppend(ions, &ion, InitIonData);
  
  return ions->dim;
  
}

void GetRecombined(int *b, int *nrec, char *name) {
  int i;
  
  *nrec = 0;
  i = 0;
  while (name[i] && name[i] != '+') i++;
  if (name[i] == '+' && name[i-1] == ' ') {
    *b = atoi(name);
    *nrec = atoi(&(name[i+1]));
  }
}

void ExtrapolateEN(int iion, ION *ion) {
  RECOMBINED *rec;
  LBLOCK blk, *blkp;
  int i, j, nr, nlev;
  int n, n0, n1, nr0;
  int t, p, q, s, k, nc;
  double a, b, c, d;
  double a0, e0, delta, gamma;

  nlev = ion->nlevels;
  for (i = 0; i < ion->recombined->dim; i++) {
    if (i > do_extrapolate) break;
    rec = (RECOMBINED *) ArrayGet(ion->recombined, i);    
    j = rec->n-1;
    if (j > 0) {
      n1 = rec->nrec[j];
      n0 = rec->nrec[j-1]+1;
      nr = rec->imax[j] - rec->imin[j] + 1;
      k = n1 - n0;
      nlev += nr*k;
    }
  }
  ion->iblock = (LBLOCK **) realloc(ion->iblock, sizeof(LBLOCK *)*nlev);
  ion->ilev = (int *) realloc(ion->ilev, sizeof(int)*nlev);
  ion->j = (int *) realloc(ion->j, sizeof(int)*nlev);
  ion->vnl = (short *) realloc(ion->vnl, sizeof(short)*nlev);
  ion->ibase = (short *) realloc(ion->ibase, sizeof(short)*nlev);
  ion->energy = (double *) realloc(ion->energy, sizeof(double)*nlev);
  
  nr0 = ion->nlevels;
  c = ion0.atom - ion->nele + 1.0;
  c = 0.5*c*c;
  for (i = 0; i < ion->recombined->dim; i++) {
    if (i > do_extrapolate) break;
    rec = (RECOMBINED *) ArrayGet(ion->recombined, i);
    j = rec->n-1;
    rec->n_ext = rec->n;
    if (j > 0) {
      CopyNComplex(blk.ncomplex, ion->iblock[rec->imin[j]]->ncomplex);
      nc = 0;
      while (blk.ncomplex[nc+1].n > 0) nc++;
      n1 = rec->nrec[j];
      n0 = rec->nrec[j-1]+1;
      nr = rec->imax[j] - rec->imin[j] + 1;
      a = -c/(n1*n1);
      a0 = -c/((n0-1.0)*(n0-1.0));
      for (n = n0, t = rec->n; n < n1; n++, t++) {
	rec->nrec[t] = n;
	rec->imin[t] = nr0;
	nr0 += nr;
	rec->imax[t] = nr0-1;
	blk.ib = blocks->dim;
	blk.iion = iion;
	blk.nlevels = nr;
	blk.n = (double *) malloc(sizeof(double)*nr);
	blk.n0 = (double *) malloc(sizeof(double)*nr);
	blk.r = (double *) malloc(sizeof(double)*nr);
	blk.total_rate = (double *) malloc(sizeof(double)*nr);
	blk.rec = rec;
	blk.irec = t;
	blk.ncomplex[nc].n = n;
	blkp = ArrayAppend(blocks, &blk, InitBlockData);
	q = -1;
	p = rec->imin[t];
	s = rec->imin[j-1];
	d = -c/(n*n);
	for (k = rec->imin[j]; k <= rec->imax[j]; k++, p++, s++) {
	  q++;
	  ion->iblock[p] = blkp;
	  ion->ilev[p] = q;
	  ion->j[p] = ion->j[k];
	  ion->vnl[p] = n*100 + ion->vnl[k]%100;
	  ion->ibase[p] = ion->ibase[k];
	  e0 = ion->energy[ion->ibase[p]];
	  delta = ion->energy[k] - e0 - a;
	  ion->energy[p] = e0 + d + delta;
	  if (s <= rec->imax[j-1]) {
	    gamma = ion->energy[s] - e0 - a0 - delta;
	    b = gamma*(n1-n)/(n1-n0+1.0);
	    ion->energy[p] += b;
	  }
	}
      }
      rec->n_ext = t;
    }
  }
  ion->nlevels = nlev;
}
  
void ExtrapolateTR(ION *ion, int inv) {
  RECOMBINED *rec;
  RATE *r, *rp, r0;
  ARRAY *rates, *rates0;
  LBLOCK *blk, *blk0;
  BLK_RATE *brts, *brts0;
  int nr;
  int n0, n1, n;
  int i, j, t, k, p, q, s;
  int imin, imax, iuta;
  double a, b, c, h;

  iuta = IsUTA();
  for (i = 0; i < ion->recombined->dim; i++) {
    if (i > do_extrapolate) break;
    rec = (RECOMBINED *) ArrayGet(ion->recombined, i);
    if (rec->n_ext == rec->n) continue;
    j = rec->n-1;
    imin = rec->imin[j];
    imax = rec->imax[j];
    blk0 = ion->iblock[imin];
    n1 = rec->nrec[j];
    a = (double) n1;
    a = a*a*a;
    k = j-1;
    n = rec->nrec[k];
    c = (double ) n;
    c = a/(c*c*c);
    for (p = 0; p < ion->tr_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr_rates, p);
      if (brts->iblock != blk0) continue;
      for (k = 0; k < ion->tr_rates->dim; k++) {
	brts0 = (BLK_RATE *) ArrayGet(ion->tr_rates, k);
	if (brts0->fblock == brts->fblock &&
	    brts0->iblock->ib == blk0->ib - 1) {
	  rates0 = brts0->rates;
	  break;
	}
      }
      blk = brts->fblock;
      rates = brts->rates;
      nr = rates->dim;
      for (t = 0; t < nr; t++) {
	r = (RATE *) ArrayGet(rates, t);
	q = -1;
	if (blk->rec) {
	  if (blk->rec->nrec[blk->irec] == n1) {
	    q = r->f - blk->rec->imin[blk->irec];
	  }
	} 
	if (q < 0) {
	  if (t < rates0->dim) {
	    rp = (RATE *) ArrayGet(rates0, t);
	    h = rp->dir/(r->dir*c);
	  } else {
	    h = 1.0;
	  }
	}
	for (k = rec->n; k < rec->n_ext; k++) {
	  n0 = rec->nrec[k];
	  r0.i = r->i - imin + rec->imin[k];
	  if (q < 0) {
	    r0.f = r->f;
	    b = 1.0/((double) n0);
	    b = b*b*b*a;
	    b *= h + (n0-n)*(1.0-h)/(n1-n);
	    r0.dir = b*r->dir;
	    r0.inv = b*r->inv;
	    AddRate(ion, ion->tr_rates, &r0, 0);
	    if (iuta) {
	      r0.dir = ion->energy[r0.i] - ion->energy[r0.f];
	      r0.inv = 0.0;
	      AddRate(ion, ion->tr_sdev, &r0, 0);
	    }
	  } else {
	    for (s = 0; s < blk->rec->n_ext; s++) {
	      if (blk->rec->nrec[s] == n0) {
		r0.f = q + blk->rec->imin[s];
		if (r0.f > blk->rec->imax[s]) break;
		r0.dir = r->dir;
		r0.inv = r->inv;
		AddRate(ion, ion->tr_rates, &r0, 0);
		if (iuta) {
		  r0.dir = ion->energy[r0.i] - ion->energy[r0.f];
		  r0.inv = 0.0;
		  AddRate(ion, ion->tr_sdev, &r0, 0);
		}
		break;
	      }
	    }
	  }
	}
      }
    }
  }
}

void ExtrapolateRR(ION *ion, int inv) {
  RECOMBINED *rec;
  RATE *r, *rp, r0;
  ARRAY *rates, *rates0;
  LBLOCK *blk0;
  BLK_RATE *brts, *brts0;
  DISTRIBUTION *dist;
  int nr;
  int n0, n1;
  int i, j, t, k, p;
  int imin, imax;
  double a, b, c, z, temp, h, d;
  double rr_extra[MAXNREC];

  dist = GetEleDist(&i);
  if (i != 0) return;
  temp = dist->params[0];
  z = ion0.atom - ion->nele + 1.0;
  for (i = 0; i < ion->recombined->dim; i++) {
    if (i > do_extrapolate) break;
    rec = (RECOMBINED *) ArrayGet(ion->recombined, i);
    if (rec->n_ext == rec->n) {
    }
    j = rec->n - 1;
    imin = rec->imin[j];
    imax = rec->imax[j];
    blk0 = ion->iblock[imin];
    n1 = rec->nrec[j];
    a = RRRateHydrogenic(temp, z, n1, &b);
    c = 0.0;
    for (p = 0; p < ion->rr_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->rr_rates, p);
      if (brts->fblock != blk0) continue;
      rates = brts->rates;
      nr = rates->dim;
      for (t = 0; t < nr; t++) {
	r = (RATE *) ArrayGet(rates, t);
	if (r->i == rec->bmin) {
	  c += r->dir;
	}
      }
    }
    c = a/c;
    rr_extra[j] = b/a;
    for (k = rec->n; k < rec->n_ext; k++) {
      n0 = rec->nrec[k];
      b = RRRateHydrogenic(temp, z, n0, NULL);
      rr_extra[k] = b/a;
    }
    k = j-1;
    n0 = rec->nrec[k];
    b = RRRateHydrogenic(temp, z, n0, NULL);
    rr_extra[k] = b/a;

    a = rr_extra[j];
    for (p = 0; p < ion->rr_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->rr_rates, p);
      if (brts->fblock != blk0) continue;
      rates = brts->rates;
      for (k = 0; k < ion->rr_rates->dim; k++) {
	brts0 = (BLK_RATE *) ArrayGet(ion->rr_rates, k);
	if (brts0->iblock == brts->iblock &&
	    brts0->fblock->ib == blk0->ib - 1) {
	  rates0 = brts0->rates;
	  break;
	}
      }
      nr = rates->dim;
      for (t = 0; t < nr; t++) {
	r = (RATE *) ArrayGet(rates, t);
	r->dir *= c;
	if (t < rates0->dim) {
	  rp = (RATE *) ArrayGet(rates0, t);
	  h = rp->dir/(r->dir*rr_extra[j-1]);
	} else {
	  h = 1.0;
	}
	for (k = rec->n; k < rec->n_ext; k++) {
	  r0.i = r->i;
	  r0.f = r->f - imin + rec->imin[k];
	  b = rr_extra[k];
	  d = (double)(rec->nrec[k]-n0)/(n1-n0);
	  b *= h + d*(1.0-h);
	  r0.dir = b*r->dir;
	  r0.inv = 0.0;
	  AddRate(ion, ion->rr_rates, &r0, 0);
	}
	r->dir *= a;
      }
    }
  }
}

void ExtrapolateAI(ION *ion, int inv) {
  RECOMBINED *rec;
  RATE *r, r0;
  LBLOCK *blk0;
  BLK_RATE *brts;
  ARRAY *rates;
  int nr;
  int n0, n1;
  int i, j, t, k, p;
  int imin, imax;
  double a, b, c, e;
  double ai_extra[MAXNREC];

  for (i = 0; i < ion->recombined->dim; i++) {
    if (i > do_extrapolate) break;
    rec = (RECOMBINED *) ArrayGet(ion->recombined, i);
    j = rec->n - 1;
    imin = rec->imin[j];
    imax = rec->imax[j];
    blk0 = ion->iblock[imin];
    n1 = rec->nrec[j];
    a = 1.0/n1;
    a = a*a*a;
    for (k = rec->n; k < rec->n_ext; k++) {
      n0 = rec->nrec[k];
      b = 1.0/n0;
      b = b*b*b;
      ai_extra[k] = b/a;
    }
    ai_extra[j] = 0.0;
    c = 0.0;
    for (k = rec->nrec[j]+1; k <= ai_extra_nmax; k++) {
      b = 1.0/k;
      b = b*b*b;
      ai_extra[j] += b;
      c += b*b;
    }
    c = (a + c/a)/(a + ai_extra[j]);
    ai_extra[j] = 1.0 + ai_extra[j]/a;
    a = ai_extra[j];
    for (p = 0; p < ion->ai_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->ai_rates, p);
      if (brts->iblock != blk0) continue;
      rates = brts->rates;
      nr = rates->dim;
      for (t = 0; t < nr; t++) {
	r = (RATE *) ArrayGet(rates, t);
	for (k = rec->n; k < rec->n_ext; k++) {
	  r0.i = r->i - imin + rec->imin[k];
	  r0.f = r->f;
	  e = ion->energy[r0.i] - ion->energy[r0.f];
	  if (e < EPS16) continue;
	  b = ai_extra[k];
	  r0.dir = b*r->dir;
	  AIRate(&(r0.dir), &(r0.inv), inv, ion->j[r->i], ion->j[r->f],
		 e, r0.dir/RATE_AU);
	  AddRate(ion, ion->ai_rates, &r0, 0);
	}
	r->inv *= a;
	r->dir *= c;
      }
    }  
  }
}
    
int SetBlocks(double ni, char *ifn) {
  ION *ion, *ion1 = NULL;
  F_HEADER fh;
  EN_HEADER h;
  EN_RECORD r, *r0, *r1;
  EN_RECORD *rionized;
  LBLOCK blk, *blkp;
  RECOMBINED *rec, rec0;
  NCOMPLEX ncomplex[MAXNCOMPLEX];
  int bmin, bmax, imin, imax, t, nrec;
  int ibase, tbase;
  FILE *f;
  int n, i, k, nb, nb0, nlevels;
  char *fn;
  int p, q = -1;
  int nionized, n0;
  int swp, sfh;

  ion0.n = ni;
  ion0.n0 = ni;
  if (ifn) {
    k = strlen(ifn);
    k += 4;
  } else {
    k = 0;
  }
  for (i = 0; i < NDB; i++) {
    if (k > 0) {
      switch (i+1) {
      case DB_EN:
	ion0.dbfiles[i] = (char *) malloc(k);
	sprintf(ion0.dbfiles[i], "%s.en", ifn);
	break;
      case DB_TR:
	ion0.dbfiles[i] = (char *) malloc(k);
	sprintf(ion0.dbfiles[i], "%s.tr", ifn);
	break;
      case DB_CE:
	ion0.dbfiles[i] = (char *) malloc(k);
	sprintf(ion0.dbfiles[i], "%s.ce", ifn);
	break;
      case DB_AI:
	ion0.dbfiles[i] = (char *) malloc(k);
	sprintf(ion0.dbfiles[i], "%s.ai", ifn);
	break;
      default:
	ion0.dbfiles[i] = NULL;
      }
    } else {
      ion0.dbfiles[i] = NULL;
    }
  }
  
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ion->n0 = ion->n;
    if (k > 0) {
      if (ion->nele != ion1->nele+1) {
	printf("ERROR: NELE for the added ions are not ");
	printf("a continuous ascending sequence\n");
	exit(1);
      }
      ifn = ion1->dbfiles[DB_EN-1];
    } else {
      ion0.nele = ion->nele - 1;
      ifn = ion0.dbfiles[DB_EN-1];
    }

    fn = ion->dbfiles[DB_EN-1];
    f = fopen(fn, "r");
    if (f == NULL) {
      printf("File %s does not exist\n", fn);
      return -1;
    }
    n = ReadFHeader(f, &fh, &swp);
    if (VersionLE((&fh), 1, 0, 8)) sfh = sizeof(F_HEADER);
    else sfh = SIZE_F_HEADER;

    if (k == 0) {
      ion0.atom = fh.atom;
      strcpy(ion0.symbol, fh.symbol);
    }

    nlevels = 0;
    nionized = 0;
    for (i = 0; i < fh.nblocks; i++) {
      n = ReadENHeader(f, &h, swp);
      nlevels += h.nlevels;
      if (h.nele == ion->nele-1) {
	nionized += h.nlevels;
      }
      fseek(f, h.length, SEEK_CUR);
    }
    ion->nlevels = nlevels;
    ion->iblock = (LBLOCK **) malloc(sizeof(LBLOCK *)*nlevels);
    for (i = 0; i < nlevels; i++) ion->iblock[i] = NULL;
    ion->ilev = (int *) malloc(sizeof(int)*nlevels);
    ion->j = (int *) malloc(sizeof(int)*nlevels);
    ion->vnl = (short *) malloc(sizeof(short)*nlevels);
    ion->ibase = (short *) malloc(sizeof(short)*nlevels);
    ion->energy = (double *) malloc(sizeof(double)*nlevels);
    rionized = (EN_RECORD *) malloc(sizeof(EN_RECORD )*nionized);
    if (k == 0 && ifn) {
      ion0.nionized = nionized;
      ion0.ionized_map[0] = (int *) malloc(sizeof(int)*nionized);
      ion0.ionized_map[1] = (int *) malloc(sizeof(int)*nionized);
      ion0.energy = (double *) malloc(sizeof(double)*nionized);
    }
    
    fseek(f, sfh, SEEK_SET);
    n0 = 0;
    nb0 = 0;
    r0 = rionized;
    if (k == 0) {
      ion0.imin[0] = 100000000;
      ion0.imin[1] = 100000000;
      ion0.imax[0] = 0;
      ion0.imax[1] = 0;
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadENHeader(f, &h, swp);
      if (h.nele == ion->nele) {
	fseek(f, h.length, SEEK_CUR);
	continue;
      } else if (h.nele == ion->nele-1) {
	for (i = 0; i < h.nlevels; i++) {
	  n = ReadENRecord(f, &r0[i], swp);
	}
	if (inner_auger) {
	  if (ion->nele >= 4 && ion->nele <= 10) {
	    GetNComplex(ncomplex, r0[0].ncomplex);
	    if (ncomplex[0].n == 1) {
	      if (ncomplex[0].nq == 1 &&
		  ncomplex[1].n == 2 &&
		  ncomplex[1].nq == ion->nele-2) {
		ion->KLN_bmin = r0[0].ilev;
		ion->KLN_bmax = r0[h.nlevels-1].ilev;
		ion->KLN_ai = (double *) malloc(sizeof(double)*h.nlevels);
		ion->KLN_nai = (int *) malloc(sizeof(int)*h.nlevels);
		for (ibase = 0; ibase < h.nlevels; ibase++) {
		  ion->KLN_ai[ibase] = 0.0;
		  ion->KLN_nai[ibase] = 0;
		}
	      } else if (ncomplex[0].nq == 2) {
		if (ncomplex[1].n > 2 ||
		    ncomplex[1].nq == ion->nele-4) {
		  ion->KLN_amin = r0[0].ilev;
		  ion->KLN_amax = r0[h.nlevels-1].ilev;
		}
	      }       
	    }
	  }
	  
	  if (ion->nele >= 12) {
	    GetNComplex(ncomplex, r0[0].ncomplex);
	    if (ncomplex[0].n == 1 && 
		ncomplex[0].nq == 2 && 
		ncomplex[1].n == 2) {
	      if (ncomplex[1].nq == 7 &&
		  ncomplex[2].n == 3 &&
		  ncomplex[2].nq == ion->nele-10) {
		ion->KLN_bmin = r0[0].ilev;
		ion->KLN_bmax= r0[h.nlevels-1].ilev;
		ion->KLN_ai = (double *) malloc(sizeof(double)*h.nlevels);
		ion->KLN_nai = (int *) malloc(sizeof(int)*h.nlevels);
		for (ibase = 0; ibase < h.nlevels; ibase++) {
		  ion->KLN_ai[ibase] = 0.0;
		  ion->KLN_nai[ibase] = 0;
		}
	      } else if (ncomplex[1].nq == 8) {
		if (ncomplex[2].n > 3 ||
		    ncomplex[2].nq == ion->nele-12) {
		  ion->KLN_amin = r0[0].ilev;
		  ion->KLN_amax = r0[h.nlevels-1].ilev;
		}
	      }       
	    }
	  }
	}

	if (ifn) {
	  r1 = (EN_RECORD *) malloc(sizeof(EN_RECORD)*h.nlevels);
	  nlevels = FindLevelBlock(h.nlevels, r0, r1, ion->nele-1, ifn); 
	  if (nlevels != h.nlevels) {
	    printf("ERROR: Ionized block %d of ion %d ", nb, ion->nele);
	    printf("does not match a block in file %s\n", ifn);
	    printf("nlevels = %d VS %d\n", h.nlevels, nlevels);
	    exit(1);
	  }
	}
	if (k > 0) {
	  for (i = 0; i < h.nlevels; i++) {
	    p = r0[i].ilev;
	    q = r1[i].ilev;
	    ion1->iblock[q]->ionized = 1;
	    ion->iblock[p] = ion1->iblock[q];
	    ion->ilev[p] = ion1->ilev[q];
	    ion->j[p] = JFromENRecord(&(r0[i]));
	    if (r0[i].p < 0) {
	      ion->vnl[p] = -r0[i].p;
	    } else {
	      ion->vnl[p] = r0[i].p;
	    }
	    ion->ibase[p] = -1;
	    ion->energy[p] = r0[i].energy;
	  }
	} else {
	  blk.ncomplex[0].n = 0;
	  blk.nlevels = 0;
	  blkp = NULL;
	  for (i = 0; i < h.nlevels; i++) {
	    GetNComplex(ncomplex, r0[i].ncomplex);
	    if (nb0 == 0 && i <= n_single_blocks) {
	      nlevels = 0;
	      blk.ib = blocks->dim;
	      blk.iion = -1;
	      blk.irec = -1;
	      blk.ionized = 1;
	      blk.rec = NULL;
	      if (i < n_single_blocks) {
		blk.nlevels = 1;
	      } else {
		blk.nlevels = h.nlevels - n_single_blocks;
	      }      
	      blk.n = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.n0 = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.r = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.total_rate = (double *) malloc(sizeof(double)*blk.nlevels);
	      CopyNComplex(blk.ncomplex, ncomplex);
	      blkp = ArrayAppend(blocks, &blk, InitBlockData);
	      q = -1;
	    } else if (CompareNComplex(ncomplex, blk.ncomplex)) {
	      if (blkp) {
		if (blkp->nlevels > nlevels) {
		  blkp->n = (double *) ReallocNew(blkp->n, 
						  sizeof(double)*nlevels);
		  blkp->n0 = (double *) ReallocNew(blkp->n0, 
						   sizeof(double)*nlevels);
		  blkp->r = (double *) ReallocNew(blkp->r, 
					       sizeof(double)*nlevels);
		  blkp->total_rate = (double *)ReallocNew(blkp->total_rate,
						       sizeof(double)*nlevels);
		  blkp->nlevels = nlevels;
		}
	      }
	      nlevels = 0;
	      blk.ib = blocks->dim;
	      blk.iion = -1;
	      blk.irec = -1;
	      blk.ionized = 1;
	      blk.rec = NULL;
	      blk.nlevels = h.nlevels;
	      blk.n = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.n0 = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.r = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.total_rate = (double *) malloc(sizeof(double)*blk.nlevels);
	      CopyNComplex(blk.ncomplex, ncomplex);
	      blkp = ArrayAppend(blocks, &blk, InitBlockData);
	      q = -1;
	    }
	    p = r0[i].ilev;
	    q++;
	    nlevels++;
	    ion->iblock[p] = blkp;
	    ion->ilev[p] = q;
	    ion->j[p] = JFromENRecord(&(r0[i]));
	    if (r0[i].p < 0) {
	      ion->vnl[p] = -r0[i].p;
	    } else {
	      ion->vnl[p] = r0[i].p;
	    }
	    ion->ibase[p] = -1;
	    ion->energy[p] = r0[i].energy;
	    if (ifn) {
	      if (r1[i].ilev > ion0.imax[0]) ion0.imax[0] = r1[i].ilev;
	      if (r0[i].ilev > ion0.imax[1]) ion0.imax[1] = r0[i].ilev;
	      if (r1[i].ilev < ion0.imin[0]) ion0.imin[0] = r1[i].ilev;
	      if (r0[i].ilev < ion0.imin[1]) ion0.imin[1] = r0[i].ilev;
	      ion0.ionized_map[0][n0] = r1[i].ilev;
	      ion0.ionized_map[1][n0] = r0[i].ilev;
	      ion0.energy[n0] = r1[i].energy;
	      n0++;
	    }
	  }
	}
	if (nb0 == 0) ion->iground = r0[0].ilev;
	if (ifn) {
	  free(r1);
	}
	nb0++;
	r0 += h.nlevels;
      } else if (h.nele == ion->nele-2) {
	fseek(f, h.length, SEEK_CUR);
	continue;	
      } else {
	printf("ERROR: Ion charge state does not match %d %d %d %d\n",
	       k, nb, h.nele, ion->nele);
	exit(1);
      }
    }

    fseek(f, sfh, SEEK_SET);
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadENHeader(f, &h, swp);
      if (h.nele != ion->nele) {
	fseek(f, h.length, SEEK_CUR);
	continue;
      }	
      blk.ncomplex[0].n = 0;
      blkp = NULL;
      nlevels = 0;
      for (i = 0; i < h.nlevels; i++) {
	n = ReadENRecord(f, &r, swp);
	GetNComplex(ncomplex, r.ncomplex);
	if (nb == 0 && i <= n_single_blocks) {
	  nlevels = 0;
	  blk.ib = blocks->dim;
	  blk.iion = k;
	  blk.irec = -1;
	  blk.ionized = 0;
	  blk.rec = NULL;
	  if (i < n_single_blocks) {
	    blk.nlevels = 1;
	  } else {
	    blk.nlevels = h.nlevels - n_single_blocks;
	  }
	  blk.n = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.n0 = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.r = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.total_rate = (double *) malloc(sizeof(double)*blk.nlevels);
	  CopyNComplex(blk.ncomplex, ncomplex);
	  blkp = ArrayAppend(blocks, &blk, InitBlockData);
	  q = -1;
	} else if (CompareNComplex(ncomplex, blk.ncomplex)) {
	  if (blkp) {
	    if (blkp->nlevels > nlevels) {
	      blkp->n = (double *) ReallocNew(blkp->n, 
					      sizeof(double)*nlevels);
	      blkp->n0 = (double *) ReallocNew(blkp->n0, 
					       sizeof(double)*nlevels);
	      blkp->r = (double *) ReallocNew(blkp->r, 
					      sizeof(double)*nlevels);
	      blkp->total_rate = (double *) ReallocNew(blkp->total_rate,
						       sizeof(double)*nlevels);
	      blkp->nlevels = nlevels;
	    }
	  }
	  nlevels = 0;
	  blk.ib = blocks->dim;
	  blk.iion = k;
	  blk.irec = -1;
	  blk.ionized = 0;
	  blk.rec = NULL;
	  blk.nlevels = h.nlevels;
	  blk.n = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.n0 = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.r = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.total_rate = (double *) malloc(sizeof(double)*blk.nlevels);
	  CopyNComplex(blk.ncomplex, ncomplex);
	  blkp = ArrayAppend(blocks, &blk, InitBlockData);
	  q = -1;
	}
	
	if (i == 0) {
	  GetRecombined(&t, &nrec, r.name);
	  if (nrec > 0) {
	    imin = r.ilev;
	    bmin = t;
	  }
	} else if (nrec > 0 && i == h.nlevels-1) {
	  GetRecombined(&t, &nrec, r.name);
	  bmax = t;
	  imax = r.ilev;
	  for (t = 0; t < ion->recombined->dim; t++) {
	    rec = (RECOMBINED *) ArrayGet(ion->recombined, t);
	    if (rec->bmin == bmin && rec->bmax == bmax) {
	      rec->imin[rec->n] = imin;
	      rec->imax[rec->n] = imax;
	      rec->nrec[rec->n] = nrec;
	      blkp->irec = rec->n;
	      blkp->rec = rec;
	      rec->n++;
	      break;
	    }
	  }
	  if (t == ion->recombined->dim) {
	    rec0.n = 1;
	    rec0.bmin = bmin;
	    rec0.bmax = bmax;
	    rec0.imin[0] = imin;
	    rec0.imax[0] = imax;
	    rec0.nrec[0] = nrec;
	    blkp->rec = ArrayAppend(ion->recombined, &rec0, NULL);
	    blkp->irec = 0;
	  }
	}
	p = r.ilev;
	q++;
	nlevels++;
	ion->iblock[p] = blkp;
	ion->ilev[p] = q;
	ion->j[p] = JFromENRecord(&(r));
	if (r.p < 0) {
	  ion->vnl[p] = -r.p;
	} else {
	  ion->vnl[p] = r.p;
	}
	ion->ibase[p] = IBaseFromENRecord(&r);
	ion->energy[p] = r.energy;
	if (ion->nele >= 4 && nrec == 0) {
	  if (ion->ibase[p] <= ion->KLN_bmax &&
	      ion->ibase[p] >= ion->KLN_bmin) {
	    ibase = ion->ibase[p] - ion->KLN_bmin;
	    if (i == 0) { 
	      ion->KLN_min = r.ilev;
	      for (tbase = ion->KLN_bmin; tbase <= ion->KLN_bmax; tbase++) {
		ion->KLN_nai[tbase-ion->KLN_bmin] = 0;
	      }
	      ion->KLN_nai[ibase]++;
	    } else if (i == h.nlevels-1) {
	      ion->KLN_max = r.ilev;
	      ion->KLN_nai[ibase]++;
	    } else {
	      ion->KLN_nai[ibase]++;
	    }
	  }
	}
      }
    }
    if (ion0.n >= 0.0) {
      ExtrapolateEN(k, ion);
    }
    ion1 = ion;
    free(rionized);
    fclose(f);
    
    /* determine the minimum ilev in each block */
    blkp = NULL;
    for (i = 0; i < ion->nlevels; i++) {
      if (ion->iblock[i] == NULL) continue;
      if (ion->iblock[i] != blkp) {
	blkp = ion->iblock[i];
	blkp->imin = i;
      }
    }
  }
  
  k = blocks->dim;
  if (bmatrix) free(bmatrix);
  if (k > 0) {
    k = 2*k*(k+1);
    bmatrix = (double *) malloc(sizeof(double)*k);
  }
  
  return 0;
}

static int CompareENRecord(const void *p0, const void *p1) {
  EN_RECORD *r0, *r1;
  
  r0 = (EN_RECORD *) p0;
  r1 = (EN_RECORD *) p1;
  
  if (r0->j < r1->j) {
    return -1;
  } else if (r0->j > r1->j) {
    return 1;
  } else {
    if (r0->p < r1->p) {
      return -1;
    } else if (r0->p > r1->p) {
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

int FindLevelBlock(int n, EN_RECORD *r0, EN_RECORD *r1, 
		   int nele, char *ifn) {
  F_HEADER fh;
  EN_HEADER h;
  EN_RECORD g;
  FILE *f;
  int i, k, j, nr, nb;
  int swp, sfh;

  f = fopen(ifn, "r");
  if (f == NULL) {
    printf("File %s does not exist\n", ifn);
    return -1;
  }
  
  nr = ReadFHeader(f, &fh, &swp);
  if (VersionLE((&fh), 1, 0, 8)) sfh = sizeof(F_HEADER);
  else sfh = SIZE_F_HEADER;

  k = 0;
  for (nb = 0; nb < fh.nblocks; nb++) {
    nr = ReadENHeader(f, &h, swp);
    if (h.nele != nele) {
      fseek(f, h.length, SEEK_CUR);
      continue;
    }
    for (i = 0; i < h.nlevels; i++) {
      nr = ReadENRecord(f, &r1[k], swp);
      if (strcmp(r1[k].ncomplex, r0[0].ncomplex) == 0) {
	k++;
	if (k == n) break;
      }
    }
    if (k == n) break;
  }

  if (k < n) {
    fseek(f, sfh, SEEK_SET);
    nr = ReadENHeader(f, &h, swp);
    j = h.nlevels;
    for (i = 0; i < h.nlevels; i++) {
      nr = ReadENRecord(f, &r1[k], swp);
      if (strcmp(r1[k].ncomplex, r0[0].ncomplex) != 0) j--;
      if (j < n) {
	k++;
	if (k == n) break;
      }
    }
  }

  if (k != n) return k;

  memcpy(&g, r0, sizeof(EN_RECORD));
  qsort(r0, n, sizeof(EN_RECORD), CompareENRecord);
  qsort(r1, k, sizeof(EN_RECORD), CompareENRecord);
  for (i = 0; i < n; i++) {
    if (g.ilev == r0[i].ilev) {
      if (i != 0) {
	memcpy(&(r0[i]), r0, sizeof(EN_RECORD));
	memcpy(r0, &g, sizeof(EN_RECORD));
	memcpy(&g, &(r1[i]), sizeof(EN_RECORD));
	memcpy(&(r1[i]), r1, sizeof(EN_RECORD));
	memcpy(r1, &g, sizeof(EN_RECORD));
      }
      break;
    }
  }

  fclose(f);
  return n;
}

int IonIndex(ION *ion, int i, int k) {
  int m;

  m = 0;
  while (m < ion->nlevels) {
    if (ion->iblock[m] == NULL) {
      m++;
    } else if (ion->iblock[m]->ib != i) {
      m += ion->iblock[m]->nlevels;
    } else {
      if (ion->ilev[m] == k) return m;
      m++;
    }
  }
  return -1;
}

int IonizedIndex(int i, int m) {
  int k;

  if (i > ion0.imax[m] || i < ion0.imin[m]) return -1;

  for (k = 0; k < ion0.nionized; k++) {
    if (ion0.ionized_map[m][k] == i) {
      return k;
    }
  }
  
  return -1;
}

int GetNComplex(NCOMPLEX *c, char *s) {
  int i, n, nq;
  char *p;
  
  i = 0;
  while (1) {
    if (i == MAXNCOMPLEX-1) {
      printf("Num of NCOMPLEX shells exceeded the limit %d\n", MAXNCOMPLEX-1);
      exit(1);
    }
    n = strtol(s, &p, 10);
    if (n == 0) {
      for (; i < MAXNCOMPLEX; i++) {
	c[i].n = 0;
	c[i].nq = 0;
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
}

int StrNComplex(char *s, NCOMPLEX *c) {
  int i;
  char a[8];
  
  i = 0;
  s[0] = '\0';
  while (i < MAXNCOMPLEX && c[i].n) {
    sprintf(a, "%d*%d ", c[i].n, c[i].nq);
    strcat(s, a);
    i++;
  }
  
  return i;
}

int CompareNComplex(NCOMPLEX *c1, NCOMPLEX *c2) {
  int i;

  i = 0;
  while (c1[i].n && c2[i].n) {
    if (c1[i].n == c2[i].n && c1[i].nq == c2[i].nq) i++;
    else return 1;
  }
  if (c1[i].n || c2[i].n) return 1;

  return 0;
}

int CopyNComplex(NCOMPLEX *d, NCOMPLEX *s) {
  int i;
  for (i = 0; i < MAXNCOMPLEX; i++) {
    d[i].n = s[i].n;
    d[i].nq = s[i].nq;
  }
  return 0;
}

int TransitionType(NCOMPLEX *ic, NCOMPLEX *fc) {
  int n1[1024];
  int n2[1024];
  int i, k, k1, k2, m1, m2;

  i = 0;
  for (i = 0; i < 1024; i++) {
    n1[i] = 0;
    n2[i] = 0;
  }
  i = 0;
  while (ic[i].n) {
    n1[ic[i].n] = ic[i].nq;
    i++;
  }
  m1 = ic[i-1].n;

  i = 0;
  while (fc[i].n) {
    n2[fc[i].n] = fc[i].nq;
    i++;
  }
  m2 = fc[i-1].n;


  k1 = 0;
  k2 = 0;
  for (i = 0; i < 1024; i++) {
    k = n1[i] - n2[i];
    if (k == 0) continue;
    else {
      if (k == -1) {
	if (k2) return -1;
	k2 = i;
      }
      else if (k == 1) {
	if (k1) return -1;
	k1 = i;
      }
      else {
	return -1;
      }
    }
  }

  if (k1 == 0) {
    if (k2 == 0) return m1*100+m2;
    else return k2;
  }
  if (m1 == m2) {
    if (k1 < m1 || k2 < m1) {
      return m1*10000+k1*100+k2;
    }
  } else {
    return k1*100+k2;
  }
}

int SetAbund(int nele, double abund) {
  ION *ion;
  int i;

  if (ion0.nele == nele) {
    ion0.n = abund;
    ion0.n0 = abund;
  } else {
    for (i = 0; i < ions->dim; i++) {
      ion = (ION *) ArrayGet(ions, i);
      if (ion->nele == nele) {
	ion->n = abund;
	ion->n0 = abund;
	break;
      }
    }
  }
  
  return 0;
}

int InitBlocks(void) {
  ION  *ion;
  RATE *r;
  BLK_RATE *brts;
  LBLOCK *blk1, *blk2;
  int k, m, i, j, p;
  double a, b;
 
  for (i = 0; i < blocks->dim; i++) {
    blk1 = (LBLOCK *) ArrayGet(blocks, i);
    blk1->nb = 1.0;
    for (k = 0; k < blk1->nlevels; k++) {
      blk1->n0[k] = 0.0;
      blk1->n[k] = 0.0;
      blk1->total_rate[k] = 0.0;
      blk1->r[k] = 0.0;
    }
    blk1->r[0] = 1.0;
  }

  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    if (ion->nele >= 4) {
      for (i = 0; i < ion->nlevels; i++) {
	blk1 = ion->iblock[i];
	if (blk1 == NULL) continue;
	if (ion->ibase[i] <= ion->KLN_bmax &&
	    ion->ibase[i] >= ion->KLN_bmin) {
	  p = ion->ibase[i] - ion->KLN_bmin;
	  if (blk1->irec >= 0) {
	    blk1->total_rate[ion->ilev[i]] = ion->KLN_ai[p];
	  }
	}
      }
    }
    if (electron_density > 0.0) {
      for (p = 0; p < ion->ce_rates->dim; p++) {
	brts = (BLK_RATE *) ArrayGet(ion->ce_rates, p);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	for (m = 0; m < brts->rates->dim; m++) {
	  r = (RATE *) ArrayGet(brts->rates, m);
	  j = ion->ilev[r->i];
	  blk1->total_rate[j] += electron_density * r->dir;
	  if (r->inv > 0.0) {
	    j = ion->ilev[r->f];
	    blk2->total_rate[j] += electron_density * r->inv;
	  }
	}
      }
    }
    for (p = 0; p < ion->tr_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr_rates, p);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	j = ion->ilev[r->i];
	blk1->total_rate[j] += r->dir;
	blk1->n[j] += r->dir;
	if (r->inv > 0.0 && photon_density > 0.0) {
	  a = photon_density * r->inv;
	  b = a * (ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
	  blk1->total_rate[j] += b;
	  j = ion->ilev[r->f];
	  blk2->total_rate[j] += a;
	}
      }
    }
    for (p = 0; p < ion->tr2_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr2_rates, p);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	j = ion->ilev[r->i];
	blk1->total_rate[j] += r->dir;
	blk1->n[j] += r->dir;
      }
    }
    for (p = 0; p < ion->rr_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->rr_rates, p);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	j = ion->ilev[r->i];
	if (electron_density > 0.0) {
	  blk1->total_rate[j] += electron_density * r->dir;
	}
	if (r->inv > 0.0 && photon_density > 0.0) {
	  j = ion->ilev[r->f];
	  blk2->total_rate[j] += photon_density * r->inv;
	}
      }
    }
    for (p = 0; p < ion->ai_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->ai_rates, p);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	j = ion->ilev[r->i];
	blk1->total_rate[j] += r->dir;
	blk1->n[j] += r->dir;
	if (r->inv > 0.0 && electron_density > 0.0) {
	  j = ion->ilev[r->f];
	  blk2->total_rate[j] += electron_density * r->inv;
	}
      }
    }
    if (electron_density > 0.0) {
      for (p = 0; p < ion->ci_rates->dim; p++) {
	brts = (BLK_RATE *) ArrayGet(ion->ci_rates, p);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	for (m = 0; m < brts->rates->dim; m++) {
	  r = (RATE *) ArrayGet(brts->rates, m);
	  j = ion->ilev[r->i];
	  blk1->total_rate[j] += electron_density * r->dir;
	  if (r->inv > 0.0) {
	    j = ion->ilev[r->f];
	    blk2->total_rate[j] += electron_density * 
	      electron_density * r->inv;
	  }
	}
      }
    }
  }

  m = -2;
  for (i = 0; i < blocks->dim; i++) {
    blk1 = (LBLOCK *) ArrayGet(blocks, i);
    for (k = 0; k < blk1->nlevels; k++) {
      if (blk1->n[k]) {
	blk1->n[k] = 0.0;
      } else {
	if (blk1->iion != m) {
	  if (blk1->nlevels > 1 && i > 0) {
	    blk1->total_rate[k] = 0.0;
	  }
	  m = blk1->iion;
	} else {
	  if (blk1->nlevels > 1) {
	    blk1->total_rate[k] = 0.0;
	  }
	}
      }
    }
  }
      
  return 0;
}

int RateTable(char *fn, int nc, char *sc[], int md) { 
  RT_RECORD rt, rt1, rt2, rt3;
  RT_HEADER rt_hdr;
  F_HEADER fhdr;
  ION *ion, *ion1;
  RATE *r;
  LBLOCK *blk, *blk1, *blk2;
  BLK_RATE *brts;
  DISTRIBUTION *edist, *pdist;
  int n, k, m, i, j, p, q, s, p1;
  double den;
  MULTI ce, tr, ci, rr, ai;
  MULTI cep, trp, cip, rrp, aip;
  int ablks[3] = {10, 5, 5};
  int index[3], index1[3];
  int *ic;
  NCOMPLEX *c, *cp;
  double *d, rtmp,  e0, abt;
  double **dce[4], **dtr[4], **drr[4], **dci[4], **dai[4];
  FILE *f;

  edist = GetEleDist(&i);
  pdist = GetPhoDist(&j);
  fhdr.type = DB_RT;
  fhdr.atom = ion0.atom;
  strcpy(fhdr.symbol, ion0.symbol);
  rt_hdr.eden = electron_density;
  rt_hdr.pden = photon_density;
  rt_hdr.iedist = i;
  rt_hdr.ipdist = j;
  rt_hdr.np_edist = edist->nparams;
  rt_hdr.np_pdist = pdist->nparams;
  rt_hdr.p_edist = edist->params;
  rt_hdr.p_pdist = pdist->params;
  f = OpenFile(fn, &fhdr);
  InitFile(f, &fhdr, &rt_hdr);

  n = blocks->dim;
  ic = (int *) malloc(sizeof(int)*n);
  if (md & 4) {
    MultiInit(&ce, sizeof(double), 3, ablks);
    MultiInit(&tr, sizeof(double), 3, ablks);
    MultiInit(&ci, sizeof(double), 3, ablks);
    MultiInit(&rr, sizeof(double), 3, ablks);
    MultiInit(&ai, sizeof(double), 3, ablks);
    MultiInit(&cep, sizeof(double), 3, ablks);
    MultiInit(&trp, sizeof(double), 3, ablks);
    MultiInit(&cip, sizeof(double), 3, ablks);
    MultiInit(&rrp, sizeof(double), 3, ablks);
    MultiInit(&aip, sizeof(double), 3, ablks);
  }
  if (nc > 0) {
    c = (NCOMPLEX *) malloc(sizeof(NCOMPLEX)*MAXNCOMPLEX*nc);
    cp = c;
  } else {
    c = NULL;
    cp = NULL;
  }
  if (nc == 1) {
    if (strcmp(sc[0], "*") == 0) {
      for (i = 0; i < n; i++) {
	ic[i] = 1;
      }
    }
  } else {
    for (i = 0; i < nc; i++) {
      GetNComplex(cp, sc[i]);
      cp += MAXNCOMPLEX;
    }
    for (i = 0; i < n; i++) {
      blk = (LBLOCK *) ArrayGet(blocks, i);
      cp = c;
      ic[i] = 0;
      for (j = 0; j < nc; j++) {
	if (CompareNComplex(blk->ncomplex, cp) == 0) {
	  ic[i] = 1;
	  break;
	}
	cp += MAXNCOMPLEX;
      }
    }
  }
  if (md & 3) {
    for (i = 0; i < 4; i++) {
      dce[i] = malloc(sizeof(double *)*n);
      dtr[i] = malloc(sizeof(double *)*n);
      drr[i] = malloc(sizeof(double *)*n);
      dci[i] = malloc(sizeof(double *)*n);
      dai[i] = malloc(sizeof(double *)*n);
      for (j = 0; j < n; j++) {
	blk = ArrayGet(blocks, j);
	if (ic[j]) m = blk->nlevels;
	else m = 1;
	dce[i][j] = malloc(sizeof(double)*m);
	dtr[i][j] = malloc(sizeof(double)*m);
	drr[i][j] = malloc(sizeof(double)*m);
	dci[i][j] = malloc(sizeof(double)*m);
	dai[i][j] = malloc(sizeof(double)*m);
	for (q = 0; q < m; q++) {
	  dce[i][j][m] = 0.0;
	  dtr[i][j][m] = 0.0;
	  drr[i][j][m] = 0.0;
	  dci[i][j][m] = 0.0;
	  dai[i][j][m] = 0.0;
	}
      }
    }
  }
  k = ions->dim - 1;
  ion = (ION *) ArrayGet(ions, k);
  e0 = 0.0;
  for (i = 0; i < ion->nlevels; i++) {
    if (ion->energy[i] < e0) e0 = ion->energy[i];
  }
  abt = ion0.nt;
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    /* store the statistical weight and energy of each level in the LBLOCK array n0
       n0,r array used in the BlockRelaxation is no longer needed, overwriting. */
    for (i = 0; i < ion->nlevels; i++) {
      blk = ion->iblock[i];
      q = ion->ilev[i];
      blk->n0[q] = ion->j[i] + 1.0;
      blk->r[q] = ion->energy[i] - e0;
    }
    abt += ion->nt;
    if (md <= 0) continue;
    for (q = 0; q < ion->ce_rates->dim; q++) {
      brts = (BLK_RATE *) ArrayGet(ion->ce_rates, q);
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	blk = ion->iblock[r->i];
	blk1 = ion->iblock[r->f];
	i = blk->ib;
	j = blk1->ib;
	if (blk == blk1 && !ic[i]) continue;
	den = blk->n[ion->ilev[r->i]];
	if (den) {
	  den *= electron_density;
	  rtmp = den*r->dir;
	  index[2] = i;
	  index[1] = j;
	  if (ic[j]) {
	    index[0] = ion->ilev[r->f];
	  } else {
	    index[0] = 0;
	  }
	  if (md & 4) {
	    d = (double *) MultiSet(&ce, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  }
	  if (md & 3) {
	    dce[0][j][index[0]] += rtmp;
	  }
	  index[2] = j;
	  index[1] = i;
	  if (ic[i]) {
	    index[0] = ion->ilev[r->i];
	  } else {
	    index[0] = 0;
	  }
	  if ((md & 4) && ic[i]) {
	    d = (double *) MultiSet(&cep, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  }
	  if (md & 3) {
	    dce[1][i][index[0]] += rtmp;
	  }
	}
	if (r->inv > 0.0) {
	  den = blk1->n[ion->ilev[r->f]];
	  if (den) {
	    den *= electron_density;
	    rtmp = den * r->inv;
	    index[2] = j;
	    index[1] = i;
	    if (ic[i]) {
	      index[0] = ion->ilev[r->i];
	    } else {
	      index[0] = 0;
	    }
	    if (md & 4) {
	      d = (double *) MultiSet(&ce, index, NULL, InitDoubleData, NULL);
	      *d += rtmp;
	    }
	    if (md & 3) {
	      dce[2][i][index[0]] += rtmp;
	    }
	    index[2] = i;
	    index[1] = j;
	    if (ic[j]) {
	      index[0] = ion->ilev[r->f];
	    } else {
	      index[0] = 0;
	    }
	    if ((md & 4) && ic[j]) {
	      d = (double *) MultiSet(&cep, index, NULL, InitDoubleData, NULL);
	      *d += rtmp;
	    }
	    if (md & 3) {
	      dce[3][j][index[0]] += rtmp;
	    }
	  }
	}
      }
    }
    for (q = 0; q < ion->tr_rates->dim; q++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr_rates, q);
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	blk = ion->iblock[r->i];
	blk1 = ion->iblock[r->f];
	i = blk->ib;
	j = blk1->ib;
	if (blk == blk1 && !ic[i]) continue;
	den = blk->n[ion->ilev[r->i]];
	if (den) {
	  rtmp = den * r->dir;
	  index[2] = i;
	  index[1] = j;
	  if (ic[j]) {
	    index[0] = ion->ilev[r->f];
	  } else {
	    index[0] = 0;
	  }
	  if (md & 4) {
	    d = (double *) MultiSet(&tr, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  }
	  if (md & 3) {
	    dtr[0][j][index[0]] += rtmp;
	  }
	  index[2] = j;
	  index[1] = i;
	  if (ic[i]) {
	    index[0] = ion->ilev[r->i];
	  } else {
	    index[0] = 0;
	  }
	  if ((md & 4) && ic[i]) {
	    d = (double *) MultiSet(&trp, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  }
	  if (md & 3) {
	    dtr[1][i][index[0]] += rtmp;
	  }
	}
	if (r->inv > 0.0 && photon_density > 0.0) {
	  den = blk1->n[ion->ilev[r->f]];
	  if (den) {
	    den *= photon_density;
	    rtmp = den * r->inv;
	    index[2] = j;
	    index[1] = i;
	    if (ic[i]) {
	      index[0] = ion->ilev[r->i];
	    } else {
	      index[0] = 0;
	    }
	    if (md & 4) {
	      d = (double *) MultiSet(&tr, index, NULL, InitDoubleData, NULL);
	      *d += rtmp;
	    }
	    if (md & 3) {
	      dtr[2][i][index[0]] += rtmp;
	    }
	    index[2] = i;
	    index[1] = j;
	    if (ic[j]) {
	      index[0] = ion->ilev[r->f];
	    } else {
	      index[0] = 0;
	    }
	    if ((md & 4) && ic[j]) {
	      d = (double *) MultiSet(&trp, index, NULL, InitDoubleData, NULL);
	      *d += rtmp;
	    }
	    if (md & 3) {
	      dtr[3][j][index[0]] += rtmp;
	    }
	  }
	}
      }
    }
    for (q = 0; q < ion->tr2_rates->dim; q++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr2_rates, q);
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	blk = ion->iblock[r->i];
	blk1 = ion->iblock[r->f];
	i = blk->ib;
	j = blk1->ib;
	if (blk == blk1 && !ic[i]) continue;
	den = blk->n[ion->ilev[r->i]];
	if (den) {
	  rtmp = den * r->dir;
	  index[2] = i;
	  index[1] = j;
	  if (ic[j]) {
	    index[0] = ion->ilev[r->f];
	  } else {
	    index[0] = 0;
	  }
	  if (md & 4) {
	    d = (double *) MultiSet(&tr, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  }
	  if (md & 3) {
	    dtr[0][j][index[0]] += rtmp;
	  }
	  index[2] = j;
	  index[1] = i;
	  if (ic[i]) {
	    index[0] = ion->ilev[r->i];
	  } else {
	    index[0] = 0;
	  }
	  if ((md & 4) && ic[i]) {
	    d = (double *) MultiSet(&trp, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  }
	  if (md & 3) {
	    dtr[1][i][index[0]] += rtmp;
	  }
	}
      }
    }
    for (q = 0; q < ion->rr_rates->dim; q++) {
      brts = (BLK_RATE *) ArrayGet(ion->rr_rates, q);
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m); 
	blk = ion->iblock[r->i];
	blk1 = ion->iblock[r->f];
	i = blk->ib;
	j = blk1->ib;
	den = blk->n[ion->ilev[r->i]];
	if (den) {
	  den *= electron_density;
	  rtmp = den * r->dir;
	  index[2] = i;
	  index[1] = j;
	  if (ic[j]) {
	    index[0] = ion->ilev[r->f];
	  } else {
	    index[0] = 0;
	  }
	  if (md & 4) {
	    d = (double *) MultiSet(&rr, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  } 
	  if (md & 3) {
	    drr[0][j][index[0]] += rtmp;
	  }
	  index[2] = j;
	  index[1] = i;
	  if (ic[i]) {
	    index[0] = ion->ilev[r->i];
	  } else {
	    index[0] = 0;
	  }
	  if ((md & 4) && ic[i]) {
	    d = (double *) MultiSet(&rrp, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  }
	  if (md & 3) {
	    drr[1][i][index[0]] += rtmp;
	  }
	}
	if (r->inv > 0.0 && photon_density > 0.0) {
	  den = blk1->n[ion->ilev[r->f]];
	  if (den) {
	    den *= photon_density;
	    rtmp = den * r->inv;
	    index[2] = j;
	    index[1] = i;
	    if (ic[i]) {
	      index[0] = ion->ilev[r->i];
	    } else {
	      index[0] = 0;
	    }
	    if (md & 4) {
	      d = (double *) MultiSet(&rr, index, NULL, InitDoubleData, NULL);
	      *d += rtmp;
	    }
	    if (md & 3) {
	      drr[2][i][index[0]] += rtmp;
	    }
	    index[2] = i;
	    index[1] = j;
	    if (ic[j]) {
	      index[0] = ion->ilev[r->f];
	    } else {
	      index[0] = 0;
	    }
	    if ((md & 4) && ic[j]) {
	      d = (double *) MultiSet(&rrp, index, NULL, InitDoubleData, NULL);
	      *d += rtmp;
	    }
	    if (md & 3) {
	      dtr[3][j][index[0]] += rtmp;
	    }
	  }
	}
      }
    }
    for (q = 0; q < ion->ai_rates->dim; q++) {
      brts = (BLK_RATE *) ArrayGet(ion->ai_rates, q);
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m); 
	blk = ion->iblock[r->i];
	blk1 = ion->iblock[r->f];
	i = blk->ib;
	j = blk1->ib;
	den = blk->n[ion->ilev[r->i]];
	if (den) {
	  rtmp = den * r->dir;
	  index[2] = i;
	  index[1] = j;
	  if (ic[j]) {
	    index[0] = ion->ilev[r->f];
	  } else {
	    index[0] = 0;
	  }
	  if (md & 4) {
	    d = (double *) MultiSet(&ai, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  } 
	  if (md & 3) {
	    dai[0][j][index[0]] += rtmp;
	  }
	  index[2] = j;
	  index[1] = i;
	  if (ic[i]) {
	    index[0] = ion->ilev[r->i];
	  } else {
	    index[0] = 0;
	  }
	  if ((md & 4) && ic[i]) {
	    d = (double *) MultiSet(&aip, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  }
	  if (md & 3) {
	    dai[1][i][index[0]] += rtmp;
	  }
	}
	if (r->inv > 0.0) {
	  den = blk1->n[ion->ilev[r->f]];
	  if (den) {
	    den *= electron_density;
	    rtmp = den * r->inv;
	    index[2] = j;
	    index[1] = i;
	    if (ic[i]) {
	      index[0] = ion->ilev[r->i];
	    } else {
	      index[0] = 0;
	    }
	    if (md & 4) {
	      d = (double *) MultiSet(&ai, index, NULL, InitDoubleData, NULL);
	      *d += rtmp;
	    }
	    if (md & 3) {
	      dai[2][i][index[0]] += rtmp;
	    }
	    index[2] = i;
	    index[1] = j;
	    if (ic[j]) {
	      index[0] = ion->ilev[r->f];
	    } else {
	      index[0] = 0;
	    }
	    if ((md & 4) && ic[j]) {
	      d = (double *) MultiSet(&aip, index, NULL, InitDoubleData, NULL);
	      *d += rtmp;
	    }
	    if (md & 3) {
	      dai[3][j][index[0]] += rtmp;
	    }
	  }
	}
      }  
    }
    for (q = 0; q < ion->ci_rates->dim; q++) {
      brts = (BLK_RATE *) ArrayGet(ion->ci_rates, q);
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m); 
	blk = ion->iblock[r->i];
	blk1 = ion->iblock[r->f];
	i = blk->ib;
	j = blk1->ib;
	den = blk->n[ion->ilev[r->i]];
	if (den) {
	  den *= electron_density;
	  rtmp = den * r->dir;
	  index[2] = i;
	  index[1] = j;
	  if (ic[j]) {
	    index[0] = ion->ilev[r->f];
	  } else {
	    index[0] = 0;
	  }
	  if (md & 4) {
	    d = (double *) MultiSet(&ci, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  } 
	  if (md & 3) {
	    dci[0][j][index[0]] += rtmp;
	  }
	  index[2] = j;
	  index[1] = i;
	  if (ic[i]) {
	    index[0] = ion->ilev[r->i];
	  } else {
	    index[0] = 0;
	  }
	  if ((md & 4) && ic[i]) {
	    d = (double *) MultiSet(&cip, index, NULL, InitDoubleData, NULL);
	    *d += rtmp;
	  }
	  if (md & 3) {
	    dci[1][i][index[0]] += rtmp;
	  }
	}
	if (r->inv > 0.0) {
	  den = blk1->n[ion->ilev[r->f]];
	  if (den) {
	    den *= electron_density*electron_density;
	    rtmp = den * r->inv;
	    index[2] = j;
	    index[1] = i;
	    if (ic[i]) {
	      index[0] = ion->ilev[r->i];
	    } else {
	      index[0] = 0;
	    }
	    if (md & 4) {
	      d = (double *) MultiSet(&ci, index, NULL, InitDoubleData, NULL);
	      *d += rtmp;
	    }
	    if (md & 3) {
	      dci[2][i][index[0]] += rtmp;
	    }
	    index[2] = i;
	    index[1] = j;
	    if (ic[j]) {
	      index[0] = ion->ilev[r->f];
	    } else {
	      index[0] = 0;
	    }
	    if ((md & 4) && ic[j]) {
	      d = (double *) MultiSet(&cip, index, NULL, InitDoubleData, NULL);
	      *d += rtmp;
	    }
	    if (md & 3) {
	      dci[3][j][index[0]] += rtmp;
	    }
	  }
	}  
      }
    }
  }

  for (i = 0; i < n; i++) {
    blk = (LBLOCK *) ArrayGet(blocks, i);
    k = blk->iion;
    if (k < 0) {
      m = ion0.nele;
      ion = (ION *) ArrayGet(ions, 0);
      rt1.nb = 0.0;
      rt2.nb = ion0.nt;
      rt3.nb = ion->nt;
    } else {
      ion = (ION *) ArrayGet(ions, k);
      m = ion->nele;
      if (k == 0) {
	rt1.nb = ion0.nt;
      } else {
	rt1.nb = ((ION *) ArrayGet(ions, k-1))->nt;
      }
      rt2.nb = ion->nt;
      if (k+1 < ions->dim) {
	rt3.nb = ((ION *) ArrayGet(ions, k+1))->nt;
      } else {
	rt3.nb = 0.0;
      }
    }
    sprintf(rt1.icomplex, "%3d", m-1);
    sprintf(rt2.icomplex, "%3d", m);
    sprintf(rt3.icomplex, "%3d", m+1);
    rt1.iblock = -1;
    rt2.iblock = -2;
    rt3.iblock = -3;
    if (ic[i]) {
      for (k = 0; k < blk->nlevels; k++) {
	StrNComplex(rt.icomplex, blk->ncomplex);
	rt.nb = blk->n[k];
	rt.iblock = i;
	if (!(blk->n[k])) continue;
	rt.dir = IonIndex(ion, i, k);
	if (blk->iion < 0 && ion0.nionized > 0) {
	  rt.dir = IonizedIndex(rt.dir, 1);
	  rt.dir = ion0.ionized_map[0][rt.dir];
	}
	rt.tr = blk->n0[k];
	rt.ce = blk->r[k];
	rt.rr = m;
	rt.ai = abt;
	rt.ci = -1.0;
	WriteRTRecord(f, &rt);
	if (md <= 0) {
	  continue;
	}	
	if ((md & 4) && (md & 1) ) {
	  index[0] = k;
	  for (j = 0; j < n; j++) {
	    rt.iblock = j;
	    blk1 = (LBLOCK *) ArrayGet(blocks, j);
	    if (abs(blk1->iion - blk->iion) > 1) continue;
	    rt.nb = blk1->nb;
	    index[2] = j;
	    index[1] = i;
	    rt.ce = 0.0;
	    rt.tr = 0.0;
	    rt.rr = 0.0;
	    rt.ai = 0.0;
	    rt.ci = 0.0;
	    rt.dir = 0;
	    if (blk1->iion == blk->iion) {
	      d = (double *) MultiGet(&ce, index);
	      if (d && *d) {
		rt.ce = *d;
	      }
	      d = (double *) MultiGet(&tr, index);
	      if (d && *d) {
		rt.tr = *d;
	      }
	    } else {
	      d = (double *) MultiGet(&rr, index);
	      if (d && *d) {
		rt.rr = *d;
	      }
	      d = (double *) MultiGet(&ai, index);
	      if (d && *d) {
		rt.ai = *d;
	      }
	      d = (double *) MultiGet(&ci, index);
	      if (d && *d) {
		rt.ci = *d;
	      }
	    }
	    if (rt.ce || rt.tr ||rt.rr || rt.ai || rt.ci) {
	      StrNComplex(rt.icomplex, blk1->ncomplex);
	      WriteRTRecord(f, &rt);
	    }
	  }
	}
	if (md & 1) {
	  rt1.dir = 0;
	  rt2.dir = 0;
	  rt3.dir = 0;
	  rt1.tr = 0.0;
	  rt1.ce = 0.0;
	  rt1.rr = drr[0][i][k];
	  rt1.ai = dai[2][i][k];
	  rt1.ci = dci[2][i][k];
	  rt2.tr = dtr[0][i][k] + dtr[2][i][k];
	  rt2.ce = dce[0][i][k] + dce[2][i][k];
	  rt2.rr = 0.0;
	  rt2.ai = 0.0;
	  rt2.ci = 0.0;
	  rt3.tr = 0.0;
	  rt3.ce = 0.0;
	  rt3.rr = drr[2][i][k];
	  rt3.ai = dai[0][i][k];
	  rt3.ci = dai[0][i][k];
	  WriteRTRecord(f, &rt1);
	  WriteRTRecord(f, &rt2);
	  WriteRTRecord(f, &rt3);
	}
	if ((md & 4) && (md & 2)) {
	  index[0] = k;
	  for (j = 0; j < n; j++) {
	    rt.iblock = j;
	    blk1 = (LBLOCK *) ArrayGet(blocks, j);
	    if (abs(blk1->iion - blk->iion) > 1) continue;
	    rt.nb = blk1->nb;
	    index[2] = j;
	    index[1] = i;
	    rt.ce = 0.0;
	    rt.tr = 0.0;
	    rt.rr = 0.0;
	    rt.ai = 0.0;
	    rt.ci = 0.0;
	    rt.dir = 1;
	    if (blk1->iion == blk->iion) {
	      d = (double *) MultiGet(&cep, index);
	      if (d && *d) {
		rt.ce = *d;
	      }
	      d = (double *) MultiGet(&trp, index);
	      if (d && *d) {
		rt.tr = *d;
	      }
	    } else {
	      d = (double *) MultiGet(&rrp, index);
	      if (d && *d) {
		rt.rr = *d;
	      }
	      d = (double *) MultiGet(&aip, index);
	      if (d && *d) {
		rt.ai = *d;
	      }
	      d = (double *) MultiGet(&cip, index);
	      if (d && *d) {
		rt.ci = *d;
	      }
	    }
	    if (rt.ce || rt.tr ||rt.rr || rt.ai || rt.ci) {
	      StrNComplex(rt.icomplex, blk1->ncomplex);
	      WriteRTRecord(f, &rt);
	    }
	  }
	}
	if (md & 2) {
	  rt1.dir = 1;
	  rt2.dir = 1;
	  rt3.dir = 1;
	  rt1.tr = 0.0;
	  rt1.ce = 0.0;
	  rt1.rr = drr[3][i][k];
	  rt1.ai = dai[1][i][k];
	  rt1.ci = dci[1][i][k];
	  rt2.tr = dtr[1][i][k] + dtr[3][i][k];
	  rt2.ce = dce[1][i][k] + dce[3][i][k];
	  rt2.rr = 0.0;
	  rt2.ai = 0.0;
	  rt2.ci = 0.0;
	  rt3.tr = 0.0;
	  rt3.ce = 0.0;
	  rt3.rr = drr[1][i][k];
	  rt3.ai = dai[3][i][k];
	  rt3.ci = dci[3][i][k];
	  WriteRTRecord(f, &rt1);
	  WriteRTRecord(f, &rt2);
	  WriteRTRecord(f, &rt3);
	}
      }
    } else {
      index[0] = 0;
      rt.nb = blk->nb;
      StrNComplex(rt.icomplex, blk->ncomplex);
      if (!(blk->nb)) continue;
      rt.dir = IonIndex(ion, i, 0);
      if (blk->iion < 0 && ion0.nionized > 0) {
	rt.dir = IonizedIndex(rt.dir, 1);
	rt.dir = ion0.ionized_map[0][rt.dir];
      }
      rt.tr = 0.0;
      rt.ce = 0.0;
      for (k = 0; k < blk->nlevels; k++) {
	rt.tr += blk->n0[k];
	rt.ce += blk->n[k]*blk->r[k];
      }
      rt.ce /= blk->nb;
      rt.dir = -(rt.dir+1);
      rt.rr = m;
      rt.ai = abt;
      rt.ci = -1.0;
      WriteRTRecord(f, &rt);
      if (md <= 0) {
	continue;
      }
      if ((md & 4) && (md & 1)) {
	for (j = 0; j < n; j++) {
	  rt.iblock = j;
	  blk1 = (LBLOCK *) ArrayGet(blocks, j);
	  if (abs(blk1->iion - blk->iion) > 1) continue;
	  rt.nb = blk1->nb;
	  index[2] = j;
	  index[1] = i;
	  rt.ce = 0.0;
	  rt.tr = 0.0;
	  rt.rr = 0.0;
	  rt.ai = 0.0;
	  rt.ci = 0.0;
	  rt.dir = 0;
	  if (blk1->iion == blk->iion) {
	    d = (double *) MultiGet(&ce, index);
	    if (d && *d) {
	      rt.ce = *d;
	    } 
	    d = (double *) MultiGet(&tr, index);
	    if (d && *d) {
	      rt.tr = *d;
	    }
	  } else {
	    d = (double *) MultiGet(&rr, index);
	    if (d && *d) {
	      rt.rr = *d;
	    }
	    d = (double *) MultiGet(&ai, index);
	    if (d && *d) {
	      rt.ai = *d;
	    }
	    d = (double *) MultiGet(&ci, index);
	    if (d && *d) {
	      rt.ci = *d;
	    }
	  }
	  if (rt.ce || rt.tr ||rt.rr || rt.ai || rt.ci) {
	    StrNComplex(rt.icomplex, blk1->ncomplex);
	    WriteRTRecord(f, &rt);	
	  }
	}
      }
      if (md & 1) {
	rt1.dir = 0;
	rt2.dir = 0;
	rt3.dir = 0;
	rt1.tr = 0.0;
	rt1.ce = 0.0;
	rt1.rr = drr[0][i][0];
	rt1.ai = dai[2][i][0];
	rt1.ci = dci[2][i][0];
	rt2.tr = dtr[0][i][0] + dtr[2][i][0];
	rt2.ce = dce[0][i][0] + dce[2][i][0];
	rt2.rr = 0.0;
	rt2.ai = 0.0;
	rt2.ci = 0.0;
	rt3.tr = 0.0;
	rt3.ce = 0.0;
	rt3.rr = drr[2][i][0];
	rt3.ai = dai[0][i][0];
	rt3.ci = dci[0][i][0];      
	WriteRTRecord(f, &rt1);
	WriteRTRecord(f, &rt2);
	WriteRTRecord(f, &rt3);
      }
      if ((md & 4) && (md & 2)) {
	for (j = 0; j < n; j++) {
	  rt.iblock = j;
	  if (abs(blk1->iion - blk->iion) > 1) continue;
	  blk1 = (LBLOCK *) ArrayGet(blocks, j);
	  rt.nb = blk1->nb;
	  index[2] = i;
	  index[1] = j;
	  rt.ce = 0.0;
	  rt.tr = 0.0;
	  rt.rr = 0.0;
	  rt.ai = 0.0;
	  rt.ci = 0.0;
	  rt.dir = 1;
	  if (blk1->iion == blk->iion) {
	    d = (double *) MultiGet(&ce, index);
	    if (d && *d) {
	      rt.ce = *d;
	    } 
	    d = (double *) MultiGet(&tr, index);
	    if (d && *d) {
	      rt.tr = *d;
	    }
	  } else {
	    d = (double *) MultiGet(&rr, index);
	    if (d && *d) {
	      rt.rr = *d;
	    }
	    d = (double *) MultiGet(&ai, index);
	    if (d && *d) {
	      rt.ai = *d;
	    }
	    d = (double *) MultiGet(&ci, index);
	    if (d && *d) {
	      rt.ci = *d;
	    }
	  }
	  if (rt.ce || rt.tr ||rt.rr || rt.ai || rt.ci) {
	    StrNComplex(rt.icomplex, blk1->ncomplex);
	    WriteRTRecord(f, &rt);
	  }
	}
      }
      if (md & 2) {
	rt1.dir = 1;
	rt2.dir = 1;
	rt3.dir = 1;
	rt1.tr = 0.0;
	rt1.ce = 0.0;
	rt1.rr = drr[3][i][0];
	rt1.ai = dai[1][i][0];
	rt1.ci = dci[1][i][0];
	rt2.tr = dtr[1][i][0] + dtr[3][i][0];
	rt2.ce = dce[1][i][0] + dce[3][i][0];
	rt2.rr = 0.0;
	rt2.ai = 0.0;
	rt2.ci = 0.0;
	rt3.tr = 0.0;
	rt3.ce = 0.0;
	rt3.rr = drr[1][i][0];
	rt3.ai = dai[3][i][0];
	rt3.ci = dci[3][i][0];
	WriteRTRecord(f, &rt1);
	WriteRTRecord(f, &rt2);
	WriteRTRecord(f, &rt3);
      }
    }
  }
  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);

  free(ic);
  if (c) free(c);
  if (md & 4) {
    MultiFree(&ce, NULL);
    MultiFree(&tr, NULL);
    MultiFree(&ci, NULL);
    MultiFree(&rr, NULL);
    MultiFree(&ai, NULL);
    MultiFree(&cep, NULL);
    MultiFree(&trp, NULL);
    MultiFree(&cip, NULL);
    MultiFree(&rrp, NULL);
    MultiFree(&aip, NULL);
  }
  if (md & 3) {
    for (i = 0; i < 4; i++) {
      for (j = 0; j < n; j++) {
	free(dtr[i][j]);
	free(dce[i][j]);
	free(drr[i][j]);
	free(dai[i][j]);
	free(dci[i][j]);
      }
      free(dtr[i]);
      free(dce[i]);
      free(drr[i]);
      free(dai[i]);
      free(dci[i]);
    }
  }

  return 0;
}

int BlockMatrix(void) {
  ION *ion;
  RATE *r;
  LBLOCK *blk1, *blk2;
  BLK_RATE *brts;
  int n, k, m, i, j, t;
  int p, q, iion, k0, k1;
  double *x, den, a;
  
  n = blocks->dim;
  for (i = 0; i < 2*n*(n+1); i++) {
    bmatrix[i] = 0.0;
  }
  x = bmatrix + n*n;

  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    if (electron_density > 0.0) {
      for (t = 0; t < ion->ce_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->ce_rates, t);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	if (blk1 == blk2) continue;
	if (rec_cascade && (blk1->rec || blk2->rec)) continue;
	for (m = 0; m < brts->rates->dim; m++) {
	  r = (RATE *) ArrayGet(brts->rates, m);
	  i = ion->iblock[r->i]->ib;
	  j = ion->iblock[r->f]->ib;
	  den = blk1->r[ion->ilev[r->i]];
	  if (den) {
	    p = i*n + j;
	    bmatrix[p] += den * electron_density * r->dir;
	  }
	  if (r->inv > 0.0) {
	    den = blk2->r[ion->ilev[r->f]];
	    if (den) {
	      p = i + j*n;
	      bmatrix[p] += den * electron_density * r->inv;
	    }
	  }
	}
      }
    }
    for (t = 0; t < ion->tr_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (blk1 == blk2) continue;
      if (rec_cascade && (blk1->rec || blk2->rec)) continue;
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	i = ion->iblock[r->i]->ib;
	j = ion->iblock[r->f]->ib;
	den = blk1->r[ion->ilev[r->i]];
	if (den) {
	  p = i*n + j;
	  bmatrix[p] += den * r->dir;
	}
	if (r->inv > 0.0 && photon_density > 0.0) {
	  if (den) {
	    a = photon_density * r->inv;
	    bmatrix[p] += den*a*(ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
	  }	  
	  den = blk2->r[ion->ilev[r->f]];
	  if (den) {
	    p = i + j*n;
	    bmatrix[p] += den * photon_density * r->inv;
	  }
	}
      }
    }
    for (t = 0; t < ion->tr2_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr2_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (blk1 == blk2) continue;
      if (rec_cascade && (blk1->rec || blk2->rec)) continue;
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	i = ion->iblock[r->i]->ib;
	j = ion->iblock[r->f]->ib;
	den = blk1->r[ion->ilev[r->i]];
	if (den) {
	  p = i*n + j;
	  bmatrix[p] += den * r->dir;
	}
      }
    }
    for (t = 0; t < ion->rr_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->rr_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (blk1 == blk2) continue;
      if (rec_cascade && (blk1->rec || blk2->rec)) continue;
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m); 
	i = ion->iblock[r->i]->ib;
	j = ion->iblock[r->f]->ib;
	den = blk1->r[ion->ilev[r->i]];
	if (den) {
	  if (electron_density > 0.0) {
	    p = i*n + j;
	    bmatrix[p] += den * electron_density * r->dir;
	  }
	}
	if (r->inv > 0.0 && photon_density > 0.0) {
	  den = blk2->r[ion->ilev[r->f]];
	  if (den) {
	    p = i + j*n;
	    bmatrix[p] += den * photon_density * r->inv;
	  }
	}
      }
    }
    for (t = 0; t < ion->ai_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->ai_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (blk1 == blk2) continue;
      if (rec_cascade && (blk1->rec || blk2->rec)) continue;
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m); 
	i = ion->iblock[r->i]->ib;
	j = ion->iblock[r->f]->ib;
	den = blk1->r[ion->ilev[r->i]];
	if (den) {
	  p = i*n + j;
	  bmatrix[p] += den * r->dir;
	}
	if (r->inv > 0.0 && electron_density > 0.0) {
	  den = blk2->r[ion->ilev[r->f]];
	  if (den) {
	    p = i + j*n;
	    bmatrix[p] += den * electron_density * r->inv;
	  }
	}
      }  
    }
    if (electron_density > 0.0) {
      for (t = 0; t < ion->ci_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->ci_rates, t);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	if (blk1 == blk2) continue;
	if (rec_cascade && (blk1->rec || blk2->rec)) continue;
	for (m = 0; m < brts->rates->dim; m++) {
	  r = (RATE *) ArrayGet(brts->rates, m); 
	  i = ion->iblock[r->i]->ib;
	  j = ion->iblock[r->f]->ib;
	  den = blk1->r[ion->ilev[r->i]];
	  if (den) {
	    p = i*n + j;
	    bmatrix[p] += den * electron_density * r->dir;
	  }
	  if (r->inv > 0.0) {
	    den = blk2->r[ion->ilev[r->f]];
	    if (den) {
	      p = i + j*n;
	      den *= electron_density;
	      bmatrix[p] += den * electron_density * r->inv;
	    }
	  }
	}  
      }
    }
  }

  for (i = 0; i < n; i++) {
    p = i*n;
    q = i + p;
    for (j = 0; j < n; j++) {
      if (j != i) bmatrix[q] += bmatrix[p];
      p++;
    }
    bmatrix[q] = - bmatrix[q];
  }

  iion = -2;
  k0 = 0;
  k1 = 0;
  for (i = 0; i < n; i++) {
    blk1 = (LBLOCK *) ArrayGet(blocks, i);
    if (blk1->iion != iion) {
      if (iion != -2) {
	k = iion;
	if (k == -1) k = 0;
	ion = (ION *) ArrayGet(ions, k);
	if (iion == -1) den = ion0.n0;
	else den = ion->n0;
	if (den > 0.0) {
	  x[k0] = den;
	  p = k0;
	  for (k = 0; k < n; k++) {
	    /*
	    if (k < k1 && k >= k0) bmatrix[p] = 1.0;
	    else bmatrix[p] = 0.0;
	    */
	    if (k == k0) bmatrix[p] = 1.0;
	    else bmatrix[p] = 0.0;
	    p += n;
	  }
	} 
      }
      iion = blk1->iion;
      k0 = k1;
    }
    k1++;
  }

  k = iion;
  ion = (ION *) ArrayGet(ions, k);
  den = ion->n0;
  if (den > 0.0) {
    x[k0] = den;
    p = k0;
    for (k = 0; k < n; k++) {
      /*
      if (k < k1 && k >= k0) bmatrix[p] = 1.0;
      else bmatrix[p] = 0.0;
      */
      if (k == k0) bmatrix[p] = 1.0;
      else bmatrix[p] = 0.0;
      p += n;
    }
  } 

  return 0;
}

int BlockPopulation(void) {
  LBLOCK *blk;
  double *b, *x;
  double *a;
  int *ipiv;
  int info;
  int n, m;
  int nrhs;
  int lda, ldb;
  int i, j, p, q;

  n = blocks->dim;
  a = bmatrix + n*n;
  x = a;
  ipiv = (int *) bmatrix;
  a = a + n;
  b = a + n*n;

  p = 0;
  q = 0;
  m = 0;
  for (i = 0; i < n; i++) {
    if (bmatrix[i+i*n] == 0) {
      x[i] = -1.0;
      p += n;
      continue;
    }
    for (j = 0; j < n; j++) {
      a[q] = bmatrix[p];
      p++;
      q++;
    }
    b[i] = x[i];
    m++;
  }

  p = 0;
  for (i = 0; i < m; i++) {
    q = p;
    for (j = 0; j < n; j++) {
      if (x[j] < 0.0) {
	continue;
      }
      a[q] = a[p+j];
      q++;
    }
    p += n;
  }

  q = 0;
  for (j = 0; j < n; j++) {
    if (x[j] < 0.0) {
      continue;
    }
    b[q] = b[j];
    q++;
  }  

  nrhs = 1;
  lda = n;
  ldb = n;

  DGESV(m, nrhs, a, lda, ipiv, b, ldb, &info);
  if (info != 0) {
    printf("Error in solving BlockMatrix\n");
    exit(1);
  }
  
  p = 0;
  for (i = 0; i < n; i++) {
    blk = (LBLOCK *) ArrayGet(blocks, i);
    if (rec_cascade && blk->rec) continue;
    if (x[i] < 0.0) {
      blk->nb = 0.0;
      for (j = 0; j < blk->nlevels; j++) {
	blk->n[j] = 0.0;
      }
    } else {
      blk->nb = b[p++];
    }
  }

  return 0;
}
  
double BlockRelaxation(int iter) {
  ION *ion;
  RATE *r;
  LBLOCK *blk1, *blk2;
  BLK_RATE *brts;
  int i, j, k, m, t;
  int p, q;
  double a, b, c, d, h, td;
  int nlevels;

  for (k = 0; k < blocks->dim; k++) {
    blk1 = (LBLOCK *) ArrayGet(blocks, k);
    if (blk1->nlevels == 1) {
      blk1->r[0] = blk1->nb;
      blk1->n[0] = 0.0;
    } else {
      for (m = 0; m < blk1->nlevels; m++) {
	blk1->r[m] = blk1->n[m];
	blk1->n[m] = 0.0;
      }
    }
  }
  
  b = 1.0-iter_stabilizer;
  c = iter_stabilizer;
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    if (electron_density > 0.0) {
      for (t = 0; t < ion->ce_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->ce_rates, t);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	if (rec_cascade && iter >= 0) {
	  if (blk1->rec || blk2->rec) continue;
	} 
	for (m = 0; m < brts->rates->dim; m++) {
	  r = (RATE *) ArrayGet(brts->rates, m);
	  i = ion->iblock[r->i]->ib;
	  p = ion->ilev[r->i];
	  j = ion->iblock[r->f]->ib;
	  q = ion->ilev[r->f];
	  if (blk1->r[p]) {
	    blk2->n[q] += blk1->r[p] * electron_density * r->dir;
	  }
	  if (r->inv > 0.0) {
	    if (blk2->r[q]) {
	      blk1->n[p] += blk2->r[q] * electron_density * r->inv;
	    }
	  }
	}
      }
    }
    for (t = 0; t < ion->tr_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (rec_cascade && iter >= 0) {
	if (blk1->rec || blk2->rec) continue;
      }
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	i = ion->iblock[r->i]->ib;
	p = ion->ilev[r->i];
	j = ion->iblock[r->f]->ib;
	q = ion->ilev[r->f];    
	if (blk1->r[p]) {
	  blk2->n[q] += blk1->r[p] * r->dir;
	}
	if (r->inv > 0.0 && photon_density > 0.0) {
	  a = photon_density * r->inv;	  
	  if (blk1->r[p]) {
	    blk2->n[q] += blk1->r[p]*a*(ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
	  }	  
	  if (blk2->r[q]) {
	    blk1->n[p] += blk2->r[q] * a;
	  }
	}
      }
    } 
    for (t = 0; t < ion->tr2_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr2_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (rec_cascade && iter >= 0) {
	if (blk1->rec || blk2->rec) continue;
      }
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	i = ion->iblock[r->i]->ib;
	p = ion->ilev[r->i];
	j = ion->iblock[r->f]->ib;
	q = ion->ilev[r->f];    
	if (blk1->r[p]) {
	  blk2->n[q] += blk1->r[p] * r->dir;
	}
      }
    }
    for (t = 0; t < ion->rr_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->rr_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (rec_cascade && iter >= 0) {
	if (blk1->rec || blk2->rec) continue;
      }
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	i = ion->iblock[r->i]->ib;
	p = ion->ilev[r->i];
	j = ion->iblock[r->f]->ib;
	q = ion->ilev[r->f];    
	if (electron_density > 0.0) {
	  if (blk1->r[p]) {
	    blk2->n[q] += blk1->r[p] * electron_density * r->dir;
	  }
	} 
	if (r->inv > 0.0 && photon_density > 0.0) {
	  if (blk2->r[q]) {
	    blk1->n[p] += blk2->r[q] * photon_density * r->inv;
	  }
	}
      }
    }
    for (t = 0; t < ion->ai_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->ai_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (rec_cascade && iter >= 0) {
	if (blk1->rec || blk2->rec) continue;
      }
      for (m = 0; m < brts->rates->dim; m++) {
	r = (RATE *) ArrayGet(brts->rates, m);
	i = ion->iblock[r->i]->ib;
	p = ion->ilev[r->i];
	j = ion->iblock[r->f]->ib;
	q = ion->ilev[r->f];    
	if (blk1->r[p]) {
	  blk2->n[q] += blk1->r[p] * r->dir;
	}
	if (r->inv > 0.0 && electron_density > 0.0) {
	  if (blk2->r[q]) {
	    blk1->n[p] += blk2->r[q] * electron_density * r->inv;
	  }
	}
      }
    }
    if (electron_density > 0.0) {
      for (t = 0; t < ion->ci_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->ci_rates, t);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	if (rec_cascade && iter >= 0) {
	  if (blk1->rec || blk2->rec) continue;
	}
	for (m = 0; m < brts->rates->dim; m++) {
	  r = (RATE *) ArrayGet(brts->rates, m);
	  i = ion->iblock[r->i]->ib;
	  p = ion->ilev[r->i];
	  j = ion->iblock[r->f]->ib;
	  q = ion->ilev[r->f];    
	  if (blk1->r[p]) {
	    blk2->n[q] += blk1->r[p] * electron_density * r->dir;
	  }
	  if (r->inv) {
	    if (blk2->r[q]) {
	      blk1->n[p] += blk2->r[q] * electron_density * 
		electron_density * r->inv;
	    }
	  }
	}
      }
    }
  }

  nlevels = 0;
  d = 0.0;
  td = 0.0;
  for (k = 0; k < blocks->dim; k++) {
    blk1 = (LBLOCK *) ArrayGet(blocks, k);
    if (rec_cascade && iter >= 0) {
      if (blk1->rec) continue;
    }
    
    if (iter < 0 || blk1->nlevels > 1) {
      a = 0.0;
      for (m = 0; m < blk1->nlevels; m++) {
	if (blk1->total_rate[m]) {
	  blk1->n[m] /= blk1->total_rate[m];
	  a += blk1->n[m];
	} else {
	  blk1->n[m] = 0.0;
	}
      }
      if (a) {
	if (iter >= 0) {
	  a = blk1->nb/a;
	  for (m = 0; m < blk1->nlevels; m++) {
	    blk1->n[m] *= a;
	  }
	} else {
	  blk1->nb = a;
	}
      }  
      if (nlevels == 1) {
	blk1->nb = blk1->n[0];
	a = blk1->nb;
      }
    } else {
      blk1->n[0] = blk1->nb;
      a = blk1->nb;
    }

    if (blk1->nb == 0.0) {
      for (m = 0; m < blk1->nlevels; m++) {
	blk1->n[m] = 0.0;
	blk1->r[m] = 0.0;
	blk1->n0[m] = 0.0;
      }
      blk1->r[0] = 1.0;
      continue;
    }
    
    if (iter >= 0) {
      if (iter > 0) {
	if (blk1->iion < 0) a = ion0.nt;
	else {
	  ion = (ION *) ArrayGet(ions, blk1->iion);
	  a = ion->nt;
	}
      } else a = 1.0;
      for (m = 0; m < blk1->nlevels; m++) {
	if (blk1->n[m]) {
	  /*d += fabs(1.0 - blk1->n0[m]/blk1->n[m]);*/
	  d += fabs((blk1->n[m]-blk1->n0[m])/blk1->n[m])*blk1->nb;
	  td += blk1->nb;
	  nlevels += 1;
	}    
	if (iter >= 2) {
	  blk1->n[m] = b*blk1->n0[m] + c*blk1->n[m];
	}
	blk1->r[m] = blk1->n[m]/blk1->nb;
	blk1->n0[m] = blk1->n[m];  
      }
    }
  }    

  q = 0;
  p = -1;
  a = 0.0;
  for (k = 0; k < blocks->dim; k++) {
    blk1 = (LBLOCK *) ArrayGet(blocks, k);
    if (blk1->iion != p) {
      if (p == -1) {
	h = ion0.n;
      } else {
	ion = (ION *) ArrayGet(ions, p);
	h = ion->n;
      }
      /*
      if (h > 0.0) {
	if (p == -1) ion0.nt = h;
	else ion->nt = h;
	if (iter < 0) {
	  a = h/a;
	  for (i = q; i < k; i++) {
	    blk2 = (LBLOCK *) ArrayGet(blocks, i);
	    blk2->nb *= a;
	    for (m = 0; m < blk2->nlevels; m++) {
	      blk2->n[m] *= a;
	    }
	  }
	}
      } else {
	if (p == -1) ion0.nt = a;
	else ion->nt = a;
      }
      */
      if (p == -1) ion0.nt = a;
      else ion->nt = a;
      /*
      if (h > 0.0) {
	blk2 = (LBLOCK *) ArrayGet(blocks, q);
	if (p == -1) ion0.n0 = blk2->nb*h/a;
	else ion->n0 = blk2->nb*h/a;
      }
      */
      p = blk1->iion;
      a = 0.0;
      q = k;
    }
    a += blk1->nb;
  }
  ion = (ION *) ArrayGet(ions, p);
  /*
  if (ion->n > 0.0) {
    ion->nt = ion->n;
    if (iter < 0) {
      a = ion->n/a;
      for (i = q; i < k; i++) {
	blk2 = (LBLOCK *) ArrayGet(blocks, i);
	blk2->nb *= a;
	for (m = 0; m < blk2->nlevels; m++) {
	  blk2->n[m] *= a;
	}
      }
    }
  } else {
    ion->nt = a;
  }
  */
  ion->nt = a;
  /*
  if (ion->n > 0.0) {
    blk2 = ArrayGet(blocks, q);
    ion->n0 = blk2->nb*ion->n/a;
  }
  */
  if (iter < 0) {
    for (k = 0; k < blocks->dim; k++) {
      blk1 = (LBLOCK *) ArrayGet(blocks, k);
      if (blk1->iion < 0) a = ion0.nt;
      else {
	ion = (ION *) ArrayGet(ions, blk1->iion);
	a = ion->nt;
      }
      for (m = 0; m < blk1->nlevels; m++) {
	if (blk1->n[m] && blk1->rec == NULL) {
	  /*d += fabs(1.0 - blk1->n0[m]/blk1->n[m]);*/
	  d += fabs((blk1->n[m] - blk1->n0[m])/blk1->n[m])*blk1->nb;
	  td += blk1->nb;
	  nlevels += 1;
	}
	if (iter <= -2) {
	  blk1->n[m] = b*blk1->n0[m] + c*blk1->n[m];
	}
	blk1->r[m] = blk1->n[m]/blk1->nb;
	blk1->n0[m] = blk1->n[m];
      }
    }
  }
  if (iter == 0) return 1.0;
  d /= td;
  return d;
}

int LevelPopulation(void) {
  int i;
  double d;

  printf("Populate Iteration:\n");
  for (i = 0; i < max_iter; i++) {
    BlockMatrix();
    BlockPopulation();
    d = BlockRelaxation(i);
    printf("%5d %11.4E\n", i, d);
    if (d < iter_accuracy) break;
  }
 
  if (i == max_iter) {
    printf("Max iteration reached\n");
  }
  return 0;
}

int Cascade(void) {
  int i;
  double d;
  
  if (!rec_cascade) return 0;
  printf("Cascade  Iteration:\n");
  d = BlockRelaxation(-1);
  for (i = 1; i <= max_iter; i++) {
    d = BlockRelaxation(-i);
    printf("%5d %11.4E\n", i, d);
    fflush(stdout);
    if (d < cas_accuracy) break;
  }
  
  if (i == max_iter) {
    printf("Max iteration reached in Cascade\n");
  }

  return 0;
}

int SpecTable(char *fn, int rrc, double strength_threshold) {
  SP_RECORD r;
  SP_EXTRA rx;
  SP_HEADER sp_hdr;
  F_HEADER fhdr;
  ION *ion;
  RATE *rt, *dev;
  LBLOCK *blk, *iblk, *fblk;
  BLK_RATE *brts, *brdev;
  int k, m, j;
  FILE *f;
  double e, a, e0;
  int i, p, q, ib, iuta;
  double smax, s;

  iuta = IsUTA();
  fhdr.type = DB_SP;
  fhdr.atom = ion0.atom;
  strcpy(fhdr.symbol, ion0.symbol);
  f = OpenFile(fn, &fhdr);

  k = ions->dim - 1;
  ion = (ION *) ArrayGet(ions, k);
  e0 = 0.0;
  for (m = 0; m < ion->nlevels; m++) {
    if (ion->energy[m] < e0) e0 = ion->energy[m];
  }
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    sp_hdr.type = 0;
    ib = -1;
    for (m = 0; m < ion->nlevels; m++) {
      blk = ion->iblock[m];
      if (blk == NULL) continue;
      if (blk->iion != k && k > 0) continue;
      i = blk->ib;
      if (i != ib) {
	if (ib >= 0) DeinitFile(f, &fhdr);
	if (blk->iion != k) {
	  sp_hdr.nele = ion->nele - 1;
	} else {
	  sp_hdr.nele = ion->nele;
	}
	sp_hdr.iblock = i;
	StrNComplex(sp_hdr.icomplex, blk->ncomplex);
	sp_hdr.fblock = 0;
	sp_hdr.fcomplex[0] = '\0';
	InitFile(f, &fhdr, &sp_hdr);
	ib = i;
      }
      p = ion->ilev[m];
      /*      if (blk->n[p]) {*/
      r.upper = m;
      r.lower = p;
      r.energy = ion->energy[m]-e0;
      r.rrate = ion->j[m]+1.0;
      r.trate = blk->total_rate[p];
      rx.sdev = 0.0;
      r.strength = blk->n[p];
      WriteSPRecord(f, &r, &rx);
	/*}*/
    }
    if (ib >= 0) DeinitFile(f, &fhdr);
    if (rrc < 0) continue;

    for (i = 0; i < ion->tr_rates->dim; i++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr_rates, i);
      if (brts->rates->dim == 0) continue;
      if (iuta) {
	brdev = (BLK_RATE *) ArrayGet(ion->tr_sdev, i);
      }
      iblk = brts->iblock;
      fblk = brts->fblock;
      if (iblk->iion != k) {
	sp_hdr.nele = ion->nele - 1;
      } else {
	sp_hdr.nele = ion->nele;
      }
      sp_hdr.iblock = iblk->ib;
      sp_hdr.fblock = fblk->ib;
      sp_hdr.type = TransitionType(iblk->ncomplex, fblk->ncomplex);
      StrNComplex(sp_hdr.icomplex, iblk->ncomplex);
      StrNComplex(sp_hdr.fcomplex, fblk->ncomplex);
      InitFile(f, &fhdr, &sp_hdr);
      smax = 0.0;
      for (m = 0; m < brts->rates->dim; m++) {
	rt = (RATE *) ArrayGet(brts->rates, m);
	if (k == 0 && 
	    ion0.nionized > 0 &&
	    (p = IonizedIndex(rt->i, 1)) >= 0 &&
	    (q = IonizedIndex(rt->f, 1)) >= 0) {
	  e = ion0.energy[p] - ion0.energy[q];
	  p = ion0.ionized_map[0][p];
	  q = ion0.ionized_map[0][q];
	} else {
	  p = rt->i;
	  q = rt->f;
	  e = ion->energy[p] - ion->energy[q];
	}
	j = ion->ilev[rt->i];
	if (iblk->n[j] > 0.0) {
	  r.lower = q;
	  r.upper = p;
	  r.energy = e;
	  r.strength = iblk->n[j] * rt->dir;
	  if (rt->inv > 0.0 && photon_density > 0.0) {
	    a = photon_density * rt->inv;
	    a *= (ion->j[rt->f]+1.0)/(ion->j[rt->i]+1.0);
	    r.strength += iblk->n[j] * a;
	  }
	  s = r.strength*e;
	  if (s < strength_threshold*smax) continue;
	  if (s > smax) smax = s;
	  if (iuta) {
	    dev = (RATE *) ArrayGet(brdev->rates, m);
	    r.energy = dev->dir;
	    rx.sdev = dev->inv;
	  }
	  r.rrate = rt->dir;
	  r.trate = iblk->total_rate[j];
	  WriteSPRecord(f, &r, &rx);
	}
      }
      DeinitFile(f, &fhdr);
    }

    if (!rrc) continue;
    for (i = 0; i < ion->rr_rates->dim; i++) {
      brts = (BLK_RATE *) ArrayGet(ion->rr_rates, i);
      if (brts->rates->dim == 0) continue;
      iblk = brts->iblock;
      fblk = brts->fblock;
      sp_hdr.nele = ion->nele;
      sp_hdr.iblock = iblk->ib;
      sp_hdr.fblock = fblk->ib;
      sp_hdr.type = TransitionType(iblk->ncomplex, fblk->ncomplex);
      StrNComplex(sp_hdr.icomplex, iblk->ncomplex);
      StrNComplex(sp_hdr.fcomplex, fblk->ncomplex);
      InitFile(f, &fhdr, &sp_hdr);
      smax = 0.0;
      for (m = 0; m < brts->rates->dim; m++) {
	rt = (RATE *) ArrayGet(brts->rates, m);
	p = rt->i;
	q = rt->f;
	e = ion->energy[p] - ion->energy[q];
	j = ion->ilev[rt->i];
	if (iblk->n[j] > 0.0) {
	  r.lower = q;
	  r.upper = p;
	  r.energy = e;
	  r.strength = electron_density * iblk->n[j] * rt->dir;
	  s = r.strength * e;
	  if (s < strength_threshold*smax) continue;
	  if (s > smax) smax = s;
	  rx.sdev = 0.0;
	  r.rrate = rt->dir*electron_density;
	  r.trate = iblk->total_rate[j];
	  WriteSPRecord(f, &r, &rx);
	}
      }
      DeinitFile(f, &fhdr);
    }
  }  
  CloseFile(f, &fhdr);
  return 0;
}

static int CompareLine(const void *p1, const void *p2) {
  double *v1, *v2;
  v1 = (double *) p1;
  v2 = (double *) p2;
  if (*v1 < *v2) return -1;
  else if (*v1 > *v2) return 1;
  else return 0;
}
      
int SelectLines(char *ifn, char *ofn, int nele, int type, 
		double emin, double emax, double fmin) {
  F_HEADER fh;
  SP_HEADER h;
  SP_RECORD r, *rp;
  SP_EXTRA rx, *rpx;
  ARRAY sp, spx;
  ARRAY linetype;
  int *tt;
  FILE *f1, *f2;
  int n, nb, i;
  int t, t0, t1, t2;
  int r0, r1;
  int low, up;
  double e, a, smax;
  int swp;  
  
  rx.sdev = 0.0;

  f1 = fopen(ifn, "r");
  if (f1 == NULL) {
    printf("ERROR: File %s does not exist\n", ifn);
    return -1;
  }
  f2 = fopen(ofn, "a");
  if (f2 == NULL) {
    printf("ERROR: Cannot open file %s\n", ofn);
    return -1;
  }

  t2 = abs(type) / 1000000;
  if (type < 0) t2 = -1;
  t = abs(type) % 1000000;
  t1 = t / 10000;
  t0 = t % 10000;
  t0 = t0/100; 

  if (fmin < 0.0) {
    low = emin;
    up = emax;
  }
  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) return -1;

  ArrayInit(&sp, sizeof(SP_RECORD), 512);
  ArrayInit(&spx, sizeof(SP_EXTRA), 512);
  ArrayInit(&linetype, sizeof(int), 512);
  smax = 0.0;
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadSPHeader(f1, &h, swp);
    if (h.ntransitions == 0) continue;
    if (h.nele != nele) goto LOOPEND;  
    r1 = h.type / 10000;
    r0 = h.type % 10000;
    r0 = r0/100;
    if (type != 0) {
      if (t2 == 0) {
	if (t != h.type) goto LOOPEND;
      } else if (t2 == 1) {
	if (r1 < t1) goto LOOPEND;
	if (h.type < 100) goto LOOPEND;
	if (t%10000 != h.type%10000) goto LOOPEND;
      } else {
	if (t < 100) {
	  if (h.type > 99) goto LOOPEND;
	  if (h.type < t) goto LOOPEND;
	} else {
	  if (r1 < t1) goto LOOPEND;
	  if (r0 < t0) goto LOOPEND;
	}
      }
    }    
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      r.energy *= HARTREE_EV;
      rx.sdev *= HARTREE_EV;
      if (fmin < 0.0) {
	if (r.lower == low && r.upper == up) {
	  ArrayAppend(&sp, &r, NULL);
	  ArrayAppend(&spx, &rx, NULL);
	  ArrayAppend(&linetype, &(h.type), NULL);
	  break;
	}
      } else {
	if (r.energy < emin || r.energy > emax) continue;
	a = r.energy * r.strength;
	if (a < smax*fmin) continue;
	if (a > smax) smax = a;
	ArrayAppend(&sp, &r, NULL);
	ArrayAppend(&spx, &rx, NULL);
	ArrayAppend(&linetype, &(h.type), NULL);
      }
    }
    continue;

  LOOPEND:
    fseek(f1, h.length, SEEK_CUR);
  }

  if (fmin < 0.0) {
    if (sp.dim > 0) {
      rp = (SP_RECORD *) ArrayGet(&sp, 0);
      rpx = (SP_EXTRA *) ArrayGet(&spx, 0);
      tt = (int *) ArrayGet(&linetype, 0);
      fprintf(f2, "%2d %6d %6d %6d %13.6E %11.4E %11.4E\n", 
	      nele, rp->lower, rp->upper, *tt, rp->energy, rpx->sdev, rp->strength);
    }
  } else {
    if (sp.dim > 0) {
      smax *= fmin;
      for (i = 0; i < sp.dim; i++) {
	rp = (SP_RECORD *) ArrayGet(&sp, i);
	rpx = (SP_EXTRA *) ArrayGet(&spx, i);
	tt = (int *) ArrayGet(&linetype, i);
	e = rp->energy;
	if (rp->strength*e > smax) {
	  fprintf(f2, "%2d %6d %6d %6d %13.6E %11.4E %11.4E\n", 
		  nele, rp->lower, rp->upper, *tt, e, rpx->sdev, rp->strength);
	}
      }
    }
  }	
  ArrayFree(&sp, NULL);
  ArrayFree(&spx, NULL);
  ArrayFree(&linetype, NULL);

  fclose(f1);
  fclose(f2);
  
  return 0;
}
    
int PlotSpec(char *ifn, char *ofn, int nele, int type, 
	     double emin, double emax, double de0, double smin) {
  F_HEADER fh;
  SP_HEADER h;
  SP_RECORD r;
  SP_EXTRA rx;
  DISTRIBUTION *dist;
  FILE *f1, *f2;
  int n, nb, i;
  int t, t0, t1, t2;
  int r0, r1;
  double e;
  int m, k, nsp;
  double *sp, *tsp, *xsp, *kernel;
  double de10, de01;
  double a, sig, factor;
  double *lines;
  double smax, de, hc=12.3984E3;
  int swp;
  int idist;

  f1 = fopen(ifn, "r");
  if (f1 == NULL) {
    printf("ERROR: File %s does not exist\n", ifn);
    return -1;
  }
  f2 = fopen(ofn, "a");
  if (f2 == NULL) {
    printf("ERROR: Cannot open file %s\n", ofn);
    return -1;
  }

  rx.sdev = 0.0;

  t2 = abs(type) / 1000000;
  if (type < 0) t2 = -1;
  t = abs(type) % 1000000;
  t1 = t / 10000;
  t0 = t % 10000;
  t0 = t0 / 100;

  de = fabs(de0);
  de01 = 0.1*de;
  de10 = 10.0*de;
  sig = de/2.35;
  factor = 1.0/(sqrt(2*PI)*sig);
  sig = 1.0/(2*sig*sig);
  kernel = (double *) malloc(sizeof(double)*128);
  e = -63.5*de01;
  for (i = 0; i < 128; i++){
    kernel[i] = factor*exp(-sig*e*e);
    e += de01;
  }

  nsp = (emax - emin)/de01;
  sp = (double *) malloc(sizeof(double)*nsp);
  xsp = (double *) malloc(sizeof(double)*nsp);
  tsp = (double *) malloc(sizeof(double)*nsp);
  sp[0] = 0.0;
  tsp[0] = 0.0;
  xsp[0] = emin;
  for (i = 1; i < nsp; i++) {
    tsp[i] = 0.0;
    sp[i] = 0.0;
    xsp[i] = xsp[i-1] + de01;
  }

  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) return -1;

  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadSPHeader(f1, &h, swp);
    if (n == 0) break;
    if (h.ntransitions == 0) continue;
    if (h.nele != nele) goto LOOPEND; 
    r1 = h.type / 10000;
    r0 = h.type % 10000;
    r0 = r0/100;
    if (type != 0) {
      if (t2 == 0) {
	if (t != h.type) goto LOOPEND;
      } else if (t2 == 1) {
	if (r1 < t1) goto LOOPEND;
	if (h.type < 100) goto LOOPEND;
	if (t%10000 != h.type%10000) goto LOOPEND;
      } else {
	if (t < 100) {
	  if (h.type > 99) goto LOOPEND;
	  if (h.type < t) goto LOOPEND;
	} else {
	  if (r1 < t1) goto LOOPEND;
	  if (r0 < t0) goto LOOPEND;
	}
      }
    }
    m = 2*h.ntransitions;
    lines = (double *) malloc(sizeof(double)*m);  
    k = 0;
    smax = 0.0;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      e = r.energy;
      a = r.strength * e;
      if (a < smax*smin) continue;
      if (a > smax) smax = a;
      e *= HARTREE_EV;
      if (de0 < 0) e = hc/e;
      lines[k++] = e;
      lines[k++] = r.strength;
    }
    m = k;
	
    qsort(lines, m/2, sizeof(double)*2, CompareLine);
    k = 0;
    i = 0;
    while (k < m && i < nsp-1) {
      if (lines[k] < xsp[i]) {
	k += 2;
      } else if (lines[k] < xsp[i+1]) {
	sp[i] += lines[k+1];
	k += 2;
      } else {
	i++;
      }
    } 
    free(lines);
    for (i = 0; i < nsp; i++) {
      if (sp[i] > 0.0) {
	for (m = i-64, k = 0; k < 128; k++, m++) {
	  if (m > 0 && m < nsp) tsp[m] += sp[i]*kernel[k];
	}
      }
      sp[i] = 0.0;
    }
    continue;
  LOOPEND:
    fseek(f1, h.length, SEEK_CUR);
  }

  if (type != 0 && t < 100 && de0 > 0) {
    dist = GetEleDist(&idist);
    sig = dist->params[0];
    m = 10*sig/de01;
    if (idist < 2 && m > 9) {
      for (i = 0; i < nsp; i++) {
	sp[i] = 0.0;
      }
      free(kernel);
      kernel = (double *) malloc(sizeof(double)*m);
      if (idist == 0) {
	kernel[0] = 0.0;
	e = de01;
	for (i = 1; i < m; i++) {
	  kernel[i] = dist->dist(e, dist->params);
	  e += de01;
	}
	for (i = 0; i < nsp; i++) {
	  if (tsp[i] > 0.0) {
	    for (k = 0; k < m; k++) {
	      r0 = i + k;
	      if (r0 >= nsp) break;
	      sp[r0] += tsp[i]*kernel[k];
	    }
	  }
	}
      } else if (idist == 1) {
	r1 = m/2;
	e = -de01*r1;
	for (i = 0; i < m; i++) {
	  kernel[i] = dist->dist(e, dist->params);
	  e += de01;
	}
	for (i = 0; i < nsp; i++) {
	  if (tsp[i] > 0.0) {
	    for (k = 0; k < m; k++) {
	      r0 = i + k - r1;
	      if (r0 < 0 || r0 >= nsp) continue;
	      sp[r0] += tsp[i]*kernel[k];
	    }
	  }
	}
      }
      for (i = 0; i < nsp; i++) {
	fprintf(f2, "%15.8E\t%15.8E\n", xsp[i], sp[i]*de01);
      }
    } else {
      for (i = 0; i < nsp; i++) {
	fprintf(f2, "%15.8E\t%15.8E\n", xsp[i], tsp[i]);
      }
    }
  } else {
    for (i = 0; i < nsp; i++) {
      fprintf(f2, "%15.8E\t%15.8E\n", xsp[i], tsp[i]);
    }
  }

  fprintf(f2, "\n\n");

  free(xsp);
  free(sp);
  free(tsp);
  free(kernel);

  fclose(f1);
  fclose(f2);

  return 0;
}

int AddRate(ION *ion, ARRAY *rts, RATE *r, int m) {
  LBLOCK *ib, *fb;
  BLK_RATE *brt, brt0;
  RATE *r0;
  int i;
  
  ib = ion->iblock[r->i];
  fb = ion->iblock[r->f];
  for (i = 0; i < rts->dim; i++) {
    brt = (BLK_RATE *) ArrayGet(rts, i);
    if (brt->iblock == ib && brt->fblock == fb) {
      break;
    }
  }
  if (i == rts->dim) {
    brt0.iblock = ib;
    brt0.fblock = fb;
    brt0.rates = (ARRAY *) malloc(sizeof(ARRAY));
    ArrayInit(brt0.rates, sizeof(RATE), RATES_BLOCK);
    ArrayAppend(brt0.rates, r, NULL);
    ArrayAppend(rts, &brt0, InitBlkRateData);
  } else {
    if (m) {
      for (i = 0; i < brt->rates->dim; i++) {
	r0 = (RATE *) ArrayGet(brt->rates, i);
	if (r0->i == r->i && r0->f == r0->f) break;
      }
      if (i == brt->rates->dim) {
	ArrayAppend(brt->rates, r, NULL);
      } else {
	r0->dir += r->dir;
	r0->inv += r->inv;
	return 1;
      }
    } else {
      ArrayAppend(brt->rates, r, NULL);
    }
  }
  return 0;
}
  
int SetCERates(int inv) {
  int nb, i, j;
  int n, m, m1, k;
  int j1, j2;
  int p, q;
  ION *ion;
  RATE rt;
  F_HEADER fh;
  CE_HEADER h;
  CE_RECORD r;
  FILE *f;
  double e;
  float *cs;
  double data[2+(1+MAXNUSR)*2];
  double *y, *x;
  double *eusr;
  int swp;
  
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  y = data + 2;
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->ce_rates, FreeBlkRateData);
    f = fopen(ion->dbfiles[DB_CE-1], "r");
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_CE-1]);
      continue;
    }
    n = ReadFHeader(f, &fh, &swp);
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadCEHeader(f, &h, swp);
      eusr = h.usr_egrid;
      if (h.nele == ion->nele-1) {
	if (k > 0 || ion0.nionized > 0) {
	  fseek(f, h.length, SEEK_CUR);
	  continue;
	}
      }
      m = h.n_usr;
      m1 = m + 1;
      x = y + m1;
      x[m] = eusr[m-1]/(h.te0+eusr[m-1]);
      data[0] = h.te0*HARTREE_EV;
      for (j = 0; j < m; j++) {
	x[j] = log((h.te0 + eusr[j])/h.te0);
      }
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadCERecord(f, &r, swp, &h);
	rt.i = r.lower;
	rt.f = r.upper;
	j1 = ion->j[r.lower];
	j2 = ion->j[r.upper];
	e = ion->energy[r.upper] - ion->energy[r.lower];
	data[1] = r.bethe;
	cs = r.strength;
	y[m] = r.born[0];
	for (j = 0; j < m; j++) {
	  y[j] = cs[j];
	}
	CERate(&(rt.dir), &(rt.inv), inv, j1, j2, e, m,
	       data, rt.i, rt.f);
	AddRate(ion, ion->ce_rates, &rt, 0);
	if (h.qk_mode == QK_FIT) free(r.params);
	free(r.strength);
      }
      free(h.tegrid);
      free(h.egrid);
      free(h.usr_egrid);
    }
    fclose(f);
    
    if (k == 0 && ion0.nionized > 0) {
      f = fopen(ion0.dbfiles[DB_CE-1], "r");
      if (f == NULL) {
	printf("File %s does not exist, skipping.\n", ion0.dbfiles[DB_CE-1]);
	continue;
      }
      n = ReadFHeader(f, &fh, &swp);
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadCEHeader(f, &h, swp);
	eusr = h.usr_egrid;
	if (h.nele != ion0.nele) {
	  fseek(f, h.length, SEEK_CUR);
	  continue;
	}
	m = h.n_usr;
	m1 = m + 1;
	x = y + m1;
	x[m] = eusr[m-1]/(h.te0+eusr[m-1]);
	data[0] = h.te0*HARTREE_EV;
        for (j = 0; j < m; j++) {
	  x[j] = log((h.te0 + eusr[j])/h.te0);
        }
	for (i = 0; i < h.ntransitions; i++) {
	  n = ReadCERecord(f, &r, swp, &h);
	  p = IonizedIndex(r.lower, 0);
	  if (p < 0) {
	    if (h.qk_mode == QK_FIT) free(r.params);
	    free(r.strength);
	    continue;
	  }
	  q = IonizedIndex(r.upper, 0);
	  if (q < 0) {
	    if (h.qk_mode == QK_FIT) free(r.params);
	    free(r.strength);
	    continue;
	  }
	  rt.i = ion0.ionized_map[1][p];
	  rt.f = ion0.ionized_map[1][q];
	  j1 = ion->j[rt.i];
	  j2 = ion->j[rt.f];
	  e = ion0.energy[q] - ion0.energy[p];
	  data[1] = r.bethe;	
	  cs = r.strength;
	  y[m] = r.born[0];
	  for (j = 0; j < m; j++) {
	    y[j] = cs[j];
	  }
	  CERate(&(rt.dir), &(rt.inv), inv, j1, j2, e, m,
		 data, rt.i, rt.f);
	  AddRate(ion, ion->ce_rates, &rt, 0);
	  if (h.qk_mode == QK_FIT) free(r.params);
	  free(r.strength);
	}
	free(h.tegrid);
	free(h.egrid);
	free(h.usr_egrid);
      }
      fclose(f);
    }
  }

  return 0;
}

int SetTRRates(int inv) {
  int nb, i;
  int n, k;
  int j1, j2;
  int p, q, m;
  ION *ion;
  RATE rt;
  F_HEADER fh;
  TR_HEADER h;
  TR_RECORD r;
  TR_EXTRA rx;
  LBLOCK *ib;
  double e, gf;
  FILE *f;  
  int swp, iuta, im;

  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->tr_rates, FreeBlkRateData);
    f = fopen(ion->dbfiles[DB_TR-1], "r");
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_TR-1]);
      continue;
    }
    n = ReadFHeader(f, &fh, &swp);
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadTRHeader(f, &h, swp);
      iuta = IsUTA();
      if (h.nele == ion->nele-1) {
	if (k > 0 || ion0.nionized > 0) {
	  fseek(f, h.length, SEEK_CUR);
	  continue;
	}
      }
      if (abs(h.multipole) == 1) m = 0;
      else m = 1;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadTRRecord(f, &r, &rx, swp);
	rt.i = r.upper;
	if (ion0.n < 0) {
	  ib = ion->iblock[r.upper];
	  if (ib->rec &&
	      ib->rec->nrec[ib->irec] > 10) {
	    continue;
	  }
	}
	rt.f = r.lower;
	j1 = ion->j[r.upper];
	j2 = ion->j[r.lower];
	e = ion->energy[r.upper] - ion->energy[r.lower];
	if (iuta) e = rx.energy;
	gf = OscillatorStrength(h.multipole, e, (double)(r.strength), NULL);
	if (iuta) gf *= rx.sci;
	TRRate(&(rt.dir), &(rt.inv), inv, j1, j2, e, (float)gf);
	im = AddRate(ion, ion->tr_rates, &rt, m);
	if (iuta && im == 0) {
	  rt.dir = rx.energy;
	  rt.inv = rx.sdev;
	  AddRate(ion, ion->tr_sdev, &rt, 0);
	}
      }
    }
    fclose(f);
    if (ion->nele == 1) {
      ArrayFree(ion->tr2_rates, FreeBlkRateData);
      rt.f = FindLevelByName(ion->dbfiles[DB_EN-1], 1, 
			     "1*1", "1s1", "1s+1(1)1");
      rt.i = FindLevelByName(ion->dbfiles[DB_EN-1], 1,
			     "2*1", "2s1", "2s+1(1)1");
      if (rt.i >= 0 && rt.f >= 0) {
	rt.dir = TwoPhotonRate(ion0.atom, 0);
	rt.inv = 0.0;
	AddRate(ion, ion->tr2_rates, &rt, 0);
      }
    } else if (ion->nele == 2) {
      ArrayFree(ion->tr2_rates, FreeBlkRateData);
      rt.f = FindLevelByName(ion->dbfiles[DB_EN-1], 2, 
			     "1*2", "1s2", "1s+2(0)0");
      rt.i = FindLevelByName(ion->dbfiles[DB_EN-1], 2,
			     "1*1 2*1", "1s1 2s1", "1s+1(1)1 2s+1(1)0");
      if (rt.i >= 0 && rt.f >= 0) {
	rt.dir = TwoPhotonRate(ion0.atom, 1);
	rt.inv = 0.0;
	AddRate(ion, ion->tr2_rates, &rt, 0);
      }
      if (k == 0 && ion0.nionized > 0.0) {
	rt.f = FindLevelByName(ion->dbfiles[DB_EN-1], 1, 
			       "1*1", "1s1", "1s+1(1)1");
	rt.i = FindLevelByName(ion->dbfiles[DB_EN-1], 1,
			       "2*1", "2s1", "2s+1(1)1");
	if (rt.i >= 0 && rt.f >= 0) {
	  rt.dir = TwoPhotonRate(ion0.atom, 0);
	  rt.inv = 0.0;
	  AddRate(ion, ion->tr2_rates, &rt, 0);
	}
      }
    }
    if (ion0.n < 0.0) continue;
    ExtrapolateTR(ion, inv);
    if (k == 0 && ion0.nionized > 0) {
      f = fopen(ion0.dbfiles[DB_TR-1], "r");
      if (f == NULL) {
	printf("File %s does not exist, skipping.\n", ion0.dbfiles[DB_TR-1]);
	continue;
      }
      n = ReadFHeader(f, &fh, &swp);
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadTRHeader(f, &h, swp);
	iuta = IsUTA();
	if (h.nele != ion0.nele) {
	  fseek(f, h.length, SEEK_CUR);
	  continue;
	}  
	if (abs(h.multipole) == 1) m = 0;
	else m = 1;
	for (i = 0; i < h.ntransitions; i++) {
	  n = ReadTRRecord(f, &r, &rx, swp);
	  p = IonizedIndex(r.lower, 0);
	  if (p < 0) {
	    continue;
	  }
	  q = IonizedIndex(r.upper, 0);
	  if (q < 0) {
	    continue;
	  }
	  rt.i = ion0.ionized_map[1][q];
	  rt.f = ion0.ionized_map[1][p];
	  j1 = ion->j[rt.i];
	  j2 = ion->j[rt.f];
	  e = ion0.energy[q] - ion0.energy[p];
	  if (iuta) e = rx.energy;	    
	  gf = OscillatorStrength(h.multipole, e, (double)(r.strength), NULL);
	  if (iuta) gf *= rx.sci;
	  TRRate(&(rt.dir), &(rt.inv), inv, j1, j2, e, (float)gf);
	  im = AddRate(ion, ion->tr_rates, &rt, m);
	  if (iuta && im == 0) {
	    rt.dir = rx.energy;
	    rt.inv = rx.sdev;
	    AddRate(ion, ion->tr_sdev, &rt, 0);
	  }
	}
      }
      fclose(f);
    }
  }
  return 0;
}

int SetCIRates(int inv) { 
  int nb, i, t;
  int n, m, k;
  int j1, j2;
  ION *ion;
  RATE rt;
  F_HEADER fh;
  CI_HEADER h;
  CI_RECORD r;
  double e;
  FILE *f;  
  int swp;

  if (ion0.n < 0.0) return 0;

  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->ci_rates, FreeBlkRateData);
    f = fopen(ion->dbfiles[DB_CI-1], "r");
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_CI-1]);
      continue;
    }
    n = ReadFHeader(f, &fh, &swp);
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadCIHeader(f, &h, swp);
      m = h.nparams;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadCIRecord(f, &r, swp, &h);
	rt.i = r.b;
	rt.f = r.f;
	j1 = ion->j[r.b];
	j2 = ion->j[r.f];
	e = ion->energy[r.f] - ion->energy[r.b];
	CIRate(&(rt.dir), &(rt.inv), inv, j1, j2, e, m, r.params,
	       rt.i, rt.f);
	AddRate(ion, ion->ci_rates, &rt, 0);
	free(r.params);
	free(r.strength);
      }
      free(h.tegrid);
      free(h.egrid);
      free(h.usr_egrid);
    }
    fclose(f);
  }
  return 0;
}

int SetRRRates(int inv) { 
  int nb, i, j;
  int n, m, k;
  int j1, j2;
  ION *ion;
  RATE rt;
  F_HEADER fh;
  RR_HEADER h;
  RR_RECORD r;
  double e;
  FILE *f;  
  int swp;
  float *cs;
  double data[1+MAXNUSR*4];
  double *eusr;
  double *x, *logx, *y, *p;

  if (ion0.n < 0.0) return 0;
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  y = data + 1;
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->rr_rates, FreeBlkRateData);
    f = fopen(ion->dbfiles[DB_RR-1], "r");
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_RR-1]);
      continue;
    }
    n = ReadFHeader(f, &fh, &swp);
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadRRHeader(f, &h, swp);
      if (h.nparams <= 0) {
	printf("RR QkMode in %s must be in QK_FIT, nb=%d\n", 
	       ion->dbfiles[DB_RR-1], nb);
	exit(1);
      }
      eusr = h.usr_egrid;
      m = h.n_usr;
      x = y + m;
      logx = x + m;
      p = logx + m;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadRRRecord(f, &r, swp, &h);
	rt.i = r.f;
	rt.f = r.b;
	j1 = ion->j[r.f];
	j2 = ion->j[r.b];
	e = ion->energy[r.f] - ion->energy[r.b];
	data[0] = 3.5 + r.kl;
	if (e < 0.0) {
	  printf("%d %d %10.3E %10.3E\n", 
		 r.f, r.b, ion->energy[r.f],ion->energy[r.b]);
	  exit(1);
	}
	cs = r.strength;
	for (j = 0; j < m; j++) {
	  x[j] = (e+eusr[j])/e;
	  logx[j] = log(x[j]);
	  y[j] = log(cs[j]);
	}
	for (j = 0; j < h.nparams; j++) {
	  p[j] = r.params[j];
	}
	p[h.nparams-1] *= HARTREE_EV;
	RRRate(&(rt.dir), &(rt.inv), inv, j1, j2, e, m, data,
	       rt.i, rt.f);
	AddRate(ion, ion->rr_rates, &rt, 0);
	free(r.params);
	free(r.strength);
      }
      free(h.tegrid);
      free(h.egrid);
      free(h.usr_egrid);
    }
    fclose(f);
    ExtrapolateRR(ion, inv);
  }
  return 0;
}

int SetAIRatesInner(char *fn) {
  int nb, n, k, i, b0, nm;
  ION *ion;
  F_HEADER fh;
  AI_HEADER h;
  AI_RECORD r;
  FILE *f;  
  int swp;
  int ibase;

  if (inner_auger != 2) {
    printf("inner_auger must be 2 for this mode\n");
    return 0;
  }
  if (ion0.n < 0.0) return 0;
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  
  f = fopen(fn, "r");
  if (f == NULL) {
    printf("File %s does not exist, skipping.\n", fn);
    return 0;
  }

  n = ReadFHeader(f, &fh, &swp);
  b0 = -1;
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadAIHeader(f, &h, swp);
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadAIRecord(f, &r, swp);
      if (b0 < 0 || r.b < b0) b0 = r.b;
    }
    free(h.egrid);
  }

  fseek(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadAIHeader(f, &h, swp);
    for (k = 0; k < ions->dim; k++) {
      ion = (ION *) ArrayGet(ions, k);
      if (ion->nele == h.nele+1) break;
    }
    nm = ion->KLN_bmax - ion->KLN_bmin;    
    if (k < ions->dim) {
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadAIRecord(f, &r, swp);
	r.rate *= RATE_AU;
	ibase = r.b - b0;
	if (ibase >= 0 && ibase <= nm) {
	  ion->KLN_ai[ibase] += r.rate;
	}
      }    
    } else {
      fseek(f, h.length, SEEK_CUR);
    }
    free(h.egrid);
  }

  fclose(f);

  return 0;
}
  
int SetAIRates(int inv) {
  int nb, i, ib;
  int n, k;
  int j1, j2;
  ION *ion, *ion1;
  RATE rt;
  F_HEADER fh;
  AI_HEADER h;
  AI_RECORD r;
  double e;
  FILE *f;  
  int swp;
  int ibase;

  if (ion0.n < 0.0) return 0;

  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  for (k = 0; k < ions->dim; k++) {
    if (k == 0) ion = (ION *) ArrayGet(ions, k);
    else ion = ion1;
    if (k < ions->dim - 1) ion1 = (ION *) ArrayGet(ions, k+1);
    else ion1 = NULL;
    ArrayFree(ion->ai_rates, FreeBlkRateData);
    f = fopen(ion->dbfiles[DB_AI-1], "r");
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_AI-1]);
      continue;
    }
    n = ReadFHeader(f, &fh, &swp);
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadAIHeader(f, &h, swp);
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadAIRecord(f, &r, swp);
	if (inner_auger == 1) {
	  if (r.b <= ion->KLN_max && 
	      r.b >= ion->KLN_min &&
	      r.f <= ion->KLN_amax &&
	      r.f >= ion->KLN_amin) {
	    ibase = ion->ibase[r.b] - ion->KLN_bmin;
	    if (ibase >= 0) {
	      ion->KLN_ai[ibase] += r.rate*RATE_AU;
	    }
	  }
	} else if (inner_auger == 3) {
	  if (h.nele == ion->nele-1 &&
	      r.b <= ion->KLN_bmax && 
	      r.b >= ion->KLN_bmin) {
	    ibase = r.b - ion->KLN_bmin;	   
	    ion->KLN_ai[ibase] += r.rate*RATE_AU;
	    continue;
	  }
	} else if (inner_auger == 4) {
	  if (ion->iblock[r.b]->ionized) {
	    ib = IonIndex(ion1, ion->iblock[r.b]->ib, ion->ilev[r.b]);
	    if (ib <= ion1->KLN_bmax && ib >= ion1->KLN_bmin) {
	      ibase = ib - ion1->KLN_bmin;
	      ion1->KLN_ai[ibase] += r.rate*RATE_AU;
	    }
	  }
	}
	if (h.nele == ion->nele - 1) continue;
	rt.i = r.b;
	rt.f = r.f;
	j1 = ion->j[r.b];
	j2 = ion->j[r.f];
	e = ion->energy[r.b] - ion->energy[r.f];
	AIRate(&(rt.dir), &(rt.inv), inv, j1, j2, e, r.rate);
	AddRate(ion, ion->ai_rates, &rt, 0);
      }
      free(h.egrid);
    }
    if (inner_auger == 1) {
      n = ion->KLN_bmax - ion->KLN_bmin + 1;
      for (ibase = 0; ibase < n; ibase++) {
	if (ion->KLN_nai[ibase]) {
	  ion->KLN_ai[ibase] /= ion->KLN_nai[ibase];
	}
      }
    }
    fclose(f);
    ExtrapolateAI(ion, inv);
    
    if (inner_auger == 4 && k == 0 && ion0.nionized > 0) {
      f = fopen(ion0.dbfiles[DB_AI-1], "r");
      if (f == NULL) {
	printf("File %s does not exist, skipping.\n", ion0.dbfiles[DB_AI-1]);
	continue;
      }
      n = ReadFHeader(f, &fh, &swp);
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadAIHeader(f, &h, swp);
	for (i = 0; i < h.ntransitions; i++) {
	  n = ReadAIRecord(f, &r, swp);
	  ib = IonizedIndex(r.b, 0);
	  if (ib >= 0) {
	    ib = ion0.ionized_map[1][ib];
	    if (ib <= ion->KLN_bmax && ib >= ion->KLN_bmin) {
	      ibase = ib - ion->KLN_bmin;
	      ion->KLN_ai[ibase] += r.rate*RATE_AU;
	    }
	  }
	}
      }
      fclose(f);
    }
  }
  return 0;
}

int DRBranch(void) {
  ION *ion;
  RATE *r;
  LBLOCK *blk1, *blk2;
  BLK_RATE *brts;
  int i, k, m, t;
  int p, q;
  double a, d;

  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }

  for (k = 0; k < blocks->dim; k++) {
    blk1 = (LBLOCK *) ArrayGet(blocks, k);
    for (m = 0; m < blk1->nlevels; m++) {
      blk1->r[m] = 1.0;
      blk1->n[m] = 0.0;
    }
  }

  for (i = 1; i <= max_iter; i++) {
    for (k = 0; k < ions->dim; k++) {
      ion = (ION *) ArrayGet(ions, k);
      for (t = 0; t < ion->tr_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->tr_rates, t);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	for (m = 0; m < brts->rates->dim; m++) {
	  r = (RATE *) ArrayGet(brts->rates, m);
	  p = ion->ilev[r->i];
	  q = ion->ilev[r->f];
	  blk1->n[p] += blk2->r[q] * r->dir;
	}
      }
    }
  
    d = 0.0;
    for (k = 0; k < blocks->dim; k++) {
      blk1 = (LBLOCK *) ArrayGet(blocks, k);
      for (m = 0; m < blk1->nlevels; m++) {
	if (blk1->total_rate[m]) {
	  blk1->n[m] /= blk1->total_rate[m];
	  a = fabs(blk1->n[m] - blk1->r[m]);
	  if (blk1->n[m]) a /= blk1->n[m];
	  if (a > d) d = a;
	  blk1->r[m] = blk1->n[m];
	}
	blk1->n[m] = 0.0;
      }
    }
    printf("%5d %11.4E\n", i, d);
    if (d < iter_accuracy) break;
  }

  if (i == max_iter) {
    printf("Max iteration reached in DRBranch\n");
  }

  return 0;
}

int DRStrength(char *fn, int nele, int mode, int ilev0) {
  ION *ion;
  RATE *r, *rp;
  LBLOCK *blk1, *blk2;
  BLK_RATE *brts, *brtsp;
  DR_RECORD r1;
  DR_HEADER hdr;
  F_HEADER fhdr;
  int k, m, t, p, n, vnl, vn, vl;
  int mp, tp;
  FILE *f;
  
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }

  fhdr.type = DB_DR;
  fhdr.atom = ion0.atom;
  strcpy(fhdr.symbol, ion0.symbol);
  f = OpenFile(fn, &fhdr);
  
  if (ilev0 >= 0) {
    for (k = 0; k < ions->dim; k++) {
      ion = (ION *) ArrayGet(ions, k);
      if (ion->nele - 1 == nele) {
	ilev0 += ion->iground;
	break;
      }
    }
  } else {
    ilev0 = -ilev0;
  }

  hdr.ilev = ilev0;
  hdr.nele = nele;
  n = -1;
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    if (ion->nele - 1 == nele) {
      hdr.energy = ion->energy[ilev0];
      hdr.j = ion->j[ilev0];
      for (t = 0; t < ion->ai_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->ai_rates, t);
	blk1 = brts->fblock;
	if (ilev0 < blk1->imin || ilev0 >= blk1->imin+blk1->nlevels) {
	  continue;
	}
	blk1 = brts->iblock;
	for (m = 0; m < brts->rates->dim; m++) {
	  r = (RATE *) ArrayGet(brts->rates, m);
	  if (r->f != ilev0) continue;
	  vnl = ion->vnl[r->i];
	  vn = vnl/100;
	  vl = vnl - vn*100;
	  if (vn != n) {
	    if (n > 0) DeinitFile(f, &fhdr);
	    hdr.vn = vn;
	    InitFile(f, &fhdr, &hdr);
	    n = vn;
	  }	
	  p = ion->ilev[r->i];
	  r1.ilev = r->i;
	  r1.ibase = ion->ibase[r->i];
	  r1.energy = ion->energy[r->i] - ion->energy[ilev0];
	  r1.j = ion->j[r->i];
	  r1.vl = vl;
	  r1.ai = r->dir;
	  r1.total_rate = blk1->total_rate[p];
	  if (mode == 0) {
	    r1.etrans = 0.0;
	    r1.flev = -1;
	    r1.fbase = -1;
	    r1.br = blk1->r[p];
	    if (!(blk1->r[p])) continue;
	    WriteDRRecord(f, &r1);
	  } else if (mode == 1) {
	    for (tp = 0; tp < ion->tr_rates->dim; tp++) {
	      brtsp = (BLK_RATE *) ArrayGet(ion->tr_rates, tp);
	      blk2 = brtsp->iblock;
	      if (r1.ilev < blk2->imin || 
		  r1.ilev >= blk2->imin+blk2->nlevels) {
		continue;
	      }
	      for (mp = 0; mp < brtsp->rates->dim; mp++) {
		rp = (RATE *) ArrayGet(brtsp->rates, mp);
		if (rp->i != r1.ilev) continue;
		r1.etrans = ion->energy[rp->i] - ion->energy[rp->f];
		r1.flev = rp->f;
		r1.fbase = ion->ibase[rp->f];
		r1.br = rp->dir/r1.total_rate;
		if (!(rp->dir)) continue;
		WriteDRRecord(f, &r1);
	      }
	    }
	  } else if (mode == 2) {
	    for (tp = 0; tp < ion->ai_rates->dim; tp++) {
	      brtsp = (BLK_RATE *) ArrayGet(ion->ai_rates, tp);
	      blk2 = brtsp->iblock;
	      if (r1.ilev < blk2->imin || 
		  r1.ilev >= blk2->imin+blk2->nlevels) {
		continue;
	      }
	      for (mp = 0; mp < brtsp->rates->dim; mp++) {
		rp = (RATE *) ArrayGet(brtsp->rates, mp);
		if (rp->i != r1.ilev || !(rp->dir)) continue;
		r1.etrans = ion->energy[rp->i] - ion->energy[rp->f];
		r1.flev = rp->f;
		r1.fbase = ion->ibase[rp->f];
		r1.br = rp->dir/r1.total_rate;
		WriteDRRecord(f, &r1);
	      }
	    }
	  }
	}
      }
      break;
    }
  }

  if (n >= 0) DeinitFile(f, &fhdr);  
  CloseFile(f, &fhdr);

  return 0;
}

int DumpRates(char *fn, int k, int m, int imax, int a) {
  FILE *f;
  int i, t, p, q;
  short nele;
  double energy;
  ION *ion;
  ARRAY *rts;
  RATE *r;
  BLK_RATE *brts;
  
  for (p = 0; p < ions->dim; p++) {
    ion = (ION *) ArrayGet(ions, p);
    if (ion->nele != k) continue;
    f = fopen(fn, "w");
    if (f == NULL) {
      printf("cannot open file %s\n", fn);
      return -1;
    }
    if (m != 0) {
      switch (m) {
      case 1:
	rts = ion->tr_rates;
	break;
      case 2:
	rts = ion->tr2_rates;
	break;
      case 3:
	rts = ion->ce_rates;
	break;
      case 4:
	rts = ion->rr_rates;
	break;
      case 5:
	rts = ion->ai_rates;
	break;
      case 6:
	rts = ion->ci_rates;
	break;
      default:
	printf("invalid mode %d\n", m);
	fclose(f);
	return -1;
      }
      for (t = 0; t < rts->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(rts, t);
	for (q = 0; q < brts->rates->dim; q++) {
	  r = (RATE *) ArrayGet(brts->rates, q);
	  if (imax < 0 || (r->i <= imax && r->f <= imax)) {
	    if (a == 0) {
	      fwrite(&(r->i), sizeof(int), 1, f);
	      fwrite(&(r->f), sizeof(int), 1, f);
	      fwrite(&(r->dir), sizeof(double), 1, f);
	      fwrite(&(r->inv), sizeof(double), 1, f);
	    } else {
	      fprintf(f, "%7d %7d %10.3E %10.3E\n", 
		      r->i, r->f, r->dir, r->inv);
	    }
	  }
	}
      }
    } else {
      for (t = 0; t < ion->nlevels; t++) {	
	if (imax < 0 || t <= imax) {
	  if (ion->iblock[t] == NULL) continue;
	  q = ion->ilev[t];
	  energy = ion->energy[t];
	  if (p == ion->iblock[t]->iion) nele = ion->nele;
	  else if (p == ion->iblock[t]->iion + 1) nele = ion->nele - 1;
	  else nele = ion->nele + 1;
	  if (a == 0) {
	    fwrite(&(nele), sizeof(short), 1, f);
	    fwrite(&t, sizeof(int), 1, f);
	    fwrite(&(ion->iblock[t]->ib), sizeof(int), 1, f);
	    fwrite(&q, sizeof(int), 1, f);
	    fwrite(&(ion->j[t]), sizeof(short), 1, f);
	    fwrite(&(ion->ibase[t]), sizeof(short), 1, f);
	    fwrite(&(ion->vnl[t]), sizeof(short), 1, f);
	    fwrite(&(ion->energy[t]), sizeof(double), 1, f);
	    fwrite(&(ion->iblock[t]->total_rate[q]), sizeof(double), 1, f);
	    fwrite(&(ion->iblock[t]->r[q]), sizeof(double), 1, f);
	  } else {
	    fprintf(f, "%2d %6d %6d %6d %2d %4d %4d %15.8E %10.3E %10.3E\n", 
		    nele, t, ion->iblock[t]->ib, q, ion->j[t],
		    ion->ibase[t], ion->vnl[t], ion->energy[t],
		    ion->iblock[t]->total_rate[q],
		    ion->iblock[t]->r[q]);
	  }
	}
      }
    }    
    fclose(f);
  }
}

static void AddSpecBB(int nx, double *xg, double *yg, double e, double s, 
		      double dv, double trate) {
  int i;
  double a, u, b;

  a = trate*4.1328e-15/(4.0*PI*dv);
  s *= e*1.6e-12;
  for (i = 0; i < nx; i++) {
    u = (xg[i] - e)/dv;
    b = s*voigt(a,u)/dv;
    yg[i] += b;
  }
}

static void AddSpecBF(int nx, double *xg, double *yg, double e, double s, 
		      double te, double alpha, double a) {
  int i;
  double b, x;
  
  s *= a*1.6e-12;
  for (i = 0; i < nx; i++) {
    if (xg[i] <= e) continue;
    x = xg[i]-e;
    b = s*pow(x, alpha)*exp(-x/te)*xg[i];
    yg[i] += b;
  }
}

static void AddSpecFF(int nx, double *xg, double *yg, 
		      double ne, double ni, int z, double te) {
  int i;
  double b;

  for (i = 0; i < nx; i++) {
    b = BremssNR(z, te, xg[i]);
    b *= ne*ni;
    yg[i] += b;
  }
}

void TabNLTE(char *fn1, char *fn2, char *fn3, char *fn,
	     double xmin, double xmax, double dx) {
  FILE *f1, *f2, *f3, *f;
  char buf[1000];
  double *ab, *stot, *scol, *spho, *saut;
  double abt, *atot, *acol, *apho, *aaut;
  double pbb, pbf, pff, pfn, eint, eint1, te, te1;
  double ne, ni, gpbf, gcbf, gaut, gpbb, gcbb, rtot;
  double zbar, m2, m3, adx, *xg, *eg, *yg[3], tmp;
  int nmaxt, *nmax, i, j, t, k, z, n, *ilev;
  int swp1, swp2, swp3, m, nx, nions, nlevs;
  NCOMPLEX cmpx[MAXNCOMPLEX];
  F_HEADER fh1, fh2, fh3;
  SP_HEADER h1;
  SP_RECORD r1;
  SP_EXTRA rx;
  RT_HEADER h2;
  RT_RECORD r2, r3;
  double dv, emin, emax, a, alpha = 0.5;

  f1 = fopen(fn1, "r");
  f2 = fopen(fn2, "r");
  f3 = fopen(fn3, "r");
 
  ReadFHeader(f1, &fh1, &swp1);
  ReadFHeader(f2, &fh2, &swp2);  
  ReadFHeader(f3, &fh3, &swp3);
 
  z = (int) fh2.atom;
  nmax = malloc(sizeof(int)*(z+1));
  ilev = malloc(sizeof(int)*(z+1));
  ab = malloc(sizeof(double)*(z+1));
  stot = malloc(sizeof(double)*(z+1));
  scol = malloc(sizeof(double)*(z+1));
  spho = malloc(sizeof(double)*(z+1));
  saut = malloc(sizeof(double)*(z+1));
  atot = malloc(sizeof(double)*(z+1));
  acol = malloc(sizeof(double)*(z+1));
  apho = malloc(sizeof(double)*(z+1));
  aaut = malloc(sizeof(double)*(z+1));
  for (i = 0; i <= z; i++) {
    ab[i] = 0;
    stot[i] = 0;
    scol[i] = 0;
    spho[i] = 0;
    saut[i] = 0;
    atot[i] = 0;
    acol[i] = 0;
    apho[i] = 0;
    aaut[i] = 0;
    nmax[i] = 0;
    ilev[i] = 0;
  }
  
  pbb = 0.0;
  pbf = 0.0;
  pff = 0.0;
  pfn = 0.0;
  eint = 0.0;
  abt = 0.0;
  for (i = 0; i < fh2.nblocks; i++) {
    n = ReadRTHeader(f2, &h2, swp2);
    if (n == 0) break;
    n = ReadRTRecord(f2, &r2, swp2);
    if (n == 0) break;
    abt = r2.ai;
    break;
  }
  fseek(f2, 0, SEEK_SET);
  ReadFHeader(f2, &fh2, &swp2);

  sprintf(buf, "%s.elev", fn);
  f = fopen(buf, "w");
  fprintf(f, "energy_levels %12d\n", 0);
  for (j = 0; j < fh2.nblocks; j++) {
    n = ReadRTHeader(f2, &h2, swp2);
    if (n == 0) break;
    te = h2.p_edist[0];
    ne = h2.eden;
    gpbf = 0.0;
    gcbf = 0.0;
    gaut = 0.0;
    gpbb = 0.0;
    gcbb = 0.0;
    for (i = 0; i < h2.ntransitions; i++) {
      n = ReadRTRecord(f2, &r2, swp2);      
      if (n == 0) break;
      if (r2.ci < 0) {
	if (i > 0) {
	  rtot = gpbf + gcbf + gpbb + gcbb + gaut;
	  gpbf /= rtot;
	  gcbf /= rtot;
	  gpbb /= rtot;
	  gcbb /= rtot;
	  gaut /= rtot;
	  rtot /= r3.nb;
	  fprintf(f, "elev  %4d %7d %12.5E %12.5E %12.5E\n",
		  k, ilev[k], r3.tr, r3.ce, r3.nb/abt);
	  fprintf(f, "      %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
		  rtot, gcbb, gpbb, gcbf, gpbf, gaut);
	  fprintf(f, "      ");
	  GetNComplex(cmpx, r3.icomplex);
	  m = 0;
	  while (m < MAXNCOMPLEX && cmpx[m].n) {
	    if (m == 0) {
	      fprintf(f, "%2d ", cmpx[m].nq);
	    } else {
	      if (cmpx[m].n <= 6) {
		for (t = cmpx[m-1].n+1; t < cmpx[m].n; t++) {
		  fprintf(f, "%2d ", 0);
		}
		fprintf(f, "%2d ", cmpx[m].nq);
	      } else {
		for (t = cmpx[m-1].n+1; t <= 6; t++) {
		  fprintf(f, "%2d ", 0);
		}
	      }
	    }      
	    m++;
	  }
	  for (t = cmpx[m-1].n+1; t <= 6; t++) {
	    fprintf(f, "%2d ", 0);
	  }
	  fprintf(f, "%2d\n", cmpx[m-1].n);
	  if (nmax[k] < cmpx[m-1].n) nmax[k] = cmpx[m-1].n;
	  eint += r3.ce*r3.nb/abt;
	  pfn += r3.tr * exp(-r3.ce/te);
	}
	gpbf = 0.0;
	gcbf = 0.0;
	gaut = 0.0;
	gpbb = 0.0;
	gcbb = 0.0;
	memcpy(&r3, &r2, sizeof(RT_RECORD));
	k = (int) r2.rr;
	ilev[k]++;
	ab[k] += r2.nb/abt;
	r3.ce *= HARTREE_EV;
	continue;
      }
      if (r2.dir == 1 && r2.iblock == -1) {
	stot[k] += r2.rr + r2.ai + r2.ci;
	scol[k] += r2.ci;
	spho[k] += r2.rr;
	saut[k] += r2.ai;
	gpbf += r2.rr;
	gaut += r2.ai;
	gcbf += r2.ci;
      } else if (r2.dir == 1 && r2.iblock == -3) {
	atot[k] += r2.rr + r2.ai + r2.ci;
	acol[k] += r2.ci;
	apho[k] += r2.rr;
	aaut[k] += r2.ai;
	gpbf += r2.rr;
	gaut += r2.ai;
	gcbf += r2.ci;
      } else if (r2.dir == 1 && r2.iblock == -2) {
	gcbb += r2.ce;
	gpbb += r2.tr;
      }
    }
    rtot = gpbf + gcbf + gpbb + gcbb + gaut;
    gpbf /= rtot;
    gcbf /= rtot;
    gpbb /= rtot;
    gcbb /= rtot;
    gaut /= rtot;
    rtot /= r3.nb;
    fprintf(f, "elev  %4d %7d %12.5E %12.5E %12.5E\n",
	    k, ilev[k], r3.tr, r3.ce, r3.nb/abt);
    fprintf(f, "      %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
	    rtot, gcbb, gpbb, gcbf, gpbf, gaut);
    fprintf(f, "      ");
    GetNComplex(cmpx, r3.icomplex);
    m = 0;
    while (m < MAXNCOMPLEX && cmpx[m].n) {
      if (m == 0) {
	fprintf(f, "%2d ", cmpx[m].nq);
      } else {
	if (cmpx[m].n <= 6) {
	  for (t = cmpx[m-1].n+1; t < cmpx[m].n; t++) {
	    fprintf(f, "%2d ", 0);
	  }
	  fprintf(f, "%2d ", cmpx[m].nq);
	} else {
	  for (t = cmpx[m-1].n+1; t <= 6; t++) {
	    fprintf(f, "%2d ", 0);
	  }
	}
      }      
      m++;
    }
    for (t = cmpx[m-1].n+1; t <= 6; t++) {
      fprintf(f, "%2d ", 0);
    }
    fprintf(f, "%2d\n", cmpx[m-1].n);
    if (nmax[k] < cmpx[m-1].n) nmax[k] = cmpx[m-1].n;
    eint += r3.ce*r3.nb/abt;
    pfn += r3.tr * exp(-r3.ce/te);
  }
  fprintf(f, "\n");
  eint1 = 0.0;
  for (i = 0; i < fh3.nblocks; i++) {
    n = ReadRTHeader(f3, &h2, swp3);
    if (n == 0) break;    
    te1 = h2.p_edist[0];
    for (t = 0; t < h2.ntransitions; t++) {
      n = ReadRTRecord(f3, &r2, swp3);
      if (r2.ci < 0) {	
	eint1 += r2.ce*HARTREE_EV*r2.nb/abt;
      }
    }
  }  
  nions = 0;
  nlevs = 0;
  nmaxt = 0;
  zbar = 0;
  for (i = 0; i <= z; i++) {
    if (ab[i] > 0) {
      nions++;
      nlevs += ilev[i];
      if (nmax[i] > nmaxt) nmaxt = nmax[i];
      zbar += ab[i]*(z - i);
      if (stot[i] > 0) {
	scol[i] /= stot[i];
	spho[i] /= stot[i];
	saut[i] /= stot[i];
	stot[i] /= ab[i]*abt;
      }
      if (atot[i] > 0) {
	acol[i] /= atot[i];
	apho[i] /= atot[i];
	aaut[i] /= atot[i];
	atot[i] /= ab[i]*abt;
      }
    }
  }
  m2 = 0.0;
  m3 = 0.0;
  for (i = 0; i <= z; i++) {
    m2 += ab[i]*pow((z-i - zbar),2);
    m3 += ab[i]*pow((z-i - zbar),3);
  }
  fseek(f, 0, SEEK_SET);
  fprintf(f, "energy_levels %12d\n", nlevs);
  fclose(f);

  sprintf(buf, "%s.spec", fn);
  f = fopen(buf, "w");
  adx = fabs(dx);
  nx = (xmax - xmin)/adx + 1;
  xg = malloc(sizeof(double)*nx);
  eg = malloc(sizeof(double)*nx);
  xg[0] = xmin;
  for (j = 1; j < nx; j++) {
    xg[j] = xg[j-1] + adx;
  }
  for (i = 0; i < 3; i++) {
    yg[i] = malloc(sizeof(double)*nx);
    for (j = 0; j < nx; j++) {
      yg[i][j] = 0.0;
    }
  }
  if (dx < 0) {
    for (i = 0; i < nx; i++) {
      eg[i] = 12.398e3/xg[i];
    }
    emin = 12.3984e3/xmax;
    emax = 12.3984e3/xmin;
  } else {
    for (i = 0; i < nx; i++) {
      eg[i] = xg[i];
    }
    emin = xmin;
    emax = xmax;
  }
  ni = ne/zbar;
  for (i = 0; i <= z; i++) {
    if (ab[i] > 0) {
      pff += BremssNR(z-i, te, -1.0)*ne*ni*ab[i];
      AddSpecFF(nx, eg, yg[2], ne, ni*ab[i], z-i, te);
    }
  }
  dv = 2.0*te/((GetAtomicMassTable())[z]*9.38272e8);
  a = DLOGAM(alpha+1.0);
  a = 1.0/(exp(a)*pow(te, alpha+1.0));
  for (m = 0; m < fh1.nblocks; m++) {
    n = ReadSPHeader(f1, &h1, swp1);
    if (n == 0) break;
    for (t = 0; t < h1.ntransitions; t++) {
      n = ReadSPRecord(f1, &r1, &rx, swp1);
      if (n == 0) break;
      if (h1.type == 0) continue;
      r1.energy *= HARTREE_EV;
      r1.strength *= ni*1e10/abt;
      if (h1.type < 100) {
	pbf += r1.strength*(r1.energy+(alpha+1.0)*te)*1.6e-12;
	if (r1.energy > emax || emin-r1.energy > 50.0*te) continue;
	AddSpecBF(nx, eg, yg[1], r1.energy, r1.strength, te, alpha, a);
      } else if (h1.type >= 100) {
	pbb += r1.strength*r1.energy*1.6e-12;
	rx.sdev *= HARTREE_EV;
	rx.sdev = sqrt(rx.sdev*rx.sdev + dv*r1.energy*r1.energy);
	if (emin-r1.energy > 3.0*dv || r1.energy-emax > 3.0*dv) continue;
	AddSpecBB(nx, eg, yg[0], r1.energy, r1.strength, rx.sdev, r1.trate);
      }
    }
  }

  if (dx < 0) {
    for (t = 0; t < 3; t++) {
      for (i = 0; i < nx; i++) {
	yg[t][i] *= eg[i]/xg[i];
      }
    }
  }
    
  fprintf(f, "spectrum\t %2s %6d\n", fh1.symbol, nx);
  for (i = 0; i < nx; i++) {
    fprintf(f, "%12.5E %12.5E %12.5E %12.5E %12.5E\n",
	    xg[i], yg[0][i], yg[1][i], yg[2][i], yg[0][i]+yg[1][i]+yg[2][i]);
  }
  fprintf(f, "\n");
  fclose(f);

  sprintf(buf, "%s.ions", fn);
  f = fopen(buf, "w");
  fprintf(f, "data\t M.F. Gu, Stanford, FAC %d.%d.%d\n", 
	  fh1.version, fh1.sversion, fh1.ssversion);
  fprintf(f, "case\t %s\n", fn);
  fprintf(f, "code\t FAC\n");
  fprintf(f, "atom\t %s %d\n", fh1.symbol, (int)fh1.atom);
  fprintf(f, "calctime 0 0\n\n");
  
  fprintf(f, "summary_quantities\n");
  fprintf(f, "plasma\t %12.5E %12.5E\n", te, ne*1e10);
  fprintf(f, "time\t %12.5E\n", 0.0);
  fprintf(f, "zbar\t %12.5E\n", zbar);
  fprintf(f, "m2\t %12.5E\n", m2);
  fprintf(f, "m3\t %12.5E\n", m3);
  fprintf(f, "eint\t %12.5E\n", eint);
  fprintf(f, "deintdt\t %12.5E\n", (eint1-eint)/(te1-te));
  fprintf(f, "pfn\t %12.5E\n", pfn);
  fprintf(f, "nmax_eff  %2d\n", nmaxt);
  fprintf(f, "ploss\t %12.5E %12.5E %12.5E %12.5E\n\n",
	  pbb, pbf, pff, pbb+pbf+pff);
  
  fprintf(f, "ion_states\t %d\n", nions);
  for (i = 0; i <= z; i++) {
    if (ab[i] == 0.0) continue;
    fprintf(f, "ion   %12d %12.5E %12d\n", i, ab[i], nmax[i]);
    fprintf(f, "      %12.5E %12.5E %12.5E %12.5E\n", 
	    stot[i], scol[i], spho[i], saut[i]);
    fprintf(f, "      %12.5E %12.5E %12.5E %12.5E\n", 
	    atot[i], acol[i], apho[i], aaut[i]);
  }
  fprintf(f, "\n");
  fclose(f);

  fclose(f1);
  fclose(f2);
  fclose(f3);
  free(nmax);
  free(ilev);
  free(ab);
  free(stot);
  free(scol);
  free(spho);
  free(saut);
  free(atot);
  free(acol);
  free(apho);
  free(aaut);
  for (i = 0; i < 3; i++) {
    free(yg[i]);
  }
  free(xg);
  free(eg);
}
