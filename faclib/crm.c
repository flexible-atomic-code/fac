#include "crm.h"
#include "grid.h"

static char *rcsid="$Id: crm.c,v 1.6 2002/01/24 03:14:30 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static IONIZED ion0;
static ARRAY *ions;
static ARRAY *blocks;
static double *bmatrix = NULL;

static int n_single_blocks = 64;

static int max_iter = 256;
static double iter_accuracy = EPS2;
static double iter_stablizer = 0.75;

static double electron_density = EPS3;
static double photon_density = 0.0;

void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);
void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipvt,
	    double *b, int *ldb, int *info);


int SetNumSingleBlocks(int n) {
  n_single_blocks = 256;
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
  if (acc >= 0) iter_accuracy = acc;
  if (s > 0.0 && s < 1.0) iter_stablizer = s;
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

static void FreeIonData(void *p) {
  ION *ion;
  int i;

  ion = (ION *) p;
  if (ion->nlevels > 0) {
    free(ion->iblock);
    free(ion->ilev);
    free(ion->j);
    free(ion->energy);
    ion->nlevels = 0;
  }
  for (i = 0; i < NDB; i++) {
    free(ion->dbfiles[i]);
  }

  ArrayFree(ion->ce_rates, NULL);
  free(ion->ce_rates);
  ArrayFree(ion->tr_rates, NULL);
  free(ion->tr_rates);
  ArrayFree(ion->ci_rates, NULL);
  free(ion->ci_rates);
  ArrayFree(ion->rr_rates, NULL);
  free(ion->rr_rates);
  ArrayFree(ion->ai_rates, NULL);
  free(ion->ai_rates);
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
  int i;

  if (m < 0) return 0;

  ReinitDBase(0);
  if (m > 0) return 0;

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
  ArrayInit(ion.ce_rates, sizeof(RATE), RATES_BLOCK);
  ion.tr_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.tr_rates, sizeof(RATE), RATES_BLOCK);
  ion.ci_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.ci_rates, sizeof(RATE), RATES_BLOCK);
  ion.rr_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.rr_rates, sizeof(RATE), RATES_BLOCK);
  ion.ai_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.ai_rates, sizeof(RATE), RATES_BLOCK);
  
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

  ArrayAppend(ions, &ion);
  
  return ions->dim;
  
}

int SetBlocks(double ni, char *ifn) {
  ION *ion, *ion1 = NULL;
  F_HEADER fh;
  EN_HEADER h;
  EN_RECORD r, *r0, *r1;
  LBLOCK blk;
  FILE *f;
  int n, i, k, nb, nlevels;
  char *fn;
  int p, q = -1;
  int nionized, n0;
  int swp;

  ion0.n = ni;
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
      default:
	ion0.dbfiles[i] = NULL;
      }
    } else {
      ion0.dbfiles[i] = NULL;
    }
  }
  
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
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
    n = fread(&fh, sizeof(F_HEADER), 1, f);
    if (CheckEndian(&fh) != CheckEndian(NULL)) {
      swp = 1;
      SwapEndianFHeader(&fh);
    } else {
      swp = 0;
    }
    if (k == 0) {
      ion0.atom = fh.atom;
      strcpy(ion0.symbol, fh.symbol);
    }

    nlevels = 0;
    nionized = 0;
    for (i = 0; i < fh.nblocks; i++) {
      n = fread(&h, sizeof(EN_HEADER), 1, f);
      if (swp) SwapEndianENHeader(&h);
      nlevels += h.nlevels;
      if (h.nele == ion->nele-1) {
	nionized += h.nlevels;
      }
      fseek(f, h.length, SEEK_CUR);
    }
    ion->nlevels = nlevels;
    ion->iblock = (int *) malloc(sizeof(int)*nlevels);
    ion->ilev = (int *) malloc(sizeof(int)*nlevels);
    ion->j = (short *) malloc(sizeof(short)*nlevels);
    ion->energy = (double *) malloc(sizeof(double)*nlevels);
 
    if (k == 0 && ifn) {
      ion0.nionized = nionized;
      ion0.ionized_map[0] = (int *) malloc(sizeof(int)*nionized);
      ion0.ionized_map[1] = (int *) malloc(sizeof(int)*nionized);
      ion0.energy = (double *) malloc(sizeof(double)*nionized);
    }
    
    fseek(f, sizeof(F_HEADER), SEEK_SET);
    n0 = 0;
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = fread(&h, sizeof(EN_HEADER), 1, f);
      if (swp) SwapEndianENHeader(&h);
      if (h.nele == ion->nele) {
	fseek(f, h.length, SEEK_CUR);
	continue;
      } else if (h.nele == ion->nele-1) {
	r0 = (EN_RECORD *) malloc(sizeof(EN_RECORD)*h.nlevels);
	for (i = 0; i < h.nlevels; i++) {
	  n = fread(&r0[i], sizeof(EN_RECORD), 1, f);
	  if (swp) SwapEndianENRecord(&(r0[i]));
	}
	if (ifn) {
	  r1 = (EN_RECORD *) malloc(sizeof(EN_RECORD)*h.nlevels);
	  nlevels = FindLevelBlock(h.nlevels, r0, r1, ion->nele-1, ifn); 
	  if (nlevels != h.nlevels) {
	    printf("ERROR: Ionized block of ion %d ", ion->nele);
	    printf("does not match a block in file %s\n", ifn);
	    exit(1);
	  }
	}
	if (k > 0) {
	  for (i = 0; i < h.nlevels; i++) {
	    p = r0[i].ilev;
	    q = r1[i].ilev;
	    ion->iblock[p] = ion1->iblock[q];
	    ion->ilev[p] = ion1->ilev[q];
	    ion->j[p] = r0[i].j;
	    ion->energy[p] = r0[i].energy;
	  }
	} else {
	  blk.iion = -1;
	  blk.nlevels = h.nlevels;
	  blk.n = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.n0 = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.r = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.total_rate = (double *) malloc(sizeof(double)*blk.nlevels);
	  ArrayAppend(blocks, &blk);
	  q = -1;
	  for (i = 0; i < h.nlevels; i++) {
	    p = r0[i].ilev;
	    q++;
	    ion->iblock[p] = blocks->dim-1;
	    ion->ilev[p] = q;
	    ion->j[p] = r0[i].j;
	    ion->energy[p] = r0[i].energy;
	    if (ifn) {
	      ion0.ionized_map[0][n0] = r1[i].ilev;
	      ion0.ionized_map[1][n0] = r0[i].ilev;
	      ion0.energy[n0] = r1[i].energy;
	      n0++;
	    }
	  }
	}
	free(r0);
	if (ifn) {
	  free(r1);
	}
      } else {
	printf("ERROR: Ion charge state does not match\n");
	exit(1);
      }
    }
  
    fseek(f, sizeof(F_HEADER), SEEK_SET);
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = fread(&h, sizeof(EN_HEADER), 1, f);
      if (swp) SwapEndianENHeader(&h);
      if (h.nele != ion->nele) {
	fseek(f, h.length, SEEK_CUR);
	continue;
      }	
      if (nb > 0) {
	blk.iion = k;
	blk.nlevels = h.nlevels;
	blk.n = (double *) malloc(sizeof(double)*blk.nlevels);
	blk.n0 = (double *) malloc(sizeof(double)*blk.nlevels);
	blk.r = (double *) malloc(sizeof(double)*blk.nlevels);
	blk.total_rate = (double *) malloc(sizeof(double)*blk.nlevels);
	ArrayAppend(blocks, &blk);
	q = -1;
      }
      for (i = 0; i < h.nlevels; i++) {
	n = fread(&r, sizeof(EN_RECORD), 1, f);
	if (swp) SwapEndianENRecord(&r);
	if (nb == 0) {
	  if (i < n_single_blocks) {
	    blk.iion = k;
	    blk.nlevels = 1;
	    blk.n = (double *) malloc(sizeof(double)*blk.nlevels);
	    blk.n0 = (double *) malloc(sizeof(double)*blk.nlevels);
	    blk.r = (double *) malloc(sizeof(double)*blk.nlevels);
	    blk.total_rate = (double *) malloc(sizeof(double)*blk.nlevels);
	    ArrayAppend(blocks, &blk);
	    q = -1;	    
	  } else if (i == n_single_blocks) {
	    blk.iion = k;
	    blk.nlevels = h.nlevels - n_single_blocks;
	    blk.n = (double *) malloc(sizeof(double)*blk.nlevels);
	    blk.n0 = (double *) malloc(sizeof(double)*blk.nlevels);
	    blk.r = (double *) malloc(sizeof(double)*blk.nlevels);
	    blk.total_rate = (double *) malloc(sizeof(double)*blk.nlevels);
	    ArrayAppend(blocks, &blk);
	    q = -1;
	  }
	}
	p = r.ilev;
	q++;
	ion->iblock[p] = blocks->dim-1;
	ion->ilev[p] = q;
	ion->j[p] = r.j;
	ion->energy[p] = r.energy;
      }
    }
    ion1 = ion;
    fclose(f);
  }
  
  if (bmatrix) free(bmatrix);
  k = blocks->dim;
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
  FILE *f;
  int i, k, nr, nb;
  int swp;

  f = fopen(ifn, "r");
  if (f == NULL) {
    printf("File %s does not exist\n", ifn);
    return -1;
  }
  
  nr = fread(&fh, sizeof(F_HEADER), 1, f);
  if (CheckEndian(&fh) != CheckEndian(NULL)) {
    swp = 1;
    SwapEndianFHeader(&fh);
  } else {
    swp = 0;
  }
  k = 0;
  for (nb = 0; nb < fh.nblocks; nb++) {
    nr = fread(&h, sizeof(EN_HEADER), 1, f);
    if (swp) SwapEndianENHeader(&h);
    if (h.nele != nele) {
      fseek(f, h.length, SEEK_CUR);
      continue;
    }
    k = 0;
    for (i = 0; i < h.nlevels; i++) {
      nr = fread(&r1[k], sizeof(EN_RECORD), 1, f);
      if (swp) SwapEndianENRecord(&(r1[k]));
      if (strcmp(r1[k].ncomplex, r0[0].ncomplex) == 0) {
	k++;
	if (k == n) break;
      }
    }
    if (k == n) break;
  }

  if (k < n) return -1;

  qsort(r0, n, sizeof(EN_RECORD), CompareENRecord);
  qsort(r1, n, sizeof(EN_RECORD), CompareENRecord);

  fclose(f);

  return n;
}

int IonizedIndex(int i, int m) {
  int k;

  for (k = 0; k < ion0.nionized; k++) {
    if (ion0.ionized_map[m][k] == i) {
      return k;
    }
  }
  
  return -1;
}

int SetAbund(int nele, double abund) {
  ION *ion;
  int i;

  if (ion0.nele == nele) ion0.n = abund;
  else {
    for (i = 0; i < ions->dim; i++) {
      ion = (ION *) ArrayGet(ions, i);
      if (ion->nele = nele) {
	ion->n = abund;
	break;
      }
    }
  }
  
  return 0;
}

int InitBlocks(void) {
  ION  *ion;
  RATE *r;
  LBLOCK *blk;
  int k, m, i, j;
  double a, b;

  for (i = 0; i < blocks->dim; i++) {
    blk = (LBLOCK *) ArrayGet(blocks, i);
    k = blk->iion;
    if (k == -1) k = 0;
    ion = (ION *) ArrayGet(ions, k);
    blk->nb = 1.0;
    for (k = 0; k < blk->nlevels; k++) {
      blk->n0[k] = 0.0;
      blk->n[k] = 0.0;
      blk->r[k] = 0.0;
      blk->total_rate[k] = 0.0;
    }
    blk->r[0] = 1.0;
  }

  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    if (electron_density > 0.0) {
      for (m = 0; m < ion->ce_rates->dim; m++) {
	r = (RATE *) ArrayGet(ion->ce_rates, m);
	i = ion->iblock[r->i];
	j = ion->ilev[r->i];
	blk = (LBLOCK *) ArrayGet(blocks, i);
	blk->total_rate[j] += electron_density * r->dir;
	if (r->inv > 0.0) {
	  i = ion->iblock[r->f];
	  j = ion->ilev[r->f];
	  blk = (LBLOCK *) ArrayGet(blocks, i);
	  blk->total_rate[j] += electron_density * r->inv;
	}
      }
    }

    for (m = 0; m < ion->tr_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->tr_rates, m);
      i = ion->iblock[r->i];
      j = ion->ilev[r->i];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      blk->total_rate[j] += r->dir;
      if (r->inv > 0.0 && photon_density > 0.0) {
	a = photon_density * r->inv;
	b = a * (ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
	blk->total_rate[j] += b;
	i = ion->iblock[r->f];
	j = ion->ilev[r->f];
	blk = (LBLOCK *) ArrayGet(blocks, i);
	blk->total_rate[j] += a;
      }
    } 
     
    for (m = 0; m < ion->rr_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->rr_rates, m);
      i = ion->iblock[r->i];
      j = ion->ilev[r->i];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      if (electron_density > 0.0) {
	blk->total_rate[j] += electron_density * r->dir;
      }
      if (r->inv > 0.0 && photon_density > 0.0) {
	i = ion->iblock[r->f];
	j = ion->ilev[r->f];
	blk = (LBLOCK *) ArrayGet(blocks, i);
	blk->total_rate[j] += photon_density * r->inv;
      }
    } 

    for (m = 0; m < ion->ai_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->ai_rates, m);
      i = ion->iblock[r->i];
      j = ion->ilev[r->i];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      blk->total_rate[j] += r->dir;
      if (r->inv > 0.0 && electron_density > 0.0) {
	i = ion->iblock[r->f];
	j = ion->ilev[r->f];
	blk = (LBLOCK *) ArrayGet(blocks, i);
	blk->total_rate[j] += electron_density * r->inv;
      }
    } 

    if (electron_density > 0.0) {
      for (m = 0; m < ion->ci_rates->dim; m++) {
	r = (RATE *) ArrayGet(ion->ci_rates, m);
	i = ion->iblock[r->i];
	j = ion->ilev[r->i];
	blk = (LBLOCK *) ArrayGet(blocks, i);
	blk->total_rate[j] += electron_density * r->dir;
	if (r->inv > 0.0) {
	  i = ion->iblock[r->f];
	  j = ion->ilev[r->f];
	  blk = (LBLOCK *) ArrayGet(blocks, i);
	  blk->total_rate[j] += electron_density * r->inv;
	}
      }
    }
  }

  return 0;
}

int RateTable(char *fn) { 
  RT_RECORD rt;
  RT_HEADER rt_hdr;
  F_HEADER fhdr;
  ION *ion;
  RATE *r;
  LBLOCK *blk, *blk1;
  DISTRIBUTION *edist, *pdist;
  int n, k, m, i, j, p;
  double den;
  double *ce, *tr, *rr, *ci, *ai;
  FILE *f;

  edist = GetEleDist(&i);
  pdist = GetPhoDist(&j);
  fhdr.type = DB_RT;
  fhdr.atom = ion0.atom;
  strcpy(fhdr.symbol, ion0.symbol);
  rt_hdr.iedist = i;
  rt_hdr.ipdist = j;
  rt_hdr.np_edist = edist->nparams;
  rt_hdr.np_pdist = pdist->nparams;
  rt_hdr.p_edist = edist->params;
  rt_hdr.p_pdist = pdist->params;

  n = blocks->dim;
  k = n*n;
  ce = (double *) malloc(sizeof(double)*k);
  tr = (double *) malloc(sizeof(double)*k);
  ci = (double *) malloc(sizeof(double)*k);
  rr = (double *) malloc(sizeof(double)*k);
  ai = (double *) malloc(sizeof(double)*k);
  for (i = 0; i < k; i++) {
    ce[i] = 0.0;
    tr[i] = 0.0;
    ci[i] = 0.0;
    rr[i] = 0.0;
    ai[i] = 0.0;
  }

  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    for (m = 0; m < ion->ce_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->ce_rates, m);
      i = ion->iblock[r->i];
      j = ion->iblock[r->f];
      if (i == j) continue;
      blk = (LBLOCK *) ArrayGet(blocks, i);
      den = blk->r[ion->ilev[r->i]];
      if (den > 0.0) {
	p = i*n + j;
	ce[p] += den * r->dir;
      }
      if (r->inv > 0.0) {
	blk = (LBLOCK *) ArrayGet(blocks, j);
	den = blk->r[ion->ilev[r->f]];
	if (den > 0.0) {
	  p = i + j*n;
	  ce[p] += den * r->inv;
	}
      }
    }

    for (m = 0; m < ion->tr_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->tr_rates, m);
      i = ion->iblock[r->i];
      j = ion->iblock[r->f];
      if (i == j) continue;
      blk = (LBLOCK *) ArrayGet(blocks, i);
      den = blk->r[ion->ilev[r->i]];
      if (den > 0.0) {
	p = i*n + j;
	tr[p] += den * r->dir;
      }
      if (r->inv > 0.0) {
	blk = (LBLOCK *) ArrayGet(blocks, j);
	den = blk->r[ion->ilev[r->f]];
	if (den > 0.0) {
	  p = i + j*n;
	  tr[p] += den * r->inv;
	}
      }
    }

    for (m = 0; m < ion->rr_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->rr_rates, m); 
      i = ion->iblock[r->i];
      j = ion->iblock[r->f];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      den = blk->r[ion->ilev[r->i]];
      if (den > 0.0) {
	p = i*n + j;
	rr[p] += den * r->dir;
      }
      if (r->inv > 0.0) {
	blk = (LBLOCK *) ArrayGet(blocks, j);
	den = blk->r[ion->ilev[r->f]];
	if (den > 0.0) {
	  p = i + j*n;
	  rr[p] += den * r->inv;
	}
      }
    }

    for (m = 0; m < ion->ai_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->ai_rates, m); 
      i = ion->iblock[r->i];
      j = ion->iblock[r->f];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      den = blk->r[ion->ilev[r->i]];
      if (den > 0.0) {
	p = i*n + j;
	ai[p] += den * r->dir;
      }
      if (r->inv > 0.0) {
	blk = (LBLOCK *) ArrayGet(blocks, j);
	den = blk->r[ion->ilev[r->f]];
	if (den > 0.0) {
	  p = i + j*n;
	  ai[p] += den * r->inv;
	}
      }
    }  

    for (m = 0; m < ion->ci_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->ci_rates, m); 
      i = ion->iblock[r->i];
      j = ion->iblock[r->f];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      den = blk->r[ion->ilev[r->i]];
      if (den > 0.0) {
	p = i*n + j;
	ci[p] += den * r->dir;
      }
      if (r->inv > 0.0) {
	blk = (LBLOCK *) ArrayGet(blocks, j);
	den = blk->r[ion->ilev[r->f]];
	if (den > 0.0) {
	  p = i + j*n;
	  ci[p] += den * r->inv;
	}
      }  
    }
  }

  p = 0;
  for (i = 0; i < n; i++) {
    blk = (LBLOCK *) ArrayGet(blocks, i);
    k = blk->iion;
    if (k < 0) {
      m = ion0.nele;
    } else {
      ion = (ION *) ArrayGet(ions, k);
      m = ion->nele;
    }
    rt_hdr.nele = m;
    rt_hdr.iblock = i;
    rt_hdr.nb = blk->nb;
    f = InitFile(fn, &fhdr, &rt_hdr);
    for (j = 0; j < n; j++) {
      blk1 = (LBLOCK *) ArrayGet(blocks, j);
      if (tr[p] + 1.0 == 1.0 &&
	  ce[p] + 1.0 == 1.0 &&
	  ci[p] + 1.0 == 1.0 &&
	  rr[p] + 1.0 == 1.0 &&
	  ai[p] + 1.0 == 1.0) {
	p++;
	continue;
      }
      rt.iblock = j;
      rt.ce = ce[p];
      rt.tr = tr[p];
      rt.rr = rr[p];
      rt.ci = ci[p];
      rt.ai = ai[p];
      WriteRTRecord(f, &rt);
      p++;
    }
    CloseFile(f, &fhdr);
  }

  free(ce);
  free(tr);
  free(ci);
  free(rr);
  free(ai);

  return 0;
}

int BlockMatrix(void) {
  ION *ion;
  RATE *r;
  LBLOCK *blk;
  int n, k, m, i, j;
  int p, q, iion, k0, k1;
  double *x, den, a;
  
  n = blocks->dim;
  for (i = 0; i < n*(n+1); i++) {
    bmatrix[i] = 0.0;
  }
  x = bmatrix + n*n;

  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    if (electron_density > 0.0) {
      for (m = 0; m < ion->ce_rates->dim; m++) {
	r = (RATE *) ArrayGet(ion->ce_rates, m);
	i = ion->iblock[r->i];
	j = ion->iblock[r->f];
	if (i == j) continue;
	blk = (LBLOCK *) ArrayGet(blocks, i);
	den = blk->r[ion->ilev[r->i]];
	if (den > 0.0) {
	  p = i*n + j;
	  bmatrix[p] += den * electron_density * r->dir;
	}
	if (r->inv > 0.0) {
	  blk = (LBLOCK *) ArrayGet(blocks, j);
	  den = blk->r[ion->ilev[r->f]];
	  if (den > 0.0) {
	    p = i + j*n;
	    bmatrix[p] += den * electron_density * r->inv;
	  }
	}
      }
    }

    for (m = 0; m < ion->tr_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->tr_rates, m);
      i = ion->iblock[r->i];
      j = ion->iblock[r->f];
      if (i == j) continue;
      blk = (LBLOCK *) ArrayGet(blocks, i);
      den = blk->r[ion->ilev[r->i]];
      if (den > 0.0) {
	p = i*n + j;
	bmatrix[p] += den * r->dir;
      }
      if (r->inv > 0.0 && photon_density > 0.0) {
	if (den > 0.0) {
	  a = photon_density * r->inv;
	  bmatrix[p] += den*a*(ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
	}
	blk = (LBLOCK *) ArrayGet(blocks, j);
	den = blk->r[ion->ilev[r->f]];
	if (den > 0.0) {
	  p = i + j*n;
	  bmatrix[p] += den * photon_density * r->inv;
	}
      }
    }

    for (m = 0; m < ion->rr_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->rr_rates, m); 
      i = ion->iblock[r->i];
      j = ion->iblock[r->f];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      den = blk->r[ion->ilev[r->i]];
      if (den > 0.0) {
	if (electron_density > 0.0) {
	  p = i*n + j;
	  bmatrix[p] += den * electron_density * r->dir;
	}
      }
      if (r->inv > 0.0 && photon_density > 0.0) {
	blk = (LBLOCK *) ArrayGet(blocks, j);
	den = blk->r[ion->ilev[r->f]];
	if (den > 0.0) {
	  p = i + j*n;
	  bmatrix[p] += den * photon_density * r->inv;
	}
      }
    }

    for (m = 0; m < ion->ai_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->ai_rates, m); 
      i = ion->iblock[r->i];
      j = ion->iblock[r->f];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      den = blk->r[ion->ilev[r->i]];
      if (den > 0.0) {
	p = i*n + j;
	bmatrix[p] += den * r->dir;
      }
      if (r->inv > 0.0 && electron_density > 0.0) {
	blk = (LBLOCK *) ArrayGet(blocks, j);
	den = blk->r[ion->ilev[r->f]];
	if (den > 0.0) {
	  p = i + j*n;
	  bmatrix[p] += den * electron_density * r->inv;
	}
      }
    }  

    if (electron_density > 0.0) {
      for (m = 0; m < ion->ci_rates->dim; m++) {
	r = (RATE *) ArrayGet(ion->ci_rates, m); 
	i = ion->iblock[r->i];
	j = ion->iblock[r->f];
	blk = (LBLOCK *) ArrayGet(blocks, i);
	den = blk->r[ion->ilev[r->i]];
	if (den > 0.0) {
	  p = i*n + j;
	  bmatrix[p] += den * electron_density * r->dir;
	}
	if (r->inv > 0.0) {
	  blk = (LBLOCK *) ArrayGet(blocks, j);
	  den = blk->r[ion->ilev[r->f]];
	  if (den > 0.0) {
	    p = i + j*n;
	    bmatrix[p] += den * electron_density * r->inv;
	  }
	}
      }  
    }
  }

  for (i = 0; i < n; i++) {
    blk = (LBLOCK *) ArrayGet(blocks, i);
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
    blk = (LBLOCK *) ArrayGet(blocks, i);
    if (blk->iion != iion) {
      if (iion != -2) {
	k = iion;
	if (k == -1) k = 0;
	ion = (ION *) ArrayGet(ions, k);
	if (iion == -1) den = ion0.n;
	else den = ion->n;
	if (den > 0.0) {
	  x[k0] = den;
	  p = k0;
	  for (k = 0; k < n; k++) {
	    if (k < k1 && k >= k0) bmatrix[p] = 1.0;
	    else bmatrix[p] = 0.0;
	    p += n;
	  }
	} 
      }
      iion = blk->iion;
      k0 = k1;
    }
    k1++;
  }

  k = iion;
  ion = (ION *) ArrayGet(ions, k);
  den = ion->n;
  if (den > 0.0) {
    x[k0] = den;
    p = k0;
    for (k = 0; k < n; k++) {
      if (k < k1 && k >= k0) bmatrix[p] = 1.0;
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
    if (1.0+bmatrix[i+i*n] == 1.0) {
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

  dgesv_(&m, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  if (info != 0) {
    printf("Error in solving BlockMatrix\n");
    exit(1);
  }
  
  p = 0;
  for (i = 0; i < n; i++) {
    blk = (LBLOCK *) ArrayGet(blocks, i);
    if (x[i] < 0.0) {
      blk->nb = 0.0;
      for (j = 0; j < blk->nlevels; j++) {
	blk->r[j] = 0.0;
	blk->n[j] = 0.0;
      }
    } else {
      blk->nb = b[p++];
      for (j = 0; j < blk->nlevels; j++) {
	blk->r[j] *= blk->nb;
	blk->n[j] = 0.0;
      }
    }
  }

  return 0;
}
  
double BlockRelaxation(int iter) {
  ION *ion;
  RATE *r;
  LBLOCK *blk1, *blk2;
  int i, j, k, m;
  int p, q;
  double a, b, c, d;
  int nlevels;

  b = 1.0-iter_stablizer;
  c = iter_stablizer;
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    if (electron_density > 0.0) {
      for (m = 0; m < ion->ce_rates->dim; m++) {
	r = (RATE *) ArrayGet(ion->ce_rates, m);
	i = ion->iblock[r->i];
	p = ion->ilev[r->i];
	j = ion->iblock[r->f];
	q = ion->ilev[r->f];
	blk1 = (LBLOCK *) ArrayGet(blocks, i);
	blk2 = (LBLOCK *) ArrayGet(blocks, j);
	if (blk1->r[p] > 0.0) {
	  blk2->n[q] += blk1->r[p] * electron_density * r->dir;
	}
	if (r->inv > 0.0) {
	  if (blk2->r[q] > 0.0) {
	    blk1->n[p] += blk2->r[q] * electron_density * r->inv;
	  }
	}
      }
    }

    for (m = 0; m < ion->tr_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->tr_rates, m);
      i = ion->iblock[r->i];
      p = ion->ilev[r->i];
      j = ion->iblock[r->f];
      q = ion->ilev[r->f];    
      blk1 = (LBLOCK *) ArrayGet(blocks, i);
      blk2 = (LBLOCK *) ArrayGet(blocks, j);
      if (blk1->r[p] > 0.0) {
	blk2->n[q] += blk1->r[p] * r->dir;
      }
      if (r->inv > 0.0 && photon_density > 0.0) {
	a = photon_density * r->inv;
	if (blk1->r[p] > 0.0) {
	  blk2->n[q] += blk1->r[p]*a*(ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
	}
	if (blk2->r[q] > 0.0) {
	  blk1->n[p] += blk2->r[q] * a;
	}
      }
    }

    for (m = 0; m < ion->rr_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->rr_rates, m);
      i = ion->iblock[r->i];
      p = ion->ilev[r->i];
      j = ion->iblock[r->f];
      q = ion->ilev[r->f];    
      blk1 = (LBLOCK *) ArrayGet(blocks, i);
      blk2 = (LBLOCK *) ArrayGet(blocks, j);
      if (electron_density > 0.0) {
	if (blk1->r[p] > 0.0) {
	  blk2->n[q] += blk1->r[p] * r->dir;
	}
      } 
      if (r->inv > 0.0 && photon_density > 0.0) {
	if (blk2->r[q] > 0.0) {
	  blk1->n[p] += blk2->r[q] * photon_density * r->inv;
	}
      }
    }

    for (m = 0; m < ion->ai_rates->dim; m++) {
      r = (RATE *) ArrayGet(ion->ai_rates, m);
      i = ion->iblock[r->i];
      p = ion->ilev[r->i];
      j = ion->iblock[r->f];
      q = ion->ilev[r->f];    
      blk1 = (LBLOCK *) ArrayGet(blocks, i);
      blk2 = (LBLOCK *) ArrayGet(blocks, j);
      if (blk1->r[p] > 0.0) {
	blk2->n[q] += blk1->r[p] * r->dir;
      }
      if (r->inv > 0.0 && electron_density > 0.0) {
	if (blk2->r[q] > 0.0) {
	  blk1->n[p] += blk2->r[q] * electron_density * r->inv;
	}
      }
    }

    if (electron_density > 0.0) {
      for (m = 0; m < ion->ci_rates->dim; m++) {
	r = (RATE *) ArrayGet(ion->ci_rates, m);
	i = ion->iblock[r->i];
	p = ion->ilev[r->i];
	j = ion->iblock[r->f];
	q = ion->ilev[r->f];    
	blk1 = (LBLOCK *) ArrayGet(blocks, i);
	blk2 = (LBLOCK *) ArrayGet(blocks, j);
	if (blk1->r[p] > 0.0) {
	  blk2->n[q] += blk1->r[p] * electron_density * r->dir;
	}
	if (r->inv) {
	  if (blk2->r[q] > 0.0) {
	    blk1->n[p] += blk2->r[q] * electron_density * r->inv;
	  }
	}
      }
    }
  }

  nlevels = 0;
  d = 0.0;
  for (k = 0; k < blocks->dim; k++) {
    blk1 = (LBLOCK *) ArrayGet(blocks, k);
    a = 0.0;
    for (m = 0; m < blk1->nlevels; m++) {
      if (1.0+blk1->total_rate[m] != 1.0) {
	blk1->n[m] /= blk1->total_rate[m];
	a += blk1->n[m];
      } else {
	blk1->n[m] = 0.0;
      }
    }

    if (a == 0.0 || blk1->nb == 0.0) {
      for (m = 0; m < blk1->nlevels; m++) {
	blk1->n[m] = 0.0;
	blk1->r[m] = 0.0;
	blk1->n0[m] = 0.0;
      }
      blk1->r[0] = 1.0;
      continue;
    }

    a = blk1->nb/a;
    for (m = 0; m < blk1->nlevels; m++) {
      blk1->n[m] *= a;
      if (blk1->n[m] > 0.0) {
	d += fabs(1.0 - blk1->n0[m]/blk1->n[m]);
      }
      if (iter >= 3) {
	blk1->n[m] = b*blk1->n0[m] + c*blk1->n[m];
      }
      blk1->r[m] = blk1->n[m]/blk1->nb;
      blk1->n0[m] = blk1->n[m];      
    }
    nlevels += blk1->nlevels;
  }    

  d /= nlevels;

  return d;
}

int LevelPopulation(void) {
  int i;
  double d;

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

int SpecTable(char *fn, double strength_threshold) {
  SP_RECORD r, *t;
  ARRAY ri;
  SP_HEADER sp_hdr;
  F_HEADER fhdr;
  ION *ion;
  RATE *rt;
  LBLOCK *blk;
  DISTRIBUTION *edist, *pdist;
  int k, m;
  FILE *f;
  double e, a;
  int i, p, q;
  int iedist, ipdist;
  double smax, s;

  edist = GetEleDist(&iedist);
  pdist = GetPhoDist(&ipdist);
  fhdr.type = DB_SP;
  fhdr.atom = ion0.atom;
  strcpy(fhdr.symbol, ion0.symbol);
  sp_hdr.eden = electron_density;
  sp_hdr.iedist = iedist;
  sp_hdr.np_edist = edist->nparams;
  sp_hdr.p_edist = edist->params;
  sp_hdr.pden = photon_density;
  sp_hdr.ipdist = ipdist;
  sp_hdr.np_pdist = pdist->nparams;
  sp_hdr.p_pdist = pdist->params;

  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    sp_hdr.nele = ion->nele;
    sp_hdr.type = 0;
    f = InitFile(fn, &fhdr, &sp_hdr);
    ArrayInit(&ri, sizeof(SP_RECORD), 512);
    for (m = 0; m < ion->nlevels; m++) {
      i = ion->iblock[m];
      p = ion->ilev[m];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      if (blk->n[p] > 0.0) {
	r.upper = m;
	r.lower = i;
	r.energy = ion->energy[m];
	r.strength = blk->n[p];
	if (blk->iion != k) {
	  ArrayAppend(&ri, &r);
	} else {
	  WriteSPRecord(f, &r);
	}
      }
    }
    CloseFile(f, &fhdr);

    if (ri.dim > 0) {
      sp_hdr.nele = ion->nele-1;
      f = InitFile(fn, &fhdr, &sp_hdr);
      for (m = 0; m < ri.dim; m++) {
	t = (SP_RECORD *) ArrayGet(&ri, m);
	WriteSPRecord(f, t);
      }
      CloseFile(f, &fhdr);
    }
    ArrayFree(&ri, NULL);

    sp_hdr.nele = ion->nele;
    sp_hdr.type = 1;
    f = InitFile(fn, &fhdr, &sp_hdr);    
    smax = 0.0;
    ArrayInit(&ri, sizeof(SP_RECORD), 512);
    for (m = 0; m < ion->tr_rates->dim; m++) {
      rt = (RATE *) ArrayGet(ion->tr_rates, m);
      if (k == 0 && 
	  ion0.nionized > 0 &&
	  (p = IonizedIndex(rt->i, 1)) >= 0 &&
	  (q = IonizedIndex(rt->f, 1)) >= 0) {
	e = ion0.energy[p] - ion0.energy[q];
      } else {
	e = ion->energy[rt->i] - ion->energy[rt->f];
      }
      i = ion->iblock[rt->i];
      p = ion->ilev[rt->i];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      if (blk->n[p] > 0.0) {
	r.lower = rt->f;
	r.upper = rt->i;
	r.energy = e;
	r.strength = blk->n[p] * rt->dir;
	if (rt->inv > 0.0 && photon_density > 0.0) {
	  a = photon_density * rt->inv;
	  a *= (ion->j[rt->f]+1.0)/(ion->j[rt->i]+1.0);
	  r.strength += blk->n[p] * a;
	}
	s = r.strength*e;
	if (s < strength_threshold*smax) continue;
	if (s > smax) smax = s;
	if (blk->iion != k) {
	  ArrayAppend(&ri, &r);
	} else {
	  WriteSPRecord(f, &r);
	}
      }
    }
    CloseFile(f, &fhdr);

    if (ri.dim > 0) {
      sp_hdr.nele = ion->nele-1;
      f = InitFile(fn, &fhdr, &sp_hdr);
      for (m = 0; m < ri.dim; m++) {
	t = (SP_RECORD *) ArrayGet(&ri, m);
	WriteSPRecord(f, t);
      }
      CloseFile(f, &fhdr);
    }
    ArrayFree(&ri, NULL);
    
    if (electron_density <= 0) continue;
    if (ion->rr_rates->dim == 0) continue;
    sp_hdr.type = ion->nele;
    sp_hdr.type = 2;
    f = InitFile(fn, &fhdr, &sp_hdr);  
    smax = 0.0;
    for (m = 0; m < ion->rr_rates->dim; m++) {
      rt = (RATE *) ArrayGet(ion->rr_rates, m);
      e = ion->energy[rt->i] - ion->energy[rt->f];
      i = ion->iblock[rt->i];
      p = ion->ilev[rt->i];
      blk = (LBLOCK *) ArrayGet(blocks, i);
      if (blk->n[p] > 0.0) {
	r.lower = rt->f;
	r.upper = rt->i;
	r.energy = e;
	r.strength = blk->n[p] * rt->dir * electron_density;
	s = r.strength*e;
	if (s < strength_threshold*smax) continue;
	if (s > smax) smax = s;
	WriteSPRecord(f, &r);
      }
    }
    CloseFile(f, &fhdr);
  }
  
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

int PlotSpec(char *ifn, char *ofn, int type, 
	     double emin, double emax, double de, double smin) {
  F_HEADER fh;
  SP_HEADER h;
  SP_RECORD r;
  DISTRIBUTION *dist;
  FILE *f1, *f2;
  int n, nb, i;
  double e;
  int m, k, nsp;
  double *sp, *xsp, *kernel;
  double de10, de01;
  double a, sig, factor;
  double *lines;
  double smax;
  int swp;

  if (type == 0 || type > 2) {
    printf("Type must be > 0 and <= 2\n");
    return -1;
  }
  
  f1 = fopen(ifn, "r");
  if (f1 == NULL) {
    printf("ERROR: File %s does not exist\n", ifn);
    return -1;
  }
  f2 = fopen(ofn, "w");
  if (f2 == NULL) {
    printf("ERROR: File %s does not exist\n", ofn);
    return -1;
  }

  dist = GetEleDist(&i);
  
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
  sp[0] = 0.0;
  xsp[0] = emin;
  for (i = 1; i < nsp; i++) {
    sp[i] = 0.0;
    xsp[i] = xsp[i-1] + de01;
  }
  
  n = fread(&fh, sizeof(F_HEADER), 1, f1);
  if (CheckEndian(&fh) != CheckEndian(NULL)) {
    swp = 1;
    SwapEndianFHeader(&fh);
  } else {
    swp = 0;
  }
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = fread(&h, sizeof(SP_HEADER), 1, f1);
    if (swp) SwapEndianSPHeader(&h);
    m = sizeof(double)*(h.np_edist + h.np_pdist);
    fseek(f1, m, SEEK_CUR);
    if (h.ntransitions == 0) continue;
    if (h.type == type) {    
      if (type == 1) {	
	m = 2*h.ntransitions;
	lines = (double *) malloc(sizeof(double)*m);  
	k = 0;
	smax = 0.0;
	for (i = 0; i < h.ntransitions; i++) {
	  n = fread(&r, sizeof(SP_RECORD), 1, f1);
	  if (swp) SwapEndianSPRecord(&r);
	  e = r.energy;
	  a = r.strength * e;
	  if (a < smax*smin) continue;
	  if (a > smax) smax = a;
	  e *= HARTREE_EV;
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
	fflush(stdout);
	free(lines);
	for (i = 0; i < nsp; i++) xsp[i] = 0.0;
	for (i = 0; i < nsp; i++) {
	  if (sp[i] > 0.0) {
	    for (m = i-64, k = 0; k < 128; k++, m++) {
	      if (m > 0 && m < nsp) xsp[m] += sp[i]*kernel[k];
	    }
	  }
	}
	e = emin;
	for (i = 0; i < nsp; i++) {
	  fprintf(f2, "%15.8E\t%15.8E\n", e, xsp[i]);
	  sp[i] = 0.0;
	  xsp[i] = e;
	  e += de01;
	}
	fprintf(f2, "\n\n");
      } else if (type == 2) {
	printf("plotting RR continuum not implemented yet\n");
      }
    } else {
      fseek(f1, h.length, SEEK_CUR);
    }
  }

  free(xsp);
  free(sp);
  free(kernel);

  fclose(f1);
  fclose(f2);

  return 0;
}

int SetCERates(int inv) {
  int nb, i, j, t;
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
  float cs[MAXNUSR];
  double data[1+(1+MAXNUSR)*4];
  double *y, *x, *x2, *logx;
  double eusr[MAXNUSR];
  int swp, endian;
  
  endian = CheckEndian(NULL);
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  y = data + 1;
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->ce_rates, NULL);
    f = fopen(ion->dbfiles[DB_CE-1], "r");
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_CE-1]);
      continue;
    }
    n = fread(&fh, sizeof(F_HEADER), 1, f);
    if (CheckEndian(&fh) != endian) {
      swp = 1;
      SwapEndianFHeader(&fh);
    } else {
      swp = 0;
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = fread(&h, sizeof(CE_HEADER), 1, f);
      if (swp) SwapEndianCEHeader(&h);
      m = h.n_tegrid + h.n_egrid;
      fseek(f, sizeof(double)*m, SEEK_CUR);
      m = h.n_usr;
      n = fread(eusr, sizeof(double), m, f);
      if (swp) {
	for (i = 0; i < m; i++) {
	  SwapEndian((char *) &(h.usr_egrid[i]), sizeof(double));
	}
      }
      if (h.nele == ion->nele-1) {
	if (k > 0 || ion0.nionized > 0) {
	  fseek(f, h.length, SEEK_CUR);
	  continue;
	}
      }
      m1 = m + 1;
      x = y + m1;
      x2 = x + m1;
      logx = x2 + m1;
      x[0] = 0.0;
      x2[0] = 0.0;
      for (i = 0; i < h.ntransitions; i++) {
	n = fread(&r, sizeof(CE_RECORD), 1, f);
	if (swp) SwapEndianCERecord(&r);
	if (h.nparams > 0) {
	  fseek(f, sizeof(float)*h.nparams, SEEK_CUR);
	}
	rt.i = r.lower;
	rt.f = r.upper;
	j1 = ion->j[r.lower];
	j2 = ion->j[r.upper];
	e = ion->energy[r.upper] - ion->energy[r.lower];
	data[0] = r.bethe;	
	n = fread(cs, sizeof(float), m, f);
	if (swp) {
	  for (j = 0; j < m; j++) {
	    SwapEndian((char *) &(cs[j]), sizeof(float));
	  }
	}
	if (r.bethe < 0.0) {
	  y[0] = 0.0;
	  for (j = 0; j < m; j++) {
	    t = m-j;
	    x[t] = e/(e + eusr[j]);
	    x2[t] = x[t]*x[t];
	    y[t] = cs[j];
	  }
	} else {
	  if (r.bethe > 0.0) {
	    for (j = 0; j < m; j++) {
	      t = m-j;
	      x[t] = e/(e + eusr[j]);
	      logx[j] = log(x[t]);
	      y[t] = cs[j] + r.bethe*logx[j];
	    }
	  } else {
	    for (j = 0; j < m; j++) {
	      t = m-j;
	      x[t] = e/(e + eusr[j]);
	      y[t] = cs[j];
	    }
	  }
	  n = 3;
	  m1 = 1;
	  uvip3p_(&n, &m, &(x[1]), &(y[1]), &m1, &(x[0]), &(y[0]));
	  if (r.bethe + 1.0 == 1.0) {
	    if (y[0] < 0.0) y[0] = 0.0;
	  }
	}
	CERate(&(rt.dir), &(rt.inv), inv, j1, j2, e, m,
	       data, rt.i, rt.f);
	ArrayAppend(ion->ce_rates, &rt);
      }
    }
    fclose(f);
    
    if (k == 0 && ion0.nionized > 0) {
      f = fopen(ion0.dbfiles[DB_CE-1], "r");
      if (f == NULL) {
	printf("File %s does not exist, skipping.\n", ion0.dbfiles[DB_CE-1]);
	continue;
      }
      n = fread(&fh, sizeof(F_HEADER), 1, f);
      if (CheckEndian(&fh) != endian) {
	swp = 1;
	SwapEndianFHeader(&fh);
      } else {
	swp = 0;
      }
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = fread(&h, sizeof(CE_HEADER), 1, f);
	if (swp) SwapEndianCEHeader(&h);
	m = h.n_tegrid + h.n_egrid;
	fseek(f, sizeof(double)*m, SEEK_CUR);
	m = h.n_usr;   
	n = fread(eusr, sizeof(double), m, f);
	if (swp) {
	  for (i = 0; i < m; i++) {
	    SwapEndian((char *) &(eusr[i]), sizeof(double));
	  }
	}
	if (h.nele != ion0.nele) {
	  fseek(f, h.length, SEEK_CUR);
	  continue;
	}
	m1 = m + 1;
	x = y + m1;
	x2 = x + m1;
	logx = x2 + m1;
	x[0] = 0.0;
	x2[0] = 0.0;
	for (i = 0; i < h.ntransitions; i++) {
	  n = fread(&r, sizeof(CE_RECORD), 1, f);
	  if (swp) SwapEndianCERecord(&r);
	  p = IonizedIndex(r.lower, 0);
	  if (p < 0) {
	    fseek(f, sizeof(float)*(h.nparams+h.n_usr), SEEK_CUR);
	    continue;
	  }
	  q = IonizedIndex(r.upper, 0);
	  if (q < 0) {
	    fseek(f, sizeof(float)*(h.nparams+h.n_usr), SEEK_CUR);
	    continue;
	  }
	  if (h.nparams > 0) {
	    fseek(f, sizeof(float)*h.nparams, SEEK_CUR);
	  }
	  rt.i = ion0.ionized_map[1][p];
	  rt.f = ion0.ionized_map[1][q];
	  j1 = ion->j[rt.i];
	  j2 = ion->j[rt.f];
	  e = ion0.energy[q] - ion0.energy[p];
	  data[0] = r.bethe;
	  n = fread(cs, sizeof(float), m, f);
	  if (swp) {
	    for (j = 0; j < m; j++) {
	      SwapEndian((char *) &(cs[j]), sizeof(float));
	    }
	  }
	  if (r.bethe < 0.0) {
	    y[0] = 0.0;
	    for (j = 0; j < m; j++) {
	      t = m-j;
	      x[t] = e/(e + eusr[j]);
	      x2[t] = x[t]*x[t];
	      y[t] = cs[j];
	    }
	  } else {
	    if (r.bethe > 0.0) {
	      for (j = 0; j < m; j++) {
		t = m-j;
		x[t] = e/(e + eusr[j]);
		logx[j] = log(x[t]);
		y[t] = cs[j] + r.bethe*logx[j];
	      }
	    } else {
	      for (j = 0; j < m; j++) {
		t = m-j;
		x[t] = e/(e + eusr[j]);
		y[t] = cs[j];
	      }
	    }
	    n = 3;
	    m1 = 1;
	    uvip3p_(&n, &m, &(x[1]), &(y[1]), &m1, &(x[0]), &(y[0]));
	    if (r.bethe + 1.0 == 1.0) {
	      if (y[0] < 0.0) y[0] = 0.0;
	    }
	  }
	  CERate(&(rt.dir), &(rt.inv), inv, j1, j2, e, m,
		 data, rt.i, rt.f);
	  ArrayAppend(ion->ce_rates, &rt);
	}
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
  int p, q;
  ION *ion;
  RATE rt;
  F_HEADER fh;
  TR_HEADER h;
  TR_RECORD r;
  double e;
  FILE *f;  
  int swp, endian;

  endian = CheckEndian(NULL);
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->tr_rates, NULL);
    f = fopen(ion->dbfiles[DB_TR-1], "r");
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_TR-1]);
      continue;
    }
    n = fread(&fh, sizeof(F_HEADER), 1, f);
    if (CheckEndian(&fh) != endian) {
      swp = 1;
      SwapEndianFHeader(&fh);
    } else {
      swp = 0;
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = fread(&h, sizeof(TR_HEADER), 1, f);
      if (swp) SwapEndianTRHeader(&h);
      if (h.nele == ion->nele-1) {
	if (k > 0 || ion0.nionized > 0) {
	  fseek(f, h.length, SEEK_CUR);
	  continue;
	}
      }
      for (i = 0; i < h.ntransitions; i++) {
	n = fread(&r, sizeof(TR_RECORD), 1, f);
	if (swp) SwapEndianTRRecord(&r);
	rt.i = r.upper;
	rt.f = r.lower;
	j1 = ion->j[r.upper];
	j2 = ion->j[r.lower];
	e = ion->energy[r.upper] - ion->energy[r.lower];
	TRRate(&(rt.dir), &(rt.inv), inv, j1, j2, e, r.strength);
	ArrayAppend(ion->tr_rates, &rt);
      }
    }
    fclose(f);

    if (k == 0 && ion0.nionized > 0) {
      f = fopen(ion0.dbfiles[DB_TR-1], "r");
      if (f == NULL) {
	printf("File %s does not exist, skipping.\n", ion0.dbfiles[DB_CE-1]);
	continue;
      }
      n = fread(&fh, sizeof(F_HEADER), 1, f);
      if (CheckEndian(&fh) != endian) {
	swp = 1;
	SwapEndianFHeader(&fh);
      } else {
	swp = 0;
      }
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = fread(&h, sizeof(TR_HEADER), 1, f);
	if (swp) SwapEndianTRHeader(&h);
	if (h.nele != ion0.nele) {
	  fseek(f, h.length, SEEK_CUR);
	  continue;
	}  
	for (i = 0; i < h.ntransitions; i++) {
	  n = fread(&r, sizeof(TR_RECORD), 1, f);
	  if (swp) SwapEndianTRRecord(&r);
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
	  TRRate(&(rt.dir), &(rt.inv), inv, j1, j2, e, r.strength);
	  ArrayAppend(ion->tr_rates, &rt);
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
  int nshells;
  ION *ion;
  RATE rt;
  F_HEADER fh;
  CI_HEADER h;
  CI_RECORD r;
  double e;
  FILE *f;  
  float *params;
  int swp, endian;

  endian = CheckEndian(NULL);

  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->ci_rates, NULL);
    f = fopen(ion->dbfiles[DB_CI-1], "r");
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_CI-1]);
      continue;
    }
    n = fread(&fh, sizeof(F_HEADER), 1, f);
    if (CheckEndian(&fh) != endian) {
      swp = 1;
      SwapEndianFHeader(&fh);
    } else {
      swp = 0;
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = fread(&h, sizeof(CI_HEADER), 1, f);
      if (swp) SwapEndianCIHeader(&h);
      m = h.n_tegrid + h.n_egrid + h.n_usr;
      fseek(f, sizeof(double)*m, SEEK_CUR);
      m = h.nparams;
      nshells = 1;
      params = (float *) malloc(sizeof(float)*m*nshells);
      for (i = 0; i < h.ntransitions; i++) {
	n = fread(&r, sizeof(CI_RECORD), 1, f);
	if (swp) SwapEndianCIRecord(&r);
	if (r.nshells > nshells) {
	  nshells = r.nshells;
	  params = (float *) realloc(params, sizeof(float)*m*nshells);
	}
	n = fread(params, sizeof(float), m*r.nshells, f);
	if (swp) {
	  for (t = 0; t < m*r.nshells; t++) {
	    SwapEndian((char *) &(params[t]), sizeof(float));
	  }
	}
	fseek(f, sizeof(float)*h.n_usr, SEEK_CUR);
	rt.i = r.b;
	rt.f = r.f;
	j1 = ion->j[r.b];
	j2 = ion->j[r.f];
	e = ion->energy[r.f] - ion->energy[r.b];
	CIRate(&(rt.dir), &(rt.inv), inv, j1, j2, e, m, nshells, params,
	       rt.i, rt.f);
	ArrayAppend(ion->ci_rates, &rt);
      }
      free(params);
    }
    fclose(f);
  }
  return 0;
}

int SetRRRates(int inv) { 
  int nb, i, t;
  int n, m, k;
  int j1, j2;
  int nshells;
  ION *ion;
  RATE rt;
  F_HEADER fh;
  RR_HEADER h;
  RR_RECORD r;
  double e;
  FILE *f;  
  float *params;
  int endian, swp;

  endian = CheckEndian(NULL);
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->rr_rates, NULL);
    f = fopen(ion->dbfiles[DB_RR-1], "r");
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_RR-1]);
      continue;
    }
    n = fread(&fh, sizeof(F_HEADER), 1, f);
    if (CheckEndian(&fh) != endian) {
      swp = 1;
      SwapEndianFHeader(&fh);
    } else {
      swp = 0;
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = fread(&h, sizeof(RR_HEADER), 1, f);
      if (swp) SwapEndianRRHeader(&h);
      m = h.n_tegrid + h.n_egrid + h.n_usr;
      fseek(f, sizeof(double)*m, SEEK_CUR);
      m = h.nparams;
      nshells = 1;
      params = (float *) malloc(sizeof(float)*m*nshells);
      for (i = 0; i < h.ntransitions; i++) {
	n = fread(&r, sizeof(RR_RECORD), 1, f);
	if (swp) SwapEndianRRRecord(&r);
	if (r.nshells > nshells) {
	  nshells = r.nshells;
	  params = (float *) realloc(params, sizeof(float)*m*nshells);
	}
	n = fread(params, sizeof(float), m*r.nshells, f);
	if (swp) {
	  for (t = 0; t < m*r.nshells; t++) {
	    SwapEndian((char *) &(params[t]), sizeof(float));
	  }
	}
	  
	fseek(f, sizeof(float)*h.n_usr, SEEK_CUR);
	rt.i = r.f;
	rt.f = r.b;
	j1 = ion->j[r.f];
	j2 = ion->j[r.b];
	e = ion->energy[r.f] - ion->energy[r.b];
	if (e < 0.0) {
	  printf("%d %d %10.3E %10.3E\n", 
		 r.f, r.b, ion->energy[r.f],ion->energy[r.b]);
	  exit(1);
	}
	RRRate(&(rt.dir), &(rt.inv), inv, j1, j2, e, m, nshells, params,
	       rt.i, rt.f);
	ArrayAppend(ion->rr_rates, &rt);
      }
      free(params);
    }
    fclose(f);
  }
  return 0;
}

int SetAIRates(int inv) {
  int nb, i;
  int n, k;
  int j1, j2;
  ION *ion;
  RATE rt;
  F_HEADER fh;
  AI_HEADER h;
  AI_RECORD r;
  double e;
  FILE *f;  
  int swp, endian;

  endian = CheckEndian(NULL);

  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->ai_rates, NULL);
    f = fopen(ion->dbfiles[DB_AI-1], "r");
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_AI-1]);
      continue;
    }
    n = fread(&fh, sizeof(F_HEADER), 1, f);
    if (CheckEndian(&fh) != endian) {
      swp = 1;
      SwapEndianFHeader(&fh);
    } else {
      swp = 0;
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = fread(&h, sizeof(AI_HEADER), 1, f);
      if (swp) SwapEndianAIHeader(&h);
      fseek(f, sizeof(double)*h.n_egrid, SEEK_CUR);
      for (i = 0; i < h.ntransitions; i++) {
	n = fread(&r, sizeof(AI_RECORD), 1, f);
	if (swp) SwapEndianAIRecord(&r);
	rt.i = r.b;
	rt.f = r.f;
	j1 = ion->j[r.b];
	j2 = ion->j[r.f];
	e = ion->energy[r.b] - ion->energy[r.f];
	AIRate(&(rt.dir), &(rt.inv), inv, j1, j2, e, r.rate);
	ArrayAppend(ion->ai_rates, &rt);
      }
    }
    fclose(f);
  }
  return 0;
}
