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

#include "crm.h"
#include "grid.h"
#include "cf77.h"
#include "mpiutil.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#define NRTB 8192

static IONIZED ion0;
static ARRAY *ions;
static ARRAY *blocks;
static double *bmatrix = NULL;

static int n_single_blocks = 64;
static int _sbnmax1 = 0;
static int _sbnmax2 = 0;

static int rec_cascade = 0;
static double cas_accuracy = EPS4;
static int max_iter = 256;
static double iter_accuracy = EPS4;
static double iter_stabilizer = 0.8;
static int  _silent = 0;
static int _krc = -1;

/* electron density in 10^10 cm-3 */
static double electron_density = 0.0;
/* photon energy density in erg cm-3 */
/* if the distribution is blackbody, then the normalization constant
 * is the dilution factor, (r0/r)^2 */
static double photon_density = 0.0;
static double cxt_density = 0.0;
static int ai_extra_nmax = 400;
static int do_extrapolate = 100;
static int inner_auger = 0;
static double ai_emin = 0.0;
static int norm_mode = 1;
static int sw_mode = 0;
static int _itol_nmax = 0;
static int _ce_nmax = 0;
static double _ce_fbr = 0.175;
static double _ce_xbr = 0.875;
static int _sp_trm = 1;
static int _rates_block = RATES_BLOCK;
static int _lblock_block = LBLOCK_BLOCK;

static double _ce_data[2+(1+MAXNUSR)*2];
static double _rr_data[1+MAXNUSR*4];

static double _starkrw = 1.0;
static double _starkqc = 1.0;
static double _starkbt = 0.0;
static int _starknp = 0;
static double *_starkzp = NULL;
static double *_starkmp = NULL;
static double *_starkwp = NULL;
static double _starkaix = 1.0;
static double _starksmx = 3.0;
static double _starkzix = -1.0;
static double _starkvg = 0.0;
static double _epstau = 0.05;
static double _reemit = 1.0;
static INTERPSP _interpsp;
static double _mfd0[10] = {0.        ,  0.97792763,  1.26508455,
			   1.22263597, -0.25197702,
			   0.15081325, -0.1970343 ,
			   0.08390483,  1.18112936,  1.19006095
};
static double _mfd1[10] = {0.56387871,  0.93997153,  1.24750821,
			   1.22705255, -0.1960854 ,
			   0.14900134, -0.08442087,  0.06137596,
			   1.13166003,  1.00234922
};
static double _mfd2[10] = {1.04218924,  0.93440751,  1.23181356,
			   1.62415867, -0.08177814,
			   0.16991271,  0.01408619, -0.02417788,
			   0.14943566,  0.26701825
};
static double _mfd3[10] = {1.11558761,  0.23828618,  1.28638484,
			   5.        , -0.03231508,
			   0.63684014,  0.09547696, -0.12634639,
			   0.02654152,  0.02242574
};
static double _mfd4[10] = {0.75575757,  0.96542404,  1.22862064,
			   0.82448354, -0.14021655,
			   0.03221276,  0.01401347,  0.24654512,
			   0.95572507,  0.1262922
};

#pragma omp threadprivate(_ce_data, _rr_data)

int NormalizeMode(int i) {
  norm_mode = i;
  return 0;
}

int SetInnerAuger(int i) {
  inner_auger = i;
  return 0;
}

int SetExtrapolate(int e) {
  do_extrapolate = e;
  return 0;
}

int SetEMinAI(double e) {
  ai_emin = e;  
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

int SetCxtDensity(double cxt) {
  if (cxt >= 0.0) cxt_density = cxt;
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

  for (i = 0; i < NDB1; i++) ion0.dbfiles[i] = NULL;
  ion0.nionized = 0;
  ion0.energy = NULL;
  ion0.atom = 0;
  ion0.atr = ion0.ace = ion0.aci = ion0.arr = ion0.aai = ion0.acx = -1;

  ions = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ions, sizeof(ION), ION_BLOCK);
  blocks = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(blocks, sizeof(LBLOCK), _lblock_block);
  bmatrix = NULL;

  _interpsp.nd = 0;
  _interpsp.nt = 0;
  _interpsp.r = NULL;
  _interpsp.xd = NULL;
  _interpsp.xt = NULL;

  SetStarkZMP(1, NULL);
  InitDBase();
  InitRates();
  InitCoulomb();

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
  if (r->rates) {
    ArrayFreeLock(r->rates, NULL);
    free(r->rates);
    r->rates = NULL;
  }
}

static void InitIonData(void *p, int n) {
  ION *ion;
  int i, k;
  
  ion = (ION *) p;
  for (k = 0; k < n; k++,ion++) {
    ion->nele = 0;
    ion->nlevels = 0;
    ion->KLN_min = 0;
    ion->KLN_max = -1;
    ion->KLN_bmin = 0;
    ion->KLN_bmax = -1;
    ion->KLN_amin = 0;
    ion->KLN_amax = -1;
    ion->ace = -1;
    ion->atr = -1;
    ion->aci = -1;
    ion->arr = -1;
    ion->aai = -1;
    ion->acx = -1;
    ion->nlevels = 0;
    ion->iblock = NULL;
    ion->ilev = NULL;
    ion->j = NULL;
    ion->vnl = NULL;
    ion->ibase = NULL;
    ion->energy = NULL;
    for (i = 0; i < NDB1; i++) {
      ion->dbfiles[i] = NULL;
    }
    for (i = 0; i < 4; i++) {
      ion->icx[i] = NULL;
    }
    ion->ce_rates = NULL;
    ion->tr_rates = NULL;
    ion->tr_sdev = NULL;
    ion->tr2_rates = NULL;
    ion->ci_rates = NULL;
    ion->rr_rates = NULL;
    ion->ai_rates = NULL;
    ion->cx_rates = NULL;
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
    free(ion->sw);
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
    ion->ace = ion->atr = ion->aci = ion->arr = ion->aai = ion->acx = -1;
  }
  for (i = 0; i < NDB1; i++) {
    if (ion->dbfiles[i]) free(ion->dbfiles[i]);
    ion->dbfiles[i] = NULL;
  }
  ArrayFreeLock(ion->ce_rates, FreeBlkRateData);
  free(ion->ce_rates);
  ion->ce_rates = NULL;
  ArrayFreeLock(ion->tr_rates, FreeBlkRateData);
  free(ion->tr_rates);
  ion->tr_rates = NULL;
  ArrayFreeLock(ion->tr_sdev, FreeBlkRateData);
  free(ion->tr_sdev);
  ion->tr_sdev = NULL;
  free(ion->tr2_rates);
  ion->tr2_rates = NULL;
  ArrayFreeLock(ion->ci_rates, FreeBlkRateData);
  free(ion->ci_rates);
  ion->ci_rates = NULL;
  ArrayFreeLock(ion->rr_rates, FreeBlkRateData);
  free(ion->rr_rates);
  ion->rr_rates = NULL;
  ArrayFreeLock(ion->ai_rates, FreeBlkRateData);
  free(ion->ai_rates);
  ion->ai_rates = NULL;
  ArrayFreeLock(ion->cx_rates, FreeBlkRateData);
  free(ion->cx_rates);
  ion->cx_rates = NULL;
  ArrayFreeLock(ion->recombined, NULL);
  free(ion->recombined);
  ion->recombined = NULL;
  for (i = 0; i < 4; i++) {
    if (ion->icx[i] != NULL) {
      free(ion->icx[i]);
      ion->icx[i] = NULL;
    }
  }
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
    free(blk->rc0);
    free(blk->rc1);
  }
  blk->nlevels = 0;
}

int ReinitCRM(int m) {
  ION *ion;
  int i, k;

  if (m < 0) return 0;

  if (m == 0) {
    ReinitDBase(0);
  } else {
    ReinitDBase(DB_SP);
    ReinitDBase(DB_RT);
    ReinitDBase(DB_DR);
  }
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
      ArrayFree(ion->cx_rates, FreeBlkRateData);
    }
    return 0;
  } else if (m == 2) {
    for (k = 0; k < ions->dim; k++) {
      ion = (ION *) ArrayGet(ions, k);
      ArrayFree(ion->ce_rates, FreeBlkRateData);
      ArrayFree(ion->ci_rates, FreeBlkRateData);
      ArrayFree(ion->rr_rates, FreeBlkRateData);
      ArrayFree(ion->ai_rates, FreeBlkRateData);
      ArrayFree(ion->cx_rates, FreeBlkRateData);
    }
    return 0;
  }

  for (i = 0; i < NDB1; i++) {
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
  if (bmatrix) {
    free(bmatrix);
  }
  bmatrix = NULL;
  
  return 0;
}

int AddIon(int nele, double n, char *pref) {
  ION ion;
  int i;
  int m;

#if USE_MPI == 2
  if (!MPIReady()) InitializeMPI(0, 1);
#endif
  ion.nlevels = 0;
  ion.ce_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInitID(ion.ce_rates, sizeof(BLK_RATE), RATES_BLOCK, "CE");  
  ion.tr_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInitID(ion.tr_rates, sizeof(BLK_RATE), RATES_BLOCK, "TR");
  ion.tr_sdev = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInitID(ion.tr_sdev, sizeof(BLK_RATE), RATES_BLOCK, "TRD");
  ion.tr2_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInitID(ion.tr2_rates, sizeof(BLK_RATE), RATES_BLOCK, "TR2");
  ion.rr_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInitID(ion.rr_rates, sizeof(BLK_RATE), RATES_BLOCK, "RR");
  ion.ci_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInitID(ion.ci_rates, sizeof(BLK_RATE), RATES_BLOCK, "CI");
  ion.ai_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInitID(ion.ai_rates, sizeof(BLK_RATE), RATES_BLOCK, "AI");
  ion.cx_rates = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInitID(ion.cx_rates, sizeof(BLK_RATE), RATES_BLOCK, "CX");

  ion.KLN_min = 0;
  ion.KLN_max = -1;
  ion.KLN_bmin = 0;
  ion.KLN_bmax = -1;
  ion.KLN_amin = 0;
  ion.KLN_amax = -1;
  ion.ace = ion.atr = ion.aci = ion.arr = ion.aai = ion.acx = -1;
  ion.KLN_ai = NULL;
  ion.KLN_nai = NULL;

  ion.recombined = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ion.recombined, sizeof(RECOMBINED), 16);
  
  ion.nele = nele;
  m = strlen(pref);
  m = m+4;
  for (i = 0; i < NDB1; i++) {
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
    case DB_RO:
      if (sw_mode >= 4) {
	sprintf(ion.dbfiles[i], "%s.sro", pref);
      } else {
	sprintf(ion.dbfiles[i], "%s.ro", pref);
      }
      break;
    case DB_CX:
      sprintf(ion.dbfiles[i], "%s.cx", pref);
      break;
    case NDB1:
      sprintf(ion.dbfiles[i], "%s.LS", pref);
      break;
    default:
      break;
    }
  }

  ion.n = n;
  for (i = 0; i < 4; i++) {
    ion.icx[i] = NULL;
  }
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
  ion->vni = (short *) realloc(ion->vni, sizeof(short)*nlev);
  ion->ibase = (short *) realloc(ion->ibase, sizeof(short)*nlev);
  ion->sw = (short *) realloc(ion->sw, sizeof(short)*nlev);
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
	blk.rc0 = (double *) malloc(sizeof(double)*nr);
	blk.rc1 = (double *) malloc(sizeof(double)*nr);
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
  
void ExtrapolateTR(ION *ion, int inv, int **irb) {
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
	    AddRate(ion, ion->tr_rates, &r0, 0, irb);
	    if (iuta) {
	      r0.dir = ion->energy[r0.i] - ion->energy[r0.f];
	      r0.inv = 0.0;
	      AddRate(ion, ion->tr_sdev, &r0, 0, irb);
	    }
	  } else {
	    for (s = 0; s < blk->rec->n_ext; s++) {
	      if (blk->rec->nrec[s] == n0) {
		r0.f = q + blk->rec->imin[s];
		if (r0.f > blk->rec->imax[s]) break;
		r0.dir = r->dir;
		r0.inv = r->inv;
		AddRate(ion, ion->tr_rates, &r0, 0, irb);
		if (iuta) {
		  r0.dir = ion->energy[r0.i] - ion->energy[r0.f];
		  r0.inv = 0.0;
		  AddRate(ion, ion->tr_sdev, &r0, 0, irb);
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

void ExtrapolateRR(ION *ion, int inv, int **irb) {
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
    if (rec->n_ext == rec->n) continue;
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
	  AddRate(ion, ion->rr_rates, &r0, 0, irb);
	}
	r->dir *= a;
      }
    }
  }
}

void ExtrapolateAI(ION *ion, int inv, int **irb) {
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
    if (rec->n_ext == rec->n) continue;
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
	  AddRate(ion, ion->ai_rates, &r0, 0, irb);
	}
	r->inv *= a;
	r->dir *= c;
      }
    }  
  }
}

void LoadSW(ION *ion) {
  FILE *f;
  int n, i, j, p, tj, ts, tk;
  double w;
  char buf[8192];

  if (sw_mode != 0 && sw_mode != 1) return;
  if (ion->dbfiles[NDB] == NULL) return;
  f = fopen(ion->dbfiles[NDB], "r");
  if (f == NULL) {
    //printf("cannot open file %s\n", ion->dbfiles[NDB]);
    return;
  }
  i = -1;
  while (1) {
    if (NULL == fgets(buf, 8192, f)) break;
    if (buf[0] == '#') {
      if (buf[1] == '#') continue;
      i = atoi(&buf[1]);
    } else {
      n = sscanf(buf, "%d %d %d %d %d %lg", &j, &p, &tj, &ts, &tk, &w);
      if (n != 6) continue;
      if (i >= 0) {
	if (i < ion->nlevels) {
	  if (sw_mode == 0) {
	    ion->sw[i] = ts;
	  } else {
	    j = ion->vnl[i]%100;
	    p = tk-j;
	    ion->sw[i] = abs(p)*100 + ts;
	    if (p < 0) ion->sw[i] = -ion->sw[i];
	  }
	  //printf("%d %d %d %d %d %d %d %g\n", i, j, p, tj, ts, tk, ion->sw[i], w);
	}
	i = -1;
      }
    }
  }
  fclose(f);
}

int SingleLevelBlock(int i, int nb, int nmx, int qmx) {
  if (n_single_blocks == 0) return 1;
  if (n_single_blocks < 0 && nb < -n_single_blocks) return 2;
  if (qmx == 1 && nmx <= _sbnmax1) return 3;
  if (qmx > 1 && nmx <= _sbnmax2) return 4;
  if (nb == 0) {
    if (i < n_single_blocks) return 5;
    if (i == n_single_blocks) return -1;
  }
  return 0;
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
  TFILE *f;
  int n, i, k, nb, nb0, nlevels;
  char *fn;
  int p, q = -1;
  int nionized, n0, imx, nmx, qmx, isb;
  int swp, sfh;
  int s3nk[10000];

  for (i = 0; i < 10000; i++) s3nk[i] = 0;
  ion0.n = ni;
  ion0.n0 = ni;
  if (ifn) {
    k = strlen(ifn);
    k += 4;
  } else {
    k = 0;
  }
  for (i = 0; i < NDB1; i++) {
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
      case NDB1:
	ion0.dbfiles[i] = (char *) malloc(k);
	sprintf(ion0.dbfiles[i], "%s.LS", ifn);
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
    f = OpenFileRO(fn, &fh, &swp);
    if (f == NULL) {
      printf("File %s does not exist\n", fn);
      return -1;
    }
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
      FSEEK(f, h.length, SEEK_CUR);
    }
    ion->nlevels = nlevels;
    ion->iblock = (LBLOCK **) malloc(sizeof(LBLOCK *)*nlevels);
    for (i = 0; i < nlevels; i++) ion->iblock[i] = NULL;
    ion->ilev = (int *) malloc(sizeof(int)*nlevels);
    ion->j = (int *) malloc(sizeof(int)*nlevels);
    ion->p = (short *) malloc(sizeof(short)*nlevels);
    ion->vnl = (short *) malloc(sizeof(short)*nlevels);
    ion->vni = (short *) malloc(sizeof(short)*nlevels);
    ion->ibase = (short *) malloc(sizeof(short)*nlevels);
    ion->sw = (short *) malloc(sizeof(short)*nlevels);
    ion->energy = (double *) malloc(sizeof(double)*nlevels);
    rionized = (EN_RECORD *) malloc(sizeof(EN_RECORD )*nionized);
    if (k == 0 && ifn) {
      ion0.nionized = nionized;
      ion0.ionized_map[0] = (int *) malloc(sizeof(int)*nionized);
      ion0.ionized_map[1] = (int *) malloc(sizeof(int)*nionized);
      ion0.energy = (double *) malloc(sizeof(double)*nionized);
    }
    
    FSEEK(f, sfh, SEEK_SET);
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
	FSEEK(f, h.length, SEEK_CUR);
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
	int nilev = h.nlevels;
	if (ifn) {
	  nilev = FindLevelBlock(h.nlevels, r0, &r1, ion->nele-1, ifn); 
	  if (nilev <= 0) {
	    printf("Ionized block %d of ion %d ", nb, ion->nele);
	    printf("does not match a block in file %s\n", ifn);
	    printf("nlevels = %d VS %d\n", h.nlevels, nilev);
	    //Abort(1);
	  }
	  printf("ionized match: %d %d %d %d %d\n", ion->nele, h.nele, nb, h.nlevels, nilev);	  
	}
	if (k > 0) {
	  for (i = 0; i < h.nlevels; i++) {
	    p = r0[i].ilev;
	    if (i >= nilev) {
	      ion->iblock[p] = NULL;
	      continue;
	    }
	    q = r1[i].ilev;
	    ion1->iblock[q]->ionized = 1;
	    ion->iblock[p] = ion1->iblock[q];
	    ion->ilev[p] = ion1->ilev[q];
	    ion->j[p] = JFromENRecord(&(r0[i]));
	    if (r0[i].p < 0) {
	      ion->vnl[p] = -r0[i].p;
	      ion->p[p] = 1;
	    } else {
	      ion->vnl[p] = r0[i].p;
	      ion->p[p] = 0;
	    }
	    ion->vni[p] = VNIFromSName(r0[i].sname);
	    ion->ibase[p] = -1;
	    ion->energy[p] = r0[i].energy;
	  }
	} else {
	  blk.ncomplex[0].n = 0;
	  blk.nlevels = 0;
	  blkp = NULL;
	  for (i = 0; i < h.nlevels; i++) {
	    imx = GetNComplex(ncomplex, r0[i].ncomplex)-1;
	    nmx = ncomplex[imx].n;
	    qmx = ncomplex[imx].nq;
	    isb = SingleLevelBlock(i, nb0, nmx, qmx);
	    if (isb != 0) {
	      nlevels = 0;
	      blk.ib = blocks->dim;
	      blk.iion = -1;
	      blk.irec = -1;
	      blk.ionized = 1;
	      blk.rec = NULL;
	      if (isb > 0) {
		blk.nlevels = 1;
	      } else {
		blk.nlevels = h.nlevels - n_single_blocks;
	      }      
	      blk.n = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.n0 = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.r = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.total_rate = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.rc0 = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.rc1 = (double *) malloc(sizeof(double)*blk.nlevels);
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
		  blkp->rc0 = (double *)ReallocNew(blkp->rc0,
						   sizeof(double)*nlevels);
		  blkp->rc1 = (double *)ReallocNew(blkp->rc1,
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
	      blk.rc0 = (double *) malloc(sizeof(double)*blk.nlevels);
	      blk.rc1 = (double *) malloc(sizeof(double)*blk.nlevels);
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
	      ion->p[p] = 1;
	    } else {
	      ion->vnl[p] = r0[i].p;
	      ion->p[p] = 0;
	    }
	    ion->vni[p] = VNIFromSName(r0[i].sname);
	    ion->ibase[p] = -1;
	    ion->energy[p] = r0[i].energy;
	    if (ifn) {
	      if (i < nilev) {
		if (r1[i].ilev > ion0.imax[0]) ion0.imax[0] = r1[i].ilev;
		if (r0[i].ilev > ion0.imax[1]) ion0.imax[1] = r0[i].ilev;
		if (r1[i].ilev < ion0.imin[0]) ion0.imin[0] = r1[i].ilev;
		if (r0[i].ilev < ion0.imin[1]) ion0.imin[1] = r0[i].ilev;
		ion0.ionized_map[0][n0] = r1[i].ilev;
		ion0.ionized_map[1][n0] = r0[i].ilev;
		ion0.energy[n0] = r1[i].energy;
		n0++;
	      } else {
		ion->iblock[p] = NULL;	      
	      }
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
	FSEEK(f, h.length, SEEK_CUR);
	continue;	
      } else {
	printf("ERROR: Ion charge state does not match %d %d %d %d\n",
	       k, nb, h.nele, ion->nele);
	FSEEK(f, h.length, SEEK_CUR);
	continue;
      }
    }
    if (k == 0) {
      ion0.nionized = n0;
    }
    FSEEK(f, sfh, SEEK_SET);
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadENHeader(f, &h, swp);
      if (h.nele != ion->nele) {
	FSEEK(f, h.length, SEEK_CUR);
	continue;
      }	
      blk.ncomplex[0].n = 0;
      blkp = NULL;
      nlevels = 0;
      for (i = 0; i < h.nlevels; i++) {
	n = ReadENRecord(f, &r, swp);
	imx = GetNComplex(ncomplex, r.ncomplex)-1;
	nmx = ncomplex[imx].n;
	qmx = ncomplex[imx].nq;
	isb = SingleLevelBlock(i, nb, nmx, qmx);
	if (isb != 0) {
	  nlevels = 0;
	  blk.ib = blocks->dim;
	  blk.iion = k;
	  blk.irec = -1;
	  blk.ionized = 0;
	  blk.rec = NULL;
	  if (isb > 0) {
	    blk.nlevels = 1;
	  } else {
	    blk.nlevels = h.nlevels - n_single_blocks;
	  }
	  blk.n = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.n0 = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.r = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.total_rate = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.rc0 = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.rc1 = (double *) malloc(sizeof(double)*blk.nlevels);
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
	      blkp->rc0 = (double *) ReallocNew(blkp->rc0,
						sizeof(double)*nlevels);
	      blkp->rc1 = (double *) ReallocNew(blkp->rc1,
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
	  blk.rc0 = (double *) malloc(sizeof(double)*blk.nlevels);
	  blk.rc1 = (double *) malloc(sizeof(double)*blk.nlevels);
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
	  ion->p[p] = 1;
	} else {
	  ion->vnl[p] = r.p;
	  ion->p[p] = 0;
	}
	ion->vni[p] = VNIFromSName(r.sname);
	ion->ibase[p] = IBaseFromENRecord(&r);
	if (ion->ibase[p] < 0 && ion->iground >= 0) {
	  if (ion->nele == 2) {
	    if (r.ncomplex[0] == '1' && r.ncomplex[1] == '*') {
	      ion->ibase[p] = ion->iground;
	    }
	  } else if (ion->nele == 1) {
	    ion->ibase[p] = ion->iground;
	  }
	}
	ion->sw[p] = 0;
	if (ion->nele == 2) {
	  if (r.ncomplex[0] == '1' && r.ncomplex[1] == '*') {
	    int nn, kk;
	    nn = ion->vnl[p]/100;
	    kk = 2*(ion->vnl[p]%100);
	    if (kk == 0) {
	      if (ion->j[p] == 0) {
		ion->sw[p] = 1;
	      } else {
		ion->sw[p] = 3;
	      }
	    } else if (ion->j[p] != kk) {
	      ion->sw[p] = 3;
	    } else if (ion->vnl[p] < 10000) {
	      if (s3nk[ion->vnl[p]] == 0) {
		ion->sw[p] = 3;
		s3nk[ion->vnl[p]] = 3;
	      } else {
		ion->sw[p] = 1;
	      }
	    }
	  }
	}
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
    LoadSW(ion);
    ion1 = ion;
    free(rionized);
    FCLOSE(f);
    
    /* determine the minimum ilev in each block */
    blkp = NULL;
    double emin = 0;
    ion->ground = 0;
    for (i = 0; i < ion->nlevels; i++) {
      if (ion->energy[i] < emin) {
	ion->ground = i;
	emin = ion->energy[i];
      }
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
    k = 2*k*(k+2);
    bmatrix = (double *) malloc(sizeof(double)*k);
  }
  
  return 0;
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

int StrNComplex(char *s, NCOMPLEX *c) {
  int i;
  char a[8];
  
  i = 0;
  s[0] = '\0';
  while (i < MAXNCOMPLEX && c[i].n) {
    sprintf(a, "%d*%d.", c[i].n, c[i].nq);
    strcat(s, a);
    i++;
  }
  int n = strlen(s);
  s[n-1] = '\0';
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
  return 0;
}

void SetRateMultiplier(int nele, int t, double a) {
  ION *ion;
  int i;
  
  for (i = 0; i < ions->dim; i++) {
    ion = (ION *) ArrayGet(ions, i);
    if (i == 0 && nele < ion->nele) {
      if (ion0.nele == nele) {
	switch (t) {
	case 0:
	  ion0.ace = a;
	  break;
	case 1:
	  ion0.atr = a;
	  break;
	case 2:
	  ion0.aci = a;
	  break;
	case 3:
	  ion0.arr = a;
	  break;
	case 4:
	  ion0.aai = a;
	  break;
	case 5:
	  ion0.acx = a;
	  break;
	default:
	  break;
	}
	break;
      }
      break;
    }
    if (ion->nele == nele) {
      switch (t) {
      case 0:
	ion->ace = a;
	break;
      case 1:
	ion->atr = a;
	break;
      case 2:
	ion->aci = a;
	break;
      case 3:
	ion->arr = a;
	break;
      case 4:
	ion->aai = a;
	break;
      case 5:
	ion->acx = a;
	break;
      default:
	break;
      }
      break;
    }
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
      blk1->rc0[k] = 0.0;
      blk1->rc1[k] = 0.0;
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
	    blk1->rc0[ion->ilev[i]] += ion->KLN_ai[p];
	  }
	}
      }
    }
    if (electron_density > 0.0) {
      ResetWidMPI();
#pragma omp parallel default(shared) private(brts, blk1, blk2, m, r, p, j)
      {
      int w = 0;
      for (p = 0; p < ion->ce_rates->dim; p++) {
	brts = (BLK_RATE *) ArrayGet(ion->ce_rates, p);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	for (m = 0; m < brts->rates->dim; m++) {
	  if (SkipWMPI(w++)) continue;
	  r = (RATE *) ArrayGet(brts->rates, m);
	  j = ion->ilev[r->i];
	  double zd = electron_density*r->dir;
	  double de = fabs(ion->energy[r->f]-ion->energy[r->i]);
	  de *= RATE_AU*_ce_fbr;
#pragma omp atomic
	  blk1->total_rate[j] += zd;
	  zd = LimitImpactWidth(zd, de);
#pragma omp atomic
	  blk1->rc1[j] += zd;
	  if (r->inv > 0.0) {
	    j = ion->ilev[r->f];
	    double zi = electron_density*r->inv;
#pragma omp atomic
	    blk2->total_rate[j] += zi;
	    zi = LimitImpactWidth(zi, de);
#pragma omp atomic
	    blk2->rc1[j] += zi;
	  }
	}	
      }
      }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(brts, blk1, blk2, m, r, p, j, a, b)
    {
    int w = 0;
    for (p = 0; p < ion->tr_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr_rates, p);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m);
	j = ion->ilev[r->i];
#pragma omp atomic
	blk1->total_rate[j] += r->dir;
#pragma omp atomic
	blk1->rc0[j] += r->dir;
#pragma omp atomic
	blk1->n[j] += r->dir;
	if (r->inv > 0.0 && photon_density > 0.0) {
	  a = photon_density * r->inv;
	  b = a * (ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
#pragma omp atomic
	  blk1->total_rate[j] += b;
	  j = ion->ilev[r->f];
#pragma omp atomic
	  blk2->total_rate[j] += a;
	}
      }
    }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(brts, blk1, blk2, m, r, p, j)
    {
    int w = 0;
    for (p = 0; p < ion->tr2_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr2_rates, p);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m);
	j = ion->ilev[r->i];
#pragma omp atomic
	blk1->total_rate[j] += r->dir;
#pragma omp atomic
	blk1->rc0[j] += r->dir;
#pragma omp atomic
	blk1->n[j] += r->dir;
      }
    }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(brts, blk1, blk2, m, r, p, j)
    {
    int w = 0;
    for (p = 0; p < ion->rr_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->rr_rates, p);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m);
	j = ion->ilev[r->i];
	if (electron_density > 0.0) {
	  double zd = electron_density*r->dir;
#pragma omp atomic
	  blk1->total_rate[j] += zd;
	  double de = fabs(ion->energy[r->f]-ion->energy[r->i]);
	  de *=  RATE_AU * _ce_fbr;
	  zd = LimitImpactWidth(zd, de);
#pragma omp atomic
	  blk1->rc1[j] += zd;
	}
	if (r->inv > 0.0 && photon_density > 0.0) {
	  j = ion->ilev[r->f];
	  double zi = photon_density * r->inv;
#pragma omp atomic
	  blk2->total_rate[j] += zi;
	}
      }
    }
    }
    if (cxt_density > 0) {
      ResetWidMPI();
#pragma omp parallel default(shared) private(brts, blk1, blk2, m, r, p, j)
      {
	int w = 0;
	for (p = 0; p < ion->cx_rates->dim; p++) {
	  brts = (BLK_RATE *) ArrayGet(ion->cx_rates, p);
	  blk1 = brts->iblock;
	  blk2 = brts->fblock;
	  for (m = 0; m < brts->rates->dim; m++) {
	    if (SkipWMPI(w++)) continue;
	    r = (RATE *) ArrayGet(brts->rates, m);
	    j = ion->ilev[r->i];
#pragma omp atomic
	    blk1->total_rate[j] += cxt_density * r->dir;
	  }
	}
      }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(brts, blk1, blk2, m, r, p, j)
    {
    int w = 0;
    for (p = 0; p < ion->ai_rates->dim; p++) {
      brts = (BLK_RATE *) ArrayGet(ion->ai_rates, p);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m);
	j = ion->ilev[r->i];
#pragma omp atomic
	blk1->total_rate[j] += r->dir;
#pragma omp atomic
	blk1->rc0[j] += r->dir;
#pragma omp atomic
	blk1->n[j] += r->dir;
	if (r->inv > 0.0 && electron_density > 0.0) {
	  j = ion->ilev[r->f];
	  double zi = electron_density*r->inv;
#pragma omp atomic
	  blk2->total_rate[j] += zi;
	  double de =  fabs(ion->energy[r->f]-ion->energy[r->i]);
	  de *= RATE_AU * _ce_fbr;
	  zi = LimitImpactWidth(zi, de);
#pragma omp  atomic
	  blk2->rc1[j] += zi;
	}
      }
    }
    }
    if (electron_density > 0.0) {
      ResetWidMPI();
#pragma omp parallel default(shared) private(brts, blk1, blk2, m, r, p, j)
      {
      int w = 0;
      for (p = 0; p < ion->ci_rates->dim; p++) {
	brts = (BLK_RATE *) ArrayGet(ion->ci_rates, p);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	for (m = 0; m < brts->rates->dim; m++) {
	  if (SkipWMPI(w++)) continue;
	  r = (RATE *) ArrayGet(brts->rates, m);
	  j = ion->ilev[r->i];
	  double zd = electron_density * r->dir;
	  double de = fabs(ion->energy[r->f]-ion->energy[r->i]);
	  de *= RATE_AU * _ce_fbr;
#pragma omp atomic
	  blk1->total_rate[j] += zd;
	  zd = LimitImpactWidth(zd, de);
#pragma omp atomic
	  blk1->rc1[j] += zd;
	  if (r->inv > 0.0) {
	    j = ion->ilev[r->f];
	    double zi = electron_density*electron_density*r->inv;	    
#pragma omp atomic
	    blk2->total_rate[j] += zi;
	  }
	}
      }
      }
    }
  }

  m = -2;
  for (i = 0; i < blocks->dim; i++) {
    blk1 = (LBLOCK *) ArrayGet(blocks, i);
    blk1->izr = 0;
    for (k = 0; k < blk1->nlevels; k++) {
      if (blk1->n[k]) {
	blk1->n[k] = 0.0;
      } else {
	if (blk1->iion != m) {
	  if (blk1->nlevels > 1 && i > 0) {
	    blk1->total_rate[k] = 0.0;
	  }
	}
      }
      if (!blk1->total_rate[k]) blk1->izr++;
      m = blk1->iion;
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
  TFILE *f;

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
    MultiInit(&ce, sizeof(double), 3, ablks, "crm_ce");
    MultiInit(&tr, sizeof(double), 3, ablks, "crm_tr");
    MultiInit(&ci, sizeof(double), 3, ablks, "crm_ci");
    MultiInit(&rr, sizeof(double), 3, ablks, "crm_rr");
    MultiInit(&ai, sizeof(double), 3, ablks, "crm_ai");
    MultiInit(&cep, sizeof(double), 3, ablks, "crm_cep");
    MultiInit(&trp, sizeof(double), 3, ablks, "crm_trp");
    MultiInit(&cip, sizeof(double), 3, ablks, "crm_cip");
    MultiInit(&rrp, sizeof(double), 3, ablks, "crm_rrp");
    MultiInit(&aip, sizeof(double), 3, ablks, "crm_aip");
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
	  dce[i][j][q] = 0.0;
	  dtr[i][j][q] = 0.0;
	  drr[i][j][q] = 0.0;
	  dci[i][j][q] = 0.0;
	  dai[i][j][q] = 0.0;
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
      if (blk == NULL) continue;
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
	if (blk == NULL || blk1 == NULL) continue;
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&ce, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&cep, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	      LOCK *lock = NULL;
	      d = (double *) MultiSet(&ce, index, NULL, &lock,
				      InitDoubleData, NULL);
	      if (lock) SetLock(lock);
	      *d += rtmp;
	      if (lock) ReleaseLock(lock);
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
	      LOCK *lock = NULL;
	      d = (double *) MultiSet(&cep, index, NULL, &lock,
				      InitDoubleData, NULL);
	      if (lock) SetLock(lock);
	      *d += rtmp;
	      if (lock) ReleaseLock(lock);
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
	if (blk == NULL || blk1 == NULL) continue;
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&tr, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&trp, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	      LOCK *lock = NULL;
	      d = (double *) MultiSet(&tr, index, NULL, &lock,
				      InitDoubleData, NULL);
	      if (lock) SetLock(lock);
	      *d += rtmp;
	      if (lock) ReleaseLock(lock);
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
	      LOCK *lock = NULL;
	      d = (double *) MultiSet(&trp, index, NULL, &lock,
				      InitDoubleData, NULL);
	      if (lock) SetLock(lock);
	      *d += rtmp;
	      if (lock) ReleaseLock(lock);
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
	if (blk == NULL || blk1 == NULL) continue;
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&tr, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&trp, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	if (blk == NULL || blk1 == NULL) continue;
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&rr, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&rrp, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	      LOCK *lock = NULL;
	      d = (double *) MultiSet(&rr, index, NULL, &lock,
				      InitDoubleData, NULL);
	      if (lock) SetLock(lock);
	      *d += rtmp;
	      if (lock) ReleaseLock(lock);
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
	      LOCK *lock = NULL;
	      d = (double *) MultiSet(&rrp, index, NULL, &lock,
				      InitDoubleData, NULL);
	      if (lock) SetLock(lock);
	      *d += rtmp;
	      if (lock) ReleaseLock(lock);
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
	if (blk == NULL || blk1 == NULL) continue;
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&ai, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&aip, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	      LOCK *lock = NULL;
	      d = (double *) MultiSet(&ai, index, NULL, &lock,
				      InitDoubleData, NULL);
	      if (lock) SetLock(lock);
	      *d += rtmp;
	      if (lock) ReleaseLock(lock);
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
	      LOCK *lock = NULL;
	      d = (double *) MultiSet(&aip, index, NULL, &lock,
				      InitDoubleData, NULL);
	      if (lock) SetLock(lock);
	      *d += rtmp;
	      if (lock) ReleaseLock(lock);
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
	if (blk == NULL || blk1 == NULL) continue;
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&ci, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	    LOCK *lock = NULL;
	    d = (double *) MultiSet(&cip, index, NULL, &lock,
				    InitDoubleData, NULL);
	    if (lock) SetLock(lock);
	    *d += rtmp;
	    if (lock) ReleaseLock(lock);
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
	      LOCK *lock = NULL;
	      d = (double *) MultiSet(&ci, index, NULL, &lock,
				      InitDoubleData, NULL);
	      if (lock) SetLock(lock);
	      *d += rtmp;
	      if (lock) ReleaseLock(lock);
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
	      LOCK *lock = NULL;
	      d = (double *) MultiSet(&cip, index, NULL, &lock,
				      InitDoubleData, NULL);
	      if (lock) SetLock(lock);
	      *d += rtmp;
	      if (lock) ReleaseLock(lock);
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
	      d = (double *) MultiGet(&ce, index, NULL);
	      if (d && *d) {
		rt.ce = *d;
	      }
	      d = (double *) MultiGet(&tr, index, NULL);
	      if (d && *d) {
		rt.tr = *d;
	      }
	    } else {
	      d = (double *) MultiGet(&rr, index, NULL);
	      if (d && *d) {
		rt.rr = *d;
	      }
	      d = (double *) MultiGet(&ai, index, NULL);
	      if (d && *d) {
		rt.ai = *d;
	      }
	      d = (double *) MultiGet(&ci, index, NULL);
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
	      d = (double *) MultiGet(&cep, index, NULL);
	      if (d && *d) {
		rt.ce = *d;
	      }
	      d = (double *) MultiGet(&trp, index, NULL);
	      if (d && *d) {
		rt.tr = *d;
	      }
	    } else {
	      d = (double *) MultiGet(&rrp, index, NULL);
	      if (d && *d) {
		rt.rr = *d;
	      }
	      d = (double *) MultiGet(&aip, index, NULL);
	      if (d && *d) {
		rt.ai = *d;
	      }
	      d = (double *) MultiGet(&cip, index, NULL);
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
      rt.iblock = i;
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
	    d = (double *) MultiGet(&ce, index, NULL);
	    if (d && *d) {
	      rt.ce = *d;
	    } 
	    d = (double *) MultiGet(&tr, index, NULL);
	    if (d && *d) {
	      rt.tr = *d;
	    }
	  } else {
	    d = (double *) MultiGet(&rr, index, NULL);
	    if (d && *d) {
	      rt.rr = *d;
	    }
	    d = (double *) MultiGet(&ai, index, NULL);
	    if (d && *d) {
	      rt.ai = *d;
	    }
	    d = (double *) MultiGet(&ci, index, NULL);
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
	  blk1 = (LBLOCK *) ArrayGet(blocks, j);
	  if (abs(blk1->iion - blk->iion) > 1) continue;
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
	    d = (double *) MultiGet(&ce, index, NULL);
	    if (d && *d) {
	      rt.ce = *d;
	    } 
	    d = (double *) MultiGet(&tr, index, NULL);
	    if (d && *d) {
	      rt.tr = *d;
	    }
	  } else {
	    d = (double *) MultiGet(&rr, index, NULL);
	    if (d && *d) {
	      rt.rr = *d;
	    }
	    d = (double *) MultiGet(&ai, index, NULL);
	    if (d && *d) {
	      rt.ai = *d;
	    }
	    d = (double *) MultiGet(&ci, index, NULL);
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
    MultiFreeLock(&ce, NULL);
    MultiFreeLock(&tr, NULL);
    MultiFreeLock(&ci, NULL);
    MultiFreeLock(&rr, NULL);
    MultiFreeLock(&ai, NULL);
    MultiFreeLock(&cep, NULL);
    MultiFreeLock(&trp, NULL);
    MultiFreeLock(&cip, NULL);
    MultiFreeLock(&rrp, NULL);
    MultiFreeLock(&aip, NULL);
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

void FixNorm(int m) {
  LBLOCK *blk1;
  ION *ion;
  int k0, k1, iion, k, p, i, n;
  double den, *x;

  n = blocks->dim;
  x = bmatrix + n*n;
  for (i = 0; i < n; i++) x[i] = 0.0;

  if (norm_mode == 3) {
    x[0] = 1.0;
    p = 0;
    for (k = 0; k < n; k++) {
      bmatrix[p] = 1.0;
      p += n;
    }
    for (k = 0; k < ions->dim; k++) {
      ion = (ION *) ArrayGet(ions, k);
      den = ion->n0;
      if (den+1 != 1) {
	blk1 = ion->iblock[0];
	if (blk1 != NULL) {
	  p = blk1->ib;
	  x[p] = den;
	}
      }
    }
  } else if (norm_mode == 2) {
    den = 0.0;
    if (ion0.n0 > 0) den += ion0.n0;
    for (i = 0; i < ions->dim; i++) {
      ion = (ION *) ArrayGet(ions, i);
      if (ion->n0 > 0) {
	den += ion->n0;
      }
    }
    x[0] = den;
    p = 0;
    for (k = 0; k < n; k++) {
      bmatrix[p] = 1.0;
      p += n;
    }
  } else {  
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
	      if (norm_mode == 1) {
		if (k < k1 && k >= k0) bmatrix[p] = 1.0;
		else bmatrix[p] = 0.0;
	      } else {
		if (k == k0) bmatrix[p] = 1.0;
		else bmatrix[p] = 0.0;
	      }
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
	if (norm_mode == 1) {
	  if (k < k1 && k >= k0) bmatrix[p] = 1.0;
	  else bmatrix[p] = 0.0;
	} else {       
	  if (k == k0) bmatrix[p] = 1.0;
	  else bmatrix[p] = 0.0;
	}
	p += n;
      }
    } 
  }
}

int BlockMatrix(void) {
  ION *ion;
  RATE *r;
  LBLOCK *blk1, *blk2;
  BLK_RATE *brts;
  int n, k, m, i, j, t, p, q;
  double a, den, *rex;
  
  n = blocks->dim;
  for (i = 0; i < 2*n*(n+2); i++) {
    bmatrix[i] = 0.0;    
  }
  rex = bmatrix + 2*n*(n+1) + n;
  
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    if (electron_density > 0.0) {
      ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, den)
      {
      int w = 0;
      for (t = 0; t < ion->ce_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->ce_rates, t);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	if (blk1 == blk2 && blk1->izr == 0) continue;
	if (rec_cascade && (blk1->rec || blk2->rec)) continue;
	for (m = 0; m < brts->rates->dim; m++) {
	  if (SkipWMPI(w++)) continue;
	  r = (RATE *) ArrayGet(brts->rates, m);
	  if (ion->iblock[r->i] == NULL ||
	      ion->iblock[r->f] == NULL) continue;	      
	  i = ion->iblock[r->i]->ib;
	  j = ion->iblock[r->f]->ib;
	  den = blk1->r[ion->ilev[r->i]];
	  if (den) {
	    p = i*n + j;
	    if (blk1 == blk2) {
	      if (blk1->total_rate[ion->ilev[r->f]] == 0) {
#pragma omp atomic
		rex[i] += den*electron_density*r->dir;
	      }
	    } else {
	      if (blk2->total_rate[ion->ilev[r->f]]) {
#pragma omp atomic
		bmatrix[p] += den * electron_density * r->dir;
	      } else {
#pragma omp atomic
		rex[i] += den*electron_density*r->dir;
	      }
	    }
	  }
	  if (r->inv > 0.0) {
	    den = blk2->r[ion->ilev[r->f]];
	    if (den) {
	      p = i + j*n;
	      if (blk1 == blk2) {
		if (blk1->total_rate[ion->ilev[r->i]] == 0) {
#pragma omp atomic
		  rex[j] += den*electron_density*r->inv;
		}
	      } else {
		if (blk1->total_rate[ion->ilev[r->i]]) {
#pragma omp atomic
		  bmatrix[p] += den * electron_density * r->inv;
		} else {
#pragma omp atomic
		  rex[j] += den*electron_density*r->inv;
		}
	      }
	    }
	  }
	}
      }
      }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, den, a)
    {
    int w = 0;
    for (t = 0; t < ion->tr_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (blk1 == blk2 && blk1->izr == 0) continue;
      if (rec_cascade && (blk1->rec || blk2->rec)) continue;
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m);
	if (ion->iblock[r->i] == NULL ||
	    ion->iblock[r->f] == NULL) continue;
	i = ion->iblock[r->i]->ib;
	j = ion->iblock[r->f]->ib;
	den = blk1->r[ion->ilev[r->i]];
	if (den) {
	  p = i*n + j;
	  if (blk1 == blk2) {
	    if (blk1->total_rate[ion->ilev[r->f]] == 0) {
#pragma omp atomic
	      rex[i] += den*r->dir;
	    }
	  } else {
	    if (blk2->total_rate[ion->ilev[r->f]]) {
#pragma omp atomic
	      bmatrix[p] += den * r->dir;
	    } else {
#pragma omp atomic
	      rex[i] += den*r->dir;
	    }
	  }
	}
	if (r->inv > 0.0 && photon_density > 0.0) {
	  if (den) {
	    a = photon_density * r->inv;
	    if (blk1 == blk2) {
	      if (blk1->total_rate[ion->ilev[r->f]] == 0) {
#pragma omp atomic
		rex[i] += den*a*(ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
	      }
	    } else {
	      if (blk2->total_rate[ion->ilev[r->f]]) {
#pragma omp atomic
		bmatrix[p] += den*a*(ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
	      } else {
#pragma omp atomic
		rex[i] += den*a*(ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
	      }
	    }
	  }	  
	  den = blk2->r[ion->ilev[r->f]];
	  if (den) {
	    p = i + j*n;
	    if (blk1 == blk2) {
	      if (blk1->total_rate[ion->ilev[r->i]] == 0) {
#pragma omp atomic
		rex[j] += den*photon_density*r->inv;
	      }
	    } else {
	      if (blk1->total_rate[ion->ilev[r->i]]) {
#pragma omp atomic
		bmatrix[p] += den * photon_density * r->inv;
	      } else {
#pragma omp atomic
		rex[j] += den*photon_density*r->inv;
	      }
	    }
	  }
	}
      }
    }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, den)
    {
    int w = 0;
    for (t = 0; t < ion->tr2_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr2_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (blk1 == blk2 && blk1->izr == 0) continue;
      if (rec_cascade && (blk1->rec || blk2->rec)) continue;
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m);
	if (ion->iblock[r->i] == NULL ||
	    ion->iblock[r->f] == NULL) continue;
	i = ion->iblock[r->i]->ib;
	j = ion->iblock[r->f]->ib;
	den = blk1->r[ion->ilev[r->i]];
	if (den) {
	  p = i*n + j;
	  if (blk1 == blk2) {
	    if (blk1->total_rate[ion->ilev[r->f]] == 0) {
#pragma omp atomic
	      rex[i] += den*r->dir;
	    }
	  } else {
	    if (blk2->total_rate[ion->ilev[r->f]]) {
#pragma omp atomic
	      bmatrix[p] += den * r->dir;
	    } else {
#pragma omp atomic
	      rex[i] += den*r->dir;
	    }
	  }
	}
      }
    }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, den)
    {
    int w = 0;
    for (t = 0; t < ion->rr_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->rr_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (blk1 == blk2 && blk1->izr == 0) continue;
      if (rec_cascade && (blk1->rec || blk2->rec)) continue;
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m); 
	if (ion->iblock[r->i] == NULL ||
	    ion->iblock[r->f] == NULL) continue;
	i = ion->iblock[r->i]->ib;
	j = ion->iblock[r->f]->ib;
	den = blk1->r[ion->ilev[r->i]];
	if (den) {
	  if (electron_density > 0.0) {
	    p = i*n + j;
	    if (blk1 == blk2) {
	      if (blk1->total_rate[ion->ilev[r->f]] == 0) {
#pragma omp atomic
		rex[i] += den*electron_density*r->dir;
	      }
	    } else {
	      if (blk2->total_rate[ion->ilev[r->f]]) {
#pragma omp atomic
		bmatrix[p] += den * electron_density * r->dir;
	      } else {
#pragma omp atomic
		rex[i] += den*electron_density*r->dir;
	      }
	    }
	  }
	}
	if (r->inv > 0.0 && photon_density > 0.0) {
	  den = blk2->r[ion->ilev[r->f]];
	  if (den) {
	    p = i + j*n;
	    if (blk1 == blk2) {
	      if (blk1->total_rate[ion->ilev[r->i]] == 0) {
#pragma omp atomic
		rex[j] += den*photon_density*r->inv;
	      }
	    } else {
	      if (blk1->total_rate[ion->ilev[r->i]]) {
#pragma omp atomic
		bmatrix[p] += den * photon_density * r->inv;
	      } else {
#pragma omp atomic
		rex[j] += den*photon_density*r->inv;
	      }
	    }
	  }
	}
      }
    }
    }
    if (cxt_density > 0) {
      ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, den)
      {
	int w = 0;
	for (t = 0; t < ion->cx_rates->dim; t++) {
	  brts = (BLK_RATE *) ArrayGet(ion->cx_rates, t);
	  blk1 = brts->iblock;
	  blk2 = brts->fblock;
	  if (blk1 == blk2 && blk1->izr == 0) continue;
	  for (m = 0; m < brts->rates->dim; m++) {
	    if (SkipWMPI(w++)) continue;
	    r = (RATE *) ArrayGet(brts->rates, m); 
	    if (ion->iblock[r->i] == NULL ||
		ion->iblock[r->f] == NULL) continue;
	    i = ion->iblock[r->i]->ib;
	    j = ion->iblock[r->f]->ib;
	    den = blk1->r[ion->ilev[r->i]];
	    if (den) {
	      p = i*n + j;
	      if (blk1 == blk2) {
		if (blk1->total_rate[ion->ilev[r->f]] == 0) {
#pragma omp atomic
		  rex[i] += den*cxt_density*r->dir;
		}
	      } else {
		if (blk2->total_rate[ion->ilev[r->f]]) {
#pragma omp atomic
		  bmatrix[p] += den * cxt_density * r->dir;
		} else {
#pragma omp atomic
		  rex[i] += den*cxt_density*r->dir;
		}
	      }
	    }
	  }
	}
      }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, den)
    {
    int w = 0;
    for (t = 0; t < ion->ai_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->ai_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (blk1 == blk2 && blk1->izr == 0) continue;
      if (rec_cascade && (blk1->rec || blk2->rec)) continue;
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m); 
	if (ion->iblock[r->i] == NULL ||
	    ion->iblock[r->f] == NULL) continue;
	i = ion->iblock[r->i]->ib;
	j = ion->iblock[r->f]->ib;
	den = blk1->r[ion->ilev[r->i]];
	if (den) {
	  p = i*n + j;
	  if (blk1 == blk2) {
	    if (blk1->total_rate[ion->ilev[r->f]] == 0) {
#pragma omp atomic
	      rex[i] += den*r->dir;
	    }
	  } else {
	    if (blk2->total_rate[ion->ilev[r->f]]) {
#pragma omp atomic
	      bmatrix[p] += den * r->dir;
	    } else {
#pragma omp atomic
	      rex[i] += den*r->dir;
	    }
	  }
	}
	if (r->inv > 0.0 && electron_density > 0.0) {
	  den = blk2->r[ion->ilev[r->f]];
	  if (den) {
	    p = i + j*n;
	    if (blk1 == blk2) {
	      if (blk1->total_rate[ion->ilev[r->i]] == 0) {
#pragma omp atomic
		rex[j] += den * electron_density * r->inv;
	      }
	    } else {
	      if (blk1->total_rate[ion->ilev[r->i]]) {
#pragma omp atomic
		bmatrix[p] += den * electron_density * r->inv;
	      } else {
#pragma omp atomic
		rex[j] += den * electron_density * r->inv;
	      }
	    }
	  }
	}
      }  
    }
    }
    if (electron_density > 0.0) {
      ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, den)
      {
      int w = 0;
      for (t = 0; t < ion->ci_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->ci_rates, t);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	if (blk1 == blk2 && blk1->izr == 0) continue;
	if (rec_cascade && (blk1->rec || blk2->rec)) continue;
	for (m = 0; m < brts->rates->dim; m++) {
	  if (SkipWMPI(w++)) continue;
	  r = (RATE *) ArrayGet(brts->rates, m); 
	  if (ion->iblock[r->i] == NULL ||
	      ion->iblock[r->f] == NULL) continue;
	  i = ion->iblock[r->i]->ib;
	  j = ion->iblock[r->f]->ib;
	  den = blk1->r[ion->ilev[r->i]];
	  if (den) {
	    p = i*n + j;
	    if (blk1 == blk2) {
	      if (blk1->total_rate[ion->ilev[r->f]] == 0) {
#pragma omp atomic
		rex[i] += den * electron_density * r->dir;
	      }
	    } else {
	      if (blk2->total_rate[ion->ilev[r->f]]) {
#pragma omp atomic
		bmatrix[p] += den * electron_density * r->dir;
	      } else {
#pragma omp atomic
		rex[i] += den * electron_density * r->dir;
	      }
	    }
	  }
	  if (r->inv > 0.0) {
	    den = blk2->r[ion->ilev[r->f]];
	    if (den) {
	      p = i + j*n;	      
	      den *= electron_density;
	      if (blk1 == blk2) {
		if (blk1->total_rate[ion->ilev[r->i]] == 0) {
#pragma omp atomic
		  rex[j] += den * electron_density * r->inv;
		}
	      } else {
		if (blk1->total_rate[ion->ilev[r->i]]) {
#pragma omp atomic
		  bmatrix[p] += den * electron_density * r->inv;
		} else {
#pragma omp atomic
		  rex[j] += den * electron_density * r->inv;
		}
	      }
	    }
	  }
	}  
      }
      }
    }
  }

  for (i = 0; i < n; i++) {
    p = i*n;
    q = i + p;
    bmatrix[q] = rex[i];
    for (j = 0; j < n; j++) {
      if (j != i) bmatrix[q] += bmatrix[p];
      p++;
    }
    bmatrix[q] = - bmatrix[q];
  }

  return 0;
}

int BlockPopulation(int miter) {
  LBLOCK *blk;
  ION *ion;
  double *b, *x;
  double *a;
  int *ipiv;
  int info;
  int n, m;
  int nrhs;
  int lda, ldb;
  int i, j, p, q, ntd, nb, niter;
  double ta, tb, td;

  n = blocks->dim;
  a = bmatrix + n*n;
  x = a;
  a = a + n;
  b = a + n*n;
  ipiv = (int *) (b+n);

  /*if (norm_mode == 0) miter = 1; use iteration to enforec normalization*/  
  miter = 1; /*disable iteration, use FixNorm to enforce normalization*/
  for (niter = 0; niter < miter; niter++) {
    FixNorm(niter);
    ta = 0.0;
    p = 0;
    q = 0;
    m = 0;
    for (i = 0; i < n; i++) {
      if (bmatrix[i+i*n] == 0) {
	x[i] = 2E50;
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
	if (x[j] > 1E50) {
	  continue;
	}
	a[q] = a[p+j];
	q++;
      }
      p += n;
    }
    
    q = 0;
    for (j = 0; j < n; j++) {
      if (x[j] > 1E50) {
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
      printf("Error in solving BlockMatrix: %d\n", info);
      exit(1);
    }

    p = 0;
    for (i = 0; i < n; i++) {
      blk = (LBLOCK *) ArrayGet(blocks, i);
      if (rec_cascade && blk->rec) continue;
      if (x[i] > 1E50) {
	blk->nb = 0.0;
	for (j = 0; j < blk->nlevels; j++) {
	  blk->n[j] = 0.0;
	}
      } else {
	blk->nb = b[p++];
	if (blk->nb < 0) blk->nb = 0.0;
      }
    }

    if (niter == miter-1) break;

    p = 0;
    q = -1;
    tb = 0.0;
    td = 0.0;
    ntd = 0;
    for (i = 0; i < n; i++) {
      blk = (LBLOCK *) ArrayGet(blocks, i);
      if (rec_cascade && blk->rec) continue;
      if (blk->iion != q) {
	if (q == -1) {
	  ion0.nt = tb;
	  if (ion0.n > 0 && miter > 1) {
	    ion0.n0 *= ion0.n/tb;
	    td += fabs(tb - ion0.n)/ion0.n;
	    ntd++;
	  }
	} else {
	  ion = (ION *) ArrayGet(ions, q);
	  ion->nt = tb;
	  if (ion->n > 0 && miter > 1) {
	    ion->n0 *= ion->n/tb;
	    td += fabs(tb - ion->n)/ion->n;
	    ntd++;
	  }
	}
	q = blk->iion;
	tb = 0.0;
      }
      tb += blk->nb;
    }
    ion = (ION *) ArrayGet(ions, q);
    ion->nt = tb;
    if (ion->n > 0 && miter > 1) {
      td += fabs(tb - ion->n)/ion->n;
      ion->n0 *= ion->n/tb;
      ntd++;
    }
    if (ntd > 0) {
      td /= ntd;
      if (td < iter_accuracy) break;
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
      ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, q)
      {
      int w = 0;
      for (t = 0; t < ion->ce_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->ce_rates, t);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	if (rec_cascade && iter >= 0) {
	  if (blk1->rec || blk2->rec) continue;
	} 
	for (m = 0; m < brts->rates->dim; m++) {
	  if (SkipWMPI(w++)) continue;
	  r = (RATE *) ArrayGet(brts->rates, m);
	  if (ion->iblock[r->i] == NULL ||
	      ion->iblock[r->f] == NULL) continue;
	  i = ion->iblock[r->i]->ib;
	  p = ion->ilev[r->i];
	  j = ion->iblock[r->f]->ib;
	  q = ion->ilev[r->f];
	  if (blk1->r[p]) {
#pragma omp atomic
	    blk2->n[q] += blk1->r[p] * electron_density * r->dir;
	  }
	  if (r->inv > 0.0) {
	    if (blk2->r[q]) {
#pragma omp atomic
	      blk1->n[p] += blk2->r[q] * electron_density * r->inv;
	    }
	  }
	}
      }
      }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, q)
    {
    int w = 0;
    for (t = 0; t < ion->tr_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (rec_cascade && iter >= 0) {
	if (blk1->rec || blk2->rec) continue;
      }
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m);
	if (ion->iblock[r->i] == NULL ||
	    ion->iblock[r->f] == NULL) continue;
	i = ion->iblock[r->i]->ib;
	p = ion->ilev[r->i];
	j = ion->iblock[r->f]->ib;
	q = ion->ilev[r->f];    
	if (blk1->r[p]) {
#pragma omp atomic
	  blk2->n[q] += blk1->r[p] * r->dir;
	}
	if (r->inv > 0.0 && photon_density > 0.0) {
	  a = photon_density * r->inv;	  
	  if (blk1->r[p]) {
#pragma omp atomic
	    blk2->n[q] += blk1->r[p]*a*(ion->j[r->f]+1.0)/(ion->j[r->i]+1.0);
	  }	  
	  if (blk2->r[q]) {
#pragma omp atomic
	    blk1->n[p] += blk2->r[q] * a;
	  }
	}
      }
    }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, q)
    {
    int w = 0;
    for (t = 0; t < ion->tr2_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->tr2_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (rec_cascade && iter >= 0) {
	if (blk1->rec || blk2->rec) continue;
      }
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m);
	if (ion->iblock[r->i] == NULL ||
	    ion->iblock[r->f] == NULL) continue;
	i = ion->iblock[r->i]->ib;
	p = ion->ilev[r->i];
	j = ion->iblock[r->f]->ib;
	q = ion->ilev[r->f];    
	if (blk1->r[p]) {
#pragma omp atomic
	  blk2->n[q] += blk1->r[p] * r->dir;
	}
      }
    }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, q)
    {
    int w = 0;
    for (t = 0; t < ion->rr_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->rr_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (rec_cascade && iter >= 0) {
	if (blk1->rec || blk2->rec) continue;
      }
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m);
	if (ion->iblock[r->i] == NULL ||
	    ion->iblock[r->f] == NULL) continue;
	i = ion->iblock[r->i]->ib;
	p = ion->ilev[r->i];
	j = ion->iblock[r->f]->ib;
	q = ion->ilev[r->f];    
	if (electron_density > 0.0) {
	  if (blk1->r[p]) {
#pragma omp atomic
	    blk2->n[q] += blk1->r[p] * electron_density * r->dir;
	  }
	} 
	if (r->inv > 0.0 && photon_density > 0.0) {
	  if (blk2->r[q]) {
#pragma omp atomic
	    blk1->n[p] += blk2->r[q] * photon_density * r->inv;
	  }
	}
      }
    }
    }
    if (cxt_density > 0) {
      ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, q)
      {
	int w = 0;
	for (t = 0; t < ion->cx_rates->dim; t++) {
	  brts = (BLK_RATE *) ArrayGet(ion->cx_rates, t);
	  blk1 = brts->iblock;
	  blk2 = brts->fblock;
	  for (m = 0; m < brts->rates->dim; m++) {
	    if (SkipWMPI(w++)) continue;
	    r = (RATE *) ArrayGet(brts->rates, m);
	    if (ion->iblock[r->i] == NULL ||
		ion->iblock[r->f] == NULL) continue;
	    i = ion->iblock[r->i]->ib;
	    p = ion->ilev[r->i];
	    j = ion->iblock[r->f]->ib;
	    q = ion->ilev[r->f];    
	    if (blk1->r[p]) {
#pragma omp atomic
	      blk2->n[q] += blk1->r[p] * cxt_density * r->dir;
	    }
	  } 
	}
      }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, q)
    {
    int w = 0;
    for (t = 0; t < ion->ai_rates->dim; t++) {
      brts = (BLK_RATE *) ArrayGet(ion->ai_rates, t);
      blk1 = brts->iblock;
      blk2 = brts->fblock;
      if (rec_cascade && iter >= 0) {
	if (blk1->rec || blk2->rec) continue;
      }
      for (m = 0; m < brts->rates->dim; m++) {
	if (SkipWMPI(w++)) continue;
	r = (RATE *) ArrayGet(brts->rates, m);
	if (ion->iblock[r->i] == NULL ||
	    ion->iblock[r->f] == NULL) continue;
	i = ion->iblock[r->i]->ib;
	p = ion->ilev[r->i];
	j = ion->iblock[r->f]->ib;
	q = ion->ilev[r->f];    
	if (blk1->r[p]) {
#pragma omp atomic
	  blk2->n[q] += blk1->r[p] * r->dir;
	}
	if (r->inv > 0.0 && electron_density > 0.0) {
	  if (blk2->r[q]) {
#pragma omp atomic
	    blk1->n[p] += blk2->r[q] * electron_density * r->inv;
	  }
	}
      }
    }
    }
    if (electron_density > 0.0) {
      ResetWidMPI();
#pragma omp parallel default(shared) private(t, brts, blk1, blk2, m, r, i, p, j, q)
      {
      int w = 0;
      for (t = 0; t < ion->ci_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->ci_rates, t);
	blk1 = brts->iblock;
	blk2 = brts->fblock;
	if (rec_cascade && iter >= 0) {
	  if (blk1->rec || blk2->rec) continue;
	}
	for (m = 0; m < brts->rates->dim; m++) {
	  if (SkipWMPI(w++)) continue;
	  r = (RATE *) ArrayGet(brts->rates, m);
	  if (ion->iblock[r->i] == NULL ||
	      ion->iblock[r->f] == NULL) continue;
	  i = ion->iblock[r->i]->ib;
	  p = ion->ilev[r->i];
	  j = ion->iblock[r->f]->ib;
	  q = ion->ilev[r->f];    
	  if (blk1->r[p]) {
#pragma omp atomic
	    blk2->n[q] += blk1->r[p] * electron_density * r->dir;
	  }
	  if (r->inv) {
	    if (blk2->r[q]) {
#pragma omp atomic
	      blk1->n[p] += blk2->r[q] * electron_density * 
		electron_density * r->inv;
	    }
	  }
	}
      }
      }
    }
  }

  nlevels = 0;
  d = 0.0;
  td = 0.0;
  double bd = 0.0;
  double tbd = 0.0;
  for (k = 0; k < blocks->dim; k++) {
    blk1 = (LBLOCK *) ArrayGet(blocks, k);
    if (rec_cascade && iter >= 0) {
      if (blk1->rec) continue;
    }
    double tbr = fabs(bmatrix[k*blocks->dim+k]);
    if (blk1->nlevels == 1 && blk1->nb) {
      blk1->r[0] = 1.0;      
      a = blk1->nb;      
      bd += fabs(a-blk1->n0[0])*tbr;
      tbd += fabs(a)*tbr;
      blk1->n[0] = a;
    } else {
      a = 0.0;
      for (m = 0; m < blk1->nlevels; m++) {
	if (blk1->total_rate[m]) {
	  blk1->n[m] /= blk1->total_rate[m];
	  a += blk1->n[m];
	} else {
	  blk1->n[m] = 0.0;
	}
      }
      bd += fabs(a - blk1->nb)*tbr;
      tbd += fabs(a)*tbr;
      if (a) {
	if (norm_mode > 1 && iter >= 0 && blk1->nb > 0) {
	  a = blk1->nb/a;
	  for (m = 0; m < blk1->nlevels; m++) {
	    blk1->n[m] *= a;
	  }
	} else {
	  blk1->nb = a;
	}
      }  
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
      /*
      if (iter > 0) {
	if (blk1->iion < 0) a = ion0.nt;
	else {
	  ion = (ION *) ArrayGet(ions, blk1->iion);
	  a = ion->nt;
	}
      } else a = 1.0;
      */
      a = 0.0;
      m = 0;
      while (blk1->ncomplex[m+1].n > 0) m++;
      int nmx = blk1->ncomplex[m].n;
      for (m = 0; m < blk1->nlevels; m++) {
	if (blk1->n[m] && blk1->rec == NULL &&
	    (_itol_nmax <= 0 || nmx <= _itol_nmax)) {
	  //d += fabs(1.0 - blk1->n0[m]/blk1->n[m]);
	  d += fabs((blk1->n[m]-blk1->n0[m])*blk1->total_rate[m]);
	  td += fabs(blk1->total_rate[m]*blk1->n[m]);
	  nlevels += 1;
	}
	if (iter >= 2) {
	  blk1->n[m] = b*blk1->n0[m] + c*blk1->n[m];
	}
	blk1->n0[m] = blk1->n[m];
	a += blk1->n[m];
      }
      if (a > 0) {
	blk1->nb = a;
	for (m = 0; m < blk1->nlevels; m++) {
	  blk1->r[m] = blk1->n[m]/blk1->nb;
	}
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
      if (p == -1) {
	ion0.nt = a;
      } else {
	ion->nt = a;
      }
      /*
      if (h > 0.0 && norm_mode != 0) {
	if (p == -1) ion0.n0 *= ion0.n/a;
	else ion->n0 *= ion->n/a;
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
  if (ion->n > 0.0 && norm_mode != 0) {
    ion->n0 *= ion->n/a;
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
	  d += fabs(((blk1->n[m] - blk1->n0[m])/blk1->n[m])*blk1->nb);
	  td += fabs(blk1->nb);
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
  bd /= tbd;
  d = Min(d, bd);
  return d;
}

int LevelPopulation(void) {
  int i, n;
  double d, c;

  if (!_silent) {
    printf("Population Iteration:\n");
  }
  d = 10.0;
  c = 1.0;
  double wt00, wt0, wt1;
  wt0 = WallTime();
  wt00 = wt0;
  for (i = 0; i < max_iter; i++) {
    BlockMatrix();
    n = max_iter;
    BlockPopulation(n);
    d = BlockRelaxation(i);
    wt1 = WallTime();
    if (!_silent) {
      if (ProcID() >= 0) {
	printf("%6d %5d %11.4E %11.4E %11.4E %11.4E\n",
	       ProcID(), i, d, wt1-wt0, wt1-wt00, TotalSize());
      } else {
	printf("%5d %11.4E %11.4E %11.4E %11.4E\n",
	       i, d, wt1-wt0, wt1-wt00, TotalSize());
      }
    }
    wt0 = wt1;
    fflush(stdout);
    if (d < iter_accuracy) break;
  }

  if (_silent < 2) {
    if (i == max_iter) {
      if (ProcID() >= 0) {
	printf("%6d max iteration reached %d\n", ProcID(), i);
      } else {
	printf("max iteration reached %d\n", i);
      }
    }
  }
  return 0;
}

int Cascade(void) {
  int i;
  double d;
  
  if (!rec_cascade) return 0;
  double wt0, wt1;
  wt0 = WallTime();
  printf("Cascade  Iteration\n");  
  d = BlockRelaxation(-1);  
  for (i = 1; i <= max_iter; i++) {
    d = BlockRelaxation(-i);
    wt1 = WallTime();
    printf("%5d %11.4E %11.4E\n", i, d, wt1-wt0);
    wt0 = wt1;
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
  TFILE *f;
  double e, a, e0;
  int i, p, q, ib, iuta;
  double smax, s, b, c;
  /* this is the factor in TRRate() *((ev/erg)/c)*1E20 */
  const double factor = 2.52433977E-4;

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
#pragma omp parallel default(shared) private(smax, m, rt, p, q, e, j, c, r, a, b, s, rx, dev)
      {
	int  w = 0;
	smax = 0.0;
	for (m = 0; m < brts->rates->dim; m++) {
	  if (SkipWMPI(w++)) continue;
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
	  if (rrc/10 == 1) {
	    j = ion->ilev[rt->f];
	    c = fblk->n[j];
	  } else {
	    j = ion->ilev[rt->i];
	    c = iblk->n[j];
	  }
	  if (c > 0.0) {
	    r.lower = q;
	    r.upper = p;
	    b = 0.0;
	    if (rrc/10 == 1) {
	      r.energy = -e;
	      a = rt->dir*(ion->j[rt->i]+1.0);
	      e *= HARTREE_EV;
	      a *= factor/(e*e*(ion->j[rt->f]+1.0));
	    } else {
	      r.energy = e;
	      a = rt->dir;
	      if (rt->inv > 0.0 && photon_density > 0.0) {
		b = photon_density * rt->inv;
		b *= (ion->j[rt->f]+1.0)/(ion->j[rt->i]+1.0);
	      }
	    }
	    r.strength = c * (a+b);	    
	    s = r.strength*e;
	    if (s < strength_threshold*smax) continue;
	    if (s > smax) smax = s;
	    if (iuta) {
	      dev = (RATE *) ArrayGet(brdev->rates, m);
	      r.energy = dev->dir;
	      rx.sdev = dev->inv;
	    }
	    if (_sp_trm == 0) {
	      r.rrate = a;
	      r.trate = iblk->total_rate[ion->ilev[rt->i]];
	    } else {
	      r.rrate = iblk->rc0[ion->ilev[rt->i]];
	      r.rrate += fblk->rc0[ion->ilev[rt->f]];
	      r.trate = iblk->rc1[ion->ilev[rt->i]];
	      r.trate += fblk->rc1[ion->ilev[rt->f]];
	    }
	    WriteSPRecord(f, &r, &rx);
	  }
	}
      }
      DeinitFile(f, &fhdr);
    }

    if (!(rrc%10)) continue;
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
#pragma omp parallel default(shared) private(smax, m, rt, p, q, e, j, r, a, b, s, rx)
      {
	int  w = 0;
	smax = 0.0;
	for (m = 0; m < brts->rates->dim; m++) {
	  if (SkipWMPI(w++)) continue;
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

double LimitImpactWidth(double zd, double de) {
  double r;
  r = zd/de;
  if (r < 1e-5) {
    return zd*(1 + 0.5*(_ce_xbr-1)*r);
  } else {
    return (1/_ce_xbr)*de*(pow(1+r,_ce_xbr)-1.0);
  }
}

int SelectLines(char *ifn, char *ofn, int nele, int type, 
		double emin, double emax, double fmin) {
  F_HEADER fh;
  SP_HEADER h;
  SP_RECORD r, *rp;
  SP_EXTRA rx, *rpx;
  ARRAY sp, spx;
  ARRAY linetype;
  TFILE *f1;
  TFILE *f2;
  int n, nb, i;
  int t, t0, t1, t2;
  int r0, r1;
  int low, up;
  double e, a, smax;
  int swp;
  LINETYPE tp, *tt;
  
  rx.sdev = 0.0;

#if USE_MPI == 2
  if (!MPIReady()) InitializeMPI(0, 1);
#endif
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("ERROR: File %s does not exist\n", ifn);
    return -1;
  }
  f2 = FOPEN(ofn, "a");
  if (f2 == NULL) {
    printf("ERROR: Cannot open file %s\n", ofn);
    return -1;
  }

  double mt = GetAtomicMassTable()[(int)(fh.atom)];
  t2 = abs(type) / 1000000;
  if (type < 0) t2 = -1;
  t = abs(type) % 1000000;
  t1 = t / 10000;
  t0 = t % 10000;
  t0 = t0/100; 

  if (fmin >= 1.0) {
    low = emin;
    up = emax;
  }

  int nbk = 1000000;
  ArrayInit(&sp, sizeof(SP_RECORD), nbk);
  ArrayInit(&spx, sizeof(SP_EXTRA), nbk);
  ArrayInit(&linetype, sizeof(LINETYPE), nbk);
  int une[256];
  for (i = 0; i < 256; i++) une[i] = 0;
  smax = 0.0;
  int ne0 = 10000;
  int ne1 = 0;
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadSPHeader(f1, &h, swp);
    if (h.ntransitions == 0) continue;
    if (nele >= 0 && h.nele != nele) goto LOOPEND;
    if (ne0 > h.nele) ne0 = h.nele;
    if (ne1 < h.nele) ne1 = h.nele;
    une[h.nele]++;
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
    } else {
      if (fmin < 0) {
	if (h.type != 0) goto LOOPEND;
      } else {
	if (h.type == 0) goto LOOPEND;
      }
    }
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      if (fmin >= 0) {
	r.energy = fabs(r.energy);
      }
      r.energy *= HARTREE_EV;
      rx.sdev *= HARTREE_EV;
      if (fmin >= 1.0) {
	if (r.lower == low && r.upper == up) {
	  ArrayAppend(&sp, &r, NULL);
	  ArrayAppend(&spx, &rx, NULL);
	  tp.type = h.type;
	  tp.nele = h.nele;
	  ArrayAppend(&linetype, &tp, NULL);
	  break;
	}
      } else {
	if (fmin >= 0) {
	  if (r.energy < emin || r.energy > emax) continue;
	  a = r.energy * r.strength;
	  if (a < smax*fmin) continue;
	  if (a > smax) smax = a;
	}
	ArrayAppend(&sp, &r, NULL);
	ArrayAppend(&spx, &rx, NULL);
	tp.type = h.type;
	tp.nele = h.nele;
	ArrayAppend(&linetype, &tp, NULL);
      }
    }
    continue;

  LOOPEND:
    FSEEK(f1, h.length, SEEK_CUR);
  }
  if (ne0 == 0) ne0 = 1;
  if (ne1 < ne0) ne1 = ne0;
  int zt = (int)(0.1+fh.atom);
  double wd = 0;
  double *wdi = NULL, *wir = NULL, *wrf = NULL, *wid=NULL;
  if (_starknp > 0) {
    wdi = malloc(sizeof(double)*_starknp);
    wir = malloc(sizeof(double)*_starknp);
    wrf = malloc(sizeof(double)*_starknp*zt);
    wid = malloc(sizeof(double)*zt);
  }
  if (_starkqc > 0 && electron_density > 0) {
    int idist;
    DISTRIBUTION *dist;
    dist = GetEleDist(&idist);
    if (idist == 0) {
      PrepStarkQC(mt, electron_density*1e10, dist->params[0], 
		  &wd, wdi, wir, zt, ne0, ne1, wrf, wid);
    }
  }
  char buf[2048];
  if (fmin >= 1.0) {
    if (sp.dim > 0) {
      rp = (SP_RECORD *) ArrayGet(&sp, 0);
      rpx = (SP_EXTRA *) ArrayGet(&spx, 0);
      tt = (LINETYPE *) ArrayGet(&linetype, 0);
      double wtr = rp->rrate;
      double wtt = rp->trate;
      double wtt0 = wtt;
      if (_sp_trm) {	
	wtr *= WCOEF;
	wtt *= WCOEF;
	wtt0 = wtt;
	wtt = CalcStarkQC(wtt0, wd, wdi, wir, wrf, wid, tt->nele, tt->type);
      }      
      sprintf(buf, "%2d %6d %6d %6d %13.6E %11.4E %15.8E %11.4E %11.4E %11.4E %11.4E\n", 
	      tt->nele, rp->lower, rp->upper, tt->type, rp->energy, rpx->sdev, rp->strength, wtr, wtt0, wtt, wd);
      FWRITE(buf, 1, strlen(buf), f2);
    }
  } else {
    if (sp.dim > 0) {
      if (fmin > 0) {
	smax *= fmin;
      } else {
	smax = 0.0;
      }
      int ie;
      for (ie = 1; ie < 256; ie++) {
	if (une[ie] == 0) continue;	
#pragma omp parallel default(shared) private(i, rp, rpx, tt, e, buf)
	{
	  int w = 0;
	  for (i = 0; i < sp.dim; i++) {
	    if (SkipWMPI(w++)) continue;
	    tt = (LINETYPE *) ArrayGet(&linetype, i);
	    if (tt->nele != ie) continue;
	    rp = (SP_RECORD *) ArrayGet(&sp, i);
	    rpx = (SP_EXTRA *) ArrayGet(&spx, i);
	    e = rp->energy;
	    if (fmin < 0 || rp->strength*e > smax) {
	      double wtr = rp->rrate;
	      double wtt = rp->trate;
	      double wtt0 = wtt;
	      if (_sp_trm) {
		wtr *= WCOEF;
		wtt *= WCOEF;
		wtt0 = wtt;
		wtt = CalcStarkQC(wtt0, wd, wdi, wir, wrf, wid, tt->nele, tt->type);
	      }      
	      sprintf(buf, "%2d %6d %6d %6d %13.6E %11.4E %15.8E %11.4E %11.4E %11.4E %11.4E\n", 
		      tt->nele, rp->lower, rp->upper, tt->type, e, rpx->sdev, 
		      rp->strength, wtr, wtt0, wtt, wd);
	      FWRITE(buf, 1, strlen(buf), f2);
	    }
	  }
	}
	FFLUSH(f2);
      }
    }
  }
  ArrayFreeLock(&sp, NULL);
  ArrayFreeLock(&spx, NULL);
  ArrayFreeLock(&linetype, NULL);

  FCLOSE(f1);
  FCLOSE(f2);

  if (_starknp > 0) {
    free(wdi);
    free(wir);
    free(wrf);
    free(wid);
  }
  return 0;
}
    
int PlotSpec(char *ifn, char *ofn, int nele, int type, 
	     double emin, double emax, double de0, double smin) {
  F_HEADER fh;
  SP_HEADER h;
  SP_RECORD r;
  SP_EXTRA rx;
  DISTRIBUTION *dist;
  TFILE *f1;
  FILE *f2;
  int n, nb, i;
  int t, t0, t1, t2;
  int r0, r1;
  double e;
  int m, k, nsp;
  double *sp, *tsp, *xsp, *kernel, *wsp;
  double de10, de01;
  double a, sig, factor;
  double *lines;
  double smax, de, hc=12.3984E3;
  int swp;
  int idist;

#if USE_MPI == 2
  if (!MPIReady()) InitializeMPI(0, 1);
#endif
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("ERROR: File %s does not exist\n", ifn);
    return -1;
  }
  f2 = fopen(ofn, "a");
  if (f2 == NULL) {
    printf("ERROR: Cannot open file %s\n", ofn);
    return -1;
  }

  int zt = (int)(0.1+fh.atom);
  double mt = GetAtomicMassTable()[(int)(fh.atom)];
  double wd = 0;
  double *wdi = NULL, *wir = NULL, *wrf = NULL, *wid = NULL;
  if (_starknp > 0) {
    wdi = malloc(sizeof(double)*_starknp);
    wir = malloc(sizeof(double)*_starknp);
    wrf = malloc(sizeof(double)*_starknp*zt);
    wid = malloc(sizeof(double)*zt);
  }
  if (_starkqc > 0 && electron_density > 0) {
    int idist;
    DISTRIBUTION *dist;
    dist = GetEleDist(&idist);
    if (idist == 0) {
      PrepStarkQC(mt, electron_density*1e10, dist->params[0], 
		  &wd, wdi, wir, zt, 1, zt, wrf, wid);
    }
  }
  
  rx.sdev = 0.0;

  t2 = abs(type) / 1000000;
  if (type < 0) t2 = -1;
  t = abs(type) % 1000000;
  t1 = t / 10000;
  t0 = t % 10000;
  t0 = t0 / 100;

  double wmin = 1e31;
  double wav = 0.0;
  double sav = 0.0;
  smax =  0.0;
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadSPHeader(f1, &h, swp);
    if (n == 0) break;
    if (h.ntransitions == 0) continue;
    if (nele >= 0 && h.nele != nele) goto LOOPEND0; 
    r1 = h.type / 10000;
    r0 = h.type % 10000;
    r0 = r0/100;
    if (type != 0) {
      if (t2 == 0) {
	if (t != h.type) goto LOOPEND0;
      } else if (t2 == 1) {
	if (r1 < t1) goto LOOPEND0;
	if (h.type < 100) goto LOOPEND0;
	if (t%10000 != h.type%10000) goto LOOPEND0;
      } else {
	if (t < 100) {
	  if (h.type > 99) goto LOOPEND0;
	  if (h.type < t) goto LOOPEND0;
	} else {
	  if (r1 < t1) goto LOOPEND0;
	  if (r0 < t0) goto LOOPEND0;
	}
      }
    }
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      e = r.energy;
      a = r.strength * e;
      if (a < smax*smin) continue;
      if (a > smax) smax = a;
      double wk = CalcStarkQC(r.trate*WCOEF, wd, wdi, wir, wrf,
			      wid, h.nele, h.type);
      wk += r.rrate*WCOEF;
      wav += wk*a;
      sav += a;
      if (wmin > wk) wmin = wk;      
    }
    continue;
  LOOPEND0:
    FSEEK(f1, h.length, SEEK_CUR);
  }
  FCLOSE(f1);
  f1 = OpenFileRO(ifn, &fh, &swp);
  de = fabs(de0);
  double ade0 = de;
  if (sav > 0) {
    wav /= sav;  
    de = sqrt(de*de + wav*wav);
  }
  de01 = 0.025*de;
  de10 = 10.0*de;
  sig = de/2.35;
  factor = 1.0/(sqrt(2*PI)*sig);
  sig = 1.0/(2*sig*sig);
  kernel = (double *) malloc(sizeof(double)*512);
  e = -256.0*de01;
  for (i = 0; i < 512; i++){
    kernel[i] = factor*exp(-sig*e*e);
    e += de01;
  }
  nsp = (emax - emin)/de01;
  sp = (double *) malloc(sizeof(double)*nsp);
  wsp = (double *) malloc(sizeof(double)*nsp);
  xsp = (double *) malloc(sizeof(double)*nsp);
  tsp = (double *) malloc(sizeof(double)*nsp);
  sp[0] = 0.0;
  tsp[0] = 0.0;
  xsp[0] = emin;
  for (i = 1; i < nsp; i++) {
    tsp[i] = 0.0;
    sp[i] = 0.0;
    wsp[i] = 0.0;
    xsp[i] = xsp[i-1] + de01;
  }

  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadSPHeader(f1, &h, swp);
    if (n == 0) break;
    if (h.ntransitions == 0) continue;
    if (nele >= 0 && h.nele != nele) goto LOOPEND; 
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
    int m0 = h.ntransitions;
    m = 3*m0;
    lines = (double *) malloc(sizeof(double)*m);  
    k = 0;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      e = r.energy;
      a = r.strength * e;
      if (a < smax*smin) continue;
      e *= HARTREE_EV;
      if (de0 < 0) e = hc/e;
      lines[k++] = e;
      lines[k++] = r.strength;
      double wk = CalcStarkQC(r.trate*WCOEF, wd, wdi, wir, wrf,
			      wid, h.nele, h.type);
      wk += r.rrate*WCOEF;
      lines[k++] = wk;
    }
    m = k;
    m0 = m/3;
    qsort(lines, m0, sizeof(double)*3, CompareLine);
    k = 0;
    i = 0;
    while (k < m && i < nsp-1) {
      if (lines[k] < xsp[i]) {
	k += 3;
      } else if (lines[k] < xsp[i+1]) {
	sp[i] += lines[k+1];
	wsp[i] += lines[k+1]*lines[k+2];
	k += 3;
      } else {
	i++;
      }
    } 
    free(lines);
    ResetWidMPI();
#pragma omp parallel default(shared) private(i, m, k)
    {
      int w = 0;
      for (i = 0; i < nsp; i++) {
	if (SkipWMPI(w++)) continue;
	if (sp[i] > 0) {
	  wsp[i] /= sp[i];
	  wsp[i] = sqrt(wsp[i]*wsp[i] + ade0*ade0);
	  wsp[i] /= de;
	  int km = (int)(wsp[i]*512);
	  double fsp = 512.0/km;
	  int km2 = km/2;
	  for (m = i-km2, k = 0; k < km; k++, m++) {
	    if (m >= 0 && m < nsp) {
	      double dk = (256+(k-km2)*fsp);
	      double fk = 0.0;
	      if (dk >= 0 && dk < 511) {
		int kk = (int)dk;
		int kk1 = kk+1;
		fk = kernel[kk]*(kk1-dk) + kernel[kk1]*(dk-kk);
#pragma omp atomic
		tsp[m] += sp[i]*fk;
	      }
	    }
	  }
	}
	sp[i] = 0.0;
	wsp[i] = 0.0;
      }
    }
    continue;
  LOOPEND:
    FSEEK(f1, h.length, SEEK_CUR);
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
	ResetWidMPI();
#pragma omp parallel default(shared) private(i, e)
	{
	  int w = 0;
	  for (i = 1; i < m; i++) {
	    if (SkipWMPI(w++)) continue;
	    e = i*de01;
	    kernel[i] = dist->dist(e, dist->params);
	  }
	}
	ResetWidMPI();
#pragma omp parallel default(shared) private(i, k, r0)
	{
	  int w = 0;
	  for (i = 0; i < nsp; i++) {
	    if (SkipWMPI(w++)) continue;
	    if (tsp[i] > 0.0) {
	      for (k = 0; k < m; k++) {
		r0 = i + k;
		if (r0 >= nsp) break;
#pragma omp atomic
		sp[r0] += tsp[i]*kernel[k];
	      }
	    }
	  }
	}
      } else if (idist == 1) {
	r1 = m/2;
	ResetWidMPI();
#pragma omp parallel default(shared) private(i, e)
	{
	  int w = 0;
	  for (i = 0; i < m; i++) {
	    if (SkipWMPI(w++)) continue;
	    e = -de01*(r1+i);
	    kernel[i] = dist->dist(e, dist->params);
	    e += de01;
	  }
	}
	ResetWidMPI();
#pragma omp parallel default(shared) private(i, k, m, r0)
	{
	  int w = 0;
	  for (i = 0; i < nsp; i++) {
	    if (SkipWMPI(w++)) continue;
	    if (tsp[i] > 0.0) {
	      for (k = 0; k < m; k++) {
		r0 = i + k - r1;
		if (r0 < 0 || r0 >= nsp) continue;
#pragma omp atomic
		sp[r0] += tsp[i]*kernel[k];
	      }
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
  free(wsp);
  free(tsp);
  free(kernel);

  FCLOSE(f1);
  fclose(f2);
  if (_starknp > 0) {
    free(wdi);
    free(wir);
    free(wrf);
    free(wid);
  }
  return 0;
}

int AddRate(ION *ion, ARRAY *rts, RATE *r, int m, int **irb) {
  LBLOCK *ib, *fb;
  BLK_RATE *brt, brt0;
  RATE *r0;
  int i, rbks;
  if (r->dir <= 0 && r->inv <= 0 && m <= 1) return 1;
  ib = ion->iblock[r->i];
  fb = ion->iblock[r->f];
  if (ib == NULL || fb == NULL) return 1;
  if (isnan(r->dir) || isinf(r->dir) || isnan(r->inv) || isinf(r->inv)) {
    MPrintf(-1, "NAN/INF rates: %s %d %d %d %g %g\n",
	    rts->id, ion->nele, r->i, r->f, r->dir, r->inv);
    return 1;
  }
  if (irb == NULL) {
    for (i = 0; i < rts->dim; i++) {
      brt = (BLK_RATE *) ArrayGet(rts, i);
      if (brt->iblock == ib && brt->fblock == fb) {
	break;
      }
    }
    if (i == rts->dim) {
      brt = NULL;
    }
  } else {
    i = irb[ib->ib][fb->ib];
    if (i >= 0 && i < rts->dim) {
      brt = (BLK_RATE *) ArrayGet(rts, i);
    } else {
      brt = NULL;
    }
  }
  if (brt == NULL) {
    brt0.iblock = ib;
    brt0.fblock = fb;
    if (irb && irb[ib->ib][fb->ib] < 0) {
      irb[ib->ib][fb->ib] = rts->dim;
    }
    brt0.rates = (ARRAY *) malloc(sizeof(ARRAY));
    rbks = ib->nlevels*fb->nlevels;
    if (rbks > _rates_block) rbks = _rates_block;
    ArrayInit(brt0.rates, sizeof(RATE), rbks);
    ArrayAppend(brt0.rates, r, NULL);
    ArrayAppend(rts, &brt0, InitBlkRateData);    
  } else {
    if (m) {
      for (i = 0; i < brt->rates->dim; i++) {
	r0 = (RATE *) ArrayGet(brt->rates, i);
	if (r0->i == r->i && r0->f == r->f) break;
      }
      if (i == brt->rates->dim) {
	if (m != 2) {
	  ArrayAppend(brt->rates, r, NULL);
	}
      } else {
	r0 = (RATE *) ArrayGet(brt->rates, i);
	if (m == 1) {
	  r0->dir += r->dir;
	  r0->inv += r->inv;
	} else if (m == 2) {
	  r0->dir *= r->dir;
	  r0->inv *= r->inv;
	} else {
	  r0->dir = r->dir;
	  r0->inv = r->inv;
	}
	return 1;
      }
    } else {
      ArrayAppend(brt->rates, r, NULL);
    }
  }
  return 0;
}

int **IdxRateBlock(int nb) {
  int **irb, i, j;
  irb = (int **) malloc(sizeof(int *)*nb);
  for (i = 0; i < nb; i++) {
    irb[i] = (int *) malloc(sizeof(int)*nb);
    for (j = 0; j < nb; j++) {
      irb[i][j] = -1;
    }
  }
  return irb;
}

void FillIdxRateBlock(int **irb, ARRAY *rts) {
  int i;
  BLK_RATE *brt;
  
  for (i = 0; i < rts->dim; i++) {
    brt = (BLK_RATE *) ArrayGet(rts, i);
    irb[brt->iblock->ib][brt->fblock->ib] = i;
  }
}

void FreeIdxRateBlock(int nb, int **irb) {
  int i;
  for (i = 0; i < nb; i++) {
    free(irb[i]);
  }
  free(irb);
}

/*
sw_mode and m0 selects which states are populated with cx.
m0%10000 specifies the n,l of the capture orbital
m = m0/10000
sw_mode=0: m = 2S+1 of the largest mixing component of the state
sw_mode=1: m = sign(L-l)*(abs(L-l)*100 + 2S+1) of the largest mixing component
sw_mode=2: m = mjl = 1 if j>l, 2 if j < l of the largest mixing component
sw_mode=3: m = sign(J-j)*(abs(J-j)*100+mjl) of any mixing component
sw_mode=4: m = 2S+1 of any mixing component of the state
sw_mode=5: m = (2J+1)*100 + (2S+1) of any mixing component
sw_mode=6: m = (2J1+1)*1000 + (2J0+1)*100 + (2S+1) of any mixing component
               such that J0 <= J <= J1
sw_mode=10: m0%10000 specifies a double capture as n2*100 + (n2-n1)
*/
int SetCXRates(int m0, char *tgt) {
  int i, k, p, ip[4], j1, j2, nn, kk, sw, jv1, jv2, m;
  int vn, vl, vl2, ix, jb, nrb, swp, nb, n, vnp, vlp;
  RATE rt, rts[NRTB];
  ION *ion;
  int **irb;
  double rcx, wt, e;
  F_HEADER fh;
  RO_HEADER h;
  RO_RECORD r[NRTB];
  CX_HEADER xh;
  CX_RECORD xr[NRTB];
  TFILE *f;
  
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  nn = -1;
  kk = -1;
  sw = 0;
  jv1 = 0;
  jv2 = 0;
  m = abs(m0);
  if (m >= 10000) {
    sw = m/10000;
    m = m%10000;
    if (sw_mode == 6 && sw > 1000) {
      jv1 = sw/100;
      sw = sw%100;
      jv2 = jv1/10;
      jv1 = jv1%10;
      if (jv2 >= 99) jv2 = 1000000;
    }
    if (m0 < 0) sw = -sw;    
  }
  if (m >= 100) {
    nn = m/100;
    kk = m%100;
    m = 1;
  } else if (m >= 10) {
    m = m%10;
    ip[3] = 1;
  } else {
    ip[3] = 0;
  }
  irb = IdxRateBlock(blocks->dim);
  KRONOS *cx;
  if (m == 1) {
    cx = KronosCX(0);
    if (kk < 0) {
      if (cx == NULL || cx->nmax <= 0) {
	printf("KRONOS CX data for H-like not setup: %d\n", m);
	return -1;
      }
      double *emass = GetAtomicMassTable();
      i = (int)ion0.atom;
      cx->pmass = emass[i-1]*AMU;
      cx->rmass = cx->pmass*cx->tmass/(cx->pmass+cx->tmass);
      rcx = 0;
      for (p = 0; p < cx->ncx; p++) {
	ip[0] = 0;
	ip[1] = 0;
	ip[2] = p;
	CXRate(&cx->rcx[p], ip, 0, 0);
	rcx += cx->rcx[p];
      }
    }
  } else if (m == 2) {
    ip[0] = 2;
    cx = KronosCX(2);
    double *emass = GetAtomicMassTable();
    i = (int)ion0.atom;
    if (i <= 0) i = 1;
    cx->pmass = emass[i-1]*AMU;
    cx->rmass = cx->pmass*cx->tmass/(cx->pmass+cx->tmass);
    cx->nmax = 0;
    cx->ilog = 3;
  }
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    int jg = ion->j[ion->ground];
    short pg = ion->p[ion->ground];
    ArrayFree(ion->cx_rates, FreeBlkRateData);
    int nbf = 0;
    if (sw_mode >= 10) {
      for (i = 0; i < ion->nlevels; i++) {
	vn = ion->vnl[i]/100;
	vl = ion->vnl[i]%100;
	vnp = ion->vni[i]/100;
	vlp = ion->vni[i]%100;
	if (vn != nn || vnp != nn-kk) continue;
	rt.i = ion->iground;
	rt.f = i;
	rt.inv = 0.0;
	rt.dir = (ion->j[i]+1.0)/((2.0*vn*vn)*(2*vnp*vnp));
	AddRate(ion, ion->cx_rates, &rt, 0, irb);
      }
      continue;
    }
    if (m == 1) {
      f = OpenFileRO(ion->dbfiles[DB_RO-1], &fh, &swp);
      if (f == NULL) {
	printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_RO-1]);
	continue;
      }
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadROHeader(f, &h, swp);
	if (h.nele != ion->nele) {
	  FSEEK(f, h.length, SEEK_CUR);
	  continue;
	}
	nrb = Min(h.ntransitions, NRTB);
	jb = 0;
	for (i = 0; i < h.ntransitions; i++) {
	  n = ReadRORecord(f, &r[jb++], swp);
	  if (jb == nrb) {
	    ResetWidMPI();
#pragma omp parallel default(shared) private(jb, j1, j2, e, p, vn, vl, ix, rcx)
	    {
	      int w = 0;
	      for (jb = 0; jb < nrb; jb++) {
		int skip = SkipWMPI(w++);
		if (skip) continue;
		rts[jb].i = r[jb].f;
		rts[jb].f = r[jb].b;
		rts[jb].dir = 0;
		rts[jb].inv = 0;
		j1 = ion->j[r[jb].f];
		j2 = ion->j[r[jb].b];
		e = ion->energy[r[jb].f] - ion->energy[r[jb].b];		
		if (e < 0) continue;
		if (sw_mode < 2 && sw != 0 && ion->sw[r[jb].b] != sw) continue;
		for (p = 0; p < r[jb].n; p++) {
		  vn = abs(r[jb].nk[p])%10000;
		  vl = vn%100;
		  vn = vn/100;
		  if (sw != 0 && sw_mode >= 2) {
		    int jw = 1;		    
		    int jv = vl*2+1;
		    if (r[jb].nk[p] < 0) {
		      jw = 2;
		      jv -= 2;
		    }
		    if (sw_mode == 3) {
		      jv = j2-jv;
		      jw += abs(jv)*100;
		      if (jv < 0) jw = -jw;
		    } else if (sw_mode == 4) {
		      jw = r[jb].nk[p]/10000;
		    } else if (sw_mode >= 5) {
		      jw = r[jb].nk[p]/10000;		      
		      jv = j2+1;
		      if (jv2 >= jv1) {
			if (jv < jv1 || jv > jv2) continue;
		      } else {
			jw += jv*100;
		      }
		    }		    
		    if (jw != sw) continue;
		  }
		  double nqr = r[jb].nq[p];
		  nqr *= nqr;
		  if (kk >= 0) {
		    if (vn != nn || vl != kk) continue;
		    rts[jb].dir += nqr/(4*vl+2.0);
		    continue;
		  }
		  if (vn > cx->nmax) continue;
		  if (m == 1) {
		    ix = cx->idn[vn-1]+vl;
		    rts[jb].dir += nqr*cx->rcx[ix]/(4*vl+2.0);
		  } else {
		    double dn = r[jb].dn[p];
		    int vn0 = (int)(vn-dn);
		    int vn1 = vn0 + 1;
		    if (vn0 == 0) vn0 = vn1;
		    if (vn1 > cx->nmax) vn1 = cx->nmax;
		    if (vn0 > cx->nmax) vn0 = cx->nmax;
		    if (vn0 == vn1) {
		      ix = cx->idn[vn-1] + vl;
		      rcx = cx->rcx[ix];
		      rts[jb].dir += nqr*rcx/(4*vl+2.0);
		    } else {
		      double dvl = vl/(double)vn;
		      int vl0 = (int)(vn0*dvl);
		      int vl1 = vl0 + 1;
		      if (vl1 >= vn0) vl1 = vn0-1;
		      int ix0 = cx->idn[vn0-1] + vl0;
		      int ix1 = cx->idn[vn0-1] + vl1;
		      double xdv = dvl*vn0 - vl0;
		      double cx0 = cx->rcx[ix0]*(1-xdv) + cx->rcx[ix1]*xdv;
		      vl0 = (int)(vn1*dvl);
		      vl1 = vl0 + 1;
		      if (vl1 >= vn1) vl1 = vn1-1;
		      ix0 = cx->idn[vn1-1] + vl0;
		      ix1 = cx->idn[vn1-1] + vl1;
		      xdv = dvl*vn1 - vl0;
		      double cx1 = cx->rcx[ix0]*(1-xdv) + cx->rcx[ix1]*xdv;
		      xdv = vn-dn - vn0;
		      rcx = cx0*(1-xdv) + cx1*xdv;
		      rts[jb].dir += nqr*rcx/(4*dvl*(vn-dn)+2.0);
		    }		    
		  }
		}
		rts[jb].dir /= (j1+1.0);
		if (ion->acx > 0) rts[jb].dir *= ion->acx;
		free(r[jb].nk);
		free(r[jb].nq);
		free(r[jb].dn);
	      }
	    }
	    for (jb = 0; jb < nrb; jb++) {
	      AddRate(ion, ion->cx_rates, &rts[jb], 0, irb);
	    }
	    nrb = h.ntransitions-i-1;
	    nrb = Min(nrb, NRTB);
	    jb = 0;
	  }
	}
      }
      continue;
    }
    if (m == 2) {
      f = OpenFileRO(ion->dbfiles[DB_CX-1], &fh, &swp);
      if (f == NULL) {
	printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_CX-1]);
	continue;
      }
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadCXHeader(f, &xh, swp);
	if (xh.nele != ion->nele ||
	    (tgt != NULL && strncmp(tgt, xh.tgts, 128))) {
	  FSEEK(f, xh.length, SEEK_CUR);
	  continue;
	}
	nbf++;
	if (xh.te0 > 0) {
	  for (i = 0; i < xh.ne0; i++) {
	    xh.e0[i] /= cx->pmass/AMU;
	  }
	}
	for (i = 0; i < xh.ne0; i++) {
	  xh.e0[i] = log(xh.e0[i]*HARTREE_EV);
	}
	nrb = Min(xh.ntransitions, NRTB);
	jb = 0;
	for (i = 0; i < xh.ntransitions; i++) {
	  n = ReadCXRecord(f, &xr[jb++], swp, &xh);
	  if (jb == nrb) {
	    ResetWidMPI();
#pragma omp parallel default(shared) private(jb, j1, j2, e, rcx)
	    {
	      int w = 0;
	      for (jb = 0; jb < nrb; jb++) {
		int skip = SkipWMPI(w++);
		if (skip) continue;
		rts[jb].i = xr[jb].f;
		rts[jb].f = xr[jb].b;
		rts[jb].dir = 0;
		rts[jb].inv = 0;
		j1 = ion->j[xr[jb].f];
		j2 = ion->j[xr[jb].b];
		e = ion->energy[xr[jb].f] - ion->energy[xr[jb].b];
		if (e < 0) continue;
		cx->nep = xh.ne0;
		cx->ep = xh.e0;
		cx->rcx = xr[jb].cx;
		int t;
		for (t = 0; t < xh.ne0; t++) {
		  if (cx->rcx[t] <= 0) {
		    cx->rcx[t] = -500.0;
		  } else {
		    cx->rcx[t] = log(cx->rcx[t]);
		  }
		}
		CXRate(&rcx, ip, xr[jb].f, xr[jb].b);
		if (ion->acx > 0) rcx *= ion->acx;
		rts[jb].dir = rcx;
		free(xr[jb].cx);
	      }
	    }
	    for (jb = 0; jb < nrb; jb++) {
	      AddRate(ion, ion->cx_rates, &rts[jb], 0, irb);
	    }
	    nrb = xh.ntransitions-i-1;
	    nrb = Min(nrb, NRTB);
	    jb = 0;
	  }
	}
	free(xh.e0);
      }
      if (nbf == 0) {
	printf("no cx data found for: %s %s\n",
	       ion->dbfiles[DB_CX-1], tgt==NULL?"NULL":tgt);
      }
      continue;
    }
    if (ion->nele < 1 || ion->nele > 2) {
      printf("CX only available for H-like and He-like ions: %d\n",
	     ion->nele);
      continue;
    }
    cx = KronosCX(ion->nele-1);
    if (cx == NULL || cx->nmax <= 0) {
      printf("KRONOS CX data not setup: %d\n", ion->nele);
      continue;
    }
    ip[0] = ion->nele-1;
    int nmax = cx->nmax;
    if (ion->icx[0] == NULL) {
      for (ix = 0; ix < 4; ix++) {
	ion->icx[ix] = malloc(sizeof(int)*cx->ncx);
	for (p = 0; p < cx->ncx; p++) {
	  ion->icx[ix][p] = -1;
	}
      }
      for (i = 0; i < ion->nlevels; i++) {
	if (ion->ibase[i] != ion->iground) {
	  continue;
	}
	vn = ion->vnl[i];
	vl = vn%100;
	vn = vn/100;
	if (vn > nmax) continue;
	p = cx->idn[vn-1]+vl;
	vl2 = 2*vl;
	if (ion->nele == 1) {
	  if (ion->j[i] > vl2) {
	    ion->icx[0][p] = i;
	  } else {
	    ion->icx[1][p] = i;
	  }
	} else {
	  if (ion->j[i] == vl2) {
	    if (ion->icx[2][p] < 0 && vl > 0) {
	      ion->icx[2][p] = i;
	    } else {
	      ion->icx[0][p] = i;
	    }
	  } else {
	    if (ion->j[i] < vl2) {
	      ion->icx[1][p] = i;
	    } else {
	      ion->icx[3][p] = i;
	    }
	  }
	}
      }
    }
    rt.inv = 0.0;
    for (p = 0; p < cx->ncx; p++) {
      ip[1] = 0;
      ip[2] = p;      
      i = ion->icx[0][p];
      if (i < 0) continue;
      rt.i = ion->iground;
      rt.f = i;
      CXRate(&rcx, ip, rt.i, rt.f);
      if (ion->acx > 0) {
	rcx *= ion->acx;
      }
      if (ion->nele == 1) {
	wt = ion->j[i]+1.0;
	vl = ion->vnl[i]%100;
	if (vl > 0) wt += wt-2;
	rt.dir = rcx*(ion->j[i]+1.0)/wt;	
	AddRate(ion, ion->cx_rates, &rt, 0, irb);
	i = ion->icx[1][p];
	if (i >= 0) {
	  rt.f = i;
	  rt.dir = rcx*(ion->j[i]+1.0)/wt;
	  AddRate(ion, ion->cx_rates, &rt, 0, irb);
	}
      } else if (ion->nele == 2) {
	rt.dir = rcx;
	AddRate(ion, ion->cx_rates, &rt, 0, irb);
	for (ix = 1; ix < 4; ix++) {
	  if (ion->icx[ix][p] >= 0) break;
	}
	if (ix == 4) continue;
	i = ion->icx[ix][p];
	wt = ion->j[i]+1.0;
	vl = ion->vnl[i]%100;
	if (vl > 0) {
	  if (ix == 1) {
	    wt = 3*wt + 6.0;
	  } else if (ix == 2) {
	    wt = 3*wt;
	  } else if (ix == 3) {
	    wt = 3*wt - 6.0;
	  }
	}
	ip[1] = 1;
	rt.i = ion->iground;
	rt.f = i;
	CXRate(&rcx, ip, rt.i, rt.f);
	if (ion->acx > 0) {
	  rcx *= ion->acx;
	}
	for (ix = 1; ix < 4; ix++) {
	  i = ion->icx[ix][p];
	  if (i < 0) continue;
	  rt.dir = rcx*(ion->j[i]+1.0)/wt;
	  rt.f = i;
	  AddRate(ion, ion->cx_rates, &rt, 0, irb);
	}
      }
    }
  }
  FreeIdxRateBlock(blocks->dim, irb);
  return 0;
}

int SetCERates(int inv) {
  int nb, i, j, ib, jb, nrb;
  int n, m, m1, k;
  int j1, j2;
  int p, q;
  ION *ion;
  RATE rt[NRTB];
  F_HEADER fh;
  CE_HEADER h;
  CE_RECORD r[NRTB];
  TFILE *f;
  double e, bte, bms;
  float *cs;
  double *data;
  double *y, *x;
  double *eusr;
  int swp;
  int **irb;

  BornFormFactorTE(&bte);
  bms = BornMass(); 
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  irb = IdxRateBlock(blocks->dim);
  for (k = 0; k < ions->dim; k++) {
    if (_krc >= 0 && k != _krc) continue;
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->ce_rates, FreeBlkRateData);
    f = OpenFileRO(ion->dbfiles[DB_CE-1], &fh, &swp);
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_CE-1]);
      continue; 
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadCEHeader(f, &h, swp);
      eusr = h.usr_egrid;
      if (h.nele == ion->nele-1) {
	if (k > 0 || ion0.nionized > 0) {
	  FSEEK(f, h.length, SEEK_CUR);
	  continue;
	}
      }
      m = h.n_usr;
      m1 = m + 1;
#pragma omp parallel default(shared) private(x, y, data, j)
      {
	data = _ce_data;
	y = data + 2;
	x = y + m1;
	data[0] = (h.te0*HARTREE_EV + bte)/bms;
	for (j = 0; j < m; j++) {
	  x[j] = log((data[0] + eusr[j]*HARTREE_EV)/data[0]);
	}
	x[m] = eusr[m-1]/(data[0]/HARTREE_EV+eusr[m-1]);
      }
      nrb = Min(NRTB, h.ntransitions);
      jb = 0;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadCERecord(f, &r[jb++], swp, &h);
	if (jb == nrb) {
	  ResetWidMPI();
#pragma omp parallel default(shared) private(ib, j1, j2, e, cs, j, data, y)
	  {
	  data = _ce_data;
	  y = data + 2;
	  int w = 0;
	  for (ib = 0; ib < nrb; ib++) {
	    int skip = SkipWMPI(w++);
	    if (skip) continue;
	    if (_ce_nmax > 0 && _ce_nmax < ion->vnl[r[ib].upper]/100) {
	      rt[ib].f = -1;
	      if (h.qk_mode == QK_FIT) free(r[ib].params);
	      free(r[ib].strength);
	      continue;
	    }	      
	    rt[ib].i = r[ib].lower;
	    rt[ib].f = r[ib].upper;
	    j1 = ion->j[r[ib].lower];
	    j2 = ion->j[r[ib].upper];
	    e = ion->energy[r[ib].upper] - ion->energy[r[ib].lower];
	    data[1] = r[ib].bethe;
	    cs = r[ib].strength;
	    y[m] = r[ib].born[0];
	    for (j = 0; j < m; j++) {
	      y[j] = cs[j];
	    }
	    CERate(&(rt[ib].dir), &(rt[ib].inv), inv, j1, j2, e, m,
		   data, rt[ib].i, rt[ib].f);
	    if (ion->ace > 0) {
	      rt[ib].dir *= ion->ace;
	      rt[ib].inv *= ion->ace;
	    }
	    //AddRate(ion, ion->ce_rates, &rt[ib], 0, irb);
	    if (h.qk_mode == QK_FIT) free(r[ib].params);
	    free(r[ib].strength);
	  }
	  }
	  for (ib = 0; ib < nrb; ib++) {
	    if (rt[ib].f < 0) continue;
	    AddRate(ion, ion->ce_rates, &rt[ib], 0, irb);
	  }
	  nrb = h.ntransitions-i-1;
	  nrb = Min(nrb, NRTB);
	  jb = 0;
	}
      }
      free(h.tegrid);
      free(h.egrid);
      free(h.usr_egrid);
    }
    FCLOSE(f);
    
    if (k == 0 && ion0.nionized > 0) {
      f = OpenFileRO(ion0.dbfiles[DB_CE-1], &fh, &swp);
      if (f == NULL) {
	printf("File %s does not exist, skipping.\n", ion0.dbfiles[DB_CE-1]);
	continue;
      }
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadCEHeader(f, &h, swp);
	eusr = h.usr_egrid;
	if (h.nele != ion0.nele) {
	  FSEEK(f, h.length, SEEK_CUR);
	  continue;
	}
	m = h.n_usr;
	m1 = m + 1;
#pragma omp parallel default(shared) private(x, y, data, j)
	{
	  data = _ce_data;
	  y = data + 2;
	  x = y + m1;
	  data[0] = (h.te0*HARTREE_EV + bte)/bms;
	  for (j = 0; j < m; j++) {
	    x[j] = log((data[0] + eusr[j]*HARTREE_EV)/data[0]);
	  }
	  x[m] = eusr[m-1]/(data[0]/HARTREE_EV+eusr[m-1]);
	}
	nrb = Min(NRTB, h.ntransitions);
	jb = 0;
	for (i = 0; i < h.ntransitions; i++) {
	  n = ReadCERecord(f, &r[jb++], swp, &h);
	  if (jb == nrb) {
	    ResetWidMPI();
#pragma omp parallel default(shared) private(ib, j1, j2, e, cs, j, p, q, data, y)
	    {	    
	    data = _ce_data;
	    y = data + 2;
	    int w = 0;
	    for (ib = 0; ib < nrb; ib++) {
	      int skip = SkipWMPI(w++);
	      if (skip) continue;
	      rt[ib].i = -1;
	      rt[ib].f = -1;
	      rt[ib].dir = 0;
	      rt[ib].inv = 0;
	      p = IonizedIndex(r[ib].lower, 0);
	      if (p < 0) {
		if (h.qk_mode == QK_FIT) free(r[ib].params);
		free(r[ib].strength);
		continue;
	      }
	      q = IonizedIndex(r[ib].upper, 0);
	      if (q < 0) {
		if (h.qk_mode == QK_FIT) free(r[ib].params);
		free(r[ib].strength);
		continue;
	      }
	      rt[ib].i = ion0.ionized_map[1][p];
	      rt[ib].f = ion0.ionized_map[1][q];
	      if (_ce_nmax > 0 && _ce_nmax < ion->vnl[rt[ib].f]/100) {
		rt[ib].f = -1;
		if (h.qk_mode == QK_FIT) free(r[ib].params);
		free(r[ib].strength);
		continue;
	      }	    
	      j1 = ion->j[rt[ib].i];
	      j2 = ion->j[rt[ib].f];
	      e = ion0.energy[q] - ion0.energy[p];
	      data[1] = r[ib].bethe;	
	      cs = r[ib].strength;
	      y[m] = r[ib].born[0];
	      for (j = 0; j < m; j++) {
		y[j] = cs[j];
	      }
	      CERate(&(rt[ib].dir), &(rt[ib].inv), inv, j1, j2, e, m,
		     data, rt[ib].i, rt[ib].f);
	      if (ion0.ace > 0) {
		rt[ib].dir *= ion0.ace;
		rt[ib].inv *= ion0.ace;
	      }
	      //AddRate(ion, ion->ce_rates, &rt, 0, irb);
	      if (h.qk_mode == QK_FIT) free(r[ib].params);
	      free(r[ib].strength);
	    }
	    }
	    for (ib = 0; ib < nrb; ib++) {
	      if (rt[ib].f < 0) continue;
	      AddRate(ion, ion->ce_rates, &rt[ib], 0, irb);
	    }
	    nrb = h.ntransitions-i-1;
	    nrb = Min(nrb, NRTB);
	    jb = 0;
	  }
	}
	free(h.tegrid);
	free(h.egrid);
	free(h.usr_egrid);
      }
      FCLOSE(f);
    }
  }

  FreeIdxRateBlock(blocks->dim, irb);
  return 0;
}

int SetTRRates(int inv) {
  int nb, i, jb, nrb;
  int n, k;
  int j1, j2;
  int p, q, m;
  ION *ion;
  RATE rt[NRTB];
  RATE rt2;
  F_HEADER fh;
  TR_HEADER h;
  TR_RECORD r[NRTB];
  TR_EXTRA rx[NRTB];
  LBLOCK *ib;
  double e, gf;
  TFILE *f;  
  int swp, iuta, im;
  int **irb, **irb2;

  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  irb = IdxRateBlock(blocks->dim);
  irb2 = IdxRateBlock(blocks->dim);
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->tr_rates, FreeBlkRateData);
    f = OpenFileRO(ion->dbfiles[DB_TR-1], &fh, &swp);
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_TR-1]);
      continue;
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadTRHeader(f, &h, swp);
      iuta = IsUTA();
      if (h.nele == ion->nele-1) {
	if (k > 0 || ion0.nionized > 0) {
	  FSEEK(f, h.length, SEEK_CUR);
	  continue;
	}
      }
      if (abs(h.multipole) <= 1) m = 0;
      else m = 1;
      nrb = Min(h.ntransitions, NRTB);
      jb = 0;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadTRRecord(f, &r[jb], &rx[jb], swp);
	jb++;
	if (jb == nrb) {
	  ResetWidMPI();
#pragma omp parallel default(shared) private(jb, j1, j2, e, gf)
	  {
	    int w = 0;
	    for (jb = 0; jb < nrb; jb++) {
	      int skip = SkipWMPI(w++);
	      if (skip) continue;
	      rt[jb].dir = 0;
	      rt[jb].inv = 0;
	      rt[jb].i = r[jb].upper;
	      if (ion0.n < 0) {
		if (ion->iblock[r[jb].upper] == NULL) continue;
		ib = ion->iblock[r[jb].upper];
		if (ib->rec &&
		    ib->rec->nrec[ib->irec] > 10) {
		  continue;
		}
	      }
	      rt[jb].f = r[jb].lower;
	      j1 = ion->j[r[jb].upper];
	      j2 = ion->j[r[jb].lower];
	      e = ion->energy[r[jb].upper] - ion->energy[r[jb].lower];
	      if (iuta) e = rx[jb].energy;
	      if (e > 0) {
		gf = OscillatorStrength(h.multipole, e,
					(double)(r[jb].strength), NULL);
		if (iuta) gf *= rx[jb].sci;
		TRRate(&(rt[jb].dir), &(rt[jb].inv), inv, j1, j2, e, (float)gf);
		if (ion->atr > 0) {
		  rt[jb].dir *= ion->atr;
		  rt[jb].inv *= ion->atr;
		}
	      }
	    }
	  }
	  for (jb = 0; jb < nrb; jb++) {
	    im = AddRate(ion, ion->tr_rates, &rt[jb], m, irb);
	    if (iuta && im == 0) {
	      rt[jb].dir = rx[jb].energy;
	      rt[jb].inv = rx[jb].sdev;
	      AddRate(ion, ion->tr_sdev, &rt[jb], 0, irb);
	    }
	  }
	  nrb = h.ntransitions-i-1;
	  nrb = Min(nrb, NRTB);
	  jb = 0;
	}
      }
    }
    FCLOSE(f);
    if (ion->nele == 1) {
      ArrayFree(ion->tr2_rates, FreeBlkRateData);
      rt2.f = FindLevelByName(ion->dbfiles[DB_EN-1], 1, 
			     "1*1", "1s1", "1s+1(1)1");
      rt2.i = FindLevelByName(ion->dbfiles[DB_EN-1], 1,
			     "2*1", "2s1", "2s+1(1)1");
      if (rt2.i >= 0 && rt2.f >= 0) {
	rt2.dir = TwoPhotonRate(ion0.atom, 0);
	rt2.inv = 0.0;
	if (ion->atr > 0) {
	  rt2.dir *= ion->atr;
	  rt2.inv *= ion->atr;
	}
	AddRate(ion, ion->tr2_rates, &rt2, 0, irb2);
      }
    } else if (ion->nele == 2) {
      ArrayFree(ion->tr2_rates, FreeBlkRateData);
      rt2.f = FindLevelByName(ion->dbfiles[DB_EN-1], 2, 
			     "1*2", "1s2", "1s+2(0)0");
      rt2.i = FindLevelByName(ion->dbfiles[DB_EN-1], 2,
			     "1*1.2*1", "1s1.2s1", "1s+1(1)1.2s+1(1)0");
      if (rt2.i >= 0 && rt2.f >= 0) {
	rt2.dir = TwoPhotonRate(ion0.atom, 1);
	rt2.inv = 0.0;
	if (ion->atr > 0) {
	  rt2.dir *= ion->atr;
	  rt2.inv *= ion->atr;
	}
	AddRate(ion, ion->tr2_rates, &rt2, 0, irb2);
      }
      if (k == 0 && ion0.nionized > 0.0) {
	rt2.f = FindLevelByName(ion->dbfiles[DB_EN-1], 1, 
			       "1*1", "1s1", "1s+1(1)1");
	rt2.i = FindLevelByName(ion->dbfiles[DB_EN-1], 1,
			       "2*1", "2s1", "2s+1(1)1");
	if (rt2.i >= 0 && rt2.f >= 0) {
	  rt2.dir = TwoPhotonRate(ion0.atom, 0);
	  rt2.inv = 0.0;
	  if (ion->atr > 0) {
	    rt2.dir *= ion->atr;
	    rt2.inv *= ion->atr;
	  }
	  AddRate(ion, ion->tr2_rates, &rt2, 0, irb2);
	}
      }
    } else if (ion->nele == 4) {
      ArrayFree(ion->tr2_rates, FreeBlkRateData);
      rt2.f = FindLevelByName(ion->dbfiles[DB_EN-1], 4,
			     "1*2.2*2", "2s2", "2s+2(0)0");
      rt2.i = FindLevelByName(ion->dbfiles[DB_EN-1], 4, 
			     "1*2.2*2", "2s1.2p1", "2s+1(1)1.2p-1(1)0");
      if (rt2.f < 0) {
	rt2.f = FindLevelByName(ion->dbfiles[DB_EN-1], 4,
				"1*2.2*2", "1s2.2s2", "1s+2(0)0.2s+2(0)0");
      }
      if (rt2.i < 0) {
	rt2.i = FindLevelByName(ion->dbfiles[DB_EN-1], 4, 
				"1*2.2*2", "1s2.2s1.2p1", "1s+2(0)0.2s+1(1)1.2p-1(1)0");
      }
      if (rt2.i >= 0 && rt2.f >= 0) {
	rt2.dir = TwoPhotonRate(ion0.atom, 2);
	rt2.inv = 0.0;
	if (ion->atr > 0) {
	  rt2.dir *= ion->atr;
	  rt2.inv *= ion->atr;
	}
	AddRate(ion, ion->tr2_rates, &rt2, 0, irb2);
      }
    } else if (ion->nele == 5) {      
      if (k == 0 && ion0.nionized > 0.0) {
	rt2.f = FindLevelByName(ion->dbfiles[DB_EN-1], 4,
				"1*2.2*2", "2s2", "2s+2(0)0");
	rt2.i = FindLevelByName(ion->dbfiles[DB_EN-1], 4, 
				"1*2.2*2", "2s1.2p1", "2s+1(1)1.2p-1(1)0");
	if (rt2.f < 0) {
	  rt2.f = FindLevelByName(ion->dbfiles[DB_EN-1], 4,
				  "1*2.2*2", "1s2.2s2", "1s+2(0)0.2s+2(0)0");
	}
	if (rt2.i < 0) {
	  rt2.i = FindLevelByName(ion->dbfiles[DB_EN-1], 4, 
				  "1*2.2*2", "1s2.2s1.2p1", "1s+2(0)0.2s+1(1)1.2p-1(1)0");
	}
	if (rt2.i >= 0 && rt2.f >= 0) {
	  rt2.dir = TwoPhotonRate(ion0.atom, 2);
	  rt2.inv = 0.0;
	  if (ion->atr > 0) {
	    rt2.dir *= ion->atr;
	    rt2.inv *= ion->atr;
	  }
	  AddRate(ion, ion->tr2_rates, &rt2, 0, irb2);
	}
      }
    }
    if (ion0.n < 0) continue;
    ExtrapolateTR(ion, inv, irb);
    if (k == 0 && ion0.nionized > 0) {
      f = OpenFileRO(ion0.dbfiles[DB_TR-1], &fh, &swp);
      if (f == NULL) {
	printf("File %s does not exist, skipping.\n", ion0.dbfiles[DB_TR-1]);
	continue;
      }
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadTRHeader(f, &h, swp);
	iuta = IsUTA();
	if (h.nele != ion0.nele) {
	  FSEEK(f, h.length, SEEK_CUR);
	  continue;
	}  
	if (abs(h.multipole) == 1) m = 0;
	else m = 1;
	nrb = Min(h.ntransitions, NRTB);
	jb = 0;
	for (i = 0; i < h.ntransitions; i++) {
	  n = ReadTRRecord(f, &r[jb], &rx[jb], swp);
	  jb++;	
	  if (jb == nrb) {
	    ResetWidMPI();
#pragma omp parallel default(shared) private(jb, j1, j2, e, gf, p, q)
	    {
	      int w = 0;
	      for (jb = 0; jb < nrb; jb++) {
		int skip = SkipWMPI(w++);
		if (skip) continue;
		rt[jb].dir = 0;
		rt[jb].inv = 0;
		p = IonizedIndex(r[jb].lower, 0);	      
		if (p < 0) {
		  continue;
		}
		q = IonizedIndex(r[jb].upper, 0);
		if (q < 0) {
		  continue;
		}
		rt[jb].i = ion0.ionized_map[1][q];
		rt[jb].f = ion0.ionized_map[1][p];
		j1 = ion->j[rt[jb].i];
		j2 = ion->j[rt[jb].f];
		e = ion0.energy[q] - ion0.energy[p];
		if (iuta) e = rx[jb].energy;	    
		if (e > 0) {
		  gf = OscillatorStrength(h.multipole, e,
					  (double)(r[jb].strength), NULL);
		  if (iuta) gf *= rx[jb].sci;
		  TRRate(&(rt[jb].dir), &(rt[jb].inv), inv, j1, j2, e, (float)gf);
		  if (ion0.atr > 0) {
		    rt[jb].dir *= ion0.atr;
		    rt[jb].inv *= ion0.atr;
		  }
		}
	      }
	    }
	    for (jb = 0; jb < nrb; jb++) {
	      im = AddRate(ion, ion->tr_rates, &rt[jb], m, irb);
	      if (iuta && im == 0) {
		rt[jb].dir = rx[jb].energy;
		rt[jb].inv = rx[jb].sdev;
		AddRate(ion, ion->tr_sdev, &rt[jb], 0, irb);
	      }
	    }
	    nrb = h.ntransitions-i-1;
	    nrb = Min(nrb, NRTB);
	    jb = 0;
	  }
	}
      }
      FCLOSE(f);
    }
  }
  FreeIdxRateBlock(blocks->dim, irb);
  FreeIdxRateBlock(blocks->dim, irb2);
  return 0;
}

int SetCIRates(int inv) { 
  int nb, i, t, jb, nrb;
  int n, m, k;
  int j1, j2;
  ION *ion;
  RATE rt[NRTB];
  F_HEADER fh;
  CI_HEADER h;
  CI_RECORD r[NRTB];
  double e;
  TFILE *f;  
  int swp;
  int **irb;
  
  if (ion0.n < 0.0) return 0;

  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }

  irb = IdxRateBlock(blocks->dim);
  for (k = 0; k < ions->dim; k++) {
    if (_krc >= 0 && k != _krc) continue;
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->ci_rates, FreeBlkRateData);
    f = OpenFileRO(ion->dbfiles[DB_CI-1], &fh, &swp);
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_CI-1]);
      continue;
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadCIHeader(f, &h, swp);
      m = h.nparams;
      if (h.nele != ion->nele) {
	FSEEK(f, h.length, SEEK_CUR);
	free(h.tegrid);
	free(h.egrid);
	free(h.usr_egrid);
	continue;
      }
      nrb = Min(h.ntransitions, NRTB);
      jb = 0;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadCIRecord(f, &r[jb++], swp, &h);
	if (jb == nrb) {	  
	  ResetWidMPI();
#pragma omp parallel default(shared) private(jb, j1, j2, e)
	  {
	    int w = 0;
	    for (jb = 0; jb < nrb; jb++) {
	      int skip = SkipWMPI(w++);
	      if (skip) continue;
	      rt[jb].i = r[jb].b;
	      rt[jb].f = r[jb].f;
	      rt[jb].dir = 0;
	      rt[jb].inv = 0;
	      j1 = ion->j[r[jb].b];
	      j2 = ion->j[r[jb].f];
	      e = ion->energy[r[jb].f] - ion->energy[r[jb].b];
	      CIRate(&(rt[jb].dir), &(rt[jb].inv), inv, j1, j2, e, m,
		     r[jb].params, rt[jb].i, rt[jb].f);
	      if (ion->aci > 0) {
		rt[jb].dir *= ion->aci;
		rt[jb].inv *= ion->aci;
	      }
	      //AddRate(ion, ion->ci_rates, &rt, 0, irb);
	      free(r[jb].params);
	      free(r[jb].strength);
	    }
	  }
	  for (jb = 0; jb < nrb; jb++) {
	    AddRate(ion, ion->ci_rates, &rt[jb], 0, irb);
	  }
	  nrb = h.ntransitions-i-1;
	  nrb = Min(nrb, NRTB);
	  jb = 0;
	}
      }
      free(h.tegrid);
      free(h.egrid);
      free(h.usr_egrid);
    }
    FCLOSE(f);
  }
  FreeIdxRateBlock(blocks->dim, irb);
  return 0;
}

int SetRRRates(int inv) { 
  int nb, i, j, jb, nrb;
  int n, m, k;
  int j1, j2;
  ION *ion;
  RATE rt[NRTB];
  F_HEADER fh;
  RR_HEADER h;
  RR_RECORD r[NRTB];
  double e;
  TFILE *f;  
  int swp;
  float *cs;
  double *data;
  double *eusr;
  double *x, *logx, *y, *p;
  int **irb;

  if (ion0.n < 0.0) return 0;
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  irb = IdxRateBlock(blocks->dim);
  for (k = 0; k < ions->dim; k++) {
    if (_krc >= 0 && k != _krc) continue;
    ion = (ION *) ArrayGet(ions, k);
    ArrayFree(ion->rr_rates, FreeBlkRateData);
    f = OpenFileRO(ion->dbfiles[DB_RR-1], &fh, &swp);
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_RR-1]);
      continue;
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadRRHeader(f, &h, swp);
      if (h.nparams <= 0) {
	printf("RR QkMode in %s must be in QK_FIT, nb=%d\n", 
	       ion->dbfiles[DB_RR-1], nb);
	exit(1);
      }
      if (h.nele != ion->nele) {
	FSEEK(f, h.length, SEEK_CUR);
	free(h.tegrid);
	free(h.egrid);
	free(h.usr_egrid);
	continue;
      }
      eusr = h.usr_egrid;
      m = h.n_usr; 
      nrb = Min(h.ntransitions, NRTB);
      jb = 0;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadRRRecord(f, &r[jb++], swp, &h);
	if (jb == nrb) {
	  ResetWidMPI();
#pragma omp parallel default(shared) private(jb, j1, j2, j, e, x, y, logx, p, data, cs)
	  {
	    data = _rr_data;
	    y = data + 1;
	    x = y + m;
	    logx = x + m;
	    p = logx + m;
	    int w = 0;
	    for (jb = 0; jb < nrb; jb++) {
	      int skip = SkipWMPI(w++);
	      if (skip) continue;
	      rt[jb].i = r[jb].f;
	      rt[jb].f = r[jb].b;
	      rt[jb].dir = 0;
	      rt[jb].inv = 0;
	      j1 = ion->j[r[jb].f];
	      j2 = ion->j[r[jb].b];
	      e = ion->energy[r[jb].f] - ion->energy[r[jb].b];
	      data[0] = 3.5 + r[jb].kl;
	      if (e < 0.0) {
		MPrintf(-1, "%d %d %10.3E %10.3E\n",
			r[jb].f, r[jb].b,
			ion->energy[r[jb].f], ion->energy[r[jb].b]);
		Abort(1);
	      }
	      cs = r[jb].strength;
	      for (j = 0; j < m; j++) {
		x[j] = (e+eusr[j])/e;
		logx[j] = log(x[j]);
		y[j] = log(cs[j]);
	      }
	      for (j = 0; j < h.nparams; j++) {
		p[j] = r[jb].params[j];
	      }
	      p[h.nparams-1] *= HARTREE_EV;
	      RRRate(&(rt[jb].dir), &(rt[jb].inv), inv, j1, j2, e, m, data,
		     rt[jb].i, rt[jb].f);
	      if (ion->arr > 0) {
		rt[jb].dir *= ion->arr;
		rt[jb].inv *= ion->arr;
	      }	    
	      free(r[jb].params);
	      free(r[jb].strength);
	    }
	  }
	  for (jb = 0; jb < nrb; jb++) {
	    AddRate(ion, ion->rr_rates, &rt[jb], 0, irb);
	  }
	  nrb = h.ntransitions-i-1;
	  nrb = Min(nrb, NRTB);
	  jb = 0;
	}
      }
      free(h.tegrid);
      free(h.egrid);
      free(h.usr_egrid);
    }
    FCLOSE(f);
    ExtrapolateRR(ion, inv, irb);
  }

  FreeIdxRateBlock(blocks->dim, irb);
  return 0;
}

int SetAIRatesInner(char *fn) {
  int nb, n, k, i, b0, nm;
  ION *ion;
  F_HEADER fh;
  AI_HEADER h;
  AI_RECORD r;
  TFILE *f;  
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
  
  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) {
    printf("File %s does not exist, skipping.\n", fn);
    return 0;
  }

  b0 = -1;
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadAIHeader(f, &h, swp);
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadAIRecord(f, &r, swp);
      if (b0 < 0 || r.b < b0) b0 = r.b;
    }
    free(h.egrid);
  }

  FSEEK(f, 0, SEEK_SET);
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
      FSEEK(f, h.length, SEEK_CUR);
    }
    free(h.egrid);
  }

  FCLOSE(f);

  return 0;
}
  
int SetAIRates(int inv) {
  int nb, i, ib, jb, nrb;
  int n, k;
  int j1, j2;
  ION *ion, *ion1;
  RATE rt[NRTB];
  F_HEADER fh;
  AI_HEADER h;
  AI_RECORD r[NRTB];
  double e;
  TFILE *f;  
  int swp;
  int ibase;
  int **irb;

  if (ion0.n < 0.0) return 0;

  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }
  irb = IdxRateBlock(blocks->dim);
  for (k = 0; k < ions->dim; k++) {
    if (_krc >= 0 && k != _krc) continue;
    if (k == 0) ion = (ION *) ArrayGet(ions, k);
    else ion = ion1;
    if (k < ions->dim - 1) ion1 = (ION *) ArrayGet(ions, k+1);
    else ion1 = NULL;
    ArrayFree(ion->ai_rates, FreeBlkRateData);
    f = OpenFileRO(ion->dbfiles[DB_AI-1], &fh, &swp);
    if (f == NULL) {
      printf("File %s does not exist, skipping.\n", ion->dbfiles[DB_AI-1]);
      continue;
    }
    for (nb = 0; nb < fh.nblocks; nb++) {
      n = ReadAIHeader(f, &h, swp);
      if (h.nele != ion->nele) {
	FSEEK(f, h.length, SEEK_CUR);
	free(h.egrid);
	continue;
      }
      nrb = Min(h.ntransitions, NRTB);
      jb = 0;
      for (i = 0; i < h.ntransitions; i++) {
	n = ReadAIRecord(f, &r[jb++], swp);
	if (jb == nrb) {
	  ResetWidMPI();
#pragma omp parallel default(shared) private(jb, j1, j2, e, ibase, ib)
	  {
	    int w = 0;
	    for (jb = 0; jb < nrb; jb++) {
	      int skip = SkipWMPI(w++);
	      if (skip) continue;
	      rt[jb].dir = 0;
	      rt[jb].inv = 0;
	      if (inner_auger == 1) {
		if (r[jb].b <= ion->KLN_max && 
		    r[jb].b >= ion->KLN_min &&
		    r[jb].f <= ion->KLN_amax &&
		    r[jb].f >= ion->KLN_amin) {
		  ibase = ion->ibase[r[jb].b] - ion->KLN_bmin;
		  if (ibase >= 0) {
#pragma omp atomic
		    ion->KLN_ai[ibase] += r[jb].rate*RATE_AU;
		  }
		}
	      } else if (inner_auger == 3) {
		if (h.nele == ion->nele-1 &&
		    r[jb].b <= ion->KLN_bmax && 
		    r[jb].b >= ion->KLN_bmin) {
		  ibase = r[jb].b - ion->KLN_bmin;
#pragma omp atomic
		  ion->KLN_ai[ibase] += r[jb].rate*RATE_AU;
		  continue;
		}
	      } else if (inner_auger == 4) {
		if (ion->iblock[r[jb].b] && ion->iblock[r[jb].b]->ionized) {
		  ib = IonIndex(ion1, ion->iblock[r[jb].b]->ib,
				ion->ilev[r[jb].b]);
		  if (ib <= ion1->KLN_bmax && ib >= ion1->KLN_bmin) {
		    ibase = ib - ion1->KLN_bmin;
#pragma omp atomic
		    ion1->KLN_ai[ibase] += r[jb].rate*RATE_AU;
		  }
		}
	      }
	      rt[jb].i = r[jb].b;
	      rt[jb].f = r[jb].f;
	      j1 = ion->j[r[jb].b];
	      j2 = ion->j[r[jb].f];
	      e = ion->energy[r[jb].b] - ion->energy[r[jb].f];
	      if (e < 0 && ion->ibase[r[jb].b] != r[jb].f) e -= ai_emin;
	      if (e > EPS16) {
		AIRate(&(rt[jb].dir), &(rt[jb].inv), inv,
		       j1, j2, e, r[jb].rate);
		if (ion->aai > 0) {
		  rt[jb].dir *= ion->aai;
		  rt[jb].inv *= ion->aai;
		}
	      }
	    }
	  }
	  for (jb = 0; jb < nrb; jb++) {
	    AddRate(ion, ion->ai_rates, &rt[jb], 0, irb);
	  }
	  nrb = h.ntransitions-i-1;
	  nrb = Min(nrb, NRTB);
	  jb = 0;
	}
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
    FCLOSE(f);
    ExtrapolateAI(ion, inv, irb);
    
    if (inner_auger == 4 && k == 0 && ion0.nionized > 0) {
      f = OpenFileRO(ion0.dbfiles[DB_AI-1], &fh, &swp);
      if (f == NULL) {
	printf("File %s does not exist, skipping.\n", ion0.dbfiles[DB_AI-1]);
	continue;
      }
      for (nb = 0; nb < fh.nblocks; nb++) {
	n = ReadAIHeader(f, &h, swp);
	for (i = 0; i < h.ntransitions; i++) {
	  n = ReadAIRecord(f, &r[0], swp);
	  ib = IonizedIndex(r[0].b, 0);
	  if (ib >= 0) {
	    ib = ion0.ionized_map[1][ib];
	    if (ib <= ion->KLN_bmax && ib >= ion->KLN_bmin) {
	      ibase = ib - ion->KLN_bmin;
	      ion->KLN_ai[ibase] += r[0].rate*RATE_AU;
	    }
	  }
	}
      }
      FCLOSE(f);
    }
  }

  FreeIdxRateBlock(blocks->dim, irb);
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
      for (t = 0; t < ion->tr2_rates->dim; t++) {
	brts = (BLK_RATE *) ArrayGet(ion->tr2_rates, t);
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

/*
 * mode = 0, total DR strength.
 * mode = 1, DR satellites.
 * mode = 2, Resonance Excitation.
 * mode = -1, total radiation branching for each level.
 * mode = -2, total autoionization branching for each level.
 */
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
  TFILE *f;
  
  if (ion0.atom <= 0) {
    printf("ERROR: Blocks not set, exitting\n");
    exit(1);
  }

  fhdr.type = DB_DR;
  fhdr.atom = ion0.atom;
  strcpy(fhdr.symbol, ion0.symbol);
  f = OpenFile(fn, &fhdr);

  if (mode >= 0) {
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
  } else {
    ilev0 = -1;
  }
  /*printf("%d %d\n", mode, ilev0);*/
  hdr.ilev = ilev0;
  hdr.nele = nele;
  n = -1;
  for (k = 0; k < ions->dim; k++) {
    ion = (ION *) ArrayGet(ions, k);
    if (ion->nele - 1 == nele) {
      if (mode < 0) {
	hdr.energy = ion->energy[ion->iground];
	hdr.j = ion->j[ion->iground];
	hdr.ilev = ion->iground;
	for (t = 0; t < ion->nlevels; t++) {
	  blk1 = ion->iblock[t];
	  if (blk1 == NULL) continue;
	  p = ion->ilev[t];
	  r1.br = blk1->r[p];
	  if (mode == -2) r1.br = 1.0 - r1.br;
	  if (!(r1.br)) continue;
	  vnl = ion->vnl[t];
	  vn = vnl/100;
	  vl = vnl - vn*100;
	  if (vn != n) {
	    if (n > 0) DeinitFile(f, &fhdr);
	    hdr.vn = vn;
	    /*printf("%d %d %d\n", hdr.ilev, hdr.j, hdr.vn);*/
	    InitFile(f, &fhdr, &hdr);
	    n = vn;
	  }
	  r1.ilev = t;
	  r1.ibase = ion->ibase[t];
	  r1.energy = ion->energy[t] - ion->energy[ion->iground];
	  r1.j = ion->j[t];
	  r1.vl = vl;
	  r1.ai = 0.0;
	  r1.total_rate = blk1->total_rate[p];
	  r1.etrans = 0.0;
	  r1.flev = -1;
	  r1.fbase = -1;
	  WriteDRRecord(f, &r1);	  
	}
      } else {
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
      }
      break;
    }
  }

  if (n >= 0) DeinitFile(f, &fhdr);  
  CloseFile(f, &fhdr);

  return 0;
}

int ModifyRates(char *fn) {
  FILE *f;
  char buf[1024];
  ION *ion;
  ARRAY *rts;
  RATE r;
  int p, n, k, m, mode;

  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return 0;
  }
  k = -1;
  while (1) {
    if (NULL == fgets(buf, 1024, f)) break;
    n = sscanf(buf, "%d %d %d", &k, &m, &mode);
    if (n == 3) break;
    k = -1;
  }
  if (k == -1) {
    printf("No line for NELE, irts, and mode\n");
    goto DONE;
  }
  for (p = 0; p < ions->dim; p++) {
    ion = (ION *) ArrayGet(ions, p);
    if (ion->nele == k) {
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
      case 7:
	rts = ion->cx_rates;
	break;
      default:
	printf("invalid mode %d\n", m);
	goto DONE;
      }
      while (1) {
	if (NULL == fgets(buf, 1024, f)) break;
	n = sscanf(buf, "%d %d %lf %lf", &(r.i), &(r.f), &(r.dir), &(r.inv));
	if (n != 4) continue;
	AddRate(ion, rts, &r, mode, NULL);
      }
      break;
    }
  }
  if (p == ions->dim) {
    printf("NELE = %d does not exist\n", k);
  }
 DONE:
  fclose(f);
  return 0;
}

int DumpRates(char *fn, int k, int m, int imax, int a) {
  TFILE *f;
  int i, t, p, q, n;
  short nele;
  double energy;
  ION *ion;
  ARRAY *rts;
  RATE *r;
  BLK_RATE *brts;
  char s[1024];
  
  if (k == 0) {
    f = FOPEN(fn, "w");
    if (f == NULL) {
      printf("cannot open file %s\n", fn);
      return -1;
    }
    n = blocks->dim*blocks->dim;
    ResetWidMPI();
#pragma omp parallel default(shared) private(p, q, t, s)
    {
      int w = 0;
      for (p = 0; p < blocks->dim; p++) {
	LBLOCK *bp = (LBLOCK *) ArrayGet(blocks, p);
	for (q = 0; q < blocks->dim; q++) {
	  if (SkipWMPI(w++)) continue;
	  LBLOCK *bq = (LBLOCK *) ArrayGet(blocks, q);
	  t = q*blocks->dim + p;	
	  sprintf(s, "%5d %5d %5d %12.5E %12.5E %12.5E %12.5E %12.5E\n",
		  p, q, bq->nlevels, bmatrix[t], bmatrix[n+p], bp->nb, bq->nb, 
		  bmatrix[t]*bq->nb);
	  FWRITE(s, 1, strlen(s), f);
	}
      }
    }
    FCLOSE(f);
    return 0;
  }
  q = 0;
  for (p = 0; p < ions->dim; p++) {
    ion = (ION *) ArrayGet(ions, p);
    if (k > 0 && ion->nele != k) continue;
    q++;
  }
  if (q == 0) return 0;
  f = FOPEN(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  for (p = 0; p < ions->dim; p++) {
    ion = (ION *) ArrayGet(ions, p);
    if (k > 0 && ion->nele != k) continue;
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
      case 7:
	rts = ion->cx_rates;
	break;
      default:
	printf("invalid mode %d\n", m);
	FCLOSE(f);
	return -1;
      }
#pragma omp parallel default(shared) private(t, q, brts, r, s)
      {
	int w = 0;
	for (t = 0; t < rts->dim; t++) {
	  brts = (BLK_RATE *) ArrayGet(rts, t);
	  for (q = 0; q < brts->rates->dim; q++) {
	    if (SkipWMPI(w++)) continue;
	    r = (RATE *) ArrayGet(brts->rates, q);
	    nele = ion->nele;
	    if (imax < 0 || (r->i <= imax && r->f <= imax)) {
	      if (a == 0) {
		char *sp = s;
		memcpy(sp, &nele, sizeof(short));
		sp += sizeof(short);
		memcpy(sp, &r->i, sizeof(int));
		sp += sizeof(int);
		memcpy(sp, &r->f, sizeof(int));
		sp += sizeof(int);
		memcpy(sp, &r->dir, sizeof(double));
		sp += sizeof(double);
		memcpy(sp, &r->inv, sizeof(double));
		sp += sizeof(double);
		FWRITE(s, 1, sp-s, f);
	      } else {
		sprintf(s, "%3d %7d %7d %10.3E %10.3E\n", 
			nele, r->i, r->f, r->dir, r->inv);
		FWRITE(s, 1, strlen(s), f);
	      }
	    }
	  }
	}
      }
      FFLUSH(f);
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
	    FWRITE(&(nele), sizeof(short), 1, f);
	    FWRITE(&t, sizeof(int), 1, f);
	    FWRITE(&(ion->iblock[t]->ib), sizeof(int), 1, f);
	    FWRITE(&q, sizeof(int), 1, f);
	    FWRITE(&(ion->j[t]), sizeof(short), 1, f);
	    FWRITE(&(ion->ibase[t]), sizeof(short), 1, f);
	    FWRITE(&(ion->vnl[t]), sizeof(short), 1, f);
	    FWRITE(&(ion->energy[t]), sizeof(double), 1, f);
	    FWRITE(&(ion->iblock[t]->n[q]), sizeof(double), 1, f);
	    FWRITE(&(ion->iblock[t]->total_rate[q]), sizeof(double), 1, f);
	    FWRITE(&(ion->iblock[t]->r[q]), sizeof(double), 1, f);
	    double x = ion->iblock[t]->rc0[q];
	    FWRITE(&x, sizeof(double), 1, f);
	    x = ion->iblock[t]->rc1[q];
	    FWRITE(&x, sizeof(double), 1, f);	    
	  } else {
	    sprintf(s, "%3d %6d %6d %6d %2d %4d %4d %15.8E %10.4E %10.4E %10.4E %10.4E %10.4E\n", 
		    nele, t, ion->iblock[t]->ib, q, ion->j[t],
		    ion->ibase[t], ion->vnl[t], ion->energy[t],
		    ion->iblock[t]->n[q],
		    ion->iblock[t]->total_rate[q],
		    ion->iblock[t]->r[q],
		    ion->iblock[t]->rc0[q],
		    ion->iblock[t]->rc1[q]);
	    FWRITE(s, 1, strlen(s), f);
	  }
	}
      }
    }
  }
  FCLOSE(f);
  return 0;
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
  TFILE *f1, *f2, *f3;
  FILE *f;
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

  f1 = OpenFileRO(fn1, &fh1, &swp1);
  f2 = OpenFileRO(fn2, &fh2, &swp2);
  if (fn3) {
    f3 = OpenFileRO(fn3, &fh3, &swp3);
  } else {
    f3 = NULL;
  }
  
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
  FSEEK(f2, 0, SEEK_SET);
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
  if (f3) {
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
  } else {
    te1 = te+1.0;
    eint1 = eint;
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

  FCLOSE(f1);
  FCLOSE(f2);
  if (f3) {
    FCLOSE(f3);
  }
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

static double vanregemoter(double y) {
  double yi[] = {0.005, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0};
  double gi[] = {1.301, 1.16, 0.977, 0.788, 0.554, 0.403, 0.29, 0.214, 0.201};
  double g;

  if (y < 0.005) {
    return 0.2757*(-0.57722 - log(y));
  } else if (y > 2.0) {
    return 0.2;
  } else {
    UVIP3P(2, 9, yi, gi, 1, &y, &g);
    return g;
  }
}

double MExpIntOne(double x) {
  double r, q2, e;

  if (x > 10) {
    return (1.-1./x+2./(x*x)-6/(x*x*x))/x;
  }
  if (x < 0.02) {
    return exp(x)*(-log(x)-0.5773+x);
  }
  if (x < 1.5) e = -0.5;
  r = (1+x)/x;
  q2 = 1./((1+x)*(1+x));
  return log(r) - (0.36+0.03*pow(x+0.01, e))*q2;
}
		   
double DRSupFactor(double z0, double d0, double t0) {
  double p0[3] = {7.966038e+00, 3.571021e-02, -1.949113e-01};
  double p[3] = {8.334987e+00, -4.561742e-02, 9.596166e-02};
  double y, n, n0, s, x0, x1, yn, a, d, t, xd, xt;

  d = exp(log(d0)-7*log(z0));
  t = exp(log(t0)-2.0*log(z0));
  xd = log(d);
  xt = log(t);

  a = pow((d/sqrt(t)),1./7.);
  y = (p[0]+p[1]*xd+p[2]*xt)/a;
  n0 = (p0[0]+p0[1]*xd+p0[2]*xt)/a;

  s = 13.605/t;
  a = s/(n0*n0);
  yn = n0*exp((log(MExpIntOne(a))-a)/7.0);
  if (yn > y) {
    x1 = n0;
    x0 = 0.5*n0;
    while (1) {
      a = s/(x0*x0);
      yn = x0*exp((log(MExpIntOne(a))-a)/7.0);
      if (yn < y) break;
      x0 *= 0.5;
    }
  } else {
    x0 = n0;
    x1 = 2.0*n0;
    while (1) {
      a = s/(x1*x1);
      yn = x1*exp((log(MExpIntOne(a))-a)/7.0);
      if (yn > y) break;
      x1 *= 2.0;
    }
  }
  while (1) {
    n0 = 0.5*(x0+x1);
    a = s/(n0*n0);
    yn = n0*exp((log(MExpIntOne(a))-a)/7.0);
    if (fabs(1-yn/y) < 1e-2) break;
    if (fabs(x1-x0) < 1e-5*n0) break;
    if (yn < y) {
      x0 = n0;
    } else {
      x1 = n0;
    }
  }
  return n0;
}

int DRSuppression(char *fn, double z, int nmax) {
  int n2, i, j, q, p, iedist, ipdist, *ipiv, info;
  double *a, *b, *c, yt, temp, trad, r, y, gy, z2, z4, ni, nj, xi, e;
  const double factor = 4.73313707E-2;
  DISTRIBUTION *edist, *pdist;
  FILE *f;

  n2 = nmax*nmax;
  if (bmatrix) free(bmatrix);
  bmatrix = (double *) malloc(sizeof(double)*n2*3);
  ipiv = (int *) malloc(sizeof(int)*nmax);
  a = bmatrix;
  b = bmatrix + n2;
  c = b + n2;
  edist = GetEleDist(&iedist);
  pdist = GetPhoDist(&ipdist);
  temp = edist->params[0];
  trad = pdist->params[0];

  for (i = 0; i < nmax; i++) {
    for (j = 0; j < nmax; j++) {
      q = i*nmax + j;
      a[q] = 0.0;
      b[q] = 0.0;
    }
  }
  z2 = z*z;
  z4 = z2*z2; 
  yt = z2*HARTREE_EV/(2.0*temp);
  for (i = 0; i < nmax; i++) {
    ni = i + 1.0;
    for (j = i+1; j < nmax; j++) {
      /* radiative transition */
      nj = j + 1.0;
      xi = 1.0/(ni*ni) - 1.0/(nj*nj); 
      r = 1.3E10*z4/(xi*ni*ni*ni)/(nj*nj*nj*nj*nj);
      q = j*nmax + i;
      p = i*nmax + j;
      a[q] += r;
      /* photo excitation */       
      if (photon_density > 0) {
	e = xi*z2*HARTREE_EV/2.0;
	y = pdist->dist(e, pdist->params)*photon_density;
	y = factor*y*r*nj*nj/(e*e*e*ni*ni);
	a[p] += y;
      }
      /* collisional excitation */
      if (electron_density > 0) {
	y = yt*xi;
	gy = vanregemoter(y);
	r = r*gy;
	r = 1.45e-6*r/(z4*z2*xi*xi*xi*sqrt(temp));      
	r = electron_density*r;
	a[q] += r;
	a[p] += r*exp(-y)*nj*nj/(ni*ni);
      }
    }
  }
  
  memcpy(c, a, sizeof(double)*n2);

  for (i = 0; i < nmax; i++) {
    ni = i + 1.0;
    q = i*nmax + i;
    for (j = 0; j < nmax; j++) {
      if (j == i) continue;
      p = i*nmax + j;
      a[q] -= a[p];
    }
    b[q] = -1.0;
    /* collisional ionization with lotz formula */
    if (i > 0 && electron_density > 0) {
      y = yt/(ni*ni);
      r = 3e4*(1.0/(temp*sqrt(temp)))*exp(-y)*FU(y)/y;
      r = electron_density*r;
      a[q] -= r;
    }
    /* photoionization with Kramer's formula */
    if (i > 0 && photon_density > 0) {
      e = z2*HARTREE_EV/(2.0*ni*ni);
      r = IntegrateRate(1, e, e, 1, &z, 0, ni, -RT_RR, PIRateKramers);
      r *= photon_density;
      p = i*nmax;
      /*printf("%f %10.3E %10.3E\n", ni, r, a[p]);*/
      a[q] -= r;
    }
  }
  
  /* ground state do not participate in the equation */
  for (i = 0; i < nmax; i++) {
    for (j = 1; j < nmax; j++) {
      q = j*nmax;
      a[q] = 0.0;
    }
    a[0] = 1E50;
  }

  /*
  for (i = 0; i < nmax; i++) {
    for (j = 0; j < nmax; j++) {
      q = i*nmax + j;
      printf("%5d %5d %10.3E %10.3E\n", i, j, a[q], b[q]);
    }
  }
  */

  DGESV(nmax, nmax, a, nmax, ipiv, b, nmax, &info);

  if (info != 0) {
    printf("Error in solving DGESV: %d\n", info);
    exit(1);
  }
  
  a[0] = 1.0;
  for (i = 1; i < nmax; i++) {
    a[i] = 0.0;
    for (j = 1; j < nmax; j++) {
      q = i*nmax + j;
      a[i] += b[q]*c[j*nmax];
      /*
      printf("%3d %3d %10.3E %10.3E\n", i, j, b[q], c[j*nmax]*b[q]);
      */
    }
  }

  if (fn) {
    f = fopen(fn, "w");
    fprintf(f, "#EDEN\t= %15.8E\n", electron_density);
    fprintf(f, "#EDIST\t= %d\n", iedist);
    fprintf(f, "#NPEDIS\t= %d\n", edist->nparams);
    for (i = 0; i < edist->nparams; i++) {    
      fprintf(f, "#\t %15.8E\n", edist->params[i]);
    }
    fprintf(f, "#PDEN\t= %15.8E\n", photon_density);
    fprintf(f, "#PDIST\t= %d\n", ipdist);
    fprintf(f, "#NPPDIS\t= %d\n", pdist->nparams);
    for (i = 0; i < pdist->nparams; i++) {
      fprintf(f, "#\t %15.8E\n", pdist->params[i]);
    }
    fprintf(f, "#Z\t = %3d\n", (int)z);
    for (i = 0; i < nmax; i++) {
      fprintf(f, "%5d %12.5E\n", i+1, a[i]);
    }
    fclose(f);
    
    /*
      for (i = 0; i < nmax; i++) {
      for (j = 0; j < nmax; j++) {
      q = i*nmax + j;
      printf("%5d %5d %10.3E\n", i, j, b[q]);
      }
      }
    */

    free(bmatrix);
    bmatrix = NULL;
  }
  free(ipiv);
  return 0;
}
  
int RydBranch(char *fn, char *ofn, int n0, int n1) {
  F_HEADER fh;
  DR_HEADER h;
  DR_RECORD r;
  TFILE *f;
  FILE *f1;
  int swp, i, j, n;
  double z, ar, br;

  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }

  if (fh.type != DB_DR) {
    printf("File type is not DB_DR\n");
    FCLOSE(f);
    return -1;
  }
   
  f1 = fopen(ofn, "w");
  
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadDRHeader(f, &h, swp);    
    if (n == 0) break;
    z = fh.atom - h.nele;
    for (j = 0; j < h.ntransitions; j++) {
      n = ReadDRRecord(f, &r, swp);
      if (n == 0) break;
      ar = 0.0;
      for (n = n0; n < h.vn; n++) {
	if (n1 >= n0 && n > n1) break;
	if (r.vl-1 < n) {
	  ar += TRRateHydrogenic(z, n, r.vl-1, h.vn, r.vl, 0);
	} 
	if (r.vl+1 < n) {
	  ar += TRRateHydrogenic(z, n, r.vl+1, h.vn, r.vl, 0);
	}
      }
      br = (r.total_rate*r.br+ar)/(r.total_rate+ar);
      fprintf(f1, "%6d %12.5E %12.5E %12.5E\n", r.ilev, ar, br, r.br);
    }
  }
  
  FCLOSE(f);
  fclose(f1);

  return 0;
}


ARRAY* _GetIons(){
/*
 * Return a reference to ions for testing purpose
 */
  return ions;
}

void SetOptionCRM(char *s, char *sp, int ip, double dp) {
  if (0 == strcmp(s, "crm:sw_mode")) {
    sw_mode = ip;
    return;
  }
  if (0 == strcmp(s, "crm:itol_nmax")) {
    _itol_nmax = ip;
    return;
  }
  if (0 == strcmp(s, "crm:sbnmax1")) {
    _sbnmax1 = ip;
    return;
  }
  if (0 == strcmp(s, "crm:sbnmax2")) {
    _sbnmax2 = ip;
    return;
  }
  if (0 == strcmp(s, "crm:silent")) {
    _silent = ip;
    return;
  }
  if (0 == strcmp(s, "crm:ce_fbr")) {
    _ce_fbr = dp;
    return;
  }
  if (0 == strcmp(s, "crm:ce_xbr")) {
    _ce_xbr = dp;
    return;
  }
  if (0 == strcmp(s, "crm:ce_nmax")) {
    _ce_nmax = ip;
    return;
  }
  if (0 == strcmp(s, "crm:epstau")) {
    _epstau = dp;
    return;
  }
  if (0 == strcmp(s, "crm:starkrw")) {
    _starkrw = dp;
    return;
  }
  if (0 == strcmp(s, "crm:starkqc")) {
    _starkqc = dp;
    return;
  }
  if (0 == strcmp(s, "crm:starkbt")) {
    _starkbt = dp;
    return;
  }
  if (0 == strcmp(s, "crm:starkaix")) {
    _starkaix = dp;
    return;
  }
  if (0 == strcmp(s, "crm:starkzix")) {
    _starkzix = dp;
    return;
  }
  if (0 == strcmp(s, "crm:starksmx")) {
    _starksmx = dp;
    return;
  }
  if (0 == strcmp(s, "crm:starkvg")) {
    _starkvg = dp;
    return;
  }
  if (0 == strcmp(s, "crm:reemit")) {
    _reemit = dp;
    return;
  }
  if (0 == strcmp(s, "crm:sp_trm")) {
    _sp_trm = ip;
    return;
  }
  if (0 == strcmp(s, "crm:rates_block")) {
    _rates_block = ip;
    return;
  }
  if (0 == strcmp(s, "crm:lblock_block")) {
    _lblock_block = ip;
    return;
  }
}

void FreeLineRec(LINEREC *r) {
  if (r->nr > 0) {
    free(r->e);
    free(r->s);
    free(r->w0);
    free(r->w);
    free(r->n0);
    free(r->n1);
    free(r->k);
    r->e = NULL;
    r->s = NULL;
    r->w0 = NULL;
    r->w = NULL;
    r->n0 = NULL;
    r->n1 = NULL;
    r->k = NULL;
  }
  r->nele = -1;
}

void PrepInterpSpec(int nd, double d0, double d1, int ds,
		    int nt, double t0, double t1, int ts,
		    double smin, double maxmem, char *fn) {
  int i, j;
  if (_interpsp.r) {
    for (i = 0; i < _interpsp.nd; i++) {
      for (j = 0; j < _interpsp.nt; j++) {
	FreeLineRec(&_interpsp.r[i][j]);
      }
      free(_interpsp.r[i]);
    }
    free(_interpsp.r);  
    if (_interpsp.xd) free(_interpsp.xd);
    if (_interpsp.xt) free(_interpsp.xt);
  }
  _interpsp.r = NULL;
  _interpsp.xd = NULL;
  _interpsp.xt = NULL;

  _interpsp.smin = smin;
  _interpsp.nd = nd;
  _interpsp.nt = nt;
  _interpsp.ds = ds;
  _interpsp.ts = ts;
  if (ds >= 0) {
    _interpsp.xd = malloc(sizeof(double)*nd);
    if (ds > 0) {
      d0 = log(d0);
      d1 = log(d1);
    }
    _interpsp.xd[0] = d0;
    if (nd > 1) {
      double dd = (d1-d0)/(nd-1);
      for (i = 1; i < nd; i++) {
	_interpsp.xd[i] = _interpsp.xd[i-1]+dd;
      }
    }
  }
  if (ts >= 0) {
    _interpsp.xt = malloc(sizeof(double)*nt);
    if (ts > 0) {
      t0 = log(t0);
      t1 = log(t1);
    }
    _interpsp.xt[0] = t0;
    if (nt > 1) {
      double dt = (t1-t0)/(nt-1);
      for (j = 1; j < nt; j++) {
	_interpsp.xt[j] = _interpsp.xt[j-1]+dt;
      }
    }
  }
  strcpy(_interpsp.fn, fn);
  _interpsp.r = malloc(sizeof(LINEREC *)*nd);
  for (i = 0; i < nd; i++) {
    _interpsp.r[i] = malloc(sizeof(LINEREC)*nt);
    for (j = 0; j < nt; j++) {
      _interpsp.r[i][j].nele = -1;
      _interpsp.r[i][j].nr = 0;
      _interpsp.r[i][j].e = NULL;
      _interpsp.r[i][j].s = NULL;
      _interpsp.r[i][j].w0 = NULL;
      _interpsp.r[i][j].w = NULL;
      _interpsp.r[i][j].n0 = NULL;
      _interpsp.r[i][j].n1 = NULL;
      _interpsp.r[i][j].k = NULL;
    }
  }
  _interpsp.tsize = 0;
  _interpsp.maxmem = maxmem*1e9;
}

void InterpSpecWF(char *fn, int nele, int type, int nmin, int nmax,
		  double c, double d, double t, double s,
		  int n, double emin, double emax) {
  FILE *f;
  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file: %s\n", fn);
    return;
  }
  double de = (emax-emin)/(n-1);
  double *x, *y;
  x = malloc(sizeof(double)*n);
  y = malloc(sizeof(double)*n);
  x[0] = emin;
  int i;
  for (i = 1; i < n; i++) {
    x[i] = x[i-1] + de;
  }
  InterpSpec(nele, type, nmin, nmax, c, d, t, s, n, x, y);
  for (i = 0; i < n; i++) {
    fprintf(f, "%17.10E %17.10E\n", x[i], y[i]);
  }
  fclose(f);
  free(x);
  free(y);
}

void InterpSpec(int nele, int type, int nmin, int nmax, double c, 
		double d, double t, double s, int n, double *x, double *y) {
  int id0, id1, it0, it1;
  double fd0, fd1, ft0, ft1;
  int i, j;
  
#if USE_MPI == 2
  if (!MPIReady()) InitializeMPI(0, 1);
#endif
  double d0 = d;
  double t0 = t;
  if (_interpsp.ds > 0) {
    d = log(d);
  }
  if (_interpsp.ts > 0) {
    t = log(t);
  }
  id0 = 0;
  id1 = 0;
  it0 = 0;
  it1 = 0;
  if (_interpsp.ds >= 0) {
    for (id1 = 0; id1 < _interpsp.nd; id1++) {
      if (d < _interpsp.xd[id1]) break;
    }  
    if (id1 == _interpsp.nd) {
      id1 = _interpsp.nd-1;
      id0 = id1;
    } else {
      id0 = id1-1;
      if (id0 < 0) id0 = 0;
    }
  }
  if (_interpsp.ts >= 0) {
    for (it1 = 0; it1 < _interpsp.nt; it1++) {
      if (t < _interpsp.xt[it1]) break;
    }
    if (it1 == _interpsp.nt) {
      it1 = _interpsp.nt-1;
      it0 = it1;
    } else {
      it0 = it1-1;
      if (it0 < 0) it0 = 0;
    }
  }
  fd0 = 1;
  ft0 = 1;
  if (id1 > id0) {
    fd0 = (_interpsp.xd[id1]-d)/(_interpsp.xd[id1]-_interpsp.xd[id0]);
    fd1 = 1-fd0;
  }
  if (it1 > it0) {
    ft0 = (_interpsp.xt[it1]-t)/(_interpsp.xt[it1]-_interpsp.xt[it0]);
    ft1 = 1-ft0;
  }
  if (fabs(fd0) < 1e-10) {
    id0 = id1;
  } else if (fabs(fd1) < 1e-10) {
    id1 = id0;
  }
  if (fabs(ft0) < 1e-10) {
    it0 = it1;
  } else if (fabs(ft1) < 1e-10) {
    it1 = it0;
  }
  for (i = 0; i < _interpsp.nd; i++) {
    for (j = 0; j < _interpsp.nt; j++) {
      _interpsp.r[i][j].ia = 0;
    }
  }
  _interpsp.r[id0][it0].ia = 1;
  _interpsp.r[id0][it1].ia = 1;
  _interpsp.r[id1][it0].ia = 1;
  _interpsp.r[id1][it1].ia = 1;
  LoadLineRec(id0, it0, nele, type, nmin, nmax);
  LoadLineRec(id0, it1, nele, type, nmin, nmax);
  LoadLineRec(id1, it0, nele, type, nmin, nmax);
  LoadLineRec(id1, it1, nele, type, nmin, nmax);
  if (_interpsp.r[id0][it0].nele < 0 ||
      _interpsp.r[id0][it1].nele < 0 ||
      _interpsp.r[id1][it0].nele < 0 ||
      _interpsp.r[id1][it1].nele < 0) {
    printf("LoadLineRec Fail: %d %d %d %d %d %d %d %d\n",
	   id0, it0, id1, it1,
	   _interpsp.r[id0][it0].nele,
	   _interpsp.r[id0][it1].nele,
	   _interpsp.r[id1][it0].nele,
	   _interpsp.r[id1][it1].nele);
    return;
  }      
  double mt = GetAtomicMassTable()[_interpsp.r[id0][it0].z];
  double dw = sqrt(t/(mt*AMU*5.11e5));
  if (id1 == id0 && it1 == it0) {
    ConvLineRec(n, x, y, mt, d0, t0, 
		s, dw, c, -1, -1, &_interpsp.r[id0][it0]);
    return;
  }
  double e0, e1, e;
  e0 = _interpsp.r[id0][it0].ae*fd0 + _interpsp.r[id1][it0].ae*fd1;
  e1 = _interpsp.r[id0][it1].ae*fd0 + _interpsp.r[id1][it1].ae*fd1;
  e = e0*ft0 + e1*ft1;
  double w0, w1, w;
  w0 = _interpsp.r[id0][it0].aw*fd0 + _interpsp.r[id1][it0].aw*fd1;
  w1 = _interpsp.r[id0][it1].aw*fd0 + _interpsp.r[id1][it1].aw*fd1;
  w = w0*ft0 + w1*ft1;
  double *y00,  *y01, *y10, *y11;
  y00 = malloc(sizeof(double)*n);
  y01 = malloc(sizeof(double)*n);
  y10 = malloc(sizeof(double)*n);
  y11 = malloc(sizeof(double)*n);
  ConvLineRec(n, x, y00, mt, d0, t0, 
	      s, dw, c, e, w, &_interpsp.r[id0][it0]);
  ConvLineRec(n, x, y01, mt, d0, t0, 
	      s, dw, c, e, w, &_interpsp.r[id0][it1]);
  ConvLineRec(n, x, y10, mt, d0, t0, 
	      s, dw, c, e, w, &_interpsp.r[id1][it0]);
  ConvLineRec(n, x, y11, mt, d0, t0, 
	      s, dw, c, e, w, &_interpsp.r[id1][it1]);
  for (i = 0; i < n; i++) {
    double d00 = log(y00[i]+1e-50);
    double d10 = log(y10[i]+1e-50);
    double d01 = log(y01[i]+1e-50);
    double d11 = log(y11[i]+1e-50);
    w0 = d00*fd0 + d10*fd1;
    w1 = d01*fd0 + d11*fd1;
    y[i] = exp(w0*ft0 + w1*ft1)-1e-50;
  }
  free(y00);
  free(y01);
  free(y10);
  free(y11);
  return;
}

double MicroFieldMode(double g, double s) {
  if (s < 1e-10) {
    double ge = 0.774 + pow(g, 0.25) + g;
    return sqrt(2/ge);
  }
  double s1 = s*exp(s);
  double bm0 = 1/(0.622 + 0.25*s1);
  double bg0 = 1/sqrt(1+(pow(g,0.25)+g)/(0.774+0.54*s1));
  return bm0*bg0;
}

double MicroFieldDist(double x, double g, double s) {
  if (s < 1e-10) {
    double x05 = sqrt(x);
    double x2 = x*x;
    double x3 = x2*x;
    double x5 = x3*x2;
    double x6 = x3*x3;
    double x35 = x3*x05;
    double x45 = x*x35;    
    double gp = g/(1+0.19*pow(g, 0.627));
    double q = 9.19 + 2.178*pow(g, 1.64);
    double z = pow((1+pow(g, 0.6)), -2.75);
    double np = exp(-gp*x05);
    double n1 = q*x3;
    double d1 = 2.25*PI*q*z+15.3*x2+1.238*q*x3+x45;
    double tn = n1*np + x6;
    double td = d1*np + x6;
    double n1d = 3*q*x2;
    double npd = -gp*np*0.5/x05;
    double x6d = 6*x5;
    double d1d = 30.6*x + 3.714*q*x2 + 4.5*x35;
    double tnd = n1d*np + n1*npd + x6d;
    double tdd = d1d*np + d1*npd + x6d;
    double r0 = tnd/td - tn*tdd/(td*td);
    if (g < 1e-10) return r0;
    double ge = 0.774 + pow(g, 0.25) + g;
    double gs = sqrt(g);
    double r1 = sqrt(2/PI)*gs*g*x2*exp(-g*x2/2);
    double a = 0.873*gs;
    return (r0 + a*r1)/(1+a);
  }
  double s2 = s*s;
  double s3 = s2*s;
  double s4 = s3*s;
  double s5 = s4*s;
  double s9 = s4*s5;
  double s016 = pow(s, 0.16);
  double s05 = sqrt(s);
  double s075 = s05*sqrt(s05);
  double s15 = s*s05;
  double s25 = s2*s05;
  double s35 = s3*s05;
  double s45 = s4*s05;
  double s55 = s5*s05;
  double s95 =  s9*s05;
  double s035 = pow(s, 0.35);
  double s18 = pow(s, 1.8);
  double s14 = s5*s9;
  double A1 = 0.59 + 2540*s4 + 3*s14;
  double A2 = (0.55 + 10*s05 + 2*s45)/(1+20*s05);
  double A3 =  2.17e-3*s5;
  double A4 = 14.8/(1+117*s35);
  double a0 = 1.15+2*s18;
  double a1 = 0.1+1.1/(1+0.145*s3);
  double a2 = 5.4/(1+20*s2) + 1.1/(1+14*s035);
  double B1 = 0.386 + 300*s2 + 1.1*s95;
  double B2 = 0.038 + 0.79*s075;
  double B3 = 3.7e-3*s55/(1+4e-3*s9);
  double b0 = (1+0.54*s25)/(1+0.07*s);
  double ga1 = 0.1 + 1.1/(1+0.174*s25);
  double ga2 = 5.4/(1+21*s15) + 1.1/(1+19*s016);
  double c = 0.097/(1+210*s25*exp(-1.3*s15));
  double g05 = sqrt(g);
  double g2 = g*g;
  double g4 = g2*g2;
  double g9 = g4*g4*g;
  double A =  A1*(1+A4*g05)/(1+A2*g2+A3*g4);
  double a = a0 + 0.5*g;
  double af = (a1 + 2*a2*g05)/(1+a2*g05);
  double B = B1/(1+B2*g2+B3*g4);
  double b = b0 + 0.25*g;
  double ga = (ga1 + 1.5*ga2*g05)/(1+ga2*g05);
  double sn = 0.0;
  if (g > 1e-10) {
    double y =  c/g9;
    double y19 = pow(y, 1.0/9.0);
    double y13 = pow(y, 1.0/3.0);
    double y23 = y13*y13;
    double y59 = y23/y19;
    double y79 = y23*y19;
    double fy = (1+0.806133*y19)/(1/240.0+0.849*y13+3.2*y59+2.43*y23+y79);
    sn = fy/(g2*g4);
  } else {
    sn = 0.806133/pow(c,0.66666667);
  }
  sn += A*exp(DLOGAM(3/af))/(af*pow(a,3/af));    
  sn += B*exp(DLOGAM(3/ga))/(ga*pow(b,3/ga));
  double x2 = x*x;
  double x4 = x2*x2;
  double x05 = sqrt(x);
  double r =  A*exp(-a*pow(x,af));
  r += B*exp(-b*pow(x,ga));
  r += exp(-g*x05)/(1+c*x4*x05);
  r *= x2/sn;
  return r;
}

double QSReduction(double g, double s) {
  double bm = MicroFieldMode(g, s);
  double b0 = -2.0;
  double b1 = 1.5;
  int nb = 1001;
  double db = (b1-b0)/(nb-1);
  double b[1001], p[1001], q[1001];
  int i;
  b[0] = b0;
  for (i = 1; i < nb; i++) {
    b[i] = b[i-1] + db;
  }
  for (i = 0; i < nb; i++) {
    b[i] = pow(10, b[i])*bm;
    p[i] = MicroFieldDist(b[i], g, s);
  }
  db = log(b[1])-log(b[0]);
  q[0] = 0.5*p[0];
  for (i = 1; i < nb; i++) {
    q[i] = q[i-1] + 0.5*(p[i-1]+p[i])*db;
  }
  double qh = 0.5*q[nb-1];
  for (i = 1; i < nb; i++) {
    if (q[i] >= qh) break;
  }
  double f = (qh-q[i-1])/(q[i]-q[i-1]);
  double xh = b[i-1]*(1-f)+b[i]*f;
  return xh/1.438;
}

double DebyeLength(double d0, double t0) {
  t0 /= HARTREE_EV;
  double t1 = sqrt(t0);
  double eta = FM1PI(d0*1.0343e-24/(t0*t1));
  double s = RBOHR*1e-8/sqrt(0.90032*t1*FM1M(eta));
  return s;
}

void ScaledSG(double s, double g, double zr, double *sn, double *gn) {
  double *a;

  if (s <= 1) a = _mfd0;
  else if (s <= 1.5) a = _mfd1;
  else if (s <= 2.0) a = _mfd2;
  else if (s <= 2.5) a = _mfd3;
  else a = _mfd4;

  double g0 = a[8]*g;
  double g1 = 1+g0;
  double zr0 = a[9]*zr;
  double zr1 = 1+zr0;
  double a0 = (a[0]+a[1]*zr0)/zr1;
  double a1 = (a[2]+a[3]*zr0)/zr1;
  *gn = g*pow(zr, (a0+a1*g0)/g1);
  a0 = (a[4]+a[5]*zr0)/zr1;
  a1 = (a[6]+a[7]*zr0)/zr1;
  *sn = s*pow(zr, (a0+a1*g0)/g1);
}

void PrepStarkQC(double mt0, double d0, double t0,
		 double *wd, double *wdi, double *wir,
		 int zt, int ne0, int ne1, double *wrf, double *wid) {
  *wd = 0;
  if (_starkqc > 0) {
    t0 /= HARTREE_EV;
    double d1 = pow(d0, ONETHIRD);
    double t1 = sqrt(t0);
    *wd = t1*2.32e-7*d1;
    int i;
    double mt;
    for (i = 0; i < _starknp; i++) {
      if (_starkmp[i] > 0 && _starkmp[i] > 0) {
	mt = _starkmp[i]*mt0/(_starkmp[i]+mt0);
	mt = sqrt(mt*AMU);
	double z1 = pow(_starkzp[i], ONETHIRD);
	wdi[i] = (*wd)/(z1*mt);
	wir[i] = _starkzp[i]*mt;
	double a = z1/(8.54e-9*d1);
	int k;
	for (k = ne0; k <= ne1; k++) {
	  double z = zt-k;
	  double g0 = _starkzp[i]*_starkzp[i]/(t0*a);
	  double eta = FM1PI(d0*1.0343e-24/(t0*t1));
	  double s0 = a*sqrt(0.90032*t1*FM1M(eta));	
	  double zix = _starkzix;
	  double g = g0, s = s0;
	  double zr = z/_starkzp[i];
	  if (fabs(zr-1)>1e-3) {
	    if (zix < 0) {
	      double zr = z/_starkzp[i];
	      ScaledSG(s0, g0, zr, &s, &g);
	    } else {
	      g = pow(zr, zix);
	    }
	  }
	  if (_starksmx >= 0 && s > _starksmx) s = _starksmx; 
	  int ik = i + (k-1)*_starknp;
	  wrf[ik] = QSReduction(g, s);
	  //printf("wrf: %d %d %g %g %g %g %g %g %g %g\n", i, k, d0, t0, a, g0, g, s0, s, wrf[ik]);
	}
      } else {
	wdi[i] = 0.0;
	wir[i] = 0.0;
	wrf[i] = 1.0;
      }
    }
    if (_starkbt > 0)  {      
      double wid0 = _starkbt*9.98e-32*(d0*d1)*HARTREE_EV*HARTREE_EV;
      for (i = ne0; i <= ne1; i++) {
	double zt2 = (zt-i+1.0);
	zt2 *= zt2;
	wid[i-1] = wid0/zt2;
      }
    } else {
      for (i = ne0; i <= ne1; i++) {
	wid[i-1] = 0.0;
      }
    }
  }
}

double CalcStarkQC(double w0, double wd, double *wdi,
		   double *wir, double *wrf, double *wid0,
		   int nele, int type) {
  double wt = w0;
  int n0, n1, dn;
  type = type%10000;
  n0 = type%100;
  n1 = type/100;
  dn = (n1-n0)%2;
  double wid = 0;
  if (wd > 0) {
    double b = 1.0;
    double r = sqrt(_starkqc*w0/wd);
    wt = wd/(_starkqc/r/r + 1/r);
    if (_starkbt > 0 && dn == 1) {
      wid = _starkqc*w0*wd/wid0[nele-1];
      b = wid/(r+wid);
      wt = 1.0/(b/wt + (1-b)/wd/1.2);
    }
    if (_starknp > 0) {
      wt = pow(wt, _starkaix);
    }
    int i;
    for (i = 0; i < _starknp; i++) {
      r = sqrt(_starkqc*w0*wir[i]/wdi[i]);
      double wti = wdi[i]/(_starkqc/r/r + 1/r);
      if (wid > 0 && dn == 1) {
	b = wid/(r+wid);
	wti = 1.0/(b/wti + (1-b)/wdi[i]/1.2);
      }
      wti *= wrf[i+(nele-1)*_starknp];
      wt += pow(wti, _starkaix)*_starkwp[i];
    }
    if (_starknp > 0) {
      wt = pow(wt, 1.0/_starkaix);
    }
  }
  return wt;
}

double NeufeldProfile(double x, double xi, double a, double t0, double ts) {
  double c0 = 0.10206;
  double c1 = 1.34308;
  double at0 = a*t0;
  double x2 = x*x;
  double xi2 = xi*xi;
  double dx = c1*fabs(x2*x - xi2*xi)/at0;
  double z = exp(-dx);
  double y;
  if (1+ts == 1) {
    y = (c0*x2/at0)*(z/(0.5*(1+z*z)));
    y *= FOUR_PI;
  } else {
    double xs = PI*ts/(2*t0);
    y = (c0*x2/at0)*cos(xs)*(z/(0.5*(1+z*z) + z*sin(xs)));
    y *= FOUR_PI*t0/(t0+ts);
  }
  return y;
}

double EscM1(double t) {
  if (t > 1e-5) {
    return (1-exp(-t))/t -1;
  }
  return -0.5*t;
}

void ConvLineRec(int n, double *x, double *y,
		 double mt, double d0, double t0,
		 double s, double dw, double c,
		 double e, double w, LINEREC *r) {
  int i, m, ny = 0;
  double de = 0;
  if (e > 0) de = e-r->ae;
  double wd;
  double *wdi = NULL, *wir = NULL, *wrf = NULL, *wid=NULL;
  int zt = r->z;
  int nele = r->nele;
  if (_starknp > 0) {
    wdi = malloc(sizeof(double)*_starknp);
    wir = malloc(sizeof(double)*_starknp);
    wrf = malloc(sizeof(double)*_starknp*zt);
    wid = malloc(sizeof(double)*zt);
  }
  PrepStarkQC(mt, d0, t0, &wd, wdi, wir, zt, nele, nele, wrf, wid);
  double rw = _starkrw;
  if (w > 0 && r->aw > 0) rw *= w/r->aw;
  double s2 = s*s;
  double sg = 1/sqrt(TWO_PI)/s;
  for (m = 0; m < n; m++) y[m] = 0.0;
  ResetWidMPI();
#pragma omp parallel default(shared) private(i, m, ny)
  {
    int wr = 0;
    for (i = 0; i < r->nr; i++) {
      if (SkipWMPI(wr++)) continue;
      double e0 = r->e[i] + de;
      double w0 = r->w[i]*rw;
      w0 = CalcStarkQC(w0, wd, wdi, wir, wrf, wid, nele, r->type);
      w0 += r->w0[i];
      //printf("wi: %d %g %g %g %g %g %g\n",  i, r->e[i], r->s[i], r->w0[i],  r->w[i], w0, wd);
      double w1 = dw*e0;
      double sw = sqrt(2*(s2 + w1*w1));
      double w2 = w0;
      double sw1 = 0.0;
      double a1 = 0.0;
      double ta = 0.0;
      double t = 0.0;
      double v0 = 0.0;
      double v = 0.0;
      double tv = 0.0;
      double ds = 0.25;
      double *yv = NULL;
      double ya = 0.0;
      double fa = 0.0;
      double y0 = 0.0;
      double x0 = 0.0;
      double ra = 0.0;
      int j;
      ny = 0;
      if (c > 0) {
	sw1 = w1*SQRT2;
	t = c*r->k[i];
	a1 = w0*0.5/sw1;
	v0 = UVoigtGauss(a1, 0.0, _starkvg);
	ta = t*v0/sw1;	
	if (ta > _epstau) {
	  if (_reemit < 0) {
	    double wg = 2.355*w1;
	    w2 = 0.5346*w0 + sqrt(0.2166*w0*w0 + wg*wg);
	    double xw = sqrt(log(ta/(log(2.0/(exp(-ta)+1.0)))));
	    xw *= (-_reemit)/0.833;
	    w2 = xw*w2 - 0.5356*w0;
	    w2 = w2*w2 - 0.2166*w0*w0;
	    if (w2 > 0) {
	      w2 = sqrt(w2)/2.355;
	      //printf("re0: %d %g %g %g %g %g %g %g %g\n", i, e0, w0, w1, w2, a1, ta, v0, xw);
	      w1 = w2;
	      sw = sqrt(2*(s2+w1*w1));
	    } else {
	      w2 = 0.0;
	    }
	  } else {
	    int dn = (int)(a1/ds);
	    dn = Max(5, dn);
	    int mn = (int)(1e3*a1/ds);
	    mn = Max(100, mn);
	    for (ny = dn; ny <= mn; ny += dn) {
	      tv = t*UVoigtGauss(a1, ny*ds, _starkvg);
	      if (tv < _epstau*sw1) break;
	    }
	    yv = malloc(sizeof(double)*ny);
	    x0 = ds;
	    ya = 0.0;
	    for (j = 0; j < ny; j++, x0 += ds) {
	      v = UVoigtGauss(a1, x0, _starkvg);
	      tv = t*v/sw1;
	      fa = EscM1(tv);
	      yv[j] = v*fa;
	      ya += yv[j];
	    }
	    y0 = v0*EscM1(ta);
	    ya = (2*ya + y0)*ds;
	    ra = 1 - _reemit*ya/(1+ya);	    
	    //printf("re1: %d %d %d %d %g %g %g %g %g %g %g %g\n", i, ny, dn, mn, e0, w0, w1, a1, y0, ta, ya, ra);
	  }
	}
      }
      double a = w0*0.5/sw;
      double b = r->s[i]/sw;
      if (ny > 0) b *= ra;
      for (m = 0; m < n; m++) {
	v = UVoigtGauss(a, (x[m]-e0)/sw, _starkvg);
#pragma omp atomic
	y[m] += b*v;
      }
      if (ny > 0) {
	x0 = ds;
	b = ds*r->s[i]*ra;
	for (j = 0; j <= ny; j++, x0 += ds) {
	  if (j < ny) {
	    ya = yv[j]*b;
	  } else {
	    ya = y0*b;
	    x0 = 0.0;
	  }
	  for (m = 0; m < n; m++) {
	    double xi = e0+x0*sw1;
	    double dx = (x[m]-xi)/s;
	    if (fabs(dx) < 10)  {
#pragma omp atomic
	      y[m] += ya*sg*exp(-0.5*dx*dx);
	    }
	    if (j < ny) {
	      xi = e0 - x0*sw1;
	      dx = (x[m]-xi)/s;
	      if (fabs(dx) < 10) {
#pragma omp atomic
		y[m] += ya*sg*exp(-0.5*dx*dx);
	      }
	    }
	  }
	}
	free(yv);
      }
    }
  }
  double dxm = 0.5*(x[1]-x[0]);
  double dxp = dxm;
  for (m = 0; m < n; m++) {
    if (m < n-1) dxp = 0.5*(x[m+1]-x[m]);
    y[m] *= (dxm+dxp);
    dxm = dxp;
  }
  if (_starknp > 0) {
    free(wdi);
    free(wir);
    free(wrf);
    free(wid);
  }
}

void LoadLineRec(int id0, int it0, int nele,
		 int type, int nmin, int nmax) {
  if (_interpsp.r[id0][it0].nele == nele &&
      (type < 0 || _interpsp.r[id0][it0].type == type) &&
      _interpsp.r[id0][it0].nmin == nmin &&
      _interpsp.r[id0][it0].nmax == nmax) return;
  int i, j;
  double ts0 = TotalSize();
  double ts1;
  if (_interpsp.maxmem > 0 && _interpsp.tsize > _interpsp.maxmem) {    
    for (i = 0; i < _interpsp.nd; i++) {
      for (j = 0; j < _interpsp.nt; j++) {
	if (!_interpsp.r[i][j].ia) {
	  FreeLineRec(&_interpsp.r[i][j]);
	}
      }
    }
    ts1 = TotalSize();
    _interpsp.tsize += ts1-ts0;
  }
  char fn[1024];
  if (_interpsp.ds < 0 && _interpsp.ts < 0) {
    strcpy(fn, _interpsp.fn);
  } else if (_interpsp.ds < 0) {
    sprintf(fn, _interpsp.fn, it0);
  } else if (_interpsp.ts < 0) {
    sprintf(fn, _interpsp.fn, id0);
  } else {
    sprintf(fn, _interpsp.fn, id0, it0);
  }
  TFILE *f1;
  F_HEADER fh;
  SP_HEADER h;
  SP_RECORD r;
  SP_EXTRA rx;
  int swp;
  
  f1 = OpenFileRO(fn, &fh, &swp);
  if (f1 == NULL) {
    printf("ERROR: File %s does not exist\n", fn);
    return;
  }
  
  LINEREC *rec = &_interpsp.r[id0][it0];
  rec->z = (int)(fh.atom);
  FreeLineRec(rec);
  rx.sdev = 0.0;
  double smax = 0.0;
  int nb, n, r0, r1, nr, imin, imax, nlev;
  int ilo0, ilo1, iup0, iup1;
  if (type < 0) {
    ilo0 = nmin%10000;
    ilo1 = nmin/10000;
    iup0 = nmax%10000;
    iup1 = nmax%10000;
    if (ilo1 == 0) ilo1 = ilo0;
    if (iup1 == 0) iup1 = iup0;
  } else {
    ilo0 = -1;
    ilo1 = -1;
    iup0 = -1;
    iup1 = -1;
  }
  nr = 0;
  nlev = 0;
  imin = -1;
  imax = -1;
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadSPHeader(f1, &h, swp);
    if (n == 0) break;
    if (h.ntransitions == 0) continue;
    if (h.nele != nele) goto LOOPEND0;
    if (h.type > 0 && type > 0) {
      r1 = h.type / 10000;
      r0 = h.type % 10000;
      if (r0 != type) goto LOOPEND0;
      if (r1 < nmin) goto LOOPEND0;
      if (r1 > nmax) goto LOOPEND0;
    }
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      if (h.type == 0) {
	if (imin < 0 || imin > r.upper) {
	  imin = r.upper;
	}
	if (imax < r.upper) {
	  imax = r.upper;
	}
      } else {
	if (type < 0) {
	  if (r.lower < ilo0 || r.lower > ilo1) continue;
	  if (r.upper < iup0 || r.upper > iup1) continue;
	}
	if (r.strength < smax*_interpsp.smin) continue;
	if (r.strength > smax) smax = r.strength;
	nr++;
      }
    }
    continue;
  LOOPEND0:
    FSEEK(f1, h.length, SEEK_CUR);
  }
  FCLOSE(f1);
  f1 = OpenFileRO(fn, &fh, &swp);
  nlev = imax-imin+1;
  rec->e = malloc(sizeof(double)*nr);
  rec->s = malloc(sizeof(double)*nr);
  rec->w0 = malloc(sizeof(double)*nr);
  rec->w = malloc(sizeof(double)*nr);
  rec->n0 = malloc(sizeof(double)*nr);
  rec->n1 = malloc(sizeof(double)*nr);
  rec->k = malloc(sizeof(double)*nr);
  double *dn = malloc(sizeof(double)*nlev);
  double *dw = malloc(sizeof(double)*nlev);
  rec->nele = nele;
  rec->type = type;
  rec->nmin = nmin;
  rec->nmax = nmax;
  rec->ae = 0.0;
  rec->aw = 0.0;
  rec->nt = 0.0;
  rec->ni = 0.0;
  rec->nr = nr;
  nr = 0;  
  for (nb = 0; nb < fh.nblocks; nb++) {
    n = ReadSPHeader(f1, &h, swp);
    if (n == 0) break;
    if (h.ntransitions == 0) continue;
    if (h.nele != nele) goto LOOPEND;
    if (h.type > 0 && type > 0) {
      r1 = h.type / 10000;
      r0 = h.type % 10000;
      if (r0 != type) goto LOOPEND;
      if (r1 < nmin) goto LOOPEND;
      if (r1 > nmax) goto LOOPEND;
    }
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadSPRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      if (h.type == 0) {
	if (h.nele == nele) {
	  rec->ni += r.strength;
	  int i1 = r.upper-imin;
	  dn[i1] = r.strength;
	  dw[i1] = r.rrate;
	}
	rec->nt += r.strength;
      } else {	
	if (type < 0) {
	  if (r.lower < ilo0 || r.lower > ilo1) continue;
	  if (r.upper < iup0 || r.upper > iup1) continue;
	}
	if (r.strength < smax*_interpsp.smin) continue;
	rec->e[nr] = r.energy;
	rec->s[nr] = r.strength;
	rec->w0[nr] = r.rrate;
	rec->w[nr] = r.trate;
	if (type < 0 && h.type > 0) {
	  rec->type = h.type%10000;
	}
	int i0 = r.lower-imin;
	int i1 = r.upper-imin;
	rec->n0[nr] = dn[i0];
	rec->n1[nr] = dn[i1];
	rec->k[nr] = 2.528e-24*(dw[i1]/dw[i0])*r.rrate*rec->n0[nr]/rec->nt;
	rec->n0[nr] /= dw[i0];
	rec->n1[nr] /= dw[i1];
	nr++;
      }
    }
    continue;
  LOOPEND:
    FSEEK(f1, h.length, SEEK_CUR);
  }
  FCLOSE(f1);
  free(dn);
  free(dw);
  if (nr < rec->nr) {
    rec->e = realloc(rec->e, sizeof(double)*nr);
    rec->s = realloc(rec->s, sizeof(double)*nr);
    rec->w0 = realloc(rec->w0, sizeof(double)*nr);
    rec->w = realloc(rec->w, sizeof(double)*nr);
    rec->n0 = realloc(rec->n0, sizeof(double)*nr);
    rec->n1 = realloc(rec->n1, sizeof(double)*nr);
    rec->k = realloc(rec->k, sizeof(double)*nr);
    rec->nr = nr;
  }

  double ts = 0.0;
  for (i = 0; i < nr; i++) {
    rec->e[i] *= HARTREE_EV;
    rec->w0[i] *= WCOEF;
    rec->w[i] *= WCOEF;
    rec->k[i] /= rec->e[i]*rec->e[i];
    rec->ae += rec->s[i]*rec->e[i];
    rec->aw += rec->s[i]*rec->w[i];
    ts += rec->s[i];
    //printf("rl: %d %d %d %g %g %g %g %g %g %g %g\n", id0, it0, i, rec->e[i], rec->w0[i], rec->w[i], rec->s[i], rec->n0[i], rec->k[i], rec->nt, rec->ni);
  }
  rec->ae /= ts;
  rec->aw /= ts;
  
  ts1 = TotalSize();
  _interpsp.tsize += ts1-ts0;
}

void SetStarkZMP(int np, double *wzm) {
  if (_starknp > 0) {
    free(_starkzp);
    free(_starkmp);
    free(_starkwp);
    _starknp = 0;
    _starkzp = NULL;
    _starkmp = NULL;
    _starkwp = NULL;
  }
  if (np <= 0) return;
  if (wzm == NULL) {
    np = 1;
  }
  _starknp = np;
  _starkzp = malloc(sizeof(double)*np);
  _starkmp = malloc(sizeof(double)*np);
  _starkwp = malloc(sizeof(double)*np);
  if (wzm == NULL) {
    _starkzp[0] = 1.0;
    _starkmp[0] = 1.0;
    _starkwp[0] = 1.0;
    return;
  }
  int i, k;
  double zt = 0.0;
  for (i = 0; i < np; i++) {
    k = 3*i;
    _starkwp[i] = wzm[k];
    _starkzp[i] = wzm[k+1];
    _starkmp[i] = wzm[k+2];
    zt += _starkwp[i];
  }
  for (i = 0; i < np; i++) {
    _starkwp[i] = _starkwp[i]/zt;
  }
}

void SortBranches(ARRAY *rts, IDXARY *idx, IDXARY *fdx, int ig, int *nr, RATE **rs) {
  int i, m, ilo, iup;
  BLK_RATE *brts;
  RATE *r;

  for (i = 0; i < idx->n; i++) {
    nr[i] = 0;
    rs[i] = NULL;
  }
  for (i = 0; i < rts->dim; i++) {
    brts = (BLK_RATE *) ArrayGet(rts, i);
    for (m = 0; m < brts->rates->dim; m++) {
      r = (RATE *) ArrayGet(brts->rates, m);
      iup = IdxGet(idx, r->i);
      ilo = IdxGet(fdx, r->f);
      if (iup >= 0 && ilo >= 0 && ilo != ig) {
	nr[iup]++;
      }
    }
  }
  for (i = 0; i < idx->n; i++) {
    if (nr[i]) {
      rs[i] = malloc(sizeof(RATE)*nr[i]);
      nr[i] = 0;
    }
  }
  for (i = 0; i < rts->dim; i++) {
    brts = (BLK_RATE *) ArrayGet(rts, i);
    for (m = 0; m < brts->rates->dim; m++) {
      r = (RATE *) ArrayGet(brts->rates, m);
      iup = IdxGet(idx, r->i);
      ilo = IdxGet(fdx, r->f);
      if (iup >= 0 && ilo >= 0 && ilo != ig) {
	rs[iup][nr[iup]].i = r->i;
	rs[iup][nr[iup]].f = r->f;
	rs[iup][nr[iup]].dir = r->dir;
	rs[iup][nr[iup]].inv = r->inv;
	nr[iup]++;
      }
    }
  }
}

void RateCoefficients(char *ofn, int k0, int k1, int nexc, int ncap0,
		      int nt, double t0, double t1,
		      int nd, double d0, double d1, int md) {
  TFILE *f;
  RC_RECORD rc;
  RC_HEADER rh;
  F_HEADER fh;
  int ms[RC_TT], ir[RC_TT];
  int *ik, *ii, *ia, *io;
  int nk, ni, na, nb, p, i, ntd, m, dj, vn, nkk, nii, nki, nbi, nr, s, n, ib;
  int it, id, ilo, iup, j0, nce, nci, nrr, ndr, nre, nea, kg, ig, n1;
  double dt, dd, rdt, rdd, *ra, *ra0, ek, ei, de, te, mp[3], br, rt, x;
  double **wr, *drs;
  int *nbai, *nbtr, ncap, mdr, mea, mdrea;
  RATE *r, **bai, **btr;
  BLK_RATE *brts;
  LBLOCK *blk;
  ION *ion;
  IDXARY iad, ibd, iid, idd;

  ncap = ncap0%100;
  mdr = ncap0/100;
  mea = mdr/10;
  mdr = mdr%10;
  mp[1] = -1.0;
  mp[2] = -1.0;
  fh.type = DB_RC;
  fh.atom = ion0.atom;
  strcpy(fh.symbol, ion0.symbol);
  f = OpenFile(ofn, &fh);
  if (f == NULL) {
    printf("cannot open file %s\n", ofn);
    return;
  }
  for (m = 1; m < RC_TT; m++) {
    ir[m] = 0;
    ms[m] = 0;
  }
  mdrea = 0;
  if (md == 0) {
    for (m = 1; m < RC_TT; m++) ms[m] = 1;
  } else {
    while (md) {
      m = md%10;
      if (m > 0 && m < RC_TT) ms[m] = 1;
      else if (m == RC_TT) {
	ms[RC_DR] = 1;
	ms[RC_EA] = 1;
	mdrea = 1;
      }
      md /= 10;
    }
  }
  ntd = nt*nd;
  rc.rc = malloc(sizeof(float)*(ntd+nt));
  dt = 0.0;
  dd = 0.0;
  if (nt > 1) dt = (log(t1)-log(t0))/(nt-1);
  if (nd > 1) dd = (log(d1)-log(d0))/(nd-1);
  rdt = exp(dt);
  rdd = exp(dd);  
  kg = 0;
  ig = 0;
  rh.nexc = nexc;
  rh.ncap = ncap;
  rh.mexc = nexc;
  rh.te0 = t0;
  rh.de0 = d0;
  rh.nte = nt;
  rh.nde = nd;
  rh.dte = dt;
  rh.dde = dd;
  
  for (p = 0; p < ions->dim; p++) {
    ion = (ION *) ArrayGet(ions, p);    
    if (ion->nele < k0 || ion->nele > k1) continue;
    nk = 0;
    ni = 0;
    na = 0;
    nb = 0;
    nbai = NULL;
    nbtr = NULL;
    bai = NULL;
    btr = NULL;
    ek = ion->energy[ion->ground];    
    if (ion->iground >= 0) {
      ei = ion->energy[ion->iground];
    } else {
      ei = 0;
    }
    j0 = ion->j[ion->ground];
    for (i = 0; i < ion->nlevels; i++) {
      vn = ion->vnl[i]/100;
      dj = abs(ion->j[i] - j0);
      if (dj%2 == 0) {
	if (vn <= nexc) {
	  if (ion->energy[i] < ei) {
	    nb++;
	  }
	  nk++;
	}
	if (vn > ncap && ion->energy[i] > ei) {
	  na++;
	}
      } else {
	if (vn <= nexc) {
	  ni++;
	}
      }
    }
    if (nk > 0) {
      ik = malloc(sizeof(int)*nk);
    } else {
      ik = NULL;
    }
    if (nb > 0) {
      io = malloc(sizeof(int)*nb);
    } else {
      io = NULL;
    }
    if (ni > 0) {
      ii = malloc(sizeof(int)*ni);
    } else {
      ii = NULL;
    }
    if (na > 0) {
      ia = malloc(sizeof(int)*na);
    } else {
      ia = NULL;
    }
    nk = 0;
    ni = 0;
    na = 0;
    nb = 0;
    n1 = 0;
    for (i = 0; i < ion->nlevels; i++) {
      vn = ion->vnl[i]/100;
      dj = abs(ion->j[i] - j0);
      if (dj%2 == 0) {
	if (vn <= nexc) {
	  if (ion->energy[i] < ei) {
	    io[nb] = i;
	    nb++;
	  }
	  ik[nk] = i;
	  nk++;
	}
	if (vn > ncap && ion->energy[i] > ei) {
	  ia[na] = i;
	  na++;
	}
	if (ion->energy[i] > ei && vn > n1) n1 = vn;
      } else {
	if (vn <= nexc) {
	  ii[ni] = i;
	  ni++;
	}
      }
    }
    if (nk > 0) {
      InitIdxAry(&iad, nk, ik);
      kg = IdxGet(&iad, ion->ground);
    }
    if (nb > 0) {
      InitIdxAry(&ibd, nb, io);
    }
    if (ni > 0) {
      InitIdxAry(&iid, ni, ii);
      ig = IdxGet(&iid, ion->iground);
    }
    if (na > 0) {      
      InitIdxAry(&idd, na, ia);
      if (ms[RC_DR]) {
	nbtr = malloc(sizeof(int)*na);
	btr = malloc(sizeof(RATE *)*na);
      }
      if (ms[RC_RE] || ms[RC_EA]) {
	nbai = malloc(sizeof(int)*na);
	bai = malloc(sizeof(RATE *)*na);
      }
    } else {
      ms[RC_DR] = 0;
      ms[RC_RE] = 0;
      ms[RC_EA] = 0;
    }
    nkk = nk*nk;
    nii = ni*ni;
    nki = nk*ni;
    nbi = nb*ni;
    nr = 0;
    if (ms[RC_CE]) {
      ir[RC_CE] = 0;
      nr += nkk;
    }
    if (ms[RC_CI]) {
      ir[RC_CI] = nr;
      nr += nki;
    }
    if (ms[RC_RR]) {
      ir[RC_RR] = nr;
      nr += nki;
    }
    if (ms[RC_DR]) {
      ir[RC_DR] += nr;
      nr += nbi;
    }
    if (ms[RC_RE]) {
      ir[RC_RE] = nr;
      nr += nii;
    }
    if (ms[RC_EA]) {
      ir[RC_EA] = nr;
      nr += nbi;
    }
    printf("irs: %d %d %d %d %d %d\n",
	   ir[RC_CE], ir[RC_CI], ir[RC_RR], ir[RC_DR], ir[RC_RE], ir[RC_EA]);
    wr = malloc(sizeof(double *)*nr);
    for (i = 0; i < nr; i++) {
      wr[i] = NULL;
    }
    drs = malloc(sizeof(double)*n1*nd);
    _krc = p;
    te = t0;
    nce = 0;
    nci = 0;
    nrr = 0;
    nre = 0;
    nea = 0;
    ndr = 0;
    for (it = 0; it < nt; it++) {
      ReinitCRM(2);
      mp[0] = te;
      printf("ion=%d te=%g nk=%d nb=%d na=%d ni=%d\n",
	     ion->nele, te, nk, nb, na, ni);
      SetEleDist(0, 3, mp);
      if (ms[RC_CE] || ms[RC_EA]) {
	SetCERates(1);
      }
      if (ms[RC_CI]) {
	SetCIRates(1);
      }
      if (ms[RC_RR]) {
	SetRRRates(0);
      }
      if (ms[RC_DR] || ms[RC_EA] || ms[RC_RE]) {
	SetEleDensity(0);
	SetAIRates(1);
	InitBlocks();
	DRBranch();
	if (it == 0) {
	  if (ms[RC_DR]) {
	    SortBranches(ion->tr_rates, &idd, &ibd, kg, nbtr, btr);
	  }
	  if (ms[RC_EA] || ms[RC_RE]) {
	    SortBranches(ion->ai_rates, &idd, &iid, ig, nbai, bai);
	  }
	}
	de = d0;
	for (id = 0; id < nd; id++) {
	  if (id == 0) {
	    for (i = 0; i < n1; i++) {
	      drs[i*nd+id] = 1.0;
	    }
	  } else {
	    SetEleDensity(de/1e10);
	    DRSuppression(NULL, (ion0.atom-ion->nele+1), n1);
	    for (i = 0; i < n1; i++) {
	      drs[i*nd+id] = bmatrix[i];
	    }
	    free(bmatrix);
	    bmatrix = NULL;
	  }
	  de *= rdd;
	}
      }
      if (ms[RC_CE]) {
	for (i = 0; i < ion->ce_rates->dim; i++) {
	  brts = (BLK_RATE *) ArrayGet(ion->ce_rates, i);
	  for (m = 0; m < brts->rates->dim; m++) {
	    r = (RATE *) ArrayGet(brts->rates, m);
	    ilo = IdxGet(&iad, r->i);
	    iup = IdxGet(&iad, r->f);
	    if (ilo >= 0 && iup >= 0) {
	      ra = wr[ir[RC_CE]+iup*nk+ilo];
	      if (ra == NULL) {
		ra = malloc(sizeof(double)*nt);
		wr[ir[RC_CE]+iup*nk+ilo] = ra;
		nce++;
	      }
	      ra[it] = r->inv;
	    }
	  }
	}
      }
      if (ms[RC_CI]) {
	for (i = 0; i < ion->ci_rates->dim; i++) {
	  brts = (BLK_RATE *) ArrayGet(ion->ci_rates, i);
	  for (m = 0; m < brts->rates->dim; m++) {
	    r = (RATE *) ArrayGet(brts->rates, m);
	    ilo = IdxGet(&iad, r->i);
	    iup = IdxGet(&iid, r->f);
	    if (ilo >= 0 && iup >= 0) {
	      ra = wr[ir[RC_CI]+iup*nk+ilo];
	      if (ra == NULL) {
		ra = malloc(sizeof(double)*nt);
		wr[ir[RC_CI]+iup*nk+ilo] = ra;
		nci++;
	      }
	      ra[it] = r->inv;
	    }
	  }
	}
      }
      if (ms[RC_RR]) {
	for (i = 0; i < ion->rr_rates->dim; i++) {
	  brts = (BLK_RATE *) ArrayGet(ion->rr_rates, i);
	  for (m = 0; m < brts->rates->dim; m++) {
	    r = (RATE *) ArrayGet(brts->rates, m);
	    iup = IdxGet(&iid, r->i);
	    ilo = IdxGet(&iad, r->f);
	    if (ilo >= 0 && iup >= 0) {
	      ra = wr[ir[RC_RR]+iup*nk+ilo];
	      if (ra == NULL) {
		ra = malloc(sizeof(double)*nt);
		wr[ir[RC_RR]+iup*nk+ilo] = ra;
		nrr++;
	      }
	      ra[it] = r->dir;
	    }
	  }
	}
      }
      if (ms[RC_DR] || ms[RC_RE]) {
	for (i = 0; i < ion->ai_rates->dim; i++) {
	  brts = (BLK_RATE *) ArrayGet(ion->ai_rates, i);
	  for (m = 0; m < brts->rates->dim; m++) {
	    r = (RATE *) ArrayGet(brts->rates, m);
	    iup = IdxGet(&iid, r->f);
	    ilo = IdxGet(&idd, r->i);
	    if (iup >= 0 && ilo >= 0) {
	      blk = ion->iblock[r->i];
	      br = blk->r[ion->ilev[r->i]];
	      rt = blk->total_rate[ion->ilev[r->i]];
	      vn = (ion->vnl[r->i]/100)-1;
	      if (ms[RC_DR]) {
		ra = wr[ir[RC_DR]+iup*nb+kg];
		if (ra == NULL) {
		  ra = malloc(sizeof(double)*ntd);
		  for (s = 0; s < ntd; s++) ra[s] = 0.0;
		  wr[ir[RC_DR]+iup*nb+kg] = ra;
		  ndr++;
		}
		de = d0;
		for (id = 0; id < nd; id++) {		  
		  ra[id*nt+it] += r->inv*br*drs[vn*nd+id];
		  de *= rdd;
		}
		if (mdr == 0) {
		  ra0 = ra;
		  for (n = 0; n < nbtr[ilo]; n++) {
		    ib = IdxGet(&ibd, btr[ilo][n].f);
		    if (ib >= 0 && ib != kg) {
		      ra = wr[ir[RC_DR]+iup*nb+ib];
		      if (ra == NULL) {
			ra = malloc(sizeof(double)*ntd);
			for (s = 0; s < ntd; s++) ra[s] = 0.0;
			wr[ir[RC_DR]+iup*nb+ib] = ra;
			ndr++;
		      }
		      br = btr[ilo][n].dir/rt;
		      de = d0;
		      for (id = 0; id < nd; id++) {
			x = r->inv*br*drs[vn*nd+id];
			ra[id*nt+it] += x;
			ra0[id*nt+it] -= x;
			de *= rdd;
		      }
		    }
		  }
		}
	      }
	      if (ms[RC_RE]) {
		for (n = 0; n < nbai[ilo]; n++) {
		  ib = IdxGet(&iid, bai[ilo][n].f);
		  if (ib >= 0 && ib > iup) {
		    ra = wr[ir[RC_RE]+ib*ni+iup];
		    if (ra == NULL) {
		      ra = malloc(sizeof(double)*nt);
		      for (s = 0; s < nt; s++) ra[s] = 0.0;
		      wr[ir[RC_RE]+ib*ni+iup] = ra;
		      nre++;
		    }
		    ra[it] += r->inv*bai[ilo][n].dir/rt;
		  }
		}
	      }
	    }
	  }
	}
      }
      if (ms[RC_EA]) {
	for (i = 0; i < ion->ce_rates->dim; i++) {
	  brts = (BLK_RATE *) ArrayGet(ion->ce_rates, i);
	  for (m = 0; m < brts->rates->dim; m++) {
	    r = (RATE *) ArrayGet(brts->rates, m);
	    ilo = IdxGet(&ibd, r->i);
	    iup = IdxGet(&idd, r->f);
	    if (ilo >= 0 && iup >= 0) {
	      blk = ion->iblock[r->f];
	      br = blk->r[ion->ilev[r->f]];
	      rt = blk->total_rate[ion->ilev[r->f]];
	      ra = wr[ir[RC_EA]+ig*nb+ilo];
	      if (ra == NULL) {
		ra = malloc(sizeof(double)*nt);
		for (s = 0; s < nt; s++) ra[s] = 0.0;
		wr[ir[RC_EA]+ig*nb+ilo] = ra;
		nea++;
	      }
	      if (br < 1) {
		ra[it] += r->dir*(1-br);
	      }
	      if (mea == 0) {
		ra0 = ra;
		for (n = 0; n < nbai[iup]; n++) {
		  ib = IdxGet(&iid, bai[iup][n].f);
		  if (ib >= 0 && ib != ig) {
		    ra = wr[ir[RC_EA]+ib*nb+ilo];
		    if (ra == NULL) {
		      ra = malloc(sizeof(double)*nt);
		      for (s = 0; s < nt; s++) ra[s] = 0.0;
		      wr[ir[RC_EA]+ib*nb+ilo] = ra;
		      nea++;
		    }
		    br = bai[iup][n].dir/rt;
		    x = r->dir*br;
		    ra[it] += x;
		    ra0[it] -= x;
		  }
		}
	      }
	    }
	  }
	}
      }
      te *= rdt;
    }

    printf("nrs: %d %d %d %d %d %d\n", nce, nci, nrr, ndr, nre, nea);
    rh.nele = ion->nele;
    if (nce) {
      rh.type = RC_CE;
      rh.nde = 1;
      InitFile(f, &fh, &rh);
      for (ilo = 0; ilo < nk; ilo++) {
	for (iup = 0; iup < nk; iup++) {
	  ra = wr[ir[RC_CE]+iup*nk+ilo];
	  if (ra) {
	    for (it = 0; it < nt; it++) {
	      rc.rc[it] = ra[it]/1e10;
	    }
	    rc.lower = iad.d[ilo];
	    rc.upper = iad.d[iup];
	    WriteRCRecord(f, &rc);
	    free(ra);
	  }
	}
      }
      DeinitFile(f, &fh);
    }
    if (nci) {
      rh.type = RC_CI;
      rh.nde = 1;
      InitFile(f, &fh, &rh);
      for (ilo = 0; ilo < nk; ilo++) {
	for (iup = 0; iup < ni; iup++) {
	  ra = wr[ir[RC_CI]+iup*nk+ilo];
	  if (ra) {
	    for (it = 0; it < nt; it++) {
	      rc.rc[it] = ra[it];
	    }
	    rc.lower = iad.d[ilo];
	    rc.upper = iid.d[iup];
	    WriteRCRecord(f, &rc);
	    free(ra);
	  }
	}
      }      
      DeinitFile(f, &fh);
    }
    if (nrr) {
      rh.type = RC_RR;
      rh.nde = 1;
      InitFile(f, &fh, &rh);
      for (ilo = 0; ilo < nk; ilo++) {
	for (iup = 0; iup < ni; iup++) {
	  ra = wr[ir[RC_RR]+iup*nk+ilo];
	  if (ra) {
	    for (it = 0; it < nt; it++) {
	      rc.rc[it] = ra[it]/1e10;
	    }
	    rc.lower = iad.d[ilo];
	    rc.upper = iid.d[iup];
	    WriteRCRecord(f, &rc);
	    free(ra);
	  }
	}
      } 
      DeinitFile(f, &fh);     
    }    
    if (ndr && !mdrea) {
      rh.type = RC_DR;
      rh.nde = nd;
      InitFile(f, &fh, &rh);
      for (ilo = 0; ilo < nb; ilo++) {
	for (iup = 0; iup < ni; iup++) {
	  ra = wr[ir[RC_DR]+iup*nb+ilo];
	  if (ra) {
	    for (id = 0; id < nd; id++) {
	      for (it = 0; it < nt; it++) {
		rc.rc[id*nt+it] = ra[id*nt+it]/1e10;
	      }
	    }
	    rc.lower = ibd.d[ilo];
	    rc.upper = iid.d[iup];
	    WriteRCRecord(f, &rc);
	    free(ra);
	  }
	}
      }
      DeinitFile(f, &fh);
    }
    if (nre) {
      rh.type = RC_RE;
      rh.nde = 1;
      InitFile(f, &fh, &rh);
      for (ilo = 0; ilo < ni; ilo++) {
	for (iup = 0; iup < ni; iup++) {
	  ra = wr[ir[RC_RE]+iup*ni+ilo];
	  if (ra) {
	    for (it = 0; it < nt; it++) {
	      rc.rc[it] = ra[it]/1e10;
	    }
	    rc.lower = iid.d[ilo];
	    rc.upper = iid.d[iup];
	    WriteRCRecord(f, &rc);
	    free(ra);
	  }
	}
      }  
      DeinitFile(f, &fh);    
    }
    if (nea && !mdrea) {
      rh.type = RC_EA;
      rh.nde = 1;
      InitFile(f, &fh, &rh);
      for (ilo = 0; ilo < nb; ilo++) {
	for (iup = 0; iup < ni; iup++) {
	  ra = wr[ir[RC_EA]+iup*nb+ilo];
	  if (ra) {
	    de = (ion->energy[iid.d[iup]]-ion->energy[ibd.d[ilo]])*HARTREE_EV;
	    te = t0;
	    for (it = 0; it < nt; it++) {
	      x = ra[it];
	      x *= ion->j[ibd.d[ilo]]+1.0;
	      x /= ion->j[iid.d[iup]]+1.0;
	      x *= exp(de/te);
	      x *= 1.64156e-12*pow(te, -1.5);
	      rc.rc[it] = x;
	      te *= rdt;
	    }
	    rc.lower = ibd.d[ilo];
	    rc.upper = iid.d[iup];
	    WriteRCRecord(f, &rc);
	    free(ra);
	  }
	}
      }  
      DeinitFile(f, &fh);       
    }

    if (mdrea && (ndr || nea)) {
      rh.type = RC_TT;
      rh.nde = nd+1;
      InitFile(f, &fh, &rh);
      for (ilo = 0; ilo < nb; ilo++) {
	for (iup = 0; iup < ni; iup++) {
	  ra = wr[ir[RC_EA]+iup*nb+ilo];
	  ra0 = wr[ir[RC_DR]+iup*nb+ilo];
	  if (!ra && !ra0) continue;
	  if (ra) {
	    de = (ion->energy[iid.d[iup]]-ion->energy[ibd.d[ilo]])*HARTREE_EV;
	    te = t0;
	    for (it = 0; it < nt; it++) {
	      x = ra[it];
	      x *= ion->j[ibd.d[ilo]]+1.0;
	      x /= ion->j[iid.d[iup]]+1.0;
	      x *= exp(de/te);
	      x *= 1.64156e-12*pow(te, -1.5);
	      rc.rc[ntd+it] = x;
	      te *= rdt;
	    }
	    free(ra);
	  } else {
	    for (it = 0; it < nt; it++) {
	      rc.rc[ntd+it] = 0.0;
	    }
	  }
	  if (ra0) {	    
	    for (id = 0; id < nd; id++) {
	      for (it = 0; it < nt; it++) {
		rc.rc[id*nt+it] = ra0[id*nt+it]/1e10;
	      }
	    }
	    free(ra0);
	  } else {
	    for (id = 0; id < nd; id++) {
	      for (it = 0; it < nt; it++) {
		rc.rc[id*nt+it] = 0.0;
	      }
	    }
	  }
	  rc.lower = ibd.d[ilo];
	  rc.upper = iid.d[iup];
	  WriteRCRecord(f, &rc);
	}
      }  
      DeinitFile(f, &fh);
    }
    
    if (nk > 0) {
      FreeIdxAry(&iad, 0);
    }
    if (nb > 0) {
      FreeIdxAry(&ibd, 0);
    }
    if (ni > 0) {
      FreeIdxAry(&iid, 0);
    }
    if (na > 0) {
      if (nbtr) {
	for (i = 0; i < na; i++) {
	  if (nbtr[i] > 0) free(btr[i]);
	}
	free(btr);
	free(nbtr);
      }
      if (nbai) {
	for (i = 0; i < na; i++) {
	  if (nbai[i] > 0) free(bai[i]);
	}
	free(bai);
	free(nbai);
      }
      FreeIdxAry(&idd, 0);
    }
    free(wr);
    free(drs);
    _krc = -1;
  }
  
  CloseFile(f, &fh);
  free(rc.rc);
}

