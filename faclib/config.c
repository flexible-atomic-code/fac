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

#include "config.h"
#include "angular.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

/*************************************************************
  Implementation of module "config".

  This module generates electron configuations and 
  carries out the angular momentum coupling. 

  Author: M. F. Gu, mfgu@space.mit.edu
**************************************************************/

/*
** VARIABLE:    cfg_groups
** TYPE:        static array
** PURPOSE:     a list of configuration groups.
** NOTE:        
*/
static CONFIG_GROUP *cfg_groups;

/*
** VARIABLE:    n_groups
** TYPE:        static int
** PURPOSE:     number of groups present.
** NOTE:        
*/
static int n_groups; 

/*
** VARIABLE:    symmetry_list
** TYPE:        static array
** PURPOSE:     a list of symmetries.
** NOTE:        the symmetry array is initialized in InitConfig.
**              the number of symmetries are fixed at 512. the i-th 
**              symmetry have the j = floor(i/2), 
**              and the parity = mod(i, 2). 
*/
static SYMMETRY *symmetry_list;

/*
** VARIABLE:    spec_symbols
** TYPE:        string.
** PURPOSE:     a string orbital angular momentum symbols.
** NOTE:        the last char "*" is not part of the symbol, 
**              rather, it represents any of the previous symbols.
*/
static char spec_symbols[MAX_SPEC_SYMBOLS+2] = "spdfghiklmnoqrtuvwxyz*"; 

#define NJQ 25
#define NJ2 10
static int _ja[25] = {1, 3, 3, 5, 5, 5, 7, 7, 7, 7, 9, 9, 9, 9, 9,
		      11, 11, 13, 13, 15, 15, 17, 17, 19, 19};
static int _qa[25] = {1, 1, 2, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 5,
		      1, 2, 1, 2, 1, 2, 1, 2};
static int _ji[NJ2] = {0, 1, 3, 6, 10, 15, 17, 19, 21, 23};
static SHELL_STATE *_csf1[NJQ];
static SHELL_STATE *_csf2[NJQ][NJQ];
static SHELL_STATE *_csf3[NJQ][NJQ][NJQ];
static SHELL_STATE *_csf4[NJQ][NJQ][NJQ][NJQ];
static SHELL_STATE *_csf5[NJQ][NJQ][NJQ][NJQ][NJQ];
static int _ncsf1[NJQ];
static int _ncsf2[NJQ][NJQ];
static int _ncsf3[NJQ][NJQ][NJQ];
static int _ncsf4[NJQ][NJQ][NJQ][NJQ];
static int _ncsf5[NJQ][NJQ][NJQ][NJQ][NJQ];
#define NCS 25
static int _isclosed[NCS][NCS];
			  
#define MAXN1 50
#define MAXN2 50
#define MAXNK1 2501
#define MAXNK2 2501
#define MAXNK12 (MAXNK1*MAXNK2)
static ARRAY _csary[MAXNK12];

int IsClosedComplex(int n, int nq) {
  if (n >= NCS) return 0;
  if (nq < 2*n*n) return 0;
  int i;
  for (i = 0; i < n; i++) {
    if (!_isclosed[n-1][i]) return 0;
  }
  return 1;
}

int IsClosedShellNR(int n, int k, int nq) {
  if (n >= NCS) return 0;
  if (nq < 2*(k*2+1)) return 0;
  return _isclosed[n-1][k];
}

int IsClosedShellFR(int n, int k, int j, int nq) {
  if (n >= NCS) return 0;
  if (nq < j+1) return 0;
  return _isclosed[n-1][k];
}

void SetClosedShellNR(int n, int k) {
  int i, j;
  if (n <= 0) {
    for (i = 0; i < NCS; i++) {
      for (j = 0; j < NCS; j++) {
	_isclosed[i][j] = 0;
      }
    }
    return;
  }
  if (n >= NCS) return;
  _isclosed[n-1][k] = 1;
}

int IndexJQ(int j, int q) {
  int i, jm, jp;  
  jp = j+1;
  if (q == 0 || q >= jp) return -1;
  if (j > 19) return -2;
  jm = (j-1)/2;
  if (q > jp/2) q = jp-q;
  if (j > 9 && q > 2) return -2;
  i = _ji[jm] + q-1;
  return i;
}
    
static void InitConfigData(void *p, int n) {
  CONFIG *d;
  int i;

  d = (CONFIG *) p;
  for (i = 0; i < n; i++) {
    d[i].n_shells = 0;
    d[i].n_csfs = 0;
    d[i].nnrs = 0;
    d[i].nrs = NULL;
    d[i].symstate = NULL;
    d[i].shells = NULL;
    d[i].csfs = NULL;
    d[i].igroup = -1;
    d[i].icfg = -1;    
  }
}
  
int ShellDegeneracy(int g, int nq) {
  if (nq == 1) {
    return g;
  } else if (nq == g) {
    return 1;
  } else {
    return (int) (exp(LnFactorial(g)-LnFactorial(nq)-LnFactorial(g-nq))+0.5);
  }
}

/* 
** FUNCTION:    DistributeElectronsShell
** PURPOSE:     distribute nq electrons among the specified shells
**              to construct all possible configurations.
** INPUT:       {CONFIG **cfg},
**              pointer to a pointer of CONFIG, which holds the
**              resulting configrations on output.
**              {int ns},
**              number of shells.
**              {SHELL *shells}
**              an array of shells.
**              {int nq},
**              number of electrons to be distributed.
**              {int *maxq},
**              maxq[i] is the maximum number of electrons allowed 
**              in all the shells except the the i-th shell. 
**              this is to determine the minimum number of electrons
**              allowed in the i-th shell, nq-maxq[i]. 
** RETURN:      {int},
**              number of configurations.
** SIDE EFFECT: 
** NOTE:        This is a static function only used in module "config".
*/
static int DistributeElectronsShell(CONFIG **cfg, int ns, SHELL *shell, 
				int nq, int *maxq) {
  CONFIG **cfg1, **cfg2;
  int *ncfg2;
  int qmin, qmax, j, q, t, k, ncfg;

  if (nq < 0) return 0;

  if (nq == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ns);
    (*cfg)->n_shells = 0;
    return 1;
  } 

  if (nq == 1){
    *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ns);
    for (t = 0; t < ns; t++) {
      (*cfg)[t].n_shells = 1;
      (*cfg)[t].shells = (SHELL *) malloc(sizeof(SHELL));
      (*cfg)[t].shells[0].n = shell[t].n;
      (*cfg)[t].shells[0].kappa = shell[t].kappa;
      (*cfg)[t].shells[0].nq = 1;
    }
    return ns;
  }

  if (ns == 1) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = ns;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    (*cfg)->shells[0].n = shell[0].n;
    (*cfg)->shells[0].kappa = shell[0].kappa;
    (*cfg)->shells[0].nq = nq;
    return 1;
  } 

  j = GetJFromKappa(shell[0].kappa);
  qmax = j+1;
  qmax = Min(qmax, nq);
  qmin = nq - maxq[0];
  qmin = Max(qmin, 0);
  ncfg = 0;
  t = qmax-qmin+1;
  cfg1 = (CONFIG **) malloc(sizeof(CONFIG *)*t);
  cfg2 = (CONFIG **) malloc(sizeof(CONFIG *)*t);
  ncfg2 = (int *) malloc(sizeof(int)*t);
  t = 0;
  for (q = qmin; q <= qmax; q++) {
    DistributeElectronsShell(cfg1+t, 1, shell, q, NULL);
    ncfg2[t] = DistributeElectronsShell(cfg2+t, ns-1, shell+1, nq-q, maxq+1);
    ncfg += ncfg2[t];
    t++;
  }
  
  *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ncfg);
  t = 0;
  k = 0;
  for (q = qmin; q <= qmax; q++) {
    for (j = 0; j < ncfg2[t]; j++) {
      (*cfg)[k].n_shells = cfg1[t]->n_shells + cfg2[t][j].n_shells;
      (*cfg)[k].shells = (SHELL *) malloc(sizeof(SHELL)*(*cfg)[k].n_shells);
      if (cfg1[t]->n_shells > 0) {
	memcpy((*cfg)[k].shells, cfg1[t]->shells, sizeof(SHELL));
      }
      if (cfg2[t][j].n_shells > 0) {
	memcpy((*cfg)[k].shells+cfg1[t]->n_shells, cfg2[t][j].shells, 
	       sizeof(SHELL)*(cfg2[t][j].n_shells));
      }
      if (cfg2[t][j].n_shells > 0) free(cfg2[t][j].shells);
      k++;
    }
    if (cfg1[t]->n_shells > 0) free(cfg1[t]->shells);
    free(cfg1[t]);
    free(cfg2[t]);
    t++;
  }
  free(cfg1);
  free(cfg2);
  free(ncfg2);
  
  return ncfg;
}

static int DistributeElectronsShellNR(CONFIG **cfg, int ns, SHELL *shell, 
				      int nq, int *maxq) {
  CONFIG **cfg1, **cfg2;
  int *ncfg2;
  int qmin, qmax, j, q, t, k, ncfg;
  
  if (nq < 0) return 0;

  if (nq == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ns);
    (*cfg)->n_shells = 0;
    return 1;
  } 

  if (nq == 1){
    *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ns);
    for (t = 0; t < ns; t++) {
      (*cfg)[t].n_shells = 1;
      (*cfg)[t].shells = (SHELL *) malloc(sizeof(SHELL));
      (*cfg)[t].shells[0].n = shell[t].n;
      (*cfg)[t].shells[0].kappa = shell[t].kappa;
      (*cfg)[t].shells[0].nq = 1;
    }
    return ns;
  }

  if (ns == 1) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = ns;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    (*cfg)->shells[0].n = shell[0].n;
    (*cfg)->shells[0].kappa = shell[0].kappa;
    (*cfg)->shells[0].nq = nq;
    return 1;
  } 

  j = shell[0].kappa;
  qmax = 2.0*(j+1);
  qmax = Min(qmax, nq);
  qmin = nq - maxq[0];
  qmin = Max(qmin, 0);
  ncfg = 0;
  t = qmax-qmin+1;
  cfg1 = (CONFIG **) malloc(sizeof(CONFIG *)*t);
  cfg2 = (CONFIG **) malloc(sizeof(CONFIG *)*t);
  ncfg2 = (int *) malloc(sizeof(int)*t);
  t = 0;
  for (q = qmin; q <= qmax; q++) {
    DistributeElectronsShellNR(cfg1+t, 1, shell, q, NULL);
    ncfg2[t] = DistributeElectronsShellNR(cfg2+t, ns-1, shell+1, nq-q, maxq+1);
    ncfg += ncfg2[t];
    t++;
  }
  
  *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ncfg);
  t = 0;
  k = 0;
  for (q = qmin; q <= qmax; q++) {
    for (j = 0; j < ncfg2[t]; j++) {
      (*cfg)[k].n_shells = cfg1[t]->n_shells + cfg2[t][j].n_shells;
      (*cfg)[k].shells = (SHELL *) malloc(sizeof(SHELL)*(*cfg)[k].n_shells);
      if (cfg1[t]->n_shells > 0) {
	memcpy((*cfg)[k].shells, cfg1[t]->shells, sizeof(SHELL));
      }
      if (cfg2[t][j].n_shells > 0) {
	memcpy((*cfg)[k].shells+cfg1[t]->n_shells, cfg2[t][j].shells, 
	       sizeof(SHELL)*(cfg2[t][j].n_shells));
      }
      if (cfg2[t][j].n_shells > 0) free(cfg2[t][j].shells);
      k++;
    }
    if (cfg1[t]->n_shells > 0) free(cfg1[t]->shells);
    free(cfg1[t]);
    free(cfg2[t]);
    t++;
  }
  free(cfg1);
  free(cfg2);
  free(ncfg2);
  
  return ncfg;
}

int ShellsFromString(char *scfg, double *dnq, SHELL **shell) {
  char token[128];
  int r, brkpos, quotepos, next;
  int nn, nkl, nkappa;
  int n[16];
  int kl[512];
  int kappa[1024];
  int i, j, t, k, kl2, ns;
  char *s;
    
  SetParserQuote("[", "]");
  SetParserBreak(spec_symbols);
  SetParserWhite("");
  SetParserEscape('\0');

  next = 0;
  r = 0;
  r = Parse(token, 512, scfg, &next, &brkpos, &quotepos);
  if (quotepos == 0) {
    nn = StrSplit(token, ',');
    if (nn > 16) {
      printf("number of n's in a single shell must be <= 16\n");
      exit(1);
    }
    s = token;
    for (i = 0; i < nn; i++) {
      while (*s == ' ' || *s == '\t') s++;
      n[i] = atoi(s);
      while (*s) s++;
      s++;
    }
  } else {
    nn = 1;
    n[0] = atoi(token);
  }
  if (brkpos < 0) {
    r = Parse(token, 512, scfg, &next, &brkpos, &quotepos);
  }
  if (brkpos >= 0) {
    kl[0] = brkpos;
    if (brkpos == MAX_SPEC_SYMBOLS) {
      if (n[nn-1] >= 512) {
	printf("not all L-terms are allowed for n >= %d\n", 512);
	exit(1);
      }
      nkl = n[nn-1];
      for (i = 0; i < nkl; i++) {
	kl[i] = i;
      }
    } else {
      nkl = 1;
      kl[0] = brkpos;
    }
    nkappa = 0;
    for (i = 0; i < nkl; i++) {
      kl2 = 2*kl[i];
      if (kl2 > 0) {
	kappa[nkappa++] = GetKappaFromJL(kl2-1, kl2);
      }
      kappa[nkappa++] = GetKappaFromJL(kl2+1, kl2);
    }
    *dnq = atof(&(scfg[next]));
    if (*dnq < 0) {
      printf("negative occupation number, use brackets to indicate j-1/2 shell\n");
      exit(1);
    }
    if (*dnq == 0 && !(isdigit(scfg[next]))) *dnq = 1;
  } else if (quotepos == 0) {
    nkl = StrSplit(token, ',');
    if (nkl > 512) {
      printf("number of L-terms must < 512\n");
      exit(1);
    }
    s = token;
    nkappa = 0;
    for (k = 0; k < nkl; k++) {
      while (*s == ' ' || *s == '\t') s++;
      GetJLFromSymbol(s, &j, &kl[k]);
      kl2 = 2*kl[k];
      if (j != 1 && kl2 > 0) {
	kappa[nkappa++] = GetKappaFromJL(kl2-1, kl2);
      }
      if (j != -1) {
	kappa[nkappa++] = GetKappaFromJL(kl2+1, kl2);
      }
      while (*s) s++;
      s++;
    }
    *dnq = atof(&(scfg[next]));
    if (*dnq == 0 && !(isdigit(scfg[next]))) *dnq = 1;
  } else {
    return -1;
  }

  *shell = (SHELL *) malloc(sizeof(SHELL)*nn*nkappa);
  t = 0;
  for (i = nn-1; i >= 0; i--) {
    for (k = nkappa-1; k >= 0; k--) {
      kl2 = GetLFromKappa(kappa[k]);
      if (kl2/2 >= n[i]) continue;
      (*shell)[t].n = n[i];
      (*shell)[t].kappa = kappa[k];
      t++;
    }
  }
  ns = t;

  if (ns == 0) {
    free(*shell);
    return -1;
  }
  
  return ns;
}

int ShellsFromStringNR(char *scfg, double *dnq, SHELL **shell) {
  char token[128];
  int r, brkpos, quotepos, next;
  int nn, nkl, nkappa;
  int n[16];
  int kl[512];
  int kappa[1024];
  int i, j, t, k, kl2, ns;
  char *s;
      
  SetParserQuote("[", "]");
  SetParserBreak(spec_symbols);
  SetParserWhite("");
  SetParserEscape('\0');

  next = 0;
  r = 0;
  r = Parse(token, 512, scfg, &next, &brkpos, &quotepos);
  if (quotepos == 0) {
    nn = StrSplit(token, ',');
    if (nn > 16) {
      printf("number of n's in a single shell must be <= 16\n");
      exit(1);
    }
    s = token;
    for (i = 0; i < nn; i++) {
      while (*s == ' ' || *s == '\t') s++;
      n[i] = atoi(s);
      while (*s) s++;
      s++;
    }
  } else {
    nn = 1;
    n[0] = atoi(token);
  }
  if (brkpos < 0) {
    r = Parse(token, 512, scfg, &next, &brkpos, &quotepos);
  }
  if (brkpos >= 0) {
    kl[0] = brkpos;
    if (brkpos == MAX_SPEC_SYMBOLS) {
      if (n[nn-1] >= 512) {
	printf("not all L-terms are allowed for n >= %d\n", 512);
	exit(1);
      }
      nkl = n[nn-1];
      for (i = 0; i < nkl; i++) {
	kl[i] = i;
      }
    } else {
      nkl = 1;
      kl[0] = brkpos;
    }
    nkappa = 0;
    for (i = 0; i < nkl; i++) {
      kl2 = 2*kl[i];
      kappa[nkappa++] = kl2;
    }
    if (scfg[next] == '+' || scfg[next] == '-') next++;
    *dnq = atof(&(scfg[next]));
    if (*dnq == 0 && !(isdigit(scfg[next]))) *dnq = 1;
  } else if (quotepos == 0) {
    nkl = StrSplit(token, ',');
    if (nkl > 512) {
      printf("number of L-terms must < 512\n");
      exit(1);
    }
    s = token;
    nkappa = 0;
    for (k = 0; k < nkl; k++) {
      while (*s == ' ' || *s == '\t') s++;
      GetJLFromSymbol(s, &j, &kl[k]);
      kl2 = 2*kl[k];
      kappa[nkappa++] = kl2;
      while (*s) s++;
      s++;
    }
    *dnq = atof(&(scfg[next]));
    if (*dnq == 0 && !(isdigit(scfg[next]))) *dnq = 1;
  } else {
    return -1;
  }

  *shell = (SHELL *) malloc(sizeof(SHELL)*nn*nkappa);
  t = 0;
  for (i = nn-1; i >= 0; i--) {
    for (k = nkappa-1; k >= 0; k--) {
      kl2 = kappa[k];
      if (kl2/2 >= n[i]) continue;
      (*shell)[t].n = n[i];
      (*shell)[t].kappa = kappa[k];
      t++;
    }
  }
  ns = t;
  
  if (ns == 0) {
    free(*shell);
    return -1;
  }

  return ns;
}

int GetRestriction(char *scfg, SHELL_RESTRICTION **sr, int m) {
  int nc, i;
  double dnq;
  char *s;

  nc = StrSplit(scfg, ';');
  nc--;

  *sr = (SHELL_RESTRICTION *) malloc(sizeof(SHELL_RESTRICTION)*nc);
 
  for (i = 0; i < nc; i++) {     
    while (*scfg) scfg++;
    scfg++;
    s = scfg;
    (*sr)[i].nq = -1;
    while (*s) {
      switch (*s) {
      case '=':
	(*sr)[i].op = 0;
	*s = '\0';
	(*sr)[i].nq = atoi(s+1);
	break;
      case '<':
	(*sr)[i].op = -1;
	*s = '\0';
	(*sr)[i].nq = atoi(s+1);
	break;
      case '>':
	(*sr)[i].op = 1;
	*s = '\0';
	(*sr)[i].nq = atoi(s+1);
	break;
      default:
	break;
      }
      s++;
    }
    if (m == 0) {
      (*sr)[i].ns = ShellsFromString(scfg, &dnq, &((*sr)[i].shells));
    } else {
      (*sr)[i].ns = ShellsFromStringNR(scfg, &dnq, &((*sr)[i].shells));
    }
    if ((*sr)[i].nq < 0) {
      (*sr)[i].op = 0;
      (*sr)[i].nq = dnq;
    }
    scfg = s;
  }
  
  return nc;
}

int ApplyRestriction(int ncfg, CONFIG *cfg, int nc, SHELL_RESTRICTION *sr) {
  int i, j, k, t, nq, c;
  
  for (k = 0; k < nc; k++) {
    for (i = 0; i < ncfg; i++) {
      nq = 0;
      for (j = 0; j < cfg[i].n_shells; j++) {
	for (t = 0; t < sr[k].ns; t++) {
	  if (cfg[i].shells[j].n == sr[k].shells[t].n &&
	      cfg[i].shells[j].kappa == sr[k].shells[t].kappa) {
	    nq += cfg[i].shells[j].nq;
	  }
	}
      }
      switch (sr[k].op) {
      case -1:
	c = (nq < sr[k].nq);
	break;
      case 0:
	c = (nq == sr[k].nq);
	break;
      case 1:
	c = (nq > sr[k].nq);
	break;
      default:
	c = 0;
	break;
      }
      if (c == 0) {
	if (cfg[i].n_shells >= 0) {
	  cfg[i].n_shells = -1;
	  free(cfg[i].shells);
	  cfg[i].shells = NULL;
	}
      }
    }
  }
  
  i = 0;
  j = 0;
  while (i < ncfg) {
    if (cfg[i].n_shells < 0) {
      i++;
    } else {
      cfg[j].n_shells = cfg[i].n_shells;
      cfg[j].shells = cfg[i].shells;
      j++;
      i++;
    }
  }

  return j;
}
  
/* 
** FUNCTION:    DistributeElectrons
** PURPOSE:     Decode a single string shell, distribute electrons
**              among all physical shells if the configurations 
**              are not the average configurations, otherwise, the
**              average configurations is constructed with all 
**              shells present, and the number of electrons returned.
** INPUT:       {CONFIG **cfg},
**              pointer to a pointer to CONFIG, which holds the 
**              resulting configurations or average configurations.
**              {double *nq},
**              pointer to double, which will hold the total number of 
**              electrons for the average configuration. it should be 
**              passed in as NULL if the actual configurations to be
**              constructed. 
**              {char *},
**              a single string shell. 
** RETURN:      {int},
**              if actual configurations to be constructed (nq == NULL),
**              return the number of configurations constructed.
**              if average configuration is to be determined,
**              return the number of shells in the average configuration.
** SIDE EFFECT: 
** NOTE:        
*/    
int DistributeElectrons(CONFIG **cfg, double *nq, char *scfg) {
  SHELL *shell;
  int ncfg, *maxq, ns, nc, i, j, inq;
  double dnq;
  SHELL_RESTRICTION *sr;

  nc = GetRestriction(scfg, &sr, 0);
  
  ns = ShellsFromString(scfg, &dnq, &shell);
  if (ns <= 0) {
    if (nc > 0) {
      for (i = 0; i < nc; i++) {
	free(sr[i].shells);
      }
      free(sr);
    }
    return ns;
  }

  if (nq) {
    *nq = dnq;
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = ns;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL)*ns);
    memcpy((*cfg)->shells, shell, sizeof(SHELL)*ns);
    free(shell);
    return ns;
  }

  maxq = (int *) malloc(sizeof(int)*ns);
  maxq[ns-1] = 0;
  for (i = ns-2; i >= 0; i--) {
    j = GetJFromKappa(shell[i+1].kappa);
    maxq[i] = maxq[i+1] + j+1;
  }

  inq = (int) dnq;
  ncfg = DistributeElectronsShell(cfg, ns, shell, inq, maxq);

  free(shell);
  free(maxq);
  if (nc > 0) {
    ncfg = ApplyRestriction(ncfg, *cfg, nc, sr);
    for (i = 0; i < nc; i++) {
      free(sr[i].shells);
    }
    free(sr);
  }
  return ncfg;
}

int DistributeElectronsNR(CONFIG **cfg, char *scfg) {
  SHELL *shell;
  int ncfg, *maxq, ns, nc, i, inq;
  double dnq;
  SHELL_RESTRICTION *sr;

  nc = GetRestriction(scfg, &sr, 1);
  
  ns = ShellsFromStringNR(scfg, &dnq, &shell);
  
  maxq = (int *) malloc(sizeof(int)*ns);
  maxq[ns-1] = 0;
  for (i = ns-2; i >= 0; i--) {
    maxq[i] = maxq[i+1] + 2*(shell[i+1].kappa + 1);
  }

  inq = (int) dnq;
  ncfg = DistributeElectronsShellNR(cfg, ns, shell, (int)dnq, maxq);

  free(shell);
  free(maxq);

  if (nc > 0) {
    ncfg = ApplyRestriction(ncfg, *cfg, nc, sr);
    for (i = 0; i < nc; i++) {
      free(sr[i].shells);
    }
    free(sr);
  }

  return ncfg;
}

/* 
** FUNCTION:    GetConfigOrAverageFromString
** PURPOSE:     decode the string representation of configurations,
**              construct all possible configurations or average 
**              configuration.
** INPUT:       {CONFIG **cfg}
**              holds the resuting configurations or average configuration.
**              {double **nq}, 
**              return fractional occupation numbers of each shell in the
**              average configuration, if it is not NULL on input. 
**              {char *scfg},
**              a string representation of configurations.
** RETURN:      {int},
**              if nq == NULL, return the number of configurations 
**              constructed.
**              if (nq != NULL, return the number of shells in the 
**              average configuration.
** SIDE EFFECT: 
** NOTE:        
*/       
int GetConfigOrAverageFromString(CONFIG **cfg, double **nq, char *scfg) {
  CONFIG **dcfg, **p1;
  double *dnq, *p2, a, b;
  char *s;
  int *isp, ncfg, *dnc;  
  int size, size_old, tmp;
  int i, t, j, k, ns;

  StrTrim(scfg, '\0');
  ns = QuotedStrSplit(scfg, ' ', '[', ']');
  if (ns == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = 1;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    PackShell((*cfg)->shells, 1, 0, 1, 0);
    return 1;
  }

  dcfg = (CONFIG **) malloc(sizeof(CONFIG *)*ns);
  dnc = (int *) malloc(sizeof(int)*ns);
  if (nq) {
    dnq = (double *) malloc(sizeof(double)*ns);
    p2 = dnq;
  } else {
    dnq = NULL;
    p2 = NULL;
  }

  s = scfg;
  isp = (int *) malloc(sizeof(int)*ns);  
  isp[0] = 0;
  for (i = 1; i < ns; i++) {
    isp[i] = isp[i-1];
    while (s[isp[i]]) isp[i]++;
    isp[i]++;
  }
  p1 = dcfg;
  t = 0;
  for (i = 0; i < ns; i++) {
    s = scfg + isp[i];
    while (*s == ' ' || *s == '\t') s++;
    dnc[t] = DistributeElectrons(p1, p2, s);
    if (dnc[t] <= 0) {
      free(dcfg);
      free(dnc);
      free(isp);
      return -1;
    }
    if (dnc[t] > 1 || (*p1)->n_shells > 0) {
      p1++;
      t++;
      if (p2) p2++;
    } else {
      free(*p1);      
    }
  }
  free(isp);

  ns = t;
  if (ns == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = 1;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    PackShell((*cfg)->shells, 1, 0, 1, 0);
    free(dcfg);
    free(dnc);
    return 1;
  }

  if (!nq) {
    ncfg = dnc[0];
    for (i = 1; i < ns; i++) ncfg *= dnc[i];
    *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ncfg);
    tmp = ncfg;
    p1 = dcfg + ns - 1;
    for (i = ns-1; i >= 0; i--) {      
      tmp /= dnc[i];
      t = 0;
      while (t < ncfg) {
	for (j = 0; j < dnc[i]; j++) {
	  for (k = 0; k < tmp; k++) {
	    if (i == ns-1) {
	      (*cfg)[t].n_shells = (*p1)[j].n_shells;
	      size = sizeof(SHELL)*(*p1)[j].n_shells;
	      (*cfg)[t].shells = (SHELL *) malloc(size);
	      memcpy((*cfg)[t].shells, (*p1)[j].shells, size);
	    } else {
	      size_old = sizeof(SHELL)*(*cfg)[t].n_shells;
	      size = sizeof(SHELL)*(*p1)[j].n_shells;
	      (*cfg)[t].shells = (SHELL *) realloc((*cfg)[t].shells, 
						   size_old+size);
	      memcpy((*cfg)[t].shells+(*cfg)[t].n_shells,
		     (*p1)[j].shells, size);
	      (*cfg)[t].n_shells += (*p1)[j].n_shells;
	    }
	    t++;
	  }
	}
      }
      p1--;
    }
  } else {
    ncfg = dnc[0];
    for (i = 1; i < ns; i++) ncfg += dnc[i];
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)[0].n_shells = ncfg;
    (*cfg)[0].shells = (SHELL *) malloc(sizeof(SHELL)*ncfg);
    *nq = (double *) malloc(sizeof(double)*ncfg);
    p1 = dcfg + ns-1;
    t = 0;
    for (i = ns-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < dnc[i]; j++) {
	a += GetJFromKappa((*p1)->shells[j].kappa) + 1.0;
      }
      for (j = 0; j < dnc[i]; j++) {
	b = GetJFromKappa((*p1)->shells[j].kappa) + 1.0;
	(*nq)[t] = dnq[i]*b/a;
	memcpy((*cfg)->shells+t, (*p1)->shells+j, sizeof(SHELL));
	t++;
      }
      p1--;
    }
  }    
  
  for (i = 0; i < ns; i++) {
    if (!nq) {
      for (j = 0; j < dnc[i]; j++) {
	if ((dcfg[i][j]).n_shells > 0) {
	  free((dcfg[i][j]).shells);
	}
      } 
    } else {
      if (dcfg[i][0].n_shells > 0) {
	free((dcfg[i][0]).shells);
      }
    }
    free(dcfg[i]);
  }
  free(dcfg);
  free(dnc);

  return ncfg;
}

int GetConfigFromStringNR(CONFIG **cfg, char *scfg) {
  CONFIG **dcfg, **p1;
  double *dnq;
  char *s;
  int *isp, ncfg, *dnc;  
  int size, size_old, tmp;
  int i, t, j, k, ns;

  StrTrim(scfg, '\0');
  ns = QuotedStrSplit(scfg, ' ', '[', ']');
  if (ns == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = 1;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    PackShell((*cfg)->shells, 1, 0, 1, 0);
    (*cfg)->shells->kappa = 0;
    return 1;
  }

  dcfg = (CONFIG **) malloc(sizeof(CONFIG *)*ns);
  dnc = (int *) malloc(sizeof(int)*ns);
  dnq = NULL;

  s = scfg;
  isp = (int *) malloc(sizeof(int)*ns);  
  isp[0] = 0;
  for (i = 1; i < ns; i++) {
    isp[i] = isp[i-1];
    while (s[isp[i]]) isp[i]++;
    isp[i]++;
  }
  p1 = dcfg;
  t = 0;
  for (i = 0; i < ns; i++) {
    s = scfg + isp[i];
    while (*s == ' ' || *s == '\t') s++;
    dnc[t] = DistributeElectronsNR(p1, s);
    if (dnc[t] <= 0) {
      free(dcfg);
      free(dnc);
      free(isp);
      return -1;
    }
    if (dnc[t] > 1 || (*p1)->n_shells > 0) {      
      p1++;
      t++;
    } else {
      free(*p1);
    }
  }
  free(isp);

  ns = t;
  if (ns == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = 1;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    PackShell((*cfg)->shells, 1, 0, 1, 0);
    (*cfg)->shells->kappa = 0;
    free(dcfg);
    free(dnc);
    return 1;
  }

  ncfg = dnc[0];
  for (i = 1; i < ns; i++) ncfg *= dnc[i];
  *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ncfg);
  tmp = ncfg;
  p1 = dcfg + ns - 1;
  for (i = ns-1; i >= 0; i--) {      
    tmp /= dnc[i];
    t = 0;
    while (t < ncfg) {
      for (j = 0; j < dnc[i]; j++) {
	for (k = 0; k < tmp; k++) {
	  if (i == ns-1) {
	    (*cfg)[t].n_shells = (*p1)[j].n_shells;
	    size = sizeof(SHELL)*(*p1)[j].n_shells;
	    (*cfg)[t].shells = (SHELL *) malloc(size);
	    memcpy((*cfg)[t].shells, (*p1)[j].shells, size);
	  } else {
	    size_old = sizeof(SHELL)*(*cfg)[t].n_shells;
	    size = sizeof(SHELL)*(*p1)[j].n_shells;
	    (*cfg)[t].shells = (SHELL *) realloc((*cfg)[t].shells, 
						 size_old+size);
	    memcpy((*cfg)[t].shells+(*cfg)[t].n_shells,
		   (*p1)[j].shells, size);
	    (*cfg)[t].n_shells += (*p1)[j].n_shells;
	  }
	  t++;
	}
      }
    }
    p1--;
  }

  for (i = 0; i < ns; i++) {
    for (j = 0; j < dnc[i]; j++) {
      if ((dcfg[i][j]).n_shells > 0) {
	free((dcfg[i][j]).shells);
      }
    }
    free(dcfg[i]);
  }
  free(dcfg);
  free(dnc);

  return ncfg;
}

/* 
** FUNCTION:    GetConfigFromString
** PURPOSE:     construct all possible cofigurations from string.
** INPUT:       {CONFIG **cfg},
**              holds the resulting configurations.
**              {char *scfg},
**              string representation of the configuraitons.
** RETURN:      {int},
**              number of the resulting configurations.
** SIDE EFFECT: 
** NOTE:        
*/
int GetConfigFromString(CONFIG **cfg, char *scfg) {
  return GetConfigOrAverageFromString(cfg, NULL, scfg);
}

/* 
** FUNCTION:    GetAverageConfigFromString
** PURPOSE:     construct the average configuration from a string.
** INPUT:       {int **n, **kappa, double **nq},
**              a list of principle quantum numbers, angular 
**              quantum numbers, and the fractional occupation
**              numbers of the resulting average configuration.
**              {char *scfg},
**              string representation of the average configuration.
** RETURN:      {int},
**              number shells in the average configuration.
** SIDE EFFECT: 
** NOTE:        
*/
int GetAverageConfigFromString(int **n, int **kappa, 
			       double **nq, char *scfg) {
  CONFIG *cfg;
  int i, ns;

  ns = GetConfigOrAverageFromString(&cfg, nq, scfg);
  if (ns <= 0) return ns;

  *n = (int *) malloc(sizeof(int)*ns);
  *kappa = (int *) malloc(sizeof(int)*ns);
  
  for (i = 0; i < ns; i++) {
    (*n)[i] = cfg->shells[i].n;
    (*kappa)[i] = cfg->shells[i].kappa;
  }

  free(cfg->shells);
  free(cfg);

  return ns;
}

/* 
** FUNCTION:    GetJLFromSymbol
** PURPOSE:     decode the spectroscopic symbol for a shell
**              and retrieve the j and L values.
** INPUT:       {char *s},
**              the spectroscopic symbol.
**              {int *j},
**              holds the total angular momentum of the shell,
**              either +1 or -1, indicates whether it's *kl+1 or
**              *kl-1.
**              {int *kl},
**              holds the orbital angular momentum of the shell.
** RETURN:      {int},
**               0: success.
**              -1: the symobel unrecoginized.
** SIDE EFFECT: 
** NOTE:        if the "+/-" sign is not present in the symbol, 
**              the j-value returned is 0, which indicates either 
**              value can be taken.
*/
int GetJLFromSymbol(char *s, int *j, int *kl) {
  int i;
  char s0[16], *p;

  strncpy(s0, s, 16);
  p = s0;
  while (*p) p++;
  p--;  
  if ((*p) == '+') {
    if (j) *j = 1;
    *p = '\0';
  } else if ((*p) == '-') {
    if (j) *j = -1;
    *p = '\0';
  } else {
    if (j) *j = 0;
  }

  if (kl) {
    if (isdigit(s0[0])) *kl = atoi(s0);
    else {
      for (i = 0; i < MAX_SPEC_SYMBOLS; i++) {
	if (spec_symbols[i] == s0[0]) {
	  *kl = i;
	  return 0;
	}
      }
      return -1;
    }
  }
  return 0;
}

/* 
** FUNCTION:    SpecSymbol
** PURPOSE:     construct the spectroscopic symbol for the 
**              non-relativistic shell.
** INPUT:       {char *s},
**              string holding the result.
**              {int kl},
**              orbital angular momentum of the shell.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        if kl >= MAX_SPEC_SYMBOLS, then the symbol is 
**              returned as "[kl]".
*/
int SpecSymbol(char *s, int kl) {
  if (kl < MAX_SPEC_SYMBOLS) {
    s[0] = spec_symbols[kl];
    s[1] = '\0';
  } else {
    sprintf(s, "[%d]", kl);
  }
  return 0;
}

/* 
** FUNCTION:    Couple
** PURPOSE:     recursively construct all possible states for a Config.
** INPUT:       {CONFIG *cfg},
**              pointer to the config. to be coupled.
** RETURN:      {int},
**               0: success.
**              <0: error.
** SIDE EFFECT: 
** NOTE:        
*/
int Couple(CONFIG *cfg) {
  CONFIG outmost, inner;
  int errcode, i;
  int *idx = NULL;

  if (cfg->n_shells == 0) {
    errcode = -1;
    goto ERROR;
  }

  if (cfg == NULL) {
    errcode = -2;
    goto ERROR; 
  }

  /* make sure that the shells are sorted in inverse order */
  for (i = 1; i < cfg->n_shells; i++) {
    if (CompareShell(&cfg->shells[i-1], &cfg->shells[i]) < 0) {
      qsort(cfg->shells, cfg->n_shells, sizeof(SHELL), CompareShellInvert);
      break;
    }
  }
  if (IsUTA()) {
    cfg->csfs = NULL;
    cfg->n_csfs = 0;
    cfg->n_electrons = 0;
    for (i = 0; i < cfg->n_shells; i++) {
      cfg->n_electrons += cfg->shells[i].nq;
    }
    return 0;
  }

  int ns, j2, nq, id[5], s;
  idx = malloc(sizeof(int)*cfg->n_shells);
  ns = 0;
  for (i = 0; i < cfg->n_shells; i++) {
    j2 = GetJ(&cfg->shells[i]);
    nq = GetNq(&cfg->shells[i]);
    idx[i] = IndexJQ(j2, nq);
    if (idx[i] == -2) {
      ns = -1;
      break;
    } else if (idx[i] >= 0) {
      if (ns == 5) {
	ns = -1;
	break;
      }
      id[ns] = i;
      ns++;
    }
  }
  SHELL_STATE *csf;
  int ncsf;
  switch (ns) {
  case 1:
    ncsf = _ncsf1[idx[id[0]]];
    csf = _csf1[idx[id[0]]];
    break;
  case 2:
    ncsf = _ncsf2[idx[id[0]]][idx[id[1]]];
    csf = _csf2[idx[id[0]]][idx[id[1]]];
    break;
  case 3:
    ncsf = _ncsf3[idx[id[0]]][idx[id[1]]][idx[id[2]]];
    csf = _csf3[idx[id[0]]][idx[id[1]]][idx[id[2]]];
    break;
  case 4:
    ncsf = _ncsf4[idx[id[0]]][idx[id[1]]][idx[id[2]]][idx[id[3]]];
    csf = _csf4[idx[id[0]]][idx[id[1]]][idx[id[2]]][idx[id[3]]];
    break;
  case 5:
    ncsf = _ncsf5[idx[id[0]]][idx[id[1]]][idx[id[2]]][idx[id[3]]][idx[id[4]]];
    csf = _csf5[idx[id[0]]][idx[id[1]]][idx[id[2]]][idx[id[3]]][idx[id[4]]];
    break;
  default:
    ncsf = 0;
    csf = NULL;
  }
  if (csf != NULL) {
    cfg->n_csfs = ncsf;
    cfg->csfs = malloc(sizeof(SHELL_STATE)*cfg->n_shells*cfg->n_csfs);
    SHELL_STATE *p1 = cfg->csfs + cfg->n_shells*cfg->n_csfs - 1;
    SHELL_STATE *p0 = csf + ns*ncsf - 1;
    for (i = 0; i < ncsf; i++) {
      for (s = cfg->n_shells-1; s >= 0; s--) {
	if (idx[s] >= 0) {
	  PackShellState(p1, p0->totalJ, p0->shellJ, p0->nu, p0->Nr);
	  p0--;
	} else {
	  if (s == cfg->n_shells-1) {
	    PackShellState(p1, 0, 0, 0, 0);
	  } else {
	    PackShellState(p1, p1[1].totalJ, 0, 0, 0);
	  }
	}
	p1--;
      }
    }
    cfg->n_electrons = 0;
    for (i = 0; i < cfg->n_shells; i++) {
      cfg->n_electrons += cfg->shells[i].nq;
    }
    return 0;
  }
  if (cfg->n_shells == 1) {
    if (GetSingleShell(cfg) < 0) {
      errcode = -3;
      goto ERROR;
    }
  } else {
    outmost.n_shells = 1;
    outmost.shells = cfg->shells;
    inner.n_shells = cfg->n_shells - 1;
    inner.shells = cfg->shells + 1;
    if (Couple(&outmost) < 0) {
      errcode = -4;
      goto ERROR;
    }
    if (Couple(&inner) < 0) {
      errcode = -5;
      free(outmost.csfs);
      goto ERROR;
    }  
    
    if (CoupleOutmost(cfg, &outmost, &inner) < 0) {
      errcode = -5;
      free(outmost.csfs);
      free(inner.csfs);
      goto ERROR;
    }

    free(outmost.csfs);
    free(inner.csfs);
  }

  if (ns > 0) {
    ncsf = cfg->n_csfs;
    csf = malloc(sizeof(SHELL_STATE)*ns*ncsf);
    SHELL_STATE *p1 = cfg->csfs;
    SHELL_STATE *p0 = csf;
    for (i = 0; i < ncsf; i++) {
      for (s = 0; s < cfg->n_shells; s++) {
	if (idx[s] >= 0) {
	  PackShellState(p0, p1->totalJ, p1->shellJ, p1->nu, p1->Nr);
	  p0++;
	}
	p1++;
      }
    }
    switch (ns) {
    case 1:
      _ncsf1[idx[id[0]]] = ncsf;
      _csf1[idx[id[0]]] = csf;
      break;
    case 2:
      _ncsf2[idx[id[0]]][idx[id[1]]] = ncsf;
      _csf2[idx[id[0]]][idx[id[1]]] = csf;
      break;
    case 3:
      _ncsf3[idx[id[0]]][idx[id[1]]][idx[id[2]]] = ncsf;
      _csf3[idx[id[0]]][idx[id[1]]][idx[id[2]]] = csf;
      break;
    case 4:
      _ncsf4[idx[id[0]]][idx[id[1]]][idx[id[2]]][idx[id[3]]] = ncsf;
      _csf4[idx[id[0]]][idx[id[1]]][idx[id[2]]][idx[id[3]]] = csf;
      break;
    case 5:
      _ncsf5[idx[id[0]]][idx[id[1]]][idx[id[2]]][idx[id[3]]][idx[id[4]]] = ncsf;
      _csf5[idx[id[0]]][idx[id[1]]][idx[id[2]]][idx[id[3]]][idx[id[4]]] = csf;
      break;
    default:
      break;
    }
  }
  if (idx) free(idx);
  return 0;

 ERROR:
  printf("****Error in Couple****\n");
  if (idx) free(idx);
  return errcode;
}

/* 
** FUNCTION:    CoupleOutmost
** PURPOSE:     constructs all possible states by coupling 
**              the outmost shell to the inner shells.
** INPUT:       {CONFIG *cfg},
**              pointer to the resulting configuaration.
**              {CONFIG *outmost, *inner},
**              pointer to the configurations of the outmost shell
**              and the inner shells.
** RETURN:      {int},
**               0: success.
**              -1: error.
** SIDE EFFECT: 
** NOTE:        both outmost shell and inner shells must 
**              have been coupled.
*/
int CoupleOutmost(CONFIG *cfg, CONFIG *outmost, CONFIG *inner) {
  int i, j, k;
  int bytes_csf, bytes_csf_inner, bytes_csf_outmost;
  int j2_min, j2_max;
  int j2_inner, j2_outmost;
  SHELL_STATE *csf, *csf_outmost, *csf_inner;

  if (outmost->n_shells != 1) goto ERROR;
  if (cfg->n_shells != 1 + inner->n_shells) goto ERROR;

  if (inner->n_shells == 0) {
    cfg->n_csfs = outmost->n_csfs;
    cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
    if (cfg->csfs == NULL) goto ERROR;
    memcpy(cfg->csfs, outmost->csfs, cfg->n_csfs * sizeof(SHELL_STATE));
    return 0;
  }

  bytes_csf_outmost = sizeof(SHELL_STATE);
  bytes_csf_inner = inner->n_shells * sizeof(SHELL_STATE);
  bytes_csf = bytes_csf_inner + bytes_csf_outmost;

  /*************************************************************************
  First calculte the total # of possible states, and allocate memory for 
  cfg->csfs.
  *************************************************************************/
  csf_outmost = outmost->csfs;
  cfg->n_csfs = 0;
  for (i = 0; i < outmost->n_csfs; i++) {
    j2_outmost = csf_outmost->totalJ;
    csf_inner = inner->csfs;
    for (j = 0; j < inner->n_csfs; j++) {
      j2_inner = csf_inner->totalJ;
      j2_min = abs(j2_outmost - j2_inner);
      j2_max = j2_outmost + j2_inner;
      cfg->n_csfs += (j2_max - j2_min) / 2 + 1;
      csf_inner += inner->n_shells;
    }
    csf_outmost ++;
  }

  cfg->csfs = malloc(cfg->n_csfs * bytes_csf);
  if (cfg->csfs == NULL) goto ERROR;
  csf = cfg->csfs;

  /** Fill in the cfg->csfs array. **/
  csf_outmost = outmost->csfs;
  for (i = 0; i < outmost->n_csfs; i++) { 
    j2_outmost = csf_outmost->totalJ;
    csf_inner = inner->csfs;
    for (j = 0; j < inner->n_csfs; j++) {
      j2_inner = csf_inner->totalJ;
      j2_min = abs(j2_outmost - j2_inner);
      j2_max = j2_outmost + j2_inner;
      for (k = j2_min; k <= j2_max; k += 2) {
	memcpy(csf, csf_outmost, bytes_csf_outmost);
	csf->totalJ = k; 
	memcpy(csf + 1, csf_inner, bytes_csf_inner);
	csf += cfg->n_shells;
      }
      csf_inner += inner->n_shells;
    }
    csf_outmost++;
  }
  
  cfg->n_electrons = outmost->n_electrons + inner->n_electrons;

  return 0;
  
 ERROR:
  printf("****Error in CoupleOutmost****\n");
  return -1;
}

/* 
** FUNCTION:    GetSingleShell
** PURPOSE:     construct all possible states for a single shell.
** INPUT:       {CONFIG *cfg},
**              pointer to the resulting configuration for the 
**              single shell.
** RETURN:      {int},
**               0: success.
**              -1: error.
** SIDE EFFECT: 
** NOTE:        for j > 9/2, no more than 2 electrons are allowed. 
**              The data were taken from "Nuclear Shell Theory" by 
**              AMOS de-SHALIT and IGAL TALMI.
*/
int GetSingleShell(CONFIG *cfg) {
  int j2, max_q;
  int occupation;
  SHELL_STATE *csf;
  int i;

  if (cfg->n_shells != 1) goto ERROR;

  j2 = GetJ(cfg->shells);
  if (!(IsOdd(j2)) || j2 < 0) goto ERROR;

  max_q = j2 + 1;
  occupation = GetNq(cfg->shells);
  cfg->n_electrons = occupation;
  if ((2 * occupation) > max_q) {
    occupation = max_q - occupation;
  }
  if (occupation < 0) goto ERROR;
  switch(occupation) {
  case 0: /** 0 occupation or closed shell **/
    cfg->n_csfs = 1;
    cfg->csfs = malloc(sizeof(SHELL_STATE));
    if (cfg->csfs == NULL) goto ERROR;
    PackShellState(cfg->csfs, 0, 0, 0, 0);
    break;

  case 1: /** 1 electron **/
    cfg->n_csfs = 1;
    cfg->csfs = malloc(sizeof(SHELL_STATE));
    if (cfg->csfs == NULL) goto ERROR;   
    PackShellState(cfg->csfs, j2, j2, 1, 0);
    break;
    
  case 2: /** 2 equivelent electrons **/
    cfg->n_csfs = (j2 + 1) / 2;
    cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
    if (cfg->csfs == NULL) goto ERROR; 
    csf = cfg->csfs;
    PackShellState(csf++, 0, 0, 0, 0);
    for (i = 2; i < j2; i += 2) {
      PackShellState(csf++, i*2, i*2, 2, 0);
    }
    break;

  case 3: /** 3 equivelant electrons **/
    switch(j2) {
    case 5: /** j = 5/2 **/
      cfg->n_csfs = 3;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR; 
      csf = cfg->csfs;
      PackShellState(csf++, 5, 5, 1, 0);
      PackShellState(csf++, 3, 3, 3, 0);
      PackShellState(csf++, 9, 9, 3, 0);
      break;
    case 7: /** j = 7/2 **/
      cfg->n_csfs = 6;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR; 
      csf = cfg->csfs;
      PackShellState(csf++, 7, 7, 1, 0);
      PackShellState(csf++, 3, 3, 3, 0);
      PackShellState(csf++, 5, 5, 3, 0);
      PackShellState(csf++, 9, 9, 3, 0);
      PackShellState(csf++, 11, 11, 3, 0);
      PackShellState(csf++, 15, 15, 3, 0);
      break;
    case 9: /** j = 9/2 **/
      cfg->n_csfs = 10;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR; 
      csf = cfg->csfs;
      PackShellState(csf++, 9, 9, 1, 0);
      PackShellState(csf++, 3, 3, 3, 0);
      PackShellState(csf++, 5, 5, 3, 0);
      PackShellState(csf++, 7, 7, 3, 0);
      PackShellState(csf++, 9, 9, 3, 0);
      PackShellState(csf++, 11, 11, 3, 0);
      PackShellState(csf++, 13, 13, 3, 0);
      PackShellState(csf++, 15, 15, 3, 0);
      PackShellState(csf++, 17, 17, 3, 0);
      PackShellState(csf++, 21, 21, 3, 0);
      break;
    default:
      goto ERROR;
    }
    break;

  case 4:
    switch(j2) {
    case 7: /** j = 7/2 **/
      cfg->n_csfs = 8;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR; 
      csf = cfg->csfs;
      PackShellState(csf++, 0, 0, 0, 0);
      PackShellState(csf++, 4, 4, 2, 0);
      PackShellState(csf++, 8, 8, 2, 0);
      PackShellState(csf++, 12, 12, 2, 0);
      PackShellState(csf++, 4, 4, 4, 0);
      PackShellState(csf++, 8, 8, 4, 0);
      PackShellState(csf++, 10, 10, 4, 0);
      PackShellState(csf++, 16, 16, 4, 0);
      break;
    case 9: /** j = 9/2 **/
      cfg->n_csfs = 18;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR; 
      csf = cfg->csfs;
      PackShellState(csf++, 0, 0, 0, 0);
      PackShellState(csf++, 4, 4, 2, 0);
      PackShellState(csf++, 8, 8, 2, 0);
      PackShellState(csf++, 12, 12, 2, 0);
      PackShellState(csf++, 16, 16, 2, 0);
      PackShellState(csf++, 0, 0, 4, 0);
      PackShellState(csf++, 4, 4, 4, 0);
      PackShellState(csf++, 6, 6, 4, 0);
      PackShellState(csf++, 8, 8, 4, 0);
      PackShellState(csf++, 8, 8, 4, 1);
      PackShellState(csf++, 10, 10, 4, 0);
      PackShellState(csf++, 12, 12, 4, 0);
      PackShellState(csf++, 12, 12, 4, 1);
      PackShellState(csf++, 14, 14, 4, 0);
      PackShellState(csf++, 16, 16, 4, 0);
      PackShellState(csf++, 18, 18, 4, 0);
      PackShellState(csf++, 20, 20, 4, 0);
      PackShellState(csf++, 24, 24, 4, 0);
      break;
    default:
      goto ERROR;
    }
    break;

  case 5:  
    switch(j2) {
    case 9: /** only j = 9/2 is allowed **/
      cfg->n_csfs = 20; 
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR; 
      csf = cfg->csfs;
      PackShellState(csf++, 9, 9, 1, 0);
      PackShellState(csf++, 3, 3, 3, 0);
      PackShellState(csf++, 5, 5, 3, 0);
      PackShellState(csf++, 7, 7, 3, 0);
      PackShellState(csf++, 9, 9, 3, 0);
      PackShellState(csf++, 11, 11, 3, 0);
      PackShellState(csf++, 13, 13, 3, 0);
      PackShellState(csf++, 15, 15, 3, 0);
      PackShellState(csf++, 17, 17, 3, 0);
      PackShellState(csf++, 21, 21, 3, 0);
      PackShellState(csf++, 1, 1, 5, 0);
      PackShellState(csf++, 5, 5, 5, 0);
      PackShellState(csf++, 7, 7, 5, 0);
      PackShellState(csf++, 9, 9, 5, 0);
      PackShellState(csf++, 11, 11, 5, 0);
      PackShellState(csf++, 13, 13, 5, 0);
      PackShellState(csf++, 15, 15, 5, 0);
      PackShellState(csf++, 17, 17, 5, 0);
      PackShellState(csf++, 19, 19, 5, 0);
      PackShellState(csf++, 25, 25, 5, 0);
      break;
    default:
      goto ERROR;
    }
    break;

  default:
    goto ERROR;

  }
  
  return 0;
  
 ERROR:
  printf("****Error in GetSingleShell****\n");
  return -1;
}

/* 
** FUNCTION:    PackShell, UnpackShell
** PURPOSE:     pack and unpack the fields of SHELL.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
void UnpackShell(SHELL *s, int *n, int *kl, int *j, int *nq) {
  *n = s->n;
  *nq = s->nq;
  *j = 2*abs(s->kappa) - 1;
  *kl = (s->kappa < 0)? (*j - 1):(*j + 1);
}

void PackShell(SHELL *s, int n, int kl, int j, int nq){
  s->n = n;
  s->nq = nq;
  s->kappa = (kl - j)*(j + 1)/2;
}

void UnpackNRShell(int *s, int *n, int *kl, int *nq) {
  *nq = (*s)&0xFF;
  *kl = 2*(((*s)>>8)&0xFF);
  *n = ((*s)>>16)&0xFF;
}

void PackNRShell(int *s, int n, int kl, int nq) {
  *s = (n<<16) | ((kl/2)<<8) | nq;
}

/* 
** FUNCTION:    GetNq, GetJ, GetL
** PURPOSE:     retrieve nq, j, and L of a shell.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int GetNq(SHELL *s){
  return s->nq;
}

int GetJ(SHELL *s){
  return 2*abs(s->kappa) - 1;
}

int GetL(SHELL *s){
  int j;
  j = 2*abs(s->kappa) - 1;
  return (s->kappa < 0)? (j - 1):(j + 1);
} 

/* 
** FUNCTION:    ShellClosed
** PURPOSE:     determine if the shell is a closed one.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int ShellClosed(SHELL *s) {
  int j;
  j = GetJ(s);
  if (s->nq < j+1) return 0;
  return 1;
}

/* 
** FUNCTION:    GetLFromKappa, GetJFromKappa, 
**              GetKappaFromJL, GetJLFromKappa
** PURPOSE:     convert between kappa and JL values.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int GetLFromKappa(int kappa) {
  int j;
  j = 2*abs(kappa) - 1;
  return (kappa < 0)? (j - 1):(j + 1);
}

int GetJFromKappa(int kappa) {
  return 2*abs(kappa) - 1;
}

int GetKappaFromJL(int j, int kl) {
  if (j <= 0 || kl < 0) return 0;
  return (kl-j)*(j+1)/2;
}

void GetJLFromKappa(int kappa, int *j, int *kl) {
  *j = 2*abs(kappa) - 1;
  *kl = (kappa < 0)? (*j - 1):(*j + 1);
}

/* 
** FUNCTION:    PackShellState
** PURPOSE:     pack fields of SHELL_STATE to the structure.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
void PackShellState(SHELL_STATE *s, int J, int j, int nu, int Nr){
  s->totalJ = J;
  s->shellJ = j;
  s->nu = nu;
  s->Nr = Nr;
}

int ShellNeedNuNr(SHELL *s, SHELL_STATE *st) {
  if (s->nq < 4) return 0;
  int j = GetJFromKappa(s->kappa);
  if (j < 7) return 0;
  if (j == 7) {
    if (st->shellJ == 4 || st->shellJ == 8) return 1;
    return 0;
  } else if (j == 9) {    
    if (s->nq == 4) {
      if (st->shellJ == 6 || st->shellJ == 10 ||
	  st->shellJ == 14 || st->shellJ == 16 || st->shellJ == 18 ||
	  st->shellJ == 20 || st->shellJ == 24) {
	return 0;
      }
      if (st->shellJ == 8 || st->shellJ == 12) return 3;
      return 1;
    }
    if (s->nq == 5) {
      if (st->shellJ == 1 ||
	  st->shellJ == 3 ||
	  st->shellJ == 19 ||
	  st->shellJ == 21 ||
	  st->shellJ == 25) {
	return 0;
      }
      return 1;
    }
  }
  return 0;
}
  
/* 
** FUNCTION:    GroupIndex
** PURPOSE:     find the index of the group by its name.
** INPUT:       {char *name},
**              the group name.
** RETURN:      {int},
**              the group index.
** SIDE EFFECT: if the group does not exist, a new one is created.
** NOTE:        
*/
int GroupIndex(char *name) {
  int i;

  for (i = n_groups - 1; i >= 0; i--) {
    if (strncmp(name, cfg_groups[i].name, GROUP_NAME_LEN) == 0) 
      break;
  }
  if (i < 0) i = AddGroup(name);
  return i;
}

/* 
** FUNCTION:    GroupExists
** PURPOSE:     determine if a group exists.
** INPUT:       {char *name},
**              the group name.
** RETURN:      {int},
**              >=0: the group index, if it exists.
**               <0: the group does not exist.
** SIDE EFFECT: 
** NOTE:        
*/
int GroupExists(char *name) {
  int i;

  for (i = n_groups - 1; i >= 0; i--) {
    if (strncmp(name, cfg_groups[i].name, GROUP_NAME_LEN) == 0) 
      break;
  }
  return i;
}

/* 
** FUNCTION:    AddGroup
** PURPOSE:     add a group to the group array.
** INPUT:       {char *name},
**              the group name.
** RETURN:      {int},
**              the index of the added group
** SIDE EFFECT: 
** NOTE:        
*/
int AddGroup(char *name) {
  if (name == NULL) return -1;
  if (n_groups == MAX_GROUPS) {
    printf("Max # groups reached\n");
    exit(1);
  }
  strncpy(cfg_groups[n_groups].name, name, GROUP_NAME_LEN);
  cfg_groups[n_groups].nmax = 0;
  n_groups++;
  return n_groups-1;
}

/* 
** FUNCTION:    GetGroup
** PURPOSE:     retrieve the pointer to the group by its index.
** INPUT:       {int k},
**              the index of the group.
** RETURN:      {CONFIG_GROUP *},
**              the pointer to the group.
** SIDE EFFECT: 
** NOTE:        
*/
CONFIG_GROUP *GetGroup(int k) {
  if (k < 0 || k >= n_groups) return NULL;
  return cfg_groups+k;
}

/* 
** FUNCTION:    GetNewGroup
** PURPOSE:     add a new group and return its pointer.
** INPUT:       
** RETURN:      {CONFIG_GROUP *},
**              pointer to the added group.
** SIDE EFFECT: 
** NOTE:        the name of the group is initialized as '_all_'.
*/
CONFIG_GROUP *GetNewGroup(void) {
  if (n_groups == MAX_GROUPS) {
    printf("Max # groups reached\n");
    exit(1);
  }
  n_groups++;
  return cfg_groups+n_groups-1;
}

/* 
** FUNCTION:    GetNumGroups
** PURPOSE:     retrieve the number of groups in the array.
** INPUT:       
** RETURN:      {int},
**              the number of groups.
** SIDE EFFECT: 
** NOTE:        
*/
int GetNumGroups(void) {
  return n_groups;
}

/* 
** FUNCTION:    GetConfig
** PURPOSE:     return a pointer to CONFIG, which a state belongs to.
** INPUT:       {STATE *s},
**              pointer to a state.
** RETURN:      {CONFIG *},
**              the configuration the state belongs to.
** SIDE EFFECT: 
** NOTE:        
*/
CONFIG *GetConfig(STATE *s) {
  CONFIG *c;
  int i, j;

  i = s->kgroup;
  j = s->kcfg;
  c = (CONFIG *) ArrayGet(&(cfg_groups[i].cfg_list), j);
  return c;
}

CONFIG *GetConfigFromGroup(int kg, int kc) {
  return (CONFIG *) ArrayGet(&(cfg_groups[kg].cfg_list), kc);
}

int ConfigToIList(CONFIG *c, int n, int *s) {
  int i, j;
  for (i = 0; i < n; i++) s[i] = 0;
  for (i = 0; i < c->n_shells; i++) {
    j = ShellToInt(c->shells[i].n, c->shells[i].kappa);
    if (j >= n) {
      printf("ConfigToIList error: %d %d\n", j, n);
      return -1;
    }
    s[j] = c->shells[i].nq;
  }
  return 0;
}

CONFIG *ConfigFromIList(int n, int *s) {
  CONFIG *c;
  int i, j;
  c = malloc(sizeof(CONFIG));
  InitConfigData(c, 1);
  for (i = 0; i < n; i++) {
    if (s[i]) c->n_shells++;
  }
  c->shells = malloc(sizeof(SHELL)*c->n_shells);
  j = 0;
  for (i = n-1; i >= 0; i--) {
    if (s[i]) {
      IntToShell(i, &c->shells[j].n, &c->shells[j].kappa);
      c->shells[j].nq = s[i];
      j++;      
    }
  }
  return c;
}

int ConfigExists(CONFIG *cfg) {
  int i, j, t, n, nele;
  CONFIG *c;
  int i0, i1, n1;

  n = cfg->shells[0].n;
  if (n > MAXN1) {
    i0 = MAXNK1-1;
  } else {
    i0 = ShellToInt(cfg->shells[0].n, cfg->shells[0].kappa);
  }
  if (cfg->shells[0].nq > 1) {
    i1 = i0;
  } else if (cfg->n_shells > 1 && cfg->shells[1].n <= MAXN2) {
    i1 = ShellToInt(cfg->shells[1].n, cfg->shells[1].kappa);
  } else {
    i1 = MAXNK2-1;
  }
  n1 = i0*MAXNK2 + i1;
  ARRAY *csa = &_csary[n1];
  if (csa->dim == 0) return 0;
  nele = 0;
  for (i = 0; i < cfg->n_shells; i++) {
    nele += cfg->shells[i].nq;
  }
  for (i = 0; i < csa->dim; i++) {
    CONFIG **p = ArrayGet(csa, i);
    c = *p;
    if (nele != c->n_electrons) continue;
    if (cfg->n_shells != c->n_shells) continue;
    for (t = 0; t < c->n_shells; t++) {
      if (cfg->shells[t].n != c->shells[t].n) break;
      if (cfg->shells[t].kappa != c->shells[t].kappa) break;
      if (cfg->shells[t].nq != c->shells[t].nq) break;
    }
    if (t == c->n_shells) return 1;
  }
  return 0;
}

/* 
** FUNCTION:    AddConfigToList
** PURPOSE:     add a configuration to the specified group,
**              and add all states to the symmetry list.
** INPUT:       {int k},
**              the group index where the config is add to.
**              {CONFIG *cfg},
**              pointer to CONFIG to be added.
** RETURN:      {int},
**               0: success.
**              -1: error.
** SIDE EFFECT: 
** NOTE:        
*/
int AddConfigToList(int k, CONFIG *cfg) {
  ARRAY *clist;  
  int n0, kl0, nq0, m, i, n, kl, j, nq;
  if (k < 0 || k >= n_groups) return -1;
  if (cfg_groups[k].n_cfgs == 0) {
    cfg_groups[k].n_electrons = cfg->n_electrons;
  } else if (cfg_groups[k].n_electrons != cfg->n_electrons) {
    printf("Error: AddConfigToList, Configurations in a group ");
    printf("must have the same number of electrons\n");
    return -1;
  }
  clist = &(cfg_groups[k].cfg_list);

  cfg->energy = 0.0;
  cfg->delta = 0.0;

  n0 = 0;
  kl0 = -1;
  nq0 = 0;
  m = 0;
  cfg->nrs = malloc(sizeof(int)*cfg->n_shells);
  cfg->sweight = 1.0;

  for (i = 0; i < cfg->n_shells; i++) {
    UnpackShell(cfg->shells+i, &n, &kl, &j, &nq);
    cfg->sweight *= ShellDegeneracy(j+1, nq);
    if (n == n0 && kl == kl0) {
      nq0 += nq;
    } else {
      if (nq0 > 0) {
	PackNRShell(cfg->nrs+m, n0, kl0, nq0);
	m++;
      }
      n0 = n;
      kl0 = kl;
      nq0 = nq;
    }
  }
  if (nq0 > 0) {
    PackNRShell(cfg->nrs+m, n0, kl0, nq0);
    m++;
  }
  
  cfg->nnrs = m;
  if (m < cfg->n_shells) {
    cfg->nrs = ReallocNew(cfg->nrs, sizeof(int)*m);
  }
  if (cfg->n_csfs > 0) {
    cfg->symstate = malloc(sizeof(int)*cfg->n_csfs);
  }
  cfg->igroup = k;
  cfg->icfg = cfg_groups[k].n_cfgs;
  CONFIG *acfg = ArrayAppend(clist, cfg, InitConfigData);
  if (acfg == NULL) return -1;
  if (cfg->n_csfs > 0) {    
    AddConfigToSymmetry(k, cfg_groups[k].n_cfgs, cfg); 
  }
  cfg_groups[k].n_cfgs++;
  if (cfg->shells[0].n > cfg_groups[k].nmax) {
    cfg_groups[k].nmax = cfg->shells[0].n;
  }
  int i0, i1, n1;
  n = acfg->shells[0].n;
  if (n > MAXN1) {
    i0 = MAXNK1-1;
  } else {
    i0 = ShellToInt(acfg->shells[0].n, acfg->shells[0].kappa);
  }
  if (acfg->shells[0].nq > 1) {
    i1 = i0;
  } else if (acfg->n_shells > 1 && acfg->shells[1].n <= MAXN2) {
    i1 = ShellToInt(acfg->shells[1].n, acfg->shells[1].kappa);
  } else {
    i1 = MAXNK2-1;
  }
  n1 = i0*MAXNK2 + i1;
  ARRAY *csa = &_csary[n1];
  ArrayAppend(csa, &acfg, InitPointerData);
  
  return 0;
}

/* 
** FUNCTION:    AddStateToSymmetry
** PURPOSE:     add a state to the symmetry list.
** INPUT:       {int kg, kc, kstate},
**              the group index, configuration index, and state index
**              of the state.
**              {int parity, j},
**              parity and total angular momentum of the state.
** RETURN:      {int},
**               0: success.
**              -1: error.
** SIDE EFFECT: 
** NOTE:        
*/
int AddStateToSymmetry(int kg, int kc, int kstate, int parity, int j) {
  int k;
  STATE s;
  ARRAY *st;
 
  k = IsEven(parity)? 2*j : (2*j+1);
  if (k >= MAX_SYMMETRIES) {
    printf("Maximum symmetry reached: %d %d\n", MAX_SYMMETRIES, k);
    int *ix = NULL;
    *ix=0;
    exit(1);
  }

  s.kgroup = kg;
  s.kcfg = kc;
  s.kstate = kstate;
  st = &(symmetry_list[k].states);
  if (ArrayAppend(st, &s, NULL) == NULL) return -1;
  symmetry_list[k].n_states++;
  return 0;
}

int ConfigParity(CONFIG *cfg) {
  int parity, i;

  parity = 0;
  for (i = 0; i < cfg->n_shells; i++) {
    parity += (cfg->shells[i].nq)*GetL(&(cfg->shells[i]));
  }
  parity /= 2;
  parity = IsOdd(parity);

  return parity;
}

int PackSymState(int s, int k) {
  return s*100000 + k;
}

void UnpackSymState(int st, int *s, int *k) {
  if (s) *s = st/100000;
  if (k) *k = st%100000;
}

/* 
** FUNCTION:    AddConfigToSymmetry
** PURPOSE:     add all states of a configuration to the symmetry list.
** INPUT:       {int kg, kc},
**              the group index and configuration index of the config.
**              a pointer of the configuration to be added.
** RETURN:      {int},
**               0: success.
**              -1: error.
** SIDE EFFECT: 
** NOTE:        
*/
int AddConfigToSymmetry(int kg, int kc, CONFIG *cfg) {
  int parity;
  int i, j, k, m;
  STATE s;
  ARRAY *st;

  parity = 0;
  for (i = 0; i < cfg->n_shells; i++) {
    parity += (cfg->shells[i].nq)*GetL(&(cfg->shells[i]));
  }
  parity /= 2;
  for (m = 0; m < cfg->n_csfs; m++) {
    i = m*cfg->n_shells;
    j = (cfg->csfs)[i].totalJ;
    k = IsEven(parity)? 2*j : (2*j+1);
    if (k >= MAX_SYMMETRIES) {
      printf("Maximum symmetry reached: %d %d\n", MAX_SYMMETRIES, k);
      int *ix = NULL;
      *ix=0;
      exit(1);
    }
    s.kgroup = kg;
    s.kcfg = kc;
    s.kstate = i;
    cfg->symstate[m] = PackSymState(k, symmetry_list[k].n_states);
    st = &(symmetry_list[k].states);
    if (ArrayAppend(st, &s, NULL) == NULL) return -1;
    symmetry_list[k].n_states++;
  }
  return 0;
}

/* 
** FUNCTION:    DecodePJ
** PURPOSE:     get the parity and J value from the symmetry index.
** INPUT:       {int i},
**              the symmetry index.
**              {int *p, int *j},
**              pointer holding the parity and J values.
** RETURN:      
** SIDE EFFECT: 
** NOTE:        if p or j is not required, pass in a NULL pointer.
*/
void DecodePJ(int i, int *p, int *j) {
  if (p) *p = IsOdd(i);
  if (j) *j = i/2;
}

/* 
** FUNCTION:    GetSymmetry
** PURPOSE:     return a pointer to the symmetry by its index.
** INPUT:       {int k},
**              the symmetry index.
** RETURN:      {SYMMETRY *},
**              pointer to the returned symmetry.
** SIDE EFFECT: 
** NOTE:        
*/
SYMMETRY *GetSymmetry(int k) {
  if (k < 0 || k >= MAX_SYMMETRIES) return NULL;
  return symmetry_list+k;
}

int ShellIndex(int n, int kappa, int ns, SHELL *s) {
  int i;

  for (i = 0; i < ns; i++) {
    if (s[i].n == n && s[i].kappa == kappa) return i;
  }
  return -1;
}

int ShellToInt(int n, int kappa) {  
  int k;
  k = (n-1)*(n-1) + 2*abs(kappa)- ((kappa>0)?1:2);
  return k;
}

void IntToShell(int i, int *n, int *kappa) {  
  int k;

  *n = ((int) sqrt(i)) + 1 ;
  k = i - ((*n)-1)*((*n)-1) + 1;
  if (IsOdd(k)) *kappa = -(k+1)/2;
  else *kappa = k/2;
}

int ConstructConfigName(char *s, int n, CONFIG *c) {
  int i, j, k, m;
  char a[16], b[32], x;

  m = 0;
  s[0] = '\0';
  for (i = c->n_shells-1; i >= 0; i--) {
    GetJLFromKappa(c->shells[i].kappa, &j, &k);
    SpecSymbol(a, k/2);
    if (j > k) x = '+';
    else x = '-';
    if (i > 0) {
      sprintf(b, "%d%s%c%d ", c->shells[i].n, a, x, c->shells[i].nq);
    } else {
      sprintf(b, "%d%s%c%d", c->shells[i].n, a, x, c->shells[i].nq);
    }
    m += strlen(b);
    if (m >= n) return -1;
    strcat(s, b);
  }
  return m;
}

void ListConfig(char *fn, int n, int *kg) {
  int i, m, j;
  CONFIG *c;
  CONFIG_GROUP *g;
  char a[2048];
  FILE *f;
  
  if (fn == NULL || strcmp(fn, "-") == 0) f = stdout;
  else f = fopen(fn, "w");

  m = 0;
  for (i = 0; i < n; i++) {
    g = GetGroup(kg[i]);
    for (j = 0; j < g->n_cfgs; j++) {
      c = GetConfigFromGroup(kg[i], j);
      ConstructConfigName(a, 2048, c);
      fprintf(f, "%32s %6d %6d %6d   %s\n",
	      g->name, kg[i], j, m, a);
      m++;
    }
  }
  
  if (f != stdout) fclose(f);
}

int ReadConfig(char *fn, char *c) {
  FILE *f;
  char buf[1024];
  char cbuf[1024];
  CONFIG *cfg;
  int t, j, ncfg;
  
  if (fn == NULL) return -1;
  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  while (1) {
    char *p = fgets(buf, 1024, f);
    if (p == NULL) break;
    p[32] = '\0';
    char *s = p;
    while(s && *s == ' ') s++;
    if (s) {
      if (c != NULL && strcmp(c, s)) continue;
      int t = GroupIndex(s);
      if (t < 0) return -1;
      char *c = &p[56];
      int i = 0;      
      while (c) {
	if (*c == '\n') {
	  break;
	}
	cbuf[i++] = *c;
	if (c[1] && c[2] && isalpha(c[1]) && (c[2] == '+' || c[2] == '-')) {
	  cbuf[i++] = '[';
	  cbuf[i++] = c[1];
	  cbuf[i++] = c[2];
	  cbuf[i++] = ']';
	  c += 2;
	}
	c++;
      }
      cbuf[i] = '\0';
      
      ncfg = GetConfigFromString(&cfg, cbuf);
      for (j = 0; j < ncfg; j++) {
	if (Couple(cfg+j) < 0) return -1;
	if (AddConfigToList(t, cfg+j) < 0) return -1;
      }
      if (ncfg > 0) free(cfg);
    }
  }
  fclose(f);
  return 0;
}

/* 
** FUNCTION:    GetAverageConfig
** PURPOSE:     determine the average configuration based on given
**              groups, the weight given to each group, and possible
**              screening orbitals.
** INPUT:       {int ng},
**              number of groups which determines the average config.
**              {int *kg},
**              ng elements array of groups indexes.
**              {double *weight},
**              weight given for each group.
**              {int n_screen},
**              number of screening orbitals.
**              {int *screened_n},
**              an array of principle quantum numbers for the 
**              screening orbitals.
**              {int screened_charge},
**              total charge to be screened by the screening orbitals.
**              {int screened_kl},
**              the orbital angular momentum used for the screening orbital.
**              -1: use the kl = 0 orbital.
**               0: use the kl = n/2 orbital.
**              +1: use the kl = n-1 orbital.
**              {AVERAGE_CONFIG *acfg},
**              pointer holding the resulting average configuration.
** RETURN:      {int},
**              >=0: success, the number of shells in the average config.
**               -1: error.
** SIDE EFFECT: 
** NOTE:        there is a limit the highest n the shells in the average
**              configuration can take. it is determined by the macro M.
**              with M = 2500, the limit is about 70, which should be 
**              more than enough.
*/
int GetAverageConfig(int ng, int *kg, int ic, double *weight,
		     int n_screen, int *screened_n, double screened_charge,
		     int screened_kl, AVERAGE_CONFIG *acfg) {
#define M 2500 /* max # of shells may be present in an average config */

  double tnq[M];
  int i, j, k, n, kappa, t;
  ARRAY *c;
  CONFIG *cfg;
  double a;

  if (ng <= 0) return -1;
  for(i = 0; i < M; i++) tnq[i] = 0.0;

  acfg->ng = ng;
  acfg->kg = malloc(sizeof(int)*ng);
  acfg->weight = malloc(sizeof(double)*ng);
  for (i = 0; i < ng; i++) {
    acfg->kg[i] = kg[i];
    if (weight != NULL) {
      acfg->weight[i] = weight[i];
    } else {
      acfg->weight[i] = 1.0;
    }
  }

  /* normalize the weight */
  a = 0.0;
  for (i = 0; i < ng; i++) {
    a += acfg->weight[i];
  }
  for (i = 0; i < ng; i++) {
    acfg->weight[i] /= a;
  }

  for (i = 0; i < ng; i++) {
    c = &(cfg_groups[kg[i]].cfg_list);
    if (ic < 0) {
      a = 1.0/cfg_groups[kg[i]].n_cfgs;
    } else {
      a = 1.0;
    }
    for (t = 0; t < cfg_groups[kg[i]].n_cfgs; t++) {
      if (ic >= 0 && t != ic) continue;
      cfg = (CONFIG *) ArrayGet(c, t);
      for (j = 0; j < cfg->n_shells; j++) {
	n = cfg->shells[j].n;
	kappa = cfg->shells[j].kappa;
	k = ShellToInt(n, kappa);
	if (k >= M) k = M-1;
	tnq[k] += (((double)(cfg->shells[j].nq)) * acfg->weight[i]*a);
      }
    }
    acfg->n_cfgs += cfg_groups[kg[i]].n_cfgs;
  }

  for (i = 0, j = 0; i < M; i++) {
    if (tnq[i] > EPS10) {
      j++;
    }
  }
  
  acfg->n_shells = j;
  acfg->n = malloc(sizeof(int)*j);
  acfg->kappa = malloc(sizeof(int)*j);
  acfg->nq = malloc(sizeof(double)*j);
  
  if (!acfg->n ||
      !acfg->kappa ||
      !acfg->nq) goto ERROR;

  for (i = 0, j = 0; i < M; i++) {
    if (tnq[i] > EPS10) {
      IntToShell(i, &n, &kappa);
      acfg->n[j] = n;
      acfg->kappa[j] = kappa;
      acfg->nq[j] = tnq[i];
#if FAC_DEBUG      
      fprintf(debug_log, "acfg: %d %d %lf\n", n, kappa, tnq[i]);
#endif
      j++;
    }
  }

  /* add in configs for screened_charge */
  if (n_screen > 0) {
    screened_charge /= (double) n_screen;
    for (i = 0; i < n_screen; i++) {
      if (screened_kl < 0) {
	t = 0;
	kappa = -1;
      } else if (screened_kl == 0) {
	t = screened_n[i];
	kappa = GetKappaFromJL(t+1, t);
      } else {
	t = screened_n[i]*2-2;
	kappa = GetKappaFromJL(t+1, t);
      }    
      for (j = 0; j < acfg->n_shells; j++) {
	k = GetLFromKappa(acfg->kappa[j]);
	if (acfg->n[j] < screened_n[i]) continue;
	if (acfg->n[j] > screened_n[i]) break;
	if (k > t) break;
	if (acfg->kappa[j] == kappa) break;
      }
      if (j < acfg->n_shells && 
	  acfg->n[j] == screened_n[i] && 
	  acfg->kappa[j] == kappa) {
	acfg->nq[j] += screened_charge; 
      } else {
	acfg->n_shells += 1;
	acfg->n = realloc(acfg->n, sizeof(int)*acfg->n_shells);
	acfg->kappa = realloc(acfg->kappa, sizeof(int)*acfg->n_shells);
	acfg->nq = realloc(acfg->nq, sizeof(double)*acfg->n_shells);
	for (k = acfg->n_shells-1; k > j; k--) {
	  acfg->n[k] = acfg->n[k-1];
	  acfg->kappa[k] = acfg->kappa[k-1];
	  acfg->nq[k] = acfg->nq[k-1];
	}
	acfg->n[j] = screened_n[i];
	acfg->kappa[j] = kappa;
	acfg->nq[j] = screened_charge;
      }
    }
  }
	  
  return j;

 ERROR:
  if (acfg->n) free(acfg->n);
  if (acfg->nq) free(acfg->nq);
  if (acfg->kappa) free(acfg->kappa);
  free(acfg->kg);
  free(acfg->weight);
  return -1;

#undef M
}

/* 
** FUNCTION:    InGroups
** PURPOSE:     determing if a group is within a list of groups.
** INPUT:       {int kg},
**              the index of the group to be tested.
**              {int ng, *kgroup}
**              the number and index of the groups in the list.
** RETURN:      {int},
**              0: group kg is not in the list.
**              1: group kg is in the list.
** SIDE EFFECT: 
** NOTE:        
*/
int InGroups(int kg, int ng, int *kgroup) {
  int i;
  if (ng < 0) {
    if (kg >= 0 && kg < n_groups) return 1;
  }

  for (i = 0; i < ng; i++) {
    if (kg == kgroup[i]) return 1;
  }

  return 0;
}

/* 
** FUNCTION:    CompareShell
** PURPOSE:     determine which of the two shells is the inner one.
** INPUT:       {const void *s1, *s2},
**              two shells in comparison.
** RETURN:      {int},
**              -1: s1 is inside s2.
**               0: s1 and s2 are the same.
**              +1: s1 is outside s2.
** SIDE EFFECT: 
** NOTE:        
*/
int CompareShell(const void *ts1, const void *ts2) {
  SHELL *s1, *s2;
  int ak1, ak2;

  s1 = (SHELL *) ts1;
  s2 = (SHELL *) ts2;
  if (s1->n > s2->n) return 1;
  else if (s1->n < s2->n) return -1;
  else {
    if (s1->kappa == s2->kappa) return 0;
    else {
      ak1 = abs(s1->kappa);
      ak2 = abs(s2->kappa);
      if (ak1 > ak2) return 1;
      else if (ak1 < ak2) return -1;
      else {
	if (s1->kappa < s2->kappa) return -1;
	else return 1;
      }
    }
  }
}

int CompareShellInvert(const void *ts1, const void *ts2) {
  return -CompareShell(ts1, ts2);
}

/* 
** FUNCTION:    InitConfig
** PURPOSE:     initialize the module "config".
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int InitConfig(void) {
  int i;

  n_groups = 0;
  cfg_groups = malloc(MAX_GROUPS*sizeof(CONFIG_GROUP));
  for (i = 0; i < MAX_GROUPS; i++) {
    strcpy(cfg_groups[i].name, "_all_");
    cfg_groups[i].n_cfgs = 0;
    ArrayInit(&(cfg_groups[i].cfg_list), sizeof(CONFIG), CONFIGS_BLOCK);
  }

  symmetry_list = malloc(MAX_SYMMETRIES*sizeof(SYMMETRY));
  for (i = 0; i < MAX_SYMMETRIES; i++) {
    symmetry_list[i].n_states = 0;
    ArrayInit(&(symmetry_list[i].states), sizeof(STATE), STATES_BLOCK);
  }
  int i1, i2, i3, i4;
  for (i = 0; i < NJQ; i++) {
    _csf1[i] = NULL;
    _ncsf1[i] = 0;
    for (i1 = 0; i1 < NJQ; i1++) {
      _csf2[i][i1] = NULL;
      _ncsf2[i][i1] = 0;
      for (i2 = 0; i2 < NJQ; i2++) {
	_csf3[i][i1][i2] = NULL;
	_ncsf3[i][i1][i2] = 0;
	for (i3 = 0; i3 < NJQ; i3++) {
	  _csf4[i][i1][i2][i3] = NULL;
	  _ncsf4[i][i1][i2][i3] = 0;
	  for (i4 = 0; i4 < NJQ; i4++) {
	    _csf5[i][i1][i2][i3][i4] = NULL;
	    _ncsf5[i][i1][i2][i3][i4] = 0;
	  }
	}
      }
    }
  }
  for (i = 0; i < MAXNK12; i++) {
    ArrayInit(&_csary[i], sizeof(CONFIG *), CONFIGS_BLOCK);
  }

  SetClosedShellNR(0, 0);
  return 0; 
}

/* 
** FUNCTION:    FreeConfigData
** PURPOSE:     free the memory in a CONFIG struct.
** INPUT:       {void *},
**              pointer to the CONFIG struct.
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
static void FreeConfigData(void *p) {
  CONFIG *c;

  c = (CONFIG *) p;
  if (c->n_shells > 0) {
    free(c->shells);
    c->n_shells = 0;
  }
  if (c->n_csfs > 0) {
    free(c->symstate);
    free(c->csfs);
    c->n_csfs = 0;
  }
  if (c->nnrs > 0) {
    free(c->nrs);
    c->nnrs = 0;
  }
}

int RemoveGroup(int k) {
  SYMMETRY *sym;
  STATE *s;
  int i, m;

  if (k != n_groups-1) {
    printf("only the last group can be removed\n");
    return -1;
  }
  ArrayFree(&(cfg_groups[k].cfg_list), FreeConfigData);
  cfg_groups[k].n_cfgs = 0;
  strcpy(cfg_groups[k].name, "_all_");
  n_groups--;

  for (i = 0; i < MAX_SYMMETRIES; i++) {
    sym = GetSymmetry(i);
    for (m = 0; m < sym->n_states; m++) {
      s = ArrayGet(&(sym->states), m);
      if (s->kgroup == k) {
	ArrayTrim(&(sym->states), m, NULL);
	sym->n_states = m;
	break;
      }
    }
  }
  
  return 0;
}
  
/* 
** FUNCTION:    ReinitConfig
** PURPOSE:     reinitialize the module "config".
** INPUT:       {int m},
**              0: do a full reinitialization.
**              -1, 1: do nothing.
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int ReinitConfig(int m) {
  int i;

  if (m) return 0;

  for (i = 0; i < n_groups; i++) {
    ArrayFree(&(cfg_groups[i].cfg_list), FreeConfigData);
    cfg_groups[i].n_cfgs = 0;
    strcpy(cfg_groups[i].name, "_all_");
  }
  n_groups = 0;

  for (i = 0; i < MAX_SYMMETRIES; i++) {
    if (symmetry_list[i].n_states > 0) {
      ArrayFree(&(symmetry_list[i].states), NULL);
      symmetry_list[i].n_states = 0;
    }
  }
  for (i = 0; i < MAXNK12; i++) {
    ArrayFree(&_csary[i], NULL);
  }

  SetClosedShellNR(0, 0);
  return 0;
}
