#include "config.h"

static char *rcsid="$Id: config.c,v 1.24 2003/04/21 02:21:35 mfgu Exp $";
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

static int nstates_partition = NSPARTITION;

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

int SetNStatesPartition(int n) {
  if (n > 0) {
    nstates_partition = n;
  } else {
    nstates_partition = NSPARTITION;
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
  char token[128];
  int r, brkpos, quotepos, next;
  int nn, nkl, nkappa;
  int n[16];
  int kl[512];
  int kappa[1024];
  int ncfg;
  double dnq;
  int *maxq;
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
      if (scfg[next] != '+' && kl2 > 0) {
	kappa[nkappa++] = GetKappaFromJL(kl2-1, kl2);
      }
      if (scfg[next] != '-') {
	kappa[nkappa++] = GetKappaFromJL(kl2+1, kl2);
      }
    }
    if (scfg[next] == '+' || scfg[next] == '-') next++;
    dnq = atof(&(scfg[next]));
    if (dnq == 0) dnq = 1;
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
    dnq = atof(&(scfg[next]));
    if (dnq == 0) dnq = 1;
  } else {
    return -1;
  }

  shell = (SHELL *) malloc(sizeof(SHELL)*nn*nkappa);
  t = 0;
  for (i = nn-1; i >= 0; i--) {
    for (k = nkappa-1; k >= 0; k--) {
      kl2 = GetLFromKappa(kappa[k]);
      if (kl2/2 >= n[i]) continue;
      shell[t].n = n[i];
      shell[t].kappa = kappa[k];
      t++;
    }
  }
  ns = t;
  
  if (ns == 0) {
    free(shell);
    return -1;
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

  ncfg = DistributeElectronsShell(cfg, ns, shell, (int)dnq, maxq);

  free(shell);
  free(maxq);

  return ncfg;
}

int DistributeElectronsNR(CONFIG **cfg, char *scfg) {
  SHELL *shell;
  char token[128];
  int r, brkpos, quotepos, next;
  int nn, nkl, nkappa;
  int n[16];
  int kl[512];
  int kappa[1024];
  int ncfg;
  double dnq;
  int *maxq;
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
    dnq = atof(&(scfg[next]));
    if (dnq == 0) dnq = 1;
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
    dnq = atof(&(scfg[next]));
    if (dnq == 0) dnq = 1;
  } else {
    return -1;
  }

  shell = (SHELL *) malloc(sizeof(SHELL)*nn*nkappa);
  t = 0;
  for (i = nn-1; i >= 0; i--) {
    for (k = nkappa-1; k >= 0; k--) {
      kl2 = kappa[k];
      if (kl2/2 >= n[i]) continue;
      shell[t].n = n[i];
      shell[t].kappa = kappa[k];
      t++;
    }
  }
  ns = t;
  
  if (ns == 0) {
    free(shell);
    return -1;
  }

  maxq = (int *) malloc(sizeof(int)*ns);
  maxq[ns-1] = 0;
  for (i = ns-2; i >= 0; i--) {
    maxq[i] = maxq[i+1] + 2*(shell[i+1].kappa + 1);
  }

  ncfg = DistributeElectronsShellNR(cfg, ns, shell, (int)dnq, maxq);

  free(shell);
  free(maxq);

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
  double *dnq, *p2;
  char *s;
  int ncfg, *dnc;  
  int size, size_old, tmp;
  int i, t, j, k, ns;

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
    dnq = (double *) malloc(sizeof(double *)*ns);
    p2 = dnq;
  } else {
    dnq = NULL;
    p2 = NULL;
  }

  s = scfg;
  p1 = dcfg;
  for (i = 0; i < ns; i++) {
    while (*s == ' ' || *s == '\t') s++;
    dnc[i] = DistributeElectrons(p1, p2, s);
    if (dnc[i] <= 0) {
      return -1;
    }
    while (*s) s++;
    s++;
    p1++;
    if (p2) p2++;
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
      for (j = 0; j < dnc[i]; j++) {
	(*nq)[t] = dnq[i]/dnc[i];
	memcpy((*cfg)->shells+t, (*p1)->shells+j, sizeof(SHELL));
	t++;
      }
      p1--;
    }
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

int GetConfigFromStringNR(CONFIG **cfg, char *scfg) {
  CONFIG **dcfg, **p1;
  double *dnq;
  char *s;
  int ncfg, *dnc;  
  int size, size_old, tmp;
  int i, t, j, k, ns;

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
  p1 = dcfg;
  for (i = 0; i < ns; i++) {
    while (*s == ' ' || *s == '\t') s++;
    dnc[i] = DistributeElectronsNR(p1, s);
    if (dnc[i] <= 0) {
      return -1;
    }
    while (*s) s++;
    s++;
    p1++;
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
**              a list of principal quantum numbers, angular 
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
  if (*p == '+') {
    if (*j) *j = 1;
    *p = '\0';
  } else if (*p == '-') {
    if (*j) *j = -1;
    *p = '\0';
  } else {
    if (*j) *j = 0;
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
  int errcode;

  if (cfg->n_shells == 0) {
    errcode = -1;
    goto ERROR;
  }

  if (cfg == NULL) {
    errcode = -2;
    goto ERROR; 
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

  return 0;

 ERROR:
  printf("****Error in Couple****\n");
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

PARTITION *GetPartition(int kg, int kp) {
  return (PARTITION *) ArrayGet(&(cfg_groups[kg].partition), kp);
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
  ARRAY *part;
  PARTITION p, *t;
  int ip, ns;

  if (cfg->n_csfs > nstates_partition) {
    printf("Error: a single configuration in group %s ", cfg_groups[k].name);
    printf("has more than %d states\n", nstates_partition);
    printf("Enlarge nstates_partition to at least %d\n", cfg->n_csfs);
    return -1;
  }
  if (k < 0 || k >= n_groups) return -1;
  if (cfg_groups[k].n_cfgs == 0) {
    cfg_groups[k].n_electrons = cfg->n_electrons;
  } else if (cfg_groups[k].n_electrons != cfg->n_electrons) {
    printf("Error: AddConfigToList, Configurations in a group ");
    printf("must have the same number of electrons\n");
    return -1;
  }
  clist = &(cfg_groups[k].cfg_list);
  part = &(cfg_groups[k].partition);

  cfg->energy = 0.0;
  cfg->delta = 0.0;

  if (part->dim == 0) {
    cfg->ipart = 0;
    cfg->npart = 0;
    p.icfg1 = clist->dim;
    p.icfg2 = p.icfg1;
    p.n_csfs = cfg->n_csfs;
    ArrayAppend(part, &p);
  } else {
    ip = part->dim - 1;
    t = ArrayGet(part, ip);
    ns = t->n_csfs + cfg->n_csfs;
    if (ns <= nstates_partition) {
      cfg->ipart = ip;
      cfg->npart = t->n_csfs;
      t->icfg2++;
      t->n_csfs = ns;
    } else {
      ip++;
      cfg->ipart = ip;
      cfg->npart = 0;
      p.icfg1 = clist->dim;
      p.icfg2 = p.icfg1;
      p.n_csfs = cfg->n_csfs;
      ArrayAppend(part, &p);
    }
  }

  if (ArrayAppend(clist, cfg) == NULL) return -1;
  AddConfigToSymmetry(k, cfg_groups[k].n_cfgs, cfg); 
  cfg_groups[k].n_cfgs++;
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
    printf("MAX_SYMMETRIES reached\n");
    exit(1);
  }

  s.kgroup = kg;
  s.kcfg = kc;
  s.kstate = kstate;
  st = &(symmetry_list[k].states);
  if (ArrayAppend(st, &s) == NULL) return -1;
  symmetry_list[k].n_states++;
  return 0;
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
  int i, j, k;
  STATE s;
  ARRAY *st;

  parity = 0;
  for (i = 0; i < cfg->n_shells; i++) {
    parity += (cfg->shells[i].nq)*GetL(&(cfg->shells[i]));
  }
  parity /= 2;
  for (i = 0; i < (cfg->n_csfs)*(cfg->n_shells); i += cfg->n_shells) {
    j = (cfg->csfs)[i].totalJ;
    k = IsEven(parity)? 2*j : (2*j+1);
    if (k >= MAX_SYMMETRIES) {
      printf("MAX_SYMMETRIES reached\n");
      exit(1);
    }

    s.kgroup = kg;
    s.kcfg = kc;
    s.kstate = i;
    st = &(symmetry_list[k].states);
    if (ArrayAppend(st, &s) == NULL) return -1;
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
**              an array of principal quantum numbers for the 
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
int GetAverageConfig(int ng, int *kg, double *weight,
		     int n_screen, int *screened_n, double screened_charge,
		     int screened_kl, AVERAGE_CONFIG *acfg) {
#define M 2500 /* max # of shells may be present in an average config */

  double tnq[M];
  int i, j, k, n, kappa, t;
  ARRAY *c;
  CONFIG *cfg;
  int weight_allocated = 0;
  double a;

  if (ng <= 0) return -1;
  for(i = 0; i < M; i++) tnq[i] = 0.0;

  if (weight == NULL) {
    weight = malloc(sizeof(double)*ng);
    if (!weight) return -1;
    for (i = 0; i < ng; i++) {
      weight[i] = 1.0;
    }
    weight_allocated = 1;
  }

  /* normalize the weight */
  a = 0.0;
  for (i = 0; i < ng; i++) {
    a += weight[i];
  }
  for (i = 0; i < ng; i++) {
    weight[i] /= a;
  }

  for (i = 0; i < ng; i++) {
    c = &(cfg_groups[kg[i]].cfg_list);
    a = 1.0/cfg_groups[kg[i]].n_cfgs;
    for (t = 0; t < cfg_groups[kg[i]].n_cfgs; t++) {
      cfg = (CONFIG *) ArrayGet(c, t);
      for (j = 0; j < cfg->n_shells; j++) {
	n = cfg->shells[j].n;
	kappa = cfg->shells[j].kappa;
	k = (n-1)*(n-1) + 2*abs(kappa)- ((kappa>0)?1:2);
	if (k >= M) k = M-1;
	tnq[k] += (((double)(cfg->shells[j].nq)) * weight[i]*a);
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
      n = ((int) sqrt(i)) + 1 ;
      k = i - (n-1)*(n-1) + 1;
      if (IsOdd(k)) kappa = -(k+1)/2;
      else kappa = k/2;
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
	  
  if (weight_allocated) {
    free(weight);
    weight = NULL;
  }

  return j;

 ERROR:
  if (acfg->n) free(acfg->n);
  if (acfg->nq) free(acfg->nq);
  if (acfg->kappa) free(acfg->kappa);
  if (weight_allocated) free(weight);
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
  if (ng == 0) {
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
** INPUT:       {SHELL *s1, *s2},
**              two shells in comparison.
** RETURN:      {int},
**              -1: s1 is inside s2.
**               0: s1 and s2 are the same.
**              +1: s1 is outside s2.
** SIDE EFFECT: 
** NOTE:        
*/
int CompareShell(SHELL *s1, SHELL *s2) {
  int ak1, ak2;
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
    ArrayInit(&(cfg_groups[i].partition), sizeof(PARTITION), CONFIGS_BLOCK);
    ArrayInit(&(cfg_groups[i].cfg_list), sizeof(CONFIG), CONFIGS_BLOCK);
  }

  symmetry_list = malloc(MAX_SYMMETRIES*sizeof(SYMMETRY));
  for (i = 0; i < MAX_SYMMETRIES; i++) {
    symmetry_list[i].n_states = 0;
    ArrayInit(&(symmetry_list[i].states), sizeof(STATE), STATES_BLOCK);
  }
  
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
void FreeConfigData(void *p) {
  CONFIG *c;

  c = (CONFIG *) p;
  if (c->n_shells > 0) {
    free(c->shells);
    c->n_shells = 0;
  }
  if (c->n_csfs > 0) {
    free(c->csfs);
    c->n_csfs = 0;
  }
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
    ArrayFree(&(cfg_groups[i].partition), NULL);
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

  return 0;
}
