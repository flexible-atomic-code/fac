#include "config.h"

static CONFIG_GROUP *cfg_groups;
static int n_groups; 

static SYMMETRY *symmetry_list;

static char spec_symbols[MAX_SPEC_SYMBOLS] = 
{'s', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'o', 'p', 'q'};

int SpecSymbol(char *s, int kl) {
  if (kl < MAX_SPEC_SYMBOLS) {
    s[0] = spec_symbols[kl];
    s[1] = '\0';
  } else {
    sprintf(s, "[%d]", kl);
  }
  return 0;
}

/** This function recursively construct all possible states for a Config. **/
int Couple(CONFIG *cfg) {
  CONFIG outmost, inner;
  SHELL *shells;
  int errcode;
  int i;
  int bytes_csf;

  if (cfg->n_shells == 0) {
    errcode = -1;
    goto ERROR;
  }

  if (cfg == NULL) {
    errcode = -2;
    goto ERROR; 
  }

  bytes_csf = cfg->n_shells * sizeof(SHELL_STATE);

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


/**************************************************************************
  This function constructs all possible states by coupling the outmost shell
  to the inner shells. both outmost shell and inner shells must have been
  coupled.
***************************************************************************/
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
  
  return 0;
  
 ERROR:
  printf("****Error in CoupleOutmost****\n");
  return -1;
}


/**************************************************************************
  GetSingleShell construct all possible states for a single shell.
  for j > 9/2, no more than 2 electrons are allowed. The data were taken
  from "Nuclear Shell Theory" by AMOS de-SHALIT and IGAL TALMI.
***************************************************************************/
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

int ShellClosed(SHELL *s) {
  int j;
  j = GetJ(s);
  if (s->nq < j+1) return 0;
  return 1;
}

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

void PackShellState(SHELL_STATE *s, int J, int j, int nu, int Nr){
  s->totalJ = J;
  s->shellJ = j;
  s->nu = nu;
  s->Nr = Nr;
}

int GroupIndex(char *name) {
  int i, j;

  for (i = n_groups - 1; i >= 0; i--) {
    if (strncmp(name, cfg_groups[i].name, GROUP_NAME_LEN) == 0) 
      break;
  }
  if (i < 0) i = AddGroup(name);
  return i;
}

int GroupExists(char *name) {
  int i;

  for (i = n_groups - 1; i >= 0; i--) {
    if (strncmp(name, cfg_groups[i].name, GROUP_NAME_LEN) == 0) 
      break;
  }
  return i;
}

int AddGroup(char *name) {
  int j;
  if (name == NULL) return -1;
  if (n_groups == MAX_GROUPS) {
    printf("Max # groups reached\n");
    abort();
  }
  strncpy(cfg_groups[n_groups].name, name, GROUP_NAME_LEN);
  n_groups++;
  return n_groups-1;
}

CONFIG_GROUP *GetGroup(int k) {
  if (k < 0 || k >= n_groups) return NULL;
  return cfg_groups+k;
}

CONFIG_GROUP *GetNewGroup() {
  int j;
  if (n_groups == MAX_GROUPS) {
    printf("Max # groups reached\n");
    abort();
  }
  n_groups++;
  return cfg_groups+n_groups-1;
}

int GetNumGroups() {
  return n_groups;
}

CONFIG *GetConfig(STATE *s) {
  CONFIG *c;
  int i, j;

  i = s->kgroup;
  j = s->kcfg;
  c = (CONFIG *) ArrayGet(&(cfg_groups[i].cfg_list), j);
  return c;
}

int AddConfigToList(int k, CONFIG *cfg) {
  ARRAY *clist;

  if (k < 0 || k >= n_groups) return -1;
  clist = &(cfg_groups[k].cfg_list);

  if (ArrayAppend(clist, cfg) == NULL) return -1;
  
  AddConfigToSymmetry(k, cfg_groups[k].n_cfgs, cfg); 
  cfg_groups[k].n_cfgs++;
  return 0;
}

int AddStateToSymmetry(int kg, int kc, int kstate, int parity, int j) {
  int k;
  STATE s;
  ARRAY *st;
 
  k = IsEven(parity)? 2*j : (2*j+1);
  if (k >= MAX_SYMMETRIES) {
    printf("MAX_SYMMETRIES reached\n");
    abort();
  }

  s.kgroup = kg;
  s.kcfg = kc;
  s.kstate = kstate;
  st = &(symmetry_list[k].states);
  if (ArrayAppend(st, &s) == NULL) return -1;
  symmetry_list[k].n_states++;
  return 0;
}

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
      abort();
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

void DecodePJ(int i, int *p, int *j) {
  if (p) *p = IsOdd(i);
  if (j) *j = i/2;
}

SYMMETRY *GetSymmetry(int k) {
  if (k < 0 || k >= MAX_SYMMETRIES) return NULL;
  return symmetry_list+k;
}

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
  int n_states_cfg, n_states_group;

  if (ng <= 0) return -1;
  for(i = 0; i < M; i++) tnq[i] = 0.0;

  acfg->n_cfgs = 0;
  acfg->n = NULL;
  acfg->nq = NULL;
  acfg->kappa = NULL;

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
    n_states_group = 0;
    for (t = 0; t < cfg_groups[kg[i]].n_cfgs; t++) {
      cfg = (CONFIG *) ArrayGet(c, t);
      n_states_cfg = 0;
      for (j = 0; j < cfg->n_csfs; j += cfg->n_shells) {
	n_states_cfg += cfg->csfs[j].totalJ + 1;
      }
      n_states_group += n_states_cfg;
    }
    for (t = 0; t < cfg_groups[kg[i]].n_cfgs; t++) {
      cfg = (CONFIG *) ArrayGet(c, t);
      n_states_cfg = 0;
      for (j = 0; j < cfg->n_csfs; j += cfg->n_shells) {
	n_states_cfg += cfg->csfs[j].totalJ + 1;
      }
      for (j = 0; j < cfg->n_shells; j++) {
	n = cfg->shells[j].n;
	kappa = cfg->shells[j].kappa;
	k = (n-1)*(n-1) + 2*abs(kappa)- ((kappa>0)?1:2);
	if (k >= M) k = M-1;
	tnq[k] += (((double)(cfg->shells[j].nq)) * weight[i] *
		   ((double)n_states_cfg)/n_states_group);
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
	kappa = GetKappaFromJL(t, t+1);
      } else {
	t = screened_n[i]*2-2;
	kappa = GetKappaFromJL(t, t+1);
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


int InitConfig() {
  int i, n, m;

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
  
  return 0; 
}


