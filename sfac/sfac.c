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

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "init.h"
#include "stoken.h"
#include "mpiutil.h"

static int PPrint(int argc, char *argv[], int argt[], ARRAY *variables) {
  int i;
  for (i = 0; i < argc; i++) {
    switch (argt[i]) {
    case NUMBER:
      printf("%s", argv[i]);
      break;
    case STRING:
      printf("%s", argv[i]);
      break;
    case LIST:
      printf("[%s]", argv[i]);
      break;
    case TUPLE:
      printf("(%s)", argv[i]);
      break;
    case KEYWORD:
      printf("%s = ", argv[i]);
    }
    if (i != argc-1 && argt[i] != KEYWORD) {
      printf(" ");
    }
  }
  if (argc > 0) printf("\n");
  fflush(stdout);
  return 0;
}

static int IntFromList(char *argv, int tp, ARRAY *variables, int **k) {
  int i;
  char *v[MAXNARGS];
  int t[MAXNARGS], n;
  
  if (tp == LIST) {
    n = DecodeArgs(argv, v, t, variables);
    if (n > 0) {
      *k = malloc(sizeof(int)*n);
      for (i = 0; i < n; i++) {
	(*k)[i] = atoi(v[i]);
	free(v[i]);
      }
    }
  } else {
    n = 1;
    *k = malloc(sizeof(int));
    (*k)[0] = atoi(argv);
  }
  
  return n;
}

static int DoubleFromList(char *argv, int tp, ARRAY *variables, double **k) {
  int i;
  char *v[MAXNARGS];
  int t[MAXNARGS], n;
  
  if (tp == LIST) {
    n = DecodeArgs(argv, v, t, variables);
    if (n > 0) {
      *k = malloc(sizeof(double)*n);
      for (i = 0; i < n; i++) {
	(*k)[i] = atof(v[i]);
	free(v[i]);
      }
    }
  } else {
    n = 1;
    *k = malloc(sizeof(double));
    (*k)[0] = atof(argv);
  }
  
  return n;
}

static int DecodeGroupArgs(int **kg, int n, int *n0, char *argv[], int argt[],
			   ARRAY *variables) {
  char *s;
  int i, k, ng;
  char *v[MAXNARGS];
  int t[MAXNARGS], nv;

  ng = n;
  nv = 0;
  if (ng > 0) {
    if (argt[0] == LIST || argt[0] == TUPLE) {
      if (ng > 1) {
	printf("there should be only one list or tuple\n");
	return -1;
      }
      ng = DecodeArgs(argv[0], v, t, variables);
      nv = ng;
    } else {
      for (i = 0; i < ng; i++) {
	v[i] = argv[i];
	t[i] = argt[i];
      }
    }
    (*kg) = malloc(sizeof(int)*ng);
    n = 0;    
    int n0p, n0q;
    if (n0) {
      n0p = *n0;
      n0q = *n0;
    } else {
      n0p = 0;
      n0q = 0;
    }
    for (i = 0; i < ng; i++) {
      if (t[i] != STRING) {
	printf("argument must be a group name\n");
	free((*kg));
	return -1;
      }
      s = v[i];
      k = GroupExists(s);      
      if (k < 0) {
	printf("group does not exist: %d %s\n", i, s);
	if (i < n0q) n0p--;
	continue;
      }
      (*kg)[n] = k;
      n++;
    }
    if (n0) *n0 = n0p;
    ng = n;
    if (ng <= 0) {
      printf("all cfg groups invalid\n");
      free(*kg);
      return ng;
    }
  } else {
    ng = GetNumGroups();
    (*kg) = malloc(sizeof(int)*ng);
    for (i = 0; i < ng; i++) (*kg)[i] = i;
  }

  for (i = 0; i < nv; i++) free(v[i]);
  
  return ng;
}

static int SelectLevels(int **t, char *argv, int argt, ARRAY *variables) {
  int n, ng, *kg, i, j, k, im, m, m0;
  int nrg, *krg, nrec;
  int ig, nlevels, iuta;
  LEVEL *lev;
  SYMMETRY *sym;
  STATE *s;
  char rgn[GROUP_NAME_LEN];
  char *v[MAXNARGS], *v1[MAXNARGS];
  int at[MAXNARGS], at1[MAXNARGS], nv, nv1, rv;

  if (argt != LIST  && argt != TUPLE) return -1;
  nv = 0; 
  nv1 = 0;
  rv = 0;

  iuta = IsUTA();
  n = DecodeArgs(argv, v, at, variables);
  nv = n;
  if (n > 0) {
    if (at[0] == STRING) {
      ng = DecodeGroupArgs(&kg, n, NULL, v, at, variables);
      if (ng <= 0) {
	rv = -1;
	goto END;
      }
      nlevels = GetNumLevels();
      (*t) = malloc(sizeof(int)*nlevels);
      k = 0;
      for (j = 0; j < nlevels; j++) {
	lev = GetLevel(j);
	if (iuta) {
	  ig = lev->iham;
	} else {
	  im = lev->pb;
	  sym = GetSymmetry(lev->pj);
	  s = (STATE *) ArrayGet(&(sym->states), im);
	  ig = s->kgroup;
	}
	if (InGroups(ig, ng, kg)) {
	  (*t)[k] = j;
	  k++;
	}
      }
      free(kg);
      (*t) = realloc(*t, k*sizeof(int));
      rv = k;
      goto END;
    } else if (at[0] == LIST) {
      if (n != 2) {
	printf("recombined states specification unrecoganized\n");
	rv = -1;
	goto END;
      }
      ng = DecodeGroupArgs(&kg, 1, NULL, v, at, variables);
      if (ng <= 0) {
	rv = -1;
	goto END;
      }
      if (at[1] == LIST) {
	m0 = 0;
	n = DecodeArgs(v[1], v1, at1, variables);
	nv1 = n;
      } else if (at[1] == NUMBER) {
	m0 = 1;
	v1[1] = v[1];
	at1[1] = at[1];
      } else {
	printf("Level specification unrecoganized\n");
	rv = -1;
	goto END;
      }
    
      nrg = ng;
      krg = malloc(sizeof(int)*nrg);
      nlevels = GetNumLevels();
      (*t) = malloc(sizeof(int)*nlevels);
      if (!(*t)) {
	rv = -1;
	goto END;
      }
      k = 0;
      for (m = m0; m < n; m++) {
	if (at1[m] != NUMBER) {
	  rv = -1;
	  goto END;
	}
	nrec = atoi(v1[m]);
	for (i = 0; i < nrg; i++) {
	  ConstructRecGroupName(rgn, GetGroup(kg[i])->name, nrec);
	  krg[i] = GroupExists(rgn);
	}
	for (j = 0; j < nlevels; j++) {
	  lev = GetLevel(j);
	  im = lev->pb;
	  sym = GetSymmetry(lev->pj);
	  s = (STATE *) ArrayGet(&(sym->states), im);
	  ig = s->kgroup;
	  if (ig < 0) {
	    if (!ValidBasis(s, ng, kg, nrec)) continue;
	    (*t)[k] = j;
	    k++;
	  } else {
	    if (InGroups(ig, nrg, krg)) {
	      (*t)[k] = j;
	      k++;
	    }
	  }
	}
      }
      
      free(krg);
      free(kg);
      (*t) = realloc(*t, k*sizeof(int));
      rv = k;
      goto END;
    } else {
      (*t) = malloc(sizeof(int)*n);
      for (i = 0; i < n; i++) {
	if (at[i] != NUMBER) return -1;
	(*t)[i] = atoi(v[i]);
      }
      rv = n;
      goto END;
    }
  }

 END:
  for (i = 0; i < nv; i++) free(v[i]);
  for (i = 0; i < nv1; i++) free(v1[i]);
  return rv;
}

static int ConfigListToC(char *clist, CONFIG **cfg, ARRAY *variables) {
  SHELL *shells;
  int i, j, k, m;
  int n_shells, n;
  char *argv[MAXNARGS], *sv[MAXNARGS];
  int argt[MAXNARGS], st[MAXNARGS], na, ns;
  
  na = 0;
  ns = 0;

  (*cfg) = (CONFIG *) malloc(sizeof(CONFIG));
  n_shells = DecodeArgs(clist, argv, argt, variables);
  na = n_shells;
  (*cfg)->n_shells = n_shells;
  (*cfg)->shells = malloc(sizeof(SHELL)*n_shells);
  shells = (*cfg)->shells;

  for (i = 0, m = n_shells-1; i < n_shells; i++, m--) {
    if (argt[i] != TUPLE) return -1;
    n = DecodeArgs(argv[i], sv, st, variables);
    ns = n;
    if (n != 4) return -1;
    shells[m].n = atoi(sv[0]);
    k = atoi(sv[1]);
    j = atoi(sv[2]);
    if (j > 0) k = -(k+1);
    shells[m].kappa = k;
    shells[m].nq = atoi(sv[3]);
  }

  for (i = 0; i < na; i++) free(argv[i]);
  for (i = 0; i < ns; i++) free(sv[i]);

  return 0;     
}  

static int PAvgConfig(int argc, char *argv[], int argt[], ARRAY *variables) {
  int ns, *n, *kappa;
  double *nq;

  if (argc != 1 || argt[0] != STRING) return -1;

  ns = GetAverageConfigFromString(&n, &kappa, &nq, argv[0]);
  if (ns <= 0) return -1;

  if (SetAverageConfig(ns, n, kappa, nq) < 0) return -1;

  free(n);
  free(kappa);
  free(nq);
  return 0;
}

static int PCheckEndian(int argc, char *argv[], int argt[], ARRAY *variables) {
  TFILE *f;
  F_HEADER fh;
  int i, swp;

  if (argc == 0) {
    i = CheckEndian(NULL);
  } else {
    f = FOPEN(argv[0], "rb");
    if (f == NULL) {
      printf("Cannot open file %s\n", argv[0]);
      return -1;
    }
    ReadFHeader(f, &fh, &swp);
    i = CheckEndian(&fh);
    FCLOSE(f);
  }

  printf("Endian: %d\n", i);
  
  return 0;
}  

static char _closed_shells[MCHSHELL] = "";
static int PClosed(int argc, char *argv[], int argt[], ARRAY *variables) {
  CONFIG *cfg;
  int i, j, kl, n, nq, ncfg;
  char *p;
  char s[16], st[16];
  int ns, k;

  if (argc == 0) _closed_shells[0] = '\0';
  for (i = 0; i < argc; i++) {
    if (argt[i] != STRING) return -1;
    ns = StrSplit(argv[i], ' ');
    p = argv[i];
    for (k = 0; k < ns; k++) {
      while (*p == ' ') p++;
      ncfg = GetConfigFromStringNR(&cfg, p);
      for (j = ncfg-1; j >= 0; j--) {
	if (cfg[j].n_shells != 1) return -1;
	n = (cfg[j].shells)[0].n;
	kl = (cfg[j].shells)[0].kappa;
	nq = 2*(kl + 1);
	kl = kl/2;
	SetClosedShellNR(n, kl);
	SpecSymbol(s, kl);
	sprintf(st, "%d%s%d ", n, s, nq);
	strcat(_closed_shells, st);
	free(cfg[j].shells);
      }
      if (ncfg > 0) free(cfg);
      while (*p) p++;
      p++;
    }
  }
  return 0;
}

static int PGetConfigNR(int argc, char *argv[], int argt[], ARRAY *variables) {
  CONFIG *cfg;
  int i, j, k, t, ncfg;
  char scfg[MCHSHELL], s[16];
  
  for (i = 0; i < argc; i++) {
    if (argt[i] != STRING) return -1;
    strncpy(scfg, _closed_shells, MCHSHELL);
    strncat(scfg, argv[i], MCHSHELL);
    ncfg = GetConfigFromStringNR(&cfg, scfg);
    for (j = 0; j < ncfg; j++) {
      scfg[0] = '\0';
      for (t = cfg[j].n_shells-1; t >= 0; t--) {
	sprintf(s, "%d", (cfg[j].shells)[t].n);
	strcat(scfg, s);
	SpecSymbol(s, (cfg[j].shells)[t].kappa/2);
	strcat(scfg, s);
	if (t == 0) {
	  sprintf(s, "%d", (cfg[j].shells)[t].nq);
	} else {
	  sprintf(s, "%d ", (cfg[j].shells)[t].nq);
	}
	strcat(scfg, s);
      }
      printf("%s\n", scfg);
    }
    if (ncfg > 0) free(cfg);
  }
  
  return 0;
}
  
static int PReadConfig(int argc, char *argv[], int argt[], ARRAY *variables) {
  if (argc < 1 || argc > 2) return -1;
  if (argt[0] != STRING) return -1;
  char *c = NULL;
  if (argc > 1) {
    if (argt[1] != STRING) return -1;
    c = argv[1];
  }
  return ReadConfig(argv[0], c);
}

static int PConfig(int argc, char *argv[], int argt[], ARRAY *variables) {
  CONFIG *cfg;
  static char gname[GROUP_NAME_LEN] = "_all_";
  int i, j, k, t, ncfg;
  char scfg[MCHSHELL], *gn1, *gn2, *s;
  int nf;
  char *v[MAXNARGS];
  int vt[MAXNARGS];

  if (argt[0] == NUMBER) {
    int ng, *kg, ngb, *kgb, n0, n1, k0, k1, m, n0d, n1d;
    double sth;
    
    m = atoi(argv[0]);
    if (m == 0) {
      if (argt[1] != STRING) return -1;
      if (argt[2] != STRING) return -1;
      sth = 0;
      ngb = 0;
      kgb = NULL;
      if (argc > 3) {
	sth = atof(argv[3]);
	if (argc > 4) {
	  if (argt[4] == NUMBER) {
	    ngb = atoi(argv[4]);
	    kgb = NULL;
	  } else {
	    ngb = DecodeGroupArgs(&kgb, 1, NULL, &argv[4], &argt[4], variables);
	  }
	}
      }
      gn1 = argv[1];
      gn2 = NULL;
      s = argv[2];
      n0 = 0;
      n1 = 0;
      n0d = 0;
      n1d = 0;
      k0 = 0;
      k1 = 0;
      t = ConfigSD(m, 0, NULL, s, gn1, gn2, n0, n1, n0d, n1d, k0, k1,
		   ngb, kgb, sth);
      if (ngb > 0 && kgb) free(kgb);
      return t;
    } else {
      if (argt[1] != STRING && argt[1] != LIST) return -1;
      if (argt[2] != LIST && argt[2] != TUPLE) return -1;
      ng = DecodeGroupArgs(&kg, 1, NULL, &argv[2], &argt[2], variables);
      if (argt[1] == STRING) {
	gn1 = argv[1];
	gn2 = NULL;
	nf = 0;
      } else {
	nf = DecodeArgs(argv[1], v, vt, variables);
	if (nf < 1 || nf > 2) return -1;
	if (vt[0] != STRING) return -1;
	gn1 = v[0];
	gn2 = NULL;
	if (nf > 1) {
	  if (vt[1] != STRING) return -1;
	  gn2 = v[1];
	}
      }
      n0 = 1;
      n1 = 1;
      k0 = 0;
      k1 = -1;
      s = NULL;
      if (argc > 3) {
	s = argv[3];
      }
      if (argc > 4) {
	n0 = atoi(argv[4]);
      }
      if (argc > 5) {
	n1 = atoi(argv[5]);
      }
      if (argc > 6) {
	k0 = atoi(argv[6]);
      }
      if (argc > 7) {
	k1 = atoi(argv[7]);
      }
      sth = 0;
      n0d = n0;
      n1d = n1;
      ngb = ng;
      kgb = kg;
      if (argc > 8) {
	n0d = atoi(argv[8]);
	if (argc > 9) {
	  n1d = atoi(argv[9]);
	  if (argc > 10) {
	    sth = atof(argv[10]);
	    if (argc > 11) {
	      if (argt[11] == NUMBER) {
		ngb = 1;
		kgb = NULL;
	      } else {
		ngb = DecodeGroupArgs(&kgb, 1, NULL,
				      &argv[11], &argt[11], variables);
	      }
	    }
	  }
	}
      }
      t = ConfigSD(m, ng, kg, s, gn1, gn2, n0, n1, n0d, n1d, k0, k1,
		   ngb, kgb, sth);
      if (nf > 0) {
	for (i = 0; i < nf; i++) {
	  free(v[i]);
	}
      }
      if (ng > 0) free(kg);
      if (ngb > 0 && kgb && kgb != kg) free(kgb);
      return t;
    }
  }
  k = -2;
  for (i = 0; i < argc; i++) {
    if (argt[i] == KEYWORD) {
      if (strcmp(argv[i], "group") != 0) {
	printf("The keyword must be group=gname\n");
	return -1;
      }
      if (i > argc-2) return -1;
      if (argt[i+1] != STRING) return -1;
      k = i;
    }
  }

  i = 0;
  
  if (k >= 0) strncpy(gname, argv[k+1], GROUP_NAME_LEN);
  else {
    if (argc == 0) return -1;
    if (argt[i] != STRING) return -1;
    strncpy(gname, argv[i], GROUP_NAME_LEN);
    i++;
  }
  
  t = GroupIndex(gname);
  if (t < 0) return -1;

  for (; i < argc; i++) {
    if (i == k || i == k+1) continue;
    if (argt[i] != STRING) return -1;
    strncpy(scfg, _closed_shells, MCHSHELL);
    strncat(scfg, argv[i], MCHSHELL);
    ncfg = GetConfigFromString(&cfg, scfg);
    for (j = 0; j < ncfg; j++) {
      if (Couple(cfg+j) < 0) return -1;
      if (AddConfigToList(t, cfg+j) < 0) return -1;
    }   
    if (ncfg > 0) free(cfg);
  }

  CONFIG_GROUP *g = GetGroup(t);
  if (g != NULL && g->n_cfgs == 0) {
    RemoveGroup(t);
  }
  
  return 0;
}      
  
static int PRemoveConfig(int argc, char *argv[], int argt[], ARRAY *variables) {
  int k, ng, *kg;
  
  if (argc <= 0) return -1;
  ng = DecodeGroupArgs(&kg, argc, NULL, argv, argt, variables);
  
  for (k = 0; k < ng; k++) {
    RemoveGroup(kg[k]);
  }
  ReinitStructure(1);
  ReinitRecouple(0);
  if (ng > 0) free(kg);

  return 0;
}
  
static int PListConfig(int argc, char *argv[], int argt[], ARRAY *variables) {
  int k, ng, *kg;
  char *s;
  
  s = NULL;
  ng = 0;
  if (argc > 0) {
    s = argv[0];
    if (argc > 1) {
      ng = DecodeGroupArgs(&kg, 1, NULL, argv+1, argt+1, variables);
    }
  }
  if (ng <= 0) {
    ng = GetNumGroups();
    kg = malloc(sizeof(int)*ng);
    for (k = 0; k < ng; k++) {
      kg[k] = k;
    }
  }

  ListConfig(s, ng, kg);

  if (ng > 0) free(kg);

  return 0;
}

static int PConfigEnergy(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int m, mr, i;
  int ng, *kg;

  if (argc == 0) return -1;
  m = atoi(argv[0]);
  
  if (argc == 1 || m != 0) {
    ConfigEnergy(m, 0, 0, NULL);
  } else {
    mr = atoi(argv[1]);
    if (argc == 2) {
      ConfigEnergy(m, mr, 0, NULL);
    } else {
      for (i = 1; i < argc; i++) {
	ng = DecodeGroupArgs(&kg, 1, NULL, argv+i, argt+i, variables);
	if (ng < 0) return -1;
	ConfigEnergy(m, mr, ng, kg);
	if (ng > 0) free(kg);
      }
    }
  }

  return 0;
}

static int PAddConfig(int argc, char *argv[], int argt[], ARRAY *variables) {
  CONFIG *cfg = NULL;
  int k;
  
  if (argc != 2) return -1;
  if (argt[0] != STRING) return -1;
  if (argt[1] != LIST) return -1;
  if (ConfigListToC(argv[1], &cfg, variables) < 0) return -1;
  if (Couple(cfg) < 0) return -1;

  k = GroupIndex(argv[0]);
  if (k < 0) return -1;
  if (AddConfigToList(k, cfg) < 0) return -1;
  free(cfg);
  
  return 0;
}

static int PAITable(int argc, char *argv[], int argt[], ARRAY *variables) {
  int nlow, *low, nup, *up;
  double c;

  if (argc != 3 && argc != 4) return -1;
  if (argt[0] != STRING) return -1;
  
  if (argc == 4) {
    if (argt[3] != NUMBER) return -1;
    c = atof(argv[3]);
  } else c = 0.0;
  
  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  SaveAI(nlow, low, nup, up, argv[0], c, 0);
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  return 0;
}

static int PAITableMSub(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int nlow, *low, nup, *up;
  double c;

  if (argc != 3 && argc != 4) return -1;
  if (argt[0] != STRING) return -1;
  
  if (argc == 4) {
    if (argt[3] != NUMBER) return -1;
    c = atof(argv[3]);
  } else c = 0.0;
  
  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  SaveAI(nlow, low, nup, up, argv[0], c, 1);
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  return 0;
}

static int PBasisTable(int argc, char *argv[], int argt[], ARRAY *variables) {
  int m, k;

  if (argc == 0) return -1;
  if (argc > 3 || argt[0] != STRING) return -1;
  m = 0;
  k = -1;
  if (argc > 1) {
    m = atoi(argv[1]);
    if (argc > 2) {
      k = atoi(argv[2]);
    }
  }

  GetBasisTable(argv[0], m, k);
  
  return 0;
}

static int PCETable(int argc, char *argv[], int argt[], ARRAY *variables) {
  int nlow, nup, *low, *up;

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;
  
  if (argc == 1) {
    if (argt[0] != STRING ) return -1;
    SaveExcitation(nlow, low, nup, up, 0, argv[0]);
  } else if (argc == 2) {
    if (argt[0] != STRING) return -1;
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    if (nlow <= 0) return -1;
    SaveExcitation(nlow, low, nlow, low, 0, argv[0]);
    free(low);
  } else if (argc == 3) {
    if (argt[0] != STRING) return -1;
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    if (nlow <= 0) return -1;
    nup = SelectLevels(&up, argv[2], argt[2], variables);
    if (nup <= 0) return -1;
    SaveExcitation(nlow, low, nup, up, 0, argv[0]);
    free(low);
    free(up);
  } else {
    return -1;
  }
    
  return 0;
}

static int PCETableMSub(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int nlow, nup, *low, *up;

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;
  
  if (argc == 1) {
    if (argt[0] != STRING ) return -1;
    SaveExcitation(nlow, low, nup, up, 1, argv[0]);
  } else if (argc == 2) {
    if (argt[0] != STRING) return -1;
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    if (nlow <= 0) return -1;
    SaveExcitation(nlow, low, nlow, low, 1, argv[0]);
    free(low);
  } else if (argc == 3) {
    if (argt[0] != STRING) return -1;
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    if (nlow <= 0) return -1;
    nup = SelectLevels(&up, argv[2], argt[2], variables);
    if (nup <= 0) return -1;
    SaveExcitation(nlow, low, nup, up, 1, argv[0]);
    free(low);
    free(up);
  } else {
    return -1;
  }
    
  return 0;
}

static int PCITable(int argc, char *argv[], int argt[], ARRAY *variables) {
  int nlow, *low, nup, *up;
  
  if (argc != 3) return -1;
  if (argt[0] != STRING) return -1;
  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  if (nlow <= 0) return -1;
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  if (nup <= 0) return -1;
  if (SaveIonization(nlow, low, nup, up, argv[0]) < 0) return -1;

  if (nlow > 0) free(low);
  if (nup > 0) free(up);
    
  return 0;
}

static int PCITableMSub(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int nlow, *low, nup, *up;
  
  if (argc != 3) return -1;
  if (argt[0] != STRING) return -1;
  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  if (nlow <= 0) return -1;
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  if (nup <= 0) return -1;
  if (SaveIonizationMSub(nlow, low, nup, up, argv[0]) < 0) return -1;

  if (nlow > 0) free(low);
  if (nup > 0) free(up);
    
  return 0;
}

static int PClearLevelTable(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  if (argc != 0) return -1;
  ClearLevelTable();
 
  return 0;
}

static int PClearOrbitalTable(int argc, char *argv[], int argt[], 
			      ARRAY *variables) {
  int m;
  
  m = 1;
  if (argc != 0 && argc != 1) return -1;
  if (argc == 1) {
    if (argt[0] != NUMBER) return -1;
    m = atoi(argv[0]);
  }

  ClearOrbitalTable(m);
  return 0;
}

static int PAdjustEnergy(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int nlevs, k, *ilevs;
  double e, *elevs;
  int i, ie, ii;
  FILE *f;
  char *iv[MAXNARGS], *ev[MAXNARGS];
  int it[MAXNARGS], et[MAXNARGS];

  ii = 0; 
  ie = 0;
  if (argc == 5) {
    if (argt[0] != STRING) {
      return -1;
    }
    f = fopen(argv[0], "r");
    
    i = 0;
    while (1) {
      if (fscanf(f, "%d%lf\n", &k, &e) == EOF) break;
      i++;
    }
    nlevs = i;    
    ilevs = (int *) malloc(sizeof(int)*nlevs);
    elevs = (double *) malloc(sizeof(double)*nlevs);
    fseek(f, 0, SEEK_SET);
    i = 0;
    while (1) {
      if (fscanf(f, "%d%lf\n", &k, &e) == EOF) break;
      e /= HARTREE_EV;
      ilevs[i] = k;
      elevs[i] = e;
      i++;
    }
    fclose(f);
    AdjustEnergy(nlevs, ilevs, elevs, argv[1], argv[2], argv[3], argv[4]);
  } else {
    if (argt[0] != LIST || argt[1] != LIST) {
      printf("The last two of three arguments ");
      printf("for CorrectEnergy must be two Lists\n");
      return -1;
    }
    ii = DecodeArgs(argv[0], iv, it, variables);
    ie = DecodeArgs(argv[1], ev, et, variables);
    if (ii != ie) return -1;
    nlevs = ii;
    ilevs = (int *) malloc(sizeof(int)*nlevs);
    elevs = (double *) malloc(sizeof(double)*nlevs);
    for (i = 0; i < nlevs; i++) {
      if (it[i] != NUMBER || et[i] != NUMBER) return -1;
      k = atoi(iv[i]);
      e = atof(ev[i]);
      e /= HARTREE_EV;
      ilevs[i] = k;
      elevs[i] = e;
    }
    AdjustEnergy(nlevs, ilevs, elevs, argv[2], argv[3], argv[4], argv[5]);
  }

  for (i = 0; i < ii; i++) free(iv[i]);
  for (i = 0; i < ie; i++) free(ev[i]);
  if (nlevs > 0) {
    free(ilevs);
    free(elevs);
  }
  
  return 0;
}

static int PCorrectEnergy(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  int n, k, kref;
  double e;
  int i, ie, ii, nmin;
  FILE *f;
  char *iv[MAXNARGS], *ev[MAXNARGS];
  int it[MAXNARGS], et[MAXNARGS];

  ii = 0; 
  ie = 0;
  kref = 0;
  if (argc == 2) {
    if (argt[0] != STRING) {
      return -1;
    }
    nmin = atoi(argv[1]);
    f = fopen(argv[0], "r");
    n = -1;
    i = 0;
    while (1) {
      if (fscanf(f, "%d%lf\n", &k, &e) == EOF) break;
      e /= HARTREE_EV;
      if (k < 0) {
	k = -k;
	kref = k;
      } else if (i == 0) {
	kref = k;
      }
      AddECorrection(kref, k, e, nmin);
      i++;
    }
    fclose(f);
  } else if (argc == 3) {
    if (argt[0] != LIST || argt[1] != LIST) {
      printf("The last two of three arguments ");
      printf("for CorrectEnergy must be two Lists\n");
      return -1;
    }
    nmin = atoi(argv[2]);
    ii = DecodeArgs(argv[0], iv, it, variables);
    ie = DecodeArgs(argv[1], ev, et, variables);
    if (ii != ie) return -1;
    n = ii;
    for (i = 0; i < n; i++) {
      if (it[i] != NUMBER || et[i] != NUMBER) return -1;
      k = atoi(iv[i]);
      e = atof(ev[i]);
      e /= HARTREE_EV;
      if (k < 0) {
	k = -k;
	kref = k;
      } else if (i == 0) {
	kref = k;
      }
      AddECorrection(kref, k, e, nmin);
    }
  } else {
    return -1;
  }

  for (i = 0; i < ii; i++) free(iv[i]);
  for (i = 0; i < ie; i++) free(ev[i]);

  return 0;
}

static int PExit(int argc, char *argv[], int argt[], ARRAY *variables) {
  if (argc != 0) return -1;
  exit(0);
}

static int PFreeExcitationQk(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  if (argc != 0) return -1;
  FreeExcitationQk();
  return 0;
}

static int PFreeIonizationQk(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  if (argc != 0) return -1;
  FreeIonizationQk();
  return 0;
}

static int PFreeMemENTable(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  
  if (argc != 0) return -1;
  FreeMemENTable();
  return 0;
}

static int PFreeMultipole(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  if (argc != 0) return -1;
  FreeMultipoleArray();
  FreeMomentsArray();
  return 0;
}

static int PFreeSlater(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  if (argc != 0) return -1;
  FreeSlaterArray();
  return 0;
}

static int PFreeResidual(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  if (argc != 0) return -1;
  FreeResidualArray();
  return 0;
}

static int PFreeRecPk(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  if (argc != 0) return -1;
  FreeRecPk();
  return 0;
}

static int PFreeRecQk(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  if (argc != 0) return -1;
  FreeRecQk();
  return 0;
}

static int PGetPotential(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  if (argc != 1 || argt[0] != STRING) return -1;
  GetPotential(argv[0]);
  return 0;
}

static int PInfo(int argc, char *argv[], int argt[], ARRAY *variables) {
  if (argc != 0) return -1;
  Info();
  return 0;
}

static int PMemENTable(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {

  if (argc != 1) return -1;
  if (argt[0] != STRING) return -1;

  MemENTable(argv[0]);

  return 0;
}

static int POptimizeRadial(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  int ng, i, k;
  int *kg;
  double z;
  double *weight;
  char *vw[MAXNARGS];
  int iw[MAXNARGS], ni;
  
  ni = 0;

  ng = argc;
  if (ng == 0) {
    ng = 0; 
    kg = NULL;
    weight = NULL;
    goto END;
  } 
  
  if (argt[0] == STRING) {
    weight = NULL;
    ng = DecodeGroupArgs(&kg, argc, NULL, argv, argt, variables);
    if (ng < 0) return -1;
  } else {
    ng = DecodeGroupArgs(&kg, 1, NULL, argv, argt, variables);
    if (ng < 0) return -1;
  
    if (argc == 1) {
      weight = NULL;
    } else {
      if (argt[1] != LIST && argt[1] != TUPLE) return -1;
      k = DecodeArgs(argv[1], vw, iw, variables);
      ni = k;
      if (k < 0 || k > ng) {
	printf("weights must be a sequence\n");
	return -1;
      } 
      weight = (double *) malloc(sizeof(double)*ng);
      z = 0.0;
      for (i = 0; i < k; i++) {
	if (iw[i] != NUMBER) {
	  return -1;
	} 
	weight[i] = atof(vw[i]);
	z += weight[i];
      }
      for (i = k; i < ng; i++) {
	if (z >= 1.0) {
	  weight[i] = weight[k-1];
	} else {
	  weight[i] = (1.0-z)/(ng-k);
	}
      }
    }
  }

 END:
  if (ng == 0) kg = NULL;  
  if (OptimizeRadial(ng, kg, -1, weight, 0) < 0) {
    if (kg) free(kg);
    if (weight) free(weight);
    return -1;
  }
  if (weight) free(weight);
  if (kg) free(kg);

  for (i = 0; i < ni; i++) free(vw[i]);

  return 0;
}

static int PRefineRadial(int argc, char *argv[], int argt[], 
		  ARRAY *variables) {
  int maxfun, msglvl;
  
  maxfun = 0;
  msglvl = 0;
  if (argc > 0) {
    maxfun = atoi(argv[0]);
    if (argc > 1) {
      msglvl = atoi(argv[1]);
    }
  }
  
  return RefineRadial(maxfun, msglvl);
}

static int PPause(int argc, char *argv[], int argt[], 
		  ARRAY *variables) {
  char s[10];

  while (1) {
    printf("Type go to continue: ");
    scanf("%s", s);
    if (strcmp(s, "go") == 0) break;
  }
  
  return 0;
}

static int PPrintMemInfo(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
#ifdef DEBUG_ARRAY
  PrintMemInfo();
#endif
 
  return 0;
}

static int PPrintTable(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int v;

  if (argc != 2 && argc != 3) return -1;
  if (argt[0] != STRING || argt[1] != STRING) return -1;
  
  v = 1;
  if (argc == 3) {
    if (argt[2] != NUMBER) return -1;
    v = atoi(argv[2]);
  }

  PrintTable(argv[0], argv[1], v);
  return 0;
}

static int PRecStates(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int ng, *kg, n;

  if (argc != 3 || argt[0] != STRING || 
      (argt[1] != LIST && argt[1] != TUPLE) ||
      argt[2] != NUMBER) 
    return -1;
  ng = DecodeGroupArgs(&kg, 1, NULL, &(argv[1]), &(argt[1]), variables);
  if (ng <= 0) return -1;
  n = atoi(argv[2]);
  if (RecStates(n, ng, kg, argv[0]) < 0) {
    printf("RecStates Error\n");
    free(kg);
    return -1;
  }

  free(kg);
  
  return 0;
}

static int PReinitConfig(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  m = atoi(argv[0]);

  ReinitConfig(m);
  _closed_shells[0] = '\0';

  return 0;
}

static int PReinitRecouple(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  m = atoi(argv[0]);

  ReinitRecouple(m);
  return 0;
}

static int PReinitRadial(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  m = atoi(argv[0]);

  ReinitRadial(m);
  return 0;
}

static int PReinitDBase(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  m = atoi(argv[0]);

  ReinitDBase(m);
  return 0;
}

static int PReinitStructure(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  m = atoi(argv[0]);

  ReinitStructure(m);
  return 0;
}

static int PReinitExcitation(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  m = atoi(argv[0]);

  ReinitExcitation(m);
  return 0;
}

static int PReinitRecombination(int argc, char *argv[], int argt[], 
				ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  m = atoi(argv[0]);

  ReinitRecombination(m);
  return 0;
}

static int PReinitIonization(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  m = atoi(argv[0]);

  ReinitIonization(m);
  return 0;
}

static int PReinit(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int i, m;  
  int m_config;
  int m_recouple;
  int m_radial;
  int m_dbase;
  int m_structure;
  int m_excitation;
  int m_recombination;
  int m_ionization;

  m_config = -1;
  m_recouple = -1;
  m_radial = -1;
  m_dbase = -1;
  m_structure = -1;
  m_excitation = -1;
  m_recombination = -1;
  m_ionization = -1;

  if (argc == 0 || argc == 1) {
    if (argc == 0) {
      m = 0;
    } else {
      if (argt[0] != NUMBER) return -1;
      m = atoi(argv[0]);
    } 
    if (m == 0) {
      m_config = 0;
      m_recouple = 0;
      m_radial = 0;
      m_dbase = 0;
      m_structure = 0;
      m_excitation = 0;
      m_recombination = 0;
      m_ionization = 0;
    } else if (m > 0) {
      m_config = -1;
      m_recouple = -1;
      m_radial = 0;
      m_dbase = 0;
      m_structure = 2;
      m_excitation = 0;
      m_recombination = 0;
      m_ionization = 0;
    } else {
      m_config = 1;
      m_recouple = 1;
      m_radial = 1;
      m_dbase = 1;
      m_structure = 1;
      m_excitation = 1;
      m_recombination = 1;
      m_ionization = 1;
    }
  } else {
    for (i = 0; i < argc; ) {
      if (argt[i] != KEYWORD) {
	i++;
	continue;
      }
      if (strcmp(argv[i], "config") == 0) {
	i++;
	if (argt[i] != NUMBER) return -1;
	m_config = atoi(argv[i]);
	i++;
      } else if (strcmp(argv[i], "recouple") == 0) {
	i++;
	if (argt[i] != NUMBER) return -1;
	m_recouple = atoi(argv[i]);
	i++;
      } else if (strcmp(argv[i], "dbase") == 0) {
	i++;
	if (argt[i] != NUMBER) return -1;
	m_dbase = atoi(argv[i]);
	i++;
      } else if (strcmp(argv[i], "structure") == 0) {
	i++;
	if (argt[i] != NUMBER) return -1;
	m_structure = atoi(argv[i]);
	i++;
      } else if (strcmp(argv[i], "excitation") == 0) {
	i++;
	if (argt[i] != NUMBER) return -1;
	m_excitation = atoi(argv[i]);
	i++;
      } else if (strcmp(argv[i], "radial") == 0) {
	i++;
	if (argt[i] != NUMBER) return -1;
	m_radial = atoi(argv[i]);
	i++;
      } else if (strcmp(argv[i], "recombination") == 0) {
	i++;
	if (argt[i] != NUMBER) return -1;
	m_recombination = atoi(argv[i]);
	i++;
      } else if (strcmp(argv[i], "ionization") == 0) {
	i++;
	if (argt[i] != NUMBER) return -1;
	m_ionization = atoi(argv[i]);
	i++;
      }
    }
  }
 
  ReinitFac(m_config, m_recouple, m_radial, m_dbase,
	    m_structure, m_excitation, m_recombination, m_ionization);
  if (m_config == 0) {
    _closed_shells[0] = '\0';
  }

  return 0;
}
  
static int PRRMultipole(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int nlow, *low, nup, *up, m;
  
  m = -1;
  if (argc != 3 && argc != 4) return -1;
  if (argt[0] != STRING) return -1;

  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  if (nlow <= 0) return -1;
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  if (nup <= 0) return -1;
  if (argc == 4) {
    if (argt[3] != NUMBER) return -1;
    m = atoi(argv[3]);
  }
  SaveRRMultipole(nlow, low, nup, up, argv[0], m);

  free(low);
  free(up);

  return 0;
}
  
static int PRRTable(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  int nlow, *low, nup, *up, m;
  
  m = -1;
  if (argc != 3 && argc != 4) return -1;
  if (argt[0] != STRING) return -1;

  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  if (nlow <= 0) return -1;
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  if (nup <= 0) return -1;
  if (argc == 4) {
    if (argt[3] != NUMBER) return -1;
    m = atoi(argv[3]);
  }
  SaveRecRR(nlow, low, nup, up, argv[0], m);

  free(low);
  free(up);

  return 0;
}

static int PSetAICut(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  double c;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  
  c = atof(argv[0]);
  SetAICut(c);

  return 0;
}

static int PSetAngZOptions(int argc, char *argv[], int argt[],
			   ARRAY *variables) {
  int n;
  double c, mc;

  n = atoi(argv[0]);
  c = EPS3;
  mc = EPS3;
  if (argc > 1) {
    mc = atof(argv[1]);
    if (argc > 2) {
      c = atof(argv[2]);
    }
  }
  SetAngZOptions(n, mc, c);
  
  return 0;
}

static int PSetCILevel(int argc, char *argv[], int argt[],
		       ARRAY *variables) {
  int i;
  
  if (argc != 1 || argt[0] != NUMBER) return -1;
  i = atoi(argv[0]);
  SetCILevel(i);

  return 0;
}

static int PSetAngZCut(int argc, char *argv[], int argt[],
		       ARRAY *variables) {
  double c;
  
  if (argc != 1 || argt[0] != NUMBER) return -1;
  c = atof(argv[0]);
  SetAngZCut(c);

  return 0;
}

static int PSetMixCut(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  double c, c2;
  
  if (argc < 1 || argc > 2) return -1;
  if (argt[0] != NUMBER) return -1;
  c = atof(argv[0]);
  c2 = -1.0;
  if (argc > 1) {
    if (argt[1] != NUMBER) return -1;
    c2 = atof(argv[1]);
  }
  SetMixCut(c, c2);
  
  return 0;
}
static int PSetExtraPotential(int argc, char *argv[], int argt[], 
			      ARRAY *variables) {
  int m;
  int n;
  double *p;

  if (argc < 1 || argc > 2) return -1;
  if (argt[0] != NUMBER) return -1;
  m = atoi(argv[0]);  
  n = 0;
  p = NULL;
  if (m >= 0 && argc == 2) {
    n = DoubleFromList(argv[1], argt[1], variables,  &p);
  }
  SetExtraPotential(m, n, p);
  return 0;
}
static int PSetAtom(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  double z, mass, rn, a, npr;

  mass = -1.0;
  z = -1.0;
  rn = -1.0;
  a = -1.0;
  npr = -1.0;
  if (argc < 1) return -1;
  if (argt[0] == STRING) {
    if (argc > 6) return -1;
    if (argc > 1) {
      z = atof(argv[1]);
      if (argc > 2) {
	mass = atof(argv[2]);
	if (argc > 3) {
	  rn = atof(argv[3]);
	  if (argc > 4) {
	    a = atof(argv[4]);
	    if (argc > 5) {
	      npr = atof(argv[5]);
	    }
	  }
	}
      }
    }
    if (SetAtom(argv[0], z, mass, rn, a, npr) < 0) return -1;
  } else if (argt[0] == NUMBER) {
    if (argc > 5) return -1;
    z = atof(argv[0]);
    if (argc > 1) {
      mass = atof(argv[1]);
      if (argc > 2) {
	rn = atof(argv[2]);
	if (argc > 3) {
	  a = atof(argv[3]);
	  if (argc > 4) {
	    npr = atof(argv[4]);
	  }
	}
      }
    }
    if (SetAtom(NULL, z, mass, rn, a, npr) < 0) return -1;
  }
  
  return 0;
}

static int PSetAvgConfig(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int ns, i, j;
  int *n, *kappa;
  double *nq;
  char *vc[MAXNARGS], *vs[MAXNARGS];
  int ic[MAXNARGS], is[MAXNARGS];

  if (argc != 1 || argt[0] != LIST) {
    return -1;
  }
  
  ns = DecodeArgs(argv[0], vc, ic, variables);
  
  n = malloc(sizeof(int)*ns);
  kappa = malloc(sizeof(int)*ns);
  nq = malloc(sizeof(double)*ns);
  
  for (i = 0; i < ns; i++) {
    if (DecodeArgs(vc[i], vs, is, variables) != 4) {
      return -1;
    }
    n[i] = atoi(vs[0]);
    kappa[i] = atoi(vs[1]);
    if (atoi(vs[2]) > 0) kappa[i] = -(kappa[i]+1);
    nq[i] = atof(vs[3]);
    for (j = 0; j < 4; j++) free(vs[j]);
  }
  
  for (j = 0; j < ns; j++) free(vc[j]);

  if (SetAverageConfig(ns, n, kappa, nq) < 0) return -1;
  free(n);
  free(kappa);
  free(nq);

  return 0;
}

static int PSetCEGrid(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int n, ng, i, err;
  double xg[MAXNE];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0; 
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetCEEGrid(ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
	xg[i] = atof(vg[i]);
	free(vg[i]);
	xg[i] /= HARTREE_EV;
      }
      err = SetCEEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetCEEGrid(ng, emin, emax, eth);
  } else {
    return -1;
  }
  
  if (err < 0) return -1;
  return 0;
}

static int PSetAngleGrid(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int n, ng, i, err, m;
  double xg[MAXNTHETA+MAXNPHI];
  double emin, emax;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  n = argc;

  if (n == 2) {
    m = atoi(argv[0]);
    if (argt[1] == NUMBER) {
      ng = atoi(argv[1]);
      if (m == 0) {
	emin = 0.0;
	emax = PI;
      } else {
	emin = 0.0;
	emax = TWO_PI;
      }
      err = SetAngleGrid(m, ng, emin, emax);
    } else if (argt[1] == LIST || argt[1] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
	xg[i] = atof(vg[i]);
	free(vg[i]);
	xg[i] /= HARTREE_EV;
      }
      err = SetAngleGridDetail(m, ng, xg);
    } else {
      return -1;
    }
  } else if (n == 4) {
    m = atoi(argv[0]);
    ng = atoi(argv[1]);
    emin = atof(argv[2]);
    emax = atof(argv[3]);
    emin *= PI/180.0;
    emax *= PI/180.0;
    err = SetAngleGrid(m, ng, emin, emax);
  } else {
    return -1;
  }
  
  if (err < 0) return -1;
  return 0;
}

static int PSetCEGridLimits(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  double emin, emax;
  int type;

  emin = -1;
  emax = -1;
  type = 0;
  
  if (argc > 0) {
    emin = atof(argv[0]);
    if (argc > 1) {
      emax = atof(argv[1]);
      if (argc > 2) {
	type = atoi(argv[2]);
      }
    }
  }
  
  SetCEEGridLimits(emin, emax, type);

  return 0;
}

static int PSetCEGridType(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetCEEGridType(atoi(argv[0]));
  return 0;
}

static int PSetTEGrid(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int i, n;
  double xg[MAXNTE];
  int ng, err;
  double emin, emax;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  n = argc;
  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetCETEGrid(ng, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
	if (ig[i] != NUMBER) return -1;
	xg[i] = atof(vg[i]);
	xg[i] /= HARTREE_EV;
	free(vg[i]);
      }
      err = SetCETEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    err = SetCETEGrid(ng, emin, emax);
  } else {
    return -1;
  }

  if (err < 0) return -1;

  return 0;
}

static int PSetCEPWOptions(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {

  int qr, max, kl_cb;
  double tol;

  qr = EXCLQR;
  max = EXCLMAX;
  kl_cb = EXCLCB;
  tol = EXCTOL;
  
  if (argc < 1 || argc > 4) return -1;
  
  tol = atof(argv[0]);
  if (argc > 1) {
    max = atoi(argv[1]);
    if (argc > 2) {
      qr = atoi(argv[2]);
      if (argc > 3) {
	kl_cb = atoi(argv[3]);
      }
    }
  }

  SetCEPWOptions(qr, max, kl_cb, tol);

  return 0;
}

static int PSetCEPWGridType(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetCEPWGridType(atoi(argv[0]));

  return 0;
}

static int PSetCEPWGrid(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int ns, i;
  int n;
  int *m, *step;
  char *v1[MAXNARGS], *v2[MAXNARGS];
  int t1[MAXNARGS], t2[MAXNARGS];

  n = argc;
  if (n == 1) {
    if (argt[0] != NUMBER) return -1;
    ns = atoi(argv[0]);
    SetCEPWGrid(-ns, NULL, NULL);
  } else {
    if (argt[0] != LIST || argt[1] != LIST) return -1;
    ns = DecodeArgs(argv[0], v1, t1, variables);
    if (ns <= 0) return -1;
    if (ns != DecodeArgs(argv[1], v2, t2, variables)) return -1;
    m = (int *) malloc(ns*sizeof(int));
    step = (int *) malloc(ns*sizeof(int));
    for (i = 0; i < ns; i++) {
      if (t1[i] != NUMBER || t2[i] != NUMBER) return -1;
      m[i] = atoi(v1[i]);
      step[i] = atoi(v2[i]);
      free(v1[i]);
      free(v2[i]);
    }
    SetCEPWGrid(ns, m, step);
    free(m);
    free(step);
  }

  return 0;
}

static int PSetCEBorn(int argc, char *argv[], int argt[],
		      ARRAY *variables) {
  double eb, x, x1, x0;

  if (argc < 1 || argc > 4) return -1;
  if (argt[0] != NUMBER) return -1;
  x0 = XBORN0;
  x1 = XBORN1;
  x = XBORN;
  if (argc > 1) {
    if (argt[1] != NUMBER) return -1;
    x = atof(argv[1]);
    if (argc > 2) {
      if (argt[2] != NUMBER) return -1;
      x1 = atof(argv[2]);
      if (argc > 3) {
        if (argt[3] != NUMBER) return -1;
        x0 = atof(argv[3]);
      }
    }
  }

  eb = atof(argv[0]);
  SetCEBorn(eb, x, x1, x0);
  
  return 0;
}

static int PSetCIBorn(int argc, char *argv[], int argt[],
		      ARRAY *variables) {
  int x;

  if (argc != 1) return -1;
  if (argt[0] != NUMBER) return -1;

  x = atoi(argv[0]);
  SetCIBorn(x);
  
  return 0;
}

static int PSetBornFormFactor(int argc, char *argv[], int argt[],
			      ARRAY *variables) {
  double te;
  char *fn;

  if (argc < 1 || argc > 2) return -1;
  if (argt[0] != NUMBER) return -1;
  te = atof(argv[0]);
  if (argc == 2) fn = argv[1];
  else fn = NULL;
  
  SetBornFormFactor(te, fn);
  
  return 0;
}

static int PSetBornMass(int argc, char *argv[], int argt[],
			ARRAY *variables) {
  double m;

  if (argc != 1) return -1;
  if (argt[0] != NUMBER) return -1;
  m = atof(argv[0]);

  SetBornMass(m);
  
  return 0;
}

static int PSetCEQkMode(int argc, char *argv[], int argt[],
			ARRAY *variables) {
  int m;
  double tol;
  
  m = QK_DEFAULT;
  tol = -1;
  
  if (argc > 2) return -1;
  if (argc > 0) {
    if (argt[0] == STRING) {
      if (strcasecmp(argv[0], "exact") == 0) m = QK_EXACT;
      else if (strcasecmp(argv[0], "interpolate") == 0) m = QK_INTERPOLATE;
      else if (strcasecmp(argv[0], "fit") == 0) m = QK_FIT;
      else return -1;
    } else if (argt[0] == NUMBER) {
      m = atoi(argv[0]);
      if (m >= QK_CB) return -1;
    }
    if (argc > 1) {
      if (argt[1] != NUMBER) return -1;
      tol = atof(argv[1]);
    }
  }

  SetCEQkMode(m, tol);
  
  return 0;
}
    
static int PSetCIEGrid(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {  int n, ng, i, err;
  double xg[MAXNE];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0; 
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetCIEGrid(ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
	xg[i] = atof(vg[i]);
	xg[i] /= HARTREE_EV;
	free(vg[i]);
      }
      err = SetCIEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetCIEGrid(ng, emin, emax, eth);
  } else {
    return -1;
  }
  
  if (err < 0) return -1;
  return 0;
}

static int PSetCIEGridLimits(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  double emin, emax;
  int type;

  emin = -1;
  emax = -1;
  type = 0;
  
  if (argc > 0) {
    emin = atof(argv[0]);
    if (argc > 1) {
      emax = atof(argv[1]);
      if (argc > 2) {
	type = atoi(argv[2]);
      }
    }
  }
  
  SetCIEGridLimits(emin, emax, type);

  return 0;
}

static int PSetIEGrid(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {  int i, n;
  double xg[MAXNTE];
  int ng, err;
  double emin, emax;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  n = argc;
  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetIEGrid(ng, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
	if (ig[i] != NUMBER) return -1;
	xg[i] = atof(vg[i]);
	xg[i] /= HARTREE_EV;
	free(vg[i]);
      }
      err = SetIEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    err = SetIEGrid(ng, emin, emax);
  } else {
    return -1;
  }

  if (err < 0) return -1;

  return 0;
}

static int PSetCIPWOptions(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  int qr, max, max_1, kl_cb;
  double tol;

  qr = IONLQR;
  max = IONLMAX;
  max_1 = IONLEJEC;
  kl_cb = IONLCB;
  tol = IONTOL;
  
  if (argc < 1 || argc > 5) return -1;
  tol = atof(argv[0]);
  if (argc > 1) {
    max = atoi(argv[1]);
    if (argc > 2) {
      max_1 = atoi(argv[2]);
      if (argc > 3) {
	qr = atoi(argv[3]);
	if (argc > 4) {
	  kl_cb = atoi(argv[4]);
	}
      }
    }
  }

  SetCIPWOptions(qr, max, max_1, kl_cb, tol);

  return 0;
}

static int PSetCIPWGrid(int argc, char *argv[], int argt[], 
			ARRAY *variables) {  int ns, i;
  int n;
  int *m, *step;
  char *v1[MAXNARGS], *v2[MAXNARGS];
  int t1[MAXNARGS], t2[MAXNARGS];

  n = argc;
  if (n == 1) {
    if (argt[0] != NUMBER) return -1;
    ns = atoi(argv[0]);
    SetCIPWGrid(-ns, NULL, NULL);
  } else {
    if (argt[0] != LIST || argt[1] != LIST) return -1;
    ns = DecodeArgs(argv[0], v1, t1, variables);
    if (ns <= 0) return -1;
    if (ns != DecodeArgs(argv[1], v2, t2, variables)) return -1;
    m = (int *) malloc(ns*sizeof(int));
    step = (int *) malloc(ns*sizeof(int));
    for (i = 0; i < ns; i++) {
      if (t1[i] != NUMBER || t2[i] != NUMBER) return -1;
      m[i] = atoi(v1[i]);
      step[i] = atoi(v2[i]);
      free(v1[i]);
      free(v2[i]);
    }
    SetCIPWGrid(ns, m, step);
    free(m);
    free(step);
  }

  return 0;
}

static int PSetCIQkMode(int argc, char *argv[], int argt[],
			ARRAY *variables) {
  int m;
  double tol;
  
  m = QK_DEFAULT;
  tol = -1;
  
  if (argc > 2) return -1;
  if (argc > 0) {
    if (argt[0] == STRING) {
      if (strcasecmp(argv[0], "cb") == 0) m = QK_CB;
      else if (strcasecmp(argv[0], "bed") == 0) m = QK_BED;
      else if (strcasecmp(argv[0],"dw") == 0) m = QK_DW;
      else return -1;
    } else if (argt[0] == NUMBER) {
      m = atoi(argv[0]);
      if (m < QK_CB) return -1;
    }
    if (argc > 1) {
      if (argt[1] != NUMBER) return -1;
      tol = atof(argv[1]);
    }
  }

  SetCIQkMode(m, tol);
  
  return 0;
}

static int PSetMaxRank(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }
  SetMaxRank(2*atoi(argv[0]));

  return 0;
}

static int PSetOptimizeControl(int argc, char *argv[], int argt[], 
			       ARRAY *variables) {
  int maxiter, iprint;
  double tol, s;

  iprint = 0;
  
  if (argc != 3 && argc != 4) return -1;
  tol = atof(argv[0]);
  s = atof(argv[1]);
  maxiter = atoi(argv[2]);
  if (argc == 4) iprint = atoi(argv[3]);
  
  SetOptimizeControl(tol, s, maxiter, iprint);

  return 0;
}

static int PSetHydrogenicNL(int argc, char *argv[], int argt[], 
			    ARRAY *variables){
  int n, k, nm, km;
  
  n = -1;
  k = -1;
  nm = -1;
  km = -1;

  if (argc > 0) {
    n = atoi(argv[0]);
    if (argc > 1) {
      k = atoi(argv[1]);
      if (argc > 2) {
	nm = atoi(argv[2]);
	if (argc > 3) {
	  km = atoi(argv[3]);
	}
      }
    }
  }
  SetHydrogenicNL(n, k, nm, km);

  return 0;
}  

static int PSetPEGrid(int argc, char *argv[], int argt[], 
		      ARRAY *variables) { int n, ng, i, err;
  double xg[MAXNE];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0; 
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetPEGrid(ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
	xg[i] = atof(vg[i]);
	xg[i] /= HARTREE_EV;
	free(vg[i]);
      }
      err = SetPEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetPEGrid(ng, emin, emax, eth);
  } else {
    return -1;
  }
  
  if (err < 0) return -1;

  return 0;
}

static int PSetPEGridLimits(int argc, char *argv[], int argt[], 
			    ARRAY *variables) { 
  double emin, emax;
  int type;

  emin = -1;
  emax = -1;
  type = 0;
  
  if (argc > 0) {
    emin = atof(argv[0]);
    if (argc > 1) {
      emax = atof(argv[1]);
      if (argc > 2) {
	type = atoi(argv[2]);
      }
    }
  }
  
  SetPEGridLimits(emin, emax, type);
  
  return 0;
}

static int PSetRecPWOptions(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  int kl_interp, max_kl;
  
  max_kl = -1;
  if (argc == 1) {
    if (argt[0] != NUMBER) return -1;
    kl_interp = atoi(argv[0]);
  } else if (argc == 2) {
    if (argt[0] != NUMBER) return -1;
    kl_interp = atoi(argv[0]);
    if (argt[1] != NUMBER) return -1;
    max_kl = atoi(argv[1]);
  } else {
    return -1;
  }
  
  if (max_kl < 0) max_kl = kl_interp;
  
  SetRecPWOptions(kl_interp, max_kl);
    
  return 0;
}

static int PSetRecQkMode(int argc, char *argv[], int argt[],
			ARRAY *variables) {
  int m;
  double tol;
  
  m = QK_DEFAULT;
  tol = -1;
  
  if (argc > 2) return -1;
  if (argc > 0) {
    if (argt[0] == STRING) {
      if (strcasecmp(argv[0], "exact") == 0) m = QK_EXACT;
      else if (strcasecmp(argv[0], "interpolate") == 0) m = QK_INTERPOLATE;
      else if (strcasecmp(argv[0], "fit") == 0) m = QK_FIT;
      else return -1;
    } else if (argt[0] == NUMBER) {
      m = atoi(argv[0]);
      if (m >= QK_CB) return -1;
    }
    if (argc > 1) {
      if (argt[1] != NUMBER) return -1;
      tol = atof(argv[1]);
    }
  }

  SetRecQkMode(m, tol);
  
  return 0;
}

static int PSolvePseudo(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int kmin, kmax, nb, nmax, nd;
  double xdf;
  if (argc < 2 || argc > 6) return -1;  
  if (argt[0] != NUMBER) return -1;
  if (argt[1] != NUMBER) return -1;
  kmax = atoi(argv[0]);
  nmax = atoi(argv[1]);
  kmin = 0;
  nb = 0;
  nd = 1;
  xdf = -1.0;
  if (argc > 2) {
    if (argt[2] != NUMBER) return -1;
    kmin = atoi(argv[2]);
    if (argc > 3) {
      if (argt[3] != NUMBER) return -1;
      nb = atoi(argv[3]);
      if (argc > 4) {
	if (argt[4] != NUMBER) return -1;
	nd = atoi(argv[4]);
	if (argc > 5) {
	  if (argt[5] != NUMBER) return -1;
	  xdf = atof(argv[5]);
	}
      }
    }
  }
  SolvePseudo(kmin, kmax, nb, nmax, nd, xdf);
  return 0;
}

static int PSetPotentialMode(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int m;
  double h, ih, h0, h1;

  h = 1E31;
  ih = 1E31;
  h0 = -1;
  h1 = -1;
  m = atoi(argv[0]);
  if (argc > 1) {
    h = atof(argv[1]);
    if (argc > 2) {
      ih = atof(argv[2]);
      if (argc > 3) {
	h0 = atof(argv[3]);
	if (argc > 4) {
	  h1 = atof(argv[4]);
	}
      }
    }
  }
  SetPotentialMode(m, h, ih, h0, h1);

  return 0;
}

static int PSetRadialGrid(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  double rmin, ratio, asym, qr;
  int maxrp;

  if (argc != 4 && argc != 5) return -1;
  
  maxrp = atoi(argv[0]);
  ratio = atof(argv[1]);
  asym = atof(argv[2]);
  rmin = atof(argv[3]);
  qr = -1;
  if (argc == 5) qr = atof(argv[4]);

  return SetRadialGrid(maxrp, ratio, asym, rmin, qr);
}

static int PSetRecPWLimits(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  int m1, m2;

  if (argc == 2) {
    if (argt[0] != NUMBER) return -1;
    m1 = atoi(argv[0]);
    if (argt[1] != NUMBER) return -1;
    m2 = atoi(argv[1]);
  } else {
    return -1;
  }
    
  SetRecPWLimits(m1, m2);

  return 0;
}

static int PSetRecSpectator(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  int n_spec, n_frozen, n_max;

  n_spec = 0;
  n_frozen = 0;
  n_max = 0;
  
  if (argc < 1 || argc > 3) return -1;
  n_spec = atoi(argv[0]);
  if (argc > 1) {
    n_frozen = atoi(argv[1]);
    if (argc > 2) {
      n_max = atoi(argv[2]);
    }
  }

  if (n_frozen == 0) n_frozen = n_spec;
  if (n_max == 0) n_max = 100;

  SetRecSpectator(n_max, n_frozen, n_spec);

  return 0;
}

static int PSetRRTEGrid(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int i, n;
  double xg[MAXNTE];
  int ng, err;
  double emin, emax;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  n = argc;
  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetRRTEGrid(ng, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
	if (ig[i] != NUMBER) return -1;
	xg[i] = atof(vg[i]);
	xg[i] /= HARTREE_EV;
	free(vg[i]);
      }
      err = SetRRTEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    err = SetRRTEGrid(ng, emin, emax);
  } else {
    return -1;
  }

  if (err < 0) return -1;

  return 0;
}

static int PSetSE(int argc, char *argv[], int argt[], 
		  ARRAY *variables) {
  int c, m, s, p;

  if (argc < 1 || argc > 4) return -1;
  c = atoi(argv[0]);
  m = -1;
  s = -1;
  p = -1;
  if (argc > 1) {
    m = atoi(argv[1]);
    if (argc > 2) {
      s = atoi(argv[2]);
      if (argc > 3) {
	p = atoi(argv[3]);
      }
    }
  }
  SetSE(c, m, s, p);
  
  return 0;
}

static int POptimizeModSE(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  double dr;
  int n, ka, ni;

  if (argc < 3 || argc > 4) return -1;
  n = atoi(argv[0]);
  ka = atoi(argv[1]);
  dr = atof(argv[2]);
  if (argc > 3) {
    ni = atoi(argv[3]);
  }
  OptimizeModSE(n, ka, dr, ni);
  return 0;
}

static int PSetModSE(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  double o0, o1, a, c0, c1, c;

  if (argc < 3 || argc > 6) return -1;
  o0 = atof(argv[0]);
  o1 = atof(argv[1]);
  a = atof(argv[2]);
  c0 = -1;
  c1 = -1;
  c = -1;
  if (argc > 3) {
    c0 = atof(argv[3]);
    if (argc > 4) {
      c1 = atof(argv[4]);
      if (argc > 5) {
	c = atof(argv[5]);
      }
    }
  }

  SetModSE(o0, o1, a, c0, c1, c);
  
  return 0;
}

static int PSetVP(int argc, char *argv[], int argt[], 
		  ARRAY *variables) {
  int c;

  if (argc != 1) return -1;
  c = atoi(argv[0]);

  SetVP(c);
  
  return 0;
}

static int PSetBreit(int argc, char *argv[], int argt[], 
		  ARRAY *variables) {
  int c, m, n, k;
  double x;

  if (argc < 1 || argc > 5) return -1;
  c = atoi(argv[0]);
  m = -1;
  n = -1;
  x = -1;
  k = 0;
  if (argc > 1) {
    m = atoi(argv[1]);
    if (argc > 2) {
      n = atoi(argv[2]);
      if (argc > 3) {
	x = atof(argv[3]);
	if (argc > 4) {
	  k = atoi(argv[4]);
	}
      }
    }
  }
  SetBreit(c, m, n, x, k);
  
  return 0;
}

static int PSetMS(int argc, char *argv[], int argt[], 
		  ARRAY *variables) {
  int c1, c2;

  if (argc != 2) return -1;
  c1 = atoi(argv[0]);
  c2 = atoi(argv[1]);

  SetMS(c1, c2);
  
  return 0;
}

static int PSetScreening(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int n_screen;
  int *screened_n = NULL;
  double screened_charge;
  int i, kl;
  char *v[MAXNARGS];
  int t[MAXNARGS];

  n_screen = 0;
  screened_charge = 1.0;
  kl = 1;
  
  if (argc < 1 || argc > 3) return -1;
  if (argt[0] != LIST && argt[0] != TUPLE) return -1;
  n_screen = DecodeArgs(argv[0], v, t, variables);
  if (argc > 1) {
    screened_charge = atof(argv[1]);
    if (argc > 2) {
      kl = atoi(argv[2]);
    }
  }
  
  screened_n = malloc(sizeof(int)*n_screen);
  for (i = 0; i < n_screen; i++) {
    if (t[i] != NUMBER) return -1;
    screened_n[i] = atoi(v[i]);
    free(v[i]);
  }
  
  SetScreening(n_screen, screened_n, screened_charge, kl);

  return 0;
}

static int PSetTransitionCut(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int c0, c;

  if (argc < 1) return -1;
  c = -1;
  c0 = atof(argv[0]);
  if (argc > 1) {
    c = atof(argv[1]);
  }

  SetTransitionCut(c0, c);
						  
  return 0;
}

static int PSetTransitionOptions(int argc, char *argv[], int argt[], 
				 ARRAY *variables) {
  int gauge, mode, max_m, max_e;

  max_e = 4;
  max_m = 4;

  if (argc < 2 || argc > 4) return -1;

  gauge = atoi(argv[0]);
  mode = atoi(argv[1]);
  if (argc > 2) {
    max_e = atoi(argv[2]);
    if (argc > 3) {
      max_m = atoi(argv[3]);
    }
  }

  SetTransitionOptions(gauge, mode, max_e, max_m);

  return 0;
}

static int PSetUsrCEGrid(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {  int n, ng, i, err;
  double xg[MAXNUSR];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0; 
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetUsrCEEGrid(ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
	xg[i] = atof(vg[i]);
	xg[i] /= HARTREE_EV;
	free(vg[i]);
      }
      err = SetUsrCEEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetUsrCEEGrid(ng, emin, emax, eth);
  } else {
    return -1;
  }
  
  if (err < 0) return -1;
  return 0;
}

static int PSetUsrCEGridType(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetUsrCEEGridType(atoi(argv[0]));
  return 0;
}

static int PSetUsrCIEGrid(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  int n, ng, i, err;
  double xg[MAXNUSR];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0; 
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetUsrCIEGrid(ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
	xg[i] = atof(vg[i]);
	xg[i] /= HARTREE_EV;
	free(vg[i]);
      }
      err = SetUsrCIEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetUsrCIEGrid(ng, emin, emax, eth);
  } else {
    return -1;
  }
  
  if (err < 0) return -1;
  return 0;
}

static int PSetUsrCIEGridType(int argc, char *argv[], int argt[], 
			      ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetUsrCIEGridType(atoi(argv[0]));
  return 0;
}

static int PSetUsrPEGrid(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int n, ng, i, err;
  double xg[MAXNUSR];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0; 
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetUsrPEGrid(ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
	xg[i] = atof(vg[i]);
	xg[i] /= HARTREE_EV;
	free(vg[i]);
      }
      err = SetUsrPEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetUsrPEGrid(ng, emin, emax, eth);
  } else {
    return -1;
  }
  
  if (err < 0) return -1;
  return 0;
}

static int PSetUsrPEGridType(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetUsrPEGridType(atoi(argv[0]));
  return 0;
}

static int PSolveBound(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int n, kappa;
  ORBITAL *orb;
  int k;
  
  if (argc != 2 || argt[0] != NUMBER || argt[1] != NUMBER) return -1;
  n = atoi(argv[0]);
  kappa = atoi(argv[1]);

  if (n <= 0) {
    printf("n must be greater than 0 for SolveBound\n");
    return -1;
  }

  k = OrbitalIndex(n, kappa, 0.0);
  if (k < 0) {
    printf("Fetal error in sloving dirac equation\n");
    return -1;
  }

  orb = GetOrbital(k);
  printf("Energy = %16.8E\n", orb->energy);
  
  return 0;
}

static int PSortLevels(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  if (argc != 0) return -1;

  SortLevels(0, 0, 0);
  
  return 0;
}

static int PTransitionMBPT(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  int m, n, nlow, *low, nup, *up;

  if (argc == 2) {
    m = atoi(argv[0]);
    n = atoi(argv[1]);
    TransitionMBPT(m, n);
  } else if (argc == 3) {
    nlow = DecodeGroupArgs(&low, 1, NULL, &(argv[1]), &(argt[1]), variables);
    nup = DecodeGroupArgs(&up, 1, NULL, &(argv[2]), &(argt[2]), variables);
    TRTableMBPT(argv[0], nlow, low, nup, up);
    if (nlow > 0) free(low);
    if (nup > 0) free(up);
  }

  return 0;
}

static int PStructureMBPT(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  int i, n, *s, n1, *ng1, n2, *ng2, nk, *nkm, kmax;
  int n3, *ng3, n4, *ng4;
  char *v[MAXNARGS], *gn;
  int t[MAXNARGS], nv;
  double d, c, e, f;

  if (argc == 1) {
    if (argt[0] != NUMBER && argt[0] != LIST) return -1;
    if (argt[0] == NUMBER) {
      f = atof(argv[0]);
      if (f < 0 || (f > 0 && f < 1)) {
	SetWarnMBPT(f, -1.0);
	return 0;
      } else {
	i = atoi(argv[0]);
	SetExtraMBPT(i);
      }
    } else {
      n1 = IntFromList(argv[0], argt[0], variables, &ng1);
      SetSymMBPT(n1, ng1);
      free(ng1);
    }
    return 0;
  }
  if (argc == 2) {
    if (argt[0] == NUMBER && argt[1] == NUMBER) {
      f = atof(argv[0]);
      d = atof(argv[1]);
      SetWarnMBPT(f, d);
      return 0;
    }
    if (argt[0] != STRING) return -1;
    if (argt[1] == NUMBER) {
      n2 = atoi(argv[1]);
      n1 = n2;
    } else if (argt[1] == LIST) {
      n3 = IntFromList(argv[1], argt[1], variables, &ng3);
      if (n3 == 0) {
	n2 = -1;
	n1 = -1;	
      } else if (n3 == 1) {
	n2 = ng3[0];
	n1 = n2;
      } else {
	n2 = ng3[0];
	n1 = ng3[1];
      }
      if (n3 > 0) free(ng3);
    }
    SetExcMBPT(n2, n1, argv[0]);
    return 0;
  }
  if ((argc == 3 || argc == 4 || argc == 5 || argc == 6)
      && argt[0] == NUMBER) {
    if (argt[1] != NUMBER) return -1;
    if (argt[2] != NUMBER) return -1;
    i = atoi(argv[0]);
    n3 = atoi(argv[1]);
    c = atof(argv[2]);
    d = -1.0;
    e = -1.0;
    f = -1.0;
    if (argc > 3 ) {
      if (argt[3] != NUMBER) return -1;
      d = atof(argv[3]);
      if (argc > 4) {
	if (argt[4] != NUMBER) return -1;
	e = atof(argv[4]);
	if (argc > 5) {
	  if (argt[5] != NUMBER) return -1;
	  f = atof(argv[5]);
	}
      }
    }      
    SetOptMBPT(i, n3, c, d, e, f);
    return 0;
  }
  /*
  if (argc == 10) {
    if (argt[1] != NUMBER) return -1;
    if (argt[2] != NUMBER) return -1;
    if (argt[3] != LIST) return -1;
    d = atof(argv[1]);
    c = atof(argv[2]);
    n = DecodeGroupArgs(&s, 1, &(argv[3]), &(argt[3]), variables);
    if (n <= 0) {
      printf("First configuration group does not exist\n");
      return -1;
    }
    kmax = atoi(argv[4]);
  
    n1 = IntFromList(argv[5], argt[5], variables, &ng1);
    n2 = IntFromList(argv[6], argt[6], variables, &ng2);
    n3 = IntFromList(argv[7], argt[7], variables, &ng3);
    n4 = IntFromList(argv[8], argt[8], variables, &ng4);
    gn = argv[9];

    StructureMBPT0(argv[0], d, c, n, s, kmax, 
		   n1, ng1, n2, ng2, n3, ng3, n4, ng4, gn);
    
    free(s);
    if (n1 > 0) free(ng1);
    if (n2 > 0) free(ng2);
    if (n3 > 0) free(ng3);
    if (n4 > 0) free(ng4);

    return 0;
  }
  */
  if (argc == 5) {
    if (argt[3] != LIST) return -1;
    if (argt[4] != NUMBER) return -1;
    n3 = atoi(argv[4]);
    n = DecodeGroupArgs(&s, 1, &n3, &(argv[3]), &(argt[3]), variables);
    if (n <= 0) return -1;
    if (argt[2] != LIST) return -1;
    n1 = DecodeArgs(argv[2], v, t, variables);
    for (i = 0; i < n1; i++) {
      if (t[i] != STRING) return -1;
    }
    if (n1 <= 0) return -1;
    StructureReadMBPT(argv[0], argv[1], n1, v, n, s, n3);
    free(s);
    for (i = 0; i < n1; i++) {
      free(v[i]);
    }
    
    return 0;
  }

  if (argc == 7 || argc == 9 || argc == 10) {
    char *hfn0, *hfn1;
    int nf = 0;
    hfn0 = NULL;
    hfn1 = argv[1];
    if (argt[1] == LIST) {
      nf = DecodeArgs(argv[1], v, t, variables);
      hfn1 = v[0];
      if (nf > 1) {
	hfn0 = v[1];
      }
    }
    if (argt[6] != NUMBER) return -1;
    n3 = atoi(argv[6]);
    if (argt[2] == LIST) {
      n = DecodeGroupArgs(&s, 1, &n3, &(argv[2]), &(argt[2]), variables);      
      if (n <= 0) {
	printf("First configuration group does not exist\n");
	return -1;
      }
    } else {
      n = 0;
      s = NULL;
    }
    if (argt[4] == LIST) {
      n1 = IntFromList(argv[4], argt[4], variables, &ng1);
    } else {
      n1 = atoi(argv[4]);
      ng1 = NULL;
    }
    if (argt[5] == LIST) {      
      n2 = IntFromList(argv[5], argt[5], variables, &ng2);
    } else {
      n2 = atoi(argv[5]);
      ng2 = NULL;
    }
    if (argt[3] == LIST) {
      nk = IntFromList(argv[3], argt[3], variables, &nkm);
    } else if (argt[3] == NUMBER) {
      nk = atoi(argv[3]) + 1;
      nkm = NULL;
    } else {
      return -1;
    }
    int icp = 0;
    int icpf = -1;
    int ncp = 0;
    if (argc >= 9) {
      ncp = atoi(argv[7]);
      icp = atoi(argv[8]);
      if (argc > 9) {
	icpf = atoi(argv[9]);
      }
    }
    StructureMBPT1(argv[0], hfn0, hfn1, n, s, nk, nkm, n1, ng1, n2, ng2, n3,
		   ncp, icp, icpf);
    free(s);
    if (nkm) free(nkm);
    for (i = 0; i < nf; i++) {
      free(v[i]);
    }
    return 0;
  } 
  
  return -1;
}

static int PCutMixing(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int nlev, n, *ilev, *kg;
  double c;
  
  nlev = 0;
  n = 0;
  if (argc < 2 || argc > 3) return -1;
  nlev = SelectLevels(&ilev, argv[0], argt[0], variables);
  if (nlev <= 0) goto DONE;
  n = DecodeGroupArgs(&kg, 1, NULL, &(argv[1]), &(argt[1]), variables);
  if (n <= 0) goto DONE;
  if (argc == 3) c = atof(argv[2]);
  else c = 0.0;

  CutMixing(nlev, ilev, n, kg, c);

 DONE:
  if (nlev > 0) free(ilev);
  if (n > 0) free(kg);

  return 0;  
}
static int PStructure(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int ng, ngp;
  int ip, i;
  int *kg, *kgp;
  int n;

  ng = 0;
  ngp = 0;
  n = argc;
  kgp = NULL;
  ip = 0;

  if (argc < 1) return -1;

  if (argt[0] == NUMBER) {      
    ip = atoi(argv[0]);
    if (argc == 1) {
      SetDiagMaxIter(ip, -1.0);
      return 0;
    }
    if (ip >= -1) {
      i = IntFromList(argv[1], argt[1], variables, &kg);
      SetSymmetry(ip, i, kg);
      free(kg);
    } else {
      if (argt[1] != NUMBER) return -1;
      double a = atof(argv[1]);
      double b = 0;
      double c = 0;
      if (ip >= -100) {
	if (argc > 2 && argt[2] == NUMBER) {
	  b = atof(argv[2]);
	}
	if (argc > 3 && argt[3] == NUMBER) {
	  c = atof(argv[3]);
	}
	SetPerturbThreshold(-ip, a, b, c);
      } else {
	SetDiagMaxIter(-ip-100, a);
      }
    }
    return 0;
  }
  char *hfn = NULL;
  if (n == 1) {
    if (argt[0] != STRING) return -1;
    ng = DecodeGroupArgs(&kg, 0, NULL, NULL, NULL, variables);
    if (ng < 0) return -1;
  } else {
    if (argt[1] == STRING) {
      hfn = argv[1];
      if (n > 5) return -1;
      if (n == 5) ip = atoi(argv[4]);
      if (argt[0] != STRING) return -1;
      if (n == 2) {
	ng = 0;
	ngp = 0;
	kg = NULL;
	kgp = NULL;
      } else {
	if (argt[2] != LIST && argt[2] != TUPLE) return -1;
	ng = DecodeGroupArgs(&kg, 1, NULL, &(argv[2]), &(argt[2]), variables);
	if (ng < 0) return -1;
	if (n >= 4) {
	  if (argt[3] != LIST && argt[3] != TUPLE) return -1;
	  ngp = DecodeGroupArgs(&kgp, 1, NULL,
				&(argv[3]), &(argt[3]), variables);
	}
      }
    } else {
      if (n > 4) return -1;
      if (n == 4) ip = atoi(argv[3]);		  
      if (argt[0] != STRING) return -1;
      if (argt[1] != LIST && argt[1] != TUPLE) return -1;
      ng = DecodeGroupArgs(&kg, 1, NULL, &(argv[1]), &(argt[1]), variables);
      if (ng < 0) return -1;
      if (n >= 3) {
	if (argt[2] != LIST && argt[2] != TUPLE) return -1;
	ngp = DecodeGroupArgs(&kgp, 1, NULL, &(argv[2]), &(argt[2]), variables);
      }
    }
  }

  return SolveStructure(argv[0], hfn, ng, kg, ngp, kgp, ip);
}

static int PSetUTA(int argc, char *argv[], int argt[], 
		   ARRAY *variables) {
  int m, mci;

  if (argc == 1) {
    if (argt[0] != NUMBER) return -1;
    m = atoi(argv[0]);
    mci = 1;
  } else if (argc == 2) {
    if (argt[0] != NUMBER || argt[1] != NUMBER) return -1;
    m = atoi(argv[0]);
    mci = atoi(argv[1]);
  }

  SetUTA(m, mci);
  
  return 0;
}

static int PSetTRF(int argc, char *argv[], int argt[], 
		   ARRAY *variables) {
  
  if (argc != 1 || argt[0] != NUMBER) return -1;
  
  SetTRF(atoi(argv[0]));
  
  return 0;
}

static int PCoulombBethe(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  double z, te, e1;

  if (argc != 4) return -1;
  z = atof(argv[1]);
  te = atof(argv[2]);
  e1 = atof(argv[3]);

  CoulombBethe(argv[0], z, te, e1);

  return 0;
}

static int PTestAngular(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  TestAngular();  
  return 0;
}

static int PTestIntegrate(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  TestIntegrate();  
  return 0;
}

static int PReportMultiStats(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  ReportMultiStats();
  return 0;
}

static int PTestMyArray(int argc, char *argv[], int argt[], 
			ARRAY *variables) { 
  ARRAY a;
  double d;
  double *b;
  MULTI ma;
  int k[3] = {101, 2550, 333};
  int block[3] = {10, 20, 5};
  int i, j, m;

  /*
  ArrayInit(&a, sizeof(double), 100);
  d = 0.1;
  m = 100000;
  printf("> ");
  scanf("%d", &i);
  for (i = 0; i < m; i++) {
    ArraySet(&a, i, &d, InitDoubleData);
  }
  printf("> ");
  scanf("%d", &i);

  b = (double *) ArrayGet(&a, 100);
  printf("%f ", *b);
  b = (double *) ArrayGet(&a, 200);
  printf("%f \n", *b);

  ArrayFree(&a, 0);
  printf("> ");
  scanf("%d", &i);
  */
  MultiInit(&ma, sizeof(double), 3, block, "test_ma");
  printf("%d %d\n", ma.esize, ma.ndim);
  for (i = 9; i < 15; i++) {
    for (j = 0; j < 30; j++) {
      k[0] = i;
      k[1] = j;
      k[2] = 20;
      b = (double *) MultiSet(&ma, k, NULL, NULL, InitDoubleData, NULL);
      *b = 0.2;
      double *b1 = (double *) MultiSet(&ma, k, NULL, NULL,
				       InitDoubleData, NULL);
    }
  }

  MultiFreeData(&ma, NULL);
  for (i = 9; i < 15; i++) {
    for (j = 0; j < 30; j++) {
      k[0] = i;
      k[1] = j;
      k[2] = 20;
      b = (double *) MultiSet(&ma, k, NULL, NULL, InitDoubleData, NULL);
      *b = 0.2;
      double *b1 = (double *) MultiSet(&ma, k, NULL, NULL,
				       InitDoubleData, NULL);
    }
  }
  
  //printf("> ");
  //scanf("%d", &i);
  MultiFree(&ma, NULL);

  //printf("> ");
  //scanf("%d", &i);

  return 0;
}

static int PPrepAngular(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int nlow, nup, *low, *up;

  nlow = 0; 
  nup = 0;
  low = NULL;
  up = NULL;

  if (argc < 1 || argc > 2) return -1;
  nlow = SelectLevels(&low, argv[0], argt[0], variables);
  if (nlow <= 0) return -1;
  if (argc == 2) {
    nup = SelectLevels(&up, argv[1], argt[1], variables);
    if (nup <= 0) {
      free(low);
      return -1;
    }
  }
  PrepAngular(nlow, low, nup, up);

  return 0;
}

static int PElectronDensity(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  int n, *ilev, t;
  if (argc < 2 || argc > 3) return -1;
  t = 1;
  n = SelectLevels(&ilev, argv[1], argt[1], variables);
  if (argc == 3) {
    t = atoi(argv[2]);
  }
  if (n > 0) {
    ElectronDensity(argv[0], n, ilev, t);
    free(ilev);
  }
  return 0;
}

static int PExpectationValue(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int n, *ilev, t;
  double a;
  if (argc < 3 || argc > 5) return -1;
  if (argt[0] != STRING) return -1;
  if (argt[1] != STRING) return -1;
  t = 1;
  a = 0;
  if (argc > 3) {
    a = atof(argv[3]);
    if (argc > 4) {
      t = atoi(argv[4]);
    }
  }
  n = SelectLevels(&ilev, argv[2], argt[2], variables);
  if (n > 0) {
    ExpectationValue(argv[0], argv[1], n, ilev, a, t);
    free(ilev);
  }
  return 0;
}

static int PTransitionTable(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  int n, m;
  int nlow, nup, *low, *up;
  
  nlow = 0; 
  nup = 0;
  low = NULL;
  up = NULL;
  m = 0;

  n = argc;

  if (n == 1) {
    if (argt[0] != STRING) return -1;
    SaveTransition(nlow, low, nup, up, argv[0], m);
  } else if (n == 2) {
    if (argt[0] != STRING || argt[1] != NUMBER) return -1;
    m = atoi(argv[1]);
    SaveTransition(nlow, low, nup, up, argv[0], m);
  } else if (n == 3) {
    if (argt[0] != STRING) return -1;
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    if (nlow <= 0) return -1;
    nup = SelectLevels(&up, argv[2], argt[2], variables);
    if (nup <= 0) return -1;
    SaveTransition(nlow, low, nup, up, argv[0], m);
    free(low);
    free(up);
  } else if (n == 4) {
    if (argt[0] != STRING || argt[3] != NUMBER) return -1;
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    if (nlow <= 0) return -1;
    nup = SelectLevels(&up, argv[2], argt[2], variables);
    if (nup <= 0) return -1;
    m = atoi(argv[3]);
    SaveTransition(nlow, low, nup, up, argv[0], m);
    free(low);
    free(up);
  } else {
    return -1;
  }
    
  return 0;
}

static int PWaveFuncTable(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  int k, n;
  double e;

  if (argc != 3 && argc != 4) return -1;
  
  n = atoi(argv[1]);
  k = atoi(argv[2]);
  if (argc == 4) {
    e = atof(argv[3]);
  } else {
    e = 0.0;
  }

  WaveFuncTable(argv[0], n, k, e);

  return 0;
}

static int PSetConfigEnergyMode(int argc, char *argv[], int argt[], 
				ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetConfigEnergyMode(m);
  return 0;
}

static int PSetOptimizeMaxIter(int argc, char *argv[], int argt[], 
			       ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetOptimizeMaxIter(m);
  return 0;
}

static int PSetOptimizeStabilizer(int argc, char *argv[], int argt[], 
				  ARRAY *variables) {
  double m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atof(argv[0]);
  SetOptimizeStabilizer(m);
  return 0;
}

static int PSetOptimizePrint(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetOptimizePrint(m);
  return 0;
}

static int PSetOptimizeTolerance(int argc, char *argv[], int argt[], 
				 ARRAY *variables) {
  double m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atof(argv[0]);
  SetOptimizeTolerance(m);
  return 0;
}

static int PSetCELQR(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCELQR(m);
  return 0;
}

static int PSetCELMax(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCELMax(m);
  return 0;
}

static int PSetCELCB(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCELCB(m);
  return 0;
}

static int PSetCETol(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  double m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atof(argv[0]);
  SetCETol(m);
  return 0;
}

static int PSetCILQR(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCILQR(m);
  return 0;
}

static int PSetCILMax(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCILMax(m);
  return 0;
}

static int PSetCILMaxEject(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCILMaxEject(m);
  return 0;
}

static int PSetCILCB(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCILCB(m);
  return 0;
}

static int PSetCITol(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  double m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atof(argv[0]);
  SetCITol(m);
  return 0;
}

static int PSetTransitionMode(int argc, char *argv[], int argt[], 
			      ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetTransitionMode(m);
  return 0;
}

static int PSetTransitionGauge(int argc, char *argv[], int argt[], 
			       ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetTransitionGauge(m);
  return 0;
}

static int PSetTransitionMaxE(int argc, char *argv[], int argt[], 
			      ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetTransitionMaxE(m);
  return 0;
}

static int PSetTransitionMaxM(int argc, char *argv[], int argt[], 
			      ARRAY *variables) {
  int m;
  
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetTransitionMaxM(m);
  return 0;
}

static int PInterpCross(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  char *ifn, *ofn;
  int i, negy, i0, i1, mp;
  double *egy;
  
  mp = 1;
  if (argc < 5 || argc > 6) return -1;
  ifn = argv[0];
  ofn = argv[1];
  i0 = atoi(argv[2]);
  i1 = atoi(argv[3]);
  negy = DoubleFromList(argv[4], argt[4], variables, &egy);
  if (argc > 5) {
    mp = atoi(argv[5]);
  }
  if (negy > 0) {
    InterpCross(ifn, ofn, i0, i1, negy, egy, mp);
    free(egy);
  }
  return 0;
}

static int PMaxwellRate(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  char *ifn, *ofn;
  int i, nt, i0, i1;
  double *temp;
  
  if (argc != 5) return -1;
  ifn = argv[0];
  ofn = argv[1];
  i0 = atoi(argv[2]);
  i1 = atoi(argv[3]);
  nt = DoubleFromList(argv[4], argt[4], variables, &temp);

  if (nt > 0) {
    MaxwellRate(ifn, ofn, i0, i1, nt, temp);
    free(temp);
  }
  return 0;
}

static int PAsymmetry(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int mx;
  
  if (argc < 2) return -1;
  if (argc == 3) mx = atoi(argv[2]);
  else mx = 1;
  
  SaveAsymmetry(argv[0], argv[1], mx);
  
  return 0;
}

static int PRadialOverlaps(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  int kappa;

  if (argc != 2 || argt[0] != STRING) return -1;

  if (argc == 2) kappa = atoi(argv[1]);
  else kappa = -1;

  RadialOverlaps(argv[0], kappa);
  
  return 0;
}

static int PFreezeOrbital(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  if (argc < 1 || argc > 2 || argt[0] != STRING) return -1;
  int m = -1;
  if (argc > 1) m = atoi(argv[1]);
  FreezeOrbital(argv[0], m);
  return 0;
}

static int PSetBoundary(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int nmax;
  double bqp, p;

  bqp = 0.0;
  p = -1.0;
  if (argc > 1) {
    p = atof(argv[1]);    
    if (argc > 2) {
      bqp = atof(argv[2]);
    }
  }

  nmax = atoi(argv[0]);
  
  return SetBoundary(nmax, p, bqp);
}

static int PRMatrixExpansion(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int m;
  double d, a, r;

  if (argc < 1) return -1;
  
  m = atoi(argv[0]);
  d = 1E-3;
  a = 1E-4;
  r = 0.0;
  if (argc > 1) {
    r = atof(argv[1]);
    if (argc > 2) {
      d = atof(argv[2]);
      if (argc > 3) {
	a = atof(argv[3]);
      }
    }
  }

  RMatrixExpansion(m, d, a, r);

  return 0;
}

static int PRMatrixNBatch(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  int m;

  if (argc != 1) return -1;
  m = atoi(argv[0]);
  RMatrixNBatch(m);
  
  return 0;
}

static int PRMatrixFMode(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int m;

  if (argc != 1) return -1;
  m = atoi(argv[0]);
  RMatrixFMode(m);
  
  return 0;
}

static int PRMatrixConvert(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int m;

  if (argc != 3) return -1;
  m = atoi(argv[2]);
  
  RMatrixConvert(argv[0], argv[1], m);
  
  return 0;
}

static int PRMatrixNMultipoles(int argc, char *argv[], int argt[], 
			       ARRAY *variables) {
  int m;

  if (argc != 1) return -1;
  m = atoi(argv[0]);
  RMatrixNMultipoles(m);
  
  return 0;
}

static int PRMatrixBoundary(int argc, char *argv[], int argt[], 
			    ARRAY *variables) {
  double r0, r1, b;
  
  if (argc != 3) return -1;
  r0 = atof(argv[0]);
  r1 = atof(argv[1]);
  b = atof(argv[2]);
  
  RMatrixBoundary(r0, r1, b);
  
  return 0;
}

static int PRMatrixBasis(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int kmax, nb;

  if (argc != 3) return -1;
  kmax = atoi(argv[1]);
  nb = atoi(argv[2]);
  
  RMatrixBasis(argv[0], kmax, nb);
  
  return 0;
}

static int PRMatrixTargets(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  int nt, *kt, nc, *kc;
  
  if (argc < 1 || argc > 2) return -1;
  if (argt[0] != LIST && argt[0] != TUPLE) return -1;
  
  nt = DecodeGroupArgs(&kt, 1, NULL, &(argv[0]), &(argt[0]), variables);
  if (nt < 0) return -1;
  nc = 0;
  kc = NULL;
  if (argc == 2) {
    nc = DecodeGroupArgs(&kc, 1, NULL, &(argv[1]), &(argt[1]), variables);
    if (nc < 0) nc = 0;
  }
    
  RMatrixTargets(nt, kt, nc, kc);

  if (nt > 0) free(kt);
  if (nc > 0) free(kc);

  return 0;
}

static int PRMatrixSurface(int argc, char *argv[], int argt[], 
			   ARRAY *variables) {
  if (argc != 1) return -1;
  if (argt[0] != STRING) return -1;

  RMatrixSurface(argv[0]);
  
  return 0;
}

static int PSetSlaterCut(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int k0, k1;

  if (argc != 2) return -1;
  
  k0 = atoi(argv[0]);
  k1 = atoi(argv[1]);
  
  SetSlaterCut(k0, k1);
  
  return 0;
}

static int PTestRMatrix(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int m;
  double e;

  if (argc != 5) return -1;
  e = atof(argv[0]);
  m = atoi(argv[1]);
  
  TestRMatrix(e, m, argv[2], argv[3], argv[4]);
  
  return 0;
}

static int PRMatrixRefine(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  int n, m;
  double r;

  if (argc < 1 || argc > 3) return -1;
  n = atoi(argv[0]);
  m = -1;
  r = -1;
  if (argc > 1) {
    m = atoi(argv[1]);
    if (argc > 2) {
      r = atof(argv[2]);
    }
  }
  RMatrixRefine(n, m, r);
  return 0;
}

static int PRMatrixCE(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  double emin, emax, de;
  int np, m, i, mb;
  char *v0[MAXNARGS], *v1[MAXNARGS];
  int t0[MAXNARGS], t1[MAXNARGS];
  
  if (argc < 6 || argc > 8) return -1;
  emin = atof(argv[3]);
  emax = atof(argv[4]);
  de = atof(argv[5]);
  m= 0;
  mb = 1;
  if (argc >= 7) {
    m = atoi(argv[6]);
    if (argc > 7) {
      mb = atoi(argv[7]);
    }
  }

  np = DecodeArgs(argv[1], v0, t0, variables);
  if (DecodeArgs(argv[2], v1, t1, variables) != np) {
    return -1;
  }
  
  RMatrixCE(argv[0], np, v0, v1, emin, emax, de, m, mb);

  for (i = 0; i < np; i++) {
    free(v0[i]);
    free(v1[i]);
  }
  return 0;
}

static int PSetCEPWFile(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  if (argc != 1) return -1;

  SetCEPWFile(argv[0]);
  
  return 0;
}

static int PPropogateDirection(int argc, char *argv[], int argt[], 
			       ARRAY *variables) {
  int m;

  if (argc != 1) return -1;
  m = atoi(argv[0]);

  PropogateDirection(m);

  return 0;
}
  
static int PAppendTable(int argc, char *argv[], int argt[], 
			ARRAY *variables) {  
  if (argc != 1) return -1;
  AppendTable(argv[0]);
  
  return 0;
}

static int PJoinTable(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  if (argc != 3) return -1;
  
  JoinTable(argv[0], argv[1], argv[2]);
  
  return 0;
}

static int PModifyTable(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  
  if (argc == 3) {  
    ModifyTable(argv[0], argv[1], argv[2], NULL);
  } else if (argc == 4) {
    ModifyTable(argv[0], argv[1], argv[2], argv[3]);
  } else {
    return -1;
  }
  
  return 0;
}

static int PLimitArray(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  int m;
  double n;

  if (argc != 2) return -1;
  m = atoi(argv[0]);
  n = atof(argv[1]);

  LimitArrayRadial(m, n);
  
  return 0;
}

static int PSetFields(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int m;
  double b, e, a;

  if (argc < 3 || argc > 4) return -1;
  m = 0;
  b = atof(argv[0]);
  e = atof(argv[1]);
  a = atof(argv[2]);
  if (argc > 3) {
    m = atoi(argv[3]);
  }

  SetFields(b, e, a, m);
  
  return 0;
}

static int PStructureEB(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int n, *ilev;

  if (argc != 2 || argt[0] != STRING || argt[1] != LIST) return -1;

  n = SelectLevels(&ilev, argv[1], argt[1], variables);
  if (n <= 0) return -1;
  
  StructureEB(argv[0], n, ilev);
  free(ilev);

  return 0;
}

static int PTransitionTableEB(int argc, char *argv[], int argt[], 
			      ARRAY *variables) {
  int m, nlow, *low, nup, *up;

  if (argc < 3 || argc > 4) return -1;
  
  if (argc == 4) m = atoi(argv[3]);
  else m = -1;
  
  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  if (nlow <= 0) {  
    printf("cannot determine levels in lower\n");
    return -1;
  }
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  if (nup <= 0) {
    printf("cannot determine levels in upper\n");
    return -1;
  }
  
  SaveTransitionEB(nlow, low, nup, up, argv[0], m);
  free(low);
  free(up);

  return 0;
}

static int PPolarizeCoeff(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  int i0, i1;

  if (argc < 2 || argc > 4) return -1;
  i0 = -1;
  i1 = -1;
  if (argc > 2) {
    i0 = atoi(argv[2]);
    if (argc > 3) {
      i1 = atoi(argv[3]);
    }
  }

  PolarizeCoeff(argv[0], argv[1], i0, i1);

  return 0;
}

static int PCETableEB(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int nlow, nup, *low, *up, m;

  if (argc < 3 || argc > 4) return -1;
  if (argc == 4) m = atoi(argv[3]);
  else m = 0;
  
  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  if (nlow <= 0) return -1;
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  if (nup <= 0) return -1;
  
  if (m == 0) {
    SaveExcitationEB(nlow, low, nup, up, argv[0]);
  } else {
    SaveExcitationEBD(nlow, low, nup, up, argv[0]);
  }
  
  free(low);
  free(up);
  
  return 0;
}

static int PCoulMultip(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  double z, te, e1;
  int k, q0, q1, m, ierr;
  char *fn;

  if (argc < 7 || argc > 8) return -1;
  fn = argv[0];
  z = atof(argv[1]);
  te = atof(argv[2]);
  e1 = atof(argv[3]);
  k = atoi(argv[4]);
  q0 = atoi(argv[5]);
  q1 = atoi(argv[6]);
  m = 1;
  if (argc > 7) m = atoi(argv[7]);
  
  ierr = CoulombMultip(fn, z, te, e1, k, q0, q1, m);
  
  return ierr;
}

static int PSlaterCoeff(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  int nlev, *ilev, na, nb, i, *n, *kappa;
  double *nq;
  SHELL *sa, *sb;
  
  if (argc != 4) return -1;
  if (argt[1] != LIST) return -1;
  if (argt[2] != STRING) return -1;
  if (argt[3] != STRING) return -1;
    
  nlev = SelectLevels(&ilev, argv[1], argt[1], variables);
  na = GetAverageConfigFromString(&n, &kappa, &nq, argv[2]);
  sa = malloc(sizeof(SHELL)*na);
  for (i = 0; i < na; i++) {
    sa[i].n = n[i];
    sa[i].kappa = kappa[i];
  }
  if (na > 0) {
    free(n);
    free(kappa);
    free(nq);
  }
  nb = GetAverageConfigFromString(&n, &kappa, &nq, argv[3]);
  sb = malloc(sizeof(SHELL)*nb);
  for (i = 0; i < nb; i++) {
    sb[i].n = n[i];
    sb[i].kappa = kappa[i];
  }
  if (nb > 0) {
    free(n);
    free(kappa);
    free(nq);
  }


  if (nlev > 0 && na > 0 && nb > 0) {
    SlaterCoeff(argv[0], nlev, ilev, na, sa, nb, sb);
    free(ilev);
    free(sa);
    free(sb);
  }
  
  return 0;
}

static int PGeneralizedMoment(int argc, char *argv[], int argt[], 
			      ARRAY *variables) {
  int n0, k0, n1, k1, m;
  double e1;
  
  m = atoi(argv[1]);
  n0 = atoi(argv[2]);
  k0 = atoi(argv[3]);
  n1 = atoi(argv[4]);
  k1 = atoi(argv[5]);
  if (argc == 7) e1 = atof(argv[6]);
  else e1 = 0.0;
  PrintGeneralizedMoments(argv[0], m, n0, k0, n1, k1, e1);
  
  return 0;
}

static int PPrintQED(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  PrintQED();
  return 0;
}

static int PPrintNucleus(int argc, char *argv[], int argt[], 
			 ARRAY *variables) {
  int m;
  char *fn;
  m = 0;
  fn = NULL;
  if (argc > 0) {
    m = atoi(argv[0]);
    if (argc > 1) {
      fn = argv[1];
    }
  }
  PrintNucleus(m, fn);
  return 0;
} 
 
static int PSavePotential(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  char *fn;
  POTENTIAL *p;

  if (argc != 1) return -1;
  fn = argv[0];
  
  p = RadialPotential();
  SavePotential(fn, p);

  return 0;
}
 
static int PRestorePotential(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  char *fn;
  POTENTIAL *p;

  if (argc != 1) return -1;
  fn = argv[0];
  
  p = RadialPotential();
  RestorePotential(fn, p);

  return 0;
} 
 
static int PModifyPotential(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  char *fn;
  POTENTIAL *p;

  if (argc != 1) return -1;
  fn = argv[0];
  
  p = RadialPotential();
  ModifyPotential(fn, p);

  return 0;
}

static int PWallTime(int argc, char *argv[], int argt[], 
		     ARRAY *variables) {
  int m = 0;
  if (argc < 1) return -1;
  if (argc > 1) {
    m = atoi(argv[1]);
  }
  PrintWallTime(argv[0], m);
  return 0;
}

static int PInitializeMPI(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
#ifdef USE_MPI
  int n = -1;
  if (argc > 0) {
    n = atoi(argv[0]);
  }
  InitializeMPI(n, 0);
#endif
  return 0;
}

static int PMPIRank(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  int n, k;
  k = MPIRank(&n);
  MPrintf(-1, "%d of %d\n", k, n);
  return 0;
}

static int PMemUsed(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  MPrintf(-1, "mem used %g\n", msize());
  return 0;
}

static int PFinalizeMPI(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
#if USE_MPI == 1
  FinalizeMPI();
#endif
  return 0;
}

static int PSetOrbMap(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int k = 0, n0 = 0, n1 = 0, n2 = 0;
  if (argc > 0) {
    k = atoi(argv[0]);
    if (argc > 1) {
      n0 = atoi(argv[1]);
      if (argc > 2) {
	n1 = atoi(argv[2]);
	if (argc > 3) {
	  n2 = atoi(argv[3]);
	}
      }
    }
  }
  SetOrbMap(k, n0, n1, n2);
  return 0;
}

static METHOD methods[] = {
  {"GeneralizedMoment", PGeneralizedMoment, METH_VARARGS},
  {"SlaterCoeff", PSlaterCoeff, METH_VARARGS},
  {"PropogateDirection", PPropogateDirection, METH_VARARGS}, 
  {"SetUTA", PSetUTA, METH_VARARGS}, 
  {"SetTRF", PSetTRF, METH_VARARGS}, 
  {"SetCEPWFile", PSetCEPWFile, METH_VARARGS}, 
  {"AppendTable", PAppendTable, METH_VARARGS}, 
  {"JoinTable", PJoinTable, METH_VARARGS}, 
  {"ModifyTable", PModifyTable, METH_VARARGS},
  {"LimitArray", PLimitArray, METH_VARARGS},
  {"RMatrixExpansion", PRMatrixExpansion, METH_VARARGS}, 
  {"RMatrixNBatch", PRMatrixNBatch, METH_VARARGS}, 
  {"RMatrixFMode", PRMatrixFMode, METH_VARARGS}, 
  {"RMatrixConvert", PRMatrixConvert, METH_VARARGS}, 
  {"RMatrixNMultipoles", PRMatrixNMultipoles, METH_VARARGS}, 
  {"TestRMatrix", PTestRMatrix, METH_VARARGS}, 
  {"RMatrixRefine", PRMatrixRefine, METH_VARARGS}, 
  {"RMatrixCE", PRMatrixCE, METH_VARARGS}, 
  {"SetSlaterCut", PSetSlaterCut, METH_VARARGS}, 
  {"RMatrixBoundary", PRMatrixBoundary, METH_VARARGS}, 
  {"RMatrixBasis", PRMatrixBasis, METH_VARARGS}, 
  {"RMatrixTargets", PRMatrixTargets, METH_VARARGS}, 
  {"RMatrixSurface", PRMatrixSurface, METH_VARARGS}, 
  {"Print", PPrint, METH_VARARGS},
  {"AddConfig", PAddConfig, METH_VARARGS},
  {"AITable", PAITable, METH_VARARGS},
  {"AITableMSub", PAITableMSub, METH_VARARGS},
  {"Asymmetry", PAsymmetry, METH_VARARGS},
  {"AvgConfig", PAvgConfig, METH_VARARGS},
  {"BasisTable", PBasisTable, METH_VARARGS},
  {"CECross", PInterpCross, METH_VARARGS},
  {"CERate", PMaxwellRate, METH_VARARGS},
  {"InterpCross", PInterpCross, METH_VARARGS},
  {"MaxwellRate", PMaxwellRate, METH_VARARGS},
  {"CETable", PCETable, METH_VARARGS},
  {"CETableMSub", PCETableMSub, METH_VARARGS},
  {"CheckEndian", PCheckEndian, METH_VARARGS},
  {"CITable", PCITable, METH_VARARGS},
  {"CITableMSub", PCITableMSub, METH_VARARGS},
  {"ClearLevelTable", PClearLevelTable, METH_VARARGS},
  {"ClearOrbitalTable", PClearOrbitalTable, METH_VARARGS},
  {"Closed", PClosed, METH_VARARGS},
  {"Config", PConfig, METH_VARARGS},
  {"ReadConfig", PReadConfig, METH_VARARGS},
  {"CutMixing", PCutMixing, METH_VARARGS},
  {"RemoveConfig", PRemoveConfig, METH_VARARGS},
  {"ListConfig", PListConfig, METH_VARARGS},
  {"GetConfigNR", PGetConfigNR, METH_VARARGS},
  {"ConfigEnergy", PConfigEnergy, METH_VARARGS},
  {"AdjustEnergy", PAdjustEnergy, METH_VARARGS},
  {"CorrectEnergy", PCorrectEnergy, METH_VARARGS},
  {"Exit", PExit, METH_VARARGS},
  {"FreeExcitationQk", PFreeExcitationQk, METH_VARARGS},
  {"FreeIonizationQk", PFreeIonizationQk, METH_VARARGS},
  {"FreeMemENTable", PFreeMemENTable, METH_VARARGS},
  {"FreeMultipole", PFreeMultipole, METH_VARARGS},
  {"FreeSlater", PFreeSlater, METH_VARARGS},
  {"FreeResidual", PFreeResidual, METH_VARARGS},
  {"FreeRecPk", PFreeRecPk, METH_VARARGS},
  {"FreeRecQk", PFreeRecQk, METH_VARARGS},
  {"GetPotential", PGetPotential, METH_VARARGS},
  {"Info", PInfo, METH_VARARGS},
  {"MemENTable", PMemENTable, METH_VARARGS},
  {"StructureMBPT", PStructureMBPT, METH_VARARGS},
  {"TransitionMBPT", PTransitionMBPT, METH_VARARGS},
  {"OptimizeRadial", POptimizeRadial, METH_VARARGS},
  {"PrepAngular", PPrepAngular, METH_VARARGS},
  {"Pause", PPause, METH_VARARGS},
  {"RadialOverlaps", PRadialOverlaps, METH_VARARGS},
  {"FreezeOrbital", PFreezeOrbital, METH_VARARGS},
  {"RefineRadial", PRefineRadial, METH_VARARGS},
  {"PrintMemInfo", PPrintMemInfo, METH_VARARGS},
  {"PrintTable", PPrintTable, METH_VARARGS},
  {"RecStates", PRecStates, METH_VARARGS},
  {"ReinitConfig", PReinitConfig, METH_VARARGS},
  {"ReinitRecouple", PReinitRecouple, METH_VARARGS},
  {"ReinitRadial", PReinitRadial, METH_VARARGS},
  {"ReinitDBase", PReinitDBase, METH_VARARGS},
  {"ReinitStructure", PReinitStructure, METH_VARARGS},
  {"ReinitExcitation", PReinitExcitation, METH_VARARGS},
  {"ReinitRecombination", PReinitRecombination, METH_VARARGS},
  {"ReinitIonization", PReinitIonization, METH_VARARGS},
  {"Reinit", PReinit, METH_VARARGS},
  {"RRTable", PRRTable, METH_VARARGS},
  {"RRMultipole", PRRMultipole, METH_VARARGS},
  {"SetAICut", PSetAICut, METH_VARARGS},
  {"SetAngZOptions", PSetAngZOptions, METH_VARARGS},
  {"SetAngZCut", PSetAngZCut, METH_VARARGS},
  {"SetCILevel", PSetCILevel, METH_VARARGS},
  {"SetBoundary", PSetBoundary, METH_VARARGS},
  {"SetMixCut", PSetMixCut, METH_VARARGS},
  {"SetAtom", PSetAtom, METH_VARARGS},
  {"SetExtraPotential", PSetExtraPotential, METH_VARARGS},
  {"SetAvgConfig", PSetAvgConfig, METH_VARARGS},
  {"SetCEGrid", PSetCEGrid, METH_VARARGS},
  {"SetTEGrid", PSetTEGrid, METH_VARARGS},
  {"SetAngleGrid", PSetAngleGrid, METH_VARARGS},
  {"SetCEBorn", PSetCEBorn, METH_VARARGS},
  {"SetCIBorn", PSetCIBorn, METH_VARARGS},
  {"SetBornFormFactor", PSetBornFormFactor, METH_VARARGS},
  {"SetBornMass", PSetBornMass, METH_VARARGS},
  {"SetCEPWOptions", PSetCEPWOptions, METH_VARARGS},
  {"SetCEPWGrid", PSetCEPWGrid, METH_VARARGS},
  {"SetCEQkMode", PSetCEQkMode, METH_VARARGS},
  {"SetCIEGrid", PSetCIEGrid, METH_VARARGS},
  {"SetCIEGridLimits", PSetCIEGridLimits, METH_VARARGS},
  {"SetIEGrid", PSetIEGrid, METH_VARARGS},
  {"SetCIQkMode", PSetCIQkMode, METH_VARARGS},
  {"SetCIPWOptions", PSetCIPWOptions, METH_VARARGS},
  {"SetCIPWGrid", PSetCIPWGrid, METH_VARARGS},
  {"SetHydrogenicNL", PSetHydrogenicNL, METH_VARARGS},
  {"SetMaxRank", PSetMaxRank, METH_VARARGS},
  {"SetOptimizeControl", PSetOptimizeControl, METH_VARARGS},
  {"SetPEGrid", PSetPEGrid, METH_VARARGS},
  {"SetPEGridLimits", PSetPEGridLimits, METH_VARARGS},
  {"SetRadialGrid", PSetRadialGrid, METH_VARARGS},
  {"SolvePseudo", PSolvePseudo, METH_VARARGS},
  {"SetPotentialMode", PSetPotentialMode, METH_VARARGS},
  {"SetRecPWLimits", PSetRecPWLimits, METH_VARARGS},
  {"SetRecPWOptions", PSetRecPWOptions, METH_VARARGS},
  {"SetRecQkMode", PSetRecQkMode, METH_VARARGS},
  {"SetRecSpectator", PSetRecSpectator, METH_VARARGS},
  {"SetRRTEGrid", PSetRRTEGrid, METH_VARARGS},
  {"SetScreening", PSetScreening, METH_VARARGS},
  {"SetSE", PSetSE, METH_VARARGS},
  {"SetModSE", PSetModSE, METH_VARARGS},
  {"OptimizeModSE", PSetModSE, METH_VARARGS},
  {"SetMS", PSetMS, METH_VARARGS},
  {"SetVP", PSetVP, METH_VARARGS},
  {"SetBreit", PSetBreit, METH_VARARGS},
  {"SetTransitionCut", PSetTransitionCut, METH_VARARGS},
  {"SetTransitionOptions", PSetTransitionOptions, METH_VARARGS},
  {"SetUsrCEGrid", PSetUsrCEGrid, METH_VARARGS},
  {"SetUsrCEGridType", PSetUsrCEGridType, METH_VARARGS},
  {"SetCEGridLimits", PSetCEGridLimits, METH_VARARGS},
  {"SetCEGridType", PSetCEGridType, METH_VARARGS},
  {"SetCEPWGridType", PSetCEPWGridType, METH_VARARGS},
  {"SetUsrCIEGrid", PSetUsrCIEGrid, METH_VARARGS},
  {"SetUsrCIEGridType", PSetUsrCIEGridType, METH_VARARGS},
  {"SetUsrPEGrid", PSetUsrPEGrid, METH_VARARGS},
  {"SetUsrPEGridType", PSetUsrPEGridType, METH_VARARGS},
  {"SolveBound", PSolveBound, METH_VARARGS},
  {"SortLevels", PSortLevels, METH_VARARGS},
  {"Structure", PStructure, METH_VARARGS},
  {"CoulombBethe", PCoulombBethe, METH_VARARGS}, 
  {"TestAngular", PTestAngular, METH_VARARGS}, 
  {"TestIntegrate", PTestIntegrate, METH_VARARGS}, 
  {"TestMyArray", PTestMyArray, METH_VARARGS},   
  {"ReportMultiStats", PReportMultiStats, METH_VARARGS},   
  {"ElectronDensity", PElectronDensity, METH_VARARGS},  
  {"ExpectationValue", PExpectationValue, METH_VARARGS},  
  {"TransitionTable", PTransitionTable, METH_VARARGS},  
  {"TRTable", PTransitionTable, METH_VARARGS},  
  {"WaveFuncTable", PWaveFuncTable, METH_VARARGS},
  {"SetConfigEnergyMode", PSetConfigEnergyMode, METH_VARARGS},
  {"SetOptimizeMaxIter", PSetOptimizeMaxIter, METH_VARARGS},
  {"SetOptimizeStabilizer", PSetOptimizeStabilizer, METH_VARARGS},
  {"SetOptimizePrint", PSetOptimizePrint, METH_VARARGS},
  {"SetOptimizeTolerance", PSetOptimizeTolerance, METH_VARARGS},
  {"SetCELQR", PSetCELQR, METH_VARARGS},
  {"SetCELMax", PSetCELMax, METH_VARARGS},
  {"SetCELCB", PSetCELCB, METH_VARARGS},
  {"SetCETol", PSetCETol, METH_VARARGS},
  {"SetCILQR", PSetCILQR, METH_VARARGS},
  {"SetCILMax", PSetCILMax, METH_VARARGS},
  {"SetCILMaxEject", PSetCILMaxEject, METH_VARARGS},
  {"SetCILCB", PSetCILCB, METH_VARARGS},
  {"SetCITol", PSetCITol, METH_VARARGS},
  {"SetTransitionMode", PSetTransitionMode, METH_VARARGS},
  {"SetTransitionGauge", PSetTransitionGauge, METH_VARARGS},
  {"SetTransitionMaxE", PSetTransitionMaxE, METH_VARARGS},
  {"SetTransitionMaxM", PSetTransitionMaxM, METH_VARARGS}, 
  {"SetFields", PSetFields, METH_VARARGS},    
  {"TRTableEB", PTransitionTableEB, METH_VARARGS}, 
  {"CETableEB", PCETableEB, METH_VARARGS},
  {"StructureEB", PStructureEB, METH_VARARGS},
  {"PolarizeCoeff", PPolarizeCoeff, METH_VARARGS}, 
  {"CoulMultipole", PCoulMultip, METH_VARARGS}, 
  {"PrintQED", PPrintQED, METH_VARARGS},
  {"PrintNucleus", PPrintNucleus, METH_VARARGS},
  {"SavePotential", PSavePotential, METH_VARARGS},
  {"RestorePotential", PRestorePotential, METH_VARARGS},
  {"ModifyPotential", PModifyPotential, METH_VARARGS},
  {"WallTime", PWallTime, METH_VARARGS},
  {"InitializeMPI", PInitializeMPI, METH_VARARGS},
  {"MPIRank", PMPIRank, METH_VARARGS},
  {"MemUsed", PMemUsed, METH_VARARGS},
  {"FinalizeMPI", PFinalizeMPI, METH_VARARGS},
  {"SetOrbMap", PSetOrbMap, METH_VARARGS},
  {"", NULL, METH_VARARGS}
};
 
int main(int argc, char *argv[]) {
  int i;
  FILE *f;

#if PMALLOC_CHECK == 1
  pmalloc_open();
#endif
  
  if (InitFac() < 0) {
    printf("initialization failed\n");
    exit(1);
  }

  SetModName("fac");

  if (argc == 1) {
    EvalFile(stdin, 1, methods);
  } else {
    for (i = 1; i < argc; i++) {
      f = fopen(argv[i], "r");
      if (!f) {
	printf("Cannot open file %s, Skipping\n", argv[i]);
	continue;
      }
      EvalFile(f, 0, methods);
    }
  }

#if PMALLOC_CHECK == 1
  pmalloc_check();
#endif

  return 0;
}


     
