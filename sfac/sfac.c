static char *rcsid="$Id: sfac.c,v 1.7 2002/01/14 23:19:51 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "init.h"
#include "stoken.h"

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
      printf(", ");
    }
  }
  if (argc > 0) printf("\n");
  fflush(stdout);
  return 0;
}

static int DecodeGroupArgs(int **kg, int n, char *argv[], int argt[],
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
    for (i = 0; i < ng; i++) {
      if (t[i] != STRING) {
	printf("argument must be a group name\n");
	return -1;
      }
      s = v[i];
      k = GroupExists(s);
      
      if (k < 0) {
	free((*kg));
	printf("group does not exist\n");
	return -1;
      }

      (*kg)[i] = k;
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
  int ig, nlevels;
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

  n = DecodeArgs(argv, v, at, variables);
  nv = n;
  if (n > 0) {
    if (at[0] == STRING) {
      ng = DecodeGroupArgs(&kg, n, v, at, variables);
      if (ng <= 0) {
	rv = -1;
	goto END;
      }
      nlevels = GetNumLevels();
      (*t) = malloc(sizeof(int)*nlevels);
      k = 0;
      for (j = 0; j < nlevels; j++) {
	lev = GetLevel(j);
	im = lev->basis[0];
	sym = GetSymmetry(lev->pj);
	s = (STATE *) ArrayGet(&(sym->states), im);
	ig = s->kgroup;
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
      ng = DecodeGroupArgs(&kg, 1, v, at, variables);
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
	  im = lev->basis[0];
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


static char _closed_shells[128] = "";
static int PClosed(int argc, char *argv[], int argt[], ARRAY *variables) {
  CONFIG *cfg;
  int i, j, kappa, jj, kl, n, nq, ncfg;
  char js, *p;
  char s[16], st[16];
  int ns, k;

  if (argc == 0) _closed_shells[0] = '\0';
  for (i = 0; i < argc; i++) {
    if (argt[i] != STRING) return -1;
    ns = StrSplit(argv[i], ' ');
    p = argv[i];
    for (k = 0; k < ns; k++) {
      while (*p == ' ') p++;
      ncfg = GetConfigFromString(&cfg, p);
      for (j = ncfg-1; j >= 0; j--) {
	if (cfg[j].n_shells != 1) return -1;
	n = (cfg[j].shells)[0].n;
	kappa = (cfg[j].shells)[0].kappa;
	GetJLFromKappa(kappa, &jj, &kl);
	nq = jj + 1;
	if (jj > kl) js = '+';
	else js = '-';
	kl = kl/2;
	SpecSymbol(s, kl);
	sprintf(st, "%d%s%c%d ", n, s, js, nq);
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
  
static int PConfig(int argc, char *argv[], int argt[], ARRAY *variables) {
  CONFIG *cfg;
  static char gname[GROUP_NAME_LEN] = "_all_";
  int i, j, k, t, ncfg;
  char scfg[1280];
  
  k = -2;
  for (i = 0; i < argc; i++) {
    if (argt[i] == KEYWORD) {
      if (strcmp(argv[i], "group") != 0) return -1;
      if (i > argc-2) return -1;
      if (argt[i+1] != STRING) return -1;
      k = i;
    }
  }

  if (k >= 0) strncpy(gname, argv[k+1], GROUP_NAME_LEN);

  for (i = 0; i < argc; i++) {
    if (i == k || i == k+1) continue;
    if (argt[i] != STRING) return -1;
    strncpy(scfg, _closed_shells, 128);
    strncat(scfg, argv[i], 1280);
    ncfg = GetConfigFromString(&cfg, scfg);
    for (j = 0; j < ncfg; j++) {
      if (Couple(cfg+j) < 0) return -1;
      t = GroupIndex(gname);
      if (t < 0) return -1;
      if (AddConfigToList(t, cfg+j) < 0) return -1;
    }   
    if (ncfg > 0) free(cfg);
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
  int nlow, *low, nup, *up, c;

  if (argc != 3 && argc != 4) return -1;
  if (argt[0] != STRING) return -1;
  
  if (argc == 4) {
    if (argt[3] != NUMBER) return -1;
    c = atoi(argv[3]);
  } else c = 0;
  
  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  SaveAI(nlow, low, nup, up, argv[0], c);
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  return 0;
}

static int PBasisTable(int argc, char *argv[], int argt[], ARRAY *variables) {
  
  if (argc != 1 || argt[0] != STRING) return -1;
  GetBasisTable(argv[0]);
  
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

  free(low);
  free(up);
    
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

static int PCorrectEnergy(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  int n, k[MAX_ENERGY_CORRECTION];
  double e[MAX_ENERGY_CORRECTION];
  int i, ie, ii;
  FILE *f;
  char *iv[MAXNARGS], *ev[MAXNARGS];
  int it[MAXNARGS], et[MAXNARGS];

  ii = 0; 
  ie = 0;
  if (argc == 1) {
    if (argt[0] != STRING) {
      printf("A single argument for CorrectEnergy must be a file name\n");
      return -1;
    }
    f = fopen(argv[0], "r");
    n = -1;
    while (1) {
      n++; 
      if (n == MAX_ENERGY_CORRECTION) {
	printf("Maximum # of levels for energy correction reached\n");
	printf("Ignoring corrections after n = %d\n", n);
	break;
      } 
      if (fscanf(f, "%d%lf\n", k+n, e+n) == EOF) break;
      e[n] /= HARTREE_EV;
    }
    fclose(f);
  } else if (argc == 2) {
    if (argt[0] != LIST || argt[1] != LIST) {
      printf("The two arguments for CorrectEnergy must be two Lists\n");
      return -1;
    }
    ii = DecodeArgs(argv[0], iv, it, variables);
    ie = DecodeArgs(argv[1], ev, et, variables);
    if (ii != ie) return -1;
    n = ii;
    for (i = 0; i < n; i++) {
      if (it[i] != NUMBER || et[i] != NUMBER) return -1;
      k[i] = atoi(iv[i]);
      e[i] = atof(ev[i]);
      e[i] /= HARTREE_EV;
    }
  } else {
    return -1;
  }

  CorrectEnergy(n, k, e);

  for (i = 0; i < ii; i++) free(iv[i]);
  for (i = 0; i < ie; i++) free(ev[i]);

  return 0;
}

static int PDRTable(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  int nf, *f, na, *a, nb, *b, ng, *g, c;
  
  if (argc == 7 && argt[6] == NUMBER) {
    c = atoi(argv[6]);
  } else if (argc == 6) {
    c = 0;
  } else {
    return -1;
  }

  if (argt[0] != STRING || argt[1] != STRING) return -1;

  nf = SelectLevels(&f, argv[2], argt[2], variables);
  na = SelectLevels(&a, argv[3], argt[3], variables);
  nb = SelectLevels(&b, argv[4], argt[4], variables);
  ng = SelectLevels(&g, argv[5], argt[5], variables);

  SaveDR(nf, f, na, a, nb, b, ng, g, argv[0], argv[1], c);

  if (nf > 0) free(f);
  if (na > 0) free(a);
  if (nb > 0) free(b);
  if (ng > 0) free(g);

  return 0;
}

static int PExit(int argc, char *argv[], int argt[], ARRAY *variables) {
  if (argc != 0) return -1;
  exit(0);
}

static int PFreeAngZ(int argc, char *argv[], int argt[], ARRAY *variables) {
  int i, m;
  int n, *kg;
  
  m = -1;
  if (argc == 0) {
    FreeAngZ(-1, m);
  } else if (argc > 0) {
    if (argc > 2) return -1;
    if (argt[0] != LIST) return -1;
    n = DecodeGroupArgs(&kg, 1, argv, argt, variables);
    if (argc == 2) {
      if (argt[1] != NUMBER) return -1;
      m = atoi(argv[1]);
    }
    for (i = 0; i < n; i++) {
      FreeAngZ(kg[i], m);
    }
  }
    
  return 0;
}

static int PFreeExcitationPk(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int ie;
  
  if (argc == 0) ie = -1;
  else if (argc == 1 && argt[0] == NUMBER) ie = atoi(argv[0]);
  else return -1;

  FreeExcitationPk(ie);
  return 0;
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

static int PFreeRecAngZ(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  if (argc != 0) return -1;
  FreeRecAngZ();
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

static int PLoadIonizationQk(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int n;
  
  if (argc == 0) {
    LoadCIRadialQkIntegrated(-1, NULL);
  } else if (argc == 1) {
    n = atoi(argv[0]);
    LoadCIRadialQkIntegrated(n, NULL);
  } else if (argc == 2) {
    n = atoi(argv[0]);
    LoadCIRadialQkIntegrated(n, argv[1]);
  } else {
    return -1;
  }
    
  return 0;
}

static int PMemENTable(int argc, char *argv[], int argt[], 
		       ARRAY *variables) {
  
  if (argc != 1 || argt[0] != STRING) return -1;
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
    ng = DecodeGroupArgs(&kg, argc, argv, argt, variables);
    if (ng < 0) return -1;
  } else {
    ng = DecodeGroupArgs(&kg, 1, argv, argt, variables);
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
  if (OptimizeRadial(ng, kg, weight) < 0) {
    if (kg) free(kg);
    if (weight) free(weight);
    return -1;
  }
  if (weight) free(weight);
  if (kg) free(kg);

  for (i = 0; i < ni; i++) free(vw[i]);

  return 0;
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

static int PPrepIonizationQk(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int nz, na, nn, nte;
  double emin, emax;
  int *n, i;
  double *z, *a;
  char *vz[MAXNARGS], *va[MAXNARGS], *vn[MAXNARGS];
  int iz[MAXNARGS], ia[MAXNARGS], in[MAXNARGS];

  nz = 0;
  na = 0;
  nn = 0;
  nte = -1;
  emin = -1.0;
  emax = -1.0;

  if (argc < 1 || argc > 7) return -1;

  if (argc > 1) {
    if (argt[1] != LIST) return -1;
    nz = DecodeArgs(argv[1], vz, iz, variables);
    if (argc > 2) {
      if (argt[2] != LIST) return -1;
      na = DecodeArgs(argv[2], va, ia, variables);
      if (argc > 3) {
	if (argt[3] != LIST) return -1;
	nn = DecodeArgs(argv[3], vn, in, variables);
	if (argc > 4) {
	  if (argt[4] != NUMBER) return -1;
	  nte = atoi(argv[4]);
	  if (argc > 5) {
	    emin = atof(argv[5]);
	    if (argc > 6) {
	      emax = atof(argv[6]);
	    }
	  }
	}
      }
    }
  }

  if (nz > 0) {
    z = (double *) malloc(sizeof(double)*nz);
    for (i = 0; i < nz; i++) {
      z[i] = atof(vz[i]);
    }
  } else {
    nz = 5;
    z = (double *) malloc(sizeof(double)*nz);
    z[0] = 10;
    z[1] = 30;
    z[2] = 50;
    z[3] = 70;
    z[4] = 90;
  }

  if (na > 0) {
    a = (double *) malloc(sizeof(double)*na);
    for (i = 0; i < na; i++) {
      a[i] = atof(va[i]);
    }
  } else { 
    na = 5;
    a = (double *) malloc(sizeof(double)*na);
    a[0] = 0.1;
    a[1] = 0.3;
    a[2] = 0.5;
    a[3] = 0.7;
    a[4] = 0.9;
  }

  if (nn > 0) {
    n = (int *) malloc(sizeof(int)*nn);
    for (i = 0; i < nn; i++) {    
      n[i] = atoi(vn[i]);
    }
  } else { 
    nn = 5;
    n = (int *) malloc(sizeof(int)*nn);
    n[0] = 1;
    n[1] = 2;
    n[2] = 3;
    n[3] = 4;
    n[4] = 5;
  }

  for (i = 0; i < nz; i++) free(vz[i]);
  for (i = 0; i < na; i++) free(va[i]);
  for (i = 0; i < nn; i++) free(vn[i]);
  
  PrepCIRadialQkIntegrated(nz, z, na, a, nn, n, nte, emin, emax, argv[0]);
 
  free(z);
  free(a);
  free(n);

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
  
  v = 0;
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
  ng = DecodeGroupArgs(&kg, 1, &(argv[1]), &(argt[1]), variables);
  if (ng <= 0) return -1;
  n = atoi(argv[2]);
  RecStates(n, ng, kg, argv[0]);
  
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

  return 0;
}

static int PSaveIonizationQk(int argc, char *argv[], int argt[], 
			     ARRAY *variables) {
  int n;
  if (argc != 2) return -1;
  
  n = atoi(argv[0]);
  SaveCIRadialQkIntegrated(n, argv[0]);
  
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
  double c;
  
  if (argc != 1 || argt[0] != NUMBER) return -1;
  c = atof(argv[0]);
  SetMixCut(c);
  
  return 0;
}

static int PSetAtom(int argc, char *argv[], int argt[], 
		    ARRAY *variables) {
  double z, mass;

  mass = 0.0;
  z = 0.0;

  if (argc < 1 || argt[0] != STRING || argc > 3) return -1;
  if (argc > 1) {
    z = atof(argv[1]);
    if (argc > 2) {
      mass = atof(argv[2]);
    }
  }
  
  if (SetAtom(argv[0], z, mass) < 0) return -1;
  
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

static int PSetCEFormat(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  if (argc != 1) return -1;
  SetCEFormat(atoi(argv[0]));

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

  qr = 0;
  max = 500;
  kl_cb = 100;
  tol = 5E-2;
  
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
    
static int PSetCIFormat(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetCIFormat(atoi(argv[0]));

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

  qr = 0;
  max = 500;
  max_1 = 8;
  kl_cb = 50;
  tol = 5E-2;
  
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
      if (strcasecmp(argv[0], "exact") == 0) m = QK_EXACT;
      else if (strcasecmp(argv[0], "interpolate") == 0) m = QK_INTERPOLATE;
      else if (strcasecmp(argv[0],"fit") == 0) m = QK_FIT;
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
  double tol;

  iprint = 0;
  
  if (argc != 2 && argc != 3) return -1;
  tol = atof(argv[0]);
  maxiter = atoi(argv[1]);
  if (argc == 3) iprint = atoi(argv[2]);
  
  SetOptimizeControl(tol, maxiter, iprint);

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

static int PSetRadialGrid(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  double rmin, rmax;

  rmin = -1.0;
  rmax = -1.0;
  
  if (argc > 2) return -1;
  
  if (argc > 0) {
    rmin = atof(argv[0]);
    if (argc > 1) {
      rmax = atof(argv[1]);
    }
  }

  SetRadialGrid(rmin, rmax);

  return 0;
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
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetTransitionCut(atof(argv[0]));
						  
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

  SortLevels(0, 0);
  
  return 0;
}

static int PStructure(int argc, char *argv[], int argt[], 
		      ARRAY *variables) {
  int i, k, ng, ngp, ns, nlevels;
  int *kg, *kgp;
  int n;

  ng = 0;
  ngp = 0;
  n = argc;
  kgp = NULL;

  if (n == 1) {
    if (argt[0] != STRING) return -1;
    ng = DecodeGroupArgs(&kg, 0, NULL, NULL, variables);
    if (ng < 0) return -1;
  } else {
    if (n > 3) return -1;
    if (argt[0] != STRING) return -1;
    if (argt[1] != LIST && argt[1] != TUPLE) return -1;
    ng = DecodeGroupArgs(&kg, 1, &(argv[1]), &(argt[1]), variables);
    if (ng < 0) return -1;
    if (n == 3) {
      if (argt[2] != LIST && argt[2] != TUPLE) return -1;
      ngp = DecodeGroupArgs(&kg, 1, &(argv[2]), &(argt[2]), variables);
    }
  }
  
  if (ngp < 0) return -1;

  nlevels = GetNumLevels();
  ns = MAX_SYMMETRIES;  
  for (i = 0; i < ns; i++) {
    k = ConstructHamilton(i, ng, kg, ngp, kgp);
    if (k < 0) continue;
    if (DiagnolizeHamilton() < 0) return -1;
    AddToLevels();
  }
  SortLevels(nlevels, -1);
  SaveLevels(argv[0], nlevels, -1);

  if (ng > 0) free(kg);
  if (ngp > 0) free(kgp);
  
  return 0;
}

static int PTestCoulomb(int argc, char *argv[], int argt[], 
			ARRAY *variables) {
  if (argc != 1 || argt[0] != STRING) return -1;
  TestCoulomb(argv[0]);

  return 0;
}

static int PTestIntegrate(int argc, char *argv[], int argt[], 
			  ARRAY *variables) {
  if (argc != 1 || argt[0] != STRING) return -1;
  TestIntegrate(argv[0]);  
  return 0;
}

static int PTestMyArray(int argc, char *argv[], int argt[], 
			ARRAY *variables) { 
  ARRAY a;
  double d;
  double *b;
  MULTI ma;
  int k[3] = {101, 2550, 333};
  int block[3] = {10, 20, 50};
  int i, j, m;
 
  ArrayInit(&a, sizeof(double), 100);
  d = 0.1;
  m = 100000;
  printf("> ");
  scanf("%d", &i);
  for (i = 0; i < m; i++) {
    ArraySet(&a, i, &d);
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

  MultiInit(&ma, sizeof(double), 3, block);
  printf("%d %d\n", ma.esize, ma.ndim);
  for (i = 9; i < 15; i++) {
    for (j = 0; j < m; j++) {
      k[0] = i;
      k[1] = j;
      k[2] = 20;	
      b = (double *) MultiSet(&ma, k, NULL);
      *b = 0.2;
      b = (double *) MultiGet(&ma, k);
    }
  }

  printf("> ");
  scanf("%d", &i);
  MultiFree(&ma, NULL);

  printf("> ");
  scanf("%d", &i);

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
  m = -1;

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

static METHOD methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"AddConfig", PAddConfig, METH_VARARGS},
  {"AITable", PAITable, METH_VARARGS},
  {"AvgConfig", PAvgConfig, METH_VARARGS},
  {"BasisTable", PBasisTable, METH_VARARGS},
  {"CETable", PCETable, METH_VARARGS},
  {"CETableMSub", PCETableMSub, METH_VARARGS},
  {"CITable", PCITable, METH_VARARGS},
  {"ClearLevelTable", PClearLevelTable, METH_VARARGS},
  {"ClearOrbitalTable", PClearOrbitalTable, METH_VARARGS},
  {"Closed", PClosed, METH_VARARGS},
  {"Config", PConfig, METH_VARARGS},
  {"CorrectEnergy", PCorrectEnergy, METH_VARARGS},
  {"DRTable", PDRTable, METH_VARARGS},
  {"Exit", PExit, METH_VARARGS},
  {"FreeAngZ", PFreeAngZ, METH_VARARGS},
  {"FreeExcitationPk", PFreeExcitationPk, METH_VARARGS},
  {"FreeExcitationQk", PFreeExcitationQk, METH_VARARGS},
  {"FreeIonizationQk", PFreeIonizationQk, METH_VARARGS},
  {"FreeMemENTable", PFreeMemENTable, METH_VARARGS},
  {"FreeMultipole", PFreeMultipole, METH_VARARGS},
  {"FreeSlater", PFreeSlater, METH_VARARGS},
  {"FreeResidual", PFreeResidual, METH_VARARGS},
  {"FreeRecPk", PFreeRecPk, METH_VARARGS},
  {"FreeRecQk", PFreeRecQk, METH_VARARGS},
  {"FreeRecAngZ", PFreeRecAngZ, METH_VARARGS},
  {"GetPotential", PGetPotential, METH_VARARGS},
  {"Info", PInfo, METH_VARARGS},
  {"LoadIonizationQk", PLoadIonizationQk, METH_VARARGS},
  {"MemENTable", PMemENTable, METH_VARARGS},
  {"OptimizeRadial", POptimizeRadial, METH_VARARGS},
  {"Pause", PPause, METH_VARARGS},
  {"PrepIonizationQk", PPrepIonizationQk, METH_VARARGS},
  {"PrintMemInfo", PPrintMemInfo, METH_VARARGS},
  {"PrintTable", PPrintTable, METH_VARARGS},
  {"RecStates", PRecStates, METH_VARARGS},
  {"Reinit", PReinit, METH_VARARGS},
  {"RRTable", PRRTable, METH_VARARGS},
  {"SaveIonizationQk", PSaveIonizationQk, METH_VARARGS},
  {"SetAICut", PSetAICut, METH_VARARGS},
  {"SetAngZOptions", PSetAngZOptions, METH_VARARGS},
  {"SetAngZCut", PSetAngZCut, METH_VARARGS},
  {"SetMixCut", PSetMixCut, METH_VARARGS},
  {"SetAtom", PSetAtom, METH_VARARGS},
  {"SetAvgConfig", PSetAvgConfig, METH_VARARGS},
  {"SetCEFormat", PSetCEFormat, METH_VARARGS},
  {"SetCEGrid", PSetCEGrid, METH_VARARGS},
  {"SetTEGrid", PSetTEGrid, METH_VARARGS},
  {"SetCEPWOptions", PSetCEPWOptions, METH_VARARGS},
  {"SetCEPWGrid", PSetCEPWGrid, METH_VARARGS},
  {"SetCEQkMode", PSetCEQkMode, METH_VARARGS},
  {"SetCIFormat", PSetCIFormat, METH_VARARGS},
  {"SetCIEGrid", PSetCIEGrid, METH_VARARGS},
  {"SetCIEGridLimits", PSetCIEGridLimits, METH_VARARGS},
  {"SetIEGrid", PSetIEGrid, METH_VARARGS},
  {"SetCIQkMode", PSetCIQkMode, METH_VARARGS},
  {"SetCIPWOptions", PSetCIPWOptions, METH_VARARGS},
  {"SetCIPWGrid", PSetCIPWGrid, METH_VARARGS},
  {"SetMaxRank", PSetMaxRank, METH_VARARGS},
  {"SetOptimizeControl", PSetOptimizeControl, METH_VARARGS},
  {"SetPEGrid", PSetPEGrid, METH_VARARGS},
  {"SetPEGridLimits", PSetPEGridLimits, METH_VARARGS},
  {"SetRadialGrid", PSetRadialGrid, METH_VARARGS},
  {"SetRecPWLimits", PSetRecPWLimits, METH_VARARGS},
  {"SetRecPWOptions", PSetRecPWOptions, METH_VARARGS},
  {"SetRecQkMode", PSetRecQkMode, METH_VARARGS},
  {"SetRecSpectator", PSetRecSpectator, METH_VARARGS},
  {"SetRRTEGrid", PSetRRTEGrid, METH_VARARGS},
  {"SetScreening", PSetScreening, METH_VARARGS},
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
  {"TestCoulomb", PTestCoulomb, METH_VARARGS}, 
  {"TestIntegrate", PTestIntegrate, METH_VARARGS}, 
  {"TestMyArray", PTestMyArray, METH_VARARGS},   
  {"TransitionTable", PTransitionTable, METH_VARARGS},  
  {"WaveFuncTable", PWaveFuncTable, METH_VARARGS},  
  {"", NULL, METH_VARARGS}
};
 
int main(int argc, char *argv[]) {
  int i;
  FILE *f;

  if (InitFac() < 0) {
    printf("initialization failed\n");
    exit(1);
  }

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

#ifdef PMALLOC_CHECK
  pmalloc_check();
#endif

  return 0;
}


     
