#include "Python.h"
#include <stdio.h>
#include <string.h>

#include "init.h"

static char *rcsid="$Id: fac.c,v 1.14 2002/01/14 23:19:47 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static PyObject *ErrorObject;
static PyObject *PFACVERSION;
static PyObject *SPECSYMBOL;
static PyObject *ATOMICSYMBOL;
static PyObject *ATOMICMASS;
static PyObject *QKMODE;

static FILE *sfac_file = NULL;

#define onError(message) {PyErr_SetString(ErrorObject, message);}

static void SFACStatement(char *func, PyObject *args, PyObject *kargs) {
  int i, n, nargs;
  PyObject *sargs;
  PyObject *klist;
  PyObject *kvar;
  PyObject *p, *q;
  char *s1, *s2;
  
  fprintf(sfac_file, "%s", func);
  nargs = PyTuple_Size(args);
  sargs = PyObject_Str(args);
  s1 = PyString_AsString(sargs);
  n = strlen(s1);
  if (nargs == 1) {
    n = n-2;
  } else {
    n = n-1;
  }
  for (i = 0; i < n; i++) {
    fprintf(sfac_file, "%c", s1[i]);
  }
  if (kargs) {
    klist = PyDict_Items(kargs);
    n = PyList_Size(klist);
    for (i = 0; i < n; i++) {
      if (nargs > 0 || i > 0) fprintf(sfac_file, ", ");
      p = PyList_GetItem(klist, i);
      q = PyTuple_GetItem(p, 0);
      s2 = PyString_AsString(q);
      fprintf(sfac_file, "%s=", s2);
      q = PyTuple_GetItem(p, 1);
      kvar = PyObject_Str(q);
      s2 = PyString_AsString(kvar);
      if (PyString_Check(q)) {
	fprintf(sfac_file, "'%s'", s2);
      } else {
	fprintf(sfac_file, "%s", s2);
      }
      Py_XDECREF(kvar);
    }
    Py_XDECREF(klist);
  }

  fprintf(sfac_file, ")\n");
    
  Py_XDECREF(sargs);

  return;
}
  
static PyObject *PConvertToSFAC(PyObject *self, PyObject *args) {
  char *fn;

  sfac_file = NULL;
  if (!PyArg_ParseTuple(args, "|s", &fn)) return NULL;

  if (fn) {
    sfac_file = fopen(fn, "w");
    if (sfac_file == NULL) return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}    
  
static PyObject *PCloseSFAC(PyObject *self, PyObject *args) {

  fclose(sfac_file);
  sfac_file = NULL;

  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PSetOptimizeControl(PyObject *self, PyObject *args) {
  int maxiter;
  double tol; 
  int iprint;
  
  if (sfac_file) {
    SFACStatement("SetOptimizeControl", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
    
  iprint = 0;
  if (!PyArg_ParseTuple(args, "di|i", &tol, &maxiter, &iprint))
    return NULL;
  SetOptimizeControl(tol, maxiter, iprint);

  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PSetScreening(PyObject *self, PyObject *args) {
  int n_screen;
  int *screened_n = NULL;
  double screened_charge;
  int i, kl;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("SetScreening", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n_screen = 0;
  screened_charge = 1.0;  
  kl = 1;
  
  if (!PyArg_ParseTuple(args, "O|di", &p, &screened_charge, &kl)) return NULL;
  if (screened_charge <= 0) {
    printf("screened charge must be positive\n");
    return NULL;
  }
  if (!PyList_Check(p) && !PyTuple_Check(p)) {
    printf("Screened n must be in a List or a Tuple\n");
    return NULL;
  }
  n_screen = PySequence_Length(p); 
  screened_n = malloc(sizeof(int)*n_screen);
  for (i = 0; i < n_screen; i++) {
    q = PySequence_GetItem(p, i);
    if (!PyInt_Check(q)) {
      printf("Screened n must be integers\n");
      free(screened_n);
      Py_DECREF(q);
      return NULL;
    }
    screened_n[i] = PyInt_AsLong(q);
    Py_DECREF(q);
  }
  
  SetScreening(n_screen, screened_n, screened_charge, kl);
  Py_INCREF(Py_None);
  return Py_None;
}  

/** Convert a Python config to a C struct **/
static int ConfigPythonToC(PyObject *python_cfg, CONFIG **cfg) {
  int n_shells;
  int i, j, k, m;
  PyObject *python_shell;

  SHELL *shells;

  shells = NULL;
  /** python_cfg should be a python list **/
  if (!PyList_Check(python_cfg)) goto ERROR;

  (*cfg) = malloc(sizeof(CONFIG));
  if ((*cfg) == NULL) goto ERROR;
  n_shells = PyList_Size(python_cfg);
  (*cfg)->n_shells = n_shells;
  (*cfg)->shells = malloc(sizeof(SHELL) * (*cfg)->n_shells);
  shells = (*cfg)->shells;
  if (shells == NULL) goto ERROR;
  
  for (i = 0, m = n_shells-1; i < n_shells; i++, m--) {
    python_shell = PyList_GetItem(python_cfg, i);
    if (!PyTuple_Check(python_shell)) goto ERROR;
    if (PyTuple_Size(python_shell) != 4) goto ERROR;
    
    shells[m].n = PyInt_AsLong(PyTuple_GetItem(python_shell, 0));
    k = PyInt_AsLong(PyTuple_GetItem(python_shell, 1));
    j = PyInt_AsLong(PyTuple_GetItem(python_shell, 2));
    if (j > 0) k = -(k+1);
    shells[m].kappa = k;
    shells[m].nq = PyInt_AsLong(PyTuple_GetItem(python_shell, 3));
  }

  return 0;
  
 ERROR:
  onError("error in conversion");
  if (shells) free(shells);
  if (cfg) free(cfg);
  return -1;
}

static char _closed_shells[128] = "";
static PyObject *PClosed(PyObject *self, PyObject *args) {
  CONFIG *cfg;
  PyObject *q;
  int i, j, kappa, jj, kl, n, nq, ncfg;
  char js, *p, argv[512];
  char s[16], st[16];
  int ns, k;
  int argc;

  if (sfac_file) {
    SFACStatement("Closed", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  argc = PyTuple_Size(args);
  if (argc == 0) _closed_shells[0] = '\0';
  for (i = 0; i < argc; i++) {
    q = PyTuple_GetItem(args, i);
    if (!PyString_Check(q)) return NULL;
    p = PyString_AsString(q);
    strncpy(argv, p, 512);
    ns = StrSplit(argv, ' ');
    p = argv;
    for (k = 0; k < ns; k++) {
      while (*p == ' ') p++;
      ncfg = GetConfigFromString(&cfg, p);
      for (j = ncfg-1; j >= 0; j--) {
	if (cfg[j].n_shells != 1) return NULL;	
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
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PConfig(PyObject *self, PyObject *args, PyObject *keywds) {
  CONFIG *cfg;
  PyObject *q;
  static char gname[GROUP_NAME_LEN] = "_all_";
  int i, j, t, ncfg;
  char scfg[1280], *p;
  int argc;

  if (sfac_file) {
    SFACStatement("Config", args, keywds);
    Py_INCREF(Py_None);
    return Py_None;
  }

  q = PyDict_GetItemString(keywds, "group");
  if (q) {
    if (!PyString_Check(q)) return NULL;
    p = PyString_AsString(q);
    strncpy(gname, p, GROUP_NAME_LEN);
  }

  argc = PyTuple_Size(args);
  
  for (i = 0; i < argc; i++) {   
    q = PyTuple_GetItem(args, i);
    if (!PyString_Check(q)) return NULL;
    p = PyString_AsString(q);
    strncpy(scfg, _closed_shells, 128);
    strncat(scfg, p, 1280);
    ncfg = GetConfigFromString(&cfg, scfg);

    for (j = 0; j < ncfg; j++) {
      if (Couple(cfg+j) < 0) return NULL;
      t = GroupIndex(gname);
      if (t < 0) return NULL;
      if (AddConfigToList(t, cfg+j) < 0) return NULL;
    }   
    if (ncfg > 0) free(cfg);
  }

  Py_INCREF(Py_None);
  return Py_None;
}  
  
static PyObject *PAvgConfig(PyObject *self, PyObject *args) {
  char *s;
  int ns, *n, *kappa;
  double *nq;
  
  if (sfac_file) {
    SFACStatement("AvgConfig", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &s)) return NULL;

  ns = GetAverageConfigFromString(&n, &kappa, &nq, s);
  if (ns <= 0) return NULL;

  if (SetAverageConfig(ns, n, kappa, nq) < 0) return NULL;

  free(n);
  free(kappa);
  free(nq);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetAvgConfig(PyObject *self, PyObject *args) {
  PyObject *acfg, *shell;
  int ns, i, m, kl, j;
  int *n, *kappa;
  double *nq, a;  

  if (sfac_file) {
    SFACStatement("SetAvgConfig", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "O", &acfg))
    return NULL;

  ns = PyList_Size(acfg);
  if (ns <= 0) return NULL;
  
  n = malloc(sizeof(int)*ns);
  kappa = malloc(sizeof(int)*ns);
  nq = malloc(sizeof(double)*ns);
  
  for (i = 0; i < ns; i++) {
    shell = PyList_GetItem(acfg, i);
    PyArg_ParseTuple(shell, "iiid", &m, &kl, &j, &a);
    n[i] = m;
    if (j < 0) kappa[i] = kl;
    else kappa[i] = -(kl+1);
    nq[i] = a;
  }

  if (SetAverageConfig(ns, n, kappa, nq) < 0) return NULL;
  free(n);
  free(kappa);
  free(nq);

  Py_INCREF(Py_None);
  return Py_None;
}    
   
/** add a configuration to the list **/
static PyObject *PAddConfig(PyObject *self, PyObject *args) {
  CONFIG *cfg = NULL;
  char *group_name;
  int k;
  PyObject *python_cfg;

  if (sfac_file) {
    SFACStatement("AddConfig", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "sO", &group_name, &python_cfg)) {
    goto ERROR;
  }
  if (ConfigPythonToC(python_cfg, &cfg) < 0) goto ERROR;

  if (Couple(cfg) < 0) goto ERROR;

  k = GroupIndex(group_name);
  if (k < 0) goto ERROR;

  if (AddConfigToList(k, cfg) < 0) goto ERROR;
  free(cfg);
  
  Py_INCREF(Py_None);
  return Py_None;

 ERROR:
  if (cfg != NULL) free(cfg);
  return NULL;
}
 

static PyObject *PSetRadialGrid(PyObject *self, PyObject *args) {
  double rmax, rmin;

  if (sfac_file) {
    SFACStatement("SetRadialGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  rmin = -1.0;
  rmax = -1.0;
  if (!PyArg_ParseTuple(args, "|dd", &rmin, &rmax))
    return NULL;
  SetRadialGrid(rmin, rmax);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetTransitionCut(PyObject *self, PyObject *args) {
  double c;

  if (sfac_file) {
    SFACStatement("SetTransitionCut", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &c))
    return NULL;
  SetTransitionCut(c);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetAICut(PyObject *self, PyObject *args) {
  double c;

  if (sfac_file) {
    SFACStatement("SetAICut", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &c))
    return NULL;
  SetAICut(c);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetMaxRank(PyObject *self, PyObject *args) {
  int k;

  if (sfac_file) {
    SFACStatement("SetMaxRank", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &k))
    return NULL;
  SetMaxRank(2*k);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetAngZOptions(PyObject *self, PyObject *args) {
  int n;
  double c;
  double mc;

  if (sfac_file) {
    SFACStatement("SetAngZOptions", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  c = EPS3;
  mc = EPS3;
  if (!PyArg_ParseTuple(args, "i|dd", &n, &mc, &c))
    return NULL;
  SetAngZOptions(n, mc, c);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetAngZCut(PyObject *self, PyObject *args) {
  double c;

  if (sfac_file) {
    SFACStatement("SetAngZCut", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &c))
    return NULL;
  SetAngZCut(c);
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PSetMixCut(PyObject *self, PyObject *args) {
  double c;

  if (sfac_file) {
    SFACStatement("SetMixCut", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &c))
    return NULL;
  SetMixCut(c);
  Py_INCREF(Py_None);
  return Py_None;
}

/** coeff. of fractional parentage **/
static PyObject *PGetCFPOld(PyObject *self, PyObject *args) {
  int j2, q, dj, dw, pj, pw;
  double coeff = 0.0;
  
  if (sfac_file) {
    SFACStatement("GetCFPOld", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "iiiiii", &j2, &q, &dj, &dw, &pj, &pw))
    return NULL;
  if (CFP(&coeff, j2, q, dj, dw, pj, pw) == -1)
    return NULL;
  return Py_BuildValue("d", coeff);
}

/** 3j symbol **/
static PyObject *PGetW3j(PyObject *self, PyObject *args) {
  int j1, j2, j3, m1, m2, m3;

  if (sfac_file) {
    SFACStatement("GetW3j", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "iiiiii", &j1, &j2, &j3, &m1, &m2, &m3))
    return NULL;
  return Py_BuildValue("d", W3j(j1, j2, j3, m1, m2, m3));
}

/** 6j symbol **/
static PyObject *PGetW6j(PyObject *self, PyObject *args) {
  int j1, j2, j3, i1, i2, i3;

  if (sfac_file) {
    SFACStatement("GetW6j", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "iiiiii", &j1, &j2, &j3, &i1, &i2, &i3))
    return NULL;
  return Py_BuildValue("d", W6j(j1, j2, j3, i1, i2, i3));
}

/** 9j symbol **/
static PyObject *PGetW9j(PyObject *self, PyObject *args) {
  int j1, j2, j3, i1, i2, i3, k1, k2, k3;

  if (sfac_file) {
    SFACStatement("GetW9j", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "iiiiiiiii", 
			&j1, &j2, &j3, &i1, &i2, &i3, &k1, &k2, &k3))
    return NULL;
  return Py_BuildValue("d", W9j(j1, j2, j3, i1, i2, i3, k1, k2, k3));
}

/** clebsch gordan coeff. **/
static PyObject *PGetCG(PyObject *self, PyObject *args) {
  int j1, j2, j3, m1, m2, m3;

  if (sfac_file) {
    SFACStatement("GetCG", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "iiiiii", &j1, &j2, &j3, &m1, &m2, &m3))
    return NULL;
  return Py_BuildValue("d", ClebschGordan(j1, j2, j3, m1, m2, m3));
}

static PyObject *PSetAtom(PyObject *self, PyObject *args) {
  char *s;
  double z, mass;

  if (sfac_file) {
    SFACStatement("SetAtom", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  mass = 0.0;
  z = 0.0;
  if (!PyArg_ParseTuple(args, "s|dd", &s, &z, &mass)) return NULL;
  if (SetAtom(s, z, mass) < 0) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetHydrogenicNL(PyObject *self, PyObject *args) {
  int n, k;
  
  if (sfac_file) {
    SFACStatement("SetHydrogenicNL", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = -1;
  k = -1;
  if (!PyArg_ParseTuple(args, "|ii", &n, &k)) return NULL;
  if (k < 0) k = n-1;

  SetHydrogenicNL(n, k);
  Py_INCREF(Py_None);
  return Py_None;
}  

static int DecodeGroupArgs(PyObject *args, int **kg) {
  PyObject *p;
  char *s;
  int i, k, ng;  

  if (args) {
    if (!PyList_Check(args) && !PyTuple_Check(args)) return -1;
    ng = PySequence_Length(args);
  } else {
    ng = 0;
  }
  if (ng > 0) {
    p = PySequence_GetItem(args, 0);
    if (PyList_Check(p) || PyTuple_Check(p)) {
      if (ng > 1) {
	onError("there should only be one list or tuple");
	return -1;
      }
      ng = PySequence_Length(p);
      args = p;
      Py_DECREF(p);
    }
    (*kg) = malloc(sizeof(int)*ng);
    if (!(*kg)) {
      onError("not enough memory");
      return -1;
    }
    for (i = 0; i < ng; i++) {
      p = PySequence_GetItem(args, i);
      if (!PyString_Check(p)) {
	free((*kg));
	onError("argument must be a group name");
	return -1;
      }
      s = PyString_AsString(p);
      k = GroupExists(s);
      Py_DECREF(p);
      
      if (k < 0) {
	free((*kg));
	onError("group does not exist");
	return -1;
      }
      (*kg)[i] = k;
    }
  } else {
    ng = GetNumGroups();
    (*kg) = malloc(sizeof(int)*ng);
    if (!(*kg)) {
      onError("not enough memory");
      return -1;
    }
    for (i = 0; i < ng; i++) (*kg)[i] = i;
  }
  return ng;
}

static PyObject *POptimizeRadial(PyObject *self, PyObject *args) {
  int ng, i, k;
  int *kg;
  double z;
  PyObject *p;
  double *weight;

  if (sfac_file) {
    SFACStatement("OptimizeRadial", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  ng = PyTuple_Size(args);
  if (ng == 0) {
    ng = 0;
    kg = NULL;
    weight = NULL;
    goto END;
  } 

  p = PyTuple_GET_ITEM(args, 0);
  if (PyString_Check(p)) {
    weight = NULL;
    ng = DecodeGroupArgs(args, &kg);
    if (ng < 0) {
      onError("the group argument format error");
      return NULL;
    }
  } else {
    ng = DecodeGroupArgs(p, &kg);
    if (ng < 0) {
      onError("the groups must be in a sequence");
      return NULL;
    }

    if (PyTuple_Size(args) == 1) {
      weight = NULL;
    } else {
      args = PyTuple_GET_ITEM(args, 1);
      k = PySequence_Length(args);
      if (k < 0 || k > ng) {
	onError("weights must be a sequence");
	return NULL;
      }
      weight = malloc(sizeof(double)*ng);
      z = 0.0;
      for (i = 0; i < k; i++) {
	p = PySequence_GetItem(args, i);
	if (!PyFloat_Check(p)) {
	  onError("weights must be float numbers");
	  if (weight) free(weight);
	  return NULL;
	}
	weight[i] = PyFloat_AsDouble(p);
	Py_DECREF(p);
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
    onError("error occured in OptimizeRadial");
    return NULL;
  }
  if (weight) free(weight);
  if (kg) free(kg);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PGetPotential(PyObject *self, PyObject *args) {
  char *s;

  if (sfac_file) {
    SFACStatement("GetPotential", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
  GetPotential(s);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSolveBound(PyObject *self, PyObject *args) {
  int n, kappa;
  ORBITAL *orb;
  int k;

  if (sfac_file) {
    SFACStatement("SolveBound", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "ii", &n, &kappa)) return NULL;
  if (n <= 0) {
    onError("n must be greater than 0 for bound states");
    return NULL;
  }
  k = OrbitalIndex(n, kappa, 0.0);
  if (k < 0) {
    onError("fatal error in solving dirac equation");
    return NULL;
  }
  orb = GetOrbital(k);
  return Py_BuildValue("d", orb->energy);
}

static PyObject *PStructure(PyObject *self, PyObject *args) {
  int i, k, ng, ns, nlevels;
  int ngp;
  int *kg, *kgp;
  char *fn;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("Structure", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  p = NULL;
  q = NULL;
  ngp = 0;
  kgp = NULL;
  
  if (!(PyArg_ParseTuple(args, "s|OO", &fn, &p, &q))) return NULL;
  
  if (p) {
    if (PyTuple_Check(p) || PyList_Check(p)) {
      ng = DecodeGroupArgs(p, &kg);
      if (ng < 0) return NULL;
      if (q) {
	if (!PyTuple_Check(q) && !PyList_Check(q)) return NULL;
	if (PySequence_Length(q) > 0) {
	  ngp = DecodeGroupArgs(q, &kgp);
	}
      }
    } else {
      return NULL;
    }
  } else {
    ng = DecodeGroupArgs(NULL, &kg);  
    if (ng < 0) return NULL;
  }

  if (ngp < 0) return NULL;
  
  nlevels = GetNumLevels();
  ns = MAX_SYMMETRIES;
  for (i = 0; i < ns; i++) {
    k = ConstructHamilton(i, ng, kg, ngp, kgp);
    if (k < 0) continue;
    if (DiagnolizeHamilton() < 0) return NULL;
    AddToLevels();
  }
  SortLevels(nlevels, -1);
  SaveLevels(fn, nlevels, -1);
  if (ng > 0) free(kg);
  if (ngp > 0) free(kgp);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PClearOrbitalTable(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("ClearOrbitalTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = 1;
  if (!PyArg_ParseTuple(args, "|i", &m)) return NULL;
  
  ClearOrbitalTable(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PClearLevelTable(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("ClearLevelTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  ClearLevelTable();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSortLevels(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("SortLevels", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  SortLevels(0, 0);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTestAngular(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("TestAngular", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  /* do nothing if the debug flag is not set in compilation */
#if FAC_DEBUG  
  TestAngular();
#else 
  printf("Turn on the FAC_DEBUG flag in compilation\n");
#endif /* FAC_DEBUG */

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetTransitionOptions(PyObject *self, PyObject *args) {
  int gauge, mode, max_m, max_e;

  if (sfac_file) {
    SFACStatement("SetTransitionOptions", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  max_e = 4;
  max_m = 4;
  if (!PyArg_ParseTuple(args, "ii|ii", &gauge, &mode, &max_e, &max_m)) 
    return NULL;
  SetTransitionOptions(gauge, mode, max_e, max_m);
  Py_INCREF(Py_None);
  return Py_None;
} 

static int SelectLevels(PyObject *p, int **t) {
  int n, ng, *kg, i, j, k, im, m, m0;
  int nrg, *krg, nrec;
  PyObject *q;
  int ig, nlevels;
  LEVEL *lev;
  SYMMETRY *sym;
  STATE *s;
  char rgn[GROUP_NAME_LEN];

  if (!PyList_Check(p) && !PyTuple_Check(p)) return 0;
  n = PySequence_Length(p);
  if (n > 0) {
    q = PySequence_GetItem(p, 0);
    if (PyString_Check(q)) {
      ng = DecodeGroupArgs(p, &kg);
      if (ng <= 0) {
	return 0;
      }
      nlevels = GetNumLevels();
      (*t) = malloc(sizeof(int)*nlevels);
      if (!(*t)) return 0;
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
      Py_DECREF(q);
      (*t) = realloc(*t, k*sizeof(int));
      return k;
    } else if (PyList_Check(q)) {
      if (n != 2) {
	printf("recombined states specification unrecoganized\n");
	return -1;
      }
      ng = DecodeGroupArgs(q, &kg);
      if (ng <= 0) return -1;
      Py_DECREF(q);
      q = PySequence_GetItem(p, 1);
      if (PyList_Check(q)) {
	p = q;
	m0 = 0;
	n = PySequence_Length(q);
      } else if (PyInt_Check(q)) {
	m0 = 1;
      } else {
	printf("Level specification unrecoganized\n");
	return -1;
      }
      nrg = ng;
      krg = malloc(sizeof(int)*nrg);
      nlevels = GetNumLevels();
      (*t) = malloc(sizeof(int)*nlevels);
      if (!(*t)) return 0;
      k = 0;
      Py_DECREF(q);
      for (m = m0; m < n; m++) {
	q = PySequence_GetItem(p, m);
	nrec = PyInt_AS_LONG(q);
	Py_DECREF(q);
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
      return k;
    } else {
      (*t) = malloc(sizeof(int)*n);
      if (!(*t)) return 0;      
      for (i = 0; i < n; i++) {
	q = PySequence_GetItem(p, i);
	if (!PyInt_Check(q)) {
	  free(*t);
	  return 0;
	}
	Py_DECREF(q);
	(*t)[i] = PyInt_AS_LONG(q);
      }
      return n;
    }
  }
  return 0;
}

static PyObject *PTransitionTable(PyObject *self, PyObject *args) {
  char *s;
  int n, m;
  int nlow, nup, *low, *up;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("TransitionTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;
  m = -1;

  n = PyTuple_Size(args); 
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
    SaveTransition(nlow, low, nup, up, s, m);
  } else if (n == 2) {
    if (!PyArg_ParseTuple(args, "si", &s, &m)) return NULL;
    SaveTransition(nlow, low, nup, up, s, m);
  } else if (n == 3) {
    if (!PyArg_ParseTuple(args, "sOO", &s, &p, &q)) return NULL;
    nlow = SelectLevels(p, &low);
    if (nlow <= 0) return NULL;
    nup = SelectLevels(q, &up);
    if (nup <= 0) return NULL;
    SaveTransition(nlow, low, nup, up, s, m);
    free(low);
    free(up);
  } else if (n == 4) {
    if (!PyArg_ParseTuple(args, "sOOi", &s, &p, &q, &m)) {
      return NULL;
    }
    nlow = SelectLevels(p, &low);
    if (nlow <= 0) return NULL;
    nup = SelectLevels(q, &up);
    if (nup <= 0) return NULL;
    SaveTransition(nlow, low, nup, up, s, m);
    free(low);
    free(up);
  } else {
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PBasisTable(PyObject *self, PyObject *args) {
  char *s;
  if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
  GetBasisTable(s);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PCETable(PyObject *self, PyObject *args) {
  char *s;
  int n;
  int nlow, nup, *low, *up;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("CETable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;

  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
    SaveExcitation(nlow, low, nup, up, 0, s);
  } else if (n == 2) {
    if (!PyArg_ParseTuple(args, "sO", &s, &p)) return NULL;
    nlow = SelectLevels(p, &low);
    if (nlow <= 0) return NULL;
    SaveExcitation(nlow, low, nlow, low, 0, s);
    free(low);
  } else if (n == 3) {
    if (!PyArg_ParseTuple(args, "sOO", &s, &p, &q)) return NULL;
    nlow = SelectLevels(p, &low);
    if (nlow <= 0) return NULL;
    nup = SelectLevels(q, &up);
    if (nup <= 0) return NULL;
    SaveExcitation(nlow, low, nup, up, 0, s);
    free(low);
    free(up);
  } else {
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PCETableMSub(PyObject *self, PyObject *args) {
  char *s;
  int n;
  int nlow, nup, *low, *up;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("CETableMSub", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;

  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
    SaveExcitation(nlow, low, nup, up, 1, s);
  } else if (n == 2) {
    if (!PyArg_ParseTuple(args, "sO", &s, &p)) return NULL;
    nlow = SelectLevels(p, &low);
    if (nlow <= 0) return NULL;
    SaveExcitation(nlow, low, nlow, low, 1, s);
    free(low);
  } else if (n == 3) {
    if (!PyArg_ParseTuple(args, "sOO", &s, &p, &q)) return NULL;
    nlow = SelectLevels(p, &low);
    if (nlow <= 0) return NULL;
    nup = SelectLevels(q, &up);
    if (nup <= 0) return NULL;
    SaveExcitation(nlow, low, nup, up, 1, s);
    free(low);
    free(up);
  } else {
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSpline(PyObject *self, PyObject *args) {
  PyObject *px, *py, *py2;
  double *x, *y, *y2, dy1, dy2;
  int n, i;

  if (sfac_file) {
    SFACStatement("Spline", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  dy1 = 1E30;
  dy2 = 1E30;
  if (!PyArg_ParseTuple(args, "OO|dd", &px, &py, &dy1, &dy2)) return NULL;
  if (!PyList_Check(px) || !PyList_Check(py)) return NULL;
  n = PyList_Size(px);
  if (PyList_Size(py) != n) return NULL;
  if (n == 0) return NULL;

  x = malloc(sizeof(double)*n);
  y = malloc(sizeof(double)*n);
  y2 = malloc(sizeof(double)*n);
  
  for (i = 0; i < n; i++) {
    x[i] = PyFloat_AsDouble(PyList_GetItem(px, i));
    y[i] = PyFloat_AsDouble(PyList_GetItem(py, i));
  }

  spline(x, y, n, dy1, dy2, y2);
  py2 = Py_BuildValue("[]");
  for (i = 0; i < n; i++) {
    PyList_Append(py2, Py_BuildValue("d", y2[i]));
  }
  free(x);
  free(y);
  free(y2);
  
  return py2;
}

static PyObject *PSplint(PyObject *self, PyObject *args) {  
  PyObject *px, *py, *py2;
  double *x, *y, *y2, x0, y0;
  int n, i;  

  if (sfac_file) {
    SFACStatement("Splint", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "OOOd", &px, &py, &py2, &x0)) return NULL;
  if (!PyList_Check(px) || !PyList_Check(py)) return NULL;
  n = PyList_Size(px);
  if (PyList_Size(py) != n) return NULL;
  if (PyList_Size(py2) != n) return NULL;
  x = malloc(sizeof(double)*n);
  y = malloc(sizeof(double)*n);
  y2 = malloc(sizeof(double)*n);  
  for (i = 0; i < n; i++) {
    x[i] = PyFloat_AsDouble(PyList_GetItem(px, i));
    y[i] = PyFloat_AsDouble(PyList_GetItem(py, i));
    y2[i] = PyFloat_AsDouble(PyList_GetItem(py2, i));
  }

  splint(x, y, y2, n, x0, &y0);

  free(x);
  free(y);
  free(y2);
  
  return Py_BuildValue("d", y0);
}  
  
static PyObject *PTestSpline(PyObject *self, PyObject *args) {
#define M 10
#define N 10
  int i, j;
  double f, ff, x1x2, xx1, xx2, x1[M], x2[N], dy[M][N], dy2[M][N];
  double *y[M], *y2[M];

  if (sfac_file) {
    SFACStatement("TestSpline", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  for (i = 0; i < M; i++) {
    x1[i] = 0.2*i;
    y[i] = dy[i];
    y2[i] = dy2[i];
  }
  for (i = 0; i < N; i++) x2[i] = 0.2*i;
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      x1x2 = x1[i]*x2[j];
      y[i][j] = x1x2*exp(-x1x2);
    }
  }
  splie2(x1, x2, y, M, N, y2);
  printf("%9s %12s %14s %12s\n","x1","x2","splin2","actual");
  for (i=0; i < 10; i++) {
    xx1=0.1*i;
    xx2=xx1*xx1;
    splin2(x1,x2,y,y2,M,N,xx1,xx2,&f);
    x1x2=xx1*xx2;
    ff=x1x2*exp(-x1x2);
    printf("%12.6f %12.6f %12.6f %12.6f\n",xx1,xx2,f,ff);
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetTEGrid(PyObject *self, PyObject *args) {
  int i, n;
  double xg[MAXNTE];
  int ng, err;
  PyObject *p, *pi;
  double emin, emax;
 
  if (sfac_file) {
    SFACStatement("SetTEGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (PyInt_Check(p)) {
      ng = PyInt_AsLong(p);
      err = SetCETEGrid(ng, -1.0, 0.0);
    } else if (!PyList_Check(p) && !PyTuple_Check(p)) {
      return NULL;
    } else {
      ng = PySequence_Length(p);      
      for (i = 0; i < ng; i++) {
	pi = PySequence_GetItem(p, i);
	xg[i] = PyFloat_AsDouble(pi)/HARTREE_EV;
	Py_DECREF(pi);
      }
      err = SetCETEGridDetail(ng, xg);
    }
  } else if (n == 3) {
    if (!PyArg_ParseTuple(args, "idd", &ng, &emin, &emax)) 
      return NULL;
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    err = SetCETEGrid(ng, emin, emax);
  } else {
    return NULL;
  }

  if (err < 0) return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}
  
static  PyObject *PSetCEQkMode(PyObject *self, PyObject *args) {
  PyObject *p;
  int m;
  double tol;

  m = QK_DEFAULT;
  tol = -1;
  if (!PyArg_ParseTuple(args, "|Od", &p, &tol)) return NULL;
  if (PyString_Check(p)) {
    p = PyDict_GetItem(QKMODE, p);
  } 
  if (PyInt_Check(p)) {
    m = PyInt_AsLong(p);
  } else {
    return NULL;
  }
  
  if (m >= QK_CB) {
    printf("CEQkMode must < %d\n", QK_CB);
    return NULL;
  }
  SetCEQkMode(m, tol);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCEFormat(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetCEFormat", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCEFormat(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCEGridType(PyObject *self, PyObject *args) {
  int type;
  
  if (sfac_file) {
    SFACStatement("SetCEGridType", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &type)) return NULL;
  SetCEEGridType(type);
  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PSetCEGridLimits(PyObject *self, PyObject *args) {
  double emin, emax;
  int type;

  if (sfac_file) {
    SFACStatement("SetCEGridLimits", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  emin = -1;
  emax = -1;
  type = 0;
  if (!PyArg_ParseTuple(args, "|ddi", &emin, &emax, &type)) return NULL;
  SetCEEGridLimits(emin, emax, type);
  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PSetCEGrid(PyObject *self, PyObject *args) {
  int n;
  double xg[MAXNE];
  int ng;
  double emin, emax, eth;
  PyObject *p, *pi;
  int i, err;
  
  if (sfac_file) {
    SFACStatement("SetCEGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  eth = 0.0;
  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (PyInt_Check(p)) {
      ng = PyInt_AsLong(p);
      err = SetCEEGrid(ng, -1.0, -1.0, 0.0);
    } else if (PyList_Check(p) || PyTuple_Check(p)) {
      ng = PySequence_Length(p);      
      for (i = 0; i < ng; i++) {
	pi = PySequence_GetItem(p, i);
	xg[i] = PyFloat_AsDouble(pi)/HARTREE_EV;
	Py_DECREF(pi);
      }
      err = SetCEEGridDetail(ng, xg);
    } else {
      return NULL;
    }
  } else if (n == 3 || n == 4) {
    if (!PyArg_ParseTuple(args, "idd|d", &ng, &emin, &emax, &eth)) 
      return NULL;
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetCEEGrid(ng, emin, emax, eth);
  } else {
    return NULL;
  }

  if (err < 0) return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PSetUsrCEGridType(PyObject *self, PyObject *args) {
  int type;
  
  if (sfac_file) {
    SFACStatement("SetUsrCEGridType", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &type)) return NULL;
  SetUsrCEEGridType(type);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetUsrCEGrid(PyObject *self, PyObject *args) {
  int n;
  double xg[MAXNUSR];
  int ng;
  double emin, emax, eth;
  PyObject *p, *pi;
  int i, err;

  if (sfac_file) {
    SFACStatement("SetUsrCEGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  eth = 0.0;
  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (PyInt_Check(p)) {
      ng = PyInt_AsLong(p);
      err = SetUsrCEEGrid(ng, -1.0, -1.0, 0.0);
    } else if (PyList_Check(p) || PyTuple_Check(p)) {
      ng = PySequence_Length(p);      
      for (i = 0; i < ng; i++) {
	pi = PySequence_GetItem(p, i);
	xg[i] = PyFloat_AsDouble(pi)/HARTREE_EV;
	Py_DECREF(pi);
      }
      if (eth > 0) eth /= HARTREE_EV;
      err = SetUsrCEEGridDetail(ng, xg);
    } else {
      return NULL;
    }
  } else if (n == 3 || n == 4) {
    if (!PyArg_ParseTuple(args, "idd|d", &ng, &emin, &emax, &eth)) 
      return NULL;
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetUsrCEEGrid(ng, emin, emax, eth);
  } else {
    return NULL;
  }

  if (err < 0) return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PSetCEPWGridType(PyObject *self, PyObject *args) {
  int type;
  
  if (sfac_file) {
    SFACStatement("SetCEPWGridType", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &type)) return NULL;
  SetCEPWGridType(type);
  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetCEPWOptions(PyObject *self, PyObject *args) {
  int qr, max, kl_cb;
  double tol;

  qr = 0;
  max = 500;
  kl_cb = 100;
  tol = 5E-2;

  if (!PyArg_ParseTuple(args, "d|iii", 
			&tol, &max, &qr, &kl_cb)) return NULL;
  SetCEPWOptions(qr, max, kl_cb, tol);
  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetCEPWGrid(PyObject *self, PyObject *args) {
  int ns, i;
  int n;
  PyObject *p, *q;
  int *m, *step;

  n = PyTuple_Size(args);
  if (n == 1) {    
    if (!PyArg_ParseTuple(args, "i", &ns)) return NULL;
    SetCEPWGrid(-ns, NULL, NULL);
  } else {
    if (!PyArg_ParseTuple(args, "OO", &p, &q)) return NULL;
    if (!PyList_Check(p) || !PyList_Check(q)) return NULL;
    ns = PyList_Size(p);
    if (ns != PyList_Size(q)) return NULL;
    if (ns == 0) return NULL;
    m = (int *) malloc(ns*sizeof(int));
    step = (int *) malloc(ns*sizeof(int));
    for (i = 0; i < ns; i++) {
      m[i] = PyInt_AsLong(PyList_GetItem(p, i));
      step[i] = PyInt_AsLong(PyList_GetItem(q, i));
    }
    SetCEPWGrid(ns, m, step);
    free(m);
    free(step);
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PWaveFuncTable(PyObject *self, PyObject *args) {
  char *s;
  int k, n;
  double e;
  
  if (sfac_file) {
    SFACStatement("WaveFuncTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  e = 0.0;
  if (!PyArg_ParseTuple(args, "sii|d", &s, &n, &k, &e)) return NULL;
  WaveFuncTable(s, n, k, e);

  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetRecQkMode(PyObject *self, PyObject *args) {
  PyObject *p;
  int m;
  double tol;

  m = QK_DEFAULT;
  tol = -1;
  if (!PyArg_ParseTuple(args, "|Od", &p, &tol)) return NULL;
  if (PyString_Check(p)) {
    p = PyDict_GetItem(QKMODE, p);
  }
  if (PyInt_Check(p)) {
    m = PyInt_AsLong(p);
  } else {
    return NULL;
  }

  if (m >= QK_CB) {
    printf("RecQkMode must < %d\n", QK_CB);
    return NULL;
  }
  SetRecQkMode(m, tol);
  Py_INCREF(Py_None);
  return Py_None;
}
  
static  PyObject *PSetRecPWOptions(PyObject *self, PyObject *args) {
  int kl_interp, max_kl;
  
  max_kl = -1;
  if (!PyArg_ParseTuple(args, "i|i", 
			&kl_interp, &max_kl)) return NULL;
  if (max_kl < 0) max_kl = kl_interp;
  SetRecPWOptions(kl_interp, max_kl);
  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetRecPWLimits(PyObject *self, PyObject *args) {
  int m1, m2;

  if (!PyArg_ParseTuple(args, "ii", &m1, &m2)) return NULL;
  SetRecPWLimits(m1, m2);
  Py_INCREF(Py_None);
  return Py_None;
}
  

static  PyObject *PSetRecSpectator(PyObject *self, PyObject *args) {
  int n_spec, n_frozen, n_max;

  n_spec = 0;
  n_frozen = 0;
  n_max = 0;
  if (!PyArg_ParseTuple(args, "i|ii", 
			&n_spec, &n_frozen, &n_max)) return NULL;
  if (n_frozen == 0) n_frozen = n_spec;
  if (n_max == 0) n_max = 100;

  SetRecSpectator(n_max, n_frozen, n_spec);
  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PRecStates(PyObject *self, PyObject *args) { 
  int ng;
  int *kg;
  int n;
  char *fn;
  PyObject *gargs;

  if (sfac_file) {
    SFACStatement("RecStates", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "sOi", &fn, &gargs, &n)) return NULL;
  ng = DecodeGroupArgs(gargs, &kg);
  RecStates(n, ng, kg, fn);
  
  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PRRTable(PyObject *self, PyObject *args) { 
  int nlow, *low, nup, *up;
  int m;
  char *s;
  PyObject *p, *q;
  
  if (sfac_file) {
    SFACStatement("RRTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = -1;
  if (!PyArg_ParseTuple(args, "sOO|i", &s, &p, &q, &m)) {
    printf("Unrecognized parameters in RRTable\n");
    return NULL;
  }
  nlow = SelectLevels(p, &low);
  nup = SelectLevels(q, &up);
  SaveRecRR(nlow, low, nup, up, s, m);
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetUsrPEGridType(PyObject *self, PyObject *args) {
  int type;
  
  if (sfac_file) {
    SFACStatement("SetUsrPEGridType", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &type)) return NULL;
  SetUsrPEGridType(type);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetUsrPEGrid(PyObject *self, PyObject *args) {
  int n;
  double xg[MAXNUSR];
  int ng;
  double emin, emax, eth;
  PyObject *p, *pi;
  int i, err;

  if (sfac_file) {
    SFACStatement("SetUsrPEGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  eth = 0.0;
  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (PyInt_Check(p)) {
      ng = PyInt_AsLong(p);
      err = SetUsrPEGrid(ng, -1.0, -1.0, 0.0);
    } else if (PyList_Check(p) || PyTuple_Check(p)) {
      ng = PySequence_Length(p);      
      for (i = 0; i < ng; i++) {
	pi = PySequence_GetItem(p, i);
	xg[i] = PyFloat_AsDouble(pi)/HARTREE_EV;
	Py_DECREF(pi);
      }
      err = SetUsrPEGridDetail(ng, xg);
    } else {
      return NULL;
    }
  } else if (n == 3 || n == 4) {
    if (!PyArg_ParseTuple(args, "idd|d", &ng, &emin, &emax, &eth)) 
      return NULL;
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetUsrPEGrid(ng, emin, emax, eth);
  } else {
    return NULL;
  }

  if (err < 0) return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PSetRRTEGrid(PyObject *self, PyObject *args) {
  int i, n;
  double xg[MAXNTE];
  int ng, err;
  PyObject *p, *pi;
  double emin, emax;
 
  if (sfac_file) {
    SFACStatement("SetRRTEGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (PyInt_Check(p)) {
      ng = PyInt_AsLong(p);
      err = SetRRTEGrid(ng, -1.0, 0.0);
    } else if (!PyList_Check(p) && !PyTuple_Check(p)) {
      return NULL;
    } else {
      ng = PySequence_Length(p);      
      for (i = 0; i < ng; i++) {
	pi = PySequence_GetItem(p, i);
	xg[i] = PyFloat_AsDouble(pi)/HARTREE_EV;
	Py_DECREF(pi);
      }
      err = SetRRTEGridDetail(ng, xg);
    }
  } else if (n == 3) {
    if (!PyArg_ParseTuple(args, "idd", &ng, &emin, &emax)) 
      return NULL;
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    err = SetRRTEGrid(ng, emin, emax);
  } else {
    return NULL;
  }

  if (err < 0) return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetPEGridLimits(PyObject *self, PyObject *args) {
  double emin, emax;
  int type;

  if (sfac_file) {
    SFACStatement("SetPEGridLimits", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  emin = -1;
  emax = -1;
  type = 0;
  if (!PyArg_ParseTuple(args, "|ddi", &emin, &emax, &type)) return NULL;
  SetPEGridLimits(emin, emax, type);
  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PSetPEGrid(PyObject *self, PyObject *args) {  
  int n;
  double xg[MAXNE];
  int ng;
  double emin, emax, eth;
  PyObject *p, *pi;
  int i, err;
  
  if (sfac_file) {
    SFACStatement("SetPEGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  eth = 0.0;
  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (PyInt_Check(p)) {
      ng = PyInt_AsLong(p);
      err = SetPEGrid(ng, -1.0, -1.0, 0.0);
    } else if (PyList_Check(p) || PyTuple_Check(p)) {
      ng = PySequence_Length(p);      
      for (i = 0; i < ng; i++) {
	pi = PySequence_GetItem(p, i);
	xg[i] = PyFloat_AsDouble(pi)/HARTREE_EV;
	Py_DECREF(pi);
      }
      err = SetPEGridDetail(ng, xg);
    } else {
      return NULL;
    }
  } else if (n == 3 || n == 4) {
    if (!PyArg_ParseTuple(args, "idd|d", &ng, &emin, &emax, &eth)) 
      return NULL;
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetPEGrid(ng, emin, emax, eth);
  } else {
    return NULL;
  }

  if (err < 0) return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PAITable(PyObject *self, PyObject *args) { 
  int nlow, *low, nup, *up, c;
  char *s;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("AITable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  c = 0;
  if (!PyArg_ParseTuple(args, "sOO|i", &s, &p, &q, &c)) return NULL;
  nlow = SelectLevels(p, &low);
  nup = SelectLevels(q, &up);
  SaveAI(nlow, low, nup, up, s, c);
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PDRTable(PyObject *self, PyObject *args) { 
  int nf, *f, na, *a, nb, *b, ng, *g, c;
  char *s1, *s2;
  PyObject *p, *q, *r, *t;

  if (sfac_file) {
    SFACStatement("DRTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "ssOOOO|i", &s1, &s2, &p, &q, &r, &t, &c)) 
    return NULL;
  nf = SelectLevels(p, &f);
  na = SelectLevels(q, &a);
  nb = SelectLevels(r, &b);
  ng = SelectLevels(t, &g);
  
  SaveDR(nf, f, na, a, nb, b, ng, g, s1, s2, c);

  if (nf > 0) free(f);
  if (na > 0) free(a);
  if (nb > 0) free(b);
  if (ng > 0) free(g);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTestMyArray(PyObject *self, PyObject *args) {
  ARRAY a;
  double d;
  double *b;
  MULTI ma;
  int k[3] = {101, 2550, 333};
  int block[3] = {10, 20, 50};
  int i, j, m;

  if (sfac_file) {
    SFACStatement("TestMyArray", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

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
  MultiFreeData(ma.array, ma.ndim, NULL);

  printf("> ");
  scanf("%d", &i); 
 
  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PSaveOrbitals(PyObject *self, PyObject *args) {  
  int n, i, norbs;
  double e;
  PyObject *p;

  if (sfac_file) {
    SFACStatement("SaveOrbitals", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  norbs = GetNumOrbitals();
  if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
  if (PyInt_Check(p)) {
    n = PyInt_AsLong(p);
    if (n <= 0) {
      if (SaveAllContinua(1) < 0) return NULL;
    } else {
      for (i = 0; i < norbs; i++) {
	if (GetOrbital(i)->n == n) {
	  if (SaveOrbital(i) < 0) return NULL;
	  FreeOrbital(i);
	}
      }
    }
  } else if (PyFloat_Check(p)) {
    e = PyFloat_AsDouble(p)/HARTREE_EV;
    if (SaveContinua(e, 1) < 0) return NULL;
  } else {
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeOrbitals(PyObject *self, PyObject *args) {  
  int n, i, norbs;
  double e;
  PyObject *p;

  if (sfac_file) {
    SFACStatement("FreeOrbitals", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  norbs = GetNumOrbitals();
  if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
  if (PyInt_Check(p)) {
    n = PyInt_AsLong(p);
    if (n <= 0) {
      for (i = 0; i < norbs; i++) {
	if (GetOrbital(i)->n <= 0) {
	  FreeOrbital(i);
	}
      }
    } else {
      for (i = 0; i < norbs; i++) {
	if (GetOrbital(i)->n == n) {
	  FreeOrbital(i);
	}
      }
    }
  } else if (PyFloat_Check(p)) {
    e = PyFloat_AsDouble(p)/HARTREE_EV;
    if (FreeContinua(e) < 0) return NULL;
  } else {
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PDROpen(PyObject *self, PyObject *args) {
  int i, n, *nlev, *n0, nop;
  PyObject *p;

  if (sfac_file) {
    SFACStatement("DROpen", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
  n = SelectLevels(p, &nlev);
  if (n == 0) return Py_BuildValue("[]");

  nop = DROpen(n, nlev, &n0); 

  p = Py_BuildValue("[]");
  for (i = nop-1; i >= 0; i--) {
    PyList_Append(p, Py_BuildValue("i", n0[i]));
  }
  free(n0);
  free(nlev);
  
  return p;
}

static PyObject *PFreeResidual(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("FreeResidual", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  FreeResidualArray();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeSlater(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("FreeSlater", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  FreeSlaterArray();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeMultipole(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("FreeMultipole", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  FreeMultipoleArray();
  FreeMomentsArray();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeRecPk(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("FreeRecPk", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  FreeRecPk();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeRecQk(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("FreeRecQk", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  FreeRecQk();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeExcitationPk(PyObject *self, PyObject *args) {
  int ie;

  if (sfac_file) {
    SFACStatement("FreeExcitationPk", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  ie = -1;
  if (!PyArg_ParseTuple(args, "|i", &ie)) return NULL;
  FreeExcitationPk(ie);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeExcitationQk(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("FreeExcitationQk", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  FreeExcitationQk();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeRecAngZ(PyObject *self, PyObject *args) { 

  if (sfac_file) {
    SFACStatement("FreeRecAngZ", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  FreeRecAngZ();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *PFreeAngZ(PyObject *self, PyObject *args) { 
  PyObject *p;
  int i, m;
  int n, *kg;

  if (sfac_file) {
    SFACStatement("FreeAngZ", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  p = NULL;
  m = -1;
  if (!PyArg_ParseTuple(args, "|Oi", &p, &m)) return NULL;
  if (p == NULL) {
    FreeAngZ(-1, m);
  } else {
    n = DecodeGroupArgs(p, &kg);
    for (i = 0; i < n; i++) {
      FreeAngZ(kg[i], m);
    }
  }

  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *PSetIEGrid(PyObject *self, PyObject *args) {
  int i, n;
  double xg[MAXNTE];
  int ng, err;
  PyObject *p, *pi;
  double emin, emax;
 
  if (sfac_file) {
    SFACStatement("SetIEGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (PyInt_Check(p)) {
      ng = PyInt_AsLong(p);
      err = SetIEGrid(ng, -1.0, 0.0);
    } else if (!PyList_Check(p) && !PyTuple_Check(p)) {
      return NULL;
    } else {
      ng = PySequence_Length(p);      
      for (i = 0; i < ng; i++) {
	pi = PySequence_GetItem(p, i);
	xg[i] = PyFloat_AsDouble(pi)/HARTREE_EV;
	Py_DECREF(pi);
      }
      err = SetIEGridDetail(ng, xg);
    }
  } else if (n == 3) {
    if (!PyArg_ParseTuple(args, "idd", &ng, &emin, &emax)) 
      return NULL;
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    err = SetIEGrid(ng, emin, emax);
  } else {
    return NULL;
  }

  if (err < 0) return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetCIQkMode(PyObject *self, PyObject *args) {
  PyObject *p;
  int m;
  double tol;

  m = QK_DEFAULT;
  tol = -1.0;
  if (!PyArg_ParseTuple(args, "|Od", &p, &tol)) return NULL;
  if (PyString_Check(p)) {
    p = PyDict_GetItem(QKMODE, p);
  }
  if (PyInt_Check(p)) {
    m = PyInt_AsLong(p);
  } else {
    return NULL;
  }
  if (m < QK_CB) {
    printf("CIQkMode must >= %d\n", QK_CB);
    return NULL;
  }
  SetCIQkMode(m, tol);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCIFormat(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetCIFormat", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCIFormat(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCIEGridLimits(PyObject *self, PyObject *args) {
  double emin, emax;
  int type;

  if (sfac_file) {
    SFACStatement("SetCIEGridLimits", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  emin = -1;
  emax = -1;
  type = 0;
  if (!PyArg_ParseTuple(args, "|ddi", &emin, &emax, &type)) return NULL;
  SetCIEGridLimits(emin, emax, type);
  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PSetCIEGrid(PyObject *self, PyObject *args) {
  int n;
  double xg[MAXNE];
  int ng;
  double emin, emax, eth;
  PyObject *p, *pi;
  int i, err;
  
  if (sfac_file) {
    SFACStatement("SetCIEGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  eth = 0.0;
  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (PyInt_Check(p)) {
      ng = PyInt_AsLong(p);
      err = SetCIEGrid(ng, -1.0, -1.0, 0.0);
    } else if (PyList_Check(p) || PyTuple_Check(p)) {
      ng = PySequence_Length(p);      
      for (i = 0; i < ng; i++) {
	pi = PySequence_GetItem(p, i);
	xg[i] = PyFloat_AsDouble(pi)/HARTREE_EV;
	Py_DECREF(pi);
      }
      err = SetCIEGridDetail(ng, xg);
    } else {
      return NULL;
    }
  } else if (n == 3 || n == 4) {
    if (!PyArg_ParseTuple(args, "idd|d", &ng, &emin, &emax, &eth)) 
      return NULL;
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetCIEGrid(ng, emin, emax, eth);
  } else {
    return NULL;
  }

  if (err < 0) return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PSetUsrCIEGridType(PyObject *self, PyObject *args) {
  int type;
  
  if (sfac_file) {
    SFACStatement("SetUsrCIEGridType", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  type = -1;
  if (!PyArg_ParseTuple(args, "i", &type)) return NULL;
  SetUsrCIEGridType(type);
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PSetUsrCIEGrid(PyObject *self, PyObject *args) {
  int n;
  double xg[MAXNUSR];
  int ng;
  double emin, emax, eth;
  PyObject *p, *pi;
  int i, err;

  if (sfac_file) {
    SFACStatement("SetUsrCIEGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  eth = 0.0;
  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "O", &p)) return NULL;
    if (PyInt_Check(p)) {
      ng = PyInt_AsLong(p);
      err = SetUsrCIEGrid(ng, -1.0, -1.0, 0.0);
    } else if (PyList_Check(p) || PyTuple_Check(p)) {
      ng = PySequence_Length(p);      
      for (i = 0; i < ng; i++) {
	pi = PySequence_GetItem(p, i);
	xg[i] = PyFloat_AsDouble(pi)/HARTREE_EV;
	Py_DECREF(pi);
      }
      err = SetUsrCIEGridDetail(ng, xg);
    } else {
      return NULL;
    }
  } else if (n == 3 || n == 4) {
    if (!PyArg_ParseTuple(args, "idd|d", &ng, &emin, &emax, &eth)) 
      return NULL;
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetUsrCIEGrid(ng, emin, emax, eth);
  } else {
    return NULL;
  }

  if (err < 0) return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PFreeIonizationQk(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("FreeIonizationQk", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  FreeIonizationQk();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSaveIonizationQk(PyObject *self, PyObject *args) {
  int n;
  char *s;
  
  if (sfac_file) {
    SFACStatement("SaveIonizationQk", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "is", &n, &s)) return NULL;
  SaveCIRadialQkIntegrated(n, s);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PLoadIonizationQk(PyObject *self, PyObject *args) {
  int n;
  char *s;
  
  if (sfac_file) {
    SFACStatement("LoadIonizationQk", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  s = NULL;
  n = -1;
  if (!PyArg_ParseTuple(args, "|is", &n, &s)) return NULL;
  if (LoadCIRadialQkIntegrated(n, s) < 0) {
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PPrepIonizationQk(PyObject *self, PyObject *args) {
  PyObject *pz, *pa, *pn, *p;
  int nz, na, nn, nte;
  double emin, emax;
  char *s;
  int *n, i;
  double *z, *a;
  
  if (sfac_file) {
    SFACStatement("PrepIonizationQk", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  pz = NULL;
  pa = NULL;
  pn = NULL;
  nte = -1;
  emin = -1.0;
  emax = -1.0;
  if (!PyArg_ParseTuple(args, "s|OOOidd", 
			&s, &pz, &pa, &pn, &nte, &emin, &emax))
    return NULL;
  if (emin < 0) emin = 0.7;
  if (emax < 0) emax = 1.3;
  if (nte <= 1) {
    nte = 1;
    emin = 1.0;
    emax = 1.0;
  }
  if (pz) {
    if (!PyList_Check(pz)) return NULL;
    nz = PyList_Size(pz);
    z = (double *) malloc(sizeof(double)*nz);
    for (i = 0; i < nz; i++) {
      p = PyList_GetItem(pz, i);
      z[i] = PyInt_AsLong(p);
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
  
  if (pa) {
    if (!PyList_Check(pa)) return NULL;
    na = PyList_Size(pa);
    a = (double *) malloc(sizeof(double)*na);
    for (i = 0; i < na; i++) {
      p = PyList_GetItem(pa, i);
      a[i] = PyFloat_AsDouble(p);
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

  if (pn) {
    if (!PyList_Check(pn)) return NULL;
    nn = PyList_Size(pn);
    n = (int *) malloc(sizeof(int)*nn);
    for (i = 0; i < nn; i++) {
      p = PyList_GetItem(pn, i);
      n[i] = PyInt_AsLong(p);
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
  PrepCIRadialQkIntegrated(nz, z, na, a, nn, n, nte, emin, emax, s);
 
  free(z);
  free(a);
  free(n);

  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetCIPWOptions(PyObject *self, PyObject *args) {
  int qr, max, max_1, kl_cb;
  double tol;

  qr = 0;
  max = 500;
  max_1 = 8;
  kl_cb = 50;
  tol = 5E-2;
  if (!PyArg_ParseTuple(args, "d|iiii", &tol, &max, &max_1, &qr, &kl_cb)) 
    return NULL;
  SetCIPWOptions(qr, max, max_1, kl_cb, tol);
  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetCIPWGrid(PyObject *self, PyObject *args) {
  int ns, i;
  int n;
  PyObject *p, *q;
  int *m, *step;

  n = PyTuple_Size(args);
  if (n == 1) {    
    if (!PyArg_ParseTuple(args, "i", &ns)) return NULL;
    SetCIPWGrid(-ns, NULL, NULL);
  } else {
    if (!PyArg_ParseTuple(args, "OO", &p, &q)) return NULL;
    if (!PyList_Check(p) || !PyList_Check(q)) return NULL;
    ns = PyList_Size(p);
    if (ns != PyList_Size(q)) return NULL;
    if (ns == 0) return NULL;
    m = (int *) malloc(ns*sizeof(int));
    step = (int *) malloc(ns*sizeof(int));
    for (i = 0; i < ns; i++) {
      m[i] = PyInt_AsLong(PyList_GetItem(p, i));
      step[i] = PyInt_AsLong(PyList_GetItem(q, i));
    }
    SetCIPWGrid(ns, m, step);
    free(m);
    free(step);
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PCITable(PyObject *self, PyObject *args) { 
  int nlow, *low, nup, *up;
  char *s;
  PyObject *p, *q;
 
  if (sfac_file) {
    SFACStatement("CITable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "sOO", &s, &p, &q)) return NULL;
  nlow = SelectLevels(p, &low);
  if (nlow <= 0) return NULL;
  nup = SelectLevels(q, &up);
  if (nup <= 0) {
    free(low);
    return NULL;
  }
  SaveIonization(nlow, low, nup, up, s);
  free(low);
  free(up);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTestIntegrate(PyObject *self, PyObject *args) { 
  char *s;

  if (sfac_file) {
    SFACStatement("TestIntegrate", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
  TestIntegrate(s);
  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject *PTestCoulomb(PyObject *self, PyObject *args) { 
  char *s;
  
  if (sfac_file) {
    SFACStatement("TestCoulomb", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
  TestCoulomb(s);
  Py_INCREF(Py_None);
  return Py_None;  
}

  
static PyObject *PCorrectEnergy(PyObject *self, PyObject *args) {
  char *s;
  PyObject *p, *q, *ip, *iq;
  int n, k[MAX_ENERGY_CORRECTION];
  double e[MAX_ENERGY_CORRECTION];
  int i;
  FILE *f;

  if (sfac_file) {
    SFACStatement("CorrectEnergy", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "s", &s)) {
      printf("A single argument for CorrectEnergy must be a file name\n");
      return NULL;
    }
    f = fopen(s, "r");
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
  } else if (n == 2) {
    if (!PyArg_ParseTuple(args, "OO", &p, &q)) return NULL;
    if (!PyList_Check(p) || !PyList_Check(q)) {
      printf("The two arguments for CorrectEnergy must be two Lists\n");
      return NULL;
    }
    n = PyList_Size(p);
    if (PyList_Size(q) != n) {
      printf("The energy list length not equal the index list length\n");
      return NULL;
    }
    for (i = 0; i < n; i++) {
      ip = PyList_GetItem(p, i);
      iq = PyList_GetItem(q, i);
      k[i] = PyInt_AsLong(ip);
      e[i] = PyFloat_AsDouble(iq);
      e[i] /= HARTREE_EV;
    }
  } else {
    printf("CorrectEnergy takes one or two arrguments\n");
    return NULL;
  }

  CorrectEnergy(n, k, e);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PInfo(PyObject *self, PyObject *args) { 

  if (sfac_file) {
    SFACStatement("Info", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  Info();

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PPrintTable(PyObject *self, PyObject *args) { 
  char *fn1, *fn2;
  int v;
  
  if (sfac_file) {
    SFACStatement("PrintTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  v = 0;
  if (!PyArg_ParseTuple(args, "ss|i", &fn1, &fn2, &v)) return NULL;
  PrintTable(fn1, fn2, v);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMemENTable(PyObject *self, PyObject *args) { 
  char *fn;
  
  if (sfac_file) {
    SFACStatement("MemENTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  MemENTable(fn);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeMemENTable(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("FreeMemENTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  FreeMemENTable();

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PReinit(PyObject *self, PyObject *args, PyObject *keywds) {
  PyObject *q;
  int m_config;
  int m_recouple;
  int m_radial;
  int m_dbase;
  int m_structure;
  int m_excitation;
  int m_recombination;
  int m_ionization;
  int m;

  if (sfac_file) {
    SFACStatement("Reinit", args, keywds);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m_config = -1;
  m_recouple = -1;
  m_radial = -1;
  m_dbase = -1;
  m_structure = -1;
  m_excitation = -1;
  m_recombination = -1;
  m_ionization = -1;

  if (PyDict_Size(keywds) == 0) {
    m = 0;
    if (!PyArg_ParseTuple(args, "|i", &m)) return NULL;
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
    q = PyDict_GetItemString(keywds, "config");
    if (q) {
      if (!PyInt_Check(q)) return NULL;
      m_config = PyInt_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "recouple");
    if (q) {
      if (!PyInt_Check(q)) return NULL;
      m_recouple = PyInt_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "dbase");
    if (q) {
      if (!PyInt_Check(q)) return NULL;
      m_dbase = PyInt_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "structure");
    if (q) {
      if (!PyInt_Check(q)) return NULL;
      m_structure = PyInt_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "excitation");
    if (q) {
      if (!PyInt_Check(q)) return NULL;
      m_excitation = PyInt_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "radial");
    if (q) {
      if (!PyInt_Check(q)) return NULL;
      m_radial = PyInt_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "recombination");
    if (q) {
      if (!PyInt_Check(q)) return NULL;
      m_recombination = PyInt_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "ionization");
    if (q) {
      if (!PyInt_Check(q)) return NULL;
      m_ionization = PyInt_AsLong(q);
    }
  }

  ReinitFac(m_config, m_recouple, m_radial, m_dbase,
	    m_structure, m_excitation, m_recombination, m_ionization);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PPrint(PyObject *self, PyObject *args) {
  PyObject *p, *q;
  char *s;
  int i, n;

  if (sfac_file) {
    SFACStatement("Print", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);

  for (i = 0; i < n; i++) {
    p = PyTuple_GetItem(args, i);
    q = PyObject_Str(p);
    s = PyString_AsString(q);
    printf("%s", s);
    if (i != n-1) {
      printf(", ");
    }
    Py_XDECREF(q);
  }
  
  if (n > 0) printf("\n");

  fflush(stdout);

  Py_INCREF(Py_None);
  return Py_None;
} 

static struct PyMethodDef fac_methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"Config", (PyCFunction) PConfig, METH_VARARGS|METH_KEYWORDS},
  {"Closed", PClosed, METH_VARARGS},
  {"AvgConfig", PAvgConfig, METH_VARARGS},
  {"AddConfig", PAddConfig, METH_VARARGS},
  {"AITable", PAITable, METH_VARARGS},
  {"BasisTable", PBasisTable, METH_VARARGS},
  {"CETable", PCETable, METH_VARARGS},
  {"CETableMSub", PCETableMSub, METH_VARARGS},
  {"CITable", PCITable, METH_VARARGS},
  {"ClearLevelTable", PClearLevelTable, METH_VARARGS},
  {"ClearOrbitalTable", PClearOrbitalTable, METH_VARARGS},
  {"CloseSFAC", PCloseSFAC, METH_VARARGS},
  {"ConvertToSFAC", PConvertToSFAC, METH_VARARGS},
  {"CorrectEnergy", PCorrectEnergy, METH_VARARGS},
  {"DROpen", PDROpen, METH_VARARGS},
  {"DRTable", PDRTable, METH_VARARGS},
  {"FreeAngZ", PFreeAngZ, METH_VARARGS},
  {"FreeExcitationPk", PFreeExcitationPk, METH_VARARGS},
  {"FreeExcitationQk", PFreeExcitationQk, METH_VARARGS},
  {"FreeIonizationQk", PFreeIonizationQk, METH_VARARGS},
  {"FreeMemENTable", PFreeMemENTable, METH_VARARGS},
  {"FreeMultipole", PFreeMultipole, METH_VARARGS},
  {"FreeOrbitals", PFreeOrbitals, METH_VARARGS},
  {"FreeSlater", PFreeSlater, METH_VARARGS},
  {"FreeResidual", PFreeResidual, METH_VARARGS},
  {"FreeRecPk", PFreeRecPk, METH_VARARGS},
  {"FreeRecQk", PFreeRecQk, METH_VARARGS},
  {"FreeRecAngZ", PFreeRecAngZ, METH_VARARGS},
  {"GetCFPOld", PGetCFPOld, METH_VARARGS},
  {"GetW3j", PGetW3j, METH_VARARGS},
  {"GetW6j", PGetW6j, METH_VARARGS},
  {"GetW9j", PGetW9j, METH_VARARGS},
  {"GetCG", PGetCG, METH_VARARGS},
  {"GetPotential", PGetPotential, METH_VARARGS},
  {"Info", PInfo, METH_VARARGS},
  {"LoadIonizationQk", PLoadIonizationQk, METH_VARARGS},
  {"MemENTable", PMemENTable, METH_VARARGS},
  {"OptimizeRadial", POptimizeRadial, METH_VARARGS},
  {"PrepIonizationQk", PPrepIonizationQk, METH_VARARGS},
  {"PrintTable", PPrintTable, METH_VARARGS},
  {"RecStates", PRecStates, METH_VARARGS},
  {"Reinit", (PyCFunction) PReinit, METH_VARARGS|METH_KEYWORDS},
  {"RRTable", PRRTable, METH_VARARGS},
  {"SaveIonizationQk", PSaveIonizationQk, METH_VARARGS},
  {"SaveOrbitals", PSaveOrbitals, METH_VARARGS},
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
  {"SetCIPWOptions", PSetCIPWOptions, METH_VARARGS},
  {"SetCIPWGrid", PSetCIPWGrid, METH_VARARGS},
  {"SetCIQkMode", PSetCIQkMode, METH_VARARGS},
  {"SetHydrogenicNL", PSetHydrogenicNL, METH_VARARGS},
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
  {"Spline", PSpline, METH_VARARGS},
  {"Splint", PSplint, METH_VARARGS},
  {"Structure", PStructure, METH_VARARGS},
  {"TestAngular", PTestAngular, METH_VARARGS},
  {"TestCoulomb", PTestCoulomb, METH_VARARGS}, 
  {"TestIntegrate", PTestIntegrate, METH_VARARGS}, 
  {"TestMyArray", PTestMyArray, METH_VARARGS},  
  {"TestSpline", PTestSpline, METH_VARARGS},     
  {"TransitionTable", PTransitionTable, METH_VARARGS},  
  {"WaveFuncTable", PWaveFuncTable, METH_VARARGS},  
  {NULL, NULL}
};


void initfac(void) {
  PyObject *m, *d;
  char v[10];
  char sp[2];
  char *ename;
  double *emass;
  int i;

  m = Py_InitModule("fac", fac_methods);
  
  d = PyModule_GetDict(m);
  ErrorObject = Py_BuildValue("s", "fac.error");
  PyDict_SetItemString(d, "error", ErrorObject);

  if (InitFac() < 0) {
    onError("initilization failed\n");
    return;
  }

  SPECSYMBOL = PyList_New(MAX_SPEC_SYMBOLS);
  sp[1] = '\0';
  for (i = 0; i < MAX_SPEC_SYMBOLS; i++) {
    SpecSymbol(sp, i);
    PyList_SetItem(SPECSYMBOL, i, Py_BuildValue("s", sp));
  }

  sprintf(v, "%d.%d.%d", VERSION, SUBVERSION, SUBSUBVERSION);
  PFACVERSION = PyString_FromString(v);

  ename = GetAtomicSymbolTable();
  emass = GetAtomicMassTable();
  ATOMICSYMBOL = PyList_New(N_ELEMENTS+1);
  ATOMICMASS = PyList_New(N_ELEMENTS+1);
  PyList_SetItem(ATOMICSYMBOL, 0, Py_BuildValue("s", ""));
  PyList_SetItem(ATOMICMASS, 0, Py_BuildValue("d", 0.0));
  
  for (i = 0; i < N_ELEMENTS; i++) {
    PyList_SetItem(ATOMICSYMBOL, i+1, Py_BuildValue("s", &(ename[i*3])));
    PyList_SetItem(ATOMICMASS, i+1, Py_BuildValue("d", emass[i]));
  }

  QKMODE = PyDict_New();
  PyDict_SetItemString(QKMODE, "DEFAULT", Py_BuildValue("i", QK_DEFAULT));
  PyDict_SetItemString(QKMODE, "EXACT", Py_BuildValue("i", QK_EXACT));
  PyDict_SetItemString(QKMODE, "INTERPOLATE", 
		       Py_BuildValue("i", QK_INTERPOLATE));
  PyDict_SetItemString(QKMODE, "FIT", Py_BuildValue("i", QK_FIT));
  PyDict_SetItemString(QKMODE, "CB", Py_BuildValue("i", QK_CB));
  PyDict_SetItemString(QKMODE, "DW", Py_BuildValue("i", QK_DW));
  PyDict_SetItemString(QKMODE, "BED", Py_BuildValue("i", QK_BED));
  PyDict_SetItemString(QKMODE, "default", Py_BuildValue("i", QK_DEFAULT));
  PyDict_SetItemString(QKMODE, "exact", Py_BuildValue("i", QK_EXACT));
  PyDict_SetItemString(QKMODE, "interpolate", 
		       Py_BuildValue("i", QK_INTERPOLATE));
  PyDict_SetItemString(QKMODE, "fit", Py_BuildValue("i", QK_FIT));
  PyDict_SetItemString(QKMODE, "cb", Py_BuildValue("i", QK_CB));
  PyDict_SetItemString(QKMODE, "dw", Py_BuildValue("i", QK_DW));
  PyDict_SetItemString(QKMODE, "bed", Py_BuildValue("i", QK_BED));

  PyDict_SetItemString(d, "VERSION", PFACVERSION);
  PyDict_SetItemString(d, "SPECSYMBOL", SPECSYMBOL);
  PyDict_SetItemString(d, "ATOMICSYMBOL", ATOMICSYMBOL);
  PyDict_SetItemString(d, "ATOMICMASS", ATOMICMASS);
  PyDict_SetItemString(d, "QKMODE", QKMODE);

  if (PyErr_Occurred()) 
    Py_FatalError("can't initialize module fac");
}



