#include "Python.h"
#include <stdio.h>
#include <string.h>

#include "array.h"
#include "coulomb.h"
#include "config.h"
#include "cfp.h"
#include "angular.h"
#include "recouple.h"
#include "radial.h"
#include "nucleus.h"
#include "structure.h"
#include "transition.h"
#include "excitation.h"
#include "recombination.h"
#include "ionization.h"

static char *rcsid="$Id: fac.c,v 1.3 2001/10/12 18:49:30 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static PyObject *ErrorObject;

#define onError(message) {PyErr_SetString(ErrorObject, message);}

static PyObject *PSetOptimizeControl(PyObject *self, PyObject *args) {
  int maxiter;
  double tol; 
  int iprint;
  
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

static PyObject *PSetAvgConfig(PyObject *self, PyObject *args) {
  PyObject *acfg, *shell;
  int ns, i, m, kl, j;
  int *n, *kappa;
  double *nq, a;  

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

  Py_INCREF(Py_None);
  return Py_None;
}    
   
/** add a configuration to the list **/
static PyObject *PAddConfig(PyObject *self, PyObject *args) {
  CONFIG *cfg = NULL;
  char *group_name;
  int k;
  PyObject *python_cfg;

  if (!PyArg_ParseTuple(args, "sO", &group_name, &python_cfg)) {
    goto ERROR;
  }
  if (ConfigPythonToC(python_cfg, &cfg) < 0) goto ERROR;

  if (Couple(cfg) < 0) goto ERROR;

  k = GroupIndex(group_name);
  if (k < 0) goto ERROR;

  if (AddConfigToList(k, cfg) < 0) goto ERROR;
  
  Py_INCREF(Py_None);
  return Py_None;

 ERROR:
  if (cfg != NULL) free(cfg);
  return NULL;
}
 

static PyObject *PSetRadialGrid(PyObject *self, PyObject *args) {
  double rmax, rmin;

  rmin = 0.0;
  rmax = 0.0;
  if (!PyArg_ParseTuple(args, "dd", &rmin, &rmax))
    return NULL;
  SetRadialGrid(rmin, rmax);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetTransitionCut(PyObject *self, PyObject *args) {
  double c;
  if (!PyArg_ParseTuple(args, "d", &c))
    return NULL;
  SetTransitionCut(c);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetAICut(PyObject *self, PyObject *args) {
  double c;
  if (!PyArg_ParseTuple(args, "d", &c))
    return NULL;
  SetAICut(c);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetMaxRank(PyObject *self, PyObject *args) {
  int k;
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

  c = EPS4;
  mc = EPS3;
  if (!PyArg_ParseTuple(args, "i|dd", &n, &mc, &c))
    return NULL;
  SetAngZOptions(n, mc, c);
  Py_INCREF(Py_None);
  return Py_None;
}

/** coeff. of fractional parentage **/
static PyObject *GetCFPOld(PyObject *self, PyObject *args) {
  int j2, q, dj, dw, pj, pw;
  double coeff = 0.0;
  
  if (!PyArg_ParseTuple(args, "iiiiii", &j2, &q, &dj, &dw, &pj, &pw))
    return NULL;
  if (CFP(&coeff, j2, q, dj, dw, pj, pw) == -1)
    return NULL;
  return Py_BuildValue("d", coeff);
}

/** 3j symbol **/
static PyObject *GetW3j(PyObject *self, PyObject *args) {
  int j1, j2, j3, m1, m2, m3;
  if (!PyArg_ParseTuple(args, "iiiiii", &j1, &j2, &j3, &m1, &m2, &m3))
    return NULL;
  return Py_BuildValue("d", W3j(j1, j2, j3, m1, m2, m3));
}

/** 6j symbol **/
static PyObject *GetW6j(PyObject *self, PyObject *args) {
  int j1, j2, j3, i1, i2, i3;
  if (!PyArg_ParseTuple(args, "iiiiii", &j1, &j2, &j3, &i1, &i2, &i3))
    return NULL;
  return Py_BuildValue("d", W6j(j1, j2, j3, i1, i2, i3));
}

/** 9j symbol **/
static PyObject *GetW9j(PyObject *self, PyObject *args) {
  int j1, j2, j3, i1, i2, i3, k1, k2, k3;
  if (!PyArg_ParseTuple(args, "iiiiiiiii", 
			&j1, &j2, &j3, &i1, &i2, &i3, &k1, &k2, &k3))
    return NULL;
  return Py_BuildValue("d", W9j(j1, j2, j3, i1, i2, i3, k1, k2, k3));
}

/** clebsch gordan coeff. **/
static PyObject *GetCG(PyObject *self, PyObject *args) {
  int j1, j2, j3, m1, m2, m3;
  if (!PyArg_ParseTuple(args, "iiiiii", &j1, &j2, &j3, &m1, &m2, &m3))
    return NULL;
  return Py_BuildValue("d", ClebschGordan(j1, j2, j3, m1, m2, m3));
}

static PyObject *PSetAtom(PyObject *self, PyObject *args) {
  char *s;
  double z, mass;

  mass = 0.0;
  z = 0.0;
  if (!PyArg_ParseTuple(args, "s|dd", &s, &z, &mass)) return NULL;
  if (SetAtom(s, z, mass) < 0) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static int DecodeGroupArgs(PyObject *args, int **kg) {
  PyObject *p;
  char *s;
  int i, k, ng;  

  if (!PyList_Check(args) && !PyTuple_Check(args)) return -1;
  ng = PySequence_Length(args);
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
      s = PyString_AS_STRING(p);
      Py_DECREF(p);
      k = GroupExists(s);
      
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
  char *s;
  int ng, i, k;
  int *kg;
  double z;
  PyObject *p;
  double *weight;

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
	onError("wieghts must be a sequence");
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

  if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
  GetPotential(s);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSolveBound(PyObject *self, PyObject *args) {
  int n, kappa;
  ORBITAL *orb;
  int k;

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
  int n;
  PyObject *p, *q;

  p = NULL;
  q = NULL;
  n = PyTuple_Size(args);
  ngp = 0;
  kgp = NULL;
  if (n == 1 || n == 2) {
    PyArg_ParseTuple(args, "O|O", &p, &q);
    if (PyTuple_Check(p) || PyList_Check(p)) {
      ng = DecodeGroupArgs(p, &kg);
      if (q) {
	if (!PyTuple_Check(q) && !PyList_Check(q)) return NULL;
	if (PySequence_Length(q) > 0) {
	  ngp = DecodeGroupArgs(q, &kgp);
	}
      }
    } else {
      ng = DecodeGroupArgs(args, &kg);      
    }
  } else {
    ng = DecodeGroupArgs(args, &kg);
  }

  if (ng < 0 || ngp < 0) return NULL;
  
  nlevels = GetNumLevels();
  ns = MAX_SYMMETRIES;
  for (i = 0; i < ns; i++) {
    k = ConstructHamilton(i, ng, kg, ngp, kgp);
    if (k < 0) continue;
    if (DiagnolizeHamilton() < 0) return NULL;
    AddToLevels();
  }
  SortLevels(nlevels, -1);
  if (ng > 0) free(kg);
  if (ngp > 0) free(kgp);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PLevelTable(PyObject *self, PyObject *args) {
  char *fn;
  int n, m;

  n = 0;
  m = 0;
  if (!PyArg_ParseTuple(args, "s|ii", &fn, &n, &m)) return NULL;
  if (SaveLevelsToAscii(fn, m, n) < 0) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PClearLevelTable(PyObject *self, PyObject *args) {
  ClearLevelTable();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSortLevels(PyObject *self, PyObject *args) {
  char *fn;
  SortLevels(0, 0);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTestAngular(PyObject *self, PyObject *args) {
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

  max_e = 4;
  max_m = 4;
  if (!PyArg_ParseTuple(args, "ii|ii", &gauge, &mode, &max_e, &max_m)) 
    return NULL;
  SetTransitionOptions(gauge, mode, max_e, max_m);
  Py_INCREF(Py_None);
  return Py_None;
} 
  
static PyObject *PTransitionAll(PyObject *self, PyObject *args) {
  char *s;
  int m;
  m = 0;
  if (!PyArg_ParseTuple(args, "s|i", &s, &m)) return NULL;
  SaveTransition(0, NULL, 0, NULL, s, m);

  Py_INCREF(Py_None);
  return Py_None;
}


static int SelectLevels(PyObject *p, int **t) {
  int n, ng, *kg, i, j, k, im, kb, m, m0;
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
	im = lev->major_component;
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
      return k;
    } else if (PyList_Check(q)) {
      ng = DecodeGroupArgs(q, &kg);
      if (ng <= 0) return -1;
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
	  im = lev->major_component;
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
    if (!PyArg_ParseTuple(args, "OOs", &p, &q, &s)) return NULL;
    nlow = SelectLevels(p, &low);
    if (nlow <= 0) return NULL;
    nup = SelectLevels(q, &up);
    if (nup <= 0) return NULL;
    SaveTransition(nlow, low, nup, up, s, m);
    free(low);
    free(up);
  } else if (n == 4) {
    if (!PyArg_ParseTuple(args, "OOsi", &p, &q, &s, &m)) {
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
  int n, m;
  int nlow, nup, *low, *up;
  PyObject *p, *q;

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;

  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
    SaveExcitation(nlow, low, nup, up, 0, s);
  } else if (n == 2) {
    if (!PyArg_ParseTuple(args, "Os", &p, &s)) return NULL;
    nlow = SelectLevels(p, &low);
    if (nlow <= 0) return NULL;
    SaveExcitation(nlow, low, nlow, low, 0, s);
    free(low);
  } else if (n == 3) {
    if (!PyArg_ParseTuple(args, "OOs", &p, &q, &s)) return NULL;
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
  int n, m;
  int nlow, nup, *low, *up;
  PyObject *p, *q;

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;

  n = PyTuple_Size(args);
  if (n == 1) {
    if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
    SaveExcitation(nlow, low, nup, up, 1, s);
  } else if (n == 2) {
    if (!PyArg_ParseTuple(args, "Os", &p, &s)) return NULL;
    nlow = SelectLevels(p, &low);
    if (nlow <= 0) return NULL;
    SaveExcitation(nlow, low, nlow, low, 1, s);
    free(low);
  } else if (n == 3) {
    if (!PyArg_ParseTuple(args, "OOs", &p, &q, &s)) return NULL;
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
  
static PyObject *PSetCEFormat(PyObject *self, PyObject *args) {
  int m;

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCEFormat(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCEGridType(PyObject *self, PyObject *args) {
  int type;
  
  if (!PyArg_ParseTuple(args, "i", &type)) return NULL;
  SetCEEGridType(type);
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
  
  e = 0.0;
  if (!PyArg_ParseTuple(args, "sii|d", &s, &n, &k, &e)) return NULL;
  WaveFuncTable(s, n, k, e);

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
  if (!PyArg_ParseTuple(args, "i|i", &m1, &m2)) return NULL;
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
  if (n_max == 0) n_max = 50;

  SetRecSpectator(n_max, n_frozen, n_spec);
  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PRecStates(PyObject *self, PyObject *args) { 
  int ng;
  int *kg;
  int n;
  PyObject *gargs;

  if (!PyArg_ParseTuple(args, "iO", &n, &gargs)) return NULL;
  ng = DecodeGroupArgs(gargs, &kg);
  RecStates(n, ng, kg);
  
  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PSetRRFormat(PyObject *self, PyObject *args) {
  int m;

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetRRFormat(m);
  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PRRTable(PyObject *self, PyObject *args) { 
  int nlow, *low, nup, *up;
  int n, m;
  char *s;
  PyObject *p, *q;
  
  m = -1;
  if (!PyArg_ParseTuple(args, "OOs|i", &p, &q, &s, &m)) {
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

static PyObject *PSetPEGrid(PyObject *self, PyObject *args) {  
  int n;
  double xg[MAXNE];
  int ng;
  double emin, emax, eth;
  PyObject *p, *pi;
  int i, err;
  
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
  c = 0;
  if (!PyArg_ParseTuple(args, "OOs|i", &p, &q, &s, &c)) return NULL;
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

  if (!PyArg_ParseTuple(args, "OOOOss|i", &p, &q, &r, &t, &s1, &s2, &c)) 
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
  
  ArrayInit(&a, sizeof(double), 100);
  d = 0.1;
  ArraySet(&a, 200, &d);
  ArraySet(&a, 100, &d);
  b = (double *) ArrayGet(&a, 100);
  printf("%f ", *b);
  b = (double *) ArrayGet(&a, 200);
  printf("%f \n", *b);
  ArrayFree(&a, 0);

  MultiInit(&ma, sizeof(double), 3, block);
  printf("%d %d\n", ma.esize, ma.ndim);
  for (i = 10; i < 15; i++) {
    for (j = 0; j < 5; j++) {
      for (m = 45; m < 46; m++) {
	k[0] = i;
	k[1] = j;
	k[2] = m;	
	b = (double *) MultiSet(&ma, k, NULL);
	*b = 0.2;
	b = (double *) MultiGet(&ma, k);
	printf("%d %d %d %f \n", i, j, m, *b);
      }
    }
  }
  MultiFree(&ma, 0);
  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PSaveOrbitals(PyObject *self, PyObject *args) {  
  int n, i, norbs;
  double e;
  PyObject *p;

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
  FreeResidualArray();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeSlater(PyObject *self, PyObject *args) {
  FreeSlaterArray();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeMultipole(PyObject *self, PyObject *args) {
  FreeMultipoleArray();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeRecPk(PyObject *self, PyObject *args) {
  FreeRecPk();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeRecQk(PyObject *self, PyObject *args) {
  FreeRecQk();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeExcitationPk(PyObject *self, PyObject *args) {
  int ie;
  ie = -1;
  if (!PyArg_ParseTuple(args, "|i", &ie)) return NULL;
  FreeExcitationPk(ie);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeExcitationQk(PyObject *self, PyObject *args) {
  FreeExcitationQk();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFreeRecAngZ(PyObject *self, PyObject *args) { 
  FreeRecAngZ();
  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *PFreeAngZ(PyObject *self, PyObject *args) { 
  PyObject *p;
  int i, m;
  int n, *kg;

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

static PyObject *PSetCIFormat(PyObject *self, PyObject *args) {
  int m;

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCIFormat(m);
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
  FreeIonizationQk();
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
 
  if (!PyArg_ParseTuple(args, "OOs", &p, &q, &s)) return NULL;
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

  if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
  TestIntegrate(s);
  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject *PTestCoulomb(PyObject *self, PyObject *args) { 
  char *s;
  
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

static struct PyMethodDef fac_methods[] = {
  {"AddConfig", PAddConfig, METH_VARARGS},
  {"AITable", PAITable, METH_VARARGS},
  {"BasisTable", PBasisTable, METH_VARARGS},
  {"CETable", PCETable, METH_VARARGS},
  {"CETableMSub", PCETableMSub, METH_VARARGS},
  {"CITable", PCITable, METH_VARARGS},
  {"ClearLevelTable", PClearLevelTable, METH_VARARGS},
  {"CorrectEnergy", PCorrectEnergy, METH_VARARGS},
  {"DROpen", PDROpen, METH_VARARGS},
  {"DRTable", PDRTable, METH_VARARGS},
  {"FreeAngZ", PFreeAngZ, METH_VARARGS},
  {"FreeExcitationPk", PFreeExcitationPk, METH_VARARGS},
  {"FreeExcitationQk", PFreeExcitationQk, METH_VARARGS},
  {"FreeIonizationQk", PFreeIonizationQk, METH_VARARGS},
  {"FreeMultipole", PFreeMultipole, METH_VARARGS},
  {"FreeOrbitals", PFreeOrbitals, METH_VARARGS},
  {"FreeSlater", PFreeSlater, METH_VARARGS},
  {"FreeResidual", PFreeResidual, METH_VARARGS},
  {"FreeRecPk", PFreeRecPk, METH_VARARGS},
  {"FreeRecQk", PFreeRecQk, METH_VARARGS},
  {"FreeRecAngZ", PFreeRecAngZ, METH_VARARGS},
  {"GetCFPOld", GetCFPOld, METH_VARARGS},
  {"Get3j", GetW3j, METH_VARARGS},
  {"Get6j", GetW6j, METH_VARARGS},
  {"Get9j", GetW9j, METH_VARARGS},
  {"GetCG", GetCG, METH_VARARGS},
  {"GetPotential", PGetPotential, METH_VARARGS},
  {"LevelTable", PLevelTable, METH_VARARGS},
  {"OptimizeRadial", POptimizeRadial, METH_VARARGS},
  {"RecStates", PRecStates, METH_VARARGS},
  {"RRTable", PRRTable, METH_VARARGS},
  {"SaveOrbitals", PSaveOrbitals, METH_VARARGS},
  {"SetAICut", PSetAICut, METH_VARARGS},
  {"SetAngZOptions", PSetAngZOptions, METH_VARARGS},
  {"SetAtom", PSetAtom, METH_VARARGS},
  {"SetAvgConfig", PSetAvgConfig, METH_VARARGS},
  {"SetCEFormat", PSetCEFormat, METH_VARARGS},
  {"SetCEGrid", PSetCEGrid, METH_VARARGS},
  {"SetTEGrid", PSetTEGrid, METH_VARARGS},
  {"SetCEPWOptions", PSetCEPWOptions, METH_VARARGS},
  {"SetCEPWGrid", PSetCEPWGrid, METH_VARARGS},
  {"SetCIFormat", PSetCIFormat, METH_VARARGS},
  {"SetCIEGrid", PSetCIEGrid, METH_VARARGS},
  {"SetIEGrid", PSetIEGrid, METH_VARARGS},
  {"SetCIPWOptions", PSetCIPWOptions, METH_VARARGS},
  {"SetCIPWGrid", PSetCIPWGrid, METH_VARARGS},
  {"SetMaxRank", PSetMaxRank, METH_VARARGS},
  {"SetOptimizeControl", PSetOptimizeControl, METH_VARARGS},
  {"SetPEGrid", PSetPEGrid, METH_VARARGS},
  {"SetRadialGrid", PSetRadialGrid, METH_VARARGS},
  {"SetRecPWLimits", PSetRecPWLimits, METH_VARARGS},
  {"SetRecPWOptions", PSetRecPWOptions, METH_VARARGS},
  {"SetRecSpectator", PSetRecSpectator, METH_VARARGS},
  {"SetRRFormat", PSetRRFormat, METH_VARARGS},
  {"SetRRTEGrid", PSetRRTEGrid, METH_VARARGS},
  {"SetScreening", PSetScreening, METH_VARARGS},
  {"SetTransitionCut", PSetTransitionCut, METH_VARARGS},
  {"SetTransitionOptions", PSetTransitionOptions, METH_VARARGS},
  {"SetUsrCEGrid", PSetUsrCEGrid, METH_VARARGS},
  {"SetUsrCEGridType", PSetUsrCEGridType, METH_VARARGS},
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
  {"TransitionAll", PTransitionAll, METH_VARARGS},  
  {"TransitionTable", PTransitionTable, METH_VARARGS},  
  {"WaveFuncTable", PWaveFuncTable, METH_VARARGS},  
  {NULL, NULL}
};


void initfac() {
  PyObject *m, *d;

  m = Py_InitModule("fac", fac_methods);
  
  d = PyModule_GetDict(m);
  ErrorObject = Py_BuildValue("s", "fac.error");
  PyDict_SetItemString(d, "error", ErrorObject);

#if FAC_DEBUG
  debug_log = fopen("debug.log", "w");
  if (!debug_log) {
    onError("can not open the debug log file");
    return;
  }
#endif
  
  if (InitConfig() < 0) {
    onError("initialize failed in InitConfig");
    return;
  }

  InitCoulomb();
  InitAngular();
  InitRecouple();

  if (InitRadial() < 0) {
    onError("initialize failed in InitRadial");
    return;
  }

  InitStructure();
  InitExcitation();
  InitRecombination();
  InitIonization();

  if (PyErr_Occurred()) 
    Py_FatalError("can't initialize module fac");
}



