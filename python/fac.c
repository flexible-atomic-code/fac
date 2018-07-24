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

#include "Python.h"
#include <stdio.h>
#include <string.h>

#include "init.h"
#include "cf77.h"
#include "mpiutil.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif
 
#if PY_MAJOR_VERSION >= 3
  #define PyUnicode_AsString(x) PyBytes_AsString(PyUnicode_AsEncodedString((x), "utf-8", "strict"))
#else
  #define PyLong_AsLong PyInt_AsLong
  #define PyLong_AS_LONG PyInt_AS_LONG
  #define PyLong_Check PyInt_Check
  #define PyUnicode_FromString PyString_FromString
  #define PyUnicode_AsString PyString_AsString
  #define PyUnicode_Check PyString_Check
#endif

#if PY_MAJOR_VERSION >= 3
  #define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
  #define GETSTATE(m) (&_state)
  static struct module_state _state;
#endif

struct module_state {
    PyObject *error;
};

static PyObject *ErrorObject;
static PyObject *PFACVERSION;
static PyObject *SPECSYMBOL;
static PyObject *ATOMICSYMBOL;
static PyObject *ATOMICMASS;
static PyObject *QKMODE;

static PyObject *_thismodule;
static PyObject *_thisdict;

static FILE *sfac_file = NULL;

//#define onError(message) {PyErr_SetString(ErrorObject, message);}
#define onError(message) {printf("%s\n", message);}

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
  s1 = PyUnicode_AsString(sargs);
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
      s2 = PyUnicode_AsString(q);
      fprintf(sfac_file, "%s=", s2);
      q = PyTuple_GetItem(p, 1);
      kvar = PyObject_Str(q);
      s2 = PyUnicode_AsString(kvar);
      if (PyUnicode_Check(q)) {
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

static PyObject *PCheckEndian(PyObject *self, PyObject *args) {
  char *fn;
  TFILE *f;
  F_HEADER fh;
  int i, swp;

  if (sfac_file) {
    SFACStatement("CheckEndian", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  fn = NULL;
  if (!PyArg_ParseTuple(args, "|s", &fn)) return NULL;
  if (fn) {
    f = OpenFileRO(fn, &fh, &swp);
    if (f == NULL) {
      printf("Cannot open file %s\n", fn);
      return NULL;
    }
    i = CheckEndian(&fh);
  } else {
    i = CheckEndian(NULL);
  }
  printf("Endian: %d\n", i);

  Py_INCREF(Py_None);
  return Py_None;
}  

static int DoubleFromList(PyObject *p, double **k) {
  int i, n;
  PyObject *q;

  if (PyList_Check(p)) {
    n = PyList_Size(p);
    if (n > 0) {
      *k = malloc(sizeof(double)*n);
      for (i = 0; i < n; i++) {
	q = PyList_GetItem(p, i);
	(*k)[i] = PyFloat_AsDouble(q);
      }
    }
  } else {
    n = 1;
    *k = malloc(sizeof(double));
    (*k)[0] = PyFloat_AsDouble(p);
  } 
  return n;
}


static int IntFromList(PyObject *p, int **k) {
  int i, n;
  PyObject *q;

  if (PyList_Check(p)) {
    n = PyList_Size(p);
    if (n > 0) {
      *k = malloc(sizeof(int)*n);
      for (i = 0; i < n; i++) {
	q = PyList_GetItem(p, i);
	(*k)[i] = PyLong_AsLong(q);
      }
    }
  } else {
    n = 1;
    *k = malloc(sizeof(int));
    (*k)[0] = PyLong_AsLong(p);
  }
  return n;
}

static int DecodeGroupArgs(PyObject *args, int **kg, int *n0) {
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
    int n0p = 0;
    int n0q = 0;
    if (n0) {
      n0p = *n0;
      n0q = *n0;
    }
    int n = 0;
    for (i = 0; i < ng; i++) {
      p = PySequence_GetItem(args, i);
      if (!PyUnicode_Check(p)) {
	free((*kg));
	onError("argument must be a group name");
	return -1;
      }
      s = PyUnicode_AsString(p);
      k = GroupExists(s);
      Py_DECREF(p);
      
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
      onError("all cfg groups invalid");
      return -1;
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

static PyObject *PGetBoundary(PyObject *self, PyObject *args) {
  int nmax, ib;
  double bqp, rmax, dr;
  
  ib = GetBoundary(&rmax, &bqp, &nmax, &dr);
  
  return Py_BuildValue("[iiddd]", ib, nmax, rmax, bqp, dr);
}

static PyObject *PSetBoundary(PyObject *self, PyObject *args) {
  int nmax, ierr;
  double bqp, p;
  
  if (sfac_file) {
    SFACStatement("SetBoundary", args, NULL);
    Py_INCREF(Py_None); 
    return Py_None;
  }
  p = -1.0;
  bqp = 0.0;
  if (!PyArg_ParseTuple(args, "i|dd", &nmax, &p, &bqp))
    return NULL;
  ierr = SetBoundary(nmax, p, bqp);
  
  if (ierr < 0) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetOptimizeControl(PyObject *self, PyObject *args) {
  int maxiter;
  double tol, s; 
  int iprint;
  
  if (sfac_file) {
    SFACStatement("SetOptimizeControl", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
    
  iprint = 0;
  if (!PyArg_ParseTuple(args, "ddi|i", &tol, &s, &maxiter, &iprint))
    return NULL;
  SetOptimizeControl(tol, s, maxiter, iprint);

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
    if (!PyLong_Check(q)) {
      printf("Screened n must be integers\n");
      free(screened_n);
      Py_DECREF(q);
      return NULL;
    }
    screened_n[i] = PyLong_AsLong(q);
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
    
    shells[m].n = PyLong_AsLong(PyTuple_GetItem(python_shell, 0));
    k = PyLong_AsLong(PyTuple_GetItem(python_shell, 1));
    j = PyLong_AsLong(PyTuple_GetItem(python_shell, 2));
    if (j > 0) k = -(k+1);
    shells[m].kappa = k;
    shells[m].nq = PyLong_AsLong(PyTuple_GetItem(python_shell, 3));
  }

  return 0;
  
 ERROR:
  onError("error in conversion");
  if (shells) free(shells);
  //if (cfg) free(cfg);
  return -1;
}

static char _closed_shells[MCHSHELL] = "";
static PyObject *PClosed(PyObject *self, PyObject *args) {
  CONFIG *cfg;
  PyObject *q;
  int i, j, kl, n, nq, ncfg;
  char *p, argv[512];
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
    if (!PyUnicode_Check(q)) return NULL;
    p = PyUnicode_AsString(q);
    strncpy(argv, p, 512);
    ns = StrSplit(argv, ' ');
    p = argv;
    for (k = 0; k < ns; k++) {
      while (*p == ' ') p++;
      ncfg = GetConfigFromStringNR(&cfg, p);
      for (j = ncfg-1; j >= 0; j--) {
	if (cfg[j].n_shells != 1) return NULL;	
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

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PGetConfigNR(PyObject *self, PyObject *args) {
  CONFIG *cfg;
  PyObject *q, *r;
  int i, j, t, ncfg;
  char scfg[MCHSHELL], *p, s[16];
  int argc;

  r = Py_BuildValue("[]");
  argc = PyTuple_Size(args);
  for (i = 0; i < argc; i++) {
    q = PyTuple_GetItem(args, i);
    if (!PyUnicode_Check(q)) return NULL;
    p = PyUnicode_AsString(q);
    strncpy(scfg, _closed_shells, MCHSHELL);
    strncat(scfg, p, MCHSHELL);
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
      PyList_Append(r, Py_BuildValue("s", scfg));
    }
    if (ncfg > 0) free(cfg);
  }

  return r;
}

static PyObject *PReadConfig(PyObject *self, PyObject *args) {
  char *s, *c;
  
  if (sfac_file) {
    SFACStatement("ReadConfig", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  c = NULL;
  if (!(PyArg_ParseTuple(args, "s|s", &s, &c))) return NULL;
  if (ReadConfig(s, c) < 0) return NULL;
  
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PConfig(PyObject *self, PyObject *args, PyObject *keywds) {
  CONFIG *cfg;
  PyObject *q, *q1, *qb;
  static char gname[GROUP_NAME_LEN] = "_all_";
  int i, j, t, ncfg;
  char scfg[MCHSHELL], *p;
  int argc;

  if (sfac_file) {
    SFACStatement("Config", args, keywds);
    Py_INCREF(Py_None);
    return Py_None;
  }

  q = PyTuple_GET_ITEM(args, 0);
  if (PyLong_Check(q)) {
    int m, r, ng, *kg, ngb, *kgb, n0, n1, k0, k1, n0d, n1d;
    double sth;
    char *gn1, *gn2, *s;
    s = NULL;
    n0 = 1;
    n1 = 1;
    k0 = 0;
    k1 = -1;
    n0d = 0;
    n1d = 0;
    sth = 0;
    qb = NULL;
    ng = 0;
    kg = NULL;
    ngb = 0;
    kgb = NULL;
    gn1 = NULL;
    gn2 = NULL;
    m = PyLong_AsLong(q);
    if (m == 0) {
      if (!(PyArg_ParseTuple(args, "iss|dO", &m, &gn1, &s, &sth, &qb))) {
	return NULL;
      }
      if (qb == NULL) sth = 0.0;
      else {
	if (PyLong_Check(qb)) {
	  ngb = PyLong_AsLong(qb);
	  kgb = NULL;
	} else {
	  ngb = DecodeGroupArgs(qb, &kgb, NULL);
	}
      }
      strncpy(scfg, s, MCHSHELL);
      r = ConfigSD(m, ng, kg, scfg, gn1, NULL,
		   n0, n1, n0d, n1d, k0, k1, ngb, kgb, sth);
      if (ngb > 0 && kgb) free(kgb);
    } else {
      if (!(PyArg_ParseTuple(args, "iOO|siiiiiidO",
			     &m, &q1, &q, &s, &n0, &n1, &k0, &k1,
			     &n0d, &n1d, &sth, &qb))) return NULL;
      if (PyList_Check(q1)) {
	int nq1 = PyList_Size(q1);
	if (nq1 < 1 || nq1 > 2) return NULL;
	PyObject *q2 = PyList_GetItem(q1, 0);
	if (!PyUnicode_Check(q2)) return NULL;
	gn1 = PyUnicode_AsString(q2);
	if (nq1 > 1) {
	  q2 = PyList_GetItem(q1, 1);
	  if (!PyUnicode_Check(q2)) return NULL;
	  gn2 = PyUnicode_AsString(q2);
	} else {
	  gn2 = NULL;
	}
      } else if (PyUnicode_Check(q1)) {
	gn1 = PyUnicode_AsString(q1);
	gn2 = NULL;
      } else {
	return NULL;
      }
      ng = DecodeGroupArgs(q, &kg, NULL);
      ngb = ng;
      kgb = kg;
      if (qb != NULL) {
	if (PyLong_Check(qb)) {
	  ngb = PyLong_AsLong(qb);
	  kgb = NULL;
	} else {
	  ngb = DecodeGroupArgs(qb, &kgb, NULL);
	}
      }      
      strncpy(scfg, s, MCHSHELL);
      r = ConfigSD(m, ng, kg, scfg, gn1, gn2,
		   n0, n1, n0d, n1d, k0, k1, ngb, kgb, sth);      
      if (ngb > 0 && kgb && kgb != kg) free(kgb);
      if (ng > 0) free(kg);
    }
    if (r < 0) {
      return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
  }
  argc = PyTuple_Size(args);
  i = 0;

  if (keywds) {
    q = PyDict_GetItemString(keywds, "group");
    if (!q || !PyUnicode_Check(q)) {
      printf("The keyword must be group=gname\n");
      return NULL;
    }
    p = PyUnicode_AsString(q);
    strncpy(gname, p, GROUP_NAME_LEN);
  } else {
    if (argc == 0) return NULL;
    q = PyTuple_GetItem(args, i);
    i++;
    if (!PyUnicode_Check(q)) return NULL;
    p = PyUnicode_AsString(q);
    strncpy(gname, p, GROUP_NAME_LEN);
  }
  
  t = GroupIndex(gname);
  if (t < 0) return NULL;
  
  for (; i < argc; i++) {   
    q = PyTuple_GetItem(args, i);
    if (!PyUnicode_Check(q)) return NULL;
    p = PyUnicode_AsString(q);
    strncpy(scfg, _closed_shells, MCHSHELL);
    strncat(scfg, p, MCHSHELL);
    ncfg = GetConfigFromString(&cfg, scfg);

    for (j = 0; j < ncfg; j++) {
      if (Couple(cfg+j) < 0) return NULL;
      if (AddConfigToList(t, cfg+j) < 0) return NULL;
    }   
    if (ncfg > 0) free(cfg);
  }

  CONFIG_GROUP *g = GetGroup(t);
  if (g != NULL && g->n_cfgs == 0) {
    RemoveGroup(t);
  }
  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PRemoveConfig(PyObject *self, PyObject *args) {
  int k, ng, *kg;
  PyObject *p;

  if (sfac_file) {
    SFACStatement("RemoveConfig", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  ng = PyTuple_Size(args);
  if (ng <= 0) return NULL;
  p = PyTuple_GET_ITEM(args, 0);
  if (PyUnicode_Check(p)) {
    ng = DecodeGroupArgs(args, &kg, NULL);
  } else {
    ng = DecodeGroupArgs(p, &kg, NULL);
  }
  
  for (k = 0; k < ng; k++) {
    RemoveGroup(kg[k]);
  }
  ReinitStructure(1);
  ReinitRecouple(0);
  if (ng > 0) free(kg);

  Py_INCREF(Py_None);
  return Py_None;
}  
 
static PyObject *PListConfig(PyObject *self, PyObject *args) {
  int k, ng, *kg;
  PyObject *p;
  char *s;

  if (sfac_file) {
    SFACStatement("ListConfig", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  s = NULL;
  p = NULL;
  if (!PyArg_ParseTuple(args, "|sO", &s, &p)) return NULL;
  if (p) {
    ng = DecodeGroupArgs(p, &kg, NULL);
  } else {
    ng = 0;
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
 
static PyObject *PSolvePseudo(PyObject *self, PyObject *args) {
  int kmin, kmax, nb, nmax, nd;
  double xdf;
  
  if (sfac_file) {
    SFACStatement("SolvePseudo", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  kmin = 0;
  nb = 0;
  nd = 1;
  xdf = -1.0;
  if (!PyArg_ParseTuple(args, "ii|iiid",
			&kmax, &nmax, &kmin, &nb, &nd, &xdf)) return NULL;
  SolvePseudo(kmin, kmax, nb, nmax, nd, xdf);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetPotentialMode(PyObject *self, PyObject *args) {
  int m;
  double h, ih, h0, h1;

  if (sfac_file) {
    SFACStatement("SetPotentialMode", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  h = 1E31;
  ih = 1E31;
  h0 = -1;
  h1 = -1;
  if (!PyArg_ParseTuple(args, "i|dddd", &m, &h, &ih, &h0, &h1))
    return NULL;

  SetPotentialMode(m, h, ih, h0, h1);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetRadialGrid(PyObject *self, PyObject *args) {
  double ratio, asym, rmin, qr;
  int maxrp;

  if (sfac_file) {
    SFACStatement("SetRadialGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  qr = -1;
  if (!PyArg_ParseTuple(args, "iddd|d", &maxrp, &ratio, &asym, &rmin, &qr))
    return NULL;
  if (-1 == SetRadialGrid(maxrp, ratio, asym, rmin, qr))
    return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetTransitionCut(PyObject *self, PyObject *args) {
  double c0, c;

  if (sfac_file) {
    SFACStatement("SetTransitionCut", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  c = -1.0;
  if (!PyArg_ParseTuple(args, "d|d", &c0, &c))
    return NULL;
  SetTransitionCut(c0, c);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetSE(PyObject *self, PyObject *args) {
  int c, m, s, p;
  
  if (sfac_file) {
    SFACStatement("SetSE", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = -1;
  s = -1;
  p = -1;
  if (!PyArg_ParseTuple(args, "i|iii", &c, &m, &s, &p))
    return NULL;
  SetSE(c, m, s, p);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetModSE(PyObject *self, PyObject *args) {
  double o0, o1, a, c0, c1, c;
  
  if (sfac_file) {
    SFACStatement("SetModSE", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  o0 = -1;
  o1 = -1;
  a = -1;
  c0 = -1;
  c1 = -1;
  c = -1;
  if (!PyArg_ParseTuple(args, "ddd|ddd", &o0, &o1, &a, &c0, &c1, &c))
    return NULL;
  SetModSE(o0, o1, a, c0, c1, c);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetVP(PyObject *self, PyObject *args) {
  int c;

  if (sfac_file) {
    SFACStatement("SetVP", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &c))
    return NULL;
  SetVP(c);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetBreit(PyObject *self, PyObject *args) {
  int c, m, n, k;
  double x;

  if (sfac_file) {
    SFACStatement("SetBreit", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = -1;
  n = -1;
  x = -1.0;
  k = 0;
  if (!PyArg_ParseTuple(args, "i|iidi", &c, &m, &n, &x, &k))
    return NULL;
  SetBreit(c, m, n, x, k);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetMS(PyObject *self, PyObject *args) {
  int c1, c2;

  if (sfac_file) {
    SFACStatement("SetMS", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "ii", &c1, &c2))
    return NULL;
  SetMS(c1, c2);
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

  c = ANGZCUT;
  mc = MIXCUT;
  if (!PyArg_ParseTuple(args, "i|dd", &n, &mc, &c))
    return NULL;
  SetAngZOptions(n, mc, c);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCILevel(PyObject *self, PyObject *args) {
  int i;

  if (sfac_file) {
    SFACStatement("SetCILevel", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &i))
    return NULL;
  SetCILevel(i);
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
  double c, c2;

  if (sfac_file) {
    SFACStatement("SetMixCut", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  c2 = -1.0;
  if (!PyArg_ParseTuple(args, "d|d", &c, &c2))
    return NULL;
  SetMixCut(c, c2);
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

  if (!PyArg_ParseTuple(args, "iiiiii", &j1, &m1, &j2, &m2, &j3, &m3))
    return NULL;
  return Py_BuildValue("d", ClebschGordan(j1, m1, j2, m2, j3, m3));
}

static PyObject *PSetExtraPotential(PyObject *self, PyObject *args) {
  int m;
  int n;
  double *p;
  PyObject *t;
  
  if (sfac_file) {
    SFACStatement("SetExtraPotential", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  t = NULL;
  n = 0;
  p = NULL;
  if (!PyArg_ParseTuple(args, "i|O", &m, &t)) {
    return NULL;
  }
  if (m >= 0 && t != NULL) {
    n = DoubleFromList(t, &p);
  }
  SetExtraPotential(m, n, p);
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PSetAtom(PyObject *self, PyObject *args) {
  char *s;
  double z, mass, rn, a, npr;
  PyObject *t;
  
  if (sfac_file) {
    SFACStatement("SetAtom", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  mass = -1.0;
  z = -1.0;
  rn = -1.0;
  a = -1.0;
  npr = -1;
  if (!PyArg_ParseTuple(args, "O|ddddd", &t, &z, &mass, &rn, &a, &npr)) {
    return NULL;
  }
  if (PyUnicode_Check(t)) {
    s = PyUnicode_AsString(t);
    if (SetAtom(s, z, mass, rn, a, npr) < 0) return NULL;
  } else if (PyFloat_Check(t)) {
    npr = a;
    a = rn;
    rn = mass;
    mass = z;
    z = PyFloat_AsDouble(t);
    if (SetAtom(NULL, z, mass, rn, a, npr) < 0) return NULL;
  } else if (PyLong_Check(t)) {
    npr = a;
    a = rn;
    rn = mass;    
    mass = z;
    z = (double) PyLong_AsLong(t);
    if (SetAtom(NULL, z, mass, rn, a, npr) < 0) return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *POptimizeModSE(PyObject *self, PyObject *args) {
  int n, ka, ni;
  double dr;

  if (sfac_file) {
    SFACStatement("OptimizeModSE", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  ni = 1000000;
  if (!PyArg_ParseTuple(args, "iid|i", &n, &ka, &dr, &ni)) return NULL;
  OptimizeModSE(n, ka, dr, ni);
  
  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PSetHydrogenicNL(PyObject *self, PyObject *args) {
  int n, k, nm, km;
  
  if (sfac_file) {
    SFACStatement("SetHydrogenicNL", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = -1;
  k = -1;
  nm = -1;
  km = -1;
  if (!PyArg_ParseTuple(args, "|iiii", &n, &k, &nm, &km)) return NULL;

  SetHydrogenicNL(n, k, nm, km);
  Py_INCREF(Py_None);
  return Py_None;
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
  if (PyUnicode_Check(p)) {
    weight = NULL;
    ng = DecodeGroupArgs(args, &kg, NULL);
    if (ng < 0) {
      onError("the group argument format error");
      return NULL;
    }
  } else {
    ng = DecodeGroupArgs(p, &kg, NULL);
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
  if (ng == 0) {
    kg = NULL;
  }
  if (OptimizeRadial(ng, kg, -1, weight, 0) < 0) {
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

static PyObject *PFreezeOrbital(PyObject *self, PyObject *args) {
  if (sfac_file) {
    SFACStatement("FreezeOrbital", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  char *s;
  int m = -1;
  if (!PyArg_ParseTuple(args, "s|i", &s, &m)) return NULL;
  FreezeOrbital(s, m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PRefineRadial(PyObject *self, PyObject *args) {
  int maxfun, msglvl;

  if (sfac_file) {
    SFACStatement("RefineRadial", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  maxfun = 0;
  msglvl = 0;
  if (!PyArg_ParseTuple(args, "|ii", &maxfun, &msglvl)) return NULL;
  
  if (RefineRadial(maxfun, msglvl)) return NULL;
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
  int ng, i;
  int ip;
  int ngp;
  int *kg, *kgp;
  char *fn;
  PyObject *p, *q, *t, *s, *r;

  if (sfac_file) {
    SFACStatement("Structure", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  p = NULL;
  q = NULL;
  s = NULL;
  r = NULL;
  ngp = 0;
  kgp = NULL;
  ip = 0;
  
  if (!(PyArg_ParseTuple(args, "O|OOOO", &t, &p, &q, &s, &r))) return NULL;
  if (r != NULL) {
    if (PyLong_Check(r)) ip = PyLong_AsLong(r);  
  } else if (s != NULL && PyLong_Check(s)) {
    ip = PyLong_AsLong(s);
  }
  if (PyLong_Check(t)) {
    ip = PyLong_AsLong(t);
    if (p == NULL) {
      SetDiagMaxIter(ip, -1.0);
      Py_INCREF(Py_None);
      return Py_None;    
    }
    if (ip >= -1) {
      i = IntFromList(p, &kg);
      SetSymmetry(ip, i, kg);
      free(kg);
    } else {
      if (!PyFloat_Check(p)) return NULL;      
      double a = PyFloat_AsDouble(p);
      double b = 0.0;
      double c = 0.0;
      if (ip >= -100) {
	if (q != NULL && PyFloat_Check(q)) {
	  b = PyFloat_AsDouble(q);
	}
	if (s != NULL) {
	  c = PyFloat_AsDouble(s);
	}
	SetPerturbThreshold(-ip, a, b, c);
      } else {
	SetDiagMaxIter(-ip-100, a);
      }
    }
    Py_INCREF(Py_None);
    return Py_None;    
  }
  
  fn = PyUnicode_AsString(t);
  char *hfn = NULL;
  if (p) {
    if (PyUnicode_Check(p)) {
      hfn = PyUnicode_AsString(p);
      if (q == NULL) {
	ng = 0;
	ngp = 0;
	kg = NULL;
	kgp = NULL;
      } else {
	if (PyTuple_Check(q) || PyList_Check(q)) {
	  ng = DecodeGroupArgs(q, &kg, NULL);
	  if (ng < 0) return NULL;
	  if (s) {
	    if (!PyTuple_Check(s) && !PyList_Check(s)) return NULL;
	    if (PySequence_Length(s) > 0) {
	      ngp = DecodeGroupArgs(s, &kgp, NULL);
	    }
	  }
	} else {
	  return NULL;
	}
      }
    } else { 
      if (PyTuple_Check(p) || PyList_Check(p)) {
	ng = DecodeGroupArgs(p, &kg, NULL);
	if (ng < 0) return NULL;
	if (q) {
	  if (!PyTuple_Check(q) && !PyList_Check(q)) return NULL;
	  if (PySequence_Length(q) > 0) {
	    ngp = DecodeGroupArgs(q, &kgp, NULL);
	  }
	}
      } else {
	return NULL;
      }
    }
  } else {
    ng = DecodeGroupArgs(NULL, &kg, NULL);  
    if (ng < 0) return NULL;
  }

  if (SolveStructure(fn, hfn, ng, kg, ngp, kgp, ip) < 0) {
    onError("Diagnolizing Hamiltonian Error");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetUTA(PyObject *self, PyObject *args) {
  int m, mci;

  if (sfac_file) {
    SFACStatement("SetUTA", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  mci = 1;
  if (!PyArg_ParseTuple(args, "i|i", &m, &mci)) return NULL;
  
  SetUTA(m, mci);

  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PSetTRF(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetTRF", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  
  SetTRF(m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTestHamilton(PyObject *self, PyObject *args) {
  
  TestHamilton();
 
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

  SortLevels(0, 0, 0);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTestAngular(PyObject *self, PyObject *args) {

  if (sfac_file) {
    SFACStatement("TestAngular", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  TestAngular();

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
  int ig, nlevels, iuta;
  LEVEL *lev;
  SYMMETRY *sym;
  STATE *s;
  char rgn[GROUP_NAME_LEN];

  iuta = IsUTA();
  if (!PyList_Check(p) && !PyTuple_Check(p)) return 0;
  n = PySequence_Length(p);
  if (n > 0) {
    q = PySequence_GetItem(p, 0);
    if (PyUnicode_Check(q)) {
      ng = DecodeGroupArgs(p, &kg, NULL);
      if (ng <= 0) {
	return 0;
      }
      nlevels = GetNumLevels();
      (*t) = malloc(sizeof(int)*nlevels);
      if (!(*t)) return 0;
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
      Py_DECREF(q);
      (*t) = realloc(*t, k*sizeof(int));
      return k;
    } else if (PyList_Check(q)) {
      if (n != 2) {
	printf("recombined states specification unrecoganized\n");
	return -1;
      }
      ng = DecodeGroupArgs(q, &kg, NULL);
      if (ng <= 0) return -1;
      Py_DECREF(q);
      q = PySequence_GetItem(p, 1);
      if (PyList_Check(q)) {
	p = q;
	m0 = 0;
	n = PySequence_Length(q);
      } else if (PyLong_Check(q)) {
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
	nrec = PyLong_AS_LONG(q);
	Py_DECREF(q);
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
      return k;
    } else {
      (*t) = malloc(sizeof(int)*n);
      if (!(*t)) return 0;      
      for (i = 0; i < n; i++) {
	q = PySequence_GetItem(p, i);
	if (!PyLong_Check(q)) {
	  free(*t);
	  return 0;
	}
	Py_DECREF(q);
	(*t)[i] = PyLong_AS_LONG(q);
      }
      return n;
    }
  }
  return 0;
}
  
static PyObject *PCutMixing(PyObject *self, PyObject *args) {
  int nlev, n, *ilev, *kg;
  double c;
  PyObject *t, *p;

  if (sfac_file) {
    SFACStatement("CutMixing", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  nlev = 0;
  n = 0;
  c = 0.0;
  if (!(PyArg_ParseTuple(args, "OO|d", &t, &p, &c))) return NULL;
  nlev = SelectLevels(t, &ilev);
  if (nlev <= 0) goto DONE;
  n = DecodeGroupArgs(p, &kg, NULL);
  if (n <= 0) goto DONE;
  
  CutMixing(nlev, ilev, n, kg, c);
  
 DONE:
  if (nlev > 0) free(ilev);
  if (n > 0) free(kg);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTransitionMBPT(PyObject *self, PyObject *args) {
  int m, n, nlow, *low, nup, *up;
  char *fn;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("TransitionMBPT", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  n = PyTuple_Size(args);
  if (n == 2) {
    if (!(PyArg_ParseTuple(args, "ii", &m, &n))) return NULL;
    TransitionMBPT(m, n);
  } else if (n == 3) {
    if (!(PyArg_ParseTuple(args, "sOO", &fn, &p, &q))) return NULL;
    nlow = DecodeGroupArgs(p, &low, NULL);
    nup = DecodeGroupArgs(q, &up, NULL);
    TRTableMBPT(fn, nlow, low, nup, up);
    if (nlow > 0) free(low);
    if (nup > 0) free(up);
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PStructureMBPT(PyObject *self, PyObject *args) {
  PyObject *p, *q, *t, *r, *x, *y;
  int i, n, n1, *ng1, n2, *ng2, *s, nk, *nkm, kmax;
  int n3, *ng3, n4, *ng4;
  char *fn, *fn1, *gn, **fn2;
  double d, c, e, f;

  if (sfac_file) {
    SFACStatement("StructureMBPT", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  
  if (n == 1) {
    if (!(PyArg_ParseTuple(args, "O", &p))) return NULL;
    if (PyFloat_Check(p) || PyLong_Check(p)) {
      f = PyFloat_AsDouble(p);
      if (f < 0 || (f > 0 && f < 1)) {
	SetWarnMBPT(f, -1.0);
	Py_INCREF(Py_None);
	return Py_None;
      }
    }
    if (PyLong_Check(p)) {
      i = PyLong_AsLong(p);
      SetExtraMBPT(i);
    } else {
      n1 = IntFromList(p, &ng1);
      SetSymMBPT(n1, ng1);
      free(ng1);
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (n == 2) {
    if (!(PyArg_ParseTuple(args, "OO", &q, &p))) return NULL;
    if ((PyFloat_Check(q) || PyLong_Check(q)) &&
	(PyFloat_Check(p) || PyLong_Check(p))) {
      f = PyFloat_AsDouble(q);
      d = PyFloat_AsDouble(p);
      SetWarnMBPT(f, d);
      Py_INCREF(Py_None);
      return Py_None;
    }
    if (!PyUnicode_Check(q)) return NULL;
    gn = PyUnicode_AsString(q);
    if (PyLong_Check(p)) {
      n2 = PyLong_AsLong(p);
      n1 = n2;
    } else {
      n3 = IntFromList(p, &ng3);
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
    SetExcMBPT(n2, n1, gn);
    Py_INCREF(Py_None);
    return Py_None;
  }

  p = PyTuple_GET_ITEM(args, 0);
  if (PyLong_Check(p)) {
    if (n == 3 || n == 4 || n == 5 || n == 6) {
      d = -1.0;
      e = -1.0;
      f = -1.0;
      if (!(PyArg_ParseTuple(args, "iid|ddd",
			     &i, &n3, &c, &d, &e, &f))) return NULL;
      SetOptMBPT(i, n3, c, d, e, f);
      Py_INCREF(Py_None);
      return Py_None;
    }
  }
  /*
  if (n == 10) { 
    c = 0.0;
    if (!(PyArg_ParseTuple(args, "sddOiOOOOs",
			   &fn, &d, &c, &p, &kmax, &q, &r, &x, &y, &gn)))
      return NULL;
    
    n = DecodeGroupArgs(p, &s, NULL);
    if (n <= 0) return NULL;
    
    n1 = IntFromList(q, &ng1);
    n2 = IntFromList(r, &ng2);
    n3 = IntFromList(x, &ng3);
    n4 = IntFromList(y, &ng4);
    
    StructureMBPT0(fn, d, c, n, s, kmax, n1, ng1, n2, ng2, n3, ng3, n4, ng4, gn);
    
    free(s);
    if (n1 > 0) free(ng1);
    if (n2 > 0) free(ng2);
    if (n3 > 0) free(ng3);
    if (n4 > 0) free(ng4);

    Py_INCREF(Py_None);
    return Py_None;
  }
  */
  if (n == 5) {
    if (!(PyArg_ParseTuple(args, "ssOOi", &fn, &fn1, &q, &p, &n3)))
      return NULL;
    n = DecodeGroupArgs(p, &s, &n3);
    if (n <= 0) return NULL;
    n1 = PyList_Size(q);
    if (n1 <= 0) return NULL;
    fn2 = malloc(sizeof(char *)*n1);
    for (i = 0; i < n1; i++) {
      t = PyList_GetItem(q, i);
      if (!PyUnicode_Check(t)) return NULL;
      fn2[i] = PyUnicode_AsString(t);
    }
    StructureReadMBPT(fn, fn1, n1, fn2, n, s, n3);
    free(s);
    free(fn2);

    Py_INCREF(Py_None);
    return Py_None;
  } 

  int icp = 0;
  int icpf = 0;
  int ncp = 0;
  
  if (n == 7 || n == 9 || n == 10) {
    PyObject *fp, *fp0, *fp1;
    if (!(PyArg_ParseTuple(args, "sOOOOOi|iii",
			   &fn, &fp, &p, &t, &q, &r, &n3, &ncp, &icp, &icpf)))
      return NULL;
    
    char *hfn0, *hfn1;
    hfn0 = NULL;
    if (PyUnicode_Check(fp)) {
      hfn1 = PyUnicode_AsString(fp);
    } else {
      if (!PyList_Check(fp) && !PyTuple_Check(fp)) return NULL;
      int nf = PySequence_Length(fp);
      if (nf == 0) return NULL;
      fp0 = PySequence_GetItem(fp, 0);      
      if (!PyUnicode_Check(fp0)) return NULL;
      hfn1 = PyUnicode_AsString(fp0);
      Py_DECREF(fp0);
      if (nf > 1) {
	fp1 = PySequence_GetItem(fp, 1);
	if (!PyUnicode_Check(fp1)) return NULL;
	hfn0 = PyUnicode_AsString(fp1);
	Py_DECREF(fp1);
      }
    }
    if (PyList_Check(p)) {
      n = DecodeGroupArgs(p, &s, &n3);
      if (n <= 0) return NULL;
    } else {
      n = 0;
      s = NULL;
    }
    if (PyList_Check(q)) {
      n1 = IntFromList(q, &ng1);
    } else {
      n1 = PyLong_AsLong(q);
      ng1 = NULL;
    }
    if (PyList_Check(r)) {      
      n2 = IntFromList(r, &ng2);
    } else {
      n2 = PyLong_AsLong(r);
      ng2 = NULL;
    }
    if (PyList_Check(t)) {
      nk = IntFromList(t, &nkm);
    } else if (PyLong_Check(t)) {
      nk = PyLong_AsLong(t)+1;
      nkm = NULL;
    } else {
      return NULL;
    }
  
    StructureMBPT1(fn, hfn0, hfn1, n, s, nk, nkm, n1, ng1, n2, ng2, n3,
		   ncp, icp, icpf);
    free(s);
    //if (n1 > 0) free(ng1);
    //if (n2 > 0) free(ng2);
    if (nkm) free(nkm);

    Py_INCREF(Py_None);
    return Py_None;
  }

  return NULL;
}  

static PyObject *PPrepAngular(PyObject *self, PyObject *args) {
  PyObject *p, *q;  
  int nlow, nup, *low, *up;

  if (sfac_file) {
    SFACStatement("PrepAngular", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;
  
  q = NULL;
  if (!PyArg_ParseTuple(args, "O|O", &p, &q)) return NULL;
  nlow = SelectLevels(p, &low);
  if (nlow <= 0) return NULL;
  if (q) {
    nup = SelectLevels(q, &up);
    if (nup <= 0) {
      free(low);
      return NULL;
    }
  }
  PrepAngular(nlow, low, nup, up);

  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PSetFields(PyObject *self, PyObject *args) {
  int m;
  double b, e, a;

  if (sfac_file) {
    SFACStatement("SetFields", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  m = 0;
  if (!(PyArg_ParseTuple(args, "ddd|i", &b, &e, &a, &m))) 
    return NULL;

  SetFields(b, e, a, m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PStructureEB(PyObject *self, PyObject *args) {
  PyObject *p;
  int n, *ilev;
  char *fn;

  if (sfac_file) {
    SFACStatement("StructureEB", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!(PyArg_ParseTuple(args, "sO", &fn, &p))) return NULL;
  
  n = SelectLevels(p, &ilev);
  if (n == 0) return NULL;
  
  StructureEB(fn, n, ilev);
  free(ilev);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTransitionTableEB(PyObject *self, PyObject *args) {
  char *s;
  int m, nlow, nup, *low, *up;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("TransitionTableEB", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;
  m = -1;

  if (!PyArg_ParseTuple(args, "sOO|i", &s, &p, &q, &m)) 
    return NULL;
  nlow = SelectLevels(p, &low);
  if (nlow <= 0) {
    printf("cannot determine levels in lower\n");
    return NULL;
  }
  nup = SelectLevels(q, &up);
  if (nup <= 0) {
    printf("cannot determine levels in upper\n");
    return NULL;
  }
  
  SaveTransitionEB(nlow, low, nup, up, s, m);
  free(low);
  free(up);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PPolarizeCoeff(PyObject *self, PyObject *args) {
  char *fn1, *fn2;
  int i0, i1;

  if (sfac_file) {
    SFACStatement("PolarizeCoeff", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  i0 = -1;
  i1 = -1;

  if (!PyArg_ParseTuple(args, "ss|ii", &fn1, &fn2, &i0, &i1)) return NULL;

  PolarizeCoeff(fn1, fn2, i0, i1);

  Py_INCREF(Py_None);
  return Py_None;
}  
  
static PyObject *PElectronDensity(PyObject *self, PyObject *args) {
  char *ofn;
  int n, *ilev, t;
  PyObject *p;
  
  if (sfac_file) {
    SFACStatement("ElectronDensity", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  t = 1;
  if (!PyArg_ParseTuple(args, "sO|i", &ofn, &p, &t)) return NULL;
  n = SelectLevels(p, &ilev);
  if (n > 0) {
    ElectronDensity(ofn, n, ilev, t);
    free(ilev);
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PExpectationValue(PyObject *self, PyObject *args) {
  char *ifn, *ofn;
  int n, *ilev, t;
  double a;
  PyObject *p;
  
  if (sfac_file) {
    SFACStatement("ExpectationValue", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  t = 1;
  if (!PyArg_ParseTuple(args, "ssO|di", &ifn, &ofn, &p, &a, &t)) return NULL;
  n = SelectLevels(p, &ilev);
  if (n > 0) {
    ExpectationValue(ifn, ofn, n, ilev, a, t);
    free(ilev);
  }
  
  Py_INCREF(Py_None);
  return Py_None;
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
  m = 0;

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
    if (nlow <= 0) {
      printf("cannot determine levels in lower\n");
      return NULL;
    }
    nup = SelectLevels(q, &up);
    if (nup <= 0) {
      printf("cannot determine levels in upper\n");
      return NULL;
    }
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
  int m, k;

  if (sfac_file) {
    SFACStatement("BasisTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = 0;
  k = -1;
  if (!PyArg_ParseTuple(args, "s|ii", &s, &m, &k)) return NULL;
  GetBasisTable(s, m, k);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PCETableEB(PyObject *self, PyObject *args) {
  char *s;
  int nlow, nup, *low, *up, m;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("CETableEB", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = 0;
  if (!PyArg_ParseTuple(args, "sOO|i", &s, &p, &q, &m)) return NULL;
  nlow = SelectLevels(p, &low);
  if (nlow <= 0) return NULL;
  nup = SelectLevels(q, &up);
  if (nup <= 0) return NULL;
  if (m == 0) {
    SaveExcitationEB(nlow, low, nup, up, s);
  } else {
    SaveExcitationEBD(nlow, low, nup, up, s);
  }

  free(low);
  free(up);  

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
    if (PyLong_Check(p)) {
      ng = PyLong_AsLong(p);
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
  
static  PyObject *PSetCEBorn(PyObject *self, PyObject *args) {
  double eb, x, x1, x0;
  
  if (sfac_file) {
    SFACStatement("SetCEBorn", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  x0 = XBORN0;
  x1 = XBORN1;
  x = XBORN;
  if (!PyArg_ParseTuple(args, "d|ddd", &eb, &x, &x1, &x0)) return NULL;

  SetCEBorn(eb, x, x1, x0);
  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetCIBorn(PyObject *self, PyObject *args) {
  int x;
  
  if (sfac_file) {
    SFACStatement("SetCIBorn", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &x)) return NULL;

  SetCIBorn(x);
  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetBornFormFactor(PyObject *self, PyObject *args) {
  double te;
  char *fn;

  if (sfac_file) {
    SFACStatement("SetBornFormFactor", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  fn = NULL;
  if (!PyArg_ParseTuple(args, "d|s", &te, &fn)) return NULL;

  SetBornFormFactor(te, fn);
  
  Py_INCREF(Py_None);
  return Py_None;
}
      
static  PyObject *PSetBornMass(PyObject *self, PyObject *args) {
  double m;

  if (sfac_file) {
    SFACStatement("SetBornMass", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "d", &m)) return NULL;

  SetBornMass(m);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetCEQkMode(PyObject *self, PyObject *args) {
  PyObject *p;
  int m;
  double tol;

  if (sfac_file) {
    SFACStatement("SetCEQkMode", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = QK_DEFAULT;
  tol = -1;
  if (!PyArg_ParseTuple(args, "|Od", &p, &tol)) return NULL;
  if (PyUnicode_Check(p)) {
    p = PyDict_GetItem(QKMODE, p);
  } 
  if (PyLong_Check(p)) {
    m = PyLong_AsLong(p);
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
    if (PyLong_Check(p)) {
      ng = PyLong_AsLong(p);
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
   
static PyObject *PSetAngleGrid(PyObject *self, PyObject *args) {
  int n, ng, m, i, err;
  double xg[MAXNPHI+MAXNTHETA];
  double emin, emax;
  PyObject *p, *pi;

  if (sfac_file) {
    SFACStatement("SetAngleGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n == 2) {
    if (!PyArg_ParseTuple(args, "iO", &m, &p)) return NULL;
    if (PyLong_Check(p)) {
      ng = PyLong_AsLong(p);
      if (m == 0) {
	emin = 0.0;
	emax = PI;
      } else {
	emin = 0.0;
	emax = TWO_PI;
      }
      err = SetAngleGrid(m, ng, emin, emax);
    } else if (PyList_Check(p) || PyTuple_Check(p)) {
      ng = PySequence_Length(p);      
      for (i = 0; i < ng; i++) {
	pi = PySequence_GetItem(p, i);
	xg[i] = PyFloat_AsDouble(pi)*PI/180.0;
	Py_DECREF(pi);
      }
      err = SetAngleGridDetail(m, ng, xg);
    } else {
      return NULL;
    }
  } else if (n == 4) {
    if (!PyArg_ParseTuple(args, "iidd", &m, &ng, &emin, &emax)) 
      return NULL;
    emin *= PI/180.0;
    emax *= PI/180.0;
    err = SetAngleGrid(m, ng, emin, emax);
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
    if (PyLong_Check(p)) {
      ng = PyLong_AsLong(p);
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

  if (sfac_file) {
    SFACStatement("SetCEPWOptions", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  qr = EXCLQR;
  max = EXCLMAX;
  kl_cb = EXCLCB;
  tol = EXCTOL;

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

  if (sfac_file) {
    SFACStatement("SetCEPWGrid", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

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
      m[i] = PyLong_AsLong(PyList_GetItem(p, i));
      step[i] = PyLong_AsLong(PyList_GetItem(q, i));
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
  
  if (sfac_file) {
    SFACStatement("SetRecQkMode", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = QK_DEFAULT;
  tol = -1;
  if (!PyArg_ParseTuple(args, "|Od", &p, &tol)) return NULL;
  if (PyUnicode_Check(p)) {
    p = PyDict_GetItem(QKMODE, p);
  }
  if (PyLong_Check(p)) {
    m = PyLong_AsLong(p);
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
   
  if (sfac_file) {
    SFACStatement("SetRecPWOptions", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
 
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

  if (sfac_file) {
    SFACStatement("SetRecPWLimits", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  if (!PyArg_ParseTuple(args, "ii", &m1, &m2)) return NULL;
  SetRecPWLimits(m1, m2);
  Py_INCREF(Py_None);
  return Py_None;
}
  

static  PyObject *PSetRecSpectator(PyObject *self, PyObject *args) {
  int n_spec, n_frozen, n_max;

  if (sfac_file) {
    SFACStatement("SetRecSpectator", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

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
  ng = DecodeGroupArgs(gargs, &kg, NULL);
  if (ng <= 0) return NULL;

  if (RecStates(n, ng, kg, fn) < 0) {
    onError("RecStates error");
    free(kg);
    return NULL;
  }

  free(kg);
  
  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PRRMultipole(PyObject *self, PyObject *args) { 
  int nlow, *low, nup, *up;
  int m;
  char *s;
  PyObject *p, *q;
  
  if (sfac_file) {
    SFACStatement("RRMultipole", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  m = -1;
  if (!PyArg_ParseTuple(args, "sOO|i", &s, &p, &q, &m)) {
    printf("Unrecognized parameters in RRMultipole\n");
    return NULL;
  }
  nlow = SelectLevels(p, &low);
  nup = SelectLevels(q, &up);
  SaveRRMultipole(nlow, low, nup, up, s, m);
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

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

static PyObject *PAsymmetry(PyObject *self, PyObject *args) {
  char *s, *fn;
  int mx;

  if (sfac_file) {
    SFACStatement("Asymmetry", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  mx = 1;
  if (!PyArg_ParseTuple(args, "ss|i", &fn, &s, &mx)) return NULL;
  
  SaveAsymmetry(fn, s, mx);
  
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
    if (PyLong_Check(p)) {
      ng = PyLong_AsLong(p);
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
    if (PyLong_Check(p)) {
      ng = PyLong_AsLong(p);
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
    if (PyLong_Check(p)) {
      ng = PyLong_AsLong(p);
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
  int nlow, *low, nup, *up;
  double emin;
  char *s;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("AITable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  emin = 0;
  if (!PyArg_ParseTuple(args, "sOO|d", &s, &p, &q, &emin)) return NULL;
  nlow = SelectLevels(p, &low);
  nup = SelectLevels(q, &up);
  SaveAI(nlow, low, nup, up, s, emin, 0);
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PAITableMSub(PyObject *self, PyObject *args) { 
  int nlow, *low, nup, *up;
  double emin;
  char *s;
  PyObject *p, *q;

  if (sfac_file) {
    SFACStatement("AITableMSub", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  emin = 0;
  if (!PyArg_ParseTuple(args, "sOO|d", &s, &p, &q, &emin)) return NULL;
  nlow = SelectLevels(p, &low);
  nup = SelectLevels(q, &up);
  SaveAI(nlow, low, nup, up, s, emin, 1);
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PReportMultiStats(PyObject *self, PyObject *args) {  
  if (sfac_file) {
    SFACStatement("ReportMultiStats", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  ReportMultiStats();
  Py_INCREF(Py_None);
  return Py_None;
}
    
static PyObject *PTestMyArray(PyObject *self, PyObject *args) {
  ARRAY a;
  double d;
  double *b;
  MULTI ma;
  int k[3] = {101, 2550, 333};
  int block[3] = {2, 2, 5};
  int i, j, m;

  if (sfac_file) {
    SFACStatement("TestMyArray", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  ArrayInit(&a, sizeof(double), 100);
  d = 0.1;
  m = 10000;

  for (i = 0; i < m; i++) {
    ArraySet(&a, i, &d, InitDoubleData);
  }

  b = (double *) ArrayGet(&a, 100);
  b = (double *) ArrayGet(&a, 200);

  printf("array set\n"); 
  ArrayFree(&a, 0);
  printf("array freed\n");

  MultiInit(&ma, sizeof(double), 3, block, "test_ma");
  for (i = 0; i < 500; i++) {
    for (j = 100; j < 4000; j++) {
      k[0] = i;
      k[1] = j;
      k[2] = 20;
      b = (double *) MultiSet(&ma, k, NULL, NULL, InitDoubleData, NULL);
      *b = 0.2;
      b = (double *) MultiGet(&ma, k, NULL);
    }
  }
  printf("set\n"); 
  MultiFreeData(&ma, NULL);
  printf("freed\n");
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
    if (PyLong_Check(p)) {
      ng = PyLong_AsLong(p);
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
  if (PyUnicode_Check(p)) {
    p = PyDict_GetItem(QKMODE, p);
  }
  if (PyLong_Check(p)) {
    m = PyLong_AsLong(p);
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
    if (PyLong_Check(p)) {
      ng = PyLong_AsLong(p);
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
    if (PyLong_Check(p)) {
      ng = PyLong_AsLong(p);
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

static  PyObject *PSetCIPWOptions(PyObject *self, PyObject *args) {
  int qr, max, max_1, kl_cb;
  double tol;

  if (sfac_file) {
    SFACStatement("SetCIPWOptions", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  qr = IONLQR;
  max = IONLMAX;
  max_1 = IONLEJEC;
  kl_cb = IONLCB;
  tol = IONTOL;
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
      m[i] = PyLong_AsLong(PyList_GetItem(p, i));
      step[i] = PyLong_AsLong(PyList_GetItem(q, i));
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
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PCITableMSub(PyObject *self, PyObject *args) { 
  int nlow, *low, nup, *up;
  char *s;
  PyObject *p, *q;
 
  if (sfac_file) {
    SFACStatement("CITableMSub", args, NULL);
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
  SaveIonizationMSub(nlow, low, nup, up, s);
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTestIntegrate(PyObject *self, PyObject *args) { 

  if (sfac_file) {
    SFACStatement("TestIntegrate", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  TestIntegrate();
  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject *PCoulombBethe(PyObject *self, PyObject *args) { 
  char *s;
  double z, te, e1;
  
  if (sfac_file) {
    SFACStatement("CoulombBethe", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "sddd", &s, &z, &te, &e1)) return NULL;
  CoulombBethe(s, z, te, e1);

  Py_INCREF(Py_None);
  return Py_None;  
}

static PyObject *PAdjustEnergy(PyObject *self, PyObject *args) {
  char *s, *efn0, *efn1, *afn0, *afn1;  
  PyObject *p, *q, *ip, *iq;
  int i, n, k, nlevs, *ilevs;
  double e, *elevs;
  FILE *f;
  
  
  if (sfac_file) {
    SFACStatement("AdjustEnergy", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n == 5) {
    if (!PyArg_ParseTuple(args, "sssss", &s, &efn0, &efn1, &afn0, &afn1)) {
      return NULL;
    }
    f = fopen(s, "r");    
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
  } else {
    if (!PyArg_ParseTuple(args, "OOssss", 
			  &p, &q, &efn0, &efn1, &afn0, &afn1)) 
      return NULL;
    if (!PyList_Check(p) || !PyList_Check(q)) {
      return NULL;
    }
    nlevs = PyList_Size(p);
    if (PyList_Size(q) != nlevs) {
      printf("The energy list length not equal the index list length\n");
      return NULL;
    }
    ilevs = (int *) malloc(sizeof(int)*nlevs);
    elevs = (double *) malloc(sizeof(double)*nlevs);
    for (i = 0; i < n; i++) {
      ip = PyList_GetItem(p, i);
      iq = PyList_GetItem(q, i);
      k = PyLong_AsLong(ip);
      e = PyFloat_AsDouble(iq);
      e /= HARTREE_EV;
      ilevs[i] = k;
      elevs[i] = e;
    }
  }

  AdjustEnergy(nlevs, ilevs, elevs, efn0, efn1, afn0, afn1);

  if (nlevs > 0) {
    free(ilevs);
    free(elevs);
  }
  Py_INCREF(Py_None);
  return Py_None;
}

  
static PyObject *PCorrectEnergy(PyObject *self, PyObject *args) {
  char *s;
  PyObject *p, *q, *ip, *iq;
  int n, k, kref;
  double e;
  int i, nmin;
  FILE *f;

  if (sfac_file) {
    SFACStatement("CorrectEnergy", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  kref = 0;
  n = PyTuple_Size(args);
  if (n == 2) {
    if (!PyArg_ParseTuple(args, "si", &s, &nmin)) {
      return NULL;
    }
    f = fopen(s, "r");
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
  } else if (n == 3) {
    if (!PyArg_ParseTuple(args, "OOi", &p, &q, &nmin)) return NULL;
    if (!PyList_Check(p) || !PyList_Check(q)) {
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
      k = PyLong_AsLong(ip);
      e = PyFloat_AsDouble(iq);
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
    return NULL;
  }

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

  v = 1;
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

static PyObject *PReinitConfig(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("ReinitConfig", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  ReinitConfig(m);
  _closed_shells[0] = '\0';

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PReinitRecouple(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("ReinitRecouple", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  ReinitRecouple(m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PReinitRadial(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("ReinitRadial", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  ReinitRadial(m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PReinitDBase(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("ReinitDBase", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  ReinitDBase(m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PReinitStructure(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("ReinitStructure", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  ReinitStructure(m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PReinitExcitation(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("ReinitExcitation", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  ReinitExcitation(m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PReinitRecombination(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("ReinitRecombination", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  ReinitRecombination(m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PReinitIonization(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("ReinitIonization", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  ReinitIonization(m);

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

  if (!keywds || PyDict_Size(keywds) == 0) {
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
      if (!PyLong_Check(q)) return NULL;
      m_config = PyLong_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "recouple");
    if (q) {
      if (!PyLong_Check(q)) return NULL;
      m_recouple = PyLong_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "dbase");
    if (q) {
      if (!PyLong_Check(q)) return NULL;
      m_dbase = PyLong_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "structure");
    if (q) {
      if (!PyLong_Check(q)) return NULL;
      m_structure = PyLong_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "excitation");
    if (q) {
      if (!PyLong_Check(q)) return NULL;
      m_excitation = PyLong_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "radial");
    if (q) {
      if (!PyLong_Check(q)) return NULL;
      m_radial = PyLong_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "recombination");
    if (q) {
      if (!PyLong_Check(q)) return NULL;
      m_recombination = PyLong_AsLong(q);
    }
    q = PyDict_GetItemString(keywds, "ionization");
    if (q) {
      if (!PyLong_Check(q)) return NULL;
      m_ionization = PyLong_AsLong(q);
    }
  }

  ReinitFac(m_config, m_recouple, m_radial, m_dbase,
	    m_structure, m_excitation, m_recombination, m_ionization);
  if (m_config == 0) _closed_shells[0] = '\0';

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
    s = PyUnicode_AsString(q);
    printf("%s", s);
    if (i != n-1) {
      printf(" ");
    }
    Py_XDECREF(q);
  }
  
  if (n > 0) printf("\n");

  fflush(stdout);

  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PConfigEnergy(PyObject *self, PyObject *args) {
  int m, mr, n, i;
  PyObject *p;
  int ng, *kg;

  if (sfac_file) {
    SFACStatement("ConfigEnergy", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n == 0) return NULL;

  p = PyTuple_GetItem(args, 0);
  m = PyLong_AsLong(p);
  if (n == 1 || m != 0) {
    ConfigEnergy(m, 0, 0, NULL);
  } else {
    p = PyTuple_GetItem(args, 1);
    mr = PyLong_AsLong(p);
    if (n == 2) {
      ConfigEnergy(m, mr, 0, NULL);
    } else {
      for (i = 1; i < n; i++) {
	p = PyTuple_GetItem(args, i);
	if (!PyList_Check(p)) return NULL;
	ng = DecodeGroupArgs(p, &kg, NULL);
	if (ng < 0) return NULL;
	ConfigEnergy(m, mr, ng, kg);
	if (ng > 0) free(kg);
      }
    }
  }

  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PTRRateH(PyObject *self, PyObject *args) {
  int os, n0, n1;
  int kl0, kl1;
  double z, r;

  os = 0;
  if (!PyArg_ParseTuple(args, "diiii|i", &z, &n0, &kl0, &n1, &kl1, &os)) 
    return NULL;
  if (n1 > 512) {
    printf("maximum NU is 512\n");
    return NULL;
  }
  if (kl0 != kl1+1 && kl0 != kl1-1) {
    r = 0.0;
  } else {
    r = TRRateHydrogenic(z, n0, kl0, n1, kl1, os);
  }
  return Py_BuildValue("d", r);
}

static PyObject *PPICrossH(PyObject *self, PyObject *args) {
  int os, n0, kl0;
  double e, z, r;

  os = 0;
  if (!PyArg_ParseTuple(args, "ddii|i", &z, &e, &n0, &kl0, &os)) 
    return NULL;
  if (n0 > 256) {
    printf("maximum NU is 256\n");
    return NULL;
  }

  e = e/HARTREE_EV;
  r = PICrossH(z, n0, kl0, e, os);

  return Py_BuildValue("d", r);
}

static PyObject *PRRCrossH(PyObject *self, PyObject *args) {
  int n0, kl0;
  double e, z, r;

  if (!PyArg_ParseTuple(args, "ddii", &z, &e, &n0, &kl0)) 
    return NULL;
  if (n0 > 256) {
    printf("maximum NU is 256\n");
    return NULL;
  }

  e = e/HARTREE_EV;
  r = RRCrossH(z, n0, kl0, e);

  return Py_BuildValue("d", r);
}

static PyObject *PRRCrossHn(PyObject *self, PyObject *args) {
  int n0;
  double e, z, r;

  if (!PyArg_ParseTuple(args, "ddi", &z, &e, &n0)) 
    return NULL;

  e = e/HARTREE_EV;
  r = AREA_AU20*RRCrossHn(z, e, n0);

  return Py_BuildValue("d", r);
}

static PyObject *PTotalRRCross(PyObject *self, PyObject *args) {
  PyObject *p, *q;
  int i, negy, ilev;
  double *egy;
  char *ifn, *ofn;
  int n0, n1, nmax, imin, imax;

  n0 = 0;
  n1 = 0;
  nmax =0;  
  imin = -1;
  imax = -1;
  if (!PyArg_ParseTuple(args, "ssiO|iiiii", 
			&ifn, &ofn, &ilev, &p, &n0, &n1, 
			&nmax, &imin, &imax))
    return NULL;
  
  if (!PyList_Check(p) && !PyTuple_Check(p)) {
    printf("Energy List must be a sequence\n");
    return NULL;
  }
  
  negy = PySequence_Length(p);
  egy = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    q = PySequence_GetItem(p, i);
    egy[i] = PyFloat_AsDouble(q);
    Py_DECREF(q);
  }
  
  TotalRRCross(ifn, ofn, ilev, negy, egy, n0, n1, nmax, imin, imax);

  free(egy);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTotalPICross(PyObject *self, PyObject *args) {
  PyObject *p, *q;
  int i, negy, ilev;
  double *egy;
  char *ifn, *ofn;
  int imin, imax;

  imin = -1;
  imax = -1;
  if (!PyArg_ParseTuple(args, "ssiO|ii", 
			&ifn, &ofn, &ilev, &p,
			&imin, &imax))
    return NULL;
  
  if (!PyList_Check(p) && !PyTuple_Check(p)) {
    printf("Energy List must be a sequence\n");
    return NULL;
  }
  
  negy = PySequence_Length(p);
  egy = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    q = PySequence_GetItem(p, i);
    egy[i] = PyFloat_AsDouble(q);
    Py_DECREF(q);
  }
  
  TotalPICross(ifn, ofn, ilev, negy, egy, imin, imax);

  free(egy);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTotalCICross(PyObject *self, PyObject *args) {
  PyObject *p, *q;
  int i, negy, ilev;
  double *egy;
  char *ifn, *ofn;
  int imin, imax;

  imin = -1;
  imax = -1;
  if (!PyArg_ParseTuple(args, "ssiO|ii", 
			&ifn, &ofn, &ilev, &p,
			&imin, &imax))
    return NULL;
  
  if (!PyList_Check(p) && !PyTuple_Check(p)) {
    printf("Energy List must be a sequence\n");
    return NULL;
  }
  
  negy = PySequence_Length(p);
  egy = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    q = PySequence_GetItem(p, i);
    egy[i] = PyFloat_AsDouble(q);
    Py_DECREF(q);
  }
  
  TotalCICross(ifn, ofn, ilev, negy, egy, imin, imax);

  free(egy);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PY5N(PyObject *self, PyObject *args) {
  double lambda, xi, eta0, x0;
  double y5, y5i, y5p, y5pi;
  int ierr;

  if (!PyArg_ParseTuple(args, "ddd", &eta0, &lambda, &x0))
    return NULL;
  
  xi = 0.0;
  Y5N(lambda, xi, eta0, x0, &y5, &y5i, &y5p, &y5pi, &ierr);
  
  return Py_BuildValue("(ddddi)", y5, y5p, y5i, y5pi, ierr);
}

static PyObject *PDiracCoulomb(PyObject *self, PyObject *args) {
  double z, e, r;
  double p, q, u, v;
  int k, ierr;

  ierr = 1;
  if (!PyArg_ParseTuple(args, "ddid|i", &z, &e, &k, &r, &ierr)) return NULL;
  if (ierr < 0) {
    e = RadialDiracCoulomb(1, &p, &q, &r, z, (int)e, k);
    return Py_BuildValue("(ddd)", p, q, e*HARTREE_EV);
  } else {
    e /= HARTREE_EV;
    DCOUL(z, e, k, r, &p, &q, &u, &v, &ierr);
    return Py_BuildValue("(ddddi)", p, q, u, v, ierr);
  }
}
  
static PyObject *PCoulombPhase(PyObject *self, PyObject *args) {
  double z, e, p;
  int k;
  
  if (!PyArg_ParseTuple(args, "ddi", &z, &e, &k)) return NULL;

  e /= HARTREE_EV;
  p = CoulombPhaseShift(z, e, k);

  while (p < 0) p += TWO_PI;
  p = p - (int)(p/TWO_PI);

  return Py_BuildValue("d", p);
}

static PyObject *PSetConfigEnergyMode(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetConfigEnergyMode", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetConfigEnergyMode(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetOptimizeMaxIter(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetOptimizeMaxIter", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetOptimizeMaxIter(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetOptimizeStabilizer(PyObject *self, PyObject *args) {
  double m;

  if (sfac_file) {
    SFACStatement("SetOptimizeStabilizer", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &m)) return NULL;
  SetOptimizeStabilizer(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetOptimizePrint(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetOptimizePrint", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetOptimizePrint(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetOptimizeTolerance(PyObject *self, PyObject *args) {
  double m;

  if (sfac_file) {
    SFACStatement("SetOptimizeTolerance", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &m)) return NULL;
  SetOptimizeTolerance(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCELQR(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetCELQR", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCELQR(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCELMax(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetCELMax", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCELMax(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCELCB(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetCELCB", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCELCB(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCETol(PyObject *self, PyObject *args) {
  double m;

  if (sfac_file) {
    SFACStatement("SetCETol", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &m)) return NULL;
  SetCETol(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCILQR(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetCILQR", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCILQR(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCILMax(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetCILMax", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCILMax(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCILMaxEject(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetCILMaxEject", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCILMaxEject(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCILCB(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetCILCB", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetCILCB(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCITol(PyObject *self, PyObject *args) {
  double m;

  if (sfac_file) {
    SFACStatement("SetCITol", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &m)) return NULL;
  SetCITol(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetTransitionMode(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetTransitionMode", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetTransitionMode(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetTransitionGauge(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetTransitionGauge", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetTransitionGauge(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetTransitionMaxE(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetTransitionMaxE", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetTransitionMaxE(m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetTransitionMaxM(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("SetTransitionMaxM", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  SetTransitionMaxM(m);
  Py_INCREF(Py_None);
  return Py_None;
}
      
static PyObject *PLevelInfor(PyObject *self, PyObject *args) { 
  char *fn;
  int i, k;
  EN_RECORD r;
  
  if (!PyArg_ParseTuple(args, "si", &fn, &i)) return NULL;

  r.ncomplex[0] = '\0';
  r.sname[0] = '\0';
  r.name[0] = '\0';
  k = LevelInfor(fn, i, &r);
  
  if (k < 0) {
    if (i >= 0) {
      return Py_BuildValue("()");
    } else {
      return Py_BuildValue("i", k);
    }
  }
  if (i >= 0) {
    return Py_BuildValue("(diisss)", r.energy, r.p, r.j,
			 r.ncomplex, r.sname, r.name);
  } else {
    return Py_BuildValue("i", k);
  }
}

static PyObject *PInterpCross(PyObject *self, PyObject *args) { 
  PyObject *p;
  int i, negy, i0, i1, mp;
  double *egy;
  char *ifn, *ofn;
  
  if (sfac_file) {
    SFACStatement("InterpCross", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  mp = 1;
  if (!PyArg_ParseTuple(args, "ssiiO|i", &ifn, &ofn, &i0, &i1, &p, &mp))
    return NULL;

  negy = DoubleFromList(p, &egy);
  if (negy > 0) {
    InterpCross(ifn, ofn, i0, i1, negy, egy, mp);
    free(egy);
  }
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PMaxwellRate(PyObject *self, PyObject *args) { 
  PyObject *p;
  int i, nt, i0, i1;
  double *temp;
  char *ifn, *ofn;
  
  if (sfac_file) {
    SFACStatement("MaxwellRate", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "ssiiO", &ifn, &ofn, &i0, &i1, &p))
    return NULL;

  nt = DoubleFromList(p, &temp);
  if (nt > 0) {
    MaxwellRate(ifn, ofn, i0, i1, nt, temp);
    free(temp);
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTRBranch(PyObject *self, PyObject *args) { 
  int i, j, nt;
  char *fn;
  double te, pa, ta;

  if (!PyArg_ParseTuple(args, "sii", &fn, &i, &j)) return NULL;
  nt = TRBranch(fn, i, j, &te, &pa, &ta);
  if (nt < 0) return NULL;

  return Py_BuildValue("(dddi)", te, pa, ta, nt);
}
  
static PyObject *PAIBranch(PyObject *self, PyObject *args) { 
  int i, j, nt;
  char *fn;
  double te, pa, ta;

  if (!PyArg_ParseTuple(args, "sii", &fn, &i, &j)) return NULL;
  nt = AIBranch(fn, i, j, &te, &pa, &ta);
  if ( nt < 0) return NULL;

  return Py_BuildValue("(dddi)", te, pa, ta, nt);
}

static PyObject *PRadialOverlaps(PyObject *self, PyObject *args) {
  char *fn;
  int kappa;

  if (sfac_file) {
    SFACStatement("RadialOverlaps", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  kappa = -1;
  if (!PyArg_ParseTuple(args, "s|i", &fn, &kappa)) return NULL;

  RadialOverlaps(fn, kappa);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PRMatrixExpansion(PyObject *self, PyObject *args) {
  int m;
  double d, a, r;

  if (sfac_file) {
    SFACStatement("RMatrixExpansion", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  d = 1E-3;
  a = 1E-4;
  r = 0.0;
  if (!PyArg_ParseTuple(args, "i|ddd", &m, &r, &d, &a)) return NULL;
  RMatrixExpansion(m, d, a, r);

  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PRMatrixNBatch(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("RMatrixNBatch", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  RMatrixNBatch(m);

  Py_INCREF(Py_None);
  return Py_None;
}  
  
static PyObject *PRMatrixNMultipoles(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("RMatrixNMultipoles", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  RMatrixNMultipoles(m);

  Py_INCREF(Py_None);
  return Py_None;
}  
  
static PyObject *PRMatrixBoundary(PyObject *self, PyObject *args) {
  double r0, r1, b;

  if (sfac_file) {
    SFACStatement("RMatrixBoundary", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "ddd", &r0, &r1, &b)) return NULL;
  RMatrixBoundary(r0, r1, b);
  
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PRMatrixBasis(PyObject *self, PyObject *args) {
  int kmax, nb;
  char *fn;

  if (sfac_file) {
    SFACStatement("RMatrixBasis", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "sii", &fn, &kmax, &nb)) return NULL;
  
  RMatrixBasis(fn, kmax, nb);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PRMatrixTargets(PyObject *self, PyObject *args) {
  PyObject *p, *q;
  int nt, *kt, nc, *kc;
  
  if (sfac_file) {
    SFACStatement("RMatrixTargets", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  nt = 0;
  nc = 0;
  q = NULL;
  kc = NULL;
  kt = NULL;
  if (!PyArg_ParseTuple(args, "O|O", &p, &q)) return NULL;
  
  if (PyTuple_Check(p) || PyList_Check(p)) {
    nt = DecodeGroupArgs(p, &kt, NULL);
    if (nt == 0) return NULL;
  }
  if (q) {
    if (PyTuple_Check(q) || PyList_Check(q)) {
      nc = DecodeGroupArgs(q, &kc, NULL);
    } else {
      nc = 0;
      kc = NULL;
    }
  }

  RMatrixTargets(nt, kt, nc, kc);

  if (nt > 0) free(kt);
  if (nc > 0) free(kc);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PRMatrixSurface(PyObject *self, PyObject *args) {
  char *fn;
  
  if (sfac_file) {
    SFACStatement("RMatrixSurface", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  RMatrixSurface(fn);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PRMatrixConvert(PyObject *self, PyObject *args) {
  char *ifn, *ofn;
  int m;

  if (sfac_file) {
    SFACStatement("RMatrixConvert", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "ssi", &ifn, &ofn, &m)) return NULL;
  RMatrixConvert(ifn, ofn, m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PRMatrixFMode(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("RMatrixFMode", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  RMatrixFMode(m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTestRMatrix(PyObject *self, PyObject *args) {
  char *fn1, *fn2, *fn3;
  double e;
  int m;
  
  if (sfac_file) {
    SFACStatement("TestRMatrix", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "disss", &e, &m, &fn1, &fn2, &fn3)) return NULL;
  
  TestRMatrix(e, m, fn1, fn2, fn3);
  
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PSetSlaterCut(PyObject *self, PyObject *args) {
  int k0, k1;
  
  if (sfac_file) {
    SFACStatement("SetSlaterCut", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "ii", &k0, &k1)) return NULL;
  
  SetSlaterCut(k0, k1);

  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PPropogateDirection(PyObject *self, PyObject *args) {
  int m;

  if (sfac_file) {
    SFACStatement("PropogateDirection", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  
  PropogateDirection(m);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PRMatrixRefine(PyObject *self, PyObject *args) {
  int n, m;
  double r;
  
  if (sfac_file) {
    SFACStatement("RMatrixRefine", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = -1;
  r = -1.0;
  if (!PyArg_ParseTuple(args, "i|id", &n, &m, &r)) return NULL;
  RMatrixRefine(n, m, r);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PRMatrixCE(PyObject *self, PyObject *args) {
  PyObject *p, *q;
  char **f1, **f2, *fn;
  double emin, emax, de;
  int i, np, m, mb;

  if (sfac_file) {
    SFACStatement("RMatrixCE", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = 0;
  mb = 1;
  if (!PyArg_ParseTuple(args, "sOOddd|ii", 
			&fn, &p, &q, &emin, &emax, &de, &m, &mb)) return NULL;
  
  if (!PyList_Check(p)) return NULL;
  if (!PyList_Check(q)) return NULL;
  np = PyList_Size(p);
  if (PyList_Size(q) != np) return NULL;
  f1 = malloc(sizeof(char *)*np);
  f2 = malloc(sizeof(char *)*np);
  for (i = 0; i < np; i++) {
    f1[i] = PyUnicode_AsString(PyList_GetItem(p, i));
    f2[i] = PyUnicode_AsString(PyList_GetItem(q, i));
  }
  RMatrixCE(fn, np, f1, f2, emin, emax, de, m, mb);
  
  free(f1);
  free(f2);

  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PSetCEPWFile(PyObject *self, PyObject *args) {
  char *fn;

  if (sfac_file) {
    SFACStatement("SetCEPWFile", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  
  SetCEPWFile(fn);

  Py_INCREF(Py_None);
  return Py_None;
}    

static PyObject *PAppendTable(PyObject *self, PyObject *args) {
  char *fn;  
  
  if (sfac_file) {
    SFACStatement("AppendTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  AppendTable(fn);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PJoinTable(PyObject *self, PyObject *args) {
  char *fn, *fn1, *fn2;  
  
  if (sfac_file) {
    SFACStatement("JoinTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  if (!PyArg_ParseTuple(args, "sss", &fn1, &fn2, &fn)) return NULL;
  JoinTable(fn1, fn2, fn);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PModifyTable(PyObject *self, PyObject *args) {
  char *fn, *fn1, *fn2, *fnm; 
  
  if (sfac_file) {
    SFACStatement("ModifyTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  fnm = NULL;
  if (!PyArg_ParseTuple(args, "sss|s", &fn, &fn1, &fn2, &fnm)) return NULL;
  ModifyTable(fn, fn1, fn2, fnm);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PLimitArray(PyObject *self, PyObject *args) {
  int m;
  double n;
  
  if (sfac_file) {
    SFACStatement("LimitArray", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  if (!PyArg_ParseTuple(args, "id", &m, &n)) return NULL;
  LimitArrayRadial(m, n);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PWignerDMatrix(PyObject *self, PyObject *args) {
  int j2, m2, n2;
  double a;

  if (!PyArg_ParseTuple(args, "diii", &a, &j2, &m2, &n2)) return NULL;

  a *= PI/180.0;
  a = WignerDMatrix(a, j2, m2, n2);
  
  return Py_BuildValue("d", a);
}

static PyObject *PCoulMultip(PyObject *self, PyObject *args) {
  double z, te, e1;
  int k, q0, q1, m, ierr;
  char *fn;

  m = 1;
  if (!PyArg_ParseTuple(args, "sdddiii|i", 
			&fn, &z, &te, &e1, &k, &q0, &q1, &m))
    return NULL;
  
  ierr = CoulombMultip(fn, z, te, e1, k, q0, q1, m);
  if (ierr) {
    return NULL;
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSlaterCoeff(PyObject *self, PyObject *args) {
  char *fn, *q1, *q2;
  PyObject *p;  
  int nlev, *ilev, na, nb, i, *n, *kappa;
  double *nq;
  SHELL *sa, *sb;
  
  if (sfac_file) {
    SFACStatement("SlaterCoeff", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!(PyArg_ParseTuple(args, "sOss", &fn, &p, &q1, &q2))) return NULL;

  nlev = SelectLevels(p, &ilev);
  na = GetAverageConfigFromString(&n, &kappa, &nq, q1);
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
  nb = GetAverageConfigFromString(&n, &kappa, &nq, q2);
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
    SlaterCoeff(fn, nlev, ilev, na, sa, nb, sb);
    free(ilev);
    free(sa);
    free(sb);
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PGeneralizedMoment(PyObject *self, PyObject *args) {
  char *fn;
  int n0, k0, n1, k1, m;
  double e1;  

  if (sfac_file) {
    SFACStatement("GeneralizedMoment", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  e1 = 0;
  if (!(PyArg_ParseTuple(args, "siiiii|d", &fn, &m, &n0, &k0, &n1, &k1, &e1))) 
    return NULL;
  
  PrintGeneralizedMoments(fn, m, n0, k0, n1, k1, e1);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PPrintQED(PyObject *self, PyObject *args) {
  if (sfac_file) {
    SFACStatement("PrintQED", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  PrintQED();
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PPrintNucleus(PyObject *self, PyObject *args) {
  if (sfac_file) {
    SFACStatement("PrintNucleus", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  int m;
  char *fn;

  m = 0;
  fn = NULL;
  if (!(PyArg_ParseTuple(args, "|is", &m, &fn))) return NULL;
  PrintNucleus(m, fn);
  
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PBreitX(PyObject *self, PyObject *args) {
  double e;
  int k0, k1, k, m;
  
  if (sfac_file) {
    SFACStatement("BreitX", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  e = -1.0;
  if (!(PyArg_ParseTuple(args, "iiii|d", &k0, &k1, &k, &m, &e)))
    return NULL;

  ORBITAL *orb0 = GetOrbital(k0);
  ORBITAL *orb1 = GetOrbital(k1);
  if (orb0 != NULL && orb1 != NULL) {
    BreitX(orb0, orb1, k, m, 0, 1, e, NULL);
  }
  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PSavePotential(PyObject *self, PyObject *args) {
  char *fn;
  POTENTIAL *p;
   
  if (sfac_file) {
    SFACStatement("SavePotential", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!(PyArg_ParseTuple(args, "s", &fn))) {
    return NULL;
  }

  p = RadialPotential();
  SavePotential(fn, p);
  
  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PRestorePotential(PyObject *self, PyObject *args) {
  char *fn;
  POTENTIAL *p;
   
  if (sfac_file) {
    SFACStatement("RestorePotential", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!(PyArg_ParseTuple(args, "s", &fn))) {
    return NULL;
  }

  p = RadialPotential();
  RestorePotential(fn, p);
  
  Py_INCREF(Py_None);
  return Py_None;
}
 
 
static PyObject *PModifyPotential(PyObject *self, PyObject *args) {
  char *fn;
  POTENTIAL *p;
   
  if (sfac_file) {
    SFACStatement("ModifyPotential", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!(PyArg_ParseTuple(args, "s", &fn))) {
    return NULL;
  }

  p = RadialPotential();
  ModifyPotential(fn, p);
  
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PWallTime(PyObject *self, PyObject *args) {
  if (sfac_file) {
    SFACStatement("WallTime", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  int m = 0;
  char *s;
  if (!(PyArg_ParseTuple(args, "s|i", &s, &m))) return NULL;
  PrintWallTime(s, m);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PInitializeMPI(PyObject *self, PyObject *args) {
  if (sfac_file) {
    SFACStatement("InitializeMPI", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

#ifdef USE_MPI
  int n = -1;
  if (!(PyArg_ParseTuple(args, "|i", &n))) {
    return NULL;
  }
  InitializeMPI(n, 0);
#endif
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMPIRank(PyObject *self, PyObject *args) {
  if (sfac_file) {
    SFACStatement("MPIRank", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  int n, k;
  k = MPIRank(&n);
  return Py_BuildValue("[ii]", k, n);
}

static PyObject *PMemUsed(PyObject *self, PyObject *args) {
  if (sfac_file) {
    SFACStatement("MemUsed", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  double m = msize();
  return Py_BuildValue("d", m);
}

static PyObject *PSetOrbMap(PyObject *self, PyObject *args) {
  if (sfac_file) {
    SFACStatement("SetOrbMap", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  int k=0, n0=0, n1=0, n2=0;
  
  if (!(PyArg_ParseTuple(args, "|iiii", &k, &n0, &n1, &n2))) return NULL;
  SetOrbMap(k, n0, n1, n2);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PFinalizeMPI(PyObject *self, PyObject *args) {
  if (sfac_file) {
    SFACStatement("FinalizeMPI", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

#if USE_MPI == 1
  FinalizeMPI();
#endif
  
  Py_INCREF(Py_None);
  return Py_None;
}
  
static struct PyMethodDef fac_methods[] = {
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
  {"RMatrixRefine", PRMatrixRefine, METH_VARARGS}, 
  {"RMatrixCE", PRMatrixCE, METH_VARARGS}, 
  {"TestRMatrix", PTestRMatrix, METH_VARARGS}, 
  {"SetSlaterCut", PSetSlaterCut, METH_VARARGS}, 
  {"RMatrixBoundary", PRMatrixBoundary, METH_VARARGS}, 
  {"RMatrixTargets", PRMatrixTargets, METH_VARARGS}, 
  {"RMatrixBasis", PRMatrixBasis, METH_VARARGS}, 
  {"RMatrixSurface", PRMatrixSurface, METH_VARARGS}, 
  {"Print", PPrint, METH_VARARGS},
  {"Asymmetry", PAsymmetry, METH_VARARGS},
  {"ReadConfig", PReadConfig, METH_VARARGS},
  {"Config", (PyCFunction) PConfig, METH_VARARGS|METH_KEYWORDS},
  {"RemoveConfig", PRemoveConfig, METH_VARARGS},
  {"ListConfig", PListConfig, METH_VARARGS},
  {"GetConfigNR", PGetConfigNR, METH_VARARGS},
  {"Closed", PClosed, METH_VARARGS},
  {"CutMixing", PCutMixing, METH_VARARGS},
  {"AvgConfig", PAvgConfig, METH_VARARGS},
  {"AddConfig", PAddConfig, METH_VARARGS},
  {"AIBranch", PAIBranch, METH_VARARGS},
  {"AITable", PAITable, METH_VARARGS},
  {"AITableMSub", PAITableMSub, METH_VARARGS},
  {"BasisTable", PBasisTable, METH_VARARGS},
  {"CECross", PInterpCross, METH_VARARGS},
  {"CERate", PMaxwellRate, METH_VARARGS},
  {"InterpCross", PInterpCross, METH_VARARGS},
  {"MaxwellRate", PMaxwellRate, METH_VARARGS},
  {"CETable", PCETable, METH_VARARGS},
  {"CETableEB", PCETableEB, METH_VARARGS},
  {"CETableMSub", PCETableMSub, METH_VARARGS},
  {"CheckEndian", PCheckEndian, METH_VARARGS},
  {"CITable", PCITable, METH_VARARGS},
  {"CITableMSub", PCITableMSub, METH_VARARGS},
  {"ClearLevelTable", PClearLevelTable, METH_VARARGS},
  {"ClearOrbitalTable", PClearOrbitalTable, METH_VARARGS},
  {"ConfigEnergy", PConfigEnergy, METH_VARARGS},
  {"CloseSFAC", PCloseSFAC, METH_VARARGS},
  {"ConvertToSFAC", PConvertToSFAC, METH_VARARGS},
  {"AdjustEnergy", PAdjustEnergy, METH_VARARGS},
  {"CorrectEnergy", PCorrectEnergy, METH_VARARGS},
  {"DROpen", PDROpen, METH_VARARGS},
  {"FreeExcitationQk", PFreeExcitationQk, METH_VARARGS},
  {"FreeIonizationQk", PFreeIonizationQk, METH_VARARGS},
  {"FreeMemENTable", PFreeMemENTable, METH_VARARGS},
  {"FreeMultipole", PFreeMultipole, METH_VARARGS},
  {"FreeSlater", PFreeSlater, METH_VARARGS},
  {"FreeResidual", PFreeResidual, METH_VARARGS},
  {"FreeRecPk", PFreeRecPk, METH_VARARGS},
  {"FreeRecQk", PFreeRecQk, METH_VARARGS},
  {"GetCFPOld", PGetCFPOld, METH_VARARGS},
  {"GetW3j", PGetW3j, METH_VARARGS},
  {"GetW6j", PGetW6j, METH_VARARGS},
  {"GetW9j", PGetW9j, METH_VARARGS},
  {"GetCG", PGetCG, METH_VARARGS},
  {"WignerDMatrix", PWignerDMatrix, METH_VARARGS},
  {"GetPotential", PGetPotential, METH_VARARGS},
  {"Info", PInfo, METH_VARARGS},
  {"StructureMBPT", PStructureMBPT, METH_VARARGS},
  {"TransitionMBPT", PTransitionMBPT, METH_VARARGS},
  {"MemENTable", PMemENTable, METH_VARARGS},
  {"LevelInfor", PLevelInfor, METH_VARARGS},
  {"LevelInfo", PLevelInfor, METH_VARARGS},
  {"OptimizeRadial", POptimizeRadial, METH_VARARGS},
  {"PrepAngular", PPrepAngular, METH_VARARGS},
  {"RadialOverlaps", PRadialOverlaps, METH_VARARGS},
  {"FreezeOrbital", PFreezeOrbital, METH_VARARGS},
  {"RefineRadial", PRefineRadial, METH_VARARGS},
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
  {"Reinit", (PyCFunction) PReinit, METH_VARARGS|METH_KEYWORDS},
  {"RRTable", PRRTable, METH_VARARGS},
  {"RRMultipole", PRRMultipole, METH_VARARGS},
  {"SetAICut", PSetAICut, METH_VARARGS},
  {"SetAngZOptions", PSetAngZOptions, METH_VARARGS},
  {"SetAngZCut", PSetAngZCut, METH_VARARGS},
  {"SetCILevel", PSetCILevel, METH_VARARGS},
  {"SetMixCut", PSetMixCut, METH_VARARGS},
  {"SetAtom", PSetAtom, METH_VARARGS},
  {"SetExtraPotential", PSetExtraPotential, METH_VARARGS},
  {"SetAvgConfig", PSetAvgConfig, METH_VARARGS},
  {"SetBoundary", PSetBoundary, METH_VARARGS},
  {"GetBoundary", PGetBoundary, METH_VARARGS},
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
  {"SetCIPWOptions", PSetCIPWOptions, METH_VARARGS},
  {"SetCIPWGrid", PSetCIPWGrid, METH_VARARGS},
  {"SetCIQkMode", PSetCIQkMode, METH_VARARGS},
  {"SetHydrogenicNL", PSetHydrogenicNL, METH_VARARGS},
  {"SetMaxRank", PSetMaxRank, METH_VARARGS},
  {"SetOptimizeControl", PSetOptimizeControl, METH_VARARGS},
  {"SetPEGrid", PSetPEGrid, METH_VARARGS},
  {"SetPEGridLimits", PSetPEGridLimits, METH_VARARGS},  
  {"SetOrbMap", PSetOrbMap, METH_VARARGS},
  {"SetRadialGrid", PSetRadialGrid, METH_VARARGS},
  {"SetPotentialMode", PSetPotentialMode, METH_VARARGS},
  {"SolvePseudo", PSolvePseudo, METH_VARARGS},
  {"SetRecPWLimits", PSetRecPWLimits, METH_VARARGS},
  {"SetRecPWOptions", PSetRecPWOptions, METH_VARARGS},
  {"SetRecQkMode", PSetRecQkMode, METH_VARARGS},
  {"SetRecSpectator", PSetRecSpectator, METH_VARARGS},
  {"SetRRTEGrid", PSetRRTEGrid, METH_VARARGS},
  {"SetScreening", PSetScreening, METH_VARARGS},
  {"SetSE", PSetSE, METH_VARARGS},
  {"SetModSE", PSetModSE, METH_VARARGS},
  {"OptimizeModSE", POptimizeModSE, METH_VARARGS},
  {"SetVP", PSetVP, METH_VARARGS},
  {"SetBreit", PSetBreit, METH_VARARGS},
  {"SetMS", PSetMS, METH_VARARGS},
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
  {"TestAngular", PTestAngular, METH_VARARGS},
  {"CoulombBethe", PCoulombBethe, METH_VARARGS}, 
  {"TestHamilton", PTestHamilton, METH_VARARGS}, 
  {"TestIntegrate", PTestIntegrate, METH_VARARGS}, 
  {"TestMyArray", PTestMyArray, METH_VARARGS},        
  {"ReportMultiStats", PReportMultiStats, METH_VARARGS},     
  {"ElectronDensity", PElectronDensity, METH_VARARGS},  
  {"ExpectationValue", PExpectationValue, METH_VARARGS},   
  {"TRTable", PTransitionTable, METH_VARARGS}, 
  {"TransitionTable", PTransitionTable, METH_VARARGS},     
  {"TRTableEB", PTransitionTableEB, METH_VARARGS},    
  {"PolarizeCoeff", PPolarizeCoeff, METH_VARARGS}, 
  {"TRBranch", PTRBranch, METH_VARARGS}, 
  {"TRRateH", PTRRateH, METH_VARARGS},  
  {"PICrossH", PPICrossH, METH_VARARGS},  
  {"RRCrossH", PRRCrossH, METH_VARARGS},  
  {"RRCrossHn", PRRCrossHn, METH_VARARGS},  
  {"TotalRRCross", PTotalRRCross, METH_VARARGS}, 
  {"TotalPICross", PTotalPICross, METH_VARARGS}, 
  {"TotalCICross", PTotalCICross, METH_VARARGS}, 
  {"WaveFuncTable", PWaveFuncTable, METH_VARARGS}, 
  {"Y5N", PY5N, METH_VARARGS},
  {"DiracCoulomb", PDiracCoulomb, METH_VARARGS},
  {"CoulombPhase", PCoulombPhase, METH_VARARGS},
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
  {"SetFields", PSetFields, METH_VARARGS},
  {"StructureEB", PStructureEB, METH_VARARGS},
  {"SetTransitionMode", PSetTransitionMode, METH_VARARGS},
  {"SetTransitionGauge", PSetTransitionGauge, METH_VARARGS},
  {"SetTransitionMaxE", PSetTransitionMaxE, METH_VARARGS},
  {"SetTransitionMaxM", PSetTransitionMaxM, METH_VARARGS},
  {"CoulMultipole", PCoulMultip, METH_VARARGS},
  {"BreitX", PBreitX, METH_VARARGS},
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
  {NULL, NULL}
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "fac",
  NULL,
  -1,
  fac_methods,
};
#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_fac(void){

#else
#define INITERROR return

void
initfac(void){
#endif
  PyObject *m, *d;
  char v[10];
  char sp[2];
  char *ename;
  double *emass;
  int i;

  #if PY_MAJOR_VERSION >= 3
    m = PyModule_Create(&moduledef);
  #else
    m = Py_InitModule("fac", fac_methods);
  #endif

  d = PyModule_GetDict(m);
  ErrorObject = Py_BuildValue("s", "fac.error");
  PyDict_SetItemString(d, "error", ErrorObject);

  if(m == NULL) INITERROR;

  if (InitFac() < 0) {
    onError("initilization failed\n");

#if PY_MAJOR_VERSION >= 3
    return m;
#else
    return;
#endif
  }

  SPECSYMBOL = PyList_New(MAX_SPEC_SYMBOLS);
  sp[1] = '\0';
  for (i = 0; i < MAX_SPEC_SYMBOLS; i++) {
    SpecSymbol(sp, i);
    PyList_SetItem(SPECSYMBOL, i, Py_BuildValue("s", sp));
  }

  sprintf(v, "%d.%d.%d", VERSION, SUBVERSION, SUBSUBVERSION);
  PFACVERSION = PyUnicode_FromString(v);

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

#if PY_MAJOR_VERSION >= 3
  return m;
#endif
}

