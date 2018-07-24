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

#include "Python.h"
#include <stdio.h>
#include <string.h>

#include "crm.h"

#if PY_MAJOR_VERSION >= 3
  #define PyUnicode_AsString(x) PyBytes_AsString(PyUnicode_AsEncodedString((x), "utf-8", "strict"))
#else
  #define PyUnicode_AsString PyString_AsString
  #define PyLong_AsLong PyInt_AsLong
  #define PyLong_Check PyInt_Check
#endif

static PyObject *ErrorObject;
#define onError(message) {PyErr_SetString(ErrorObject, message);}

static FILE *scrm_file = NULL;

static void SCRMStatement(char *func, PyObject *args, PyObject *kargs) {
  int i, n, nargs;
  PyObject *sargs;
  PyObject *klist;
  PyObject *kvar;
  PyObject *p, *q;
  char *s1, *s2;
  
  fprintf(scrm_file, "%s", func);
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
    fprintf(scrm_file, "%c", s1[i]);
  }
  if (kargs) {
    klist = PyDict_Items(kargs);
    n = PyList_Size(klist);
    for (i = 0; i < n; i++) {
      if (nargs > 0 || i > 0) fprintf(scrm_file, ", ");
      p = PyList_GetItem(klist, i);
      q = PyTuple_GetItem(p, 0);
      s2 = PyUnicode_AsString(q);
      fprintf(scrm_file, "%s=", s2);
      q = PyTuple_GetItem(p, 1);
      kvar = PyObject_Str(q);
      s2 = PyUnicode_AsString(kvar);
      if (PyUnicode_Check(q)) {
	fprintf(scrm_file, "'%s'", s2);
      } else {
	fprintf(scrm_file, "%s", s2);
      }
      Py_XDECREF(kvar);
    }
    Py_XDECREF(klist);
  }

  fprintf(scrm_file, ")\n");
    
  Py_XDECREF(sargs);

  return;
}

static PyObject *PConvertToSCRM(PyObject *self, PyObject *args) {
  char *fn;

  scrm_file = NULL;
  if (!PyArg_ParseTuple(args, "|s", &fn)) return NULL;

  if (fn) {
    scrm_file = fopen(fn, "w");
    if (scrm_file == NULL) return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}    
  
static PyObject *PCloseSCRM(PyObject *self, PyObject *args) {

  fclose(scrm_file);
  scrm_file = NULL;

  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PCheckEndian(PyObject *self, PyObject *args) {
  char *fn;
  TFILE *f;
  F_HEADER fh;
  int i, swp;

  if (scrm_file) {
    SCRMStatement("CheckEndian", args, NULL);
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

static PyObject *PPrint(PyObject *self, PyObject *args) {
  PyObject *p, *q;
  char *s;
  int i, n;

  if (scrm_file) {
    SCRMStatement("Print", args, NULL);
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

static PyObject *PEleDist(PyObject *self, PyObject *args) {
  int n;
  char *fn;
  
  if (scrm_file) {
    SCRMStatement("EleDist", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "si", &fn, &n)) return NULL;
  EleDist(fn, n);

  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetEleDist(PyObject *self, PyObject *args) {
  PyObject *p;
  int i, n, k, np;
  double *par = NULL;
  char *fn;

  if (scrm_file) {
    SCRMStatement("SetEleDist", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n < 1) return NULL;
  p = PyTuple_GetItem(args, 0);
  i = PyLong_AsLong(p);
  if (i == -1) {
    p = PyTuple_GetItem(args, 1);
    fn = PyUnicode_AsString(p);
    np = DistFromFile(fn, &par);
    if (np <= 0) return NULL;
  } else {
    np = n-1;
    if (np > 0) {
      par = (double *) malloc(sizeof(double)*np);
    }  
    for (k = 1; k < n; k++) {
      p = PyTuple_GetItem(args, k);
      par[k-1] = PyFloat_AsDouble(p);
    }
  }
  if (SetEleDist(i, np, par) < 0) return NULL;
  if (np > 0) free(par);
 
  Py_INCREF(Py_None);
  return Py_None;
} 
 
static PyObject *PPhoDist(PyObject *self, PyObject *args) {
  int n;
  char *fn;
  
  if (scrm_file) {
    SCRMStatement("PhoDist", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "si", &fn, &n)) return NULL;
  PhoDist(fn, n);

  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetPhoDist(PyObject *self, PyObject *args) {
  PyObject *p;
  int i, n, k, np;
  double *par = NULL;
  char *fn;

  if (scrm_file) {
    SCRMStatement("SetPhoDist", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n < 1) return NULL;
  p = PyTuple_GetItem(args, 0);
  i = PyLong_AsLong(p);
  if (i == -1) {
    p = PyTuple_GetItem(args, 1);
    fn = PyUnicode_AsString(p);
    np = DistFromFile(fn, &par);
    if (np <= 0) return NULL;
  } else {
    np = n-1;
    if (np > 0) {
      par = (double *) malloc(sizeof(double)*np);
    }
    for (k = 1; k < n; k++) {
      p = PyTuple_GetItem(args, k);
      par[k-1] = PyFloat_AsDouble(p);
    }
  }
  if (SetPhoDist(i, np, par) < 0) return NULL;
  if (np > 0) free(par);
 
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetNumSingleBlocks(PyObject *self, PyObject *args) {
  int n;
  
  if (scrm_file) {
    SCRMStatement("SetNumSingleBlocks", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &n)) return NULL;
  SetNumSingleBlocks(n);
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetExtrapolate(PyObject *self, PyObject *args) {
  int n;
  
  if (scrm_file) {
    SCRMStatement("SetExtrapolate", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &n)) return NULL;
  SetExtrapolate(n);

  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetInnerAuger(PyObject *self, PyObject *args) {
  int n;
  
  if (scrm_file) {
    SCRMStatement("SetInnerAuger", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &n)) return NULL;
  SetInnerAuger(n);

  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetEMinAI(PyObject *self, PyObject *args) {
  double e;
  
  if (scrm_file) {
    SCRMStatement("SetEMinAI", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &e)) return NULL;
  SetEMinAI(e);

  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetEleDensity(PyObject *self, PyObject *args) {
  double den;
  
  if (scrm_file) {
    SCRMStatement("SetEleDensity", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &den)) return NULL;
  SetEleDensity(den);
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetPhoDensity(PyObject *self, PyObject *args) {
  double den;
  
  if (scrm_file) {
    SCRMStatement("SetPhoDensity", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &den)) return NULL;
  SetPhoDensity(den);
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetRateAccuracy(PyObject *self, PyObject *args) {
  double epsrel, epsabs;

  if (scrm_file) {
    SCRMStatement("SetRateAccuracy", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  epsabs = -1.0;
  if (!PyArg_ParseTuple(args, "d|d", &epsrel, &epsabs)) return NULL;
  SetRateAccuracy(epsrel, epsabs);
  
  Py_INCREF(Py_None);
  return Py_None;
}  
    
static PyObject *PSetCascade(PyObject *self, PyObject *args) {
  int c;
  double a;
  
  if (scrm_file) {
    SCRMStatement("SetCascade", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  a = 0.0;
  if (!PyArg_ParseTuple(args, "i|d", &c, &a)) return NULL;
  SetCascade(c, a);

  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetIteration(PyObject *self, PyObject *args) {
  int maxiter;
  double a, s;
  
  if (scrm_file) {
    SCRMStatement("SetIteration", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  s = -1.0;
  maxiter = -1;
  if (!PyArg_ParseTuple(args, "d|di", &a, &s, &maxiter)) return NULL;
  SetIteration(a, s, maxiter);
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetBlocks(PyObject *self, PyObject *args) {
  char *ifn;
  double n;

  if (scrm_file) {
    SCRMStatement("SetBlocks", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = 0.0;
  ifn = NULL;
  if (!PyArg_ParseTuple(args, "|ds", &n, &ifn)) return NULL;
  SetBlocks(n, ifn);

  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PInitBlocks(PyObject *self, PyObject *args) {
  
  if (scrm_file) {
    SCRMStatement("InitBlocks", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  InitBlocks();
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PLevelPopulation(PyObject *self, PyObject *args) {
  
  if (scrm_file) {
    SCRMStatement("LevelPopulation", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  LevelPopulation();
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PCascade(PyObject *self, PyObject *args) {
  
  if (scrm_file) {
    SCRMStatement("Cascade", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  Cascade();
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSpecTable(PyObject *self, PyObject *args) {
  char *fn;
  double smin;
  int rrc;

  if (scrm_file) {
    SCRMStatement("SpecTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  smin = EPS10;
  rrc = 0;
  if (!PyArg_ParseTuple(args, "s|id", &fn, &rrc, &smin)) return NULL;
  
  SpecTable(fn, rrc, smin);
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PPlotSpec(PyObject *self, PyObject *args) {
  char *fn1, *fn2;
  double emin, emax, de, smin;
  int nele, type;

  if (scrm_file) {
    SCRMStatement("PlotSpec", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  smin = EPS6;
  if (!PyArg_ParseTuple(args, "ssiiddd|d", 
			&fn1, &fn2, &nele, &type, &emin, &emax, &de, &smin))
    return NULL;
      
  PlotSpec(fn1, fn2, nele, type, emin, emax, de, smin);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PTabNLTE(PyObject *self, PyObject *args) {
  char *fn1, *fn2, *fn3, *fn;
  double xmin, xmax, dx;

  if (scrm_file) {
    SCRMStatement("TabNLTE", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  fn3 = NULL;
  if (!PyArg_ParseTuple(args, "sssddd|s", &fn1, &fn2, &fn, 
			&xmin, &xmax, &dx, &fn3)) 
    return NULL;

  TabNLTE(fn1, fn2, fn3, fn, xmin, xmax, dx);

  Py_INCREF(Py_None);
  return Py_None;
}    
  
static PyObject *PAddIon(PyObject *self, PyObject *args) {
  char *fn;
  int i;
  double n;

  if (scrm_file) {
    SCRMStatement("AddIon", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "ids", &i, &n, &fn))
    return NULL;
      
  AddIon(i, n, fn);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCERates(PyObject *self, PyObject *args) {
  int inv;

  if (scrm_file) {
    SCRMStatement("SetCERates", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &inv)) return NULL;
  SetCERates(inv);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetTRRates(PyObject *self, PyObject *args) {
  int inv;

  if (scrm_file) {
    SCRMStatement("SetTRRates", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &inv)) return NULL;
  SetTRRates(inv);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetCIRates(PyObject *self, PyObject *args) {
  int inv;

  if (scrm_file) {
    SCRMStatement("SetCIRates", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &inv)) return NULL;
  SetCIRates(inv);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetRRRates(PyObject *self, PyObject *args) {
  int inv;
  
  if (scrm_file) {
    SCRMStatement("SetRRRates", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &inv)) return NULL;
  SetRRRates(inv);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetAIRates(PyObject *self, PyObject *args) {
  int inv;
  
  if (scrm_file) {
    SCRMStatement("SetAIRates", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "i", &inv)) return NULL;
  SetAIRates(inv);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetAIRatesInner(PyObject *self, PyObject *args) {
  char *fn;

  if (scrm_file) {
    SCRMStatement("SetAIRatesInner", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  SetAIRatesInner(fn);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PPrintTable(PyObject *self, PyObject *args) { 
  char *fn1, *fn2;
  int v;
  
  if (scrm_file) {
    SCRMStatement("PrintTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  v = 1;
  if (!PyArg_ParseTuple(args, "ss|i", &fn1, &fn2, &v)) return NULL;
  PrintTable(fn1, fn2, v);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PReinitCRM(PyObject *self, PyObject *args) { 
  int m;

  if (scrm_file) {
    SCRMStatement("ReinitCRM", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  m = 0;
  if (!PyArg_ParseTuple(args, "|i", &m)) return NULL;
  ReinitCRM(m);
 
  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PRateTable(PyObject *self, PyObject *args) { 
  PyObject *p;
  int i, nc, md;
  char **sc;
  char *fn;
    
  if (scrm_file) {
    SCRMStatement("RateTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  p = NULL;
  md = 0;
  if (!PyArg_ParseTuple(args, "s|Oi", &fn, &p, &md)) return NULL;
  
  nc = 0;
  sc = NULL;
  if (p) {
    if (!PyList_Check(p)) return NULL;
    nc = PyList_Size(p);
    if (nc > 0) {
      sc = (char **) malloc(sizeof(char *)*nc);
      for (i = 0; i < nc; i++) {
	sc[i] = PyUnicode_AsString(PyList_GetItem(p, i));
      }
    }
  }
  
  RateTable(fn, nc, sc, md);

  if (nc > 0) {
    free(sc);
  }

  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PSetAbund(PyObject *self, PyObject *args) { 
  int nele;
  double a;
    
  if (scrm_file) {
    SCRMStatement("SetAbund", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "id", &nele, &a)) return NULL;
  SetAbund(nele, a);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSetRateMultiplier(PyObject *self, PyObject *args) { 
  int nele;
  int t;
  double a;
    
  if (scrm_file) {
    SCRMStatement("SetRateMultiplier", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "iid", &nele, &t, &a)) return NULL;
  SetRateMultiplier(nele, t, a);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PSelectLines(PyObject *self, PyObject *args) {
  char *ifn, *ofn;
  double emin, emax, fmin;
  int nele, type;

  if (scrm_file) {
    SCRMStatement("SelectLines", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  fmin = 0.0;
  if (!PyArg_ParseTuple(args, "ssiidd|d", 
			&ifn, &ofn, &nele, &type, &emin, &emax, &fmin)) 
    return NULL;
  
  SelectLines(ifn, ofn, nele, type, emin, emax, fmin);
  
  Py_INCREF(Py_None);
  return Py_None;
}  

static int ShellIndexFromString(PyObject *s) {
  int ks;
  char *c;

  c = PyUnicode_AsString(s);
  if (strcmp(c, "1s") == 0) ks = 1;
  else if (strcmp(c, "2s") == 0) ks = 2;
  else if (strcmp(c, "2p") == 0) ks = 3;
  else if (strcmp(c, "3s") == 0) ks = 4;
  else if (strcmp(c, "3p") == 0) ks = 5;
  else if (strcmp(c, "3d") == 0) ks = 6;
  else if (strcmp(c, "4s") == 0) ks = 7;
  else ks = 0;
  return ks;
}

static PyObject *PBremss(PyObject *self, PyObject *args) {
  int z;
  double te, e, a;
  
  e = -1.0;
  if (!PyArg_ParseTuple(args, "id|d", &z, &te, &e)) return NULL;
  
  a = BremssNR(z, te, e);
  
  return Py_BuildValue("d", a);
}

static PyObject *PEPhFit(PyObject *self, PyObject *args) {
  PyObject *s;
  int z, nele, ks;
  double e;
  
  if (!PyArg_ParseTuple(args, "iiO", &z, &nele, &s)) return NULL;
  ks = -1;
  if (PyUnicode_Check(s)) {
    ks = ShellIndexFromString(s);
  } else if (PyLong_Check(s)) {
    ks = PyLong_AsLong(s);
  }
  if (ks < 1) return NULL;

  e = EPhFit2(z, nele, ks);
  return Py_BuildValue("d", e);
}

static PyObject *PPhFit(PyObject *self, PyObject *args) {
  PyObject *s;
  int z, nele, ks;
  double e, r;
  
  if (!PyArg_ParseTuple(args, "iidO", &z, &nele, &e, &s)) return NULL;
  ks = -1;
  if (PyUnicode_Check(s)) {
    ks = ShellIndexFromString(s);
  } else if (PyLong_Check(s)) {
    ks = PyLong_AsLong(s);
  }
  if (ks < 1) return NULL;

  r = PhFit2(z, nele, ks, e);
  return Py_BuildValue("d", r);
}

static PyObject *PNRRFit(PyObject *self, PyObject *args) {
  int z, nele;
  double t, r;
  
  if (!PyArg_ParseTuple(args, "iid", &z, &nele, &t)) return NULL;
  r = NRRFit(z, nele, t);
  return Py_BuildValue("d", r);
}

static PyObject *PNDRFit(PyObject *self, PyObject *args) {
  int z, nele;
  double t, r;
  
  if (!PyArg_ParseTuple(args, "iid", &z, &nele, &t)) return NULL;
  r = NDRFit(z, nele, t);
  return Py_BuildValue("d", r);
}

static PyObject *PRRFit(PyObject *self, PyObject *args) {
  int z, nele;
  double t, r;
  
  if (!PyArg_ParseTuple(args, "iid", &z, &nele, &t)) return NULL;
  r = RRFit(z, nele, t);
  return Py_BuildValue("d", r);
}

static PyObject *PDRFit(PyObject *self, PyObject *args) {
  int z, nele;
  double t, r;
  
  if (!PyArg_ParseTuple(args, "iid", &z, &nele, &t)) return NULL;
  r = DRFit(z, nele, t);
  return Py_BuildValue("d", r);
}

static PyObject *PRBeli(PyObject *self, PyObject *args) {
  int z, nele;
  double t, d, a, r;

  if (!PyArg_ParseTuple(args, "iid", &z, &nele, &t)) return NULL;
  r = RBeli(z, nele, t, &a, &d);
  return Py_BuildValue("(ddd)", r, a, d);
}

static PyObject *PCBeli(PyObject *self, PyObject *args) {
  int z, nele;
  double ene, d, a, err, r;

  if (!PyArg_ParseTuple(args, "iid", &z, &nele, &ene)) return NULL;
  r = CBeli(z, nele, ene, &a, &d, &err);
  return Py_BuildValue("(dddd)", r, a, d, err);
}

static PyObject *PCFit(PyObject *self, PyObject *args) {
  int z, nele;
  double t, r, a, d;
  
  if (!PyArg_ParseTuple(args, "iid", &z, &nele, &t)) return NULL;
  r = CFit(z, nele, t, &a, &d);
  return Py_BuildValue("(ddd)", r, a, d);
}

static PyObject *PColFit(PyObject *self, PyObject *args) {
  PyObject *s;
  int z, nele, ks;
  double t, r, a, d;
  
  s = NULL;
  if (!PyArg_ParseTuple(args, "iid|O", &z, &nele, &t, &s)) return NULL;
  if (s == NULL) {
    ks = 0;
  } else if (PyUnicode_Check(s)) {
    ks = ShellIndexFromString(s);
  } else if (PyLong_Check(s)) {
    ks = PyLong_AsLong(s);
  } else {
    return NULL;
  }
  if (ks < 0) return NULL;
  r = ColFit(z, nele, ks, t, &a, &d);
  return Py_BuildValue("(ddd)", r, a, d);
}

static PyObject *PCColFit(PyObject *self, PyObject *args) {
  PyObject *s;
  int z, nele, ks;
  double t, r, a, d;
  
  s = NULL;
  if (!PyArg_ParseTuple(args, "iid|O", &z, &nele, &t, &s)) return NULL;
  if (s == NULL) {
    ks = 0;
  } else if (PyUnicode_Check(s)) {
    ks = ShellIndexFromString(s);
  } else if (PyLong_Check(s)) {
    ks = PyLong_AsLong(s);
  } else {
    return NULL;
  }
  if (ks < 0) return NULL;
  r = CColFit(z, nele, ks, t, &a, &d);
  return Py_BuildValue("(ddd)", r, a, d);
}

static PyObject *PEColFit(PyObject *self, PyObject *args) {
  PyObject *s;
  int z, nele, ks;
  double e;
  
  s = NULL;
  if (!PyArg_ParseTuple(args, "ii|O", &z, &nele, &s)) return NULL;
  if (s == NULL) {
    ks = 0;
  } else if (PyUnicode_Check(s)) {
    ks = ShellIndexFromString(s);
  } else if (PyLong_Check(s)) {
    ks = PyLong_AsLong(s);
  } else {
    return NULL;
  }
  if (ks < 0) return NULL;
  e = EColFit(z, nele, ks);
  return Py_BuildValue("d", e);
}

static PyObject *PEBeli(PyObject *self, PyObject *args) {
  int z, nele;
  double e;
  
  if (!PyArg_ParseTuple(args, "ii", &z, &nele)) return NULL;
  e = EBeli(z, nele);
  return Py_BuildValue("d", e);
}

static PyObject *PIonis(PyObject *self, PyObject *args) {
  int z, nele, m;
  double t, total, a, d;

  m = 1;
  if (!PyArg_ParseTuple(args, "iid|i", &z, &nele, &t, &m)) return NULL;
  
  total = Ionis(z, nele, t, &a, &d, m);
  
  return Py_BuildValue("(ddd)", total, a, d);
}

static PyObject *PRRRateH(PyObject *self, PyObject *args) {
  int n;
  double z, t, r, top;

  if (!PyArg_ParseTuple(args, "did", &z, &n, &t)) return NULL;
  
  r = RRRateHydrogenic(t, z, n, &top);
  
  return Py_BuildValue("(dd)", r, top);
}

static PyObject *PRecomb(PyObject *self, PyObject *args) {
  int z, nele, m;
  double t, total, r, d;

  m = 1;
  if (!PyArg_ParseTuple(args, "iid|i", &z, &nele, &t, &m)) return NULL;
  
  total = Recomb(z, nele, t, &r, &d, m);
  
  return Py_BuildValue("(ddd)", total, r, d);
}

static PyObject *PFracAbund(PyObject *self, PyObject *args) {
  PyObject *pa;
  int z, i, im, rm;
  double t, *a;

  im = 1;
  rm = 1;
  if (!PyArg_ParseTuple(args, "id|ii", &z, &t, &im, &rm)) return NULL;
  
  a = (double *) malloc(sizeof(double)*(z+1));
  FracAbund(z, t, a, im, rm);
  pa = Py_BuildValue("[]");
  for (i = 0; i <= z; i++) {
    PyList_Append(pa, Py_BuildValue("d", a[i]));
  }
  free(a);

  return pa;
}

static PyObject *PMaxAbund(PyObject *self, PyObject *args) {
  PyObject *pa;
  int i, z, nele, im, rm;
  double eps, *a, tmax;

  im = 1;
  rm = 1;
  eps = 1E-4;
  if (!PyArg_ParseTuple(args, "ii|iid", &z, &nele, &im, &rm, &eps)) 
    return NULL;
  
  a = (double *) malloc(sizeof(double)*(z+1));
  tmax = MaxAbund(z, nele, a, eps, im, rm);
  pa = Py_BuildValue("[]");
  for (i = 0; i <= z; i++) {
    PyList_Append(pa, Py_BuildValue("d", a[i]));
  }
  free(a);
  return Py_BuildValue("(dO)", tmax, pa);
}  

static PyObject *PDRBranch(PyObject *self, PyObject *args) {
  
  if (scrm_file) {
    SCRMStatement("DRBranch", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  DRBranch();

  Py_INCREF(Py_None);  
  return Py_None;
}

static PyObject *PDRStrength(PyObject *self, PyObject *args) {
  int n, i, m;
  char *s;
  
  if (scrm_file) {
    SCRMStatement("DRStrength", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  i = 0;
  m = 0;
  if (!PyArg_ParseTuple(args, "si|ii", &s, &n, &m, &i)) 
    return NULL;

  DRStrength(s, n, m, i);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PRydBranch(PyObject *self, PyObject *args) {
  char *fn, *ofn;
  int n0, n1;

  if (scrm_file) {
    SCRMStatement("RydBranch", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  n1 = -1;
  if (!PyArg_ParseTuple(args, "ssi|i", &fn, &ofn, &n0, &n1)) 
    return NULL;
  
  RydBranch(fn, ofn, n0, n1);

  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PTwoPhoton(PyObject *self, PyObject *args) {
  int t;
  double z;
  
  if (!PyArg_ParseTuple(args, "di", &z, &t))
    return NULL;

  z = TwoPhotonRate(z, t);
  
  return Py_BuildValue("d", z);
}

static PyObject *PDumpRates(PyObject *self, PyObject *args) {
  int k, m, imax, a;
  char *fn;

  if (scrm_file) {
    SCRMStatement("DumpRates", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  imax = -1;
  a = 0;
  if (!PyArg_ParseTuple(args, "sii|ii", &fn, &k, &m, &imax, &a)) return NULL;
  
  DumpRates(fn, k, m, imax, a);
  
  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PModifyRates(PyObject *self, PyObject *args) {
  char *fn;

  if (scrm_file) {
    SCRMStatement("ModifyRates", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  
  ModifyRates(fn);
  
  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PIonDensity(PyObject *self, PyObject *args) {
  char *fn;
  double d;
  int k;

  if (!PyArg_ParseTuple(args, "si", &fn, &k)) return NULL;
  
  d = IonDensity(fn, k);
  
  return Py_BuildValue("d", d);
}

static PyObject *PIonRadiation(PyObject *self, PyObject *args) {
  char *fn;
  double d;
  int k, m;

  m = 0;
  if (!PyArg_ParseTuple(args, "si|i", &fn, &k, &m)) return NULL;
  
  d = IonRadiation(fn, k, m);
  
  return Py_BuildValue("d", d);
}

static PyObject *PSetUTA(PyObject *self, PyObject *args) {
  int m, mci;

  if (scrm_file) {
    SCRMStatement("SetUTA", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  mci = 1;
  if (!PyArg_ParseTuple(args, "i|i", &m, &mci)) return NULL;
  
  SetUTA(m, mci);

  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PDRSuppression(PyObject *self, PyObject *args) {
  char *fn;
  int nmax;
  double z;

  if (scrm_file) {
    SCRMStatement("DRSuppression", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "sdi", &fn, &z, &nmax)) return NULL;
  
  DRSuppression(fn, z, nmax);
      
  Py_INCREF(Py_None);
  return Py_None;
}
  
static PyObject *PNormalizeMode(PyObject *self, PyObject *args) {
  int m;

  if (scrm_file) {
    SCRMStatement("NormalizeMode", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  
  NormalizeMode(m);
      
  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetBornFormFactor(PyObject *self, PyObject *args) {
  double te;
  char *fn;

  if (scrm_file) {
    SCRMStatement("SetBornFormFactor", args, NULL);
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

  if (scrm_file) {
    SCRMStatement("SetBornMass", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "d", &m)) return NULL;

  SetBornMass(m);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *PSetGamma3B(PyObject *self, PyObject *args) {
  double m;

  if (scrm_file) {
    SCRMStatement("SetGamma3B", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "d", &m)) return NULL;

  SetGamma3B(m);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PWallTime(PyObject *self, PyObject *args) {
  if (scrm_file) {
    SCRMStatement("WallTime", args, NULL);
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
  if (scrm_file) {
    SCRMStatement("InitializeMPI", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

#ifdef USE_MPI
  int n = -1;
  if (!(PyArg_ParseTuple(args, "|i", &n))) {
    return NULL;
  }
  InitializeMPI(n, 1);
#endif
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PMPIRank(PyObject *self, PyObject *args) {
  if (scrm_file) {
    SCRMStatement("MPIRank", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  int n, k;
  k = MPIRank(&n);
  return Py_BuildValue("[ii]", k, n);
}

static PyObject *PMemUsed(PyObject *self, PyObject *args) {
  if (scrm_file) {
    SCRMStatement("MemUsed", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  double m = msize();
  return Py_BuildValue("d", m);
}

static PyObject *PFinalizeMPI(PyObject *self, PyObject *args) {
  if (scrm_file) {
    SCRMStatement("FinalizeMPI", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

#if USE_MPI == 1
  FinalizeMPI();
#endif
  
  Py_INCREF(Py_None);
  return Py_None;
}

static struct PyMethodDef crm_methods[] = {
  {"Print", PPrint, METH_VARARGS}, 
  {"SetUTA", PSetUTA, METH_VARARGS}, 
  {"CloseSCRM", PCloseSCRM, METH_VARARGS},
  {"ConvertToSCRM", PConvertToSCRM, METH_VARARGS},
  {"CheckEndian", PCheckEndian, METH_VARARGS},
  {"DRSuppression", PDRSuppression, METH_VARARGS},
  {"EleDist", PEleDist, METH_VARARGS},
  {"IonDensity", PIonDensity, METH_VARARGS},
  {"IonRadiation", PIonRadiation, METH_VARARGS},
  {"PhoDist", PPhoDist, METH_VARARGS},
  {"SetEleDist", PSetEleDist, METH_VARARGS},
  {"SetPhoDist", PSetPhoDist, METH_VARARGS},
  {"SetNumSingleBlocks", PSetNumSingleBlocks, METH_VARARGS},
  {"SetExtrapolate", PSetExtrapolate, METH_VARARGS},
  {"SetEMinAI", PSetEMinAI, METH_VARARGS},
  {"SetInnerAuger", PSetInnerAuger, METH_VARARGS},
  {"SetEleDensity", PSetEleDensity, METH_VARARGS},
  {"SetPhoDensity", PSetPhoDensity, METH_VARARGS},
  {"SetCascade", PSetCascade, METH_VARARGS},
  {"SetIteration", PSetIteration, METH_VARARGS},
  {"SetRateAccuracy", PSetRateAccuracy, METH_VARARGS},
  {"SetBlocks", PSetBlocks, METH_VARARGS},
  {"RateTable", PRateTable, METH_VARARGS},
  {"AddIon", PAddIon, METH_VARARGS},
  {"SetCERates", PSetCERates, METH_VARARGS},
  {"SetTRRates", PSetTRRates, METH_VARARGS},
  {"SetCIRates", PSetCIRates, METH_VARARGS},
  {"SetRRRates", PSetRRRates, METH_VARARGS},
  {"SetAIRates", PSetAIRates, METH_VARARGS},
  {"SetAIRatesInner", PSetAIRatesInner, METH_VARARGS},
  {"SetAbund", PSetAbund, METH_VARARGS},
  {"SetRateMultiplier", PSetRateMultiplier, METH_VARARGS},
  {"InitBlocks", PInitBlocks, METH_VARARGS},
  {"LevelPopulation", PLevelPopulation, METH_VARARGS},
  {"Cascade", PCascade, METH_VARARGS},
  {"SpecTable", PSpecTable, METH_VARARGS},
  {"SelectLines", PSelectLines, METH_VARARGS},
  {"PlotSpec", PPlotSpec, METH_VARARGS},
  {"TabNLTE", PTabNLTE, METH_VARARGS},
  {"PrintTable", PPrintTable, METH_VARARGS}, 
  {"ReinitCRM", PReinitCRM, METH_VARARGS},
  {"Bremss", PBremss, METH_VARARGS},
  {"PhFit", PPhFit, METH_VARARGS},
  {"EPhFit", PEPhFit, METH_VARARGS},
  {"RRFit", PRRFit, METH_VARARGS},
  {"DRFit", PDRFit, METH_VARARGS},
  {"NRRFit", PNRRFit, METH_VARARGS},
  {"NDRFit", PNDRFit, METH_VARARGS},
  {"CBeli", PCBeli, METH_VARARGS},
  {"CFit", PCFit, METH_VARARGS},
  {"ColFit", PColFit, METH_VARARGS},
  {"CColFit", PCColFit, METH_VARARGS},
  {"EColFit", PEColFit, METH_VARARGS},
  {"Ionis", PIonis, METH_VARARGS},
  {"RBeli", PRBeli, METH_VARARGS},
  {"EBeli", PEBeli, METH_VARARGS},
  {"Recomb", PRecomb, METH_VARARGS},
  {"RRRateH", PRRRateH, METH_VARARGS},
  {"TwoPhoton", PTwoPhoton, METH_VARARGS},
  {"FracAbund", PFracAbund, METH_VARARGS},
  {"MaxAbund", PMaxAbund, METH_VARARGS},
  {"DRBranch", PDRBranch, METH_VARARGS},
  {"DRStrength", PDRStrength, METH_VARARGS},
  {"DumpRates", PDumpRates, METH_VARARGS},
  {"ModifyRates", PModifyRates, METH_VARARGS},
  {"RydBranch", PRydBranch, METH_VARARGS},
  {"NormalizeMode", PNormalizeMode, METH_VARARGS},
  {"SetBornFormFactor", PSetBornFormFactor, METH_VARARGS},
  {"SetBornMass", PSetBornMass, METH_VARARGS},
  {"SetGamma3B", PSetGamma3B, METH_VARARGS},
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
  "crm",
  NULL,
  -1,
  crm_methods,
};
#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_crm(void){

#else
#define INITERROR return

void
initcrm(void) {
#endif

  PyObject *m, *d;

  #if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
  #else
  m = Py_InitModule("crm", crm_methods);
  #endif
  d = PyModule_GetDict(m);
  ErrorObject = Py_BuildValue("s", "crm.error");
  PyDict_SetItemString(d, "error", ErrorObject);
  InitCRM();
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module crm");

#if PY_MAJOR_VERSION >= 3
  return m;
#else
  return;
#endif
}
