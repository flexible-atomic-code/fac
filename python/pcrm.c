
static char *rcsid="$Id: pcrm.c,v 1.28 2003/06/09 14:49:48 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "Python.h"
#include <stdio.h>
#include <string.h>

#include "crm.h"

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
  s1 = PyString_AsString(sargs);
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
      s2 = PyString_AsString(q);
      fprintf(scrm_file, "%s=", s2);
      q = PyTuple_GetItem(p, 1);
      kvar = PyObject_Str(q);
      s2 = PyString_AsString(kvar);
      if (PyString_Check(q)) {
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
  FILE *f;
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
    f = fopen(fn, "rb");
    if (f == NULL) {
      printf("Cannot open file %s\n", fn);
      return NULL;
    }
    ReadFHeader(f, &fh, &swp);
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

static PyObject *PSetEleDist(PyObject *self, PyObject *args) {
  PyObject *p;
  int i, n, k, np;
  double *par = NULL;

  if (scrm_file) {
    SCRMStatement("SetEleDist", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n < 1) return NULL;
  p = PyTuple_GetItem(args, 0);
  i = PyInt_AsLong(p);
  np = n-1;
  if (np > 0) 
    par = (double *) malloc(sizeof(double)*np);
  for (k = 1; k < n; k++) {
    p = PyTuple_GetItem(args, k);
    par[k-1] = PyFloat_AsDouble(p);
  }
  if (SetEleDist(i, np, par) < 0) return NULL;
  if (np > 0) free(par);
 
  Py_INCREF(Py_None);
  return Py_None;
} 
 
static PyObject *PSetPhoDist(PyObject *self, PyObject *args) {
  PyObject *p;
  int i, n, k, np;
  double *par = NULL;

  if (scrm_file) {
    SCRMStatement("SetPhoDist", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  n = PyTuple_Size(args);
  if (n < 1) return NULL;
  p = PyTuple_GetItem(args, 0);
  i = PyInt_AsLong(p);
  np = n-1;
  if (np > 0) 
    par = (double *) malloc(sizeof(double)*np);
  for (k = 1; k < n; k++) {
    p = PyTuple_GetItem(args, k);
    par[k-1] = PyFloat_AsDouble(p);
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

static PyObject *PPrintTable(PyObject *self, PyObject *args) { 
  char *fn1, *fn2;
  int v;
  
  if (scrm_file) {
    SCRMStatement("PrintTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  v = 0;
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
  int i, nc;
  char **sc;
    
  if (scrm_file) {
    SCRMStatement("RateTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  nc = PyTuple_Size(args);
  if (nc == 0) return NULL;
  sc = (char **) malloc(sizeof(char *)*nc);
  for (i = 0; i < nc; i++) {
    sc[i] = PyString_AsString(PyTuple_GetItem(args, i));
  }
  RateTable(sc[0], nc-1, &(sc[1]));

  free(sc);
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

  c = PyString_AsString(s);
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

static PyObject *PEPhFit(PyObject *self, PyObject *args) {
  PyObject *s;
  int z, nele, ks;
  double e;
  
  if (!PyArg_ParseTuple(args, "iiO", &z, &nele, &s)) return NULL;
  ks = -1;
  if (PyString_Check(s)) {
    ks = ShellIndexFromString(s);
  } else if (PyInt_Check(s)) {
    ks = PyInt_AsLong(s);
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
  if (PyString_Check(s)) {
    ks = ShellIndexFromString(s);
  } else if (PyInt_Check(s)) {
    ks = PyInt_AsLong(s);
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
  } else if (PyString_Check(s)) {
    ks = ShellIndexFromString(s);
  } else if (PyInt_Check(s)) {
    ks = PyInt_AsLong(s);
  } else {
    return NULL;
  }
  if (ks < 0) return NULL;
  r = ColFit(z, nele, ks, t, &a, &d);
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
  } else if (PyString_Check(s)) {
    ks = ShellIndexFromString(s);
  } else if (PyInt_Check(s)) {
    ks = PyInt_AsLong(s);
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
  
  DRBranch();

  Py_INCREF(Py_None);  
  return Py_None;
}

static PyObject *PDRStrength(PyObject *self, PyObject *args) {
  int n, i, m;
  char *s;
  
  i = 0;
  m = 0;
  if (!PyArg_ParseTuple(args, "si|ii", &s, &n, &m, &i)) 
    return NULL;

  DRStrength(s, n, m, i);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static struct PyMethodDef crm_methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"CloseSCRM", PCloseSCRM, METH_VARARGS},
  {"ConvertToSCRM", PConvertToSCRM, METH_VARARGS},
  {"CheckEndian", PCheckEndian, METH_VARARGS},
  {"SetEleDist", PSetEleDist, METH_VARARGS},
  {"SetPhoDist", PSetPhoDist, METH_VARARGS},
  {"SetNumSingleBlocks", PSetNumSingleBlocks, METH_VARARGS},
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
  {"SetAbund", PSetAbund, METH_VARARGS},
  {"InitBlocks", PInitBlocks, METH_VARARGS},
  {"LevelPopulation", PLevelPopulation, METH_VARARGS},
  {"Cascade", PCascade, METH_VARARGS},
  {"SpecTable", PSpecTable, METH_VARARGS},
  {"SelectLines", PSelectLines, METH_VARARGS},
  {"PlotSpec", PPlotSpec, METH_VARARGS},
  {"PrintTable", PPrintTable, METH_VARARGS}, 
  {"ReinitCRM", PReinitCRM, METH_VARARGS},
  {"PhFit", PPhFit, METH_VARARGS},
  {"EPhFit", PEPhFit, METH_VARARGS},
  {"RRFit", PRRFit, METH_VARARGS},
  {"DRFit", PDRFit, METH_VARARGS},
  {"NRRFit", PNRRFit, METH_VARARGS},
  {"NDRFit", PNDRFit, METH_VARARGS},
  {"CBeli", PCBeli, METH_VARARGS},
  {"CFit", PCFit, METH_VARARGS},
  {"ColFit", PColFit, METH_VARARGS},
  {"EColFit", PEColFit, METH_VARARGS},
  {"Ionis", PIonis, METH_VARARGS},
  {"RBeli", PRBeli, METH_VARARGS},
  {"EBeli", PEBeli, METH_VARARGS},
  {"Recomb", PRecomb, METH_VARARGS},
  {"RRRateH", PRRRateH, METH_VARARGS},
  {"FracAbund", PFracAbund, METH_VARARGS},
  {"MaxAbund", PMaxAbund, METH_VARARGS},
  {"DRBranch", PDRBranch, METH_VARARGS},
  {"DRStrength", PDRStrength, METH_VARARGS},
  {NULL, NULL}
};

void initcrm(void) {
  PyObject *m, *d;
  
  m = Py_InitModule("crm", crm_methods);
  
  d = PyModule_GetDict(m);
  ErrorObject = Py_BuildValue("s", "crm.error");
  PyDict_SetItemString(d, "error", ErrorObject);

  InitCRM();

  if (PyErr_Occurred()) 
    Py_FatalError("can't initialize module crm");
}
