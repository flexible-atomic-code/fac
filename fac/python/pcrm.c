static char *rcsid="$Id: pcrm.c,v 1.8 2002/02/12 20:32:17 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "Python.h"
#include <stdio.h>
#include <string.h>

#include "interpolation.h"
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
  int i;

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
    fread(&fh, sizeof(F_HEADER), 1, f);
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
  SetEleDist(i, np, par);
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
  SetPhoDist(i, np, par);
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

  if (scrm_file) {
    SCRMStatement("SpecTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  smin = EPS10;
  if (!PyArg_ParseTuple(args, "s|d", &fn, &smin)) return NULL;
  
  SpecTable(fn, smin);
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

static PyObject *PFreeMemENTable(PyObject *self, PyObject *args) {

  if (scrm_file) {
    SCRMStatement("FreeMemENTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  FreeMemENTable();

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *PLevelName(PyObject *self, PyObject *args) { 
  char *fn;
  int i, k;
  char cname[LNCOMPLEX];
  char sname[LSNAME];
  char name[LNAME];
  
  if (!PyArg_ParseTuple(args, "si", &fn, &i)) return NULL;

  k = LevelName(fn, i, cname, sname, name);
  if (k < 0) return NULL;
  if (k > 0) {
    cname[0] = '\0';
    sname[0] = '\0';
    name[0] = '\0';
  }
  return Py_BuildValue("(sss)", cname, sname, name);
}
  
static PyObject *PMemENTable(PyObject *self, PyObject *args) { 
  char *fn;
  
  if (scrm_file) {
    SCRMStatement("MemENTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  MemENTable(fn);

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
  char *fn;
    
  if (scrm_file) {
    SCRMStatement("RateTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;

  RateTable(fn);
 
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

static PyObject *PSpline(PyObject *self, PyObject *args) {
  PyObject *px, *py, *py2;
  double *x, *y, *y2, dy1, dy2;
  int n, i;

  if (scrm_file) {
    printf("SCRM does not support Spline\n");
    return NULL;
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

  if (scrm_file) {
    printf("SCRM does not support Splint\n");
    return NULL;
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

static PyObject *PSelectLines(PyObject *self, PyObject *args) {
  char *ifn, *ofn;
  double emin, emax, fmin;
  int nele, type;

  if (scrm_file) {
    SCRMStatement("SelectLines", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  fmin = EPS6;
  if (!PyArg_ParseTuple(args, "ssiidd|d", 
			&ifn, &ofn, &nele, &type, &emin, &emax, &fmin)) 
    return NULL;
  
  SelectLines(ifn, ofn, nele, type, emin, emax, fmin);
  
  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PIonis(PyObject *self, PyObject *args) {
  int z, nele;
  double t, total, a, d;

  if (!PyArg_ParseTuple(args, "iid", &z, &nele, &t)) return NULL;
  
  total = Ionis(z, nele, t, &a, &d);
  
  return Py_BuildValue("(ddd)", total, a, d);
}

static PyObject *PRecomb(PyObject *self, PyObject *args) {
  int z, nele;
  double t, total, r, d;

  if (!PyArg_ParseTuple(args, "iid", &z, &nele, &t)) return NULL;
  
  total = Recomb(z, nele, t, &r, &d);
  
  return Py_BuildValue("(ddd)", total, r, d);
}

static PyObject *PFracAbund(PyObject *self, PyObject *args) {
  PyObject *pa;
  int z, i;
  double t, *a;

  if (!PyArg_ParseTuple(args, "id", &z, &t)) return NULL;
  
  a = (double *) malloc(sizeof(double)*(z+1));
  FracAbund(z, t, a);
  pa = Py_BuildValue("[]");
  for (i = 0; i <= z; i++) {
    PyList_Append(pa, Py_BuildValue("d", a[i]));
  }
  free(a);

  return pa;
}

static PyObject *PMaxAbund(PyObject *self, PyObject *args) {
  PyObject *pa;
  int i, z, nele;
  double eps, *a, tmax;

  eps = 1E-4;
  if (!PyArg_ParseTuple(args, "ii|d", &z, &nele, &eps)) return NULL;
  
  a = (double *) malloc(sizeof(double)*(z+1));
  tmax = MaxAbund(z, nele, a, eps);
  pa = Py_BuildValue("[]");
  for (i = 0; i <= z; i++) {
    PyList_Append(pa, Py_BuildValue("d", a[i]));
  }
  free(a);
  return Py_BuildValue("(dO)", tmax, pa);
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
  {"FreeMemENTable", PFreeMemENTable, METH_VARARGS},
  {"LevelName", PLevelName, METH_VARARGS},
  {"MemENTable", PMemENTable, METH_VARARGS},
  {"PrintTable", PPrintTable, METH_VARARGS}, 
  {"ReinitCRM", PReinitCRM, METH_VARARGS},
  {"Spline", PSpline, METH_VARARGS},
  {"Splint", PSplint, METH_VARARGS},
  {"Ionis", PIonis, METH_VARARGS},
  {"Recomb", PRecomb, METH_VARARGS},
  {"FracAbund", PFracAbund, METH_VARARGS},
  {"MaxAbund", PMaxAbund, METH_VARARGS},
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
