static char *rcsid="$Id: pcrm.c,v 1.1 2002/01/14 23:19:48 mfgu Exp $";
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

static PyObject *PSpecTable(PyObject *self, PyObject *args) {
  char *fn;

  if (scrm_file) {
    SCRMStatement("SpecTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  
  SpecTable(fn);
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PPlotSpec(PyObject *self, PyObject *args) {
  char *fn1, *fn2;
  double emin, emax, de;
  int type;

  if (scrm_file) {
    SCRMStatement("PlotSpec", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "ssiddd", 
			&fn1, &fn2, &type, &emin, &emax, &de))
    return NULL;
      
  PlotSpec(fn1, fn2, type, emin, emax, de);
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
    
  if (scrm_file) {
    SCRMStatement("ReinitCRM", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  ReinitCRM();
 
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
 
static struct PyMethodDef crm_methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"CloseSCRM", PCloseSCRM, METH_VARARGS},
  {"ConvertToSCRM", PConvertToSCRM, METH_VARARGS},
  {"SetEleDist", PSetEleDist, METH_VARARGS},
  {"SetPhoDist", PSetPhoDist, METH_VARARGS},
  {"SetNumSingleBlocks", PSetNumSingleBlocks, METH_VARARGS},
  {"SetEleDensity", PSetEleDensity, METH_VARARGS},
  {"SetPhoDensity", PSetPhoDensity, METH_VARARGS},
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
  {"InitBlocks", PInitBlocks, METH_VARARGS},
  {"LevelPopulation", PLevelPopulation, METH_VARARGS},
  {"SpecTable", PSpecTable, METH_VARARGS},
  {"PlotSpec", PPlotSpec, METH_VARARGS},
  {"FreeMemENTable", PFreeMemENTable, METH_VARARGS},
  {"MemENTable", PMemENTable, METH_VARARGS},
  {"PrintTable", PPrintTable, METH_VARARGS}, 
  {"ReinitCRM", PReinitCRM, METH_VARARGS},
  {"", NULL, METH_VARARGS}
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
