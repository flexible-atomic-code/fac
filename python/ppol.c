
static char *rcsid="$Id: ppol.c,v 1.4 2003/08/05 16:25:59 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "Python.h"
#include "polarization.h"

static PyObject *ErrorObject;
#define onError(message) {PyErr_SetString(ErrorObject, message);}

static FILE *spol_file = NULL;

static void SPOLStatement(char *func, PyObject *args, PyObject *kargs) {
  int i, n, nargs;
  PyObject *sargs;
  PyObject *klist;
  PyObject *kvar;
  PyObject *p, *q;
  char *s1, *s2;
  
  fprintf(spol_file, "%s", func);
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
    fprintf(spol_file, "%c", s1[i]);
  }
  if (kargs) {
    klist = PyDict_Items(kargs);
    n = PyList_Size(klist);
    for (i = 0; i < n; i++) {
      if (nargs > 0 || i > 0) fprintf(spol_file, ", ");
      p = PyList_GetItem(klist, i);
      q = PyTuple_GetItem(p, 0);
      s2 = PyString_AsString(q);
      fprintf(spol_file, "%s=", s2);
      q = PyTuple_GetItem(p, 1);
      kvar = PyObject_Str(q);
      s2 = PyString_AsString(kvar);
      if (PyString_Check(q)) {
	fprintf(spol_file, "'%s'", s2);
      } else {
	fprintf(spol_file, "%s", s2);
      }
      Py_XDECREF(kvar);
    }
    Py_XDECREF(klist);
  }

  fprintf(spol_file, ")\n");
    
  Py_XDECREF(sargs);

  return;
}

static PyObject *PConvertToSPOL(PyObject *self, PyObject *args) {
  char *fn;

  spol_file = NULL;
  if (!PyArg_ParseTuple(args, "|s", &fn)) return NULL;

  if (fn) {
    spol_file = fopen(fn, "w");
    if (spol_file == NULL) return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}    
  
static PyObject *PCloseSPOL(PyObject *self, PyObject *args) {

  fclose(spol_file);
  spol_file = NULL;

  Py_INCREF(Py_None);
  return Py_None;
}  

static PyObject *PPrint(PyObject *self, PyObject *args) {
  PyObject *p, *q;
  char *s;
  int i, n;

  if (spol_file) {
    SPOLStatement("Print", args, NULL);
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

static PyObject *PSetMaxLevels(PyObject *self, PyObject *args) {
  int m;

  if (spol_file) {
    SPOLStatement("SetMaxLevels", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (!PyArg_ParseTuple(args, "i", &m)) return NULL;
  
  SetMaxLevels(m);

  Py_INCREF(Py_None);
  return Py_None;
}   

static PyObject *PSetMIteration(PyObject *self, PyObject *args) {
  int m;
  double a;

  if (spol_file) {
    SPOLStatement("SetIteration", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  m = 0;
  if (!PyArg_ParseTuple(args, "d|i", &a, &m)) return NULL;
  
  SetMIteration(a, m);

  Py_INCREF(Py_None);
  return Py_None;
}   

static PyObject *PSetIDR(PyObject *self, PyObject *args) {
  int idr, ndr, i;
  double *pdr;
  PyObject *p, *q;

  if (spol_file) {
    SPOLStatement("SetIDR", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  p = NULL;
  if (!PyArg_ParseTuple(args, "i|O", &idr, &p)) return NULL;

  ndr = 0;
  pdr = NULL;
  if (p) {
    if (!PyList_Check(p) && !PyTuple_Check(p)) return NULL;
    ndr = PySequence_Length(p);
    pdr = (double *) malloc(sizeof(double)*ndr);
    for (i = 0; i < ndr; i++) {
      q = PySequence_GetItem(p, i);
      pdr[i] = PyFloat_AsDouble(q);
      Py_DECREF(q);
    }
  }
  
  if (SetIDR(idr, ndr, pdr) < 0) return NULL;

  Py_INCREF(Py_None);
  return Py_None;
} 
  
static PyObject *PSetEnergy(PyObject *self, PyObject *args) {
  double e, es;

  if (spol_file) {
    SPOLStatement("SetEnergy", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  es = 0.0;
  if (!PyArg_ParseTuple(args, "d|d", &e, &es)) return NULL;

  SetEnergy(e, es);

  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject *PSetDensity(PyObject *self, PyObject *args) {
  double d;

  if (spol_file) {
    SPOLStatement("SetDensity", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "d", &d)) return NULL;

  SetDensity(d);

  Py_INCREF(Py_None);
  return Py_None;
}     
  
static PyObject *PSetMLevels(PyObject *self, PyObject *args) {
  char *efn, *tfn;
  
  if (spol_file) {
    SPOLStatement("SetMLevels", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "ss", &efn, &tfn)) return NULL;
  
  if (SetMLevels(efn, tfn) < 0) return NULL;
  
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetMCERates(PyObject *self, PyObject *args) {
  char *fn;
  
  if (spol_file) {
    SPOLStatement("SetMCERates", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  
  if (SetMCERates(fn) < 0) return NULL;
  
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PSetMAIRates(PyObject *self, PyObject *args) {
  char *fn;
  
  if (spol_file) {
    SPOLStatement("SetMAIRates", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  
  if (SetMAIRates(fn) < 0) return NULL;
  
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PPopulationTable(PyObject *self, PyObject *args) {
  char *fn;
  
  if (spol_file) {
    SPOLStatement("PopulationTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  
  if (PopulationTable(fn) < 0) return NULL;
  
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PPolarizationTable(PyObject *self, PyObject *args) {
  char *fn;
  
  if (spol_file) {
    SPOLStatement("PolarizationTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (!PyArg_ParseTuple(args, "s", &fn)) return NULL;
  
  if (PolarizationTable(fn) < 0) return NULL;
  
  Py_INCREF(Py_None);
  return Py_None;
} 

static struct PyMethodDef pol_methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"ConvertToSPOL", PConvertToSPOL, METH_VARARGS},
  {"CloseSPOL", PCloseSPOL, METH_VARARGS},
  {"SetIDR", PSetIDR, METH_VARARGS},
  {"SetMaxLevels", PSetMaxLevels, METH_VARARGS},
  {"SetMIteration", PSetMIteration, METH_VARARGS},
  {"SetEnergy", PSetEnergy, METH_VARARGS},
  {"SetDensity", PSetDensity, METH_VARARGS},
  {"SetMLevels", PSetMLevels, METH_VARARGS},
  {"SetMCERates", PSetMCERates, METH_VARARGS},
  {"SetMAIRates", PSetMAIRates, METH_VARARGS},
  {"PopulationTable", PPopulationTable, METH_VARARGS}, 
  {"PolarizationTable", PPolarizationTable, METH_VARARGS},  
  {"", NULL, METH_VARARGS}
};

void initpol(void) {
  PyObject *m, *d;
    
  m = Py_InitModule("pol", pol_methods);
  
  d = PyModule_GetDict(m);
  ErrorObject = Py_BuildValue("s", "pol.error");
  PyDict_SetItemString(d, "error", ErrorObject);

  InitPolarization();

  if (PyErr_Occurred()) 
    Py_FatalError("can't initialize module pol");
}
