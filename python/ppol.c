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


static char *rcsid="$Id: ppol.c,v 1.10 2003/11/25 21:00:33 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "Python.h"
#include "polarization.h"

#if PY_MAJOR_VERSION >= 3
  #define PyUnicode_AsString(x) PyBytes_AsString(PyUnicode_AsEncodedString((x), "utf-8", "strict"))
#else
  #define PyUnicode_AsString PyString_AsString
  #define PyLong_AsLong PyInt_AsLong
  #define PyLong_Check PyInt_Check
#endif

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
  s1 = PyUnicode_AsString(sargs);
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
      s2 = PyUnicode_AsString(q);
      fprintf(spol_file, "%s=", s2);
      q = PyTuple_GetItem(p, 1);
      kvar = PyObject_Str(q);
      s2 = PyUnicode_AsString(kvar);
      if (PyUnicode_Check(q)) {
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
    SPOLStatement("SetMIteration", args, NULL);
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

static PyObject *POrientation(PyObject *self, PyObject *args) {
  double e;
  char *fn;
  
  if (spol_file) {
    SPOLStatement("Orientation", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  fn = NULL;
  e = 0.0;
  if (!PyArg_ParseTuple(args, "|ds", &e, &fn)) return NULL;
  
  if (Orientation(fn, e) < 0) return NULL;
  
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PPolarizationTable(PyObject *self, PyObject *args) {
  PyObject *p, *q;
  char *fn;
  char *ifn;
  char **sc;
  int n, i;
  
  if (spol_file) {
    SPOLStatement("PolarizationTable", args, NULL);
    Py_INCREF(Py_None);
    return Py_None;
  }

  p = NULL;
  if (!PyArg_ParseTuple(args, "s|O", &fn, &p)) return NULL;

  ifn = NULL;
  n = 0;
  sc = NULL;
  if (p) {
    if (PyUnicode_Check(p)) {
      ifn = PyUnicode_AsString(p);
    } else if (PyList_Check(p)) {
      n = PyList_Size(p);
      if (n > 0) {
	sc = malloc(sizeof(char *)*n);
	for (i = 0; i < n; i++) {
	  q = PyList_GetItem(p, i);
	  if (!PyUnicode_Check(q)) {
	    free(sc);
	    return NULL;
	  }
	  sc[i] = PyUnicode_AsString(q);
	}
      }
    }
  }
  
  i = PolarizationTable(fn, ifn, n, sc);
  if (n > 0) {
    free(sc);
  }

  if (i < 0) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
} 

static PyObject *PWallTime(PyObject *self, PyObject *args) {
  if (spol_file) {
    SPOLStatement("WallTime", args, NULL);
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

static struct PyMethodDef pol_methods[] = {
  {"Print", PPrint, METH_VARARGS},
  {"ConvertToSPOL", PConvertToSPOL, METH_VARARGS},
  {"CloseSPOL", PCloseSPOL, METH_VARARGS},
  {"Orientation", POrientation, METH_VARARGS},
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
  {"WallTime", PWallTime, METH_VARARGS}, 
  {NULL, NULL, METH_VARARGS}
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "pol",
  NULL,
  -1,
  pol_methods,
};
#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_pol(void){

#else
#define INITERROR return

void
initpol(void) {
#endif

  PyObject *m, *d;

  #if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
  #else
  m = Py_InitModule("pol", pol_methods);
  #endif
  d = PyModule_GetDict(m);
  ErrorObject = Py_BuildValue("s", "pol.error");
  PyDict_SetItemString(d, "error", ErrorObject);
  InitPolarization();

  if (PyErr_Occurred())
    Py_FatalError("can't initialize module pol");

#if PY_MAJOR_VERSION >= 3
  return m;
#else
  return;
#endif
}
