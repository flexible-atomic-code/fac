static char *rcsid="$Id: util.c,v 1.1 2002/04/25 16:22:29 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "Python.h"
#include <stdio.h>

#include "interpolation.h"

void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);

static PyObject *ErrorObject;
#define onError(message) {PyErr_SetString(ErrorObject, message);}

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

static PyObject *PUVIP3P(PyObject *self, PyObject *args) {
  PyObject *px, *py, *px0, *py0;
  double *x, *y, *x0, *y0;
  int n, m, np, i, f;  

  if (!PyArg_ParseTuple(args, "OOO", &px, &py, &px0)) return NULL;
  if (!PyList_Check(px)) return NULL;
  if (!PyList_Check(py)) return NULL;
  n = PyList_Size(px);
  np = 3;
  if (PyList_Size(py) != n) return NULL;
  if (!PyList_Check(px0)) {
    if (!PyFloat_Check(px0)) return NULL;
    x0 = (double *) malloc(sizeof(double));
    y0 = (double *) malloc(sizeof(double));
    x0[0] = PyFloat_AsDouble(px0);
    m = 1;
    f = 1;
  } else {
    m = PyList_Size(px0);
    x0 = (double *) malloc(sizeof(double)*m);
    y0 = (double *) malloc(sizeof(double)*m);
    for (i = 0; i < m; i++) {
      x0[i] = PyFloat_AsDouble(PyList_GetItem(px0, i));
    }
    f = 0;
  }
  x = (double *) malloc(sizeof(double)*n);
  y = (double *) malloc(sizeof(double)*n);
  for (i = 0; i < n; i++) {
    x[i] = PyFloat_AsDouble(PyList_GetItem(px, i));
    y[i] = PyFloat_AsDouble(PyList_GetItem(py, i));
  }
  uvip3p_(&np, &n, x, y, &m, x0, y0);
  if (f) {
    py0 = Py_BuildValue("d", y0[0]);
  } else {
    py0 = Py_BuildValue("[]");
    for (i = 0; i < m; i++) {
      PyList_Append(py0, Py_BuildValue("d", y0[i]));
    }
  }
  free(x);
  free(y);
  free(x0);
  free(y0);

  return py0;
}

static struct PyMethodDef util_methods[] = {
  {"Spline", PSpline, METH_VARARGS},
  {"Splint", PSplint, METH_VARARGS},
  {"UVIP3P", PUVIP3P, METH_VARARGS},
  {NULL, NULL}
};

void initutil(void) {
  PyObject *m, *d;
  m = Py_InitModule("util", util_methods);
  
  d = PyModule_GetDict(m);
  ErrorObject = Py_BuildValue("s", "util.error");
  PyDict_SetItemString(d, "error", ErrorObject);

  if (PyErr_Occurred()) 
    Py_FatalError("can't initialize module crm");
}

  
