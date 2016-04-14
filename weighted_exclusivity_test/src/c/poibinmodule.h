// Requirements
#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// Function declarations
double pmf_recursion(int k, int j, double *ps, double **cache);
double pmf(int k, int N, double *ps);
PyObject *py_pmf(PyObject *self, PyObject *args);
