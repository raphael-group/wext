#include "poibinmodule.h"

// Poisson-Binomial PMF

// Recursive function for PMF that includes cache
double pmf_recursion(int k, int j, double *ps, double **cache){
  // Boundary cases
  if (k == 0 && j == 0){
    return 1.0;
  } else if (k == -1 || k == j + 1) {
    return 0.0;
  // Check to see if this value has been cached
  } else if (cache[k][j] != -1.0){
    return cache[k][j];
  // Recursion
  } else {
    cache[k][j] = (1.-ps[j-1]) * pmf_recursion(k, j-1, ps, cache) +
                  ps[j-1]      * pmf_recursion(k-1, j-1, ps, cache);
    return cache[k][j];
  }
}

// Wrapper for PMF recursion that creates and uses its own cache
double pmf(int k, int N, double *ps){
    // Initialize a cache
    int i, j;
    double mass, **cache;
    cache = malloc(sizeof(double*) * N+1);
    for (i = 0; i < N+1; i++){
        cache[i] = malloc(sizeof(double) * N+1);
        for (j = 0; j < N+1; j++){
            cache[i][j] = (double) -1.0;
        }
    }

    // Compute the mass
    mass = pmf_recursion( k, N, ps, cache);

    // Clean up and return
    for (i = 0; i < 2; i++){
      free(cache[i]);
    }
    free(cache);

    return mass;
}

PyObject *py_pmf(PyObject *self, PyObject *args){
  // Parameters
  int i, k, N;
  double result, *ps;
  PyObject *py_ps;

  // Parse Python arguments
  if (! PyArg_ParseTuple( args, "iO!", &k, &PyList_Type, &py_ps )){
    return NULL;
  }

  N  = PyList_Size(py_ps);
  ps = malloc(sizeof(double) * N);
  for (i = 0; i < N; i ++){
    ps[i] = (double) PyFloat_AsDouble(PyList_GetItem(py_ps, i));
  }

  // Call the PMF
  result = pmf(k, N, ps);

  // Free memory
  free(ps);

  return Py_BuildValue("d", result);
}

// Register the functions we want to be accessible from Python
PyMethodDef poibinMethods[] = {
    {"pmf", py_pmf, METH_VARARGS, "Poisson-Binomial PMF"}
};

// Note that the suffix of init has to match the name of the module,
// both here and in the setup.py file
PyMODINIT_FUNC initcpoibin(void) {
    PyObject *m = Py_InitModule("cpoibin", poibinMethods);
    if (m == NULL) {
        return;
    }
}
