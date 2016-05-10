#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "poibinmodule.c"

// HELPERS
int max(int x, int y){
  if (x > y) return x;
  else return y;
}

int min(int x, int y){
  if (x < y) return x;
  else return y;
}

// Joint probability (see http://goo.gl/QLDTUF)
double joint_mass(int n, int z, int x, int y, double *p_x, double *p_y, double ****cache){
    int i;
    // Base cases (top of page 7)
    if (0 > min(min(z, x), y) || z > min(x, y) || max(x, y) > n){
        return 0.0;
    } else if (n == 0 && z == 0 && x == 0 && y == 0){
        return 1.0;
    // Check the cache
    } else if (cache[n][z][x][y] != -1.0){
        return cache[n][z][x][y];
    // Recursive case
    } else{
        // Equation at bottom of page 6
		i = n-1; // subtract one because of zero-based indexing
		cache[n][z][x][y] =  p_x[i]      * p_y[i]      * joint_mass(n-1, z-1, x-1, y-1, p_x, p_y, cache) +
			                 p_x[i]      * (1.-p_y[i]) * joint_mass(n-1, z, x-1, y, p_x, p_y, cache) +
			                (1.-p_x[i]) * p_y[i]      * joint_mass(n-1, z, x, y-1, p_x, p_y, cache) +
			                (1.-p_x[i]) * (1.-p_y[i]) * joint_mass(n-1, z, x, y, p_x, p_y, cache);
    }
    return cache[n][z][x][y];
}

PyObject *py_conditional(PyObject *self, PyObject *args){
    // Parameters
    int i, j, i2, j2, N, x, y, *zs, num_zs;
    double *p_x, *p_y, joint_marginal, mass, ****cache;
    PyObject *py_zs, *py_p_x, *py_p_y, *results;

    // Parse Python arguments
    if (! PyArg_ParseTuple( args, "iO!iiO!O!", &N, &PyList_Type, &py_zs,
                            &x, &y, &PyList_Type, &py_p_x, &PyList_Type,
                            &py_p_y)){
        return NULL;
    }

    p_x = malloc(sizeof(double) * N);
    p_y = malloc(sizeof(double) * N);
    for (i = 0; i < N; i ++){
        p_x[i] = (double) PyFloat_AsDouble(PyList_GetItem(py_p_x, i));
        p_y[i] = (double) PyFloat_AsDouble(PyList_GetItem(py_p_y, i));
    }

    num_zs = PyList_Size(py_zs);
    zs  = malloc(sizeof(int) * num_zs);
    for (i = 0; i < num_zs; i++){
      zs[i]  = (int) PyLong_AsLong(PyList_GetItem(py_zs, i));
    }

    // Initialize the cache
    cache = malloc(sizeof(double ***) * (N+1));
    for (i = 0; i < N+1; i++){
        cache[i] = malloc(sizeof(double **) * (N+1));
        for (j = 0; j < N+1; j++){
            cache[i][j] = malloc(sizeof(double *) * (x+1));
            for ( i2 = 0; i2 < x+1; i2++){
                cache[i][j][i2] = malloc(sizeof(double) * (y+1));
                for (j2 = 0; j2 < y+1; j2++){
                    cache[i][j][i2][j2] = (double) -1.0;
                }
            }
        }
    }

    // Create a list of the masses of each
    results        = PyList_New(num_zs);
    joint_marginal = pmf(x, N, p_x) * pmf(y, N, p_y);
    for (i = 0; i < num_zs; i++){
        mass = joint_mass( N, zs[i], x, y, p_x, p_y, cache );
        PyList_SetItem(results, i, Py_BuildValue("f", mass / joint_marginal));
    }

    // Free memory
    free(zs);
    free(p_x);
    free(p_y);
    for (i = 0; i < N+1; i++){
        for (j = 0; j < N+1; j++){
            for (i2 = 0; i2 < x+1; i2++){
                free(cache[i][j][i2]);
            }
            free(cache[i][j]);
        }
        free(cache[i]);
    }
    free(cache);

    return results;
}


////////////////////////////////////////////////////////////////////////////////
// TEST FOR TRIPLES
////////////////////////////////////////////////////////////////////////////////
double P(int n, int t, int w, int x, int y, double **p, double *****cache){
    int i;

    // Base cases
    if (t < 0 || w < 0 || x < 0 || y < 0 || t > w + x + y || w > n || x > n || y > n){
        return 0.0;
    } else if (t == 0 && n == 0 && w == 0 && x == 0 && y == 0){
        return 1.0;
    // Use the cache when possible
    } else if (cache[n][t][w][x][y] != -1.0){
        return cache[n][t][w][x][y];
    // Recurse
    } else {
        i = n - 1;
		cache[n][t][w][x][y] = (1. - p[0][i]) * (1. - p[1][i]) * (1. - p[2][i]) * P(n-1, t,   w,   x,   y,   p, cache) +
							   (1. - p[0][i]) * (1. - p[1][i]) * p[2][i]        * P(n-1, t-1, w,   x,   y-1, p, cache) +
							   (1. - p[0][i]) * p[1][i]        * (1. - p[2][i]) * P(n-1, t-1, w,   x-1, y,   p, cache) +
							   (1. - p[0][i]) * p[1][i]        * p[2][i]        * P(n-1, t,   w,   x-1, y-1, p, cache) +
							   p[0][i]        * (1. - p[1][i]) * (1. - p[2][i]) * P(n-1, t-1, w-1, x,   y,   p, cache) +
							   p[0][i]        * (1. - p[1][i]) * p[2][i]        * P(n-1, t,   w-1, x,   y-1, p, cache) +
							   p[0][i]        * p[1][i]        * (1. - p[2][i]) * P(n-1, t,   w-1, x-1, y,   p, cache) +
							   p[0][i]        * p[1][i]        * p[2][i]        * P(n-1, t,   w-1, x-1, y-1, p, cache);
    }
	return cache[n][t][w][x][y];
}

PyObject *triple_exact_test(PyObject *self, PyObject *args){
    // Parameters
    int i, j, i2, j2, i3, N, w, x, y, t, T;
    double **p, marginals, joint, result, *****cache;
    PyObject *py_p;

    // Parse Python arguments
    if (! PyArg_ParseTuple( args, "iiiiiO!", &N, &T, &w, &x, &y, &PyList_Type, &py_p)){
        return NULL;
    }

    // Load the success probability matrix
    p = malloc(sizeof(double *)  * 3);
    for (i = 0; i < 3; i++){
        p[i] = malloc(sizeof(double) * N);
        for (j = 0; j < N; j ++){
            p[i][j] = (double) PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_p, i), j));
        }
    }

    // Initialize the cache
    cache = malloc(sizeof(double ****) * (N+1));
    for (i = 0; i < N+1; i++){
        cache[i] = malloc(sizeof(double ***) * (N+1));
        for (j = 0; j < N+1; j++){
            cache[i][j] = malloc(sizeof(double **) * (w+1));
            for (i2 = 0; i2 < w+1; i2++){
                cache[i][j][i2] = malloc(sizeof(double *) * (x+1));
                for (j2 = 0; j2 < x+1; j2++){
                    cache[i][j][i2][j2] = malloc(sizeof(double) * (y+1));
                    for (i3 = 0; i3 < y+1; i3++){
                        cache[i][j][i2][j2][i3] = (double) -1.0;
                    }
                }
            }
        }
    }

    // Create a list of the masses of each
    result = 0.0;
    marginals = pmf(w, N, p[0]) * pmf(x, N, p[1]) * pmf(y, N, p[2]);
    for (t = T; t < min(N, w + x + y)+1; t++){
        joint = P( N, t, w, x, y, p, cache );
        result += joint/marginals;
    }

    // Free memory
    for (i = 0; i < 3; i++) free(p[i]);
    free(p);
    for (i = 0; i < N+1; i++){
        for (j = 0; j < N+1; j++){
            for (i2 = 0; i2 < w+1; i2++){
                for (j2 = 0; j2 < x+1; j2++){
                    free(cache[i][j][i2][j2]);
                }
                free(cache[i][j][i2]);
            }
            free(cache[i][j]);
        }
        free(cache[i]);
    }
    free(cache);

    return Py_BuildValue("f", result);
}

////////////////////////////////////////////////////////////////////////////////
// Register the functions we want to be accessible from Python
PyMethodDef weightedEnrichmentMethods[] = {
    {"conditional", py_conditional, METH_VARARGS, "Weighted enrichment test conditional PMF for pairs"},
    {"triple_exact_test", triple_exact_test, METH_VARARGS, "Weighted enrichment test for triples"}
};

// Note that the suffix of init has to match the name of the module,
// both here and in the setup.py file
PyMODINIT_FUNC initwext_exact_test(void) {
    PyObject *m = Py_InitModule("wext_exact_test", weightedEnrichmentMethods);
    if (m == NULL) {
        return;
    }
}
