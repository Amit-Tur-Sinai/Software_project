# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "symnmf.h"



double** convertPyToC(PyObject* py_list, int N, int cols) {
    double** mat = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        PyObject* row = PyList_GetItem(py_list, i);
        mat[i] = (double*)malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++) {
            PyObject* item = PyList_GetItem(row, j);
            mat[i][j] = PyFloat_AsDouble(item);    
        }
    }
    return mat;
}

PyObject* convertCToPy(double** mat, int N, int cols) {
    PyObject* py_list = PyList_New(N);
    for (int i = 0; i < N; i++) {
        PyObject* py_row = PyList_New(cols);
        for (int j = 0; j < cols; j++) {
            PyObject* py_value = PyFloat_FromDouble(mat[i][j]);
            PyList_SetItem(py_row, j, py_value); 
        }
        /* Set the inner list in the outer list */
        PyList_SetItem(py_list, i, py_row); 
    }
    return py_list;
}


static PyObject* ddg_func(PyObject *self, PyObject *args){
    PyObject *X, *result;
    double **X_mat, **D_mat, **A_mat;
    int N, cols;

    if(!PyArg_ParseTuple(args, "O", &X)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    N = PyList_Size(X); /* get the number of data points in the X matrix */
    cols = PyList_Size( PyList_GetItem(X, 0) ); /* get the number of columns (the dimension) of each data point */
    X_mat = convertPyToC(X, N, cols); 
    A_mat = sym(X_mat, N, cols); /* calls the sym function from the symnmf.c */
    D_mat = ddg(A_mat, N); /* calls the ddg function from the symnmf.c */
    result = convertCToPy(D_mat, N, N);
    
    freeMat(X_mat, N);
    freeMat(A_mat, N);
    freeMat(D_mat, N);
    return result;
}

static PyObject* sym_func(PyObject *self, PyObject *args){
    PyObject *X, *result;
    double **X_mat, **A_mat;
    int N, cols;

    if(!PyArg_ParseTuple(args, "O", &X)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a PyObject* so it is used to signal that an error has occurred. */
    }

    N = PyList_Size(X); /* get the number of data points in the X matrix */
    cols = PyList_Size( PyList_GetItem(X, 0) ); /* get the number of columns (the dimension) of each data point */
    X_mat = convertPyToC(X, N, cols); 
    A_mat = sym(X_mat, N, cols); /* calls the sym function from the symnmf.c */
    result = convertCToPy(A_mat, N, N);
    
    freeMat(X_mat, N);
    freeMat(A_mat, N);
    return result;
}

static PyObject* norm_func(PyObject *self, PyObject *args){
    PyObject *X, *result;
    double **X_mat, **A_mat, **D_mat, **W_mat;
    int N, cols;

    if(!PyArg_ParseTuple(args, "O", &X)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    N = PyList_Size(X); /* get the number of data points in the X matrix */
    cols = PyList_Size( PyList_GetItem(X, 0) ); /* get the number of columns (the dimension) of each data point */
    X_mat = convertPyToC(X, N, cols); 
    A_mat = sym(X_mat, N, cols); /* calls the sym function from the symnmf.c*/
    D_mat = ddg(A_mat, N); /* calls the ddg function from the symnmf.c */
    W_mat = norm(A_mat, D_mat, N);
    result = convertCToPy(W_mat, N, N);
    
    freeMat(X_mat, N);
    freeMat(A_mat, N);
    freeMat(D_mat, N);
    freeMat(W_mat, N);
    return result;
}


static PyObject* symnmf_func(PyObject *self, PyObject *args){
    PyObject *H, *W, *result;
    int N, K;
    double **W_mat, **H_mat;

    if(!PyArg_ParseTuple(args, "OOii", &H, &W, &N, &K)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    
    W_mat = convertPyToC(W, N, N); 
    H_mat = convertPyToC(H, N, K);
    H_mat = symnmf(H_mat, W_mat, N, K);
    result = convertCToPy(H_mat, N, K);
    freeMat(H_mat, N);
    freeMat(W_mat, N);

    return result;
}

/* module's function table */
static PyMethodDef Symnmf_FunctionsTable[] = {
    {"symnmf",              
      (PyCFunction) symnmf_func, 
      METH_VARARGS,           
      PyDoc_STR("Perform full the symNMF as described and output H")}, 
    {"norm",              
      (PyCFunction) norm_func, 
      METH_VARARGS,          
      PyDoc_STR("Calculate and output the normalized similarity matrix as described")},
    {"ddg",              
      (PyCFunction) ddg_func, 
      METH_VARARGS,          
      PyDoc_STR("Calculate and output the Diagonal Degree Matrix as described")},
    {"sym",                   
      (PyCFunction) sym_func, 
      METH_VARARGS,         
      PyDoc_STR("Calculate and output the similarity matrix as described")},
    {NULL, NULL, 0, NULL} 
};

/* modules definition */
static struct PyModuleDef Symnmf_Module = {
    PyModuleDef_HEAD_INIT,
    "symnmf",  /* name of module exposed to Python */
    NULL, /* module documentation */
    -1,
    Symnmf_FunctionsTable
};

PyMODINIT_FUNC PyInit_symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&Symnmf_Module);
    if (!m) {
        return NULL;
    }
    return m;
}