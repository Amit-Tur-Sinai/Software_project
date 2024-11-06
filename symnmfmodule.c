# define PY_SSIZE_T_CLEAN
# include <Python.h>

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "symnmf.h"


static PyObject* ddg_func(PyObject *self, PyObject *args){
    PyObject *X;

    if(!PyArg_ParseTuple(args, "O", &X)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    return (ddg_implementation(X));
}

static PyObject* sym_func(PyObject *self, PyObject *args){
    PyObject *X;

    if(!PyArg_ParseTuple(args, "O", &X)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    return (sym_implementation(X));
}

static PyObject* norm_func(PyObject *self, PyObject *args){
    PyObject *X;

    if(!PyArg_ParseTuple(args, "O", &X)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    return (norm_implementation(X));
}


static PyObject* symnmf_func(PyObject *self, PyObject *args){
    PyObject *H, *W;
    int N,K;

    if(!PyArg_ParseTuple(args, "OOii", &H, &W, &N, &K)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    return (symnmf_implementation(W,H,n,K));
}

// module's function table
static PyMethodDef Symnmf_FunctionsTable[] = {
    {"symnmf",              
      (PyCFunction) symnmf_func, 
      METH_VARARGS,           
      PyDoc_STR("ADD----------------------")}, 
    {"norm",              
      (PyCFunction) norm_func, 
      METH_VARARGS,          
      PyDoc_STR("ADD----------------------")},
    {"ddg",              
      (PyCFunction) ddg_func, 
      METH_VARARGS,          
      PyDoc_STR("ADD----------------------")},
    {"sym",                   
      (PyCFunction) sym_func, 
      METH_VARARGS,         
      PyDoc_STR("ADD----------------------")},
    {NULL, NULL, 0, NULL} 
};

// modules definition
static struct PyModuleDef Symnmf_Module = {
    PyModuleDef_HEAD_INIT,
    "symnmf",     // name of module exposed to Python
    NULL, // module documentation
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