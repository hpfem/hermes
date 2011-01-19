# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

# Compile this with:
# $ cython --cplus python_engine.pyx
# and it will generate the python_engine.cpp and python_engine_api.h files,
# that are included with the project. Do this whenever you modify this file.

#-----------------------------------------------------------------------
# Common C++ <-> Python+NumPy conversion tools:

#ctypedef double complex cplx

cdef extern from *:
    ctypedef char* char_p       "char*"
    ctypedef char* const_char_p "const char*"

cdef extern from "utilities.h":
    void throw_exception(char *msg)

import sys
import traceback
from numpy cimport (ndarray, npy_intp, import_array, PyArray_SimpleNewFromData,
        NPY_INT, NPY_DOUBLE, NPY_COMPLEX128)
from libc.string cimport memcpy

# this is important to be called here, otherwise we can't use the NumPy C/API:
import_array()

cdef api object namespace_create():
    return {"verbose": False}

cdef api void namespace_push(object namespace, const_char_p name, object o):
    namespace.update({name: o})

cdef api void namespace_print(object namespace):
    print "-"*80
    print "namespace:"
    print namespace

cdef api object namespace_pull(object namespace, const_char_p name):
    return namespace.get(name)

cdef api object c2py_int(int i):
    return i

cdef api int py2c_int(object i):
    try:
        return i
    except:
        etype, value, tb = sys.exc_info()
        s = "".join(traceback.format_exception(etype, value, tb))
        s = "Exception raised in the Python code:\n" + s
        throw_exception(s)

cdef api object c2py_double(double i):
    return i

cdef api double py2c_double(object i):
    try:
        return i
    except:
        etype, value, tb = sys.exc_info()
        s = "".join(traceback.format_exception(etype, value, tb))
        s = "Exception raised in the Python code:\n" + s
        throw_exception(s)

cdef api object c2py_str(const_char_p s):
    return s

cdef api char* py2c_str(object s):
    try:
        return s
    except:
        etype, value, tb = sys.exc_info()
        s = "".join(traceback.format_exception(etype, value, tb))
        s = "Exception raised in the Python code:\n" + s
        throw_exception(s)

cdef api object c2numpy_int(int *A, int len):
    """
    Construct the integer NumPy array by copying the data.
    """
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="int32")
    cdef int *pvec = <int *>vec.data
    memcpy(pvec, A, len*sizeof(int))
    return vec

cdef api object c2numpy_int_inplace(int *A, int len):
    """
    Construct the integer NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_INT, A)

cdef api object c2numpy_double(double *A, int len):
    """
    Construct the double NumPy array by copying the data.
    """
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="double")
    cdef double *pvec = <double *>vec.data
    memcpy(pvec, A, len*sizeof(double))
    return vec

cdef api object c2numpy_double_inplace(double *A, int len):
    """
    Construct the double NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, A)

#cdef api object c2numpy_double_complex_inplace(cplx *A, int len):
#    """
#    Construct the double NumPy array inplace (don't copy any data).
#    """
#    cdef npy_intp dim = len
#    return PyArray_SimpleNewFromData(1, &dim, NPY_COMPLEX128, A)

_AA = None

cdef api void numpy2c_int_inplace(object A_n, int **A_c, int *n):
    """
    Returns the C array, that points to the numpy array (inplace).

    Only if strides != sizeof(int), the data get copied first.

    Important note: you need to use the A_c array immediately after calling
    this function in C, otherwise numpy could deallocate the array, especially
    if the _AA global variable was deallocated.
    """
    cdef ndarray A = A_n
    if not (A.ndim == 1 and A.strides[0] == sizeof(int)):
        from numpy import array
        A = array(A.flat, dtype="int32")
        # this is needed so that numpy doesn't dealocate the arrays
        global _AA
        _AA = A
    n[0] = len(A)
    A_c[0] = <int *>(A.data)

cdef api void numpy2c_double_inplace(object A_n, double **A_c, int *n):
    """
    Returns the C array, that points to the numpy array (inplace).

    Only if strides != sizeof(double), the data get copied first.

    Important note: you need to use the A_c array immediately after calling
    this function in C, otherwise numpy could deallocate the array, especially
    if the _AA global variable was deallocated.
    """
    cdef ndarray A = A_n
    if not (A.ndim == 1 and A.strides[0] == sizeof(double)):
        from numpy import array
        A = array(A.flat, dtype="double")
        # this is needed so that numpy doesn't dealocate the arrays
        global _AA
        _AA = A
    n[0] = len(A)
    A_c[0] = <double *>(A.data)

#cdef api void numpy2c_double_complex_inplace(object A_n, cplx **A_c, int *n):
#    """
#    Returns the C array, that points to the numpy array (inplace).
#
#    Only if strides != sizeof(double), the data get copied first.
#
#    Important note: you need to use the A_c array immediately after calling
#    this function in C, otherwise numpy could deallocate the array, especially
#    if the _AA global variable was deallocated.
#    """
#    cdef ndarray A = A_n
#    if not (A.nd == 1 and A.strides[0] == sizeof(cplx)):
#        from numpy import array
#        A = array(A.flat, dtype="complex128")
#        # this is needed so that numpy doesn't dealocate the arrays
#        global _AA
#        _AA = A
#    n[0] = len(A)
#    A_c[0] = <cplx *>(A.data)

cdef api void run_cmd(const_char_p text, object namespace):
    try:
        verbose = namespace.get("verbose")
        if verbose:
            print "got a text:", text
        if verbose:
            print "evaluting in the namespace:"
            print namespace
        code = compile(text, "", "exec")
        eval(code, {}, namespace)
        if verbose:
            print "new namespace:"
            print namespace
    except SystemExit, e:
        try:
            exit_code = int(e)
        except:
            exit_code = -1
        exit(exit_code)
    except:
        etype, value, tb = sys.exc_info()
        s = "".join(traceback.format_exception(etype, value, tb))
        s = "Exception raised in the Python code:\n" + s
        throw_exception(s)
