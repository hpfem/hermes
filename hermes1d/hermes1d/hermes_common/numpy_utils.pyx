from libc.string cimport memcpy

from numpy import empty
from numpy cimport npy_intp, PyArray_SimpleNewFromData, NPY_INT, NPY_DOUBLE

cdef ndarray c2numpy_int(int *A, int len):
    """
    Construct the integer NumPy array by copying the data.
    """
    cdef ndarray vec = empty([len], dtype="int32")
    cdef int *pvec = <int *>vec.data
    memcpy(pvec, A, len*sizeof(int))
    return vec

cdef ndarray c2numpy_double(double *A, int len):
    """
    Construct the double NumPy array by copying the data.
    """
    cdef ndarray vec = empty([len], dtype="double")
    cdef double *pvec = <double *>vec.data
    memcpy(pvec, A, len*sizeof(double))
    return vec

cdef ndarray c2numpy_int_inplace(int *A, int len):
    """
    Construct the integer NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_INT, A)

cdef ndarray c2numpy_double_inplace(double *A, int len):
    """
    Construct the double NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, A)
