from libc.string cimport memcpy

from numpy import empty
from numpy cimport (npy_intp, PyArray_SimpleNewFromData, NPY_INT, NPY_DOUBLE,
        import_array)
# Initialize NumPy
import_array()

cdef ndarray c2numpy_int(int *A, int n):
    """
    Construct the integer NumPy array by copying the data.
    """
    cdef ndarray[int, mode="c"] a = empty(n, dtype="int32")
    memcpy(&a[0], A, n*sizeof(int))
    return a

cdef ndarray c2numpy_double(double *A, int n):
    """
    Construct the double NumPy array by copying the data.
    """
    cdef ndarray[double, mode="c"] a = empty([n], dtype="double")
    memcpy(&a[0], A, n*sizeof(double))
    return a

cdef ndarray c2numpy_int_inplace(int *A, int n):
    """
    Construct the integer NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = n
    return PyArray_SimpleNewFromData(1, &dim, NPY_INT, A)

cdef ndarray c2numpy_double_inplace(double *A, int n):
    """
    Construct the double NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = n
    return PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, A)
