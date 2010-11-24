from libc.string cimport memcpy

from numpy import empty

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
