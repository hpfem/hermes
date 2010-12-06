from cpython.object cimport PyObject

cdef extern from "Python.h":
    ctypedef struct RichPyTypeObject "PyTypeObject":
        object tp_new(RichPyTypeObject*, object, PyObject*)

cdef inline PY_NEW(T):
    # The line below is roughly equivalent to "return T()", except that the
    # "__init__()" function is not being called:
    return (<RichPyTypeObject*>T).tp_new((<RichPyTypeObject*>T), (), NULL)
