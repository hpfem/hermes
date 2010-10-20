from ref cimport PyObject

cdef extern from *:
    ctypedef char* char_p       "char*"
    ctypedef char* const_char_p "const char*"

    ctypedef struct RichPyTypeObject "PyTypeObject":
        object tp_new(RichPyTypeObject*, object, PyObject*)

cdef extern from "cpp_api.h":
    cdef struct c_CppCallback "CppCallback":
        void event(const_char_p msg)
    c_CppCallback *new_CppCallback "new CppCallback" ()
    void delete(...)
