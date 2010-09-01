cdef class CppCallback:
    cdef c_CppCallback *thisptr

    def __init__(self):
        self.thisptr = new_CppCallback()

    def __del__(self):
        delete(self.thisptr)

    def event(self, msg):
        self.thisptr.event(msg)

cdef inline PY_NEW(T):
    # The line below is roughly equivalent to "return T()", except that the
    # "__init__()" function is not being called:
    return (<RichPyTypeObject*>T).tp_new((<RichPyTypeObject*>T), (), NULL)

cdef api object c2py_CppCallback(c_CppCallback *m):
    cdef CppCallback c
    c = <CppCallback>PY_NEW(CppCallback)
    c.thisptr = <c_CppCallback *>m
    return c

