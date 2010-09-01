cdef extern from "Python.h":
	ctypedef void PyObject
	void Py_INCREF(PyObject *x)
	void Py_DECREF(PyObject *x)


cdef extern from "mesh.h":
	cdef struct c_Mesh "Mesh":
		pass
	c_Mesh *new_Mesh "new Mesh" ()
	void del_Mesh "delete" (c_Mesh *o)
