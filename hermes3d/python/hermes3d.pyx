cdef class Mesh:
	cdef c_Mesh *thisptr

	def __cinit__(self):
		self.thisptr = new_Mesh()
