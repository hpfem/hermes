cdef extern from "../../hermes_common/matrix.h":

    cdef cppclass Matrix:
        pass

    cdef cppclass SparseMatrix(Matrix):
        pass

    cdef cppclass Vector:
        pass
