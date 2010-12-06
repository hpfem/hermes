cdef extern from "../../hermes_common/matrix.h":

    cdef cppclass Matrix:
        int get_size()

    cdef cppclass SparseMatrix(Matrix):
        pass

    cdef cppclass Vector:
        pass

cdef extern from "../../hermes_common/solver/umfpack_solver.h":

    cdef cppclass UMFPackMatrix(SparseMatrix):
        int *get_Ap()
        int *get_Ai()
        double *get_Ax()
        int get_nnz()

    cdef cppclass UMFPackVector(Vector):
        pass

