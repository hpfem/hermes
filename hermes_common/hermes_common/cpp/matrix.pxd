cdef extern from "matrix.h":

    cdef cppclass Matrix:
        int get_size()

    cdef cppclass SparseMatrix(Matrix):
        pass

    cdef cppclass Vector:
        int length()

cdef extern from "solver/umfpack_solver.h":

    cdef cppclass CSCMatrix(SparseMatrix):
        int *get_Ap()
        int *get_Ai()
        double *get_Ax()
        int get_nnz()

    cdef cppclass UMFPackVector(Vector):
        double *get_c_array()

    cdef cppclass UMFPackVector(CSCMatrix):
        pass
