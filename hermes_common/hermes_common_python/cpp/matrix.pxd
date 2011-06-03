cdef extern from "matrix.h":

    cdef cppclass Matrix[double]:
        int get_size()

    cdef cppclass SparseMatrix[double](Matrix):
        pass

    cdef cppclass Vector[double]:
        int length()

cdef extern from "solver/umfpack_solver.h":

    cdef cppclass CSCMatrix[double](SparseMatrix):
        int *get_Ap()
        int *get_Ai()
        double *get_Ax()
        int get_nnz()

    cdef cppclass UMFPackVector[double](Vector):
        double *get_c_array()

    cdef cppclass UMFPackVector[double](CSCMatrix):
        pass
