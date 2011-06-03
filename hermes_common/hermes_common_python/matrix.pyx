from numpy_utils cimport c2numpy_double_inplace, c2numpy_int_inplace

cdef class Matrix:

    cpdef int get_size(self):
        return self.thisptr.get_size()

cdef class SparseMatrix(Matrix):

    cdef cSparseMatrix* as_SparseMatrix(self):
        return <cSparseMatrix*> self.thisptr

cdef class CSCMatrix(SparseMatrix):
    """
    Represents a CSC matrix.
    """

    def __cinit__(self):
        self.thisptr = new cCSCMatrix()

    def __dealloc__(self):
        del self.thisptr

    cdef cCSCMatrix* as_CSCMatrix(self):
        return <cCSCMatrix*> self.thisptr

    @property
    def IA(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef cCSCMatrix *this = self.as_CSCMatrix()
        return c2numpy_int_inplace(this.get_Ai(), this.get_nnz())

    @property
    def JA(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef cCSCMatrix *this = self.as_CSCMatrix()
        return c2numpy_int_inplace(this.get_Ap(), self.get_size()+1)

    @property
    def A(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef cCSCMatrix *this = self.as_CSCMatrix()
        return c2numpy_double_inplace(this.get_Ax(), this.get_nnz())

    def to_scipy_csc(self):
        """
        Converts itself to the scipy sparse CSC format.
        """
        from scipy.sparse import csc_matrix
        n = self.get_size()
        return csc_matrix((self.A, self.IA, self.JA), shape=(n, n))

    def __str__(self):
        return str(self.to_scipy_csc())

cdef class Vector:
    pass

cdef class AVector(Vector):
    """
    A Vector, represented by a C array internally.
    """

    def __cinit__(self):
        self.thisptr = new cUMFPackVector()

    def __dealloc__(self):
        del self.thisptr

    cdef cUMFPackVector* as_UMFPackVector(self):
        return <cUMFPackVector*> self.thisptr

    def to_numpy(self):
        cdef cUMFPackVector *this = self.as_UMFPackVector()
        return c2numpy_double_inplace(this.get_c_array(), this.length())
