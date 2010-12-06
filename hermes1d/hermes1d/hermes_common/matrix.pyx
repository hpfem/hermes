from numpy_utils cimport c2numpy_double_inplace, c2numpy_int_inplace

cdef class Matrix:

    def get_size(self):
        return self.thisptr.get_size()

cdef class SparseMatrix(Matrix):

    cdef matrix.SparseMatrix* as_SparseMatrix(self):
        return <matrix.SparseMatrix*> self.thisptr

cdef class CSCMatrix(SparseMatrix):
    """
    Represents a CSC matrix.
    """

    def __cinit__(self):
        self.thisptr = new matrix.UMFPackMatrix()

    def __dealloc__(self):
        del self.thisptr

    cdef matrix.UMFPackMatrix* as_UMFPackMatrix(self):
        return <matrix.UMFPackMatrix*> self.thisptr

    @property
    def IA(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef matrix.UMFPackMatrix *this = self.as_UMFPackMatrix()
        return c2numpy_int_inplace(this.get_Ai(), this.get_nnz())

    @property
    def JA(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef matrix.UMFPackMatrix *this = self.as_UMFPackMatrix()
        return c2numpy_int_inplace(this.get_Ap(), self.get_size()+1)

    @property
    def A(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef matrix.UMFPackMatrix *this = self.as_UMFPackMatrix()
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
        self.thisptr = new matrix.UMFPackVector()

    def __dealloc__(self):
        del self.thisptr
