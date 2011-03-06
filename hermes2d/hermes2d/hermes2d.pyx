from numpy cimport ndarray

from hermes_common.numpy_utils cimport c2numpy_int, c2numpy_double
cimport hermes2d_defs

cdef class MeshFunction:
    pass

cdef class Solution(MeshFunction):

    def __cinit__(self):
        self.thisptr = new hermes2d_defs.Solution()

    def __dealloc__(self):
        del self.thisptr

    cdef hermes2d_defs.Solution *getptr(self):
        return <hermes2d_defs.Solution *>(self.thisptr)

cdef class Linearizer:
    """
    Linearizes the solution.

    It returns the triangles and vertices and you can then use it to visualize
    the solution.

    Example::

        In [40]: l = Linearizer()

        In [44]: l.process_solution(sln)
        Linearizer: 5519 verts, 10847 tris in 0.05 sec

        In [45]: l.get_vertices()
        Out[45]:
        array([[  0.00000000e+00,  -1.00000000e+00,  -2.22396971e-17],
               [  1.00000000e+00,  -1.00000000e+00,  -1.64798730e-17],
               [ -1.00000000e+00,   0.00000000e+00,   8.09899023e-17],
               ...,
               [  1.48437500e-01,  -1.56250000e-02,   1.62359362e-01],
               [  1.32812500e-01,   0.00000000e+00,   1.56012622e-01],
               [  1.32812500e-01,  -1.56250000e-02,   1.50562411e-01]])

        In [46]: l.get_triangles()
        Out[46]:
        array([[   3, 5448,   29],
               [  27, 5445,   28],
               [  29,   28,   26],
               ...,
               [5499, 5498, 5479],
               [5510, 5493, 5491],
               [5513, 5508, 5491]], dtype=int32)

    """
    cdef hermes2d_defs.Linearizer *thisptr

    def __cinit__(self):
        self.thisptr = new hermes2d_defs.Linearizer()

    def __dealloc__(self):
        del self.thisptr

    def process_solution(self, MeshFunction sln):
        self.thisptr.process_solution(sln.thisptr)

    def get_vertices(self):
        """
        Returns the list of vertices.

        It's a list of triples, where each triple contains (x, y, val), where
        x, y are the "x" and "y" coordinates of the vertex and "val" is the
        value of the solution at the vertex.

        Example::

            In [45]: l.get_vertices()
            Out[45]:
            array([[  0.00000000e+00,  -1.00000000e+00,  -2.22396971e-17],
                   [  1.00000000e+00,  -1.00000000e+00,  -1.64798730e-17],
                   [ -1.00000000e+00,   0.00000000e+00,   8.09899023e-17],
                   ...,
                   [  1.48437500e-01,  -1.56250000e-02,   1.62359362e-01],
                   [  1.32812500e-01,   0.00000000e+00,   1.56012622e-01],
                   [  1.32812500e-01,  -1.56250000e-02,   1.50562411e-01]])

        """
        cdef hermes2d_defs.double3 *vert = self.thisptr.get_vertices()
        cdef int nvert = self.thisptr.get_num_vertices()
        cdef ndarray vec = c2numpy_double(<double *>vert, 3*nvert)
        return vec.reshape((nvert, 3))

    def get_num_vertices(self):
        return self.thisptr.get_num_vertices()

    def get_triangles(self):
        """
        Returns a list of triangles.

        The list contains triples of vertices IDs. Use get_vertices() to obtain
        vertices coordinates.

        Example::

            In [46]: l.get_triangles()
            Out[46]:
            array([[   3, 5448,   29],
                   [  27, 5445,   28],
                   [  29,   28,   26],
                   ...,
                   [5499, 5498, 5479],
                   [5510, 5493, 5491],
                   [5513, 5508, 5491]], dtype=int32)

        """
        cdef hermes2d_defs.int3 *tri = self.thisptr.get_triangles()
        cdef int ntri = self.thisptr.get_num_triangles()
        cdef ndarray vec = c2numpy_int(<int *>tri, 3*ntri)
        return vec.reshape((ntri, 3))

    def get_num_triangles(self):
        return self.thisptr.get_num_triangles()

    def get_edges(self):
        """
        Returns a list of edges.

        The list contains triples of vertices IDs. Use get_vertices() to obtain
        vertices coordinates.

        Example::

            In [47]: l.get_edges()
            Out[47]:
            array([[   3,   27,    0],
                   [  27,   24,    0],
                   [  24,   30,    0],
                   ...,
                   [5339, 5070,    4],
                   [5070, 5077,    4],
                   [5077,   11,    4]], dtype=int32)

        """
        cdef hermes2d_defs.int3 *edges = self.thisptr.get_edges()
        cdef int nedges = self.thisptr.get_num_edges()
        cdef ndarray vec = c2numpy_int(<int *>edges, 3*nedges)
        return vec.reshape((nedges, 3))

    def get_num_edges(self):
        return self.thisptr.get_num_edges()

    def get_min_value(self):
        return self.thisptr.get_min_value()

    def get_max_value(self):
        return self.thisptr.get_max_value()

    def save_data(self, char* filename):
        self.thisptr.save_data(filename)

    def load_data(self, char* filename):
        self.thisptr.load_data(filename)
