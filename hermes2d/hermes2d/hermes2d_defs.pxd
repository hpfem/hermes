cdef extern from "hermes2d.h":

    ctypedef double scalar
    ctypedef double double4[4]
    ctypedef double double3[3]
    ctypedef int int3[3]
    ctypedef int int2[2]

    cdef cppclass Function:
        pass

    cdef cppclass MeshFunction(Function):
        #RefMap* get_refmap()
        #c_Mesh* get_mesh()
        scalar get_pt_value(double x, double y)

    cdef cppclass Solution(MeshFunction):
        #void set_zero(c_Mesh *m)
        #void set_const(c_Mesh *m, scalar c)
        void copy(Solution *s)
        #void set_fe_solution(c_H1Space *s, c_PrecalcShapeset *pss, scalar *vec)
        #void get_fe_solution(int *Ylen, scalar **Y)

    cdef cppclass Linearizer:
        void process_solution(MeshFunction* sln, ...)
        void lock_data()
        void unlock_data()
        double3* get_vertices()
        int get_num_vertices()
        int3* get_triangles()
        int get_num_triangles()
        int3* get_edges()
        int get_num_edges()
        double get_min_value()
        double get_max_value()
        void save_data(char* filename)
        void load_data(char* filename)
