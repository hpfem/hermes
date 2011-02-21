cdef extern from "math.h":

    double c_sqr "sqr"(double x)
    double c_sqrt "sqrt"(double x)
    double c_atan "atan"(double x)
    double c_pi "M_PI"

cdef extern from "stdlib.h":

    ctypedef unsigned long size_t
    void *malloc (size_t size)
    void free(void *mem)
    void *memcpy(void *dst, void *src, long n)

    void exit(int exit_code)

cdef extern from "arrayobject.h":

    ctypedef int intp

    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef int flags

cdef extern from "Python.h":
    ctypedef void PyObject
    void Py_INCREF(PyObject *x)
    void Py_DECREF(PyObject *x)

cdef extern from "stdcython.h":
    void init_global_empty_tuple()
    object PY_NEW(object t)

cdef extern from "hermes2d.h":

    # This is just the C++ "delete" statement
    void delete(...)

    void hermes2d_initialize(int* argc, char* argv[])
    void hermes2d_finalize()
    void throw_exception(char *msg)
    #void finish_glut_main_loop(int force_quit)
    int BC_ESSENTIAL "BC_ESSENTIAL"
    int BC_NATURAL "BC_NATURAL"
    int BC_NONE "BC_NONE"
    int H2D_ANTISYM "H2D_ANTISYM"
    int H2D_UNSYM "H2D_UNSYM"
    int H2D_SYM "H2D_SYM"
    int H2D_ANY "H2D_ANY"
    int c_H2D_FN_VAL "H2D_FN_VAL"
    int c_H2D_FN_VAL_0 "H2D_FN_VAL_0"
    int c_H2D_FN_DX "H2D_FN_DX"
    int c_H2D_FN_DY "H2D_FN_DY"
    int c_H2D_FN_DX_0 "H2D_FN_DX_0"
    int c_H2D_FN_DY_0 "H2D_FN_DY_0"    
    int c_H2D_FN_DXX "H2D_FN_DXX"
    int c_H2D_FN_DYY "H2D_FN_DYY"
    int c_H2D_FN_DXY "H2D_FN_DXY"
    int c_H2D_FN_DEFAULT "H2D_FN_DEFAULT"
    int c_H2D_FN_ALL "H2D_FN_ALL"
    int c_H2D_EPS_LOW "H2D_EPS_LOW"
    int c_H2D_EPS_NORMAL "H2D_EPS_NORMAL"
    int c_H2D_EPS_HIGH "H2D_EPS_HIGH"
    int c_verbose_mode "__h2d_report_verbose"
    int c_info_mode "__h2d_report_info"
    int c_warn_integration "__h2d_report_warn_intr"

    int c_H2D_P_ISO "RefinementSelectors::H2D_P_ISO"
    int c_H2D_P_ANISO "RefinementSelectors::H2D_P_ANISO"
    int c_H2D_H_ISO "RefinementSelectors::H2D_H_ISO"
    int c_H2D_H_ANISO "RefinementSelectors::H2D_H_ANISO"
    int c_H2D_HP_ISO "RefinementSelectors::H2D_HP_ISO"
    int c_H2D_HP_ANISO_H "RefinementSelectors::H2D_HP_ANISO_H"
    int c_H2D_HP_ANISO_P "RefinementSelectors::H2D_HP_ANISO_P"
    int c_H2D_HP_ANISO "RefinementSelectors::H2D_HP_ANISO"

    int c_H2D_PREFER_SYMMETRIC_MESH "RefinementSelectors::H2D_PREFER_SYMMETRIC_MESH"
    int c_H2D_APPLY_CONV_EXP_DOF "RefinementSelectors::H2D_APPLY_CONV_EXP_DOF"

    int c_H2D_TOTAL_ERROR_REL "H2D_TOTAL_ERROR_REL"
    int c_H2D_TOTAL_ERROR_ABS "H2D_TOTAL_ERROR_ABS"
    int c_H2D_ELEMENT_ERROR_REL "H2D_ELEMENT_ERROR_REL"
    int c_H2D_ELEMENT_ERROR_ABS "H2D_ELEMENT_ERROR_ABS"

    ctypedef double double4[4]
    ctypedef double double3[3]
    ctypedef double double2[2]
    ctypedef int int3[3]
    ctypedef int int2[2]

    cdef struct c_Nurbs "Nurbs":
        int degree
        int np
        double3 *pt
        int nk
        double *kv
        int ref

    cdef struct c_CurvMap "CurvMap":
        int toplevel
        int order
        c_Nurbs* nurbs[4]

    cdef struct c_Node "Node":
        int id
        unsigned ref
        unsigned type
        unsigned bnd
        unsigned marker
        unsigned used
        double x, y

    cdef struct c_Element "Element":
        int id
        unsigned nvert
        unsigned active
        unsigned used
        int marker
        int userdata
        int iro_cache
        c_Node* vn[4]
        c_Node* en[4]
        c_Element* sons[4]
        c_CurvMap* cm
        double get_area()
        double get_diameter()

    cdef struct c_Mesh "Mesh":
        void load(char* filename)
        void load_str(char* mesh)
        void save(char* filename)
        void copy(c_Mesh *m)
        void refine_element_id(int id, int refinement)
        void refine_all_elements()
        void refine_towards_boundary(int marker, int depth)
        void refine_towards_vertex(int marker, int depth)
        c_Element* get_element(int id)
        int get_num_elements()
        int get_num_base_elements()
        int get_num_active_elements()
        int get_max_element_id()
        void convert_triangles_to_quads()
        void convert_quads_to_triangles()
    c_Mesh *new_Mesh "new Mesh" ()

    cdef struct c_H2DReader "H2DReader":
        void load(char* filename, c_Mesh *m)
        void load_str(char* mesh, c_Mesh *m)

    ctypedef struct c_H1Shapeset "H1Shapeset"
    c_H1Shapeset *new_H1Shapeset "new H1Shapeset" ()

    ctypedef struct c_L2Shapeset "L2Shapeset"
    c_L2Shapeset *new_L2Shapeset "new L2Shapeset" ()

    cdef struct c_PrecalcShapeset "PrecalcShapeset"
    c_PrecalcShapeset *new_PrecalcShapeset "new PrecalcShapeset" (c_H1Shapeset *s)

    ctypedef double scalar
    cdef struct EdgePos "EdgePos":
        int marker

    ctypedef enum c_BCType "BCType":
        pass

    cdef struct c_H1Space "H1Space":
        void set_uniform_order(int tri_order)
        int assign_dofs(int first_dof, int stride)
        void copy_orders(c_H1Space *s, int inc)
        int get_element_order(int id)
        int get_num_dofs()
        void set_bc_types(c_BCType (*bc_type_callback)(int marker))
        void set_essential_bc_values(scalar (*bc_value_callback_by_coord)(int ess_bdy_marker,
            double x, double y))
        void set_essential_bc_values_edge "set_essential_bc_values"(scalar (*bc_value_callback_by_edge)(EdgePos *ep))
    c_H1Space *new_H1Space "new H1Space" (c_Mesh *m, c_BCType (*bc_type_callback)(int marker), scalar (*bc_value_callback_by_coord)(int, double, double), int p_init, c_H1Shapeset *h)

    cdef struct c_L2Space "L2Space":
        void set_uniform_order(int tri_order)
        int assign_dofs(int first_dof, int stride)
        void copy_orders(c_L2Space *s, int inc)
        int get_element_order(int id)
        int get_num_dofs()
        void set_bc_types(c_BCType (*bc_type_callback)(int marker))
        void set_essential_bc_values(scalar (*bc_value_callback_by_coord)(int ess_bdy_marker,
            double x, double y))
        void set_essential_bc_values_edge "set_essential_bc_values"(scalar (*bc_value_callback_by_edge)(EdgePos *ep))
    c_L2Space *new_L2Space "new L2Space" (c_Mesh *m)

    ctypedef struct RealFunction "Function<double>":
        c_Element* get_active_element()
    cdef struct RefMap "RefMap"
    ctypedef struct c_ScalarFunction "Function<scalar>"

    ctypedef struct c_Ord "Ord":
        int get_order()
        c_Ord mul_double "operator*"(double right)
    c_Ord create_Ord "Ord"(int n)

    ctypedef struct FuncReal "Func<double>":
        double *val
        double *dx, *dy 
    ctypedef struct GeomReal "Geom<double>":
        int marker
        double *x, *y
    ctypedef struct ExtDataReal "ExtData<double>":
        FuncReal **fn
    ctypedef struct FuncOrd "Func<Ord>":
        c_Ord *val
        c_Ord *dx, *dy  
    ctypedef struct GeomOrd "Geom<Ord>":
        int marker
        c_Ord *x, *y    
    ctypedef struct ExtDataOrd "ExtData<Ord>":
        pass

    ctypedef scalar (*BiFormFnVol)(RealFunction *fu, RealFunction *fv,
            RefMap *ru, RefMap *rv)

    ctypedef scalar (*LiFormFnVol)(RealFunction *fv, RefMap *rv)

    cdef BiFormFnVol BF_SYM = <BiFormFnVol> 1
    cdef BiFormFnVol BF_ANTISYM = <BiFormFnVol> 2

    cdef struct c_Function "MeshFunction":
        pass

    cdef struct c_MeshFunction "MeshFunction":
        RefMap* get_refmap()
        c_Mesh* get_mesh()
        #scalar get_pt_value(double x, double y, int item)
        scalar get_pt_value(double x, double y)

    cdef struct c_Solution "Solution":
        void set_zero(c_Mesh *m)
        void set_const(c_Mesh *m, scalar c)
        void copy(c_Solution *s)
        void set_fe_solution(c_H1Space *s, c_PrecalcShapeset *pss, scalar *vec)
        void get_fe_solution(int *Ylen, scalar **Y)
    c_Solution *new_Solution "new Solution" ()

    cdef struct c_VonMisesFilter "VonMisesFilter"
    c_VonMisesFilter *new_VonMisesFilter "new VonMisesFilter" (c_MeshFunction *sln1, c_MeshFunction *sln2, double l, double m)

    cdef struct c_MagFilter "MagFilter"
    c_MagFilter *new_MagFilter "new MagFilter" (c_MeshFunction *sln1,
            c_MeshFunction *sln2, int item1, int item2)

    cdef struct c_DiffFilter "DiffFilter"
    c_DiffFilter *new_DiffFilter "new DiffFilter" (c_MeshFunction *sln1,
            c_MeshFunction *sln2, int item1, int item2)

    cdef struct c_SumFilter "SumFilter"
    c_SumFilter *new_SumFilter "new SumFilter" (c_MeshFunction *sln1,
            c_MeshFunction *sln2, int item1, int item2)

    cdef struct c_SquareFilter "SquareFilter"
    c_SquareFilter *new_SquareFilter "new SquareFilter" (c_MeshFunction *sln1,
            int item1)

    cdef struct c_WeakForm "WeakForm":
        void add_matrix_form(int i, int j, ...)
        void add_matrix_form_surf(int i, int j, ...)
        void add_vector_form(int i, ...)
        void add_vector_form_data(int i, void *data)
        void add_vector_form_surf(int i, ...)
    c_WeakForm *new_WeakForm "new WeakForm" (int neq)

    cdef struct c_CommonSolver "CommonSolver":
        pass

    cdef struct c_LinSystem "LinSystem":
        void set_spaces(int n, ...)
        void set_spaces2(int n, ...)
        void set_pss(int n, ...)
        void set_pss2(int n, ...)
        c_H1Space *get_space(int n)
        c_PrecalcShapeset *get_pss(int n)
        void copy(c_LinSystem *sys)
        void assemble()
        int get_num_dofs()
        int solve(int n, ...)
        int solve2(int n, ...)        
        void save_matrix_matlab(char *filename, char *varname)
        void get_matrix(int *Ap, int *Ai, scalar *Ax, int size)
        void get_rhs(scalar *RHS, int size)
    #c_LinSystem *new_LinSystem "new LinSystem" (c_WeakForm *wf,
    #        c_CommonSolver *solver)
    c_LinSystem *new_LinSystem "new LinSystem" (c_WeakForm *wf)#,
          #  c_H1Space *sp)
            
    cdef struct c_RefSystem "RefSystem":
        void assemble()
        c_H1Space *get_ref_space(int eq)
    c_RefSystem *new_RefSystem "new RefSystem" (c_LinSystem *ls)

    #cdef struct c_DiscreteProblem "DiscreteProblem":
    #    void set_num_equations(int neq)
    #    void set_external_fns(int n, ...)
    #    void set_bilinear_form(int i, int j, ...)
    #    void set_linear_form(int i, ...)
    #    void create_matrix()
    #    void set_quiet(int quiet)
    #    void assemble_matrix_and_rhs()
    #    void solve_system(int n, ...)
    #    void save_matrix_coo(char *filename)
    #    void free_matrix()
    #c_DiscreteProblem *new_DiscreteProblem "new DiscreteProblem" ()

    ctypedef enum c_CandList "RefinementSelectors::CandList":
        pass

    ctypedef enum c_SelOption "RefinementSelectors::SelOption":
        pass

    ctypedef struct c_H1SpaceTuple "Hermes::Tuple<Space*>":
        void (* push_back)(c_H1Space*)

    ctypedef struct c_L2SpaceTuple "Hermes::Tuple<Space*>":
        void (* push_back)(c_L2Space*)

    ctypedef struct c_SolutionTuple "Hermes::Tuple<Solution*>":
        void (* push_back)(c_Solution*)

    ctypedef struct c_ProjBasedSelector "RefinementSelectors::ProjBasedSelector":
        void set_error_weights(double, double, double)
        void set_option(c_SelOption, int)

    ctypedef struct c_H1ProjBasedSelector "RefinementSelectors::H1ProjBasedSelector":
        pass

    c_H1ProjBasedSelector *new_H1ProjBasedSelector "new RefinementSelectors::H1ProjBasedSelector" (c_CandList, double, int, c_H1Shapeset*)

    ctypedef struct c_L2ProjBasedSelector "RefinementSelectors::L2ProjBasedSelector":
        pass
        
    c_L2ProjBasedSelector *new_L2ProjBasedSelector "new RefinementSelectors::L2ProjBasedSelector" (c_CandList, double, int, c_L2Shapeset*)

    ctypedef struct c_Adapt "Adapt":
        void set_solutions(c_SolutionTuple, c_SolutionTuple)
        double calc_error(int)
        void set_error_form(int i, int j, ...)
        int adapt(c_ProjBasedSelector*, double, int, int, double)

    ctypedef struct c_H1Adapt "H1Adapt":
        pass

     #c_H1Adapt *new_H1Adapt "new H1Adapt" (c_H1SpaceTuple)
    c_H1Adapt *new_H1Adapt "new H1Adapt" (c_LinSystem*)

    ctypedef struct c_L2Adapt "L2Adapt":
        pass

    #c_L2Adapt *new_L2Adapt "new L2Adapt" (c_L2SpaceTuple)
    c_L2Adapt *new_L2Adapt "new L2Adapt" (c_LinSystem*)

    cdef struct c_Linearizer "Linearizer":
        void process_solution(c_MeshFunction* sln, ...)
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
    c_Linearizer *new_Linearizer "new Linearizer" ()

    cdef struct c_Vectorizer "Vectorizer":
        void process_solution(c_MeshFunction* xsln, ...)
        double4* get_vertices()
        int2* get_dashes()
        int get_num_dashes()
    c_Vectorizer *new_Vectorizer "new Vectorizer" ()


    double int_u "int_u<double, double>"(int n, double *wt, FuncReal *u)
    double int_l2_norm(RealFunction* fu, RefMap* ru)
    double l2_norm(c_MeshFunction* fu)
    double h1_norm(c_MeshFunction* fu)
    double integrate(c_MeshFunction *sln)
    double int_grad_u_grad_v "int_grad_u_grad_v<double, double>"(int n,
            double *wt, FuncReal *u, FuncReal *v)
    c_Ord int_grad_u_grad_v_ord "int_grad_u_grad_v<Ord, Ord>"(int n,
            double *wt, FuncOrd *u, FuncOrd *v)
    double int_x_grad_u_grad_v "int_x_grad_u_grad_v<double, double>"(int n,
            double *wt, FuncReal *u, FuncReal *v)
    double int_u_dvdx "int_u_dvdx<double, double>"(...)
    double int_u_dvdy "int_u_dvdy<double, double>"(...)
    double int_u_v "int_u_v<double, double>"(int n, double *wt, FuncReal *u,
            FuncReal *v)
    c_Ord int_u_v_ord "int_u_v<Ord, Ord>"(int n, double *wt, FuncOrd *u,
            FuncOrd *v)
    double int_F_u_v(...)
    double int_F_v "int_F_v<double, double>"(...)
    c_Ord int_F_v_ord "int_F_v<Ord, Ord>"(...)
    double int_v "int_v<double, double>"(int n, double *wt, FuncReal *v)
    c_Ord int_v_ord "int_v<Ord, Ord>"(int n, double *wt, FuncOrd *v)
    double int_w_nabla_u_v "int_w_nabla_u_v<double, double>"(...)
    #double surf_int_G_v(RealFunction *fv, RefMap *rv, EdgePos *ep)
    #double surf_int_v(RealFunction *fv, RefMap *rv, EdgePos *ep)
    #double surf_int_u_v(RealFunction *fu, RealFunction *fv,
    #        RefMap *ru, RefMap *rv, EdgePos *ep)
    #double int_a_dudx_dvdx_b_dudy_dvdy(double a, RealFunction *fu, double b,
    #        RealFunction *fv, RefMap *ru, RefMap *rv)
    #double int_a_dudx_dvdy_b_dudy_dvdx(double a, RealFunction *fu, double b,
    #        RealFunction *fv, RefMap *ru, RefMap *rv)
    double int_dudx_dvdx "int_dudx_dvdx<double, double>"(...)
    double int_dudx_dvdy "int_dudx_dvdy<double, double>"(...)
    double int_dudy_dvdy "int_dudy_dvdy<double, double>"(...)
    double int_dudy_dvdx "int_dudy_dvdx<double, double>"(...)


    cdef struct c_View "View":
        void set_title(char *title)
        void set_min_max_range(double min, double max)

    void View_wait "View::wait"()

    cdef struct c_ScalarView "ScalarView":
        void show(c_MeshFunction *s, ...)
        void show_mesh(int show)
        void show_scale(int show)
        void set_title(char *title)
        void set_min_max_range(double min, double max)
        void wait()
    c_ScalarView *new_ScalarView "new ScalarView" (char *title, ...)

    cdef struct c_BaseView "BaseView":
        void show(c_H1Space *s)
    c_BaseView *new_BaseView "new BaseView" (char *title, ...)

    cdef struct c_MeshView "MeshView":
        void show(c_Mesh *s)
    c_MeshView *new_MeshView "new MeshView" (char *title, ...)

    cdef struct c_OrderView "OrderView":
        void show(c_H1Space *s)
    c_OrderView *new_OrderView "new OrderView" (char *title, ...)

    cdef struct c_VectorView "VectorView":
        void show(c_MeshFunction *s1, c_MeshFunction *s2, double eps)
        void show_scale(int show)
        void set_min_max_range(double min, double max)
    c_VectorView *new_VectorView "new VectorView" (char *title, ...)

    #cdef struct c_MatrixView "MatrixView":
    #    pass
        #void show(c_DiscreteProblem *ep)
    #c_MatrixView *new_MatrixView "new MatrixView" (char *title, ...)

    double2 *transform(c_Element *e)
    void element_polygonal_boundary(c_Element *e, double2 **tp, int *n)

cdef extern from "dummy_solver.h":

    cdef struct c_DummySolver "DummySolver":
        pass
    c_DummySolver *new_DummySolver "new DummySolver" ()


cdef class LinSystem:
    cdef c_LinSystem *thisptr
    cdef object _spaces
    cdef object _pss

cdef class RefSystem(LinSystem):
    pass

cdef class CommonSolver:
    cdef c_CommonSolver *thisptr

cdef class DummySolver(CommonSolver):
    pass

cdef class WeakForm:
    cdef c_WeakForm *thisptr

cdef class H1Space:
    cdef c_H1Space *thisptr

cdef class L2Space:
    cdef c_L2Space *thisptr

cdef class Transformable:
    pass

cdef class Function(Transformable):
    cdef c_Function *thisptr

cdef class ScalarFunction(Function):
    pass

cdef class MeshFunction(ScalarFunction):
    pass

cdef class Solution(MeshFunction):
    pass

cdef class Linearizer:
    cdef c_Linearizer *thisptr

cdef class Vectorizer(Linearizer):
    pass

cdef class ProjBasedSelector:
    cdef c_ProjBasedSelector *thisptr

cdef class H1ProjBasedSelector(ProjBasedSelector):
    pass

cdef class L2ProjBasedSelector(ProjBasedSelector):
    pass

cdef class Adapt:
    cdef c_Adapt *thisptr

cdef class H1Adapt(Adapt):
    pass

cdef class L2Adapt(Adapt):
    pass

cdef class Mesh:
    cdef c_Mesh *thisptr
