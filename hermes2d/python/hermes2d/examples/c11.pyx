from hermes2d._hermes2d cimport scalar, H1Space, BC_ESSENTIAL, BC_NATURAL, c_BCType, int_u_v, int_grad_u_grad_v, int_v, \
    FuncReal, GeomReal, ExtDataReal, WeakForm, c_Ord, create_Ord, FuncOrd, GeomOrd, ExtDataOrd, Solution, H2D_ANY, H2D_SYM, \
    int_dudx_dvdx, int_dudy_dvdy, int_dudy_dvdx, int_dudx_dvdy, Adapt, int_v_ord


# Problem constants
cdef double E  = 200e9  # Young modulus for steel: 200 GPa
cdef double nu = 0.3    # Poisson ratio
cdef double f  = 1e3    # load force: 10^3 N
cdef double lamda = (E * nu) / ((1 + nu) * (1 - 2*nu))
cdef double mu = E / (2*(1 + nu))

# Boundary markers
cdef int marker_left = 1
cdef int marker_top = 2

# Boundary condition types
cdef c_BCType bc_types(int marker):
    return <c_BCType>BC_NATURAL

# function values for essential(Dirichlet) boundary markers
# (if the return value is zero, this can be omitted)
cdef scalar essential_bc_values(int ess_bdy_marker, double x, double y):
    return 0
    

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_types)
    space.thisptr.set_essential_bc_values(&essential_bc_values)


cdef scalar bilinear_form_0_0(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return (lamda + 2*mu) * int_dudx_dvdx(n, wt, u, v) + mu * int_dudy_dvdy(n, wt, u, v)


cdef scalar bilinear_form_0_1(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return lamda * int_dudy_dvdx(n, wt, u, v) + mu * int_dudx_dvdy(n, wt, u, v)


cdef scalar bilinear_form_1_0(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return mu * int_dudy_dvdx(n, wt, u, v) + lamda * int_dudx_dvdy(n, wt, u, v)


cdef scalar bilinear_form_1_1(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return mu * int_dudx_dvdx(n, wt, u, v) + (lamda + 2*mu) * int_dudy_dvdy(n, wt, u, v)


cdef scalar linear_form_surf_1(int n, double *wt, FuncReal **t, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return -f * int_v(n, wt, v)

cdef c_Ord _order_bf(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    return create_Ord(20)

cdef c_Ord _order_lf(int n, double *wt, FuncOrd **t, FuncOrd *u, GeomOrd *e, ExtDataOrd *ext):
    #return int_v_ord(n, wt, u).mul_double(-f)
    return create_Ord(20)

def set_wf_forms(WeakForm wf):
    wf.thisptr.add_matrix_form(0, 0, &bilinear_form_0_0, &_order_bf, H2D_SYM)
    wf.thisptr.add_matrix_form(0, 1, &bilinear_form_0_1, &_order_bf, H2D_SYM)
    wf.thisptr.add_matrix_form(1, 1, &bilinear_form_1_1, &_order_bf, H2D_SYM)
    wf.thisptr.add_vector_form_surf(1, &linear_form_surf_1, &_order_lf, marker_top)


def set_hp_forms(Adapt hp):
    hp.thisptr.set_error_form(0, 0, &bilinear_form_0_0, &_order_bf)
    hp.thisptr.set_error_form(0, 1, &bilinear_form_0_1, &_order_bf)
    hp.thisptr.set_error_form(1, 0, &bilinear_form_1_0, &_order_bf)
    hp.thisptr.set_error_form(1, 1, &bilinear_form_1_1, &_order_bf)
