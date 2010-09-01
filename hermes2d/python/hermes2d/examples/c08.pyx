from hermes2d._hermes2d cimport scalar, WeakForm, H1Space, EdgePos, \
        FuncReal, GeomReal, ExtDataReal, BC_ESSENTIAL, c_BCType, \
        BC_NATURAL, int_v, c_Ord, create_Ord, FuncOrd, GeomOrd, ExtDataOrd, \
        int_dudy_dvdy, int_dudy_dvdx, int_dudx_dvdx, int_dudx_dvdy, H2D_SYM, int_v_ord

# Problem constants
cdef double E  = 200e9                               # Young modulus (steel)
cdef double nu = 0.3                                 # Poisson ration
cdef double l = (E * nu) / ((1 + nu) * (1 - 2*nu))   # external force in x-direction
cdef double mu = E / (2*(1 + nu))                    # external force in y-direction
cdef double f_1 = 1e4                                # first Lame constant
cdef double f_0 = 0                                  # second Lame constant

# Boundary marker (external force).
cdef int GAMMA_3_BDY = 3

# Boundary condition types
cdef c_BCType bc_type(int marker):
    if marker == 1:
        return <c_BCType>BC_ESSENTIAL
    return <c_BCType>BC_NATURAL

# Function values for essential(Dirichlet) boundary conditions    
cdef scalar essential_bc_values(int ess_bdy_marker, double x, double y):
    return 0.0

# Bilinear forms
cdef scalar bilinear_form_0_0(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return (l +2*mu) * int_dudx_dvdx(n, wt, u, v) + mu * int_dudy_dvdy(n, wt, u, v)

cdef scalar bilinear_form_0_1(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return l * int_dudy_dvdx(n, wt, u, v) + mu * int_dudx_dvdy(n, wt, u, v)

cdef scalar bilinear_form_1_1(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return mu * int_dudx_dvdx(n, wt, u, v) + (l + 2 * mu) * int_dudy_dvdy(n, wt, u, v)

# Linear forms
cdef scalar linear_form_surf_0(int n, double *wt, FuncReal **t, FuncReal *u, GeomReal *e, ExtDataReal *ext):
    return f_0 * int_v(n, wt, u)

cdef scalar linear_form_surf_1(int n, double *wt, FuncReal **t, FuncReal *u, GeomReal *e, ExtDataReal *ext):
    return f_1 * int_v(n, wt, u)

cdef c_Ord _order_bf(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    # XXX: with 9 it doesn't shout about the integration order, but gives wrong
    # results...
    return create_Ord(20)

cdef c_Ord _order_lf(int n, double *wt, FuncOrd **t, FuncOrd *u, GeomOrd *e, ExtDataOrd *ext):
    return int_v_ord(n, wt, u).mul_double(f_1)

def set_forms(WeakForm dp):
    dp.thisptr.add_matrix_form(0, 0, &bilinear_form_0_0, &_order_bf, H2D_SYM)
    dp.thisptr.add_matrix_form(0, 1, &bilinear_form_0_1, &_order_bf, H2D_SYM)
    dp.thisptr.add_matrix_form(1, 1, &bilinear_form_1_1, &_order_bf, H2D_SYM)
    dp.thisptr.add_vector_form_surf(0, &linear_form_surf_0, &_order_lf, GAMMA_3_BDY)
    dp.thisptr.add_vector_form_surf(1, &linear_form_surf_1, &_order_lf, GAMMA_3_BDY)

def set_bc(H1Space xdisp, H1Space ydisp):
    xdisp.thisptr.set_bc_types(&bc_type)
    ydisp.thisptr.set_bc_types(&bc_type)

    xdisp.thisptr.set_essential_bc_values(&essential_bc_values)
    ydisp.thisptr.set_essential_bc_values(&essential_bc_values)
