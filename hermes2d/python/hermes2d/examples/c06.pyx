from hermes2d._hermes2d cimport scalar, FuncReal, GeomReal, WeakForm, \
        int_grad_u_grad_v, int_grad_u_grad_v_ord, int_v, int_v_ord, malloc, ExtDataReal, c_Ord, create_Ord, \
        FuncOrd, GeomOrd, ExtDataOrd, int_u_v, int_u_v_ord, BC_NATURAL, BC_ESSENTIAL, c_BCType, \
        H1Space
# Problem parameters
cdef double T1 = 30.0                # prescribed temperature on Gamma_3
cdef double T0 = 20.0                # outer temperature on Gamma_1
cdef double H  = 0.05                # Heat flux on Gamma_1

# Boundary condition types
cdef c_BCType bc_type_06(int marker):
    if marker == 3:
        return <c_BCType>BC_ESSENTIAL
    else:
        return <c_BCType>BC_NATURAL

# Function values for essential(Dirichlet) boundary markers
cdef scalar essential_bc_values(int ess_bdy_marker, double x, double y):
    if ess_bdy_marker == 3:
        return 100.
    return T1

cdef scalar bilinear_form(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal
        *e, ExtDataReal *ext):
    return int_grad_u_grad_v(n, wt, u, v)

cdef scalar bilinear_form_surf_06(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v,
        GeomReal *e, ExtDataReal *ext):
    if e.marker != 1:
        return 0.
    return H * int_u_v(n, wt, u, v)

cdef scalar linear_form_surf_06(int n, double *wt, FuncReal **t, FuncReal *v, GeomReal *e,
        ExtDataReal *ext):
    if e.marker != 1:
        return 0.
    return T0 * H * int_v(n, wt, v)

cdef c_Ord _order_bf(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v, GeomOrd
        *e, ExtDataOrd *ext):
    return int_grad_u_grad_v_ord(n, wt, u, v)

cdef c_Ord _order_bf_surf_06(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v,
        GeomOrd *e, ExtDataOrd *ext):
    return int_u_v_ord(n, wt, u, v)

cdef c_Ord _order_lf(int n, double *wt, FuncOrd **t, FuncOrd *u, GeomOrd
        *e, ExtDataOrd *ext):
    return int_v_ord(n, wt, u)

def set_forms(WeakForm dp):
    dp.thisptr.add_matrix_form(0, 0, &bilinear_form, &_order_bf)
    dp.thisptr.add_matrix_form_surf(0, 0, &bilinear_form_surf_06, &_order_bf_surf_06)
    dp.thisptr.add_vector_form_surf(0, &linear_form_surf_06, &_order_lf)

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_type_06)
    space.thisptr.set_essential_bc_values(&essential_bc_values)
