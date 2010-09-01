from hermes2d._hermes2d cimport scalar, H1Space, BC_ESSENTIAL, BC_NATURAL, c_BCType, int_u_v, int_grad_u_grad_v, int_v, int_grad_u_grad_v_ord, \
    FuncReal, GeomReal, ExtDataReal, WeakForm, c_Ord, create_Ord, FuncOrd, GeomOrd, ExtDataOrd, Solution, H2D_ANY, int_u_v_ord, int_v_ord

import math

# Problem parameters.
cdef double HEATCAP = 1e6
cdef double RHO = 3000
cdef double LAMBDA = 1e5
cdef double ALPHA = 10
cdef double T_INIT = 10
cdef double FINAL_TIME = 86400

# Global time variable
cdef double TIME = 0

cdef double TAU = 300.0

# Boundary markers
cdef int marker_ground = 1
cdef int marker_air = 2

# Boundary condition types
cdef c_BCType bc_types(int marker):
    if marker == marker_ground:
        return <c_BCType>BC_ESSENTIAL
    else:
        return <c_BCType>BC_NATURAL

# Function values for essential(Dirichlet) boundary markers
cdef scalar essential_bc_values(int ess_bdy_marker, double x, double y):
    return T_INIT

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_types)
    space.thisptr.set_essential_bc_values(&essential_bc_values)

# Time-dependent exterior temperature
def temp_ext(t):
    return T_INIT + 10. * math.sin(2*math.pi*t/FINAL_TIME);

cdef scalar bilinear_form(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return HEATCAP * RHO * int_u_v(n, wt, u, v) / TAU + LAMBDA * int_grad_u_grad_v(n, wt, u, v)

cdef c_Ord bilinear_form_ord(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    #return create_Ord(20)
    return int_grad_u_grad_v_ord(n, wt, u, v)

cdef scalar bilinear_form_surf(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return LAMBDA * ALPHA * int_u_v(n, wt, u, v)

cdef c_Ord bilinear_form_surf_ord(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    #return create_Ord(20)
    return int_u_v_ord(n, wt, u, v)

cdef scalar linear_form(int n, double *wt, FuncReal **t, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return HEATCAP * RHO * int_u_v(n, wt, ext.fn[0], v) / TAU;

cdef c_Ord linear_form_ord(int n, double *wt, FuncOrd **t, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    return create_Ord(20)
#return int_v_ord(n, wt, v)

cdef scalar linear_form_surf(int n, double *wt, FuncReal **t, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return LAMBDA * ALPHA * temp_ext(TIME) * int_v(n, wt, v)

def update_time(t):
    global TIME
    TIME = t

cdef c_Ord linear_form_surf_ord(int n, double *wt, FuncOrd **t, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    return create_Ord(20)
#return int_v_ord(n, wt, v)

def set_forms(WeakForm wf):
    wf.thisptr.add_matrix_form(0, 0, &bilinear_form, &bilinear_form_ord)
    wf.thisptr.add_matrix_form_surf(0, 0, &bilinear_form_surf, &bilinear_form_surf_ord, marker_air)

    #wf.thisptr.add_vector_form(0, &linear_form, &linear_form_ord, H2D_ANY, 1, s.thisptr)
    wf.thisptr.add_vector_form_surf(0, &linear_form_surf, &linear_form_surf_ord, marker_air)
