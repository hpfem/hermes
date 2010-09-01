from hermes2d._hermes2d cimport scalar, FuncReal, GeomReal, WeakForm, \
        int_grad_u_grad_v, int_v, malloc, ExtDataReal, c_Ord, create_Ord, \
        FuncOrd, GeomOrd, ExtDataOrd, H1Space, BC_ESSENTIAL, BC_NATURAL, c_BCType, \
        c_atan, c_pi, c_sqrt, c_sqr, int_F_v, int_grad_u_grad_v_ord
from hermes2d._hermes2d cimport int_F_v_ord

cdef double fn(double x, double y):
    return c_atan(60 * (c_sqrt(c_sqr(x-1.25) + c_sqr(y+0.25)) - c_pi/3))

cdef double rhs(double x, double y):
    cdef double t1 = c_sqrt(16*(x*x + y*y) - 40*x + 8*y + 26)
    cdef double t2 = 3600*(x*x + y*y) - 9000*x + 1800*y
    return -(240 * (t2 + 5849 - 400*c_pi*c_pi)) / (t1 * c_sqr(5851 + t2 - \
              600*t1*c_pi + 400*c_pi*c_pi))

cdef scalar bilinear_form(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal
        *e, ExtDataReal *ext):
    return int_grad_u_grad_v(n, wt, u, v)

cdef scalar linear_form(int n, double *wt, FuncReal **t, FuncReal *u, GeomReal
        *e, ExtDataReal *ext):
    return -int_F_v(n, wt, rhs, u, e)

cdef c_BCType bc_type(int marker):
    return <c_BCType>BC_ESSENTIAL

cdef scalar essential_bc_values(int ess_bdy_marker, double x, double y):
    return fn(x, y)

cdef c_Ord _order_bf(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v, GeomOrd
        *e, ExtDataOrd *ext):
    return int_grad_u_grad_v_ord(n, wt, u, v)

cdef c_Ord _order_lf(int n, double *wt, FuncOrd **t, FuncOrd *u, GeomOrd
            *e, ExtDataOrd *ext):
    # this doesn't work, unless we redefine rhs using Ord:
    #return int_F_v_ord(n, wt, rhs, u, e).mul_double(-1)
    return create_Ord(20)

def set_forms(WeakForm dp):
    dp.thisptr.add_matrix_form(0, 0, &bilinear_form, &_order_bf)
    dp.thisptr.add_vector_form(0, &linear_form, &_order_lf)

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_type)
    space.thisptr.set_essential_bc_values(&essential_bc_values)
