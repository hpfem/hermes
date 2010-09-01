from hermes2d._hermes2d cimport scalar, WeakForm, H1Space, EdgePos, \
        FuncReal, GeomReal, ExtDataReal, BC_ESSENTIAL, c_BCType, \
        BC_NATURAL, c_Ord, create_Ord, FuncOrd, GeomOrd, \
        ExtDataOrd, ExtDataReal, FuncReal, GeomReal, H2D_SYM, int_F_v

import math

# Problem parameters

def a_11(x, y):
    if (y > 0):
        return 1 + x*x + y*y
    else:
        return 1

def a_22(x, y):
    if (y > 0):
        return 1
    else:
        return 1 + x*x + y*y

def a_12(x, y):
    return 1

def a_21(x, y):
  return 1

def a_1(x, y):
    return 0.0

def a_2(x, y):
    return 0.0

def a_0(x, y):
  return 0.0

cdef double rhs(double x, double y):
    return 1 + x*x + y*y

def g_D(x, y):
    return -math.cos(math.pi*x)

cdef double g_N(double x, double y):
    return 0

# Boundary condition types
cdef c_BCType bc_types(int marker):
    if marker == 1:
        return <c_BCType>BC_ESSENTIAL
    else:
        return <c_BCType>BC_NATURAL

# Dirichlet boundary condition values
cdef scalar essential_bc_values(int marker, double x, double y):
    return g_D(x, y)

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_types)
    space.thisptr.set_essential_bc_values(&essential_bc_values)

def set_forms(WeakForm wf):
    wf.thisptr.add_matrix_form(0, 0, &bilinear_form, &bilinear_form_ord, H2D_SYM);
    wf.thisptr.add_vector_form(0, &linear_form, &linear_form_ord);
    wf.thisptr.add_vector_form_surf(0, &linear_form_surf, &linear_form_surf_ord, 2);

# (Volumetric) bilinear form
cdef scalar bilinear_form(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    result = 0
    for i in range(n):
        x = e.x[i]
        y = e.y[i]
        result += (a_11(x, y)*u.dx[i]*v.dx[i] + a_12(x, y)*u.dy[i]*v.dx[i] + a_21(x, y)*u.dx[i]*v.dy[i] + a_22(x, y)*u.dy[i]*v.dy[i] + a_1(x, y)*u.dx[i]*v.val[i] + a_2(x, y)*u.dy[i]*v.val[i] + a_0(x, y)*u.val[i]*v.val[i]) * wt[i]

    return result

# Integration order for the bilinear form
cdef c_Ord bilinear_form_ord(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    #return u.val[0] * v.val[0] * e.x[0] * e.x[0]
    return create_Ord(20)

# Surface linear form (natural boundary conditions)
cdef scalar linear_form_surf(int n, double *wt, FuncReal **t, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return int_F_v(n, wt, g_N, v, e)

# Integration order for surface linear form
cdef c_Ord linear_form_surf_ord(int n, double *wt, FuncOrd **t, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    #return v.val[0] * e.x[0] * e.x[0]
    return create_Ord(20)

# Volumetric linear form (right-hand side)
cdef scalar linear_form(int n, double *wt, FuncReal **t, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return int_F_v(n, wt, rhs, v, e)

# Integration order for the volumetric linear form
cdef c_Ord linear_form_ord(int n, double *wt, FuncOrd **t, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    #return v.val[0] * e.x[0] * e.x[0]
    return create_Ord(20)
