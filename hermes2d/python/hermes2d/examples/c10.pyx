from hermes2d._hermes2d cimport scalar, H1Space, BC_ESSENTIAL, BC_NATURAL, c_BCType, int_u_v, int_grad_u_grad_v, int_v, \
    FuncReal, GeomReal, ExtDataReal, WeakForm, c_Ord, create_Ord, FuncOrd, GeomOrd, ExtDataOrd, Solution, H2D_ANY, H2D_SYM

# Problem constants
cdef int OMEGA_1 = 1
cdef int OMEGA_2 = 2
cdef int STATOR_BDY = 2
cdef double VOLTAGE = 50.0      # Voltage on the stator.
cdef double EPS1 = 1.0          # Relative electric permittivity in Omega_1.
cdef double EPS2 = 10.0         # Relative electric permittivity in Omega_2.

# Boundary condition types
cdef c_BCType bc_types(int marker):
    return <c_BCType>BC_ESSENTIAL

# Dirichlet boundary condition values
cdef scalar essential_bc_values(int ess_bdy_marker, double x, double y):
    if ess_bdy_marker == STATOR_BDY:
        return VOLTAGE
    else:
        return 0.0

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_types)
    space.thisptr.set_essential_bc_values(&essential_bc_values)


cdef scalar biform1(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return EPS1 * int_grad_u_grad_v(n, wt, u, v)

cdef scalar biform2(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return EPS2 * int_grad_u_grad_v(n, wt, u, v);

cdef c_Ord _order_bf(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    return create_Ord(20)


def set_forms(WeakForm wf):
    wf.thisptr.add_matrix_form(0, 0, &biform1, &_order_bf, H2D_SYM, OMEGA_1)
    wf.thisptr.add_matrix_form(0, 0, &biform2, &_order_bf, H2D_SYM, OMEGA_2)
