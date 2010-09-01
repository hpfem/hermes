from hermes2d._hermes2d cimport scalar, WeakForm, H1Space, EdgePos, \
        FuncReal, GeomReal, ExtDataReal, BC_ESSENTIAL, c_BCType, \
        BC_NATURAL, int_v, int_grad_u_grad_v, c_Ord, create_Ord, FuncOrd, GeomOrd, ExtDataOrd, \
        int_v_ord
#from hermes2d._hermes2d.forms import bilinear_form, bilinear_form_ord
        
CONST_F = -1.0                       # right-hand side
CONST_GAMMA = [-0.5, 1.0, -0.5]

# Boundary condition types
# Note: natural means Neumann, Newton, or any other type of condition
# where the solution value is not prescribed.
cdef c_BCType bc_type_05(int marker):
    if marker == 4:
        return <c_BCType>BC_ESSENTIAL
    else:
        return <c_BCType>BC_NATURAL
        
# Function values for essential(Dirichlet) boundary markers
cdef scalar essential_bc_values(int marker, double x, double y):
    return 0.0

cdef scalar bilinear_form(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal
        *e, ExtDataReal *ext):
    return int_grad_u_grad_v(n, wt, u, v)
    
cdef c_Ord bilinear_form_ord(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v, GeomOrd
        *e, ExtDataOrd *ext):
    # surprisingly, uncommenting this makes it a lot slower:
    #print int_grad_u_grad_v_ord(n, wt, u, v).get_order(),
    #return int_grad_u_grad_v_ord(n, wt, u, v)
    return create_Ord(20)

cdef scalar linear_form_05(int n, double *wt, FuncReal **t, FuncReal *u, GeomReal *e, ExtDataReal *ext):
    return CONST_F*int_v(n, wt, u)

cdef c_Ord _order_lf(int n, double *wt, FuncOrd **t, FuncOrd *u, GeomOrd
        *e, ExtDataOrd *ext):
    return int_v_ord(n, wt, u).mul_double(CONST_GAMMA[e.marker-1])

cdef scalar linear_form_surf_05(int n, double *wt, FuncReal **t, FuncReal *u, GeomReal
        *e, ExtDataReal *ext):
    return CONST_GAMMA[e.marker-1] * int_v(n, wt, u)

def set_forms(WeakForm wf):
    wf.thisptr.add_matrix_form(0, 0, &bilinear_form, &bilinear_form_ord);
    wf.thisptr.add_vector_form(0, &linear_form_05, &_order_lf); 
    wf.thisptr.add_vector_form_surf(0, &linear_form_surf_05, &_order_lf);

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_type_05)
    space.thisptr.set_essential_bc_values(&essential_bc_values)
