from hermes2d._hermes2d cimport scalar, FuncReal, GeomReal, WeakForm, \
        int_grad_u_grad_v, int_grad_u_grad_v_ord, int_v, malloc, ExtDataReal, \
        c_Ord, create_Ord, FuncOrd, GeomOrd, ExtDataOrd, int_v_ord

cdef scalar bilinear_form(int n, double *wt, FuncReal **t, FuncReal *u, FuncReal *v, GeomReal
        *e, ExtDataReal *ext):
    return int_grad_u_grad_v(n, wt, u, v)

cdef c_Ord bilinear_form_ord(int n, double *wt, FuncOrd **t, FuncOrd *u, FuncOrd *v, GeomOrd
        *e, ExtDataOrd *ext):
    # surprisingly, uncommenting this makes it a lot slower:
    #print int_grad_u_grad_v_ord(n, wt, u, v).get_order(),
    #return int_grad_u_grad_v_ord(n, wt, u, v)
    return create_Ord(20)

cdef scalar linear_form_p2(int n, double *wt, FuncReal **t, FuncReal *u, GeomReal
        *e, ExtDataReal *ext):
    return 2*int_v(n, wt, u)

cdef scalar linear_form_m1(int n, double *wt, FuncReal **t, FuncReal *u, GeomReal
        *e, ExtDataReal *ext):
    return -int_v(n, wt, u)

cdef scalar linear_form_m4(int n, double *wt, FuncReal **t, FuncReal *u, GeomReal
        *e, ExtDataReal *ext):
    return -4*int_v(n, wt, u)

cdef c_Ord _order_lf(int n, double *wt, FuncOrd **t, FuncOrd *u, GeomOrd
        *e, ExtDataOrd *ext):
    return int_v_ord(n, wt, u)

def set_forms(WeakForm wf, int k=2):
    wf.thisptr.add_matrix_form(0, 0, &bilinear_form, &bilinear_form_ord)
    
    if k == 2:
        wf.thisptr.add_vector_form(0, &linear_form_p2, &_order_lf)
    elif k == -1:
        wf.thisptr.add_vector_form(0, &linear_form_m1, &_order_lf)
    elif k == -4:
        wf.thisptr.add_vector_form(0, &linear_form_m4, &_order_lf)
    else:
        raise NotImplementedError()
