# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

from libc.math cimport sin, cos, sqrt
import traceback

from numpy import empty, array
from numpy cimport ndarray

from hermes_common.matrix cimport SparseMatrix, Vector
cimport hermes1d
from hermes_common.cython_utils cimport PY_NEW
from hermes_common.numpy_utils cimport c2numpy_double
from hermes1d.fekete._fekete cimport get_gauss_points_phys, int_f2, \
        int_f2_f2

cdef class Element:
    cdef hermes1d.Element *thisptr

    @property
    def p(self):
        return self.thisptr.p

    @property
    def x1(self):
        return self.thisptr.x1

    @property
    def x2(self):
        return self.thisptr.x2

    def get_coeffs(self, int comp=0):
        cdef ndarray[double, mode="c"] coeffs = \
                empty((self.p+1,), dtype="double")
        self.thisptr.get_coeffs(0, comp, &coeffs[0])
        return coeffs

cdef api object c2py_Element(hermes1d.Element *h):
    cdef Element n
    n = <Element>PY_NEW(Element)
    n.thisptr = h
    return n

cdef class Mesh:

    def __init__(self, *args):
        cdef double a, b
        cdef int n_elem, p_init, eq_num, n
        cdef ndarray[double, mode="c"] pts
        cdef ndarray[int, mode="c"] p
        cdef ndarray[int, mode="c"] m
        cdef ndarray[int, mode="c"] div
        if len(args) == 5:
            a, b, n_elem, p_init, eq_num = args
            self.thisptr = new hermes1d.Space(a, b, n_elem, p_init, eq_num,
                    1, 0)
        elif len(args) == 4:
            a, b, n_elem, p_init = args
            eq_num = 1
            self.thisptr = new hermes1d.Space(a, b, n_elem, p_init, eq_num,
                    1, 0)
        elif len(args) == 2:
            pts = array(args[0], dtype="double")
            p = array(args[1], dtype="int32")
            if not (len(pts) == len(p) + 1):
                raise ValueError("len(pts) must be equal to len(p) + 1")
            n_elem = len(p)
            m = array([1]*len(p), dtype="int32")
            div = array([1]*len(p), dtype="int32")
            eq_num = 1
            self.thisptr = new hermes1d.Space(n_elem, &pts[0], &p[0],
                    &m[0], &div[0], eq_num, 1, 0)
        else:
            raise ValueError("Don't understand the arguments")
        self.delptr = True

    def __dealloc__(self):
        if self.delptr:
            del self.thisptr

    def set_coeff_vector(self, ndarray[double, mode="c"] sol not None,
            int comp):
        self.thisptr.set_coeff_vector(&sol[0], comp)

    def assign_dofs(self):
        return self.thisptr.assign_dofs()

    def plot_to_file(self, filename):
        self.thisptr.plot(filename)

    def get_mesh_data(self):
        """
        Returns the list of node coordinates and element orders.
        """
        pts = []
        p = []
        I = Iterator(self)
        cdef hermes1d.Element *e = I._next_active_element()
        while e != NULL:
            if len(pts) == 0:
                pts.append(e.x1)
            pts.append(e.x2)
            p.append(e.p)
            e = I._next_active_element()
        return array(pts), array(p)

    def set_bc_left_dirichlet(self, eq_n, val):
        self.thisptr.set_bc_left_dirichlet(eq_n, val)

    def set_bc_right_dirichlet(self, eq_n, val):
        self.thisptr.set_bc_right_dirichlet(eq_n, val)

    def replicate(self):
        cdef Mesh m = c2py_Mesh(self.thisptr.replicate())
        m.delptr = True
        return m

    def get_n_active_elem(self):
        return self.thisptr.get_n_active_elem()

    def get_num_dofs(self):
        return self.thisptr.get_num_dofs()

    def _reference_refinement(self, a, b):
        """
        Refines the mesh in place.
        """
        self.thisptr.reference_refinement(a, b)

    def reference_refinement(self, start_elem_id=None, num_to_ref=None):
        """
        Returns a new instance of refined mesh.
        """
        mesh_ref = self.replicate()
        if start_elem_id is None:
            start_elem_id = 0
        if num_to_ref is None:
            num_to_ref = mesh_ref.get_n_active_elem()
        mesh_ref._reference_refinement(start_elem_id, num_to_ref)
        return mesh_ref

    def copy(self):
        pts, orders = self.get_mesh_data()
        return Mesh(pts, orders)


cdef api object c2py_Mesh(hermes1d.Space *h):
    cdef Mesh n
    n = <Mesh>PY_NEW(Mesh)
    n.thisptr = h
    n.delptr = False
    return n

cdef class Linearizer:
    cdef hermes1d.Linearizer *thisptr
    cdef Mesh mesh

    def __cinit__(self, Mesh mesh):
        self.thisptr = new hermes1d.Linearizer(mesh.thisptr)
        self.mesh = mesh

    def plot_solution(self, out_filename, plotting_elem_subdivision):
        #cdef double *A
        cdef int n
        #numpy2c_double_inplace(y_prev, &A, &n)
        self.thisptr.plot_solution(out_filename, plotting_elem_subdivision)

    def get_xy(self, sol, int comp, int plotting_elem_subdivision):
        """
        Returns (x, y), where x, y are arrays of points.

        sol ... is the solution to linearize
        comp ... solution component
        """
        cdef double *x
        cdef double *y
        cdef int n
        self.mesh.set_coeff_vector(sol, comp)
        self.thisptr.get_xy_space(comp, plotting_elem_subdivision,
                &x, &y, &n)
        x_numpy = c2numpy_double(x, n)
        y_numpy = c2numpy_double(y, n)
        return x_numpy, y_numpy

cdef class Iterator:
    cdef hermes1d.Iterator *thisptr

    def __cinit__(self, Mesh mesh):
        self.thisptr = new hermes1d.Iterator(mesh.thisptr)

    def __dealloc__(self):
        del self.thisptr

    cdef hermes1d.Element* _first_active_element(self):
        return self.thisptr.first_active_element()

    def first_active_element(self):
        return c2py_Element(self._first_active_element())

    cdef hermes1d.Element* _next_active_element(self):
        return self.thisptr.next_active_element()

    def next_active_element(self):
        return c2py_Element(self._next_active_element())

    cdef hermes1d.Element* _last_active_element(self):
        return self.thisptr.last_active_element()

    def last_active_element(self):
        return c2py_Element(self._last_active_element())

class FESolution:
    """
    Represents an FE solution on an hp-mesh.

    It represents all the components (for multiple equations).
    """

    def __init__(self, mesh, coefs):
        self._mesh = mesh
        self._coefs = coefs

    def value(self, x, int comp=0):
        """
        Returns the value of the solution at a point 'x'.
        """
        self._mesh.set_coeff_vector(self._coefs, comp)

        pts = []
        p = []
        I = Iterator(self._mesh)
        cdef hermes1d.Element *e = I._next_active_element()
        while e != NULL:
            if e.x1 <= x and x <= e.x2:
                return e.get_solution_value(x, comp)
            e = I._next_active_element()

    def deriv(self, x, int comp=0):
        """
        Returns the derivative of the solution at a point 'x'.
        """
        self._mesh.set_coeff_vector(self._coefs, comp)

        pts = []
        p = []
        I = Iterator(self._mesh)
        cdef hermes1d.Element *e = I._next_active_element()
        while e != NULL:
            if e.x1 <= x and x <= e.x2:
                return e.get_solution_deriv(x, comp)
            e = I._next_active_element()

    def l2_norm(self, int comp=0):
        """
        Returns the L2 norm of the solution.
        """
        self._mesh.set_coeff_vector(self._coefs, comp)

        pts = []
        p = []
        I = Iterator(self._mesh)
        cdef hermes1d.Element *e = I._next_active_element()
        l2_norm_squared = 0
        while e != NULL:
            x, w = get_gauss_points_phys(e.x1, e.x2, e.p+1)
            vals = array([e.get_solution_value(_x, comp) for _x in x])
            l2_norm_squared += int_f2(w, vals)
            e = I._next_active_element()
        return sqrt(l2_norm_squared)

    def h1_norm(self, int comp=0):
        """
        Returns the H1 norm of the solution.
        """
        self._mesh.set_coeff_vector(self._coefs, comp)

        pts = []
        p = []
        I = Iterator(self._mesh)
        cdef hermes1d.Element *e = I._next_active_element()
        h1_norm_squared = 0
        while e != NULL:
            x, w = get_gauss_points_phys(e.x1, e.x2, e.p+1)
            vals = array([e.get_solution_value(_x, comp) for _x in x])
            derivs = array([e.get_solution_deriv(_x, comp) for _x in x])
            h1_norm_squared += int_f2_f2(w, vals, derivs)
            e = I._next_active_element()
        return sqrt(h1_norm_squared)

    def to_discrete_function(self, int comp=0):
        """
        Returns a Function instance, that uses Fekete points to represent the
        solution.
        """
        from hermes1d.fekete.fekete import Function, Mesh1D
        pts, orders = self._mesh.get_mesh_data()
        m = Mesh1D(pts, orders)
        # This can be sped up by evaluating at fekete points at each element at
        # once:
        self._component = comp
        return Function(self, m)

    def __call__(self, x):
        return self.value(x, self._component)

    def get_element_coeffs(self, int comp=0):
        """
        Returns a list of FE coefficients for each element, corresponding to
        the solution component 'comp'.
        """
        self._mesh.set_coeff_vector(self._coefs, comp)
        coeffs = []
        I = Iterator(self._mesh)
        cdef hermes1d.Element *e = I._next_active_element()
        while e != NULL:
            coeffs.append(c2py_Element(e).get_coeffs(comp))
            e = I._next_active_element()
        return coeffs

def calc_err_est(int norm, Mesh mesh, Mesh mesh_ref, int sln=0):
    cdef ndarray[double] err_array = empty(mesh.get_n_active_elem())
    err_total = hermes1d.calc_err_est(norm, mesh.thisptr,
            mesh_ref.thisptr, &(err_array[0]), sln)
    return err_total, err_array

def calc_solution_norm(int norm, Mesh mesh):
    return hermes1d.calc_solution_norm(norm, mesh.thisptr)

def adapt(norm, adapt_type, threshold,
        ndarray[double, mode="c"] err_squared_array, Mesh mesh, Mesh mesh_ref):
    hermes1d.adapt(norm, adapt_type, threshold, &err_squared_array[0],
            mesh.thisptr, mesh_ref.thisptr)


cdef void fn_sin(int n, double x[], double f[], double dfdx[]):
    for i in range(n):
        f[i] = sin(x[i])
        if dfdx != NULL:
            dfdx[i] = cos(x[i])

cdef class Function:

    cpdef double eval_f(self, double x) except *:
        return 0

    cpdef double eval_dfdx(self, double x) except *:
        return 0

cdef int fn(int n, double x[], double f[], double dfdx[], void *data):
    cdef Function func = <Function> data
    try:
        for i in range(n):
            f[i] = func.eval_f(x[i])
            if dfdx != NULL:
                dfdx[i] = func.eval_dfdx(x[i])
    except:
        traceback.print_exc()
        return 1
    return 0

def assemble_projection_matrix_rhs(Mesh mesh, SparseMatrix A, Vector rhs,
    f, projection_type=None):
    cdef int prj_type
    if projection_type == "L2":
        prj_type = hermes1d.H1D_L2_ortho_global
    elif projection_type == "H1":
        prj_type = hermes1d.H1D_H1_ortho_global
    else:
        raise ValueError("Unknown projection type")
    cdef Function func = f
    hermes1d.assemble_projection_matrix_rhs(mesh.thisptr, A.as_SparseMatrix(),
            rhs.thisptr, &fn, prj_type, <void *>func)
