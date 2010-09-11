"""
Module for handling Fekete points approximations.
"""

from math import pi, sin, log, sqrt, exp, cos
import logging
import datetime
import time

from numpy import empty, arange, array, ndarray, zeros, real
from numpy.linalg import solve
from scipy.integrate import quadrature, fixed_quad
from sympy import vectorize

from hermes1d.h1d_wrapper import h1d_wrapper

from hydrogen import R_nl_numeric
import _fekete

_logger_Function = logging.getLogger("hermes1d.Function")

@vectorize(0, 1)
def feq(a, b, max_relative_error=1e-12, max_absolute_error=1e-12):
    a = float(a)
    b = float(b)
    # if the numbers are close enough (absolutely), then they are equal
    if abs(a-b) < max_absolute_error:
        return True
    # if not, they can still be equal if their relative error is small
    if abs(b) > abs(a):
        relative_error = abs((a-b)/b)
    else:
        relative_error = abs((a-b)/a)
    return relative_error <= max_relative_error

def generate_candidates(a, b, order):
    def cand(divisions, orders):
        if len(divisions) == 0:
            assert len(orders) == 1
            return Mesh1D((a, b), (order + orders[0],))
        elif len(divisions) == 1:
            assert len(orders) == 2
            return Mesh1D((a, _fekete.get_x_phys(divisions[0], a, b), b),
                    (order + orders[0], order + orders[1]))
        else:
            raise NotImplementedError()
    cands = [
            cand([], [1]),
            cand([], [2]),
            cand([0], [0, 0]),
            cand([0], [1, 0]),
            cand([0], [0, 1]),
            cand([0], [1, 1]),
            cand([0], [2, 0]),
            cand([0], [2, 1]),
            cand([0], [1, 2]),
            cand([0], [0, 2]),
            cand([0], [2, 2]),
            ]
    if order > 1:
        cands.extend([
            cand([0], [-1, 0]),
            cand([0], [0, -1]),
            cand([0], [-1, -1]),
            ])
    if order > 2:
        cands.extend([
            cand([0], [-2, 0]),
            cand([0], [-2, -1]),
            cand([0], [-1, -2]),
            cand([0], [0, -2]),
            cand([0], [-2, -2]),
            ])
    return cands

class Mesh1D(object):

    def __init__(self, points, orders):
        if not (len(points) == len(orders) + 1):
            raise Exception("points vs order mismatch")
        self._points = tuple(points)
        self._orders = tuple(orders)
        self.logger = logging.getLogger("hermes1d.Mesh1D")

    def iter_elems(self, elems_id=None):
        if elems_id is None:
            elems_id = range(len(self._orders))
        for i in elems_id:
            yield (self._points[i], self._points[i+1], self._orders[i])

    def __str__(self):
        return "<%s, orders: %s>" % (self._points, self._orders)

    def get_mesh_data(self):
        return self._points, self._orders

    def plot(self, call_show=True):
        try:
            from jsplot import plot, show
        except ImportError:
            from pylab import plot, show
        odd = False
        for a, b, order in self.iter_elems():
            fekete_points = _fekete.get_fekete_points_phys(order, a, b)
            if odd:
                format = "y-"
            else:
                format = "k-"
            odd = not odd
            plot([a, a, b, b], [0, order, order, 0], format, lw=2)
        if call_show:
            show()

    def element_at_point(self, x):
        """
        Returns the element id at a point "x".

        """
        for n, (a, b, order) in enumerate(self.iter_elems()):
            if b < x:
                continue
            return n
        raise Exception("Element at the point '%f' not found" % x)

    def get_element_by_id(self, id):
        return (self._points[id], self._points[id+1], self._orders[id])

    def get_node_id_by_coord(self, x):
        for i, node in enumerate(self._points):
            if feq(node, x):
                return i
        raise ValueError("Node not found")

    def restrict_to_elements(self, n1, n2):
        """
        Returns the submesh of elements [n1, n2] (in the slice notation).
        """
        return Mesh1D(self._points[n1:n2+1], self._orders[n1:n2])

    def restrict_to_interval(self, A, B):
        assert B > A
        n1 = self.element_at_point(A)
        n2 = self.element_at_point(B)
        points = []
        orders = []

        # first element:
        a, b, order = self.get_element_by_id(n1)
        if feq(b, A):
            pass
        else:
            if feq(a, A):
                pass
            elif a < A:
                a = A
            else:
                raise NotImplementedError()
            points.append(a)
            orders.append(order)

        #middle elements
        for n in range(n1+1, n2):
            a, b, order = self.get_element_by_id(n)
            points.append(a)
            orders.append(order)

        # last element:
        a, b, order = self.get_element_by_id(n2)
        if feq(a, A):
            pass
        elif a < A:
            a = A
        if feq(b, B):
            pass
        elif B < b:
            b = B
        else:
            raise NotImplementedError()
        if len(points) == 0 or not feq(points[-1], a):
            points.append(a)
            orders.append(order)
        points.append(b)

        return Mesh1D(points, orders)

    def union(self, o):
        p1 = self._points
        p2 = o._points
        p = list(p1)
        p.extend(p2)
        p.sort()
        points = [p[0]]
        for point in p[1:]:
            if feq(points[-1], point):
                continue
            points.append(point)
        # points now contains the sorted list of union points
        orders = []
        for n, p in enumerate(points[1:]):
            p1 = points[n]
            p2 = p
            mid = (p1+p2)/2.
            o1 = self._orders[self.element_at_point(mid)]
            o2 = o._orders[o.element_at_point(mid)]
            orders.append(max(o1, o2))

        return Mesh1D(points, orders)

    def __eq__(self, o):
        if isinstance(o, Mesh1D):
            if self._orders == o._orders:
                d = array(self._points) - array(o._points)
                if array(feq(d, 0.0)).all():
                    return True
        return False

    def __ne__(self, o):
        return not self.__eq__(o)

    def use_candidate(self, cand):
        n1 = self.get_node_id_by_coord(cand._points[0])
        n2 = self.get_node_id_by_coord(cand._points[-1])
        points = self._points[:n1] + cand._points + self._points[n2+1:]
        orders = self._orders[:n1] + cand._orders + self._orders[n2:]
        return Mesh1D(points, orders)

    def dofs(self):
        """
        Returns the number of DOFs needed to represent the function using H1
        elements.

        It assumes no Dirichlet BC at either end. If you have a Dirichlet BC at
        one end, substract 1, if you have Dirichlet BC at both ends, substract
        2.

        Example::

        >>> from fekete import Mesh1D, Function
        >>> from math import sin
        >>> m = Mesh1D((0, 10, 20, 30, 40), (6, 6, 6, 6))
        >>> f = Function(sin, m)
        >>> f.dofs()
        25

        If you prescribe values at both ends using Dirichlet, then f.dofs()
        returns 25 (as it doesn't know about any FEM), but the number of DOFs
        that you get in FEM is 23.

        """
        n = 1
        for a, b, order in self.iter_elems():
            n += order
        return n

    def adapt(self, exact_fns):
        """
        Adapts the mesh using the set of exact (or reference) solutions
        'exact_fns'.
        """
        g_mesh = self
        g_fns = [f.project_onto(g_mesh) for f in exact_fns]
        errors = [(g-f).l2_norm() for f, g in zip(exact_fns, g_fns)]
        graph = []
        for i in range(100000):
            self.logger.info("-"*80)
            self.logger.info("Adaptivity step: %s", i)
            self.logger.info("Current errors: %s", errors)
            self.logger.info("Current DOFs  : %s", g_mesh.dofs())
            graph.append((g_mesh.dofs(), errors))
            self.logger.info("Current mesh (and orders):")
            self.logger.info("   %s", g_mesh._points)
            self.logger.info("   %s", g_mesh._orders)
            self.logger.info("Determining which elements to refine...")
            #elems_to_refine = [0]
            elems_to_refine = None
            self.logger.info("Will refine the following elements: %s",
                    elems_to_refine)
            self.logger.info("    Done.")
            self.logger.info("Calculating errors for all candidates...")
            cands = []
            for n, f, g in zip(range(len(exact_fns)), exact_fns, g_fns):
                self.logger.info("  Considering mesh %d:" % n)
                c = g.get_candidates_with_errors(f, elems_to_refine)
                cands.extend(c)
            cands.sort(key=lambda x: x[1])
            #cands = cands[:4+len(cands)/3]
            cands = cands[:1]
            self.logger.info("    Done.")
            self.logger.info("Will use the following candidates:")
            for m, err in cands:
                self.logger.info("   %s %s %s", err, m._points, m._orders)
            self.logger.info("Refining mesh...")
            for m, _ in cands:
                g_mesh = g_mesh.use_candidate(m)
            self.logger.info("    Done.")
            self.logger.info("Projecting onto the new mesh (and calculating the error)...")
            g_fns = [f.project_onto(g_mesh) for f in exact_fns]
            errors = [(g-f).l2_norm() for f, g in zip(exact_fns, g_fns)]
            self.logger.info("    Done.")
            if max(errors) < 1e-9:
                break
        graph.append((g.dofs(), errors))
        return g_mesh, errors

    def refine_all_elements(self):
        pts = []
        orders = []
        for a, b, order in self.iter_elems():
            mid = (a+b)/2.
            pts.append(a)
            pts.append(mid)
            orders.append(order)
            orders.append(order)
            last = b
        pts.append(last)
        return Mesh1D(pts, orders)

    def increase_order(self):
        return Mesh1D(self._points, array(self._orders)+1)

class Function(h1d_wrapper.Function):
    """
    Represents a function on a mesh.

    The values are given in the Fekete points.
    """

    def __init__(self, obj, mesh=None):
        if not isinstance(mesh, Mesh1D):
            raise Exception("You need to specify a mesh.")
        self._mesh = mesh
        if isinstance(obj, (tuple, list, ndarray)):
            self._values = obj
        else:
            self._values = []
            for a, b, order in mesh.iter_elems():
                fekete_points = _fekete.get_fekete_points_phys(order, a, b)
                # Note: this is not a projection (it only evaluates obj in
                # fekete points), so the result is not the best
                # approximation possible:
                elem_values = [obj(p) for p in fekete_points]
                self._values.append(array(elem_values))

        self._poly_coeffs = {}
        self._fe_sol = None

        self.logger = _logger_Function

    def get_polynomial(self, x, values, a, b):
        """
        Returns the interpolating polynomial's coeffs.

        The len(values) specifies the order and we work in the element <a, b>
        """
        # Note: the version in _fekete is 2.6x faster
        n = len(values)
        A = empty((n, n), dtype="double")
        y = empty((n,), dtype="double")
        assert len(x) == n
        for i in range(n):
            for j in range(n):
                A[i, j] = _fekete.get_x_phys(x[i], a, b)**(n-j-1)
            y[i] = values[i]
        a = solve(A, y)
        return a

    def restrict_to_interval(self, A, B):
        """
        Returns the same function, with the mesh (domain) restricted to the
        interval (A, B).
        """
        m = self._mesh.restrict_to_interval(A, B)
        return Function(self, m)


    def eval_polynomial(self, coeffs, x):
        # This is about 15x faster
        return _fekete.eval_polynomial(coeffs, x)
        # than this:
        r = 0
        n = len(coeffs)
        for i, a in enumerate(coeffs):
            r += a*x**(n-i-1)
        return r

    def eval_polynomial_array(self, coeffs, x):
        # This is about 6x faster
        return _fekete.eval_polynomial_array(coeffs, x)
        # than this:
        r = zeros(len(x))
        n = len(coeffs)
        for i, a in enumerate(coeffs):
            r += a*x**(n-i-1)
        return r

    def __call__(self, x):
        for n, (a, b, order) in enumerate(self._mesh.iter_elems()):
            if b < x:
                continue
            y = _fekete.eval_poly(array([float(x)]), self._values[n],
                    a, b)[0]
            return y

    def eval_f(self, x):
        return self(x)

    def eval_dfdx(self, x):
        self.calculate_FE_coeffs()
        return self._fe_sol.deriv(x)

    def calculate_FE_coeffs(self):
        if self._fe_sol is None:
            from hermes1d.h1d_wrapper.h1d_wrapper import \
                    (assemble_projection_matrix_rhs, Mesh, FESolution)
            from hermes_common._hermes_common import CooMatrix
            pts, orders = self._mesh.get_mesh_data()
            m = Mesh(pts, orders)
            n_dof = m.assign_dofs()
            A = CooMatrix(n_dof)
            rhs = empty(n_dof)
            assemble_projection_matrix_rhs(m, A, rhs, self,
                    projection_type="L2")
            coeffs = solve(A.to_scipy_coo().todense(), rhs)
            self._fe_sol = FESolution(m, coeffs)

    def get_values_in_element(self, n, x):
        """
        Return the values in points 'x' in the element 'n'.

        'x' is a numpy array of points
        'n' is an element id

        It returns a numpy array of values of the function in points 'x'.
        All points 'x' must be in the element 'n'.
        """
        a, b, order = self._mesh.get_element_by_id(n)
        assert (a<=x).all()
        assert (x<=b).all()
        return _fekete.eval_poly(x, self._values[n], a, b)

    def get_polynomial_coeffs(self, n, values, a, b):
        if n not in self._poly_coeffs:
            vals = array(values)
            x = array(points[len(vals)-1])
            self._poly_coeffs[n] = _fekete.get_polynomial(x, vals, a, b)
        return self._poly_coeffs[n]

    def project_onto(self, mesh, proj_type="Fekete"):
        """
        Projects 'self' onto the 'mesh' using the 'proj_type' projection.

        proj_type == "Fekete"/"L2"/"H1"
        """
        if mesh == self._mesh:
            return self
        if proj_type == "Fekete":
            return Function(self, mesh)
        elif proj_type in ["L2", "H1"]:
            from hermes1d.h1d_wrapper.h1d_wrapper import \
                    (assemble_projection_matrix_rhs, Mesh, FESolution)
            from hermes_common._hermes_common import CooMatrix
            pts, orders = mesh.get_mesh_data()
            m = Mesh(pts, orders)
            n_dof = m.assign_dofs()
            A = CooMatrix(n_dof)
            rhs = empty(n_dof)
            assemble_projection_matrix_rhs(m, A, rhs, self,
                    projection_type=proj_type)
            coeffs = solve(A.to_scipy_coo().todense(), rhs)
            return FESolution(m, coeffs).to_discrete_function()
        else:
            raise ValueError("Unknown projection type")

    def project_onto_union(self, mesh):
        """
        The same as project_onto, only "mesh" is a subset of self._mesh.
        """
        if mesh == self._mesh:
            return self
        values = []
        n = 0
        for a, b, order in mesh.iter_elems():
            if a >= self._mesh._points[n+1]:
                n += 1
            fekete_points = _fekete.get_fekete_points_phys(order, a, b)
            elem_values = []
            # Note: this is not a projection (it only evaluates obj in
            # fekete points), so the result is not the best
            # approximation possible:
            elem_values = self.get_values_in_element(n, fekete_points)
            values.append(elem_values)
        return Function(values, mesh)

    def plot(self, call_show=True):
        try:
            from jsplot import plot, show
        except ImportError:
            from pylab import plot, show
        odd = False
        for n, (a, b, order) in enumerate(self._mesh.iter_elems()):
            fekete_points = _fekete.get_fekete_points_phys(order, a, b)
            vals = self._values[n]
            assert len(vals) == len(fekete_points)
            x = arange(a, b, 0.1)
            y = [self(_x) for _x in x]
            if odd:
                format = "g-"
            else:
                format = "r-"
            odd = not odd
            plot(x, y, format)
            plot(fekete_points, vals, "ko")
        if call_show:
            show()

    def __eq__(self, o):
        if isinstance(o, Function):
            for a, b, order in self._mesh.iter_elems():
                fekete_points = _fekete.get_fekete_points_phys(order, a, b)
                for p in fekete_points:
                    if not feq(self(p), o(p)):
                        return False
            for a, b, order in o._mesh.iter_elems():
                fekete_points = _fekete.get_fekete_points_phys(order, a, b)
                for p in fekete_points:
                    if not feq(self(p), o(p)):
                        return False
            return True
        else:
            return False

    def __ne__(self, o):
        return not self.__eq__(o)

    def __add__(self, o):
        if self._mesh == o._mesh:
            values = [array(x)+array(y) for x, y in zip(self._values,
                o._values)]
            return Function(values, self._mesh)
        else:
            union_mesh = self._mesh.union(o._mesh)
            return self.project_onto_union(union_mesh) + o.project_onto_union(union_mesh)

    def __sub__(self, o):
        return self + (-o)

    def __neg__(self):
        values = [-x for x in self._values]
        return Function(values, self._mesh)

    def __pow__(self, o):
        if isinstance(o, (int, long)):
            pts = self._mesh._points
            orders = empty(len(self._mesh._orders), dtype="int")
            values = []
            for n, (a, b, order) in enumerate(self._mesh.iter_elems()):
                order = o*order
                x = _fekete.get_fekete_points_phys(order, a, b)
                vals = _fekete.eval_poly(x, self._values[n], a, b)**o
                values.append(vals)
                orders[n] = order
            mesh = Mesh1D(pts, orders)
            return Function(values, mesh)
        else:
            return NotImplemented

    def get_mesh_adapt(self, max_order=12):
        return self._mesh

    def l2_norm(self, method="Fekete"):
        """
        Calculates the L2 norm of the function.

        method == "Fekete" or "FE": Use Legendre interpolation, or first
            project to a FE basis and then calculate the L2 norm
        """
        if method=="Fekete":
            r = 0
            for n, (a, b, order) in enumerate(self._mesh.iter_elems()):
                x, w = _fekete.get_gauss_points_phys(a, b, order+1)
                vals = _fekete.eval_poly(x, self._values[n], a, b)
                r += _fekete.int_f2(w, vals)
            return sqrt(r)
        elif method=="FE":
            self.calculate_FE_coeffs()
            return self._fe_sol.l2_norm()

    def h1_norm(self):
        """
        Calculates the H1 norm of the function.
        """
        self.calculate_FE_coeffs()
        return self._fe_sol.h1_norm()

    def get_candidates_with_errors(self, f, elems=None):
        """
        Returns a sorted list of all candidates and their errors.

        The best candidate is first, the worst candidate is last.

        The "f" is the reference function which we want to approximate using
        "self".
        """
        cand_with_errors = []
        orig = f.project_onto(self._mesh)
        dof_orig = orig.dofs()
        err_orig = (f - orig).l2_norm()
        for a, b, order in self._mesh.iter_elems(elems):
            self.logger.info("Considering element: %s %s %s", a, b, order)
            cands = generate_candidates(a, b, order)
            #self.logger.debug("Candidates %s", cands)
            #print "-"*40
            #print a, b, order
            best_crit = 1e10
            for m in cands:
                #self.logger.debug("Candidate: %s", m)
                #self.logger.debug("  create_mesh...")
                cand_mesh = self._mesh.use_candidate(m)
                #self.logger.debug("  project...")
                cand = f.project_onto(cand_mesh)
                dof_cand = cand.dofs()
                #self.logger.debug("  l2_norm...")
                # This is slow, because we are calculating the l2_norm over the
                # whole mesh and we are doing the union mesh as well. All we
                # need to do is integrate over the element we are interested
                # and cache the rest of the elements.
                diff = f - cand
                err_cand = diff.l2_norm()
                #self.logger.debug("  Choose...")
                if dof_cand == dof_orig:
                    if err_cand < err_orig:
                        # if this happens, it means that we can get better
                        # approximation with the same DOFs, so we definitely take
                        # this candidate:
                        print "DOF_cand == DOF_orig:", dof_cand, dof_orig, err_cand, err_orig
                        crit = -1e10
                    else:
                        crit = 1e10 # forget this candidate
                elif dof_cand > dof_orig:
                    # if DOF rises, we take the candidate that goes the steepest in
                    # the log/sqrt error/DOFs convergence graph
                    # we want 'crit' as negative as possible:
                    crit = (log(err_cand) - log(err_orig)) / \
                            sqrt(dof_cand - dof_orig)
                    #print crit, err_cand, err_orig, dof_cand, dof_orig
                else:
                    if err_cand < err_orig:
                        # if this happens, it means that we can get better
                        # approximation with less DOFs, so we definitely take
                        # this candidate:
                        print "Nice!", dof_cand, dof_orig, err_cand, err_orig
                        crit = -1e10
                    else:
                        crit = 1e10 # forget this candidate
                if crit < best_crit:
                    best_m = m
                    best_crit = crit
            cand_with_errors.append((best_m, best_crit))
                #if crit < -1e9:
                #    cand_with_errors.sort(key=lambda x: x[1])
                #    return cand_with_errors

        cand_with_errors.sort(key=lambda x: x[1])
        return cand_with_errors

    def dofs(self):
        return self._mesh.dofs()

def main():
    logger = logging.getLogger("hermes1d")

    pts = (0, 2, 4, 6, 8, 10, 50, 80, 100, 120, 150)
    orders = tuple([12]*(len(pts)-1))
    f_mesh = Mesh1D(pts, orders)
    logger.info("Projecting the hydrogen wavefunctions into a very fine mesh...")
    exact_fns = [
        Function(lambda x: R_nl_numeric(1, 0, x), f_mesh),
        Function(lambda x: R_nl_numeric(2, 0, x), f_mesh),
        Function(lambda x: R_nl_numeric(3, 0, x), f_mesh),
        Function(lambda x: R_nl_numeric(4, 0, x), f_mesh),
        ]
    logger.info("    Done.")
    g_mesh = Mesh1D((0, pts[-1]), (1,))
    g_mesh, errors = g_mesh.adapt(exact_fns)
    logger.info("Adaptivity converged.")
    logger.info("Final errors: %s", errors)
    logger.info("Final DOFs  : %s", g_mesh.dofs())
    logger.info("DOFs used to approximate the exact functions:    %s", \
        f_mesh.dofs())
    logger.info("Final mesh (and orders):")
    logger.info("   %s", g_mesh._points)
    logger.info("   %s", g_mesh._orders)
    #f.plot(False)
    #g.plot(False)
    #g._mesh.plot(False)
    #from pylab import figure, show, semilogy, grid
    #figure()
    #xx = [x[0] for x in graph]
    #yy = [x[1] for x in graph]
    #semilogy(xx, yy)
    #grid()
    #show()

if __name__ == "__main__":
    class TimeFormatter(logging.Formatter):
        """
        Uses datetime to format the time.

        Also counts the time interval from creating this instance.

        Otherwise it's the same as Formatter.  This is useful if you need to
        print milliseconds, that Formatter doesn't support.
        """
        def __init__(self, *args, **kwargs):
            self._created = time.time()
            logging.Formatter.__init__(self, *args, **kwargs)

        def formatTime(self, record, datefmt=None):
            delta = record.created - self._created
            return "%8.3f" % delta

    logging.basicConfig(level=logging.DEBUG, filename="fekete.log",
            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
            filemode="w")
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    formatter = TimeFormatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger("").addHandler(console)
    main()
