"""
Module for handling Fekete points approximations.
"""

from math import pi, sin, log, sqrt, exp, cos

from numpy import empty, arange, array, ndarray
from numpy.linalg import solve

from scipy.integrate import quadrature

from gauss_lobatto_points import points

def get_x_phys(x_ref, a, b):
    return (a+b)/2. + x_ref*(b-a)/2.;

def generate_candidates(a, b, order):
    def cand(divisions, orders):
        if len(divisions) == 0:
            assert len(orders) == 1
            return Mesh1D((a, b), (order + orders[0],))
        elif len(divisions) == 1:
            assert len(orders) == 2
            return Mesh1D((a, get_x_phys(divisions[0], a, b), b),
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

    def iter_elems(self):
        for i in range(len(self._orders)):
            yield (self._points[i], self._points[i+1], self._orders[i])

    def plot(self, call_show=True):
        try:
            from jsplot import plot, show
        except ImportError:
            from pylab import plot, show
        odd = False
        for a, b, order in self.iter_elems():
            fekete_points = points[order]
            fekete_points = [get_x_phys(x, a, b) for x in fekete_points]
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

    def get_element_by_id(self, id):
        return list(self.iter_elems())[id]

    def get_node_id_by_coord(self, x):
        eps = 1e-10
        for i, node in enumerate(self._points):
            if abs(node-x) < eps:
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
        eps = 1e-12
        if abs(b-A) < eps:
            pass
        else:
            if abs(a-A) < eps:
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
        eps = 1e-12
        if abs(a-A) < eps:
            pass
        elif a < A:
            a = A
        if abs(b-B) < eps:
            pass
        elif B < b:
            b = B
        else:
            raise NotImplementedError()
        if len(points) == 0 or not (abs(points[-1] - a) < eps):
            points.append(a)
            orders.append(order)
        points.append(b)

        return Mesh1D(points, orders)

    def union(self, o):
        eps = 1e-12
        p1 = self._points
        p2 = o._points
        p = list(p1)
        p.extend(p2)
        p.sort()
        points = [p[0]]
        for point in p[1:]:
            if abs(points[-1] - point) < eps:
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
        eps = 1e-12
        if isinstance(o, Mesh1D):
            if self._orders == o._orders:
                d = array(self._points) - array(o._points)
                if (abs(d) < eps).all():
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

class Function(object):
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
                if order not in points:
                    raise ValueError("order '%d' not implememented" % order)
                fekete_points = points[order]
                elem_values = []
                # Note: this is not a projection (it only evaluates obj in
                # fekete points), so the result is not the best
                # approximation possible:
                for p in fekete_points:
                    p = get_x_phys(p, a, b)
                    val = obj(p)
                    elem_values.append(val)
                self._values.append(elem_values)

    def get_polynomial(self, values, a, b):
        """
        Returns the interpolating polynomial's coeffs.

        The len(values) specifies the order and we work in the element <a, b>
        """
        n = len(values)
        A = empty((n, n), dtype="double")
        y = empty((n,), dtype="double")
        x = points[n-1]
        assert len(x) == n
        for i in range(n):
            for j in range(n):
                A[i, j] = get_x_phys(x[i], a, b)**(n-j-1)
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
        r = 0
        n = len(coeffs)
        for i, a in enumerate(coeffs):
            r += a*x**(n-i-1)
        return r

    def __call__(self, x):
        for n, (a, b, order) in enumerate(self._mesh.iter_elems()):
            if b < x:
                continue
            # This can be made faster by using Lagrange interpolation
            # polynomials (no need to invert a matrix in order to get the
            # polynomial below). The results are however identical.
            coeffs = self.get_polynomial(self._values[n], a, b)
            return self.eval_polynomial(coeffs, x)

    def project_onto(self, mesh):
        # This is not a true projection, only some approximation:
        return Function(self, mesh)

    def plot(self, call_show=True):
        try:
            from jsplot import plot, show
        except ImportError:
            from pylab import plot, show
        odd = False
        for n, (a, b, order) in enumerate(self._mesh.iter_elems()):
            fekete_points = points[order]
            vals = self._values[n]
            assert len(vals) == len(fekete_points)
            fekete_points = [get_x_phys(x, a, b) for x in fekete_points]
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
        eps = 1e-12
        if isinstance(o, Function):
            for a, b, order in self._mesh.iter_elems():
                fekete_points = points[order]
                fekete_points = [get_x_phys(x, a, b) for x in fekete_points]
                for p in fekete_points:
                    if abs(self(p) - o(p)) > eps:
                        return False
            for a, b, order in o._mesh.iter_elems():
                fekete_points = points[order]
                fekete_points = [get_x_phys(x, a, b) for x in fekete_points]
                for p in fekete_points:
                    if abs(self(p) - o(p)) > eps:
                        return False
            return True
        else:
            return False

    def __ne__(self, o):
        return not self.__eq__(o)

    def __add__(self, o):
        if self._mesh == o._mesh:
            values = array(self._values) + array(o._values)
            return Function(values, self._mesh)
        else:
            union_mesh = self._mesh.union(o._mesh)
            return self.project_onto(union_mesh) + o.project_onto(union_mesh)

    def __sub__(self, o):
        return self + (-o)

    def __neg__(self):
        values = [-array(x) for x in self._values]
        return Function(values, self._mesh)


    def get_mesh_adapt(self, max_order=12):
        return self._mesh

    def l2_norm(self):
        """
        Returns the L2 norm of the function.
        """
        i = 0
        def f(x):
            return [self(_)**2 for _ in x]
        for a, b, order in self._mesh.iter_elems():
            val, err = quadrature(f, a, b)
            i += val
        return i

    def dofs(self):
        """
        Returns the number of DOFs needed to represent the function.
        """
        n = 1
        for a, b, order in self._mesh.iter_elems():
            n += order
        return n

    def get_candidates_with_errors(self, f):
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
        for a, b, order in self._mesh.iter_elems():
            cands = generate_candidates(a, b, order)
            #print "-"*40
            #print a, b, order
            for m in cands:
                cand_mesh = self._mesh.use_candidate(m)
                cand = f.project_onto(cand_mesh)
                dof_cand = cand.dofs()
                err_cand = (f - cand).l2_norm()
                if dof_cand == dof_orig:
                    if err_cand < err_orig:
                        # if this happens, it means that we can get better
                        # approximation with the same DOFs, so we definitely take
                        # this candidate:
                        print "XXX", dof_cand, dof_orig, err_cand, err_orig
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
                cand_with_errors.append((m, crit))
        cand_with_errors.sort(key=lambda x: x[1])
        return cand_with_errors

def main():
    f_mesh = Mesh1D((0, 2, 4, 6, 8, 10), (12, 12, 12, 12, 12))
    f = Function(lambda x: exp(-x), f_mesh)
    #mesh = f.get_mesh_adapt(max_order=1)
    g_mesh = Mesh1D((0, 10), (1,))
    #mesh.plot(False)
    g = f.project_onto(g_mesh)
    error = (g-f).l2_norm()
    graph = []
    while error > 1e-9:
        print error, g.dofs()
        graph.append((g.dofs(), error))
        print "  ", g._mesh._points
        print "  ", g._mesh._orders
        cand_with_errors = g.get_candidates_with_errors(f)
        #for m, err in cand_with_errors[:10]:
        #    print "   ", err, m._points, m._orders
        m, _ = cand_with_errors[0]
        g_mesh = g_mesh.use_candidate(m)
        g = f.project_onto(g_mesh)
        error = (g-f).l2_norm()
    graph.append((g.dofs(), error))
    print "Done.", error
    error = (g - f)
    print "error:     ", error.l2_norm()
    print "f dofs:    ", f.dofs()
    print "g dofs:    ", g.dofs()
    f.plot(False)
    g.plot(False)
    g._mesh.plot(False)
    from pylab import figure, show, semilogy, grid
    figure()
    xx = [x[0] for x in graph]
    yy = [x[1] for x in graph]
    semilogy(xx, yy)
    grid()
    show()

if __name__ == "__main__":
    main()
