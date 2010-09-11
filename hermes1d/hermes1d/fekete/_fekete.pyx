cimport cython

from numpy cimport ndarray
from numpy import zeros, empty, real, array, pi, arange, cos
from numpy.linalg import solve
from scipy.special.orthogonal import p_roots
from gauss_lobatto_points import points

def eval_polynomial_orig(coeffs, x):
    r = 0
    n = len(coeffs)
    for i, a in enumerate(coeffs):
        r += a*x**(n-i-1)
    return r

@cython.boundscheck(False)
def eval_polynomial(ndarray[double, mode="c"] coeffs not None, double x):
    cdef double r=0
    cdef unsigned n = len(coeffs)
    cdef unsigned i
    for i in range(n):
        r += coeffs[i]*x**(n-i-1)
    return r

def get_x_phys_orig(x_ref, a, b):
    return (a+b)/2. + x_ref*(b-a)/2.

cpdef double get_x_phys(double x_ref, double a, double b):
    return (a+b)/2. + x_ref*(b-a)/2.

def eval_polynomial_array_orig(coeffs, x):
    r = zeros(len(x))
    n = len(coeffs)
    for i, a in enumerate(coeffs):
        r += a*x**(n-i-1)
    return r

@cython.boundscheck(False)
cpdef eval_polynomial_array(ndarray[double, mode="c"] coeffs, ndarray[double, mode="c"] x):
    cdef unsigned n_coeffs = len(coeffs)
    cdef unsigned n_x = len(x)
    cdef ndarray[double, mode="c"] r=zeros(n_x)
    cdef unsigned i, j
    for j in range(n_x):
        for i in range(n_coeffs):
            r[j] += coeffs[i]*x[j]**(n_coeffs-i-1)
    return r

def get_polynomial_orig(x, values, a, b):
    """
    Returns the interpolating polynomial's coeffs.

    The len(values) specifies the order and we work in the element <a, b>
    """
    n = len(values)
    A = empty((n, n), dtype="double")
    y = empty((n,), dtype="double")
    assert len(x) == n
    for i in range(n):
        for j in range(n):
            A[i, j] = get_x_phys(x[i], a, b)**(n-j-1)
        y[i] = values[i]
    a = solve(A, y)
    return a

@cython.boundscheck(False)
cpdef ndarray[double, mode="c"] get_polynomial(ndarray[double, mode="c"] x,
        ndarray[double, mode="c"] values,
        double a, double b):
    """
    Returns the interpolating polynomial's coeffs.

    The len(values) specifies the order and we work in the element <a, b>
    """
    cdef unsigned n = len(values)
    cdef ndarray[double, ndim=2, mode="c"] A = empty((n, n), dtype="double")
    cdef ndarray[double, mode="c"] y = empty((n,), dtype="double")
    assert len(x) == n
    cdef unsigned i, j
    for i in range(n):
        for j in range(n):
            A[i, j] = get_x_phys(x[i], a, b)**(n-j-1)
        y[i] = values[i]
    r = solve(A, y)
    return r

@cython.boundscheck(False)
cpdef double int_f2(ndarray[double, mode="c"] w,
        ndarray[double, mode="c"] values):
    cdef double r=0
    cdef unsigned n = len(w)
    cdef unsigned i
    for i in range(n):
        r += (values[i]**2)*w[i]
    return r

@cython.boundscheck(False)
cpdef double int_f2_f2(ndarray[double, mode="c"] w,
        ndarray[double, mode="c"] values1,
        ndarray[double, mode="c"] values2,
        ):
    cdef double r=0
    cdef unsigned n = len(w)
    cdef unsigned i
    for i in range(n):
        r += (values1[i]**2 + values2[i]**2)*w[i]
    return r

# gauss points+weights on the reference domain
cdef list _gauss_points_reference = None

def init_gauss_points(N=50):
    print "Precalculating gauss points..."
    global _gauss_points_reference
    _gauss_points_reference = range(N)
    for n in range(1, N):
        x, w = p_roots(n)
        _gauss_points_reference[n] = (real(x), w)
    print "    Done."

@cython.boundscheck(False)
cpdef get_gauss_points_phys(double a, double b, int n):
    cdef double J = (b-a)/2.0
    cdef unsigned n_points, i
    cdef ndarray[double, mode="c"] x, w
    if _gauss_points_reference is None:
        init_gauss_points()
    x, w = _gauss_points_reference[n]
    n_points = len(x)
    cdef ndarray[double, mode="c"] x_phys=empty(n_points), w_phys=empty(n_points)
    for i in range(n_points):
        x_phys[i] = J*(x[i]+1) + a
        w_phys[i] = w[i]*J
    return x_phys, w_phys

cpdef ndarray[double, mode="c"] get_fekete_points_phys(int order, double a,
        double b):
    """
    Returns the array of fekete points in physical domain.

    'order' ... is the polynomial order of the element
    (a, b)  ... physical interval where we want the fekete points

    There are "order+1" fekete points.

    Examples::

    >>> from hermes1d.fekete._fekete import get_fekete_points_phys
    >>> get_fekete_points_phys(1, 5, 6)
    array([ 5.,  6.])
    >>> get_fekete_points_phys(2, 5, 6)
    array([ 5. ,  5.5,  6. ])
    >>> get_fekete_points_phys(3, 5, 6)
    array([ 5.       ,  5.2763932,  5.7236068,  6.       ])
    >>> get_fekete_points_phys(4, 5, 6)
    array([ 5.        ,  5.17267316,  5.5       ,  5.82732684,  6.        ])

    """
    if order not in points:
        return get_chebyshev_points_phys(order, a, b)
    cdef ndarray[double, mode="c"] fekete_points = array(points[order])
    cdef double J = (b-a)/2.
    for i in range(order+1):
        fekete_points[i] = J*(fekete_points[i]+1) + a
    return fekete_points

cpdef ndarray[double, mode="c"] get_chebyshev_points_phys(int order, double a,
        double b):
    """
    Returns the array of Chebyshev points in physical domain.

    'order' ... is the polynomial order of the element
    (a, b)  ... physical interval where we want the Chebyshev points

    There are "order+1" Chebyshev points.

    Examples::

    >>> from hermes1d.fekete._fekete import get_chebyshev_points_phys
    >>> get_chebyshev_points_phys(1, 5, 6)
    array([ 5.,  6.])
    >>> get_chebyshev_points_phys(2, 5, 6)
    array([ 5. ,  5.5,  6. ])
    >>> get_chebyshev_points_phys(3, 5, 6)
    array([ 5.  ,  5.25,  5.75,  6.  ])
    >>> get_chebyshev_points_phys(4, 5, 6)
    array([ 5.        ,  5.14644661,  5.5       ,  5.85355339,  6.        ])

    """
    cdef unsigned N = order+1
    k = arange(N)
    cdef ndarray[double, mode="c"] points = -cos(pi*k/(N-1))
    cdef double J = (b-a)/2.
    for i in range(N):
        points[i] = J*(points[i]+1) + a
    return points

@cython.boundscheck(False)
def eval_poly(ndarray[double, mode="c"] x not None,
        ndarray[double, mode="c"] values not None, double a, double b):
    """
    Evaluates polynomial at points 'x'.

    The polynomial is defined by 'values' in fekete points, on the interval
    (a, b).
    """
    cdef unsigned n = len(x)
    cdef unsigned n_fekete = len(values)
    cdef unsigned i
    cdef ndarray[double, mode="c"] fekete_points = \
            get_fekete_points_phys(n_fekete-1, a, b)
    assert len(fekete_points) == n_fekete
    cdef ndarray[double, mode="c"] y = empty(n)
    for i in range(n):
        y[i] = lagrange_interpolating_polynomial(&(fekete_points[0]),
                &(values[0]), n_fekete, x[i])
    return y

cdef double lagrange_interpolating_polynomial(double *pos, double *val,
        int degree, double x):
    cdef double weight
    cdef double y = 0
    for i in range(degree):
        weight = 1
        for j in range(degree):
            # The i-th term has to be skipped
            if j != i:
                weight *= (x - pos[j]) / (pos[i] - pos[j])
        y += weight * val[i]
    return y
