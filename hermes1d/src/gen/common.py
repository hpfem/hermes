"""
Common functions for all basis generators.
"""

from sympy import Symbol, legendre, sqrt, integrate, ccode, Poly
from sympy.mpmath import polyroots, nstr, workdps

def ccode_pow2(s):
    s = ccode(s)
    s = s.replace("pow(x,2)", "((x)*(x))")
    return s

def legendre_int(i, x):
    """
    Returns the normalized integrated Legendre polynomial.
    """
    y = Symbol("y", dummy=True)
    f = legendre(i, y)
    n = sqrt(integrate(f**2, (y, -1, 1)))
    return integrate(f, (y, -1, x))/n

def legendre_norm(i, x):
    """
    Returns the normalized integrated Legendre polynomial.
    """
    f = legendre(i, x)
    n = sqrt(integrate(f**2, (x, -1, 1)))
    return f/n

def legendre_shape_function(i, x):
    """
    Returns the Lobatto shape function.
    """
    f = legendre_norm(i, x)
    return f.expand()

def lobatto(i, x):
    """
    Returns the Lobatto shape function.
    """
    if i == 0:
        return 1-lobatto(1, x)
    f = legendre_int(i-1, x)
    if i == 1:
        f /= sqrt(2)
    return f.expand()

def horner_scheme(p, x, factor_const=False):
    """
    Rewrites the polynomial using the Horner scheme.
    """
    a = p.subs(x, 0)
    if p == a:
        return p
    rest = ((p-a)/x).expand()
    if factor_const:
        const = rest.subs(x, 0)
        if const != 0:
            rest = (rest/const).expand()
        else:
            const = 1
    else:
        const = 1
    return const*x*horner_scheme(rest, x, factor_const)+a

def gauss_lobatto_points(p):
    """
    Returns the list of Gauss-Lobatto points of the order 'p'.
    """
    x = Symbol("x")
    e = (1-x**2)*legendre(p, x).diff(x)
    e = Poly(e, x)
    if e == 0:
        return []
    with workdps(40):
        r, err = polyroots(e.all_coeffs(), error=True)
    if err > 1e-40:
        raise Exception("Internal Error: Root is not precise")
    p = []
    for x in r:
        if abs(x) < 1e-40: x = 0
        p.append(str(x))
    return p
