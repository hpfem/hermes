"""
SymPy orthogonalization module. This should become part of SymPy, after better
docstrings (+doctests) as well as tests are written.

It should probably go into sympy/matrices/
"""

from sympy import sqrt, integrate, conjugate, Symbol, sympify

def l2_inner_product(a, b, lim):
    """
    Calculates the L2 inner product (a, b) over the domain lim.
    """
    return integrate(a*b, lim)

def h1_inner_product(a, b, lim):
    """
    Calculates the L2 inner product (a, b) over the domain lim.
    """
    x, _, _ = lim
    a = sympify(a)
    b = sympify(b)
    return integrate(a*b + a.diff(x)*b.diff(x), lim)


def norm_induced(v, inner_product):
    return sqrt(inner_product(v, v))

def normalize(v, inner_product):
    return v/norm_induced(v, inner_product)

def projection(v, basis, inner_product):
    """
    Projects the function v on the basis using the inner_product.
    """
    r = 0
    for b in basis:
        r += inner_product(v, b) * b
    return r

def gram_schmidt(basis, inner_product=None):
    """
    Orthonormalizes the "basis" using the Gram-Schmidt process.

    basis ........... an arbitrary set of linearly independent SymPy expressions
    inner_product ... the inner product to use for projections

    Example:
    >>> l2_gram_schmidt([1, x, x**2], (x, -1, 1)]
    [1/2*2**(1/2), x*6**(1/2)/2, -3*10**(1/2)*(1/3 - x**2)/4]

    """
    if inner_product is None:
        x = basis[-1].atoms(Symbol).pop()
        inner_product = lambda f, g: l2_inner_product(f, g, (x, -1, 1))

    r = []
    for a in basis:
        if r == []:
            v = a
        else:
            v = a - projection(a, r, inner_product)
        v = normalize(v, inner_product)
        r.append(v)
    return r
