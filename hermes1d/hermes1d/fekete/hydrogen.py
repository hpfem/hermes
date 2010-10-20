"""
This module will eventually end up in SymPy itself.

For now we'll keep it here, until we polish the interface and get the job done.
"""

from sympy import var, factorial, sqrt, exp, S

def laguerre_poly(alpha, n, x):
    r = 0
    for m in range(n+1):
        c = 1
        for i in range(m+1, n+1):
            c *= alpha+i
        r += (-1)**m * c * x**m/(factorial(m)*factorial(n-m))
    return r

def R_nl(n, l, r, a):
    """
    Returns the radial Hydrogen wavefunction, as a function of "r".

    "a" is the Bohr radius.
    """
    return S(1)/n * sqrt(factorial(n-l-1)/factorial(n+l)/a**3) * (2*r/n/a)**l * \
            laguerre_poly(2*l+1, n-l-1, 2*r/n/a)*exp(-r/n/a)

def R_nl_numeric(n, l, x):
    """
    Returns the exact floating point value of the R_nl at the (float) point
    'x'.
    """
    return float(R_nl(n, l, x, 1).evalf())

if __name__ == "__main__":
    var("alpha z")
    print laguerre_poly(alpha, 3, z)
    var("r")
    for n in range(1, 7):
        for l in range(0, n):
            print "(%d, %d): %s" % (n, l, R_nl(n, l, r, 1))

    print
    print R_nl(4, 0, S(1)/2, 1).evalf()
    print R_nl_numeric(4, 0, 0.5)
    #from pylab import plot, show
    #from numpy import arange
    #x = arange(0, 20, 0.1)
    #y = [R_nl(4, 0, r, 1) for r in x]
    #plot(x, y)
    #show()
