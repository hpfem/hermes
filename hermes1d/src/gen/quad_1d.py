# Generates 1D Gauss quadrature points and 
# weights in (-1,1).

from scipy.special.orthogonal import legendre
from scipy.special.orthogonal import p_roots
from scipy.integrate import quadrature

# This script prints Gauss quadrature points and
# weights for a given number n of points.
 
n = 101
print "n =", n
print "orders: %d, %d" % (2*n - 2, 2*n - 1)
roots, weights = p_roots(n)
print "{"
for root, weight in zip(roots, weights):
    print "  { %.20f, %.20f }," % (root, weight)
