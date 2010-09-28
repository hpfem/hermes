from sympy import var, sqrt, sin, exp, cos, pprint, pi, simplify, ccode
var("x y EPSILON")

def laplace(f):
    return f.diff(x, 2) + f.diff(y, 2)

# exact solution (given)
u = (1 - exp(-(1-x)/EPSILON)) * (1 - exp(-(1-y)/EPSILON)) * cos(pi * (x + y))
# right-hand-side (to be calculated)
rhs = -EPSILON * laplace(u) + 2*u.diff(x, 1) + u.diff(y, 1)
print "rhs, as a formula:"
pprint(rhs)
print
print "-"*60
print "rhs, as C code:"
print ccode(rhs)
# x-derivative of u (to be calculated)
dudx = u.diff(x, 1)
print "dudx, as C code:"
print ccode(dudx)
# y-derivative of u (to be calculated)
dudy = u.diff(y, 1)
print "dudy, as C code:"
print ccode(dudy)
