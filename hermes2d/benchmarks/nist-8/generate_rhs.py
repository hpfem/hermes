from sympy import var, sqrt, sin, pprint, simplify, ccode
var("x y ALPHA")

def laplace(f):
    return f.diff(x, 2) + f.diff(y, 2)

r = sqrt(x**2+y**2)
u = sin(1/(ALPHA+r))
lhs = -laplace(u) - u/(ALPHA+r)**4
print "lhs, as a formula:"
pprint(lhs)
print
print "-"*60
print "lhs, as C code:"
print ccode(lhs)
