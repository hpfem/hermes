from sympy import var, sqrt, sin, pprint, simplify, ccode, atan2, atan, exp, S
var("x y alpha_w alpha_p omega pi x_w y_w x_p y_p epsilon r0")

def laplace(f):
    return f.diff(x, 2) + f.diff(y, 2)

r = sqrt(x**2+y**2)
r_w = sqrt((x-x_w)**2+(y-y_w)**2)
r_p_squared = (x-x_p)**2+(y-y_p)**2
theta = atan(y/x)
u = r**(pi/omega) * sin(theta * pi / omega) + atan(alpha_w * (r_w - r0)) + exp(-alpha_p*r_p_squared) + exp(-(1+y)/epsilon)
params = {
    omega: 3*pi/2,
    x_w: 0,
    y_w: -S(3)/4,
    r0: S(3)/4,
    alpha_p: 1000,
    alpha_w: 200,
    x_p: sqrt(5)/4,
    y_p: -S(1)/4,
    epsilon: S(1)/100,
    }
u = u.subs(params)
lhs = laplace(u)
	
print "lhs, as a formula:"
pprint(lhs)
print
print "-"*60
print "lhs, as C code:"
print ccode(lhs)
