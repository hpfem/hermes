#! /usr/bin/env python

import os

from jinja2 import Environment, FileSystemLoader

from sympy import var, pprint, ccode
from orthogonalization import (gram_schmidt, l2_inner_product,
        h1_inner_product, integrate)

N = 20
precision = 25

def check(basis):
    print "orthonormality matrix:"
    for f in basis:
        for g in basis:
            print integrate(f*g+f.diff(x)*g.diff(x), (x, -1, 1)),
        print

var("x")
l2_ip = lambda f, g: l2_inner_product(f, g, (x, -1, 1))
h1_ip = lambda f, g: h1_inner_product(f, g, (x, -1, 1))

polys = [x**i for i in range(N)]
print "Applying Gram Schmidt process..."
h1_basis = gram_schmidt(polys, inner_product=h1_ip)
print "Orthonormal basis:"
#print h1_basis
#check(h1_basis)

functions = []
for i, v in enumerate(h1_basis):
    f = v.n(precision)
    f_diff = v.diff(x).n(precision)
    print "i:", f
    functions.append({"id": i,
        "expr": ccode(f),
        "expr_diff": ccode(f_diff),
        })

print "Generating the C file..."
template = "h1_polys.cpp"
env = Environment(loader=FileSystemLoader('.'))
t = env.get_template(template)
open(os.path.join("..", template), "w").write(t.render({
    "functions": functions,
    }))
