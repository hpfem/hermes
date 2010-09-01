#! /usr/bin/env python

print "Importing..."
import os
from jinja2 import Environment, FileSystemLoader
from sympy import Symbol

from common import lobatto, horner_scheme, ccode_pow2

n_functions = 50
precision = 25
factor_const = True

x = Symbol("x")
env = Environment(loader=FileSystemLoader('.'))

functions = []
print "Calculating shape functions..."
for i in range(n_functions):
    print "  i=%d" % i
    lob = lobatto(i, x)
    lob_diff = lob.diff(x)
    lob = horner_scheme(lob, x, factor_const=factor_const)
    print lob
    lob_diff = horner_scheme(lob_diff, x, factor_const=factor_const)
    functions.append({"id": i,
        "expr": ccode_pow2(lob),
        "expr_diff": ccode_pow2(lob_diff),
        })

print "Generating the C file..."
template = "lobatto.cpp"
t = env.get_template(template)
open(os.path.join("..", template), "w").write(t.render({
    "functions": functions,
    }))
