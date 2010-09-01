#! /usr/bin/env python

print "Importing..."
import os
from jinja2 import Environment, FileSystemLoader

from common import lobatto, horner_scheme, ccode_pow2
from common import gauss_lobatto_points

N = 13

print "Generating points"
points = []
for i in range(N):
    print i
    p = gauss_lobatto_points(i)
    points.append({
        "p": i,
        "points": p
        })

print points

print "Generating the Python file..."
template = "gauss_lobatto_points.py"
env = Environment(loader=FileSystemLoader('.'))
t = env.get_template(template)
open(os.path.join("..", template), "w").write(t.render({
    "data": points,
    }))
