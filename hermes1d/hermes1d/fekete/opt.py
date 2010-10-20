"""
A simple script for optimizing and profiling fekete.py.
"""
from numpy import array
from math import exp
from fekete import Function, Mesh1D

pts = array([   0.,           0.5859375,    1.171875,     2.34375,      4.6875,
    7.03125, 9.375,       11.71875,     14.0625,      18.75,        37.5, 75.,
    112.5,        150.       ])
orders = [12,  9, 16, 16,  6,  5,  5,  4,  1,  1,  1,  1,  1]

@profile
def p(a, b):
    diff = a-b
    return diff.l2_norm()

print "mesh"
m = Mesh1D(pts, orders)
m_ref = m.refine_all_elements()
m_ref = m_ref.increase_order()
f = Function(lambda x: exp(-x)*(x*x+1), m)
f_ref = Function(lambda x: exp(-x)*(x*x+1), m_ref)
print f
print f_ref
print "loop"
for i in range(100):
    l2 = p(f, f_ref)
print l2
