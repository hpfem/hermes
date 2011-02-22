#!/usr/bin/env python
from numpy import *
import sys
def list_prod(*args):
    " returns the product of an arbitrary number of lists"
    ndim=len(args)
    tmp= [(i,"args[%d]" % i) for i in range(ndim)]
    fs=range(ndim)+list(reduce(lambda x,y:x+y,tmp))
    fs=tuple(fs)
    ex=("[["+ndim*"x%d,"+"]"+ndim*"for x%d in %s "+"]") %  fs
    return eval(ex)
hex_order=[
[0,  0,  0],
[1,  0,  0],
[1,  1,  0],
[0,  1,  0],
[0,  0,  1],
[1,  0,  1],
[1,  1,  1],
[0,  1,  1]
]
standard_quads=[
[1, 2, 3, 4],
[1, 2, 6, 5],
[2, 3, 7, 6],
[3, 4, 8, 7],
[4, 1, 5, 8],
[5, 6, 7, 8]
]
xx=range(3)
xx[0],xx[1],xx[2]=map(eval,sys.argv[1:4])
nx,ny,nz=map(len,xx)
x=range(3)
x[0]=range(nx)
x[1]=range(ny)
x[2]=range(nz)
def bdcon(point):
    x,y,z=point
    if  ((x == 0) or ( x == nx-1) or \
             (y == 0) or ( y == ny-1) or \
             (z == 0 ) or ( z == nz-1)): 
        return 1
    else:
        return 0
vert=list_prod(x[0],x[1],x[2])
v0=list_prod(x[0][:-1],x[1][:-1],x[2][:-1])
print "# vertices"
print len(vert)
for v in vert:
    print "%25.16f %25.16f %25.16f" % (xx[0][v[0]],xx[1][v[1]],xx[2][v[2]])
print 
print "# tetras"
print 0
print 
print "# hexes"
print len(v0)
quads=[]
oquads=[]
for v in v0:
    p=range(8)
    for i in [0,1]:
        for j in [0,1]:
            for k in [0,1]:
                l=hex_order.index([i,j,k])
                p[l]=array(v)+array([i,j,k])
    hi=[]
    for pt in p:
        tmp=pt.tolist()
        ind=vert.index(tmp)+1
        hi.append(ind)
    print "%d %d %d %d %d %d %d %d" % tuple(hi)
    for [i1,i2,i3,i4] in standard_quads:
        q=[hi[i1-1],hi[i2-1],hi[i3-1],hi[i4-1]]
        qtmp=q[:]
        qtmp.sort()
        if qtmp not in oquads:
            quads.append(q)
            oquads.append(qtmp)
print\
"""
# prisms
0 

# tris
0 

#quads
"""
print len(quads)
for q in quads:
    # check boundary condition for all four points of quad
    p1,p2,p3,p4=vert[q[0]-1],vert[q[1]-1],vert[q[2]-1],vert[q[3]-1]
    # logical and gives final marker
    bm=bdcon(p1)*bdcon(p2)*bdcon(p3)*bdcon(p4)
    print "%d %d %d %d " % tuple(q),bm 


