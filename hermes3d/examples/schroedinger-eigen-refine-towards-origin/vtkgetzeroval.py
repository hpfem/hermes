#!/usr/bin/env python
from numpy import *
from sys import argv
vtkfname=argv[1]
f=open(vtkfname,"r")
h=f.readline()
f.readline()
type=f.readline()
f.readline()
dset=f.readline()
l=f.readline()
npts=eval(l.split(" ")[1])
pts=[]
for i in range(npts):
    l=f.readline()
    point=tuple(map(eval,l.split(" ")))
    pts.append(point)
line="xxxxxxxxxxxxxxx"
while(line[:12] != "LOOKUP_TABLE"):
    line = f.readline()
F={}
i=0
for l in f:
    F[pts[i]]=eval(l)
    i=i+1
f.close()
print F[(0.0,0.0,0.0)]

