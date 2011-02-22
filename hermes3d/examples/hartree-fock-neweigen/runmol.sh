#!/bin/bash
molcode.py < $1.dat > $1.cpp
\rm mol.cpp
\rm mol.mesh3d
ln -s $1.cpp mol.cpp
ln -s $1.mesh3d mol.mesh3d
make 
molecule3d-eigen 





