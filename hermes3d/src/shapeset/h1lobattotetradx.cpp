// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

//
// h1lobbatotetradx.cc
//
//

#include "../h3dconfig.h"
#include "common.h"
#include "lobatto.h"

#ifdef WITH_TETRA

// Derivatives of the shape functions (dx)

// DEGREE 1
//------------

// Vertex shape functions, degree 1

double lobatto_dx_f0(double x, double y, double z) {
	return lambda1dx(x, y, z);
}

double lobatto_dx_f1(double x, double y, double z) {
	return lambda2dx(x, y, z);
}

double lobatto_dx_f2(double x, double y, double z) {
	return lambda0dx(x, y, z);
}

double lobatto_dx_f3(double x, double y, double z) {
	return lambda3dx(x, y, z);
}

// DEGREE 2
//------------

// Edge shape functions, degree 2

// edge 0
double lobatto_dx_f4(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z));
}

// edge 1
double lobatto_dx_f5(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z));
}

// edge 2
double lobatto_dx_f6(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z));
}

// edge 3
double lobatto_dx_f7(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z));
}

// edge 4
double lobatto_dx_f8(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z));
}

// edge 5
double lobatto_dx_f9(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z));
}


// DEGREE 3
//------------

// Edge shape functions, degree 3

// edge 0
double lobatto_dx_f10_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z));
}

double lobatto_dx_f10_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)));
}

// edge 1
double lobatto_dx_f11_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z));
}

double lobatto_dx_f11_1(double x, double y, double z) {
	return -(
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)));
}

// edge 2
double lobatto_dx_f12_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z));
}

double lobatto_dx_f12_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)));
}

// edge 3
double lobatto_dx_f13_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f13_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)));
}

// edge 4
double lobatto_dx_f14_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f14_1(double x, double y, double z) {
	return -(
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)));
}

// edge 5
double lobatto_dx_f15_0(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f15_1(double x, double y, double z) {
	return -(
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)));
}


// Face shape functions, degree 3

// face 0
double lobatto_dx_f16_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f16_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f16_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f16_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f16_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f16_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 1
double lobatto_dx_f17_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f17_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f17_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f17_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f17_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f17_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}


// face 2
double lobatto_dx_f18_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f18_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f18_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f18_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f18_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f18_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 3
double lobatto_dx_f19_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f19_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f19_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f19_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f19_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f19_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}


// DEGREE 4
//------------

// Edge shape functions, degree 4

// edge 0
double lobatto_dx_f20(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z));
}

// edge 1
double lobatto_dx_f21(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z));
}

// edge 2
double lobatto_dx_f22(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z));
}

// edge 3
double lobatto_dx_f23(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z));
}

// edge 4
double lobatto_dx_f24(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z));
}

// edge 5
double lobatto_dx_f25(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z));
}


// Face shape functions, degree 4

// face 0
double lobatto_dx_f26_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f26_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f26_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f26_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f26_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f26_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f27_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f27_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f27_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f27_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f27_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f27_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 1
double lobatto_dx_f28_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f28_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f28_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f28_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f28_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f28_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f29_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f29_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f29_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f29_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f29_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f29_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}


// face 2
double lobatto_dx_f30_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f30_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f30_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f30_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f30_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f30_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f31_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f31_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f31_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f31_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f31_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f31_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 3
double lobatto_dx_f32_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f32_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f32_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f32_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f32_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f32_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f33_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f33_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f33_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f33_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f33_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f33_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}


// Bubble shape functions, degree 4

double lobatto_dx_f34(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

// DEGREE 5
//------------

// Edge shape functions, degree 5

// edge 0
double lobatto_dx_f35_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z));
}

double lobatto_dx_f35_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)));
}

// edge 1
double lobatto_dx_f36_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z));
}

double lobatto_dx_f36_1(double x, double y, double z) {
	return -(
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)));
}

// edge 2
double lobatto_dx_f37_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z));
}

double lobatto_dx_f37_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)));
}

// edge 3
double lobatto_dx_f38_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f38_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)));
}

// edge 4
double lobatto_dx_f39_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f39_1(double x, double y, double z) {
	return -(
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)));
}

// edge 5
double lobatto_dx_f40_0(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f40_1(double x, double y, double z) {
	return -(
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)));
}


// Face shape functions, degree 5

// face 0
double lobatto_dx_f41_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f41_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f41_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f41_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f41_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f41_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f42_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f42_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f42_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f42_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f42_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f42_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f43_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f43_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f43_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f43_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f43_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f43_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 1
double lobatto_dx_f44_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f44_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f44_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f44_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f44_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f44_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f45_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f45_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f45_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f45_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f45_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f45_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f46_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f46_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f46_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f46_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f46_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f46_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}


// face 2
double lobatto_dx_f47_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f47_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f47_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f47_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f47_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f47_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f48_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f48_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f48_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f48_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f48_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f48_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f49_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f49_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f49_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f49_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f49_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f49_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 3
double lobatto_dx_f50_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f50_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f50_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f50_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f50_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f50_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f51_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f51_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f51_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f51_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f51_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f51_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f52_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f52_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f52_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f52_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f52_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f52_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}


// Bubble shape functions, degree 5

double lobatto_dx_f53(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f54(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f55(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

// DEGREE 6
//------------

// Edge shape functions, degree 6

// edge 0
double lobatto_dx_f56(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z));
}

// edge 1
double lobatto_dx_f57(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z));
}

// edge 2
double lobatto_dx_f58(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z));
}

// edge 3
double lobatto_dx_f59(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z));
}

// edge 4
double lobatto_dx_f60(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z));
}

// edge 5
double lobatto_dx_f61(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z));
}


// Face shape functions, degree 6

// face 0
double lobatto_dx_f62_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f62_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f62_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f62_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f62_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f62_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f63_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f63_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f63_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f63_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f63_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f63_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f64_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f64_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f64_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f64_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f64_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f64_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f65_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f65_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f65_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f65_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f65_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f65_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 1
double lobatto_dx_f66_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f66_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f66_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f66_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f66_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f66_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f67_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f67_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f67_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f67_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f67_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f67_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f68_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f68_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f68_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f68_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f68_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f68_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f69_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f69_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f69_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f69_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f69_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f69_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}


// face 2
double lobatto_dx_f70_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f70_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f70_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f70_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f70_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f70_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f71_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f71_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f71_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f71_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f71_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f71_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f72_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f72_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f72_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f72_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f72_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f72_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f73_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f73_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f73_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f73_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f73_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f73_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 3
double lobatto_dx_f74_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f74_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f74_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f74_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f74_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f74_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f75_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f75_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f75_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f75_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f75_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f75_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f76_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f76_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f76_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f76_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f76_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f76_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f77_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f77_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f77_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f77_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f77_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f77_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}


// Bubble shape functions, degree 6

double lobatto_dx_f78(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f79(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f80(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f81(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f82(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f83(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

// DEGREE 7
//------------

// Edge shape functions, degree 7

// edge 0
double lobatto_dx_f84_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z));
}

double lobatto_dx_f84_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)));
}

// edge 1
double lobatto_dx_f85_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z));
}

double lobatto_dx_f85_1(double x, double y, double z) {
	return -(
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)));
}

// edge 2
double lobatto_dx_f86_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z));
}

double lobatto_dx_f86_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)));
}

// edge 3
double lobatto_dx_f87_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f87_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)));
}

// edge 4
double lobatto_dx_f88_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f88_1(double x, double y, double z) {
	return -(
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)));
}

// edge 5
double lobatto_dx_f89_0(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f89_1(double x, double y, double z) {
	return -(
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)));
}


// Face shape functions, degree 7

// face 0
double lobatto_dx_f90_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f90_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f90_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f90_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f90_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f90_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f91_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f91_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f91_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f91_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f91_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f91_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f92_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f92_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f92_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f92_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f92_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f92_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f93_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f93_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f93_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f93_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f93_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f93_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f94_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f94_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f94_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f94_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f94_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f94_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 1
double lobatto_dx_f95_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f95_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f95_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f95_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f95_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f95_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f96_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f96_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f96_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f96_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f96_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f96_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f97_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f97_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f97_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f97_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f97_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f97_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f98_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f98_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f98_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f98_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f98_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f98_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f99_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f99_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f99_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f99_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f99_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f99_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}


// face 2
double lobatto_dx_f100_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f100_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f100_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f100_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f100_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f100_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f101_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f101_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f101_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f101_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f101_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f101_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f102_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f102_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f102_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f102_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f102_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f102_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f103_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f103_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f103_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f103_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f103_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f103_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f104_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f104_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f104_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f104_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f104_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f104_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 3
double lobatto_dx_f105_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f105_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f105_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f105_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f105_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f105_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f106_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f106_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f106_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f106_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f106_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f106_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f107_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f107_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f107_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f107_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f107_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f107_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f108_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f108_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f108_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f108_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f108_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f108_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f109_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f109_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f109_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f109_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f109_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f109_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}


// Bubble shape functions, degree 7

double lobatto_dx_f110(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f111(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f112(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f113(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f114(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f115(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f116(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f117(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f118(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f119(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

// DEGREE 8
//------------

// Edge shape functions, degree 8

// edge 0
double lobatto_dx_f120(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi6dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z));
}

// edge 1
double lobatto_dx_f121(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi6dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z));
}

// edge 2
double lobatto_dx_f122(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi6dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z));
}

// edge 3
double lobatto_dx_f123(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi6dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z));
}

// edge 4
double lobatto_dx_f124(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi6dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z));
}

// edge 5
double lobatto_dx_f125(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi6dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z));
}


// Face shape functions, degree 8

// face 0
double lobatto_dx_f126_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f126_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f126_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f126_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f126_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f126_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f127_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f127_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f127_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f127_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f127_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f127_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f128_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f128_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f128_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f128_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f128_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f128_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f129_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f129_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f129_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f129_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f129_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f129_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f130_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f130_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f130_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f130_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f130_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f130_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f131_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f131_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f131_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f131_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f131_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f131_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 1
double lobatto_dx_f132_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f132_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f132_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f132_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f132_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f132_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f133_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f133_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f133_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f133_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f133_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f133_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f134_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f134_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f134_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f134_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f134_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f134_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f135_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f135_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f135_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f135_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f135_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f135_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f136_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f136_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f136_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f136_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f136_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f136_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f137_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f137_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f137_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f137_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f137_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f137_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}


// face 2
double lobatto_dx_f138_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f138_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f138_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f138_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f138_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f138_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f139_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f139_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f139_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f139_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f139_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f139_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f140_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f140_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f140_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f140_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f140_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f140_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f141_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f141_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f141_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f141_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f141_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f141_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f142_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f142_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f142_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f142_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f142_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f142_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f143_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f143_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f143_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f143_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f143_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f143_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 3
double lobatto_dx_f144_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f144_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f144_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f144_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f144_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f144_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f145_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f145_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f145_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f145_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f145_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f145_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f146_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f146_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f146_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f146_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f146_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f146_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f147_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f147_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f147_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f147_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f147_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f147_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f148_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f148_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f148_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f148_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f148_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f148_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f149_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f149_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f149_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f149_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f149_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f149_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}


// Bubble shape functions, degree 8

double lobatto_dx_f150(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f151(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f152(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f153(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f154(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f155(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f156(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f157(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f158(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f159(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f160(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f161(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f162(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f163(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f164(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

// DEGREE 9
//------------

// Edge shape functions, degree 9

// edge 0
double lobatto_dx_f165_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi7dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z));
}

double lobatto_dx_f165_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi7dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)));
}

// edge 1
double lobatto_dx_f166_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi7dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z));
}

double lobatto_dx_f166_1(double x, double y, double z) {
	return -(
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi7dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)));
}

// edge 2
double lobatto_dx_f167_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi7dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z));
}

double lobatto_dx_f167_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi7dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)));
}

// edge 3
double lobatto_dx_f168_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi7dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f168_1(double x, double y, double z) {
	return -(
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi7dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)));
}

// edge 4
double lobatto_dx_f169_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi7dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f169_1(double x, double y, double z) {
	return -(
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi7dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)));
}

// edge 5
double lobatto_dx_f170_0(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi7dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z));
}

double lobatto_dx_f170_1(double x, double y, double z) {
	return -(
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi7dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)));
}


// Face shape functions, degree 9

// face 0
double lobatto_dx_f171_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f171_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f171_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f171_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f171_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f171_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f172_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f172_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f172_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f172_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f172_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f172_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f173_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f173_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f173_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f173_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f173_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f173_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f174_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f174_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f174_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f174_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f174_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f174_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f175_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f175_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f175_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f175_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f175_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f175_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f176_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f176_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f176_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f176_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f176_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f176_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f177_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f177_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f177_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f177_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f177_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f177_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 1
double lobatto_dx_f178_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f178_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f178_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f178_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f178_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f178_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f179_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f179_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f179_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f179_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f179_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f179_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f180_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f180_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f180_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f180_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f180_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f180_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f181_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f181_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f181_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f181_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f181_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f181_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f182_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f182_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f182_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f182_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f182_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f182_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f183_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f183_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f183_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f183_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f183_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f183_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f184_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f184_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f184_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f184_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f184_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f184_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}


// face 2
double lobatto_dx_f185_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f185_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f185_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f185_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f185_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f185_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f186_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f186_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f186_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f186_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f186_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f186_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f187_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f187_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f187_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f187_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f187_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f187_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f188_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f188_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f188_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f188_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f188_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f188_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f189_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f189_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f189_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f189_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f189_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f189_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f190_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f190_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f190_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f190_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f190_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f190_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f191_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f191_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f191_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f191_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f191_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f191_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 3
double lobatto_dx_f192_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f192_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f192_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f192_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f192_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f192_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f193_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f193_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f193_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f193_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f193_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f193_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f194_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f194_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f194_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f194_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f194_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f194_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f195_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f195_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f195_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f195_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f195_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f195_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f196_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f196_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f196_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f196_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f196_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f196_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f197_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f197_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f197_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f197_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f197_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f197_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f198_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f198_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f198_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f198_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f198_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f198_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}


// Bubble shape functions, degree 9

double lobatto_dx_f199(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f200(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f201(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f202(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f203(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f204(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f205(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f206(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f207(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f208(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f209(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f210(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f211(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f212(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f213(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f214(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f215(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f216(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f217(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f218(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f219(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

// DEGREE 10
//------------

// Edge shape functions, degree 10

// edge 0
double lobatto_dx_f220(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * phi8(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * phi8(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * phi8dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z));
}

// edge 1
double lobatto_dx_f221(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * phi8(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * phi8(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * phi8dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z));
}

// edge 2
double lobatto_dx_f222(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * phi8(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * phi8(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * phi8dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z));
}

// edge 3
double lobatto_dx_f223(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * phi8(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * phi8(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * phi8dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z));
}

// edge 4
double lobatto_dx_f224(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * phi8(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * phi8(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * phi8dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z));
}

// edge 5
double lobatto_dx_f225(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * phi8(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * phi8(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * phi8dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z));
}


// Face shape functions, degree 10

// face 0
double lobatto_dx_f226_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi7dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f226_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi7dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f226_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi7dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f226_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi7dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f226_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi7dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f226_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi7dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f227_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f227_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f227_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f227_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f227_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f227_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f228_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f228_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f228_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f228_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f228_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f228_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f229_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f229_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f229_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f229_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f229_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f229_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f230_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f230_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f230_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f230_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f230_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f230_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f231_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f231_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f231_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f231_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f231_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f231_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f232_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f232_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f232_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f232_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f232_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f232_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f233_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi7dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f233_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi7dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f233_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi7dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f233_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi7dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f233_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi7dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f233_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi7dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 1
double lobatto_dx_f234_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi7dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f234_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi7dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f234_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) * phi7dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f234_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) * phi7dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f234_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi7dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f234_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi7dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f235_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f235_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f235_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f235_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f235_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f235_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f236_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f236_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f236_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f236_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f236_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f236_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f237_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f237_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f237_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f237_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f237_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f237_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f238_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f238_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f238_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f238_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f238_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f238_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f239_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f239_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f239_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f239_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f239_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f239_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f240_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f240_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f240_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f240_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f240_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f240_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f241_0(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi7dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda3(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f241_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda2(x, y, z) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2dx(x, y, z) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi7dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda2(x, y, z) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f241_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi7dx(lambda2(x, y, z) - lambda3(x, y, z)) * (lambda2dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi7(lambda2(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f241_3(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi7dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi7(lambda3(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f241_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi7dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f241_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi7dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda2(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda2(x, y, z)) * (lambda3dx(x, y, z) - lambda2dx(x, y, z))
;
}


// face 2
double lobatto_dx_f242_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi7dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f242_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) * phi7dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f242_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) * phi7dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f242_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) * phi7dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f242_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi7dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f242_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) * phi7dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f243_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f243_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f243_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f243_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f243_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f243_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f244_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f244_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f244_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f244_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f244_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f244_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f245_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f245_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f245_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f245_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f245_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f245_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f246_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f246_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi4(lambda3(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f246_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f246_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f246_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f246_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda3(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f247_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f247_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi5(lambda3(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f247_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f247_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f247_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f247_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda3(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f248_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f248_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi6(lambda3(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f248_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f248_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f248_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f248_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda3(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f249_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda3(x, y, z) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3dx(x, y, z) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi7dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda3(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda3(x, y, z) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f249_1(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3dx(x, y, z) * lambda1(x, y, z) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1dx(x, y, z) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi7dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda3(x, y, z) * lambda1(x, y, z) * phi7(lambda3(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f249_2(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi7dx(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda0(x, y, z)) +
		lambda3(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi7(lambda1(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda0(x, y, z)) * (lambda3dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f249_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3dx(x, y, z) * lambda0(x, y, z) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0dx(x, y, z) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi7dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda3(x, y, z) * lambda0(x, y, z) * phi7(lambda3(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f249_4(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda3(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3dx(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi7dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda3(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda3(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z))
;
}

double lobatto_dx_f249_5(double x, double y, double z) {
	return
		lambda3dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi7dx(lambda0(x, y, z) - lambda3(x, y, z)) * (lambda0dx(x, y, z) - lambda3dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda3(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi7(lambda0(x, y, z) - lambda3(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z))
;
}


// face 3
double lobatto_dx_f250_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi7dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f250_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) * phi7dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f250_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) * phi7dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f250_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi7dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f250_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) * phi7dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f250_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) * phi7dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f251_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f251_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f251_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f251_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f251_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f251_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) * phi6dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f252_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f252_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f252_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f252_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f252_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f252_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f253_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f253_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f253_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f253_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f253_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f253_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f254_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f254_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi4(lambda0(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f254_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi4(lambda1(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f254_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f254_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi4(lambda1(x, y, z) - lambda2(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f254_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi4(lambda2(x, y, z) - lambda0(x, y, z)) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f255_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f255_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi5(lambda0(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f255_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi5(lambda1(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f255_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f255_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi5(lambda1(x, y, z) - lambda2(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f255_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi5(lambda2(x, y, z) - lambda0(x, y, z)) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f256_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f256_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi6(lambda0(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f256_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi6(lambda1(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f256_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f256_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi6(lambda1(x, y, z) - lambda2(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f256_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi6(lambda2(x, y, z) - lambda0(x, y, z)) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f257_0(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2dx(x, y, z) * lambda0(x, y, z) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0dx(x, y, z) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi7dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda0(x, y, z)) +
		lambda1(x, y, z) * lambda2(x, y, z) * lambda0(x, y, z) * phi7(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f257_1(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0dx(x, y, z) * lambda1(x, y, z) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1dx(x, y, z) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi7dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) +
		lambda2(x, y, z) * lambda0(x, y, z) * lambda1(x, y, z) * phi7(lambda0(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z))
;
}

double lobatto_dx_f257_2(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi7dx(lambda1(x, y, z) - lambda0(x, y, z)) * (lambda1dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda2(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * phi7(lambda1(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda2(x, y, z)) * (lambda0dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f257_3(double x, double y, double z) {
	return
		lambda1dx(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0dx(x, y, z) * lambda2(x, y, z) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2dx(x, y, z) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi7dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) +
		lambda1(x, y, z) * lambda0(x, y, z) * lambda2(x, y, z) * phi7(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z))
;
}

double lobatto_dx_f257_4(double x, double y, double z) {
	return
		lambda2dx(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1dx(x, y, z) * lambda0(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0dx(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi7dx(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1dx(x, y, z) - lambda2dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda0(x, y, z)) +
		lambda2(x, y, z) * lambda1(x, y, z) * lambda0(x, y, z) * phi7(lambda1(x, y, z) - lambda2(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z))
;
}

double lobatto_dx_f257_5(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2dx(x, y, z) * lambda1(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1dx(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi7dx(lambda2(x, y, z) - lambda0(x, y, z)) * (lambda2dx(x, y, z) - lambda0dx(x, y, z)) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda2(x, y, z) * lambda1(x, y, z) * phi7(lambda2(x, y, z) - lambda0(x, y, z)) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z))
;
}


// Bubble shape functions, degree 10

double lobatto_dx_f258(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f259(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f260(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f261(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f262(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f263(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f264(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi6(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f265(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f266(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f267(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f268(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f269(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f270(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi1(lambda0(x, y, z) - lambda1(x, y, z)) * phi5(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f271(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f272(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f273(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f274(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f275(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi2(lambda0(x, y, z) - lambda1(x, y, z)) * phi4(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f276(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f277(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f278(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f279(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi3(lambda0(x, y, z) - lambda1(x, y, z)) * phi3(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f280(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f281(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f282(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi4(lambda0(x, y, z) - lambda1(x, y, z)) * phi2(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f283(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f284(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi5(lambda0(x, y, z) - lambda1(x, y, z)) * phi1(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}

double lobatto_dx_f285(double x, double y, double z) {
	return
		lambda0dx(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1dx(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2dx(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3dx(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2dx(x, y, z) - lambda1dx(x, y, z)) * phi0(lambda3(x, y, z) - lambda1(x, y, z)) +
		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi6(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi0dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3dx(x, y, z) - lambda1dx(x, y, z));
}


shape_fn_t lobatto_tetra_dx[] = {
	lobatto_dx_f0, lobatto_dx_f1, lobatto_dx_f2, lobatto_dx_f3, lobatto_dx_f4, lobatto_dx_f5, lobatto_dx_f6,
	lobatto_dx_f7, lobatto_dx_f8, lobatto_dx_f9, lobatto_dx_f10_0, lobatto_dx_f10_1, lobatto_dx_f11_0, lobatto_dx_f11_1,
	lobatto_dx_f12_0, lobatto_dx_f12_1, lobatto_dx_f13_0, lobatto_dx_f13_1, lobatto_dx_f14_0, lobatto_dx_f14_1, lobatto_dx_f15_0,
	lobatto_dx_f15_1, lobatto_dx_f16_0, lobatto_dx_f16_1, lobatto_dx_f16_2, lobatto_dx_f16_3, lobatto_dx_f16_4, lobatto_dx_f16_5,
	lobatto_dx_f17_0, lobatto_dx_f17_1, lobatto_dx_f17_2, lobatto_dx_f17_3, lobatto_dx_f17_4, lobatto_dx_f17_5, lobatto_dx_f18_0,
	lobatto_dx_f18_1, lobatto_dx_f18_2, lobatto_dx_f18_3, lobatto_dx_f18_4, lobatto_dx_f18_5, lobatto_dx_f19_0, lobatto_dx_f19_1,
	lobatto_dx_f19_2, lobatto_dx_f19_3, lobatto_dx_f19_4, lobatto_dx_f19_5, lobatto_dx_f20, lobatto_dx_f21, lobatto_dx_f22,
	lobatto_dx_f23, lobatto_dx_f24, lobatto_dx_f25, lobatto_dx_f26_0, lobatto_dx_f26_1, lobatto_dx_f26_2, lobatto_dx_f26_3,
	lobatto_dx_f26_4, lobatto_dx_f26_5, lobatto_dx_f27_0, lobatto_dx_f27_1, lobatto_dx_f27_2, lobatto_dx_f27_3, lobatto_dx_f27_4,
	lobatto_dx_f27_5, lobatto_dx_f28_0, lobatto_dx_f28_1, lobatto_dx_f28_2, lobatto_dx_f28_3, lobatto_dx_f28_4, lobatto_dx_f28_5,
	lobatto_dx_f29_0, lobatto_dx_f29_1, lobatto_dx_f29_2, lobatto_dx_f29_3, lobatto_dx_f29_4, lobatto_dx_f29_5, lobatto_dx_f30_0,
	lobatto_dx_f30_1, lobatto_dx_f30_2, lobatto_dx_f30_3, lobatto_dx_f30_4, lobatto_dx_f30_5, lobatto_dx_f31_0, lobatto_dx_f31_1,
	lobatto_dx_f31_2, lobatto_dx_f31_3, lobatto_dx_f31_4, lobatto_dx_f31_5, lobatto_dx_f32_0, lobatto_dx_f32_1, lobatto_dx_f32_2,
	lobatto_dx_f32_3, lobatto_dx_f32_4, lobatto_dx_f32_5, lobatto_dx_f33_0, lobatto_dx_f33_1, lobatto_dx_f33_2, lobatto_dx_f33_3,
	lobatto_dx_f33_4, lobatto_dx_f33_5, lobatto_dx_f34, lobatto_dx_f35_0, lobatto_dx_f35_1, lobatto_dx_f36_0, lobatto_dx_f36_1,
	lobatto_dx_f37_0, lobatto_dx_f37_1, lobatto_dx_f38_0, lobatto_dx_f38_1, lobatto_dx_f39_0, lobatto_dx_f39_1, lobatto_dx_f40_0,
	lobatto_dx_f40_1, lobatto_dx_f41_0, lobatto_dx_f41_1, lobatto_dx_f41_2, lobatto_dx_f41_3, lobatto_dx_f41_4, lobatto_dx_f41_5,
	lobatto_dx_f42_0, lobatto_dx_f42_1, lobatto_dx_f42_2, lobatto_dx_f42_3, lobatto_dx_f42_4, lobatto_dx_f42_5, lobatto_dx_f43_0,
	lobatto_dx_f43_1, lobatto_dx_f43_2, lobatto_dx_f43_3, lobatto_dx_f43_4, lobatto_dx_f43_5, lobatto_dx_f44_0, lobatto_dx_f44_1,
	lobatto_dx_f44_2, lobatto_dx_f44_3, lobatto_dx_f44_4, lobatto_dx_f44_5, lobatto_dx_f45_0, lobatto_dx_f45_1, lobatto_dx_f45_2,
	lobatto_dx_f45_3, lobatto_dx_f45_4, lobatto_dx_f45_5, lobatto_dx_f46_0, lobatto_dx_f46_1, lobatto_dx_f46_2, lobatto_dx_f46_3,
	lobatto_dx_f46_4, lobatto_dx_f46_5, lobatto_dx_f47_0, lobatto_dx_f47_1, lobatto_dx_f47_2, lobatto_dx_f47_3, lobatto_dx_f47_4,
	lobatto_dx_f47_5, lobatto_dx_f48_0, lobatto_dx_f48_1, lobatto_dx_f48_2, lobatto_dx_f48_3, lobatto_dx_f48_4, lobatto_dx_f48_5,
	lobatto_dx_f49_0, lobatto_dx_f49_1, lobatto_dx_f49_2, lobatto_dx_f49_3, lobatto_dx_f49_4, lobatto_dx_f49_5, lobatto_dx_f50_0,
	lobatto_dx_f50_1, lobatto_dx_f50_2, lobatto_dx_f50_3, lobatto_dx_f50_4, lobatto_dx_f50_5, lobatto_dx_f51_0, lobatto_dx_f51_1,
	lobatto_dx_f51_2, lobatto_dx_f51_3, lobatto_dx_f51_4, lobatto_dx_f51_5, lobatto_dx_f52_0, lobatto_dx_f52_1, lobatto_dx_f52_2,
	lobatto_dx_f52_3, lobatto_dx_f52_4, lobatto_dx_f52_5, lobatto_dx_f53, lobatto_dx_f54, lobatto_dx_f55, lobatto_dx_f56,
	lobatto_dx_f57, lobatto_dx_f58, lobatto_dx_f59, lobatto_dx_f60, lobatto_dx_f61, lobatto_dx_f62_0, lobatto_dx_f62_1,
	lobatto_dx_f62_2, lobatto_dx_f62_3, lobatto_dx_f62_4, lobatto_dx_f62_5, lobatto_dx_f63_0, lobatto_dx_f63_1, lobatto_dx_f63_2,
	lobatto_dx_f63_3, lobatto_dx_f63_4, lobatto_dx_f63_5, lobatto_dx_f64_0, lobatto_dx_f64_1, lobatto_dx_f64_2, lobatto_dx_f64_3,
	lobatto_dx_f64_4, lobatto_dx_f64_5, lobatto_dx_f65_0, lobatto_dx_f65_1, lobatto_dx_f65_2, lobatto_dx_f65_3, lobatto_dx_f65_4,
	lobatto_dx_f65_5, lobatto_dx_f66_0, lobatto_dx_f66_1, lobatto_dx_f66_2, lobatto_dx_f66_3, lobatto_dx_f66_4, lobatto_dx_f66_5,
	lobatto_dx_f67_0, lobatto_dx_f67_1, lobatto_dx_f67_2, lobatto_dx_f67_3, lobatto_dx_f67_4, lobatto_dx_f67_5, lobatto_dx_f68_0,
	lobatto_dx_f68_1, lobatto_dx_f68_2, lobatto_dx_f68_3, lobatto_dx_f68_4, lobatto_dx_f68_5, lobatto_dx_f69_0, lobatto_dx_f69_1,
	lobatto_dx_f69_2, lobatto_dx_f69_3, lobatto_dx_f69_4, lobatto_dx_f69_5, lobatto_dx_f70_0, lobatto_dx_f70_1, lobatto_dx_f70_2,
	lobatto_dx_f70_3, lobatto_dx_f70_4, lobatto_dx_f70_5, lobatto_dx_f71_0, lobatto_dx_f71_1, lobatto_dx_f71_2, lobatto_dx_f71_3,
	lobatto_dx_f71_4, lobatto_dx_f71_5, lobatto_dx_f72_0, lobatto_dx_f72_1, lobatto_dx_f72_2, lobatto_dx_f72_3, lobatto_dx_f72_4,
	lobatto_dx_f72_5, lobatto_dx_f73_0, lobatto_dx_f73_1, lobatto_dx_f73_2, lobatto_dx_f73_3, lobatto_dx_f73_4, lobatto_dx_f73_5,
	lobatto_dx_f74_0, lobatto_dx_f74_1, lobatto_dx_f74_2, lobatto_dx_f74_3, lobatto_dx_f74_4, lobatto_dx_f74_5, lobatto_dx_f75_0,
	lobatto_dx_f75_1, lobatto_dx_f75_2, lobatto_dx_f75_3, lobatto_dx_f75_4, lobatto_dx_f75_5, lobatto_dx_f76_0, lobatto_dx_f76_1,
	lobatto_dx_f76_2, lobatto_dx_f76_3, lobatto_dx_f76_4, lobatto_dx_f76_5, lobatto_dx_f77_0, lobatto_dx_f77_1, lobatto_dx_f77_2,
	lobatto_dx_f77_3, lobatto_dx_f77_4, lobatto_dx_f77_5, lobatto_dx_f78, lobatto_dx_f79, lobatto_dx_f80, lobatto_dx_f81,
	lobatto_dx_f82, lobatto_dx_f83, lobatto_dx_f84_0, lobatto_dx_f84_1, lobatto_dx_f85_0, lobatto_dx_f85_1, lobatto_dx_f86_0,
	lobatto_dx_f86_1, lobatto_dx_f87_0, lobatto_dx_f87_1, lobatto_dx_f88_0, lobatto_dx_f88_1, lobatto_dx_f89_0, lobatto_dx_f89_1,
	lobatto_dx_f90_0, lobatto_dx_f90_1, lobatto_dx_f90_2, lobatto_dx_f90_3, lobatto_dx_f90_4, lobatto_dx_f90_5, lobatto_dx_f91_0,
	lobatto_dx_f91_1, lobatto_dx_f91_2, lobatto_dx_f91_3, lobatto_dx_f91_4, lobatto_dx_f91_5, lobatto_dx_f92_0, lobatto_dx_f92_1,
	lobatto_dx_f92_2, lobatto_dx_f92_3, lobatto_dx_f92_4, lobatto_dx_f92_5, lobatto_dx_f93_0, lobatto_dx_f93_1, lobatto_dx_f93_2,
	lobatto_dx_f93_3, lobatto_dx_f93_4, lobatto_dx_f93_5, lobatto_dx_f94_0, lobatto_dx_f94_1, lobatto_dx_f94_2, lobatto_dx_f94_3,
	lobatto_dx_f94_4, lobatto_dx_f94_5, lobatto_dx_f95_0, lobatto_dx_f95_1, lobatto_dx_f95_2, lobatto_dx_f95_3, lobatto_dx_f95_4,
	lobatto_dx_f95_5, lobatto_dx_f96_0, lobatto_dx_f96_1, lobatto_dx_f96_2, lobatto_dx_f96_3, lobatto_dx_f96_4, lobatto_dx_f96_5,
	lobatto_dx_f97_0, lobatto_dx_f97_1, lobatto_dx_f97_2, lobatto_dx_f97_3, lobatto_dx_f97_4, lobatto_dx_f97_5, lobatto_dx_f98_0,
	lobatto_dx_f98_1, lobatto_dx_f98_2, lobatto_dx_f98_3, lobatto_dx_f98_4, lobatto_dx_f98_5, lobatto_dx_f99_0, lobatto_dx_f99_1,
	lobatto_dx_f99_2, lobatto_dx_f99_3, lobatto_dx_f99_4, lobatto_dx_f99_5, lobatto_dx_f100_0, lobatto_dx_f100_1, lobatto_dx_f100_2,
	lobatto_dx_f100_3, lobatto_dx_f100_4, lobatto_dx_f100_5, lobatto_dx_f101_0, lobatto_dx_f101_1, lobatto_dx_f101_2, lobatto_dx_f101_3,
	lobatto_dx_f101_4, lobatto_dx_f101_5, lobatto_dx_f102_0, lobatto_dx_f102_1, lobatto_dx_f102_2, lobatto_dx_f102_3, lobatto_dx_f102_4,
	lobatto_dx_f102_5, lobatto_dx_f103_0, lobatto_dx_f103_1, lobatto_dx_f103_2, lobatto_dx_f103_3, lobatto_dx_f103_4, lobatto_dx_f103_5,
	lobatto_dx_f104_0, lobatto_dx_f104_1, lobatto_dx_f104_2, lobatto_dx_f104_3, lobatto_dx_f104_4, lobatto_dx_f104_5, lobatto_dx_f105_0,
	lobatto_dx_f105_1, lobatto_dx_f105_2, lobatto_dx_f105_3, lobatto_dx_f105_4, lobatto_dx_f105_5, lobatto_dx_f106_0, lobatto_dx_f106_1,
	lobatto_dx_f106_2, lobatto_dx_f106_3, lobatto_dx_f106_4, lobatto_dx_f106_5, lobatto_dx_f107_0, lobatto_dx_f107_1, lobatto_dx_f107_2,
	lobatto_dx_f107_3, lobatto_dx_f107_4, lobatto_dx_f107_5, lobatto_dx_f108_0, lobatto_dx_f108_1, lobatto_dx_f108_2, lobatto_dx_f108_3,
	lobatto_dx_f108_4, lobatto_dx_f108_5, lobatto_dx_f109_0, lobatto_dx_f109_1, lobatto_dx_f109_2, lobatto_dx_f109_3, lobatto_dx_f109_4,
	lobatto_dx_f109_5, lobatto_dx_f110, lobatto_dx_f111, lobatto_dx_f112, lobatto_dx_f113, lobatto_dx_f114, lobatto_dx_f115,
	lobatto_dx_f116, lobatto_dx_f117, lobatto_dx_f118, lobatto_dx_f119, lobatto_dx_f120, lobatto_dx_f121, lobatto_dx_f122,
	lobatto_dx_f123, lobatto_dx_f124, lobatto_dx_f125, lobatto_dx_f126_0, lobatto_dx_f126_1, lobatto_dx_f126_2, lobatto_dx_f126_3,
	lobatto_dx_f126_4, lobatto_dx_f126_5, lobatto_dx_f127_0, lobatto_dx_f127_1, lobatto_dx_f127_2, lobatto_dx_f127_3, lobatto_dx_f127_4,
	lobatto_dx_f127_5, lobatto_dx_f128_0, lobatto_dx_f128_1, lobatto_dx_f128_2, lobatto_dx_f128_3, lobatto_dx_f128_4, lobatto_dx_f128_5,
	lobatto_dx_f129_0, lobatto_dx_f129_1, lobatto_dx_f129_2, lobatto_dx_f129_3, lobatto_dx_f129_4, lobatto_dx_f129_5, lobatto_dx_f130_0,
	lobatto_dx_f130_1, lobatto_dx_f130_2, lobatto_dx_f130_3, lobatto_dx_f130_4, lobatto_dx_f130_5, lobatto_dx_f131_0, lobatto_dx_f131_1,
	lobatto_dx_f131_2, lobatto_dx_f131_3, lobatto_dx_f131_4, lobatto_dx_f131_5, lobatto_dx_f132_0, lobatto_dx_f132_1, lobatto_dx_f132_2,
	lobatto_dx_f132_3, lobatto_dx_f132_4, lobatto_dx_f132_5, lobatto_dx_f133_0, lobatto_dx_f133_1, lobatto_dx_f133_2, lobatto_dx_f133_3,
	lobatto_dx_f133_4, lobatto_dx_f133_5, lobatto_dx_f134_0, lobatto_dx_f134_1, lobatto_dx_f134_2, lobatto_dx_f134_3, lobatto_dx_f134_4,
	lobatto_dx_f134_5, lobatto_dx_f135_0, lobatto_dx_f135_1, lobatto_dx_f135_2, lobatto_dx_f135_3, lobatto_dx_f135_4, lobatto_dx_f135_5,
	lobatto_dx_f136_0, lobatto_dx_f136_1, lobatto_dx_f136_2, lobatto_dx_f136_3, lobatto_dx_f136_4, lobatto_dx_f136_5, lobatto_dx_f137_0,
	lobatto_dx_f137_1, lobatto_dx_f137_2, lobatto_dx_f137_3, lobatto_dx_f137_4, lobatto_dx_f137_5, lobatto_dx_f138_0, lobatto_dx_f138_1,
	lobatto_dx_f138_2, lobatto_dx_f138_3, lobatto_dx_f138_4, lobatto_dx_f138_5, lobatto_dx_f139_0, lobatto_dx_f139_1, lobatto_dx_f139_2,
	lobatto_dx_f139_3, lobatto_dx_f139_4, lobatto_dx_f139_5, lobatto_dx_f140_0, lobatto_dx_f140_1, lobatto_dx_f140_2, lobatto_dx_f140_3,
	lobatto_dx_f140_4, lobatto_dx_f140_5, lobatto_dx_f141_0, lobatto_dx_f141_1, lobatto_dx_f141_2, lobatto_dx_f141_3, lobatto_dx_f141_4,
	lobatto_dx_f141_5, lobatto_dx_f142_0, lobatto_dx_f142_1, lobatto_dx_f142_2, lobatto_dx_f142_3, lobatto_dx_f142_4, lobatto_dx_f142_5,
	lobatto_dx_f143_0, lobatto_dx_f143_1, lobatto_dx_f143_2, lobatto_dx_f143_3, lobatto_dx_f143_4, lobatto_dx_f143_5, lobatto_dx_f144_0,
	lobatto_dx_f144_1, lobatto_dx_f144_2, lobatto_dx_f144_3, lobatto_dx_f144_4, lobatto_dx_f144_5, lobatto_dx_f145_0, lobatto_dx_f145_1,
	lobatto_dx_f145_2, lobatto_dx_f145_3, lobatto_dx_f145_4, lobatto_dx_f145_5, lobatto_dx_f146_0, lobatto_dx_f146_1, lobatto_dx_f146_2,
	lobatto_dx_f146_3, lobatto_dx_f146_4, lobatto_dx_f146_5, lobatto_dx_f147_0, lobatto_dx_f147_1, lobatto_dx_f147_2, lobatto_dx_f147_3,
	lobatto_dx_f147_4, lobatto_dx_f147_5, lobatto_dx_f148_0, lobatto_dx_f148_1, lobatto_dx_f148_2, lobatto_dx_f148_3, lobatto_dx_f148_4,
	lobatto_dx_f148_5, lobatto_dx_f149_0, lobatto_dx_f149_1, lobatto_dx_f149_2, lobatto_dx_f149_3, lobatto_dx_f149_4, lobatto_dx_f149_5,
	lobatto_dx_f150, lobatto_dx_f151, lobatto_dx_f152, lobatto_dx_f153, lobatto_dx_f154, lobatto_dx_f155, lobatto_dx_f156,
	lobatto_dx_f157, lobatto_dx_f158, lobatto_dx_f159, lobatto_dx_f160, lobatto_dx_f161, lobatto_dx_f162, lobatto_dx_f163,
	lobatto_dx_f164, lobatto_dx_f165_0, lobatto_dx_f165_1, lobatto_dx_f166_0, lobatto_dx_f166_1, lobatto_dx_f167_0, lobatto_dx_f167_1,
	lobatto_dx_f168_0, lobatto_dx_f168_1, lobatto_dx_f169_0, lobatto_dx_f169_1, lobatto_dx_f170_0, lobatto_dx_f170_1, lobatto_dx_f171_0,
	lobatto_dx_f171_1, lobatto_dx_f171_2, lobatto_dx_f171_3, lobatto_dx_f171_4, lobatto_dx_f171_5, lobatto_dx_f172_0, lobatto_dx_f172_1,
	lobatto_dx_f172_2, lobatto_dx_f172_3, lobatto_dx_f172_4, lobatto_dx_f172_5, lobatto_dx_f173_0, lobatto_dx_f173_1, lobatto_dx_f173_2,
	lobatto_dx_f173_3, lobatto_dx_f173_4, lobatto_dx_f173_5, lobatto_dx_f174_0, lobatto_dx_f174_1, lobatto_dx_f174_2, lobatto_dx_f174_3,
	lobatto_dx_f174_4, lobatto_dx_f174_5, lobatto_dx_f175_0, lobatto_dx_f175_1, lobatto_dx_f175_2, lobatto_dx_f175_3, lobatto_dx_f175_4,
	lobatto_dx_f175_5, lobatto_dx_f176_0, lobatto_dx_f176_1, lobatto_dx_f176_2, lobatto_dx_f176_3, lobatto_dx_f176_4, lobatto_dx_f176_5,
	lobatto_dx_f177_0, lobatto_dx_f177_1, lobatto_dx_f177_2, lobatto_dx_f177_3, lobatto_dx_f177_4, lobatto_dx_f177_5, lobatto_dx_f178_0,
	lobatto_dx_f178_1, lobatto_dx_f178_2, lobatto_dx_f178_3, lobatto_dx_f178_4, lobatto_dx_f178_5, lobatto_dx_f179_0, lobatto_dx_f179_1,
	lobatto_dx_f179_2, lobatto_dx_f179_3, lobatto_dx_f179_4, lobatto_dx_f179_5, lobatto_dx_f180_0, lobatto_dx_f180_1, lobatto_dx_f180_2,
	lobatto_dx_f180_3, lobatto_dx_f180_4, lobatto_dx_f180_5, lobatto_dx_f181_0, lobatto_dx_f181_1, lobatto_dx_f181_2, lobatto_dx_f181_3,
	lobatto_dx_f181_4, lobatto_dx_f181_5, lobatto_dx_f182_0, lobatto_dx_f182_1, lobatto_dx_f182_2, lobatto_dx_f182_3, lobatto_dx_f182_4,
	lobatto_dx_f182_5, lobatto_dx_f183_0, lobatto_dx_f183_1, lobatto_dx_f183_2, lobatto_dx_f183_3, lobatto_dx_f183_4, lobatto_dx_f183_5,
	lobatto_dx_f184_0, lobatto_dx_f184_1, lobatto_dx_f184_2, lobatto_dx_f184_3, lobatto_dx_f184_4, lobatto_dx_f184_5, lobatto_dx_f185_0,
	lobatto_dx_f185_1, lobatto_dx_f185_2, lobatto_dx_f185_3, lobatto_dx_f185_4, lobatto_dx_f185_5, lobatto_dx_f186_0, lobatto_dx_f186_1,
	lobatto_dx_f186_2, lobatto_dx_f186_3, lobatto_dx_f186_4, lobatto_dx_f186_5, lobatto_dx_f187_0, lobatto_dx_f187_1, lobatto_dx_f187_2,
	lobatto_dx_f187_3, lobatto_dx_f187_4, lobatto_dx_f187_5, lobatto_dx_f188_0, lobatto_dx_f188_1, lobatto_dx_f188_2, lobatto_dx_f188_3,
	lobatto_dx_f188_4, lobatto_dx_f188_5, lobatto_dx_f189_0, lobatto_dx_f189_1, lobatto_dx_f189_2, lobatto_dx_f189_3, lobatto_dx_f189_4,
	lobatto_dx_f189_5, lobatto_dx_f190_0, lobatto_dx_f190_1, lobatto_dx_f190_2, lobatto_dx_f190_3, lobatto_dx_f190_4, lobatto_dx_f190_5,
	lobatto_dx_f191_0, lobatto_dx_f191_1, lobatto_dx_f191_2, lobatto_dx_f191_3, lobatto_dx_f191_4, lobatto_dx_f191_5, lobatto_dx_f192_0,
	lobatto_dx_f192_1, lobatto_dx_f192_2, lobatto_dx_f192_3, lobatto_dx_f192_4, lobatto_dx_f192_5, lobatto_dx_f193_0, lobatto_dx_f193_1,
	lobatto_dx_f193_2, lobatto_dx_f193_3, lobatto_dx_f193_4, lobatto_dx_f193_5, lobatto_dx_f194_0, lobatto_dx_f194_1, lobatto_dx_f194_2,
	lobatto_dx_f194_3, lobatto_dx_f194_4, lobatto_dx_f194_5, lobatto_dx_f195_0, lobatto_dx_f195_1, lobatto_dx_f195_2, lobatto_dx_f195_3,
	lobatto_dx_f195_4, lobatto_dx_f195_5, lobatto_dx_f196_0, lobatto_dx_f196_1, lobatto_dx_f196_2, lobatto_dx_f196_3, lobatto_dx_f196_4,
	lobatto_dx_f196_5, lobatto_dx_f197_0, lobatto_dx_f197_1, lobatto_dx_f197_2, lobatto_dx_f197_3, lobatto_dx_f197_4, lobatto_dx_f197_5,
	lobatto_dx_f198_0, lobatto_dx_f198_1, lobatto_dx_f198_2, lobatto_dx_f198_3, lobatto_dx_f198_4, lobatto_dx_f198_5, lobatto_dx_f199,
	lobatto_dx_f200, lobatto_dx_f201, lobatto_dx_f202, lobatto_dx_f203, lobatto_dx_f204, lobatto_dx_f205, lobatto_dx_f206,
	lobatto_dx_f207, lobatto_dx_f208, lobatto_dx_f209, lobatto_dx_f210, lobatto_dx_f211, lobatto_dx_f212, lobatto_dx_f213,
	lobatto_dx_f214, lobatto_dx_f215, lobatto_dx_f216, lobatto_dx_f217, lobatto_dx_f218, lobatto_dx_f219, lobatto_dx_f220,
	lobatto_dx_f221, lobatto_dx_f222, lobatto_dx_f223, lobatto_dx_f224, lobatto_dx_f225, lobatto_dx_f226_0, lobatto_dx_f226_1,
	lobatto_dx_f226_2, lobatto_dx_f226_3, lobatto_dx_f226_4, lobatto_dx_f226_5, lobatto_dx_f227_0, lobatto_dx_f227_1, lobatto_dx_f227_2,
	lobatto_dx_f227_3, lobatto_dx_f227_4, lobatto_dx_f227_5, lobatto_dx_f228_0, lobatto_dx_f228_1, lobatto_dx_f228_2, lobatto_dx_f228_3,
	lobatto_dx_f228_4, lobatto_dx_f228_5, lobatto_dx_f229_0, lobatto_dx_f229_1, lobatto_dx_f229_2, lobatto_dx_f229_3, lobatto_dx_f229_4,
	lobatto_dx_f229_5, lobatto_dx_f230_0, lobatto_dx_f230_1, lobatto_dx_f230_2, lobatto_dx_f230_3, lobatto_dx_f230_4, lobatto_dx_f230_5,
	lobatto_dx_f231_0, lobatto_dx_f231_1, lobatto_dx_f231_2, lobatto_dx_f231_3, lobatto_dx_f231_4, lobatto_dx_f231_5, lobatto_dx_f232_0,
	lobatto_dx_f232_1, lobatto_dx_f232_2, lobatto_dx_f232_3, lobatto_dx_f232_4, lobatto_dx_f232_5, lobatto_dx_f233_0, lobatto_dx_f233_1,
	lobatto_dx_f233_2, lobatto_dx_f233_3, lobatto_dx_f233_4, lobatto_dx_f233_5, lobatto_dx_f234_0, lobatto_dx_f234_1, lobatto_dx_f234_2,
	lobatto_dx_f234_3, lobatto_dx_f234_4, lobatto_dx_f234_5, lobatto_dx_f235_0, lobatto_dx_f235_1, lobatto_dx_f235_2, lobatto_dx_f235_3,
	lobatto_dx_f235_4, lobatto_dx_f235_5, lobatto_dx_f236_0, lobatto_dx_f236_1, lobatto_dx_f236_2, lobatto_dx_f236_3, lobatto_dx_f236_4,
	lobatto_dx_f236_5, lobatto_dx_f237_0, lobatto_dx_f237_1, lobatto_dx_f237_2, lobatto_dx_f237_3, lobatto_dx_f237_4, lobatto_dx_f237_5,
	lobatto_dx_f238_0, lobatto_dx_f238_1, lobatto_dx_f238_2, lobatto_dx_f238_3, lobatto_dx_f238_4, lobatto_dx_f238_5, lobatto_dx_f239_0,
	lobatto_dx_f239_1, lobatto_dx_f239_2, lobatto_dx_f239_3, lobatto_dx_f239_4, lobatto_dx_f239_5, lobatto_dx_f240_0, lobatto_dx_f240_1,
	lobatto_dx_f240_2, lobatto_dx_f240_3, lobatto_dx_f240_4, lobatto_dx_f240_5, lobatto_dx_f241_0, lobatto_dx_f241_1, lobatto_dx_f241_2,
	lobatto_dx_f241_3, lobatto_dx_f241_4, lobatto_dx_f241_5, lobatto_dx_f242_0, lobatto_dx_f242_1, lobatto_dx_f242_2, lobatto_dx_f242_3,
	lobatto_dx_f242_4, lobatto_dx_f242_5, lobatto_dx_f243_0, lobatto_dx_f243_1, lobatto_dx_f243_2, lobatto_dx_f243_3, lobatto_dx_f243_4,
	lobatto_dx_f243_5, lobatto_dx_f244_0, lobatto_dx_f244_1, lobatto_dx_f244_2, lobatto_dx_f244_3, lobatto_dx_f244_4, lobatto_dx_f244_5,
	lobatto_dx_f245_0, lobatto_dx_f245_1, lobatto_dx_f245_2, lobatto_dx_f245_3, lobatto_dx_f245_4, lobatto_dx_f245_5, lobatto_dx_f246_0,
	lobatto_dx_f246_1, lobatto_dx_f246_2, lobatto_dx_f246_3, lobatto_dx_f246_4, lobatto_dx_f246_5, lobatto_dx_f247_0, lobatto_dx_f247_1,
	lobatto_dx_f247_2, lobatto_dx_f247_3, lobatto_dx_f247_4, lobatto_dx_f247_5, lobatto_dx_f248_0, lobatto_dx_f248_1, lobatto_dx_f248_2,
	lobatto_dx_f248_3, lobatto_dx_f248_4, lobatto_dx_f248_5, lobatto_dx_f249_0, lobatto_dx_f249_1, lobatto_dx_f249_2, lobatto_dx_f249_3,
	lobatto_dx_f249_4, lobatto_dx_f249_5, lobatto_dx_f250_0, lobatto_dx_f250_1, lobatto_dx_f250_2, lobatto_dx_f250_3, lobatto_dx_f250_4,
	lobatto_dx_f250_5, lobatto_dx_f251_0, lobatto_dx_f251_1, lobatto_dx_f251_2, lobatto_dx_f251_3, lobatto_dx_f251_4, lobatto_dx_f251_5,
	lobatto_dx_f252_0, lobatto_dx_f252_1, lobatto_dx_f252_2, lobatto_dx_f252_3, lobatto_dx_f252_4, lobatto_dx_f252_5, lobatto_dx_f253_0,
	lobatto_dx_f253_1, lobatto_dx_f253_2, lobatto_dx_f253_3, lobatto_dx_f253_4, lobatto_dx_f253_5, lobatto_dx_f254_0, lobatto_dx_f254_1,
	lobatto_dx_f254_2, lobatto_dx_f254_3, lobatto_dx_f254_4, lobatto_dx_f254_5, lobatto_dx_f255_0, lobatto_dx_f255_1, lobatto_dx_f255_2,
	lobatto_dx_f255_3, lobatto_dx_f255_4, lobatto_dx_f255_5, lobatto_dx_f256_0, lobatto_dx_f256_1, lobatto_dx_f256_2, lobatto_dx_f256_3,
	lobatto_dx_f256_4, lobatto_dx_f256_5, lobatto_dx_f257_0, lobatto_dx_f257_1, lobatto_dx_f257_2, lobatto_dx_f257_3, lobatto_dx_f257_4,
	lobatto_dx_f257_5, lobatto_dx_f258, lobatto_dx_f259, lobatto_dx_f260, lobatto_dx_f261, lobatto_dx_f262, lobatto_dx_f263,
	lobatto_dx_f264, lobatto_dx_f265, lobatto_dx_f266, lobatto_dx_f267, lobatto_dx_f268, lobatto_dx_f269, lobatto_dx_f270,
	lobatto_dx_f271, lobatto_dx_f272, lobatto_dx_f273, lobatto_dx_f274, lobatto_dx_f275, lobatto_dx_f276, lobatto_dx_f277,
	lobatto_dx_f278, lobatto_dx_f279, lobatto_dx_f280, lobatto_dx_f281, lobatto_dx_f282, lobatto_dx_f283, lobatto_dx_f284,
	lobatto_dx_f285
};

#endif
