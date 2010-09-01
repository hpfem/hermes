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

#ifndef _SHAPESET_COMMON_H_
#define _SHAPESET_COMMON_H_

#include "../common.h"

struct QuadPt3D;

// shape function
typedef double (*shape_fn_t)(double, double, double);

// shape function that is calculated on-the-fly
typedef void (*shape_fn_deleg_t)(int, int, QuadPt3D *, int, double *);

// shape function in 1D (used by product-based shapesets)
typedef double (*shape_fn_1d_t)(double);


/// Affine coordinates

#define lambda0(x,y,z) (((y) + 1) / 2)
#define lambda1(x,y,z) (-(1 + (x) + (y) + (z)) / 2)
#define lambda2(x,y,z) (((x) + 1) / 2)
#define lambda3(x,y,z) (((z) + 1) / 2)

/// X derivatives of affine coordinates

#define lambda0dx(x,y,z) (0.0)
#define lambda1dx(x,y,z) (-1.0 / 2.0)
#define lambda2dx(x,y,z) (1.0 / 2.0)
#define lambda3dx(x,y,z) (0.0)

/// Y derivatives of affine coordinates

#define lambda0dy(x,y,z) (1.0 / 2.0)
#define lambda1dy(x,y,z) (-1.0 / 2.0)
#define lambda2dy(x,y,z) (0.0)
#define lambda3dy(x,y,z) (0.0)

/// Z derivatives of affine coordinates

#define lambda0dz(x,y,z) (0.0)
#define lambda1dz(x,y,z) (-1.0 / 2.0)
#define lambda2dz(x,y,z) (0.0)
#define lambda3dz(x,y,z) (1.0 / 2.0)


// validation macros
#define CHECK_ORDER(o)  			assert((o) >= 0 && (o) <= max_order)
#define CHECK_INDEX(index)     		assert(index >= 0 && index <= max_index)
#define CHECK_COMPONENT(comp)	 	assert(comp >= 0 && comp < num_components)

#endif
