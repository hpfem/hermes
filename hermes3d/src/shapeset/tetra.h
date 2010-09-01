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

#ifndef _SHAPESET_TETRA_H_
#define _SHAPESET_TETRA_H_

// tetra specific macros

// validation macros
#define CHECK_VERTEX(vertex)    	assert(vertex >= 0 && vertex < 4)

#define CHECK_EDGE(edge)	      	assert(edge >= 0 && edge < 6)
#define CHECK_EDGE_ORI(ori)	      	assert(ori == 0 || ori == 1)

#define CHECK_FACE(face)			assert(face >= 0 && face < 4)
#define CHECK_FACE_MODE(mode)		assert(mode == MODE_TRIANGLE)
#define CHECK_FACE_ORI(ori)	      	assert(ori >= 0 || ori <= 6)

//#define CHECK_PART(p)				assert(p >= 0)

#endif /* SHAPESET_TETRA_H_ */
