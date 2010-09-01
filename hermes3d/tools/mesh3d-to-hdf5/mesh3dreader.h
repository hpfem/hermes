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

#ifndef MESH3D_READER_H_
#define MESH3D_READER_H_

class Mesh;

/// Mesh loader from Mesh3D format
///
/// @ingroup meshloaders
class Mesh3DLoader {
public:
	Mesh3DLoader();
	virtual ~Mesh3DLoader();

	virtual bool load(const char *file_name, Mesh *mesh);
};

#endif /*MESH3D_READER_H_*/
