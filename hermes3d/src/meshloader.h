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

#ifndef _MESHLOADER_H_
#define _MESHLOADER_H_

#include "mesh.h"

/// @defgroup meshloaders Mesh loaders

/// Abstract class for mesh loaders
///
/// @ingroup meshloaders
class MeshLoader {
public:
	virtual ~MeshLoader() { }

	/// Loads the mesh from a file. Aborts the program on error.
	/// @param filename [in] The name of the file.
	/// @param mesh [out] The mesh.
	virtual bool load(const char *file_name, Mesh *mesg) = 0;

	/// Saves the mesh, including all refinements, to a file.
	/// Caution: never overwrite hand-created meshes with this function --
	/// all comments in the original file will be lost.
	/// @param filename [in] The name of the file.
	/// @param mesh [out] The mesh.
	virtual bool save(const char *file_name, Mesh *mesh) = 0;
};

#endif
