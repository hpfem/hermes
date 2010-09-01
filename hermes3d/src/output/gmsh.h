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

#ifndef _GMSH_OUTPUT_ENGINE_H_
#define _GMSH_OUTPUT_ENGINE_H_

#include "../output.h"
#include "../matrix.h"

/// GMSH output engine.
///
///
///
/// @ingroup visualization
class GmshOutputEngine : public OutputEngine {
public:
	GmshOutputEngine(FILE *file);
	virtual ~GmshOutputEngine();

	/// Run the output with specified output engine
	///
	/// @return true if ok
	/// @param[in] fn A function that will be visualized
	virtual void out(MeshFunction *fn, const char *name, int item = FN_VAL);
	virtual void out(MeshFunction *fn1, MeshFunction *fn2, MeshFunction *fn3, const char *name, int item = FN_VAL_0);
	virtual void out(Mesh *mesh);
	virtual void out_bc(Mesh *mesh, const char *name = "BCs");

	virtual void out_orders(Space *space, const char *name = "orders");

	virtual void out(Matrix *mat);

protected:
	/// file into which the output is done
	FILE *out_file;

	void dump_scalars(int mode, int num_pts, Point3D *pts, double *value);
	void dump_vectors(int mode, int num_pts, Point3D *pts, double *v0, double *v1, double *v2);
	void dump_mesh(Mesh *mesh);
};

#endif
