// This file is part of Hermes2D
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, see <http://www.gnu.prg/licenses/>.

#ifndef _H2D_READER_H_
#define _H2D_READER_H_

#include "mesh_loader.h"

/// Mesh loader from Hermes2D format
///
/// @ingroup meshloaders
class H2D_API H2DReader : public MeshLoader
{
public:
  H2DReader();
  virtual ~H2DReader();

  virtual bool load(const char *file_name, Mesh *mesh);
  virtual bool save(const char *file_name, Mesh *mesh);

  void load_old(const char* filename, Mesh *mesh);
  void load_str(char* mesh_str, Mesh *mesh);
  void load_stream(FILE *f, Mesh *mesh);

protected:
	Nurbs* load_nurbs_old(Mesh *mesh, FILE* f, Node** en, int &p1, int &p2);
  Nurbs* load_nurbs(Mesh *mesh, MItem* curve, int id, Node** en, int &p1, int &p2);

	void save_refinements(Mesh *mesh, FILE* f, Element* e, int id, bool& first);
	void save_nurbs(Mesh *mesh, FILE* f, int p1, int p2, Nurbs* nurbs);
};

#endif

