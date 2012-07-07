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

#ifndef _MESH_READER_H2D_H_
#define _MESH_READER_H2D_H_

#include "mesh_reader.h"
#include "mesh_data.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Mesh reader from Hermes2D format
    ///
    /// @ingroup mesh_readers
    class HERMES_API MeshReaderH2D : public MeshReader
    {
    public:
      MeshReaderH2D();
      virtual ~MeshReaderH2D();

      virtual bool load(const char *filename, Mesh *mesh);
      virtual bool save(const char *filename, Mesh *mesh);

    protected:
      Nurbs* load_nurbs(Mesh *mesh, MeshData *m, int id, Node** en, int &p1, int &p2);

      void save_refinements(Mesh *mesh, FILE* f, Element* e, int id, bool& first);
      void save_nurbs(Mesh *mesh, FILE* f, int p1, int p2, Nurbs* nurbs);
    };
  }
}
#endif