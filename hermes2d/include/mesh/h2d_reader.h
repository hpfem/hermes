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
#include "mesh_data.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Mesh loader from Hermes2D format
    ///
    /// @ingroup meshloaders
    /// \todo Change the name to H2DMeshLoader
    class HERMES_API H2DReader : public MeshLoader
    {
    public:
      H2DReader();
      virtual ~H2DReader();

      virtual bool load(const char *file_name, Mesh *mesh);
      virtual bool save(const char *file_name, Mesh *mesh);

      void load_str(const char* mesh_str, Mesh *mesh);
      bool load_stream(std::istream &is, Mesh *mesh, const char *filename);

    protected:
      Nurbs* load_nurbs_old(Mesh *mesh, FILE* f, Node** en, int &p1, int &p2);
      /// Load NURBS.
      Nurbs* load_nurbs(Mesh *mesh, MeshData *m, int id, Node** en, int &p1, int &p2);

      void save_refinements(Mesh *mesh, FILE* f, Element* e, int id, bool& first);
      void save_nurbs(Mesh *mesh, FILE* f, int p1, int p2, Nurbs* nurbs);
    };
  }
}
#endif

