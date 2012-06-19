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

#ifndef _MESHLOADER_H_
#define _MESHLOADER_H_

#include "mesh.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @defgroup mesh_readers Mesh readers

    /// Abstract class for mesh readers
    ///
    /// @ingroup mesh_readers
    class HERMES_API MeshReader : public Hermes::Mixins::Loggable
    {
    public:
      virtual ~MeshReader() { }

      /// Loads the mesh from a file. Aborts the program on error.
      /// @param filename [in] The name of the file.
      /// @param mesh [out] The mesh.
      virtual bool load(const char *filename, Mesh *mesg) = 0;

    protected:
      static bool is_twin_nurbs(Element* e, int i)
      {
        // on internal edges, where there are two Nurbs', we only save one of them
        return e->cm->nurbs[i]->twin && e->en[i]->ref == 2;
      }
    };
  }
}
#endif
