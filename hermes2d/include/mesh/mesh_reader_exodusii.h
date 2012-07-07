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

#ifndef _EXODUSII_LOADER_H_
#define _EXODUSII_LOADER_H_

#include "mesh_reader.h"
namespace Hermes
{
  namespace Hermes2D
  {
    /// Mesh loader from EXODUSII format
    ///
    /// @ingroup mesh_readers
    class HERMES_API MeshReaderExodusII : public MeshReader
    {
    public:
      MeshReaderExodusII();
      virtual ~MeshReaderExodusII();

      virtual bool load(const char *file_name, Mesh *mesh);
    };
  }
}
#endif