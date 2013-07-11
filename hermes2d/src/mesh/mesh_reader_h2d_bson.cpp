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

#include "mesh.h"
#include "api2d.h"
#include "mesh_reader_h2d_bson.h"

using namespace std;

namespace Hermes
{
  namespace Hermes2D
  {
    MeshReaderH2DBSON::MeshReaderH2DBSON()
    {
    }

    MeshReaderH2DBSON::~MeshReaderH2DBSON()
    {
    }

    bool MeshReaderH2DBSON::load(const char *filename, MeshSharedPtr mesh)
    {
      throw Exceptions::MethodNotOverridenException("MeshReaderH2DBSON::load");

      return true;
    }

    bool MeshReaderH2DBSON::save(const char *filename, MeshSharedPtr mesh)
    {
      throw Exceptions::MethodNotOverridenException("MeshReaderH2DBSON::save");

      return true;
    }

    bool MeshReaderH2DBSON::load(const char *filename, Hermes::vector<MeshSharedPtr > meshes)
    {
      throw Exceptions::MethodNotOverridenException("MeshReaderH2DBSON::load");

      return true;
    }

    bool MeshReaderH2DBSON::save(const char *filename, Hermes::vector<MeshSharedPtr > meshes)
    {
      throw Exceptions::MethodNotOverridenException("MeshReaderH2DBSON::save");

      return true;
    }
  }
}
