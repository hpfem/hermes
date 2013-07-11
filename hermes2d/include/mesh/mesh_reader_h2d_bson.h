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

#ifndef _MESH_READER_H2D_BSON_H_
#define _MESH_READER_H2D_BSON_H_

#include "mesh_reader.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Mesh reader from BSON format
    ///
    /// @ingroup mesh_readers
    /// Typical usage:
    /// MeshSharedPtr mesh;
    /// Hermes::Hermes2D::MeshReaderH2DBSON mloader;
    /// try
    /// {
    ///&nbsp;mloader.load("mesh.bson", &mesh);
    /// }
    /// catch(Exceptions::MeshLoadFailureException& e)
    /// {
    ///&nbsp;e.print_msg();
    ///&nbsp;return -1;
    /// }
    /// 
    class HERMES_API MeshReaderH2DBSON : public MeshReader
    {
    public:
      MeshReaderH2DBSON();
      virtual ~MeshReaderH2DBSON();

      /// This method loads a single mesh from a file.
      virtual bool load(const char *filename, MeshSharedPtr mesh);

      /// This method saves a single mesh to a file.
      bool save(const char *filename, MeshSharedPtr mesh);

      /// This method loads multiple meshes according to subdomains described in the meshfile.
      /// \param[in] meshes Meshes to be loaded, the number must correspond to the subdomains described in the file.
      ///&nbsp;         also the order is determined by the order in the file.
      bool load(const char *filename, Hermes::vector<MeshSharedPtr > meshes);

      /// This method saves multiple meshes according to subdomains in the vector meshes.
      bool save(const char *filename, Hermes::vector<MeshSharedPtr > meshes);
    };
  }
}
#endif
