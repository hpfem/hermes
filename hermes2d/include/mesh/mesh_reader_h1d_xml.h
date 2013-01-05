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

#ifndef _MESH_READER_H1D_XML_H_
#define _MESH_READER_H1D_XML_H_

#include "mesh_reader.h"

// This is here mainly because XSD uses its own error, therefore it has to be undefined here.
#ifdef error
#undef error
#endif

#include "mesh_h1d_xml.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Mesh reader from Hermes2D format
    ///
    /// @ingroup mesh_readers
    class HERMES_API MeshReaderH1DXML : public MeshReader, public Hermes::Hermes2D::Mixins::XMLParsing
    {
    public:
      MeshReaderH1DXML();
      virtual ~MeshReaderH1DXML();

      /// This method loads a single mesh from a file.
      virtual bool load(const char *filename, Mesh *mesh);

      /// This method saves a single mesh to a file.
      bool save(const char *filename, Mesh *mesh);
    };
  }
}
#endif
