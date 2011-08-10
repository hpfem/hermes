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

#ifndef _MESH_READER_H2D_XML_H_
#define _MESH_READER_H2D_XML_H_

#include "mesh_reader.h"
#ifdef error(...)
#undef error(...)
#endif
#include "mesh_h2d_xml.h"

// This is here mainly because XSD uses its own error, therefore it had to be undefined previously.
#ifndef error(...)
#define error(...) hermes_exit_if(hermes_log_message_if(true, HERMES_BUILD_LOG_INFO(HERMES_EC_ERROR), __VA_ARGS__))
#endif

namespace Hermes
{
  namespace Hermes2D
  {
    /// Mesh reader from Hermes2D format
    ///
    /// @ingroup mesh_readers
    class HERMES_API MeshReaderH2DXML : public MeshReader
    {
    public:
      MeshReaderH2DXML();
      virtual ~MeshReaderH2DXML();

      virtual bool load(const char *file_name, Mesh *mesh);

    protected:
      Nurbs* load_nurbs(Mesh *mesh, std::auto_ptr<mesh_h2d> & m, int id, Node** en, int &p1, int &p2);
    };
  }
}
#endif

