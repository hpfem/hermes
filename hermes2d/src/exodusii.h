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

#include "mesh_loader.h"

/// Mesh loader from EXODUSII format
///
/// @ingroup meshloaders
class H2D_API ExodusIIReader : public MeshLoader
{
public:
  ExodusIIReader();
  virtual ~ExodusIIReader();

  virtual bool load(const char *file_name, Mesh *mesh);
};

#endif
