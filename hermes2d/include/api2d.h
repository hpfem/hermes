// This file is part of Hermes2D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file api.h
\brief Main Hermes API
*/
#ifndef __HERMES_API_2D_H_
#define __HERMES_API_2D_H_

#include "compat.h"
#include "hermes_common.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// API Class containing settings for the whole Hermes.
    class HERMES_API Api2D : public Hermes::Api
    {
    public:
      Api2D();
      virtual void init();
    };

    // Global declarations.
    extern HERMES_API Hermes::Hermes2D::Api2D Hermes2DApi;
  }
}
#endif