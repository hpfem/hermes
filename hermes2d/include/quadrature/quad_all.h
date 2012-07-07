// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_QUAD_ALL_H
#define __H2D_QUAD_ALL_H

// This is a common header for all available 1D and 2D quadrature tables

#include "quad.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// 1D quadrature points on the standard reference domain (-1,1)
    class HERMES_API Quad1DStd : public Quad1D
    {
    public: Quad1DStd();

            virtual void dummy_fn() {}
    };

    /// 2D quadrature points on the standard reference domains (-1,1)^2
    class HERMES_API Quad2DStd : public Quad2D
    {
    public:  Quad2DStd();
             ~Quad2DStd();

             virtual void dummy_fn() {}
    };

    extern HERMES_API Quad1DStd g_quad_1d_std;
    extern HERMES_API Quad2DStd g_quad_2d_std;
  }
}
#endif