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
             virtual unsigned char get_id()
             {
               return 1;
             };

             virtual void dummy_fn() {}
    };

    extern HERMES_API Quad1DStd g_quad_1d_std;
    extern HERMES_API Quad2DStd g_quad_2d_std;

    //// linearization "quadrature" ////////////////////////////////////////////////////////////////////

    /// The tables with index zero are for obtaining solution values at the element
    /// vertices. Index one tables serve for the retrieval of interior values. Index one tables
    /// are used for adaptive approximation of the solution by transforming their points to sub-elements.
    /// Actually, the tables contain two levels of refinement -- this is an optimization to reduce
    /// the number of calls to sln->get_values().
    extern double3 lin_pts_0_tri[];

    extern double3 lin_pts_0_quad[];

    extern double3 lin_pts_1_tri[12];

    extern double3 lin_pts_1_quad[21];

    extern unsigned short quad_indices[9][5];

    extern unsigned short tri_indices[5][3];

    extern unsigned char lin_np_tri[2];
    extern unsigned char lin_np_quad[2];
    extern unsigned char* lin_np[2];

    extern double3*  lin_tables_tri[2];
    extern double3*  lin_tables_quad[2];
    extern double3** lin_tables[2];

    class Quad2DLin : public Quad2D
    {
    public:
      Quad2DLin();
      virtual unsigned char get_id()
      {
        return 2;
      };
    };

    extern HERMES_API Quad2DLin g_quad_lin;
  }
}
#endif