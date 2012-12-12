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

#ifndef __H2D_SHAPESET_L2_ALL
#define __H2D_SHAPESET_L2_ALL

// This file is a common header for all L2 shapesets.

#include "shapeset.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// L2 shapeset - products of legendre polynomials
    /// @ingroup spaces
    class HERMES_API L2ShapesetLegendre : public Shapeset
    {
    public:
      L2ShapesetLegendre();
      virtual Shapeset* clone() { return new L2ShapesetLegendre(*this); };
    protected:
      virtual int get_id() const { return 30; }
      virtual SpaceType get_space_type() const { return HERMES_L2_SPACE; }
      template<typename Scalar> friend class DiscreteProblem;
      template<typename Scalar> friend class Solution;
      friend class CurvMap; friend class RefMap;
      template<typename Scalar> friend class RefinementSelectors::H1ProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::L2ProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::HcurlProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::OptimumSelector; friend class PrecalcShapeset;
    };

    /// This is the default shapeset typedef
    typedef L2ShapesetLegendre L2Shapeset;
  }
}
#endif