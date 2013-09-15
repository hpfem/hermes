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

#ifndef __H2D_SHAPESET_HC_ALL_H
#define __H2D_SHAPESET_HC_ALL_H

#include "shapeset.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// H(curl) shapeset based on Legendre polynomials.
    /// @ingroup spaces
    class HERMES_API HcurlShapesetLegendre : public Shapeset
    {
    public:
      HcurlShapesetLegendre();
      virtual Shapeset* clone() { return new HcurlShapesetLegendre(*this); };
      virtual SpaceType get_space_type() const { return HERMES_HCURL_SPACE; }
      virtual int get_max_index(ElementMode2D mode);
      virtual int get_id() const { return 10; }
      
      static const int max_index[H2D_NUM_MODES];
    };

    /*
    // Experimental.
    /// @ingroup spaces
    class HERMES_API HcurlShapesetEigen2 : public Shapeset
    {
    public:
      HcurlShapesetEigen2();
      virtual Shapeset* clone() { return new HcurlShapesetEigen2(*this); };
      virtual int get_max_index(ElementMode2D mode);
    protected:
      template<typename Scalar> friend class VectorForm;
      template<typename Scalar> friend class MatrixForm;
      virtual int get_id() const { return 11; }
      virtual SpaceType get_space_type() const { return HERMES_HCURL_SPACE; }
      template<typename Scalar> friend class DiscreteProblem; template<typename Scalar> friend class Solution; friend class CurvMap; friend class RefMap; template<typename Scalar> friend class RefinementSelectors::H1ProjBasedSelector; template<typename Scalar> friend class RefinementSelectors::L2ProjBasedSelector; template<typename Scalar> friend class RefinementSelectors::HcurlProjBasedSelector; template<typename Scalar> friend class RefinementSelectors::OptimumSelector; friend class PrecalcShapeset;
      static const int max_index[2];
    };
    */

    /*
    /// Experimental.
    /// @ingroup spaces
    class HERMES_API HcurlShapesetGradEigen : public Shapeset
    {
    public:
      HcurlShapesetGradEigen();
      virtual Shapeset* clone() { return new HcurlShapesetGradEigen(*this); };
      virtual int get_max_index(ElementMode2D mode);
    protected:
      template<typename Scalar> friend class VectorForm;
      template<typename Scalar> friend class MatrixForm;
      virtual int get_id() const { return 12; }
      virtual SpaceType get_space_type() const { return HERMES_HCURL_SPACE; }
      template<typename Scalar> friend class DiscreteProblem; template<typename Scalar> friend class Solution; friend class CurvMap; friend class RefMap; template<typename Scalar> friend class RefinementSelectors::H1ProjBasedSelector; template<typename Scalar> friend class RefinementSelectors::L2ProjBasedSelector; template<typename Scalar> friend class RefinementSelectors::HcurlProjBasedSelector; template<typename Scalar> friend class RefinementSelectors::OptimumSelector; friend class PrecalcShapeset;
      static const int max_index[2];
    };
    */

    /// H(curl) shapeset with Legendre bubbles and gradients of H1 functions as edges
    /// @ingroup spaces
    class HERMES_API HcurlShapesetGradLeg : public Shapeset
    {
    public:
      HcurlShapesetGradLeg();
      virtual Shapeset* clone() { return new HcurlShapesetGradLeg(*this); };
      virtual int get_max_index(ElementMode2D mode);
      virtual int get_id() const { return 13; }
      virtual SpaceType get_space_type() const { return HERMES_HCURL_SPACE; }
      
      static const int max_index[H2D_NUM_MODES];
    };

    /// This is the default Hcurl shapeset typedef.
    typedef HcurlShapesetGradLeg HcurlShapeset;
  }
}
#endif