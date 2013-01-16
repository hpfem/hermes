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
// You should have received a copy of the GNU General Public Licenserix
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_ELEMWISE_PARAMETER_MESH_FUNC_H
#define __H2D_ELEMWISE_PARAMETER_MESH_FUNC_H

#include "elemwise_parameter_func.h"
#include "../function/mesh_function.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar> class MeshFunction;

    /// Mesh function - based parameter.
    template<typename Scalar>
    class HERMES_API ElemwiseParameterMeshFunc : public ElemwiseParameterFunc<Scalar>, public Hermes2D::Mixins::StateQueryable
    {
    public:
      /// Constructor.
      ElemwiseParameterMeshFunc(MeshFunction<Scalar>* mesh_function);
      ElemwiseParameterMeshFunc();

      /// Destructor.
      virtual ~ElemwiseParameterMeshFunc();

      /// Ask if the instance is fine.
      virtual bool isOkay() const;

        /// Get class name, for the purpose of messaging.
      virtual std::string getClassName() const;

    protected:
      virtual Scalar get_value(Element* e);

      template<typename T> friend class DiscreteProblemLinear;

    private:
      MeshFunction<Scalar>* mesh_function;
    };
  }
}
#endif
