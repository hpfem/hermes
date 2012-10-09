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

#ifndef __H2D_LOCALPROJECTION_H
#define __H2D_LOCALPROJECTION_H

#include "../function/solution.h"
#include "../forms.h"
#include "../weakform/weakform.h"
#include "../views/scalar_view.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup projections
    template<typename Scalar>
    class HERMES_API LocalProjection
    {
    public:
      LocalProjection();

      // Main functionality.
      static void project_local(const Space<Scalar>* space, MeshFunction<Scalar>* meshfn,
          Scalar* target_vec, ProjNormType proj_norm = HERMES_UNSET_NORM);

      // Wrapper that delivers a Solution instead of coefficient vector.
      static void project_local(const Space<Scalar>* space,
    Solution<Scalar>* source_sln, Solution<Scalar>* target_sln, ProjNormType proj_norm = HERMES_UNSET_NORM);

      // Wrapper that takes multiple MeshFunctions.
      static void project_local(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<MeshFunction<Scalar>*> meshfns,
          Scalar* target_vec, Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());

      // Wrapper that takes multiple Solutions.
      static void project_local(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> slns,
          Scalar* target_vec, Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());

      // Wrapper that delivers Solutions instead of a coefficient vector.
      static void project_local(Hermes::vector<const Space<Scalar>*> spaces,
          Hermes::vector<Solution<Scalar>*> source_slns, Hermes::vector<Solution<Scalar>*> target_slns,
          Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>(), bool delete_old_mesh = false);

    protected:
      static int ndof;
    };
  }
}
#endif