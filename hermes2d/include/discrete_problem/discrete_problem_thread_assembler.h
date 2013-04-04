/// This file is part of Hermes2D.
///
/// Hermes2D is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 2 of the License, or
/// (at your option) any later version.
///
/// Hermes2D is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY;without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Hermes2D. If not, see <http:///www.gnu.org/licenses/>.

#ifndef __H2D_DISCRETE_PROBLEM_ASSEMBLY_DATA_H
#define __H2D_DISCRETE_PROBLEM_ASSEMBLY_DATA_H

#include "../weakform/weakform.h"
#include "../shapeset/precalc.h"
#include "../function/solution.h"
#include "discrete_problem_helpers.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup inner
    /// DiscreteProblemAssemblyData class.
    ///
    /// This class does assembling into external matrix / vector structures.
    ///
    template<typename Scalar>
    class HERMES_API DiscreteProblemThreadAssembler : 
      public Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemRungeKutta<Scalar>
    {
    private:
      DiscreteProblemThreadAssembler();
      ~DiscreteProblemThreadAssembler();

      void init_assembling();
      void init_assembling_one_state(Traverse::State* current_state);

      void init_spaces(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      void set_weak_formulation(WeakForm<Scalar>* wf);
      void init_u_ext(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Solution<Scalar>** u_ext_sln);

      void free();
      void free_spaces();
      void free_weak_formulation();
      void free_u_ext();

      PrecalcShapeset** pss;
      RefMap** refmaps;
      Solution<Scalar>** u_ext;
      AsmList<Scalar>** als;
      AsmList<Scalar>*** alsSurface;

      Hermes::vector<Transformable *> fns;  

      int spaces_size;
      friend class DiscreteProblem<Scalar>;
    };
  }
}
#endif
