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

#ifndef __H2D_DISCRETE_PROBLEM_LINEAR_H
#define __H2D_DISCRETE_PROBLEM_LINEAR_H

#include "discrete_problem.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// SIMPLIFICATION of Discrete problem class for LINEAR PROBLEMS.
    ///
    /// This class does assembling into external matrix / vector structures.
    ///
    template<typename Scalar>
    class HERMES_API DiscreteProblemLinear : public DiscreteProblem<Scalar>
    {
    public:
      /// Constructor for multiple components / equations.
      DiscreteProblemLinear(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar> *> spaces);

      /// Constructor for one equation.
      DiscreteProblemLinear(const WeakForm<Scalar>* wf, const Space<Scalar>* space);

      /// Destuctor.
      virtual ~DiscreteProblemLinear();

      /// Assembling.
      /// Light version, linear problems.
      virtual void assemble(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL, bool force_diagonal_blocks = false,
        Table* block_weights = NULL);

    protected:
      /// Methods different to those of the parent class.
      /// Matrix forms.
      virtual void assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state);

      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class NewtonSolver;
      template<typename T> friend class PicardSolver;
      template<typename T> friend class RungeKutta;
    };
  }
}
#endif
