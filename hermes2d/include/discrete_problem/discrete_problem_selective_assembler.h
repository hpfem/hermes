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

#ifndef __H2D_DISCRETE_PROBLEM_SELECTIVE_ASSEMBLER_H
#define __H2D_DISCRETE_PROBLEM_SELECTIVE_ASSEMBLER_H

#include "hermes_common.h"
#include "forms.h"
#include "weakform/weakform.h"
#include "mesh/traverse.h"
#include "space/space.h"
#include "mixins2d.h"
#include "discrete_problem_helpers.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class PrecalcShapeset;
    /// @ingroup inner
    /// Discrete problem class.
    ///
    /// This class does assembling into external matrix / vector structures.
    ///
    template<typename Scalar>
    class HERMES_API DiscreteProblemSelectiveAssembler : 
      public Hermes::Hermes2D::Mixins::DiscreteProblemRungeKutta<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>
    {
    public:
      DiscreteProblemSelectiveAssembler();
      ~DiscreteProblemSelectiveAssembler();
      
      /// Preassembling.
      /// Precalculate matrix sparse structure.
      /// If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems.
      /// Returns false if there are no states to assemble.
      bool prepare_sparse_structure(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State**& states, int& num_states);
      
      /// Sets new spaces for the instance.
      void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      void alloc_recalculation_tables_spaces_settings(typename WeakForm<Scalar>::FormIntegrationDimension dimension, typename WeakForm<Scalar>::FormEquationSide equation_side, int new_markers_count);
      
      /// Set the weak forms.
      void set_weak_formulation(WeakForm<Scalar>* wf);
      void alloc_recalculation_tables_weakform_settings(typename WeakForm<Scalar>::FormIntegrationDimension dimension, typename WeakForm<Scalar>::FormEquationSide equation_side);

      /// The form will be assembled.
      bool form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormDG<Scalar>* form, Traverse::State* current_state);

      bool form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormDG<Scalar>* form, Traverse::State* current_state);

      /// Recalculation storages.
      bool* state_reuse_kept[2][2];
      int markers_size[2][2];

    protected:
      /// Spaces.
      int spaces_size;

      /// Seq numbers of Space instances in spaces.
      int* sp_seq;

      /// Matrix structure can be reused.
      /// If other conditions apply.
      bool matrix_structure_reusable;
      bool vector_structure_reusable;

      friend class DiscreteProblem<Scalar>;
    };
  }
}
#endif
