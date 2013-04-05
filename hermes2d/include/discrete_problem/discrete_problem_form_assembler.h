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

#ifndef __H2D_DISCRETE_PROBLEM_FORM_ASSEMBLER_H
#define __H2D_DISCRETE_PROBLEM_FORM_ASSEMBLER_H

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
    class HERMES_API DiscreteProblemFormAssembler : 
      public Hermes::Hermes2D::Mixins::StateQueryable,
      public Hermes::Hermes2D::Mixins::DiscreteProblemRungeKutta<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemMatrixVector<Scalar>
    {
    public:

      /// State querying helpers.
      virtual bool isOkay() const { return true; }
      virtual inline std::string getClassName() const { return "DiscreteProblem"; }


      /// The form will be assembled.
      bool form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state);

      /// Utilities.
      void init_u_ext();
      void init_ext();
      
      void assemble_forms(int order, Func<double>** base_fns, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext,
        AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);

      /// Matrix volumetric forms - assemble the form.
      void assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext,
        AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);

      /// Vector volumetric forms - assemble the form.
      void assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext, 
        AsmList<Scalar>* current_als, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);

      Func<double>*** test_fns;

      friend class DiscreteProblemThreadAssembler<Scalar>;
      friend class DiscreteProblemDGAssembler<Scalar>;
    };
  }
}
#endif
