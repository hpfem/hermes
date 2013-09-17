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

#ifndef __H2D_DISCRETE_PROBLEM_INTEGRATION_ORDER_CALCULATOR_H
#define __H2D_DISCRETE_PROBLEM_INTEGRATION_ORDER_CALCULATOR_H

#include "hermes_common.h"
#include "adapt/adapt.h"
#include "graph.h"
#include "forms.h"
#include "weakform/weakform.h"
#include "function/function.h"
#include "neighbor_search.h"
#include "refinement_selectors/selector.h"
#include "exceptions.h"
#include "mixins2d.h"
#include "discrete_problem_helpers.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class PrecalcShapeset;
    /// @ingroup inner
    /// DiscreteProblemIntegrationOrderCalculator class.
    /// \brief Provides methods of integration order calculation.
    template<typename Scalar>
    class HERMES_API DiscreteProblemIntegrationOrderCalculator :
      public Hermes::Hermes2D::Mixins::DiscreteProblemRungeKutta<Scalar>
    {
    private:
      DiscreteProblemIntegrationOrderCalculator(DiscreteProblemSelectiveAssembler<Scalar>* selectiveAssembler);

      /// Adjusts order to refmaps.
      void adjust_order_to_refmaps(Form<Scalar> *form, int& order, Hermes::Ord* o, RefMap** current_refmaps);

      /// Matrix volumetric forms - calculate the integration order.
      int calc_order_matrix_form(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, MatrixForm<Scalar>* mfv, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, Traverse::State* current_state);

      /// Vector volumetric forms - calculate the integration order.
      int calc_order_vector_form(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, VectorForm<Scalar>* mfv, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, Traverse::State* current_state);

      /// Order calculation.
      int calculate_order(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State* state, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, WeakForm<Scalar>* current_wf);
      
      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Calculates orders for external functions.
      void init_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, Func<Hermes::Ord>** oext, Solution<Scalar>** current_u_ext, Traverse::State* current_state);
      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Cleans up after init_ext_orders.
      template<typename FormType>
      void deinit_ext_orders(Form<Scalar> *form, FormType** oi, FormType** oext);

      /// Calculates integration order for DG matrix forms.
      int calc_order_dg_matrix_form(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State* state, MatrixFormDG<Scalar>* mfDG, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, bool neighbor_supp_u, bool neighbor_supp_v, NeighborSearch<Scalar>** neighbor_searches);

      /// Calculates integration order for DG vector forms.
      int calc_order_dg_vector_form(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State* state, VectorFormDG<Scalar>* vfDG, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, bool neighbor_supp_v, NeighborSearch<Scalar>** neighbor_searches);

      /// One external function for DG.
      DiscontinuousFunc<Hermes::Ord>* init_ext_fn_ord(NeighborSearch<Scalar>* ns, MeshFunctionSharedPtr<Scalar> fu);

      /// Initialize orders of external functions for DG forms.
      DiscontinuousFunc<Hermes::Ord>** init_ext_fns_ord(Hermes::vector<MeshFunctionSharedPtr<Scalar> > &ext,
        NeighborSearch<Scalar>** neighbor_searches);

      /// For selective assembling.
      DiscreteProblemSelectiveAssembler<Scalar>* selectiveAssembler;

      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemThreadAssembler;
    };
  }
}
#endif
