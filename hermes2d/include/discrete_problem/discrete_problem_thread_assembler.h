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
#include "discrete_problem_integration_order_calculator.h"
#include "discrete_problem_selective_assembler.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup inner
    /// Discrete problem thread assembler class
    /// \brief This class is a one-thread (non-DG) assembly worker.
    ///
    template<typename Scalar>
    class HERMES_API DiscreteProblemThreadAssembler : 
      public Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemRungeKutta<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemMatrixVector<Scalar>
    {
    private:
      DiscreteProblemThreadAssembler(DiscreteProblemSelectiveAssembler<Scalar>* selectiveAssembler);
      ~DiscreteProblemThreadAssembler();

      /// Initialization of all structures concerning space - assembly lists, precalculated shapesets, ..
      void init_spaces(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      /// Initialization of the weak formulation.
      void set_weak_formulation(WeakForm<Scalar>* wf);
      /// Initialization of previous iterations for non-linear solvers.
      void init_u_ext(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Solution<Scalar>** u_ext_sln);

      /// Initializes the Transformable array for doing transformations.
      void init_assembling(Solution<Scalar>** u_ext_sln, const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, bool nonlinear, bool add_dirichlet_lift);

      /// Initialize Func storages.
      void init_funcs();
      /// De-initialize Func storages.
      void deinit_funcs();
      /// Initializitation of u-ext values into Funcs
      void init_u_ext_values(int order);
      /// Initializitation of ext values into Funcs
      void init_ext_values(Func<Scalar>** target_array, Hermes::vector<MeshFunctionSharedPtr<Scalar> >& ext, Hermes::vector<UExtFunctionSharedPtr<Scalar> >& u_ext_fns, int order, Func<Scalar>** u_ext_func, Geom<double>* geometry);

      /// Sets active elements & transformations
      void init_assembling_one_state(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State* current_state);
      /// Assemble the state.
      void assemble_one_state();
      /// Matrix volumetric forms - assemble the form.
      void assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns,
        AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);
      /// Vector volumetric forms - assemble the form.
      void assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, AsmList<Scalar>* current_als,
        int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);
      /// De-initialization of 1 state assembly
      void deinit_assembling_one_state();
      
      /// De-initialization.
      void deinit_assembling();
      
      /// Free all data.
      void free();
      /// Free space-related data.
      void free_spaces();
      /// Free weak formulation data.
      void free_weak_formulation();
      /// Free nonlinearities-related data.
      void free_u_ext();

      PrecalcShapeset** pss;
      RefMap** refmaps;
      RefMap* rep_refmap;
      Solution<Scalar>** u_ext;
      Hermes::vector<Transformable *> fns;
      
      /// For selective reassembling.
      DiscreteProblemSelectiveAssembler<Scalar>* selectiveAssembler;
      
      /// Currently assembled state.
      Traverse::State* current_state;

      /// Integration orders for the currently assembled state.
      /// - calculator
      DiscreteProblemIntegrationOrderCalculator<Scalar> integrationOrderCalculator;
      /// - volumetric
      int order;
      /// - surface
      int orderSurface[H2D_MAX_NUMBER_EDGES];

      /// Holding values formerly held by cache record.
      void init_calculation_variables();
      void deinit_calculation_variables();
      Func<double>* funcs[H2D_MAX_COMPONENTS][H2D_MAX_LOCAL_BASIS_SIZE];
      Func<double>* funcsSurface[H2D_MAX_NUMBER_EDGES][H2D_MAX_COMPONENTS][H2D_MAX_LOCAL_BASIS_SIZE];
      Geom<double>* geometry;
      Geom<double>* geometrySurface[H2D_MAX_NUMBER_EDGES];
      double jacobian_x_weights[H2D_MAX_INTEGRATION_POINTS_COUNT];
      double jacobian_x_weightsSurface[H2D_MAX_NUMBER_EDGES][H2D_MAX_INTEGRATION_POINTS_COUNT];
      int n_quadrature_points;
      int n_quadrature_pointsSurface[H2D_MAX_NUMBER_EDGES];

      /// Ext values storage.
      Func<Scalar>* u_ext_funcs[H2D_MAX_COMPONENTS];
      Func<Scalar>** ext_funcs;
      int ext_funcs_allocated_size;
      Func<Scalar>** ext_funcs_local;
      int ext_funcs_local_allocated_size;

      AsmList<Scalar> als[H2D_MAX_COMPONENTS];
      AsmList<Scalar> alsSurface[H2D_MAX_NUMBER_EDGES][H2D_MAX_COMPONENTS];
      int spaces_size;
      bool nonlinear, add_dirichlet_lift;

      friend class DiscreteProblem<Scalar>;
      friend class DiscreteProblemDGAssembler<Scalar>;
    };
  }
}
#endif
