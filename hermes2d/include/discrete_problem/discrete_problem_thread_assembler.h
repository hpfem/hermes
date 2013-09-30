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
#include "discrete_problem_cache.h"
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
      public Hermes::Hermes2D::Mixins::DiscreteProblemCacheSettings,
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

      /// Sets active elements & transformations
      void init_assembling_one_state(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State* current_state);
      /// Takes the global cache and according to the current_state, it
      /// either gets or gets a record.
      void handle_cache(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, DiscreteProblemCache<Scalar>* cache);
      /// Assemble the state.
      void assemble_one_state();
      /// Matrix volumetric forms - assemble the form.
      void assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext,
        AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);
      /// Vector volumetric forms - assemble the form.
      void assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext, 
        AsmList<Scalar>* current_als, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);

      /// Delete the cache record if do_not_use_cache etc.
      void deinit_assembling_one_state();
      /// Delete the cache record if do_not_use_cache etc.
      void deinit_assembling();
      
      /// Free all data.
      void free();
      /// Free space-related data.
      void free_spaces();
      /// Free weak formulation data.
      void free_weak_formulation();
      /// Free nonlinearities-related data.
      void free_u_ext();

      Func<Scalar>** init_u_ext_values(int order);
      void deinit_u_ext_values(Func<Scalar>** u_ext_func);

      Func<Scalar>** init_ext_values(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& ext, int order, Func<Scalar>** u_ext_func);
      void deinit_ext_values(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& ext, Func<Scalar>** ext_func);

      PrecalcShapeset** pss;
      RefMap** refmaps;
      Solution<Scalar>** u_ext;
      AsmList<Scalar>** als;
      AsmList<Scalar>*** alsSurface;
      Hermes::vector<Transformable *> fns;  
      int spaces_size;
      bool nonlinear, add_dirichlet_lift;

      Traverse::State* current_state;
      typename DiscreteProblemCache<Scalar>::CacheRecord* current_cache_record;
      
      /// For selective reassembling.
      DiscreteProblemSelectiveAssembler<Scalar>* selectiveAssembler;

      /// Integration order calculator.
      DiscreteProblemIntegrationOrderCalculator<Scalar> integrationOrderCalculator;

      friend class DiscreteProblem<Scalar>;
      friend class DiscreteProblemDGAssembler<Scalar>;
    };
  }
}
#endif
