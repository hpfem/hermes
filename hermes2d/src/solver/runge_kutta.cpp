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

#include "solver/runge_kutta.h"
#include "discrete_problem.h"
#include "projections/ogprojection.h"
#include "projections/localprojection.h"
#include "weakform_library/weakforms_hcurl.h"
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    RungeKutta<Scalar>::RungeKutta(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> > spaces, ButcherTable* bt)
      : wf(wf), bt(bt), num_stages(bt->get_size()), stage_wf_right(bt->get_size() * spaces.size()),
      stage_wf_left(spaces.size()), start_from_zero_K_vector(false), block_diagonal_jacobian(false), residual_as_vector(true), iteration(0),
      freeze_jacobian(false), newton_tol(1e-6), newton_max_iter(20), newton_damping_coeff(1.0), newton_max_allowed_residual_norm(1e10)
    {
      for(unsigned int i = 0; i < spaces.size(); i++)
      {
        this->spaces.push_back(spaces.at(i));
        this->spaces_seqs.push_back(spaces.at(i)->get_seq());
      }
      for(unsigned int i = 0; i < spaces.size(); i++)
        this->spaces_mutable.push_back(spaces.at(i));

      if(bt==NULL)
        throw Exceptions::NullException(2);

      do_global_projections = true;

      matrix_right = create_matrix<Scalar>();
      matrix_left = create_matrix<Scalar>();
      vector_right = create_vector<Scalar>();
      // Create matrix solver.
      solver = create_linear_solver(matrix_right, vector_right);

      // Vector K_vector of length num_stages * ndof. will represent
      // the 'K_i' vectors in the usual R-K notation.
      K_vector = new Scalar[num_stages * Space<Scalar>::get_num_dofs(this->spaces)];

      // Vector u_ext_vec will represent h \sum_{j = 1}^s a_{ij} K_i.
      u_ext_vec = new Scalar[num_stages * Space<Scalar>::get_num_dofs(this->spaces)];

      // Vector for the left part of the residual.
      vector_left = new Scalar[num_stages*  Space<Scalar>::get_num_dofs(this->spaces)];

      this->stage_dp_left = NULL;
      this->stage_dp_right = NULL;
    }

    template<typename Scalar>
    RungeKutta<Scalar>::RungeKutta(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar> space, ButcherTable* bt)
      : wf(wf), bt(bt), num_stages(bt->get_size()), stage_wf_right(bt->get_size() * 1),
      stage_wf_left(1), start_from_zero_K_vector(false), block_diagonal_jacobian(false), residual_as_vector(true), iteration(0),
      freeze_jacobian(false), newton_tol(1e-6), newton_max_iter(20), newton_damping_coeff(1.0), newton_max_allowed_residual_norm(1e10)
    {
      this->spaces.push_back(space);
      this->spaces_seqs.push_back(space->get_seq());
      this->spaces_mutable.push_back(space);

      if(bt==NULL) throw Exceptions::NullException(2);

      do_global_projections = true;

      matrix_right = create_matrix<Scalar>();
      matrix_left = create_matrix<Scalar>();
      vector_right = create_vector<Scalar>();
      // Create matrix solver.
      solver = create_linear_solver(matrix_right, vector_right);

      // Vector K_vector of length num_stages * ndof. will represent
      // the 'K_i' vectors in the usual R-K notation.
      K_vector = new Scalar[num_stages * Space<Scalar>::get_num_dofs(spaces)];

      // Vector u_ext_vec will represent h \sum_{j = 1}^s a_{ij} K_i.
      u_ext_vec = new Scalar[num_stages * Space<Scalar>::get_num_dofs(spaces)];

      // Vector for the left part of the residual.
      vector_left = new Scalar[num_stages*  Space<Scalar>::get_num_dofs(spaces)];

      this->stage_dp_left = NULL;
      this->stage_dp_right = NULL;
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces)
    {
      bool delete_K_vector = false;
      for(unsigned int i = 0; i < spaces.size(); i++)
      {
        if(spaces[i]->get_seq() != this->spaces_seqs[i])
          delete_K_vector = true;
      }

      this->spaces = spaces;
      this->spaces_seqs.clear();
      for(unsigned int i = 0; i < spaces.size(); i++)
        this->spaces_seqs.push_back(spaces.at(i)->get_seq());
      this->spaces_mutable.clear();
      for(unsigned int i = 0; i < this->spaces.size(); i++)
        this->spaces_mutable.push_back(this->spaces.at(i));

      if(delete_K_vector)
      {
        delete [] K_vector;
        K_vector = new Scalar[num_stages * Space<Scalar>::get_num_dofs(this->spaces)];
        this->info("\tRunge-Kutta: K vectors are being set to zero, as the spaces changed during computation.");
        memset(K_vector, 0, num_stages * Space<Scalar>::get_num_dofs(this->spaces) * sizeof(Scalar));
      }
      delete [] u_ext_vec;
      u_ext_vec = new Scalar[num_stages * Space<Scalar>::get_num_dofs(this->spaces)];
      delete [] vector_left;
      vector_left = new Scalar[num_stages*  Space<Scalar>::get_num_dofs(this->spaces)];

      if(this->stage_dp_left != NULL)
        this->stage_dp_left->set_spaces(this->spaces);
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::set_space(SpaceSharedPtr<Scalar>& space)
    {
      bool delete_K_vector = false;
      if(space->get_seq() != this->spaces_seqs[0])
        delete_K_vector = true;

      this->spaces.clear();
			this->spaces.push_back(space);
      this->spaces_seqs.clear();
      this->spaces_seqs.push_back(space->get_seq());
      this->spaces_mutable.clear();
      this->spaces_mutable.push_back(space);

      if(delete_K_vector)
      {
        delete [] K_vector;
        K_vector = new Scalar[num_stages * Space<Scalar>::get_num_dofs(this->spaces)];
        this->info("\tRunge-Kutta: K vector is being set to zero, as the spaces changed during computation.");
        memset(K_vector, 0, num_stages * Space<Scalar>::get_num_dofs(this->spaces) * sizeof(Scalar));
      }
      delete [] u_ext_vec;
      u_ext_vec = new Scalar[num_stages * Space<Scalar>::get_num_dofs(this->spaces)];
      delete [] vector_left;
      vector_left = new Scalar[num_stages*  Space<Scalar>::get_num_dofs(this->spaces)];

      if(this->stage_dp_left != NULL)
        this->stage_dp_left->set_space(space);
    }

    template<typename Scalar>
    Hermes::vector<SpaceSharedPtr<Scalar> >& RungeKutta<Scalar>::get_spaces()
    {
      return this->spaces;
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::init()
    {
      this->create_stage_wf(spaces.size(), block_diagonal_jacobian);

      if(this->get_verbose_output())
      {
        this->stage_wf_left.set_verbose_output(true);
        this->stage_wf_right.set_verbose_output(true);
      }
      else
      {
        this->stage_wf_left.set_verbose_output(false);
        this->stage_wf_right.set_verbose_output(false);
      }

      // The tensor discrete problem is created in two parts. First, matrix_left is the Jacobian
      // matrix of the term coming from the left-hand side of the RK formula k_i = f(...). This is
      // a block-diagonal mass matrix. The corresponding part of the residual is obtained by multiplying
      // this block mass matrix with the tensor vector K. Next, matrix_right and vector_right are the Jacobian
      // matrix and residula vector coming from the function f(...). Of course the RK equation is assumed
      // in a form suitable for the Newton's method: k_i - f(...) = 0. At the end, matrix_left and vector_left
      // are added to matrix_right and vector_right, respectively.
      this->stage_dp_left = new DiscreteProblem<Scalar>(&stage_wf_left, spaces);

      // All Spaces of the problem.
      Hermes::vector<SpaceSharedPtr<Scalar> > stage_spaces_vector;

      // Create spaces for stage solutions K_i. This is necessary
      // to define a num_stages x num_stages block weak formulation.
      for (unsigned int i = 0; i < num_stages; i++)
        for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
          stage_spaces_vector.push_back(spaces[space_i]);

      this->stage_dp_right = new DiscreteProblem<Scalar>(&stage_wf_right, stage_spaces_vector);

      // Prepare residuals of stage solutions.
      if(!residual_as_vector)
        for (unsigned int i = 0; i < num_stages; i++)
          for(unsigned int sln_i = 0; sln_i < spaces.size(); sln_i++)
            residuals_vector.push_back(new Solution<Scalar>(spaces[sln_i]->get_mesh()));
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::set_start_from_zero_K_vector()
    {
      this->start_from_zero_K_vector = true;
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::set_residual_as_solutions()
    {
      this->residual_as_vector = false;
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::set_block_diagonal_jacobian()
    {
      this->block_diagonal_jacobian = true;
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::set_freeze_jacobian()
    {
      this->freeze_jacobian = true;
    }
    template<typename Scalar>
    void RungeKutta<Scalar>::set_tolerance(double newton_tol)
    {
      this->newton_tol = newton_tol;
    }
    template<typename Scalar>
    void RungeKutta<Scalar>::set_max_allowed_iterations(int newton_max_iter)
    {
      this->newton_max_iter = newton_max_iter;
    }
    template<typename Scalar>
    void RungeKutta<Scalar>::set_newton_damping_coeff(double newton_damping_coeff)
    {
      this->newton_damping_coeff = newton_damping_coeff;
    }
    template<typename Scalar>
    void RungeKutta<Scalar>::set_newton_max_allowed_residual_norm(double newton_max_allowed_residual_norm)
    {
      this->newton_max_allowed_residual_norm = newton_max_allowed_residual_norm;
    }

    template<typename Scalar>
    RungeKutta<Scalar>::~RungeKutta()
    {
      if(stage_dp_left != NULL)
        delete stage_dp_left;
      if(stage_dp_right != NULL)
        delete stage_dp_right;
      delete solver;
      delete matrix_right;
      delete matrix_left;
      delete vector_right;
      delete [] K_vector;
      delete [] u_ext_vec;
      delete [] vector_left;
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::use_local_projections()
    {
      do_global_projections = false;
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::multiply_as_diagonal_block_matrix(SparseMatrix<Scalar>* matrix, int num_blocks,
      Scalar* source_vec, Scalar* target_vec)
    {
      int size = matrix->get_size();
      for (int i = 0; i < num_blocks; i++)
      {
        matrix->multiply_with_vector(source_vec + i*size, target_vec + i*size);
      }
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::rk_time_step_newton(MeshFunctionSharedPtr<Scalar>  sln_time_prev,
                                          MeshFunctionSharedPtr<Scalar>  sln_time_new, MeshFunctionSharedPtr<Scalar>  error_fn)
    {
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns_time_prev = Hermes::vector<MeshFunctionSharedPtr<Scalar> >();
      slns_time_prev.push_back(sln_time_prev);
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns_time_new  = Hermes::vector<MeshFunctionSharedPtr<Scalar> >();
      slns_time_new.push_back(sln_time_new);
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > error_fns      = Hermes::vector<MeshFunctionSharedPtr<Scalar> >();
      error_fns.push_back(error_fn);
      return rk_time_step_newton(slns_time_prev, slns_time_new,
        error_fns);
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::rk_time_step_newton(Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns_time_prev,
                                          Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns_time_new,
                                          Hermes::vector<MeshFunctionSharedPtr<Scalar> > error_fns)
    {
      this->tick();

      int ndof = Space<Scalar>::get_num_dofs(spaces);

      if(this->stage_dp_left == NULL)
        this->init();

      // Creates the stage weak formulation.
      update_stage_wf(slns_time_prev);

      // Check whether the user provided a nonzero B2-row if he wants temporal error estimation.
      if(error_fns != Hermes::vector<MeshFunctionSharedPtr<Scalar> >() && bt->is_embedded() == false)
        throw Hermes::Exceptions::Exception("rk_time_step_newton(): R-K method must be embedded if temporal error estimate is requested.");

      info("\tRunge-Kutta: time step, time: %f, time step: %f", this->time, this->time_step);

      // Set the correct time to the essential boundary conditions.
      for (unsigned int stage_i = 0; stage_i < num_stages; stage_i++)
        Space<Scalar>::update_essential_bc_values(spaces_mutable, this->time + bt->get_C(stage_i)*this->time_step);

      // All Spaces of the problem.
      Hermes::vector<SpaceSharedPtr<Scalar> > stage_spaces_vector;

      // Create spaces for stage solutions K_i. This is necessary
      // to define a num_stages x num_stages block weak formulation.
      for (unsigned int i = 0; i < num_stages; i++)
        for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
        {
          typename Space<Scalar>::ReferenceSpaceCreator ref_space_creator(spaces[space_i], spaces[space_i]->get_mesh(), 0);
          stage_spaces_vector.push_back(ref_space_creator.create_ref_space());
        }

      Space<Scalar>::assign_dofs(stage_spaces_vector);
      this->stage_dp_right->set_spaces(stage_spaces_vector);

      // Zero utility vectors.
      if(start_from_zero_K_vector || !iteration)
        memset(K_vector, 0, num_stages * ndof * sizeof(Scalar));
      memset(u_ext_vec, 0, num_stages * ndof * sizeof(Scalar));
      memset(vector_left, 0, num_stages * ndof * sizeof(Scalar));

      // Assemble the block-diagonal mass matrix M of size ndof times ndof.
      // The corresponding part of the global residual vector is obtained
      // just by multiplication with the stage vector K.
      // FIXME: This should not be repeated if spaces have not changed.
      stage_dp_left->assemble(matrix_left);

      // The Newton's loop.
      double residual_norm;
      int it = 1;
      while (true)
      {
        // Prepare vector h\sum_{j = 1}^s a_{ij} K_j.
        prepare_u_ext_vec();

        // Reinitialize filters.
        if(this->filters_to_reinit.size() > 0)
        {
          Solution<Scalar>::vector_to_solutions(u_ext_vec, spaces, slns_time_new);

          for(unsigned int filters_i = 0; filters_i < this->filters_to_reinit.size(); filters_i++)
            filters_to_reinit.at(filters_i)->reinit();
        }

        // Residual corresponding to the stage derivatives k_i in the equation k_i - f(...) = 0.
        multiply_as_diagonal_block_matrix(matrix_left, num_stages, K_vector, vector_left);

        // Assemble the block Jacobian matrix of the stationary residual F.
        // Diagonal blocks are created even if empty, so that matrix_left can be added later.
        bool force_diagonal_blocks = true;
        stage_dp_right->set_RK(spaces.size(), force_diagonal_blocks, this->bt);
        stage_dp_right->assemble(u_ext_vec, NULL, vector_right);

        // Finalizing the residual vector.
        vector_right->add_vector(vector_left);

        // Multiply the residual vector with -1 since the matrix
        // equation reads J(Y^n) \deltaY^{n + 1} = -F(Y^n).
        vector_right->change_sign();
        if(this->output_rhsOn && (this->output_rhsIterations == -1 || this->output_rhsIterations >= it))
        {
          char* fileName = new char[this->RhsFilename.length() + 5];
          if(this->RhsFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
            sprintf(fileName, "%s%i.m", this->RhsFilename.c_str(), it);
          else
            sprintf(fileName, "%s%i", this->RhsFilename.c_str(), it);
          FILE* rhs_file = fopen(fileName, "w+");
          vector_right->dump(rhs_file, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
          fclose(rhs_file);
        }

        // Measure the residual norm.
        if(residual_as_vector)
          // Calculate the l2-norm of residual vector.
          residual_norm = Global<Scalar>::get_l2_norm(vector_right);
        else
        {
          // Translate residual vector into residual functions.
          Hermes::vector<bool> add_dir_lift_vector;
          add_dir_lift_vector.reserve(1);
          add_dir_lift_vector.push_back(false);
          Solution<Scalar>::vector_to_solutions(vector_right, stage_dp_right->get_spaces(), residuals_vector, false);
          
          Hermes::vector<MeshFunction<Scalar>*> meshFns;
          for(int i = 0; i < residuals_vector.size(); i++)
            meshFns.push_back(residuals_vector[i].get());
          residual_norm = Global<Scalar>::calc_norms(meshFns);
        }

        // Info for the user.
        if(it == 1)
          this->info("\tRunge-Kutta: Newton initial residual norm: %g", residual_norm);
        else
          this->info("\tRunge-Kutta: Newton iteration %d, residual norm: %g", it-1, residual_norm);

        // If maximum allowed residual norm is exceeded, fail.
        if(residual_norm > newton_max_allowed_residual_norm)
        {
          throw Exceptions::ValueException("residual norm", residual_norm, newton_max_allowed_residual_norm);
        }

        // If residual norm is within tolerance, or the maximum number
        // of iteration has been reached, or the problem is linear, then quit.
        if((residual_norm < newton_tol || it > newton_max_iter) && it > 1)
          break;

        bool rhs_only = (freeze_jacobian && it > 1);
        if(!rhs_only)
        {
          // Assemble the block Jacobian matrix of the stationary residual F
          // Diagonal blocks are created even if empty, so that matrix_left
          // can be added later.
          stage_dp_right->set_RK(spaces.size(), force_diagonal_blocks);
          stage_dp_right->assemble(u_ext_vec, matrix_right, NULL);

          // Adding the block mass matrix M to matrix_right. This completes the
          // resulting tensor Jacobian.
          matrix_right->add_sparse_to_diagonal_blocks(num_stages, matrix_left);

          if(this->output_matrixOn && (this->output_matrixIterations == -1 || this->output_matrixIterations >= it))
          {
            char* fileName = new char[this->matrixFilename.length() + 5];
            if(this->matrixFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
              sprintf(fileName, "%s%i.m", this->matrixFilename.c_str(), it);
            else
              sprintf(fileName, "%s%i", this->matrixFilename.c_str(), it);
            FILE* matrix_file = fopen(fileName, "w+");

            matrix_right->dump(matrix_file, this->matrixVarname.c_str(), this->matrixFormat, this->matrix_number_format);
            fclose(matrix_file);
          }

          matrix_right->finish();
        }
        else
          solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

        // Solve the linear system.
        if(!solver->solve())
          throw Exceptions::LinearMatrixSolverException();

        // Add \deltaK^{n + 1} to K^n.
        for (unsigned int i = 0; i < num_stages*ndof; i++)
          K_vector[i] += newton_damping_coeff * solver->get_sln_vector()[i];

        // Increase iteration counter.
        it++;
      }

      // If max number of iterations was exceeded, fail.
      if(it >= newton_max_iter)
      {
        this->tick();
        this->info("\tRunge-Kutta: time step duration: %f s.\n", this->last());
        throw Exceptions::ValueException("Newton iterations", it, newton_max_iter);
      }

      // Project previous time level solution on the stage space,
      // to be able to add them together. The result of the projection
      // will be stored in the vector coeff_vec.
      // FIXME - this projection is not needed when the
      //         spaces are the same (if spatial adaptivity is not used).
      Scalar* coeff_vec = new Scalar[ndof];
      if(do_global_projections)
      {
        OGProjection<Scalar> ogProjection;
        ogProjection.project_global(spaces, slns_time_prev, coeff_vec);
      }
      else
      {
        LocalProjection<Scalar> ogProjection;
        ogProjection.project_local(spaces, slns_time_prev, coeff_vec);
      }

      if(do_global_projections)
      {
        OGProjection<Scalar> ogProjection;
        ogProjection.project_global(spaces, slns_time_prev, coeff_vec);
      }
      else
      {
        LocalProjection<Scalar> ogProjection;
        ogProjection.project_local(spaces, slns_time_prev, coeff_vec);
      }

      // Calculate new time level solution in the stage space (u_{n + 1} = u_n + h \sum_{j = 1}^s b_j k_j).
      for (int i = 0; i < ndof; i++)
        for (unsigned int j = 0; j < num_stages; j++)
          coeff_vec[i] += this->time_step * bt->get_B(j) * K_vector[j * ndof + i];

      Solution<Scalar>::vector_to_solutions(coeff_vec, spaces, slns_time_new);

      // If error_fn is not NULL, use the B2-row in the Butcher's
      // table to calculate the temporal error estimate.
      if(error_fns != Hermes::vector<MeshFunctionSharedPtr<Scalar> >())
      {
        for (int i = 0; i < ndof; i++)
        {
          coeff_vec[i] = 0.;
          for (unsigned int j = 0; j < num_stages; j++)
            coeff_vec[i] += (bt->get_B(j) - bt->get_B2(j)) * K_vector[j * ndof + i];
          coeff_vec[i] *= this->time_step;
        }
        Solution<Scalar>::vector_to_solutions_common_dir_lift(coeff_vec, spaces, error_fns);
      }

      // Clean up.
      delete [] coeff_vec;

      iteration++;
      this->tick();
      this->info("\tRunge-Kutta: time step duration: %f s.\n", this->last());
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::rk_time_step_newton(Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns_time_prev,
                                          Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns_time_new)
    {
      return rk_time_step_newton(slns_time_prev, slns_time_new,
                          Hermes::vector<MeshFunctionSharedPtr<Scalar> >());
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::rk_time_step_newton(MeshFunctionSharedPtr<Scalar>  sln_time_prev, MeshFunctionSharedPtr<Scalar>  sln_time_new)
    {
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns_time_prev;;
      slns_time_prev.push_back(sln_time_prev);
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns_time_new;
      slns_time_new.push_back(sln_time_new);
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > error_fns;
      return rk_time_step_newton(slns_time_prev, slns_time_new, error_fns);
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::set_filters_to_reinit(Hermes::vector<Filter<Scalar>*> filters_to_reinit)
    {
      for(int i = 0; i < filters_to_reinit.size(); i++)
        this->filters_to_reinit.push_back(filters_to_reinit.at(i));
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::create_stage_wf(unsigned int size, bool block_diagonal_jacobian)
    {
      // Clear the WeakForms.
      stage_wf_left.delete_all();
      stage_wf_right.delete_all();

      // First let's do the mass matrix (only one block ndof times ndof).
      for(unsigned int component_i = 0; component_i < size; component_i++)
      {
        if(spaces[component_i]->get_type() == HERMES_H1_SPACE
           || spaces[component_i]->get_type() == HERMES_L2_SPACE)
        {
          MatrixFormVolL2<Scalar>* proj_form = new MatrixFormVolL2<Scalar>(component_i, component_i);
          proj_form->areas.push_back(HERMES_ANY);
          proj_form->scaling_factor = 1.0;
          proj_form->u_ext_offset = 0;
          stage_wf_left.add_matrix_form(proj_form);
        }
        if(spaces[component_i]->get_type() == HERMES_HDIV_SPACE
           || spaces[component_i]->get_type() == HERMES_HCURL_SPACE)
        {
          MatrixFormVolHCurl<Scalar>* proj_form = new MatrixFormVolHCurl<Scalar>(component_i, component_i);
          proj_form->areas.push_back(HERMES_ANY);
          proj_form->scaling_factor = 1.0;
          proj_form->u_ext_offset = 0;
          stage_wf_left.add_matrix_form(proj_form);
        }
      }

      // In the rest we will take the stationary jacobian and residual forms
      // (right-hand side) and use them to create a block Jacobian matrix of
      // size (num_stages*ndof times num_stages*ndof) and a block residual
      // vector of length num_stages*ndof.

      // Extracting volume and surface matrix and vector forms from the
      // original weak formulation.
      Hermes::vector<MatrixFormVol<Scalar> *> mfvol_base = wf->mfvol;
      Hermes::vector<MatrixFormSurf<Scalar> *> mfsurf_base = wf->mfsurf;
      Hermes::vector<VectorFormVol<Scalar> *> vfvol_base = wf->vfvol;
      Hermes::vector<VectorFormSurf<Scalar> *> vfsurf_base = wf->vfsurf;

      // Duplicate matrix volume forms, scale them according
      // to the Butcher's table, enhance them with additional
      // external solutions, and anter them as blocks to the
      // new stage Jacobian. If block_diagonal_jacobian = true
      // then only diagonal blocks are considered.
      for (unsigned int m = 0; m < mfvol_base.size(); m++)
      {
        for (unsigned int i = 0; i < num_stages; i++)
        {
          for (unsigned int j = 0; j < num_stages; j++)
          {
            if(block_diagonal_jacobian && i != j) continue;

            MatrixFormVol<Scalar>* mfv_ij = mfvol_base[m]->clone();

            mfv_ij->i = mfv_ij->i + i * spaces.size();
            mfv_ij->j = mfv_ij->j + j * spaces.size();

            mfv_ij->u_ext_offset = i * spaces.size();

            // Add the matrix form to the corresponding block of the
            // stage Jacobian matrix.
            stage_wf_right.add_matrix_form(mfv_ij);
          }
        }
      }

      // Duplicate matrix surface forms, enhance them with
      // additional external solutions, and anter them as
      // blocks of the stage Jacobian.
      for (unsigned int m = 0; m < mfsurf_base.size(); m++)
      {
        for (unsigned int i = 0; i < num_stages; i++)
        {
          for (unsigned int j = 0; j < num_stages; j++)
          {
            if(block_diagonal_jacobian && i != j) continue;

            MatrixFormSurf<Scalar>* mfs_ij = mfsurf_base[m]->clone();

            mfs_ij->i = mfs_ij->i + i * spaces.size();
            mfs_ij->j = mfs_ij->j + j * spaces.size();

            mfs_ij->u_ext_offset = i * spaces.size();

            // Add the matrix form to the corresponding block of the
            // stage Jacobian matrix.
            stage_wf_right.add_matrix_form_surf(mfs_ij);
          }
        }
      }

      // Duplicate vector volume forms, enhance them with
      // additional external solutions, and anter them as
      // blocks of the stage residual.
      for (unsigned int m = 0; m < vfvol_base.size(); m++)
      {
        for (unsigned int i = 0; i < num_stages; i++)
        {
          VectorFormVol<Scalar>* vfv_i = vfvol_base[m]->clone();

          vfv_i->i = vfv_i->i + i * spaces.size();

          vfv_i->scaling_factor = -1.0;
          vfv_i->u_ext_offset = i * spaces.size();

          // Add the matrix form to the corresponding block of the
          // stage Jacobian matrix.
          stage_wf_right.add_vector_form(vfv_i);
        }
      }

      // Duplicate vector surface forms, enhance them with
      // additional external solutions, and anter them as
      // blocks of the stage residual.
      for (unsigned int m = 0; m < vfsurf_base.size(); m++)
      {
        for (unsigned int i = 0; i < num_stages; i++)
        {
          VectorFormSurf<Scalar>* vfs_i = vfsurf_base[m]->clone();

          vfs_i->i = vfs_i->i + i * spaces.size();

          vfs_i->scaling_factor = -1.0;
          vfs_i->u_ext_offset = i * spaces.size();

          // Add the matrix form to the corresponding block of the
          // stage Jacobian matrix.
          stage_wf_right.add_vector_form_surf(vfs_i);
        }
      }
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::update_stage_wf(Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns_time_prev)
    {
      if(this->wf->global_integration_order_set)
      {
        this->stage_wf_left.set_global_integration_order(this->wf->global_integration_order);
        this->stage_wf_right.set_global_integration_order(this->wf->global_integration_order);
      }

      // Extracting volume and surface matrix and vector forms from the
      // 'right' weak formulation.
      Hermes::vector<MatrixFormVol<Scalar> *> mfvol = stage_wf_right.mfvol;
      Hermes::vector<MatrixFormSurf<Scalar> *> mfsurf = stage_wf_right.mfsurf;
      Hermes::vector<VectorFormVol<Scalar> *> vfvol = stage_wf_right.vfvol;
      Hermes::vector<VectorFormSurf<Scalar> *> vfsurf = stage_wf_right.vfsurf;

      stage_wf_right.ext.clear();

      for(unsigned int slns_time_prev_i = 0; slns_time_prev_i < slns_time_prev.size(); slns_time_prev_i++)
          stage_wf_right.ext.push_back(slns_time_prev[slns_time_prev_i]);

      // Duplicate matrix volume forms, scale them according
      // to the Butcher's table, enhance them with additional
      // external solutions, and anter them as blocks to the
      // new stage Jacobian. If block_diagonal_jacobian = true
      // then only diagonal blocks are considered.
      for (unsigned int m = 0; m < mfvol.size(); m++)
      {
        MatrixFormVol<Scalar> *mfv_ij = mfvol[m];
        mfv_ij->scaling_factor = -this->time_step * bt->get_A(mfv_ij->i / spaces.size(), mfv_ij->j / spaces.size());
        mfv_ij->set_current_stage_time(this->time + bt->get_C(mfv_ij->i / spaces.size()) * this->time_step);
      }

      // Duplicate matrix surface forms, enhance them with
      // additional external solutions, and anter them as
      // blocks of the stage Jacobian.
      for (unsigned int m = 0; m < mfsurf.size(); m++)
      {
        MatrixFormSurf<Scalar> *mfs_ij = mfsurf[m];
        mfs_ij->scaling_factor = -this->time_step * bt->get_A(mfs_ij->i / spaces.size(), mfs_ij->j / spaces.size());
        mfs_ij->set_current_stage_time(this->time + bt->get_C(mfs_ij->i / spaces.size()) * this->time_step);
      }

      // Duplicate vector volume forms, enhance them with
      // additional external solutions, and anter them as
      // blocks of the stage residual.
      for (unsigned int m = 0; m < vfvol.size(); m++)
      {
        VectorFormVol<Scalar>* vfv_i = vfvol[m];
        vfv_i->set_current_stage_time(this->time + bt->get_C(vfv_i->i / spaces.size())*this->time_step);
      }

      // Duplicate vector surface forms, enhance them with
      // additional external solutions, and anter them as
      // blocks of the stage residual.
      for (unsigned int m = 0; m < vfsurf.size(); m++)
      {
        VectorFormSurf<Scalar>* vfs_i = vfsurf[m];
        vfs_i->set_current_stage_time(this->time + bt->get_C(vfs_i->i / spaces.size())*this->time_step);
      }
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::prepare_u_ext_vec()
    {
      unsigned int ndof = Space<Scalar>::get_num_dofs(spaces);
      for (unsigned int stage_i = 0; stage_i < num_stages; stage_i++)
      {
        unsigned int running_space_ndofs = 0;
        for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
        {
          for (int idx = 0; idx < spaces[space_i]->get_num_dofs(); idx++)
          {
            Scalar increment = 0;
            for (unsigned int stage_j = 0; stage_j < num_stages; stage_j++)
              increment += bt->get_A(stage_i, stage_j) * K_vector[stage_j * ndof + running_space_ndofs + idx];
            u_ext_vec[stage_i * ndof + running_space_ndofs + idx] = this->time_step * increment;
          }
          running_space_ndofs += spaces[space_i]->get_num_dofs();
        }
      }
    }
    template class HERMES_API RungeKutta<double>;
    template class HERMES_API RungeKutta<std::complex<double> >;
  }
}