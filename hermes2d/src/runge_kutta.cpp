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

#include <typeinfo>
#include "runge_kutta.h"
#include "discrete_problem.h"
#include "ogprojection.h"
#include "weakform_library/weakforms_hcurl.h"
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    RungeKutta<Scalar>::RungeKutta(DiscreteProblem<Scalar>* dp, ButcherTable* bt, Hermes::MatrixSolverType matrix_solver_type, bool start_from_zero_K_vector, bool residual_as_vector)
      : dp(dp), bt(bt), num_stages(bt->get_size()), stage_wf_right(bt->get_size() * dp->get_spaces().size()),
      stage_wf_left(dp->get_spaces().size()), start_from_zero_K_vector(start_from_zero_K_vector), residual_as_vector(residual_as_vector), iteration(0) , matrix_solver_type(matrix_solver_type)
    {
      _F_
      if (dp==NULL) throw Exceptions::NullException(1);
      if (bt==NULL) throw Exceptions::NullException(2);

      matrix_right = create_matrix<Scalar>(matrix_solver_type);
      matrix_left = create_matrix<Scalar>(matrix_solver_type);
      vector_right = create_vector<Scalar>(matrix_solver_type);
      // Create matrix solver.
      solver = create_linear_solver(matrix_solver_type, matrix_right, vector_right);

      // Vector K_vector of length num_stages * ndof. will represent
      // the 'K_i' vectors in the usual R-K notation.
      K_vector = new Scalar[num_stages * dp->get_num_dofs()];

      // Vector u_ext_vec will represent h \sum_{j = 1}^s a_{ij} K_i.
      u_ext_vec = new Scalar[num_stages * dp->get_num_dofs()];

      // Vector for the left part of the residual.
      vector_left = new Scalar[num_stages*  dp->get_num_dofs()];
    }

    template<typename Scalar>
    RungeKutta<Scalar>::~RungeKutta()
    {
      delete solver;
      delete matrix_right;
      delete matrix_left;
      delete vector_right;
      delete [] K_vector;
      delete [] u_ext_vec;
      delete [] vector_left;
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
    void RungeKutta<Scalar>::rk_time_step_newton(double current_time, double time_step, Solution<Scalar>* sln_time_prev,
                                          Solution<Scalar>* sln_time_new, Solution<Scalar>* error_fn,
                                          bool freeze_jacobian, bool block_diagonal_jacobian,
                                          bool verbose, double newton_tol, int newton_max_iter,
                                          double newton_damping_coeff, double newton_max_allowed_residual_norm)
    {
      Hermes::vector<Solution<Scalar>*> slns_time_prev = Hermes::vector<Solution<Scalar>*>();
      slns_time_prev.push_back(sln_time_prev);
      Hermes::vector<Solution<Scalar>*> slns_time_new  = Hermes::vector<Solution<Scalar>*>();
      slns_time_new.push_back(sln_time_new);
      Hermes::vector<Solution<Scalar>*> error_fns      = Hermes::vector<Solution<Scalar>*>();
      error_fns.push_back(error_fn);
      return rk_time_step_newton(current_time, time_step, slns_time_prev, slns_time_new,
        error_fns, freeze_jacobian, block_diagonal_jacobian, verbose, newton_tol, newton_max_iter,
        newton_damping_coeff, newton_max_allowed_residual_norm);
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::rk_time_step_newton(double current_time, double time_step,
                                          Hermes::vector<Solution<Scalar>*> slns_time_prev,
                                          Hermes::vector<Solution<Scalar>*> slns_time_new,
                                          Hermes::vector<Solution<Scalar>*> error_fns,
                                          bool freeze_jacobian, bool block_diagonal_jacobian,
                                          bool verbose, double newton_tol,
                                          int newton_max_iter, double newton_damping_coeff,
                                          double newton_max_allowed_residual_norm)
    {
      // Check whether the user provided a nonzero B2-row if he wants temporal error estimation.
      if(error_fns != Hermes::vector<Solution<Scalar>*>() && bt->is_embedded() == false)
        error("rk_time_step_newton(): R-K method must be embedded if temporal error estimate is requested.");

      // All Spaces of the problem.
      Hermes::vector<Space<Scalar>*> stage_spaces_vector;

      // Create spaces for stage solutions K_i. This is necessary
      // to define a num_stages x num_stages block weak formulation.
      for (unsigned int i = 0; i < num_stages; i++)
        for(unsigned int space_i = 0; space_i < dp->get_spaces().size(); space_i++)
          stage_spaces_vector.push_back(dp->get_space(space_i)->dup(dp->get_space(space_i)->get_mesh()));

      int ndof = dp->get_num_dofs();

      // Creates the stage weak formulation.
      create_stage_wf(dp->get_spaces().size(), current_time, time_step,
                      slns_time_prev, block_diagonal_jacobian);

      // Set the correct time to the essential boundary conditions.
      for (unsigned int stage_i = 0; stage_i < num_stages; stage_i++)
        Space<Scalar>::update_essential_bc_values(dp->get_spaces(), current_time + bt->get_C(stage_i)*time_step);

      // The tensor discrete problem is created in two parts. First, matrix_left is the Jacobian
      // matrix of the term coming from the left-hand side of the RK formula k_i = f(...). This is
      // a block-diagonal mass matrix. The corresponding part of the residual is obtained by multiplying
      // this block mass matrix with the tensor vector K. Next, matrix_right and vector_right are the Jacobian
      // matrix and residula vector coming from the function f(...). Of course the RK equation is assumed
      // in a form suitable for the Newton's method: k_i - f(...) = 0. At the end, matrix_left and vector_left
      // are added to matrix_right and vector_right, respectively.
      DiscreteProblem<Scalar> stage_dp_left(&stage_wf_left, dp->get_spaces());
      DiscreteProblem<Scalar> stage_dp_right(&stage_wf_right, stage_spaces_vector);
      stage_dp_right.set_RK(dp->get_spaces().size());

      // Prepare residuals of stage solutions.
      Hermes::vector<Solution<Scalar>*> residuals_vector;
      // A technical workaround.
      Hermes::vector<bool> add_dir_lift;
      for (unsigned int i = 0; i < num_stages; i++)
      {
        for(unsigned int sln_i = 0; sln_i < dp->get_spaces().size(); sln_i++)
        {
          residuals_vector.push_back(new Solution<Scalar>(dp->get_space(sln_i)->get_mesh()));
          add_dir_lift.push_back(false);
        }
      }

      // Zero utility vectors.
      if(start_from_zero_K_vector || !iteration)
        memset(K_vector, 0, num_stages * ndof * sizeof(Scalar));
      memset(u_ext_vec, 0, num_stages * ndof * sizeof(Scalar));
      memset(vector_left, 0, num_stages * ndof * sizeof(Scalar));

      // Assemble the block-diagonal mass matrix M of size ndof times ndof.
      // The corresponding part of the global residual vector is obtained
      // just by multiplication with the stage vector K.
      // FIXME: This should not be repeated if spaces have not changed.
      stage_dp_left.assemble(matrix_left, NULL);

      // The Newton's loop.
      double residual_norm;
      int it = 1;
      while (true)
      {
        // Prepare vector h\sum_{j = 1}^s a_{ij} K_j.
        prepare_u_ext_vec(time_step);

        // Residual corresponding to the stage derivatives k_i in the equation k_i - f(...) = 0.
        multiply_as_diagonal_block_matrix(matrix_left, num_stages, K_vector, vector_left);

        // Assemble the block Jacobian matrix of the stationary residual F.
        // Diagonal blocks are created even if empty, so that matrix_left can be added later.
        bool force_diagonal_blocks = true;
        bool add_dir_lift = false;
        stage_dp_right.assemble(u_ext_vec, NULL, vector_right, force_diagonal_blocks, add_dir_lift);

        // Finalizing the residual vector.
        vector_right->add_vector(vector_left);

        // Multiply the residual vector with -1 since the matrix
        // equation reads J(Y^n) \deltaY^{n + 1} = -F(Y^n).
        vector_right->change_sign();

        // Measure the residual norm.
        if (residual_as_vector)
          // Calculate the l2-norm of residual vector.
          residual_norm = Global<Scalar>::get_l2_norm(vector_right);
        else
        {
          // Translate residual vector into residual functions.
          Hermes::vector<bool> add_dir_lift_vector;
          add_dir_lift_vector.reserve(1);
          add_dir_lift_vector.push_back(add_dir_lift);
          Solution<Scalar>::vector_to_solutions(vector_right, stage_dp_right.get_spaces(),
            residuals_vector, add_dir_lift_vector);
          residual_norm = Global<Scalar>::calc_norms(residuals_vector);
        }

        // Info for the user.
        if (it == 1)
        {
          if(verbose)
            info("---- Newton initial residual norm: %g", residual_norm);
        }
        else
          if(verbose)
            info("---- Newton iter %d, residual norm: %g", it-1, residual_norm);

        // If maximum allowed residual norm is exceeded, fail.
        if (residual_norm > newton_max_allowed_residual_norm)
        {
          throw Exceptions::ValueException("residual norm", residual_norm, newton_max_allowed_residual_norm);
        }

        // If residual norm is within tolerance, or the maximum number
        // of iteration has been reached, or the problem is linear, then quit.
        if ((residual_norm < newton_tol || it > newton_max_iter) && it > 1)
          break;

        bool rhs_only = (freeze_jacobian && it > 1);
        if (!rhs_only)
        {
          // Assemble the block Jacobian matrix of the stationary residual F
          // Diagonal blocks are created even if empty, so that matrix_left
          // can be added later.
          stage_dp_right.assemble(u_ext_vec, matrix_right, NULL, force_diagonal_blocks, add_dir_lift);

          // Adding the block mass matrix M to matrix_right. This completes the
          // resulting tensor Jacobian.
          matrix_right->add_sparse_to_diagonal_blocks(num_stages, matrix_left);
          matrix_right->finish();
        }
        else
          solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

        // Solve the linear system.
        if(!solver->solve())
          throw Exceptions::LinearSolverException();

        // Add \deltaK^{n + 1} to K^n.
        for (unsigned int i = 0; i < num_stages*ndof; i++)
          K_vector[i] += newton_damping_coeff * solver->get_sln_vector()[i];

        // Increase iteration counter.
        it++;
      }

      // If max number of iterations was exceeded, fail.
      if (it >= newton_max_iter)
      {
        throw Exceptions::ValueException("newton iterations", it, newton_max_iter);
      }

      // Project previous time level solution on the stage space,
      // to be able to add them together. The result of the projection
      // will be stored in the vector coeff_vec.
      // FIXME - this projection is slow and it is not needed when the
      //         spaces are the same (if spatial adaptivity does not take place).
      Scalar* coeff_vec = new Scalar[ndof];
      OGProjection<Scalar>::project_global(dp->get_spaces(), slns_time_prev, coeff_vec);

      // Calculate new time level solution in the stage space (u_{n + 1} = u_n + h \sum_{j = 1}^s b_j k_j).
      for (int i = 0; i < ndof; i++)
        for (unsigned int j = 0; j < num_stages; j++)
          coeff_vec[i] += time_step * bt->get_B(j) * K_vector[j * ndof + i];

      Solution<Scalar>::vector_to_solutions(coeff_vec, dp->get_spaces(), slns_time_new);

      // If error_fn is not NULL, use the B2-row in the Butcher's
      // table to calculate the temporal error estimate.
      if (error_fns != Hermes::vector<Solution<Scalar>*>())
      {
        for (int i = 0; i < ndof; i++)
        {
          coeff_vec[i] = 0;
          for (unsigned int j = 0; j < num_stages; j++)
            coeff_vec[i] += (bt->get_B(j) - bt->get_B2(j)) * K_vector[j * ndof + i];
          coeff_vec[i] *= time_step;
        }
        Solution<Scalar>::vector_to_solutions(coeff_vec, dp->get_spaces(), error_fns, add_dir_lift);
      }

      // Delete stage spaces.
      for (unsigned int i = 0; i < num_stages; i++)
        delete stage_spaces_vector[i];

      // Delete all residuals.
      for (unsigned int i = 0; i < num_stages; i++)
        delete residuals_vector[i];

      // Clean up.
      delete [] coeff_vec;

      iteration++;
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::rk_time_step_newton(double current_time, double time_step,
                                          Hermes::vector<Solution<Scalar>*> slns_time_prev,
                                          Hermes::vector<Solution<Scalar>*> slns_time_new,
                                          bool freeze_jacobian, bool block_diagonal_jacobian,
                                          bool verbose, double newton_tol, int newton_max_iter,
                                          double newton_damping_coeff,
                                          double newton_max_allowed_residual_norm)
    {
      return rk_time_step_newton(current_time, time_step, slns_time_prev, slns_time_new,
                          Hermes::vector<Solution<Scalar>*>(), freeze_jacobian, block_diagonal_jacobian,
                          verbose, newton_tol, newton_max_iter, newton_damping_coeff,
                          newton_max_allowed_residual_norm);
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::rk_time_step_newton(double current_time, double time_step, Solution<Scalar>* sln_time_prev,
                                          Solution<Scalar>* sln_time_new,
                                          bool freeze_jacobian, bool block_diagonal_jacobian,
                                          bool verbose, double newton_tol, int newton_max_iter,
                                          double newton_damping_coeff, double newton_max_allowed_residual_norm)
    {
      Hermes::vector<Solution<Scalar>*> slns_time_prev = Hermes::vector<Solution<Scalar>*>();
      slns_time_prev.push_back(sln_time_prev);
      Hermes::vector<Solution<Scalar>*> slns_time_new  = Hermes::vector<Solution<Scalar>*>();
      slns_time_new.push_back(sln_time_new);
      Hermes::vector<Solution<Scalar>*> error_fns      = Hermes::vector<Solution<Scalar>*>();
      return rk_time_step_newton(current_time, time_step, slns_time_prev, slns_time_new,
                          error_fns, freeze_jacobian, block_diagonal_jacobian, verbose, newton_tol,
                          newton_max_iter, newton_damping_coeff, newton_max_allowed_residual_norm);
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::create_stage_wf(unsigned int size, double current_time, double time_step,
                                             Hermes::vector<Solution<Scalar>*> slns_time_prev,
                                             bool block_diagonal_jacobian)
    {
      // Clear the WeakForms.
      stage_wf_left.delete_all();
      stage_wf_right.delete_all();

      // First let's do the mass matrix (only one block ndof times ndof).
      for(unsigned int component_i = 0; component_i < size; component_i++)
      {
        if(dp->get_spaces()[component_i]->get_type() == HERMES_H1_SPACE
           || dp->get_spaces()[component_i]->get_type() == HERMES_L2_SPACE)
        {
          MatrixFormVolL2<Scalar>* proj_form = new MatrixFormVolL2<Scalar>(component_i, component_i);
          proj_form->areas.push_back(HERMES_ANY);
          proj_form->scaling_factor = 1.0;
          proj_form->u_ext_offset = 0;
          stage_wf_left.add_matrix_form(proj_form);
        }
        if(dp->get_spaces()[component_i]->get_type() == HERMES_HDIV_SPACE
           || dp->get_spaces()[component_i]->get_type() == HERMES_HCURL_SPACE)
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

      // Original weak formulation.
      WeakForm<Scalar>* wf = dp->get_weak_formulation();

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
            if (block_diagonal_jacobian && i != j) continue;

            MatrixFormVol<Scalar>* mfv_ij = mfvol_base[m]->clone();

            mfv_ij->i = mfv_ij->i + i * dp->get_spaces().size();
            mfv_ij->j = mfv_ij->j + j * dp->get_spaces().size();

            mfv_ij->scaling_factor = -time_step * bt->get_A(i, j);

            mfv_ij->u_ext_offset = i * dp->get_spaces().size();

            for(unsigned int slns_time_prev_i = 0; slns_time_prev_i < slns_time_prev.size(); slns_time_prev_i++)
              mfv_ij->ext.push_back(slns_time_prev[slns_time_prev_i]);

            mfv_ij->set_current_stage_time(current_time + bt->get_C(i)*time_step);

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
            if (block_diagonal_jacobian && i != j) continue;

            MatrixFormSurf<Scalar>* mfs_ij = mfsurf_base[m]->clone();

            mfs_ij->i = mfs_ij->i + i * dp->get_spaces().size();
            mfs_ij->j = mfs_ij->j + j * dp->get_spaces().size();

            mfs_ij->scaling_factor = -time_step * bt->get_A(i, j);

            mfs_ij->u_ext_offset = i * dp->get_spaces().size();

            for(unsigned int slns_time_prev_i = 0; slns_time_prev_i < slns_time_prev.size(); slns_time_prev_i++)
              mfs_ij->ext.push_back(slns_time_prev[slns_time_prev_i]);

            mfs_ij->set_current_stage_time(current_time + bt->get_C(i)*time_step);

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

          vfv_i->i = vfv_i->i + i * dp->get_spaces().size();

          vfv_i->scaling_factor = -1.0;
          vfv_i->u_ext_offset = i * dp->get_spaces().size();

          for(unsigned int slns_time_prev_i = 0; slns_time_prev_i < slns_time_prev.size(); slns_time_prev_i++)
            vfv_i->ext.push_back(slns_time_prev[slns_time_prev_i]);

          vfv_i->set_current_stage_time(current_time + bt->get_C(i)*time_step);

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

          vfs_i->i = vfs_i->i + i * dp->get_spaces().size();

          vfs_i->scaling_factor = -1.0;
          vfs_i->u_ext_offset = i * dp->get_spaces().size();

          for(unsigned int slns_time_prev_i = 0; slns_time_prev_i < slns_time_prev.size(); slns_time_prev_i++)
            vfs_i->ext.push_back(slns_time_prev[slns_time_prev_i]);

          vfs_i->set_current_stage_time(current_time + bt->get_C(i)*time_step);

          // Add the matrix form to the corresponding block of the
          // stage Jacobian matrix.
          stage_wf_right.add_vector_form_surf(vfs_i);
        }
      }
    }

    template<typename Scalar>
    void RungeKutta<Scalar>::prepare_u_ext_vec(double time_step)
    {
      unsigned int ndof = dp->get_num_dofs();
      for (unsigned int stage_i = 0; stage_i < num_stages; stage_i++)
      {
        unsigned int running_space_ndofs = 0;
        for(unsigned int space_i = 0; space_i < dp->get_spaces().size(); space_i++)
        {
          for (int idx = 0; idx < dp->get_space(space_i)->get_num_dofs(); idx++)
          {
            Scalar increment = 0;
            for (unsigned int stage_j = 0; stage_j < num_stages; stage_j++)
              increment += bt->get_A(stage_i, stage_j) * K_vector[stage_j * ndof + running_space_ndofs + idx];
            u_ext_vec[stage_i * ndof + running_space_ndofs + idx] = time_step * increment;
          }
          running_space_ndofs += dp->get_space(space_i)->get_num_dofs();
        }
      }
    }
    template class HERMES_API RungeKutta<double>;
    template class HERMES_API RungeKutta<std::complex<double> >;
  }
}
