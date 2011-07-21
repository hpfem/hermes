// This file is part of Hermes2D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#include "newton_solver.h"
#include "hermes_common.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    double NewtonSolver<Scalar>::max_allowed_residual_norm = 1E6;

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(DiscreteProblem<Scalar>* dp) : NonlinearSolver<Scalar>(dp), kept_jacobian(NULL)
    {
      init_linear_solver();
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(DiscreteProblem<Scalar>* dp, Hermes::MatrixSolverType matrix_solver_type) : NonlinearSolver<Scalar>(dp, matrix_solver_type), kept_jacobian(NULL)
    {
      init_linear_solver();
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::init_linear_solver()
    {
      // Set up the solver, jacobian, and residual according to the solver selection.
      jacobian = create_matrix<Scalar>(this->matrix_solver_type);
      residual = create_vector<Scalar>(this->matrix_solver_type);
      linear_solver = create_linear_solver<Scalar>(this->matrix_solver_type, jacobian, residual);
      reset_times();
      this->timer = NULL;
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::~NewtonSolver()
    {
      if(kept_jacobian != NULL)
        delete kept_jacobian;
      delete jacobian;
      delete residual;
      delete linear_solver;
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      return solve(coeff_vec, 1E-8, 100, false);
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::solve(Scalar* coeff_vec, bool residual_as_function)
    {
      return solve(coeff_vec, 1E-8, 100, residual_as_function);
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::solve(Scalar* coeff_vec, double newton_tol, int newton_max_iter, bool residual_as_function)
    {      
      // Delete the old solution vector, if there is any.
      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }

      // Obtain the number of degrees of freedom.
      int ndof = this->dp->get_num_dofs();

      // The Newton's loop.
      double residual_norm;
      int it = 1;
      
      bool delete_timer = false;
      if (this->timer == NULL)
      {
        this->timer = new TimePeriod;
        delete_timer = true;
      }
            
      this->timer->tick();
      setup_time += this->timer->last();
      
      while (1)
      {        
        // Assemble the residual vector.
        this->dp->assemble(coeff_vec, residual);
        
        this->timer->tick();
        assemble_time += this->timer->last();

        // Measure the residual norm.
        if (residual_as_function)
        {
          // Prepare solutions for measuring residual norm.
          Hermes::vector<Solution<Scalar>*> solutions;
          Hermes::vector<bool> dir_lift_false;
          for (unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++) {
            solutions.push_back(new Solution<Scalar>());
            dir_lift_false.push_back(false);
          }
          Solution<Scalar>::vector_to_solutions(residual, static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces(), solutions, dir_lift_false);

          // Calculate the norm.
          residual_norm = Global<Scalar>::calc_norms(solutions);

          // Clean up.
          for (unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++)
            delete solutions[i];
        }
        else 
        {
          // Calculate the l2-norm of residual vector, this is the traditional way.
          residual_norm = Global<Scalar>::get_l2_norm(residual);
        }

        // Info for the user.
        if(it == 1) {
          if(this->verbose_output)
            info("---- Newton initial residual norm: %g", residual_norm);
        }
        else 
          if(this->verbose_output)
            info("---- Newton iter %d, residual norm: %g", it - 1, residual_norm);

        // If maximum allowed residual norm is exceeded, fail.
        if (residual_norm > max_allowed_residual_norm)
        {
          if (this->verbose_output)
          {
            info("Current residual norm: %g", residual_norm);
            info("Maximum allowed residual norm: %g", max_allowed_residual_norm);
            info("Newton solve not successful, returning false.");
          }
          break;
        }

        // If residual norm is within tolerance, return 'true'.
        // This is the only correct way of ending.
        if (residual_norm < newton_tol && it > 1) {
          // We want to return the solution in a different structure.
          this->sln_vector = new Scalar[ndof];
          for (int i = 0; i < ndof; i++)
            this->sln_vector[i] = coeff_vec[i];

          this->timer->tick();
          solve_time += this->timer->last();
          
          if (delete_timer)
          {
            delete this->timer;
            this->timer = NULL;
          }
          
          return true;
        }
        
        this->timer->tick();
        solve_time += this->timer->last();

        // Assemble the jacobian.
        this->dp->assemble(coeff_vec, jacobian);
        
        this->timer->tick();
        assemble_time += this->timer->last();

        // Multiply the residual vector with -1 since the matrix
        // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
        residual->change_sign();
        
        // Solve the linear system.
        if(!linear_solver->solve()) {
          if (this->verbose_output) 
            info ("Matrix<Scalar> solver failed. Returning false.\n");
          break;
        }
        
        // Add \deltaY^{n+1} to Y^n.
        for (int i = 0; i < ndof; i++)
          coeff_vec[i] += linear_solver->get_sln_vector()[i];

        // Increase the number of iterations and test if we are still under the limit.
        if (it++ >= newton_max_iter)
        {
          if (this->verbose_output) 
            info("Maximum allowed number of Newton iterations exceeded, returning false.");
          break;
        }
        
        this->timer->tick();
        solve_time += this->timer->last();
      }
      // Return false.
      // All 'bad' situations end here. 
      
      if (delete_timer)
      {
        delete this->timer;
        this->timer = NULL;
      }
      
      return false;
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::solve_keep_jacobian(Scalar* coeff_vec, bool residual_as_function)
    {
      return solve(coeff_vec, 1E-8, 100, residual_as_function);
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::solve_keep_jacobian(Scalar* coeff_vec, double newton_tol, int newton_max_iter, bool residual_as_function)
    {
      // Obtain the number of degrees of freedom.
      int ndof = this->dp->get_num_dofs();

      // The Newton's loop.
      double residual_norm;
      int it = 1;
      while (1)
      {
        // Assemble the residual vector.
        this->dp->assemble(coeff_vec, residual);

        // Measure the residual norm.
        if (residual_as_function)
        {
          // Prepare solutions for measuring residual norm.
          Hermes::vector<Solution<Scalar>*> solutions;
          Hermes::vector<bool> dir_lift_false;
          for (unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++) {
            solutions.push_back(new Solution<Scalar>());
            dir_lift_false.push_back(false);
          }
          Solution<Scalar>::vector_to_solutions(residual, static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces(), solutions, dir_lift_false);

          // Calculate the norm.
          residual_norm = Global<Scalar>::calc_norms(solutions);

          // Clean up.
          for (unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++)
            delete solutions[i];
        }
        else 
        {
          // Calculate the l2-norm of residual vector, this is the traditional way.
          residual_norm = Global<Scalar>::get_l2_norm(residual);
        }

        // Info for the user.
        if(it == 1) 
        {
          if(this->verbose_output)
            info("---- Newton initial residual norm: %g", residual_norm);
        }
        else 
          if(this->verbose_output)
            info("---- Newton iter %d, residual norm: %g", it - 1, residual_norm);

        // If maximum allowed residual norm is exceeded, fail.
        if (residual_norm > max_allowed_residual_norm)
        {
          if (this->verbose_output)
          {
            info("Current residual norm: %g", residual_norm);
            info("Maximum allowed residual norm: %g", max_allowed_residual_norm);
            info("Newton solve not successful, returning false.");
          }
          break;
        }

        // If residual norm is within tolerance, return 'true'.
        // This is the only correct way of ending.
        if (residual_norm < newton_tol && it > 1) {
          // We want to return the solution in a different structure.
          this->sln_vector = new Scalar[ndof];
          for (int i = 0; i < ndof; i++)
            this->sln_vector[i] = coeff_vec[i];

          return true;
        }

        // Assemble and keep the jacobian if this has not been done before.
        if(kept_jacobian == NULL) {
          kept_jacobian = create_matrix<Scalar>(this->matrix_solver_type);
          this->dp->assemble(coeff_vec, kept_jacobian);
        }

        // Multiply the residual vector with -1 since the matrix
        // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
        residual->change_sign();

        // Solve the linear system.
        if(!linear_solver->solve()) {
          if (this->verbose_output) 
            info ("Matrix<Scalar> solver failed. Returning false.\n");
          break;
        }

        // Add \deltaY^{n+1} to Y^n.
        for (int i = 0; i < ndof; i++)
          coeff_vec[i] += linear_solver->get_sln_vector()[i];

        // Increase the number of iterations and test if we are still under the limit.
        if (it++ >= newton_max_iter)
        {
          if (this->verbose_output) 
            info("Maximum allowed number of Newton iterations exceeded, returning false.");
          break;
        }
      }

      // Return false.
      // All 'bad' situations end here.
      return false;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_iterative_method(const char* iterative_method_name)
    {
      NonlinearSolver<Scalar>::set_iterative_method(iterative_method_name);
      // Set iterative method and preconditioner in case of iterative solver AztecOO.
#ifdef HAVE_AZTECOO
      if(this->matrix_solver_type != SOLVER_AZTECOO)
      {
        warning("Trying to set iterative method for a different solver than AztecOO.");
        return;
      }
      if (this->matrix_solver_type == SOLVER_AZTECOO)
        dynamic_cast<Hermes::Solvers::AztecOOSolver<Scalar>*>(linear_solver)->set_solver(iterative_method_name);
#else
      warning("Trying to set iterative method without AztecOO present.");
#endif
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_preconditioner(const char* preconditioner_name)
    {
      NonlinearSolver<Scalar>::set_preconditioner(preconditioner_name);
      // Set iterative method and preconditioner in case of iterative solver AztecOO.
#ifdef HAVE_AZTECOO
      if(this->matrix_solver_type != SOLVER_AZTECOO)
      {
        warning("Trying to set iterative method for a different solver than AztecOO.");
        return;
      }
      if (this->matrix_solver_type == SOLVER_AZTECOO)
        dynamic_cast<Hermes::Solvers::AztecOOSolver<Scalar> *>(linear_solver)->set_precond(preconditioner_name);
#else
      warning("Trying to set iterative method without AztecOO present.");
#endif
    }

    template class HERMES_API NewtonSolver<double>;
    template class HERMES_API NewtonSolver<std::complex<double> >;
  }
}
