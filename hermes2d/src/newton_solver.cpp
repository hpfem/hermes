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
#include "projections/ogprojection.h"
#include "hermes_common.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver() : NonlinearSolver<Scalar>(new DiscreteProblem<Scalar>()), own_dp(true), kept_jacobian(NULL)
    {
      init_attributes();
      init_linear_solver();
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(DiscreteProblem<Scalar>* dp) : NonlinearSolver<Scalar>(dp), own_dp(false), kept_jacobian(NULL)
    {
      init_attributes();
      init_linear_solver();
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(const WeakForm<Scalar>* wf, const Space<Scalar>* space) : NonlinearSolver<Scalar>(new DiscreteProblem<Scalar>(wf, space)), own_dp(true), kept_jacobian(NULL)
    {
      init_attributes();
      init_linear_solver();
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar> *> spaces) : NonlinearSolver<Scalar>(new DiscreteProblem<Scalar>(wf, spaces)), own_dp(true), kept_jacobian(NULL)
    {
      init_attributes();
      init_linear_solver();
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::isOkay() const
    {
      if(static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_weak_formulation() == NULL)
        return false;
      if(static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size() == 0)
        return false;
      return true;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::init_attributes()
    {
      this->newton_tol = 1e-8;
      this->newton_max_iter = 15;
      this->residual_as_function = false;
      this->max_allowed_residual_norm = 1E9;
      this->min_allowed_damping_coeff = 1E-4;
      this->currentDampingCofficient = 1.0;
      this->manual_damping = false;
      this->auto_damping_ratio = 2.0;
      this->initial_auto_damping_ratio = 1.0;
      this->sufficient_improvement_factor = 0.95;
      this->necessary_successful_steps_to_increase = 1;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_newton_tol(double newton_tol)
    {
      this->newton_tol = newton_tol;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_weak_formulation(const WeakForm<Scalar>* wf)
    {
      (static_cast<DiscreteProblem<Scalar>*>(this->dp))->set_weak_formulation(wf);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_newton_max_iter(int newton_max_iter)
    {
      if(newton_max_iter < 1)
        throw Exceptions::ValueException("newton_max_iter", newton_max_iter, 1);
      this->newton_max_iter = newton_max_iter;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set)
    {
      if(max_allowed_residual_norm_to_set <= 0.0)
        throw Exceptions::ValueException("max_allowed_residual_norm_to_set", max_allowed_residual_norm_to_set, 0.0);
      this->max_allowed_residual_norm = max_allowed_residual_norm_to_set;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_min_allowed_damping_coeff(double min_allowed_damping_coeff_to_set)
    {
      if(min_allowed_damping_coeff_to_set <= 0.0)
        throw Exceptions::ValueException("min_allowed_damping_coeff_to_set", min_allowed_damping_coeff_to_set, 0.0);
      this->min_allowed_damping_coeff = min_allowed_damping_coeff_to_set;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_residual_as_function()
    {
      this->residual_as_function = true;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_time(double time)
    {
      Hermes::vector<Space<Scalar>*> spaces;
      for(unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++)
        spaces.push_back(const_cast<Space<Scalar>*>(static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(i)));

      Space<Scalar>::update_essential_bc_values(spaces, time);
      const_cast<WeakForm<Scalar>*>(static_cast<DiscreteProblem<Scalar>*>(this->dp)->wf)->set_current_time(time);
    }
      
    template<typename Scalar>
    void NewtonSolver<Scalar>::set_time_step(double time_step)
    {
      const_cast<WeakForm<Scalar>*>(static_cast<DiscreteProblem<Scalar>*>(this->dp)->wf)->set_current_time_step(time_step);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_spaces(Hermes::vector<const Space<Scalar>*> spaces)
    {
      static_cast<DiscreteProblem<Scalar>*>(this->dp)->set_spaces(spaces);
      if(kept_jacobian != NULL)
        delete kept_jacobian;
      kept_jacobian = NULL;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_space(const Space<Scalar>* space)
    {
      static_cast<DiscreteProblem<Scalar>*>(this->dp)->set_space(space);
      if(kept_jacobian != NULL)
        delete kept_jacobian;
      kept_jacobian = NULL;
    }
    
    template<typename Scalar>
    Hermes::vector<const Space<Scalar>*> NewtonSolver<Scalar>::get_spaces() const
    {
      return static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces();
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::init_linear_solver()
    {
      // Set up the solver, jacobian, and residual according to the solver selection.
      jacobian = create_matrix<Scalar>();
      residual = create_vector<Scalar>();
      linear_solver = create_linear_solver<Scalar>(jacobian, residual);
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::~NewtonSolver()
    {
      if(kept_jacobian != NULL)
        delete kept_jacobian;
      delete jacobian;
      delete residual;
      delete linear_solver;
      if(own_dp)
        delete this->dp;
      else
        static_cast<DiscreteProblem<Scalar>*>(this->dp)->have_matrix = false;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_manual_damping_coeff(bool onOff, double coeff)
    {
      if(coeff <= 0.0 || coeff > 1.0)
        throw Exceptions::ValueException("coeff", coeff, 0.0, 1.0);
      if(onOff)
      {
        this->manual_damping = true;
        this->currentDampingCofficient = coeff;
      }
      else
        this->manual_damping = false;
    }
      
    template<typename Scalar>
    void NewtonSolver<Scalar>::set_initial_auto_damping_coeff(double coeff)
    {
      if(coeff <= 0.0 || coeff > 1.0)
        throw Exceptions::ValueException("coeff", coeff, 0.0, 1.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      else
        this->currentDampingCofficient = coeff;
    }
      
    template<typename Scalar>
    void NewtonSolver<Scalar>::set_auto_damping_ratio(double ratio)
    {
      if(ratio <= 0.0)
        throw Exceptions::ValueException("ratio", ratio, 0.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
        this->auto_damping_ratio = ratio;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_sufficient_improvement_factor(double ratio)
    {
      if(ratio <= 0.0)
        throw Exceptions::ValueException("ratio", ratio, 0.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      this->sufficient_improvement_factor = ratio;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_necessary_successful_steps_to_increase(unsigned int steps)
    {
      if(steps < 1)
        throw Exceptions::ValueException("steps", steps, 0.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      this->necessary_successful_steps_to_increase = steps;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve(Solution<Scalar>* initial_guess)
    {
      Hermes::vector<Solution<Scalar>*> vectorToPass;
      vectorToPass.push_back(initial_guess);
      this->solve(vectorToPass);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve(Hermes::vector<Solution<Scalar>*> initial_guess)
    {
      int ndof = this->dp->get_num_dofs();
      Scalar* coeff_vec = new Scalar[ndof];
      OGProjection<Scalar> ogProjection;
      ogProjection.project_global(static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces(), initial_guess, coeff_vec);
      this->solve(coeff_vec);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      this->check();

      this->tick();

      // Obtain the number of degrees of freedom.
      int ndof = this->dp->get_num_dofs();

      // If coeff_vec == NULL then create a zero vector.
      bool delete_coeff_vec = false;
      if(coeff_vec == NULL)
      {
        coeff_vec = new Scalar[ndof];
        memset(coeff_vec, 0, ndof*sizeof(Scalar));
        delete_coeff_vec = true;
      }

      // Vector for restarting steps.
      Scalar* coeff_vec_back = new Scalar[ndof];
      memset(coeff_vec_back, 0, ndof*sizeof(Scalar));

      // Delete the old solution vector, if there is any.
      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }

      // The Newton's loop.
      double residual_norm;
      double last_residual_norm;
      int it = 1;
      int successfulSteps = 0;

      this->on_initialization();

      while (true)
      {
        this->on_step_begin();

        // Assemble just the residual vector.
        this->dp->assemble(coeff_vec, residual);
        if(this->output_rhsOn && (this->output_rhsIterations == -1 || this->output_rhsIterations >= it))
        {
          char* fileName = new char[this->RhsFilename.length() + 5];
          if(this->RhsFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
            sprintf(fileName, "%s%i.m", this->RhsFilename.c_str(), it);
          else
            sprintf(fileName, "%s%i", this->RhsFilename.c_str(), it);
          FILE* rhs_file = fopen(fileName, "w+");
          residual->dump(rhs_file, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
          fclose(rhs_file);
          delete [] fileName;
        }
        
        Element* e;
        for(unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++)
          for_all_active_elements(e, static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces()[i]->get_mesh())
            static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces()[i]->edata[e->id].changed_in_last_adaptation = false;

        // Measure the residual norm.
        if(residual_as_function)
        {
          // Prepare solutions for measuring residual norm.
          Hermes::vector<Solution<Scalar>*> solutions;
          Hermes::vector<bool> dir_lift_false;
          for (unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++) 
          {
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
          if(it == 1)
            last_residual_norm = residual_norm = Global<Scalar>::get_l2_norm(residual);
          else
          {
            last_residual_norm = residual_norm;
            residual_norm = Global<Scalar>::get_l2_norm(residual);
          }
        }

        // Info for the user.
        if(it == 1)
          this->info("\tNewton: initial residual norm: %g", residual_norm);
        else
        {
          this->info("\tNewton: iteration %d, residual norm: %g", it - 1, residual_norm);
          if(!this->manual_damping && !((residual_norm > max_allowed_residual_norm) || (residual_norm < newton_tol && it > 1)))
          {
            if(residual_norm < last_residual_norm * this->sufficient_improvement_factor)
            {
              if(++successfulSteps >= this->necessary_successful_steps_to_increase)
                this->currentDampingCofficient = std::min(this->initial_auto_damping_ratio, 2 * this->currentDampingCofficient);
              if(residual_norm < last_residual_norm)
                this->info("\t Newton: step successful, calculation continues with damping coefficient: %g.", this->currentDampingCofficient);
            }
            else
            {
              successfulSteps = 0;
              if(this->currentDampingCofficient < this->min_allowed_damping_coeff)
              {
                this->warn("\t Newton: results NOT improved, current damping coefficient is at the minimum possible level: %g.", min_allowed_damping_coeff);
                this->info("\t  If you want to decrease the minimum level, use the method set_min_allowed_damping_coeff()");
                if(delete_coeff_vec)
                {
                  delete [] coeff_vec;
                  coeff_vec = NULL;
                }

                delete [] coeff_vec_back;

                throw Exceptions::Exception("Newton NOT converged because of damping coefficient could not be decreased anymore to possibly handle non-converging process.");
              }
              else
              {
                this->currentDampingCofficient = (1 / this->auto_damping_ratio) * this->currentDampingCofficient;
                this->warn(" Newton: results NOT improved, step restarted with damping coefficient: %g.", this->currentDampingCofficient);
              }
            }
          }
        }

        // If maximum allowed residual norm is exceeded, fail.
        if(residual_norm > max_allowed_residual_norm)
        {
          this->tick();
          this->info("\tNewton: solution duration: %f s.\n", this->last());
          this->on_finish();
          
          if(delete_coeff_vec)
          {
            delete [] coeff_vec;
            coeff_vec = NULL;
          }

          delete [] coeff_vec_back;
          
          throw Exceptions::ValueException("residual norm", residual_norm, max_allowed_residual_norm);
        }

        // If residual norm is within tolerance, return 'true'.
        // This is the only correct way of ending.
        if(residual_norm < newton_tol && it > 1)
        {
          // We want to return the solution in a different structure.
          this->sln_vector = new Scalar[ndof];
          for (int i = 0; i < ndof; i++)
            this->sln_vector[i] = coeff_vec[i];

          if(delete_coeff_vec)
          {
            delete [] coeff_vec;
            coeff_vec = NULL;
          }

          delete [] coeff_vec_back;

          this->on_finish();

          this->tick();
          this->info("\tNewton: solution duration: %f s.\n", this->last());

          return;
        }

        // Assemble just the jacobian.
        this->dp->assemble(coeff_vec, jacobian);
        if(this->output_matrixOn && (this->output_matrixIterations == -1 || this->output_matrixIterations >= it))
        {
          char* fileName = new char[this->matrixFilename.length() + 5];
          if(this->matrixFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
            sprintf(fileName, "%s%i.m", this->matrixFilename.c_str(), it);
          else
            sprintf(fileName, "%s%i", this->matrixFilename.c_str(), it);
          FILE* matrix_file = fopen(fileName, "w+");

          jacobian->dump(matrix_file, this->matrixVarname.c_str(), this->matrixFormat, this->matrix_number_format);
          fclose(matrix_file);
          delete [] fileName;
        }

        this->on_step_end();

        // Multiply the residual vector with -1 since the matrix
        // equation reads J(Y^n) \deltaY^{n + 1} = -F(Y^n).
        residual->change_sign();

        // Solve the linear system.
        if(!linear_solver->solve())
          throw Exceptions::LinearMatrixSolverException();

        // Add \deltaY^{n + 1} to Y^n.
        // The good case.
        if(residual_norm < last_residual_norm * this->sufficient_improvement_factor || this->manual_damping || it == 1)
        {
          memcpy(coeff_vec_back, coeff_vec, sizeof(Scalar)*ndof);
          for (int i = 0; i < ndof; i++)
            coeff_vec[i] += currentDampingCofficient * linear_solver->get_sln_vector()[i];
        }
        else
        {
          for (int i = 0; i < ndof; i++)
            coeff_vec[i] = coeff_vec_back[i] + currentDampingCofficient * (coeff_vec[i] - coeff_vec_back[i]);
        }

        // Increase the number of iterations and test if we are still under the limit.
        if(it++ >= newton_max_iter)
        {
          // We want to return the solution in a different structure.
          this->sln_vector = new Scalar[ndof];
          for (int i = 0; i < ndof; i++)
            this->sln_vector[i] = coeff_vec[i];

          if(delete_coeff_vec)
          {
            delete [] coeff_vec;
            coeff_vec = NULL;
          }

          delete [] coeff_vec_back;
          coeff_vec_back = NULL;

          this->tick();
          this->info("\tNewton: solution duration: %f s.\n", this->last());

          this->on_finish();
          throw Exceptions::ValueException("iterations", it, newton_max_iter);
        }
      }
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve_keep_jacobian(Solution<Scalar>* initial_guess)
    {
      Hermes::vector<Solution<Scalar>*> vectorToPass;
      vectorToPass.push_back(initial_guess);
      this->solve_keep_jacobian(vectorToPass);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve_keep_jacobian(Hermes::vector<Solution<Scalar>*> initial_guess)
    {
      int ndof = this->dp->get_num_dofs();
      Scalar* coeff_vec = new Scalar[ndof];
      OGProjection<Scalar> ogProjection;
      ogProjection.project_global(static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces(), initial_guess, coeff_vec);
      this->solve_keep_jacobian(coeff_vec);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve_keep_jacobian(Scalar* coeff_vec)
    {
      this->check();

      this->tick();

      // Obtain the number of degrees of freedom.
      int ndof = this->dp->get_num_dofs();

      // If coeff_vec == NULL then create a zero vector.
      bool delete_coeff_vec = false;
      if(coeff_vec == NULL)
      {
        coeff_vec = new Scalar[ndof];
        memset(coeff_vec, 0, ndof*sizeof(Scalar));
        delete_coeff_vec = true;
      }

      // Vector for restarting steps.
      Scalar* coeff_vec_back = new Scalar[ndof];
      memset(coeff_vec_back, 0, ndof*sizeof(Scalar));

      
      // The Newton's loop.
      double residual_norm;
      double last_residual_norm;
      int it = 1;
      int successfulSteps = 0;

      this->on_initialization();

      while (true)
      {
        this->on_step_begin();

        // Assemble the residual vector.
        this->dp->assemble(coeff_vec, residual);
        if(this->output_rhsOn && (this->output_rhsIterations == -1 || this->output_rhsIterations >= it))
        {
          char* fileName = new char[this->RhsFilename.length() + 5];
          if(this->RhsFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
            sprintf(fileName, "%s%i.m", this->RhsFilename.c_str(), it);
          else
            sprintf(fileName, "%s%i", this->RhsFilename.c_str(), it);
          FILE* rhs_file = fopen(fileName, "w+");
          residual->dump(rhs_file, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
          fclose(rhs_file);
        }

        Element* e;
        for(unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++)
          for_all_active_elements(e, static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces()[i]->get_mesh())
            static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces()[i]->edata[e->id].changed_in_last_adaptation = false;

        // Measure the residual norm.
        if(residual_as_function)
        {
          // Prepare solutions for measuring residual norm.
          Hermes::vector<Solution<Scalar>* > solutions;
          Hermes::vector<bool> dir_lift_false;
          for (unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++) {
            solutions.push_back(new Solution<Scalar>());
            dir_lift_false.push_back(false);
          }
          Solution<Scalar>::vector_to_solutions(residual,
            static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces(), solutions, dir_lift_false);

          // Calculate the norm.
          residual_norm = Global<Scalar>::calc_norms(solutions);

          // Clean up.
          for (unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++)
            delete solutions[i];
        }
        else
        {
          // Calculate the l2-norm of residual vector, this is the traditional way.
          if(it == 1)
            last_residual_norm = residual_norm = Global<Scalar>::get_l2_norm(residual);
          else
          {
            last_residual_norm = residual_norm;
            residual_norm = Global<Scalar>::get_l2_norm(residual);
          }
        }

        // Info for the user.
        if(it == 1)
          this->info("\tNewton: initial residual norm: %g", residual_norm);
        else
        {
          this->info("\tNewton: iteration %d, residual norm: %g", it - 1, residual_norm);
          if(!this->manual_damping && !((residual_norm > max_allowed_residual_norm) || (residual_norm < newton_tol && it > 1)))
          {
            if(residual_norm < last_residual_norm * this->sufficient_improvement_factor)
            {
              if(++successfulSteps >= this->necessary_successful_steps_to_increase)
                this->currentDampingCofficient = std::min(this->initial_auto_damping_ratio, 2 * this->currentDampingCofficient);
              if(residual_norm < last_residual_norm)
                this->info("\t Newton: step successful, calculation continues with damping coefficient: %g.", this->currentDampingCofficient);
            }
            else
            {
              successfulSteps = 0;
              if(this->currentDampingCofficient < this->min_allowed_damping_coeff)
              {
                this->warn("\t Newton: results NOT improved, current damping coefficient is at the minimum possible level: %g.", min_allowed_damping_coeff);
                this->info("\t  If you want to decrease the minimum level, use the method set_min_allowed_damping_coeff()");
                
                if(delete_coeff_vec)
                {
                  delete [] coeff_vec;
                  coeff_vec = NULL;
                }

                delete [] coeff_vec_back;
          
                throw Exceptions::Exception("Newton NOT converged because of damping coefficient could not be decreased anymore to possibly handle non-converging process.");
              }
              else
              {
                this->currentDampingCofficient = (1 / this->auto_damping_ratio) * this->currentDampingCofficient;
                this->warn("\t Newton: results NOT improved, step restarted with damping coefficient: %g.", this->currentDampingCofficient);
              }
            }
          }
        }

        // If maximum allowed residual norm is exceeded, fail.
        if(residual_norm > max_allowed_residual_norm)
        {
          this->tick();
          this->info("\tNewton: solution duration: %f s.", this->last());

          if(delete_coeff_vec)
          {
            delete [] coeff_vec;
            coeff_vec = NULL;
          }

          delete [] coeff_vec_back;

          this->on_finish();

          throw Exceptions::ValueException("residual norm", residual_norm, max_allowed_residual_norm);
        }

        // If residual norm is within tolerance, return 'true'.
        // This is the only correct way of ending.
        if(residual_norm < newton_tol && it > 1) 
        {
          // We want to return the solution in a different structure.
          this->sln_vector = new Scalar[ndof];
          for (int i = 0; i < ndof; i++)
            this->sln_vector[i] = coeff_vec[i];

          if(delete_coeff_vec)
          {
            delete [] coeff_vec;
            coeff_vec = NULL;
          }

          delete [] coeff_vec_back;

          this->tick();
          this->info("\tNewton: solution duration: %f s.", this->last());

          this->on_finish();
          return;
        }

        // Assemble and keep the jacobian if this has not been done before.
        // Also declare that LU-factorization in case of a direct solver will be done only once and reused afterwards.
        if(kept_jacobian == NULL || !(static_cast<DiscreteProblem<Scalar>*>(this->dp))->have_matrix) 
        {
          if(kept_jacobian != NULL)
            delete kept_jacobian;

          kept_jacobian = create_matrix<Scalar>();

          // Give the matrix solver the correct Jacobian. NOTE: It would be cleaner if the whole decision whether to keep
          // Jacobian or not was made in the constructor.
          //
          // Delete the matrix solver created in the constructor.
          delete linear_solver;
          // Create new matrix solver with correct matrix.
          linear_solver = create_linear_solver<Scalar>(kept_jacobian, residual);

          this->dp->assemble(coeff_vec, kept_jacobian);

          if(this->output_matrixOn && (this->output_matrixIterations == -1 || this->output_matrixIterations >= it))
          {
            char* fileName = new char[this->matrixFilename.length() + 5];
            if(this->matrixFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
              sprintf(fileName, "%s%i.m", this->matrixFilename.c_str(), it);
            else
              sprintf(fileName, "%s%i", this->matrixFilename.c_str(), it);
            FILE* matrix_file = fopen(fileName, "w+");

            kept_jacobian->dump(matrix_file, this->matrixVarname.c_str(), this->matrixFormat, this->matrix_number_format);
            fclose(matrix_file);
          }

          linear_solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
        }

        this->on_step_end();

        // Multiply the residual vector with -1 since the matrix
        // equation reads J(Y^n) \deltaY^{n + 1} = -F(Y^n).
        residual->change_sign();

        // Solve the linear system.
        if(!linear_solver->solve()) 
        {
          throw Exceptions::LinearMatrixSolverException();
        }

         // Add \deltaY^{n + 1} to Y^n.
        // The good case.
        if(residual_norm < last_residual_norm * this->sufficient_improvement_factor || this->manual_damping || it == 1)
        {
          memcpy(coeff_vec_back, coeff_vec, sizeof(Scalar)*ndof);
          for (int i = 0; i < ndof; i++)
            coeff_vec[i] += currentDampingCofficient * linear_solver->get_sln_vector()[i];
        }
        else
        {
          for (int i = 0; i < ndof; i++)
            coeff_vec[i] = coeff_vec_back[i] + currentDampingCofficient * (coeff_vec[i] - coeff_vec_back[i]);
        }

        // Increase the number of iterations and test if we are still under the limit.
        if(it++ >= newton_max_iter)
        {
          // We want to return the solution in a different structure.
          this->sln_vector = new Scalar[ndof];
          for (int i = 0; i < ndof; i++)
            this->sln_vector[i] = coeff_vec[i];

          if(delete_coeff_vec)
          {
            delete [] coeff_vec;
            coeff_vec = NULL;
          }

          delete [] coeff_vec_back;

          this->tick();
          this->info("\tNewton: solution duration: %f s.\n", this->last());

          this->on_finish();
          throw Exceptions::ValueException("iterations", it, newton_max_iter);
        }
      }
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_iterative_method(const char* iterative_method_name)
    {
      NonlinearSolver<Scalar>::set_iterative_method(iterative_method_name);
      // Set iterative method in case of iterative solver AztecOO.
#ifdef HAVE_AZTECOO
      dynamic_cast<Hermes::Solvers::AztecOOSolver<Scalar>*>(linear_solver)->set_solver(iterative_method_name);
#else
      this->warn("Trying to set iterative method without AztecOO present.");
#endif
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_preconditioner(const char* preconditioner_name)
    {
      NonlinearSolver<Scalar>::set_preconditioner(preconditioner_name);
      // Set preconditioner in case of iterative solver AztecOO.
#ifdef HAVE_AZTECOO
      dynamic_cast<Hermes::Solvers::AztecOOSolver<Scalar> *>(linear_solver)->set_precond(preconditioner_name);
#else
      this->warn("Trying to set iterative method without AztecOO present.");
#endif
    }

    template class HERMES_API NewtonSolver<double>;
    template class HERMES_API NewtonSolver<std::complex<double> >;
  }
}
