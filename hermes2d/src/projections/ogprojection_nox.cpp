// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.
/*! \file ogprojection_nox.cpp
\brief Orthogonal projection via NOX (matrix-free).
*/
#include "projections/ogprojection_nox.h"
#include "space.h"
#include "discrete_problem.h"
#include "newton_solver_nox.h"
#if(defined HAVE_NOX && defined HAVE_EPETRA && defined HAVE_TEUCHOS)

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    OGProjectionNOX<Scalar>::OGProjectionNOX() : ndof(0)
    {
    }

    template<typename Scalar>
    void OGProjectionNOX<Scalar>::project_internal(const Space<Scalar>* space, WeakForm<Scalar>* wf,
      Scalar* target_vec, double newton_tol, int newton_max_iter)
    {
        // Sanity check.
        if(space == NULL)
          throw Hermes::Exceptions::Exception("this->space == NULL in project_internal().");

      // Get dimension of the space.
      int ndof = space->get_num_dofs();

      // Initialize DiscreteProblem.
      DiscreteProblem<Scalar> dp(wf, space);

      // Initial coefficient vector for the Newton's method.
      Scalar* coeff_vec = new Scalar[ndof];
      memset(coeff_vec, 0, ndof*sizeof(Scalar));

      const char* iterative_method = "GMRES";           // Name of the iterative method employed by AztecOO (ignored
      // by the other solvers).
      // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
      const char* preconditioner = "New Ifpack";           // Name of the preconditioner employed by AztecOO
      // Possibilities: None" - No preconditioning.
      // "AztecOO" - AztecOO internal preconditioner.
      // "New Ifpack" - Ifpack internal preconditioner.
      // "ML" - Multi level preconditione
      unsigned message_type = NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration | NOX::Utils::InnerIteration | NOX::Utils::Parameters | NOX::Utils::LinearSolverDetails;
      // NOX error messages, see NOX_Utils.h.
      double ls_tolerance = 1e-5;                       // Tolerance for linear system.
      unsigned flag_absresid = 0;                       // Flag for absolute value of the residuum.
      double abs_resid = 1.0e-3;                        // Tolerance for absolute value of the residuum.
      unsigned flag_relresid = 1;                       // Flag for relative value of the residuum.
      double rel_resid = newton_tol;                        // Tolerance for relative value of the residuum.
      int max_iters = newton_max_iter;                  // Max number of iterations.

      // Initialize NOX.
      NewtonSolverNOX<Scalar> newton_nox(&dp);

      // Set NOX parameters.
      newton_nox.set_verbose_output(false);
      newton_nox.set_output_flags(message_type);
      newton_nox.set_ls_type(iterative_method);
      newton_nox.set_ls_tolerance(ls_tolerance);
      newton_nox.set_conv_iters(max_iters);
      if(flag_absresid)
        newton_nox.set_conv_abs_resid(abs_resid);
      if(flag_relresid)
        newton_nox.set_conv_rel_resid(rel_resid);
      newton_nox.set_precond(preconditioner);
      newton_nox.set_precond_reuse("Rebuild");

      // Perform Newton's iteration via NOX
      newton_nox.solve(coeff_vec);

      delete [] coeff_vec;

      if(target_vec != NULL)
        for (int i = 0; i < ndof; i++)
          target_vec[i] = newton_nox.get_sln_vector()[i];
    }

    template<typename Scalar>
    void OGProjectionNOX<Scalar>::project_global(const Space<Scalar>* space,
      MatrixFormVol<Scalar>* custom_projection_jacobian,
      VectorFormVol<Scalar>* custom_projection_residual,
      Scalar* target_vec,
      double newton_tol, int newton_max_iter)
    {
      // Define projection weak form.
      WeakForm<Scalar>* proj_wf = new WeakForm<Scalar>(1);
      proj_wf->add_matrix_form(custom_projection_jacobian);
      proj_wf->add_vector_form(custom_projection_residual);

      // Call the main function.
      project_internal(space, proj_wf, target_vec, newton_tol, newton_max_iter);

      // Clean up.
      delete proj_wf;
    }

    template<typename Scalar>
    void OGProjectionNOX<Scalar>::project_global(const Space<Scalar>* space,
      MeshFunction<Scalar>* source_meshfn, Scalar* target_vec,
      ProjNormType proj_norm,
      double newton_tol, int newton_max_iter)
    {
      bool sln_vector_loaded = true;

      if(dynamic_cast<Solution<Scalar>*>(source_meshfn) != NULL && dynamic_cast<Solution<Scalar>*>(source_meshfn)->get_type() == HERMES_SLN)
        if(dynamic_cast<Solution<Scalar>*>(source_meshfn)->get_space() != NULL)
          if(dynamic_cast<Solution<Scalar>*>(source_meshfn)->get_space_seq() == space->get_seq()
            && dynamic_cast<Solution<Scalar>*>(source_meshfn)->get_sln_vector() != NULL)
            for(int j = 0; j < space->get_num_dofs(); j++)
              target_vec[j] = dynamic_cast<Solution<Scalar>*>(source_meshfn)->get_sln_vector()[j];
          else
            sln_vector_loaded = false;
        else
          sln_vector_loaded = false;
      else
        sln_vector_loaded = false;

      if(sln_vector_loaded)
        return;

      // Sanity checks.
      if(target_vec == NULL) throw Exceptions::NullException(3);

      // If projection norm is not provided, set it
      // to match the type of the space.
      ProjNormType norm = HERMES_UNSET_NORM;
      if(proj_norm == HERMES_UNSET_NORM)
      {
        SpaceType space_type = space->get_type();
        switch (space_type)
        {
        case HERMES_H1_SPACE: norm = HERMES_H1_NORM; break;
        case HERMES_HCURL_SPACE: norm = HERMES_HCURL_NORM; break;
        case HERMES_HDIV_SPACE: norm = HERMES_HDIV_NORM; break;
        case HERMES_L2_SPACE: norm = HERMES_L2_NORM; break;
        default: throw Hermes::Exceptions::Exception("Unknown space type in OGProjectionNOX<Scalar>::project_global().");
        }
      }
      else norm = proj_norm;

      // Define temporary projection weak form.
      WeakForm<Scalar>* proj_wf = new WeakForm<Scalar>(1);
      // Add Jacobian.
      proj_wf->add_matrix_form(new ProjectionMatrixFormVol(0, 0, norm));
      // Add Residual.
      proj_wf->add_vector_form(new ProjectionVectorFormVol(0, source_meshfn, norm));

      // Call main function.
      project_internal(space, proj_wf, target_vec, newton_tol, newton_max_iter);

      // Clean up.
      delete proj_wf;
    }

    template<typename Scalar>
    void OGProjectionNOX<Scalar>::project_global(const Space<Scalar>* space,
      Solution<Scalar>* source_sln, Solution<Scalar>* target_sln,
      ProjNormType proj_norm,
      double newton_tol, int newton_max_iter)
    {
      if(proj_norm == HERMES_UNSET_NORM)
      {
        SpaceType space_type = space->get_type();
        switch (space_type)
        {
        case HERMES_H1_SPACE: proj_norm = HERMES_H1_NORM; break;
        case HERMES_HCURL_SPACE: proj_norm = HERMES_HCURL_NORM; break;
        case HERMES_HDIV_SPACE: proj_norm = HERMES_HDIV_NORM; break;
        case HERMES_L2_SPACE: proj_norm = HERMES_L2_NORM; break;
        default: throw Hermes::Exceptions::Exception("Unknown space type in OGProjectionNOX<Scalar>::project_global().");
        }
      }

      // Calculate the coefficient vector.
      int ndof = space->get_num_dofs();
      Scalar* target_vec = new Scalar[ndof];
      project_global(space, source_sln, target_vec, proj_norm, newton_tol, newton_max_iter);

      // Translate coefficient vector into a Solution.
      Solution<Scalar>::vector_to_solution(target_vec, space, target_sln);

      // Clean up.
      delete [] target_vec;
    }

    template<typename Scalar>
    void OGProjectionNOX<Scalar>::project_global(Hermes::vector<const Space<Scalar>*> spaces,
      Hermes::vector<MeshFunction<Scalar>*> source_meshfns,
      Scalar* target_vec, Hermes::vector<ProjNormType> proj_norms,
      double newton_tol, int newton_max_iter)
    {
      int n = spaces.size();

      // Sanity checks.
      if(n != source_meshfns.size()) throw Exceptions::LengthException(1, 2, n, source_meshfns.size());
      if(target_vec == NULL) throw Exceptions::NullException(3);
      if(!proj_norms.empty() && n != proj_norms.size()) throw Exceptions::LengthException(1, 5, n, proj_norms.size());

      int start_index = 0;
      for (int i = 0; i < n; i++)
      {
        if(proj_norms.empty())
          project_global(spaces[i], source_meshfns[i], target_vec + start_index, HERMES_UNSET_NORM, newton_tol, newton_max_iter);
        else
          project_global(spaces[i], source_meshfns[i], target_vec + start_index, proj_norms[i], newton_tol, newton_max_iter);
        spaces[i]->assign_dofs(start_index);
        start_index += spaces[i]->get_num_dofs();
      }
    }

    template<typename Scalar>
    void OGProjectionNOX<Scalar>::project_global(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> source_slns,
      Scalar* target_vec, Hermes::vector<ProjNormType> proj_norms,
      double newton_tol, int newton_max_iter)
    {
      int n = spaces.size();

      // Sanity checks.
      if(n != source_slns.size()) throw Exceptions::LengthException(1, 2, n, source_slns.size());
      if(target_vec == NULL) throw Exceptions::NullException(3);
      if(!proj_norms.empty() && n != proj_norms.size()) throw Exceptions::LengthException(1, 5, n, proj_norms.size());

      int start_index = 0;
      for (int i = 0; i < n; i++)
      {
        if(proj_norms.empty())
          project_global(spaces[i], source_slns[i], target_vec + start_index, HERMES_UNSET_NORM, newton_tol, newton_max_iter);
        else
          project_global(spaces[i], source_slns[i], target_vec + start_index, proj_norms[i], newton_tol, newton_max_iter);
        start_index += spaces[i]->get_num_dofs();
      }
    }

    template<typename Scalar>
    void OGProjectionNOX<Scalar>::project_global(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> source_slns,
      Hermes::vector<Solution<Scalar>*> target_slns,
      Hermes::vector<ProjNormType> proj_norms, bool delete_old_meshes,
      double newton_tol, int newton_max_iter)
    {
      int n = spaces.size();

      // Sanity checks.
      if(n != source_slns.size()) throw Exceptions::LengthException(1, 2, n, source_slns.size());
      if(n != target_slns.size()) throw Exceptions::LengthException(1, 2, n, target_slns.size());
      if(!proj_norms.empty() && n != proj_norms.size()) throw Exceptions::LengthException(1, 5, n, proj_norms.size());

      int start_index = 0;
      for (int i = 0; i < n; i++)
      {
        if(proj_norms.empty())
          project_global(spaces[i], source_slns[i], target_slns[i], HERMES_UNSET_NORM, newton_tol, newton_max_iter);
        else
          project_global(spaces[i], source_slns[i], target_slns[i], proj_norms[i], newton_tol, newton_max_iter);
        start_index += spaces[i]->get_num_dofs();
      }
    }

    template class HERMES_API OGProjectionNOX<double>;
    // template class HERMES_API OGProjectionNOX<std::complex<double> >; //complex version of nox solver is not implemented
  }
}

#endif