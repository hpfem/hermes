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

#include "projections/ogprojection.h"
#include "space.h"
#include "linear_solver.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    OGProjection<Scalar>::OGProjection() : Hermes::Mixins::Loggable(), ndof(0)
    {
    }

    template<typename Scalar>
    void OGProjection<Scalar>::project_internal(const Space<Scalar>* space, WeakForm<Scalar>* wf,
  Scalar* target_vec, double newton_tol, int newton_max_iter)
    {
      // Sanity check.
      if(space == NULL)
        throw Hermes::Exceptions::Exception("this->space == NULL in project_internal().");

      // Initialize DiscreteProblem.
      DiscreteProblemLinear<Scalar> dp(wf, space);
      dp.setDoNotUseCache();

      // Initialize linear solver.
      Hermes::Hermes2D::LinearSolver<Scalar> linear_solver(&dp);

      // Perform Newton iteration.
      try
      {
        linear_solver.solve();
      }
      catch(Hermes::Exceptions::Exception e)
      {
        e.printMsg();
      }

      if(target_vec != NULL)
        for (int i = 0; i < space->get_num_dofs(); i++)
          target_vec[i] = linear_solver.get_sln_vector()[i];
    }

    template<typename Scalar>
    void OGProjection<Scalar>::project_global(const Space<Scalar>* space,
        MatrixFormVol<Scalar>* custom_projection_jacobian,
        VectorFormVol<Scalar>* custom_projection_residual,
        Scalar* target_vec, double newton_tol, int newton_max_iter)
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
    void OGProjection<Scalar>::project_global(const Space<Scalar>* space,
        MeshFunction<Scalar>* source_meshfn, Scalar* target_vec,
        ProjNormType proj_norm, double newton_tol, int newton_max_iter)
    {
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
          default: throw Hermes::Exceptions::Exception("Unknown space type in OGProjection<Scalar>::project_global().");
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
    void OGProjection<Scalar>::project_global(const Space<Scalar>* space,
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
          default: throw Hermes::Exceptions::Exception("Unknown space type in OGProjection<Scalar>::project_global().");
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
    void OGProjection<Scalar>::project_global(Hermes::vector<const Space<Scalar>*> spaces,
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
    void OGProjection<Scalar>::project_global(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> source_slns,
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
    void OGProjection<Scalar>::project_global(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> source_slns,
        Hermes::vector<Solution<Scalar>*> target_slns, Hermes::vector<ProjNormType> proj_norms, bool delete_old_meshes,
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

    template class HERMES_API OGProjection<double>;
    template class HERMES_API OGProjection<std::complex<double> >;
  }
}