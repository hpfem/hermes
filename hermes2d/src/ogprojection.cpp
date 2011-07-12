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

#include "ogprojection.h"
#include "space.h"
#include "discrete_problem.h"
#include "newton_solver.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    void OGProjection<Scalar>::project_internal(Hermes::vector<Space<Scalar>*> spaces, WeakForm<Scalar>* wf,
      Scalar* target_vec, Hermes::MatrixSolverType matrix_solver_type)
    {
      _F_;
      unsigned int n = spaces.size();

      // sanity checks
      for (unsigned int i = 0; i < n; i++) 
        if(spaces[i] == NULL) 
          error("this->spaces[%d] == NULL in project_internal().", i);

      // this is needed since spaces may have their DOFs enumerated only locally.
      int ndof = Space<Scalar>::assign_dofs(spaces);

      // Initialize DiscreteProblem.
      DiscreteProblem<Scalar> dp(wf, spaces);

      // Initial coefficient vector for the Newton's method.  
      Scalar* coeff_vec = new Scalar[ndof];
      memset(coeff_vec, 0, ndof*sizeof(Scalar));

      // Perform Newton's iteration.
      NewtonSolver<Scalar> newton(&dp, matrix_solver_type);
      // No output for the Newton's loop.
      newton.set_verbose_output(false);
      if (!newton.solve(coeff_vec)) 
        error("Newton's iteration failed.");

      delete [] coeff_vec;

      if (target_vec != NULL)
        for (int i = 0; i < ndof; i++) 
          target_vec[i] = newton.get_sln_vector()[i];
    }

    template<typename Scalar>
    void OGProjection<Scalar>::project_global(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<MeshFunction<Scalar>*> source_meshfns,
      Scalar* target_vec, Hermes::MatrixSolverType matrix_solver_type, Hermes::vector<ProjNormType> proj_norms)
    {
      _F_;
      int n = spaces.size();

      // define temporary projection weak form
      WeakForm<Scalar>* proj_wf = new WeakForm<Scalar>(n);
      
      // For error detection.
      bool* found = new bool[n];
      for (int i = 0; i < n; i++) 
        found[i] = false;

      for (int i = 0; i < n; i++)
      {
        ProjNormType norm = HERMES_UNSET_NORM;
        if (proj_norms == Hermes::vector<ProjNormType>()) 
        {
          SpaceType space_type = spaces[i]->get_type();
          switch (space_type) 
          {
          case HERMES_H1_SPACE: norm = HERMES_H1_NORM; break;
          case HERMES_HCURL_SPACE: norm = HERMES_HCURL_NORM; break;
          case HERMES_HDIV_SPACE: norm = HERMES_HDIV_NORM; break;
          case HERMES_L2_SPACE: norm = HERMES_L2_NORM; break;
          default: error("Unknown space type in OGProjection<Scalar>::project_global().");
          }
        }
        else norm = proj_norms[i];

        // FIXME - memory leak - create Projection class and encapsulate this function project_global(...)
        // maybe in more general form
        found[i] = true;
        // Jacobian.
        proj_wf->add_matrix_form(new ProjectionMatrixFormVol(i, i, norm));
        // Residual.
        proj_wf->add_vector_form(new ProjectionVectorFormVol(i, source_meshfns[i], norm));
      }
      for (int i=0; i < n; i++)
      {
        if (!found[i])
        {
          warn("Index of component: %d\n", i);
          error("Wrong projection norm in project_global().");
        }
      }

      project_internal(spaces, proj_wf, target_vec, matrix_solver_type);
    }

    template<typename Scalar>
    void OGProjection<Scalar>::project_global(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> source_sols,
      Scalar* target_vec, Hermes::MatrixSolverType matrix_solver_type, Hermes::vector<ProjNormType> proj_norms)
    {
      Hermes::vector<MeshFunction<Scalar>*> mesh_fns;
      for(unsigned int i = 0; i < source_sols.size(); i++)
        mesh_fns.push_back(source_sols[i]);
      project_global(spaces, mesh_fns, target_vec, matrix_solver_type, proj_norms);
    }

    template<typename Scalar>
    void OGProjection<Scalar>::project_global(Space<Scalar>* space, MeshFunction<Scalar>* source_meshfn,
      Scalar* target_vec, Hermes::MatrixSolverType matrix_solver_type,
      ProjNormType proj_norm)
    {
      Hermes::vector<Space<Scalar>*> spaces;
      spaces.push_back(space);
      Hermes::vector<MeshFunction<Scalar>*> source_meshfns;
      source_meshfns.push_back(source_meshfn);
      Hermes::vector<ProjNormType> proj_norms;
      proj_norms.push_back(proj_norm);
      project_global(spaces, source_meshfns, target_vec, matrix_solver_type, proj_norms);
    }

    template<typename Scalar>
    void OGProjection<Scalar>::project_global(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> sols_src,
      Hermes::vector<Solution<Scalar>*> sols_dest, Hermes::MatrixSolverType matrix_solver_type,
      Hermes::vector<ProjNormType> proj_norms, bool delete_old_meshes)
    {
      _F_;

      Scalar* target_vec = new Scalar[Space<Scalar>::get_num_dofs(spaces)];
      Hermes::vector<MeshFunction<Scalar>*> ref_slns_mf;
      for (unsigned int i = 0; i < sols_src.size(); i++)
        ref_slns_mf.push_back(static_cast<MeshFunction<Scalar>*>(sols_src[i]));

      OGProjection<Scalar>::project_global(spaces, ref_slns_mf, target_vec, matrix_solver_type, proj_norms);

      if(delete_old_meshes)
        for(unsigned int i = 0; i < sols_src.size(); i++) 
        {
          delete sols_src[i]->get_mesh();
          sols_src[i]->own_mesh = false;
        }

        Solution<Scalar>::vector_to_solutions(target_vec, spaces, sols_dest);

        delete [] target_vec;
    }

    template<typename Scalar>
    void OGProjection<Scalar>::project_global(Space<Scalar>* space,
      Solution<Scalar>* sol_src, Solution<Scalar>* sol_dest,
      Hermes::MatrixSolverType matrix_solver_type,
      ProjNormType proj_norm)
    {
      Hermes::vector<Space<Scalar>*> spaces;
      spaces.push_back(space);
      Hermes::vector<Solution<Scalar>*> sols_src;
      sols_src.push_back(sol_src);
      Hermes::vector<Solution<Scalar>*> sols_dest;
      sols_dest.push_back(sol_dest);
      Hermes::vector<ProjNormType> proj_norms;
      if(proj_norm != HERMES_UNSET_NORM)
        proj_norms.push_back(proj_norm);

      project_global(spaces, sols_src, sols_dest, matrix_solver_type, proj_norms);
    }

    template<typename Scalar>
    void OGProjection<Scalar>::project_global(Hermes::vector<Space<Scalar>*> spaces,
      Hermes::vector<MatrixFormVol<Scalar> *> mfvol,
      Hermes::vector<VectorFormVol<Scalar> *> vfvol,
      Hermes::vector<MeshFunction<Scalar>*> source_meshfns,
      Scalar* target_vec, Hermes::MatrixSolverType matrix_solver_type)
    {
      _F_;
      unsigned int n = spaces.size();
      unsigned int n_biforms = mfvol.size();
      if (n_biforms == 0)
        error("Please use the simpler version of project_global with the argument Hermes::vector<ProjNormType> proj_norms if you do not provide your own projection norm.");
      if (n_biforms != vfvol.size())
        error("Mismatched numbers of projection forms in project_global().");
      if (n != n_biforms)
        error("Mismatched numbers of projected functions and projection forms in project_global().");

      // This is needed since spaces may have their DOFs enumerated only locally
      // when they come here.
      int ndof = Space<Scalar>::assign_dofs(spaces);

      // Define projection weak form.
      WeakForm<Scalar>* proj_wf = new WeakForm<Scalar>(n);
      for (unsigned int i = 0; i < n; i++) 
      {
        proj_wf->add_matrix_form(mfvol[i]);
        // FIXME
        // proj_wf->add_vector_form(i, proj_liforms[i].first, proj_liforms[i].second, HERMES_ANY, source_meshfns[i]);
      }

      project_internal(spaces, proj_wf, target_vec, matrix_solver_type);
    }

    template class HERMES_API OGProjection<double>;
    template class HERMES_API OGProjection<std::complex<double> >;
  }
}