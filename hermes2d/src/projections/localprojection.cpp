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

#include "projections/localprojection.h"
#include "space.h"
#include "discrete_problem.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    int LocalProjection<Scalar>::ndof = 0;

    template<typename Scalar>
    LocalProjection<Scalar>::LocalProjection()
    {
    }

    template<typename Scalar>
    void LocalProjection<Scalar>::project_local(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<MeshFunction<Scalar>*> source_meshfns,
      Scalar* target_vec, Hermes::MatrixSolverType matrix_solver, Hermes::vector<ProjNormType> proj_norms)
    {
      _F_;
      int n = spaces.size();

      if (n!=source_meshfns.size()) throw Exceptions::LengthException(1, 2, n, source_meshfns.size());
      if (target_vec==NULL) throw Exceptions::NullException(3);
      if (!proj_norms.empty() && n!=proj_norms.size()) throw Exceptions::LengthException(1, 5, n, proj_norms.size());

      // Get ndof of each space.
      int* ndof_space = new int[n];
      for (int i=0; i < n; i++) ndof_space[i] = spaces[i]->get_num_dofs();

      // Create a zero target vector for each space. 
      Scalar** target_vec_space = new Scalar*[n];
      for (int i=0; i < n; i++) 
      {
        int ndof_i = ndof_space[i];
        target_vec_space[i] = new Scalar[ndof_i];
        memset(target_vec_space[i], 0, ndof_i*sizeof(Scalar));
      }

      // For each function source_meshfns[i], dump into the first part of target_vec_space[i] the values 
      // of active vertex dofs, then add values of active edge dofs, and finally also values 
      // of active bubble dofs.
      // Start with active vertex dofs.
      for (int i = 0; i < n; i++) 
      { 
        Mesh* mesh_i = spaces[i]->get_mesh();
	Element* e;
        // Go through all active elements in mesh[i]
        // to collect active vertex DOF.
	for_all_active_elements(e, mesh_i)
        {
          int order = spaces[i]->get_element_order(e->id);
          if (order > 0)
          {
            for (unsigned int j = 0; j < e->get_nvert(); j++)
            {
              // FIXME - the same vertex is visited several times since it 
              // belongs to multiple elements!
              Node* vn = e->vn[j];
              double x = e->vn[j]->x;
              double y = e->vn[j]->y;
              //info("Probing vertex %g %g\n", x, y);
              typename Space<Scalar>::NodeData* nd = spaces[i]->ndata + vn->id;
              if (!vn->is_constrained_vertex() && nd->dof >= 0)
              {
                int dof_num = nd->dof;
                // FIXME: If this is a Solution, the it would be MUCH faster to just 
                // retrieve the value from the coefficient vector stored in the Solution.

                // FIXME: Retrieving the value through get_pt_value() is slow and this 
                // should be only done if we are dealing with MeshFunction (not a Solution).
                Scalar val = source_meshfns[i]->get_pt_value(x, y);
                //printf("Found active vertex %g %g, val = %g, dof_num = %d\n", x, y, std::abs(val), dof_num);
                target_vec_space[i][dof_num] = val;
              }
            }
          }
        }

        // TODO: Calculate coefficients of edge functions and copy into target_vec_space[i] as well

        // TODO: Calculate coefficients of bubble functions and copy into target_vec_space[i] as well.

      }

      /*
      // debug
      Solution<double> debug_sln;
      Solution<double>::vector_to_solution(target_vec_space[3], spaces[3], &debug_sln);
      Hermes::Hermes2D::Views::ScalarView sv;
      sv.show(&debug_sln);
      Hermes::Hermes2D::Views::View::wait();
      */

      // Copy all the vectors target_vec_space[i] into target_vec.
      int counter = 0;
      for (int i = 0; i < n; i++) 
      {
        for (int j = 0; j < ndof_space[i]; j++) 
	{
          target_vec[counter + j] = target_vec_space[i][j];
	}
        counter += ndof_space[i];
      }

      // debug
      //for (int i=0; i<counter; i++) 
      //{
      //  printf("dof %d val %g\n", i, std::abs(target_vec[i]));
      //}

      // Clean up.
      for (int i = 0; i < n; i++) delete [] target_vec_space[i];

      return;
    }

    template<typename Scalar>
    void LocalProjection<Scalar>::project_local(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> source_sols,
      Scalar* target_vec, Hermes::MatrixSolverType matrix_solver, Hermes::vector<ProjNormType> proj_norms)
    {
      int n = spaces.size();

      // Sanity checks.
      if (n != source_sols.size()) throw Exceptions::LengthException(1, 2, n, source_sols.size());

      Hermes::vector<MeshFunction<Scalar>*> mesh_fns = Hermes::vector<MeshFunction<Scalar>*>();
      for(unsigned int i = 0; i < source_sols.size(); i++) mesh_fns.push_back(source_sols[i]);
      project_local(spaces, mesh_fns, target_vec, matrix_solver, proj_norms);
    }

    template<typename Scalar>
    void LocalProjection<Scalar>::project_local(Space<Scalar>* space, MeshFunction<Scalar>* source_meshfn,
      Scalar* target_vec, Hermes::MatrixSolverType matrix_solver,
      ProjNormType proj_norm)
    {
      if (proj_norm == HERMES_UNSET_NORM) 
      { 
        SpaceType space_type = space->get_type();
        switch (space_type)
        {
          case HERMES_H1_SPACE: proj_norm = HERMES_H1_NORM; break;
          case HERMES_HCURL_SPACE: proj_norm = HERMES_HCURL_NORM; break;
          case HERMES_HDIV_SPACE: proj_norm = HERMES_HDIV_NORM; break;
          case HERMES_L2_SPACE: proj_norm = HERMES_L2_NORM; break;
          default: error("Unknown space type in OGProjection<Scalar>::project_global().");
        }
      }

      Hermes::vector<Space<Scalar>*> spaces;
      spaces.push_back(space);
      Hermes::vector<MeshFunction<Scalar>*> source_meshfns;
      source_meshfns.push_back(source_meshfn);
      Hermes::vector<ProjNormType> proj_norms;
      proj_norms.push_back(proj_norm);
      project_local(spaces, source_meshfns, target_vec, matrix_solver, proj_norms);
    }

    template<typename Scalar>
    void LocalProjection<Scalar>::project_local(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> sols_src,
      Hermes::vector<Solution<Scalar>*> sols_dest, Hermes::MatrixSolverType matrix_solver,
      Hermes::vector<ProjNormType> proj_norms, bool delete_old_meshes)
    {
      _F_;
      int n = spaces.size();

      // Sanity checks.
      if (n != sols_src.size()) throw Exceptions::LengthException(1, 2, n, sols_src.size());
      if (!proj_norms.empty() && n != proj_norms.size()) throw Exceptions::LengthException(1, 5, n, proj_norms.size());

      // This vector will contain projection norms. 
      Hermes::vector<ProjNormType> proj_norms_new = Hermes::vector<ProjNormType>();

      // If projection norms are not provided, set them
      // to match the type of each space.
      for (int i = 0; i < n; i++)
      {
        ProjNormType norm_i = HERMES_UNSET_NORM;
        if (proj_norms == Hermes::vector<ProjNormType>())
        {
          SpaceType space_type_i = spaces[i]->get_type();
          switch (space_type_i)
          {
            case HERMES_H1_SPACE: norm_i = HERMES_H1_NORM; break;
            case HERMES_HCURL_SPACE: norm_i = HERMES_HCURL_NORM; break;
            case HERMES_HDIV_SPACE: norm_i = HERMES_HDIV_NORM; break;
            case HERMES_L2_SPACE: norm_i = HERMES_L2_NORM; break;
            default: error("Unknown space type in OGProjection<Scalar>::project_global().");
          }
        }
        else norm_i = proj_norms[i];

        proj_norms_new.push_back(norm_i);
      }

      Scalar* target_vec = new Scalar[Space<Scalar>::get_num_dofs(spaces)];
      Hermes::vector<MeshFunction<Scalar>*> ref_slns_mf;
      for (unsigned int i = 0; i < sols_src.size(); i++)
        ref_slns_mf.push_back(static_cast<MeshFunction<Scalar>*>(sols_src[i]));

      LocalProjection<Scalar>::project_local(spaces, ref_slns_mf, target_vec, matrix_solver, proj_norms_new);

      if(delete_old_meshes) 
      {
        for(unsigned int i = 0; i < sols_src.size(); i++)
        {
          delete sols_src[i]->get_mesh();
          sols_src[i]->own_mesh = false;
        }
      }
      
      Solution<Scalar>::vector_to_solutions(target_vec, spaces, sols_dest);

      delete [] target_vec;
    }

    template<typename Scalar>
    void LocalProjection<Scalar>::project_local(Space<Scalar>* space,
        Solution<Scalar>* sol_src, Solution<Scalar>* sol_dest,
        Hermes::MatrixSolverType matrix_solver,
        ProjNormType proj_norm)
    {
      if (proj_norm == HERMES_UNSET_NORM) 
      { 
        SpaceType space_type = space->get_type();
        switch (space_type)
        {
          case HERMES_H1_SPACE: proj_norm = HERMES_H1_NORM; break;
          case HERMES_HCURL_SPACE: proj_norm = HERMES_HCURL_NORM; break;
          case HERMES_HDIV_SPACE: proj_norm = HERMES_HDIV_NORM; break;
          case HERMES_L2_SPACE: proj_norm = HERMES_L2_NORM; break;
          default: error("Unknown space type in OGProjection<Scalar>::project_global().");
        }
      }

      Hermes::vector<Space<Scalar>*> spaces;
      spaces.push_back(space);
      Hermes::vector<Solution<Scalar>*> sols_src;
      sols_src.push_back(sol_src);
      Hermes::vector<Solution<Scalar>*> sols_dest;
      sols_dest.push_back(sol_dest);
      Hermes::vector<ProjNormType> proj_norms;
      proj_norms.push_back(proj_norm);

      project_local(spaces, sols_src, sols_dest, matrix_solver, proj_norms);
    }

    template class HERMES_API LocalProjection<double>;
    template class HERMES_API LocalProjection<std::complex<double> >;
  }
}
