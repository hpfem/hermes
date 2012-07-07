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
    void LocalProjection<Scalar>::project_local(const Space<Scalar>* space, MeshFunction<Scalar>* meshfn,
      Scalar* target_vec, ProjNormType proj_norm)
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

      // Get dimension of the space.
      int ndof = space->get_num_dofs();

      // Erase the target vector.
      memset(target_vec, 0, ndof*sizeof(Scalar));

      // Dump into the first part of target_vec the values of active vertex dofs, then add values
      // of active edge dofs, and finally also values of active bubble dofs.
      // Start with active vertex dofs.
      Mesh* mesh = space->get_mesh();
      Element* e;
      // Go through all active elements in mesh to collect active vertex DOF.
      for_all_active_elements(e, mesh)
      {
        int order = space->get_element_order(e->id);
        if(order > 0)
        {
          for (unsigned int j = 0; j < e->get_nvert(); j++)
          {
            // FIXME - the same vertex is visited several times since it
            // belongs to multiple elements!
            Node* vn = e->vn[j];
            double x = e->vn[j]->x;
            double y = e->vn[j]->y;
            //this->info("Probing vertex %g %g\n", x, y);
            typename Space<Scalar>::NodeData* nd = space->ndata + vn->id;
            if(!vn->is_constrained_vertex() && nd->dof >= 0)
            {
              int dof_num = nd->dof;
              // FIXME: If this is a Solution, the it would be MUCH faster to just
              // retrieve the value from the coefficient vector stored in the Solution.

              // FIXME: Retrieving the value through get_pt_value() is slow and this
              // should be only done if we are dealing with MeshFunction (not a Solution).
              Scalar val = meshfn->get_pt_value(x, y);
              //printf("Found active vertex %g %g, val = %g, dof_num = %d\n", x, y, std::abs(val), dof_num);
              target_vec[dof_num] = val;
            }
          }
        }

        // TODO: Calculate coefficients of edge functions and copy into target_vec_space[i] as well

        // TODO: Calculate coefficients of bubble functions and copy into target_vec_space[i] as well.
      }
    }

    template<typename Scalar>
    void LocalProjection<Scalar>::project_local(const Space<Scalar>* space,
        Solution<Scalar>* source_sln, Solution<Scalar>* target_sln,
        ProjNormType proj_norm)
    {
      int ndof = space->get_num_dofs();
      Scalar* coeff_vec = new Scalar[ndof];
      project_local(space, source_sln, coeff_vec, proj_norm);
      Solution<Scalar>::vector_to_solution(coeff_vec, space, target_sln, proj_norm);
      delete [] coeff_vec;
    }

    template<typename Scalar>
    void LocalProjection<Scalar>::project_local(Hermes::vector<const Space<Scalar>*> spaces,
        Hermes::vector<MeshFunction<Scalar>*> meshfns, Scalar* target_vec,
        Hermes::vector<ProjNormType> proj_norms)
    {
      int n = spaces.size();

      if(n != meshfns.size()) throw Exceptions::LengthException(1, 2, n, meshfns.size());
      if(target_vec == NULL) throw Exceptions::NullException(3);
      if(!proj_norms.empty() && n!=proj_norms.size()) throw Exceptions::LengthException(1, 5, n, proj_norms.size());

      int start_index = 0;
      for (int i = 0; i < n; i++)
      {
        if(proj_norms.empty())
          project_local(spaces[i], meshfns[i], target_vec + start_index);
        else
          project_local(spaces[i], meshfns[i], target_vec + start_index, proj_norms[i]);
        start_index += spaces[i]->get_num_dofs();
      }
    }

    template<typename Scalar>
    void LocalProjection<Scalar>::project_local(Hermes::vector<const Space<Scalar>*> spaces,
        Hermes::vector<Solution<Scalar>*> slns, Scalar* target_vec,
        Hermes::vector<ProjNormType> proj_norms)
    {
      int n = spaces.size();

      // Sanity checks.
      if(n != slns.size()) throw Exceptions::LengthException(1, 2, n, slns.size());
      if(target_vec == NULL) throw Exceptions::NullException(3);
      if(!proj_norms.empty() && n!=proj_norms.size()) throw Exceptions::LengthException(1, 5, n, proj_norms.size());

      int start_index = 0;
      for (int i = 0; i < n; i++)
      {
        if(proj_norms.empty())
          project_local(spaces[i], slns[i], target_vec + start_index);
        else
          project_local(spaces[i], slns[i], target_vec + start_index, proj_norms[i]);
        start_index += spaces[i]->get_num_dofs();
      }
    }

    template<typename Scalar>
    void LocalProjection<Scalar>::project_local(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> source_slns,
      Hermes::vector<Solution<Scalar>*> target_slns, Hermes::vector<ProjNormType> proj_norms, bool delete_old_meshes)
    {
      int n = spaces.size();

      // Sanity checks.
      if(n != source_slns.size()) throw Exceptions::LengthException(1, 2, n, source_slns.size());
      if(n != target_slns.size()) throw Exceptions::LengthException(1, 2, n, target_slns.size());
      if(!proj_norms.empty() && n != proj_norms.size()) throw Exceptions::LengthException(1, 5, n, proj_norms.size());

      for (int i = 0; i < n; i++)
      {
        if(proj_norms.empty())
          project_local(spaces[i], source_slns[i], target_slns[i]);
        else
          project_local(spaces[i], source_slns[i], target_slns[i], proj_norms[i]);
      }
    }

    template class HERMES_API LocalProjection<double>;
    template class HERMES_API LocalProjection<std::complex<double> >;
  }
}