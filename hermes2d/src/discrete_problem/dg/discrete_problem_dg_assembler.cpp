//#define DEBUG_DG_ASSEMBLING
//#define DEBUG_DG_ASSEMBLING_ELEMENT 44
//#define DEBUG_DG_ASSEMBLING_ISURF 3
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

#include "discrete_problem/dg/discrete_problem_dg_assembler.h"
#include "discrete_problem/discrete_problem_thread_assembler.h"

using namespace Hermes::Algebra::DenseMatrixOperations;

namespace Hermes
{
  namespace Hermes2D
  {
    static const std::string H2D_DG_INNER_EDGE = "-1234567";

    template<typename Scalar>
    unsigned int DiscreteProblemDGAssembler<Scalar>::dg_order = 20;

    template<typename Scalar>
    DiscreteProblemDGAssembler<Scalar>::DiscreteProblemDGAssembler(DiscreteProblemThreadAssembler<Scalar>* threadAssembler, const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Hermes::vector<MeshSharedPtr>& meshes)
      : pss(threadAssembler->pss),
      refmaps(threadAssembler->refmaps),
      u_ext(threadAssembler->u_ext),
      als(threadAssembler->als),
      fns(threadAssembler->fns),
      wf(threadAssembler->wf),
      spaces_size(threadAssembler->spaces_size),
      nonlinear(threadAssembler->nonlinear),
      current_mat(threadAssembler->current_mat),
      current_rhs(threadAssembler->current_rhs),
      current_state(NULL),
      selectiveAssembler(threadAssembler->selectiveAssembler),
      do_not_use_cache(threadAssembler->do_not_use_cache),
      spaces(spaces),
      meshes(meshes)
    {
      this->DG_matrix_forms_present = false;
      this->DG_vector_forms_present = false;
      if(this->wf)
      {
        if(!this->wf->mfDG.empty())
          this->DG_matrix_forms_present = true;

        if(!this->wf->vfDG.empty())
          this->DG_vector_forms_present = true;
      }

      if(DG_matrix_forms_present)
      {
        npss = new PrecalcShapeset*[spaces_size];
        nrefmaps = new RefMap*[spaces_size];

        for (unsigned int j = 0; j < spaces_size; j++)
        {
          npss[j] = new PrecalcShapeset(spaces[j]->shapeset);
          nrefmaps[j] = new RefMap();
        }
      }
    }

    template<typename Scalar>
    NeighborSearch<Scalar>* DiscreteProblemDGAssembler<Scalar>::get_neighbor_search_ext(WeakForm<Scalar>* wf, NeighborSearch<Scalar>** neighbor_searches, int index)
    {
      return neighbor_searches[index + wf->get_neq()];
    }

    template<typename Scalar>
    DiscreteProblemDGAssembler<Scalar>::~DiscreteProblemDGAssembler()
    {
      if(DG_matrix_forms_present)
      {
        for (unsigned int j = 0; j < spaces_size; j++)
        {
          delete npss[j];
          delete nrefmaps[j];
        }
        delete [] npss;
        delete [] nrefmaps;
      }
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::init_assembling_one_state(Traverse::State* current_state_)
    {
      this->current_state = current_state_;

      this->neighbor_searches = new NeighborSearch<Scalar>**[this->current_state->rep->nvert];
      for(int i = 0; i < this->current_state->rep->nvert; i++)
        this->neighbor_searches[i] = new NeighborSearch<Scalar>*[this->current_state->num];
      this->num_neighbors = new int[this->current_state->rep->nvert];
      processed = new bool*[current_state->rep->nvert];

      if(DG_matrix_forms_present)
      {
        for (unsigned int i = 0; i < spaces_size; i++)
        {
          npss[i]->set_quad_2d(&g_quad_2d_std);
          nrefmaps[i]->set_quad_2d(&g_quad_2d_std);
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::assemble_one_state()
    {
#pragma omp critical (DG)
      {
        for(unsigned int i = 0; i < current_state->num; i++)
          current_state->e[i]->visited = true;

        for(current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
        {
          if(!current_state->bnd[current_state->isurf])
          {
            // If this edge is an inter-element one on all meshes.
            if(!init_neighbors(neighbor_searches[current_state->isurf], current_state))
              continue;

            // Create a multimesh tree;
            MultimeshDGNeighborTree<Scalar>::process_edge(neighbor_searches[current_state->isurf], this->current_state->num, this->num_neighbors[current_state->isurf], this->processed[current_state->isurf]);
          }
        }
        for(current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
        {
          if(!current_state->bnd[current_state->isurf])
          {
#ifdef DEBUG_DG_ASSEMBLING
            debug();
#endif
            for(unsigned int neighbor_i = 0; neighbor_i < num_neighbors[current_state->isurf]; neighbor_i++)
            {
              if(!DG_vector_forms_present && processed[current_state->isurf][neighbor_i])
                continue;

              // DG-inner-edge-wise parameters for WeakForm.
              wf->set_active_DG_state(current_state->e, current_state->isurf);

              assemble_one_neighbor(processed[current_state->isurf][neighbor_i], neighbor_i, neighbor_searches[current_state->isurf]);
            }

            deinit_neighbors(neighbor_searches[current_state->isurf], current_state);
          }
          else
            processed[current_state->isurf] = NULL;
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::deinit_assembling_one_state()
    {
      for(int i = 0; i < this->current_state->rep->nvert; i++)
      {
        delete [] neighbor_searches[i];
        if(processed[i])
          delete [] processed[i];
      }
      delete [] neighbor_searches;
      delete [] num_neighbors;
      delete [] processed;
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::init_assembling_one_neighbor()
    {
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::assemble_one_neighbor(bool edge_processed, unsigned int neighbor_i, NeighborSearch<Scalar>** current_neighbor_searches)
    {
      // Set the active segment in all NeighborSearches
      for(unsigned int i = 0; i < this->current_state->num; i++)
      {
        NeighborSearch<Scalar>* ns = current_neighbor_searches[i];
        ns->active_segment = neighbor_i;
        ns->neighb_el = ns->neighbors[neighbor_i];
        ns->neighbor_edge = ns->neighbor_edges[neighbor_i];
      }

      // Push all the necessary transformations to all functions of this stage.
      // The important thing is that the transformations to the current subelement are already there.
      for(unsigned int fns_i = 0; fns_i < current_state->num; fns_i++)
      {
        NeighborSearch<Scalar>* ns = current_neighbor_searches[fns_i];
        if(neighbor_i < ns->central_transformations_alloc_size && ns->central_transformations[neighbor_i])
          ns->central_transformations[neighbor_i]->apply_on(fns[fns_i]);
      }

      // For neighbor psss.
      if(current_mat && DG_matrix_forms_present && !edge_processed)
      {
        for(unsigned int idx_i = 0; idx_i < spaces_size; idx_i++)
        {
          NeighborSearch<Scalar>* ns = current_neighbor_searches[idx_i];
          npss[idx_i]->set_active_element((*ns->get_neighbors())[neighbor_i]);
          if(neighbor_i < ns->neighbor_transformations_alloc_size && ns->neighbor_transformations[neighbor_i])
            ns->neighbor_transformations[neighbor_i]->apply_on(npss[idx_i]);
        }
      }

      // Also push the transformations to the refmaps.
      for (unsigned int i = 0; i < spaces_size; i++)
      {
        refmaps[i]->force_transform(pss[i]->get_transform(), pss[i]->get_ctm());

        // Neighbor.
        if(current_mat && DG_matrix_forms_present && !edge_processed)
        {
          nrefmaps[i]->set_active_element(npss[i]->get_active_element());
          nrefmaps[i]->force_transform(npss[i]->get_transform(), npss[i]->get_ctm());
        }
      }

      /***/
      // The computation takes place here.
      typename NeighborSearch<Scalar>::ExtendedShapeset** ext_asmlist = new typename NeighborSearch<Scalar>::ExtendedShapeset*[this->spaces_size];
      int n_quadrature_points;
      Geom<double>** geometry = new Geom<double>*[this->spaces_size];
      double** jacobian_x_weights = new double*[this->spaces_size];
      Geom<double>** e = new Geom<double>*[this->spaces_size];
      DiscontinuousFunc<double>*** testFunctions = new DiscontinuousFunc<double>**[this->spaces_size];

      // Create the extended shapeset on the union of the central element and its current neighbor.
      int order = DiscreteProblemDGAssembler<Scalar>::dg_order;
      int order_base = DiscreteProblemDGAssembler<Scalar>::dg_order;
      for (unsigned int i = 0; i < this->spaces_size; i++)
      {
        current_neighbor_searches[i]->set_quad_order(order);
        order_base = order;
        n_quadrature_points = init_surface_geometry_points(refmaps, this->spaces_size, order_base, current_state->isurf, current_state->rep->marker, geometry[i], jacobian_x_weights[i]);
        e[i] = new InterfaceGeom<double>(geometry[i], current_neighbor_searches[i]->neighb_el->marker, current_neighbor_searches[i]->neighb_el->id, current_neighbor_searches[i]->neighb_el->get_diameter());

        if(current_mat && DG_matrix_forms_present && !edge_processed)
        {
          ext_asmlist[i] = current_neighbor_searches[i]->create_extended_asmlist(spaces[i], als[i]);
          testFunctions[i] = new DiscontinuousFunc<double>*[ext_asmlist[i]->cnt];
          for (int func_i = 0; func_i < ext_asmlist[i]->cnt; func_i++)
          {
            if(ext_asmlist[i]->dof[func_i] < 0)
              continue;

            // Choose the correct shapeset for the test function.
            if(ext_asmlist[i]->has_support_on_neighbor(func_i))
            {
              npss[i]->set_active_shape(ext_asmlist[i]->neighbor_al->idx[func_i - ext_asmlist[i]->central_al->cnt]);
              testFunctions[i][func_i] = new DiscontinuousFunc<double>(init_fn(npss[i], nrefmaps[i], current_neighbor_searches[i]->get_quad_eo(true)), true, current_neighbor_searches[i]->neighbor_edge.orientation);
            }
            else
            {
              pss[i]->set_active_shape(ext_asmlist[i]->central_al->idx[func_i]);
              testFunctions[i][func_i] = new DiscontinuousFunc<double>(init_fn(pss[i], refmaps[i], current_neighbor_searches[i]->get_quad_eo(false)), false, current_neighbor_searches[i]->neighbor_edge.orientation);
            }
          }
        }
      }

      DiscontinuousFunc<Scalar>** ext = init_ext_fns(wf->ext, current_neighbor_searches, order);

      DiscontinuousFunc<Scalar>** u_ext_func = new DiscontinuousFunc<Scalar>*[this->spaces_size];
      if(this->nonlinear)
      {
        if(u_ext)
        {
          for(int u_ext_func_i = 0; u_ext_func_i < this->spaces_size; u_ext_func_i++)
            if(u_ext[u_ext_func_i])
            {
              current_neighbor_searches[u_ext_func_i]->set_quad_order(order);
              u_ext_func[u_ext_func_i]  = current_neighbor_searches[u_ext_func_i]->init_ext_fn(u_ext[u_ext_func_i]);
            }
            else
              u_ext_func[u_ext_func_i] = NULL;
        }
        else
          for(int u_ext_func_i = 0; u_ext_func_i < this->spaces_size; u_ext_func_i++)
            u_ext_func[u_ext_func_i] = NULL;
      }

      if(current_mat && DG_matrix_forms_present && !edge_processed)
      {
        for(int current_mfsurf_i = 0; current_mfsurf_i < wf->mfDG.size(); current_mfsurf_i++)
        {
          if(!this->selectiveAssembler->form_to_be_assembled((MatrixForm<Scalar>*)wf->mfDG[current_mfsurf_i], current_state))
            continue;

          MatrixFormDG<Scalar>* mfs = wf->mfDG[current_mfsurf_i];

          int m = mfs->i;
          int n = mfs->j;

          // Precalc shapeset and refmaps used for the evaluation.
          bool support_neigh_u, support_neigh_v;
          typename NeighborSearch<Scalar>::ExtendedShapeset* ext_asmlist_u = ext_asmlist[n];
          typename NeighborSearch<Scalar>::ExtendedShapeset* ext_asmlist_v = ext_asmlist[m];

          Scalar **local_stiffness_matrix = new_matrix<Scalar>(std::max(ext_asmlist_u->cnt, ext_asmlist_v->cnt));
          for (int i = 0; i < ext_asmlist_v->cnt; i++)
          {
            if(ext_asmlist_v->dof[i] < 0)
              continue;

            support_neigh_v = ext_asmlist_v->has_support_on_neighbor(i);

            for (int j = 0; j < ext_asmlist_u->cnt; j++)
            {
              if(ext_asmlist_u->dof[j] >= 0)
              {
                // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
                DiscontinuousFunc<double>* u = testFunctions[n][j];
                DiscontinuousFunc<double>* v = testFunctions[m][i];

                Scalar res = mfs->value(n_quadrature_points, jacobian_x_weights[n], u_ext_func, u, v, e[n], ext) * mfs->scaling_factor;

                support_neigh_u = ext_asmlist_u->has_support_on_neighbor(j);

                Scalar val = 0.5 * res * (support_neigh_u ? ext_asmlist_u->neighbor_al->coef[j - ext_asmlist_u->central_al->cnt]: ext_asmlist_u->central_al->coef[j])
                  * (support_neigh_v ? ext_asmlist_v->neighbor_al->coef[i - ext_asmlist_v->central_al->cnt]: ext_asmlist_v->central_al->coef[i]);

                local_stiffness_matrix[i][j] = val;
              }
            }
          }

          current_mat->add(ext_asmlist_v->cnt, ext_asmlist_u->cnt, local_stiffness_matrix, ext_asmlist_v->dof, ext_asmlist_u->dof);

          delete [] local_stiffness_matrix;
        }
      }

      if(current_mat && DG_matrix_forms_present && !edge_processed)
      {
        for(int i = 0; i < this->spaces_size; i++)
        {
          for (int func_i = 0; func_i < ext_asmlist[i]->cnt; func_i++)
          {
            if(ext_asmlist[i]->dof[func_i] < 0)
              continue;
            testFunctions[i][func_i]->free_fn();
            delete testFunctions[i][func_i];
          }        
          delete ext_asmlist[i];
          delete [] testFunctions[i];
        }
      }
      delete [] testFunctions;
      delete [] ext_asmlist;

      if(current_rhs && DG_vector_forms_present)
      {
        for (unsigned int ww = 0; ww < wf->vfDG.size(); ww++)
        {
          VectorFormDG<Scalar>* vfs = wf->vfDG[ww];
          if(vfs->areas[0] != H2D_DG_INNER_EDGE)
            continue;

          int n = vfs->i;

          if(!this->selectiveAssembler->form_to_be_assembled((VectorForm<Scalar>*)vfs, current_state))
            continue;

          NeighborSearch<Scalar>* current_neighbor_searches_v = current_neighbor_searches[n];

          // Here we use the standard pss, possibly just transformed by NeighborSearch.
          for (unsigned int dof_i = 0; dof_i < als[n]->cnt; dof_i++)
          {
            if(als[n]->dof[dof_i] < 0)
              continue;
            pss[n]->set_active_shape(als[n]->idx[dof_i]);

            Func<double>* v = init_fn(pss[n], refmaps[n], current_neighbor_searches_v->get_quad_eo());

            current_rhs->add(als[n]->dof[dof_i], 0.5 * vfs->value(n_quadrature_points, jacobian_x_weights[n], u_ext_func, v, e[n], ext) * vfs->scaling_factor * als[n]->coef[dof_i]);

            v->free_fn();
            delete v;
          }
        }
      }

      if(ext)
      {
        for(unsigned int i = 0; i < wf->ext.size(); i++)
        {
          ext[i]->free_fn();
          delete ext[i];
        }
        delete [] ext;
      }

      if(this->nonlinear)
      {
        if(u_ext)
        {
          for(int u_ext_i = 0; u_ext_i < this->spaces_size; u_ext_i++)
            if(u_ext[u_ext_i])
            {
              u_ext_func[u_ext_i]->free_fn();
              delete u_ext_func[u_ext_i];
            }
        }
      }

      delete [] u_ext_func;


      for(int i = 0; i < this->spaces_size; i++)
      {
        if(this->spaces[i]->get_type() != HERMES_L2_SPACE)
          continue;
        delete [] jacobian_x_weights[i];
        e[i]->free();
        delete e[i];
      }

      delete [] geometry;
      delete [] jacobian_x_weights;
      delete [] e;

      // This is just cleaning after ourselves.
      // Clear the transformations from the RefMaps and all functions.
      for(unsigned int fns_i = 0; fns_i < current_state->num; fns_i++)
      {
        fns[fns_i]->set_transform(current_neighbor_searches[fns_i]->original_central_el_transform);
      }

      // Also clear the transformations from the slave psss and refmaps.
      for (unsigned int i = 0; i < spaces_size; i++)
        refmaps[i]->force_transform(pss[i]->get_transform(), pss[i]->get_ctm());
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::deinit_assembling_one_neighbor()
    {
    }

    template<typename Scalar>
    DiscontinuousFunc<Scalar>** DiscreteProblemDGAssembler<Scalar>::init_ext_fns(Hermes::vector<MeshFunctionSharedPtr<Scalar> > ext,
      NeighborSearch<Scalar>** current_neighbor_searches, int order)
    {
      DiscontinuousFunc<Scalar>** ext_fns = new DiscontinuousFunc<Scalar>*[ext.size()];
      for(unsigned int j = 0; j < ext.size(); j++)
      {
        NeighborSearch<Scalar>* ns = get_neighbor_search_ext(this->wf, current_neighbor_searches, j);
        ns->set_quad_order(order);
        ext_fns[j] = ns->init_ext_fn(ext[j].get());
      }

      return ext_fns;
    }

    template<typename Scalar>
    bool DiscreteProblemDGAssembler<Scalar>::init_neighbors(NeighborSearch<Scalar>** current_neighbor_searches, Traverse::State* current_state)
    {
      // Initialize the NeighborSearches.
      bool DG_intra = false;
      for(unsigned int i = 0; i < current_state->num; i++)
      {
        bool existing_ns = false;
        for(int j = i - 1; j >= 0; j--)
          if(current_state->e[i] == current_state->e[j])
          {
            current_neighbor_searches[i] = current_neighbor_searches[j];
            existing_ns = true;
            break;
          }
          if(!existing_ns)
          {
            NeighborSearch<Scalar>* ns = new NeighborSearch<Scalar>(current_state->e[i], this->meshes[i]);
            ns->original_central_el_transform = current_state->sub_idx[i];
            current_neighbor_searches[i] = ns;
            if(current_neighbor_searches[i]->set_active_edge_multimesh(current_state->isurf) && (i >= this->spaces_size || spaces[i]->get_type() == HERMES_L2_SPACE))
              DG_intra = true;
            current_neighbor_searches[i]->clear_initial_sub_idx();
          }
      }

      return DG_intra;
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::deinit_neighbors(NeighborSearch<Scalar>** current_neighbor_searches, Traverse::State* current_state)
    {
      for(unsigned int i = 0; i < current_state->num; i++)
      {
        bool existing_ns = false;
        for(int j = i - 1; j >= 0; j--)
          if(current_state->e[i] == current_state->e[j])
          {
            existing_ns = true;
            break;
          }
          if(!existing_ns)
            delete current_neighbor_searches[i];
      }
    }

#ifdef DEBUG_DG_ASSEMBLING
    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::debug()
    {
#pragma omp critical (debug_DG)
      {
        int id = 0;
        bool pass = true;
        if(DEBUG_DG_ASSEMBLING_ELEMENT != -1)
        {
          for(unsigned int i = 0; i < this->current_state->num; i++)
            if(neighbor_searches[current_state->isurf][i]->central_el->id == DEBUG_DG_ASSEMBLING_ELEMENT)
              pass = false;
        }
        else
          pass = false;

        if(!pass)
          if(DEBUG_DG_ASSEMBLING_ISURF != -1)
            if(current_state->isurf != DEBUG_DG_ASSEMBLING_ISURF)
              pass = true;

        if(!pass)
        {
          for(unsigned int i = 0; i < this->current_state->num; i++)
          {
            NeighborSearch<Scalar>* ns = neighbor_searches[current_state->isurf][i];
            std::cout << "The " << ++id << "-th Neighbor search: " << ns->n_neighbors << " neighbors." << std::endl;
            std::cout << "\tCentral element: " << ns->central_el->id << ", Isurf: " << current_state->isurf << ", Original sub_idx: " << ns->original_central_el_transform << std::endl;
            for(int j = 0; j < ns->n_neighbors; j++)
            {
              std::cout << '\t' << "The " << j << "-th neighbor element: " << ns->neighbors[j]->id << std::endl;
              if(ns->central_transformations[j])
              {
                std::cout << '\t' << "Central transformations: " << std::endl;
                for(int k = 0; k < ns->central_transformations[j]->num_levels; k++)
                  std::cout << '\t' << '\t' << ns->central_transformations[j]->transf[k] << std::endl;
              }
              if(ns->neighbor_transformations[j])
              {
                std::cout << '\t' << "Neighbor transformations: " << std::endl;
                for(int k = 0; k < ns->neighbor_transformations[j]->num_levels; k++)
                  std::cout << '\t' << '\t' << ns->neighbor_transformations[j]->transf[k] << std::endl;
              }
            }
          }
        }
      }
    }
#endif

    template class HERMES_API DiscreteProblemDGAssembler<double>;
    template class HERMES_API DiscreteProblemDGAssembler<std::complex<double> >;
  }
}
