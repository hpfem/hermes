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

#include "discrete_problem.h"
#include "function/exact_solution.h"
#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "global.h"
#include "integrals/h1.h"
#include "quadrature/limit_order.h"
#include "mesh/traverse.h"
#include "space/space.h"
#include "shapeset/precalc.h"
#include "mesh/refmap.h"
#include "function/solution.h"
#include "neighbor.h"
#include "api2d.h"

using namespace Hermes::Algebra::DenseMatrixOperations;

namespace Hermes
{
  namespace Hermes2D
  {
    static const std::string H2D_DG_INNER_EDGE = "-1234567";

    template<typename Scalar>
    DiscreteProblemDGAssembler<Scalar>::DiscreteProblemDGAssembler(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> > spaces)
      : Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>(wf)
    {
      this->set_spaces(spaces);
      init();
    }

    template<typename Scalar>
    DiscreteProblemDGAssembler<Scalar>::DiscreteProblemDGAssembler(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar> space)
      : Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>(wf)
    {

      this->set_space(space);
      init();
    }

    template<typename Scalar>
    DiscreteProblemDGAssembler<Scalar>::DiscreteProblemDGAssembler()
    {
      init();
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::init()
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
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::init_assembling()
    {
    }

    template<typename Scalar>
    bool DiscreteProblemDGAssembler<Scalar>::isOkay() const
    {
      if(!this->wf)
        return false;
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> > spacesToSet)
    {
      for(unsigned int i = 0; i < spacesToSet.size(); i++)
      {
        if(!spacesToSet[i])
          throw Exceptions::NullException(0, i);

        spacesToSet[i]->check();
      }

      if(this->spaces.size() != spacesToSet.size() && this->spaces.size() > 0)
        throw Hermes::Exceptions::LengthException(0, spacesToSet.size(), this->spaces.size());

      this->spaces = spacesToSet;
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::set_space(SpaceSharedPtr<Scalar> space)
    {
      Hermes::vector<SpaceSharedPtr<Scalar> > spaces;
      spaces.push_back(space);
      this->set_spaces(spaces);
    }
    
    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::assemble_one_DG_state(PrecalcShapeset** current_pss, RefMap** current_refmaps,  Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als,
      Traverse::State* current_state, Hermes::vector<MatrixFormDG<Scalar>*> current_mfDG, Hermes::vector<VectorFormDG<Scalar>*> current_vfDG, Transformable** fn, WeakForm<Scalar>* current_wf)
    {
      // Determine the minimum mesh seq.
      unsigned int min_dg_mesh_seq = 0;
      for(unsigned int i = 0; i < spaces.size(); i++)
        if(spaces[i]->get_mesh()->get_seq() < min_dg_mesh_seq || i == 0)
          min_dg_mesh_seq = spaces[i]->get_mesh()->get_seq();

      // Create neighbor psss, refmaps.
      std::map<unsigned int, PrecalcShapeset *> npss;
      std::map<unsigned int, PrecalcShapeset *> nspss;
      std::map<unsigned int, RefMap *> nrefmap;

      // Initialize neighbor precalc shapesets and refmaps.
      // This is only needed when there are matrix DG forms present.
      if(DG_matrix_forms_present)
      {
        for (unsigned int i = 0; i < spaces.size(); i++)
        {
          PrecalcShapeset* new_ps = new PrecalcShapeset(spaces[i]->shapeset);
          npss.insert(std::pair<unsigned int, PrecalcShapeset*>(i, new_ps));
          PrecalcShapeset* new_pss = new PrecalcShapeset(npss[i]);
          nspss.insert(std::pair<unsigned int, PrecalcShapeset*>(i, new_pss));
          RefMap* new_rm = new RefMap();
          new_rm->set_quad_2d(&g_quad_2d_std);
          nrefmap.insert(std::pair<unsigned int, RefMap*>(i, new_rm));
        }
      }

      bool** processed = new bool*[current_state->rep->nvert];
      LightArray<NeighborSearch<Scalar>*>** neighbor_searches = new LightArray<NeighborSearch<Scalar>*>*[current_state->rep->nvert];
      (5);
      unsigned int* num_neighbors = new unsigned int[current_state->rep->nvert];

      bool intra_edge_passed_DG[H2D_MAX_NUMBER_VERTICES];
      for(int a = 0; a < H2D_MAX_NUMBER_VERTICES; a++)
        intra_edge_passed_DG[a] = false;

#pragma omp critical (DG)
      {
        for(unsigned int i = 0; i < current_state->num; i++)
          current_state->e[i]->visited = true;

        for(current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
        {
          bool inner_edge_for_dg = false;
          for(int i = 0; i < this->spaces.size(); i++)
            if(current_state->e[i]->en[current_state->isurf]->marker == 0)
              inner_edge_for_dg = true;
          if(inner_edge_for_dg)
          {
            neighbor_searches[current_state->isurf] = new LightArray<NeighborSearch<Scalar>*>(5);

            if(!init_neighbors((*neighbor_searches[current_state->isurf]), current_state, min_dg_mesh_seq))
            {
              intra_edge_passed_DG[current_state->isurf] = true;
              continue;
            }
            // Create a multimesh tree;
            MultimeshDGNeighborTreeNode* root = new MultimeshDGNeighborTreeNode(NULL, 0);
            build_multimesh_tree(root, (*neighbor_searches[current_state->isurf]));

#ifdef DEBUG_DG_ASSEMBLING
#pragma omp critical (debug_DG)
            {
              int id = 0;
              bool pass = true;
              if(DEBUG_DG_ASSEMBLING_ELEMENT != -1)
              {
                for(unsigned int i = 0; i < (*neighbor_searches[current_state->isurf]).get_size(); i++)
                  if((*neighbor_searches[current_state->isurf]).present(i))
                    if((*neighbor_searches[current_state->isurf]).get(i)->central_el->id == DEBUG_DG_ASSEMBLING_ELEMENT)
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
                for(unsigned int i = 0; i < (*neighbor_searches[current_state->isurf]).get_size(); i++)
                {
                  if((*neighbor_searches[current_state->isurf]).present(i))
                  {
                    NeighborSearch<Scalar>* ns = (*neighbor_searches[current_state->isurf]).get(i);
                    std::cout << (std::string)"The " << ++id << (std::string)"-th Neighbor search:: " << (std::string)"Central element: " << ns->central_el->id << (std::string)", Isurf: " << current_state->isurf << (std::string)", Original sub_idx: " << ns->original_central_el_transform << std::endl;
                    for(int j = 0; j < ns->n_neighbors; j++)
                    {
                      std::cout << '\t' << (std::string)"The " << j << (std::string)"-th neighbor element: " << ns->neighbors[j]->id << std::endl;
                      if(ns->central_transformations.present(j))
                      {
                        std::cout << '\t' << (std::string)"Central transformations: " << std::endl;
                        for(int k = 0; k < ns->central_transformations.get(j)->num_levels; k++)
                          std::cout << '\t' << '\t' << ns->central_transformations.get(j)->transf[k] << std::endl;
                      }
                      if(ns->neighbor_transformations.present(j))
                      {
                        std::cout << '\t' << (std::string)"Neighbor transformations: " << std::endl;
                        for(int k = 0; k < ns->neighbor_transformations.get(j)->num_levels; k++)
                          std::cout << '\t' << '\t' << ns->neighbor_transformations.get(j)->transf[k] << std::endl;
                      }
                    }
                  }
                }
              }
            }
#endif

            // Update all NeighborSearches according to the multimesh tree.
            // After this, all NeighborSearches in neighbor_searches should have the same count
            // of neighbors and proper set of transformations
            // for the central and the neighbor element(s) alike.
            // Also check that every NeighborSearch has the same number of neighbor elements.
            num_neighbors[current_state->isurf] = 0;
            for(unsigned int i = 0; i < (*neighbor_searches[current_state->isurf]).get_size(); i++)
            {
              if((*neighbor_searches[current_state->isurf]).present(i))
              {
                NeighborSearch<Scalar>* ns = (*neighbor_searches[current_state->isurf]).get(i);
                update_neighbor_search(ns, root);
                if(num_neighbors[current_state->isurf] == 0)
                  num_neighbors[current_state->isurf] = ns->n_neighbors;
                if(ns->n_neighbors != num_neighbors[current_state->isurf])
                  throw Hermes::Exceptions::Exception("Num_neighbors of different NeighborSearches not matching in assemble_one_DG_state().");
              }
            }

            // Delete the multimesh tree;
            delete root;

            processed[current_state->isurf] = new bool[num_neighbors[current_state->isurf]];

            for(unsigned int neighbor_i = 0; neighbor_i < num_neighbors[current_state->isurf]; neighbor_i++)
            {
              // If the active segment has already been processed (when the neighbor element was assembled), it is skipped.
              // We test all neighbor searches, because in the case of intra-element edge, the neighboring (the same as central) element
              // will be marked as visited, even though the edge was not calculated.
              processed[current_state->isurf][neighbor_i] = true;
              for(unsigned int i = 0; i < (*neighbor_searches[current_state->isurf]).get_size(); i++)
              {
                if((*neighbor_searches[current_state->isurf]).present(i))
                {
                  if(!(*neighbor_searches[current_state->isurf]).get(i)->neighbors.at(neighbor_i)->visited)
                  {
                    processed[current_state->isurf][neighbor_i] = false;
                    break;
                  }
                }
              }
            }
          }
        }
      }

      for(current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
      {
        if(intra_edge_passed_DG[current_state->isurf])
          continue;

        bool inner_edge_for_dg = false;
        for(int i = 0; i < this->spaces.size(); i++)
          if(current_state->e[i]->en[current_state->isurf]->marker == 0)
            inner_edge_for_dg = true;
        if(!inner_edge_for_dg)
          continue;

        for(unsigned int neighbor_i = 0; neighbor_i < num_neighbors[current_state->isurf]; neighbor_i++)
        {
          if(!DG_vector_forms_present && processed[current_state->isurf][neighbor_i])
            continue;

          // DG-inner-edge-wise parameters for WeakForm.
          (const_cast<WeakForm<Scalar>*>(current_wf))->set_active_DG_state(current_state->e, current_state->isurf);

          assemble_DG_one_neighbor(processed[current_state->isurf][neighbor_i], neighbor_i, current_pss, current_refmaps, current_u_ext, current_als,
            current_state, current_mfDG, current_vfDG, fn, npss, nspss, nrefmap, (*neighbor_searches[current_state->isurf]), min_dg_mesh_seq, current_wf);
        }

        // Delete the neighbor_searches array.
        for(unsigned int i = 0; i < (*neighbor_searches[current_state->isurf]).get_size(); i++)
          if((*neighbor_searches[current_state->isurf]).present(i))
            delete (*neighbor_searches[current_state->isurf]).get(i);
        delete neighbor_searches[current_state->isurf];
        delete [] processed[current_state->isurf];
      }

      delete [] processed;
      delete [] neighbor_searches;
      delete [] num_neighbors;

      // Deinitialize neighbor pss's, refmaps.
      if(DG_matrix_forms_present)
      {
        for(std::map<unsigned int, PrecalcShapeset *>::iterator it = nspss.begin(); it != nspss.end(); it++)
          delete it->second;
        for(std::map<unsigned int, PrecalcShapeset *>::iterator it = npss.begin(); it != npss.end(); it++)
          delete it->second;
        for(std::map<unsigned int, RefMap *>::iterator it = nrefmap.begin(); it != nrefmap.end(); it++)
          delete it->second;
      }
    }

    template<typename Scalar>
    bool DiscreteProblemDGAssembler<Scalar>::form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state)
    {
      if(current_state->e[form->i] && current_state->e[form->j])
      {
        if(fabs(form->scaling_factor) < 1e-12)
          return false;

        // If a block scaling table is provided, and if the scaling coefficient
        // A_mn for this block is zero, then the form does not need to be assembled.
        return true;
      }
      return false;
    }

    template<typename Scalar>
    bool DiscreteProblemDGAssembler<Scalar>::form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((MatrixForm<Scalar>*)form, current_state))
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    bool DiscreteProblemDGAssembler<Scalar>::form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((MatrixForm<Scalar>*)form, current_state))
        return false;

      if(current_state->rep->en[current_state->isurf]->marker == 0)
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->en[current_state->isurf]->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    bool DiscreteProblemDGAssembler<Scalar>::form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state)
    {
      if(!current_state->e[form->i])
        return false;
      if(fabs(form->scaling_factor) < 1e-12)
        return false;

      return true;
    }

    template<typename Scalar>
    bool DiscreteProblemDGAssembler<Scalar>::form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((VectorForm<Scalar>*)form, current_state))
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    bool DiscreteProblemDGAssembler<Scalar>::form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((VectorForm<Scalar>*)form, current_state))
        return false;

      if(current_state->rep->en[current_state->isurf]->marker == 0)
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->en[current_state->isurf]->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::assemble_DG_one_neighbor(bool edge_processed, unsigned int neighbor_i,
      PrecalcShapeset** current_pss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als,
      Traverse::State* current_state, Hermes::vector<MatrixFormDG<Scalar>*> current_mfDG, Hermes::vector<VectorFormDG<Scalar>*> current_vfDG, Transformable** fn,
      std::map<unsigned int, PrecalcShapeset *> npss, std::map<unsigned int, PrecalcShapeset *> nspss, std::map<unsigned int, RefMap *> nrefmap,
      LightArray<NeighborSearch<Scalar>*>& neighbor_searches, unsigned int min_dg_mesh_seq, WeakForm<Scalar>* current_wf)
    {
      // Set the active segment in all NeighborSearches
      for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
      {
        if(neighbor_searches.present(i))
        {
          NeighborSearch<Scalar>* ns = neighbor_searches.get(i);
          ns->active_segment = neighbor_i;
          ns->neighb_el = ns->neighbors[neighbor_i];
          ns->neighbor_edge = ns->neighbor_edges[neighbor_i];
        }
      }

      // Push all the necessary transformations to all functions of this stage.
      // The important thing is that the transformations to the current subelement are already there.
      for(unsigned int fns_i = 0; fns_i < current_state->num; fns_i++)
      {
        MeshSharedPtr mesh_i;
        if(fns_i < this->spaces.size())
          mesh_i = spaces[fns_i]->get_mesh();
        else
          mesh_i = (static_cast<MeshFunction<Scalar>* >(fn[fns_i]))->get_mesh();
        NeighborSearch<Scalar>* ns = neighbor_searches.get(mesh_i->get_seq() - min_dg_mesh_seq);
        if(ns->central_transformations[neighbor_i])
          ns->central_transformations[neighbor_i]->apply_on(fn[fns_i]);
      }

      // For neighbor psss.
      if(current_mat && DG_matrix_forms_present && !edge_processed)
      {
        for(unsigned int idx_i = 0; idx_i < spaces.size(); idx_i++)
        {
          NeighborSearch<Scalar>* ns = neighbor_searches.get(spaces[idx_i]->get_mesh()->get_seq() - min_dg_mesh_seq);
          npss[idx_i]->set_active_element((*ns->get_neighbors())[neighbor_i]);
          if(ns->neighbor_transformations[neighbor_i])
            ns->neighbor_transformations[neighbor_i]->apply_on(npss[idx_i]);
        }
      }

      // Also push the transformations to the slave psss and refmaps.
      for (unsigned int i = 0; i < spaces.size(); i++)
      {
        current_refmaps[i]->force_transform(current_pss[i]->get_transform(), current_pss[i]->get_ctm());

        // Neighbor.
        if(current_mat && DG_matrix_forms_present && !edge_processed)
        {
          nrefmap[i]->set_active_element(npss[i]->get_active_element());
          nrefmap[i]->force_transform(npss[i]->get_transform(), npss[i]->get_ctm());
        }
      }

      /***/
      // The computation takes place here.
      NeighborSearch<Scalar>** nbs = new NeighborSearch<Scalar>*[this->spaces.size()];
      typename NeighborSearch<Scalar>::ExtendedShapeset** ext_asmlist = new typename NeighborSearch<Scalar>::ExtendedShapeset*[this->spaces.size()];
      int n_quadrature_points;
      Geom<double>** geometry = new Geom<double>*[this->spaces.size()];
      double** jacobian_x_weights = new double*[this->spaces.size()];
      Geom<double>** e = new Geom<double>*[this->spaces.size()];
      DiscontinuousFunc<double>*** testFunctions = new DiscontinuousFunc<double>**[this->spaces.size()];

      // Create the extended shapeset on the union of the central element and its current neighbor.
      int order = 20;
      int order_base = 20;
      for (unsigned int i = 0; i < this->spaces.size(); i++)
      {
        if(this->spaces[i]->get_type() != HERMES_L2_SPACE)
          continue;

        nbs[i] = neighbor_searches.get(spaces[i]->get_mesh()->get_seq() - min_dg_mesh_seq);
        ext_asmlist[i] = nbs[i]->create_extended_asmlist(spaces[i], current_als[i]);
        nbs[i]->set_quad_order(order);
        order_base = order;
        n_quadrature_points = init_surface_geometry_points(current_refmaps[i], order_base, i, current_state->rep->marker, geometry[i], jacobian_x_weights[i]);
        e[i] = new InterfaceGeom<double>(geometry[i], nbs[i]->neighb_el->marker, nbs[i]->neighb_el->id, nbs[i]->neighb_el->get_diameter());

        testFunctions[i] = new DiscontinuousFunc<double>*[ext_asmlist[i]->cnt];
        for (int func_i = 0; func_i < ext_asmlist[i]->cnt; func_i++)
        {
          if(ext_asmlist[i]->dof[func_i] < 0)
            continue;

          // Choose the correct shapeset for the test function.
          if(!ext_asmlist[i]->has_support_on_neighbor(func_i))
          {
            current_pss[i]->set_active_shape(ext_asmlist[i]->central_al->idx[func_i]);
            PrecalcShapeset* func = current_pss[i];
            RefMap* refmap = current_refmaps[i];
            testFunctions[i][func_i] = new DiscontinuousFunc<double>(init_fn(func, refmap, nbs[i]->get_quad_eo(false)), false, nbs[i]->neighbor_edge.orientation);
          }
          else
          {
            npss[i]->set_active_shape(ext_asmlist[i]->neighbor_al->idx[func_i - ext_asmlist[i]->central_al->cnt]);
            PrecalcShapeset* func = npss[i];
            RefMap* refmap = nrefmap[i];
            testFunctions[i][func_i] = new DiscontinuousFunc<double>(init_fn(func, refmap, nbs[i]->get_quad_eo(true)), true, nbs[i]->neighbor_edge.orientation);
          }
        }
      }

      DiscontinuousFunc<Scalar>** ext = init_ext_fns(current_wf->ext, neighbor_searches, order, min_dg_mesh_seq);

      int prevNewtonSize = this->wf->get_neq();
      DiscontinuousFunc<Scalar>** u_ext = new DiscontinuousFunc<Scalar>*[prevNewtonSize];
      if(!this->is_linear)
      {
        if(current_u_ext)
          for(int u_ext_i = 0; u_ext_i < prevNewtonSize; u_ext_i++)
            if(current_u_ext[u_ext_i])
            {
              neighbor_searches.get(current_u_ext[u_ext_i]->get_mesh()->get_seq() - min_dg_mesh_seq)->set_quad_order(order);
              u_ext[u_ext_i]  = neighbor_searches.get(current_u_ext[u_ext_i]->get_mesh()->get_seq() - min_dg_mesh_seq)->init_ext_fn(current_u_ext[u_ext_i]);
            }
            else
              u_ext[u_ext_i] = NULL;
        else
          for(int u_ext_i = 0; u_ext_i < prevNewtonSize; u_ext_i++)
            u_ext[u_ext_i] = NULL;
      }



      if(current_mat && DG_matrix_forms_present && !edge_processed)
      {
        for(int current_mfsurf_i = 0; current_mfsurf_i < wf->mfDG.size(); current_mfsurf_i++)
        {
          if(!form_to_be_assembled((MatrixForm<Scalar>*)current_mfDG[current_mfsurf_i], current_state))
            continue;

          MatrixFormDG<Scalar>* mfs = current_mfDG[current_mfsurf_i];

          double block_scaling_coefficient = block_scaling_coeff(mfs);

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

                Scalar res = mfs->value(n_quadrature_points, jacobian_x_weights[n], u_ext, u, v, e[n], ext) * mfs->scaling_factor;

                support_neigh_u = ext_asmlist_u->has_support_on_neighbor(j);

                Scalar val = block_scaling_coefficient * 0.5 * res * (support_neigh_u ? ext_asmlist_u->neighbor_al->coef[j - ext_asmlist_u->central_al->cnt]: ext_asmlist_u->central_al->coef[j])
                  * (support_neigh_v ? ext_asmlist_v->neighbor_al->coef[i - ext_asmlist_v->central_al->cnt]: ext_asmlist_v->central_al->coef[i]);

                local_stiffness_matrix[i][j] = val;
              }
            }
          }

          current_mat->add(ext_asmlist_v->cnt, ext_asmlist_u->cnt, local_stiffness_matrix, ext_asmlist_v->dof, ext_asmlist_u->dof);

          delete [] local_stiffness_matrix;
        }
      }

      for(int i = 0; i < this->spaces.size(); i++)
      {
        if(this->spaces[i]->get_type() != HERMES_L2_SPACE)
          continue;
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

      delete [] testFunctions;
      delete [] ext_asmlist;

      if(current_rhs && DG_vector_forms_present)
      {
        for (unsigned int ww = 0; ww < wf->vfDG.size(); ww++)
        {
          VectorFormDG<Scalar>* vfs = current_vfDG[ww];
          if(vfs->areas[0] != H2D_DG_INNER_EDGE)
            continue;
          int n = vfs->i;

          if(!form_to_be_assembled((VectorForm<Scalar>*)vfs, current_state))
            continue;

          NeighborSearch<Scalar>* nbs_v = nbs[n];

          // Here we use the standard pss, possibly just transformed by NeighborSearch.
          for (unsigned int dof_i = 0; dof_i < current_als[n]->cnt; dof_i++)
          {
            if(current_als[n]->dof[dof_i] < 0)
              continue;
            current_pss[n]->set_active_shape(current_als[n]->idx[dof_i]);

            Func<double>* v = init_fn(current_pss[n], current_refmaps[n], nbs_v->get_quad_eo());

            current_rhs->add(current_als[n]->dof[dof_i], 0.5 * vfs->value(n_quadrature_points, jacobian_x_weights[n], u_ext, v, e[n], ext) * vfs->scaling_factor * current_als[n]->coef[dof_i]);

            v->free_fn();
            delete v;
          }
        }
      }

      if(ext)
      {
        for(unsigned int i = 0; i < current_wf->ext.size(); i++)
        {
          ext[i]->free_fn();
          delete ext[i];
        }
        delete [] ext;
      }

      if(!this->is_linear)
      {
        if(current_u_ext)
        {
          for(int u_ext_i = 0; u_ext_i < prevNewtonSize; u_ext_i++)
            if(current_u_ext[u_ext_i])
            {
              u_ext[u_ext_i]->free_fn();
              delete u_ext[u_ext_i];
            }
        }
      }

      delete [] u_ext;


      for(int i = 0; i < this->spaces.size(); i++)
      {
        if(this->spaces[i]->get_type() != HERMES_L2_SPACE)
          continue;
        delete [] jacobian_x_weights[i];
        e[i]->free();
        delete e[i];
      }

      delete [] nbs;
      delete [] geometry;
      delete [] jacobian_x_weights;
      delete [] e;

      // This is just cleaning after ourselves.
      // Clear the transformations from the RefMaps and all functions.
      for(unsigned int fns_i = 0; fns_i < current_state->num; fns_i++)
      {
        MeshSharedPtr mesh_i;
        if(fns_i < this->spaces.size())
          mesh_i = spaces[fns_i]->get_mesh();
        else
          mesh_i = (static_cast<MeshFunction<Scalar>* >(fn[fns_i]))->get_mesh();

        fn[fns_i]->set_transform(neighbor_searches.get(mesh_i->get_seq() - min_dg_mesh_seq)->original_central_el_transform);
      }

      // Also clear the transformations from the slave psss and refmaps.
      for (unsigned int i = 0; i < spaces.size(); i++)
        current_refmaps[i]->force_transform(current_pss[i]->get_transform(), current_pss[i]->get_ctm());
    }

    template<typename Scalar>
    DiscontinuousFunc<Scalar>** DiscreteProblemDGAssembler<Scalar>::init_ext_fns(Hermes::vector<MeshFunctionSharedPtr<Scalar> > ext,
      LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int order, unsigned int min_dg_mesh_seq)
    {
      DiscontinuousFunc<Scalar>** ext_fns = new DiscontinuousFunc<Scalar>*[ext.size()];
      for(unsigned int j = 0; j < ext.size(); j++)
      {
        neighbor_searches.get(ext[j]->get_mesh()->get_seq() - min_dg_mesh_seq)->set_quad_order(order);
        ext_fns[j] = neighbor_searches.get(ext[j]->get_mesh()->get_seq() - min_dg_mesh_seq)->init_ext_fn(ext[j].get());
      }

      return ext_fns;
    }

    template<typename Scalar>
    bool DiscreteProblemDGAssembler<Scalar>::init_neighbors(LightArray<NeighborSearch<Scalar>*>& neighbor_searches,
      Traverse::State* current_state, unsigned int min_dg_mesh_seq)
    {
      // Initialize the NeighborSearches.
      for(unsigned int i = 0; i < spaces.size(); i++)
      {
        if(i > 0 && spaces[i - 1]->get_mesh()->get_seq() == spaces[i]->get_mesh()->get_seq())
          continue;
        else
          if(!neighbor_searches.present(spaces[i]->get_mesh()->get_seq() - min_dg_mesh_seq))
          {
            NeighborSearch<Scalar>* ns = new NeighborSearch<Scalar>(current_state->e[i], spaces[i]->get_mesh());
            ns->original_central_el_transform = current_state->sub_idx[i];
            neighbor_searches.add(ns, spaces[i]->get_mesh()->get_seq() - min_dg_mesh_seq);
          }
      }

      // Calculate respective neighbors.
      // Also clear the initial_sub_idxs from the central element transformations
      // of NeighborSearches with multiple neighbors.
      // If all DG meshes have this edge as intra-edge, pass.
      bool DG_intra = false;
      for(unsigned int i = 0; i < spaces.size(); i++)
      {
        if(!(i > 0 && spaces[i]->get_mesh()->get_seq() - min_dg_mesh_seq == spaces[i-1]->get_mesh()->get_seq() - min_dg_mesh_seq))
        {
          if(neighbor_searches.get(spaces[i]->get_mesh()->get_seq() - min_dg_mesh_seq)->set_active_edge_multimesh(current_state->isurf) && spaces[i]->get_type() == HERMES_L2_SPACE)
            DG_intra = true;
          neighbor_searches.get(spaces[i]->get_mesh()->get_seq() - min_dg_mesh_seq)->clear_initial_sub_idx();
        }
      }
      return DG_intra;
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::build_multimesh_tree(MultimeshDGNeighborTreeNode* root,
      LightArray<NeighborSearch<Scalar>*>& neighbor_searches)
    {
      for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
        if(neighbor_searches.present(i))
        {
          NeighborSearch<Scalar>* ns = neighbor_searches.get(i);
          if(ns->n_neighbors == 1 &&
            (ns->central_transformations_size == 0 || ns->central_transformations[0]->num_levels == 0))
            continue;
          for(unsigned int j = 0; j < ns->n_neighbors; j++)
            if(ns->central_transformations[j])
              insert_into_multimesh_tree(root, ns->central_transformations[j]->transf, ns->central_transformations[j]->num_levels);
        }
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::insert_into_multimesh_tree(MultimeshDGNeighborTreeNode* node,
      unsigned int* transformations,
      unsigned int transformation_count)
    {
      // If we are already in the leaf.
      if(transformation_count == 0)
        return;
      // Both sons are null. We have to add a new Node. Let us do it for the left sone of node.
      if(node->get_left_son() == NULL && node->get_right_son() == NULL)
      {
        node->set_left_son(new MultimeshDGNeighborTreeNode(node, transformations[0]));
        insert_into_multimesh_tree(node->get_left_son(), transformations + 1, transformation_count - 1);
      }
      // At least the left son is not null (it is impossible only for the right one to be not null, because
      // the left one always gets into the tree first, as seen above).
      else
      {
        // The existing left son is the right one to continue through.
        if(node->get_left_son()->get_transformation() == transformations[0])
          insert_into_multimesh_tree(node->get_left_son(), transformations + 1, transformation_count - 1);
        // The right one also exists, check that it is the right one, or return an error.
        else if(node->get_right_son())
        {
          if(node->get_right_son()->get_transformation() == transformations[0])
            insert_into_multimesh_tree(node->get_right_son(), transformations + 1, transformation_count - 1);
          else
            throw Hermes::Exceptions::Exception("More than two possible sons in insert_into_multimesh_tree().");
        }
        // If the right one does not exist and the left one was not correct, create a right son and continue this way.
        else
        {
          node->set_right_son(new MultimeshDGNeighborTreeNode(node, transformations[0]));
          insert_into_multimesh_tree(node->get_right_son(), transformations + 1, transformation_count - 1);
        }
      }
    }

    template<typename Scalar>
    Hermes::vector<Hermes::vector<unsigned int>*> DiscreteProblemDGAssembler<Scalar>::get_multimesh_neighbors_transformations(MultimeshDGNeighborTreeNode* multimesh_tree)
    {
      // Initialize the vector.
      Hermes::vector<Hermes::vector<unsigned int>*> running_transformations;
      // Prepare the first neighbor's vector.
      running_transformations.push_back(new Hermes::vector<unsigned int>);
      // Fill the vector.
      traverse_multimesh_tree(multimesh_tree, running_transformations);
      return running_transformations;
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::traverse_multimesh_tree(MultimeshDGNeighborTreeNode* node,
      Hermes::vector<Hermes::vector<unsigned int>*>& running_transformations)
    {
      // If we are in the root.
      if(node->get_transformation() == 0)
      {
        if(node->get_left_son())
          traverse_multimesh_tree(node->get_left_son(), running_transformations);
        if(node->get_right_son())
          traverse_multimesh_tree(node->get_right_son(), running_transformations);
        // Delete the vector prepared by the last accessed leaf.
        delete running_transformations.back();
        running_transformations.pop_back();
        return;
      }
      // If we are in a leaf.
      if(node->get_left_son() == NULL && node->get_right_son() == NULL)
      {
        // Create a vector for the new neighbor.
        Hermes::vector<unsigned int>* new_neighbor_transformations = new Hermes::vector<unsigned int>;
        // Copy there the whole path except for this leaf.
        for(unsigned int i = 0; i < running_transformations.back()->size(); i++)
          new_neighbor_transformations->push_back((*running_transformations.back())[i]);
        // Insert this leaf into the current running transformation, thus complete it.
        running_transformations.back()->push_back(node->get_transformation());
        // Make the new_neighbor_transformations the current running transformation.
        running_transformations.push_back(new_neighbor_transformations);
        return;
      }
      else
      {
        running_transformations.back()->push_back(node->get_transformation());
        if(node->get_left_son())
          traverse_multimesh_tree(node->get_left_son(), running_transformations);
        if(node->get_right_son())
          traverse_multimesh_tree(node->get_right_son(), running_transformations);
        running_transformations.back()->pop_back();
        return;
      }
      return;
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::update_neighbor_search(NeighborSearch<Scalar>* ns, MultimeshDGNeighborTreeNode* multimesh_tree)
    {
      // This has to be done, because we pass ns by reference and the number of neighbors is changing.
      unsigned int num_neighbors = ns->get_num_neighbors();

      for(int i = 0; i < num_neighbors; i++)
      {
        // Find the node corresponding to this neighbor in the tree.
        MultimeshDGNeighborTreeNode* node;
        if(ns->central_transformations[i])
          node = find_node(ns->central_transformations[i]->transf, ns->central_transformations[i]->num_levels, multimesh_tree);
        else
          node = multimesh_tree;

        // Update the NeighborSearch.
        int added = update_ns_subtree(ns, node, i);
        i -= added;
        num_neighbors -= added;
      }
    }

    template<typename Scalar>
    MultimeshDGNeighborTreeNode* DiscreteProblemDGAssembler<Scalar>::find_node(unsigned int* transformations,
      unsigned int transformation_count,
      MultimeshDGNeighborTreeNode* node)
    {
      // If there are no transformations left.
      if(transformation_count == 0)
        return node;
      else
      {
        if(node->get_left_son())
        {
          if(node->get_left_son()->get_transformation() == transformations[0])
            return find_node(transformations + 1, transformation_count - 1, node->get_left_son());
        }
        if(node->get_right_son())
        {
          if(node->get_right_son()->get_transformation() == transformations[0])
            return find_node(transformations + 1, transformation_count - 1, node->get_right_son());
        }
      }
      // We always should be able to empty the transformations array.
      throw
        Hermes::Exceptions::Exception("Transformation of a central element not found in the multimesh tree.");
      return NULL;
    }

    template<typename Scalar>
    int DiscreteProblemDGAssembler<Scalar>::update_ns_subtree(NeighborSearch<Scalar>* ns,
      MultimeshDGNeighborTreeNode* node, unsigned int ith_neighbor)
    {
      int current_count = ns->get_num_neighbors();

      // No subtree => no work.
      // Also check the assertion that if one son is null, then the other too.
      if(node->get_left_son() == NULL)
      {
        if(node->get_right_son())
          throw Hermes::Exceptions::Exception("Only one son (right) not null in DiscreteProblemDGAssembler<Scalar>::update_ns_subtree.");
        return 0;
      }

      // Key part.
      // Begin with storing the info about the current neighbor.
      Element* neighbor = ns->neighbors[ith_neighbor];
      typename NeighborSearch<Scalar>::NeighborEdgeInfo edge_info = ns->neighbor_edges[ith_neighbor];

      // Initialize the vector for central transformations->
      Hermes::vector<Hermes::vector<unsigned int>*> running_central_transformations;
      // Prepare the first new neighbor's vector. Push back the current transformations (in case of GO_DOWN neighborhood).
      running_central_transformations.push_back(new Hermes::vector<unsigned int>);
      if(ns->central_transformations[ith_neighbor])
        ns->central_transformations[ith_neighbor]->copy_to(running_central_transformations.back());

      // Initialize the vector for neighbor transformations->
      Hermes::vector<Hermes::vector<unsigned int>*> running_neighbor_transformations;
      // Prepare the first new neighbor's vector. Push back the current transformations (in case of GO_UP/NO_TRF neighborhood).
      running_neighbor_transformations.push_back(new Hermes::vector<unsigned int>);
      if(ns->neighbor_transformations[ith_neighbor])
        ns->neighbor_transformations[ith_neighbor]->copy_to(running_neighbor_transformations.back());

      // Delete the current neighbor.
      ns->delete_neighbor(ith_neighbor);

      // Move down the subtree.
      if(node->get_left_son())
        traverse_multimesh_subtree(node->get_left_son(), running_central_transformations,
        running_neighbor_transformations, edge_info, ns->active_edge,
        ns->central_el->get_mode());
      if(node->get_right_son())
        traverse_multimesh_subtree(node->get_right_son(), running_central_transformations,
        running_neighbor_transformations, edge_info, ns->active_edge,
        ns->central_el->get_mode());

      // Delete the last neighbors' info (this is a dead end, caused by the function traverse_multimesh_subtree.
      delete running_central_transformations.back();
      running_central_transformations.pop_back();
      delete running_neighbor_transformations.back();
      running_neighbor_transformations.pop_back();

      // Insert new neighbors.
      for(unsigned int i = 0; i < running_central_transformations.size(); i++)
      {
        ns->neighbors.push_back(neighbor);
        ns->neighbor_edges.push_back(edge_info);

        if(!ns->central_transformations[ns->n_neighbors])
          ns->add_central_transformations(new typename NeighborSearch<Scalar>::Transformations, ns->n_neighbors);

        if(!ns->neighbor_transformations[ns->n_neighbors])
          ns->add_neighbor_transformations(new typename NeighborSearch<Scalar>::Transformations, ns->n_neighbors);

        ns->central_transformations[ns->n_neighbors]->copy_from(*running_central_transformations[i]);
        ns->neighbor_transformations[ns->n_neighbors]->copy_from(*running_neighbor_transformations[i]);

        ns->n_neighbors++;
      }

      for(unsigned int i = 0; i < running_central_transformations.size(); i++)
        delete running_central_transformations[i];
      for(unsigned int i = 0; i < running_neighbor_transformations.size(); i++)
        delete running_neighbor_transformations[i];

      // Return the number of neighbors added/deleted.
      return ns->get_num_neighbors() - current_count;
    }

    template<typename Scalar>
    void DiscreteProblemDGAssembler<Scalar>::traverse_multimesh_subtree(MultimeshDGNeighborTreeNode* node,
      Hermes::vector<Hermes::vector<unsigned int>*>& running_central_transformations,
      Hermes::vector<Hermes::vector<unsigned int>*>& running_neighbor_transformations,
      const typename NeighborSearch<Scalar>::NeighborEdgeInfo& edge_info, const int& active_edge, const int& mode)
    {
      // If we are in a leaf.
      if(node->get_left_son() == NULL && node->get_right_son() == NULL)
      {
        // Create vectors for the new neighbor.
        Hermes::vector<unsigned int>* new_neighbor_central_transformations = new Hermes::vector<unsigned int>;
        Hermes::vector<unsigned int>* new_neighbor_neighbor_transformations = new Hermes::vector<unsigned int>;

        // Copy there the whole path except for this leaf.
        for(unsigned int i = 0; i < running_central_transformations.back()->size(); i++)
          new_neighbor_central_transformations->push_back((*running_central_transformations.back())[i]);
        for(unsigned int i = 0; i < running_neighbor_transformations.back()->size(); i++)
          new_neighbor_neighbor_transformations->push_back((*running_neighbor_transformations.back())[i]);

        // Insert this leaf into the current running central transformation, thus complete it.
        running_central_transformations.back()->push_back(node->get_transformation());

        // Make the new_neighbor_central_transformations the current running central transformation.
        running_central_transformations.push_back(new_neighbor_central_transformations);

        // Take care of the neighbor transformation.
        // Insert appropriate info from this leaf into the current running neighbor transformation, thus complete it.
        if(mode == HERMES_MODE_TRIANGLE)
          if((active_edge == 0 && node->get_transformation() == 0) ||
            (active_edge == 1 && node->get_transformation() == 1) ||
            (active_edge == 2 && node->get_transformation() == 2))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
        // Quads.
        else
          if((active_edge == 0 && (node->get_transformation() == 0 || node->get_transformation() == 6)) ||
            (active_edge == 1 && (node->get_transformation() == 1 || node->get_transformation() == 4)) ||
            (active_edge == 2 && (node->get_transformation() == 2 || node->get_transformation() == 7)) ||
            (active_edge == 3 && (node->get_transformation() == 3 || node->get_transformation() == 5)))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % H2D_MAX_NUMBER_EDGES));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % H2D_MAX_NUMBER_EDGES));

        // Make the new_neighbor_neighbor_transformations the current running neighbor transformation.
        running_neighbor_transformations.push_back(new_neighbor_neighbor_transformations);

        return;
      }
      else
      {
        // Insert this leaf into the current running central transformation, thus complete it.
        running_central_transformations.back()->push_back(node->get_transformation());

        // Insert appropriate info from this leaf into the current running neighbor transformation, thus complete it.
        // Triangles.
        if(mode == HERMES_MODE_TRIANGLE)
          if((active_edge == 0 && node->get_transformation() == 0) ||
            (active_edge == 1 && node->get_transformation() == 1) ||
            (active_edge == 2 && node->get_transformation() == 2))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
        // Quads.
        else
          if((active_edge == 0 && (node->get_transformation() == 0 || node->get_transformation() == 6)) ||
            (active_edge == 1 && (node->get_transformation() == 1 || node->get_transformation() == 4)) ||
            (active_edge == 2 && (node->get_transformation() == 2 || node->get_transformation() == 7)) ||
            (active_edge == 3 && (node->get_transformation() == 3 || node->get_transformation() == 5)))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % H2D_MAX_NUMBER_EDGES));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % H2D_MAX_NUMBER_EDGES));

        // Move down.
        if(node->get_left_son())
          traverse_multimesh_subtree(node->get_left_son(), running_central_transformations, running_neighbor_transformations,
          edge_info, active_edge, mode);
        if(node->get_right_son())
          traverse_multimesh_subtree(node->get_right_son(), running_central_transformations, running_neighbor_transformations,
          edge_info, active_edge, mode);

        // Take this transformation out.
        running_central_transformations.back()->pop_back();
        running_neighbor_transformations.back()->pop_back();
        return;
      }
      return;
    }

    template class HERMES_API DiscreteProblemDGAssembler<double>;
    template class HERMES_API DiscreteProblemDGAssembler<std::complex<double> >;
  }
}