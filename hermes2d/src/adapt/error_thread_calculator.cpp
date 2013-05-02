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

#include "adapt/error_thread_calculator.h"
#include "discrete_problem/discrete_problem_helpers.h"
#include "discrete_problem/dg/multimesh_dg_neighbor_tree.h"
#include "neighbor_search.h"
#include "mesh/refmap.h"
#include "forms.h"
#include <stdlib.h>

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ErrorThreadCalculator<Scalar>::ErrorThreadCalculator(ErrorCalculator<Scalar>* errorCalculator) :
      errorCalculator(errorCalculator)
    {
      slns = new Solution<Scalar>*[this->errorCalculator->component_count];
      rslns = new Solution<Scalar>*[this->errorCalculator->component_count];

      for (int j = 0; j < this->errorCalculator->component_count; j++)
      {
        slns[j] = static_cast<Solution<Scalar>* >(errorCalculator->coarse_solutions[j]->clone());
        rslns[j] = static_cast<Solution<Scalar>* >(errorCalculator->fine_solutions[j]->clone());
      }
    }

    template<typename Scalar>
    ErrorThreadCalculator<Scalar>::~ErrorThreadCalculator()
    {
      for (int j = 0; j < this->errorCalculator->component_count; j++)
      {
        delete slns[j];
        delete rslns[j];
      }

      delete [] slns;
      delete [] rslns;
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::evaluate_one_state(Traverse::State* current_state_)
    {
      this->current_state = current_state_;

      // Initialization.
      for(int i = 0; i < this->errorCalculator->component_count; i++)
      {
        slns[i]->set_active_element(current_state->e[i]);
        slns[i]->set_transform(current_state->sub_idx[i]);

        rslns[i]->set_active_element(current_state->e[this->errorCalculator->component_count + i]);
        rslns[i]->set_transform(current_state->sub_idx[this->errorCalculator->component_count + i]);
      }

      // Max order imposement.
      int order = g_quad_2d_std.get_max_order(current_state->rep->get_mode());

      // Volumetric forms.
      this->evaluate_volumetric_forms(current_state, order);

      // Surface forms.
      if(current_state->isBnd && this->errorCalculator->mfsurf.size() > 0)
      {
        for (current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
        {
          if(!current_state->bnd[current_state->isurf])
            continue;

          this->evaluate_surface_forms_one_edge(current_state, order);
        }
      }

      // DG forms.
      if(this->errorCalculator->mfDG.size() > 0)
      {
        for (current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
        {
          if(current_state->bnd[current_state->isurf])
            continue;

          ErrorThreadCalculator<Scalar>::DGErrorCalculator dGErrorCalculator(this);
          dGErrorCalculator.assemble_one_edge();
        }
      }
    }

    template<typename Scalar>
    ErrorThreadCalculator<Scalar>::DGErrorCalculator::DGErrorCalculator(ErrorThreadCalculator* errorThreadCalculator) : errorThreadCalculator(errorThreadCalculator), current_state(errorThreadCalculator->current_state)
    {
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::DGErrorCalculator::assemble_one_edge()
    {
      this->neighbor_searches = new NeighborSearch<Scalar>*[this->current_state->num];
      this->processed = NULL;

#pragma omp critical (DG)
      {
        for(unsigned int i = 0; i < current_state->num; i++)
          current_state->e[i]->visited = true;

        // If this edge is an inter-element one on all meshes.
        if(init_neighbors())
        {

          // Create a multimesh tree;
          MultimeshDGNeighborTree<Scalar>::process_edge(this->neighbor_searches, this->current_state->num, this->num_neighbors, this->processed);

          for(unsigned int neighbor_i = 0; neighbor_i < num_neighbors; neighbor_i++)
          {
            if(processed[neighbor_i])
              continue;

            this->assemble_one_neighbor(processed[neighbor_i], neighbor_i);
          }
        }
      }

      delete [] neighbor_searches;
      if(processed)
        delete [] processed;
    }

    template<typename Scalar>
    bool ErrorThreadCalculator<Scalar>::DGErrorCalculator::init_neighbors()
    {
      // Initialize the NeighborSearches.
      bool DG_intra = false;
      for(unsigned int i = 0; i < current_state->num; i++)
      {
        bool existing_ns = false;
        for(int j = i - 1; j >= 0; j--)
          if(current_state->e[i] == current_state->e[j])
          {
            neighbor_searches[i] = neighbor_searches[j];
            existing_ns = true;
            break;
          }
        if(!existing_ns)
        {
          NeighborSearch<Scalar>* ns;
          if(i < this->errorThreadCalculator->errorCalculator->component_count)
            ns = new NeighborSearch<Scalar>(current_state->e[i], this->errorThreadCalculator->slns[i]->get_mesh());
          else
            ns = new NeighborSearch<Scalar>(current_state->e[i], this->errorThreadCalculator->rslns[i - this->errorThreadCalculator->errorCalculator->component_count]->get_mesh());

          ns->original_central_el_transform = current_state->sub_idx[i];
          neighbor_searches[i] = ns;
          if(neighbor_searches[i]->set_active_edge_multimesh(current_state->isurf))
            DG_intra = true;
          neighbor_searches[i]->clear_initial_sub_idx();
        }
      }

      return DG_intra;
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::DGErrorCalculator::assemble_one_neighbor(bool edge_processed, unsigned int neighbor_i)
    {
      if(edge_processed)
        return;

      // Set the active segment in all NeighborSearches
      for(unsigned int i = 0; i < this->current_state->num; i++)
      {
        NeighborSearch<Scalar>* ns = neighbor_searches[i];
        ns->active_segment = neighbor_i;
        ns->neighb_el = ns->neighbors[neighbor_i];
        ns->neighbor_edge = ns->neighbor_edges[neighbor_i];
      }

      // Push all the necessary transformations to all functions of this stage.
      // The important thing is that the transformations to the current subelement are already there.
      for(unsigned int fns_i = 0; fns_i < this->errorThreadCalculator->errorCalculator->component_count; fns_i++)
      {
        NeighborSearch<Scalar>* ns = neighbor_searches[fns_i];
        if(ns->central_transformations[neighbor_i])
          ns->central_transformations[neighbor_i]->apply_on(this->errorThreadCalculator->slns[fns_i]);

        ns = neighbor_searches[fns_i + this->errorThreadCalculator->errorCalculator->component_count];
        if(ns->central_transformations[neighbor_i])
          ns->central_transformations[neighbor_i]->apply_on(this->errorThreadCalculator->rslns[fns_i]);
      }

      /***/
      // The computation takes place here.
      // Create the extended shapeset on the union of the central element and its current neighbor.
      int order = 20;
      int order_base = 20;

      for(int i = 0; i < this->errorThreadCalculator->errorCalculator->component_count; i++)
      {
        neighbor_searches[i]->set_quad_order(order);
        neighbor_searches[i + this->errorThreadCalculator->errorCalculator->component_count]->set_quad_order(order);
      }
      order_base = order;

      RefMap** refmaps = new RefMap*[this->errorThreadCalculator->errorCalculator->component_count];
      for(int i = 0; i < this->errorThreadCalculator->errorCalculator->component_count; i++)
        refmaps[i] = this->errorThreadCalculator->slns[i]->get_refmap();
      this->errorThreadCalculator->n_quadrature_points = init_surface_geometry_points(refmaps, this->errorThreadCalculator->errorCalculator->component_count, order_base, current_state->isurf, current_state->rep->marker, this->errorThreadCalculator->geometry, this->errorThreadCalculator->jacobian_x_weights);
      delete [] refmaps;

      DiscontinuousFunc<Scalar>** difference_funcs = new DiscontinuousFunc<Scalar>*[this->errorThreadCalculator->errorCalculator->component_count];
      DiscontinuousFunc<Scalar>** rsln_funcs = new DiscontinuousFunc<Scalar>*[this->errorThreadCalculator->errorCalculator->component_count];
      for(int i = 0; i < this->errorThreadCalculator->errorCalculator->component_count; i++)
      {
        difference_funcs[i] = neighbor_searches[i]->init_ext_fn(this->errorThreadCalculator->slns[i]);
        rsln_funcs[i] = neighbor_searches[i + this->errorThreadCalculator->errorCalculator->component_count]->init_ext_fn(this->errorThreadCalculator->rslns[i]);
        difference_funcs[i]->subtract(*rsln_funcs[i]);
      }

      for(int current_mfDG_i = 0; current_mfDG_i < this->errorThreadCalculator->errorCalculator->mfDG.size(); current_mfDG_i++)
      {
        NormFormDG<Scalar>* mfs = this->errorThreadCalculator->errorCalculator->mfDG[current_mfDG_i];

        double* error = &this->errorThreadCalculator->errorCalculator->errors[mfs->i][current_state->e[mfs->i]->id];
        double* norm = &this->errorThreadCalculator->errorCalculator->norms[mfs->i][current_state->e[mfs->i]->id];
        
        this->errorThreadCalculator->evaluate_DG_form(mfs, difference_funcs[mfs->i], difference_funcs[mfs->j], rsln_funcs[mfs->i], rsln_funcs[mfs->j], error, norm);
      }

      // deinitialize Funcs
      for(int i = 0; i < this->errorThreadCalculator->errorCalculator->component_count; i++)
      {
        difference_funcs[i]->free_fn();
        delete difference_funcs[i];
        rsln_funcs[i]->free_fn();
        delete rsln_funcs[i];
      }

      delete [] this->errorThreadCalculator->jacobian_x_weights;
      this->errorThreadCalculator->geometry->free();
      delete this->errorThreadCalculator->geometry;

      // This is just cleaning after ourselves.
      // Clear the transformations from the RefMaps and all functions.
      for(unsigned int fns_i = 0; fns_i < this->errorThreadCalculator->errorCalculator->component_count; fns_i++)
      {
        this->errorThreadCalculator->slns[fns_i]->set_transform(neighbor_searches[fns_i]->original_central_el_transform);
        this->errorThreadCalculator->rslns[fns_i]->set_transform(neighbor_searches[this->errorThreadCalculator->errorCalculator->component_count + fns_i]->original_central_el_transform);
      }
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::evaluate_volumetric_forms(Traverse::State* current_state, int order)
    {
      // initialize points & geometry & jacobian times weights.
      RefMap** refmaps = new RefMap*[this->errorCalculator->component_count];
      for(int i = 0; i < this->errorCalculator->component_count; i++)
        refmaps[i] = slns[i]->get_refmap();
      this->n_quadrature_points = init_geometry_points(refmaps, this->errorCalculator->component_count, order, this->geometry, this->jacobian_x_weights);
      delete [] refmaps;
      
      // initialize Funcs
      Func<Scalar>** difference_funcs = new Func<Scalar>*[this->errorCalculator->component_count];
      Func<Scalar>** rsln_funcs = new Func<Scalar>*[this->errorCalculator->component_count];
      for(int i = 0; i < this->errorCalculator->component_count; i++)
      {
        difference_funcs[i] = init_fn(slns[i], order);
        rsln_funcs[i] = init_fn(rslns[i], order);
        difference_funcs[i]->subtract(rsln_funcs[i]);
      }

      for(int i = 0; i < this->errorCalculator->mfvol.size(); i++)
      {
        NormFormVol<Scalar>* form = this->errorCalculator->mfvol[i];

        double* error = &this->errorCalculator->errors[form->i][current_state->e[form->i]->id];
        double* norm = &this->errorCalculator->norms[form->i][current_state->e[form->i]->id];

        this->evaluate_volumetric_form(form, difference_funcs[form->i], difference_funcs[form->j], rsln_funcs[form->i], rsln_funcs[form->j], error, norm);
      }

      // deinitialize Funcs
      for(int i = 0; i < this->errorCalculator->component_count; i++)
      {
        difference_funcs[i]->free_fn();
        delete difference_funcs[i];
        rsln_funcs[i]->free_fn();
        delete rsln_funcs[i];
      }
      delete [] difference_funcs;
      delete [] rsln_funcs;

      // deinitialize points & geometry & jacobian times weights
      geometry->free();
      delete geometry;
      delete [] this->jacobian_x_weights;
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::evaluate_surface_forms_one_edge(Traverse::State* current_state, int order)
    {
      RefMap** refmaps = new RefMap*[this->errorCalculator->component_count];
      for(int i = 0; i < this->errorCalculator->component_count; i++)
        refmaps[i] = slns[i]->get_refmap();
      this->n_quadrature_points = init_surface_geometry_points(refmaps, this->errorCalculator->component_count, order, current_state->isurf, current_state->rep->marker, this->geometry, this->jacobian_x_weights);
      delete [] refmaps;

      // initialize Funcs
      Func<Scalar>** difference_funcs = new Func<Scalar>*[this->errorCalculator->component_count];
      Func<Scalar>** rsln_funcs = new Func<Scalar>*[this->errorCalculator->component_count];
      for(int i = 0; i < this->errorCalculator->component_count; i++)
      {
        difference_funcs[i] = init_fn(slns[i], order);
        rsln_funcs[i] = init_fn(rslns[i], order);
        difference_funcs[i]->subtract(rsln_funcs[i]);
      }

      for(int i = 0; i < this->errorCalculator->mfsurf.size(); i++)
      {
        NormFormSurf<Scalar>* form = this->errorCalculator->mfsurf[i];
     
        double* error = &this->errorCalculator->errors[form->i][current_state->e[form->i]->id];
        double* norm = &this->errorCalculator->norms[form->i][current_state->e[form->i]->id];

        this->evaluate_surface_form(form, difference_funcs[form->i], difference_funcs[form->j], rsln_funcs[form->i], rsln_funcs[form->j], error, norm);
      }

      // deinitialize Funcs
      for(int i = 0; i < this->errorCalculator->component_count; i++)
      {
        difference_funcs[i]->free_fn();
        delete difference_funcs[i];
        rsln_funcs[i]->free_fn();
        delete rsln_funcs[i];
      }
      delete [] difference_funcs;
      delete [] rsln_funcs;

      // deinitialize points & geometry & jacobian times weights
      geometry->free();
      delete geometry;
      delete [] this->jacobian_x_weights;
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::evaluate_volumetric_form(NormFormVol<Scalar>* form, Func<Scalar>* difference_func_i, Func<Scalar>* difference_func_j, Func<Scalar>* rsln_i, Func<Scalar>* rsln_j, double* error, double* norm)
    {
      double error_value = std::abs(form->value(this->n_quadrature_points, this->jacobian_x_weights, difference_func_i, difference_func_j, this->geometry));
#pragma omp atomic
      (*error) += error_value;

      double norm_value = std::abs(form->value(this->n_quadrature_points, this->jacobian_x_weights, rsln_i, rsln_j, this->geometry));
        
#pragma omp atomic
      (*norm) += norm_value;
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::evaluate_surface_form(NormFormSurf<Scalar>* form, Func<Scalar>* difference_func_i, Func<Scalar>* difference_func_j, Func<Scalar>* rsln_i, Func<Scalar>* rsln_j, double* error, double* norm)
    {
      double error_value = std::abs(form->value(this->n_quadrature_points, this->jacobian_x_weights, difference_func_i, difference_func_j, this->geometry));
      
      // 1D quadrature has the weights summed to 2.
      error_value *= 0.5;

#pragma omp atomic
      (*error) += error_value;

      double norm_value = std::abs(form->value(this->n_quadrature_points, this->jacobian_x_weights, rsln_i, rsln_j, this->geometry));

      // 1D quadrature has the weights summed to 2.
      norm_value *= 0.5;

#pragma omp atomic
      (*norm) += norm_value;
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::evaluate_DG_form(NormFormDG<Scalar>* form, DiscontinuousFunc<Scalar>* difference_func_i, DiscontinuousFunc<Scalar>* difference_func_j, DiscontinuousFunc<Scalar>* rsln_i, DiscontinuousFunc<Scalar>* rsln_j, double* error, double* norm)
    {
      double error_value = std::abs(form->value(this->n_quadrature_points, this->jacobian_x_weights, difference_func_i, difference_func_j, this->geometry));

      // 1D quadrature has the weights summed to 2.
      error_value *= 0.5;

#pragma omp atomic
      (*error) += error_value;

      double norm_value = std::abs(form->value(this->n_quadrature_points, this->jacobian_x_weights, rsln_i, rsln_j, this->geometry));

#pragma omp atomic
      (*norm) += norm_value;
    }

    template HERMES_API class ErrorThreadCalculator<double>;
    template HERMES_API class ErrorThreadCalculator<std::complex<double> >;
  }
}