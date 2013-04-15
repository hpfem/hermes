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
#include <stdlib.h>

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ErrorThreadCalculator<Scalar>::ErrorThreadCalculator(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& coarse_solutions, Hermes::vector<MeshFunctionSharedPtr<Scalar> >& fine_solutions, ErrorCalculator<Scalar>* errorCalculator) :
      errorCalculator(errorCalculator)
    {
      slns = new Solution<Scalar>*[this->errorCalculator->component_count];
      rslns = new Solution<Scalar>*[this->errorCalculator->component_count];

      for (int j = 0; j < this->errorCalculator->component_count; j++)
      {
        slns[j] = static_cast<Solution<Scalar>* >(coarse_solutions[j]->clone());
        rslns[j] = static_cast<Solution<Scalar>* >(fine_solutions[j]->clone());
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
    void ErrorThreadCalculator<Scalar>::evaluate_one_state(Traverse::State* current_state)
    {
      int running_count = 0;
      for(int i = 0; i < this->errorCalculator->component_count; i++)
      {
        slns[i]->set_active_element(current_state->e[i]);
        slns[i]->set_transform(current_state->sub_idx[i]);

        rslns[i]->set_active_element(current_state->e[this->errorCalculator->component_count + i]);
        rslns[i]->set_transform(current_state->sub_idx[this->errorCalculator->component_count + i]);

        if(!this->errorCalculator->element_references[running_count + current_state->e[i]->id])
          this->errorCalculator->element_references[running_count + current_state->e[i]->id] = new ElementReference(i, current_state->e[i]->id, &this->errorCalculator->errors[i][current_state->e[i]->id]);
        running_count += slns[i]->get_mesh()->get_num_active_elements();
      }

      int order = g_quad_2d_std.get_max_order(current_state->rep->get_mode());

      this->evaluate_volumetric_forms(current_state, order);

      if(current_state->isBnd && this->errorCalculator->mfsurf.size() > 0)
      {
        for (current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
        {
          if(!current_state->bnd[current_state->isurf])
            continue;

          this->evaluate_surface_forms_one_edge(current_state, order);
        }
      }

      if(this->errorCalculator->mfDG.size() > 0)
      {
        for (current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
        {
          if(current_state->bnd[current_state->isurf])
            continue;

          this->evaluate_DG_forms_one_edge(current_state, order);
        }
      }
    }

    template<typename Scalar>
    ErrorThreadCalculator<Scalar>::DGErrorCalculator::DGErrorCalculator(ErrorThreadCalculator* errorThreadCalculator) : DiscreteProblemDGAssembler<Scalar>(errorThreadCalculator->errorCalculator->mfDG), errorThreadCalculator(errorThreadCalculator)
    {
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::DGErrorCalculator::assemble_one_neighbor(bool edge_processed, unsigned int neighbor_i, NeighborSearch<Scalar>** current_neighbor_searches)
    {
      if(edge_processed)
        return;

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
      for(unsigned int fns_i = 0; fns_i < this->errorThreadCalculator->errorCalculator->component_count; fns_i++)
      {
        NeighborSearch<Scalar>* ns = current_neighbor_searches[fns_i];
        if(ns->central_transformations[neighbor_i])
          ns->central_transformations[neighbor_i]->apply_on(this->errorThreadCalculator->slns[fns_i]);

        ns = current_neighbor_searches[fns_i + this->errorThreadCalculator->errorCalculator->component_count];
        if(ns->central_transformations[neighbor_i])
          ns->central_transformations[neighbor_i]->apply_on(this->errorThreadCalculator->rslns[fns_i]);
      }

      /***/
      // The computation takes place here.
      int n_quadrature_points;
      Geom<double>* geometry;
      double* jacobian_x_weights;
      Geom<double>* e;

      // Create the extended shapeset on the union of the central element and its current neighbor.
      int order = 20;
      int order_base = 20;

      for(int i = 0; i < this->errorThreadCalculator->errorCalculator->component_count; i++)
        current_neighbor_searches[i]->set_quad_order(order);
      order_base = order;
      n_quadrature_points = init_surface_geometry_points(this->errorThreadCalculator->rslns[0]->get_refmap(), order_base, current_state->isurf, current_state->rep->marker, geometry, jacobian_x_weights);
      e = new InterfaceGeom<double>(geometry[i], current_neighbor_searches[0]->neighb_el->marker, current_neighbor_searches[0]->neighb_el->id, current_neighbor_searches[0]->neighb_el->get_diameter());

      DiscontinuousFunc<Scalar>** difference_funcs = new DiscontinuousFunc<Scalar>*[this->errorThreadCalculator->errorCalculator->component_count];
      DiscontinuousFunc<Scalar>** rsln_funcs = new DiscontinuousFunc<Scalar>*[this->errorThreadCalculator->errorCalculator->component_count];
      for(int i = 0; i < this->errorThreadCalculator->errorCalculator->component_count; i++)
      {
        difference_funcs[i] = current_neighbor_searches[i]->init_ext_fn(slns[i]);
        rsln_funcs[i] = current_neighbor_searches[i + this->errorThreadCalculator->errorCalculator->component_count]->init_ext_fn(rslns[i]);
        difference_funcs[i]->subtract(rsln_funcs[i]);
      }

      for(int current_mfDG_i = 0; current_mfDG_i < this->mfDG.size(); current_mfDG_i++)
      {
        MatrixFormDG<Scalar>* mfs = this->mfDG[current_mfDG_i];

        this->evaluate_form(form, difference_funcs[form->i], difference_funcs[form->j], rsln_funcs[form->i], rsln_funcs[form->j], &this->errors[form->i][current_state->e[form->i]->id], (this->strategy == RelativeError ? &this->norms[form->i][current_state->e[form->i]->id] : NULL));
      }

      // deinitialize Funcs
      for(int i = 0; i < this->errorThreadCalculator->errorCalculator->component_count; i++)
      {
        difference_funcs[i]->free_fn();
        delete difference_funcs[i];
        rsln_funcs[i]->free_fn();
        delete rsln_funcs[i];
      }

      delete [] jacobian_x_weights;
      e->free();
      delete e;

      // This is just cleaning after ourselves.
      // Clear the transformations from the RefMaps and all functions.
      for(unsigned int fns_i = 0; fns_i < this->errorThreadCalculator->errorCalculator->component_count; fns_i++)
      {
        this->errorThreadCalculator->slns[fns_i]->set_transform(current_neighbor_searches[fns_i]->original_central_el_transform);
        this->errorThreadCalculator->rslns[fns_i]->set_transform(current_neighbor_searches[this->errorThreadCalculator->errorCalculator->component_count + fns_i]->original_central_el_transform);
      }
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::evaluate_volumetric_forms(Traverse::State* current_state, int order)
    {
      // initialize points & geometry & jacobian times weights
      this->n_quadrature_points = init_geometry_points(rslns[0]->get_refmap(), order, this->geometry, this->jacobian_x_weights);

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
        MatrixFormVol<Scalar>* form = this->errorCalculator->mfvol[i];
        this->evaluate_form(form, difference_funcs[form->i], difference_funcs[form->j], rsln_funcs[form->i], rsln_funcs[form->j], &this->errorCalculator->errors[form->i][current_state->e[form->i]->id], (this->errorCalculator->strategy == RelativeError ? &this->errorCalculator->norms[form->i][current_state->e[form->i]->id] : NULL));
      }

      // deinitialize Funcs
      for(int i = 0; i < this->errorCalculator->component_count; i++)
      {
        difference_funcs[i]->free_fn();
        delete difference_funcs[i];
        rsln_funcs[i]->free_fn();
        delete rsln_funcs[i];
      }

      // deinitialize points & geometry & jacobian times weights
      geometry->free();
      delete geometry;
      delete [] this->jacobian_x_weights;
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::evaluate_surface_forms_one_edge(Traverse::State* current_state, int order)
    {
      this->n_quadrature_points = init_surface_geometry_points(rslns[0]->get_refmap(), order, current_state->isurf, current_state->rep->marker, this->geometry, this->jacobian_x_weights);

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
        MatrixFormSurf<Scalar>* form = this->errorCalculator->mfsurf[i];
        this->evaluate_form(form, difference_funcs[form->i], difference_funcs[form->j], rsln_funcs[form->i], rsln_funcs[form->j], &this->errorCalculator->errors[form->i][current_state->e[form->i]->id], (this->errorCalculator->strategy == RelativeError ? &this->errorCalculator->norms[form->i][current_state->e[form->i]->id] : NULL));
      }

      // deinitialize Funcs
      for(int i = 0; i < this->errorCalculator->component_count; i++)
      {
        difference_funcs[i]->free_fn();
        delete difference_funcs[i];
        rsln_funcs[i]->free_fn();
        delete rsln_funcs[i];
      }

      // deinitialize points & geometry & jacobian times weights
      geometry->free();
      delete geometry;
      delete [] this->jacobian_x_weights;
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::evaluate_DG_forms_one_edge(Traverse::State* current_state, int order)
    {
      ErrorThreadCalculator<Scalar>::DGErrorCalculator dGErrorCalculator(this);
      dGErrorCalculator.init_assembling_one_state(current_state);
      dGErrorCalculator.assemble_one_state();
      dGErrorCalculator.deinit_assembling_one_state();
    }

    template<typename Scalar>
    void ErrorThreadCalculator<Scalar>::evaluate_form(MatrixForm<Scalar>* form, Func<Scalar>* difference_func_i, Func<Scalar>* difference_func_j, Func<Scalar>* rsln_i, Func<Scalar>* rsln_j, double* error, double* norm)
    {
      bool volumetric_form = (dynamic_cast<MatrixFormVol<Scalar>*>(form));

      // Calculation
      Scalar error_value = (volumetric_form ? 1.0 : 0.5) * std::abs(form->value(this->n_quadrature_points, this->jacobian_x_weights, NULL, difference_func_i, difference_func_j, this->geometry, NULL));

#pragma omp atomic
      (*error) += error_value;

      if(norm)
      {
        Scalar norm_value = std::abs(form->value(this->n_quadrature_points, this->jacobian_x_weights, NULL, rsln_i, rsln_j, this->geometry, NULL));

#pragma omp atomic
        (*norm) += norm_value;
      }
    }
  }
}