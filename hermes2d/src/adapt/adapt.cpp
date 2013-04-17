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

#include "adapt.h"
#include "projections/ogprojection.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    Adapt<Scalar>::Adapt(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces_, ErrorCalculator<Scalar>* errorCalculator) : errorCalculator(errorCalculator)
    {
      for(unsigned int i = 0; i < spaces_.size(); i++)
        spaces.push_back(spaces_[i]);

      this->init();
    }

    template<typename Scalar>
    Adapt<Scalar>::Adapt(SpaceSharedPtr<Scalar>& space, ErrorCalculator<Scalar>* errorCalculator) : errorCalculator(errorCalculator)
    {
      spaces.push_back(space);
      this->init();
    }

    template<typename Scalar>
    Adapt<Scalar>::~Adapt()
    {
    }

    template<typename Scalar>
    void Adapt<Scalar>::init()
    {
      if(!this->errorCalculator)
        throw Exceptions::Exception("Error calculator must not be NULL in Adapt::Adapt().");

      this->num = spaces.size();

      if(this->num > H2D_MAX_COMPONENTS)
        throw Exceptions::ValueException("components", this->num, 0, H2D_MAX_COMPONENTS);

      for(unsigned int i = 0; i < this->num; i++)
      {
        if(!spaces[i])
          throw Exceptions::NullException(0, i);
        spaces[i]->check();
      }

      regularization = -1;

      this->set_strategy(AdaptStoppingCriterionCumulative);
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_strategy(AdaptivityStoppingCriterion strategy_, double threshold_)
    {
      this->strategy = strategy_;
      this->threshold = threshold_;
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_regularization_level(int regularization_)
    {
      this->regularization = regularization_;
    }

    template<typename Scalar>
    bool Adapt<Scalar>::add_refinement(double processed_error_squared, double max_error_squared, double error_sum_squared, typename ErrorCalculator<Scalar>::ElementReference* element_reference)
    {
      switch(this->strategy)
      {
      case AdaptStoppingCriterionCumulative:
        if(processed_error_squared > (threshold*threshold) * error_sum_squared)
          return false;
        else
          return true;
        break;
      case AdaptStoppingCriterionSingleElement:
        if(*(element_reference->error) > (threshold*threshold) * max_error_squared)
          return true;
        else
          return false;
        break;
      default:
        throw Exceptions::Exception("Unknown strategy in Adaptivity.");
      }

      return false;
    }

    template<typename Scalar>
    void Adapt<Scalar>::init_adapt(Hermes::vector<RefinementSelectors::Selector<Scalar>*>& refinement_selectors, int** element_refinement_location, MeshSharedPtr* meshes)
    {
      // Start time measurement.
      this->tick();

      // Important, sets the current caughtException to NULL.
      this->caughtException = NULL;

      // Checks.
      if(!this->errorCalculator->elements_stored)
        throw Exceptions::Exception("element errors have to be calculated first, call ErrorCalculator::calculate_errors().");
      if(refinement_selectors.empty())
        throw Exceptions::NullException(1);
      if(spaces.size() != refinement_selectors.size())
        throw Exceptions::LengthException(1, refinement_selectors.size(), spaces.size());

      // Get meshes
      for (int j = 0; j < this->num; j++)
      {
        meshes[j] = this->spaces[j]->get_mesh();
        element_refinement_location[j] = (int*)calloc(meshes[j]->get_max_element_id(), sizeof(int));
      }

      for (int j = 0; j < this->num; j++)
      {
        Solution<Scalar>* temp = dynamic_cast<Solution<Scalar>*>(this->errorCalculator->fine_solutions[j].get());
        if(temp)
          temp->enable_transform(false);
      }
    }

    template<typename Scalar>
    int Adapt<Scalar>::calculate_attempted_element_refinements_count()
    {
      unsigned int attempted_element_refinements_count = 0;

      // Processed error so far.
      double processed_error_squared = 0.0;
      // Maximum error - the first one in the error calculator's array.
      double max_error_squared = *(this->errorCalculator->element_references[0].error);
      // Total error.
      double sum_error_squared = this->errorCalculator->get_total_error_squared();

      // For ALL elements on ALL meshes.
      // The stopping condition for this loop is the stopping condition for adaptivity.
      for(int element_inspected_i = 0; element_inspected_i < this->errorCalculator->num_act_elems; element_inspected_i++)
      {
        // Get the element info from the error calculator.
        ErrorCalculator<Scalar>::ElementReference element_reference = this->errorCalculator->element_references[element_inspected_i];

        // Ask the strategy if we should add this refinement or break the loop.
        if(!this->add_refinement(processed_error_squared, max_error_squared, sum_error_squared, &element_reference))
          break;

        processed_error_squared += *(element_reference.error);

        attempted_element_refinements_count++;
      }

      return attempted_element_refinements_count;
    }

    template<typename Scalar>
    bool Adapt<Scalar>::adapt(Hermes::vector<RefinementSelectors::Selector<Scalar> *>& refinement_selectors)
    {
      // Initialize.
      MeshSharedPtr meshes[H2D_MAX_COMPONENTS];
      int* element_refinement_location[H2D_MAX_COMPONENTS];
      this->init_adapt(refinement_selectors, element_refinement_location, meshes);

      // This is the number of refinements attempted.
      int attempted_element_refinements_count = calculate_attempted_element_refinements_count();

      // Time measurement.
      this->tick();
      this->info("\tAdaptivity: data preparation duration: %f s.", this->last());

      // List of indices of elements that are going to be processed
      ElementToRefine* elements_to_refine = new ElementToRefine[attempted_element_refinements_count];

      // Projected solutions obtaining.
      MeshFunctionSharedPtr<Scalar>* rslns = new MeshFunctionSharedPtr<Scalar>[this->num];
      OGProjection<Scalar> ogProjection;

      for(unsigned int i = 0; i < this->num; i++)
      {
        rslns[i] = MeshFunctionSharedPtr<Scalar>(new Solution<Scalar>());

        typename Mesh::ReferenceMeshCreator ref_mesh_creator(this->spaces[i]->get_mesh());
        MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
        typename Space<Scalar>::ReferenceSpaceCreator ref_space_creator(this->spaces[i], ref_mesh);
        SpaceSharedPtr<Scalar> ref_space = ref_space_creator.create_ref_space();

        ogProjection.project_global(ref_space, this->errorCalculator->fine_solutions[i], rslns[i]);
      }

      // Parallel section
#pragma omp parallel num_threads(this->num_threads_used)
      {
        int thread_number = omp_get_thread_num();
        int start = (attempted_element_refinements_count / this->num_threads_used) * thread_number;
        int end = (attempted_element_refinements_count / this->num_threads_used) * (thread_number + 1);
        if(thread_number == this->num_threads_used - 1)
          end = attempted_element_refinements_count;

        // rslns cloning.
        MeshFunctionSharedPtr<Scalar>* current_rslns = new MeshFunctionSharedPtr<Scalar>[this->num];
        for(unsigned int i = 0; i < this->num; i++)
          current_rslns[i] = rslns[i]->clone();

        for(int id_to_refine = start; id_to_refine < end; id_to_refine++)
        {
          try
          {
            // Get the appropriate element reference from the error calculator.
            ErrorCalculator<Scalar>::ElementReference element_reference = this->errorCalculator->element_references[id_to_refine];
            int element_id = element_reference.element_id;
            int component = element_reference.comp;
            int current_order = this->spaces[component]->get_element_order(element_id);

            // Get refinement suggestion.
            ElementToRefine elem_ref(element_id, component);

            // Rsln[comp] may be unset if refinement_selectors[comp] == HOnlySelector or POnlySelector
            refinement_selectors[component]->select_refinement(meshes[component]->get_element(element_id), current_order, current_rslns[component].get(), elem_ref, this->errorCalculator->errorType);

            // Put this refinement to the storage.
            elements_to_refine[id_to_refine] = elem_ref;

          }
          catch(Hermes::Exceptions::Exception& exception)
          {
            if(this->caughtException == NULL)   
              this->caughtException = exception.clone();
          }
          catch(std::exception& exception)
          {
            if(this->caughtException == NULL)
              this->caughtException = new std::exception(exception);
          }
        }
      
        delete [] current_rslns;
      }

      delete [] rslns;

      if(this->caughtException)
      {
        this->deinit_adapt(element_refinement_location);
        throw *(this->caughtException);
        return false;
      }

      // Time measurement.
      this->tick();
      this->info("\tAdaptivity: refinement selection duration: %f s.", this->last());

      // Before applying, fix the shared mesh refinements.
      fix_shared_mesh_refinements(meshes, elements_to_refine, attempted_element_refinements_count, element_refinement_location, &refinement_selectors.front());

      // Apply refinements
      apply_refinements(elements_to_refine, attempted_element_refinements_count);

      // in singlemesh case, impose same orders across meshes
      homogenize_shared_mesh_orders(meshes);

      this->adapt_postprocess(meshes, attempted_element_refinements_count);

      delete [] elements_to_refine;

      this->deinit_adapt(element_refinement_location);

      return false;
    }

    template<typename Scalar>
    void Adapt<Scalar>::deinit_adapt(int** element_refinement_location)
    {
      // Get meshes
      for (int j = 0; j < this->num; j++)
        ::free(element_refinement_location[j]);

      for (int j = 0; j < this->num; j++)
      {
        Solution<Scalar>* temp = dynamic_cast<Solution<Scalar>*>(this->errorCalculator->fine_solutions[j].get());
        if(temp)
          temp->enable_transform(true);
      }
    }

    template<typename Scalar>
    void Adapt<Scalar>::adapt_postprocess(MeshSharedPtr* meshes, int element_refinements_count)
    {
      // mesh regularization
      if(this->regularization >= 0)
      {
        if(this->regularization == 0)
        {
          this->regularization = 1;
          this->warn("Total mesh regularization is not supported in adaptivity. 1-irregular mesh is used instead.");
        }
        for (int i = 0; i < this->num; i++)
        {
          int* parents;
          parents = meshes[i]->regularize(this->regularization);
          this->spaces[i]->distribute_orders(meshes[i], parents);
          ::free(parents);
        }
      }

      // since space changed, assign dofs:
      for(unsigned int i = 0; i < this->spaces.size(); i++)
        this->spaces[i]->assign_dofs();

      Element* e;
      for (int i = 0; i < this->num; i++)
      {
        for_all_active_elements(e, this->spaces[i]->get_mesh())
          this->spaces[i]->edata[e->id].changed_in_last_adaptation = false;
      }

      // EXTREMELY important - set the changed_in_last_adaptation to changed elements.
      for(int i = 0; i < element_refinements_count; i++)
      {
        ErrorCalculator<Scalar>::ElementReference element_reference = this->errorCalculator->element_references[i];
        int element_id = element_reference.element_id;
        int component = element_reference.comp;
        this->spaces[component]->edata[element_id].changed_in_last_adaptation = true;
      }
    }

    template<typename Scalar>
    bool Adapt<Scalar>::adapt(RefinementSelectors::Selector<Scalar>* refinement_selector)
    {
      if(!refinement_selector)
        throw Exceptions::NullException(1);
      Hermes::vector<RefinementSelectors::Selector<Scalar> *> refinement_selectors;
      refinement_selectors.push_back(refinement_selector);
      return adapt(refinement_selectors);
    }

    template<typename Scalar>
    void Adapt<Scalar>::fix_shared_mesh_refinements(MeshSharedPtr* meshes, ElementToRefine*& elems_to_refine, int num_elem_to_proc, int** idx, RefinementSelectors::Selector<Scalar>** refinement_selectors)
    {
      // For additions.
      std::vector<ElementToRefine> new_elems_to_refine;

      for(int inx = 0; inx < num_elem_to_proc; inx++)
      {
        ElementToRefine& elem_ref = elems_to_refine[inx];
        int current_quad_order = this->spaces[elem_ref.comp]->get_element_order(elem_ref.id);
        Element* current_elem = meshes[elem_ref.comp]->get_element(elem_ref.id);

        //select a refinement used by all components that share a mesh which is about to be refined
        int selected_refinement = elem_ref.split;
        for (int j = 0; j < this->num; j++)
        {
          if(selected_refinement == H2D_REFINEMENT_H) break; // iso refinement is max what can be recieved
          if(j != elem_ref.comp && meshes[j] == meshes[elem_ref.comp])
          { // if a mesh is shared
            int ii = idx[elem_ref.id][j];
            if(ii >= 0) { // and the sample element is about to be refined by another compoment
              const ElementToRefine& elem_ref_ii = elems_to_refine[ii];
              if(elem_ref_ii.split != selected_refinement && elem_ref_ii.split != H2D_REFINEMENT_P) { //select more complicated refinement
                if((elem_ref_ii.split == H2D_REFINEMENT_ANISO_H || elem_ref_ii.split == H2D_REFINEMENT_ANISO_V) && selected_refinement == H2D_REFINEMENT_P)
                  selected_refinement = elem_ref_ii.split;
                else
                  selected_refinement = H2D_REFINEMENT_H;
              }
            }
          }
        }

        //fix other refinements according to the selected refinement
        if(selected_refinement != H2D_REFINEMENT_P)
        {
          //get suggested orders for the selected refinement
          const int* suggested_orders = NULL;
          if(selected_refinement == H2D_REFINEMENT_H)
            suggested_orders = elem_ref.q;

          //update orders
          for (int j = 0; j < this->num; j++)
          {
            if(j != elem_ref.comp && meshes[j] == meshes[elem_ref.comp]) { // if components share the mesh
              // change currently processed refinement
              if(elem_ref.split != selected_refinement)
              {
                elem_ref.split = selected_refinement;
                refinement_selectors[j]->generate_shared_mesh_orders(current_elem, current_quad_order, elem_ref.split, elem_ref.p, suggested_orders);
              }

              // change other refinements
              int ii = idx[elem_ref.id][j];
              if(ii >= 0)
              {
                ElementToRefine& elem_ref_ii = elems_to_refine[ii];
                if(elem_ref_ii.split != selected_refinement)
                {
                  elem_ref_ii.split = selected_refinement;
                  refinement_selectors[j]->generate_shared_mesh_orders(current_elem, current_quad_order, elem_ref_ii.split, elem_ref_ii.p, suggested_orders);
                }
              }
              else
              { // element (of the other comp.) not refined at all: assign refinement
                ElementToRefine elem_ref_new(elem_ref.id, j);
                elem_ref_new.split = selected_refinement;
                refinement_selectors[j]->generate_shared_mesh_orders(current_elem, current_quad_order, elem_ref_new.split, elem_ref_new.p, suggested_orders);
                new_elems_to_refine.push_back(elem_ref_new);
              }
            }
          }
        }
      }

      // Adding the additions.
      if(new_elems_to_refine.size() > 0)
        elems_to_refine = (ElementToRefine*)realloc(elems_to_refine, (num_elem_to_proc + new_elems_to_refine.size()) * sizeof(ElementToRefine));
      for(int inx = 0; inx < new_elems_to_refine.size(); inx++)
        elems_to_refine[num_elem_to_proc + inx] = new_elems_to_refine[inx];
    }

    template<typename Scalar>
    void Adapt<Scalar>::homogenize_shared_mesh_orders(MeshSharedPtr* meshes)
    {
      Element* e;
      for (int i = 0; i < this->num; i++)
      {
        for_all_active_elements(e, meshes[i])
        {
          int current_quad_order = this->spaces[i]->get_element_order(e->id);
          int current_order_h = H2D_GET_H_ORDER(current_quad_order), current_order_v = H2D_GET_V_ORDER(current_quad_order);

          for (int j = 0; j < this->num; j++)
            if((j != i) && (meshes[j] == meshes[i])) // components share the mesh
            {
              int quad_order = this->spaces[j]->get_element_order(e->id);
              current_order_h = std::max(current_order_h, H2D_GET_H_ORDER(quad_order));
              current_order_v = std::max(current_order_v, H2D_GET_V_ORDER(quad_order));
            }

            this->spaces[i]->set_element_order_internal(e->id, H2D_MAKE_QUAD_ORDER(current_order_h, current_order_v));
        }
      }
    }

    template<typename Scalar>
    void Adapt<Scalar>::apply_refinements(ElementToRefine* elems_to_refine, int num_elem_to_process)
    {
      for (int i = 0; i < num_elem_to_process; i++)
        apply_refinement(elems_to_refine[i]);
    }

    template<typename Scalar>
    void Adapt<Scalar>::apply_refinement(const ElementToRefine& elem_ref)
    {
      SpaceSharedPtr<Scalar>& space = this->spaces[elem_ref.comp];
      MeshSharedPtr& mesh = space->get_mesh();

      Element* e = mesh->get_element(elem_ref.id);

      if(elem_ref.split == H2D_REFINEMENT_P)
      {
        space->set_element_order_internal(elem_ref.id, elem_ref.p[0]);
        space->edata[elem_ref.id].changed_in_last_adaptation = true;
      }
      else if(elem_ref.split == H2D_REFINEMENT_H)
      {
        if(e->active)
          mesh->refine_element_id(elem_ref.id);
        for (int j = 0; j < 4; j++)
        {
          space->set_element_order_internal(e->sons[j]->id, elem_ref.p[j]);
          space->edata[e->sons[j]->id].changed_in_last_adaptation = true;
        }
      }
      else
      {
        if(e->active)
          mesh->refine_element_id(elem_ref.id, elem_ref.split);
        for (int j = 0; j < 2; j++)
        {
          space->set_element_order_internal(e->sons[ (elem_ref.split == 1) ? j : j + 2 ]->id, elem_ref.p[j]);
          space->edata[e->sons[ (elem_ref.split == 1) ? j : j + 2 ]->id].changed_in_last_adaptation = true;
        }
      }
    }

    template HERMES_API class Adapt<double>;
    template HERMES_API class Adapt<std::complex<double> >;
  }
}
