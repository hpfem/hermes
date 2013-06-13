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
#include "refinement_selectors/candidates.h"
#include "function/exact_solution.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    AdaptStoppingCriterionCumulative<Scalar>::AdaptStoppingCriterionCumulative(double threshold) : threshold(threshold)
    {
    }

    template<typename Scalar>
    bool AdaptStoppingCriterionCumulative<Scalar>::add_refinement(ErrorCalculator<Scalar>* error_calculator, double processed_error_squared, double max_error_squared, int element_inspected_i)
    {
      if(processed_error_squared > (threshold*threshold) * error_calculator->get_total_error_squared())
        return false;
      else
        return true;
    }

    template<typename Scalar>
    AdaptStoppingCriterionSingleElement<Scalar>::AdaptStoppingCriterionSingleElement(double threshold) : threshold(threshold)
    {
    }

    template<typename Scalar>
    bool AdaptStoppingCriterionSingleElement<Scalar>::add_refinement(ErrorCalculator<Scalar>* error_calculator, double processed_error_squared, double max_error_squared, int element_inspected_i)
    {
      const typename ErrorCalculator<Scalar>::ElementReference& element_reference = error_calculator->get_element_reference(element_inspected_i);
      if(*(element_reference.error) > (threshold*threshold) * max_error_squared)
        return true;
      else
        return false;
    }

    template<typename Scalar>
    AdaptStoppingCriterionLevels<Scalar>::AdaptStoppingCriterionLevels(double threshold) : threshold(threshold)
    {
    }

    template<typename Scalar>
    bool AdaptStoppingCriterionLevels<Scalar>::add_refinement(ErrorCalculator<Scalar>* error_calculator, double processed_error_squared, double max_error_squared, int element_inspected_i)
    {
      const typename ErrorCalculator<Scalar>::ElementReference& element_reference = error_calculator->get_element_reference(element_inspected_i);
      if(element_inspected_i == 0)
        return true;
      else
      {
        const typename ErrorCalculator<Scalar>::ElementReference previous_element_reference = error_calculator->get_element_reference(element_inspected_i - 1);
        if(*(element_reference.error) > (threshold*threshold) * *((previous_element_reference).error))
          return true;
        else
          return false;
      }
    }

    template<typename Scalar>
    Adapt<Scalar>::Adapt(ErrorCalculator<Scalar>* errorCalculator, AdaptivityStoppingCriterion<Scalar>* strategy) : errorCalculator(errorCalculator), strategy(strategy)
    {
      this->init();
      this->set_defaults();
    }

    template<typename Scalar>
    Adapt<Scalar>::Adapt(Hermes::vector<SpaceSharedPtr<Scalar> > spaces_, ErrorCalculator<Scalar>* errorCalculator, AdaptivityStoppingCriterion<Scalar>* strategy) : errorCalculator(errorCalculator), spaces(spaces_), strategy(strategy)
    {
      this->init();
      this->set_defaults();
    }

    template<typename Scalar>
    Adapt<Scalar>::Adapt(SpaceSharedPtr<Scalar> space, ErrorCalculator<Scalar>* errorCalculator, AdaptivityStoppingCriterion<Scalar>* strategy) : errorCalculator(errorCalculator), strategy(strategy)
    {
      spaces.push_back(space);
      this->init();
      this->set_defaults();
    }

    template<typename Scalar>
    Adapt<Scalar>::~Adapt()
    {
      if(elements_to_refine)
        delete [] elements_to_refine;
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> > spaces)
    {
      this->spaces = spaces;
      this->num = spaces.size();
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_space(SpaceSharedPtr<Scalar> space)
    {
      this->spaces.clear();
      this->spaces.push_back(space);
      this->num = 1;
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_defaults()
    {
      regularization = -1;
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

      elements_to_refine = NULL;
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_strategy(AdaptivityStoppingCriterion<Scalar>* strategy_)
    {
      this->strategy = strategy_;
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_regularization_level(int regularization_)
    {
      this->regularization = regularization_;
    }

    template<typename Scalar>
    bool Adapt<Scalar>::isOkay() const
    {
      if(!this->spaces.size())
      {
        this->info("\tAdaptivity: spaces are missing.");
        return false;
      }

      if(!this->strategy)
      {
        this->info("\tAdaptivity: strategy is missing.");
        return false;
      }
      return true;
    }

    template<typename Scalar>
    void Adapt<Scalar>::init_adapt(Hermes::vector<RefinementSelectors::Selector<Scalar>*>& refinement_selectors, ElementToRefine*** element_refinement_location, MeshSharedPtr* meshes)
    {
      // Start time measurement.
      this->tick();

      // Check.
      this->check();

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
        element_refinement_location[j] = (ElementToRefine**)calloc(meshes[j]->get_max_element_id() + 1, sizeof(ElementToRefine*));
      }

      // Clearing.
      if(elements_to_refine)
        delete [] elements_to_refine;

      // Also handle the refinementInfoMeshFunction.
      for(int i = 0; i < this->num; i++)
        if(this->refinementInfoMeshFunction[i])
          this->refinementInfoMeshFunction[i].reset();
      
      // Init the caught parallel exception message.
      this->exceptionMessageCaughtInParallelBlock.clear();
    }

    template<typename Scalar>
    int Adapt<Scalar>::calculate_attempted_element_refinements_count()
    {
      // Processed error so far.
      double processed_error_squared = 0.0;
      // Maximum error - the first one in the error calculator's array.
      double max_error_squared = *(this->errorCalculator->get_element_reference(0).error);

      unsigned int attempted_element_refinements_count = 0;

      // For ALL elements on ALL meshes.
      // The stopping condition for this loop is the stopping condition for adaptivity.
      for(int element_inspected_i = 0; element_inspected_i < this->errorCalculator->num_act_elems; element_inspected_i++)
      {
        // Get the element info from the error calculator.
        typename ErrorCalculator<Scalar>::ElementReference element_reference = this->errorCalculator->get_element_reference(element_inspected_i);

        // Ask the strategy if we should add this refinement or break the loop.
        if(!this->strategy->add_refinement(this->errorCalculator, processed_error_squared, max_error_squared, element_inspected_i))
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
      ElementToRefine** element_refinement_location[H2D_MAX_COMPONENTS];
      this->init_adapt(refinement_selectors, element_refinement_location, meshes);

      // This is the number of refinements attempted.
      int attempted_element_refinements_count = calculate_attempted_element_refinements_count();

      // Time measurement.
      this->tick();
      this->info("\tAdaptivity: data preparation duration: %f s.", this->last());

      // List of indices of elements that are going to be processed
      this->elements_to_refine_count = attempted_element_refinements_count;
      this->elements_to_refine = new ElementToRefine[elements_to_refine_count];

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
            typename ErrorCalculator<Scalar>::ElementReference element_reference = this->errorCalculator->get_element_reference(id_to_refine);
            int element_id = element_reference.element_id;
            int component = element_reference.comp;
            int current_order = this->spaces[component]->get_element_order(element_id);

            // Get refinement suggestion.
            ElementToRefine elem_ref(element_id, component);

            // Rsln[comp] may be unset if refinement_selectors[comp] == HOnlySelector or POnlySelector
            if(refinement_selectors[component]->select_refinement(meshes[component]->get_element(element_id), current_order, current_rslns[component].get(), elem_ref))
            {
              // Put this refinement to the storage.
              elements_to_refine[id_to_refine] = elem_ref;
              element_refinement_location[component][element_id] = &elements_to_refine[id_to_refine];
            }
            else
              elements_to_refine[id_to_refine] = ElementToRefine(-1, -1);
          }
          catch(std::exception& exception)
          {
            this->exceptionMessageCaughtInParallelBlock = exception.what();
          }
        }

        delete [] current_rslns;
      }

      delete [] rslns;

      if(!this->exceptionMessageCaughtInParallelBlock.empty())
      {
        this->deinit_adapt(element_refinement_location);
        throw Hermes::Exceptions::Exception(this->exceptionMessageCaughtInParallelBlock.c_str());
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

      this->deinit_adapt(element_refinement_location);

      this->adapt_postprocess(meshes, attempted_element_refinements_count);

      return false;
    }

    template<typename Scalar>
    void Adapt<Scalar>::deinit_adapt(ElementToRefine*** element_refinement_location)
    {
      // Free data.
      for (int j = 0; j < this->num; j++)
        ::free(element_refinement_location[j]);
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
          for(int j = 0; j < this->num; j++)
            if(this->spaces[i]->get_mesh()->get_seq() == this->spaces[j]->get_mesh()->get_seq())
              this->spaces[j]->distribute_orders(meshes[i], parents);
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
        typename ErrorCalculator<Scalar>::ElementReference element_reference = this->errorCalculator->get_element_reference(i);
        int element_id = element_reference.element_id;
        int component = element_reference.comp;
        this->spaces[component]->edata[element_id].changed_in_last_adaptation = true;
      }
    }

    template<typename Scalar>
    MeshFunctionSharedPtr<double> Adapt<Scalar>::get_refinementInfoMeshFunction(int component)
    {
      if(component >= this->num)
        throw Exceptions::ValueException("component", component, this->num);

      // The value is ready to be returned if it has been initialized and no other adaptivity run has been performed since.
      if(this->refinementInfoMeshFunction[component])
        return this->refinementInfoMeshFunction[component];
      else
      {
        int* info_array = new int[this->spaces[component]->get_mesh()->get_num_elements()];
        memset(info_array, 0, sizeof(int) * this->spaces[component]->get_mesh()->get_num_elements());
        for(int i = 0; i < this->elements_to_refine_count; i++)
        {
          if(this->elements_to_refine[i].comp == component)
          {
            if(this->elements_to_refine[i].split == 0)
              info_array[this->elements_to_refine[i].id] = 2;
            else
            {
              int id;
              for(int sons_i = 0; sons_i < H2D_MAX_ELEMENT_SONS; sons_i++)
                if(this->spaces[component]->get_mesh()->get_element(this->elements_to_refine[i].id)->sons[sons_i])
                  info_array[this->spaces[component]->get_mesh()->get_element(this->elements_to_refine[i].id)->sons[sons_i]->id] = this->elements_to_refine[i].split == 1 ? 1 : 2;
            }
          }
        }
        this->refinementInfoMeshFunction[component].reset(new ExactSolutionConstantArray<double, int>(spaces[component]->get_mesh(), info_array, true));
        return this->refinementInfoMeshFunction[component];
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
    void Adapt<Scalar>::fix_shared_mesh_refinements(MeshSharedPtr* meshes, ElementToRefine*& elems_to_refine, int& num_elem_to_proc, ElementToRefine*** idx, RefinementSelectors::Selector<Scalar>** refinement_selectors)
    {
      // Simple returns.
      if(this->num == 1)
        return;

      // For additions.
      std::vector<ElementToRefine> new_elems_to_refine;

      for(int inx = 0; inx < num_elem_to_proc; inx++)
      {
        ElementToRefine& elem_ref = elems_to_refine[inx];
        if(elem_ref.id == -1)
          continue;

        //select a refinement used by all components that share a mesh which is about to be refined
        int selected_refinement = elem_ref.split;
        for (int j = 0; j < this->num; j++)
        {
          if(selected_refinement == H2D_REFINEMENT_H)
            break; // iso refinement is max what can be recieved

          // if a mesh is shared
          if(j != elem_ref.comp && meshes[j]->get_seq() == meshes[elem_ref.comp]->get_seq())
          {
            // and the sample element is about to be refined by another compoment
            if(idx[j][elem_ref.id])
            { 
              ElementToRefine* elem_ref_ii = idx[j][elem_ref.id];
              //select more complicated refinement
              if((elem_ref_ii->split != selected_refinement) && (elem_ref_ii->split != H2D_REFINEMENT_P))
              { 
                if((elem_ref_ii->split == H2D_REFINEMENT_ANISO_H || elem_ref_ii->split == H2D_REFINEMENT_ANISO_V) && selected_refinement == H2D_REFINEMENT_P)
                  selected_refinement = elem_ref_ii->split;
                else
                  selected_refinement = H2D_REFINEMENT_H;
              }
            }
          }
        }

        //fix other refinements according to the selected refinement
        if(selected_refinement != H2D_REFINEMENT_P)
        {
          // change currently processed refinement
          if(elem_ref.split != selected_refinement)
          {
            elem_ref.split = selected_refinement;
            ElementToRefine::copy_orders(elem_ref.refinement_polynomial_order, elem_ref.best_refinement_polynomial_order_type[selected_refinement]);
          }

          //update orders
          for (int j = 0; j < this->num; j++)
          {
            // if components share the mesh
            if(j != elem_ref.comp && meshes[j]->get_seq() == meshes[elem_ref.comp]->get_seq())
            { 
              // change other refinements
              if(idx[j][elem_ref.id])
              {
                ElementToRefine* elem_ref_ii = idx[j][elem_ref.id];
                if(elem_ref_ii->split != selected_refinement)
                {
                  elem_ref_ii->split = selected_refinement;
                  if(elem_ref_ii->best_refinement_polynomial_order_type[selected_refinement])
                    ElementToRefine::copy_orders(elem_ref_ii->refinement_polynomial_order, elem_ref_ii->best_refinement_polynomial_order_type[selected_refinement]);
                  else
                  {
                    // This should occur only if the original refinement was a p-refinement.
#ifdef _DEBUG
                    this->warn("The best refinement poly degree is missing in fix_shared_mesh_refinements.");
#endif 
                    elem_ref_ii->refinement_polynomial_order[3] = elem_ref_ii->refinement_polynomial_order[2] = elem_ref_ii->refinement_polynomial_order[1] = elem_ref_ii->refinement_polynomial_order[0];
                  }
                }
              }
              else
              { // element (of the other comp.) not refined at all: assign refinement
                ElementToRefine elem_ref_new(elem_ref.id, j);
                elem_ref_new.split = selected_refinement;
                ElementToRefine::copy_orders(elem_ref_new.refinement_polynomial_order, elem_ref.refinement_polynomial_order);
                new_elems_to_refine.push_back(elem_ref_new);
              }
            }
          }
        }
      }

      // Adding the additions.
      if(new_elems_to_refine.size() > 0)
      {
        ElementToRefine* new_elems_to_refine_array = new ElementToRefine[num_elem_to_proc + new_elems_to_refine.size()];
        memcpy(new_elems_to_refine_array, elems_to_refine, num_elem_to_proc * sizeof(ElementToRefine));
        delete [] elems_to_refine;
        elems_to_refine = new_elems_to_refine_array;

        for(int inx = 0; inx < new_elems_to_refine.size(); inx++)
          elems_to_refine[num_elem_to_proc + inx] = new_elems_to_refine[inx];
        num_elem_to_proc += new_elems_to_refine.size();
      }
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
          {
            if((j != i) && (meshes[j]->get_seq() == meshes[i]->get_seq()))
            {
              int quad_order = this->spaces[j]->get_element_order(e->id);
              current_order_h = std::max(current_order_h, H2D_GET_H_ORDER(quad_order));
              current_order_v = std::max(current_order_v, H2D_GET_V_ORDER(quad_order));
            }
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
      if(elem_ref.id == -1)
        return;

      SpaceSharedPtr<Scalar>& space = this->spaces[elem_ref.comp];

      Element* e = space->get_mesh()->get_element(elem_ref.id);

      if(elem_ref.split == H2D_REFINEMENT_P)
      {
        space->set_element_order_internal(elem_ref.id, elem_ref.refinement_polynomial_order[0]);
        space->edata[elem_ref.id].changed_in_last_adaptation = true;
      }
      else if(elem_ref.split == H2D_REFINEMENT_H)
      {
        if(e->active)
          space->get_mesh()->refine_element_id(elem_ref.id);
        for (int j = 0; j < 4; j++)
        {
          space->set_element_order_internal(e->sons[j]->id, elem_ref.refinement_polynomial_order[j]);
          space->edata[e->sons[j]->id].changed_in_last_adaptation = true;
        }
      }
      else
      {
        if(e->active)
        {
          space->get_mesh()->refine_element_id(elem_ref.id, (elem_ref.split == H2D_REFINEMENT_ANISO_H ? 1 : 2));
        }
        for (int j = 0; j < 2; j++)
        {
          space->set_element_order_internal(e->sons[ (elem_ref.split == H2D_REFINEMENT_ANISO_H) ? j : j + 2 ]->id, elem_ref.refinement_polynomial_order[j]);
          space->edata[e->sons[ (elem_ref.split == H2D_REFINEMENT_ANISO_H) ? j : j + 2 ]->id].changed_in_last_adaptation = true;
        }
      }
    }

    template HERMES_API class AdaptStoppingCriterionCumulative<double>;
    template HERMES_API class AdaptStoppingCriterionCumulative<std::complex<double> >;
    template HERMES_API class AdaptStoppingCriterionSingleElement<double>;
    template HERMES_API class AdaptStoppingCriterionSingleElement<std::complex<double> >;
    template HERMES_API class AdaptStoppingCriterionLevels<double>;
    template HERMES_API class AdaptStoppingCriterionLevels<std::complex<double> >;

    template HERMES_API class Adapt<double>;
    template HERMES_API class Adapt<std::complex<double> >;
  }
}
