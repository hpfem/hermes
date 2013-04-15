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

#include "error_calculator.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    Adapt<Scalar>::Adapt(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces_, ErrorCalculator<Scalar>* error_calculator) : error_calculator(error_calculator)
    {
      for(unsigned int i = 0; i < spaces_.size(); i++)
        spaces.push_back(spaces_[i]);

      this->init();
    }

    template<typename Scalar>
    Adapt<Scalar>::Adapt(SpaceSharedPtr<Scalar>& space) : error_calculator(error_calculator)
    {
      spaces.push_back(space);
      this->init();
    }

    template<typename Scalar>
    Adapt<Scalar>::~Adapt()
    {
      for (int i = 0; i < this->num; i++)
        delete [] errors[i];

      for (int i = 0; i < this->num; i++)
        for (int j = 0; j < this->num; j++)
          if(error_form[i][j] && own_forms[i][j])
          {
            delete error_form[i][j];
            own_forms[i][j] = false;
          }

      for(int i = 0; i < H2D_MAX_COMPONENTS; i++)
        delete [] own_forms[i];
      delete [] own_forms;
    }

    template<typename Scalar>
    void Adapt<Scalar>::init()
    {
      if(!this->error_calculator)
        throw Exceptions::Exception("Error calculator must not be NULL in Adapt::Adapt().");

      this->num = spaces.size();
      
      if(this->num > H2D_MAX_COMPONENTS)
        throw Exceptions::ValueException("components", this->num, 0, H2D_MAX_COMPONENTS);

      for(unsigned int i = 0; i < this->num; i++)
      {
        if(spaces[i].empty())
          throw Exceptions::NullException(0, i);
        spaces[i]->check();
      }

      regularization = -1;
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_strategy(typename Adapt<Scalar>::StoppingCriterionStrategy strategy_, double threshold_)
    {
      this->strategy = strategy_;
      this->threshold = threshold_;
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_regularization_level(int regularization-)
    {
      this->regularization = regularization_;
    }

    template<typename Scalar>
    bool Adapt<Scalar>::add_refinement(double processed_error_squared, double max_error, double threshold, int refinement_count, double current_error)
    {
      switch(this->strategy)
      {
      case Cumulative:
        if(processed_error_squared > (threshold*threshold) * errors_squared_sum)
          return false;
        else
          return true;
        break;
      case SingleElement:
        if(current_error > threshold * max_error)
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
    bool Adapt<Scalar>::adapt(Hermes::vector<RefinementSelectors::Selector<Scalar> *> refinement_selectors, double threshold, int regularize)
    {
      this->tick();

      // Important, sets the current caughtException to NULL.
      this->caughtException = NULL;

      if(!this->error_calculator->elements_stored)
        throw Exceptions::Exception("element errors have to be calculated first, call ErrorCalculator::calculate_errors().");

      if(refinement_selectors.empty())
        throw Exceptions::NullException(1);

      if(spaces.size() != refinement_selectors.size())
        throw Exceptions::LengthException(1, refinement_selectors.size(), spaces.size());

      //get meshes
      int max_id = -1;
      MeshSharedPtr meshes[H2D_MAX_COMPONENTS];
      for (int j = 0; j < this->num; j++)
      {
        meshes[j] = this->spaces[j]->get_mesh();
        if(rsln[j])
        {
          rsln[j]->set_quad_2d(&g_quad_2d_std);
          rsln[j]->enable_transform(false);
        }
        if(meshes[j]->get_max_element_id() > max_id)
          max_id = meshes[j]->get_max_element_id();
      }

      //reset element refinement info
      int** idx = new int*[max_id];
      for(int i = 0; i < max_id; i++)
        idx[i] = new int[num];

      Element* e;
      for(int j = 0; j < max_id; j++)
        for(int l = 0; l < this->num; l++)
          idx[j][l] = -1; // element not refined

      double processed_error_squared = 0.0;

      Hermes::vector<ElementToRefine> elem_inx_to_proc; //list of indices of elements that are going to be processed
      elem_inx_to_proc.reserve(this->error_calculator->num_act_elems);

      double max_error;
      int num_ignored_elem = 0; //a number of ignored elements
      int num_not_changed = 0; //a number of element that were not changed
      int num_priority_elem = 0; //a number of elements that were processed using priority queue

      // Structures traversed in reality using strategies.
      Hermes::vector<int> ids;
      Hermes::vector<int> components;
      Hermes::vector<int> current_orders;
      bool first_regular_element = true; // true if first regular element was not processed yet
      bool error_level_reached = false;

      for(int element_inspected_i = 0; element_inspected_i < this->error_calculator->num_act_elems; element_inspected_i++)
      {
        // Set the maximum error
        if(element_inspected_i == 0)
          max_error = thr * err_squared;

        ErrorCalculator<Scalar>::ElementReference* element_reference = this->error_calculator->element_references[element_inspected_i];

        if(this->add_refinement(processed_error_squared, max_error, this->threshold, element_inspected_i, (*element_reference->error))
        {
          err0_squared = err_squared;
          processed_error_squared += err_squared;
          ids.push_back(id);
          components.push_back(comp);
          current_orders.push_back(this->spaces[comp]->get_element_order(id));
          spaces[comp]->edata[id].changed_in_last_adaptation = true;
        }
        else
          break;
      }

      if(ids.empty())
      {
        this->warn("None of the elements selected for refinement was refined. Adaptivity step successful, returning 'true'.");
        return true;
      }

      // RefinementSelectors cloning.
      int num_threads_used = Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads);
      RefinementSelectors::Selector<Scalar>*** global_refinement_selectors = new RefinementSelectors::Selector<Scalar>**[num_threads_used];

      for(unsigned int i = 0; i < num_threads_used; i++)
      {
        global_refinement_selectors[i] = new RefinementSelectors::Selector<Scalar>*[refinement_selectors.size()];
        for (unsigned int j = 0; j < refinement_selectors.size(); j++)
        {
          if(i == 0)
            global_refinement_selectors[i][j] = refinement_selectors[j];
          else
          {
            global_refinement_selectors[i][j] = refinement_selectors[j]->clone();
            RefinementSelectors::ProjBasedSelector<Scalar>* proj_based_selector_i_j = dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(global_refinement_selectors[i][j]);
            RefinementSelectors::ProjBasedSelector<Scalar>* proj_based_selector_j = dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(refinement_selectors[j]);
            RefinementSelectors::ProjBasedSelector<Scalar>* optimum_selector_i_j = dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(global_refinement_selectors[i][j]);
            RefinementSelectors::ProjBasedSelector<Scalar>* optimum_selector_j = dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(refinement_selectors[j]);
            if(proj_based_selector_i_j)
            {
              proj_based_selector_i_j->cached_shape_vals_valid = proj_based_selector_j->cached_shape_vals_valid;
              proj_based_selector_i_j->cached_shape_ortho_vals = proj_based_selector_j->cached_shape_ortho_vals;
              proj_based_selector_i_j->cached_shape_vals = proj_based_selector_j->cached_shape_vals;
            }
            if(optimum_selector_i_j)
              optimum_selector_i_j->num_shapes = optimum_selector_j->num_shapes;
          }
        }
      }

      // Solution cloning.
      Solution<Scalar>*** rslns = new Solution<Scalar>**[num_threads_used];

      for(unsigned int i = 0; i < num_threads_used; i++)
      {
        rslns[i] = new Solution<Scalar>*[this->num];
        for (int j = 0; j < this->num; j++)
        {
          if(rsln[j])
            rslns[i][j] = static_cast<Solution<Scalar>* >(rsln[j]->clone());
        }
      }

      this->tick();
      this->info("\tAdaptivity: data preparation duration: %f s.", this->last());

      // For statistics.
      int num_elements_for_refinenement = ids.size();
      int* numberOfCandidates = new int[num_elements_for_refinenement];

      // Parallel section
#pragma omp parallel num_threads(num_threads_used)
      {
        int thread_number = omp_get_thread_num();
        int start = (num_elements_for_refinenement / num_threads_used) * thread_number;
        int end = (num_elements_for_refinenement / num_threads_used) * (thread_number + 1);
        if(thread_number == num_threads_used - 1)
          end = num_elements_for_refinenement;
        for(int id_to_refine = start; id_to_refine < end; id_to_refine++)
        {
          try
          {
            RefinementSelectors::Selector<Scalar>** current_refinement_selectors = global_refinement_selectors[thread_number];
            Solution<Scalar>** current_rslns = rslns[thread_number];

            // Get refinement suggestion
            ElementToRefine elem_ref(ids[id_to_refine], components[id_to_refine]);

            // rsln[comp] may be unset if refinement_selectors[comp] == HOnlySelector or POnlySelector
            current_refinement_selectors[components[id_to_refine]]->select_refinement(meshes[components[id_to_refine]]->get_element(ids[id_to_refine]), current_orders[id_to_refine], current_rslns[components[id_to_refine]], elem_ref);
            
            //add to a list of elements that are going to be refined
#pragma omp critical (elem_ref_being_pushed_back)
            {
              idx[ids[id_to_refine]][components[id_to_refine]] = elem_inx_to_proc.size();
              elem_inx_to_proc.push_back(elem_ref);
            }

            if(this->get_verbose_output())
            {
						  if(dynamic_cast<Hermes::Hermes2D::RefinementSelectors::OptimumSelector<Scalar>*>(current_refinement_selectors[components[id_to_refine]]))
								  numberOfCandidates[id_to_refine] = dynamic_cast<Hermes::Hermes2D::RefinementSelectors::OptimumSelector<Scalar>*>(current_refinement_selectors[components[id_to_refine]])->get_candidates().size();
              else
								  numberOfCandidates[id_to_refine] = 0;
            }
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
      }

      if(this->caughtException == NULL)
      {
        if(this->get_verbose_output())
        {
          int averageNumberOfCandidates = 0;
          for(int i = 0; i < num_elements_for_refinenement; i++)
              averageNumberOfCandidates += numberOfCandidates[i];
          averageNumberOfCandidates = averageNumberOfCandidates / num_elements_for_refinenement;

          this->info("\tAdaptivity: total number of refined Elements: %i.", num_elements_for_refinenement);
          this->info("\tAdaptivity: average number of candidates per refined Element: %i.", averageNumberOfCandidates);
        }
      }

      delete [] numberOfCandidates;

      this->tick();
      this->info("\tAdaptivity: refinement selection duration: %f s.", this->last());

      if(this->caughtException == NULL)
        fix_shared_mesh_refinements(meshes, elem_inx_to_proc, idx, global_refinement_selectors);

      for(unsigned int i = 0; i < num_threads_used; i++)
      {
        if(i > 0)
          for (unsigned int j = 0; j < refinement_selectors.size(); j++)
            delete global_refinement_selectors[i][j];
        delete [] global_refinement_selectors[i];
      }
      delete [] global_refinement_selectors;

      for(unsigned int i = 0; i < num_threads_used; i++)
      {
        if(rslns[i])
        {
          for (unsigned int j = 0; j < this->num; j++)
            if(rsln[j])
              delete rslns[i][j];
          delete [] rslns[i];
        }
      }
      delete [] rslns;

      for(int i = 0; i < max_id; i++)
        delete [] idx[i];
      delete [] idx;

      if(this->caughtException)
        throw *(this->caughtException);
      
      //apply refinements
      apply_refinements(elem_inx_to_proc);

      // in singlemesh case, impose same orders across meshes
      homogenize_shared_mesh_orders(meshes);

      // mesh regularization
      if(regularize >= 0)
      {
        if(regularize == 0)
        {
          regularize = 1;
          this->warn("Total mesh regularization is not supported in adaptivity. 1-irregular mesh is used instead.");
        }
        for (int i = 0; i < this->num; i++)
        {
          int* parents;
          parents = meshes[i]->regularize(regularize);
          this->spaces[i]->distribute_orders(meshes[i], parents);
          ::free(parents);
        }
      }

      for (int j = 0; j < this->num; j++)
        if(rsln[j])
          rsln[j]->enable_transform(true);

      //store for the user to retrieve
      last_refinements.swap(elem_inx_to_proc);

      have_errors = false;
      if(strat == 2)
        have_errors = true; // space without changes

      // since space changed, assign dofs:
      for(unsigned int i = 0; i < this->spaces.size(); i++)
        this->spaces[i]->assign_dofs();

      for (int i = 0; i < this->num; i++)
      {
        for_all_active_elements(e, this->spaces[i]->get_mesh())
          this->spaces[i]->edata[e->id].changed_in_last_adaptation = false;
        for(int id_to_refine = 0; id_to_refine < ids.size(); id_to_refine++)
          this->spaces[i]->edata[ids[id_to_refine]].changed_in_last_adaptation = false;
      }

      return false;
    }

    template<typename Scalar>
    bool Adapt<Scalar>::adapt(RefinementSelectors::Selector<Scalar>* refinement_selector, double thr, int strat,
      int regularize, double to_be_processed)
    {
      if(refinement_selector==NULL)
        throw Exceptions::NullException(1);
      Hermes::vector<RefinementSelectors::Selector<Scalar> *> refinement_selectors;
      refinement_selectors.push_back(refinement_selector);
      return adapt(refinement_selectors, thr, strat, regularize, to_be_processed);
    }

    template<typename Scalar>
    void Adapt<Scalar>::fix_shared_mesh_refinements(MeshSharedPtr* meshes, std::vector<ElementToRefine>& elems_to_refine,
      int** idx, RefinementSelectors::Selector<Scalar> *** refinement_selectors)
    {
      int num_elem_to_proc = elems_to_refine.size();

      RefinementSelectors::Selector<Scalar>** current_refinement_selectors;

      for(int inx = 0; inx < num_elem_to_proc; inx++)
      {
        current_refinement_selectors = refinement_selectors[omp_get_thread_num()];
        ElementToRefine& elem_ref = elems_to_refine[inx];
        int current_quad_order = this->spaces[elem_ref.comp]->get_element_order(elem_ref.id);
        Element* current_elem = meshes[elem_ref.comp]->get_element(elem_ref.id);

        //select a refinement used by all components that share a mesh which is about to be refined
        int selected_refinement = elem_ref.split;
        for (int j = 0; j < this->num; j++)
        {
          if(selected_refinement == H2D_REFINEMENT_H) break; // iso refinement is max what can be recieved
          if(j != elem_ref.comp && meshes[j] == meshes[elem_ref.comp]) { // if a mesh is shared
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
                current_refinement_selectors[j]->generate_shared_mesh_orders(current_elem, current_quad_order, elem_ref.split, elem_ref.p, suggested_orders);
              }

              // change other refinements
              int ii = idx[elem_ref.id][j];
              if(ii >= 0)
              {
                ElementToRefine& elem_ref_ii = elems_to_refine[ii];
                if(elem_ref_ii.split != selected_refinement)
                {
                  elem_ref_ii.split = selected_refinement;
                  current_refinement_selectors[j]->generate_shared_mesh_orders(current_elem, current_quad_order, elem_ref_ii.split, elem_ref_ii.p, suggested_orders);
                }
              }
              else
              { // element (of the other comp.) not refined at all: assign refinement
                ElementToRefine elem_ref_new(elem_ref.id, j);
                elem_ref_new.split = selected_refinement;
                current_refinement_selectors[j]->generate_shared_mesh_orders(current_elem, current_quad_order, elem_ref_new.split, elem_ref_new.p, suggested_orders);
                elems_to_refine.push_back(elem_ref_new);
              }
            }
          }
        }
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
    void Adapt<Scalar>::apply_refinements(std::vector<ElementToRefine>& elems_to_refine)
    {
      for (std::vector<ElementToRefine>::const_iterator elem_ref = elems_to_refine.begin();
        elem_ref != elems_to_refine.end(); elem_ref++)
        apply_refinement(*elem_ref);
    }

    template<typename Scalar>
    void Adapt<Scalar>::apply_refinement(const ElementToRefine& elem_ref)
    {
      SpaceSharedPtr<Scalar> space = this->spaces[elem_ref.comp];
      MeshSharedPtr mesh = space->get_mesh();

      Element* e;
      e = mesh->get_element(elem_ref.id);

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
