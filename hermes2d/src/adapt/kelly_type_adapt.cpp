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

/*
#include "kelly_type_adapt.h"

namespace Hermes
{
  namespace Hermes2D
  {
    static const int H2D_DG_INNER_EDGE_INT = -1234567;
    static const std::string H2D_DG_INNER_EDGE = "-1234567";

    template<typename Scalar>
    BasicKellyAdapt<Scalar>::ErrorEstimatorFormKelly::ErrorEstimatorFormKelly(int i, double const_by_laplacian) : KellyTypeAdapt<Scalar>::ErrorEstimatorForm(i, H2D_DG_INNER_EDGE), const_by_laplacian(const_by_laplacian)
    {}

    template<typename Scalar>
    void KellyTypeAdapt<Scalar>::ErrorEstimatorForm::setAsInterface() 
    { 
      this->set_area(H2D_DG_INNER_EDGE);
    }

    template<typename Scalar>
    KellyTypeAdapt<Scalar>::KellyTypeAdapt(Hermes::vector< SpaceSharedPtr<Scalar> >& spaces_,
                                           bool ignore_visited_segments_,
                                           Hermes::vector<const InterfaceEstimatorScalingFunction*> interface_scaling_fns_,
                                           Hermes::vector<NormType > norms_)
      : Adapt<Scalar>(spaces_, norms_)
    {
      error_estimators_surf.reserve(this->num);
      error_estimators_vol.reserve(this->num);

      if(interface_scaling_fns_.size() == 0)
      {
        interface_scaling_fns_.reserve(this->num);
        for (int i = 0; i < this->num; i++)
          interface_scaling_fns_.push_back(new ScaleByElementDiameter);
      }
      use_aposteriori_interface_scaling = true;
      interface_scaling_fns = interface_scaling_fns_;
      interface_scaling_const = boundary_scaling_const = volumetric_scaling_const = 1.0;
      ignore_visited_segments = ignore_visited_segments_;

      element_markers_conversion = spaces_[0]->get_mesh()->element_markers_conversion;
      boundary_markers_conversion = spaces_[0]->get_mesh()->boundary_markers_conversion;
    }

    template<typename Scalar>
    KellyTypeAdapt<Scalar>::KellyTypeAdapt(SpaceSharedPtr<Scalar> space_,
                                           bool ignore_visited_segments_,
                                           const InterfaceEstimatorScalingFunction* interface_scaling_fn_,
                                           NormType norm_)
      : Adapt<Scalar>(space_, norm_)
    {
      if(interface_scaling_fn_ == NULL)
        interface_scaling_fns.push_back(new ScaleByElementDiameter);
      else
        interface_scaling_fns.push_back(interface_scaling_fn_);

      use_aposteriori_interface_scaling = true;

      interface_scaling_const = boundary_scaling_const = volumetric_scaling_const = 1.0;
      ignore_visited_segments = ignore_visited_segments_;

      element_markers_conversion = space_->get_mesh()->element_markers_conversion;
      boundary_markers_conversion = space_->get_mesh()->boundary_markers_conversion;
    }

    template<typename Scalar>
    bool KellyTypeAdapt<Scalar>::adapt(double thr, int strat, int regularize, double to_be_processed)
    {
      Hermes::vector<RefinementSelectors::Selector<Scalar> *> refinement_selectors;
      RefinementSelectors::HOnlySelector<Scalar> selector;
      for (int i = 0; i < this->num; i++)
        refinement_selectors.push_back(&selector);

      return Adapt<Scalar>::adapt(refinement_selectors, thr, strat, regularize, to_be_processed);
    }

    template<typename Scalar>
    void KellyTypeAdapt<Scalar>::add_error_estimator_vol(typename KellyTypeAdapt<Scalar>::ErrorEstimatorForm* form)
    {
      if(form->i < 0 || form->i >= this->num)
        throw Exceptions::ValueException("component number", form->i, 0, this->num);

      form->adapt = this;
      this->error_estimators_vol.push_back(form);
    }

    template<typename Scalar>
    void KellyTypeAdapt<Scalar>::add_error_estimator_surf(typename KellyTypeAdapt<Scalar>::ErrorEstimatorForm* form)
    {
      if(form->i < 0 || form->i >= this->num)
        throw Exceptions::ValueException("component number", form->i, 0, this->num);

      form->adapt = this;
      this->error_estimators_surf.push_back(form);
    }

    template<typename Scalar>
    double KellyTypeAdapt<Scalar>::calc_err_internal(Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns,
                                                     Hermes::vector<double>* component_errors,
                                                     unsigned int error_flags)
    {
      if(slns.size() != this->num)
        throw Hermes::Exceptions::LengthException(0, slns.size(), this->num);

      this->tick();

      for (int i = 0; i < this->num; i++)
      {
        Solution<Scalar>* solution = dynamic_cast<Solution<Scalar>*>(slns[i].get());
        if(solution == NULL)
          throw Exceptions::Exception("Passed solution is in fact not a Solution instance in KellyTypeAdapt::calc_err_*().");

        this->sln[i] = solution;
        this->sln[i]->set_quad_2d(&g_quad_2d_std);
      }

      this->have_coarse_solutions = true;

      MeshSharedPtr* meshes = new MeshSharedPtr[this->num];
      Transformable** fns = new Transformable*[this->num];

      this->num_act_elems = 0;
      for (int i = 0; i < this->num; i++)
      {
        meshes[i] = (this->sln[i]->get_mesh());
        fns[i] = (this->sln[i]);

        this->num_act_elems += meshes[i]->get_num_active_elements();
        int max = meshes[i]->get_max_element_id();

        if(this->errors[i] != NULL) delete [] this->errors[i];
        this->errors[i] = new double[max];
        memset(this->errors[i], 0, sizeof(double) * max);
      }

      double total_norm = 0.0;

      bool calc_norm = false;
      if((error_flags & this->HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL ||
        (error_flags & this->HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL) calc_norm = true;

      double *norms = NULL;
      if(calc_norm)
      {
        norms = new double[this->num];
        memset(norms, 0, this->num * sizeof(double));
      }

      double *errors_components = new double[this->num];
      memset(errors_components, 0, this->num * sizeof(double));
      this->errors_squared_sum = 0.0;
      double total_error = 0.0;

      bool bnd[H2D_MAX_NUMBER_EDGES];
      SurfPos surf_pos[H2D_MAX_NUMBER_EDGES];
      Traverse::State *ee;
      Traverse trav(this->num);

      // Reset the e->visited status of each element of each mesh (most likely it will be set to true from
      // the latest assembling procedure).
      if(ignore_visited_segments)
      {
        for (int i = 0; i < this->num; i++)
        {
          Element* e;
          for_all_active_elements(e, meshes[i])
            e->visited = false;
        }
      }

      // Begin the multimesh traversal.
      trav.begin(this->num, meshes, fns);
      while ((ee = trav.get_next_state()) != NULL)
      {
        bnd[0] = ee->bnd[0];
        bnd[1] = ee->bnd[1];
        bnd[2] = ee->bnd[2];
        bnd[3] = ee->bnd[3];

        surf_pos[0].marker = ee->rep->en[0]->marker;
        surf_pos[1].marker = ee->rep->en[1]->marker;
        surf_pos[2].marker = ee->rep->en[2]->marker;
        surf_pos[3].marker = ee->rep->en[3]->marker;

        surf_pos[0].surf_num = 0;
        surf_pos[1].surf_num = 1;
        surf_pos[2].surf_num = 2;
        surf_pos[3].surf_num = 3;

        // Go through all solution components.
        for (int i = 0; i < this->num; i++)
        {
          if(ee->e[i] == NULL)
            continue;

          // Set maximum integration order for use in integrals, see limit_order()
          update_limit_table(ee->e[i]->get_mode());

          RefMap *rm = this->sln[i]->get_refmap();

          double err = 0.0;

          // Go through all volumetric error estimators.
          for (unsigned int iest = 0; iest < error_estimators_vol.size(); iest++)
          {
            // Skip current error estimator if it is assigned to a different component or geometric area
            // different from that of the current active element.

            if(error_estimators_vol[iest]->i != i)
              continue;

            if(error_estimators_vol[iest]->area != HERMES_ANY)
              if(!element_markers_conversion.get_internal_marker(error_estimators_vol[iest]->area).valid || element_markers_conversion.get_internal_marker(error_estimators_vol[iest]->area).marker != ee->e[i]->marker)
                continue;

            err += eval_volumetric_estimator(error_estimators_vol[iest], rm);
          }

          // Go through all surface error estimators (includes both interface and boundary est's).
          for (unsigned int iest = 0; iest < error_estimators_surf.size(); iest++)
          {
            if(error_estimators_surf[iest]->i != i)
              continue;

            for (int isurf = 0; isurf < ee->e[i]->get_nvert(); isurf++)
            {
              if(bnd[isurf])   // Boundary
              {
                if(error_estimators_surf[iest]->area != HERMES_ANY)
                {
                  if(!boundary_markers_conversion.get_internal_marker(error_estimators_surf[iest]->area).valid)
                    continue;
                  int imarker = boundary_markers_conversion.get_internal_marker(error_estimators_surf[iest]->area).marker;

                  if(imarker == H2D_DG_INNER_EDGE_INT)
                    continue;
                  if(imarker != surf_pos[isurf].marker)
                    continue;
                }

                err += eval_boundary_estimator(error_estimators_surf[iest], rm, &surf_pos[isurf]);
              }
              else              // Interface
              {
                if(error_estimators_surf[iest]->area != H2D_DG_INNER_EDGE)
                  continue;

                // BEGIN COPY FROM DISCRETE_PROBLEM.CPP

                // 5 is for bits per page in the array.
                LightArray<NeighborSearch<Scalar>*> neighbor_searches(5);
                unsigned int num_neighbors = 0;
                NeighborNode* root;
                int ns_index;

                // Determine the minimum mesh seq in this stage.
                unsigned int min_dg_mesh_seq = 0;
                for(unsigned int j = 0; j < this->spaces.size(); j++)
                  if(this->spaces[j]->get_mesh()->get_seq() < min_dg_mesh_seq || j == 0)
                    min_dg_mesh_seq = this->spaces[j]->get_mesh()->get_seq();

                ns_index = meshes[i]->get_seq() - min_dg_mesh_seq; // = 0 for single mesh

                // Initialize the NeighborSearches.
                this->dp.init_neighbors(neighbor_searches, ee, min_dg_mesh_seq);

                // Create a multimesh tree;
                root = new NeighborNode(NULL, 0);
                this->dp.build_multimesh_tree(root, neighbor_searches);

                // Update all NeighborSearches according to the multimesh tree.
                // After this, all NeighborSearches in neighbor_searches should have the same count
                // of neighbors and proper set of transformations
                // for the central and the neighbor element(s) alike.
                // Also check that every NeighborSearch has the same number of neighbor elements.
                for(unsigned int j = 0; j < neighbor_searches.get_size(); j++)
                {
                  if(neighbor_searches.present(j))
                  {
                    NeighborSearch<Scalar>* ns = neighbor_searches.get(j);
                    this->dp.update_neighbor_search(ns, root);
                    if(num_neighbors == 0)
                      num_neighbors = ns->n_neighbors;
                    if(ns->n_neighbors != num_neighbors)
                      throw Hermes::Exceptions::Exception("Num_neighbors of different NeighborSearches not matching in KellyTypeAdapt<Scalar>::calc_err_internal.");
                  }
                }

                // Go through all segments of the currently processed interface (segmentation is caused
                // by hanging nodes on the other side of the interface).
                for (unsigned int neighbor = 0; neighbor < num_neighbors; neighbor++)
                {
                  if(ignore_visited_segments)
                  {
                    bool processed = true;
                    for(unsigned int j = 0; j < neighbor_searches.get_size(); j++)
                      if(neighbor_searches.present(j))
                        if(!neighbor_searches.get(j)->neighbors.at(neighbor)->visited)
                        {
                          processed = false;
                          break;
                        }

                    if(processed) continue;
                  }

                  // We do not use cache_e and cache_jwt here.

                  // Set the active segment in all NeighborSearches
                  for(unsigned int j = 0; j < neighbor_searches.get_size(); j++)
                  {
                    if(neighbor_searches.present(j))
                    {
                      neighbor_searches.get(j)->active_segment = neighbor;
                      neighbor_searches.get(j)->neighb_el = neighbor_searches.get(j)->neighbors[neighbor];
                      neighbor_searches.get(j)->neighbor_edge = neighbor_searches.get(j)->neighbor_edges[neighbor];
                    }
                  }

                  // Push all the necessary transformations to all functions of this stage.
                  // The important thing is that the transformations to the current subelement are already there.
                  // Also store the current neighbor element and neighbor edge in neighb_el, neighbor_edge.
                  for(unsigned int fns_i = 0; fns_i < this->num; fns_i++)
                  {
                    NeighborSearch<Scalar> *ns = neighbor_searches.get(meshes[fns_i]->get_seq() - min_dg_mesh_seq);
                    if(ns->central_transformations[neighbor])
                      ns->central_transformations[neighbor]->apply_on(fns[fns_i]);
                  }

                  // END COPY FROM DISCRETE_PROBLEM.CPP
                  rm->force_transform(this->sln[i]->get_transform(), this->sln[i]->get_ctm());

                  // The estimate is multiplied by 0.5 in order to distribute the error equally onto
                  // the two neighboring elements.
                  double central_err = 0.5 * eval_interface_estimator(error_estimators_surf[iest],
                                                                      rm, &surf_pos[isurf], neighbor_searches,
                                                                      ns_index);
                  double neighb_err = central_err;

                  // Scale the error estimate by the scaling function dependent on the element diameter
                  // (use the central element's diameter).
                  if(use_aposteriori_interface_scaling && interface_scaling_fns[i])
                    if(!element_markers_conversion.get_user_marker(ee->e[i]->marker).valid)
                      throw Hermes::Exceptions::Exception("Marker not valid.");
                    else
                      central_err *= interface_scaling_fns[i]->value(ee->e[i]->get_diameter(), element_markers_conversion.get_user_marker(ee->e[i]->marker).marker);

                  // In the case this edge will be ignored when calculating the error for the element on
                  // the other side, add the now computed error to that element as well.
                  if(ignore_visited_segments)
                  {
                    Element *neighb = neighbor_searches.get(ns_index)->neighb_el;

                    // Scale the error estimate by the scaling function dependent on the element diameter
                    // (use the diameter of the element on the other side).
                    if(use_aposteriori_interface_scaling && interface_scaling_fns[i])
                      if(!element_markers_conversion.get_user_marker(neighb->marker).valid)
                      throw Hermes::Exceptions::Exception("Marker not valid.");
                    else
                      neighb_err *= interface_scaling_fns[i]->value(neighb->get_diameter(), element_markers_conversion.get_user_marker(neighb->marker).marker);

                    errors_components[i] += central_err + neighb_err;
                    total_error += central_err + neighb_err;
                    this->errors[i][ee->e[i]->id] += central_err;
                    this->errors[i][neighb->id] += neighb_err;
                  }
                  else
                    err += central_err;

                  // BEGIN COPY FROM DISCRETE_PROBLEM.CPP

                  // Clear the transformations from the RefMaps and all functions.
                  for(unsigned int fns_i = 0; fns_i < this->num; fns_i++)
                    fns[fns_i]->set_transform(neighbor_searches.get(meshes[fns_i]->get_seq() - min_dg_mesh_seq)->original_central_el_transform);

                  rm->set_transform(neighbor_searches.get(ns_index)->original_central_el_transform);

                  // END COPY FROM DISCRETE_PROBLEM.CPP
                }

                // BEGIN COPY FROM DISCRETE_PROBLEM.CPP

                // Delete the multimesh tree;
                delete root;

                // Delete the neighbor_searches array.
                for(unsigned int j = 0; j < neighbor_searches.get_size(); j++)
                  if(neighbor_searches.present(j))
                    delete neighbor_searches.get(j);

                // END COPY FROM DISCRETE_PROBLEM.CPP
              }
            }
          }

          if(calc_norm)
          {
            double nrm = eval_solution_norm(this->norm_form[i][i], rm, this->sln[i]);
            norms[i] += nrm;
            total_norm += nrm;
          }

          errors_components[i] += err;
          total_error += err;
          this->errors[i][ee->e[i]->id] += err;

          ee->e[i]->visited = true;
        }
      }
      trav.finish();

      // Store the calculation for each solution component separately.
      if(component_errors != NULL)
      {
        component_errors->clear();
        for (int i = 0; i < this->num; i++)
        {
          if((error_flags & this->HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_ABS)
            component_errors->push_back(sqrt(errors_components[i]));
          else if((error_flags & this->HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL)
            component_errors->push_back(sqrt(errors_components[i]/norms[i]));
          else
          {
            throw Hermes::Exceptions::Exception("Unknown total error type (0x%x).", error_flags & this->HERMES_TOTAL_ERROR_MASK);
            return -1.0;
          }
        }
      }

      this->tick();
      this->error_time = this->accumulated();

      // Make the error relative if needed.
      if((error_flags & this->HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL)
      {
        for (int i = 0; i < this->num; i++)
        {
          Element* e;
          for_all_active_elements(e, meshes[i])
            this->errors[i][e->id] /= norms[i];
        }
      }

      this->errors_squared_sum = total_error;

      // Element error mask is used here, because this variable is used in the adapt()
      // function, where the processed error (sum of errors of processed element errors)
      // is matched to this variable.
      if((error_flags & this->HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL)
        this->errors_squared_sum /= total_norm;

      // Prepare an ordered list of elements according to an error.
      this->fill_regular_queue(meshes);
      this->have_errors = true;

      if(calc_norm)
        delete [] norms;
      delete [] errors_components;

      // Return error value.
      if((error_flags & this->HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_ABS)
        return sqrt(total_error);
      else if((error_flags & this->HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL)
        return sqrt(total_error / total_norm);
      else
      {
        throw Hermes::Exceptions::Exception("Unknown total error type (0x%x).", error_flags & this->HERMES_TOTAL_ERROR_MASK);
        return -1.0;
      }
    }

    template<typename Scalar>
    double KellyTypeAdapt<Scalar>::eval_solution_norm(typename Adapt<Scalar>::MatrixFormVolError* form,
                                                      RefMap *rm, MeshFunctionSharedPtr<Scalar> sln)
    {
      // Determine the integration order.
      int inc = (sln->get_num_components() == 2) ? 1 : 0;
      Func<Hermes::Ord>* ou = init_fn_ord(sln->get_fn_order() + inc);

      double fake_wt = 1.0;
      Geom<Hermes::Ord>* fake_e = init_geom_ord();
      Hermes::Ord o = form->ord(1, &fake_wt, NULL, ou, ou, fake_e, NULL);
      int order = rm->get_inv_ref_order();
      order += o.get_order();

      limit_order(order, rm->get_active_element()->get_mode());

      ou->free_ord(); delete ou;
      delete fake_e;

      // Evaluate the form.
      Quad2D* quad = sln->get_quad_2d();
      double3* pt = quad->get_points(order, sln->get_active_element()->get_mode());
      int np = quad->get_num_points(order, sln->get_active_element()->get_mode());

      // Initialize geometry and jacobian*weights.
      Geom<double>* e = init_geom_vol(rm, order);
      double* jac = rm->get_jacobian(order);
      double* jwt = new double[np];
      for(int i = 0; i < np; i++)
        jwt[i] = pt[i][2] * jac[i];

      // Function values.
      Func<Scalar>* u = init_fn(sln.get(), order);
      Scalar res = form->value(np, jwt, NULL, u, u, e, NULL);

      e->free(); delete e;
      delete [] jwt;
      u->free_fn(); delete u;

      return std::abs(res);
    }

    template<typename Scalar>
    double KellyTypeAdapt<Scalar>::eval_volumetric_estimator(typename KellyTypeAdapt<Scalar>::ErrorEstimatorForm* err_est_form,
                                                             RefMap *rm)
    {
      // Determine the integration order.
      int inc = (this->sln[err_est_form->i]->get_num_components() == 2) ? 1 : 0;

      Func<Hermes::Ord>** oi = new Func<Hermes::Ord>*[this->num];
      for (int i = 0; i < this->num; i++)
        oi[i] = init_fn_ord(this->sln[i]->get_fn_order() + inc);

      // Polynomial order of additional external functions.
      Func<Hermes::Ord>** fake_ext_fn = new Func<Hermes::Ord>*[err_est_form->ext.size()];
      for (int i = 0; i < err_est_form->ext.size(); i++)
        fake_ext_fn[i] = init_fn_ord(err_est_form->ext[i]->get_fn_order());

      double fake_wt = 1.0;
      Geom<Hermes::Ord>* fake_e = init_geom_ord();

      DiscontinuousFunc<Hermes::Ord> oi_i(oi[err_est_form->i], false, false);

      Hermes::Ord o = err_est_form->ord(1, &fake_wt, oi, &oi_i, fake_e, fake_ext_fn);
      int order = rm->get_inv_ref_order();
      order += o.get_order();

      limit_order(order, rm->get_active_element()->get_mode());

      // Clean up.
      for (int i = 0; i < this->num; i++)
      {
        if(oi[i] != NULL)
        {
          oi[i]->free_ord();
          delete oi[i];
        }
      }
      delete [] oi;
      delete fake_e;
      for(int i = 0; i < err_est_form->ext.size(); i++)
        fake_ext_fn[i]->free_ord();
      delete [] fake_ext_fn;

      // eval the form
      Quad2D* quad = this->sln[err_est_form->i]->get_quad_2d();
      double3* pt = quad->get_points(order, rm->get_active_element()->get_mode());
      int np = quad->get_num_points(order, rm->get_active_element()->get_mode());

      // Initialize geometry and jacobian*weights
      Geom<double>* e = init_geom_vol(rm, order);
      double* jac = rm->get_jacobian(order);
      double* jwt = new double[np];
      for(int i = 0; i < np; i++)
        jwt[i] = pt[i][2] * jac[i];

      // Function values.
      Func<Scalar>** ui = new Func<Scalar>*[this->num];

      for (int i = 0; i < this->num; i++)
        ui[i] = init_fn(this->sln[i], order);

      Func<Scalar>** ext_fn = new Func<Scalar>*[err_est_form->ext.size()];
      for (unsigned i = 0; i < err_est_form->ext.size(); i++)
      {
        if(err_est_form->ext[i] != NULL)
          ext_fn[i] = init_fn(err_est_form->ext[i].get(), order);
        else
          ext_fn[i] = NULL;
      }

      DiscontinuousFunc<Scalar> ui_i(ui[err_est_form->i], false, false);

      Scalar res = volumetric_scaling_const * err_est_form->value(np, jwt, ui, &ui_i, e, ext_fn);

      for (int i = 0; i < this->num; i++)
      {
        if(ui[i] != NULL)
        {
          ui[i]->free_fn();
          delete ui[i];
        }
      }
      delete [] ui;

      for(int i = 0; i < err_est_form->ext.size(); i++)
        fake_ext_fn[i]->free_fn();
      delete [] fake_ext_fn;

      e->free();
      delete e;

      delete [] jwt;

      return std::abs(res);
    }

    template<typename Scalar>
    double KellyTypeAdapt<Scalar>::eval_boundary_estimator(typename KellyTypeAdapt<Scalar>::ErrorEstimatorForm* err_est_form,
                                                           RefMap *rm, SurfPos* surf_pos)
    {
      // Determine the integration order.
      int inc = (this->sln[err_est_form->i]->get_num_components() == 2) ? 1 : 0;
      Func<Hermes::Ord>** oi = new Func<Hermes::Ord>*[this->num];
      for (int i = 0; i < this->num; i++)
        oi[i] = init_fn_ord(this->sln[i]->get_edge_fn_order(surf_pos->surf_num) + inc);

      // Polynomial order of additional external functions.
      Func<Hermes::Ord>** fake_ext_fn = new Func<Hermes::Ord>*[err_est_form->ext.size()];
      for (int i = 0; i < err_est_form->ext.size(); i++)
        fake_ext_fn[i] = init_fn_ord(err_est_form->ext[i]->get_fn_order());

      double fake_wt = 1.0;
      Geom<Hermes::Ord>* fake_e = init_geom_ord();
      DiscontinuousFunc<Hermes::Ord> oi_i(oi[err_est_form->i], false, false);
      Hermes::Ord o = err_est_form->ord(1, &fake_wt, oi, &oi_i, fake_e, fake_ext_fn);
      int order = rm->get_inv_ref_order();
      order += o.get_order();

      limit_order(order, rm->get_active_element()->get_mode());

      // Clean up.
      for (int i = 0; i < this->num; i++)
        if(oi[i] != NULL)
        {
          oi[i]->free_ord();
          delete oi[i];
        }

      delete [] oi;
      delete fake_e;
      for(int i = 0; i < err_est_form->ext.size(); i++)
        fake_ext_fn[i]->free_ord();
      delete [] fake_ext_fn;

      // Evaluate the form.
      Quad2D* quad = this->sln[err_est_form->i]->get_quad_2d();
      int eo = quad->get_edge_points(surf_pos->surf_num, order, rm->get_active_element()->get_mode());
      double3* pt = quad->get_points(eo, rm->get_active_element()->get_mode());
      int np = quad->get_num_points(eo, rm->get_active_element()->get_mode());

      // Initialize geometry and jacobian*weights.
      double3* tan;
      Geom<double>* e = init_geom_surf(rm, surf_pos->surf_num, surf_pos->marker, eo, tan);
      double* jwt = new double[np];
      for(int i = 0; i < np; i++)
        jwt[i] = pt[i][2] * tan[i][2];

      // Function values
      Func<Scalar>** ui = new Func<Scalar>*[this->num];
      for (int i = 0; i < this->num; i++)
        ui[i] = init_fn(this->sln[i], eo);

      Func<Scalar>** ext_fn = new Func<Scalar>*[err_est_form->ext.size()];
      for (unsigned i = 0; i < err_est_form->ext.size(); i++)
      {
        if(err_est_form->ext[i] != NULL)
          ext_fn[i] = init_fn(err_est_form->ext[i].get(), order);
        else
          ext_fn[i] = NULL;
      }

      DiscontinuousFunc<Scalar> ui_i(ui[err_est_form->i], false, false);

      Scalar res = boundary_scaling_const *
        err_est_form->value(np, jwt, ui, &ui_i, e, ext_fn);

      for (int i = 0; i < this->num; i++)
        if(ui[i] != NULL)
        {
          ui[i]->free_fn();
          delete ui[i];
        }

      delete [] ui;
      for(int i = 0; i < err_est_form->ext.size(); i++)
        fake_ext_fn[i]->free_fn();
      delete [] fake_ext_fn;

      e->free();
      delete e;

      delete [] jwt;

      return std::abs(0.5*res);   // Edges are parameterized from 0 to 1 while integration weights
                                  // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                                  // the weights.
    }

    template<typename Scalar>
    double KellyTypeAdapt<Scalar>::eval_interface_estimator(typename KellyTypeAdapt<Scalar>::ErrorEstimatorForm* err_est_form,
                                                            RefMap *rm, SurfPos* surf_pos,
                                                            LightArray<NeighborSearch<Scalar>*>& neighbor_searches,
                                                            int neighbor_index)
    {
      NeighborSearch<Scalar>* nbs = neighbor_searches.get(neighbor_index);
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > slns;
      for (int i = 0; i < this->num; i++)
        slns.push_back(this->sln[i]);

      // Determine integration order.
      Func<Hermes::Ord>** fake_ext_fns = new Func<Hermes::Ord>*[err_est_form->ext.size()];
      for (unsigned int j = 0; j < err_est_form->ext.size(); j++)
      {
        int inc = (err_est_form->ext[j]->get_num_components() == 2) ? 1 : 0;
        int central_order = err_est_form->ext[j]->get_edge_fn_order(neighbor_searches.get(err_est_form->ext[j]->get_mesh()->get_seq())->active_edge) + inc;
        int neighbor_order = err_est_form->ext[j]->get_edge_fn_order(neighbor_searches.get(err_est_form->ext[j]->get_mesh()->get_seq())->neighbor_edge.local_num_of_edge) + inc;
        fake_ext_fns[j] = new DiscontinuousFunc<Hermes::Ord>(init_fn_ord(central_order), init_fn_ord(neighbor_order));
      }

      // Polynomial order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
      Geom<Hermes::Ord>* fake_e = new InterfaceGeom<Hermes::Ord>(init_geom_ord(), nbs->neighb_el->marker, nbs->neighb_el->id, Hermes::Ord(nbs->neighb_el->get_diameter()));
      double fake_wt = 1.0;
      DiscontinuousFunc<Hermes::Ord> fake_ext_fns_i(fake_ext_fns[err_est_form->i], false, false);
      
      Hermes::Ord o = err_est_form->ord(1, &fake_wt, fake_ext_fns, &fake_ext_fns_i, fake_e, NULL);

      int order = rm->get_inv_ref_order();
      order += o.get_order();

      limit_order(order, rm->get_active_element()->get_mode());

      // Clean up.
      for (int i = 0; i < this->num; i++)
      {
        fake_ext_fns[i]->free_ord();
        delete fake_ext_fns[i];
      }
      delete [] fake_ext_fns;

      fake_e->free_ord();
      delete fake_e;

      //delete fake_ext;

      Quad2D* quad = this->sln[err_est_form->i]->get_quad_2d();
      int eo = quad->get_edge_points(surf_pos->surf_num, order, rm->get_active_element()->get_mode());
      int np = quad->get_num_points(eo, rm->get_active_element()->get_mode());
      double3* pt = quad->get_points(eo, rm->get_active_element()->get_mode());

      // Initialize geometry and jacobian*weights (do not use the NeighborSearch caching mechanism).
      double3* tan;
      Geom<double>* e = new InterfaceGeom<double>(init_geom_surf(rm, surf_pos->surf_num, surf_pos->marker, eo, tan),
                                                  nbs->neighb_el->marker,
                                                  nbs->neighb_el->id,
                                                  nbs->neighb_el->get_diameter());

      double* jwt = new double[np];
      for(int i = 0; i < np; i++)
        jwt[i] = pt[i][2] * tan[i][2];

      // Function values.
      DiscontinuousFunc<Scalar>** ui = this->dp.init_ext_fns(slns, neighbor_searches, order, 0);

      Scalar res = interface_scaling_const *
        err_est_form->value(np, jwt, NULL, ui[err_est_form->i], e, NULL);

      if(ui != NULL)
      {
        for(unsigned int i = 0; i < slns.size(); i++)
          ui[i]->free_fn();
        delete [] ui;
      }

      e->free();
      delete e;

      delete [] jwt;

      return std::abs(0.5*res);   // Edges are parameterized from 0 to 1 while integration weights
                                  // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                                  // the weights.
    }

    // #endif
    template HERMES_API class KellyTypeAdapt<double>;
    template HERMES_API class KellyTypeAdapt<std::complex<double> >;
    template HERMES_API class BasicKellyAdapt<double>;
    template HERMES_API class BasicKellyAdapt<std::complex<double> >;
  }
}
*/