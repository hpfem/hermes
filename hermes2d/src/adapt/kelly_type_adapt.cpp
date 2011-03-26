#include "kelly_type_adapt.h"

// #ifdef KELLY_TYPE_ADAPT_H_IS_REWORKED

KellyTypeAdapt::KellyTypeAdapt(Hermes::vector< Space* > spaces_,
                               Hermes::vector< ProjNormType > norms_,
                               bool ignore_visited_segments_,
                               Hermes::vector<interface_estimator_scaling_fn_t> interface_scaling_fns_)
  : Adapt(spaces_, norms_)
{
  error_estimators_surf.reserve(num);
  error_estimators_vol.reserve(num);

  if (interface_scaling_fns_.size() == 0)
  {
    interface_scaling_fns_.reserve(num);
    for (int i = 0; i < num; i++)
      interface_scaling_fns_.push_back(scale_by_element_diameter);
  }
  use_aposteriori_interface_scaling = true;
  interface_scaling_fns = interface_scaling_fns_;
  interface_scaling_const = boundary_scaling_const = volumetric_scaling_const = 1.0;
  ignore_visited_segments = ignore_visited_segments_;

  element_markers_conversion = spaces_[0]->get_mesh()->element_markers_conversion;
  boundary_markers_conversion = spaces_[0]->get_mesh()->boundary_markers_conversion;
}

KellyTypeAdapt::KellyTypeAdapt(Space* space_, 
                               ProjNormType norm_, 
                               bool ignore_visited_segments_, 
                               interface_estimator_scaling_fn_t interface_scaling_fn_) : Adapt(space_, norm_)
{ 
  if (interface_scaling_fn_ == NULL)
    interface_scaling_fns.push_back(scale_by_element_diameter);
  else
    interface_scaling_fns.push_back(interface_scaling_fn_);
  
  use_aposteriori_interface_scaling = true;

  interface_scaling_const = boundary_scaling_const = volumetric_scaling_const = 1.0;
  ignore_visited_segments = ignore_visited_segments_;
  
  element_markers_conversion = space_->get_mesh()->element_markers_conversion;
  boundary_markers_conversion = space_->get_mesh()->boundary_markers_conversion;
}

bool KellyTypeAdapt::adapt(double thr, int strat, int regularize, double to_be_processed)
{
  Hermes::vector<RefinementSelectors::Selector *> refinement_selectors;
  RefinementSelectors::HOnlySelector selector;
  for (int i = 0; i < this->num; i++)
    refinement_selectors.push_back(&selector);

  return Adapt::adapt(refinement_selectors, thr, strat, regularize, to_be_processed);
}

void KellyTypeAdapt::add_error_estimator_vol(KellyTypeAdapt::ErrorEstimatorForm* form)
{
  error_if(form->i < 0 || form->i >= this->num,
           "Invalid component number (%d), max. supported components: %d", form->i, H2D_MAX_COMPONENTS);

  form->adapt = this;
  this->error_estimators_vol.push_back(form);
}

void KellyTypeAdapt::add_error_estimator_surf(KellyTypeAdapt::ErrorEstimatorForm* form)
{
  error_if (form->i < 0 || form->i >= this->num,
            "Invalid component number (%d), max. supported components: %d", form->i, H2D_MAX_COMPONENTS);

  form->adapt = this;
  this->error_estimators_surf.push_back(form);
}

double KellyTypeAdapt::calc_err_internal(Hermes::vector<Solution *> slns,
                                         Hermes::vector<double>* component_errors,
                                         unsigned int error_flags)
{    
  int n = slns.size();
  error_if (n != this->num,
            "Wrong number of solutions.");

  TimePeriod tmr;

  for (int i = 0; i < n; i++)
  {
    this->sln[i] = slns[i];
    sln[i]->set_quad_2d(&g_quad_2d_std);
  }

  have_coarse_solutions = true;
  
  WeakForm::Stage stage;

  num_act_elems = 0;
  for (int i = 0; i < num; i++)
  {
    stage.meshes.push_back(sln[i]->get_mesh());
    stage.fns.push_back(sln[i]);

    num_act_elems += stage.meshes[i]->get_num_active_elements();
    int max = stage.meshes[i]->get_max_element_id();

    if (errors[i] != NULL) delete [] errors[i];
    errors[i] = new double[max];
    memset(errors[i], 0.0, sizeof(double) * max);
  }
/*
  for (unsigned int i = 0; i < error_estimators_vol.size(); i++)
    trset.insert(error_estimators_vol[i].ext.begin(), error_estimators_vol[i].ext.end());
  for (unsigned int i = 0; i < error_estimators_surf.size(); i++)
    trset.insert(error_estimators_surf[i].ext.begin(), error_estimators_surf[i].ext.end());
*/

  double total_norm = 0.0;

  bool calc_norm = false;
  if ((error_flags & HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL ||
      (error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL) calc_norm = true;

  double *norms = NULL;
  if (calc_norm)
  {
    norms = new double[num];
    memset(norms, 0.0, num * sizeof(double));
  }

  double *errors_components = new double[num];
  memset(errors_components, 0.0, num * sizeof(double));
  this->errors_squared_sum = 0.0;
  double total_error = 0.0;

  bool bnd[4];          // FIXME: magic number - maximal possible number of element surfaces
  SurfPos surf_pos[4];
  Element **ee;
  Traverse trav;

  // Reset the e->visited status of each element of each mesh (most likely it will be set to true from
  // the latest assembling procedure).
  if (ignore_visited_segments)
  {
    for (int i = 0; i < num; i++)
    {
      Element* e;
      for_all_active_elements(e, stage.meshes[i])
        e->visited = false;
    }
  }
  
  //WARNING: AD HOC debugging parameter.
  bool multimesh = false;

  // Begin the multimesh traversal.
  trav.begin(num, &(stage.meshes.front()), &(stage.fns.front()));
  while ((ee = trav.get_next_state(bnd, surf_pos)) != NULL)
  {   
    // Go through all solution components.
    for (int i = 0; i < num; i++)
    {
      if (ee[i] == NULL)
        continue;
      
      // Set maximum integration order for use in integrals, see limit_order()
      update_limit_table(ee[i]->get_mode());

      RefMap *rm = sln[i]->get_refmap();

      double err = 0.0;

      // Go through all volumetric error estimators.
      for (unsigned int iest = 0; iest < error_estimators_vol.size(); iest++)
      {
        // Skip current error estimator if it is assigned to a different component or geometric area
        // different from that of the current active element.
        if (error_estimators_vol[iest]->i != i)
          continue;
        /*
        if (error_estimators_vol[iest].area != ee[i]->marker)
          continue;
          */
        else if (error_estimators_vol[iest]->area != HERMES_ANY)
          continue;

        err += eval_volumetric_estimator(error_estimators_vol[iest], rm);
      }

      // Go through all surface error estimators (includes both interface and boundary est's).
      for (unsigned int iest = 0; iest < error_estimators_surf.size(); iest++)
      {
        if (error_estimators_surf[iest]->i != i)
          continue;

        for (int isurf = 0; isurf < ee[i]->get_num_surf(); isurf++)
        {
            /*
          if (error_estimators_surf[iest].area > 0 &&
              error_estimators_surf[iest].area != surf_pos[isurf].marker) continue;
          */
          if (bnd[isurf])   // Boundary
          {
            if (error_estimators_surf[iest]->area == H2D_DG_INNER_EDGE) continue;
            
            /*
            if (boundary_markers_conversion.get_internal_marker(error_estimators_surf[iest].area) < 0 &&
                error_estimators_surf[iest].area != HERMES_ANY) continue;
            */    
            
            err += eval_boundary_estimator(error_estimators_surf[iest], rm, surf_pos);
          }
          else              // Interface
          {
            if (error_estimators_surf[iest]->area != H2D_DG_INNER_EDGE) continue;

            /* BEGIN COPY FROM DISCRETE_PROBLEM.CPP */
            
            // 5 is for bits per page in the array.
            LightArray<NeighborSearch*> neighbor_searches(5);
            unsigned int num_neighbors = 0;
            DiscreteProblem::NeighborNode* root;
            int ns_index;
            
            dp.min_dg_mesh_seq = 0;
            for(int j = 0; j < num; j++)
              if(stage.meshes[j]->get_seq() < dp.min_dg_mesh_seq || j == 0)
                dp.min_dg_mesh_seq = stage.meshes[j]->get_seq();
            
            ns_index = stage.meshes[i]->get_seq() - dp.min_dg_mesh_seq; // = 0 for single mesh
            
            // Determine the minimum mesh seq in this stage.
            if (multimesh) 
            {              
              // Initialize the NeighborSearches.
              dp.init_neighbors(neighbor_searches, stage, isurf);
              
              // Create a multimesh tree;
              root = new DiscreteProblem::NeighborNode(NULL, 0);
              dp.build_multimesh_tree(root, neighbor_searches);
              
              // Update all NeighborSearches according to the multimesh tree.
              // After this, all NeighborSearches in neighbor_searches should have the same count 
              // of neighbors and proper set of transformations
              // for the central and the neighbor element(s) alike.
              // Also check that every NeighborSearch has the same number of neighbor elements.
              for(unsigned int j = 0; j < neighbor_searches.get_size(); j++)
                if(neighbor_searches.present(j)) {
                  NeighborSearch* ns = neighbor_searches.get(j);
                  dp.update_neighbor_search(ns, root);
                  if(num_neighbors == 0)
                    num_neighbors = ns->n_neighbors;
                  if(ns->n_neighbors != num_neighbors)
                    error("Num_neighbors of different NeighborSearches not matching in KellyTypeAdapt::calc_err_internal.");
                }
            }
            else
            {
              NeighborSearch *ns = new NeighborSearch(ee[i], stage.meshes[i]);
              ns->original_central_el_transform = stage.fns[i]->get_transform();
              ns->set_active_edge(isurf);
              ns->clear_initial_sub_idx();
              num_neighbors = ns->n_neighbors;
              neighbor_searches.add(ns, ns_index);
            }

            // Go through all segments of the currently processed interface (segmentation is caused
            // by hanging nodes on the other side of the interface).
            for (unsigned int neighbor = 0; neighbor < num_neighbors; neighbor++)
            {              
              if (ignore_visited_segments) {
                bool processed = true;
                for(unsigned int j = 0; j < neighbor_searches.get_size(); j++)
                  if(neighbor_searches.present(j))
                    if(!neighbor_searches.get(j)->neighbors.at(neighbor)->visited) {
                      processed = false;
                      break;
                    }
                if (processed) continue;
              }
              
              // Set the active segment in all NeighborSearches
              for(unsigned int j = 0; j < neighbor_searches.get_size(); j++)
                if(neighbor_searches.present(j)) {
                  neighbor_searches.get(j)->active_segment = neighbor;
                  neighbor_searches.get(j)->neighb_el = neighbor_searches.get(j)->neighbors[neighbor];
                  neighbor_searches.get(j)->neighbor_edge = neighbor_searches.get(j)->neighbor_edges[neighbor];
                }
                
              // Push all the necessary transformations to all functions of this stage.
              // The important thing is that the transformations to the current subelement are already there.
              // Also store the current neighbor element and neighbor edge in neighb_el, neighbor_edge.
              if (multimesh) 
              {
                for(unsigned int fns_i = 0; fns_i < stage.fns.size(); fns_i++)
                  for(unsigned int trf_i = 0; trf_i < neighbor_searches.get(stage.meshes[fns_i]->get_seq() - dp.min_dg_mesh_seq)->central_n_trans[neighbor]; trf_i++)
                    stage.fns[fns_i]->push_transform(neighbor_searches.get(stage.meshes[fns_i]->get_seq() - dp.min_dg_mesh_seq)->central_transformations[neighbor][trf_i]);
              }
              else
              {            
                // Push the transformations only to the solution on the current mesh
                for(unsigned int trf_i = 0; trf_i < neighbor_searches.get(ns_index)->central_n_trans[neighbor]; trf_i++)
                  stage.fns[i]->push_transform(neighbor_searches.get(ns_index)->central_transformations[neighbor][trf_i]);
              }
              /* END COPY FROM DISCRETE_PROBLEM.CPP */
              rm->force_transform(this->sln[i]->get_transform(), this->sln[i]->get_ctm());
              
              // The estimate is multiplied by 0.5 in order to distribute the error equally onto
              // the two neighboring elements.
              double central_err = 0.5 * eval_interface_estimator(error_estimators_surf[iest],
                                                                  rm, surf_pos, neighbor_searches, 
                                                                  ns_index);
              double neighb_err = central_err;

              // Scale the error estimate by the scaling function dependent on the element diameter
              // (use the central element's diameter).
              if (use_aposteriori_interface_scaling && interface_scaling_fns[i])
                central_err *= interface_scaling_fns[i](ee[i]->get_diameter());

              // In the case this edge will be ignored when calculating the error for the element on
              // the other side, add the now computed error to that element as well.
              if (ignore_visited_segments)
              {
                Element *neighb = neighbor_searches.get(i)->neighb_el;

                // Scale the error estimate by the scaling function dependent on the element diameter
                // (use the diameter of the element on the other side).
                if (use_aposteriori_interface_scaling && interface_scaling_fns[i])
                  neighb_err *= interface_scaling_fns[i](neighb->get_diameter());

                errors_components[i] += central_err + neighb_err;
                total_error += central_err + neighb_err;
                errors[i][ee[i]->id] += central_err;
                errors[i][neighb->id] += neighb_err;
              }
              else
                err += central_err;
              
              /* BEGIN COPY FROM DISCRETE_PROBLEM.CPP */
              
              // Clear the transformations from the RefMaps and all functions.
              if (multimesh)
                for(unsigned int fns_i = 0; fns_i < stage.fns.size(); fns_i++)
                  stage.fns[fns_i]->set_transform(neighbor_searches.get(stage.meshes[fns_i]->get_seq() - dp.min_dg_mesh_seq)->original_central_el_transform);
              else
                stage.fns[i]->set_transform(neighbor_searches.get(ns_index)->original_central_el_transform);

              rm->set_transform(neighbor_searches.get(ns_index)->original_central_el_transform);

              
              /* END COPY FROM DISCRETE_PROBLEM.CPP */
            }
            
            /* BEGIN COPY FROM DISCRETE_PROBLEM.CPP */
            
            if (multimesh)
              // Delete the multimesh tree;
              delete root;
            
            // Delete the neighbor_searches array.
            for(unsigned int j = 0; j < neighbor_searches.get_size(); j++) 
              if(neighbor_searches.present(j))
                delete neighbor_searches.get(j);
              
            /* END COPY FROM DISCRETE_PROBLEM.CPP */
            
          }
        }
      }

      if (calc_norm)
      {
        double nrm = eval_solution_norm(error_form[i][i], rm, sln[i]);
        norms[i] += nrm;
        total_norm += nrm;
      }

      errors_components[i] += err;
      total_error += err;
      errors[i][ee[i]->id] += err;

      ee[i]->visited = true;
    }
  }
  trav.finish();

  // Store the calculation for each solution component separately.
  if(component_errors != NULL)
  {
    component_errors->clear();
    for (int i = 0; i < num; i++)
    {
      if((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_ABS)
        component_errors->push_back(sqrt(errors_components[i]));
      else if ((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL)
        component_errors->push_back(sqrt(errors_components[i]/norms[i]));
      else
      {
        error("Unknown total error type (0x%x).", error_flags & HERMES_TOTAL_ERROR_MASK);
        return -1.0;
      }
    }
  }

  tmr.tick();
  error_time = tmr.accumulated();

  // Make the error relative if needed.
  if ((error_flags & HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL)
  {
    for (int i = 0; i < num; i++)
    {
      Element* e;
      for_all_active_elements(e, stage.meshes[i])
        errors[i][e->id] /= norms[i];
    }
  }

  this->errors_squared_sum = total_error;

  // Element error mask is used here, because this variable is used in the adapt()
  // function, where the processed error (sum of errors of processed element errors)
  // is matched to this variable.
  if ((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL)
    errors_squared_sum /= total_norm;

  // Prepare an ordered list of elements according to an error.
  fill_regular_queue(&(stage.meshes.front()));
  have_errors = true;

  if (calc_norm)
    delete [] norms;
  delete [] errors_components;

  // Return error value.
  if ((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_ABS)
    return sqrt(total_error);
  else if ((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL)
    return sqrt(total_error / total_norm);
  else
  {
    error("Unknown total error type (0x%x).", error_flags & HERMES_TOTAL_ERROR_MASK);
    return -1.0;
  }
}

double KellyTypeAdapt::eval_solution_norm(Adapt::MatrixFormVolError* form, RefMap *rm, MeshFunction* sln)
{
  // determine the integration order
  int inc = (sln->get_num_components() == 2) ? 1 : 0;
  Func<Ord>* ou = init_fn_ord(sln->get_fn_order() + inc);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = form->ord(1, &fake_wt, NULL, ou, ou, fake_e, NULL);
  int order = rm->get_inv_ref_order();
  order += o.get_order();

  Solution *sol = static_cast<Solution *>(sln);
  if(sol && sol->get_type() == HERMES_EXACT) {
    limit_order_nowarn(order);
  }
  else {
    limit_order(order);
  }

  ou->free_ord(); delete ou;
  delete fake_e;

  // eval the form
  Quad2D* quad = sln->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // init geometry and jacobian*weights
  Geom<double>* e = init_geom_vol(rm, order);
  double* jac = rm->get_jacobian(order);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * jac[i];

  // function values
  Func<scalar>* u = init_fn(sln, order);
  scalar res = form->value(np, jwt, NULL, u, u, e, NULL);

  e->free(); delete e;
  delete [] jwt;
  u->free_fn(); delete u;

  return std::abs(res);
}

double KellyTypeAdapt::eval_volumetric_estimator(KellyTypeAdapt::ErrorEstimatorForm* err_est_form, RefMap *rm)
{
  // determine the integration order
  int inc = (this->sln[err_est_form->i]->get_num_components() == 2) ? 1 : 0;

  Func<Ord>** oi = new Func<Ord>* [num];
  for (int i = 0; i < num; i++)
    oi[i] = init_fn_ord(this->sln[i]->get_fn_order() + inc);

  // Order of additional external functions.
  ExtData<Ord>* fake_ext = dp.init_ext_fns_ord(err_est_form->ext);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = err_est_form->ord(1, &fake_wt, oi, oi[err_est_form->i], fake_e, fake_ext);
  int order = rm->get_inv_ref_order();
  order += o.get_order();

  limit_order(order);

  // Clean up.
  for (int i = 0; i < this->num; i++)
    if (oi[i] != NULL) { oi[i]->free_ord(); delete oi[i]; }
  delete [] oi;
  delete fake_e;
  delete fake_ext;

  // eval the form
  Quad2D* quad = this->sln[err_est_form->i]->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // init geometry and jacobian*weights
  Geom<double>* e = init_geom_vol(rm, order);
  double* jac = rm->get_jacobian(order);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * jac[i];

  // function values
  Func<scalar>** ui = new Func<scalar>* [num];
  
  for (int i = 0; i < num; i++)
    ui[i] = init_fn(this->sln[i], order);
  
  ExtData<scalar>* ext = dp.init_ext_fns(err_est_form->ext, rm, order);

  scalar res = volumetric_scaling_const *
                err_est_form->value(np, jwt, ui, ui[err_est_form->i], e, ext);

  for (int i = 0; i < this->num; i++)
    if (ui[i] != NULL) { ui[i]->free_fn(); delete ui[i]; }
  delete [] ui;
  if (ext != NULL) { ext->free(); delete ext; }
  e->free(); delete e;
  delete [] jwt;

  return std::abs(res);
}

double KellyTypeAdapt::eval_boundary_estimator(KellyTypeAdapt::ErrorEstimatorForm* err_est_form, RefMap *rm, SurfPos* surf_pos)
{
  // determine the integration order
  int inc = (this->sln[err_est_form->i]->get_num_components() == 2) ? 1 : 0;
  Func<Ord>** oi = new Func<Ord>* [num];
  for (int i = 0; i < num; i++)
    oi[i] = init_fn_ord(this->sln[i]->get_edge_fn_order(surf_pos->surf_num) + inc);

  // Order of additional external functions.
  ExtData<Ord>* fake_ext = dp.init_ext_fns_ord(err_est_form->ext, surf_pos->surf_num);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = err_est_form->ord(1, &fake_wt, oi, oi[err_est_form->i], fake_e, fake_ext);
  int order = rm->get_inv_ref_order();
  order += o.get_order();

  limit_order(order);

  // Clean up.
  for (int i = 0; i < this->num; i++)
    if (oi[i] != NULL) { oi[i]->free_ord(); delete oi[i]; }
  delete [] oi;
  delete fake_e;
  delete fake_ext;

  // eval the form
  Quad2D* quad = this->sln[err_est_form->i]->get_quad_2d();
  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // init geometry and jacobian*weights
  Geom<double>* e = init_geom_surf(rm, surf_pos, eo);
  double3* tan = rm->get_tangent(surf_pos->surf_num, eo);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * tan[i][2];

  // function values
  Func<scalar>** ui = new Func<scalar>* [num];
  for (int i = 0; i < num; i++)
    ui[i] = init_fn(this->sln[i], eo);
  ExtData<scalar>* ext = dp.init_ext_fns(err_est_form->ext, rm, eo);

  scalar res = boundary_scaling_const *
                err_est_form->value(np, jwt, ui, ui[err_est_form->i], e, ext);

  for (int i = 0; i < this->num; i++)
    if (ui[i] != NULL) { ui[i]->free_fn(); delete ui[i]; }
  delete [] ui;
  if (ext != NULL) { ext->free(); delete ext; }
  e->free(); delete e;
  delete [] jwt;

  return std::abs(0.5*res);   // Edges are parameterized from 0 to 1 while integration weights
                              // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                              // the weights.
}

double KellyTypeAdapt::eval_interface_estimator(KellyTypeAdapt::ErrorEstimatorForm* err_est_form,
                                                RefMap *rm, SurfPos* surf_pos,
                                                LightArray<NeighborSearch*>& neighbor_searches, int neighbor_index)
{
  NeighborSearch* nbs = neighbor_searches.get(neighbor_index);
  Hermes::vector<MeshFunction*> slns;
  for (int i = 0; i < num; i++)
    slns.push_back(this->sln[i]);
  
  // Determine integration order.
  ExtData<Ord>* fake_ui = dp.init_ext_fns_ord(slns, neighbor_searches);
  
  // Order of additional external functions.
  // ExtData<Ord>* fake_ext = dp.init_ext_fns_ord(err_est_form->ext, nbs);

  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  Geom<Ord>* fake_e = new InterfaceGeom<Ord>(init_geom_ord(), nbs->neighb_el->marker, nbs->neighb_el->id, nbs->neighb_el->get_diameter());
  double fake_wt = 1.0;
  Ord o = err_est_form->ord(1, &fake_wt, fake_ui->fn, fake_ui->fn[err_est_form->i], fake_e, NULL);

  int order = rm->get_inv_ref_order();
  order += o.get_order();

  limit_order(order);

  // Clean up.
  if (fake_ui != NULL)
  {
    for (int i = 0; i < num; i++)
      delete fake_ui->fn[i];
    fake_ui->free_ord();
    delete fake_ui;
  }
  
  delete fake_e;
  
  //delete fake_ext;
  
  Quad2D* quad = this->sln[err_est_form->i]->get_quad_2d();
  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  int np = quad->get_num_points(eo);
  double3* pt = quad->get_points(eo);
  
  // Init geometry and jacobian*weights (do not use the NeighborSearch caching mechanism).
  double3* tan = rm->get_tangent(surf_pos->surf_num, eo);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * tan[i][2];
  
  Geom<double>* e = new InterfaceGeom<double>(init_geom_surf(rm, surf_pos, eo), 
                                              nbs->neighb_el->marker, 
                                              nbs->neighb_el->id, 
                                              nbs->neighb_el->get_diameter());
    
  // function values
  ExtData<scalar>* ui = dp.init_ext_fns(slns, neighbor_searches, order);
  //ExtData<scalar>* ext = dp.init_ext_fns(err_est_form->ext, nbs);

  scalar res = interface_scaling_const *
                err_est_form->value(np, jwt, ui->fn, ui->fn[err_est_form->i], e, NULL);

  if (ui != NULL) { ui->free(); delete ui; }
  //if (ext != NULL) { ext->free(); delete ext; }
  e->free(); delete e;
  delete [] jwt;

  return std::abs(0.5*res);   // Edges are parameterized from 0 to 1 while integration weights
                              // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                              // the weights.
}

// #endif
