#include "kelly_type_adapt.h"

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
}

bool KellyTypeAdapt::adapt(double thr, int strat, int regularize, double to_be_processed)
{
  Hermes::vector<RefinementSelectors::Selector *> refinement_selectors;
  RefinementSelectors::HOnlySelector selector;
  for (int i = 0; i < this->num; i++)
    refinement_selectors.push_back(&selector);

  return Adapt::adapt(refinement_selectors, thr, strat, regularize, to_be_processed);
}

void KellyTypeAdapt::add_error_estimator_vol( int i,
                                              WeakForm::error_vector_form_val_t vfv, WeakForm::error_vector_form_ord_t vfo,
                                              int area,
                                              Hermes::vector<MeshFunction*> ext )
{
  error_if(i < 0 || i >= this->num,
           "Invalid component number (%d), max. supported components: %d", i, H2D_MAX_COMPONENTS);
  error_if(area != HERMES_ANY && area < 0,
           "Invalid area number.");

  KellyTypeAdapt::ErrorEstimatorForm form = { i, area, vfv, vfo, ext };
  this->error_estimators_vol.push_back(form);
}

void KellyTypeAdapt::add_error_estimator_surf(int i,
                                              WeakForm::error_vector_form_val_t vfv, WeakForm::error_vector_form_ord_t vfo,
                                              int area,
                                              Hermes::vector<MeshFunction*> ext)
{
  error_if (i < 0 || i >= this->num,
            "Invalid equation number.");
  error_if (area != HERMES_ANY && area != H2D_DG_INNER_EDGE && area < 0,
            "Invalid area number.");

  KellyTypeAdapt::ErrorEstimatorForm form = { i, area, vfv, vfo, ext };
  this->error_estimators_surf.push_back(form);
}

double KellyTypeAdapt::calc_err_internal(Hermes::vector< Solution* > slns,
                                         Hermes::vector< double >* component_errors,
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

  // Prepare mesh traversal and error arrays.
  Mesh **meshes = new Mesh* [num];
  std::set<Transformable*> trset;

  num_act_elems = 0;
  for (int i = 0; i < num; i++)
  {
    meshes[i] = sln[i]->get_mesh();
    trset.insert(sln[i]);

    num_act_elems += meshes[i]->get_num_active_elements();
    int max = meshes[i]->get_max_element_id();

    if (errors[i] != NULL) delete [] errors[i];
    errors[i] = new double[max];
    memset(errors[i], 0.0, sizeof(double) * max);
  }

  for (unsigned int i = 0; i < error_estimators_vol.size(); i++)
    trset.insert(error_estimators_vol[i].ext.begin(), error_estimators_vol[i].ext.end());
  for (unsigned int i = 0; i < error_estimators_surf.size(); i++)
    trset.insert(error_estimators_surf[i].ext.begin(), error_estimators_surf[i].ext.end());

  Transformable **tr = new Transformable* [trset.size()];
  std::copy(trset.begin(), trset.end(), tr);
  trset.clear();

  double total_norm = 0.0;

  bool calc_norm = false;
  if ((error_flags & HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL ||
      (error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL) calc_norm = true;

  double *norms;
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
      for_all_active_elements(e, meshes[i])
        e->visited = false;
    }
  }

  // Begin the multimesh traversal.
  trav.begin(num, meshes, tr);
  while ((ee = trav.get_next_state(bnd, surf_pos)) != NULL)
  {
    // Go through all solution components.
    // WARNING: This may not work yet because of the untested multimesh behavior of NeighborSearch.
    for (int i = 0; i < num; i++)
    {
      // Set maximum integration order for use in integrals, see limit_order()
      update_limit_table(ee[i]->get_mode());

      NeighborSearch *nbs = new NeighborSearch(ee[i], meshes[i]);

      RefMap *rm = sln[i]->get_refmap();
      nbs->attach_rm(rm); // This is needed for querying the geometric and quadrature data.

      double err = 0.0;

      // Go through all volumetric error estimators.
      for (unsigned int iest = 0; iest < error_estimators_vol.size(); iest++)
      {
        // Skip current error estimator if it is assigned to a different component or geometric area
        // different from that of the current active element.
        if (error_estimators_vol[iest].i != i)
          continue;
        if (error_estimators_vol[iest].area >= 0 && error_estimators_vol[iest].area != ee[i]->marker)
          continue;
        else if (error_estimators_vol[iest].area != HERMES_ANY)
          continue;

        err += eval_volumetric_estimator(&error_estimators_vol[iest], rm);
      }

      // Go through all surface error estimators (includes both interface and boundary est's).
      for (unsigned int iest = 0; iest < error_estimators_surf.size(); iest++)
      {
        if (error_estimators_surf[iest].i != i)
          continue;

        for (int isurf = 0; isurf < ee[i]->get_num_surf(); isurf++)
        {
          if (error_estimators_surf[iest].area > 0 &&
              error_estimators_surf[iest].area != surf_pos[isurf].marker) continue;

          if (bnd[isurf])   // Boundary
          {
            if (error_estimators_surf[iest].area < 0 &&
                error_estimators_surf[iest].area != HERMES_ANY) continue;

            err += eval_boundary_estimator(&error_estimators_surf[iest], rm, surf_pos);
          }
          else              // Interface
          {
            if (error_estimators_surf[iest].area != H2D_DG_INNER_EDGE) continue;

            nbs->set_active_edge(isurf, ignore_visited_segments);

            // Go through all segments of the currently processed interface (segmentation is caused
            // by hanging nodes on the other side of the interface).
            for (int neighbor = 0; neighbor < nbs->get_num_neighbors(); neighbor++)
            {
              bool use_extended_shapeset = false;
              bool needs_processing = nbs->set_active_segment(neighbor, use_extended_shapeset);
              if (!needs_processing) continue;

              // The estimate is multiplied by 0.5 in order to distribute the error equally onto
              // the two neighboring elements.
              double central_err = 0.5 * eval_interface_estimator(&error_estimators_surf[iest],
                                                                  rm, surf_pos, nbs);
              double neighb_err = central_err;

              // Scale the error estimate by the scaling function dependent on the element diameter
              // (use the central element's diameter).
              if (use_aposteriori_interface_scaling && interface_scaling_fns[i])
                central_err *= interface_scaling_fns[i](ee[i]->get_diameter());

              // In the case this edge will be ignored when calculating the error for the element on
              // the other side, add the now computed error to that element as well.
              if (ignore_visited_segments)
              {
                Element *neighb = nbs->get_current_neighbor_element();

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
            }
          }
        }
      }

      if (calc_norm)
      {
        double nrm = eval_solution_norm(error_form[i][i], error_ord[i][i], rm, sln[i]);
        norms[i] += nrm;
        total_norm += nrm;
      }

      errors_components[i] += err;
      total_error += err;
      errors[i][ee[i]->id] += err;

      ee[i]->visited = true;

      delete nbs;

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
      for_all_active_elements(e, meshes[i])
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
  fill_regular_queue(meshes);
  have_errors = true;

  delete [] meshes;
  delete [] tr;
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

double KellyTypeAdapt::eval_solution_norm(error_matrix_form_val_t val, error_matrix_form_ord_t ord, RefMap *rm, MeshFunction* sln)
{
  // determine the integration order
  int inc = (sln->get_num_components() == 2) ? 1 : 0;
  Func<Ord>* ou = init_fn_ord(sln->get_fn_order() + inc);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = ord(1, &fake_wt, NULL, ou, ou, fake_e, NULL);
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
  Func<scalar>* u = init_fn(sln, rm, order);
  scalar res = val(np, jwt, NULL, u, u, e, NULL);

  e->free(); delete e;
  delete [] jwt;
  u->free_fn(); delete u;

  return std::abs(res);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// COPIED FROM discrete_problem.cpp

// Initialize integration order for external functions
ExtData<Ord>* init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext)
{
  _F_
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = init_fn_ord(ext[i]->get_fn_order());
  fake_ext->fn = fake_ext_fn;

  return fake_ext;
}

// Initialize external functions (obtain values, derivatives,...)
ExtData<scalar>* init_ext_fns(Hermes::vector<MeshFunction *> &ext, RefMap *rm, const int order)
{
  _F_
  ExtData<scalar>* ext_data = new ExtData<scalar>;
  Func<scalar>** ext_fn = new Func<scalar>*[ext.size()];
  for (unsigned i = 0; i < ext.size(); i++) {
    if (ext[i] != NULL) ext_fn[i] = init_fn(ext[i], rm, order);
    else ext_fn[i] = NULL;
  }
  ext_data->nf = ext.size();
  ext_data->fn = ext_fn;

  return ext_data;
}

// Initialize integration order on a given edge for external functions
ExtData<Ord>* init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext, int edge)
{
  _F_
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = init_fn_ord(ext[i]->get_edge_fn_order(edge));
  fake_ext->fn = fake_ext_fn;

  return fake_ext;
}

// Initialize discontinuous external functions (obtain values, derivatives,... on both sides of the
// supplied NeighborSearch's active edge).
ExtData<scalar>* init_ext_fns(Hermes::vector<MeshFunction *> &ext, NeighborSearch* nbs)
{
  Func<scalar>** ext_fns = new Func<scalar>*[ext.size()];
  for(unsigned int j = 0; j < ext.size(); j++)
    ext_fns[j] = nbs->init_ext_fn(ext[j]);

  ExtData<scalar>* ext_data = new ExtData<scalar>;
  ext_data->fn = ext_fns;
  ext_data->nf = ext.size();

  return ext_data;
}

// Initialize integration order for discontinuous external functions.
ExtData<Ord>* init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext, NeighborSearch* nbs)
{
  Func<Ord>** fake_ext_fns = new Func<Ord>*[ext.size()];
  for (unsigned int j = 0; j < ext.size(); j++)
    fake_ext_fns[j] = nbs->init_ext_fn_ord(ext[j]);

  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->fn = fake_ext_fns;
  fake_ext->nf = ext.size();

  return fake_ext;
}

// END COPIED FROM discrete_problem.cpp
//////////////////////////////////////////////////////////////////////////////////////////////////////////

double KellyTypeAdapt::eval_volumetric_estimator(KellyTypeAdapt::ErrorEstimatorForm* err_est_form, RefMap *rm)
{
  // determine the integration order
  int inc = (this->sln[err_est_form->i]->get_num_components() == 2) ? 1 : 0;

  Func<Ord>** oi = new Func<Ord>* [num];
  for (int i = 0; i < num; i++)
    oi[i] = init_fn_ord(this->sln[i]->get_fn_order() + inc);

  // Order of additional external functions.
  ExtData<Ord>* fake_ext = init_ext_fns_ord(err_est_form->ext);

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
    ui[i] = init_fn(this->sln[i], rm, order);
  ExtData<scalar>* ext = init_ext_fns(err_est_form->ext, rm, order);

  scalar res = volumetric_scaling_const *
                err_est_form->fn(np, jwt, ui, ui[err_est_form->i], e, ext);

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
  ExtData<Ord>* fake_ext = init_ext_fns_ord(err_est_form->ext, surf_pos->surf_num);

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
    ui[i] = init_fn(this->sln[i], rm, eo);
  ExtData<scalar>* ext = init_ext_fns(err_est_form->ext, rm, eo);

  scalar res = boundary_scaling_const *
                err_est_form->fn(np, jwt, ui, ui[err_est_form->i], e, ext);

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
                                                RefMap *rm, SurfPos* surf_pos, NeighborSearch* nbs)
{
  // Determine integration order.
  Func<Ord>** oi = new Func<Ord>* [num];
  for (int i = 0; i < num; i++)
    oi[i] = nbs->init_ext_fn_ord(this->sln[i]);

  // Order of additional external functions.
  ExtData<Ord>* fake_ext = init_ext_fns_ord(err_est_form->ext, nbs);

  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  Element *neighb_el = nbs->get_current_neighbor_element();
  Geom<Ord>* fake_e = new InterfaceGeom<Ord>(init_geom_ord(), neighb_el->marker, neighb_el->id, neighb_el->get_diameter());
  double fake_wt = 1.0;
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
  nbs->set_quad_order(order);

  // Init geometry and jacobian*weights (do not use the NeighborSearch caching mechanism).
  Geom<double>* e = nbs->init_geometry(NULL, surf_pos);
  double* jwt = nbs->init_jwt(NULL);

  // function values
  Func<scalar>** ui = new Func<scalar>* [num];
  for (int i = 0; i < num; i++)
    ui[i] = nbs->init_ext_fn(this->sln[i]);
  ExtData<scalar>* ext = init_ext_fns(err_est_form->ext, nbs);

  scalar res = interface_scaling_const *
                err_est_form->fn(nbs->get_quad_np(), jwt, ui, ui[err_est_form->i], e, ext);

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
