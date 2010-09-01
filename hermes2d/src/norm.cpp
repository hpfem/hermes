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

#include "common.h"
#include "norm.h"
#include "limit_order.h"
#include "discrete_problem.h"
#include "refmap.h"
#include "integrals_h1.h"
#include "traverse.h"
#include "feproblem.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////

double calc_abs_error(MeshFunction* sln1, MeshFunction* sln2, int norm_type)
{
  double error;
  switch (norm_type) {
  case H2D_L2_NORM: 
    error = calc_abs_error(error_fn_l2, sln1, sln2);
    break;
  case H2D_H1_NORM: 
    error = calc_abs_error(error_fn_h1, sln1, sln2);
    break;
  case H2D_HCURL_NORM: 
    error = calc_abs_error(error_fn_hc, sln1, sln2);
    break;
  case H2D_HDIV_NORM: 
    error = calc_abs_error(error_fn_hdiv, sln1, sln2);
    break;
  default: error("Unknown norm in calc_error().");
  }

  return error;
}

double calc_norm(MeshFunction* ref_sln, int norm_type)
{
  double norm;
  switch (norm_type) {
  case H2D_L2_NORM: 
    norm = calc_norm(norm_fn_l2, ref_sln);
    break;
  case H2D_H1_NORM: 
    norm = calc_norm(norm_fn_h1, ref_sln);
    break;
  case H2D_HCURL_NORM: 
    norm = calc_norm(norm_fn_hc, ref_sln);
    break;
  case H2D_HDIV_NORM: 
    norm = calc_norm(norm_fn_hdiv, ref_sln);
    break;
  default: error("Unknown norm in calc_norm().");
  }

  return norm;
}

/// Calculates the absolute error between sln1 and sln2 using function fn
double calc_abs_error(double (*fn)(MeshFunction*, MeshFunction*, RefMap*, RefMap*), MeshFunction* sln1, 
                      MeshFunction* sln2)
{
  // sanity checks
  if (fn == NULL) error("error norm function is NULL in calc_abs_error().");
  if (sln1 == NULL) error("sln1 is NULL in calc_abs_error().");
  if (sln2 == NULL) error("sln2 is NULL in calc_abs_error().");

  Quad2D* quad = &g_quad_2d_std;
  sln1->set_quad_2d(quad);
  sln2->set_quad_2d(quad);

  Mesh* meshes[2] = { sln1->get_mesh(), sln2->get_mesh() };
  Transformable* tr[2] = { sln1, sln2 };
  Traverse trav;
  trav.begin(2, meshes, tr);

  double error = 0.0;
  Element** ee;
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
  {
    update_limit_table(ee[0]->get_mode());

    RefMap* ru = sln1->get_refmap();
    RefMap* rv = sln2->get_refmap();

    error += fn(sln1, sln2, ru, rv);
  }
  trav.finish();
  return sqrt(error);
}


/// Calculates the norm of sln using function fn
double calc_norm(double (*fn)(MeshFunction*, RefMap*), MeshFunction* sln)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);

  double norm = 0.0;
  Element* e;
  Mesh* mesh = sln->get_mesh();

  for_all_active_elements(e, mesh)
  {
    // set maximum integration order for use in integrals, see limit_order()
    update_limit_table(e->get_mode());

    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();

    norm += fn(sln, ru);
  }
  return sqrt(norm);
}

double calc_rel_error(MeshFunction* sln, MeshFunction* ref_sln, int norm_type)
{
  double error = calc_abs_error(sln, ref_sln, norm_type);
  double norm = calc_norm(ref_sln, norm_type);
  
  return error/norm;
}


//// H1 space //////////////////////////////////////////////////////////////////////////////////////

// function used to calculate error in H1 norm
double error_fn_h1(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2*std::max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o);
  sln2->set_quad_order(o);

  scalar *uval, *vval, *dudx, *dudy, *dvdx, *dvdy;
  uval = sln1->get_fn_values();
  vval = sln2->get_fn_values();
  sln1->get_dx_dy_values(dudx, dudy);
  sln2->get_dx_dy_values(dvdx, dvdy);

  double result = 0.0;
  h1_integrate_expression(sqr(uval[i] - vval[i]) +
                          sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]));
  return result;
}

// function used to calculate H1 norm of the solution
double norm_fn_h1(MeshFunction* sln, RefMap* ru)
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  scalar *uval, *dudx, *dudy;
  uval = sln->get_fn_values();
  sln->get_dx_dy_values(dudx, dudy);

  double result = 0.0;
  h1_integrate_expression(sqr(uval[i]) + sqr(dudx[i]) + sqr(dudy[i]));
  return result;
}

// function used to calculate error in L2 norm
double error_fn_l2(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2*std::max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o, H2D_FN_VAL);
  sln2->set_quad_order(o, H2D_FN_VAL);

  scalar *uval, *vval;
  uval = sln1->get_fn_values();
  vval = sln2->get_fn_values();

  double result = 0.0;
  h1_integrate_expression(sqr(uval[i] - vval[i]));
  return result;
}


// function used to calculate L2 norm of the solution
double norm_fn_l2(MeshFunction* sln, RefMap* ru)
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 *sln->get_fn_order() + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o, H2D_FN_VAL);

  scalar* uval = sln->get_fn_values();

  double result = 0.0;
  h1_integrate_expression(sqr(uval[i]));
  return result;
}

//// Hcurl space ///////////////////////////////////////////////////////////////////////////////////

// function used to calculate error in Hcurl norm
double error_fn_hc(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2 * std::max(sln1->get_fn_order(), sln2->get_fn_order()) + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o);
  sln2->set_quad_order(o);


  scalar *uval0 = sln1->get_fn_values(0), *uval1 = sln1->get_fn_values(1);
  scalar *udx1  = sln1->get_dx_values(1), *udy0  = sln1->get_dy_values(0);
  scalar *vval0 = sln2->get_fn_values(0), *vval1 = sln2->get_fn_values(1);
  scalar *vdx1  = sln2->get_dx_values(1), *vdy0  = sln2->get_dy_values(0);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i] - vval0[i]) + sqr(uval1[i] - vval1[i]) +
                          sqr((udx1[i] - udy0[i]) - (vdx1[i] - vdy0[i])));
  return result;
}


// function used to calculate Hcurl norm
double norm_fn_hc(MeshFunction* sln, RefMap* ru)
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  scalar *uval0 = sln->get_fn_values(0), *uval1 = sln->get_fn_values(1);
  scalar *udx1  = sln->get_dx_values(1), *udy0  = sln->get_dy_values(0);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i]) + sqr(uval1[i]) + sqr(udx1[i] - udy0[i]));
  return result;
}

// function used to calculate error in Hcurl norm
double error_fn_hcl2(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2 * std::max(sln1->get_fn_order(), sln2->get_fn_order()) + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o);
  sln2->set_quad_order(o);


  scalar *uval0 = sln1->get_fn_values(0), *uval1 = sln1->get_fn_values(1);
  scalar *vval0 = sln2->get_fn_values(0), *vval1 = sln2->get_fn_values(1);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i] - vval0[i]) + sqr(uval1[i] - vval1[i]));
  return result;
}


// function used to calculate Hcurl norm
double norm_fn_hcl2(MeshFunction* sln, RefMap* ru)
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  scalar *uval0 = sln->get_fn_values(0), *uval1 = sln->get_fn_values(1);
  scalar *udx1  = sln->get_dx_values(1), *udy0  = sln->get_dy_values(0);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i]) + sqr(uval1[i]));
  return result;
}


double hcurl_l2error(MeshFunction* sln1, MeshFunction* sln2)
{
  double error = calc_abs_error(error_fn_hcl2, sln1, sln2);
  double norm = calc_norm(norm_fn_hcl2, sln2);
  return error / norm;
}


double hcurl_l2norm(MeshFunction* sln)
{
  return calc_norm(norm_fn_hcl2, sln);
}
//// Hdiv space ///////////////////////////////////////////////////////////////////////////////////

// function used to calculate error in Hcurl norm
double error_fn_hdiv(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv)
{
  error("error_fn_hdiv() not implemented yet.");

  // Hcurl code
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2 * std::max(sln1->get_fn_order(), sln2->get_fn_order()) + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o);
  sln2->set_quad_order(o);


  scalar *uval0 = sln1->get_fn_values(0), *uval1 = sln1->get_fn_values(1);
  scalar *udx1  = sln1->get_dx_values(1), *udy0  = sln1->get_dy_values(0);
  scalar *vval0 = sln2->get_fn_values(0), *vval1 = sln2->get_fn_values(1);
  scalar *vdx1  = sln2->get_dx_values(1), *vdy0  = sln2->get_dy_values(0);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i] - vval0[i]) + sqr(uval1[i] - vval1[i]) +
                          sqr((udx1[i] - udy0[i]) - (vdx1[i] - vdy0[i])));
  return result;
}


// function used to calculate Hcurl norm
double norm_fn_hdiv(MeshFunction* sln, RefMap* ru)
{
  error("norm_fn_hdiv() not implemented yet.");

  // Hcurl code
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  scalar *uval0 = sln->get_fn_values(0), *uval1 = sln->get_fn_values(1);
  scalar *udx1  = sln->get_dx_values(1), *udy0  = sln->get_dy_values(0);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i]) + sqr(uval1[i]) + sqr(udx1[i] - udy0[i]));
  return result;
}

