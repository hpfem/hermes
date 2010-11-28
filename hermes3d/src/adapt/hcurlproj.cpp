// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "../function.h"
#include "../solution.h"
#include "hcurlproj.h"
#include "../quad.h"
#include "../refdomain.h"
#include "../transform.h"
#include "../shapeset/common.h"
#include "../shapeset/lobatto.h"
#include "../../../hermes_common/callstack.h"
#include "../../../hermes_common/matrix.h"

#ifdef DEBUG
	#define PRINTF			printf
#else
	#define PRINTF(...)
#endif


std::map<HCurlProjection::Key, double*, HCurlProjection::Compare> HCurlProjection::proj_mat_entries;

HCurlProjection::HCurlProjection(Solution *afn, Element *e, Shapeset *ss) : Projection(afn, e, ss)
{
}

double HCurlProjection::get_error(int split, int son, const Ord3 &order)
{
	_F_
	sln->enable_transform(false);

	Ord3 order_rhs = order;

	calc_projection(split, son + 1, order_rhs);

	QuadPt3D *pt = quad->get_points(order_rhs);
	int np = quad->get_num_points(order_rhs);

	double error = 0.0;
	for (int i = 0; i < int_ns[split]; i++) {
		Trf *tr = get_trf(int_trf[split][i]);

    unsigned int son_idx = base_elem->get_son(int_son[son + 1][i]);
	  sln->set_active_element(mesh->elements[son_idx]);
	  sln->precalculate(np, pt, FN_DEFAULT);
	  
    scalar *rval0 = sln->get_fn_values(0);
    scalar *rval1 = sln->get_fn_values(1);
    scalar *rval2 = sln->get_fn_values(2);

    // Derivatives for the calculation of curl values.
    scalar *drdx0, *drdy0, *drdz0;
    scalar *drdx1, *drdy1, *drdz1;
    scalar *drdx2, *drdy2, *drdz2;
    
    sln->get_dx_dy_dz_values(drdx0, drdy0, drdz0, 0);
    sln->get_dx_dy_dz_values(drdx1, drdy1, drdz1, 1);
    sln->get_dx_dy_dz_values(drdx2, drdy2, drdz2, 2);

    scalar *rcurl0 = new scalar[np];
    scalar *rcurl1 = new scalar[np];
    scalar *rcurl2 = new scalar[np];

    for (int point_num = 0; point_num < np; point_num++)
    {
      rcurl0[point_num] = drdy2[point_num] - drdz1[point_num];
      rcurl1[point_num] = drdz0[point_num] - drdx2[point_num];
      rcurl2[point_num] = drdx1[point_num] - drdy0[point_num];
    }

	  QuadPt3D *tpt = new QuadPt3D[np];
	  transform_points(np, pt, tr, tpt);
	  scalar *prfn0 = new scalar[np];
    scalar *prdx0 = new scalar[np];
    scalar *prdy0 = new scalar[np];
    scalar *prdz0 = new scalar[np];

	  memset(prfn0, 0, np * sizeof(double));
	  memset(prdx0, 0, np * sizeof(double));
	  memset(prdy0, 0, np * sizeof(double));
	  memset(prdz0, 0, np * sizeof(double));

    scalar *prfn1 = new scalar[np];
    scalar *prdx1 = new scalar[np];
    scalar *prdy1 = new scalar[np];
    scalar *prdz1 = new scalar[np];

	  memset(prfn1, 0, np * sizeof(double));
	  memset(prdx1, 0, np * sizeof(double));
	  memset(prdy1, 0, np * sizeof(double));
	  memset(prdz1, 0, np * sizeof(double));

    scalar *prfn2 = new scalar[np];
    scalar *prdx2 = new scalar[np];
    scalar *prdy2 = new scalar[np];
    scalar *prdz2 = new scalar[np];

	  memset(prfn2, 0, np * sizeof(double));
	  memset(prdx2, 0, np * sizeof(double));
	  memset(prdy2, 0, np * sizeof(double));
	  memset(prdz2, 0, np * sizeof(double));

    scalar *prcurl0 = new scalar[np];
    scalar *prcurl1 = new scalar[np];
    scalar *prcurl2 = new scalar[np];

	  for (int j = 0; j < n_fns; j++) {
#ifndef H3D_COMPLEX
		  double *tmp = new double[np];
		  ss->get_fn_values(fn_idx[j], np, tpt, 0, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prfn0, 1);
		  ss->get_dx_values(fn_idx[j], np, tpt, 0, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prdx0, 1);
		  ss->get_dy_values(fn_idx[j], np, tpt, 0, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prdy0, 1);
		  ss->get_dz_values(fn_idx[j], np, tpt, 0, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prdz0, 1);
#else
		  double *tmp = new double[np];
		  scalar *sctmp = new scalar[np];
		  ss->get_fn_values(fn_idx[j], np, tpt, 0, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prfn0, 1);
		  ss->get_dx_values(fn_idx[j], np, tpt, 0, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prdx0, 1);
		  ss->get_dy_values(fn_idx[j], np, tpt, 0, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prdy0, 1);
		  ss->get_dz_values(fn_idx[j], np, tpt, 0, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prdz0, 1);
#endif

#ifndef H3D_COMPLEX
		  memset(tmp, 0, np*sizeof(double));
		  ss->get_fn_values(fn_idx[j], np, tpt, 1, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prfn1, 1);
		  ss->get_dx_values(fn_idx[j], np, tpt, 1, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prdx1, 1);
		  ss->get_dy_values(fn_idx[j], np, tpt, 1, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prdy1, 1);
		  ss->get_dz_values(fn_idx[j], np, tpt, 1, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prdz1, 1);
#else
		  memset(tmp, 0, np*sizeof(double));
		  memset(sctmp, 0, np*sizeof(scalar));
		  ss->get_fn_values(fn_idx[j], np, tpt, 1, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prfn1, 1);
		  ss->get_dx_values(fn_idx[j], np, tpt, 1, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prdx1, 1);
		  ss->get_dy_values(fn_idx[j], np, tpt, 1, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prdy1, 1);
		  ss->get_dz_values(fn_idx[j], np, tpt, 1, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prdz1, 1);
#endif

#ifndef H3D_COMPLEX
		  memset(tmp, 0, np*sizeof(double));
		  ss->get_fn_values(fn_idx[j], np, tpt, 2, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prfn2, 1);
		  ss->get_dx_values(fn_idx[j], np, tpt, 2, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prdx2, 1);
		  ss->get_dy_values(fn_idx[j], np, tpt, 2, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prdy2, 1);
		  ss->get_dz_values(fn_idx[j], np, tpt, 2, tmp);
		  blas_axpy(np, proj_coef[j], tmp, 1, prdz2, 1);
#else
		  memset(tmp, 0, np*sizeof(double));
		  memset(sctmp, 0, np*sizeof(scalar));
		  ss->get_fn_values(fn_idx[j], np, tpt, 2, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prfn2, 1);
		  ss->get_dx_values(fn_idx[j], np, tpt, 2, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prdx2, 1);
		  ss->get_dy_values(fn_idx[j], np, tpt, 2, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prdy2, 1);
		  ss->get_dz_values(fn_idx[j], np, tpt, 2, tmp);
		  for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
		  blas_axpy(np, proj_coef[j], sctmp, 1, prdz2, 1);
      delete [] tmp;
      delete [] sctmp;
#endif
	  }
	  for (int k = 0; k < np; k++)
    {
      prcurl0[k] = prdy2[k] - prdz1[k];
      prcurl1[k] = prdz0[k] - prdx2[k];
      prcurl2[k] = prdx1[k] - prdy0[k];

		  error += pt[k].w *
			  (
          sqr(rcurl0[k] - prcurl0[k]) + 
          sqr(rcurl1[k] - prcurl1[k]) + 
          sqr(rcurl2[k] - prcurl2[k]) + 
          sqr(rval0[k] - prfn0[k]) +
          sqr(rval1[k] - prfn1[k]) +
          sqr(rval2[k] - prfn2[k])
         );
    }
    delete [] rcurl0;
    delete [] rcurl1;
    delete [] rcurl2;
    
    delete [] tpt;

    delete [] prfn0;
    delete [] prfn1;
    delete [] prfn2;

    delete [] prdx0;
    delete [] prdx1;
    delete [] prdx2;
    delete [] prdy0;
    delete [] prdy1;
    delete [] prdy2;
    delete [] prdz0;
    delete [] prdz1;
    delete [] prdz2;

    delete [] prcurl0;
    delete [] prcurl1;
    delete [] prcurl2;
	}

  sln->enable_transform(true);
	return error;
}

void HCurlProjection::calc_projection(int split, int son, Ord3 &order)
{
	_F_
  
  // Max order used for the use of the get_error() function.
  Ord3 max_order_used = order;

  // Number of edge functions.
  int n_edge_fns = 0;
  for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++)
    if (ss->get_num_edge_fns(order.get_edge_order(iedge)) > 0)
      n_edge_fns += ss->get_num_edge_fns(order.get_edge_order(iedge));
  
  // Number of face functions.
  int n_face_fns = 0;
  for (int iface = 0; iface < Hex::NUM_FACES; iface++)
	  if (ss->get_num_face_fns(order.get_face_order(iface)) > 0)
      n_face_fns +=ss->get_num_face_fns(order.get_face_order(iface));

  // Number of bubble functions.
  int n_bubble_fns = ss->get_num_bubble_fns(order);

  // Total and array reallocation.
  this->n_fns = n_bubble_fns + n_face_fns + n_edge_fns;
 	delete [] fn_idx;
	fn_idx = new int [n_fns];

  // Running index.
	int mm = 0;
	
	// edge functions
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
    const int *edge_fn_idx = ss->get_edge_indices(iedge, base_elem->get_edge_orientation(iedge), order.get_edge_order(iedge));
			for (int i = 0; i < ss->get_num_edge_fns(order.get_edge_order(iedge)); i++, mm++)
				fn_idx[mm] = edge_fn_idx[i];
	}

	// face functions
	for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
    const int *face_fn_idx = ss->get_face_indices(iface, base_elem->get_face_orientation(iface), order.get_face_order(iface));
			for (int i = 0; i < ss->get_num_face_fns(order.get_face_order(iface)); i++, mm++)
				fn_idx[mm] = face_fn_idx[i];
	}

  // bubble functions
  if (n_bubble_fns > 0) {
	  const int *bubble_fn_idx = ss->get_bubble_indices(order);
	  for (int i = 0; i < n_bubble_fns; i++, mm++)
		  fn_idx[mm] = bubble_fn_idx[i];
  }

	double **proj_mat = new_matrix<double>(n_fns, n_fns);
	scalar *proj_rhs = new scalar[n_fns];
	memset(proj_rhs, 0, sizeof(scalar) * n_fns);

  double d;
	int * iperm = new int[n_fns];

  // proj matrix
  for (int i = 0; i < n_fns; i++) {
    for (int j = 0; j < n_fns; j++) {

      // Looking for a value in the cache.
      Key key(fn_idx[i], fn_idx[j]);
      bool calculation_needed = true;
      if (proj_mat_entries[key] != NULL) {
        proj_mat[i][j] += *proj_mat_entries[key];
        calculation_needed = false;
      }

      // Also inverted key can be used.
      if(calculation_needed)
      {
        Key key_inverted(fn_idx[j], fn_idx[i]);
        if (proj_mat_entries[key_inverted] != NULL) {
          proj_mat[i][j] += *proj_mat_entries[key_inverted];
          calculation_needed = false;
        }
      }

      if(calculation_needed) {
        fu->set_active_shape(fn_idx[i]);

        // Calculating of necessary order + updating max order used.
        Ord3 order_rhs = ss->get_order(fn_idx[i]) + ss->get_order(fn_idx[j]) + order;
        if(order_rhs.x > max_order_used.x || order_rhs.y > max_order_used.y || order_rhs.z > max_order_used.z)
          max_order_used = order_rhs;

        QuadPt3D *pt = quad->get_points(order_rhs);
        int np = quad->get_num_points(order_rhs);

        fu->precalculate(np, pt, FN_DEFAULT);

        // Function values.
        double *uval0 = fu->get_fn_values(0);
        double *uval1 = fu->get_fn_values(1);
        double *uval2 = fu->get_fn_values(2);
        
        // Derivatives for the calculation of curl values.
        double *dudx0, *dudy0, *dudz0;
        double *dudx1, *dudy1, *dudz1;
        double *dudx2, *dudy2, *dudz2;

        fu->get_dx_dy_dz_values(dudx0, dudy0, dudz0, 0);
        fu->get_dx_dy_dz_values(dudx1, dudy1, dudz1, 1);
        fu->get_dx_dy_dz_values(dudx2, dudy2, dudz2, 2);

        double *ucurl0 = new double[np];
        double *ucurl1 = new double[np];
        double *ucurl2 = new double[np];

        fv->set_active_shape(fn_idx[j]);

        // Function values.
	      fv->precalculate(np, pt, FN_DEFAULT);
        double *vval0 = fv->get_fn_values(0);
	      double *vval1 = fv->get_fn_values(1);
	      double *vval2 = fv->get_fn_values(2);

        // Derivatives for the calculation of curl values.
        double *dvdx0, *dvdy0, *dvdz0;
        double *dvdx1, *dvdy1, *dvdz1;
        double *dvdx2, *dvdy2, *dvdz2;
	      
        fv->get_dx_dy_dz_values(dvdx0, dvdy0, dvdz0, 0);
        fv->get_dx_dy_dz_values(dvdx1, dvdy1, dvdz1, 1);
        fv->get_dx_dy_dz_values(dvdx2, dvdy2, dvdz2, 2);
        
        double *vcurl0 = new double[np];
        double *vcurl1 = new double[np];
        double *vcurl2 = new double[np];

        for (int point_num = 0; point_num < np; point_num++)
        {
          ucurl0[point_num] = dudy2[point_num] - dudz1[point_num];
          ucurl1[point_num] = dudz0[point_num] - dudx2[point_num];
          ucurl2[point_num] = dudx1[point_num] - dudy0[point_num];

          vcurl0[point_num] = dvdy2[point_num] - dvdz1[point_num];
          vcurl1[point_num] = dvdz0[point_num] - dvdx2[point_num];
          vcurl2[point_num] = dvdx1[point_num] - dvdy0[point_num];
        }

	      double value = 0.0;
	      for (int k = 0; k < np; k++) {
		      value += pt[k].w * 
            (
              (ucurl0[k] * conj(vcurl0[k]) + ucurl1[k] * conj(vcurl1[k]) + ucurl2[k] * conj(vcurl2[k])) + 
              (uval0[k] * conj(vval0[k]) + uval1[k] * conj(vval1[k]) + uval2[k] * conj(vval2[k]))
            );
        }
	      proj_mat[i][j] += value;
        
        proj_mat_entries[key] = new double(value);

        delete [] vcurl0;
        delete [] vcurl1;
        delete [] vcurl2;

        delete [] ucurl0;
        delete [] ucurl1;
        delete [] ucurl2;
      }
    }
  }

  // rhs
  for (int e = 0; e < int_ns[split]; e++) {
	  sln->set_active_element(mesh->elements[base_elem->get_son(int_son[son][e])]);
    Trf *tr = get_trf(int_trf[split][e]);
	  for (int i = 0; i < n_fns; i++) {
		  int iidx = fn_idx[i];
		  fu->set_active_shape(iidx);

		  Ord3 order_rhs = ss->get_order(iidx) + order;
		  QuadPt3D *pt = quad->get_points(order_rhs);
		  int np = quad->get_num_points(order_rhs);

		  if (int_trf[split][e] != -1) fu->push_transform(int_trf[split][e]);

      // Function values.
	    fu->precalculate(np, pt, FN_DEFAULT);
	    sln->precalculate(np, pt, FN_DEFAULT);

	    double *uval0 = fu->get_fn_values(0);
	    double *uval1 = fu->get_fn_values(1);
	    double *uval2 = fu->get_fn_values(2);
      scalar *rval0 = sln->get_fn_values(0);
	    scalar *rval1 = sln->get_fn_values(1);
	    scalar *rval2 = sln->get_fn_values(2);

      // Derivatives for the calculation of curl values.
	    double *dudx0, *dudy0, *dudz0;
      double *dudx1, *dudy1, *dudz1;
      double *dudx2, *dudy2, *dudz2;
      scalar *drdx0, *drdy0, *drdz0;
      scalar *drdx1, *drdy1, *drdz1;
      scalar *drdx2, *drdy2, *drdz2;
	    
      fu->get_dx_dy_dz_values(dudx0, dudy0, dudz0, 0);
      fu->get_dx_dy_dz_values(dudx1, dudy1, dudz1, 1);
      fu->get_dx_dy_dz_values(dudx2, dudy2, dudz2, 2);
      sln->get_dx_dy_dz_values(drdx0, drdy0, drdz0, 0);
      sln->get_dx_dy_dz_values(drdx1, drdy1, drdz1, 1);
      sln->get_dx_dy_dz_values(drdx2, drdy2, drdz2, 2);

      double *ucurl0 = new double [np];
      double *ucurl1 = new double [np];
      double *ucurl2 = new double [np];

      scalar *rcurl0 = new scalar [np];
      scalar *rcurl1 = new scalar [np];
      scalar *rcurl2 = new scalar [np];

      for (int point_num = 0; point_num < np; point_num++)
      {
        ucurl0[point_num] = dudy2[point_num] - dudz1[point_num];
        ucurl1[point_num] = dudz0[point_num] - dudx2[point_num];
        ucurl2[point_num] = dudx1[point_num] - dudy0[point_num];
      
        rcurl0[point_num] = drdy2[point_num] - drdz1[point_num];
        rcurl1[point_num] = drdz0[point_num] - drdx2[point_num];
        rcurl2[point_num] = drdx1[point_num] - drdy0[point_num];
      }

		  scalar value = 0.0;
		  for (int k = 0; k < np; k++) {
			   value += pt[k].w * 
            (
              (ucurl0[k] * conj(rcurl0[k]) + ucurl1[k] * conj(rcurl1[k]) + ucurl2[k] * conj(rcurl2[k])) + 
              (uval0[k] * conj(rval0[k]) + uval1[k] * conj(rval1[k]) + uval2[k] * conj(rval2[k]))
            );
		  }
		  proj_rhs[i] += value * (1 / (double) int_ns[split]);

		  if (int_trf[split][e] != -1) fu->pop_transform();

      delete [] ucurl0;
      delete [] ucurl1;
      delete [] ucurl2;
      delete [] rcurl0;
      delete [] rcurl1;
      delete [] rcurl2;
	  }
  }

  ludcmp(proj_mat, n_fns, iperm, &d);
  lubksb(proj_mat, n_fns, iperm, proj_rhs);

  proj_coef = new double [n_fns];
  memcpy(proj_coef, proj_rhs, n_fns * sizeof(double));

  delete [] proj_mat;
	delete [] proj_rhs;
  delete [] iperm;

  order = max_order_used;
  return;
}
