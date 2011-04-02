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

#include "../h2d_common.h"
#include "../shapeset/shapeset_h1_all.h"
#include "../shapeset/shapeset_common.h"
#include "../shapeset/precalc.h"
#include "curved.h"
#include "mesh.h"
#include "../quadrature/quad_all.h"
#include "../../../hermes_common/matrix.h"
  
H1ShapesetJacobi CurvMap::ref_map_shapeset;
PrecalcShapeset CurvMap::ref_map_pss(&ref_map_shapeset);

double** CurvMap::edge_proj_matrix;
double** CurvMap::bubble_proj_matrix_tri;
double** CurvMap::bubble_proj_matrix_quad;

double* CurvMap::edge_p;
double* CurvMap::bubble_tri_p;
double* CurvMap::bubble_quad_p;

Quad1DStd CurvMap::quad1d;
Quad2DStd CurvMap::quad2d;

Trf CurvMap::ctm;

//// NURBS //////////////////////////////////////////////////////////////////////////////////////////
// recursive calculation of the basis function N_i,k
bool CurvMap::warning_issued = false;

double CurvMap::nurbs_basis_fn(int i, int k, double t, double* knot)
{
  _F_
  if (k == 0)
  {
    return (t >= knot[i] && t <= knot[i+1] && knot[i] < knot[i+1]) ? 1.0 : 0.0;
  }
  else
  {
    double N1 = nurbs_basis_fn(i, k-1, t, knot);
    double N2 = nurbs_basis_fn(i+1, k-1, t, knot);

    double result = 0.0;
    if (knot[i+k] != knot[i])
    {
      result += ((t - knot[i]) / (knot[i+k] - knot[i])) * N1;
    }
    if (knot[i+k+1] != knot[i+1])
    {
      result += ((knot[i+k+1] - t) / (knot[i+k+1] - knot[i+1])) * N2;
    }
    return result;
  }
}

// Nurbs curve: t goes from -1 to 1, function returns x, y coordinates in plane
// as well as the unit normal and unit tangential vectors.
void CurvMap::nurbs_edge(Element* e, Nurbs* nurbs, int edge, double t, double& x, 
                         double& y, double& n_x, double& n_y, double& t_x, double& t_y)
{
  _F_
  t = (t + 1) / 2.0; // nurbs curves are parametrized from 0 to 1
  if (nurbs == NULL)
  {
    double2 v;
    v[0] = e->vn[e->next_vert(edge)]->x - e->vn[edge]->x;
    v[1] = e->vn[e->next_vert(edge)]->y - e->vn[edge]->y;
    x = e->vn[edge]->x + t * v[0];
    y = e->vn[edge]->y + t * v[1];
    double abs_v = sqrt(sqr(v[0]) + sqr(v[1]));
    t_x = v[0] / abs_v;
    t_y = v[1] / abs_v;
    n_x = t_y;
    n_y = -t_x;
  }
  else
  {
    double3* cp = nurbs->pt;
    x = y = 0.0;
    double sum = 0.0;  // sum of basis fns and weights

    for (int i = 0; i < nurbs->np; i++)
    {
      double basis = nurbs_basis_fn(i, nurbs->degree, t, nurbs->kv);
      sum += cp[i][2] * basis;
      x   += cp[i][2] * basis * cp[i][0];
      y   += cp[i][2] * basis * cp[i][1];
    }

    sum = 1.0 / sum;
    x *= sum;
    y *= sum;

    if(!warning_issued) {
      printf("FIXME: IMPLEMENT CALCULATION OF n_x, n_y, t_x, t_y in nurbs_edge() !!!\n");
      warning_issued = true;
    }
    n_x = 0;
    n_y = 0;
    t_x = 0;
    t_y = 0;
  }
}

//// non-polynomial reference map //////////////////////////////////////////////////////////////////////////////////
const double2 CurvMap::ref_vert[2][4] = {
    { { -1.0, -1.0 }, { 1.0, -1.0 }, { -1.0, 1.0 }, {  0.0, 0.0 } },
    { { -1.0, -1.0 }, { 1.0, -1.0 }, {  1.0, 1.0 }, { -1.0, 1.0 } }
  };

// subtraction of straight edge and nurbs curve
void CurvMap::nurbs_edge_0(Element* e, Nurbs* nurbs, int edge, double t, double& x, double& y, double& n_x, double& n_y, double& t_x, double& t_y)
{
  int va = edge;
  int vb = e->next_vert(edge);
  nurbs_edge(e, nurbs, edge, t, x, y, n_x, n_y, t_x, t_y);

  x -= 0.5 * ((1-t) * (e->vn[va]->x) + (1+t) * (e->vn[vb]->x));
  y -= 0.5 * ((1-t) * (e->vn[va]->y) + (1+t) * (e->vn[vb]->y));

  double k = 4.0 / ((1-t) * (1+t));
  x *= k;
  y *= k;
}

// calculation of nonpolynomial reference mapping on curved element
void CurvMap::calc_ref_map_tri(Element* e, Nurbs** nurbs, double xi_1, double xi_2, double& x, double& y)
{
  _F_
  double  fx,  fy;
  x = y = 0.0;

  for (unsigned int j = 0; j < e->nvert; j++)
  {
    int va = j;
    int vb = e->next_vert(j);
    double l_a = 0;
    double l_b = 0;
    switch(va) {
    case 0:
      l_a = lambda_0(xi_1, xi_2);
      break;
    case 1:
      l_a = lambda_1(xi_1, xi_2);
      break;
    case 2:
      l_a = lambda_2(xi_1, xi_2);
      break;
    }

    switch(vb) {
    case 0:
      l_b = lambda_0(xi_1, xi_2);
      break;
    case 1:
      l_b = lambda_1(xi_1, xi_2);
      break;
    case 2:
      l_b = lambda_2(xi_1, xi_2);
      break;
    }

    // vertex part
    x += e->vn[j]->x * l_a;
    y += e->vn[j]->y * l_a;

    if (!(((ref_vert[0][va][0] == xi_1) && (ref_vert[0][va][1] == xi_2)) ||
          ((ref_vert[0][vb][0] == xi_1) && (ref_vert[0][vb][1] == xi_2))))
    {
      // edge part
      double t = l_b - l_a;
      double n_x, n_y, t_x, t_y;
      nurbs_edge_0(e, nurbs[j], j, t, fx, fy, n_x, n_y, t_x, t_y);
      x += fx * l_a  * l_b;
      y += fy * l_a  * l_b;
    }
  }
}


void CurvMap::calc_ref_map_quad(Element* e, Nurbs** nurbs, double xi_1, double xi_2,
                              double& x, double& y)
{
  _F_
  double ex[4], ey[4];

  double n_x, n_y, t_x, t_y;
  nurbs_edge(e, nurbs[0], 0,  xi_1, ex[0], ey[0], n_x, n_y, t_x, t_y);
  nurbs_edge(e, nurbs[1], 1,  xi_2, ex[1], ey[1], n_x, n_y, t_x, t_y);
  nurbs_edge(e, nurbs[2], 2, -xi_1, ex[2], ey[2], n_x, n_y, t_x, t_y);
  nurbs_edge(e, nurbs[3], 3, -xi_2, ex[3], ey[3], n_x, n_y, t_x, t_y);

  x = (1-xi_2)/2.0 * ex[0] + (1+xi_1)/2.0 * ex[1] +
      (1+xi_2)/2.0 * ex[2] + (1-xi_1)/2.0 * ex[3] -
      (1-xi_1)*(1-xi_2)/4.0 * e->vn[0]->x - (1+xi_1)*(1-xi_2)/4.0 * e->vn[1]->x -
      (1+xi_1)*(1+xi_2)/4.0 * e->vn[2]->x - (1-xi_1)*(1+xi_2)/4.0 * e->vn[3]->x;

  y = (1-xi_2)/2.0 * ey[0] + (1+xi_1)/2.0 * ey[1] +
      (1+xi_2)/2.0 * ey[2] + (1-xi_1)/2.0 * ey[3] -
      (1-xi_1)*(1-xi_2)/4.0 * e->vn[0]->y - (1+xi_1)*(1-xi_2)/4.0 * e->vn[1]->y -
      (1+xi_1)*(1+xi_2)/4.0 * e->vn[2]->y - (1-xi_1)*(1+xi_2)/4.0 * e->vn[3]->y;
}


void CurvMap::calc_ref_map(Element* e, Nurbs** nurbs, double xi_1, double xi_2, double2& f)
{
  _F_
  if (e->get_mode() == HERMES_MODE_QUAD)
    calc_ref_map_quad(e, nurbs, xi_1, xi_2, f[0], f[1]);
  else
    calc_ref_map_tri(e, nurbs, xi_1, xi_2, f[0], f[1]);
}


//// projection based interpolation ////////////////////////////////////////////////////////////////

// preparation of projection matrices, Cholesky factorization
void CurvMap::precalculate_cholesky_projection_matrix_edge()
{
  _F_
  int order = ref_map_shapeset.get_max_order();
  int n = order - 1; // number of edge basis functions
  edge_proj_matrix = new_matrix<double>(n, n);

  // calculate projection matrix of maximum order
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      int o = i + j + 4;
      double2* pt = quad1d.get_points(o);
      double val = 0.0;
      for (int k = 0; k < quad1d.get_num_points(o); k++)
      {
        double fi = 0;
        double fj = 0;
        double x = pt[k][0];
        switch(i+2) {
        case 0:
          fi = lob0(x);
          break;
        case 1:
          fi = lob1(x);
          break;
        case 2:
          fi = lob2(x);
          break;
        case 3:
          fi = lob3(x);
          break;
        case 4:
          fi = lob4(x);
          break;
        case 5:
          fi = lob5(x);
          break;
        case 6:
          fi = lob6(x);
          break;
        case 7:
          fi = lob7(x);
          break;
        case 8:
          fi = lob8(x);
          break;
        case 9:
          fi = lob9(x);
          break;
        case 10:
          fi = lob10(x);
          break;
        case 11:
          fi = lob11(x);
          break;
        }
        switch(j+2) {
        case 0:
          fj = lob0(x);
          break;
        case 1:
          fj = lob1(x);
          break;
        case 2:
          fj = lob2(x);
          break;
        case 3:
          fj = lob3(x);
          break;
        case 4:
          fj = lob4(x);
          break;
        case 5:
          fj = lob5(x);
          break;
        case 6:
          fj = lob6(x);
          break;
        case 7:
          fj = lob7(x);
          break;
        case 8:
          fj = lob8(x);
          break;
        case 9:
          fj = lob9(x);
          break;
        case 10:
          fj = lob10(x);
          break;
        case 11:
          fj = lob11(x);
          break;
        }
        val += pt[k][1] * (fi * fj);
      }
      edge_proj_matrix[i][j] = edge_proj_matrix[j][i] = val;
    }
  }

  // Cholesky factorization of the matrix
  edge_p = new double[n];
  choldc(edge_proj_matrix, n, edge_p);
}

// calculate the H1 seminorm products (\phi_i, \phi_j) for all 0 <= i,j < n, n is the number of bubble functions
double** CurvMap::calculate_bubble_projection_matrix(int nb, int* indices)
{
  _F_
  double** mat = new_matrix<double>(nb, nb);

  for (int i = 0; i < nb; i++)
  {
    for (int j = i; j < nb; j++)
    {
      int ii = indices[i], ij = indices[j];
      int o = ref_map_shapeset.get_order(ii) + ref_map_shapeset.get_order(ij);
      o = std::max(H2D_GET_V_ORDER(o), H2D_GET_H_ORDER(o));

      ref_map_pss.set_active_shape(ii);
      ref_map_pss.set_quad_order(o);
      double* fni = ref_map_pss.get_fn_values();

      ref_map_pss.set_active_shape(ij);
      ref_map_pss.set_quad_order(o);
      double* fnj = ref_map_pss.get_fn_values();

      double3* pt = quad2d.get_points(o);
      double val = 0.0;
      for (int k = 0; k < quad2d.get_num_points(o); k++)
        val += pt[k][2] * (fni[k] * fnj[k]);

      mat[i][j] = mat[j][i] = val;
    }
  }

  return mat;
}


void CurvMap::precalculate_cholesky_projection_matrices_bubble()
{
  _F_
  // *** triangles ***
  ref_map_pss.set_mode(HERMES_MODE_TRIANGLE);
  int order = ref_map_shapeset.get_max_order();

  // calculate projection matrix of maximum order
  int nb = ref_map_shapeset.get_num_bubbles(order);
  int* indices = ref_map_shapeset.get_bubble_indices(order);
  bubble_proj_matrix_tri = calculate_bubble_projection_matrix(nb, indices);

  // cholesky factorization of the matrix
  bubble_tri_p = new double[nb];
  choldc(bubble_proj_matrix_tri, nb, bubble_tri_p);

  // *** quads ***
  ref_map_pss.set_mode(HERMES_MODE_QUAD);
  order = ref_map_shapeset.get_max_order();
  order = H2D_MAKE_QUAD_ORDER(order, order);

  // calculate projection matrix of maximum order
  nb = ref_map_shapeset.get_num_bubbles(order);
  indices = ref_map_shapeset.get_bubble_indices(order);
  bubble_proj_matrix_quad = calculate_bubble_projection_matrix(nb, indices);

  // cholesky factorization of the matrix
  bubble_quad_p = new double[nb];
  choldc(bubble_proj_matrix_quad, nb, bubble_quad_p);
}


//// edge part of projection based interpolation ///////////////////////////////////////////////////

// compute point (x,y) in reference element, edge vector (v1, v2)
void CurvMap::edge_coord(Element* e, int edge, double t, double2& x, double2& v)
{
  _F_
  int mode = e->get_mode();
  double2 a, b;
  a[0] = ctm.m[0] * ref_vert[mode][edge][0] + ctm.t[0];
  a[1] = ctm.m[1] * ref_vert[mode][edge][1] + ctm.t[1];
  b[0] = ctm.m[0] * ref_vert[mode][e->next_vert(edge)][0] + ctm.t[0];
  b[1] = ctm.m[1] * ref_vert[mode][e->next_vert(edge)][1] + ctm.t[1];

  for (int i = 0; i < 2; i++)
  {
    v[i] = b[i] - a[i];
    x[i] = a[i] + (t+1.0)/2.0 * v[i];
  }
  double lenght = sqrt(v[0] * v[0] + v[1] * v[1]);
  v[0] /= lenght; v[1] /= lenght;
}

void CurvMap::calc_edge_projection(Element* e, int edge, Nurbs** nurbs, int order, double2* proj)
{
  _F_
  ref_map_pss.set_active_element(e);

  int i, j, k;
  int mo1 = quad1d.get_max_order();
  int np = quad1d.get_num_points(mo1);
  int ne = order - 1;
  int mode = e->get_mode();

  assert(np <= 15 && ne <= 10);
  double2 fn[15];
  double rhside[2][10];
  memset(fn, 0, sizeof(double2) * np);
  memset(rhside[0], 0, sizeof(double) * ne);
  memset(rhside[1], 0, sizeof(double) * ne);

  double a_1, a_2, b_1, b_2;
  a_1 = ctm.m[0] * ref_vert[mode][edge][0] + ctm.t[0];
  a_2 = ctm.m[1] * ref_vert[mode][edge][1] + ctm.t[1];
  b_1 = ctm.m[0] * ref_vert[mode][e->next_vert(edge)][0] + ctm.t[0];
  b_2 = ctm.m[1] * ref_vert[mode][e->next_vert(edge)][1] + ctm.t[1];

  // values of nonpolynomial function in two vertices
  double2 fa, fb;
  calc_ref_map(e, nurbs, a_1, a_2, fa);
  calc_ref_map(e, nurbs, b_1, b_2, fb);

  double2* pt = quad1d.get_points(mo1);
  for (j = 0; j < np; j++) // over all integration points
  {
    double2 x, v;
    double t = pt[j][0];
    edge_coord(e, edge, t, x, v);
    calc_ref_map(e, nurbs, x[0], x[1], fn[j]);

    for (k = 0; k < 2; k++)
      fn[j][k] = fn[j][k] - (fa[k] + (t+1)/2.0 * (fb[k] - fa[k]));
  }

  double2* result = proj + e->nvert + edge * (order - 1);
  for (k = 0; k < 2; k++)
  {
    for (i = 0; i < ne; i++)
    {
      for (j = 0; j < np; j++)
      {
        double t = pt[j][0];
        double fi = 0;
        switch(i+2) {
        case 0:
          fi = lob0(t);
          break;
        case 1:
          fi = lob1(t);
          break;
        case 2:
          fi = lob2(t);
          break;
        case 3:
          fi = lob3(t);
          break;
        case 4:
          fi = lob4(t);
          break;
        case 5:
          fi = lob5(t);
          break;
        case 6:
          fi = lob6(t);
          break;
        case 7:
          fi = lob7(t);
          break;
        case 8:
          fi = lob8(t);
          break;
        case 9:
          fi = lob9(t);
          break;
        case 10:
          fi = lob10(t);
          break;
        case 11:
          fi = lob11(t);
          break;
        }
        rhside[k][i] += pt[j][1] * (fi * fn[j][k]);
      }
    }
    // solve
    cholsl(edge_proj_matrix, ne, edge_p, rhside[k], rhside[k]);
    for (i = 0; i < ne; i++)
      result[i][k] = rhside[k][i];
  }
}

//// bubble part of projection based interpolation /////////////////////////////////////////////////

void CurvMap::old_projection(Element* e, int order, double2* proj, double* old[2])
{
  _F_
  int mo2 = quad2d.get_max_order();
  int np = quad2d.get_num_points(mo2);

  for (unsigned int k = 0; k < e->nvert; k++) // loop over vertices
  {
    // vertex basis functions in all integration points
    double* vd;
    int index_v = ref_map_shapeset.get_vertex_index(k);
    ref_map_pss.set_active_shape(index_v);
    ref_map_pss.set_quad_order(mo2);
    vd = ref_map_pss.get_fn_values();

    for (int m = 0; m < 2; m++)   // part 0 or 1
      for (int j = 0; j < np; j++)
        old[m][j] += proj[k][m] * vd[j];

    for (int ii = 0; ii < order - 1; ii++)
    {
      // edge basis functions in all integration points
      double* ed;
      int index_e = ref_map_shapeset.get_edge_index(k,0,ii+2);
      ref_map_pss.set_active_shape(index_e);
      ref_map_pss.set_quad_order(mo2);
      ed = ref_map_pss.get_fn_values();

      for (int m = 0; m < 2; m++)  //part 0 or 1
        for (int j = 0; j < np; j++)
          old[m][j] += proj[e->nvert + k * (order-1) + ii][m] * ed[j];
    }
  }
}

void CurvMap::calc_bubble_projection(Element* e, Nurbs** nurbs, int order, double2* proj)
{
  _F_
  ref_map_pss.set_active_element(e);

  int i, j, k;
  int mo2 = quad2d.get_max_order();
  int np = quad2d.get_num_points(mo2);
  int qo = e->is_quad() ? H2D_MAKE_QUAD_ORDER(order, order) : order;
  int nb = ref_map_shapeset.get_num_bubbles(qo);

  double2* fn = new double2[np];
  memset(fn, 0, np * sizeof(double2));

  double* rhside[2];
  double* old[2];
  for (i = 0; i < 2; i++) {
    rhside[i] = new double[nb];
    old[i] = new double[np];
    memset(rhside[i], 0, sizeof(double) * nb);
    memset(old[i], 0, sizeof(double) * np);
  }

  // compute known part of projection (vertex and edge part)
  old_projection(e, order, proj, old);

  // fn values of both components of nonpolynomial function
  double3* pt = quad2d.get_points(mo2);
  for (j = 0; j < np; j++)  // over all integration points
  {
    double2 a;
    a[0] = ctm.m[0] * pt[j][0] + ctm.t[0];
    a[1] = ctm.m[1] * pt[j][1] + ctm.t[1];
    calc_ref_map(e, nurbs, a[0], a[1], fn[j]);
  }

  double2* result = proj + e->nvert + e->nvert * (order - 1);
  for (k = 0; k < 2; k++)
  {
    for (i = 0; i < nb; i++) // loop over bubble basis functions
    {
      // bubble basis functions in all integration points
      double *bfn;
      int index_i = ref_map_shapeset.get_bubble_indices(qo)[i];
      ref_map_pss.set_active_shape(index_i);
      ref_map_pss.set_quad_order(mo2);
      bfn = ref_map_pss.get_fn_values();

      for (j = 0; j < np; j++) // over all integration points
        rhside[k][i] += pt[j][2] * (bfn[j] * (fn[j][k] - old[k][j]));
    }

    // solve
    if (e->nvert == 3)
      cholsl(bubble_proj_matrix_tri, nb, bubble_tri_p, rhside[k], rhside[k]);
    else
      cholsl(bubble_proj_matrix_quad, nb, bubble_quad_p, rhside[k], rhside[k]);

    for (i = 0; i < nb; i++)
      result[i][k] = rhside[k][i];
  }

  for (i = 0; i < 2; i++) {
    delete [] rhside[i];
    delete [] old[i];
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void CurvMap::ref_map_projection(Element* e, Nurbs** nurbs, int order, double2* proj)
{
  _F_
  // vertex part
  for (unsigned int i = 0; i < e->nvert; i++)
  {
    proj[i][0] = e->vn[i]->x;
    proj[i][1] = e->vn[i]->y;
  }

  if (e->cm->toplevel == false)
    e = e->cm->parent;

  // edge part
  for (int edge = 0; edge < (int)e->nvert; edge++)
    calc_edge_projection(e, edge, nurbs, order, proj);

  //bubble part
  calc_bubble_projection(e, nurbs, order, proj);
}


void CurvMap::update_refmap_coeffs(Element* e)
{
  _F_
  ref_map_pss.set_quad_2d(&quad2d);
  //ref_map_pss.set_active_element(e);

  // calculation of projection matrices
  if (edge_proj_matrix == NULL) precalculate_cholesky_projection_matrix_edge();
  if (bubble_proj_matrix_tri == NULL) precalculate_cholesky_projection_matrices_bubble();

  ref_map_pss.set_mode(e->get_mode());
  ref_map_shapeset.set_mode(e->get_mode());

  // allocate projection coefficients
  int nv = e->nvert;
  int ne = order - 1;
  int qo = e->is_quad() ? H2D_MAKE_QUAD_ORDER(order, order) : order;
  int nb = ref_map_shapeset.get_num_bubbles(qo);
  nc = nv + nv*ne + nb;
  if (coeffs != NULL) {
    delete [] coeffs;
    coeffs = NULL;
  }
  coeffs = new double2[nc];

  // WARNING: do not change the format of the array 'coeffs'. If it changes,
  // RefMap::set_active_element() has to be changed too.

  Nurbs** nurbs;
  if (toplevel == false)
  {
    ref_map_pss.set_active_element(e);
    ref_map_pss.set_transform(part);
    nurbs = parent->cm->nurbs;
  }
  else
  {
    ref_map_pss.reset_transform();
    nurbs = e->cm->nurbs;
  }
  ctm = *(ref_map_pss.get_ctm());
  ref_map_pss.reset_transform(); // fixme - do we need this?

  // calculation of new projection coefficients
  ref_map_projection(e, nurbs, order, coeffs);
}

void CurvMap::get_mid_edge_points(Element* e, double2* pt, int n)
{
  _F_
  Nurbs** nurbs = this->nurbs;
  Transformable tran;
  tran.set_active_element(e);

  if (toplevel == false)
  {
    tran.set_transform(part);
    e = e->cm->parent;
    nurbs = e->cm->nurbs;
  }

  ctm = *(tran.get_ctm());
  double xi_1, xi_2;
  for (int i = 0; i < n; i++)
  {
    xi_1 = ctm.m[0] * pt[i][0] + ctm.t[0];
    xi_2 = ctm.m[1] * pt[i][1] + ctm.t[1];
    calc_ref_map(e, nurbs, xi_1, xi_2, pt[i]);
  }
}

void Nurbs::unref()
{
  _F_
  if (!--ref) // fixme: possible leak, we need ~Nurbs too
  {
    delete [] pt;
    delete [] kv;
    delete this;
  }
}


CurvMap::CurvMap(CurvMap* cm)
{
  _F_
  memcpy(this, cm, sizeof(CurvMap));
  coeffs = new double2[nc];
  memcpy(coeffs, cm->coeffs, sizeof(double2) * nc);

  if (toplevel)
    for (int i = 0; i < 4; i++)
      if (nurbs[i] != NULL)
        nurbs[i]->ref++;
}

CurvMap::~CurvMap()
{
  _F_
  if (coeffs != NULL) {
    delete [] coeffs;
    coeffs = NULL;
  }
  if (toplevel)
    for (int i = 0; i < 4; i++)
      if (nurbs[i] != NULL)
        nurbs[i]->unref();
}
