#include "common.h"
#include "shapeset_h1_all.h"
#include "shapeset_common.h"
#include "precalc.h"
#include "curved.h"
#include "mesh.h"
#include "quad_all.h"
#include "matrix_old.h"


// defined in refmap.cpp
extern H1ShapesetOrtho ref_map_shapeset;
extern PrecalcShapeset ref_map_pss;

static double** edge_proj_matrix = NULL;  //projection matrix for each edge is the same
static double** bubble_proj_matrix_tri = NULL; //projection matrix for bubble
static double** bubble_proj_matrix_quad = NULL; //projection matrix for bubble

static double* edge_p = NULL;  // diagonal vector in choleski factorization
static double* bubble_tri_p = NULL; // diagonal vector in choleski factorization
static double* bubble_quad_p = NULL; // diagonal vector in choleski factorization

static Quad1DStd quad1d;
static Quad2DStd quad2d;


//// NURBS //////////////////////////////////////////////////////////////////////////////////////////

// recursive calculation of the basis function N_i,k and their derivatives
double nurbs_basis_fn(int i, int k, double t, double* knot, double& fn, double& d_fn)
{
  if (k == 0)
  {
    fn = (t >= knot[i] && t <= knot[i+1] && knot[i] < knot[i+1]) ? 1.0 : 0.0;
    d_fn = 0.0;
  }
  else
  {
    double N1, N2, dN1, dN2;
    nurbs_basis_fn(i, k-1, t, knot, N1, dN1);
    nurbs_basis_fn(i+1, k-1, t, knot, N2, dN2);

    fn = d_fn = 0.0;
    if (knot[i+k] != knot[i])
    {
      fn += ((t - knot[i]) / (knot[i+k] - knot[i])) * N1;
      d_fn += N1 / (knot[i+k] - knot[i]) + ((t - knot[i])/(knot[i+k] - knot[i])) * dN1;
    }
    if (knot[i+k+1] != knot[i+1])
    {
      fn += ((knot[i+k+1] - t) / (knot[i+k+1] - knot[i+1])) * N2;
      d_fn += -N2 / (knot[i+k+1] - knot[i+1]) + ((knot[i+k+1] - t)/(knot[i+k+1] - knot[i+1])) * dN2;
    }
  }
}


// nurbs curve: t goes from -1 to 1, function returns coordinates x, y in plane and derivatives
void nurbs_edge(Element* e, Nurbs* nurbs, int edge, double t, double& x, double& y, double& d_x, double& d_y)
{
  if (nurbs == NULL)
  {
    double2 v;
    v[0] = e->vn[e->next_vert(edge)]->x - e->vn[edge]->x;
    v[1] = e->vn[e->next_vert(edge)]->y - e->vn[edge]->y;
    x = e->vn[edge]->x + (t+1)/2.0 * v[0];
    y = e->vn[edge]->y + (t+1)/2.0 * v[1];
    d_x = 0.5 * v[0];
    d_y = 0.5 * v[1];
  }
  else
  {
    double3* cp = nurbs->pt;

    int i;
    double basis;
    double d_basis;

    double sum = 0.0;   // sum of basis fn and weights
    double d_sum = 0.0; // sum of derivative of basis fn and weights

    x = y = 0.0;
    d_x = d_y = 0.0;
    t = (t + 1) / 2.0; // nurbs curves are parametrized from 0 to 1
    for (i = 0; i < nurbs->np; i++)
    {
      nurbs_basis_fn(i, nurbs->degree, t, nurbs->kv, basis, d_basis);
      sum    += cp[i][2] * basis;
      d_sum  += cp[i][2] * d_basis;

      x   += cp[i][2] * basis * cp[i][0];
      y   += cp[i][2] * basis * cp[i][1];
      d_x += cp[i][2] * d_basis * cp[i][0];
      d_y += cp[i][2] * d_basis * cp[i][1];
    }

    // coordinates x and y
    // derivative of components x and y - derivative of fraction
    d_x  = d_x * sum - x * d_sum;
    d_x /= sum * sum;
    d_y  = d_y * sum - y * d_sum;
    d_y /= sum * sum;

    x /= sum;
    y /= sum;

  }
}


//// non-polynomial reference map //////////////////////////////////////////////////////////////////////////////////

// definition of vertex basis functions for triangle
double lambda_0(double x, double y) { return -0.5 * (x + y); }
double lambda_1(double x, double y) { return  0.5 * (x + 1); }
double lambda_2(double x, double y) { return  0.5 * (y + 1); }

double (*lambda[3])(double, double) = { lambda_0, lambda_1, lambda_2 };
double lambda_dx[3] = { -0.5, 0.5, 0.0 };
double lambda_dy[3] = { -0.5, 0.0, 0.5 };

// 1D Lobatto functions
double lob0(double x)  { return l0(x); }
double lob1(double x)  { return l1(x); }
double lob2(double x)  { return l2(x); }
double lob3(double x)  { return l3(x); }
double lob4(double x)  { return l4(x); }
double lob5(double x)  { return l5(x); }
double lob6(double x)  { return l6(x); }
double lob7(double x)  { return l7(x); }
double lob8(double x)  { return l8(x); }
double lob9(double x)  { return l9(x); }
double lob10(double x) { return l10(x); }
double lob11(double x) { return l11(x); }

// derivatives of Lobatto functions
double dlob0(double x)  { return dl0(x); }
double dlob1(double x)  { return dl1(x); }
double dlob2(double x)  { return dl2(x); }
double dlob3(double x)  { return dl3(x); }
double dlob4(double x)  { return dl4(x); }
double dlob5(double x)  { return dl5(x); }
double dlob6(double x)  { return dl6(x); }
double dlob7(double x)  { return dl7(x); }
double dlob8(double x)  { return dl8(x); }
double dlob9(double x)  { return dl9(x); }
double dlob10(double x) { return dl10(x); }
double dlob11(double x) { return dl11(x); }

double     (*lob[12])(double) = {  lob0,  lob1,  lob2,  lob3,  lob4,  lob5,  lob6,  lob7,  lob8,  lob9,  lob10,  lob11 };
double (*der_lob[12])(double) = { dlob0, dlob1, dlob2, dlob3, dlob4, dlob5, dlob6, dlob7, dlob8, dlob9, dlob10, dlob11 };

static double2 ref_vert[2][4] =
{
  { { -1.0, -1.0 }, { 1.0, -1.0 }, { -1.0, 1.0 }, {  0.0, 0.0 } },
  { { -1.0, -1.0 }, { 1.0, -1.0 }, {  1.0, 1.0 }, { -1.0, 1.0 } }
};

// subtraction of straight edge and nurbs curve
void nurbs_edge_0(Element* e, Nurbs* nurbs, int edge, double t, double& x, double& y, double& dx, double& dy)
{
  int va = edge;
  int vb = e->next_vert(edge);
  double k = 1.0/4.0 * (1-t) * (1+t);

  nurbs_edge(e, nurbs, edge, t, x, y, dx, dy);

  x = x - 0.5 * ((1-t) * (e->vn[va]->x) + (1+t) * (e->vn[vb]->x));
  y = y - 0.5 * ((1-t) * (e->vn[va]->y) + (1+t) * (e->vn[vb]->y));

  dx += 0.5 * (e->vn[va]->x - e->vn[vb]->x);
  dy += 0.5 * (e->vn[va]->y - e->vn[vb]->y);

  dx = (dx * k - x * (-0.5) * t) / (k * k); // derivative of fraction
  dy = (dy * k - y * (-0.5) * t) / (k * k);

  x = x / k;
  y = y / k;

}


// calculation of nonpolynomial reference mapping on curved element
void calc_ref_map_tri(Element* e, Nurbs** nurbs, double xi_1, double xi_2,
                      double& x, double& y, double& x_dx, double& x_dy, double& y_dx, double& y_dy)
{
  int i, j;
  double  fx,  fy;
  double dfx, dfy;

  x = y = 0.0;
  x_dx = x_dy = 0.0;
  y_dx = y_dy = 0.0;

  for (j = 0; j < e->nvert; j++)
  {
    int va = j;
    int vb = e->next_vert(j);
    double l_a  = lambda[va](xi_1, xi_2);
    double l_b  = lambda[vb](xi_1, xi_2);
    double lx_a = lambda_dx[va], lx_b = lambda_dx[vb];
    double ly_a = lambda_dy[va], ly_b = lambda_dy[vb];

    // vertex part
    x    += e->vn[j]->x * l_a;
    y    += e->vn[j]->y * l_a;
    x_dx += e->vn[j]->x * lx_a;
    x_dy += e->vn[j]->x * ly_a;
    y_dx += e->vn[j]->y * lx_a;
    y_dy += e->vn[j]->y * ly_a;

    if (!(((ref_vert[0][va][0] == xi_1) && (ref_vert[0][va][1] == xi_2)) ||
          ((ref_vert[0][vb][0] == xi_1) && (ref_vert[0][vb][1] == xi_2))))
    {
      // edge part
      double t = l_b - l_a;
      nurbs_edge_0(e, nurbs[j], j, t, fx, fy, dfx, dfy);
      x    += fx * l_a  * l_b;
      y    += fy * l_a  * l_b;

      x_dx += fx * lx_a * l_b  +  fx * l_a * lx_b  +  dfx * l_a * l_b * (lx_b - lx_a);
      x_dy += fx * ly_a * l_b  +  fx * l_a * ly_b  +  dfx * l_a * l_b * (ly_b - ly_a);
      y_dx += fy * lx_a * l_b  +  fy * l_a * lx_b  +  dfy * l_a * l_b * (lx_b - lx_a);
      y_dy += fy * ly_a * l_b  +  fy * l_a * ly_b  +  dfy * l_a * l_b * (ly_b - ly_a);
    }
  }
}


void calc_ref_map_quad(Element* e, Nurbs** nurbs, double xi_1, double xi_2,
                       double& x, double& y, double& x_dx, double& x_dy, double& y_dx, double& y_dy)
{
  int i, j;
  double ex[4], ey[4];
  double dex[4], dey[4];

  nurbs_edge(e, nurbs[0], 0,  xi_1, ex[0], ey[0], dex[0], dey[0]);
  nurbs_edge(e, nurbs[1], 1,  xi_2, ex[1], ey[1], dex[1], dey[1]);
  nurbs_edge(e, nurbs[2], 2, -xi_1, ex[2], ey[2], dex[2], dey[2]); dex[2] = -dex[2]; dey[2] = -dey[2];
  nurbs_edge(e, nurbs[3], 3, -xi_2, ex[3], ey[3], dex[3], dey[3]); dex[3] = -dex[3]; dey[3] = -dey[3];

  x = (1-xi_2)/2.0 * ex[0] + (1+xi_1)/2.0 * ex[1] +
      (1+xi_2)/2.0 * ex[2] + (1-xi_1)/2.0 * ex[3] -
      (1-xi_1)*(1-xi_2)/4.0 * e->vn[0]->x - (1+xi_1)*(1-xi_2)/4.0 * e->vn[1]->x -
      (1+xi_1)*(1+xi_2)/4.0 * e->vn[2]->x - (1-xi_1)*(1+xi_2)/4.0 * e->vn[3]->x;
  y = (1-xi_2)/2.0 * ey[0] + (1+xi_1)/2.0 * ey[1] +
      (1+xi_2)/2.0 * ey[2] + (1-xi_1)/2.0 * ey[3] -
      (1-xi_1)*(1-xi_2)/4.0 * e->vn[0]->y - (1+xi_1)*(1-xi_2)/4.0 * e->vn[1]->y -
      (1+xi_1)*(1+xi_2)/4.0 * e->vn[2]->y - (1-xi_1)*(1+xi_2)/4.0 * e->vn[3]->y;

  x_dx = (1-xi_2)/2.0 * dex[0] + (0.5) * ex[1] + (1+xi_2)/2.0 * dex[2] + (-0.5) * ex[3] +
         (1-xi_2)/4.0 * e->vn[0]->x - (1-xi_2)/4.0 * e->vn[1]->x - (1+xi_2)/4.0 * e->vn[2]->x + (1+xi_2)/4.0 * e->vn[3]->x;
  x_dy = -0.5 * ex[0] + (1+xi_1)/2.0 * dex[1] + (0.5) * ex[2] + (1-xi_1)/2.0 * dex[3] +
         (1-xi_1)/4.0 * e->vn[0]->x + (1+xi_1)/4.0 * e->vn[1]->x - (1+xi_1)/4.0 * e->vn[2]->x - (1-xi_1)/4.0 * e->vn[3]->x;
  y_dx = (1-xi_2)/2.0 * dey[0] + (0.5) * ey[1] + (1+xi_2)/2.0 * dey[2] + (-0.5) * ey[3] +
         (1-xi_2)/4.0 * e->vn[0]->y - (1-xi_2)/4.0 * e->vn[1]->y - (1+xi_2)/4.0 * e->vn[2]->y + (1+xi_2)/4.0 * e->vn[3]->y;
  y_dy = -0.5 * ey[0] + (1+xi_1)/2.0 * dey[1] + (0.5) * ey[2] + (1-xi_1)/2.0 * dey[3] +
         (1-xi_1)/4.0 * e->vn[0]->y + (1+xi_1)/4.0 * e->vn[1]->y - (1+xi_1)/4.0 * e->vn[2]->y - (1-xi_1)/4.0 * e->vn[3]->y;

}


void calc_ref_map(Element* e, Nurbs** nurbs, double xi_1, double xi_2,
                  double2& f, double2& f_dx, double2& f_dy)
{
  if (e->get_mode() == MODE_QUAD)
    calc_ref_map_quad(e, nurbs, xi_1, xi_2, f[0], f[1], f_dx[0], f_dy[0], f_dx[1], f_dy[1]);
  else
    calc_ref_map_tri( e, nurbs, xi_1, xi_2, f[0], f[1], f_dx[0], f_dy[0], f_dx[1], f_dy[1]);
}


//// projection based interpolation ////////////////////////////////////////////////////////////////

// preparation of projection matrices, Cholesky factorization

// calculate the H1 products over edges (\phi_i, \phi_j) for all 0 <= i,j < n, n is number of edge function
/*double** calculate_edge_projection_matrix_h1(int order)
{
  int n = order - 1; // number of edge basis functions
  double** mat = new_matrix<double>(n, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      int o = i + j + 4;
      double2* pt = quad1d.get_points(o);
      double val = 0.0;
      for (int k = 0; k < quad1d.get_num_points(o); k++)
      {
        double x = pt[k][0];
        double fi = lob[i+2](x);
        double fj = lob[j+2](x);
        double fi_dx = der_lob[i+2](x);
        double fj_dx = der_lob[j+2](x);
        val += pt[k][1] * (fi * fj   +   fi_dx * fj_dx);
      }
      mat[i][j] = mat[j][i] = val;
    }
  }

  return mat;
}*/


void precalculate_cholesky_projection_matrix_edge()
{
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
        double x = pt[k][0];
        double fi = lob[i+2](x);
        double fj = lob[j+2](x);
        double fi_dx = der_lob[i+2](x);
        double fj_dx = der_lob[j+2](x);
        //val += pt[k][1] * (fi * fj   +   fi_dx * fj_dx);
        val += pt[k][1] * (fi * fj);
      }
      edge_proj_matrix[i][j] = edge_proj_matrix[j][i] = val;
    }
  }

  // Cholesky factorization of matrix
  edge_p = new double[n];
  choldc(edge_proj_matrix, n, edge_p);
}


// calculate the H1 seminorm products (\phi_i, \phi_j) for all 0 <= i,j < n, n is number of bubble functions
double** calculate_bubble_projection_matrix_h1(int order)
{
  int n = ref_map_shapeset.get_num_bubbles(order); // number of bubble basis functions
  double** mat = new_matrix<double>(n, n);

  double *i_dx, *i_dy, *j_dx, *j_dy;
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      int index_i = ref_map_shapeset.get_bubble_indices(order)[i];
      int index_j = ref_map_shapeset.get_bubble_indices(order)[j];
      int o = ref_map_shapeset.get_order(index_i) + ref_map_shapeset.get_order(index_j);

      ref_map_pss.set_active_shape(index_i);
      ref_map_pss.set_quad_order(o);
      ref_map_pss.get_dx_dy_values(i_dx, i_dy);
      //i_dx = ref_map_pss.get_fn_values();
      ref_map_pss.set_active_shape(index_j);
      ref_map_pss.set_quad_order(o);
      ref_map_pss.get_dx_dy_values(j_dx, j_dy);
      //j_dx = ref_map_pss.get_fn_values();

      double3* pt = quad2d.get_points(o);
      double val = 0.0;
      for (int k = 0; k < quad2d.get_num_points(o); k++)
      {
        val += pt[k][2] * (i_dx[k] * j_dx[k]  +  i_dy[k] * j_dy[k]);
        //val += pt[k][2] * (i_dx[k] * j_dx[k]);
      }
      mat[i][j] = mat[j][i] = val;
    }
  }

  return mat;
}


void precalculate_cholesky_projection_matrices_bubble()
{
  // triangles
  ref_map_pss.set_mode(MODE_TRIANGLE);
  int order = ref_map_shapeset.get_max_order();
  int n = ref_map_shapeset.get_num_bubbles(order);

  // calculate projection matrix of maximum order
  bubble_proj_matrix_tri = calculate_bubble_projection_matrix_h1(order);

  // Cholesky factorization of matrix
  bubble_tri_p = new double[n];
  choldc(bubble_proj_matrix_tri, n, bubble_tri_p);

  // quads
  ref_map_pss.set_mode(MODE_QUAD);
  order = ref_map_shapeset.get_max_order();
  n = ref_map_shapeset.get_num_bubbles(make_quad_order(order, order));

  // calculate projection matrix of maximum order
  bubble_proj_matrix_quad = calculate_bubble_projection_matrix_h1(make_quad_order(order, order));

  // Cholesky factorization of the matrix
  bubble_quad_p = new double[n];
  choldc(bubble_proj_matrix_quad, n, bubble_quad_p);
}


//// edge part of projection based interpolation ///////////////////////////////////////////////////

// compute point (x,y) in reference element, edge vector (v1, v2)
void edge_coord(Element* e, int edge, double t, double2& x, double2& v)
{
  int mode = e->get_mode();
  for (int i = 0; i < 2; i++)
  {
    v[i] = ref_vert[mode][e->next_vert(edge)][i] - ref_vert[mode][edge][i];
    x[i] = ref_vert[mode][edge][i] + (t+1.0)/2.0 * v[i];
  }
  double lenght = sqrt(v[0] * v[0] + v[1] * v[1]);
  v[0] /= lenght; v[1] /= lenght;
}


void calc_edge_projection(Element* e, int edge, Nurbs** nurbs, int order, double2* proj)
{
  double2* rhs = proj + e->nvert + edge * (order - 1);
  int mo1 = quad1d.get_max_order();
  int np = quad1d.get_num_points(mo1);
  int ne = order - 1;
  int mode = e->get_mode();

  double2* fn = new double2[np];
  double2* fn_d = new double2[np];
  memset(fn, 0, sizeof(double2) * np);
  memset(fn_d, 0, sizeof(double2) * np);

  double* rhside[2];
  for (int i = 0; i < 2; i++) {
    rhside[i] = new double[ne];
    memset(rhside[i], 0, sizeof(double) * ne);
  }

/*  double* old[2][2];
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
      old[i][j] = new double[np];
      memset(old[i][j], 0, sizeof(double) * np);
    }
*/
  // values of nonpolynomial function in two vertices
  double2 fa, fa_x, fa_y, fb, fb_x, fb_y;
  calc_ref_map(e, nurbs, ref_vert[mode][edge][0], ref_vert[mode][edge][1], fa, fa_x, fa_y);
  calc_ref_map(e, nurbs, ref_vert[mode][e->next_vert(edge)][0], ref_vert[mode][e->next_vert(edge)][1], fb, fb_x, fb_y);

  double2* pt = quad1d.get_points(mo1);
  for (int j = 0; j < np; j++)  // over all integration points
  {

    double2 x, v, fn_x, fn_y;
    double t = pt[j][0];
    edge_coord(e, edge, t, x, v);

    calc_ref_map(e, nurbs, x[0], x[1], fn[j], fn_x, fn_y);

    for (int part = 0; part < 2; part++)
    {
      fn[j][part] = fn[j][part] - (fa[part] + (t+1)/2.0 * (fb[part] - fa[part]));
      fn_d[j][part] = (fn_x[part] * v[0] + fn_y[part] * v[1]) - (0.5 * (fb[part] - fa[part]));
    }
  }

  for (int part = 0; part < 2; part++)
  {
    for (int i = 0; i < ne; i++)
    {
      for (int j = 0; j < np; j++)
      {
        double t = pt[j][0];
        double fi = lob[i+2](t);
        double fi_d = der_lob[i+2](t);

        rhside[part][i] += pt[j][1] * (fi * fn[j][part]  +  fi_d * fn_d[j][part]);
        //rhside[part][i] += pt[j][1] * (fi_d * fn_d[j][part]);
        //rhside[part][i] += pt[j][1] * (fi * fn[j][part]);
      }
    }

    // solve
    cholsl(edge_proj_matrix, ne, edge_p, rhside[part], rhside[part]);
    for (int i = 0; i < ne; i++)
      rhs[i][part] = rhside[part][i];
  }
}

void calc_edge_projection2(Element* e, int edge, Nurbs** nurbs, int order, double2* proj)
{
  double2* rhs = proj + e->nvert + edge * (order - 1);
  int mo1 = quad1d.get_max_order();
  int np = quad1d.get_num_points(mo1);
  int ne = order - 1;
  int mode = e->get_mode();

  double2* fn = new double2[np];
  double2* fn_d = new double2[np];
  memset(fn, 0, sizeof(double2) * np);
  memset(fn_d, 0, sizeof(double2) * np);

  double* rhside[2];
  for (int i = 0; i < 2; i++) {
    rhside[i] = new double[ne];
    memset(rhside[i], 0, sizeof(double) * ne);
  }

/*  double* old[2][2];
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
      old[i][j] = new double[np];
      memset(old[i][j], 0, sizeof(double) * np);
    }
*/
  // values of nonpolynomial function in two vertices
  double2 fa, fa_x, fa_y, fb, fb_x, fb_y;
  calc_ref_map(e, nurbs, ref_vert[mode][edge][0], ref_vert[mode][edge][1], fa, fa_x, fa_y);
  calc_ref_map(e, nurbs, ref_vert[mode][e->next_vert(edge)][0], ref_vert[mode][e->next_vert(edge)][1], fb, fb_x, fb_y);

  double2* pt = quad1d.get_points(mo1);
  for (int j = 0; j < np; j++)  // over all integration points
  {

    double2 x, v, fn_x, fn_y;
    double t = pt[j][0];
    edge_coord(e, edge, t, x, v);

    calc_ref_map(e, nurbs, x[0], x[1], fn[j], fn_x, fn_y);

    for (int part = 0; part < 2; part++)
    {
      fn[j][part] = fn[j][part] - (fa[part] + (t+1)/2.0 * (fb[part] - fa[part]));
      fn_d[j][part] = (fn_x[part] * v[0] + fn_y[part] * v[1]) - (0.5 * (fb[part] - fa[part]));
    }
  }

  for (int part = 0; part < 2; part++)
  {
    for (int i = 0; i < ne; i++)
    {
      for (int j = 0; j < np; j++)
      {
        double t = pt[j][0];
        double fi = lob[i+2](t);
        double fi_d = der_lob[i+2](t);

        //rhside[part][i] += pt[j][1] * (fi * fn[j][part]  +  fi_d * fn_d[j][part]);
        //rhside[part][i] += pt[j][1] * (fi_d * fn_d[j][part]);
        rhside[part][i] += pt[j][1] * (fi * fn[j][part]);
      }
    }

    // solve
    cholsl(edge_proj_matrix, ne, edge_p, rhside[part], rhside[part]);
    for (int i = 0; i < ne; i++)
      rhs[i][part] = rhside[part][i];
  }
}


//// bubble part of projection based interpolation /////////////////////////////////////////////////

void old_projection(Element* e, int order, double2* proj, double* old[2][2])
{

  int mo2 = quad2d.get_max_order();
  int np = quad2d.get_num_points(mo2);

  for (int k = 0; k < e->nvert; k++) // loop over vertices
  {
    // derivatives of vertex basis functions in all integration points
    double* vd[2];
    int index_v = ref_map_shapeset.get_vertex_index(k);
    ref_map_pss.set_active_shape(index_v);
    ref_map_pss.set_quad_order(mo2);
    ref_map_pss.get_dx_dy_values(vd[0], vd[1]);
    //vd[0] = ref_map_pss.get_fn_values();

    for (int m = 0; m < 2; m++)   // part 0 or 1
      for (int l = 0; l < 2; l++) // derivative - x or y
        for (int j = 0; j < np; j++)
          old[m][l][j] += proj[k][m] * vd[l][j];

    for (int ii = 0; ii < order - 1; ii++)
    {
      // derivatives of edge basis functions in all integration points
      double* ed[2];
      int index_e = ref_map_shapeset.get_edge_index(k,0,ii+2);
      ref_map_pss.set_active_shape(index_e);
      ref_map_pss.set_quad_order(mo2);
      ref_map_pss.get_dx_dy_values(ed[0], ed[1]);
      //ed[0] = ref_map_pss.get_fn_values();

      for (int m = 0; m < 2; m++)  //part 0 or 1
        for (int l = 0; l < 2; l++)  //derivative  x or y
          for (int j = 0; j < np; j++)
            old[m][l][j] += proj[e->nvert + k * (order-1) + ii][m] * ed[l][j];
    }
  }
}

void old_projection2(Element* e, int order, double2* proj, double* old[2])
{

  int mo2 = quad2d.get_max_order();
  int np = quad2d.get_num_points(mo2);

  for (int k = 0; k < e->nvert; k++) // loop over vertices
  {
    // derivatives of vertex basis functions in all integration points
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
      // derivatives of edge basis functions in all integration points
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


void calc_bubble_projection(Element* e, Nurbs** nurbs, int order, double2* proj)
{
  ref_map_pss.set_active_element(e);

  double2* rhs = proj + e->nvert + e->nvert * (order - 1);
  int mo2 = quad2d.get_max_order();
  int np = quad2d.get_num_points(mo2);
  int qo = e->is_quad() ? make_quad_order(order, order) : order;
  int nb = ref_map_shapeset.get_num_bubbles(qo);

  double2 fn[np];
  double2 fn_dx[np];
  double2 fn_dy[np];
  memset(fn, 0, sizeof(fn));
  memset(fn_dx, 0, sizeof(fn_dx));
  memset(fn_dy, 0, sizeof(fn_dy));

  double* rhside[2];
  for (int i = 0; i < 2; i++) {
    rhside[i] = new double[nb];
    memset(rhside[i], 0, sizeof(double) * nb);
  }

  double* old[2][2];
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
      old[i][j] = new double[np];
      memset(old[i][j], 0, sizeof(double) * np);
    }

  // compute known part of projection (vertex and edge part)
  old_projection(e, order, proj, old);

  // fn values and derivatives of both components of nonpolynomial function
  double3* pt = quad2d.get_points(mo2);
  for (int j = 0; j < np; j++)  // over all integration points
    calc_ref_map(e, nurbs, pt[j][0], pt[j][1], fn[j], fn_dx[j], fn_dy[j]);

  for (int part = 0; part < 2; part++)
  {
    for (int i = 0; i < nb; i++) // loop over bubble basis functions
    {
      // derivatives of bubble basis functions in all integration points
      double *dx, *dy;
      int index_i = ref_map_shapeset.get_bubble_indices(qo)[i];
      ref_map_pss.set_active_shape(index_i);
      ref_map_pss.set_quad_order(mo2);
      ref_map_pss.get_dx_dy_values(dx, dy);
      //dx = ref_map_pss.get_fn_values();

      for (int j = 0; j < np; j++) // over all integration points
        rhside[part][i] += pt[j][2] * (dx[j] * (fn_dx[j][part] - old[part][0][j]) + dy[j] * (fn_dy[j][part] - old[part][1][j]));
        //rhside[part][i] += pt[j][2] * (dx[j] * (fn[j][part] - old[part][0][j]));
    }

    // solve
    if (e->nvert == 3)
      cholsl(bubble_proj_matrix_tri, nb, bubble_tri_p, rhside[part], rhside[part]);
    else
      cholsl(bubble_proj_matrix_quad, nb, bubble_quad_p, rhside[part], rhside[part]);

    for (int i = 0; i < nb; i++)
      rhs[i][part] = rhside[part][i];
  }
}

void calc_bubble_projection2(Element* e, Nurbs** nurbs, int order, double2* proj)
{
  ref_map_pss.set_active_element(e);

  double2* rhs = proj + e->nvert + e->nvert * (order - 1);
  int mo2 = quad2d.get_max_order();
  int np = quad2d.get_num_points(mo2);
  int qo = e->is_quad() ? make_quad_order(order, order) : order;
  int nb = ref_map_shapeset.get_num_bubbles(qo);

  double2 fn[np];
  double2 fn_dx[np];
  double2 fn_dy[np];
  memset(fn, 0, sizeof(fn));
  memset(fn_dx, 0, sizeof(fn_dx));
  memset(fn_dy, 0, sizeof(fn_dy));

  double* rhside[2];
  for (int i = 0; i < 2; i++) {
    rhside[i] = new double[nb];
    memset(rhside[i], 0, sizeof(double) * nb);
  }

  double* old[2];
  for (int i = 0; i < 2; i++)
  {
    old[i] = new double[np];
    memset(old[i], 0, sizeof(double) * np);
  }

  // compute known part of projection (vertex and edge part)
  old_projection2(e, order, proj, old);

  // fn values and derivatives of both components of nonpolynomial function
  double3* pt = quad2d.get_points(mo2);
  for (int j = 0; j < np; j++)  // over all integration points
    calc_ref_map(e, nurbs, pt[j][0], pt[j][1], fn[j], fn_dx[j], fn_dy[j]);

  for (int part = 0; part < 2; part++)
  {
    for (int i = 0; i < nb; i++) // loop over bubble basis functions
    {
      // derivatives of bubble basis functions in all integration points
      double *bfn;
      int index_i = ref_map_shapeset.get_bubble_indices(qo)[i];
      ref_map_pss.set_active_shape(index_i);
      ref_map_pss.set_quad_order(mo2);
      bfn = ref_map_pss.get_fn_values();

      for (int j = 0; j < np; j++) // over all integration points
        rhside[part][i] += pt[j][2] * (bfn[j] * (fn[j][part] - old[part][j]));
    }

    // solve
    if (e->nvert == 3)
      cholsl(bubble_proj_matrix_tri, nb, bubble_tri_p, rhside[part], rhside[part]);
    else
      cholsl(bubble_proj_matrix_quad, nb, bubble_quad_p, rhside[part], rhside[part]);

    for (int i = 0; i < nb; i++)
      rhs[i][part] = rhside[part][i];
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void ref_map_projection(Element* e, Nurbs** nurbs, int order, double2* proj)
{
  // vertex part
  for (int i = 0; i < e->nvert; i++)
  {
    proj[i][0] = e->vn[i]->x;
    proj[i][1] = e->vn[i]->y;
  }

  // edge part
  for (int edge = 0; edge < e->nvert; edge++)
    calc_edge_projection2(e, edge, nurbs, order, proj);

  //bubble part
  calc_bubble_projection(e, nurbs, order, proj);
}


void CurvMap::update_refmap_coefs(Element* e, int order)
{
  ref_map_pss.set_quad_2d(&quad2d);

  // calculation of projection matrices
  if (edge_proj_matrix == NULL) precalculate_cholesky_projection_matrix_edge();
  if (bubble_proj_matrix_tri == NULL) precalculate_cholesky_projection_matrices_bubble();

  ref_map_pss.set_mode(e->get_mode());
  ref_map_shapeset.set_mode(e->get_mode());

  int nv = e->nvert;
  int ne = order - 1;
  int qo = e->is_quad() ? make_quad_order(order, order) : order;
  int nb = ref_map_shapeset.get_num_bubbles(qo);
  nc = nv + nv*ne + nb;
  // allocate projection coefficients
  if (coefs != NULL) delete [] coefs;
  coefs = new double2[nc];

/*  if (toplevel == false)
  {
    ref_map_pss.set_transform(part);
    Nurbs* nurbs = parent->nurbs;
  }*/

  // calculation of new projection coefficients
  ref_map_projection(e, nurbs, order, coefs);

}
