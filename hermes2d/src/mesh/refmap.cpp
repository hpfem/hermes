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

#include "global.h"
#include "mesh.h"
#include "refmap.h"

namespace Hermes
{
  namespace Hermes2D
  {
    RefMap::RefMap() : ref_map_shapeset(H1ShapesetJacobi()), ref_map_pss(PrecalcShapeset(&ref_map_shapeset))
    {
      quad_2d = NULL;
      num_tables = 0;
      cur_node = NULL;
      overflow = NULL;
      set_quad_2d(&g_quad_2d_std); // default quadrature
    }

    RefMap::~RefMap() { free(); }

    /// Sets the quadrature points in which the reference map will be evaluated.
    /// \param quad_2d[in] The quadrature points.
    void set_quad_2d(Quad2D* quad_2d);

    /// Returns the current quadrature points.
    Quad2D* RefMap::get_quad_2d() const
    {
      return quad_2d;
    }

    /// Returns the 1D quadrature for use in surface integrals.
    const Quad1D* RefMap::get_quad_1d() const
    {
      return &quad_1d;
    }

    /// Returns true if the jacobian of the reference map is constant (which
    /// is the case for non-curvilinear triangular elements), false otherwise.
    bool RefMap::is_jacobian_const() const
    {
      return is_const;
    }

    /// Returns the increase in the integration order due to the reference map.
    int RefMap::get_inv_ref_order() const
    {
      return inv_ref_order;
    }

    /// If the jacobian of the reference map is constant, this is the fast
    /// way to obtain it.
    double RefMap::get_const_jacobian() const
    {
      return const_jacobian;
    }

    /// If the reference map is constant, this is the fast way to obtain
    /// its inverse matrix.
    double2x2* RefMap::get_const_inv_ref_map()
    {
      return &const_inv_ref_map;
    }

    /// Returns the jacobian of the reference map precalculated at the integration
    /// points of the specified order. Intended for non-constant jacobian elements.
    double* RefMap::get_jacobian(int order)
    {
      if(cur_node == NULL)
        throw Hermes::Exceptions::Exception("Cur_node == NULL in RefMap - inner algorithms failed");
      if(cur_node->inv_ref_map[order] == NULL)
        calc_inv_ref_map(order);
      return cur_node->jacobian[order];
    }

    /// Returns the inverse matrices of the reference map precalculated at the
    /// integration points of the specified order. Intended for non-constant
    /// jacobian elements.
    double2x2* RefMap::get_inv_ref_map(int order)
    {
      if(cur_node == NULL)
        throw Hermes::Exceptions::Exception("Cur_node == NULL in RefMap - inner algorithms failed");
      if(cur_node->inv_ref_map[order] == NULL)
        calc_inv_ref_map(order);
      return cur_node->inv_ref_map[order];
    }

    /// Returns coefficients for weak forms with second derivatives.
    double3x2* RefMap::get_second_ref_map(int order)
    {
      if(cur_node == NULL)
        throw Hermes::Exceptions::Exception("Cur_node == NULL in RefMap - inner algorithms failed");
      if(cur_node->second_ref_map[order] == NULL) calc_second_ref_map(order);
      return cur_node->second_ref_map[order];
    }

    /// Returns the x-coordinates of the integration points transformed to the
    /// physical domain of the element. Intended for integrals containing spatial
    /// variables.
    double* RefMap::get_phys_x(int order)
    {
      if(cur_node == NULL)
        throw Hermes::Exceptions::Exception("Cur_node == NULL in RefMap - inner algorithms failed");
      if(cur_node->phys_x[order] == NULL) calc_phys_x(order);
      return cur_node->phys_x[order];
    }

    /// Returns he y-coordinates of the integration points transformed to the
    /// physical domain of the element. Intended for integrals containing spatial
    /// variables.
    double* RefMap::get_phys_y(int order)
    {
      if(cur_node == NULL)
        throw Hermes::Exceptions::Exception("Cur_node == NULL in RefMap - inner algorithms failed");
      if(cur_node->phys_y[order] == NULL) calc_phys_y(order);
      return cur_node->phys_y[order];
    }

    /// Returns the triples[x, y, norm] of the tangent to the specified (possibly
    /// curved) edge at the 1D integration points along the edge. The maximum
    /// 1D quadrature rule is used by default, but the user may specify his own
    /// order. In this case, the edge pseudo-order is expected (as returned by
    /// Quad2D::get_edge_points).
    double3* RefMap::get_tangent(int edge, int order)
    {
      if(quad_2d == NULL)
        throw Hermes::Exceptions::Exception("2d quadrature wasn't set.");
      if(order == -1)
        order = quad_2d->get_edge_points(edge, quad_2d->get_max_order(element->get_mode()), element->get_mode());

      // NOTE: Hermes::Order-based caching of geometric data is already employed in DiscreteProblem.
      if(cur_node->tan[edge] != NULL)
      {
        delete [] cur_node->tan[edge];
        cur_node->tan[edge] = NULL;
      }
      calc_tangent(edge, order);

      return cur_node->tan[edge];
    }

    void RefMap::set_quad_2d(Quad2D* quad_2d)
    {
      free();
      this->quad_2d = quad_2d;
      ref_map_pss.set_quad_2d(quad_2d);
    }

    void RefMap::set_active_element(Element* e)
    {
      if(e != element) free();

      ref_map_pss.set_active_element(e);
      num_tables = quad_2d->get_num_tables(e->get_mode());
      assert(num_tables <= H2D_MAX_TABLES);

      if(e == element) return;
      Transformable::set_active_element(e);

      update_cur_node();

      is_const = !element->is_curved() &&
        (element->is_triangle() || is_parallelogram());

      // prepare the shapes and coefficients of the reference map
      int j, k = 0;
      for (unsigned int i = 0; i < e->get_num_surf(); i++)
        indices[k++] = ref_map_shapeset.get_vertex_index(i, e->get_mode());

      // straight-edged element
      if(e->cm == NULL)
      {
        for (unsigned int i = 0; i < e->get_num_surf(); i++)
        {
          lin_coeffs[i][0] = e->vn[i]->x;
          lin_coeffs[i][1] = e->vn[i]->y;
        }
        coeffs = lin_coeffs;
        nc = e->get_num_surf();
      }
      else // curvilinear element - edge and bubble shapes
      {
        int o = e->cm->order;
        for (unsigned int i = 0; i < e->get_num_surf(); i++)
          for (j = 2; j <= o; j++)
            indices[k++] = ref_map_shapeset.get_edge_index(i, 0, j, e->get_mode());

        if(e->is_quad()) o = H2D_MAKE_QUAD_ORDER(o, o);
        memcpy(indices + k, ref_map_shapeset.get_bubble_indices(o, e->get_mode()),
          ref_map_shapeset.get_num_bubbles(o, e->get_mode()) * sizeof(int));

        coeffs = e->cm->coeffs;
        nc = e->cm->nc;
      }

      // calculate the order of the inverse reference map
      if(element->iro_cache == -1 && quad_2d->get_max_order(e->get_mode()) > 1)
      {
        element->iro_cache = is_const ? 0 : calc_inv_ref_order();
      }
      inv_ref_order = element->iro_cache;

      // constant inverse reference map
      if(is_const) calc_const_inv_ref_map(); else const_jacobian = 0.0;
    }

    void RefMap::push_transform(int son)
    {
      Transformable::push_transform(son);
      update_cur_node();
      const_jacobian *= 0.25;
    }

    void RefMap::pop_transform()
    {
      Transformable::pop_transform();
      update_cur_node();
      const_jacobian *= 4;
    }

    void RefMap::calc_inv_ref_map(int order)
    {
      assert(quad_2d != NULL);
      int i, j, np = quad_2d->get_num_points(order, element->get_mode());

      // construct jacobi matrices of the direct reference map for all integration points

      double2x2* m = new double2x2[np];
      memset(m, 0, np * sizeof(double2x2));
      ref_map_pss.force_transform(sub_idx, ctm);
      for (i = 0; i < nc; i++)
      {
        double *dx, *dy;
        ref_map_pss.set_active_shape(indices[i]);
        ref_map_pss.set_quad_order(order);
        ref_map_pss.get_dx_dy_values(dx, dy);
        for (j = 0; j < np; j++)
        {
          m[j][0][0] += coeffs[i][0] * dx[j];
          m[j][0][1] += coeffs[i][0] * dy[j];
          m[j][1][0] += coeffs[i][1] * dx[j];
          m[j][1][1] += coeffs[i][1] * dy[j];
        }
      }

      // calculate the jacobian and inverted matrix
      double trj = get_transform_jacobian();
      double2x2* irm = cur_node->inv_ref_map[order] = new double2x2[np];
      double* jac = cur_node->jacobian[order] = new double[np];
      for (i = 0; i < np; i++)
      {
        jac[i] = (m[i][0][0] * m[i][1][1] - m[i][0][1] * m[i][1][0]);
        double ij = 1.0 / jac[i];
        if(!finite(ij))
        {
          delete [] m;
          throw Hermes::Exceptions::Exception("1/jac[%d] is infinity when calculating inv. ref. map for order %d (jac = %g)", i, order);
        }
        if(ij != ij)
        {
          delete [] m;
          throw Hermes::Exceptions::Exception("1/jac[%d] is NaN when calculating inv. ref. map for order %d (jac = %g)", i, order);
        }

        // invert and transpose the matrix
        irm[i][0][0] =  m[i][1][1] * ij;
        irm[i][0][1] = -m[i][1][0] * ij;
        irm[i][1][0] = -m[i][0][1] * ij;
        irm[i][1][1] =  m[i][0][0] * ij;

        jac[i] *= trj;
      }

      delete [] m;
    }

    void RefMap::calc_second_ref_map(int order)
    {
      assert(quad_2d != NULL);
      int i, j, np = quad_2d->get_num_points(order, element->get_mode());

      double3x2* k = new double3x2[np];
      memset(k, 0, np * sizeof(double3x2));
      ref_map_pss.force_transform(sub_idx, ctm);
      for (i = 0; i < nc; i++)
      {
        double *dxy, *dxx, *dyy;
        ref_map_pss.set_active_shape(indices[i]);
        ref_map_pss.set_quad_order(order, H2D_FN_ALL);
        dxx = ref_map_pss.get_dxx_values();
        dyy = ref_map_pss.get_dyy_values();
        dxy = ref_map_pss.get_dxy_values();
        for (j = 0; j < np; j++)
        {
          k[j][0][0] += coeffs[i][0] * dxx[j];
          k[j][0][1] += coeffs[i][1] * dxx[j];
          k[j][1][0] += coeffs[i][0] * dxy[j];
          k[j][1][1] += coeffs[i][1] * dxy[j];
          k[j][2][0] += coeffs[i][0] * dyy[j];
          k[j][2][1] += coeffs[i][1] * dyy[j];
        }
      }

      double3x2* mm = cur_node->second_ref_map[order] = new double3x2[np];
      double2x2* m = get_inv_ref_map(order);
      for (j = 0; j < np; j++)
      {
        double a, b;
        // coefficients in second derivative with respect to xx
        a = sqr(m[j][0][0])*k[j][0][0] + 2*m[j][0][0]*m[j][0][1]*k[j][1][0] + sqr(m[j][0][1])*k[j][2][0];
        b = sqr(m[j][0][0])*k[j][0][1] + 2*m[j][0][0]*m[j][0][1]*k[j][1][1] + sqr(m[j][0][1])*k[j][2][1];
        mm[j][0][0] = -(a * m[j][0][0] + b * m[j][1][0]); // du/dx
        mm[j][0][1] = -(a * m[j][0][1] + b * m[j][1][1]); // du/dy

        // coefficients in second derivative with respect to xy
        a = m[j][0][0]*m[j][1][0]*k[j][0][0] + (m[j][0][1]*m[j][1][0] + m[j][0][0]*m[j][1][1])*k[j][1][0] + m[j][0][1]*m[j][1][1]*k[j][2][0];
        b = m[j][0][0]*m[j][1][0]*k[j][0][1] + (m[j][0][1]*m[j][1][0] + m[j][0][0]*m[j][1][1])*k[j][1][1] + m[j][0][1]*m[j][1][1]*k[j][2][1];
        mm[j][1][0] = -(a * m[j][0][0] + b * m[j][1][0]); // du/dx
        mm[j][1][1] = -(a * m[j][0][1] + b * m[j][1][1]); // du/dy

        // coefficients in second derivative with respect to yy
        a = sqr(m[j][1][0])*k[j][0][0] + 2*m[j][1][0]*m[j][1][1]*k[j][1][0] + sqr(m[j][1][1])*k[j][2][0];
        b = sqr(m[j][1][0])*k[j][0][1] + 2*m[j][1][0]*m[j][1][1]*k[j][1][1] + sqr(m[j][1][1])*k[j][2][1];
        mm[j][2][0] = -(a * m[j][0][0] + b * m[j][1][0]); // du/dx
        mm[j][2][1] = -(a * m[j][0][1] + b * m[j][1][1]); // du/dy
      }

      delete [] k;
    }

    bool RefMap::is_parallelogram()
    {
      const double eps = 1e-14;
      Element* e = element;
      assert(e->is_quad());
      return fabs(e->vn[2]->x - (e->vn[1]->x + e->vn[3]->x - e->vn[0]->x)) < eps &&
        fabs(e->vn[2]->y - (e->vn[1]->y + e->vn[3]->y - e->vn[0]->y)) < eps;
    }

    void RefMap::calc_const_inv_ref_map()
    {
      if(element == NULL)
        throw Hermes::Exceptions::Exception("The element variable must not be NULL.");
      int k = element->is_triangle() ? 2 : 3;
      double m[2][2] = { { element->vn[1]->x - element->vn[0]->x,  element->vn[k]->x - element->vn[0]->x },
      { element->vn[1]->y - element->vn[0]->y,  element->vn[k]->y - element->vn[0]->y } };

      const_jacobian = 0.25 * (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
      if(const_jacobian <= 0.0)
        throw Hermes::Exceptions::Exception("Element #%d is concave or badly oriented.", element->id);

      double ij = 0.5 / const_jacobian;

      const_inv_ref_map[0][0] =  m[1][1] * ij;
      const_inv_ref_map[1][0] = -m[0][1] * ij;
      const_inv_ref_map[0][1] = -m[1][0] * ij;
      const_inv_ref_map[1][1] =  m[0][0] * ij;

      const_jacobian *= get_transform_jacobian();
    }

    void RefMap::calc_phys_x(int order)
    {
      // transform all x coordinates of the integration points
      int i, j, np = quad_2d->get_num_points(order, element->get_mode());
      double* x = cur_node->phys_x[order] = new double[np];
      memset(x, 0, np * sizeof(double));
      ref_map_pss.force_transform(sub_idx, ctm);
      for (i = 0; i < nc; i++)
      {
        ref_map_pss.set_active_shape(indices[i]);
        ref_map_pss.set_quad_order(order);
        double* fn = ref_map_pss.get_fn_values();
        for (j = 0; j < np; j++)
          x[j] += coeffs[i][0] * fn[j];
      }
    }

    void RefMap::calc_phys_y(int order)
    {
      // transform all y coordinates of the integration points
      int i, j, np = quad_2d->get_num_points(order, element->get_mode());
      double* y = cur_node->phys_y[order] = new double[np];
      memset(y, 0, np * sizeof(double));
      ref_map_pss.force_transform(sub_idx, ctm);
      for (i = 0; i < nc; i++)
      {
        ref_map_pss.set_active_shape(indices[i]);
        ref_map_pss.set_quad_order(order);
        double* fn = ref_map_pss.get_fn_values();
        for (j = 0; j < np; j++)
          y[j] += coeffs[i][1] * fn[j];
      }
    }

    void RefMap::calc_tangent(int edge, int eo)
    {
      int i, j;
      int np = quad_2d->get_num_points(eo, element->get_mode());
      double3* tan = cur_node->tan[edge] = new double3[np];
      int a = edge, b = element->next_vert(edge);

      if(!element->is_curved())
      {
        // straight edges: the tangent at each point is just the edge length
        tan[0][0] = element->vn[b]->x - element->vn[a]->x;
        tan[0][1] = element->vn[b]->y - element->vn[a]->y;
        tan[0][2] = sqrt(sqr(tan[0][0]) + sqr(tan[0][1]));
        double inorm = 1.0 / tan[0][2];
        tan[0][0] *= inorm;
        tan[0][1] *= inorm;
        tan[0][2] *= (edge == 0 || edge == 2) ? ctm->m[0] : ctm->m[1];

        for (i = 1; i < np; i++)
          memcpy(tan + i, tan, sizeof(double3));
      }
      else
      {
        // construct jacobi matrices of the direct reference map at integration points along the edge
        double2x2 m[15];
        assert(np <= 15);
        memset(m, 0, np*sizeof(double2x2));
        ref_map_pss.force_transform(sub_idx, ctm);
        for (i = 0; i < nc; i++)
        {
          double *dx, *dy;
          ref_map_pss.set_active_shape(indices[i]);
          ref_map_pss.set_quad_order(eo);
          ref_map_pss.get_dx_dy_values(dx, dy);
          for (j = 0; j < np; j++)
          {
            m[j][0][0] += coeffs[i][0] * dx[j];
            m[j][0][1] += coeffs[i][0] * dy[j];
            m[j][1][0] += coeffs[i][1] * dx[j];
            m[j][1][1] += coeffs[i][1] * dy[j];
          }
        }

        // multiply them by the vector of the reference edge
        double2* v1 = ref_map_shapeset.get_ref_vertex(a, element->get_mode());
        double2* v2 = ref_map_shapeset.get_ref_vertex(b, element->get_mode());

        double ex = (*v2)[0] - (*v1)[0];
        double ey = (*v2)[1] - (*v1)[1];
        for (i = 0; i < np; i++)
        {
          double3& t = tan[i];
          t[0] = m[i][0][0]*ex + m[i][0][1]*ey;
          t[1] = m[i][1][0]*ex + m[i][1][1]*ey;
          t[2] = sqrt(sqr(t[0]) + sqr(t[1]));
          double inorm = 1.0 / t[2];
          t[0] *= inorm;
          t[1] *= inorm;
          t[2] *= (edge == 0 || edge == 2) ? ctm->m[0] : ctm->m[1];
        }
      }
    }

    int RefMap::calc_inv_ref_order()
    {
      Quad2D* quad = get_quad_2d();
      int i, o, mo = quad->get_max_order(element->get_mode());

      // check first the positivity of the jacobian
      double3* pt = quad->get_points(mo, element->get_mode());
      double2x2* m = get_inv_ref_map(mo);
      double* jac = get_jacobian(mo);
      for (i = 0; i < quad->get_num_points(mo, element->get_mode()); i++)
        if(jac[i] <= 0.0)
          throw Hermes::Exceptions::Exception("Element #%d is concave or badly oriented.", element->id);

      // next, estimate the "exact" value of the typical integral int_grad_u_grad_v
      // (with grad_u == grad_v == (1, 1)) using the maximum integration rule
      double exact1 = 0.0;
      double exact2 = 0.0;
      for (i = 0; i < quad->get_num_points(mo, element->get_mode()); i++, m++)
      {
        exact1 += pt[i][2] * jac[i] * (sqr((*m)[0][0] + (*m)[0][1]) + sqr((*m)[1][0] + (*m)[1][1]));
        exact2 += pt[i][2] / jac[i];
      }
      // find sufficient quadrature degree
      for (o = 0; o < mo; o++)
      {
        pt = quad->get_points(o, element->get_mode());
        m = get_inv_ref_map(o);
        jac = get_jacobian(o);
        double result1 = 0.0;
        double result2 = 0.0;
        for (i = 0; i < quad->get_num_points(o, element->get_mode()); i++, m++)
        {
          result1 += pt[i][2] * jac[i] * (sqr((*m)[0][0] + (*m)[0][1]) + sqr((*m)[1][0] + (*m)[1][1]));
          result2 += pt[i][2] / jac[i] ;
        }
        if((fabs((exact1 - result1) / exact1) < 1e-8) &&
          (fabs((exact2 - result2) / exact2) < 1e-8)) break;
      }
      if(o >= 10)
      {
        this->warn("Element #%d is too distorted (iro ~ %d).", element->id, o);
      }
      return o;
    }

    void RefMap::inv_ref_map_at_point(double xi1, double xi2, double& x, double& y, double2x2& m)
    {
      double2x2 tmp;
      memset(tmp, 0, sizeof(double2x2));
      x = y = 0;
      for (int i = 0; i < nc; i++)
      {
        double val = ref_map_shapeset.get_fn_value(indices[i], xi1, xi2, 0, element->get_mode());
        x += coeffs[i][0] * val;
        y += coeffs[i][1] * val;

        double dx =  ref_map_shapeset.get_dx_value(indices[i], xi1, xi2, 0, element->get_mode());
        double dy =  ref_map_shapeset.get_dy_value(indices[i], xi1, xi2, 0, element->get_mode());
        tmp[0][0] += coeffs[i][0] * dx;
        tmp[0][1] += coeffs[i][0] * dy;
        tmp[1][0] += coeffs[i][1] * dx;
        tmp[1][1] += coeffs[i][1] * dy;
      }

      // inverse matrix
      double jac = tmp[0][0] * tmp[1][1] - tmp[0][1] * tmp[1][0];
      m[0][0] =  tmp[1][1] / jac;
      m[0][1] = -tmp[1][0] / jac;
      m[1][0] = -tmp[0][1] / jac;
      m[1][1] =  tmp[0][0] / jac;
    }

    void RefMap::second_ref_map_at_point(double xi1, double xi2, double& x, double& y, double3x2& mm)
    {
      double3x2 k;
      memset(k, 0, sizeof(double3x2));
      x = y = 0;
      for (int i = 0; i < nc; i++)
      {
        double val = ref_map_shapeset.get_fn_value(indices[i], xi1, xi2, 0, element->get_mode());
        x += coeffs[i][0] * val;
        y += coeffs[i][1] * val;

        double dxy, dxx, dyy;

        dxx = ref_map_shapeset.get_dxx_value(indices[i], xi1, xi2, 0, element->get_mode());
        dxy = ref_map_shapeset.get_dxy_value(indices[i], xi1, xi2, 0, element->get_mode());
        dyy = ref_map_shapeset.get_dxy_value(indices[i], xi1, xi2, 0, element->get_mode());

        k[0][0] += coeffs[i][0] * dxx;
        k[0][1] += coeffs[i][1] * dxx;
        k[1][0] += coeffs[i][0] * dxy;
        k[1][1] += coeffs[i][1] * dxy;
        k[2][0] += coeffs[i][0] * dyy;
        k[2][1] += coeffs[i][1] * dyy;
      }

      double2x2 m;
      this->inv_ref_map_at_point(xi1, xi2, x, y, m);
      double a, b;

      // coefficients in second derivative with respect to xx
      a = sqr(m[0][0])*k[0][0] + 2*m[0][0]*m[0][1]*k[1][0] + sqr(m[0][1])*k[2][0];
      b = sqr(m[0][0])*k[0][1] + 2*m[0][0]*m[0][1]*k[1][1] + sqr(m[0][1])*k[2][1];
      mm[0][0] = -(a * m[0][0] + b * m[1][0]); // du/dx
      mm[0][1] = -(a * m[0][1] + b * m[1][1]); // du/dy

      // coefficients in second derivative with respect to xy
      a = m[0][0]*m[1][0]*k[0][0] + (m[0][1]*m[1][0] + m[0][0]*m[1][1])*k[1][0] + m[0][1]*m[1][1]*k[2][0];
      b = m[0][0]*m[1][0]*k[0][1] + (m[0][1]*m[1][0] + m[0][0]*m[1][1])*k[1][1] + m[0][1]*m[1][1]*k[2][1];
      mm[1][0] = -(a * m[0][0] + b * m[1][0]); // du/dx
      mm[1][1] = -(a * m[0][1] + b * m[1][1]); // du/dy

      // coefficients in second derivative with respect to yy
      a = sqr(m[1][0])*k[0][0] + 2*m[1][0]*m[1][1]*k[1][0] + sqr(m[1][1])*k[2][0];
      b = sqr(m[1][0])*k[0][1] + 2*m[1][0]*m[1][1]*k[1][1] + sqr(m[1][1])*k[2][1];
      mm[2][0] = -(a * m[0][0] + b * m[1][0]); // du/dx
      mm[2][1] = -(a * m[0][1] + b * m[1][1]); // du/dy
    }

    void RefMap::untransform(Element* e, double x, double y, double& xi1, double& xi2)
    {
      const double TOL = 1e-12;

      if(is_const)
      {
        double dx = e->vn[0]->x - x;
        double dy = e->vn[0]->y - y;
        xi1 = -1.0 - (const_inv_ref_map[0][0] * dx + const_inv_ref_map[1][0] * dy);
        xi2 = -1.0 - (const_inv_ref_map[0][1] * dx + const_inv_ref_map[1][1] * dy);
      }
      else // Newton Method
      {
        double xi1_old = 0.0, xi2_old = 0.0;
        double vx, vy;
        double2x2 m;
        int it = 0; // number of Newton iterations
        while (1)
        {
          inv_ref_map_at_point(xi1_old, xi2_old, vx, vy, m);
          xi1 = xi1_old - (m[0][0] * (vx - x) + m[1][0] * (vy - y));
          xi2 = xi2_old - (m[0][1] * (vx - x) + m[1][1] * (vy - y));
          if(fabs(xi1 - xi1_old) < TOL && fabs(xi2 - xi2_old) < TOL) return;
          if(it > 1 && (xi1 > 1.5 || xi2 > 1.5 || xi1 < -1.5 || xi2 < -1.5)) return;
          if(it > 100) { this->warn("Could not find reference coordinates - Newton method did not converge."); return; }
          xi1_old = xi1;
          xi2_old = xi2;
          it++;
        }
      }
    }

    void RefMap::init_node(Node* pp)
    {
      memset(pp->inv_ref_map, 0, num_tables * sizeof(double2x2*));
      memset(pp->second_ref_map, 0, num_tables * sizeof(double3x2*));
      memset(pp->phys_x, 0, num_tables * sizeof(double*));
      memset(pp->phys_y, 0, num_tables * sizeof(double*));
      memset(pp->tan, 0, sizeof(pp->tan));
    }

    void RefMap::free_node(Node* node)
    {
      // destroy all precalculated tables
      for (int i = 0; i < num_tables; i++)
      {
        if(node->inv_ref_map[i] != NULL)
        {
          delete [] node->inv_ref_map[i];
          delete [] node->jacobian[i];
        }
        if(node->second_ref_map[i] != NULL) delete [] node->second_ref_map[i];
        if(node->phys_x[i] != NULL) delete [] node->phys_x[i];
        if(node->phys_y[i] != NULL) delete [] node->phys_y[i];
      }

      for (int i = 0; i < 4; i++)
        if(node->tan[i] != NULL)
          delete [] node->tan[i];

      delete node;
    }

    void RefMap::free()
    {
      std::map<uint64_t, Node*>::iterator it;

      for (it = nodes.begin(); it != nodes.end(); ++it)
        free_node(it->second);
      nodes.clear();
      if(overflow != NULL)
      {
        free_node(overflow); overflow = NULL;
      }
    }

    RefMap::Node* RefMap::handle_overflow()
    {
      if(overflow != NULL)
        free_node(overflow);
      overflow = new Node;
      init_node(overflow);
      return overflow;
    }
  }
}