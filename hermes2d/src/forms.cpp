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

#include "forms.h"
#include <complex>

namespace Hermes
{
  namespace Hermes2D
  {
    // Explicit template specializations are needed here, general template<T> T DiscontinuousFunc<T>::zero = T(0) doesn't work.
    template<> Hermes::Ord DiscontinuousFunc<Hermes::Ord>::zero = Hermes::Ord(0);
    template<> double DiscontinuousFunc<double>::zero = 0.0;
    template<> std::complex<double> DiscontinuousFunc<std::complex<double> >::zero = std::complex<double>(0);

    /// Integration order for coordinates, normals and tangents is one.
    Geom<Hermes::Ord>* init_geom_ord()
    {
      Geom<Hermes::Ord>* e = new Geom<Hermes::Ord>;
      static Hermes::Ord x[] = { Hermes::Ord(1) };
      static Hermes::Ord y[] = { Hermes::Ord(1) };

      static Hermes::Ord nx[] = { Hermes::Ord(1) };
      static Hermes::Ord ny[] = { Hermes::Ord(1) };

      static Hermes::Ord tx[] = { Hermes::Ord(1) };
      static Hermes::Ord ty[] = { Hermes::Ord(1) };

      static Hermes::Ord diam = Hermes::Ord(1);

      e->x = x; e->y = y;
      e->nx = nx; e->ny = ny;
      e->tx = tx; e->ty = ty;
      e->diam = diam;
      e->edge_marker = -8888;
      e->elem_marker = -9999;
      return e;
    }

    /// Initialize element marker and coordinates.
    Geom<double>* init_geom_vol(RefMap *rm, const int order)
    {
      Geom<double>* e = new Geom<double>;
      e->diam = rm->get_active_element()->get_diameter();
      e->id = rm->get_active_element()->id;
      e->elem_marker = rm->get_active_element()->marker;
      e->x = rm->get_phys_x(order);
      e->y = rm->get_phys_y(order);
      return e;
    }

    /// Initialize edge marker, coordinates, tangent and normals.
    Geom<double>* init_geom_surf(RefMap *rm, SurfPos* surf_pos, const int order)
    {
      Geom<double>* e = new Geom<double>;
      e->edge_marker = surf_pos->marker;
      e->elem_marker = rm->get_active_element()->marker;
      e->diam = rm->get_active_element()->get_diameter();
      e->id = rm->get_active_element()->en[surf_pos->surf_num]->id;
      e->x = rm->get_phys_x(order);
      e->y = rm->get_phys_y(order);
      double3 *tan;
      tan = rm->get_tangent(surf_pos->surf_num, order);

      Quad2D* quad = rm->get_quad_2d();
      int np = quad->get_num_points(order);
      e->tx = new double [np];
      e->ty = new double [np];
      e->nx = new double [np];
      e->ny = new double [np];
      for (int i = 0; i < np; i++) 
      {
        e->tx[i] = tan[i][0];  e->ty[i] =   tan[i][1];
        e->nx[i] = tan[i][1];  e->ny[i] = - tan[i][0];
      }
      e->orientation = rm->get_active_element()->get_edge_orientation(surf_pos->surf_num);
      return e;
    }

    /// Initialize integration order for function values and derivatives.
    Func<Hermes::Ord>* init_fn_ord(const int order)
    {
      Hermes::Ord *d = new Hermes::Ord(order);

      Func<Hermes::Ord>* f = new Func<Hermes::Ord>(1, 2);
      f->val = d;
      f->dx = f->dy = d;
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
      f->laplace = d;
#endif
      f->val0 = f->val1 = d;
      f->dx0 = f->dx1 = d;
      f->dy0 = f->dy1 = d;
      f->curl = d;
      f->div = d;
      return f;
    }

    /// Transformation of shape functions using reference mapping.
    Func<double>* init_fn(PrecalcShapeset *fu, RefMap *rm, const int order)
    {
      int nc = fu->get_num_components();
      SpaceType space_type = fu->get_space_type();
      Quad2D* quad = fu->get_quad_2d();
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
      if (space_type == HERMES_H1_SPACE)
        fu->set_quad_order(order, H2D_FN_ALL);
      else
#endif
        fu->set_quad_order(order);
      double3* pt = quad->get_points(order);
      int np = quad->get_num_points(order);
      Func<double>* u = new Func<double>(np, nc);

      // H1 space.
      if (space_type == HERMES_H1_SPACE) 
      {
        u->val = new double [np];
        u->dx  = new double [np];
        u->dy  = new double [np];
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
        u->laplace = new double [np];
#endif
        double *fn = fu->get_fn_values();
        double *dx = fu->get_dx_values();
        double *dy = fu->get_dy_values();
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
        double *dxx = fu->get_dxx_values();
        double *dxy = fu->get_dxy_values();
        double *dyy = fu->get_dyy_values();
#endif

        double2x2 *m;
        if(rm->is_jacobian_const()) 
        {
          m = new double2x2[np];
          double2x2 const_inv_ref_map;

          const_inv_ref_map[0][0] = rm->get_const_inv_ref_map()[0][0][0];
          const_inv_ref_map[0][1] = rm->get_const_inv_ref_map()[0][0][1];
          const_inv_ref_map[1][0] = rm->get_const_inv_ref_map()[0][1][0];
          const_inv_ref_map[1][1] = rm->get_const_inv_ref_map()[0][1][1];

          for(int i = 0; i < np; i++) 
          {
            m[i][0][0] = const_inv_ref_map[0][0];
            m[i][0][1] = const_inv_ref_map[0][1];
            m[i][1][0] = const_inv_ref_map[1][0];
            m[i][1][1] = const_inv_ref_map[1][1];
          }

        }
        else
          m = rm->get_inv_ref_map(order);

#ifdef H2D_SECOND_DERIVATIVES_ENABLED
        double3x2 *mm = rm->get_second_ref_map(order);
#endif

#ifdef H2D_SECOND_DERIVATIVES_ENABLED
        for (int i = 0; i < np; i++, m++, mm++)
#else
        for (int i = 0; i < np; i++, m++)
#endif
        {
          u->val[i] = fn[i];
          u->dx[i] = (dx[i] * (*m)[0][0] + dy[i] * (*m)[0][1]);
          u->dy[i] = (dx[i] * (*m)[1][0] + dy[i] * (*m)[1][1]);

#ifdef H2D_SECOND_DERIVATIVES_ENABLED
          double axx = (Hermes::sqr((*m)[0][0]) + Hermes::sqr((*m)[1][0]));
          double ayy = (Hermes::sqr((*m)[0][1]) + Hermes::sqr((*m)[1][1]));
          double axy = 2.0 * ((*m)[0][0]*(*m)[0][1] + (*m)[1][0]*(*m)[1][1]);
          double ax = (*mm)[0][0] + (*mm)[2][0];
          double ay = (*mm)[0][1] + (*mm)[2][1];
          u->laplace[i] = ( dx[i] * ax + dy[i] * ay + dxx[i] * axx + dxy[i] * axy + dyy[i] * ayy );
#endif
        }

        m -= np;
        if(rm->is_jacobian_const())
          delete [] m;
      }
      // Hcurl space.
      else if (space_type == HERMES_HCURL_SPACE) 
      {
        u->val0 = new double [np];
        u->val1 = new double [np];
        u->curl = new double [np];

        double *fn0 = fu->get_fn_values(0);
        double *fn1 = fu->get_fn_values(1);
        double *dx1 = fu->get_dx_values(1);
        double *dy0 = fu->get_dy_values(0);
        double2x2 *m;
        if(rm->is_jacobian_const()) 
        {
          m = new double2x2[np];
          double2x2 const_inv_ref_map;

          const_inv_ref_map[0][0] = rm->get_const_inv_ref_map()[0][0][0];
          const_inv_ref_map[0][1] = rm->get_const_inv_ref_map()[0][0][1];
          const_inv_ref_map[1][0] = rm->get_const_inv_ref_map()[0][1][0];
          const_inv_ref_map[1][1] = rm->get_const_inv_ref_map()[0][1][1];

          for(int i = 0; i < np; i++) 
          {
            m[i][0][0] = const_inv_ref_map[0][0];
            m[i][0][1] = const_inv_ref_map[0][1];
            m[i][1][0] = const_inv_ref_map[1][0];
            m[i][1][1] = const_inv_ref_map[1][1];
          }
        }
        else
          m = rm->get_inv_ref_map(order);
        for (int i = 0; i < np; i++, m++) 
        {
          u->val0[i] = (fn0[i] * (*m)[0][0] + fn1[i] * (*m)[0][1]);
          u->val1[i] = (fn0[i] * (*m)[1][0] + fn1[i] * (*m)[1][1]);
          u->curl[i] = ((*m)[0][0] * (*m)[1][1] - (*m)[1][0] * (*m)[0][1]) * (dx1[i] - dy0[i]);
        }

        m -= np;
        if(rm->is_jacobian_const())
          delete [] m;
      }
      // Hdiv space.
      else if (space_type == HERMES_HDIV_SPACE) 
      {
        u->val0 = new double [np];
        u->val1 = new double [np];

        double *fn0 = fu->get_fn_values(0);
        double *fn1 = fu->get_fn_values(1);
        double *dx0 = fu->get_dx_values(0);
        double *dy1 = fu->get_dy_values(1);
        double2x2 *m;
        if(rm->is_jacobian_const()) 
        {
          m = new double2x2[np];
          double2x2 const_inv_ref_map;

          const_inv_ref_map[0][0] = rm->get_const_inv_ref_map()[0][0][0];
          const_inv_ref_map[0][1] = rm->get_const_inv_ref_map()[0][0][1];
          const_inv_ref_map[1][0] = rm->get_const_inv_ref_map()[0][1][0];
          const_inv_ref_map[1][1] = rm->get_const_inv_ref_map()[0][1][1];

          for(int i = 0; i < np; i++) 
          {
            m[i][0][0] = const_inv_ref_map[0][0];
            m[i][0][1] = const_inv_ref_map[0][1];
            m[i][1][0] = const_inv_ref_map[1][0];
            m[i][1][1] = const_inv_ref_map[1][1];
          }
        }
        else
          m = rm->get_inv_ref_map(order);
        for (int i = 0; i < np; i++, m++) 
        {
          u->val0[i] = (  fn0[i] * (*m)[1][1] - fn1[i] * (*m)[1][0]);
          u->val1[i] = (- fn0[i] * (*m)[0][1] + fn1[i] * (*m)[0][0]);
          u->div[i] = ((*m)[0][0] * (*m)[1][1] - (*m)[1][0] * (*m)[0][1]) * (dx0[i] + dy1[i]);
        }
        m -= np;
        if(rm->is_jacobian_const())
          delete [] m;
      }
      // L2 Space.
      else if (space_type == HERMES_L2_SPACE) 
      {
        // Same as for H1, except that we currently do not have
        // second derivatives of L2 shape functions for triangles.
        u->val = new double [np];
        u->dx  = new double [np];
        u->dy  = new double [np];

        double *fn = fu->get_fn_values();
        double *dx = fu->get_dx_values();
        double *dy = fu->get_dy_values();

        double2x2 *m;
        if(rm->is_jacobian_const()) 
        {
          m = new double2x2[np];
          double2x2 const_inv_ref_map;

          const_inv_ref_map[0][0] = rm->get_const_inv_ref_map()[0][0][0];
          const_inv_ref_map[0][1] = rm->get_const_inv_ref_map()[0][0][1];
          const_inv_ref_map[1][0] = rm->get_const_inv_ref_map()[0][1][0];
          const_inv_ref_map[1][1] = rm->get_const_inv_ref_map()[0][1][1];

          for(int i = 0; i < np; i++) 
          {
            m[i][0][0] = const_inv_ref_map[0][0];
            m[i][0][1] = const_inv_ref_map[0][1];
            m[i][1][0] = const_inv_ref_map[1][0];
            m[i][1][1] = const_inv_ref_map[1][1];
          }
        }
        else
          m = rm->get_inv_ref_map(order);

        for (int i = 0; i < np; i++, m++) 
        {
          u->val[i] = fn[i];
          u->dx[i] = (dx[i] * (*m)[0][0] + dy[i] * (*m)[0][1]);
          u->dy[i] = (dx[i] * (*m)[1][0] + dy[i] * (*m)[1][1]);
        }

        m -= np;
        if(rm->is_jacobian_const())
          delete [] m;
      }
      else
        error("Wrong space type - space has to be either H1, Hcurl, Hdiv or L2");

      return u;
    }

    /// Preparation of mesh functions.
    template<typename Scalar>
    Func<Scalar>* init_fn(MeshFunction<Scalar>*fu, const int order)
    {
      // Sanity checks.
      if (fu == NULL) error("NULL MeshFunction in Func<Scalar>*::init_fn().");
      if (fu->get_mesh() == NULL) error("Uninitialized MeshFunction used.");

      int nc = fu->get_num_components();
      Quad2D* quad = fu->get_quad_2d();
      fu->set_quad_order(order);
      double3* pt = quad->get_points(order);
      int np = quad->get_num_points(order);
      Func<Scalar>* u = new Func<Scalar>(np, nc);

      if (u->nc == 1) 
      {
        u->val = new Scalar [np];
        u->dx  = new Scalar [np];
        u->dy  = new Scalar [np];
        memcpy(u->val, fu->get_fn_values(), np * sizeof(Scalar));
        memcpy(u->dx, fu->get_dx_values(), np * sizeof(Scalar));
        memcpy(u->dy, fu->get_dy_values(), np * sizeof(Scalar));
      }
      else if (u->nc == 2) 
      {
        u->val0 = new Scalar [np];
        u->val1 = new Scalar [np];
        u->curl = new Scalar [np];
        u->div = new Scalar [np];

        memcpy(u->val0, fu->get_fn_values(0), np * sizeof(Scalar));
        memcpy(u->val1, fu->get_fn_values(1), np * sizeof(Scalar));

        Scalar *dx1 = fu->get_dx_values(1);
        Scalar *dy0 = fu->get_dy_values(0);
        for (int i = 0; i < np; i++) u->curl[i] = dx1[i] - dy0[i];

        Scalar *dx0 = fu->get_dx_values(0);
        Scalar *dy1 = fu->get_dy_values(1);
        for (int i = 0; i < np; i++) u->div[i] = dx0[i] + dy1[i];
      }
      return u;
    }

    /// Preparation of solutions.
    template<typename Scalar>
    Func<Scalar>* init_fn(Solution<Scalar>*fu, const int order)
    {
      // Sanity checks.
      if (fu == NULL) error("NULL MeshFunction in Func<Scalar>*::init_fn().");
      if (fu->get_mesh() == NULL) error("Uninitialized MeshFunction used.");

      SpaceType space_type = fu->get_space_type();
      SolutionType sln_type = fu->get_type();

      int nc = fu->get_num_components();
      Quad2D* quad = fu->get_quad_2d();
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
      if (space_type == HERMES_H1_SPACE && sln_type != HERMES_EXACT)
        fu->set_quad_order(order, H2D_FN_ALL);
      else
#endif
        fu->set_quad_order(order);

      double3* pt = quad->get_points(order);
      int np = quad->get_num_points(order);
      Func<Scalar>* u = new Func<Scalar>(np, nc);

      if (u->nc == 1) 
      {
        u->val = new Scalar [np];
        u->dx  = new Scalar [np];
        u->dy  = new Scalar [np];
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
        if (space_type == HERMES_H1_SPACE && sln_type != HERMES_EXACT)
          u->laplace = new Scalar [np];
#endif
        memcpy(u->val, fu->get_fn_values(), np * sizeof(Scalar));
        memcpy(u->dx, fu->get_dx_values(), np * sizeof(Scalar));
        memcpy(u->dy, fu->get_dy_values(), np * sizeof(Scalar));
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
        if (space_type == HERMES_H1_SPACE) 
        {
          if(sln_type == HERMES_SLN) 
          {
            Scalar *dxx = fu->get_dxx_values();
            Scalar *dyy = fu->get_dyy_values();
            for (int i = 0; i < np; i++)
              u->laplace[i] = dxx[i] + dyy[i];
          }
          else if (sln_type == HERMES_CONST)
            memset(u->laplace, 0, np * sizeof(Scalar));
        }
#endif
      }
      else if (u->nc == 2) 
      {
        u->val0 = new Scalar [np];
        u->val1 = new Scalar [np];
        u->curl = new Scalar [np];
        u->div = new Scalar [np];

        memcpy(u->val0, fu->get_fn_values(0), np * sizeof(Scalar));
        memcpy(u->val1, fu->get_fn_values(1), np * sizeof(Scalar));

        Scalar *dx1 = fu->get_dx_values(1);
        Scalar *dy0 = fu->get_dy_values(0);
        for (int i = 0; i < np; i++) u->curl[i] = dx1[i] - dy0[i];

        Scalar *dx0 = fu->get_dx_values(0);
        Scalar *dy1 = fu->get_dy_values(1);
        for (int i = 0; i < np; i++) u->div[i] = dx0[i] + dy1[i];
      }
      return u;
    }

    template Func<double>* init_fn(MeshFunction<double>*fu, const int order);
    template Func<std::complex<double> >* init_fn(MeshFunction<std::complex<double> >*fu, const int order);

    template HERMES_API Func<double>* init_fn(Solution<double>*fu, const int order);
    template HERMES_API Func<std::complex<double> >* init_fn(Solution<std::complex<double> >*fu, const int order);
  }
}