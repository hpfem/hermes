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

/// \file This file contains definition of methods of classes for form evaluation (hence the name) Geom, and Func.

#include "forms.h"
#include <complex>

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename T>
    Func<T>::Func(int num_gip, int num_comps) : num_gip(num_gip), nc(num_comps)
    {
      val = val0 = val1 = NULL;
      dx = dx0 = dx1 = NULL;
      dy = dy0 = dy1 = NULL;
      curl = NULL;
      div = NULL;
      laplace = NULL;
      sub_idx = 99999;
    }

    template<typename T>
    void Func<T>::subtract(const Func<T>& func)
    {
      if(num_gip != func.num_gip)
        throw Hermes::Exceptions::Exception("Unable to subtract a function due to a different number of integration points (this: %d, other: %d)", num_gip, func.num_gip);
      if(nc != func.nc)
        throw Hermes::Exceptions::Exception("Unable to subtract a function due to a different number of components (this: %d, other: %d)", nc, func.nc);

      subtract(this->val, func.val);
      subtract(this->dx, func.dx);
      subtract(this->dy, func.dy);
      // The following function checks that laplace != NULL.
      subtract(this->laplace, func.laplace);
      if(nc > 1)
      {
        subtract(this->val0, func.val0);
        subtract(this->val1, func.val1);
        subtract(this->dx0, func.dx0);
        subtract(this->dx1, func.dx1);
        subtract(this->dy0, func.dy0);
        subtract(this->dy1, func.dy1);
        subtract(this->curl, func.curl);
        subtract(this->div, func.div);
      }
    };

    template<typename T>
    int Func<T>::get_num_gip() const
    {
      return num_gip;
    }

    template<typename T>
    void Func<T>::subtract(T* attribute, T* other_attribute)
    {
      if(attribute != NULL && other_attribute != NULL)
      {
        for(int i = 0; i < num_gip; i++)
          attribute[i] -= other_attribute[i];
      }
    }

    template<typename T>
    void Func<T>::add(const Func<T>& func)
    {
      if(num_gip != func.num_gip)
        throw Hermes::Exceptions::Exception("Unable to add a function due to a different number of integration points (this: %d, other: %d)", num_gip, func.num_gip);
      if(nc != func.nc)
        throw Hermes::Exceptions::Exception("Unable to add a function due to a different number of components (this: %d, other: %d)", nc, func.nc);

      add(this->val, func.val);
      add(this->dx, func.dx);
      add(this->dy, func.dy);
      // The following function checks that laplace != NULL.
      add(this->laplace, func.laplace);
      if(nc > 1)
      {
        add(this->val0, func.val0);
        add(this->val1, func.val1);
        add(this->dx0, func.dx0);
        add(this->dx1, func.dx1);
        add(this->dy0, func.dy0);
        add(this->dy1, func.dy1);
        add(this->curl, func.curl);
        add(this->div, func.div);
      }
    };

    template<typename T>
    void Func<T>::add(T* attribute, T* other_attribute)
    {
      if(attribute != NULL && other_attribute != NULL)
      {
        if(other_attribute == NULL)
          throw Hermes::Exceptions::Exception("Unable to add a function expansion, the desired attribute is NULL in the other function.");
        for(int i = 0; i < num_gip; i++)
          attribute[i] += other_attribute[i];
      }
    }

    template<typename T>
    void Func<T>::free_ord()
    {
      delete val;
      val = val0 = val1 = NULL;
      dx = dx0 = dx1 = NULL;
      dy = dy0 = dy1 = NULL;
      curl = NULL;
      div = NULL;
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
      laplace = NULL;
#endif
    }

    template<typename T>
    void Func<T>::free_fn()
    {
      delete [] val; val = NULL;
      delete [] dx; dx = NULL;
      delete [] dy; dy = NULL;
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
      delete [] laplace; laplace = NULL;
#endif

      delete [] val0; delete [] val1; val0 = val1 = NULL;
      delete [] dx0;  delete [] dx1; dx0 = dx1 = NULL;
      delete [] dy0;  delete [] dy1; dy0 = dy1 = NULL;
      delete [] curl; curl = NULL;
      delete [] div; div = NULL;
    }

    template<typename T>
    DiscontinuousFunc<T>::DiscontinuousFunc(Func<T>* fn, bool support_on_neighbor, bool reverse) :
    Func<T>(fn->num_gip, fn->nc),
      fn_central(NULL),
      fn_neighbor(NULL),
      reverse_neighbor_side(reverse)
    {
      this->sub_idx = fn->sub_idx;
      if(fn == NULL)
        throw Hermes::Exceptions::Exception("Invalid arguments to DiscontinuousFunc constructor.");
      if(support_on_neighbor)
        fn_neighbor = fn;
      else
        fn_central = fn;
    }

    template<typename T>
    DiscontinuousFunc<T>::DiscontinuousFunc(Func<T>* fn_c, Func<T>* fn_n, bool reverse) :
    Func<T>(fn_c->num_gip, fn_c->nc),
      fn_central(fn_c),
      fn_neighbor(fn_n),
      reverse_neighbor_side(reverse)
    {
      if(fn_c == NULL)
        throw Hermes::Exceptions::NullException(0);
      if(fn_n == NULL)
        throw Hermes::Exceptions::NullException(1);
      if(fn_c->num_gip != fn_n->num_gip || fn_c->nc != fn_n->nc)
        throw Hermes::Exceptions::Exception("DiscontinuousFunc must be formed by two Func's with same number of integration points and components.");
    }

    // Explicit template specializations are needed here, general template<T> T DiscontinuousFunc<T>::zero = T(0) doesn't work.
    template<> Hermes::Ord DiscontinuousFunc<Hermes::Ord>::zero = Hermes::Ord(0);
    template<> double DiscontinuousFunc<double>::zero = 0.0;
    template<> std::complex<double> DiscontinuousFunc<std::complex<double> >::zero = std::complex<double>(0);

    template<typename T>
    T& DiscontinuousFunc<T>::get_val_central(int k) const { return (fn_central != NULL) ? fn_central->val[k] : zero; }

    template<typename T>
    T& DiscontinuousFunc<T>::get_val_neighbor(int k) const { return (fn_neighbor != NULL) ? fn_neighbor->val[ reverse_neighbor_side ? fn_neighbor->num_gip-k-1 : k ] : zero; }

    template<typename T>
    T& DiscontinuousFunc<T>::get_dx_central(int k) const { return (fn_central != NULL) ? fn_central->dx[k] : zero; }

    template<typename T>
    T& DiscontinuousFunc<T>::get_dx_neighbor(int k) const { return (fn_neighbor != NULL) ? fn_neighbor->dx[ reverse_neighbor_side ? fn_neighbor->num_gip-k-1 : k ] : zero; }

    template<typename T>
    T& DiscontinuousFunc<T>::get_dy_central(int k) const { return (fn_central != NULL) ? fn_central->dy[k] : zero; }

    template<typename T>
    T& DiscontinuousFunc<T>::get_dy_neighbor(int k) const { return (fn_neighbor != NULL) ? fn_neighbor->dy[ reverse_neighbor_side ? fn_neighbor->num_gip-k-1 : k ] : zero; }

    template<typename T>
    T& DiscontinuousFunc<T>::get_laplace_central(int k) { return (fn_central != NULL) ? fn_central->laplace[k] : zero; }

    template<typename T>
    T& DiscontinuousFunc<T>::get_laplace_neighbor(int k) { return (fn_neighbor != NULL) ? fn_neighbor->laplace[ reverse_neighbor_side ? fn_neighbor->num_gip-k-1 : k ] : zero; }

    template<typename T>
    void DiscontinuousFunc<T>::subtract(const DiscontinuousFunc<T>& func)
    {
      if(fn_central != NULL && func.fn_central != NULL)
        fn_central->subtract(*func.fn_central);
      if(fn_neighbor != NULL && func.fn_neighbor != NULL)
        fn_neighbor->subtract(*func.fn_neighbor);
    }

    template<typename T>
    void DiscontinuousFunc<T>::free_fn()
    {
      if(fn_central != NULL)
      {
        fn_central->free_fn();
        delete fn_central;
        fn_central = NULL;
      }
      if(fn_neighbor != NULL)
      {
        fn_neighbor->free_fn();
        delete fn_neighbor;
        fn_neighbor = NULL;
      }
    }

    template<typename T>
    void DiscontinuousFunc<T>::free_ord()
    {
      if(fn_central != NULL)
      {
        fn_central->free_ord();
        delete fn_central;
        fn_central = NULL;
      }
      if(fn_neighbor != NULL)
      {
        fn_neighbor->free_ord();
        delete fn_neighbor;
        fn_neighbor = NULL;
      }
    }

    template<typename T>
    Geom<T>::Geom()
    {
      elem_marker = -1;
      edge_marker = -1;
      id = 0;
      isurf = 4;
      x = y = NULL;
      nx = ny = NULL;
      tx = ty = NULL;
    }

    template<typename T>
    void Geom<T>::free()
    {
      delete [] tx;    delete [] ty;
      delete [] nx;    delete [] ny;
    }

    template<typename T>
    ExtData<T>::ExtData()
    {
      nf = 0;
      fn = NULL;
    }

    template<typename T>
    void ExtData<T>::free()
    {
      for (int i = 0; i < nf; i++)
      {
        fn[i]->free_fn();
        delete fn[i];
      }
      delete [] fn;
    }

    template<typename T>
    void ExtData<T>::free_ord()
    {
      delete [] fn;
    }

    template<typename T>
    InterfaceGeom<T>::InterfaceGeom(Geom<T>* geom, int n_marker, int n_id, T n_diam) :
    Geom<T>(), neighb_marker(n_marker), neighb_id(n_id), neighb_diam(n_diam)
    {
      // Let this class expose the standard Geom interface.
      this->edge_marker = geom->edge_marker;
      this->elem_marker = geom->elem_marker;
      this->id = geom->id;
      this->isurf = geom->isurf;
      this->diam = geom->diam;
      this->area = geom->area;
      this->x = geom->x;
      this->y = geom->y;
      this->tx = geom->tx;
      this->ty = geom->ty;
      this->nx = geom->nx;
      this->ny = geom->ny;
      this->orientation = geom->orientation;
      this->wrapped_geom = geom;
    }

    template<typename T>
    void InterfaceGeom<T>::free()
    {
      wrapped_geom->free();
      delete wrapped_geom;
    }

    template<typename T>
    void InterfaceGeom<T>::free_ord()
    {
      delete wrapped_geom;
    }

    template<typename T>
    int InterfaceGeom<T>::get_neighbor_marker() const
    {
      return neighb_marker;
    }

    template<typename T>
    int InterfaceGeom<T>::get_neighbor_id() const
    {
      return neighb_id;
    }

    template<typename T>
    T InterfaceGeom<T>::get_neighbor_diam() const
    {
      return neighb_diam;
    }

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

      return e;
    }

    Geom<double>* init_geom_vol(RefMap *rm, const int order)
    {
      Geom<double>* e = new Geom<double>;
      e->diam = rm->get_active_element()->get_diameter();
      e->area = rm->get_active_element()->get_area();
      e->id = rm->get_active_element()->id;
      e->elem_marker = rm->get_active_element()->marker;
      e->x = rm->get_phys_x(order);
      e->y = rm->get_phys_y(order);
      return e;
    }

    Geom<double>* init_geom_surf(RefMap *rm, int isurf, int marker, const int order, double3*& tan)
    {
      Geom<double>* e = new Geom<double>;
      e->edge_marker = marker;
      e->elem_marker = rm->get_active_element()->marker;
      e->diam = rm->get_active_element()->get_diameter();
      e->area = rm->get_active_element()->get_area();
      e->id = rm->get_active_element()->id;
      e->isurf = isurf;
      e->x = rm->get_phys_x(order);
      e->y = rm->get_phys_y(order);

      tan = rm->get_tangent(isurf, order);

      Quad2D* quad = rm->get_quad_2d();
      int np = quad->get_num_points(order, rm->get_active_element()->get_mode());
      e->tx = new double[np];
      e->ty = new double[np];
      e->nx = new double[np];
      e->ny = new double[np];
      for (int i = 0; i < np; i++)
      {
        e->tx[i] = tan[i][0];  e->ty[i] =   tan[i][1];
        e->nx[i] = tan[i][1];  e->ny[i] = - tan[i][0];
      }
      e->orientation = rm->get_active_element()->get_edge_orientation(isurf);
      return e;
    }

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

    Func<double>* init_fn(PrecalcShapeset *fu, RefMap *rm, const int order)
    {
      int nc = fu->get_num_components();
      SpaceType space_type = fu->get_space_type();
      Quad2D* quad = fu->get_quad_2d();
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
      if(space_type == HERMES_H1_SPACE)
        fu->set_quad_order(order, H2D_FN_ALL);
      else
#endif
        fu->set_quad_order(order);
      double3* pt = quad->get_points(order, rm->get_active_element()->get_mode());
      int np = quad->get_num_points(order, rm->get_active_element()->get_mode());
      Func<double>* u = new Func<double>(np, nc);

      // H1 space.
      if(space_type == HERMES_H1_SPACE)
      {
        u->val = new double[np];
        u->dx  = new double[np];
        u->dy  = new double[np];
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
        u->laplace = new double[np];
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
      else if(space_type == HERMES_HCURL_SPACE)
      {
        u->val0 = new double[np];
        u->val1 = new double[np];
        u->curl = new double[np];

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
      else if(space_type == HERMES_HDIV_SPACE)
      {
        u->val0 = new double[np];
        u->val1 = new double[np];

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
      else if(space_type == HERMES_L2_SPACE)
      {
        // Same as for H1, except that we currently do not have
        // second derivatives of L2 shape functions for triangles.
        u->val = new double[np];
        u->dx  = new double[np];
        u->dy  = new double[np];

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
        throw Hermes::Exceptions::Exception("Wrong space type - space has to be either H1, Hcurl, Hdiv or L2");

      u->sub_idx = rm->get_transform();
      return u;
    }

    template<typename Scalar>
    Func<Scalar>* init_fn(MeshFunction<Scalar>*fu, const int order)
    {
      // Sanity checks.
      if(fu == NULL) throw Hermes::Exceptions::Exception("NULL MeshFunction in Func<Scalar>*::init_fn().");
      if(fu->get_mesh() == NULL) throw Hermes::Exceptions::Exception("Uninitialized MeshFunction used.");

      int nc = fu->get_num_components();
      Quad2D* quad = fu->get_quad_2d();
      fu->set_quad_order(order);
      double3* pt = quad->get_points(order, fu->get_active_element()->get_mode());
      int np = quad->get_num_points(order, fu->get_active_element()->get_mode());
      Func<Scalar>* u = new Func<Scalar>(np, nc);

      if(u->nc == 1)
      {
        u->val = new Scalar[np];
        u->dx  = new Scalar[np];
        u->dy  = new Scalar[np];
        memcpy(u->val, fu->get_fn_values(), np * sizeof(Scalar));
        memcpy(u->dx, fu->get_dx_values(), np * sizeof(Scalar));
        memcpy(u->dy, fu->get_dy_values(), np * sizeof(Scalar));
      }
      else if(u->nc == 2)
      {
        u->val0 = new Scalar[np];
        u->val1 = new Scalar[np];
        u->curl = new Scalar[np];
        u->div = new Scalar[np];

        memcpy(u->val0, fu->get_fn_values(0), np * sizeof(Scalar));
        memcpy(u->val1, fu->get_fn_values(1), np * sizeof(Scalar));

        Scalar *dx1 = fu->get_dx_values(1);
        Scalar *dy0 = fu->get_dy_values(0);
        for (int i = 0; i < np; i++) u->curl[i] = dx1[i] - dy0[i];

        Scalar *dx0 = fu->get_dx_values(0);
        Scalar *dy1 = fu->get_dy_values(1);
        for (int i = 0; i < np; i++) u->div[i] = dx0[i] + dy1[i];
      }
      u->sub_idx = fu->get_transform();
      return u;
    }

    template<typename Scalar>
    Func<Scalar>* init_fn(Solution<Scalar>*fu, const int order)
    {
      // Sanity checks.
      if(fu == NULL) throw Hermes::Exceptions::Exception("NULL MeshFunction in Func<Scalar>*::init_fn().");
      if(fu->get_mesh() == NULL) throw Hermes::Exceptions::Exception("Uninitialized MeshFunction used.");

      SpaceType space_type = fu->get_space_type();
      SolutionType sln_type = fu->get_type();

      int nc = fu->get_num_components();
      Quad2D* quad = fu->get_quad_2d();
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
      if(space_type == HERMES_H1_SPACE && sln_type != HERMES_EXACT)
        fu->set_quad_order(order, H2D_FN_ALL);
      else
#endif
        fu->set_quad_order(order);

      double3* pt = quad->get_points(order, fu->get_active_element()->get_mode());
      int np = quad->get_num_points(order, fu->get_active_element()->get_mode());
      Func<Scalar>* u = new Func<Scalar>(np, nc);

      if(u->nc == 1)
      {
        u->val = new Scalar[np];
        u->dx  = new Scalar[np];
        u->dy  = new Scalar[np];
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
        if(space_type == HERMES_H1_SPACE && sln_type != HERMES_EXACT)
          u->laplace = new Scalar[np];
#endif
        memcpy(u->val, fu->get_fn_values(), np * sizeof(Scalar));
        memcpy(u->dx, fu->get_dx_values(), np * sizeof(Scalar));
        memcpy(u->dy, fu->get_dy_values(), np * sizeof(Scalar));
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
        if(space_type == HERMES_H1_SPACE)
        {
          if(sln_type == HERMES_SLN)
          {
            Scalar *dxx = fu->get_dxx_values();
            Scalar *dyy = fu->get_dyy_values();
            for (int i = 0; i < np; i++)
              u->laplace[i] = dxx[i] + dyy[i];
          }
        }
#endif
      }
      else if(u->nc == 2)
      {
        u->val0 = new Scalar[np];
        u->val1 = new Scalar[np];
        u->curl = new Scalar[np];
        u->div = new Scalar[np];

        memcpy(u->val0, fu->get_fn_values(0), np * sizeof(Scalar));
        memcpy(u->val1, fu->get_fn_values(1), np * sizeof(Scalar));

        Scalar *dx1 = fu->get_dx_values(1);
        Scalar *dy0 = fu->get_dy_values(0);
        for (int i = 0; i < np; i++) u->curl[i] = dx1[i] - dy0[i];

        Scalar *dx0 = fu->get_dx_values(0);
        Scalar *dy1 = fu->get_dy_values(1);
        for (int i = 0; i < np; i++) u->div[i] = dx0[i] + dy1[i];
      }
      u->sub_idx = fu->get_transform();
      return u;
    }

    template Func<double>* init_fn(MeshFunction<double>*fu, const int order);
    template Func<std::complex<double> >* init_fn(MeshFunction<std::complex<double> >*fu, const int order);

    template HERMES_API Func<double>* init_fn(Solution<double>*fu, const int order);
    template HERMES_API Func<std::complex<double> >* init_fn(Solution<std::complex<double> >*fu, const int order);

    template class HERMES_API Func<Hermes::Ord>;
    template class HERMES_API Func<double>;
    template class HERMES_API Func<std::complex<double> >;
    template class HERMES_API DiscontinuousFunc<Hermes::Ord>;
    template class HERMES_API DiscontinuousFunc<double>;
    template class HERMES_API DiscontinuousFunc<std::complex<double> >;
    template class HERMES_API Geom<Hermes::Ord>;
    template class HERMES_API Geom<double>;
    template class HERMES_API InterfaceGeom<Hermes::Ord>;
    template class HERMES_API InterfaceGeom<double>;
    template class HERMES_API ExtData<Hermes::Ord>;
    template class HERMES_API ExtData<double>;
    template class HERMES_API ExtData<std::complex<double> >;
  }
}