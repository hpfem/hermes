// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "weakforms_elasticity.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace WeakFormsElasticity 
    {
      template<typename Scalar>
      DefaultJacobianElasticity_0_0<Scalar>::DefaultJacobianElasticity_0_0
        (unsigned int i, unsigned int j, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      DefaultJacobianElasticity_0_0<Scalar>::DefaultJacobianElasticity_0_0
        (unsigned int i, unsigned int j, std::string area, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j, area, HERMES_SYM), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      Scalar DefaultJacobianElasticity_0_0<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const 
      {
        return (lambda + 2*mu) * int_dudx_dvdx<double, Scalar>(n, wt, u, v) +
          mu * int_dudy_dvdy<double, Scalar>(n, wt, u, v);
      }

      template<typename Scalar>
      Ord DefaultJacobianElasticity_0_0<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        return (lambda + 2*mu) * int_dudx_dvdx<Ord, Ord>(n, wt, u, v) +
          mu * int_dudy_dvdy<Ord, Ord>(n, wt, u, v);
      }

      template<typename Scalar>
      DefaultJacobianElasticity_0_1<Scalar>::DefaultJacobianElasticity_0_1
        (unsigned int i, unsigned int j, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      DefaultJacobianElasticity_0_1<Scalar>::DefaultJacobianElasticity_0_1
        (unsigned int i, unsigned int j, std::string area, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j, area, HERMES_SYM), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      Scalar DefaultJacobianElasticity_0_1<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const
      {
        return lambda * int_dudy_dvdx<double, Scalar>(n, wt, u, v) +
          mu * int_dudx_dvdy<double, Scalar>(n, wt, u, v);
      }

      template<typename Scalar>
      Ord DefaultJacobianElasticity_0_1<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        return lambda * int_dudy_dvdx<Ord, Ord>(n, wt, u, v) +
          mu * int_dudx_dvdy<Ord, Ord>(n, wt, u, v);
      }

      template<typename Scalar>
      DefaultResidualElasticity_0_0<Scalar>::DefaultResidualElasticity_0_0
        (unsigned int i, double lambda, double mu)
        : VectorFormVol<Scalar>(i), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      DefaultResidualElasticity_0_0<Scalar>::DefaultResidualElasticity_0_0
        (unsigned int i, std::string area, double lambda, double mu)
        : VectorFormVol<Scalar>(i,  area), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      Scalar DefaultResidualElasticity_0_0<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<Scalar> *ext) const
      {
        return (2*mu + lambda) * int_dudx_dvdx<Scalar, double>(n, wt, u_ext[0], v) +
          mu * int_dudy_dvdy<Scalar, double>(n, wt, u_ext[0], v);
      }

      template<typename Scalar>
      Ord DefaultResidualElasticity_0_0<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
        Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        return (2*mu + lambda) * int_dudx_dvdx<Ord, Ord>(n, wt, u_ext[0], v) +
          mu * int_dudy_dvdy<Ord, Ord>(n, wt, u_ext[0], v);
      }

      template<typename Scalar>
      DefaultResidualElasticity_0_1<Scalar>::DefaultResidualElasticity_0_1
        (unsigned int i, double lambda, double mu)
        : VectorFormVol<Scalar>(i), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      DefaultResidualElasticity_0_1<Scalar>::DefaultResidualElasticity_0_1
        (unsigned int i, std::string area, double lambda, double mu)
        : VectorFormVol<Scalar>(i,  area), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      Scalar DefaultResidualElasticity_0_1<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<Scalar> *ext) const
      {
        return lambda * int_dudy_dvdx<double, Scalar>(n, wt, u_ext[1], v) +
          mu * int_dudx_dvdy<double, Scalar>(n, wt, u_ext[1], v);
      }

      template<typename Scalar>
      Ord DefaultResidualElasticity_0_1<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
        Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        return lambda * int_dudy_dvdx<Ord, Ord>(n, wt, u_ext[1], v) +
          mu * int_dudx_dvdy<Ord, Ord>(n, wt, u_ext[1], v);
      }

      template<typename Scalar>
      DefaultResidualElasticity_1_0<Scalar>::DefaultResidualElasticity_1_0
        (unsigned int i, double lambda, double mu)
        : VectorFormVol<Scalar>(i), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      DefaultResidualElasticity_1_0<Scalar>::DefaultResidualElasticity_1_0
        (unsigned int i, std::string area, double lambda, double mu)
        : VectorFormVol<Scalar>(i,  area), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      Scalar DefaultResidualElasticity_1_0<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<Scalar> *ext) const
      {
        return mu * int_dudy_dvdx<double, Scalar>(n, wt, u_ext[0], v) +
          lambda * int_dudx_dvdy<double, Scalar>(n, wt, u_ext[0], v);
      }

      template<typename Scalar>
      Ord DefaultResidualElasticity_1_0<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
        Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        return mu * int_dudy_dvdx<Ord, Ord>(n, wt, u_ext[0], v) +
          lambda * int_dudx_dvdy<Ord, Ord>(n, wt, u_ext[0], v);
      }


      template<typename Scalar>
      DefaultResidualElasticity_1_1<Scalar>::DefaultResidualElasticity_1_1
        (unsigned int i, double lambda, double mu)
        : VectorFormVol<Scalar>(i), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      DefaultResidualElasticity_1_1<Scalar>::DefaultResidualElasticity_1_1
        (unsigned int i, std::string area, double lambda, double mu)
        : VectorFormVol<Scalar>(i,  area), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      Scalar DefaultResidualElasticity_1_1<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<Scalar> *ext) const
      {
        return (2*mu + lambda) * int_dudy_dvdy<double, Scalar>(n, wt, u_ext[1], v) +
          mu * int_dudx_dvdx<double, Scalar>(n, wt, u_ext[1], v);
      }

      template<typename Scalar>
      Ord DefaultResidualElasticity_1_1<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
        Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        return (2*mu + lambda) * int_dudy_dvdy<Ord, Ord>(n, wt, u_ext[1], v) +
          mu * int_dudx_dvdx<Ord, Ord>(n, wt, u_ext[1], v);
      }

      template<typename Scalar>
      DefaultJacobianElasticity_1_1<Scalar>::DefaultJacobianElasticity_1_1
        (unsigned int i, unsigned int j, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      DefaultJacobianElasticity_1_1<Scalar>::DefaultJacobianElasticity_1_1
        (unsigned int i, unsigned int j, std::string area, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j, area, HERMES_SYM), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      Scalar DefaultJacobianElasticity_1_1<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const
      {
        return  mu * int_dudx_dvdx<double, Scalar>(n, wt, u, v) +
          (lambda + 2*mu) * int_dudy_dvdy<double, Scalar>(n, wt, u, v);
      }

      template<typename Scalar>
      Ord DefaultJacobianElasticity_1_1<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        return  mu * int_dudx_dvdx<Ord, Ord>(n, wt, u, v) +
          (lambda + 2*mu) * int_dudy_dvdy<Ord, Ord>(n, wt, u, v);
      }

      template<typename Scalar>
      DefaultJacobianElasticity_00_11<Scalar>::DefaultJacobianElasticity_00_11
        (Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, double lambda, double mu)
        : MultiComponentMatrixFormVol<Scalar>(coordinates, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      DefaultJacobianElasticity_00_11<Scalar>::DefaultJacobianElasticity_00_11
        (Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, std::string area, double lambda, double mu)
        : MultiComponentMatrixFormVol<Scalar>(coordinates, area, HERMES_SYM), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      void DefaultJacobianElasticity_00_11<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const 
      {
        Scalar result_0_0 = 0;
        Scalar result_1_1 = 0;
        for (int i = 0; i < n; i++) 
        {
          result_0_0 += wt[i] * ((lambda + 2*mu) * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
          result_1_1 += wt[i] * (mu * u->dx[i] * v->dx[i] + (lambda + 2*mu) * u->dy[i] * v->dy[i]);
        }
        result.push_back(result_0_0);
        result.push_back(result_1_1);
      }

      template<typename Scalar>
      Ord DefaultJacobianElasticity_00_11<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Hermes::vector<Ord> result;
        Ord result_0_0(0);
        Ord result_1_1(0);
        for (int i = 0; i < n; i++) 
        {
          result_0_0 += wt[i] * ((lambda + 2*mu) * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
          result_1_1 += wt[i] * (mu * u->dx[i] * v->dx[i] + (lambda + 2*mu) * u->dy[i] * v->dy[i]);
        }
        result.push_back(result_0_0);
        result.push_back(result_1_1);

        // Choose the maximum order.
        Ord to_return = result[0];
        if(result[1] > to_return)
          to_return = result[1];
        return to_return;
      }


      template<typename Scalar>
      DefaultResidualElasticity_00_11<Scalar>::DefaultResidualElasticity_00_11
        (Hermes::vector<unsigned int> coordinates, double lambda, double mu)
        : MultiComponentVectorFormVol<Scalar>(coordinates), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      DefaultResidualElasticity_00_11<Scalar>::DefaultResidualElasticity_00_11
        (Hermes::vector<unsigned int> coordinates, std::string area, double lambda, double mu)
        : MultiComponentVectorFormVol<Scalar>(coordinates, area), lambda(lambda), mu(mu) 
      {
      }

      template<typename Scalar>
      void DefaultResidualElasticity_00_11<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const 
      {
        Scalar result_0 = 0;
        Scalar result_1 = 0;
        for (int i = 0; i < n; i++) 
        {
          result_0 += wt[i] * ((lambda + 2*mu) * u_ext[0]->dx[i] * v->dx[i] + mu * u_ext[0]->dy[i] * v->dy[i]);
          result_1 += wt[i] * (mu * u_ext[1]->dx[i] * v->dx[i] + (lambda + 2*mu) * u_ext[1]->dy[i] * v->dy[i]);
        }
        result.push_back(result_0);
        result.push_back(result_1);
      }

      template<typename Scalar>
      Ord DefaultResidualElasticity_00_11<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Hermes::vector<Ord> result;
        Ord result_0(0);
        Ord result_1(0);
        for (int i = 0; i < n; i++) 
        {
          result_0 += wt[i] * ((lambda + 2*mu) * u_ext[0]->dx[i] * v->dx[i] + mu * u_ext[0]->dy[i] * v->dy[i]);
          result_1 += wt[i] * (mu * u_ext[1]->dx[i] * v->dx[i] + (lambda + 2*mu) * u_ext[1]->dy[i] * v->dy[i]);
        }
        result.push_back(result_0);
        result.push_back(result_1);

        // Choose the maximum order.
        Ord to_return = result[0];
        if(result[1] > to_return)
          to_return = result[1];
        return to_return;
      }

      template class HERMES_API DefaultJacobianElasticity_0_0<double>;
      template class HERMES_API DefaultJacobianElasticity_0_1<double>;
      template class HERMES_API DefaultResidualElasticity_0_0<double>;
      template class HERMES_API DefaultResidualElasticity_0_1<double>;
      template class HERMES_API DefaultResidualElasticity_1_0<double>;
      template class HERMES_API DefaultResidualElasticity_1_1<double>;
      template class HERMES_API DefaultJacobianElasticity_1_1<double>;
      template class HERMES_API DefaultJacobianElasticity_00_11<double>;
      template class HERMES_API DefaultResidualElasticity_00_11<double>;
    };
  }
}