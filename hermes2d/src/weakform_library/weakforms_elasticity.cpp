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
        : MatrixFormVol<Scalar>(i, j), lambda(lambda), mu(mu)
      {
        this->setSymFlag(HERMES_SYM);
      }

      template<typename Scalar>
      DefaultJacobianElasticity_0_0<Scalar>::DefaultJacobianElasticity_0_0
        (unsigned int i, unsigned int j, std::string area, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j), lambda(lambda), mu(mu)
      {
        this->setSymFlag(HERMES_SYM);
        this->set_area(area);
      }

      template<typename Scalar>
      Scalar DefaultJacobianElasticity_0_0<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const
      {
        return (lambda + 2*mu) * int_dudx_dvdx<double, Scalar>(n, wt, u, v) +
          mu * int_dudy_dvdy<double, Scalar>(n, wt, u, v);
      }

      template<typename Scalar>
      Ord DefaultJacobianElasticity_0_0<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        return (lambda + 2*mu) * int_dudx_dvdx<Ord, Ord>(n, wt, u, v) +
          mu * int_dudy_dvdy<Ord, Ord>(n, wt, u, v);
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultJacobianElasticity_0_0<Scalar>::clone() const
      {
        return new DefaultJacobianElasticity_0_0<Scalar>(this->i, this->j, this->areas[0], this->lambda, this->mu);
      }

      template<typename Scalar>
      DefaultJacobianElasticity_0_1<Scalar>::DefaultJacobianElasticity_0_1
        (unsigned int i, unsigned int j, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j), lambda(lambda), mu(mu)
      {
        this->setSymFlag(HERMES_SYM);
      }

      template<typename Scalar>
      DefaultJacobianElasticity_0_1<Scalar>::DefaultJacobianElasticity_0_1
        (unsigned int i, unsigned int j, std::string area, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j), lambda(lambda), mu(mu)
      {
        this->setSymFlag(HERMES_SYM);
        this->set_area(area);
      }

      template<typename Scalar>
      Scalar DefaultJacobianElasticity_0_1<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const
      {
        return lambda * int_dudy_dvdx<double, Scalar>(n, wt, u, v) +
          mu * int_dudx_dvdy<double, Scalar>(n, wt, u, v);
      }

      template<typename Scalar>
      Ord DefaultJacobianElasticity_0_1<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        return lambda * int_dudy_dvdx<Ord, Ord>(n, wt, u, v) +
          mu * int_dudx_dvdy<Ord, Ord>(n, wt, u, v);
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultJacobianElasticity_0_1<Scalar>::clone() const
      {
        return new DefaultJacobianElasticity_0_1<Scalar>(this->i, this->j, this->areas[0], this->lambda, this->mu);
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
        : VectorFormVol<Scalar>(i), lambda(lambda), mu(mu)
      {
        this->set_area(area);
      }

      template<typename Scalar>
      Scalar DefaultResidualElasticity_0_0<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        return (2*mu + lambda) * int_dudx_dvdx<Scalar, double>(n, wt, u_ext[0], v) +
          mu * int_dudy_dvdy<Scalar, double>(n, wt, u_ext[0], v);
      }

      template<typename Scalar>
      Ord DefaultResidualElasticity_0_0<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        return (2*mu + lambda) * int_dudx_dvdx<Ord, Ord>(n, wt, u_ext[0], v) +
          mu * int_dudy_dvdy<Ord, Ord>(n, wt, u_ext[0], v);
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualElasticity_0_0<Scalar>::clone() const
      {
        return new DefaultResidualElasticity_0_0<Scalar>(this->i, this->areas[0], this->lambda, this->mu);
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
        : VectorFormVol<Scalar>(i), lambda(lambda), mu(mu)
      {
        this->set_area(area);
      }

      template<typename Scalar>
      Scalar DefaultResidualElasticity_0_1<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        return lambda * int_dudy_dvdx<double, Scalar>(n, wt, u_ext[1], v) +
          mu * int_dudx_dvdy<double, Scalar>(n, wt, u_ext[1], v);
      }

      template<typename Scalar>
      Ord DefaultResidualElasticity_0_1<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        return lambda * int_dudy_dvdx<Ord, Ord>(n, wt, u_ext[1], v) +
          mu * int_dudx_dvdy<Ord, Ord>(n, wt, u_ext[1], v);
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualElasticity_0_1<Scalar>::clone() const
      {
        return new DefaultResidualElasticity_0_1<Scalar>(this->i, this->areas[0], this->lambda, this->mu);
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
        : VectorFormVol<Scalar>(i), lambda(lambda), mu(mu)
      {
        this->set_area(area);
      }

      template<typename Scalar>
      Scalar DefaultResidualElasticity_1_0<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        return mu * int_dudy_dvdx<double, Scalar>(n, wt, u_ext[0], v) +
          lambda * int_dudx_dvdy<double, Scalar>(n, wt, u_ext[0], v);
      }

      template<typename Scalar>
      Ord DefaultResidualElasticity_1_0<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        return mu * int_dudy_dvdx<Ord, Ord>(n, wt, u_ext[0], v) +
          lambda * int_dudx_dvdy<Ord, Ord>(n, wt, u_ext[0], v);
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualElasticity_1_0<Scalar>::clone() const
      {
        return new DefaultResidualElasticity_1_0<Scalar>(this->i, this->areas[0], this->lambda, this->mu);
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
        : VectorFormVol<Scalar>(i), lambda(lambda), mu(mu)
      {
        this->set_area(area);
      }

      template<typename Scalar>
      Scalar DefaultResidualElasticity_1_1<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        return (2*mu + lambda) * int_dudy_dvdy<double, Scalar>(n, wt, u_ext[1], v) +
          mu * int_dudx_dvdx<double, Scalar>(n, wt, u_ext[1], v);
      }

      template<typename Scalar>
      Ord DefaultResidualElasticity_1_1<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        return (2*mu + lambda) * int_dudy_dvdy<Ord, Ord>(n, wt, u_ext[1], v) +
          mu * int_dudx_dvdx<Ord, Ord>(n, wt, u_ext[1], v);
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualElasticity_1_1<Scalar>::clone() const
      {
        return new DefaultResidualElasticity_1_1<Scalar>(this->i, this->areas[0], this->lambda, this->mu);
      }

      template<typename Scalar>
      DefaultJacobianElasticity_1_1<Scalar>::DefaultJacobianElasticity_1_1
        (unsigned int i, unsigned int j, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j), lambda(lambda), mu(mu)
      {
      }

      template<typename Scalar>
      DefaultJacobianElasticity_1_1<Scalar>::DefaultJacobianElasticity_1_1
        (unsigned int i, unsigned int j, std::string area, double lambda, double mu)
        : MatrixFormVol<Scalar>(i, j), lambda(lambda), mu(mu)
      {
        this->set_area(area);
      }

      template<typename Scalar>
      Scalar DefaultJacobianElasticity_1_1<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const
      {
        return  mu * int_dudx_dvdx<double, Scalar>(n, wt, u, v) +
          (lambda + 2*mu) * int_dudy_dvdy<double, Scalar>(n, wt, u, v);
      }

      template<typename Scalar>
      Ord DefaultJacobianElasticity_1_1<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        return  mu * int_dudx_dvdx<Ord, Ord>(n, wt, u, v) +
          (lambda + 2*mu) * int_dudy_dvdy<Ord, Ord>(n, wt, u, v);
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultJacobianElasticity_1_1<Scalar>::clone() const
      {
        return new DefaultJacobianElasticity_1_1<Scalar>(this->i, this->j, this->areas[0], this->lambda, this->mu);
      }

      template class HERMES_API DefaultJacobianElasticity_0_0<double>;
      template class HERMES_API DefaultJacobianElasticity_0_1<double>;
      template class HERMES_API DefaultResidualElasticity_0_0<double>;
      template class HERMES_API DefaultResidualElasticity_0_1<double>;
      template class HERMES_API DefaultResidualElasticity_1_0<double>;
      template class HERMES_API DefaultResidualElasticity_1_1<double>;
      template class HERMES_API DefaultJacobianElasticity_1_1<double>;
    };
  }
}