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

#ifndef __H2D_ELASTICITY_WEAK_FORMS_H
#define __H2D_ELASTICITY_WEAK_FORMS_H

#include "../integrals/h1.h"

/* Default weak form for linear elasticity (Lame equations)
   with Dirichlet and/or zero Neumann BC

   Nonzero Neumann and Newton boundary conditions can be enabled
   by creating a descendant and adding surface forms to it.
*/

namespace WeakFormsElasticity {

  /* Single-component version -- to be used for multimesh assembling */

  class DefaultJacobianElasticity_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianElasticity_0_0(unsigned int i, unsigned int j, double lambda, double mu)
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) { }
    DefaultJacobianElasticity_0_0(unsigned int i, unsigned int j, 
                                          std::string area, double lambda, double mu)
      : WeakForm::MatrixFormVol(i, j, area, HERMES_SYM), lambda(lambda), mu(mu) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      return (lambda + 2*mu) * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                          mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

  private:
      double lambda, mu;
  };

  class DefaultJacobianElasticity_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianElasticity_0_1(unsigned int i, unsigned int j, double lambda, double mu)
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) {}
    DefaultJacobianElasticity_0_1(unsigned int i, unsigned int j, 
                                          std::string area, double lambda, double mu)
      : WeakForm::MatrixFormVol(i, j, area, HERMES_SYM), lambda(lambda), mu(mu) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
    {
      return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
                 mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
            Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
       return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

  private:
    double lambda, mu;
  };

  class DefaultResidualElasticity_0_0 : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualElasticity_0_0(unsigned int i, double lambda, double mu)
      : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu) {}
    DefaultResidualElasticity_0_0(unsigned int i, std::string area, 
                                  double lambda, double mu)
      : WeakForm::VectorFormVol(i,  area), lambda(lambda), mu(mu) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const
    {
      return (2*mu + lambda) * int_dudx_dvdx<Real, Scalar>(n, wt, u_ext[0], v) +
                 mu * int_dudy_dvdy<Real, Scalar>(n, wt, u_ext[0], v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const
    {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const
    {
       return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

  private:
    double lambda, mu;
  };

  class DefaultResidualElasticity_0_1 : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualElasticity_0_1(unsigned int i, double lambda, double mu)
      : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu) {}
    DefaultResidualElasticity_0_1(unsigned int i, std::string area, 
                                            double lambda, double mu)
      : WeakForm::VectorFormVol(i,  area), lambda(lambda), mu(mu) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const
    {
      return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u_ext[1], v) +
                 mu * int_dudx_dvdy<Real, Scalar>(n, wt, u_ext[1], v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const
    {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const
    {
       return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

  private:
    double lambda, mu;
  };

  class DefaultResidualElasticity_1_0 : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualElasticity_1_0(unsigned int i, double lambda, double mu)
      : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu) {}
    DefaultResidualElasticity_1_0(unsigned int i, std::string area, 
                                            double lambda, double mu)
      : WeakForm::VectorFormVol(i,  area), lambda(lambda), mu(mu) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const
    {
      return mu * int_dudy_dvdx<Real, Scalar>(n, wt, u_ext[0], v) +
             lambda * int_dudx_dvdy<Real, Scalar>(n, wt, u_ext[0], v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const
    {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const
    {
       return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

  private:
    double lambda, mu;
  };

  class DefaultResidualElasticity_1_1 : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualElasticity_1_1(unsigned int i, double lambda, double mu)
      : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu) {}
    DefaultResidualElasticity_1_1(unsigned int i, std::string area, 
                                  double lambda, double mu)
      : WeakForm::VectorFormVol(i,  area), lambda(lambda), mu(mu) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const
    {
      return (2*mu + lambda) * int_dudy_dvdy<Real, Scalar>(n, wt, u_ext[1], v) +
                 mu * int_dudx_dvdx<Real, Scalar>(n, wt, u_ext[1], v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const
    {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const
    {
       return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

  private:
    double lambda, mu;
  };

  class DefaultJacobianElasticity_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianElasticity_1_1(unsigned int i, unsigned int j, double lambda, double mu)
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) { }
    DefaultJacobianElasticity_1_1(unsigned int i, unsigned int j, 
                                          std::string area, double lambda, double mu)
      : WeakForm::MatrixFormVol(i, j, area, HERMES_SYM), lambda(lambda), mu(mu) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
    {
      return              mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
             (lambda + 2*mu) * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const
    {
       return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

  private:
    double lambda, mu;
  };

  /* Multicomponent version */

  class DefaultJacobianElasticity_00_11 
    : public WeakForm::MultiComponentMatrixFormVol
  {
  public:
    DefaultJacobianElasticity_00_11(
      Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, double lambda, double mu)
      : WeakForm::MultiComponentMatrixFormVol(coordinates, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) { }
    DefaultJacobianElasticity_00_11(
      Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, std::string area, double lambda, double mu)
      : WeakForm::MultiComponentMatrixFormVol(coordinates, area, HERMES_SYM), lambda(lambda), mu(mu) { }

    template<typename Real, typename Scalar>
    void matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                     Geom<Real> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const {
      Scalar result_0_0 = 0;
      Scalar result_1_1 = 0;
      for (int i = 0; i < n; i++) {
        result_0_0 += wt[i] * ((lambda + 2*mu) * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
        result_1_1 += wt[i] * (mu * u->dx[i] * v->dx[i] + (lambda + 2*mu) * u->dy[i] * v->dy[i]);
      }
      result.push_back(result_0_0);
      result.push_back(result_1_1);
    }

    virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                        Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const {
      matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext, result);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const {
      Hermes::vector<Ord> result;
      matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext, result);

      // Choose the maximum order.
      Ord to_return = result[0];
      if(result[1] > to_return)
        to_return = result[1];
      return to_return;
    }

  private:
    double lambda, mu;
  };

  class DefaultResidualElasticity_00_11 : public WeakForm::MultiComponentVectorFormVol
  {
  public:
    DefaultResidualElasticity_00_11(
      Hermes::vector<unsigned int> coordinates, double lambda, double mu)
      : WeakForm::MultiComponentVectorFormVol(coordinates), lambda(lambda), mu(mu) { }
    DefaultResidualElasticity_00_11(
      Hermes::vector<unsigned int> coordinates, std::string area, double lambda, double mu)
      : WeakForm::MultiComponentVectorFormVol(coordinates, area), lambda(lambda), mu(mu) { }

    template<typename Real, typename Scalar>
    void vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                        Geom<Real> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const {
      Scalar result_0 = 0;
      Scalar result_1 = 0;
      for (int i = 0; i < n; i++) {
        result_0 += wt[i] * ((lambda + 2*mu) * u_ext[0]->dx[i] * v->dx[i] + mu * u_ext[0]->dy[i] * v->dy[i]);
        result_1 += wt[i] * (mu * u_ext[1]->dx[i] * v->dx[i] + (lambda + 2*mu) * u_ext[1]->dy[i] * v->dy[i]);
      }
      result.push_back(result_0);
      result.push_back(result_1);
    }

    virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const {
      vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext, result);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const {
      Hermes::vector<Ord> result;
      vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext, result);

      // Choose the maximum order.
      Ord to_return = result[0];
      if(result[1] > to_return)
        to_return = result[1];
      return to_return;
    }

  private:
    double lambda, mu;
  };
};
#endif
