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

#ifndef __H2D_SAMPLE_WEAK_FORMS_H
#define __H2D_SAMPLE_WEAK_FORMS_H

#include "../integrals/integrals_h1.h"

/* Default vector form corresponding to constant right-hand side */

class VectorFormConstant : public WeakForm::VectorFormVol
{
public:
  VectorFormConstant(int i, double const_f) : WeakForm::VectorFormVol(i), const_f(const_f) { }
  VectorFormConstant(int i, std::string area, double const_f) : WeakForm::VectorFormVol(i, area), const_f(const_f) { }

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
               Geom<double> *e, ExtData<scalar> *ext) {
    return const_f * int_v<scalar, scalar>(n, wt, v);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) {
    return v->val[0]; // Make the parser return the polynomial degree of 'v'.
  }

private:
  // Constant right hand side.
  double const_f;
};

/* Default weak form for the Laplace equation -Laplace u = 0 
   with constant Dirichlet or Neumann BC

   Nonconstant Neumann and Newton boundary conditions can be enabled
   by creating a descendant and adding volume and surface forms to it.
*/

class WeakFormLaplace : public WeakForm
{
public:
  WeakFormLaplace(unsigned int neq = 1) : WeakForm(neq) { }

  // Constant matrix volume form
  class MatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVol(int i, int j, double const_f = 1.0) : WeakForm::MatrixFormVol(i, j, HERMES_SYM), const_f(const_f) { }
    MatrixFormVol(int i, int j, std::string area, double const_f = 1.0) : WeakForm::MatrixFormVol(i, j, HERMES_SYM, area), const_f(const_f) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return const_f * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

  private:
    double const_f;
  };

  // Constant Newton matrix surface form
  class MatrixFormSurfNewton : public WeakForm::MatrixFormSurf
  {
  public:
    MatrixFormSurfNewton(int i, int j, double const_f) : WeakForm::MatrixFormSurf(i, j), const_f(const_f) { }
    MatrixFormSurfNewton(int i, int j, std::string area, double const_f) : WeakForm::MatrixFormSurf(i, j, area), const_f(const_f) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return const_f * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
        return matrix_form_surf<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
        return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

  private:
      double const_f;
  };

  // Constant Newton vector surface form
  class VectorFormSurfNewton : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurfNewton(int i, double const_f) : WeakForm::VectorFormSurf(i), const_f(const_f) { }
    VectorFormSurfNewton(int i, std::string area, double const_f) : WeakForm::VectorFormSurf(i, area), const_f(const_f) { }

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return const_f * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
        return vector_form_surf<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
        return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

  private:
      double const_f;
  };

  // Constant Neumann vector surface form
  class VectorFormSurfNeumann : public VectorFormSurfNewton
  {
  public:
    VectorFormSurfNeumann(int i, double const_f) : VectorFormSurfNewton(i, const_f) { }
    VectorFormSurfNeumann(int i, std::string area, double const_f) : VectorFormSurfNewton(i, area, const_f) { }
  };
};

/* Default weak form for the Poisson equation -Laplace u = f
   with constant Dirichlet and Neumann BC

   Nonconstant Neumann and Newton boundary conditions can be enabled
   by creating a descendant and adding volume and surface forms to it.
*/

class WeakFormPoisson : public WeakFormLaplace
{
public:
  WeakFormPoisson(unsigned int neq = 1) : WeakFormLaplace(neq) { }

  // Constant right side
  class VectorFormVol : public VectorFormConstant
  {
  public:
    VectorFormVol(int i, double const_f) : VectorFormConstant(i, const_f) { }
    VectorFormVol(int i, std::string area, double const_f) : VectorFormConstant(i, area, const_f)  { }
  };
};

/* Default (vector-valued) weak form for linear elasticity (Lame equations)  
   with Dirichlet and/or zero Neumann BC (just volumetric forms).

   Nonzero Neumann and Newton boundary conditions can be enabled 
   by creating a descendant and adding surface forms to it. 
*/

class WeakFormLinearElasticity : public WeakForm
{
public:
  WeakFormLinearElasticity(double E, double nu, double rho_g) : WeakForm(2)
  {
    double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
    double mu = E / (2*(1 + nu));
    
    add_matrix_form(new MatrixFormVolLinearElasticity_0_0(lambda, mu));
    add_matrix_form(new MatrixFormVolLinearElasticity_0_1(lambda, mu)); 
    add_matrix_form(new MatrixFormVolLinearElasticity_1_1(lambda, mu));
    add_vector_form(new VectorFormGravity(rho_g));                   // gravity loading
  }

private:
  class MatrixFormVolLinearElasticity_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolLinearElasticity_0_0(double lambda, double mu) 
      : WeakForm::MatrixFormVol(0, 0, HERMES_SYM), lambda(lambda), mu(mu) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return (lambda + 2*mu) * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                          mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext)
    {
       return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Members.
    double lambda, mu;
  };

  class MatrixFormVolLinearElasticity_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolLinearElasticity_0_1(double lambda, double mu) 
            : WeakForm::MatrixFormVol(0, 1, HERMES_SYM), lambda(lambda), mu(mu) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
                 mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
            Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
       return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Members.
    double lambda, mu;
  };

  class MatrixFormVolLinearElasticity_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolLinearElasticity_1_1(double lambda, double mu) 
            : WeakForm::MatrixFormVol(1, 1, HERMES_SYM), lambda(lambda), mu(mu) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return              mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
             (lambda + 2*mu) * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext)
    {
       return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Members.
    double lambda, mu;
  };

  class VectorFormGravity : public WeakForm::VectorFormVol
  {
  public:
    VectorFormGravity(double rho_g) : WeakForm::VectorFormVol(1), rho_g(rho_g) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) {
      return rho_g * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);;
    }

    // Member.
    double rho_g;
  };
};




#endif
