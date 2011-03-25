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

#ifndef __H2D_LAPLACE_WEAK_FORMS_H
#define __H2D_LAPLACE_WEAK_FORMS_H

#include "../integrals/integrals_h1.h"

/* Default volumetric matrix form \int_{area} coeff \nabla u \cdot \nabla v d\bfx 
   coeff... constant number
*/

class DefaultMatrixFormStiffness : public WeakForm::MatrixFormVol
{
public:
  DefaultMatrixFormStiffness(int i, int j, double coeff = 1.0) 
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM), coeff(coeff) { }
  DefaultMatrixFormStiffness(int i, int j, std::string area, double coeff = 1.0) 
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM, area), coeff(coeff) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    return coeff * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
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
    double coeff;
};

/* Default volumetric matrix form \int_{area} coeff u v d\bfx 
   coeff... constant number
*/

class DefaultMatrixFormMass : public WeakForm::MatrixFormVol
{
public:
  DefaultMatrixFormMass(int i, int j, double coeff = 1.0) 
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM), coeff(coeff) { }
  DefaultMatrixFormMass(int i, int j, std::string area, double coeff = 1.0) 
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM, area), coeff(coeff) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    return coeff * int_u_v<Real, Scalar>(n, wt, u, v);
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
    double coeff;
};

/* Default volumetric matrix form \int_{area} (coeff1, coeff2) \cdot \nabla u vd\bfx 
   coeff1, coeff2... constant number
*/

class DefaultMatrixFormAdvection : public WeakForm::MatrixFormVol
{
public:
 DefaultMatrixFormAdvection(int i, int j, double coeff1, double coeff2) 
   : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM), coeff1(coeff1), coeff2(coeff2) { }
 DefaultMatrixFormAdvection(int i, int j, std::string area, double coeff1, double coeff2) 
   : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM, area), coeff1(coeff1), coeff2(coeff2) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    return   coeff1 * int_dudx_v<Real, Scalar>(n, wt, u, v)
           + coeff2 * int_dudy_v<Real, Scalar>(n, wt, u, v);
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
  double coeff1, coeff2;
};

/* Default volumetric vector form \int_{area} coeff v d\bfx 
   coeff... constant number
*/

class DefaultVectorFormVolConst : public WeakForm::VectorFormVol
{
public:
  DefaultVectorFormVolConst(int i, double coeff) 
               : WeakForm::VectorFormVol(i), coeff(coeff) { }
  DefaultVectorFormVolConst(int i, std::string area, double coeff) 
               : WeakForm::VectorFormVol(i, area), coeff(coeff) { }

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<scalar> *ext) {
    return coeff * int_v<scalar, scalar>(n, wt, v);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) {
    return int_v<Ord, Ord>(n, wt, v);
  }

private:
  double coeff;
};

/* Default volumetric vector form \int_{area} rhs(x, y) v d\bfx 
   rhs(x, y)... non-constant right-hand side
*/

// Generic class for non-constant right-hand side. 
class DefaultNonConstRightHandSide
{
public:
  DefaultNonConstRightHandSide() { };

  virtual scalar value(double x, double y) const = 0;
  virtual Ord ord(Ord x, Ord y) const = 0;
};

class DefaultVectorFormVolNonConst : public WeakForm::VectorFormVol
{
public:
  DefaultVectorFormVolNonConst(int i, DefaultNonConstRightHandSide* rhs) 
               : WeakForm::VectorFormVol(i), rhs(rhs) { }
  DefaultVectorFormVolNonConst(int i, std::string area, DefaultNonConstRightHandSide* rhs) 
               : WeakForm::VectorFormVol(i, area), rhs(rhs) { }

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
               Geom<double> *e, ExtData<scalar> *ext) {
    scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (rhs->value(e->x[i], e->y[i]) * v->val[i]);
    return result;
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) {
    Ord result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (rhs->ord(e->x[i], e->y[i]) * v->val[i]);
    return result;
  }

private:
  DefaultNonConstRightHandSide* rhs;
};

/* Default surface matrix form \int_{area} coeff u v dS
   coeff... constant number
*/

class DefaultMatrixFormSurf : public WeakForm::MatrixFormSurf
{
public:
  DefaultMatrixFormSurf(int i, int j, double coeff) 
        : WeakForm::MatrixFormSurf(i, j), coeff(coeff) { }
  DefaultMatrixFormSurf(int i, int j, std::string area, double coeff) 
        : WeakForm::MatrixFormSurf(i, j, area), coeff(coeff) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                          Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    return coeff * int_u_v<Real, Scalar>(n, wt, u, v);
  }

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
               Geom<double> *e, ExtData<scalar> *ext) {
    return matrix_form_surf<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
    return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  private:
    double coeff;
};


/* Default surface vector form \int_{area} coeff v dS
   coeff... constant number
*/

class DefaultVectorFormSurf : public WeakForm::VectorFormSurf
{
public:
  DefaultVectorFormSurf(int i, double coeff) 
         : WeakForm::VectorFormSurf(i), coeff(coeff) { }
  DefaultVectorFormSurf(int i, std::string area, double coeff) 
         : WeakForm::VectorFormSurf(i, area), coeff(coeff) { }

  template<typename Real, typename Scalar>
  Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], 
                          Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    return coeff * int_v<Real, Scalar>(n, wt, v);
  }

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
               Geom<double> *e, ExtData<scalar> *ext) {
    return vector_form_surf<scalar, scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
    return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }

private:
  double coeff;
};

/* Default weak form for the Laplace equation -Laplace u = 0
*/

class DefaultWeakFormLaplace : public WeakForm
{
public:
  DefaultWeakFormLaplace() : WeakForm(1)
  {
    add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
  };
};

/* Default non-constant Dirichlet boundary condition 
   based on an exact solution
*/

class DefaultEssentialBCNonConst : public EssentialBC
{
public:
  DefaultEssentialBCNonConst(Hermes::vector<std::string> markers_, 
                             ExactSolutionScalar* exact_solution) : 
        EssentialBC(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    for (unsigned int i=0; i < markers.size(); i++) markers.push_back(markers_[i]);
  };
  DefaultEssentialBCNonConst(std::string marker, ExactSolutionScalar* exact_solution) : 
        EssentialBC(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    markers.push_back(marker);
  };
  
  ~DefaultEssentialBCNonConst() {};

  virtual EssentialBCValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar function(double x, double y) const {
    return exact_solution->value(x, y);
  };

  ExactSolutionScalar* exact_solution;
};

#endif
