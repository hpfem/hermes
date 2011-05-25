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

#ifndef __H2D_H1_WEAK_FORMS_H
#define __H2D_H1_WEAK_FORMS_H

#include "../integrals/h1.h"

namespace WeakFormsH1 
{
  /* Default volumetric matrix form \int_{area} const_coeff * function_coeff(x, y) * u * v \bfx
     const_coeff... constant number
     function_coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    DefaultMatrixFormVol(int i, int j, std::string area = HERMES_ANY,
                         scalar const_coeff = 1.0, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                         SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    DefaultMatrixFormVol(int i, int j, Hermes::vector<std::string> areas,
                         scalar const_coeff = 1.0, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                         SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    ~DefaultMatrixFormVol();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    private:
      scalar const_coeff;
      DefaultFunction* function_coeff;
      GeomType gt;
  };

  /* Default volumetric matrix form \int_{area} const_coeff * spline_coeff'(u_ext[0]) u \nabla u_ext[0] \cdot \nabla v
     + const_coeff * spline_coeff(u_ext[0]) * \nabla u \cdot \nabla v d\bfx
     const_coeff... constant number
     spline_coeff... nonconstant parameter given by cubic spline
  */

  class HERMES_API DefaultJacobianDiffusion : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianDiffusion(int i, int j, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                             CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                             SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    DefaultJacobianDiffusion(int i, int j, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                             CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                             SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    ~DefaultJacobianDiffusion();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    private:
      int idx_j;
      scalar const_coeff;
      CubicSpline* spline_coeff;
      GeomType gt;
  };

  /* Default volumetric matrix form
     \int_{area} spline_coeff1`(u_ext[0]) * u * u_ext[0]->dx * v
     + spline_coeff1(u_ext[0]) * u->dx * v
     + spline_coeff2`(u_ext[0]) * u * u_ext[0]->dy * v
     + spline_coeff2(u_ext[0]) * u->dy * v d\bfx.
     spline_coeff1, spline_coeff2... non-constant parameters given by cubic splines
  */

  class HERMES_API DefaultJacobianAdvection : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianAdvection(int i, int j, std::string area = HERMES_ANY, 
                             scalar const_coeff1 = 1.0, scalar const_coeff2 = 1.0,
                             CubicSpline* c_spline1 = HERMES_DEFAULT_SPLINE,
                             CubicSpline* c_spline2 = HERMES_DEFAULT_SPLINE, GeomType gt = HERMES_PLANAR);

   DefaultJacobianAdvection(int i, int j, Hermes::vector<std::string> areas, 
                            scalar const_coeff1 = 1.0, scalar const_coeff2 = 1.0,
                            CubicSpline* c_spline1 = HERMES_DEFAULT_SPLINE,
                            CubicSpline* c_spline2 = HERMES_DEFAULT_SPLINE,
                            GeomType gt = HERMES_PLANAR);

   ~DefaultJacobianAdvection();

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    private:
      int idx_j;
      scalar const_coeff1, const_coeff2;
      CubicSpline* spline_coeff1, *spline_coeff2;
      GeomType gt;
  };

  /* Default volumetric vector form \int_{area} const_coeff * function_coeff(x, y) * v d\bfx
     const_coeff... constant number
     function_coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    DefaultVectorFormVol(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                         DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                         GeomType gt = HERMES_PLANAR);

    DefaultVectorFormVol(int i, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                         DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                         GeomType gt = HERMES_PLANAR);

    ~DefaultVectorFormVol();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    private:
      scalar const_coeff;
      DefaultFunction* function_coeff;
      GeomType gt;
  };

  /* Default volumetric vector form \int_{area} const_coeff * function_coeff(x, y) * u_ext[0] * v d\bfx
     const_coeff... constant number
     function_coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultResidualVol : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualVol(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                       DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                       GeomType gt = HERMES_PLANAR);
    DefaultResidualVol(int i, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                       DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                       GeomType gt = HERMES_PLANAR);

    ~DefaultResidualVol();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    private:
      int idx_i;
      scalar const_coeff;
      DefaultFunction* function_coeff;
      GeomType gt;
  };

  /* Default volumetric vector form \int_{area} const_coeff * spline_coeff(u_ext[0]) *
     \nabla u_ext[0] \cdot \nabla v d\bfx
     const_coeff... constant number
     spline_coeff... non-constant parameter given by a cubic spline
  */

  class HERMES_API DefaultResidualDiffusion : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualDiffusion(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                             CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                             GeomType gt = HERMES_PLANAR);

    DefaultResidualDiffusion(int i, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                             CubicSpline* c_spline = HERMES_DEFAULT_SPLINE, 
                             GeomType gt = HERMES_PLANAR);

    ~DefaultResidualDiffusion();

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    private:
      int idx_i;
      scalar const_coeff;
      CubicSpline* spline_coeff;
      GeomType gt;
  };

  /* Default volumetric vector form \int_{area} spline_coeff1(u_ext[0]) * u->dx * v->val
     + spline_coeff2(u_ext[0]) * u->dy * v->val d\bfx
     spline_coeff1, spline_coeff2... non-constant parameters given by cubic splines
  */

  class HERMES_API DefaultResidualAdvection : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualAdvection(int i, std::string area = HERMES_ANY, 
                             scalar const_coeff1 = 1.0, scalar const_coeff2 = 1.0, 
                             CubicSpline* c_spline1 = HERMES_DEFAULT_SPLINE,
                             CubicSpline* c_spline2 = HERMES_DEFAULT_SPLINE,
                             GeomType gt = HERMES_PLANAR);
    DefaultResidualAdvection(int i, Hermes::vector<std::string> areas,\
                             scalar const_coeff1 = 1.0, scalar const_coeff2 = 1.0,
                             CubicSpline* c_spline1 = HERMES_DEFAULT_SPLINE,
                             CubicSpline* c_spline2 = HERMES_DEFAULT_SPLINE, GeomType gt = HERMES_PLANAR);

    ~DefaultResidualAdvection();

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    private:
      int idx_i;
      scalar const_coeff1, const_coeff2;
      CubicSpline* spline_coeff1, *spline_coeff2;
      GeomType gt;
  };

  /* Default surface matrix form \int_{area} const_coeff * function_coeff(x, y) * u * v dS
     const_coeff... constant number
     function_coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultMatrixFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    DefaultMatrixFormSurf(int i, int j, std::string area = HERMES_ANY,
                          scalar const_coeff = 1.0, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                          GeomType gt = HERMES_PLANAR);

    DefaultMatrixFormSurf(int i, int j, Hermes::vector<std::string> areas,
                          scalar const_coeff = 1.0, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                          GeomType gt = HERMES_PLANAR);

    ~DefaultMatrixFormSurf();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;
      
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormSurf* clone();

    private:
      scalar const_coeff;
      DefaultFunction* function_coeff;
      GeomType gt;
  };

  /* Default surface matrix form \int_{area} const_coeff * spline_coeff'(u_ext[0]) * u_ext[0] * u * v
     + const_coeff * spline_coeff(u_ext[0]) * u * v dS
     spline_coeff... non-constant parameter given by a spline
  */

  class HERMES_API DefaultJacobianFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    DefaultJacobianFormSurf(int i, int j, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                            CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                            GeomType gt = HERMES_PLANAR);
    DefaultJacobianFormSurf(int i, int j, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                            CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                            GeomType gt = HERMES_PLANAR);

    ~DefaultJacobianFormSurf();
      
    template<typename Real, typename Scalar>
    Scalar matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormSurf* clone();

    private:
      int idx_j;
      scalar const_coeff;
      CubicSpline* spline_coeff;
      GeomType gt;
  };

  /* Default surface vector form \int_{area} const_coeff * function_coeff(x, y) * v dS
     const_coeff... constant number
     function_coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultVectorFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    DefaultVectorFormSurf(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                          DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                          GeomType gt = HERMES_PLANAR);
    DefaultVectorFormSurf(int i, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                          DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                          GeomType gt = HERMES_PLANAR);

    ~DefaultVectorFormSurf();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormSurf* clone();

    private:
      scalar const_coeff;
      DefaultFunction* function_coeff;
      GeomType gt;
  };

  class HERMES_API DefaultMultiComponentVectorFormSurf : public WeakForm::MultiComponentVectorFormSurf
  {
  public:
    DefaultMultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates,
                                        std::string area = HERMES_ANY,
                                        Hermes::vector<scalar> coeffs = Hermes::vector<scalar>(1.0),
                                        GeomType gt = HERMES_PLANAR);
    DefaultMultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates,
                                        Hermes::vector<std::string> areas,
                                        Hermes::vector<scalar> coeffs, GeomType gt = HERMES_PLANAR);

    virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    /* FIXME
    virtual WeakForm::VectorFormSurf* clone() {
      return new DefaultMultiComponentVectorFormSurf(*this);
    }
    */

    private:
      Hermes::vector<scalar> coeffs;
      GeomType gt;
  };

  /* Default surface vector form \int_{area} const_coeff * function_coeff(x, y) * u_ext[0] v dS
     const_coeff... constant number
     function_coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultResidualSurf : public WeakForm::VectorFormSurf
  {
  public:
    DefaultResidualSurf(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                        DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                        GeomType gt = HERMES_PLANAR);
    DefaultResidualSurf(int i, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                        DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                        GeomType gt = HERMES_PLANAR);

    ~DefaultResidualSurf();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormSurf* clone();

    private:
      int idx_i;
      scalar const_coeff;
      DefaultFunction* function_coeff;
      GeomType gt;
  };

  /* Default weak form for the Laplace equation -div(const_coeff spline_coeff(u) grad u) = 0. */

  class HERMES_API DefaultWeakFormLaplace : public WeakForm
  {
  public:
    DefaultWeakFormLaplace(std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                           CubicSpline* spline_coeff = HERMES_DEFAULT_SPLINE,
                           GeomType gt = HERMES_PLANAR);
  };


  /* Default weak form for the Poisson equation -div(const_coeff spline_coeff(u) grad u) - rhs = 0. */

  class HERMES_API DefaultWeakFormPoisson : public WeakForm
  {
  public:
  DefaultWeakFormPoisson(DefaultFunction* rhs = HERMES_DEFAULT_FUNCTION,
                         std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                         CubicSpline* spline_coeff = HERMES_DEFAULT_SPLINE,
                         GeomType gt = HERMES_PLANAR);
  };
};

#endif
