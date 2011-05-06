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

#ifndef __H2D_HCURL_WEAK_FORMS_H
#define __H2D_HCURL_WEAK_FORMS_H

#include "../integrals/hcurl.h"


namespace WeakFormsHcurl {

  /* Default volumetric matrix form \int_{area} coeff \curl E \curl F d\bfx
    coeff... constant number
  */

  class DefaultLinearCurlCurl : public WeakForm::MatrixFormVol
  {
  public:
    DefaultLinearCurlCurl(int i, int j, scalar coeff = 1.0, SymFlag sym = HERMES_SYM);
    DefaultLinearCurlCurl(int i, int j, std::string area, scalar coeff = 1.0, SymFlag sym = HERMES_SYM);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                          Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

  private:
      scalar coeff;
  };

  /* Default volumetric matrix form \int_{area} coeff E \cdot F d\bfx
      coeff... constant number
  */

  class DefaultLinearMass : public WeakForm::MatrixFormVol
  {
  public:
    DefaultLinearMass(int i, int j, scalar coeff = 1.0, SymFlag sym = HERMES_SYM);
    DefaultLinearMass(int i, int j, std::string area, scalar coeff = 1.0, SymFlag sym = HERMES_SYM);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                          Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

  private:
      scalar coeff;
  };

  /* Default surface matrix form \int_{area} coeff e tau f tau dS
      coeff... constant number
  */

  class DefaultMatrixFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    DefaultMatrixFormSurf(int i, int j, scalar coeff);
    DefaultMatrixFormSurf(int i, int j, std::string area, scalar coeff);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                          Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormSurf* clone();

  private:
      scalar coeff;
  };

  /* Default volumetric vector form \int_{area} (coeff0, coeff1) \cdot E d\bfx
      coeff... constant number
  */

  class DefaultVectorFormConst : public WeakForm::VectorFormVol
  {
  public:
    DefaultVectorFormConst(int i, scalar coeff0, scalar coeff1);
    DefaultVectorFormConst(int i, std::string area, scalar coeff0, scalar coeff1);

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                          Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

  private:
    scalar coeff0, coeff1;
  };
};
#endif
