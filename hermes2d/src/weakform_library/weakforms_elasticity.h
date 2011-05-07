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
#ifndef H2D_COMPLEX
namespace WeakFormsElasticity {

  /* Single-component version -- to be used for multimesh assembling */

  class HERMES_API DefaultJacobianElasticity_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianElasticity_0_0(unsigned int i, unsigned int j, double lambda, double mu);
    DefaultJacobianElasticity_0_0(unsigned int i, unsigned int j, std::string area, double lambda, double mu);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
  
  private:
      double lambda, mu;
  };

  class HERMES_API DefaultJacobianElasticity_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianElasticity_0_1(unsigned int i, unsigned int j, double lambda, double mu);
    DefaultJacobianElasticity_0_1(unsigned int i, unsigned int j, std::string area, double lambda, double mu);
    
    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
    
    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
            Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
  
  private:
    double lambda, mu;
  };

  class HERMES_API DefaultResidualElasticity_0_0 : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualElasticity_0_0(unsigned int i, double lambda, double mu);
    DefaultResidualElasticity_0_0(unsigned int i, std::string area, double lambda, double mu);

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double lambda, mu;
  };

  class HERMES_API DefaultResidualElasticity_0_1 : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualElasticity_0_1(unsigned int i, double lambda, double mu);
    DefaultResidualElasticity_0_1(unsigned int i, std::string area, double lambda, double mu);

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double lambda, mu;
  };

  class HERMES_API DefaultResidualElasticity_1_0 : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualElasticity_1_0(unsigned int i, double lambda, double mu);
    DefaultResidualElasticity_1_0(unsigned int i, std::string area, double lambda, double mu);

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double lambda, mu;
  };

  class HERMES_API DefaultResidualElasticity_1_1 : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualElasticity_1_1(unsigned int i, double lambda, double mu);
    DefaultResidualElasticity_1_1(unsigned int i, std::string area, double lambda, double mu);

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double lambda, mu;
  };

  class HERMES_API DefaultJacobianElasticity_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianElasticity_1_1(unsigned int i, unsigned int j, double lambda, double mu);
    DefaultJacobianElasticity_1_1(unsigned int i, unsigned int j, std::string area, double lambda, double mu);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double lambda, mu;
  };

  class HERMES_API DefaultJacobianElasticity_00_11 
    : public WeakForm::MultiComponentMatrixFormVol
  {
  public:
    DefaultJacobianElasticity_00_11
      (Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, double lambda, double mu);
    DefaultJacobianElasticity_00_11
      (Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, std::string area, double lambda, double mu);

    template<typename Real, typename Scalar>
    void matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                     Geom<Real> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const;

    virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                        Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double lambda, mu;
  };

  class HERMES_API DefaultResidualElasticity_00_11 : public WeakForm::MultiComponentVectorFormVol
  {
  public:
    DefaultResidualElasticity_00_11
      (Hermes::vector<unsigned int> coordinates, double lambda, double mu);
    DefaultResidualElasticity_00_11
      (Hermes::vector<unsigned int> coordinates, std::string area, double lambda, double mu);

    template<typename Real, typename Scalar>
    void vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                        Geom<Real> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const;

    virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double lambda, mu;
  };
};
#endif

#endif
