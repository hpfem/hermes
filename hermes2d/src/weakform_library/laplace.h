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

namespace Laplace {
  namespace VolumetricMatrixForms {
    class DefaultMatrixFormStiffness : public MatrixFormVol<double>
    {
    public:
      DefaultMatrixFormStiffness(int i, int j, double coeff = 1.0, SymFlag sym = HERMES_SYM) 
            : MatrixFormVol<double>(i, j, sym), coeff(coeff) { }
      DefaultMatrixFormStiffness(int i, int j, std::string area, double coeff = 1.0, SymFlag sym = HERMES_SYM) 
            : MatrixFormVol<double>(i, j, sym, area), coeff(coeff) { }

      template<typename real, typename scalar>
      scalar matrix_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, 
                         Func<real> *v, Geom<real> *e, ExtData<scalar> *ext) {
        return coeff * int_grad_u_grad_v<real, scalar>(n, wt, u, v);
      }

      double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                   Func<double> *v, Geom<double> *e, ExtData<double> *ext) {
        return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
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

    class DefaultMatrixFormMass : public MatrixFormVol<double>
    {
    public:
      DefaultMatrixFormMass(int i, int j, double coeff = 1.0, SymFlag sym = HERMES_SYM) 
            : MatrixFormVol<double>(i, j, sym), coeff(coeff) { }
      DefaultMatrixFormMass(int i, int j, std::string area, double coeff = 1.0, SymFlag sym = HERMES_SYM) 
            : MatrixFormVol<double>(i, j, sym, area), coeff(coeff) { }

      template<typename real, typename scalar>
      scalar matrix_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, 
                         Func<real> *v, Geom<real> *e, ExtData<scalar> *ext) {
        return coeff * int_u_v<real, scalar>(n, wt, u, v);
      }

      double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                   Func<double> *v, Geom<double> *e, ExtData<double> *ext) {
        return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
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

    class DefaultMatrixFormAdvection : public MatrixFormVol<double>
    {
    public:
     DefaultMatrixFormAdvection(int i, int j, double coeff1, double coeff2) 
       : MatrixFormVol<double>(i, j, HERMES_NONSYM), coeff1(coeff1), coeff2(coeff2) { }
     DefaultMatrixFormAdvection(int i, int j, std::string area, double coeff1, double coeff2) 
       : MatrixFormVol<double>(i, j, HERMES_NONSYM, area), coeff1(coeff1), coeff2(coeff2) { }

      template<typename real, typename scalar>
      scalar matrix_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, 
                         Func<real> *v, Geom<real> *e, ExtData<scalar> *ext) {
        return   coeff1 * int_dudx_v<real, scalar>(n, wt, u, v)
               + coeff2 * int_dudy_v<real, scalar>(n, wt, u, v);
      }

      double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                   Func<double> *v, Geom<double> *e, ExtData<double> *ext) {
        return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
      }

      Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
              Geom<Ord> *e, ExtData<Ord> *ext) {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      private:
      double coeff1, coeff2;
    };
  }

  namespace RightHandSides {
    // Generic class for non-constant right-hand side. 
    class DefaultNonConstRightHandSide
    {
    public:
      DefaultNonConstRightHandSide() { };

      virtual double value(double x, double y) const = 0;
      virtual Ord ord(Ord x, Ord y) const = 0;
    };
  }

  namespace VolumetricVectorForms {
      /* Default volumetric vector form \int_{area} coeff v d\bfx 
         coeff... constant number
      */

      class DefaultVectorFormConst : public VectorFormVol<double>
      {
      public:
        DefaultVectorFormConst(int i, double coeff) 
                     : VectorFormVol<double>(i), coeff(coeff) { }
        DefaultVectorFormConst(int i, std::string area, double coeff) 
                     : VectorFormVol<double>(i, area), coeff(coeff) { }

        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                             Geom<double> *e, ExtData<double> *ext) {
          return coeff * int_v<double, double>(n, wt, v);
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
      class DefaultVectorFormNonConst : public VectorFormVol<double>
      {
      public:
        DefaultVectorFormNonConst(int i, RightHandSides::DefaultNonConstRightHandSide* rhs) 
                     : VectorFormVol<double>(i), rhs(rhs) { }
        DefaultVectorFormNonConst(int i, std::string area, RightHandSides::DefaultNonConstRightHandSide* rhs) 
                     : VectorFormVol<double>(i, area), rhs(rhs) { }

        double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                     Geom<double> *e, ExtData<double> *ext) {
          double result = 0;
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
        RightHandSides::DefaultNonConstRightHandSide* rhs;
      };
  }

  namespace SurfaceMatrixForms {
    /* Default surface matrix form \int_{area} coeff u v dS
       coeff... constant number
    */

    class DefaultSurfaceMatrixForm : public MatrixFormSurf<double>
    {
    public:
      DefaultSurfaceMatrixForm(int i, int j, double coeff) 
            : MatrixFormSurf<double>(i, j), coeff(coeff) { }
      DefaultSurfaceMatrixForm(int i, int j, std::string area, double coeff) 
            : MatrixFormSurf<double>(i, j, area), coeff(coeff) { }

      template<typename real, typename scalar>
      scalar matrix_form_surf(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, 
                              Func<real> *v, Geom<real> *e, ExtData<scalar> *ext) {
        return coeff * int_u_v<real, scalar>(n, wt, u, v);
      }

      double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                   Geom<double> *e, ExtData<double> *ext) {
        return matrix_form_surf<double, double>(n, wt, u_ext, u, v, e, ext);
      }

      Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
        return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      private:
        double coeff;
    };
  }

  namespace SurfaceVectorForms {
    /* Default surface vector form \int_{area} coeff v dS
       coeff... constant number
    */

    class DefaultSurfaceVectorForm : public VectorFormSurf<double>
    {
    public:
      DefaultSurfaceVectorForm(int i, double coeff) 
             : VectorFormSurf<double>(i), coeff(coeff) { }
      DefaultSurfaceVectorForm(int i, std::string area, double coeff) 
             : VectorFormSurf<double>(i, area), coeff(coeff) { }

      template<typename real, typename scalar>
      scalar vector_form_surf(int n, double *wt, Func<scalar> *u_ext[], 
                              Func<real> *v, Geom<real> *e, ExtData<scalar> *ext) {
        return coeff * int_v<real, scalar>(n, wt, v);
      }

      double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                   Geom<double> *e, ExtData<double> *ext) {
        return vector_form_surf<double, double>(n, wt, u_ext, v, e, ext);
      }

      Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
        return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }

    private:
      double coeff;
    };
  }

  namespace WeakForms {
    /* Default weak form for the Laplace equation -Laplace u = 0
    */

    class DefaultWeakFormLaplace : public WeakForm<double>
    {
    public:
      DefaultWeakFormLaplace() : WeakForm<double>(1)
      {
        add_matrix_form(new VolumetricMatrixForms::DefaultMatrixFormStiffness(0, 0));
      };
    };
  }
}
#endif
