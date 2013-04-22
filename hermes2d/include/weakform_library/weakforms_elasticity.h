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

#include "../weakform/weakform.h"
#include "../spline.h"

/* Default weak form for linear elasticity (Lame equations)
with Dirichlet and/or zero Neumann BC

Nonzero Neumann and Newton boundary conditions can be enabled
by creating a descendant and adding surface forms to it.
*/
namespace Hermes
{
  namespace Hermes2D
  {
    namespace WeakFormsElasticity {
      /* Single-component version -- to be used for multimesh assembling */

      template<typename Scalar>
      class HERMES_API DefaultJacobianElasticity_0_0 : public MatrixFormVol<Scalar>
      {
      public:
        DefaultJacobianElasticity_0_0(unsigned int i, unsigned int j, double lambda, double mu);
        DefaultJacobianElasticity_0_0(unsigned int i, unsigned int j, std::string area, double lambda, double mu);

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
          Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      private:
        double lambda, mu;
      };

      template<typename Scalar>
      class HERMES_API DefaultJacobianElasticity_0_1 : public MatrixFormVol<Scalar>
      {
      public:
        DefaultJacobianElasticity_0_1(unsigned int i, unsigned int j, double lambda, double mu);
        DefaultJacobianElasticity_0_1(unsigned int i, unsigned int j, std::string area, double lambda, double mu);

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
          Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      private:
        double lambda, mu;
      };

      template<typename Scalar>
      class HERMES_API DefaultResidualElasticity_0_0 : public VectorFormVol<Scalar>
      {
      public:
        DefaultResidualElasticity_0_0(unsigned int i, double lambda, double mu);
        DefaultResidualElasticity_0_0(unsigned int i, std::string area, double lambda, double mu);

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormVol<Scalar>* clone() const;
      private:
        double lambda, mu;
      };

      template<typename Scalar>
      class HERMES_API DefaultResidualElasticity_0_1 : public VectorFormVol<Scalar>
      {
      public:
        DefaultResidualElasticity_0_1(unsigned int i, double lambda, double mu);
        DefaultResidualElasticity_0_1(unsigned int i, std::string area, double lambda, double mu);

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormVol<Scalar>* clone() const;
      private:
        double lambda, mu;
      };

      template<typename Scalar>
      class HERMES_API DefaultResidualElasticity_1_0 : public VectorFormVol<Scalar>
      {
      public:
        DefaultResidualElasticity_1_0(unsigned int i, double lambda, double mu);
        DefaultResidualElasticity_1_0(unsigned int i, std::string area, double lambda, double mu);

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormVol<Scalar>* clone() const;
      private:
        double lambda, mu;
      };

      template<typename Scalar>
      class HERMES_API DefaultResidualElasticity_1_1 : public VectorFormVol<Scalar>
      {
      public:
        DefaultResidualElasticity_1_1(unsigned int i, double lambda, double mu);
        DefaultResidualElasticity_1_1(unsigned int i, std::string area, double lambda, double mu);

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormVol<Scalar>* clone() const;
      private:
        double lambda, mu;
      };

      template<typename Scalar>
      class HERMES_API DefaultJacobianElasticity_1_1 : public MatrixFormVol<Scalar>
      {
      public:
        DefaultJacobianElasticity_1_1(unsigned int i, unsigned int j, double lambda, double mu);
        DefaultJacobianElasticity_1_1(unsigned int i, unsigned int j, std::string area, double lambda, double mu);

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
          Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      private:
        double lambda, mu;
      };
    };
  }
}
#endif