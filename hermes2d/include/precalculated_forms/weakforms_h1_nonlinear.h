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

#ifndef __H2D_NONLINEAR_H1_WEAK_FORMS_H
#define __H2D_NONLINEAR_H1_WEAK_FORMS_H

#include "../integrals/h1.h"
#include "../weakform/weakform.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace NonlinearElemwiseWeakFormsH1
    {
      template<typename Scalar>
      class HERMES_API NonlinearElemwiseDiffusionDerivativeForm : public MatrixFormVol<Scalar>
      {
      public:
        NonlinearElemwiseDiffusionDerivativeForm(int i, int j, std::string area = HERMES_ANY);

        NonlinearElemwiseDiffusionDerivativeForm(int i, int j, Hermes::vector<std::string> areas);

        ~NonlinearElemwiseDiffusionDerivativeForm();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      };

      template<typename Scalar>
      class HERMES_API NonlinearElemwiseDiffusionValueForm : public MatrixFormVol<Scalar>
      {
      public:
        NonlinearElemwiseDiffusionValueForm(int i, int j, std::string area = HERMES_ANY);

        NonlinearElemwiseDiffusionValueForm(int i, int j, Hermes::vector<std::string> areas);

        void init_tables();

        ~NonlinearElemwiseDiffusionValueForm();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      };
    }
  }
}
#endif
  
