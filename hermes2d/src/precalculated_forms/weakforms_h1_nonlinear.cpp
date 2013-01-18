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

#include "weakforms_h1_nonlinear.h"
#include "api2d.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace NonlinearElemwiseWeakFormsH1
    {
      template<typename Scalar>
      NonlinearElemwiseDiffusionDerivativeForm<Scalar>::NonlinearElemwiseDiffusionDerivativeForm(int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->previous_iteration_space_index = j;
      }

      template<typename Scalar>
      NonlinearElemwiseDiffusionDerivativeForm<Scalar>::NonlinearElemwiseDiffusionDerivativeForm(int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->previous_iteration_space_index = j;
      }

      template<typename Scalar>
      NonlinearElemwiseDiffusionDerivativeForm<Scalar>::~NonlinearElemwiseDiffusionDerivativeForm()
      {
      };

      template<typename Scalar>
      Scalar NonlinearElemwiseDiffusionDerivativeForm<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * (u_ext[this->previous_iteration_space_index]->dx[i] * v->dx[i] + u_ext[this->previous_iteration_space_index]->dy[i] * v->dy[i]);
        return result;
      }

      template<typename Scalar>
      Ord NonlinearElemwiseDiffusionDerivativeForm<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * (u_ext[this->previous_iteration_space_index]->dx[i] * v->dx[i] + u_ext[this->previous_iteration_space_index]->dy[i] * v->dy[i]);
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* NonlinearElemwiseDiffusionDerivativeForm<Scalar>::clone() const
      {
        return new NonlinearElemwiseDiffusionDerivativeForm<Scalar>(*this);
      }

      template<typename Scalar>
      NonlinearElemwiseDiffusionValueForm<Scalar>::NonlinearElemwiseDiffusionValueForm
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
      }

      template<typename Scalar>
      NonlinearElemwiseDiffusionValueForm<Scalar>::NonlinearElemwiseDiffusionValueForm
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
      }

      template<typename Scalar>
      NonlinearElemwiseDiffusionValueForm<Scalar>::~NonlinearElemwiseDiffusionValueForm()
      {
      };

      template<typename Scalar>
      Scalar NonlinearElemwiseDiffusionValueForm<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
            result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        return result;
      }

      template<typename Scalar>
      Ord NonlinearElemwiseDiffusionValueForm<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* NonlinearElemwiseDiffusionValueForm<Scalar>::clone() const
      {
        return new NonlinearElemwiseDiffusionValueForm<Scalar>(*this);
      }

      template class HERMES_API NonlinearElemwiseDiffusionDerivativeForm<double>;
      template class HERMES_API NonlinearElemwiseDiffusionDerivativeForm<std::complex<double> >;
      template class HERMES_API NonlinearElemwiseDiffusionValueForm<double>;
      template class HERMES_API NonlinearElemwiseDiffusionValueForm<std::complex<double> >;
    }
  }
}
