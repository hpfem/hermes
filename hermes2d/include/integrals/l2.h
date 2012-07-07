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

#ifndef __H2D_INTEGRALS_L2_H
#define __H2D_INTEGRALS_L2_H

#include "../quadrature/limit_order.h"
#include "../weakform/weakform.h"
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    class MatrixFormVolL2 : public MatrixFormVol<Scalar>
    {
    public:
      // One area.
      MatrixFormVolL2(int i, int j, std::string area = HERMES_ANY,
        SymFlag sym = HERMES_NONSYM) : MatrixFormVol<Scalar>(i, j, area, sym) {}
      // Multiple areas.
      MatrixFormVolL2(int i, int j, Hermes::vector<std::string> areas,
        SymFlag sym = HERMES_NONSYM) : MatrixFormVol<Scalar>(i, j, areas, sym) {}

      virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * conj(v->val[i]));
        return result;
      }

      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const
      {
        Hermes::Ord result = Hermes::Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * conj(v->val[i]));
        return result;
      }

      MatrixFormVol<Scalar>* clone()
      {
        return new MatrixFormVolL2<Scalar>(*this);
      }
    };
  }
}
#endif