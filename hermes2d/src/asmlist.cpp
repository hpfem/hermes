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

#include "asmlist.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    AsmList<Scalar>::AsmList(const AsmList<Scalar> & other)
    {
      this->cnt = other.cnt;
      this->cap = other.cap;
      this->idx = (int*) malloc(sizeof(int) * cap);
      this->dof = (int*) malloc(sizeof(int) * cap);
      this->coef = (Scalar*) malloc(sizeof(Scalar) * cap);
      for(unsigned int i = 0; i < cnt; i++) {
        coef[i] = other.coef[i];
        dof[i] = other.dof[i];
        idx[i] = other.idx[i];
      }
    }

    template<typename Scalar>
    AsmList<Scalar>::AsmList()
    {
      cnt = 0;
      cap = 128;
      idx = (int*) malloc(sizeof(int) * cap);
      dof = (int*) malloc(sizeof(int) * cap);
      coef = (Scalar*) malloc(sizeof(Scalar) * cap);
    }

    template<typename Scalar>
    AsmList<Scalar>::~AsmList()
    {
      free(idx);
      free(dof);
      free(coef);
    }

    template<typename Scalar>
    int* AsmList<Scalar>::get_idx()
    {
      return this->idx;
    }

    template<typename Scalar>
    int* AsmList<Scalar>::get_dof()
    {
      return this->dof;
    }

    template<typename Scalar>
    unsigned int AsmList<Scalar>::get_cnt()
    {
      return this->cnt;
    }

    template<typename Scalar>
    Scalar* AsmList<Scalar>::get_coef()
    {
      return this->coef;
    }

    template<typename Scalar>
    void AsmList<Scalar>::add_triplet(int i, int d, Scalar c)
    {
      if(cnt >= cap)
        enlarge();
      idx[cnt] = i;
      dof[cnt] = d;
      coef[cnt++] = c;
    }

    template<typename Scalar>
    void AsmList<Scalar>::enlarge()
    {
      cap = !cap ? 128 : cap * 2;
      idx = (int*) realloc(idx, sizeof(int) * cap);
      dof = (int*) realloc(dof, sizeof(int) * cap);
      coef = (Scalar*) realloc(coef, sizeof(Scalar) * cap);
    }

    template HERMES_API class AsmList<double>;
    template HERMES_API class AsmList<std::complex<double> >;
  }
}