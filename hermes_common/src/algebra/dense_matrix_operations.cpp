// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file dense_matrix_operations.cpp
\brief Dense (small) simply stored matrix operations.
*/
#include "dense_matrix_operations.h"
#include "exceptions.h"

namespace Hermes
{
  namespace Algebra
  {
    namespace DenseMatrixOperations
    {
      template<typename T>
      void ludcmp(T **a, int n, int *indx, double *d)
      {
        int i, imax = 0, j, k;
        T big, dum, sum, temp;
        T *vv = new T[n];

        *d = 1.0;
        for (i = 0; i < n; i++)
        {
          big = 0.0;
          for (j = 0; j < n; j++)
          {
            temp = a[i][j];
            if(std::abs(temp) > std::abs(big))
              big = temp;
          }
          if(big == 0.0)
          {
            delete [] vv;
            throw Exceptions::Exception("Singular matrix in routine LUDCMP!");
          }
          vv[i] = 1.0 / big;
        }
        for (j = 0; j < n; j++)
        {
          for (i = 0; i < j; i++)
          {
            sum = a[i][j];
            for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
          }
          big = 0.0;
          for (i = j; i < n; i++)
          {
            sum = a[i][j];
            for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            dum = vv[i]*std::abs(sum);
            if(std::abs(dum) >= std::abs(big))
            {
              big = dum;
              imax = i;
            }
          }
          if(j != imax)
          {
            for (k = 0; k < n; k++)
            {
              dum = a[imax][k];
              a[imax][k] = a[j][k];
              a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
          }
          indx[j] = imax;
          if(a[j][j] == 0.0) a[j][j] = 1.0e-20;
          if(j != n-1)
          {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i < n; i++) a[i][j] *= dum;
          }
        }
        delete [] vv;
      }

      template<typename T, typename S>
      void lubksb(T **a, int n, int *indx, S *b)
      {
        int i, ip, j;
        S sum;

        for (i = 0; i < n; i++)
        {
          ip = indx[i];
          sum = b[ip];
          b[ip] = b[i];
          for (j = 0; j < i; j++)
            sum = sum - a[i][j]*b[j];
          b[i] = sum;
        }
        for (i = n-1; i >= 0; i--)
        {
          sum = b[i];
          for (j = i + 1; j < n; j++)
            sum -= a[i][j]*b[j];
          b[i] = sum / a[i][i];
        }
      }

      template<typename T>
      void choldc(T **a, int n, T p[])
      {
        int i, j, k;
        for (i = 0; i < n; i++)
        {
          for (j = i; j < n; j++)
          {
            T sum = a[i][j];
            k = i;
            while (--k >= 0) sum -= a[i][k] * a[j][k];
            if(i == j)
            {
              if((std::complex<double>(sum)).real() <= 0.0)
                throw Exceptions::Exception("CHOLDC failed!");
              else p[i] = sqrt(sum);
            }
            else a[j][i] = sum / p[i];
          }
        }
      }

      template HERMES_API void ludcmp<double>(double **a, int n, int *indx, double *d);
      template HERMES_API void ludcmp<std::complex<double> >(std::complex<double> **a, int n, int *indx, double *d);
      template HERMES_API void lubksb<double, std::complex<double> >(double **a, int n, int *indx, std::complex<double> *d);
      template HERMES_API void lubksb<double, double>(double **a, int n, int *indx, double *d);
      template HERMES_API void lubksb<std::complex<double>, std::complex<double> >(std::complex<double> **a, int n, int *indx, std::complex<double> *d);
      template HERMES_API void choldc<double>(double **a, int n, double p[]);
      template HERMES_API void choldc<std::complex<double> >(std::complex<double> **a, int n, std::complex<double> p[]);
    }
  }
}
