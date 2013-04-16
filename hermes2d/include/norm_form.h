// This file is part of Hermes2D
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
/*! \file norm_form.h
\brief Containc class that calculates the norm.
*/
#ifndef __H2D_NORM_FORM_H_
#define __H2D_NORM_FORM_H_

#include "weakform.h"

namespace Hermes
{
  /// Namespace containing definitions specific for Hermes2D.
  namespace Hermes2D
  {
    class NormForm
    {
    public:
      NormForm(int i, int j);

    protected:
      int i, j;
    };

    template<typename Scalar>
    class ContinuousNormForm : public NormForm
    {
    public:
      ContinuousNormForm(int i, int j);

      void set_marker(int marker);

      virtual Scalar value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const = 0;
    protected:
      int marker;
    };

    template<typename Scalar>
    class DiscontinuousNormForm : public NormForm
    {
    public:
      DiscontinuousNormForm(int i, int j);

      virtual Scalar value(int n, double *wt, DiscontinuousFunc<Scalar> *u, DiscontinuousFunc<Scalar> *v, Geom<double> *e) const = 0;
    };

    template<typename Scalar>
    class MatrixNormFormVol : public MatrixFormVol<Scalar>
    {
    public:
      MatrixNormFormVol(int i, int j, ContinuousNormForm<Scalar>* normForm) : MatrixFormVol<Scalar>(i, j), normForm(normForm) {};

      virtual double value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const
      {
        return this->normForm->value(n, wt, u, v);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, Func<Ord> **ext) const
      {
        return u->val[0] * v->val[0];
      }
    
      MatrixFormVol<double>* clone() const { return new MatrixNormFormVol(*this); };

      protected:
        ContinuousNormForm<Scalar>* normForm;
    };

    template <typename Scalar, NormType normType>
    class DefaultContinuousNormForm : public ContinuousNormForm<Scalar>
    {
    public:
      DefaultContinuousNormForm(int i, int j);
      virtual Scalar value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const;
    };

    template <typename Scalar>
    class DefaultContinuousNormForm<Scalar, HERMES_L2_NORM> : public ContinuousNormForm<Scalar>
    {
    public:
      DefaultContinuousNormForm(int i, int j);
      virtual Scalar value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const;
    };

    template <typename Scalar>
    class DefaultContinuousNormForm<Scalar, HERMES_H1_NORM> : public ContinuousNormForm<Scalar>
    {
    public:
      DefaultContinuousNormForm(int i, int j);
      virtual Scalar value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const;
    };

    template <typename Scalar>
    class DefaultContinuousNormForm<Scalar, HERMES_H1_SEMINORM> : public ContinuousNormForm<Scalar>
    {
    public:
      DefaultContinuousNormForm(int i, int j);
      virtual Scalar value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const;
    };

    template <typename Scalar>
    class DefaultContinuousNormForm<Scalar, HERMES_HCURL_NORM> : public ContinuousNormForm<Scalar>
    {
    public:
      DefaultContinuousNormForm(int i, int j);
      virtual Scalar value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const;
    };

    template <typename Scalar>
    class DefaultContinuousNormForm<Scalar, HERMES_HDIV_NORM> : public ContinuousNormForm<Scalar>
    {
    public:
      DefaultContinuousNormForm(int i, int j);
      virtual Scalar value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const;
    };
  }
}
#endif