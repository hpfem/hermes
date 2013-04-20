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

#include "weakform/weakform.h"

namespace Hermes
{
  /// Namespace containing definitions specific for Hermes2D.
  namespace Hermes2D
  {
    class HERMES_API NormForm
    {
    public:
      NormForm(int i, int j);

      /// set area
      /// \todo use this - it is not used so far.
      void set_area(std::string area);

      /// Coordinates.
      int i, j;

    protected:
      std::string area;
    };

    template<typename Scalar>
    class HERMES_API NormFormVol : public NormForm
    {
    public:
      NormFormVol(int i, int j);

      virtual Scalar value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const = 0;
    };

    template<typename Scalar>
    class HERMES_API NormFormSurf : public NormForm
    {
    public:
      NormFormSurf(int i, int j);

      virtual Scalar value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const = 0;
    };

    template<typename Scalar>
    class HERMES_API NormFormDG : public NormForm
    {
    public:
      NormFormDG(int i, int j);

      virtual Scalar value(int n, double *wt, DiscontinuousFunc<Scalar> *u, DiscontinuousFunc<Scalar> *v, Geom<double> *e) const = 0;
    };

    template <typename Scalar>
    class HERMES_API DefaultNormFormVol : public NormFormVol<Scalar>
    {
    public:
      DefaultNormFormVol(int i, int j, NormType normType);
      
      Scalar value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const;

    protected:
      NormType normType;
    };

    template<typename Scalar>
    class HERMES_API MatrixDefaultNormFormVol : public MatrixFormVol<Scalar>
    {
    public:
      MatrixDefaultNormFormVol(int i, int j, NormType normType);

      Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const;

      Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, Func<Ord> **ext) const;
    
      MatrixFormVol<Scalar>* clone() const;

      protected:
      NormType normType;
    };

    template<typename Scalar>
    class HERMES_API VectorDefaultNormFormVol : public VectorFormVol<Scalar>
    {
    public:
      VectorDefaultNormFormVol(int i, NormType normType);

      Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const;

      Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;
    
      VectorFormVol<Scalar>* clone() const;

      protected:
      NormType normType;
    };
  }
}
#endif