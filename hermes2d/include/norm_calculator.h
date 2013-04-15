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
/*! \file norm_calculator.h
\brief Containc class that calculates the norm.
*/
#ifndef __H2D_DISTANCE_FUNCTION_CALCULATOR_H_
#define __H2D_DISTANCE_FUNCTION_CALCULATOR_H_

#include "weakform.h"
#include "global.h"

namespace Hermes
{
  /// Namespace containing definitions specific for Hermes2D.
  namespace Hermes2D
  {
    class HERMES_API NormCalculator : public Hermes::Mixins::Loggable
    {
    public:
      template<typename FormType, typename Scalar, NormType normType>
      class NormForm : public FormType
      {
      public:
        NormForm(int i, int j, Hermes::vector<std::string> areas) : FormType(i, j, areas) {};
        NormForm(int i, int j, std::string area) : FormType(i, j, area) {};

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[],
          Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
          Func<Scalar> **ext) const = 0;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
          Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
          Func<Ord> **ext) const = 0;

        virtual FormType* clone() const { return this; };
      };

      template <typename FormType, typename Scalar>
      class NormForm<FormType, Func<Scalar>, HERMES_L2_NORM>
      {
      public:
        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[],
          Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
          Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
          Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
          Func<Ord> **ext) const;
      };

      template <typename FormType, typename Scalar>
      class NormForm<FormType, Func<Scalar>, HERMES_H1_NORM>
      {
      public:
        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[],
          Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
          Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
          Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
          Func<Ord> **ext) const;
      };

      template <typename FormType, typename Scalar>
      class NormForm<FormType, Func<Scalar>, HERMES_H1_SEMINORM>
      {
      public:
        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[],
          Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
          Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
          Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
          Func<Ord> **ext) const;
      };

      template <typename FormType, typename Scalar>
      class NormForm<FormType, Func<Scalar>, HERMES_HCURL_NORM>
      {
      public:
        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[],
          Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
          Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
          Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
          Func<Ord> **ext) const;
      };

      template <typename FormType, typename Scalar>
      class NormForm<FormType, Func<Scalar>, HERMES_HDIV_NORM>
      {
      public:
        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[],
          Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
          Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
          Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
          Func<Ord> **ext) const;
      };
    };
  }
}
#endif