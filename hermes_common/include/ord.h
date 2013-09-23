// This file is part of HermesCommon
//
// Hermes is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes; if not, see <http://www.gnu.prg/licenses/>.
/*! \file ord.h
\brief Contains class Ord for calculation of integration order.
*/
#include <complex>
#include "util/compat.h"

#ifndef __HERMES_COMMON_ORD_H_
#define __HERMES_COMMON_ORD_H_

namespace Hermes
{
  /// Base type for orders of functions.
  ///
  /// We defined a special arithmetics with this type to be able to analyze forms
  /// and determine the necessary integration order.  This works for forms, but it also
  /// works for user-defined functions.
  class HERMES_API Ord
  {
  public:

    Ord();
    explicit Ord(int o);
    explicit Ord(double o);

    int get_order() const;

    static Ord get_max_order();

    Ord operator + (const Ord &o);
    Ord operator + (double d);
    Ord operator + (std::complex<double> d);
    Ord operator-(const Ord &o);
    Ord operator-(double d);
    Ord operator-(std::complex<double> d);
    Ord operator*(const Ord &o);
    Ord operator*(double d);
    Ord operator*(std::complex<double> d);
    Ord operator/(const Ord &o);
    Ord operator/(double d);
    Ord operator/(std::complex<double> d);

    Ord operator+=(const Ord &o);
    Ord operator-=(const Ord &o);

    Ord operator+=(const double &d);
    Ord operator+=(const std::complex<double> &d);
    Ord operator-=(const double &d);
    Ord operator-=(const std::complex<double> &d);
    Ord operator*=(const double &d);
    Ord operator*=(const std::complex<double> &d);
    Ord operator/=(const double &d);
    Ord operator/=(const std::complex<double> &d);

    bool operator<(double d);
    bool operator<(std::complex<double> d);
    bool operator>(double d);
    bool operator>(std::complex<double> d);
    bool operator<(const Ord &o);
    bool operator>(const Ord &o);
    
    friend std::ostream & operator<< (std::ostream& os, const Ord& ord)
    {
      os << "Integration order: " << ord.get_order() << std::endl;
      return os;
    }

  protected:
    int order;
  };

  HERMES_API Ord operator/(const double &a, const Ord &b);
  HERMES_API Ord operator*(const double &a, const Ord &b);
  HERMES_API Ord operator + (const double &a, const Ord &b);
  HERMES_API Ord operator-(const double &a, const Ord &b);
  HERMES_API Ord operator/(const std::complex<double> &a, const Ord &b);
  HERMES_API Ord operator*(const std::complex<double> &a, const Ord &b);
  HERMES_API Ord operator + (const std::complex<double> &a, const Ord &b);
  HERMES_API Ord operator-(const std::complex<double> &a, const Ord &b);
  HERMES_API Ord operator-(const Ord &a);

  HERMES_API Ord pow(const Ord &a, const double &b);
  HERMES_API Ord sqrt(const Ord &a);
  HERMES_API Ord sqr(const Ord &a);
  HERMES_API Ord conj(const Ord &a);
  HERMES_API Ord abs(const Ord &a);
  HERMES_API Ord magn(const Ord &a);

  HERMES_API Ord atan2(const Ord &a, const Ord &b);
  HERMES_API Ord atan(const Ord &a);
  HERMES_API Ord sin(const Ord &a);
  HERMES_API Ord cos(const Ord &a);
  HERMES_API Ord log(const Ord &a);
  HERMES_API Ord exp(const Ord &a);
}

#endif