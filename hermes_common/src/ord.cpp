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
/*! \file ord.cpp
\brief Contains class Ord for calculation of integration order.
*/
#include "ord.h"
#include <algorithm>

namespace Hermes
{
  Ord::Ord(): order(0) {}
  Ord::Ord(int o): order(o) {}
  Ord::Ord(double o): order(0) {}

  int Ord::get_order() const { return order; }

  Ord Ord::get_max_order() {return Ord(30);}

  Ord Ord::operator+ (const Ord &o) { return Ord(std::max(this->order, o.order)); }
  Ord Ord::operator+ (double d) { return *this; }
  Ord Ord::operator+ (std::complex<double> d) { return *this; }
  Ord Ord::operator-(const Ord &o) { return Ord(std::max(this->order, o.order)); }
  Ord Ord::operator-(double d) { return *this; }
  Ord Ord::operator-(std::complex<double> d) { return *this; }
  Ord Ord::operator*(const Ord &o) { return Ord(this->order + o.order); }
  Ord Ord::operator*(double d) { return *this; }
  Ord Ord::operator*(std::complex<double> d) { return *this; }
  Ord Ord::operator/(const Ord &o) { return Ord::get_max_order(); }
  Ord Ord::operator/(double d) { return *this; }
  Ord Ord::operator/(std::complex<double> d) { return *this; }

  Ord Ord::operator+=(const Ord &o) { this->order = std::max(this->order, o.order); return *this; }
  Ord Ord::operator-=(const Ord &o) { this->order = std::max(this->order, o.order); return *this; }

  Ord Ord::operator+=(const double &d) { return *this; }
  Ord Ord::operator+=(const std::complex<double> &d) { return *this; }
  Ord Ord::operator-=(const double &d) { return *this; }
  Ord Ord::operator-=(const std::complex<double> &d) { return *this; }
  Ord Ord::operator*=(const double &d) { return *this; }
  Ord Ord::operator*=(const std::complex<double> &d) { return *this; }
  Ord Ord::operator/=(const double &d) { return *this; }
  Ord Ord::operator/=(const std::complex<double> &d) { return *this; }

  bool Ord::operator<(double d) { return true; }
  bool Ord::operator<(std::complex<double> d) { return true; }
  bool Ord::operator>(double d) { return false; }
  bool Ord::operator>(std::complex<double> d) { return false; }
  bool Ord::operator<(const Ord &o) { return this->order < o.order; }
  bool Ord::operator>(const Ord &o) { return this->order > o.order; }

  Ord operator/(const double &a, const Ord &b) { return Ord::get_max_order(); }
  Ord operator*(const double &a, const Ord &b) { return b; }
  Ord operator + (const double &a, const Ord &b) { return b; }
  Ord operator-(const double &a, const Ord &b) { return b; }
  Ord operator/(const std::complex<double> &a, const Ord &b) { return Ord::get_max_order(); }
  Ord operator*(const std::complex<double> &a, const Ord &b) { return b; }
  Ord operator + (const std::complex<double> &a, const Ord &b) { return b; }
  Ord operator-(const std::complex<double> &a, const Ord &b) { return b; }
  Ord operator-(const Ord &a) { return a; }

  Ord pow(const Ord &a, const double &b) { return Ord((int) ceil(fabs(b)) * a.get_order()); }
  Ord sqrt(const Ord &a) { return a; }
  Ord sqr(const Ord &a) { return Ord(2 * a.get_order()); }
  Ord conj(const Ord &a) { return a; }
  Ord abs(const Ord &a) { return a; }
  Ord magn(const Ord &a) { return a; }

  Ord atan2(const Ord &a, const Ord &b) { return Ord::get_max_order(); }
  Ord atan(const Ord &a) { return Ord::get_max_order(); }
  Ord sin(const Ord &a) { return Ord::get_max_order(); }
  Ord cos(const Ord &a) { return Ord::get_max_order(); }
  Ord log(const Ord &a) { return Ord::get_max_order(); }
  Ord exp(const Ord &a) { return Ord(3 * a.get_order()); }
}