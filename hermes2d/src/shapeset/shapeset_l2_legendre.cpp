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

#include "global.h"
#include "shapeset.h"
#include "shapeset_common.h"
#include "shapeset_l2_all.h"

namespace Hermes
{
  namespace Hermes2D
  {
    static double leg_quad_l0_l0(double x, double y)
    {
      return Legendre0(x) * Legendre0(y);
    }

    static double leg_quad_l0_l0x(double x, double y)
    {
      return Legendre0x(x) * Legendre0(y);
    }

    static double leg_quad_l0_l0y(double x, double y)
    {
      return Legendre0(x) * Legendre0x(y);
    }

    static double leg_quad_l0_l0xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre0(y);
    }

    static double leg_quad_l0_l0xy(double x, double y)
    {
      return Legendre0x(x) * Legendre0x(y);
    }

    static double leg_quad_l0_l0yy(double x, double y)
    {
      return  Legendre0(x) * Legendre0xx(y);
    }

    static double leg_quad_l0_l1(double x, double y)
    {
      return Legendre0(x) * Legendre1(y);
    }

    static double leg_quad_l0_l1x(double x, double y)
    {
      return Legendre0x(x) * Legendre1(y);
    }

    static double leg_quad_l0_l1y(double x, double y)
    {
      return Legendre0(x) * Legendre1x(y);
    }

    static double leg_quad_l0_l1xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre1(y);
    }

    static double leg_quad_l0_l1xy(double x, double y)
    {
      return Legendre0x(x) * Legendre1x(y);
    }

    static double leg_quad_l0_l1yy(double x, double y)
    {
      return  Legendre0(x) * Legendre1xx(y);
    }

    static double leg_quad_l0_l2(double x, double y)
    {
      return Legendre0(x) * Legendre2(y);
    }

    static double leg_quad_l0_l2x(double x, double y)
    {
      return Legendre0x(x) * Legendre2(y);
    }

    static double leg_quad_l0_l2y(double x, double y)
    {
      return Legendre0(x) * Legendre2x(y);
    }

    static double leg_quad_l0_l2xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre2(y);
    }

    static double leg_quad_l0_l2xy(double x, double y)
    {
      return Legendre0x(x) * Legendre2x(y);
    }

    static double leg_quad_l0_l2yy(double x, double y)
    {
      return  Legendre0(x) * Legendre2xx(y);
    }

    static double leg_quad_l0_l3(double x, double y)
    {
      return Legendre0(x) * Legendre3(y);
    }

    static double leg_quad_l0_l3x(double x, double y)
    {
      return Legendre0x(x) * Legendre3(y);
    }

    static double leg_quad_l0_l3y(double x, double y)
    {
      return Legendre0(x) * Legendre3x(y);
    }

    static double leg_quad_l0_l3xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre3(y);
    }

    static double leg_quad_l0_l3xy(double x, double y)
    {
      return Legendre0x(x) * Legendre3x(y);
    }

    static double leg_quad_l0_l3yy(double x, double y)
    {
      return  Legendre0(x) * Legendre3xx(y);
    }

    static double leg_quad_l0_l4(double x, double y)
    {
      return Legendre0(x) * Legendre4(y);
    }

    static double leg_quad_l0_l4x(double x, double y)
    {
      return Legendre0x(x) * Legendre4(y);
    }

    static double leg_quad_l0_l4y(double x, double y)
    {
      return Legendre0(x) * Legendre4x(y);
    }

    static double leg_quad_l0_l4xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre4(y);
    }

    static double leg_quad_l0_l4xy(double x, double y)
    {
      return Legendre0x(x) * Legendre4x(y);
    }

    static double leg_quad_l0_l4yy(double x, double y)
    {
      return  Legendre0(x) * Legendre4xx(y);
    }

    static double leg_quad_l0_l5(double x, double y)
    {
      return Legendre0(x) * Legendre5(y);
    }

    static double leg_quad_l0_l5x(double x, double y)
    {
      return Legendre0x(x) * Legendre5(y);
    }

    static double leg_quad_l0_l5y(double x, double y)
    {
      return Legendre0(x) * Legendre5x(y);
    }

    static double leg_quad_l0_l5xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre5(y);
    }

    static double leg_quad_l0_l5xy(double x, double y)
    {
      return Legendre0x(x) * Legendre5x(y);
    }

    static double leg_quad_l0_l5yy(double x, double y)
    {
      return  Legendre0(x) * Legendre5xx(y);
    }

    static double leg_quad_l0_l6(double x, double y)
    {
      return Legendre0(x) * Legendre6(y);
    }

    static double leg_quad_l0_l6x(double x, double y)
    {
      return Legendre0x(x) * Legendre6(y);
    }

    static double leg_quad_l0_l6y(double x, double y)
    {
      return Legendre0(x) * Legendre6x(y);
    }

    static double leg_quad_l0_l6xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre6(y);
    }

    static double leg_quad_l0_l6xy(double x, double y)
    {
      return Legendre0x(x) * Legendre6x(y);
    }

    static double leg_quad_l0_l6yy(double x, double y)
    {
      return  Legendre0(x) * Legendre6xx(y);
    }

    static double leg_quad_l0_l7(double x, double y)
    {
      return Legendre0(x) * Legendre7(y);
    }

    static double leg_quad_l0_l7x(double x, double y)
    {
      return Legendre0x(x) * Legendre7(y);
    }

    static double leg_quad_l0_l7y(double x, double y)
    {
      return Legendre0(x) * Legendre7x(y);
    }

    static double leg_quad_l0_l7xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre7(y);
    }

    static double leg_quad_l0_l7xy(double x, double y)
    {
      return Legendre0x(x) * Legendre7x(y);
    }

    static double leg_quad_l0_l7yy(double x, double y)
    {
      return  Legendre0(x) * Legendre7xx(y);
    }

    static double leg_quad_l0_l8(double x, double y)
    {
      return Legendre0(x) * Legendre8(y);
    }

    static double leg_quad_l0_l8x(double x, double y)
    {
      return Legendre0x(x) * Legendre8(y);
    }

    static double leg_quad_l0_l8y(double x, double y)
    {
      return Legendre0(x) * Legendre8x(y);
    }

    static double leg_quad_l0_l8xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre8(y);
    }

    static double leg_quad_l0_l8xy(double x, double y)
    {
      return Legendre0x(x) * Legendre8x(y);
    }

    static double leg_quad_l0_l8yy(double x, double y)
    {
      return  Legendre0(x) * Legendre8xx(y);
    }

    static double leg_quad_l0_l9(double x, double y)
    {
      return Legendre0(x) * Legendre9(y);
    }

    static double leg_quad_l0_l9x(double x, double y)
    {
      return Legendre0x(x) * Legendre9(y);
    }

    static double leg_quad_l0_l9y(double x, double y)
    {
      return Legendre0(x) * Legendre9x(y);
    }

    static double leg_quad_l0_l9xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre9(y);
    }

    static double leg_quad_l0_l9xy(double x, double y)
    {
      return Legendre0x(x) * Legendre9x(y);
    }

    static double leg_quad_l0_l9yy(double x, double y)
    {
      return  Legendre0(x) * Legendre9xx(y);
    }

    static double leg_quad_l0_l10(double x, double y)
    {
      return Legendre0(x) * Legendre10(y);
    }

    static double leg_quad_l0_l10x(double x, double y)
    {
      return Legendre0x(x) * Legendre10(y);
    }

    static double leg_quad_l0_l10y(double x, double y)
    {
      return Legendre0(x) * Legendre10x(y);
    }

    static double leg_quad_l0_l10xx(double x, double y)
    {
      return Legendre0xx(x) * Legendre10(y);
    }

    static double leg_quad_l0_l10xy(double x, double y)
    {
      return Legendre0x(x) * Legendre10x(y);
    }

    static double leg_quad_l0_l10yy(double x, double y)
    {
      return  Legendre0(x) * Legendre10xx(y);
    }

    static double leg_quad_l1_l0(double x, double y)
    {
      return Legendre1(x) * Legendre0(y);
    }

    static double leg_quad_l1_l0x(double x, double y)
    {
      return Legendre1x(x) * Legendre0(y);
    }

    static double leg_quad_l1_l0y(double x, double y)
    {
      return Legendre1(x) * Legendre0x(y);
    }

    static double leg_quad_l1_l0xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre0(y);
    }

    static double leg_quad_l1_l0xy(double x, double y)
    {
      return Legendre1x(x) * Legendre0x(y);
    }

    static double leg_quad_l1_l0yy(double x, double y)
    {
      return  Legendre1(x) * Legendre0xx(y);
    }

    static double leg_quad_l1_l1(double x, double y)
    {
      return Legendre1(x) * Legendre1(y);
    }

    static double leg_quad_l1_l1x(double x, double y)
    {
      return Legendre1x(x) * Legendre1(y);
    }

    static double leg_quad_l1_l1y(double x, double y)
    {
      return Legendre1(x) * Legendre1x(y);
    }

    static double leg_quad_l1_l1xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre1(y);
    }

    static double leg_quad_l1_l1xy(double x, double y)
    {
      return Legendre1x(x) * Legendre1x(y);
    }

    static double leg_quad_l1_l1yy(double x, double y)
    {
      return  Legendre1(x) * Legendre1xx(y);
    }

    static double leg_quad_l1_l2(double x, double y)
    {
      return Legendre1(x) * Legendre2(y);
    }

    static double leg_quad_l1_l2x(double x, double y)
    {
      return Legendre1x(x) * Legendre2(y);
    }

    static double leg_quad_l1_l2y(double x, double y)
    {
      return Legendre1(x) * Legendre2x(y);
    }

    static double leg_quad_l1_l2xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre2(y);
    }

    static double leg_quad_l1_l2xy(double x, double y)
    {
      return Legendre1x(x) * Legendre2x(y);
    }

    static double leg_quad_l1_l2yy(double x, double y)
    {
      return  Legendre1(x) * Legendre2xx(y);
    }

    static double leg_quad_l1_l3(double x, double y)
    {
      return Legendre1(x) * Legendre3(y);
    }

    static double leg_quad_l1_l3x(double x, double y)
    {
      return Legendre1x(x) * Legendre3(y);
    }

    static double leg_quad_l1_l3y(double x, double y)
    {
      return Legendre1(x) * Legendre3x(y);
    }

    static double leg_quad_l1_l3xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre3(y);
    }

    static double leg_quad_l1_l3xy(double x, double y)
    {
      return Legendre1x(x) * Legendre3x(y);
    }

    static double leg_quad_l1_l3yy(double x, double y)
    {
      return  Legendre1(x) * Legendre3xx(y);
    }

    static double leg_quad_l1_l4(double x, double y)
    {
      return Legendre1(x) * Legendre4(y);
    }

    static double leg_quad_l1_l4x(double x, double y)
    {
      return Legendre1x(x) * Legendre4(y);
    }

    static double leg_quad_l1_l4y(double x, double y)
    {
      return Legendre1(x) * Legendre4x(y);
    }

    static double leg_quad_l1_l4xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre4(y);
    }

    static double leg_quad_l1_l4xy(double x, double y)
    {
      return Legendre1x(x) * Legendre4x(y);
    }

    static double leg_quad_l1_l4yy(double x, double y)
    {
      return  Legendre1(x) * Legendre4xx(y);
    }

    static double leg_quad_l1_l5(double x, double y)
    {
      return Legendre1(x) * Legendre5(y);
    }

    static double leg_quad_l1_l5x(double x, double y)
    {
      return Legendre1x(x) * Legendre5(y);
    }

    static double leg_quad_l1_l5y(double x, double y)
    {
      return Legendre1(x) * Legendre5x(y);
    }

    static double leg_quad_l1_l5xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre5(y);
    }

    static double leg_quad_l1_l5xy(double x, double y)
    {
      return Legendre1x(x) * Legendre5x(y);
    }

    static double leg_quad_l1_l5yy(double x, double y)
    {
      return  Legendre1(x) * Legendre5xx(y);
    }

    static double leg_quad_l1_l6(double x, double y)
    {
      return Legendre1(x) * Legendre6(y);
    }

    static double leg_quad_l1_l6x(double x, double y)
    {
      return Legendre1x(x) * Legendre6(y);
    }

    static double leg_quad_l1_l6y(double x, double y)
    {
      return Legendre1(x) * Legendre6x(y);
    }

    static double leg_quad_l1_l6xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre6(y);
    }

    static double leg_quad_l1_l6xy(double x, double y)
    {
      return Legendre1x(x) * Legendre6x(y);
    }

    static double leg_quad_l1_l6yy(double x, double y)
    {
      return  Legendre1(x) * Legendre6xx(y);
    }

    static double leg_quad_l1_l7(double x, double y)
    {
      return Legendre1(x) * Legendre7(y);
    }

    static double leg_quad_l1_l7x(double x, double y)
    {
      return Legendre1x(x) * Legendre7(y);
    }

    static double leg_quad_l1_l7y(double x, double y)
    {
      return Legendre1(x) * Legendre7x(y);
    }

    static double leg_quad_l1_l7xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre7(y);
    }

    static double leg_quad_l1_l7xy(double x, double y)
    {
      return Legendre1x(x) * Legendre7x(y);
    }

    static double leg_quad_l1_l7yy(double x, double y)
    {
      return  Legendre1(x) * Legendre7xx(y);
    }

    static double leg_quad_l1_l8(double x, double y)
    {
      return Legendre1(x) * Legendre8(y);
    }

    static double leg_quad_l1_l8x(double x, double y)
    {
      return Legendre1x(x) * Legendre8(y);
    }

    static double leg_quad_l1_l8y(double x, double y)
    {
      return Legendre1(x) * Legendre8x(y);
    }

    static double leg_quad_l1_l8xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre8(y);
    }

    static double leg_quad_l1_l8xy(double x, double y)
    {
      return Legendre1x(x) * Legendre8x(y);
    }

    static double leg_quad_l1_l8yy(double x, double y)
    {
      return  Legendre1(x) * Legendre8xx(y);
    }

    static double leg_quad_l1_l9(double x, double y)
    {
      return Legendre1(x) * Legendre9(y);
    }

    static double leg_quad_l1_l9x(double x, double y)
    {
      return Legendre1x(x) * Legendre9(y);
    }

    static double leg_quad_l1_l9y(double x, double y)
    {
      return Legendre1(x) * Legendre9x(y);
    }

    static double leg_quad_l1_l9xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre9(y);
    }

    static double leg_quad_l1_l9xy(double x, double y)
    {
      return Legendre1x(x) * Legendre9x(y);
    }

    static double leg_quad_l1_l9yy(double x, double y)
    {
      return  Legendre1(x) * Legendre9xx(y);
    }

    static double leg_quad_l1_l10(double x, double y)
    {
      return Legendre1(x) * Legendre10(y);
    }

    static double leg_quad_l1_l10x(double x, double y)
    {
      return Legendre1x(x) * Legendre10(y);
    }

    static double leg_quad_l1_l10y(double x, double y)
    {
      return Legendre1(x) * Legendre10x(y);
    }

    static double leg_quad_l1_l10xx(double x, double y)
    {
      return Legendre1xx(x) * Legendre10(y);
    }

    static double leg_quad_l1_l10xy(double x, double y)
    {
      return Legendre1x(x) * Legendre10x(y);
    }

    static double leg_quad_l1_l10yy(double x, double y)
    {
      return  Legendre1(x) * Legendre10xx(y);
    }

    static double leg_quad_l2_l0(double x, double y)
    {
      return Legendre2(x) * Legendre0(y);
    }

    static double leg_quad_l2_l0x(double x, double y)
    {
      return Legendre2x(x) * Legendre0(y);
    }

    static double leg_quad_l2_l0y(double x, double y)
    {
      return Legendre2(x) * Legendre0x(y);
    }

    static double leg_quad_l2_l0xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre0(y);
    }

    static double leg_quad_l2_l0xy(double x, double y)
    {
      return Legendre2x(x) * Legendre0x(y);
    }

    static double leg_quad_l2_l0yy(double x, double y)
    {
      return  Legendre2(x) * Legendre0xx(y);
    }

    static double leg_quad_l2_l1(double x, double y)
    {
      return Legendre2(x) * Legendre1(y);
    }

    static double leg_quad_l2_l1x(double x, double y)
    {
      return Legendre2x(x) * Legendre1(y);
    }

    static double leg_quad_l2_l1y(double x, double y)
    {
      return Legendre2(x) * Legendre1x(y);
    }

    static double leg_quad_l2_l1xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre1(y);
    }

    static double leg_quad_l2_l1xy(double x, double y)
    {
      return Legendre2x(x) * Legendre1x(y);
    }

    static double leg_quad_l2_l1yy(double x, double y)
    {
      return  Legendre2(x) * Legendre1xx(y);
    }

    static double leg_quad_l2_l2(double x, double y)
    {
      return Legendre2(x) * Legendre2(y);
    }

    static double leg_quad_l2_l2x(double x, double y)
    {
      return Legendre2x(x) * Legendre2(y);
    }

    static double leg_quad_l2_l2y(double x, double y)
    {
      return Legendre2(x) * Legendre2x(y);
    }

    static double leg_quad_l2_l2xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre2(y);
    }

    static double leg_quad_l2_l2xy(double x, double y)
    {
      return Legendre2x(x) * Legendre2x(y);
    }

    static double leg_quad_l2_l2yy(double x, double y)
    {
      return  Legendre2(x) * Legendre2xx(y);
    }

    static double leg_quad_l2_l3(double x, double y)
    {
      return Legendre2(x) * Legendre3(y);
    }

    static double leg_quad_l2_l3x(double x, double y)
    {
      return Legendre2x(x) * Legendre3(y);
    }

    static double leg_quad_l2_l3y(double x, double y)
    {
      return Legendre2(x) * Legendre3x(y);
    }

    static double leg_quad_l2_l3xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre3(y);
    }

    static double leg_quad_l2_l3xy(double x, double y)
    {
      return Legendre2x(x) * Legendre3x(y);
    }

    static double leg_quad_l2_l3yy(double x, double y)
    {
      return  Legendre2(x) * Legendre3xx(y);
    }

    static double leg_quad_l2_l4(double x, double y)
    {
      return Legendre2(x) * Legendre4(y);
    }

    static double leg_quad_l2_l4x(double x, double y)
    {
      return Legendre2x(x) * Legendre4(y);
    }

    static double leg_quad_l2_l4y(double x, double y)
    {
      return Legendre2(x) * Legendre4x(y);
    }

    static double leg_quad_l2_l4xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre4(y);
    }

    static double leg_quad_l2_l4xy(double x, double y)
    {
      return Legendre2x(x) * Legendre4x(y);
    }

    static double leg_quad_l2_l4yy(double x, double y)
    {
      return  Legendre2(x) * Legendre4xx(y);
    }

    static double leg_quad_l2_l5(double x, double y)
    {
      return Legendre2(x) * Legendre5(y);
    }

    static double leg_quad_l2_l5x(double x, double y)
    {
      return Legendre2x(x) * Legendre5(y);
    }

    static double leg_quad_l2_l5y(double x, double y)
    {
      return Legendre2(x) * Legendre5x(y);
    }

    static double leg_quad_l2_l5xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre5(y);
    }

    static double leg_quad_l2_l5xy(double x, double y)
    {
      return Legendre2x(x) * Legendre5x(y);
    }

    static double leg_quad_l2_l5yy(double x, double y)
    {
      return  Legendre2(x) * Legendre5xx(y);
    }

    static double leg_quad_l2_l6(double x, double y)
    {
      return Legendre2(x) * Legendre6(y);
    }

    static double leg_quad_l2_l6x(double x, double y)
    {
      return Legendre2x(x) * Legendre6(y);
    }

    static double leg_quad_l2_l6y(double x, double y)
    {
      return Legendre2(x) * Legendre6x(y);
    }

    static double leg_quad_l2_l6xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre6(y);
    }

    static double leg_quad_l2_l6xy(double x, double y)
    {
      return Legendre2x(x) * Legendre6x(y);
    }

    static double leg_quad_l2_l6yy(double x, double y)
    {
      return  Legendre2(x) * Legendre6xx(y);
    }

    static double leg_quad_l2_l7(double x, double y)
    {
      return Legendre2(x) * Legendre7(y);
    }

    static double leg_quad_l2_l7x(double x, double y)
    {
      return Legendre2x(x) * Legendre7(y);
    }

    static double leg_quad_l2_l7y(double x, double y)
    {
      return Legendre2(x) * Legendre7x(y);
    }

    static double leg_quad_l2_l7xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre7(y);
    }

    static double leg_quad_l2_l7xy(double x, double y)
    {
      return Legendre2x(x) * Legendre7x(y);
    }

    static double leg_quad_l2_l7yy(double x, double y)
    {
      return  Legendre2(x) * Legendre7xx(y);
    }

    static double leg_quad_l2_l8(double x, double y)
    {
      return Legendre2(x) * Legendre8(y);
    }

    static double leg_quad_l2_l8x(double x, double y)
    {
      return Legendre2x(x) * Legendre8(y);
    }

    static double leg_quad_l2_l8y(double x, double y)
    {
      return Legendre2(x) * Legendre8x(y);
    }

    static double leg_quad_l2_l8xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre8(y);
    }

    static double leg_quad_l2_l8xy(double x, double y)
    {
      return Legendre2x(x) * Legendre8x(y);
    }

    static double leg_quad_l2_l8yy(double x, double y)
    {
      return  Legendre2(x) * Legendre8xx(y);
    }

    static double leg_quad_l2_l9(double x, double y)
    {
      return Legendre2(x) * Legendre9(y);
    }

    static double leg_quad_l2_l9x(double x, double y)
    {
      return Legendre2x(x) * Legendre9(y);
    }

    static double leg_quad_l2_l9y(double x, double y)
    {
      return Legendre2(x) * Legendre9x(y);
    }

    static double leg_quad_l2_l9xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre9(y);
    }

    static double leg_quad_l2_l9xy(double x, double y)
    {
      return Legendre2x(x) * Legendre9x(y);
    }

    static double leg_quad_l2_l9yy(double x, double y)
    {
      return  Legendre2(x) * Legendre9xx(y);
    }

    static double leg_quad_l2_l10(double x, double y)
    {
      return Legendre2(x) * Legendre10(y);
    }

    static double leg_quad_l2_l10x(double x, double y)
    {
      return Legendre2x(x) * Legendre10(y);
    }

    static double leg_quad_l2_l10y(double x, double y)
    {
      return Legendre2(x) * Legendre10x(y);
    }

    static double leg_quad_l2_l10xx(double x, double y)
    {
      return Legendre2xx(x) * Legendre10(y);
    }

    static double leg_quad_l2_l10xy(double x, double y)
    {
      return Legendre2x(x) * Legendre10x(y);
    }

    static double leg_quad_l2_l10yy(double x, double y)
    {
      return  Legendre2(x) * Legendre10xx(y);
    }

    static double leg_quad_l3_l0(double x, double y)
    {
      return Legendre3(x) * Legendre0(y);
    }

    static double leg_quad_l3_l0x(double x, double y)
    {
      return Legendre3x(x) * Legendre0(y);
    }

    static double leg_quad_l3_l0y(double x, double y)
    {
      return Legendre3(x) * Legendre0x(y);
    }

    static double leg_quad_l3_l0xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre0(y);
    }

    static double leg_quad_l3_l0xy(double x, double y)
    {
      return Legendre3x(x) * Legendre0x(y);
    }

    static double leg_quad_l3_l0yy(double x, double y)
    {
      return  Legendre3(x) * Legendre0xx(y);
    }

    static double leg_quad_l3_l1(double x, double y)
    {
      return Legendre3(x) * Legendre1(y);
    }

    static double leg_quad_l3_l1x(double x, double y)
    {
      return Legendre3x(x) * Legendre1(y);
    }

    static double leg_quad_l3_l1y(double x, double y)
    {
      return Legendre3(x) * Legendre1x(y);
    }

    static double leg_quad_l3_l1xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre1(y);
    }

    static double leg_quad_l3_l1xy(double x, double y)
    {
      return Legendre3x(x) * Legendre1x(y);
    }

    static double leg_quad_l3_l1yy(double x, double y)
    {
      return  Legendre3(x) * Legendre1xx(y);
    }

    static double leg_quad_l3_l2(double x, double y)
    {
      return Legendre3(x) * Legendre2(y);
    }

    static double leg_quad_l3_l2x(double x, double y)
    {
      return Legendre3x(x) * Legendre2(y);
    }

    static double leg_quad_l3_l2y(double x, double y)
    {
      return Legendre3(x) * Legendre2x(y);
    }

    static double leg_quad_l3_l2xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre2(y);
    }

    static double leg_quad_l3_l2xy(double x, double y)
    {
      return Legendre3x(x) * Legendre2x(y);
    }

    static double leg_quad_l3_l2yy(double x, double y)
    {
      return  Legendre3(x) * Legendre2xx(y);
    }

    static double leg_quad_l3_l3(double x, double y)
    {
      return Legendre3(x) * Legendre3(y);
    }

    static double leg_quad_l3_l3x(double x, double y)
    {
      return Legendre3x(x) * Legendre3(y);
    }

    static double leg_quad_l3_l3y(double x, double y)
    {
      return Legendre3(x) * Legendre3x(y);
    }

    static double leg_quad_l3_l3xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre3(y);
    }

    static double leg_quad_l3_l3xy(double x, double y)
    {
      return Legendre3x(x) * Legendre3x(y);
    }

    static double leg_quad_l3_l3yy(double x, double y)
    {
      return  Legendre3(x) * Legendre3xx(y);
    }

    static double leg_quad_l3_l4(double x, double y)
    {
      return Legendre3(x) * Legendre4(y);
    }

    static double leg_quad_l3_l4x(double x, double y)
    {
      return Legendre3x(x) * Legendre4(y);
    }

    static double leg_quad_l3_l4y(double x, double y)
    {
      return Legendre3(x) * Legendre4x(y);
    }

    static double leg_quad_l3_l4xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre4(y);
    }

    static double leg_quad_l3_l4xy(double x, double y)
    {
      return Legendre3x(x) * Legendre4x(y);
    }

    static double leg_quad_l3_l4yy(double x, double y)
    {
      return  Legendre3(x) * Legendre4xx(y);
    }

    static double leg_quad_l3_l5(double x, double y)
    {
      return Legendre3(x) * Legendre5(y);
    }

    static double leg_quad_l3_l5x(double x, double y)
    {
      return Legendre3x(x) * Legendre5(y);
    }

    static double leg_quad_l3_l5y(double x, double y)
    {
      return Legendre3(x) * Legendre5x(y);
    }

    static double leg_quad_l3_l5xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre5(y);
    }

    static double leg_quad_l3_l5xy(double x, double y)
    {
      return Legendre3x(x) * Legendre5x(y);
    }

    static double leg_quad_l3_l5yy(double x, double y)
    {
      return  Legendre3(x) * Legendre5xx(y);
    }

    static double leg_quad_l3_l6(double x, double y)
    {
      return Legendre3(x) * Legendre6(y);
    }

    static double leg_quad_l3_l6x(double x, double y)
    {
      return Legendre3x(x) * Legendre6(y);
    }

    static double leg_quad_l3_l6y(double x, double y)
    {
      return Legendre3(x) * Legendre6x(y);
    }

    static double leg_quad_l3_l6xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre6(y);
    }

    static double leg_quad_l3_l6xy(double x, double y)
    {
      return Legendre3x(x) * Legendre6x(y);
    }

    static double leg_quad_l3_l6yy(double x, double y)
    {
      return  Legendre3(x) * Legendre6xx(y);
    }

    static double leg_quad_l3_l7(double x, double y)
    {
      return Legendre3(x) * Legendre7(y);
    }

    static double leg_quad_l3_l7x(double x, double y)
    {
      return Legendre3x(x) * Legendre7(y);
    }

    static double leg_quad_l3_l7y(double x, double y)
    {
      return Legendre3(x) * Legendre7x(y);
    }

    static double leg_quad_l3_l7xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre7(y);
    }

    static double leg_quad_l3_l7xy(double x, double y)
    {
      return Legendre3x(x) * Legendre7x(y);
    }

    static double leg_quad_l3_l7yy(double x, double y)
    {
      return  Legendre3(x) * Legendre7xx(y);
    }

    static double leg_quad_l3_l8(double x, double y)
    {
      return Legendre3(x) * Legendre8(y);
    }

    static double leg_quad_l3_l8x(double x, double y)
    {
      return Legendre3x(x) * Legendre8(y);
    }

    static double leg_quad_l3_l8y(double x, double y)
    {
      return Legendre3(x) * Legendre8x(y);
    }

    static double leg_quad_l3_l8xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre8(y);
    }

    static double leg_quad_l3_l8xy(double x, double y)
    {
      return Legendre3x(x) * Legendre8x(y);
    }

    static double leg_quad_l3_l8yy(double x, double y)
    {
      return  Legendre3(x) * Legendre8xx(y);
    }

    static double leg_quad_l3_l9(double x, double y)
    {
      return Legendre3(x) * Legendre9(y);
    }

    static double leg_quad_l3_l9x(double x, double y)
    {
      return Legendre3x(x) * Legendre9(y);
    }

    static double leg_quad_l3_l9y(double x, double y)
    {
      return Legendre3(x) * Legendre9x(y);
    }

    static double leg_quad_l3_l9xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre9(y);
    }

    static double leg_quad_l3_l9xy(double x, double y)
    {
      return Legendre3x(x) * Legendre9x(y);
    }

    static double leg_quad_l3_l9yy(double x, double y)
    {
      return  Legendre3(x) * Legendre9xx(y);
    }

    static double leg_quad_l3_l10(double x, double y)
    {
      return Legendre3(x) * Legendre10(y);
    }

    static double leg_quad_l3_l10x(double x, double y)
    {
      return Legendre3x(x) * Legendre10(y);
    }

    static double leg_quad_l3_l10y(double x, double y)
    {
      return Legendre3(x) * Legendre10x(y);
    }

    static double leg_quad_l3_l10xx(double x, double y)
    {
      return Legendre3xx(x) * Legendre10(y);
    }

    static double leg_quad_l3_l10xy(double x, double y)
    {
      return Legendre3x(x) * Legendre10x(y);
    }

    static double leg_quad_l3_l10yy(double x, double y)
    {
      return  Legendre3(x) * Legendre10xx(y);
    }

    static double leg_quad_l4_l0(double x, double y)
    {
      return Legendre4(x) * Legendre0(y);
    }

    static double leg_quad_l4_l0x(double x, double y)
    {
      return Legendre4x(x) * Legendre0(y);
    }

    static double leg_quad_l4_l0y(double x, double y)
    {
      return Legendre4(x) * Legendre0x(y);
    }

    static double leg_quad_l4_l0xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre0(y);
    }

    static double leg_quad_l4_l0xy(double x, double y)
    {
      return Legendre4x(x) * Legendre0x(y);
    }

    static double leg_quad_l4_l0yy(double x, double y)
    {
      return  Legendre4(x) * Legendre0xx(y);
    }

    static double leg_quad_l4_l1(double x, double y)
    {
      return Legendre4(x) * Legendre1(y);
    }

    static double leg_quad_l4_l1x(double x, double y)
    {
      return Legendre4x(x) * Legendre1(y);
    }

    static double leg_quad_l4_l1y(double x, double y)
    {
      return Legendre4(x) * Legendre1x(y);
    }

    static double leg_quad_l4_l1xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre1(y);
    }

    static double leg_quad_l4_l1xy(double x, double y)
    {
      return Legendre4x(x) * Legendre1x(y);
    }

    static double leg_quad_l4_l1yy(double x, double y)
    {
      return  Legendre4(x) * Legendre1xx(y);
    }

    static double leg_quad_l4_l2(double x, double y)
    {
      return Legendre4(x) * Legendre2(y);
    }

    static double leg_quad_l4_l2x(double x, double y)
    {
      return Legendre4x(x) * Legendre2(y);
    }

    static double leg_quad_l4_l2y(double x, double y)
    {
      return Legendre4(x) * Legendre2x(y);
    }

    static double leg_quad_l4_l2xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre2(y);
    }

    static double leg_quad_l4_l2xy(double x, double y)
    {
      return Legendre4x(x) * Legendre2x(y);
    }

    static double leg_quad_l4_l2yy(double x, double y)
    {
      return  Legendre4(x) * Legendre2xx(y);
    }

    static double leg_quad_l4_l3(double x, double y)
    {
      return Legendre4(x) * Legendre3(y);
    }

    static double leg_quad_l4_l3x(double x, double y)
    {
      return Legendre4x(x) * Legendre3(y);
    }

    static double leg_quad_l4_l3y(double x, double y)
    {
      return Legendre4(x) * Legendre3x(y);
    }

    static double leg_quad_l4_l3xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre3(y);
    }

    static double leg_quad_l4_l3xy(double x, double y)
    {
      return Legendre4x(x) * Legendre3x(y);
    }

    static double leg_quad_l4_l3yy(double x, double y)
    {
      return  Legendre4(x) * Legendre3xx(y);
    }

    static double leg_quad_l4_l4(double x, double y)
    {
      return Legendre4(x) * Legendre4(y);
    }

    static double leg_quad_l4_l4x(double x, double y)
    {
      return Legendre4x(x) * Legendre4(y);
    }

    static double leg_quad_l4_l4y(double x, double y)
    {
      return Legendre4(x) * Legendre4x(y);
    }

    static double leg_quad_l4_l4xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre4(y);
    }

    static double leg_quad_l4_l4xy(double x, double y)
    {
      return Legendre4x(x) * Legendre4x(y);
    }

    static double leg_quad_l4_l4yy(double x, double y)
    {
      return  Legendre4(x) * Legendre4xx(y);
    }

    static double leg_quad_l4_l5(double x, double y)
    {
      return Legendre4(x) * Legendre5(y);
    }

    static double leg_quad_l4_l5x(double x, double y)
    {
      return Legendre4x(x) * Legendre5(y);
    }

    static double leg_quad_l4_l5y(double x, double y)
    {
      return Legendre4(x) * Legendre5x(y);
    }

    static double leg_quad_l4_l5xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre5(y);
    }

    static double leg_quad_l4_l5xy(double x, double y)
    {
      return Legendre4x(x) * Legendre5x(y);
    }

    static double leg_quad_l4_l5yy(double x, double y)
    {
      return  Legendre4(x) * Legendre5xx(y);
    }

    static double leg_quad_l4_l6(double x, double y)
    {
      return Legendre4(x) * Legendre6(y);
    }

    static double leg_quad_l4_l6x(double x, double y)
    {
      return Legendre4x(x) * Legendre6(y);
    }

    static double leg_quad_l4_l6y(double x, double y)
    {
      return Legendre4(x) * Legendre6x(y);
    }

    static double leg_quad_l4_l6xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre6(y);
    }

    static double leg_quad_l4_l6xy(double x, double y)
    {
      return Legendre4x(x) * Legendre6x(y);
    }

    static double leg_quad_l4_l6yy(double x, double y)
    {
      return  Legendre4(x) * Legendre6xx(y);
    }

    static double leg_quad_l4_l7(double x, double y)
    {
      return Legendre4(x) * Legendre7(y);
    }

    static double leg_quad_l4_l7x(double x, double y)
    {
      return Legendre4x(x) * Legendre7(y);
    }

    static double leg_quad_l4_l7y(double x, double y)
    {
      return Legendre4(x) * Legendre7x(y);
    }

    static double leg_quad_l4_l7xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre7(y);
    }

    static double leg_quad_l4_l7xy(double x, double y)
    {
      return Legendre4x(x) * Legendre7x(y);
    }

    static double leg_quad_l4_l7yy(double x, double y)
    {
      return  Legendre4(x) * Legendre7xx(y);
    }

    static double leg_quad_l4_l8(double x, double y)
    {
      return Legendre4(x) * Legendre8(y);
    }

    static double leg_quad_l4_l8x(double x, double y)
    {
      return Legendre4x(x) * Legendre8(y);
    }

    static double leg_quad_l4_l8y(double x, double y)
    {
      return Legendre4(x) * Legendre8x(y);
    }

    static double leg_quad_l4_l8xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre8(y);
    }

    static double leg_quad_l4_l8xy(double x, double y)
    {
      return Legendre4x(x) * Legendre8x(y);
    }

    static double leg_quad_l4_l8yy(double x, double y)
    {
      return  Legendre4(x) * Legendre8xx(y);
    }

    static double leg_quad_l4_l9(double x, double y)
    {
      return Legendre4(x) * Legendre9(y);
    }

    static double leg_quad_l4_l9x(double x, double y)
    {
      return Legendre4x(x) * Legendre9(y);
    }

    static double leg_quad_l4_l9y(double x, double y)
    {
      return Legendre4(x) * Legendre9x(y);
    }

    static double leg_quad_l4_l9xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre9(y);
    }

    static double leg_quad_l4_l9xy(double x, double y)
    {
      return Legendre4x(x) * Legendre9x(y);
    }

    static double leg_quad_l4_l9yy(double x, double y)
    {
      return  Legendre4(x) * Legendre9xx(y);
    }

    static double leg_quad_l4_l10(double x, double y)
    {
      return Legendre4(x) * Legendre10(y);
    }

    static double leg_quad_l4_l10x(double x, double y)
    {
      return Legendre4x(x) * Legendre10(y);
    }

    static double leg_quad_l4_l10y(double x, double y)
    {
      return Legendre4(x) * Legendre10x(y);
    }

    static double leg_quad_l4_l10xx(double x, double y)
    {
      return Legendre4xx(x) * Legendre10(y);
    }

    static double leg_quad_l4_l10xy(double x, double y)
    {
      return Legendre4x(x) * Legendre10x(y);
    }

    static double leg_quad_l4_l10yy(double x, double y)
    {
      return  Legendre4(x) * Legendre10xx(y);
    }

    static double leg_quad_l5_l0(double x, double y)
    {
      return Legendre5(x) * Legendre0(y);
    }

    static double leg_quad_l5_l0x(double x, double y)
    {
      return Legendre5x(x) * Legendre0(y);
    }

    static double leg_quad_l5_l0y(double x, double y)
    {
      return Legendre5(x) * Legendre0x(y);
    }

    static double leg_quad_l5_l0xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre0(y);
    }

    static double leg_quad_l5_l0xy(double x, double y)
    {
      return Legendre5x(x) * Legendre0x(y);
    }

    static double leg_quad_l5_l0yy(double x, double y)
    {
      return  Legendre5(x) * Legendre0xx(y);
    }

    static double leg_quad_l5_l1(double x, double y)
    {
      return Legendre5(x) * Legendre1(y);
    }

    static double leg_quad_l5_l1x(double x, double y)
    {
      return Legendre5x(x) * Legendre1(y);
    }

    static double leg_quad_l5_l1y(double x, double y)
    {
      return Legendre5(x) * Legendre1x(y);
    }

    static double leg_quad_l5_l1xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre1(y);
    }

    static double leg_quad_l5_l1xy(double x, double y)
    {
      return Legendre5x(x) * Legendre1x(y);
    }

    static double leg_quad_l5_l1yy(double x, double y)
    {
      return  Legendre5(x) * Legendre1xx(y);
    }

    static double leg_quad_l5_l2(double x, double y)
    {
      return Legendre5(x) * Legendre2(y);
    }

    static double leg_quad_l5_l2x(double x, double y)
    {
      return Legendre5x(x) * Legendre2(y);
    }

    static double leg_quad_l5_l2y(double x, double y)
    {
      return Legendre5(x) * Legendre2x(y);
    }

    static double leg_quad_l5_l2xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre2(y);
    }

    static double leg_quad_l5_l2xy(double x, double y)
    {
      return Legendre5x(x) * Legendre2x(y);
    }

    static double leg_quad_l5_l2yy(double x, double y)
    {
      return  Legendre5(x) * Legendre2xx(y);
    }

    static double leg_quad_l5_l3(double x, double y)
    {
      return Legendre5(x) * Legendre3(y);
    }

    static double leg_quad_l5_l3x(double x, double y)
    {
      return Legendre5x(x) * Legendre3(y);
    }

    static double leg_quad_l5_l3y(double x, double y)
    {
      return Legendre5(x) * Legendre3x(y);
    }

    static double leg_quad_l5_l3xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre3(y);
    }

    static double leg_quad_l5_l3xy(double x, double y)
    {
      return Legendre5x(x) * Legendre3x(y);
    }

    static double leg_quad_l5_l3yy(double x, double y)
    {
      return  Legendre5(x) * Legendre3xx(y);
    }

    static double leg_quad_l5_l4(double x, double y)
    {
      return Legendre5(x) * Legendre4(y);
    }

    static double leg_quad_l5_l4x(double x, double y)
    {
      return Legendre5x(x) * Legendre4(y);
    }

    static double leg_quad_l5_l4y(double x, double y)
    {
      return Legendre5(x) * Legendre4x(y);
    }

    static double leg_quad_l5_l4xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre4(y);
    }

    static double leg_quad_l5_l4xy(double x, double y)
    {
      return Legendre5x(x) * Legendre4x(y);
    }

    static double leg_quad_l5_l4yy(double x, double y)
    {
      return  Legendre5(x) * Legendre4xx(y);
    }

    static double leg_quad_l5_l5(double x, double y)
    {
      return Legendre5(x) * Legendre5(y);
    }

    static double leg_quad_l5_l5x(double x, double y)
    {
      return Legendre5x(x) * Legendre5(y);
    }

    static double leg_quad_l5_l5y(double x, double y)
    {
      return Legendre5(x) * Legendre5x(y);
    }

    static double leg_quad_l5_l5xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre5(y);
    }

    static double leg_quad_l5_l5xy(double x, double y)
    {
      return Legendre5x(x) * Legendre5x(y);
    }

    static double leg_quad_l5_l5yy(double x, double y)
    {
      return  Legendre5(x) * Legendre5xx(y);
    }

    static double leg_quad_l5_l6(double x, double y)
    {
      return Legendre5(x) * Legendre6(y);
    }

    static double leg_quad_l5_l6x(double x, double y)
    {
      return Legendre5x(x) * Legendre6(y);
    }

    static double leg_quad_l5_l6y(double x, double y)
    {
      return Legendre5(x) * Legendre6x(y);
    }

    static double leg_quad_l5_l6xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre6(y);
    }

    static double leg_quad_l5_l6xy(double x, double y)
    {
      return Legendre5x(x) * Legendre6x(y);
    }

    static double leg_quad_l5_l6yy(double x, double y)
    {
      return  Legendre5(x) * Legendre6xx(y);
    }

    static double leg_quad_l5_l7(double x, double y)
    {
      return Legendre5(x) * Legendre7(y);
    }

    static double leg_quad_l5_l7x(double x, double y)
    {
      return Legendre5x(x) * Legendre7(y);
    }

    static double leg_quad_l5_l7y(double x, double y)
    {
      return Legendre5(x) * Legendre7x(y);
    }

    static double leg_quad_l5_l7xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre7(y);
    }

    static double leg_quad_l5_l7xy(double x, double y)
    {
      return Legendre5x(x) * Legendre7x(y);
    }

    static double leg_quad_l5_l7yy(double x, double y)
    {
      return  Legendre5(x) * Legendre7xx(y);
    }

    static double leg_quad_l5_l8(double x, double y)
    {
      return Legendre5(x) * Legendre8(y);
    }

    static double leg_quad_l5_l8x(double x, double y)
    {
      return Legendre5x(x) * Legendre8(y);
    }

    static double leg_quad_l5_l8y(double x, double y)
    {
      return Legendre5(x) * Legendre8x(y);
    }

    static double leg_quad_l5_l8xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre8(y);
    }

    static double leg_quad_l5_l8xy(double x, double y)
    {
      return Legendre5x(x) * Legendre8x(y);
    }

    static double leg_quad_l5_l8yy(double x, double y)
    {
      return  Legendre5(x) * Legendre8xx(y);
    }

    static double leg_quad_l5_l9(double x, double y)
    {
      return Legendre5(x) * Legendre9(y);
    }

    static double leg_quad_l5_l9x(double x, double y)
    {
      return Legendre5x(x) * Legendre9(y);
    }

    static double leg_quad_l5_l9y(double x, double y)
    {
      return Legendre5(x) * Legendre9x(y);
    }

    static double leg_quad_l5_l9xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre9(y);
    }

    static double leg_quad_l5_l9xy(double x, double y)
    {
      return Legendre5x(x) * Legendre9x(y);
    }

    static double leg_quad_l5_l9yy(double x, double y)
    {
      return  Legendre5(x) * Legendre9xx(y);
    }

    static double leg_quad_l5_l10(double x, double y)
    {
      return Legendre5(x) * Legendre10(y);
    }

    static double leg_quad_l5_l10x(double x, double y)
    {
      return Legendre5x(x) * Legendre10(y);
    }

    static double leg_quad_l5_l10y(double x, double y)
    {
      return Legendre5(x) * Legendre10x(y);
    }

    static double leg_quad_l5_l10xx(double x, double y)
    {
      return Legendre5xx(x) * Legendre10(y);
    }

    static double leg_quad_l5_l10xy(double x, double y)
    {
      return Legendre5x(x) * Legendre10x(y);
    }

    static double leg_quad_l5_l10yy(double x, double y)
    {
      return  Legendre5(x) * Legendre10xx(y);
    }

    static double leg_quad_l6_l0(double x, double y)
    {
      return Legendre6(x) * Legendre0(y);
    }

    static double leg_quad_l6_l0x(double x, double y)
    {
      return Legendre6x(x) * Legendre0(y);
    }

    static double leg_quad_l6_l0y(double x, double y)
    {
      return Legendre6(x) * Legendre0x(y);
    }

    static double leg_quad_l6_l0xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre0(y);
    }

    static double leg_quad_l6_l0xy(double x, double y)
    {
      return Legendre6x(x) * Legendre0x(y);
    }

    static double leg_quad_l6_l0yy(double x, double y)
    {
      return  Legendre6(x) * Legendre0xx(y);
    }

    static double leg_quad_l6_l1(double x, double y)
    {
      return Legendre6(x) * Legendre1(y);
    }

    static double leg_quad_l6_l1x(double x, double y)
    {
      return Legendre6x(x) * Legendre1(y);
    }

    static double leg_quad_l6_l1y(double x, double y)
    {
      return Legendre6(x) * Legendre1x(y);
    }

    static double leg_quad_l6_l1xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre1(y);
    }

    static double leg_quad_l6_l1xy(double x, double y)
    {
      return Legendre6x(x) * Legendre1x(y);
    }

    static double leg_quad_l6_l1yy(double x, double y)
    {
      return  Legendre6(x) * Legendre1xx(y);
    }

    static double leg_quad_l6_l2(double x, double y)
    {
      return Legendre6(x) * Legendre2(y);
    }

    static double leg_quad_l6_l2x(double x, double y)
    {
      return Legendre6x(x) * Legendre2(y);
    }

    static double leg_quad_l6_l2y(double x, double y)
    {
      return Legendre6(x) * Legendre2x(y);
    }

    static double leg_quad_l6_l2xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre2(y);
    }

    static double leg_quad_l6_l2xy(double x, double y)
    {
      return Legendre6x(x) * Legendre2x(y);
    }

    static double leg_quad_l6_l2yy(double x, double y)
    {
      return  Legendre6(x) * Legendre2xx(y);
    }

    static double leg_quad_l6_l3(double x, double y)
    {
      return Legendre6(x) * Legendre3(y);
    }

    static double leg_quad_l6_l3x(double x, double y)
    {
      return Legendre6x(x) * Legendre3(y);
    }

    static double leg_quad_l6_l3y(double x, double y)
    {
      return Legendre6(x) * Legendre3x(y);
    }

    static double leg_quad_l6_l3xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre3(y);
    }

    static double leg_quad_l6_l3xy(double x, double y)
    {
      return Legendre6x(x) * Legendre3x(y);
    }

    static double leg_quad_l6_l3yy(double x, double y)
    {
      return  Legendre6(x) * Legendre3xx(y);
    }

    static double leg_quad_l6_l4(double x, double y)
    {
      return Legendre6(x) * Legendre4(y);
    }

    static double leg_quad_l6_l4x(double x, double y)
    {
      return Legendre6x(x) * Legendre4(y);
    }

    static double leg_quad_l6_l4y(double x, double y)
    {
      return Legendre6(x) * Legendre4x(y);
    }

    static double leg_quad_l6_l4xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre4(y);
    }

    static double leg_quad_l6_l4xy(double x, double y)
    {
      return Legendre6x(x) * Legendre4x(y);
    }

    static double leg_quad_l6_l4yy(double x, double y)
    {
      return  Legendre6(x) * Legendre4xx(y);
    }

    static double leg_quad_l6_l5(double x, double y)
    {
      return Legendre6(x) * Legendre5(y);
    }

    static double leg_quad_l6_l5x(double x, double y)
    {
      return Legendre6x(x) * Legendre5(y);
    }

    static double leg_quad_l6_l5y(double x, double y)
    {
      return Legendre6(x) * Legendre5x(y);
    }

    static double leg_quad_l6_l5xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre5(y);
    }

    static double leg_quad_l6_l5xy(double x, double y)
    {
      return Legendre6x(x) * Legendre5x(y);
    }

    static double leg_quad_l6_l5yy(double x, double y)
    {
      return  Legendre6(x) * Legendre5xx(y);
    }

    static double leg_quad_l6_l6(double x, double y)
    {
      return Legendre6(x) * Legendre6(y);
    }

    static double leg_quad_l6_l6x(double x, double y)
    {
      return Legendre6x(x) * Legendre6(y);
    }

    static double leg_quad_l6_l6y(double x, double y)
    {
      return Legendre6(x) * Legendre6x(y);
    }

    static double leg_quad_l6_l6xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre6(y);
    }

    static double leg_quad_l6_l6xy(double x, double y)
    {
      return Legendre6x(x) * Legendre6x(y);
    }

    static double leg_quad_l6_l6yy(double x, double y)
    {
      return  Legendre6(x) * Legendre6xx(y);
    }

    static double leg_quad_l6_l7(double x, double y)
    {
      return Legendre6(x) * Legendre7(y);
    }

    static double leg_quad_l6_l7x(double x, double y)
    {
      return Legendre6x(x) * Legendre7(y);
    }

    static double leg_quad_l6_l7y(double x, double y)
    {
      return Legendre6(x) * Legendre7x(y);
    }

    static double leg_quad_l6_l7xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre7(y);
    }

    static double leg_quad_l6_l7xy(double x, double y)
    {
      return Legendre6x(x) * Legendre7x(y);
    }

    static double leg_quad_l6_l7yy(double x, double y)
    {
      return  Legendre6(x) * Legendre7xx(y);
    }

    static double leg_quad_l6_l8(double x, double y)
    {
      return Legendre6(x) * Legendre8(y);
    }

    static double leg_quad_l6_l8x(double x, double y)
    {
      return Legendre6x(x) * Legendre8(y);
    }

    static double leg_quad_l6_l8y(double x, double y)
    {
      return Legendre6(x) * Legendre8x(y);
    }

    static double leg_quad_l6_l8xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre8(y);
    }

    static double leg_quad_l6_l8xy(double x, double y)
    {
      return Legendre6x(x) * Legendre8x(y);
    }

    static double leg_quad_l6_l8yy(double x, double y)
    {
      return  Legendre6(x) * Legendre8xx(y);
    }

    static double leg_quad_l6_l9(double x, double y)
    {
      return Legendre6(x) * Legendre9(y);
    }

    static double leg_quad_l6_l9x(double x, double y)
    {
      return Legendre6x(x) * Legendre9(y);
    }

    static double leg_quad_l6_l9y(double x, double y)
    {
      return Legendre6(x) * Legendre9x(y);
    }

    static double leg_quad_l6_l9xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre9(y);
    }

    static double leg_quad_l6_l9xy(double x, double y)
    {
      return Legendre6x(x) * Legendre9x(y);
    }

    static double leg_quad_l6_l9yy(double x, double y)
    {
      return  Legendre6(x) * Legendre9xx(y);
    }

    static double leg_quad_l6_l10(double x, double y)
    {
      return Legendre6(x) * Legendre10(y);
    }

    static double leg_quad_l6_l10x(double x, double y)
    {
      return Legendre6x(x) * Legendre10(y);
    }

    static double leg_quad_l6_l10y(double x, double y)
    {
      return Legendre6(x) * Legendre10x(y);
    }

    static double leg_quad_l6_l10xx(double x, double y)
    {
      return Legendre6xx(x) * Legendre10(y);
    }

    static double leg_quad_l6_l10xy(double x, double y)
    {
      return Legendre6x(x) * Legendre10x(y);
    }

    static double leg_quad_l6_l10yy(double x, double y)
    {
      return  Legendre6(x) * Legendre10xx(y);
    }

    static double leg_quad_l7_l0(double x, double y)
    {
      return Legendre7(x) * Legendre0(y);
    }

    static double leg_quad_l7_l0x(double x, double y)
    {
      return Legendre7x(x) * Legendre0(y);
    }

    static double leg_quad_l7_l0y(double x, double y)
    {
      return Legendre7(x) * Legendre0x(y);
    }

    static double leg_quad_l7_l0xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre0(y);
    }

    static double leg_quad_l7_l0xy(double x, double y)
    {
      return Legendre7x(x) * Legendre0x(y);
    }

    static double leg_quad_l7_l0yy(double x, double y)
    {
      return  Legendre7(x) * Legendre0xx(y);
    }

    static double leg_quad_l7_l1(double x, double y)
    {
      return Legendre7(x) * Legendre1(y);
    }

    static double leg_quad_l7_l1x(double x, double y)
    {
      return Legendre7x(x) * Legendre1(y);
    }

    static double leg_quad_l7_l1y(double x, double y)
    {
      return Legendre7(x) * Legendre1x(y);
    }

    static double leg_quad_l7_l1xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre1(y);
    }

    static double leg_quad_l7_l1xy(double x, double y)
    {
      return Legendre7x(x) * Legendre1x(y);
    }

    static double leg_quad_l7_l1yy(double x, double y)
    {
      return  Legendre7(x) * Legendre1xx(y);
    }

    static double leg_quad_l7_l2(double x, double y)
    {
      return Legendre7(x) * Legendre2(y);
    }

    static double leg_quad_l7_l2x(double x, double y)
    {
      return Legendre7x(x) * Legendre2(y);
    }

    static double leg_quad_l7_l2y(double x, double y)
    {
      return Legendre7(x) * Legendre2x(y);
    }

    static double leg_quad_l7_l2xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre2(y);
    }

    static double leg_quad_l7_l2xy(double x, double y)
    {
      return Legendre7x(x) * Legendre2x(y);
    }

    static double leg_quad_l7_l2yy(double x, double y)
    {
      return  Legendre7(x) * Legendre2xx(y);
    }

    static double leg_quad_l7_l3(double x, double y)
    {
      return Legendre7(x) * Legendre3(y);
    }

    static double leg_quad_l7_l3x(double x, double y)
    {
      return Legendre7x(x) * Legendre3(y);
    }

    static double leg_quad_l7_l3y(double x, double y)
    {
      return Legendre7(x) * Legendre3x(y);
    }

    static double leg_quad_l7_l3xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre3(y);
    }

    static double leg_quad_l7_l3xy(double x, double y)
    {
      return Legendre7x(x) * Legendre3x(y);
    }

    static double leg_quad_l7_l3yy(double x, double y)
    {
      return  Legendre7(x) * Legendre3xx(y);
    }

    static double leg_quad_l7_l4(double x, double y)
    {
      return Legendre7(x) * Legendre4(y);
    }

    static double leg_quad_l7_l4x(double x, double y)
    {
      return Legendre7x(x) * Legendre4(y);
    }

    static double leg_quad_l7_l4y(double x, double y)
    {
      return Legendre7(x) * Legendre4x(y);
    }

    static double leg_quad_l7_l4xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre4(y);
    }

    static double leg_quad_l7_l4xy(double x, double y)
    {
      return Legendre7x(x) * Legendre4x(y);
    }

    static double leg_quad_l7_l4yy(double x, double y)
    {
      return  Legendre7(x) * Legendre4xx(y);
    }

    static double leg_quad_l7_l5(double x, double y)
    {
      return Legendre7(x) * Legendre5(y);
    }

    static double leg_quad_l7_l5x(double x, double y)
    {
      return Legendre7x(x) * Legendre5(y);
    }

    static double leg_quad_l7_l5y(double x, double y)
    {
      return Legendre7(x) * Legendre5x(y);
    }

    static double leg_quad_l7_l5xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre5(y);
    }

    static double leg_quad_l7_l5xy(double x, double y)
    {
      return Legendre7x(x) * Legendre5x(y);
    }

    static double leg_quad_l7_l5yy(double x, double y)
    {
      return  Legendre7(x) * Legendre5xx(y);
    }

    static double leg_quad_l7_l6(double x, double y)
    {
      return Legendre7(x) * Legendre6(y);
    }

    static double leg_quad_l7_l6x(double x, double y)
    {
      return Legendre7x(x) * Legendre6(y);
    }

    static double leg_quad_l7_l6y(double x, double y)
    {
      return Legendre7(x) * Legendre6x(y);
    }

    static double leg_quad_l7_l6xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre6(y);
    }

    static double leg_quad_l7_l6xy(double x, double y)
    {
      return Legendre7x(x) * Legendre6x(y);
    }

    static double leg_quad_l7_l6yy(double x, double y)
    {
      return  Legendre7(x) * Legendre6xx(y);
    }

    static double leg_quad_l7_l7(double x, double y)
    {
      return Legendre7(x) * Legendre7(y);
    }

    static double leg_quad_l7_l7x(double x, double y)
    {
      return Legendre7x(x) * Legendre7(y);
    }

    static double leg_quad_l7_l7y(double x, double y)
    {
      return Legendre7(x) * Legendre7x(y);
    }

    static double leg_quad_l7_l7xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre7(y);
    }

    static double leg_quad_l7_l7xy(double x, double y)
    {
      return Legendre7x(x) * Legendre7x(y);
    }

    static double leg_quad_l7_l7yy(double x, double y)
    {
      return  Legendre7(x) * Legendre7xx(y);
    }

    static double leg_quad_l7_l8(double x, double y)
    {
      return Legendre7(x) * Legendre8(y);
    }

    static double leg_quad_l7_l8x(double x, double y)
    {
      return Legendre7x(x) * Legendre8(y);
    }

    static double leg_quad_l7_l8y(double x, double y)
    {
      return Legendre7(x) * Legendre8x(y);
    }

    static double leg_quad_l7_l8xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre8(y);
    }

    static double leg_quad_l7_l8xy(double x, double y)
    {
      return Legendre7x(x) * Legendre8x(y);
    }

    static double leg_quad_l7_l8yy(double x, double y)
    {
      return  Legendre7(x) * Legendre8xx(y);
    }

    static double leg_quad_l7_l9(double x, double y)
    {
      return Legendre7(x) * Legendre9(y);
    }

    static double leg_quad_l7_l9x(double x, double y)
    {
      return Legendre7x(x) * Legendre9(y);
    }

    static double leg_quad_l7_l9y(double x, double y)
    {
      return Legendre7(x) * Legendre9x(y);
    }

    static double leg_quad_l7_l9xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre9(y);
    }

    static double leg_quad_l7_l9xy(double x, double y)
    {
      return Legendre7x(x) * Legendre9x(y);
    }

    static double leg_quad_l7_l9yy(double x, double y)
    {
      return  Legendre7(x) * Legendre9xx(y);
    }

    static double leg_quad_l7_l10(double x, double y)
    {
      return Legendre7(x) * Legendre10(y);
    }

    static double leg_quad_l7_l10x(double x, double y)
    {
      return Legendre7x(x) * Legendre10(y);
    }

    static double leg_quad_l7_l10y(double x, double y)
    {
      return Legendre7(x) * Legendre10x(y);
    }

    static double leg_quad_l7_l10xx(double x, double y)
    {
      return Legendre7xx(x) * Legendre10(y);
    }

    static double leg_quad_l7_l10xy(double x, double y)
    {
      return Legendre7x(x) * Legendre10x(y);
    }

    static double leg_quad_l7_l10yy(double x, double y)
    {
      return  Legendre7(x) * Legendre10xx(y);
    }

    static double leg_quad_l8_l0(double x, double y)
    {
      return Legendre8(x) * Legendre0(y);
    }

    static double leg_quad_l8_l0x(double x, double y)
    {
      return Legendre8x(x) * Legendre0(y);
    }

    static double leg_quad_l8_l0y(double x, double y)
    {
      return Legendre8(x) * Legendre0x(y);
    }

    static double leg_quad_l8_l0xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre0(y);
    }

    static double leg_quad_l8_l0xy(double x, double y)
    {
      return Legendre8x(x) * Legendre0x(y);
    }

    static double leg_quad_l8_l0yy(double x, double y)
    {
      return  Legendre8(x) * Legendre0xx(y);
    }

    static double leg_quad_l8_l1(double x, double y)
    {
      return Legendre8(x) * Legendre1(y);
    }

    static double leg_quad_l8_l1x(double x, double y)
    {
      return Legendre8x(x) * Legendre1(y);
    }

    static double leg_quad_l8_l1y(double x, double y)
    {
      return Legendre8(x) * Legendre1x(y);
    }

    static double leg_quad_l8_l1xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre1(y);
    }

    static double leg_quad_l8_l1xy(double x, double y)
    {
      return Legendre8x(x) * Legendre1x(y);
    }

    static double leg_quad_l8_l1yy(double x, double y)
    {
      return  Legendre8(x) * Legendre1xx(y);
    }

    static double leg_quad_l8_l2(double x, double y)
    {
      return Legendre8(x) * Legendre2(y);
    }

    static double leg_quad_l8_l2x(double x, double y)
    {
      return Legendre8x(x) * Legendre2(y);
    }

    static double leg_quad_l8_l2y(double x, double y)
    {
      return Legendre8(x) * Legendre2x(y);
    }

    static double leg_quad_l8_l2xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre2(y);
    }

    static double leg_quad_l8_l2xy(double x, double y)
    {
      return Legendre8x(x) * Legendre2x(y);
    }

    static double leg_quad_l8_l2yy(double x, double y)
    {
      return  Legendre8(x) * Legendre2xx(y);
    }

    static double leg_quad_l8_l3(double x, double y)
    {
      return Legendre8(x) * Legendre3(y);
    }

    static double leg_quad_l8_l3x(double x, double y)
    {
      return Legendre8x(x) * Legendre3(y);
    }

    static double leg_quad_l8_l3y(double x, double y)
    {
      return Legendre8(x) * Legendre3x(y);
    }

    static double leg_quad_l8_l3xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre3(y);
    }

    static double leg_quad_l8_l3xy(double x, double y)
    {
      return Legendre8x(x) * Legendre3x(y);
    }

    static double leg_quad_l8_l3yy(double x, double y)
    {
      return  Legendre8(x) * Legendre3xx(y);
    }

    static double leg_quad_l8_l4(double x, double y)
    {
      return Legendre8(x) * Legendre4(y);
    }

    static double leg_quad_l8_l4x(double x, double y)
    {
      return Legendre8x(x) * Legendre4(y);
    }

    static double leg_quad_l8_l4y(double x, double y)
    {
      return Legendre8(x) * Legendre4x(y);
    }

    static double leg_quad_l8_l4xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre4(y);
    }

    static double leg_quad_l8_l4xy(double x, double y)
    {
      return Legendre8x(x) * Legendre4x(y);
    }

    static double leg_quad_l8_l4yy(double x, double y)
    {
      return  Legendre8(x) * Legendre4xx(y);
    }

    static double leg_quad_l8_l5(double x, double y)
    {
      return Legendre8(x) * Legendre5(y);
    }

    static double leg_quad_l8_l5x(double x, double y)
    {
      return Legendre8x(x) * Legendre5(y);
    }

    static double leg_quad_l8_l5y(double x, double y)
    {
      return Legendre8(x) * Legendre5x(y);
    }

    static double leg_quad_l8_l5xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre5(y);
    }

    static double leg_quad_l8_l5xy(double x, double y)
    {
      return Legendre8x(x) * Legendre5x(y);
    }

    static double leg_quad_l8_l5yy(double x, double y)
    {
      return  Legendre8(x) * Legendre5xx(y);
    }

    static double leg_quad_l8_l6(double x, double y)
    {
      return Legendre8(x) * Legendre6(y);
    }

    static double leg_quad_l8_l6x(double x, double y)
    {
      return Legendre8x(x) * Legendre6(y);
    }

    static double leg_quad_l8_l6y(double x, double y)
    {
      return Legendre8(x) * Legendre6x(y);
    }

    static double leg_quad_l8_l6xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre6(y);
    }

    static double leg_quad_l8_l6xy(double x, double y)
    {
      return Legendre8x(x) * Legendre6x(y);
    }

    static double leg_quad_l8_l6yy(double x, double y)
    {
      return  Legendre8(x) * Legendre6xx(y);
    }

    static double leg_quad_l8_l7(double x, double y)
    {
      return Legendre8(x) * Legendre7(y);
    }

    static double leg_quad_l8_l7x(double x, double y)
    {
      return Legendre8x(x) * Legendre7(y);
    }

    static double leg_quad_l8_l7y(double x, double y)
    {
      return Legendre8(x) * Legendre7x(y);
    }

    static double leg_quad_l8_l7xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre7(y);
    }

    static double leg_quad_l8_l7xy(double x, double y)
    {
      return Legendre8x(x) * Legendre7x(y);
    }

    static double leg_quad_l8_l7yy(double x, double y)
    {
      return  Legendre8(x) * Legendre7xx(y);
    }

    static double leg_quad_l8_l8(double x, double y)
    {
      return Legendre8(x) * Legendre8(y);
    }

    static double leg_quad_l8_l8x(double x, double y)
    {
      return Legendre8x(x) * Legendre8(y);
    }

    static double leg_quad_l8_l8y(double x, double y)
    {
      return Legendre8(x) * Legendre8x(y);
    }

    static double leg_quad_l8_l8xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre8(y);
    }

    static double leg_quad_l8_l8xy(double x, double y)
    {
      return Legendre8x(x) * Legendre8x(y);
    }

    static double leg_quad_l8_l8yy(double x, double y)
    {
      return  Legendre8(x) * Legendre8xx(y);
    }

    static double leg_quad_l8_l9(double x, double y)
    {
      return Legendre8(x) * Legendre9(y);
    }

    static double leg_quad_l8_l9x(double x, double y)
    {
      return Legendre8x(x) * Legendre9(y);
    }

    static double leg_quad_l8_l9y(double x, double y)
    {
      return Legendre8(x) * Legendre9x(y);
    }

    static double leg_quad_l8_l9xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre9(y);
    }

    static double leg_quad_l8_l9xy(double x, double y)
    {
      return Legendre8x(x) * Legendre9x(y);
    }

    static double leg_quad_l8_l9yy(double x, double y)
    {
      return  Legendre8(x) * Legendre9xx(y);
    }

    static double leg_quad_l8_l10(double x, double y)
    {
      return Legendre8(x) * Legendre10(y);
    }

    static double leg_quad_l8_l10x(double x, double y)
    {
      return Legendre8x(x) * Legendre10(y);
    }

    static double leg_quad_l8_l10y(double x, double y)
    {
      return Legendre8(x) * Legendre10x(y);
    }

    static double leg_quad_l8_l10xx(double x, double y)
    {
      return Legendre8xx(x) * Legendre10(y);
    }

    static double leg_quad_l8_l10xy(double x, double y)
    {
      return Legendre8x(x) * Legendre10x(y);
    }

    static double leg_quad_l8_l10yy(double x, double y)
    {
      return  Legendre8(x) * Legendre10xx(y);
    }

    static double leg_quad_l9_l0(double x, double y)
    {
      return Legendre9(x) * Legendre0(y);
    }

    static double leg_quad_l9_l0x(double x, double y)
    {
      return Legendre9x(x) * Legendre0(y);
    }

    static double leg_quad_l9_l0y(double x, double y)
    {
      return Legendre9(x) * Legendre0x(y);
    }

    static double leg_quad_l9_l0xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre0(y);
    }

    static double leg_quad_l9_l0xy(double x, double y)
    {
      return Legendre9x(x) * Legendre0x(y);
    }

    static double leg_quad_l9_l0yy(double x, double y)
    {
      return  Legendre9(x) * Legendre0xx(y);
    }

    static double leg_quad_l9_l1(double x, double y)
    {
      return Legendre9(x) * Legendre1(y);
    }

    static double leg_quad_l9_l1x(double x, double y)
    {
      return Legendre9x(x) * Legendre1(y);
    }

    static double leg_quad_l9_l1y(double x, double y)
    {
      return Legendre9(x) * Legendre1x(y);
    }

    static double leg_quad_l9_l1xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre1(y);
    }

    static double leg_quad_l9_l1xy(double x, double y)
    {
      return Legendre9x(x) * Legendre1x(y);
    }

    static double leg_quad_l9_l1yy(double x, double y)
    {
      return  Legendre9(x) * Legendre1xx(y);
    }

    static double leg_quad_l9_l2(double x, double y)
    {
      return Legendre9(x) * Legendre2(y);
    }

    static double leg_quad_l9_l2x(double x, double y)
    {
      return Legendre9x(x) * Legendre2(y);
    }

    static double leg_quad_l9_l2y(double x, double y)
    {
      return Legendre9(x) * Legendre2x(y);
    }

    static double leg_quad_l9_l2xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre2(y);
    }

    static double leg_quad_l9_l2xy(double x, double y)
    {
      return Legendre9x(x) * Legendre2x(y);
    }

    static double leg_quad_l9_l2yy(double x, double y)
    {
      return  Legendre9(x) * Legendre2xx(y);
    }

    static double leg_quad_l9_l3(double x, double y)
    {
      return Legendre9(x) * Legendre3(y);
    }

    static double leg_quad_l9_l3x(double x, double y)
    {
      return Legendre9x(x) * Legendre3(y);
    }

    static double leg_quad_l9_l3y(double x, double y)
    {
      return Legendre9(x) * Legendre3x(y);
    }

    static double leg_quad_l9_l3xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre3(y);
    }

    static double leg_quad_l9_l3xy(double x, double y)
    {
      return Legendre9x(x) * Legendre3x(y);
    }

    static double leg_quad_l9_l3yy(double x, double y)
    {
      return  Legendre9(x) * Legendre3xx(y);
    }

    static double leg_quad_l9_l4(double x, double y)
    {
      return Legendre9(x) * Legendre4(y);
    }

    static double leg_quad_l9_l4x(double x, double y)
    {
      return Legendre9x(x) * Legendre4(y);
    }

    static double leg_quad_l9_l4y(double x, double y)
    {
      return Legendre9(x) * Legendre4x(y);
    }

    static double leg_quad_l9_l4xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre4(y);
    }

    static double leg_quad_l9_l4xy(double x, double y)
    {
      return Legendre9x(x) * Legendre4x(y);
    }

    static double leg_quad_l9_l4yy(double x, double y)
    {
      return  Legendre9(x) * Legendre4xx(y);
    }

    static double leg_quad_l9_l5(double x, double y)
    {
      return Legendre9(x) * Legendre5(y);
    }

    static double leg_quad_l9_l5x(double x, double y)
    {
      return Legendre9x(x) * Legendre5(y);
    }

    static double leg_quad_l9_l5y(double x, double y)
    {
      return Legendre9(x) * Legendre5x(y);
    }

    static double leg_quad_l9_l5xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre5(y);
    }

    static double leg_quad_l9_l5xy(double x, double y)
    {
      return Legendre9x(x) * Legendre5x(y);
    }

    static double leg_quad_l9_l5yy(double x, double y)
    {
      return  Legendre9(x) * Legendre5xx(y);
    }

    static double leg_quad_l9_l6(double x, double y)
    {
      return Legendre9(x) * Legendre6(y);
    }

    static double leg_quad_l9_l6x(double x, double y)
    {
      return Legendre9x(x) * Legendre6(y);
    }

    static double leg_quad_l9_l6y(double x, double y)
    {
      return Legendre9(x) * Legendre6x(y);
    }

    static double leg_quad_l9_l6xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre6(y);
    }

    static double leg_quad_l9_l6xy(double x, double y)
    {
      return Legendre9x(x) * Legendre6x(y);
    }

    static double leg_quad_l9_l6yy(double x, double y)
    {
      return  Legendre9(x) * Legendre6xx(y);
    }

    static double leg_quad_l9_l7(double x, double y)
    {
      return Legendre9(x) * Legendre7(y);
    }

    static double leg_quad_l9_l7x(double x, double y)
    {
      return Legendre9x(x) * Legendre7(y);
    }

    static double leg_quad_l9_l7y(double x, double y)
    {
      return Legendre9(x) * Legendre7x(y);
    }

    static double leg_quad_l9_l7xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre7(y);
    }

    static double leg_quad_l9_l7xy(double x, double y)
    {
      return Legendre9x(x) * Legendre7x(y);
    }

    static double leg_quad_l9_l7yy(double x, double y)
    {
      return  Legendre9(x) * Legendre7xx(y);
    }

    static double leg_quad_l9_l8(double x, double y)
    {
      return Legendre9(x) * Legendre8(y);
    }

    static double leg_quad_l9_l8x(double x, double y)
    {
      return Legendre9x(x) * Legendre8(y);
    }

    static double leg_quad_l9_l8y(double x, double y)
    {
      return Legendre9(x) * Legendre8x(y);
    }

    static double leg_quad_l9_l8xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre8(y);
    }

    static double leg_quad_l9_l8xy(double x, double y)
    {
      return Legendre9x(x) * Legendre8x(y);
    }

    static double leg_quad_l9_l8yy(double x, double y)
    {
      return  Legendre9(x) * Legendre8xx(y);
    }

    static double leg_quad_l9_l9(double x, double y)
    {
      return Legendre9(x) * Legendre9(y);
    }

    static double leg_quad_l9_l9x(double x, double y)
    {
      return Legendre9x(x) * Legendre9(y);
    }

    static double leg_quad_l9_l9y(double x, double y)
    {
      return Legendre9(x) * Legendre9x(y);
    }

    static double leg_quad_l9_l9xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre9(y);
    }

    static double leg_quad_l9_l9xy(double x, double y)
    {
      return Legendre9x(x) * Legendre9x(y);
    }

    static double leg_quad_l9_l9yy(double x, double y)
    {
      return  Legendre9(x) * Legendre9xx(y);
    }

    static double leg_quad_l9_l10(double x, double y)
    {
      return Legendre9(x) * Legendre10(y);
    }

    static double leg_quad_l9_l10x(double x, double y)
    {
      return Legendre9x(x) * Legendre10(y);
    }

    static double leg_quad_l9_l10y(double x, double y)
    {
      return Legendre9(x) * Legendre10x(y);
    }

    static double leg_quad_l9_l10xx(double x, double y)
    {
      return Legendre9xx(x) * Legendre10(y);
    }

    static double leg_quad_l9_l10xy(double x, double y)
    {
      return Legendre9x(x) * Legendre10x(y);
    }

    static double leg_quad_l9_l10yy(double x, double y)
    {
      return  Legendre9(x) * Legendre10xx(y);
    }

    static double leg_quad_l10_l0(double x, double y)
    {
      return Legendre10(x) * Legendre0(y);
    }

    static double leg_quad_l10_l0x(double x, double y)
    {
      return Legendre10x(x) * Legendre0(y);
    }

    static double leg_quad_l10_l0y(double x, double y)
    {
      return Legendre10(x) * Legendre0x(y);
    }

    static double leg_quad_l10_l0xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre0(y);
    }

    static double leg_quad_l10_l0xy(double x, double y)
    {
      return Legendre10x(x) * Legendre0x(y);
    }

    static double leg_quad_l10_l0yy(double x, double y)
    {
      return  Legendre10(x) * Legendre0xx(y);
    }

    static double leg_quad_l10_l1(double x, double y)
    {
      return Legendre10(x) * Legendre1(y);
    }

    static double leg_quad_l10_l1x(double x, double y)
    {
      return Legendre10x(x) * Legendre1(y);
    }

    static double leg_quad_l10_l1y(double x, double y)
    {
      return Legendre10(x) * Legendre1x(y);
    }

    static double leg_quad_l10_l1xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre1(y);
    }

    static double leg_quad_l10_l1xy(double x, double y)
    {
      return Legendre10x(x) * Legendre1x(y);
    }

    static double leg_quad_l10_l1yy(double x, double y)
    {
      return  Legendre10(x) * Legendre1xx(y);
    }

    static double leg_quad_l10_l2(double x, double y)
    {
      return Legendre10(x) * Legendre2(y);
    }

    static double leg_quad_l10_l2x(double x, double y)
    {
      return Legendre10x(x) * Legendre2(y);
    }

    static double leg_quad_l10_l2y(double x, double y)
    {
      return Legendre10(x) * Legendre2x(y);
    }

    static double leg_quad_l10_l2xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre2(y);
    }

    static double leg_quad_l10_l2xy(double x, double y)
    {
      return Legendre10x(x) * Legendre2x(y);
    }

    static double leg_quad_l10_l2yy(double x, double y)
    {
      return  Legendre10(x) * Legendre2xx(y);
    }

    static double leg_quad_l10_l3(double x, double y)
    {
      return Legendre10(x) * Legendre3(y);
    }

    static double leg_quad_l10_l3x(double x, double y)
    {
      return Legendre10x(x) * Legendre3(y);
    }

    static double leg_quad_l10_l3y(double x, double y)
    {
      return Legendre10(x) * Legendre3x(y);
    }

    static double leg_quad_l10_l3xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre3(y);
    }

    static double leg_quad_l10_l3xy(double x, double y)
    {
      return Legendre10x(x) * Legendre3x(y);
    }

    static double leg_quad_l10_l3yy(double x, double y)
    {
      return  Legendre10(x) * Legendre3xx(y);
    }

    static double leg_quad_l10_l4(double x, double y)
    {
      return Legendre10(x) * Legendre4(y);
    }

    static double leg_quad_l10_l4x(double x, double y)
    {
      return Legendre10x(x) * Legendre4(y);
    }

    static double leg_quad_l10_l4y(double x, double y)
    {
      return Legendre10(x) * Legendre4x(y);
    }

    static double leg_quad_l10_l4xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre4(y);
    }

    static double leg_quad_l10_l4xy(double x, double y)
    {
      return Legendre10x(x) * Legendre4x(y);
    }

    static double leg_quad_l10_l4yy(double x, double y)
    {
      return  Legendre10(x) * Legendre4xx(y);
    }

    static double leg_quad_l10_l5(double x, double y)
    {
      return Legendre10(x) * Legendre5(y);
    }

    static double leg_quad_l10_l5x(double x, double y)
    {
      return Legendre10x(x) * Legendre5(y);
    }

    static double leg_quad_l10_l5y(double x, double y)
    {
      return Legendre10(x) * Legendre5x(y);
    }

    static double leg_quad_l10_l5xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre5(y);
    }

    static double leg_quad_l10_l5xy(double x, double y)
    {
      return Legendre10x(x) * Legendre5x(y);
    }

    static double leg_quad_l10_l5yy(double x, double y)
    {
      return  Legendre10(x) * Legendre5xx(y);
    }

    static double leg_quad_l10_l6(double x, double y)
    {
      return Legendre10(x) * Legendre6(y);
    }

    static double leg_quad_l10_l6x(double x, double y)
    {
      return Legendre10x(x) * Legendre6(y);
    }

    static double leg_quad_l10_l6y(double x, double y)
    {
      return Legendre10(x) * Legendre6x(y);
    }

    static double leg_quad_l10_l6xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre6(y);
    }

    static double leg_quad_l10_l6xy(double x, double y)
    {
      return Legendre10x(x) * Legendre6x(y);
    }

    static double leg_quad_l10_l6yy(double x, double y)
    {
      return  Legendre10(x) * Legendre6xx(y);
    }

    static double leg_quad_l10_l7(double x, double y)
    {
      return Legendre10(x) * Legendre7(y);
    }

    static double leg_quad_l10_l7x(double x, double y)
    {
      return Legendre10x(x) * Legendre7(y);
    }

    static double leg_quad_l10_l7y(double x, double y)
    {
      return Legendre10(x) * Legendre7x(y);
    }

    static double leg_quad_l10_l7xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre7(y);
    }

    static double leg_quad_l10_l7xy(double x, double y)
    {
      return Legendre10x(x) * Legendre7x(y);
    }

    static double leg_quad_l10_l7yy(double x, double y)
    {
      return  Legendre10(x) * Legendre7xx(y);
    }

    static double leg_quad_l10_l8(double x, double y)
    {
      return Legendre10(x) * Legendre8(y);
    }

    static double leg_quad_l10_l8x(double x, double y)
    {
      return Legendre10x(x) * Legendre8(y);
    }

    static double leg_quad_l10_l8y(double x, double y)
    {
      return Legendre10(x) * Legendre8x(y);
    }

    static double leg_quad_l10_l8xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre8(y);
    }

    static double leg_quad_l10_l8xy(double x, double y)
    {
      return Legendre10x(x) * Legendre8x(y);
    }

    static double leg_quad_l10_l8yy(double x, double y)
    {
      return  Legendre10(x) * Legendre8xx(y);
    }

    static double leg_quad_l10_l9(double x, double y)
    {
      return Legendre10(x) * Legendre9(y);
    }

    static double leg_quad_l10_l9x(double x, double y)
    {
      return Legendre10x(x) * Legendre9(y);
    }

    static double leg_quad_l10_l9y(double x, double y)
    {
      return Legendre10(x) * Legendre9x(y);
    }

    static double leg_quad_l10_l9xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre9(y);
    }

    static double leg_quad_l10_l9xy(double x, double y)
    {
      return Legendre10x(x) * Legendre9x(y);
    }

    static double leg_quad_l10_l9yy(double x, double y)
    {
      return  Legendre10(x) * Legendre9xx(y);
    }

    static double leg_quad_l10_l10(double x, double y)
    {
      return Legendre10(x) * Legendre10(y);
    }

    static double leg_quad_l10_l10x(double x, double y)
    {
      return Legendre10x(x) * Legendre10(y);
    }

    static double leg_quad_l10_l10y(double x, double y)
    {
      return Legendre10(x) * Legendre10x(y);
    }

    static double leg_quad_l10_l10xx(double x, double y)
    {
      return Legendre10xx(x) * Legendre10(y);
    }

    static double leg_quad_l10_l10xy(double x, double y)
    {
      return Legendre10x(x) * Legendre10x(y);
    }

    static double leg_quad_l10_l10yy(double x, double y)
    {
      return  Legendre10(x) * Legendre10xx(y);
    }

    static Shapeset::shape_fn_t leg_quad_fn[] =
    {
      leg_quad_l0_l0,   leg_quad_l0_l1,   leg_quad_l0_l2,   leg_quad_l0_l3,   leg_quad_l0_l4,
      leg_quad_l0_l5,   leg_quad_l0_l6,   leg_quad_l0_l7,   leg_quad_l0_l8,   leg_quad_l0_l9,
      leg_quad_l0_l10,   leg_quad_l1_l0,   leg_quad_l1_l1,   leg_quad_l1_l2,   leg_quad_l1_l3,
      leg_quad_l1_l4,   leg_quad_l1_l5,   leg_quad_l1_l6,   leg_quad_l1_l7,   leg_quad_l1_l8,
      leg_quad_l1_l9,   leg_quad_l1_l10,   leg_quad_l2_l0,   leg_quad_l2_l1,   leg_quad_l2_l2,
      leg_quad_l2_l3,   leg_quad_l2_l4,   leg_quad_l2_l5,   leg_quad_l2_l6,   leg_quad_l2_l7,
      leg_quad_l2_l8,   leg_quad_l2_l9,   leg_quad_l2_l10,   leg_quad_l3_l0,   leg_quad_l3_l1,
      leg_quad_l3_l2,   leg_quad_l3_l3,   leg_quad_l3_l4,   leg_quad_l3_l5,   leg_quad_l3_l6,
      leg_quad_l3_l7,   leg_quad_l3_l8,   leg_quad_l3_l9,   leg_quad_l3_l10,   leg_quad_l4_l0,
      leg_quad_l4_l1,   leg_quad_l4_l2,   leg_quad_l4_l3,   leg_quad_l4_l4,   leg_quad_l4_l5,
      leg_quad_l4_l6,   leg_quad_l4_l7,   leg_quad_l4_l8,   leg_quad_l4_l9,   leg_quad_l4_l10,
      leg_quad_l5_l0,   leg_quad_l5_l1,   leg_quad_l5_l2,   leg_quad_l5_l3,   leg_quad_l5_l4,
      leg_quad_l5_l5,   leg_quad_l5_l6,   leg_quad_l5_l7,   leg_quad_l5_l8,   leg_quad_l5_l9,
      leg_quad_l5_l10,   leg_quad_l6_l0,   leg_quad_l6_l1,   leg_quad_l6_l2,   leg_quad_l6_l3,
      leg_quad_l6_l4,   leg_quad_l6_l5,   leg_quad_l6_l6,   leg_quad_l6_l7,   leg_quad_l6_l8,
      leg_quad_l6_l9,   leg_quad_l6_l10,   leg_quad_l7_l0,   leg_quad_l7_l1,   leg_quad_l7_l2,
      leg_quad_l7_l3,   leg_quad_l7_l4,   leg_quad_l7_l5,   leg_quad_l7_l6,   leg_quad_l7_l7,
      leg_quad_l7_l8,   leg_quad_l7_l9,   leg_quad_l7_l10,   leg_quad_l8_l0,   leg_quad_l8_l1,
      leg_quad_l8_l2,   leg_quad_l8_l3,   leg_quad_l8_l4,   leg_quad_l8_l5,   leg_quad_l8_l6,
      leg_quad_l8_l7,   leg_quad_l8_l8,   leg_quad_l8_l9,   leg_quad_l8_l10,   leg_quad_l9_l0,
      leg_quad_l9_l1,   leg_quad_l9_l2,   leg_quad_l9_l3,   leg_quad_l9_l4,   leg_quad_l9_l5,
      leg_quad_l9_l6,   leg_quad_l9_l7,   leg_quad_l9_l8,   leg_quad_l9_l9,   leg_quad_l9_l10,
      leg_quad_l10_l0,   leg_quad_l10_l1,   leg_quad_l10_l2,   leg_quad_l10_l3,   leg_quad_l10_l4,
      leg_quad_l10_l5,   leg_quad_l10_l6,   leg_quad_l10_l7,   leg_quad_l10_l8,   leg_quad_l10_l9,
      leg_quad_l10_l10,
    };
    static Shapeset::shape_fn_t leg_quad_fn_dx[] =
    {
      leg_quad_l0_l0x,   leg_quad_l0_l1x,   leg_quad_l0_l2x,   leg_quad_l0_l3x,   leg_quad_l0_l4x,
      leg_quad_l0_l5x,   leg_quad_l0_l6x,   leg_quad_l0_l7x,   leg_quad_l0_l8x,   leg_quad_l0_l9x,
      leg_quad_l0_l10x,   leg_quad_l1_l0x,   leg_quad_l1_l1x,   leg_quad_l1_l2x,   leg_quad_l1_l3x,
      leg_quad_l1_l4x,   leg_quad_l1_l5x,   leg_quad_l1_l6x,   leg_quad_l1_l7x,   leg_quad_l1_l8x,
      leg_quad_l1_l9x,   leg_quad_l1_l10x,   leg_quad_l2_l0x,   leg_quad_l2_l1x,   leg_quad_l2_l2x,
      leg_quad_l2_l3x,   leg_quad_l2_l4x,   leg_quad_l2_l5x,   leg_quad_l2_l6x,   leg_quad_l2_l7x,
      leg_quad_l2_l8x,   leg_quad_l2_l9x,   leg_quad_l2_l10x,   leg_quad_l3_l0x,   leg_quad_l3_l1x,
      leg_quad_l3_l2x,   leg_quad_l3_l3x,   leg_quad_l3_l4x,   leg_quad_l3_l5x,   leg_quad_l3_l6x,
      leg_quad_l3_l7x,   leg_quad_l3_l8x,   leg_quad_l3_l9x,   leg_quad_l3_l10x,   leg_quad_l4_l0x,
      leg_quad_l4_l1x,   leg_quad_l4_l2x,   leg_quad_l4_l3x,   leg_quad_l4_l4x,   leg_quad_l4_l5x,
      leg_quad_l4_l6x,   leg_quad_l4_l7x,   leg_quad_l4_l8x,   leg_quad_l4_l9x,   leg_quad_l4_l10x,
      leg_quad_l5_l0x,   leg_quad_l5_l1x,   leg_quad_l5_l2x,   leg_quad_l5_l3x,   leg_quad_l5_l4x,
      leg_quad_l5_l5x,   leg_quad_l5_l6x,   leg_quad_l5_l7x,   leg_quad_l5_l8x,   leg_quad_l5_l9x,
      leg_quad_l5_l10x,   leg_quad_l6_l0x,   leg_quad_l6_l1x,   leg_quad_l6_l2x,   leg_quad_l6_l3x,
      leg_quad_l6_l4x,   leg_quad_l6_l5x,   leg_quad_l6_l6x,   leg_quad_l6_l7x,   leg_quad_l6_l8x,
      leg_quad_l6_l9x,   leg_quad_l6_l10x,   leg_quad_l7_l0x,   leg_quad_l7_l1x,   leg_quad_l7_l2x,
      leg_quad_l7_l3x,   leg_quad_l7_l4x,   leg_quad_l7_l5x,   leg_quad_l7_l6x,   leg_quad_l7_l7x,
      leg_quad_l7_l8x,   leg_quad_l7_l9x,   leg_quad_l7_l10x,   leg_quad_l8_l0x,   leg_quad_l8_l1x,
      leg_quad_l8_l2x,   leg_quad_l8_l3x,   leg_quad_l8_l4x,   leg_quad_l8_l5x,   leg_quad_l8_l6x,
      leg_quad_l8_l7x,   leg_quad_l8_l8x,   leg_quad_l8_l9x,   leg_quad_l8_l10x,   leg_quad_l9_l0x,
      leg_quad_l9_l1x,   leg_quad_l9_l2x,   leg_quad_l9_l3x,   leg_quad_l9_l4x,   leg_quad_l9_l5x,
      leg_quad_l9_l6x,   leg_quad_l9_l7x,   leg_quad_l9_l8x,   leg_quad_l9_l9x,   leg_quad_l9_l10x,
      leg_quad_l10_l0x,   leg_quad_l10_l1x,   leg_quad_l10_l2x,   leg_quad_l10_l3x,   leg_quad_l10_l4x,
      leg_quad_l10_l5x,   leg_quad_l10_l6x,   leg_quad_l10_l7x,   leg_quad_l10_l8x,   leg_quad_l10_l9x,
      leg_quad_l10_l10x,
    };
    static Shapeset::shape_fn_t leg_quad_fn_dy[] =
    {
      leg_quad_l0_l0y,   leg_quad_l0_l1y,   leg_quad_l0_l2y,   leg_quad_l0_l3y,   leg_quad_l0_l4y,
      leg_quad_l0_l5y,   leg_quad_l0_l6y,   leg_quad_l0_l7y,   leg_quad_l0_l8y,   leg_quad_l0_l9y,
      leg_quad_l0_l10y,   leg_quad_l1_l0y,   leg_quad_l1_l1y,   leg_quad_l1_l2y,   leg_quad_l1_l3y,
      leg_quad_l1_l4y,   leg_quad_l1_l5y,   leg_quad_l1_l6y,   leg_quad_l1_l7y,   leg_quad_l1_l8y,
      leg_quad_l1_l9y,   leg_quad_l1_l10y,   leg_quad_l2_l0y,   leg_quad_l2_l1y,   leg_quad_l2_l2y,
      leg_quad_l2_l3y,   leg_quad_l2_l4y,   leg_quad_l2_l5y,   leg_quad_l2_l6y,   leg_quad_l2_l7y,
      leg_quad_l2_l8y,   leg_quad_l2_l9y,   leg_quad_l2_l10y,   leg_quad_l3_l0y,   leg_quad_l3_l1y,
      leg_quad_l3_l2y,   leg_quad_l3_l3y,   leg_quad_l3_l4y,   leg_quad_l3_l5y,   leg_quad_l3_l6y,
      leg_quad_l3_l7y,   leg_quad_l3_l8y,   leg_quad_l3_l9y,   leg_quad_l3_l10y,   leg_quad_l4_l0y,
      leg_quad_l4_l1y,   leg_quad_l4_l2y,   leg_quad_l4_l3y,   leg_quad_l4_l4y,   leg_quad_l4_l5y,
      leg_quad_l4_l6y,   leg_quad_l4_l7y,   leg_quad_l4_l8y,   leg_quad_l4_l9y,   leg_quad_l4_l10y,
      leg_quad_l5_l0y,   leg_quad_l5_l1y,   leg_quad_l5_l2y,   leg_quad_l5_l3y,   leg_quad_l5_l4y,
      leg_quad_l5_l5y,   leg_quad_l5_l6y,   leg_quad_l5_l7y,   leg_quad_l5_l8y,   leg_quad_l5_l9y,
      leg_quad_l5_l10y,   leg_quad_l6_l0y,   leg_quad_l6_l1y,   leg_quad_l6_l2y,   leg_quad_l6_l3y,
      leg_quad_l6_l4y,   leg_quad_l6_l5y,   leg_quad_l6_l6y,   leg_quad_l6_l7y,   leg_quad_l6_l8y,
      leg_quad_l6_l9y,   leg_quad_l6_l10y,   leg_quad_l7_l0y,   leg_quad_l7_l1y,   leg_quad_l7_l2y,
      leg_quad_l7_l3y,   leg_quad_l7_l4y,   leg_quad_l7_l5y,   leg_quad_l7_l6y,   leg_quad_l7_l7y,
      leg_quad_l7_l8y,   leg_quad_l7_l9y,   leg_quad_l7_l10y,   leg_quad_l8_l0y,   leg_quad_l8_l1y,
      leg_quad_l8_l2y,   leg_quad_l8_l3y,   leg_quad_l8_l4y,   leg_quad_l8_l5y,   leg_quad_l8_l6y,
      leg_quad_l8_l7y,   leg_quad_l8_l8y,   leg_quad_l8_l9y,   leg_quad_l8_l10y,   leg_quad_l9_l0y,
      leg_quad_l9_l1y,   leg_quad_l9_l2y,   leg_quad_l9_l3y,   leg_quad_l9_l4y,   leg_quad_l9_l5y,
      leg_quad_l9_l6y,   leg_quad_l9_l7y,   leg_quad_l9_l8y,   leg_quad_l9_l9y,   leg_quad_l9_l10y,
      leg_quad_l10_l0y,   leg_quad_l10_l1y,   leg_quad_l10_l2y,   leg_quad_l10_l3y,   leg_quad_l10_l4y,
      leg_quad_l10_l5y,   leg_quad_l10_l6y,   leg_quad_l10_l7y,   leg_quad_l10_l8y,   leg_quad_l10_l9y,
      leg_quad_l10_l10y,
    };
    static Shapeset::shape_fn_t leg_quad_fn_dxx[] =
    {
      leg_quad_l0_l0xx,   leg_quad_l0_l1xx,   leg_quad_l0_l2xx,   leg_quad_l0_l3xx,   leg_quad_l0_l4xx,
      leg_quad_l0_l5xx,   leg_quad_l0_l6xx,   leg_quad_l0_l7xx,   leg_quad_l0_l8xx,   leg_quad_l0_l9xx,
      leg_quad_l0_l10xx,   leg_quad_l1_l0xx,   leg_quad_l1_l1xx,   leg_quad_l1_l2xx,   leg_quad_l1_l3xx,
      leg_quad_l1_l4xx,   leg_quad_l1_l5xx,   leg_quad_l1_l6xx,   leg_quad_l1_l7xx,   leg_quad_l1_l8xx,
      leg_quad_l1_l9xx,   leg_quad_l1_l10xx,   leg_quad_l2_l0xx,   leg_quad_l2_l1xx,   leg_quad_l2_l2xx,
      leg_quad_l2_l3xx,   leg_quad_l2_l4xx,   leg_quad_l2_l5xx,   leg_quad_l2_l6xx,   leg_quad_l2_l7xx,
      leg_quad_l2_l8xx,   leg_quad_l2_l9xx,   leg_quad_l2_l10xx,   leg_quad_l3_l0xx,   leg_quad_l3_l1xx,
      leg_quad_l3_l2xx,   leg_quad_l3_l3xx,   leg_quad_l3_l4xx,   leg_quad_l3_l5xx,   leg_quad_l3_l6xx,
      leg_quad_l3_l7xx,   leg_quad_l3_l8xx,   leg_quad_l3_l9xx,   leg_quad_l3_l10xx,   leg_quad_l4_l0xx,
      leg_quad_l4_l1xx,   leg_quad_l4_l2xx,   leg_quad_l4_l3xx,   leg_quad_l4_l4xx,   leg_quad_l4_l5xx,
      leg_quad_l4_l6xx,   leg_quad_l4_l7xx,   leg_quad_l4_l8xx,   leg_quad_l4_l9xx,   leg_quad_l4_l10xx,
      leg_quad_l5_l0xx,   leg_quad_l5_l1xx,   leg_quad_l5_l2xx,   leg_quad_l5_l3xx,   leg_quad_l5_l4xx,
      leg_quad_l5_l5xx,   leg_quad_l5_l6xx,   leg_quad_l5_l7xx,   leg_quad_l5_l8xx,   leg_quad_l5_l9xx,
      leg_quad_l5_l10xx,   leg_quad_l6_l0xx,   leg_quad_l6_l1xx,   leg_quad_l6_l2xx,   leg_quad_l6_l3xx,
      leg_quad_l6_l4xx,   leg_quad_l6_l5xx,   leg_quad_l6_l6xx,   leg_quad_l6_l7xx,   leg_quad_l6_l8xx,
      leg_quad_l6_l9xx,   leg_quad_l6_l10xx,   leg_quad_l7_l0xx,   leg_quad_l7_l1xx,   leg_quad_l7_l2xx,
      leg_quad_l7_l3xx,   leg_quad_l7_l4xx,   leg_quad_l7_l5xx,   leg_quad_l7_l6xx,   leg_quad_l7_l7xx,
      leg_quad_l7_l8xx,   leg_quad_l7_l9xx,   leg_quad_l7_l10xx,   leg_quad_l8_l0xx,   leg_quad_l8_l1xx,
      leg_quad_l8_l2xx,   leg_quad_l8_l3xx,   leg_quad_l8_l4xx,   leg_quad_l8_l5xx,   leg_quad_l8_l6xx,
      leg_quad_l8_l7xx,   leg_quad_l8_l8xx,   leg_quad_l8_l9xx,   leg_quad_l8_l10xx,   leg_quad_l9_l0xx,
      leg_quad_l9_l1xx,   leg_quad_l9_l2xx,   leg_quad_l9_l3xx,   leg_quad_l9_l4xx,   leg_quad_l9_l5xx,
      leg_quad_l9_l6xx,   leg_quad_l9_l7xx,   leg_quad_l9_l8xx,   leg_quad_l9_l9xx,   leg_quad_l9_l10xx,
      leg_quad_l10_l0xx,   leg_quad_l10_l1xx,   leg_quad_l10_l2xx,   leg_quad_l10_l3xx,   leg_quad_l10_l4xx,
      leg_quad_l10_l5xx,   leg_quad_l10_l6xx,   leg_quad_l10_l7xx,   leg_quad_l10_l8xx,   leg_quad_l10_l9xx,
      leg_quad_l10_l10xx,
    };
    static Shapeset::shape_fn_t leg_quad_fn_dxy[] =
    {
      leg_quad_l0_l0xy,   leg_quad_l0_l1xy,   leg_quad_l0_l2xy,   leg_quad_l0_l3xy,   leg_quad_l0_l4xy,
      leg_quad_l0_l5xy,   leg_quad_l0_l6xy,   leg_quad_l0_l7xy,   leg_quad_l0_l8xy,   leg_quad_l0_l9xy,
      leg_quad_l0_l10xy,   leg_quad_l1_l0xy,   leg_quad_l1_l1xy,   leg_quad_l1_l2xy,   leg_quad_l1_l3xy,
      leg_quad_l1_l4xy,   leg_quad_l1_l5xy,   leg_quad_l1_l6xy,   leg_quad_l1_l7xy,   leg_quad_l1_l8xy,
      leg_quad_l1_l9xy,   leg_quad_l1_l10xy,   leg_quad_l2_l0xy,   leg_quad_l2_l1xy,   leg_quad_l2_l2xy,
      leg_quad_l2_l3xy,   leg_quad_l2_l4xy,   leg_quad_l2_l5xy,   leg_quad_l2_l6xy,   leg_quad_l2_l7xy,
      leg_quad_l2_l8xy,   leg_quad_l2_l9xy,   leg_quad_l2_l10xy,   leg_quad_l3_l0xy,   leg_quad_l3_l1xy,
      leg_quad_l3_l2xy,   leg_quad_l3_l3xy,   leg_quad_l3_l4xy,   leg_quad_l3_l5xy,   leg_quad_l3_l6xy,
      leg_quad_l3_l7xy,   leg_quad_l3_l8xy,   leg_quad_l3_l9xy,   leg_quad_l3_l10xy,   leg_quad_l4_l0xy,
      leg_quad_l4_l1xy,   leg_quad_l4_l2xy,   leg_quad_l4_l3xy,   leg_quad_l4_l4xy,   leg_quad_l4_l5xy,
      leg_quad_l4_l6xy,   leg_quad_l4_l7xy,   leg_quad_l4_l8xy,   leg_quad_l4_l9xy,   leg_quad_l4_l10xy,
      leg_quad_l5_l0xy,   leg_quad_l5_l1xy,   leg_quad_l5_l2xy,   leg_quad_l5_l3xy,   leg_quad_l5_l4xy,
      leg_quad_l5_l5xy,   leg_quad_l5_l6xy,   leg_quad_l5_l7xy,   leg_quad_l5_l8xy,   leg_quad_l5_l9xy,
      leg_quad_l5_l10xy,   leg_quad_l6_l0xy,   leg_quad_l6_l1xy,   leg_quad_l6_l2xy,   leg_quad_l6_l3xy,
      leg_quad_l6_l4xy,   leg_quad_l6_l5xy,   leg_quad_l6_l6xy,   leg_quad_l6_l7xy,   leg_quad_l6_l8xy,
      leg_quad_l6_l9xy,   leg_quad_l6_l10xy,   leg_quad_l7_l0xy,   leg_quad_l7_l1xy,   leg_quad_l7_l2xy,
      leg_quad_l7_l3xy,   leg_quad_l7_l4xy,   leg_quad_l7_l5xy,   leg_quad_l7_l6xy,   leg_quad_l7_l7xy,
      leg_quad_l7_l8xy,   leg_quad_l7_l9xy,   leg_quad_l7_l10xy,   leg_quad_l8_l0xy,   leg_quad_l8_l1xy,
      leg_quad_l8_l2xy,   leg_quad_l8_l3xy,   leg_quad_l8_l4xy,   leg_quad_l8_l5xy,   leg_quad_l8_l6xy,
      leg_quad_l8_l7xy,   leg_quad_l8_l8xy,   leg_quad_l8_l9xy,   leg_quad_l8_l10xy,   leg_quad_l9_l0xy,
      leg_quad_l9_l1xy,   leg_quad_l9_l2xy,   leg_quad_l9_l3xy,   leg_quad_l9_l4xy,   leg_quad_l9_l5xy,
      leg_quad_l9_l6xy,   leg_quad_l9_l7xy,   leg_quad_l9_l8xy,   leg_quad_l9_l9xy,   leg_quad_l9_l10xy,
      leg_quad_l10_l0xy,   leg_quad_l10_l1xy,   leg_quad_l10_l2xy,   leg_quad_l10_l3xy,   leg_quad_l10_l4xy,
      leg_quad_l10_l5xy,   leg_quad_l10_l6xy,   leg_quad_l10_l7xy,   leg_quad_l10_l8xy,   leg_quad_l10_l9xy,
      leg_quad_l10_l10xy,
    };
    static Shapeset::shape_fn_t leg_quad_fn_dyy[] =
    {
      leg_quad_l0_l0yy,   leg_quad_l0_l1yy,   leg_quad_l0_l2yy,   leg_quad_l0_l3yy,   leg_quad_l0_l4yy,
      leg_quad_l0_l5yy,   leg_quad_l0_l6yy,   leg_quad_l0_l7yy,   leg_quad_l0_l8yy,   leg_quad_l0_l9yy,
      leg_quad_l0_l10yy,   leg_quad_l1_l0yy,   leg_quad_l1_l1yy,   leg_quad_l1_l2yy,   leg_quad_l1_l3yy,
      leg_quad_l1_l4yy,   leg_quad_l1_l5yy,   leg_quad_l1_l6yy,   leg_quad_l1_l7yy,   leg_quad_l1_l8yy,
      leg_quad_l1_l9yy,   leg_quad_l1_l10yy,   leg_quad_l2_l0yy,   leg_quad_l2_l1yy,   leg_quad_l2_l2yy,
      leg_quad_l2_l3yy,   leg_quad_l2_l4yy,   leg_quad_l2_l5yy,   leg_quad_l2_l6yy,   leg_quad_l2_l7yy,
      leg_quad_l2_l8yy,   leg_quad_l2_l9yy,   leg_quad_l2_l10yy,   leg_quad_l3_l0yy,   leg_quad_l3_l1yy,
      leg_quad_l3_l2yy,   leg_quad_l3_l3yy,   leg_quad_l3_l4yy,   leg_quad_l3_l5yy,   leg_quad_l3_l6yy,
      leg_quad_l3_l7yy,   leg_quad_l3_l8yy,   leg_quad_l3_l9yy,   leg_quad_l3_l10yy,   leg_quad_l4_l0yy,
      leg_quad_l4_l1yy,   leg_quad_l4_l2yy,   leg_quad_l4_l3yy,   leg_quad_l4_l4yy,   leg_quad_l4_l5yy,
      leg_quad_l4_l6yy,   leg_quad_l4_l7yy,   leg_quad_l4_l8yy,   leg_quad_l4_l9yy,   leg_quad_l4_l10yy,
      leg_quad_l5_l0yy,   leg_quad_l5_l1yy,   leg_quad_l5_l2yy,   leg_quad_l5_l3yy,   leg_quad_l5_l4yy,
      leg_quad_l5_l5yy,   leg_quad_l5_l6yy,   leg_quad_l5_l7yy,   leg_quad_l5_l8yy,   leg_quad_l5_l9yy,
      leg_quad_l5_l10yy,   leg_quad_l6_l0yy,   leg_quad_l6_l1yy,   leg_quad_l6_l2yy,   leg_quad_l6_l3yy,
      leg_quad_l6_l4yy,   leg_quad_l6_l5yy,   leg_quad_l6_l6yy,   leg_quad_l6_l7yy,   leg_quad_l6_l8yy,
      leg_quad_l6_l9yy,   leg_quad_l6_l10yy,   leg_quad_l7_l0yy,   leg_quad_l7_l1yy,   leg_quad_l7_l2yy,
      leg_quad_l7_l3yy,   leg_quad_l7_l4yy,   leg_quad_l7_l5yy,   leg_quad_l7_l6yy,   leg_quad_l7_l7yy,
      leg_quad_l7_l8yy,   leg_quad_l7_l9yy,   leg_quad_l7_l10yy,   leg_quad_l8_l0yy,   leg_quad_l8_l1yy,
      leg_quad_l8_l2yy,   leg_quad_l8_l3yy,   leg_quad_l8_l4yy,   leg_quad_l8_l5yy,   leg_quad_l8_l6yy,
      leg_quad_l8_l7yy,   leg_quad_l8_l8yy,   leg_quad_l8_l9yy,   leg_quad_l8_l10yy,   leg_quad_l9_l0yy,
      leg_quad_l9_l1yy,   leg_quad_l9_l2yy,   leg_quad_l9_l3yy,   leg_quad_l9_l4yy,   leg_quad_l9_l5yy,
      leg_quad_l9_l6yy,   leg_quad_l9_l7yy,   leg_quad_l9_l8yy,   leg_quad_l9_l9yy,   leg_quad_l9_l10yy,
      leg_quad_l10_l0yy,   leg_quad_l10_l1yy,   leg_quad_l10_l2yy,   leg_quad_l10_l3yy,   leg_quad_l10_l4yy,
      leg_quad_l10_l5yy,   leg_quad_l10_l6yy,   leg_quad_l10_l7yy,   leg_quad_l10_l8yy,   leg_quad_l10_l9yy,
      leg_quad_l10_l10yy,
    };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table[1]     = { leg_quad_fn };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table_dx[1]  = { leg_quad_fn_dx };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table_dy[1]  = { leg_quad_fn_dy };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table_dxx[1] = { leg_quad_fn_dxx };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table_dxy[1] = { leg_quad_fn_dxy };
    Shapeset::shape_fn_t* leg_quad_shape_fn_table_dyy[1] = { leg_quad_fn_dyy };

    static int qb_0_0[] = { 0, };
    static int qb_0_1[] = { 0, 1, };
    static int qb_0_2[] = { 0, 1, 2, };
    static int qb_0_3[] = { 0, 1, 2, 3, };
    static int qb_0_4[] = { 0, 1, 2, 3, 4, };
    static int qb_0_5[] = { 0, 1, 2, 3, 4, 5, };
    static int qb_0_6[] = { 0, 1, 2, 3, 4, 5, 6, };
    static int qb_0_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, };
    static int qb_0_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, };
    static int qb_0_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, };
    static int qb_0_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, };
    static int qb_1_0[] = { 0, 11, };
    static int qb_1_1[] = { 0, 1, 11, 12, };
    static int qb_1_2[] = { 0, 1, 2, 11, 12, 13, };
    static int qb_1_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, };
    static int qb_1_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, };
    static int qb_1_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, };
    static int qb_1_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, };
    static int qb_1_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, };
    static int qb_1_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, };
    static int qb_1_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, };
    static int qb_1_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, };
    static int qb_2_0[] = { 0, 11, 22, };
    static int qb_2_1[] = { 0, 1, 11, 12, 22, 23, };
    static int qb_2_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, };
    static int qb_2_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, };
    static int qb_2_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, };
    static int qb_2_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, };
    static int qb_2_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, };
    static int qb_2_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, };
    static int qb_2_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, };
    static int qb_2_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, };
    static int qb_2_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, };
    static int qb_3_0[] = { 0, 11, 22, 33, };
    static int qb_3_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, };
    static int qb_3_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, };
    static int qb_3_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, };
    static int qb_3_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, };
    static int qb_3_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, };
    static int qb_3_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, };
    static int qb_3_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, };
    static int qb_3_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, };
    static int qb_3_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, };
    static int qb_3_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, };
    static int qb_4_0[] = { 0, 11, 22, 33, 44, };
    static int qb_4_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, };
    static int qb_4_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, };
    static int qb_4_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, };
    static int qb_4_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, };
    static int qb_4_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, };
    static int qb_4_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, };
    static int qb_4_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, };
    static int qb_4_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, };
    static int qb_4_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, };
    static int qb_4_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, };
    static int qb_5_0[] = { 0, 11, 22, 33, 44, 55, };
    static int qb_5_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, };
    static int qb_5_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, };
    static int qb_5_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, };
    static int qb_5_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, };
    static int qb_5_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, };
    static int qb_5_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, };
    static int qb_5_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, };
    static int qb_5_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, };
    static int qb_5_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, };
    static int qb_5_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, };
    static int qb_6_0[] = { 0, 11, 22, 33, 44, 55, 66, };
    static int qb_6_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, 66, 67, };
    static int qb_6_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, 66, 67, 68, };
    static int qb_6_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, 66, 67, 68, 69, };
    static int qb_6_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, };
    static int qb_6_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, };
    static int qb_6_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, };
    static int qb_6_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, };
    static int qb_6_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 72, 73, 74, };
    static int qb_6_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, };
    static int qb_6_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, };
    static int qb_7_0[] = { 0, 11, 22, 33, 44, 55, 66, 77, };
    static int qb_7_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, 66, 67, 77, 78, };
    static int qb_7_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, 66, 67, 68, 77, 78, 79, };
    static int qb_7_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, 66, 67, 68, 69, 77, 78, 79, 80, };
    static int qb_7_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, 77, 78, 79, 80, 81, };
    static int qb_7_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, 77, 78, 79, 80, 81, 82, };
    static int qb_7_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, 77, 78, 79, 80, 81, 82, 83, };
    static int qb_7_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 77, 78, 79, 80, 81, 82, 83, 84, };
    static int qb_7_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 72, 73, 74, 77, 78, 79, 80, 81, 82, 83, 84, 85, };
    static int qb_7_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, };
    static int qb_7_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, };
    static int qb_8_0[] = { 0, 11, 22, 33, 44, 55, 66, 77, 88, };
    static int qb_8_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, 66, 67, 77, 78, 88, 89, };
    static int qb_8_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, 66, 67, 68, 77, 78, 79, 88, 89, 90, };
    static int qb_8_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, 66, 67, 68, 69, 77, 78, 79, 80, 88, 89, 90, 91, };
    static int qb_8_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, 77, 78, 79, 80, 81, 88, 89, 90, 91, 92, };
    static int qb_8_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, 77, 78, 79, 80, 81, 82, 88, 89, 90, 91, 92, 93, };
    static int qb_8_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, 77, 78, 79, 80, 81, 82, 83, 88, 89, 90, 91, 92, 93, 94, };
    static int qb_8_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 77, 78, 79, 80, 81, 82, 83, 84, 88, 89, 90, 91, 92, 93, 94, 95, };
    static int qb_8_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 72, 73, 74, 77, 78, 79, 80, 81, 82, 83, 84, 85, 88, 89, 90, 91, 92, 93, 94, 95, 96, };
    static int qb_8_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, };
    static int qb_8_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, };
    static int qb_9_0[] = { 0, 11, 22, 33, 44, 55, 66, 77, 88, 99, };
    static int qb_9_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, 66, 67, 77, 78, 88, 89, 99, 100, };
    static int qb_9_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, 66, 67, 68, 77, 78, 79, 88, 89, 90, 99, 100, 101, };
    static int qb_9_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, 66, 67, 68, 69, 77, 78, 79, 80, 88, 89, 90, 91, 99, 100, 101, 102, };
    static int qb_9_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, 77, 78, 79, 80, 81, 88, 89, 90, 91, 92, 99, 100, 101, 102, 103, };
    static int qb_9_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, 77, 78, 79, 80, 81, 82, 88, 89, 90, 91, 92, 93, 99, 100, 101, 102, 103, 104, };
    static int qb_9_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, 77, 78, 79, 80, 81, 82, 83, 88, 89, 90, 91, 92, 93, 94, 99, 100, 101, 102, 103, 104, 105, };
    static int qb_9_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 77, 78, 79, 80, 81, 82, 83, 84, 88, 89, 90, 91, 92, 93, 94, 95, 99, 100, 101, 102, 103, 104, 105, 106, };
    static int qb_9_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 72, 73, 74, 77, 78, 79, 80, 81, 82, 83, 84, 85, 88, 89, 90, 91, 92, 93, 94, 95, 96, 99, 100, 101, 102, 103, 104, 105, 106, 107, };
    static int qb_9_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, };
    static int qb_9_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, };
    static int qb_10_0[] = { 0, 11, 22, 33, 44, 55, 66, 77, 88, 99, 110, };
    static int qb_10_1[] = { 0, 1, 11, 12, 22, 23, 33, 34, 44, 45, 55, 56, 66, 67, 77, 78, 88, 89, 99, 100, 110, 111, };
    static int qb_10_2[] = { 0, 1, 2, 11, 12, 13, 22, 23, 24, 33, 34, 35, 44, 45, 46, 55, 56, 57, 66, 67, 68, 77, 78, 79, 88, 89, 90, 99, 100, 101, 110, 111, 112, };
    static int qb_10_3[] = { 0, 1, 2, 3, 11, 12, 13, 14, 22, 23, 24, 25, 33, 34, 35, 36, 44, 45, 46, 47, 55, 56, 57, 58, 66, 67, 68, 69, 77, 78, 79, 80, 88, 89, 90, 91, 99, 100, 101, 102, 110, 111, 112, 113, };
    static int qb_10_4[] = { 0, 1, 2, 3, 4, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, 77, 78, 79, 80, 81, 88, 89, 90, 91, 92, 99, 100, 101, 102, 103, 110, 111, 112, 113, 114, };
    static int qb_10_5[] = { 0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, 77, 78, 79, 80, 81, 82, 88, 89, 90, 91, 92, 93, 99, 100, 101, 102, 103, 104, 110, 111, 112, 113, 114, 115, };
    static int qb_10_6[] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, 77, 78, 79, 80, 81, 82, 83, 88, 89, 90, 91, 92, 93, 94, 99, 100, 101, 102, 103, 104, 105, 110, 111, 112, 113, 114, 115, 116, };
    static int qb_10_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 77, 78, 79, 80, 81, 82, 83, 84, 88, 89, 90, 91, 92, 93, 94, 95, 99, 100, 101, 102, 103, 104, 105, 106, 110, 111, 112, 113, 114, 115, 116, 117, };
    static int qb_10_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 72, 73, 74, 77, 78, 79, 80, 81, 82, 83, 84, 85, 88, 89, 90, 91, 92, 93, 94, 95, 96, 99, 100, 101, 102, 103, 104, 105, 106, 107, 110, 111, 112, 113, 114, 115, 116, 117, 118, };
    static int qb_10_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, };
    static int qb_10_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, };

    #define NULL16 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL

    int* leg_quad_bubble_indices[] =
    {
      qb_0_0,  qb_0_1,  qb_0_2,  qb_0_3,  qb_0_4,  qb_0_5,  qb_0_6,  qb_0_7,  qb_0_8,  qb_0_9,  qb_0_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_1_0,  qb_1_1,  qb_1_2,  qb_1_3,  qb_1_4,  qb_1_5,  qb_1_6,  qb_1_7,  qb_1_8,  qb_1_9,  qb_1_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_2_0,  qb_2_1,  qb_2_2,  qb_2_3,  qb_2_4,  qb_2_5,  qb_2_6,  qb_2_7,  qb_2_8,  qb_2_9,  qb_2_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_3_0,  qb_3_1,  qb_3_2,  qb_3_3,  qb_3_4,  qb_3_5,  qb_3_6,  qb_3_7,  qb_3_8,  qb_3_9,  qb_3_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_4_0,  qb_4_1,  qb_4_2,  qb_4_3,  qb_4_4,  qb_4_5,  qb_4_6,  qb_4_7,  qb_4_8,  qb_4_9,  qb_4_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_5_0,  qb_5_1,  qb_5_2,  qb_5_3,  qb_5_4,  qb_5_5,  qb_5_6,  qb_5_7,  qb_5_8,  qb_5_9,  qb_5_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_6_0,  qb_6_1,  qb_6_2,  qb_6_3,  qb_6_4,  qb_6_5,  qb_6_6,  qb_6_7,  qb_6_8,  qb_6_9,  qb_6_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_7_0,  qb_7_1,  qb_7_2,  qb_7_3,  qb_7_4,  qb_7_5,  qb_7_6,  qb_7_7,  qb_7_8,  qb_7_9,  qb_7_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_8_0,  qb_8_1,  qb_8_2,  qb_8_3,  qb_8_4,  qb_8_5,  qb_8_6,  qb_8_7,  qb_8_8,  qb_8_9,  qb_8_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_9_0,  qb_9_1,  qb_9_2,  qb_9_3,  qb_9_4,  qb_9_5,  qb_9_6,  qb_9_7,  qb_9_8,  qb_9_9,  qb_9_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_10_0,  qb_10_1,  qb_10_2,  qb_10_3,  qb_10_4,  qb_10_5,  qb_10_6,  qb_10_7,  qb_10_8,  qb_10_9,  qb_10_10,   NULL, NULL, NULL, NULL, NULL, NULL16,
    };

    #define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    int leg_quad_bubble_count[] =
    {
      1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  0,  0,  0,  0,  0, zero16,
      2,  4,  6,  8,  10,  12,  14,  16,  18,  20,  22,  0,  0,  0,  0,  0, zero16,
      3,  6,  9,  12,  15,  18,  21,  24,  27,  30,  33,  0,  0,  0,  0,  0, zero16,
      4,  8,  12,  16,  20,  24,  28,  32,  36,  40,  44,  0,  0,  0,  0,  0, zero16,
      5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  0,  0,  0,  0,  0, zero16,
      6,  12,  18,  24,  30,  36,  42,  48,  54,  60,  66,  0,  0,  0,  0,  0, zero16,
      7,  14,  21,  28,  35,  42,  49,  56,  63,  70,  77,  0,  0,  0,  0,  0, zero16,
      8,  16,  24,  32,  40,  48,  56,  64,  72,  80,  88,  0,  0,  0,  0,  0, zero16,
      9,  18,  27,  36,  45,  54,  63,  72,  81,  90,  99,  0,  0,  0,  0,  0, zero16,
      10,  20,  30,  40,  50,  60,  70,  80,  90,  100,  110,  0,  0,  0,  0,  0, zero16,
      11,  22,  33,  44,  55,  66,  77,  88,  99,  110,  121,  0,  0,  0,  0,  0, zero16,
    };

    int leg_quad_vertex_indices[4] = { -1, -1, -1, -1 };

    static int leg_quad_edge_indices_0[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_quad_edge_indices_1[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_quad_edge_indices_2[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_quad_edge_indices_3[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };

    int* leg_quad_edge_indices[4] =
    {
      leg_quad_edge_indices_0,
      leg_quad_edge_indices_1,
      leg_quad_edge_indices_2,
      leg_quad_edge_indices_3
    };

    #define oo H2D_MAKE_QUAD_ORDER
    #define XX(a, b) oo(a, b), oo(a, b)

    int leg_quad_index_to_order[] =
    {
      oo(0, 0),   oo(0, 1),   oo(0, 2),   oo(0, 3),   oo(0, 4),   oo(0, 5),   oo(0, 6),   oo(0, 7),   oo(0, 8),   oo(0, 9),   oo(0, 10),
      oo(1, 0),   oo(1, 1),   oo(1, 2),   oo(1, 3),   oo(1, 4),   oo(1, 5),   oo(1, 6),   oo(1, 7),   oo(1, 8),   oo(1, 9),   oo(1, 10),
      oo(2, 0),   oo(2, 1),   oo(2, 2),   oo(2, 3),   oo(2, 4),   oo(2, 5),   oo(2, 6),   oo(2, 7),   oo(2, 8),   oo(2, 9),   oo(2, 10),
      oo(3, 0),   oo(3, 1),   oo(3, 2),   oo(3, 3),   oo(3, 4),   oo(3, 5),   oo(3, 6),   oo(3, 7),   oo(3, 8),   oo(3, 9),   oo(3, 10),
      oo(4, 0),   oo(4, 1),   oo(4, 2),   oo(4, 3),   oo(4, 4),   oo(4, 5),   oo(4, 6),   oo(4, 7),   oo(4, 8),   oo(4, 9),   oo(4, 10),
      oo(5, 0),   oo(5, 1),   oo(5, 2),   oo(5, 3),   oo(5, 4),   oo(5, 5),   oo(5, 6),   oo(5, 7),   oo(5, 8),   oo(5, 9),   oo(5, 10),
      oo(6, 0),   oo(6, 1),   oo(6, 2),   oo(6, 3),   oo(6, 4),   oo(6, 5),   oo(6, 6),   oo(6, 7),   oo(6, 8),   oo(6, 9),   oo(6, 10),
      oo(7, 0),   oo(7, 1),   oo(7, 2),   oo(7, 3),   oo(7, 4),   oo(7, 5),   oo(7, 6),   oo(7, 7),   oo(7, 8),   oo(7, 9),   oo(7, 10),
      oo(8, 0),   oo(8, 1),   oo(8, 2),   oo(8, 3),   oo(8, 4),   oo(8, 5),   oo(8, 6),   oo(8, 7),   oo(8, 8),   oo(8, 9),   oo(8, 10),
      oo(9, 0),   oo(9, 1),   oo(9, 2),   oo(9, 3),   oo(9, 4),   oo(9, 5),   oo(9, 6),   oo(9, 7),   oo(9, 8),   oo(9, 9),   oo(9, 10),
      oo(10, 0),   oo(10, 1),   oo(10, 2),   oo(10, 3),   oo(10, 4),   oo(10, 5),   oo(10, 6),   oo(10, 7),   oo(10, 8),   oo(10, 9),   oo(10, 10),
    };

    //// triangle legendre shapeset /////////////////////////////////////////////////////////////////

    static double leg_tri_l0_l0(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l0_l1(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre1(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l1x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre1x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre1x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l1y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre1x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre1x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l1_l0(double x, double y)
    {
      return Legendre1(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l1_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l1_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l0_l2(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre2(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l2x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre2x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre2x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l2y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre2x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre2x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l2_l0(double x, double y)
    {
      return Legendre2(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l2_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l2_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l0_l3(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre3(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l3x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre3x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre3x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l3y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre3x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre3x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l3_l0(double x, double y)
    {
      return Legendre3(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l3_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l3_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l0_l4(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre4(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l4x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre4x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre4x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l4y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre4x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre4x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l4_l0(double x, double y)
    {
      return Legendre4(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l4_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l4_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l0_l5(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre5(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l5x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre5x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre5x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l5y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre5x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre5x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l5_l0(double x, double y)
    {
      return Legendre5(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l5_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l5_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l0_l6(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre6(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l6x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre6x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre6x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l6y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre6x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre6x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l6_l0(double x, double y)
    {
      return Legendre6(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l6_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre6(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre6(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l6_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre6(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre6(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l0_l7(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre7(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l7x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre7x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre7x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l7y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre7x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre7x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l7_l0(double x, double y)
    {
      return Legendre7(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l7_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre7(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre7(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l7_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre7(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre7(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l0_l8(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre8(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l8x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre8x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre8x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l8y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre8x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre8x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l8_l0(double x, double y)
    {
      return Legendre8(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l8_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre8(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre8(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l8_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre8(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre8(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l0_l9(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre9(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l9x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre9x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre9x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l9y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre9x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre9x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l9_l0(double x, double y)
    {
      return Legendre9(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l9_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre9(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre9(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l9_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre9(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre9(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l0_l10(double x, double y)
    {
      return Legendre0(lambda3(x, y) - lambda2(x, y)) * Legendre10(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l0_l10x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1x = Legendre10x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2x = Legendre10x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l0_l10y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre0(l3 - l2), L1y = Legendre10x(l3 - l2);
      double L2 = Legendre0(l2 - l1), L2y = Legendre10x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l10_l0(double x, double y)
    {
      return Legendre10(lambda3(x, y) - lambda2(x, y)) * Legendre0(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l10_l0x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre10(l3 - l2), L1x = Legendre0x(l3 - l2);
      double L2 = Legendre10(l2 - l1), L2x = Legendre0x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l10_l0y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre10(l3 - l2), L1y = Legendre0x(l3 - l2);
      double L2 = Legendre10(l2 - l1), L2y = Legendre0x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l1_l1(double x, double y)
    {
      return Legendre1(lambda3(x, y) - lambda2(x, y)) * Legendre1(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l1_l1x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1x = Legendre1x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2x = Legendre1x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l1_l1y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1y = Legendre1x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2y = Legendre1x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l1_l2(double x, double y)
    {
      return Legendre1(lambda3(x, y) - lambda2(x, y)) * Legendre2(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l1_l2x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1x = Legendre2x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2x = Legendre2x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l1_l2y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1y = Legendre2x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2y = Legendre2x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l2_l1(double x, double y)
    {
      return Legendre2(lambda3(x, y) - lambda2(x, y)) * Legendre1(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l2_l1x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1x = Legendre1x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2x = Legendre1x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l2_l1y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1y = Legendre1x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2y = Legendre1x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l1_l3(double x, double y)
    {
      return Legendre1(lambda3(x, y) - lambda2(x, y)) * Legendre3(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l1_l3x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1x = Legendre3x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2x = Legendre3x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l1_l3y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1y = Legendre3x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2y = Legendre3x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l3_l1(double x, double y)
    {
      return Legendre3(lambda3(x, y) - lambda2(x, y)) * Legendre1(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l3_l1x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1x = Legendre1x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2x = Legendre1x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l3_l1y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1y = Legendre1x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2y = Legendre1x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l1_l4(double x, double y)
    {
      return Legendre1(lambda3(x, y) - lambda2(x, y)) * Legendre4(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l1_l4x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1x = Legendre4x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2x = Legendre4x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l1_l4y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1y = Legendre4x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2y = Legendre4x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l4_l1(double x, double y)
    {
      return Legendre4(lambda3(x, y) - lambda2(x, y)) * Legendre1(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l4_l1x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1x = Legendre1x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2x = Legendre1x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l4_l1y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1y = Legendre1x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2y = Legendre1x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l1_l5(double x, double y)
    {
      return Legendre1(lambda3(x, y) - lambda2(x, y)) * Legendre5(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l1_l5x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1x = Legendre5x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2x = Legendre5x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l1_l5y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1y = Legendre5x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2y = Legendre5x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l5_l1(double x, double y)
    {
      return Legendre5(lambda3(x, y) - lambda2(x, y)) * Legendre1(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l5_l1x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1x = Legendre1x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2x = Legendre1x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l5_l1y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1y = Legendre1x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2y = Legendre1x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l1_l6(double x, double y)
    {
      return Legendre1(lambda3(x, y) - lambda2(x, y)) * Legendre6(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l1_l6x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1x = Legendre6x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2x = Legendre6x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l1_l6y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1y = Legendre6x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2y = Legendre6x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l6_l1(double x, double y)
    {
      return Legendre6(lambda3(x, y) - lambda2(x, y)) * Legendre1(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l6_l1x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre6(l3 - l2), L1x = Legendre1x(l3 - l2);
      double L2 = Legendre6(l2 - l1), L2x = Legendre1x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l6_l1y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre6(l3 - l2), L1y = Legendre1x(l3 - l2);
      double L2 = Legendre6(l2 - l1), L2y = Legendre1x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l1_l7(double x, double y)
    {
      return Legendre1(lambda3(x, y) - lambda2(x, y)) * Legendre7(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l1_l7x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1x = Legendre7x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2x = Legendre7x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l1_l7y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1y = Legendre7x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2y = Legendre7x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l7_l1(double x, double y)
    {
      return Legendre7(lambda3(x, y) - lambda2(x, y)) * Legendre1(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l7_l1x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre7(l3 - l2), L1x = Legendre1x(l3 - l2);
      double L2 = Legendre7(l2 - l1), L2x = Legendre1x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l7_l1y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre7(l3 - l2), L1y = Legendre1x(l3 - l2);
      double L2 = Legendre7(l2 - l1), L2y = Legendre1x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l1_l8(double x, double y)
    {
      return Legendre1(lambda3(x, y) - lambda2(x, y)) * Legendre8(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l1_l8x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1x = Legendre8x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2x = Legendre8x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l1_l8y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1y = Legendre8x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2y = Legendre8x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l8_l1(double x, double y)
    {
      return Legendre8(lambda3(x, y) - lambda2(x, y)) * Legendre1(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l8_l1x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre8(l3 - l2), L1x = Legendre1x(l3 - l2);
      double L2 = Legendre8(l2 - l1), L2x = Legendre1x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l8_l1y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre8(l3 - l2), L1y = Legendre1x(l3 - l2);
      double L2 = Legendre8(l2 - l1), L2y = Legendre1x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l1_l9(double x, double y)
    {
      return Legendre1(lambda3(x, y) - lambda2(x, y)) * Legendre9(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l1_l9x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1x = Legendre9x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2x = Legendre9x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l1_l9y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre1(l3 - l2), L1y = Legendre9x(l3 - l2);
      double L2 = Legendre1(l2 - l1), L2y = Legendre9x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l9_l1(double x, double y)
    {
      return Legendre9(lambda3(x, y) - lambda2(x, y)) * Legendre1(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l9_l1x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre9(l3 - l2), L1x = Legendre1x(l3 - l2);
      double L2 = Legendre9(l2 - l1), L2x = Legendre1x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l9_l1y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre9(l3 - l2), L1y = Legendre1x(l3 - l2);
      double L2 = Legendre9(l2 - l1), L2y = Legendre1x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l2_l2(double x, double y)
    {
      return Legendre2(lambda3(x, y) - lambda2(x, y)) * Legendre2(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l2_l2x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1x = Legendre2x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2x = Legendre2x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l2_l2y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1y = Legendre2x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2y = Legendre2x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l2_l3(double x, double y)
    {
      return Legendre2(lambda3(x, y) - lambda2(x, y)) * Legendre3(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l2_l3x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1x = Legendre3x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2x = Legendre3x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l2_l3y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1y = Legendre3x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2y = Legendre3x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l3_l2(double x, double y)
    {
      return Legendre3(lambda3(x, y) - lambda2(x, y)) * Legendre2(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l3_l2x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1x = Legendre2x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2x = Legendre2x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l3_l2y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1y = Legendre2x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2y = Legendre2x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l2_l4(double x, double y)
    {
      return Legendre2(lambda3(x, y) - lambda2(x, y)) * Legendre4(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l2_l4x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1x = Legendre4x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2x = Legendre4x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l2_l4y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1y = Legendre4x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2y = Legendre4x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l4_l2(double x, double y)
    {
      return Legendre4(lambda3(x, y) - lambda2(x, y)) * Legendre2(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l4_l2x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1x = Legendre2x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2x = Legendre2x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l4_l2y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1y = Legendre2x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2y = Legendre2x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l2_l5(double x, double y)
    {
      return Legendre2(lambda3(x, y) - lambda2(x, y)) * Legendre5(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l2_l5x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1x = Legendre5x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2x = Legendre5x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l2_l5y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1y = Legendre5x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2y = Legendre5x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l5_l2(double x, double y)
    {
      return Legendre5(lambda3(x, y) - lambda2(x, y)) * Legendre2(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l5_l2x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1x = Legendre2x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2x = Legendre2x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l5_l2y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1y = Legendre2x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2y = Legendre2x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l2_l6(double x, double y)
    {
      return Legendre2(lambda3(x, y) - lambda2(x, y)) * Legendre6(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l2_l6x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1x = Legendre6x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2x = Legendre6x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l2_l6y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1y = Legendre6x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2y = Legendre6x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l6_l2(double x, double y)
    {
      return Legendre6(lambda3(x, y) - lambda2(x, y)) * Legendre2(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l6_l2x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre6(l3 - l2), L1x = Legendre2x(l3 - l2);
      double L2 = Legendre6(l2 - l1), L2x = Legendre2x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l6_l2y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre6(l3 - l2), L1y = Legendre2x(l3 - l2);
      double L2 = Legendre6(l2 - l1), L2y = Legendre2x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l2_l7(double x, double y)
    {
      return Legendre2(lambda3(x, y) - lambda2(x, y)) * Legendre7(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l2_l7x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1x = Legendre7x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2x = Legendre7x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l2_l7y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1y = Legendre7x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2y = Legendre7x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l7_l2(double x, double y)
    {
      return Legendre7(lambda3(x, y) - lambda2(x, y)) * Legendre2(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l7_l2x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre7(l3 - l2), L1x = Legendre2x(l3 - l2);
      double L2 = Legendre7(l2 - l1), L2x = Legendre2x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l7_l2y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre7(l3 - l2), L1y = Legendre2x(l3 - l2);
      double L2 = Legendre7(l2 - l1), L2y = Legendre2x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l2_l8(double x, double y)
    {
      return Legendre2(lambda3(x, y) - lambda2(x, y)) * Legendre8(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l2_l8x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1x = Legendre8x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2x = Legendre8x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l2_l8y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre2(l3 - l2), L1y = Legendre8x(l3 - l2);
      double L2 = Legendre2(l2 - l1), L2y = Legendre8x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l8_l2(double x, double y)
    {
      return Legendre8(lambda3(x, y) - lambda2(x, y)) * Legendre2(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l8_l2x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre8(l3 - l2), L1x = Legendre2x(l3 - l2);
      double L2 = Legendre8(l2 - l1), L2x = Legendre2x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l8_l2y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre8(l3 - l2), L1y = Legendre2x(l3 - l2);
      double L2 = Legendre8(l2 - l1), L2y = Legendre2x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l3_l3(double x, double y)
    {
      return Legendre3(lambda3(x, y) - lambda2(x, y)) * Legendre3(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l3_l3x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1x = Legendre3x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2x = Legendre3x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l3_l3y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1y = Legendre3x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2y = Legendre3x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l3_l4(double x, double y)
    {
      return Legendre3(lambda3(x, y) - lambda2(x, y)) * Legendre4(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l3_l4x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1x = Legendre4x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2x = Legendre4x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l3_l4y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1y = Legendre4x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2y = Legendre4x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l4_l3(double x, double y)
    {
      return Legendre4(lambda3(x, y) - lambda2(x, y)) * Legendre3(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l4_l3x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1x = Legendre3x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2x = Legendre3x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l4_l3y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1y = Legendre3x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2y = Legendre3x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l3_l5(double x, double y)
    {
      return Legendre3(lambda3(x, y) - lambda2(x, y)) * Legendre5(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l3_l5x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1x = Legendre5x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2x = Legendre5x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l3_l5y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1y = Legendre5x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2y = Legendre5x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l5_l3(double x, double y)
    {
      return Legendre5(lambda3(x, y) - lambda2(x, y)) * Legendre3(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l5_l3x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1x = Legendre3x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2x = Legendre3x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l5_l3y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1y = Legendre3x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2y = Legendre3x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l3_l6(double x, double y)
    {
      return Legendre3(lambda3(x, y) - lambda2(x, y)) * Legendre6(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l3_l6x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1x = Legendre6x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2x = Legendre6x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l3_l6y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1y = Legendre6x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2y = Legendre6x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l6_l3(double x, double y)
    {
      return Legendre6(lambda3(x, y) - lambda2(x, y)) * Legendre3(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l6_l3x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre6(l3 - l2), L1x = Legendre3x(l3 - l2);
      double L2 = Legendre6(l2 - l1), L2x = Legendre3x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l6_l3y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre6(l3 - l2), L1y = Legendre3x(l3 - l2);
      double L2 = Legendre6(l2 - l1), L2y = Legendre3x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l3_l7(double x, double y)
    {
      return Legendre3(lambda3(x, y) - lambda2(x, y)) * Legendre7(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l3_l7x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1x = Legendre7x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2x = Legendre7x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l3_l7y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre3(l3 - l2), L1y = Legendre7x(l3 - l2);
      double L2 = Legendre3(l2 - l1), L2y = Legendre7x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l7_l3(double x, double y)
    {
      return Legendre7(lambda3(x, y) - lambda2(x, y)) * Legendre3(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l7_l3x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre7(l3 - l2), L1x = Legendre3x(l3 - l2);
      double L2 = Legendre7(l2 - l1), L2x = Legendre3x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l7_l3y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre7(l3 - l2), L1y = Legendre3x(l3 - l2);
      double L2 = Legendre7(l2 - l1), L2y = Legendre3x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l4_l4(double x, double y)
    {
      return Legendre4(lambda3(x, y) - lambda2(x, y)) * Legendre4(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l4_l4x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1x = Legendre4x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2x = Legendre4x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l4_l4y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1y = Legendre4x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2y = Legendre4x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l4_l5(double x, double y)
    {
      return Legendre4(lambda3(x, y) - lambda2(x, y)) * Legendre5(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l4_l5x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1x = Legendre5x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2x = Legendre5x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l4_l5y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1y = Legendre5x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2y = Legendre5x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l5_l4(double x, double y)
    {
      return Legendre5(lambda3(x, y) - lambda2(x, y)) * Legendre4(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l5_l4x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1x = Legendre4x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2x = Legendre4x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l5_l4y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1y = Legendre4x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2y = Legendre4x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l4_l6(double x, double y)
    {
      return Legendre4(lambda3(x, y) - lambda2(x, y)) * Legendre6(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l4_l6x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1x = Legendre6x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2x = Legendre6x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l4_l6y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre4(l3 - l2), L1y = Legendre6x(l3 - l2);
      double L2 = Legendre4(l2 - l1), L2y = Legendre6x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l6_l4(double x, double y)
    {
      return Legendre6(lambda3(x, y) - lambda2(x, y)) * Legendre4(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l6_l4x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre6(l3 - l2), L1x = Legendre4x(l3 - l2);
      double L2 = Legendre6(l2 - l1), L2x = Legendre4x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l6_l4y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre6(l3 - l2), L1y = Legendre4x(l3 - l2);
      double L2 = Legendre6(l2 - l1), L2y = Legendre4x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static double leg_tri_l5_l5(double x, double y)
    {
      return Legendre5(lambda3(x, y) - lambda2(x, y)) * Legendre5(lambda2(x, y) - lambda1(x, y));
    }

    static double leg_tri_l5_l5x(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1x = Legendre5x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2x = Legendre5x(l2 - l1);
      return L1x * (lambda3x(x, y) - lambda2x(x, y)) * L2 + L1 * L2x * (lambda2x(x, y) - lambda1x(x, y));
    }

    static double leg_tri_l5_l5y(double x, double y)
    {
      double l1 = lambda1(x, y), l2 = lambda2(x, y), l3 = lambda3(x, y);
      double L1 = Legendre5(l3 - l2), L1y = Legendre5x(l3 - l2);
      double L2 = Legendre5(l2 - l1), L2y = Legendre5x(l2 - l1);
      return L1y * (lambda3y(x, y) - lambda2y(x, y)) * L2 + L1 * L2y * (lambda2y(x, y) - lambda1y(x, y));
    }

    static Shapeset::shape_fn_t leg_tri_fn[] =
    {
      leg_tri_l0_l0,   leg_tri_l0_l1,   leg_tri_l1_l0,   leg_tri_l0_l2,   leg_tri_l2_l0,
      leg_tri_l0_l3,   leg_tri_l3_l0,   leg_tri_l0_l4,   leg_tri_l4_l0,   leg_tri_l0_l5,
      leg_tri_l5_l0,   leg_tri_l0_l6,   leg_tri_l6_l0,   leg_tri_l0_l7,   leg_tri_l7_l0,
      leg_tri_l0_l8,   leg_tri_l8_l0,   leg_tri_l0_l9,   leg_tri_l9_l0,   leg_tri_l0_l10,
      leg_tri_l10_l0,   leg_tri_l1_l1,   leg_tri_l1_l2,   leg_tri_l2_l1,   leg_tri_l1_l3,
      leg_tri_l3_l1,   leg_tri_l1_l4,   leg_tri_l4_l1,   leg_tri_l1_l5,   leg_tri_l5_l1,
      leg_tri_l1_l6,   leg_tri_l6_l1,   leg_tri_l1_l7,   leg_tri_l7_l1,   leg_tri_l1_l8,
      leg_tri_l8_l1,   leg_tri_l1_l9,   leg_tri_l9_l1,   leg_tri_l2_l2,   leg_tri_l2_l3,
      leg_tri_l3_l2,   leg_tri_l2_l4,   leg_tri_l4_l2,   leg_tri_l2_l5,   leg_tri_l5_l2,
      leg_tri_l2_l6,   leg_tri_l6_l2,   leg_tri_l2_l7,   leg_tri_l7_l2,   leg_tri_l2_l8,
      leg_tri_l8_l2,   leg_tri_l3_l3,   leg_tri_l3_l4,   leg_tri_l4_l3,   leg_tri_l3_l5,
      leg_tri_l5_l3,   leg_tri_l3_l6,   leg_tri_l6_l3,   leg_tri_l3_l7,   leg_tri_l7_l3,
      leg_tri_l4_l4,   leg_tri_l4_l5,   leg_tri_l5_l4,   leg_tri_l4_l6,   leg_tri_l6_l4,
      leg_tri_l5_l5,
    };
    static Shapeset::shape_fn_t leg_tri_fn_dx[] =
    {
      leg_tri_l0_l0x,   leg_tri_l0_l1x,   leg_tri_l1_l0x,   leg_tri_l0_l2x,   leg_tri_l2_l0x,
      leg_tri_l0_l3x,   leg_tri_l3_l0x,   leg_tri_l0_l4x,   leg_tri_l4_l0x,   leg_tri_l0_l5x,
      leg_tri_l5_l0x,   leg_tri_l0_l6x,   leg_tri_l6_l0x,   leg_tri_l0_l7x,   leg_tri_l7_l0x,
      leg_tri_l0_l8x,   leg_tri_l8_l0x,   leg_tri_l0_l9x,   leg_tri_l9_l0x,   leg_tri_l0_l10x,
      leg_tri_l10_l0x,   leg_tri_l1_l1x,   leg_tri_l1_l2x,   leg_tri_l2_l1x,   leg_tri_l1_l3x,
      leg_tri_l3_l1x,   leg_tri_l1_l4x,   leg_tri_l4_l1x,   leg_tri_l1_l5x,   leg_tri_l5_l1x,
      leg_tri_l1_l6x,   leg_tri_l6_l1x,   leg_tri_l1_l7x,   leg_tri_l7_l1x,   leg_tri_l1_l8x,
      leg_tri_l8_l1x,   leg_tri_l1_l9x,   leg_tri_l9_l1x,   leg_tri_l2_l2x,   leg_tri_l2_l3x,
      leg_tri_l3_l2x,   leg_tri_l2_l4x,   leg_tri_l4_l2x,   leg_tri_l2_l5x,   leg_tri_l5_l2x,
      leg_tri_l2_l6x,   leg_tri_l6_l2x,   leg_tri_l2_l7x,   leg_tri_l7_l2x,   leg_tri_l2_l8x,
      leg_tri_l8_l2x,   leg_tri_l3_l3x,   leg_tri_l3_l4x,   leg_tri_l4_l3x,   leg_tri_l3_l5x,
      leg_tri_l5_l3x,   leg_tri_l3_l6x,   leg_tri_l6_l3x,   leg_tri_l3_l7x,   leg_tri_l7_l3x,
      leg_tri_l4_l4x,   leg_tri_l4_l5x,   leg_tri_l5_l4x,   leg_tri_l4_l6x,   leg_tri_l6_l4x,
      leg_tri_l5_l5x,
    };
    static Shapeset::shape_fn_t leg_tri_fn_dy[] =
    {
      leg_tri_l0_l0y,   leg_tri_l0_l1y,   leg_tri_l1_l0y,   leg_tri_l0_l2y,   leg_tri_l2_l0y,
      leg_tri_l0_l3y,   leg_tri_l3_l0y,   leg_tri_l0_l4y,   leg_tri_l4_l0y,   leg_tri_l0_l5y,
      leg_tri_l5_l0y,   leg_tri_l0_l6y,   leg_tri_l6_l0y,   leg_tri_l0_l7y,   leg_tri_l7_l0y,
      leg_tri_l0_l8y,   leg_tri_l8_l0y,   leg_tri_l0_l9y,   leg_tri_l9_l0y,   leg_tri_l0_l10y,
      leg_tri_l10_l0y,   leg_tri_l1_l1y,   leg_tri_l1_l2y,   leg_tri_l2_l1y,   leg_tri_l1_l3y,
      leg_tri_l3_l1y,   leg_tri_l1_l4y,   leg_tri_l4_l1y,   leg_tri_l1_l5y,   leg_tri_l5_l1y,
      leg_tri_l1_l6y,   leg_tri_l6_l1y,   leg_tri_l1_l7y,   leg_tri_l7_l1y,   leg_tri_l1_l8y,
      leg_tri_l8_l1y,   leg_tri_l1_l9y,   leg_tri_l9_l1y,   leg_tri_l2_l2y,   leg_tri_l2_l3y,
      leg_tri_l3_l2y,   leg_tri_l2_l4y,   leg_tri_l4_l2y,   leg_tri_l2_l5y,   leg_tri_l5_l2y,
      leg_tri_l2_l6y,   leg_tri_l6_l2y,   leg_tri_l2_l7y,   leg_tri_l7_l2y,   leg_tri_l2_l8y,
      leg_tri_l8_l2y,   leg_tri_l3_l3y,   leg_tri_l3_l4y,   leg_tri_l4_l3y,   leg_tri_l3_l5y,
      leg_tri_l5_l3y,   leg_tri_l3_l6y,   leg_tri_l6_l3y,   leg_tri_l3_l7y,   leg_tri_l7_l3y,
      leg_tri_l4_l4y,   leg_tri_l4_l5y,   leg_tri_l5_l4y,   leg_tri_l4_l6y,   leg_tri_l6_l4y,
      leg_tri_l5_l5y,
    };
    Shapeset::shape_fn_t* leg_tri_shape_fn_table[1]     = { leg_tri_fn };
    Shapeset::shape_fn_t* leg_tri_shape_fn_table_dx[1]  = { leg_tri_fn_dx };
    Shapeset::shape_fn_t* leg_tri_shape_fn_table_dy[1]  = { leg_tri_fn_dy };

    static int qb_0[] = { 0, };
    static int qb_1[] = { 0, 1, 2, };
    static int qb_2[] = { 0, 1, 2, 3, 4, 21, };
    static int qb_3[] = { 0, 1, 2, 3, 4, 5, 6, 21, 22, 23, };
    static int qb_4[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 21, 22, 23, 24, 25, 38, };
    static int qb_5[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 21, 22, 23, 24, 25, 26, 27, 38, 39, 40, };
    static int qb_6[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 26, 27, 28, 29, 38, 39, 40, 41, 42, 51, };
    static int qb_7[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 38, 39, 40, 41, 42, 43, 44, 51, 52, 53, };
    static int qb_8[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 38, 39, 40, 41, 42, 43, 44, 45, 46, 51, 52, 53, 54, 55, 60, };
    static int qb_9[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 51, 52, 53, 54, 55, 56, 57, 60, 61, 62, };
    static int qb_10[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, };

    int* leg_tri_bubble_indices[11] = {  qb_0,   qb_1,   qb_2,   qb_3,   qb_4,   qb_5,   qb_6,   qb_7,   qb_8,   qb_9,   qb_10 };

    int leg_tri_bubble_count[11] = { 1,  3,  6,  10,  15,  21,  28,  36,  45,  55,  66 };

    int leg_tri_vertex_indices[4] = { -1, -1, -1, -1 };

    static int leg_tri_edge_indices_0[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_tri_edge_indices_1[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_tri_edge_indices_2[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int leg_tri_edge_indices_3[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };

    int* leg_tri_edge_indices[4] =
    {
      leg_tri_edge_indices_0,
      leg_tri_edge_indices_1,
      leg_tri_edge_indices_2,
      leg_tri_edge_indices_3
    };

    int leg_tri_index_to_order[] =
    {
      0,   1,   1,   2,   2,   3,   3,   4,   4,   5,  5,   6,   6,   7,   7,   8,   8,   9,   9,   10,  10,
      2,   3,   3,   4,   4,   5,  5,   6,   6,   7,   7,   8,   8,   9,   9,   10,  10,
      4,   5,  5,   6,   6,   7,   7,   8,   8,   9,   9,   10,  10,
      6,   7,   7,   8,   8,   9,   9,   10,  10,
      8,   9,   9,   10,  10,
      10,
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Shapeset::shape_fn_t** leg_shape_fn_table[2] =
    {
      leg_tri_shape_fn_table,
      leg_quad_shape_fn_table
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_dx[2] =
    {
      leg_tri_shape_fn_table_dx,
      leg_quad_shape_fn_table_dx
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_dy[2] =
    {
      leg_tri_shape_fn_table_dy,
      leg_quad_shape_fn_table_dy
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_dxx[2] =
    {
      NULL,
      leg_quad_shape_fn_table_dxx
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_dyy[2] =
    {
      NULL,
      leg_quad_shape_fn_table_dyy
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_dxy[2] =
    {
      NULL,
      leg_quad_shape_fn_table_dxy
    };

    static int* leg_vertex_indices[2] =
    {
      leg_tri_vertex_indices,
      leg_quad_vertex_indices
    };

    static int** leg_edge_indices[2] =
    {
      leg_tri_edge_indices,
      leg_quad_edge_indices
    };

    static int** leg_bubble_indices[2] =
    {
      leg_tri_bubble_indices,
      leg_quad_bubble_indices
    };

    static int* leg_bubble_count[2] =
    {
      leg_tri_bubble_count,
      leg_quad_bubble_count
    };

    static int* leg_index_to_order[2] =
    {
      leg_tri_index_to_order,
      leg_quad_index_to_order
    };

    L2ShapesetLegendre::L2ShapesetLegendre()
    {
      shape_table[0] = leg_shape_fn_table;
      shape_table[1] = leg_shape_fn_table_dx;
      shape_table[2] = leg_shape_fn_table_dy;
      shape_table[3] = leg_shape_fn_table_dxx;
      shape_table[4] = leg_shape_fn_table_dyy;
      shape_table[5] = leg_shape_fn_table_dxy;

      vertex_indices = leg_vertex_indices;
      edge_indices = leg_edge_indices;
      bubble_indices = leg_bubble_indices;
      bubble_count = leg_bubble_count;
      index_to_order = leg_index_to_order;

      ref_vert[0][0][0] = -1.0;
      ref_vert[0][0][1] = -1.0;
      ref_vert[0][1][0] =  1.0;
      ref_vert[0][1][1] = -1.0;
      ref_vert[0][2][0] = -1.0;
      ref_vert[0][2][1] =  1.0;

      ref_vert[1][0][0] = -1.0;
      ref_vert[1][0][1] = -1.0;
      ref_vert[1][1][0] =  1.0;
      ref_vert[1][1][1] = -1.0;
      ref_vert[1][2][0] =  1.0;
      ref_vert[1][2][1] =  1.0;
      ref_vert[1][3][0] = -1.0;
      ref_vert[1][3][1] =  1.0;

      max_order = 10;
      num_components = 1;

      max_index[0] = 66;
      max_index[1] = 121;

      ebias = 2;

      comb_table = NULL;
    }
  }
}