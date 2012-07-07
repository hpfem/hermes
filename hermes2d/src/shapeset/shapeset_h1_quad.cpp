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

//// quad simple lobatto shapeset /////////////////////////////////////////////////////////////////
namespace Hermes
{
  namespace Hermes2D
  {
    static double simple_quad_l0_l0(double x, double y)
    {
      return   l0(x) * l0(y);
    }

    static double simple_quad_l0_l0x(double x, double y)
    {
      return   dl0(x) * l0(y);
    }

    static double simple_quad_l0_l0y(double x, double y)
    {
      return   l0(x) * dl0(y);
    }

    static double simple_quad_l0_l0xx(double x, double y)
    {
      return   d2l0(x) * l0(y);
    }

    static double simple_quad_l0_l0xy(double x, double y)
    {
      return   dl0(x) * dl0(y);
    }

    static double simple_quad_l0_l0yy(double x, double y)
    {
      return   l0(x) * d2l0(y);
    }

    static double simple_quad_l0_l1(double x, double y)
    {
      return   l0(x) * l1(y);
    }

    static double simple_quad_l0_l1x(double x, double y)
    {
      return   dl0(x) * l1(y);
    }

    static double simple_quad_l0_l1y(double x, double y)
    {
      return   l0(x) * dl1(y);
    }

    static double simple_quad_l0_l1xx(double x, double y)
    {
      return   d2l0(x) * l1(y);
    }

    static double simple_quad_l0_l1xy(double x, double y)
    {
      return   dl0(x) * dl1(y);
    }

    static double simple_quad_l0_l1yy(double x, double y)
    {
      return   l0(x) * d2l1(y);
    }

    static double simple_quad_l0_l2(double x, double y)
    {
      return   l0(x) * l2(y);
    }

    static double simple_quad_l0_l2x(double x, double y)
    {
      return   dl0(x) * l2(y);
    }

    static double simple_quad_l0_l2y(double x, double y)
    {
      return   l0(x) * dl2(y);
    }

    static double simple_quad_l0_l2xx(double x, double y)
    {
      return   d2l0(x) * l2(y);
    }

    static double simple_quad_l0_l2xy(double x, double y)
    {
      return   dl0(x) * dl2(y);
    }

    static double simple_quad_l0_l2yy(double x, double y)
    {
      return   l0(x) * d2l2(y);
    }

    static double simple_quad_l0_l3_0(double x, double y)
    {
      return - l0(x) * l3(y);
    }

    static double simple_quad_l0_l3x_0(double x, double y)
    {
      return - dl0(x) * l3(y);
    }

    static double simple_quad_l0_l3y_0(double x, double y)
    {
      return - l0(x) * dl3(y);
    }

    static double simple_quad_l0_l3xx_0(double x, double y)
    {
      return - d2l0(x) * l3(y);
    }

    static double simple_quad_l0_l3xy_0(double x, double y)
    {
      return - dl0(x) * dl3(y);
    }

    static double simple_quad_l0_l3yy_0(double x, double y)
    {
      return - l0(x) * d2l3(y);
    }

    static double simple_quad_l0_l3_1(double x, double y)
    {
      return -(- l0(x) * l3(y));
    }

    static double simple_quad_l0_l3x_1(double x, double y)
    {
      return -(- dl0(x) * l3(y));
    }

    static double simple_quad_l0_l3y_1(double x, double y)
    {
      return -(- l0(x) * dl3(y));
    }

    static double simple_quad_l0_l3xx_1(double x, double y)
    {
      return -(- d2l0(x) * l3(y));
    }

    static double simple_quad_l0_l3xy_1(double x, double y)
    {
      return -(- dl0(x) * dl3(y));
    }

    static double simple_quad_l0_l3yy_1(double x, double y)
    {
      return -(- l0(x) * d2l3(y));
    }

    static double simple_quad_l0_l4(double x, double y)
    {
      return   l0(x) * l4(y);
    }

    static double simple_quad_l0_l4x(double x, double y)
    {
      return   dl0(x) * l4(y);
    }

    static double simple_quad_l0_l4y(double x, double y)
    {
      return   l0(x) * dl4(y);
    }

    static double simple_quad_l0_l4xx(double x, double y)
    {
      return   d2l0(x) * l4(y);
    }

    static double simple_quad_l0_l4xy(double x, double y)
    {
      return   dl0(x) * dl4(y);
    }

    static double simple_quad_l0_l4yy(double x, double y)
    {
      return   l0(x) * d2l4(y);
    }

    static double simple_quad_l0_l5_0(double x, double y)
    {
      return - l0(x) * l5(y);
    }

    static double simple_quad_l0_l5x_0(double x, double y)
    {
      return - dl0(x) * l5(y);
    }

    static double simple_quad_l0_l5y_0(double x, double y)
    {
      return - l0(x) * dl5(y);
    }

    static double simple_quad_l0_l5xx_0(double x, double y)
    {
      return - d2l0(x) * l5(y);
    }

    static double simple_quad_l0_l5xy_0(double x, double y)
    {
      return - dl0(x) * dl5(y);
    }

    static double simple_quad_l0_l5yy_0(double x, double y)
    {
      return - l0(x) * d2l5(y);
    }

    static double simple_quad_l0_l5_1(double x, double y)
    {
      return -(- l0(x) * l5(y));
    }

    static double simple_quad_l0_l5x_1(double x, double y)
    {
      return -(- dl0(x) * l5(y));
    }

    static double simple_quad_l0_l5y_1(double x, double y)
    {
      return -(- l0(x) * dl5(y));
    }

    static double simple_quad_l0_l5xx_1(double x, double y)
    {
      return -(- d2l0(x) * l5(y));
    }

    static double simple_quad_l0_l5xy_1(double x, double y)
    {
      return -(- dl0(x) * dl5(y));
    }

    static double simple_quad_l0_l5yy_1(double x, double y)
    {
      return -(- l0(x) * d2l5(y));
    }

    static double simple_quad_l0_l6(double x, double y)
    {
      return   l0(x) * l6(y);
    }

    static double simple_quad_l0_l6x(double x, double y)
    {
      return   dl0(x) * l6(y);
    }

    static double simple_quad_l0_l6y(double x, double y)
    {
      return   l0(x) * dl6(y);
    }

    static double simple_quad_l0_l6xx(double x, double y)
    {
      return   d2l0(x) * l6(y);
    }

    static double simple_quad_l0_l6xy(double x, double y)
    {
      return   dl0(x) * dl6(y);
    }

    static double simple_quad_l0_l6yy(double x, double y)
    {
      return   l0(x) * d2l6(y);
    }

    static double simple_quad_l0_l7_0(double x, double y)
    {
      return - l0(x) * l7(y);
    }

    static double simple_quad_l0_l7x_0(double x, double y)
    {
      return - dl0(x) * l7(y);
    }

    static double simple_quad_l0_l7y_0(double x, double y)
    {
      return - l0(x) * dl7(y);
    }

    static double simple_quad_l0_l7xx_0(double x, double y)
    {
      return - d2l0(x) * l7(y);
    }

    static double simple_quad_l0_l7xy_0(double x, double y)
    {
      return - dl0(x) * dl7(y);
    }

    static double simple_quad_l0_l7yy_0(double x, double y)
    {
      return - l0(x) * d2l7(y);
    }

    static double simple_quad_l0_l7_1(double x, double y)
    {
      return -(- l0(x) * l7(y));
    }

    static double simple_quad_l0_l7x_1(double x, double y)
    {
      return -(- dl0(x) * l7(y));
    }

    static double simple_quad_l0_l7y_1(double x, double y)
    {
      return -(- l0(x) * dl7(y));
    }

    static double simple_quad_l0_l7xx_1(double x, double y)
    {
      return -(- d2l0(x) * l7(y));
    }

    static double simple_quad_l0_l7xy_1(double x, double y)
    {
      return -(- dl0(x) * dl7(y));
    }

    static double simple_quad_l0_l7yy_1(double x, double y)
    {
      return -(- l0(x) * d2l7(y));
    }

    static double simple_quad_l0_l8(double x, double y)
    {
      return   l0(x) * l8(y);
    }

    static double simple_quad_l0_l8x(double x, double y)
    {
      return   dl0(x) * l8(y);
    }

    static double simple_quad_l0_l8y(double x, double y)
    {
      return   l0(x) * dl8(y);
    }

    static double simple_quad_l0_l8xx(double x, double y)
    {
      return   d2l0(x) * l8(y);
    }

    static double simple_quad_l0_l8xy(double x, double y)
    {
      return   dl0(x) * dl8(y);
    }

    static double simple_quad_l0_l8yy(double x, double y)
    {
      return   l0(x) * d2l8(y);
    }

    static double simple_quad_l0_l9_0(double x, double y)
    {
      return - l0(x) * l9(y);
    }

    static double simple_quad_l0_l9x_0(double x, double y)
    {
      return - dl0(x) * l9(y);
    }

    static double simple_quad_l0_l9y_0(double x, double y)
    {
      return - l0(x) * dl9(y);
    }

    static double simple_quad_l0_l9xx_0(double x, double y)
    {
      return - d2l0(x) * l9(y);
    }

    static double simple_quad_l0_l9xy_0(double x, double y)
    {
      return - dl0(x) * dl9(y);
    }

    static double simple_quad_l0_l9yy_0(double x, double y)
    {
      return - l0(x) * d2l9(y);
    }

    static double simple_quad_l0_l9_1(double x, double y)
    {
      return -(- l0(x) * l9(y));
    }

    static double simple_quad_l0_l9x_1(double x, double y)
    {
      return -(- dl0(x) * l9(y));
    }

    static double simple_quad_l0_l9y_1(double x, double y)
    {
      return -(- l0(x) * dl9(y));
    }

    static double simple_quad_l0_l9xx_1(double x, double y)
    {
      return -(- d2l0(x) * l9(y));
    }

    static double simple_quad_l0_l9xy_1(double x, double y)
    {
      return -(- dl0(x) * dl9(y));
    }

    static double simple_quad_l0_l9yy_1(double x, double y)
    {
      return -(- l0(x) * d2l9(y));
    }

    static double simple_quad_l0_l10(double x, double y)
    {
      return   l0(x) * l10(y);
    }

    static double simple_quad_l0_l10x(double x, double y)
    {
      return   dl0(x) * l10(y);
    }

    static double simple_quad_l0_l10y(double x, double y)
    {
      return   l0(x) * dl10(y);
    }

    static double simple_quad_l0_l10xx(double x, double y)
    {
      return   d2l0(x) * l10(y);
    }

    static double simple_quad_l0_l10xy(double x, double y)
    {
      return   dl0(x) * dl10(y);
    }

    static double simple_quad_l0_l10yy(double x, double y)
    {
      return   l0(x) * d2l10(y);
    }

    static double simple_quad_l1_l0(double x, double y)
    {
      return   l1(x) * l0(y);
    }

    static double simple_quad_l1_l0x(double x, double y)
    {
      return   dl1(x) * l0(y);
    }

    static double simple_quad_l1_l0y(double x, double y)
    {
      return   l1(x) * dl0(y);
    }

    static double simple_quad_l1_l0xx(double x, double y)
    {
      return   d2l1(x) * l0(y);
    }

    static double simple_quad_l1_l0xy(double x, double y)
    {
      return   dl1(x) * dl0(y);
    }

    static double simple_quad_l1_l0yy(double x, double y)
    {
      return   l1(x) * d2l0(y);
    }

    static double simple_quad_l1_l1(double x, double y)
    {
      return   l1(x) * l1(y);
    }

    static double simple_quad_l1_l1x(double x, double y)
    {
      return   dl1(x) * l1(y);
    }

    static double simple_quad_l1_l1y(double x, double y)
    {
      return   l1(x) * dl1(y);
    }

    static double simple_quad_l1_l1xx(double x, double y)
    {
      return   d2l1(x) * l1(y);
    }

    static double simple_quad_l1_l1xy(double x, double y)
    {
      return   dl1(x) * dl1(y);
    }

    static double simple_quad_l1_l1yy(double x, double y)
    {
      return   l1(x) * d2l1(y);
    }

    static double simple_quad_l1_l2(double x, double y)
    {
      return   l1(x) * l2(y);
    }

    static double simple_quad_l1_l2x(double x, double y)
    {
      return   dl1(x) * l2(y);
    }

    static double simple_quad_l1_l2y(double x, double y)
    {
      return   l1(x) * dl2(y);
    }

    static double simple_quad_l1_l2xx(double x, double y)
    {
      return   d2l1(x) * l2(y);
    }

    static double simple_quad_l1_l2xy(double x, double y)
    {
      return   dl1(x) * dl2(y);
    }

    static double simple_quad_l1_l2yy(double x, double y)
    {
      return   l1(x) * d2l2(y);
    }

    static double simple_quad_l1_l3_0(double x, double y)
    {
      return   l1(x) * l3(y);
    }

    static double simple_quad_l1_l3x_0(double x, double y)
    {
      return   dl1(x) * l3(y);
    }

    static double simple_quad_l1_l3y_0(double x, double y)
    {
      return   l1(x) * dl3(y);
    }

    static double simple_quad_l1_l3xx_0(double x, double y)
    {
      return   d2l1(x) * l3(y);
    }

    static double simple_quad_l1_l3xy_0(double x, double y)
    {
      return   dl1(x) * dl3(y);
    }

    static double simple_quad_l1_l3yy_0(double x, double y)
    {
      return   l1(x) * d2l3(y);
    }

    static double simple_quad_l1_l3_1(double x, double y)
    {
      return -(  l1(x) * l3(y));
    }

    static double simple_quad_l1_l3x_1(double x, double y)
    {
      return -(  dl1(x) * l3(y));
    }

    static double simple_quad_l1_l3y_1(double x, double y)
    {
      return -(  l1(x) * dl3(y));
    }

    static double simple_quad_l1_l3xx_1(double x, double y)
    {
      return -(  d2l1(x) * l3(y));
    }

    static double simple_quad_l1_l3xy_1(double x, double y)
    {
      return -(  dl1(x) * dl3(y));
    }

    static double simple_quad_l1_l3yy_1(double x, double y)
    {
      return -(  l1(x) * d2l3(y));
    }

    static double simple_quad_l1_l4(double x, double y)
    {
      return   l1(x) * l4(y);
    }

    static double simple_quad_l1_l4x(double x, double y)
    {
      return   dl1(x) * l4(y);
    }

    static double simple_quad_l1_l4y(double x, double y)
    {
      return   l1(x) * dl4(y);
    }

    static double simple_quad_l1_l4xx(double x, double y)
    {
      return   d2l1(x) * l4(y);
    }

    static double simple_quad_l1_l4xy(double x, double y)
    {
      return   dl1(x) * dl4(y);
    }

    static double simple_quad_l1_l4yy(double x, double y)
    {
      return   l1(x) * d2l4(y);
    }

    static double simple_quad_l1_l5_0(double x, double y)
    {
      return   l1(x) * l5(y);
    }

    static double simple_quad_l1_l5x_0(double x, double y)
    {
      return   dl1(x) * l5(y);
    }

    static double simple_quad_l1_l5y_0(double x, double y)
    {
      return   l1(x) * dl5(y);
    }

    static double simple_quad_l1_l5xx_0(double x, double y)
    {
      return   d2l1(x) * l5(y);
    }

    static double simple_quad_l1_l5xy_0(double x, double y)
    {
      return   dl1(x) * dl5(y);
    }

    static double simple_quad_l1_l5yy_0(double x, double y)
    {
      return   l1(x) * d2l5(y);
    }

    static double simple_quad_l1_l5_1(double x, double y)
    {
      return -(  l1(x) * l5(y));
    }

    static double simple_quad_l1_l5x_1(double x, double y)
    {
      return -(  dl1(x) * l5(y));
    }

    static double simple_quad_l1_l5y_1(double x, double y)
    {
      return -(  l1(x) * dl5(y));
    }

    static double simple_quad_l1_l5xx_1(double x, double y)
    {
      return -(  d2l1(x) * l5(y));
    }

    static double simple_quad_l1_l5xy_1(double x, double y)
    {
      return -(  dl1(x) * dl5(y));
    }

    static double simple_quad_l1_l5yy_1(double x, double y)
    {
      return -(  l1(x) * d2l5(y));
    }

    static double simple_quad_l1_l6(double x, double y)
    {
      return   l1(x) * l6(y);
    }

    static double simple_quad_l1_l6x(double x, double y)
    {
      return   dl1(x) * l6(y);
    }

    static double simple_quad_l1_l6y(double x, double y)
    {
      return   l1(x) * dl6(y);
    }

    static double simple_quad_l1_l6xx(double x, double y)
    {
      return   d2l1(x) * l6(y);
    }

    static double simple_quad_l1_l6xy(double x, double y)
    {
      return   dl1(x) * dl6(y);
    }

    static double simple_quad_l1_l6yy(double x, double y)
    {
      return   l1(x) * d2l6(y);
    }

    static double simple_quad_l1_l7_0(double x, double y)
    {
      return   l1(x) * l7(y);
    }

    static double simple_quad_l1_l7x_0(double x, double y)
    {
      return   dl1(x) * l7(y);
    }

    static double simple_quad_l1_l7y_0(double x, double y)
    {
      return   l1(x) * dl7(y);
    }

    static double simple_quad_l1_l7xx_0(double x, double y)
    {
      return   d2l1(x) * l7(y);
    }

    static double simple_quad_l1_l7xy_0(double x, double y)
    {
      return   dl1(x) * dl7(y);
    }

    static double simple_quad_l1_l7yy_0(double x, double y)
    {
      return   l1(x) * d2l7(y);
    }

    static double simple_quad_l1_l7_1(double x, double y)
    {
      return -(  l1(x) * l7(y));
    }

    static double simple_quad_l1_l7x_1(double x, double y)
    {
      return -(  dl1(x) * l7(y));
    }

    static double simple_quad_l1_l7y_1(double x, double y)
    {
      return -(  l1(x) * dl7(y));
    }

    static double simple_quad_l1_l7xx_1(double x, double y)
    {
      return -(  d2l1(x) * l7(y));
    }

    static double simple_quad_l1_l7xy_1(double x, double y)
    {
      return -(  dl1(x) * dl7(y));
    }

    static double simple_quad_l1_l7yy_1(double x, double y)
    {
      return -(  l1(x) * d2l7(y));
    }

    static double simple_quad_l1_l8(double x, double y)
    {
      return   l1(x) * l8(y);
    }

    static double simple_quad_l1_l8x(double x, double y)
    {
      return   dl1(x) * l8(y);
    }

    static double simple_quad_l1_l8y(double x, double y)
    {
      return   l1(x) * dl8(y);
    }

    static double simple_quad_l1_l8xx(double x, double y)
    {
      return   d2l1(x) * l8(y);
    }

    static double simple_quad_l1_l8xy(double x, double y)
    {
      return   dl1(x) * dl8(y);
    }

    static double simple_quad_l1_l8yy(double x, double y)
    {
      return   l1(x) * d2l8(y);
    }

    static double simple_quad_l1_l9_0(double x, double y)
    {
      return   l1(x) * l9(y);
    }

    static double simple_quad_l1_l9x_0(double x, double y)
    {
      return   dl1(x) * l9(y);
    }

    static double simple_quad_l1_l9y_0(double x, double y)
    {
      return   l1(x) * dl9(y);
    }

    static double simple_quad_l1_l9xx_0(double x, double y)
    {
      return   d2l1(x) * l9(y);
    }

    static double simple_quad_l1_l9xy_0(double x, double y)
    {
      return   dl1(x) * dl9(y);
    }

    static double simple_quad_l1_l9yy_0(double x, double y)
    {
      return   l1(x) * d2l9(y);
    }

    static double simple_quad_l1_l9_1(double x, double y)
    {
      return -(  l1(x) * l9(y));
    }

    static double simple_quad_l1_l9x_1(double x, double y)
    {
      return -(  dl1(x) * l9(y));
    }

    static double simple_quad_l1_l9y_1(double x, double y)
    {
      return -(  l1(x) * dl9(y));
    }

    static double simple_quad_l1_l9xx_1(double x, double y)
    {
      return -(  d2l1(x) * l9(y));
    }

    static double simple_quad_l1_l9xy_1(double x, double y)
    {
      return -(  dl1(x) * dl9(y));
    }

    static double simple_quad_l1_l9yy_1(double x, double y)
    {
      return -(  l1(x) * d2l9(y));
    }

    static double simple_quad_l1_l10(double x, double y)
    {
      return   l1(x) * l10(y);
    }

    static double simple_quad_l1_l10x(double x, double y)
    {
      return   dl1(x) * l10(y);
    }

    static double simple_quad_l1_l10y(double x, double y)
    {
      return   l1(x) * dl10(y);
    }

    static double simple_quad_l1_l10xx(double x, double y)
    {
      return   d2l1(x) * l10(y);
    }

    static double simple_quad_l1_l10xy(double x, double y)
    {
      return   dl1(x) * dl10(y);
    }

    static double simple_quad_l1_l10yy(double x, double y)
    {
      return   l1(x) * d2l10(y);
    }

    static double simple_quad_l2_l0(double x, double y)
    {
      return   l2(x) * l0(y);
    }

    static double simple_quad_l2_l0x(double x, double y)
    {
      return   dl2(x) * l0(y);
    }

    static double simple_quad_l2_l0y(double x, double y)
    {
      return   l2(x) * dl0(y);
    }

    static double simple_quad_l2_l0xx(double x, double y)
    {
      return   d2l2(x) * l0(y);
    }

    static double simple_quad_l2_l0xy(double x, double y)
    {
      return   dl2(x) * dl0(y);
    }

    static double simple_quad_l2_l0yy(double x, double y)
    {
      return   l2(x) * d2l0(y);
    }

    static double simple_quad_l2_l1(double x, double y)
    {
      return   l2(x) * l1(y);
    }

    static double simple_quad_l2_l1x(double x, double y)
    {
      return   dl2(x) * l1(y);
    }

    static double simple_quad_l2_l1y(double x, double y)
    {
      return   l2(x) * dl1(y);
    }

    static double simple_quad_l2_l1xx(double x, double y)
    {
      return   d2l2(x) * l1(y);
    }

    static double simple_quad_l2_l1xy(double x, double y)
    {
      return   dl2(x) * dl1(y);
    }

    static double simple_quad_l2_l1yy(double x, double y)
    {
      return   l2(x) * d2l1(y);
    }

    static double simple_quad_l2_l2(double x, double y)
    {
      return   l2(x) * l2(y);
    }

    static double simple_quad_l2_l2x(double x, double y)
    {
      return   dl2(x) * l2(y);
    }

    static double simple_quad_l2_l2y(double x, double y)
    {
      return   l2(x) * dl2(y);
    }

    static double simple_quad_l2_l2xx(double x, double y)
    {
      return   d2l2(x) * l2(y);
    }

    static double simple_quad_l2_l2xy(double x, double y)
    {
      return   dl2(x) * dl2(y);
    }

    static double simple_quad_l2_l2yy(double x, double y)
    {
      return   l2(x) * d2l2(y);
    }

    static double simple_quad_l2_l3(double x, double y)
    {
      return   l2(x) * l3(y);
    }

    static double simple_quad_l2_l3x(double x, double y)
    {
      return   dl2(x) * l3(y);
    }

    static double simple_quad_l2_l3y(double x, double y)
    {
      return   l2(x) * dl3(y);
    }

    static double simple_quad_l2_l3xx(double x, double y)
    {
      return   d2l2(x) * l3(y);
    }

    static double simple_quad_l2_l3xy(double x, double y)
    {
      return   dl2(x) * dl3(y);
    }

    static double simple_quad_l2_l3yy(double x, double y)
    {
      return   l2(x) * d2l3(y);
    }

    static double simple_quad_l2_l4(double x, double y)
    {
      return   l2(x) * l4(y);
    }

    static double simple_quad_l2_l4x(double x, double y)
    {
      return   dl2(x) * l4(y);
    }

    static double simple_quad_l2_l4y(double x, double y)
    {
      return   l2(x) * dl4(y);
    }

    static double simple_quad_l2_l4xx(double x, double y)
    {
      return   d2l2(x) * l4(y);
    }

    static double simple_quad_l2_l4xy(double x, double y)
    {
      return   dl2(x) * dl4(y);
    }

    static double simple_quad_l2_l4yy(double x, double y)
    {
      return   l2(x) * d2l4(y);
    }

    static double simple_quad_l2_l5(double x, double y)
    {
      return   l2(x) * l5(y);
    }

    static double simple_quad_l2_l5x(double x, double y)
    {
      return   dl2(x) * l5(y);
    }

    static double simple_quad_l2_l5y(double x, double y)
    {
      return   l2(x) * dl5(y);
    }

    static double simple_quad_l2_l5xx(double x, double y)
    {
      return   d2l2(x) * l5(y);
    }

    static double simple_quad_l2_l5xy(double x, double y)
    {
      return   dl2(x) * dl5(y);
    }

    static double simple_quad_l2_l5yy(double x, double y)
    {
      return   l2(x) * d2l5(y);
    }

    static double simple_quad_l2_l6(double x, double y)
    {
      return   l2(x) * l6(y);
    }

    static double simple_quad_l2_l6x(double x, double y)
    {
      return   dl2(x) * l6(y);
    }

    static double simple_quad_l2_l6y(double x, double y)
    {
      return   l2(x) * dl6(y);
    }

    static double simple_quad_l2_l6xx(double x, double y)
    {
      return   d2l2(x) * l6(y);
    }

    static double simple_quad_l2_l6xy(double x, double y)
    {
      return   dl2(x) * dl6(y);
    }

    static double simple_quad_l2_l6yy(double x, double y)
    {
      return   l2(x) * d2l6(y);
    }

    static double simple_quad_l2_l7(double x, double y)
    {
      return   l2(x) * l7(y);
    }

    static double simple_quad_l2_l7x(double x, double y)
    {
      return   dl2(x) * l7(y);
    }

    static double simple_quad_l2_l7y(double x, double y)
    {
      return   l2(x) * dl7(y);
    }

    static double simple_quad_l2_l7xx(double x, double y)
    {
      return   d2l2(x) * l7(y);
    }

    static double simple_quad_l2_l7xy(double x, double y)
    {
      return   dl2(x) * dl7(y);
    }

    static double simple_quad_l2_l7yy(double x, double y)
    {
      return   l2(x) * d2l7(y);
    }

    static double simple_quad_l2_l8(double x, double y)
    {
      return   l2(x) * l8(y);
    }

    static double simple_quad_l2_l8x(double x, double y)
    {
      return   dl2(x) * l8(y);
    }

    static double simple_quad_l2_l8y(double x, double y)
    {
      return   l2(x) * dl8(y);
    }

    static double simple_quad_l2_l8xx(double x, double y)
    {
      return   d2l2(x) * l8(y);
    }

    static double simple_quad_l2_l8xy(double x, double y)
    {
      return   dl2(x) * dl8(y);
    }

    static double simple_quad_l2_l8yy(double x, double y)
    {
      return   l2(x) * d2l8(y);
    }

    static double simple_quad_l2_l9(double x, double y)
    {
      return   l2(x) * l9(y);
    }

    static double simple_quad_l2_l9x(double x, double y)
    {
      return   dl2(x) * l9(y);
    }

    static double simple_quad_l2_l9y(double x, double y)
    {
      return   l2(x) * dl9(y);
    }

    static double simple_quad_l2_l9xx(double x, double y)
    {
      return   d2l2(x) * l9(y);
    }

    static double simple_quad_l2_l9xy(double x, double y)
    {
      return   dl2(x) * dl9(y);
    }

    static double simple_quad_l2_l9yy(double x, double y)
    {
      return   l2(x) * d2l9(y);
    }

    static double simple_quad_l2_l10(double x, double y)
    {
      return   l2(x) * l10(y);
    }

    static double simple_quad_l2_l10x(double x, double y)
    {
      return   dl2(x) * l10(y);
    }

    static double simple_quad_l2_l10y(double x, double y)
    {
      return   l2(x) * dl10(y);
    }

    static double simple_quad_l2_l10xx(double x, double y)
    {
      return   d2l2(x) * l10(y);
    }

    static double simple_quad_l2_l10xy(double x, double y)
    {
      return   dl2(x) * dl10(y);
    }

    static double simple_quad_l2_l10yy(double x, double y)
    {
      return   l2(x) * d2l10(y);
    }

    static double simple_quad_l3_l0_0(double x, double y)
    {
      return   l3(x) * l0(y);
    }

    static double simple_quad_l3_l0x_0(double x, double y)
    {
      return   dl3(x) * l0(y);
    }

    static double simple_quad_l3_l0y_0(double x, double y)
    {
      return   l3(x) * dl0(y);
    }

    static double simple_quad_l3_l0xx_0(double x, double y)
    {
      return   d2l3(x) * l0(y);
    }

    static double simple_quad_l3_l0xy_0(double x, double y)
    {
      return   dl3(x) * dl0(y);
    }

    static double simple_quad_l3_l0yy_0(double x, double y)
    {
      return   l3(x) * d2l0(y);
    }

    static double simple_quad_l3_l0_1(double x, double y)
    {
      return -(  l3(x) * l0(y));
    }

    static double simple_quad_l3_l0x_1(double x, double y)
    {
      return -(  dl3(x) * l0(y));
    }

    static double simple_quad_l3_l0y_1(double x, double y)
    {
      return -(  l3(x) * dl0(y));
    }

    static double simple_quad_l3_l0xx_1(double x, double y)
    {
      return -(  d2l3(x) * l0(y));
    }

    static double simple_quad_l3_l0xy_1(double x, double y)
    {
      return -(  dl3(x) * dl0(y));
    }

    static double simple_quad_l3_l0yy_1(double x, double y)
    {
      return -(  l3(x) * d2l0(y));
    }

    static double simple_quad_l3_l1_0(double x, double y)
    {
      return - l3(x) * l1(y);
    }

    static double simple_quad_l3_l1x_0(double x, double y)
    {
      return - dl3(x) * l1(y);
    }

    static double simple_quad_l3_l1y_0(double x, double y)
    {
      return - l3(x) * dl1(y);
    }

    static double simple_quad_l3_l1xx_0(double x, double y)
    {
      return - d2l3(x) * l1(y);
    }

    static double simple_quad_l3_l1xy_0(double x, double y)
    {
      return - dl3(x) * dl1(y);
    }

    static double simple_quad_l3_l1yy_0(double x, double y)
    {
      return - l3(x) * d2l1(y);
    }

    static double simple_quad_l3_l1_1(double x, double y)
    {
      return -(- l3(x) * l1(y));
    }

    static double simple_quad_l3_l1x_1(double x, double y)
    {
      return -(- dl3(x) * l1(y));
    }

    static double simple_quad_l3_l1y_1(double x, double y)
    {
      return -(- l3(x) * dl1(y));
    }

    static double simple_quad_l3_l1xx_1(double x, double y)
    {
      return -(- d2l3(x) * l1(y));
    }

    static double simple_quad_l3_l1xy_1(double x, double y)
    {
      return -(- dl3(x) * dl1(y));
    }

    static double simple_quad_l3_l1yy_1(double x, double y)
    {
      return -(- l3(x) * d2l1(y));
    }

    static double simple_quad_l3_l2(double x, double y)
    {
      return   l3(x) * l2(y);
    }

    static double simple_quad_l3_l2x(double x, double y)
    {
      return   dl3(x) * l2(y);
    }

    static double simple_quad_l3_l2y(double x, double y)
    {
      return   l3(x) * dl2(y);
    }

    static double simple_quad_l3_l2xx(double x, double y)
    {
      return   d2l3(x) * l2(y);
    }

    static double simple_quad_l3_l2xy(double x, double y)
    {
      return   dl3(x) * dl2(y);
    }

    static double simple_quad_l3_l2yy(double x, double y)
    {
      return   l3(x) * d2l2(y);
    }

    static double simple_quad_l3_l3(double x, double y)
    {
      return   l3(x) * l3(y);
    }

    static double simple_quad_l3_l3x(double x, double y)
    {
      return   dl3(x) * l3(y);
    }

    static double simple_quad_l3_l3y(double x, double y)
    {
      return   l3(x) * dl3(y);
    }

    static double simple_quad_l3_l3xx(double x, double y)
    {
      return   d2l3(x) * l3(y);
    }

    static double simple_quad_l3_l3xy(double x, double y)
    {
      return   dl3(x) * dl3(y);
    }

    static double simple_quad_l3_l3yy(double x, double y)
    {
      return   l3(x) * d2l3(y);
    }

    static double simple_quad_l3_l4(double x, double y)
    {
      return   l3(x) * l4(y);
    }

    static double simple_quad_l3_l4x(double x, double y)
    {
      return   dl3(x) * l4(y);
    }

    static double simple_quad_l3_l4y(double x, double y)
    {
      return   l3(x) * dl4(y);
    }

    static double simple_quad_l3_l4xx(double x, double y)
    {
      return   d2l3(x) * l4(y);
    }

    static double simple_quad_l3_l4xy(double x, double y)
    {
      return   dl3(x) * dl4(y);
    }

    static double simple_quad_l3_l4yy(double x, double y)
    {
      return   l3(x) * d2l4(y);
    }

    static double simple_quad_l3_l5(double x, double y)
    {
      return   l3(x) * l5(y);
    }

    static double simple_quad_l3_l5x(double x, double y)
    {
      return   dl3(x) * l5(y);
    }

    static double simple_quad_l3_l5y(double x, double y)
    {
      return   l3(x) * dl5(y);
    }

    static double simple_quad_l3_l5xx(double x, double y)
    {
      return   d2l3(x) * l5(y);
    }

    static double simple_quad_l3_l5xy(double x, double y)
    {
      return   dl3(x) * dl5(y);
    }

    static double simple_quad_l3_l5yy(double x, double y)
    {
      return   l3(x) * d2l5(y);
    }

    static double simple_quad_l3_l6(double x, double y)
    {
      return   l3(x) * l6(y);
    }

    static double simple_quad_l3_l6x(double x, double y)
    {
      return   dl3(x) * l6(y);
    }

    static double simple_quad_l3_l6y(double x, double y)
    {
      return   l3(x) * dl6(y);
    }

    static double simple_quad_l3_l6xx(double x, double y)
    {
      return   d2l3(x) * l6(y);
    }

    static double simple_quad_l3_l6xy(double x, double y)
    {
      return   dl3(x) * dl6(y);
    }

    static double simple_quad_l3_l6yy(double x, double y)
    {
      return   l3(x) * d2l6(y);
    }

    static double simple_quad_l3_l7(double x, double y)
    {
      return   l3(x) * l7(y);
    }

    static double simple_quad_l3_l7x(double x, double y)
    {
      return   dl3(x) * l7(y);
    }

    static double simple_quad_l3_l7y(double x, double y)
    {
      return   l3(x) * dl7(y);
    }

    static double simple_quad_l3_l7xx(double x, double y)
    {
      return   d2l3(x) * l7(y);
    }

    static double simple_quad_l3_l7xy(double x, double y)
    {
      return   dl3(x) * dl7(y);
    }

    static double simple_quad_l3_l7yy(double x, double y)
    {
      return   l3(x) * d2l7(y);
    }

    static double simple_quad_l3_l8(double x, double y)
    {
      return   l3(x) * l8(y);
    }

    static double simple_quad_l3_l8x(double x, double y)
    {
      return   dl3(x) * l8(y);
    }

    static double simple_quad_l3_l8y(double x, double y)
    {
      return   l3(x) * dl8(y);
    }

    static double simple_quad_l3_l8xx(double x, double y)
    {
      return   d2l3(x) * l8(y);
    }

    static double simple_quad_l3_l8xy(double x, double y)
    {
      return   dl3(x) * dl8(y);
    }

    static double simple_quad_l3_l8yy(double x, double y)
    {
      return   l3(x) * d2l8(y);
    }

    static double simple_quad_l3_l9(double x, double y)
    {
      return   l3(x) * l9(y);
    }

    static double simple_quad_l3_l9x(double x, double y)
    {
      return   dl3(x) * l9(y);
    }

    static double simple_quad_l3_l9y(double x, double y)
    {
      return   l3(x) * dl9(y);
    }

    static double simple_quad_l3_l9xx(double x, double y)
    {
      return   d2l3(x) * l9(y);
    }

    static double simple_quad_l3_l9xy(double x, double y)
    {
      return   dl3(x) * dl9(y);
    }

    static double simple_quad_l3_l9yy(double x, double y)
    {
      return   l3(x) * d2l9(y);
    }

    static double simple_quad_l3_l10(double x, double y)
    {
      return   l3(x) * l10(y);
    }

    static double simple_quad_l3_l10x(double x, double y)
    {
      return   dl3(x) * l10(y);
    }

    static double simple_quad_l3_l10y(double x, double y)
    {
      return   l3(x) * dl10(y);
    }

    static double simple_quad_l3_l10xx(double x, double y)
    {
      return   d2l3(x) * l10(y);
    }

    static double simple_quad_l3_l10xy(double x, double y)
    {
      return   dl3(x) * dl10(y);
    }

    static double simple_quad_l3_l10yy(double x, double y)
    {
      return   l3(x) * d2l10(y);
    }

    static double simple_quad_l4_l0(double x, double y)
    {
      return   l4(x) * l0(y);
    }

    static double simple_quad_l4_l0x(double x, double y)
    {
      return   dl4(x) * l0(y);
    }

    static double simple_quad_l4_l0y(double x, double y)
    {
      return   l4(x) * dl0(y);
    }

    static double simple_quad_l4_l0xx(double x, double y)
    {
      return   d2l4(x) * l0(y);
    }

    static double simple_quad_l4_l0xy(double x, double y)
    {
      return   dl4(x) * dl0(y);
    }

    static double simple_quad_l4_l0yy(double x, double y)
    {
      return   l4(x) * d2l0(y);
    }

    static double simple_quad_l4_l1(double x, double y)
    {
      return   l4(x) * l1(y);
    }

    static double simple_quad_l4_l1x(double x, double y)
    {
      return   dl4(x) * l1(y);
    }

    static double simple_quad_l4_l1y(double x, double y)
    {
      return   l4(x) * dl1(y);
    }

    static double simple_quad_l4_l1xx(double x, double y)
    {
      return   d2l4(x) * l1(y);
    }

    static double simple_quad_l4_l1xy(double x, double y)
    {
      return   dl4(x) * dl1(y);
    }

    static double simple_quad_l4_l1yy(double x, double y)
    {
      return   l4(x) * d2l1(y);
    }

    static double simple_quad_l4_l2(double x, double y)
    {
      return   l4(x) * l2(y);
    }

    static double simple_quad_l4_l2x(double x, double y)
    {
      return   dl4(x) * l2(y);
    }

    static double simple_quad_l4_l2y(double x, double y)
    {
      return   l4(x) * dl2(y);
    }

    static double simple_quad_l4_l2xx(double x, double y)
    {
      return   d2l4(x) * l2(y);
    }

    static double simple_quad_l4_l2xy(double x, double y)
    {
      return   dl4(x) * dl2(y);
    }

    static double simple_quad_l4_l2yy(double x, double y)
    {
      return   l4(x) * d2l2(y);
    }

    static double simple_quad_l4_l3(double x, double y)
    {
      return   l4(x) * l3(y);
    }

    static double simple_quad_l4_l3x(double x, double y)
    {
      return   dl4(x) * l3(y);
    }

    static double simple_quad_l4_l3y(double x, double y)
    {
      return   l4(x) * dl3(y);
    }

    static double simple_quad_l4_l3xx(double x, double y)
    {
      return   d2l4(x) * l3(y);
    }

    static double simple_quad_l4_l3xy(double x, double y)
    {
      return   dl4(x) * dl3(y);
    }

    static double simple_quad_l4_l3yy(double x, double y)
    {
      return   l4(x) * d2l3(y);
    }

    static double simple_quad_l4_l4(double x, double y)
    {
      return   l4(x) * l4(y);
    }

    static double simple_quad_l4_l4x(double x, double y)
    {
      return   dl4(x) * l4(y);
    }

    static double simple_quad_l4_l4y(double x, double y)
    {
      return   l4(x) * dl4(y);
    }

    static double simple_quad_l4_l4xx(double x, double y)
    {
      return   d2l4(x) * l4(y);
    }

    static double simple_quad_l4_l4xy(double x, double y)
    {
      return   dl4(x) * dl4(y);
    }

    static double simple_quad_l4_l4yy(double x, double y)
    {
      return   l4(x) * d2l4(y);
    }

    static double simple_quad_l4_l5(double x, double y)
    {
      return   l4(x) * l5(y);
    }

    static double simple_quad_l4_l5x(double x, double y)
    {
      return   dl4(x) * l5(y);
    }

    static double simple_quad_l4_l5y(double x, double y)
    {
      return   l4(x) * dl5(y);
    }

    static double simple_quad_l4_l5xx(double x, double y)
    {
      return   d2l4(x) * l5(y);
    }

    static double simple_quad_l4_l5xy(double x, double y)
    {
      return   dl4(x) * dl5(y);
    }

    static double simple_quad_l4_l5yy(double x, double y)
    {
      return   l4(x) * d2l5(y);
    }

    static double simple_quad_l4_l6(double x, double y)
    {
      return   l4(x) * l6(y);
    }

    static double simple_quad_l4_l6x(double x, double y)
    {
      return   dl4(x) * l6(y);
    }

    static double simple_quad_l4_l6y(double x, double y)
    {
      return   l4(x) * dl6(y);
    }

    static double simple_quad_l4_l6xx(double x, double y)
    {
      return   d2l4(x) * l6(y);
    }

    static double simple_quad_l4_l6xy(double x, double y)
    {
      return   dl4(x) * dl6(y);
    }

    static double simple_quad_l4_l6yy(double x, double y)
    {
      return   l4(x) * d2l6(y);
    }

    static double simple_quad_l4_l7(double x, double y)
    {
      return   l4(x) * l7(y);
    }

    static double simple_quad_l4_l7x(double x, double y)
    {
      return   dl4(x) * l7(y);
    }

    static double simple_quad_l4_l7y(double x, double y)
    {
      return   l4(x) * dl7(y);
    }

    static double simple_quad_l4_l7xx(double x, double y)
    {
      return   d2l4(x) * l7(y);
    }

    static double simple_quad_l4_l7xy(double x, double y)
    {
      return   dl4(x) * dl7(y);
    }

    static double simple_quad_l4_l7yy(double x, double y)
    {
      return   l4(x) * d2l7(y);
    }

    static double simple_quad_l4_l8(double x, double y)
    {
      return   l4(x) * l8(y);
    }

    static double simple_quad_l4_l8x(double x, double y)
    {
      return   dl4(x) * l8(y);
    }

    static double simple_quad_l4_l8y(double x, double y)
    {
      return   l4(x) * dl8(y);
    }

    static double simple_quad_l4_l8xx(double x, double y)
    {
      return   d2l4(x) * l8(y);
    }

    static double simple_quad_l4_l8xy(double x, double y)
    {
      return   dl4(x) * dl8(y);
    }

    static double simple_quad_l4_l8yy(double x, double y)
    {
      return   l4(x) * d2l8(y);
    }

    static double simple_quad_l4_l9(double x, double y)
    {
      return   l4(x) * l9(y);
    }

    static double simple_quad_l4_l9x(double x, double y)
    {
      return   dl4(x) * l9(y);
    }

    static double simple_quad_l4_l9y(double x, double y)
    {
      return   l4(x) * dl9(y);
    }

    static double simple_quad_l4_l9xx(double x, double y)
    {
      return   d2l4(x) * l9(y);
    }

    static double simple_quad_l4_l9xy(double x, double y)
    {
      return   dl4(x) * dl9(y);
    }

    static double simple_quad_l4_l9yy(double x, double y)
    {
      return   l4(x) * d2l9(y);
    }

    static double simple_quad_l4_l10(double x, double y)
    {
      return   l4(x) * l10(y);
    }

    static double simple_quad_l4_l10x(double x, double y)
    {
      return   dl4(x) * l10(y);
    }

    static double simple_quad_l4_l10y(double x, double y)
    {
      return   l4(x) * dl10(y);
    }

    static double simple_quad_l4_l10xx(double x, double y)
    {
      return   d2l4(x) * l10(y);
    }

    static double simple_quad_l4_l10xy(double x, double y)
    {
      return   dl4(x) * dl10(y);
    }

    static double simple_quad_l4_l10yy(double x, double y)
    {
      return   l4(x) * d2l10(y);
    }

    static double simple_quad_l5_l0_0(double x, double y)
    {
      return   l5(x) * l0(y);
    }

    static double simple_quad_l5_l0x_0(double x, double y)
    {
      return   dl5(x) * l0(y);
    }

    static double simple_quad_l5_l0y_0(double x, double y)
    {
      return   l5(x) * dl0(y);
    }

    static double simple_quad_l5_l0xx_0(double x, double y)
    {
      return   d2l5(x) * l0(y);
    }

    static double simple_quad_l5_l0xy_0(double x, double y)
    {
      return   dl5(x) * dl0(y);
    }

    static double simple_quad_l5_l0yy_0(double x, double y)
    {
      return   l5(x) * d2l0(y);
    }

    static double simple_quad_l5_l0_1(double x, double y)
    {
      return -(  l5(x) * l0(y));
    }

    static double simple_quad_l5_l0x_1(double x, double y)
    {
      return -(  dl5(x) * l0(y));
    }

    static double simple_quad_l5_l0y_1(double x, double y)
    {
      return -(  l5(x) * dl0(y));
    }

    static double simple_quad_l5_l0xx_1(double x, double y)
    {
      return -(  d2l5(x) * l0(y));
    }

    static double simple_quad_l5_l0xy_1(double x, double y)
    {
      return -(  dl5(x) * dl0(y));
    }

    static double simple_quad_l5_l0yy_1(double x, double y)
    {
      return -(  l5(x) * d2l0(y));
    }

    static double simple_quad_l5_l1_0(double x, double y)
    {
      return - l5(x) * l1(y);
    }

    static double simple_quad_l5_l1x_0(double x, double y)
    {
      return - dl5(x) * l1(y);
    }

    static double simple_quad_l5_l1y_0(double x, double y)
    {
      return - l5(x) * dl1(y);
    }

    static double simple_quad_l5_l1xx_0(double x, double y)
    {
      return - d2l5(x) * l1(y);
    }

    static double simple_quad_l5_l1xy_0(double x, double y)
    {
      return - dl5(x) * dl1(y);
    }

    static double simple_quad_l5_l1yy_0(double x, double y)
    {
      return - l5(x) * d2l1(y);
    }

    static double simple_quad_l5_l1_1(double x, double y)
    {
      return -(- l5(x) * l1(y));
    }

    static double simple_quad_l5_l1x_1(double x, double y)
    {
      return -(- dl5(x) * l1(y));
    }

    static double simple_quad_l5_l1y_1(double x, double y)
    {
      return -(- l5(x) * dl1(y));
    }

    static double simple_quad_l5_l1xx_1(double x, double y)
    {
      return -(- d2l5(x) * l1(y));
    }

    static double simple_quad_l5_l1xy_1(double x, double y)
    {
      return -(- dl5(x) * dl1(y));
    }

    static double simple_quad_l5_l1yy_1(double x, double y)
    {
      return -(- l5(x) * d2l1(y));
    }

    static double simple_quad_l5_l2(double x, double y)
    {
      return   l5(x) * l2(y);
    }

    static double simple_quad_l5_l2x(double x, double y)
    {
      return   dl5(x) * l2(y);
    }

    static double simple_quad_l5_l2y(double x, double y)
    {
      return   l5(x) * dl2(y);
    }

    static double simple_quad_l5_l2xx(double x, double y)
    {
      return   d2l5(x) * l2(y);
    }

    static double simple_quad_l5_l2xy(double x, double y)
    {
      return   dl5(x) * dl2(y);
    }

    static double simple_quad_l5_l2yy(double x, double y)
    {
      return   l5(x) * d2l2(y);
    }

    static double simple_quad_l5_l3(double x, double y)
    {
      return   l5(x) * l3(y);
    }

    static double simple_quad_l5_l3x(double x, double y)
    {
      return   dl5(x) * l3(y);
    }

    static double simple_quad_l5_l3y(double x, double y)
    {
      return   l5(x) * dl3(y);
    }

    static double simple_quad_l5_l3xx(double x, double y)
    {
      return   d2l5(x) * l3(y);
    }

    static double simple_quad_l5_l3xy(double x, double y)
    {
      return   dl5(x) * dl3(y);
    }

    static double simple_quad_l5_l3yy(double x, double y)
    {
      return   l5(x) * d2l3(y);
    }

    static double simple_quad_l5_l4(double x, double y)
    {
      return   l5(x) * l4(y);
    }

    static double simple_quad_l5_l4x(double x, double y)
    {
      return   dl5(x) * l4(y);
    }

    static double simple_quad_l5_l4y(double x, double y)
    {
      return   l5(x) * dl4(y);
    }

    static double simple_quad_l5_l4xx(double x, double y)
    {
      return   d2l5(x) * l4(y);
    }

    static double simple_quad_l5_l4xy(double x, double y)
    {
      return   dl5(x) * dl4(y);
    }

    static double simple_quad_l5_l4yy(double x, double y)
    {
      return   l5(x) * d2l4(y);
    }

    static double simple_quad_l5_l5(double x, double y)
    {
      return   l5(x) * l5(y);
    }

    static double simple_quad_l5_l5x(double x, double y)
    {
      return   dl5(x) * l5(y);
    }

    static double simple_quad_l5_l5y(double x, double y)
    {
      return   l5(x) * dl5(y);
    }

    static double simple_quad_l5_l5xx(double x, double y)
    {
      return   d2l5(x) * l5(y);
    }

    static double simple_quad_l5_l5xy(double x, double y)
    {
      return   dl5(x) * dl5(y);
    }

    static double simple_quad_l5_l5yy(double x, double y)
    {
      return   l5(x) * d2l5(y);
    }

    static double simple_quad_l5_l6(double x, double y)
    {
      return   l5(x) * l6(y);
    }

    static double simple_quad_l5_l6x(double x, double y)
    {
      return   dl5(x) * l6(y);
    }

    static double simple_quad_l5_l6y(double x, double y)
    {
      return   l5(x) * dl6(y);
    }

    static double simple_quad_l5_l6xx(double x, double y)
    {
      return   d2l5(x) * l6(y);
    }

    static double simple_quad_l5_l6xy(double x, double y)
    {
      return   dl5(x) * dl6(y);
    }

    static double simple_quad_l5_l6yy(double x, double y)
    {
      return   l5(x) * d2l6(y);
    }

    static double simple_quad_l5_l7(double x, double y)
    {
      return   l5(x) * l7(y);
    }

    static double simple_quad_l5_l7x(double x, double y)
    {
      return   dl5(x) * l7(y);
    }

    static double simple_quad_l5_l7y(double x, double y)
    {
      return   l5(x) * dl7(y);
    }

    static double simple_quad_l5_l7xx(double x, double y)
    {
      return   d2l5(x) * l7(y);
    }

    static double simple_quad_l5_l7xy(double x, double y)
    {
      return   dl5(x) * dl7(y);
    }

    static double simple_quad_l5_l7yy(double x, double y)
    {
      return   l5(x) * d2l7(y);
    }

    static double simple_quad_l5_l8(double x, double y)
    {
      return   l5(x) * l8(y);
    }

    static double simple_quad_l5_l8x(double x, double y)
    {
      return   dl5(x) * l8(y);
    }

    static double simple_quad_l5_l8y(double x, double y)
    {
      return   l5(x) * dl8(y);
    }

    static double simple_quad_l5_l8xx(double x, double y)
    {
      return   d2l5(x) * l8(y);
    }

    static double simple_quad_l5_l8xy(double x, double y)
    {
      return   dl5(x) * dl8(y);
    }

    static double simple_quad_l5_l8yy(double x, double y)
    {
      return   l5(x) * d2l8(y);
    }

    static double simple_quad_l5_l9(double x, double y)
    {
      return   l5(x) * l9(y);
    }

    static double simple_quad_l5_l9x(double x, double y)
    {
      return   dl5(x) * l9(y);
    }

    static double simple_quad_l5_l9y(double x, double y)
    {
      return   l5(x) * dl9(y);
    }

    static double simple_quad_l5_l9xx(double x, double y)
    {
      return   d2l5(x) * l9(y);
    }

    static double simple_quad_l5_l9xy(double x, double y)
    {
      return   dl5(x) * dl9(y);
    }

    static double simple_quad_l5_l9yy(double x, double y)
    {
      return   l5(x) * d2l9(y);
    }

    static double simple_quad_l5_l10(double x, double y)
    {
      return   l5(x) * l10(y);
    }

    static double simple_quad_l5_l10x(double x, double y)
    {
      return   dl5(x) * l10(y);
    }

    static double simple_quad_l5_l10y(double x, double y)
    {
      return   l5(x) * dl10(y);
    }

    static double simple_quad_l5_l10xx(double x, double y)
    {
      return   d2l5(x) * l10(y);
    }

    static double simple_quad_l5_l10xy(double x, double y)
    {
      return   dl5(x) * dl10(y);
    }

    static double simple_quad_l5_l10yy(double x, double y)
    {
      return   l5(x) * d2l10(y);
    }

    static double simple_quad_l6_l0(double x, double y)
    {
      return   l6(x) * l0(y);
    }

    static double simple_quad_l6_l0x(double x, double y)
    {
      return   dl6(x) * l0(y);
    }

    static double simple_quad_l6_l0y(double x, double y)
    {
      return   l6(x) * dl0(y);
    }

    static double simple_quad_l6_l0xx(double x, double y)
    {
      return   d2l6(x) * l0(y);
    }

    static double simple_quad_l6_l0xy(double x, double y)
    {
      return   dl6(x) * dl0(y);
    }

    static double simple_quad_l6_l0yy(double x, double y)
    {
      return   l6(x) * d2l0(y);
    }

    static double simple_quad_l6_l1(double x, double y)
    {
      return   l6(x) * l1(y);
    }

    static double simple_quad_l6_l1x(double x, double y)
    {
      return   dl6(x) * l1(y);
    }

    static double simple_quad_l6_l1y(double x, double y)
    {
      return   l6(x) * dl1(y);
    }

    static double simple_quad_l6_l1xx(double x, double y)
    {
      return   d2l6(x) * l1(y);
    }

    static double simple_quad_l6_l1xy(double x, double y)
    {
      return   dl6(x) * dl1(y);
    }

    static double simple_quad_l6_l1yy(double x, double y)
    {
      return   l6(x) * d2l1(y);
    }

    static double simple_quad_l6_l2(double x, double y)
    {
      return   l6(x) * l2(y);
    }

    static double simple_quad_l6_l2x(double x, double y)
    {
      return   dl6(x) * l2(y);
    }

    static double simple_quad_l6_l2y(double x, double y)
    {
      return   l6(x) * dl2(y);
    }

    static double simple_quad_l6_l2xx(double x, double y)
    {
      return   d2l6(x) * l2(y);
    }

    static double simple_quad_l6_l2xy(double x, double y)
    {
      return   dl6(x) * dl2(y);
    }

    static double simple_quad_l6_l2yy(double x, double y)
    {
      return   l6(x) * d2l2(y);
    }

    static double simple_quad_l6_l3(double x, double y)
    {
      return   l6(x) * l3(y);
    }

    static double simple_quad_l6_l3x(double x, double y)
    {
      return   dl6(x) * l3(y);
    }

    static double simple_quad_l6_l3y(double x, double y)
    {
      return   l6(x) * dl3(y);
    }

    static double simple_quad_l6_l3xx(double x, double y)
    {
      return   d2l6(x) * l3(y);
    }

    static double simple_quad_l6_l3xy(double x, double y)
    {
      return   dl6(x) * dl3(y);
    }

    static double simple_quad_l6_l3yy(double x, double y)
    {
      return   l6(x) * d2l3(y);
    }

    static double simple_quad_l6_l4(double x, double y)
    {
      return   l6(x) * l4(y);
    }

    static double simple_quad_l6_l4x(double x, double y)
    {
      return   dl6(x) * l4(y);
    }

    static double simple_quad_l6_l4y(double x, double y)
    {
      return   l6(x) * dl4(y);
    }

    static double simple_quad_l6_l4xx(double x, double y)
    {
      return   d2l6(x) * l4(y);
    }

    static double simple_quad_l6_l4xy(double x, double y)
    {
      return   dl6(x) * dl4(y);
    }

    static double simple_quad_l6_l4yy(double x, double y)
    {
      return   l6(x) * d2l4(y);
    }

    static double simple_quad_l6_l5(double x, double y)
    {
      return   l6(x) * l5(y);
    }

    static double simple_quad_l6_l5x(double x, double y)
    {
      return   dl6(x) * l5(y);
    }

    static double simple_quad_l6_l5y(double x, double y)
    {
      return   l6(x) * dl5(y);
    }

    static double simple_quad_l6_l5xx(double x, double y)
    {
      return   d2l6(x) * l5(y);
    }

    static double simple_quad_l6_l5xy(double x, double y)
    {
      return   dl6(x) * dl5(y);
    }

    static double simple_quad_l6_l5yy(double x, double y)
    {
      return   l6(x) * d2l5(y);
    }

    static double simple_quad_l6_l6(double x, double y)
    {
      return   l6(x) * l6(y);
    }

    static double simple_quad_l6_l6x(double x, double y)
    {
      return   dl6(x) * l6(y);
    }

    static double simple_quad_l6_l6y(double x, double y)
    {
      return   l6(x) * dl6(y);
    }

    static double simple_quad_l6_l6xx(double x, double y)
    {
      return   d2l6(x) * l6(y);
    }

    static double simple_quad_l6_l6xy(double x, double y)
    {
      return   dl6(x) * dl6(y);
    }

    static double simple_quad_l6_l6yy(double x, double y)
    {
      return   l6(x) * d2l6(y);
    }

    static double simple_quad_l6_l7(double x, double y)
    {
      return   l6(x) * l7(y);
    }

    static double simple_quad_l6_l7x(double x, double y)
    {
      return   dl6(x) * l7(y);
    }

    static double simple_quad_l6_l7y(double x, double y)
    {
      return   l6(x) * dl7(y);
    }

    static double simple_quad_l6_l7xx(double x, double y)
    {
      return   d2l6(x) * l7(y);
    }

    static double simple_quad_l6_l7xy(double x, double y)
    {
      return   dl6(x) * dl7(y);
    }

    static double simple_quad_l6_l7yy(double x, double y)
    {
      return   l6(x) * d2l7(y);
    }

    static double simple_quad_l6_l8(double x, double y)
    {
      return   l6(x) * l8(y);
    }

    static double simple_quad_l6_l8x(double x, double y)
    {
      return   dl6(x) * l8(y);
    }

    static double simple_quad_l6_l8y(double x, double y)
    {
      return   l6(x) * dl8(y);
    }

    static double simple_quad_l6_l8xx(double x, double y)
    {
      return   d2l6(x) * l8(y);
    }

    static double simple_quad_l6_l8xy(double x, double y)
    {
      return   dl6(x) * dl8(y);
    }

    static double simple_quad_l6_l8yy(double x, double y)
    {
      return   l6(x) * d2l8(y);
    }

    static double simple_quad_l6_l9(double x, double y)
    {
      return   l6(x) * l9(y);
    }

    static double simple_quad_l6_l9x(double x, double y)
    {
      return   dl6(x) * l9(y);
    }

    static double simple_quad_l6_l9y(double x, double y)
    {
      return   l6(x) * dl9(y);
    }

    static double simple_quad_l6_l9xx(double x, double y)
    {
      return   d2l6(x) * l9(y);
    }

    static double simple_quad_l6_l9xy(double x, double y)
    {
      return   dl6(x) * dl9(y);
    }

    static double simple_quad_l6_l9yy(double x, double y)
    {
      return   l6(x) * d2l9(y);
    }

    static double simple_quad_l6_l10(double x, double y)
    {
      return   l6(x) * l10(y);
    }

    static double simple_quad_l6_l10x(double x, double y)
    {
      return   dl6(x) * l10(y);
    }

    static double simple_quad_l6_l10y(double x, double y)
    {
      return   l6(x) * dl10(y);
    }

    static double simple_quad_l6_l10xx(double x, double y)
    {
      return   d2l6(x) * l10(y);
    }

    static double simple_quad_l6_l10xy(double x, double y)
    {
      return   dl6(x) * dl10(y);
    }

    static double simple_quad_l6_l10yy(double x, double y)
    {
      return   l6(x) * d2l10(y);
    }

    static double simple_quad_l7_l0_0(double x, double y)
    {
      return   l7(x) * l0(y);
    }

    static double simple_quad_l7_l0x_0(double x, double y)
    {
      return   dl7(x) * l0(y);
    }

    static double simple_quad_l7_l0y_0(double x, double y)
    {
      return   l7(x) * dl0(y);
    }

    static double simple_quad_l7_l0xx_0(double x, double y)
    {
      return   d2l7(x) * l0(y);
    }

    static double simple_quad_l7_l0xy_0(double x, double y)
    {
      return   dl7(x) * dl0(y);
    }

    static double simple_quad_l7_l0yy_0(double x, double y)
    {
      return   l7(x) * d2l0(y);
    }

    static double simple_quad_l7_l0_1(double x, double y)
    {
      return -(  l7(x) * l0(y));
    }

    static double simple_quad_l7_l0x_1(double x, double y)
    {
      return -(  dl7(x) * l0(y));
    }

    static double simple_quad_l7_l0y_1(double x, double y)
    {
      return -(  l7(x) * dl0(y));
    }

    static double simple_quad_l7_l0xx_1(double x, double y)
    {
      return -(  d2l7(x) * l0(y));
    }

    static double simple_quad_l7_l0xy_1(double x, double y)
    {
      return -(  dl7(x) * dl0(y));
    }

    static double simple_quad_l7_l0yy_1(double x, double y)
    {
      return -(  l7(x) * d2l0(y));
    }

    static double simple_quad_l7_l1_0(double x, double y)
    {
      return - l7(x) * l1(y);
    }

    static double simple_quad_l7_l1x_0(double x, double y)
    {
      return - dl7(x) * l1(y);
    }

    static double simple_quad_l7_l1y_0(double x, double y)
    {
      return - l7(x) * dl1(y);
    }

    static double simple_quad_l7_l1xx_0(double x, double y)
    {
      return - d2l7(x) * l1(y);
    }

    static double simple_quad_l7_l1xy_0(double x, double y)
    {
      return - dl7(x) * dl1(y);
    }

    static double simple_quad_l7_l1yy_0(double x, double y)
    {
      return - l7(x) * d2l1(y);
    }

    static double simple_quad_l7_l1_1(double x, double y)
    {
      return -(- l7(x) * l1(y));
    }

    static double simple_quad_l7_l1x_1(double x, double y)
    {
      return -(- dl7(x) * l1(y));
    }

    static double simple_quad_l7_l1y_1(double x, double y)
    {
      return -(- l7(x) * dl1(y));
    }

    static double simple_quad_l7_l1xx_1(double x, double y)
    {
      return -(- d2l7(x) * l1(y));
    }

    static double simple_quad_l7_l1xy_1(double x, double y)
    {
      return -(- dl7(x) * dl1(y));
    }

    static double simple_quad_l7_l1yy_1(double x, double y)
    {
      return -(- l7(x) * d2l1(y));
    }

    static double simple_quad_l7_l2(double x, double y)
    {
      return   l7(x) * l2(y);
    }

    static double simple_quad_l7_l2x(double x, double y)
    {
      return   dl7(x) * l2(y);
    }

    static double simple_quad_l7_l2y(double x, double y)
    {
      return   l7(x) * dl2(y);
    }

    static double simple_quad_l7_l2xx(double x, double y)
    {
      return   d2l7(x) * l2(y);
    }

    static double simple_quad_l7_l2xy(double x, double y)
    {
      return   dl7(x) * dl2(y);
    }

    static double simple_quad_l7_l2yy(double x, double y)
    {
      return   l7(x) * d2l2(y);
    }

    static double simple_quad_l7_l3(double x, double y)
    {
      return   l7(x) * l3(y);
    }

    static double simple_quad_l7_l3x(double x, double y)
    {
      return   dl7(x) * l3(y);
    }

    static double simple_quad_l7_l3y(double x, double y)
    {
      return   l7(x) * dl3(y);
    }

    static double simple_quad_l7_l3xx(double x, double y)
    {
      return   d2l7(x) * l3(y);
    }

    static double simple_quad_l7_l3xy(double x, double y)
    {
      return   dl7(x) * dl3(y);
    }

    static double simple_quad_l7_l3yy(double x, double y)
    {
      return   l7(x) * d2l3(y);
    }

    static double simple_quad_l7_l4(double x, double y)
    {
      return   l7(x) * l4(y);
    }

    static double simple_quad_l7_l4x(double x, double y)
    {
      return   dl7(x) * l4(y);
    }

    static double simple_quad_l7_l4y(double x, double y)
    {
      return   l7(x) * dl4(y);
    }

    static double simple_quad_l7_l4xx(double x, double y)
    {
      return   d2l7(x) * l4(y);
    }

    static double simple_quad_l7_l4xy(double x, double y)
    {
      return   dl7(x) * dl4(y);
    }

    static double simple_quad_l7_l4yy(double x, double y)
    {
      return   l7(x) * d2l4(y);
    }

    static double simple_quad_l7_l5(double x, double y)
    {
      return   l7(x) * l5(y);
    }

    static double simple_quad_l7_l5x(double x, double y)
    {
      return   dl7(x) * l5(y);
    }

    static double simple_quad_l7_l5y(double x, double y)
    {
      return   l7(x) * dl5(y);
    }

    static double simple_quad_l7_l5xx(double x, double y)
    {
      return   d2l7(x) * l5(y);
    }

    static double simple_quad_l7_l5xy(double x, double y)
    {
      return   dl7(x) * dl5(y);
    }

    static double simple_quad_l7_l5yy(double x, double y)
    {
      return   l7(x) * d2l5(y);
    }

    static double simple_quad_l7_l6(double x, double y)
    {
      return   l7(x) * l6(y);
    }

    static double simple_quad_l7_l6x(double x, double y)
    {
      return   dl7(x) * l6(y);
    }

    static double simple_quad_l7_l6y(double x, double y)
    {
      return   l7(x) * dl6(y);
    }

    static double simple_quad_l7_l6xx(double x, double y)
    {
      return   d2l7(x) * l6(y);
    }

    static double simple_quad_l7_l6xy(double x, double y)
    {
      return   dl7(x) * dl6(y);
    }

    static double simple_quad_l7_l6yy(double x, double y)
    {
      return   l7(x) * d2l6(y);
    }

    static double simple_quad_l7_l7(double x, double y)
    {
      return   l7(x) * l7(y);
    }

    static double simple_quad_l7_l7x(double x, double y)
    {
      return   dl7(x) * l7(y);
    }

    static double simple_quad_l7_l7y(double x, double y)
    {
      return   l7(x) * dl7(y);
    }

    static double simple_quad_l7_l7xx(double x, double y)
    {
      return   d2l7(x) * l7(y);
    }

    static double simple_quad_l7_l7xy(double x, double y)
    {
      return   dl7(x) * dl7(y);
    }

    static double simple_quad_l7_l7yy(double x, double y)
    {
      return   l7(x) * d2l7(y);
    }

    static double simple_quad_l7_l8(double x, double y)
    {
      return   l7(x) * l8(y);
    }

    static double simple_quad_l7_l8x(double x, double y)
    {
      return   dl7(x) * l8(y);
    }

    static double simple_quad_l7_l8y(double x, double y)
    {
      return   l7(x) * dl8(y);
    }

    static double simple_quad_l7_l8xx(double x, double y)
    {
      return   d2l7(x) * l8(y);
    }

    static double simple_quad_l7_l8xy(double x, double y)
    {
      return   dl7(x) * dl8(y);
    }

    static double simple_quad_l7_l8yy(double x, double y)
    {
      return   l7(x) * d2l8(y);
    }

    static double simple_quad_l7_l9(double x, double y)
    {
      return   l7(x) * l9(y);
    }

    static double simple_quad_l7_l9x(double x, double y)
    {
      return   dl7(x) * l9(y);
    }

    static double simple_quad_l7_l9y(double x, double y)
    {
      return   l7(x) * dl9(y);
    }

    static double simple_quad_l7_l9xx(double x, double y)
    {
      return   d2l7(x) * l9(y);
    }

    static double simple_quad_l7_l9xy(double x, double y)
    {
      return   dl7(x) * dl9(y);
    }

    static double simple_quad_l7_l9yy(double x, double y)
    {
      return   l7(x) * d2l9(y);
    }

    static double simple_quad_l7_l10(double x, double y)
    {
      return   l7(x) * l10(y);
    }

    static double simple_quad_l7_l10x(double x, double y)
    {
      return   dl7(x) * l10(y);
    }

    static double simple_quad_l7_l10y(double x, double y)
    {
      return   l7(x) * dl10(y);
    }

    static double simple_quad_l7_l10xx(double x, double y)
    {
      return   d2l7(x) * l10(y);
    }

    static double simple_quad_l7_l10xy(double x, double y)
    {
      return   dl7(x) * dl10(y);
    }

    static double simple_quad_l7_l10yy(double x, double y)
    {
      return   l7(x) * d2l10(y);
    }

    static double simple_quad_l8_l0(double x, double y)
    {
      return   l8(x) * l0(y);
    }

    static double simple_quad_l8_l0x(double x, double y)
    {
      return   dl8(x) * l0(y);
    }

    static double simple_quad_l8_l0y(double x, double y)
    {
      return   l8(x) * dl0(y);
    }

    static double simple_quad_l8_l0xx(double x, double y)
    {
      return   d2l8(x) * l0(y);
    }

    static double simple_quad_l8_l0xy(double x, double y)
    {
      return   dl8(x) * dl0(y);
    }

    static double simple_quad_l8_l0yy(double x, double y)
    {
      return   l8(x) * d2l0(y);
    }

    static double simple_quad_l8_l1(double x, double y)
    {
      return   l8(x) * l1(y);
    }

    static double simple_quad_l8_l1x(double x, double y)
    {
      return   dl8(x) * l1(y);
    }

    static double simple_quad_l8_l1y(double x, double y)
    {
      return   l8(x) * dl1(y);
    }

    static double simple_quad_l8_l1xx(double x, double y)
    {
      return   d2l8(x) * l1(y);
    }

    static double simple_quad_l8_l1xy(double x, double y)
    {
      return   dl8(x) * dl1(y);
    }

    static double simple_quad_l8_l1yy(double x, double y)
    {
      return   l8(x) * d2l1(y);
    }

    static double simple_quad_l8_l2(double x, double y)
    {
      return   l8(x) * l2(y);
    }

    static double simple_quad_l8_l2x(double x, double y)
    {
      return   dl8(x) * l2(y);
    }

    static double simple_quad_l8_l2y(double x, double y)
    {
      return   l8(x) * dl2(y);
    }

    static double simple_quad_l8_l2xx(double x, double y)
    {
      return   d2l8(x) * l2(y);
    }

    static double simple_quad_l8_l2xy(double x, double y)
    {
      return   dl8(x) * dl2(y);
    }

    static double simple_quad_l8_l2yy(double x, double y)
    {
      return   l8(x) * d2l2(y);
    }

    static double simple_quad_l8_l3(double x, double y)
    {
      return   l8(x) * l3(y);
    }

    static double simple_quad_l8_l3x(double x, double y)
    {
      return   dl8(x) * l3(y);
    }

    static double simple_quad_l8_l3y(double x, double y)
    {
      return   l8(x) * dl3(y);
    }

    static double simple_quad_l8_l3xx(double x, double y)
    {
      return   d2l8(x) * l3(y);
    }

    static double simple_quad_l8_l3xy(double x, double y)
    {
      return   dl8(x) * dl3(y);
    }

    static double simple_quad_l8_l3yy(double x, double y)
    {
      return   l8(x) * d2l3(y);
    }

    static double simple_quad_l8_l4(double x, double y)
    {
      return   l8(x) * l4(y);
    }

    static double simple_quad_l8_l4x(double x, double y)
    {
      return   dl8(x) * l4(y);
    }

    static double simple_quad_l8_l4y(double x, double y)
    {
      return   l8(x) * dl4(y);
    }

    static double simple_quad_l8_l4xx(double x, double y)
    {
      return   d2l8(x) * l4(y);
    }

    static double simple_quad_l8_l4xy(double x, double y)
    {
      return   dl8(x) * dl4(y);
    }

    static double simple_quad_l8_l4yy(double x, double y)
    {
      return   l8(x) * d2l4(y);
    }

    static double simple_quad_l8_l5(double x, double y)
    {
      return   l8(x) * l5(y);
    }

    static double simple_quad_l8_l5x(double x, double y)
    {
      return   dl8(x) * l5(y);
    }

    static double simple_quad_l8_l5y(double x, double y)
    {
      return   l8(x) * dl5(y);
    }

    static double simple_quad_l8_l5xx(double x, double y)
    {
      return   d2l8(x) * l5(y);
    }

    static double simple_quad_l8_l5xy(double x, double y)
    {
      return   dl8(x) * dl5(y);
    }

    static double simple_quad_l8_l5yy(double x, double y)
    {
      return   l8(x) * d2l5(y);
    }

    static double simple_quad_l8_l6(double x, double y)
    {
      return   l8(x) * l6(y);
    }

    static double simple_quad_l8_l6x(double x, double y)
    {
      return   dl8(x) * l6(y);
    }

    static double simple_quad_l8_l6y(double x, double y)
    {
      return   l8(x) * dl6(y);
    }

    static double simple_quad_l8_l6xx(double x, double y)
    {
      return   d2l8(x) * l6(y);
    }

    static double simple_quad_l8_l6xy(double x, double y)
    {
      return   dl8(x) * dl6(y);
    }

    static double simple_quad_l8_l6yy(double x, double y)
    {
      return   l8(x) * d2l6(y);
    }

    static double simple_quad_l8_l7(double x, double y)
    {
      return   l8(x) * l7(y);
    }

    static double simple_quad_l8_l7x(double x, double y)
    {
      return   dl8(x) * l7(y);
    }

    static double simple_quad_l8_l7y(double x, double y)
    {
      return   l8(x) * dl7(y);
    }

    static double simple_quad_l8_l7xx(double x, double y)
    {
      return   d2l8(x) * l7(y);
    }

    static double simple_quad_l8_l7xy(double x, double y)
    {
      return   dl8(x) * dl7(y);
    }

    static double simple_quad_l8_l7yy(double x, double y)
    {
      return   l8(x) * d2l7(y);
    }

    static double simple_quad_l8_l8(double x, double y)
    {
      return   l8(x) * l8(y);
    }

    static double simple_quad_l8_l8x(double x, double y)
    {
      return   dl8(x) * l8(y);
    }

    static double simple_quad_l8_l8y(double x, double y)
    {
      return   l8(x) * dl8(y);
    }

    static double simple_quad_l8_l8xx(double x, double y)
    {
      return   d2l8(x) * l8(y);
    }

    static double simple_quad_l8_l8xy(double x, double y)
    {
      return   dl8(x) * dl8(y);
    }

    static double simple_quad_l8_l8yy(double x, double y)
    {
      return   l8(x) * d2l8(y);
    }

    static double simple_quad_l8_l9(double x, double y)
    {
      return   l8(x) * l9(y);
    }

    static double simple_quad_l8_l9x(double x, double y)
    {
      return   dl8(x) * l9(y);
    }

    static double simple_quad_l8_l9y(double x, double y)
    {
      return   l8(x) * dl9(y);
    }

    static double simple_quad_l8_l9xx(double x, double y)
    {
      return   d2l8(x) * l9(y);
    }

    static double simple_quad_l8_l9xy(double x, double y)
    {
      return   dl8(x) * dl9(y);
    }

    static double simple_quad_l8_l9yy(double x, double y)
    {
      return   l8(x) * d2l9(y);
    }

    static double simple_quad_l8_l10(double x, double y)
    {
      return   l8(x) * l10(y);
    }

    static double simple_quad_l8_l10x(double x, double y)
    {
      return   dl8(x) * l10(y);
    }

    static double simple_quad_l8_l10y(double x, double y)
    {
      return   l8(x) * dl10(y);
    }

    static double simple_quad_l8_l10xx(double x, double y)
    {
      return   d2l8(x) * l10(y);
    }

    static double simple_quad_l8_l10xy(double x, double y)
    {
      return   dl8(x) * dl10(y);
    }

    static double simple_quad_l8_l10yy(double x, double y)
    {
      return   l8(x) * d2l10(y);
    }

    static double simple_quad_l9_l0_0(double x, double y)
    {
      return   l9(x) * l0(y);
    }

    static double simple_quad_l9_l0x_0(double x, double y)
    {
      return   dl9(x) * l0(y);
    }

    static double simple_quad_l9_l0y_0(double x, double y)
    {
      return   l9(x) * dl0(y);
    }

    static double simple_quad_l9_l0xx_0(double x, double y)
    {
      return   d2l9(x) * l0(y);
    }

    static double simple_quad_l9_l0xy_0(double x, double y)
    {
      return   dl9(x) * dl0(y);
    }

    static double simple_quad_l9_l0yy_0(double x, double y)
    {
      return   l9(x) * d2l0(y);
    }

    static double simple_quad_l9_l0_1(double x, double y)
    {
      return -(  l9(x) * l0(y));
    }

    static double simple_quad_l9_l0x_1(double x, double y)
    {
      return -(  dl9(x) * l0(y));
    }

    static double simple_quad_l9_l0y_1(double x, double y)
    {
      return -(  l9(x) * dl0(y));
    }

    static double simple_quad_l9_l0xx_1(double x, double y)
    {
      return -(  d2l9(x) * l0(y));
    }

    static double simple_quad_l9_l0xy_1(double x, double y)
    {
      return -(  dl9(x) * dl0(y));
    }

    static double simple_quad_l9_l0yy_1(double x, double y)
    {
      return -(  l9(x) * d2l0(y));
    }

    static double simple_quad_l9_l1_0(double x, double y)
    {
      return - l9(x) * l1(y);
    }

    static double simple_quad_l9_l1x_0(double x, double y)
    {
      return - dl9(x) * l1(y);
    }

    static double simple_quad_l9_l1y_0(double x, double y)
    {
      return - l9(x) * dl1(y);
    }

    static double simple_quad_l9_l1xx_0(double x, double y)
    {
      return - d2l9(x) * l1(y);
    }

    static double simple_quad_l9_l1xy_0(double x, double y)
    {
      return - dl9(x) * dl1(y);
    }

    static double simple_quad_l9_l1yy_0(double x, double y)
    {
      return - l9(x) * d2l1(y);
    }

    static double simple_quad_l9_l1_1(double x, double y)
    {
      return -(- l9(x) * l1(y));
    }

    static double simple_quad_l9_l1x_1(double x, double y)
    {
      return -(- dl9(x) * l1(y));
    }

    static double simple_quad_l9_l1y_1(double x, double y)
    {
      return -(- l9(x) * dl1(y));
    }

    static double simple_quad_l9_l1xx_1(double x, double y)
    {
      return -(- d2l9(x) * l1(y));
    }

    static double simple_quad_l9_l1xy_1(double x, double y)
    {
      return -(- dl9(x) * dl1(y));
    }

    static double simple_quad_l9_l1yy_1(double x, double y)
    {
      return -(- l9(x) * d2l1(y));
    }

    static double simple_quad_l9_l2(double x, double y)
    {
      return   l9(x) * l2(y);
    }

    static double simple_quad_l9_l2x(double x, double y)
    {
      return   dl9(x) * l2(y);
    }

    static double simple_quad_l9_l2y(double x, double y)
    {
      return   l9(x) * dl2(y);
    }

    static double simple_quad_l9_l2xx(double x, double y)
    {
      return   d2l9(x) * l2(y);
    }

    static double simple_quad_l9_l2xy(double x, double y)
    {
      return   dl9(x) * dl2(y);
    }

    static double simple_quad_l9_l2yy(double x, double y)
    {
      return   l9(x) * d2l2(y);
    }

    static double simple_quad_l9_l3(double x, double y)
    {
      return   l9(x) * l3(y);
    }

    static double simple_quad_l9_l3x(double x, double y)
    {
      return   dl9(x) * l3(y);
    }

    static double simple_quad_l9_l3y(double x, double y)
    {
      return   l9(x) * dl3(y);
    }

    static double simple_quad_l9_l3xx(double x, double y)
    {
      return   d2l9(x) * l3(y);
    }

    static double simple_quad_l9_l3xy(double x, double y)
    {
      return   dl9(x) * dl3(y);
    }

    static double simple_quad_l9_l3yy(double x, double y)
    {
      return   l9(x) * d2l3(y);
    }

    static double simple_quad_l9_l4(double x, double y)
    {
      return   l9(x) * l4(y);
    }

    static double simple_quad_l9_l4x(double x, double y)
    {
      return   dl9(x) * l4(y);
    }

    static double simple_quad_l9_l4y(double x, double y)
    {
      return   l9(x) * dl4(y);
    }

    static double simple_quad_l9_l4xx(double x, double y)
    {
      return   d2l9(x) * l4(y);
    }

    static double simple_quad_l9_l4xy(double x, double y)
    {
      return   dl9(x) * dl4(y);
    }

    static double simple_quad_l9_l4yy(double x, double y)
    {
      return   l9(x) * d2l4(y);
    }

    static double simple_quad_l9_l5(double x, double y)
    {
      return   l9(x) * l5(y);
    }

    static double simple_quad_l9_l5x(double x, double y)
    {
      return   dl9(x) * l5(y);
    }

    static double simple_quad_l9_l5y(double x, double y)
    {
      return   l9(x) * dl5(y);
    }

    static double simple_quad_l9_l5xx(double x, double y)
    {
      return   d2l9(x) * l5(y);
    }

    static double simple_quad_l9_l5xy(double x, double y)
    {
      return   dl9(x) * dl5(y);
    }

    static double simple_quad_l9_l5yy(double x, double y)
    {
      return   l9(x) * d2l5(y);
    }

    static double simple_quad_l9_l6(double x, double y)
    {
      return   l9(x) * l6(y);
    }

    static double simple_quad_l9_l6x(double x, double y)
    {
      return   dl9(x) * l6(y);
    }

    static double simple_quad_l9_l6y(double x, double y)
    {
      return   l9(x) * dl6(y);
    }

    static double simple_quad_l9_l6xx(double x, double y)
    {
      return   d2l9(x) * l6(y);
    }

    static double simple_quad_l9_l6xy(double x, double y)
    {
      return   dl9(x) * dl6(y);
    }

    static double simple_quad_l9_l6yy(double x, double y)
    {
      return   l9(x) * d2l6(y);
    }

    static double simple_quad_l9_l7(double x, double y)
    {
      return   l9(x) * l7(y);
    }

    static double simple_quad_l9_l7x(double x, double y)
    {
      return   dl9(x) * l7(y);
    }

    static double simple_quad_l9_l7y(double x, double y)
    {
      return   l9(x) * dl7(y);
    }

    static double simple_quad_l9_l7xx(double x, double y)
    {
      return   d2l9(x) * l7(y);
    }

    static double simple_quad_l9_l7xy(double x, double y)
    {
      return   dl9(x) * dl7(y);
    }

    static double simple_quad_l9_l7yy(double x, double y)
    {
      return   l9(x) * d2l7(y);
    }

    static double simple_quad_l9_l8(double x, double y)
    {
      return   l9(x) * l8(y);
    }

    static double simple_quad_l9_l8x(double x, double y)
    {
      return   dl9(x) * l8(y);
    }

    static double simple_quad_l9_l8y(double x, double y)
    {
      return   l9(x) * dl8(y);
    }

    static double simple_quad_l9_l8xx(double x, double y)
    {
      return   d2l9(x) * l8(y);
    }

    static double simple_quad_l9_l8xy(double x, double y)
    {
      return   dl9(x) * dl8(y);
    }

    static double simple_quad_l9_l8yy(double x, double y)
    {
      return   l9(x) * d2l8(y);
    }

    static double simple_quad_l9_l9(double x, double y)
    {
      return   l9(x) * l9(y);
    }

    static double simple_quad_l9_l9x(double x, double y)
    {
      return   dl9(x) * l9(y);
    }

    static double simple_quad_l9_l9y(double x, double y)
    {
      return   l9(x) * dl9(y);
    }

    static double simple_quad_l9_l9xx(double x, double y)
    {
      return   d2l9(x) * l9(y);
    }

    static double simple_quad_l9_l9xy(double x, double y)
    {
      return   dl9(x) * dl9(y);
    }

    static double simple_quad_l9_l9yy(double x, double y)
    {
      return   l9(x) * d2l9(y);
    }

    static double simple_quad_l9_l10(double x, double y)
    {
      return   l9(x) * l10(y);
    }

    static double simple_quad_l9_l10x(double x, double y)
    {
      return   dl9(x) * l10(y);
    }

    static double simple_quad_l9_l10y(double x, double y)
    {
      return   l9(x) * dl10(y);
    }

    static double simple_quad_l9_l10xx(double x, double y)
    {
      return   d2l9(x) * l10(y);
    }

    static double simple_quad_l9_l10xy(double x, double y)
    {
      return   dl9(x) * dl10(y);
    }

    static double simple_quad_l9_l10yy(double x, double y)
    {
      return   l9(x) * d2l10(y);
    }

    static double simple_quad_l10_l0(double x, double y)
    {
      return   l10(x) * l0(y);
    }

    static double simple_quad_l10_l0x(double x, double y)
    {
      return   dl10(x) * l0(y);
    }

    static double simple_quad_l10_l0y(double x, double y)
    {
      return   l10(x) * dl0(y);
    }

    static double simple_quad_l10_l0xx(double x, double y)
    {
      return   d2l10(x) * l0(y);
    }

    static double simple_quad_l10_l0xy(double x, double y)
    {
      return   dl10(x) * dl0(y);
    }

    static double simple_quad_l10_l0yy(double x, double y)
    {
      return   l10(x) * d2l0(y);
    }

    static double simple_quad_l10_l1(double x, double y)
    {
      return   l10(x) * l1(y);
    }

    static double simple_quad_l10_l1x(double x, double y)
    {
      return   dl10(x) * l1(y);
    }

    static double simple_quad_l10_l1y(double x, double y)
    {
      return   l10(x) * dl1(y);
    }

    static double simple_quad_l10_l1xx(double x, double y)
    {
      return   d2l10(x) * l1(y);
    }

    static double simple_quad_l10_l1xy(double x, double y)
    {
      return   dl10(x) * dl1(y);
    }

    static double simple_quad_l10_l1yy(double x, double y)
    {
      return   l10(x) * d2l1(y);
    }

    static double simple_quad_l10_l2(double x, double y)
    {
      return   l10(x) * l2(y);
    }

    static double simple_quad_l10_l2x(double x, double y)
    {
      return   dl10(x) * l2(y);
    }

    static double simple_quad_l10_l2y(double x, double y)
    {
      return   l10(x) * dl2(y);
    }

    static double simple_quad_l10_l2xx(double x, double y)
    {
      return   d2l10(x) * l2(y);
    }

    static double simple_quad_l10_l2xy(double x, double y)
    {
      return   dl10(x) * dl2(y);
    }

    static double simple_quad_l10_l2yy(double x, double y)
    {
      return   l10(x) * d2l2(y);
    }

    static double simple_quad_l10_l3(double x, double y)
    {
      return   l10(x) * l3(y);
    }

    static double simple_quad_l10_l3x(double x, double y)
    {
      return   dl10(x) * l3(y);
    }

    static double simple_quad_l10_l3y(double x, double y)
    {
      return   l10(x) * dl3(y);
    }

    static double simple_quad_l10_l3xx(double x, double y)
    {
      return   d2l10(x) * l3(y);
    }

    static double simple_quad_l10_l3xy(double x, double y)
    {
      return   dl10(x) * dl3(y);
    }

    static double simple_quad_l10_l3yy(double x, double y)
    {
      return   l10(x) * d2l3(y);
    }

    static double simple_quad_l10_l4(double x, double y)
    {
      return   l10(x) * l4(y);
    }

    static double simple_quad_l10_l4x(double x, double y)
    {
      return   dl10(x) * l4(y);
    }

    static double simple_quad_l10_l4y(double x, double y)
    {
      return   l10(x) * dl4(y);
    }

    static double simple_quad_l10_l4xx(double x, double y)
    {
      return   d2l10(x) * l4(y);
    }

    static double simple_quad_l10_l4xy(double x, double y)
    {
      return   dl10(x) * dl4(y);
    }

    static double simple_quad_l10_l4yy(double x, double y)
    {
      return   l10(x) * d2l4(y);
    }

    static double simple_quad_l10_l5(double x, double y)
    {
      return   l10(x) * l5(y);
    }

    static double simple_quad_l10_l5x(double x, double y)
    {
      return   dl10(x) * l5(y);
    }

    static double simple_quad_l10_l5y(double x, double y)
    {
      return   l10(x) * dl5(y);
    }

    static double simple_quad_l10_l5xx(double x, double y)
    {
      return   d2l10(x) * l5(y);
    }

    static double simple_quad_l10_l5xy(double x, double y)
    {
      return   dl10(x) * dl5(y);
    }

    static double simple_quad_l10_l5yy(double x, double y)
    {
      return   l10(x) * d2l5(y);
    }

    static double simple_quad_l10_l6(double x, double y)
    {
      return   l10(x) * l6(y);
    }

    static double simple_quad_l10_l6x(double x, double y)
    {
      return   dl10(x) * l6(y);
    }

    static double simple_quad_l10_l6y(double x, double y)
    {
      return   l10(x) * dl6(y);
    }

    static double simple_quad_l10_l6xx(double x, double y)
    {
      return   d2l10(x) * l6(y);
    }

    static double simple_quad_l10_l6xy(double x, double y)
    {
      return   dl10(x) * dl6(y);
    }

    static double simple_quad_l10_l6yy(double x, double y)
    {
      return   l10(x) * d2l6(y);
    }

    static double simple_quad_l10_l7(double x, double y)
    {
      return   l10(x) * l7(y);
    }

    static double simple_quad_l10_l7x(double x, double y)
    {
      return   dl10(x) * l7(y);
    }

    static double simple_quad_l10_l7y(double x, double y)
    {
      return   l10(x) * dl7(y);
    }

    static double simple_quad_l10_l7xx(double x, double y)
    {
      return   d2l10(x) * l7(y);
    }

    static double simple_quad_l10_l7xy(double x, double y)
    {
      return   dl10(x) * dl7(y);
    }

    static double simple_quad_l10_l7yy(double x, double y)
    {
      return   l10(x) * d2l7(y);
    }

    static double simple_quad_l10_l8(double x, double y)
    {
      return   l10(x) * l8(y);
    }

    static double simple_quad_l10_l8x(double x, double y)
    {
      return   dl10(x) * l8(y);
    }

    static double simple_quad_l10_l8y(double x, double y)
    {
      return   l10(x) * dl8(y);
    }

    static double simple_quad_l10_l8xx(double x, double y)
    {
      return   d2l10(x) * l8(y);
    }

    static double simple_quad_l10_l8xy(double x, double y)
    {
      return   dl10(x) * dl8(y);
    }

    static double simple_quad_l10_l8yy(double x, double y)
    {
      return   l10(x) * d2l8(y);
    }

    static double simple_quad_l10_l9(double x, double y)
    {
      return   l10(x) * l9(y);
    }

    static double simple_quad_l10_l9x(double x, double y)
    {
      return   dl10(x) * l9(y);
    }

    static double simple_quad_l10_l9y(double x, double y)
    {
      return   l10(x) * dl9(y);
    }

    static double simple_quad_l10_l9xx(double x, double y)
    {
      return   d2l10(x) * l9(y);
    }

    static double simple_quad_l10_l9xy(double x, double y)
    {
      return   dl10(x) * dl9(y);
    }

    static double simple_quad_l10_l9yy(double x, double y)
    {
      return   l10(x) * d2l9(y);
    }

    static double simple_quad_l10_l10(double x, double y)
    {
      return   l10(x) * l10(y);
    }

    static double simple_quad_l10_l10x(double x, double y)
    {
      return   dl10(x) * l10(y);
    }

    static double simple_quad_l10_l10y(double x, double y)
    {
      return   l10(x) * dl10(y);
    }

    static double simple_quad_l10_l10xx(double x, double y)
    {
      return   d2l10(x) * l10(y);
    }

    static double simple_quad_l10_l10xy(double x, double y)
    {
      return   dl10(x) * dl10(y);
    }

    static double simple_quad_l10_l10yy(double x, double y)
    {
      return   l10(x) * d2l10(y);
    }

    static Shapeset::shape_fn_t simple_quad_fn[] =
    {
      simple_quad_l0_l0,   simple_quad_l0_l1,   simple_quad_l0_l2,   simple_quad_l0_l3_0, simple_quad_l0_l3_1,
      simple_quad_l0_l4,   simple_quad_l0_l5_0, simple_quad_l0_l5_1, simple_quad_l0_l6,   simple_quad_l0_l7_0,
      simple_quad_l0_l7_1, simple_quad_l0_l8,   simple_quad_l0_l9_0, simple_quad_l0_l9_1, simple_quad_l0_l10,
      simple_quad_l1_l0,   simple_quad_l1_l1,   simple_quad_l1_l2,   simple_quad_l1_l3_0, simple_quad_l1_l3_1,
      simple_quad_l1_l4,   simple_quad_l1_l5_0, simple_quad_l1_l5_1, simple_quad_l1_l6,   simple_quad_l1_l7_0,
      simple_quad_l1_l7_1, simple_quad_l1_l8,   simple_quad_l1_l9_0, simple_quad_l1_l9_1, simple_quad_l1_l10,
      simple_quad_l2_l0,   simple_quad_l2_l1,   simple_quad_l2_l2,   simple_quad_l2_l3,   simple_quad_l2_l4,
      simple_quad_l2_l5,   simple_quad_l2_l6,   simple_quad_l2_l7,   simple_quad_l2_l8,   simple_quad_l2_l9,
      simple_quad_l2_l10,   simple_quad_l3_l0_0, simple_quad_l3_l0_1, simple_quad_l3_l1_0, simple_quad_l3_l1_1,
      simple_quad_l3_l2,   simple_quad_l3_l3,   simple_quad_l3_l4,   simple_quad_l3_l5,   simple_quad_l3_l6,
      simple_quad_l3_l7,   simple_quad_l3_l8,   simple_quad_l3_l9,   simple_quad_l3_l10,   simple_quad_l4_l0,
      simple_quad_l4_l1,   simple_quad_l4_l2,   simple_quad_l4_l3,   simple_quad_l4_l4,   simple_quad_l4_l5,
      simple_quad_l4_l6,   simple_quad_l4_l7,   simple_quad_l4_l8,   simple_quad_l4_l9,   simple_quad_l4_l10,
      simple_quad_l5_l0_0, simple_quad_l5_l0_1, simple_quad_l5_l1_0, simple_quad_l5_l1_1, simple_quad_l5_l2,
      simple_quad_l5_l3,   simple_quad_l5_l4,   simple_quad_l5_l5,   simple_quad_l5_l6,   simple_quad_l5_l7,
      simple_quad_l5_l8,   simple_quad_l5_l9,   simple_quad_l5_l10,   simple_quad_l6_l0,   simple_quad_l6_l1,
      simple_quad_l6_l2,   simple_quad_l6_l3,   simple_quad_l6_l4,   simple_quad_l6_l5,   simple_quad_l6_l6,
      simple_quad_l6_l7,   simple_quad_l6_l8,   simple_quad_l6_l9,   simple_quad_l6_l10,   simple_quad_l7_l0_0,
      simple_quad_l7_l0_1, simple_quad_l7_l1_0, simple_quad_l7_l1_1, simple_quad_l7_l2,   simple_quad_l7_l3,
      simple_quad_l7_l4,   simple_quad_l7_l5,   simple_quad_l7_l6,   simple_quad_l7_l7,   simple_quad_l7_l8,
      simple_quad_l7_l9,   simple_quad_l7_l10,   simple_quad_l8_l0,   simple_quad_l8_l1,   simple_quad_l8_l2,
      simple_quad_l8_l3,   simple_quad_l8_l4,   simple_quad_l8_l5,   simple_quad_l8_l6,   simple_quad_l8_l7,
      simple_quad_l8_l8,   simple_quad_l8_l9,   simple_quad_l8_l10,   simple_quad_l9_l0_0, simple_quad_l9_l0_1,
      simple_quad_l9_l1_0, simple_quad_l9_l1_1, simple_quad_l9_l2,   simple_quad_l9_l3,   simple_quad_l9_l4,
      simple_quad_l9_l5,   simple_quad_l9_l6,   simple_quad_l9_l7,   simple_quad_l9_l8,   simple_quad_l9_l9,
      simple_quad_l9_l10,   simple_quad_l10_l0,   simple_quad_l10_l1,   simple_quad_l10_l2,   simple_quad_l10_l3,
      simple_quad_l10_l4,   simple_quad_l10_l5,   simple_quad_l10_l6,   simple_quad_l10_l7,   simple_quad_l10_l8,
      simple_quad_l10_l9,   simple_quad_l10_l10,
    };
    static Shapeset::shape_fn_t simple_quad_fn_dx[] =
    {
      simple_quad_l0_l0x,   simple_quad_l0_l1x,   simple_quad_l0_l2x,   simple_quad_l0_l3x_0, simple_quad_l0_l3x_1,
      simple_quad_l0_l4x,   simple_quad_l0_l5x_0, simple_quad_l0_l5x_1, simple_quad_l0_l6x,   simple_quad_l0_l7x_0,
      simple_quad_l0_l7x_1, simple_quad_l0_l8x,   simple_quad_l0_l9x_0, simple_quad_l0_l9x_1, simple_quad_l0_l10x,
      simple_quad_l1_l0x,   simple_quad_l1_l1x,   simple_quad_l1_l2x,   simple_quad_l1_l3x_0, simple_quad_l1_l3x_1,
      simple_quad_l1_l4x,   simple_quad_l1_l5x_0, simple_quad_l1_l5x_1, simple_quad_l1_l6x,   simple_quad_l1_l7x_0,
      simple_quad_l1_l7x_1, simple_quad_l1_l8x,   simple_quad_l1_l9x_0, simple_quad_l1_l9x_1, simple_quad_l1_l10x,
      simple_quad_l2_l0x,   simple_quad_l2_l1x,   simple_quad_l2_l2x,   simple_quad_l2_l3x,   simple_quad_l2_l4x,
      simple_quad_l2_l5x,   simple_quad_l2_l6x,   simple_quad_l2_l7x,   simple_quad_l2_l8x,   simple_quad_l2_l9x,
      simple_quad_l2_l10x,   simple_quad_l3_l0x_0, simple_quad_l3_l0x_1, simple_quad_l3_l1x_0, simple_quad_l3_l1x_1,
      simple_quad_l3_l2x,   simple_quad_l3_l3x,   simple_quad_l3_l4x,   simple_quad_l3_l5x,   simple_quad_l3_l6x,
      simple_quad_l3_l7x,   simple_quad_l3_l8x,   simple_quad_l3_l9x,   simple_quad_l3_l10x,   simple_quad_l4_l0x,
      simple_quad_l4_l1x,   simple_quad_l4_l2x,   simple_quad_l4_l3x,   simple_quad_l4_l4x,   simple_quad_l4_l5x,
      simple_quad_l4_l6x,   simple_quad_l4_l7x,   simple_quad_l4_l8x,   simple_quad_l4_l9x,   simple_quad_l4_l10x,
      simple_quad_l5_l0x_0, simple_quad_l5_l0x_1, simple_quad_l5_l1x_0, simple_quad_l5_l1x_1, simple_quad_l5_l2x,
      simple_quad_l5_l3x,   simple_quad_l5_l4x,   simple_quad_l5_l5x,   simple_quad_l5_l6x,   simple_quad_l5_l7x,
      simple_quad_l5_l8x,   simple_quad_l5_l9x,   simple_quad_l5_l10x,   simple_quad_l6_l0x,   simple_quad_l6_l1x,
      simple_quad_l6_l2x,   simple_quad_l6_l3x,   simple_quad_l6_l4x,   simple_quad_l6_l5x,   simple_quad_l6_l6x,
      simple_quad_l6_l7x,   simple_quad_l6_l8x,   simple_quad_l6_l9x,   simple_quad_l6_l10x,   simple_quad_l7_l0x_0,
      simple_quad_l7_l0x_1, simple_quad_l7_l1x_0, simple_quad_l7_l1x_1, simple_quad_l7_l2x,   simple_quad_l7_l3x,
      simple_quad_l7_l4x,   simple_quad_l7_l5x,   simple_quad_l7_l6x,   simple_quad_l7_l7x,   simple_quad_l7_l8x,
      simple_quad_l7_l9x,   simple_quad_l7_l10x,   simple_quad_l8_l0x,   simple_quad_l8_l1x,   simple_quad_l8_l2x,
      simple_quad_l8_l3x,   simple_quad_l8_l4x,   simple_quad_l8_l5x,   simple_quad_l8_l6x,   simple_quad_l8_l7x,
      simple_quad_l8_l8x,   simple_quad_l8_l9x,   simple_quad_l8_l10x,   simple_quad_l9_l0x_0, simple_quad_l9_l0x_1,
      simple_quad_l9_l1x_0, simple_quad_l9_l1x_1, simple_quad_l9_l2x,   simple_quad_l9_l3x,   simple_quad_l9_l4x,
      simple_quad_l9_l5x,   simple_quad_l9_l6x,   simple_quad_l9_l7x,   simple_quad_l9_l8x,   simple_quad_l9_l9x,
      simple_quad_l9_l10x,   simple_quad_l10_l0x,   simple_quad_l10_l1x,   simple_quad_l10_l2x,   simple_quad_l10_l3x,
      simple_quad_l10_l4x,   simple_quad_l10_l5x,   simple_quad_l10_l6x,   simple_quad_l10_l7x,   simple_quad_l10_l8x,
      simple_quad_l10_l9x,   simple_quad_l10_l10x,
    };
    static Shapeset::shape_fn_t simple_quad_fn_dy[] =
    {
      simple_quad_l0_l0y,   simple_quad_l0_l1y,   simple_quad_l0_l2y,   simple_quad_l0_l3y_0, simple_quad_l0_l3y_1,
      simple_quad_l0_l4y,   simple_quad_l0_l5y_0, simple_quad_l0_l5y_1, simple_quad_l0_l6y,   simple_quad_l0_l7y_0,
      simple_quad_l0_l7y_1, simple_quad_l0_l8y,   simple_quad_l0_l9y_0, simple_quad_l0_l9y_1, simple_quad_l0_l10y,
      simple_quad_l1_l0y,   simple_quad_l1_l1y,   simple_quad_l1_l2y,   simple_quad_l1_l3y_0, simple_quad_l1_l3y_1,
      simple_quad_l1_l4y,   simple_quad_l1_l5y_0, simple_quad_l1_l5y_1, simple_quad_l1_l6y,   simple_quad_l1_l7y_0,
      simple_quad_l1_l7y_1, simple_quad_l1_l8y,   simple_quad_l1_l9y_0, simple_quad_l1_l9y_1, simple_quad_l1_l10y,
      simple_quad_l2_l0y,   simple_quad_l2_l1y,   simple_quad_l2_l2y,   simple_quad_l2_l3y,   simple_quad_l2_l4y,
      simple_quad_l2_l5y,   simple_quad_l2_l6y,   simple_quad_l2_l7y,   simple_quad_l2_l8y,   simple_quad_l2_l9y,
      simple_quad_l2_l10y,   simple_quad_l3_l0y_0, simple_quad_l3_l0y_1, simple_quad_l3_l1y_0, simple_quad_l3_l1y_1,
      simple_quad_l3_l2y,   simple_quad_l3_l3y,   simple_quad_l3_l4y,   simple_quad_l3_l5y,   simple_quad_l3_l6y,
      simple_quad_l3_l7y,   simple_quad_l3_l8y,   simple_quad_l3_l9y,   simple_quad_l3_l10y,   simple_quad_l4_l0y,
      simple_quad_l4_l1y,   simple_quad_l4_l2y,   simple_quad_l4_l3y,   simple_quad_l4_l4y,   simple_quad_l4_l5y,
      simple_quad_l4_l6y,   simple_quad_l4_l7y,   simple_quad_l4_l8y,   simple_quad_l4_l9y,   simple_quad_l4_l10y,
      simple_quad_l5_l0y_0, simple_quad_l5_l0y_1, simple_quad_l5_l1y_0, simple_quad_l5_l1y_1, simple_quad_l5_l2y,
      simple_quad_l5_l3y,   simple_quad_l5_l4y,   simple_quad_l5_l5y,   simple_quad_l5_l6y,   simple_quad_l5_l7y,
      simple_quad_l5_l8y,   simple_quad_l5_l9y,   simple_quad_l5_l10y,   simple_quad_l6_l0y,   simple_quad_l6_l1y,
      simple_quad_l6_l2y,   simple_quad_l6_l3y,   simple_quad_l6_l4y,   simple_quad_l6_l5y,   simple_quad_l6_l6y,
      simple_quad_l6_l7y,   simple_quad_l6_l8y,   simple_quad_l6_l9y,   simple_quad_l6_l10y,   simple_quad_l7_l0y_0,
      simple_quad_l7_l0y_1, simple_quad_l7_l1y_0, simple_quad_l7_l1y_1, simple_quad_l7_l2y,   simple_quad_l7_l3y,
      simple_quad_l7_l4y,   simple_quad_l7_l5y,   simple_quad_l7_l6y,   simple_quad_l7_l7y,   simple_quad_l7_l8y,
      simple_quad_l7_l9y,   simple_quad_l7_l10y,   simple_quad_l8_l0y,   simple_quad_l8_l1y,   simple_quad_l8_l2y,
      simple_quad_l8_l3y,   simple_quad_l8_l4y,   simple_quad_l8_l5y,   simple_quad_l8_l6y,   simple_quad_l8_l7y,
      simple_quad_l8_l8y,   simple_quad_l8_l9y,   simple_quad_l8_l10y,   simple_quad_l9_l0y_0, simple_quad_l9_l0y_1,
      simple_quad_l9_l1y_0, simple_quad_l9_l1y_1, simple_quad_l9_l2y,   simple_quad_l9_l3y,   simple_quad_l9_l4y,
      simple_quad_l9_l5y,   simple_quad_l9_l6y,   simple_quad_l9_l7y,   simple_quad_l9_l8y,   simple_quad_l9_l9y,
      simple_quad_l9_l10y,   simple_quad_l10_l0y,   simple_quad_l10_l1y,   simple_quad_l10_l2y,   simple_quad_l10_l3y,
      simple_quad_l10_l4y,   simple_quad_l10_l5y,   simple_quad_l10_l6y,   simple_quad_l10_l7y,   simple_quad_l10_l8y,
      simple_quad_l10_l9y,   simple_quad_l10_l10y,
    };
    static Shapeset::shape_fn_t simple_quad_fn_dxx[] =
    {
      simple_quad_l0_l0xx,   simple_quad_l0_l1xx,   simple_quad_l0_l2xx,   simple_quad_l0_l3xx_0, simple_quad_l0_l3xx_1,
      simple_quad_l0_l4xx,   simple_quad_l0_l5xx_0, simple_quad_l0_l5xx_1, simple_quad_l0_l6xx,   simple_quad_l0_l7xx_0,
      simple_quad_l0_l7xx_1, simple_quad_l0_l8xx,   simple_quad_l0_l9xx_0, simple_quad_l0_l9xx_1, simple_quad_l0_l10xx,
      simple_quad_l1_l0xx,   simple_quad_l1_l1xx,   simple_quad_l1_l2xx,   simple_quad_l1_l3xx_0, simple_quad_l1_l3xx_1,
      simple_quad_l1_l4xx,   simple_quad_l1_l5xx_0, simple_quad_l1_l5xx_1, simple_quad_l1_l6xx,   simple_quad_l1_l7xx_0,
      simple_quad_l1_l7xx_1, simple_quad_l1_l8xx,   simple_quad_l1_l9xx_0, simple_quad_l1_l9xx_1, simple_quad_l1_l10xx,
      simple_quad_l2_l0xx,   simple_quad_l2_l1xx,   simple_quad_l2_l2xx,   simple_quad_l2_l3xx,   simple_quad_l2_l4xx,
      simple_quad_l2_l5xx,   simple_quad_l2_l6xx,   simple_quad_l2_l7xx,   simple_quad_l2_l8xx,   simple_quad_l2_l9xx,
      simple_quad_l2_l10xx,   simple_quad_l3_l0xx_0, simple_quad_l3_l0xx_1, simple_quad_l3_l1xx_0, simple_quad_l3_l1xx_1,
      simple_quad_l3_l2xx,   simple_quad_l3_l3xx,   simple_quad_l3_l4xx,   simple_quad_l3_l5xx,   simple_quad_l3_l6xx,
      simple_quad_l3_l7xx,   simple_quad_l3_l8xx,   simple_quad_l3_l9xx,   simple_quad_l3_l10xx,   simple_quad_l4_l0xx,
      simple_quad_l4_l1xx,   simple_quad_l4_l2xx,   simple_quad_l4_l3xx,   simple_quad_l4_l4xx,   simple_quad_l4_l5xx,
      simple_quad_l4_l6xx,   simple_quad_l4_l7xx,   simple_quad_l4_l8xx,   simple_quad_l4_l9xx,   simple_quad_l4_l10xx,
      simple_quad_l5_l0xx_0, simple_quad_l5_l0xx_1, simple_quad_l5_l1xx_0, simple_quad_l5_l1xx_1, simple_quad_l5_l2xx,
      simple_quad_l5_l3xx,   simple_quad_l5_l4xx,   simple_quad_l5_l5xx,   simple_quad_l5_l6xx,   simple_quad_l5_l7xx,
      simple_quad_l5_l8xx,   simple_quad_l5_l9xx,   simple_quad_l5_l10xx,   simple_quad_l6_l0xx,   simple_quad_l6_l1xx,
      simple_quad_l6_l2xx,   simple_quad_l6_l3xx,   simple_quad_l6_l4xx,   simple_quad_l6_l5xx,   simple_quad_l6_l6xx,
      simple_quad_l6_l7xx,   simple_quad_l6_l8xx,   simple_quad_l6_l9xx,   simple_quad_l6_l10xx,   simple_quad_l7_l0xx_0,
      simple_quad_l7_l0xx_1, simple_quad_l7_l1xx_0, simple_quad_l7_l1xx_1, simple_quad_l7_l2xx,   simple_quad_l7_l3xx,
      simple_quad_l7_l4xx,   simple_quad_l7_l5xx,   simple_quad_l7_l6xx,   simple_quad_l7_l7xx,   simple_quad_l7_l8xx,
      simple_quad_l7_l9xx,   simple_quad_l7_l10xx,   simple_quad_l8_l0xx,   simple_quad_l8_l1xx,   simple_quad_l8_l2xx,
      simple_quad_l8_l3xx,   simple_quad_l8_l4xx,   simple_quad_l8_l5xx,   simple_quad_l8_l6xx,   simple_quad_l8_l7xx,
      simple_quad_l8_l8xx,   simple_quad_l8_l9xx,   simple_quad_l8_l10xx,   simple_quad_l9_l0xx_0, simple_quad_l9_l0xx_1,
      simple_quad_l9_l1xx_0, simple_quad_l9_l1xx_1, simple_quad_l9_l2xx,   simple_quad_l9_l3xx,   simple_quad_l9_l4xx,
      simple_quad_l9_l5xx,   simple_quad_l9_l6xx,   simple_quad_l9_l7xx,   simple_quad_l9_l8xx,   simple_quad_l9_l9xx,
      simple_quad_l9_l10xx,   simple_quad_l10_l0xx,   simple_quad_l10_l1xx,   simple_quad_l10_l2xx,   simple_quad_l10_l3xx,
      simple_quad_l10_l4xx,   simple_quad_l10_l5xx,   simple_quad_l10_l6xx,   simple_quad_l10_l7xx,   simple_quad_l10_l8xx,
      simple_quad_l10_l9xx,   simple_quad_l10_l10xx,
    };
    static Shapeset::shape_fn_t simple_quad_fn_dxy[] =
    {
      simple_quad_l0_l0xy,   simple_quad_l0_l1xy,   simple_quad_l0_l2xy,   simple_quad_l0_l3xy_0, simple_quad_l0_l3xy_1,
      simple_quad_l0_l4xy,   simple_quad_l0_l5xy_0, simple_quad_l0_l5xy_1, simple_quad_l0_l6xy,   simple_quad_l0_l7xy_0,
      simple_quad_l0_l7xy_1, simple_quad_l0_l8xy,   simple_quad_l0_l9xy_0, simple_quad_l0_l9xy_1, simple_quad_l0_l10xy,
      simple_quad_l1_l0xy,   simple_quad_l1_l1xy,   simple_quad_l1_l2xy,   simple_quad_l1_l3xy_0, simple_quad_l1_l3xy_1,
      simple_quad_l1_l4xy,   simple_quad_l1_l5xy_0, simple_quad_l1_l5xy_1, simple_quad_l1_l6xy,   simple_quad_l1_l7xy_0,
      simple_quad_l1_l7xy_1, simple_quad_l1_l8xy,   simple_quad_l1_l9xy_0, simple_quad_l1_l9xy_1, simple_quad_l1_l10xy,
      simple_quad_l2_l0xy,   simple_quad_l2_l1xy,   simple_quad_l2_l2xy,   simple_quad_l2_l3xy,   simple_quad_l2_l4xy,
      simple_quad_l2_l5xy,   simple_quad_l2_l6xy,   simple_quad_l2_l7xy,   simple_quad_l2_l8xy,   simple_quad_l2_l9xy,
      simple_quad_l2_l10xy,   simple_quad_l3_l0xy_0, simple_quad_l3_l0xy_1, simple_quad_l3_l1xy_0, simple_quad_l3_l1xy_1,
      simple_quad_l3_l2xy,   simple_quad_l3_l3xy,   simple_quad_l3_l4xy,   simple_quad_l3_l5xy,   simple_quad_l3_l6xy,
      simple_quad_l3_l7xy,   simple_quad_l3_l8xy,   simple_quad_l3_l9xy,   simple_quad_l3_l10xy,   simple_quad_l4_l0xy,
      simple_quad_l4_l1xy,   simple_quad_l4_l2xy,   simple_quad_l4_l3xy,   simple_quad_l4_l4xy,   simple_quad_l4_l5xy,
      simple_quad_l4_l6xy,   simple_quad_l4_l7xy,   simple_quad_l4_l8xy,   simple_quad_l4_l9xy,   simple_quad_l4_l10xy,
      simple_quad_l5_l0xy_0, simple_quad_l5_l0xy_1, simple_quad_l5_l1xy_0, simple_quad_l5_l1xy_1, simple_quad_l5_l2xy,
      simple_quad_l5_l3xy,   simple_quad_l5_l4xy,   simple_quad_l5_l5xy,   simple_quad_l5_l6xy,   simple_quad_l5_l7xy,
      simple_quad_l5_l8xy,   simple_quad_l5_l9xy,   simple_quad_l5_l10xy,   simple_quad_l6_l0xy,   simple_quad_l6_l1xy,
      simple_quad_l6_l2xy,   simple_quad_l6_l3xy,   simple_quad_l6_l4xy,   simple_quad_l6_l5xy,   simple_quad_l6_l6xy,
      simple_quad_l6_l7xy,   simple_quad_l6_l8xy,   simple_quad_l6_l9xy,   simple_quad_l6_l10xy,   simple_quad_l7_l0xy_0,
      simple_quad_l7_l0xy_1, simple_quad_l7_l1xy_0, simple_quad_l7_l1xy_1, simple_quad_l7_l2xy,   simple_quad_l7_l3xy,
      simple_quad_l7_l4xy,   simple_quad_l7_l5xy,   simple_quad_l7_l6xy,   simple_quad_l7_l7xy,   simple_quad_l7_l8xy,
      simple_quad_l7_l9xy,   simple_quad_l7_l10xy,   simple_quad_l8_l0xy,   simple_quad_l8_l1xy,   simple_quad_l8_l2xy,
      simple_quad_l8_l3xy,   simple_quad_l8_l4xy,   simple_quad_l8_l5xy,   simple_quad_l8_l6xy,   simple_quad_l8_l7xy,
      simple_quad_l8_l8xy,   simple_quad_l8_l9xy,   simple_quad_l8_l10xy,   simple_quad_l9_l0xy_0, simple_quad_l9_l0xy_1,
      simple_quad_l9_l1xy_0, simple_quad_l9_l1xy_1, simple_quad_l9_l2xy,   simple_quad_l9_l3xy,   simple_quad_l9_l4xy,
      simple_quad_l9_l5xy,   simple_quad_l9_l6xy,   simple_quad_l9_l7xy,   simple_quad_l9_l8xy,   simple_quad_l9_l9xy,
      simple_quad_l9_l10xy,   simple_quad_l10_l0xy,   simple_quad_l10_l1xy,   simple_quad_l10_l2xy,   simple_quad_l10_l3xy,
      simple_quad_l10_l4xy,   simple_quad_l10_l5xy,   simple_quad_l10_l6xy,   simple_quad_l10_l7xy,   simple_quad_l10_l8xy,
      simple_quad_l10_l9xy,   simple_quad_l10_l10xy,
    };
    static Shapeset::shape_fn_t simple_quad_fn_dyy[] =
    {
      simple_quad_l0_l0yy,   simple_quad_l0_l1yy,   simple_quad_l0_l2yy,   simple_quad_l0_l3yy_0, simple_quad_l0_l3yy_1,
      simple_quad_l0_l4yy,   simple_quad_l0_l5yy_0, simple_quad_l0_l5yy_1, simple_quad_l0_l6yy,   simple_quad_l0_l7yy_0,
      simple_quad_l0_l7yy_1, simple_quad_l0_l8yy,   simple_quad_l0_l9yy_0, simple_quad_l0_l9yy_1, simple_quad_l0_l10yy,
      simple_quad_l1_l0yy,   simple_quad_l1_l1yy,   simple_quad_l1_l2yy,   simple_quad_l1_l3yy_0, simple_quad_l1_l3yy_1,
      simple_quad_l1_l4yy,   simple_quad_l1_l5yy_0, simple_quad_l1_l5yy_1, simple_quad_l1_l6yy,   simple_quad_l1_l7yy_0,
      simple_quad_l1_l7yy_1, simple_quad_l1_l8yy,   simple_quad_l1_l9yy_0, simple_quad_l1_l9yy_1, simple_quad_l1_l10yy,
      simple_quad_l2_l0yy,   simple_quad_l2_l1yy,   simple_quad_l2_l2yy,   simple_quad_l2_l3yy,   simple_quad_l2_l4yy,
      simple_quad_l2_l5yy,   simple_quad_l2_l6yy,   simple_quad_l2_l7yy,   simple_quad_l2_l8yy,   simple_quad_l2_l9yy,
      simple_quad_l2_l10yy,   simple_quad_l3_l0yy_0, simple_quad_l3_l0yy_1, simple_quad_l3_l1yy_0, simple_quad_l3_l1yy_1,
      simple_quad_l3_l2yy,   simple_quad_l3_l3yy,   simple_quad_l3_l4yy,   simple_quad_l3_l5yy,   simple_quad_l3_l6yy,
      simple_quad_l3_l7yy,   simple_quad_l3_l8yy,   simple_quad_l3_l9yy,   simple_quad_l3_l10yy,   simple_quad_l4_l0yy,
      simple_quad_l4_l1yy,   simple_quad_l4_l2yy,   simple_quad_l4_l3yy,   simple_quad_l4_l4yy,   simple_quad_l4_l5yy,
      simple_quad_l4_l6yy,   simple_quad_l4_l7yy,   simple_quad_l4_l8yy,   simple_quad_l4_l9yy,   simple_quad_l4_l10yy,
      simple_quad_l5_l0yy_0, simple_quad_l5_l0yy_1, simple_quad_l5_l1yy_0, simple_quad_l5_l1yy_1, simple_quad_l5_l2yy,
      simple_quad_l5_l3yy,   simple_quad_l5_l4yy,   simple_quad_l5_l5yy,   simple_quad_l5_l6yy,   simple_quad_l5_l7yy,
      simple_quad_l5_l8yy,   simple_quad_l5_l9yy,   simple_quad_l5_l10yy,   simple_quad_l6_l0yy,   simple_quad_l6_l1yy,
      simple_quad_l6_l2yy,   simple_quad_l6_l3yy,   simple_quad_l6_l4yy,   simple_quad_l6_l5yy,   simple_quad_l6_l6yy,
      simple_quad_l6_l7yy,   simple_quad_l6_l8yy,   simple_quad_l6_l9yy,   simple_quad_l6_l10yy,   simple_quad_l7_l0yy_0,
      simple_quad_l7_l0yy_1, simple_quad_l7_l1yy_0, simple_quad_l7_l1yy_1, simple_quad_l7_l2yy,   simple_quad_l7_l3yy,
      simple_quad_l7_l4yy,   simple_quad_l7_l5yy,   simple_quad_l7_l6yy,   simple_quad_l7_l7yy,   simple_quad_l7_l8yy,
      simple_quad_l7_l9yy,   simple_quad_l7_l10yy,   simple_quad_l8_l0yy,   simple_quad_l8_l1yy,   simple_quad_l8_l2yy,
      simple_quad_l8_l3yy,   simple_quad_l8_l4yy,   simple_quad_l8_l5yy,   simple_quad_l8_l6yy,   simple_quad_l8_l7yy,
      simple_quad_l8_l8yy,   simple_quad_l8_l9yy,   simple_quad_l8_l10yy,   simple_quad_l9_l0yy_0, simple_quad_l9_l0yy_1,
      simple_quad_l9_l1yy_0, simple_quad_l9_l1yy_1, simple_quad_l9_l2yy,   simple_quad_l9_l3yy,   simple_quad_l9_l4yy,
      simple_quad_l9_l5yy,   simple_quad_l9_l6yy,   simple_quad_l9_l7yy,   simple_quad_l9_l8yy,   simple_quad_l9_l9yy,
      simple_quad_l9_l10yy,   simple_quad_l10_l0yy,   simple_quad_l10_l1yy,   simple_quad_l10_l2yy,   simple_quad_l10_l3yy,
      simple_quad_l10_l4yy,   simple_quad_l10_l5yy,   simple_quad_l10_l6yy,   simple_quad_l10_l7yy,   simple_quad_l10_l8yy,
      simple_quad_l10_l9yy,   simple_quad_l10_l10yy,
    };
    Shapeset::shape_fn_t* simple_quad_shape_fn_table[1]     = { simple_quad_fn };
    Shapeset::shape_fn_t* simple_quad_shape_fn_table_dx[1]  = { simple_quad_fn_dx };
    Shapeset::shape_fn_t* simple_quad_shape_fn_table_dy[1]  = { simple_quad_fn_dy };
    Shapeset::shape_fn_t* simple_quad_shape_fn_table_dxx[1] = { simple_quad_fn_dxx };
    Shapeset::shape_fn_t* simple_quad_shape_fn_table_dxy[1] = { simple_quad_fn_dxy };
    Shapeset::shape_fn_t* simple_quad_shape_fn_table_dyy[1] = { simple_quad_fn_dyy };

    static int qb_2_2[] = { 32, };
    static int qb_2_3[] = { 32, 33, };
    static int qb_2_4[] = { 32, 33, 34, };
    static int qb_2_5[] = { 32, 33, 34, 35, };
    static int qb_2_6[] = { 32, 33, 34, 35, 36, };
    static int qb_2_7[] = { 32, 33, 34, 35, 36, 37, };
    static int qb_2_8[] = { 32, 33, 34, 35, 36, 37, 38, };
    static int qb_2_9[] = { 32, 33, 34, 35, 36, 37, 38, 39, };
    static int qb_2_10[] = { 32, 33, 34, 35, 36, 37, 38, 39, 40, };
    static int qb_3_2[] = { 32, 45, };
    static int qb_3_3[] = { 32, 33, 45, 46, };
    static int qb_3_4[] = { 32, 33, 34, 45, 46, 47, };
    static int qb_3_5[] = { 32, 33, 34, 35, 45, 46, 47, 48, };
    static int qb_3_6[] = { 32, 33, 34, 35, 36, 45, 46, 47, 48, 49, };
    static int qb_3_7[] = { 32, 33, 34, 35, 36, 37, 45, 46, 47, 48, 49, 50, };
    static int qb_3_8[] = { 32, 33, 34, 35, 36, 37, 38, 45, 46, 47, 48, 49, 50, 51, };
    static int qb_3_9[] = { 32, 33, 34, 35, 36, 37, 38, 39, 45, 46, 47, 48, 49, 50, 51, 52, };
    static int qb_3_10[] = { 32, 33, 34, 35, 36, 37, 38, 39, 40, 45, 46, 47, 48, 49, 50, 51, 52, 53, };
    static int qb_4_2[] = { 32, 45, 56, };
    static int qb_4_3[] = { 32, 33, 45, 46, 56, 57, };
    static int qb_4_4[] = { 32, 33, 34, 45, 46, 47, 56, 57, 58, };
    static int qb_4_5[] = { 32, 33, 34, 35, 45, 46, 47, 48, 56, 57, 58, 59, };
    static int qb_4_6[] = { 32, 33, 34, 35, 36, 45, 46, 47, 48, 49, 56, 57, 58, 59, 60, };
    static int qb_4_7[] = { 32, 33, 34, 35, 36, 37, 45, 46, 47, 48, 49, 50, 56, 57, 58, 59, 60, 61, };
    static int qb_4_8[] = { 32, 33, 34, 35, 36, 37, 38, 45, 46, 47, 48, 49, 50, 51, 56, 57, 58, 59, 60, 61, 62, };
    static int qb_4_9[] = { 32, 33, 34, 35, 36, 37, 38, 39, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 62, 63, };
    static int qb_4_10[] = { 32, 33, 34, 35, 36, 37, 38, 39, 40, 45, 46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 60, 61, 62, 63, 64, };
    static int qb_5_2[] = { 32, 45, 56, 69, };
    static int qb_5_3[] = { 32, 33, 45, 46, 56, 57, 69, 70, };
    static int qb_5_4[] = { 32, 33, 34, 45, 46, 47, 56, 57, 58, 69, 70, 71, };
    static int qb_5_5[] = { 32, 33, 34, 35, 45, 46, 47, 48, 56, 57, 58, 59, 69, 70, 71, 72, };
    static int qb_5_6[] = { 32, 33, 34, 35, 36, 45, 46, 47, 48, 49, 56, 57, 58, 59, 60, 69, 70, 71, 72, 73, };
    static int qb_5_7[] = { 32, 33, 34, 35, 36, 37, 45, 46, 47, 48, 49, 50, 56, 57, 58, 59, 60, 61, 69, 70, 71, 72, 73, 74, };
    static int qb_5_8[] = { 32, 33, 34, 35, 36, 37, 38, 45, 46, 47, 48, 49, 50, 51, 56, 57, 58, 59, 60, 61, 62, 69, 70, 71, 72, 73, 74, 75, };
    static int qb_5_9[] = { 32, 33, 34, 35, 36, 37, 38, 39, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 62, 63, 69, 70, 71, 72, 73, 74, 75, 76, };
    static int qb_5_10[] = { 32, 33, 34, 35, 36, 37, 38, 39, 40, 45, 46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 60, 61, 62, 63, 64, 69, 70, 71, 72, 73, 74, 75, 76, 77, };
    static int qb_6_2[] = { 32, 45, 56, 69, 80, };
    static int qb_6_3[] = { 32, 33, 45, 46, 56, 57, 69, 70, 80, 81, };
    static int qb_6_4[] = { 32, 33, 34, 45, 46, 47, 56, 57, 58, 69, 70, 71, 80, 81, 82, };
    static int qb_6_5[] = { 32, 33, 34, 35, 45, 46, 47, 48, 56, 57, 58, 59, 69, 70, 71, 72, 80, 81, 82, 83, };
    static int qb_6_6[] = { 32, 33, 34, 35, 36, 45, 46, 47, 48, 49, 56, 57, 58, 59, 60, 69, 70, 71, 72, 73, 80, 81, 82, 83, 84, };
    static int qb_6_7[] = { 32, 33, 34, 35, 36, 37, 45, 46, 47, 48, 49, 50, 56, 57, 58, 59, 60, 61, 69, 70, 71, 72, 73, 74, 80, 81, 82, 83, 84, 85, };
    static int qb_6_8[] = { 32, 33, 34, 35, 36, 37, 38, 45, 46, 47, 48, 49, 50, 51, 56, 57, 58, 59, 60, 61, 62, 69, 70, 71, 72, 73, 74, 75, 80, 81, 82, 83, 84, 85, 86, };
    static int qb_6_9[] = { 32, 33, 34, 35, 36, 37, 38, 39, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 62, 63, 69, 70, 71, 72, 73, 74, 75, 76, 80, 81, 82, 83, 84, 85, 86, 87, };
    static int qb_6_10[] = { 32, 33, 34, 35, 36, 37, 38, 39, 40, 45, 46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 60, 61, 62, 63, 64, 69, 70, 71, 72, 73, 74, 75, 76, 77, 80, 81, 82, 83, 84, 85, 86, 87, 88, };
    static int qb_7_2[] = { 32, 45, 56, 69, 80, 93, };
    static int qb_7_3[] = { 32, 33, 45, 46, 56, 57, 69, 70, 80, 81, 93, 94, };
    static int qb_7_4[] = { 32, 33, 34, 45, 46, 47, 56, 57, 58, 69, 70, 71, 80, 81, 82, 93, 94, 95, };
    static int qb_7_5[] = { 32, 33, 34, 35, 45, 46, 47, 48, 56, 57, 58, 59, 69, 70, 71, 72, 80, 81, 82, 83, 93, 94, 95, 96, };
    static int qb_7_6[] = { 32, 33, 34, 35, 36, 45, 46, 47, 48, 49, 56, 57, 58, 59, 60, 69, 70, 71, 72, 73, 80, 81, 82, 83, 84, 93, 94, 95, 96, 97, };
    static int qb_7_7[] = { 32, 33, 34, 35, 36, 37, 45, 46, 47, 48, 49, 50, 56, 57, 58, 59, 60, 61, 69, 70, 71, 72, 73, 74, 80, 81, 82, 83, 84, 85, 93, 94, 95, 96, 97, 98, };
    static int qb_7_8[] = { 32, 33, 34, 35, 36, 37, 38, 45, 46, 47, 48, 49, 50, 51, 56, 57, 58, 59, 60, 61, 62, 69, 70, 71, 72, 73, 74, 75, 80, 81, 82, 83, 84, 85, 86, 93, 94, 95, 96, 97, 98, 99, };
    static int qb_7_9[] = { 32, 33, 34, 35, 36, 37, 38, 39, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 62, 63, 69, 70, 71, 72, 73, 74, 75, 76, 80, 81, 82, 83, 84, 85, 86, 87, 93, 94, 95, 96, 97, 98, 99, 100, };
    static int qb_7_10[] = { 32, 33, 34, 35, 36, 37, 38, 39, 40, 45, 46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 60, 61, 62, 63, 64, 69, 70, 71, 72, 73, 74, 75, 76, 77, 80, 81, 82, 83, 84, 85, 86, 87, 88, 93, 94, 95, 96, 97, 98, 99, 100, 101, };
    static int qb_8_2[] = { 32, 45, 56, 69, 80, 93, 104, };
    static int qb_8_3[] = { 32, 33, 45, 46, 56, 57, 69, 70, 80, 81, 93, 94, 104, 105, };
    static int qb_8_4[] = { 32, 33, 34, 45, 46, 47, 56, 57, 58, 69, 70, 71, 80, 81, 82, 93, 94, 95, 104, 105, 106, };
    static int qb_8_5[] = { 32, 33, 34, 35, 45, 46, 47, 48, 56, 57, 58, 59, 69, 70, 71, 72, 80, 81, 82, 83, 93, 94, 95, 96, 104, 105, 106, 107, };
    static int qb_8_6[] = { 32, 33, 34, 35, 36, 45, 46, 47, 48, 49, 56, 57, 58, 59, 60, 69, 70, 71, 72, 73, 80, 81, 82, 83, 84, 93, 94, 95, 96, 97, 104, 105, 106, 107, 108, };
    static int qb_8_7[] = { 32, 33, 34, 35, 36, 37, 45, 46, 47, 48, 49, 50, 56, 57, 58, 59, 60, 61, 69, 70, 71, 72, 73, 74, 80, 81, 82, 83, 84, 85, 93, 94, 95, 96, 97, 98, 104, 105, 106, 107, 108, 109, };
    static int qb_8_8[] = { 32, 33, 34, 35, 36, 37, 38, 45, 46, 47, 48, 49, 50, 51, 56, 57, 58, 59, 60, 61, 62, 69, 70, 71, 72, 73, 74, 75, 80, 81, 82, 83, 84, 85, 86, 93, 94, 95, 96, 97, 98, 99, 104, 105, 106, 107, 108, 109, 110, };
    static int qb_8_9[] = { 32, 33, 34, 35, 36, 37, 38, 39, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 62, 63, 69, 70, 71, 72, 73, 74, 75, 76, 80, 81, 82, 83, 84, 85, 86, 87, 93, 94, 95, 96, 97, 98, 99, 100, 104, 105, 106, 107, 108, 109, 110, 111, };
    static int qb_8_10[] = { 32, 33, 34, 35, 36, 37, 38, 39, 40, 45, 46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 60, 61, 62, 63, 64, 69, 70, 71, 72, 73, 74, 75, 76, 77, 80, 81, 82, 83, 84, 85, 86, 87, 88, 93, 94, 95, 96, 97, 98, 99, 100, 101, 104, 105, 106, 107, 108, 109, 110, 111, 112, };
    static int qb_9_2[] = { 32, 45, 56, 69, 80, 93, 104, 117, };
    static int qb_9_3[] = { 32, 33, 45, 46, 56, 57, 69, 70, 80, 81, 93, 94, 104, 105, 117, 118, };
    static int qb_9_4[] = { 32, 33, 34, 45, 46, 47, 56, 57, 58, 69, 70, 71, 80, 81, 82, 93, 94, 95, 104, 105, 106, 117, 118, 119, };
    static int qb_9_5[] = { 32, 33, 34, 35, 45, 46, 47, 48, 56, 57, 58, 59, 69, 70, 71, 72, 80, 81, 82, 83, 93, 94, 95, 96, 104, 105, 106, 107, 117, 118, 119, 120, };
    static int qb_9_6[] = { 32, 33, 34, 35, 36, 45, 46, 47, 48, 49, 56, 57, 58, 59, 60, 69, 70, 71, 72, 73, 80, 81, 82, 83, 84, 93, 94, 95, 96, 97, 104, 105, 106, 107, 108, 117, 118, 119, 120, 121, };
    static int qb_9_7[] = { 32, 33, 34, 35, 36, 37, 45, 46, 47, 48, 49, 50, 56, 57, 58, 59, 60, 61, 69, 70, 71, 72, 73, 74, 80, 81, 82, 83, 84, 85, 93, 94, 95, 96, 97, 98, 104, 105, 106, 107, 108, 109, 117, 118, 119, 120, 121, 122, };
    static int qb_9_8[] = { 32, 33, 34, 35, 36, 37, 38, 45, 46, 47, 48, 49, 50, 51, 56, 57, 58, 59, 60, 61, 62, 69, 70, 71, 72, 73, 74, 75, 80, 81, 82, 83, 84, 85, 86, 93, 94, 95, 96, 97, 98, 99, 104, 105, 106, 107, 108, 109, 110, 117, 118, 119, 120, 121, 122, 123, };
    static int qb_9_9[] = { 32, 33, 34, 35, 36, 37, 38, 39, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 62, 63, 69, 70, 71, 72, 73, 74, 75, 76, 80, 81, 82, 83, 84, 85, 86, 87, 93, 94, 95, 96, 97, 98, 99, 100, 104, 105, 106, 107, 108, 109, 110, 111, 117, 118, 119, 120, 121, 122, 123, 124, };
    static int qb_9_10[] = { 32, 33, 34, 35, 36, 37, 38, 39, 40, 45, 46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 60, 61, 62, 63, 64, 69, 70, 71, 72, 73, 74, 75, 76, 77, 80, 81, 82, 83, 84, 85, 86, 87, 88, 93, 94, 95, 96, 97, 98, 99, 100, 101, 104, 105, 106, 107, 108, 109, 110, 111, 112, 117, 118, 119, 120, 121, 122, 123, 124, 125, };
    static int qb_10_2[] = { 32, 45, 56, 69, 80, 93, 104, 117, 128, };
    static int qb_10_3[] = { 32, 33, 45, 46, 56, 57, 69, 70, 80, 81, 93, 94, 104, 105, 117, 118, 128, 129, };
    static int qb_10_4[] = { 32, 33, 34, 45, 46, 47, 56, 57, 58, 69, 70, 71, 80, 81, 82, 93, 94, 95, 104, 105, 106, 117, 118, 119, 128, 129, 130, };
    static int qb_10_5[] = { 32, 33, 34, 35, 45, 46, 47, 48, 56, 57, 58, 59, 69, 70, 71, 72, 80, 81, 82, 83, 93, 94, 95, 96, 104, 105, 106, 107, 117, 118, 119, 120, 128, 129, 130, 131, };
    static int qb_10_6[] = { 32, 33, 34, 35, 36, 45, 46, 47, 48, 49, 56, 57, 58, 59, 60, 69, 70, 71, 72, 73, 80, 81, 82, 83, 84, 93, 94, 95, 96, 97, 104, 105, 106, 107, 108, 117, 118, 119, 120, 121, 128, 129, 130, 131, 132, };
    static int qb_10_7[] = { 32, 33, 34, 35, 36, 37, 45, 46, 47, 48, 49, 50, 56, 57, 58, 59, 60, 61, 69, 70, 71, 72, 73, 74, 80, 81, 82, 83, 84, 85, 93, 94, 95, 96, 97, 98, 104, 105, 106, 107, 108, 109, 117, 118, 119, 120, 121, 122, 128, 129, 130, 131, 132, 133, };
    static int qb_10_8[] = { 32, 33, 34, 35, 36, 37, 38, 45, 46, 47, 48, 49, 50, 51, 56, 57, 58, 59, 60, 61, 62, 69, 70, 71, 72, 73, 74, 75, 80, 81, 82, 83, 84, 85, 86, 93, 94, 95, 96, 97, 98, 99, 104, 105, 106, 107, 108, 109, 110, 117, 118, 119, 120, 121, 122, 123, 128, 129, 130, 131, 132, 133, 134, };
    static int qb_10_9[] = { 32, 33, 34, 35, 36, 37, 38, 39, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 62, 63, 69, 70, 71, 72, 73, 74, 75, 76, 80, 81, 82, 83, 84, 85, 86, 87, 93, 94, 95, 96, 97, 98, 99, 100, 104, 105, 106, 107, 108, 109, 110, 111, 117, 118, 119, 120, 121, 122, 123, 124, 128, 129, 130, 131, 132, 133, 134, 135, };
    static int qb_10_10[] = { 32, 33, 34, 35, 36, 37, 38, 39, 40, 45, 46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 60, 61, 62, 63, 64, 69, 70, 71, 72, 73, 74, 75, 76, 77, 80, 81, 82, 83, 84, 85, 86, 87, 88, 93, 94, 95, 96, 97, 98, 99, 100, 101, 104, 105, 106, 107, 108, 109, 110, 111, 112, 117, 118, 119, 120, 121, 122, 123, 124, 125, 128, 129, 130, 131, 132, 133, 134, 135, 136, };

#define NULL16 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL

    int* simple_quad_bubble_indices[] =
    {
      NULL, NULL, NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,     NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL, NULL, NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,     NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL, NULL, qb_2_2,  qb_2_3,  qb_2_4,  qb_2_5,  qb_2_6,  qb_2_7,  qb_2_8,  qb_2_9,  qb_2_10,  NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL, NULL, qb_3_2,  qb_3_3,  qb_3_4,  qb_3_5,  qb_3_6,  qb_3_7,  qb_3_8,  qb_3_9,  qb_3_10,  NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL, NULL, qb_4_2,  qb_4_3,  qb_4_4,  qb_4_5,  qb_4_6,  qb_4_7,  qb_4_8,  qb_4_9,  qb_4_10,  NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL, NULL, qb_5_2,  qb_5_3,  qb_5_4,  qb_5_5,  qb_5_6,  qb_5_7,  qb_5_8,  qb_5_9,  qb_5_10,  NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL, NULL, qb_6_2,  qb_6_3,  qb_6_4,  qb_6_5,  qb_6_6,  qb_6_7,  qb_6_8,  qb_6_9,  qb_6_10,  NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL, NULL, qb_7_2,  qb_7_3,  qb_7_4,  qb_7_5,  qb_7_6,  qb_7_7,  qb_7_8,  qb_7_9,  qb_7_10,  NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL, NULL, qb_8_2,  qb_8_3,  qb_8_4,  qb_8_5,  qb_8_6,  qb_8_7,  qb_8_8,  qb_8_9,  qb_8_10,  NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL, NULL, qb_9_2,  qb_9_3,  qb_9_4,  qb_9_5,  qb_9_6,  qb_9_7,  qb_9_8,  qb_9_9,  qb_9_10,  NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL, NULL, qb_10_2, qb_10_3, qb_10_4, qb_10_5, qb_10_6, qb_10_7, qb_10_8, qb_10_9, qb_10_10,  NULL, NULL, NULL, NULL, NULL, NULL16,
    };

#define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    int simple_quad_bubble_count[] =
    {
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16,
      0,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  0,  0,  0,  0,  0, zero16,
      0,  0,  2,  4,  6,  8, 10, 12, 14, 16, 18,  0,  0,  0,  0,  0, zero16,
      0,  0,  3,  6,  9, 12, 15, 18, 21, 24, 27,  0,  0,  0,  0,  0, zero16,
      0,  0,  4,  8, 12, 16, 20, 24, 28, 32, 36,  0,  0,  0,  0,  0, zero16,
      0,  0,  5, 10, 15, 20, 25, 30, 35, 40, 45,  0,  0,  0,  0,  0, zero16,
      0,  0,  6, 12, 18, 24, 30, 36, 42, 48, 54,  0,  0,  0,  0,  0, zero16,
      0,  0,  7, 14, 21, 28, 35, 42, 49, 56, 63,  0,  0,  0,  0,  0, zero16,
      0,  0,  8, 16, 24, 32, 40, 48, 56, 64, 72,  0,  0,  0,  0,  0, zero16,
      0,  0,  9, 18, 27, 36, 45, 54, 63, 72, 81,  0,  0,  0,  0,  0, zero16,
    };

    int simple_quad_vertex_indices[4] = { 0, 15, 16, 1 };

    static int simple_quad_edge_indices_0[22] =  {  0, 15, 15, 0,  30, 30, 41, 42, 54, 54, 65, 66, 78, 78, 89, 90, 102, 102, 113, 114, 126, 126 };
    static int simple_quad_edge_indices_1[22] =  { 15, 16, 16, 15, 17, 17, 18, 19, 20, 20, 21, 22, 23, 23, 24, 25, 26,   26,  27,  28,  29,  29 };
    static int simple_quad_edge_indices_2[22] =  { 16,  1,  1, 16, 31, 31, 43, 44, 55, 55, 67, 68, 79, 79, 91, 92, 103, 103, 115, 116, 127, 127 };
    static int simple_quad_edge_indices_3[22] =  {  1,  0,  0,  1,  2,  2,  3,  4,  5,  5,  6,  7,  8,  8,  9, 10,  11,  11,  12,  13,  14,  14 };

    int* simple_quad_edge_indices[4] =
    {
      simple_quad_edge_indices_0,
      simple_quad_edge_indices_1,
      simple_quad_edge_indices_2,
      simple_quad_edge_indices_3
    };

#define oo H2D_MAKE_QUAD_ORDER
#define XX(a, b) oo(a, b), oo(a, b)

    int simple_quad_index_to_order[] =
    {
      oo(1, 1),   oo(1, 1),   oo(1, 2),   XX(1, 3),   oo(1, 4),   XX(1, 5),   oo(1, 6),   XX(1, 7),   oo(1, 8),   XX(1, 9),   oo(1, 10),
      oo(1, 1),   oo(1, 1),   oo(1, 2),   XX(1, 3),   oo(1, 4),   XX(1, 5),   oo(1, 6),   XX(1, 7),   oo(1, 8),   XX(1, 9),   oo(1, 10),
      oo(2, 1),   oo(2, 1),   oo(2, 2),   oo(2, 3),   oo(2, 4),   oo(2, 5),   oo(2, 6),   oo(2, 7),   oo(2, 8),   oo(2, 9),   oo(2, 10),
      XX(3, 1),   XX(3, 1),   oo(3, 2),   oo(3, 3),   oo(3, 4),   oo(3, 5),   oo(3, 6),   oo(3, 7),   oo(3, 8),   oo(3, 9),   oo(3, 10),
      oo(4, 1),   oo(4, 1),   oo(4, 2),   oo(4, 3),   oo(4, 4),   oo(4, 5),   oo(4, 6),   oo(4, 7),   oo(4, 8),   oo(4, 9),   oo(4, 10),
      XX(5, 1),   XX(5, 1),   oo(5, 2),   oo(5, 3),   oo(5, 4),   oo(5, 5),   oo(5, 6),   oo(5, 7),   oo(5, 8),   oo(5, 9),   oo(5, 10),
      oo(6, 1),   oo(6, 1),   oo(6, 2),   oo(6, 3),   oo(6, 4),   oo(6, 5),   oo(6, 6),   oo(6, 7),   oo(6, 8),   oo(6, 9),   oo(6, 10),
      XX(7, 1),   XX(7, 1),   oo(7, 2),   oo(7, 3),   oo(7, 4),   oo(7, 5),   oo(7, 6),   oo(7, 7),   oo(7, 8),   oo(7, 9),   oo(7, 10),
      oo(8, 1),   oo(8, 1),   oo(8, 2),   oo(8, 3),   oo(8, 4),   oo(8, 5),   oo(8, 6),   oo(8, 7),   oo(8, 8),   oo(8, 9),   oo(8, 10),
      XX(9, 1),   XX(9, 1),   oo(9, 2),   oo(9, 3),   oo(9, 4),   oo(9, 5),   oo(9, 6),   oo(9, 7),   oo(9, 8),   oo(9, 9),   oo(9, 10),
      oo(10, 1),  oo(10, 1),  oo(10, 2),  oo(10, 3),  oo(10, 4),  oo(10, 5),  oo(10, 6),  oo(10, 7),  oo(10, 8),  oo(10, 9),  oo(10, 10),
    };
  }
}