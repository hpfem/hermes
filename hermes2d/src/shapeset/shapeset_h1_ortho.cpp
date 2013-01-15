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
#include "shapeset_common.h"
#include "shapeset_h1_all.h"

//// \file triangle ortho2 shapeset /////////////////////////////////////////////////////////////////////
namespace Hermes
{
  namespace Hermes2D
  {
    // ORDER 1

    // vertex functions, order 1

    // number 1
    inline double
    ortho2_f1 (double x, double y)
    {
      return lambda2 (x, y);
    }

    inline double
    ortho2_f1x (double x, double y)
    {
      return lambda2x (x, y);
    }

    inline double
    ortho2_f1y (double x, double y)
    {
      return lambda2y (x, y);
    }

    // number 2
    inline double
    ortho2_f2 (double x, double y)
    {
      return lambda3 (x, y);
    }

    inline double
    ortho2_f2x (double x, double y)
    {
      return lambda3x (x, y);
    }

    inline double
    ortho2_f2y (double x, double y)
    {
      return lambda3y (x, y);
    }

    // number 3
    inline double
    ortho2_f3 (double x, double y)
    {
      return lambda1 (x, y);
    }

    inline double
    ortho2_f3x (double x, double y)
    {
      return lambda1x (x, y);
    }

    inline double
    ortho2_f3y (double x, double y)
    {
      return lambda1y (x, y);
    }

    // ORDER 2

    // Edge functions, order 2

    // number 4
    inline double
    ortho2_f4 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l2 * l3 * phi0 (l3 - l2);
    }

    inline double
    ortho2_f4x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l2x * l3 + l2 * l3x) * phi0 (l3 - l2) + l2 * l3 * phi0x (l3 - l2) * (l3x - l2x);
    }

    inline double
    ortho2_f4y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l2y * l3 + l2 * l3y) * phi0 (l3 - l2) + l2 * l3 * phi0x (l3 - l2) * (l3y - l2y);
    }

    // number 5
    inline double
    ortho2_f5 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l3 * l1 * phi0 (l1 - l3);
    }

    inline double
    ortho2_f5x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l3x * l1 + l3 * l1x) * phi0 (l1 - l3) + l3 * l1 * phi0x (l1 - l3) *  (l1x - l3x);
    }

    inline double
    ortho2_f5y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l3y * l1 + l3 * l1y) * phi0 (l1 - l3) + l3 * l1 * phi0x (l1 - l3) * (l1y - l3y);
    }

    // number 6
    inline double
    ortho2_f6 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l1 * l2 * phi0 (l2 - l1);
    }

    inline double
    ortho2_f6x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l1x * l2 + l1 * l2x) * phi0 (l2 - l1) + l1 * l2 * phi0x (l2 - l1) * (l2x - l1x);
    }

    inline double
    ortho2_f6y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l1y * l2 + l1 * l2y) * phi0 (l2 - l1) + l1 * l2 * phi0x (l2 - l1) * (l2y - l1y);
    }

    // ORDER 3

    // Edge functions, order 3

    // number 7
    inline double
    ortho2_f7_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l2 * l3 * phi1 (l3 - l2);
    }

    inline double
    ortho2_f7_1 (double x, double y)
    {
      return -ortho2_f7_0(x, y);
    }

    inline double
    ortho2_f7x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l2x * l3 + l2 * l3x) * phi1 (l3 - l2) + l2 * l3 * phi1x (l3 - l2) * (l3x - l2x);
    }

    inline double
    ortho2_f7x_1 (double x, double y)
    {
      return -ortho2_f7x_0(x, y);
    }

    inline double
    ortho2_f7y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l2y * l3 + l2 * l3y) * phi1 (l3 - l2) + l2 * l3 * phi1x (l3 - l2) * (l3y - l2y);
    }

    inline double
    ortho2_f7y_1 (double x, double y)
    {
      return -ortho2_f7y_0(x, y);
    }

    // number 8
    inline double
    ortho2_f8_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l3 * l1 * phi1 (l1 - l3);
    }

    inline double
    ortho2_f8_1 (double x, double y)
    {
      return -ortho2_f8_0(x, y);
    }

    inline double
    ortho2_f8x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l3x * l1 + l3 * l1x) * phi1 (l1 - l3) + l3 * l1 * phi1x (l1 - l3) * (l1x - l3x);
    }

    inline double
    ortho2_f8x_1 (double x, double y)
    {
      return -ortho2_f8x_0(x, y);
    }

    inline double
    ortho2_f8y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l3y * l1 + l3 * l1y) * phi1 (l1 - l3) + l3 * l1 * phi1x (l1 - l3) * (l1y - l3y);
    }

    inline double
    ortho2_f8y_1 (double x, double y)
    {
      return -ortho2_f8y_0(x, y);
    }

    // number 9
    inline double
    ortho2_f9_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l1 * l2 * phi1 (l2 - l1);
    }

    inline double
    ortho2_f9_1 (double x, double y)
    {
      return -ortho2_f9_0(x, y);
    }

    inline double
    ortho2_f9x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l1x * l2 + l1 * l2x) * phi1 (l2 - l1) + l1 * l2 * phi1x (l2 - l1) * (l2x - l1x);
    }

    inline double
    ortho2_f9x_1 (double x, double y)
    {
      return -ortho2_f9x_0(x, y);
    }

    inline double
    ortho2_f9y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l1y * l2 + l1 * l2y) * phi1 (l2 - l1) + l1 * l2 * phi1x (l2 - l1) * (l2y - l1y);
    }

    inline double
    ortho2_f9y_1 (double x, double y)
    {
      return -ortho2_f9y_0(x, y);
    }

    // Bubble functions, order 3

    // number 10
    // * f10 *********************************************************************

    double
    ortho2_f10 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t8 = 0.9486832980505138E1 * (0.5 + t3) * (-t3 - t1) * (0.5 + t1);
      return t8;
    }

    double
    ortho2_f10x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t10 =
        -0.4743416490252569E1 * (0.5 + t3) * t2 + 0.4743416490252569E1 * (-t3 - t1) * t2;
      return t10;
    }

    double
    ortho2_f10y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t4 = 0.5 + t1;
      double t10 =
        0.4743416490252569E1 * t4 * (-t1 - t2) -
        0.4743416490252569E1 * t4 * (0.5 + t2);
      return t10;
    }

    // ORDER 4

    // Edge functions, order 4

    // number 11
    inline double
    ortho2_f11 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l2 * l3 * phi2 (l3 - l2);
    }

    inline double
    ortho2_f11x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l2x * l3 + l2 * l3x) * phi2 (l3 - l2) + l2 * l3 * phi2x (l3 - l2) * (l3x - l2x);
    }

    inline double
    ortho2_f11y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l2y * l3 + l2 * l3y) * phi2 (l3 - l2) + l2 * l3 * phi2x (l3 - l2) * (l3y - l2y);
    }

    // number 12
    inline double
    ortho2_f12 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l3 * l1 * phi2 (l1 - l3);
    }

    inline double
    ortho2_f12x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l3x * l1 + l3 * l1x) * phi2 (l1 - l3) + l3 * l1 * phi2x (l1 - l3) * (l1x - l3x);
    }

    inline double
    ortho2_f12y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l3y * l1 + l3 * l1y) * phi2 (l1 - l3) + l3 * l1 * phi2x (l1 - l3) * (l1y - l3y);
    }

    // number 13
    inline double
    ortho2_f13 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l1 * l2 * phi2 (l2 - l1);
    }

    inline double
    ortho2_f13x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l1x * l2 + l1 * l2x) * phi2 (l2 - l1) + l1 * l2 * phi2x (l2 - l1) * (l2x - l1x);
    }

    inline double
    ortho2_f13y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l1y * l2 + l1 * l2y) * phi2 (l2 - l1) + l1 * l2 * phi2x (l2 - l1) * (l2y - l1y);
    }

    // Bubble functions, order 4

    // * f14 *********************************************************************

    double
    ortho2_f14 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t15 =
        -0.7202940575985371E1 * (0.158113883008419E1 * x +
               0.3162277660168379E1 * y +
               0.158113883008419E1) * t6 * t5 +
        0.1626978433639921E1 * t6 * t5;
      return t15;
    }

    double
    ortho2_f14x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t8 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t12 = (-t3 - t1) * t2;
      double t19 =
        0.3601470287992686E1 * t8 * t5 - 0.3601470287992686E1 * t8 * t12 -
        0.1138884903547945E2 * t4 * t12 - 0.8134892168199606 * t5 +
        0.8134892168199606 * t12;
      return t19;
    }

    double
    ortho2_f14y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t8 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t11 = 0.5 + t2;
      double t12 = t4 * t11;
      double t20 =
        -0.3601470287992686E1 * t8 * t5 + 0.3601470287992686E1 * t8 * t12 -
        0.227776980709589E2 * t4 * t3 * t11 + 0.8134892168199606 * t5 -
        0.8134892168199606 * t12;
      return t20;
    }

    // * f15 *********************************************************************

    double
    ortho2_f15 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t21 =
        -0.8906926143924925E1 * (-0.158113883008419E1 - 0.3162277660168379E1 * x -
               0.158113883008419E1 * y) * t6 * t5 -
        0.828416869579514 * t6 * t5 -
        0.5239368319955838E1 * (0.158113883008419E1 * x +
              0.3162277660168379E1 * y +
              0.158113883008419E1) * t6 * t5;
      return t21;
    }

    double
    ortho2_f15x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t8 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t12 = (-t3 - t1) * t2;
      double t21 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t26 =
        0.4453463071962462E1 * t8 * t5 - 0.4453463071962462E1 * t8 * t12 +
        0.1988200486990834E2 * t4 * t12 + 0.414208434789757 * t5 -
        0.414208434789757 * t12 + 0.2619684159977919E1 * t21 * t5 -
        0.2619684159977919E1 * t21 * t12;
      return t26;
    }

    double
    ortho2_f15y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t8 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 = 0.5 + t2;
      double t12 = t4 * t11;
      double t22 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t27 =
        -0.4453463071962462E1 * t8 * t5 + 0.4453463071962462E1 * t8 * t12 -
        0.2485250608738542E1 * t4 * t3 * t11 - 0.414208434789757 * t5 +
        0.414208434789757 * t12 - 0.2619684159977919E1 * t22 * t5 +
        0.2619684159977919E1 * t22 * t12;
      return t27;
    }

    // ORDER 5

    // Edge functions, order 5

    // number 16
    inline double
    ortho2_f16_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l2 * l3 * phi3 (l3 - l2);
    }

    inline double
    ortho2_f16_1 (double x, double y)
    {
      return -ortho2_f16_0(x, y);
    }

    inline double
    ortho2_f16x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l2x * l3 + l2 * l3x) * phi3 (l3 - l2) + l2 * l3 * phi3x (l3 -
                       l2) *
        (l3x - l2x);
    }

    inline double
    ortho2_f16x_1 (double x, double y)
    {
      return -ortho2_f16x_0(x, y);
    }

    inline double
    ortho2_f16y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l2y * l3 + l2 * l3y) * phi3 (l3 - l2) + l2 * l3 * phi3x (l3 -
                       l2) *
        (l3y - l2y);
    }

    inline double
    ortho2_f16y_1 (double x, double y)
    {
      return -ortho2_f16y_0(x, y);
    }

    // number 17
    inline double
    ortho2_f17_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l3 * l1 * phi3 (l1 - l3);
    }

    inline double
    ortho2_f17_1 (double x, double y)
    {
      return -ortho2_f17_0(x, y);
    }

    inline double
    ortho2_f17x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l3x * l1 + l3 * l1x) * phi3 (l1 - l3) + l3 * l1 * phi3x (l1 -
                       l3) *
        (l1x - l3x);
    }

    inline double
    ortho2_f17x_1 (double x, double y)
    {
      return -ortho2_f17x_0(x, y);
    }

    inline double
    ortho2_f17y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l3y * l1 + l3 * l1y) * phi3 (l1 - l3) + l3 * l1 * phi3x (l1 -
                       l3) *
        (l1y - l3y);
    }

    inline double
    ortho2_f17y_1 (double x, double y)
    {
      return -ortho2_f17y_0(x, y);
    }

    // number 18
    inline double
    ortho2_f18_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l1 * l2 * phi3 (l2 - l1);
    }

    inline double
    ortho2_f18_1 (double x, double y)
    {
      return -ortho2_f18_0(x, y);
    }

    inline double
    ortho2_f18x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l1x * l2 + l1 * l2x) * phi3 (l2 - l1) + l1 * l2 * phi3x (l2 -
                       l1) *
        (l2x - l1x);
    }

    inline double
    ortho2_f18x_1 (double x, double y)
    {
      return -ortho2_f18x_0(x, y);
    }

    inline double
    ortho2_f18y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l1y * l2 + l1 * l2y) * phi3 (l2 - l1) + l1 * l2 * phi3x (l2 -
                       l1) *
        (l2y - l1y);
    }

    inline double
    ortho2_f18y_1 (double x, double y)
    {
      return -ortho2_f18y_0(x, y);
    }

    // Bubble functions, order 5

    // * f19 *********************************************************************

    double
    ortho2_f19 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * y;
      double t28 =
        0.5071884024840934E2 * (-t3 - t7 - 0.9472135954999579) * (-t3 - t7 -
                        0.5278640450004206E-1)
        * t6 * t5 + 0.5682388583386602E1 * t6 * t5 -
        0.2153344726179916E1 * (0.158113883008419E1 * x +
              0.3162277660168379E1 * y +
              0.158113883008419E1) * t6 * t5 -
        0.1856331660499927E1 * (-0.158113883008419E1 - 0.3162277660168379E1 * x -
              0.158113883008419E1 * y) * t6 * t5;
      return t28;
    }

    double
    ortho2_f19x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * y;
      double t7 = -t3 - t6 - 0.5278640450004206E-1;
      double t8 = -t3 - t6 - 0.9472135954999579;
      double t9 = t8 * t7;
      double t13 = (-t3 - t1) * t2;
      double t26 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t35 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t40 =
        -0.2535942012420467E2 * t9 * t5 + 0.2535942012420467E2 * t9 * t13 -
        0.2535942012420467E2 * t8 * t4 * t13 -
        0.2535942012420467E2 * t7 * t4 * t13 - 0.2841194291693301E1 * t5 +
        0.2841194291693301E1 * t13 + 0.1076672363089958E1 * t26 * t5 -
        0.1076672363089958E1 * t26 * t13 + 0.2465499178742121E1 * t4 * t13 +
        0.9281658302499636 * t35 * t5 - 0.9281658302499636 * t35 * t13;
      return t40;
    }

    double
    ortho2_f19y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t6 = 0.1E1 * y;
      double t7 = -t1 - t6 - 0.5278640450004206E-1;
      double t8 = -t1 - t6 - 0.9472135954999579;
      double t9 = t8 * t7;
      double t12 = 0.5 + t2;
      double t13 = t4 * t12;
      double t16 = t3 * t12;
      double t27 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t36 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t41 =
        0.2535942012420467E2 * t9 * t5 - 0.2535942012420467E2 * t9 * t13 -
        0.5071884024840934E2 * t8 * t4 * t16 -
        0.5071884024840934E2 * t7 * t4 * t16 + 0.2841194291693301E1 * t5 -
        0.2841194291693301E1 * t13 - 0.1076672363089958E1 * t27 * t5 +
        0.1076672363089958E1 * t27 * t13 - 0.3874355852309047E1 * t4 * t16 -
        0.9281658302499636 * t36 * t5 + 0.9281658302499636 * t36 * t13;
      return t41;
    }

    // * f20 *********************************************************************

    double
    ortho2_f20 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t10 =
        (-0.158113883008419E1 - 0.3162277660168379E1 * x -
         0.158113883008419E1 * y) * t6;
      double t13 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t24 = 0.1E1 * y;
      double t31 =
        0.8791860079241937E1 * t13 * t10 * t5 + 0.8739923287658513E1 * t6 * t5 -
        0.2335472262979457E1 * t13 * t6 * t5 - 0.153398811084586E1 * t10 * t5 +
        0.4894523238821834E2 * (-t3 - t24 - 0.9472135954999579) * (-t3 - t24 -
                         0.5278640450004206E-1)
        * t6 * t5;
      return t31;
    }

    double
    ortho2_f20x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t8 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t12 = t11 * t8;
      double t16 = (-t3 - t1) * t2;
      double t37 = 0.1E1 * y;
      double t38 = -t3 - t37 - 0.5278640450004206E-1;
      double t39 = -t3 - t37 - 0.9472135954999579;
      double t40 = t39 * t38;
      double t51 =
        -0.4395930039620968E1 * t12 * t5 + 0.4395930039620968E1 * t12 * t16 -
        0.2780230271991297E2 * t11 * t4 * t16 +
        0.1390115135995649E2 * t8 * t4 * t16 - 0.4369961643829257E1 * t5 +
        0.4369961643829257E1 * t16 + 0.1167736131489729E1 * t11 * t5 -
        0.1167736131489729E1 * t11 * t16 + 0.1158190452310345E1 * t4 * t16 +
        0.76699405542293 * t8 * t5 - 0.76699405542293 * t8 * t16 -
        0.2447261619410917E2 * t40 * t5 + 0.2447261619410917E2 * t40 * t16 -
        0.2447261619410917E2 * t39 * t4 * t16 -
        0.2447261619410917E2 * t38 * t4 * t16;
      return t51;
    }

    double
    ortho2_f20y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t8 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t12 = t11 * t8;
      double t15 = 0.5 + t2;
      double t16 = t4 * t15;
      double t19 = t3 * t15;
      double t38 = 0.1E1 * y;
      double t39 = -t1 - t38 - 0.5278640450004206E-1;
      double t40 = -t1 - t38 - 0.9472135954999579;
      double t41 = t40 * t39;
      double t52 =
        0.4395930039620968E1 * t12 * t5 - 0.4395930039620968E1 * t12 * t16 -
        0.1390115135995649E2 * t11 * t4 * t19 +
        0.2780230271991297E2 * t8 * t4 * t19 + 0.4369961643829257E1 * t5 -
        0.4369961643829257E1 * t16 - 0.1167736131489729E1 * t11 * t5 +
        0.1167736131489729E1 * t11 * t16 - 0.4959963596216948E1 * t4 * t19 -
        0.76699405542293 * t8 * t5 + 0.76699405542293 * t8 * t16 +
        0.2447261619410917E2 * t41 * t5 - 0.2447261619410917E2 * t41 * t16 -
        0.4894523238821834E2 * t40 * t4 * t19 -
        0.4894523238821834E2 * t39 * t4 * t19;
      return t52;
    }

    // * f21 *********************************************************************

    double
    ortho2_f21 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * x;
      double t18 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t25 =
        (-0.158113883008419E1 - 0.3162277660168379E1 * x -
         0.158113883008419E1 * y) * t6;
      double t28 = 0.1E1 * y;
      double t38 =
        0.7468713169834236E2 * (0.5278640450004206E-1 + t7 +
              t1) * (0.9472135954999579 + t7 + t1) * t6 * t5 +
        0.1432355950379169E2 * t6 * t5 + 0.143280330532747E1 * t18 * t6 * t5 +
        0.2518960649688617E1 * t25 * t5 + 0.2470083220551831E2 * (-t3 - t28 -
                        0.9472135954999579)
        * (-t3 - t28 - 0.5278640450004206E-1) * t6 * t5 +
        0.8608751579319697E1 * t18 * t25 * t5;
      return t38;
    }

    double
    ortho2_f21x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * x;
      double t7 = 0.9472135954999579 + t6 + t1;
      double t8 = 0.5278640450004206E-1 + t6 + t1;
      double t9 = t8 * t7;
      double t13 = (-t3 - t1) * t2;
      double t26 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t35 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t40 = 0.1E1 * y;
      double t41 = -t3 - t40 - 0.5278640450004206E-1;
      double t42 = -t3 - t40 - 0.9472135954999579;
      double t43 = t42 * t41;
      double t54 = t26 * t35;
      double t65 =
        -0.3734356584917118E2 * t9 * t5 + 0.3734356584917118E2 * t9 * t13 +
        0.7468713169834236E2 * t8 * t4 * t13 +
        0.7468713169834236E2 * t7 * t4 * t13 - 0.7161779751895843E1 * t5 +
        0.7161779751895843E1 * t13 - 0.716401652663735 * t26 * t5 +
        0.716401652663735 * t26 * t13 - 0.5700192047427303E1 * t4 * t13 -
        0.1259480324844308E1 * t35 * t5 + 0.1259480324844308E1 * t35 * t13 -
        0.1235041610275916E2 * t43 * t5 + 0.1235041610275916E2 * t43 * t13 -
        0.1235041610275916E2 * t42 * t4 * t13 -
        0.1235041610275916E2 * t41 * t4 * t13 - 0.4304375789659848E1 * t54 * t5 +
        0.4304375789659848E1 * t54 * t13 - 0.2722326280122193E2 * t26 * t4 * t13 +
        0.1361163140061097E2 * t35 * t4 * t13;
      return t65;
    }

    double
    ortho2_f21y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t6 = 0.1E1 * x;
      double t7 = 0.9472135954999579 + t6 + t2;
      double t8 = 0.5278640450004206E-1 + t6 + t2;
      double t9 = t8 * t7;
      double t12 = 0.5 + t2;
      double t13 = t4 * t12;
      double t16 = t3 * t12;
      double t27 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t36 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t41 = 0.1E1 * y;
      double t42 = -t1 - t41 - 0.5278640450004206E-1;
      double t43 = -t1 - t41 - 0.9472135954999579;
      double t44 = t43 * t42;
      double t55 = t27 * t36;
      double t66 =
        0.3734356584917118E2 * t9 * t5 - 0.3734356584917118E2 * t9 * t13 +
        0.3734356584917118E2 * t8 * t4 * t16 +
        0.3734356584917118E2 * t7 * t4 * t16 + 0.7161779751895843E1 * t5 -
        0.7161779751895843E1 * t13 + 0.716401652663735 * t27 * t5 -
        0.716401652663735 * t27 * t13 + 0.5480953891757022 * t4 * t16 +
        0.1259480324844308E1 * t36 * t5 - 0.1259480324844308E1 * t36 * t13 +
        0.1235041610275916E2 * t44 * t5 - 0.1235041610275916E2 * t44 * t13 -
        0.2470083220551831E2 * t43 * t4 * t16 -
        0.2470083220551831E2 * t42 * t4 * t16 + 0.4304375789659848E1 * t55 * t5 -
        0.4304375789659848E1 * t55 * t13 - 0.1361163140061097E2 * t27 * t4 * t16 +
        0.2722326280122193E2 * t36 * t4 * t16;
      return t66;
    }

    // ORDER 6

    // Edge functions, order 6

    // number 22
    inline double
    ortho2_f22 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l2 * l3 * phi4 (l3 - l2);
    }

    inline double
    ortho2_f22x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l2x * l3 + l2 * l3x) * phi4 (l3 - l2) + l2 * l3 * phi4x (l3 -
                       l2) *
        (l3x - l2x);
    }

    inline double
    ortho2_f22y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l2y * l3 + l2 * l3y) * phi4 (l3 - l2) + l2 * l3 * phi4x (l3 -
                       l2) *
        (l3y - l2y);
    }

    // number 23
    inline double
    ortho2_f23 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l3 * l1 * phi4 (l1 - l3);
    }

    inline double
    ortho2_f23x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l3x * l1 + l3 * l1x) * phi4 (l1 - l3) + l3 * l1 * phi4x (l1 -
                       l3) *
        (l1x - l3x);
    }

    inline double
    ortho2_f23y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l3y * l1 + l3 * l1y) * phi4 (l1 - l3) + l3 * l1 * phi4x (l1 -
                       l3) *
        (l1y - l3y);
    }

    // number 24
    inline double
    ortho2_f24 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l1 * l2 * phi4 (l2 - l1);
    }

    inline double
    ortho2_f24x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l1x * l2 + l1 * l2x) * phi4 (l2 - l1) + l1 * l2 * phi4x (l2 -
                       l1) *
        (l2x - l1x);
    }

    inline double
    ortho2_f24y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l1y * l2 + l1 * l2y) * phi4 (l2 - l1) + l1 * l2 * phi4x (l2 -
                       l1) *
        (l2y - l1y);
    }

    // Bubble functions, order 6

    // * f25 *********************************************************************

    double
    ortho2_f25 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t19 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t26 =
        (-0.158113883008419E1 - 0.3162277660168379E1 * x -
         0.158113883008419E1 * y) * t6;
      double t38 = 0.1E1 * x;
      double t45 =
        0.1102119728831204E3 * (-t3 - t8 - 0.1154653670707977E1) * (-t3 - t8 -
                    0.5) * (-t3 -
                      t8 +
                      0.1546536707079771)
        * t7 + 0.4359261447370632E1 * t7 - 0.8084773402927861E1 * t19 * t6 * t5 -
        0.1630506946111695 * t26 * t5 + 0.2537043288590704E2 * (-t3 - t8 -
                      0.9472135954999579)
        * (-t3 - t8 - 0.5278640450004206E-1) * t6 * t5 +
        0.3339756754033952E1 * t19 * t26 * t5 -
        0.2343688950199264 * (0.5278640450004206E-1 + t38 +
            t1) * (0.9472135954999579 + t38 + t1) * t6 * t5;
      return t45;
    }

    double
    ortho2_f25x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t8 = (-t3 - t1) * t2;
      double t10 = 0.1E1 * y;
      double t11 = -t3 - t10 - 0.5278640450004206E-1;
      double t12 = -t3 - t10 - 0.9472135954999579;
      double t13 = t12 * t11;
      double t26 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t32 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t36 = t26 * t32;
      double t45 = -0.2179630723685316E1 * t5 + 0.2179630723685316E1 * t8
        - 0.1268521644295352E2 * t13 * t5 + 0.1268521644295352E2 * t13 * t8 -
        0.1268521644295352E2 * t12 * t4 * t8 -
        0.1268521644295352E2 * t11 * t4 * t8 -
        0.1056123817367803E2 * t26 * t4 * t8 +
        0.5280619086839013E1 * t32 * t4 * t8 - 0.1669878377016976E1 * t36 * t5 +
        0.1669878377016976E1 * t36 * t8 + 0.8152534730558477E-1 * t32 * t5 -
        0.8152534730558477E-1 * t32 * t8;
      double t50 = -t3 - t10 - 0.5;
      double t52 = -t3 - t10 - 0.1154653670707977E1;
      double t56 = -t3 - t10 + 0.1546536707079771;
      double t57 = t56 * t4;
      double t62 = t52 * t50 * t56;
      double t72 = 0.1E1 * x;
      double t73 = 0.9472135954999579 + t72 + t1;
      double t74 = 0.5278640450004206E-1 + t72 + t1;
      double t75 = t74 * t73;
      double t86 = 0.404238670146393E1 * t26 * t5 - 0.404238670146393E1 * t26 * t8
        - 0.551059864415602E2 * t52 * t50 * t4 * t8 -
        0.551059864415602E2 * t52 * t57 * t8 - 0.551059864415602E2 * t62 * t5 +
        0.551059864415602E2 * t62 * t8 - 0.1226753759075729E2 * t4 * t8 -
        0.551059864415602E2 * t50 * t57 * t8 + 0.1171844475099632 * t75 * t5 -
        0.2343688950199264 * t73 * t4 * t8 - 0.2343688950199264 * t74 * t4 * t8 -
        0.1171844475099632 * t75 * t8;
      double t87 = t45 + t86;
      return t87;
    }

    double
    ortho2_f25y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t6 = 0.1E1 * y;
      double t7 = -t1 - t6 + 0.1546536707079771;
      double t8 = -t1 - t6 - 0.5;
      double t10 = -t1 - t6 - 0.1154653670707977E1;
      double t11 = t10 * t8 * t7;
      double t14 = 0.1E1 * x;
      double t15 = 0.9472135954999579 + t14 + t2;
      double t16 = 0.5278640450004206E-1 + t14 + t2;
      double t17 = t16 * t15;
      double t20 = 0.5 + t2;
      double t21 = t4 * t20;
      double t24 = -t1 - t6 - 0.5278640450004206E-1;
      double t25 = -t1 - t6 - 0.9472135954999579;
      double t26 = t25 * t24;
      double t29 = t3 * t20;
      double t40 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t46 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t50 = t40 * t46;
      double t55 = 0.551059864415602E2 * t11 * t5 - 0.1171844475099632 * t17 * t5
        - 0.2179630723685316E1 * t21 + 0.2179630723685316E1 * t5 -
        0.1268521644295352E2 * t26 * t21 - 0.2537043288590704E2 * t25 * t4 * t29 -
        0.2537043288590704E2 * t24 * t4 * t29 + 0.1268521644295352E2 * t26 * t5 -
        0.5280619086839013E1 * t40 * t4 * t29 +
        0.1056123817367803E2 * t46 * t4 * t29 - 0.1669878377016976E1 * t50 * t21 +
        0.8152534730558477E-1 * t46 * t21;
      double t64 = t7 * t4;
      double t87 =
        0.404238670146393E1 * t40 * t21 - 0.8152534730558477E-1 * t46 * t5 -
        0.1102119728831204E3 * t10 * t4 * t8 * t29 -
        0.1102119728831204E3 * t10 * t64 * t29 - 0.551059864415602E2 * t11 * t21 -
        0.2530849253508034E2 * t4 * t29 + 0.1669878377016976E1 * t50 * t5 -
        0.1102119728831204E3 * t8 * t64 * t29 - 0.404238670146393E1 * t40 * t5 +
        0.1171844475099632 * t17 * t21 - 0.1171844475099632 * t15 * t4 * t29 -
        0.1171844475099632 * t16 * t4 * t29;
      double t88 = t55 + t87;
      return t88;
    }

    // * f26 *********************************************************************

    double
    ortho2_f26 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t10 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 = 0.1E1 * y;
      double t12 = -t3 - t11 - 0.5278640450004206E-1;
      double t14 = -t3 - t11 - 0.9472135954999579;
      double t18 = 0.1E1 * x;
      double t29 = t10 * t6;
      double t32 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t49 =
        -0.750062486130387E2 * t14 * t12 * t10 * t7 +
        0.2431428729906529E2 * (0.5278640450004206E-1 + t18 +
              t1) * (0.9472135954999579 + t18 + t1) * t6 * t5 +
        0.4475799965898882E2 * t14 * t12 * t6 * t5 +
        0.4156173881104336E1 * t32 * t29 * t5 -
        0.1408679299827772E2 * t32 * t6 * t5 - 0.9829176222510025E1 * t29 * t5 +
        0.1292473852643748E3 * (-t3 - t11 - 0.1154653670707977E1) * (-t3 - t11 -
                     0.5) * (-t3 -
                       t11 +
                       0.1546536707079771)
        * t7 + 0.118332076609073E2 * t7;
      return t49;
    }

    double
    ortho2_f26x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t8 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t12 = t11 * t8;
      double t15 = 0.1E1 * y;
      double t16 = -t3 - t15 - 0.5278640450004206E-1;
      double t17 = -t3 - t15 - 0.9472135954999579;
      double t18 = t17 * t16;
      double t22 = t17 * t16 * t8;
      double t26 = (-t3 - t1) * t2;
      double t29 = t16 * t4;
      double t39 = t8 * t4;
      double t53 =
        -0.2078086940552168E1 * t12 * t5 - 0.2237899982949441E2 * t18 * t5 +
        0.3750312430651935E2 * t22 * t5 + 0.2078086940552168E1 * t12 * t26 -
        0.2237899982949441E2 * t29 * t26 - 0.2237899982949441E2 * t17 * t4 * t26 -
        0.3750312430651935E2 * t22 * t26 + 0.2237899982949441E2 * t18 * t26 +
        0.6571487907995775E1 * t39 * t26 - 0.1314297581599155E2 * t11 * t4 * t26 +
        0.3750312430651935E2 * t17 * t39 * t26 +
        0.3750312430651935E2 * t16 * t39 * t26 - 0.591660383045365E1 * t5 +
        0.591660383045365E1 * t26;
      double t54 = -t3 - t15 - 0.5;
      double t56 = -t3 - t15 - 0.1154653670707977E1;
      double t60 = -t3 - t15 + 0.1546536707079771;
      double t61 = t60 * t4;
      double t66 = t56 * t54 * t60;
      double t69 = 0.1E1 * x;
      double t70 = 0.9472135954999579 + t69 + t1;
      double t71 = 0.5278640450004206E-1 + t69 + t1;
      double t72 = t71 * t70;
      double t101 =
        -0.6462369263218741E2 * t56 * t54 * t4 * t26 -
        0.6462369263218741E2 * t56 * t61 * t26 - 0.6462369263218741E2 * t66 * t5 +
        0.1215714364953265E2 * t72 * t26 - 0.4914588111255012E1 * t8 * t26 +
        0.4914588111255012E1 * t8 * t5 + 0.2371905843620478E3 * t17 * t29 * t26 +
        0.7043396499138861E1 * t11 * t5 + 0.8809408985366679E1 * t4 * t26 -
        0.7043396499138861E1 * t11 * t26 + 0.6462369263218741E2 * t66 * t26 -
        0.6462369263218741E2 * t54 * t61 * t26 +
        0.2431428729906529E2 * t71 * t4 * t26 - 0.1215714364953265E2 * t72 * t5 +
        0.2431428729906529E2 * t70 * t4 * t26;
      double t102 = t53 + t101;
      return t102;
    }

    double
    ortho2_f26y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = 0.1E1 * y;
      double t12 = -t3 - t11 - 0.9472135954999579;
      double t16 = -t3 - t11 - 0.5278640450004206E-1;
      double t20 = t6 * t2;
      double t23 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t26 = t6 * t4;
      double t32 = -t3 - t11 - 0.5;
      double t34 = -t3 - t11 - 0.1154653670707977E1;
      double t39 = t12 * t16 * t9;
      double t42 = -t3 - t11 + 0.1546536707079771;
      double t43 = t42 * t6;
      double t50 = t12 * t16;
      double t53 = 0.1E1 * x;
      double t54 = 0.9472135954999579 + t53 + t1;
      double t55 = 0.5278640450004206E-1 + t53 + t1;
      double t56 = t55 * t54;
      double t60 = t34 * t32 * t42;
      double t63 = t16 * t6;
      double t66 =
        0.750062486130387E2 * t12 * t10 * t5 +
        0.750062486130387E2 * t16 * t10 * t5 + 0.7043396499138861E1 * t23 * t20 -
        0.4914588111255012E1 * t9 * t26 - 0.2900505860871915E2 * t6 * t5 +
        0.591660383045365E1 * t26 - 0.1292473852643748E3 * t34 * t32 * t6 * t5 -
        0.3750312430651935E2 * t39 * t26 - 0.1292473852643748E3 * t34 * t43 * t5 -
        0.4475799965898882E2 * t12 * t6 * t5 - 0.2237899982949441E2 * t50 * t20 +
        0.1215714364953265E2 * t56 * t26 - 0.6462369263218741E2 * t60 * t20 -
        0.4475799965898882E2 * t63 * t5;
      double t67 = t23 * t9;
      double t102 =
        -0.2078086940552168E1 * t67 * t20 + 0.1314297581599155E2 * t10 * t5 -
        0.6571487907995775E1 * t23 * t6 * t5 + 0.2237899982949441E2 * t50 * t26 +
        0.2078086940552168E1 * t67 * t26 + 0.1215714364953265E2 * t55 * t6 * t5 +
        0.1215714364953265E2 * t54 * t6 * t5 - 0.1215714364953265E2 * t56 * t20 -
        0.591660383045365E1 * t20 - 0.7043396499138861E1 * t23 * t26 +
        0.3750312430651935E2 * t39 * t20 + 0.4914588111255012E1 * t9 * t20 -
        0.1292473852643748E3 * t32 * t43 * t5 + 0.6462369263218741E2 * t60 * t26 +
        0.1185952921810239E3 * t12 * t63 * t5;
      double t103 = t66 + t102;
      return t103;
    }

    // * f27 *********************************************************************

    double
    ortho2_f27 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t10 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 = 0.1E1 * y;
      double t12 = -t3 - t11 - 0.5278640450004206E-1;
      double t14 = -t3 - t11 - 0.9472135954999579;
      double t18 = 0.1E1 * x;
      double t19 = 0.9472135954999579 + t18 + t1;
      double t21 = 0.5278640450004206E-1 + t18 + t1;
      double t29 = t10 * t6;
      double t32 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t53 =
        -0.9654367571766702E2 * t14 * t12 * t10 * t7 +
        0.1401882261448316E2 * t21 * t19 * t6 * t5 +
        0.1009339240923249E2 * t14 * t12 * t6 * t5 +
        0.2580465905695021E1 * t32 * t29 * t5 -
        0.2331196977901257E2 * t32 * t6 * t5 - 0.1788640850123086E2 * t29 * t5 -
        0.8723366372133481E2 * t32 * t21 * t19 * t7 +
        0.8417401345620565E2 * (-t3 - t11 - 0.1154653670707977E1) * (-t3 - t11 -
                     0.5) * (-t3 -
                       t11 +
                       0.1546536707079771)
        * t7 + 0.4163606342625211E1 * t7;
      return t53;
    }

    double
    ortho2_f27x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t7 = 0.5 + t3;
      double t8 = t7 * t2;
      double t12 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t13 = 0.1E1 * y;
      double t14 = -t3 - t13 - 0.5278640450004206E-1;
      double t16 = -t3 - t13 - 0.9472135954999579;
      double t17 = t16 * t14 * t12;
      double t20 = t16 * t14;
      double t25 = t12 * t7;
      double t30 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t34 = t14 * t7;
      double t40 = t30 * t12;
      double t45 = 0.1E1 * x;
      double t46 = 0.9472135954999579 + t45 + t1;
      double t47 = 0.5278640450004206E-1 + t45 + t1;
      double t48 = t47 * t46;
      double t53 = t47 * t7;
      double t56 = t46 * t7;
      double t61 = t30 * t48;
      double t64 = 0.2081803171312606E1 * t5 - 0.2081803171312606E1 * t8 +
        0.4827183785883351E2 * t17 * t8 - 0.5046696204616243E1 * t20 * t8 +
        0.5046696204616243E1 * t20 * t5 + 0.4080074843202764E1 * t25 * t5 -
        0.8160149686405528E1 * t30 * t7 * t5 - 0.5046696204616243E1 * t34 * t5 -
        0.5046696204616243E1 * t16 * t7 * t5 + 0.129023295284751E1 * t40 * t5 -
        0.129023295284751E1 * t40 * t8 + 0.700941130724158E1 * t48 * t5 -
        0.4827183785883351E2 * t17 * t5 + 0.1401882261448316E2 * t53 * t5 +
        0.1401882261448316E2 * t56 * t5 - 0.700941130724158E1 * t48 * t8 +
        0.436168318606674E2 * t61 * t8;
      double t89 = -t3 - t13 - 0.5;
      double t91 = -t3 - t13 - 0.1154653670707977E1;
      double t95 = -t3 - t13 + 0.1546536707079771;
      double t96 = t95 * t7;
      double t101 = t91 * t89 * t95;
      double t115 =
        -0.436168318606674E2 * t61 * t5 + 0.4827183785883351E2 * t16 * t25 * t5 +
        0.4827183785883351E2 * t14 * t25 * t5 -
        0.8723366372133481E2 * t30 * t53 * t5 -
        0.8723366372133481E2 * t30 * t56 * t5 + 0.1165598488950628E2 * t30 * t8 -
        0.894320425061543E1 * t12 * t5 + 0.894320425061543E1 * t12 * t8 -
        0.1165598488950628E2 * t30 * t5 + 0.1970232940074221E2 * t7 * t5 -
        0.4208700672810283E2 * t91 * t89 * t7 * t5 -
        0.4208700672810283E2 * t91 * t96 * t5 - 0.4208700672810283E2 * t101 * t8 +
        0.4208700672810283E2 * t101 * t5 - 0.4208700672810283E2 * t89 * t96 * t5 -
        0.1379285330003089E3 * t47 * t56 * t5 +
        0.3052979089525188E3 * t16 * t34 * t5;
      double t116 = t64 + t115;
      return t116;
    }

    double
    ortho2_f27y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t7 = 0.5 + t2;
      double t8 = t4 * t7;
      double t10 = t3 * t7;
      double t11 = 0.1E1 * y;
      double t12 = -t1 - t11 - 0.5278640450004206E-1;
      double t13 = t12 * t4;
      double t16 = -t1 - t11 - 0.9472135954999579;
      double t20 = t16 * t12;
      double t23 = 0.1E1 * x;
      double t24 = 0.9472135954999579 + t23 + t2;
      double t25 = 0.5278640450004206E-1 + t23 + t2;
      double t26 = t25 * t24;
      double t31 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t39 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t40 = t31 * t39;
      double t43 = t39 * t4;
      double t48 = t25 * t4;
      double t51 = t24 * t4;
      double t56 = -t1 - t11 + 0.1546536707079771;
      double t57 = -t1 - t11 - 0.5;
      double t59 = -t1 - t11 - 0.1154653670707977E1;
      double t60 = t59 * t57 * t56;
      double t68 = 0.2081803171312606E1 * t5 - 0.2081803171312606E1 * t8
        - 0.1009339240923249E2 * t13 * t10 - 0.1009339240923249E2 * t16 * t4 * t10
        - 0.5046696204616243E1 * t20 * t8 + 0.700941130724158E1 * t26 * t5 -
        0.4080074843202764E1 * t31 * t4 * t10 + 0.5046696204616243E1 * t20 * t5 -
        0.129023295284751E1 * t40 * t8 + 0.8160149686405528E1 * t43 * t10 +
        0.129023295284751E1 * t40 * t5 + 0.700941130724158E1 * t48 * t10 +
        0.700941130724158E1 * t51 * t10 - 0.700941130724158E1 * t26 * t8 +
        0.4208700672810283E2 * t60 * t5 - 0.2758570660006179E3 * t25 * t51 * t10 -
        0.4543802623464777E2 * t4 * t10;
      double t73 = t16 * t12 * t39;
      double t86 = t31 * t26;
      double t93 = t56 * t4;
      double t116 =
        0.1526489544762594E3 * t16 * t13 * t10 - 0.4827183785883351E2 * t73 * t5 +
        0.1165598488950628E2 * t31 * t8 - 0.894320425061543E1 * t39 * t5 -
        0.436168318606674E2 * t31 * t51 * t10 -
        0.436168318606674E2 * t31 * t48 * t10 + 0.436168318606674E2 * t86 * t8 +
        0.894320425061543E1 * t39 * t8 - 0.4208700672810283E2 * t60 * t8 -
        0.8417401345620565E2 * t59 * t93 * t10 -
        0.8417401345620565E2 * t59 * t57 * t4 * t10 -
        0.1165598488950628E2 * t31 * t5 - 0.436168318606674E2 * t86 * t5 -
        0.8417401345620565E2 * t57 * t93 * t10 +
        0.9654367571766702E2 * t12 * t43 * t10 +
        0.9654367571766702E2 * t16 * t43 * t10 + 0.4827183785883351E2 * t73 * t8;
      double t117 = t68 + t116;
      return t117;
    }

    // * f28 *********************************************************************

    double
    ortho2_f28 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t10 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 = 0.1E1 * y;
      double t12 = -t3 - t11 - 0.5278640450004206E-1;
      double t14 = -t3 - t11 - 0.9472135954999579;
      double t18 = 0.1E1 * x;
      double t19 = 0.9472135954999579 + t18 + t1;
      double t21 = 0.5278640450004206E-1 + t18 + t1;
      double t29 = t10 * t6;
      double t32 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t60 =
        -0.5516494482649984E2 * t14 * t12 * t10 * t7 -
        0.466228674418177E2 * t21 * t19 * t6 * t5 -
        0.1458234147916518E2 * t14 * t12 * t6 * t5 -
        0.5004661733822727E1 * t32 * t29 * t5 -
        0.1911261620960439E2 * t32 * t6 * t5 - 0.2749506194739649E2 * t29 * t5 -
        0.1018826663546988E3 * t32 * t21 * t19 * t7 +
        0.192256519358649E3 * (-0.1546536707079771 + t18 + t1) * (0.5 + t18 +
                        t1) *
        (0.1154653670707977E1 + t18 + t1) * t7 + 0.3202716237991868E2 * (-t3 -
                         t11 -
                         0.1154653670707977E1)
        * (-t3 - t11 - 0.5) * (-t3 - t11 + 0.1546536707079771) * t7 -
        0.1052382140289826E2 * t7;
      return t60;
    }

    double
    ortho2_f28x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t7 = 0.5 + t3;
      double t8 = t7 * t2;
      double t10 = 0.1E1 * y;
      double t11 = -t3 - t10 - 0.5;
      double t13 = -t3 - t10 - 0.1154653670707977E1;
      double t17 = -t3 - t10 + 0.1546536707079771;
      double t18 = t17 * t7;
      double t22 = -t3 - t10 - 0.5278640450004206E-1;
      double t23 = -t3 - t10 - 0.9472135954999579;
      double t24 = t23 * t22;
      double t30 = t13 * t11 * t17;
      double t33 = t22 * t7;
      double t41 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t42 = t41 * t7;
      double t49 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t53 = t49 * t41;
      double t63 = 0.1E1 * x;
      double t64 = 0.9472135954999579 + t63 + t1;
      double t65 = 0.5278640450004206E-1 + t63 + t1;
      double t66 = t65 * t64;
      double t69 = t65 * t7;
      double t72 = t64 * t7;
      double t75 = -0.5261910701449129E1 * t5 + 0.5261910701449129E1 * t8
        - 0.1601358118995934E2 * t13 * t11 * t7 * t5 -
        0.1601358118995934E2 * t13 * t18 * t5 - 0.729117073958259E1 * t24 * t5 +
        0.729117073958259E1 * t24 * t8 - 0.1601358118995934E2 * t30 * t8 +
        0.729117073958259E1 * t33 * t5 + 0.729117073958259E1 * t23 * t7 * t5 -
        0.7913064998783579E1 * t42 * t5 + 0.1601358118995934E2 * t30 * t5 +
        0.1582612999756716E2 * t49 * t7 * t5 - 0.2502330866911363E1 * t53 * t5 +
        0.2502330866911363E1 * t53 * t8 - 0.1601358118995934E2 * t11 * t18 * t5 +
        0.5672732052769561E2 * t7 * t5 - 0.2331143372090885E2 * t66 * t5 -
        0.466228674418177E2 * t69 * t5 - 0.466228674418177E2 * t72 * t5;
      double t78 = 0.1154653670707977E1 + t63 + t1;
      double t79 = 0.5 + t63 + t1;
      double t81 = -0.1546536707079771 + t63 + t1;
      double t82 = t81 * t79 * t78;
      double t97 = t78 * t7;
      double t108 = t23 * t22 * t41;
      double t122 = t49 * t66;
      double t135 =
        0.2331143372090885E2 * t66 * t8 + 0.9612825967932449E2 * t82 * t5 -
        0.1374753097369824E2 * t41 * t5 + 0.1374753097369824E2 * t41 * t8 -
        0.9556308104802193E1 * t49 * t5 + 0.9556308104802193E1 * t49 * t8 +
        0.192256519358649E3 * t81 * t79 * t7 * t5 +
        0.192256519358649E3 * t81 * t97 * t5 +
        0.192256519358649E3 * t79 * t97 * t5 -
        0.1610906398859262E3 * t65 * t72 * t5 + 0.2758247241324992E2 * t108 * t8 +
        0.1744468726492617E3 * t23 * t33 * t5 - 0.2758247241324992E2 * t108 * t5 +
        0.2758247241324992E2 * t23 * t42 * t5 +
        0.2758247241324992E2 * t22 * t42 * t5 + 0.5094133317734938E2 * t122 * t8 -
        0.5094133317734938E2 * t122 * t5 - 0.1018826663546988E3 * t49 * t69 * t5 -
        0.1018826663546988E3 * t49 * t72 * t5 - 0.9612825967932449E2 * t82 * t8;
      double t136 = t75 + t135;
      return t136;
    }

    double
    ortho2_f28y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * x;
      double t8 = 0.5278640450004206E-1 + t7 + t1;
      double t9 = t8 * t6;
      double t12 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t16 = t6 * t4;
      double t17 = 0.1154653670707977E1 + t7 + t1;
      double t18 = 0.5 + t7 + t1;
      double t20 = -0.1546536707079771 + t7 + t1;
      double t21 = t20 * t18 * t17;
      double t24 = t6 * t2;
      double t25 = 0.1E1 * y;
      double t26 = -t3 - t25 + 0.1546536707079771;
      double t27 = -t3 - t25 - 0.5;
      double t29 = -t3 - t25 - 0.1154653670707977E1;
      double t30 = t29 * t27 * t26;
      double t33 = 0.9472135954999579 + t7 + t1;
      double t34 = t33 * t6;
      double t38 = t17 * t6;
      double t44 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t45 = -t3 - t25 - 0.5278640450004206E-1;
      double t47 = -t3 - t25 - 0.9472135954999579;
      double t48 = t47 * t45 * t44;
      double t53 = t8 * t33;
      double t54 = t12 * t53;
      double t62 = t44 * t6;
      double t72 = t26 * t6;
      double t76 = t45 * t6;
      double t89 =
        -0.5094133317734938E2 * t12 * t9 * t5 + 0.9612825967932449E2 * t21 * t16 -
        0.1601358118995934E2 * t30 * t24 - 0.3221812797718524E3 * t8 * t34 * t5 +
        0.9612825967932449E2 * t18 * t38 * t5 + 0.2758247241324992E2 * t48 * t24 +
        0.1601358118995934E2 * t30 * t16 - 0.5094133317734938E2 * t54 * t16 -
        0.5094133317734938E2 * t12 * t34 * t5 - 0.9556308104802193E1 * t12 * t16 +
        0.5516494482649984E2 * t47 * t62 * t5 + 0.9556308104802193E1 * t12 * t24 +
        0.1374753097369824E2 * t44 * t24 - 0.2758247241324992E2 * t48 * t16 -
        0.3202716237991868E2 * t27 * t72 * t5 +
        0.8722343632463083E2 * t47 * t76 * t5 +
        0.5516494482649984E2 * t45 * t62 * t5 +
        0.9612825967932449E2 * t20 * t18 * t6 * t5 +
        0.5094133317734938E2 * t54 * t24;
      double t103 = t47 * t45;
      double t115 = t12 * t44;
      double t136 =
        -0.1374753097369824E2 * t44 * t16 - 0.1696588918640519E2 * t6 * t5 +
        0.5261910701449129E1 * t24 - 0.9612825967932449E2 * t21 * t24 +
        0.9612825967932449E2 * t20 * t38 * t5 +
        0.1458234147916518E2 * t47 * t6 * t5 + 0.729117073958259E1 * t103 * t24 -
        0.2331143372090885E2 * t53 * t16 + 0.7913064998783579E1 * t12 * t6 * t5 -
        0.729117073958259E1 * t103 * t16 + 0.1458234147916518E2 * t76 * t5 +
        0.2502330866911363E1 * t115 * t24 - 0.1582612999756716E2 * t62 * t5 -
        0.2502330866911363E1 * t115 * t16 -
        0.3202716237991868E2 * t29 * t27 * t6 * t5 +
        0.2331143372090885E2 * t53 * t24 - 0.2331143372090885E2 * t34 * t5 -
        0.5261910701449129E1 * t16 - 0.3202716237991868E2 * t29 * t72 * t5 -
        0.2331143372090885E2 * t9 * t5;
      double t137 = t89 + t136;
      return t137;
    }

    // ORDER 7

    // Edge functions, order 7

    // number 29
    inline double
    ortho2_f29_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l2 * l3 * phi5 (l3 - l2);
    }

    inline double
    ortho2_f29_1 (double x, double y)
    {
      return -ortho2_f29_0(x, y);
    }

    inline double
    ortho2_f29x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l2x * l3 + l2 * l3x) * phi5 (l3 - l2) + l2 * l3 * phi5x (l3 -
                       l2) *
        (l3x - l2x);
    }

    inline double
    ortho2_f29x_1 (double x, double y)
    {
      return -ortho2_f29x_0(x, y);
    }

    inline double
    ortho2_f29y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l2y * l3 + l2 * l3y) * phi5 (l3 - l2) + l2 * l3 * phi5x (l3 -
                       l2) *
        (l3y - l2y);
    }

    inline double
    ortho2_f29y_1 (double x, double y)
    {
      return -ortho2_f29y_0(x, y);
    }

    // number 30
    inline double
    ortho2_f30_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l3 * l1 * phi5 (l1 - l3);
    }

    inline double
    ortho2_f30_1 (double x, double y)
    {
      return -ortho2_f30_0(x, y);
    }

    inline double
    ortho2_f30x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l3x * l1 + l3 * l1x) * phi5 (l1 - l3) + l3 * l1 * phi5x (l1 -
                       l3) *
        (l1x - l3x);
    }

    inline double
    ortho2_f30x_1 (double x, double y)
    {
      return -ortho2_f30x_0(x, y);
    }

    inline double
    ortho2_f30y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l3y * l1 + l3 * l1y) * phi5 (l1 - l3) + l3 * l1 * phi5x (l1 -
                       l3) *
        (l1y - l3y);
    }

    inline double
    ortho2_f30y_1 (double x, double y)
    {
      return -ortho2_f30y_0(x, y);
    }

    // number 31
    inline double
    ortho2_f31_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l1 * l2 * phi5 (l2 - l1);
    }

    inline double
    ortho2_f31_1 (double x, double y)
    {
      return -ortho2_f31_0(x, y);
    }

    inline double
    ortho2_f31x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l1x * l2 + l1 * l2x) * phi5 (l2 - l1) + l1 * l2 * phi5x (l2 -
                       l1) *
        (l2x - l1x);
    }

    inline double
    ortho2_f31x_1 (double x, double y)
    {
      return -ortho2_f31x_0(x, y);
    }

    inline double
    ortho2_f31y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l1y * l2 + l1 * l2y) * phi5 (l2 - l1) + l1 * l2 * phi5x (l2 -
                       l1) *
        (l2y - l1y);
    }

    inline double
    ortho2_f31y_1 (double x, double y)
    {
      return -ortho2_f31y_0(x, y);
    }

    // Bubble functions, order 7

    // * f32 *********************************************************************

    double
    ortho2_f32 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t20 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t21 = -t3 - t8 - 0.5278640450004206E-1;
      double t23 = -t3 - t8 - 0.9472135954999579;
      double t27 = 0.1E1 * x;
      double t28 = 0.9472135954999579 + t27 + t1;
      double t30 = 0.5278640450004206E-1 + t27 + t1;
      double t38 = t20 * t6;
      double t41 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t69 =
        0.2362831568572143E3 * (-t3 - t8 - 0.1265055323929465E1) * (-t3 - t8 -
                    0.7852315164806451)
        * (-t3 - t8 - 0.2147684835193549) * (-t3 - t8 + 0.2650553239294647) * t7 -
        0.399012080959134E2 * t23 * t21 * t20 * t7 +
        0.6556143495302777E1 * t30 * t28 * t6 * t5 +
        0.8720207566480212E2 * t23 * t21 * t6 * t5 +
        0.1044871155641011E1 * t41 * t38 * t5 -
        0.8582914604564961E1 * t41 * t6 * t5 - 0.5289268356110648E1 * t38 * t5 +
        0.4340933853627001 * t41 * t30 * t28 * t7 -
        0.7759664772434117 * (-0.1546536707079771 + t27 + t1) * (0.5 + t27 +
                       t1) *
        (0.1154653670707977E1 + t27 + t1) * t7 + 0.8196043575638882E2 * (-t3 -
                         t8 -
                         0.1154653670707977E1)
        * (-t3 - t8 - 0.5) * (-t3 - t8 + 0.1546536707079771) * t7 +
        0.1001403753788686E2 * t7;
      return t69;
    }

    double
    ortho2_f32x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * y;
      double t7 = -t3 - t6 - 0.5278640450004206E-1;
      double t8 = -t3 - t6 - 0.9472135954999579;
      double t9 = t8 * t7;
      double t13 = (-t3 - t1) * t2;
      double t16 = t4 * t13;
      double t20 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t24 = 0.1E1 * x;
      double t25 = 0.9472135954999579 + t24 + t1;
      double t26 = t25 * t4;
      double t27 = 0.5278640450004206E-1 + t24 + t1;
      double t31 = t7 * t4;
      double t36 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t48 =
        -0.4360103783240106E2 * t9 * t5 + 0.4360103783240106E2 * t9 * t13 +
        0.3155355604589882E1 * t16 - 0.3304172713237888E1 * t20 * t4 * t13 +
        0.6863619074796649 * t27 * t26 * t13 - 0.4360103783240106E2 * t31 * t13 +
        0.2644634178055324E1 * t36 * t5 - 0.4360103783240106E2 * t8 * t4 * t13 -
        0.429145730228248E1 * t20 * t13 + 0.429145730228248E1 * t20 * t5 -
        0.2644634178055324E1 * t36 * t13;
      double t49 = t20 * t36;
      double t54 = t36 * t4;
      double t57 = t27 * t4;
      double t62 = t27 * t25;
      double t70 = t20 * t62;
      double t78 =
        0.5224355778205057 * t49 * t13 - 0.5224355778205057 * t49 * t5 +
        0.1652086356618944E1 * t54 * t13 + 0.6556143495302777E1 * t57 * t13 +
        0.6556143495302777E1 * t26 * t13 - 0.3278071747651388E1 * t62 * t5 +
        0.3278071747651388E1 * t62 * t13 + 0.1261786989754366E3 * t8 * t31 * t13 -
        0.21704669268135 * t70 * t5 + 0.21704669268135 * t70 * t13 +
        0.4340933853627001 * t20 * t57 * t13;
      double t83 = -t3 - t6 + 0.2650553239294647;
      double t85 = -t3 - t6 - 0.2147684835193549;
      double t86 = -t3 - t6 - 0.7852315164806451;
      double t88 = -t3 - t6 - 0.1265055323929465E1;
      double t89 = t88 * t86 * t85;
      double t94 = 0.1154653670707977E1 + t24 + t1;
      double t95 = 0.5 + t24 + t1;
      double t97 = -0.1546536707079771 + t24 + t1;
      double t98 = t97 * t95 * t94;
      double t110 = t85 * t83;
      double t117 = -t3 - t6 - 0.5;
      double t119 = -t3 - t6 - 0.1154653670707977E1;
      double t123 =
        0.4340933853627001 * t20 * t26 * t13 -
        0.1181415784286071E3 * t89 * t83 * t5 + 0.5007018768943429E1 * t13 -
        0.5007018768943429E1 * t5 + 0.3879832386217059 * t98 * t5 -
        0.1181415784286071E3 * t89 * t16 +
        0.1181415784286071E3 * t89 * t83 * t13 -
        0.1181415784286071E3 * t88 * t86 * t83 * t16 -
        0.1181415784286071E3 * t88 * t110 * t16 -
        0.1181415784286071E3 * t86 * t110 * t16 -
        0.4098021787819441E2 * t119 * t117 * t4 * t13;
      double t124 = -t3 - t6 + 0.1546536707079771;
      double t125 = t124 * t4;
      double t130 = t119 * t117 * t124;
      double t140 = t94 * t4;
      double t152 = t8 * t7 * t36;
      double t163 =
        -0.4098021787819441E2 * t119 * t125 * t13 -
        0.4098021787819441E2 * t130 * t5 + 0.4098021787819441E2 * t130 * t13 -
        0.4098021787819441E2 * t117 * t125 * t13 -
        0.3879832386217059 * t98 * t13 - 0.7759664772434117 * t97 * t140 * t13 -
        0.7759664772434117 * t97 * t95 * t4 * t13 -
        0.7759664772434117 * t95 * t140 * t13 + 0.199506040479567E2 * t152 * t5 -
        0.199506040479567E2 * t152 * t13 + 0.199506040479567E2 * t8 * t54 * t13 +
        0.199506040479567E2 * t7 * t54 * t13;
      double t165 = t48 + t78 + t123 + t163;
      return t165;
    }

    double
    ortho2_f32y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * x;
      double t8 = 0.5278640450004206E-1 + t7 + t1;
      double t9 = t8 * t6;
      double t12 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t16 = 0.9472135954999579 + t7 + t1;
      double t17 = t16 * t6;
      double t21 = t6 * t2;
      double t24 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t27 = t6 * t4;
      double t30 = 0.1E1 * y;
      double t31 = -t3 - t30 + 0.2650553239294647;
      double t33 = -t3 - t30 - 0.2147684835193549;
      double t34 = -t3 - t30 - 0.7852315164806451;
      double t36 = -t3 - t30 - 0.1265055323929465E1;
      double t37 = t36 * t34 * t33;
      double t42 = 0.1154653670707977E1 + t7 + t1;
      double t43 = 0.5 + t7 + t1;
      double t45 = -0.1546536707079771 + t7 + t1;
      double t46 = t45 * t43 * t42;
      double t51 = t6 * t5;
      double t58 =
        0.21704669268135 * t12 * t9 * t5 + 0.21704669268135 * t12 * t17 * t5 +
        0.2644634178055324E1 * t24 * t21 + 0.5007018768943429E1 * t27 -
        0.5007018768943429E1 * t21 - 0.1181415784286071E3 * t37 * t31 * t21 +
        0.429145730228248E1 * t12 * t21 + 0.3879832386217059 * t46 * t21 -
        0.2644634178055324E1 * t24 * t27 - 0.2362831568572143E3 * t37 * t51 -
        0.2362831568572143E3 * t36 * t34 * t31 * t51;
      double t59 = t33 * t31;
      double t66 = -t3 - t30 - 0.5;
      double t68 = -t3 - t30 - 0.1154653670707977E1;
      double t72 = -t3 - t30 + 0.1546536707079771;
      double t73 = t72 * t6;
      double t79 = t68 * t66 * t72;
      double t84 = t8 * t16;
      double t90 = -t3 - t30 - 0.5278640450004206E-1;
      double t91 = t90 * t6;
      double t94 = -t3 - t30 - 0.9472135954999579;
      double t98 =
        -0.2362831568572143E3 * t36 * t59 * t51 -
        0.2362831568572143E3 * t34 * t59 * t51 -
        0.8196043575638882E2 * t68 * t66 * t6 * t5 -
        0.8196043575638882E2 * t68 * t73 * t5 - 0.1877849153256658E2 * t51 -
        0.4098021787819441E2 * t79 * t21 - 0.429145730228248E1 * t12 * t27 +
        0.3278071747651388E1 * t84 * t27 - 0.8196043575638882E2 * t66 * t73 * t5 -
        0.8720207566480212E2 * t91 * t5 - 0.8720207566480212E2 * t94 * t6 * t5;
      double t100 = t94 * t90;
      double t105 = t12 * t24;
      double t108 = t24 * t6;
      double t126 = t42 * t6;
      double t130 =
        -0.4360103783240106E2 * t100 * t21 + 0.4360103783240106E2 * t100 * t27 -
        0.5224355778205057 * t105 * t21 + 0.3304172713237888E1 * t108 * t5 -
        0.1652086356618944E1 * t12 * t6 * t5 + 0.3278071747651388E1 * t17 * t5 -
        0.3278071747651388E1 * t84 * t21 + 0.5224355778205057 * t105 * t27 +
        0.3278071747651388E1 * t9 * t5 -
        0.3879832386217059 * t45 * t43 * t6 * t5 -
        0.3879832386217059 * t45 * t126 * t5;
      double t134 = t12 * t84;
      double t145 = t94 * t90 * t24;
      double t164 =
        -0.3879832386217059 * t43 * t126 * t5 + 0.21704669268135 * t134 * t27 +
        0.4098021787819441E2 * t79 * t27 +
        0.1181415784286071E3 * t37 * t31 * t27 - 0.3879832386217059 * t46 * t27 +
        0.199506040479567E2 * t145 * t21 + 0.399012080959134E2 * t94 * t108 * t5 +
        0.137272381495933E1 * t8 * t17 * t5 +
        0.630893494877183E2 * t94 * t91 * t5 +
        0.399012080959134E2 * t90 * t108 * t5 - 0.199506040479567E2 * t145 * t27 -
        0.21704669268135 * t134 * t21;
      double t166 = t58 + t98 + t130 + t164;
      return t166;
    }

    // * f33 *********************************************************************

    double
    ortho2_f33 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t20 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t21 = -t3 - t8 - 0.5278640450004206E-1;
      double t23 = -t3 - t8 - 0.9472135954999579;
      double t27 = 0.1E1 * x;
      double t28 = 0.9472135954999579 + t27 + t1;
      double t30 = 0.5278640450004206E-1 + t27 + t1;
      double t38 = t20 * t6;
      double t41 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t61 = -t3 - t8 + 0.1546536707079771;
      double t62 = -t3 - t8 - 0.5;
      double t64 = -t3 - t8 - 0.1154653670707977E1;
      double t74 =
        0.3233244967814759E3 * (-t3 - t8 - 0.1265055323929465E1) * (-t3 - t8 -
                    0.7852315164806451)
        * (-t3 - t8 - 0.2147684835193549) * (-t3 - t8 + 0.2650553239294647) * t7 -
        0.6678724813609661E2 * t23 * t21 * t20 * t7 +
        0.3334698673272348E1 * t30 * t28 * t6 * t5 +
        0.1684811178666115E3 * t23 * t21 * t6 * t5 +
        0.1759563970786134E2 * t41 * t38 * t5 -
        0.2724006919222775E2 * t41 * t6 * t5 - 0.1222514087066895E2 * t38 * t5 -
        0.5057000636728777E2 * t41 * t30 * t28 * t7 -
        0.1015537544238751E1 * (-0.1546536707079771 + t27 + t1) * (0.5 + t27 +
                         t1) *
        (0.1154653670707977E1 + t27 + t1) * t7 +
        0.173747923689801E3 * t64 * t62 * t61 * t7 -
        0.1900540576406206E3 * t64 * t62 * t61 * t20 * t7 +
        0.1921620998468437E2 * t7;
      return t74;
    }

    double
    ortho2_f33x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = 0.1E1 * y;
      double t12 = -t3 - t11 - 0.5278640450004206E-1;
      double t16 = t6 * t2;
      double t17 = t9 * t16;
      double t18 = -t3 - t11 + 0.1546536707079771;
      double t19 = -t3 - t11 - 0.5;
      double t21 = -t3 - t11 - 0.1154653670707977E1;
      double t22 = t21 * t19 * t18;
      double t25 = t9 * t5;
      double t28 = t6 * t5;
      double t33 = t18 * t9;
      double t40 = -t3 - t11 - 0.9472135954999579;
      double t41 = t40 * t12;
      double t49 = 0.1E1 * x;
      double t50 = 0.9472135954999579 + t49 + t1;
      double t51 = t50 * t6;
      double t52 = 0.5278640450004206E-1 + t49 + t1;
      double t56 = t12 * t6;
      double t61 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t62 = t61 * t9;
      double t65 =
        0.3339362406804831E2 * t12 * t10 * t5 + 0.9502702882031028E2 * t22 * t17 -
        0.9502702882031028E2 * t22 * t25 +
        0.9502702882031028E2 * t21 * t19 * t9 * t28 +
        0.9502702882031028E2 * t21 * t33 * t28 +
        0.9502702882031028E2 * t19 * t33 * t28 + 0.8424055893330573E2 * t41 * t5 -
        0.8424055893330573E2 * t41 * t16 - 0.8424055893330573E2 * t40 * t6 * t5 -
        0.799582007049234E2 * t52 * t51 * t5 - 0.8424055893330573E2 * t56 * t5 +
        0.8797819853930669E1 * t62 * t5;
      double t73 = t52 * t50;
      double t76 = t52 * t6;
      double t86 = t61 * t73;
      double t99 =
        -0.8797819853930669E1 * t62 * t16 + 0.2782114918227079E2 * t10 * t5 -
        0.5564229836454158E2 * t61 * t6 * t5 + 0.1667349336636174E1 * t73 * t5 +
        0.3334698673272348E1 * t76 * t5 + 0.3334698673272348E1 * t51 * t5 -
        0.1667349336636174E1 * t73 * t16 + 0.2111998227649005E3 * t40 * t56 * t5 +
        0.2528500318364388E2 * t86 * t16 - 0.2528500318364388E2 * t86 * t5 +
        0.6010037007012879E3 * t22 * t28 - 0.5057000636728777E2 * t61 * t76 * t5 -
        0.5057000636728777E2 * t61 * t51 * t5;
      double t102 = -t3 - t11 + 0.2650553239294647;
      double t104 = -t3 - t11 - 0.2147684835193549;
      double t105 = -t3 - t11 - 0.7852315164806451;
      double t107 = -t3 - t11 - 0.1265055323929465E1;
      double t108 = t107 * t105 * t104;
      double t117 = 0.1154653670707977E1 + t49 + t1;
      double t118 = 0.5 + t49 + t1;
      double t120 = -0.1546536707079771 + t49 + t1;
      double t121 = t120 * t118 * t117;
      double t133 = t104 * t102;
      double t144 =
        0.6112570435334473E1 * t17 - 0.161662248390738E3 * t108 * t102 * t16 -
        0.6112570435334473E1 * t25 - 0.1362003459611388E2 * t61 * t5 +
        0.1362003459611388E2 * t61 * t16 - 0.4411041266283546E1 * t28 +
        0.5077687721193755 * t121 * t16 + 0.161662248390738E3 * t108 * t102 * t5 -
        0.161662248390738E3 * t108 * t28 -
        0.161662248390738E3 * t107 * t105 * t102 * t28 -
        0.161662248390738E3 * t107 * t133 * t28 -
        0.161662248390738E3 * t105 * t133 * t28 -
        0.8687396184490052E2 * t21 * t19 * t6 * t5;
      double t145 = t18 * t6;
      double t162 = t117 * t6;
      double t172 = t40 * t12 * t9;
      double t180 =
        -0.8687396184490052E2 * t21 * t145 * t5 -
        0.8687396184490052E2 * t22 * t16 + 0.8687396184490052E2 * t22 * t5 -
        0.8687396184490052E2 * t19 * t145 * t5 - 0.5077687721193755 * t121 * t5 -
        0.1015537544238751E1 * t120 * t118 * t6 * t5 -
        0.1015537544238751E1 * t120 * t162 * t5 + 0.9608104992342185E1 * t5 -
        0.1015537544238751E1 * t118 * t162 * t5 - 0.9608104992342185E1 * t16 +
        0.3339362406804831E2 * t172 * t16 - 0.3339362406804831E2 * t172 * t5 +
        0.3339362406804831E2 * t40 * t10 * t5;
      double t182 = t65 + t99 + t144 + t180;
      return t182;
    }

    double
    ortho2_f33y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * y;
      double t8 = -t3 - t7 + 0.1546536707079771;
      double t9 = t8 * t6;
      double t10 = -t3 - t7 - 0.5;
      double t14 = t6 * t5;
      double t17 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t18 = t8 * t17;
      double t22 = 0.1E1 * x;
      double t23 = 0.5 + t22 + t1;
      double t25 = -0.1546536707079771 + t22 + t1;
      double t29 = t6 * t2;
      double t32 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t35 = t6 * t4;
      double t36 = t17 * t35;
      double t41 = t17 * t6;
      double t42 = -t3 - t7 - 0.5278640450004206E-1;
      double t47 = -t3 - t7 - 0.1154653670707977E1;
      double t48 = t47 * t10 * t8;
      double t54 = 0.5278640450004206E-1 + t22 + t1;
      double t55 = t54 * t6;
      double t61 =
        -0.173747923689801E3 * t10 * t9 * t5 +
        0.1900540576406206E3 * t10 * t18 * t14 -
        0.5077687721193755 * t25 * t23 * t6 * t5 +
        0.1362003459611388E2 * t32 * t29 - 0.6112570435334473E1 * t36 -
        0.6681101733415882E2 * t14 - 0.1362003459611388E2 * t32 * t35 +
        0.6678724813609661E2 * t42 * t41 * t5 - 0.8687396184490052E2 * t48 * t29 -
        0.173747923689801E3 * t47 * t9 * t5 -
        0.2528500318364388E2 * t32 * t55 * t5 - 0.9502702882031028E2 * t48 * t36;
      double t62 = t42 * t6;
      double t63 = -t3 - t7 - 0.9472135954999579;
      double t68 = 0.9472135954999579 + t22 + t1;
      double t69 = t68 * t6;
      double t73 = -t3 - t7 + 0.2650553239294647;
      double t74 = -t3 - t7 - 0.2147684835193549;
      double t75 = t74 * t73;
      double t76 = -t3 - t7 - 0.1265055323929465E1;
      double t80 = 0.1154653670707977E1 + t22 + t1;
      double t81 = t80 * t6;
      double t85 = t54 * t68;
      double t86 = t32 * t85;
      double t94 = t17 * t29;
      double t97 = -t3 - t7 - 0.7852315164806451;
      double t99 = t76 * t97 * t74;
      double t107 = t25 * t23 * t80;
      double t110 =
        0.1055999113824503E3 * t63 * t62 * t5 - 0.9608104992342185E1 * t29 -
        0.2528500318364388E2 * t32 * t69 * t5 -
        0.3233244967814759E3 * t76 * t75 * t14 -
        0.5077687721193755 * t25 * t81 * t5 + 0.2528500318364388E2 * t86 * t29 +
        0.300501850350644E3 * t48 * t14 - 0.5077687721193755 * t23 * t81 * t5 +
        0.6112570435334473E1 * t94 - 0.161662248390738E3 * t99 * t73 * t29 +
        0.9608104992342185E1 * t35 + 0.161662248390738E3 * t99 * t73 * t35 -
        0.5077687721193755 * t107 * t35;
      double t113 = t63 * t42 * t17;
      double t149 =
        -0.3339362406804831E2 * t113 * t35 + 0.8687396184490052E2 * t48 * t35 +
        0.9502702882031028E2 * t48 * t94 -
        0.3233244967814759E3 * t97 * t75 * t14 -
        0.2528500318364388E2 * t86 * t35 -
        0.173747923689801E3 * t47 * t10 * t6 * t5 +
        0.5077687721193755 * t107 * t29 - 0.3233244967814759E3 * t99 * t14 -
        0.3233244967814759E3 * t76 * t97 * t73 * t14 +
        0.1900540576406206E3 * t47 * t10 * t17 * t14 +
        0.1900540576406206E3 * t47 * t18 * t14 +
        0.1667349336636174E1 * t85 * t35 - 0.2782114918227079E2 * t32 * t6 * t5;
      double t150 = t63 * t42;
      double t160 = t32 * t17;
      double t181 =
        0.8424055893330573E2 * t150 * t35 - 0.1684811178666115E3 * t62 * t5 -
        0.1684811178666115E3 * t63 * t6 * t5 - 0.8424055893330573E2 * t150 * t29 -
        0.8797819853930669E1 * t160 * t29 + 0.5564229836454158E2 * t41 * t5 +
        0.8797819853930669E1 * t160 * t35 + 0.3339362406804831E2 * t113 * t29 +
        0.1667349336636174E1 * t55 * t5 + 0.1667349336636174E1 * t69 * t5 -
        0.1667349336636174E1 * t85 * t29 - 0.1599164014098468E3 * t54 * t69 * t5 +
        0.6678724813609661E2 * t63 * t41 * t5;
      double t183 = t61 + t110 + t149 + t181;
      return t183;
    }

    // * f34 *********************************************************************

    double
    ortho2_f34 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t20 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t21 = -t3 - t8 - 0.5278640450004206E-1;
      double t23 = -t3 - t8 - 0.9472135954999579;
      double t27 = 0.1E1 * x;
      double t28 = 0.9472135954999579 + t27 + t1;
      double t30 = 0.5278640450004206E-1 + t27 + t1;
      double t38 = t20 * t6;
      double t41 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t45 = t30 * t28;
      double t65 = -t3 - t8 + 0.1546536707079771;
      double t66 = -t3 - t8 - 0.5;
      double t68 = -t3 - t8 - 0.1154653670707977E1;
      double t78 =
        0.2511246756682848E3 * (-t3 - t8 - 0.1265055323929465E1) * (-t3 - t8 -
                    0.7852315164806451)
        * (-t3 - t8 - 0.2147684835193549) * (-t3 - t8 + 0.2650553239294647) * t7 -
        0.831812204547802E2 * t23 * t21 * t20 * t7 +
        0.1176036009327637E3 * t30 * t28 * t6 * t5 +
        0.2880344953802456E3 * t23 * t21 * t6 * t5 +
        0.329692714047373E2 * t41 * t38 * t5 +
        0.8561461356421576E3 * t23 * t21 * t45 * t7 -
        0.2260724429175527E2 * t41 * t6 * t5 - 0.2380663034055424E2 * t38 * t5 -
        0.567181428190866E2 * t41 * t45 * t7 +
        0.7856907955363914E2 * (-0.1546536707079771 + t27 + t1) * (0.5 + t27 +
                         t1) *
        (0.1154653670707977E1 + t27 + t1) * t7 +
        0.103024103259328E3 * t68 * t66 * t65 * t7 -
        0.293859669769857E3 * t68 * t66 * t65 * t20 * t7 +
        0.3753527718161051E2 * t7;
      return t78;
    }

    double
    ortho2_f34x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t8 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t9 = 0.1E1 * y;
      double t10 = -t3 - t9 - 0.5278640450004206E-1;
      double t12 = -t3 - t9 - 0.9472135954999579;
      double t13 = t12 * t10 * t8;
      double t17 = (-t3 - t1) * t2;
      double t20 = t8 * t4;
      double t27 = t8 * t5;
      double t28 = -t3 - t9 + 0.1546536707079771;
      double t29 = -t3 - t9 - 0.5;
      double t31 = -t3 - t9 - 0.1154653670707977E1;
      double t32 = t31 * t29 * t28;
      double t35 = t8 * t17;
      double t38 = t4 * t17;
      double t43 = t28 * t8;
      double t50 = 0.1E1 * x;
      double t51 = 0.9472135954999579 + t50 + t1;
      double t53 = 0.5278640450004206E-1 + t50 + t1;
      double t55 = t12 * t10 * t53;
      double t67 = t53 * t51;
      double t71 =
        0.415906102273901E2 * t13 * t5 - 0.415906102273901E2 * t13 * t17 +
        0.415906102273901E2 * t12 * t20 * t17 +
        0.415906102273901E2 * t10 * t20 * t17 + 0.1469298348849285E3 * t32 * t27 -
        0.1469298348849285E3 * t32 * t35 +
        0.1469298348849285E3 * t31 * t29 * t8 * t38 +
        0.1469298348849285E3 * t31 * t43 * t38 +
        0.1469298348849285E3 * t29 * t43 * t38 -
        0.4280730678210788E3 * t55 * t51 * t5 +
        0.4280730678210788E3 * t55 * t51 * t17 +
        0.8561461356421576E3 * t55 * t38 +
        0.8561461356421576E3 * t12 * t10 * t51 * t38 -
        0.4280730678210788E3 * t12 * t67 * t38;
      double t75 = t12 * t10;
      double t82 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t83 = t82 * t8;
      double t93 = t10 * t4;
      double t101 = t53 * t4;
      double t104 = t51 * t4;
      double t110 =
        -0.4280730678210788E3 * t10 * t67 * t38 +
        0.1440172476901228E3 * t75 * t17 - 0.1440172476901228E3 * t75 * t5 +
        0.1648463570236865E2 * t83 * t17 - 0.1648463570236865E2 * t83 * t5 +
        0.5212899521761446E2 * t20 * t17 - 0.1042579904352289E3 * t82 * t4 * t17 -
        0.1440172476901228E3 * t93 * t17 - 0.1440172476901228E3 * t12 * t4 * t17 +
        0.5880180046638183E2 * t67 * t17 + 0.1176036009327637E3 * t101 * t17 +
        0.1176036009327637E3 * t104 * t17 - 0.5880180046638183E2 * t67 * t5 -
        0.1876763859080525E2 * t5;
      double t126 = t82 * t67;
      double t139 = -t3 - t9 + 0.2650553239294647;
      double t141 = -t3 - t9 - 0.2147684835193549;
      double t142 = -t3 - t9 - 0.7852315164806451;
      double t144 = -t3 - t9 - 0.1265055323929465E1;
      double t145 = t144 * t142 * t141;
      double t148 = 0.1876763859080525E2 * t17 - 0.1190331517027712E2 * t35 +
        0.1190331517027712E2 * t27 - 0.8967925798151857E2 * t53 * t104 * t17 -
        0.1130362214587764E2 * t82 * t17 + 0.1130362214587764E2 * t82 * t5 +
        0.2630421151896925E3 * t12 * t93 * t17 + 0.39537983498928E2 * t38 +
        0.283590714095433E2 * t126 * t5 - 0.283590714095433E2 * t126 * t17 +
        0.9292658689376761E3 * t32 * t38 -
        0.567181428190866E2 * t82 * t101 * t17 -
        0.567181428190866E2 * t82 * t104 * t17 -
        0.1255623378341424E3 * t145 * t139 * t5;
      double t149 = 0.1154653670707977E1 + t50 + t1;
      double t150 = 0.5 + t50 + t1;
      double t152 = -0.1546536707079771 + t50 + t1;
      double t153 = t152 * t150 * t149;
      double t165 = t141 * t139;
      double t176 = t28 * t4;
      double t193 = t149 * t4;
      double t200 =
        -0.3928453977681957E2 * t153 * t5 +
        0.1255623378341424E3 * t145 * t139 * t17 -
        0.1255623378341424E3 * t145 * t38 -
        0.1255623378341424E3 * t144 * t142 * t139 * t38 -
        0.1255623378341424E3 * t144 * t165 * t38 -
        0.1255623378341424E3 * t142 * t165 * t38 -
        0.5151205162966402E2 * t31 * t29 * t4 * t17 -
        0.5151205162966402E2 * t31 * t176 * t17 -
        0.5151205162966402E2 * t32 * t5 + 0.5151205162966402E2 * t32 * t17 -
        0.5151205162966402E2 * t29 * t176 * t17 +
        0.3928453977681957E2 * t153 * t17 +
        0.7856907955363914E2 * t152 * t150 * t4 * t17 +
        0.7856907955363914E2 * t152 * t193 * t17 +
        0.7856907955363914E2 * t150 * t193 * t17;
      double t202 = t71 + t110 + t148 + t200;
      return t202;
    }

    double
    ortho2_f34y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.9472135954999579 + t8 + t1;
      double t10 = 0.5278640450004206E-1 + t8 + t1;
      double t11 = t10 * t9;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 - 0.9472135954999579;
      double t17 = t6 * t2;
      double t18 = 0.1154653670707977E1 + t8 + t1;
      double t19 = 0.5 + t8 + t1;
      double t21 = -0.1546536707079771 + t8 + t1;
      double t22 = t21 * t19 * t18;
      double t25 = t10 * t6;
      double t28 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t32 = -t3 - t12 + 0.2650553239294647;
      double t33 = -t3 - t12 - 0.7852315164806451;
      double t35 = -t3 - t12 - 0.1265055323929465E1;
      double t39 = t6 * t4;
      double t42 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t43 = t42 * t39;
      double t45 = t9 * t6;
      double t52 = -t3 - t12 + 0.1546536707079771;
      double t53 = -t3 - t12 - 0.5;
      double t55 = -t3 - t12 - 0.1154653670707977E1;
      double t56 = t55 * t53 * t52;
      double t59 = t28 * t42;
      double t75 =
        -0.8561461356421576E3 * t13 * t11 * t7 -
        0.3928453977681957E2 * t22 * t17 - 0.283590714095433E2 * t28 * t25 * t5 -
        0.2511246756682848E3 * t35 * t33 * t32 * t7 - 0.1190331517027712E2 * t43 -
        0.1793585159630371E3 * t10 * t45 * t5 - 0.1130362214587764E2 * t28 * t39 -
        0.3384879593687611E2 * t7 + 0.5151205162966402E2 * t56 * t39 +
        0.1648463570236865E2 * t59 * t39 +
        0.293859669769857E3 * t55 * t53 * t42 * t7 -
        0.283590714095433E2 * t28 * t45 * t5 -
        0.103024103259328E3 * t55 * t53 * t6 * t5 +
        0.3928453977681957E2 * t22 * t39;
      double t77 = -t3 - t12 - 0.2147684835193549;
      double t79 = t35 * t33 * t77;
      double t82 = -t3 - t12 - 0.5278640450004206E-1;
      double t84 = t13 * t82 * t42;
      double t87 = t18 * t6;
      double t91 = t52 * t42;
      double t95 = t77 * t32;
      double t102 = t52 * t6;
      double t111 = t13 * t82 * t10;
      double t119 = t28 * t11;
      double t124 =
        0.1255623378341424E3 * t79 * t32 * t39 - 0.415906102273901E2 * t84 * t39 +
        0.3928453977681957E2 * t19 * t87 * t5 +
        0.293859669769857E3 * t55 * t91 * t7 -
        0.2511246756682848E3 * t35 * t95 * t7 -
        0.8561461356421576E3 * t82 * t11 * t7 -
        0.103024103259328E3 * t55 * t102 * t5 -
        0.103024103259328E3 * t53 * t102 * t5 +
        0.4280730678210788E3 * t111 * t9 * t39 + 0.1876763859080525E2 * t39 -
        0.1876763859080525E2 * t17 + 0.293859669769857E3 * t53 * t91 * t7 -
        0.283590714095433E2 * t119 * t39 + 0.1130362214587764E2 * t28 * t17;
      double t134 = t42 * t6;
      double t153 = t42 * t17;
      double t162 =
        0.415906102273901E2 * t84 * t17 + 0.283590714095433E2 * t119 * t17 -
        0.5151205162966402E2 * t56 * t17 - 0.1469298348849285E3 * t56 * t43 +
        0.831812204547802E2 * t13 * t134 * t5 -
        0.2511246756682848E3 * t33 * t95 * t7 -
        0.4280730678210788E3 * t111 * t9 * t17 +
        0.831812204547802E2 * t82 * t134 * t5 - 0.5880180046638183E2 * t11 * t17 +
        0.5880180046638183E2 * t45 * t5 + 0.5880180046638183E2 * t25 * t5 +
        0.1190331517027712E2 * t153 + 0.3928453977681957E2 * t21 * t19 * t6 * t5 -
        0.1255623378341424E3 * t79 * t32 * t17;
      double t165 = t82 * t6;
      double t181 = t13 * t82;
      double t201 =
        0.4280730678210788E3 * t111 * t7 +
        0.1315210575948462E3 * t13 * t165 * t5 + 0.464632934468838E3 * t56 * t7 +
        0.3928453977681957E2 * t21 * t87 * t5 + 0.5880180046638183E2 * t11 * t39 -
        0.2511246756682848E3 * t79 * t7 - 0.2880344953802456E3 * t13 * t6 * t5 -
        0.1440172476901228E3 * t181 * t17 + 0.1469298348849285E3 * t56 * t153 -
        0.2880344953802456E3 * t165 * t5 + 0.1440172476901228E3 * t181 * t39 -
        0.1648463570236865E2 * t59 * t17 + 0.1042579904352289E3 * t134 * t5 -
        0.5212899521761446E2 * t28 * t6 * t5 +
        0.4280730678210788E3 * t13 * t82 * t9 * t7;
      double t203 = t75 + t124 + t162 + t201;
      return t203;
    }

    // * f35 *********************************************************************

    double
    ortho2_f35 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t20 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t21 = -t3 - t8 - 0.5278640450004206E-1;
      double t23 = -t3 - t8 - 0.9472135954999579;
      double t27 = 0.1E1 * x;
      double t30 = (0.5 + t27 + t1) * (0.1154653670707977E1 + t27 + t1);
      double t31 = -0.1546536707079771 + t27 + t1;
      double t34 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t39 = 0.9472135954999579 + t27 + t1;
      double t41 = 0.5278640450004206E-1 + t27 + t1;
      double t49 = t20 * t6;
      double t53 = t41 * t39;
      double t69 = -t3 - t8 + 0.1546536707079771;
      double t70 = -t3 - t8 - 0.5;
      double t72 = -t3 - t8 - 0.1154653670707977E1;
      double t82 =
        0.1223215991420282E3 * (-t3 - t8 - 0.1265055323929465E1) * (-t3 - t8 -
                    0.7852315164806451)
        * (-t3 - t8 - 0.2147684835193549) * (-t3 - t8 + 0.2650553239294647) * t7 +
        0.6459440914577556E-2 * t23 * t21 * t20 * t7 -
        0.2542303516446689E3 * t34 * t31 * t30 * t7 +
        0.2458707861302519E3 * t41 * t39 * t6 * t5 +
        0.2777517978656213E3 * t23 * t21 * t6 * t5 +
        0.5141622074373635E2 * t34 * t49 * t5 +
        0.1314495925277894E4 * t23 * t21 * t53 * t7 +
        0.3721523473287692 * t34 * t6 * t5 - 0.3864494903587631E1 * t49 * t5 -
        0.9028021086561E1 * t34 * t53 * t7 +
        0.3904616157666349E2 * t31 * t30 * t7 -
        0.8119934108954854E1 * t72 * t70 * t69 * t7 -
        0.2199997587150551E3 * t72 * t70 * t69 * t20 * t7 +
        0.4714215058969815E2 * t7;
      return t82;
    }

    double
    ortho2_f35x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.9472135954999579 + t8 + t1;
      double t10 = 0.5278640450004206E-1 + t8 + t1;
      double t11 = t10 * t9;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 - 0.5278640450004206E-1;
      double t17 = t6 * t2;
      double t20 = t9 * t6;
      double t24 = t13 * t6;
      double t25 = -t3 - t12 - 0.9472135954999579;
      double t29 = t25 * t13;
      double t36 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t46 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t47 = t46 * t6;
      double t53 = t46 * t5;
      double t55 = t46 * t17;
      double t58 =
        -0.6572479626389472E3 * t13 * t11 * t7 - 0.2357107529484908E2 * t17 +
        0.2357107529484908E2 * t5 - 0.1427455469878045E2 * t10 * t20 * t5 -
        0.2042654570134621E-1 * t25 * t24 * t5 -
        0.1388758989328107E3 * t29 * t17 + 0.1388758989328107E3 * t29 * t5 +
        0.1860761736643846 * t36 * t5 - 0.1388758989328107E3 * t25 * t6 * t5 -
        0.1388758989328107E3 * t24 * t5 + 0.8129618311410174E2 * t47 * t5 -
        0.1625923662282035E3 * t36 * t6 * t5 - 0.1932247451793815E1 * t53 +
        0.1932247451793815E1 * t55 + 0.1280903042851821E2 * t7;
      double t59 = t36 * t46;
      double t70 = t10 * t6;
      double t75 = 0.1154653670707977E1 + t8 + t1;
      double t76 = 0.5 + t8 + t1;
      double t77 = t76 * t75;
      double t78 = -0.1546536707079771 + t8 + t1;
      double t79 = t78 * t77;
      double t85 = t36 * t11;
      double t90 = -t3 - t12 + 0.1546536707079771;
      double t91 = -t3 - t12 - 0.5;
      double t93 = -t3 - t12 - 0.1154653670707977E1;
      double t94 = t93 * t91 * t90;
      double t103 = -t3 - t12 + 0.2650553239294647;
      double t105 = -t3 - t12 - 0.2147684835193549;
      double t106 = -t3 - t12 - 0.7852315164806451;
      double t108 = -t3 - t12 - 0.1265055323929465E1;
      double t109 = t108 * t106 * t105;
      double t114 =
        0.2570811037186818E2 * t59 * t5 - 0.1860761736643846 * t36 * t17 -
        0.2570811037186818E2 * t59 * t17 - 0.122935393065126E3 * t11 * t17 +
        0.122935393065126E3 * t11 * t5 + 0.2458707861302519E3 * t70 * t5 +
        0.2458707861302519E3 * t20 * t5 - 0.4019734807713439E3 * t79 * t7 -
        0.2542303516446689E3 * t36 * t77 * t7 + 0.45140105432805E1 * t85 * t17 -
        0.45140105432805E1 * t85 * t5 + 0.6957003222270525E3 * t94 * t7 -
        0.9028021086561E1 * t36 * t70 * t5 - 0.9028021086561E1 * t36 * t20 * t5 -
        0.6116079957101409E2 * t109 * t103 * t17 -
        0.1952308078833175E2 * t79 * t17;
      double t117 = t36 * t78 * t76;
      double t133 = t105 * t103;
      double t144 = t90 * t6;
      double t161 = t75 * t6;
      double t168 =
        -0.2542303516446689E3 * t117 * t7 -
        0.2542303516446689E3 * t36 * t78 * t75 * t7 +
        0.6116079957101409E2 * t109 * t103 * t5 -
        0.6116079957101409E2 * t109 * t7 -
        0.6116079957101409E2 * t108 * t106 * t103 * t7 -
        0.6116079957101409E2 * t108 * t133 * t7 -
        0.6116079957101409E2 * t106 * t133 * t7 +
        0.4059967054477427E1 * t93 * t91 * t6 * t5 +
        0.4059967054477427E1 * t93 * t144 * t5 +
        0.4059967054477427E1 * t94 * t17 - 0.4059967054477427E1 * t94 * t5 +
        0.4059967054477427E1 * t91 * t144 * t5 + 0.1952308078833175E2 * t79 * t5 +
        0.3904616157666349E2 * t78 * t76 * t6 * t5 +
        0.3904616157666349E2 * t78 * t161 * t5 +
        0.3904616157666349E2 * t76 * t161 * t5;
      double t170 = t25 * t13 * t46;
      double t195 = t90 * t46;
      double t204 = t25 * t13 * t10;
      double t219 =
        -0.3229720457288778E-2 * t170 * t17 + 0.3229720457288778E-2 * t170 * t5 -
        0.3229720457288778E-2 * t25 * t47 * t5 -
        0.3229720457288778E-2 * t13 * t47 * t5 +
        0.1271151758223344E3 * t117 * t75 * t17 -
        0.1271151758223344E3 * t117 * t75 * t5 +
        0.1099998793575276E3 * t94 * t55 - 0.1099998793575276E3 * t94 * t53 +
        0.1099998793575276E3 * t93 * t91 * t46 * t7 +
        0.1099998793575276E3 * t93 * t195 * t7 +
        0.1099998793575276E3 * t91 * t195 * t7 -
        0.6572479626389472E3 * t204 * t9 * t17 +
        0.6572479626389472E3 * t204 * t9 * t5 + 0.1314495925277894E4 * t204 * t7 +
        0.1314495925277894E4 * t25 * t13 * t9 * t7 -
        0.6572479626389472E3 * t25 * t11 * t7;
      double t221 = t58 + t114 + t168 + t219;
      return t221;
    }

    double
    ortho2_f35y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t9 = -t3 - t8 + 0.2650553239294647;
      double t10 = -t3 - t8 - 0.7852315164806451;
      double t12 = -t3 - t8 - 0.1265055323929465E1;
      double t16 = 0.1E1 * x;
      double t17 = 0.5278640450004206E-1 + t16 + t1;
      double t18 = -t3 - t8 - 0.5278640450004206E-1;
      double t20 = -t3 - t8 - 0.9472135954999579;
      double t21 = t20 * t18 * t17;
      double t24 = -t3 - t8 - 0.2147684835193549;
      double t25 = t24 * t9;
      double t29 = 0.5 + t16 + t1;
      double t31 = -0.1546536707079771 + t16 + t1;
      double t35 = t6 * t4;
      double t38 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t41 = 0.1154653670707977E1 + t16 + t1;
      double t42 = t41 * t6;
      double t49 = 0.9472135954999579 + t16 + t1;
      double t55 = t17 * t49;
      double t56 = t38 * t55;
      double t59 = t6 * t2;
      double t61 = t49 * t6;
      double t70 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t71 = t70 * t59;
      double t76 =
        -0.1223215991420282E3 * t12 * t10 * t9 * t7 +
        0.6572479626389472E3 * t21 * t7 - 0.1223215991420282E3 * t12 * t25 * t7 +
        0.1952308078833175E2 * t31 * t29 * t6 * t5 +
        0.1860761736643846 * t38 * t35 + 0.1952308078833175E2 * t31 * t42 * t5 +
        0.1952308078833175E2 * t29 * t42 * t5 +
        0.6572479626389472E3 * t20 * t18 * t49 * t7 + 0.2357107529484908E2 * t35 -
        0.45140105432805E1 * t56 * t35 - 0.2357107529484908E2 * t59 -
        0.2854910939756091E2 * t17 * t61 * t5 -
        0.1314495925277894E4 * t20 * t55 * t7 + 0.1932247451793815E1 * t71 -
        0.1314495925277894E4 * t18 * t55 * t7;
      double t84 = t12 * t10 * t24;
      double t90 = -t3 - t8 - 0.5;
      double t92 = -t3 - t8 - 0.1154653670707977E1;
      double t96 = t29 * t41;
      double t97 = t31 * t96;
      double t101 = t20 * t18 * t70;
      double t104 = -t3 - t8 + 0.1546536707079771;
      double t105 = t104 * t6;
      double t109 = t70 * t6;
      double t113 = t18 * t6;
      double t132 =
        -0.1860761736643846 * t38 * t59 - 0.1223215991420282E3 * t10 * t25 * t7 +
        0.6116079957101409E2 * t84 * t9 * t35 -
        0.6116079957101409E2 * t84 * t9 * t59 +
        0.8119934108954854E1 * t92 * t90 * t6 * t5 +
        0.1952308078833175E2 * t97 * t35 - 0.3229720457288778E-2 * t101 * t59 +
        0.8119934108954854E1 * t92 * t105 * t5 -
        0.6459440914577556E-2 * t20 * t109 * t5 -
        0.1021327285067311E-1 * t20 * t113 * t5 -
        0.8039469615426877E3 * t97 * t7 -
        0.6459440914577556E-2 * t18 * t109 * t5 +
        0.8119934108954854E1 * t90 * t105 * t5 + 0.728715200486185E1 * t7 +
        0.6572479626389472E3 * t21 * t49 * t35 -
        0.1271151758223344E3 * t38 * t96 * t7;
      double t138 = t38 * t31 * t29;
      double t147 = t92 * t90 * t104;
      double t161 = t104 * t70;
      double t176 = t70 * t35;
      double t181 =
        -0.1952308078833175E2 * t97 * t59 +
        0.1271151758223344E3 * t138 * t41 * t59 -
        0.1271151758223344E3 * t138 * t7 -
        0.1271151758223344E3 * t138 * t41 * t35 +
        0.1099998793575276E3 * t147 * t71 - 0.45140105432805E1 * t38 * t61 * t5 +
        0.45140105432805E1 * t56 * t59 +
        0.2199997587150551E3 * t92 * t90 * t70 * t7 -
        0.4059967054477427E1 * t147 * t35 +
        0.2199997587150551E3 * t92 * t161 * t7 +
        0.4059967054477427E1 * t147 * t59 -
        0.1271151758223344E3 * t38 * t31 * t41 * t7 +
        0.2199997587150551E3 * t90 * t161 * t7 +
        0.3478501611135263E3 * t147 * t7 - 0.1099998793575276E3 * t147 * t176 +
        0.3229720457288778E-2 * t101 * t35;
      double t182 = t17 * t6;
      double t188 = t20 * t18;
      double t203 = t38 * t70;
      double t220 =
        -0.45140105432805E1 * t38 * t182 * t5 + 0.122935393065126E3 * t55 * t35 -
        0.1388758989328107E3 * t188 * t59 - 0.2777517978656213E3 * t113 * t5 -
        0.2777517978656213E3 * t20 * t6 * t5 + 0.1625923662282035E3 * t109 * t5 -
        0.8129618311410174E2 * t38 * t6 * t5 + 0.1388758989328107E3 * t188 * t35 -
        0.2570811037186818E2 * t203 * t59 + 0.2570811037186818E2 * t203 * t35 +
        0.122935393065126E3 * t182 * t5 + 0.122935393065126E3 * t61 * t5 -
        0.122935393065126E3 * t55 * t59 - 0.6572479626389472E3 * t21 * t49 * t59 -
        0.1932247451793815E1 * t176 - 0.1223215991420282E3 * t84 * t7;
      double t222 = t76 + t132 + t181 + t220;
      return t222;
    }

    // * f36 *********************************************************************

    double
    ortho2_f36 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t20 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t21 = -t3 - t8 - 0.5278640450004206E-1;
      double t23 = -t3 - t8 - 0.9472135954999579;
      double t27 = 0.1E1 * x;
      double t30 = (0.5 + t27 + t1) * (0.1154653670707977E1 + t27 + t1);
      double t31 = -0.1546536707079771 + t27 + t1;
      double t34 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t48 = 0.9472135954999579 + t27 + t1;
      double t50 = 0.5278640450004206E-1 + t27 + t1;
      double t58 = t20 * t6;
      double t62 = t50 * t48;
      double t78 = -t3 - t8 + 0.1546536707079771;
      double t79 = -t3 - t8 - 0.5;
      double t81 = -t3 - t8 - 0.1154653670707977E1;
      double MapleGenVar1 =
        0.3733028101803165E2 * (-t3 - t8 - 0.1265055323929465E1) * (-t3 - t8 -
                    0.7852315164806451)
        * (-t3 - t8 - 0.2147684835193549) * (-t3 - t8 + 0.2650553239294647) * t7 +
        0.567712971738462E2 * t23 * t21 * t20 * t7 -
        0.3379495617142792E3 * t34 * t31 * t30 * t7 +
        0.4896817197532036E3 * (-0.2650553239294647 + t27 +
              t1) * (0.2147684835193549 + t27 +
               t1) * (0.7852315164806451 + t27 +
                t1) * (0.1265055323929465E1 + t27 +
                 t1) * t7 +
        0.3823270589331786E3 * t50 * t48 * t6 * t5 +
        0.1562405607566647E3 * t23 * t21 * t6 * t5 +
        0.4532913346829967E2 * t34 * t58 * t5;
      double t91 =
        MapleGenVar1 + 0.8307493512880589E3 * t23 * t21 * t62 * t7 +
        0.2135883494598093E2 * t34 * t6 * t5 + 0.3199989401092891E2 * t58 * t5 +
        0.101679750743881E3 * t34 * t62 * t7 -
        0.202516380602846E3 * t31 * t30 * t7 -
        0.3356289097373485E2 * t81 * t79 * t78 * t7 -
        0.8926944484340175E2 * t81 * t79 * t78 * t20 * t7 +
        0.5255152376132534E2 * t7;
      return t91;
    }

    double
    ortho2_f36x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * x;
      double t7 = 0.1154653670707977E1 + t6 + t1;
      double t8 = 0.5 + t6 + t1;
      double t9 = t8 * t7;
      double t10 = -0.1546536707079771 + t6 + t1;
      double t11 = t10 * t9;
      double t15 = (-t3 - t1) * t2;
      double t16 = t4 * t15;
      double t20 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t21 = t20 * t10 * t8;
      double t28 = 0.1E1 * y;
      double t29 = -t3 - t28 + 0.2650553239294647;
      double t31 = -t3 - t28 - 0.2147684835193549;
      double t32 = -t3 - t28 - 0.7852315164806451;
      double t34 = -t3 - t28 - 0.1265055323929465E1;
      double t35 = t34 * t32 * t31;
      double t44 = t31 * t29;
      double t51 = -t3 - t28 - 0.5;
      double t53 = -t3 - t28 - 0.1154653670707977E1;
      double t57 = -t3 - t28 + 0.1546536707079771;
      double t58 = t57 * t4;
      double t63 = t53 * t51 * t57;
      double t77 = t7 * t4;
      double t84 =
        0.101258190301423E3 * t11 * t5 - 0.3379495617142792E3 * t21 * t16 -
        0.3379495617142792E3 * t20 * t10 * t7 * t16 +
        0.1866514050901582E2 * t35 * t29 * t15 -
        0.1866514050901582E2 * t35 * t16 -
        0.1866514050901582E2 * t34 * t32 * t29 * t16 -
        0.1866514050901582E2 * t34 * t44 * t16 -
        0.1866514050901582E2 * t32 * t44 * t16 +
        0.1678144548686743E2 * t53 * t51 * t4 * t15 +
        0.1678144548686743E2 * t53 * t58 * t15 + 0.1678144548686743E2 * t63 * t5 -
        0.1678144548686743E2 * t63 * t15 +
        0.1678144548686743E2 * t51 * t58 * t15 - 0.101258190301423E3 * t11 * t15 -
        0.202516380602846E3 * t10 * t8 * t4 * t15 -
        0.202516380602846E3 * t10 * t77 * t15 -
        0.202516380602846E3 * t8 * t77 * t15;
      double t85 = -t3 - t28 - 0.5278640450004206E-1;
      double t86 = -t3 - t28 - 0.9472135954999579;
      double t87 = t86 * t85;
      double t94 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t96 = t86 * t85 * t94;
      double t101 = t94 * t4;
      double t105 = t94 * t5;
      double t119 = t85 * t4;
      double t122 = t20 * t94;
      double t127 = 0.5278640450004206E-1 + t6 + t1;
      double t128 = t127 * t4;
      double t131 = 0.9472135954999579 + t6 + t1;
      double t132 = t131 * t4;
      double t135 = t127 * t131;
      double t138 =
        -0.7812028037833235E2 * t87 * t5 + 0.7812028037833235E2 * t87 * t15 -
        0.283856485869231E2 * t96 * t5 + 0.283856485869231E2 * t96 * t15 -
        0.283856485869231E2 * t86 * t101 * t15 - 0.1599994700546445E2 * t105 -
        0.283856485869231E2 * t85 * t101 * t15 - 0.6742126666006681E2 * t16 -
        0.7812028037833235E2 * t86 * t4 * t15 +
        0.7167165306079742E2 * t101 * t15 -
        0.1433433061215948E3 * t20 * t4 * t15 -
        0.7812028037833235E2 * t119 * t15 + 0.2266456673414983E2 * t122 * t15 -
        0.2266456673414983E2 * t122 * t5 + 0.3823270589331786E3 * t128 * t15 +
        0.3823270589331786E3 * t132 * t15 - 0.1911635294665893E3 * t135 * t5;
      double t143 = t94 * t15;
      double t162 = t57 * t94;
      double t171 = t86 * t85 * t127;
      double t189 =
        0.1689747808571396E3 * t21 * t7 * t5 + 0.1599994700546445E2 * t143 +
        0.1911635294665893E3 * t135 * t15 - 0.1067941747299047E2 * t20 * t5 +
        0.1067941747299047E2 * t20 * t15 - 0.1689747808571396E3 * t21 * t15 * t7 +
        0.4463472242170087E2 * t63 * t105 - 0.4463472242170087E2 * t63 * t143 +
        0.4463472242170087E2 * t53 * t51 * t94 * t16 +
        0.4463472242170087E2 * t53 * t162 * t16 +
        0.4463472242170087E2 * t51 * t162 * t16 -
        0.4153746756440294E3 * t171 * t131 * t5 +
        0.4153746756440294E3 * t171 * t131 * t15 +
        0.8307493512880589E3 * t171 * t16 +
        0.8307493512880589E3 * t86 * t85 * t131 * t16 -
        0.4153746756440294E3 * t86 * t135 * t16 -
        0.4153746756440294E3 * t85 * t135 * t16;
      double t190 = 0.1265055323929465E1 + t6 + t1;
      double t192 = 0.7852315164806451 + t6 + t1;
      double t193 = 0.2147684835193549 + t6 + t1;
      double t195 = -0.2650553239294647 + t6 + t1;
      double t196 = t195 * t193 * t192;
      double t210 = t192 * t190;
      double t228 = t20 * t135;
      double t244 =
        -0.2448408598766018E3 * t196 * t190 * t5 +
        0.2448408598766018E3 * t196 * t190 * t15 +
        0.4896817197532036E3 * t196 * t16 - 0.2627576188066267E2 * t5 +
        0.4896817197532036E3 * t195 * t193 * t190 * t16 +
        0.2627576188066267E2 * t15 + 0.4896817197532036E3 * t195 * t210 * t16 +
        0.4896817197532036E3 * t193 * t210 * t16 +
        0.160769802134432E3 * t127 * t132 * t15 -
        0.1795266047916341E3 * t86 * t119 * t15 - 0.53434517463638E3 * t11 * t16 -
        0.3379495617142792E3 * t20 * t9 * t16 - 0.508398753719405E2 * t228 * t5 +
        0.508398753719405E2 * t228 * t15 + 0.2822947711639227E3 * t63 * t16 +
        0.101679750743881E3 * t20 * t128 * t15 +
        0.101679750743881E3 * t20 * t132 * t15 -
        0.1866514050901582E2 * t35 * t29 * t5;
      double t246 = t84 + t138 + t189 + t244;
      return t246;
    }

    double
    ortho2_f36y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t10 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 = 0.1E1 * y;
      double t12 = -t3 - t11 + 0.1546536707079771;
      double t13 = t12 * t10;
      double t14 = -t3 - t11 - 0.1154653670707977E1;
      double t18 = t6 * t2;
      double t19 = 0.1E1 * x;
      double t20 = 0.9472135954999579 + t19 + t1;
      double t22 = 0.5278640450004206E-1 + t19 + t1;
      double t23 = -t3 - t11 - 0.5278640450004206E-1;
      double t25 = -t3 - t11 - 0.9472135954999579;
      double t26 = t25 * t23 * t22;
      double t31 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t34 = t22 * t20;
      double t38 = t10 * t18;
      double t39 = -t3 - t11 - 0.5;
      double t41 = t14 * t39 * t12;
      double t47 = 0.1265055323929465E1 + t19 + t1;
      double t48 = 0.2147684835193549 + t19 + t1;
      double t50 = -0.2650553239294647 + t19 + t1;
      double t57 = t25 * t23 * t10;
      double t60 = -t3 - t11 + 0.2650553239294647;
      double t62 = -t3 - t11 - 0.2147684835193549;
      double t63 = -t3 - t11 - 0.7852315164806451;
      double t65 = -t3 - t11 - 0.1265055323929465E1;
      double t66 = t65 * t63 * t62;
      double t69 = t20 * t6;
      double t73 = t6 * t4;
      double t74 = t34 * t31;
      double t79 = 0.1154653670707977E1 + t19 + t1;
      double t80 = 0.5 + t19 + t1;
      double t81 = t80 * t79;
      double t82 = -0.1546536707079771 + t19 + t1;
      double t83 = t82 * t81;
      double t89 = t31 * t82 * t80;
      double t92 =
        0.8926944484340175E2 * t14 * t13 * t7 -
        0.4153746756440294E3 * t26 * t20 * t18 -
        0.1067941747299047E2 * t31 * t18 - 0.8307493512880589E3 * t25 * t34 * t7 +
        0.4463472242170087E2 * t41 * t38 - 0.8307493512880589E3 * t23 * t34 * t7 +
        0.2448408598766018E3 * t50 * t48 * t47 * t7 +
        0.4153746756440294E3 * t26 * t7 - 0.283856485869231E2 * t57 * t18 -
        0.1866514050901582E2 * t66 * t60 * t18 +
        0.321539604268864E3 * t22 * t69 * t5 + 0.508398753719405E2 * t74 * t73 -
        0.1678144548686743E2 * t41 * t73 + 0.101258190301423E3 * t83 * t18 +
        0.2627576188066267E2 * t73 - 0.2627576188066267E2 * t18 -
        0.1689747808571396E3 * t89 * t7;
      double t93 = 0.7852315164806451 + t19 + t1;
      double t94 = t93 * t47;
      double t98 = t10 * t6;
      double t132 = t50 * t48 * t93;
      double t137 = t62 * t60;
      double t149 =
        0.2448408598766018E3 * t50 * t94 * t7 -
        0.567712971738462E2 * t25 * t98 * t5 +
        0.4153746756440294E3 * t26 * t20 * t73 -
        0.1689747808571396E3 * t31 * t82 * t79 * t7 +
        0.4153746756440294E3 * t25 * t23 * t20 * t7 + 0.16946291617641E2 * t7 -
        0.567712971738462E2 * t23 * t98 * t5 +
        0.1689747808571396E3 * t89 * t79 * t18 - 0.3733028101803165E2 * t66 * t7 -
        0.3733028101803165E2 * t65 * t63 * t60 * t7 +
        0.8926944484340175E2 * t14 * t39 * t10 * t7 -
        0.2448408598766018E3 * t132 * t47 * t18 +
        0.1067941747299047E2 * t31 * t73 -
        0.3733028101803165E2 * t65 * t137 * t7 +
        0.8926944484340175E2 * t39 * t13 * t7 -
        0.3733028101803165E2 * t63 * t137 * t7 - 0.106869034927276E4 * t83 * t7;
      double t161 = t23 * t6;
      double t165 = t10 * t73;
      double t178 = t12 * t6;
      double t185 = t22 * t6;
      double t197 =
        -0.1689747808571396E3 * t31 * t81 * t7 +
        0.3356289097373485E2 * t14 * t39 * t6 * t5 +
        0.2448408598766018E3 * t48 * t94 * t7 -
        0.8976330239581705E2 * t25 * t161 * t5 + 0.1599994700546445E2 * t165 -
        0.508398753719405E2 * t74 * t18 - 0.4463472242170087E2 * t41 * t165 +
        0.1411473855819613E3 * t41 * t7 + 0.1866514050901582E2 * t66 * t60 * t73 +
        0.1678144548686743E2 * t41 * t18 +
        0.3356289097373485E2 * t14 * t178 * t5 - 0.101258190301423E3 * t83 * t73 -
        0.1599994700546445E2 * t38 + 0.508398753719405E2 * t31 * t185 * t5 +
        0.283856485869231E2 * t57 * t73 + 0.3356289097373485E2 * t39 * t178 * t5 +
        0.2448408598766018E3 * t132 * t47 * t73;
      double t207 = t79 * t6;
      double t213 = t25 * t23;
      double t234 = t31 * t10;
      double t245 =
        0.2448408598766018E3 * t132 * t7 + 0.508398753719405E2 * t31 * t69 * t5 -
        0.101258190301423E3 * t82 * t80 * t6 * t5 -
        0.101258190301423E3 * t82 * t207 * t5 + 0.1911635294665893E3 * t34 * t73 -
        0.7812028037833235E2 * t213 * t18 -
        0.101258190301423E3 * t80 * t207 * t5 -
        0.1562405607566647E3 * t25 * t6 * t5 -
        0.1689747808571396E3 * t89 * t79 * t73 -
        0.1562405607566647E3 * t161 * t5 + 0.1433433061215948E3 * t98 * t5 -
        0.7167165306079742E2 * t31 * t6 * t5 + 0.7812028037833235E2 * t213 * t73 -
        0.2266456673414983E2 * t234 * t18 + 0.2266456673414983E2 * t234 * t73 -
        0.1911635294665893E3 * t34 * t18 + 0.1911635294665893E3 * t185 * t5 +
        0.1911635294665893E3 * t69 * t5;
      double t247 = t92 + t149 + t197 + t245;
      return t247;
    }

    // ORDER 8

    // Edge functions, order 8

    // number 37
    inline double
    ortho2_f37 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l2 * l3 * phi6 (l3 - l2);
    }

    inline double
    ortho2_f37x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l2x * l3 + l2 * l3x) * phi6 (l3 - l2) + l2 * l3 * phi6x (l3 -
                       l2) *
        (l3x - l2x);
    }

    inline double
    ortho2_f37y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l2y * l3 + l2 * l3y) * phi6 (l3 - l2) + l2 * l3 * phi6x (l3 -
                       l2) *
        (l3y - l2y);
    }

    // number 38
    inline double
    ortho2_f38 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l3 * l1 * phi6 (l1 - l3);
    }

    inline double
    ortho2_f38x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l3x * l1 + l3 * l1x) * phi6 (l1 - l3) + l3 * l1 * phi6x (l1 -
                       l3) *
        (l1x - l3x);
    }

    inline double
    ortho2_f38y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l3y * l1 + l3 * l1y) * phi6 (l1 - l3) + l3 * l1 * phi6x (l1 -
                       l3) *
        (l1y - l3y);
    }

    // number 39
    inline double
    ortho2_f39 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l1 * l2 * phi6 (l2 - l1);
    }

    inline double
    ortho2_f39x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l1x * l2 + l1 * l2x) * phi6 (l2 - l1) + l1 * l2 * phi6x (l2 -
                       l1) *
        (l2x - l1x);
    }

    inline double
    ortho2_f39y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l1y * l2 + l1 * l2y) * phi6 (l2 - l1) + l1 * l2 * phi6x (l2 -
                       l1) *
        (l2y - l1y);
    }

    // Bubble functions, order 8

    // * f40 *********************************************************************
    double
    ortho2_f40 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t20 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t21 = -t3 - t8 - 0.5278640450004206E-1;
      double t23 = -t3 - t8 - 0.9472135954999579;
      double t27 = 0.1E1 * x;
      double t30 = (0.5 + t27 + t1) * (0.1154653670707977E1 + t27 + t1);
      double t31 = -0.1546536707079771 + t27 + t1;
      double t34 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t48 = 0.9472135954999579 + t27 + t1;
      double t50 = 0.5278640450004206E-1 + t27 + t1;
      double t58 = t20 * t6;
      double t62 = t50 * t48;
      double t78 = -t3 - t8 + 0.1546536707079771;
      double t79 = -t3 - t8 - 0.5;
      double t81 = -t3 - t8 - 0.1154653670707977E1;
      double MapleGenVar1 =
        0.238446312838284E3 * (-t3 - t8 - 0.1265055323929465E1) * (-t3 - t8 -
                         0.7852315164806451)
        * (-t3 - t8 - 0.2147684835193549) * (-t3 - t8 + 0.2650553239294647) * t7 -
        0.2182504088190724E2 * t23 * t21 * t20 * t7 +
        0.2166031926508729E1 * t34 * t31 * t30 * t7 -
        0.103233637463148E1 * (-0.2650553239294647 + t27 +
             t1) * (0.2147684835193549 + t27 +
              t1) * (0.7852315164806451 + t27 +
               t1) * (0.1265055323929465E1 + t27 +
                t1) * t7 -
        0.1372528991489462E1 * t50 * t48 * t6 * t5 +
        0.1143842569064458E3 * t23 * t21 * t6 * t5 +
        0.1063094007467139E2 * t34 * t58 * t5 -
        0.8582957314521049E1 * t23 * t21 * t62 * t7;
      double t102 =
        MapleGenVar1 - 0.2029846685550216E2 * t34 * t6 * t5 -
        0.4054932721507037E1 * t58 * t5 - 0.166921527234769E2 * t34 * t62 * t7 +
        0.2605141925241498 * t31 * t30 * t7 +
        0.2645344493657566E3 * t81 * t79 * t78 * t7 -
        0.1248150209615594E3 * t81 * t79 * t78 * t20 * t7 +
        0.1205309328170986E2 * t7 + 0.5020599445851945E3 * (-t3 - t8 -
                  0.1330223896278567E1)
        * (-t3 - t8 - 0.9688487934707142) * t79 * (-t3 - t8 -
                     0.3115120652928579E-1) * (-t3 -
                       t8 +
                       0.3302238962785669)
        * t6 * t5;
      return t102;
    }

    double
    ortho2_f40x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * x;
      double t7 = 0.1265055323929465E1 + t6 + t1;
      double t9 = 0.7852315164806451 + t6 + t1;
      double t10 = 0.2147684835193549 + t6 + t1;
      double t12 = -0.2650553239294647 + t6 + t1;
      double t13 = t12 * t10 * t9;
      double t17 = (-t3 - t1) * t2;
      double t21 = t4 * t17;
      double t28 = t9 * t7;
      double t35 = 0.1E1 * y;
      double t36 = -t3 - t35 - 0.3115120652928579E-1;
      double t37 = -t3 - t35 - 0.5;
      double t39 = -t3 - t35 - 0.9688487934707142;
      double t40 = -t3 - t35 - 0.1330223896278567E1;
      double t41 = t40 * t39;
      double t42 = t41 * t37 * t36;
      double t45 = -t3 - t35 + 0.3302238962785669;
      double t56 = t36 * t45;
      double t68 = -t3 - t35 - 0.5278640450004206E-1;
      double t69 = -t3 - t35 - 0.9472135954999579;
      double t70 = t69 * t68;
      double t75 = t68 * t4;
      double t83 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t86 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t87 = t86 * t83;
      double t90 = t83 * t4;
      double t93 =
        0.5161681873157401 * t13 * t7 * t5 - 0.5161681873157401 * t13 * t7 * t17 -
        0.103233637463148E1 * t13 * t21 -
        0.103233637463148E1 * t12 * t10 * t7 * t21 -
        0.103233637463148E1 * t12 * t28 * t21 -
        0.103233637463148E1 * t10 * t28 * t21 - 0.2510299722925972E3 * t42 * t21 +
        0.2510299722925972E3 * t42 * t45 * t17 -
        0.2510299722925972E3 * t42 * t45 * t5 -
        0.2510299722925972E3 * t41 * t37 * t45 * t21 -
        0.2510299722925972E3 * t41 * t56 * t21 -
        0.2510299722925972E3 * t39 * t37 * t56 * t21 -
        0.2510299722925972E3 * t40 * t37 * t56 * t21 +
        0.5719212845322289E2 * t70 * t17 - 0.5719212845322289E2 * t70 * t5 -
        0.5719212845322289E2 * t75 * t17 - 0.5719212845322289E2 * t69 * t4 * t17 -
        0.5315470037335697E1 * t87 * t5 + 0.1680899215236106E2 * t90 * t17;
      double t99 = 0.9472135954999579 + t6 + t1;
      double t100 = 0.5278640450004206E-1 + t6 + t1;
      double t101 = t100 * t99;
      double t104 = t100 * t4;
      double t107 = t99 * t4;
      double t115 = t83 * t5;
      double t117 = t83 * t17;
      double t129 = 0.1154653670707977E1 + t6 + t1;
      double t130 = 0.5 + t6 + t1;
      double t131 = t130 * t129;
      double t132 = -0.1546536707079771 + t6 + t1;
      double t133 = t132 * t131;
      double t139 = t86 * t101;
      double t144 =
        -0.3361798430472211E2 * t86 * t4 * t17 +
        0.5315470037335697E1 * t87 * t17 - 0.6862644957447309 * t101 * t17 -
        0.1372528991489462E1 * t104 * t17 - 0.1372528991489462E1 * t107 * t17 +
        0.6862644957447309 * t101 * t5 - 0.1014923342775108E2 * t86 * t17 -
        0.1927187097770391E2 * t21 + 0.2027466360753519E1 * t115 -
        0.2027466360753519E1 * t117 + 0.1014923342775108E2 * t86 * t5 -
        0.2639261082878488E2 * t100 * t107 * t17 +
        0.6901683921311685E2 * t69 * t75 * t17 + 0.6026546640854932E1 * t17 -
        0.6026546640854932E1 * t5 + 0.3424797186205015E1 * t133 * t21 +
        0.2166031926508729E1 * t86 * t131 * t21 +
        0.8346076361738448E1 * t139 * t5 - 0.8346076361738448E1 * t139 * t17;
      double t146 = -t3 - t35 + 0.1546536707079771;
      double t148 = -t3 - t35 - 0.1154653670707977E1;
      double t149 = t148 * t37 * t146;
      double t158 = -t3 - t35 + 0.2650553239294647;
      double t160 = -t3 - t35 - 0.2147684835193549;
      double t161 = -t3 - t35 - 0.7852315164806451;
      double t163 = -t3 - t35 - 0.1265055323929465E1;
      double t164 = t163 * t161 * t160;
      double t170 = t86 * t132 * t130;
      double t186 = t160 * t158;
      double t197 = t146 * t4;
      double t214 =
        0.3946997524401873E3 * t149 * t21 -
        0.166921527234769E2 * t86 * t104 * t17 -
        0.166921527234769E2 * t86 * t107 * t17 -
        0.119223156419142E3 * t164 * t158 * t5 - 0.1302570962620749 * t133 * t5 +
        0.2166031926508729E1 * t170 * t21 +
        0.2166031926508729E1 * t86 * t132 * t129 * t21 +
        0.119223156419142E3 * t164 * t158 * t17 -
        0.119223156419142E3 * t164 * t21 -
        0.119223156419142E3 * t163 * t161 * t158 * t21 -
        0.119223156419142E3 * t163 * t186 * t21 -
        0.119223156419142E3 * t161 * t186 * t21 -
        0.1322672246828783E3 * t148 * t37 * t4 * t17 -
        0.1322672246828783E3 * t148 * t197 * t17 -
        0.1322672246828783E3 * t149 * t5 + 0.1322672246828783E3 * t149 * t17 -
        0.1322672246828783E3 * t37 * t197 * t17 +
        0.1302570962620749 * t133 * t17 +
        0.2605141925241498 * t132 * t130 * t4 * t17;
      double t215 = t129 * t4;
      double t223 = t69 * t68 * t83;
      double t248 = t146 * t83;
      double t257 = t69 * t68 * t100;
      double t275 =
        0.2605141925241498 * t132 * t215 * t17 +
        0.2605141925241498 * t130 * t215 * t17 +
        0.1091252044095362E2 * t223 * t5 - 0.1091252044095362E2 * t223 * t17 +
        0.1091252044095362E2 * t69 * t90 * t17 +
        0.1091252044095362E2 * t68 * t90 * t17 -
        0.1083015963254364E1 * t170 * t129 * t5 +
        0.1083015963254364E1 * t170 * t129 * t17 +
        0.624075104807797E2 * t149 * t115 - 0.624075104807797E2 * t149 * t117 +
        0.624075104807797E2 * t148 * t37 * t83 * t21 +
        0.624075104807797E2 * t148 * t248 * t21 +
        0.624075104807797E2 * t37 * t248 * t21 +
        0.4291478657260525E1 * t257 * t99 * t5 -
        0.4291478657260525E1 * t257 * t99 * t17 -
        0.8582957314521049E1 * t257 * t21 -
        0.8582957314521049E1 * t69 * t68 * t99 * t21 +
        0.4291478657260525E1 * t69 * t101 * t21 +
        0.4291478657260525E1 * t68 * t101 * t21;
      double t277 = t93 + t144 + t214 + t275;
      return t277;
    }

    double
    ortho2_f40y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t10 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 = 0.1E1 * y;
      double t12 = -t3 - t11 - 0.5;
      double t14 = -t3 - t11 - 0.1154653670707977E1;
      double t18 = 0.1E1 * x;
      double t19 = 0.7852315164806451 + t18 + t1;
      double t20 = 0.2147684835193549 + t18 + t1;
      double t22 = -0.2650553239294647 + t18 + t1;
      double t23 = t22 * t20 * t19;
      double t26 = t6 * t2;
      double t27 = 0.1154653670707977E1 + t18 + t1;
      double t29 = 0.5 + t18 + t1;
      double t30 = -0.1546536707079771 + t18 + t1;
      double t34 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t35 = t34 * t30 * t29;
      double t38 = t10 * t26;
      double t39 = -t3 - t11 + 0.1546536707079771;
      double t41 = t14 * t12 * t39;
      double t44 = t39 * t10;
      double t48 = t6 * t4;
      double t49 = 0.1265055323929465E1 + t18 + t1;
      double t53 = t19 * t49;
      double t57 = t10 * t48;
      double t60 = t10 * t6;
      double t61 = -t3 - t11 - 0.9472135954999579;
      double t65 = 0.9472135954999579 + t18 + t1;
      double t67 = 0.5278640450004206E-1 + t18 + t1;
      double t68 = -t3 - t11 - 0.5278640450004206E-1;
      double t70 = t61 * t68 * t67;
      double t73 = t67 * t65;
      double t82 = t29 * t27;
      double t83 = t30 * t82;
      double t86 = -t3 - t11 + 0.3302238962785669;
      double t87 = -t3 - t11 - 0.3115120652928579E-1;
      double t88 = t87 * t86;
      double t89 = -t3 - t11 - 0.9688487934707142;
      double t90 = -t3 - t11 - 0.1330223896278567E1;
      double t91 = t90 * t89;
      double t102 = t34 * t73;
      double t105 =
        0.1248150209615594E3 * t14 * t12 * t10 * t7 -
        0.5161681873157401 * t23 * t7 - 0.1083015963254364E1 * t35 * t27 * t26 +
        0.624075104807797E2 * t41 * t38 + 0.1248150209615594E3 * t12 * t44 * t7 -
        0.5161681873157401 * t23 * t49 * t48 -
        0.5161681873157401 * t22 * t53 * t7 - 0.624075104807797E2 * t41 * t57 +
        0.2182504088190724E2 * t61 * t60 * t5 +
        0.4291478657260525E1 * t70 * t65 * t26 +
        0.8582957314521049E1 * t61 * t73 * t7 +
        0.8582957314521049E1 * t68 * t73 * t7 - 0.4291478657260525E1 * t70 * t7 +
        0.684959437241003E1 * t83 * t7 - 0.5020599445851945E3 * t91 * t88 * t7 +
        0.1083015963254364E1 * t34 * t82 * t7 + 0.2027466360753519E1 * t38 +
        0.5161681873157401 * t23 * t49 * t26 + 0.8346076361738448E1 * t102 * t26;
      double t125 = t67 * t6;
      double t129 = t65 * t6;
      double t148 = -t3 - t11 + 0.2650553239294647;
      double t150 = -t3 - t11 - 0.2147684835193549;
      double t151 = -t3 - t11 - 0.7852315164806451;
      double t153 = -t3 - t11 - 0.1265055323929465E1;
      double t154 = t153 * t151 * t150;
      double t160 = t91 * t12 * t87;
      double t171 =
        0.1973498762200937E3 * t41 * t7 + 0.2182504088190724E2 * t68 * t60 * t5 +
        0.1248150209615594E3 * t14 * t44 * t7 -
        0.5020599445851945E3 * t90 * t12 * t88 * t7 -
        0.5161681873157401 * t20 * t53 * t7 -
        0.5161681873157401 * t22 * t20 * t49 * t7 -
        0.8346076361738448E1 * t34 * t125 * t5 -
        0.8346076361738448E1 * t34 * t129 * t5 -
        0.4291478657260525E1 * t61 * t68 * t65 * t7 +
        0.1083015963254364E1 * t35 * t27 * t48 -
        0.5020599445851945E3 * t89 * t12 * t88 * t7 +
        0.1014923342775108E2 * t34 * t26 + 0.1322672246828783E3 * t41 * t48 -
        0.119223156419142E3 * t154 * t148 * t26 - 0.2027466360753519E1 * t57 +
        0.2510299722925972E3 * t160 * t86 * t48 - 0.1302570962620749 * t83 * t26 +
        0.1083015963254364E1 * t35 * t7 +
        0.1083015963254364E1 * t34 * t30 * t27 * t7;
      double t190 = t150 * t148;
      double t200 = t68 * t6;
      double t204 = t61 * t68;
      double t216 = t39 * t6;
      double t222 =
        -0.238446312838284E3 * t154 * t7 - 0.6862644957447309 * t73 * t48 -
        0.6026546640854932E1 * t26 - 0.5777797669346903E2 * t7 -
        0.6862644957447309 * t125 * t5 + 0.6026546640854932E1 * t48 -
        0.238446312838284E3 * t153 * t151 * t148 * t7 -
        0.1014923342775108E2 * t34 * t48 + 0.6862644957447309 * t73 * t26 -
        0.238446312838284E3 * t153 * t190 * t7 -
        0.5278522165756976E2 * t67 * t129 * t5 -
        0.238446312838284E3 * t151 * t190 * t7 +
        0.3450841960655843E2 * t61 * t200 * t5 -
        0.5719212845322289E2 * t204 * t26 -
        0.4291478657260525E1 * t70 * t65 * t48 -
        0.2645344493657566E3 * t14 * t12 * t6 * t5 -
        0.6862644957447309 * t129 * t5 - 0.2645344493657566E3 * t14 * t216 * t5 -
        0.1322672246828783E3 * t41 * t26;
      double t227 = t61 * t68 * t10;
      double t246 = t27 * t6;
      double t250 = t34 * t10;
      double t276 =
        -0.1143842569064458E3 * t61 * t6 * t5 -
        0.1091252044095362E2 * t227 * t48 -
        0.2645344493657566E3 * t12 * t216 * t5 -
        0.1143842569064458E3 * t200 * t5 -
        0.2510299722925972E3 * t160 * t86 * t26 +
        0.1302570962620749 * t30 * t29 * t6 * t5 +
        0.5719212845322289E2 * t204 * t48 + 0.3361798430472211E2 * t60 * t5 +
        0.1302570962620749 * t30 * t246 * t5 - 0.5315470037335697E1 * t250 * t26 -
        0.5020599445851945E3 * t160 * t7 + 0.5315470037335697E1 * t250 * t48 +
        0.1302570962620749 * t29 * t246 * t5 - 0.8346076361738448E1 * t102 * t48 -
        0.5020599445851945E3 * t91 * t12 * t86 * t7 -
        0.1680899215236106E2 * t34 * t6 * t5 +
        0.119223156419142E3 * t154 * t148 * t48 + 0.1302570962620749 * t83 * t48 +
        0.1091252044095362E2 * t227 * t26;
      double t278 = t105 + t171 + t222 + t276;
      return t278;
    }

    // * f41 *********************************************************************

    double
    ortho2_f41 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = t10 * t5;
      double t12 = 0.1E1 * y;
      double t19 =
        (-t3 - t12 - 0.1265055323929465E1) * (-t3 - t12 -
                0.7852315164806451) * (-t3 - t12 -
                     0.2147684835193549)
        * (-t3 - t12 + 0.2650553239294647);
      double t22 = t6 * t5;
      double t25 = -t3 - t12 - 0.5278640450004206E-1;
      double t27 = -t3 - t12 - 0.9472135954999579;
      double t31 = 0.1E1 * x;
      double t34 = (0.5 + t31 + t1) * (0.1154653670707977E1 + t31 + t1);
      double t35 = -0.1546536707079771 + t31 + t1;
      double t38 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t52 = 0.9472135954999579 + t31 + t1;
      double t54 = 0.5278640450004206E-1 + t31 + t1;
      double t65 = t54 * t52;
      double t80 = -t3 - t12 + 0.1546536707079771;
      double t81 = -t3 - t12 - 0.5;
      double t83 = -t3 - t12 - 0.1154653670707977E1;
      double MapleGenVar1 =
        -0.4640828799519412E3 * t19 * t11 + 0.5661821975206096E3 * t19 * t22 -
        0.2317957249892383E3 * t27 * t25 * t9 * t22 +
        0.3033443353082674E1 * t38 * t35 * t34 * t22 -
        0.3049013103392047E1 * (-0.2650553239294647 + t31 +
              t1) * (0.2147684835193549 + t31 +
               t1) * (0.7852315164806451 + t31 +
                t1) * (0.1265055323929465E1 + t31 +
                 t1) * t22 +
        0.1000636590337293E3 * t54 * t52 * t6 * t5 +
        0.4020717491243223E3 * t27 * t25 * t6 * t5 +
        0.2958674505593556E2 * t38 * t10 * t5;
      double t104 =
        MapleGenVar1 + 0.6888869135866671E3 * t27 * t25 * t65 * t22 -
        0.4451610833825012E2 * t38 * t6 * t5 - 0.3027495223509753E2 * t11 -
        0.2268937223381916E2 * t38 * t65 * t22 +
        0.3145273981003103E2 * t35 * t34 * t22 +
        0.5644101376801965E3 * t83 * t81 * t80 * t22 -
        0.2625126183009276E3 * t83 * t81 * t80 * t9 * t22 +
        0.4889912996252441E2 * t22 + 0.7824599923319738E3 * (-t3 - t12 -
                   0.1330223896278567E1)
        * (-t3 - t12 - 0.9688487934707142) * t81 * (-t3 - t12 -
                0.3115120652928579E-1) *
        (-t3 - t12 + 0.3302238962785669) * t6 * t5;
      return t104;
    }

    double
    ortho2_f41x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t9 = -t3 - t8 + 0.2650553239294647;
      double t10 = -t3 - t8 - 0.2147684835193549;
      double t11 = t10 * t9;
      double t12 = -t3 - t8 - 0.7852315164806451;
      double t13 = -t3 - t8 - 0.1265055323929465E1;
      double t14 = t13 * t12;
      double t15 = t14 * t11;
      double t18 = 0.1E1 * x;
      double t19 = 0.9472135954999579 + t18 + t1;
      double t20 = t19 * t6;
      double t21 = 0.5278640450004206E-1 + t18 + t1;
      double t25 = -t3 - t8 - 0.5278640450004206E-1;
      double t26 = t25 * t6;
      double t27 = -t3 - t8 - 0.9472135954999579;
      double t31 = t6 * t2;
      double t32 = -t3 - t8 + 0.3302238962785669;
      double t34 = -t3 - t8 - 0.3115120652928579E-1;
      double t35 = -t3 - t8 - 0.5;
      double t37 = -t3 - t8 - 0.9688487934707142;
      double t38 = -t3 - t8 - 0.1330223896278567E1;
      double t39 = t38 * t37;
      double t40 = t39 * t35 * t34;
      double t52 = t34 * t32;
      double t56 = 0.1154653670707977E1 + t18 + t1;
      double t57 = 0.5 + t18 + t1;
      double t58 = t57 * t56;
      double t59 = -0.1546536707079771 + t18 + t1;
      double t60 = t59 * t58;
      double t65 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t69 = t21 * t19;
      double t70 = t65 * t69;
      double t75 = -t3 - t8 + 0.1546536707079771;
      double t77 = -t3 - t8 - 0.1154653670707977E1;
      double t78 = t77 * t35 * t75;
      double t88 = t21 * t6;
      double t94 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t95 = t94 * t31;
      double t102 = t12 * t10;
      double t103 = t13 * t102;
      double t108 =
        0.1467558923738627E4 * t15 * t7 - 0.3587504746912553E2 * t21 * t20 * t5 +
        0.7330024428560015E3 * t27 * t26 * t5 -
        0.3912299961659869E3 * t40 * t32 * t31 +
        0.3912299961659869E3 * t40 * t32 * t5 - 0.3912299961659869E3 * t40 * t7 -
        0.3912299961659869E3 * t39 * t35 * t32 * t7 -
        0.3912299961659869E3 * t39 * t52 * t7 + 0.4796295074419801E1 * t60 * t7 +
        0.3033443353082674E1 * t65 * t58 * t7 - 0.1134468611690958E2 * t70 * t5 +
        0.1134468611690958E2 * t70 * t31 + 0.8301377883653322E3 * t78 * t7 -
        0.3912299961659869E3 * t38 * t35 * t52 * t7 -
        0.2268937223381916E2 * t65 * t20 * t5 -
        0.2268937223381916E2 * t65 * t88 * t5 + 0.1513747611754877E2 * t95 -
        0.3912299961659869E3 * t37 * t35 * t52 * t7 -
        0.2830910987603048E3 * t103 * t9 * t31 - 0.1572636990501552E2 * t60 * t31;
      double t110 = t65 * t59 * t57;
      double t136 = t75 * t6;
      double t149 = t56 * t6;
      double t168 = t94 * t5;
      double t171 = t9 * t94;
      double t175 =
        -0.2822050688400982E3 * t78 * t31 + 0.2822050688400982E3 * t78 * t5 -
        0.2822050688400982E3 * t35 * t136 * t5 +
        0.3145273981003103E2 * t59 * t149 * t5 +
        0.3145273981003103E2 * t59 * t57 * t6 * t5 +
        0.1572636990501552E2 * t60 * t5 + 0.3145273981003103E2 * t57 * t149 * t5 +
        0.2320414399759706E3 * t15 * t95 +
        0.2320414399759706E3 * t14 * t10 * t94 * t7 -
        0.2320414399759706E3 * t15 * t168 +
        0.2320414399759706E3 * t14 * t171 * t7;
      double t186 = t27 * t25 * t94;
      double t191 = t94 * t6;
      double t213 = t75 * t94;
      double t222 = t27 * t25 * t21;
      double t240 = 0.1265055323929465E1 + t18 + t1;
      double t242 = 0.7852315164806451 + t18 + t1;
      double t243 = 0.2147684835193549 + t18 + t1;
      double t245 = -0.2650553239294647 + t18 + t1;
      double t246 = t245 * t243 * t242;
      double t251 =
        0.1312563091504638E3 * t77 * t35 * t94 * t7 +
        0.1312563091504638E3 * t35 * t213 * t7 +
        0.1312563091504638E3 * t77 * t213 * t7 -
        0.3444434567933336E3 * t222 * t19 * t31 +
        0.3444434567933336E3 * t222 * t19 * t5 +
        0.6888869135866671E3 * t222 * t7 - 0.3444434567933336E3 * t27 * t69 * t7 +
        0.6888869135866671E3 * t27 * t25 * t19 * t7 -
        0.3444434567933336E3 * t25 * t69 * t7 +
        0.1524506551696024E1 * t246 * t240 * t31 -
        0.3049013103392047E1 * t246 * t7;
      double t256 = t242 * t240;
      double t267 = t27 * t25;
      double t283 = t65 * t94;
      double t303 =
        -0.2225805416912506E2 * t65 * t5 - 0.1479337252796778E2 * t283 * t31 +
        0.4678075146374114E2 * t191 * t5 - 0.9356150292748228E2 * t65 * t6 * t5 +
        0.1479337252796778E2 * t283 * t5 - 0.5003182951686463E2 * t69 * t31 +
        0.1000636590337293E3 * t88 * t5 + 0.1000636590337293E3 * t20 * t5 +
        0.5003182951686463E2 * t69 * t5 - 0.2444956498126221E2 * t31 +
        0.2444956498126221E2 * t5;
      double MapleGenVar1 =
        t175 - 0.1513747611754877E2 * t168 - 0.2830910987603048E3 * t103 * t7 +
        0.1158978624946191E3 * t186 * t31 - 0.1158978624946191E3 * t186 * t5 -
        0.1312563091504638E3 * t78 * t168 + 0.1312563091504638E3 * t78 * t95 +
        0.3033443353082674E1 * t110 * t7 + 0.2010358745621611E3 * t267 * t5 +
        0.2225805416912506E2 * t65 * t31 - 0.2010358745621611E3 * t26 * t5 +
        t108 - 0.3049013103392047E1 * t243 * t256 * t7 +
        0.2320414399759706E3 * t13 * t10 * t171 * t7 + 0.2535165765787183E2 * t7 +
        0.2320414399759706E3 * t102 * t171 * t7 +
        0.1158978624946191E3 * t25 * t191 * t5;
      double t306 =
        MapleGenVar1 + 0.1158978624946191E3 * t27 * t191 * t5 -
        0.1516721676541337E1 * t110 * t56 * t31 +
        0.1516721676541337E1 * t110 * t56 * t5 + t251 +
        0.3033443353082674E1 * t65 * t59 * t56 * t7 + t303 -
        0.2830910987603048E3 * t13 * t12 * t9 * t7 -
        0.2822050688400982E3 * t77 * t35 * t6 * t5 -
        0.2010358745621611E3 * t27 * t6 * t5 - 0.2010358745621611E3 * t267 * t31 +
        0.2830910987603048E3 * t103 * t9 * t5 -
        0.2830910987603048E3 * t13 * t11 * t7 -
        0.1524506551696024E1 * t246 * t240 * t5 -
        0.2830910987603048E3 * t12 * t11 * t7 -
        0.3049013103392047E1 * t245 * t243 * t240 * t7 -
        0.2822050688400982E3 * t77 * t136 * t5 -
        0.3049013103392047E1 * t245 * t256 * t7;
      return t306;
    }

    double
    ortho2_f41y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.9472135954999579 + t8 + t1;
      double t10 = 0.1E1 * y;
      double t11 = -t3 - t10 - 0.5278640450004206E-1;
      double t13 = -t3 - t10 - 0.9472135954999579;
      double t17 = 0.5278640450004206E-1 + t8 + t1;
      double t18 = t17 * t9;
      double t22 = t11 * t6;
      double t25 = -t3 - t10 - 0.5;
      double t27 = -t3 - t10 - 0.1154653670707977E1;
      double t33 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t34 = t33 * t6;
      double t39 = t13 * t11 * t17;
      double t42 = -t3 - t10 + 0.3302238962785669;
      double t44 = -t3 - t10 - 0.9688487934707142;
      double t45 = -t3 - t10 - 0.1330223896278567E1;
      double t46 = t45 * t44;
      double t50 = -t3 - t10 + 0.1546536707079771;
      double t51 = t50 * t33;
      double t55 = t6 * t4;
      double t59 = 0.1154653670707977E1 + t8 + t1;
      double t60 = t59 * t6;
      double t61 = 0.5 + t8 + t1;
      double t71 = t13 * t11 * t33;
      double t75 = t27 * t25 * t50;
      double t78 = t9 * t6;
      double t83 = t6 * t2;
      double t85 = -t3 - t10 - 0.3115120652928579E-1;
      double t87 = t46 * t25 * t85;
      double t90 = -t3 - t10 + 0.2650553239294647;
      double t91 = -t3 - t10 - 0.2147684835193549;
      double t92 = t91 * t90;
      double t93 = -t3 - t10 - 0.7852315164806451;
      double t94 = -t3 - t10 - 0.1265055323929465E1;
      double t95 = t94 * t93;
      double t96 = t95 * t92;
      double t100 = t13 * t11;
      double t103 =
        0.3444434567933336E3 * t13 * t11 * t9 * t7 -
        0.6888869135866671E3 * t13 * t18 * t7 - 0.4020717491243223E3 * t22 * t5 -
        0.5644101376801965E3 * t27 * t25 * t6 * t5 +
        0.2317957249892383E3 * t11 * t34 * t5 + 0.3444434567933336E3 * t39 * t7 -
        0.7824599923319738E3 * t46 * t25 * t42 * t7 +
        0.2625126183009276E3 * t27 * t51 * t7 +
        0.3444434567933336E3 * t39 * t9 * t55 +
        0.1572636990501552E2 * t61 * t60 * t5 + 0.5003182951686463E2 * t18 * t55 +
        0.3665012214280007E3 * t13 * t22 * t5 - 0.1158978624946191E3 * t71 * t55 +
        0.2822050688400982E3 * t75 * t55 - 0.7175009493825106E2 * t17 * t78 * t5 +
        0.2444956498126221E2 * t55 - 0.3912299961659869E3 * t87 * t42 * t83 +
        0.7337794618693137E3 * t96 * t7 - 0.2444956498126221E2 * t83 -
        0.2010358745621611E3 * t100 * t83;
      double t104 = t90 * t33;
      double t108 = t61 * t59;
      double t109 = -0.1546536707079771 + t8 + t1;
      double t110 = t109 * t108;
      double t117 = t50 * t6;
      double t121 = t85 * t42;
      double t127 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t131 = t127 * t18;
      double t134 = t93 * t91;
      double t135 = t94 * t134;
      double t138 = t33 * t55;
      double t152 = t33 * t83;
      double t177 =
        -0.7824599923319738E3 * t45 * t25 * t121 * t7 +
        0.4150688941826661E3 * t75 * t7 - 0.1134468611690958E2 * t131 * t55 +
        0.2320414399759706E3 * t96 * t152 -
        0.5661821975206096E3 * t94 * t93 * t90 * t7 -
        0.5661821975206096E3 * t94 * t92 * t7 - 0.2822050688400982E3 * t75 * t83 -
        0.5661821975206096E3 * t93 * t92 * t7 +
        0.4640828799519412E3 * t134 * t104 * t7 +
        0.3912299961659869E3 * t87 * t42 * t55 +
        0.1516721676541337E1 * t127 * t109 * t59 * t7;
      double t186 = t17 * t6;
      double t200 = t127 * t109 * t61;
      double t203 = t127 * t33;
      double t237 =
        -0.5003182951686463E2 * t18 * t83 + 0.5003182951686463E2 * t78 * t5 +
        0.5003182951686463E2 * t186 * t5 - 0.1572636990501552E2 * t110 * t83 -
        0.2830910987603048E3 * t135 * t90 * t83 -
        0.4020717491243223E3 * t13 * t6 * t5 +
        0.4640828799519412E3 * t95 * t91 * t33 * t7 +
        0.1572636990501552E2 * t109 * t61 * t6 * t5 +
        0.2830910987603048E3 * t135 * t90 * t55 +
        0.1572636990501552E2 * t109 * t60 * t5 +
        0.1572636990501552E2 * t110 * t55;
      double t247 = 0.1265055323929465E1 + t8 + t1;
      double t249 = 0.7852315164806451 + t8 + t1;
      double t250 = 0.2147684835193549 + t8 + t1;
      double t252 = -0.2650553239294647 + t8 + t1;
      double t253 = t252 * t250 * t249;
      double t283 = t249 * t247;
      double t304 =
        0.1524506551696024E1 * t253 * t247 * t83 -
        0.1524506551696024E1 * t253 * t7 -
        0.1524506551696024E1 * t252 * t250 * t247 * t7 -
        0.1524506551696024E1 * t252 * t283 * t7 -
        0.1524506551696024E1 * t250 * t283 * t7 -
        0.6888869135866671E3 * t11 * t18 * t7 -
        0.2320414399759706E3 * t96 * t138 + 0.2010358745621611E3 * t100 * t55 -
        0.4678075146374114E2 * t127 * t6 * t5 + 0.9356150292748228E2 * t34 * t5 -
        0.1479337252796778E2 * t203 * t83;
      double MapleGenVar1 =
        t177 - 0.7824599923319738E3 * t46 * t121 * t7 + t103 +
        0.1513747611754877E2 * t152 + 0.9592590148839602E1 * t110 * t7 + t237 +
        t304 - 0.1513747611754877E2 * t138 - 0.9290339235782683E2 * t7 -
        0.3444434567933336E3 * t39 * t9 * t83 -
        0.1516721676541337E1 * t200 * t59 * t83 +
        0.4640828799519412E3 * t94 * t91 * t104 * t7 +
        0.1516721676541337E1 * t127 * t108 * t7 -
        0.5644101376801965E3 * t27 * t117 * t5 +
        0.4640828799519412E3 * t95 * t104 * t7 -
        0.5644101376801965E3 * t25 * t117 * t5 +
        0.1134468611690958E2 * t131 * t83;
      double t307 =
        MapleGenVar1 - 0.5661821975206096E3 * t135 * t7 +
        0.2225805416912506E2 * t127 * t83 - 0.2225805416912506E2 * t127 * t55 -
        0.7824599923319738E3 * t87 * t7 + 0.1516721676541337E1 * t200 * t7 +
        0.1479337252796778E2 * t203 * t55 + 0.1158978624946191E3 * t71 * t83 -
        0.1312563091504638E3 * t75 * t138 + 0.1312563091504638E3 * t75 * t152 -
        0.1134468611690958E2 * t127 * t186 * t5 -
        0.1134468611690958E2 * t127 * t78 * t5 +
        0.2317957249892383E3 * t13 * t34 * t5 -
        0.1524506551696024E1 * t253 * t247 * t55 +
        0.1516721676541337E1 * t200 * t59 * t55 +
        0.2625126183009276E3 * t25 * t51 * t7 +
        0.2625126183009276E3 * t27 * t25 * t33 * t7 -
        0.7824599923319738E3 * t44 * t25 * t121 * t7;
      return t307;
    }

    // * f42 *********************************************************************

    double
    ortho2_f42 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = t10 * t5;
      double t12 = 0.1E1 * y;
      double t19 =
        (-t3 - t12 - 0.1265055323929465E1) * (-t3 - t12 -
                0.7852315164806451) * (-t3 - t12 -
                     0.2147684835193549)
        * (-t3 - t12 + 0.2650553239294647);
      double t22 = 0.1E1 * x;
      double t23 = 0.9472135954999579 + t22 + t1;
      double t24 = t23 * t6;
      double t26 = 0.5278640450004206E-1 + t22 + t1;
      double t27 = -t3 - t12 + 0.1546536707079771;
      double t29 = -t3 - t12 - 0.5;
      double t30 = -t3 - t12 - 0.1154653670707977E1;
      double t31 = t30 * t29;
      double t35 = t6 * t5;
      double t38 = -t3 - t12 - 0.5278640450004206E-1;
      double t40 = -t3 - t12 - 0.9472135954999579;
      double t46 = (0.5 + t22 + t1) * (0.1154653670707977E1 + t22 + t1);
      double t47 = -0.1546536707079771 + t22 + t1;
      double t50 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t74 = t26 * t23;
      double MapleGenVar1 =
        -0.8297514121605843E3 * t19 * t11 +
        0.2447460819436755E4 * t31 * t27 * t26 * t24 * t5 +
        0.4784957998111317E3 * t19 * t35 -
        0.4530409344362302E3 * t40 * t38 * t9 * t35 -
        0.1859379908129035E3 * t50 * t47 * t46 * t35 -
        0.1416473927093414E1 * (-0.2650553239294647 + t22 +
              t1) * (0.2147684835193549 + t22 +
               t1) * (0.7852315164806451 + t22 +
                t1) * (0.1265055323929465E1 + t22 +
                 t1) * t35 +
        0.2065065364391397E3 * t26 * t24 * t5 +
        0.4660632883609144E3 * t40 * t38 * t6 * t5 +
        0.6852464824289748E2 * t50 * t10 * t5;
      double t109 =
        MapleGenVar1 + 0.1099160133578124E4 * t40 * t38 * t74 * t35 -
        0.8785967965229223E2 * t50 * t6 * t5 - 0.5308071599948685E2 * t11 -
        0.2454458174580278E3 * t50 * t74 * t35 +
        0.1276659334611126E2 * t47 * t46 * t35 +
        0.1002384876260218E4 * t30 * t29 * t27 * t35 -
        0.4141189246796004E3 * t31 * t27 * t9 * t35 + 0.6834852260962291E2 * t35 +
        0.700454232477025E3 * (-t3 - t12 - 0.1330223896278567E1) * (-t3 - t12 -
                    0.9688487934707142)
        * t29 * (-t3 - t12 - 0.3115120652928579E-1) * (-t3 - t12 +
                   0.3302238962785669) * t6 *
        t5;
      return t109;
    }

    double
    ortho2_f42x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.1E1 * y;
      double t7 = -t3 - t6 + 0.3302238962785669;
      double t9 = -t3 - t6 - 0.3115120652928579E-1;
      double t10 = -t3 - t6 - 0.5;
      double t12 = -t3 - t6 - 0.9688487934707142;
      double t13 = -t3 - t6 - 0.1330223896278567E1;
      double t14 = t13 * t12;
      double t15 = t14 * t10 * t9;
      double t18 = 0.5 + t3;
      double t19 = t18 * t2;
      double t23 = t18 * t5;
      double t24 = -t3 - t6 + 0.2650553239294647;
      double t25 = -t3 - t6 - 0.2147684835193549;
      double t26 = t25 * t24;
      double t27 = -t3 - t6 - 0.7852315164806451;
      double t28 = -t3 - t6 - 0.1265055323929465E1;
      double t29 = t28 * t27;
      double t30 = t29 * t26;
      double t39 = 0.1E1 * x;
      double t40 = 0.1154653670707977E1 + t39 + t1;
      double t41 = 0.5 + t39 + t1;
      double t42 = t41 * t40;
      double t43 = -0.1546536707079771 + t39 + t1;
      double t44 = t43 * t42;
      double t47 = t9 * t7;
      double t53 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t57 = 0.9472135954999579 + t39 + t1;
      double t58 = 0.5278640450004206E-1 + t39 + t1;
      double t59 = t58 * t57;
      double t60 = t53 * t59;
      double t65 = -t3 - t6 + 0.1546536707079771;
      double t66 = t10 * t65;
      double t67 = -t3 - t6 - 0.1154653670707977E1;
      double t68 = t67 * t66;
      double t71 =
        0.3502271162385125E3 * t15 * t7 * t5 -
        0.3502271162385125E3 * t15 * t7 * t19 + 0.2623904354168581E4 * t30 * t23 -
        0.3502271162385125E3 * t14 * t10 * t7 * t23 -
        0.3502271162385125E3 * t15 * t23 - 0.293993777262119E3 * t44 * t23 -
        0.3502271162385125E3 * t14 * t47 * t23 -
        0.1859379908129035E3 * t53 * t42 * t23 +
        0.1227229087290139E3 * t60 * t19 - 0.1227229087290139E3 * t60 * t5 +
        0.1309559024167252E4 * t68 * t23;
      double t76 = t58 * t18;
      double t80 = t57 * t18;
      double t89 = t27 * t25;
      double t90 = t28 * t89;
      double t96 = t53 * t43 * t41;
      double t109 =
        -0.3502271162385125E3 * t13 * t10 * t47 * t23 -
        0.2454458174580278E3 * t53 * t76 * t5 -
        0.2454458174580278E3 * t53 * t80 * t5 -
        0.3502271162385125E3 * t12 * t10 * t47 * t23 -
        0.2392478999055658E3 * t90 * t24 * t19 -
        0.6383296673055628E1 * t44 * t19 - 0.1859379908129035E3 * t96 * t23 -
        0.1859379908129035E3 * t53 * t43 * t40 * t23 +
        0.2392478999055658E3 * t90 * t24 * t5 + 0.289376112939225E2 * t23 -
        0.2392478999055658E3 * t90 * t23;
      double t120 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t121 = t120 * t19;
      double t130 = t120 * t5;
      double t136 = t65 * t18;
      double t144 =
        -0.2392478999055658E3 * t28 * t27 * t24 * t23 -
        0.2392478999055658E3 * t28 * t26 * t23 + 0.2654035799974342E2 * t121 -
        0.4392983982614612E2 * t53 * t5 - 0.2392478999055658E3 * t27 * t26 * t23 +
        0.4392983982614612E2 * t53 * t19 - 0.2654035799974342E2 * t130 -
        0.5011924381301092E3 * t67 * t10 * t18 * t5 -
        0.5011924381301092E3 * t67 * t136 * t5 -
        0.5011924381301092E3 * t68 * t19 + 0.5011924381301092E3 * t68 * t5;
      double t154 = t40 * t18;
      double t169 = t24 * t120;
      double t180 = -t3 - t6 - 0.5278640450004206E-1;
      double t182 = -t3 - t6 - 0.9472135954999579;
      double t183 = t182 * t180 * t120;
      double t186 =
        -0.5011924381301092E3 * t10 * t136 * t5 +
        0.6383296673055628E1 * t44 * t5 +
        0.1276659334611126E2 * t43 * t41 * t18 * t5 +
        0.1276659334611126E2 * t43 * t154 * t5 +
        0.1276659334611126E2 * t41 * t154 * t5 +
        0.4148757060802921E3 * t30 * t121 - 0.4148757060802921E3 * t30 * t130 +
        0.4148757060802921E3 * t29 * t25 * t120 * t23 +
        0.4148757060802921E3 * t29 * t169 * t23 +
        0.4148757060802921E3 * t28 * t25 * t169 * t23 +
        0.4148757060802921E3 * t89 * t169 * t23 +
        0.2265204672181151E3 * t183 * t19;
      double t191 = t120 * t18;
      double t204 = t182 * t180;
      double t217 = t65 * t120;
      double t221 =
        -0.2265204672181151E3 * t183 * t5 +
        0.2265204672181151E3 * t182 * t191 * t5 +
        0.2265204672181151E3 * t180 * t191 * t5 +
        0.9296899540645173E2 * t96 * t40 * t19 -
        0.9296899540645173E2 * t96 * t40 * t5 -
        0.2330316441804572E3 * t204 * t19 + 0.2330316441804572E3 * t204 * t5 +
        0.2070594623398002E3 * t68 * t121 - 0.2070594623398002E3 * t68 * t130 +
        0.2070594623398002E3 * t67 * t10 * t120 * t23 +
        0.2070594623398002E3 * t67 * t217 * t23;
      double t228 = t57 * t19;
      double t230 = t67 * t10;
      double t231 = t230 * t65 * t58;
      double t234 = t57 * t5;
      double t239 = t180 * t18;
      double t247 = t53 * t120;
      double t253 = t182 * t180 * t58;
      double t256 =
        0.2070594623398002E3 * t10 * t217 * t23 -
        0.2330316441804572E3 * t182 * t18 * t5 -
        0.1223730409718377E4 * t231 * t228 + 0.1223730409718377E4 * t231 * t234 +
        0.2447460819436755E4 * t231 * t23 - 0.2330316441804572E3 * t239 * t5 +
        0.1083469821547055E3 * t191 * t5 - 0.2166939643094111E3 * t53 * t18 * t5 +
        0.3426232412144874E2 * t247 * t5 - 0.3426232412144874E2 * t247 * t19 -
        0.5495800667890619E3 * t253 * t228;
      double t288 =
        0.2065065364391397E3 * t80 * t5 - 0.1032532682195698E3 * t59 * t19 +
        0.1032532682195698E3 * t59 * t5 + 0.2065065364391397E3 * t76 * t5 +
        0.5495800667890619E3 * t253 * t234 +
        0.1099160133578124E4 * t182 * t180 * t57 * t23 +
        0.1099160133578124E4 * t253 * t23 -
        0.5495800667890619E3 * t182 * t59 * t23 +
        0.2447460819436755E4 * t230 * t65 * t57 * t23 -
        0.1223730409718377E4 * t230 * t59 * t23 -
        0.1223730409718377E4 * t67 * t65 * t59 * t23;
      double t295 = 0.1265055323929465E1 + t39 + t1;
      double t297 = 0.7852315164806451 + t39 + t1;
      double t298 = 0.2147684835193549 + t39 + t1;
      double t300 = -0.2650553239294647 + t39 + t1;
      double t301 = t300 * t298 * t297;
      double t313 = t297 * t295;
      double t328 =
        -0.1223730409718377E4 * t66 * t59 * t23 -
        0.5495800667890619E3 * t180 * t59 * t23 +
        0.708236963546707 * t301 * t295 * t19 -
        0.708236963546707 * t301 * t295 * t5 - 0.1416473927093414E1 * t301 * t23 -
        0.1416473927093414E1 * t300 * t298 * t295 * t23 -
        0.1416473927093414E1 * t300 * t313 * t23 -
        0.1416473927093414E1 * t298 * t313 * t23 - 0.3417426130481145E2 * t19 +
        0.3417426130481145E2 * t5 - 0.3880839126646436E3 * t58 * t80 * t5 +
        0.1432641226109498E4 * t182 * t239 * t5;
      double t331 = t71 + t109 + t144 + t186 + t221 + t256 + t288 + t328;
      return t331;
    }

    double
    ortho2_f42y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t9 = -t3 - t8 - 0.2147684835193549;
      double t10 = -t3 - t8 - 0.7852315164806451;
      double t11 = t10 * t9;
      double t12 = -t3 - t8 - 0.1265055323929465E1;
      double t13 = t12 * t11;
      double t16 = t6 * t2;
      double t17 = 0.1E1 * x;
      double t18 = 0.1154653670707977E1 + t17 + t1;
      double t20 = 0.5 + t17 + t1;
      double t21 = -0.1546536707079771 + t17 + t1;
      double t25 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t26 = t25 * t21 * t20;
      double t29 = 0.1265055323929465E1 + t17 + t1;
      double t31 = 0.7852315164806451 + t17 + t1;
      double t32 = 0.2147684835193549 + t17 + t1;
      double t34 = -0.2650553239294647 + t17 + t1;
      double t35 = t34 * t32 * t31;
      double t44 = -t3 - t8 + 0.2650553239294647;
      double t49 = t31 * t29;
      double t56 = t9 * t44;
      double t62 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t63 = t44 * t62;
      double t68 = t6 * t4;
      double t72 =
        -0.4784957998111317E3 * t13 * t7 +
        0.9296899540645173E2 * t26 * t18 * t16 +
        0.708236963546707 * t35 * t29 * t16 - 0.708236963546707 * t35 * t7 -
        0.708236963546707 * t34 * t32 * t29 * t7 -
        0.4784957998111317E3 * t12 * t10 * t44 * t7 -
        0.708236963546707 * t34 * t49 * t7 - 0.708236963546707 * t32 * t49 * t7 -
        0.4784957998111317E3 * t12 * t56 * t7 +
        0.8297514121605843E3 * t12 * t9 * t63 * t7 -
        0.708236963546707 * t35 * t29 * t68;
      double t82 = t62 * t16;
      double t83 = -t3 - t8 + 0.1546536707079771;
      double t84 = -t3 - t8 - 0.5;
      double t85 = t84 * t83;
      double t86 = -t3 - t8 - 0.1154653670707977E1;
      double t87 = t86 * t85;
      double t90 = 0.9472135954999579 + t17 + t1;
      double t91 = 0.5278640450004206E-1 + t17 + t1;
      double t92 = t91 * t90;
      double t95 = t90 * t6;
      double t98 = t91 * t6;
      double t101 = t90 * t68;
      double t103 = t86 * t84;
      double t104 = t103 * t83 * t91;
      double t111 = t83 * t6;
      double t117 =
        -0.9296899540645173E2 * t26 * t18 * t68 +
        0.8297514121605843E3 * t11 * t63 * t7 -
        0.4784957998111317E3 * t10 * t56 * t7 + 0.2070594623398002E3 * t87 * t82 -
        0.1032532682195698E3 * t92 * t16 + 0.1032532682195698E3 * t95 * t5 +
        0.1032532682195698E3 * t98 * t5 + 0.1223730409718377E4 * t104 * t101 -
        0.1002384876260218E4 * t86 * t84 * t6 * t5 -
        0.1002384876260218E4 * t86 * t111 * t5 + 0.5011924381301092E3 * t87 * t68;
      double t123 = t83 * t62;
      double t127 = -t3 - t8 + 0.3302238962785669;
      double t129 = -t3 - t8 - 0.3115120652928579E-1;
      double t131 = -t3 - t8 - 0.9688487934707142;
      double t132 = -t3 - t8 - 0.1330223896278567E1;
      double t133 = t132 * t131;
      double t134 = t133 * t84 * t129;
      double t142 = t12 * t10;
      double t143 = t142 * t56;
      double t146 = t62 * t68;
      double t155 = -t3 - t8 - 0.5278640450004206E-1;
      double t156 = t155 * t6;
      double t157 = -t3 - t8 - 0.9472135954999579;
      double t161 = t90 * t16;
      double t164 =
        0.4141189246796004E3 * t86 * t84 * t62 * t7 +
        0.4141189246796004E3 * t86 * t123 * t7 +
        0.3502271162385125E3 * t134 * t127 * t68 -
        0.5011924381301092E3 * t87 * t16 +
        0.4141189246796004E3 * t84 * t123 * t7 +
        0.1311952177084291E4 * t143 * t7 - 0.2070594623398002E3 * t87 * t146 -
        0.7761678253292872E3 * t91 * t95 * t5 +
        0.2392478999055658E3 * t13 * t44 * t68 +
        0.7163206130547491E3 * t157 * t156 * t5 -
        0.1223730409718377E4 * t104 * t161;
      double t169 = t157 * t155 * t91;
      double t174 = t157 * t155;
      double t178 = t157 * t155 * t62;
      double t181 = t20 * t18;
      double t182 = t21 * t181;
      double t202 =
        -0.1002384876260218E4 * t84 * t111 * t5 +
        0.5495800667890619E3 * t169 * t101 + 0.1223730409718377E4 * t104 * t7 -
        0.2330316441804572E3 * t174 * t16 - 0.2265204672181151E3 * t178 * t68 +
        0.6383296673055628E1 * t182 * t68 -
        0.3502271162385125E3 * t134 * t127 * t16 -
        0.4660632883609144E3 * t157 * t6 * t5 - 0.4660632883609144E3 * t156 * t5 +
        0.2330316441804572E3 * t174 * t68 - 0.1083469821547055E3 * t25 * t6 * t5 +
        0.6383296673055628E1 * t21 * t20 * t6 * t5;
      double t211 = t18 * t6;
      double t226 = t62 * t6;
      double t229 = t25 * t62;
      double t232 = 0.3417426130481145E2 * t68 - 0.3417426130481145E2 * t16 +
        0.1032532682195698E3 * t92 * t68 - 0.700454232477025E3 * t134 * t7 +
        0.6383296673055628E1 * t21 * t211 * t5 -
        0.5495800667890619E3 * t169 * t161 + 0.5495800667890619E3 * t169 * t7 -
        0.700454232477025E3 * t133 * t84 * t127 * t7 +
        0.6383296673055628E1 * t20 * t211 * t5 +
        0.2166939643094111E3 * t226 * t5 - 0.3426232412144874E2 * t229 * t16;
      double t237 = t25 * t92;
      double t247 = t129 * t127;
      double t263 =
        0.5495800667890619E3 * t157 * t155 * t90 * t7 -
        0.1227229087290139E3 * t237 * t68 - 0.587987554524238E3 * t182 * t7 -
        0.1099160133578124E4 * t157 * t92 * t7 +
        0.2265204672181151E3 * t178 * t16 -
        0.700454232477025E3 * t133 * t247 * t7 +
        0.1223730409718377E4 * t103 * t83 * t90 * t7 - 0.1939087209985343E3 * t7 +
        0.4530409344362302E3 * t157 * t226 * t5 -
        0.9296899540645173E2 * t25 * t181 * t7 + 0.2654035799974342E2 * t82;
      double t293 =
        0.1227229087290139E3 * t237 * t16 -
        0.2447460819436755E4 * t103 * t92 * t7 +
        0.4392983982614612E2 * t25 * t16 - 0.2654035799974342E2 * t146 -
        0.4392983982614612E2 * t25 * t68 + 0.6547795120836261E3 * t87 * t7 +
        0.4148757060802921E3 * t143 * t82 -
        0.2447460819436755E4 * t86 * t83 * t92 * t7 -
        0.700454232477025E3 * t132 * t84 * t247 * t7 -
        0.1227229087290139E3 * t25 * t98 * t5 -
        0.1227229087290139E3 * t25 * t95 * t5;
      double t329 =
        0.8297514121605843E3 * t142 * t9 * t62 * t7 -
        0.700454232477025E3 * t131 * t84 * t247 * t7 -
        0.2392478999055658E3 * t13 * t44 * t16 +
        0.4530409344362302E3 * t155 * t226 * t5 -
        0.2447460819436755E4 * t85 * t92 * t7 -
        0.6383296673055628E1 * t182 * t16 - 0.9296899540645173E2 * t26 * t7 -
        0.1099160133578124E4 * t155 * t92 * t7 -
        0.9296899540645173E2 * t25 * t21 * t18 * t7 +
        0.8297514121605843E3 * t142 * t63 * t7 -
        0.4148757060802921E3 * t143 * t146 + 0.3426232412144874E2 * t229 * t68;
      double t332 = t72 + t117 + t164 + t202 + t232 + t263 + t293 + t329;
      return t332;
    }

    // * f43 *********************************************************************

    double
    ortho2_f43 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = t10 * t5;
      double t12 = 0.1E1 * y;
      double t19 =
        (-t3 - t12 - 0.1265055323929465E1) * (-t3 - t12 -
                0.7852315164806451) * (-t3 - t12 -
                     0.2147684835193549)
        * (-t3 - t12 + 0.2650553239294647);
      double t22 = 0.1E1 * x;
      double t23 = 0.9472135954999579 + t22 + t1;
      double t24 = t23 * t6;
      double t26 = 0.5278640450004206E-1 + t22 + t1;
      double t27 = -t3 - t12 + 0.1546536707079771;
      double t29 = -t3 - t12 - 0.5;
      double t30 = -t3 - t12 - 0.1154653670707977E1;
      double t31 = t30 * t29;
      double t35 = t6 * t5;
      double t38 = -t3 - t12 - 0.5278640450004206E-1;
      double t40 = -t3 - t12 - 0.9472135954999579;
      double t44 = 0.1154653670707977E1 + t22 + t1;
      double t45 = 0.5 + t22 + t1;
      double t46 = t45 * t44;
      double t47 = -0.1546536707079771 + t22 + t1;
      double t50 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t74 = t26 * t23;
      double t75 = t40 * t38;
      double MapleGenVar1 =
        -0.7427114999812399E3 * t19 * t11 +
        0.4476149005487822E4 * t31 * t27 * t26 * t24 * t5 +
        0.1319342091877083E3 * t19 * t35 -
        0.7187648907364357E3 * t40 * t38 * t9 * t35 -
        0.2090949573729648E3 * t50 * t47 * t46 * t35 +
        0.2305189355770199E3 * (-0.2650553239294647 + t22 +
              t1) * (0.2147684835193549 + t22 +
               t1) * (0.7852315164806451 + t22 +
                t1) * (0.1265055323929465E1 + t22 +
                 t1) * t35 +
        0.3155586322274129E3 * t26 * t24 * t5 +
        0.2503489132304689E3 * t40 * t38 * t6 * t5 +
        0.5060803494957414E2 * t50 * t10 * t5;
      double t115 =
        MapleGenVar1 + 0.101640355235855E4 * t75 * t74 * t35 -
        0.1090214776952226E3 * t50 * t6 * t5 - 0.9138384797240655E2 * t11 -
        0.5035075741297269E3 * t50 * t74 * t35 +
        0.3804118857440016E3 * t47 * t46 * t35 +
        0.1069570146684975E4 * t30 * t29 * t27 * t35 -
        0.213034906419206E3 * t31 * t27 * t9 * t35 +
        0.277287385932578E4 * t75 * t47 * t45 * t44 * t6 * t5 +
        0.5634363474712685E2 * t35 + 0.4064062919510089E3 * (-t3 - t12 -
                   0.1330223896278567E1)
        * (-t3 - t12 - 0.9688487934707142) * t29 * (-t3 - t12 -
                0.3115120652928579E-1) *
        (-t3 - t12 + 0.3302238962785669) * t6 * t5;
      return t115;
    }

    double
    ortho2_f43x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * x;
      double t7 = 0.1154653670707977E1 + t6 + t1;
      double t8 = t7 * t5;
      double t9 = 0.5 + t6 + t1;
      double t10 = -0.1546536707079771 + t6 + t1;
      double t11 = t10 * t9;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 - 0.5278640450004206E-1;
      double t14 = -t3 - t12 - 0.9472135954999579;
      double t15 = t14 * t13;
      double t16 = t15 * t11;
      double t20 = (-t3 - t1) * t2;
      double t21 = t4 * t20;
      double t22 = -t3 - t12 + 0.2650553239294647;
      double t23 = -t3 - t12 - 0.2147684835193549;
      double t24 = t23 * t22;
      double t25 = -t3 - t12 - 0.7852315164806451;
      double t26 = -t3 - t12 - 0.1265055323929465E1;
      double t27 = t26 * t25;
      double t28 = t27 * t24;
      double t31 = -t3 - t12 + 0.3302238962785669;
      double t33 = -t3 - t12 - 0.3115120652928579E-1;
      double t34 = -t3 - t12 - 0.5;
      double t36 = -t3 - t12 - 0.9688487934707142;
      double t37 = -t3 - t12 - 0.1330223896278567E1;
      double t38 = t37 * t36;
      double t39 = t38 * t34 * t33;
      double t45 = t9 * t7;
      double t57 = t10 * t7;
      double t63 = t7 * t20;
      double t72 =
        -0.138643692966289E4 * t16 * t8 + 0.2348659984340823E4 * t28 * t21 +
        0.2032031459755044E3 * t39 * t31 * t20 -
        0.2032031459755044E3 * t39 * t31 * t5 -
        0.138643692966289E4 * t13 * t10 * t45 * t21 -
        0.138643692966289E4 * t14 * t10 * t45 * t21 +
        0.277287385932578E4 * t15 * t45 * t21 +
        0.277287385932578E4 * t15 * t57 * t21 + 0.277287385932578E4 * t16 * t21 +
        0.138643692966289E4 * t16 * t63 -
        0.2032031459755044E3 * t38 * t34 * t31 * t21 -
        0.2032031459755044E3 * t39 * t21;
      double t73 = t33 * t31;
      double t86 = t25 * t23;
      double t87 = t26 * t86;
      double t103 = -t3 - t12 - 0.1154653670707977E1;
      double t107 = -t3 - t12 + 0.1546536707079771;
      double t108 = t107 * t4;
      double t112 = t34 * t107;
      double t113 = t103 * t112;
      double t118 =
        -0.2032031459755044E3 * t38 * t73 * t21 -
        0.2032031459755044E3 * t36 * t34 * t73 * t21 -
        0.2032031459755044E3 * t37 * t34 * t73 * t21 +
        0.6596710459385414E2 * t87 * t22 * t20 -
        0.6596710459385414E2 * t87 * t21 -
        0.6596710459385414E2 * t26 * t25 * t22 * t21 -
        0.6596710459385414E2 * t26 * t24 * t21 -
        0.6596710459385414E2 * t25 * t24 * t21 -
        0.5347850733424873E3 * t103 * t34 * t4 * t20 -
        0.5347850733424873E3 * t103 * t108 * t20 -
        0.5347850733424873E3 * t113 * t5 + 0.5347850733424873E3 * t113 * t20;
      double t123 = t10 * t45;
      double t130 = t7 * t4;
      double t139 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t140 = t139 * t5;
      double t143 = t139 * t20;
      double t150 = t22 * t139;
      double t162 =
        -0.5347850733424873E3 * t34 * t108 * t20 +
        0.1902059428720008E3 * t123 * t20 +
        0.3804118857440016E3 * t10 * t9 * t4 * t20 +
        0.3804118857440016E3 * t10 * t130 * t20 +
        0.3804118857440016E3 * t9 * t130 * t20 + 0.37135574999062E3 * t28 * t140 -
        0.37135574999062E3 * t28 * t143 +
        0.37135574999062E3 * t27 * t23 * t139 * t21 +
        0.37135574999062E3 * t27 * t150 * t21 +
        0.37135574999062E3 * t26 * t23 * t150 * t21 +
        0.37135574999062E3 * t86 * t150 * t21 + 0.4569192398620327E2 * t140;
      double t165 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t173 = t14 * t13 * t139;
      double t178 = t139 * t4;
      double t182 = t165 * t11;
      double t194 =
        -0.5451073884761128E2 * t165 * t20 + 0.5451073884761128E2 * t165 * t5 -
        0.4569192398620327E2 * t143 + 0.1166030092462909E3 * t21 +
        0.3593824453682178E3 * t173 * t5 - 0.3593824453682178E3 * t173 * t20 +
        0.3593824453682178E3 * t14 * t178 * t20 +
        0.1045474786864824E3 * t182 * t8 +
        0.3593824453682178E3 * t13 * t178 * t20 -
        0.1045474786864824E3 * t182 * t63 + 0.106517453209603E3 * t113 * t140 -
        0.106517453209603E3 * t113 * t143;
      double t205 = t107 * t139;
      double t212 = 0.9472135954999579 + t6 + t1;
      double t213 = t212 * t5;
      double t214 = 0.5278640450004206E-1 + t6 + t1;
      double t216 = t103 * t34;
      double t217 = t216 * t107 * t214;
      double t220 = t212 * t20;
      double t226 = t14 * t13 * t214;
      double t236 =
        -0.1251744566152345E3 * t15 * t5 +
        0.106517453209603E3 * t103 * t34 * t139 * t21 +
        0.1251744566152345E3 * t15 * t20 +
        0.106517453209603E3 * t103 * t205 * t21 +
        0.106517453209603E3 * t34 * t205 * t21 -
        0.2238074502743911E4 * t217 * t213 + 0.2238074502743911E4 * t217 * t220 +
        0.4476149005487822E4 * t217 * t21 - 0.508201776179275E3 * t226 * t213 +
        0.508201776179275E3 * t226 * t220 + 0.101640355235855E4 * t226 * t21 -
        0.1251744566152345E3 * t14 * t4 * t20;
      double t240 = t13 * t4;
      double t245 = t165 * t139;
      double t254 = t212 * t4;
      double t257 = t214 * t212;
      double t262 = t214 * t4;
      double t272 =
        -0.1600366583460589E3 * t165 * t4 * t20 -
        0.1251744566152345E3 * t240 * t20 + 0.8001832917302944E2 * t178 * t20 -
        0.2530401747478707E2 * t245 * t5 +
        0.101640355235855E4 * t14 * t13 * t212 * t21 +
        0.2530401747478707E2 * t245 * t20 + 0.3155586322274129E3 * t254 * t20 -
        0.1577793161137065E3 * t257 * t5 + 0.1577793161137065E3 * t257 * t20 +
        0.3155586322274129E3 * t262 * t20 -
        0.508201776179275E3 * t14 * t257 * t21 +
        0.4476149005487822E4 * t216 * t107 * t212 * t21;
      double t287 = 0.1265055323929465E1 + t6 + t1;
      double t289 = 0.7852315164806451 + t6 + t1;
      double t290 = 0.2147684835193549 + t6 + t1;
      double t292 = -0.2650553239294647 + t6 + t1;
      double t293 = t292 * t290 * t289;
      double t305 = t289 * t287;
      double t314 =
        -0.2238074502743911E4 * t216 * t257 * t21 -
        0.2238074502743911E4 * t103 * t107 * t257 * t21 -
        0.2238074502743911E4 * t112 * t257 * t21 -
        0.508201776179275E3 * t13 * t257 * t21 -
        0.11525946778851E3 * t293 * t287 * t5 +
        0.11525946778851E3 * t293 * t287 * t20 +
        0.2305189355770199E3 * t293 * t21 +
        0.2305189355770199E3 * t292 * t290 * t287 * t21 +
        0.2305189355770199E3 * t292 * t305 * t21 +
        0.2305189355770199E3 * t290 * t305 * t21 - 0.2817181737356343E2 * t5 +
        0.2817181737356343E2 * t20;
      double t326 = t165 * t257;
      double t349 =
        -0.7961153766980048E3 * t214 * t254 * t20 +
        0.2272934156889197E4 * t14 * t240 * t20 -
        0.3306081562771931E3 * t123 * t21 -
        0.2090949573729648E3 * t165 * t45 * t21 +
        0.2517537870648634E3 * t326 * t5 - 0.2517537870648634E3 * t326 * t20 +
        0.6736755254055165E3 * t113 * t21 -
        0.5035075741297269E3 * t165 * t262 * t20 -
        0.5035075741297269E3 * t165 * t254 * t20 -
        0.6596710459385414E2 * t87 * t22 * t5 - 0.1902059428720008E3 * t123 * t5 -
        0.2090949573729648E3 * t182 * t21 -
        0.2090949573729648E3 * t165 * t57 * t21;
      double t352 = t72 + t118 + t162 + t194 + t236 + t272 + t314 + t349;

      return t352;
    }

    double
    ortho2_f43y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t6 = 0.1E1 * y;
      double t7 = -t1 - t6 + 0.2650553239294647;
      double t9 = -t1 - t6 - 0.2147684835193549;
      double t10 = -t1 - t6 - 0.7852315164806451;
      double t11 = t10 * t9;
      double t12 = -t1 - t6 - 0.1265055323929465E1;
      double t13 = t12 * t11;
      double t16 = 0.5 + t2;
      double t17 = t3 * t16;
      double t18 = t4 * t17;
      double t19 = t9 * t7;
      double t23 = t4 * t16;
      double t25 = -t1 - t6 + 0.3302238962785669;
      double t27 = -t1 - t6 - 0.3115120652928579E-1;
      double t28 = -t1 - t6 - 0.5;
      double t30 = -t1 - t6 - 0.9688487934707142;
      double t31 = -t1 - t6 - 0.1330223896278567E1;
      double t32 = t31 * t30;
      double t33 = t32 * t28 * t27;
      double t37 = -t1 - t6 - 0.1154653670707977E1;
      double t41 = 0.1E1 * x;
      double t42 = 0.1154653670707977E1 + t41 + t2;
      double t43 = 0.5 + t41 + t2;
      double t44 = t43 * t42;
      double t45 = -0.1546536707079771 + t41 + t2;
      double t46 = t45 * t44;
      double t49 = 0.1265055323929465E1 + t41 + t2;
      double t51 = 0.7852315164806451 + t41 + t2;
      double t52 = 0.2147684835193549 + t41 + t2;
      double t54 = -0.2650553239294647 + t41 + t2;
      double t55 = t54 * t52 * t51;
      double t58 = -t1 - t6 + 0.1546536707079771;
      double t59 = t58 * t4;
      double t63 = t12 * t10;
      double t64 = t63 * t19;
      double t67 = 0.9472135954999579 + t41 + t2;
      double t68 = t67 * t23;
      double t69 = 0.5278640450004206E-1 + t41 + t2;
      double t71 = t37 * t28;
      double t72 = t71 * t58 * t69;
      double t79 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t80 = -t1 - t6 - 0.5278640450004206E-1;
      double t82 = -t1 - t6 - 0.9472135954999579;
      double t83 = t82 * t80 * t79;
      double t86 =
        0.6596710459385414E2 * t13 * t7 * t5 -
        0.1319342091877083E3 * t10 * t19 * t18 - 0.2817181737356343E2 * t23 +
        0.2032031459755044E3 * t33 * t25 * t5 -
        0.1069570146684975E4 * t37 * t28 * t4 * t17 +
        0.1902059428720008E3 * t46 * t5 - 0.11525946778851E3 * t55 * t49 * t23 -
        0.1069570146684975E4 * t37 * t59 * t17 +
        0.1174329992170411E4 * t64 * t18 - 0.2238074502743911E4 * t72 * t68 +
        0.11525946778851E3 * t55 * t18 - 0.3593824453682178E3 * t83 * t5;
      double t89 = t28 * t58;
      double t90 = t37 * t89;
      double t97 = t42 * t23;
      double t98 = t45 * t43;
      double t99 = t82 * t80;
      double t100 = t99 * t98;
      double t105 = t51 * t49;
      double t111 = t79 * t4;
      double t121 = t45 * t42;
      double t126 = t82 * t80 * t69;
      double t129 =
        0.3593824453682178E3 * t83 * t23 - 0.5347850733424873E3 * t90 * t23 +
        0.11525946778851E3 * t54 * t52 * t49 * t18 -
        0.138643692966289E4 * t100 * t97 + 0.2238074502743911E4 * t72 * t18 +
        0.11525946778851E3 * t54 * t105 * t18 + 0.138643692966289E4 * t100 * t18 +
        0.7187648907364357E3 * t82 * t111 * t17 +
        0.11525946778851E3 * t52 * t105 * t18 +
        0.7187648907364357E3 * t80 * t111 * t17 +
        0.138643692966289E4 * t99 * t121 * t18 - 0.508201776179275E3 * t126 * t68;
      double t136 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t137 = t136 * t79;
      double t140 = t69 * t67;
      double t143 = t67 * t4;
      double t146 = t69 * t4;
      double t164 = t80 * t4;
      double t170 =
        0.138643692966289E4 * t99 * t44 * t18 + 0.2530401747478707E2 * t137 * t5 -
        0.1577793161137065E3 * t140 * t23 + 0.1577793161137065E3 * t143 * t17 +
        0.1577793161137065E3 * t146 * t17 -
        0.277287385932578E4 * t82 * t45 * t44 * t18 -
        0.277287385932578E4 * t80 * t45 * t44 * t18 +
        0.508201776179275E3 * t126 * t18 - 0.1251744566152345E3 * t99 * t23 -
        0.2503489132304689E3 * t82 * t4 * t17 -
        0.2503489132304689E3 * t164 * t17 -
        0.1069570146684975E4 * t28 * t59 * t17;
      double t175 = t136 * t98;
      double t184 = t67 * t5;
      double t197 = t42 * t4;
      double t210 =
        0.508201776179275E3 * t82 * t80 * t67 * t18 +
        0.1045474786864824E3 * t175 * t97 -
        0.159223075339601E4 * t69 * t143 * t17 +
        0.1136467078444598E4 * t82 * t164 * t17 +
        0.508201776179275E3 * t126 * t184 -
        0.101640355235855E4 * t82 * t140 * t18 +
        0.1902059428720008E3 * t45 * t43 * t4 * t17 -
        0.2032031459755044E3 * t33 * t25 * t23 +
        0.1902059428720008E3 * t45 * t197 * t17 +
        0.2238074502743911E4 * t71 * t58 * t67 * t18 +
        0.1251744566152345E3 * t99 * t5 - 0.8001832917302944E2 * t136 * t4 * t17;
      double t226 = t136 * t140;
      double t238 = t27 * t25;
      double t242 = t79 * t23;
      double t244 =
        0.1600366583460589E3 * t111 * t17 - 0.2530401747478707E2 * t137 * t23 -
        0.4064062919510089E3 * t33 * t18 +
        0.1902059428720008E3 * t43 * t197 * t17 -
        0.4064062919510089E3 * t32 * t28 * t25 * t18 -
        0.2517537870648634E3 * t226 * t5 -
        0.4476149005487822E4 * t71 * t140 * t18 +
        0.11525946778851E3 * t55 * t49 * t5 - 0.6612163125543862E3 * t46 * t18 -
        0.2002656329224652E3 * t18 - 0.4064062919510089E3 * t32 * t238 * t18 +
        0.4569192398620327E2 * t242;
      double t249 = t79 * t5;
      double t253 = t42 * t5;
      double t278 =
        0.5451073884761128E2 * t136 * t23 + 0.37135574999062E3 * t64 * t242 -
        0.4569192398620327E2 * t249 - 0.5451073884761128E2 * t136 * t5 -
        0.1045474786864824E3 * t175 * t253 -
        0.4476149005487822E4 * t37 * t58 * t140 * t18 +
        0.2517537870648634E3 * t226 * t23 -
        0.1045474786864824E3 * t136 * t44 * t18 +
        0.3368377627027583E3 * t90 * t18 -
        0.4064062919510089E3 * t31 * t28 * t238 * t18 +
        0.7427114999812399E3 * t63 * t9 * t79 * t18 -
        0.2517537870648634E3 * t136 * t146 * t17;
      double t282 = t7 * t79;
      double t309 = t58 * t79;
      double t319 =
        0.106517453209603E3 * t90 * t242 +
        0.7427114999812399E3 * t63 * t282 * t18 -
        0.4476149005487822E4 * t89 * t140 * t18 -
        0.2517537870648634E3 * t136 * t143 * t17 -
        0.4064062919510089E3 * t30 * t28 * t238 * t18 +
        0.7427114999812399E3 * t12 * t9 * t282 * t18 -
        0.101640355235855E4 * t80 * t140 * t18 +
        0.138643692966289E4 * t100 * t253 +
        0.213034906419206E3 * t37 * t28 * t79 * t18 +
        0.213034906419206E3 * t37 * t309 * t18 -
        0.6596710459385414E2 * t13 * t7 * t23 +
        0.213034906419206E3 * t28 * t309 * t18;
      double t350 =
        0.7427114999812399E3 * t11 * t282 * t18 -
        0.1902059428720008E3 * t46 * t23 + 0.1577793161137065E3 * t140 * t5 +
        0.2238074502743911E4 * t72 * t184 - 0.1045474786864824E3 * t175 * t18 -
        0.1045474786864824E3 * t136 * t121 * t18 -
        0.37135574999062E3 * t64 * t249 - 0.106517453209603E3 * t90 * t249 -
        0.1319342091877083E3 * t13 * t18 -
        0.1319342091877083E3 * t12 * t10 * t7 * t18 +
        0.5347850733424873E3 * t90 * t5 + 0.2817181737356343E2 * t5 -
        0.1319342091877083E3 * t12 * t19 * t18;
      double t353 = t86 + t129 + t170 + t210 + t244 + t278 + t319 + t350;
      return t353;
    }

    // * f44 *********************************************************************

    double
    ortho2_f44 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = t10 * t5;
      double t12 = 0.1E1 * y;
      double t19 =
        (-t3 - t12 - 0.1265055323929465E1) * (-t3 - t12 -
                0.7852315164806451) * (-t3 - t12 -
                     0.2147684835193549)
        * (-t3 - t12 + 0.2650553239294647);
      double t22 = 0.1E1 * x;
      double t23 = 0.9472135954999579 + t22 + t1;
      double t24 = t23 * t6;
      double t26 = 0.5278640450004206E-1 + t22 + t1;
      double t27 = -t3 - t12 + 0.1546536707079771;
      double t29 = -t3 - t12 - 0.5;
      double t30 = -t3 - t12 - 0.1154653670707977E1;
      double t31 = t30 * t29;
      double t35 = t6 * t5;
      double t38 = -t3 - t12 - 0.5278640450004206E-1;
      double t40 = -t3 - t12 - 0.9472135954999579;
      double t44 = 0.1154653670707977E1 + t22 + t1;
      double t45 = 0.5 + t22 + t1;
      double t46 = t45 * t44;
      double t47 = -0.1546536707079771 + t22 + t1;
      double t50 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t55 = 0.1265055323929465E1 + t22 + t1;
      double t56 = 0.7852315164806451 + t22 + t1;
      double t58 = 0.2147684835193549 + t22 + t1;
      double t59 = -0.2650553239294647 + t22 + t1;
      double t74 = t26 * t23;
      double t75 = t40 * t38;
      double MapleGenVar1 =
        -0.3881050652464114E3 * t19 * t11 +
        0.3714143415349624E4 * t31 * t27 * t26 * t24 * t5 -
        0.7443964996605715E2 * t19 * t35 -
        0.7233599916796882E3 * t40 * t38 * t9 * t35 +
        0.5369046204632096E2 * t50 * t47 * t46 * t35 +
        0.1032641494677201E3 * t59 * t58 * t56 * t55 * t35 -
        0.3475490625224816E2 * t26 * t24 * t5 -
        0.1283721139761255E3 * t40 * t38 * t6 * t5 -
        0.1999834990545062E2 * t50 * t10 * t5 -
        0.4288238606240546E3 * t75 * t74 * t35;
      double t122 =
        MapleGenVar1 - 0.118290851565856E3 * t50 * t6 * t5 -
        0.1203624724285125E3 * t11 - 0.789327558665435E3 * t50 * t74 * t35 +
        0.910294384281492E3 * t47 * t46 * t35 -
        0.7139746658567991E3 * t50 * t59 * t58 * t56 * t55 * t6 * t5 +
        0.7182343948002611E3 * t30 * t29 * t27 * t35 +
        0.9910605273148391E2 * t31 * t27 * t9 * t35 +
        0.4847453898874358E4 * t75 * t47 * t45 * t44 * t6 * t5 -
        0.1641670879611324E2 * t35 + 0.1524238661659677E3 * (-t3 - t12 -
                   0.1330223896278567E1)
        * (-t3 - t12 - 0.9688487934707142) * t29 * (-t3 - t12 -
                0.3115120652928579E-1) *
        (-t3 - t12 + 0.3302238962785669) * t6 * t5;
      return t122;
    }

    double
    ortho2_f44x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.5 + t8 + t1;
      double t10 = -0.1546536707079771 + t8 + t1;
      double t11 = t10 * t9;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 - 0.5278640450004206E-1;
      double t14 = -t3 - t12 - 0.9472135954999579;
      double t15 = t14 * t13;
      double t16 = t15 * t11;
      double t19 = 0.1154653670707977E1 + t8 + t1;
      double t20 = t19 * t5;
      double t23 = t6 * t2;
      double t24 = t19 * t23;
      double t27 = -t3 - t12 + 0.2650553239294647;
      double t28 = -t3 - t12 - 0.2147684835193549;
      double t29 = t28 * t27;
      double t30 = -t3 - t12 - 0.7852315164806451;
      double t31 = -t3 - t12 - 0.1265055323929465E1;
      double t32 = t31 * t30;
      double t33 = t32 * t29;
      double t36 = t10 * t19;
      double t40 = t9 * t19;
      double t48 = -t3 - t12 + 0.3302238962785669;
      double t50 = -t3 - t12 - 0.3115120652928579E-1;
      double t51 = -t3 - t12 - 0.5;
      double t53 = -t3 - t12 - 0.9688487934707142;
      double t54 = -t3 - t12 - 0.1330223896278567E1;
      double t55 = t54 * t53;
      double t56 = t55 * t51 * t50;
      double t59 = 0.1265055323929465E1 + t8 + t1;
      double t60 = 0.7852315164806451 + t8 + t1;
      double t61 = t60 * t59;
      double t62 = 0.2147684835193549 + t8 + t1;
      double t63 = -0.2650553239294647 + t8 + t1;
      double t81 =
        0.4847453898874358E4 * t16 * t7 + 0.2423726949437179E4 * t16 * t20 -
        0.2423726949437179E4 * t16 * t24 + 0.1227295977626918E4 * t33 * t7 +
        0.4847453898874358E4 * t15 * t36 * t7 -
        0.2423726949437179E4 * t14 * t10 * t40 * t7 +
        0.4847453898874358E4 * t15 * t40 * t7 -
        0.7621193308298383E2 * t56 * t48 * t23 -
        0.112889306788257E4 * t63 * t62 * t61 * t7 -
        0.2423726949437179E4 * t13 * t10 * t40 * t7 -
        0.7621193308298383E2 * t56 * t7 + 0.7621193308298383E2 * t56 * t48 * t5 -
        0.7621193308298383E2 * t55 * t51 * t48 * t7;
      double t82 = t50 * t48;
      double t96 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t97 = t96 * t5;
      double t100 = t96 * t23;
      double t103 = t27 * t96;
      double t115 = t30 * t28;
      double t119 = 0.9472135954999579 + t8 + t1;
      double t120 = t119 * t23;
      double t121 = 0.5278640450004206E-1 + t8 + t1;
      double t122 = -t3 - t12 + 0.1546536707079771;
      double t124 = -t3 - t12 - 0.1154653670707977E1;
      double t125 = t124 * t51;
      double t126 = t125 * t122 * t121;
      double t129 = t119 * t5;
      double t138 =
        -0.7621193308298383E2 * t55 * t82 * t7 -
        0.7621193308298383E2 * t53 * t51 * t82 * t7 -
        0.7621193308298383E2 * t54 * t51 * t82 * t7 -
        0.1940525326232057E3 * t33 * t97 + 0.1940525326232057E3 * t33 * t100 +
        0.1940525326232057E3 * t31 * t28 * t103 * t7 +
        0.1940525326232057E3 * t32 * t103 * t7 +
        0.1940525326232057E3 * t32 * t28 * t96 * t7 +
        0.1940525326232057E3 * t115 * t103 * t7 -
        0.1857071707674812E4 * t126 * t120 + 0.1857071707674812E4 * t126 * t129 +
        0.3714143415349624E4 * t126 * t7 +
        0.3714143415349624E4 * t125 * t122 * t119 * t7;
      double t140 = t121 * t119;
      double t144 = t51 * t122;
      double t152 = t59 * t23;
      double t153 = t62 * t60;
      double t156 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t157 = t156 * t63;
      double t158 = t157 * t153;
      double t161 = t59 * t5;
      double t167 = t62 * t59;
      double t182 =
        -0.1857071707674812E4 * t125 * t140 * t7 -
        0.1857071707674812E4 * t144 * t140 * t7 -
        0.1857071707674812E4 * t124 * t122 * t140 * t7 +
        0.3569873329283995E3 * t158 * t152 - 0.3569873329283995E3 * t158 * t161 -
        0.7139746658567991E3 * t157 * t61 * t7 -
        0.7139746658567991E3 * t157 * t167 * t7 -
        0.7139746658567991E3 * t158 * t7 -
        0.7139746658567991E3 * t156 * t62 * t61 * t7 -
        0.5914542578292801E2 * t156 * t5 + 0.6018123621425624E2 * t100 +
        0.8208354398056619E1 * t23 - 0.8208354398056619E1 * t5;
      double t187 = t119 * t6;
      double t191 = t13 * t6;
      double t195 = t10 * t40;
      double t205 = t156 * t140;
      double t216 = 0.1935852990288174E3 * t7 + 0.5914542578292801E2 * t156 * t23
        - 0.6018123621425624E2 * t97 - 0.1248036452661476E4 * t121 * t187 * t5 +
        0.2287465141948263E4 * t14 * t191 * t5 + 0.848920743465995E2 * t195 * t7
        - 0.6418605698806277E2 * t15 * t5 + 0.6418605698806277E2 * t15 * t23 +
        0.5369046204632096E2 * t156 * t40 * t7 +
        0.3946637793327175E3 * t205 * t23 - 0.3946637793327175E3 * t205 * t5 +
        0.6418605698806277E2 * t14 * t6 * t5 +
        0.632403351462369E2 * t156 * t6 * t5;
      double t221 = t156 * t96;
      double t224 = t96 * t6;
      double t231 = t124 * t144;
      double t234 = t121 * t6;
      double t248 = t31 * t115;
      double t253 =
        0.6418605698806277E2 * t191 * t5 + 0.9999174952725308E1 * t221 * t23 -
        0.3162016757311845E2 * t224 * t5 - 0.9999174952725308E1 * t221 * t5 +
        0.1737745312612408E2 * t140 * t23 - 0.313400856540241E3 * t231 * t7 -
        0.3475490625224816E2 * t234 * t5 - 0.3475490625224816E2 * t187 * t5 -
        0.1737745312612408E2 * t140 * t5 -
        0.789327558665435E3 * t156 * t234 * t5 -
        0.789327558665435E3 * t156 * t187 * t5 +
        0.3721982498302858E2 * t248 * t27 * t23 -
        0.455147192140746E3 * t195 * t23;
      double t254 = t156 * t11;
      double t279 = t122 * t6;
      double t292 =
        0.5369046204632096E2 * t254 * t7 +
        0.5369046204632096E2 * t156 * t36 * t7 -
        0.3721982498302858E2 * t248 * t27 * t5 +
        0.3721982498302858E2 * t248 * t7 +
        0.3721982498302858E2 * t31 * t30 * t27 * t7 +
        0.3721982498302858E2 * t31 * t29 * t7 +
        0.3721982498302858E2 * t30 * t29 * t7 -
        0.3591171974001305E3 * t124 * t51 * t6 * t5 -
        0.3591171974001305E3 * t124 * t279 * t5 -
        0.3591171974001305E3 * t231 * t23 + 0.3591171974001305E3 * t231 * t5 -
        0.3591171974001305E3 * t51 * t279 * t5 + 0.455147192140746E3 * t195 * t5;
      double t298 = t19 * t6;
      double t306 = t14 * t13 * t96;
      double t329 = t122 * t96;
      double t333 =
        0.910294384281492E3 * t10 * t9 * t6 * t5 +
        0.910294384281492E3 * t10 * t298 * t5 +
        0.910294384281492E3 * t9 * t298 * t5 + 0.3616799958398441E3 * t306 * t23 -
        0.3616799958398441E3 * t306 * t5 +
        0.3616799958398441E3 * t14 * t224 * t5 +
        0.3616799958398441E3 * t13 * t224 * t5 -
        0.2684523102316048E2 * t254 * t24 + 0.2684523102316048E2 * t254 * t20 -
        0.4955302636574196E2 * t231 * t100 + 0.4955302636574196E2 * t231 * t97 -
        0.4955302636574196E2 * t124 * t51 * t96 * t7 -
        0.4955302636574196E2 * t124 * t329 * t7;
      double t338 = t14 * t13 * t121;
      double t355 = t63 * t153;
      double t371 =
        -0.4955302636574196E2 * t51 * t329 * t7 +
        0.2144119303120273E3 * t338 * t120 - 0.2144119303120273E3 * t338 * t129 -
        0.4288238606240546E3 * t338 * t7 -
        0.4288238606240546E3 * t14 * t13 * t119 * t7 +
        0.2144119303120273E3 * t14 * t140 * t7 +
        0.2144119303120273E3 * t13 * t140 * t7 -
        0.5163207473386007E2 * t355 * t152 + 0.5163207473386007E2 * t355 * t161 +
        0.1032641494677201E3 * t355 * t7 +
        0.1032641494677201E3 * t63 * t167 * t7 +
        0.1032641494677201E3 * t63 * t61 * t7 +
        0.1032641494677201E3 * t62 * t61 * t7;
      double t374 = t81 + t138 + t182 + t216 + t253 + t292 + t333 + t371;

      return t374;
    }

    double
    ortho2_f44y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * y;
      double t8 = -t3 - t7 + 0.1546536707079771;
      double t9 = t8 * t6;
      double t10 = -t3 - t7 - 0.1154653670707977E1;
      double t16 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t17 = t16 * t6;
      double t18 = -t3 - t7 - 0.9472135954999579;
      double t22 = 0.1E1 * x;
      double t23 = 0.1154653670707977E1 + t22 + t1;
      double t24 = t23 * t6;
      double t25 = -0.1546536707079771 + t22 + t1;
      double t29 = -t3 - t7 - 0.5;
      double t34 = t6 * t2;
      double t35 = -t3 - t7 - 0.5278640450004206E-1;
      double t37 = t18 * t35 * t16;
      double t40 = t6 * t5;
      double t41 = -t3 - t7 + 0.2650553239294647;
      double t42 = t41 * t16;
      double t43 = -t3 - t7 - 0.7852315164806451;
      double t44 = -t3 - t7 - 0.1265055323929465E1;
      double t45 = t44 * t43;
      double t49 = t23 * t34;
      double t50 = 0.5 + t22 + t1;
      double t51 = t25 * t50;
      double t54 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t55 = t54 * t51;
      double t58 = t6 * t4;
      double t59 = t50 * t23;
      double t60 = t25 * t59;
      double t63 = t54 * t16;
      double t66 = 0.9472135954999579 + t22 + t1;
      double t67 = t66 * t6;
      double t70 = 0.5278640450004206E-1 + t22 + t1;
      double t71 = t70 * t6;
      double t74 = t16 * t34;
      double t75 = -t3 - t7 - 0.2147684835193549;
      double t76 = t75 * t41;
      double t77 = t45 * t76;
      double t80 = t70 * t66;
      double t83 =
        -0.7182343948002611E3 * t10 * t9 * t5 +
        0.7233599916796882E3 * t18 * t17 * t5 +
        0.455147192140746E3 * t25 * t24 * t5 -
        0.7182343948002611E3 * t10 * t29 * t6 * t5 +
        0.3616799958398441E3 * t37 * t34 +
        0.3881050652464114E3 * t45 * t42 * t40 -
        0.2684523102316048E2 * t55 * t49 + 0.455147192140746E3 * t60 * t58 -
        0.9999174952725308E1 * t63 * t58 - 0.1737745312612408E2 * t67 * t5 -
        0.1737745312612408E2 * t71 * t5 + 0.1940525326232057E3 * t77 * t74 +
        0.1737745312612408E2 * t80 * t34;
      double t84 = t29 * t8;
      double t85 = t10 * t84;
      double t89 = t43 * t75;
      double t90 = t44 * t89;
      double t93 = t54 * t80;
      double t96 = 0.1265055323929465E1 + t22 + t1;
      double t97 = t96 * t58;
      double t98 = 0.7852315164806451 + t22 + t1;
      double t99 = 0.2147684835193549 + t22 + t1;
      double t100 = t99 * t98;
      double t101 = -0.2650553239294647 + t22 + t1;
      double t102 = t101 * t100;
      double t105 = t23 * t58;
      double t127 = t8 * t16;
      double t134 =
        -0.3591171974001305E3 * t85 * t34 -
        0.3721982498302858E2 * t90 * t41 * t58 -
        0.3946637793327175E3 * t93 * t58 + 0.5163207473386007E2 * t102 * t97 +
        0.2684523102316048E2 * t55 * t105 +
        0.7443964996605715E2 * t43 * t76 * t40 +
        0.455147192140746E3 * t50 * t24 * t5 - 0.4955302636574196E2 * t85 * t74 +
        0.3881050652464114E3 * t45 * t75 * t16 * t40 -
        0.7182343948002611E3 * t29 * t9 * t5 -
        0.9910605273148391E2 * t10 * t29 * t16 * t40 -
        0.9910605273148391E2 * t10 * t127 * t40 +
        0.7233599916796882E3 * t35 * t17 * t5;
      double t150 = t16 * t58;
      double t153 = t66 * t34;
      double t155 = t10 * t29;
      double t156 = t155 * t8 * t70;
      double t168 = t18 * t35 * t70;
      double t179 =
        -0.9910605273148391E2 * t29 * t127 * t40 +
        0.7443964996605715E2 * t44 * t76 * t40 +
        0.455147192140746E3 * t25 * t50 * t6 * t5 +
        0.3881050652464114E3 * t44 * t75 * t42 * t40 +
        0.4955302636574196E2 * t85 * t150 - 0.1857071707674812E4 * t156 * t153 +
        0.7443964996605715E2 * t44 * t43 * t41 * t40 +
        0.7443964996605715E2 * t90 * t40 + 0.1857071707674812E4 * t156 * t40 +
        0.2144119303120273E3 * t168 * t153 - 0.2144119303120273E3 * t168 * t40 -
        0.1737745312612408E2 * t80 * t58 -
        0.2144119303120273E3 * t18 * t35 * t66 * t40;
      double t183 = t66 * t58;
      double t192 = t54 * t101;
      double t193 = t192 * t100;
      double t196 = -t3 - t7 + 0.3302238962785669;
      double t198 = -t3 - t7 - 0.3115120652928579E-1;
      double t200 = -t3 - t7 - 0.9688487934707142;
      double t201 = -t3 - t7 - 0.1330223896278567E1;
      double t202 = t201 * t200;
      double t203 = t202 * t29 * t198;
      double t223 =
        0.4288238606240546E3 * t18 * t80 * t40 +
        0.1857071707674812E4 * t156 * t183 + 0.3591171974001305E3 * t85 * t58 +
        0.1857071707674812E4 * t155 * t8 * t66 * t40 -
        0.3569873329283995E3 * t193 * t97 +
        0.7621193308298383E2 * t203 * t196 * t58 -
        0.3714143415349624E4 * t155 * t80 * t40 -
        0.3714143415349624E4 * t10 * t8 * t80 * t40 - 0.8208354398056619E1 * t58 +
        0.3881050652464114E3 * t89 * t42 * t40 +
        0.6136479888134591E3 * t77 * t40 -
        0.3714143415349624E4 * t84 * t80 * t40 + 0.8208354398056619E1 * t34;
      double t229 = t18 * t35;
      double t230 = t229 * t51;
      double t235 = t25 * t23;
      double t252 = t96 * t34;
      double t255 = t98 * t96;
      double t264 =
        0.4288238606240546E3 * t35 * t80 * t40 -
        0.2423726949437179E4 * t230 * t49 + 0.2423726949437179E4 * t230 * t40 +
        0.2423726949437179E4 * t229 * t235 * t40 +
        0.2423726949437179E4 * t229 * t59 * t40 -
        0.1940525326232057E3 * t77 * t150 -
        0.4847453898874358E4 * t18 * t25 * t59 * t40 -
        0.4847453898874358E4 * t35 * t25 * t59 * t40 -
        0.5163207473386007E2 * t102 * t252 -
        0.2257786135765139E4 * t101 * t99 * t255 * t40 -
        0.1837587384673415E3 * t40 + 0.6018123621425624E2 * t74 -
        0.5914542578292801E2 * t54 * t58;
      double t270 = t35 * t6;
      double t279 = t99 * t96;
      double t298 =
        -0.2496072905322951E4 * t70 * t67 * t5 +
        0.5163207473386007E2 * t102 * t40 +
        0.1143732570974131E4 * t18 * t270 * t5 -
        0.2144119303120273E3 * t168 * t183 + 0.5914542578292801E2 * t54 * t34 -
        0.6018123621425624E2 * t150 + 0.5163207473386007E2 * t101 * t279 * t40 -
        0.3616799958398441E3 * t37 * t58 +
        0.5163207473386007E2 * t101 * t255 * t40 -
        0.7621193308298383E2 * t203 * t196 * t34 +
        0.1283721139761255E3 * t270 * t5 - 0.6418605698806277E2 * t229 * t58 +
        0.3162016757311845E2 * t54 * t6 * t5;
      double t313 = t198 * t196;
      double t334 =
        -0.632403351462369E2 * t17 * t5 - 0.1524238661659677E3 * t203 * t40 -
        0.1524238661659677E3 * t202 * t29 * t196 * t40 +
        0.5163207473386007E2 * t99 * t255 * t40 +
        0.169784148693199E3 * t60 * t40 -
        0.1524238661659677E3 * t202 * t313 * t40 +
        0.3569873329283995E3 * t193 * t252 +
        0.2684523102316048E2 * t54 * t59 * t40 -
        0.3569873329283995E3 * t193 * t40 -
        0.3569873329283995E3 * t192 * t279 * t40 +
        0.3946637793327175E3 * t93 * t34 -
        0.3569873329283995E3 * t192 * t255 * t40 -
        0.1567004282701205E3 * t85 * t40;
      double t372 =
        -0.1524238661659677E3 * t201 * t29 * t313 * t40 -
        0.3569873329283995E3 * t54 * t99 * t255 * t40 +
        0.9999174952725308E1 * t63 * t34 - 0.3946637793327175E3 * t54 * t71 * t5 -
        0.3946637793327175E3 * t54 * t67 * t5 +
        0.6418605698806277E2 * t229 * t34 + 0.1283721139761255E3 * t18 * t6 * t5 -
        0.1524238661659677E3 * t200 * t29 * t313 * t40 +
        0.2423726949437179E4 * t230 * t105 +
        0.3721982498302858E2 * t90 * t41 * t34 - 0.455147192140746E3 * t60 * t34 +
        0.2684523102316048E2 * t55 * t40 +
        0.2684523102316048E2 * t54 * t235 * t40;
      double t375 = t83 + t134 + t179 + t223 + t264 + t298 + t334 + t372;
      return t375;
    }

    // * f45 *********************************************************************

    double
    ortho2_f45 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = t10 * t5;
      double t12 = 0.1E1 * y;
      double t19 =
        (-t3 - t12 - 0.1265055323929465E1) * (-t3 - t12 -
                0.7852315164806451) * (-t3 - t12 -
                     0.2147684835193549)
        * (-t3 - t12 + 0.2650553239294647);
      double t22 = 0.1E1 * x;
      double t23 = 0.9472135954999579 + t22 + t1;
      double t24 = t23 * t6;
      double t26 = 0.5278640450004206E-1 + t22 + t1;
      double t27 = -t3 - t12 + 0.1546536707079771;
      double t29 = -t3 - t12 - 0.5;
      double t30 = -t3 - t12 - 0.1154653670707977E1;
      double t31 = t30 * t29;
      double t35 = t6 * t5;
      double t38 = -t3 - t12 - 0.5278640450004206E-1;
      double t40 = -t3 - t12 - 0.9472135954999579;
      double t48 = 0.5 + t22 + t1;
      double t56 = 0.1154653670707977E1 + t22 + t1;
      double t57 = t48 * t56;
      double t58 = -0.1546536707079771 + t22 + t1;
      double t61 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t66 = 0.1265055323929465E1 + t22 + t1;
      double t67 = 0.7852315164806451 + t22 + t1;
      double t69 = 0.2147684835193549 + t22 + t1;
      double t70 = -0.2650553239294647 + t22 + t1;
      double t86 = t26 * t23;
      double t87 = t40 * t38;
      double t134 =
        -0.1288167888880889E4 * t87 * t86 * t35 -
        0.9969787481987889E2 * t61 * t6 * t5 - 0.13907695463633E3 * t11 -
        0.7479865155940853E3 * t61 * t86 * t35 +
        0.1474763450130956E4 * t58 * t57 * t35 -
        0.1054064757613177E4 * t61 * t70 * t69 * t67 * t66 * t6 * t5 +
        0.310282078438805E3 * t30 * t29 * t27 * t35 +
        0.1477324329797148E3 * t31 * t27 * t9 * t35 +
        0.3349716319662684E4 * t87 * t58 * t48 * t56 * t6 * t5 -
        0.8777090031616546E2 * t35 + 0.3975476104111114E2 * (-t3 - t12 -
                   0.1330223896278567E1)
        * (-t3 - t12 - 0.9688487934707142) * t29 * (-t3 - t12 -
                0.3115120652928579E-1) *
        (-t3 - t12 + 0.3302238962785669) * t6 * t5;
      double t135 =
        -0.1212686992021287E3 * t19 * t11 +
        0.1608229607943578E4 * t31 * t27 * t26 * t24 * t5 -
        0.6077066017834362E2 * t19 * t35 -
        0.4446959185129003E3 * t40 * t38 * t9 * t35 +
        0.1242949520089937E4 * (-0.3302238962785669 + t22 +
              t1) * (0.3115120652928579E-1 + t22 +
               t1) * t48 * (0.9688487934707142 + t22 +
                t1) * (0.1330223896278567E1 +
                       t22 + t1) * t6 * t5 +
        0.5018233792718522E3 * t61 * t58 * t57 * t35 -
        0.7544556267199322E3 * t70 * t69 * t67 * t66 * t35 -
        0.6198258990089902E3 * t26 * t24 * t5 -
        0.2579039585178655E3 * t40 * t38 * t6 * t5 -
        0.7328687740830976E2 * t61 * t10 * t5 + t134;
      return t135;
    }

    double
    ortho2_f45x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * x;
      double t7 = 0.1154653670707977E1 + t6 + t1;
      double t8 = t7 * t5;
      double t9 = 0.5 + t6 + t1;
      double t10 = -0.1546536707079771 + t6 + t1;
      double t11 = t10 * t9;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 - 0.5278640450004206E-1;
      double t14 = -t3 - t12 - 0.9472135954999579;
      double t15 = t14 * t13;
      double t16 = t15 * t11;
      double t20 = (-t3 - t1) * t2;
      double t21 = t4 * t20;
      double t22 = -t3 - t12 + 0.2650553239294647;
      double t23 = -t3 - t12 - 0.2147684835193549;
      double t24 = t23 * t22;
      double t25 = -t3 - t12 - 0.7852315164806451;
      double t26 = -t3 - t12 - 0.1265055323929465E1;
      double t27 = t26 * t25;
      double t28 = t27 * t24;
      double t31 = t7 * t20;
      double t36 = t9 * t7;
      double t44 = t10 * t7;
      double t48 = 0.1265055323929465E1 + t6 + t1;
      double t49 = 0.7852315164806451 + t6 + t1;
      double t50 = t49 * t48;
      double t51 = 0.2147684835193549 + t6 + t1;
      double t52 = -0.2650553239294647 + t6 + t1;
      double t61 = -t3 - t12 + 0.3302238962785669;
      double t63 = -t3 - t12 - 0.3115120652928579E-1;
      double t64 = -t3 - t12 - 0.5;
      double t66 = -t3 - t12 - 0.9688487934707142;
      double t67 = -t3 - t12 - 0.1330223896278567E1;
      double t68 = t67 * t66;
      double t69 = t68 * t64 * t63;
      double t81 =
        -0.1674858159831342E4 * t16 * t8 + 0.3834852983645706E3 * t28 * t21 +
        0.1674858159831342E4 * t16 * t31 + 0.3349716319662684E4 * t16 * t21 -
        0.1674858159831342E4 * t14 * t10 * t36 * t21 +
        0.3349716319662684E4 * t15 * t36 * t21 +
        0.3349716319662684E4 * t15 * t44 * t21 -
        0.1666622717685473E4 * t52 * t51 * t50 * t21 -
        0.1674858159831342E4 * t13 * t10 * t36 * t21 -
        0.1987738052055557E2 * t69 * t61 * t5 +
        0.1987738052055557E2 * t69 * t61 * t20 -
        0.1987738052055557E2 * t68 * t64 * t61 * t21 -
        0.1987738052055557E2 * t69 * t21;
      double t82 = t63 * t61;
      double t86 = 0.1330223896278567E1 + t6 + t1;
      double t88 = 0.9688487934707142 + t6 + t1;
      double t90 = 0.3115120652928579E-1 + t6 + t1;
      double t91 = -0.3302238962785669 + t6 + t1;
      double t92 = t91 * t90;
      double t93 = t92 * t9 * t88;
      double t105 = t88 * t86;
      double t127 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t128 = t22 * t127;
      double t136 = t127 * t20;
      double t139 = t127 * t5;
      double t142 =
        -0.1987738052055557E2 * t68 * t82 * t21 -
        0.6214747600449685E3 * t93 * t86 * t5 +
        0.1242949520089937E4 * t92 * t9 * t86 * t21 +
        0.1242949520089937E4 * t93 * t21 +
        0.6214747600449685E3 * t93 * t86 * t20 +
        0.1242949520089937E4 * t92 * t105 * t21 +
        0.1242949520089937E4 * t91 * t9 * t105 * t21 +
        0.1242949520089937E4 * t90 * t9 * t105 * t21 -
        0.1987738052055557E2 * t66 * t64 * t82 * t21 -
        0.1987738052055557E2 * t67 * t64 * t82 * t21 +
        0.6063434960106435E2 * t27 * t128 * t21 +
        0.6063434960106435E2 * t27 * t23 * t127 * t21 -
        0.6063434960106435E2 * t28 * t136 + 0.6063434960106435E2 * t28 * t139;
      double t148 = t25 * t23;
      double t153 = t14 * t13 * t127;
      double t158 = t127 * t4;
      double t167 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t168 = t167 * t11;
      double t173 = -t3 - t12 + 0.1546536707079771;
      double t174 = t64 * t173;
      double t175 = -t3 - t12 - 0.1154653670707977E1;
      double t176 = t175 * t174;
      double t187 = t173 * t127;
      double t191 =
        0.6063434960106435E2 * t26 * t23 * t128 * t21 +
        0.6063434960106435E2 * t148 * t128 * t21 +
        0.2223479592564502E3 * t153 * t5 - 0.2223479592564502E3 * t153 * t20 +
        0.2223479592564502E3 * t14 * t158 * t20 +
        0.2223479592564502E3 * t13 * t158 * t20 -
        0.2509116896359261E3 * t168 * t8 + 0.2509116896359261E3 * t168 * t31 -
        0.738662164898574E2 * t176 * t139 + 0.738662164898574E2 * t176 * t136 +
        0.2821637655361341E3 * t21 + 0.6953847731816499E2 * t139 -
        0.738662164898574E2 * t175 * t64 * t127 * t21 -
        0.738662164898574E2 * t175 * t187 * t21;
      double t195 = 0.9472135954999579 + t6 + t1;
      double t196 = t195 * t5;
      double t197 = 0.5278640450004206E-1 + t6 + t1;
      double t199 = t175 * t64;
      double t200 = t199 * t173 * t197;
      double t208 = t195 * t20;
      double t214 = t14 * t13 * t197;
      double t225 = t197 * t195;
      double t236 =
        -0.738662164898574E2 * t64 * t187 * t21 -
        0.8041148039717888E3 * t200 * t196 - 0.6953847731816499E2 * t136 -
        0.4984893740993945E2 * t167 * t20 + 0.4984893740993945E2 * t167 * t5 +
        0.8041148039717888E3 * t200 * t208 + 0.1608229607943578E4 * t200 * t21 +
        0.6440839444404447E3 * t214 * t196 - 0.6440839444404447E3 * t214 * t208 -
        0.1288167888880889E4 * t214 * t21 -
        0.1288167888880889E4 * t14 * t13 * t195 * t21 +
        0.6440839444404447E3 * t14 * t225 * t21 +
        0.1608229607943578E4 * t199 * t173 * t195 * t21 -
        0.8041148039717888E3 * t199 * t225 * t21;
      double t253 = t48 * t5;
      double t254 = t51 * t49;
      double t255 = t52 * t254;
      double t261 = t48 * t20;
      double t264 = t13 * t4;
      double t267 = t167 * t127;
      double t277 = t195 * t4;
      double t280 =
        -0.1289519792589328E3 * t15 * t20 + 0.1289519792589328E3 * t15 * t5 -
        0.8041148039717888E3 * t175 * t173 * t225 * t21 -
        0.8041148039717888E3 * t174 * t225 * t21 +
        0.6440839444404447E3 * t13 * t225 * t21 +
        0.3772278133599661E3 * t255 * t253 +
        0.1289519792589328E3 * t14 * t4 * t20 -
        0.3772278133599661E3 * t255 * t261 + 0.1289519792589328E3 * t264 * t20 +
        0.3664343870415488E2 * t267 * t5 - 0.1158767276058983E3 * t158 * t20 +
        0.2317534552117966E3 * t167 * t4 * t20 -
        0.3664343870415488E2 * t267 * t20 - 0.6198258990089902E3 * t277 * t20;
      double t287 = t197 * t4;
      double t290 = t51 * t48;
      double t300 = t167 * t52;
      double t301 = t300 * t254;
      double t319 =
        0.3099129495044951E3 * t225 * t5 - 0.3099129495044951E3 * t225 * t20 -
        0.7544556267199322E3 * t255 * t21 - 0.6198258990089902E3 * t287 * t20 -
        0.7544556267199322E3 * t52 * t290 * t21 -
        0.7544556267199322E3 * t52 * t50 * t21 -
        0.7544556267199322E3 * t51 * t50 * t21 +
        0.5270323788065883E3 * t301 * t253 - 0.5270323788065883E3 * t301 * t261 -
        0.1054064757613177E4 * t300 * t290 * t21 -
        0.1054064757613177E4 * t301 * t21 -
        0.1054064757613177E4 * t300 * t50 * t21 -
        0.1054064757613177E4 * t167 * t51 * t50 * t21 + 0.4388545015808273E2 * t5;
      double t328 = t10 * t36;
      double t334 = t167 * t225;
      double t348 = t26 * t148;
      double t358 =
        -0.4388545015808273E2 * t20 - 0.1182670524185181E4 * t197 * t277 * t20 +
        0.1406251968681403E4 * t14 * t264 * t20 +
        0.793452430810791E3 * t328 * t21 +
        0.5018233792718522E3 * t167 * t36 * t21 +
        0.3739932577970426E3 * t334 * t5 - 0.3739932577970426E3 * t334 * t20 -
        0.4671709724940744E3 * t176 * t21 -
        0.7479865155940853E3 * t167 * t287 * t20 -
        0.7479865155940853E3 * t167 * t277 * t20 +
        0.3038533008917181E2 * t348 * t22 * t5 - 0.737381725065478E3 * t328 * t5 +
        0.5018233792718522E3 * t168 * t21 +
        0.5018233792718522E3 * t167 * t44 * t21;
      double t378 = t173 * t4;
      double t395 = t7 * t4;
      double t402 =
        -0.3038533008917181E2 * t348 * t22 * t20 +
        0.3038533008917181E2 * t348 * t21 +
        0.3038533008917181E2 * t26 * t25 * t22 * t21 +
        0.3038533008917181E2 * t26 * t24 * t21 +
        0.3038533008917181E2 * t25 * t24 * t21 -
        0.1551410392194025E3 * t175 * t64 * t4 * t20 -
        0.1551410392194025E3 * t175 * t378 * t20 -
        0.1551410392194025E3 * t176 * t5 + 0.1551410392194025E3 * t176 * t20 -
        0.1551410392194025E3 * t64 * t378 * t20 +
        0.737381725065478E3 * t328 * t20 +
        0.1474763450130956E4 * t10 * t9 * t4 * t20 +
        0.1474763450130956E4 * t10 * t395 * t20 +
        0.1474763450130956E4 * t9 * t395 * t20;
      double t405 = t81 + t142 + t191 + t236 + t280 + t319 + t358 + t402;
      return t405;
    }

    double
    ortho2_f45y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t6 = 0.1E1 * y;
      double t7 = -t1 - t6 - 0.5278640450004206E-1;
      double t8 = -t1 - t6 - 0.9472135954999579;
      double t9 = t8 * t7;
      double t12 = 0.5 + t2;
      double t13 = t3 * t12;
      double t16 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t22 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t23 = t22 * t4;
      double t26 = t4 * t13;
      double t27 = -t1 - t6 + 0.2650553239294647;
      double t28 = -t1 - t6 - 0.7852315164806451;
      double t30 = -t1 - t6 - 0.1265055323929465E1;
      double t34 = -t1 - t6 - 0.2147684835193549;
      double t35 = t34 * t27;
      double t42 = -t1 - t6 - 0.5;
      double t44 = -t1 - t6 - 0.1154653670707977E1;
      double t48 = -t1 - t6 + 0.1546536707079771;
      double t49 = t48 * t4;
      double t53 = 0.1E1 * x;
      double t54 = 0.1265055323929465E1 + t53 + t2;
      double t55 = 0.7852315164806451 + t53 + t2;
      double t56 = t55 * t54;
      double t57 = 0.2147684835193549 + t53 + t2;
      double t58 = -0.2650553239294647 + t53 + t2;
      double t63 = t16 * t22;
      double t66 = t4 * t12;
      double t67 = 0.9472135954999579 + t53 + t2;
      double t68 = 0.5278640450004206E-1 + t53 + t2;
      double t69 = t68 * t67;
      double t72 = t42 * t48;
      double t73 = t44 * t72;
      double t76 = t67 * t4;
      double t80 =
        -0.1289519792589328E3 * t9 * t5 + 0.1158767276058983E3 * t16 * t4 * t13 -
        0.2317534552117966E3 * t23 * t13 +
        0.6077066017834362E2 * t30 * t28 * t27 * t26 +
        0.6077066017834362E2 * t30 * t35 * t26 +
        0.6077066017834362E2 * t28 * t35 * t26 -
        0.310282078438805E3 * t44 * t42 * t4 * t13 -
        0.310282078438805E3 * t44 * t49 * t13 -
        0.3333245435370946E4 * t58 * t57 * t56 * t26 -
        0.3664343870415488E2 * t63 * t5 + 0.3099129495044951E3 * t69 * t66 -
        0.1551410392194025E3 * t73 * t66 - 0.2365341048370363E4 * t68 * t76 * t13;
      double t81 = t7 * t4;
      double t90 = t67 * t5;
      double t92 = t8 * t7 * t68;
      double t95 = 0.5 + t53 + t2;
      double t97 = -0.1546536707079771 + t53 + t2;
      double t101 = 0.1154653670707977E1 + t53 + t2;
      double t102 = t101 * t4;
      double t111 = t95 * t101;
      double t116 = t16 * t69;
      double t120 = t8 * t7 * t22;
      double t127 = t22 * t66;
      double t128 = t30 * t28;
      double t129 = t128 * t35;
      double t132 = -t1 - t6 + 0.3302238962785669;
      double t134 = -t1 - t6 - 0.3115120652928579E-1;
      double t136 = -t1 - t6 - 0.9688487934707142;
      double t137 = -t1 - t6 - 0.1330223896278567E1;
      double t138 = t137 * t136;
      double t139 = t138 * t42 * t134;
      double t142 =
        0.7031259843407014E3 * t8 * t81 * t13 + 0.2579039585178655E3 * t81 * t13 -
        0.310282078438805E3 * t42 * t49 * t13 - 0.6440839444404447E3 * t92 * t90 +
        0.737381725065478E3 * t97 * t95 * t4 * t13 +
        0.737381725065478E3 * t97 * t102 * t13 - 0.3099129495044951E3 * t69 * t5 +
        0.737381725065478E3 * t95 * t102 * t13 -
        0.3349716319662684E4 * t7 * t97 * t111 * t26 -
        0.3739932577970426E3 * t116 * t5 - 0.2223479592564502E3 * t120 * t5 -
        0.3349716319662684E4 * t8 * t97 * t111 * t26 +
        0.6063434960106435E2 * t129 * t127 -
        0.1987738052055557E2 * t139 * t132 * t66;
      double t148 = t27 * t22;
      double t156 = t28 * t34;
      double t162 = 0.1330223896278567E1 + t53 + t2;
      double t164 = 0.9688487934707142 + t53 + t2;
      double t166 = 0.3115120652928579E-1 + t53 + t2;
      double t167 = -0.3302238962785669 + t53 + t2;
      double t168 = t167 * t166;
      double t169 = t168 * t95 * t164;
      double t174 = t68 * t4;
      double t178 = t30 * t156;
      double t181 = t97 * t111;
      double t196 =
        0.1212686992021287E3 * t128 * t34 * t22 * t26 +
        0.1212686992021287E3 * t128 * t148 * t26 +
        0.1212686992021287E3 * t30 * t34 * t148 * t26 +
        0.1212686992021287E3 * t156 * t148 * t26 -
        0.3975476104111114E2 * t139 * t26 +
        0.6214747600449685E3 * t169 * t162 * t5 -
        0.3099129495044951E3 * t76 * t13 - 0.3099129495044951E3 * t174 * t13 -
        0.3038533008917181E2 * t178 * t27 * t5 + 0.737381725065478E3 * t181 * t5 -
        0.3975476104111114E2 * t138 * t42 * t132 * t26 +
        0.2223479592564502E3 * t120 * t66 +
        0.4446959185129003E3 * t8 * t23 * t13 +
        0.4446959185129003E3 * t7 * t23 * t13;
      double t197 = t101 * t66;
      double t198 = t97 * t95;
      double t199 = t16 * t198;
      double t202 = t54 * t5;
      double t203 = t57 * t55;
      double t204 = t58 * t203;
      double t207 = t101 * t5;
      double t211 = t44 * t42;
      double t212 = t211 * t48 * t68;
      double t221 = t48 * t22;
      double t232 = t22 * t5;
      double t235 = t67 * t66;
      double t242 =
        -0.2509116896359261E3 * t199 * t197 - 0.3772278133599661E3 * t204 * t202 +
        0.2509116896359261E3 * t199 * t207 + 0.8041148039717888E3 * t212 * t90 -
        0.738662164898574E2 * t73 * t127 -
        0.1477324329797148E3 * t44 * t42 * t22 * t26 -
        0.1477324329797148E3 * t44 * t221 * t26 +
        0.1551410392194025E3 * t73 * t5 + 0.1586904861621582E4 * t181 * t26 -
        0.1477324329797148E3 * t42 * t221 * t26 +
        0.738662164898574E2 * t73 * t232 - 0.8041148039717888E3 * t212 * t235 +
        0.6077066017834362E2 * t178 * t26 + 0.8041148039717888E3 * t212 * t26;
      double t245 = t134 * t132;
      double t260 = t16 * t58;
      double t261 = t260 * t203;
      double t287 = t164 * t162;
      double t291 =
        -0.3975476104111114E2 * t138 * t245 * t26 +
        0.6440839444404447E3 * t92 * t235 - 0.6440839444404447E3 * t92 * t26 +
        0.2509116896359261E3 * t16 * t111 * t26 -
        0.6440839444404447E3 * t8 * t7 * t67 * t26 -
        0.5270323788065883E3 * t261 * t202 +
        0.1288167888880889E4 * t8 * t69 * t26 +
        0.8041148039717888E3 * t211 * t48 * t67 * t26 -
        0.6214747600449685E3 * t169 * t162 * t66 -
        0.1608229607943578E4 * t211 * t69 * t26 +
        0.6214747600449685E3 * t169 * t26 +
        0.6214747600449685E3 * t168 * t95 * t162 * t26 -
        0.1608229607943578E4 * t44 * t48 * t69 * t26 +
        0.6214747600449685E3 * t168 * t287 * t26;
      double t312 = t54 * t66;
      double t325 =
        -0.1608229607943578E4 * t72 * t69 * t26 +
        0.6214747600449685E3 * t167 * t95 * t287 * t26 +
        0.1288167888880889E4 * t7 * t69 * t26 +
        0.6214747600449685E3 * t166 * t95 * t287 * t26 -
        0.6063434960106435E2 * t129 * t232 - 0.9537238896380791E2 * t26 +
        0.6953847731816499E2 * t127 + 0.3739932577970426E3 * t116 * t66 +
        0.3772278133599661E3 * t204 * t312 - 0.6953847731816499E2 * t232 -
        0.3772278133599661E3 * t204 * t26 +
        0.1987738052055557E2 * t139 * t132 * t5 +
        0.4984893740993945E2 * t16 * t66 - 0.4984893740993945E2 * t16 * t5;
      double t329 = t57 * t54;
      double t367 =
        -0.2335854862470372E3 * t73 * t26 -
        0.3772278133599661E3 * t58 * t329 * t26 -
        0.3772278133599661E3 * t58 * t56 * t26 -
        0.3772278133599661E3 * t57 * t56 * t26 +
        0.3664343870415488E2 * t63 * t66 -
        0.3975476104111114E2 * t137 * t42 * t245 * t26 +
        0.5270323788065883E3 * t261 * t312 - 0.5270323788065883E3 * t261 * t26 -
        0.5270323788065883E3 * t260 * t329 * t26 -
        0.5270323788065883E3 * t260 * t56 * t26 -
        0.5270323788065883E3 * t16 * t57 * t56 * t26 -
        0.3739932577970426E3 * t16 * t174 * t13 -
        0.3739932577970426E3 * t16 * t76 * t13 +
        0.1917426491822853E3 * t129 * t26;
      double t372 = t9 * t198;
      double t388 = t97 * t101;
      double t403 =
        -0.3975476104111114E2 * t136 * t42 * t245 * t26 +
        0.1674858159831342E4 * t372 * t207 +
        0.3038533008917181E2 * t178 * t27 * t66 - 0.4388545015808273E2 * t5 +
        0.4388545015808273E2 * t66 - 0.1674858159831342E4 * t372 * t197 +
        0.1674858159831342E4 * t372 * t26 - 0.737381725065478E3 * t181 * t66 +
        0.2509116896359261E3 * t199 * t26 +
        0.2509116896359261E3 * t16 * t388 * t26 +
        0.1674858159831342E4 * t9 * t388 * t26 +
        0.2579039585178655E3 * t8 * t4 * t13 +
        0.1674858159831342E4 * t9 * t111 * t26 + 0.1289519792589328E3 * t9 * t66;
      double t406 = t80 + t142 + t196 + t242 + t291 + t325 + t367 + t403;
      return t406;
    }

    // ORDER 9

    // Edge functions, order 9

    // number 46
    inline double
    ortho2_f46_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l2 * l3 * phi7 (l3 - l2);
    }

    inline double
    ortho2_f46_1 (double x, double y)
    {
      return -ortho2_f46_0(x, y);
    }

    inline double
    ortho2_f46x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l2x * l3 + l2 * l3x) * phi7 (l3 - l2) + l2 * l3 * phi7x (l3 -
                       l2) *
        (l3x - l2x);
    }

    inline double
    ortho2_f46x_1 (double x, double y)
    {
      return -ortho2_f46x_0(x, y);
    }

    inline double
    ortho2_f46y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l2y * l3 + l2 * l3y) * phi7 (l3 - l2) + l2 * l3 * phi7x (l3 -
                       l2) *
        (l3y - l2y);
    }

    inline double
    ortho2_f46y_1 (double x, double y)
    {
      return -ortho2_f46y_0(x, y);
    }

    // number 47
    inline double
    ortho2_f47_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l3 * l1 * phi7 (l1 - l3);
    }

    inline double
    ortho2_f47_1 (double x, double y)
    {
      return -ortho2_f47_0(x, y);
    }

    inline double
    ortho2_f47x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l3x * l1 + l3 * l1x) * phi7 (l1 - l3) + l3 * l1 * phi7x (l1 -
                       l3) *
        (l1x - l3x);
    }

    inline double
    ortho2_f47x_1 (double x, double y)
    {
      return -ortho2_f47x_0(x, y);
    }

    inline double
    ortho2_f47y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l3y * l1 + l3 * l1y) * phi7 (l1 - l3) + l3 * l1 * phi7x (l1 -
                       l3) *
        (l1y - l3y);
    }

    inline double
    ortho2_f47y_1 (double x, double y)
    {
      return -ortho2_f47y_0(x, y);
    }

    // number 48
    inline double
    ortho2_f48_0 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l1 * l2 * phi7 (l2 - l1);
    }

    inline double
    ortho2_f48_1 (double x, double y)
    {
      return -ortho2_f48_0(x, y);
    }

    inline double
    ortho2_f48x_0 (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l1x * l2 + l1 * l2x) * phi7 (l2 - l1) + l1 * l2 * phi7x (l2 -
                       l1) *
        (l2x - l1x);
    }

    inline double
    ortho2_f48x_1 (double x, double y)
    {
      return -ortho2_f48x_0(x, y);
    }

    inline double
    ortho2_f48y_0 (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l1y * l2 + l1 * l2y) * phi7 (l2 - l1) + l1 * l2 * phi7x (l2 -
                       l1) *
        (l2y - l1y);
    }

    inline double
    ortho2_f48y_1 (double x, double y)
    {
      return -ortho2_f48y_0(x, y);
    }

    // Bubble functions, order 9

    // number 49
    // * f49 *********************************************************************

    double
    ortho2_f49 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = t10 * t5;
      double t12 = 0.1E1 * y;
      double t19 =
        (-t3 - t12 - 0.1265055323929465E1) * (-t3 - t12 -
                0.7852315164806451) * (-t3 - t12 -
                     0.2147684835193549)
        * (-t3 - t12 + 0.2650553239294647);
      double t22 = 0.1E1 * x;
      double t23 = 0.9472135954999579 + t22 + t1;
      double t24 = t23 * t6;
      double t26 = 0.5278640450004206E-1 + t22 + t1;
      double t27 = -t3 - t12 + 0.1546536707079771;
      double t29 = -t3 - t12 - 0.5;
      double t30 = -t3 - t12 - 0.1154653670707977E1;
      double t31 = t30 * t29;
      double t35 = t6 * t5;
      double t38 = -t3 - t12 - 0.5278640450004206E-1;
      double t40 = -t3 - t12 - 0.9472135954999579;
      double t48 = 0.5 + t22 + t1;
      double t70 = 0.1154653670707977E1 + t22 + t1;
      double t71 = t48 * t70;
      double t72 = -0.1546536707079771 + t22 + t1;
      double t75 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t80 = 0.1265055323929465E1 + t22 + t1;
      double t81 = 0.7852315164806451 + t22 + t1;
      double t83 = 0.2147684835193549 + t22 + t1;
      double t84 = -0.2650553239294647 + t22 + t1;
      double MapleGenVar1 =
        -0.3523053470186945E3 * t19 * t11 -
        0.3550572246380408E2 * t31 * t27 * t26 * t24 * t5 +
        0.7488021844033388E3 * t19 * t35 -
        0.1453285249712096E3 * t40 * t38 * t9 * t35 -
        0.5443986854044022 * (-0.3302238962785669 + t22 +
            t1) * (0.3115120652928579E-1 + t22 +
             t1) * t48 * (0.9688487934707142 + t22 +
                    t1) * (0.1330223896278567E1 +
                     t22 + t1) * t6 * t5;
      double t99 =
        MapleGenVar1 + 0.105971340598455E4 * (-t3 - t12 -
                0.1371740148509607E1) * (-t3 - t12 -
                       0.1091700181433142E1)
        * (-t3 - t12 - 0.7092992179024789) * (-t3 - t12 -
                0.2907007820975211) * (-t3 - t12 +
                     0.917001814331423E-1)
        * (-t3 - t12 + 0.3717401485096066) * t6 * t5 +
        0.8329302166120903 * t75 * t72 * t71 * t35 -
        0.1193848984509831E1 * t84 * t83 * t81 * t80 * t35 +
        0.3746701779680405E2 * t26 * t24 * t5 +
        0.2819810433360698E3 * t40 * t38 * t6 * t5 +
        0.1112756411871235E2 * t75 * t10 * t5;
      double t100 = t26 * t23;
      double t101 = t40 * t38;
      double t148 =
        0.2560564644646672E3 * t101 * t100 * t35 -
        0.2705062432907437E2 * t75 * t6 * t5 - 0.1571530000316444E2 * t11 -
        0.2513479686250514 * t75 * t100 * t35 +
        0.1999736322485727E1 * t72 * t71 * t35 +
        0.197969751926724E1 * t75 * t84 * t83 * t81 * t80 * t6 * t5 +
        0.4045851542419182E3 * t30 * t29 * t27 * t35 -
        0.9865198545674906E2 * t31 * t27 * t9 * t35 -
        0.3008401330350343E2 * t101 * t72 * t48 * t70 * t6 * t5 +
        0.2910785337424417E2 * t35 + 0.6466598509071278E3 * (-t3 - t12 -
                   0.1330223896278567E1)
        * (-t3 - t12 - 0.9688487934707142) * t29 * (-t3 - t12 -
                0.3115120652928579E-1) *
        (-t3 - t12 + 0.3302238962785669) * t6 * t5;
      double t149 = t99 + t148;
      return t149;
    }

    double
    ortho2_f49x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t8 = (-t3 - t1) * t2;
      double t10 = 0.1E1 * x;
      double t11 = 0.9472135954999579 + t10 + t1;
      double t12 = 0.5278640450004206E-1 + t10 + t1;
      double t13 = t12 * t11;
      double t16 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t17 = t16 * t13;
      double t22 = t4 * t8;
      double t23 = 0.1E1 * y;
      double t24 = -t3 - t23 + 0.1546536707079771;
      double t25 = -t3 - t23 - 0.5;
      double t26 = t25 * t24;
      double t27 = -t3 - t23 - 0.1154653670707977E1;
      double t28 = t27 * t26;
      double t31 = t12 * t4;
      double t35 = t11 * t4;
      double t39 = -t3 - t23 + 0.2650553239294647;
      double t41 = -t3 - t23 - 0.2147684835193549;
      double t42 = -t3 - t23 - 0.7852315164806451;
      double t43 = t42 * t41;
      double t44 = -t3 - t23 - 0.1265055323929465E1;
      double t45 = t44 * t43;
      double t48 = 0.1330223896278567E1 + t10 + t1;
      double t50 = 0.9688487934707142 + t10 + t1;
      double t51 = 0.5 + t10 + t1;
      double t53 = 0.3115120652928579E-1 + t10 + t1;
      double t54 = -0.3302238962785669 + t10 + t1;
      double t55 = t54 * t53;
      double t56 = t55 * t51 * t50;
      double t59 = -t3 - t23 + 0.3302238962785669;
      double t60 = -t3 - t23 - 0.3115120652928579E-1;
      double t61 = t60 * t59;
      double t62 = -t3 - t23 - 0.9688487934707142;
      double t63 = -t3 - t23 - 0.1330223896278567E1;
      double t64 = t63 * t62;
      double t73 = t64 * t25 * t60;
      double t82 =
        -0.1455392668712209E2 * t5 + 0.1455392668712209E2 * t8 +
        0.1256739843125257 * t17 * t5 - 0.1256739843125257 * t17 * t8 +
        0.3119649697411334E3 * t28 * t22 - 0.2513479686250514 * t16 * t31 * t8 -
        0.2513479686250514 * t16 * t35 * t8 -
        0.3744010922016694E3 * t45 * t39 * t5 +
        0.2721993427022011 * t56 * t48 * t5 -
        0.3233299254535639E3 * t64 * t61 * t22 -
        0.3233299254535639E3 * t64 * t25 * t59 * t22 -
        0.3233299254535639E3 * t73 * t22 + 0.3233299254535639E3 * t73 * t59 * t8 -
        0.3233299254535639E3 * t73 * t59 * t5;
      double t83 = 0.1265055323929465E1 + t10 + t1;
      double t84 = 0.7852315164806451 + t10 + t1;
      double t85 = t84 * t83;
      double t86 = 0.2147684835193549 + t10 + t1;
      double t87 = -0.2650553239294647 + t10 + t1;
      double t92 = 0.1154653670707977E1 + t10 + t1;
      double t93 = t51 * t92;
      double t94 = -0.1546536707079771 + t10 + t1;
      double t95 = -t3 - t23 - 0.5278640450004206E-1;
      double t100 = -t3 - t23 - 0.9472135954999579;
      double t105 = t100 * t95;
      double t109 = t94 * t92;
      double t113 = t94 * t51;
      double t114 = t105 * t113;
      double t117 = t92 * t8;
      double t120 = t92 * t5;
      double t123 = t41 * t39;
      double t124 = t44 * t42;
      double t125 = t124 * t123;
      double t137 = t50 * t48;
      double t149 =
        0.3130176619534776E1 * t87 * t86 * t85 * t22 +
        0.1504200665175171E2 * t95 * t94 * t93 * t22 +
        0.1504200665175171E2 * t100 * t94 * t93 * t22 -
        0.3008401330350343E2 * t105 * t93 * t22 -
        0.3008401330350343E2 * t105 * t109 * t22 -
        0.3008401330350343E2 * t114 * t22 - 0.1504200665175171E2 * t114 * t117 +
        0.1504200665175171E2 * t114 * t120 + 0.1114087328435086E4 * t125 * t22 -
        0.5443986854044022 * t55 * t51 * t48 * t22 -
        0.5443986854044022 * t56 * t22 - 0.2721993427022011 * t56 * t48 * t8 -
        0.5443986854044022 * t55 * t137 * t22 -
        0.5443986854044022 * t54 * t51 * t137 * t22 -
        0.5443986854044022 * t53 * t51 * t137 * t22;
      double t157 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t158 = t157 * t5;
      double t165 = t39 * t157;
      double t180 = t157 * t8;
      double t183 = t11 * t5;
      double t185 = t27 * t25;
      double t186 = t185 * t24 * t12;
      double t189 = t11 * t8;
      double t201 = t94 * t93;
      double t208 =
        -0.3233299254535639E3 * t63 * t25 * t61 * t22 +
        0.1761526735093473E3 * t125 * t158 -
        0.3233299254535639E3 * t62 * t25 * t61 * t22 +
        0.1761526735093473E3 * t43 * t165 * t22 +
        0.1761526735093473E3 * t44 * t41 * t165 * t22 +
        0.1761526735093473E3 * t124 * t165 * t22 +
        0.1761526735093473E3 * t124 * t41 * t157 * t22 -
        0.1761526735093473E3 * t125 * t180 + 0.1775286123190204E2 * t186 * t183 -
        0.1775286123190204E2 * t186 * t189 - 0.3550572246380408E2 * t186 * t22 +
        0.1775286123190204E2 * t185 * t13 * t22 -
        0.3550572246380408E2 * t185 * t24 * t11 * t22 -
        0.9998681612428635 * t201 * t5 +
        0.1775286123190204E2 * t27 * t24 * t13 * t22;
      double t212 = t16 * t113;
      double t218 = t83 * t5;
      double t219 = t86 * t84;
      double t220 = t16 * t87;
      double t221 = t220 * t219;
      double t227 = t86 * t83;
      double t233 = t83 * t8;
      double t259 =
        0.1775286123190204E2 * t26 * t13 * t22 + 0.8329302166120903 * t212 * t22 +
        0.8329302166120903 * t16 * t109 * t22 - 0.9898487596336198 * t221 * t218 +
        0.197969751926724E1 * t220 * t85 * t22 +
        0.197969751926724E1 * t220 * t227 * t22 +
        0.197969751926724E1 * t221 * t22 + 0.9898487596336198 * t221 * t233 +
        0.197969751926724E1 * t16 * t86 * t85 * t22 +
        0.3744010922016694E3 * t45 * t39 * t8 - 0.3744010922016694E3 * t45 * t22 -
        0.3744010922016694E3 * t44 * t42 * t39 * t22 -
        0.3744010922016694E3 * t44 * t123 * t22 -
        0.3744010922016694E3 * t42 * t123 * t22 -
        0.2022925771209591E3 * t27 * t25 * t4 * t8;
      double t262 = t24 * t4;
      double t273 = -t3 - t23 + 0.3717401485096066;
      double t274 = -t3 - t23 + 0.917001814331423E-1;
      double t275 = t274 * t273;
      double t277 = -t3 - t23 - 0.2907007820975211;
      double t278 = -t3 - t23 - 0.7092992179024789;
      double t280 = -t3 - t23 - 0.1091700181433142E1;
      double t281 = -t3 - t23 - 0.1371740148509607E1;
      double t282 = t281 * t280;
      double t283 = t282 * t278 * t277;
      double t308 = t273 * t4 * t8;
      double t311 = t16 * t157;
      double t314 =
        -0.2022925771209591E3 * t27 * t262 * t8 -
        0.2022925771209591E3 * t28 * t5 + 0.2022925771209591E3 * t28 * t8 -
        0.2022925771209591E3 * t25 * t262 * t8 +
        0.5298567029922751E3 * t283 * t275 * t8 -
        0.5298567029922751E3 * t283 * t275 * t5 -
        0.5298567029922751E3 * t283 * t274 * t4 * t8 +
        0.9998681612428635 * t201 * t8 +
        0.1999736322485727E1 * t94 * t51 * t4 * t8 +
        0.1873350889840203E2 * t13 * t8 + 0.3746701779680405E2 * t31 * t8 +
        0.3746701779680405E2 * t35 * t8 - 0.1873350889840203E2 * t13 * t5 -
        0.5298567029922751E3 * t283 * t308 + 0.5563782059356174E1 * t311 * t8;
      double t317 = t157 * t4;
      double t323 = t95 * t4;
      double t333 = t277 * t274;
      double t349 = t92 * t4;
      double t358 =
        -0.5563782059356174E1 * t311 * t5 + 0.1759422371234765E2 * t317 * t8 -
        0.351884474246953E2 * t16 * t4 * t8 - 0.1409905216680349E3 * t323 * t8 -
        0.1409905216680349E3 * t100 * t4 * t8 + 0.1409905216680349E3 * t105 * t8 -
        0.1409905216680349E3 * t105 * t5 -
        0.5298567029922751E3 * t281 * t278 * t333 * t308 -
        0.5298567029922751E3 * t282 * t333 * t308 -
        0.5298567029922751E3 * t282 * t278 * t274 * t308 -
        0.5298567029922751E3 * t280 * t278 * t333 * t308 +
        0.1999736322485727E1 * t94 * t349 * t8 +
        0.1999736322485727E1 * t51 * t349 * t8 + 0.692534961813138E1 * t22 -
        0.7857650001582218E1 * t180;
      double t366 = t100 * t95 * t157;
      double t389 = t24 * t157;
      double t397 = t100 * t95 * t12;
      double t400 =
        -0.1352531216453719E2 * t16 * t8 + 0.1352531216453719E2 * t16 * t5 +
        0.7857650001582218E1 * t158 + 0.7266426248560478E2 * t366 * t5 -
        0.7266426248560478E2 * t366 * t8 +
        0.7266426248560478E2 * t100 * t317 * t8 +
        0.7266426248560478E2 * t95 * t317 * t8 -
        0.4164651083060452 * t212 * t120 + 0.4164651083060452 * t212 * t117 +
        0.4932599272837453E2 * t28 * t158 - 0.4932599272837453E2 * t28 * t180 +
        0.4932599272837453E2 * t27 * t25 * t157 * t22 +
        0.4932599272837453E2 * t27 * t389 * t22 +
        0.4932599272837453E2 * t25 * t389 * t22 -
        0.1280282322323336E3 * t397 * t183;
      double t415 = t87 * t219;
      double t442 =
        0.1280282322323336E3 * t397 * t189 + 0.2560564644646672E3 * t397 * t22 +
        0.2560564644646672E3 * t100 * t95 * t11 * t22 -
        0.1280282322323336E3 * t100 * t13 * t22 -
        0.1280282322323336E3 * t95 * t13 * t22 +
        0.5969244922549153 * t415 * t218 - 0.5969244922549153 * t415 * t233 -
        0.1193848984509831E1 * t415 * t22 -
        0.1193848984509831E1 * t87 * t227 * t22 -
        0.1193848984509831E1 * t87 * t85 * t22 -
        0.1193848984509831E1 * t86 * t85 * t22 -
        0.3974160330558515 * t12 * t35 * t8 +
        0.4595691479016784E3 * t100 * t323 * t8 +
        0.1316978308235811E1 * t201 * t22 + 0.8329302166120903 * t16 * t93 * t22;
      double t445 = t82 + t149 + t208 + t259 + t314 + t358 + t400 + t442;
      return t445;
    }

    double
    ortho2_f49y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t9 = -t3 - t8 + 0.2650553239294647;
      double t10 = -t3 - t8 - 0.7852315164806451;
      double t12 = -t3 - t8 - 0.1265055323929465E1;
      double t18 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t19 = t9 * t18;
      double t20 = -t3 - t8 - 0.2147684835193549;
      double t25 = t6 * t2;
      double t26 = 0.1E1 * x;
      double t27 = 0.1265055323929465E1 + t26 + t1;
      double t28 = t27 * t25;
      double t29 = 0.7852315164806451 + t26 + t1;
      double t30 = 0.2147684835193549 + t26 + t1;
      double t31 = t30 * t29;
      double t32 = -0.2650553239294647 + t26 + t1;
      double t33 = t32 * t31;
      double t36 = t6 * t4;
      double t39 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t42 = t18 * t25;
      double t44 = 0.9472135954999579 + t26 + t1;
      double t45 = t44 * t25;
      double t46 = 0.5278640450004206E-1 + t26 + t1;
      double t47 = -t3 - t8 + 0.1546536707079771;
      double t49 = -t3 - t8 - 0.5;
      double t50 = -t3 - t8 - 0.1154653670707977E1;
      double t51 = t50 * t49;
      double t52 = t51 * t47 * t46;
      double t55 = -t3 - t8 + 0.3717401485096066;
      double t57 = t55 * t6 * t5;
      double t58 = -t3 - t8 + 0.917001814331423E-1;
      double t59 = -t3 - t8 - 0.7092992179024789;
      double t61 = -t3 - t8 - 0.1091700181433142E1;
      double t62 = -t3 - t8 - 0.1371740148509607E1;
      double t63 = t62 * t61;
      double t68 = t18 * t36;
      double t72 = 0.1330223896278567E1 + t26 + t1;
      double t74 = 0.9688487934707142 + t26 + t1;
      double t75 = 0.5 + t26 + t1;
      double t77 = 0.3115120652928579E-1 + t26 + t1;
      double t78 = -0.3302238962785669 + t26 + t1;
      double t79 = t78 * t77;
      double t80 = t79 * t75 * t74;
      double t83 = t39 * t18;
      double t88 = -t3 - t8 - 0.2907007820975211;
      double t90 = t63 * t59 * t88;
      double t95 =
        -0.7488021844033388E3 * t12 * t10 * t9 * t7 +
        0.3523053470186945E3 * t12 * t20 * t19 * t7 +
        0.5969244922549153 * t33 * t28 - 0.1352531216453719E2 * t39 * t36 +
        0.7857650001582218E1 * t42 + 0.1775286123190204E2 * t52 * t45 -
        0.105971340598455E4 * t63 * t59 * t58 * t57 - 0.6069351394801367E2 * t7 -
        0.7857650001582218E1 * t68 + 0.1352531216453719E2 * t39 * t25 +
        0.2721993427022011 * t80 * t72 * t25 - 0.5563782059356174E1 * t83 * t25 -
        0.105971340598455E4 * t90 * t58 * t6 * t5 - 0.2721993427022011 * t80 * t7;
      double t98 = 0.1154653670707977E1 + t26 + t1;
      double t99 = t75 * t98;
      double t100 = -0.1546536707079771 + t26 + t1;
      double t101 = -t3 - t8 - 0.5278640450004206E-1;
      double t110 = -t3 - t8 - 0.9472135954999579;
      double t116 = t10 * t20;
      double t120 = t18 * t6;
      double t124 = t100 * t99;
      double t128 = t110 * t101 * t46;
      double t131 = t74 * t72;
      double t135 = t30 * t27;
      double t139 = t20 * t9;
      double t148 = t110 * t101;
      double t151 =
        -0.1775286123190204E2 * t52 * t7 +
        0.3008401330350343E2 * t101 * t100 * t99 * t7 -
        0.2721993427022011 * t79 * t75 * t72 * t7 -
        0.2819810433360698E3 * t110 * t6 * t5 - 0.5969244922549153 * t33 * t7 +
        0.3523053470186945E3 * t116 * t19 * t7 +
        0.1453285249712096E3 * t110 * t120 * t5 +
        0.2633956616471622E1 * t124 * t7 - 0.1280282322323336E3 * t128 * t45 -
        0.2721993427022011 * t79 * t131 * t7 -
        0.5969244922549153 * t32 * t135 * t7 -
        0.7488021844033388E3 * t12 * t139 * t7 + 0.351884474246953E2 * t120 * t5 -
        0.1759422371234765E2 * t39 * t6 * t5 + 0.1409905216680349E3 * t148 * t36;
      double t153 = t101 * t6;
      double t158 = t46 * t6;
      double t161 = t44 * t6;
      double t164 = t46 * t44;
      double t167 = t27 * t36;
      double t168 = t39 * t32;
      double t169 = t168 * t31;
      double t174 = t29 * t27;
      double t193 = t39 * t164;
      double t204 =
        -0.2819810433360698E3 * t153 * t5 + 0.1280282322323336E3 * t128 * t7 +
        0.1873350889840203E2 * t158 * t5 + 0.1873350889840203E2 * t161 * t5 -
        0.1873350889840203E2 * t164 * t25 + 0.9898487596336198 * t169 * t167 +
        0.5563782059356174E1 * t83 * t36 - 0.5969244922549153 * t32 * t174 * t7 +
        0.9998681612428635 * t100 * t75 * t6 * t5 +
        0.1453285249712096E3 * t101 * t120 * t5 -
        0.2721993427022011 * t78 * t75 * t131 * t7 +
        0.1280282322323336E3 * t110 * t101 * t44 * t7 +
        0.1256739843125257 * t193 * t25 -
        0.2721993427022011 * t77 * t75 * t131 * t7 +
        0.6260353239069551E1 * t32 * t30 * t174 * t7;
      double t205 = t98 * t6;
      double t224 = t49 * t47;
      double t225 = t50 * t224;
      double t231 = t98 * t25;
      double t232 = t100 * t75;
      double t233 = t39 * t232;
      double t241 = -t3 - t8 + 0.3302238962785669;
      double t242 = -t3 - t8 - 0.3115120652928579E-1;
      double t243 = t242 * t241;
      double t244 = -t3 - t8 - 0.1330223896278567E1;
      double t259 =
        0.9998681612428635 * t100 * t205 * t5 -
        0.7488021844033388E3 * t10 * t139 * t7 -
        0.2721993427022011 * t80 * t72 * t36 +
        0.9998681612428635 * t75 * t205 * t5 -
        0.7948320661117029 * t46 * t161 * t5 -
        0.5969244922549153 * t30 * t174 * t7 + 0.1559824848705667E3 * t225 * t7 +
        0.2297845739508392E3 * t110 * t153 * t5 -
        0.4164651083060452 * t233 * t231 -
        0.2560564644646672E3 * t110 * t164 * t7 -
        0.9898487596336198 * t169 * t28 -
        0.6466598509071278E3 * t244 * t49 * t243 * t7 -
        0.4045851542419182E3 * t50 * t49 * t6 * t5 -
        0.1775286123190204E2 * t51 * t47 * t44 * t7 +
        0.9898487596336198 * t169 * t7;
      double t267 = t44 * t36;
      double t280 = -t3 - t8 - 0.9688487934707142;
      double t281 = t244 * t280;
      double t282 = t281 * t49 * t242;
      double t287 = t47 * t6;
      double t294 = t98 * t36;
      double t300 = t12 * t10;
      double t301 = t300 * t139;
      double t310 =
        -0.1256739843125257 * t39 * t158 * t5 - 0.5969244922549153 * t33 * t167 +
        0.1280282322323336E3 * t128 * t267 -
        0.1256739843125257 * t39 * t161 * t5 - 0.1256739843125257 * t193 * t36 +
        0.9898487596336198 * t168 * t135 * t7 +
        0.3233299254535639E3 * t282 * t241 * t36 +
        0.1873350889840203E2 * t164 * t36 -
        0.4045851542419182E3 * t50 * t287 * t5 +
        0.3550572246380408E2 * t51 * t164 * t7 +
        0.4164651083060452 * t233 * t294 -
        0.6466598509071278E3 * t281 * t243 * t7 +
        0.1761526735093473E3 * t301 * t42 - 0.6466598509071278E3 * t282 * t7 +
        0.3550572246380408E2 * t50 * t47 * t164 * t7;
      double t320 = t148 * t232;
      double t329 = t12 * t116;
      double t334 = t88 * t58;
      double t354 =
        -0.2022925771209591E3 * t225 * t25 -
        0.6466598509071278E3 * t280 * t49 * t243 * t7 +
        0.3550572246380408E2 * t224 * t164 * t7 -
        0.1504200665175171E2 * t320 * t294 + 0.4932599272837453E2 * t225 * t42 +
        0.9898487596336198 * t168 * t174 * t7 -
        0.3744010922016694E3 * t329 * t9 * t25 -
        0.1775286123190204E2 * t52 * t267 -
        0.105971340598455E4 * t63 * t334 * t57 +
        0.3523053470186945E3 * t300 * t20 * t18 * t7 +
        0.4164651083060452 * t39 * t99 * t7 - 0.9998681612428635 * t124 * t25 +
        0.9865198545674906E2 * t50 * t49 * t18 * t7 + 0.1455392668712209E2 * t36 +
        0.5570436642175431E3 * t301 * t7;
      double t367 = t110 * t101 * t18;
      double t378 = t100 * t98;
      double t386 = t47 * t18;
      double t401 =
        -0.1409905216680349E3 * t148 * t25 +
        0.9898487596336198 * t39 * t30 * t174 * t7 -
        0.105971340598455E4 * t62 * t59 * t334 * t57 -
        0.7266426248560478E2 * t367 * t36 + 0.4164651083060452 * t233 * t7 -
        0.1455392668712209E2 * t25 + 0.2022925771209591E3 * t225 * t36 -
        0.2560564644646672E3 * t101 * t164 * t7 +
        0.4164651083060452 * t39 * t378 * t7 -
        0.6466598509071278E3 * t281 * t49 * t241 * t7 +
        0.9865198545674906E2 * t50 * t386 * t7 +
        0.1504200665175171E2 * t320 * t231 +
        0.3744010922016694E3 * t329 * t9 * t36 -
        0.1504200665175171E2 * t320 * t7 -
        0.105971340598455E4 * t61 * t59 * t334 * t57;
      double t408 = t58 * t55;
      double t443 =
        0.3523053470186945E3 * t300 * t19 * t7 -
        0.1504200665175171E2 * t148 * t378 * t7 +
        0.5298567029922751E3 * t90 * t408 * t36 +
        0.9865198545674906E2 * t49 * t386 * t7 -
        0.4045851542419182E3 * t49 * t287 * t5 -
        0.1504200665175171E2 * t148 * t99 * t7 +
        0.3008401330350343E2 * t110 * t100 * t99 * t7 -
        0.5298567029922751E3 * t90 * t408 * t25 -
        0.105971340598455E4 * t90 * t57 + 0.9998681612428635 * t124 * t36 -
        0.1761526735093473E3 * t301 * t68 -
        0.3233299254535639E3 * t282 * t241 * t25 -
        0.7488021844033388E3 * t329 * t7 + 0.7266426248560478E2 * t367 * t25 -
        0.4932599272837453E2 * t225 * t68;
      double t446 = t95 + t151 + t204 + t259 + t310 + t354 + t401 + t443;
      return t446;
    }

    // * f50 *********************************************************************

    double
    ortho2_f50 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = t10 * t5;
      double t12 = 0.1E1 * y;
      double t19 =
        (-t3 - t12 - 0.1265055323929465E1) * (-t3 - t12 -
                0.7852315164806451) * (-t3 - t12 -
                     0.2147684835193549)
        * (-t3 - t12 + 0.2650553239294647);
      double t22 = 0.1E1 * x;
      double t23 = 0.9472135954999579 + t22 + t1;
      double t24 = t23 * t6;
      double t26 = 0.5278640450004206E-1 + t22 + t1;
      double t27 = -t3 - t12 + 0.1546536707079771;
      double t29 = -t3 - t12 - 0.5;
      double t30 = -t3 - t12 - 0.1154653670707977E1;
      double t31 = t30 * t29;
      double t35 = t6 * t5;
      double t38 = -t3 - t12 - 0.5278640450004206E-1;
      double t40 = -t3 - t12 - 0.9472135954999579;
      double t48 = 0.5 + t22 + t1;
      double t70 = 0.1154653670707977E1 + t22 + t1;
      double t71 = t48 * t70;
      double t72 = -0.1546536707079771 + t22 + t1;
      double t75 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t80 = 0.1265055323929465E1 + t22 + t1;
      double t81 = 0.7852315164806451 + t22 + t1;
      double t83 = 0.2147684835193549 + t22 + t1;
      double t84 = -0.2650553239294647 + t22 + t1;
      double MapleGenVar1 =
        -0.8889061863602187E3 * t19 * t11 +
        0.2433663888227475E4 * t31 * t27 * t26 * t24 * t5 +
        0.1758998262238863E4 * t19 * t35 -
        0.4666763467016183E3 * t40 * t38 * t9 * t35 -
        0.4148207212977512E1 * (-0.3302238962785669 + t22 +
              t1) * (0.3115120652928579E-1 + t22 +
               t1) * t48 * (0.9688487934707142 + t22 +
                t1) * (0.1330223896278567E1 +
                       t22 + t1) * t6 * t5;
      double t99 =
        MapleGenVar1 + 0.1850321001886696E4 * (-t3 - t12 -
                 0.1371740148509607E1) * (-t3 -
                        t12 -
                        0.1091700181433142E1)
        * (-t3 - t12 - 0.7092992179024789) * (-t3 - t12 -
                0.2907007820975211) * (-t3 - t12 +
                     0.917001814331423E-1)
        * (-t3 - t12 + 0.3717401485096066) * t6 * t5 -
        0.9225687774179959E2 * t75 * t72 * t71 * t35 +
        0.1205290876568007E1 * t84 * t83 * t81 * t80 * t35 +
        0.1044928509097744E3 * t26 * t24 * t5 +
        0.7083627997564722E3 * t40 * t38 * t6 * t5 +
        0.7544593661307843E2 * t75 * t10 * t5;
      double t100 = t26 * t23;
      double t101 = t40 * t38;
      double t137 = -t3 - t12 + 0.3302238962785669;
      double t140 = -t3 - t12 - 0.3115120652928579E-1;
      double t142 = -t3 - t12 - 0.9688487934707142;
      double t143 = -t3 - t12 - 0.1330223896278567E1;
      double t154 =
        0.549804259066764E3 * t101 * t100 * t35 -
        0.1175929080145098E3 * t75 * t6 * t5 - 0.5212211172635059E2 * t11 -
        0.2328101041304904E3 * t75 * t100 * t35 -
        0.9898064554760503E1 * t72 * t71 * t35 +
        0.9857756960509461E1 * t75 * t84 * t83 * t81 * t80 * t6 * t5 +
        0.1534633544851E4 * t30 * t29 * t27 * t35 -
        0.8422150708862237E3 * t31 * t27 * t9 * t35 -
        0.6112075447357898E2 * t101 * t72 * t48 * t70 * t6 * t5 +
        0.7844325172297547E2 * t35 +
        0.167522671187423E4 * t143 * t142 * t29 * t140 * t137 * t6 * t5 -
        0.1104965769548309E4 * t143 * t142 * t29 * t140 * t137 * t11;
      double t155 = t99 + t154;
      return t155;
    }

    double
    ortho2_f50x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t9 = -t3 - t8 + 0.3302238962785669;
      double t10 = -t3 - t8 - 0.3115120652928579E-1;
      double t11 = t10 * t9;
      double t12 = -t3 - t8 - 0.5;
      double t13 = -t3 - t8 - 0.9688487934707142;
      double t15 = t13 * t12 * t11;
      double t18 = -t3 - t8 - 0.1330223896278567E1;
      double t20 = t18 * t12 * t11;
      double t23 = 0.1E1 * x;
      double t24 = 0.1330223896278567E1 + t23 + t1;
      double t25 = 0.9688487934707142 + t23 + t1;
      double t26 = t25 * t24;
      double t27 = 0.5 + t23 + t1;
      double t28 = 0.3115120652928579E-1 + t23 + t1;
      double t33 = -0.3302238962785669 + t23 + t1;
      double t38 = t33 * t28;
      double t47 = t38 * t27 * t25;
      double t53 = t6 * t2;
      double t57 = t18 * t13;
      double t58 = t57 * t11;
      double t62 = t57 * t12 * t9;
      double t66 = t57 * t12 * t10;
      double t75 = 0.1265055323929465E1 + t23 + t1;
      double t76 = 0.7852315164806451 + t23 + t1;
      double t77 = t76 * t75;
      double t78 = 0.2147684835193549 + t23 + t1;
      double t79 = -0.2650553239294647 + t23 + t1;
      double t84 =
        -0.8376133559371152E3 * t15 * t7 - 0.8376133559371152E3 * t20 * t7 -
        0.4148207212977512E1 * t28 * t27 * t26 * t7 -
        0.4148207212977512E1 * t33 * t27 * t26 * t7 -
        0.4148207212977512E1 * t38 * t26 * t7 -
        0.4148207212977512E1 * t38 * t27 * t24 * t7 -
        0.4148207212977512E1 * t47 * t7 - 0.2074103606488756E1 * t47 * t24 * t5 +
        0.2074103606488756E1 * t47 * t24 * t53 - 0.8376133559371152E3 * t58 * t7 -
        0.8376133559371152E3 * t62 * t7 - 0.8376133559371152E3 * t66 * t7 +
        0.8376133559371152E3 * t66 * t9 * t5 -
        0.8376133559371152E3 * t66 * t9 * t53 +
        0.1558648230779421E2 * t79 * t78 * t77 * t7;
      double t85 = 0.1154653670707977E1 + t23 + t1;
      double t86 = t27 * t85;
      double t87 = -0.1546536707079771 + t23 + t1;
      double t88 = -t3 - t8 - 0.5278640450004206E-1;
      double t93 = -t3 - t8 - 0.9472135954999579;
      double t98 = t93 * t88;
      double t102 = t87 * t85;
      double t106 = t87 * t27;
      double t107 = t98 * t106;
      double t110 = t85 * t5;
      double t113 = t85 * t53;
      double t116 = -t3 - t8 + 0.2650553239294647;
      double t117 = -t3 - t8 - 0.2147684835193549;
      double t118 = t117 * t116;
      double t119 = -t3 - t8 - 0.7852315164806451;
      double t120 = -t3 - t8 - 0.1265055323929465E1;
      double t121 = t120 * t119;
      double t122 = t121 * t118;
      double t127 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t128 = t127 * t5;
      double t131 = t127 * t53;
      double t134 = t116 * t127;
      double t146 = t119 * t117;
      double t150 = 0.9472135954999579 + t23 + t1;
      double t151 = t150 * t5;
      double t152 = 0.5278640450004206E-1 + t23 + t1;
      double t153 = -t3 - t8 + 0.1546536707079771;
      double t155 = -t3 - t8 - 0.1154653670707977E1;
      double t156 = t155 * t12;
      double t157 = t156 * t153 * t152;
      double t160 = t150 * t53;
      double t163 =
        0.3056037723678949E2 * t88 * t87 * t86 * t7 +
        0.3056037723678949E2 * t93 * t87 * t86 * t7 -
        0.6112075447357898E2 * t98 * t86 * t7 -
        0.6112075447357898E2 * t98 * t102 * t7 -
        0.6112075447357898E2 * t107 * t7 - 0.3056037723678949E2 * t107 * t110 +
        0.3056037723678949E2 * t107 * t113 + 0.281096817511239E4 * t122 * t7 -
        0.4444530931801093E3 * t122 * t128 + 0.4444530931801093E3 * t122 * t131 +
        0.4444530931801093E3 * t121 * t134 * t7 +
        0.4444530931801093E3 * t121 * t117 * t127 * t7 +
        0.4444530931801093E3 * t120 * t117 * t134 * t7 +
        0.4444530931801093E3 * t146 * t134 * t7 +
        0.1216831944113737E4 * t157 * t151 - 0.1216831944113737E4 * t157 * t160;
      double t167 = t152 * t150;
      double t168 = t12 * t153;
      double t185 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t190 = t185 * t79;
      double t194 = t78 * t75;
      double t198 = t78 * t76;
      double t199 = t190 * t198;
      double t202 = t75 * t5;
      double t205 = t75 * t53;
      double t210 = t185 * t127;
      double t213 = t127 * t6;
      double t214 = t213 * t5;
      double t219 = t88 * t6;
      double t222 =
        0.2433663888227475E4 * t157 * t7 -
        0.1216831944113737E4 * t168 * t167 * t7 -
        0.1216831944113737E4 * t155 * t153 * t167 * t7 -
        0.1216831944113737E4 * t156 * t167 * t7 +
        0.2433663888227475E4 * t156 * t153 * t150 * t7 +
        0.9857756960509461E1 * t185 * t78 * t77 * t7 +
        0.9857756960509461E1 * t190 * t77 * t7 +
        0.9857756960509461E1 * t190 * t194 * t7 +
        0.9857756960509461E1 * t199 * t7 + 0.492887848025473E1 * t199 * t202 -
        0.492887848025473E1 * t199 * t205 - 0.3541813998782361E3 * t98 * t53 -
        0.3772296830653921E2 * t210 * t53 + 0.1192904999510088E3 * t214 -
        0.2385809999020175E3 * t185 * t6 * t5 - 0.3541813998782361E3 * t219 * t5;
      double t232 = t152 * t6;
      double t235 = t150 * t6;
      double t247 = -t3 - t8 + 0.3717401485096066;
      double t248 = -t3 - t8 + 0.917001814331423E-1;
      double t249 = t248 * t247;
      double t251 = -t3 - t8 - 0.2907007820975211;
      double t252 = -t3 - t8 - 0.7092992179024789;
      double t254 = -t3 - t8 - 0.1091700181433142E1;
      double t255 = -t3 - t8 - 0.1371740148509607E1;
      double t256 = t255 * t254;
      double t257 = t256 * t252 * t251;
      double t271 =
        -0.3541813998782361E3 * t93 * t6 * t5 + 0.3541813998782361E3 * t98 * t5 +
        0.3772296830653921E2 * t210 * t5 + 0.5224642545488721E2 * t167 * t5 +
        0.1044928509097744E3 * t232 * t5 + 0.1044928509097744E3 * t235 * t5 -
        0.5224642545488721E2 * t167 * t53 + 0.5879645400725488E2 * t185 * t53 -
        0.211061234912209E2 * t7 - 0.2606105586317529E2 * t128 +
        0.2606105586317529E2 * t131 - 0.5879645400725488E2 * t185 * t5 +
        0.9251605009433481E3 * t257 * t249 * t5 -
        0.9251605009433481E3 * t257 * t249 * t53 +
        0.349420856829338E4 * t66 * t9 * t6 * t5 -
        0.9251605009433481E3 * t257 * t248 * t6 * t5;
      double t275 = t247 * t6 * t5;
      double t276 = t251 * t248;
      double t294 = t9 * t127;
      double t319 =
        -0.9251605009433481E3 * t256 * t276 * t275 -
        0.9251605009433481E3 * t256 * t252 * t248 * t275 -
        0.9251605009433481E3 * t257 * t275 -
        0.9251605009433481E3 * t254 * t252 * t276 * t275 -
        0.9251605009433481E3 * t255 * t252 * t276 * t275 +
        0.5524828847741546E3 * t66 * t294 * t53 -
        0.5524828847741546E3 * t66 * t294 * t5 +
        0.5524828847741546E3 * t66 * t214 + 0.5524828847741546E3 * t58 * t214 +
        0.5524828847741546E3 * t62 * t214 + 0.5524828847741546E3 * t15 * t214 +
        0.5524828847741546E3 * t20 * t214 -
        0.368105095676662E3 * t152 * t235 * t5 +
        0.1475760185703521E4 * t93 * t219 * t5 - 0.3922162586148773E2 * t53 +
        0.3922162586148773E2 * t5;
      double t320 = t87 * t86;
      double t326 = t185 * t167;
      double t331 = t155 * t168;
      double t341 = t120 * t146;
      double t344 = t185 * t106;
      double t367 =
        -0.1458709317398891E3 * t320 * t7 -
        0.9225687774179959E2 * t185 * t86 * t7 +
        0.1164050520652452E3 * t326 * t53 - 0.1164050520652452E3 * t326 * t5 +
        0.2663317903720633E4 * t331 * t7 -
        0.2328101041304904E3 * t185 * t232 * t5 -
        0.2328101041304904E3 * t185 * t235 * t5 -
        0.8794991311194317E3 * t341 * t116 * t53 -
        0.9225687774179959E2 * t344 * t7 + 0.4949032277380252E1 * t320 * t53 +
        0.8794991311194317E3 * t341 * t116 * t5 -
        0.9225687774179959E2 * t185 * t102 * t7 -
        0.8794991311194317E3 * t341 * t7 -
        0.8794991311194317E3 * t120 * t118 * t7 -
        0.8794991311194317E3 * t120 * t119 * t116 * t7 -
        0.8794991311194317E3 * t119 * t118 * t7;
      double t373 = t153 * t6;
      double t390 = t85 * t6;
      double t398 = t93 * t88 * t127;
      double t415 =
        -0.7673167724255E3 * t155 * t12 * t6 * t5 -
        0.7673167724255E3 * t155 * t373 * t5 - 0.7673167724255E3 * t331 * t53 +
        0.7673167724255E3 * t331 * t5 - 0.7673167724255E3 * t12 * t373 * t5 -
        0.4949032277380252E1 * t320 * t5 -
        0.9898064554760503E1 * t87 * t27 * t6 * t5 -
        0.9898064554760503E1 * t87 * t390 * t5 -
        0.9898064554760503E1 * t27 * t390 * t5 +
        0.2333381733508091E3 * t398 * t53 - 0.2333381733508091E3 * t398 * t5 +
        0.2333381733508091E3 * t93 * t213 * t5 +
        0.461284388708998E2 * t344 * t113 +
        0.2333381733508091E3 * t88 * t213 * t5 -
        0.461284388708998E2 * t344 * t110 + 0.4211075354431118E3 * t331 * t131;
      double t416 = t153 * t127;
      double t430 = t93 * t88 * t152;
      double t447 = t79 * t198;
      double t463 =
        0.4211075354431118E3 * t155 * t416 * t7 +
        0.4211075354431118E3 * t155 * t12 * t127 * t7 -
        0.4211075354431118E3 * t331 * t128 +
        0.4211075354431118E3 * t12 * t416 * t7 -
        0.274902129533382E3 * t430 * t160 + 0.274902129533382E3 * t430 * t151 +
        0.549804259066764E3 * t430 * t7 +
        0.549804259066764E3 * t93 * t88 * t150 * t7 -
        0.274902129533382E3 * t93 * t167 * t7 -
        0.274902129533382E3 * t88 * t167 * t7 - 0.6026454382840033 * t447 * t205 +
        0.6026454382840033 * t447 * t202 + 0.1205290876568007E1 * t447 * t7 +
        0.1205290876568007E1 * t79 * t194 * t7 +
        0.1205290876568007E1 * t78 * t77 * t7 +
        0.1205290876568007E1 * t79 * t77 * t7;
      double t466 = t84 + t163 + t222 + t271 + t319 + t367 + t415 + t463;
      return t466;
    }

    double
    ortho2_f50y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = t10 * t5;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 + 0.3302238962785669;
      double t14 = -t3 - t12 - 0.5;
      double t16 = -t3 - t12 - 0.9688487934707142;
      double t17 = -t3 - t12 - 0.1330223896278567E1;
      double t18 = t17 * t16;
      double t19 = t18 * t14 * t13;
      double t22 = t6 * t2;
      double t23 = t9 * t22;
      double t25 = t6 * t5;
      double t27 = t6 * t4;
      double t30 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t35 = t9 * t27;
      double t37 = -t3 - t12 - 0.5278640450004206E-1;
      double t39 = -t3 - t12 - 0.9472135954999579;
      double t40 = t39 * t37 * t9;
      double t43 = 0.1E1 * x;
      double t44 = 0.1154653670707977E1 + t43 + t1;
      double t45 = t44 * t6;
      double t46 = 0.5 + t43 + t1;
      double t50 = 0.7852315164806451 + t43 + t1;
      double t51 = 0.2147684835193549 + t43 + t1;
      double t52 = t51 * t50;
      double t53 = -0.2650553239294647 + t43 + t1;
      double t54 = t53 * t52;
      double t57 = -t3 - t12 - 0.3115120652928579E-1;
      double t58 = t57 * t13;
      double t59 = t18 * t58;
      double t62 = 0.1330223896278567E1 + t43 + t1;
      double t63 = 0.9688487934707142 + t43 + t1;
      double t64 = t63 * t62;
      double t65 = 0.3115120652928579E-1 + t43 + t1;
      double t66 = -0.3302238962785669 + t43 + t1;
      double t67 = t66 * t65;
      double t73 = t67 * t46 * t63;
      double t76 = t44 * t27;
      double t77 = -0.1546536707079771 + t43 + t1;
      double t78 = t77 * t46;
      double t79 = t39 * t37;
      double t80 = t79 * t78;
      double t83 = t46 * t44;
      double t90 = t18 * t14 * t57;
      double t93 = 0.1104965769548309E4 * t19 * t11 + 0.2606105586317529E2 * t23
        - 0.289449131252E3 * t25 - 0.5879645400725488E2 * t30 * t27 +
        0.5879645400725488E2 * t30 * t22 - 0.2606105586317529E2 * t35 -
        0.2333381733508091E3 * t40 * t27 - 0.4949032277380252E1 * t46 * t45 * t5 +
        0.6026454382840033 * t54 * t25 + 0.1104965769548309E4 * t59 * t11 -
        0.2074103606488756E1 * t67 * t64 * t25 -
        0.2074103606488756E1 * t73 * t62 * t27 -
        0.3056037723678949E2 * t80 * t76 +
        0.6112075447357898E2 * t39 * t77 * t83 * t25 -
        0.8376133559371152E3 * t90 * t13 * t22;
      double t94 = -t3 - t12 + 0.1546536707079771;
      double t95 = t14 * t94;
      double t96 = -t3 - t12 - 0.1154653670707977E1;
      double t97 = t96 * t95;
      double t100 = t94 * t9;
      double t104 = -t3 - t12 + 0.2650553239294647;
      double t106 = -t3 - t12 - 0.2147684835193549;
      double t107 = -t3 - t12 - 0.7852315164806451;
      double t108 = t107 * t106;
      double t109 = -t3 - t12 - 0.1265055323929465E1;
      double t110 = t109 * t108;
      double t113 = 0.1265055323929465E1 + t43 + t1;
      double t114 = t50 * t113;
      double t115 = t30 * t53;
      double t119 = 0.9472135954999579 + t43 + t1;
      double t120 = 0.5278640450004206E-1 + t43 + t1;
      double t121 = t120 * t119;
      double t122 = t30 * t121;
      double t125 = t13 * t9;
      double t133 = t106 * t104;
      double t137 = t77 * t83;
      double t140 = t109 * t107;
      double t141 = t140 * t133;
      double t144 = t119 * t22;
      double t146 = t96 * t14;
      double t147 = t146 * t94 * t120;
      double t153 = -t3 - t12 + 0.3717401485096066;
      double t155 = t153 * t6 * t5;
      double t156 = -t3 - t12 - 0.2907007820975211;
      double t157 = -t3 - t12 - 0.7092992179024789;
      double t159 = -t3 - t12 - 0.1091700181433142E1;
      double t160 = -t3 - t12 - 0.1371740148509607E1;
      double t161 = t160 * t159;
      double t162 = t161 * t157 * t156;
      double t173 = t94 * t6;
      double t177 =
        -0.4211075354431118E3 * t97 * t35 +
        0.8422150708862237E3 * t96 * t100 * t25 -
        0.8794991311194317E3 * t110 * t104 * t22 +
        0.492887848025473E1 * t115 * t114 * t25 -
        0.1164050520652452E3 * t122 * t27 -
        0.5524828847741546E3 * t90 * t125 * t27 +
        0.492887848025473E1 * t30 * t51 * t114 * t25 -
        0.1758998262238863E4 * t107 * t133 * t25 +
        0.4949032277380252E1 * t137 * t22 + 0.1405484087556195E4 * t141 * t25 -
        0.1216831944113737E4 * t147 * t144 +
        0.8794991311194317E3 * t110 * t104 * t27 -
        0.1850321001886696E4 * t162 * t155 -
        0.2074103606488756E1 * t66 * t46 * t64 * t25 -
        0.1534633544851E4 * t96 * t14 * t6 * t5 -
        0.1534633544851E4 * t96 * t173 * t5;
      double t181 = t30 * t78;
      double t190 = t51 * t113;
      double t196 = t37 * t6;
      double t202 = t30 * t9;
      double t207 = t77 * t44;
      double t219 = t120 * t6;
      double t222 =
        0.5224642545488721E2 * t121 * t27 - 0.461284388708998E2 * t181 * t25 +
        0.6112075447357898E2 * t37 * t77 * t83 * t25 -
        0.4949032277380252E1 * t137 * t27 +
        0.6026454382840033 * t53 * t190 * t25 + 0.3541813998782361E3 * t79 * t27 -
        0.7083627997564722E3 * t196 * t5 - 0.7083627997564722E3 * t39 * t6 * t5 +
        0.3772296830653921E2 * t202 * t27 - 0.167522671187423E4 * t90 * t25 -
        0.461284388708998E2 * t30 * t207 * t25 -
        0.3772296830653921E2 * t202 * t22 + 0.2385809999020175E3 * t11 +
        0.4444530931801093E3 * t141 * t23 - 0.1192904999510088E3 * t30 * t6 * t5 +
        0.5224642545488721E2 * t219 * t5;
      double t223 = t119 * t6;
      double t240 = t39 * t37 * t120;
      double t253 = t17 * t14 * t58;
      double t266 =
        0.5224642545488721E2 * t223 * t5 - 0.5224642545488721E2 * t121 * t22 -
        0.2074103606488756E1 * t65 * t46 * t64 * t25 +
        0.1216831944113737E4 * t147 * t25 - 0.7673167724255E3 * t97 * t22 -
        0.3922162586148773E2 * t22 + 0.1164050520652452E3 * t122 * t22 -
        0.274902129533382E3 * t240 * t144 + 0.2333381733508091E3 * t40 * t22 +
        0.274902129533382E3 * t240 * t25 + 0.4666763467016183E3 * t39 * t10 * t5 +
        0.1331658951860317E4 * t97 * t25 - 0.167522671187423E4 * t253 * t25 +
        0.274902129533382E3 * t39 * t37 * t119 * t25 +
        0.8376133559371152E3 * t90 * t13 * t27 +
        0.4666763467016183E3 * t37 * t10 * t5;
      double t291 = t44 * t22;
      double t297 = -t3 - t12 + 0.917001814331423E-1;
      double t320 =
        -0.167522671187423E4 * t19 * t25 -
        0.549804259066764E3 * t39 * t121 * t25 -
        0.1534633544851E4 * t14 * t173 * t5 +
        0.3117296461558841E2 * t53 * t51 * t114 * t25 +
        0.1216831944113737E4 * t146 * t94 * t119 * t25 +
        0.6026454382840033 * t53 * t114 * t25 -
        0.2433663888227475E4 * t146 * t121 * t25 +
        0.3056037723678949E2 * t80 * t291 +
        0.6026454382840033 * t51 * t114 * t25 -
        0.1850321001886696E4 * t161 * t157 * t297 * t155 +
        0.461284388708998E2 * t181 * t291 +
        0.8889061863602187E3 * t140 * t106 * t9 * t25 -
        0.2433663888227475E4 * t96 * t94 * t121 * t25 -
        0.7362101913533239E3 * t120 * t223 * t5 -
        0.2917418634797783E3 * t137 * t25 -
        0.2433663888227475E4 * t95 * t121 * t25;
      double t321 = t297 * t153;
      double t338 = t156 * t297;
      double t342 = t119 * t27;
      double t345 = t113 * t27;
      double t348 = t104 * t9;
      double t370 =
        -0.9251605009433481E3 * t162 * t321 * t22 -
        0.549804259066764E3 * t37 * t121 * t25 +
        0.2074103606488756E1 * t73 * t62 * t22 +
        0.7378800928517604E3 * t39 * t196 * t5 -
        0.1850321001886696E4 * t162 * t297 * t6 * t5 -
        0.1850321001886696E4 * t161 * t338 * t155 +
        0.1216831944113737E4 * t147 * t342 + 0.6026454382840033 * t54 * t345 +
        0.8889061863602187E3 * t140 * t348 * t25 -
        0.1758998262238863E4 * t110 * t25 + 0.7673167724255E3 * t97 * t27 -
        0.1850321001886696E4 * t157 * t160 * t338 * t155 -
        0.3056037723678949E2 * t80 * t25 - 0.3541813998782361E3 * t79 * t22 +
        0.174710428414669E4 * t90 * t13 * t6 * t5 +
        0.274902129533382E3 * t240 * t342;
      double t374 = t113 * t22;
      double t375 = t115 * t52;
      double t420 =
        -0.461284388708998E2 * t181 * t76 - 0.492887848025473E1 * t375 * t374 +
        0.8889061863602187E3 * t109 * t106 * t348 * t25 -
        0.1758998262238863E4 * t109 * t107 * t104 * t25 +
        0.4211075354431118E3 * t97 * t23 -
        0.1164050520652452E3 * t30 * t219 * t5 -
        0.1850321001886696E4 * t159 * t157 * t338 * t155 +
        0.9251605009433481E3 * t162 * t321 * t27 +
        0.8422150708862237E3 * t14 * t100 * t25 -
        0.4949032277380252E1 * t77 * t46 * t6 * t5 +
        0.1104965769548309E4 * t253 * t11 +
        0.5524828847741546E3 * t90 * t125 * t22 -
        0.3056037723678949E2 * t79 * t207 * t25 +
        0.8889061863602187E3 * t108 * t348 * t25 +
        0.492887848025473E1 * t375 * t25 - 0.2074103606488756E1 * t73 * t25;
      double t422 = t16 * t14 * t58;
      double t464 =
        0.1104965769548309E4 * t422 * t11 -
        0.1164050520652452E3 * t30 * t223 * t5 - 0.167522671187423E4 * t59 * t25 +
        0.3922162586148773E2 * t27 -
        0.2074103606488756E1 * t67 * t46 * t62 * t25 +
        0.1104965769548309E4 * t90 * t11 - 0.461284388708998E2 * t30 * t83 * t25 -
        0.4444530931801093E3 * t141 * t35 +
        0.8422150708862237E3 * t96 * t14 * t9 * t25 -
        0.4949032277380252E1 * t77 * t45 * t5 - 0.6026454382840033 * t54 * t374 +
        0.492887848025473E1 * t115 * t190 * t25 -
        0.1758998262238863E4 * t109 * t133 * t25 -
        0.167522671187423E4 * t422 * t25 -
        0.3056037723678949E2 * t79 * t83 * t25 +
        0.492887848025473E1 * t375 * t345;
      double t467 = t93 + t177 + t222 + t266 + t320 + t370 + t420 + t464;
      return t467;
    }

    // * f51 *********************************************************************

    double
    ortho2_f51 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = t10 * t5;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 + 0.2650553239294647;
      double t14 = -t3 - t12 - 0.2147684835193549;
      double t16 = -t3 - t12 - 0.7852315164806451;
      double t17 = -t3 - t12 - 0.1265055323929465E1;
      double t19 = t17 * t16 * t14 * t13;
      double t22 = 0.1E1 * x;
      double t23 = 0.9472135954999579 + t22 + t1;
      double t24 = t23 * t6;
      double t25 = t24 * t5;
      double t26 = 0.5278640450004206E-1 + t22 + t1;
      double t27 = -t3 - t12 + 0.1546536707079771;
      double t29 = -t3 - t12 - 0.5;
      double t30 = -t3 - t12 - 0.1154653670707977E1;
      double t31 = t30 * t29;
      double t35 = t6 * t5;
      double t38 = -t3 - t12 - 0.5278640450004206E-1;
      double t40 = -t3 - t12 - 0.9472135954999579;
      double t48 = 0.5 + t22 + t1;
      double t70 = 0.1154653670707977E1 + t22 + t1;
      double t71 = t48 * t70;
      double t72 = -0.1546536707079771 + t22 + t1;
      double t75 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t80 = 0.1265055323929465E1 + t22 + t1;
      double t81 = 0.7852315164806451 + t22 + t1;
      double t83 = 0.2147684835193549 + t22 + t1;
      double t84 = -0.2650553239294647 + t22 + t1;
      double t99 = t26 * t23;
      double t100 = t40 * t38;
      double MapleGenVar1 =
        -0.1618543974138404E4 * t19 * t11 +
        0.4941463028333687E4 * t31 * t27 * t26 * t25 +
        0.3259684839629149E4 * t19 * t35 -
        0.1202181340260288E4 * t40 * t38 * t9 * t35 -
        0.8624392878061751E1 * (-0.3302238962785669 + t22 +
              t1) * (0.3115120652928579E-1 + t22 +
               t1) * t48 * (0.9688487934707142 + t22 +
                t1) * (0.1330223896278567E1 +
                       t22 + t1) * t6 * t5 +
        0.1870935748025166E4 * (-t3 - t12 - 0.1371740148509607E1) * (-t3 - t12 -
                     0.1091700181433142E1)
        * (-t3 - t12 - 0.7092992179024789) * (-t3 - t12 -
                0.2907007820975211) * (-t3 - t12 +
                     0.917001814331423E-1)
        * (-t3 - t12 + 0.3717401485096066) * t6 * t5;
      double t104 =
        MapleGenVar1 - 0.1070075643657266E3 * t75 * t72 * t71 * t35 +
        0.1154193062952326E3 * t84 * t83 * t81 * t80 * t35 +
        0.5217306959738636E3 * t26 * t24 * t5 +
        0.1451178479755592E4 * t40 * t38 * t6 * t5 +
        0.1437749442579747E3 * t75 * t10 * t5 +
        0.3770257664377148E4 * t100 * t99 * t35;
      double t143 = -t3 - t12 + 0.3302238962785669;
      double t146 = -t3 - t12 - 0.3115120652928579E-1;
      double t148 = -t3 - t12 - 0.9688487934707142;
      double t149 = -t3 - t12 - 0.1330223896278567E1;
      double t160 =
        0.6648285590338501E4 * t17 * t16 * t14 * t13 * t26 * t25 -
        0.1981091017543117E3 * t75 * t6 * t5 - 0.1512821197015158E3 * t11 -
        0.5753867520189039E3 * t75 * t99 * t35 +
        0.4287430975168431E3 * t72 * t71 * t35 +
        0.972140815863854E1 * t75 * t84 * t83 * t81 * t80 * t6 * t5 +
        0.2187016565401361E4 * t30 * t29 * t27 * t35 -
        0.1740910333200721E4 * t31 * t27 * t9 * t35 +
        0.2822891849735354E4 * t100 * t72 * t48 * t70 * t6 * t5 +
        0.1679004141043073E3 * t35 +
        0.1754475502216867E4 * t149 * t148 * t29 * t146 * t143 * t6 * t5 -
        0.2234318048354632E4 * t149 * t148 * t29 * t146 * t143 * t11;
      double t161 = t104 + t160;
      return t161;
    }

    double
    ortho2_f51x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.1E1 * x;
      double t7 = 0.9472135954999579 + t6 + t1;
      double t8 = t7 * t5;
      double t9 = 0.5278640450004206E-1 + t6 + t1;
      double t10 = 0.1E1 * y;
      double t11 = -t3 - t10 - 0.5278640450004206E-1;
      double t13 = -t3 - t10 - 0.9472135954999579;
      double t14 = t13 * t11 * t9;
      double t17 = 0.5 + t3;
      double t18 = t17 * t5;
      double t21 = t17 * t2;
      double t24 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t25 = -t3 - t10 + 0.3302238962785669;
      double t26 = t25 * t24;
      double t28 = -t3 - t10 - 0.3115120652928579E-1;
      double t29 = -t3 - t10 - 0.5;
      double t31 = -t3 - t10 - 0.9688487934707142;
      double t32 = -t3 - t10 - 0.1330223896278567E1;
      double t33 = t32 * t31;
      double t34 = t33 * t29 * t28;
      double t41 = t9 * t7;
      double t45 = -t3 - t10 + 0.1546536707079771;
      double t47 = -t3 - t10 - 0.1154653670707977E1;
      double t48 = t47 * t29;
      double t52 = -t3 - t10 + 0.3717401485096066;
      double t53 = -t3 - t10 + 0.917001814331423E-1;
      double t54 = t53 * t52;
      double t56 = -t3 - t10 - 0.2907007820975211;
      double t57 = -t3 - t10 - 0.7092992179024789;
      double t59 = -t3 - t10 - 0.1091700181433142E1;
      double t60 = -t3 - t10 - 0.1371740148509607E1;
      double t61 = t60 * t59;
      double t62 = t61 * t57 * t56;
      double t65 = t7 * t21;
      double t67 = t48 * t45 * t9;
      double t73 = t29 * t45;
      double t87 = t24 * t17;
      double t88 = t87 * t5;
      double t89 = t28 * t25;
      double t91 = t32 * t29 * t89;
      double t95 = t31 * t29 * t89;
      double t102 =
        0.1885128832188574E4 * t14 * t8 + 0.3770257664377148E4 * t14 * t18 +
        0.1117159024177316E4 * t34 * t26 * t21 +
        0.3770257664377148E4 * t13 * t11 * t7 * t18 -
        0.1885128832188574E4 * t13 * t41 * t18 +
        0.4941463028333687E4 * t48 * t45 * t7 * t18 -
        0.9354678740125828E3 * t62 * t54 * t21 -
        0.2470731514166843E4 * t67 * t65 -
        0.2470731514166843E4 * t48 * t41 * t18 -
        0.2470731514166843E4 * t73 * t41 * t18 +
        0.9354678740125828E3 * t62 * t54 * t5 -
        0.2470731514166843E4 * t47 * t45 * t41 * t18 -
        0.1885128832188574E4 * t11 * t41 * t18 +
        0.1117159024177316E4 * t91 * t88 + 0.1117159024177316E4 * t95 * t88 -
        0.9354678740125828E3 * t62 * t53 * t17 * t5;
      double t103 = 0.1265055323929465E1 + t6 + t1;
      double t104 = t103 * t21;
      double t105 = 0.7852315164806451 + t6 + t1;
      double t106 = 0.2147684835193549 + t6 + t1;
      double t107 = t106 * t105;
      double t108 = -0.2650553239294647 + t6 + t1;
      double t109 = t108 * t107;
      double t112 = t103 * t5;
      double t117 = t106 * t103;
      double t121 = t105 * t103;
      double t134 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t135 = t134 * t108;
      double t136 = t135 * t107;
      double t149 = 0.1154653670707977E1 + t6 + t1;
      double t150 = 0.5 + t6 + t1;
      double t151 = t150 * t149;
      double t152 = -0.1546536707079771 + t6 + t1;
      double t153 = t152 * t151;
      double t160 = t41 * t21;
      double t161 = -t3 - t10 + 0.2650553239294647;
      double t162 = -t3 - t10 - 0.2147684835193549;
      double t163 = t162 * t161;
      double t164 = -t3 - t10 - 0.7852315164806451;
      double t165 = -t3 - t10 - 0.1265055323929465E1;
      double t166 = t165 * t164;
      double t167 = t166 * t163;
      double t170 = t9 * t17;
      double t171 = t170 * t5;
      double t174 = t41 * t5;
      double t177 =
        -0.577096531476163E2 * t109 * t104 + 0.577096531476163E2 * t109 * t112 +
        0.1154193062952326E3 * t109 * t18 +
        0.1154193062952326E3 * t108 * t117 * t18 +
        0.1154193062952326E3 * t108 * t121 * t18 +
        0.1154193062952326E3 * t106 * t121 * t18 +
        0.7065534050022867E4 * t34 * t25 * t17 * t5 -
        0.486070407931927E1 * t136 * t104 + 0.486070407931927E1 * t136 * t112 +
        0.972140815863854E1 * t136 * t18 +
        0.972140815863854E1 * t135 * t121 * t18 +
        0.972140815863854E1 * t135 * t117 * t18 +
        0.2143715487584216E3 * t153 * t5 +
        0.972140815863854E1 * t134 * t106 * t121 * t18 -
        0.3324142795169251E4 * t167 * t160 + 0.6648285590338501E4 * t167 * t171 +
        0.3324142795169251E4 * t167 * t174;
      double t179 = t7 * t17;
      double t180 = t179 * t5;
      double t191 = t149 * t17;
      double t195 = t161 * t9;
      double t199 = t165 * t162;
      double t206 = t164 * t162;
      double t210 = t24 * t21;
      double t218 = t24 * t5;
      double t226 = t161 * t24;
      double t233 =
        0.6648285590338501E4 * t167 * t180 +
        0.4287430975168431E3 * t152 * t150 * t17 * t5 -
        0.3324142795169251E4 * t166 * t162 * t9 * t180 +
        0.4287430975168431E3 * t152 * t191 * t5 -
        0.3324142795169251E4 * t166 * t195 * t180 -
        0.3324142795169251E4 * t199 * t195 * t180 +
        0.4287430975168431E3 * t150 * t191 * t5 -
        0.3324142795169251E4 * t206 * t195 * t180 +
        0.8092719870692021E3 * t167 * t210 + 0.7564105985075791E2 * t210 -
        0.9905455087715583E2 * t134 * t5 + 0.9905455087715583E2 * t134 * t21 -
        0.7564105985075791E2 * t218 - 0.8092719870692021E3 * t167 * t218 +
        0.8092719870692021E3 * t166 * t162 * t24 * t18 +
        0.8092719870692021E3 * t166 * t226 * t18 +
        0.8092719870692021E3 * t199 * t226 * t18;
      double t241 = t13 * t11 * t24;
      double t256 = t13 * t11;
      double t259 = t11 * t17;
      double t269 = t134 * t24;
      double t274 = t149 * t21;
      double t275 = t152 * t150;
      double t276 = t134 * t275;
      double t279 =
        0.8092719870692021E3 * t206 * t226 * t18 -
        0.1117159024177316E4 * t34 * t26 * t5 +
        0.6010906701301438E3 * t241 * t21 + 0.2470731514166843E4 * t67 * t8 -
        0.6010906701301438E3 * t241 * t5 + 0.6010906701301438E3 * t13 * t87 * t5 +
        0.1117159024177316E4 * t34 * t88 + 0.6010906701301438E3 * t11 * t87 * t5 +
        0.7255892398777961E3 * t256 * t5 - 0.7255892398777961E3 * t259 * t5 -
        0.7255892398777961E3 * t13 * t17 * t5 -
        0.4546562943189475E3 * t134 * t17 * t5 + 0.2273281471594737E3 * t88 +
        0.7188747212898736E2 * t269 * t5 + 0.5217306959738636E3 * t171 +
        0.2608653479869318E3 * t174 + 0.5350378218286328E2 * t276 * t274;
      double t284 = t149 * t5;
      double t289 = t47 * t73;
      double t293 = t33 * t29 * t25;
      double t298 = t256 * t275;
      double t310 = t152 * t149;
      double t322 = t45 * t24;
      double t330 =
        0.4941463028333687E4 * t67 * t18 - 0.5350378218286328E2 * t276 * t284 +
        0.5118285451418023E4 * t167 * t18 + 0.8704551666003607E3 * t289 * t210 +
        0.1117159024177316E4 * t293 * t88 - 0.8704551666003607E3 * t289 * t218 -
        0.1411445924867677E4 * t298 * t274 + 0.1411445924867677E4 * t298 * t284 +
        0.2822891849735354E4 * t298 * t18 + 0.8395020705215363E2 * t5 +
        0.8704551666003607E3 * t47 * t29 * t24 * t18 +
        0.2822891849735354E4 * t256 * t310 * t18 +
        0.2822891849735354E4 * t256 * t151 * t18 - 0.8395020705215363E2 * t21 -
        0.1411445924867677E4 * t13 * t152 * t151 * t18 +
        0.8704551666003607E3 * t47 * t322 * t18 -
        0.1411445924867677E4 * t11 * t152 * t151 * t18;
      double t352 = t52 * t17 * t5;
      double t363 = t33 * t89;
      double t368 = t56 * t53;
      double t377 = 0.1330223896278567E1 + t6 + t1;
      double t379 = 0.9688487934707142 + t6 + t1;
      double t381 = 0.3115120652928579E-1 + t6 + t1;
      double t382 = -0.3302238962785669 + t6 + t1;
      double t383 = t382 * t381;
      double t384 = t383 * t150 * t379;
      double t387 =
        0.1537089592272064E2 * t108 * t106 * t121 * t18 -
        0.9097663359331114E3 * t9 * t179 * t5 -
        0.7255892398777961E3 * t256 * t21 +
        0.3801631195776389E4 * t13 * t259 * t5 -
        0.8772377511084337E3 * t34 * t25 * t21 +
        0.8772377511084337E3 * t34 * t25 * t5 - 0.1885128832188574E4 * t14 * t65 -
        0.9354678740125828E3 * t62 * t352 - 0.8772377511084337E3 * t34 * t18 -
        0.8772377511084337E3 * t293 * t18 -
        0.9354678740125828E3 * t61 * t57 * t53 * t352 +
        0.1117159024177316E4 * t363 * t88 - 0.1691938151313835E3 * t153 * t18 -
        0.9354678740125828E3 * t61 * t368 * t352 -
        0.8772377511084337E3 * t363 * t18 -
        0.1070075643657266E3 * t134 * t151 * t18 +
        0.4312196439030875E1 * t384 * t377 * t21;
      double t398 = t379 * t377;
      double t410 = t134 * t41;
      double t430 = t165 * t206;
      double t439 =
        -0.8624392878061751E1 * t384 * t18 -
        0.4312196439030875E1 * t384 * t377 * t5 -
        0.8624392878061751E1 * t383 * t150 * t377 * t18 -
        0.8624392878061751E1 * t383 * t398 * t18 -
        0.8624392878061751E1 * t381 * t150 * t398 * t18 -
        0.8624392878061751E1 * t382 * t150 * t398 * t18 +
        0.2876933760094519E3 * t410 * t21 - 0.7188747212898736E2 * t269 * t21 -
        0.2876933760094519E3 * t410 * t5 + 0.5505241855036931E4 * t289 * t18 -
        0.8772377511084337E3 * t91 * t18 -
        0.5753867520189039E3 * t134 * t170 * t5 -
        0.5753867520189039E3 * t134 * t179 * t5 -
        0.8772377511084337E3 * t95 * t18 -
        0.1629842419814574E4 * t430 * t161 * t21 -
        0.9354678740125828E3 * t60 * t57 * t368 * t352 -
        0.2143715487584216E3 * t153 * t21;
      double t470 = t45 * t17;
      double t485 =
        0.5217306959738636E3 * t180 - 0.1070075643657266E3 * t276 * t18 -
        0.1070075643657266E3 * t134 * t310 * t18 + 0.1651580741381801E3 * t18 +
        0.1629842419814574E4 * t430 * t161 * t5 -
        0.1629842419814574E4 * t430 * t18 -
        0.9354678740125828E3 * t59 * t57 * t368 * t352 -
        0.1629842419814574E4 * t165 * t164 * t161 * t18 -
        0.1629842419814574E4 * t165 * t163 * t18 -
        0.1629842419814574E4 * t164 * t163 * t18 -
        0.109350828270068E4 * t47 * t29 * t17 * t5 -
        0.109350828270068E4 * t47 * t470 * t5 - 0.109350828270068E4 * t289 * t21 +
        0.8704551666003607E3 * t29 * t322 * t18 +
        0.109350828270068E4 * t289 * t5 - 0.109350828270068E4 * t29 * t470 * t5 -
        0.2608653479869318E3 * t160;
      double t488 = t102 + t177 + t233 + t279 + t330 + t387 + t439 + t485;
      return t488;
    }

    double
    ortho2_f51y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * y;
      double t7 = -t3 - t6 + 0.1546536707079771;
      double t8 = -t3 - t6 - 0.5;
      double t9 = t8 * t7;
      double t10 = -t3 - t6 - 0.1154653670707977E1;
      double t11 = t10 * t9;
      double t14 = -t3 - t1;
      double t15 = t14 * t2;
      double t16 = t4 * t15;
      double t20 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t21 = t20 * t5;
      double t23 = t4 * t14;
      double t26 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t29 = 0.1E1 * x;
      double t30 = 0.1154653670707977E1 + t29 + t1;
      double t31 = t30 * t5;
      double t32 = 0.5 + t29 + t1;
      double t33 = -0.1546536707079771 + t29 + t1;
      double t34 = t33 * t32;
      double t35 = t26 * t34;
      double t38 = t32 * t30;
      double t39 = -t3 - t6 - 0.5278640450004206E-1;
      double t40 = -t3 - t6 - 0.9472135954999579;
      double t41 = t40 * t39;
      double t45 = t20 * t23;
      double t49 = 0.1330223896278567E1 + t29 + t1;
      double t51 = 0.9688487934707142 + t29 + t1;
      double t53 = 0.3115120652928579E-1 + t29 + t1;
      double t54 = -0.3302238962785669 + t29 + t1;
      double t55 = t54 * t53;
      double t56 = t55 * t32 * t51;
      double t59 = -t3 - t6 + 0.3717401485096066;
      double t61 = t59 * t4 * t15;
      double t62 = -t3 - t6 + 0.917001814331423E-1;
      double t63 = -t3 - t6 - 0.7092992179024789;
      double t65 = -t3 - t6 - 0.1091700181433142E1;
      double t66 = -t3 - t6 - 0.1371740148509607E1;
      double t67 = t66 * t65;
      double t71 = t20 * t4;
      double t72 = t71 * t15;
      double t73 = -t3 - t6 + 0.3302238962785669;
      double t74 = -t3 - t6 - 0.3115120652928579E-1;
      double t75 = t74 * t73;
      double t76 = -t3 - t6 - 0.9688487934707142;
      double t77 = -t3 - t6 - 0.1330223896278567E1;
      double t78 = t77 * t76;
      double t79 = t78 * t75;
      double t82 = -t3 - t6 + 0.2650553239294647;
      double t83 = -t3 - t6 - 0.2147684835193549;
      double t84 = t83 * t82;
      double t85 = -t3 - t6 - 0.7852315164806451;
      double t86 = -t3 - t6 - 0.1265055323929465E1;
      double t87 = t86 * t85;
      double t88 = t87 * t84;
      double t94 = t78 * t8 * t74;
      double t97 = 0.5278640450004206E-1 + t29 + t1;
      double t98 = t97 * t4;
      double t99 = t98 * t15;
      double t106 = -0.109350828270068E4 * t11 * t5 - 0.387277952996173E3 * t16 +
        0.7564105985075791E2 * t21 - 0.9905455087715583E2 * t26 * t23 +
        0.5350378218286328E2 * t35 * t31 +
        0.1411445924867677E4 * t41 * t38 * t16 - 0.7564105985075791E2 * t45 +
        0.9905455087715583E2 * t26 * t5 + 0.4312196439030875E1 * t56 * t49 * t5 -
        0.1870935748025166E4 * t67 * t63 * t62 * t61 +
        0.2234318048354632E4 * t79 * t72 - 0.8092719870692021E3 * t88 * t45 -
        0.8704551666003607E3 * t11 * t45 - 0.1754475502216867E4 * t94 * t16 +
        0.3324142795169251E4 * t88 * t99 -
        0.2822891849735354E4 * t40 * t33 * t38 * t16;
      double t107 = t33 * t30;
      double t111 = 0.9472135954999579 + t29 + t1;
      double t112 = t111 * t5;
      double t114 = t10 * t8;
      double t115 = t114 * t7 * t97;
      double t118 = 0.1265055323929465E1 + t29 + t1;
      double t119 = t118 * t5;
      double t120 = 0.7852315164806451 + t29 + t1;
      double t121 = 0.2147684835193549 + t29 + t1;
      double t122 = t121 * t120;
      double t123 = -0.2650553239294647 + t29 + t1;
      double t124 = t123 * t122;
      double t127 = t111 * t4;
      double t128 = t127 * t15;
      double t132 = t78 * t8 * t73;
      double t141 = -t3 - t6 - 0.2907007820975211;
      double t142 = t141 * t62;
      double t152 = t62 * t59;
      double t155 = t67 * t63 * t141;
      double t163 = t40 * t39 * t20;
      double t168 = t51 * t49;
      double t176 = t118 * t23;
      double t179 =
        -0.5350378218286328E2 * t26 * t107 * t16 -
        0.2470731514166843E4 * t115 * t112 - 0.577096531476163E2 * t124 * t119 +
        0.3324142795169251E4 * t88 * t128 - 0.1754475502216867E4 * t132 * t16 -
        0.2822891849735354E4 * t39 * t33 * t38 * t16 -
        0.4312196439030875E1 * t56 * t16 -
        0.1870935748025166E4 * t67 * t142 * t61 +
        0.2470731514166843E4 * t115 * t16 -
        0.6648285590338501E4 * t87 * t83 * t97 * t128 +
        0.9354678740125828E3 * t155 * t152 * t23 -
        0.4312196439030875E1 * t55 * t32 * t49 * t16 -
        0.6010906701301438E3 * t163 * t23 + 0.577096531476163E2 * t124 * t16 -
        0.4312196439030875E1 * t55 * t168 * t16 -
        0.4312196439030875E1 * t54 * t32 * t168 * t16 +
        0.577096531476163E2 * t124 * t176;
      double t183 = t97 * t111;
      double t184 = t183 * t23;
      double t186 = t82 * t97;
      double t194 = t40 * t39 * t97;
      double t203 = t26 * t20;
      double t213 = t85 * t83;
      double t214 = t86 * t213;
      double t217 = t39 * t4;
      double t224 = t26 * t183;
      double t229 = t30 * t23;
      double t232 = 0.2559142725709011E4 * t88 * t16 + 0.2608653479869318E3 * t184
        - 0.6648285590338501E4 * t87 * t186 * t128 -
        0.1451178479755592E4 * t40 * t4 * t15 -
        0.1885128832188574E4 * t194 * t112 - 0.7255892398777961E3 * t41 * t5 -
        0.4312196439030875E1 * t53 * t32 * t168 * t16 -
        0.7188747212898736E2 * t203 * t5 + 0.4546562943189475E3 * t72 -
        0.2273281471594737E3 * t26 * t4 * t15 + 0.7255892398777961E3 * t41 * t23 +
        0.1629842419814574E4 * t214 * t82 * t23 -
        0.1451178479755592E4 * t217 * t15 -
        0.1870935748025166E4 * t66 * t63 * t142 * t61 +
        0.2876933760094519E3 * t224 * t5 + 0.7188747212898736E2 * t203 * t23 -
        0.5350378218286328E2 * t35 * t229;
      double t235 = t183 * t5;
      double t240 = t120 * t118;
      double t246 = t7 * t4;
      double t256 = t77 * t8 * t75;
      double t262 = t121 * t118;
      double t276 = t33 * t38;
      double t282 = 0.2752620927518466E4 * t11 * t16 - 0.2608653479869318E3 * t235
        - 0.4941463028333687E4 * t9 * t183 * t16 +
        0.3074179184544127E2 * t123 * t121 * t240 * t16 +
        0.2608653479869318E3 * t99 - 0.2187016565401361E4 * t8 * t246 * t15 +
        0.2608653479869318E3 * t128 -
        0.1870935748025166E4 * t65 * t63 * t142 * t61 -
        0.1754475502216867E4 * t256 * t16 +
        0.577096531476163E2 * t123 * t240 * t16 +
        0.577096531476163E2 * t123 * t262 * t16 -
        0.8772377511084337E3 * t94 * t73 * t5 -
        0.3259684839629149E4 * t214 * t16 +
        0.577096531476163E2 * t121 * t240 * t16 +
        0.1885128832188574E4 * t194 * t16 + 0.2143715487584216E3 * t276 * t23 -
        0.3770257664377148E4 * t39 * t183 * t16;
      double t288 = t7 * t20;
      double t317 = t86 * t83;
      double t325 = t26 * t123;
      double t326 = t325 * t122;
      double t330 = t76 * t8 * t75;
      double t334 = t41 * t34;
      double t337 =
        -0.2876933760094519E3 * t26 * t98 * t15 +
        0.1740910333200721E4 * t8 * t288 * t16 -
        0.2876933760094519E3 * t26 * t127 * t15 +
        0.1885128832188574E4 * t40 * t39 * t111 * t16 +
        0.8704551666003607E3 * t11 * t21 -
        0.3259684839629149E4 * t86 * t85 * t82 * t16 +
        0.6010906701301438E3 * t163 * t5 -
        0.9354678740125828E3 * t155 * t152 * t5 + 0.8395020705215363E2 * t23 -
        0.3259684839629149E4 * t86 * t84 * t16 -
        0.1819532671866223E4 * t97 * t127 * t15 -
        0.6648285590338501E4 * t317 * t186 * t128 +
        0.1740910333200721E4 * t10 * t8 * t20 * t16 +
        0.486070407931927E1 * t326 * t176 - 0.1754475502216867E4 * t330 * t16 -
        0.8395020705215363E2 * t5 - 0.1411445924867677E4 * t334 * t31;
      double t352 = t73 * t20;
      double t386 =
        0.1411445924867677E4 * t334 * t229 -
        0.3770257664377148E4 * t40 * t183 * t16 -
        0.486070407931927E1 * t326 * t119 -
        0.1629842419814574E4 * t214 * t82 * t5 +
        0.2470731514166843E4 * t114 * t7 * t111 * t16 +
        0.1117159024177316E4 * t94 * t352 * t5 -
        0.6648285590338501E4 * t213 * t186 * t128 +
        0.1411445924867677E4 * t334 * t16 -
        0.1870935748025166E4 * t155 * t62 * t4 * t15 +
        0.3324142795169251E4 * t88 * t184 -
        0.3259684839629149E4 * t85 * t84 * t16 +
        0.1900815597888194E4 * t40 * t217 * t15 +
        0.486070407931927E1 * t326 * t16 - 0.338387630262767E3 * t276 * t16 +
        0.2234318048354632E4 * t256 * t72 +
        0.1202181340260288E4 * t40 * t71 * t15 +
        0.3532767025011433E4 * t94 * t73 * t4 * t15;
      double t388 = t111 * t23;
      double t397 = t30 * t4;
      double t436 = t82 * t20;
      double t440 =
        0.2470731514166843E4 * t115 * t388 + 0.2234318048354632E4 * t94 * t72 +
        0.2143715487584216E3 * t33 * t32 * t4 * t15 +
        0.2143715487584216E3 * t33 * t397 * t15 -
        0.1754475502216867E4 * t79 * t16 +
        0.486070407931927E1 * t325 * t262 * t16 +
        0.2143715487584216E3 * t32 * t397 * t15 +
        0.8772377511084337E3 * t94 * t73 * t23 -
        0.2187016565401361E4 * t10 * t8 * t4 * t15 +
        0.486070407931927E1 * t325 * t240 * t16 -
        0.2876933760094519E3 * t224 * t23 +
        0.1740910333200721E4 * t10 * t288 * t16 -
        0.2187016565401361E4 * t10 * t246 * t15 +
        0.8092719870692021E3 * t88 * t21 +
        0.1618543974138404E4 * t87 * t83 * t20 * t16 -
        0.1117159024177316E4 * t94 * t352 * t23 +
        0.1618543974138404E4 * t87 * t436 * t16;
      double t486 =
        0.1885128832188574E4 * t194 * t388 +
        0.1618543974138404E4 * t317 * t436 * t16 -
        0.2143715487584216E3 * t276 * t5 + 0.2234318048354632E4 * t330 * t72 +
        0.109350828270068E4 * t11 * t23 - 0.1870935748025166E4 * t155 * t61 +
        0.1618543974138404E4 * t213 * t436 * t16 -
        0.4941463028333687E4 * t114 * t183 * t16 -
        0.5350378218286328E2 * t26 * t38 * t16 +
        0.1411445924867677E4 * t41 * t107 * t16 -
        0.4312196439030875E1 * t56 * t49 * t23 -
        0.4941463028333687E4 * t10 * t7 * t183 * t16 -
        0.5350378218286328E2 * t35 * t16 +
        0.486070407931927E1 * t26 * t121 * t240 * t16 +
        0.2234318048354632E4 * t132 * t72 +
        0.1202181340260288E4 * t39 * t71 * t15 -
        0.3324142795169251E4 * t88 * t235;
      double t489 = t106 + t179 + t232 + t282 + t337 + t386 + t440 + t486;
      return t489;
    }

    // * f52 *********************************************************************

    double
    ortho2_f52 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t9 = -t3 - t8 + 0.1546536707079771;
      double t10 = -t3 - t8 - 0.5;
      double t12 = -t3 - t8 - 0.1154653670707977E1;
      double t13 = t12 * t10 * t9;
      double t16 = 0.1E1 * x;
      double t17 = 0.1154653670707977E1 + t16 + t1;
      double t18 = 0.5 + t16 + t1;
      double t19 = t18 * t17;
      double t20 = -0.1546536707079771 + t16 + t1;
      double t24 = 0.9472135954999579 + t16 + t1;
      double t25 = 0.5278640450004206E-1 + t16 + t1;
      double t26 = t25 * t24;
      double t29 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t35 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t36 = -t3 - t8 - 0.5278640450004206E-1;
      double t38 = -t3 - t8 - 0.9472135954999579;
      double t43 = t24 * t6;
      double t51 = t35 * t6;
      double t69 = t51 * t5;
      double t70 = -t3 - t8 + 0.3302238962785669;
      double t71 = -t3 - t8 - 0.3115120652928579E-1;
      double t73 = -t3 - t8 - 0.9688487934707142;
      double t75 = -t3 - t8 - 0.1330223896278567E1;
      double t81 = t17 * t6 * t5;
      double t82 = t20 * t18;
      double t86 = t43 * t5;
      double t87 = -t3 - t8 + 0.2650553239294647;
      double t89 = -t3 - t8 - 0.2147684835193549;
      double t90 = -t3 - t8 - 0.7852315164806451;
      double t92 = -t3 - t8 - 0.1265055323929465E1;
      double t97 =
        0.1820068381573848E4 * t13 * t7 + 0.8809646492467098E3 * t20 * t19 * t7 -
        0.1101072487729106E4 * t29 * t26 * t7 -
        0.1297500910708003E4 * t38 * t36 * t35 * t7 + 0.2082207423944037E3 * t7 +
        0.8981261799019269E3 * t25 * t43 * t5 +
        0.1848506869817779E4 * t38 * t36 * t6 * t5 +
        0.2409173318021205E3 * t29 * t51 * t5 + 0.1250923670061809E4 * (-t3 - t8 -
                        0.1371740148509607E1)
        * (-t3 - t8 - 0.1091700181433142E1) * (-t3 - t8 -
                 0.7092992179024789) * (-t3 - t8 -
                      0.2907007820975211)
        * (-t3 - t8 + 0.917001814331423E-1) * (-t3 - t8 +
                 0.3717401485096066) * t6 * t5 -
        0.2300613577316539E4 * t75 * t73 * t10 * t71 * t70 * t69 +
        0.8710409471020156E4 * t13 * t82 * t81 +
        0.139327173855154E5 * t92 * t90 * t89 * t87 * t25 * t86;
      double t98 = t38 * t36;
      double t102 = 0.1265055323929465E1 + t16 + t1;
      double t103 = 0.7852315164806451 + t16 + t1;
      double t105 = 0.2147684835193549 + t16 + t1;
      double t106 = -0.2650553239294647 + t16 + t1;
      double t117 = t92 * t90 * t89 * t87;
      double t121 = t12 * t10;
      double t163 =
        0.7664414936237447E4 * t98 * t26 * t7 +
        0.3996818745482428E2 * t106 * t105 * t103 * t102 * t7 -
        0.8855282029741453E3 * t29 * t20 * t19 * t7 +
        0.3773488637506035E4 * t117 * t7 -
        0.2816878062162366E4 * t121 * t9 * t35 * t7 -
        0.2267874668403667E3 * t29 * t6 * t5 - 0.1928061474454036E3 * t69 +
        0.194591173823294E1 * (-0.3302238962785669 + t16 +
             t1) * (0.3115120652928579E-1 + t16 +
              t1) * t18 * (0.9688487934707142 + t16 +
                     t1) * (0.1330223896278567E1 +
                      t16 + t1) * t6 * t5 +
        0.6454549858382595E4 * t121 * t9 * t25 * t86 -
        0.130681915980417E4 * t117 * t69 -
        0.605806373058894E3 * t29 * t106 * t105 * t103 * t102 * t6 * t5 +
        0.4621342875099631E4 * t98 * t82 * t81 +
        0.8586065532091051E3 * t75 * t73 * t10 * t71 * t70 * t6 * t5;
      double t164 = t97 + t163;
      return t164;
    }

    double
    ortho2_f52x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.1154653670707977E1 + t8 + t1;
      double t10 = 0.5 + t8 + t1;
      double t11 = t10 * t9;
      double t14 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t18 = -0.1546536707079771 + t8 + t1;
      double t19 = t18 * t11;
      double t22 = 0.1E1 * y;
      double t23 = -t3 - t22 + 0.1546536707079771;
      double t24 = -t3 - t22 - 0.5;
      double t25 = t24 * t23;
      double t26 = -t3 - t22 - 0.1154653670707977E1;
      double t27 = t26 * t25;
      double t30 = t6 * t2;
      double t31 = -t3 - t22 + 0.2650553239294647;
      double t33 = -t3 - t22 - 0.2147684835193549;
      double t34 = -t3 - t22 - 0.7852315164806451;
      double t35 = t34 * t33;
      double t36 = -t3 - t22 - 0.1265055323929465E1;
      double t37 = t36 * t35;
      double t49 = t18 * t9;
      double t53 = t18 * t10;
      double t54 = t14 * t53;
      double t57 = t33 * t31;
      double t66 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t67 = -t3 - t22 - 0.5278640450004206E-1;
      double t69 = -t3 - t22 - 0.9472135954999579;
      double t70 = t69 * t67 * t66;
      double t75 = t66 * t6;
      double t79 = t9 * t30;
      double t85 = t9 * t5;
      double t88 =
        -0.8855282029741453E3 * t14 * t11 * t7 - 0.1400143026857095E4 * t19 * t7 +
        0.8907750567394444E4 * t27 * t7 - 0.1886744318753018E4 * t37 * t31 * t30 -
        0.1886744318753018E4 * t36 * t34 * t31 * t7 -
        0.1886744318753018E4 * t37 * t7 + 0.1886744318753018E4 * t37 * t31 * t5 -
        0.8855282029741453E3 * t14 * t49 * t7 - 0.8855282029741453E3 * t54 * t7 -
        0.1886744318753018E4 * t34 * t57 * t7 -
        0.1886744318753018E4 * t36 * t57 * t7 + 0.6487504553540017E3 * t70 * t30 -
        0.6487504553540017E3 * t70 * t5 + 0.6487504553540017E3 * t69 * t75 * t5 +
        0.4427641014870726E3 * t54 * t79 + 0.6487504553540017E3 * t67 * t75 * t5 -
        0.4427641014870726E3 * t54 * t85;
      double t89 = t66 * t30;
      double t92 = t66 * t5;
      double t99 = t23 * t66;
      double t106 = 0.9472135954999579 + t8 + t1;
      double t107 = t106 * t30;
      double t108 = 0.5278640450004206E-1 + t8 + t1;
      double t110 = t69 * t67 * t108;
      double t113 = t106 * t5;
      double t122 = t108 * t106;
      double t129 = 0.1265055323929465E1 + t8 + t1;
      double t130 = t129 * t30;
      double t131 = 0.7852315164806451 + t8 + t1;
      double t132 = 0.2147684835193549 + t8 + t1;
      double t133 = t132 * t131;
      double t134 = -0.2650553239294647 + t8 + t1;
      double t135 = t134 * t133;
      double t138 = t129 * t5;
      double t143 = t131 * t129;
      double t147 = t132 * t129;
      double t154 = t36 * t34;
      double t155 = t154 * t57;
      double t158 =
        0.1408439031081183E4 * t27 * t89 - 0.1408439031081183E4 * t27 * t92 +
        0.1408439031081183E4 * t26 * t24 * t66 * t7 +
        0.1408439031081183E4 * t26 * t99 * t7 +
        0.1408439031081183E4 * t24 * t99 * t7 -
        0.3832207468118723E4 * t110 * t107 + 0.3832207468118723E4 * t110 * t113 +
        0.7664414936237447E4 * t69 * t67 * t106 * t7 +
        0.7664414936237447E4 * t110 * t7 -
        0.3832207468118723E4 * t69 * t122 * t7 -
        0.3832207468118723E4 * t67 * t122 * t7 -
        0.1998409372741214E2 * t135 * t130 + 0.1998409372741214E2 * t135 * t138 +
        0.3996818745482428E2 * t135 * t7 +
        0.3996818745482428E2 * t134 * t143 * t7 +
        0.3996818745482428E2 * t134 * t147 * t7 +
        0.3996818745482428E2 * t132 * t143 * t7 +
        0.4132525034928739E4 * t155 * t7;
      double t160 = t69 * t67;
      double t161 = t160 * t53;
      double t182 = -t3 - t22 + 0.3302238962785669;
      double t184 = -t3 - t22 - 0.9688487934707142;
      double t185 = -t3 - t22 - 0.1330223896278567E1;
      double t186 = t185 * t184;
      double t187 = t186 * t24 * t182;
      double t190 = -t3 - t22 - 0.3115120652928579E-1;
      double t192 = t186 * t24 * t190;
      double t205 = t190 * t182;
      double t206 = t186 * t205;
      double t209 = 0.1330223896278567E1 + t8 + t1;
      double t210 = 0.9688487934707142 + t8 + t1;
      double t211 = t210 * t209;
      double t212 = 0.3115120652928579E-1 + t8 + t1;
      double t217 = -0.3302238962785669 + t8 + t1;
      double t222 = t217 * t212;
      double t231 = t222 * t10 * t210;
      double t234 =
        0.4621342875099631E4 * t161 * t7 + 0.2310671437549815E4 * t161 * t85 -
        0.2310671437549815E4 * t161 * t79 +
        0.4621342875099631E4 * t160 * t49 * t7 +
        0.4621342875099631E4 * t160 * t11 * t7 -
        0.2310671437549815E4 * t67 * t18 * t11 * t7 -
        0.2310671437549815E4 * t69 * t18 * t11 * t7 -
        0.4293032766045525E3 * t187 * t7 - 0.4293032766045525E3 * t192 * t7 +
        0.4293032766045525E3 * t192 * t182 * t5 -
        0.4293032766045525E3 * t192 * t182 * t30 -
        0.9578639799558858E3 * t134 * t132 * t143 * t7 -
        0.4293032766045525E3 * t206 * t7 +
        0.194591173823294E1 * t212 * t10 * t211 * t7 +
        0.194591173823294E1 * t217 * t10 * t211 * t7 +
        0.194591173823294E1 * t222 * t211 * t7 +
        0.194591173823294E1 * t222 * t10 * t209 * t7 +
        0.194591173823294E1 * t231 * t7;
      double t242 = t184 * t24 * t205;
      double t246 = t185 * t24 * t205;
      double t249 = t31 * t66;
      double t253 = t36 * t33;
      double t269 = t26 * t24;
      double t270 = t269 * t23 * t108;
      double t280 = t26 * t23;
      double t295 =
        0.9729558691164699 * t231 * t209 * t5 -
        0.9729558691164699 * t231 * t209 * t30 -
        0.4293032766045525E3 * t242 * t7 - 0.4293032766045525E3 * t246 * t7 +
        0.6534095799020851E3 * t35 * t249 * t7 +
        0.6534095799020851E3 * t253 * t249 * t7 +
        0.6534095799020851E3 * t154 * t249 * t7 +
        0.6534095799020851E3 * t154 * t33 * t66 * t7 -
        0.6534095799020851E3 * t155 * t92 + 0.6534095799020851E3 * t155 * t89 +
        0.6454549858382595E4 * t270 * t7 + 0.3227274929191298E4 * t270 * t113 -
        0.3227274929191298E4 * t270 * t107 -
        0.3227274929191298E4 * t25 * t122 * t7 -
        0.3227274929191298E4 * t280 * t122 * t7 -
        0.3227274929191298E4 * t269 * t122 * t7 +
        0.6454549858382595E4 * t269 * t23 * t106 * t7 -
        0.605806373058894E3 * t14 * t132 * t143 * t7;
      double t298 = t14 * t134;
      double t305 = t298 * t133;
      double t321 = t269 * t23 * t18;
      double t327 = t9 * t6;
      double t328 = t327 * t5;
      double t331 = t10 * t6;
      double t335 = -t3 - t22 + 0.3717401485096066;
      double t336 = -t3 - t22 + 0.917001814331423E-1;
      double t337 = t336 * t335;
      double t339 = -t3 - t22 - 0.2907007820975211;
      double t340 = -t3 - t22 - 0.7092992179024789;
      double t342 = -t3 - t22 - 0.1091700181433142E1;
      double t343 = -t3 - t22 - 0.1371740148509607E1;
      double t344 = t343 * t342;
      double t345 = t344 * t340 * t339;
      double t359 =
        -0.605806373058894E3 * t298 * t143 * t7 -
        0.605806373058894E3 * t298 * t147 * t7 - 0.605806373058894E3 * t305 * t7 -
        0.302903186529447E3 * t305 * t138 + 0.302903186529447E3 * t305 * t130 +
        0.1133937334201833E3 * t14 * t30 - 0.9640307372270178E2 * t92 +
        0.251124102811996E3 * t7 - 0.1133937334201833E3 * t14 * t5 +
        0.9640307372270178E2 * t89 - 0.4355204735510078E4 * t321 * t11 * t30 +
        0.4355204735510078E4 * t321 * t11 * t5 +
        0.8710409471020156E4 * t321 * t328 +
        0.8710409471020156E4 * t321 * t331 * t5 +
        0.6254618350309046E3 * t345 * t337 * t5 -
        0.6254618350309046E3 * t345 * t337 * t30 +
        0.8710409471020156E4 * t269 * t23 * t10 * t328 +
        0.727517892022815E4 * t192 * t182 * t6 * t5;
      double t374 = t335 * t6 * t5;
      double t381 = t339 * t336;
      double t393 = t182 * t66;
      double t400 = t75 * t5;
      double t411 = t14 * t66;
      double t416 =
        -0.6254618350309046E3 * t345 * t336 * t6 * t5 -
        0.4355204735510078E4 * t25 * t53 * t328 -
        0.4355204735510078E4 * t280 * t53 * t328 -
        0.4355204735510078E4 * t269 * t53 * t328 -
        0.6254618350309046E3 * t344 * t340 * t336 * t374 -
        0.6254618350309046E3 * t345 * t374 -
        0.6254618350309046E3 * t344 * t381 * t374 -
        0.6254618350309046E3 * t342 * t340 * t381 * t374 -
        0.6254618350309046E3 * t343 * t340 * t381 * t374 -
        0.115030678865827E4 * t192 * t393 * t5 +
        0.115030678865827E4 * t192 * t393 * t30 +
        0.115030678865827E4 * t206 * t400 + 0.115030678865827E4 * t187 * t400 +
        0.115030678865827E4 * t192 * t400 + 0.115030678865827E4 * t242 * t400 +
        0.115030678865827E4 * t246 * t400 + 0.1204586659010603E3 * t411 * t5 -
        0.1204586659010603E3 * t411 * t30;
      double t422 = t67 * t6;
      double t432 = t122 * t5;
      double t435 = t122 * t30;
      double t438 = t106 * t6;
      double t439 = t438 * t5;
      double t446 = t108 * t6;
      double t447 = t446 * t5;
      double t451 = t31 * t108;
      double t464 =
        0.3809237481526094E3 * t400 - 0.7618474963052188E3 * t14 * t6 * t5 -
        0.9242534349088893E3 * t422 * t5 - 0.9242534349088893E3 * t69 * t6 * t5 +
        0.9242534349088893E3 * t160 * t5 - 0.9242534349088893E3 * t160 * t30 +
        0.6966358692757702E4 * t155 * t432 - 0.6966358692757702E4 * t155 * t435 -
        0.6966358692757702E4 * t154 * t33 * t108 * t439 +
        0.139327173855154E5 * t155 * t439 + 0.139327173855154E5 * t155 * t447 -
        0.4490630899509634E3 * t435 - 0.6966358692757702E4 * t35 * t451 * t439 +
        0.4490630899509634E3 * t432 + 0.8981261799019269E3 * t447 -
        0.6966358692757702E4 * t253 * t451 * t439 -
        0.6966358692757702E4 * t154 * t451 * t439 + 0.8981261799019269E3 * t439;
      double t473 = t14 * t122;
      double t490 = t23 * t6;
      double t512 = -0.1041103711972018E3 * t30 + 0.1041103711972018E3 * t5
        - 0.1740948465085887E4 * t108 * t438 * t5 +
        0.4103058143980046E4 * t69 * t422 * t5 +
        0.5505362438645531E3 * t473 * t30 - 0.5505362438645531E3 * t473 * t5 -
        0.1101072487729106E4 * t14 * t446 * t5 -
        0.1101072487729106E4 * t14 * t438 * t5 -
        0.4404823246233549E3 * t19 * t30 -
        0.910034190786924E3 * t26 * t24 * t6 * t5 -
        0.910034190786924E3 * t26 * t490 * t5 - 0.910034190786924E3 * t27 * t30 +
        0.910034190786924E3 * t27 * t5 - 0.910034190786924E3 * t24 * t490 * t5 +
        0.4404823246233549E3 * t19 * t5 + 0.8809646492467098E3 * t18 * t331 * t5 +
        0.8809646492467098E3 * t18 * t327 * t5 +
        0.8809646492467098E3 * t10 * t327 * t5;
      double t515 = t88 + t158 + t234 + t295 + t359 + t416 + t464 + t512;
      return t515;
    }

    double
    ortho2_f52y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t8 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t9 = t8 * t5;
      double t11 = -t3 - t1;
      double t12 = t11 * t2;
      double t13 = t4 * t12;
      double t15 = t4 * t11;
      double t18 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t23 = t8 * t15;
      double t25 = 0.1E1 * x;
      double t26 = 0.1265055323929465E1 + t25 + t1;
      double t27 = t26 * t15;
      double t28 = 0.7852315164806451 + t25 + t1;
      double t29 = 0.2147684835193549 + t25 + t1;
      double t30 = t29 * t28;
      double t31 = -0.2650553239294647 + t25 + t1;
      double t32 = t18 * t31;
      double t33 = t32 * t30;
      double t36 = 0.1E1 * y;
      double t37 = -t3 - t36 + 0.3302238962785669;
      double t39 = -t3 - t36 - 0.3115120652928579E-1;
      double t40 = -t3 - t36 - 0.5;
      double t42 = -t3 - t36 - 0.9688487934707142;
      double t43 = -t3 - t36 - 0.1330223896278567E1;
      double t44 = t43 * t42;
      double t45 = t44 * t40 * t39;
      double t48 = 0.9472135954999579 + t25 + t1;
      double t49 = 0.5278640450004206E-1 + t25 + t1;
      double t50 = t49 * t48;
      double t51 = t50 * t15;
      double t52 = -t3 - t36 + 0.2650553239294647;
      double t53 = -t3 - t36 - 0.2147684835193549;
      double t54 = t53 * t52;
      double t55 = -t3 - t36 - 0.7852315164806451;
      double t56 = -t3 - t36 - 0.1265055323929465E1;
      double t57 = t56 * t55;
      double t58 = t57 * t54;
      double t61 = t39 * t37;
      double t63 = t42 * t40 * t61;
      double t66 = 0.1154653670707977E1 + t25 + t1;
      double t67 = t66 * t4;
      double t68 = t67 * t12;
      double t69 = 0.5 + t25 + t1;
      double t70 = -0.1546536707079771 + t25 + t1;
      double t71 = t70 * t69;
      double t72 = -t3 - t36 - 0.1154653670707977E1;
      double t73 = t72 * t40;
      double t77 = -t3 - t36 + 0.1546536707079771;
      double t82 = t40 * t77;
      double t83 = t72 * t82;
      double t88 = t28 * t26;
      double t94 = -t3 - t36 - 0.5278640450004206E-1;
      double t95 = -t3 - t36 - 0.9472135954999579;
      double t96 = t95 * t94;
      double t102 = 0.9640307372270178E2 * t9 - 0.4123116535906034E3 * t13
        - 0.1133937334201833E3 * t18 * t15 + 0.1133937334201833E3 * t18 * t5 -
        0.9640307372270178E2 * t23 - 0.302903186529447E3 * t33 * t27 -
        0.4293032766045525E3 * t45 * t37 * t5 + 0.6966358692757702E4 * t58 * t51 -
        0.8586065532091051E3 * t63 * t13 -
        0.8710409471020156E4 * t73 * t71 * t68 +
        0.4355204735510078E4 * t73 * t77 * t69 * t68 -
        0.1408439031081183E4 * t83 * t23 + 0.910034190786924E3 * t83 * t15 -
        0.302903186529447E3 * t18 * t29 * t88 * t13 + 0.4490630899509634E3 * t51 -
        0.9242534349088893E3 * t96 * t5 - 0.1848506869817779E4 * t95 * t4 * t12;
      double t103 = t94 * t4;
      double t111 = t8 * t4;
      double t112 = t111 * t12;
      double t114 = t18 * t8;
      double t117 = t31 * t30;
      double t120 = t55 * t53;
      double t121 = t56 * t120;
      double t129 = t66 * t15;
      double t130 = t18 * t71;
      double t134 = t43 * t40 * t61;
      double t137 = t50 * t5;
      double t144 = t52 * t8;
      double t149 = t44 * t40 * t37;
      double t152 = t72 * t77;
      double t156 = t96 * t71;
      double t162 =
        -0.1848506869817779E4 * t103 * t12 + 0.9242534349088893E3 * t96 * t15 -
        0.3809237481526094E3 * t18 * t4 * t12 + 0.7618474963052188E3 * t112 -
        0.1204586659010603E3 * t114 * t5 + 0.1998409372741214E2 * t117 * t27 -
        0.3773488637506035E4 * t121 * t13 - 0.8586065532091051E3 * t45 * t13 -
        0.7664414936237447E4 * t94 * t50 * t13 -
        0.4427641014870726E3 * t130 * t129 + 0.2300613577316539E4 * t134 * t112 -
        0.6966358692757702E4 * t58 * t137 -
        0.3773488637506035E4 * t56 * t55 * t52 * t13 +
        0.130681915980417E4 * t120 * t144 * t13 -
        0.8586065532091051E3 * t149 * t13 -
        0.8710409471020156E4 * t152 * t71 * t68 +
        0.2310671437549815E4 * t156 * t129 -
        0.3773488637506035E4 * t56 * t54 * t13;
      double t174 = t48 * t5;
      double t176 = t73 * t77 * t49;
      double t186 = t69 * t66;
      double t187 = t70 * t186;
      double t191 = t48 * t4;
      double t192 = t191 * t12;
      double t194 = t49 * t4;
      double t195 = t194 * t12;
      double t197 = 0.1330223896278567E1 + t25 + t1;
      double t198 = 0.9688487934707142 + t25 + t1;
      double t199 = t198 * t197;
      double t200 = 0.3115120652928579E-1 + t25 + t1;
      double t201 = -0.3302238962785669 + t25 + t1;
      double t202 = t201 * t200;
      double t217 = t202 * t69 * t198;
      double t220 = t52 * t49;
      double t224 =
        -0.1886744318753018E4 * t121 * t52 * t5 +
        0.1204586659010603E3 * t114 * t15 + 0.1408439031081183E4 * t83 * t9 +
        0.4293032766045525E3 * t45 * t37 * t15 -
        0.3227274929191298E4 * t176 * t174 + 0.2300613577316539E4 * t63 * t112 -
        0.8710409471020156E4 * t82 * t71 * t68 +
        0.3227274929191298E4 * t176 * t13 - 0.280028605371419E4 * t187 * t13 -
        0.4490630899509634E3 * t137 + 0.4490630899509634E3 * t192 +
        0.4490630899509634E3 * t195 + 0.9729558691164699 * t202 * t199 * t13 +
        0.2816878062162366E4 * t72 * t40 * t8 * t13 -
        0.3773488637506035E4 * t55 * t54 * t13 -
        0.4404823246233549E3 * t187 * t5 +
        0.9729558691164699 * t217 * t197 * t15 -
        0.139327173855154E5 * t120 * t220 * t192;
      double t225 = t44 * t61;
      double t228 = t37 * t8;
      double t234 = t73 * t77 * t70;
      double t250 = t70 * t66;
      double t256 = t77 * t8;
      double t261 = t95 * t94 * t49;
      double t267 = t26 * t5;
      double t270 = -t3 - t36 + 0.3717401485096066;
      double t272 = t270 * t4 * t12;
      double t273 = -t3 - t36 + 0.917001814331423E-1;
      double t274 = -t3 - t36 - 0.7092992179024789;
      double t276 = -t3 - t36 - 0.1091700181433142E1;
      double t277 = -t3 - t36 - 0.1371740148509607E1;
      double t278 = t277 * t276;
      double t282 = -t3 - t36 - 0.2907007820975211;
      double t284 = t278 * t274 * t282;
      double t291 = t77 * t4;
      double t297 =
        -0.8586065532091051E3 * t225 * t13 -
        0.115030678865827E4 * t45 * t228 * t15 +
        0.4355204735510078E4 * t234 * t186 * t15 -
        0.1820068381573848E4 * t72 * t40 * t4 * t12 -
        0.1915727959911772E4 * t31 * t29 * t88 * t13 -
        0.4427641014870726E3 * t130 * t13 -
        0.3481896930171775E4 * t49 * t191 * t12 -
        0.4427641014870726E3 * t18 * t250 * t13 -
        0.6534095799020851E3 * t58 * t23 +
        0.2816878062162366E4 * t72 * t256 * t13 -
        0.3832207468118723E4 * t261 * t174 -
        0.4427641014870726E3 * t18 * t186 * t13 -
        0.1998409372741214E2 * t117 * t267 -
        0.1250923670061809E4 * t278 * t274 * t273 * t272 -
        0.1250923670061809E4 * t284 * t272 +
        0.9729558691164699 * t201 * t69 * t199 * t13 -
        0.1820068381573848E4 * t72 * t291 * t12 +
        0.3832207468118723E4 * t261 * t13;
      double t304 = t69 * t4;
      double t330 = t18 * t50;
      double t335 = t273 * t270;
      double t339 = t282 * t273;
      double t343 = t48 * t15;
      double t353 =
        0.3832207468118723E4 * t95 * t94 * t48 * t13 +
        0.4404823246233549E3 * t70 * t304 * t12 +
        0.4404823246233549E3 * t70 * t67 * t12 +
        0.6966358692757702E4 * t58 * t195 + 0.6966358692757702E4 * t58 * t192 +
        0.2051529071990023E4 * t95 * t103 * t12 +
        0.4404823246233549E3 * t69 * t67 * t12 - 0.910034190786924E3 * t83 * t5 +
        0.9729558691164699 * t200 * t69 * t199 * t13 -
        0.4355204735510078E4 * t234 * t186 * t5 -
        0.5505362438645531E3 * t330 * t15 + 0.5505362438645531E3 * t330 * t5 +
        0.6254618350309046E3 * t284 * t335 * t15 -
        0.1250923670061809E4 * t278 * t339 * t272 +
        0.3832207468118723E4 * t261 * t343 + 0.2066262517464369E4 * t58 * t13 +
        0.1886744318753018E4 * t121 * t52 * t15 + 0.6534095799020851E3 * t58 * t9;
      double t365 = t29 * t26;
      double t405 = t95 * t94 * t8;
      double t411 =
        -0.1820068381573848E4 * t40 * t291 * t12 -
        0.7664414936237447E4 * t95 * t50 * t13 +
        0.4355204735510078E4 * t234 * t304 * t12 +
        0.1998409372741214E2 * t117 * t13 +
        0.1998409372741214E2 * t31 * t365 * t13 -
        0.1250923670061809E4 * t277 * t274 * t339 * t272 +
        0.4453875283697222E4 * t83 * t13 +
        0.1998409372741214E2 * t31 * t88 * t13 +
        0.4404823246233549E3 * t187 * t15 +
        0.3227274929191298E4 * t73 * t77 * t48 * t13 -
        0.9729558691164699 * t217 * t197 * t5 -
        0.1250923670061809E4 * t276 * t274 * t339 * t272 -
        0.6454549858382595E4 * t73 * t50 * t13 -
        0.139327173855154E5 * t57 * t53 * t49 * t192 +
        0.2816878062162366E4 * t40 * t256 * t13 +
        0.1998409372741214E2 * t29 * t88 * t13 +
        0.6487504553540017E3 * t405 * t5 +
        0.1297500910708003E4 * t95 * t111 * t12;
      double t420 = t66 * t5;
      double t462 =
        -0.6254618350309046E3 * t284 * t335 * t5 +
        0.4355204735510078E4 * t234 * t68 - 0.8586065532091051E3 * t134 * t13 -
        0.2310671437549815E4 * t156 * t420 -
        0.139327173855154E5 * t57 * t220 * t192 +
        0.302903186529447E3 * t33 * t267 + 0.2310671437549815E4 * t156 * t13 -
        0.302903186529447E3 * t33 * t13 + 0.115030678865827E4 * t45 * t228 * t5 +
        0.2300613577316539E4 * t45 * t112 -
        0.6454549858382595E4 * t152 * t50 * t13 -
        0.1250923670061809E4 * t284 * t273 * t4 * t12 +
        0.130681915980417E4 * t57 * t53 * t8 * t13 +
        0.2310671437549815E4 * t96 * t250 * t13 -
        0.5505362438645531E3 * t18 * t194 * t12 +
        0.2300613577316539E4 * t149 * t112 +
        0.1297500910708003E4 * t94 * t111 * t12 -
        0.302903186529447E3 * t32 * t365 * t13;
      double t486 = t56 * t53;
      double t513 =
        -0.5505362438645531E3 * t18 * t191 * t12 +
        0.4427641014870726E3 * t130 * t420 + 0.2300613577316539E4 * t225 * t112 +
        0.2310671437549815E4 * t96 * t186 * t13 + 0.1041103711972018E3 * t15 -
        0.6454549858382595E4 * t82 * t50 * t13 +
        0.130681915980417E4 * t57 * t144 * t13 + 0.9729558691164699 * t217 * t13 +
        0.9729558691164699 * t202 * t69 * t197 * t13 -
        0.139327173855154E5 * t486 * t220 * t192 -
        0.302903186529447E3 * t32 * t88 * t13 -
        0.4621342875099631E4 * t95 * t70 * t186 * t13 -
        0.1041103711972018E3 * t5 + 0.3227274929191298E4 * t176 * t343 +
        0.3637589460114075E4 * t45 * t37 * t4 * t12 +
        0.130681915980417E4 * t486 * t144 * t13 -
        0.4621342875099631E4 * t94 * t70 * t186 * t13 -
        0.6487504553540017E3 * t405 * t15;
      double t516 = t102 + t162 + t224 + t297 + t353 + t411 + t462 + t513;
      return t516;
    }

    // * f53 *********************************************************************

    double
    ortho2_f53 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t13 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t17 = 0.1E1 * y;
      double t18 = -t3 - t17 - 0.5278640450004206E-1;
      double t20 = -t3 - t17 - 0.9472135954999579;
      double t24 = 0.1E1 * x;
      double t25 = 0.9472135954999579 + t24 + t1;
      double t26 = t25 * t6;
      double t27 = 0.5278640450004206E-1 + t24 + t1;
      double t31 = t6 * t5;
      double t32 = t27 * t25;
      double t40 = -t3 - t17 + 0.1546536707079771;
      double t41 = -t3 - t17 - 0.5;
      double t43 = -t3 - t17 - 0.1154653670707977E1;
      double t44 = t43 * t41 * t40;
      double t47 = 0.1154653670707977E1 + t24 + t1;
      double t48 = 0.5 + t24 + t1;
      double t49 = t48 * t47;
      double t50 = -0.1546536707079771 + t24 + t1;
      double t56 = t43 * t41;
      double t60 = t20 * t18;
      double t64 = 0.1265055323929465E1 + t24 + t1;
      double t65 = 0.7852315164806451 + t24 + t1;
      double t67 = 0.2147684835193549 + t24 + t1;
      double t68 = -0.2650553239294647 + t24 + t1;
      double t77 = -t3 - t17 + 0.2650553239294647;
      double t78 = -t3 - t17 - 0.2147684835193549;
      double t80 = -t3 - t17 - 0.7852315164806451;
      double t81 = -t3 - t17 - 0.1265055323929465E1;
      double t83 = t81 * t80 * t78 * t77;
      double t86 =
        0.2982166719335347E3 * t13 * t10 * t5 +
        0.1946100011820659E4 * t20 * t18 * t6 * t5 +
        0.1476313010788564E4 * t27 * t26 * t5 -
        0.6861248656465258E3 * t13 * t32 * t31 -
        0.4658436012012737E3 * t20 * t18 * t9 * t31 +
        0.3774290410247878E3 * t44 * t31 +
        0.1094544143263798E4 * t50 * t49 * t31 + 0.2320554691692617E3 * t31 -
        0.3044206802479674E4 * t56 * t40 * t9 * t31 +
        0.1199633016720218E5 * t60 * t32 * t31 +
        0.1130005678384205E4 * t68 * t67 * t65 * t64 * t31 -
        0.2008298687123097E4 * t13 * t50 * t49 * t31 +
        0.2836542630529092E4 * t83 * t31;
      double t88 = t47 * t6 * t5;
      double t89 = t50 * t48;
      double t93 = t26 * t5;
      double t114 = t10 * t5;
      double t115 = -t3 - t17 + 0.3302238962785669;
      double t116 = -t3 - t17 - 0.3115120652928579E-1;
      double t118 = -t3 - t17 - 0.9688487934707142;
      double t120 = -t3 - t17 - 0.1330223896278567E1;
      double t126 = t64 * t6 * t5;
      double t127 = t67 * t65;
      double MapleGenVar1 =
        0.1805304996254271E5 * t44 * t89 * t88 +
        0.1380528777907254E5 * t81 * t80 * t78 * t77 * t27 * t93 +
        0.5528469356472047E3 * (-t3 - t17 - 0.1371740148509607E1) * (-t3 - t17 -
                     0.1091700181433142E1)
        * (-t3 - t17 - 0.7092992179024789) * (-t3 - t17 -
                0.2907007820975211) * (-t3 - t17 +
                     0.917001814331423E-1)
        * (-t3 - t17 + 0.3717401485096066) * t6 * t5 -
        0.1436332603686921E4 * t120 * t118 * t41 * t116 * t115 * t114 +
        0.8475034230283761E4 * t20 * t18 * t68 * t127 * t126 -
        0.1177251798820042E3 * t114;
      double t168 =
        MapleGenVar1 - 0.9363683116433664E2 * t13 * t6 * t5 +
        0.6439929362750927E3 * (-0.3302238962785669 + t24 +
              t1) * (0.3115120652928579E-1 + t24 +
               t1) * t48 * (0.9688487934707142 + t24 +
                t1) * (0.1330223896278567E1 +
                       t24 + t1) * t6 * t5 +
        0.2327982157053983E4 * t56 * t40 * t27 * t93 -
        0.1920773965135279E3 * t83 * t114 -
        0.698768915145629E3 * t13 * t68 * t127 * t126 +
        0.3156600791226379E4 * t60 * t89 * t88 +
        0.3799026905676787E1 * t120 * t118 * t41 * t116 * t115 * t6 * t5;
      double t169 = t86 + t168;
      return t169;
    }

    double
    ortho2_f53x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.1154653670707977E1 + t8 + t1;
      double t10 = 0.5 + t8 + t1;
      double t11 = t10 * t9;
      double t14 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t18 = -0.1546536707079771 + t8 + t1;
      double t19 = t18 * t11;
      double t22 = t18 * t9;
      double t26 = t18 * t10;
      double t27 = t14 * t26;
      double t30 = t6 * t2;
      double t31 = 0.1E1 * y;
      double t32 = -t3 - t31 + 0.2650553239294647;
      double t34 = -t3 - t31 - 0.2147684835193549;
      double t35 = -t3 - t31 - 0.7852315164806451;
      double t36 = t35 * t34;
      double t37 = -t3 - t31 - 0.1265055323929465E1;
      double t38 = t37 * t36;
      double t41 = -t3 - t31 + 0.1546536707079771;
      double t42 = -t3 - t31 - 0.5;
      double t43 = t42 * t41;
      double t44 = -t3 - t31 - 0.1154653670707977E1;
      double t45 = t44 * t43;
      double t48 = t34 * t32;
      double t66 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t67 = t66 * t5;
      double t70 = t66 * t30;
      double t73 = t9 * t5;
      double t76 = t9 * t30;
      double t79 = t41 * t66;
      double t87 = 0.9472135954999579 + t8 + t1;
      double t88 = t87 * t30;
      double t89 = 0.5278640450004206E-1 + t8 + t1;
      double t90 = -t3 - t31 - 0.5278640450004206E-1;
      double t92 = -t3 - t31 - 0.9472135954999579;
      double t93 = t92 * t90 * t89;
      double t96 =
        -0.2008298687123097E4 * t14 * t11 * t7 - 0.3175399036617427E4 * t19 * t7 -
        0.2008298687123097E4 * t14 * t22 * t7 - 0.2008298687123097E4 * t27 * t7 -
        0.1418271315264546E4 * t38 * t32 * t30 + 0.9626627164414088E4 * t45 * t7 -
        0.1418271315264546E4 * t35 * t48 * t7 -
        0.1418271315264546E4 * t37 * t48 * t7 -
        0.1418271315264546E4 * t37 * t35 * t32 * t7 -
        0.1418271315264546E4 * t38 * t7 + 0.1418271315264546E4 * t38 * t32 * t5 -
        0.1522103401239837E4 * t45 * t67 + 0.1522103401239837E4 * t45 * t70 -
        0.1004149343561548E4 * t27 * t73 + 0.1004149343561548E4 * t27 * t76 +
        0.1522103401239837E4 * t44 * t79 * t7 +
        0.1522103401239837E4 * t44 * t42 * t66 * t7 -
        0.5998165083601091E4 * t93 * t88;
      double t100 = t89 * t87;
      double t113 = t87 * t5;
      double t116 = 0.1265055323929465E1 + t8 + t1;
      double t117 = 0.7852315164806451 + t8 + t1;
      double t118 = t117 * t116;
      double t119 = 0.2147684835193549 + t8 + t1;
      double t123 = -0.2650553239294647 + t8 + t1;
      double t127 = t119 * t116;
      double t131 = t119 * t117;
      double t132 = t123 * t131;
      double t135 = t116 * t5;
      double t138 = t116 * t30;
      double t141 = t37 * t35;
      double t142 = t141 * t48;
      double t145 = t92 * t90;
      double t146 = t145 * t26;
      double t156 = -t3 - t31 + 0.3302238962785669;
      double t158 = -t3 - t31 - 0.3115120652928579E-1;
      double t160 = -t3 - t31 - 0.9688487934707142;
      double t161 = -t3 - t31 - 0.1330223896278567E1;
      double t162 = t161 * t160;
      double t163 = t162 * t42 * t158;
      double t169 =
        0.1522103401239837E4 * t42 * t79 * t7 -
        0.5998165083601091E4 * t90 * t100 * t7 -
        0.5998165083601091E4 * t92 * t100 * t7 +
        0.1199633016720218E5 * t92 * t90 * t87 * t7 +
        0.1199633016720218E5 * t93 * t7 + 0.5998165083601091E4 * t93 * t113 +
        0.1130005678384205E4 * t119 * t118 * t7 +
        0.1130005678384205E4 * t123 * t118 * t7 +
        0.1130005678384205E4 * t123 * t127 * t7 +
        0.1130005678384205E4 * t132 * t7 + 0.5650028391921024E3 * t132 * t135 -
        0.5650028391921024E3 * t132 * t138 + 0.6074020600180329E3 * t142 * t7 +
        0.157830039561319E4 * t146 * t73 - 0.157830039561319E4 * t146 * t76 +
        0.3156600791226379E4 * t145 * t22 * t7 +
        0.3156600791226379E4 * t146 * t7 +
        0.1899513452838393E1 * t163 * t156 * t5 -
        0.1899513452838393E1 * t163 * t156 * t30;
      double t171 = t123 * t119;
      double t186 = 0.1330223896278567E1 + t8 + t1;
      double t188 = 0.9688487934707142 + t8 + t1;
      double t190 = 0.3115120652928579E-1 + t8 + t1;
      double t191 = -0.3302238962785669 + t8 + t1;
      double t192 = t191 * t190;
      double t193 = t192 * t10 * t188;
      double t196 = t158 * t156;
      double t197 = t162 * t196;
      double t201 = t162 * t42 * t156;
      double t207 = t160 * t42 * t196;
      double t211 = t161 * t42 * t196;
      double t214 = t188 * t186;
      double t243 =
        -0.1104850664992558E4 * t171 * t118 * t7 -
        0.157830039561319E4 * t90 * t18 * t11 * t7 -
        0.157830039561319E4 * t92 * t18 * t11 * t7 +
        0.3156600791226379E4 * t145 * t11 * t7 -
        0.3219964681375463E3 * t193 * t186 * t30 -
        0.1899513452838393E1 * t197 * t7 - 0.1899513452838393E1 * t201 * t7 -
        0.1899513452838393E1 * t163 * t7 - 0.1899513452838393E1 * t207 * t7 -
        0.1899513452838393E1 * t211 * t7 +
        0.6439929362750927E3 * t190 * t10 * t214 * t7 +
        0.6439929362750927E3 * t191 * t10 * t214 * t7 +
        0.6439929362750927E3 * t192 * t214 * t7 +
        0.6439929362750927E3 * t192 * t10 * t186 * t7 +
        0.6439929362750927E3 * t193 * t7 +
        0.3219964681375463E3 * t193 * t186 * t5 +
        0.9603869825676393E2 * t141 * t34 * t66 * t7 -
        0.9603869825676393E2 * t142 * t67 + 0.9603869825676393E2 * t142 * t70;
      double t244 = t32 * t66;
      double t248 = t37 * t34;
      double t255 = t44 * t41;
      double t259 = t44 * t42;
      double t268 = t259 * t41 * t89;
      double t282 = t14 * t123;
      double t289 = t282 * t131;
      double t301 =
        0.9603869825676393E2 * t36 * t244 * t7 +
        0.9603869825676393E2 * t248 * t244 * t7 +
        0.9603869825676393E2 * t141 * t244 * t7 -
        0.1163991078526991E4 * t255 * t100 * t7 -
        0.1163991078526991E4 * t259 * t100 * t7 +
        0.2327982157053983E4 * t259 * t41 * t87 * t7 +
        0.2327982157053983E4 * t268 * t7 + 0.1163991078526991E4 * t268 * t113 -
        0.1163991078526991E4 * t268 * t88 -
        0.1163991078526991E4 * t43 * t100 * t7 -
        0.698768915145629E3 * t14 * t119 * t118 * t7 -
        0.698768915145629E3 * t282 * t118 * t7 -
        0.698768915145629E3 * t282 * t127 * t7 - 0.698768915145629E3 * t289 * t7 -
        0.3493844575728145E3 * t289 * t135 + 0.3493844575728145E3 * t289 * t138 -
        0.5886258994100211E2 * t67 - 0.4681841558216832E2 * t14 * t5 +
        0.4681841558216832E2 * t14 * t30;
      double t316 = t90 * t6;
      double t319 = t14 * t66;
      double t324 = t66 * t6;
      double t325 = t324 * t5;
      double t329 = t259 * t41 * t18;
      double t335 = t9 * t6;
      double t336 = t335 * t5;
      double t343 = t10 * t6;
      double t347 = t100 * t5;
      double t349 = t89 * t6;
      double t350 = t349 * t5;
      double t352 = t87 * t6;
      double t353 = t352 * t5;
      double t355 = t100 * t30;
      double t357 = 0.5886258994100211E2 * t70 + 0.2242268767001959E3 * t7
        - 0.9730500059103294E3 * t145 * t30 -
        0.9730500059103294E3 * t92 * t6 * t5 + 0.9730500059103294E3 * t145 * t5 -
        0.9430439195451794E3 * t14 * t6 * t5 - 0.9730500059103294E3 * t316 * t5 +
        0.1491083359667674E3 * t319 * t5 - 0.1491083359667674E3 * t319 * t30 +
        0.4715219597725897E3 * t325 + 0.9026524981271357E4 * t329 * t11 * t5 -
        0.9026524981271357E4 * t329 * t11 * t30 +
        0.1805304996254271E5 * t259 * t41 * t10 * t336 +
        0.1805304996254271E5 * t329 * t336 +
        0.1805304996254271E5 * t329 * t343 * t5 + 0.7381565053942821E3 * t347 +
        0.1476313010788564E4 * t350 + 0.1476313010788564E4 * t353 -
        0.7381565053942821E3 * t355;
      double t368 = -t3 - t31 + 0.917001814331423E-1;
      double t371 = -t3 - t31 - 0.2907007820975211;
      double t372 = -t3 - t31 - 0.7092992179024789;
      double t374 = -t3 - t31 - 0.1091700181433142E1;
      double t375 = -t3 - t31 - 0.1371740148509607E1;
      double t376 = t375 * t374;
      double t377 = t376 * t372 * t371;
      double t380 = -t3 - t31 + 0.3717401485096066;
      double t381 = t368 * t380;
      double t392 = t145 * t171;
      double t403 = t116 * t6 * t5;
      double t421 = t156 * t66;
      double t426 = t380 * t6 * t5;
      double t427 = t371 * t368;
      double t439 =
        -0.9026524981271357E4 * t255 * t26 * t336 -
        0.9026524981271357E4 * t259 * t26 * t336 +
        0.4542082505210631E4 * t163 * t156 * t6 * t5 -
        0.2764234678236024E3 * t377 * t368 * t6 * t5 +
        0.2764234678236024E3 * t377 * t381 * t5 -
        0.2764234678236024E3 * t377 * t381 * t30 -
        0.9026524981271357E4 * t43 * t26 * t336 -
        0.4237517115141881E4 * t392 * t118 * t30 +
        0.8475034230283761E4 * t392 * t117 * t6 * t5 +
        0.4237517115141881E4 * t392 * t118 * t5 +
        0.8475034230283761E4 * t145 * t123 * t117 * t403 +
        0.8475034230283761E4 * t392 * t403 -
        0.4237517115141881E4 * t92 * t123 * t131 * t403 +
        0.8475034230283761E4 * t145 * t131 * t403 -
        0.4237517115141881E4 * t90 * t123 * t131 * t403 +
        0.7181663018434603E3 * t163 * t421 * t30 -
        0.2764234678236024E3 * t374 * t372 * t427 * t426 -
        0.2764234678236024E3 * t375 * t372 * t427 * t426 -
        0.2764234678236024E3 * t376 * t427 * t426;
      double t472 = t32 * t89;
      double t487 =
        -0.2764234678236024E3 * t376 * t372 * t368 * t426 -
        0.2764234678236024E3 * t377 * t426 + 0.7181663018434603E3 * t163 * t325 -
        0.7181663018434603E3 * t163 * t421 * t5 +
        0.7181663018434603E3 * t197 * t325 + 0.7181663018434603E3 * t201 * t325 +
        0.7181663018434603E3 * t211 * t325 + 0.7181663018434603E3 * t207 * t325 -
        0.6902643889536271E4 * t141 * t34 * t89 * t353 +
        0.1380528777907254E5 * t142 * t353 + 0.1380528777907254E5 * t142 * t350 +
        0.6902643889536271E4 * t142 * t347 - 0.6902643889536271E4 * t142 * t355 -
        0.6902643889536271E4 * t141 * t472 * t353 -
        0.6902643889536271E4 * t248 * t472 * t353 -
        0.6902643889536271E4 * t36 * t472 * t353 - 0.1160277345846308E3 * t30 +
        0.1160277345846308E3 * t5 - 0.108485866736002E4 * t89 * t352 * t5;
      double t491 = t14 * t100;
      double t508 = t41 * t6;
      double t531 = t92 * t90 * t66;
      double t542 =
        0.1473126813211175E4 * t92 * t316 * t5 +
        0.3430624328232629E3 * t491 * t30 - 0.3430624328232629E3 * t491 * t5 -
        0.6861248656465258E3 * t14 * t349 * t5 -
        0.6861248656465258E3 * t14 * t352 * t5 - 0.547272071631899E3 * t19 * t30 -
        0.1887145205123939E3 * t44 * t42 * t6 * t5 -
        0.1887145205123939E3 * t44 * t508 * t5 -
        0.1887145205123939E3 * t45 * t30 + 0.1887145205123939E3 * t45 * t5 -
        0.1887145205123939E3 * t42 * t508 * t5 + 0.547272071631899E3 * t19 * t5 +
        0.1094544143263798E4 * t18 * t343 * t5 +
        0.1094544143263798E4 * t18 * t335 * t5 +
        0.1094544143263798E4 * t10 * t335 * t5 +
        0.2329218006006368E3 * t531 * t30 - 0.2329218006006368E3 * t531 * t5 +
        0.2329218006006368E3 * t90 * t324 * t5 +
        0.2329218006006368E3 * t92 * t324 * t5;
      double t545 = t96 + t169 + t243 + t301 + t357 + t439 + t487 + t542;
      return t545;
    }

    double
    ortho2_f53y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * y;
      double t7 = -t3 - t6 + 0.2650553239294647;
      double t9 = -t3 - t6 - 0.2147684835193549;
      double t10 = -t3 - t6 - 0.7852315164806451;
      double t11 = t10 * t9;
      double t12 = -t3 - t6 - 0.1265055323929465E1;
      double t13 = t12 * t11;
      double t16 = -t3 - t1;
      double t17 = t16 * t2;
      double t18 = t4 * t17;
      double t19 = 0.1E1 * x;
      double t20 = 0.1154653670707977E1 + t19 + t1;
      double t21 = 0.5 + t19 + t1;
      double t22 = t21 * t20;
      double t23 = -0.1546536707079771 + t19 + t1;
      double t24 = t23 * t22;
      double t30 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t35 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t36 = t35 * t5;
      double t38 = t4 * t16;
      double t39 = t35 * t38;
      double t43 = -t3 - t6 - 0.5278640450004206E-1;
      double t45 = -t3 - t6 - 0.9472135954999579;
      double t46 = t45 * t43 * t35;
      double t53 = -t3 - t6 + 0.3302238962785669;
      double t54 = -t3 - t6 - 0.3115120652928579E-1;
      double t55 = t54 * t53;
      double t56 = -t3 - t6 - 0.9688487934707142;
      double t57 = -t3 - t6 - 0.1330223896278567E1;
      double t58 = t57 * t56;
      double t59 = t58 * t55;
      double t62 = 0.5278640450004206E-1 + t19 + t1;
      double t64 = t45 * t43 * t62;
      double t67 = 0.9472135954999579 + t19 + t1;
      double t68 = t67 * t38;
      double t69 = -t3 - t6 + 0.1546536707079771;
      double t71 = -t3 - t6 - 0.5;
      double t72 = -t3 - t6 - 0.1154653670707977E1;
      double t73 = t72 * t71;
      double t74 = t73 * t69 * t62;
      double t80 = t58 * t71 * t54;
      double t86 = t9 * t7;
      double t87 = t12 * t10;
      double t88 = t87 * t86;
      double t91 = 0.1265055323929465E1 + t19 + t1;
      double t93 = t91 * t4 * t17;
      double t94 = 0.2147684835193549 + t19 + t1;
      double t95 = -0.2650553239294647 + t19 + t1;
      double t96 = t95 * t94;
      double t97 = t45 * t43;
      double t98 = t97 * t96;
      double t101 = t20 * t5;
      double t102 = t23 * t21;
      double t103 = t30 * t102;
      double t106 = t67 * t4;
      double t107 = t106 * t17;
      double t108 = t7 * t62;
      double t112 =
        -0.1418271315264546E4 * t13 * t7 * t5 - 0.6350798073234854E4 * t24 * t18 -
        0.1099658061698571E3 * t18 + 0.4681841558216832E2 * t30 * t5 +
        0.5886258994100211E2 * t36 - 0.5886258994100211E2 * t39 -
        0.4681841558216832E2 * t30 * t38 - 0.2329218006006368E3 * t46 * t38 -
        0.2836542630529092E4 * t12 * t10 * t7 * t18 -
        0.3799026905676787E1 * t59 * t18 + 0.5998165083601091E4 * t64 * t18 +
        0.1163991078526991E4 * t74 * t68 +
        0.2271041252605316E4 * t80 * t53 * t4 * t17 -
        0.1004149343561548E4 * t30 * t22 * t18 +
        0.9603869825676393E2 * t88 * t36 + 0.4237517115141881E4 * t98 * t93 +
        0.1004149343561548E4 * t103 * t101 -
        0.1380528777907254E5 * t11 * t108 * t107;
      double t116 = t91 * t38;
      double t117 = 0.7852315164806451 + t19 + t1;
      double t118 = t94 * t117;
      double t119 = t95 * t118;
      double t122 = t20 * t38;
      double t125 = t71 * t69;
      double t126 = t72 * t125;
      double t129 = t62 * t67;
      double t130 = t129 * t5;
      double t133 = t20 * t4;
      double t134 = t133 * t17;
      double t138 = t72 * t69;
      double t146 = t7 * t35;
      double t147 = t12 * t9;
      double t151 = t69 * t35;
      double t155 = t30 * t95;
      double t156 = t155 * t118;
      double t159 = t23 * t20;
      double t163 = t62 * t4;
      double t164 = t163 * t17;
      double t169 = -t3 - t6 + 0.3717401485096066;
      double t171 = t169 * t4 * t17;
      double t172 = -t3 - t6 - 0.2907007820975211;
      double t173 = -t3 - t6 - 0.7092992179024789;
      double t175 = -t3 - t6 - 0.1091700181433142E1;
      double t176 = -t3 - t6 - 0.1371740148509607E1;
      double t177 = t176 * t175;
      double t178 = t177 * t173 * t172;
      double t181 = 0.1330223896278567E1 + t19 + t1;
      double t183 = 0.9688487934707142 + t19 + t1;
      double t185 = 0.3115120652928579E-1 + t19 + t1;
      double t186 = -0.3302238962785669 + t19 + t1;
      double t187 = t186 * t185;
      double t188 = t187 * t21 * t183;
      double t197 = t53 * t35;
      double t201 =
        -0.2836542630529092E4 * t12 * t86 * t18 +
        0.5650028391921024E3 * t119 * t116 - 0.1004149343561548E4 * t103 * t122 +
        0.1522103401239837E4 * t126 * t36 - 0.6902643889536271E4 * t88 * t130 -
        0.1805304996254271E5 * t73 * t102 * t134 -
        0.1805304996254271E5 * t138 * t102 * t134 +
        0.3044206802479674E4 * t72 * t71 * t35 * t18 +
        0.1920773965135279E3 * t147 * t146 * t18 +
        0.3044206802479674E4 * t72 * t151 * t18 -
        0.3493844575728145E3 * t156 * t116 +
        0.157830039561319E4 * t97 * t159 * t18 +
        0.6902643889536271E4 * t88 * t164 + 0.6902643889536271E4 * t88 * t107 -
        0.5528469356472047E3 * t178 * t171 -
        0.3219964681375463E3 * t188 * t181 * t5 -
        0.2836542630529092E4 * t10 * t86 * t18 +
        0.157830039561319E4 * t97 * t22 * t18 -
        0.7181663018434603E3 * t80 * t197 * t38;
      double t213 = -t3 - t6 + 0.917001814331423E-1;
      double t225 = t69 * t4;
      double t235 = t129 * t38;
      double t242 = t172 * t213;
      double t256 = t73 * t69 * t23;
      double t259 = t21 * t4;
      double t267 =
        0.1920773965135279E3 * t11 * t146 * t18 -
        0.1899513452838393E1 * t80 * t53 * t5 - 0.547272071631899E3 * t24 * t5 +
        0.1887145205123939E3 * t126 * t38 -
        0.5528469356472047E3 * t177 * t173 * t213 * t171 -
        0.3774290410247878E3 * t72 * t71 * t4 * t17 -
        0.1805304996254271E5 * t125 * t102 * t134 -
        0.3774290410247878E3 * t72 * t225 * t17 -
        0.1887145205123939E3 * t126 * t5 +
        0.4237517115141881E4 * t97 * t95 * t117 * t93 +
        0.6902643889536271E4 * t88 * t235 +
        0.5998165083601091E4 * t45 * t43 * t67 * t18 -
        0.5528469356472047E3 * t177 * t242 * t171 +
        0.3219964681375463E3 * t188 * t18 -
        0.1199633016720218E5 * t45 * t129 * t18 +
        0.3219964681375463E3 * t188 * t181 * t38 +
        0.9026524981271357E4 * t256 * t22 * t38 +
        0.547272071631899E3 * t23 * t259 * t17 +
        0.1163991078526991E4 * t73 * t69 * t67 * t18;
      double t279 = t30 * t35;
      double t288 = t43 * t4;
      double t296 = t35 * t4;
      double t297 = t296 * t17;
      double t304 = t183 * t181;
      double t320 =
        -0.2327982157053983E4 * t73 * t129 * t18 +
        0.3219964681375463E3 * t187 * t21 * t181 * t18 -
        0.5528469356472047E3 * t176 * t173 * t242 * t171 +
        0.1491083359667674E3 * t279 * t38 - 0.7381565053942821E3 * t130 +
        0.7381565053942821E3 * t107 + 0.7381565053942821E3 * t164 -
        0.1946100011820659E4 * t45 * t4 * t17 -
        0.1946100011820659E4 * t288 * t17 + 0.9730500059103294E3 * t97 * t38 -
        0.4715219597725897E3 * t30 * t4 * t17 + 0.9430439195451794E3 * t297 -
        0.1491083359667674E3 * t279 * t5 -
        0.2327982157053983E4 * t138 * t129 * t18 +
        0.3219964681375463E3 * t187 * t304 * t18 -
        0.1199633016720218E5 * t43 * t129 * t18 -
        0.2327982157053983E4 * t125 * t129 * t18 +
        0.547272071631899E3 * t23 * t133 * t17 +
        0.4237517115141881E4 * t97 * t118 * t93;
      double t343 = t58 * t71 * t53;
      double t372 = t56 * t71 * t55;
      double t376 = t57 * t71 * t55;
      double t381 =
        -0.8475034230283761E4 * t45 * t95 * t118 * t93 +
        0.547272071631899E3 * t21 * t133 * t17 -
        0.8475034230283761E4 * t43 * t95 * t118 * t93 -
        0.1004149343561548E4 * t103 * t18 +
        0.3219964681375463E3 * t186 * t21 * t304 * t18 +
        0.1436332603686921E4 * t80 * t297 + 0.1436332603686921E4 * t343 * t297 +
        0.3219964681375463E3 * t185 * t21 * t304 * t18 -
        0.1160277345846308E3 * t5 - 0.3156600791226379E4 * t45 * t23 * t22 * t18 +
        0.1436332603686921E4 * t59 * t297 + 0.7381565053942821E3 * t235 -
        0.3774290410247878E3 * t71 * t225 * t17 -
        0.3156600791226379E4 * t43 * t23 * t22 * t18 +
        0.1418271315264546E4 * t13 * t38 * t7 -
        0.1004149343561548E4 * t30 * t159 * t18 +
        0.1436332603686921E4 * t372 * t297 + 0.1436332603686921E4 * t376 * t297 -
        0.3799026905676787E1 * t80 * t18;
      double t385 = t30 * t129;
      double t388 = t117 * t91;
      double t411 = t213 * t169;
      double t430 = t91 * t5;
      double t435 =
        0.1899513452838393E1 * t80 * t53 * t38 +
        0.3430624328232629E3 * t385 * t5 -
        0.2209701329985117E4 * t96 * t388 * t18 -
        0.1522103401239837E4 * t126 * t39 -
        0.9026524981271357E4 * t256 * t22 * t5 + 0.547272071631899E3 * t24 * t38 -
        0.9603869825676393E2 * t88 * t39 +
        0.1920773965135279E3 * t87 * t9 * t35 * t18 +
        0.2329218006006368E3 * t46 * t5 - 0.3799026905676787E1 * t18 * t343 -
        0.3430624328232629E3 * t385 * t38 -
        0.2764234678236024E3 * t178 * t411 * t5 -
        0.2169717334720039E4 * t62 * t106 * t17 +
        0.7365634066055877E3 * t45 * t288 * t17 +
        0.4658436012012737E3 * t45 * t296 * t17 -
        0.5528469356472047E3 * t175 * t173 * t242 * t171 +
        0.4813313582207044E4 * t126 * t18 - 0.5650028391921024E3 * t119 * t430 +
        0.5650028391921024E3 * t119 * t18;
      double t476 = t94 * t91;
      double t493 =
        0.9026524981271357E4 * t256 * t259 * t17 -
        0.3799026905676787E1 * t376 * t18 -
        0.1380528777907254E5 * t87 * t9 * t62 * t107 +
        0.5998165083601091E4 * t64 * t68 -
        0.1380528777907254E5 * t87 * t108 * t107 -
        0.4237517115141881E4 * t98 * t388 * t5 +
        0.3044206802479674E4 * t71 * t151 * t18 +
        0.9026524981271357E4 * t256 * t134 -
        0.3430624328232629E3 * t30 * t163 * t17 -
        0.5528469356472047E3 * t178 * t213 * t4 * t17 +
        0.4237517115141881E4 * t98 * t117 * t4 * t17 +
        0.1920773965135279E3 * t87 * t146 * t18 +
        0.7181663018434603E3 * t80 * t197 * t5 +
        0.5650028391921024E3 * t95 * t476 * t18 +
        0.3037010300090164E3 * t88 * t18 - 0.9730500059103294E3 * t97 * t5 +
        0.4658436012012737E3 * t43 * t296 * t17 +
        0.4237517115141881E4 * t98 * t388 * t38 -
        0.3430624328232629E3 * t30 * t106 * t17;
      double t516 = t67 * t5;
      double t526 = t97 * t102;
      double t543 =
        0.3493844575728145E3 * t156 * t430 +
        0.5650028391921024E3 * t94 * t388 * t18 +
        0.5650028391921024E3 * t95 * t388 * t18 +
        0.9026524981271357E4 * t73 * t69 * t21 * t134 -
        0.3493844575728145E3 * t30 * t94 * t388 * t18 -
        0.3493844575728145E3 * t155 * t388 * t18 -
        0.3493844575728145E3 * t155 * t476 * t18 -
        0.1163991078526991E4 * t74 * t516 - 0.3493844575728145E3 * t156 * t18 +
        0.2764234678236024E3 * t178 * t411 * t38 -
        0.3799026905676787E1 * t372 * t18 + 0.157830039561319E4 * t526 * t122 -
        0.1380528777907254E5 * t147 * t108 * t107 -
        0.5998165083601091E4 * t64 * t516 + 0.1163991078526991E4 * t74 * t18 +
        0.157830039561319E4 * t526 * t18 - 0.157830039561319E4 * t526 * t101 +
        0.1160277345846308E3 * t38 - 0.2836542630529092E4 * t13 * t18;
      double t546 = t112 + t201 + t267 + t320 + t381 + t435 + t493 + t543;

      return t546;
    }

    // * f54 *********************************************************************

    double
    ortho2_f54 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t13 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t17 = 0.1E1 * y;
      double t18 = -t3 - t17 - 0.5278640450004206E-1;
      double t20 = -t3 - t17 - 0.9472135954999579;
      double t24 = 0.1E1 * x;
      double t25 = 0.9472135954999579 + t24 + t1;
      double t26 = t25 * t6;
      double t27 = 0.5278640450004206E-1 + t24 + t1;
      double t31 = t6 * t5;
      double t33 = -t3 - t17 + 0.1546536707079771;
      double t34 = -t3 - t17 - 0.5;
      double t36 = -t3 - t17 - 0.1154653670707977E1;
      double t37 = t36 * t34 * t33;
      double t40 = 0.1154653670707977E1 + t24 + t1;
      double t41 = 0.5 + t24 + t1;
      double t42 = t41 * t40;
      double t43 = -0.1546536707079771 + t24 + t1;
      double t47 = t27 * t25;
      double t55 = -t3 - t17 + 0.2650553239294647;
      double t56 = -t3 - t17 - 0.2147684835193549;
      double t58 = -t3 - t17 - 0.7852315164806451;
      double t59 = -t3 - t17 - 0.1265055323929465E1;
      double t61 = t59 * t58 * t56 * t55;
      double t68 = 0.1265055323929465E1 + t24 + t1;
      double t69 = 0.7852315164806451 + t24 + t1;
      double t71 = 0.2147684835193549 + t24 + t1;
      double t72 = -0.2650553239294647 + t24 + t1;
      double t78 = t36 * t34;
      double t82 = t20 * t18;
      double t86 =
        0.3348659386389576E3 * t13 * t10 * t5 +
        0.1720339219500472E4 * t20 * t18 * t6 * t5 +
        0.2093078722891669E4 * t27 * t26 * t5 + 0.2675916433772571E3 * t31 -
        0.7757360586169351E3 * t37 * t31 -
        0.5317640379405782E3 * t43 * t42 * t31 +
        0.6697229743186177E3 * t13 * t47 * t31 +
        0.6929602906904184E3 * t20 * t18 * t9 * t31 +
        0.1436624516341072E4 * t61 * t31 -
        0.3286703354692339E4 * t13 * t43 * t42 * t31 +
        0.3097119891778271E4 * t72 * t71 * t69 * t68 * t31 -
        0.218046329738323E4 * t78 * t33 * t9 * t31 +
        0.1276149101998645E5 * t82 * t47 * t31;
      double t103 = (0.1330223896278567E1 + t24 + t1) * t6 * t5;
      double t105 = t41 * (0.9688487934707142 + t24 + t1);
      double t108 =
        (-0.3302238962785669 + t24 + t1) * (0.3115120652928579E-1 + t24 + t1);
      double t113 = t26 * t5;
      double t121 = t40 * t6 * t5;
      double t122 = t43 * t41;
      double t126 = t10 * t5;
      double t127 = -t3 - t17 + 0.3302238962785669;
      double t128 = -t3 - t17 - 0.3115120652928579E-1;
      double t130 = -t3 - t17 - 0.9688487934707142;
      double t132 = -t3 - t17 - 0.1330223896278567E1;
      double t138 = t68 * t6 * t5;
      double t139 = t71 * t69;
      double MapleGenVar1 =
        0.162287634880332E3 * (-t3 - t17 - 0.1371740148509607E1) * (-t3 - t17 -
                    0.1091700181433142E1)
        * (-t3 - t17 - 0.7092992179024789) * (-t3 - t17 -
                0.2907007820975211) * (-t3 - t17 +
                     0.917001814331423E-1)
        * (-t3 - t17 + 0.3717401485096066) * t6 * t5 -
        0.1963591413892188E4 * t13 * t108 * t105 * t103 +
        0.7704596907161649E4 * t59 * t58 * t56 * t55 * t27 * t113 +
        0.1636896355275561E5 * t37 * t122 * t121 -
        0.5530210100784448E3 * t132 * t130 * t34 * t128 * t127 * t126 +
        0.1645844409636311E5 * t20 * t18 * t72 * t139 * t138 +
        0.2652950810378689E3 * t108 * t105 * t103;
      double t172 =
        MapleGenVar1 - 0.3342890675852115E4 * t78 * t33 * t27 * t113 +
        0.4423616526402223E3 * t61 * t126 +
        0.4374147456031232E3 * t13 * t72 * t139 * t138 +
        0.1183431720422361E3 * t13 * t6 * t5 -
        0.3596184470413046E4 * t82 * t122 * t121 + 0.1040167628112854E3 * t126 -
        0.2106844111412697E3 * t132 * t130 * t34 * t128 * t127 * t6 * t5;
      double t173 = t86 + t172;
      return t173;
    }

    double
    ortho2_f54x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.1154653670707977E1 + t8 + t1;
      double t10 = 0.5 + t8 + t1;
      double t11 = t10 * t9;
      double t14 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t18 = -0.1546536707079771 + t8 + t1;
      double t19 = t18 * t11;
      double t22 = 0.1E1 * y;
      double t23 = -t3 - t22 + 0.1546536707079771;
      double t24 = -t3 - t22 - 0.5;
      double t25 = t24 * t23;
      double t26 = -t3 - t22 - 0.1154653670707977E1;
      double t27 = t26 * t25;
      double t30 = t18 * t9;
      double t34 = t18 * t10;
      double t35 = t14 * t34;
      double t38 = t6 * t2;
      double t39 = -t3 - t22 + 0.2650553239294647;
      double t41 = -t3 - t22 - 0.2147684835193549;
      double t42 = -t3 - t22 - 0.7852315164806451;
      double t43 = t42 * t41;
      double t44 = -t3 - t22 - 0.1265055323929465E1;
      double t45 = t44 * t43;
      double t51 = t41 * t39;
      double t64 = t9 * t38;
      double t67 = t9 * t5;
      double t70 = 0.1265055323929465E1 + t8 + t1;
      double t71 = 0.7852315164806451 + t8 + t1;
      double t72 = t71 * t70;
      double t73 = 0.2147684835193549 + t8 + t1;
      double t77 = -0.2650553239294647 + t8 + t1;
      double t81 = t73 * t70;
      double t85 = t73 * t71;
      double t86 = t77 * t85;
      double t89 = t70 * t5;
      double t92 = t70 * t38;
      double t95 =
        -0.3286703354692339E4 * t14 * t11 * t7 - 0.5196734297072026E4 * t19 * t7 +
        0.6895230374132071E4 * t27 * t7 - 0.3286703354692339E4 * t14 * t30 * t7 -
        0.3286703354692339E4 * t35 * t7 - 0.7183122581705362E3 * t45 * t39 * t38 +
        0.7183122581705362E3 * t45 * t39 * t5 -
        0.7183122581705362E3 * t42 * t51 * t7 -
        0.7183122581705362E3 * t44 * t51 * t7 -
        0.7183122581705362E3 * t44 * t42 * t39 * t7 -
        0.7183122581705362E3 * t45 * t7 + 0.1643351677346169E4 * t35 * t64 -
        0.1643351677346169E4 * t35 * t67 + 0.3097119891778271E4 * t73 * t72 * t7 +
        0.3097119891778271E4 * t77 * t72 * t7 +
        0.3097119891778271E4 * t77 * t81 * t7 + 0.3097119891778271E4 * t86 * t7 +
        0.1548559945889136E4 * t86 * t89 - 0.1548559945889136E4 * t86 * t92;
      double t96 = 0.9472135954999579 + t8 + t1;
      double t97 = 0.5278640450004206E-1 + t8 + t1;
      double t98 = t97 * t96;
      double t99 = -t3 - t22 - 0.5278640450004206E-1;
      double t103 = -t3 - t22 - 0.9472135954999579;
      double t112 = t103 * t99 * t97;
      double t115 = t96 * t5;
      double t118 = t96 * t38;
      double t123 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t124 = t23 * t123;
      double t135 = t123 * t5;
      double t138 = t123 * t38;
      double t141 = -t3 - t22 + 0.3302238962785669;
      double t143 = -t3 - t22 - 0.9688487934707142;
      double t144 = -t3 - t22 - 0.1330223896278567E1;
      double t145 = t144 * t143;
      double t146 = t145 * t24 * t141;
      double t149 = -t3 - t22 - 0.3115120652928579E-1;
      double t151 = t145 * t24 * t149;
      double t160 = t77 * t73;
      double t172 = t103 * t99;
      double MapleGenVar1 =
        -0.6380745509993226E4 * t99 * t98 * t7 -
        0.6380745509993226E4 * t103 * t98 * t7 +
        0.1276149101998645E5 * t103 * t99 * t96 * t7 +
        0.1276149101998645E5 * t112 * t7 + 0.6380745509993226E4 * t112 * t115 -
        0.6380745509993226E4 * t112 * t118 +
        0.1090231648691615E4 * t24 * t124 * t7 +
        0.1090231648691615E4 * t26 * t124 * t7 +
        0.1090231648691615E4 * t26 * t24 * t123 * t7 -
        0.1090231648691615E4 * t27 * t135;
      double t179 =
        MapleGenVar1 + 0.1090231648691615E4 * t27 * t138 +
        0.1053422055706348E3 * t146 * t7 + 0.1053422055706348E3 * t151 * t7 -
        0.1053422055706348E3 * t151 * t141 * t5 +
        0.1053422055706348E3 * t151 * t141 * t38 +
        0.6916134391244957E3 * t160 * t72 * t7 +
        0.1798092235206523E4 * t99 * t18 * t11 * t7 +
        0.1798092235206523E4 * t103 * t18 * t11 * t7 -
        0.3596184470413046E4 * t172 * t11 * t7 -
        0.3596184470413046E4 * t172 * t30 * t7;
      double t181 = t172 * t34;
      double t188 = t44 * t42;
      double t189 = t188 * t51;
      double t192 = t149 * t141;
      double t193 = t145 * t192;
      double t201 = t143 * t24 * t192;
      double t205 = t144 * t24 * t192;
      double t208 = 0.1330223896278567E1 + t8 + t1;
      double t209 = 0.9688487934707142 + t8 + t1;
      double t210 = t209 * t208;
      double t211 = 0.3115120652928579E-1 + t8 + t1;
      double t212 = t211 * t10;
      double t216 = -0.3302238962785669 + t8 + t1;
      double t221 = t216 * t211;
      double t229 = t10 * t209;
      double t230 = t221 * t229;
      double t243 = t39 * t123;
      double t247 = t44 * t41;
      double t254 =
        -0.3596184470413046E4 * t181 * t7 - 0.1798092235206523E4 * t181 * t67 +
        0.1798092235206523E4 * t181 * t64 - 0.139887037185934E4 * t189 * t7 +
        0.1053422055706348E3 * t193 * t7 + 0.2211808263201112E3 * t189 * t135 -
        0.2211808263201112E3 * t189 * t138 + 0.1053422055706348E3 * t201 * t7 +
        0.1053422055706348E3 * t205 * t7 +
        0.2652950810378689E3 * t212 * t210 * t7 +
        0.2652950810378689E3 * t216 * t10 * t210 * t7 +
        0.2652950810378689E3 * t221 * t210 * t7 +
        0.2652950810378689E3 * t221 * t10 * t208 * t7 +
        0.2652950810378689E3 * t230 * t7 +
        0.1326475405189345E3 * t230 * t208 * t5 -
        0.1326475405189345E3 * t230 * t208 * t38 -
        0.2211808263201112E3 * t188 * t41 * t123 * t7 -
        0.2211808263201112E3 * t43 * t243 * t7 -
        0.2211808263201112E3 * t247 * t243 * t7 -
        0.2211808263201112E3 * t188 * t243 * t7;
      double t256 = t26 * t24;
      double t257 = t256 * t23 * t97;
      double t262 = t26 * t23;
      double t279 = t14 * t77;
      double t286 = t279 * t85;
      double t305 = t14 * t216;
      double t306 = t305 * t212;
      double t312 =
        -0.1671445337926057E4 * t257 * t115 + 0.1671445337926057E4 * t257 * t118 +
        0.1671445337926057E4 * t262 * t98 * t7 +
        0.1671445337926057E4 * t256 * t98 * t7 -
        0.3342890675852115E4 * t256 * t23 * t96 * t7 -
        0.3342890675852115E4 * t257 * t7 +
        0.4374147456031232E3 * t14 * t73 * t72 * t7 +
        0.4374147456031232E3 * t279 * t72 * t7 +
        0.4374147456031232E3 * t279 * t81 * t7 +
        0.4374147456031232E3 * t286 * t7 + 0.2187073728015616E3 * t286 * t89 -
        0.2187073728015616E3 * t286 * t92 +
        0.1671445337926057E4 * t25 * t98 * t7 - 0.520083814056427E2 * t138 -
        0.1418129007298477E3 * t7 + 0.520083814056427E2 * t135 +
        0.5917158602111804E2 * t14 * t5 - 0.5917158602111804E2 * t14 * t38 -
        0.1963591413892188E4 * t306 * t209 * t6 * t5 -
        0.981795706946094E3 * t306 * t210 * t5;
      double t319 = t208 * t6 * t5;
      double t335 = t96 * t6;
      double t339 = t99 * t6;
      double t345 = t256 * t23 * t18;
      double t351 = t10 * t6;
      double t355 = t9 * t6;
      double t356 = t355 * t5;
      double t363 = -t3 - t22 + 0.3717401485096066;
      double t364 = -t3 - t22 + 0.917001814331423E-1;
      double t365 = t364 * t363;
      double t367 = -t3 - t22 - 0.2907007820975211;
      double t368 = -t3 - t22 - 0.7092992179024789;
      double t370 = -t3 - t22 - 0.1091700181433142E1;
      double t371 = -t3 - t22 - 0.1371740148509607E1;
      double t372 = t371 * t370;
      double t373 = t372 * t368 * t367;
      double t394 = t70 * t6 * t5;
      double t395 = t172 * t160;
      MapleGenVar1 =
        0.981795706946094E3 * t306 * t210 * t38 -
        0.1963591413892188E4 * t305 * t229 * t319 -
        0.1963591413892188E4 * t305 * t211 * t209 * t319 -
        0.1963591413892188E4 * t306 * t319 -
        0.1963591413892188E4 * t14 * t211 * t229 * t319 -
        0.3104710630924854E4 * t230 * t319 +
        0.1058925000094643E4 * t97 * t335 * t5 -
        0.2191332846634096E4 * t103 * t339 * t5 +
        0.8184481776377804E4 * t345 * t11 * t5 -
        0.8184481776377804E4 * t345 * t11 * t38;
      double t398 =
        MapleGenVar1 + 0.1636896355275561E5 * t345 * t351 * t5 +
        0.1636896355275561E5 * t256 * t23 * t10 * t356 +
        0.1636896355275561E5 * t345 * t356 -
        0.8114381744016598E2 * t373 * t365 * t38 -
        0.8184481776377804E4 * t256 * t34 * t356 +
        0.1748805985774818E4 * t151 * t141 * t6 * t5 -
        0.8114381744016598E2 * t373 * t364 * t6 * t5 +
        0.8114381744016598E2 * t373 * t365 * t5 -
        0.8184481776377804E4 * t262 * t34 * t356 +
        0.1645844409636311E5 * t395 * t394;
      double t428 = t363 * t6 * t5;
      double t435 = t367 * t364;
      double t447 = t14 * t98;
      double t452 = t141 * t123;
      double t456 = t123 * t6;
      double t457 = t456 * t5;
      MapleGenVar1 =
        0.1645844409636311E5 * t395 * t71 * t6 * t5 +
        0.8229222048181555E4 * t395 * t72 * t5 -
        0.8229222048181555E4 * t395 * t72 * t38 -
        0.8184481776377804E4 * t25 * t34 * t356 -
        0.8229222048181555E4 * t99 * t77 * t85 * t394 -
        0.8229222048181555E4 * t103 * t77 * t85 * t394 +
        0.1645844409636311E5 * t172 * t85 * t394 +
        0.1645844409636311E5 * t172 * t77 * t71 * t394 -
        0.8114381744016598E2 * t372 * t368 * t364 * t428 -
        0.8114381744016598E2 * t373 * t428;
      double t467 =
        MapleGenVar1 - 0.8114381744016598E2 * t371 * t368 * t435 * t428 -
        0.8114381744016598E2 * t372 * t435 * t428 -
        0.8114381744016598E2 * t370 * t368 * t435 * t428 -
        0.3348614871593089E3 * t447 * t38 + 0.3348614871593089E3 * t447 * t5 +
        0.2765105050392224E3 * t151 * t452 * t38 +
        0.2765105050392224E3 * t151 * t457 -
        0.2765105050392224E3 * t151 * t452 * t5 +
        0.2765105050392224E3 * t193 * t457 + 0.2765105050392224E3 * t146 * t457;
      double t471 = t97 * t6;
      double t480 = t98 * t38;
      double t483 = t471 * t5;
      double t486 = t98 * t5;
      double t489 = t335 * t5;
      double t498 = t39 * t97;
      double t508 = t23 * t6;
      double t525 = t14 * t123;
      MapleGenVar1 =
        0.2765105050392224E3 * t205 * t457 +
        0.6697229743186177E3 * t14 * t471 * t5 +
        0.2765105050392224E3 * t201 * t457 +
        0.6697229743186177E3 * t14 * t335 * t5 -
        0.3852298453580824E4 * t189 * t480 + 0.7704596907161649E4 * t189 * t483 +
        0.3852298453580824E4 * t189 * t486 -
        0.3852298453580824E4 * t188 * t41 * t97 * t489 +
        0.7704596907161649E4 * t189 * t489 + 0.2658820189702891E3 * t19 * t38;
      double t528 =
        MapleGenVar1 - 0.3852298453580824E4 * t188 * t498 * t489 -
        0.3852298453580824E4 * t43 * t498 * t489 -
        0.3852298453580824E4 * t247 * t498 * t489 +
        0.3878680293084675E3 * t26 * t508 * t5 +
        0.3878680293084675E3 * t26 * t24 * t6 * t5 +
        0.3878680293084675E3 * t27 * t38 - 0.3878680293084675E3 * t27 * t5 +
        0.3878680293084675E3 * t24 * t508 * t5 +
        0.8601696097502358E3 * t172 * t5 + 0.1674329693194788E3 * t525 * t5;
      double t556 = t103 * t99 * t123;
      double t571 =
        -0.1674329693194788E3 * t525 * t38 + 0.5294695384546455E3 * t457 -
        0.1058939076909291E4 * t14 * t6 * t5 - 0.8601696097502358E3 * t339 * t5 -
        0.8601696097502358E3 * t103 * t6 * t5 - 0.1046539361445834E4 * t480 +
        0.2093078722891669E4 * t489 + 0.1046539361445834E4 * t486 +
        0.2093078722891669E4 * t483 - 0.2658820189702891E3 * t19 * t5 -
        0.5317640379405782E3 * t10 * t355 * t5 -
        0.5317640379405782E3 * t18 * t355 * t5 -
        0.5317640379405782E3 * t18 * t351 * t5 +
        0.3464801453452092E3 * t556 * t5 - 0.3464801453452092E3 * t556 * t38 -
        0.3464801453452092E3 * t103 * t456 * t5 -
        0.3464801453452092E3 * t99 * t456 * t5 -
        0.8601696097502358E3 * t172 * t38 + 0.1337958216886285E3 * t5 -
        0.1337958216886285E3 * t38;
      double t574 = t95 + t179 + t254 + t312 + t398 + t467 + t528 + t571;
      return t574;
    }

    double
    ortho2_f54y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t8 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t9 = t8 * t5;
      double t11 = -t3 - t1;
      double t12 = t11 * t2;
      double t13 = t4 * t12;
      double t15 = t4 * t11;
      double t16 = t8 * t15;
      double t20 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t23 = 0.1E1 * x;
      double t24 = 0.1154653670707977E1 + t23 + t1;
      double t25 = 0.5 + t23 + t1;
      double t26 = t25 * t24;
      double t30 = 0.1330223896278567E1 + t23 + t1;
      double t32 = 0.9688487934707142 + t23 + t1;
      double t33 = t25 * t32;
      double t34 = 0.3115120652928579E-1 + t23 + t1;
      double t35 = -0.3302238962785669 + t23 + t1;
      double t36 = t35 * t34;
      double t37 = t36 * t33;
      double t40 = -0.1546536707079771 + t23 + t1;
      double t41 = t40 * t26;
      double t44 = t24 * t4;
      double t45 = t44 * t12;
      double t46 = t40 * t25;
      double t47 = 0.1E1 * y;
      double t48 = -t3 - t47 + 0.1546536707079771;
      double t49 = -t3 - t47 - 0.1154653670707977E1;
      double t50 = t49 * t48;
      double t56 = -t3 - t47 - 0.5;
      double t61 = 0.1265055323929465E1 + t23 + t1;
      double t62 = 0.7852315164806451 + t23 + t1;
      double t63 = t62 * t61;
      double t64 = 0.2147684835193549 + t23 + t1;
      double t68 = t61 * t5;
      double t69 = t64 * t62;
      double t70 = -0.2650553239294647 + t23 + t1;
      double t71 = t20 * t70;
      double t72 = t71 * t69;
      double t76 = t61 * t4 * t12;
      double t78 = -t3 - t47 - 0.5278640450004206E-1;
      double t79 = -t3 - t47 - 0.9472135954999579;
      double t80 = t79 * t78;
      double t89 = t25 * t34;
      double t90 = t20 * t35;
      double t91 = t90 * t89;
      double t94 = t61 * t15;
      double t97 = 0.9472135954999579 + t23 + t1;
      double t98 = 0.5278640450004206E-1 + t23 + t1;
      double t99 = t98 * t97;
      double t100 = t49 * t56;
      double t109 = t30 * t4 * t12;
      double t112 = -0.520083814056427E2 * t9 + 0.2097690265220458E3 * t13 +
        0.520083814056427E2 * t16 - 0.5917158602111804E2 * t20 * t5 -
        0.1643351677346169E4 * t20 * t26 * t13 -
        0.1326475405189345E3 * t37 * t30 * t5 + 0.2658820189702891E3 * t41 * t5 -
        0.1636896355275561E5 * t50 * t46 * t45 +
        0.1326475405189345E3 * t37 * t13 +
        0.218046329738323E4 * t49 * t56 * t8 * t13 +
        0.1548559945889136E4 * t64 * t63 * t13 -
        0.2187073728015616E3 * t72 * t68 +
        0.8229222048181555E4 * t80 * t70 * t62 * t76 +
        0.8229222048181555E4 * t80 * t69 * t76 -
        0.981795706946094E3 * t91 * t32 * t4 * t12 +
        0.2187073728015616E3 * t72 * t94 +
        0.3342890675852115E4 * t100 * t99 * t13 +
        0.1326475405189345E3 * t36 * t25 * t30 * t13 -
        0.981795706946094E3 * t91 * t109;
      double t116 = t56 * t48;
      double t124 = t20 * t46;
      double t127 = -t3 - t47 + 0.2650553239294647;
      double t128 = -t3 - t47 - 0.7852315164806451;
      double t130 = -t3 - t47 - 0.1265055323929465E1;
      double t135 = t100 * t48 * t98;
      double t138 = t32 * t30;
      double t143 = t79 * t78 * t8;
      double t150 = t48 * t8;
      double t154 = -t3 - t47 + 0.3717401485096066;
      double t156 = t154 * t4 * t12;
      double t157 = -t3 - t47 - 0.2907007820975211;
      double t158 = -t3 - t47 - 0.7092992179024789;
      double t160 = -t3 - t47 - 0.1091700181433142E1;
      double t161 = -t3 - t47 - 0.1371740148509607E1;
      double t162 = t161 * t160;
      double t163 = t162 * t158 * t157;
      double t169 = -t3 - t47 - 0.2147684835193549;
      double t170 = t169 * t127;
      double t178 = t97 * t5;
      double t180 = t79 * t78 * t98;
      double t187 = t8 * t4;
      double t188 = t187 * t12;
      double t189 = -t3 - t47 + 0.3302238962785669;
      double t190 = -t3 - t47 - 0.3115120652928579E-1;
      double t191 = t190 * t189;
      double t192 = -t3 - t47 - 0.1330223896278567E1;
      double t194 = t192 * t56 * t191;
      double t197 = t40 * t24;
      double t204 = -t3 - t47 - 0.9688487934707142;
      double t206 = t204 * t56 * t191;
      double t209 =
        0.3342890675852115E4 * t50 * t99 * t13 +
        0.3342890675852115E4 * t116 * t99 * t13 -
        0.981795706946094E3 * t90 * t34 * t32 * t109 -
        0.1643351677346169E4 * t124 * t13 -
        0.1436624516341072E4 * t130 * t128 * t127 * t13 -
        0.1671445337926057E4 * t135 * t13 +
        0.1326475405189345E3 * t36 * t138 * t13 -
        0.3464801453452092E3 * t143 * t5 +
        0.1326475405189345E3 * t35 * t25 * t138 * t13 +
        0.218046329738323E4 * t49 * t150 * t13 -
        0.162287634880332E3 * t163 * t156 -
        0.1276149101998645E5 * t78 * t99 * t13 -
        0.1436624516341072E4 * t130 * t170 * t13 +
        0.3596184470413046E4 * t79 * t40 * t26 * t13 -
        0.6380745509993226E4 * t180 * t178 -
        0.1645844409636311E5 * t79 * t70 * t69 * t76 +
        0.5530210100784448E3 * t194 * t188 -
        0.1643351677346169E4 * t20 * t197 * t13 -
        0.1636896355275561E5 * t116 * t46 * t45 +
        0.5530210100784448E3 * t206 * t188;
      double t214 = t99 * t15;
      double t215 = t130 * t128;
      double t216 = t215 * t170;
      double t221 = t100 * t48 * t40;
      double t228 = t204 * t192;
      double t229 = t228 * t56 * t190;
      double t232 = t49 * t116;
      double t235 = t97 * t15;
      double t250 = t64 * t61;
      double t254 = t25 * t4;
      double t261 = t127 * t8;
      double t276 = t130 * t169;
      double t281 = t228 * t56 * t189;
      double MapleGenVar1 =
        -0.6929602906904184E3 * t79 * t187 * t12 +
        0.3852298453580824E4 * t216 * t214 +
        0.8184481776377804E4 * t221 * t26 * t15 +
        0.2187073728015616E3 * t72 * t13 -
        0.1053422055706348E3 * t229 * t189 * t15 -
        0.3878680293084675E3 * t232 * t15 - 0.1671445337926057E4 * t135 * t235 -
        0.1645844409636311E5 * t78 * t70 * t69 * t76 -
        0.4423616526402223E3 * t215 * t169 * t8 * t13 +
        0.3596184470413046E4 * t78 * t40 * t26 * t13;
      double t284 =
        MapleGenVar1 + 0.2187073728015616E3 * t71 * t250 * t13 -
        0.2658820189702891E3 * t40 * t254 * t12 +
        0.2187073728015616E3 * t71 * t63 * t13 -
        0.4423616526402223E3 * t215 * t261 * t13 +
        0.1326475405189345E3 * t89 * t138 * t13 +
        0.6380745509993226E4 * t180 * t13 +
        0.6380745509993226E4 * t79 * t78 * t97 * t13 +
        0.2106844111412697E3 * t229 * t13 -
        0.4423616526402223E3 * t276 * t261 * t13 +
        0.2106844111412697E3 * t281 * t13;
      double t291 = -t3 - t47 + 0.917001814331423E-1;
      double t299 = t24 * t5;
      double t302 = t20 * t99;
      double t316 = t99 * t5;
      double t337 = t70 * t69;
      double t340 = t48 * t4;
      double t346 =
        -0.8184481776377804E4 * t221 * t26 * t5 -
        0.981795706946094E3 * t90 * t33 * t109 -
        0.162287634880332E3 * t162 * t158 * t291 * t156 -
        0.6929602906904184E3 * t78 * t187 * t12 +
        0.1643351677346169E4 * t124 * t299 - 0.3348614871593089E3 * t302 * t5 -
        0.1436624516341072E4 * t128 * t170 * t13 +
        0.7757360586169351E3 * t49 * t56 * t4 * t12 +
        0.2187073728015616E3 * t20 * t64 * t63 * t13 -
        0.3852298453580824E4 * t216 * t316 + 0.3447615187066035E4 * t232 * t13 +
        0.2211808263201112E3 * t216 * t16 -
        0.981795706946094E3 * t20 * t34 * t33 * t109 -
        0.2658820189702891E3 * t40 * t44 * t12 +
        0.8184481776377804E4 * t221 * t45 +
        0.8184481776377804E4 * t221 * t254 * t12 +
        0.2106844111412697E3 * t194 * t13 + 0.1548559945889136E4 * t337 * t13 +
        0.7757360586169351E3 * t49 * t340 * t12 +
        0.3878680293084675E3 * t232 * t5;
      double t352 = t157 * t291;
      double t361 = t98 * t4;
      double t370 = t361 * t12;
      double t373 = t291 * t154;
      double t377 = t128 * t169;
      double t392 = t70 * t64;
      double t393 = t80 * t392;
      double t396 = t97 * t4;
      double t397 = t396 * t12;
      double t405 = t189 * t8;
      double t411 =
        -0.6994351859296698E3 * t216 * t13 + 0.1337958216886285E3 * t15 -
        0.162287634880332E3 * t162 * t352 * t156 -
        0.2658820189702891E3 * t25 * t44 * t12 +
        0.3348614871593089E3 * t302 * t15 +
        0.3348614871593089E3 * t20 * t361 * t12 +
        0.7757360586169351E3 * t56 * t340 * t12 -
        0.6209421261849708E4 * t37 * t109 + 0.3852298453580824E4 * t216 * t370 +
        0.8114381744016598E2 * t163 * t373 * t15 -
        0.4423616526402223E3 * t377 * t261 * t13 -
        0.162287634880332E3 * t161 * t158 * t352 * t156 -
        0.162287634880332E3 * t160 * t158 * t352 * t156 -
        0.1090231648691615E4 * t232 * t16 +
        0.8229222048181555E4 * t393 * t63 * t15 +
        0.3852298453580824E4 * t216 * t397 + 0.3464801453452092E3 * t143 * t15 +
        0.218046329738323E4 * t56 * t150 * t13 +
        0.2765105050392224E3 * t229 * t405 * t5 -
        0.1039346859414405E5 * t41 * t13;
      double t445 = t24 * t15;
      double t466 = t127 * t98;
      MapleGenVar1 =
        -0.8114381744016598E2 * t163 * t373 * t5 -
        0.981795706946094E3 * t91 * t138 * t15 +
        0.8184481776377804E4 * t100 * t48 * t25 * t45 +
        0.3348614871593089E3 * t20 * t396 * t12 -
        0.7704596907161649E4 * t215 * t169 * t98 * t397 +
        0.5530210100784448E3 * t229 * t188 + 0.5530210100784448E3 * t281 * t188 -
        0.8229222048181555E4 * t393 * t63 * t5 +
        0.1326475405189345E3 * t37 * t30 * t15 +
        0.1548559945889136E4 * t337 * t94;
      double t470 =
        MapleGenVar1 - 0.1337958216886285E3 * t5 +
        0.1053422055706348E3 * t229 * t189 * t5 -
        0.1643351677346169E4 * t124 * t445 +
        0.1548559945889136E4 * t70 * t250 * t13 +
        0.1548559945889136E4 * t70 * t63 * t13 +
        0.1383226878248991E4 * t392 * t63 * t13 +
        0.8229222048181555E4 * t393 * t62 * t4 * t12 -
        0.1548559945889136E4 * t337 * t68 +
        0.2117850000189286E4 * t98 * t396 * t12 -
        0.7704596907161649E4 * t215 * t466 * t397;
      double t474 = t228 * t191;
      double t484 = t80 * t46;
      double t505 = t130 * t377;
      double t514 = t78 * t4;
      double t525 = 0.1046539361445834E4 * t397 + 0.1046539361445834E4 * t370 +
        0.5530210100784448E3 * t474 * t188 -
        0.7704596907161649E4 * t276 * t466 * t397 +
        0.1090231648691615E4 * t232 * t9 + 0.8229222048181555E4 * t393 * t76 +
        0.1798092235206523E4 * t484 * t299 + 0.1671445337926057E4 * t135 * t178 -
        0.162287634880332E3 * t163 * t291 * t4 * t12 +
        0.2106844111412697E3 * t206 * t13 - 0.1798092235206523E4 * t484 * t445 -
        0.1798092235206523E4 * t484 * t13 -
        0.1798092235206523E4 * t80 * t197 * t13 +
        0.2106844111412697E3 * t474 * t13 -
        0.7183122581705362E3 * t505 * t127 * t5 -
        0.1276149101998645E5 * t79 * t99 * t13 -
        0.7704596907161649E4 * t377 * t466 * t397 -
        0.1095666423317048E4 * t79 * t514 * t12 +
        0.874402992887409E3 * t229 * t189 * t4 * t12 -
        0.1798092235206523E4 * t80 * t26 * t13;
      double t552 = t20 * t8;
      double t572 =
        -0.1671445337926057E4 * t100 * t48 * t97 * t13 -
        0.2211808263201112E3 * t216 * t9 + 0.5917158602111804E2 * t20 * t15 -
        0.2765105050392224E3 * t229 * t405 * t15 +
        0.7183122581705362E3 * t505 * t127 * t15 -
        0.2658820189702891E3 * t41 * t15 - 0.8601696097502358E3 * t80 * t5 -
        0.1720339219500472E4 * t514 * t12 + 0.8601696097502358E3 * t80 * t15 -
        0.5294695384546455E3 * t20 * t4 * t12 + 0.1058939076909291E4 * t188 -
        0.1674329693194788E3 * t552 * t5 + 0.1674329693194788E3 * t552 * t15 -
        0.1046539361445834E4 * t316 - 0.1720339219500472E4 * t79 * t4 * t12 +
        0.6380745509993226E4 * t180 * t235 + 0.1046539361445834E4 * t214 +
        0.981795706946094E3 * t91 * t138 * t5 -
        0.1436624516341072E4 * t505 * t13 -
        0.1636896355275561E5 * t100 * t46 * t45;
      double t575 = t112 + t209 + t284 + t346 + t411 + t470 + t525 + t572;
      return t575;
    }

    // * f55 *********************************************************************

    double
    ortho2_f55 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * y;
      double t8 = -t3 - t7 - 0.5278640450004206E-1;
      double t10 = -t3 - t7 - 0.9472135954999579;
      double t14 = 0.1E1 * x;
      double t15 = 0.9472135954999579 + t14 + t1;
      double t16 = t15 * t6;
      double t17 = 0.5278640450004206E-1 + t14 + t1;
      double t23 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t24 = t23 * t6;
      double t27 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t31 = t6 * t5;
      double t32 = -t3 - t7 + 0.1546536707079771;
      double t33 = -t3 - t7 - 0.5;
      double t35 = -t3 - t7 - 0.1154653670707977E1;
      double t36 = t35 * t33 * t32;
      double t39 = 0.1154653670707977E1 + t14 + t1;
      double t40 = 0.5 + t14 + t1;
      double t41 = t40 * t39;
      double t42 = -0.1546536707079771 + t14 + t1;
      double t46 = t17 * t15;
      double t59 = -t3 - t7 + 0.2650553239294647;
      double t60 = -t3 - t7 - 0.2147684835193549;
      double t62 = -t3 - t7 - 0.7852315164806451;
      double t63 = -t3 - t7 - 0.1265055323929465E1;
      double t65 = t63 * t62 * t60 * t59;
      double t68 = t10 * t8;
      double t72 = 0.1265055323929465E1 + t14 + t1;
      double t73 = 0.7852315164806451 + t14 + t1;
      double t75 = 0.2147684835193549 + t14 + t1;
      double t76 = -0.2650553239294647 + t14 + t1;
      double t82 = t35 * t33;
      double MapleGenVar1 =
        0.1139640535344328E4 * t10 * t8 * t6 * t5 +
        0.2549157879614935E4 * t17 * t16 * t5 +
        0.3023413771542384E3 * t27 * t24 * t5 - 0.7228434489313551E3 * t36 * t31 -
        0.3072520170666446E4 * t42 * t41 * t31 +
        0.1583653772747944E4 * t27 * t46 * t31 +
        0.9798094737081272E3 * t10 * t8 * t23 * t31;
      double t100 =
        MapleGenVar1 + 0.2963461994798898E3 * t31 -
        0.3336873502349874E4 * t27 * t42 * t41 * t31 +
        0.5115703995204265E3 * t65 * t31 +
        0.8557537406959445E4 * t68 * t46 * t31 +
        0.5292175668942183E4 * t76 * t75 * t73 * t72 * t31 -
        0.1033579062670333E4 * t82 * t32 * t23 * t31 +
        0.3153661601350457E4 * (-0.3717401485096066 + t14 +
              t1) * (-0.917001814331423E-1 + t14 +
               t1) * (0.2907007820975211 + t14 +
                t1) * (0.7092992179024789 + t14 +
                 t1) * (0.1091700181433142E1 +
                  t14 +
                  t1) *
        (0.1371740148509607E1 + t14 + t1) * t6 * t5;
      double t116 = t39 * t6 * t5;
      double t117 = t42 * t40;
      double t123 = (0.1330223896278567E1 + t14 + t1) * t6 * t5;
      double t125 = t40 * (0.9688487934707142 + t14 + t1);
      double t128 =
        (-0.3302238962785669 + t14 + t1) * (0.3115120652928579E-1 + t14 + t1);
      double t133 = t16 * t5;
      double t141 = t72 * t6 * t5;
      double t142 = t75 * t73;
      double t148 = t24 * t5;
      double t149 = -t3 - t7 + 0.3302238962785669;
      double t150 = -t3 - t7 - 0.3115120652928579E-1;
      double t152 = -t3 - t7 - 0.9688487934707142;
      double t154 = -t3 - t7 - 0.1330223896278567E1;
      MapleGenVar1 =
        0.4139529848926616E2 * (-t3 - t7 - 0.1371740148509607E1) * (-t3 - t7 -
                    0.1091700181433142E1)
        * (-t3 - t7 - 0.7092992179024789) * (-t3 - t7 -
               0.2907007820975211) * (-t3 - t7 +
                    0.917001814331423E-1)
        * (-t3 - t7 + 0.3717401485096066) * t6 * t5 +
        0.7561406648554379E4 * t36 * t117 * t116 -
        0.3168514940764076E4 * t27 * t128 * t125 * t123 +
        0.2494366772317908E4 * t63 * t62 * t60 * t59 * t17 * t133 +
        0.1231999908755893E5 * t10 * t8 * t76 * t142 * t141 -
        0.1426368092188121E3 * t154 * t152 * t33 * t150 * t149 * t148 -
        0.2570313390091366E4 * t128 * t125 * t123;
      double t186 =
        MapleGenVar1 - 0.3736999873970488E4 * t82 * t32 * t17 * t133 +
        0.3001693060286986E3 * t65 * t148 -
        0.7094786670845755E4 * t68 * t117 * t116 +
        0.212546623830242E4 * t27 * t76 * t142 * t141 -
        0.8584650399287748E2 * t154 * t152 * t33 * t150 * t149 * t6 * t5 +
        0.2930929349370947E3 * t148 + 0.2135270163889491E3 * t27 * t6 * t5;
      double t187 = t100 + t186;
      return t187;
    }

    double
    ortho2_f55x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.1154653670707977E1 + t8 + t1;
      double t10 = 0.5 + t8 + t1;
      double t11 = t10 * t9;
      double t14 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t18 = -0.1546536707079771 + t8 + t1;
      double t19 = t18 * t11;
      double t22 = t6 * t2;
      double t23 = 0.1E1 * y;
      double t24 = -t3 - t23 + 0.2650553239294647;
      double t26 = -t3 - t23 - 0.2147684835193549;
      double t27 = -t3 - t23 - 0.7852315164806451;
      double t28 = t27 * t26;
      double t29 = -t3 - t23 - 0.1265055323929465E1;
      double t30 = t29 * t28;
      double t33 = -t3 - t23 + 0.1546536707079771;
      double t34 = -t3 - t23 - 0.5;
      double t35 = t34 * t33;
      double t36 = -t3 - t23 - 0.1154653670707977E1;
      double t37 = t36 * t35;
      double t40 = t18 * t9;
      double t44 = t18 * t10;
      double t45 = t14 * t44;
      double t48 = t26 * t24;
      double t66 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t67 = t66 * t22;
      double t70 = t9 * t5;
      double t73 = t9 * t22;
      double t76 = t66 * t5;
      double t79 = t33 * t66;
      double t90 = 0.9472135954999579 + t8 + t1;
      double t91 = 0.5278640450004206E-1 + t8 + t1;
      double t92 = t91 * t90;
      double t93 = -t3 - t23 - 0.9472135954999579;
      double t97 = -t3 - t23 - 0.5278640450004206E-1;
      double t102 =
        -0.3336873502349874E4 * t14 * t11 * t7 - 0.5276060265644412E4 * t19 * t7 -
        0.2557851997602132E3 * t30 * t24 * t22 + 0.3268463979900167E4 * t37 * t7 -
        0.3336873502349874E4 * t14 * t40 * t7 - 0.3336873502349874E4 * t45 * t7 -
        0.2557851997602132E3 * t27 * t48 * t7 -
        0.2557851997602132E3 * t29 * t48 * t7 -
        0.2557851997602132E3 * t29 * t27 * t24 * t7 -
        0.2557851997602132E3 * t30 * t7 + 0.2557851997602132E3 * t30 * t24 * t5 +
        0.5167895313351665E3 * t37 * t67 - 0.1668436751174937E4 * t45 * t70 +
        0.1668436751174937E4 * t45 * t73 - 0.5167895313351665E3 * t37 * t76 +
        0.5167895313351665E3 * t34 * t79 * t7 +
        0.5167895313351665E3 * t36 * t79 * t7 +
        0.5167895313351665E3 * t36 * t34 * t66 * t7 -
        0.4278768703479722E4 * t93 * t92 * t7 +
        0.8557537406959445E4 * t93 * t97 * t90 * t7;
      double t104 = t93 * t97 * t91;
      double t107 = t90 * t5;
      double t110 = t90 * t22;
      double t113 = 0.1265055323929465E1 + t8 + t1;
      double t114 = 0.7852315164806451 + t8 + t1;
      double t115 = t114 * t113;
      double t116 = 0.2147684835193549 + t8 + t1;
      double t120 = -0.2650553239294647 + t8 + t1;
      double t124 = t116 * t113;
      double t128 = t116 * t114;
      double t129 = t120 * t128;
      double t132 = t113 * t5;
      double t135 = t113 * t22;
      double t142 = t120 * t116;
      double t154 = t93 * t97;
      double t161 = t154 * t44;
      double t168 = t29 * t27;
      double t169 = t168 * t48;
      double t172 = -t3 - t23 + 0.3302238962785669;
      double t174 = -t3 - t23 - 0.9688487934707142;
      double t175 = -t3 - t23 - 0.1330223896278567E1;
      double t176 = t175 * t174;
      double t177 = t176 * t34 * t172;
      double t180 = -t3 - t23 - 0.3115120652928579E-1;
      double t182 = t176 * t34 * t180;
      double t185 =
        0.3360657201412931E4 * t142 * t115 * t7 +
        0.3547393335422877E4 * t97 * t18 * t11 * t7 +
        0.3547393335422877E4 * t93 * t18 * t11 * t7 -
        0.7094786670845755E4 * t154 * t11 * t7 -
        0.7094786670845755E4 * t154 * t40 * t7 -
        0.7094786670845755E4 * t161 * t7 - 0.3547393335422877E4 * t161 * t70 +
        0.3547393335422877E4 * t161 * t73 - 0.9492186907227991E3 * t169 * t7 +
        0.4292325199643874E2 * t177 * t7 + 0.4292325199643874E2 * t182 * t7;
      double t194 = 0.1330223896278567E1 + t8 + t1;
      double t196 = 0.9688487934707142 + t8 + t1;
      double t197 = t10 * t196;
      double t198 = 0.3115120652928579E-1 + t8 + t1;
      double t199 = -0.3302238962785669 + t8 + t1;
      double t200 = t199 * t198;
      double t201 = t200 * t197;
      double t204 = t180 * t172;
      double t205 = t176 * t204;
      double t208 = t196 * t194;
      double t225 = t198 * t10;
      double t231 = t174 * t34 * t204;
      double t235 = t175 * t34 * t204;
      double t238 = t24 * t66;
      double t242 = t29 * t26;
      double t258 = t36 * t34;
      double t259 = t258 * t33 * t91;
      double t266 =
        0.4292325199643874E2 * t231 * t7 + 0.4292325199643874E2 * t235 * t7 -
        0.1500846530143493E3 * t28 * t238 * t7 -
        0.1500846530143493E3 * t242 * t238 * t7 -
        0.1500846530143493E3 * t168 * t238 * t7 -
        0.1500846530143493E3 * t168 * t26 * t66 * t7 +
        0.1500846530143493E3 * t169 * t76 - 0.1500846530143493E3 * t169 * t67 -
        0.1868499936985244E4 * t259 * t107 + 0.1868499936985244E4 * t259 * t110 -
        0.3736999873970488E4 * t259 * t7;
      double t278 = t36 * t33;
      double t282 = t14 * t120;
      double t283 = t282 * t128;
      double t310 = t14 * t199;
      double t311 = t310 * t225;
      double t321 = t194 * t6 * t5;
      double t332 = 0.1067635081944745E3 * t14 * t5 - 0.1465464674685473E3 * t67
        - 0.5892253836201683E3 * t7 - 0.1067635081944745E3 * t14 * t22 +
        0.1465464674685473E3 * t76 -
        0.3168514940764076E4 * t311 * t196 * t6 * t5 -
        0.1584257470382038E4 * t311 * t208 * t5 +
        0.1584257470382038E4 * t311 * t208 * t22 -
        0.3168514940764076E4 * t311 * t321 -
        0.3168514940764076E4 * t310 * t198 * t196 * t321 -
        0.3168514940764076E4 * t14 * t198 * t197 * t321;
      double t342 = 0.1371740148509607E1 + t8 + t1;
      double t344 = t342 * t6 * t5;
      double t345 = 0.1091700181433142E1 + t8 + t1;
      double t346 = 0.7092992179024789 + t8 + t1;
      double t347 = t346 * t345;
      double t348 = -0.917001814331423E-1 + t8 + t1;
      double t349 = -0.3717401485096066 + t8 + t1;
      double t350 = t349 * t348;
      double t354 = 0.2907007820975211 + t8 + t1;
      double t360 = t350 * t354 * t346;
      double t367 = t345 * t342;
      double t376 = t258 * t33 * t18;
      double t384 = t9 * t6;
      double t385 = t384 * t5;
      double t392 = t10 * t6;
      double t396 = -t3 - t23 + 0.3717401485096066;
      double t397 = -t3 - t23 + 0.917001814331423E-1;
      double t398 = t397 * t396;
      double t400 = -t3 - t23 - 0.2907007820975211;
      double t401 = -t3 - t23 - 0.7092992179024789;
      double t403 = -t3 - t23 - 0.1091700181433142E1;
      double t404 = -t3 - t23 - 0.1371740148509607E1;
      double t405 = t404 * t403;
      double t406 = t405 * t401 * t400;
      double t426 =
        0.378070332427719E4 * t376 * t11 * t5 + 0.1481730997399449E3 * t5 +
        0.7561406648554379E4 * t258 * t33 * t10 * t385 +
        0.7561406648554379E4 * t376 * t385 +
        0.7561406648554379E4 * t376 * t392 * t5 -
        0.2069764924463308E2 * t406 * t398 * t22 +
        0.2069764924463308E2 * t406 * t398 * t5 +
        0.4510571953103487E3 * t182 * t172 * t6 * t5 -
        0.2069764924463308E2 * t406 * t397 * t6 * t5 -
        0.378070332427719E4 * t35 * t44 * t385 -
        0.378070332427719E4 * t278 * t44 * t385;
      double t433 = t154 * t142;
      double t451 = t113 * t6 * t5;
      double t470 = t90 * t6;
      double t474 = t97 * t6;
      double t479 = t396 * t6 * t5;
      double t486 = t400 * t397;
      double t498 = t172 * t66;
      double t505 = t66 * t6;
      double t506 = t505 * t5;
      double t509 =
        -0.6159999543779463E4 * t97 * t120 * t128 * t451 +
        0.2503976473501098E4 * t91 * t470 * t5 -
        0.3098429609928548E4 * t93 * t474 * t5 -
        0.2069764924463308E2 * t405 * t401 * t397 * t479 -
        0.2069764924463308E2 * t406 * t479 -
        0.2069764924463308E2 * t405 * t486 * t479 -
        0.2069764924463308E2 * t404 * t401 * t486 * t479 -
        0.2069764924463308E2 * t403 * t401 * t486 * t479 -
        0.7131840460940606E2 * t182 * t498 * t5 +
        0.7131840460940606E2 * t182 * t498 * t22 +
        0.7131840460940606E2 * t177 * t506;
      double t520 = t92 * t22;
      double t523 = t92 * t5;
      double t526 = t91 * t6;
      double t527 = t526 * t5;
      double t530 = t470 * t5;
      double t531 = t24 * t91;
      double t556 = t14 * t66;
      double t567 =
        -0.1247183386158954E4 * t28 * t531 * t530 -
        0.1247183386158954E4 * t242 * t531 * t530 -
        0.5698202676721638E3 * t154 * t22 - 0.5698202676721638E3 * t93 * t6 * t5 +
        0.5698202676721638E3 * t154 * t5 - 0.1274578939807468E4 * t520 +
        0.1511706885771192E3 * t556 * t5 - 0.1511706885771192E3 * t556 * t22 +
        0.4780436913596953E3 * t506 - 0.9560873827193906E3 * t14 * t6 * t5 -
        0.5698202676721638E3 * t474 * t5;
      double t572 = t14 * t92;
      double t589 = t33 * t6;
      double t613 = t93 * t97 * t66;
      double t624 =
        0.3614217244656775E3 * t37 * t22 - 0.3614217244656775E3 * t37 * t5 +
        0.3614217244656775E3 * t34 * t589 * t5 - 0.1536260085333223E4 * t19 * t5 -
        0.3072520170666446E4 * t18 * t392 * t5 -
        0.3072520170666446E4 * t18 * t384 * t5 -
        0.3072520170666446E4 * t10 * t384 * t5 -
        0.4899047368540636E3 * t613 * t22 -
        0.4899047368540636E3 * t97 * t505 * t5 -
        0.4899047368540636E3 * t93 * t505 * t5 + 0.4899047368540636E3 * t613 * t5;
      double MapleGenVar2 =
        t102 + t185 - 0.1481730997399449E3 * t22 +
        0.1231999908755893E5 * t433 * t451 + 0.106273311915121E4 * t283 * t132 -
        0.106273311915121E4 * t283 * t135 +
        0.1583653772747944E4 * t14 * t470 * t5 +
        0.1231999908755893E5 * t154 * t128 * t451 + t509 + t266 +
        0.6159999543779463E4 * t433 * t115 * t5 +
        0.7131840460940606E2 * t182 * t506 -
        0.6159999543779463E4 * t433 * t115 * t22 + t624 +
        0.1583653772747944E4 * t14 * t526 * t5 -
        0.1576830800675228E4 * t360 * t367 * t22 -
        0.3168514940764076E4 * t310 * t197 * t321 +
        0.7131840460940606E2 * t205 * t506 +
        0.3153661601350457E4 * t350 * t347 * t344;
      double MapleGenVar3 =
        MapleGenVar2 + 0.3614217244656775E3 * t36 * t589 * t5 +
        0.5292175668942183E4 * t120 * t115 * t7 +
        0.5292175668942183E4 * t116 * t115 * t7 -
        0.378070332427719E4 * t258 * t44 * t385 +
        0.1576830800675228E4 * t360 * t367 * t5 -
        0.378070332427719E4 * t376 * t11 * t22 +
        0.5292175668942183E4 * t120 * t124 * t7 +
        0.212546623830242E4 * t282 * t124 * t7 -
        0.1247183386158954E4 * t168 * t531 * t530;
      double MapleGenVar1 =
        MapleGenVar3 - 0.4278768703479722E4 * t97 * t92 * t7 +
        0.7131840460940606E2 * t231 * t506 -
        0.4292325199643874E2 * t182 * t172 * t5 +
        0.4292325199643874E2 * t182 * t172 * t22 +
        0.212546623830242E4 * t282 * t115 * t7 +
        0.1285156695045683E4 * t201 * t194 * t22 +
        0.1868499936985244E4 * t278 * t92 * t7 -
        0.2570313390091366E4 * t200 * t208 * t7 +
        0.1868499936985244E4 * t35 * t92 * t7 +
        0.7131840460940606E2 * t235 * t506 -
        0.1285156695045683E4 * t201 * t194 * t5;
      MapleGenVar2 =
        MapleGenVar1 + 0.1868499936985244E4 * t258 * t92 * t7 -
        0.2570313390091366E4 * t225 * t208 * t7 + t426 -
        0.1247183386158954E4 * t169 * t520 + 0.1274578939807468E4 * t523 +
        0.1247183386158954E4 * t169 * t523 + 0.2549157879614935E4 * t527 +
        0.2494366772317908E4 * t169 * t527 + 0.2549157879614935E4 * t530 +
        0.2494366772317908E4 * t169 * t530 - 0.2570313390091366E4 * t201 * t7 +
        0.3153661601350457E4 * t360 * t344 - 0.5009862006543987E4 * t201 * t321 +
        0.4292325199643874E2 * t205 * t7 + t567 + t332 -
        0.7918268863739721E3 * t572 * t22 + 0.7918268863739721E3 * t572 * t5 +
        0.1536260085333223E4 * t19 * t22;
      MapleGenVar3 =
        MapleGenVar2 + 0.212546623830242E4 * t283 * t7 -
        0.2646087834471091E4 * t129 * t135 -
        0.2570313390091366E4 * t200 * t10 * t194 * t7 -
        0.2570313390091366E4 * t199 * t10 * t208 * t7 +
        0.2646087834471091E4 * t129 * t132 +
        0.212546623830242E4 * t14 * t116 * t115 * t7 -
        0.3736999873970488E4 * t258 * t33 * t90 * t7 +
        0.3153661601350457E4 * t360 * t345 * t6 * t5 +
        0.3153661601350457E4 * t350 * t354 * t345 * t344;
      double t628 =
        MapleGenVar3 + 0.1231999908755893E5 * t433 * t114 * t6 * t5 +
        0.3153661601350457E4 * t349 * t354 * t347 * t344 +
        0.3153661601350457E4 * t348 * t354 * t347 * t344 +
        0.1231999908755893E5 * t154 * t120 * t114 * t451 +
        0.5292175668942183E4 * t129 * t7 -
        0.6159999543779463E4 * t93 * t120 * t128 * t451 +
        0.3614217244656775E3 * t36 * t34 * t6 * t5 -
        0.1247183386158954E4 * t168 * t26 * t91 * t530 -
        0.4278768703479722E4 * t104 * t110 + 0.4278768703479722E4 * t104 * t107 +
        0.8557537406959445E4 * t104 * t7;
      return t628;
    }

    double
    ortho2_f55y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * y;
      double t8 = -t3 - t7 + 0.3717401485096066;
      double t10 = t8 * t6 * t5;
      double t11 = -t3 - t7 + 0.917001814331423E-1;
      double t12 = -t3 - t7 - 0.2907007820975211;
      double t13 = t12 * t11;
      double t14 = -t3 - t7 - 0.7092992179024789;
      double t15 = -t3 - t7 - 0.1371740148509607E1;
      double t20 = -t3 - t7 - 0.1091700181433142E1;
      double t25 = t6 * t5;
      double t28 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t29 = -t3 - t7 + 0.1546536707079771;
      double t30 = t29 * t28;
      double t31 = -t3 - t7 - 0.5;
      double t35 = -t3 - t7 - 0.2147684835193549;
      double t36 = -t3 - t7 - 0.7852315164806451;
      double t37 = t36 * t35;
      double t38 = -t3 - t7 - 0.1265055323929465E1;
      double t39 = t38 * t37;
      double t42 = t6 * t4;
      double t43 = 0.1E1 * x;
      double t44 = 0.1265055323929465E1 + t43 + t1;
      double t45 = t44 * t42;
      double t46 = 0.7852315164806451 + t43 + t1;
      double t47 = 0.2147684835193549 + t43 + t1;
      double t48 = t47 * t46;
      double t49 = -0.2650553239294647 + t43 + t1;
      double t52 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t53 = t52 * t49;
      double t54 = t53 * t48;
      double t57 = t6 * t2;
      double t58 = 0.9472135954999579 + t43 + t1;
      double t59 = t58 * t57;
      double t60 = 0.5278640450004206E-1 + t43 + t1;
      double t61 = -t3 - t7 - 0.5278640450004206E-1;
      double t63 = -t3 - t7 - 0.9472135954999579;
      double t64 = t63 * t61 * t60;
      double t69 = -t3 - t7 + 0.2650553239294647;
      double t74 = t35 * t69;
      double t78 = t49 * t48;
      double t81 = t28 * t6;
      double t85 = 0.1154653670707977E1 + t43 + t1;
      double t86 = t85 * t57;
      double t87 = 0.5 + t43 + t1;
      double t88 = -0.1546536707079771 + t43 + t1;
      double t89 = t88 * t87;
      double t90 = t52 * t89;
      double t95 = t60 * t58;
      double t96 = t95 * t57;
      double t98 = t58 * t6;
      double t99 = t98 * t5;
      double t101 = t60 * t6;
      double t102 = t101 * t5;
      double t104 = -t3 - t7 + 0.3302238962785669;
      double t106 = -t3 - t7 - 0.3115120652928579E-1;
      double t108 = -t3 - t7 - 0.9688487934707142;
      double t109 = -t3 - t7 - 0.1330223896278567E1;
      double t110 = t109 * t108;
      double t111 = t110 * t31 * t106;
      double t114 = t63 * t61;
      double t115 = t114 * t89;
      double t121 = t52 * t95;
      double t124 =
        -0.4139529848926616E2 * t15 * t14 * t13 * t10 -
        0.4139529848926616E2 * t20 * t14 * t13 * t10 +
        0.1033579062670333E4 * t31 * t30 * t25 -
        0.5115703995204265E3 * t39 * t25 + 0.106273311915121E4 * t54 * t45 -
        0.4278768703479722E4 * t64 * t59 + 0.4278768703479722E4 * t64 * t25 -
        0.5115703995204265E3 * t38 * t36 * t69 * t25 -
        0.5115703995204265E3 * t38 * t74 * t25 +
        0.2646087834471091E4 * t78 * t25 - 0.9798094737081272E3 * t61 * t81 * t5 +
        0.1668436751174937E4 * t90 * t86 + 0.2646087834471091E4 * t78 * t45 -
        0.1274578939807468E4 * t96 + 0.1274578939807468E4 * t99 +
        0.1274578939807468E4 * t102 - 0.4292325199643874E2 * t111 * t104 * t42 +
        0.3547393335422877E4 * t115 * t86 -
        0.2557851997602132E3 * t39 * t69 * t57 +
        0.7918268863739721E3 * t121 * t42;
      double t125 = t85 * t6;
      double t132 = t28 * t42;
      double t134 = t87 * t85;
      double t135 = t88 * t134;
      double t146 = t44 * t6 * t5;
      double t155 = t63 * t61 * t28;
      double t164 = t104 * t28;
      double t168 = t46 * t44;
      double t169 = t49 * t47;
      double t176 = t61 * t6;
      double t180 = 0.1330223896278567E1 + t43 + t1;
      double t182 = 0.9688487934707142 + t43 + t1;
      double t183 = t87 * t182;
      double t184 = 0.3115120652928579E-1 + t43 + t1;
      double t185 = -0.3302238962785669 + t43 + t1;
      double t186 = t185 * t184;
      double t187 = t186 * t183;
      double t191 = t114 * t169;
      double t194 = t88 * t85;
      double t198 = t95 * t42;
      double t202 =
        -0.1668436751174937E4 * t52 * t134 * t25 -
        0.1668436751174937E4 * t90 * t25 +
        0.7131840460940606E2 * t111 * t164 * t57 +
        0.6721314402825863E4 * t169 * t168 * t25 +
        0.5007952947002195E4 * t60 * t98 * t5 -
        0.1549214804964274E4 * t63 * t176 * t5 +
        0.1285156695045683E4 * t187 * t180 * t57 -
        0.6159999543779463E4 * t191 * t168 * t57 -
        0.3547393335422877E4 * t114 * t194 * t25 + 0.1274578939807468E4 * t198 -
        0.5698202676721638E3 * t114 * t57;
      double t209 = -t3 - t7 - 0.1154653670707977E1;
      double t213 = t29 * t6;
      double t224 = t38 * t36;
      double t225 = t224 * t74;
      double t228 = t58 * t42;
      double t231 = 0.1371740148509607E1 + t43 + t1;
      double t233 = t231 * t6 * t5;
      double t234 = 0.1091700181433142E1 + t43 + t1;
      double t235 = 0.2907007820975211 + t43 + t1;
      double t237 = -0.917001814331423E-1 + t43 + t1;
      double t238 = -0.3717401485096066 + t43 + t1;
      double t239 = t238 * t237;
      double t243 = t69 * t60;
      double t244 = t38 * t35;
      double t252 = t11 * t8;
      double t255 = t15 * t20;
      double t256 = t255 * t14 * t12;
      double t259 = t31 * t29;
      double t260 = t209 * t259;
      double t263 = t85 * t42;
      double t266 = t28 * t57;
      double t282 = t81 * t5;
      double t285 =
        -0.2494366772317908E4 * t37 * t243 * t99 -
        0.2069764924463308E2 * t256 * t252 * t57 +
        0.1634231989950084E4 * t260 * t25 - 0.1668436751174937E4 * t90 * t263 +
        0.5167895313351665E3 * t260 * t266 - 0.1067635081944745E3 * t52 * t57 -
        0.1465464674685473E3 * t266 + 0.1067635081944745E3 * t52 * t42 +
        0.1033579062670333E4 * t209 * t31 * t28 * t25 -
        0.4139529848926616E2 * t256 * t11 * t6 * t5 +
        0.1426368092188121E3 * t111 * t282;
      double t288 = t110 * t31 * t104;
      double t291 = 0.7092992179024789 + t43 + t1;
      double t292 = t291 * t234;
      double t297 = t47 * t44;
      double t349 = t182 * t180;
      double t351 = t184 * t87;
      double t352 = t52 * t185;
      double t353 = t352 * t351;
      double t360 = t234 * t231;
      double t363 = t239 * t235 * t291;
      double t366 =
        -0.1285156695045683E4 * t187 * t180 * t42 -
        0.4139529848926616E2 * t256 * t10 +
        0.1033579062670333E4 * t209 * t30 * t25 +
        0.2646087834471091E4 * t49 * t168 * t25 +
        0.2646087834471091E4 * t47 * t168 * t25 -
        0.8557537406959445E4 * t63 * t95 * t25 +
        0.6159999543779463E4 * t191 * t146 +
        0.1576830800675228E4 * t239 * t292 * t233 +
        0.1584257470382038E4 * t353 * t349 * t57 -
        0.1584257470382038E4 * t353 * t182 * t6 * t5 +
        0.1576830800675228E4 * t363 * t360 * t42;
      double t370 = t44 * t57;
      double t378 = t87 * t6;
      double t381 = t209 * t31;
      double t382 = t381 * t29 * t88;
      double t385 = t125 * t5;
      double t391 = t106 * t104;
      double t393 = t109 * t31 * t391;
      double t401 = t180 * t6 * t5;
      double t430 = t110 * t391;
      double t443 =
        -0.1285156695045683E4 * t185 * t87 * t349 * t25 -
        0.7561406648554379E4 * t381 * t89 * t385 +
        0.106273311915121E4 * t54 * t25 -
        0.3001693060286986E3 * t224 * t35 * t28 * t25 +
        0.7094786670845755E4 * t63 * t88 * t134 * t25 -
        0.1500846530143493E3 * t225 * t266 +
        0.2069764924463308E2 * t256 * t252 * t42 +
        0.1426368092188121E3 * t430 * t282 - 0.5167895313351665E3 * t260 * t132 +
        0.7094786670845755E4 * t61 * t88 * t134 * t25 -
        0.1868499936985244E4 * t381 * t29 * t58 * t25;
      double t453 = t209 * t29;
      double t474 = t69 * t28;
      double t493 = t52 * t28;
      double t504 =
        -0.3001693060286986E3 * t224 * t474 * t25 +
        0.378070332427719E4 * t381 * t29 * t87 * t385 -
        0.1139640535344328E4 * t63 * t6 * t5 - 0.1139640535344328E4 * t176 * t5 +
        0.5698202676721638E3 * t114 * t42 - 0.4780436913596953E3 * t52 * t6 * t5 +
        0.9560873827193906E3 * t282 - 0.1511706885771192E3 * t493 * t57 +
        0.1511706885771192E3 * t493 * t42 +
        0.106273311915121E4 * t52 * t47 * t168 * t25 -
        0.1247183386158954E4 * t225 * t96;
      double t512 = t381 * t29 * t60;
      double t539 = t108 * t31 * t391;
      double t566 =
        -0.1584257470382038E4 * t353 * t349 * t42 +
        0.1426368092188121E3 * t393 * t282 + 0.1426368092188121E3 * t539 * t282 +
        0.3614217244656775E3 * t260 * t57 +
        0.7918268863739721E3 * t52 * t101 * t5 -
        0.7561406648554379E4 * t259 * t89 * t385 -
        0.7561406648554379E4 * t453 * t89 * t385 -
        0.2494366772317908E4 * t224 * t35 * t60 * t99 -
        0.2494366772317908E4 * t224 * t243 * t99 +
        0.378070332427719E4 * t382 * t134 * t42 -
        0.1536260085333223E4 * t88 * t378 * t5;
      double t625 =
        -0.1584257470382038E4 * t352 * t184 * t182 * t401 -
        0.1584257470382038E4 * t352 * t183 * t401 -
        0.1584257470382038E4 * t52 * t184 * t183 * t401 -
        0.3547393335422877E4 * t115 * t263 + 0.8584650399287748E2 * t539 * t25 -
        0.1576830800675228E4 * t363 * t360 * t57 +
        0.1576830800675228E4 * t363 * t234 * t6 * t5 +
        0.1576830800675228E4 * t363 * t233 +
        0.6159999543779463E4 * t191 * t168 * t42 +
        0.1868499936985244E4 * t512 * t59 - 0.1868499936985244E4 * t512 * t25;
      double MapleGenVar2 =
        -0.7918268863739721E3 * t121 * t57 + 0.1500846530143493E3 * t225 * t132 -
        0.2646087834471091E4 * t78 * t370 + 0.8584650399287748E2 * t430 * t25 +
        0.1247183386158954E4 * t225 * t102 - 0.1481730997399449E3 * t57 +
        0.1247183386158954E4 * t225 * t198 + 0.4278768703479722E4 * t64 * t228 +
        0.3736999873970488E4 * t259 * t95 * t25 -
        0.7131840460940606E2 * t111 * t164 * t42 -
        0.1285156695045683E4 * t351 * t349 * t25 -
        0.1536260085333223E4 * t88 * t125 * t5 -
        0.1536260085333223E4 * t87 * t125 * t5 + t366 -
        0.106273311915121E4 * t54 * t370 + t202 -
        0.9798094737081272E3 * t63 * t81 * t5 +
        0.2646087834471091E4 * t49 * t297 * t25 -
        0.5115703995204265E3 * t36 * t74 * t25;
      double MapleGenVar1 =
        MapleGenVar2 + 0.7228434489313551E3 * t209 * t213 * t5 -
        0.3547393335422877E4 * t114 * t134 * t25 +
        0.4292325199643874E2 * t111 * t104 * t57 +
        0.7918268863739721E3 * t52 * t98 * t5 +
        0.1576830800675228E4 * t238 * t235 * t292 * t233 +
        0.1576830800675228E4 * t237 * t235 * t292 * t233 + t504 -
        0.1536260085333223E4 * t135 * t42 - 0.4899047368540636E3 * t155 * t57 -
        0.1584257470382038E4 * t353 * t401 -
        0.1668436751174937E4 * t52 * t194 * t25 -
        0.8557537406959445E4 * t61 * t95 * t25 +
        0.378070332427719E4 * t382 * t378 * t5 +
        0.4899047368540636E3 * t155 * t42 +
        0.6159999543779463E4 * t114 * t49 * t46 * t146 +
        0.7228434489313551E3 * t209 * t31 * t6 * t5 +
        0.8584650399287748E2 * t393 * t25 - 0.3614217244656775E3 * t260 * t42 -
        0.1055212053128882E5 * t135 * t25 + 0.1247183386158954E4 * t225 * t99;
      MapleGenVar2 =
        MapleGenVar1 - 0.1285156695045683E4 * t186 * t349 * t25 +
        0.6159999543779463E4 * t114 * t48 * t146 -
        0.3001693060286986E3 * t244 * t474 * t25 +
        0.6159999543779463E4 * t191 * t46 * t6 * t5 -
        0.1285156695045683E4 * t186 * t87 * t180 * t25 -
        0.1231999908755893E5 * t63 * t49 * t48 * t146 -
        0.1231999908755893E5 * t61 * t49 * t48 * t146 +
        0.2255285976551744E3 * t111 * t104 * t6 * t5 -
        0.4139529848926616E2 * t255 * t14 * t11 * t10 +
        0.2118110935168016E3 * t25 + 0.378070332427719E4 * t382 * t385 + t124 +
        0.1465464674685473E3 * t132 + 0.7228434489313551E3 * t31 * t213 * t5 +
        0.4278768703479722E4 * t63 * t61 * t58 * t25 + t566 + t285 +
        0.1426368092188121E3 * t288 * t282 -
        0.3001693060286986E3 * t37 * t474 * t25;
      double t629 =
        MapleGenVar2 + 0.106273311915121E4 * t53 * t168 * t25 -
        0.2494366772317908E4 * t244 * t243 * t99 +
        0.1536260085333223E4 * t135 * t57 - 0.3547393335422877E4 * t115 * t25 -
        0.1001972401308797E5 * t187 * t401 + t443 +
        0.8584650399287748E2 * t288 * t25 - 0.1285156695045683E4 * t187 * t25 -
        0.4746093453613996E3 * t225 * t25 + 0.1481730997399449E3 * t42 +
        0.1576830800675228E4 * t239 * t235 * t234 * t233 + t625 -
        0.1868499936985244E4 * t512 * t228 +
        0.106273311915121E4 * t53 * t297 * t25 -
        0.4139529848926616E2 * t255 * t13 * t10 +
        0.3736999873970488E4 * t381 * t95 * t25 +
        0.2557851997602132E3 * t39 * t69 * t42 -
        0.378070332427719E4 * t382 * t134 * t57 +
        0.8584650399287748E2 * t111 * t25 +
        0.3736999873970488E4 * t453 * t95 * t25;
      return t629;
    }

    // ORDER 10

    // Edge functions, order 10

    // number 56
    inline double
    ortho2_f56 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l2 * l3 * phi8 (l3 - l2);
    }

    inline double
    ortho2_f56x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l2x * l3 + l2 * l3x) * phi8 (l3 - l2) + l2 * l3 * phi8x (l3 -
                       l2) *
        (l3x - l2x);
    }

    inline double
    ortho2_f56y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l2y * l3 + l2 * l3y) * phi8 (l3 - l2) + l2 * l3 * phi8x (l3 -
                       l2) *
        (l3y - l2y);
    }

    // number 57
    inline double
    ortho2_f57 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l3 * l1 * phi8 (l1 - l3);
    }

    inline double
    ortho2_f57x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l3x * l1 + l3 * l1x) * phi8 (l1 - l3) + l3 * l1 * phi8x (l1 -
                       l3) *
        (l1x - l3x);
    }

    inline double
    ortho2_f57y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l3y * l1 + l3 * l1y) * phi8 (l1 - l3) + l3 * l1 * phi8x (l1 -
                       l3) *
        (l1y - l3y);
    }

    // number 58
    inline double
    ortho2_f58 (double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      return l1 * l2 * phi8 (l2 - l1);
    }

    inline double
    ortho2_f58x (double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1x = lambda1x (x, y);
      l2x = lambda2x (x, y);
      l3x = lambda3x (x, y);
      return (l1x * l2 + l1 * l2x) * phi8 (l2 - l1) + l1 * l2 * phi8x (l2 -
                       l1) *
        (l2x - l1x);
    }

    inline double
    ortho2_f58y (double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y;
      l1 = lambda1 (x, y);
      l2 = lambda2 (x, y);
      l3 = lambda3 (x, y);
      l1y = lambda1y (x, y);
      l2y = lambda2y (x, y);
      l3y = lambda3y (x, y);
      return (l1y * l2 + l1 * l2y) * phi8 (l2 - l1) + l1 * l2 * phi8x (l2 -
                       l1) *
        (l2y - l1y);
    }

    // Bubble functions, order 10

    // * f59 *********************************************************************

    double
    ortho2_f59 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * y;
      double t14 = -t3 - t7 - 0.5;
      double t24 = t6 * t5;
      double t25 = -t3 - t7 + 0.1546536707079771;
      double t27 = -t3 - t7 - 0.1154653670707977E1;
      double t28 = t27 * t14 * t25;
      double t31 = 0.1E1 * x;
      double t32 = 0.1154653670707977E1 + t31 + t1;
      double t33 = 0.5 + t31 + t1;
      double t34 = t33 * t32;
      double t35 = -0.1546536707079771 + t31 + t1;
      double t39 = 0.9472135954999579 + t31 + t1;
      double t40 = 0.5278640450004206E-1 + t31 + t1;
      double t41 = t40 * t39;
      double t44 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t50 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t51 = -t3 - t7 - 0.5278640450004206E-1;
      double t53 = -t3 - t7 - 0.9472135954999579;
      double t57 = t50 * t6;
      double t65 = t39 * t6;
      double t70 = -t3 - t7 + 0.2650553239294647;
      double t71 = -t3 - t7 - 0.2147684835193549;
      double t73 = -t3 - t7 - 0.7852315164806451;
      double t74 = -t3 - t7 - 0.1265055323929465E1;
      double t76 = t74 * t73 * t71 * t70;
      double t79 = t53 * t51;
      double t83 = 0.1265055323929465E1 + t31 + t1;
      double t84 = 0.7852315164806451 + t31 + t1;
      double t86 = 0.2147684835193549 + t31 + t1;
      double t87 = -0.2650553239294647 + t31 + t1;
      double t97 = t27 * t14;
      double MapleGenVar1 =
        0.22250594477178E4 * (-t3 - t7 - 0.139975799541146E1) * (-t3 - t7 -
                       0.1177186279510738E1)
        * (-t3 - t7 - 0.8631174638261782) * t14 * (-t3 - t7 -
                     0.1368825361738218) * (-t3 -
                          t7 +
                          0.1771862795107378)
        * (-t3 - t7 + 0.3997579954114602) * t6 * t5 +
        0.1040372256714634E4 * t28 * t24 -
        0.4426934589560384E1 * t35 * t34 * t24 -
        0.923693451666489E2 * t44 * t41 * t24 -
        0.1892762020444952E3 * t53 * t51 * t50 * t24 +
        0.3695738687939205E2 * t44 * t57 * t5 +
        0.4115914087765948E3 * t53 * t51 * t6 * t5;
      double t101 =
        MapleGenVar1 + 0.1640452312932599E2 * t40 * t65 * t5 +
        0.3972171851120012E2 * t24 + 0.1291046144090737E4 * t76 * t24 +
        0.6219545471548346E2 * t79 * t41 * t24 -
        0.6367288425307034E1 * t87 * t86 * t84 * t83 * t24 -
        0.1071108804605722E2 * t44 * t35 * t34 * t24 -
        0.5281691920541455E3 * t97 * t25 * t50 * t24;
      double t103 = t32 * t6 * t5;
      double t104 = t35 * t33;
      double t110 = (0.1330223896278567E1 + t31 + t1) * t6 * t5;
      double t112 = t33 * (0.9688487934707142 + t31 + t1);
      double t115 =
        (-0.3302238962785669 + t31 + t1) * (0.3115120652928579E-1 + t31 + t1);
      double t120 = t65 * t5;
      double t156 = t83 * t6 * t5;
      double t157 = t86 * t84;
      double t163 = t57 * t5;
      double t164 = -t3 - t7 + 0.3302238962785669;
      double t165 = -t3 - t7 - 0.3115120652928579E-1;
      double t167 = -t3 - t7 - 0.9688487934707142;
      double t169 = -t3 - t7 - 0.1330223896278567E1;
      double MapleGenVar2 =
        -0.1045764578632491E3 * t28 * t104 * t103 +
        0.2330296513480324E1 * t44 * t115 * t112 * t110 -
        0.1195880576620191E3 * t74 * t73 * t71 * t70 * t40 * t120;
      MapleGenVar1 =
        MapleGenVar2 - 0.1294797561569409E1 * (-0.3717401485096066 + t31 +
                 t1) * (-0.917001814331423E-1 +
                  t31 +
                  t1) * (0.2907007820975211 +
                   t31 +
                   t1) *
        (0.7092992179024789 + t31 + t1) * (0.1091700181433142E1 + t31 +
                   t1) * (0.1371740148509607E1 + t31 +
                    t1) * t6 * t5 +
        0.1669421501547021E4 * (-t3 - t7 - 0.1371740148509607E1) * (-t3 - t7 -
                    0.1091700181433142E1)
        * (-t3 - t7 - 0.7092992179024789) * (-t3 - t7 -
               0.2907007820975211) * (-t3 - t7 +
                    0.917001814331423E-1)
        * (-t3 - t7 + 0.3717401485096066) * t6 * t5 -
        0.289813189885209E2 * t53 * t51 * t87 * t157 * t156 -
        0.9322085311855808E3 * t169 * t167 * t14 * t165 * t164 * t163;
      double t201 =
        MapleGenVar1 - 0.3680091408067806E3 * t76 * t163 +
        0.1808378606046524 * t115 * t112 * t110 +
        0.9894140496190919E3 * t97 * t25 * t40 * t120 +
        0.3419671758443138E1 * t44 * t87 * t157 * t156 -
        0.273362156439465E2 * t79 * t104 * t103 +
        0.2022426084311444E4 * t169 * t167 * t14 * t165 * t164 * t6 * t5 -
        0.6708680441856585E2 * t44 * t6 * t5 - 0.2079672579801747E2 * t163;
      double t202 = t101 + t201;
      return t202;
    }

    double
    ortho2_f59x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.1154653670707977E1 + t8 + t1;
      double t10 = 0.5 + t8 + t1;
      double t11 = t10 * t9;
      double t14 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t18 = -0.1546536707079771 + t8 + t1;
      double t19 = t18 * t11;
      double t22 = 0.1E1 * y;
      double t23 = -t3 - t22 + 0.1546536707079771;
      double t24 = -t3 - t22 - 0.5;
      double t25 = t24 * t23;
      double t26 = -t3 - t22 - 0.1154653670707977E1;
      double t27 = t26 * t25;
      double t30 = t6 * t2;
      double t31 = -t3 - t22 + 0.2650553239294647;
      double t33 = -t3 - t22 - 0.2147684835193549;
      double t34 = -t3 - t22 - 0.7852315164806451;
      double t35 = t34 * t33;
      double t36 = -t3 - t22 - 0.1265055323929465E1;
      double t37 = t36 * t35;
      double t40 = t18 * t10;
      double t41 = t14 * t40;
      double t47 = t18 * t9;
      double t57 = t33 * t31;
      double t64 =
        -0.1071108804605722E2 * t14 * t11 * t7 - 0.1693571722207165E2 * t19 * t7 +
        0.1670217636822007E4 * t27 * t7 - 0.6455230720453683E3 * t37 * t31 * t30 -
        0.1071108804605722E2 * t41 * t7 + 0.6455230720453683E3 * t37 * t31 * t5 -
        0.1071108804605722E2 * t14 * t47 * t7 -
        0.6455230720453683E3 * t36 * t34 * t31 * t7 -
        0.6455230720453683E3 * t37 * t7 - 0.6455230720453683E3 * t34 * t57 * t7 -
        0.6455230720453683E3 * t36 * t57 * t7;
      double t67 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t68 = t23 * t67;
      double t76 = t67 * t5;
      double t79 = t67 * t30;
      double t82 = t9 * t5;
      double t85 = t9 * t30;
      double t88 = 0.5278640450004206E-1 + t8 + t1;
      double t89 = -t3 - t22 - 0.5278640450004206E-1;
      double t91 = -t3 - t22 - 0.9472135954999579;
      double t92 = t91 * t89 * t88;
      double t95 = 0.9472135954999579 + t8 + t1;
      double t96 = t95 * t5;
      double t99 = t95 * t30;
      double t109 =
        0.2640845960270727E3 * t26 * t68 * t7 +
        0.2640845960270727E3 * t26 * t24 * t67 * t7 -
        0.2640845960270727E3 * t27 * t76 + 0.2640845960270727E3 * t27 * t79 -
        0.5355544023028608E1 * t41 * t82 + 0.5355544023028608E1 * t41 * t85 +
        0.6219545471548346E2 * t92 * t7 + 0.3109772735774173E2 * t92 * t96 -
        0.3109772735774173E2 * t92 * t99 + 0.2640845960270727E3 * t24 * t68 * t7 +
        0.6219545471548346E2 * t91 * t89 * t95 * t7;
      double t111 = 0.1265055323929465E1 + t8 + t1;
      double t112 = 0.7852315164806451 + t8 + t1;
      double t113 = t112 * t111;
      double t114 = 0.2147684835193549 + t8 + t1;
      double t118 = -0.2650553239294647 + t8 + t1;
      double t122 = t114 * t111;
      double t126 = t114 * t112;
      double t127 = t118 * t126;
      double t130 = t111 * t5;
      double t133 = t111 * t30;
      double t136 = t88 * t95;
      double t143 = t91 * t89;
      double t144 = t143 * t40;
      double t147 = t36 * t34;
      double t148 = t147 * t57;
      double t151 = -t3 - t22 - 0.3115120652928579E-1;
      double t153 = -t3 - t22 - 0.9688487934707142;
      double t154 = -t3 - t22 - 0.1330223896278567E1;
      double t155 = t154 * t153;
      double t156 = t155 * t24 * t151;
      double t159 =
        -0.6367288425307034E1 * t114 * t113 * t7 -
        0.6367288425307034E1 * t118 * t113 * t7 -
        0.6367288425307034E1 * t118 * t122 * t7 -
        0.6367288425307034E1 * t127 * t7 - 0.3183644212653517E1 * t127 * t130 +
        0.3183644212653517E1 * t127 * t133 -
        0.3109772735774173E2 * t89 * t136 * t7 -
        0.3109772735774173E2 * t91 * t136 * t7 +
        0.1366810782197325E2 * t144 * t85 + 0.1163747084711042E4 * t148 * t7 -
        0.1011213042155722E4 * t156 * t7;
      double t160 = -t3 - t22 + 0.3302238962785669;
      double t167 = t118 * t114;
      double t189 = 0.1330223896278567E1 + t8 + t1;
      double t191 = 0.9688487934707142 + t8 + t1;
      double t192 = t10 * t191;
      double t193 = 0.3115120652928579E-1 + t8 + t1;
      double t194 = -0.3302238962785669 + t8 + t1;
      double t195 = t194 * t193;
      double t196 = t195 * t192;
      double t199 = t151 * t160;
      double t200 = t155 * t199;
      double t203 =
        0.1011213042155722E4 * t156 * t160 * t5 -
        0.1011213042155722E4 * t156 * t160 * t30 +
        0.5406975803416726E1 * t167 * t113 * t7 +
        0.1366810782197325E2 * t89 * t18 * t11 * t7 +
        0.1366810782197325E2 * t91 * t18 * t11 * t7 -
        0.273362156439465E2 * t143 * t11 * t7 -
        0.273362156439465E2 * t143 * t47 * t7 - 0.273362156439465E2 * t144 * t7 -
        0.1366810782197325E2 * t144 * t82 -
        0.9041893030232619E-1 * t196 * t189 * t30 -
        0.1011213042155722E4 * t200 * t7;
      double t207 = t155 * t24 * t160;
      double t220 = t153 * t24 * t199;
      double t224 = t154 * t24 * t199;
      double t227 = t191 * t189;
      double t228 = t193 * t10;
      double t239 = t31 * t67;
      double t243 = t36 * t33;
      double t247 =
        -0.1011213042155722E4 * t207 * t7 +
        0.1808378606046524 * t195 * t10 * t189 * t7 +
        0.1808378606046524 * t196 * t7 +
        0.9041893030232619E-1 * t196 * t189 * t5 -
        0.1011213042155722E4 * t220 * t7 - 0.1011213042155722E4 * t224 * t7 +
        0.1808378606046524 * t228 * t227 * t7 +
        0.1808378606046524 * t194 * t10 * t227 * t7 +
        0.1808378606046524 * t195 * t227 * t7 +
        0.1840045704033903E3 * t35 * t239 * t7 +
        0.1840045704033903E3 * t243 * t239 * t7;
      double t260 = t26 * t24;
      double t261 = t260 * t23 * t88;
      double t275 = t26 * t23;
      double t282 =
        0.1840045704033903E3 * t147 * t239 * t7 +
        0.1840045704033903E3 * t147 * t33 * t67 * t7 -
        0.1840045704033903E3 * t148 * t76 + 0.1840045704033903E3 * t148 * t79 +
        0.9894140496190919E3 * t261 * t7 + 0.4947070248095459E3 * t261 * t96 -
        0.4947070248095459E3 * t261 * t99 +
        0.9894140496190919E3 * t260 * t23 * t95 * t7 -
        0.4947070248095459E3 * t25 * t136 * t7 -
        0.4947070248095459E3 * t275 * t136 * t7 -
        0.4947070248095459E3 * t260 * t136 * t7;
      double t284 = t14 * t118;
      double t285 = t284 * t126;
      double t302 = -t3 - t22 + 0.3997579954114602;
      double t304 = t302 * t6 * t5;
      double t305 = -t3 - t22 + 0.1771862795107378;
      double t307 = -t3 - t22 - 0.8631174638261782;
      double t308 = -t3 - t22 - 0.1177186279510738E1;
      double t310 = -t3 - t22 - 0.139975799541146E1;
      double t311 = t310 * t308 * t307;
      double t315 = -t3 - t22 - 0.1368825361738218;
      double t317 = t311 * t24 * t315;
      double t324 = t305 * t302;
      double t331 =
        0.1709835879221569E1 * t285 * t130 - 0.1709835879221569E1 * t285 * t133 +
        0.3419671758443138E1 * t14 * t114 * t113 * t7 +
        0.3419671758443138E1 * t284 * t113 * t7 +
        0.3419671758443138E1 * t284 * t122 * t7 +
        0.3419671758443138E1 * t285 * t7 -
        0.11125297238589E4 * t311 * t24 * t305 * t304 -
        0.11125297238589E4 * t317 * t304 -
        0.11125297238589E4 * t317 * t305 * t6 * t5 +
        0.11125297238589E4 * t317 * t324 * t5 -
        0.11125297238589E4 * t317 * t324 * t30;
      double t332 = t315 * t305;
      double t341 = t307 * t24;
      double t351 = t14 * t194;
      double t352 = t351 * t228;
      double t356 = t189 * t6 * t5;
      double t374 =
        -0.11125297238589E4 * t310 * t308 * t24 * t332 * t304 -
        0.11125297238589E4 * t311 * t332 * t304 -
        0.11125297238589E4 * t308 * t341 * t332 * t304 -
        0.11125297238589E4 * t310 * t341 * t332 * t304 -
        0.1165148256740162E1 * t352 * t227 * t30 +
        0.2330296513480324E1 * t352 * t356 - 0.4030853005674001E2 * t7 +
        0.2330296513480324E1 * t352 * t191 * t6 * t5 +
        0.1165148256740162E1 * t352 * t227 * t5 +
        0.2330296513480324E1 * t351 * t193 * t191 * t356 +
        0.2330296513480324E1 * t351 * t192 * t356;
      double t390 = 0.1371740148509607E1 + t8 + t1;
      double t391 = 0.1091700181433142E1 + t8 + t1;
      double t392 = t391 * t390;
      double t394 = 0.7092992179024789 + t8 + t1;
      double t395 = 0.2907007820975211 + t8 + t1;
      double t397 = -0.917001814331423E-1 + t8 + t1;
      double t398 = -0.3717401485096066 + t8 + t1;
      double t399 = t398 * t397;
      double t400 = t399 * t395 * t394;
      double t411 = t390 * t6 * t5;
      double t414 = t394 * t391;
      double t418 =
        0.1039836289900873E2 * t79 +
        0.2330296513480324E1 * t14 * t193 * t192 * t356 -
        0.3354340220928293E2 * t14 * t5 + 0.3354340220928293E2 * t14 * t30 -
        0.1039836289900873E2 * t76 + 0.3684522303073546E1 * t196 * t356 +
        0.6473987807847046 * t400 * t392 * t30 -
        0.1294797561569409E1 * t400 * t391 * t6 * t5 -
        0.6473987807847046 * t400 * t392 * t5 -
        0.1294797561569409E1 * t400 * t411 -
        0.1294797561569409E1 * t399 * t414 * t411;
      double t425 = t89 * t6;
      double t433 = t88 * t6;
      double t434 = t433 * t5;
      double t436 = t95 * t6;
      double t437 = t436 * t5;
      double t439 = t136 * t30;
      double t441 = t14 * t67;
      double t446 = t67 * t6;
      double t447 = t446 * t5;
      double t449 =
        -0.1294797561569409E1 * t399 * t395 * t391 * t411 -
        0.2057957043882974E3 * t143 * t30 - 0.2057957043882974E3 * t425 * t5 -
        0.2057957043882974E3 * t91 * t6 * t5 + 0.2057957043882974E3 * t143 * t5 +
        0.1640452312932599E2 * t434 + 0.1640452312932599E2 * t437 -
        0.8202261564662995E1 * t439 + 0.1847869343969603E2 * t441 * t5 -
        0.1847869343969603E2 * t441 * t30 + 0.5843475945345073E2 * t447;
      double t456 = t260 * t23 * t18;
      double t462 = t136 * t5;
      double t464 = t9 * t6;
      double t465 = t464 * t5;
      double t472 = t10 * t6;
      double t476 = -t3 - t22 + 0.3717401485096066;
      double t477 = -t3 - t22 + 0.917001814331423E-1;
      double t478 = t477 * t476;
      double t480 = -t3 - t22 - 0.2907007820975211;
      double t481 = -t3 - t22 - 0.7092992179024789;
      double t483 = -t3 - t22 - 0.1091700181433142E1;
      double t484 = -t3 - t22 - 0.1371740148509607E1;
      double t485 = t484 * t483;
      double t486 = t485 * t481 * t480;
      double t499 =
        -0.1168695189069015E3 * t14 * t6 * t5 -
        0.5228822893162457E2 * t456 * t11 * t5 +
        0.5228822893162457E2 * t456 * t11 * t30 + 0.8202261564662995E1 * t462 -
        0.1045764578632491E3 * t260 * t23 * t10 * t465 -
        0.1045764578632491E3 * t456 * t465 -
        0.1045764578632491E3 * t456 * t472 * t5 -
        0.8347107507735107E3 * t486 * t478 * t30 -
        0.8347107507735107E3 * t486 * t477 * t6 * t5 +
        0.8347107507735107E3 * t486 * t478 * t5 +
        0.5228822893162457E2 * t25 * t40 * t465;
      double t515 = t143 * t167;
      double t526 = t111 * t6 * t5;
      double t541 =
        0.5228822893162457E2 * t275 * t40 * t465 +
        0.5228822893162457E2 * t260 * t40 * t465 +
        0.294790221278654E4 * t156 * t160 * t6 * t5 -
        0.1294797561569409E1 * t398 * t395 * t414 * t411 -
        0.1449065949426045E2 * t515 * t113 * t5 +
        0.1449065949426045E2 * t515 * t113 * t30 -
        0.1294797561569409E1 * t397 * t395 * t414 * t411 -
        0.289813189885209E2 * t515 * t526 -
        0.289813189885209E2 * t515 * t112 * t6 * t5 +
        0.1449065949426045E2 * t89 * t118 * t126 * t526 +
        0.1449065949426045E2 * t91 * t118 * t126 * t526;
      double t552 = t476 * t6 * t5;
      double t559 = t160 * t67;
      double t563 = t480 * t477;
      double t582 =
        -0.289813189885209E2 * t143 * t126 * t526 -
        0.289813189885209E2 * t143 * t118 * t112 * t526 -
        0.8347107507735107E3 * t485 * t481 * t477 * t552 -
        0.8347107507735107E3 * t486 * t552 +
        0.4661042655927904E3 * t156 * t559 * t30 -
        0.8347107507735107E3 * t483 * t481 * t563 * t552 -
        0.8347107507735107E3 * t484 * t481 * t563 * t552 -
        0.8347107507735107E3 * t485 * t563 * t552 +
        0.4661042655927904E3 * t207 * t447 + 0.4661042655927904E3 * t156 * t447 -
        0.4661042655927904E3 * t156 * t559 * t5;
      double t607 = t31 * t88;
      double t611 =
        0.4661042655927904E3 * t200 * t447 + 0.4661042655927904E3 * t224 * t447 +
        0.4661042655927904E3 * t220 * t447 +
        0.5985439053268237E3 * t91 * t425 * t5 -
        0.1460487583524379E3 * t88 * t436 * t5 +
        0.5979402883100954E2 * t148 * t439 +
        0.5979402883100954E2 * t147 * t33 * t88 * t437 -
        0.1195880576620191E3 * t148 * t437 - 0.1195880576620191E3 * t148 * t434 -
        0.5979402883100954E2 * t148 * t462 +
        0.5979402883100954E2 * t35 * t607 * t437;
      double t619 = t14 * t136;
      double t636 = t23 * t6;
      double t644 =
        0.5979402883100954E2 * t243 * t607 * t437 +
        0.5979402883100954E2 * t147 * t607 * t437 -
        0.4618467258332445E2 * t619 * t5 + 0.4618467258332445E2 * t619 * t30 -
        0.923693451666489E2 * t14 * t436 * t5 -
        0.923693451666489E2 * t14 * t433 * t5 + 0.2213467294780192E1 * t19 * t30 -
        0.5201861283573172E3 * t26 * t24 * t6 * t5 -
        0.5201861283573172E3 * t26 * t636 * t5 -
        0.5201861283573172E3 * t27 * t30 + 0.5201861283573172E3 * t27 * t5;
      double t660 = t91 * t89 * t67;
      double t673 =
        -0.5201861283573172E3 * t24 * t636 * t5 -
        0.2213467294780192E1 * t19 * t5 - 0.4426934589560384E1 * t10 * t464 * t5 -
        0.4426934589560384E1 * t18 * t464 * t5 -
        0.4426934589560384E1 * t18 * t472 * t5 +
        0.946381010222476E2 * t660 * t30 - 0.946381010222476E2 * t660 * t5 +
        0.946381010222476E2 * t91 * t446 * t5 +
        0.946381010222476E2 * t89 * t446 * t5 - 0.1986085925560006E2 * t30 +
        0.1986085925560006E2 * t5;
      double t677 =
        t64 + t109 + t159 + t203 + t247 + t282 + t331 + t374 + t418 + t449 +
        t499 + t541 + t582 + t611 + t644 + t673;
      return t677;
    }

    double
    ortho2_f59y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * x;
      double t7 = 0.1154653670707977E1 + t6 + t1;
      double t8 = 0.5 + t6 + t1;
      double t9 = t8 * t7;
      double t11 = -0.1546536707079771 + t6 + t1;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 + 0.1546536707079771;
      double t15 = -t3 - t12 - 0.5;
      double t16 = -t3 - t12 - 0.1154653670707977E1;
      double t17 = t16 * t15;
      double t18 = t17 * t13 * t11;
      double t21 = -t3 - t1;
      double t22 = t21 * t2;
      double t23 = t8 * t4;
      double t27 = t7 * t4;
      double t28 = t27 * t22;
      double t31 = t4 * t21;
      double t32 = -t3 - t12 + 0.3302238962785669;
      double t34 = -t3 - t12 - 0.3115120652928579E-1;
      double t36 = -t3 - t12 - 0.9688487934707142;
      double t37 = -t3 - t12 - 0.1330223896278567E1;
      double t38 = t37 * t36;
      double t39 = t38 * t15 * t34;
      double t42 = 0.1330223896278567E1 + t6 + t1;
      double t44 = t42 * t4 * t22;
      double t45 = 0.9688487934707142 + t6 + t1;
      double t46 = t8 * t45;
      double t47 = -0.3302238962785669 + t6 + t1;
      double t50 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t51 = t50 * t47;
      double t55 = 0.3115120652928579E-1 + t6 + t1;
      double t60 = -t3 - t12 - 0.5278640450004206E-1;
      double t61 = -t3 - t12 - 0.9472135954999579;
      double t62 = t61 * t60;
      double t65 = t47 * t55;
      double t66 = t65 * t46;
      double t69 = t15 * t13;
      double t70 = t16 * t69;
      double t73 = 0.9472135954999579 + t6 + t1;
      double t74 = t73 * t4;
      double t80 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t81 = t80 * t4;
      double t82 = t81 * t22;
      double t84 = t38 * t15 * t32;
      double t87 =
        0.5228822893162457E2 * t18 * t9 * t5 -
        0.5228822893162457E2 * t18 * t23 * t22 -
        0.5228822893162457E2 * t18 * t28 +
        0.1011213042155722E4 * t39 * t32 * t31 +
        0.1165148256740162E1 * t51 * t46 * t44 +
        0.1165148256740162E1 * t50 * t55 * t46 * t44 -
        0.2057957043882974E3 * t62 * t5 + 0.7369044606147092E1 * t66 * t44 +
        0.5201861283573172E3 * t70 * t31 -
        0.4618467258332445E2 * t50 * t74 * t22 + 0.9322085311855808E3 * t84 * t82;
      double t88 = t34 * t32;
      double t89 = t38 * t88;
      double t96 = t4 * t22;
      double t97 = -t3 - t12 + 0.2650553239294647;
      double t98 = t97 * t80;
      double t99 = -t3 - t12 - 0.7852315164806451;
      double t100 = -t3 - t12 - 0.1265055323929465E1;
      double t101 = t100 * t99;
      double t106 = t36 * t15 * t88;
      double t109 = t7 * t31;
      double t110 = t11 * t8;
      double t111 = t62 * t110;
      double t115 = -t3 - t12 - 0.2147684835193549;
      double t116 = t99 * t115;
      double t117 = t100 * t116;
      double t127 = t16 * t13;
      double t131 = t80 * t31;
      double t134 = -t3 - t12 + 0.3717401485096066;
      double t136 = t134 * t4 * t22;
      double t137 = -t3 - t12 + 0.917001814331423E-1;
      double t138 = -t3 - t12 - 0.7092992179024789;
      double t140 = -t3 - t12 - 0.1091700181433142E1;
      double t141 = -t3 - t12 - 0.1371740148509607E1;
      double t142 = t141 * t140;
      double t146 =
        0.9322085311855808E3 * t89 * t82 -
        0.5228822893162457E2 * t17 * t13 * t8 * t28 +
        0.3680091408067806E3 * t101 * t98 * t96 -
        0.2022426084311444E4 * t106 * t96 - 0.1366810782197325E2 * t111 * t109 -
        0.6455230720453683E3 * t117 * t97 * t5 +
        0.147395110639327E4 * t39 * t32 * t4 * t22 +
        0.1045764578632491E3 * t17 * t110 * t28 +
        0.1045764578632491E3 * t127 * t110 * t28 -
        0.2640845960270727E3 * t70 * t131 -
        0.1669421501547021E4 * t142 * t138 * t137 * t136;
      double t148 = -t3 - t12 - 0.2907007820975211;
      double t149 = t148 * t137;
      double t153 = t100 * t115;
      double t166 = 0.7852315164806451 + t6 + t1;
      double t167 = 0.2147684835193549 + t6 + t1;
      double t168 = t167 * t166;
      double t169 = -0.2650553239294647 + t6 + t1;
      double t170 = t50 * t169;
      double t171 = t170 * t168;
      double t174 = 0.1265055323929465E1 + t6 + t1;
      double t175 = t167 * t174;
      double t179 = t166 * t174;
      double t195 =
        -0.1669421501547021E4 * t142 * t149 * t136 +
        0.3680091408067806E3 * t153 * t98 * t96 -
        0.1291046144090737E4 * t117 * t96 -
        0.1291046144090737E4 * t100 * t99 * t97 * t96 +
        0.3680091408067806E3 * t116 * t98 * t96 +
        0.1709835879221569E1 * t171 * t96 +
        0.1709835879221569E1 * t170 * t175 * t96 +
        0.1709835879221569E1 * t170 * t179 * t96 +
        0.1709835879221569E1 * t50 * t167 * t179 * t96 -
        0.1669421501547021E4 * t141 * t138 * t149 * t136 -
        0.1669421501547021E4 * t140 * t138 * t149 * t136;
      double t196 = 0.5278640450004206E-1 + t6 + t1;
      double t197 = t196 * t73;
      double t198 = t50 * t197;
      double t201 = t115 * t97;
      double t202 = t101 * t201;
      double t205 = t174 * t5;
      double t206 = t169 * t168;
      double t223 = t137 * t134;
      double t226 = t142 * t138 * t148;
      double t231 =
        0.4618467258332445E2 * t198 * t5 - 0.1840045704033903E3 * t202 * t131 +
        0.3183644212653517E1 * t206 * t205 - 0.3183644212653517E1 * t206 * t96 -
        0.2022426084311444E4 * t39 * t96 - 0.2022426084311444E4 * t84 * t96 +
        0.1045764578632491E3 * t69 * t110 * t28 +
        0.9041893030232619E-1 * t66 * t42 * t31 -
        0.5201861283573172E3 * t70 * t5 +
        0.8347107507735107E3 * t226 * t223 * t31 +
        0.8351088184110033E3 * t70 * t96;
      double t234 = t45 * t42;
      double t236 = t55 * t8;
      double t237 = t51 * t236;
      double t240 = t169 * t167;
      double t250 = t197 * t5;
      double t253 = t196 * t4;
      double t254 = t253 * t22;
      double t259 = t60 * t4;
      double t269 = t73 * t31;
      double t271 = t61 * t60 * t196;
      double t274 =
        0.1165148256740162E1 * t237 * t234 * t31 +
        0.1081395160683345E2 * t240 * t179 * t96 -
        0.2920975167048759E3 * t196 * t74 * t22 -
        0.5228822893162457E2 * t18 * t9 * t31 +
        0.5979402883100954E2 * t202 * t250 - 0.5979402883100954E2 * t202 * t254 +
        0.581873542355521E3 * t202 * t96 +
        0.2992719526634118E3 * t61 * t259 * t22 -
        0.2213467294780192E1 * t11 * t23 * t22 -
        0.2213467294780192E1 * t11 * t27 * t22 +
        0.3109772735774173E2 * t271 * t269;
      double t275 = t13 * t4;
      double t280 = 0.1371740148509607E1 + t6 + t1;
      double t281 = 0.1091700181433142E1 + t6 + t1;
      double t282 = t281 * t280;
      double t284 = 0.7092992179024789 + t6 + t1;
      double t285 = 0.2907007820975211 + t6 + t1;
      double t287 = -0.917001814331423E-1 + t6 + t1;
      double t288 = -0.3717401485096066 + t6 + t1;
      double t289 = t288 * t287;
      double t290 = t289 * t285 * t284;
      double t296 = t74 * t22;
      double t300 = t62 * t240;
      double t303 = t73 * t5;
      double t305 = t17 * t13 * t196;
      double t319 =
        -0.1040372256714634E4 * t15 * t275 * t22 - 0.1039836289900873E2 * t131 +
        0.6473987807847046 * t290 * t282 * t5 -
        0.8347107507735107E3 * t226 * t223 * t5 -
        0.5979402883100954E2 * t202 * t296 -
        0.1449065949426045E2 * t300 * t179 * t31 -
        0.4947070248095459E3 * t305 * t303 + 0.4947070248095459E3 * t305 * t96 -
        0.2213467294780192E1 * t8 * t27 * t22 -
        0.3109772735774173E2 * t271 * t303 -
        0.6473987807847046 * t290 * t281 * t4 * t22;
      double t322 = t280 * t4 * t22;
      double t332 = t11 * t9;
      double t335 = t50 * t110;
      double t338 = t11 * t7;
      double t347 = t7 * t5;
      double t356 =
        -0.6473987807847046 * t290 * t322 -
        0.1165148256740162E1 * t237 * t234 * t5 +
        0.1165148256740162E1 * t237 * t45 * t4 * t22 +
        0.2213467294780192E1 * t332 * t5 - 0.5355544023028608E1 * t335 * t96 -
        0.5355544023028608E1 * t50 * t338 * t96 -
        0.4618467258332445E2 * t198 * t31 +
        0.1449065949426045E2 * t300 * t179 * t5 +
        0.1366810782197325E2 * t111 * t347 +
        0.6455230720453683E3 * t117 * t97 * t31 -
        0.1291046144090737E4 * t100 * t201 * t96;
      double t375 = t13 * t80;
      double t387 =
        -0.1291046144090737E4 * t99 * t201 * t96 +
        0.1165148256740162E1 * t237 * t44 + 0.3109772735774173E2 * t271 * t96 +
        0.3109772735774173E2 * t61 * t60 * t73 * t96 -
        0.6219545471548346E2 * t61 * t197 * t96 -
        0.3387143444414331E2 * t332 * t96 - 0.2022426084311444E4 * t89 * t96 +
        0.5281691920541455E3 * t15 * t375 * t96 -
        0.3183644212653517E1 * t169 * t175 * t96 -
        0.3183644212653517E1 * t169 * t179 * t96 -
        0.1366810782197325E2 * t111 * t96;
      double t395 = -t3 - t12 + 0.3997579954114602;
      double t396 = -t3 - t12 + 0.1771862795107378;
      double t397 = t396 * t395;
      double t399 = -t3 - t12 - 0.1368825361738218;
      double t401 = -t3 - t12 - 0.8631174638261782;
      double t402 = -t3 - t12 - 0.1177186279510738E1;
      double t404 = -t3 - t12 - 0.139975799541146E1;
      double t405 = t404 * t402 * t401;
      double t406 = t405 * t15 * t399;
      double t424 = t61 * t60 * t80;
      double t433 =
        -0.1366810782197325E2 * t62 * t338 * t96 + 0.1986085925560006E2 * t31 -
        0.11125297238589E4 * t406 * t397 * t5 +
        0.1165148256740162E1 * t51 * t55 * t45 * t44 +
        0.1195880576620191E3 * t101 * t115 * t196 * t296 -
        0.1366810782197325E2 * t62 * t9 * t96 - 0.1986085925560006E2 * t5 -
        0.2213467294780192E1 * t332 * t31 + 0.946381010222476E2 * t424 * t5 -
        0.5355544023028608E1 * t50 * t9 * t96 -
        0.9041893030232619E-1 * t66 * t42 * t5;
      double t447 = t395 * t4 * t22;
      double t454 = t284 * t281;
      double t466 = t97 * t196;
      double t476 =
        -0.6473987807847046 * t289 * t285 * t281 * t322 -
        0.22250594477178E4 * t406 * t396 * t4 * t22 +
        0.4947070248095459E3 * t17 * t13 * t73 * t96 -
        0.22250594477178E4 * t406 * t447 -
        0.22250594477178E4 * t405 * t15 * t396 * t447 -
        0.6473987807847046 * t288 * t285 * t454 * t322 +
        0.1892762020444952E3 * t61 * t81 * t22 -
        0.1669421501547021E4 * t226 * t137 * t4 * t22 +
        0.1195880576620191E3 * t101 * t466 * t296 +
        0.1195880576620191E3 * t153 * t466 * t296 +
        0.1195880576620191E3 * t116 * t466 * t296;
      double t487 = t174 * t31;
      double t501 = t174 * t4 * t22;
      double t504 = t80 * t5;
      double t514 =
        -0.6473987807847046 * t287 * t285 * t454 * t322 +
        0.9041893030232619E-1 * t66 * t96 -
        0.3183644212653517E1 * t167 * t179 * t96 +
        0.1709835879221569E1 * t171 * t487 -
        0.9894140496190919E3 * t17 * t197 * t96 -
        0.9894140496190919E3 * t127 * t197 * t96 -
        0.1449065949426045E2 * t300 * t166 * t4 * t22 -
        0.1449065949426045E2 * t300 * t501 + 0.1840045704033903E3 * t202 * t504 +
        0.9041893030232619E-1 * t65 * t8 * t42 * t96 -
        0.9894140496190919E3 * t69 * t197 * t96;
      double t539 = t197 * t31;
      double t546 =
        0.1892762020444952E3 * t60 * t81 * t22 +
        0.11125297238589E4 * t406 * t397 * t31 +
        0.273362156439465E2 * t61 * t11 * t9 * t96 +
        0.273362156439465E2 * t60 * t11 * t9 * t96 +
        0.5355544023028608E1 * t335 * t347 +
        0.9041893030232619E-1 * t65 * t234 * t96 -
        0.946381010222476E2 * t424 * t31 - 0.1011213042155722E4 * t39 * t32 * t5 -
        0.5979402883100954E2 * t202 * t539 + 0.1039836289900873E2 * t504 -
        0.6473987807847046 * t289 * t454 * t322;
      double t550 = t37 * t15 * t88;
      double t553 = t32 * t80;
      double t569 = t399 * t396;
      double t578 = t401 * t15;
      double t589 =
        -0.2022426084311444E4 * t550 * t96 +
        0.4661042655927904E3 * t39 * t553 * t5 -
        0.3183644212653517E1 * t206 * t487 -
        0.1449065949426045E2 * t62 * t169 * t166 * t501 -
        0.1449065949426045E2 * t62 * t168 * t501 -
        0.4661042655927904E3 * t39 * t553 * t31 -
        0.22250594477178E4 * t405 * t569 * t447 -
        0.22250594477178E4 * t404 * t402 * t15 * t569 * t447 -
        0.22250594477178E4 * t404 * t578 * t569 * t447 -
        0.4618467258332445E2 * t50 * t253 * t22 -
        0.6219545471548346E2 * t60 * t197 * t96;
      double t620 =
        0.9041893030232619E-1 * t47 * t8 * t234 * t96 -
        0.3354340220928293E2 * t50 * t31 +
        0.9041893030232619E-1 * t236 * t234 * t96 -
        0.1040372256714634E4 * t16 * t15 * t4 * t22 -
        0.1040372256714634E4 * t16 * t275 * t22 -
        0.5355544023028608E1 * t335 * t109 + 0.2640845960270727E3 * t70 * t504 -
        0.1792645922070571E3 * t96 + 0.3354340220928293E2 * t50 * t5 -
        0.22250594477178E4 * t402 * t578 * t569 * t447 -
        0.4115914087765948E3 * t61 * t4 * t22;
      double t630 = t50 * t80;
      double t646 =
        -0.4115914087765948E3 * t259 * t22 + 0.2057957043882974E3 * t62 * t31 -
        0.5843475945345073E2 * t50 * t4 * t22 + 0.1168695189069015E3 * t82 -
        0.1847869343969603E2 * t630 * t5 + 0.1847869343969603E2 * t630 * t31 -
        0.8202261564662995E1 * t250 + 0.8202261564662995E1 * t296 +
        0.8202261564662995E1 * t254 +
        0.289813189885209E2 * t61 * t169 * t168 * t501 +
        0.289813189885209E2 * t60 * t169 * t168 * t501;
      double t674 =
        0.5281691920541455E3 * t16 * t15 * t80 * t96 +
        0.5281691920541455E3 * t16 * t375 * t96 +
        0.4947070248095459E3 * t305 * t269 +
        0.3680091408067806E3 * t101 * t115 * t80 * t96 -
        0.1669421501547021E4 * t226 * t136 -
        0.6473987807847046 * t290 * t282 * t31 -
        0.1709835879221569E1 * t171 * t205 + 0.9322085311855808E3 * t39 * t82 +
        0.8202261564662995E1 * t539 + 0.9322085311855808E3 * t550 * t82 +
        0.9322085311855808E3 * t106 * t82;
      double t678 =
        t87 + t146 + t195 + t231 + t274 + t319 + t356 + t387 + t433 + t476 +
        t514 + t546 + t589 + t620 + t646 + t674;
      return t678;
    }

    // * f60 *********************************************************************

    double
    ortho2_f60 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t10 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 = 0.1E1 * y;
      double t12 = -t3 - t11 - 0.5278640450004206E-1;
      double t14 = -t3 - t11 - 0.9472135954999579;
      double t18 = -t3 - t11 + 0.1546536707079771;
      double t19 = -t3 - t11 - 0.5;
      double t21 = -t3 - t11 - 0.1154653670707977E1;
      double t22 = t21 * t19 * t18;
      double t25 = 0.1E1 * x;
      double t26 = 0.1154653670707977E1 + t25 + t1;
      double t27 = 0.5 + t25 + t1;
      double t28 = t27 * t26;
      double t29 = -0.1546536707079771 + t25 + t1;
      double t33 = 0.9472135954999579 + t25 + t1;
      double t34 = 0.5278640450004206E-1 + t25 + t1;
      double t35 = t34 * t33;
      double t38 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t46 = t33 * t6;
      double t50 = t10 * t6;
      double t70 = -t3 - t11 + 0.3717401485096066;
      double t81 =
        (-t3 - t11 - 0.1371740148509607E1) * (-t3 - t11 -
                0.1091700181433142E1) * (-t3 - t11 -
                       0.7092992179024789)
        * (-t3 - t11 - 0.2907007820975211) * (-t3 - t11 + 0.917001814331423E-1);
      double t84 = -t3 - t11 + 0.2650553239294647;
      double t85 = -t3 - t11 - 0.2147684835193549;
      double t87 = -t3 - t11 - 0.7852315164806451;
      double t88 = -t3 - t11 - 0.1265055323929465E1;
      double t90 = t88 * t87 * t85 * t84;
      double t94 = t21 * t19;
      double t98 = t14 * t12;
      double t102 = 0.1265055323929465E1 + t25 + t1;
      double t103 = 0.7852315164806451 + t25 + t1;
      double t105 = 0.2147684835193549 + t25 + t1;
      double t106 = -0.2650553239294647 + t25 + t1;
      double MapleGenVar1 =
        -0.1251478649981627E4 * t14 * t12 * t10 * t7 +
        0.2960694596759533E4 * t22 * t7 + 0.2496294901559572E3 * t29 * t28 * t7 -
        0.3344643134161015E3 * t38 * t35 * t7 +
        0.1889846853861325E4 * t14 * t12 * t6 * t5 +
        0.4225517595916278E3 * t34 * t46 * t5 +
        0.1375756782543202E3 * t38 * t50 * t5;
      double t115 =
        MapleGenVar1 + 0.1907241707285146E3 * t7 + 0.4299672921430364E4 * (-t3 -
                           t11 -
                           0.139975799541146E1)
        * (-t3 - t11 - 0.1177186279510738E1) * (-t3 - t11 -
                  0.8631174638261782) * t19 * (-t3 -
                       t11 -
                       0.1368825361738218)
        * (-t3 - t11 + 0.1771862795107378) * (-t3 - t11 +
                0.3997579954114602) * t6 * t5 -
        0.2581740783477578E4 * t81 * t70 * t50 * t5 +
        0.5215417796841247E4 * t90 * t7 -
        0.1911546075590841E4 * t94 * t18 * t10 * t7 +
        0.3652744885876795E4 * t98 * t35 * t7 +
        0.1013349632260505E2 * t106 * t105 * t103 * t102 * t7 +
        0.1190984908687336E1 * t38 * t29 * t28 * t7;
      double t118 = (0.1330223896278567E1 + t25 + t1) * t6 * t5;
      double t120 = t27 * (0.9688487934707142 + t25 + t1);
      double t123 =
        (-0.3302238962785669 + t25 + t1) * (0.3115120652928579E-1 + t25 + t1);
      double t128 = t46 * t5;
      double t154 = t26 * t6 * t5;
      double t155 = t29 * t27;
      double t159 = t50 * t5;
      double t160 = -t3 - t11 + 0.3302238962785669;
      double t161 = -t3 - t11 - 0.3115120652928579E-1;
      double t163 = -t3 - t11 - 0.9688487934707142;
      double t165 = -t3 - t11 - 0.1330223896278567E1;
      double t171 = t102 * t6 * t5;
      double t172 = t105 * t103;
      MapleGenVar1 = 0.1055352582198434E2 * t38 * t123 * t120 * t118 +
        0.7666365967377942E4 * t88 * t87 * t85 * t84 * t34 * t128 -
        0.2895696814322679E1 * (-0.3717401485096066 + t25 +
              t1) * (-0.917001814331423E-1 + t25 +
               t1) * (0.2907007820975211 + t25 +
                t1) * (0.7092992179024789 + t25 +
                 t1) * (0.1091700181433142E1 +
                  t25 +
                  t1) *
        (0.1371740148509607E1 + t25 + t1) * t6 * t5 +
        0.4658972010987756E4 * t81 * t70 * t6 * t5 -
        0.2778496433973706E3 * t22 * t155 * t154 -
        0.2736529982271208E4 * t165 * t163 * t19 * t161 * t160 * t159 -
        0.1656416798924622E3 * t14 * t12 * t106 * t172 * t171;
      double t205 =
        MapleGenVar1 + 0.2847622329768398E4 * t94 * t18 * t34 * t128 -
        0.280584537891879E4 * t90 * t159 -
        0.4436219744433057E1 * t123 * t120 * t118 +
        0.3862040588646039E1 * t38 * t106 * t172 * t171 +
        0.1606012365029544E4 * t98 * t155 * t154 +
        0.5204320085351685E4 * t165 * t163 * t19 * t161 * t160 * t6 * t5 -
        0.1383286906474909E3 * t159 - 0.209416434036252E3 * t38 * t6 * t5;
      double t206 = t115 + t205;
      return t206;
    }

    double
    ortho2_f60x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t11 = t10 * t5;
      double t13 = t6 * t2;
      double t14 = 0.1E1 * y;
      double t15 = -t3 - t14 + 0.2650553239294647;
      double t17 = -t3 - t14 - 0.2147684835193549;
      double t18 = -t3 - t14 - 0.7852315164806451;
      double t19 = t18 * t17;
      double t20 = -t3 - t14 - 0.1265055323929465E1;
      double t21 = t20 * t19;
      double t24 = t6 * t5;
      double t25 = -t3 - t14 + 0.1546536707079771;
      double t26 = -t3 - t14 - 0.5;
      double t27 = t26 * t25;
      double t28 = -t3 - t14 - 0.1154653670707977E1;
      double t29 = t28 * t27;
      double t32 = 0.1E1 * x;
      double t33 = 0.1154653670707977E1 + t32 + t1;
      double t34 = 0.5 + t32 + t1;
      double t35 = t34 * t33;
      double t38 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t42 = -0.1546536707079771 + t32 + t1;
      double t43 = t42 * t35;
      double t49 = t42 * t33;
      double t53 = t42 * t34;
      double t54 = t38 * t53;
      double t57 = t17 * t15;
      double t67 =
        0.2175262469630748E3 * t11 - 0.2607708898420623E4 * t21 * t15 * t13 +
        0.6044839451223451E4 * t29 * t24 +
        0.1190984908687336E1 * t38 * t35 * t24 + 0.188311248516982E1 * t43 * t24 +
        0.2607708898420623E4 * t21 * t15 * t5 +
        0.1190984908687336E1 * t38 * t49 * t24 +
        0.1190984908687336E1 * t54 * t24 -
        0.2607708898420623E4 * t20 * t57 * t24 -
        0.2607708898420623E4 * t20 * t18 * t15 * t24 -
        0.2607708898420623E4 * t21 * t24;
      double t71 = 0.9472135954999579 + t32 + t1;
      double t72 = 0.5278640450004206E-1 + t32 + t1;
      double t73 = t72 * t71;
      double t74 = -t3 - t14 - 0.9472135954999579;
      double t78 = -t3 - t14 - 0.5278640450004206E-1;
      double t84 = t74 * t78 * t72;
      double t87 = t71 * t5;
      double t90 = t71 * t13;
      double t93 = t25 * t9;
      double t104 = t9 * t5;
      double t107 = t9 * t13;
      double t110 = t33 * t5;
      double t113 =
        -0.2607708898420623E4 * t18 * t57 * t24 -
        0.1826372442938397E4 * t74 * t73 * t24 +
        0.3652744885876795E4 * t74 * t78 * t71 * t24 +
        0.3652744885876795E4 * t84 * t24 + 0.1826372442938397E4 * t84 * t87 -
        0.1826372442938397E4 * t84 * t90 +
        0.9557730377954203E3 * t26 * t93 * t24 +
        0.9557730377954203E3 * t28 * t93 * t24 +
        0.9557730377954203E3 * t28 * t26 * t9 * t24 -
        0.9557730377954203E3 * t29 * t104 + 0.9557730377954203E3 * t29 * t107 +
        0.595492454343668 * t54 * t110;
      double t115 = t33 * t13;
      double t121 = 0.1265055323929465E1 + t32 + t1;
      double t122 = 0.7852315164806451 + t32 + t1;
      double t123 = t122 * t121;
      double t124 = -0.2650553239294647 + t32 + t1;
      double t128 = 0.2147684835193549 + t32 + t1;
      double t129 = t128 * t121;
      double t133 = t128 * t122;
      double t134 = t124 * t133;
      double t137 = t121 * t5;
      double t140 = t121 * t13;
      double t146 = t124 * t128;
      double t158 =
        -0.595492454343668 * t54 * t115 - 0.1826372442938397E4 * t78 * t73 * t24 +
        0.1013349632260505E2 * t124 * t123 * t24 +
        0.1013349632260505E2 * t124 * t129 * t24 +
        0.1013349632260505E2 * t134 * t24 + 0.5066748161302524E1 * t134 * t137 -
        0.5066748161302524E1 * t134 * t140 +
        0.1013349632260505E2 * t128 * t123 * t24 +
        0.6106422338069453E1 * t146 * t123 * t24 -
        0.8030061825147721E3 * t78 * t42 * t35 * t24 -
        0.8030061825147721E3 * t74 * t42 * t35 * t24;
      double t159 = t74 * t78;
      double t166 = t159 * t53;
      double t173 = t20 * t18;
      double t174 = t173 * t57;
      double t177 = -t3 - t14 + 0.3302238962785669;
      double t179 = -t3 - t14 - 0.3115120652928579E-1;
      double t181 = -t3 - t14 - 0.9688487934707142;
      double t182 = -t3 - t14 - 0.1330223896278567E1;
      double t183 = t182 * t181;
      double t184 = t183 * t26 * t179;
      double t191 = t183 * t26 * t177;
      double t196 = 0.1330223896278567E1 + t32 + t1;
      double t197 = 0.9688487934707142 + t32 + t1;
      double t198 = t197 * t196;
      double t199 = 0.3115120652928579E-1 + t32 + t1;
      double t200 = -0.3302238962785669 + t32 + t1;
      double t201 = t200 * t199;
      double t209 =
        0.1606012365029544E4 * t159 * t35 * t24 +
        0.1606012365029544E4 * t159 * t49 * t24 +
        0.1606012365029544E4 * t166 * t24 + 0.8030061825147721E3 * t166 * t110 -
        0.8030061825147721E3 * t166 * t115 + 0.887286215964157E4 * t174 * t24 +
        0.2602160042675842E4 * t184 * t177 * t5 -
        0.2602160042675842E4 * t184 * t177 * t13 -
        0.2602160042675842E4 * t191 * t24 - 0.2602160042675842E4 * t184 * t24 -
        0.4436219744433057E1 * t201 * t198 * t24 -
        0.4436219744433057E1 * t201 * t34 * t196 * t24;
      double t212 = t34 * t197;
      double t213 = t201 * t212;
      double t222 = t179 * t177;
      double t223 = t183 * t222;
      double t226 = t199 * t34;
      double t234 = t15 * t9;
      double t238 = t20 * t17;
      double t251 =
        -0.4436219744433057E1 * t213 * t24 -
        0.2218109872216528E1 * t213 * t196 * t5 +
        0.2218109872216528E1 * t213 * t196 * t13 -
        0.2602160042675842E4 * t223 * t24 -
        0.4436219744433057E1 * t226 * t198 * t24 -
        0.4436219744433057E1 * t200 * t34 * t198 * t24 +
        0.1402922689459395E4 * t19 * t234 * t24 +
        0.1402922689459395E4 * t238 * t234 * t24 +
        0.1402922689459395E4 * t173 * t234 * t24 +
        0.1402922689459395E4 * t173 * t17 * t9 * t24 -
        0.1402922689459395E4 * t174 * t104;
      double t255 = t181 * t26 * t222;
      double t259 = t182 * t26 * t222;
      double t262 = t28 * t26;
      double t271 = t262 * t25 * t72;
      double t281 = t28 * t25;
      double t285 = t38 * t124;
      double t286 = t285 * t133;
      double t291 =
        0.1402922689459395E4 * t174 * t107 - 0.2602160042675842E4 * t255 * t24 -
        0.2602160042675842E4 * t259 * t24 -
        0.1423811164884199E4 * t262 * t73 * t24 +
        0.2847622329768398E4 * t262 * t25 * t71 * t24 +
        0.2847622329768398E4 * t271 * t24 + 0.1423811164884199E4 * t271 * t87 -
        0.1423811164884199E4 * t271 * t90 -
        0.1423811164884199E4 * t27 * t73 * t24 -
        0.1423811164884199E4 * t281 * t73 * t24 -
        0.1931020294323019E1 * t286 * t140 + 0.3862040588646039E1 * t286 * t24;
      double t305 = -t3 - t14 + 0.3717401485096066;
      double t306 = t305 * t9;
      double t308 = -t3 - t14 + 0.917001814331423E-1;
      double t309 = -t3 - t14 - 0.2907007820975211;
      double t310 = t309 * t308;
      double t311 = -t3 - t14 - 0.7092992179024789;
      double t312 = -t3 - t14 - 0.1091700181433142E1;
      double t313 = t312 * t311;
      double t314 = -t3 - t14 - 0.1371740148509607E1;
      double t315 = t314 * t313;
      double t316 = t315 * t310;
      double t324 = t308 * t305;
      double t332 = t311 * t309;
      double t342 =
        0.1931020294323019E1 * t286 * t137 +
        0.3862040588646039E1 * t285 * t123 * t24 +
        0.3862040588646039E1 * t285 * t129 * t24 +
        0.3862040588646039E1 * t38 * t128 * t123 * t24 -
        0.1290870391738789E4 * t316 * t306 * t5 +
        0.1290870391738789E4 * t316 * t306 * t13 +
        0.1290870391738789E4 * t316 * t11 +
        0.1290870391738789E4 * t315 * t324 * t11 +
        0.1290870391738789E4 * t315 * t309 * t305 * t11 +
        0.1290870391738789E4 * t314 * t332 * t324 * t11 +
        0.1290870391738789E4 * t314 * t312 * t309 * t324 * t11;
      double t347 = -t3 - t14 + 0.3997579954114602;
      double t348 = -t3 - t14 + 0.1771862795107378;
      double t349 = t348 * t347;
      double t351 = -t3 - t14 - 0.1368825361738218;
      double t353 = -t3 - t14 - 0.8631174638261782;
      double t354 = -t3 - t14 - 0.1177186279510738E1;
      double t356 = -t3 - t14 - 0.139975799541146E1;
      double t357 = t356 * t354 * t353;
      double t358 = t357 * t26 * t351;
      double t365 = t305 * t6 * t5;
      double t373 = t347 * t6 * t5;
      double t374 = t351 * t348;
      double t389 = t353 * t26;
      double t398 = t78 * t6;
      double t402 =
        0.1290870391738789E4 * t312 * t332 * t324 * t11 +
        0.2149836460715182E4 * t358 * t349 * t5 -
        0.2149836460715182E4 * t358 * t349 * t13 +
        0.8164181203936754E4 * t316 * t365 -
        0.2149836460715182E4 * t358 * t348 * t6 * t5 -
        0.2149836460715182E4 * t356 * t354 * t26 * t374 * t373 -
        0.2149836460715182E4 * t357 * t374 * t373 -
        0.2149836460715182E4 * t357 * t26 * t348 * t373 -
        0.2149836460715182E4 * t358 * t373 -
        0.2149836460715182E4 * t354 * t389 * t374 * t373 -
        0.2149836460715182E4 * t356 * t389 * t374 * t373 +
        0.3957522977014582E4 * t74 * t398 * t5;
      double t406 = t71 * t6;
      double t413 = t196 * t6 * t5;
      double t415 = t38 * t200;
      double t419 = t415 * t226;
      double t441 =
        -0.5288345132196465E3 * t72 * t406 * t5 - 0.104708217018126E3 * t38 * t5 +
        0.1055352582198434E2 * t415 * t199 * t197 * t413 +
        0.1055352582198434E2 * t419 * t413 +
        0.1055352582198434E2 * t419 * t197 * t6 * t5 +
        0.5276762910992171E1 * t419 * t198 * t5 -
        0.5276762910992171E1 * t419 * t198 * t13 +
        0.1055352582198434E2 * t38 * t199 * t212 * t413 +
        0.1055352582198434E2 * t415 * t212 * t413 + 0.1063172726824206E3 * t24 +
        0.6916434532374543E2 * t107;
      double t445 = 0.1371740148509607E1 + t32 + t1;
      double t447 = t445 * t6 * t5;
      double t448 = 0.7092992179024789 + t32 + t1;
      double t449 = 0.2907007820975211 + t32 + t1;
      double t451 = -0.917001814331423E-1 + t32 + t1;
      double t452 = -0.3717401485096066 + t32 + t1;
      double t453 = t452 * t451;
      double t454 = t453 * t449 * t448;
      double t457 = 0.1091700181433142E1 + t32 + t1;
      double t462 = t457 * t445;
      double t473 = t448 * t457;
      double t479 = t262 * t25 * t42;
      double t482 = t33 * t6;
      double t483 = t482 * t5;
      double t490 = t34 * t6;
      double t494 =
        0.1668658947143561E2 * t213 * t413 - 0.6916434532374543E2 * t104 -
        0.2895696814322679E1 * t454 * t447 -
        0.2895696814322679E1 * t454 * t457 * t6 * t5 -
        0.144784840716134E1 * t454 * t462 * t5 +
        0.144784840716134E1 * t454 * t462 * t13 -
        0.2895696814322679E1 * t453 * t449 * t457 * t447 -
        0.2895696814322679E1 * t453 * t473 * t447 +
        0.1389248216986853E3 * t479 * t35 * t13 -
        0.2778496433973706E3 * t262 * t25 * t34 * t483 -
        0.2778496433973706E3 * t479 * t483 -
        0.2778496433973706E3 * t479 * t490 * t5;
      double t501 = t314 * t312;
      double t502 = t501 * t332;
      double t512 = t159 * t146;
      double t535 =
        -0.1389248216986853E3 * t479 * t35 * t5 -
        0.2329486005493878E4 * t502 * t308 * t6 * t5 +
        0.2329486005493878E4 * t502 * t324 * t5 -
        0.2329486005493878E4 * t502 * t324 * t13 +
        0.8282083994623109E2 * t512 * t123 * t13 +
        0.1389248216986853E3 * t27 * t53 * t483 +
        0.1389248216986853E3 * t281 * t53 * t483 +
        0.1389248216986853E3 * t262 * t53 * t483 +
        0.8653667629317213E4 * t184 * t177 * t6 * t5 -
        0.2895696814322679E1 * t452 * t449 * t473 * t447 -
        0.8282083994623109E2 * t512 * t123 * t5;
      double t537 = t121 * t6 * t5;
      double t569 = t177 * t9;
      double t580 =
        0.8282083994623109E2 * t78 * t124 * t133 * t537 +
        0.8282083994623109E2 * t74 * t124 * t133 * t537 -
        0.1656416798924622E3 * t159 * t133 * t537 -
        0.1656416798924622E3 * t159 * t124 * t122 * t537 -
        0.1656416798924622E3 * t512 * t537 -
        0.1656416798924622E3 * t512 * t122 * t6 * t5 -
        0.2895696814322679E1 * t451 * t449 * t473 * t447 -
        0.2329486005493878E4 * t502 * t365 -
        0.2329486005493878E4 * t501 * t311 * t308 * t365 +
        0.1368264991135604E4 * t184 * t569 * t13 -
        0.2329486005493878E4 * t313 * t310 * t365 -
        0.2329486005493878E4 * t314 * t311 * t310 * t365;
      double t588 = t38 * t73;
      double t604 = t406 * t5;
      double t611 =
        -0.2329486005493878E4 * t501 * t310 * t365 +
        0.1368264991135604E4 * t223 * t11 - 0.1672321567080507E3 * t588 * t5 +
        0.1672321567080507E3 * t588 * t13 + 0.1368264991135604E4 * t191 * t11 +
        0.1368264991135604E4 * t184 * t11 -
        0.1368264991135604E4 * t184 * t569 * t5 +
        0.1368264991135604E4 * t259 * t11 + 0.1368264991135604E4 * t255 * t11 -
        0.3833182983688971E4 * t173 * t17 * t72 * t604 +
        0.7666365967377942E4 * t174 * t604;
      double t612 = t72 * t6;
      double t613 = t612 * t5;
      double t616 = t73 * t5;
      double t619 = t73 * t13;
      double t622 = t15 * t72;
      double t647 =
        0.7666365967377942E4 * t174 * t613 + 0.3833182983688971E4 * t174 * t616 -
        0.3833182983688971E4 * t174 * t619 -
        0.3833182983688971E4 * t238 * t622 * t604 -
        0.3833182983688971E4 * t173 * t622 * t604 -
        0.3833182983688971E4 * t19 * t622 * t604 -
        0.3344643134161015E3 * t38 * t406 * t5 -
        0.3344643134161015E3 * t38 * t612 * t5 -
        0.9449234269306626E3 * t159 * t13 - 0.1248147450779786E3 * t43 * t13 +
        0.9449234269306626E3 * t159 * t5 - 0.4350524939261495E3 * t38 * t6 * t5;
      double t661 = t25 * t6;
      double t675 =
        -0.9449234269306626E3 * t398 * t5 - 0.9449234269306626E3 * t74 * t6 * t5 -
        0.2112758797958139E3 * t619 - 0.1480347298379766E4 * t28 * t26 * t6 * t5 +
        0.4225517595916278E3 * t613 + 0.4225517595916278E3 * t604 -
        0.1480347298379766E4 * t28 * t661 * t5 + 0.2112758797958139E3 * t616 -
        0.1480347298379766E4 * t26 * t661 * t5 + 0.1480347298379766E4 * t29 * t5 -
        0.1480347298379766E4 * t29 * t13 + 0.1248147450779786E3 * t43 * t5;
      double t686 = t74 * t78 * t9;
      double t697 = t38 * t9;
      double t706 =
        0.2496294901559572E3 * t42 * t482 * t5 +
        0.2496294901559572E3 * t42 * t490 * t5 +
        0.2496294901559572E3 * t34 * t482 * t5 +
        0.6257393249908136E3 * t686 * t13 - 0.6257393249908136E3 * t686 * t5 +
        0.6257393249908136E3 * t74 * t10 * t5 +
        0.6257393249908136E3 * t78 * t10 * t5 + 0.6878783912716011E2 * t697 * t5 +
        0.104708217018126E3 * t38 * t13 - 0.9536208536425731E2 * t13 +
        0.9536208536425731E2 * t5 - 0.6878783912716011E2 * t697 * t13;
      double t710 =
        t67 + t113 + t158 + t209 + t251 + t291 + t342 + t402 + t441 + t494 +
        t535 + t580 + t611 + t647 + t675 + t706;
      return t710;
    }

    double
    ortho2_f60y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.9472135954999579 + t8 + t1;
      double t10 = 0.5278640450004206E-1 + t8 + t1;
      double t11 = t10 * t9;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 - 0.9472135954999579;
      double t17 = -t3 - t12 + 0.1546536707079771;
      double t19 = -t3 - t12 - 0.5;
      double t20 = -t3 - t12 - 0.1154653670707977E1;
      double t21 = t20 * t19;
      double t25 = 0.1371740148509607E1 + t8 + t1;
      double t27 = t25 * t6 * t5;
      double t28 = 0.1091700181433142E1 + t8 + t1;
      double t29 = 0.7092992179024789 + t8 + t1;
      double t30 = t29 * t28;
      double t31 = -0.917001814331423E-1 + t8 + t1;
      double t32 = -0.3717401485096066 + t8 + t1;
      double t33 = t32 * t31;
      double t37 = 0.1154653670707977E1 + t8 + t1;
      double t38 = t37 * t6;
      double t39 = -0.1546536707079771 + t8 + t1;
      double t43 = 0.5 + t8 + t1;
      double t47 = t6 * t4;
      double t50 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t51 = t50 * t11;
      double t54 = -t3 - t12 + 0.2650553239294647;
      double t55 = -t3 - t12 - 0.7852315164806451;
      double t57 = -t3 - t12 - 0.1265055323929465E1;
      double t61 = -t3 - t12 - 0.2147684835193549;
      double t62 = t61 * t54;
      double t66 = t9 * t47;
      double t68 = t21 * t17 * t10;
      double t74 = t9 * t6;
      double t75 = t74 * t5;
      double t76 = t54 * t10;
      double t77 = t57 * t55;
      double t81 =
        -0.3652744885876795E4 * t13 * t11 * t7 +
        0.1423811164884199E4 * t21 * t17 * t9 * t7 -
        0.144784840716134E1 * t33 * t30 * t27 +
        0.1248147450779786E3 * t39 * t38 * t5 +
        0.1248147450779786E3 * t43 * t38 * t5 - 0.1672321567080507E3 * t51 * t47 -
        0.5215417796841247E4 * t57 * t55 * t54 * t7 -
        0.5215417796841247E4 * t57 * t62 * t7 + 0.1423811164884199E4 * t68 * t66 -
        0.5215417796841247E4 * t55 * t62 * t7 -
        0.7666365967377942E4 * t77 * t76 * t75;
      double t82 = -t3 - t12 + 0.3302238962785669;
      double t83 = -t3 - t12 - 0.3115120652928579E-1;
      double t84 = t83 * t82;
      double t85 = -t3 - t12 - 0.1330223896278567E1;
      double t87 = t85 * t19 * t84;
      double t90 = 0.1330223896278567E1 + t8 + t1;
      double t92 = t90 * t6 * t5;
      double t93 = 0.9688487934707142 + t8 + t1;
      double t94 = t43 * t93;
      double t95 = 0.3115120652928579E-1 + t8 + t1;
      double t104 = t17 * t6;
      double t108 = t43 * t37;
      double t112 = t6 * t2;
      double t114 = -0.3302238962785669 + t8 + t1;
      double t115 = t114 * t95;
      double t116 = t115 * t94;
      double t120 = -t3 - t12 - 0.9688487934707142;
      double t121 = t85 * t120;
      double t122 = t121 * t19 * t83;
      double t126 = t121 * t19 * t82;
      double t129 = t19 * t17;
      double t130 = t20 * t129;
      double t135 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t136 = t82 * t135;
      double t140 = t135 * t6;
      double t141 = t140 * t5;
      double t146 =
        -0.5204320085351685E4 * t87 * t7 +
        0.5276762910992171E1 * t50 * t95 * t94 * t92 -
        0.2960694596759533E4 * t20 * t19 * t6 * t5 -
        0.2960694596759533E4 * t20 * t104 * t5 +
        0.595492454343668 * t50 * t108 * t7 +
        0.2218109872216528E1 * t116 * t90 * t112 -
        0.5204320085351685E4 * t122 * t7 - 0.5204320085351685E4 * t126 * t7 -
        0.1480347298379766E4 * t130 * t112 +
        0.1368264991135604E4 * t122 * t136 * t112 +
        0.2736529982271208E4 * t122 * t141 + 0.104708217018126E3 * t50 * t112;
      double t150 = -t3 - t12 - 0.5278640450004206E-1;
      double t151 = t13 * t150;
      double t154 = t57 * t61;
      double t158 = t55 * t61;
      double t165 = -t3 - t12 + 0.3717401485096066;
      double t166 = t165 * t135;
      double t168 = -t3 - t12 + 0.917001814331423E-1;
      double t169 = -t3 - t12 - 0.2907007820975211;
      double t170 = t169 * t168;
      double t171 = -t3 - t12 - 0.7092992179024789;
      double t172 = -t3 - t12 - 0.1091700181433142E1;
      double t173 = t172 * t171;
      double t174 = -t3 - t12 - 0.1371740148509607E1;
      double t175 = t174 * t173;
      double t176 = t175 * t170;
      double t179 = 0.1265055323929465E1 + t8 + t1;
      double t180 = 0.7852315164806451 + t8 + t1;
      double t181 = t180 * t179;
      double t183 = 0.2147684835193549 + t8 + t1;
      double t184 = -0.2650553239294647 + t8 + t1;
      double t185 = t184 * t183;
      double t186 = t151 * t185;
      double t189 = t179 * t112;
      double t190 = t183 * t180;
      double t191 = t50 * t184;
      double t192 = t191 * t190;
      double t198 = t20 * t17;
      double t205 =
        -0.104708217018126E3 * t50 * t47 - 0.9449234269306626E3 * t151 * t112 -
        0.7666365967377942E4 * t154 * t76 * t75 -
        0.7666365967377942E4 * t158 * t76 * t75 -
        0.2847622329768398E4 * t21 * t11 * t7 +
        0.1290870391738789E4 * t176 * t166 * t112 +
        0.8282083994623109E2 * t186 * t181 * t112 -
        0.1931020294323019E1 * t192 * t189 -
        0.2960694596759533E4 * t19 * t104 * t5 -
        0.2847622329768398E4 * t198 * t11 * t7 -
        0.2847622329768398E4 * t129 * t11 * t7;
      double t215 = t93 * t90;
      double t226 = t179 * t6 * t5;
      double t235 = t54 * t135;
      double t241 = t168 * t165;
      double t243 = t171 * t169;
      double t244 = t174 * t172;
      double t245 = t244 * t243;
      double t249 =
        -0.3652744885876795E4 * t150 * t11 * t7 -
        0.2218109872216528E1 * t116 * t7 -
        0.2218109872216528E1 * t115 * t43 * t90 * t7 -
        0.2218109872216528E1 * t115 * t215 * t7 +
        0.2581740783477578E4 * t176 * t141 +
        0.2581740783477578E4 * t175 * t169 * t165 * t141 +
        0.1656416798924622E3 * t150 * t184 * t190 * t226 +
        0.280584537891879E4 * t77 * t61 * t135 * t7 +
        0.280584537891879E4 * t77 * t235 * t7 +
        0.2736529982271208E4 * t126 * t141 -
        0.2329486005493878E4 * t245 * t241 * t112 + 0.9536208536425731E2 * t47;
      double t255 = t120 * t19 * t84;
      double t262 = t10 * t6;
      double t266 = t95 * t43;
      double t281 = 0.2907007820975211 + t8 + t1;
      double t286 = t121 * t84;
      double t289 = t135 * t47;
      double t292 =
        0.2736529982271208E4 * t87 * t141 + 0.2736529982271208E4 * t255 * t141 -
        0.2218109872216528E1 * t114 * t43 * t215 * t7 -
        0.1672321567080507E3 * t50 * t262 * t5 -
        0.2218109872216528E1 * t266 * t215 * t7 +
        0.2581740783477578E4 * t175 * t241 * t141 +
        0.2581740783477578E4 * t174 * t172 * t169 * t241 * t141 -
        0.1672321567080507E3 * t50 * t74 * t5 -
        0.144784840716134E1 * t32 * t281 * t30 * t27 +
        0.2736529982271208E4 * t286 * t141 - 0.9557730377954203E3 * t130 * t289;
      double t293 = t77 * t62;
      double t311 = t183 * t179;
      double t328 =
        0.4436431079820785E4 * t293 * t7 -
        0.4658972010987756E4 * t245 * t168 * t6 * t5 +
        0.1931020294323019E1 * t192 * t7 +
        0.280584537891879E4 * t154 * t235 * t7 +
        0.280584537891879E4 * t158 * t235 * t7 -
        0.1368264991135604E4 * t122 * t136 * t47 +
        0.1931020294323019E1 * t191 * t311 * t7 +
        0.1931020294323019E1 * t191 * t181 * t7 -
        0.5204320085351685E4 * t255 * t7 +
        0.1931020294323019E1 * t50 * t183 * t181 * t7 +
        0.1672321567080507E3 * t51 * t112 + 0.3022419725611726E4 * t130 * t7;
      double t332 = t37 * t112;
      double t333 = t39 * t43;
      double t334 = t151 * t333;
      double t337 = t179 * t47;
      double t340 = t37 * t47;
      double t344 = t57 * t158;
      double t364 =
        0.1480347298379766E4 * t130 * t47 - 0.8030061825147721E3 * t334 * t332 +
        0.1931020294323019E1 * t192 * t337 + 0.8030061825147721E3 * t334 * t340 -
        0.2607708898420623E4 * t344 * t54 * t112 -
        0.2218109872216528E1 * t116 * t90 * t47 - 0.4435160469275132E3 * t7 -
        0.8282083994623109E2 * t186 * t181 * t47 -
        0.144784840716134E1 * t31 * t281 * t30 * t27 -
        0.8282083994623109E2 * t186 * t180 * t6 * t5 +
        0.8030061825147721E3 * t334 * t7;
      double t370 = t150 * t6;
      double t379 = t50 * t135;
      double t384 = t11 * t112;
      double t387 = t262 * t5;
      double t389 = t39 * t108;
      double t392 =
        -0.8282083994623109E2 * t186 * t226 -
        0.1889846853861325E4 * t13 * t6 * t5 - 0.1889846853861325E4 * t370 * t5 +
        0.9449234269306626E3 * t151 * t47 - 0.2175262469630748E3 * t50 * t6 * t5 +
        0.4350524939261495E3 * t141 - 0.6878783912716011E2 * t379 * t112 +
        0.6878783912716011E2 * t379 * t47 - 0.2112758797958139E3 * t384 +
        0.2112758797958139E3 * t75 + 0.2112758797958139E3 * t387 +
        0.376622497033964E1 * t389 * t7;
      double t398 = t184 * t190;
      double t412 = t165 * t6 * t5;
      double t420 = t39 * t37;
      double t425 = t50 * t114;
      double t426 = t425 * t266;
      double t432 =
        -0.1402922689459395E4 * t293 * t289 - 0.5066748161302524E1 * t398 * t189 +
        0.5066748161302524E1 * t398 * t7 +
        0.2581740783477578E4 * t174 * t243 * t241 * t141 +
        0.2581740783477578E4 * t172 * t243 * t241 * t141 +
        0.4082090601968377E4 * t176 * t412 + 0.3337317894287122E2 * t116 * t92 +
        0.2607708898420623E4 * t344 * t54 * t47 +
        0.8030061825147721E3 * t151 * t420 * t7 +
        0.5276762910992171E1 * t426 * t215 * t47 +
        0.1221284467613891E2 * t185 * t181 * t7;
      double t441 = t135 * t112;
      double t444 = t9 * t112;
      double t448 = t13 * t150 * t135;
      double t457 = t50 * t333;
      double t467 =
        -0.1057669026439293E4 * t74 * t10 * t5 +
        0.1248147450779786E3 * t389 * t47 +
        0.1978761488507291E4 * t13 * t370 * t5 + 0.6916434532374543E2 * t441 -
        0.6916434532374543E2 * t289 - 0.1423811164884199E4 * t68 * t444 +
        0.6257393249908136E3 * t448 * t112 +
        0.1251478649981627E4 * t13 * t140 * t5 +
        0.1251478649981627E4 * t150 * t140 * t5 -
        0.595492454343668 * t457 * t332 +
        0.8030061825147721E3 * t151 * t108 * t7 -
        0.1606012365029544E4 * t13 * t39 * t108 * t7;
      double t488 = t13 * t150 * t10;
      double t499 =
        0.2602160042675842E4 * t122 * t82 * t47 +
        0.2329486005493878E4 * t245 * t241 * t47 +
        0.5066748161302524E1 * t184 * t311 * t7 +
        0.1423811164884199E4 * t68 * t7 - 0.1248147450779786E3 * t389 * t112 +
        0.595492454343668 * t457 * t7 + 0.595492454343668 * t50 * t420 * t7 +
        0.1826372442938397E4 * t488 * t66 - 0.3833182983688971E4 * t293 * t384 +
        0.5066748161302524E1 * t398 * t337 +
        0.4326833814658606E4 * t122 * t82 * t6 * t5;
      double t500 = t38 * t5;
      double t516 = -t3 - t12 + 0.3997579954114602;
      double t518 = t516 * t6 * t5;
      double t519 = -t3 - t12 + 0.1771862795107378;
      double t521 = -t3 - t12 - 0.8631174638261782;
      double t522 = -t3 - t12 - 0.1177186279510738E1;
      double t524 = -t3 - t12 - 0.139975799541146E1;
      double t525 = t524 * t522 * t521;
      double t531 = -t3 - t12 - 0.1368825361738218;
      double t533 = t525 * t19 * t531;
      double t538 = t519 * t516;
      double t545 = t17 * t135;
      double t549 =
        0.2778496433973706E3 * t21 * t333 * t500 +
        0.1911546075590841E4 * t20 * t19 * t135 * t7 +
        0.595492454343668 * t457 * t340 + 0.9557730377954203E3 * t130 * t441 +
        0.3833182983688971E4 * t293 * t387 + 0.3833182983688971E4 * t293 * t75 -
        0.4299672921430364E4 * t525 * t19 * t519 * t518 -
        0.4299672921430364E4 * t533 * t519 * t6 * t5 -
        0.4299672921430364E4 * t533 * t518 -
        0.2149836460715182E4 * t533 * t538 * t112 -
        0.1290870391738789E4 * t176 * t166 * t47 +
        0.1911546075590841E4 * t20 * t545 * t7;
      double t566 = t21 * t17 * t39;
      double t569 = t43 * t6;
      double t580 = t531 * t519;
      double t581 = t521 * t19;
      double t591 =
        -0.4658972010987756E4 * t245 * t412 +
        0.2778496433973706E3 * t198 * t333 * t500 -
        0.6257393249908136E3 * t448 * t47 -
        0.2602160042675842E4 * t122 * t82 * t112 +
        0.1402922689459395E4 * t293 * t441 +
        0.1389248216986853E3 * t566 * t108 * t112 -
        0.1389248216986853E3 * t566 * t569 * t5 +
        0.5066748161302524E1 * t184 * t181 * t7 -
        0.1606012365029544E4 * t150 * t39 * t108 * t7 -
        0.4299672921430364E4 * t524 * t581 * t580 * t518 -
        0.4299672921430364E4 * t524 * t522 * t19 * t580 * t518;
      double t592 = t11 * t47;
      double t613 = t28 * t25;
      double t616 = t33 * t281 * t29;
      double t631 =
        0.3833182983688971E4 * t293 * t592 -
        0.4299672921430364E4 * t525 * t580 * t518 -
        0.1826372442938397E4 * t488 * t444 + 0.1826372442938397E4 * t488 * t7 +
        0.1826372442938397E4 * t13 * t150 * t9 * t7 +
        0.5066748161302524E1 * t183 * t181 * t7 -
        0.4299672921430364E4 * t522 * t581 * t580 * t518 +
        0.144784840716134E1 * t616 * t613 * t112 -
        0.5204320085351685E4 * t286 * t7 -
        0.4658972010987756E4 * t244 * t171 * t168 * t412 -
        0.4658972010987756E4 * t244 * t170 * t412 -
        0.5276762910992171E1 * t426 * t215 * t112;
      double t666 =
        0.5276762910992171E1 * t426 * t93 * t6 * t5 +
        0.5276762910992171E1 * t426 * t92 - 0.9536208536425731E2 * t112 -
        0.4658972010987756E4 * t174 * t171 * t170 * t412 +
        0.2112758797958139E3 * t592 - 0.1389248216986853E3 * t566 * t500 -
        0.1389248216986853E3 * t21 * t17 * t43 * t500 -
        0.4658972010987756E4 * t173 * t170 * t412 +
        0.1911546075590841E4 * t19 * t545 * t7 -
        0.144784840716134E1 * t616 * t28 * t6 * t5 -
        0.144784840716134E1 * t616 * t27 +
        0.2778496433973706E3 * t129 * t333 * t500;
      double t707 =
        -0.1389248216986853E3 * t566 * t108 * t47 -
        0.8282083994623109E2 * t151 * t184 * t180 * t226 -
        0.8282083994623109E2 * t151 * t190 * t226 -
        0.144784840716134E1 * t33 * t281 * t28 * t27 +
        0.1248147450779786E3 * t39 * t569 * t5 +
        0.2149836460715182E4 * t533 * t538 * t47 -
        0.144784840716134E1 * t616 * t613 * t47 +
        0.1656416798924622E3 * t13 * t184 * t190 * t226 -
        0.5215417796841247E4 * t344 * t7 -
        0.7666365967377942E4 * t77 * t61 * t10 * t75 +
        0.5276762910992171E1 * t425 * t95 * t93 * t92 +
        0.5276762910992171E1 * t425 * t94 * t92;
      double t711 =
        t81 + t146 + t205 + t249 + t292 + t328 + t364 + t392 + t432 + t467 +
        t499 + t549 + t591 + t631 + t666 + t707;
      return t711;
    }

    // * f61 *********************************************************************

    double
    ortho2_f61 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * x;
      double t8 = 0.9472135954999579 + t7 + t1;
      double t9 = t8 * t6;
      double t10 = 0.5278640450004206E-1 + t7 + t1;
      double t12 = t10 * t9 * t5;
      double t14 = 0.1E1 * y;
      double t15 = -t3 - t14 - 0.5278640450004206E-1;
      double t17 = -t3 - t14 - 0.9472135954999579;
      double t23 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t24 = t23 * t6;
      double t27 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t31 = t6 * t5;
      double t33 = t10 * t8;
      double t41 = -t3 - t14 + 0.1546536707079771;
      double t42 = -t3 - t14 - 0.5;
      double t44 = -t3 - t14 - 0.1154653670707977E1;
      double t45 = t44 * t42 * t41;
      double t48 = 0.1154653670707977E1 + t7 + t1;
      double t49 = 0.5 + t7 + t1;
      double t50 = t49 * t48;
      double t51 = -0.1546536707079771 + t7 + t1;
      double t55 = -t3 - t14 + 0.2650553239294647;
      double t56 = -t3 - t14 - 0.2147684835193549;
      double t58 = -t3 - t14 - 0.7852315164806451;
      double t59 = -t3 - t14 - 0.1265055323929465E1;
      double t61 = t59 * t58 * t56 * t55;
      double t64 = t17 * t15;
      double t68 = 0.1265055323929465E1 + t7 + t1;
      double t69 = 0.7852315164806451 + t7 + t1;
      double t71 = 0.2147684835193549 + t7 + t1;
      double t72 = -0.2650553239294647 + t7 + t1;
      double t82 = t44 * t42;
      double t101 = -t3 - t14 + 0.3717401485096066;
      double t112 =
        (-t3 - t14 - 0.1371740148509607E1) * (-t3 - t14 -
                0.1091700181433142E1) * (-t3 - t14 -
                       0.7092992179024789)
        * (-t3 - t14 - 0.2907007820975211) * (-t3 - t14 + 0.917001814331423E-1);
      double MapleGenVar1 =
        0.1179889359572706E4 * t12 + 0.3644895809114196E4 * t17 * t15 * t6 * t5 +
        0.4239974851768125E3 * t27 * t24 * t5 + 0.3864320576256994E3 * t31 -
        0.1529327582048436E4 * t27 * t33 * t31 -
        0.2646255087966301E4 * t17 * t15 * t23 * t31 +
        0.6394431503163879E4 * t45 * t31;
      double t115 =
        MapleGenVar1 + 0.5665161428070201E3 * t51 * t50 * t31 +
        0.8640439762651027E4 * t61 * t31 +
        0.1027134453831671E5 * t64 * t33 * t31 -
        0.383793845985807E2 * t72 * t71 * t69 * t68 * t31 -
        0.1118221465928798E4 * t27 * t51 * t50 * t31 -
        0.5316822559346911E4 * t82 * t41 * t23 * t31 +
        0.4842289924427221E4 * (-t3 - t14 - 0.139975799541146E1) * (-t3 - t14 -
                    0.1177186279510738E1)
        * (-t3 - t14 - 0.8631174638261782) * t42 * (-t3 - t14 -
                0.1368825361738218) * (-t3 -
                           t14 +
                           0.1771862795107378)
        * (-t3 - t14 + 0.3997579954114602) * t6 * t5 -
        0.5816387269908822E4 * t112 * t101 * t24 * t5;
      double t116 = -t3 - t14 + 0.3302238962785669;
      double t117 = -t3 - t14 - 0.3115120652928579E-1;
      double t119 = -t3 - t14 - 0.9688487934707142;
      double t121 = -t3 - t14 - 0.1330223896278567E1;
      double t123 = t121 * t119 * t42 * t117 * t116;
      double t145 = t48 * t6 * t5;
      double t146 = t51 * t49;
      double t152 = (0.1330223896278567E1 + t7 + t1) * t6 * t5;
      double t154 = t49 * (0.9688487934707142 + t7 + t1);
      double t157 =
        (-0.3302238962785669 + t7 + t1) * (0.3115120652928579E-1 + t7 + t1);
      double t162 = t9 * t5;
      double t170 = t68 * t6 * t5;
      double t171 = t71 * t69;
      double t177 = t24 * t5;
      MapleGenVar1 =
        0.1742158675177348E5 * t123 * t12 -
        0.1129686939771449E2 * (-0.3717401485096066 + t7 +
              t1) * (-0.917001814331423E-1 + t7 +
               t1) * (0.2907007820975211 + t7 +
                t1) * (0.7092992179024789 + t7 +
                 t1) * (0.1091700181433142E1 +
                  t7 +
                  t1) *
        (0.1371740148509607E1 + t7 + t1) * t6 * t5 +
        0.5698761970291791E4 * t112 * t101 * t6 * t5 +
        0.1101065341415125E5 * t45 * t146 * t145 +
        0.3178691657686378E2 * t27 * t157 * t154 * t152 +
        0.1872758552635813E5 * t59 * t58 * t56 * t55 * t10 * t162 -
        0.247599194517437E3 * t17 * t15 * t72 * t171 * t170 -
        0.5554261285327401E4 * t123 * t177;
      double t207 =
        MapleGenVar1 + 0.158414173059413E5 * t82 * t41 * t10 * t162 -
        0.6182991517419212E4 * t61 * t177 +
        0.4239259297140644E1 * t157 * t154 * t152 -
        0.3795249324381826E3 * t27 * t72 * t171 * t170 +
        0.2947695329493827E4 * t64 * t146 * t145 +
        0.1009757207040712E5 * t121 * t119 * t42 * t117 * t116 * t6 * t5 -
        0.2998967244915291E3 * t177 - 0.4858762230760809E3 * t27 * t6 * t5;
      double t208 = t115 + t207;
      return t208;
    }

    double
    ortho2_f61x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * x;
      double t7 = 0.1154653670707977E1 + t6 + t1;
      double t8 = t7 * t5;
      double t9 = 0.5 + t6 + t1;
      double t10 = -0.1546536707079771 + t6 + t1;
      double t11 = t10 * t9;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 - 0.5278640450004206E-1;
      double t14 = -t3 - t12 - 0.9472135954999579;
      double t15 = t14 * t13;
      double t16 = t15 * t11;
      double t20 = (-t3 - t1) * t2;
      double t21 = t4 * t20;
      double t22 = -t3 - t12 + 0.2650553239294647;
      double t23 = -t3 - t12 - 0.2147684835193549;
      double t24 = t23 * t22;
      double t25 = -t3 - t12 - 0.7852315164806451;
      double t26 = -t3 - t12 - 0.1265055323929465E1;
      double t27 = t26 * t25;
      double t28 = t27 * t24;
      double t33 = t7 * t20;
      double t36 = -t3 - t12 + 0.3302238962785669;
      double t37 = -t3 - t12 - 0.3115120652928579E-1;
      double t38 = t37 * t36;
      double t39 = -t3 - t12 - 0.9688487934707142;
      double t40 = -t3 - t12 - 0.1330223896278567E1;
      double t41 = t40 * t39;
      double t42 = t41 * t38;
      double t45 = -t3 - t12 - 0.5;
      double t47 = t41 * t45 * t36;
      double t50 = t45 * t37;
      double t51 = t41 * t50;
      double t60 = 0.1265055323929465E1 + t6 + t1;
      double t61 = 0.7852315164806451 + t6 + t1;
      double t62 = t61 * t60;
      double t63 = 0.2147684835193549 + t6 + t1;
      double t64 = -0.2650553239294647 + t6 + t1;
      double t65 = t64 * t63;
      double t69 = t9 * t7;
      double t78 =
        -0.1473847664746913E4 * t16 * t8 + 0.1955233594854536E5 * t28 * t21 +
        0.2947695329493827E4 * t16 * t21 + 0.1473847664746913E4 * t16 * t33 -
        0.5048786035203561E4 * t42 * t21 - 0.5048786035203561E4 * t47 * t21 -
        0.5048786035203561E4 * t51 * t21 +
        0.5048786035203561E4 * t51 * t36 * t20 -
        0.5048786035203561E4 * t51 * t36 * t5 -
        0.6000816076630892E3 * t65 * t62 * t21 -
        0.1473847664746913E4 * t13 * t10 * t69 * t21 -
        0.1473847664746913E4 * t14 * t10 * t69 * t21;
      double t82 = t10 * t7;
      double t86 = 0.9472135954999579 + t6 + t1;
      double t87 = 0.5278640450004206E-1 + t6 + t1;
      double t88 = t87 * t86;
      double t89 = -t3 - t12 + 0.1546536707079771;
      double t90 = t45 * t89;
      double t94 = -t3 - t12 - 0.1154653670707977E1;
      double t95 = t94 * t89;
      double t99 = t94 * t45;
      double t108 = t99 * t89 * t87;
      double t111 = t86 * t20;
      double t114 = t86 * t5;
      double t119 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t120 = t22 * t119;
      double t121 = t25 * t23;
      double t125 = t26 * t23;
      double t132 =
        0.2947695329493827E4 * t15 * t69 * t21 +
        0.2947695329493827E4 * t15 * t82 * t21 -
        0.792070865297065E4 * t90 * t88 * t21 -
        0.792070865297065E4 * t95 * t88 * t21 -
        0.792070865297065E4 * t99 * t88 * t21 +
        0.158414173059413E5 * t99 * t89 * t86 * t21 +
        0.158414173059413E5 * t108 * t21 + 0.792070865297065E4 * t108 * t111 -
        0.792070865297065E4 * t108 * t114 +
        0.3091495758709606E4 * t121 * t120 * t21 +
        0.3091495758709606E4 * t125 * t120 * t21 +
        0.3091495758709606E4 * t27 * t120 * t21;
      double t138 = t119 * t20;
      double t141 = t119 * t5;
      double t144 = t39 * t45;
      double t145 = t144 * t38;
      double t149 = t40 * t45 * t38;
      double t152 = 0.1330223896278567E1 + t6 + t1;
      double t153 = 0.9688487934707142 + t6 + t1;
      double t154 = t153 * t152;
      double t155 = 0.3115120652928579E-1 + t6 + t1;
      double t156 = t155 * t9;
      double t160 = -0.3302238962785669 + t6 + t1;
      double t165 = t160 * t155;
      double t173 = t9 * t153;
      double t174 = t165 * t173;
      double t183 =
        0.3091495758709606E4 * t27 * t23 * t119 * t21 -
        0.3091495758709606E4 * t28 * t138 + 0.3091495758709606E4 * t28 * t141 -
        0.5048786035203561E4 * t145 * t21 - 0.5048786035203561E4 * t149 * t21 +
        0.4239259297140644E1 * t156 * t154 * t21 +
        0.4239259297140644E1 * t160 * t9 * t154 * t21 +
        0.4239259297140644E1 * t165 * t154 * t21 +
        0.4239259297140644E1 * t165 * t9 * t152 * t21 +
        0.4239259297140644E1 * t174 * t21 +
        0.2119629648570322E1 * t174 * t152 * t20 -
        0.2119629648570322E1 * t174 * t152 * t5;
      double t184 = t10 * t69;
      double t187 = t63 * t61;
      double t190 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t191 = t190 * t64;
      double t192 = t191 * t187;
      double t195 = t60 * t20;
      double t198 = t60 * t5;
      double t204 = t190 * t11;
      double t208 = t26 * t121;
      double t211 = t94 * t90;
      double t221 = t63 * t60;
      double t228 =
        -0.1768063380413687E4 * t184 * t21 - 0.3795249324381826E3 * t192 * t21 -
        0.1897624662190913E3 * t192 * t195 + 0.1897624662190913E3 * t192 * t198 -
        0.1118221465928798E4 * t190 * t82 * t21 -
        0.1118221465928798E4 * t204 * t21 -
        0.4320219881325513E4 * t208 * t22 * t5 +
        0.1681326920250201E5 * t211 * t21 -
        0.3795249324381826E3 * t190 * t63 * t62 * t21 -
        0.3795249324381826E3 * t191 * t62 * t21 -
        0.3795249324381826E3 * t191 * t221 * t21 -
        0.1118221465928798E4 * t190 * t69 * t21;
      double t258 = t89 * t119;
      double t265 =
        -0.4320219881325513E4 * t25 * t24 * t21 -
        0.4320219881325513E4 * t26 * t24 * t21 -
        0.4320219881325513E4 * t26 * t25 * t22 * t21 -
        0.4320219881325513E4 * t208 * t21 +
        0.4320219881325513E4 * t208 * t22 * t20 +
        0.2658411279673456E4 * t94 * t45 * t119 * t21 -
        0.2658411279673456E4 * t211 * t138 + 0.2658411279673456E4 * t211 * t141 -
        0.559110732964399E3 * t204 * t33 + 0.559110732964399E3 * t204 * t8 +
        0.2658411279673456E4 * t94 * t258 * t21 +
        0.2658411279673456E4 * t45 * t258 * t21;
      double t269 = t64 * t187;
      double t283 = t14 * t13 * t87;
      double t297 =
        -0.383793845985807E2 * t64 * t221 * t21 -
        0.383793845985807E2 * t269 * t21 - 0.1918969229929035E2 * t269 * t195 +
        0.1918969229929035E2 * t269 * t198 -
        0.5135672269158357E4 * t13 * t88 * t21 -
        0.5135672269158357E4 * t14 * t88 * t21 +
        0.1027134453831671E5 * t283 * t21 + 0.5135672269158357E4 * t283 * t111 -
        0.5135672269158357E4 * t283 * t114 -
        0.383793845985807E2 * t64 * t62 * t21 -
        0.383793845985807E2 * t63 * t62 * t21 + 0.1801189492969945E3 * t21;
      double t305 = -t3 - t12 + 0.3717401485096066;
      double t306 = t305 * t119;
      double t308 = -t3 - t12 + 0.917001814331423E-1;
      double t309 = -t3 - t12 - 0.2907007820975211;
      double t310 = t309 * t308;
      double t311 = -t3 - t12 - 0.7092992179024789;
      double t312 = -t3 - t12 - 0.1091700181433142E1;
      double t313 = t312 * t311;
      double t314 = -t3 - t12 - 0.1371740148509607E1;
      double t315 = t314 * t313;
      double t316 = t315 * t310;
      double t319 = t119 * t4;
      double t320 = t319 * t20;
      double t330 = t308 * t305;
      double t339 = t311 * t309;
      double t348 =
        0.1499483622457646E3 * t141 - 0.2429381115380405E3 * t190 * t20 +
        0.2429381115380405E3 * t190 * t5 - 0.1499483622457646E3 * t138 +
        0.2908193634954411E4 * t316 * t306 * t5 +
        0.2908193634954411E4 * t315 * t309 * t305 * t320 +
        0.2908193634954411E4 * t316 * t320 -
        0.2908193634954411E4 * t316 * t306 * t20 +
        0.2908193634954411E4 * t314 * t312 * t309 * t330 * t320 +
        0.2908193634954411E4 * t315 * t330 * t320 +
        0.2908193634954411E4 * t312 * t339 * t330 * t320 +
        0.2908193634954411E4 * t314 * t339 * t330 * t320;
      double t349 = -t3 - t12 + 0.3997579954114602;
      double t350 = -t3 - t12 + 0.1771862795107378;
      double t351 = t350 * t349;
      double t353 = -t3 - t12 - 0.1368825361738218;
      double t355 = -t3 - t12 - 0.8631174638261782;
      double t356 = -t3 - t12 - 0.1177186279510738E1;
      double t358 = -t3 - t12 - 0.139975799541146E1;
      double t359 = t358 * t356 * t355;
      double t360 = t359 * t45 * t353;
      double t367 = t305 * t4 * t20;
      double t375 = t349 * t4 * t20;
      double t378 = t353 * t350;
      double t379 = t355 * t45;
      double t400 = t86 * t4;
      double t401 = t400 * t20;
      double t402 = t36 * t87;
      double t409 = t190 * t160;
      double t410 = t409 * t156;
      double t413 = t40 * t144;
      double t417 =
        0.242114496221361E4 * t360 * t351 * t20 -
        0.242114496221361E4 * t360 * t351 * t5 +
        0.1839303152652042E5 * t316 * t367 -
        0.242114496221361E4 * t360 * t350 * t4 * t20 -
        0.242114496221361E4 * t360 * t375 -
        0.242114496221361E4 * t358 * t379 * t378 * t375 -
        0.242114496221361E4 * t358 * t356 * t45 * t378 * t375 -
        0.242114496221361E4 * t359 * t378 * t375 -
        0.242114496221361E4 * t359 * t45 * t350 * t375 -
        0.242114496221361E4 * t356 * t379 * t378 * t375 -
        0.8710793375886742E4 * t40 * t39 * t37 * t402 * t401 -
        0.1589345828843189E2 * t410 * t154 * t5 -
        0.8710793375886742E4 * t413 * t402 * t401;
      double t425 = t413 * t38;
      double t428 = t87 * t4;
      double t429 = t428 * t20;
      double t432 = t88 * t20;
      double t435 = t88 * t5;
      double t447 = t152 * t4 * t20;
      double t465 =
        -0.8710793375886742E4 * t413 * t37 * t87 * t401 +
        0.1742158675177348E5 * t425 * t401 + 0.1742158675177348E5 * t425 * t429 +
        0.8710793375886742E4 * t425 * t432 - 0.8710793375886742E4 * t425 * t435 -
        0.8710793375886742E4 * t39 * t50 * t402 * t401 -
        0.8710793375886742E4 * t40 * t50 * t402 * t401 +
        0.3178691657686378E2 * t190 * t155 * t173 * t447 +
        0.3178691657686378E2 * t409 * t173 * t447 +
        0.3178691657686378E2 * t409 * t155 * t153 * t447 +
        0.3178691657686378E2 * t410 * t447 +
        0.3178691657686378E2 * t410 * t153 * t4 * t20;
      double t471 = 0.1371740148509607E1 + t6 + t1;
      double t473 = t471 * t4 * t20;
      double t474 = 0.7092992179024789 + t6 + t1;
      double t475 = 0.2907007820975211 + t6 + t1;
      double t477 = -0.917001814331423E-1 + t6 + t1;
      double t478 = -0.3717401485096066 + t6 + t1;
      double t479 = t478 * t477;
      double t480 = t479 * t475 * t474;
      double t483 = 0.1091700181433142E1 + t6 + t1;
      double t488 = t483 * t471;
      double t499 = t474 * t483;
      double t505 = t99 * t89 * t10;
      double t508 = t7 * t4;
      double t509 = t508 * t20;
      double t512 = t9 * t4;
      double t519 =
        0.1589345828843189E2 * t410 * t154 * t20 +
        0.5025952808832613E2 * t174 * t447 - 0.1129686939771449E2 * t480 * t473 -
        0.1129686939771449E2 * t480 * t483 * t4 * t20 -
        0.5648434698857245E1 * t480 * t488 * t20 +
        0.5648434698857245E1 * t480 * t488 * t5 -
        0.1129686939771449E2 * t479 * t475 * t483 * t473 -
        0.1129686939771449E2 * t479 * t499 * t473 -
        0.5505326707075624E4 * t505 * t69 * t5 +
        0.1101065341415125E5 * t505 * t509 +
        0.1101065341415125E5 * t505 * t512 * t20 +
        0.5505326707075624E4 * t505 * t69 * t20;
      double t526 = t314 * t312;
      double t527 = t526 * t339;
      double t551 = t15 * t65;
      double t562 = t60 * t4 * t20;
      double t565 =
        0.1101065341415125E5 * t99 * t89 * t9 * t509 -
        0.2849380985145896E4 * t527 * t330 * t5 -
        0.2849380985145896E4 * t527 * t308 * t4 * t20 +
        0.2849380985145896E4 * t527 * t330 * t20 -
        0.5505326707075624E4 * t99 * t11 * t509 +
        0.1756411638132895E5 * t51 * t36 * t4 * t20 -
        0.5505326707075624E4 * t95 * t11 * t509 -
        0.5505326707075624E4 * t90 * t11 * t509 -
        0.1237995972587185E3 * t551 * t62 * t20 +
        0.1237995972587185E3 * t551 * t62 * t5 -
        0.1129686939771449E2 * t478 * t475 * t499 * t473 -
        0.247599194517437E3 * t551 * t562;
      double t604 = t36 * t119;
      double t608 =
        -0.247599194517437E3 * t551 * t61 * t4 * t20 -
        0.1129686939771449E2 * t477 * t475 * t499 * t473 +
        0.1237995972587185E3 * t14 * t64 * t187 * t562 -
        0.247599194517437E3 * t15 * t187 * t562 -
        0.247599194517437E3 * t15 * t64 * t61 * t562 +
        0.1237995972587185E3 * t13 * t64 * t187 * t562 -
        0.2849380985145896E4 * t527 * t367 -
        0.2849380985145896E4 * t526 * t311 * t308 * t367 +
        0.1822447904557098E4 * t15 * t20 -
        0.2849380985145896E4 * t526 * t310 * t367 -
        0.2849380985145896E4 * t314 * t311 * t310 * t367 +
        0.2777130642663701E4 * t51 * t604 * t5;
      double t631 = t22 * t87;
      double t642 =
        -0.2849380985145896E4 * t313 * t310 * t367 -
        0.2777130642663701E4 * t51 * t604 * t20 +
        0.2777130642663701E4 * t51 * t320 + 0.2777130642663701E4 * t47 * t320 +
        0.2777130642663701E4 * t42 * t320 + 0.2777130642663701E4 * t145 * t320 +
        0.2777130642663701E4 * t149 * t320 + 0.9363792763179067E4 * t28 * t432 -
        0.9363792763179067E4 * t28 * t435 -
        0.9363792763179067E4 * t125 * t631 * t401 -
        0.9363792763179067E4 * t27 * t631 * t401 -
        0.9363792763179067E4 * t27 * t23 * t87 * t401;
      double t653 = t13 * t4;
      double t659 = t190 * t119;
      double t668 =
        0.1872758552635813E5 * t28 * t401 + 0.1872758552635813E5 * t28 * t429 -
        0.9363792763179067E4 * t121 * t631 * t401 -
        0.1822447904557098E4 * t14 * t4 * t20 -
        0.1822447904557098E4 * t653 * t20 -
        0.1340797775342208E4 * t190 * t4 * t20 +
        0.2119987425884062E3 * t659 * t20 - 0.2119987425884062E3 * t659 * t5 +
        0.6703988876711038E3 * t320 - 0.589944679786353E3 * t435 +
        0.1179889359572706E4 * t429 + 0.1179889359572706E4 * t401;
      double t685 = t190 * t88;
      double t698 =
        0.589944679786353E3 * t432 +
        0.1027134453831671E5 * t14 * t13 * t86 * t21 -
        0.1822447904557098E4 * t15 * t5 -
        0.2418079223895547E4 * t87 * t400 * t20 +
        0.8368193347782743E4 * t14 * t653 * t20 + 0.1932160288128497E3 * t20 -
        0.1932160288128497E3 * t5 + 0.7646637910242182E3 * t685 * t5 -
        0.7646637910242182E3 * t685 * t20 -
        0.1529327582048436E4 * t190 * t428 * t20 -
        0.1529327582048436E4 * t190 * t400 * t20 - 0.28325807140351E3 * t184 * t5;
      double t703 = t89 * t4;
      double t726 = t14 * t13 * t119;
      double t737 =
        -0.319721575158194E4 * t94 * t45 * t4 * t20 -
        0.319721575158194E4 * t94 * t703 * t20 - 0.319721575158194E4 * t211 * t5 +
        0.319721575158194E4 * t211 * t20 -
        0.319721575158194E4 * t45 * t703 * t20 + 0.28325807140351E3 * t184 * t20 +
        0.5665161428070201E3 * t10 * t512 * t20 +
        0.5665161428070201E3 * t10 * t508 * t20 +
        0.5665161428070201E3 * t9 * t508 * t20 +
        0.1323127543983151E4 * t726 * t5 - 0.1323127543983151E4 * t726 * t20 +
        0.1323127543983151E4 * t14 * t319 * t20 +
        0.1323127543983151E4 * t13 * t319 * t20;
      double t741 =
        t78 + t132 + t183 + t228 + t265 + t297 + t348 + t417 + t465 + t519 +
        t565 + t608 + t642 + t668 + t698 + t737;
      return t741;
    }

    double
    ortho2_f61y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.1265055323929465E1 + t8 + t1;
      double t10 = 0.7852315164806451 + t8 + t1;
      double t11 = t10 * t9;
      double t12 = 0.2147684835193549 + t8 + t1;
      double t13 = -0.2650553239294647 + t8 + t1;
      double t14 = t13 * t12;
      double t18 = 0.9472135954999579 + t8 + t1;
      double t19 = t18 * t6;
      double t20 = 0.5278640450004206E-1 + t8 + t1;
      double t24 = t6 * t4;
      double t25 = 0.1154653670707977E1 + t8 + t1;
      double t26 = t25 * t24;
      double t27 = 0.5 + t8 + t1;
      double t28 = -0.1546536707079771 + t8 + t1;
      double t29 = t28 * t27;
      double t32 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t33 = t32 * t29;
      double t36 = t6 * t2;
      double t39 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t40 = t39 * t36;
      double t41 = 0.1E1 * y;
      double t42 = -t3 - t41 + 0.1546536707079771;
      double t43 = -t3 - t41 - 0.5;
      double t44 = t43 * t42;
      double t45 = -t3 - t41 - 0.1154653670707977E1;
      double t46 = t45 * t44;
      double t50 = t9 * t6 * t5;
      double t51 = -t3 - t41 - 0.5278640450004206E-1;
      double t52 = -t3 - t41 - 0.9472135954999579;
      double t53 = t52 * t51;
      double t54 = t53 * t14;
      double t57 = -t3 - t41 + 0.2650553239294647;
      double t58 = -t3 - t41 - 0.2147684835193549;
      double t59 = t58 * t57;
      double t60 = -t3 - t41 - 0.7852315164806451;
      double t61 = -t3 - t41 - 0.1265055323929465E1;
      double t62 = t61 * t60;
      double t63 = t62 * t59;
      double t70 = t12 * t10;
      double t74 = t39 * t6;
      double t75 = t74 * t5;
      double t76 = -t3 - t41 + 0.3302238962785669;
      double t77 = -t3 - t41 - 0.3115120652928579E-1;
      double t78 = t77 * t76;
      double t79 = -t3 - t41 - 0.9688487934707142;
      double t80 = t79 * t43;
      double t81 = t80 * t78;
      double t84 = t39 * t24;
      double t87 = t42 * t6;
      double t95 =
        -0.1200163215326178E4 * t14 * t11 * t7 -
        0.4836158447791095E4 * t20 * t19 * t5 - 0.559110732964399E3 * t33 * t26 +
        0.2658411279673456E4 * t46 * t40 - 0.1237995972587185E3 * t54 * t50 +
        0.3091495758709606E4 * t63 * t40 -
        0.1237995972587185E3 * t53 * t13 * t10 * t50 -
        0.1237995972587185E3 * t53 * t70 * t50 +
        0.5554261285327401E4 * t81 * t75 - 0.3091495758709606E4 * t63 * t84 -
        0.6394431503163879E4 * t43 * t87 * t5 +
        0.247599194517437E3 * t52 * t13 * t70 * t50;
      double t96 = 0.1091700181433142E1 + t8 + t1;
      double t99 = 0.7092992179024789 + t8 + t1;
      double t100 = 0.2907007820975211 + t8 + t1;
      double t102 = -0.917001814331423E-1 + t8 + t1;
      double t103 = -0.3717401485096066 + t8 + t1;
      double t104 = t103 * t102;
      double t105 = t104 * t100 * t99;
      double t116 = t20 * t18;
      double t117 = t116 * t36;
      double t118 = -t3 - t41 - 0.1330223896278567E1;
      double t119 = t118 * t80;
      double t120 = t119 * t78;
      double t123 = t20 * t6;
      double t124 = t123 * t5;
      double t127 = t43 * t77;
      double t128 = t118 * t79;
      double t129 = t128 * t127;
      double t132 = 0.1371740148509607E1 + t8 + t1;
      double t134 = t132 * t6 * t5;
      double t141 = t19 * t5;
      double t144 = t9 * t36;
      double t145 = t13 * t70;
      double t152 = t76 * t20;
      double t156 =
        -0.5648434698857245E1 * t105 * t96 * t6 * t5 +
        0.247599194517437E3 * t51 * t13 * t70 * t50 +
        0.6182991517419212E4 * t62 * t58 * t39 * t7 -
        0.8710793375886742E4 * t120 * t117 + 0.8710793375886742E4 * t120 * t124 -
        0.1009757207040712E5 * t129 * t7 - 0.5648434698857245E1 * t105 * t134 +
        0.5316822559346911E4 * t45 * t43 * t39 * t7 +
        0.8710793375886742E4 * t120 * t141 + 0.1918969229929035E2 * t145 * t144 -
        0.1742158675177348E5 * t119 * t77 * t20 * t141 -
        0.1742158675177348E5 * t119 * t152 * t141;
      double t158 = -t3 - t41 + 0.3717401485096066;
      double t159 = t158 * t39;
      double t161 = -t3 - t41 + 0.917001814331423E-1;
      double t162 = -t3 - t41 - 0.2907007820975211;
      double t163 = t162 * t161;
      double t164 = -t3 - t41 - 0.7092992179024789;
      double t165 = -t3 - t41 - 0.1091700181433142E1;
      double t166 = t165 * t164;
      double t167 = -t3 - t41 - 0.1371740148509607E1;
      double t168 = t167 * t166;
      double t169 = t168 * t163;
      double t172 = t128 * t78;
      double t180 = t42 * t39;
      double t185 = t158 * t6 * t5;
      double t186 = t164 * t162;
      double t187 = t167 * t165;
      double t188 = t187 * t186;
      double t197 = t32 * t39;
      double t200 = t18 * t36;
      double t202 = t45 * t43;
      double t203 = t202 * t42 * t20;
      double t207 = t60 * t58;
      double t208 = t61 * t207;
      double t211 =
        -0.2908193634954411E4 * t169 * t159 * t24 +
        0.5554261285327401E4 * t172 * t75 - 0.2658411279673456E4 * t46 * t84 -
        0.1237995972587185E3 * t54 * t11 * t24 +
        0.5316822559346911E4 * t45 * t180 * t7 -
        0.5698761970291791E4 * t188 * t185 +
        0.2908193634954411E4 * t169 * t159 * t36 +
        0.5816387269908822E4 * t169 * t75 + 0.1932160288128497E3 * t24 -
        0.2119987425884062E3 * t197 * t36 - 0.792070865297065E4 * t203 * t200 -
        0.4320219881325513E4 * t208 * t57 * t36;
      double t212 = t27 * t25;
      double t213 = t28 * t212;
      double t220 = t128 * t43 * t76;
      double t223 = t25 * t6;
      double t224 = t223 * t5;
      double t228 = t57 * t39;
      double t229 = t61 * t58;
      double t236 = 0.1330223896278567E1 + t8 + t1;
      double t237 = 0.9688487934707142 + t8 + t1;
      double t238 = t237 * t236;
      double t240 = 0.3115120652928579E-1 + t8 + t1;
      double t241 = t240 * t27;
      double t242 = -0.3302238962785669 + t8 + t1;
      double t243 = t32 * t242;
      double t244 = t243 * t241;
      double t256 = t99 * t96;
      double t260 =
        -0.28325807140351E3 * t213 * t36 - 0.1499483622457646E3 * t84 +
        0.9776167974272681E4 * t63 * t7 - 0.1009757207040712E5 * t220 * t7 -
        0.1101065341415125E5 * t202 * t29 * t224 +
        0.6182991517419212E4 * t229 * t228 * t7 +
        0.6182991517419212E4 * t62 * t228 * t7 -
        0.1589345828843189E2 * t244 * t238 * t36 +
        0.1589345828843189E2 * t244 * t237 * t6 * t5 -
        0.1932160288128497E3 * t36 -
        0.5648434698857245E1 * t104 * t100 * t96 * t134 -
        0.5648434698857245E1 * t104 * t256 * t134;
      double t263 = t51 * t6;
      double t267 = t18 * t24;
      double t269 = t52 * t51 * t20;
      double t274 = t96 * t132;
      double t278 = t32 * t13;
      double t279 = t278 * t70;
      double t286 = t12 * t9;
      double t293 = t28 * t25;
      double t298 = t236 * t6 * t5;
      double t306 =
        0.4184096673891372E4 * t52 * t263 * t5 +
        0.5135672269158357E4 * t269 * t267 - 0.559110732964399E3 * t33 * t7 -
        0.5648434698857245E1 * t105 * t274 * t24 +
        0.1897624662190913E3 * t279 * t144 - 0.1897624662190913E3 * t279 * t7 -
        0.1918969229929035E2 * t145 * t7 -
        0.1918969229929035E2 * t13 * t286 * t7 -
        0.1918969229929035E2 * t13 * t11 * t7 -
        0.559110732964399E3 * t32 * t293 * t7 +
        0.1589345828843189E2 * t244 * t298 -
        0.1742158675177348E5 * t118 * t79 * t77 * t152 * t141;
      double t317 = t45 * t42;
      double t326 = t202 * t42 * t28;
      double t329 = t27 * t6;
      double t340 = t116 * t24;
      double t342 =
        -0.1742158675177348E5 * t118 * t127 * t152 * t141 +
        0.792070865297065E4 * t203 * t7 - 0.5135672269158357E4 * t269 * t200 +
        0.5135672269158357E4 * t269 * t7 -
        0.1101065341415125E5 * t317 * t29 * t224 -
        0.1101065341415125E5 * t44 * t29 * t224 +
        0.5505326707075624E4 * t326 * t212 * t24 +
        0.28325807140351E3 * t28 * t329 * t5 +
        0.28325807140351E3 * t28 * t223 * t5 - 0.589944679786353E3 * t117 -
        0.1897624662190913E3 * t278 * t286 * t7 + 0.589944679786353E3 * t340;
      double t351 = t32 * t116;
      double t360 = t161 * t158;
      double t368 = t25 * t36;
      double t369 = t53 * t29;
      double t383 = t118 * t43 * t78;
      double t386 =
        0.5135672269158357E4 * t52 * t51 * t18 * t7 +
        0.28325807140351E3 * t27 * t223 * t5 - 0.7646637910242182E3 * t351 * t24 -
        0.5505326707075624E4 * t326 * t212 * t36 +
        0.1237995972587185E3 * t54 * t11 * t36 -
        0.2849380985145896E4 * t188 * t360 * t36 -
        0.5698761970291791E4 * t188 * t161 * t6 * t5 -
        0.1473847664746913E4 * t369 * t368 +
        0.5816387269908822E4 * t168 * t360 * t75 +
        0.5816387269908822E4 * t168 * t162 * t158 * t75 -
        0.1897624662190913E3 * t278 * t11 * t7 - 0.1009757207040712E5 * t383 * t7;
      double t401 = t27 * t237;
      double t421 =
        0.5505326707075624E4 * t326 * t329 * t5 +
        0.1473847664746913E4 * t369 * t7 +
        0.1473847664746913E4 * t53 * t293 * t7 +
        0.1589345828843189E2 * t243 * t240 * t237 * t298 -
        0.1822447904557098E4 * t53 * t36 +
        0.1589345828843189E2 * t243 * t401 * t298 -
        0.3536126760827375E4 * t213 * t7 - 0.1009757207040712E5 * t172 * t7 -
        0.1918969229929035E2 * t12 * t11 * t7 - 0.2429381115380405E3 * t32 * t24 +
        0.2429381115380405E3 * t32 * t36 - 0.559110732964399E3 * t32 * t212 * t7 +
        0.5505326707075624E4 * t326 * t224;
      double t445 = t242 * t240;
      double t446 = t445 * t401;
      double t460 =
        0.5505326707075624E4 * t202 * t42 * t27 * t224 +
        0.1589345828843189E2 * t32 * t240 * t401 * t298 -
        0.5698761970291791E4 * t187 * t164 * t161 * t185 -
        0.5698761970291791E4 * t187 * t163 * t185 -
        0.5698761970291791E4 * t167 * t164 * t163 * t185 -
        0.2119629648570322E1 * t446 * t236 * t36 -
        0.3644895809114196E4 * t263 * t5 + 0.1822447904557098E4 * t53 * t24 -
        0.6703988876711038E3 * t32 * t6 * t5 + 0.1340797775342208E4 * t75 +
        0.2119987425884062E3 * t197 * t24 + 0.589944679786353E3 * t141;
      double t504 =
        0.5816387269908822E4 * t167 * t165 * t162 * t360 * t75 +
        0.5816387269908822E4 * t167 * t186 * t360 * t75 -
        0.1027134453831671E5 * t52 * t116 * t7 +
        0.1473847664746913E4 * t53 * t212 * t7 -
        0.1742158675177348E5 * t79 * t127 * t152 * t141 +
        0.792070865297065E4 * t202 * t42 * t18 * t7 -
        0.5648434698857245E1 * t103 * t100 * t256 * t134 -
        0.1897624662190913E3 * t32 * t12 * t11 * t7 -
        0.9363792763179067E4 * t63 * t117 + 0.9363792763179067E4 * t63 * t124 -
        0.5648434698857245E1 * t102 * t100 * t256 * t134 -
        0.1237995972587185E3 * t54 * t10 * t6 * t5;
      double t512 = -t3 - t41 + 0.3997579954114602;
      double t513 = -t3 - t41 + 0.1771862795107378;
      double t514 = t513 * t512;
      double t516 = -t3 - t41 - 0.1368825361738218;
      double t518 = -t3 - t41 - 0.8631174638261782;
      double t519 = -t3 - t41 - 0.1177186279510738E1;
      double t521 = -t3 - t41 - 0.139975799541146E1;
      double t522 = t521 * t519 * t518;
      double t523 = t522 * t43 * t516;
      double t526 = t76 * t39;
      double t533 = t9 * t24;
      double t552 =
        -0.5698761970291791E4 * t166 * t163 * t185 +
        0.5316822559346911E4 * t43 * t180 * t7 +
        0.242114496221361E4 * t523 * t514 * t24 +
        0.2777130642663701E4 * t129 * t526 * t36 +
        0.6182991517419212E4 * t207 * t228 * t7 -
        0.1897624662190913E3 * t279 * t533 + 0.2119629648570322E1 * t446 * t7 -
        0.1062297169731862E4 * t7 - 0.7646637910242182E3 * t32 * t123 * t5 +
        0.5816387269908822E4 * t165 * t186 * t360 * t75 +
        0.1005190561766523E3 * t446 * t298 +
        0.2119629648570322E1 * t445 * t27 * t236 * t7;
      double t561 = t512 * t6 * t5;
      double t593 =
        -0.242114496221361E4 * t523 * t514 * t36 -
        0.4842289924427221E4 * t523 * t513 * t6 * t5 -
        0.4842289924427221E4 * t523 * t561 - 0.8640439762651027E4 * t208 * t7 -
        0.8640439762651027E4 * t61 * t60 * t57 * t7 +
        0.2119629648570322E1 * t445 * t238 * t7 +
        0.2119629648570322E1 * t242 * t27 * t238 * t7 +
        0.2119629648570322E1 * t446 * t236 * t24 -
        0.8640439762651027E4 * t60 * t59 * t7 -
        0.6394431503163879E4 * t45 * t43 * t6 * t5 -
        0.8640439762651027E4 * t61 * t59 * t7 -
        0.158414173059413E5 * t202 * t116 * t7;
      double t623 = t52 * t51 * t39;
      double t632 =
        -0.2947695329493827E4 * t52 * t28 * t212 * t7 -
        0.2947695329493827E4 * t51 * t28 * t212 * t7 -
        0.4842289924427221E4 * t522 * t43 * t513 * t561 +
        0.589944679786353E3 * t124 + 0.792070865297065E4 * t203 * t267 +
        0.5648434698857245E1 * t105 * t274 * t36 -
        0.7646637910242182E3 * t32 * t19 * t5 + 0.319721575158194E4 * t46 * t24 -
        0.6394431503163879E4 * t45 * t87 * t5 -
        0.1323127543983151E4 * t623 * t24 -
        0.5048786035203561E4 * t129 * t76 * t36 -
        0.158414173059413E5 * t317 * t116 * t7;
      double t643 = t516 * t513;
      double t664 =
        -0.158414173059413E5 * t44 * t116 * t7 +
        0.4320219881325513E4 * t208 * t57 * t24 +
        0.28325807140351E3 * t213 * t24 + 0.9363792763179067E4 * t63 * t340 -
        0.4842289924427221E4 * t522 * t643 * t561 -
        0.1009757207040712E5 * t81 * t7 + 0.1499483622457646E3 * t40 -
        0.3644895809114196E4 * t52 * t6 * t5 + 0.9363792763179067E4 * t63 * t141 +
        0.1323127543983151E4 * t623 * t36 +
        0.2646255087966301E4 * t52 * t74 * t5 -
        0.1872758552635813E5 * t62 * t58 * t20 * t141;
      double t666 = t57 * t20;
      double t698 =
        -0.1872758552635813E5 * t62 * t666 * t141 +
        0.2646255087966301E4 * t51 * t74 * t5 + 0.559110732964399E3 * t33 * t368 -
        0.319721575158194E4 * t46 * t36 +
        0.2119629648570322E1 * t241 * t238 * t7 +
        0.9196515763260209E4 * t169 * t185 +
        0.1589345828843189E2 * t244 * t238 * t24 +
        0.8710793375886742E4 * t120 * t340 + 0.1473847664746913E4 * t369 * t26 +
        0.5554261285327401E4 * t129 * t75 + 0.5554261285327401E4 * t220 * t75 -
        0.4842289924427221E4 * t521 * t519 * t43 * t643 * t561;
      double t699 = t518 * t43;
      double t738 =
        -0.4842289924427221E4 * t521 * t699 * t643 * t561 -
        0.4842289924427221E4 * t519 * t699 * t643 * t561 +
        0.5048786035203561E4 * t129 * t76 * t24 -
        0.1918969229929035E2 * t145 * t533 +
        0.2849380985145896E4 * t188 * t360 * t24 +
        0.8782058190664475E4 * t129 * t76 * t6 * t5 +
        0.7646637910242182E3 * t351 * t36 + 0.8406634601251003E4 * t46 * t7 -
        0.1872758552635813E5 * t229 * t666 * t141 -
        0.1027134453831671E5 * t51 * t116 * t7 +
        0.5554261285327401E4 * t383 * t75 -
        0.1872758552635813E5 * t207 * t666 * t141 -
        0.2777130642663701E4 * t129 * t526 * t24;
      double t742 =
        t95 + t156 + t211 + t260 + t306 + t342 + t386 + t421 + t460 + t504 +
        t552 + t593 + t632 + t664 + t698 + t738;
      return t742;
    }

    // * f62 *********************************************************************

    double
    ortho2_f62 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * y;
      double t8 = -t3 - t7 - 0.5278640450004206E-1;
      double t10 = -t3 - t7 - 0.9472135954999579;
      double t14 = 0.1E1 * x;
      double t15 = 0.9472135954999579 + t14 + t1;
      double t16 = t15 * t6;
      double t17 = 0.5278640450004206E-1 + t14 + t1;
      double t19 = t17 * t16 * t5;
      double t23 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t24 = t23 * t6;
      double t27 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t31 = t6 * t5;
      double t33 = -t3 - t7 + 0.2650553239294647;
      double t34 = -t3 - t7 - 0.2147684835193549;
      double t36 = -t3 - t7 - 0.7852315164806451;
      double t37 = -t3 - t7 - 0.1265055323929465E1;
      double t39 = t37 * t36 * t34 * t33;
      double t42 = 0.1154653670707977E1 + t14 + t1;
      double t43 = 0.5 + t14 + t1;
      double t44 = t43 * t42;
      double t45 = -0.1546536707079771 + t14 + t1;
      double t50 = t17 * t15;
      double t51 = t10 * t8;
      double t55 = 0.1265055323929465E1 + t14 + t1;
      double t56 = 0.7852315164806451 + t14 + t1;
      double t58 = 0.2147684835193549 + t14 + t1;
      double t59 = -0.2650553239294647 + t14 + t1;
      double t64 = -t3 - t7 + 0.1546536707079771;
      double t66 = -t3 - t7 - 0.5;
      double t67 = -t3 - t7 - 0.1154653670707977E1;
      double t68 = t67 * t66;
      double t77 = t67 * t66 * t64;
      double t101 = t42 * t6;
      double t106 = t37 * t36 * t34;
      double t110 = -t3 - t7 + 0.3302238962785669;
      double t111 = -t3 - t7 - 0.3115120652928579E-1;
      double t113 = -t3 - t7 - 0.9688487934707142;
      double t115 = -t3 - t7 - 0.1330223896278567E1;
      double t117 = t115 * t113 * t66 * t111 * t110;
      double MapleGenVar1 =
        0.4966359698523725E4 * t10 * t8 * t6 * t5 + 0.2838141049686084E4 * t19 +
        0.6457003867411289E3 * t27 * t24 * t5 + 0.5867497943922963E3 * t31 +
        0.8994429964242004E4 * t39 * t31 -
        0.2726237759905083E4 * t27 * t45 * t44 * t31 +
        0.2223180731003101E5 * t51 * t50 * t31 +
        0.156324097610188E4 * t59 * t58 * t56 * t55 * t31;
      double t120 =
        MapleGenVar1 - 0.7063898286376077E4 * t68 * t64 * t23 * t31 -
        0.4525128512731688E4 * t10 * t8 * t23 * t31 +
        0.8608730879084671E4 * t77 * t31 +
        0.2189331658876131E4 * t45 * t44 * t31 -
        0.2740504385244097E4 * t27 * t50 * t31 + 0.3656576713920043E4 * (-t3 -
                         t7 -
                         0.139975799541146E1)
        * (-t3 - t7 - 0.1177186279510738E1) * (-t3 - t7 -
                 0.8631174638261782) * t66 * (-t3 -
                      t7 -
                      0.1368825361738218)
        * (-t3 - t7 + 0.1771862795107378) * (-t3 - t7 +
               0.3997579954114602) * t6 * t5 +
        0.2578858190533426E5 * t106 * t33 * t45 * t43 * t101 * t5 +
        0.4095641972506376E5 * t117 * t19;
      double t121 = -t3 - t7 + 0.3717401485096066;
      double t132 =
        (-t3 - t7 - 0.1371740148509607E1) * (-t3 - t7 -
               0.1091700181433142E1) * (-t3 - t7 -
                      0.7092992179024789)
        * (-t3 - t7 - 0.2907007820975211) * (-t3 - t7 + 0.917001814331423E-1);
      double t156 = t24 * t5;
      double t160 = (0.1330223896278567E1 + t14 + t1) * t6 * t5;
      double t162 = t43 * (0.9688487934707142 + t14 + t1);
      double t165 =
        (-0.3302238962785669 + t14 + t1) * (0.3115120652928579E-1 + t14 + t1);
      double t170 = t16 * t5;
      double t178 = t55 * t6 * t5;
      double t179 = t58 * t56;
      double t185 = t101 * t5;
      double t186 = t45 * t43;
      MapleGenVar1 = -0.6741048997173348E4 * t132 * t121 * t24 * t5
        - 0.1998987179858919E2 * (-0.3717401485096066 + t14 +
                t1) * (-0.917001814331423E-1 + t14 +
                 t1) * (0.2907007820975211 + t14 +
                  t1) * (0.7092992179024789 + t14 +
                   t1) *
        (0.1091700181433142E1 + t14 + t1) * (0.1371740148509607E1 + t14 +
               t1) * t6 * t5 +
        0.3715834555594522E4 * t132 * t121 * t6 * t5 -
        0.6546972432240988E3 * t27 * t6 * t5 - 0.5209369827736014E3 * t156 +
        0.1994791224891356E2 * t27 * t165 * t162 * t160 +
        0.2948754515800036E5 * t106 * t33 * t17 * t170 -
        0.5672739849408487E4 * t117 * t156;
      double t213 =
        MapleGenVar1 + 0.1006826808160069E5 * t10 * t8 * t59 * t179 * t178 +
        0.2302483456293612E5 * t77 * t186 * t185 -
        0.1029319606251424E5 * t39 * t156 +
        0.3698591049065542E3 * t165 * t162 * t160 +
        0.3271146426942721E5 * t68 * t64 * t17 * t170 -
        0.4158894483907378E3 * t27 * t59 * t179 * t178 +
        0.125704235001869E5 * t115 * t113 * t66 * t111 * t110 * t6 * t5 +
        0.1537690234755264E5 * t51 * t186 * t185;
      double t214 = t120 + t213;
      return t214;
    }

    double
    ortho2_f62x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.1265055323929465E1 + t8 + t1;
      double t10 = 0.7852315164806451 + t8 + t1;
      double t11 = t10 * t9;
      double t12 = 0.2147684835193549 + t8 + t1;
      double t13 = -0.2650553239294647 + t8 + t1;
      double t14 = t13 * t12;
      double t18 = 0.1154653670707977E1 + t8 + t1;
      double t19 = 0.5 + t8 + t1;
      double t20 = t19 * t18;
      double t21 = -0.1546536707079771 + t8 + t1;
      double t22 = 0.1E1 * y;
      double t23 = -t3 - t22 - 0.5278640450004206E-1;
      double t28 = -t3 - t22 - 0.9472135954999579;
      double t33 = t28 * t23;
      double t37 = t21 * t18;
      double t41 = t21 * t19;
      double t42 = t33 * t41;
      double t45 = t18 * t5;
      double t48 = t6 * t2;
      double t49 = t18 * t48;
      double t52 = -t3 - t22 + 0.2650553239294647;
      double t53 = -t3 - t22 - 0.2147684835193549;
      double t54 = t53 * t52;
      double t55 = -t3 - t22 - 0.7852315164806451;
      double t56 = -t3 - t22 - 0.1265055323929465E1;
      double t57 = t56 * t55;
      double t58 = t57 * t54;
      double t61 = -t3 - t22 + 0.3302238962785669;
      double t62 = -t3 - t22 - 0.5;
      double t64 = -t3 - t22 - 0.9688487934707142;
      double t65 = -t3 - t22 - 0.1330223896278567E1;
      double t66 = t65 * t64;
      double t67 = t66 * t62 * t61;
      double t70 = -t3 - t22 - 0.3115120652928579E-1;
      double t71 = t62 * t70;
      double t72 = t66 * t71;
      double t78 =
        -0.6575789558728902E3 * t14 * t11 * t7 -
        0.7688451173776322E4 * t23 * t21 * t20 * t7 -
        0.7688451173776322E4 * t28 * t21 * t20 * t7 +
        0.1537690234755264E5 * t33 * t20 * t7 +
        0.1537690234755264E5 * t33 * t37 * t7 + 0.1537690234755264E5 * t42 * t7 +
        0.7688451173776322E4 * t42 * t45 - 0.7688451173776322E4 * t42 * t49 +
        0.3254994396022192E5 * t58 * t7 - 0.6285211750093451E4 * t67 * t7 -
        0.6285211750093451E4 * t72 * t7 + 0.6285211750093451E4 * t72 * t61 * t5;
      double t82 = t70 * t61;
      double t83 = t66 * t82;
      double t86 = 0.1330223896278567E1 + t8 + t1;
      double t88 = 0.9688487934707142 + t8 + t1;
      double t89 = t19 * t88;
      double t90 = 0.3115120652928579E-1 + t8 + t1;
      double t91 = -0.3302238962785669 + t8 + t1;
      double t92 = t91 * t90;
      double t93 = t92 * t89;
      double t105 = t88 * t86;
      double t109 = t90 * t19;
      double t117 = t64 * t62;
      double t118 = t117 * t82;
      double t122 = t65 * t62 * t82;
      double t127 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t128 = t52 * t127;
      double t129 = t55 * t53;
      double t133 = t56 * t53;
      double t137 =
        -0.6285211750093451E4 * t72 * t61 * t48 -
        0.6285211750093451E4 * t83 * t7 - 0.1849295524532771E3 * t93 * t86 * t48 +
        0.1849295524532771E3 * t93 * t86 * t5 + 0.3698591049065542E3 * t93 * t7 +
        0.3698591049065542E3 * t92 * t19 * t86 * t7 +
        0.3698591049065542E3 * t92 * t105 * t7 +
        0.3698591049065542E3 * t109 * t105 * t7 +
        0.3698591049065542E3 * t91 * t19 * t105 * t7 -
        0.6285211750093451E4 * t118 * t7 - 0.6285211750093451E4 * t122 * t7 +
        0.5146598031257122E4 * t129 * t128 * t7 +
        0.5146598031257122E4 * t133 * t128 * t7;
      double t146 = t127 * t5;
      double t149 = t127 * t48;
      double t152 = 0.9472135954999579 + t8 + t1;
      double t153 = -t3 - t22 + 0.1546536707079771;
      double t155 = -t3 - t22 - 0.1154653670707977E1;
      double t156 = t155 * t62;
      double t160 = 0.5278640450004206E-1 + t8 + t1;
      double t162 = t156 * t153 * t160;
      double t165 = t152 * t5;
      double t168 = t152 * t48;
      double t171 = t160 * t152;
      double t172 = t155 * t153;
      double t181 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t182 = t181 * t13;
      double t186 = t12 * t9;
      double t190 =
        0.5146598031257122E4 * t57 * t128 * t7 +
        0.5146598031257122E4 * t57 * t53 * t127 * t7 -
        0.5146598031257122E4 * t58 * t146 + 0.5146598031257122E4 * t58 * t149 +
        0.3271146426942721E5 * t156 * t153 * t152 * t7 +
        0.3271146426942721E5 * t162 * t7 + 0.1635573213471361E5 * t165 * t162 -
        0.1635573213471361E5 * t162 * t168 -
        0.1635573213471361E5 * t172 * t171 * t7 -
        0.1635573213471361E5 * t156 * t171 * t7 -
        0.4158894483907378E3 * t182 * t11 * t7 -
        0.4158894483907378E3 * t182 * t186 * t7;
      double t191 = t12 * t10;
      double t192 = t182 * t191;
      double t195 = t9 * t5;
      double t198 = t9 * t48;
      double t201 = t62 * t153;
      double t209 = t155 * t201;
      double t215 = t21 * t20;
      double t218 = t181 * t41;
      double t222 = t56 * t129;
      double t225 = t55 * t54;
      double t228 = t56 * t54;
      double t232 = t56 * t55 * t52;
      double t235 =
        -0.4158894483907378E3 * t192 * t7 - 0.2079447241953689E3 * t192 * t195 +
        0.2079447241953689E3 * t192 * t198 -
        0.1635573213471361E5 * t201 * t171 * t7 -
        0.4158894483907378E3 * t181 * t12 * t11 * t7 +
        0.2233800774470877E5 * t209 * t7 -
        0.2726237759905083E4 * t181 * t20 * t7 -
        0.4310560382227665E4 * t215 * t7 - 0.2726237759905083E4 * t218 * t7 -
        0.4497214982121002E4 * t222 * t52 * t48 -
        0.4497214982121002E4 * t225 * t7 - 0.4497214982121002E4 * t228 * t7 -
        0.4497214982121002E4 * t232 * t7;
      double t247 = t28 * t23 * t160;
      double t250 = t153 * t127;
      double t272 =
        -0.4497214982121002E4 * t222 * t7 +
        0.4497214982121002E4 * t222 * t52 * t5 -
        0.2726237759905083E4 * t181 * t37 * t7 -
        0.111159036550155E5 * t247 * t168 +
        0.3531949143188038E4 * t62 * t250 * t7 +
        0.3531949143188038E4 * t155 * t250 * t7 +
        0.3531949143188038E4 * t155 * t62 * t127 * t7 -
        0.3531949143188038E4 * t209 * t146 + 0.3531949143188038E4 * t209 * t149 -
        0.1363118879952542E4 * t218 * t45 + 0.1363118879952542E4 * t218 * t49 +
        0.156324097610188E4 * t12 * t11 * t7;
      double t279 = t13 * t191;
      double t304 =
        0.156324097610188E4 * t13 * t11 * t7 +
        0.156324097610188E4 * t13 * t186 * t7 + 0.156324097610188E4 * t279 * t7 +
        0.7816204880509401E3 * t279 * t195 - 0.7816204880509401E3 * t279 * t198 -
        0.111159036550155E5 * t23 * t171 * t7 -
        0.111159036550155E5 * t28 * t171 * t7 +
        0.2223180731003101E5 * t28 * t23 * t152 * t7 +
        0.2223180731003101E5 * t247 * t7 + 0.111159036550155E5 * t247 * t165 +
        0.3273486216120494E3 * t181 * t48 + 0.6121801497697838E3 * t7 +
        0.2604684913868007E3 * t149;
      double t306 = t127 * t6;
      double t307 = t306 * t5;
      double t308 = -t3 - t22 + 0.917001814331423E-1;
      double t309 = -t3 - t22 - 0.2907007820975211;
      double t310 = t309 * t308;
      double t311 = -t3 - t22 - 0.7092992179024789;
      double t312 = -t3 - t22 - 0.1091700181433142E1;
      double t313 = t312 * t311;
      double t314 = -t3 - t22 - 0.1371740148509607E1;
      double t315 = t314 * t313;
      double t316 = t315 * t310;
      double t319 = -t3 - t22 + 0.3717401485096066;
      double t320 = t319 * t127;
      double t327 = t308 * t319;
      double t328 = t311 * t309;
      double t351 = t319 * t6 * t5;
      double t354 = -t3 - t22 + 0.3997579954114602;
      double t356 = t354 * t6 * t5;
      double t357 = -t3 - t22 + 0.1771862795107378;
      double t358 = -t3 - t22 - 0.1368825361738218;
      double t359 = t358 * t357;
      double t360 = -t3 - t22 - 0.1177186279510738E1;
      double t362 = -t3 - t22 - 0.139975799541146E1;
      double t369 = -t3 - t22 - 0.8631174638261782;
      double t371 = t362 * t360 * t369;
      double t375 =
        0.3370524498586674E4 * t316 * t307 -
        0.3370524498586674E4 * t316 * t320 * t5 +
        0.3370524498586674E4 * t316 * t320 * t48 +
        0.3370524498586674E4 * t312 * t328 * t327 * t307 -
        0.2604684913868007E3 * t146 +
        0.3370524498586674E4 * t314 * t328 * t327 * t307 +
        0.3370524498586674E4 * t314 * t312 * t309 * t327 * t307 +
        0.3370524498586674E4 * t315 * t327 * t307 +
        0.3370524498586674E4 * t315 * t309 * t319 * t307 +
        0.2131706864986174E5 * t316 * t351 -
        0.1828288356960021E4 * t362 * t360 * t62 * t359 * t356 -
        0.3273486216120494E3 * t181 * t5 -
        0.1828288356960021E4 * t371 * t359 * t356;
      double t381 = t371 * t62 * t358;
      double t388 = t357 * t354;
      double t395 = t369 * t62;
      double t404 = t18 * t6;
      double t405 = t404 * t5;
      double t414 = t222 * t52 * t21;
      double t417 = t19 * t6;
      double t418 = t417 * t5;
      double t421 = t20 * t5;
      double t424 = t20 * t48;
      double t427 =
        -0.1828288356960021E4 * t371 * t62 * t357 * t356 -
        0.1828288356960021E4 * t381 * t356 -
        0.1828288356960021E4 * t381 * t357 * t6 * t5 +
        0.1828288356960021E4 * t381 * t388 * t5 -
        0.1828288356960021E4 * t381 * t388 * t48 -
        0.1828288356960021E4 * t360 * t395 * t359 * t356 -
        0.1828288356960021E4 * t362 * t395 * t359 * t356 -
        0.1289429095266713E5 * t222 * t41 * t405 +
        0.2578858190533426E5 * t222 * t52 * t19 * t405 +
        0.2578858190533426E5 * t414 * t405 + 0.2578858190533426E5 * t414 * t418 +
        0.1289429095266713E5 * t414 * t421 - 0.1289429095266713E5 * t414 * t424;
      double t440 = t171 * t48;
      double t441 = t65 * t117;
      double t442 = t441 * t82;
      double t445 = t152 * t6;
      double t446 = t445 * t5;
      double t447 = t61 * t160;
      double t462 = t160 * t6;
      double t463 = t462 * t5;
      double t466 = t171 * t5;
      double t477 =
        -0.1289429095266713E5 * t228 * t41 * t405 -
        0.1289429095266713E5 * t232 * t41 * t405 -
        0.1289429095266713E5 * t225 * t41 * t405 -
        0.2047820986253188E5 * t442 * t440 -
        0.2047820986253188E5 * t65 * t64 * t70 * t447 * t446 -
        0.2047820986253188E5 * t441 * t447 * t446 -
        0.2047820986253188E5 * t441 * t70 * t160 * t446 +
        0.4095641972506376E5 * t442 * t446 + 0.4095641972506376E5 * t442 * t463 +
        0.2047820986253188E5 * t442 * t466 -
        0.2047820986253188E5 * t64 * t71 * t447 * t446 -
        0.2047820986253188E5 * t65 * t71 * t447 * t446;
      double t479 = t181 * t91;
      double t480 = t479 * t109;
      double t491 = t86 * t6 * t5;
      double t510 = 0.1371740148509607E1 + t8 + t1;
      double t512 = t510 * t6 * t5;
      double t513 = 0.7092992179024789 + t8 + t1;
      double t514 = 0.2907007820975211 + t8 + t1;
      double t516 = -0.917001814331423E-1 + t8 + t1;
      double t517 = -0.3717401485096066 + t8 + t1;
      double t518 = t517 * t516;
      double t519 = t518 * t514 * t513;
      double t522 = 0.1091700181433142E1 + t8 + t1;
      double t527 = t522 * t510;
      double t534 =
        -0.997395612445678E1 * t480 * t105 * t48 +
        0.1994791224891356E2 * t480 * t88 * t6 * t5 +
        0.997395612445678E1 * t480 * t105 * t5 +
        0.1994791224891356E2 * t181 * t90 * t89 * t491 +
        0.1994791224891356E2 * t479 * t89 * t491 +
        0.1994791224891356E2 * t479 * t90 * t88 * t491 +
        0.1994791224891356E2 * t480 * t491 + 0.3154041863586926E2 * t93 * t491 -
        0.2483179849261862E4 * t28 * t6 * t5 -
        0.1998987179858919E2 * t519 * t512 -
        0.1998987179858919E2 * t519 * t522 * t6 * t5 -
        0.9994935899294593E1 * t519 * t527 * t5 +
        0.9994935899294593E1 * t519 * t527 * t48;
      double t540 = t513 * t522;
      double t545 = t156 * t153 * t21;
      double t570 = t314 * t312;
      double t571 = t570 * t328;
      double t580 =
        -0.1998987179858919E2 * t518 * t514 * t522 * t512 -
        0.1998987179858919E2 * t518 * t540 * t512 -
        0.1151241728146806E5 * t545 * t424 +
        0.2302483456293612E5 * t156 * t153 * t19 * t405 +
        0.2302483456293612E5 * t545 * t405 + 0.2302483456293612E5 * t545 * t418 +
        0.1151241728146806E5 * t545 * t421 -
        0.1151241728146806E5 * t172 * t41 * t405 -
        0.1151241728146806E5 * t156 * t41 * t405 +
        0.1793877849773139E5 * t72 * t61 * t6 * t5 -
        0.1857917277797261E4 * t571 * t308 * t6 * t5 +
        0.1857917277797261E4 * t571 * t327 * t5 -
        0.1857917277797261E4 * t571 * t327 * t48;
      double t582 = t33 * t14;
      double t589 = t9 * t6 * t5;
      double t622 = t61 * t127;
      double t629 =
        -0.5034134040800346E4 * t582 * t11 * t48 -
        0.1151241728146806E5 * t201 * t41 * t405 -
        0.5034134040800346E4 * t28 * t13 * t191 * t589 +
        0.1006826808160069E5 * t33 * t191 * t589 +
        0.1006826808160069E5 * t33 * t13 * t10 * t589 +
        0.1006826808160069E5 * t582 * t589 +
        0.1006826808160069E5 * t582 * t10 * t6 * t5 -
        0.1998987179858919E2 * t516 * t514 * t540 * t512 -
        0.1998987179858919E2 * t517 * t514 * t540 * t512 +
        0.5034134040800346E4 * t582 * t11 * t5 -
        0.5034134040800346E4 * t23 * t13 * t191 * t589 +
        0.2836369924704243E4 * t72 * t622 * t48 -
        0.1857917277797261E4 * t313 * t310 * t351;
      double t656 = t23 * t6;
      double t659 = t181 * t127;
      double t664 =
        -0.1857917277797261E4 * t314 * t311 * t310 * t351 -
        0.1857917277797261E4 * t570 * t310 * t351 -
        0.1857917277797261E4 * t570 * t311 * t308 * t351 -
        0.1857917277797261E4 * t571 * t351 -
        0.2836369924704243E4 * t72 * t622 * t5 + 0.2483179849261862E4 * t33 * t5 +
        0.2836369924704243E4 * t72 * t307 + 0.2836369924704243E4 * t67 * t307 +
        0.2836369924704243E4 * t83 * t307 - 0.2483179849261862E4 * t656 * t5 +
        0.3228501933705645E3 * t659 * t5 - 0.3228501933705645E3 * t659 * t48;
      double t689 =
        0.1020941954076777E4 * t307 - 0.2041883908153555E4 * t181 * t6 * t5 +
        0.2836369924704243E4 * t122 * t307 + 0.2836369924704243E4 * t118 * t307 -
        0.1419070524843042E4 * t440 + 0.2838141049686084E4 * t446 +
        0.1419070524843042E4 * t466 + 0.2838141049686084E4 * t463 +
        0.2948754515800036E5 * t58 * t463 + 0.1474377257900018E5 * t58 * t466 -
        0.1474377257900018E5 * t58 * t440 -
        0.1474377257900018E5 * t57 * t53 * t160 * t446 +
        0.2948754515800036E5 * t58 * t446;
      double t691 = t52 * t160;
      double t709 = t181 * t171;
      double t726 = t153 * t6;
      double t730 =
        -0.1474377257900018E5 * t129 * t691 * t446 -
        0.1474377257900018E5 * t133 * t691 * t446 -
        0.1474377257900018E5 * t57 * t691 * t446 -
        0.2483179849261862E4 * t33 * t48 -
        0.4333117897525444E4 * t160 * t445 * t5 +
        0.1430971280520238E5 * t28 * t656 * t5 +
        0.1370252192622049E4 * t709 * t48 - 0.1370252192622049E4 * t709 * t5 -
        0.2740504385244097E4 * t181 * t462 * t5 -
        0.2740504385244097E4 * t181 * t445 * t5 -
        0.1094665829438066E4 * t215 * t48 -
        0.4304365439542336E4 * t155 * t62 * t6 * t5 -
        0.4304365439542336E4 * t155 * t726 * t5;
      double t750 = t28 * t23 * t127;
      double t763 =
        -0.4304365439542336E4 * t209 * t48 + 0.4304365439542336E4 * t209 * t5 -
        0.4304365439542336E4 * t62 * t726 * t5 +
        0.1094665829438066E4 * t215 * t5 +
        0.2189331658876131E4 * t21 * t417 * t5 +
        0.2189331658876131E4 * t21 * t404 * t5 +
        0.2189331658876131E4 * t19 * t404 * t5 +
        0.2262564256365844E4 * t750 * t48 - 0.2262564256365844E4 * t750 * t5 +
        0.2262564256365844E4 * t28 * t306 * t5 +
        0.2262564256365844E4 * t23 * t306 * t5 + 0.2933748971961482E3 * t5 -
        0.2933748971961482E3 * t48;
      double t767 =
        t78 + t137 + t190 + t235 + t272 + t304 + t375 + t427 + t477 + t534 +
        t580 + t629 + t664 + t689 + t730 + t763;
      return t767;
    }

    double
    ortho2_f62y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t9 = -t3 - t8 + 0.3302238962785669;
      double t10 = -t3 - t8 - 0.3115120652928579E-1;
      double t11 = t10 * t9;
      double t12 = -t3 - t8 - 0.5;
      double t13 = -t3 - t8 - 0.9688487934707142;
      double t14 = t13 * t12;
      double t15 = t14 * t11;
      double t18 = 0.1E1 * x;
      double t19 = 0.1371740148509607E1 + t18 + t1;
      double t21 = t19 * t6 * t5;
      double t22 = 0.1091700181433142E1 + t18 + t1;
      double t23 = 0.7092992179024789 + t18 + t1;
      double t24 = t23 * t22;
      double t25 = 0.2907007820975211 + t18 + t1;
      double t26 = -0.3717401485096066 + t18 + t1;
      double t31 = 0.5278640450004206E-1 + t18 + t1;
      double t32 = t31 * t6;
      double t33 = t32 * t5;
      double t35 = -t3 - t8 + 0.3997579954114602;
      double t37 = t35 * t6 * t5;
      double t38 = -t3 - t8 + 0.1771862795107378;
      double t39 = -t3 - t8 - 0.1368825361738218;
      double t40 = t39 * t38;
      double t41 = -t3 - t8 - 0.8631174638261782;
      double t42 = t41 * t12;
      double t43 = -t3 - t8 - 0.139975799541146E1;
      double t48 = -t3 - t8 - 0.1177186279510738E1;
      double t53 = t6 * t4;
      double t54 = 0.1154653670707977E1 + t18 + t1;
      double t55 = t54 * t53;
      double t56 = 0.5 + t18 + t1;
      double t57 = -0.1546536707079771 + t18 + t1;
      double t58 = t57 * t56;
      double t59 = -t3 - t8 - 0.5278640450004206E-1;
      double t60 = -t3 - t8 - 0.9472135954999579;
      double t61 = t60 * t59;
      double t62 = t61 * t58;
      double t65 = t12 * t10;
      double t66 = -t3 - t8 - 0.1330223896278567E1;
      double t67 = t66 * t13;
      double t68 = t67 * t65;
      double t71 = -t3 - t8 + 0.1546536707079771;
      double t72 = t71 * t6;
      double t73 = -t3 - t8 - 0.1154653670707977E1;
      double t77 = t6 * t2;
      double t78 = t12 * t71;
      double t79 = t73 * t78;
      double t82 = t54 * t6;
      double t83 = t82 * t5;
      double t84 = -t3 - t8 + 0.2650553239294647;
      double t85 = -t3 - t8 - 0.2147684835193549;
      double t86 = t85 * t84;
      double t87 = -t3 - t8 - 0.1265055323929465E1;
      double t88 = t87 * t86;
      double t92 = -0.917001814331423E-1 + t18 + t1;
      double t97 = 0.7852315164806451 + t18 + t1;
      double t100 = 0.2147684835193549 + t18 + t1;
      double t101 = -0.2650553239294647 + t18 + t1;
      double t102 = t101 * t100;
      double t103 = t61 * t102;
      double t106 =
        -0.125704235001869E5 * t15 * t7 -
        0.9994935899294593E1 * t26 * t25 * t24 * t21 +
        0.1419070524843042E4 * t33 -
        0.3656576713920043E4 * t43 * t42 * t40 * t37 -
        0.3656576713920043E4 * t48 * t42 * t40 * t37 +
        0.7688451173776322E4 * t62 * t55 - 0.125704235001869E5 * t68 * t7 -
        0.8608730879084671E4 * t73 * t72 * t5 - 0.4304365439542336E4 * t79 * t77 -
        0.2578858190533426E5 * t88 * t58 * t83 -
        0.9994935899294593E1 * t92 * t25 * t24 * t21 +
        0.5034134040800346E4 * t103 * t97 * t6 * t5;
      double t107 = t22 * t19;
      double t110 = t26 * t92;
      double t111 = t110 * t25 * t23;
      double t118 = 0.9472135954999579 + t18 + t1;
      double t119 = t118 * t77;
      double t121 = t60 * t59 * t31;
      double t130 = t118 * t6;
      double t131 = t130 * t5;
      double t133 = 0.1265055323929465E1 + t18 + t1;
      double t135 = t133 * t6 * t5;
      double t140 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t141 = t140 * t77;
      double t142 = -t3 - t8 - 0.7852315164806451;
      double t143 = t87 * t142;
      double t144 = t143 * t86;
      double t148 = t142 * t85;
      double t149 = t87 * t148;
      double t159 = t59 * t6;
      double t162 =
        0.9994935899294593E1 * t111 * t107 * t77 -
        0.9994935899294593E1 * t111 * t22 * t6 * t5 -
        0.111159036550155E5 * t121 * t119 + 0.111159036550155E5 * t121 * t7 +
        0.111159036550155E5 * t60 * t59 * t118 * t7 +
        0.1419070524843042E4 * t131 + 0.5034134040800346E4 * t103 * t135 +
        0.5146598031257122E4 * t144 * t141 -
        0.4497214982121002E4 * t149 * t84 * t77 -
        0.9994935899294593E1 * t111 * t21 - 0.2483179849261862E4 * t61 * t77 -
        0.4966359698523725E4 * t60 * t6 * t5 - 0.4966359698523725E4 * t159 * t5;
      double t168 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t172 = t140 * t6;
      double t173 = t172 * t5;
      double t175 = t168 * t140;
      double t180 = t100 * t133;
      double t184 = t97 * t133;
      double t188 = t56 * t54;
      double t189 = t57 * t188;
      double t192 = 0.1330223896278567E1 + t18 + t1;
      double t194 = t192 * t6 * t5;
      double t195 = 0.9688487934707142 + t18 + t1;
      double t196 = t56 * t195;
      double t197 = 0.3115120652928579E-1 + t18 + t1;
      double t202 = t168 * t58;
      double t205 = -t3 - t8 + 0.3717401485096066;
      double t207 = t205 * t6 * t5;
      double t208 = -t3 - t8 + 0.917001814331423E-1;
      double t209 = -t3 - t8 - 0.2907007820975211;
      double t210 = t209 * t208;
      double t211 = -t3 - t8 - 0.1091700181433142E1;
      double t212 = -t3 - t8 - 0.1371740148509607E1;
      double t213 = t212 * t211;
      double t217 = t133 * t53;
      double t218 = t100 * t97;
      double t219 = t168 * t101;
      double t220 = t219 * t218;
      double t223 =
        0.2483179849261862E4 * t61 * t53 - 0.1020941954076777E4 * t168 * t6 * t5 +
        0.2041883908153555E4 * t173 - 0.3228501933705645E3 * t175 * t77 +
        0.3228501933705645E3 * t175 * t53 +
        0.7816204880509401E3 * t101 * t180 * t7 +
        0.7816204880509401E3 * t101 * t184 * t7 -
        0.1094665829438066E4 * t189 * t77 +
        0.997395612445678E1 * t168 * t197 * t196 * t194 -
        0.1363118879952542E4 * t202 * t55 -
        0.3715834555594522E4 * t213 * t210 * t207 -
        0.2079447241953689E3 * t220 * t217;
      double t229 = t57 * t54;
      double t234 = t67 * t12 * t9;
      double t237 = -t3 - t8 - 0.7092992179024789;
      double t242 = t192 * t195;
      double t244 = t197 * t56;
      double t245 = -0.3302238962785669 + t18 + t1;
      double t246 = t168 * t245;
      double t247 = t246 * t244;
      double t250 = t31 * t118;
      double t254 = t9 * t31;
      double t259 = t73 * t71;
      double t270 = t101 * t218;
      double t273 =
        0.2604684913868007E3 * t141 + 0.3273486216120494E3 * t168 * t77 -
        0.1363118879952542E4 * t202 * t7 -
        0.1363118879952542E4 * t168 * t229 * t7 -
        0.125704235001869E5 * t234 * t7 -
        0.3715834555594522E4 * t212 * t237 * t210 * t207 -
        0.997395612445678E1 * t247 * t242 * t77 -
        0.2223180731003101E5 * t60 * t250 * t7 -
        0.4095641972506376E5 * t13 * t65 * t254 * t131 -
        0.2302483456293612E5 * t259 * t58 * t83 -
        0.2302483456293612E5 * t78 * t58 * t83 +
        0.5034134040800346E4 * t61 * t101 * t97 * t135 +
        0.7816204880509401E3 * t270 * t217;
      double t281 = t208 * t205;
      double t283 = t237 * t209;
      double t284 = t213 * t283;
      double t288 = t211 * t237;
      double t289 = t212 * t288;
      double t297 = t73 * t12;
      double t311 = t188 * t77;
      double t313 = t297 * t71 * t57;
      double t316 = t250 * t77;
      double t317 = t66 * t14;
      double t318 = t317 * t11;
      double t325 = t43 * t48 * t41;
      double t326 = t325 * t12 * t39;
      double t329 =
        0.3531949143188038E4 * t79 * t141 -
        0.8608730879084671E4 * t12 * t72 * t5 -
        0.1857917277797261E4 * t284 * t281 * t77 +
        0.6741048997173348E4 * t289 * t209 * t205 * t173 +
        0.6741048997173348E4 * t289 * t281 * t173 +
        0.1635573213471361E5 * t297 * t71 * t118 * t7 +
        0.5034134040800346E4 * t61 * t218 * t135 -
        0.1006826808160069E5 * t60 * t101 * t218 * t135 +
        0.7816204880509401E3 * t100 * t184 * t7 -
        0.1151241728146806E5 * t313 * t311 - 0.2047820986253188E5 * t318 * t316 -
        0.3656576713920043E4 * t326 * t38 * t6 * t5;
      double t330 = t38 * t35;
      double t337 = t71 * t140;
      double t341 = t56 * t6;
      double t342 = t341 * t5;
      double t370 = t245 * t197;
      double t371 = t370 * t196;
      double t374 = t142 * t86;
      double t378 =
        -0.1828288356960021E4 * t326 * t330 * t77 -
        0.3715834555594522E4 * t288 * t210 * t207 +
        0.7063898286376077E4 * t12 * t337 * t7 +
        0.1151241728146806E5 * t313 * t342 - 0.3656576713920043E4 * t326 * t37 -
        0.3656576713920043E4 * t325 * t12 * t38 * t37 -
        0.9994935899294593E1 * t110 * t25 * t22 * t21 +
        0.1627497198011096E5 * t144 * t7 +
        0.6741048997173348E4 * t212 * t283 * t281 * t173 +
        0.6741048997173348E4 * t212 * t211 * t209 * t281 * t173 +
        0.6741048997173348E4 * t211 * t283 * t281 * t173 +
        0.6308083727173852E2 * t371 * t194 -
        0.2578858190533426E5 * t374 * t58 * t83;
      double t384 = t188 * t53;
      double t387 = t84 * t140;
      double t391 = t87 * t85;
      double t400 = t54 * t77;
      double t407 = t149 * t84 * t57;
      double t420 =
        0.1029319606251424E5 * t143 * t85 * t140 * t7 +
        0.1151241728146806E5 * t313 * t384 +
        0.1029319606251424E5 * t143 * t387 * t7 +
        0.1029319606251424E5 * t391 * t387 * t7 +
        0.2047820986253188E5 * t318 * t33 -
        0.3271146426942721E5 * t297 * t250 * t7 -
        0.7688451173776322E4 * t62 * t400 +
        0.1029319606251424E5 * t148 * t387 * t7 +
        0.1289429095266713E5 * t407 * t384 -
        0.2079447241953689E3 * t168 * t100 * t184 * t7 +
        0.1828288356960021E4 * t326 * t330 * t53 -
        0.8621120764455331E4 * t189 * t7 - 0.2933748971961482E3 * t77;
      double t426 = t66 * t12 * t11;
      double t429 = t289 * t210;
      double t438 = t67 * t11;
      double t444 = t118 * t53;
      double t446 = t297 * t71 * t31;
      double t460 =
        -0.1419070524843042E4 * t316 - 0.9994935899294593E1 * t110 * t24 * t21 -
        0.125704235001869E5 * t426 * t7 + 0.1065853432493087E5 * t429 * t207 -
        0.1006826808160069E5 * t59 * t101 * t218 * t135 -
        0.1289429095266713E5 * t407 * t311 - 0.125704235001869E5 * t438 * t7 -
        0.1363118879952542E4 * t168 * t188 * t7 +
        0.1635573213471361E5 * t446 * t444 + 0.4304365439542336E4 * t79 * t53 +
        0.6285211750093451E4 * t68 * t9 * t53 -
        0.3271146426942721E5 * t259 * t250 * t7 -
        0.3271146426942721E5 * t78 * t250 * t7;
      double t492 = t60 * t59 * t140;
      double t500 =
        -0.2223180731003101E5 * t59 * t250 * t7 +
        0.7688451173776322E4 * t62 * t7 + 0.7688451173776322E4 * t61 * t229 * t7 +
        0.2047820986253188E5 * t318 * t131 -
        0.4095641972506376E5 * t317 * t10 * t31 * t131 -
        0.4095641972506376E5 * t317 * t254 * t131 -
        0.3715834555594522E4 * t284 * t208 * t6 * t5 +
        0.1151241728146806E5 * t313 * t83 +
        0.1151241728146806E5 * t297 * t71 * t56 * t83 -
        0.2262564256365844E4 * t492 * t53 -
        0.6285211750093451E4 * t68 * t9 * t77 - 0.8994429964242004E4 * t149 * t7;
      double t501 = t9 * t140;
      double t507 = t250 * t53;
      double t520 = t87 * t142 * t84;
      double t525 = t140 * t53;
      double t534 =
        0.2836369924704243E4 * t68 * t501 * t77 +
        0.5672739849408487E4 * t68 * t173 + 0.1474377257900018E5 * t144 * t507 -
        0.3656576713920043E4 * t325 * t40 * t37 -
        0.1849295524532771E3 * t371 * t192 * t77 +
        0.7688451173776322E4 * t61 * t188 * t7 -
        0.8994429964242004E4 * t520 * t7 - 0.3273486216120494E3 * t168 * t53 -
        0.2604684913868007E3 * t525 - 0.1474377257900018E5 * t144 * t316 +
        0.1474377257900018E5 * t144 * t33 + 0.1289429095266713E5 * t407 * t342 +
        0.1419070524843042E4 * t507;
      double t572 =
        0.1849295524532771E3 * t371 * t192 * t53 +
        0.1849295524532771E3 * t371 * t7 +
        0.1849295524532771E3 * t370 * t56 * t192 * t7 -
        0.8994429964242004E4 * t88 * t7 - 0.8994429964242004E4 * t374 * t7 +
        0.1289429095266713E5 * t407 * t83 +
        0.1289429095266713E5 * t149 * t84 * t56 * t83 +
        0.7063898286376077E4 * t73 * t12 * t140 * t7 +
        0.7063898286376077E4 * t73 * t337 * t7 +
        0.5672739849408487E4 * t426 * t173 +
        0.4497214982121002E4 * t149 * t84 * t53 +
        0.1094665829438066E4 * t189 * t53 +
        0.1849295524532771E3 * t370 * t242 * t7;
      double t580 = t133 * t77;
      double t608 =
        0.5672739849408487E4 * t234 * t173 + 0.5672739849408487E4 * t438 * t173 -
        0.9994935899294593E1 * t111 * t107 * t53 +
        0.2079447241953689E3 * t220 * t580 - 0.2079447241953689E3 * t220 * t7 +
        0.5672739849408487E4 * t15 * t173 +
        0.1849295524532771E3 * t245 * t56 * t242 * t7 -
        0.1537690234755264E5 * t60 * t57 * t188 * t7 -
        0.1537690234755264E5 * t59 * t57 * t188 * t7 -
        0.3531949143188038E4 * t79 * t525 +
        0.5034134040800346E4 * t103 * t184 * t53 -
        0.1635573213471361E5 * t446 * t119 + 0.1635573213471361E5 * t446 * t7;
      double t647 =
        -0.1370252192622049E4 * t168 * t32 * t5 -
        0.1370252192622049E4 * t168 * t130 * t5 +
        0.1474377257900018E5 * t144 * t131 + 0.2262564256365844E4 * t492 * t77 +
        0.997395612445678E1 * t247 * t195 * t6 * t5 +
        0.997395612445678E1 * t247 * t194 -
        0.4095641972506376E5 * t66 * t13 * t10 * t254 * t131 -
        0.4095641972506376E5 * t66 * t65 * t254 * t131 +
        0.1849295524532771E3 * t244 * t242 * t7 +
        0.997395612445678E1 * t246 * t197 * t195 * t194 +
        0.2933748971961482E3 * t53 + 0.1094665829438066E4 * t57 * t341 * t5;
      double t658 = t168 * t250;
      double t681 =
        0.1094665829438066E4 * t57 * t82 * t5 +
        0.8969389248865697E4 * t68 * t9 * t6 * t5 +
        0.1094665829438066E4 * t56 * t82 * t5 +
        0.1370252192622049E4 * t658 * t77 - 0.5146598031257122E4 * t144 * t525 -
        0.1246660774931152E4 * t7 + 0.1116900387235438E5 * t79 * t7 +
        0.997395612445678E1 * t247 * t242 * t53 -
        0.7816204880509401E3 * t270 * t580 + 0.7816204880509401E3 * t270 * t7 +
        0.2047820986253188E5 * t318 * t507 +
        0.4525128512731688E4 * t60 * t172 * t5 +
        0.4525128512731688E4 * t59 * t172 * t5;
      double t701 = t84 * t31;
      double t717 = t205 * t140;
      double t723 =
        0.1363118879952542E4 * t202 * t400 +
        0.1857917277797261E4 * t284 * t281 * t53 -
        0.131515791174578E4 * t102 * t184 * t7 -
        0.8666235795050887E4 * t31 * t130 * t5 +
        0.997395612445678E1 * t246 * t196 * t194 -
        0.2948754515800036E5 * t143 * t85 * t31 * t131 -
        0.2948754515800036E5 * t143 * t701 * t131 -
        0.2079447241953689E3 * t219 * t180 * t7 -
        0.2079447241953689E3 * t219 * t184 * t7 -
        0.2578858190533426E5 * t149 * t58 * t83 -
        0.2578858190533426E5 * t520 * t58 * t83 +
        0.3370524498586674E4 * t429 * t717 * t77 +
        0.6741048997173348E4 * t429 * t173;
      double t764 =
        0.715485640260119E4 * t60 * t159 * t5 +
        0.111159036550155E5 * t121 * t444 -
        0.3656576713920043E4 * t43 * t48 * t12 * t40 * t37 -
        0.8608730879084671E4 * t73 * t12 * t6 * t5 -
        0.2302483456293612E5 * t297 * t58 * t83 -
        0.1370252192622049E4 * t658 * t53 -
        0.5034134040800346E4 * t103 * t184 * t77 -
        0.2948754515800036E5 * t391 * t701 * t131 -
        0.2948754515800036E5 * t148 * t701 * t131 -
        0.2836369924704243E4 * t68 * t501 * t53 -
        0.3715834555594522E4 * t284 * t207 -
        0.3715834555594522E4 * t213 * t237 * t208 * t207 -
        0.3370524498586674E4 * t429 * t717 * t53;
      double t768 =
        t106 + t162 + t223 + t273 + t329 + t378 + t420 + t460 + t500 + t534 +
        t572 + t608 + t647 + t681 + t723 + t764;
      return t768;
    }

    // * f63 *********************************************************************

    double
    ortho2_f63 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * y;
      double t8 = -t3 - t7 - 0.5278640450004206E-1;
      double t10 = -t3 - t7 - 0.9472135954999579;
      double t14 = 0.1E1 * x;
      double t15 = 0.9472135954999579 + t14 + t1;
      double t16 = t15 * t6;
      double t17 = 0.5278640450004206E-1 + t14 + t1;
      double t19 = t17 * t16 * t5;
      double t23 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t24 = t23 * t6;
      double t27 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t31 = t6 * t5;
      double t33 = t17 * t15;
      double t41 = 0.1154653670707977E1 + t14 + t1;
      double t42 = 0.5 + t14 + t1;
      double t43 = t42 * t41;
      double t44 = -0.1546536707079771 + t14 + t1;
      double t48 = -t3 - t7 + 0.1546536707079771;
      double t49 = -t3 - t7 - 0.5;
      double t51 = -t3 - t7 - 0.1154653670707977E1;
      double t52 = t51 * t49 * t48;
      double t55 = -t3 - t7 + 0.2650553239294647;
      double t56 = -t3 - t7 - 0.2147684835193549;
      double t58 = -t3 - t7 - 0.7852315164806451;
      double t59 = -t3 - t7 - 0.1265055323929465E1;
      double t61 = t59 * t58 * t56 * t55;
      double t68 = 0.1265055323929465E1 + t14 + t1;
      double t69 = 0.7852315164806451 + t14 + t1;
      double t71 = 0.2147684835193549 + t14 + t1;
      double t72 = -0.2650553239294647 + t14 + t1;
      double t73 = t72 * t71;
      double t77 = t10 * t8;
      double t82 = t51 * t49;
      double t86 = -t3 - t7 + 0.3717401485096066;
      double t97 =
        (-t3 - t7 - 0.1371740148509607E1) * (-t3 - t7 -
               0.1091700181433142E1) * (-t3 - t7 -
                      0.7092992179024789)
        * (-t3 - t7 - 0.2907007820975211) * (-t3 - t7 + 0.917001814331423E-1);
      double t115 = t68 * t6;
      double MapleGenVar1 =
        0.3665670843231257E4 * t10 * t8 * t6 * t5 + 0.3363458299515579E4 * t19 +
        0.6193382776831578E3 * t27 * t24 * t5 + 0.5156621548470918E3 * t31 -
        0.4203367765597256E4 * t27 * t33 * t31 -
        0.549754198888656E4 * t10 * t8 * t23 * t31 +
        0.3801812734468263E4 * t44 * t43 * t31 + 0.9078373652226388E4 * t52 * t31;
      double t121 =
        MapleGenVar1 + 0.4632145837726292E4 * t61 * t31 -
        0.4552100788113374E4 * t27 * t44 * t43 * t31 +
        0.3303211643947558E4 * t73 * t69 * t68 * t31 +
        0.2221938296640474E5 * t77 * t33 * t31 -
        0.4957386626598052E4 * t82 * t48 * t23 * t31 -
        0.4859572376279952E4 * t97 * t86 * t24 * t5 + 0.186642201717685E4 * (-t3 -
                       t7 -
                       0.139975799541146E1)
        * (-t3 - t7 - 0.1177186279510738E1) * (-t3 - t7 -
                 0.8631174638261782) * t49 * (-t3 -
                      t7 -
                      0.1368825361738218)
        * (-t3 - t7 + 0.1771862795107378) * (-t3 - t7 +
               0.3997579954114602) * t6 * t5 +
        0.2873448180739235E5 * t52 * t73 * t69 * t115 * t5;
      double t122 = -t3 - t7 + 0.3302238962785669;
      double t123 = -t3 - t7 - 0.3115120652928579E-1;
      double t125 = -t3 - t7 - 0.9688487934707142;
      double t127 = -t3 - t7 - 0.1330223896278567E1;
      double t129 = t127 * t125 * t49 * t123 * t122;
      double t132 = t41 * t6;
      double t137 = t59 * t58 * t56;
      double t144 = t24 * t5;
      double t164 = t115 * t5;
      double t165 = t71 * t69;
      double t171 = t132 * t5;
      double t172 = t44 * t42;
      double t178 = (0.1330223896278567E1 + t14 + t1) * t6 * t5;
      double t180 = t42 * (0.9688487934707142 + t14 + t1);
      double t183 =
        (-0.3302238962785669 + t14 + t1) * (0.3115120652928579E-1 + t14 + t1);
      double t188 = t16 * t5;
      MapleGenVar1 =
        0.4656720864519461E5 * t129 * t19 +
        0.6085105985839932E5 * t137 * t55 * t44 * t42 * t132 * t5 -
        0.7083045568091604E3 * t27 * t6 * t5 - 0.6047447301921191E3 * t144 +
        0.9041017619025774E3 * t97 * t86 * t6 * t5 +
        0.2017737155627241E2 * (-0.3717401485096066 + t14 +
              t1) * (-0.917001814331423E-1 + t14 +
               t1) * (0.2907007820975211 + t14 +
                t1) * (0.7092992179024789 + t14 +
                 t1) * (0.1091700181433142E1 +
                  t14 +
                  t1) *
        (0.1371740148509607E1 + t14 + t1) * t6 * t5 +
        0.1716815978136522E5 * t10 * t8 * t72 * t165 * t164 +
        0.2599768958793846E5 * t52 * t172 * t171;
      double t218 =
        MapleGenVar1 - 0.1846714705976828E4 * t27 * t183 * t180 * t178 +
        0.2128758599163483E5 * t137 * t55 * t17 * t188 -
        0.1176556539100669E5 * t61 * t144 +
        0.5116126402262862E5 * t82 * t48 * t17 * t188 -
        0.2439531448743868E4 * t129 * t144 +
        0.1109737843191367E3 * t183 * t180 * t178 -
        0.284993719277665E4 * t27 * t72 * t165 * t164 +
        0.3294753734167192E5 * t77 * t172 * t171 +
        0.1037757810531379E5 * t127 * t125 * t49 * t123 * t122 * t6 * t5;
      double t219 = t121 + t218;
      return t219;
    }

    double
    ortho2_f63x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * y;
      double t9 = -t3 - t8 + 0.1546536707079771;
      double t10 = -t3 - t8 - 0.5;
      double t11 = t10 * t9;
      double t12 = -t3 - t8 - 0.1154653670707977E1;
      double t13 = t12 * t11;
      double t16 = 0.1E1 * x;
      double t17 = 0.1154653670707977E1 + t16 + t1;
      double t18 = 0.5 + t16 + t1;
      double t19 = t18 * t17;
      double t22 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t26 = -0.1546536707079771 + t16 + t1;
      double t27 = t26 * t19;
      double t30 = t26 * t18;
      double t31 = t22 * t30;
      double t34 = t6 * t2;
      double t35 = -t3 - t8 + 0.2650553239294647;
      double t37 = -t3 - t8 - 0.2147684835193549;
      double t38 = -t3 - t8 - 0.7852315164806451;
      double t39 = t38 * t37;
      double t40 = -t3 - t8 - 0.1265055323929465E1;
      double t41 = t40 * t39;
      double t44 = t37 * t35;
      double t45 = t40 * t44;
      double t49 = t40 * t38 * t35;
      double t57 = t26 * t17;
      double t61 = t38 * t44;
      double t64 = t17 * t5;
      double t67 = t17 * t34;
      double t70 =
        0.156766329821085E5 * t13 * t7 - 0.4552100788113374E4 * t22 * t19 * t7 -
        0.7197503314542898E4 * t27 * t7 - 0.4552100788113374E4 * t31 * t7 -
        0.2316072918863146E4 * t41 * t35 * t34 - 0.2316072918863146E4 * t45 * t7 -
        0.2316072918863146E4 * t49 * t7 - 0.2316072918863146E4 * t41 * t7 +
        0.2316072918863146E4 * t41 * t35 * t5 -
        0.4552100788113374E4 * t22 * t57 * t7 - 0.2316072918863146E4 * t61 * t7 -
        0.2276050394056687E4 * t31 * t64 + 0.2276050394056687E4 * t31 * t67;
      double t73 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t74 = t73 * t34;
      double t77 = t73 * t5;
      double t84 = 0.9472135954999579 + t16 + t1;
      double t85 = 0.5278640450004206E-1 + t16 + t1;
      double t86 = t85 * t84;
      double t87 = -t3 - t8 - 0.9472135954999579;
      double t91 = -t3 - t8 - 0.5278640450004206E-1;
      double t92 = t87 * t91;
      double t96 = t92 * t30;
      double t106 = t87 * t91 * t85;
      double t109 = t84 * t5;
      double t112 = t84 * t34;
      double t115 = t9 * t73;
      double t121 =
        0.2478693313299026E4 * t13 * t74 - 0.2478693313299026E4 * t13 * t77 +
        0.2478693313299026E4 * t12 * t10 * t73 * t7 -
        0.1110969148320237E5 * t87 * t86 * t7 +
        0.3294753734167192E5 * t92 * t57 * t7 + 0.3294753734167192E5 * t96 * t7 +
        0.1647376867083596E5 * t96 * t64 +
        0.2221938296640474E5 * t87 * t91 * t84 * t7 +
        0.2221938296640474E5 * t106 * t7 + 0.1110969148320237E5 * t106 * t109 -
        0.1110969148320237E5 * t106 * t112 +
        0.2478693313299026E4 * t10 * t115 * t7 - 0.1647376867083596E5 * t96 * t67;
      double t126 = t40 * t38;
      double t127 = t126 * t44;
      double t130 = 0.1265055323929465E1 + t16 + t1;
      double t131 = t130 * t34;
      double t132 = 0.7852315164806451 + t16 + t1;
      double t133 = 0.2147684835193549 + t16 + t1;
      double t134 = t133 * t132;
      double t135 = -0.2650553239294647 + t16 + t1;
      double t136 = t135 * t134;
      double t142 = t130 * t5;
      double t147 = t133 * t130;
      double t162 = -t3 - t8 + 0.3302238962785669;
      double t164 = -t3 - t8 - 0.3115120652928579E-1;
      double t165 = t10 * t164;
      double t166 = -t3 - t8 - 0.9688487934707142;
      double t167 = -t3 - t8 - 0.1330223896278567E1;
      double t168 = t167 * t166;
      double t169 = t168 * t165;
      double t175 = t132 * t130;
      double t176 = t135 * t133;
      double t180 =
        0.2478693313299026E4 * t12 * t115 * t7 + 0.372059845952307E5 * t127 * t7 -
        0.1651605821973779E4 * t136 * t131 -
        0.1110969148320237E5 * t91 * t86 * t7 +
        0.1651605821973779E4 * t136 * t142 + 0.3303211643947558E4 * t136 * t7 +
        0.3303211643947558E4 * t135 * t147 * t7 -
        0.1647376867083596E5 * t91 * t26 * t19 * t7 -
        0.1647376867083596E5 * t87 * t26 * t19 * t7 +
        0.3294753734167192E5 * t92 * t19 * t7 +
        0.5188789052656896E4 * t169 * t162 * t5 -
        0.5188789052656896E4 * t169 * t162 * t34 -
        0.4506146358800292E4 * t176 * t175 * t7;
      double t188 = t168 * t10 * t162;
      double t193 = t164 * t162;
      double t194 = t168 * t193;
      double t197 = 0.1330223896278567E1 + t16 + t1;
      double t199 = 0.9688487934707142 + t16 + t1;
      double t200 = t18 * t199;
      double t201 = 0.3115120652928579E-1 + t16 + t1;
      double t202 = -0.3302238962785669 + t16 + t1;
      double t203 = t202 * t201;
      double t204 = t203 * t200;
      double t210 = t199 * t197;
      double t211 = t201 * t18;
      double t229 = t167 * t10 * t193;
      double t232 = t10 * t166;
      double t233 = t232 * t193;
      double t236 =
        0.3303211643947558E4 * t133 * t175 * t7 +
        0.3303211643947558E4 * t135 * t175 * t7 -
        0.5188789052656896E4 * t188 * t7 - 0.5188789052656896E4 * t169 * t7 -
        0.5188789052656896E4 * t194 * t7 +
        0.5548689215956835E2 * t204 * t197 * t5 -
        0.5548689215956835E2 * t204 * t197 * t34 +
        0.1109737843191367E3 * t211 * t210 * t7 +
        0.1109737843191367E3 * t202 * t18 * t210 * t7 +
        0.1109737843191367E3 * t203 * t210 * t7 +
        0.1109737843191367E3 * t203 * t18 * t197 * t7 +
        0.1109737843191367E3 * t204 * t7 - 0.5188789052656896E4 * t229 * t7 -
        0.5188789052656896E4 * t233 * t7;
      double t243 = t35 * t73;
      double t247 = t40 * t37;
      double t259 = t12 * t10;
      double t260 = t259 * t9 * t85;
      double t270 = t12 * t9;
      double t281 =
        0.5882782695503345E4 * t127 * t74 - 0.5882782695503345E4 * t127 * t77 +
        0.5882782695503345E4 * t39 * t243 * t7 +
        0.5882782695503345E4 * t247 * t243 * t7 +
        0.5882782695503345E4 * t126 * t243 * t7 +
        0.5882782695503345E4 * t126 * t37 * t73 * t7 -
        0.2558063201131431E5 * t260 * t112 + 0.5116126402262862E5 * t260 * t7 +
        0.2558063201131431E5 * t260 * t109 -
        0.2558063201131431E5 * t11 * t86 * t7 -
        0.2558063201131431E5 * t270 * t86 * t7 -
        0.2558063201131431E5 * t259 * t86 * t7 +
        0.5116126402262862E5 * t259 * t9 * t84 * t7;
      double t282 = t22 * t135;
      double t283 = t282 * t134;
      double t308 = t130 * t6 * t5;
      double t309 = t9 * t135;
      double t318 =
        -0.1424968596388325E4 * t283 * t142 + 0.1424968596388325E4 * t283 * t131 -
        0.284993719277665E4 * t283 * t7 - 0.284993719277665E4 * t282 * t175 * t7 -
        0.284993719277665E4 * t282 * t147 * t7 -
        0.284993719277665E4 * t22 * t133 * t175 * t7 -
        0.3541522784045802E3 * t22 * t5 + 0.3023723650960595E3 * t74 +
        0.3541522784045802E3 * t22 * t34 - 0.3023723650960595E3 * t77 +
        0.7924429120945558E3 * t7 -
        0.1436724090369617E5 * t12 * t309 * t134 * t308 -
        0.1436724090369617E5 * t10 * t309 * t134 * t308;
      double t320 = t73 * t6;
      double t321 = t320 * t5;
      double t322 = -t3 - t8 + 0.917001814331423E-1;
      double t323 = -t3 - t8 - 0.2907007820975211;
      double t324 = t323 * t322;
      double t325 = -t3 - t8 - 0.7092992179024789;
      double t326 = -t3 - t8 - 0.1091700181433142E1;
      double t327 = t326 * t325;
      double t328 = -t3 - t8 - 0.1371740148509607E1;
      double t329 = t328 * t327;
      double t330 = t329 * t324;
      double t333 = -t3 - t8 + 0.3717401485096066;
      double t334 = t333 * t73;
      double t342 = t333 * t6 * t5;
      double t345 = t322 * t333;
      double t346 = t325 * t323;
      double t367 = -t3 - t8 + 0.3997579954114602;
      double t368 = -t3 - t8 + 0.1771862795107378;
      double t369 = t368 * t367;
      double t371 = -t3 - t8 - 0.1368825361738218;
      double t373 = -t3 - t8 - 0.8631174638261782;
      double t374 = -t3 - t8 - 0.1177186279510738E1;
      double t376 = -t3 - t8 - 0.139975799541146E1;
      double t377 = t376 * t374 * t373;
      double t378 = t377 * t10 * t371;
      double t381 = t175 * t34;
      double t382 = t13 * t176;
      double t392 =
        0.2429786188139976E4 * t330 * t321 -
        0.2429786188139976E4 * t330 * t334 * t5 +
        0.2429786188139976E4 * t330 * t334 * t34 +
        0.1536731716348146E5 * t330 * t342 +
        0.2429786188139976E4 * t326 * t346 * t345 * t321 +
        0.2429786188139976E4 * t328 * t346 * t345 * t321 +
        0.2429786188139976E4 * t328 * t326 * t323 * t345 * t321 +
        0.2429786188139976E4 * t329 * t345 * t321 +
        0.2429786188139976E4 * t329 * t323 * t333 * t321 -
        0.9332110085884252E3 * t378 * t369 * t34 -
        0.1436724090369617E5 * t382 * t381 -
        0.9332110085884252E3 * t378 * t368 * t6 * t5 +
        0.9332110085884252E3 * t378 * t369 * t5;
      double t394 = t132 * t6 * t5;
      double t398 = t367 * t6 * t5;
      double t399 = t371 * t368;
      double t400 = t373 * t10;
      double t417 = t175 * t5;
      double t426 = t19 * t5;
      double t428 = t41 * t35 * t26;
      double t431 = t19 * t34;
      double t434 = t17 * t6;
      double t435 = t434 * t5;
      double t442 = t18 * t6;
      double t443 = t442 * t5;
      double t449 =
        0.2873448180739235E5 * t382 * t394 -
        0.9332110085884252E3 * t374 * t400 * t399 * t398 -
        0.9332110085884252E3 * t376 * t400 * t399 * t398 -
        0.9332110085884252E3 * t376 * t374 * t10 * t399 * t398 -
        0.9332110085884252E3 * t377 * t399 * t398 +
        0.1436724090369617E5 * t382 * t417 -
        0.9332110085884252E3 * t377 * t10 * t368 * t398 -
        0.9332110085884252E3 * t378 * t398 + 0.3042552992919966E5 * t428 * t426 -
        0.3042552992919966E5 * t428 * t431 +
        0.6085105985839932E5 * t41 * t35 * t18 * t435 +
        0.6085105985839932E5 * t428 * t435 + 0.6085105985839932E5 * t428 * t443 -
        0.3042552992919966E5 * t45 * t30 * t435;
      double t459 = t135 * t132;
      double t471 = t86 * t34;
      double t472 = t167 * t232;
      double t473 = t472 * t193;
      double t476 = t85 * t6;
      double t477 = t476 * t5;
      double t480 = t86 * t5;
      double t488 = t84 * t6;
      double t489 = t488 * t5;
      double t490 = t162 * t85;
      double t504 =
        -0.3042552992919966E5 * t49 * t30 * t435 -
        0.3042552992919966E5 * t41 * t30 * t435 +
        0.2873448180739235E5 * t13 * t459 * t308 +
        0.2873448180739235E5 * t382 * t308 -
        0.3042552992919966E5 * t61 * t30 * t435 +
        0.2873448180739235E5 * t13 * t134 * t308 -
        0.232836043225973E5 * t473 * t471 + 0.4656720864519461E5 * t473 * t477 +
        0.232836043225973E5 * t473 * t480 -
        0.1436724090369617E5 * t12 * t10 * t135 * t134 * t308 -
        0.232836043225973E5 * t166 * t165 * t490 * t489 -
        0.232836043225973E5 * t167 * t165 * t490 * t489 -
        0.232836043225973E5 * t167 * t166 * t164 * t490 * t489;
      double t515 = t22 * t202;
      double t516 = t515 * t211;
      double t527 = t197 * t6 * t5;
      double t543 = 0.1371740148509607E1 + t16 + t1;
      double t544 = 0.1091700181433142E1 + t16 + t1;
      double t545 = t544 * t543;
      double t547 = 0.7092992179024789 + t16 + t1;
      double t548 = 0.2907007820975211 + t16 + t1;
      double t550 = -0.917001814331423E-1 + t16 + t1;
      double t551 = -0.3717401485096066 + t16 + t1;
      double t552 = t551 * t550;
      double t553 = t552 * t548 * t547;
      double t559 =
        -0.232836043225973E5 * t472 * t490 * t489 -
        0.232836043225973E5 * t472 * t164 * t85 * t489 +
        0.4656720864519461E5 * t473 * t489 +
        0.9233573529884142E3 * t516 * t210 * t34 -
        0.9233573529884142E3 * t516 * t210 * t5 -
        0.1846714705976828E4 * t516 * t199 * t6 * t5 -
        0.1846714705976828E4 * t516 * t527 -
        0.1846714705976828E4 * t515 * t200 * t527 -
        0.1846714705976828E4 * t515 * t201 * t199 * t527 -
        0.1846714705976828E4 * t22 * t201 * t200 * t527 -
        0.2919912329707471E4 * t204 * t527 +
        0.100886857781362E2 * t553 * t545 * t5 -
        0.100886857781362E2 * t553 * t545 * t34;
      double t562 = t543 * t6 * t5;
      double t563 = t547 * t544;
      double t578 = t259 * t9 * t26;
      double t592 = t328 * t326;
      double t593 = t592 * t346;
      double t605 =
        0.2017737155627241E2 * t552 * t563 * t562 +
        0.2017737155627241E2 * t552 * t548 * t544 * t562 +
        0.2017737155627241E2 * t553 * t562 +
        0.2017737155627241E2 * t553 * t544 * t6 * t5 +
        0.2599768958793846E5 * t578 * t443 + 0.1299884479396923E5 * t578 * t426 -
        0.1299884479396923E5 * t578 * t431 +
        0.2599768958793846E5 * t259 * t9 * t18 * t435 +
        0.2599768958793846E5 * t578 * t435 -
        0.4520508809512887E3 * t593 * t345 * t34 -
        0.4520508809512887E3 * t593 * t322 * t6 * t5 +
        0.4520508809512887E3 * t593 * t345 * t5 -
        0.1832835421615628E4 * t92 * t34;
      double t619 = t92 * t176;
      double t642 = t91 * t6;
      double t648 =
        -0.1299884479396923E5 * t11 * t30 * t435 -
        0.1299884479396923E5 * t270 * t30 * t435 -
        0.1299884479396923E5 * t259 * t30 * t435 +
        0.7714475801640934E4 * t169 * t162 * t6 * t5 -
        0.8584079890682611E4 * t619 * t381 + 0.1716815978136522E5 * t619 * t308 +
        0.1716815978136522E5 * t619 * t394 +
        0.2017737155627241E2 * t550 * t548 * t563 * t562 +
        0.2017737155627241E2 * t551 * t548 * t563 * t562 +
        0.8584079890682611E4 * t619 * t417 +
        0.1716815978136522E5 * t92 * t459 * t308 -
        0.195851959960461E4 * t22 * t6 * t5 - 0.1832835421615628E4 * t642 * t5 -
        0.1832835421615628E4 * t87 * t6 * t5;
      double t657 = t22 * t73;
      double t681 =
        0.1832835421615628E4 * t92 * t5 + 0.9792597998023052E3 * t321 +
        0.1716815978136522E5 * t92 * t134 * t308 +
        0.3096691388415789E3 * t657 * t5 - 0.3096691388415789E3 * t657 * t34 -
        0.8584079890682611E4 * t91 * t135 * t134 * t308 -
        0.8584079890682611E4 * t87 * t135 * t134 * t308 +
        0.3363458299515579E4 * t489 - 0.1681729149757789E4 * t471 +
        0.1681729149757789E4 * t480 + 0.3363458299515579E4 * t477 -
        0.4520508809512887E3 * t592 * t324 * t342 -
        0.4520508809512887E3 * t592 * t325 * t322 * t342;
      double t695 = t162 * t73;
      double t716 =
        -0.4520508809512887E3 * t593 * t342 -
        0.4520508809512887E3 * t328 * t325 * t324 * t342 -
        0.4520508809512887E3 * t327 * t324 * t342 +
        0.1219765724371934E4 * t188 * t321 + 0.1219765724371934E4 * t169 * t321 -
        0.1219765724371934E4 * t169 * t695 * t5 +
        0.1219765724371934E4 * t169 * t695 * t34 +
        0.1219765724371934E4 * t194 * t321 + 0.1219765724371934E4 * t233 * t321 +
        0.1219765724371934E4 * t229 * t321 - 0.1064379299581742E5 * t127 * t471 +
        0.1064379299581742E5 * t127 * t480 -
        0.1064379299581742E5 * t126 * t37 * t85 * t489;
      double t722 = t35 * t85;
      double t738 = t22 * t86;
      double t755 =
        0.2128758599163483E5 * t127 * t489 + 0.2128758599163483E5 * t127 * t477 -
        0.1064379299581742E5 * t247 * t722 * t489 -
        0.1064379299581742E5 * t126 * t722 * t489 -
        0.1064379299581742E5 * t39 * t722 * t489 -
        0.664610799131004E4 * t85 * t488 * t5 +
        0.1738475421729361E5 * t87 * t642 * t5 +
        0.2101683882798628E4 * t738 * t34 - 0.2101683882798628E4 * t738 * t5 -
        0.4203367765597256E4 * t22 * t476 * t5 -
        0.4203367765597256E4 * t22 * t488 * t5 -
        0.1900906367234132E4 * t27 * t34 -
        0.4539186826113194E4 * t12 * t10 * t6 * t5;
      double t758 = t9 * t6;
      double t779 = t87 * t91 * t73;
      double t792 =
        -0.4539186826113194E4 * t13 * t34 -
        0.4539186826113194E4 * t12 * t758 * t5 + 0.4539186826113194E4 * t13 * t5 -
        0.4539186826113194E4 * t10 * t758 * t5 + 0.1900906367234132E4 * t27 * t5 +
        0.3801812734468263E4 * t26 * t442 * t5 +
        0.3801812734468263E4 * t26 * t434 * t5 +
        0.3801812734468263E4 * t18 * t434 * t5 +
        0.274877099444328E4 * t779 * t34 - 0.274877099444328E4 * t779 * t5 +
        0.274877099444328E4 * t87 * t320 * t5 +
        0.274877099444328E4 * t91 * t320 * t5 - 0.2578310774235459E3 * t34 +
        0.2578310774235459E3 * t5;
      double t796 =
        t70 + t121 + t180 + t236 + t281 + t318 + t392 + t449 + t504 + t559 +
        t605 + t648 + t681 + t716 + t755 + t792;
      return t796;
    }

    double
    ortho2_f63y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * x;
      double t8 = 0.1265055323929465E1 + t7 + t1;
      double t10 = t8 * t6 * t5;
      double t11 = 0.7852315164806451 + t7 + t1;
      double t12 = 0.2147684835193549 + t7 + t1;
      double t13 = t12 * t11;
      double t14 = -0.2650553239294647 + t7 + t1;
      double t15 = 0.1E1 * y;
      double t16 = -t3 - t15 - 0.9472135954999579;
      double t21 = t6 * t5;
      double t22 = 0.9472135954999579 + t7 + t1;
      double t23 = 0.5278640450004206E-1 + t7 + t1;
      double t24 = t23 * t22;
      double t25 = -t3 - t15 + 0.1546536707079771;
      double t26 = -t3 - t15 - 0.5;
      double t27 = t26 * t25;
      double t31 = -t3 - t15 - 0.5278640450004206E-1;
      double t37 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t39 = -t3 - t15 - 0.1154653670707977E1;
      double t43 = t11 * t8;
      double t46 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t47 = t46 * t14;
      double t55 = t6 * t2;
      double t56 = t24 * t55;
      double t57 = -t3 - t15 + 0.2650553239294647;
      double t58 = -t3 - t15 - 0.2147684835193549;
      double t59 = t58 * t57;
      double t60 = -t3 - t15 - 0.7852315164806451;
      double t61 = -t3 - t15 - 0.1265055323929465E1;
      double t62 = t61 * t60;
      double t63 = t62 * t59;
      double t66 = 0.5 + t7 + t1;
      double t67 = t66 * t6;
      double t68 = -0.1546536707079771 + t7 + t1;
      double t72 = 0.1154653670707977E1 + t7 + t1;
      double t73 = t72 * t6;
      double t77 = t37 * t6;
      double t78 = t77 * t5;
      double t79 = -t3 - t15 + 0.3302238962785669;
      double t80 = -t3 - t15 - 0.3115120652928579E-1;
      double t81 = t80 * t79;
      double t82 = -t3 - t15 - 0.1330223896278567E1;
      double t84 = t82 * t26 * t81;
      double t87 = 0.9688487934707142 + t7 + t1;
      double t88 = t66 * t87;
      double t89 = 0.3115120652928579E-1 + t7 + t1;
      double t90 = -0.3302238962785669 + t7 + t1;
      double t91 = t90 * t89;
      double t92 = t91 * t88;
      double t95 = 0.1330223896278567E1 + t7 + t1;
      double t100 = -t3 - t15 + 0.3717401485096066;
      double t101 = -t3 - t15 + 0.917001814331423E-1;
      double t102 = t101 * t100;
      double t103 = -t3 - t15 - 0.2907007820975211;
      double t104 = -t3 - t15 - 0.7092992179024789;
      double t105 = t104 * t103;
      double t106 = -t3 - t15 - 0.1371740148509607E1;
      double t111 =
        -0.1716815978136522E5 * t16 * t14 * t13 * t10 -
        0.5116126402262862E5 * t27 * t24 * t21 -
        0.2221938296640474E5 * t31 * t24 * t21 +
        0.4957386626598052E4 * t39 * t26 * t37 * t21 -
        0.1424968596388325E4 * t47 * t43 * t21 -
        0.1424968596388325E4 * t46 * t12 * t43 * t21 -
        0.1064379299581742E5 * t63 * t56 + 0.1900906367234132E4 * t68 * t67 * t5 +
        0.1900906367234132E4 * t68 * t73 * t5 + 0.2439531448743868E4 * t84 * t78 +
        0.5548689215956835E2 * t92 * t21 +
        0.5548689215956835E2 * t91 * t66 * t95 * t21 +
        0.4859572376279952E4 * t106 * t105 * t102 * t78;
      double t112 = -t3 - t15 - 0.1091700181433142E1;
      double t117 = t61 * t59;
      double t120 = t60 * t59;
      double t127 = -t3 - t15 - 0.9688487934707142;
      double t128 = t127 * t26;
      double t129 = t128 * t81;
      double t132 = t6 * t4;
      double t134 = t106 * t112;
      double t135 = t134 * t105;
      double t138 = t37 * t132;
      double t145 = t8 * t55;
      double t146 = t14 * t13;
      double t149 = t25 * t37;
      double t154 = t100 * t6 * t5;
      double t157 =
        0.4859572376279952E4 * t112 * t105 * t102 * t78 -
        0.4632145837726292E4 * t117 * t21 - 0.4632145837726292E4 * t120 * t21 -
        0.9078373652226388E4 * t39 * t26 * t6 * t5 -
        0.1037757810531379E5 * t129 * t21 +
        0.4520508809512887E3 * t135 * t102 * t132 - 0.3023723650960595E3 * t138 -
        0.1283670301397526E4 * t21 + 0.2439531448743868E4 * t129 * t78 -
        0.5882782695503345E4 * t63 * t138 - 0.1651605821973779E4 * t146 * t145 +
        0.4957386626598052E4 * t39 * t149 * t21 -
        0.9041017619025774E3 * t135 * t154;
      double t162 = t39 * t27;
      double t166 = t72 * t132;
      double t167 = t68 * t66;
      double t168 = t16 * t31;
      double t169 = t168 * t167;
      double t173 = t60 * t58;
      double t174 = t61 * t173;
      double t177 = t87 * t95;
      double t181 = -t3 - t15 + 0.1771862795107378;
      double t184 = -t3 - t15 - 0.1368825361738218;
      double t186 = -t3 - t15 - 0.8631174638261782;
      double t187 = -t3 - t15 - 0.1177186279510738E1;
      double t189 = -t3 - t15 - 0.139975799541146E1;
      double t190 = t189 * t187 * t186;
      double t191 = t190 * t26 * t184;
      double t203 = t22 * t55;
      double t205 = t16 * t31 * t23;
      double t208 = t25 * t6;
      double t212 = t82 * t128;
      double t213 = t212 * t81;
      double t217 = t95 * t6 * t5;
      double t220 =
        0.1900906367234132E4 * t66 * t73 * t5 +
        0.1436724090369617E5 * t162 * t13 * t10 +
        0.1647376867083596E5 * t169 * t166 -
        0.2316072918863146E4 * t174 * t57 * t55 +
        0.5548689215956835E2 * t91 * t177 * t21 -
        0.186642201717685E4 * t191 * t181 * t6 * t5 +
        0.2316072918863146E4 * t174 * t57 * t132 +
        0.1860299229761535E5 * t63 * t21 -
        0.1716815978136522E5 * t31 * t14 * t13 * t10 -
        0.1110969148320237E5 * t205 * t203 -
        0.9078373652226388E4 * t39 * t208 * t5 -
        0.232836043225973E5 * t213 * t56 - 0.5839824659414941E4 * t92 * t217;
      double t221 = t103 * t101;
      double t222 = t112 * t104;
      double t223 = t106 * t222;
      double t224 = t223 * t221;
      double t229 = 0.1371740148509607E1 + t7 + t1;
      double t231 = t229 * t6 * t5;
      double t232 = 0.7092992179024789 + t7 + t1;
      double t233 = 0.2907007820975211 + t7 + t1;
      double t235 = -0.917001814331423E-1 + t7 + t1;
      double t236 = -0.3717401485096066 + t7 + t1;
      double t237 = t236 * t235;
      double t238 = t237 * t233 * t232;
      double t245 = t37 * t55;
      double t247 = t66 * t72;
      double t248 = t68 * t247;
      double t252 = t46 * t90;
      double t260 = t89 * t66;
      double t261 = t252 * t260;
      double t264 = t46 * t167;
      double t267 = t46 * t24;
      double t270 = t43 * t55;
      double t271 = t14 * t12;
      double t272 = t168 * t271;
      double t277 = 0.1091700181433142E1 + t7 + t1;
      double t278 = t232 * t277;
      double t283 =
        0.7683658581740728E4 * t224 * t154 - 0.1832835421615628E4 * t168 * t55 +
        0.100886857781362E2 * t238 * t231 -
        0.9041017619025774E3 * t134 * t104 * t101 * t154 +
        0.3023723650960595E3 * t245 - 0.1900906367234132E4 * t248 * t55 -
        0.9233573529884142E3 * t252 * t89 * t87 * t217 -
        0.9233573529884142E3 * t252 * t88 * t217 -
        0.9233573529884142E3 * t261 * t177 * t132 -
        0.2276050394056687E4 * t264 * t21 - 0.2101683882798628E4 * t267 * t132 -
        0.8584079890682611E4 * t272 * t270 + 0.1900906367234132E4 * t248 * t132 +
        0.100886857781362E2 * t236 * t233 * t278 * t231;
      double t286 = t24 * t132;
      double t292 = t68 * t72;
      double t296 = t100 * t37;
      double t300 = t22 * t6;
      double t304 = -t3 - t15 + 0.3997579954114602;
      double t306 = t304 * t6 * t5;
      double t318 = t16 * t31 * t37;
      double t321 = t72 * t55;
      double t326 = t162 * t271;
      double t331 =
        0.232836043225973E5 * t213 * t286 -
        0.9012292717600583E4 * t271 * t43 * t21 -
        0.2276050394056687E4 * t46 * t292 * t21 -
        0.2429786188139976E4 * t224 * t296 * t132 -
        0.1329221598262008E5 * t23 * t300 * t5 -
        0.186642201717685E4 * t191 * t306 -
        0.186642201717685E4 * t190 * t26 * t181 * t306 +
        0.5548689215956835E2 * t90 * t66 * t177 * t21 -
        0.274877099444328E4 * t318 * t132 - 0.1647376867083596E5 * t169 * t321 +
        0.1647376867083596E5 * t169 * t21 + 0.1436724090369617E5 * t326 * t10 +
        0.274877099444328E4 * t318 * t55;
      double t332 = t247 * t55;
      double t334 = t39 * t26;
      double t335 = t334 * t25 * t68;
      double t339 = t174 * t57 * t68;
      double t342 = t67 * t5;
      double t345 = t73 * t5;
      double t357 = t12 * t8;
      double t361 = t23 * t6;
      double t362 = t361 * t5;
      double t371 = t8 * t132;
      double t372 = t47 * t13;
      double t379 =
        -0.1299884479396923E5 * t335 * t332 - 0.3042552992919966E5 * t339 * t332 +
        0.3042552992919966E5 * t339 * t342 + 0.3042552992919966E5 * t339 * t345 -
        0.9041017619025774E3 * t134 * t221 * t154 -
        0.9041017619025774E3 * t106 * t104 * t221 * t154 +
        0.1651605821973779E4 * t146 * t21 +
        0.1651605821973779E4 * t14 * t357 * t21 +
        0.232836043225973E5 * t213 * t362 +
        0.1651605821973779E4 * t14 * t43 * t21 +
        0.1651605821973779E4 * t12 * t43 * t21 -
        0.1424968596388325E4 * t372 * t371 +
        0.100886857781362E2 * t235 * t233 * t278 * t231;
      double t382 = t11 * t6 * t5;
      double t388 = t26 * t80;
      double t389 = t82 * t127;
      double t390 = t389 * t388;
      double t401 = t31 * t6;
      double t410 = t247 * t132;
      double t413 = t79 * t37;
      double t424 =
        0.8584079890682611E4 * t272 * t382 - 0.4539186826113194E4 * t162 * t55 +
        0.5188789052656896E4 * t390 * t79 * t132 -
        0.9041017619025774E3 * t222 * t221 * t154 +
        0.5548689215956835E2 * t260 * t177 * t21 +
        0.2101683882798628E4 * t267 * t55 +
        0.8692377108646805E4 * t16 * t401 * t5 +
        0.7838316491054252E4 * t162 * t21 +
        0.4957386626598052E4 * t26 * t149 * t21 +
        0.3042552992919966E5 * t339 * t410 +
        0.1219765724371934E4 * t390 * t413 * t55 +
        0.3042552992919966E5 * t174 * t57 * t66 * t345 -
        0.6085105985839932E5 * t174 * t167 * t345;
      double t439 = t300 * t5;
      double t446 = t181 * t304;
      double t450 = t277 * t229;
      double t458 = t25 * t14;
      double t465 =
        0.8584079890682611E4 * t272 * t10 +
        0.1647376867083596E5 * t168 * t292 * t21 +
        0.1647376867083596E5 * t168 * t247 * t21 - 0.2578310774235459E3 * t55 +
        0.1832835421615628E4 * t168 * t132 -
        0.5188789052656896E4 * t390 * t79 * t55 +
        0.232836043225973E5 * t213 * t439 -
        0.4656720864519461E5 * t212 * t80 * t23 * t439 +
        0.9332110085884252E3 * t191 * t446 * t132 +
        0.100886857781362E2 * t238 * t450 * t132 +
        0.1424968596388325E4 * t372 * t145 + 0.1064379299581742E5 * t63 * t362 -
        0.2873448180739235E5 * t26 * t458 * t13 * t10 -
        0.1436724090369617E5 * t326 * t270;
      double t477 = t57 * t23;
      double t508 =
        -0.1424968596388325E4 * t372 * t21 + 0.1064379299581742E5 * t63 * t439 -
        0.2128758599163483E5 * t62 * t58 * t23 * t439 -
        0.2128758599163483E5 * t62 * t477 * t439 +
        0.549754198888656E4 * t16 * t77 * t5 +
        0.549754198888656E4 * t31 * t77 * t5 +
        0.1299884479396923E5 * t335 * t342 + 0.1299884479396923E5 * t335 * t345 -
        0.2873448180739235E5 * t39 * t458 * t13 * t10 +
        0.1110969148320237E5 * t205 * t21 +
        0.1110969148320237E5 * t16 * t31 * t22 * t21 -
        0.2221938296640474E5 * t16 * t24 * t21 -
        0.9233573529884142E3 * t46 * t89 * t88 * t217;
      double t519 = t43 * t132;
      double t534 = t61 * t58;
      double t550 =
        -0.9078373652226388E4 * t26 * t208 * t5 -
        0.4520508809512887E3 * t135 * t102 * t55 +
        0.100886857781362E2 * t237 * t233 * t277 * t231 +
        0.1436724090369617E5 * t326 * t519 -
        0.3294753734167192E5 * t16 * t68 * t247 * t21 -
        0.3294753734167192E5 * t31 * t68 * t247 * t21 +
        0.1299884479396923E5 * t334 * t25 * t66 * t345 -
        0.2128758599163483E5 * t534 * t477 * t439 +
        0.100886857781362E2 * t237 * t278 * t231 -
        0.1037757810531379E5 * t84 * t21 +
        0.2429786188139976E4 * t224 * t296 * t55 +
        0.4859572376279952E4 * t224 * t78 - 0.3541522784045802E3 * t46 * t132;
      double t558 = t79 * t23;
      double t564 = t389 * t26 * t79;
      double t570 = t61 * t60 * t57;
      double t583 = t184 * t181;
      double t592 =
        -0.4632145837726292E4 * t174 * t21 + 0.2439531448743868E4 * t390 * t78 -
        0.1037757810531379E5 * t390 * t21 -
        0.4656720864519461E5 * t212 * t558 * t439 + 0.2578310774235459E3 * t132 +
        0.2439531448743868E4 * t564 * t78 + 0.1064379299581742E5 * t63 * t286 -
        0.6085105985839932E5 * t570 * t167 * t345 -
        0.6085105985839932E5 * t117 * t167 * t345 -
        0.6085105985839932E5 * t120 * t167 * t345 -
        0.1424968596388325E4 * t47 * t357 * t21 -
        0.186642201717685E4 * t190 * t583 * t306 -
        0.186642201717685E4 * t189 * t187 * t26 * t583 * t306;
      double t593 = t22 * t132;
      double t599 = t389 * t81;
      double t628 = t334 * t25 * t23;
      double t636 = t39 * t25;
      double t640 =
        0.1110969148320237E5 * t205 * t593 +
        0.9233573529884142E3 * t261 * t177 * t55 +
        0.2439531448743868E4 * t599 * t78 - 0.2478693313299026E4 * t162 * t138 -
        0.4656720864519461E5 * t82 * t127 * t80 * t558 * t439 -
        0.4656720864519461E5 * t82 * t388 * t558 * t439 -
        0.4656720864519461E5 * t127 * t388 * t558 * t439 -
        0.9041017619025774E3 * t135 * t101 * t6 * t5 +
        0.3857237900820467E4 * t390 * t79 * t6 * t5 -
        0.4632145837726292E4 * t570 * t21 + 0.2558063201131431E5 * t628 * t593 +
        0.8584079890682611E4 * t272 * t519 -
        0.2599768958793846E5 * t334 * t167 * t345 -
        0.2599768958793846E5 * t636 * t167 * t345;
      double t647 = t186 * t26;
      double t658 = t14 * t11;
      double t688 =
        0.4859572376279952E4 * t223 * t103 * t100 * t78 -
        0.186642201717685E4 * t189 * t647 * t583 * t306 -
        0.186642201717685E4 * t187 * t647 * t583 * t306 +
        0.5882782695503345E4 * t63 * t245 +
        0.8584079890682611E4 * t168 * t658 * t10 +
        0.2558063201131431E5 * t334 * t25 * t22 * t21 -
        0.5116126402262862E5 * t334 * t24 * t21 +
        0.4859572376279952E4 * t223 * t102 * t78 +
        0.4859572376279952E4 * t106 * t112 * t103 * t102 * t78 +
        0.8584079890682611E4 * t168 * t13 * t10 +
        0.1176556539100669E5 * t62 * t58 * t37 * t21 +
        0.2276050394056687E4 * t264 * t321 + 0.1651605821973779E4 * t146 * t371;
      double t691 = t57 * t37;
      double t724 =
        -0.1037757810531379E5 * t564 * t21 +
        0.1176556539100669E5 * t62 * t691 * t21 + 0.1681729149757789E4 * t286 -
        0.100886857781362E2 * t238 * t450 * t55 +
        0.100886857781362E2 * t238 * t277 * t6 * t5 -
        0.2558063201131431E5 * t628 * t203 + 0.2558063201131431E5 * t628 * t21 -
        0.2101683882798628E4 * t46 * t361 * t5 -
        0.2599768958793846E5 * t27 * t167 * t345 +
        0.1299884479396923E5 * t335 * t410 +
        0.1176556539100669E5 * t534 * t691 * t21 -
        0.9233573529884142E3 * t261 * t87 * t6 * t5 -
        0.9233573529884142E3 * t261 * t217;
      double t745 = t46 * t37;
      double t753 =
        0.3541522784045802E3 * t46 * t55 - 0.3665670843231257E4 * t16 * t6 * t5 -
        0.9332110085884252E3 * t191 * t446 * t55 +
        0.4539186826113194E4 * t162 * t132 -
        0.2128758599163483E5 * t173 * t477 * t439 -
        0.3665670843231257E4 * t401 * t5 - 0.9792597998023052E3 * t46 * t6 * t5 +
        0.195851959960461E4 * t78 - 0.3096691388415789E3 * t745 * t55 +
        0.3096691388415789E3 * t745 * t132 - 0.1681729149757789E4 * t56 +
        0.1681729149757789E4 * t439 + 0.1681729149757789E4 * t362;
      double t793 =
        0.1436724090369617E5 * t326 * t382 -
        0.2873448180739235E5 * t39 * t26 * t14 * t13 * t10 -
        0.1219765724371934E4 * t390 * t413 * t132 -
        0.143950066290858E5 * t248 * t21 -
        0.5116126402262862E5 * t636 * t24 * t21 -
        0.1037757810531379E5 * t599 * t21 -
        0.2276050394056687E4 * t46 * t247 * t21 +
        0.1176556539100669E5 * t173 * t691 * t21 +
        0.5548689215956835E2 * t92 * t95 * t132 +
        0.1436724090369617E5 * t162 * t658 * t10 -
        0.2276050394056687E4 * t264 * t166 + 0.2478693313299026E4 * t162 * t245 -
        0.5548689215956835E2 * t92 * t95 * t55 -
        0.2101683882798628E4 * t46 * t300 * t5;
      double t797 =
        t111 + t157 + t220 + t283 + t331 + t379 + t424 + t465 + t508 + t550 +
        t592 + t640 + t688 + t724 + t753 + t793;
      return t797;
    }

    // * f64 *********************************************************************

    double
    ortho2_f64 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t9 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t10 = t9 * t6;
      double t13 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t17 = 0.1E1 * y;
      double t18 = -t3 - t17 - 0.5278640450004206E-1;
      double t20 = -t3 - t17 - 0.9472135954999579;
      double t24 = 0.1E1 * x;
      double t25 = 0.9472135954999579 + t24 + t1;
      double t26 = t25 * t6;
      double t27 = 0.5278640450004206E-1 + t24 + t1;
      double t29 = t27 * t26 * t5;
      double t31 = t6 * t5;
      double t33 = -t3 - t17 + 0.2650553239294647;
      double t34 = -t3 - t17 - 0.2147684835193549;
      double t36 = -t3 - t17 - 0.7852315164806451;
      double t37 = -t3 - t17 - 0.1265055323929465E1;
      double t39 = t37 * t36 * t34 * t33;
      double t42 = 0.1265055323929465E1 + t24 + t1;
      double t43 = 0.7852315164806451 + t24 + t1;
      double t45 = 0.2147684835193549 + t24 + t1;
      double t46 = -0.2650553239294647 + t24 + t1;
      double t47 = t46 * t45;
      double t51 = 0.1154653670707977E1 + t24 + t1;
      double t52 = 0.5 + t24 + t1;
      double t53 = t52 * t51;
      double t54 = -0.1546536707079771 + t24 + t1;
      double t59 = t27 * t25;
      double t60 = t20 * t18;
      double t64 = -t3 - t17 + 0.1546536707079771;
      double t66 = -t3 - t17 - 0.5;
      double t67 = -t3 - t17 - 0.1154653670707977E1;
      double t68 = t67 * t66;
      double t80 = t67 * t66 * t64;
      double t86 = t42 * t6;
      double t93 = (0.1330223896278567E1 + t24 + t1) * t6;
      double t94 = 0.9688487934707142 + t24 + t1;
      double t97 = 0.3115120652928579E-1 + t24 + t1;
      double t99 = -0.3302238962785669 + t24 + t1;
      double t120 = -t3 - t17 + 0.3717401485096066;
      double t131 =
        (-t3 - t17 - 0.1371740148509607E1) * (-t3 - t17 -
                0.1091700181433142E1) * (-t3 - t17 -
                       0.7092992179024789)
        * (-t3 - t17 - 0.2907007820975211) * (-t3 - t17 + 0.917001814331423E-1);
      double MapleGenVar1 =
        0.1343458479500793E3 * t13 * t10 * t5 +
        0.134693390869038E3 * t20 * t18 * t6 * t5 + 0.1424330719807454E4 * t29 +
        0.122265673540252E3 * t31 - 0.8585939423551394E3 * t39 * t31 +
        0.3415937934684501E4 * t47 * t43 * t42 * t31 -
        0.2203551587560495E4 * t13 * t54 * t53 * t31 +
        0.3706919888485688E4 * t60 * t59 * t31;
      double t134 =
        MapleGenVar1 + 0.1605221334142759E3 * t68 * t64 * t9 * t31 -
        0.5755979150724977E4 * t20 * t18 * t9 * t31 -
        0.5373181954143113E4 * t13 * t59 * t31 + 0.790303613912945E4 * t80 * t31 +
        0.6312762553650661E4 * t54 * t53 * t31 +
        0.6594944612030961E5 * t80 * t47 * t43 * t86 * t5 +
        0.2500586384583724E5 * t20 * t18 * t99 * t97 * t52 * t94 * t93 * t5 +
        0.6227107348316808E3 * (-t3 - t17 - 0.139975799541146E1) * (-t3 - t17 -
                    0.1177186279510738E1)
        * (-t3 - t17 - 0.8631174638261782) * t66 * (-t3 - t17 -
                0.1368825361738218) * (-t3 -
                           t17 +
                           0.1771862795107378)
        * (-t3 - t17 + 0.3997579954114602) * t6 * t5 -
        0.2211504185494865E4 * t131 * t120 * t10 * t5;
      double t135 = -t3 - t17 + 0.3302238962785669;
      double t136 = -t3 - t17 - 0.3115120652928579E-1;
      double t138 = -t3 - t17 - 0.9688487934707142;
      double t140 = -t3 - t17 - 0.1330223896278567E1;
      double t142 = t140 * t138 * t66 * t136 * t135;
      double t145 = t51 * t6;
      double t150 = t37 * t36 * t34;
      double t154 = t93 * t5;
      double t155 = t52 * t94;
      double t156 = t99 * t97;
      double t160 = t26 * t5;
      double t165 = t10 * t5;
      double t168 = t86 * t5;
      double t169 = t45 * t43;
      double t181 = t145 * t5;
      double t182 = t54 * t52;
      MapleGenVar1 =
        0.3099229373788809E5 * t142 * t29 +
        0.6566894860589722E5 * t150 * t33 * t54 * t52 * t145 * t5 +
        0.3194277810660077E4 * t156 * t155 * t154 +
        0.5699405234470576E5 * t68 * t64 * t27 * t160 -
        0.9168362198145701E4 * t39 * t165 -
        0.7209656522802525E4 * t13 * t46 * t169 * t168 +
        0.5834536221727242E4 * t140 * t138 * t66 * t136 * t135 * t6 * t5 +
        0.5297503343516344E5 * t60 * t182 * t181;
      double t226 =
        MapleGenVar1 - 0.705276823212651E3 * t13 * t6 * t5 -
        0.6624570099446339E3 * t165 -
        0.4566291315721297E3 * t131 * t120 * t6 * t5 +
        0.1748873503842613E4 * (-0.3717401485096066 + t24 +
              t1) * (-0.917001814331423E-1 + t24 +
               t1) * (0.2907007820975211 + t24 +
                t1) * (0.7092992179024789 + t24 +
                 t1) * (0.1091700181433142E1 +
                  t24 +
                  t1) *
        (0.1371740148509607E1 + t24 + t1) * t6 * t5 -
        0.2196498431742222E4 * t13 * t156 * t155 * t154 -
        0.1535661541073266E4 * t150 * t33 * t27 * t160 +
        0.6399173559048418E3 * t142 * t165 +
        0.8417870238879249E4 * t20 * t18 * t46 * t169 * t168 +
        0.3863456908359779E4 * t80 * t182 * t181;
      double t227 = t134 + t226;
      return t227;
    }

    double
    ortho2_f64x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = 0.5 + t3;
      double t5 = t4 * t2;
      double t6 = 0.1E1 * x;
      double t7 = 0.1154653670707977E1 + t6 + t1;
      double t8 = t7 * t5;
      double t9 = 0.5 + t6 + t1;
      double t10 = -0.1546536707079771 + t6 + t1;
      double t11 = t10 * t9;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 - 0.5278640450004206E-1;
      double t14 = -t3 - t12 - 0.9472135954999579;
      double t15 = t14 * t13;
      double t16 = t15 * t11;
      double t20 = (-t3 - t1) * t2;
      double t21 = t4 * t20;
      double t22 = -t3 - t12 + 0.2650553239294647;
      double t23 = -t3 - t12 - 0.2147684835193549;
      double t24 = t23 * t22;
      double t25 = -t3 - t12 - 0.7852315164806451;
      double t26 = -t3 - t12 - 0.1265055323929465E1;
      double t27 = t26 * t25;
      double t28 = t27 * t24;
      double t31 = t7 * t20;
      double t36 = t10 * t7;
      double t40 = t9 * t7;
      double t48 = -t3 - t12 + 0.3302238962785669;
      double t50 = -t3 - t12 - 0.3115120652928579E-1;
      double t51 = -t3 - t12 - 0.5;
      double t52 = t51 * t50;
      double t53 = -t3 - t12 - 0.9688487934707142;
      double t54 = -t3 - t12 - 0.1330223896278567E1;
      double t55 = t54 * t53;
      double t56 = t55 * t52;
      double t62 = 0.1265055323929465E1 + t6 + t1;
      double t63 = 0.7852315164806451 + t6 + t1;
      double t64 = t63 * t62;
      double t65 = 0.2147684835193549 + t6 + t1;
      double t66 = -0.2650553239294647 + t6 + t1;
      double t67 = t66 * t65;
      double t76 = t55 * t51 * t48;
      double t81 =
        -0.2648751671758172E5 * t16 * t8 + 0.2899290695952841E5 * t28 * t21 +
        0.2648751671758172E5 * t16 * t31 + 0.5297503343516344E5 * t16 * t21 +
        0.5297503343516344E5 * t15 * t36 * t21 -
        0.2648751671758172E5 * t14 * t10 * t40 * t21 +
        0.5297503343516344E5 * t15 * t40 * t21 +
        0.2917268110863621E4 * t56 * t48 * t20 -
        0.2917268110863621E4 * t56 * t48 * t5 -
        0.1139946787977283E5 * t67 * t64 * t21 -
        0.2648751671758172E5 * t13 * t10 * t40 * t21 -
        0.2917268110863621E4 * t76 * t21 - 0.2917268110863621E4 * t56 * t21;
      double t82 = t50 * t48;
      double t83 = t55 * t82;
      double t86 = 0.1330223896278567E1 + t6 + t1;
      double t87 = 0.9688487934707142 + t6 + t1;
      double t88 = t87 * t86;
      double t89 = -0.3302238962785669 + t6 + t1;
      double t94 = 0.3115120652928579E-1 + t6 + t1;
      double t95 = t89 * t94;
      double t103 = t9 * t87;
      double t104 = t95 * t103;
      double t113 = t94 * t9;
      double t117 = t53 * t51;
      double t118 = t117 * t82;
      double t122 = t54 * t51 * t82;
      double t127 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t128 = t127 * t5;
      double t131 = t22 * t127;
      double t139 = t127 * t20;
      double t142 =
        -0.2917268110863621E4 * t83 * t21 +
        0.3194277810660077E4 * t89 * t9 * t88 * t21 +
        0.3194277810660077E4 * t95 * t88 * t21 +
        0.3194277810660077E4 * t95 * t9 * t86 * t21 +
        0.3194277810660077E4 * t104 * t21 +
        0.1597138905330039E4 * t104 * t86 * t20 -
        0.1597138905330039E4 * t104 * t86 * t5 +
        0.3194277810660077E4 * t113 * t88 * t21 -
        0.2917268110863621E4 * t118 * t21 - 0.2917268110863621E4 * t122 * t21 +
        0.458418109907285E4 * t28 * t128 +
        0.458418109907285E4 * t27 * t131 * t21 +
        0.458418109907285E4 * t27 * t23 * t127 * t21 -
        0.458418109907285E4 * t28 * t139;
      double t144 = t25 * t23;
      double t148 = t26 * t23;
      double t152 = 0.9472135954999579 + t6 + t1;
      double t153 = t152 * t20;
      double t154 = 0.5278640450004206E-1 + t6 + t1;
      double t155 = -t3 - t12 + 0.1546536707079771;
      double t157 = -t3 - t12 - 0.1154653670707977E1;
      double t158 = t157 * t51;
      double t159 = t158 * t155 * t154;
      double t162 = t152 * t5;
      double t171 = t154 * t152;
      double t172 = t51 * t155;
      double t176 = t157 * t155;
      double t183 = t62 * t5;
      double t184 = t65 * t63;
      double t187 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t188 = t187 * t66;
      double t189 = t188 * t184;
      double t192 = t65 * t62;
      double t198 = t62 * t20;
      double t201 = t10 * t40;
      double t204 =
        0.458418109907285E4 * t144 * t131 * t21 +
        0.458418109907285E4 * t148 * t131 * t21 +
        0.2849702617235288E5 * t159 * t153 - 0.2849702617235288E5 * t159 * t162 +
        0.5699405234470576E5 * t159 * t21 +
        0.5699405234470576E5 * t158 * t155 * t152 * t21 -
        0.2849702617235288E5 * t172 * t171 * t21 -
        0.2849702617235288E5 * t176 * t171 * t21 -
        0.2849702617235288E5 * t158 * t171 * t21 +
        0.3604828261401262E4 * t189 * t183 -
        0.7209656522802525E4 * t188 * t192 * t21 -
        0.7209656522802525E4 * t189 * t21 - 0.3604828261401262E4 * t189 * t198 -
        0.348412097918556E4 * t201 * t21;
      double t215 = t157 * t172;
      double t219 = t26 * t144;
      double t222 = t187 * t11;
      double t225 = t26 * t24;
      double t229 = t26 * t25 * t22;
      double t240 = t25 * t24;
      double t247 =
        -0.2203551587560495E4 * t187 * t40 * t21 -
        0.7209656522802525E4 * t188 * t64 * t21 -
        0.7209656522802525E4 * t187 * t65 * t64 * t21 -
        0.5076155564585329E3 * t215 * t21 +
        0.4292969711775697E3 * t219 * t22 * t5 -
        0.2203551587560495E4 * t222 * t21 + 0.4292969711775697E3 * t225 * t21 +
        0.4292969711775697E3 * t229 * t21 + 0.4292969711775697E3 * t219 * t21 -
        0.4292969711775697E3 * t219 * t22 * t20 -
        0.2203551587560495E4 * t187 * t36 * t21 +
        0.4292969711775697E3 * t240 * t21 + 0.1101775793780247E4 * t222 * t8 -
        0.8026106670713796E2 * t215 * t128;
      double t252 = t155 * t127;
      double t266 = t14 * t13 * t154;
      double t282 = t66 * t184;
      double t287 =
        -0.1101775793780247E4 * t222 * t31 -
        0.8026106670713796E2 * t51 * t252 * t21 -
        0.8026106670713796E2 * t157 * t252 * t21 -
        0.8026106670713796E2 * t157 * t51 * t127 * t21 +
        0.8026106670713796E2 * t215 * t139 + 0.1853459944242844E4 * t266 * t153 -
        0.1853459944242844E4 * t266 * t162 + 0.3706919888485688E4 * t266 * t21 +
        0.3415937934684501E4 * t65 * t64 * t21 +
        0.3415937934684501E4 * t66 * t64 * t21 +
        0.3415937934684501E4 * t66 * t192 * t21 +
        0.3415937934684501E4 * t282 * t21 + 0.170796896734225E4 * t282 * t198;
      double t308 = t62 * t4 * t20;
      double t309 = t155 * t66;
      double t314 = t88 * t5;
      double t316 = t14 * t13 * t89;
      double t317 = t316 * t113;
      double t320 = t88 * t20;
      double t324 = t86 * t4 * t20;
      double t325 = t94 * t87;
      double t331 =
        -0.170796896734225E4 * t282 * t183 -
        0.1853459944242844E4 * t13 * t171 * t21 -
        0.1853459944242844E4 * t14 * t171 * t21 +
        0.3706919888485688E4 * t14 * t13 * t152 * t21 -
        0.331228504972317E3 * t139 - 0.3526384116063255E3 * t187 * t20 +
        0.3526384116063255E3 * t187 * t5 + 0.9797324322299129E3 * t21 +
        0.331228504972317E3 * t128 -
        0.329747230601548E5 * t157 * t309 * t184 * t308 -
        0.1250293192291862E5 * t317 * t314 + 0.1250293192291862E5 * t317 * t320 +
        0.2500586384583724E5 * t316 * t325 * t324 +
        0.2500586384583724E5 * t317 * t324;
      double t334 = t87 * t4 * t20;
      double t353 = t127 * t4;
      double t354 = t353 * t20;
      double t355 = -t3 - t12 + 0.3717401485096066;
      double t356 = -t3 - t12 + 0.917001814331423E-1;
      double t357 = t356 * t355;
      double t358 = -t3 - t12 - 0.7092992179024789;
      double t359 = -t3 - t12 - 0.1091700181433142E1;
      double t360 = t359 * t358;
      double t361 = -t3 - t12 - 0.1371740148509607E1;
      double t362 = t361 * t360;
      double t366 = -t3 - t12 - 0.2907007820975211;
      double t371 = t366 * t356;
      double t372 = t362 * t371;
      double t375 = t355 * t127;
      double t382 = -t3 - t12 + 0.3997579954114602;
      double t383 = -t3 - t12 + 0.1771862795107378;
      double t384 = t383 * t382;
      double t386 = -t3 - t12 - 0.1368825361738218;
      double t388 = -t3 - t12 - 0.8631174638261782;
      double t389 = -t3 - t12 - 0.1177186279510738E1;
      double t391 = -t3 - t12 - 0.139975799541146E1;
      double t392 = t391 * t389 * t388;
      double t393 = t392 * t51 * t386;
      double t399 = t64 * t5;
      double t400 = t215 * t67;
      double t404 = t355 * t4 * t20;
      double t407 =
        0.2500586384583724E5 * t317 * t334 +
        0.2500586384583724E5 * t316 * t103 * t324 -
        0.1250293192291862E5 * t13 * t95 * t103 * t324 -
        0.1250293192291862E5 * t14 * t95 * t103 * t324 +
        0.2500586384583724E5 * t14 * t13 * t94 * t103 * t324 +
        0.1105752092747433E4 * t362 * t357 * t354 +
        0.1105752092747433E4 * t362 * t366 * t355 * t354 +
        0.1105752092747433E4 * t372 * t354 -
        0.1105752092747433E4 * t372 * t375 * t20 +
        0.1105752092747433E4 * t372 * t375 * t5 +
        0.3113553674158404E3 * t393 * t384 * t20 -
        0.3113553674158404E3 * t393 * t384 * t5 -
        0.329747230601548E5 * t400 * t399 + 0.6993390281159279E4 * t372 * t404;
      double t408 = t358 * t366;
      double t423 = t382 * t4 * t20;
      double t435 = t63 * t4 * t20;
      double t438 = t386 * t383;
      double t439 = t388 * t51;
      double t456 = t64 * t20;
      double t459 = t40 * t5;
      double t461 = t219 * t22 * t10;
      double t464 = t7 * t4;
      double t465 = t464 * t20;
      double t468 =
        0.1105752092747433E4 * t359 * t408 * t357 * t354 +
        0.1105752092747433E4 * t361 * t408 * t357 * t354 +
        0.1105752092747433E4 * t361 * t359 * t366 * t357 * t354 -
        0.3113553674158404E3 * t392 * t51 * t383 * t423 -
        0.3113553674158404E3 * t393 * t423 -
        0.3113553674158404E3 * t393 * t383 * t4 * t20 +
        0.6594944612030961E5 * t400 * t435 -
        0.3113553674158404E3 * t389 * t439 * t438 * t423 -
        0.3113553674158404E3 * t391 * t439 * t438 * t423 -
        0.3113553674158404E3 * t391 * t389 * t51 * t438 * t423 -
        0.3113553674158404E3 * t392 * t438 * t423 +
        0.329747230601548E5 * t400 * t456 - 0.3283447430294861E5 * t461 * t459 +
        0.6566894860589722E5 * t461 * t465;
      double t472 = t9 * t4;
      double t473 = t472 * t20;
      double t476 = t40 * t20;
      double t497 = t66 * t63;
      double t501 = t187 * t89;
      double t502 = t501 * t113;
      double t512 =
        0.6566894860589722E5 * t461 * t473 + 0.3283447430294861E5 * t461 * t476 +
        0.6566894860589722E5 * t219 * t22 * t9 * t465 +
        0.6594944612030961E5 * t400 * t308 -
        0.3283447430294861E5 * t240 * t11 * t465 -
        0.3283447430294861E5 * t225 * t11 * t465 -
        0.3283447430294861E5 * t229 * t11 * t465 -
        0.3283447430294861E5 * t219 * t11 * t465 +
        0.6594944612030961E5 * t215 * t497 * t308 -
        0.2196498431742222E4 * t502 * t334 - 0.1098249215871111E4 * t502 * t320 +
        0.1098249215871111E4 * t502 * t314 -
        0.2196498431742222E4 * t501 * t103 * t324;
      double t522 = t154 * t4;
      double t523 = t522 * t20;
      double t524 = t54 * t117;
      double t525 = t524 * t82;
      double t528 = t171 * t20;
      double t533 = t171 * t5;
      double t539 = t152 * t4;
      double t540 = t539 * t20;
      double t543 = t48 * t154;
      double t564 =
        -0.2196498431742222E4 * t501 * t325 * t324 -
        0.2196498431742222E4 * t502 * t324 -
        0.2196498431742222E4 * t187 * t94 * t103 * t324 +
        0.3099229373788809E5 * t525 * t523 + 0.1549614686894405E5 * t525 * t528 -
        0.3472968960646655E4 * t104 * t324 - 0.1549614686894405E5 * t525 * t533 +
        0.6594944612030961E5 * t215 * t184 * t308 +
        0.3099229373788809E5 * t525 * t540 -
        0.1549614686894405E5 * t53 * t52 * t543 * t540 -
        0.1549614686894405E5 * t54 * t52 * t543 * t540 -
        0.1549614686894405E5 * t54 * t53 * t50 * t543 * t540 -
        0.1549614686894405E5 * t524 * t543 * t540 -
        0.1549614686894405E5 * t524 * t50 * t154 * t540;
      double t571 = 0.1371740148509607E1 + t6 + t1;
      double t572 = 0.1091700181433142E1 + t6 + t1;
      double t573 = t572 * t571;
      double t575 = 0.7092992179024789 + t6 + t1;
      double t576 = 0.2907007820975211 + t6 + t1;
      double t578 = -0.917001814331423E-1 + t6 + t1;
      double t579 = -0.3717401485096066 + t6 + t1;
      double t580 = t579 * t578;
      double t581 = t580 * t576 * t575;
      double t592 = t571 * t4 * t20;
      double t599 = t575 * t572;
      double t604 = t158 * t155 * t10;
      double t618 = t361 * t359;
      double t619 = t618 * t408;
      double t626 =
        -0.329747230601548E5 * t157 * t51 * t66 * t184 * t308 +
        0.8744367519213065E3 * t581 * t573 * t20 -
        0.8744367519213065E3 * t581 * t573 * t5 +
        0.1748873503842613E4 * t581 * t572 * t4 * t20 +
        0.1748873503842613E4 * t581 * t592 +
        0.1748873503842613E4 * t580 * t576 * t572 * t592 +
        0.1748873503842613E4 * t580 * t599 * t592 +
        0.3863456908359779E4 * t604 * t473 + 0.193172845417989E4 * t604 * t476 -
        0.193172845417989E4 * t604 * t459 + 0.3863456908359779E4 * t604 * t465 +
        0.3863456908359779E4 * t158 * t155 * t9 * t465 +
        0.2283145657860649E3 * t619 * t357 * t5 +
        0.2283145657860649E3 * t619 * t356 * t4 * t20;
      double t643 = t15 * t67;
      double t670 =
        -0.2283145657860649E3 * t619 * t357 * t20 -
        0.2023596358931899E4 * t56 * t48 * t4 * t20 -
        0.193172845417989E4 * t176 * t11 * t465 -
        0.193172845417989E4 * t158 * t11 * t465 -
        0.193172845417989E4 * t172 * t11 * t465 -
        0.4208935119439625E4 * t643 * t399 + 0.4208935119439625E4 * t643 * t456 +
        0.1748873503842613E4 * t579 * t576 * t599 * t592 +
        0.8417870238879249E4 * t643 * t435 +
        0.1748873503842613E4 * t578 * t576 * t599 * t592 +
        0.8417870238879249E4 * t643 * t308 +
        0.8417870238879249E4 * t15 * t184 * t308 +
        0.8417870238879249E4 * t15 * t497 * t308 -
        0.4208935119439625E4 * t14 * t66 * t184 * t308;
      double t699 = t48 * t127;
      double t712 =
        -0.4208935119439625E4 * t13 * t66 * t184 * t308 +
        0.2283145657860649E3 * t619 * t404 +
        0.2283145657860649E3 * t360 * t371 * t404 +
        0.2283145657860649E3 * t361 * t358 * t371 * t404 +
        0.2283145657860649E3 * t618 * t371 * t404 +
        0.2283145657860649E3 * t618 * t358 * t356 * t404 -
        0.3199586779524209E3 * t83 * t354 - 0.3199586779524209E3 * t76 * t354 -
        0.3199586779524209E3 * t56 * t354 +
        0.3199586779524209E3 * t56 * t699 * t20 -
        0.3199586779524209E3 * t56 * t699 * t5 -
        0.3199586779524209E3 * t118 * t354 - 0.3199586779524209E3 * t122 * t354 -
        0.7678307705366328E3 * t28 * t528;
      double t715 = t22 * t154;
      double t738 = t13 * t4;
      double t750 =
        0.7678307705366328E3 * t28 * t533 +
        0.7678307705366328E3 * t148 * t715 * t540 +
        0.7678307705366328E3 * t27 * t715 * t540 +
        0.7678307705366328E3 * t27 * t23 * t154 * t540 -
        0.1535661541073266E4 * t28 * t540 - 0.1535661541073266E4 * t28 * t523 +
        0.7678307705366328E3 * t144 * t715 * t540 -
        0.6734669543451902E2 * t15 * t5 -
        0.8495746628803322E4 * t154 * t539 * t20 +
        0.1820200428073255E5 * t14 * t738 * t20 +
        0.6734669543451902E2 * t15 * t20 - 0.6734669543451902E2 * t738 * t20 -
        0.6734669543451902E2 * t14 * t4 * t20 + 0.2124194368544568E3 * t354;
      double t755 = t187 * t127;
      double t764 = t187 * t171;
      double t783 =
        -0.4248388737089136E3 * t187 * t4 * t20 -
        0.6717292397503965E2 * t755 * t5 + 0.6717292397503965E2 * t755 * t20 +
        0.7121653599037271E3 * t528 + 0.1424330719807454E4 * t523 +
        0.1424330719807454E4 * t540 - 0.7121653599037271E3 * t533 +
        0.2686590977071556E4 * t764 * t5 - 0.2686590977071556E4 * t764 * t20 -
        0.5373181954143113E4 * t187 * t522 * t20 -
        0.5373181954143113E4 * t187 * t539 * t20 -
        0.3156381276825331E4 * t201 * t5 -
        0.3951518069564725E4 * t157 * t51 * t4 * t20 +
        0.3951518069564725E4 * t215 * t20;
      double t786 = t155 * t4;
      double t805 = t14 * t13 * t127;
      double t822 =
        -0.3951518069564725E4 * t215 * t5 -
        0.3951518069564725E4 * t157 * t786 * t20 -
        0.3951518069564725E4 * t51 * t786 * t20 +
        0.3156381276825331E4 * t201 * t20 +
        0.6312762553650661E4 * t10 * t472 * t20 +
        0.6312762553650661E4 * t10 * t464 * t20 +
        0.6312762553650661E4 * t9 * t464 * t20 +
        0.2877989575362488E4 * t805 * t5 - 0.2877989575362488E4 * t805 * t20 +
        0.2877989575362488E4 * t14 * t353 * t20 +
        0.2877989575362488E4 * t13 * t353 * t20 -
        0.329747230601548E5 * t51 * t309 * t184 * t308 -
        0.6113283677012599E2 * t5 + 0.6113283677012599E2 * t20;
      double t826 =
        t81 + t142 + t204 + t247 + t287 + t331 + t407 + t468 + t512 + t564 +
        t626 + t670 + t712 + t750 + t783 + t822;
      return t826;
    }

    double
    ortho2_f64y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t6 = 0.1E1 * y;
      double t7 = -t1 - t6 + 0.3997579954114602;
      double t8 = -t1 - t6 + 0.1771862795107378;
      double t9 = t8 * t7;
      double t11 = -t1 - t6 - 0.1368825361738218;
      double t12 = -t1 - t6 - 0.5;
      double t14 = -t1 - t6 - 0.8631174638261782;
      double t15 = -t1 - t6 - 0.1177186279510738E1;
      double t17 = -t1 - t6 - 0.139975799541146E1;
      double t18 = t17 * t15 * t14;
      double t19 = t18 * t12 * t11;
      double t22 = 0.1E1 * x;
      double t23 = 0.1371740148509607E1 + t22 + t2;
      double t24 = 0.1091700181433142E1 + t22 + t2;
      double t25 = t24 * t23;
      double t27 = 0.7092992179024789 + t22 + t2;
      double t28 = 0.2907007820975211 + t22 + t2;
      double t30 = -0.917001814331423E-1 + t22 + t2;
      double t31 = -0.3717401485096066 + t22 + t2;
      double t32 = t31 * t30;
      double t33 = t32 * t28 * t27;
      double t36 = 0.5 + t2;
      double t37 = t3 * t36;
      double t38 = 0.1330223896278567E1 + t22 + t2;
      double t40 = t38 * t4 * t37;
      double t41 = 0.9688487934707142 + t22 + t2;
      double t42 = 0.5 + t22 + t2;
      double t43 = t42 * t41;
      double t44 = 0.3115120652928579E-1 + t22 + t2;
      double t45 = -0.3302238962785669 + t22 + t2;
      double t46 = t45 * t44;
      double t47 = -t1 - t6 - 0.9472135954999579;
      double t52 = t4 * t37;
      double t53 = -t1 - t6 + 0.3302238962785669;
      double t54 = -t1 - t6 - 0.3115120652928579E-1;
      double t55 = t54 * t53;
      double t56 = -t1 - t6 - 0.9688487934707142;
      double t57 = t56 * t12;
      double t58 = t57 * t55;
      double t61 = 0.1154653670707977E1 + t22 + t2;
      double t62 = t61 * t5;
      double t63 = -0.1546536707079771 + t22 + t2;
      double t64 = t63 * t42;
      double t65 = -t1 - t6 - 0.5278640450004206E-1;
      double t66 = t47 * t65;
      double t67 = t66 * t64;
      double t70 = 0.1265055323929465E1 + t22 + t2;
      double t72 = t70 * t4 * t37;
      double t73 = 0.7852315164806451 + t22 + t2;
      double t74 = 0.2147684835193549 + t22 + t2;
      double t75 = t74 * t73;
      double t76 = -t1 - t6 + 0.1546536707079771;
      double t77 = t12 * t76;
      double t78 = -t1 - t6 - 0.1154653670707977E1;
      double t79 = t78 * t77;
      double t83 = t4 * t36;
      double t84 = 0.9472135954999579 + t22 + t2;
      double t85 = 0.5278640450004206E-1 + t22 + t2;
      double t86 = t85 * t84;
      double t87 = t86 * t83;
      double t88 = -t1 - t6 - 0.1330223896278567E1;
      double t89 = t88 * t57;
      double t90 = t89 * t55;
      double t93 = t61 * t4;
      double t94 = t93 * t37;
      double t96 = t78 * t12;
      double t97 = t96 * t76 * t63;
      double t104 = t70 * t83;
      double t105 = -0.2650553239294647 + t22 + t2;
      double t108 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t109 = t108 * t105;
      double t110 = t109 * t75;
      double t117 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t118 = -t1 - t6 + 0.2650553239294647;
      double t119 = t118 * t117;
      double t120 = -t1 - t6 - 0.2147684835193549;
      double t121 = -t1 - t6 - 0.7852315164806451;
      double t122 = t121 * t120;
      double t126 = t117 * t4;
      double t127 = t126 * t37;
      double t129 = t88 * t56;
      double t130 = t129 * t12 * t53;
      double t133 =
        0.3113553674158404E3 * t19 * t9 * t5 +
        0.8744367519213065E3 * t33 * t25 * t5 -
        0.2500586384583724E5 * t47 * t46 * t43 * t40 -
        0.5834536221727242E4 * t58 * t52 + 0.2648751671758172E5 * t67 * t62 +
        0.329747230601548E5 * t79 * t75 * t72 - 0.1549614686894405E5 * t90 * t87 +
        0.193172845417989E4 * t97 * t94 +
        0.193172845417989E4 * t96 * t76 * t42 * t94 +
        0.3604828261401262E4 * t104 * t110 - 0.3604828261401262E4 * t110 * t52 +
        0.9168362198145701E4 * t122 * t119 * t52 -
        0.6399173559048418E3 * t130 * t127;
      double t134 = t129 * t55;
      double t137 = t78 * t76;
      double t143 = t12 * t54;
      double t144 = t129 * t143;
      double t151 = -t1 - t6 - 0.1265055323929465E1;
      double t152 = t151 * t122;
      double t155 = t42 * t61;
      double t156 = t63 * t155;
      double t162 = t73 * t4 * t37;
      double t163 = t105 * t74;
      double t164 = t79 * t163;
      double t167 = t74 * t70;
      double t171 = t66 * t163;
      double t178 = -t1 - t6 + 0.3717401485096066;
      double t179 = t178 * t117;
      double t181 = -t1 - t6 + 0.917001814331423E-1;
      double t182 = -t1 - t6 - 0.2907007820975211;
      double t183 = t182 * t181;
      double t184 = -t1 - t6 - 0.7092992179024789;
      double t185 = -t1 - t6 - 0.1091700181433142E1;
      double t186 = t185 * t184;
      double t187 = -t1 - t6 - 0.1371740148509607E1;
      double t188 = t187 * t186;
      double t189 = t188 * t183;
      double t192 = t86 * t5;
      double t194 =
        -0.6399173559048418E3 * t134 * t127 -
        0.5699405234470576E5 * t137 * t86 * t52 -
        0.101179817946595E4 * t144 * t53 * t4 * t37 -
        0.3863456908359779E4 * t96 * t64 * t94 +
        0.4292969711775697E3 * t152 * t118 * t83 -
        0.3156381276825331E4 * t156 * t83 + 0.8585939423551394E3 * t152 * t52 +
        0.329747230601548E5 * t164 * t162 -
        0.3604828261401262E4 * t109 * t167 * t52 +
        0.4208935119439625E4 * t171 * t162 - 0.6968241958371119E4 * t156 * t52 -
        0.5834536221727242E4 * t134 * t52 -
        0.1105752092747433E4 * t189 * t179 * t5 + 0.7121653599037271E3 * t192;
      double t203 = t7 * t4 * t37;
      double t204 = t11 * t8;
      double t210 = t117 * t83;
      double t211 = t120 * t118;
      double t212 = t151 * t121;
      double t213 = t212 * t211;
      double t216 = t105 * t73;
      double t220 = t108 * t64;
      double t224 = t151 * t121 * t118;
      double t227 = t151 * t211;
      double t244 =
        -0.5699405234470576E5 * t77 * t86 * t52 -
        0.3706919888485688E4 * t65 * t86 * t52 -
        0.6227107348316808E3 * t18 * t204 * t203 +
        0.4208935119439625E4 * t171 * t72 + 0.458418109907285E4 * t213 * t210 +
        0.329747230601548E5 * t79 * t216 * t72 -
        0.1101775793780247E4 * t220 * t62 + 0.8585939423551394E3 * t224 * t52 +
        0.8585939423551394E3 * t227 * t52 -
        0.2500586384583724E5 * t65 * t46 * t43 * t40 -
        0.8026106670713796E2 * t79 * t210 -
        0.1605221334142759E3 * t78 * t12 * t117 * t52 -
        0.3863456908359779E4 * t137 * t64 * t94 + 0.6113283677012599E2 * t5;
      double t245 = t73 * t70;
      double t249 = t65 * t4;
      double t253 = t245 * t83;
      double t256 = t84 * t4;
      double t257 = t256 * t37;
      double t260 = t85 * t4;
      double t261 = t260 * t37;
      double t265 = t88 * t12 * t55;
      double t274 = t117 * t5;
      double t277 = t245 * t5;
      double t280 = t181 * t178;
      double t281 = t184 * t182;
      double t286 = t46 * t43;
      double t293 = t178 * t4 * t37;
      double t296 =
        -0.3604828261401262E4 * t109 * t245 * t52 +
        0.9101002140366277E4 * t47 * t249 * t37 -
        0.329747230601548E5 * t164 * t253 + 0.1549614686894405E5 * t90 * t257 +
        0.1549614686894405E5 * t90 * t261 - 0.6399173559048418E3 * t265 * t127 -
        0.6399173559048418E3 * t58 * t127 -
        0.3604828261401262E4 * t108 * t74 * t245 * t52 +
        0.8026106670713796E2 * t79 * t274 + 0.4208935119439625E4 * t171 * t277 +
        0.2211504185494865E4 * t185 * t281 * t280 * t127 -
        0.694593792129331E4 * t286 * t40 +
        0.2917268110863621E4 * t144 * t53 * t5 +
        0.349669514057964E4 * t189 * t293;
      double t300 = t47 * t65 * t45;
      double t317 = t155 * t5;
      double t320 = t42 * t4;
      double t324 = t61 * t83;
      double t331 = t41 * t38;
      double t332 = t331 * t5;
      double t333 = t44 * t42;
      double t334 = t108 * t45;
      double t335 = t334 * t333;
      double t345 =
        0.1250293192291862E5 * t300 * t43 * t40 +
        0.1250293192291862E5 * t47 * t65 * t44 * t43 * t40 -
        0.6227107348316808E3 * t17 * t15 * t12 * t204 * t203 -
        0.3863456908359779E4 * t77 * t64 * t94 +
        0.193172845417989E4 * t97 * t317 +
        0.3156381276825331E4 * t63 * t320 * t37 -
        0.2648751671758172E5 * t67 * t324 + 0.2648751671758172E5 * t67 * t52 +
        0.3526384116063255E3 * t108 * t83 - 0.1098249215871111E4 * t335 * t332 +
        0.1549614686894405E5 * t90 * t192 -
        0.2279893575954566E5 * t163 * t245 * t52 +
        0.7678307705366328E3 * t213 * t87;
      double t354 = t331 * t83;
      double t363 = t300 * t333;
      double t366 = t63 * t61;
      double t371 = t41 * t4 * t37;
      double t378 = t84 * t83;
      double t380 = t96 * t76 * t85;
      double t383 = t121 * t211;
      double t390 =
        -0.7678307705366328E3 * t213 * t261 +
        0.4208935119439625E4 * t66 * t216 * t72 +
        0.4208935119439625E4 * t66 * t75 * t72 +
        0.1098249215871111E4 * t335 * t354 - 0.1182844640594961E4 * t52 -
        0.134693390869038E3 * t249 * t37 -
        0.1101775793780247E4 * t108 * t155 * t52 +
        0.1250293192291862E5 * t363 * t332 +
        0.2648751671758172E5 * t66 * t366 * t52 -
        0.1098249215871111E4 * t335 * t371 -
        0.8417870238879249E4 * t47 * t105 * t75 * t72 -
        0.2849702617235288E5 * t380 * t378 + 0.8585939423551394E3 * t383 * t52 -
        0.790303613912945E4 * t78 * t12 * t4 * t37;
      double t392 = t76 * t117;
      double t399 = t47 * t65 * t85;
      double t402 = t76 * t4;
      double t419 = t53 * t85;
      double t423 = t84 * t5;
      double t436 =
        -0.1605221334142759E3 * t78 * t392 * t52 +
        0.2849702617235288E5 * t380 * t52 - 0.1853459944242844E4 * t399 * t378 -
        0.790303613912945E4 * t78 * t402 * t37 +
        0.3156381276825331E4 * t63 * t93 * t37 +
        0.3156381276825331E4 * t42 * t93 * t37 +
        0.2648751671758172E5 * t66 * t155 * t52 -
        0.3099229373788809E5 * t89 * t54 * t85 * t257 -
        0.3099229373788809E5 * t89 * t419 * t257 +
        0.1853459944242844E4 * t399 * t423 - 0.6734669543451902E2 * t66 * t83 -
        0.1597138905330039E4 * t286 * t38 * t83 +
        0.1597138905330039E4 * t286 * t52 -
        0.8744367519213065E3 * t33 * t25 * t83;
      double t451 = t108 * t86;
      double t454 = t187 * t185;
      double t455 = t454 * t281;
      double t476 =
        -0.1250293192291862E5 * t363 * t354 + 0.331228504972317E3 * t210 -
        0.3526384116063255E3 * t108 * t5 - 0.331228504972317E3 * t274 -
        0.5297503343516344E5 * t47 * t63 * t155 * t52 -
        0.5297503343516344E5 * t65 * t63 * t155 * t52 -
        0.2686590977071556E4 * t451 * t5 + 0.4566291315721297E3 * t455 * t293 +
        0.4566291315721297E3 * t454 * t184 * t181 * t293 -
        0.8417870238879249E4 * t65 * t105 * t75 * t72 -
        0.1098249215871111E4 * t335 * t40 +
        0.4566291315721297E3 * t454 * t183 * t293 -
        0.1101775793780247E4 * t220 * t52 +
        0.1597138905330039E4 * t286 * t38 * t5;
      double t494 = t155 * t83;
      double t496 = t152 * t118 * t63;
      double t503 = t105 * t75;
      double t508 = t320 * t37;
      double t519 =
        0.4566291315721297E3 * t187 * t184 * t183 * t293 +
        0.1597138905330039E4 * t46 * t42 * t38 * t52 +
        0.1597138905330039E4 * t46 * t331 * t52 -
        0.4292969711775697E3 * t152 * t118 * t5 -
        0.3283447430294861E5 * t496 * t494 + 0.6734669543451902E2 * t66 * t5 -
        0.458418109907285E4 * t213 * t274 - 0.170796896734225E4 * t503 * t104 -
        0.3951518069564725E4 * t79 * t83 + 0.3283447430294861E5 * t496 * t508 +
        0.3283447430294861E5 * t496 * t94 + 0.170796896734225E4 * t503 * t52 +
        0.3283447430294861E5 * t152 * t118 * t42 * t94;
      double t531 = t23 * t4 * t37;
      double t547 = t76 * t105;
      double t553 = t47 * t65 * t117;
      double t559 = t44 * t41;
      double t566 =
        0.170796896734225E4 * t105 * t167 * t52 +
        0.170796896734225E4 * t105 * t245 * t52 +
        0.1597138905330039E4 * t45 * t42 * t331 * t52 +
        0.8744367519213065E3 * t33 * t531 + 0.3156381276825331E4 * t156 * t5 +
        0.329747230601548E5 * t164 * t72 - 0.7678307705366328E3 * t213 * t257 +
        0.1535661541073266E4 * t212 * t120 * t85 * t257 -
        0.1101775793780247E4 * t108 * t366 * t52 -
        0.6594944612030961E5 * t12 * t547 * t75 * t72 +
        0.2877989575362488E4 * t553 * t83 +
        0.5755979150724977E4 * t47 * t126 * t37 -
        0.1098249215871111E4 * t334 * t559 * t40 -
        0.1098249215871111E4 * t334 * t43 * t40;
      double t578 = t27 * t24;
      double t583 = t14 * t12;
      double t598 = t70 * t5;
      double t616 =
        -0.3113553674158404E3 * t19 * t9 * t83 -
        0.3099229373788809E5 * t88 * t56 * t54 * t419 * t257 -
        0.4208935119439625E4 * t171 * t253 +
        0.8744367519213065E3 * t31 * t28 * t578 * t531 -
        0.6227107348316808E3 * t17 * t583 * t204 * t203 -
        0.6227107348316808E3 * t15 * t583 * t204 * t203 -
        0.790303613912945E4 * t12 * t402 * t37 +
        0.2283145657860649E3 * t455 * t280 * t83 -
        0.3604828261401262E4 * t110 * t598 +
        0.5755979150724977E4 * t65 * t126 * t37 -
        0.6227107348316808E3 * t19 * t8 * t4 * t37 -
        0.6227107348316808E3 * t19 * t203 -
        0.6594944612030961E5 * t78 * t547 * t75 * t72 +
        0.329747230601548E5 * t164 * t277;
      double t621 = t108 * t117;
      double t648 =
        -0.2124194368544568E3 * t108 * t4 * t37 + 0.4248388737089136E3 * t127 +
        0.6717292397503965E2 * t621 * t5 + 0.7121653599037271E3 * t257 +
        0.7121653599037271E3 * t261 + 0.1853459944242844E4 * t399 * t52 +
        0.1853459944242844E4 * t47 * t65 * t84 * t52 +
        0.1597138905330039E4 * t333 * t331 * t52 +
        0.2686590977071556E4 * t451 * t83 -
        0.1098249215871111E4 * t108 * t44 * t43 * t40 -
        0.6113283677012599E2 * t83 + 0.144964534797642E5 * t213 * t52 +
        0.2849702617235288E5 * t380 * t423 + 0.3951518069564725E4 * t79 * t5;
      double t694 =
        -0.6227107348316808E3 * t18 * t12 * t8 * t203 -
        0.2283145657860649E3 * t455 * t280 * t5 +
        0.1250293192291862E5 * t363 * t40 +
        0.4566291315721297E3 * t186 * t183 * t293 -
        0.1605221334142759E3 * t12 * t392 * t52 -
        0.6566894860589722E5 * t152 * t64 * t94 -
        0.6566894860589722E5 * t224 * t64 * t94 +
        0.1250293192291862E5 * t300 * t559 * t40 -
        0.193172845417989E4 * t97 * t494 +
        0.8744367519213065E3 * t30 * t28 * t578 * t531 +
        0.8744367519213065E3 * t32 * t28 * t24 * t531 +
        0.8744367519213065E3 * t32 * t578 * t531 -
        0.6566894860589722E5 * t227 * t64 * t94 -
        0.6566894860589722E5 * t383 * t64 * t94;
      double t731 = t151 * t120;
      double t739 = t118 * t85;
      double t743 =
        -0.3099229373788809E5 * t88 * t143 * t419 * t257 -
        0.3099229373788809E5 * t56 * t143 * t419 * t257 +
        0.1105752092747433E4 * t189 * t179 * t83 +
        0.2211504185494865E4 * t189 * t127 +
        0.2211504185494865E4 * t188 * t182 * t178 * t127 +
        0.9168362198145701E4 * t212 * t120 * t117 * t52 -
        0.7121653599037271E3 * t87 + 0.2211504185494865E4 * t188 * t280 * t127 +
        0.2211504185494865E4 * t187 * t185 * t182 * t280 * t127 +
        0.170796896734225E4 * t74 * t245 * t52 +
        0.9168362198145701E4 * t212 * t119 * t52 +
        0.9168362198145701E4 * t731 * t119 * t52 +
        0.2211504185494865E4 * t187 * t281 * t280 * t127 +
        0.1535661541073266E4 * t212 * t739 * t257;
      double t788 =
        0.1535661541073266E4 * t731 * t739 * t257 +
        0.1535661541073266E4 * t122 * t739 * t257 -
        0.5834536221727242E4 * t144 * t52 -
        0.1699149325760664E5 * t85 * t256 * t37 +
        0.8744367519213065E3 * t33 * t24 * t4 * t37 +
        0.193172845417989E4 * t97 * t508 -
        0.3706919888485688E4 * t47 * t86 * t52 +
        0.2849702617235288E5 * t96 * t76 * t84 * t52 -
        0.5699405234470576E5 * t96 * t86 * t52 -
        0.6594944612030961E5 * t78 * t12 * t105 * t75 * t72 +
        0.1250293192291862E5 * t363 * t371 - 0.2877989575362488E4 * t553 * t5 -
        0.134693390869038E3 * t47 * t4 * t37 +
        0.4566291315721297E3 * t455 * t181 * t4 * t37;
      double t795 = t53 * t117;
      double t823 =
        -0.5834536221727242E4 * t265 * t52 - 0.2538077782292664E3 * t79 * t52 -
        0.6717292397503965E2 * t621 * t83 +
        0.3199586779524209E3 * t144 * t795 * t5 -
        0.5834536221727242E4 * t130 * t52 -
        0.2917268110863621E4 * t144 * t53 * t83 -
        0.7678307705366328E3 * t213 * t192 - 0.6399173559048418E3 * t144 * t127 -
        0.2686590977071556E4 * t108 * t260 * t37 -
        0.2686590977071556E4 * t108 * t256 * t37 +
        0.3283447430294861E5 * t496 * t317 -
        0.3199586779524209E3 * t144 * t795 * t83 +
        0.1101775793780247E4 * t220 * t324 + 0.170796896734225E4 * t503 * t598;
      double t827 =
        t133 + t194 + t244 + t296 + t345 + t390 + t436 + t476 + t519 + t566 +
        t616 + t648 + t694 + t743 + t788 + t823;
      return t827;
    }

    // * f65 *********************************************************************

    double
    ortho2_f65 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = 0.1E1 * x;
      double t8 = 0.9472135954999579 + t7 + t1;
      double t9 = t8 * t6;
      double t10 = 0.5278640450004206E-1 + t7 + t1;
      double t12 = t10 * t9 * t5;
      double t16 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t17 = t16 * t6;
      double t20 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t24 = 0.1E1 * y;
      double t25 = -t3 - t24 - 0.5278640450004206E-1;
      double t27 = -t3 - t24 - 0.9472135954999579;
      double t31 = t6 * t5;
      double t37 = -t3 - t24 + 0.1546536707079771;
      double t38 = -t3 - t24 - 0.5;
      double t40 = -t3 - t24 - 0.1154653670707977E1;
      double t41 = t40 * t38 * t37;
      double t44 = 0.1154653670707977E1 + t7 + t1;
      double t45 = 0.5 + t7 + t1;
      double t46 = t45 * t44;
      double t47 = -0.1546536707079771 + t7 + t1;
      double t51 = t10 * t8;
      double t55 = -t3 - t24 + 0.2650553239294647;
      double t56 = -t3 - t24 - 0.2147684835193549;
      double t58 = -t3 - t24 - 0.7852315164806451;
      double t59 = -t3 - t24 - 0.1265055323929465E1;
      double t61 = t59 * t58 * t56 * t55;
      double t68 = t27 * t25;
      double t72 = 0.1265055323929465E1 + t7 + t1;
      double t73 = 0.7852315164806451 + t7 + t1;
      double t75 = 0.2147684835193549 + t7 + t1;
      double t76 = -0.2650553239294647 + t7 + t1;
      double t77 = t76 * t75;
      double t82 = t40 * t38;
      double t87 = (0.1330223896278567E1 + t7 + t1) * t6;
      double t88 = 0.9688487934707142 + t7 + t1;
      double t91 = 0.3115120652928579E-1 + t7 + t1;
      double t93 = -0.3302238962785669 + t7 + t1;
      double t99 = t72 * t6;
      double t105 = -t3 - t24 + 0.3302238962785669;
      double t106 = -t3 - t24 - 0.3115120652928579E-1;
      double t108 = -t3 - t24 - 0.9688487934707142;
      double t110 = -t3 - t24 - 0.1330223896278567E1;
      double t112 = t110 * t108 * t38 * t106 * t105;
      double t116 = (0.1371740148509607E1 + t7 + t1) * t6;
      double t117 = 0.1091700181433142E1 + t7 + t1;
      double t120 = 0.7092992179024789 + t7 + t1;
      double t121 = 0.2907007820975211 + t7 + t1;
      double t123 = -0.917001814331423E-1 + t7 + t1;
      double t124 = -0.3717401485096066 + t7 + t1;
      double t130 =
        -0.2988816868777699E4 * t12 - 0.5281227320252506E3 * t20 * t17 * t5 -
        0.2892514555814678E4 * t27 * t25 * t6 * t5 - 0.4275353643066911E3 * t31 -
        0.544183830106299E4 * t27 * t25 * t16 * t31 +
        0.5693275451437552E4 * t41 * t31 +
        0.9763906215069606E4 * t47 * t46 * t31 -
        0.6484667054292267E4 * t20 * t51 * t31 -
        0.2808604154251295E4 * t61 * t31 +
        0.4468541289188007E4 * t20 * t47 * t46 * t31 -
        0.1926653845960513E5 * t68 * t51 * t31 -
        0.3375151029277206E4 * t77 * t73 * t72 * t31 +
        0.3679824669625454E4 * t82 * t37 * t16 * t31 +
        0.5311034339944347E5 * t27 * t25 * t93 * t91 * t45 * t88 * t87 * t5 +
        0.6475518832663362E5 * t41 * t77 * t73 * t99 * t5 +
        0.1235611895559165E5 * t112 * t12 -
        0.5332552025469556E4 * t20 * t124 * t123 * t121 * t120 * t117 * t116 * t5;
      double t131 = -t3 - t24 + 0.3717401485096066;
      double t142 =
        (-t3 - t24 - 0.1371740148509607E1) * (-t3 - t24 -
                0.1091700181433142E1) * (-t3 - t24 -
                       0.7092992179024789)
        * (-t3 - t24 - 0.2907007820975211) * (-t3 - t24 + 0.917001814331423E-1);
      double t160 = t44 * t6;
      double t165 = t59 * t58 * t56;
      double t169 = t87 * t5;
      double t170 = t45 * t88;
      double t171 = t93 * t91;
      double t175 = t9 * t5;
      double t180 = t17 * t5;
      double t183 = t99 * t5;
      double t184 = t75 * t73;
      double t189 = t160 * t5;
      double t190 = t47 * t45;
      double MapleGenVar1 =
        -0.6333884631802828E3 * t142 * t131 * t17 * t5 +
        0.14992522414727E3 * (-t3 - t24 - 0.139975799541146E1) * (-t3 - t24 -
                        0.1177186279510738E1)
        * (-t3 - t24 - 0.8631174638261782) * t38 * (-t3 - t24 -
                0.1368825361738218) * (-t3 -
                           t24 +
                           0.1771862795107378)
        * (-t3 - t24 + 0.3997579954114602) * t6 * t5 +
        0.3909813047090412E5 * t165 * t55 * t47 * t45 * t160 * t5 +
        0.100096466702681E5 * t171 * t170 * t169 +
        0.4401842632932012E5 * t82 * t37 * t10 * t175 -
        0.5001594966430308E4 * t61 * t180 -
        0.1256966354416985E5 * t20 * t76 * t184 * t183 +
        0.6010073026537092E5 * t68 * t190 * t189 +
        0.2324464911146762E4 * t110 * t108 * t38 * t106 * t105 * t6 * t5;
      double t234 =
        MapleGenVar1 - 0.3934913214320119E3 * t142 * t131 * t6 * t5 +
        0.1206349486912768E4 * t112 * t180 -
        0.2009361265752683E5 * t27 * t25 * t76 * t184 * t183 -
        0.2314872646412975E5 * t41 * t190 * t189 +
        0.2112721710477496E4 * t20 * t171 * t170 * t169 -
        0.1345822630669334E5 * t165 * t55 * t10 * t175 +
        0.6706658638066913E3 * t124 * t123 * t121 * t120 * t117 * t116 * t5 -
        0.8309007801355711E3 * t180 - 0.7622982032747365E3 * t20 * t6 * t5;
      double t235 = t130 + t234;

      return t235;
    }

    double
    ortho2_f65x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.1154653670707977E1 + t8 + t1;
      double t10 = 0.5 + t8 + t1;
      double t11 = t10 * t9;
      double t14 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t18 = -0.1546536707079771 + t8 + t1;
      double t19 = t18 * t11;
      double t22 = 0.1E1 * y;
      double t23 = -t3 - t22 + 0.2650553239294647;
      double t24 = -t3 - t22 - 0.2147684835193549;
      double t25 = t24 * t23;
      double t26 = -t3 - t22 - 0.7852315164806451;
      double t27 = t26 * t25;
      double t30 = -t3 - t22 - 0.1265055323929465E1;
      double t31 = t30 * t25;
      double t35 = t30 * t26 * t23;
      double t38 = t26 * t24;
      double t39 = t30 * t38;
      double t45 = t18 * t9;
      double t49 = t18 * t10;
      double t50 = t14 * t49;
      double t53 = t6 * t2;
      double t57 = -t3 - t22 + 0.1546536707079771;
      double t58 = -t3 - t22 - 0.5;
      double t59 = t58 * t57;
      double t60 = -t3 - t22 - 0.1154653670707977E1;
      double t61 = t60 * t59;
      double t66 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t67 = t66 * t5;
      double t70 = t66 * t53;
      double t73 = t9 * t5;
      double t76 =
        0.4468541289188007E4 * t14 * t11 * t7 + 0.7065384146169622E4 * t19 * t7 +
        0.1404302077125648E4 * t27 * t7 + 0.1404302077125648E4 * t31 * t7 +
        0.1404302077125648E4 * t35 * t7 + 0.1404302077125648E4 * t39 * t7 -
        0.1404302077125648E4 * t39 * t23 * t5 +
        0.4468541289188007E4 * t14 * t45 * t7 + 0.4468541289188007E4 * t50 * t7 +
        0.1404302077125648E4 * t39 * t23 * t53 - 0.1163662734609306E5 * t61 * t7 +
        0.1839912334812727E4 * t61 * t67 - 0.1839912334812727E4 * t61 * t70 +
        0.2234270644594004E4 * t50 * t73;
      double t77 = t9 * t53;
      double t84 = t57 * t66;
      double t91 = 0.5278640450004206E-1 + t8 + t1;
      double t92 = -t3 - t22 - 0.5278640450004206E-1;
      double t94 = -t3 - t22 - 0.9472135954999579;
      double t95 = t94 * t92 * t91;
      double t98 = 0.9472135954999579 + t8 + t1;
      double t99 = t98 * t5;
      double t102 = t98 * t53;
      double t105 = t91 * t98;
      double t116 = t94 * t92;
      double t117 = t116 * t49;
      double t120 = 0.1265055323929465E1 + t8 + t1;
      double t121 = t120 * t5;
      double t122 = 0.7852315164806451 + t8 + t1;
      double t123 = 0.2147684835193549 + t8 + t1;
      double t124 = t123 * t122;
      double t125 = -0.2650553239294647 + t8 + t1;
      double t126 = t125 * t124;
      double t131 = t30 * t26;
      double t132 = t131 * t25;
      double t135 =
        -0.2234270644594004E4 * t50 * t77 -
        0.1839912334812727E4 * t60 * t58 * t66 * t7 -
        0.1839912334812727E4 * t60 * t84 * t7 -
        0.1839912334812727E4 * t58 * t84 * t7 - 0.1926653845960513E5 * t95 * t7 -
        0.9633269229802566E4 * t95 * t99 + 0.9633269229802566E4 * t95 * t102 +
        0.9633269229802566E4 * t92 * t105 * t7 +
        0.9633269229802566E4 * t94 * t105 * t7 -
        0.1926653845960513E5 * t94 * t92 * t98 * t7 +
        0.3005036513268546E5 * t117 * t73 - 0.1687575514638603E4 * t126 * t121 -
        0.3005036513268546E5 * t117 * t77 + 0.1581643202755318E5 * t132 * t7;
      double t137 = t120 * t53;
      double t140 = t122 * t120;
      double t147 = t123 * t120;
      double t165 = -t3 - t22 + 0.3302238962785669;
      double t167 = -t3 - t22 - 0.3115120652928579E-1;
      double t168 = t58 * t167;
      double t169 = -t3 - t22 - 0.9688487934707142;
      double t170 = -t3 - t22 - 0.1330223896278567E1;
      double t171 = t170 * t169;
      double t172 = t171 * t168;
      double t175 = t125 * t123;
      double t183 = t167 * t165;
      double t184 = t171 * t183;
      double t188 = t171 * t58 * t165;
      double t191 =
        0.1687575514638603E4 * t126 * t137 -
        0.3375151029277206E4 * t123 * t140 * t7 -
        0.3375151029277206E4 * t125 * t140 * t7 -
        0.3375151029277206E4 * t125 * t147 * t7 -
        0.3375151029277206E4 * t126 * t7 + 0.6010073026537092E5 * t117 * t7 +
        0.6010073026537092E5 * t116 * t45 * t7 +
        0.6010073026537092E5 * t116 * t11 * t7 -
        0.3005036513268546E5 * t94 * t18 * t11 * t7 -
        0.1162232455573381E4 * t172 * t165 * t53 -
        0.198743831107806E5 * t175 * t140 * t7 -
        0.3005036513268546E5 * t92 * t18 * t11 * t7 -
        0.1162232455573381E4 * t184 * t7 - 0.1162232455573381E4 * t188 * t7;
      double t197 = 0.9688487934707142 + t8 + t1;
      double t198 = t10 * t197;
      double t199 = 0.3115120652928579E-1 + t8 + t1;
      double t200 = -0.3302238962785669 + t8 + t1;
      double t201 = t200 * t199;
      double t202 = t201 * t198;
      double t205 = 0.1330223896278567E1 + t8 + t1;
      double t212 = t197 * t205;
      double t224 = t199 * t10;
      double t229 = t170 * t58 * t183;
      double t232 = t169 * t58;
      double t233 = t232 * t183;
      double t244 = t23 * t66;
      double t245 = t30 * t24;
      double t249 =
        -0.1162232455573381E4 * t172 * t7 +
        0.1162232455573381E4 * t172 * t165 * t5 +
        0.100096466702681E5 * t202 * t7 +
        0.5004823335134051E4 * t202 * t205 * t5 -
        0.5004823335134051E4 * t202 * t205 * t53 +
        0.100096466702681E5 * t201 * t212 * t7 +
        0.100096466702681E5 * t201 * t10 * t205 * t7 +
        0.100096466702681E5 * t200 * t10 * t212 * t7 +
        0.100096466702681E5 * t224 * t212 * t7 -
        0.1162232455573381E4 * t229 * t7 - 0.1162232455573381E4 * t233 * t7 +
        0.2500797483215154E4 * t132 * t70 - 0.2500797483215154E4 * t132 * t67 +
        0.2500797483215154E4 * t131 * t24 * t66 * t7 +
        0.2500797483215154E4 * t245 * t244 * t7;
      double t259 = t60 * t58;
      double t260 = t259 * t57 * t91;
      double t274 = t60 * t57;
      double t285 = t14 * t125;
      double t292 = t285 * t124;
      double t297 =
        0.2500797483215154E4 * t131 * t244 * t7 +
        0.2500797483215154E4 * t38 * t244 * t7 +
        0.4401842632932012E5 * t260 * t7 + 0.2200921316466006E5 * t260 * t99 -
        0.2200921316466006E5 * t260 * t102 -
        0.2200921316466006E5 * t259 * t105 * t7 +
        0.4401842632932012E5 * t259 * t57 * t98 * t7 -
        0.2200921316466006E5 * t274 * t105 * t7 -
        0.2200921316466006E5 * t59 * t105 * t7 -
        0.1256966354416985E5 * t14 * t123 * t140 * t7 -
        0.1256966354416985E5 * t285 * t140 * t7 -
        0.1256966354416985E5 * t285 * t147 * t7 -
        0.1256966354416985E5 * t292 * t7 - 0.6284831772084924E4 * t292 * t121;
      double t308 = t120 * t6 * t5;
      double t309 = t57 * t125;
      double t318 = t212 * t53;
      double t320 = t94 * t92 * t200;
      double t321 = t320 * t224;
      double t325 = t205 * t6 * t5;
      double t328 = t212 * t5;
      double t331 = t199 * t197;
      double t336 = t197 * t6 * t5;
      double t347 =
        0.6284831772084924E4 * t292 * t137 + 0.4154503900677855E3 * t70 -
        0.4154503900677855E3 * t67 + 0.1422239685538098E4 * t7 +
        0.3811491016373682E3 * t14 * t53 - 0.3811491016373682E3 * t14 * t5 -
        0.3237759416331681E5 * t60 * t309 * t124 * t308 -
        0.3237759416331681E5 * t58 * t309 * t124 * t308 -
        0.2655517169972174E5 * t321 * t318 + 0.5311034339944347E5 * t321 * t325 +
        0.2655517169972174E5 * t321 * t328 +
        0.5311034339944347E5 * t320 * t331 * t325 +
        0.5311034339944347E5 * t321 * t336 +
        0.5311034339944347E5 * t320 * t198 * t325 +
        0.5311034339944347E5 * t94 * t92 * t199 * t198 * t325;
      double t357 = -t3 - t22 + 0.3717401485096066;
      double t358 = t357 * t66;
      double t360 = -t3 - t22 + 0.917001814331423E-1;
      double t361 = -t3 - t22 - 0.2907007820975211;
      double t362 = t361 * t360;
      double t363 = -t3 - t22 - 0.7092992179024789;
      double t364 = -t3 - t22 - 0.1091700181433142E1;
      double t365 = t364 * t363;
      double t366 = -t3 - t22 - 0.1371740148509607E1;
      double t367 = t366 * t365;
      double t368 = t367 * t362;
      double t374 = t66 * t6;
      double t375 = t374 * t5;
      double t376 = t360 * t357;
      double t391 = t363 * t361;
      double t396 = t140 * t53;
      double t397 = t61 * t175;
      double t401 = t357 * t6 * t5;
      double t408 = -t3 - t22 + 0.3997579954114602;
      double t409 = -t3 - t22 + 0.1771862795107378;
      double t410 = t409 * t408;
      double t412 = -t3 - t22 - 0.1368825361738218;
      double t414 = -t3 - t22 - 0.8631174638261782;
      double t415 = -t3 - t22 - 0.1177186279510738E1;
      double t417 = -t3 - t22 - 0.139975799541146E1;
      double t418 = t417 * t415 * t414;
      double t419 = t418 * t58 * t412;
      double t425 =
        -0.2655517169972174E5 * t94 * t201 * t198 * t325 -
        0.2655517169972174E5 * t92 * t201 * t198 * t325 -
        0.3166942315901414E3 * t368 * t358 * t5 +
        0.3166942315901414E3 * t368 * t358 * t53 +
        0.3166942315901414E3 * t367 * t376 * t375 +
        0.3166942315901414E3 * t367 * t361 * t357 * t375 +
        0.3166942315901414E3 * t368 * t375 +
        0.3166942315901414E3 * t366 * t364 * t361 * t376 * t375 +
        0.3166942315901414E3 * t366 * t391 * t376 * t375 -
        0.3237759416331681E5 * t397 * t396 + 0.200295018732339E4 * t368 * t401 +
        0.3166942315901414E3 * t364 * t391 * t376 * t375 -
        0.7496261207363501E2 * t419 * t410 * t53 +
        0.7496261207363501E2 * t419 * t410 * t5;
      double t427 = t408 * t6 * t5;
      double t438 = t412 * t409;
      double t439 = t414 * t58;
      double t452 = t140 * t5;
      double t456 = t122 * t6 * t5;
      double t463 = 0.1371740148509607E1 + t8 + t1;
      double t465 = t463 * t6 * t5;
      double t466 = 0.1091700181433142E1 + t8 + t1;
      double t467 = 0.7092992179024789 + t8 + t1;
      double t468 = t467 * t466;
      double t469 = 0.2907007820975211 + t8 + t1;
      double t470 = -0.917001814331423E-1 + t8 + t1;
      double t471 = t470 * t469;
      double t472 = -0.3717401485096066 + t8 + t1;
      double t477 = t11 * t5;
      double t479 = t39 * t23 * t18;
      double t482 = t11 * t53;
      double t485 = t9 * t6;
      double t486 = t485 * t5;
      double t489 = t10 * t6;
      double t490 = t489 * t5;
      double t496 =
        -0.7496261207363501E2 * t418 * t58 * t409 * t427 -
        0.7496261207363501E2 * t419 * t427 -
        0.7496261207363501E2 * t419 * t409 * t6 * t5 -
        0.7496261207363501E2 * t417 * t439 * t438 * t427 -
        0.7496261207363501E2 * t417 * t415 * t58 * t438 * t427 -
        0.7496261207363501E2 * t418 * t438 * t427 +
        0.3237759416331681E5 * t397 * t452 + 0.6475518832663362E5 * t397 * t456 -
        0.7496261207363501E2 * t415 * t439 * t438 * t427 -
        0.843150507091401E4 * t472 * t471 * t468 * t465 +
        0.1954906523545206E5 * t479 * t477 - 0.1954906523545206E5 * t479 * t482 +
        0.3909813047090412E5 * t479 * t486 + 0.3909813047090412E5 * t479 * t490 -
        0.1954906523545206E5 * t39 * t49 * t486;
      double t504 = t14 * t200;
      double t505 = t504 * t224;
      double t519 = t125 * t122;
      double t541 =
        0.3909813047090412E5 * t39 * t23 * t10 * t486 -
        0.1056360855238748E4 * t505 * t318 -
        0.1954906523545206E5 * t27 * t49 * t486 -
        0.1954906523545206E5 * t31 * t49 * t486 -
        0.1954906523545206E5 * t35 * t49 * t486 +
        0.1056360855238748E4 * t505 * t328 +
        0.6475518832663362E5 * t61 * t519 * t308 +
        0.6475518832663362E5 * t397 * t308 + 0.2112721710477496E4 * t505 * t336 +
        0.2112721710477496E4 * t505 * t325 +
        0.2112721710477496E4 * t14 * t199 * t198 * t325 +
        0.2112721710477496E4 * t504 * t198 * t325 +
        0.2112721710477496E4 * t504 * t331 * t325 +
        0.3340506333597856E4 * t202 * t325;
      double t542 = t98 * t6;
      double t543 = t542 * t5;
      double t544 = t170 * t232;
      double t545 = t544 * t183;
      double t548 = t91 * t6;
      double t549 = t548 * t5;
      double t552 = t105 * t5;
      double t555 = t105 * t53;
      double t561 = t165 * t91;
      double t570 = t466 * t6 * t5;
      double t571 = t469 * t467;
      double t572 = t472 * t470;
      double t573 = t14 * t572;
      double t574 = t573 * t571;
      double t577 = t466 * t463;
      double t578 = t577 * t5;
      double t581 = t577 * t53;
      double t599 =
        0.1235611895559165E5 * t545 * t543 + 0.1235611895559165E5 * t545 * t549 +
        0.6178059477795823E4 * t545 * t552 - 0.6178059477795823E4 * t545 * t555 +
        0.6475518832663362E5 * t61 * t124 * t308 -
        0.6178059477795823E4 * t544 * t561 * t543 -
        0.6178059477795823E4 * t544 * t167 * t91 * t543 -
        0.5332552025469556E4 * t574 * t570 - 0.2666276012734778E4 * t574 * t578 +
        0.2666276012734778E4 * t574 * t581 -
        0.6178059477795823E4 * t169 * t168 * t561 * t543 -
        0.6178059477795823E4 * t170 * t168 * t561 * t543 -
        0.6178059477795823E4 * t170 * t169 * t167 * t561 * t543 -
        0.5332552025469556E4 * t574 * t465;
      double t601 = t469 * t466;
      double t610 = t472 * t469;
      double t625 = t572 * t571;
      double t632 = t92 * t6;
      double t645 = t259 * t57 * t18;
      double t648 =
        -0.5332552025469556E4 * t573 * t601 * t465 -
        0.3237759416331681E5 * t60 * t58 * t125 * t124 * t308 -
        0.5332552025469556E4 * t14 * t610 * t468 * t465 -
        0.5332552025469556E4 * t573 * t468 * t465 -
        0.5332552025469556E4 * t14 * t471 * t468 * t465 -
        0.1025315887970916E5 * t91 * t542 * t5 +
        0.6706658638066913E3 * t625 * t570 + 0.3353329319033456E3 * t625 * t578 -
        0.3353329319033456E3 * t625 * t581 +
        0.1720860368970014E5 * t94 * t632 * t5 +
        0.6706658638066913E3 * t572 * t601 * t465 +
        0.6706658638066913E3 * t625 * t465 +
        0.6706658638066913E3 * t572 * t468 * t465 -
        0.1157436323206488E5 * t645 * t477;
      double t660 = t366 * t364;
      double t661 = t660 * t391;
      double t684 = t116 * t175;
      double t693 =
        0.1157436323206488E5 * t645 * t482 - 0.2314872646412975E5 * t645 * t490 -
        0.2314872646412975E5 * t259 * t57 * t10 * t486 -
        0.2314872646412975E5 * t645 * t486 +
        0.1967456607160059E3 * t661 * t376 * t53 -
        0.1967456607160059E3 * t661 * t376 * t5 -
        0.3814812032819833E4 * t172 * t165 * t6 * t5 +
        0.1967456607160059E3 * t661 * t360 * t6 * t5 +
        0.1157436323206488E5 * t259 * t49 * t486 +
        0.1157436323206488E5 * t274 * t49 * t486 +
        0.1157436323206488E5 * t59 * t49 * t486 -
        0.1004680632876341E5 * t684 * t452 + 0.1004680632876341E5 * t684 * t396 -
        0.2009361265752683E5 * t684 * t308 - 0.2009361265752683E5 * t684 * t456;
      double t716 = t165 * t66;
      double t741 =
        0.6706658638066913E3 * t471 * t468 * t465 +
        0.6706658638066913E3 * t610 * t468 * t465 +
        0.1004680632876341E5 * t94 * t125 * t124 * t308 -
        0.2009361265752683E5 * t116 * t124 * t308 -
        0.2009361265752683E5 * t116 * t519 * t308 +
        0.1004680632876341E5 * t92 * t125 * t124 * t308 -
        0.603174743456384E3 * t172 * t716 * t53 +
        0.1967456607160059E3 * t365 * t362 * t401 +
        0.1967456607160059E3 * t366 * t363 * t362 * t401 +
        0.1967456607160059E3 * t660 * t362 * t401 +
        0.1967456607160059E3 * t660 * t363 * t360 * t401 +
        0.1967456607160059E3 * t661 * t401 +
        0.603174743456384E3 * t172 * t716 * t5 -
        0.603174743456384E3 * t172 * t375;
      double t758 = t23 * t91;
      double t772 = t14 * t105;
      double t780 =
        -0.603174743456384E3 * t184 * t375 - 0.603174743456384E3 * t188 * t375 -
        0.603174743456384E3 * t233 * t375 - 0.603174743456384E3 * t229 * t375 -
        0.6729113153346672E4 * t132 * t552 + 0.6729113153346672E4 * t132 * t555 -
        0.1345822630669334E5 * t132 * t549 - 0.1345822630669334E5 * t132 * t543 +
        0.6729113153346672E4 * t131 * t758 * t543 +
        0.6729113153346672E4 * t131 * t24 * t91 * t543 +
        0.6729113153346672E4 * t38 * t758 * t543 +
        0.6729113153346672E4 * t245 * t758 * t543 +
        0.3242333527146133E4 * t772 * t53 - 0.3242333527146133E4 * t772 * t5 -
        0.6484667054292267E4 * t14 * t548 * t5;
      double t796 = t57 * t6;
      double t805 = t14 * t66;
      double t817 =
        -0.6484667054292267E4 * t14 * t542 * t5 -
        0.4881953107534803E4 * t19 * t53 -
        0.2846637725718776E4 * t60 * t58 * t6 * t5 +
        0.1446257277907339E4 * t94 * t6 * t5 - 0.2846637725718776E4 * t61 * t53 -
        0.2846637725718776E4 * t60 * t796 * t5 + 0.2846637725718776E4 * t61 * t5 -
        0.2846637725718776E4 * t58 * t796 * t5 -
        0.2640613660126253E3 * t805 * t5 + 0.2640613660126253E3 * t805 * t53 -
        0.8350353586552708E3 * t375 + 0.1670070717310542E4 * t14 * t6 * t5 +
        0.1446257277907339E4 * t632 * t5 - 0.2988816868777699E4 * t543;
      double t833 = t94 * t92 * t66;
      double t850 = 0.149440843438885E4 * t555 - 0.149440843438885E4 * t552
        - 0.2988816868777699E4 * t549 + 0.4881953107534803E4 * t19 * t5 +
        0.9763906215069606E4 * t18 * t489 * t5 +
        0.9763906215069606E4 * t10 * t485 * t5 +
        0.9763906215069606E4 * t18 * t485 * t5 +
        0.2720919150531495E4 * t833 * t53 - 0.2720919150531495E4 * t833 * t5 +
        0.2720919150531495E4 * t94 * t374 * t5 +
        0.2720919150531495E4 * t92 * t374 * t5 +
        0.1446257277907339E4 * t116 * t53 - 0.1446257277907339E4 * t116 * t5 -
        0.2137676821533456E3 * t5 + 0.2137676821533456E3 * t53;
      double t854 =
        t76 + t135 + t191 + t249 + t297 + t347 + t425 + t496 + t541 + t599 +
        t648 + t693 + t741 + t780 + t817 + t850;
      return t854;
    }

    double
    ortho2_f65y (double x, double y)
    {
      double t1 = 0.5 * x;
      double t2 = 0.5 * y;
      double t3 = -t1 - t2;
      double t4 = 0.5 + t1;
      double t5 = t4 * t3;
      double t6 = 0.1E1 * x;
      double t7 = 0.9472135954999579 + t6 + t2;
      double t8 = 0.5278640450004206E-1 + t6 + t2;
      double t9 = t8 * t7;
      double t10 = t9 * t5;
      double t11 = 0.1E1 * y;
      double t12 = -t1 - t11 + 0.3302238962785669;
      double t13 = -t1 - t11 - 0.3115120652928579E-1;
      double t14 = t13 * t12;
      double t15 = -t1 - t11 - 0.5;
      double t16 = -t1 - t11 - 0.9688487934707142;
      double t17 = t16 * t15;
      double t18 = -t1 - t11 - 0.1330223896278567E1;
      double t19 = t18 * t17;
      double t20 = t19 * t14;
      double t23 = 0.1154653670707977E1 + t6 + t2;
      double t24 = t23 * t5;
      double t25 = 0.5 + t6 + t2;
      double t26 = -0.1546536707079771 + t6 + t2;
      double t27 = t26 * t25;
      double t28 = -t1 - t11 - 0.5278640450004206E-1;
      double t29 = -t1 - t11 - 0.9472135954999579;
      double t30 = t29 * t28;
      double t31 = t30 * t27;
      double t35 = 0.5 + t2;
      double t36 = t4 * t35;
      double t38 = t3 * t35;
      double t39 = -t1 - t11 + 0.3717401485096066;
      double t41 = t39 * t4 * t38;
      double t42 = -t1 - t11 + 0.917001814331423E-1;
      double t43 = -t1 - t11 - 0.2907007820975211;
      double t44 = t43 * t42;
      double t45 = -t1 - t11 - 0.7092992179024789;
      double t46 = -t1 - t11 - 0.1091700181433142E1;
      double t47 = t46 * t45;
      double t48 = -t1 - t11 - 0.1371740148509607E1;
      double t49 = t48 * t47;
      double t50 = t49 * t44;
      double t53 = 0.1265055323929465E1 + t6 + t2;
      double t54 = 0.7852315164806451 + t6 + t2;
      double t55 = t54 * t53;
      double t56 = t55 * t36;
      double t57 = 0.2147684835193549 + t6 + t2;
      double t58 = -0.2650553239294647 + t6 + t2;
      double t59 = t58 * t57;
      double t60 = -t1 - t11 + 0.1546536707079771;
      double t61 = t15 * t60;
      double t62 = -t1 - t11 - 0.1154653670707977E1;
      double t63 = t62 * t61;
      double t64 = t63 * t59;
      double t67 = t53 * t5;
      double t68 = t57 * t54;
      double t71 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t72 = t71 * t58;
      double t73 = t72 * t68;
      double t76 = -t1 - t11 + 0.3997579954114602;
      double t78 = t76 * t4 * t38;
      double t79 = -t1 - t11 - 0.1368825361738218;
      double t81 = -t1 - t11 - 0.8631174638261782;
      double t82 = -t1 - t11 - 0.1177186279510738E1;
      double t84 = -t1 - t11 - 0.139975799541146E1;
      double t85 = t84 * t82 * t81;
      double t86 = t85 * t15 * t79;
      double t91 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t92 = t91 * t5;
      double t94 = -t1 - t11 + 0.1771862795107378;
      double t99 = 0.1330223896278567E1 + t6 + t2;
      double t100 = 0.9688487934707142 + t6 + t2;
      double t101 = t100 * t99;
      double t102 = t101 * t5;
      double t103 = 0.3115120652928579E-1 + t6 + t2;
      double t104 = t103 * t25;
      double t105 = -0.3302238962785669 + t6 + t2;
      double t106 = t71 * t105;
      double t107 = t106 * t104;
      double t110 = t94 * t76;
      double t116 = -t1 - t11 + 0.2650553239294647;
      double t117 = -t1 - t11 - 0.2147684835193549;
      double t118 = t117 * t116;
      double t119 = -t1 - t11 - 0.7852315164806451;
      double t120 = -t1 - t11 - 0.1265055323929465E1;
      double t121 = t120 * t119;
      double t122 = t121 * t118;
      double t125 =
        0.6178059477795823E4 * t20 * t10 + 0.3005036513268546E5 * t31 * t24 -
        0.2137676821533456E3 * t5 + 0.2137676821533456E3 * t36 +
        0.1001475093661695E4 * t50 * t41 - 0.3237759416331681E5 * t64 * t56 -
        0.6284831772084924E4 * t73 * t67 - 0.14992522414727E3 * t86 * t78 -
        0.4154503900677855E3 * t92 - 0.14992522414727E3 * t86 * t94 * t4 * t38 +
        0.1056360855238748E4 * t107 * t102 -
        0.7496261207363501E2 * t86 * t110 * t36 +
        0.3811491016373682E3 * t71 * t36 - 0.2500797483215154E4 * t122 * t92;
      double t126 = t91 * t36;
      double t131 = t15 * t13;
      double t132 = t18 * t16;
      double t133 = t132 * t131;
      double t136 = t7 * t4;
      double t144 = t4 * t38;
      double t146 = t101 * t36;
      double t148 = t29 * t28 * t105;
      double t149 = t148 * t104;
      double t152 = 0.1371740148509607E1 + t6 + t2;
      double t154 = t152 * t4 * t38;
      double t155 = 0.1091700181433142E1 + t6 + t2;
      double t156 = 0.7092992179024789 + t6 + t2;
      double t157 = t156 * t155;
      double t158 = 0.2907007820975211 + t6 + t2;
      double t159 = -0.917001814331423E-1 + t6 + t2;
      double t160 = t159 * t158;
      double t161 = -0.3717401485096066 + t6 + t2;
      double t167 = t54 * t4 * t38;
      double t170 = t53 * t36;
      double t171 = t58 * t68;
      double t174 = t79 * t94;
      double t175 = t81 * t15;
      double t192 =
        0.4154503900677855E3 * t126 - 0.3811491016373682E3 * t71 * t5 +
        0.1162232455573381E4 * t133 * t12 * t5 -
        0.2050631775941833E5 * t8 * t136 * t38 -
        0.14992522414727E3 * t85 * t15 * t94 * t78 - 0.1096829091182596E4 * t144 -
        0.2655517169972174E5 * t149 * t146 -
        0.1686301014182802E5 * t161 * t160 * t157 * t154 +
        0.3237759416331681E5 * t64 * t167 + 0.1687575514638603E4 * t171 * t170 -
        0.14992522414727E3 * t82 * t175 * t174 * t78 -
        0.14992522414727E3 * t84 * t175 * t174 * t78 -
        0.14992522414727E3 * t84 * t82 * t15 * t174 * t78 -
        0.14992522414727E3 * t85 * t174 * t78;
      double t194 = t42 * t39;
      double t196 = t45 * t43;
      double t197 = t48 * t46;
      double t198 = t197 * t196;
      double t201 = t7 * t5;
      double t203 = t29 * t28 * t8;
      double t208 = t23 * t36;
      double t209 = t71 * t27;
      double t214 = t28 * t4;
      double t221 = t23 * t4;
      double t222 = t221 * t38;
      double t223 = t120 * t118;
      double t227 = t155 * t152;
      double t228 = t227 * t36;
      double t229 = t158 * t156;
      double t230 = t161 * t159;
      double t231 = t71 * t230;
      double t232 = t231 * t229;
      double t238 = t53 * t4 * t38;
      double t243 = t58 * t54;
      double t247 = t91 * t4;
      double t248 = t247 * t38;
      double t253 =
        -0.1967456607160059E3 * t198 * t194 * t5 -
        0.9633269229802566E4 * t203 * t201 - 0.1687575514638603E4 * t171 * t144 -
        0.2234270644594004E4 * t209 * t208 + 0.7908216013776589E4 * t122 * t144 +
        0.860430184485007E4 * t29 * t214 * t38 +
        0.1967456607160059E3 * t198 * t194 * t36 -
        0.3909813047090412E5 * t223 * t27 * t222 +
        0.2666276012734778E4 * t232 * t228 - 0.1687575514638603E4 * t171 * t67 +
        0.2009361265752683E5 * t28 * t58 * t68 * t238 +
        0.3237759416331681E5 * t63 * t243 * t238 - 0.1670070717310542E4 * t248 +
        0.8350353586552708E3 * t71 * t4 * t38;
      double t258 = t57 * t53;
      double t262 = t25 * t23;
      double t263 = t26 * t262;
      double t276 = t136 * t38;
      double t281 = t8 * t4;
      double t282 = t281 * t38;
      double t284 = t71 * t91;
      double t290 = t155 * t4 * t38;
      double t298 =
        -0.1446257277907339E4 * t30 * t5 + 0.2892514555814678E4 * t214 * t38 -
        0.1687575514638603E4 * t58 * t258 * t144 -
        0.4881953107534803E4 * t263 * t36 +
        0.2892514555814678E4 * t29 * t4 * t38 + 0.1446257277907339E4 * t30 * t36 -
        0.149440843438885E4 * t10 + 0.3934913214320119E3 * t198 * t42 * t4 * t38 -
        0.1235611895559165E5 * t19 * t13 * t8 * t276 -
        0.149440843438885E4 * t282 - 0.2640613660126253E3 * t284 * t5 +
        0.2640613660126253E3 * t284 * t36 - 0.2666276012734778E4 * t232 * t290 -
        0.1687575514638603E4 * t58 * t55 * t144 +
        0.2234270644594004E4 * t209 * t24;
      double t301 = t71 * t9;
      double t304 = t119 * t118;
      double t323 = t62 * t15;
      double t327 = t62 * t60;
      double t331 = t26 * t23;
      double t338 = t12 * t91;
      double t344 =
        0.3242333527146133E4 * t301 * t36 -
        0.3909813047090412E5 * t304 * t27 * t222 -
        0.1687575514638603E4 * t57 * t55 * t144 -
        0.1907406016409916E4 * t133 * t12 * t4 * t38 -
        0.3005036513268546E5 * t31 * t208 + 0.2234270644594004E4 * t209 * t144 -
        0.581831367304653E4 * t63 * t144 - 0.1839912334812727E4 * t63 * t126 +
        0.2314872646412975E5 * t323 * t27 * t222 +
        0.2314872646412975E5 * t327 * t27 * t222 +
        0.2234270644594004E4 * t71 * t331 * t144 +
        0.7496261207363501E2 * t86 * t110 * t5 +
        0.603174743456384E3 * t133 * t338 * t5 -
        0.1056360855238748E4 * t107 * t146;
      double t356 = t119 * t117;
      double t357 = t120 * t356;
      double t358 = t357 * t116 * t26;
      double t363 = t227 * t5;
      double t364 = t230 * t229;
      double t371 = t100 * t4 * t38;
      double t375 = t99 * t4 * t38;
      double t376 = t103 * t100;
      double t383 = t262 * t5;
      double t385 = t323 * t60 * t26;
      double t390 = t12 * t8;
      double t394 = t60 * t91;
      double t401 =
        -0.3679824669625454E4 * t62 * t15 * t91 * t144 +
        0.3005036513268546E5 * t31 * t144 +
        0.5001594966430308E4 * t121 * t117 * t91 * t144 +
        0.1954906523545206E5 * t358 * t222 - 0.2666276012734778E4 * t232 * t154 +
        0.3353329319033456E3 * t364 * t363 +
        0.2314872646412975E5 * t61 * t27 * t222 +
        0.1056360855238748E4 * t107 * t371 +
        0.1056360855238748E4 * t106 * t376 * t375 +
        0.3005036513268546E5 * t30 * t331 * t144 -
        0.1157436323206488E5 * t385 * t383 + 0.1056360855238748E4 * t107 * t375 -
        0.1235611895559165E5 * t19 * t390 * t276 -
        0.3679824669625454E4 * t62 * t394 * t144 +
        0.3005036513268546E5 * t30 * t262 * t144;
      double t403 = t25 * t100;
      double t413 = t116 * t91;
      double t424 = t25 * t4;
      double t445 =
        0.1056360855238748E4 * t106 * t403 * t375 +
        0.3237759416331681E5 * t64 * t238 + 0.6284831772084924E4 * t73 * t170 -
        0.6729113153346672E4 * t122 * t10 +
        0.5001594966430308E4 * t121 * t413 * t144 +
        0.3934913214320119E3 * t198 * t41 +
        0.2655517169972174E5 * t148 * t376 * t375 +
        0.2655517169972174E5 * t149 * t375 +
        0.4881953107534803E4 * t26 * t424 * t38 -
        0.6284831772084924E4 * t73 * t144 +
        0.1954906523545206E5 * t357 * t116 * t25 * t222 +
        0.3934913214320119E3 * t197 * t45 * t42 * t41 -
        0.6010073026537092E5 * t29 * t26 * t262 * t144 +
        0.3934913214320119E3 * t197 * t44 * t41;
      double t457 = t39 * t91;
      double t475 = t120 * t117;
      double t495 =
        -0.4401842632932012E5 * t61 * t9 * t144 +
        0.4881953107534803E4 * t26 * t221 * t38 -
        0.6475518832663362E5 * t62 * t15 * t58 * t68 * t238 -
        0.3166942315901414E3 * t50 * t457 * t5 -
        0.6010073026537092E5 * t28 * t26 * t262 * t144 +
        0.1056360855238748E4 * t71 * t103 * t403 * t375 +
        0.2655517169972174E5 * t149 * t371 - 0.149440843438885E4 * t276 -
        0.6284831772084924E4 * t72 * t258 * t144 +
        0.5001594966430308E4 * t475 * t413 * t144 +
        0.3934913214320119E3 * t48 * t45 * t44 * t41 -
        0.397487662215612E5 * t59 * t55 * t144 +
        0.3934913214320119E3 * t47 * t44 * t41 +
        0.2655517169972174E5 * t148 * t403 * t375 +
        0.4881953107534803E4 * t25 * t221 * t38;
      double t502 = t158 * t155;
      double t511 = t262 * t36;
      double t514 = t105 * t103;
      double t536 = t29 * t28 * t91;
      double t545 =
        -0.3679824669625454E4 * t15 * t394 * t144 -
        0.2666276012734778E4 * t231 * t502 * t154 +
        0.2808604154251295E4 * t357 * t144 -
        0.3242333527146133E4 * t71 * t281 * t38 -
        0.1954906523545206E5 * t358 * t511 -
        0.5311034339944347E5 * t28 * t514 * t403 * t375 -
        0.5311034339944347E5 * t29 * t514 * t403 * t375 +
        0.2655517169972174E5 * t29 * t28 * t103 * t403 * t375 -
        0.6284831772084924E4 * t72 * t55 * t144 +
        0.1954906523545206E5 * t358 * t383 - 0.3242333527146133E4 * t301 * t5 -
        0.2720919150531495E4 * t536 * t5 -
        0.603174743456384E3 * t133 * t338 * t36 +
        0.5001594966430308E4 * t356 * t413 * t144;
      double t549 = t120 * t119 * t116;
      double t556 = t132 * t15 * t12;
      double t559 = t30 * t59;
      double t576 = t132 * t14;
      double t588 =
        -0.1206349486912768E4 * t133 * t248 + 0.2808604154251295E4 * t549 * t144 +
        0.3237759416331681E5 * t63 * t68 * t238 -
        0.1206349486912768E4 * t556 * t248 + 0.1004680632876341E5 * t559 * t56 -
        0.1235611895559165E5 * t18 * t16 * t13 * t390 * t276 +
        0.2808604154251295E4 * t223 * t144 -
        0.6284831772084924E4 * t71 * t57 * t55 * t144 -
        0.1162232455573381E4 * t133 * t12 * t36 -
        0.1206349486912768E4 * t576 * t248 - 0.2666276012734778E4 * t232 * t363 -
        0.2324464911146762E4 * t133 * t144 +
        0.1926653845960513E5 * t28 * t9 * t144 -
        0.2324464911146762E4 * t556 * t144;
      double t592 = t161 * t158;
      double t601 = t9 * t36;
      double t623 = t55 * t5;
      double t626 =
        0.2655517169972174E5 * t149 * t102 +
        0.3353329319033456E3 * t592 * t157 * t154 -
        0.3242333527146133E4 * t71 * t136 * t38 +
        0.1839912334812727E4 * t63 * t92 + 0.6729113153346672E4 * t122 * t601 +
        0.2808604154251295E4 * t304 * t144 + 0.149440843438885E4 * t601 -
        0.3353329319033456E3 * t364 * t228 + 0.2720919150531495E4 * t536 * t36 +
        0.3353329319033456E3 * t364 * t290 -
        0.5693275451437552E4 * t62 * t15 * t4 * t38 +
        0.3353329319033456E3 * t230 * t157 * t154 +
        0.3353329319033456E3 * t230 * t502 * t154 -
        0.1004680632876341E5 * t559 * t623;
      double t634 = t7 * t36;
      double t636 = t323 * t60 * t8;
      double t648 = t60 * t4;
      double t662 = t18 * t15 * t14;
      double t669 =
        0.3353329319033456E3 * t364 * t154 - 0.6178059477795823E4 * t20 * t601 -
        0.2666276012734778E4 * t231 * t157 * t154 -
        0.2200921316466006E5 * t636 * t634 - 0.2324464911146762E4 * t576 * t144 +
        0.3353329319033456E3 * t160 * t157 * t154 +
        0.2200921316466006E5 * t636 * t144 - 0.6729113153346672E4 * t122 * t282 -
        0.5693275451437552E4 * t62 * t648 * t38 +
        0.9633269229802566E4 * t203 * t634 -
        0.3909813047090412E5 * t357 * t27 * t222 -
        0.1235611895559165E5 * t18 * t131 * t390 * t276 -
        0.1206349486912768E4 * t662 * t248 - 0.9633269229802566E4 * t203 * t144 -
        0.1004680632876341E5 * t559 * t167;
      double t686 = t514 * t403;
      double t710 =
        -0.2846637725718776E4 * t63 * t36 - 0.6729113153346672E4 * t122 * t276 +
        0.544183830106299E4 * t29 * t247 * t38 -
        0.9633269229802566E4 * t29 * t28 * t7 * t144 +
        0.6178059477795823E4 * t20 * t282 -
        0.5004823335134051E4 * t686 * t99 * t36 +
        0.5004823335134051E4 * t686 * t99 * t5 +
        0.1926653845960513E5 * t29 * t9 * t144 -
        0.1004680632876341E5 * t559 * t238 + 0.1157436323206488E5 * t385 * t511 +
        0.5004823335134051E4 * t686 * t144 -
        0.2666276012734778E4 * t71 * t592 * t157 * t154 +
        0.3166942315901414E3 * t50 * t457 * t36 +
        0.2500797483215154E4 * t122 * t126;
      double t717 = t424 * t38;
      double t731 = t116 * t8;
      double t755 = t17 * t14;
      double t761 =
        0.1345822630669334E5 * t121 * t117 * t8 * t276 -
        0.1157436323206488E5 * t385 * t222 - 0.1157436323206488E5 * t385 * t717 -
        0.3909813047090412E5 * t549 * t27 * t222 -
        0.1157436323206488E5 * t323 * t60 * t25 * t222 +
        0.5004823335134051E4 * t514 * t25 * t99 * t144 +
        0.1345822630669334E5 * t121 * t731 * t276 +
        0.2200921316466006E5 * t323 * t60 * t7 * t144 -
        0.1004680632876341E5 * t30 * t243 * t238 +
        0.6333884631802828E3 * t50 * t248 +
        0.5004823335134051E4 * t514 * t101 * t144 +
        0.5004823335134051E4 * t105 * t25 * t101 * t144 +
        0.2009361265752683E5 * t29 * t58 * t68 * t238 -
        0.1206349486912768E4 * t755 * t248 -
        0.4401842632932012E5 * t323 * t9 * t144;
      double t806 =
        0.1345822630669334E5 * t475 * t731 * t276 +
        0.5004823335134051E4 * t104 * t101 * t144 +
        0.6333884631802828E3 * t49 * t43 * t39 * t248 +
        0.1954906523545206E5 * t358 * t717 - 0.2324464911146762E4 * t662 * t144 -
        0.1404302077125648E4 * t357 * t116 * t5 +
        0.6333884631802828E3 * t49 * t194 * t248 +
        0.6333884631802828E3 * t48 * t46 * t43 * t194 * t248 -
        0.2324464911146762E4 * t755 * t144 +
        0.1345822630669334E5 * t356 * t731 * t276 -
        0.1235611895559165E5 * t16 * t131 * t390 * t276 +
        0.1413076829233924E5 * t263 * t144 -
        0.5693275451437552E4 * t15 * t648 * t38 +
        0.6333884631802828E3 * t48 * t196 * t194 * t248;
      double t822 = t60 * t58;
      double t851 =
        0.6178059477795823E4 * t20 * t276 -
        0.1004680632876341E5 * t30 * t68 * t238 +
        0.6333884631802828E3 * t46 * t196 * t194 * t248 +
        0.2234270644594004E4 * t71 * t262 * t144 +
        0.1404302077125648E4 * t357 * t116 * t36 -
        0.6475518832663362E5 * t15 * t822 * t68 * t238 -
        0.4401842632932012E5 * t327 * t9 * t144 -
        0.6475518832663362E5 * t62 * t822 * t68 * t238 +
        0.3237759416331681E5 * t64 * t623 + 0.2200921316466006E5 * t636 * t201 -
        0.2666276012734778E4 * t71 * t160 * t157 * t154 +
        0.2846637725718776E4 * t63 * t5 + 0.544183830106299E4 * t28 * t247 * t38 +
        0.6681012667195712E4 * t686 * t375 + 0.4881953107534803E4 * t263 * t5;
      double t855 =
        t125 + t192 + t253 + t298 + t344 + t401 + t445 + t495 + t545 + t588 +
        t626 + t669 + t710 + t761 + t806 + t851;
      return t855;
    }

    // * f66 *********************************************************************

    double
    ortho2_f66 (double x, double y)
    {
      double t1 = 0.5 * y;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * (0.5 + t1);
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t9 = 0.1E1 * y;
      double t10 = -t3 - t9 - 0.5278640450004206E-1;
      double t12 = -t3 - t9 - 0.9472135954999579;
      double t16 = 0.1E1 * x;
      double t17 = 0.9472135954999579 + t16 + t1;
      double t18 = t17 * t6;
      double t19 = 0.5278640450004206E-1 + t16 + t1;
      double t21 = t19 * t18 * t5;
      double t25 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t26 = t25 * t6;
      double t29 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t37 = t19 * t17;
      double t41 = 0.1154653670707977E1 + t16 + t1;
      double t42 = 0.5 + t16 + t1;
      double t43 = t42 * t41;
      double t44 = -0.1546536707079771 + t16 + t1;
      double t48 = -t3 - t9 + 0.2650553239294647;
      double t49 = -t3 - t9 - 0.2147684835193549;
      double t51 = -t3 - t9 - 0.7852315164806451;
      double t52 = -t3 - t9 - 0.1265055323929465E1;
      double t54 = t52 * t51 * t49 * t48;
      double t57 = -t3 - t9 + 0.1546536707079771;
      double t58 = -t3 - t9 - 0.5;
      double t60 = -t3 - t9 - 0.1154653670707977E1;
      double t61 = t60 * t58 * t57;
      double t64 = t12 * t10;
      double t68 = 0.1265055323929465E1 + t16 + t1;
      double t69 = 0.7852315164806451 + t16 + t1;
      double t71 = 0.2147684835193549 + t16 + t1;
      double t72 = -0.2650553239294647 + t16 + t1;
      double t73 = t72 * t71;
      double t82 = t60 * t58;
      double t87 = (0.1330223896278567E1 + t16 + t1) * t6;
      double t88 = 0.9688487934707142 + t16 + t1;
      double t91 = 0.3115120652928579E-1 + t16 + t1;
      double t93 = -0.3302238962785669 + t16 + t1;
      double t114 = t68 * t6;
      double t121 = (0.1371740148509607E1 + t16 + t1) * t6;
      double t122 = 0.1091700181433142E1 + t16 + t1;
      double t125 = 0.7092992179024789 + t16 + t1;
      double t126 = 0.2907007820975211 + t16 + t1;
      double t128 = -0.917001814331423E-1 + t16 + t1;
      double t129 = -0.3717401485096066 + t16 + t1;
      double t135 = -t3 - t9 + 0.3717401485096066;
      double t146 =
        (-t3 - t9 - 0.1371740148509607E1) * (-t3 - t9 -
               0.1091700181433142E1) * (-t3 - t9 -
                      0.7092992179024789)
        * (-t3 - t9 - 0.2907007820975211) * (-t3 - t9 + 0.917001814331423E-1);
      double MapleGenVar1 =
        -0.7321104826850749E3 * t7 - 0.3039978717550737E4 * t12 * t10 * t6 * t5 -
        0.6391691498246073E4 * t21 - 0.7753863286709051E3 * t29 * t26 * t5 -
        0.3937610814400307E4 * t12 * t10 * t25 * t7 -
        0.6244535848343604E4 * t29 * t37 * t7 +
        0.1252982454838416E5 * t44 * t43 * t7 - 0.1584292313575895E4 * t54 * t7 +
        0.3110666290793476E4 * t61 * t7;
      double t149 =
        MapleGenVar1 - 0.233004131924886E5 * t64 * t37 * t7 -
        0.1343873984047569E5 * t73 * t69 * t68 * t7 +
        0.866754541388572E4 * t29 * t44 * t43 * t7 +
        0.3018539765005329E4 * t82 * t57 * t25 * t7 +
        0.427419277015135E5 * t12 * t10 * t93 * t91 * t42 * t88 * t87 * t5 +
        0.8009529691663972E4 * (-0.3997579954114602 + t16 +
              t1) * (-0.1771862795107378 + t16 +
               t1) * (0.1368825361738218 + t16 +
                t1) * t42 * (0.8631174638261782 +
                       t16 +
                       t1) *
        (0.1177186279510738E1 + t16 + t1) * (0.139975799541146E1 + t16 +
               t1) * t6 * t5 +
        0.3187110592391705E5 * t61 * t73 * t69 * t114 * t5 -
        0.9296948523935838E4 * t29 * t129 * t128 * t126 * t125 * t122 * t121 *
        t5 - 0.1565236214672771E3 * t146 * t135 * t26 * t5;
      double t165 = t41 * t6;
      double t170 = t52 * t51 * t49;
      double t174 = -t3 - t9 + 0.3302238962785669;
      double t175 = -t3 - t9 - 0.3115120652928579E-1;
      double t177 = -t3 - t9 - 0.9688487934707142;
      double t179 = -t3 - t9 - 0.1330223896278567E1;
      double t181 = t179 * t177 * t58 * t175 * t174;
      double t184 = t26 * t5;
      double t187 = t18 * t5;
      double t192 = t87 * t5;
      double t193 = t42 * t88;
      double t194 = t93 * t91;
      double t205 = t165 * t5;
      double t206 = t44 * t42;
      double t210 = t114 * t5;
      double t211 = t71 * t69;
      MapleGenVar1 =
        0.4885436707797789E2 * (-t3 - t9 - 0.139975799541146E1) * (-t3 - t9 -
                         0.1177186279510738E1)
        * (-t3 - t9 - 0.8631174638261782) * t58 * (-t3 - t9 -
                     0.1368825361738218) * (-t3 -
                          t9 +
                          0.1771862795107378)
        * (-t3 - t9 + 0.3997579954114602) * t6 * t5 +
        0.1322968769260516E5 * t170 * t48 * t44 * t42 * t165 * t5 +
        0.3198955594291514E4 * t181 * t21 - 0.1959192364970552E4 * t54 * t184 +
        0.2284393995262139E5 * t82 * t57 * t19 * t187 +
        0.1811889849856858E5 * t194 * t193 * t192 +
        0.7411664677842599E3 * t179 * t177 * t58 * t175 * t174 * t6 * t5 +
        0.4363063669972135E5 * t64 * t206 * t205 -
        0.136056785568102E5 * t29 * t72 * t211 * t210;
      double t249 =
        MapleGenVar1 -
        0.8273751188185237E4 * t129 * t128 * t126 * t125 * t122 * t121 * t5 -
        0.8488811532974639E2 * t146 * t135 * t6 * t5 -
        0.8424990827982554E4 * t170 * t48 * t19 * t187 -
        0.2277620371046502E5 * t61 * t206 * t205 +
        0.8168906203729521E4 * t29 * t194 * t193 * t192 -
        0.3340404201829962E5 * t12 * t10 * t72 * t211 * t210 +
        0.4770905375469706E3 * t181 * t184 - 0.9621911560865351E3 * t184 -
        0.713011053183776E3 * t29 * t6 * t5;
      double t250 = t149 + t249;
      return t250;
    }

    double
    ortho2_f66x (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t5 = (-t3 - t1) * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t8 = 0.1E1 * x;
      double t9 = 0.1154653670707977E1 + t8 + t1;
      double t10 = 0.5 + t8 + t1;
      double t11 = t10 * t9;
      double t12 = 0.1E1 * y;
      double t13 = -t3 - t12 - 0.5278640450004206E-1;
      double t14 = -t3 - t12 - 0.9472135954999579;
      double t15 = t14 * t13;
      double t19 = -0.1546536707079771 + t8 + t1;
      double t20 = t19 * t9;
      double t24 = t19 * t10;
      double t25 = t15 * t24;
      double t28 = t9 * t5;
      double t31 = t6 * t2;
      double t32 = t9 * t31;
      double t35 = -t3 - t12 + 0.2650553239294647;
      double t36 = -t3 - t12 - 0.2147684835193549;
      double t37 = t36 * t35;
      double t38 = -t3 - t12 - 0.7852315164806451;
      double t39 = -t3 - t12 - 0.1265055323929465E1;
      double t40 = t39 * t38;
      double t41 = t40 * t37;
      double t52 = -t3 - t12 + 0.3302238962785669;
      double t54 = -t3 - t12 - 0.3115120652928579E-1;
      double t55 = -t3 - t12 - 0.5;
      double t56 = t55 * t54;
      double t57 = -t3 - t12 - 0.9688487934707142;
      double t58 = -t3 - t12 - 0.1330223896278567E1;
      double t59 = t58 * t57;
      double t60 = t59 * t56;
      double t63 = 0.1265055323929465E1 + t8 + t1;
      double t64 = 0.7852315164806451 + t8 + t1;
      double t65 = t64 * t63;
      double t66 = 0.2147684835193549 + t8 + t1;
      double t67 = -0.2650553239294647 + t8 + t1;
      double t68 = t67 * t66;
      double t78 = t59 * t55 * t52;
      double t81 = 0.1330223896278567E1 + t8 + t1;
      double t83 = 0.3115120652928579E-1 + t8 + t1;
      double t84 = -0.3302238962785669 + t8 + t1;
      double t85 = t84 * t83;
      double t89 =
        0.4363063669972135E5 * t15 * t11 * t7 +
        0.4363063669972135E5 * t15 * t20 * t7 + 0.4363063669972135E5 * t25 * t7 +
        0.2181531834986068E5 * t25 * t28 - 0.2181531834986068E5 * t25 * t32 +
        0.6195510247718829E4 * t41 * t7 -
        0.2181531834986068E5 * t14 * t19 * t11 * t7 -
        0.2181531834986068E5 * t13 * t19 * t11 * t7 -
        0.3705832338921299E3 * t60 * t52 * t31 -
        0.2151246667581642E5 * t68 * t65 * t7 - 0.3705832338921299E3 * t60 * t7 +
        0.3705832338921299E3 * t60 * t52 * t5 - 0.3705832338921299E3 * t78 * t7 +
        0.1811889849856858E5 * t85 * t10 * t81 * t7;
      double t90 = 0.9688487934707142 + t8 + t1;
      double t91 = t10 * t90;
      double t92 = t85 * t91;
      double t101 = t54 * t52;
      double t102 = t59 * t101;
      double t105 = t90 * t81;
      double t113 = t83 * t10;
      double t118 = t58 * t55 * t101;
      double t121 = t57 * t55;
      double t122 = t121 * t101;
      double t127 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t132 = t127 * t5;
      double t135 = t127 * t31;
      double t138 = t35 * t127;
      double t142 = t39 * t36;
      double t146 = t38 * t36;
      double t150 =
        0.1811889849856858E5 * t92 * t7 + 0.905944924928429E4 * t81 * t5 * t92 -
        0.905944924928429E4 * t92 * t81 * t31 - 0.3705832338921299E3 * t102 * t7 +
        0.1811889849856858E5 * t85 * t105 * t7 +
        0.1811889849856858E5 * t84 * t10 * t105 * t7 +
        0.1811889849856858E5 * t113 * t105 * t7 -
        0.3705832338921299E3 * t118 * t7 - 0.3705832338921299E3 * t122 * t7 +
        0.9795961824852758E3 * t40 * t36 * t127 * t7 -
        0.9795961824852758E3 * t41 * t132 + 0.9795961824852758E3 * t41 * t135 +
        0.9795961824852758E3 * t40 * t138 * t7 +
        0.9795961824852758E3 * t142 * t138 * t7 +
        0.9795961824852758E3 * t146 * t138 * t7;
      double t153 = t39 * t146;
      double t158 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t162 = 0.5278640450004206E-1 + t8 + t1;
      double t163 = -t3 - t12 + 0.1546536707079771;
      double t165 = -t3 - t12 - 0.1154653670707977E1;
      double t166 = t165 * t55;
      double t167 = t166 * t163 * t162;
      double t170 = t158 * t24;
      double t173 = 0.9472135954999579 + t8 + t1;
      double t174 = t173 * t5;
      double t180 = t173 * t31;
      double t183 = t55 * t163;
      double t184 = t165 * t183;
      double t190 = t19 * t11;
      double t193 = t38 * t37;
      double t196 = t39 * t37;
      double t200 = t39 * t38 * t35;
      double t205 = t162 * t173;
      double t209 =
        -0.7921461567879473E3 * t153 * t35 * t5 +
        0.866754541388572E4 * t158 * t20 * t7 + 0.2284393995262139E5 * t167 * t7 +
        0.866754541388572E4 * t170 * t7 + 0.1142196997631069E5 * t167 * t174 +
        0.7921461567879473E3 * t153 * t35 * t31 -
        0.1142196997631069E5 * t167 * t180 - 0.9545460865206262E4 * t184 * t7 +
        0.866754541388572E4 * t158 * t11 * t7 + 0.1370459261541285E5 * t190 * t7 +
        0.7921461567879473E3 * t193 * t7 + 0.7921461567879473E3 * t196 * t7 +
        0.7921461567879473E3 * t200 * t7 + 0.7921461567879473E3 * t153 * t7 -
        0.1142196997631069E5 * t166 * t205 * t7;
      double t217 = t165 * t163;
      double t225 = t158 * t67;
      double t229 = t66 * t63;
      double t233 = t66 * t64;
      double t234 = t225 * t233;
      double t237 = t63 * t5;
      double t240 = t63 * t31;
      double t256 = t14 * t13 * t162;
      double t259 =
        0.2284393995262139E5 * t166 * t163 * t173 * t7 -
        0.1142196997631069E5 * t183 * t205 * t7 -
        0.1142196997631069E5 * t217 * t205 * t7 -
        0.136056785568102E5 * t158 * t66 * t65 * t7 -
        0.136056785568102E5 * t225 * t65 * t7 -
        0.136056785568102E5 * t225 * t229 * t7 - 0.136056785568102E5 * t234 * t7 -
        0.6802839278405098E4 * t234 * t237 + 0.6802839278405098E4 * t234 * t240 -
        0.433377270694286E4 * t170 * t32 - 0.1509269882502665E4 * t184 * t135 +
        0.433377270694286E4 * t170 * t28 -
        0.1509269882502665E4 * t165 * t55 * t127 * t7 +
        0.1509269882502665E4 * t184 * t132 + 0.116502065962443E5 * t256 * t180;
      double t262 = t163 * t127;
      double t289 = t67 * t233;
      double t301 =
        -0.1509269882502665E4 * t55 * t262 * t7 -
        0.1509269882502665E4 * t165 * t262 * t7 -
        0.116502065962443E5 * t256 * t174 +
        0.116502065962443E5 * t14 * t205 * t7 -
        0.233004131924886E5 * t14 * t13 * t173 * t7 -
        0.233004131924886E5 * t256 * t7 - 0.1343873984047569E5 * t66 * t65 * t7 -
        0.1343873984047569E5 * t67 * t65 * t7 -
        0.1343873984047569E5 * t67 * t229 * t7 -
        0.1343873984047569E5 * t289 * t7 - 0.6719369920237843E4 * t289 * t237 +
        0.6719369920237843E4 * t289 * t240 +
        0.116502065962443E5 * t13 * t205 * t7 - 0.4810955780432676E3 * t132 +
        0.4810955780432676E3 * t135;
      double t307 = t105 * t31;
      double t308 = t158 * t84;
      double t309 = t308 * t113;
      double t312 = t105 * t5;
      double t316 = t90 * t6 * t5;
      double t320 = t81 * t6 * t5;
      double t323 = t83 * t90;
      double t335 = t63 * t6 * t5;
      double t336 = t163 * t67;
      double t348 = t14 * t13 * t84;
      double t349 = t348 * t113;
      double t354 = -0.356505526591888E3 * t158 * t5 + 0.1915346135235945E4 * t7 +
        0.356505526591888E3 * t158 * t31 - 0.408445310186476E4 * t309 * t307 +
        0.408445310186476E4 * t309 * t312 + 0.8168906203729521E4 * t309 * t316 +
        0.8168906203729521E4 * t309 * t320 +
        0.8168906203729521E4 * t308 * t323 * t320 +
        0.8168906203729521E4 * t308 * t91 * t320 +
        0.8168906203729521E4 * t158 * t83 * t91 * t320 -
        0.1593555296195853E5 * t165 * t336 * t233 * t335 -
        0.1593555296195853E5 * t55 * t336 * t233 * t335 +
        0.1291617479803237E5 * t92 * t320 - 0.2137096385075675E5 * t349 * t307 +
        0.427419277015135E5 * t349 * t316;
      double t379 = -t3 - t12 + 0.3717401485096066;
      double t380 = t379 * t127;
      double t382 = -t3 - t12 + 0.917001814331423E-1;
      double t383 = -t3 - t12 - 0.2907007820975211;
      double t384 = t383 * t382;
      double t385 = -t3 - t12 - 0.7092992179024789;
      double t386 = -t3 - t12 - 0.1091700181433142E1;
      double t387 = t386 * t385;
      double t388 = -t3 - t12 - 0.1371740148509607E1;
      double t389 = t388 * t387;
      double t390 = t389 * t384;
      double t396 = t127 * t6;
      double t397 = t396 * t5;
      double t404 = t382 * t379;
      double t413 = t385 * t383;
      double t422 =
        0.2137096385075675E5 * t349 * t312 +
        0.427419277015135E5 * t348 * t323 * t320 +
        0.427419277015135E5 * t349 * t320 +
        0.427419277015135E5 * t14 * t13 * t83 * t91 * t320 +
        0.427419277015135E5 * t348 * t91 * t320 -
        0.2137096385075675E5 * t14 * t85 * t91 * t320 -
        0.2137096385075675E5 * t13 * t85 * t91 * t320 -
        0.7826181073363854E2 * t390 * t380 * t5 +
        0.7826181073363854E2 * t390 * t380 * t31 +
        0.7826181073363854E2 * t390 * t397 +
        0.7826181073363854E2 * t389 * t383 * t379 * t397 +
        0.7826181073363854E2 * t389 * t404 * t397 +
        0.7826181073363854E2 * t388 * t386 * t383 * t404 * t397 +
        0.7826181073363854E2 * t386 * t413 * t404 * t397 +
        0.7826181073363854E2 * t388 * t413 * t404 * t397;
      double t423 = 0.1091700181433142E1 + t8 + t1;
      double t425 = t423 * t6 * t5;
      double t426 = 0.7092992179024789 + t8 + t1;
      double t427 = 0.2907007820975211 + t8 + t1;
      double t428 = t427 * t426;
      double t429 = -0.917001814331423E-1 + t8 + t1;
      double t430 = -0.3717401485096066 + t8 + t1;
      double t431 = t430 * t429;
      double t432 = t431 * t428;
      double t435 = 0.1371740148509607E1 + t8 + t1;
      double t436 = t423 * t435;
      double t437 = t436 * t5;
      double t440 = t436 * t31;
      double t444 = t379 * t6 * t5;
      double t448 = t435 * t6 * t5;
      double t451 = -t3 - t12 + 0.3997579954114602;
      double t452 = -t3 - t12 + 0.1771862795107378;
      double t453 = t452 * t451;
      double t455 = -t3 - t12 - 0.1368825361738218;
      double t457 = -t3 - t12 - 0.8631174638261782;
      double t458 = -t3 - t12 - 0.1177186279510738E1;
      double t460 = -t3 - t12 - 0.139975799541146E1;
      double t461 = t460 * t458 * t457;
      double t462 = t461 * t55 * t455;
      double t465 = t427 * t423;
      double t469 = t65 * t31;
      double t470 = t184 * t68;
      double t473 = t426 * t423;
      double t481 = t451 * t6 * t5;
      double t482 = t455 * t452;
      double t486 = t65 * t5;
      double t499 =
        -0.8273751188185237E4 * t432 * t425 - 0.4136875594092619E4 * t432 * t437 +
        0.4136875594092619E4 * t432 * t440 + 0.4949711514546221E3 * t390 * t444 -
        0.8273751188185237E4 * t432 * t448 -
        0.2442718353898894E2 * t462 * t453 * t31 -
        0.8273751188185237E4 * t431 * t465 * t448 -
        0.1593555296195853E5 * t470 * t469 -
        0.8273751188185237E4 * t431 * t473 * t448 +
        0.2442718353898894E2 * t462 * t453 * t5 -
        0.2442718353898894E2 * t461 * t482 * t481 +
        0.1593555296195853E5 * t470 * t486 -
        0.2442718353898894E2 * t461 * t55 * t452 * t481 -
        0.2442718353898894E2 * t462 * t481 -
        0.2442718353898894E2 * t462 * t452 * t6 * t5;
      double t503 = t10 * t6;
      double t504 = t503 * t5;
      double t506 = t166 * t163 * t19;
      double t509 = t11 * t5;
      double t512 = t11 * t31;
      double t520 = t9 * t6;
      double t521 = t520 * t5;
      double t524 = t457 * t55;
      double t534 = t64 * t6 * t5;
      double t541 = t429 * t427;
      double t548 = t388 * t386;
      double t549 = t548 * t413;
      double t565 =
        -0.2277620371046502E5 * t506 * t504 - 0.1138810185523251E5 * t506 * t509 +
        0.1138810185523251E5 * t506 * t512 -
        0.2442718353898894E2 * t460 * t458 * t55 * t482 * t481 -
        0.2277620371046502E5 * t506 * t521 -
        0.2442718353898894E2 * t460 * t524 * t482 * t481 -
        0.2277620371046502E5 * t166 * t163 * t10 * t521 +
        0.3187110592391705E5 * t470 * t534 -
        0.2442718353898894E2 * t458 * t524 * t482 * t481 -
        0.1469976631248884E5 * t430 * t541 * t473 * t448 +
        0.4244405766487319E2 * t549 * t382 * t6 * t5 -
        0.4244405766487319E2 * t549 * t404 * t5 +
        0.4244405766487319E2 * t549 * t404 * t31 -
        0.1508692748762508E4 * t60 * t52 * t6 * t5 +
        0.1138810185523251E5 * t166 * t24 * t521;
      double t572 = t15 * t68;
      double t577 = t430 * t427;
      double t588 = t67 * t64;
      double t612 =
        0.1138810185523251E5 * t183 * t24 * t521 +
        0.1138810185523251E5 * t217 * t24 * t521 -
        0.1670202100914981E5 * t572 * t486 + 0.1670202100914981E5 * t572 * t469 -
        0.8273751188185237E4 * t577 * t473 * t448 -
        0.3340404201829962E5 * t572 * t335 - 0.3340404201829962E5 * t572 * t534 -
        0.8273751188185237E4 * t541 * t473 * t448 -
        0.3340404201829962E5 * t15 * t588 * t335 -
        0.3340404201829962E5 * t15 * t233 * t335 +
        0.1670202100914981E5 * t13 * t67 * t233 * t335 +
        0.1670202100914981E5 * t14 * t67 * t233 * t335 -
        0.6614843846302581E4 * t193 * t24 * t521 -
        0.6614843846302581E4 * t196 * t24 * t521 -
        0.6614843846302581E4 * t200 * t24 * t521;
      double t622 = t153 * t35 * t19;
      double t633 = t52 * t127;
      double t643 = 0.139975799541146E1 + t8 + t1;
      double t645 = t643 * t6 * t5;
      double t646 = 0.8631174638261782 + t8 + t1;
      double t648 = 0.1368825361738218 + t8 + t1;
      double t649 = -0.1771862795107378 + t8 + t1;
      double t651 = -0.3997579954114602 + t8 + t1;
      double t652 = t651 * t649 * t648;
      double t653 = t652 * t10 * t646;
      double t660 = 0.1177186279510738E1 + t8 + t1;
      double t672 =
        -0.6614843846302581E4 * t153 * t24 * t521 +
        0.1322968769260516E5 * t153 * t35 * t10 * t521 +
        0.1322968769260516E5 * t622 * t521 + 0.1322968769260516E5 * t622 * t504 +
        0.6614843846302581E4 * t622 * t509 - 0.6614843846302581E4 * t622 * t512 +
        0.3187110592391705E5 * t470 * t335 +
        0.2385452687734853E3 * t60 * t633 * t5 -
        0.2385452687734853E3 * t60 * t633 * t31 +
        0.4244405766487319E2 * t387 * t384 * t444 +
        0.8009529691663972E4 * t653 * t645 +
        0.4244405766487319E2 * t388 * t385 * t384 * t444 +
        0.8009529691663972E4 * t653 * t660 * t6 * t5 +
        0.4244405766487319E2 * t548 * t384 * t444 +
        0.4244405766487319E2 * t548 * t385 * t382 * t444;
      double t673 = t660 * t643;
      double t693 = t646 * t660;
      double t694 = t648 * t10;
      double t717 = t162 * t6;
      double t718 = t717 * t5;
      double t721 =
        0.4004764845831986E4 * t653 * t673 * t5 +
        0.4244405766487319E2 * t549 * t444 -
        0.4004764845831986E4 * t653 * t673 * t31 +
        0.3187110592391705E5 * t184 * t588 * t335 -
        0.2385452687734853E3 * t102 * t397 - 0.2385452687734853E3 * t78 * t397 -
        0.2385452687734853E3 * t60 * t397 - 0.2385452687734853E3 * t118 * t397 +
        0.8009529691663972E4 * t649 * t694 * t693 * t645 +
        0.8009529691663972E4 * t651 * t694 * t693 * t645 +
        0.8009529691663972E4 * t651 * t649 * t10 * t693 * t645 +
        0.8009529691663972E4 * t652 * t693 * t645 +
        0.8009529691663972E4 * t652 * t10 * t660 * t645 -
        0.2385452687734853E3 * t122 * t397 - 0.8424990827982554E4 * t41 * t718;
      double t724 = t205 * t5;
      double t727 = t205 * t31;
      double t730 = t173 * t6;
      double t731 = t730 * t5;
      double t738 = t58 * t121;
      double t739 = t738 * t101;
      double t747 = t35 * t162;
      double t765 = t52 * t162;
      double t774 =
        -0.4212495413991277E4 * t41 * t724 + 0.4212495413991277E4 * t41 * t727 +
        0.4212495413991277E4 * t40 * t36 * t162 * t731 -
        0.8424990827982554E4 * t41 * t731 - 0.1599477797145757E4 * t739 * t727 +
        0.3187110592391705E5 * t184 * t233 * t335 +
        0.1599477797145757E4 * t739 * t724 +
        0.4212495413991277E4 * t146 * t747 * t731 +
        0.4212495413991277E4 * t142 * t747 * t731 -
        0.1599477797145757E4 * t738 * t54 * t162 * t731 +
        0.3198955594291514E4 * t739 * t731 + 0.3198955594291514E4 * t739 * t718 +
        0.4212495413991277E4 * t40 * t747 * t731 -
        0.1599477797145757E4 * t57 * t56 * t765 * t731 -
        0.1599477797145757E4 * t58 * t56 * t765 * t731;
      double t783 = t158 * t431;
      double t784 = t783 * t428;
      double t815 = t13 * t6;
      double t824 =
        -0.1599477797145757E4 * t58 * t57 * t54 * t765 * t731 -
        0.1599477797145757E4 * t738 * t765 * t731 -
        0.9296948523935838E4 * t784 * t425 - 0.4648474261967919E4 * t784 * t437 +
        0.4648474261967919E4 * t784 * t440 - 0.9296948523935838E4 * t784 * t448 -
        0.9296948523935838E4 * t783 * t465 * t448 -
        0.1593555296195853E5 * t165 * t55 * t67 * t233 * t335 -
        0.9296948523935838E4 * t783 * t473 * t448 -
        0.9296948523935838E4 * t158 * t577 * t473 * t448 -
        0.9296948523935838E4 * t158 * t541 * t473 * t448 -
        0.9873478105668788E4 * t162 * t730 * t5 +
        0.1245181871281551E5 * t14 * t815 * t5 +
        0.1519989358775368E4 * t14 * t6 * t5 + 0.1519989358775368E4 * t815 * t5;
      double t829 = t158 * t127;
      double t839 = t158 * t205;
      double t856 = t163 * t6;
      double t860 =
        0.245198686515598E4 * t158 * t6 * t5 + 0.3876931643354526E3 * t829 * t31 -
        0.122599343257799E4 * t397 - 0.3876931643354526E3 * t829 * t5 +
        0.3195845749123036E4 * t727 - 0.6391691498246073E4 * t731 -
        0.6391691498246073E4 * t718 - 0.3195845749123036E4 * t724 +
        0.3122267924171802E4 * t839 * t31 - 0.3122267924171802E4 * t839 * t5 -
        0.6244535848343604E4 * t158 * t717 * t5 -
        0.6244535848343604E4 * t158 * t730 * t5 -
        0.626491227419208E4 * t190 * t31 -
        0.1555333145396738E4 * t165 * t55 * t6 * t5 -
        0.1555333145396738E4 * t165 * t856 * t5;
      double t880 = t14 * t13 * t127;
      double t897 =
        -0.1555333145396738E4 * t184 * t31 + 0.1555333145396738E4 * t184 * t5 -
        0.1555333145396738E4 * t55 * t856 * t5 + 0.626491227419208E4 * t190 * t5 +
        0.1252982454838416E5 * t19 * t520 * t5 +
        0.1252982454838416E5 * t19 * t503 * t5 +
        0.1252982454838416E5 * t10 * t520 * t5 +
        0.1968805407200153E4 * t880 * t31 - 0.1968805407200153E4 * t880 * t5 +
        0.1968805407200153E4 * t14 * t396 * t5 +
        0.1968805407200153E4 * t13 * t396 * t5 - 0.1519989358775368E4 * t15 * t5 +
        0.1519989358775368E4 * t15 * t31 + 0.3660552413425374E3 * t31 -
        0.3660552413425374E3 * t5;
      double t901 =
        t89 + t150 + t209 + t259 + t301 + t354 + t422 + t499 + t565 + t612 +
        t672 + t721 + t774 + t824 + t860 + t897;
      return t901;
    }

    double
    ortho2_f66y (double x, double y)
    {
      double t1 = 0.5 * y;
      double t2 = 0.5 + t1;
      double t3 = 0.5 * x;
      double t4 = -t3 - t1;
      double t5 = t4 * t2;
      double t6 = 0.5 + t3;
      double t7 = t6 * t5;
      double t10 =
        -0.158113883008419E1 - 0.3162277660168379E1 * x - 0.158113883008419E1 * y;
      double t11 = 0.1E1 * y;
      double t12 = -t3 - t11 + 0.2650553239294647;
      double t13 = t12 * t10;
      double t14 = -t3 - t11 - 0.7852315164806451;
      double t15 = -t3 - t11 - 0.1265055323929465E1;
      double t16 = t15 * t14;
      double t20 = -t3 - t11 - 0.2147684835193549;
      double t25 = 0.1E1 * x;
      double t26 = 0.1265055323929465E1 + t25 + t1;
      double t27 = 0.7852315164806451 + t25 + t1;
      double t28 = t27 * t26;
      double t29 = 0.2147684835193549 + t25 + t1;
      double t30 = -0.2650553239294647 + t25 + t1;
      double t31 = t30 * t29;
      double t35 = 0.1371740148509607E1 + t25 + t1;
      double t37 = t35 * t6 * t5;
      double t38 = 0.1091700181433142E1 + t25 + t1;
      double t39 = 0.2907007820975211 + t25 + t1;
      double t40 = t39 * t38;
      double t41 = -0.917001814331423E-1 + t25 + t1;
      double t42 = -0.3717401485096066 + t25 + t1;
      double t43 = t42 * t41;
      double t47 = t29 * t26;
      double t51 = 0.9472135954999579 + t25 + t1;
      double t52 = t51 * t6;
      double t53 = 0.5278640450004206E-1 + t25 + t1;
      double t57 = t15 * t20;
      double t61 = t14 * t20;
      double t65 = t6 * t4;
      double t66 = -t3 - t11 - 0.5278640450004206E-1;
      double t68 = -t3 - t11 - 0.9472135954999579;
      double t69 = t68 * t66 * t10;
      double t72 = 0.139975799541146E1 + t25 + t1;
      double t74 = t72 * t6 * t5;
      double t75 = 0.1177186279510738E1 + t25 + t1;
      double t76 = 0.8631174638261782 + t25 + t1;
      double t77 = t76 * t75;
      double t78 = 0.5 + t25 + t1;
      double t79 = -0.1771862795107378 + t25 + t1;
      double t81 = -0.3997579954114602 + t25 + t1;
      double t86 = t66 * t6;
      double t90 = t53 * t6;
      double t91 = t90 * t5;
      double t92 = t20 * t12;
      double t93 = t16 * t92;
      double t97 = t26 * t6 * t5;
      double t98 = t30 * t27;
      double t99 = -t3 - t11 + 0.1546536707079771;
      double t100 = -t3 - t11 - 0.5;
      double t101 = t100 * t99;
      double t102 = -t3 - t11 - 0.1154653670707977E1;
      double t103 = t102 * t101;
      double t107 = t103 * t31;
      double t110 =
        0.1959192364970552E4 * t16 * t13 * t7 +
        0.1959192364970552E4 * t16 * t20 * t10 * t7 -
        0.4302493335163284E5 * t31 * t28 * t7 -
        0.4136875594092619E4 * t43 * t40 * t37 -
        0.6719369920237843E4 * t30 * t47 * t7 -
        0.1974695621133758E5 * t53 * t52 * t5 +
        0.1959192364970552E4 * t57 * t13 * t7 +
        0.1959192364970552E4 * t61 * t13 * t7 - 0.1968805407200153E4 * t69 * t65 +
        0.4004764845831986E4 * t81 * t79 * t78 * t77 * t74 +
        0.6225909356407754E4 * t68 * t86 * t5 - 0.4212495413991277E4 * t93 * t91 +
        0.1593555296195853E5 * t103 * t98 * t97 +
        0.1593555296195853E5 * t107 * t97;
      double t111 = 0.1330223896278567E1 + t25 + t1;
      double t112 = 0.9688487934707142 + t25 + t1;
      double t113 = t112 * t111;
      double t114 = t113 * t65;
      double t115 = 0.3115120652928579E-1 + t25 + t1;
      double t116 = t115 * t78;
      double t117 = -0.3302238962785669 + t25 + t1;
      double t119 = t68 * t66 * t117;
      double t120 = t119 * t116;
      double t123 = 0.1154653670707977E1 + t25 + t1;
      double t124 = t123 * t6;
      double t125 = t124 * t5;
      double t126 = -0.1546536707079771 + t25 + t1;
      double t127 = t126 * t78;
      double t128 = t14 * t92;
      double t132 = t15 * t92;
      double t136 = -t3 - t11 + 0.3717401485096066;
      double t138 = t136 * t6 * t5;
      double t139 = -t3 - t11 + 0.917001814331423E-1;
      double t140 = -t3 - t11 - 0.2907007820975211;
      double t141 = t140 * t139;
      double t142 = -t3 - t11 - 0.7092992179024789;
      double t143 = -t3 - t11 - 0.1091700181433142E1;
      double t144 = t143 * t142;
      double t145 = -t3 - t11 - 0.1371740148509607E1;
      double t146 = t145 * t144;
      double t147 = t146 * t141;
      double t151 = t15 * t14 * t12;
      double t155 = t15 * t61;
      double t164 = t155 * t12 * t126;
      double t167 = t78 * t6;
      double t168 = t167 * t5;
      double t171 = t6 * t2;
      double t172 = t78 * t123;
      double t173 = t172 * t171;
      double t176 = 0.7092992179024789 + t25 + t1;
      double t177 = t176 * t38;
      double t183 =
        0.158113883008419E1 * x + 0.3162277660168379E1 * y + 0.158113883008419E1;
      double t184 = t183 * t117;
      double t185 = t184 * t116;
      double t189 = t111 * t6 * t5;
      double t190 = t78 * t112;
      double t191 = t117 * t115;
      double t192 = t191 * t190;
      double t195 = t53 * t51;
      double t196 = t183 * t195;
      double t199 = -t3 - t11 + 0.3302238962785669;
      double t201 = -t3 - t11 - 0.3115120652928579E-1;
      double t202 = t100 * t201;
      double t203 = -t3 - t11 - 0.9688487934707142;
      double t204 = -t3 - t11 - 0.1330223896278567E1;
      double t205 = t204 * t203;
      double t206 = t205 * t202;
      double t209 =
        0.2137096385075675E5 * t120 * t114 -
        0.1322968769260516E5 * t128 * t127 * t125 -
        0.1322968769260516E5 * t132 * t127 * t125 +
        0.2474855757273111E3 * t147 * t138 -
        0.1322968769260516E5 * t151 * t127 * t125 -
        0.1322968769260516E5 * t155 * t127 * t125 +
        0.6614843846302581E4 * t155 * t12 * t78 * t125 +
        0.6614843846302581E4 * t164 * t125 + 0.6614843846302581E4 * t164 * t168 -
        0.6614843846302581E4 * t164 * t173 -
        0.4136875594092619E4 * t43 * t177 * t37 +
        0.408445310186476E4 * t185 * t114 + 0.2583234959606475E5 * t192 * t189 +
        0.3122267924171802E4 * t196 * t171 -
        0.3705832338921299E3 * t206 * t199 * t171;
      double t211 = t126 * t172;
      double t221 = 0.1368825361738218 + t25 + t1;
      double t223 = t81 * t79 * t221;
      double t224 = t223 * t78 * t76;
      double t227 = t51 * t65;
      double t229 = t68 * t66 * t53;
      double t236 = t75 * t72;
      double t245 = t52 * t5;
      double t248 = t10 * t171;
      double t255 = t183 * t10;
      double t258 = t10 * t6;
      double t259 = t258 * t5;
      double t261 = t201 * t199;
      double t263 = t204 * t100 * t261;
      double t266 =
        -0.626491227419208E4 * t211 * t171 -
        0.3122267924171802E4 * t183 * t52 * t5 -
        0.6719369920237843E4 * t30 * t28 * t7 +
        0.4004764845831986E4 * t224 * t74 - 0.116502065962443E5 * t229 * t227 +
        0.4004764845831986E4 * t224 * t75 * t6 * t5 -
        0.4004764845831986E4 * t224 * t236 * t171 -
        0.3122267924171802E4 * t183 * t90 * t5 + 0.274091852308257E5 * t211 * t7 -
        0.4212495413991277E4 * t93 * t245 - 0.1509269882502665E4 * t103 * t248 +
        0.8424990827982554E4 * t16 * t20 * t53 * t245 +
        0.3876931643354526E3 * t255 * t171 - 0.245198686515598E4 * t259 -
        0.4770905375469706E3 * t263 * t259;
      double t270 = t68 * t66;
      double t285 = t195 * t171;
      double t289 = t12 * t53;
      double t293 = t99 * t10;
      double t304 = t142 * t140;
      double t305 = t145 * t143;
      double t306 = t305 * t304;
      double t309 =
        0.122599343257799E4 * t183 * t6 * t5 - 0.1519989358775368E4 * t270 * t65 -
        0.3018539765005329E4 * t102 * t100 * t10 * t7 +
        0.3039978717550737E4 * t86 * t5 + 0.3039978717550737E4 * t68 * t6 * t5 +
        0.1519989358775368E4 * t270 * t171 - 0.3195845749123036E4 * t245 +
        0.3195845749123036E4 * t285 - 0.3876931643354526E3 * t255 * t65 +
        0.8424990827982554E4 * t16 * t289 * t245 -
        0.3018539765005329E4 * t102 * t293 * t7 +
        0.433377270694286E4 * t183 * t172 * t7 -
        0.6719369920237843E4 * t29 * t28 * t7 - 0.3195845749123036E4 * t91 +
        0.8488811532974639E2 * t306 * t138;
      double t333 = t126 * t123;
      double t340 = t183 * t127;
      double t344 = t205 * t100 * t199;
      double t347 = -t3 - t11 + 0.3997579954114602;
      double t348 = -t3 - t11 + 0.1771862795107378;
      double t349 = t348 * t347;
      double t351 = -t3 - t11 - 0.1368825361738218;
      double t353 = -t3 - t11 - 0.8631174638261782;
      double t354 = -t3 - t11 - 0.1177186279510738E1;
      double t356 = -t3 - t11 - 0.139975799541146E1;
      double t357 = t356 * t354 * t353;
      double t358 = t357 * t100 * t351;
      double t371 =
        -0.4772730432603131E4 * t103 * t7 +
        0.8488811532974639E2 * t305 * t142 * t139 * t138 -
        0.7411664677842599E3 * t206 * t7 -
        0.3110666290793476E4 * t102 * t100 * t6 * t5 +
        0.905944924928429E4 * t192 * t111 * t65 +
        0.8488811532974639E2 * t305 * t141 * t138 +
        0.7921461567879473E3 * t155 * t12 * t171 +
        0.433377270694286E4 * t183 * t333 * t7 +
        0.8424990827982554E4 * t57 * t289 * t245 +
        0.433377270694286E4 * t340 * t7 - 0.7411664677842599E3 * t344 * t7 +
        0.2442718353898894E2 * t358 * t349 * t65 +
        0.8488811532974639E2 * t145 * t142 * t141 * t138 -
        0.7921461567879473E3 * t155 * t12 * t65 +
        0.8424990827982554E4 * t61 * t289 * t245;
      double t372 = t99 * t6;
      double t376 = t203 * t100;
      double t377 = t376 * t261;
      double t386 = t102 * t100;
      double t387 = t386 * t99 * t126;
      double t390 = t102 * t99;
      double t401 = t205 * t261;
      double t404 = t199 * t10;
      double t412 = t38 * t35;
      double t413 = t412 * t65;
      double t414 = t39 * t176;
      double t415 = t43 * t414;
      double t420 =
        -0.3110666290793476E4 * t102 * t372 * t5 -
        0.4770905375469706E3 * t377 * t259 +
        0.8488811532974639E2 * t144 * t141 * t138 +
        0.626491227419208E4 * t211 * t65 + 0.1138810185523251E5 * t387 * t173 -
        0.2284393995262139E5 * t390 * t195 * t7 +
        0.1584292313575895E4 * t151 * t7 -
        0.3018539765005329E4 * t100 * t293 * t7 +
        0.1584292313575895E4 * t155 * t7 - 0.7411664677842599E3 * t401 * t7 +
        0.2385452687734853E3 * t206 * t404 * t65 -
        0.1138810185523251E5 * t387 * t125 - 0.1138810185523251E5 * t387 * t168 -
        0.4136875594092619E4 * t415 * t413 + 0.1584292313575895E4 * t128 * t7;
      double t429 = t29 * t27;
      double t433 = t195 * t65;
      double t434 = t204 * t376;
      double t435 = t434 * t261;
      double t438 = t172 * t65;
      double t447 = t10 * t65;
      double t450 = t221 * t78;
      double t467 =
        -0.1555333145396738E4 * t103 * t171 + 0.1584292313575895E4 * t132 * t7 +
        0.4004764845831986E4 * t224 * t236 * t65 +
        0.1593555296195853E5 * t103 * t429 * t97 +
        0.1599477797145757E4 * t435 * t433 + 0.6614843846302581E4 * t164 * t438 -
        0.1138810185523251E5 * t386 * t99 * t78 * t125 -
        0.1599477797145757E4 * t435 * t285 - 0.9795961824852758E3 * t93 * t447 +
        0.4004764845831986E4 * t81 * t450 * t77 * t74 +
        0.1599477797145757E4 * t435 * t91 -
        0.2442718353898894E2 * t358 * t349 * t171 +
        0.1968805407200153E4 * t69 * t171 -
        0.905944924928429E4 * t192 * t111 * t171 -
        0.4212495413991277E4 * t93 * t433;
      double t468 = t28 * t171;
      double t473 = t26 * t171;
      double t474 = t183 * t30;
      double t475 = t474 * t429;
      double t496 = t199 * t53;
      double t504 = t347 * t6 * t5;
      double t511 = t139 * t136;
      double t519 =
        -0.1593555296195853E5 * t107 * t468 + 0.905944924928429E4 * t192 * t7 +
        0.6802839278405098E4 * t475 * t473 + 0.1599477797145757E4 * t435 * t245 -
        0.4885436707797789E2 * t358 * t348 * t6 * t5 -
        0.2385452687734853E3 * t206 * t404 * t171 -
        0.3198955594291514E4 * t434 * t201 * t53 * t245 +
        0.2137096385075675E5 * t68 * t66 * t115 * t190 * t189 -
        0.3198955594291514E4 * t434 * t496 * t245 +
        0.3937610814400307E4 * t68 * t258 * t5 -
        0.4885436707797789E2 * t357 * t100 * t348 * t504 -
        0.4885436707797789E2 * t358 * t504 +
        0.4244405766487319E2 * t306 * t511 * t171 -
        0.6802839278405098E4 * t475 * t7 - 0.4770905375469706E3 * t206 * t259;
      double t567 = t351 * t348;
      double t576 =
        -0.3198955594291514E4 * t204 * t202 * t496 * t245 +
        0.8488811532974639E2 * t306 * t139 * t6 * t5 +
        0.905944924928429E4 * t191 * t78 * t111 * t7 -
        0.3198955594291514E4 * t204 * t203 * t201 * t496 * t245 +
        0.2277620371046502E5 * t390 * t127 * t125 -
        0.3110666290793476E4 * t100 * t372 * t5 -
        0.4770905375469706E3 * t344 * t259 +
        0.3937610814400307E4 * t66 * t258 * t5 +
        0.2277620371046502E5 * t386 * t127 * t125 -
        0.7543463743812542E3 * t206 * t199 * t6 * t5 +
        0.2277620371046502E5 * t101 * t127 * t125 -
        0.3198955594291514E4 * t203 * t202 * t496 * t245 -
        0.1138810185523251E5 * t387 * t438 -
        0.4885436707797789E2 * t356 * t354 * t100 * t567 * t504 +
        0.905944924928429E4 * t191 * t113 * t7;
      double t580 = t123 * t171;
      double t587 = t353 * t100;
      double t594 = t115 * t112;
      double t601 = t42 * t39;
      double t605 = t270 * t31;
      double t612 = t27 * t6 * t5;
      double t624 = t412 * t171;
      double t625 = t183 * t43;
      double t626 = t625 * t414;
      double t629 =
        -0.4885436707797789E2 * t357 * t567 * t504 -
        0.433377270694286E4 * t340 * t580 +
        0.905944924928429E4 * t117 * t78 * t113 * t7 -
        0.4885436707797789E2 * t356 * t587 * t567 * t504 +
        0.408445310186476E4 * t185 * t189 +
        0.408445310186476E4 * t184 * t594 * t189 +
        0.408445310186476E4 * t184 * t190 * t189 -
        0.4136875594092619E4 * t601 * t177 * t37 +
        0.1670202100914981E5 * t605 * t468 +
        0.905944924928429E4 * t116 * t113 * t7 +
        0.1593555296195853E5 * t107 * t612 -
        0.4885436707797789E2 * t354 * t587 * t567 * t504 -
        0.4770905375469706E3 * t401 * t259 -
        0.2284393995262139E5 * t101 * t195 * t7 +
        0.4648474261967919E4 * t626 * t624;
      double t631 = t41 * t39;
      double t646 = t30 * t429;
      double t649 = t26 * t65;
      double t660 = t38 * t6 * t5;
      double t676 =
        -0.4136875594092619E4 * t631 * t177 * t37 -
        0.1670202100914981E5 * t605 * t612 + 0.1565236214672771E3 * t147 * t259 +
        0.1565236214672771E3 * t146 * t140 * t136 * t259 -
        0.6802839278405098E4 * t474 * t47 * t7 +
        0.6719369920237843E4 * t646 * t473 - 0.6719369920237843E4 * t646 * t649 -
        0.1670202100914981E5 * t605 * t97 -
        0.1670202100914981E5 * t270 * t98 * t97 +
        0.1509269882502665E4 * t103 * t447 - 0.4648474261967919E4 * t626 * t660 -
        0.4648474261967919E4 * t626 * t37 -
        0.2939953262497769E5 * t42 * t631 * t177 * t37 -
        0.7411664677842599E3 * t263 * t7 -
        0.3187110592391705E5 * t102 * t100 * t30 * t429 * t97;
      double t694 = t28 * t65;
      double t700 = t123 * t65;
      double t722 = t270 * t127;
      double t725 =
        -0.1670202100914981E5 * t270 * t429 * t97 +
        0.626491227419208E4 * t126 * t167 * t5 -
        0.6802839278405098E4 * t474 * t28 * t7 +
        0.3340404201829962E5 * t68 * t30 * t429 * t97 +
        0.4004764845831986E4 * t79 * t450 * t77 * t74 -
        0.1670202100914981E5 * t605 * t694 -
        0.4648474261967919E4 * t625 * t40 * t37 +
        0.433377270694286E4 * t340 * t700 -
        0.4648474261967919E4 * t625 * t177 * t37 -
        0.7411664677842599E3 * t377 * t7 +
        0.3340404201829962E5 * t66 * t30 * t429 * t97 +
        0.408445310186476E4 * t183 * t115 * t190 * t189 +
        0.3097755123859415E4 * t93 * t7 -
        0.4648474261967919E4 * t183 * t601 * t177 * t37 +
        0.2181531834986068E5 * t722 * t700;
      double t728 = t51 * t171;
      double t730 = t386 * t99 * t53;
      double t738 = t112 * t6 * t5;
      double t757 = t113 * t171;
      double t778 =
        -0.1142196997631069E5 * t730 * t728 - 0.6719369920237843E4 * t646 * t7 +
        0.1142196997631069E5 * t730 * t7 + 0.408445310186476E4 * t185 * t738 +
        0.1565236214672771E3 * t146 * t511 * t259 +
        0.1565236214672771E3 * t145 * t143 * t140 * t511 * t259 +
        0.1565236214672771E3 * t145 * t304 * t511 * t259 -
        0.6802839278405098E4 * t183 * t29 * t28 * t7 -
        0.2137096385075675E5 * t120 * t757 + 0.116502065962443E5 * t229 * t728 -
        0.4648474261967919E4 * t183 * t631 * t177 * t37 +
        0.1565236214672771E3 * t143 * t304 * t511 * t259 -
        0.116502065962443E5 * t229 * t7 +
        0.1142196997631069E5 * t386 * t99 * t51 * t7 +
        0.4212495413991277E4 * t93 * t285;
      double t783 = t99 * t30;
      double t827 =
        0.4004764845831986E4 * t223 * t78 * t75 * t74 -
        0.3187110592391705E5 * t100 * t783 * t429 * t97 -
        0.3187110592391705E5 * t102 * t783 * t429 * t97 -
        0.6802839278405098E4 * t475 * t649 -
        0.427419277015135E5 * t68 * t191 * t190 * t189 -
        0.116502065962443E5 * t68 * t66 * t51 * t7 +
        0.626491227419208E4 * t126 * t124 * t5 -
        0.2181531834986068E5 * t722 * t580 +
        0.4004764845831986E4 * t223 * t77 * t74 +
        0.233004131924886E5 * t66 * t195 * t7 + 0.2181531834986068E5 * t722 * t7 +
        0.2181531834986068E5 * t270 * t333 * t7 +
        0.2181531834986068E5 * t270 * t172 * t7 +
        0.4136875594092619E4 * t415 * t624 -
        0.427419277015135E5 * t66 * t191 * t190 * t189;
      double t860 = t136 * t10;
      double t865 = 0.4810955780432676E3 * t248 + 0.3660552413425374E3 * t171 +
        0.626491227419208E4 * t78 * t124 * t5 + 0.9795961824852758E3 * t93 * t248
        - 0.4363063669972135E5 * t68 * t126 * t172 * t7 -
        0.4363063669972135E5 * t66 * t126 * t172 * t7 +
        0.1593555296195853E5 * t107 * t694 + 0.1142196997631069E5 * t730 * t227 -
        0.3122267924171802E4 * t196 * t65 + 0.1555333145396738E4 * t103 * t65 +
        0.3705832338921299E3 * t206 * t199 * t65 -
        0.4244405766487319E2 * t306 * t511 * t65 -
        0.4136875594092619E4 * t415 * t660 +
        0.7826181073363854E2 * t147 * t860 * t171 - 0.3660552413425374E3 * t65;
      double t898 =
        0.2137096385075675E5 * t120 * t189 - 0.4136875594092619E4 * t415 * t37 -
        0.3195845749123036E4 * t433 - 0.733381126084165E3 * t7 -
        0.4810955780432676E3 * t447 + 0.356505526591888E3 * t183 * t171 +
        0.2137096385075675E5 * t119 * t594 * t189 +
        0.2137096385075675E5 * t120 * t738 - 0.4648474261967919E4 * t626 * t413 -
        0.408445310186476E4 * t185 * t757 -
        0.2284393995262139E5 * t386 * t195 * t7 +
        0.233004131924886E5 * t68 * t195 * t7 - 0.356505526591888E3 * t183 * t65 -
        0.7826181073363854E2 * t147 * t860 * t65 +
        0.2137096385075675E5 * t119 * t190 * t189;
      double t902 =
        t110 + t209 + t266 + t309 + t371 + t420 + t467 + t519 + t576 + t629 +
        t676 + t725 + t778 + t827 + t865 + t898;
      return t902;
    }

    static Shapeset::shape_fn_t ortho2_tri_fn[] =
    {
      ortho2_f1,    ortho2_f2,    ortho2_f3,    ortho2_f4,    ortho2_f5,    ortho2_f6,    ortho2_f7_0,
      ortho2_f7_1,  ortho2_f8_0,  ortho2_f8_1,  ortho2_f9_0,  ortho2_f9_1,  ortho2_f10,   ortho2_f11,
      ortho2_f12,   ortho2_f13,   ortho2_f14,   ortho2_f15,   ortho2_f16_0, ortho2_f16_1, ortho2_f17_0,
      ortho2_f17_1, ortho2_f18_0, ortho2_f18_1, ortho2_f19,   ortho2_f20,   ortho2_f21,   ortho2_f22,
      ortho2_f23,   ortho2_f24,   ortho2_f25,   ortho2_f26,   ortho2_f27,   ortho2_f28,   ortho2_f29_0,
      ortho2_f29_1, ortho2_f30_0, ortho2_f30_1, ortho2_f31_0, ortho2_f31_1, ortho2_f32,   ortho2_f33,
      ortho2_f34,   ortho2_f35,   ortho2_f36,   ortho2_f37,   ortho2_f38,   ortho2_f39,   ortho2_f40,
      ortho2_f41,   ortho2_f42,   ortho2_f43,   ortho2_f44,   ortho2_f45,   ortho2_f46_0, ortho2_f46_1,
      ortho2_f47_0, ortho2_f47_1, ortho2_f48_0, ortho2_f48_1, ortho2_f49,   ortho2_f50,   ortho2_f51,
      ortho2_f52,   ortho2_f53,   ortho2_f54,   ortho2_f55,   ortho2_f56,   ortho2_f57,   ortho2_f58,
      ortho2_f59,   ortho2_f60,   ortho2_f61,   ortho2_f62,   ortho2_f63,   ortho2_f64,   ortho2_f65,
      ortho2_f66
    };

    static Shapeset::shape_fn_t ortho2_tri_fn_dx[] =
    {
     ortho2_f1x,    ortho2_f2x,    ortho2_f3x,    ortho2_f4x,    ortho2_f5x,    ortho2_f6x,    ortho2_f7x_0,
      ortho2_f7x_1,  ortho2_f8x_0,  ortho2_f8x_1,  ortho2_f9x_0,  ortho2_f9x_1,  ortho2_f10x,   ortho2_f11x,
      ortho2_f12x,   ortho2_f13x,   ortho2_f14x,   ortho2_f15x,   ortho2_f16x_0, ortho2_f16x_1, ortho2_f17x_0,
      ortho2_f17x_1, ortho2_f18x_0, ortho2_f18x_1, ortho2_f19x,   ortho2_f20x,   ortho2_f21x,   ortho2_f22x,
      ortho2_f23x,   ortho2_f24x,   ortho2_f25x,   ortho2_f26x,   ortho2_f27x,   ortho2_f28x,   ortho2_f29x_0,
      ortho2_f29x_1, ortho2_f30x_0, ortho2_f30x_1, ortho2_f31x_0, ortho2_f31x_1, ortho2_f32x,   ortho2_f33x,
      ortho2_f34x,   ortho2_f35x,   ortho2_f36x,   ortho2_f37x,   ortho2_f38x,   ortho2_f39x,   ortho2_f40x,
      ortho2_f41x,   ortho2_f42x,   ortho2_f43x,   ortho2_f44x,   ortho2_f45x,   ortho2_f46x_0, ortho2_f46x_1,
      ortho2_f47x_0, ortho2_f47x_1, ortho2_f48x_0, ortho2_f48x_1, ortho2_f49x,   ortho2_f50x,   ortho2_f51x,
      ortho2_f52x,   ortho2_f53x,   ortho2_f54x,   ortho2_f55x,   ortho2_f56x,   ortho2_f57x,   ortho2_f58x,
      ortho2_f59x,   ortho2_f60x,   ortho2_f61x,   ortho2_f62x,   ortho2_f63x,   ortho2_f64x,   ortho2_f65x,
      ortho2_f66x
    };

    static Shapeset::shape_fn_t ortho2_tri_fn_dy[] =
    {
     ortho2_f1y,    ortho2_f2y,    ortho2_f3y,    ortho2_f4y,    ortho2_f5y,    ortho2_f6y,    ortho2_f7y_0,
      ortho2_f7y_1,  ortho2_f8y_0,  ortho2_f8y_1,  ortho2_f9y_0,  ortho2_f9y_1,  ortho2_f10y,   ortho2_f11y,
      ortho2_f12y,   ortho2_f13y,   ortho2_f14y,   ortho2_f15y,   ortho2_f16y_0, ortho2_f16y_1, ortho2_f17y_0,
      ortho2_f17y_1, ortho2_f18y_0, ortho2_f18y_1, ortho2_f19y,   ortho2_f20y,   ortho2_f21y,   ortho2_f22y,
      ortho2_f23y,   ortho2_f24y,   ortho2_f25y,   ortho2_f26y,   ortho2_f27y,   ortho2_f28y,   ortho2_f29y_0,
      ortho2_f29y_1, ortho2_f30y_0, ortho2_f30y_1, ortho2_f31y_0, ortho2_f31y_1, ortho2_f32y,   ortho2_f33y,
      ortho2_f34y,   ortho2_f35y,   ortho2_f36y,   ortho2_f37y,   ortho2_f38y,   ortho2_f39y,   ortho2_f40y,
      ortho2_f41y,   ortho2_f42y,   ortho2_f43y,   ortho2_f44y,   ortho2_f45y,   ortho2_f46y_0, ortho2_f46y_1,
      ortho2_f47y_0, ortho2_f47y_1, ortho2_f48y_0, ortho2_f48y_1, ortho2_f49y,   ortho2_f50y,   ortho2_f51y,
      ortho2_f52y,   ortho2_f53y,   ortho2_f54y,   ortho2_f55y,   ortho2_f56y,   ortho2_f57y,   ortho2_f58y,
      ortho2_f59y,   ortho2_f60y,   ortho2_f61y,   ortho2_f62y,   ortho2_f63y,   ortho2_f64y,   ortho2_f65y,
      ortho2_f66y
    };

    static int ortho2_tri_bubble_indices_all_orders[] =
    {
      12,
      16, 17,
      24, 25, 26,
      30, 31, 32, 33,
      40, 41, 42, 43, 44,
      48, 49, 50, 51, 52, 53,
      60, 61, 62, 63, 64, 65, 66,
      70, 71, 72, 73, 74, 75, 76, 77
    };

    static int* ortho2_tri_bubble_indices[11] =
    {
      NULL, NULL, NULL,
      ortho2_tri_bubble_indices_all_orders,
      ortho2_tri_bubble_indices_all_orders,
      ortho2_tri_bubble_indices_all_orders,
      ortho2_tri_bubble_indices_all_orders,
      ortho2_tri_bubble_indices_all_orders,
      ortho2_tri_bubble_indices_all_orders,
      ortho2_tri_bubble_indices_all_orders,
      ortho2_tri_bubble_indices_all_orders
    };

    static int ortho2_tri_bubble_count[11] = { 0, 0, 0, 1, 3, 6, 10, 15, 21, 28, 36 };

    static int ortho2_tri_edge_indices_0[22] =  { 0, 1, 1, 0, 3, 3, 6,  7,  13, 13, 18, 19, 27, 27, 34, 35, 45, 45, 54, 55, 67, 67 };
    static int ortho2_tri_edge_indices_1[22] =  { 1, 2, 2, 1, 4, 4, 8,  9,  14, 14, 20, 21, 28, 28, 36, 37, 46, 46, 56, 57, 68, 68 };
    static int ortho2_tri_edge_indices_2[22] =  { 2, 0, 0, 2, 5, 5, 10, 11, 15, 15, 22, 23, 29, 29, 38, 39, 47, 47, 58, 59, 69, 69 };

    static int* ortho2_tri_edge_indices[3] =
    {
      ortho2_tri_edge_indices_0,
      ortho2_tri_edge_indices_1,
      ortho2_tri_edge_indices_2
    };

    static int ortho2_tri_vertex_indices[3] = { 0, 1, 2 };

    static int ortho2_tri_index_to_order[78] =
    {
      1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6,
      7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
      10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
    };

    static Shapeset::shape_fn_t* ortho2_tri_shape_fn_table[1] =
    {
      ortho2_tri_fn
    };

    static Shapeset::shape_fn_t* ortho2_tri_shape_fn_table_dx[1] =
    {
      ortho2_tri_fn_dx
    };

    static Shapeset::shape_fn_t* ortho2_tri_shape_fn_table_dy[1] =
    {
      ortho2_tri_fn_dy
    };

    //// triangle and quad tables and class constructor ///////////////////////////////////////////////

    #include "shapeset_h1_quad.h"

    static Shapeset::shape_fn_t** ortho2_shape_fn_table[2] =
    {
      ortho2_tri_shape_fn_table,
      simple_quad_shape_fn_table
    };

    static Shapeset::shape_fn_t** ortho2_shape_fn_table_dx[2] =
    {
      ortho2_tri_shape_fn_table_dx,
      simple_quad_shape_fn_table_dx
    };

    static Shapeset::shape_fn_t** ortho2_shape_fn_table_dy[2] =
    {
      ortho2_tri_shape_fn_table_dy,
      simple_quad_shape_fn_table_dy
    };

    static int* ortho2_vertex_indices[2] =
    {
      ortho2_tri_vertex_indices,
      simple_quad_vertex_indices
    };

    static int** ortho2_edge_indices[2] =
    {
      ortho2_tri_edge_indices,
      simple_quad_edge_indices
    };

    static int** ortho2_bubble_indices[2] =
    {
      ortho2_tri_bubble_indices,
      simple_quad_bubble_indices
    };

    static int* ortho2_bubble_count[2] =
    {
      ortho2_tri_bubble_count,
      simple_quad_bubble_count
    };

    static int* ortho2_index_to_order[2] =
    {
      ortho2_tri_index_to_order,
      simple_quad_index_to_order
    };

    int H1ShapesetOrtho::get_max_index(ElementMode2D mode) { return max_index[mode]; }

    H1ShapesetOrtho::H1ShapesetOrtho()
    {
      shape_table[0] = ortho2_shape_fn_table;
      shape_table[1] = ortho2_shape_fn_table_dx;
      shape_table[2] = ortho2_shape_fn_table_dy;
      shape_table[3] = NULL;
      shape_table[4] = NULL;
      shape_table[5] = NULL;

      vertex_indices = ortho2_vertex_indices;
      edge_indices = ortho2_edge_indices;
      bubble_indices = ortho2_bubble_indices;
      bubble_count = ortho2_bubble_count;
      index_to_order = ortho2_index_to_order;

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

      ebias = 2;

      comb_table = NULL;
    }

    const int H1ShapesetOrtho::max_index[2] = { 77, 136 };
  }
}