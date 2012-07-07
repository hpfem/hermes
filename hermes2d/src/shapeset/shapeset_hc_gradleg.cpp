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
#include "shapeset_hc_all.h"
namespace Hermes
{
  namespace Hermes2D
  {
    // Shape functions for the curl operator for triangles,

    //Shape functions

    ///////////////////////////////// ORDER 0 //////////////////////////////////

    /* Whitney fns - constant tangential component */

    // EDGE 1

    static double gradleg_tri_p0_e1_a_0(double x, double y)
    {
      return psi0e1_1(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e1_a_1(double x, double y)
    {
      return -(gradleg_tri_p0_e1_a_0(x, y));
    }

    static double gradleg_tri_p0_e1_b_0(double x, double y)
    {
      return psi0e1_2(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e1_b_1(double x, double y)
    {
      return -(gradleg_tri_p0_e1_b_0(x, y));
    }

    static double gradleg_tri_p0_e1_ax_0(double x, double y)
    {
      return psi0e1x_1(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e1_ax_1(double x, double y)
    {
      return -(gradleg_tri_p0_e1_ax_0(x, y));
    }

    static double gradleg_tri_p0_e1_bx_0(double x, double y)
    {
      return psi0e1x_2(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e1_bx_1(double x, double y)
    {
      return -(gradleg_tri_p0_e1_bx_0(x, y));
    }

    static double gradleg_tri_p0_e1_ay_0(double x, double y)
    {
      return psi0e1y_1(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e1_ay_1(double x, double y)
    {
      return -(gradleg_tri_p0_e1_ay_0(x, y));
    }

    static double gradleg_tri_p0_e1_by_0(double x, double y)
    {
      return psi0e1y_2(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e1_by_1(double x, double y)
    {
      return -(gradleg_tri_p0_e1_by_0(x, y));
    }

    // EDGE 2

    static double gradleg_tri_p0_e2_a_0(double x, double y)
    {
      return psi0e2_1(x, y) / 1.4142135623731;
    }

    static double gradleg_tri_p0_e2_a_1(double x, double y)
    {
      return -(gradleg_tri_p0_e2_a_0(x, y));
    }

    static double gradleg_tri_p0_e2_b_0(double x, double y)
    {
      return psi0e2_2(x, y) / 1.4142135623731;
    }

    static double gradleg_tri_p0_e2_b_1(double x, double y)
    {
      return -(gradleg_tri_p0_e2_b_0(x, y));
    }

    static double gradleg_tri_p0_e2_ax_0(double x, double y)
    {
      return psi0e2x_1(x, y) / 1.4142135623731;
    }

    static double gradleg_tri_p0_e2_ax_1(double x, double y)
    {
      return -(gradleg_tri_p0_e2_ax_0(x, y));
    }

    static double gradleg_tri_p0_e2_bx_0(double x, double y)
    {
      return psi0e2x_2(x, y) / 1.4142135623731;
    }

    static double gradleg_tri_p0_e2_bx_1(double x, double y)
    {
      return -(gradleg_tri_p0_e2_bx_0(x, y));
    }

    static double gradleg_tri_p0_e2_ay_0(double x, double y)
    {
      return psi0e2y_1(x, y) / 1.4142135623731;
    }

    static double gradleg_tri_p0_e2_ay_1(double x, double y)
    {
      return -(gradleg_tri_p0_e2_ay_0(x, y));
    }

    static double gradleg_tri_p0_e2_by_0(double x, double y)
    {
      return psi0e2y_2(x, y) / 1.4142135623731;
    }

    static double gradleg_tri_p0_e2_by_1(double x, double y)
    {
      return -(gradleg_tri_p0_e2_by_0(x, y));
    }

    // EDGE 3

    static double gradleg_tri_p0_e3_a_0(double x, double y)
    {
      return psi0e3_1(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e3_a_1(double x, double y)
    {
      return -(gradleg_tri_p0_e3_a_0(x, y));
    }

    static double gradleg_tri_p0_e3_b_0(double x, double y)
    {
      return psi0e3_2(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e3_b_1(double x, double y)
    {
      return -(gradleg_tri_p0_e3_b_0(x, y));
    }

    static double gradleg_tri_p0_e3_ax_0(double x, double y)
    {
      return psi0e3x_1(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e3_ax_1(double x, double y)
    {
      return -(gradleg_tri_p0_e3_ax_0(x, y));
    }

    static double gradleg_tri_p0_e3_bx_0(double x, double y)
    {
      return psi0e3x_2(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e3_bx_1(double x, double y)
    {
      return -(gradleg_tri_p0_e3_bx_0(x, y));
    }

    static double gradleg_tri_p0_e3_ay_0(double x, double y)
    {
      return psi0e3y_1(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e3_ay_1(double x, double y)
    {
      return -(gradleg_tri_p0_e3_ay_0(x, y));
    }

    static double gradleg_tri_p0_e3_by_0(double x, double y)
    {
      return psi0e3y_2(x, y) / 1.0000000000000;
    }

    static double gradleg_tri_p0_e3_by_1(double x, double y)
    {
      return -(gradleg_tri_p0_e3_by_0(x, y));
    }

    ///////////////////////////////// ORDER 1 //////////////////////////////////

    /* EDGE FUNCTIONS - order 1*/

     /* EDGE 1 */

    static double gradleg_tri_p1_e1_a(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return  (l3x * l2 * phi0(l3 - l2) + l3 * l2x * phi0(l3 - l2) + l3 * l2 * phi0x(l3 - l2) * (l3x - l2x)) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e1_b(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return  (l3y * l2 * phi0(l3 - l2) + l3 * l2y * phi0(l3 - l2) + l3 * l2 * phi0x(l3 - l2) * (l3y - l2y)) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     double ker, kerx, kerxx;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     ker = phi0(l3 - l2); kerx = phi0x(l3 - l2) * (l3x - l2x); kerxx = phi0xx(l3 - l2) * sqr(l3x - l2x);
     return  (2.0 * l3x * l2x * ker + 2.0 * l3x * l2 * kerx + 2.0 * l3 * l2x * kerx + l3 * l2 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e1_ay(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi0(l3 - l2); kerx = phi0x(l3 - l2) * (l3x - l2x);
     kery = phi0x(l3 - l2) * (l3y - l2y); kerxy = phi0xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e1_bx(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi0(l3 - l2); kerx = phi0x(l3 - l2) * (l3x - l2x);
     kery = phi0x(l3 - l2) * (l3y - l2y); kerxy = phi0xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     double ker, kery, keryy;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     ker = phi0(l3 - l2); kery = phi0x(l3 - l2) * (l3y - l2y); keryy = phi0xx(l3 - l2) * sqr(l3y - l2y);
     return  (2.0 * l3y * l2y * ker + 2.0 * l3y * l2 * kery + 2.0 * l3 * l2y * kery + l3 * l2 * keryy) / 1.0000000000000;
    }

     /* EDGE 2 */

    static double gradleg_tri_p1_e2_a(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return  (l1x * l3 * phi0(l1 - l3) + l1 * l3x * phi0(l1 - l3) + l1 * l3 * phi0x(l1 - l3) * (l1x - l3x)) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e2_b(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return  (l1y * l3 * phi0(l1 - l3) + l1 * l3y * phi0(l1 - l3) + l1 * l3 * phi0x(l1 - l3) * (l1y - l3y)) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     double ker, kerx, kerxx;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     ker = phi0(l1 - l3); kerx = phi0x(l1 - l3) * (l1x - l3x); kerxx = phi0xx(l1 - l3) * sqr(l1x - l3x);
     return  (2.0 * l1x * l3x * ker + 2.0 * l1x * l3 * kerx + 2.0 * l1 * l3x * kerx + l1 * l3 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e2_ay(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi0(l1 - l3); kerx = phi0x(l1 - l3) * (l1x - l3x);
     kery = phi0x(l1 - l3) * (l1y - l3y); kerxy = phi0xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e2_bx(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi0(l1 - l3); kerx = phi0x(l1 - l3) * (l1x - l3x);
     kery = phi0x(l1 - l3) * (l1y - l3y); kerxy = phi0xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     double ker, kery, keryy;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     ker = phi0(l1 - l3); kery = phi0x(l1 - l3) * (l1y - l3y); keryy = phi0xx(l1 - l3) * sqr(l1y - l3y);
     return  (2.0 * l1y * l3y * ker + 2.0 * l1y * l3 * kery + 2.0 * l1 * l3y * kery + l1 * l3 * keryy) / 1.0000000000000;
    }

     /* EDGE 3 */

    static double gradleg_tri_p1_e3_a(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return  (l2x * l1 * phi0(l2 - l1) + l2 * l1x * phi0(l2 - l1) + l2 * l1 * phi0x(l2 - l1) * (l2x - l1x)) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e3_b(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return  (l2y * l1 * phi0(l2 - l1) + l2 * l1y * phi0(l2 - l1) + l2 * l1 * phi0x(l2 - l1) * (l2y - l1y)) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     double ker, kerx, kerxx;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     ker = phi0(l2 - l1); kerx = phi0x(l2 - l1) * (l2x - l1x); kerxx = phi0xx(l2 - l1) * sqr(l2x - l1x);
     return  (2.0 * l2x * l1x * ker + 2.0 * l2x * l1 * kerx + 2.0 * l2 * l1x * kerx + l2 * l1 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e3_ay(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi0(l2 - l1); kerx = phi0x(l2 - l1) * (l2x - l1x);
     kery = phi0x(l2 - l1) * (l2y - l1y); kerxy = phi0xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e3_bx(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi0(l2 - l1); kerx = phi0x(l2 - l1) * (l2x - l1x);
     kery = phi0x(l2 - l1) * (l2y - l1y); kerxy = phi0xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p1_e3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     double ker, kery, keryy;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     ker = phi0(l2 - l1); kery = phi0x(l2 - l1) * (l2y - l1y); keryy = phi0xx(l2 - l1) * sqr(l2y - l1y);
     return  (2.0 * l2y * l1y * ker + 2.0 * l2y * l1 * kery + 2.0 * l2 * l1y * kery + l2 * l1 * keryy) / 1.0000000000000;
    }

    ///////////////////////////////// ORDER 2 //////////////////////////////////

    /* EDGE FUNCTIONS - order 2*/

     /* EDGE 1 */

    static double gradleg_tri_p2_e1_a_0(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return  (l3x * l2 * phi1(l3 - l2) + l3 * l2x * phi1(l3 - l2) + l3 * l2 * phi1x(l3 - l2) * (l3x - l2x)) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e1_a_1(double x, double y)
    {
     return -(gradleg_tri_p2_e1_a_0(x, y));
    }

    static double gradleg_tri_p2_e1_b_0(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return  (l3y * l2 * phi1(l3 - l2) + l3 * l2y * phi1(l3 - l2) + l3 * l2 * phi1x(l3 - l2) * (l3y - l2y)) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e1_b_1(double x, double y)
    {
     return -(gradleg_tri_p2_e1_b_0(x, y));
    }

    static double gradleg_tri_p2_e1_ax_0(double x, double y)
    {
     double l3, l3x, l2, l2x;
     double ker, kerx, kerxx;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     ker = phi1(l3 - l2); kerx = phi1x(l3 - l2) * (l3x - l2x); kerxx = phi1xx(l3 - l2) * sqr(l3x - l2x);
     return  (2.0 * l3x * l2x * ker + 2.0 * l3x * l2 * kerx + 2.0 * l3 * l2x * kerx + l3 * l2 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e1_ax_1(double x, double y)
    {
     return -(gradleg_tri_p2_e1_ax_0(x, y));
    }

    static double gradleg_tri_p2_e1_ay_0(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi1(l3 - l2); kerx = phi1x(l3 - l2) * (l3x - l2x);
     kery = phi1x(l3 - l2) * (l3y - l2y); kerxy = phi1xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e1_ay_1(double x, double y)
    {
     return -(gradleg_tri_p2_e1_ay_0(x, y));
    }

    static double gradleg_tri_p2_e1_bx_0(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi1(l3 - l2); kerx = phi1x(l3 - l2) * (l3x - l2x);
     kery = phi1x(l3 - l2) * (l3y - l2y); kerxy = phi1xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e1_bx_1(double x, double y)
    {
     return -(gradleg_tri_p2_e1_bx_0(x, y));
    }

    static double gradleg_tri_p2_e1_by_0(double x, double y)
    {
     double l3, l3y, l2, l2y;
     double ker, kery, keryy;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     ker = phi1(l3 - l2); kery = phi1x(l3 - l2) * (l3y - l2y); keryy = phi1xx(l3 - l2) * sqr(l3y - l2y);
     return  (2.0 * l3y * l2y * ker + 2.0 * l3y * l2 * kery + 2.0 * l3 * l2y * kery + l3 * l2 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e1_by_1(double x, double y)
    {
     return -(gradleg_tri_p2_e1_by_0(x, y));
    }

     /* EDGE 2 */

    static double gradleg_tri_p2_e2_a_0(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return  (l1x * l3 * phi1(l1 - l3) + l1 * l3x * phi1(l1 - l3) + l1 * l3 * phi1x(l1 - l3) * (l1x - l3x)) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e2_a_1(double x, double y)
    {
     return -(gradleg_tri_p2_e2_a_0(x, y));
    }

    static double gradleg_tri_p2_e2_b_0(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return  (l1y * l3 * phi1(l1 - l3) + l1 * l3y * phi1(l1 - l3) + l1 * l3 * phi1x(l1 - l3) * (l1y - l3y)) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e2_b_1(double x, double y)
    {
     return -(gradleg_tri_p2_e2_b_0(x, y));
    }

    static double gradleg_tri_p2_e2_ax_0(double x, double y)
    {
     double l1, l1x, l3, l3x;
     double ker, kerx, kerxx;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     ker = phi1(l1 - l3); kerx = phi1x(l1 - l3) * (l1x - l3x); kerxx = phi1xx(l1 - l3) * sqr(l1x - l3x);
     return  (2.0 * l1x * l3x * ker + 2.0 * l1x * l3 * kerx + 2.0 * l1 * l3x * kerx + l1 * l3 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e2_ax_1(double x, double y)
    {
     return -(gradleg_tri_p2_e2_ax_0(x, y));
    }

    static double gradleg_tri_p2_e2_ay_0(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi1(l1 - l3); kerx = phi1x(l1 - l3) * (l1x - l3x);
     kery = phi1x(l1 - l3) * (l1y - l3y); kerxy = phi1xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e2_ay_1(double x, double y)
    {
     return -(gradleg_tri_p2_e2_ay_0(x, y));
    }

    static double gradleg_tri_p2_e2_bx_0(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi1(l1 - l3); kerx = phi1x(l1 - l3) * (l1x - l3x);
     kery = phi1x(l1 - l3) * (l1y - l3y); kerxy = phi1xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e2_bx_1(double x, double y)
    {
     return -(gradleg_tri_p2_e2_bx_0(x, y));
    }

    static double gradleg_tri_p2_e2_by_0(double x, double y)
    {
     double l1, l1y, l3, l3y;
     double ker, kery, keryy;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     ker = phi1(l1 - l3); kery = phi1x(l1 - l3) * (l1y - l3y); keryy = phi1xx(l1 - l3) * sqr(l1y - l3y);
     return  (2.0 * l1y * l3y * ker + 2.0 * l1y * l3 * kery + 2.0 * l1 * l3y * kery + l1 * l3 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e2_by_1(double x, double y)
    {
     return -(gradleg_tri_p2_e2_by_0(x, y));
    }

     /* EDGE 3 */

    static double gradleg_tri_p2_e3_a_0(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return  (l2x * l1 * phi1(l2 - l1) + l2 * l1x * phi1(l2 - l1) + l2 * l1 * phi1x(l2 - l1) * (l2x - l1x)) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e3_a_1(double x, double y)
    {
     return -(gradleg_tri_p2_e3_a_0(x, y));
    }

    static double gradleg_tri_p2_e3_b_0(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return  (l2y * l1 * phi1(l2 - l1) + l2 * l1y * phi1(l2 - l1) + l2 * l1 * phi1x(l2 - l1) * (l2y - l1y)) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e3_b_1(double x, double y)
    {
     return -(gradleg_tri_p2_e3_b_0(x, y));
    }

    static double gradleg_tri_p2_e3_ax_0(double x, double y)
    {
     double l2, l2x, l1, l1x;
     double ker, kerx, kerxx;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     ker = phi1(l2 - l1); kerx = phi1x(l2 - l1) * (l2x - l1x); kerxx = phi1xx(l2 - l1) * sqr(l2x - l1x);
     return  (2.0 * l2x * l1x * ker + 2.0 * l2x * l1 * kerx + 2.0 * l2 * l1x * kerx + l2 * l1 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e3_ax_1(double x, double y)
    {
     return -(gradleg_tri_p2_e3_ax_0(x, y));
    }

    static double gradleg_tri_p2_e3_ay_0(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi1(l2 - l1); kerx = phi1x(l2 - l1) * (l2x - l1x);
     kery = phi1x(l2 - l1) * (l2y - l1y); kerxy = phi1xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e3_ay_1(double x, double y)
    {
     return -(gradleg_tri_p2_e3_ay_0(x, y));
    }

    static double gradleg_tri_p2_e3_bx_0(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi1(l2 - l1); kerx = phi1x(l2 - l1) * (l2x - l1x);
     kery = phi1x(l2 - l1) * (l2y - l1y); kerxy = phi1xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e3_bx_1(double x, double y)
    {
     return -(gradleg_tri_p2_e3_bx_0(x, y));
    }

    static double gradleg_tri_p2_e3_by_0(double x, double y)
    {
     double l2, l2y, l1, l1y;
     double ker, kery, keryy;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     ker = phi1(l2 - l1); kery = phi1x(l2 - l1) * (l2y - l1y); keryy = phi1xx(l2 - l1) * sqr(l2y - l1y);
     return  (2.0 * l2y * l1y * ker + 2.0 * l2y * l1 * kery + 2.0 * l2 * l1y * kery + l2 * l1 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p2_e3_by_1(double x, double y)
    {
     return -(gradleg_tri_p2_e3_by_0(x, y));
    }

    /* BUBBLE */

    /* Edge-based BUBBLE - order 2 */

     // EDGE 1
    static double gradleg_tri_p2_b1_a(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n11 * (l3 * l2 * Legendre0(l3 - l2));
    }

    static double gradleg_tri_p2_b1_b(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n12 * (l3 * l2 * Legendre0(l3 - l2));
    }

    static double gradleg_tri_p2_b1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n11 * (l3x * l2 * Legendre0(l3 - l2) + l3 * l2x * Legendre0(l3 - l2) + l3 * l2 * Legendre0x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p2_b1_bx(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n12 * (l3x * l2 * Legendre0(l3 - l2) + l3 * l2x * Legendre0(l3 - l2) + l3 * l2 * Legendre0x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p2_b1_ay(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n11 * (l3y * l2 * Legendre0(l3 - l2) + l3 * l2y * Legendre0(l3 - l2) + l3 * l2 * Legendre0x(l3 - l2) * (l3y - l2y));
    }

    static double gradleg_tri_p2_b1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n12 * (l3y * l2 * Legendre0(l3 - l2) + l3 * l2y * Legendre0(l3 - l2) + l3 * l2 * Legendre0x(l3 - l2) * (l3y - l2y));
    }

     // EDGE 2
    static double gradleg_tri_p2_b2_a(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n21 * (l1 * l3 * Legendre0(l1 - l3));
    }

    static double gradleg_tri_p2_b2_b(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n22 * (l1 * l3 * Legendre0(l1 - l3));
    }

    static double gradleg_tri_p2_b2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n21 * (l1x * l3 * Legendre0(l1 - l3) + l1 * l3x * Legendre0(l1 - l3) + l1 * l3 * Legendre0x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p2_b2_bx(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n22 * (l1x * l3 * Legendre0(l1 - l3) + l1 * l3x * Legendre0(l1 - l3) + l1 * l3 * Legendre0x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p2_b2_ay(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n21 * (l1y * l3 * Legendre0(l1 - l3) + l1 * l3y * Legendre0(l1 - l3) + l1 * l3 * Legendre0x(l1 - l3) * (l1y - l3y));
    }

    static double gradleg_tri_p2_b2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n22 * (l1y * l3 * Legendre0(l1 - l3) + l1 * l3y * Legendre0(l1 - l3) + l1 * l3 * Legendre0x(l1 - l3) * (l1y - l3y));
    }

     // EDGE 3
    static double gradleg_tri_p2_b3_a(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n31 * (l2 * l1 * Legendre0(l2 - l1));
    }

    static double gradleg_tri_p2_b3_b(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n32 * (l2 * l1 * Legendre0(l2 - l1));
    }

    static double gradleg_tri_p2_b3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n31 * (l2x * l1 * Legendre0(l2 - l1) + l2 * l1x * Legendre0(l2 - l1) + l2 * l1 * Legendre0x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p2_b3_bx(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n32 * (l2x * l1 * Legendre0(l2 - l1) + l2 * l1x * Legendre0(l2 - l1) + l2 * l1 * Legendre0x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p2_b3_ay(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n31 * (l2y * l1 * Legendre0(l2 - l1) + l2 * l1y * Legendre0(l2 - l1) + l2 * l1 * Legendre0x(l2 - l1) * (l2y - l1y));
    }

    static double gradleg_tri_p2_b3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n32 * (l2y * l1 * Legendre0(l2 - l1) + l2 * l1y * Legendre0(l2 - l1) + l2 * l1 * Legendre0x(l2 - l1) * (l2y - l1y));
    }

    ///////////////////////////////// ORDER 3 //////////////////////////////////

    /* EDGE FUNCTIONS - order 3*/

     /* EDGE 1 */

    static double gradleg_tri_p3_e1_a(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return  (l3x * l2 * phi2(l3 - l2) + l3 * l2x * phi2(l3 - l2) + l3 * l2 * phi2x(l3 - l2) * (l3x - l2x)) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e1_b(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return  (l3y * l2 * phi2(l3 - l2) + l3 * l2y * phi2(l3 - l2) + l3 * l2 * phi2x(l3 - l2) * (l3y - l2y)) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     double ker, kerx, kerxx;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     ker = phi2(l3 - l2); kerx = phi2x(l3 - l2) * (l3x - l2x); kerxx = phi2xx(l3 - l2) * sqr(l3x - l2x);
     return  (2.0 * l3x * l2x * ker + 2.0 * l3x * l2 * kerx + 2.0 * l3 * l2x * kerx + l3 * l2 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e1_ay(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi2(l3 - l2); kerx = phi2x(l3 - l2) * (l3x - l2x);
     kery = phi2x(l3 - l2) * (l3y - l2y); kerxy = phi2xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e1_bx(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi2(l3 - l2); kerx = phi2x(l3 - l2) * (l3x - l2x);
     kery = phi2x(l3 - l2) * (l3y - l2y); kerxy = phi2xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     double ker, kery, keryy;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     ker = phi2(l3 - l2); kery = phi2x(l3 - l2) * (l3y - l2y); keryy = phi2xx(l3 - l2) * sqr(l3y - l2y);
     return  (2.0 * l3y * l2y * ker + 2.0 * l3y * l2 * kery + 2.0 * l3 * l2y * kery + l3 * l2 * keryy) / 1.0000000000000;
    }

     /* EDGE 2 */

    static double gradleg_tri_p3_e2_a(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return  (l1x * l3 * phi2(l1 - l3) + l1 * l3x * phi2(l1 - l3) + l1 * l3 * phi2x(l1 - l3) * (l1x - l3x)) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e2_b(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return  (l1y * l3 * phi2(l1 - l3) + l1 * l3y * phi2(l1 - l3) + l1 * l3 * phi2x(l1 - l3) * (l1y - l3y)) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     double ker, kerx, kerxx;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     ker = phi2(l1 - l3); kerx = phi2x(l1 - l3) * (l1x - l3x); kerxx = phi2xx(l1 - l3) * sqr(l1x - l3x);
     return  (2.0 * l1x * l3x * ker + 2.0 * l1x * l3 * kerx + 2.0 * l1 * l3x * kerx + l1 * l3 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e2_ay(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi2(l1 - l3); kerx = phi2x(l1 - l3) * (l1x - l3x);
     kery = phi2x(l1 - l3) * (l1y - l3y); kerxy = phi2xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e2_bx(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi2(l1 - l3); kerx = phi2x(l1 - l3) * (l1x - l3x);
     kery = phi2x(l1 - l3) * (l1y - l3y); kerxy = phi2xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     double ker, kery, keryy;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     ker = phi2(l1 - l3); kery = phi2x(l1 - l3) * (l1y - l3y); keryy = phi2xx(l1 - l3) * sqr(l1y - l3y);
     return  (2.0 * l1y * l3y * ker + 2.0 * l1y * l3 * kery + 2.0 * l1 * l3y * kery + l1 * l3 * keryy) / 1.0000000000000;
    }

     /* EDGE 3 */

    static double gradleg_tri_p3_e3_a(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return  (l2x * l1 * phi2(l2 - l1) + l2 * l1x * phi2(l2 - l1) + l2 * l1 * phi2x(l2 - l1) * (l2x - l1x)) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e3_b(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return  (l2y * l1 * phi2(l2 - l1) + l2 * l1y * phi2(l2 - l1) + l2 * l1 * phi2x(l2 - l1) * (l2y - l1y)) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     double ker, kerx, kerxx;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     ker = phi2(l2 - l1); kerx = phi2x(l2 - l1) * (l2x - l1x); kerxx = phi2xx(l2 - l1) * sqr(l2x - l1x);
     return  (2.0 * l2x * l1x * ker + 2.0 * l2x * l1 * kerx + 2.0 * l2 * l1x * kerx + l2 * l1 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e3_ay(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi2(l2 - l1); kerx = phi2x(l2 - l1) * (l2x - l1x);
     kery = phi2x(l2 - l1) * (l2y - l1y); kerxy = phi2xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e3_bx(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi2(l2 - l1); kerx = phi2x(l2 - l1) * (l2x - l1x);
     kery = phi2x(l2 - l1) * (l2y - l1y); kerxy = phi2xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p3_e3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     double ker, kery, keryy;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     ker = phi2(l2 - l1); kery = phi2x(l2 - l1) * (l2y - l1y); keryy = phi2xx(l2 - l1) * sqr(l2y - l1y);
     return  (2.0 * l2y * l1y * ker + 2.0 * l2y * l1 * kery + 2.0 * l2 * l1y * kery + l2 * l1 * keryy) / 1.0000000000000;
    }

    /* BUBBLE */

    /* Edge-based BUBBLE - order 3 */

     // EDGE 1
    static double gradleg_tri_p3_b1_a(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n11 * (l3 * l2 * Legendre1(l3 - l2));
    }

    static double gradleg_tri_p3_b1_b(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n12 * (l3 * l2 * Legendre1(l3 - l2));
    }

    static double gradleg_tri_p3_b1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n11 * (l3x * l2 * Legendre1(l3 - l2) + l3 * l2x * Legendre1(l3 - l2) + l3 * l2 * Legendre1x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p3_b1_bx(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n12 * (l3x * l2 * Legendre1(l3 - l2) + l3 * l2x * Legendre1(l3 - l2) + l3 * l2 * Legendre1x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p3_b1_ay(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n11 * (l3y * l2 * Legendre1(l3 - l2) + l3 * l2y * Legendre1(l3 - l2) + l3 * l2 * Legendre1x(l3 - l2) * (l3y - l2y));
    }

    static double gradleg_tri_p3_b1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n12 * (l3y * l2 * Legendre1(l3 - l2) + l3 * l2y * Legendre1(l3 - l2) + l3 * l2 * Legendre1x(l3 - l2) * (l3y - l2y));
    }

     // EDGE 2
    static double gradleg_tri_p3_b2_a(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n21 * (l1 * l3 * Legendre1(l1 - l3));
    }

    static double gradleg_tri_p3_b2_b(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n22 * (l1 * l3 * Legendre1(l1 - l3));
    }

    static double gradleg_tri_p3_b2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n21 * (l1x * l3 * Legendre1(l1 - l3) + l1 * l3x * Legendre1(l1 - l3) + l1 * l3 * Legendre1x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p3_b2_bx(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n22 * (l1x * l3 * Legendre1(l1 - l3) + l1 * l3x * Legendre1(l1 - l3) + l1 * l3 * Legendre1x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p3_b2_ay(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n21 * (l1y * l3 * Legendre1(l1 - l3) + l1 * l3y * Legendre1(l1 - l3) + l1 * l3 * Legendre1x(l1 - l3) * (l1y - l3y));
    }

    static double gradleg_tri_p3_b2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n22 * (l1y * l3 * Legendre1(l1 - l3) + l1 * l3y * Legendre1(l1 - l3) + l1 * l3 * Legendre1x(l1 - l3) * (l1y - l3y));
    }

     // EDGE 3
    static double gradleg_tri_p3_b3_a(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n31 * (l2 * l1 * Legendre1(l2 - l1));
    }

    static double gradleg_tri_p3_b3_b(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n32 * (l2 * l1 * Legendre1(l2 - l1));
    }

    static double gradleg_tri_p3_b3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n31 * (l2x * l1 * Legendre1(l2 - l1) + l2 * l1x * Legendre1(l2 - l1) + l2 * l1 * Legendre1x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p3_b3_bx(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n32 * (l2x * l1 * Legendre1(l2 - l1) + l2 * l1x * Legendre1(l2 - l1) + l2 * l1 * Legendre1x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p3_b3_ay(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n31 * (l2y * l1 * Legendre1(l2 - l1) + l2 * l1y * Legendre1(l2 - l1) + l2 * l1 * Legendre1x(l2 - l1) * (l2y - l1y));
    }

    static double gradleg_tri_p3_b3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n32 * (l2y * l1 * Legendre1(l2 - l1) + l2 * l1y * Legendre1(l2 - l1) + l2 * l1 * Legendre1x(l2 - l1) * (l2y - l1y));
    }

    /* Genuine BUBBLE - order 3 */

    static double gradleg_tri_b1_b1_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b1_b1_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b1_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b1_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b1_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b1_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b1_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b1_b1_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b1_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b1_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b1_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b1_2_ay(double x, double y)
    {
     return 0.0;
    }

    ///////////////////////////////// ORDER 4 //////////////////////////////////

    /* EDGE FUNCTIONS - order 4*/

     /* EDGE 1 */

    static double gradleg_tri_p4_e1_a_0(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return  (l3x * l2 * phi3(l3 - l2) + l3 * l2x * phi3(l3 - l2) + l3 * l2 * phi3x(l3 - l2) * (l3x - l2x)) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e1_a_1(double x, double y)
    {
     return -(gradleg_tri_p4_e1_a_0(x, y));
    }

    static double gradleg_tri_p4_e1_b_0(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return  (l3y * l2 * phi3(l3 - l2) + l3 * l2y * phi3(l3 - l2) + l3 * l2 * phi3x(l3 - l2) * (l3y - l2y)) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e1_b_1(double x, double y)
    {
     return -(gradleg_tri_p4_e1_b_0(x, y));
    }

    static double gradleg_tri_p4_e1_ax_0(double x, double y)
    {
     double l3, l3x, l2, l2x;
     double ker, kerx, kerxx;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     ker = phi3(l3 - l2); kerx = phi3x(l3 - l2) * (l3x - l2x); kerxx = phi3xx(l3 - l2) * sqr(l3x - l2x);
     return  (2.0 * l3x * l2x * ker + 2.0 * l3x * l2 * kerx + 2.0 * l3 * l2x * kerx + l3 * l2 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e1_ax_1(double x, double y)
    {
     return -(gradleg_tri_p4_e1_ax_0(x, y));
    }

    static double gradleg_tri_p4_e1_ay_0(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi3(l3 - l2); kerx = phi3x(l3 - l2) * (l3x - l2x);
     kery = phi3x(l3 - l2) * (l3y - l2y); kerxy = phi3xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e1_ay_1(double x, double y)
    {
     return -(gradleg_tri_p4_e1_ay_0(x, y));
    }

    static double gradleg_tri_p4_e1_bx_0(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi3(l3 - l2); kerx = phi3x(l3 - l2) * (l3x - l2x);
     kery = phi3x(l3 - l2) * (l3y - l2y); kerxy = phi3xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e1_bx_1(double x, double y)
    {
     return -(gradleg_tri_p4_e1_bx_0(x, y));
    }

    static double gradleg_tri_p4_e1_by_0(double x, double y)
    {
     double l3, l3y, l2, l2y;
     double ker, kery, keryy;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     ker = phi3(l3 - l2); kery = phi3x(l3 - l2) * (l3y - l2y); keryy = phi3xx(l3 - l2) * sqr(l3y - l2y);
     return  (2.0 * l3y * l2y * ker + 2.0 * l3y * l2 * kery + 2.0 * l3 * l2y * kery + l3 * l2 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e1_by_1(double x, double y)
    {
     return -(gradleg_tri_p4_e1_by_0(x, y));
    }

     /* EDGE 2 */

    static double gradleg_tri_p4_e2_a_0(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return  (l1x * l3 * phi3(l1 - l3) + l1 * l3x * phi3(l1 - l3) + l1 * l3 * phi3x(l1 - l3) * (l1x - l3x)) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e2_a_1(double x, double y)
    {
     return -(gradleg_tri_p4_e2_a_0(x, y));
    }

    static double gradleg_tri_p4_e2_b_0(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return  (l1y * l3 * phi3(l1 - l3) + l1 * l3y * phi3(l1 - l3) + l1 * l3 * phi3x(l1 - l3) * (l1y - l3y)) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e2_b_1(double x, double y)
    {
     return -(gradleg_tri_p4_e2_b_0(x, y));
    }

    static double gradleg_tri_p4_e2_ax_0(double x, double y)
    {
     double l1, l1x, l3, l3x;
     double ker, kerx, kerxx;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     ker = phi3(l1 - l3); kerx = phi3x(l1 - l3) * (l1x - l3x); kerxx = phi3xx(l1 - l3) * sqr(l1x - l3x);
     return  (2.0 * l1x * l3x * ker + 2.0 * l1x * l3 * kerx + 2.0 * l1 * l3x * kerx + l1 * l3 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e2_ax_1(double x, double y)
    {
     return -(gradleg_tri_p4_e2_ax_0(x, y));
    }

    static double gradleg_tri_p4_e2_ay_0(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi3(l1 - l3); kerx = phi3x(l1 - l3) * (l1x - l3x);
     kery = phi3x(l1 - l3) * (l1y - l3y); kerxy = phi3xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e2_ay_1(double x, double y)
    {
     return -(gradleg_tri_p4_e2_ay_0(x, y));
    }

    static double gradleg_tri_p4_e2_bx_0(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi3(l1 - l3); kerx = phi3x(l1 - l3) * (l1x - l3x);
     kery = phi3x(l1 - l3) * (l1y - l3y); kerxy = phi3xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e2_bx_1(double x, double y)
    {
     return -(gradleg_tri_p4_e2_bx_0(x, y));
    }

    static double gradleg_tri_p4_e2_by_0(double x, double y)
    {
     double l1, l1y, l3, l3y;
     double ker, kery, keryy;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     ker = phi3(l1 - l3); kery = phi3x(l1 - l3) * (l1y - l3y); keryy = phi3xx(l1 - l3) * sqr(l1y - l3y);
     return  (2.0 * l1y * l3y * ker + 2.0 * l1y * l3 * kery + 2.0 * l1 * l3y * kery + l1 * l3 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e2_by_1(double x, double y)
    {
     return -(gradleg_tri_p4_e2_by_0(x, y));
    }

     /* EDGE 3 */

    static double gradleg_tri_p4_e3_a_0(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return  (l2x * l1 * phi3(l2 - l1) + l2 * l1x * phi3(l2 - l1) + l2 * l1 * phi3x(l2 - l1) * (l2x - l1x)) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e3_a_1(double x, double y)
    {
     return -(gradleg_tri_p4_e3_a_0(x, y));
    }

    static double gradleg_tri_p4_e3_b_0(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return  (l2y * l1 * phi3(l2 - l1) + l2 * l1y * phi3(l2 - l1) + l2 * l1 * phi3x(l2 - l1) * (l2y - l1y)) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e3_b_1(double x, double y)
    {
     return -(gradleg_tri_p4_e3_b_0(x, y));
    }

    static double gradleg_tri_p4_e3_ax_0(double x, double y)
    {
     double l2, l2x, l1, l1x;
     double ker, kerx, kerxx;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     ker = phi3(l2 - l1); kerx = phi3x(l2 - l1) * (l2x - l1x); kerxx = phi3xx(l2 - l1) * sqr(l2x - l1x);
     return  (2.0 * l2x * l1x * ker + 2.0 * l2x * l1 * kerx + 2.0 * l2 * l1x * kerx + l2 * l1 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e3_ax_1(double x, double y)
    {
     return -(gradleg_tri_p4_e3_ax_0(x, y));
    }

    static double gradleg_tri_p4_e3_ay_0(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi3(l2 - l1); kerx = phi3x(l2 - l1) * (l2x - l1x);
     kery = phi3x(l2 - l1) * (l2y - l1y); kerxy = phi3xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e3_ay_1(double x, double y)
    {
     return -(gradleg_tri_p4_e3_ay_0(x, y));
    }

    static double gradleg_tri_p4_e3_bx_0(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi3(l2 - l1); kerx = phi3x(l2 - l1) * (l2x - l1x);
     kery = phi3x(l2 - l1) * (l2y - l1y); kerxy = phi3xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e3_bx_1(double x, double y)
    {
     return -(gradleg_tri_p4_e3_bx_0(x, y));
    }

    static double gradleg_tri_p4_e3_by_0(double x, double y)
    {
     double l2, l2y, l1, l1y;
     double ker, kery, keryy;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     ker = phi3(l2 - l1); kery = phi3x(l2 - l1) * (l2y - l1y); keryy = phi3xx(l2 - l1) * sqr(l2y - l1y);
     return  (2.0 * l2y * l1y * ker + 2.0 * l2y * l1 * kery + 2.0 * l2 * l1y * kery + l2 * l1 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p4_e3_by_1(double x, double y)
    {
     return -(gradleg_tri_p4_e3_by_0(x, y));
    }

    /* BUBBLE */

    /* Edge-based BUBBLE - order 4 */

     // EDGE 1
    static double gradleg_tri_p4_b1_a(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n11 * (l3 * l2 * Legendre2(l3 - l2));
    }

    static double gradleg_tri_p4_b1_b(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n12 * (l3 * l2 * Legendre2(l3 - l2));
    }

    static double gradleg_tri_p4_b1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n11 * (l3x * l2 * Legendre2(l3 - l2) + l3 * l2x * Legendre2(l3 - l2) + l3 * l2 * Legendre2x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p4_b1_bx(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n12 * (l3x * l2 * Legendre2(l3 - l2) + l3 * l2x * Legendre2(l3 - l2) + l3 * l2 * Legendre2x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p4_b1_ay(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n11 * (l3y * l2 * Legendre2(l3 - l2) + l3 * l2y * Legendre2(l3 - l2) + l3 * l2 * Legendre2x(l3 - l2) * (l3y - l2y));
    }

    static double gradleg_tri_p4_b1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n12 * (l3y * l2 * Legendre2(l3 - l2) + l3 * l2y * Legendre2(l3 - l2) + l3 * l2 * Legendre2x(l3 - l2) * (l3y - l2y));
    }

     // EDGE 2
    static double gradleg_tri_p4_b2_a(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n21 * (l1 * l3 * Legendre2(l1 - l3));
    }

    static double gradleg_tri_p4_b2_b(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n22 * (l1 * l3 * Legendre2(l1 - l3));
    }

    static double gradleg_tri_p4_b2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n21 * (l1x * l3 * Legendre2(l1 - l3) + l1 * l3x * Legendre2(l1 - l3) + l1 * l3 * Legendre2x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p4_b2_bx(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n22 * (l1x * l3 * Legendre2(l1 - l3) + l1 * l3x * Legendre2(l1 - l3) + l1 * l3 * Legendre2x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p4_b2_ay(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n21 * (l1y * l3 * Legendre2(l1 - l3) + l1 * l3y * Legendre2(l1 - l3) + l1 * l3 * Legendre2x(l1 - l3) * (l1y - l3y));
    }

    static double gradleg_tri_p4_b2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n22 * (l1y * l3 * Legendre2(l1 - l3) + l1 * l3y * Legendre2(l1 - l3) + l1 * l3 * Legendre2x(l1 - l3) * (l1y - l3y));
    }

     // EDGE 3
    static double gradleg_tri_p4_b3_a(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n31 * (l2 * l1 * Legendre2(l2 - l1));
    }

    static double gradleg_tri_p4_b3_b(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n32 * (l2 * l1 * Legendre2(l2 - l1));
    }

    static double gradleg_tri_p4_b3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n31 * (l2x * l1 * Legendre2(l2 - l1) + l2 * l1x * Legendre2(l2 - l1) + l2 * l1 * Legendre2x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p4_b3_bx(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n32 * (l2x * l1 * Legendre2(l2 - l1) + l2 * l1x * Legendre2(l2 - l1) + l2 * l1 * Legendre2x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p4_b3_ay(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n31 * (l2y * l1 * Legendre2(l2 - l1) + l2 * l1y * Legendre2(l2 - l1) + l2 * l1 * Legendre2x(l2 - l1) * (l2y - l1y));
    }

    static double gradleg_tri_p4_b3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n32 * (l2y * l1 * Legendre2(l2 - l1) + l2 * l1y * Legendre2(l2 - l1) + l2 * l1 * Legendre2x(l2 - l1) * (l2y - l1y));
    }

    /* Genuine BUBBLE - order 4 */

    static double gradleg_tri_b1_b2_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b1_b2_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b2_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b2_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b2_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b2_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b2_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b1_b2_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b2_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b2_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b2_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b2_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b1_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b2_b1_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b1_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b1_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b1_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b1_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b1_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b2_b1_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b1_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b1_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b1_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b1_2_ay(double x, double y)
    {
     return 0.0;
    }

    ///////////////////////////////// ORDER 5 //////////////////////////////////

    /* EDGE FUNCTIONS - order 5*/

     /* EDGE 1 */

    static double gradleg_tri_p5_e1_a(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return  (l3x * l2 * phi4(l3 - l2) + l3 * l2x * phi4(l3 - l2) + l3 * l2 * phi4x(l3 - l2) * (l3x - l2x)) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e1_b(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return  (l3y * l2 * phi4(l3 - l2) + l3 * l2y * phi4(l3 - l2) + l3 * l2 * phi4x(l3 - l2) * (l3y - l2y)) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     double ker, kerx, kerxx;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     ker = phi4(l3 - l2); kerx = phi4x(l3 - l2) * (l3x - l2x); kerxx = phi4xx(l3 - l2) * sqr(l3x - l2x);
     return  (2.0 * l3x * l2x * ker + 2.0 * l3x * l2 * kerx + 2.0 * l3 * l2x * kerx + l3 * l2 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e1_ay(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi4(l3 - l2); kerx = phi4x(l3 - l2) * (l3x - l2x);
     kery = phi4x(l3 - l2) * (l3y - l2y); kerxy = phi4xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e1_bx(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi4(l3 - l2); kerx = phi4x(l3 - l2) * (l3x - l2x);
     kery = phi4x(l3 - l2) * (l3y - l2y); kerxy = phi4xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     double ker, kery, keryy;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     ker = phi4(l3 - l2); kery = phi4x(l3 - l2) * (l3y - l2y); keryy = phi4xx(l3 - l2) * sqr(l3y - l2y);
     return  (2.0 * l3y * l2y * ker + 2.0 * l3y * l2 * kery + 2.0 * l3 * l2y * kery + l3 * l2 * keryy) / 1.0000000000000;
    }

     /* EDGE 2 */

    static double gradleg_tri_p5_e2_a(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return  (l1x * l3 * phi4(l1 - l3) + l1 * l3x * phi4(l1 - l3) + l1 * l3 * phi4x(l1 - l3) * (l1x - l3x)) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e2_b(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return  (l1y * l3 * phi4(l1 - l3) + l1 * l3y * phi4(l1 - l3) + l1 * l3 * phi4x(l1 - l3) * (l1y - l3y)) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     double ker, kerx, kerxx;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     ker = phi4(l1 - l3); kerx = phi4x(l1 - l3) * (l1x - l3x); kerxx = phi4xx(l1 - l3) * sqr(l1x - l3x);
     return  (2.0 * l1x * l3x * ker + 2.0 * l1x * l3 * kerx + 2.0 * l1 * l3x * kerx + l1 * l3 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e2_ay(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi4(l1 - l3); kerx = phi4x(l1 - l3) * (l1x - l3x);
     kery = phi4x(l1 - l3) * (l1y - l3y); kerxy = phi4xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e2_bx(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi4(l1 - l3); kerx = phi4x(l1 - l3) * (l1x - l3x);
     kery = phi4x(l1 - l3) * (l1y - l3y); kerxy = phi4xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     double ker, kery, keryy;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     ker = phi4(l1 - l3); kery = phi4x(l1 - l3) * (l1y - l3y); keryy = phi4xx(l1 - l3) * sqr(l1y - l3y);
     return  (2.0 * l1y * l3y * ker + 2.0 * l1y * l3 * kery + 2.0 * l1 * l3y * kery + l1 * l3 * keryy) / 1.0000000000000;
    }

     /* EDGE 3 */

    static double gradleg_tri_p5_e3_a(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return  (l2x * l1 * phi4(l2 - l1) + l2 * l1x * phi4(l2 - l1) + l2 * l1 * phi4x(l2 - l1) * (l2x - l1x)) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e3_b(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return  (l2y * l1 * phi4(l2 - l1) + l2 * l1y * phi4(l2 - l1) + l2 * l1 * phi4x(l2 - l1) * (l2y - l1y)) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     double ker, kerx, kerxx;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     ker = phi4(l2 - l1); kerx = phi4x(l2 - l1) * (l2x - l1x); kerxx = phi4xx(l2 - l1) * sqr(l2x - l1x);
     return  (2.0 * l2x * l1x * ker + 2.0 * l2x * l1 * kerx + 2.0 * l2 * l1x * kerx + l2 * l1 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e3_ay(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi4(l2 - l1); kerx = phi4x(l2 - l1) * (l2x - l1x);
     kery = phi4x(l2 - l1) * (l2y - l1y); kerxy = phi4xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e3_bx(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi4(l2 - l1); kerx = phi4x(l2 - l1) * (l2x - l1x);
     kery = phi4x(l2 - l1) * (l2y - l1y); kerxy = phi4xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p5_e3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     double ker, kery, keryy;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     ker = phi4(l2 - l1); kery = phi4x(l2 - l1) * (l2y - l1y); keryy = phi4xx(l2 - l1) * sqr(l2y - l1y);
     return  (2.0 * l2y * l1y * ker + 2.0 * l2y * l1 * kery + 2.0 * l2 * l1y * kery + l2 * l1 * keryy) / 1.0000000000000;
    }

    /* BUBBLE */

    /* Edge-based BUBBLE - order 5 */

     // EDGE 1
    static double gradleg_tri_p5_b1_a(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n11 * (l3 * l2 * Legendre3(l3 - l2));
    }

    static double gradleg_tri_p5_b1_b(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n12 * (l3 * l2 * Legendre3(l3 - l2));
    }

    static double gradleg_tri_p5_b1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n11 * (l3x * l2 * Legendre3(l3 - l2) + l3 * l2x * Legendre3(l3 - l2) + l3 * l2 * Legendre3x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p5_b1_bx(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n12 * (l3x * l2 * Legendre3(l3 - l2) + l3 * l2x * Legendre3(l3 - l2) + l3 * l2 * Legendre3x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p5_b1_ay(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n11 * (l3y * l2 * Legendre3(l3 - l2) + l3 * l2y * Legendre3(l3 - l2) + l3 * l2 * Legendre3x(l3 - l2) * (l3y - l2y));
    }

    static double gradleg_tri_p5_b1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n12 * (l3y * l2 * Legendre3(l3 - l2) + l3 * l2y * Legendre3(l3 - l2) + l3 * l2 * Legendre3x(l3 - l2) * (l3y - l2y));
    }

     // EDGE 2
    static double gradleg_tri_p5_b2_a(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n21 * (l1 * l3 * Legendre3(l1 - l3));
    }

    static double gradleg_tri_p5_b2_b(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n22 * (l1 * l3 * Legendre3(l1 - l3));
    }

    static double gradleg_tri_p5_b2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n21 * (l1x * l3 * Legendre3(l1 - l3) + l1 * l3x * Legendre3(l1 - l3) + l1 * l3 * Legendre3x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p5_b2_bx(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n22 * (l1x * l3 * Legendre3(l1 - l3) + l1 * l3x * Legendre3(l1 - l3) + l1 * l3 * Legendre3x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p5_b2_ay(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n21 * (l1y * l3 * Legendre3(l1 - l3) + l1 * l3y * Legendre3(l1 - l3) + l1 * l3 * Legendre3x(l1 - l3) * (l1y - l3y));
    }

    static double gradleg_tri_p5_b2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n22 * (l1y * l3 * Legendre3(l1 - l3) + l1 * l3y * Legendre3(l1 - l3) + l1 * l3 * Legendre3x(l1 - l3) * (l1y - l3y));
    }

     // EDGE 3
    static double gradleg_tri_p5_b3_a(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n31 * (l2 * l1 * Legendre3(l2 - l1));
    }

    static double gradleg_tri_p5_b3_b(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n32 * (l2 * l1 * Legendre3(l2 - l1));
    }

    static double gradleg_tri_p5_b3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n31 * (l2x * l1 * Legendre3(l2 - l1) + l2 * l1x * Legendre3(l2 - l1) + l2 * l1 * Legendre3x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p5_b3_bx(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n32 * (l2x * l1 * Legendre3(l2 - l1) + l2 * l1x * Legendre3(l2 - l1) + l2 * l1 * Legendre3x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p5_b3_ay(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n31 * (l2y * l1 * Legendre3(l2 - l1) + l2 * l1y * Legendre3(l2 - l1) + l2 * l1 * Legendre3x(l2 - l1) * (l2y - l1y));
    }

    static double gradleg_tri_p5_b3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n32 * (l2y * l1 * Legendre3(l2 - l1) + l2 * l1y * Legendre3(l2 - l1) + l2 * l1 * Legendre3x(l2 - l1) * (l2y - l1y));
    }

    /* Genuine BUBBLE - order 5 */

    static double gradleg_tri_b1_b3_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b1_b3_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b3_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b3_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b3_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b3_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b3_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b1_b3_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b3_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b3_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b3_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b3_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b2_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b2_b2_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b2_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b2_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b2_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b2_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b2_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b2_b2_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b2_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b2_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b2_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b2_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b1_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b3_b1_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b1_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b1_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b1_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b1_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b1_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b3_b1_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b1_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b1_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b1_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b1_2_ay(double x, double y)
    {
     return 0.0;
    }

    ///////////////////////////////// ORDER 6 //////////////////////////////////

    /* EDGE FUNCTIONS - order 6*/

     /* EDGE 1 */

    static double gradleg_tri_p6_e1_a_0(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return  (l3x * l2 * phi5(l3 - l2) + l3 * l2x * phi5(l3 - l2) + l3 * l2 * phi5x(l3 - l2) * (l3x - l2x)) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e1_a_1(double x, double y)
    {
     return -(gradleg_tri_p6_e1_a_0(x, y));
    }

    static double gradleg_tri_p6_e1_b_0(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return  (l3y * l2 * phi5(l3 - l2) + l3 * l2y * phi5(l3 - l2) + l3 * l2 * phi5x(l3 - l2) * (l3y - l2y)) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e1_b_1(double x, double y)
    {
     return -(gradleg_tri_p6_e1_b_0(x, y));
    }

    static double gradleg_tri_p6_e1_ax_0(double x, double y)
    {
     double l3, l3x, l2, l2x;
     double ker, kerx, kerxx;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     ker = phi5(l3 - l2); kerx = phi5x(l3 - l2) * (l3x - l2x); kerxx = phi5xx(l3 - l2) * sqr(l3x - l2x);
     return  (2.0 * l3x * l2x * ker + 2.0 * l3x * l2 * kerx + 2.0 * l3 * l2x * kerx + l3 * l2 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e1_ax_1(double x, double y)
    {
     return -(gradleg_tri_p6_e1_ax_0(x, y));
    }

    static double gradleg_tri_p6_e1_ay_0(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi5(l3 - l2); kerx = phi5x(l3 - l2) * (l3x - l2x);
     kery = phi5x(l3 - l2) * (l3y - l2y); kerxy = phi5xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e1_ay_1(double x, double y)
    {
     return -(gradleg_tri_p6_e1_ay_0(x, y));
    }

    static double gradleg_tri_p6_e1_bx_0(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi5(l3 - l2); kerx = phi5x(l3 - l2) * (l3x - l2x);
     kery = phi5x(l3 - l2) * (l3y - l2y); kerxy = phi5xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e1_bx_1(double x, double y)
    {
     return -(gradleg_tri_p6_e1_bx_0(x, y));
    }

    static double gradleg_tri_p6_e1_by_0(double x, double y)
    {
     double l3, l3y, l2, l2y;
     double ker, kery, keryy;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     ker = phi5(l3 - l2); kery = phi5x(l3 - l2) * (l3y - l2y); keryy = phi5xx(l3 - l2) * sqr(l3y - l2y);
     return  (2.0 * l3y * l2y * ker + 2.0 * l3y * l2 * kery + 2.0 * l3 * l2y * kery + l3 * l2 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e1_by_1(double x, double y)
    {
     return -(gradleg_tri_p6_e1_by_0(x, y));
    }

     /* EDGE 2 */

    static double gradleg_tri_p6_e2_a_0(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return  (l1x * l3 * phi5(l1 - l3) + l1 * l3x * phi5(l1 - l3) + l1 * l3 * phi5x(l1 - l3) * (l1x - l3x)) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e2_a_1(double x, double y)
    {
     return -(gradleg_tri_p6_e2_a_0(x, y));
    }

    static double gradleg_tri_p6_e2_b_0(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return  (l1y * l3 * phi5(l1 - l3) + l1 * l3y * phi5(l1 - l3) + l1 * l3 * phi5x(l1 - l3) * (l1y - l3y)) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e2_b_1(double x, double y)
    {
     return -(gradleg_tri_p6_e2_b_0(x, y));
    }

    static double gradleg_tri_p6_e2_ax_0(double x, double y)
    {
     double l1, l1x, l3, l3x;
     double ker, kerx, kerxx;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     ker = phi5(l1 - l3); kerx = phi5x(l1 - l3) * (l1x - l3x); kerxx = phi5xx(l1 - l3) * sqr(l1x - l3x);
     return  (2.0 * l1x * l3x * ker + 2.0 * l1x * l3 * kerx + 2.0 * l1 * l3x * kerx + l1 * l3 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e2_ax_1(double x, double y)
    {
     return -(gradleg_tri_p6_e2_ax_0(x, y));
    }

    static double gradleg_tri_p6_e2_ay_0(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi5(l1 - l3); kerx = phi5x(l1 - l3) * (l1x - l3x);
     kery = phi5x(l1 - l3) * (l1y - l3y); kerxy = phi5xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e2_ay_1(double x, double y)
    {
     return -(gradleg_tri_p6_e2_ay_0(x, y));
    }

    static double gradleg_tri_p6_e2_bx_0(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi5(l1 - l3); kerx = phi5x(l1 - l3) * (l1x - l3x);
     kery = phi5x(l1 - l3) * (l1y - l3y); kerxy = phi5xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e2_bx_1(double x, double y)
    {
     return -(gradleg_tri_p6_e2_bx_0(x, y));
    }

    static double gradleg_tri_p6_e2_by_0(double x, double y)
    {
     double l1, l1y, l3, l3y;
     double ker, kery, keryy;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     ker = phi5(l1 - l3); kery = phi5x(l1 - l3) * (l1y - l3y); keryy = phi5xx(l1 - l3) * sqr(l1y - l3y);
     return  (2.0 * l1y * l3y * ker + 2.0 * l1y * l3 * kery + 2.0 * l1 * l3y * kery + l1 * l3 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e2_by_1(double x, double y)
    {
     return -(gradleg_tri_p6_e2_by_0(x, y));
    }

     /* EDGE 3 */

    static double gradleg_tri_p6_e3_a_0(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return  (l2x * l1 * phi5(l2 - l1) + l2 * l1x * phi5(l2 - l1) + l2 * l1 * phi5x(l2 - l1) * (l2x - l1x)) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e3_a_1(double x, double y)
    {
     return -(gradleg_tri_p6_e3_a_0(x, y));
    }

    static double gradleg_tri_p6_e3_b_0(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return  (l2y * l1 * phi5(l2 - l1) + l2 * l1y * phi5(l2 - l1) + l2 * l1 * phi5x(l2 - l1) * (l2y - l1y)) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e3_b_1(double x, double y)
    {
     return -(gradleg_tri_p6_e3_b_0(x, y));
    }

    static double gradleg_tri_p6_e3_ax_0(double x, double y)
    {
     double l2, l2x, l1, l1x;
     double ker, kerx, kerxx;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     ker = phi5(l2 - l1); kerx = phi5x(l2 - l1) * (l2x - l1x); kerxx = phi5xx(l2 - l1) * sqr(l2x - l1x);
     return  (2.0 * l2x * l1x * ker + 2.0 * l2x * l1 * kerx + 2.0 * l2 * l1x * kerx + l2 * l1 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e3_ax_1(double x, double y)
    {
     return -(gradleg_tri_p6_e3_ax_0(x, y));
    }

    static double gradleg_tri_p6_e3_ay_0(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi5(l2 - l1); kerx = phi5x(l2 - l1) * (l2x - l1x);
     kery = phi5x(l2 - l1) * (l2y - l1y); kerxy = phi5xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e3_ay_1(double x, double y)
    {
     return -(gradleg_tri_p6_e3_ay_0(x, y));
    }

    static double gradleg_tri_p6_e3_bx_0(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi5(l2 - l1); kerx = phi5x(l2 - l1) * (l2x - l1x);
     kery = phi5x(l2 - l1) * (l2y - l1y); kerxy = phi5xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e3_bx_1(double x, double y)
    {
     return -(gradleg_tri_p6_e3_bx_0(x, y));
    }

    static double gradleg_tri_p6_e3_by_0(double x, double y)
    {
     double l2, l2y, l1, l1y;
     double ker, kery, keryy;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     ker = phi5(l2 - l1); kery = phi5x(l2 - l1) * (l2y - l1y); keryy = phi5xx(l2 - l1) * sqr(l2y - l1y);
     return  (2.0 * l2y * l1y * ker + 2.0 * l2y * l1 * kery + 2.0 * l2 * l1y * kery + l2 * l1 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p6_e3_by_1(double x, double y)
    {
     return -(gradleg_tri_p6_e3_by_0(x, y));
    }

    /* BUBBLE */

    /* Edge-based BUBBLE - order 6 */

     // EDGE 1
    static double gradleg_tri_p6_b1_a(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n11 * (l3 * l2 * Legendre4(l3 - l2));
    }

    static double gradleg_tri_p6_b1_b(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n12 * (l3 * l2 * Legendre4(l3 - l2));
    }

    static double gradleg_tri_p6_b1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n11 * (l3x * l2 * Legendre4(l3 - l2) + l3 * l2x * Legendre4(l3 - l2) + l3 * l2 * Legendre4x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p6_b1_bx(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n12 * (l3x * l2 * Legendre4(l3 - l2) + l3 * l2x * Legendre4(l3 - l2) + l3 * l2 * Legendre4x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p6_b1_ay(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n11 * (l3y * l2 * Legendre4(l3 - l2) + l3 * l2y * Legendre4(l3 - l2) + l3 * l2 * Legendre4x(l3 - l2) * (l3y - l2y));
    }

    static double gradleg_tri_p6_b1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n12 * (l3y * l2 * Legendre4(l3 - l2) + l3 * l2y * Legendre4(l3 - l2) + l3 * l2 * Legendre4x(l3 - l2) * (l3y - l2y));
    }

     // EDGE 2
    static double gradleg_tri_p6_b2_a(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n21 * (l1 * l3 * Legendre4(l1 - l3));
    }

    static double gradleg_tri_p6_b2_b(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n22 * (l1 * l3 * Legendre4(l1 - l3));
    }

    static double gradleg_tri_p6_b2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n21 * (l1x * l3 * Legendre4(l1 - l3) + l1 * l3x * Legendre4(l1 - l3) + l1 * l3 * Legendre4x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p6_b2_bx(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n22 * (l1x * l3 * Legendre4(l1 - l3) + l1 * l3x * Legendre4(l1 - l3) + l1 * l3 * Legendre4x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p6_b2_ay(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n21 * (l1y * l3 * Legendre4(l1 - l3) + l1 * l3y * Legendre4(l1 - l3) + l1 * l3 * Legendre4x(l1 - l3) * (l1y - l3y));
    }

    static double gradleg_tri_p6_b2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n22 * (l1y * l3 * Legendre4(l1 - l3) + l1 * l3y * Legendre4(l1 - l3) + l1 * l3 * Legendre4x(l1 - l3) * (l1y - l3y));
    }

     // EDGE 3
    static double gradleg_tri_p6_b3_a(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n31 * (l2 * l1 * Legendre4(l2 - l1));
    }

    static double gradleg_tri_p6_b3_b(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n32 * (l2 * l1 * Legendre4(l2 - l1));
    }

    static double gradleg_tri_p6_b3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n31 * (l2x * l1 * Legendre4(l2 - l1) + l2 * l1x * Legendre4(l2 - l1) + l2 * l1 * Legendre4x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p6_b3_bx(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n32 * (l2x * l1 * Legendre4(l2 - l1) + l2 * l1x * Legendre4(l2 - l1) + l2 * l1 * Legendre4x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p6_b3_ay(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n31 * (l2y * l1 * Legendre4(l2 - l1) + l2 * l1y * Legendre4(l2 - l1) + l2 * l1 * Legendre4x(l2 - l1) * (l2y - l1y));
    }

    static double gradleg_tri_p6_b3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n32 * (l2y * l1 * Legendre4(l2 - l1) + l2 * l1y * Legendre4(l2 - l1) + l2 * l1 * Legendre4x(l2 - l1) * (l2y - l1y));
    }

    /* Genuine BUBBLE - order 6 */

    static double gradleg_tri_b1_b4_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre3(l2 - l1);
    }

    static double gradleg_tri_b1_b4_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b4_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre3(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre3x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b4_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b4_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre3(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre3x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b4_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b4_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre3(l2 - l1);
    }

    static double gradleg_tri_b1_b4_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b4_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre3(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre3x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b4_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b4_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre3(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre3x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b4_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b3_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b2_b3_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b3_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b3_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b3_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b3_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b3_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b2_b3_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b3_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b3_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b3_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b3_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b2_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b3_b2_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b2_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b2_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b2_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b2_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b2_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b3_b2_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b2_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b2_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b2_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b2_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b1_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b4_b1_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b1_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre3x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b4_b1_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b1_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre3x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b4_b1_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b1_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b4_b1_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b1_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre3x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b4_b1_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b1_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre3x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b4_b1_2_ay(double x, double y)
    {
     return 0.0;
    }

    ///////////////////////////////// ORDER 7 //////////////////////////////////

    /* EDGE FUNCTIONS - order 7*/

     /* EDGE 1 */

    static double gradleg_tri_p7_e1_a(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return  (l3x * l2 * phi6(l3 - l2) + l3 * l2x * phi6(l3 - l2) + l3 * l2 * phi6x(l3 - l2) * (l3x - l2x)) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e1_b(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return  (l3y * l2 * phi6(l3 - l2) + l3 * l2y * phi6(l3 - l2) + l3 * l2 * phi6x(l3 - l2) * (l3y - l2y)) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     double ker, kerx, kerxx;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     ker = phi6(l3 - l2); kerx = phi6x(l3 - l2) * (l3x - l2x); kerxx = phi6xx(l3 - l2) * sqr(l3x - l2x);
     return  (2.0 * l3x * l2x * ker + 2.0 * l3x * l2 * kerx + 2.0 * l3 * l2x * kerx + l3 * l2 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e1_ay(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi6(l3 - l2); kerx = phi6x(l3 - l2) * (l3x - l2x);
     kery = phi6x(l3 - l2) * (l3y - l2y); kerxy = phi6xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e1_bx(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi6(l3 - l2); kerx = phi6x(l3 - l2) * (l3x - l2x);
     kery = phi6x(l3 - l2) * (l3y - l2y); kerxy = phi6xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     double ker, kery, keryy;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     ker = phi6(l3 - l2); kery = phi6x(l3 - l2) * (l3y - l2y); keryy = phi6xx(l3 - l2) * sqr(l3y - l2y);
     return  (2.0 * l3y * l2y * ker + 2.0 * l3y * l2 * kery + 2.0 * l3 * l2y * kery + l3 * l2 * keryy) / 1.0000000000000;
    }

     /* EDGE 2 */

    static double gradleg_tri_p7_e2_a(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return  (l1x * l3 * phi6(l1 - l3) + l1 * l3x * phi6(l1 - l3) + l1 * l3 * phi6x(l1 - l3) * (l1x - l3x)) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e2_b(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return  (l1y * l3 * phi6(l1 - l3) + l1 * l3y * phi6(l1 - l3) + l1 * l3 * phi6x(l1 - l3) * (l1y - l3y)) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     double ker, kerx, kerxx;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     ker = phi6(l1 - l3); kerx = phi6x(l1 - l3) * (l1x - l3x); kerxx = phi6xx(l1 - l3) * sqr(l1x - l3x);
     return  (2.0 * l1x * l3x * ker + 2.0 * l1x * l3 * kerx + 2.0 * l1 * l3x * kerx + l1 * l3 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e2_ay(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi6(l1 - l3); kerx = phi6x(l1 - l3) * (l1x - l3x);
     kery = phi6x(l1 - l3) * (l1y - l3y); kerxy = phi6xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e2_bx(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi6(l1 - l3); kerx = phi6x(l1 - l3) * (l1x - l3x);
     kery = phi6x(l1 - l3) * (l1y - l3y); kerxy = phi6xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     double ker, kery, keryy;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     ker = phi6(l1 - l3); kery = phi6x(l1 - l3) * (l1y - l3y); keryy = phi6xx(l1 - l3) * sqr(l1y - l3y);
     return  (2.0 * l1y * l3y * ker + 2.0 * l1y * l3 * kery + 2.0 * l1 * l3y * kery + l1 * l3 * keryy) / 1.0000000000000;
    }

     /* EDGE 3 */

    static double gradleg_tri_p7_e3_a(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return  (l2x * l1 * phi6(l2 - l1) + l2 * l1x * phi6(l2 - l1) + l2 * l1 * phi6x(l2 - l1) * (l2x - l1x)) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e3_b(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return  (l2y * l1 * phi6(l2 - l1) + l2 * l1y * phi6(l2 - l1) + l2 * l1 * phi6x(l2 - l1) * (l2y - l1y)) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     double ker, kerx, kerxx;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     ker = phi6(l2 - l1); kerx = phi6x(l2 - l1) * (l2x - l1x); kerxx = phi6xx(l2 - l1) * sqr(l2x - l1x);
     return  (2.0 * l2x * l1x * ker + 2.0 * l2x * l1 * kerx + 2.0 * l2 * l1x * kerx + l2 * l1 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e3_ay(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi6(l2 - l1); kerx = phi6x(l2 - l1) * (l2x - l1x);
     kery = phi6x(l2 - l1) * (l2y - l1y); kerxy = phi6xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e3_bx(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi6(l2 - l1); kerx = phi6x(l2 - l1) * (l2x - l1x);
     kery = phi6x(l2 - l1) * (l2y - l1y); kerxy = phi6xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p7_e3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     double ker, kery, keryy;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     ker = phi6(l2 - l1); kery = phi6x(l2 - l1) * (l2y - l1y); keryy = phi6xx(l2 - l1) * sqr(l2y - l1y);
     return  (2.0 * l2y * l1y * ker + 2.0 * l2y * l1 * kery + 2.0 * l2 * l1y * kery + l2 * l1 * keryy) / 1.0000000000000;
    }

    /* BUBBLE */

    /* Edge-based BUBBLE - order 7 */

     // EDGE 1
    static double gradleg_tri_p7_b1_a(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n11 * (l3 * l2 * Legendre5(l3 - l2));
    }

    static double gradleg_tri_p7_b1_b(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n12 * (l3 * l2 * Legendre5(l3 - l2));
    }

    static double gradleg_tri_p7_b1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n11 * (l3x * l2 * Legendre5(l3 - l2) + l3 * l2x * Legendre5(l3 - l2) + l3 * l2 * Legendre5x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p7_b1_bx(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n12 * (l3x * l2 * Legendre5(l3 - l2) + l3 * l2x * Legendre5(l3 - l2) + l3 * l2 * Legendre5x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p7_b1_ay(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n11 * (l3y * l2 * Legendre5(l3 - l2) + l3 * l2y * Legendre5(l3 - l2) + l3 * l2 * Legendre5x(l3 - l2) * (l3y - l2y));
    }

    static double gradleg_tri_p7_b1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n12 * (l3y * l2 * Legendre5(l3 - l2) + l3 * l2y * Legendre5(l3 - l2) + l3 * l2 * Legendre5x(l3 - l2) * (l3y - l2y));
    }

     // EDGE 2
    static double gradleg_tri_p7_b2_a(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n21 * (l1 * l3 * Legendre5(l1 - l3));
    }

    static double gradleg_tri_p7_b2_b(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n22 * (l1 * l3 * Legendre5(l1 - l3));
    }

    static double gradleg_tri_p7_b2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n21 * (l1x * l3 * Legendre5(l1 - l3) + l1 * l3x * Legendre5(l1 - l3) + l1 * l3 * Legendre5x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p7_b2_bx(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n22 * (l1x * l3 * Legendre5(l1 - l3) + l1 * l3x * Legendre5(l1 - l3) + l1 * l3 * Legendre5x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p7_b2_ay(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n21 * (l1y * l3 * Legendre5(l1 - l3) + l1 * l3y * Legendre5(l1 - l3) + l1 * l3 * Legendre5x(l1 - l3) * (l1y - l3y));
    }

    static double gradleg_tri_p7_b2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n22 * (l1y * l3 * Legendre5(l1 - l3) + l1 * l3y * Legendre5(l1 - l3) + l1 * l3 * Legendre5x(l1 - l3) * (l1y - l3y));
    }

     // EDGE 3
    static double gradleg_tri_p7_b3_a(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n31 * (l2 * l1 * Legendre5(l2 - l1));
    }

    static double gradleg_tri_p7_b3_b(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n32 * (l2 * l1 * Legendre5(l2 - l1));
    }

    static double gradleg_tri_p7_b3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n31 * (l2x * l1 * Legendre5(l2 - l1) + l2 * l1x * Legendre5(l2 - l1) + l2 * l1 * Legendre5x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p7_b3_bx(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n32 * (l2x * l1 * Legendre5(l2 - l1) + l2 * l1x * Legendre5(l2 - l1) + l2 * l1 * Legendre5x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p7_b3_ay(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n31 * (l2y * l1 * Legendre5(l2 - l1) + l2 * l1y * Legendre5(l2 - l1) + l2 * l1 * Legendre5x(l2 - l1) * (l2y - l1y));
    }

    static double gradleg_tri_p7_b3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n32 * (l2y * l1 * Legendre5(l2 - l1) + l2 * l1y * Legendre5(l2 - l1) + l2 * l1 * Legendre5x(l2 - l1) * (l2y - l1y));
    }

    /* Genuine BUBBLE - order 7 */

    static double gradleg_tri_b1_b5_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre4(l2 - l1);
    }

    static double gradleg_tri_b1_b5_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b5_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre4(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre4x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b5_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b5_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre4(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre4x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b5_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b5_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre4(l2 - l1);
    }

    static double gradleg_tri_b1_b5_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b5_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre4(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre4x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b5_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b5_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre4(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre4x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b5_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b4_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre3(l2 - l1);
    }

    static double gradleg_tri_b2_b4_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b4_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre3(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre3x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b4_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b4_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre3(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre3x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b4_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b4_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre3(l2 - l1);
    }

    static double gradleg_tri_b2_b4_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b4_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre3(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre3x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b4_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b4_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre3(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre3x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b4_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b3_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b3_b3_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b3_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b3_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b3_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b3_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b3_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b3_b3_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b3_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b3_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b3_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b3_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b2_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b4_b2_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b2_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre3x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b4_b2_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b2_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre3x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b4_b2_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b2_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b4_b2_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b2_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre3x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b4_b2_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b2_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre3x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b4_b2_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b1_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b5_b1_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b1_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre4x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b5_b1_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b1_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre4x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b5_b1_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b1_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b5_b1_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b1_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre4x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b5_b1_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b1_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre4x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b5_b1_2_ay(double x, double y)
    {
     return 0.0;
    }

    ///////////////////////////////// ORDER 8 //////////////////////////////////

    /* EDGE FUNCTIONS - order 8*/

     /* EDGE 1 */

    static double gradleg_tri_p8_e1_a_0(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return  (l3x * l2 * phi7(l3 - l2) + l3 * l2x * phi7(l3 - l2) + l3 * l2 * phi7x(l3 - l2) * (l3x - l2x)) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e1_a_1(double x, double y)
    {
     return -(gradleg_tri_p8_e1_a_0(x, y));
    }

    static double gradleg_tri_p8_e1_b_0(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return  (l3y * l2 * phi7(l3 - l2) + l3 * l2y * phi7(l3 - l2) + l3 * l2 * phi7x(l3 - l2) * (l3y - l2y)) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e1_b_1(double x, double y)
    {
     return -(gradleg_tri_p8_e1_b_0(x, y));
    }

    static double gradleg_tri_p8_e1_ax_0(double x, double y)
    {
     double l3, l3x, l2, l2x;
     double ker, kerx, kerxx;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     ker = phi7(l3 - l2); kerx = phi7x(l3 - l2) * (l3x - l2x); kerxx = phi7xx(l3 - l2) * sqr(l3x - l2x);
     return  (2.0 * l3x * l2x * ker + 2.0 * l3x * l2 * kerx + 2.0 * l3 * l2x * kerx + l3 * l2 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e1_ax_1(double x, double y)
    {
     return -(gradleg_tri_p8_e1_ax_0(x, y));
    }

    static double gradleg_tri_p8_e1_ay_0(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi7(l3 - l2); kerx = phi7x(l3 - l2) * (l3x - l2x);
     kery = phi7x(l3 - l2) * (l3y - l2y); kerxy = phi7xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e1_ay_1(double x, double y)
    {
     return -(gradleg_tri_p8_e1_ay_0(x, y));
    }

    static double gradleg_tri_p8_e1_bx_0(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi7(l3 - l2); kerx = phi7x(l3 - l2) * (l3x - l2x);
     kery = phi7x(l3 - l2) * (l3y - l2y); kerxy = phi7xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e1_bx_1(double x, double y)
    {
     return -(gradleg_tri_p8_e1_bx_0(x, y));
    }

    static double gradleg_tri_p8_e1_by_0(double x, double y)
    {
     double l3, l3y, l2, l2y;
     double ker, kery, keryy;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     ker = phi7(l3 - l2); kery = phi7x(l3 - l2) * (l3y - l2y); keryy = phi7xx(l3 - l2) * sqr(l3y - l2y);
     return  (2.0 * l3y * l2y * ker + 2.0 * l3y * l2 * kery + 2.0 * l3 * l2y * kery + l3 * l2 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e1_by_1(double x, double y)
    {
     return -(gradleg_tri_p8_e1_by_0(x, y));
    }

     /* EDGE 2 */

    static double gradleg_tri_p8_e2_a_0(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return  (l1x * l3 * phi7(l1 - l3) + l1 * l3x * phi7(l1 - l3) + l1 * l3 * phi7x(l1 - l3) * (l1x - l3x)) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e2_a_1(double x, double y)
    {
     return -(gradleg_tri_p8_e2_a_0(x, y));
    }

    static double gradleg_tri_p8_e2_b_0(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return  (l1y * l3 * phi7(l1 - l3) + l1 * l3y * phi7(l1 - l3) + l1 * l3 * phi7x(l1 - l3) * (l1y - l3y)) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e2_b_1(double x, double y)
    {
     return -(gradleg_tri_p8_e2_b_0(x, y));
    }

    static double gradleg_tri_p8_e2_ax_0(double x, double y)
    {
     double l1, l1x, l3, l3x;
     double ker, kerx, kerxx;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     ker = phi7(l1 - l3); kerx = phi7x(l1 - l3) * (l1x - l3x); kerxx = phi7xx(l1 - l3) * sqr(l1x - l3x);
     return  (2.0 * l1x * l3x * ker + 2.0 * l1x * l3 * kerx + 2.0 * l1 * l3x * kerx + l1 * l3 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e2_ax_1(double x, double y)
    {
     return -(gradleg_tri_p8_e2_ax_0(x, y));
    }

    static double gradleg_tri_p8_e2_ay_0(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi7(l1 - l3); kerx = phi7x(l1 - l3) * (l1x - l3x);
     kery = phi7x(l1 - l3) * (l1y - l3y); kerxy = phi7xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e2_ay_1(double x, double y)
    {
     return -(gradleg_tri_p8_e2_ay_0(x, y));
    }

    static double gradleg_tri_p8_e2_bx_0(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi7(l1 - l3); kerx = phi7x(l1 - l3) * (l1x - l3x);
     kery = phi7x(l1 - l3) * (l1y - l3y); kerxy = phi7xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e2_bx_1(double x, double y)
    {
     return -(gradleg_tri_p8_e2_bx_0(x, y));
    }

    static double gradleg_tri_p8_e2_by_0(double x, double y)
    {
     double l1, l1y, l3, l3y;
     double ker, kery, keryy;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     ker = phi7(l1 - l3); kery = phi7x(l1 - l3) * (l1y - l3y); keryy = phi7xx(l1 - l3) * sqr(l1y - l3y);
     return  (2.0 * l1y * l3y * ker + 2.0 * l1y * l3 * kery + 2.0 * l1 * l3y * kery + l1 * l3 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e2_by_1(double x, double y)
    {
     return -(gradleg_tri_p8_e2_by_0(x, y));
    }

     /* EDGE 3 */

    static double gradleg_tri_p8_e3_a_0(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return  (l2x * l1 * phi7(l2 - l1) + l2 * l1x * phi7(l2 - l1) + l2 * l1 * phi7x(l2 - l1) * (l2x - l1x)) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e3_a_1(double x, double y)
    {
     return -(gradleg_tri_p8_e3_a_0(x, y));
    }

    static double gradleg_tri_p8_e3_b_0(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return  (l2y * l1 * phi7(l2 - l1) + l2 * l1y * phi7(l2 - l1) + l2 * l1 * phi7x(l2 - l1) * (l2y - l1y)) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e3_b_1(double x, double y)
    {
     return -(gradleg_tri_p8_e3_b_0(x, y));
    }

    static double gradleg_tri_p8_e3_ax_0(double x, double y)
    {
     double l2, l2x, l1, l1x;
     double ker, kerx, kerxx;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     ker = phi7(l2 - l1); kerx = phi7x(l2 - l1) * (l2x - l1x); kerxx = phi7xx(l2 - l1) * sqr(l2x - l1x);
     return  (2.0 * l2x * l1x * ker + 2.0 * l2x * l1 * kerx + 2.0 * l2 * l1x * kerx + l2 * l1 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e3_ax_1(double x, double y)
    {
     return -(gradleg_tri_p8_e3_ax_0(x, y));
    }

    static double gradleg_tri_p8_e3_ay_0(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi7(l2 - l1); kerx = phi7x(l2 - l1) * (l2x - l1x);
     kery = phi7x(l2 - l1) * (l2y - l1y); kerxy = phi7xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e3_ay_1(double x, double y)
    {
     return -(gradleg_tri_p8_e3_ay_0(x, y));
    }

    static double gradleg_tri_p8_e3_bx_0(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi7(l2 - l1); kerx = phi7x(l2 - l1) * (l2x - l1x);
     kery = phi7x(l2 - l1) * (l2y - l1y); kerxy = phi7xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e3_bx_1(double x, double y)
    {
     return -(gradleg_tri_p8_e3_bx_0(x, y));
    }

    static double gradleg_tri_p8_e3_by_0(double x, double y)
    {
     double l2, l2y, l1, l1y;
     double ker, kery, keryy;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     ker = phi7(l2 - l1); kery = phi7x(l2 - l1) * (l2y - l1y); keryy = phi7xx(l2 - l1) * sqr(l2y - l1y);
     return  (2.0 * l2y * l1y * ker + 2.0 * l2y * l1 * kery + 2.0 * l2 * l1y * kery + l2 * l1 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p8_e3_by_1(double x, double y)
    {
     return -(gradleg_tri_p8_e3_by_0(x, y));
    }

    /* BUBBLE */

    /* Edge-based BUBBLE - order 8 */

     // EDGE 1
    static double gradleg_tri_p8_b1_a(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n11 * (l3 * l2 * Legendre6(l3 - l2));
    }

    static double gradleg_tri_p8_b1_b(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n12 * (l3 * l2 * Legendre6(l3 - l2));
    }

    static double gradleg_tri_p8_b1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n11 * (l3x * l2 * Legendre6(l3 - l2) + l3 * l2x * Legendre6(l3 - l2) + l3 * l2 * Legendre6x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p8_b1_bx(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n12 * (l3x * l2 * Legendre6(l3 - l2) + l3 * l2x * Legendre6(l3 - l2) + l3 * l2 * Legendre6x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p8_b1_ay(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n11 * (l3y * l2 * Legendre6(l3 - l2) + l3 * l2y * Legendre6(l3 - l2) + l3 * l2 * Legendre6x(l3 - l2) * (l3y - l2y));
    }

    static double gradleg_tri_p8_b1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n12 * (l3y * l2 * Legendre6(l3 - l2) + l3 * l2y * Legendre6(l3 - l2) + l3 * l2 * Legendre6x(l3 - l2) * (l3y - l2y));
    }

     // EDGE 2
    static double gradleg_tri_p8_b2_a(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n21 * (l1 * l3 * Legendre6(l1 - l3));
    }

    static double gradleg_tri_p8_b2_b(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n22 * (l1 * l3 * Legendre6(l1 - l3));
    }

    static double gradleg_tri_p8_b2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n21 * (l1x * l3 * Legendre6(l1 - l3) + l1 * l3x * Legendre6(l1 - l3) + l1 * l3 * Legendre6x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p8_b2_bx(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n22 * (l1x * l3 * Legendre6(l1 - l3) + l1 * l3x * Legendre6(l1 - l3) + l1 * l3 * Legendre6x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p8_b2_ay(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n21 * (l1y * l3 * Legendre6(l1 - l3) + l1 * l3y * Legendre6(l1 - l3) + l1 * l3 * Legendre6x(l1 - l3) * (l1y - l3y));
    }

    static double gradleg_tri_p8_b2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n22 * (l1y * l3 * Legendre6(l1 - l3) + l1 * l3y * Legendre6(l1 - l3) + l1 * l3 * Legendre6x(l1 - l3) * (l1y - l3y));
    }

     // EDGE 3
    static double gradleg_tri_p8_b3_a(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n31 * (l2 * l1 * Legendre6(l2 - l1));
    }

    static double gradleg_tri_p8_b3_b(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n32 * (l2 * l1 * Legendre6(l2 - l1));
    }

    static double gradleg_tri_p8_b3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n31 * (l2x * l1 * Legendre6(l2 - l1) + l2 * l1x * Legendre6(l2 - l1) + l2 * l1 * Legendre6x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p8_b3_bx(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n32 * (l2x * l1 * Legendre6(l2 - l1) + l2 * l1x * Legendre6(l2 - l1) + l2 * l1 * Legendre6x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p8_b3_ay(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n31 * (l2y * l1 * Legendre6(l2 - l1) + l2 * l1y * Legendre6(l2 - l1) + l2 * l1 * Legendre6x(l2 - l1) * (l2y - l1y));
    }

    static double gradleg_tri_p8_b3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n32 * (l2y * l1 * Legendre6(l2 - l1) + l2 * l1y * Legendre6(l2 - l1) + l2 * l1 * Legendre6x(l2 - l1) * (l2y - l1y));
    }

    /* Genuine BUBBLE - order 8 */

    static double gradleg_tri_b1_b6_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre5(l2 - l1);
    }

    static double gradleg_tri_b1_b6_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b6_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre5(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre5x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b6_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b6_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre5(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre5x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b6_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b6_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre5(l2 - l1);
    }

    static double gradleg_tri_b1_b6_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b6_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre5(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre5x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b6_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b6_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre5(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre5x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b6_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b5_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre4(l2 - l1);
    }

    static double gradleg_tri_b2_b5_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b5_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre4(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre4x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b5_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b5_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre4(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre4x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b5_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b5_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre4(l2 - l1);
    }

    static double gradleg_tri_b2_b5_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b5_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre4(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre4x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b5_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b5_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre4(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre4x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b5_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b4_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre3(l2 - l1);
    }

    static double gradleg_tri_b3_b4_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b4_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre3(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre3x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b4_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b4_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre3(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre3x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b4_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b4_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre3(l2 - l1);
    }

    static double gradleg_tri_b3_b4_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b4_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre3(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre3x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b4_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b4_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre3(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre3x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b4_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b3_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b4_b3_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b3_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre3x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b4_b3_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b3_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre3x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b4_b3_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b3_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b4_b3_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b3_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre3x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b4_b3_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b3_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre3x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b4_b3_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b2_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b5_b2_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b2_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre4x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b5_b2_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b2_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre4x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b5_b2_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b2_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b5_b2_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b2_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre4x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b5_b2_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b2_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre4x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b5_b2_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b1_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b6_b1_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b1_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre5x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b6_b1_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b1_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre5x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b6_b1_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b1_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b6_b1_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b1_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre5x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b6_b1_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b1_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre5x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b6_b1_2_ay(double x, double y)
    {
     return 0.0;
    }

    ///////////////////////////////// ORDER 9 //////////////////////////////////

    /* EDGE FUNCTIONS - order 9*/

     /* EDGE 1 */

    static double gradleg_tri_p9_e1_a(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return  (l3x * l2 * phi8(l3 - l2) + l3 * l2x * phi8(l3 - l2) + l3 * l2 * phi8x(l3 - l2) * (l3x - l2x)) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e1_b(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return  (l3y * l2 * phi8(l3 - l2) + l3 * l2y * phi8(l3 - l2) + l3 * l2 * phi8x(l3 - l2) * (l3y - l2y)) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     double ker, kerx, kerxx;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     ker = phi8(l3 - l2); kerx = phi8x(l3 - l2) * (l3x - l2x); kerxx = phi8xx(l3 - l2) * sqr(l3x - l2x);
     return  (2.0 * l3x * l2x * ker + 2.0 * l3x * l2 * kerx + 2.0 * l3 * l2x * kerx + l3 * l2 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e1_ay(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi8(l3 - l2); kerx = phi8x(l3 - l2) * (l3x - l2x);
     kery = phi8x(l3 - l2) * (l3y - l2y); kerxy = phi8xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e1_bx(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi8(l3 - l2); kerx = phi8x(l3 - l2) * (l3x - l2x);
     kery = phi8x(l3 - l2) * (l3y - l2y); kerxy = phi8xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     double ker, kery, keryy;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     ker = phi8(l3 - l2); kery = phi8x(l3 - l2) * (l3y - l2y); keryy = phi8xx(l3 - l2) * sqr(l3y - l2y);
     return  (2.0 * l3y * l2y * ker + 2.0 * l3y * l2 * kery + 2.0 * l3 * l2y * kery + l3 * l2 * keryy) / 1.0000000000000;
    }

     /* EDGE 2 */

    static double gradleg_tri_p9_e2_a(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return  (l1x * l3 * phi8(l1 - l3) + l1 * l3x * phi8(l1 - l3) + l1 * l3 * phi8x(l1 - l3) * (l1x - l3x)) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e2_b(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return  (l1y * l3 * phi8(l1 - l3) + l1 * l3y * phi8(l1 - l3) + l1 * l3 * phi8x(l1 - l3) * (l1y - l3y)) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     double ker, kerx, kerxx;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     ker = phi8(l1 - l3); kerx = phi8x(l1 - l3) * (l1x - l3x); kerxx = phi8xx(l1 - l3) * sqr(l1x - l3x);
     return  (2.0 * l1x * l3x * ker + 2.0 * l1x * l3 * kerx + 2.0 * l1 * l3x * kerx + l1 * l3 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e2_ay(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi8(l1 - l3); kerx = phi8x(l1 - l3) * (l1x - l3x);
     kery = phi8x(l1 - l3) * (l1y - l3y); kerxy = phi8xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e2_bx(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi8(l1 - l3); kerx = phi8x(l1 - l3) * (l1x - l3x);
     kery = phi8x(l1 - l3) * (l1y - l3y); kerxy = phi8xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     double ker, kery, keryy;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     ker = phi8(l1 - l3); kery = phi8x(l1 - l3) * (l1y - l3y); keryy = phi8xx(l1 - l3) * sqr(l1y - l3y);
     return  (2.0 * l1y * l3y * ker + 2.0 * l1y * l3 * kery + 2.0 * l1 * l3y * kery + l1 * l3 * keryy) / 1.0000000000000;
    }

     /* EDGE 3 */

    static double gradleg_tri_p9_e3_a(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return  (l2x * l1 * phi8(l2 - l1) + l2 * l1x * phi8(l2 - l1) + l2 * l1 * phi8x(l2 - l1) * (l2x - l1x)) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e3_b(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return  (l2y * l1 * phi8(l2 - l1) + l2 * l1y * phi8(l2 - l1) + l2 * l1 * phi8x(l2 - l1) * (l2y - l1y)) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     double ker, kerx, kerxx;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     ker = phi8(l2 - l1); kerx = phi8x(l2 - l1) * (l2x - l1x); kerxx = phi8xx(l2 - l1) * sqr(l2x - l1x);
     return  (2.0 * l2x * l1x * ker + 2.0 * l2x * l1 * kerx + 2.0 * l2 * l1x * kerx + l2 * l1 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e3_ay(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi8(l2 - l1); kerx = phi8x(l2 - l1) * (l2x - l1x);
     kery = phi8x(l2 - l1) * (l2y - l1y); kerxy = phi8xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e3_bx(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi8(l2 - l1); kerx = phi8x(l2 - l1) * (l2x - l1x);
     kery = phi8x(l2 - l1) * (l2y - l1y); kerxy = phi8xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p9_e3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     double ker, kery, keryy;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     ker = phi8(l2 - l1); kery = phi8x(l2 - l1) * (l2y - l1y); keryy = phi8xx(l2 - l1) * sqr(l2y - l1y);
     return  (2.0 * l2y * l1y * ker + 2.0 * l2y * l1 * kery + 2.0 * l2 * l1y * kery + l2 * l1 * keryy) / 1.0000000000000;
    }

    /* BUBBLE */

    /* Edge-based BUBBLE - order 9 */

     // EDGE 1
    static double gradleg_tri_p9_b1_a(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n11 * (l3 * l2 * Legendre7(l3 - l2));
    }

    static double gradleg_tri_p9_b1_b(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n12 * (l3 * l2 * Legendre7(l3 - l2));
    }

    static double gradleg_tri_p9_b1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n11 * (l3x * l2 * Legendre7(l3 - l2) + l3 * l2x * Legendre7(l3 - l2) + l3 * l2 * Legendre7x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p9_b1_bx(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n12 * (l3x * l2 * Legendre7(l3 - l2) + l3 * l2x * Legendre7(l3 - l2) + l3 * l2 * Legendre7x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p9_b1_ay(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n11 * (l3y * l2 * Legendre7(l3 - l2) + l3 * l2y * Legendre7(l3 - l2) + l3 * l2 * Legendre7x(l3 - l2) * (l3y - l2y));
    }

    static double gradleg_tri_p9_b1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n12 * (l3y * l2 * Legendre7(l3 - l2) + l3 * l2y * Legendre7(l3 - l2) + l3 * l2 * Legendre7x(l3 - l2) * (l3y - l2y));
    }

     // EDGE 2
    static double gradleg_tri_p9_b2_a(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n21 * (l1 * l3 * Legendre7(l1 - l3));
    }

    static double gradleg_tri_p9_b2_b(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n22 * (l1 * l3 * Legendre7(l1 - l3));
    }

    static double gradleg_tri_p9_b2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n21 * (l1x * l3 * Legendre7(l1 - l3) + l1 * l3x * Legendre7(l1 - l3) + l1 * l3 * Legendre7x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p9_b2_bx(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n22 * (l1x * l3 * Legendre7(l1 - l3) + l1 * l3x * Legendre7(l1 - l3) + l1 * l3 * Legendre7x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p9_b2_ay(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n21 * (l1y * l3 * Legendre7(l1 - l3) + l1 * l3y * Legendre7(l1 - l3) + l1 * l3 * Legendre7x(l1 - l3) * (l1y - l3y));
    }

    static double gradleg_tri_p9_b2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n22 * (l1y * l3 * Legendre7(l1 - l3) + l1 * l3y * Legendre7(l1 - l3) + l1 * l3 * Legendre7x(l1 - l3) * (l1y - l3y));
    }

     // EDGE 3
    static double gradleg_tri_p9_b3_a(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n31 * (l2 * l1 * Legendre7(l2 - l1));
    }

    static double gradleg_tri_p9_b3_b(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n32 * (l2 * l1 * Legendre7(l2 - l1));
    }

    static double gradleg_tri_p9_b3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n31 * (l2x * l1 * Legendre7(l2 - l1) + l2 * l1x * Legendre7(l2 - l1) + l2 * l1 * Legendre7x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p9_b3_bx(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n32 * (l2x * l1 * Legendre7(l2 - l1) + l2 * l1x * Legendre7(l2 - l1) + l2 * l1 * Legendre7x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p9_b3_ay(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n31 * (l2y * l1 * Legendre7(l2 - l1) + l2 * l1y * Legendre7(l2 - l1) + l2 * l1 * Legendre7x(l2 - l1) * (l2y - l1y));
    }

    static double gradleg_tri_p9_b3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n32 * (l2y * l1 * Legendre7(l2 - l1) + l2 * l1y * Legendre7(l2 - l1) + l2 * l1 * Legendre7x(l2 - l1) * (l2y - l1y));
    }

    /* Genuine BUBBLE - order 9 */

    static double gradleg_tri_b1_b7_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre6(l2 - l1);
    }

    static double gradleg_tri_b1_b7_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b7_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre6(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre6x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b7_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b7_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre6(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre6x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b7_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b7_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre6(l2 - l1);
    }

    static double gradleg_tri_b1_b7_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b7_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre6(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre6x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b7_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b7_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre6(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre6x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b7_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b6_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre5(l2 - l1);
    }

    static double gradleg_tri_b2_b6_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b6_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre5(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre5x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b6_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b6_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre5(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre5x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b6_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b6_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre5(l2 - l1);
    }

    static double gradleg_tri_b2_b6_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b6_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre5(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre5x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b6_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b6_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre5(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre5x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b6_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b5_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre4(l2 - l1);
    }

    static double gradleg_tri_b3_b5_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b5_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre4(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre4x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b5_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b5_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre4(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre4x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b5_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b5_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre4(l2 - l1);
    }

    static double gradleg_tri_b3_b5_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b5_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre4(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre4x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b5_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b5_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre4(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre4x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b5_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b4_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre3(l2 - l1);
    }

    static double gradleg_tri_b4_b4_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b4_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre3(l2 - l1);
     L1x = Legendre3x(l3 - l2) * (l3x - l2x); L2x = Legendre3x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b4_b4_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b4_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre3(l2 - l1);
     L1y = Legendre3x(l3 - l2) * (l3y - l2y); L2y = Legendre3x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b4_b4_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b4_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre3(l2 - l1);
    }

    static double gradleg_tri_b4_b4_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b4_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre3(l2 - l1);
     L1x = Legendre3x(l3 - l2) * (l3x - l2x); L2x = Legendre3x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b4_b4_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b4_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre3(l2 - l1);
     L1y = Legendre3x(l3 - l2) * (l3y - l2y); L2y = Legendre3x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b4_b4_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b3_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b5_b3_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b3_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre4x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b5_b3_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b3_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre4x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b5_b3_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b3_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b5_b3_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b3_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre4x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b5_b3_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b3_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre4x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b5_b3_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b2_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b6_b2_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b2_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre5x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b6_b2_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b2_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre5x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b6_b2_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b2_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b6_b2_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b2_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre5x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b6_b2_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b2_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre5x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b6_b2_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b1_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre6(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b7_b1_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b1_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre6(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre6x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b7_b1_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b1_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre6(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre6x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b7_b1_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b1_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre6(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b7_b1_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b1_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre6(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre6x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b7_b1_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b1_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre6(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre6x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b7_b1_2_ay(double x, double y)
    {
     return 0.0;
    }

    ///////////////////////////////// ORDER 10 //////////////////////////////////

    /* EDGE FUNCTIONS - order 10*/

     /* EDGE 1 */

    static double gradleg_tri_p10_e1_a_0(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return  (l3x * l2 * phi9(l3 - l2) + l3 * l2x * phi9(l3 - l2) + l3 * l2 * phi9x(l3 - l2) * (l3x - l2x)) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e1_a_1(double x, double y)
    {
     return -(gradleg_tri_p10_e1_a_0(x, y));
    }

    static double gradleg_tri_p10_e1_b_0(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return  (l3y * l2 * phi9(l3 - l2) + l3 * l2y * phi9(l3 - l2) + l3 * l2 * phi9x(l3 - l2) * (l3y - l2y)) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e1_b_1(double x, double y)
    {
     return -(gradleg_tri_p10_e1_b_0(x, y));
    }

    static double gradleg_tri_p10_e1_ax_0(double x, double y)
    {
     double l3, l3x, l2, l2x;
     double ker, kerx, kerxx;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     ker = phi9(l3 - l2); kerx = phi9x(l3 - l2) * (l3x - l2x); kerxx = phi9xx(l3 - l2) * sqr(l3x - l2x);
     return  (2.0 * l3x * l2x * ker + 2.0 * l3x * l2 * kerx + 2.0 * l3 * l2x * kerx + l3 * l2 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e1_ax_1(double x, double y)
    {
     return -(gradleg_tri_p10_e1_ax_0(x, y));
    }

    static double gradleg_tri_p10_e1_ay_0(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi9(l3 - l2); kerx = phi9x(l3 - l2) * (l3x - l2x);
     kery = phi9x(l3 - l2) * (l3y - l2y); kerxy = phi9xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e1_ay_1(double x, double y)
    {
     return -(gradleg_tri_p10_e1_ay_0(x, y));
    }

    static double gradleg_tri_p10_e1_bx_0(double x, double y)
    {
     double l3, l3x, l3y, l2, l2x, l2y;
     double ker, kerx, kery, kerxy;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
      l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
     ker = phi9(l3 - l2); kerx = phi9x(l3 - l2) * (l3x - l2x);
     kery = phi9x(l3 - l2) * (l3y - l2y); kerxy = phi9xx(l3 - l2) * (l3x - l2x) * (l3y - l2y);
     return  (l3x * l2y * ker + l3y * l2x * ker + l3x * l2 * kery + l3 * l2x * kery + l3y * l2 * kerx + l3 * l2y * kerx +  l3 * l2 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e1_bx_1(double x, double y)
    {
     return -(gradleg_tri_p10_e1_bx_0(x, y));
    }

    static double gradleg_tri_p10_e1_by_0(double x, double y)
    {
     double l3, l3y, l2, l2y;
     double ker, kery, keryy;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     ker = phi9(l3 - l2); kery = phi9x(l3 - l2) * (l3y - l2y); keryy = phi9xx(l3 - l2) * sqr(l3y - l2y);
     return  (2.0 * l3y * l2y * ker + 2.0 * l3y * l2 * kery + 2.0 * l3 * l2y * kery + l3 * l2 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e1_by_1(double x, double y)
    {
     return -(gradleg_tri_p10_e1_by_0(x, y));
    }

     /* EDGE 2 */

    static double gradleg_tri_p10_e2_a_0(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return  (l1x * l3 * phi9(l1 - l3) + l1 * l3x * phi9(l1 - l3) + l1 * l3 * phi9x(l1 - l3) * (l1x - l3x)) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e2_a_1(double x, double y)
    {
     return -(gradleg_tri_p10_e2_a_0(x, y));
    }

    static double gradleg_tri_p10_e2_b_0(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return  (l1y * l3 * phi9(l1 - l3) + l1 * l3y * phi9(l1 - l3) + l1 * l3 * phi9x(l1 - l3) * (l1y - l3y)) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e2_b_1(double x, double y)
    {
     return -(gradleg_tri_p10_e2_b_0(x, y));
    }

    static double gradleg_tri_p10_e2_ax_0(double x, double y)
    {
     double l1, l1x, l3, l3x;
     double ker, kerx, kerxx;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     ker = phi9(l1 - l3); kerx = phi9x(l1 - l3) * (l1x - l3x); kerxx = phi9xx(l1 - l3) * sqr(l1x - l3x);
     return  (2.0 * l1x * l3x * ker + 2.0 * l1x * l3 * kerx + 2.0 * l1 * l3x * kerx + l1 * l3 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e2_ax_1(double x, double y)
    {
     return -(gradleg_tri_p10_e2_ax_0(x, y));
    }

    static double gradleg_tri_p10_e2_ay_0(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi9(l1 - l3); kerx = phi9x(l1 - l3) * (l1x - l3x);
     kery = phi9x(l1 - l3) * (l1y - l3y); kerxy = phi9xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e2_ay_1(double x, double y)
    {
     return -(gradleg_tri_p10_e2_ay_0(x, y));
    }

    static double gradleg_tri_p10_e2_bx_0(double x, double y)
    {
     double l1, l1x, l1y, l3, l3x, l3y;
     double ker, kerx, kery, kerxy;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
      l3 = lambda3(x, y); l3x = lambda3x(x, y); l3y = lambda3y(x, y);
     ker = phi9(l1 - l3); kerx = phi9x(l1 - l3) * (l1x - l3x);
     kery = phi9x(l1 - l3) * (l1y - l3y); kerxy = phi9xx(l1 - l3) * (l1x - l3x) * (l1y - l3y);
     return  (l1x * l3y * ker + l1y * l3x * ker + l1x * l3 * kery + l1 * l3x * kery + l1y * l3 * kerx + l1 * l3y * kerx +  l1 * l3 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e2_bx_1(double x, double y)
    {
     return -(gradleg_tri_p10_e2_bx_0(x, y));
    }

    static double gradleg_tri_p10_e2_by_0(double x, double y)
    {
     double l1, l1y, l3, l3y;
     double ker, kery, keryy;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     ker = phi9(l1 - l3); kery = phi9x(l1 - l3) * (l1y - l3y); keryy = phi9xx(l1 - l3) * sqr(l1y - l3y);
     return  (2.0 * l1y * l3y * ker + 2.0 * l1y * l3 * kery + 2.0 * l1 * l3y * kery + l1 * l3 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e2_by_1(double x, double y)
    {
     return -(gradleg_tri_p10_e2_by_0(x, y));
    }

     /* EDGE 3 */

    static double gradleg_tri_p10_e3_a_0(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return  (l2x * l1 * phi9(l2 - l1) + l2 * l1x * phi9(l2 - l1) + l2 * l1 * phi9x(l2 - l1) * (l2x - l1x)) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e3_a_1(double x, double y)
    {
     return -(gradleg_tri_p10_e3_a_0(x, y));
    }

    static double gradleg_tri_p10_e3_b_0(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return  (l2y * l1 * phi9(l2 - l1) + l2 * l1y * phi9(l2 - l1) + l2 * l1 * phi9x(l2 - l1) * (l2y - l1y)) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e3_b_1(double x, double y)
    {
     return -(gradleg_tri_p10_e3_b_0(x, y));
    }

    static double gradleg_tri_p10_e3_ax_0(double x, double y)
    {
     double l2, l2x, l1, l1x;
     double ker, kerx, kerxx;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     ker = phi9(l2 - l1); kerx = phi9x(l2 - l1) * (l2x - l1x); kerxx = phi9xx(l2 - l1) * sqr(l2x - l1x);
     return  (2.0 * l2x * l1x * ker + 2.0 * l2x * l1 * kerx + 2.0 * l2 * l1x * kerx + l2 * l1 * kerxx) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e3_ax_1(double x, double y)
    {
     return -(gradleg_tri_p10_e3_ax_0(x, y));
    }

    static double gradleg_tri_p10_e3_ay_0(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi9(l2 - l1); kerx = phi9x(l2 - l1) * (l2x - l1x);
     kery = phi9x(l2 - l1) * (l2y - l1y); kerxy = phi9xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e3_ay_1(double x, double y)
    {
     return -(gradleg_tri_p10_e3_ay_0(x, y));
    }

    static double gradleg_tri_p10_e3_bx_0(double x, double y)
    {
     double l2, l2x, l2y, l1, l1x, l1y;
     double ker, kerx, kery, kerxy;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l2y = lambda2y(x, y);
      l1 = lambda1(x, y); l1x = lambda1x(x, y); l1y = lambda1y(x, y);
     ker = phi9(l2 - l1); kerx = phi9x(l2 - l1) * (l2x - l1x);
     kery = phi9x(l2 - l1) * (l2y - l1y); kerxy = phi9xx(l2 - l1) * (l2x - l1x) * (l2y - l1y);
     return  (l2x * l1y * ker + l2y * l1x * ker + l2x * l1 * kery + l2 * l1x * kery + l2y * l1 * kerx + l2 * l1y * kerx +  l2 * l1 * kerxy) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e3_bx_1(double x, double y)
    {
     return -(gradleg_tri_p10_e3_bx_0(x, y));
    }

    static double gradleg_tri_p10_e3_by_0(double x, double y)
    {
     double l2, l2y, l1, l1y;
     double ker, kery, keryy;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     ker = phi9(l2 - l1); kery = phi9x(l2 - l1) * (l2y - l1y); keryy = phi9xx(l2 - l1) * sqr(l2y - l1y);
     return  (2.0 * l2y * l1y * ker + 2.0 * l2y * l1 * kery + 2.0 * l2 * l1y * kery + l2 * l1 * keryy) / 1.0000000000000;
    }

    static double gradleg_tri_p10_e3_by_1(double x, double y)
    {
     return -(gradleg_tri_p10_e3_by_0(x, y));
    }

    /* BUBBLE */

    /* Edge-based BUBBLE - order 10 */

     // EDGE 1
    static double gradleg_tri_p10_b1_a(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n11 * (l3 * l2 * Legendre8(l3 - l2));
    }

    static double gradleg_tri_p10_b1_b(double x, double y)
    {
     double l3, l2;
     l3 = lambda3(x, y); l2 = lambda2(x, y);
     return n12 * (l3 * l2 * Legendre8(l3 - l2));
    }

    static double gradleg_tri_p10_b1_ax(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n11 * (l3x * l2 * Legendre8(l3 - l2) + l3 * l2x * Legendre8(l3 - l2) + l3 * l2 * Legendre8x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p10_b1_bx(double x, double y)
    {
     double l3, l3x, l2, l2x;
     l3 = lambda3(x, y); l3x = lambda3x(x, y); l2 = lambda2(x, y); l2x = lambda2x(x, y);
     return n12 * (l3x * l2 * Legendre8(l3 - l2) + l3 * l2x * Legendre8(l3 - l2) + l3 * l2 * Legendre8x(l3 - l2) * (l3x - l2x));
    }

    static double gradleg_tri_p10_b1_ay(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n11 * (l3y * l2 * Legendre8(l3 - l2) + l3 * l2y * Legendre8(l3 - l2) + l3 * l2 * Legendre8x(l3 - l2) * (l3y - l2y));
    }

    static double gradleg_tri_p10_b1_by(double x, double y)
    {
     double l3, l3y, l2, l2y;
     l3 = lambda3(x, y); l3y = lambda3y(x, y); l2 = lambda2(x, y); l2y = lambda2y(x, y);
     return n12 * (l3y * l2 * Legendre8(l3 - l2) + l3 * l2y * Legendre8(l3 - l2) + l3 * l2 * Legendre8x(l3 - l2) * (l3y - l2y));
    }

     // EDGE 2
    static double gradleg_tri_p10_b2_a(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n21 * (l1 * l3 * Legendre8(l1 - l3));
    }

    static double gradleg_tri_p10_b2_b(double x, double y)
    {
     double l1, l3;
     l1 = lambda1(x, y); l3 = lambda3(x, y);
     return n22 * (l1 * l3 * Legendre8(l1 - l3));
    }

    static double gradleg_tri_p10_b2_ax(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n21 * (l1x * l3 * Legendre8(l1 - l3) + l1 * l3x * Legendre8(l1 - l3) + l1 * l3 * Legendre8x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p10_b2_bx(double x, double y)
    {
     double l1, l1x, l3, l3x;
     l1 = lambda1(x, y); l1x = lambda1x(x, y); l3 = lambda3(x, y); l3x = lambda3x(x, y);
     return n22 * (l1x * l3 * Legendre8(l1 - l3) + l1 * l3x * Legendre8(l1 - l3) + l1 * l3 * Legendre8x(l1 - l3) * (l1x - l3x));
    }

    static double gradleg_tri_p10_b2_ay(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n21 * (l1y * l3 * Legendre8(l1 - l3) + l1 * l3y * Legendre8(l1 - l3) + l1 * l3 * Legendre8x(l1 - l3) * (l1y - l3y));
    }

    static double gradleg_tri_p10_b2_by(double x, double y)
    {
     double l1, l1y, l3, l3y;
     l1 = lambda1(x, y); l1y = lambda1y(x, y); l3 = lambda3(x, y); l3y = lambda3y(x, y);
     return n22 * (l1y * l3 * Legendre8(l1 - l3) + l1 * l3y * Legendre8(l1 - l3) + l1 * l3 * Legendre8x(l1 - l3) * (l1y - l3y));
    }

     // EDGE 3
    static double gradleg_tri_p10_b3_a(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n31 * (l2 * l1 * Legendre8(l2 - l1));
    }

    static double gradleg_tri_p10_b3_b(double x, double y)
    {
     double l2, l1;
     l2 = lambda2(x, y); l1 = lambda1(x, y);
     return n32 * (l2 * l1 * Legendre8(l2 - l1));
    }

    static double gradleg_tri_p10_b3_ax(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n31 * (l2x * l1 * Legendre8(l2 - l1) + l2 * l1x * Legendre8(l2 - l1) + l2 * l1 * Legendre8x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p10_b3_bx(double x, double y)
    {
     double l2, l2x, l1, l1x;
     l2 = lambda2(x, y); l2x = lambda2x(x, y); l1 = lambda1(x, y); l1x = lambda1x(x, y);
     return n32 * (l2x * l1 * Legendre8(l2 - l1) + l2 * l1x * Legendre8(l2 - l1) + l2 * l1 * Legendre8x(l2 - l1) * (l2x - l1x));
    }

    static double gradleg_tri_p10_b3_ay(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n31 * (l2y * l1 * Legendre8(l2 - l1) + l2 * l1y * Legendre8(l2 - l1) + l2 * l1 * Legendre8x(l2 - l1) * (l2y - l1y));
    }

    static double gradleg_tri_p10_b3_by(double x, double y)
    {
     double l2, l2y, l1, l1y;
     l2 = lambda2(x, y); l2y = lambda2y(x, y); l1 = lambda1(x, y); l1y = lambda1y(x, y);
     return n32 * (l2y * l1 * Legendre8(l2 - l1) + l2 * l1y * Legendre8(l2 - l1) + l2 * l1 * Legendre8x(l2 - l1) * (l2y - l1y));
    }

    /* Genuine BUBBLE - order 10 */

    static double gradleg_tri_b1_b8_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre7(l2 - l1);
    }

    static double gradleg_tri_b1_b8_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b8_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre7(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre7x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b8_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b8_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre7(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre7x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b8_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b8_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre7(l2 - l1);
    }

    static double gradleg_tri_b1_b8_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b8_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre7(l2 - l1);
     L1x = Legendre0x(l3 - l2) * (l3x - l2x); L2x = Legendre7x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b1_b8_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b1_b8_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre0(l3 - l2); L2 = Legendre7(l2 - l1);
     L1y = Legendre0x(l3 - l2) * (l3y - l2y); L2y = Legendre7x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b1_b8_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b7_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre6(l2 - l1);
    }

    static double gradleg_tri_b2_b7_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b7_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre6(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre6x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b7_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b7_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre6(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre6x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b7_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b7_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre6(l2 - l1);
    }

    static double gradleg_tri_b2_b7_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b7_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre6(l2 - l1);
     L1x = Legendre1x(l3 - l2) * (l3x - l2x); L2x = Legendre6x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b2_b7_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b2_b7_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre1(l3 - l2); L2 = Legendre6(l2 - l1);
     L1y = Legendre1x(l3 - l2) * (l3y - l2y); L2y = Legendre6x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b2_b7_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b6_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre5(l2 - l1);
    }

    static double gradleg_tri_b3_b6_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b6_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre5(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre5x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b6_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b6_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre5(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre5x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b6_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b6_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre5(l2 - l1);
    }

    static double gradleg_tri_b3_b6_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b6_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre5(l2 - l1);
     L1x = Legendre2x(l3 - l2) * (l3x - l2x); L2x = Legendre5x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b3_b6_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b3_b6_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre2(l3 - l2); L2 = Legendre5(l2 - l1);
     L1y = Legendre2x(l3 - l2) * (l3y - l2y); L2y = Legendre5x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b3_b6_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b5_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre4(l2 - l1);
    }

    static double gradleg_tri_b4_b5_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b5_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre4(l2 - l1);
     L1x = Legendre3x(l3 - l2) * (l3x - l2x); L2x = Legendre4x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b4_b5_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b5_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre4(l2 - l1);
     L1y = Legendre3x(l3 - l2) * (l3y - l2y); L2y = Legendre4x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b4_b5_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b5_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre4(l2 - l1);
    }

    static double gradleg_tri_b4_b5_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b5_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre4(l2 - l1);
     L1x = Legendre3x(l3 - l2) * (l3x - l2x); L2x = Legendre4x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b4_b5_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b4_b5_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre3(l3 - l2); L2 = Legendre4(l2 - l1);
     L1y = Legendre3x(l3 - l2) * (l3y - l2y); L2y = Legendre4x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b4_b5_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b4_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre3(l2 - l1);
    }

    static double gradleg_tri_b5_b4_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b4_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre3(l2 - l1);
     L1x = Legendre4x(l3 - l2) * (l3x - l2x); L2x = Legendre3x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b5_b4_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b4_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre3(l2 - l1);
     L1y = Legendre4x(l3 - l2) * (l3y - l2y); L2y = Legendre3x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b5_b4_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b4_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre3(l2 - l1);
    }

    static double gradleg_tri_b5_b4_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b4_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre3(l2 - l1);
     L1x = Legendre4x(l3 - l2) * (l3x - l2x); L2x = Legendre3x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b5_b4_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b5_b4_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre4(l3 - l2); L2 = Legendre3(l2 - l1);
     L1y = Legendre4x(l3 - l2) * (l3y - l2y); L2y = Legendre3x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b5_b4_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b3_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b6_b3_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b3_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre5x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b6_b3_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b3_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre5x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b6_b3_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b3_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre2(l2 - l1);
    }

    static double gradleg_tri_b6_b3_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b3_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre2(l2 - l1);
     L1x = Legendre5x(l3 - l2) * (l3x - l2x); L2x = Legendre2x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b6_b3_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b6_b3_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre5(l3 - l2); L2 = Legendre2(l2 - l1);
     L1y = Legendre5x(l3 - l2) * (l3y - l2y); L2y = Legendre2x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b6_b3_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b2_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre6(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b7_b2_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b2_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre6(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre6x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b7_b2_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b2_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre6(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre6x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b7_b2_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b2_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre6(l3 - l2) * Legendre1(l2 - l1);
    }

    static double gradleg_tri_b7_b2_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b2_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre6(l3 - l2); L2 = Legendre1(l2 - l1);
     L1x = Legendre6x(l3 - l2) * (l3x - l2x); L2x = Legendre1x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b7_b2_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b7_b2_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre6(l3 - l2); L2 = Legendre1(l2 - l1);
     L1y = Legendre6x(l3 - l2) * (l3y - l2y); L2y = Legendre1x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b7_b2_2_ay(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b8_b1_1_a(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre7(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b8_b1_1_b(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b8_b1_1_ax(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre7(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre7x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b8_b1_1_bx(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b8_b1_1_ay(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre7(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre7x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b8_b1_1_by(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b8_b1_2_b(double x, double y)
    {
     double l1, l2, l3;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     return l1 * l2 * l3 * Legendre7(l3 - l2) * Legendre0(l2 - l1);
    }

    static double gradleg_tri_b8_b1_2_a(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b8_b1_2_bx(double x, double y)
    {
     double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
     L1 = Legendre7(l3 - l2); L2 = Legendre0(l2 - l1);
     L1x = Legendre7x(l3 - l2) * (l3x - l2x); L2x = Legendre0x(l2 - l1) * (l2x - l1x);
     return l1x * l2 *  l3 *  L1 *  L2 +          l1 * l2x * l3 *  L1 *  L2 +          l1 * l2 *  l3x * L1 *  L2 +          l1 * l2 *  l3 *  L1x * L2 +          l1 * l2 *  l3 *  L1 *  L2x;
    }

    static double gradleg_tri_b8_b1_2_ax(double x, double y)
    {
     return 0.0;
    }

    static double gradleg_tri_b8_b1_2_by(double x, double y)
    {
     double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;
     l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
     l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
     L1 = Legendre7(l3 - l2); L2 = Legendre0(l2 - l1);
     L1y = Legendre7x(l3 - l2) * (l3y - l2y); L2y = Legendre0x(l2 - l1) * (l2y - l1y);
     return l1y * l2 *  l3 *  L1 *  L2 +          l1 * l2y * l3 *  L1 *  L2 +          l1 * l2 *  l3y * L1 *  L2 +          l1 * l2 *  l3 *  L1y * L2 +          l1 * l2 *  l3 *  L1 *  L2y;
    }

    static double gradleg_tri_b8_b1_2_ay(double x, double y)
    {
     return 0.0;
    }

    static Shapeset::shape_fn_t gradleg_tri_fn_a[] =
    {
     gradleg_tri_p0_e1_a_0, gradleg_tri_p0_e1_a_1, gradleg_tri_p0_e2_a_0,  gradleg_tri_p0_e2_a_1, gradleg_tri_p0_e3_a_0, gradleg_tri_p0_e3_a_1,
     gradleg_tri_p1_e1_a,   gradleg_tri_p1_e1_a,   gradleg_tri_p1_e2_a,    gradleg_tri_p1_e2_a,   gradleg_tri_p1_e3_a,   gradleg_tri_p1_e3_a,
     gradleg_tri_p2_e1_a_0, gradleg_tri_p2_e1_a_1, gradleg_tri_p2_e2_a_0,  gradleg_tri_p2_e2_a_1, gradleg_tri_p2_e3_a_0, gradleg_tri_p2_e3_a_1,
     gradleg_tri_p3_e1_a,   gradleg_tri_p3_e1_a,   gradleg_tri_p3_e2_a,    gradleg_tri_p3_e2_a,   gradleg_tri_p3_e3_a,   gradleg_tri_p3_e3_a,
     gradleg_tri_p4_e1_a_0, gradleg_tri_p4_e1_a_1, gradleg_tri_p4_e2_a_0,  gradleg_tri_p4_e2_a_1, gradleg_tri_p4_e3_a_0, gradleg_tri_p4_e3_a_1,
     gradleg_tri_p5_e1_a,   gradleg_tri_p5_e1_a,   gradleg_tri_p5_e2_a,    gradleg_tri_p5_e2_a,   gradleg_tri_p5_e3_a,   gradleg_tri_p5_e3_a,
     gradleg_tri_p6_e1_a_0, gradleg_tri_p6_e1_a_1, gradleg_tri_p6_e2_a_0,  gradleg_tri_p6_e2_a_1, gradleg_tri_p6_e3_a_0, gradleg_tri_p6_e3_a_1,
     gradleg_tri_p7_e1_a,   gradleg_tri_p7_e1_a,   gradleg_tri_p7_e2_a,    gradleg_tri_p7_e2_a,   gradleg_tri_p7_e3_a,   gradleg_tri_p7_e3_a,
     gradleg_tri_p8_e1_a_0, gradleg_tri_p8_e1_a_1, gradleg_tri_p8_e2_a_0,  gradleg_tri_p8_e2_a_1, gradleg_tri_p8_e3_a_0, gradleg_tri_p8_e3_a_1,
     gradleg_tri_p9_e1_a,   gradleg_tri_p9_e1_a,   gradleg_tri_p9_e2_a,    gradleg_tri_p9_e2_a,   gradleg_tri_p9_e3_a,   gradleg_tri_p9_e3_a,
     gradleg_tri_p10_e1_a_0, gradleg_tri_p10_e1_a_1, gradleg_tri_p10_e2_a_0,  gradleg_tri_p10_e2_a_1, gradleg_tri_p10_e3_a_0, gradleg_tri_p10_e3_a_1,

      gradleg_tri_p2_b1_a,   gradleg_tri_p2_b2_a,   gradleg_tri_p2_b3_a,   gradleg_tri_p3_b1_a,   gradleg_tri_p3_b2_a,   gradleg_tri_p3_b3_a,   gradleg_tri_p4_b1_a,   gradleg_tri_p4_b2_a,   gradleg_tri_p4_b3_a,   gradleg_tri_p5_b1_a,   gradleg_tri_p5_b2_a,   gradleg_tri_p5_b3_a,   gradleg_tri_p6_b1_a,   gradleg_tri_p6_b2_a,   gradleg_tri_p6_b3_a,   gradleg_tri_p7_b1_a,   gradleg_tri_p7_b2_a,   gradleg_tri_p7_b3_a,   gradleg_tri_p8_b1_a,   gradleg_tri_p8_b2_a,   gradleg_tri_p8_b3_a,   gradleg_tri_p9_b1_a,   gradleg_tri_p9_b2_a,   gradleg_tri_p9_b3_a,   gradleg_tri_p10_b1_a,   gradleg_tri_p10_b2_a,   gradleg_tri_p10_b3_a,

      gradleg_tri_b1_b1_1_a, gradleg_tri_b1_b1_2_a,  gradleg_tri_b1_b2_1_a, gradleg_tri_b1_b2_2_a,  gradleg_tri_b2_b1_1_a, gradleg_tri_b2_b1_2_a,  gradleg_tri_b1_b3_1_a, gradleg_tri_b1_b3_2_a,  gradleg_tri_b2_b2_1_a, gradleg_tri_b2_b2_2_a,  gradleg_tri_b3_b1_1_a, gradleg_tri_b3_b1_2_a,  gradleg_tri_b1_b4_1_a, gradleg_tri_b1_b4_2_a,  gradleg_tri_b2_b3_1_a, gradleg_tri_b2_b3_2_a,  gradleg_tri_b3_b2_1_a, gradleg_tri_b3_b2_2_a,  gradleg_tri_b4_b1_1_a, gradleg_tri_b4_b1_2_a,  gradleg_tri_b1_b5_1_a, gradleg_tri_b1_b5_2_a,  gradleg_tri_b2_b4_1_a, gradleg_tri_b2_b4_2_a,  gradleg_tri_b3_b3_1_a, gradleg_tri_b3_b3_2_a,  gradleg_tri_b4_b2_1_a, gradleg_tri_b4_b2_2_a,  gradleg_tri_b5_b1_1_a, gradleg_tri_b5_b1_2_a,  gradleg_tri_b1_b6_1_a, gradleg_tri_b1_b6_2_a,  gradleg_tri_b2_b5_1_a, gradleg_tri_b2_b5_2_a,  gradleg_tri_b3_b4_1_a, gradleg_tri_b3_b4_2_a,  gradleg_tri_b4_b3_1_a, gradleg_tri_b4_b3_2_a,  gradleg_tri_b5_b2_1_a, gradleg_tri_b5_b2_2_a,  gradleg_tri_b6_b1_1_a, gradleg_tri_b6_b1_2_a,  gradleg_tri_b1_b7_1_a, gradleg_tri_b1_b7_2_a,  gradleg_tri_b2_b6_1_a, gradleg_tri_b2_b6_2_a,  gradleg_tri_b3_b5_1_a, gradleg_tri_b3_b5_2_a,  gradleg_tri_b4_b4_1_a, gradleg_tri_b4_b4_2_a,  gradleg_tri_b5_b3_1_a, gradleg_tri_b5_b3_2_a,  gradleg_tri_b6_b2_1_a, gradleg_tri_b6_b2_2_a,  gradleg_tri_b7_b1_1_a, gradleg_tri_b7_b1_2_a,  gradleg_tri_b1_b8_1_a, gradleg_tri_b1_b8_2_a,  gradleg_tri_b2_b7_1_a, gradleg_tri_b2_b7_2_a,  gradleg_tri_b3_b6_1_a, gradleg_tri_b3_b6_2_a,  gradleg_tri_b4_b5_1_a, gradleg_tri_b4_b5_2_a,  gradleg_tri_b5_b4_1_a, gradleg_tri_b5_b4_2_a,  gradleg_tri_b6_b3_1_a, gradleg_tri_b6_b3_2_a,  gradleg_tri_b7_b2_1_a, gradleg_tri_b7_b2_2_a,  gradleg_tri_b8_b1_1_a, gradleg_tri_b8_b1_2_a,
    };

    static Shapeset::shape_fn_t gradleg_tri_fn_b[] =
    {
     gradleg_tri_p0_e1_b_0, gradleg_tri_p0_e1_b_1, gradleg_tri_p0_e2_b_0,  gradleg_tri_p0_e2_b_1, gradleg_tri_p0_e3_b_0, gradleg_tri_p0_e3_b_1,
     gradleg_tri_p1_e1_b,   gradleg_tri_p1_e1_b,   gradleg_tri_p1_e2_b,    gradleg_tri_p1_e2_b,   gradleg_tri_p1_e3_b,   gradleg_tri_p1_e3_b,
     gradleg_tri_p2_e1_b_0, gradleg_tri_p2_e1_b_1, gradleg_tri_p2_e2_b_0,  gradleg_tri_p2_e2_b_1, gradleg_tri_p2_e3_b_0, gradleg_tri_p2_e3_b_1,
     gradleg_tri_p3_e1_b,   gradleg_tri_p3_e1_b,   gradleg_tri_p3_e2_b,    gradleg_tri_p3_e2_b,   gradleg_tri_p3_e3_b,   gradleg_tri_p3_e3_b,
     gradleg_tri_p4_e1_b_0, gradleg_tri_p4_e1_b_1, gradleg_tri_p4_e2_b_0,  gradleg_tri_p4_e2_b_1, gradleg_tri_p4_e3_b_0, gradleg_tri_p4_e3_b_1,
     gradleg_tri_p5_e1_b,   gradleg_tri_p5_e1_b,   gradleg_tri_p5_e2_b,    gradleg_tri_p5_e2_b,   gradleg_tri_p5_e3_b,   gradleg_tri_p5_e3_b,
     gradleg_tri_p6_e1_b_0, gradleg_tri_p6_e1_b_1, gradleg_tri_p6_e2_b_0,  gradleg_tri_p6_e2_b_1, gradleg_tri_p6_e3_b_0, gradleg_tri_p6_e3_b_1,
     gradleg_tri_p7_e1_b,   gradleg_tri_p7_e1_b,   gradleg_tri_p7_e2_b,    gradleg_tri_p7_e2_b,   gradleg_tri_p7_e3_b,   gradleg_tri_p7_e3_b,
     gradleg_tri_p8_e1_b_0, gradleg_tri_p8_e1_b_1, gradleg_tri_p8_e2_b_0,  gradleg_tri_p8_e2_b_1, gradleg_tri_p8_e3_b_0, gradleg_tri_p8_e3_b_1,
     gradleg_tri_p9_e1_b,   gradleg_tri_p9_e1_b,   gradleg_tri_p9_e2_b,    gradleg_tri_p9_e2_b,   gradleg_tri_p9_e3_b,   gradleg_tri_p9_e3_b,
     gradleg_tri_p10_e1_b_0, gradleg_tri_p10_e1_b_1, gradleg_tri_p10_e2_b_0,  gradleg_tri_p10_e2_b_1, gradleg_tri_p10_e3_b_0, gradleg_tri_p10_e3_b_1,

      gradleg_tri_p2_b1_b,   gradleg_tri_p2_b2_b,   gradleg_tri_p2_b3_b,   gradleg_tri_p3_b1_b,   gradleg_tri_p3_b2_b,   gradleg_tri_p3_b3_b,   gradleg_tri_p4_b1_b,   gradleg_tri_p4_b2_b,   gradleg_tri_p4_b3_b,   gradleg_tri_p5_b1_b,   gradleg_tri_p5_b2_b,   gradleg_tri_p5_b3_b,   gradleg_tri_p6_b1_b,   gradleg_tri_p6_b2_b,   gradleg_tri_p6_b3_b,   gradleg_tri_p7_b1_b,   gradleg_tri_p7_b2_b,   gradleg_tri_p7_b3_b,   gradleg_tri_p8_b1_b,   gradleg_tri_p8_b2_b,   gradleg_tri_p8_b3_b,   gradleg_tri_p9_b1_b,   gradleg_tri_p9_b2_b,   gradleg_tri_p9_b3_b,   gradleg_tri_p10_b1_b,   gradleg_tri_p10_b2_b,   gradleg_tri_p10_b3_b,

      gradleg_tri_b1_b1_1_b, gradleg_tri_b1_b1_2_b,  gradleg_tri_b1_b2_1_b, gradleg_tri_b1_b2_2_b,  gradleg_tri_b2_b1_1_b, gradleg_tri_b2_b1_2_b,  gradleg_tri_b1_b3_1_b, gradleg_tri_b1_b3_2_b,  gradleg_tri_b2_b2_1_b, gradleg_tri_b2_b2_2_b,  gradleg_tri_b3_b1_1_b, gradleg_tri_b3_b1_2_b,  gradleg_tri_b1_b4_1_b, gradleg_tri_b1_b4_2_b,  gradleg_tri_b2_b3_1_b, gradleg_tri_b2_b3_2_b,  gradleg_tri_b3_b2_1_b, gradleg_tri_b3_b2_2_b,  gradleg_tri_b4_b1_1_b, gradleg_tri_b4_b1_2_b,  gradleg_tri_b1_b5_1_b, gradleg_tri_b1_b5_2_b,  gradleg_tri_b2_b4_1_b, gradleg_tri_b2_b4_2_b,  gradleg_tri_b3_b3_1_b, gradleg_tri_b3_b3_2_b,  gradleg_tri_b4_b2_1_b, gradleg_tri_b4_b2_2_b,  gradleg_tri_b5_b1_1_b, gradleg_tri_b5_b1_2_b,  gradleg_tri_b1_b6_1_b, gradleg_tri_b1_b6_2_b,  gradleg_tri_b2_b5_1_b, gradleg_tri_b2_b5_2_b,  gradleg_tri_b3_b4_1_b, gradleg_tri_b3_b4_2_b,  gradleg_tri_b4_b3_1_b, gradleg_tri_b4_b3_2_b,  gradleg_tri_b5_b2_1_b, gradleg_tri_b5_b2_2_b,  gradleg_tri_b6_b1_1_b, gradleg_tri_b6_b1_2_b,  gradleg_tri_b1_b7_1_b, gradleg_tri_b1_b7_2_b,  gradleg_tri_b2_b6_1_b, gradleg_tri_b2_b6_2_b,  gradleg_tri_b3_b5_1_b, gradleg_tri_b3_b5_2_b,  gradleg_tri_b4_b4_1_b, gradleg_tri_b4_b4_2_b,  gradleg_tri_b5_b3_1_b, gradleg_tri_b5_b3_2_b,  gradleg_tri_b6_b2_1_b, gradleg_tri_b6_b2_2_b,  gradleg_tri_b7_b1_1_b, gradleg_tri_b7_b1_2_b,  gradleg_tri_b1_b8_1_b, gradleg_tri_b1_b8_2_b,  gradleg_tri_b2_b7_1_b, gradleg_tri_b2_b7_2_b,  gradleg_tri_b3_b6_1_b, gradleg_tri_b3_b6_2_b,  gradleg_tri_b4_b5_1_b, gradleg_tri_b4_b5_2_b,  gradleg_tri_b5_b4_1_b, gradleg_tri_b5_b4_2_b,  gradleg_tri_b6_b3_1_b, gradleg_tri_b6_b3_2_b,  gradleg_tri_b7_b2_1_b, gradleg_tri_b7_b2_2_b,  gradleg_tri_b8_b1_1_b, gradleg_tri_b8_b1_2_b,
    };

    static Shapeset::shape_fn_t gradleg_tri_fn_ax[] =
    {
     gradleg_tri_p0_e1_ax_0, gradleg_tri_p0_e1_ax_1, gradleg_tri_p0_e2_ax_0,  gradleg_tri_p0_e2_ax_1, gradleg_tri_p0_e3_ax_0, gradleg_tri_p0_e3_ax_1,
     gradleg_tri_p1_e1_ax,   gradleg_tri_p1_e1_ax,   gradleg_tri_p1_e2_ax,    gradleg_tri_p1_e2_ax,   gradleg_tri_p1_e3_ax,   gradleg_tri_p1_e3_ax,
     gradleg_tri_p2_e1_ax_0, gradleg_tri_p2_e1_ax_1, gradleg_tri_p2_e2_ax_0,  gradleg_tri_p2_e2_ax_1, gradleg_tri_p2_e3_ax_0, gradleg_tri_p2_e3_ax_1,
     gradleg_tri_p3_e1_ax,   gradleg_tri_p3_e1_ax,   gradleg_tri_p3_e2_ax,    gradleg_tri_p3_e2_ax,   gradleg_tri_p3_e3_ax,   gradleg_tri_p3_e3_ax,
     gradleg_tri_p4_e1_ax_0, gradleg_tri_p4_e1_ax_1, gradleg_tri_p4_e2_ax_0,  gradleg_tri_p4_e2_ax_1, gradleg_tri_p4_e3_ax_0, gradleg_tri_p4_e3_ax_1,
     gradleg_tri_p5_e1_ax,   gradleg_tri_p5_e1_ax,   gradleg_tri_p5_e2_ax,    gradleg_tri_p5_e2_ax,   gradleg_tri_p5_e3_ax,   gradleg_tri_p5_e3_ax,
     gradleg_tri_p6_e1_ax_0, gradleg_tri_p6_e1_ax_1, gradleg_tri_p6_e2_ax_0,  gradleg_tri_p6_e2_ax_1, gradleg_tri_p6_e3_ax_0, gradleg_tri_p6_e3_ax_1,
     gradleg_tri_p7_e1_ax,   gradleg_tri_p7_e1_ax,   gradleg_tri_p7_e2_ax,    gradleg_tri_p7_e2_ax,   gradleg_tri_p7_e3_ax,   gradleg_tri_p7_e3_ax,
     gradleg_tri_p8_e1_ax_0, gradleg_tri_p8_e1_ax_1, gradleg_tri_p8_e2_ax_0,  gradleg_tri_p8_e2_ax_1, gradleg_tri_p8_e3_ax_0, gradleg_tri_p8_e3_ax_1,
     gradleg_tri_p9_e1_ax,   gradleg_tri_p9_e1_ax,   gradleg_tri_p9_e2_ax,    gradleg_tri_p9_e2_ax,   gradleg_tri_p9_e3_ax,   gradleg_tri_p9_e3_ax,
     gradleg_tri_p10_e1_ax_0, gradleg_tri_p10_e1_ax_1, gradleg_tri_p10_e2_ax_0,  gradleg_tri_p10_e2_ax_1, gradleg_tri_p10_e3_ax_0, gradleg_tri_p10_e3_ax_1,

      gradleg_tri_p2_b1_ax,   gradleg_tri_p2_b2_ax,   gradleg_tri_p2_b3_ax,   gradleg_tri_p3_b1_ax,   gradleg_tri_p3_b2_ax,   gradleg_tri_p3_b3_ax,   gradleg_tri_p4_b1_ax,   gradleg_tri_p4_b2_ax,   gradleg_tri_p4_b3_ax,   gradleg_tri_p5_b1_ax,   gradleg_tri_p5_b2_ax,   gradleg_tri_p5_b3_ax,   gradleg_tri_p6_b1_ax,   gradleg_tri_p6_b2_ax,   gradleg_tri_p6_b3_ax,   gradleg_tri_p7_b1_ax,   gradleg_tri_p7_b2_ax,   gradleg_tri_p7_b3_ax,   gradleg_tri_p8_b1_ax,   gradleg_tri_p8_b2_ax,   gradleg_tri_p8_b3_ax,   gradleg_tri_p9_b1_ax,   gradleg_tri_p9_b2_ax,   gradleg_tri_p9_b3_ax,   gradleg_tri_p10_b1_ax,   gradleg_tri_p10_b2_ax,   gradleg_tri_p10_b3_ax,

      gradleg_tri_b1_b1_1_ax, gradleg_tri_b1_b1_2_ax,  gradleg_tri_b1_b2_1_ax, gradleg_tri_b1_b2_2_ax,  gradleg_tri_b2_b1_1_ax, gradleg_tri_b2_b1_2_ax,  gradleg_tri_b1_b3_1_ax, gradleg_tri_b1_b3_2_ax,  gradleg_tri_b2_b2_1_ax, gradleg_tri_b2_b2_2_ax,  gradleg_tri_b3_b1_1_ax, gradleg_tri_b3_b1_2_ax,  gradleg_tri_b1_b4_1_ax, gradleg_tri_b1_b4_2_ax,  gradleg_tri_b2_b3_1_ax, gradleg_tri_b2_b3_2_ax,  gradleg_tri_b3_b2_1_ax, gradleg_tri_b3_b2_2_ax,  gradleg_tri_b4_b1_1_ax, gradleg_tri_b4_b1_2_ax,  gradleg_tri_b1_b5_1_ax, gradleg_tri_b1_b5_2_ax,  gradleg_tri_b2_b4_1_ax, gradleg_tri_b2_b4_2_ax,  gradleg_tri_b3_b3_1_ax, gradleg_tri_b3_b3_2_ax,  gradleg_tri_b4_b2_1_ax, gradleg_tri_b4_b2_2_ax,  gradleg_tri_b5_b1_1_ax, gradleg_tri_b5_b1_2_ax,  gradleg_tri_b1_b6_1_ax, gradleg_tri_b1_b6_2_ax,  gradleg_tri_b2_b5_1_ax, gradleg_tri_b2_b5_2_ax,  gradleg_tri_b3_b4_1_ax, gradleg_tri_b3_b4_2_ax,  gradleg_tri_b4_b3_1_ax, gradleg_tri_b4_b3_2_ax,  gradleg_tri_b5_b2_1_ax, gradleg_tri_b5_b2_2_ax,  gradleg_tri_b6_b1_1_ax, gradleg_tri_b6_b1_2_ax,  gradleg_tri_b1_b7_1_ax, gradleg_tri_b1_b7_2_ax,  gradleg_tri_b2_b6_1_ax, gradleg_tri_b2_b6_2_ax,  gradleg_tri_b3_b5_1_ax, gradleg_tri_b3_b5_2_ax,  gradleg_tri_b4_b4_1_ax, gradleg_tri_b4_b4_2_ax,  gradleg_tri_b5_b3_1_ax, gradleg_tri_b5_b3_2_ax,  gradleg_tri_b6_b2_1_ax, gradleg_tri_b6_b2_2_ax,  gradleg_tri_b7_b1_1_ax, gradleg_tri_b7_b1_2_ax,  gradleg_tri_b1_b8_1_ax, gradleg_tri_b1_b8_2_ax,  gradleg_tri_b2_b7_1_ax, gradleg_tri_b2_b7_2_ax,  gradleg_tri_b3_b6_1_ax, gradleg_tri_b3_b6_2_ax,  gradleg_tri_b4_b5_1_ax, gradleg_tri_b4_b5_2_ax,  gradleg_tri_b5_b4_1_ax, gradleg_tri_b5_b4_2_ax,  gradleg_tri_b6_b3_1_ax, gradleg_tri_b6_b3_2_ax,  gradleg_tri_b7_b2_1_ax, gradleg_tri_b7_b2_2_ax,  gradleg_tri_b8_b1_1_ax, gradleg_tri_b8_b1_2_ax,
    };

    static Shapeset::shape_fn_t gradleg_tri_fn_bx[] =
    {
     gradleg_tri_p0_e1_bx_0, gradleg_tri_p0_e1_bx_1, gradleg_tri_p0_e2_bx_0,  gradleg_tri_p0_e2_bx_1, gradleg_tri_p0_e3_bx_0, gradleg_tri_p0_e3_bx_1,
     gradleg_tri_p1_e1_bx,   gradleg_tri_p1_e1_bx,   gradleg_tri_p1_e2_bx,    gradleg_tri_p1_e2_bx,   gradleg_tri_p1_e3_bx,   gradleg_tri_p1_e3_bx,
     gradleg_tri_p2_e1_bx_0, gradleg_tri_p2_e1_bx_1, gradleg_tri_p2_e2_bx_0,  gradleg_tri_p2_e2_bx_1, gradleg_tri_p2_e3_bx_0, gradleg_tri_p2_e3_bx_1,
     gradleg_tri_p3_e1_bx,   gradleg_tri_p3_e1_bx,   gradleg_tri_p3_e2_bx,    gradleg_tri_p3_e2_bx,   gradleg_tri_p3_e3_bx,   gradleg_tri_p3_e3_bx,
     gradleg_tri_p4_e1_bx_0, gradleg_tri_p4_e1_bx_1, gradleg_tri_p4_e2_bx_0,  gradleg_tri_p4_e2_bx_1, gradleg_tri_p4_e3_bx_0, gradleg_tri_p4_e3_bx_1,
     gradleg_tri_p5_e1_bx,   gradleg_tri_p5_e1_bx,   gradleg_tri_p5_e2_bx,    gradleg_tri_p5_e2_bx,   gradleg_tri_p5_e3_bx,   gradleg_tri_p5_e3_bx,
     gradleg_tri_p6_e1_bx_0, gradleg_tri_p6_e1_bx_1, gradleg_tri_p6_e2_bx_0,  gradleg_tri_p6_e2_bx_1, gradleg_tri_p6_e3_bx_0, gradleg_tri_p6_e3_bx_1,
     gradleg_tri_p7_e1_bx,   gradleg_tri_p7_e1_bx,   gradleg_tri_p7_e2_bx,    gradleg_tri_p7_e2_bx,   gradleg_tri_p7_e3_bx,   gradleg_tri_p7_e3_bx,
     gradleg_tri_p8_e1_bx_0, gradleg_tri_p8_e1_bx_1, gradleg_tri_p8_e2_bx_0,  gradleg_tri_p8_e2_bx_1, gradleg_tri_p8_e3_bx_0, gradleg_tri_p8_e3_bx_1,
     gradleg_tri_p9_e1_bx,   gradleg_tri_p9_e1_bx,   gradleg_tri_p9_e2_bx,    gradleg_tri_p9_e2_bx,   gradleg_tri_p9_e3_bx,   gradleg_tri_p9_e3_bx,
     gradleg_tri_p10_e1_bx_0, gradleg_tri_p10_e1_bx_1, gradleg_tri_p10_e2_bx_0,  gradleg_tri_p10_e2_bx_1, gradleg_tri_p10_e3_bx_0, gradleg_tri_p10_e3_bx_1,

      gradleg_tri_p2_b1_bx,   gradleg_tri_p2_b2_bx,   gradleg_tri_p2_b3_bx,   gradleg_tri_p3_b1_bx,   gradleg_tri_p3_b2_bx,   gradleg_tri_p3_b3_bx,   gradleg_tri_p4_b1_bx,   gradleg_tri_p4_b2_bx,   gradleg_tri_p4_b3_bx,   gradleg_tri_p5_b1_bx,   gradleg_tri_p5_b2_bx,   gradleg_tri_p5_b3_bx,   gradleg_tri_p6_b1_bx,   gradleg_tri_p6_b2_bx,   gradleg_tri_p6_b3_bx,   gradleg_tri_p7_b1_bx,   gradleg_tri_p7_b2_bx,   gradleg_tri_p7_b3_bx,   gradleg_tri_p8_b1_bx,   gradleg_tri_p8_b2_bx,   gradleg_tri_p8_b3_bx,   gradleg_tri_p9_b1_bx,   gradleg_tri_p9_b2_bx,   gradleg_tri_p9_b3_bx,   gradleg_tri_p10_b1_bx,   gradleg_tri_p10_b2_bx,   gradleg_tri_p10_b3_bx,

      gradleg_tri_b1_b1_1_bx, gradleg_tri_b1_b1_2_bx,  gradleg_tri_b1_b2_1_bx, gradleg_tri_b1_b2_2_bx,  gradleg_tri_b2_b1_1_bx, gradleg_tri_b2_b1_2_bx,  gradleg_tri_b1_b3_1_bx, gradleg_tri_b1_b3_2_bx,  gradleg_tri_b2_b2_1_bx, gradleg_tri_b2_b2_2_bx,  gradleg_tri_b3_b1_1_bx, gradleg_tri_b3_b1_2_bx,  gradleg_tri_b1_b4_1_bx, gradleg_tri_b1_b4_2_bx,  gradleg_tri_b2_b3_1_bx, gradleg_tri_b2_b3_2_bx,  gradleg_tri_b3_b2_1_bx, gradleg_tri_b3_b2_2_bx,  gradleg_tri_b4_b1_1_bx, gradleg_tri_b4_b1_2_bx,  gradleg_tri_b1_b5_1_bx, gradleg_tri_b1_b5_2_bx,  gradleg_tri_b2_b4_1_bx, gradleg_tri_b2_b4_2_bx,  gradleg_tri_b3_b3_1_bx, gradleg_tri_b3_b3_2_bx,  gradleg_tri_b4_b2_1_bx, gradleg_tri_b4_b2_2_bx,  gradleg_tri_b5_b1_1_bx, gradleg_tri_b5_b1_2_bx,  gradleg_tri_b1_b6_1_bx, gradleg_tri_b1_b6_2_bx,  gradleg_tri_b2_b5_1_bx, gradleg_tri_b2_b5_2_bx,  gradleg_tri_b3_b4_1_bx, gradleg_tri_b3_b4_2_bx,  gradleg_tri_b4_b3_1_bx, gradleg_tri_b4_b3_2_bx,  gradleg_tri_b5_b2_1_bx, gradleg_tri_b5_b2_2_bx,  gradleg_tri_b6_b1_1_bx, gradleg_tri_b6_b1_2_bx,  gradleg_tri_b1_b7_1_bx, gradleg_tri_b1_b7_2_bx,  gradleg_tri_b2_b6_1_bx, gradleg_tri_b2_b6_2_bx,  gradleg_tri_b3_b5_1_bx, gradleg_tri_b3_b5_2_bx,  gradleg_tri_b4_b4_1_bx, gradleg_tri_b4_b4_2_bx,  gradleg_tri_b5_b3_1_bx, gradleg_tri_b5_b3_2_bx,  gradleg_tri_b6_b2_1_bx, gradleg_tri_b6_b2_2_bx,  gradleg_tri_b7_b1_1_bx, gradleg_tri_b7_b1_2_bx,  gradleg_tri_b1_b8_1_bx, gradleg_tri_b1_b8_2_bx,  gradleg_tri_b2_b7_1_bx, gradleg_tri_b2_b7_2_bx,  gradleg_tri_b3_b6_1_bx, gradleg_tri_b3_b6_2_bx,  gradleg_tri_b4_b5_1_bx, gradleg_tri_b4_b5_2_bx,  gradleg_tri_b5_b4_1_bx, gradleg_tri_b5_b4_2_bx,  gradleg_tri_b6_b3_1_bx, gradleg_tri_b6_b3_2_bx,  gradleg_tri_b7_b2_1_bx, gradleg_tri_b7_b2_2_bx,  gradleg_tri_b8_b1_1_bx, gradleg_tri_b8_b1_2_bx,
    };

    static Shapeset::shape_fn_t gradleg_tri_fn_ay[] =
    {
     gradleg_tri_p0_e1_ay_0, gradleg_tri_p0_e1_ay_1, gradleg_tri_p0_e2_ay_0,  gradleg_tri_p0_e2_ay_1, gradleg_tri_p0_e3_ay_0, gradleg_tri_p0_e3_ay_1,
     gradleg_tri_p1_e1_ay,   gradleg_tri_p1_e1_ay,   gradleg_tri_p1_e2_ay,    gradleg_tri_p1_e2_ay,   gradleg_tri_p1_e3_ay,   gradleg_tri_p1_e3_ay,
     gradleg_tri_p2_e1_ay_0, gradleg_tri_p2_e1_ay_1, gradleg_tri_p2_e2_ay_0,  gradleg_tri_p2_e2_ay_1, gradleg_tri_p2_e3_ay_0, gradleg_tri_p2_e3_ay_1,
     gradleg_tri_p3_e1_ay,   gradleg_tri_p3_e1_ay,   gradleg_tri_p3_e2_ay,    gradleg_tri_p3_e2_ay,   gradleg_tri_p3_e3_ay,   gradleg_tri_p3_e3_ay,
     gradleg_tri_p4_e1_ay_0, gradleg_tri_p4_e1_ay_1, gradleg_tri_p4_e2_ay_0,  gradleg_tri_p4_e2_ay_1, gradleg_tri_p4_e3_ay_0, gradleg_tri_p4_e3_ay_1,
     gradleg_tri_p5_e1_ay,   gradleg_tri_p5_e1_ay,   gradleg_tri_p5_e2_ay,    gradleg_tri_p5_e2_ay,   gradleg_tri_p5_e3_ay,   gradleg_tri_p5_e3_ay,
     gradleg_tri_p6_e1_ay_0, gradleg_tri_p6_e1_ay_1, gradleg_tri_p6_e2_ay_0,  gradleg_tri_p6_e2_ay_1, gradleg_tri_p6_e3_ay_0, gradleg_tri_p6_e3_ay_1,
     gradleg_tri_p7_e1_ay,   gradleg_tri_p7_e1_ay,   gradleg_tri_p7_e2_ay,    gradleg_tri_p7_e2_ay,   gradleg_tri_p7_e3_ay,   gradleg_tri_p7_e3_ay,
     gradleg_tri_p8_e1_ay_0, gradleg_tri_p8_e1_ay_1, gradleg_tri_p8_e2_ay_0,  gradleg_tri_p8_e2_ay_1, gradleg_tri_p8_e3_ay_0, gradleg_tri_p8_e3_ay_1,
     gradleg_tri_p9_e1_ay,   gradleg_tri_p9_e1_ay,   gradleg_tri_p9_e2_ay,    gradleg_tri_p9_e2_ay,   gradleg_tri_p9_e3_ay,   gradleg_tri_p9_e3_ay,
     gradleg_tri_p10_e1_ay_0, gradleg_tri_p10_e1_ay_1, gradleg_tri_p10_e2_ay_0,  gradleg_tri_p10_e2_ay_1, gradleg_tri_p10_e3_ay_0, gradleg_tri_p10_e3_ay_1,

      gradleg_tri_p2_b1_ay,   gradleg_tri_p2_b2_ay,   gradleg_tri_p2_b3_ay,   gradleg_tri_p3_b1_ay,   gradleg_tri_p3_b2_ay,   gradleg_tri_p3_b3_ay,   gradleg_tri_p4_b1_ay,   gradleg_tri_p4_b2_ay,   gradleg_tri_p4_b3_ay,   gradleg_tri_p5_b1_ay,   gradleg_tri_p5_b2_ay,   gradleg_tri_p5_b3_ay,   gradleg_tri_p6_b1_ay,   gradleg_tri_p6_b2_ay,   gradleg_tri_p6_b3_ay,   gradleg_tri_p7_b1_ay,   gradleg_tri_p7_b2_ay,   gradleg_tri_p7_b3_ay,   gradleg_tri_p8_b1_ay,   gradleg_tri_p8_b2_ay,   gradleg_tri_p8_b3_ay,   gradleg_tri_p9_b1_ay,   gradleg_tri_p9_b2_ay,   gradleg_tri_p9_b3_ay,   gradleg_tri_p10_b1_ay,   gradleg_tri_p10_b2_ay,   gradleg_tri_p10_b3_ay,

      gradleg_tri_b1_b1_1_ay, gradleg_tri_b1_b1_2_ay,  gradleg_tri_b1_b2_1_ay, gradleg_tri_b1_b2_2_ay,  gradleg_tri_b2_b1_1_ay, gradleg_tri_b2_b1_2_ay,  gradleg_tri_b1_b3_1_ay, gradleg_tri_b1_b3_2_ay,  gradleg_tri_b2_b2_1_ay, gradleg_tri_b2_b2_2_ay,  gradleg_tri_b3_b1_1_ay, gradleg_tri_b3_b1_2_ay,  gradleg_tri_b1_b4_1_ay, gradleg_tri_b1_b4_2_ay,  gradleg_tri_b2_b3_1_ay, gradleg_tri_b2_b3_2_ay,  gradleg_tri_b3_b2_1_ay, gradleg_tri_b3_b2_2_ay,  gradleg_tri_b4_b1_1_ay, gradleg_tri_b4_b1_2_ay,  gradleg_tri_b1_b5_1_ay, gradleg_tri_b1_b5_2_ay,  gradleg_tri_b2_b4_1_ay, gradleg_tri_b2_b4_2_ay,  gradleg_tri_b3_b3_1_ay, gradleg_tri_b3_b3_2_ay,  gradleg_tri_b4_b2_1_ay, gradleg_tri_b4_b2_2_ay,  gradleg_tri_b5_b1_1_ay, gradleg_tri_b5_b1_2_ay,  gradleg_tri_b1_b6_1_ay, gradleg_tri_b1_b6_2_ay,  gradleg_tri_b2_b5_1_ay, gradleg_tri_b2_b5_2_ay,  gradleg_tri_b3_b4_1_ay, gradleg_tri_b3_b4_2_ay,  gradleg_tri_b4_b3_1_ay, gradleg_tri_b4_b3_2_ay,  gradleg_tri_b5_b2_1_ay, gradleg_tri_b5_b2_2_ay,  gradleg_tri_b6_b1_1_ay, gradleg_tri_b6_b1_2_ay,  gradleg_tri_b1_b7_1_ay, gradleg_tri_b1_b7_2_ay,  gradleg_tri_b2_b6_1_ay, gradleg_tri_b2_b6_2_ay,  gradleg_tri_b3_b5_1_ay, gradleg_tri_b3_b5_2_ay,  gradleg_tri_b4_b4_1_ay, gradleg_tri_b4_b4_2_ay,  gradleg_tri_b5_b3_1_ay, gradleg_tri_b5_b3_2_ay,  gradleg_tri_b6_b2_1_ay, gradleg_tri_b6_b2_2_ay,  gradleg_tri_b7_b1_1_ay, gradleg_tri_b7_b1_2_ay,  gradleg_tri_b1_b8_1_ay, gradleg_tri_b1_b8_2_ay,  gradleg_tri_b2_b7_1_ay, gradleg_tri_b2_b7_2_ay,  gradleg_tri_b3_b6_1_ay, gradleg_tri_b3_b6_2_ay,  gradleg_tri_b4_b5_1_ay, gradleg_tri_b4_b5_2_ay,  gradleg_tri_b5_b4_1_ay, gradleg_tri_b5_b4_2_ay,  gradleg_tri_b6_b3_1_ay, gradleg_tri_b6_b3_2_ay,  gradleg_tri_b7_b2_1_ay, gradleg_tri_b7_b2_2_ay,  gradleg_tri_b8_b1_1_ay, gradleg_tri_b8_b1_2_ay,
    };

    static Shapeset::shape_fn_t gradleg_tri_fn_by[] =
    {
     gradleg_tri_p0_e1_by_0, gradleg_tri_p0_e1_by_1, gradleg_tri_p0_e2_by_0,  gradleg_tri_p0_e2_by_1, gradleg_tri_p0_e3_by_0, gradleg_tri_p0_e3_by_1,
     gradleg_tri_p1_e1_by,   gradleg_tri_p1_e1_by,   gradleg_tri_p1_e2_by,    gradleg_tri_p1_e2_by,   gradleg_tri_p1_e3_by,   gradleg_tri_p1_e3_by,
     gradleg_tri_p2_e1_by_0, gradleg_tri_p2_e1_by_1, gradleg_tri_p2_e2_by_0,  gradleg_tri_p2_e2_by_1, gradleg_tri_p2_e3_by_0, gradleg_tri_p2_e3_by_1,
     gradleg_tri_p3_e1_by,   gradleg_tri_p3_e1_by,   gradleg_tri_p3_e2_by,    gradleg_tri_p3_e2_by,   gradleg_tri_p3_e3_by,   gradleg_tri_p3_e3_by,
     gradleg_tri_p4_e1_by_0, gradleg_tri_p4_e1_by_1, gradleg_tri_p4_e2_by_0,  gradleg_tri_p4_e2_by_1, gradleg_tri_p4_e3_by_0, gradleg_tri_p4_e3_by_1,
     gradleg_tri_p5_e1_by,   gradleg_tri_p5_e1_by,   gradleg_tri_p5_e2_by,    gradleg_tri_p5_e2_by,   gradleg_tri_p5_e3_by,   gradleg_tri_p5_e3_by,
     gradleg_tri_p6_e1_by_0, gradleg_tri_p6_e1_by_1, gradleg_tri_p6_e2_by_0,  gradleg_tri_p6_e2_by_1, gradleg_tri_p6_e3_by_0, gradleg_tri_p6_e3_by_1,
     gradleg_tri_p7_e1_by,   gradleg_tri_p7_e1_by,   gradleg_tri_p7_e2_by,    gradleg_tri_p7_e2_by,   gradleg_tri_p7_e3_by,   gradleg_tri_p7_e3_by,
     gradleg_tri_p8_e1_by_0, gradleg_tri_p8_e1_by_1, gradleg_tri_p8_e2_by_0,  gradleg_tri_p8_e2_by_1, gradleg_tri_p8_e3_by_0, gradleg_tri_p8_e3_by_1,
     gradleg_tri_p9_e1_by,   gradleg_tri_p9_e1_by,   gradleg_tri_p9_e2_by,    gradleg_tri_p9_e2_by,   gradleg_tri_p9_e3_by,   gradleg_tri_p9_e3_by,
     gradleg_tri_p10_e1_by_0, gradleg_tri_p10_e1_by_1, gradleg_tri_p10_e2_by_0,  gradleg_tri_p10_e2_by_1, gradleg_tri_p10_e3_by_0, gradleg_tri_p10_e3_by_1,

      gradleg_tri_p2_b1_by,   gradleg_tri_p2_b2_by,   gradleg_tri_p2_b3_by,   gradleg_tri_p3_b1_by,   gradleg_tri_p3_b2_by,   gradleg_tri_p3_b3_by,   gradleg_tri_p4_b1_by,   gradleg_tri_p4_b2_by,   gradleg_tri_p4_b3_by,   gradleg_tri_p5_b1_by,   gradleg_tri_p5_b2_by,   gradleg_tri_p5_b3_by,   gradleg_tri_p6_b1_by,   gradleg_tri_p6_b2_by,   gradleg_tri_p6_b3_by,   gradleg_tri_p7_b1_by,   gradleg_tri_p7_b2_by,   gradleg_tri_p7_b3_by,   gradleg_tri_p8_b1_by,   gradleg_tri_p8_b2_by,   gradleg_tri_p8_b3_by,   gradleg_tri_p9_b1_by,   gradleg_tri_p9_b2_by,   gradleg_tri_p9_b3_by,   gradleg_tri_p10_b1_by,   gradleg_tri_p10_b2_by,   gradleg_tri_p10_b3_by,

      gradleg_tri_b1_b1_1_by, gradleg_tri_b1_b1_2_by,  gradleg_tri_b1_b2_1_by, gradleg_tri_b1_b2_2_by,  gradleg_tri_b2_b1_1_by, gradleg_tri_b2_b1_2_by,  gradleg_tri_b1_b3_1_by, gradleg_tri_b1_b3_2_by,  gradleg_tri_b2_b2_1_by, gradleg_tri_b2_b2_2_by,  gradleg_tri_b3_b1_1_by, gradleg_tri_b3_b1_2_by,  gradleg_tri_b1_b4_1_by, gradleg_tri_b1_b4_2_by,  gradleg_tri_b2_b3_1_by, gradleg_tri_b2_b3_2_by,  gradleg_tri_b3_b2_1_by, gradleg_tri_b3_b2_2_by,  gradleg_tri_b4_b1_1_by, gradleg_tri_b4_b1_2_by,  gradleg_tri_b1_b5_1_by, gradleg_tri_b1_b5_2_by,  gradleg_tri_b2_b4_1_by, gradleg_tri_b2_b4_2_by,  gradleg_tri_b3_b3_1_by, gradleg_tri_b3_b3_2_by,  gradleg_tri_b4_b2_1_by, gradleg_tri_b4_b2_2_by,  gradleg_tri_b5_b1_1_by, gradleg_tri_b5_b1_2_by,  gradleg_tri_b1_b6_1_by, gradleg_tri_b1_b6_2_by,  gradleg_tri_b2_b5_1_by, gradleg_tri_b2_b5_2_by,  gradleg_tri_b3_b4_1_by, gradleg_tri_b3_b4_2_by,  gradleg_tri_b4_b3_1_by, gradleg_tri_b4_b3_2_by,  gradleg_tri_b5_b2_1_by, gradleg_tri_b5_b2_2_by,  gradleg_tri_b6_b1_1_by, gradleg_tri_b6_b1_2_by,  gradleg_tri_b1_b7_1_by, gradleg_tri_b1_b7_2_by,  gradleg_tri_b2_b6_1_by, gradleg_tri_b2_b6_2_by,  gradleg_tri_b3_b5_1_by, gradleg_tri_b3_b5_2_by,  gradleg_tri_b4_b4_1_by, gradleg_tri_b4_b4_2_by,  gradleg_tri_b5_b3_1_by, gradleg_tri_b5_b3_2_by,  gradleg_tri_b6_b2_1_by, gradleg_tri_b6_b2_2_by,  gradleg_tri_b7_b1_1_by, gradleg_tri_b7_b1_2_by,  gradleg_tri_b1_b8_1_by, gradleg_tri_b1_b8_2_by,  gradleg_tri_b2_b7_1_by, gradleg_tri_b2_b7_2_by,  gradleg_tri_b3_b6_1_by, gradleg_tri_b3_b6_2_by,  gradleg_tri_b4_b5_1_by, gradleg_tri_b4_b5_2_by,  gradleg_tri_b5_b4_1_by, gradleg_tri_b5_b4_2_by,  gradleg_tri_b6_b3_1_by, gradleg_tri_b6_b3_2_by,  gradleg_tri_b7_b2_1_by, gradleg_tri_b7_b2_2_by,  gradleg_tri_b8_b1_1_by, gradleg_tri_b8_b1_2_by,
    };

    static int gradleg_tri_bubble_indices_all_orders[] =
    {
     66, 67, 68,
     69, 70, 71, 93, 94,
     72, 73, 74, 95, 96, 97, 98,
     75, 76, 77, 99, 100, 101, 102, 103, 104,
     78, 79, 80, 105, 106, 107, 108, 109, 110, 111, 112,
     81, 82, 83, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
     84, 85, 86, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134,
     87, 88, 89, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148,
     90, 91, 92, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
    };

    static int* gradleg_tri_bubble_indices[11] =
    {
      NULL, NULL,
      gradleg_tri_bubble_indices_all_orders,
      gradleg_tri_bubble_indices_all_orders,
      gradleg_tri_bubble_indices_all_orders,
      gradleg_tri_bubble_indices_all_orders,
      gradleg_tri_bubble_indices_all_orders,
      gradleg_tri_bubble_indices_all_orders,
      gradleg_tri_bubble_indices_all_orders,
      gradleg_tri_bubble_indices_all_orders,
      gradleg_tri_bubble_indices_all_orders
    };

    static int gradleg_tri_bubble_count[11] = { 0, 0, 3, 8, 15, 24, 35, 48, 63, 80, 99, };

    static int gradleg_tri_edge_indices_0[22] =  {  0, 1, 6, 7, 12, 13, 18, 19, 24, 25, 30, 31, 36, 37, 42, 43, 48, 49, 54, 55, 60, 61, };
    static int gradleg_tri_edge_indices_1[22] =  {  2, 3, 8, 9, 14, 15, 20, 21, 26, 27, 32, 33, 38, 39, 44, 45, 50, 51, 56, 57, 62, 63, };
    static int gradleg_tri_edge_indices_2[22] =  {  4, 5, 10, 11, 16, 17, 22, 23, 28, 29, 34, 35, 40, 41, 46, 47, 52, 53, 58, 59, 64, 65, };

    static int* gradleg_tri_edge_indices[3] =
    {
      gradleg_tri_edge_indices_0,
      gradleg_tri_edge_indices_1,
      gradleg_tri_edge_indices_2,
    };

    static int gradleg_tri_vertex_indices[3] = { -1, -1, -1 };

    static int gradleg_tri_index_to_order[] =
    {
     0, 0, 0, 0, 0, 0,
     1, 1, 1, 1, 1, 1,
     2, 2, 2, 2, 2, 2,
     3, 3, 3, 3, 3, 3,
     4, 4, 4, 4, 4, 4,
     5, 5, 5, 5, 5, 5,
     6, 6, 6, 6, 6, 6,
     7, 7, 7, 7, 7, 7,
     8, 8, 8, 8, 8, 8,
     9, 9, 9, 9, 9, 9,
     10, 10, 10, 10, 10, 10,
     2, 2, 2,
     3, 3, 3,
     4, 4, 4,
     5, 5, 5,
     6, 6, 6,
     7, 7, 7,
     8, 8, 8,
     9, 9, 9,
     10, 10, 10,
     3, 3,
     4, 4, 4, 4,
     5, 5, 5, 5, 5, 5,
     6, 6, 6, 6, 6, 6, 6, 6,
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
     9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
     10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    };

    static Shapeset::shape_fn_t* gradleg_tri_shape_fn_table[2] =
    {
      gradleg_tri_fn_a,
      gradleg_tri_fn_b
    };

    static Shapeset::shape_fn_t* gradleg_tri_shape_fn_table_x[2] =
    {
      gradleg_tri_fn_ax,
      gradleg_tri_fn_bx
    };

    static Shapeset::shape_fn_t* gradleg_tri_shape_fn_table_y[2] =
    {
      gradleg_tri_fn_ay,
      gradleg_tri_fn_by
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // QUADS

    /* Whitney fns - constant tangential component */

    static double gradleg_quad_p0_e1_a_0(double x, double y)
    {
      return Legendre0(x) * l0(y);
    }

    static double gradleg_quad_p0_e1_a_1(double x, double y)
    {
      return -(Legendre0(x) * l0(y));
    }

    static double gradleg_quad_p0_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e1_ax_0(double x, double y)
    {
      return Legendre0x(x) * l0(y);
    }

    static double gradleg_quad_p0_e1_ax_1(double x, double y)
    {
      return -(Legendre0x(x) * l0(y));
    }

    static double gradleg_quad_p0_e1_ay_0(double x, double y)
    {
      return Legendre0(x) * dl0(y);
    }

    static double gradleg_quad_p0_e1_ay_1(double x, double y)
    {
      return -(Legendre0(x) * dl0(y));
    }

    static double gradleg_quad_p0_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e2_b_0(double x, double y)
    {
      return l1(x) * Legendre0(y);
    }

    static double gradleg_quad_p0_e2_b_1(double x, double y)
    {
      return -(l1(x) * Legendre0(y));
    }

    static double gradleg_quad_p0_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e2_bx_0(double x, double y)
    {
      return dl1(x) * Legendre0(y);
    }

    static double gradleg_quad_p0_e2_bx_1(double x, double y)
    {
      return -(dl1(x) * Legendre0(y));
    }

    static double gradleg_quad_p0_e2_by_0(double x, double y)
    {
      return l1(x) * Legendre0x(y);
    }

    static double gradleg_quad_p0_e2_by_1(double x, double y)
    {
      return -(l1(x) * Legendre0x(y));
    }

    static double gradleg_quad_p0_e3_a_1(double x, double y)
    {
      return Legendre0(x) * l1(y);
    }

    static double gradleg_quad_p0_e3_a_0(double x, double y)
    {
      return -(Legendre0(x) * l1(y));
    }

    static double gradleg_quad_p0_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e3_ax_1(double x, double y)
    {
      return Legendre0x(x) * l1(y);
    }

    static double gradleg_quad_p0_e3_ax_0(double x, double y)
    {
      return -(Legendre0x(x) * l1(y));
    }

    static double gradleg_quad_p0_e3_ay_1(double x, double y)
    {
      return Legendre0(x) * dl1(y);
    }

    static double gradleg_quad_p0_e3_ay_0(double x, double y)
    {
      return -(Legendre0(x) * dl1(y));
    }

    static double gradleg_quad_p0_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e4_b_1(double x, double y)
    {
      return l0(x) * Legendre0(y);
    }

    static double gradleg_quad_p0_e4_b_0(double x, double y)
    {
      return -(l0(x) * Legendre0(y));
    }

    static double gradleg_quad_p0_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e4_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0_e4_bx_1(double x, double y)
    {
      return dl0(x) * Legendre0(y);
    }

    static double gradleg_quad_p0_e4_bx_0(double x, double y)
    {
      return -(dl0(x) * Legendre0(y));
    }

    static double gradleg_quad_p0_e4_by_1(double x, double y)
    {
      return l0(x) * Legendre0x(y);
    }

    static double gradleg_quad_p0_e4_by_0(double x, double y)
    {
      return -(l0(x) * Legendre0x(y));
    }

    /* Edge fns - gradients of Scalar lobatto edge functions */

    static double gradleg_quad_l0_l2_a(double x, double y)
    {
     return dl0(x) * l2(y);
    }

    static double gradleg_quad_l0_l2_b(double x, double y)
    {
     return l0(x) * dl2(y);
    }

    static double gradleg_quad_l0_l2_ax(double x, double y)
    {
     return d2l0(x) * l2(y);
    }

    static double gradleg_quad_l0_l2_bx(double x, double y)
    {
     return dl0(x) * dl2(y);
    }

    static double gradleg_quad_l0_l2_ay(double x, double y)
    {
     return dl0(x) * dl2(y);
    }

    static double gradleg_quad_l0_l2_by(double x, double y)
    {
     return l0(x) * d2l2(y);
    }

    static double gradleg_quad_l0_l3_a_0(double x, double y)
    {
     return -dl0(x) * l3(y);
    }

    static double gradleg_quad_l0_l3_a_1(double x, double y)
    {
     return  dl0(x) * l3(y);
    }

    static double gradleg_quad_l0_l3_b_0(double x, double y)
    {
     return -l0(x) * dl3(y);
    }

    static double gradleg_quad_l0_l3_b_1(double x, double y)
    {
     return  l0(x) * dl3(y);
    }

    static double gradleg_quad_l0_l3_ax_0(double x, double y)
    {
     return -d2l0(x) * l3(y);
    }

    static double gradleg_quad_l0_l3_ax_1(double x, double y)
    {
     return  d2l0(x) * l3(y);
    }

    static double gradleg_quad_l0_l3_bx_0(double x, double y)
    {
     return -dl0(x) * dl3(y);
    }

    static double gradleg_quad_l0_l3_bx_1(double x, double y)
    {
     return  dl0(x) * dl3(y);
    }

    static double gradleg_quad_l0_l3_ay_0(double x, double y)
    {
     return -dl0(x) * dl3(y);
    }

    static double gradleg_quad_l0_l3_ay_1(double x, double y)
    {
     return  dl0(x) * dl3(y);
    }

    static double gradleg_quad_l0_l3_by_0(double x, double y)
    {
     return -l0(x) * d2l3(y);
    }

    static double gradleg_quad_l0_l3_by_1(double x, double y)
    {
     return  l0(x) * d2l3(y);
    }

    static double gradleg_quad_l0_l4_a(double x, double y)
    {
     return dl0(x) * l4(y);
    }

    static double gradleg_quad_l0_l4_b(double x, double y)
    {
     return l0(x) * dl4(y);
    }

    static double gradleg_quad_l0_l4_ax(double x, double y)
    {
     return d2l0(x) * l4(y);
    }

    static double gradleg_quad_l0_l4_bx(double x, double y)
    {
     return dl0(x) * dl4(y);
    }

    static double gradleg_quad_l0_l4_ay(double x, double y)
    {
     return dl0(x) * dl4(y);
    }

    static double gradleg_quad_l0_l4_by(double x, double y)
    {
     return l0(x) * d2l4(y);
    }

    static double gradleg_quad_l0_l5_a_0(double x, double y)
    {
     return -dl0(x) * l5(y);
    }

    static double gradleg_quad_l0_l5_a_1(double x, double y)
    {
     return  dl0(x) * l5(y);
    }

    static double gradleg_quad_l0_l5_b_0(double x, double y)
    {
     return -l0(x) * dl5(y);
    }

    static double gradleg_quad_l0_l5_b_1(double x, double y)
    {
     return  l0(x) * dl5(y);
    }

    static double gradleg_quad_l0_l5_ax_0(double x, double y)
    {
     return -d2l0(x) * l5(y);
    }

    static double gradleg_quad_l0_l5_ax_1(double x, double y)
    {
     return  d2l0(x) * l5(y);
    }

    static double gradleg_quad_l0_l5_bx_0(double x, double y)
    {
     return -dl0(x) * dl5(y);
    }

    static double gradleg_quad_l0_l5_bx_1(double x, double y)
    {
     return  dl0(x) * dl5(y);
    }

    static double gradleg_quad_l0_l5_ay_0(double x, double y)
    {
     return -dl0(x) * dl5(y);
    }

    static double gradleg_quad_l0_l5_ay_1(double x, double y)
    {
     return  dl0(x) * dl5(y);
    }

    static double gradleg_quad_l0_l5_by_0(double x, double y)
    {
     return -l0(x) * d2l5(y);
    }

    static double gradleg_quad_l0_l5_by_1(double x, double y)
    {
     return  l0(x) * d2l5(y);
    }

    static double gradleg_quad_l0_l6_a(double x, double y)
    {
     return dl0(x) * l6(y);
    }

    static double gradleg_quad_l0_l6_b(double x, double y)
    {
     return l0(x) * dl6(y);
    }

    static double gradleg_quad_l0_l6_ax(double x, double y)
    {
     return d2l0(x) * l6(y);
    }

    static double gradleg_quad_l0_l6_bx(double x, double y)
    {
     return dl0(x) * dl6(y);
    }

    static double gradleg_quad_l0_l6_ay(double x, double y)
    {
     return dl0(x) * dl6(y);
    }

    static double gradleg_quad_l0_l6_by(double x, double y)
    {
     return l0(x) * d2l6(y);
    }

    static double gradleg_quad_l0_l7_a_0(double x, double y)
    {
     return -dl0(x) * l7(y);
    }

    static double gradleg_quad_l0_l7_a_1(double x, double y)
    {
     return  dl0(x) * l7(y);
    }

    static double gradleg_quad_l0_l7_b_0(double x, double y)
    {
     return -l0(x) * dl7(y);
    }

    static double gradleg_quad_l0_l7_b_1(double x, double y)
    {
     return  l0(x) * dl7(y);
    }

    static double gradleg_quad_l0_l7_ax_0(double x, double y)
    {
     return -d2l0(x) * l7(y);
    }

    static double gradleg_quad_l0_l7_ax_1(double x, double y)
    {
     return  d2l0(x) * l7(y);
    }

    static double gradleg_quad_l0_l7_bx_0(double x, double y)
    {
     return -dl0(x) * dl7(y);
    }

    static double gradleg_quad_l0_l7_bx_1(double x, double y)
    {
     return  dl0(x) * dl7(y);
    }

    static double gradleg_quad_l0_l7_ay_0(double x, double y)
    {
     return -dl0(x) * dl7(y);
    }

    static double gradleg_quad_l0_l7_ay_1(double x, double y)
    {
     return  dl0(x) * dl7(y);
    }

    static double gradleg_quad_l0_l7_by_0(double x, double y)
    {
     return -l0(x) * d2l7(y);
    }

    static double gradleg_quad_l0_l7_by_1(double x, double y)
    {
     return  l0(x) * d2l7(y);
    }

    static double gradleg_quad_l0_l8_a(double x, double y)
    {
     return dl0(x) * l8(y);
    }

    static double gradleg_quad_l0_l8_b(double x, double y)
    {
     return l0(x) * dl8(y);
    }

    static double gradleg_quad_l0_l8_ax(double x, double y)
    {
     return d2l0(x) * l8(y);
    }

    static double gradleg_quad_l0_l8_bx(double x, double y)
    {
     return dl0(x) * dl8(y);
    }

    static double gradleg_quad_l0_l8_ay(double x, double y)
    {
     return dl0(x) * dl8(y);
    }

    static double gradleg_quad_l0_l8_by(double x, double y)
    {
     return l0(x) * d2l8(y);
    }

    static double gradleg_quad_l0_l9_a_0(double x, double y)
    {
     return -dl0(x) * l9(y);
    }

    static double gradleg_quad_l0_l9_a_1(double x, double y)
    {
     return  dl0(x) * l9(y);
    }

    static double gradleg_quad_l0_l9_b_0(double x, double y)
    {
     return -l0(x) * dl9(y);
    }

    static double gradleg_quad_l0_l9_b_1(double x, double y)
    {
     return  l0(x) * dl9(y);
    }

    static double gradleg_quad_l0_l9_ax_0(double x, double y)
    {
     return -d2l0(x) * l9(y);
    }

    static double gradleg_quad_l0_l9_ax_1(double x, double y)
    {
     return  d2l0(x) * l9(y);
    }

    static double gradleg_quad_l0_l9_bx_0(double x, double y)
    {
     return -dl0(x) * dl9(y);
    }

    static double gradleg_quad_l0_l9_bx_1(double x, double y)
    {
     return  dl0(x) * dl9(y);
    }

    static double gradleg_quad_l0_l9_ay_0(double x, double y)
    {
     return -dl0(x) * dl9(y);
    }

    static double gradleg_quad_l0_l9_ay_1(double x, double y)
    {
     return  dl0(x) * dl9(y);
    }

    static double gradleg_quad_l0_l9_by_0(double x, double y)
    {
     return -l0(x) * d2l9(y);
    }

    static double gradleg_quad_l0_l9_by_1(double x, double y)
    {
     return  l0(x) * d2l9(y);
    }

    static double gradleg_quad_l0_l10_a(double x, double y)
    {
     return dl0(x) * l10(y);
    }

    static double gradleg_quad_l0_l10_b(double x, double y)
    {
     return l0(x) * dl10(y);
    }

    static double gradleg_quad_l0_l10_ax(double x, double y)
    {
     return d2l0(x) * l10(y);
    }

    static double gradleg_quad_l0_l10_bx(double x, double y)
    {
     return dl0(x) * dl10(y);
    }

    static double gradleg_quad_l0_l10_ay(double x, double y)
    {
     return dl0(x) * dl10(y);
    }

    static double gradleg_quad_l0_l10_by(double x, double y)
    {
     return l0(x) * d2l10(y);
    }

    static double gradleg_quad_l0_l11_a_0(double x, double y)
    {
     return -dl0(x) * l11(y);
    }

    static double gradleg_quad_l0_l11_a_1(double x, double y)
    {
     return  dl0(x) * l11(y);
    }

    static double gradleg_quad_l0_l11_b_0(double x, double y)
    {
     return -l0(x) * dl11(y);
    }

    static double gradleg_quad_l0_l11_b_1(double x, double y)
    {
     return  l0(x) * dl11(y);
    }

    static double gradleg_quad_l0_l11_ax_0(double x, double y)
    {
     return -d2l0(x) * l11(y);
    }

    static double gradleg_quad_l0_l11_ax_1(double x, double y)
    {
     return  d2l0(x) * l11(y);
    }

    static double gradleg_quad_l0_l11_bx_0(double x, double y)
    {
     return -dl0(x) * dl11(y);
    }

    static double gradleg_quad_l0_l11_bx_1(double x, double y)
    {
     return  dl0(x) * dl11(y);
    }

    static double gradleg_quad_l0_l11_ay_0(double x, double y)
    {
     return -dl0(x) * dl11(y);
    }

    static double gradleg_quad_l0_l11_ay_1(double x, double y)
    {
     return  dl0(x) * dl11(y);
    }

    static double gradleg_quad_l0_l11_by_0(double x, double y)
    {
     return -l0(x) * d2l11(y);
    }

    static double gradleg_quad_l0_l11_by_1(double x, double y)
    {
     return  l0(x) * d2l11(y);
    }

    static double gradleg_quad_l1_l2_a(double x, double y)
    {
     return dl1(x) * l2(y);
    }

    static double gradleg_quad_l1_l2_b(double x, double y)
    {
     return l1(x) * dl2(y);
    }

    static double gradleg_quad_l1_l2_ax(double x, double y)
    {
     return d2l1(x) * l2(y);
    }

    static double gradleg_quad_l1_l2_bx(double x, double y)
    {
     return dl1(x) * dl2(y);
    }

    static double gradleg_quad_l1_l2_ay(double x, double y)
    {
     return dl1(x) * dl2(y);
    }

    static double gradleg_quad_l1_l2_by(double x, double y)
    {
     return l1(x) * d2l2(y);
    }

    static double gradleg_quad_l1_l3_a_0(double x, double y)
    {
     return  dl1(x) * l3(y);
    }

    static double gradleg_quad_l1_l3_a_1(double x, double y)
    {
     return -dl1(x) * l3(y);
    }

    static double gradleg_quad_l1_l3_b_0(double x, double y)
    {
     return  l1(x) * dl3(y);
    }

    static double gradleg_quad_l1_l3_b_1(double x, double y)
    {
     return -l1(x) * dl3(y);
    }

    static double gradleg_quad_l1_l3_ax_0(double x, double y)
    {
     return  d2l1(x) * l3(y);
    }

    static double gradleg_quad_l1_l3_ax_1(double x, double y)
    {
     return -d2l1(x) * l3(y);
    }

    static double gradleg_quad_l1_l3_bx_0(double x, double y)
    {
     return  dl1(x) * dl3(y);
    }

    static double gradleg_quad_l1_l3_bx_1(double x, double y)
    {
     return -dl1(x) * dl3(y);
    }

    static double gradleg_quad_l1_l3_ay_0(double x, double y)
    {
     return  dl1(x) * dl3(y);
    }

    static double gradleg_quad_l1_l3_ay_1(double x, double y)
    {
     return -dl1(x) * dl3(y);
    }

    static double gradleg_quad_l1_l3_by_0(double x, double y)
    {
     return  l1(x) * d2l3(y);
    }

    static double gradleg_quad_l1_l3_by_1(double x, double y)
    {
     return -l1(x) * d2l3(y);
    }

    static double gradleg_quad_l1_l4_a(double x, double y)
    {
     return dl1(x) * l4(y);
    }

    static double gradleg_quad_l1_l4_b(double x, double y)
    {
     return l1(x) * dl4(y);
    }

    static double gradleg_quad_l1_l4_ax(double x, double y)
    {
     return d2l1(x) * l4(y);
    }

    static double gradleg_quad_l1_l4_bx(double x, double y)
    {
     return dl1(x) * dl4(y);
    }

    static double gradleg_quad_l1_l4_ay(double x, double y)
    {
     return dl1(x) * dl4(y);
    }

    static double gradleg_quad_l1_l4_by(double x, double y)
    {
     return l1(x) * d2l4(y);
    }

    static double gradleg_quad_l1_l5_a_0(double x, double y)
    {
     return  dl1(x) * l5(y);
    }

    static double gradleg_quad_l1_l5_a_1(double x, double y)
    {
     return -dl1(x) * l5(y);
    }

    static double gradleg_quad_l1_l5_b_0(double x, double y)
    {
     return  l1(x) * dl5(y);
    }

    static double gradleg_quad_l1_l5_b_1(double x, double y)
    {
     return -l1(x) * dl5(y);
    }

    static double gradleg_quad_l1_l5_ax_0(double x, double y)
    {
     return  d2l1(x) * l5(y);
    }

    static double gradleg_quad_l1_l5_ax_1(double x, double y)
    {
     return -d2l1(x) * l5(y);
    }

    static double gradleg_quad_l1_l5_bx_0(double x, double y)
    {
     return  dl1(x) * dl5(y);
    }

    static double gradleg_quad_l1_l5_bx_1(double x, double y)
    {
     return -dl1(x) * dl5(y);
    }

    static double gradleg_quad_l1_l5_ay_0(double x, double y)
    {
     return  dl1(x) * dl5(y);
    }

    static double gradleg_quad_l1_l5_ay_1(double x, double y)
    {
     return -dl1(x) * dl5(y);
    }

    static double gradleg_quad_l1_l5_by_0(double x, double y)
    {
     return  l1(x) * d2l5(y);
    }

    static double gradleg_quad_l1_l5_by_1(double x, double y)
    {
     return -l1(x) * d2l5(y);
    }

    static double gradleg_quad_l1_l6_a(double x, double y)
    {
     return dl1(x) * l6(y);
    }

    static double gradleg_quad_l1_l6_b(double x, double y)
    {
     return l1(x) * dl6(y);
    }

    static double gradleg_quad_l1_l6_ax(double x, double y)
    {
     return d2l1(x) * l6(y);
    }

    static double gradleg_quad_l1_l6_bx(double x, double y)
    {
     return dl1(x) * dl6(y);
    }

    static double gradleg_quad_l1_l6_ay(double x, double y)
    {
     return dl1(x) * dl6(y);
    }

    static double gradleg_quad_l1_l6_by(double x, double y)
    {
     return l1(x) * d2l6(y);
    }

    static double gradleg_quad_l1_l7_a_0(double x, double y)
    {
     return  dl1(x) * l7(y);
    }

    static double gradleg_quad_l1_l7_a_1(double x, double y)
    {
     return -dl1(x) * l7(y);
    }

    static double gradleg_quad_l1_l7_b_0(double x, double y)
    {
     return  l1(x) * dl7(y);
    }

    static double gradleg_quad_l1_l7_b_1(double x, double y)
    {
     return -l1(x) * dl7(y);
    }

    static double gradleg_quad_l1_l7_ax_0(double x, double y)
    {
     return  d2l1(x) * l7(y);
    }

    static double gradleg_quad_l1_l7_ax_1(double x, double y)
    {
     return -d2l1(x) * l7(y);
    }

    static double gradleg_quad_l1_l7_bx_0(double x, double y)
    {
     return  dl1(x) * dl7(y);
    }

    static double gradleg_quad_l1_l7_bx_1(double x, double y)
    {
     return -dl1(x) * dl7(y);
    }

    static double gradleg_quad_l1_l7_ay_0(double x, double y)
    {
     return  dl1(x) * dl7(y);
    }

    static double gradleg_quad_l1_l7_ay_1(double x, double y)
    {
     return -dl1(x) * dl7(y);
    }

    static double gradleg_quad_l1_l7_by_0(double x, double y)
    {
     return  l1(x) * d2l7(y);
    }

    static double gradleg_quad_l1_l7_by_1(double x, double y)
    {
     return -l1(x) * d2l7(y);
    }

    static double gradleg_quad_l1_l8_a(double x, double y)
    {
     return dl1(x) * l8(y);
    }

    static double gradleg_quad_l1_l8_b(double x, double y)
    {
     return l1(x) * dl8(y);
    }

    static double gradleg_quad_l1_l8_ax(double x, double y)
    {
     return d2l1(x) * l8(y);
    }

    static double gradleg_quad_l1_l8_bx(double x, double y)
    {
     return dl1(x) * dl8(y);
    }

    static double gradleg_quad_l1_l8_ay(double x, double y)
    {
     return dl1(x) * dl8(y);
    }

    static double gradleg_quad_l1_l8_by(double x, double y)
    {
     return l1(x) * d2l8(y);
    }

    static double gradleg_quad_l1_l9_a_0(double x, double y)
    {
     return  dl1(x) * l9(y);
    }

    static double gradleg_quad_l1_l9_a_1(double x, double y)
    {
     return -dl1(x) * l9(y);
    }

    static double gradleg_quad_l1_l9_b_0(double x, double y)
    {
     return  l1(x) * dl9(y);
    }

    static double gradleg_quad_l1_l9_b_1(double x, double y)
    {
     return -l1(x) * dl9(y);
    }

    static double gradleg_quad_l1_l9_ax_0(double x, double y)
    {
     return  d2l1(x) * l9(y);
    }

    static double gradleg_quad_l1_l9_ax_1(double x, double y)
    {
     return -d2l1(x) * l9(y);
    }

    static double gradleg_quad_l1_l9_bx_0(double x, double y)
    {
     return  dl1(x) * dl9(y);
    }

    static double gradleg_quad_l1_l9_bx_1(double x, double y)
    {
     return -dl1(x) * dl9(y);
    }

    static double gradleg_quad_l1_l9_ay_0(double x, double y)
    {
     return  dl1(x) * dl9(y);
    }

    static double gradleg_quad_l1_l9_ay_1(double x, double y)
    {
     return -dl1(x) * dl9(y);
    }

    static double gradleg_quad_l1_l9_by_0(double x, double y)
    {
     return  l1(x) * d2l9(y);
    }

    static double gradleg_quad_l1_l9_by_1(double x, double y)
    {
     return -l1(x) * d2l9(y);
    }

    static double gradleg_quad_l1_l10_a(double x, double y)
    {
     return dl1(x) * l10(y);
    }

    static double gradleg_quad_l1_l10_b(double x, double y)
    {
     return l1(x) * dl10(y);
    }

    static double gradleg_quad_l1_l10_ax(double x, double y)
    {
     return d2l1(x) * l10(y);
    }

    static double gradleg_quad_l1_l10_bx(double x, double y)
    {
     return dl1(x) * dl10(y);
    }

    static double gradleg_quad_l1_l10_ay(double x, double y)
    {
     return dl1(x) * dl10(y);
    }

    static double gradleg_quad_l1_l10_by(double x, double y)
    {
     return l1(x) * d2l10(y);
    }

    static double gradleg_quad_l1_l11_a_0(double x, double y)
    {
     return  dl1(x) * l11(y);
    }

    static double gradleg_quad_l1_l11_a_1(double x, double y)
    {
     return -dl1(x) * l11(y);
    }

    static double gradleg_quad_l1_l11_b_0(double x, double y)
    {
     return  l1(x) * dl11(y);
    }

    static double gradleg_quad_l1_l11_b_1(double x, double y)
    {
     return -l1(x) * dl11(y);
    }

    static double gradleg_quad_l1_l11_ax_0(double x, double y)
    {
     return  d2l1(x) * l11(y);
    }

    static double gradleg_quad_l1_l11_ax_1(double x, double y)
    {
     return -d2l1(x) * l11(y);
    }

    static double gradleg_quad_l1_l11_bx_0(double x, double y)
    {
     return  dl1(x) * dl11(y);
    }

    static double gradleg_quad_l1_l11_bx_1(double x, double y)
    {
     return -dl1(x) * dl11(y);
    }

    static double gradleg_quad_l1_l11_ay_0(double x, double y)
    {
     return  dl1(x) * dl11(y);
    }

    static double gradleg_quad_l1_l11_ay_1(double x, double y)
    {
     return -dl1(x) * dl11(y);
    }

    static double gradleg_quad_l1_l11_by_0(double x, double y)
    {
     return  l1(x) * d2l11(y);
    }

    static double gradleg_quad_l1_l11_by_1(double x, double y)
    {
     return -l1(x) * d2l11(y);
    }

    static double gradleg_quad_l2_l0_a(double x, double y)
    {
     return dl2(x) * l0(y);
    }

    static double gradleg_quad_l2_l0_b(double x, double y)
    {
     return l2(x) * dl0(y);
    }

    static double gradleg_quad_l2_l0_ax(double x, double y)
    {
     return d2l2(x) * l0(y);
    }

    static double gradleg_quad_l2_l0_bx(double x, double y)
    {
     return dl2(x) * dl0(y);
    }

    static double gradleg_quad_l2_l0_ay(double x, double y)
    {
     return dl2(x) * dl0(y);
    }

    static double gradleg_quad_l2_l0_by(double x, double y)
    {
     return l2(x) * d2l0(y);
    }

    static double gradleg_quad_l3_l0_a_0(double x, double y)
    {
     return  dl3(x) * l0(y);
    }

    static double gradleg_quad_l3_l0_a_1(double x, double y)
    {
     return -dl3(x) * l0(y);
    }

    static double gradleg_quad_l3_l0_b_0(double x, double y)
    {
     return  l3(x) * dl0(y);
    }

    static double gradleg_quad_l3_l0_b_1(double x, double y)
    {
     return -l3(x) * dl0(y);
    }

    static double gradleg_quad_l3_l0_ax_0(double x, double y)
    {
     return  d2l3(x) * l0(y);
    }

    static double gradleg_quad_l3_l0_ax_1(double x, double y)
    {
     return -d2l3(x) * l0(y);
    }

    static double gradleg_quad_l3_l0_bx_0(double x, double y)
    {
     return  dl3(x) * dl0(y);
    }

    static double gradleg_quad_l3_l0_bx_1(double x, double y)
    {
     return -dl3(x) * dl0(y);
    }

    static double gradleg_quad_l3_l0_ay_0(double x, double y)
    {
     return  dl3(x) * dl0(y);
    }

    static double gradleg_quad_l3_l0_ay_1(double x, double y)
    {
     return -dl3(x) * dl0(y);
    }

    static double gradleg_quad_l3_l0_by_0(double x, double y)
    {
     return  l3(x) * d2l0(y);
    }

    static double gradleg_quad_l3_l0_by_1(double x, double y)
    {
     return -l3(x) * d2l0(y);
    }

    static double gradleg_quad_l4_l0_a(double x, double y)
    {
     return dl4(x) * l0(y);
    }

    static double gradleg_quad_l4_l0_b(double x, double y)
    {
     return l4(x) * dl0(y);
    }

    static double gradleg_quad_l4_l0_ax(double x, double y)
    {
     return d2l4(x) * l0(y);
    }

    static double gradleg_quad_l4_l0_bx(double x, double y)
    {
     return dl4(x) * dl0(y);
    }

    static double gradleg_quad_l4_l0_ay(double x, double y)
    {
     return dl4(x) * dl0(y);
    }

    static double gradleg_quad_l4_l0_by(double x, double y)
    {
     return l4(x) * d2l0(y);
    }

    static double gradleg_quad_l5_l0_a_0(double x, double y)
    {
     return  dl5(x) * l0(y);
    }

    static double gradleg_quad_l5_l0_a_1(double x, double y)
    {
     return -dl5(x) * l0(y);
    }

    static double gradleg_quad_l5_l0_b_0(double x, double y)
    {
     return  l5(x) * dl0(y);
    }

    static double gradleg_quad_l5_l0_b_1(double x, double y)
    {
     return -l5(x) * dl0(y);
    }

    static double gradleg_quad_l5_l0_ax_0(double x, double y)
    {
     return  d2l5(x) * l0(y);
    }

    static double gradleg_quad_l5_l0_ax_1(double x, double y)
    {
     return -d2l5(x) * l0(y);
    }

    static double gradleg_quad_l5_l0_bx_0(double x, double y)
    {
     return  dl5(x) * dl0(y);
    }

    static double gradleg_quad_l5_l0_bx_1(double x, double y)
    {
     return -dl5(x) * dl0(y);
    }

    static double gradleg_quad_l5_l0_ay_0(double x, double y)
    {
     return  dl5(x) * dl0(y);
    }

    static double gradleg_quad_l5_l0_ay_1(double x, double y)
    {
     return -dl5(x) * dl0(y);
    }

    static double gradleg_quad_l5_l0_by_0(double x, double y)
    {
     return  l5(x) * d2l0(y);
    }

    static double gradleg_quad_l5_l0_by_1(double x, double y)
    {
     return -l5(x) * d2l0(y);
    }

    static double gradleg_quad_l6_l0_a(double x, double y)
    {
     return dl6(x) * l0(y);
    }

    static double gradleg_quad_l6_l0_b(double x, double y)
    {
     return l6(x) * dl0(y);
    }

    static double gradleg_quad_l6_l0_ax(double x, double y)
    {
     return d2l6(x) * l0(y);
    }

    static double gradleg_quad_l6_l0_bx(double x, double y)
    {
     return dl6(x) * dl0(y);
    }

    static double gradleg_quad_l6_l0_ay(double x, double y)
    {
     return dl6(x) * dl0(y);
    }

    static double gradleg_quad_l6_l0_by(double x, double y)
    {
     return l6(x) * d2l0(y);
    }

    static double gradleg_quad_l7_l0_a_0(double x, double y)
    {
     return  dl7(x) * l0(y);
    }

    static double gradleg_quad_l7_l0_a_1(double x, double y)
    {
     return -dl7(x) * l0(y);
    }

    static double gradleg_quad_l7_l0_b_0(double x, double y)
    {
     return  l7(x) * dl0(y);
    }

    static double gradleg_quad_l7_l0_b_1(double x, double y)
    {
     return -l7(x) * dl0(y);
    }

    static double gradleg_quad_l7_l0_ax_0(double x, double y)
    {
     return  d2l7(x) * l0(y);
    }

    static double gradleg_quad_l7_l0_ax_1(double x, double y)
    {
     return -d2l7(x) * l0(y);
    }

    static double gradleg_quad_l7_l0_bx_0(double x, double y)
    {
     return  dl7(x) * dl0(y);
    }

    static double gradleg_quad_l7_l0_bx_1(double x, double y)
    {
     return -dl7(x) * dl0(y);
    }

    static double gradleg_quad_l7_l0_ay_0(double x, double y)
    {
     return  dl7(x) * dl0(y);
    }

    static double gradleg_quad_l7_l0_ay_1(double x, double y)
    {
     return -dl7(x) * dl0(y);
    }

    static double gradleg_quad_l7_l0_by_0(double x, double y)
    {
     return  l7(x) * d2l0(y);
    }

    static double gradleg_quad_l7_l0_by_1(double x, double y)
    {
     return -l7(x) * d2l0(y);
    }

    static double gradleg_quad_l8_l0_a(double x, double y)
    {
     return dl8(x) * l0(y);
    }

    static double gradleg_quad_l8_l0_b(double x, double y)
    {
     return l8(x) * dl0(y);
    }

    static double gradleg_quad_l8_l0_ax(double x, double y)
    {
     return d2l8(x) * l0(y);
    }

    static double gradleg_quad_l8_l0_bx(double x, double y)
    {
     return dl8(x) * dl0(y);
    }

    static double gradleg_quad_l8_l0_ay(double x, double y)
    {
     return dl8(x) * dl0(y);
    }

    static double gradleg_quad_l8_l0_by(double x, double y)
    {
     return l8(x) * d2l0(y);
    }

    static double gradleg_quad_l9_l0_a_0(double x, double y)
    {
     return  dl9(x) * l0(y);
    }

    static double gradleg_quad_l9_l0_a_1(double x, double y)
    {
     return -dl9(x) * l0(y);
    }

    static double gradleg_quad_l9_l0_b_0(double x, double y)
    {
     return  l9(x) * dl0(y);
    }

    static double gradleg_quad_l9_l0_b_1(double x, double y)
    {
     return -l9(x) * dl0(y);
    }

    static double gradleg_quad_l9_l0_ax_0(double x, double y)
    {
     return  d2l9(x) * l0(y);
    }

    static double gradleg_quad_l9_l0_ax_1(double x, double y)
    {
     return -d2l9(x) * l0(y);
    }

    static double gradleg_quad_l9_l0_bx_0(double x, double y)
    {
     return  dl9(x) * dl0(y);
    }

    static double gradleg_quad_l9_l0_bx_1(double x, double y)
    {
     return -dl9(x) * dl0(y);
    }

    static double gradleg_quad_l9_l0_ay_0(double x, double y)
    {
     return  dl9(x) * dl0(y);
    }

    static double gradleg_quad_l9_l0_ay_1(double x, double y)
    {
     return -dl9(x) * dl0(y);
    }

    static double gradleg_quad_l9_l0_by_0(double x, double y)
    {
     return  l9(x) * d2l0(y);
    }

    static double gradleg_quad_l9_l0_by_1(double x, double y)
    {
     return -l9(x) * d2l0(y);
    }

    static double gradleg_quad_l10_l0_a(double x, double y)
    {
     return dl10(x) * l0(y);
    }

    static double gradleg_quad_l10_l0_b(double x, double y)
    {
     return l10(x) * dl0(y);
    }

    static double gradleg_quad_l10_l0_ax(double x, double y)
    {
     return d2l10(x) * l0(y);
    }

    static double gradleg_quad_l10_l0_bx(double x, double y)
    {
     return dl10(x) * dl0(y);
    }

    static double gradleg_quad_l10_l0_ay(double x, double y)
    {
     return dl10(x) * dl0(y);
    }

    static double gradleg_quad_l10_l0_by(double x, double y)
    {
     return l10(x) * d2l0(y);
    }

    static double gradleg_quad_l11_l0_a_0(double x, double y)
    {
     return  dl11(x) * l0(y);
    }

    static double gradleg_quad_l11_l0_a_1(double x, double y)
    {
     return -dl11(x) * l0(y);
    }

    static double gradleg_quad_l11_l0_b_0(double x, double y)
    {
     return  l11(x) * dl0(y);
    }

    static double gradleg_quad_l11_l0_b_1(double x, double y)
    {
     return -l11(x) * dl0(y);
    }

    static double gradleg_quad_l11_l0_ax_0(double x, double y)
    {
     return  d2l11(x) * l0(y);
    }

    static double gradleg_quad_l11_l0_ax_1(double x, double y)
    {
     return -d2l11(x) * l0(y);
    }

    static double gradleg_quad_l11_l0_bx_0(double x, double y)
    {
     return  dl11(x) * dl0(y);
    }

    static double gradleg_quad_l11_l0_bx_1(double x, double y)
    {
     return -dl11(x) * dl0(y);
    }

    static double gradleg_quad_l11_l0_ay_0(double x, double y)
    {
     return  dl11(x) * dl0(y);
    }

    static double gradleg_quad_l11_l0_ay_1(double x, double y)
    {
     return -dl11(x) * dl0(y);
    }

    static double gradleg_quad_l11_l0_by_0(double x, double y)
    {
     return  l11(x) * d2l0(y);
    }

    static double gradleg_quad_l11_l0_by_1(double x, double y)
    {
     return -l11(x) * d2l0(y);
    }

    static double gradleg_quad_l2_l1_a(double x, double y)
    {
     return dl2(x) * l1(y);
    }

    static double gradleg_quad_l2_l1_b(double x, double y)
    {
     return l2(x) * dl1(y);
    }

    static double gradleg_quad_l2_l1_ax(double x, double y)
    {
     return d2l2(x) * l1(y);
    }

    static double gradleg_quad_l2_l1_bx(double x, double y)
    {
     return dl2(x) * dl1(y);
    }

    static double gradleg_quad_l2_l1_ay(double x, double y)
    {
     return dl2(x) * dl1(y);
    }

    static double gradleg_quad_l2_l1_by(double x, double y)
    {
     return l2(x) * d2l1(y);
    }

    static double gradleg_quad_l3_l1_a_0(double x, double y)
    {
     return -dl3(x) * l1(y);
    }

    static double gradleg_quad_l3_l1_a_1(double x, double y)
    {
     return  dl3(x) * l1(y);
    }

    static double gradleg_quad_l3_l1_b_0(double x, double y)
    {
     return -l3(x) * dl1(y);
    }

    static double gradleg_quad_l3_l1_b_1(double x, double y)
    {
     return  l3(x) * dl1(y);
    }

    static double gradleg_quad_l3_l1_ax_0(double x, double y)
    {
     return -d2l3(x) * l1(y);
    }

    static double gradleg_quad_l3_l1_ax_1(double x, double y)
    {
     return  d2l3(x) * l1(y);
    }

    static double gradleg_quad_l3_l1_bx_0(double x, double y)
    {
     return -dl3(x) * dl1(y);
    }

    static double gradleg_quad_l3_l1_bx_1(double x, double y)
    {
     return  dl3(x) * dl1(y);
    }

    static double gradleg_quad_l3_l1_ay_0(double x, double y)
    {
     return -dl3(x) * dl1(y);
    }

    static double gradleg_quad_l3_l1_ay_1(double x, double y)
    {
     return  dl3(x) * dl1(y);
    }

    static double gradleg_quad_l3_l1_by_0(double x, double y)
    {
     return -l3(x) * d2l1(y);
    }

    static double gradleg_quad_l3_l1_by_1(double x, double y)
    {
     return  l3(x) * d2l1(y);
    }

    static double gradleg_quad_l4_l1_a(double x, double y)
    {
     return dl4(x) * l1(y);
    }

    static double gradleg_quad_l4_l1_b(double x, double y)
    {
     return l4(x) * dl1(y);
    }

    static double gradleg_quad_l4_l1_ax(double x, double y)
    {
     return d2l4(x) * l1(y);
    }

    static double gradleg_quad_l4_l1_bx(double x, double y)
    {
     return dl4(x) * dl1(y);
    }

    static double gradleg_quad_l4_l1_ay(double x, double y)
    {
     return dl4(x) * dl1(y);
    }

    static double gradleg_quad_l4_l1_by(double x, double y)
    {
     return l4(x) * d2l1(y);
    }

    static double gradleg_quad_l5_l1_a_0(double x, double y)
    {
     return -dl5(x) * l1(y);
    }

    static double gradleg_quad_l5_l1_a_1(double x, double y)
    {
     return  dl5(x) * l1(y);
    }

    static double gradleg_quad_l5_l1_b_0(double x, double y)
    {
     return -l5(x) * dl1(y);
    }

    static double gradleg_quad_l5_l1_b_1(double x, double y)
    {
     return  l5(x) * dl1(y);
    }

    static double gradleg_quad_l5_l1_ax_0(double x, double y)
    {
     return -d2l5(x) * l1(y);
    }

    static double gradleg_quad_l5_l1_ax_1(double x, double y)
    {
     return  d2l5(x) * l1(y);
    }

    static double gradleg_quad_l5_l1_bx_0(double x, double y)
    {
     return -dl5(x) * dl1(y);
    }

    static double gradleg_quad_l5_l1_bx_1(double x, double y)
    {
     return  dl5(x) * dl1(y);
    }

    static double gradleg_quad_l5_l1_ay_0(double x, double y)
    {
     return -dl5(x) * dl1(y);
    }

    static double gradleg_quad_l5_l1_ay_1(double x, double y)
    {
     return  dl5(x) * dl1(y);
    }

    static double gradleg_quad_l5_l1_by_0(double x, double y)
    {
     return -l5(x) * d2l1(y);
    }

    static double gradleg_quad_l5_l1_by_1(double x, double y)
    {
     return  l5(x) * d2l1(y);
    }

    static double gradleg_quad_l6_l1_a(double x, double y)
    {
     return dl6(x) * l1(y);
    }

    static double gradleg_quad_l6_l1_b(double x, double y)
    {
     return l6(x) * dl1(y);
    }

    static double gradleg_quad_l6_l1_ax(double x, double y)
    {
     return d2l6(x) * l1(y);
    }

    static double gradleg_quad_l6_l1_bx(double x, double y)
    {
     return dl6(x) * dl1(y);
    }

    static double gradleg_quad_l6_l1_ay(double x, double y)
    {
     return dl6(x) * dl1(y);
    }

    static double gradleg_quad_l6_l1_by(double x, double y)
    {
     return l6(x) * d2l1(y);
    }

    static double gradleg_quad_l7_l1_a_0(double x, double y)
    {
     return -dl7(x) * l1(y);
    }

    static double gradleg_quad_l7_l1_a_1(double x, double y)
    {
     return  dl7(x) * l1(y);
    }

    static double gradleg_quad_l7_l1_b_0(double x, double y)
    {
     return -l7(x) * dl1(y);
    }

    static double gradleg_quad_l7_l1_b_1(double x, double y)
    {
     return  l7(x) * dl1(y);
    }

    static double gradleg_quad_l7_l1_ax_0(double x, double y)
    {
     return -d2l7(x) * l1(y);
    }

    static double gradleg_quad_l7_l1_ax_1(double x, double y)
    {
     return  d2l7(x) * l1(y);
    }

    static double gradleg_quad_l7_l1_bx_0(double x, double y)
    {
     return -dl7(x) * dl1(y);
    }

    static double gradleg_quad_l7_l1_bx_1(double x, double y)
    {
     return  dl7(x) * dl1(y);
    }

    static double gradleg_quad_l7_l1_ay_0(double x, double y)
    {
     return -dl7(x) * dl1(y);
    }

    static double gradleg_quad_l7_l1_ay_1(double x, double y)
    {
     return  dl7(x) * dl1(y);
    }

    static double gradleg_quad_l7_l1_by_0(double x, double y)
    {
     return -l7(x) * d2l1(y);
    }

    static double gradleg_quad_l7_l1_by_1(double x, double y)
    {
     return  l7(x) * d2l1(y);
    }

    static double gradleg_quad_l8_l1_a(double x, double y)
    {
     return dl8(x) * l1(y);
    }

    static double gradleg_quad_l8_l1_b(double x, double y)
    {
     return l8(x) * dl1(y);
    }

    static double gradleg_quad_l8_l1_ax(double x, double y)
    {
     return d2l8(x) * l1(y);
    }

    static double gradleg_quad_l8_l1_bx(double x, double y)
    {
     return dl8(x) * dl1(y);
    }

    static double gradleg_quad_l8_l1_ay(double x, double y)
    {
     return dl8(x) * dl1(y);
    }

    static double gradleg_quad_l8_l1_by(double x, double y)
    {
     return l8(x) * d2l1(y);
    }

    static double gradleg_quad_l9_l1_a_0(double x, double y)
    {
     return -dl9(x) * l1(y);
    }

    static double gradleg_quad_l9_l1_a_1(double x, double y)
    {
     return  dl9(x) * l1(y);
    }

    static double gradleg_quad_l9_l1_b_0(double x, double y)
    {
     return -l9(x) * dl1(y);
    }

    static double gradleg_quad_l9_l1_b_1(double x, double y)
    {
     return  l9(x) * dl1(y);
    }

    static double gradleg_quad_l9_l1_ax_0(double x, double y)
    {
     return -d2l9(x) * l1(y);
    }

    static double gradleg_quad_l9_l1_ax_1(double x, double y)
    {
     return  d2l9(x) * l1(y);
    }

    static double gradleg_quad_l9_l1_bx_0(double x, double y)
    {
     return -dl9(x) * dl1(y);
    }

    static double gradleg_quad_l9_l1_bx_1(double x, double y)
    {
     return  dl9(x) * dl1(y);
    }

    static double gradleg_quad_l9_l1_ay_0(double x, double y)
    {
     return -dl9(x) * dl1(y);
    }

    static double gradleg_quad_l9_l1_ay_1(double x, double y)
    {
     return  dl9(x) * dl1(y);
    }

    static double gradleg_quad_l9_l1_by_0(double x, double y)
    {
     return -l9(x) * d2l1(y);
    }

    static double gradleg_quad_l9_l1_by_1(double x, double y)
    {
     return  l9(x) * d2l1(y);
    }

    static double gradleg_quad_l10_l1_a(double x, double y)
    {
     return dl10(x) * l1(y);
    }

    static double gradleg_quad_l10_l1_b(double x, double y)
    {
     return l10(x) * dl1(y);
    }

    static double gradleg_quad_l10_l1_ax(double x, double y)
    {
     return d2l10(x) * l1(y);
    }

    static double gradleg_quad_l10_l1_bx(double x, double y)
    {
     return dl10(x) * dl1(y);
    }

    static double gradleg_quad_l10_l1_ay(double x, double y)
    {
     return dl10(x) * dl1(y);
    }

    static double gradleg_quad_l10_l1_by(double x, double y)
    {
     return l10(x) * d2l1(y);
    }

    static double gradleg_quad_l11_l1_a_0(double x, double y)
    {
     return -dl11(x) * l1(y);
    }

    static double gradleg_quad_l11_l1_a_1(double x, double y)
    {
     return  dl11(x) * l1(y);
    }

    static double gradleg_quad_l11_l1_b_0(double x, double y)
    {
     return -l11(x) * dl1(y);
    }

    static double gradleg_quad_l11_l1_b_1(double x, double y)
    {
     return  l11(x) * dl1(y);
    }

    static double gradleg_quad_l11_l1_ax_0(double x, double y)
    {
     return -d2l11(x) * l1(y);
    }

    static double gradleg_quad_l11_l1_ax_1(double x, double y)
    {
     return  d2l11(x) * l1(y);
    }

    static double gradleg_quad_l11_l1_bx_0(double x, double y)
    {
     return -dl11(x) * dl1(y);
    }

    static double gradleg_quad_l11_l1_bx_1(double x, double y)
    {
     return  dl11(x) * dl1(y);
    }

    static double gradleg_quad_l11_l1_ay_0(double x, double y)
    {
     return -dl11(x) * dl1(y);
    }

    static double gradleg_quad_l11_l1_ay_1(double x, double y)
    {
     return  dl11(x) * dl1(y);
    }

    static double gradleg_quad_l11_l1_by_0(double x, double y)
    {
     return -l11(x) * d2l1(y);
    }

    static double gradleg_quad_l11_l1_by_1(double x, double y)
    {
     return  l11(x) * d2l1(y);
    }

    /* BUBBLE */

    /* BUBBLE ( 1 , 0 ) */

    static double gradleg_quad_p0p2_b1_a(double x, double y)
    {
      return Legendre0(x) * l2(y);
    }

    static double gradleg_quad_p0p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p2_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l2(y);
    }

    static double gradleg_quad_p0p2_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl2(y);
    }

    static double gradleg_quad_p0p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p3_b1_a(double x, double y)
    {
      return Legendre0(x) * l3(y);
    }

    static double gradleg_quad_p0p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p3_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l3(y);
    }

    static double gradleg_quad_p0p3_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl3(y);
    }

    static double gradleg_quad_p0p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p4_b1_a(double x, double y)
    {
      return Legendre0(x) * l4(y);
    }

    static double gradleg_quad_p0p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p4_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l4(y);
    }

    static double gradleg_quad_p0p4_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl4(y);
    }

    static double gradleg_quad_p0p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p5_b1_a(double x, double y)
    {
      return Legendre0(x) * l5(y);
    }

    static double gradleg_quad_p0p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p5_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l5(y);
    }

    static double gradleg_quad_p0p5_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl5(y);
    }

    static double gradleg_quad_p0p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p6_b1_a(double x, double y)
    {
      return Legendre0(x) * l6(y);
    }

    static double gradleg_quad_p0p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p6_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l6(y);
    }

    static double gradleg_quad_p0p6_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl6(y);
    }

    static double gradleg_quad_p0p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p7_b1_a(double x, double y)
    {
      return Legendre0(x) * l7(y);
    }

    static double gradleg_quad_p0p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p7_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l7(y);
    }

    static double gradleg_quad_p0p7_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl7(y);
    }

    static double gradleg_quad_p0p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p8_b1_a(double x, double y)
    {
      return Legendre0(x) * l8(y);
    }

    static double gradleg_quad_p0p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p8_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l8(y);
    }

    static double gradleg_quad_p0p8_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl8(y);
    }

    static double gradleg_quad_p0p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p9_b1_a(double x, double y)
    {
      return Legendre0(x) * l9(y);
    }

    static double gradleg_quad_p0p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p9_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l9(y);
    }

    static double gradleg_quad_p0p9_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl9(y);
    }

    static double gradleg_quad_p0p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p10_b1_a(double x, double y)
    {
      return Legendre0(x) * l10(y);
    }

    static double gradleg_quad_p0p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p10_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l10(y);
    }

    static double gradleg_quad_p0p10_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl10(y);
    }

    static double gradleg_quad_p0p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p11_b1_a(double x, double y)
    {
      return Legendre0(x) * l11(y);
    }

    static double gradleg_quad_p0p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p11_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l11(y);
    }

    static double gradleg_quad_p0p11_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl11(y);
    }

    static double gradleg_quad_p0p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p0p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p2_b1_a(double x, double y)
    {
      return Legendre1(x) * l2(y);
    }

    static double gradleg_quad_p1p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p2_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l2(y);
    }

    static double gradleg_quad_p1p2_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl2(y);
    }

    static double gradleg_quad_p1p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p3_b1_a(double x, double y)
    {
      return Legendre1(x) * l3(y);
    }

    static double gradleg_quad_p1p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p3_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l3(y);
    }

    static double gradleg_quad_p1p3_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl3(y);
    }

    static double gradleg_quad_p1p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p4_b1_a(double x, double y)
    {
      return Legendre1(x) * l4(y);
    }

    static double gradleg_quad_p1p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p4_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l4(y);
    }

    static double gradleg_quad_p1p4_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl4(y);
    }

    static double gradleg_quad_p1p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p5_b1_a(double x, double y)
    {
      return Legendre1(x) * l5(y);
    }

    static double gradleg_quad_p1p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p5_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l5(y);
    }

    static double gradleg_quad_p1p5_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl5(y);
    }

    static double gradleg_quad_p1p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p6_b1_a(double x, double y)
    {
      return Legendre1(x) * l6(y);
    }

    static double gradleg_quad_p1p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p6_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l6(y);
    }

    static double gradleg_quad_p1p6_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl6(y);
    }

    static double gradleg_quad_p1p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p7_b1_a(double x, double y)
    {
      return Legendre1(x) * l7(y);
    }

    static double gradleg_quad_p1p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p7_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l7(y);
    }

    static double gradleg_quad_p1p7_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl7(y);
    }

    static double gradleg_quad_p1p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p8_b1_a(double x, double y)
    {
      return Legendre1(x) * l8(y);
    }

    static double gradleg_quad_p1p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p8_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l8(y);
    }

    static double gradleg_quad_p1p8_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl8(y);
    }

    static double gradleg_quad_p1p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p9_b1_a(double x, double y)
    {
      return Legendre1(x) * l9(y);
    }

    static double gradleg_quad_p1p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p9_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l9(y);
    }

    static double gradleg_quad_p1p9_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl9(y);
    }

    static double gradleg_quad_p1p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p10_b1_a(double x, double y)
    {
      return Legendre1(x) * l10(y);
    }

    static double gradleg_quad_p1p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p10_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l10(y);
    }

    static double gradleg_quad_p1p10_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl10(y);
    }

    static double gradleg_quad_p1p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p11_b1_a(double x, double y)
    {
      return Legendre1(x) * l11(y);
    }

    static double gradleg_quad_p1p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p11_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l11(y);
    }

    static double gradleg_quad_p1p11_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl11(y);
    }

    static double gradleg_quad_p1p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p1p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p2_b1_a(double x, double y)
    {
      return Legendre2(x) * l2(y);
    }

    static double gradleg_quad_p2p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p2_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l2(y);
    }

    static double gradleg_quad_p2p2_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl2(y);
    }

    static double gradleg_quad_p2p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p3_b1_a(double x, double y)
    {
      return Legendre2(x) * l3(y);
    }

    static double gradleg_quad_p2p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p3_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l3(y);
    }

    static double gradleg_quad_p2p3_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl3(y);
    }

    static double gradleg_quad_p2p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p4_b1_a(double x, double y)
    {
      return Legendre2(x) * l4(y);
    }

    static double gradleg_quad_p2p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p4_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l4(y);
    }

    static double gradleg_quad_p2p4_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl4(y);
    }

    static double gradleg_quad_p2p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p5_b1_a(double x, double y)
    {
      return Legendre2(x) * l5(y);
    }

    static double gradleg_quad_p2p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p5_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l5(y);
    }

    static double gradleg_quad_p2p5_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl5(y);
    }

    static double gradleg_quad_p2p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p6_b1_a(double x, double y)
    {
      return Legendre2(x) * l6(y);
    }

    static double gradleg_quad_p2p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p6_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l6(y);
    }

    static double gradleg_quad_p2p6_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl6(y);
    }

    static double gradleg_quad_p2p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p7_b1_a(double x, double y)
    {
      return Legendre2(x) * l7(y);
    }

    static double gradleg_quad_p2p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p7_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l7(y);
    }

    static double gradleg_quad_p2p7_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl7(y);
    }

    static double gradleg_quad_p2p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p8_b1_a(double x, double y)
    {
      return Legendre2(x) * l8(y);
    }

    static double gradleg_quad_p2p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p8_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l8(y);
    }

    static double gradleg_quad_p2p8_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl8(y);
    }

    static double gradleg_quad_p2p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p9_b1_a(double x, double y)
    {
      return Legendre2(x) * l9(y);
    }

    static double gradleg_quad_p2p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p9_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l9(y);
    }

    static double gradleg_quad_p2p9_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl9(y);
    }

    static double gradleg_quad_p2p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p10_b1_a(double x, double y)
    {
      return Legendre2(x) * l10(y);
    }

    static double gradleg_quad_p2p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p10_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l10(y);
    }

    static double gradleg_quad_p2p10_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl10(y);
    }

    static double gradleg_quad_p2p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p11_b1_a(double x, double y)
    {
      return Legendre2(x) * l11(y);
    }

    static double gradleg_quad_p2p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p11_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l11(y);
    }

    static double gradleg_quad_p2p11_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl11(y);
    }

    static double gradleg_quad_p2p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p2_b1_a(double x, double y)
    {
      return Legendre3(x) * l2(y);
    }

    static double gradleg_quad_p3p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p2_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l2(y);
    }

    static double gradleg_quad_p3p2_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl2(y);
    }

    static double gradleg_quad_p3p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p3_b1_a(double x, double y)
    {
      return Legendre3(x) * l3(y);
    }

    static double gradleg_quad_p3p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p3_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l3(y);
    }

    static double gradleg_quad_p3p3_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl3(y);
    }

    static double gradleg_quad_p3p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p4_b1_a(double x, double y)
    {
      return Legendre3(x) * l4(y);
    }

    static double gradleg_quad_p3p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p4_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l4(y);
    }

    static double gradleg_quad_p3p4_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl4(y);
    }

    static double gradleg_quad_p3p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p5_b1_a(double x, double y)
    {
      return Legendre3(x) * l5(y);
    }

    static double gradleg_quad_p3p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p5_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l5(y);
    }

    static double gradleg_quad_p3p5_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl5(y);
    }

    static double gradleg_quad_p3p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p6_b1_a(double x, double y)
    {
      return Legendre3(x) * l6(y);
    }

    static double gradleg_quad_p3p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p6_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l6(y);
    }

    static double gradleg_quad_p3p6_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl6(y);
    }

    static double gradleg_quad_p3p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p7_b1_a(double x, double y)
    {
      return Legendre3(x) * l7(y);
    }

    static double gradleg_quad_p3p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p7_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l7(y);
    }

    static double gradleg_quad_p3p7_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl7(y);
    }

    static double gradleg_quad_p3p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p8_b1_a(double x, double y)
    {
      return Legendre3(x) * l8(y);
    }

    static double gradleg_quad_p3p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p8_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l8(y);
    }

    static double gradleg_quad_p3p8_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl8(y);
    }

    static double gradleg_quad_p3p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p9_b1_a(double x, double y)
    {
      return Legendre3(x) * l9(y);
    }

    static double gradleg_quad_p3p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p9_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l9(y);
    }

    static double gradleg_quad_p3p9_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl9(y);
    }

    static double gradleg_quad_p3p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p10_b1_a(double x, double y)
    {
      return Legendre3(x) * l10(y);
    }

    static double gradleg_quad_p3p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p10_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l10(y);
    }

    static double gradleg_quad_p3p10_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl10(y);
    }

    static double gradleg_quad_p3p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p11_b1_a(double x, double y)
    {
      return Legendre3(x) * l11(y);
    }

    static double gradleg_quad_p3p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p11_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l11(y);
    }

    static double gradleg_quad_p3p11_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl11(y);
    }

    static double gradleg_quad_p3p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p2_b1_a(double x, double y)
    {
      return Legendre4(x) * l2(y);
    }

    static double gradleg_quad_p4p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p2_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l2(y);
    }

    static double gradleg_quad_p4p2_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl2(y);
    }

    static double gradleg_quad_p4p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p3_b1_a(double x, double y)
    {
      return Legendre4(x) * l3(y);
    }

    static double gradleg_quad_p4p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p3_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l3(y);
    }

    static double gradleg_quad_p4p3_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl3(y);
    }

    static double gradleg_quad_p4p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p4_b1_a(double x, double y)
    {
      return Legendre4(x) * l4(y);
    }

    static double gradleg_quad_p4p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p4_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l4(y);
    }

    static double gradleg_quad_p4p4_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl4(y);
    }

    static double gradleg_quad_p4p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p5_b1_a(double x, double y)
    {
      return Legendre4(x) * l5(y);
    }

    static double gradleg_quad_p4p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p5_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l5(y);
    }

    static double gradleg_quad_p4p5_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl5(y);
    }

    static double gradleg_quad_p4p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p6_b1_a(double x, double y)
    {
      return Legendre4(x) * l6(y);
    }

    static double gradleg_quad_p4p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p6_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l6(y);
    }

    static double gradleg_quad_p4p6_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl6(y);
    }

    static double gradleg_quad_p4p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p7_b1_a(double x, double y)
    {
      return Legendre4(x) * l7(y);
    }

    static double gradleg_quad_p4p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p7_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l7(y);
    }

    static double gradleg_quad_p4p7_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl7(y);
    }

    static double gradleg_quad_p4p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p8_b1_a(double x, double y)
    {
      return Legendre4(x) * l8(y);
    }

    static double gradleg_quad_p4p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p8_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l8(y);
    }

    static double gradleg_quad_p4p8_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl8(y);
    }

    static double gradleg_quad_p4p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p9_b1_a(double x, double y)
    {
      return Legendre4(x) * l9(y);
    }

    static double gradleg_quad_p4p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p9_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l9(y);
    }

    static double gradleg_quad_p4p9_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl9(y);
    }

    static double gradleg_quad_p4p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p10_b1_a(double x, double y)
    {
      return Legendre4(x) * l10(y);
    }

    static double gradleg_quad_p4p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p10_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l10(y);
    }

    static double gradleg_quad_p4p10_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl10(y);
    }

    static double gradleg_quad_p4p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p11_b1_a(double x, double y)
    {
      return Legendre4(x) * l11(y);
    }

    static double gradleg_quad_p4p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p11_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l11(y);
    }

    static double gradleg_quad_p4p11_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl11(y);
    }

    static double gradleg_quad_p4p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p2_b1_a(double x, double y)
    {
      return Legendre5(x) * l2(y);
    }

    static double gradleg_quad_p5p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p2_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l2(y);
    }

    static double gradleg_quad_p5p2_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl2(y);
    }

    static double gradleg_quad_p5p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p3_b1_a(double x, double y)
    {
      return Legendre5(x) * l3(y);
    }

    static double gradleg_quad_p5p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p3_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l3(y);
    }

    static double gradleg_quad_p5p3_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl3(y);
    }

    static double gradleg_quad_p5p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p4_b1_a(double x, double y)
    {
      return Legendre5(x) * l4(y);
    }

    static double gradleg_quad_p5p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p4_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l4(y);
    }

    static double gradleg_quad_p5p4_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl4(y);
    }

    static double gradleg_quad_p5p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p5_b1_a(double x, double y)
    {
      return Legendre5(x) * l5(y);
    }

    static double gradleg_quad_p5p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p5_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l5(y);
    }

    static double gradleg_quad_p5p5_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl5(y);
    }

    static double gradleg_quad_p5p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p6_b1_a(double x, double y)
    {
      return Legendre5(x) * l6(y);
    }

    static double gradleg_quad_p5p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p6_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l6(y);
    }

    static double gradleg_quad_p5p6_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl6(y);
    }

    static double gradleg_quad_p5p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p7_b1_a(double x, double y)
    {
      return Legendre5(x) * l7(y);
    }

    static double gradleg_quad_p5p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p7_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l7(y);
    }

    static double gradleg_quad_p5p7_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl7(y);
    }

    static double gradleg_quad_p5p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p8_b1_a(double x, double y)
    {
      return Legendre5(x) * l8(y);
    }

    static double gradleg_quad_p5p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p8_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l8(y);
    }

    static double gradleg_quad_p5p8_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl8(y);
    }

    static double gradleg_quad_p5p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p9_b1_a(double x, double y)
    {
      return Legendre5(x) * l9(y);
    }

    static double gradleg_quad_p5p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p9_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l9(y);
    }

    static double gradleg_quad_p5p9_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl9(y);
    }

    static double gradleg_quad_p5p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p10_b1_a(double x, double y)
    {
      return Legendre5(x) * l10(y);
    }

    static double gradleg_quad_p5p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p10_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l10(y);
    }

    static double gradleg_quad_p5p10_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl10(y);
    }

    static double gradleg_quad_p5p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p11_b1_a(double x, double y)
    {
      return Legendre5(x) * l11(y);
    }

    static double gradleg_quad_p5p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p11_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l11(y);
    }

    static double gradleg_quad_p5p11_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl11(y);
    }

    static double gradleg_quad_p5p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p2_b1_a(double x, double y)
    {
      return Legendre6(x) * l2(y);
    }

    static double gradleg_quad_p6p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p2_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l2(y);
    }

    static double gradleg_quad_p6p2_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl2(y);
    }

    static double gradleg_quad_p6p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p3_b1_a(double x, double y)
    {
      return Legendre6(x) * l3(y);
    }

    static double gradleg_quad_p6p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p3_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l3(y);
    }

    static double gradleg_quad_p6p3_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl3(y);
    }

    static double gradleg_quad_p6p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p4_b1_a(double x, double y)
    {
      return Legendre6(x) * l4(y);
    }

    static double gradleg_quad_p6p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p4_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l4(y);
    }

    static double gradleg_quad_p6p4_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl4(y);
    }

    static double gradleg_quad_p6p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p5_b1_a(double x, double y)
    {
      return Legendre6(x) * l5(y);
    }

    static double gradleg_quad_p6p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p5_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l5(y);
    }

    static double gradleg_quad_p6p5_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl5(y);
    }

    static double gradleg_quad_p6p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p6_b1_a(double x, double y)
    {
      return Legendre6(x) * l6(y);
    }

    static double gradleg_quad_p6p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p6_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l6(y);
    }

    static double gradleg_quad_p6p6_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl6(y);
    }

    static double gradleg_quad_p6p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p7_b1_a(double x, double y)
    {
      return Legendre6(x) * l7(y);
    }

    static double gradleg_quad_p6p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p7_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l7(y);
    }

    static double gradleg_quad_p6p7_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl7(y);
    }

    static double gradleg_quad_p6p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p8_b1_a(double x, double y)
    {
      return Legendre6(x) * l8(y);
    }

    static double gradleg_quad_p6p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p8_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l8(y);
    }

    static double gradleg_quad_p6p8_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl8(y);
    }

    static double gradleg_quad_p6p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p9_b1_a(double x, double y)
    {
      return Legendre6(x) * l9(y);
    }

    static double gradleg_quad_p6p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p9_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l9(y);
    }

    static double gradleg_quad_p6p9_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl9(y);
    }

    static double gradleg_quad_p6p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p10_b1_a(double x, double y)
    {
      return Legendre6(x) * l10(y);
    }

    static double gradleg_quad_p6p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p10_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l10(y);
    }

    static double gradleg_quad_p6p10_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl10(y);
    }

    static double gradleg_quad_p6p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p11_b1_a(double x, double y)
    {
      return Legendre6(x) * l11(y);
    }

    static double gradleg_quad_p6p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p11_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l11(y);
    }

    static double gradleg_quad_p6p11_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl11(y);
    }

    static double gradleg_quad_p6p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p2_b1_a(double x, double y)
    {
      return Legendre7(x) * l2(y);
    }

    static double gradleg_quad_p7p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p2_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l2(y);
    }

    static double gradleg_quad_p7p2_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl2(y);
    }

    static double gradleg_quad_p7p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p3_b1_a(double x, double y)
    {
      return Legendre7(x) * l3(y);
    }

    static double gradleg_quad_p7p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p3_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l3(y);
    }

    static double gradleg_quad_p7p3_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl3(y);
    }

    static double gradleg_quad_p7p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p4_b1_a(double x, double y)
    {
      return Legendre7(x) * l4(y);
    }

    static double gradleg_quad_p7p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p4_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l4(y);
    }

    static double gradleg_quad_p7p4_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl4(y);
    }

    static double gradleg_quad_p7p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p5_b1_a(double x, double y)
    {
      return Legendre7(x) * l5(y);
    }

    static double gradleg_quad_p7p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p5_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l5(y);
    }

    static double gradleg_quad_p7p5_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl5(y);
    }

    static double gradleg_quad_p7p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p6_b1_a(double x, double y)
    {
      return Legendre7(x) * l6(y);
    }

    static double gradleg_quad_p7p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p6_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l6(y);
    }

    static double gradleg_quad_p7p6_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl6(y);
    }

    static double gradleg_quad_p7p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p7_b1_a(double x, double y)
    {
      return Legendre7(x) * l7(y);
    }

    static double gradleg_quad_p7p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p7_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l7(y);
    }

    static double gradleg_quad_p7p7_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl7(y);
    }

    static double gradleg_quad_p7p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p8_b1_a(double x, double y)
    {
      return Legendre7(x) * l8(y);
    }

    static double gradleg_quad_p7p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p8_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l8(y);
    }

    static double gradleg_quad_p7p8_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl8(y);
    }

    static double gradleg_quad_p7p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p9_b1_a(double x, double y)
    {
      return Legendre7(x) * l9(y);
    }

    static double gradleg_quad_p7p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p9_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l9(y);
    }

    static double gradleg_quad_p7p9_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl9(y);
    }

    static double gradleg_quad_p7p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p10_b1_a(double x, double y)
    {
      return Legendre7(x) * l10(y);
    }

    static double gradleg_quad_p7p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p10_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l10(y);
    }

    static double gradleg_quad_p7p10_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl10(y);
    }

    static double gradleg_quad_p7p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p11_b1_a(double x, double y)
    {
      return Legendre7(x) * l11(y);
    }

    static double gradleg_quad_p7p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p11_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l11(y);
    }

    static double gradleg_quad_p7p11_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl11(y);
    }

    static double gradleg_quad_p7p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p2_b1_a(double x, double y)
    {
      return Legendre8(x) * l2(y);
    }

    static double gradleg_quad_p8p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p2_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l2(y);
    }

    static double gradleg_quad_p8p2_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl2(y);
    }

    static double gradleg_quad_p8p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p3_b1_a(double x, double y)
    {
      return Legendre8(x) * l3(y);
    }

    static double gradleg_quad_p8p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p3_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l3(y);
    }

    static double gradleg_quad_p8p3_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl3(y);
    }

    static double gradleg_quad_p8p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p4_b1_a(double x, double y)
    {
      return Legendre8(x) * l4(y);
    }

    static double gradleg_quad_p8p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p4_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l4(y);
    }

    static double gradleg_quad_p8p4_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl4(y);
    }

    static double gradleg_quad_p8p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p5_b1_a(double x, double y)
    {
      return Legendre8(x) * l5(y);
    }

    static double gradleg_quad_p8p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p5_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l5(y);
    }

    static double gradleg_quad_p8p5_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl5(y);
    }

    static double gradleg_quad_p8p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p6_b1_a(double x, double y)
    {
      return Legendre8(x) * l6(y);
    }

    static double gradleg_quad_p8p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p6_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l6(y);
    }

    static double gradleg_quad_p8p6_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl6(y);
    }

    static double gradleg_quad_p8p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p7_b1_a(double x, double y)
    {
      return Legendre8(x) * l7(y);
    }

    static double gradleg_quad_p8p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p7_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l7(y);
    }

    static double gradleg_quad_p8p7_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl7(y);
    }

    static double gradleg_quad_p8p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p8_b1_a(double x, double y)
    {
      return Legendre8(x) * l8(y);
    }

    static double gradleg_quad_p8p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p8_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l8(y);
    }

    static double gradleg_quad_p8p8_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl8(y);
    }

    static double gradleg_quad_p8p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p9_b1_a(double x, double y)
    {
      return Legendre8(x) * l9(y);
    }

    static double gradleg_quad_p8p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p9_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l9(y);
    }

    static double gradleg_quad_p8p9_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl9(y);
    }

    static double gradleg_quad_p8p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p10_b1_a(double x, double y)
    {
      return Legendre8(x) * l10(y);
    }

    static double gradleg_quad_p8p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p10_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l10(y);
    }

    static double gradleg_quad_p8p10_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl10(y);
    }

    static double gradleg_quad_p8p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p11_b1_a(double x, double y)
    {
      return Legendre8(x) * l11(y);
    }

    static double gradleg_quad_p8p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p11_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l11(y);
    }

    static double gradleg_quad_p8p11_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl11(y);
    }

    static double gradleg_quad_p8p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p2_b1_a(double x, double y)
    {
      return Legendre9(x) * l2(y);
    }

    static double gradleg_quad_p9p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p2_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l2(y);
    }

    static double gradleg_quad_p9p2_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl2(y);
    }

    static double gradleg_quad_p9p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p3_b1_a(double x, double y)
    {
      return Legendre9(x) * l3(y);
    }

    static double gradleg_quad_p9p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p3_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l3(y);
    }

    static double gradleg_quad_p9p3_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl3(y);
    }

    static double gradleg_quad_p9p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p4_b1_a(double x, double y)
    {
      return Legendre9(x) * l4(y);
    }

    static double gradleg_quad_p9p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p4_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l4(y);
    }

    static double gradleg_quad_p9p4_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl4(y);
    }

    static double gradleg_quad_p9p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p5_b1_a(double x, double y)
    {
      return Legendre9(x) * l5(y);
    }

    static double gradleg_quad_p9p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p5_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l5(y);
    }

    static double gradleg_quad_p9p5_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl5(y);
    }

    static double gradleg_quad_p9p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p6_b1_a(double x, double y)
    {
      return Legendre9(x) * l6(y);
    }

    static double gradleg_quad_p9p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p6_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l6(y);
    }

    static double gradleg_quad_p9p6_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl6(y);
    }

    static double gradleg_quad_p9p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p7_b1_a(double x, double y)
    {
      return Legendre9(x) * l7(y);
    }

    static double gradleg_quad_p9p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p7_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l7(y);
    }

    static double gradleg_quad_p9p7_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl7(y);
    }

    static double gradleg_quad_p9p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p8_b1_a(double x, double y)
    {
      return Legendre9(x) * l8(y);
    }

    static double gradleg_quad_p9p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p8_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l8(y);
    }

    static double gradleg_quad_p9p8_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl8(y);
    }

    static double gradleg_quad_p9p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p9_b1_a(double x, double y)
    {
      return Legendre9(x) * l9(y);
    }

    static double gradleg_quad_p9p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p9_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l9(y);
    }

    static double gradleg_quad_p9p9_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl9(y);
    }

    static double gradleg_quad_p9p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p10_b1_a(double x, double y)
    {
      return Legendre9(x) * l10(y);
    }

    static double gradleg_quad_p9p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p10_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l10(y);
    }

    static double gradleg_quad_p9p10_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl10(y);
    }

    static double gradleg_quad_p9p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p11_b1_a(double x, double y)
    {
      return Legendre9(x) * l11(y);
    }

    static double gradleg_quad_p9p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p11_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l11(y);
    }

    static double gradleg_quad_p9p11_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl11(y);
    }

    static double gradleg_quad_p9p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p2_b1_a(double x, double y)
    {
      return Legendre10(x) * l2(y);
    }

    static double gradleg_quad_p10p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p2_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l2(y);
    }

    static double gradleg_quad_p10p2_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl2(y);
    }

    static double gradleg_quad_p10p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p3_b1_a(double x, double y)
    {
      return Legendre10(x) * l3(y);
    }

    static double gradleg_quad_p10p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p3_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l3(y);
    }

    static double gradleg_quad_p10p3_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl3(y);
    }

    static double gradleg_quad_p10p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p4_b1_a(double x, double y)
    {
      return Legendre10(x) * l4(y);
    }

    static double gradleg_quad_p10p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p4_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l4(y);
    }

    static double gradleg_quad_p10p4_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl4(y);
    }

    static double gradleg_quad_p10p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p5_b1_a(double x, double y)
    {
      return Legendre10(x) * l5(y);
    }

    static double gradleg_quad_p10p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p5_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l5(y);
    }

    static double gradleg_quad_p10p5_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl5(y);
    }

    static double gradleg_quad_p10p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p6_b1_a(double x, double y)
    {
      return Legendre10(x) * l6(y);
    }

    static double gradleg_quad_p10p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p6_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l6(y);
    }

    static double gradleg_quad_p10p6_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl6(y);
    }

    static double gradleg_quad_p10p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p7_b1_a(double x, double y)
    {
      return Legendre10(x) * l7(y);
    }

    static double gradleg_quad_p10p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p7_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l7(y);
    }

    static double gradleg_quad_p10p7_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl7(y);
    }

    static double gradleg_quad_p10p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p8_b1_a(double x, double y)
    {
      return Legendre10(x) * l8(y);
    }

    static double gradleg_quad_p10p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p8_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l8(y);
    }

    static double gradleg_quad_p10p8_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl8(y);
    }

    static double gradleg_quad_p10p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p9_b1_a(double x, double y)
    {
      return Legendre10(x) * l9(y);
    }

    static double gradleg_quad_p10p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p9_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l9(y);
    }

    static double gradleg_quad_p10p9_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl9(y);
    }

    static double gradleg_quad_p10p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p10_b1_a(double x, double y)
    {
      return Legendre10(x) * l10(y);
    }

    static double gradleg_quad_p10p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p10_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l10(y);
    }

    static double gradleg_quad_p10p10_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl10(y);
    }

    static double gradleg_quad_p10p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p11_b1_a(double x, double y)
    {
      return Legendre10(x) * l11(y);
    }

    static double gradleg_quad_p10p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p11_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l11(y);
    }

    static double gradleg_quad_p10p11_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl11(y);
    }

    static double gradleg_quad_p10p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    /* BUBBLE ( 0 , 1 ) */

    static double gradleg_quad_p2p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p0_b2_b(double x, double y)
    {
      return l2(x) * Legendre0(y);
    }

    static double gradleg_quad_p2p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p0_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre0(y);
    }

    static double gradleg_quad_p2p0_b2_by(double x, double y)
    {
      return l2(x) * Legendre0x(y);
    }

    static double gradleg_quad_p2p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p1_b2_b(double x, double y)
    {
      return l2(x) * Legendre1(y);
    }

    static double gradleg_quad_p2p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p1_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre1(y);
    }

    static double gradleg_quad_p2p1_b2_by(double x, double y)
    {
      return l2(x) * Legendre1x(y);
    }

    static double gradleg_quad_p2p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p2_b2_b(double x, double y)
    {
      return l2(x) * Legendre2(y);
    }

    static double gradleg_quad_p2p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p2_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre2(y);
    }

    static double gradleg_quad_p2p2_b2_by(double x, double y)
    {
      return l2(x) * Legendre2x(y);
    }

    static double gradleg_quad_p2p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p3_b2_b(double x, double y)
    {
      return l2(x) * Legendre3(y);
    }

    static double gradleg_quad_p2p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p3_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre3(y);
    }

    static double gradleg_quad_p2p3_b2_by(double x, double y)
    {
      return l2(x) * Legendre3x(y);
    }

    static double gradleg_quad_p2p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p4_b2_b(double x, double y)
    {
      return l2(x) * Legendre4(y);
    }

    static double gradleg_quad_p2p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p4_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre4(y);
    }

    static double gradleg_quad_p2p4_b2_by(double x, double y)
    {
      return l2(x) * Legendre4x(y);
    }

    static double gradleg_quad_p2p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p5_b2_b(double x, double y)
    {
      return l2(x) * Legendre5(y);
    }

    static double gradleg_quad_p2p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p5_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre5(y);
    }

    static double gradleg_quad_p2p5_b2_by(double x, double y)
    {
      return l2(x) * Legendre5x(y);
    }

    static double gradleg_quad_p2p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p6_b2_b(double x, double y)
    {
      return l2(x) * Legendre6(y);
    }

    static double gradleg_quad_p2p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p6_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre6(y);
    }

    static double gradleg_quad_p2p6_b2_by(double x, double y)
    {
      return l2(x) * Legendre6x(y);
    }

    static double gradleg_quad_p2p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p7_b2_b(double x, double y)
    {
      return l2(x) * Legendre7(y);
    }

    static double gradleg_quad_p2p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p7_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre7(y);
    }

    static double gradleg_quad_p2p7_b2_by(double x, double y)
    {
      return l2(x) * Legendre7x(y);
    }

    static double gradleg_quad_p2p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p8_b2_b(double x, double y)
    {
      return l2(x) * Legendre8(y);
    }

    static double gradleg_quad_p2p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p8_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre8(y);
    }

    static double gradleg_quad_p2p8_b2_by(double x, double y)
    {
      return l2(x) * Legendre8x(y);
    }

    static double gradleg_quad_p2p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p9_b2_b(double x, double y)
    {
      return l2(x) * Legendre9(y);
    }

    static double gradleg_quad_p2p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p9_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre9(y);
    }

    static double gradleg_quad_p2p9_b2_by(double x, double y)
    {
      return l2(x) * Legendre9x(y);
    }

    static double gradleg_quad_p2p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p10_b2_b(double x, double y)
    {
      return l2(x) * Legendre10(y);
    }

    static double gradleg_quad_p2p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p2p10_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre10(y);
    }

    static double gradleg_quad_p2p10_b2_by(double x, double y)
    {
      return l2(x) * Legendre10x(y);
    }

    static double gradleg_quad_p3p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p0_b2_b(double x, double y)
    {
      return l3(x) * Legendre0(y);
    }

    static double gradleg_quad_p3p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p0_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre0(y);
    }

    static double gradleg_quad_p3p0_b2_by(double x, double y)
    {
      return l3(x) * Legendre0x(y);
    }

    static double gradleg_quad_p3p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p1_b2_b(double x, double y)
    {
      return l3(x) * Legendre1(y);
    }

    static double gradleg_quad_p3p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p1_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre1(y);
    }

    static double gradleg_quad_p3p1_b2_by(double x, double y)
    {
      return l3(x) * Legendre1x(y);
    }

    static double gradleg_quad_p3p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p2_b2_b(double x, double y)
    {
      return l3(x) * Legendre2(y);
    }

    static double gradleg_quad_p3p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p2_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre2(y);
    }

    static double gradleg_quad_p3p2_b2_by(double x, double y)
    {
      return l3(x) * Legendre2x(y);
    }

    static double gradleg_quad_p3p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p3_b2_b(double x, double y)
    {
      return l3(x) * Legendre3(y);
    }

    static double gradleg_quad_p3p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p3_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre3(y);
    }

    static double gradleg_quad_p3p3_b2_by(double x, double y)
    {
      return l3(x) * Legendre3x(y);
    }

    static double gradleg_quad_p3p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p4_b2_b(double x, double y)
    {
      return l3(x) * Legendre4(y);
    }

    static double gradleg_quad_p3p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p4_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre4(y);
    }

    static double gradleg_quad_p3p4_b2_by(double x, double y)
    {
      return l3(x) * Legendre4x(y);
    }

    static double gradleg_quad_p3p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p5_b2_b(double x, double y)
    {
      return l3(x) * Legendre5(y);
    }

    static double gradleg_quad_p3p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p5_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre5(y);
    }

    static double gradleg_quad_p3p5_b2_by(double x, double y)
    {
      return l3(x) * Legendre5x(y);
    }

    static double gradleg_quad_p3p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p6_b2_b(double x, double y)
    {
      return l3(x) * Legendre6(y);
    }

    static double gradleg_quad_p3p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p6_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre6(y);
    }

    static double gradleg_quad_p3p6_b2_by(double x, double y)
    {
      return l3(x) * Legendre6x(y);
    }

    static double gradleg_quad_p3p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p7_b2_b(double x, double y)
    {
      return l3(x) * Legendre7(y);
    }

    static double gradleg_quad_p3p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p7_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre7(y);
    }

    static double gradleg_quad_p3p7_b2_by(double x, double y)
    {
      return l3(x) * Legendre7x(y);
    }

    static double gradleg_quad_p3p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p8_b2_b(double x, double y)
    {
      return l3(x) * Legendre8(y);
    }

    static double gradleg_quad_p3p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p8_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre8(y);
    }

    static double gradleg_quad_p3p8_b2_by(double x, double y)
    {
      return l3(x) * Legendre8x(y);
    }

    static double gradleg_quad_p3p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p9_b2_b(double x, double y)
    {
      return l3(x) * Legendre9(y);
    }

    static double gradleg_quad_p3p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p9_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre9(y);
    }

    static double gradleg_quad_p3p9_b2_by(double x, double y)
    {
      return l3(x) * Legendre9x(y);
    }

    static double gradleg_quad_p3p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p10_b2_b(double x, double y)
    {
      return l3(x) * Legendre10(y);
    }

    static double gradleg_quad_p3p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p3p10_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre10(y);
    }

    static double gradleg_quad_p3p10_b2_by(double x, double y)
    {
      return l3(x) * Legendre10x(y);
    }

    static double gradleg_quad_p4p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p0_b2_b(double x, double y)
    {
      return l4(x) * Legendre0(y);
    }

    static double gradleg_quad_p4p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p0_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre0(y);
    }

    static double gradleg_quad_p4p0_b2_by(double x, double y)
    {
      return l4(x) * Legendre0x(y);
    }

    static double gradleg_quad_p4p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p1_b2_b(double x, double y)
    {
      return l4(x) * Legendre1(y);
    }

    static double gradleg_quad_p4p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p1_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre1(y);
    }

    static double gradleg_quad_p4p1_b2_by(double x, double y)
    {
      return l4(x) * Legendre1x(y);
    }

    static double gradleg_quad_p4p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p2_b2_b(double x, double y)
    {
      return l4(x) * Legendre2(y);
    }

    static double gradleg_quad_p4p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p2_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre2(y);
    }

    static double gradleg_quad_p4p2_b2_by(double x, double y)
    {
      return l4(x) * Legendre2x(y);
    }

    static double gradleg_quad_p4p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p3_b2_b(double x, double y)
    {
      return l4(x) * Legendre3(y);
    }

    static double gradleg_quad_p4p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p3_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre3(y);
    }

    static double gradleg_quad_p4p3_b2_by(double x, double y)
    {
      return l4(x) * Legendre3x(y);
    }

    static double gradleg_quad_p4p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p4_b2_b(double x, double y)
    {
      return l4(x) * Legendre4(y);
    }

    static double gradleg_quad_p4p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p4_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre4(y);
    }

    static double gradleg_quad_p4p4_b2_by(double x, double y)
    {
      return l4(x) * Legendre4x(y);
    }

    static double gradleg_quad_p4p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p5_b2_b(double x, double y)
    {
      return l4(x) * Legendre5(y);
    }

    static double gradleg_quad_p4p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p5_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre5(y);
    }

    static double gradleg_quad_p4p5_b2_by(double x, double y)
    {
      return l4(x) * Legendre5x(y);
    }

    static double gradleg_quad_p4p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p6_b2_b(double x, double y)
    {
      return l4(x) * Legendre6(y);
    }

    static double gradleg_quad_p4p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p6_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre6(y);
    }

    static double gradleg_quad_p4p6_b2_by(double x, double y)
    {
      return l4(x) * Legendre6x(y);
    }

    static double gradleg_quad_p4p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p7_b2_b(double x, double y)
    {
      return l4(x) * Legendre7(y);
    }

    static double gradleg_quad_p4p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p7_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre7(y);
    }

    static double gradleg_quad_p4p7_b2_by(double x, double y)
    {
      return l4(x) * Legendre7x(y);
    }

    static double gradleg_quad_p4p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p8_b2_b(double x, double y)
    {
      return l4(x) * Legendre8(y);
    }

    static double gradleg_quad_p4p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p8_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre8(y);
    }

    static double gradleg_quad_p4p8_b2_by(double x, double y)
    {
      return l4(x) * Legendre8x(y);
    }

    static double gradleg_quad_p4p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p9_b2_b(double x, double y)
    {
      return l4(x) * Legendre9(y);
    }

    static double gradleg_quad_p4p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p9_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre9(y);
    }

    static double gradleg_quad_p4p9_b2_by(double x, double y)
    {
      return l4(x) * Legendre9x(y);
    }

    static double gradleg_quad_p4p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p10_b2_b(double x, double y)
    {
      return l4(x) * Legendre10(y);
    }

    static double gradleg_quad_p4p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p4p10_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre10(y);
    }

    static double gradleg_quad_p4p10_b2_by(double x, double y)
    {
      return l4(x) * Legendre10x(y);
    }

    static double gradleg_quad_p5p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p0_b2_b(double x, double y)
    {
      return l5(x) * Legendre0(y);
    }

    static double gradleg_quad_p5p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p0_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre0(y);
    }

    static double gradleg_quad_p5p0_b2_by(double x, double y)
    {
      return l5(x) * Legendre0x(y);
    }

    static double gradleg_quad_p5p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p1_b2_b(double x, double y)
    {
      return l5(x) * Legendre1(y);
    }

    static double gradleg_quad_p5p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p1_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre1(y);
    }

    static double gradleg_quad_p5p1_b2_by(double x, double y)
    {
      return l5(x) * Legendre1x(y);
    }

    static double gradleg_quad_p5p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p2_b2_b(double x, double y)
    {
      return l5(x) * Legendre2(y);
    }

    static double gradleg_quad_p5p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p2_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre2(y);
    }

    static double gradleg_quad_p5p2_b2_by(double x, double y)
    {
      return l5(x) * Legendre2x(y);
    }

    static double gradleg_quad_p5p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p3_b2_b(double x, double y)
    {
      return l5(x) * Legendre3(y);
    }

    static double gradleg_quad_p5p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p3_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre3(y);
    }

    static double gradleg_quad_p5p3_b2_by(double x, double y)
    {
      return l5(x) * Legendre3x(y);
    }

    static double gradleg_quad_p5p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p4_b2_b(double x, double y)
    {
      return l5(x) * Legendre4(y);
    }

    static double gradleg_quad_p5p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p4_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre4(y);
    }

    static double gradleg_quad_p5p4_b2_by(double x, double y)
    {
      return l5(x) * Legendre4x(y);
    }

    static double gradleg_quad_p5p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p5_b2_b(double x, double y)
    {
      return l5(x) * Legendre5(y);
    }

    static double gradleg_quad_p5p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p5_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre5(y);
    }

    static double gradleg_quad_p5p5_b2_by(double x, double y)
    {
      return l5(x) * Legendre5x(y);
    }

    static double gradleg_quad_p5p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p6_b2_b(double x, double y)
    {
      return l5(x) * Legendre6(y);
    }

    static double gradleg_quad_p5p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p6_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre6(y);
    }

    static double gradleg_quad_p5p6_b2_by(double x, double y)
    {
      return l5(x) * Legendre6x(y);
    }

    static double gradleg_quad_p5p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p7_b2_b(double x, double y)
    {
      return l5(x) * Legendre7(y);
    }

    static double gradleg_quad_p5p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p7_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre7(y);
    }

    static double gradleg_quad_p5p7_b2_by(double x, double y)
    {
      return l5(x) * Legendre7x(y);
    }

    static double gradleg_quad_p5p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p8_b2_b(double x, double y)
    {
      return l5(x) * Legendre8(y);
    }

    static double gradleg_quad_p5p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p8_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre8(y);
    }

    static double gradleg_quad_p5p8_b2_by(double x, double y)
    {
      return l5(x) * Legendre8x(y);
    }

    static double gradleg_quad_p5p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p9_b2_b(double x, double y)
    {
      return l5(x) * Legendre9(y);
    }

    static double gradleg_quad_p5p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p9_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre9(y);
    }

    static double gradleg_quad_p5p9_b2_by(double x, double y)
    {
      return l5(x) * Legendre9x(y);
    }

    static double gradleg_quad_p5p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p10_b2_b(double x, double y)
    {
      return l5(x) * Legendre10(y);
    }

    static double gradleg_quad_p5p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p5p10_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre10(y);
    }

    static double gradleg_quad_p5p10_b2_by(double x, double y)
    {
      return l5(x) * Legendre10x(y);
    }

    static double gradleg_quad_p6p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p0_b2_b(double x, double y)
    {
      return l6(x) * Legendre0(y);
    }

    static double gradleg_quad_p6p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p0_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre0(y);
    }

    static double gradleg_quad_p6p0_b2_by(double x, double y)
    {
      return l6(x) * Legendre0x(y);
    }

    static double gradleg_quad_p6p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p1_b2_b(double x, double y)
    {
      return l6(x) * Legendre1(y);
    }

    static double gradleg_quad_p6p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p1_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre1(y);
    }

    static double gradleg_quad_p6p1_b2_by(double x, double y)
    {
      return l6(x) * Legendre1x(y);
    }

    static double gradleg_quad_p6p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p2_b2_b(double x, double y)
    {
      return l6(x) * Legendre2(y);
    }

    static double gradleg_quad_p6p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p2_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre2(y);
    }

    static double gradleg_quad_p6p2_b2_by(double x, double y)
    {
      return l6(x) * Legendre2x(y);
    }

    static double gradleg_quad_p6p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p3_b2_b(double x, double y)
    {
      return l6(x) * Legendre3(y);
    }

    static double gradleg_quad_p6p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p3_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre3(y);
    }

    static double gradleg_quad_p6p3_b2_by(double x, double y)
    {
      return l6(x) * Legendre3x(y);
    }

    static double gradleg_quad_p6p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p4_b2_b(double x, double y)
    {
      return l6(x) * Legendre4(y);
    }

    static double gradleg_quad_p6p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p4_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre4(y);
    }

    static double gradleg_quad_p6p4_b2_by(double x, double y)
    {
      return l6(x) * Legendre4x(y);
    }

    static double gradleg_quad_p6p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p5_b2_b(double x, double y)
    {
      return l6(x) * Legendre5(y);
    }

    static double gradleg_quad_p6p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p5_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre5(y);
    }

    static double gradleg_quad_p6p5_b2_by(double x, double y)
    {
      return l6(x) * Legendre5x(y);
    }

    static double gradleg_quad_p6p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p6_b2_b(double x, double y)
    {
      return l6(x) * Legendre6(y);
    }

    static double gradleg_quad_p6p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p6_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre6(y);
    }

    static double gradleg_quad_p6p6_b2_by(double x, double y)
    {
      return l6(x) * Legendre6x(y);
    }

    static double gradleg_quad_p6p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p7_b2_b(double x, double y)
    {
      return l6(x) * Legendre7(y);
    }

    static double gradleg_quad_p6p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p7_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre7(y);
    }

    static double gradleg_quad_p6p7_b2_by(double x, double y)
    {
      return l6(x) * Legendre7x(y);
    }

    static double gradleg_quad_p6p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p8_b2_b(double x, double y)
    {
      return l6(x) * Legendre8(y);
    }

    static double gradleg_quad_p6p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p8_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre8(y);
    }

    static double gradleg_quad_p6p8_b2_by(double x, double y)
    {
      return l6(x) * Legendre8x(y);
    }

    static double gradleg_quad_p6p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p9_b2_b(double x, double y)
    {
      return l6(x) * Legendre9(y);
    }

    static double gradleg_quad_p6p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p9_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre9(y);
    }

    static double gradleg_quad_p6p9_b2_by(double x, double y)
    {
      return l6(x) * Legendre9x(y);
    }

    static double gradleg_quad_p6p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p10_b2_b(double x, double y)
    {
      return l6(x) * Legendre10(y);
    }

    static double gradleg_quad_p6p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p6p10_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre10(y);
    }

    static double gradleg_quad_p6p10_b2_by(double x, double y)
    {
      return l6(x) * Legendre10x(y);
    }

    static double gradleg_quad_p7p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p0_b2_b(double x, double y)
    {
      return l7(x) * Legendre0(y);
    }

    static double gradleg_quad_p7p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p0_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre0(y);
    }

    static double gradleg_quad_p7p0_b2_by(double x, double y)
    {
      return l7(x) * Legendre0x(y);
    }

    static double gradleg_quad_p7p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p1_b2_b(double x, double y)
    {
      return l7(x) * Legendre1(y);
    }

    static double gradleg_quad_p7p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p1_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre1(y);
    }

    static double gradleg_quad_p7p1_b2_by(double x, double y)
    {
      return l7(x) * Legendre1x(y);
    }

    static double gradleg_quad_p7p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p2_b2_b(double x, double y)
    {
      return l7(x) * Legendre2(y);
    }

    static double gradleg_quad_p7p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p2_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre2(y);
    }

    static double gradleg_quad_p7p2_b2_by(double x, double y)
    {
      return l7(x) * Legendre2x(y);
    }

    static double gradleg_quad_p7p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p3_b2_b(double x, double y)
    {
      return l7(x) * Legendre3(y);
    }

    static double gradleg_quad_p7p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p3_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre3(y);
    }

    static double gradleg_quad_p7p3_b2_by(double x, double y)
    {
      return l7(x) * Legendre3x(y);
    }

    static double gradleg_quad_p7p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p4_b2_b(double x, double y)
    {
      return l7(x) * Legendre4(y);
    }

    static double gradleg_quad_p7p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p4_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre4(y);
    }

    static double gradleg_quad_p7p4_b2_by(double x, double y)
    {
      return l7(x) * Legendre4x(y);
    }

    static double gradleg_quad_p7p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p5_b2_b(double x, double y)
    {
      return l7(x) * Legendre5(y);
    }

    static double gradleg_quad_p7p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p5_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre5(y);
    }

    static double gradleg_quad_p7p5_b2_by(double x, double y)
    {
      return l7(x) * Legendre5x(y);
    }

    static double gradleg_quad_p7p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p6_b2_b(double x, double y)
    {
      return l7(x) * Legendre6(y);
    }

    static double gradleg_quad_p7p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p6_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre6(y);
    }

    static double gradleg_quad_p7p6_b2_by(double x, double y)
    {
      return l7(x) * Legendre6x(y);
    }

    static double gradleg_quad_p7p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p7_b2_b(double x, double y)
    {
      return l7(x) * Legendre7(y);
    }

    static double gradleg_quad_p7p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p7_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre7(y);
    }

    static double gradleg_quad_p7p7_b2_by(double x, double y)
    {
      return l7(x) * Legendre7x(y);
    }

    static double gradleg_quad_p7p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p8_b2_b(double x, double y)
    {
      return l7(x) * Legendre8(y);
    }

    static double gradleg_quad_p7p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p8_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre8(y);
    }

    static double gradleg_quad_p7p8_b2_by(double x, double y)
    {
      return l7(x) * Legendre8x(y);
    }

    static double gradleg_quad_p7p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p9_b2_b(double x, double y)
    {
      return l7(x) * Legendre9(y);
    }

    static double gradleg_quad_p7p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p9_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre9(y);
    }

    static double gradleg_quad_p7p9_b2_by(double x, double y)
    {
      return l7(x) * Legendre9x(y);
    }

    static double gradleg_quad_p7p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p10_b2_b(double x, double y)
    {
      return l7(x) * Legendre10(y);
    }

    static double gradleg_quad_p7p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p7p10_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre10(y);
    }

    static double gradleg_quad_p7p10_b2_by(double x, double y)
    {
      return l7(x) * Legendre10x(y);
    }

    static double gradleg_quad_p8p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p0_b2_b(double x, double y)
    {
      return l8(x) * Legendre0(y);
    }

    static double gradleg_quad_p8p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p0_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre0(y);
    }

    static double gradleg_quad_p8p0_b2_by(double x, double y)
    {
      return l8(x) * Legendre0x(y);
    }

    static double gradleg_quad_p8p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p1_b2_b(double x, double y)
    {
      return l8(x) * Legendre1(y);
    }

    static double gradleg_quad_p8p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p1_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre1(y);
    }

    static double gradleg_quad_p8p1_b2_by(double x, double y)
    {
      return l8(x) * Legendre1x(y);
    }

    static double gradleg_quad_p8p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p2_b2_b(double x, double y)
    {
      return l8(x) * Legendre2(y);
    }

    static double gradleg_quad_p8p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p2_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre2(y);
    }

    static double gradleg_quad_p8p2_b2_by(double x, double y)
    {
      return l8(x) * Legendre2x(y);
    }

    static double gradleg_quad_p8p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p3_b2_b(double x, double y)
    {
      return l8(x) * Legendre3(y);
    }

    static double gradleg_quad_p8p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p3_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre3(y);
    }

    static double gradleg_quad_p8p3_b2_by(double x, double y)
    {
      return l8(x) * Legendre3x(y);
    }

    static double gradleg_quad_p8p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p4_b2_b(double x, double y)
    {
      return l8(x) * Legendre4(y);
    }

    static double gradleg_quad_p8p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p4_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre4(y);
    }

    static double gradleg_quad_p8p4_b2_by(double x, double y)
    {
      return l8(x) * Legendre4x(y);
    }

    static double gradleg_quad_p8p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p5_b2_b(double x, double y)
    {
      return l8(x) * Legendre5(y);
    }

    static double gradleg_quad_p8p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p5_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre5(y);
    }

    static double gradleg_quad_p8p5_b2_by(double x, double y)
    {
      return l8(x) * Legendre5x(y);
    }

    static double gradleg_quad_p8p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p6_b2_b(double x, double y)
    {
      return l8(x) * Legendre6(y);
    }

    static double gradleg_quad_p8p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p6_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre6(y);
    }

    static double gradleg_quad_p8p6_b2_by(double x, double y)
    {
      return l8(x) * Legendre6x(y);
    }

    static double gradleg_quad_p8p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p7_b2_b(double x, double y)
    {
      return l8(x) * Legendre7(y);
    }

    static double gradleg_quad_p8p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p7_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre7(y);
    }

    static double gradleg_quad_p8p7_b2_by(double x, double y)
    {
      return l8(x) * Legendre7x(y);
    }

    static double gradleg_quad_p8p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p8_b2_b(double x, double y)
    {
      return l8(x) * Legendre8(y);
    }

    static double gradleg_quad_p8p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p8_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre8(y);
    }

    static double gradleg_quad_p8p8_b2_by(double x, double y)
    {
      return l8(x) * Legendre8x(y);
    }

    static double gradleg_quad_p8p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p9_b2_b(double x, double y)
    {
      return l8(x) * Legendre9(y);
    }

    static double gradleg_quad_p8p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p9_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre9(y);
    }

    static double gradleg_quad_p8p9_b2_by(double x, double y)
    {
      return l8(x) * Legendre9x(y);
    }

    static double gradleg_quad_p8p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p10_b2_b(double x, double y)
    {
      return l8(x) * Legendre10(y);
    }

    static double gradleg_quad_p8p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p8p10_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre10(y);
    }

    static double gradleg_quad_p8p10_b2_by(double x, double y)
    {
      return l8(x) * Legendre10x(y);
    }

    static double gradleg_quad_p9p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p0_b2_b(double x, double y)
    {
      return l9(x) * Legendre0(y);
    }

    static double gradleg_quad_p9p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p0_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre0(y);
    }

    static double gradleg_quad_p9p0_b2_by(double x, double y)
    {
      return l9(x) * Legendre0x(y);
    }

    static double gradleg_quad_p9p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p1_b2_b(double x, double y)
    {
      return l9(x) * Legendre1(y);
    }

    static double gradleg_quad_p9p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p1_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre1(y);
    }

    static double gradleg_quad_p9p1_b2_by(double x, double y)
    {
      return l9(x) * Legendre1x(y);
    }

    static double gradleg_quad_p9p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p2_b2_b(double x, double y)
    {
      return l9(x) * Legendre2(y);
    }

    static double gradleg_quad_p9p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p2_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre2(y);
    }

    static double gradleg_quad_p9p2_b2_by(double x, double y)
    {
      return l9(x) * Legendre2x(y);
    }

    static double gradleg_quad_p9p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p3_b2_b(double x, double y)
    {
      return l9(x) * Legendre3(y);
    }

    static double gradleg_quad_p9p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p3_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre3(y);
    }

    static double gradleg_quad_p9p3_b2_by(double x, double y)
    {
      return l9(x) * Legendre3x(y);
    }

    static double gradleg_quad_p9p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p4_b2_b(double x, double y)
    {
      return l9(x) * Legendre4(y);
    }

    static double gradleg_quad_p9p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p4_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre4(y);
    }

    static double gradleg_quad_p9p4_b2_by(double x, double y)
    {
      return l9(x) * Legendre4x(y);
    }

    static double gradleg_quad_p9p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p5_b2_b(double x, double y)
    {
      return l9(x) * Legendre5(y);
    }

    static double gradleg_quad_p9p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p5_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre5(y);
    }

    static double gradleg_quad_p9p5_b2_by(double x, double y)
    {
      return l9(x) * Legendre5x(y);
    }

    static double gradleg_quad_p9p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p6_b2_b(double x, double y)
    {
      return l9(x) * Legendre6(y);
    }

    static double gradleg_quad_p9p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p6_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre6(y);
    }

    static double gradleg_quad_p9p6_b2_by(double x, double y)
    {
      return l9(x) * Legendre6x(y);
    }

    static double gradleg_quad_p9p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p7_b2_b(double x, double y)
    {
      return l9(x) * Legendre7(y);
    }

    static double gradleg_quad_p9p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p7_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre7(y);
    }

    static double gradleg_quad_p9p7_b2_by(double x, double y)
    {
      return l9(x) * Legendre7x(y);
    }

    static double gradleg_quad_p9p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p8_b2_b(double x, double y)
    {
      return l9(x) * Legendre8(y);
    }

    static double gradleg_quad_p9p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p8_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre8(y);
    }

    static double gradleg_quad_p9p8_b2_by(double x, double y)
    {
      return l9(x) * Legendre8x(y);
    }

    static double gradleg_quad_p9p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p9_b2_b(double x, double y)
    {
      return l9(x) * Legendre9(y);
    }

    static double gradleg_quad_p9p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p9_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre9(y);
    }

    static double gradleg_quad_p9p9_b2_by(double x, double y)
    {
      return l9(x) * Legendre9x(y);
    }

    static double gradleg_quad_p9p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p10_b2_b(double x, double y)
    {
      return l9(x) * Legendre10(y);
    }

    static double gradleg_quad_p9p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p9p10_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre10(y);
    }

    static double gradleg_quad_p9p10_b2_by(double x, double y)
    {
      return l9(x) * Legendre10x(y);
    }

    static double gradleg_quad_p10p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p0_b2_b(double x, double y)
    {
      return l10(x) * Legendre0(y);
    }

    static double gradleg_quad_p10p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p0_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre0(y);
    }

    static double gradleg_quad_p10p0_b2_by(double x, double y)
    {
      return l10(x) * Legendre0x(y);
    }

    static double gradleg_quad_p10p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p1_b2_b(double x, double y)
    {
      return l10(x) * Legendre1(y);
    }

    static double gradleg_quad_p10p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p1_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre1(y);
    }

    static double gradleg_quad_p10p1_b2_by(double x, double y)
    {
      return l10(x) * Legendre1x(y);
    }

    static double gradleg_quad_p10p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p2_b2_b(double x, double y)
    {
      return l10(x) * Legendre2(y);
    }

    static double gradleg_quad_p10p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p2_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre2(y);
    }

    static double gradleg_quad_p10p2_b2_by(double x, double y)
    {
      return l10(x) * Legendre2x(y);
    }

    static double gradleg_quad_p10p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p3_b2_b(double x, double y)
    {
      return l10(x) * Legendre3(y);
    }

    static double gradleg_quad_p10p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p3_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre3(y);
    }

    static double gradleg_quad_p10p3_b2_by(double x, double y)
    {
      return l10(x) * Legendre3x(y);
    }

    static double gradleg_quad_p10p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p4_b2_b(double x, double y)
    {
      return l10(x) * Legendre4(y);
    }

    static double gradleg_quad_p10p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p4_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre4(y);
    }

    static double gradleg_quad_p10p4_b2_by(double x, double y)
    {
      return l10(x) * Legendre4x(y);
    }

    static double gradleg_quad_p10p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p5_b2_b(double x, double y)
    {
      return l10(x) * Legendre5(y);
    }

    static double gradleg_quad_p10p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p5_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre5(y);
    }

    static double gradleg_quad_p10p5_b2_by(double x, double y)
    {
      return l10(x) * Legendre5x(y);
    }

    static double gradleg_quad_p10p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p6_b2_b(double x, double y)
    {
      return l10(x) * Legendre6(y);
    }

    static double gradleg_quad_p10p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p6_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre6(y);
    }

    static double gradleg_quad_p10p6_b2_by(double x, double y)
    {
      return l10(x) * Legendre6x(y);
    }

    static double gradleg_quad_p10p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p7_b2_b(double x, double y)
    {
      return l10(x) * Legendre7(y);
    }

    static double gradleg_quad_p10p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p7_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre7(y);
    }

    static double gradleg_quad_p10p7_b2_by(double x, double y)
    {
      return l10(x) * Legendre7x(y);
    }

    static double gradleg_quad_p10p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p8_b2_b(double x, double y)
    {
      return l10(x) * Legendre8(y);
    }

    static double gradleg_quad_p10p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p8_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre8(y);
    }

    static double gradleg_quad_p10p8_b2_by(double x, double y)
    {
      return l10(x) * Legendre8x(y);
    }

    static double gradleg_quad_p10p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p9_b2_b(double x, double y)
    {
      return l10(x) * Legendre9(y);
    }

    static double gradleg_quad_p10p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p9_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre9(y);
    }

    static double gradleg_quad_p10p9_b2_by(double x, double y)
    {
      return l10(x) * Legendre9x(y);
    }

    static double gradleg_quad_p10p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p10_b2_b(double x, double y)
    {
      return l10(x) * Legendre10(y);
    }

    static double gradleg_quad_p10p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p10p10_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre10(y);
    }

    static double gradleg_quad_p10p10_b2_by(double x, double y)
    {
      return l10(x) * Legendre10x(y);
    }

    static double gradleg_quad_p11p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p0_b2_b(double x, double y)
    {
      return l11(x) * Legendre0(y);
    }

    static double gradleg_quad_p11p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p0_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre0(y);
    }

    static double gradleg_quad_p11p0_b2_by(double x, double y)
    {
      return l11(x) * Legendre0x(y);
    }

    static double gradleg_quad_p11p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p1_b2_b(double x, double y)
    {
      return l11(x) * Legendre1(y);
    }

    static double gradleg_quad_p11p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p1_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre1(y);
    }

    static double gradleg_quad_p11p1_b2_by(double x, double y)
    {
      return l11(x) * Legendre1x(y);
    }

    static double gradleg_quad_p11p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p2_b2_b(double x, double y)
    {
      return l11(x) * Legendre2(y);
    }

    static double gradleg_quad_p11p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p2_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre2(y);
    }

    static double gradleg_quad_p11p2_b2_by(double x, double y)
    {
      return l11(x) * Legendre2x(y);
    }

    static double gradleg_quad_p11p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p3_b2_b(double x, double y)
    {
      return l11(x) * Legendre3(y);
    }

    static double gradleg_quad_p11p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p3_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre3(y);
    }

    static double gradleg_quad_p11p3_b2_by(double x, double y)
    {
      return l11(x) * Legendre3x(y);
    }

    static double gradleg_quad_p11p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p4_b2_b(double x, double y)
    {
      return l11(x) * Legendre4(y);
    }

    static double gradleg_quad_p11p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p4_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre4(y);
    }

    static double gradleg_quad_p11p4_b2_by(double x, double y)
    {
      return l11(x) * Legendre4x(y);
    }

    static double gradleg_quad_p11p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p5_b2_b(double x, double y)
    {
      return l11(x) * Legendre5(y);
    }

    static double gradleg_quad_p11p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p5_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre5(y);
    }

    static double gradleg_quad_p11p5_b2_by(double x, double y)
    {
      return l11(x) * Legendre5x(y);
    }

    static double gradleg_quad_p11p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p6_b2_b(double x, double y)
    {
      return l11(x) * Legendre6(y);
    }

    static double gradleg_quad_p11p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p6_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre6(y);
    }

    static double gradleg_quad_p11p6_b2_by(double x, double y)
    {
      return l11(x) * Legendre6x(y);
    }

    static double gradleg_quad_p11p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p7_b2_b(double x, double y)
    {
      return l11(x) * Legendre7(y);
    }

    static double gradleg_quad_p11p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p7_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre7(y);
    }

    static double gradleg_quad_p11p7_b2_by(double x, double y)
    {
      return l11(x) * Legendre7x(y);
    }

    static double gradleg_quad_p11p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p8_b2_b(double x, double y)
    {
      return l11(x) * Legendre8(y);
    }

    static double gradleg_quad_p11p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p8_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre8(y);
    }

    static double gradleg_quad_p11p8_b2_by(double x, double y)
    {
      return l11(x) * Legendre8x(y);
    }

    static double gradleg_quad_p11p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p9_b2_b(double x, double y)
    {
      return l11(x) * Legendre9(y);
    }

    static double gradleg_quad_p11p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p9_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre9(y);
    }

    static double gradleg_quad_p11p9_b2_by(double x, double y)
    {
      return l11(x) * Legendre9x(y);
    }

    static double gradleg_quad_p11p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p10_b2_b(double x, double y)
    {
      return l11(x) * Legendre10(y);
    }

    static double gradleg_quad_p11p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double gradleg_quad_p11p10_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre10(y);
    }

    static double gradleg_quad_p11p10_b2_by(double x, double y)
    {
      return l11(x) * Legendre10x(y);
    }

    static Shapeset::shape_fn_t gradleg_quad_fn_a[] =
    {
      gradleg_quad_p0_e1_a_0, gradleg_quad_p0_e1_a_1, gradleg_quad_p0_e2_a, gradleg_quad_p0_e2_a, gradleg_quad_p0_e3_a_0, gradleg_quad_p0_e3_a_1, gradleg_quad_p0_e4_a, gradleg_quad_p0_e4_a,
      gradleg_quad_l2_l0_a, gradleg_quad_l2_l0_a, gradleg_quad_l1_l2_a, gradleg_quad_l1_l2_a, gradleg_quad_l2_l1_a, gradleg_quad_l2_l1_a, gradleg_quad_l0_l2_a, gradleg_quad_l0_l2_a,
      gradleg_quad_l3_l0_a_0, gradleg_quad_l3_l0_a_1, gradleg_quad_l1_l3_a_0, gradleg_quad_l1_l3_a_1,  gradleg_quad_l3_l1_a_0, gradleg_quad_l3_l1_a_1, gradleg_quad_l0_l3_a_0, gradleg_quad_l0_l3_a_1,
      gradleg_quad_l4_l0_a, gradleg_quad_l4_l0_a, gradleg_quad_l1_l4_a, gradleg_quad_l1_l4_a, gradleg_quad_l4_l1_a, gradleg_quad_l4_l1_a, gradleg_quad_l0_l4_a, gradleg_quad_l0_l4_a,
      gradleg_quad_l5_l0_a_0, gradleg_quad_l5_l0_a_1, gradleg_quad_l1_l5_a_0, gradleg_quad_l1_l5_a_1,  gradleg_quad_l5_l1_a_0, gradleg_quad_l5_l1_a_1, gradleg_quad_l0_l5_a_0, gradleg_quad_l0_l5_a_1,
      gradleg_quad_l6_l0_a, gradleg_quad_l6_l0_a, gradleg_quad_l1_l6_a, gradleg_quad_l1_l6_a, gradleg_quad_l6_l1_a, gradleg_quad_l6_l1_a, gradleg_quad_l0_l6_a, gradleg_quad_l0_l6_a,
      gradleg_quad_l7_l0_a_0, gradleg_quad_l7_l0_a_1, gradleg_quad_l1_l7_a_0, gradleg_quad_l1_l7_a_1,  gradleg_quad_l7_l1_a_0, gradleg_quad_l7_l1_a_1, gradleg_quad_l0_l7_a_0, gradleg_quad_l0_l7_a_1,
      gradleg_quad_l8_l0_a, gradleg_quad_l8_l0_a, gradleg_quad_l1_l8_a, gradleg_quad_l1_l8_a, gradleg_quad_l8_l1_a, gradleg_quad_l8_l1_a, gradleg_quad_l0_l8_a, gradleg_quad_l0_l8_a,
      gradleg_quad_l9_l0_a_0, gradleg_quad_l9_l0_a_1, gradleg_quad_l1_l9_a_0, gradleg_quad_l1_l9_a_1,  gradleg_quad_l9_l1_a_0, gradleg_quad_l9_l1_a_1, gradleg_quad_l0_l9_a_0, gradleg_quad_l0_l9_a_1,
      gradleg_quad_l10_l0_a, gradleg_quad_l10_l0_a, gradleg_quad_l1_l10_a, gradleg_quad_l1_l10_a, gradleg_quad_l10_l1_a, gradleg_quad_l10_l1_a, gradleg_quad_l0_l10_a, gradleg_quad_l0_l10_a,
      gradleg_quad_l11_l0_a_0, gradleg_quad_l11_l0_a_1, gradleg_quad_l1_l11_a_0, gradleg_quad_l1_l11_a_1,  gradleg_quad_l11_l1_a_0, gradleg_quad_l11_l1_a_1, gradleg_quad_l0_l11_a_0, gradleg_quad_l0_l11_a_1,

      gradleg_quad_p0p2_b1_a,   gradleg_quad_p0p3_b1_a,   gradleg_quad_p0p4_b1_a,   gradleg_quad_p0p5_b1_a,   gradleg_quad_p0p6_b1_a,   gradleg_quad_p0p7_b1_a,   gradleg_quad_p0p8_b1_a,   gradleg_quad_p0p9_b1_a,   gradleg_quad_p0p10_b1_a,   gradleg_quad_p0p11_b1_a,   gradleg_quad_p1p2_b1_a,   gradleg_quad_p1p3_b1_a,   gradleg_quad_p1p4_b1_a,   gradleg_quad_p1p5_b1_a,   gradleg_quad_p1p6_b1_a,   gradleg_quad_p1p7_b1_a,   gradleg_quad_p1p8_b1_a,   gradleg_quad_p1p9_b1_a,   gradleg_quad_p1p10_b1_a,   gradleg_quad_p1p11_b1_a,   gradleg_quad_p2p2_b1_a,   gradleg_quad_p2p3_b1_a,   gradleg_quad_p2p4_b1_a,   gradleg_quad_p2p5_b1_a,   gradleg_quad_p2p6_b1_a,   gradleg_quad_p2p7_b1_a,   gradleg_quad_p2p8_b1_a,   gradleg_quad_p2p9_b1_a,   gradleg_quad_p2p10_b1_a,   gradleg_quad_p2p11_b1_a,   gradleg_quad_p3p2_b1_a,   gradleg_quad_p3p3_b1_a,   gradleg_quad_p3p4_b1_a,   gradleg_quad_p3p5_b1_a,   gradleg_quad_p3p6_b1_a,   gradleg_quad_p3p7_b1_a,   gradleg_quad_p3p8_b1_a,   gradleg_quad_p3p9_b1_a,   gradleg_quad_p3p10_b1_a,   gradleg_quad_p3p11_b1_a,   gradleg_quad_p4p2_b1_a,   gradleg_quad_p4p3_b1_a,   gradleg_quad_p4p4_b1_a,   gradleg_quad_p4p5_b1_a,   gradleg_quad_p4p6_b1_a,   gradleg_quad_p4p7_b1_a,   gradleg_quad_p4p8_b1_a,   gradleg_quad_p4p9_b1_a,   gradleg_quad_p4p10_b1_a,   gradleg_quad_p4p11_b1_a,   gradleg_quad_p5p2_b1_a,   gradleg_quad_p5p3_b1_a,   gradleg_quad_p5p4_b1_a,   gradleg_quad_p5p5_b1_a,   gradleg_quad_p5p6_b1_a,   gradleg_quad_p5p7_b1_a,   gradleg_quad_p5p8_b1_a,   gradleg_quad_p5p9_b1_a,   gradleg_quad_p5p10_b1_a,   gradleg_quad_p5p11_b1_a,   gradleg_quad_p6p2_b1_a,   gradleg_quad_p6p3_b1_a,   gradleg_quad_p6p4_b1_a,   gradleg_quad_p6p5_b1_a,   gradleg_quad_p6p6_b1_a,   gradleg_quad_p6p7_b1_a,   gradleg_quad_p6p8_b1_a,   gradleg_quad_p6p9_b1_a,   gradleg_quad_p6p10_b1_a,   gradleg_quad_p6p11_b1_a,   gradleg_quad_p7p2_b1_a,   gradleg_quad_p7p3_b1_a,   gradleg_quad_p7p4_b1_a,   gradleg_quad_p7p5_b1_a,   gradleg_quad_p7p6_b1_a,   gradleg_quad_p7p7_b1_a,   gradleg_quad_p7p8_b1_a,   gradleg_quad_p7p9_b1_a,   gradleg_quad_p7p10_b1_a,   gradleg_quad_p7p11_b1_a,   gradleg_quad_p8p2_b1_a,   gradleg_quad_p8p3_b1_a,   gradleg_quad_p8p4_b1_a,   gradleg_quad_p8p5_b1_a,   gradleg_quad_p8p6_b1_a,   gradleg_quad_p8p7_b1_a,   gradleg_quad_p8p8_b1_a,   gradleg_quad_p8p9_b1_a,   gradleg_quad_p8p10_b1_a,   gradleg_quad_p8p11_b1_a,   gradleg_quad_p9p2_b1_a,   gradleg_quad_p9p3_b1_a,   gradleg_quad_p9p4_b1_a,   gradleg_quad_p9p5_b1_a,   gradleg_quad_p9p6_b1_a,   gradleg_quad_p9p7_b1_a,   gradleg_quad_p9p8_b1_a,   gradleg_quad_p9p9_b1_a,   gradleg_quad_p9p10_b1_a,   gradleg_quad_p9p11_b1_a,   gradleg_quad_p10p2_b1_a,   gradleg_quad_p10p3_b1_a,   gradleg_quad_p10p4_b1_a,   gradleg_quad_p10p5_b1_a,   gradleg_quad_p10p6_b1_a,   gradleg_quad_p10p7_b1_a,   gradleg_quad_p10p8_b1_a,   gradleg_quad_p10p9_b1_a,   gradleg_quad_p10p10_b1_a,   gradleg_quad_p10p11_b1_a,   gradleg_quad_p2p0_b2_a,   gradleg_quad_p2p1_b2_a,   gradleg_quad_p2p2_b2_a,   gradleg_quad_p2p3_b2_a,   gradleg_quad_p2p4_b2_a,   gradleg_quad_p2p5_b2_a,   gradleg_quad_p2p6_b2_a,   gradleg_quad_p2p7_b2_a,   gradleg_quad_p2p8_b2_a,   gradleg_quad_p2p9_b2_a,   gradleg_quad_p2p10_b2_a,   gradleg_quad_p3p0_b2_a,   gradleg_quad_p3p1_b2_a,   gradleg_quad_p3p2_b2_a,   gradleg_quad_p3p3_b2_a,   gradleg_quad_p3p4_b2_a,   gradleg_quad_p3p5_b2_a,   gradleg_quad_p3p6_b2_a,   gradleg_quad_p3p7_b2_a,   gradleg_quad_p3p8_b2_a,   gradleg_quad_p3p9_b2_a,   gradleg_quad_p3p10_b2_a,   gradleg_quad_p4p0_b2_a,   gradleg_quad_p4p1_b2_a,   gradleg_quad_p4p2_b2_a,   gradleg_quad_p4p3_b2_a,   gradleg_quad_p4p4_b2_a,   gradleg_quad_p4p5_b2_a,   gradleg_quad_p4p6_b2_a,   gradleg_quad_p4p7_b2_a,   gradleg_quad_p4p8_b2_a,   gradleg_quad_p4p9_b2_a,   gradleg_quad_p4p10_b2_a,   gradleg_quad_p5p0_b2_a,   gradleg_quad_p5p1_b2_a,   gradleg_quad_p5p2_b2_a,   gradleg_quad_p5p3_b2_a,   gradleg_quad_p5p4_b2_a,   gradleg_quad_p5p5_b2_a,   gradleg_quad_p5p6_b2_a,   gradleg_quad_p5p7_b2_a,   gradleg_quad_p5p8_b2_a,   gradleg_quad_p5p9_b2_a,   gradleg_quad_p5p10_b2_a,   gradleg_quad_p6p0_b2_a,   gradleg_quad_p6p1_b2_a,   gradleg_quad_p6p2_b2_a,   gradleg_quad_p6p3_b2_a,   gradleg_quad_p6p4_b2_a,   gradleg_quad_p6p5_b2_a,   gradleg_quad_p6p6_b2_a,   gradleg_quad_p6p7_b2_a,   gradleg_quad_p6p8_b2_a,   gradleg_quad_p6p9_b2_a,   gradleg_quad_p6p10_b2_a,   gradleg_quad_p7p0_b2_a,   gradleg_quad_p7p1_b2_a,   gradleg_quad_p7p2_b2_a,   gradleg_quad_p7p3_b2_a,   gradleg_quad_p7p4_b2_a,   gradleg_quad_p7p5_b2_a,   gradleg_quad_p7p6_b2_a,   gradleg_quad_p7p7_b2_a,   gradleg_quad_p7p8_b2_a,   gradleg_quad_p7p9_b2_a,   gradleg_quad_p7p10_b2_a,   gradleg_quad_p8p0_b2_a,   gradleg_quad_p8p1_b2_a,   gradleg_quad_p8p2_b2_a,   gradleg_quad_p8p3_b2_a,   gradleg_quad_p8p4_b2_a,   gradleg_quad_p8p5_b2_a,   gradleg_quad_p8p6_b2_a,   gradleg_quad_p8p7_b2_a,   gradleg_quad_p8p8_b2_a,   gradleg_quad_p8p9_b2_a,   gradleg_quad_p8p10_b2_a,   gradleg_quad_p9p0_b2_a,   gradleg_quad_p9p1_b2_a,   gradleg_quad_p9p2_b2_a,   gradleg_quad_p9p3_b2_a,   gradleg_quad_p9p4_b2_a,   gradleg_quad_p9p5_b2_a,   gradleg_quad_p9p6_b2_a,   gradleg_quad_p9p7_b2_a,   gradleg_quad_p9p8_b2_a,   gradleg_quad_p9p9_b2_a,   gradleg_quad_p9p10_b2_a,   gradleg_quad_p10p0_b2_a,   gradleg_quad_p10p1_b2_a,   gradleg_quad_p10p2_b2_a,   gradleg_quad_p10p3_b2_a,   gradleg_quad_p10p4_b2_a,   gradleg_quad_p10p5_b2_a,   gradleg_quad_p10p6_b2_a,   gradleg_quad_p10p7_b2_a,   gradleg_quad_p10p8_b2_a,   gradleg_quad_p10p9_b2_a,   gradleg_quad_p10p10_b2_a,   gradleg_quad_p11p0_b2_a,   gradleg_quad_p11p1_b2_a,   gradleg_quad_p11p2_b2_a,   gradleg_quad_p11p3_b2_a,   gradleg_quad_p11p4_b2_a,   gradleg_quad_p11p5_b2_a,   gradleg_quad_p11p6_b2_a,   gradleg_quad_p11p7_b2_a,   gradleg_quad_p11p8_b2_a,   gradleg_quad_p11p9_b2_a,   gradleg_quad_p11p10_b2_a, };

    static Shapeset::shape_fn_t gradleg_quad_fn_b[] =
    {
      gradleg_quad_p0_e1_b, gradleg_quad_p0_e1_b, gradleg_quad_p0_e2_b_0, gradleg_quad_p0_e2_b_1,  gradleg_quad_p0_e3_b, gradleg_quad_p0_e3_b, gradleg_quad_p0_e4_b_0, gradleg_quad_p0_e4_b_1,
      gradleg_quad_l2_l0_b, gradleg_quad_l2_l0_b, gradleg_quad_l1_l2_b, gradleg_quad_l1_l2_b, gradleg_quad_l2_l1_b, gradleg_quad_l2_l1_b, gradleg_quad_l0_l2_b,  gradleg_quad_l0_l2_b,
      gradleg_quad_l3_l0_b_0, gradleg_quad_l3_l0_b_1, gradleg_quad_l1_l3_b_0, gradleg_quad_l1_l3_b_1, gradleg_quad_l3_l1_b_0, gradleg_quad_l3_l1_b_1, gradleg_quad_l0_l3_b_0, gradleg_quad_l0_l3_b_1,
      gradleg_quad_l4_l0_b, gradleg_quad_l4_l0_b, gradleg_quad_l1_l4_b, gradleg_quad_l1_l4_b, gradleg_quad_l4_l1_b, gradleg_quad_l4_l1_b, gradleg_quad_l0_l4_b,  gradleg_quad_l0_l4_b,
      gradleg_quad_l5_l0_b_0, gradleg_quad_l5_l0_b_1, gradleg_quad_l1_l5_b_0, gradleg_quad_l1_l5_b_1, gradleg_quad_l5_l1_b_0, gradleg_quad_l5_l1_b_1, gradleg_quad_l0_l5_b_0, gradleg_quad_l0_l5_b_1,
      gradleg_quad_l6_l0_b, gradleg_quad_l6_l0_b, gradleg_quad_l1_l6_b, gradleg_quad_l1_l6_b, gradleg_quad_l6_l1_b, gradleg_quad_l6_l1_b, gradleg_quad_l0_l6_b,  gradleg_quad_l0_l6_b,
      gradleg_quad_l7_l0_b_0, gradleg_quad_l7_l0_b_1, gradleg_quad_l1_l7_b_0, gradleg_quad_l1_l7_b_1, gradleg_quad_l7_l1_b_0, gradleg_quad_l7_l1_b_1, gradleg_quad_l0_l7_b_0, gradleg_quad_l0_l7_b_1,
      gradleg_quad_l8_l0_b, gradleg_quad_l8_l0_b, gradleg_quad_l1_l8_b, gradleg_quad_l1_l8_b, gradleg_quad_l8_l1_b, gradleg_quad_l8_l1_b, gradleg_quad_l0_l8_b,  gradleg_quad_l0_l8_b,
      gradleg_quad_l9_l0_b_0, gradleg_quad_l9_l0_b_1, gradleg_quad_l1_l9_b_0, gradleg_quad_l1_l9_b_1, gradleg_quad_l9_l1_b_0, gradleg_quad_l9_l1_b_1, gradleg_quad_l0_l9_b_0, gradleg_quad_l0_l9_b_1,
      gradleg_quad_l10_l0_b, gradleg_quad_l10_l0_b, gradleg_quad_l1_l10_b, gradleg_quad_l1_l10_b, gradleg_quad_l10_l1_b, gradleg_quad_l10_l1_b, gradleg_quad_l0_l10_b,  gradleg_quad_l0_l10_b,
      gradleg_quad_l11_l0_b_0, gradleg_quad_l11_l0_b_1, gradleg_quad_l1_l11_b_0, gradleg_quad_l1_l11_b_1, gradleg_quad_l11_l1_b_0, gradleg_quad_l11_l1_b_1, gradleg_quad_l0_l11_b_0, gradleg_quad_l0_l11_b_1,

      gradleg_quad_p0p2_b1_b,   gradleg_quad_p0p3_b1_b,   gradleg_quad_p0p4_b1_b,   gradleg_quad_p0p5_b1_b,   gradleg_quad_p0p6_b1_b,   gradleg_quad_p0p7_b1_b,   gradleg_quad_p0p8_b1_b,   gradleg_quad_p0p9_b1_b,   gradleg_quad_p0p10_b1_b,   gradleg_quad_p0p11_b1_b,   gradleg_quad_p1p2_b1_b,   gradleg_quad_p1p3_b1_b,   gradleg_quad_p1p4_b1_b,   gradleg_quad_p1p5_b1_b,   gradleg_quad_p1p6_b1_b,   gradleg_quad_p1p7_b1_b,   gradleg_quad_p1p8_b1_b,   gradleg_quad_p1p9_b1_b,   gradleg_quad_p1p10_b1_b,   gradleg_quad_p1p11_b1_b,   gradleg_quad_p2p2_b1_b,   gradleg_quad_p2p3_b1_b,   gradleg_quad_p2p4_b1_b,   gradleg_quad_p2p5_b1_b,   gradleg_quad_p2p6_b1_b,   gradleg_quad_p2p7_b1_b,   gradleg_quad_p2p8_b1_b,   gradleg_quad_p2p9_b1_b,   gradleg_quad_p2p10_b1_b,   gradleg_quad_p2p11_b1_b,   gradleg_quad_p3p2_b1_b,   gradleg_quad_p3p3_b1_b,   gradleg_quad_p3p4_b1_b,   gradleg_quad_p3p5_b1_b,   gradleg_quad_p3p6_b1_b,   gradleg_quad_p3p7_b1_b,   gradleg_quad_p3p8_b1_b,   gradleg_quad_p3p9_b1_b,   gradleg_quad_p3p10_b1_b,   gradleg_quad_p3p11_b1_b,   gradleg_quad_p4p2_b1_b,   gradleg_quad_p4p3_b1_b,   gradleg_quad_p4p4_b1_b,   gradleg_quad_p4p5_b1_b,   gradleg_quad_p4p6_b1_b,   gradleg_quad_p4p7_b1_b,   gradleg_quad_p4p8_b1_b,   gradleg_quad_p4p9_b1_b,   gradleg_quad_p4p10_b1_b,   gradleg_quad_p4p11_b1_b,   gradleg_quad_p5p2_b1_b,   gradleg_quad_p5p3_b1_b,   gradleg_quad_p5p4_b1_b,   gradleg_quad_p5p5_b1_b,   gradleg_quad_p5p6_b1_b,   gradleg_quad_p5p7_b1_b,   gradleg_quad_p5p8_b1_b,   gradleg_quad_p5p9_b1_b,   gradleg_quad_p5p10_b1_b,   gradleg_quad_p5p11_b1_b,   gradleg_quad_p6p2_b1_b,   gradleg_quad_p6p3_b1_b,   gradleg_quad_p6p4_b1_b,   gradleg_quad_p6p5_b1_b,   gradleg_quad_p6p6_b1_b,   gradleg_quad_p6p7_b1_b,   gradleg_quad_p6p8_b1_b,   gradleg_quad_p6p9_b1_b,   gradleg_quad_p6p10_b1_b,   gradleg_quad_p6p11_b1_b,   gradleg_quad_p7p2_b1_b,   gradleg_quad_p7p3_b1_b,   gradleg_quad_p7p4_b1_b,   gradleg_quad_p7p5_b1_b,   gradleg_quad_p7p6_b1_b,   gradleg_quad_p7p7_b1_b,   gradleg_quad_p7p8_b1_b,   gradleg_quad_p7p9_b1_b,   gradleg_quad_p7p10_b1_b,   gradleg_quad_p7p11_b1_b,   gradleg_quad_p8p2_b1_b,   gradleg_quad_p8p3_b1_b,   gradleg_quad_p8p4_b1_b,   gradleg_quad_p8p5_b1_b,   gradleg_quad_p8p6_b1_b,   gradleg_quad_p8p7_b1_b,   gradleg_quad_p8p8_b1_b,   gradleg_quad_p8p9_b1_b,   gradleg_quad_p8p10_b1_b,   gradleg_quad_p8p11_b1_b,   gradleg_quad_p9p2_b1_b,   gradleg_quad_p9p3_b1_b,   gradleg_quad_p9p4_b1_b,   gradleg_quad_p9p5_b1_b,   gradleg_quad_p9p6_b1_b,   gradleg_quad_p9p7_b1_b,   gradleg_quad_p9p8_b1_b,   gradleg_quad_p9p9_b1_b,   gradleg_quad_p9p10_b1_b,   gradleg_quad_p9p11_b1_b,   gradleg_quad_p10p2_b1_b,   gradleg_quad_p10p3_b1_b,   gradleg_quad_p10p4_b1_b,   gradleg_quad_p10p5_b1_b,   gradleg_quad_p10p6_b1_b,   gradleg_quad_p10p7_b1_b,   gradleg_quad_p10p8_b1_b,   gradleg_quad_p10p9_b1_b,   gradleg_quad_p10p10_b1_b,   gradleg_quad_p10p11_b1_b,   gradleg_quad_p2p0_b2_b,   gradleg_quad_p2p1_b2_b,   gradleg_quad_p2p2_b2_b,   gradleg_quad_p2p3_b2_b,   gradleg_quad_p2p4_b2_b,   gradleg_quad_p2p5_b2_b,   gradleg_quad_p2p6_b2_b,   gradleg_quad_p2p7_b2_b,   gradleg_quad_p2p8_b2_b,   gradleg_quad_p2p9_b2_b,   gradleg_quad_p2p10_b2_b,   gradleg_quad_p3p0_b2_b,   gradleg_quad_p3p1_b2_b,   gradleg_quad_p3p2_b2_b,   gradleg_quad_p3p3_b2_b,   gradleg_quad_p3p4_b2_b,   gradleg_quad_p3p5_b2_b,   gradleg_quad_p3p6_b2_b,   gradleg_quad_p3p7_b2_b,   gradleg_quad_p3p8_b2_b,   gradleg_quad_p3p9_b2_b,   gradleg_quad_p3p10_b2_b,   gradleg_quad_p4p0_b2_b,   gradleg_quad_p4p1_b2_b,   gradleg_quad_p4p2_b2_b,   gradleg_quad_p4p3_b2_b,   gradleg_quad_p4p4_b2_b,   gradleg_quad_p4p5_b2_b,   gradleg_quad_p4p6_b2_b,   gradleg_quad_p4p7_b2_b,   gradleg_quad_p4p8_b2_b,   gradleg_quad_p4p9_b2_b,   gradleg_quad_p4p10_b2_b,   gradleg_quad_p5p0_b2_b,   gradleg_quad_p5p1_b2_b,   gradleg_quad_p5p2_b2_b,   gradleg_quad_p5p3_b2_b,   gradleg_quad_p5p4_b2_b,   gradleg_quad_p5p5_b2_b,   gradleg_quad_p5p6_b2_b,   gradleg_quad_p5p7_b2_b,   gradleg_quad_p5p8_b2_b,   gradleg_quad_p5p9_b2_b,   gradleg_quad_p5p10_b2_b,   gradleg_quad_p6p0_b2_b,   gradleg_quad_p6p1_b2_b,   gradleg_quad_p6p2_b2_b,   gradleg_quad_p6p3_b2_b,   gradleg_quad_p6p4_b2_b,   gradleg_quad_p6p5_b2_b,   gradleg_quad_p6p6_b2_b,   gradleg_quad_p6p7_b2_b,   gradleg_quad_p6p8_b2_b,   gradleg_quad_p6p9_b2_b,   gradleg_quad_p6p10_b2_b,   gradleg_quad_p7p0_b2_b,   gradleg_quad_p7p1_b2_b,   gradleg_quad_p7p2_b2_b,   gradleg_quad_p7p3_b2_b,   gradleg_quad_p7p4_b2_b,   gradleg_quad_p7p5_b2_b,   gradleg_quad_p7p6_b2_b,   gradleg_quad_p7p7_b2_b,   gradleg_quad_p7p8_b2_b,   gradleg_quad_p7p9_b2_b,   gradleg_quad_p7p10_b2_b,   gradleg_quad_p8p0_b2_b,   gradleg_quad_p8p1_b2_b,   gradleg_quad_p8p2_b2_b,   gradleg_quad_p8p3_b2_b,   gradleg_quad_p8p4_b2_b,   gradleg_quad_p8p5_b2_b,   gradleg_quad_p8p6_b2_b,   gradleg_quad_p8p7_b2_b,   gradleg_quad_p8p8_b2_b,   gradleg_quad_p8p9_b2_b,   gradleg_quad_p8p10_b2_b,   gradleg_quad_p9p0_b2_b,   gradleg_quad_p9p1_b2_b,   gradleg_quad_p9p2_b2_b,   gradleg_quad_p9p3_b2_b,   gradleg_quad_p9p4_b2_b,   gradleg_quad_p9p5_b2_b,   gradleg_quad_p9p6_b2_b,   gradleg_quad_p9p7_b2_b,   gradleg_quad_p9p8_b2_b,   gradleg_quad_p9p9_b2_b,   gradleg_quad_p9p10_b2_b,   gradleg_quad_p10p0_b2_b,   gradleg_quad_p10p1_b2_b,   gradleg_quad_p10p2_b2_b,   gradleg_quad_p10p3_b2_b,   gradleg_quad_p10p4_b2_b,   gradleg_quad_p10p5_b2_b,   gradleg_quad_p10p6_b2_b,   gradleg_quad_p10p7_b2_b,   gradleg_quad_p10p8_b2_b,   gradleg_quad_p10p9_b2_b,   gradleg_quad_p10p10_b2_b,   gradleg_quad_p11p0_b2_b,   gradleg_quad_p11p1_b2_b,   gradleg_quad_p11p2_b2_b,   gradleg_quad_p11p3_b2_b,   gradleg_quad_p11p4_b2_b,   gradleg_quad_p11p5_b2_b,   gradleg_quad_p11p6_b2_b,   gradleg_quad_p11p7_b2_b,   gradleg_quad_p11p8_b2_b,   gradleg_quad_p11p9_b2_b,   gradleg_quad_p11p10_b2_b, };

    static Shapeset::shape_fn_t gradleg_quad_fn_ax[] =
    {
      gradleg_quad_p0_e1_ax_0, gradleg_quad_p0_e1_ax_1, gradleg_quad_p0_e2_ax, gradleg_quad_p0_e2_ax, gradleg_quad_p0_e3_ax_0, gradleg_quad_p0_e3_ax_1, gradleg_quad_p0_e4_ax, gradleg_quad_p0_e4_ax,
      gradleg_quad_l2_l0_ax, gradleg_quad_l2_l0_ax, gradleg_quad_l1_l2_ax, gradleg_quad_l1_l2_ax, gradleg_quad_l2_l1_ax, gradleg_quad_l2_l1_ax, gradleg_quad_l0_l2_ax, gradleg_quad_l0_l2_ax,
      gradleg_quad_l3_l0_ax_0, gradleg_quad_l3_l0_ax_1, gradleg_quad_l1_l3_ax_0, gradleg_quad_l1_l3_ax_1,  gradleg_quad_l3_l1_ax_0, gradleg_quad_l3_l1_ax_1, gradleg_quad_l0_l3_ax_0, gradleg_quad_l0_l3_ax_1,
      gradleg_quad_l4_l0_ax, gradleg_quad_l4_l0_ax, gradleg_quad_l1_l4_ax, gradleg_quad_l1_l4_ax, gradleg_quad_l4_l1_ax, gradleg_quad_l4_l1_ax, gradleg_quad_l0_l4_ax, gradleg_quad_l0_l4_ax,
      gradleg_quad_l5_l0_ax_0, gradleg_quad_l5_l0_ax_1, gradleg_quad_l1_l5_ax_0, gradleg_quad_l1_l5_ax_1,  gradleg_quad_l5_l1_ax_0, gradleg_quad_l5_l1_ax_1, gradleg_quad_l0_l5_ax_0, gradleg_quad_l0_l5_ax_1,
      gradleg_quad_l6_l0_ax, gradleg_quad_l6_l0_ax, gradleg_quad_l1_l6_ax, gradleg_quad_l1_l6_ax, gradleg_quad_l6_l1_ax, gradleg_quad_l6_l1_ax, gradleg_quad_l0_l6_ax, gradleg_quad_l0_l6_ax,
      gradleg_quad_l7_l0_ax_0, gradleg_quad_l7_l0_ax_1, gradleg_quad_l1_l7_ax_0, gradleg_quad_l1_l7_ax_1,  gradleg_quad_l7_l1_ax_0, gradleg_quad_l7_l1_ax_1, gradleg_quad_l0_l7_ax_0, gradleg_quad_l0_l7_ax_1,
      gradleg_quad_l8_l0_ax, gradleg_quad_l8_l0_ax, gradleg_quad_l1_l8_ax, gradleg_quad_l1_l8_ax, gradleg_quad_l8_l1_ax, gradleg_quad_l8_l1_ax, gradleg_quad_l0_l8_ax, gradleg_quad_l0_l8_ax,
      gradleg_quad_l9_l0_ax_0, gradleg_quad_l9_l0_ax_1, gradleg_quad_l1_l9_ax_0, gradleg_quad_l1_l9_ax_1,  gradleg_quad_l9_l1_ax_0, gradleg_quad_l9_l1_ax_1, gradleg_quad_l0_l9_ax_0, gradleg_quad_l0_l9_ax_1,
      gradleg_quad_l10_l0_ax, gradleg_quad_l10_l0_ax, gradleg_quad_l1_l10_ax, gradleg_quad_l1_l10_ax, gradleg_quad_l10_l1_ax, gradleg_quad_l10_l1_ax, gradleg_quad_l0_l10_ax, gradleg_quad_l0_l10_ax,
      gradleg_quad_l11_l0_ax_0, gradleg_quad_l11_l0_ax_1, gradleg_quad_l1_l11_ax_0, gradleg_quad_l1_l11_ax_1,  gradleg_quad_l11_l1_ax_0, gradleg_quad_l11_l1_ax_1, gradleg_quad_l0_l11_ax_0, gradleg_quad_l0_l11_ax_1,

      gradleg_quad_p0p2_b1_ax,   gradleg_quad_p0p3_b1_ax,   gradleg_quad_p0p4_b1_ax,   gradleg_quad_p0p5_b1_ax,   gradleg_quad_p0p6_b1_ax,   gradleg_quad_p0p7_b1_ax,   gradleg_quad_p0p8_b1_ax,   gradleg_quad_p0p9_b1_ax,   gradleg_quad_p0p10_b1_ax,   gradleg_quad_p0p11_b1_ax,   gradleg_quad_p1p2_b1_ax,   gradleg_quad_p1p3_b1_ax,   gradleg_quad_p1p4_b1_ax,   gradleg_quad_p1p5_b1_ax,   gradleg_quad_p1p6_b1_ax,   gradleg_quad_p1p7_b1_ax,   gradleg_quad_p1p8_b1_ax,   gradleg_quad_p1p9_b1_ax,   gradleg_quad_p1p10_b1_ax,   gradleg_quad_p1p11_b1_ax,   gradleg_quad_p2p2_b1_ax,   gradleg_quad_p2p3_b1_ax,   gradleg_quad_p2p4_b1_ax,   gradleg_quad_p2p5_b1_ax,   gradleg_quad_p2p6_b1_ax,   gradleg_quad_p2p7_b1_ax,   gradleg_quad_p2p8_b1_ax,   gradleg_quad_p2p9_b1_ax,   gradleg_quad_p2p10_b1_ax,   gradleg_quad_p2p11_b1_ax,   gradleg_quad_p3p2_b1_ax,   gradleg_quad_p3p3_b1_ax,   gradleg_quad_p3p4_b1_ax,   gradleg_quad_p3p5_b1_ax,   gradleg_quad_p3p6_b1_ax,   gradleg_quad_p3p7_b1_ax,   gradleg_quad_p3p8_b1_ax,   gradleg_quad_p3p9_b1_ax,   gradleg_quad_p3p10_b1_ax,   gradleg_quad_p3p11_b1_ax,   gradleg_quad_p4p2_b1_ax,   gradleg_quad_p4p3_b1_ax,   gradleg_quad_p4p4_b1_ax,   gradleg_quad_p4p5_b1_ax,   gradleg_quad_p4p6_b1_ax,   gradleg_quad_p4p7_b1_ax,   gradleg_quad_p4p8_b1_ax,   gradleg_quad_p4p9_b1_ax,   gradleg_quad_p4p10_b1_ax,   gradleg_quad_p4p11_b1_ax,   gradleg_quad_p5p2_b1_ax,   gradleg_quad_p5p3_b1_ax,   gradleg_quad_p5p4_b1_ax,   gradleg_quad_p5p5_b1_ax,   gradleg_quad_p5p6_b1_ax,   gradleg_quad_p5p7_b1_ax,   gradleg_quad_p5p8_b1_ax,   gradleg_quad_p5p9_b1_ax,   gradleg_quad_p5p10_b1_ax,   gradleg_quad_p5p11_b1_ax,   gradleg_quad_p6p2_b1_ax,   gradleg_quad_p6p3_b1_ax,   gradleg_quad_p6p4_b1_ax,   gradleg_quad_p6p5_b1_ax,   gradleg_quad_p6p6_b1_ax,   gradleg_quad_p6p7_b1_ax,   gradleg_quad_p6p8_b1_ax,   gradleg_quad_p6p9_b1_ax,   gradleg_quad_p6p10_b1_ax,   gradleg_quad_p6p11_b1_ax,   gradleg_quad_p7p2_b1_ax,   gradleg_quad_p7p3_b1_ax,   gradleg_quad_p7p4_b1_ax,   gradleg_quad_p7p5_b1_ax,   gradleg_quad_p7p6_b1_ax,   gradleg_quad_p7p7_b1_ax,   gradleg_quad_p7p8_b1_ax,   gradleg_quad_p7p9_b1_ax,   gradleg_quad_p7p10_b1_ax,   gradleg_quad_p7p11_b1_ax,   gradleg_quad_p8p2_b1_ax,   gradleg_quad_p8p3_b1_ax,   gradleg_quad_p8p4_b1_ax,   gradleg_quad_p8p5_b1_ax,   gradleg_quad_p8p6_b1_ax,   gradleg_quad_p8p7_b1_ax,   gradleg_quad_p8p8_b1_ax,   gradleg_quad_p8p9_b1_ax,   gradleg_quad_p8p10_b1_ax,   gradleg_quad_p8p11_b1_ax,   gradleg_quad_p9p2_b1_ax,   gradleg_quad_p9p3_b1_ax,   gradleg_quad_p9p4_b1_ax,   gradleg_quad_p9p5_b1_ax,   gradleg_quad_p9p6_b1_ax,   gradleg_quad_p9p7_b1_ax,   gradleg_quad_p9p8_b1_ax,   gradleg_quad_p9p9_b1_ax,   gradleg_quad_p9p10_b1_ax,   gradleg_quad_p9p11_b1_ax,   gradleg_quad_p10p2_b1_ax,   gradleg_quad_p10p3_b1_ax,   gradleg_quad_p10p4_b1_ax,   gradleg_quad_p10p5_b1_ax,   gradleg_quad_p10p6_b1_ax,   gradleg_quad_p10p7_b1_ax,   gradleg_quad_p10p8_b1_ax,   gradleg_quad_p10p9_b1_ax,   gradleg_quad_p10p10_b1_ax,   gradleg_quad_p10p11_b1_ax,   gradleg_quad_p2p0_b2_ax,   gradleg_quad_p2p1_b2_ax,   gradleg_quad_p2p2_b2_ax,   gradleg_quad_p2p3_b2_ax,   gradleg_quad_p2p4_b2_ax,   gradleg_quad_p2p5_b2_ax,   gradleg_quad_p2p6_b2_ax,   gradleg_quad_p2p7_b2_ax,   gradleg_quad_p2p8_b2_ax,   gradleg_quad_p2p9_b2_ax,   gradleg_quad_p2p10_b2_ax,   gradleg_quad_p3p0_b2_ax,   gradleg_quad_p3p1_b2_ax,   gradleg_quad_p3p2_b2_ax,   gradleg_quad_p3p3_b2_ax,   gradleg_quad_p3p4_b2_ax,   gradleg_quad_p3p5_b2_ax,   gradleg_quad_p3p6_b2_ax,   gradleg_quad_p3p7_b2_ax,   gradleg_quad_p3p8_b2_ax,   gradleg_quad_p3p9_b2_ax,   gradleg_quad_p3p10_b2_ax,   gradleg_quad_p4p0_b2_ax,   gradleg_quad_p4p1_b2_ax,   gradleg_quad_p4p2_b2_ax,   gradleg_quad_p4p3_b2_ax,   gradleg_quad_p4p4_b2_ax,   gradleg_quad_p4p5_b2_ax,   gradleg_quad_p4p6_b2_ax,   gradleg_quad_p4p7_b2_ax,   gradleg_quad_p4p8_b2_ax,   gradleg_quad_p4p9_b2_ax,   gradleg_quad_p4p10_b2_ax,   gradleg_quad_p5p0_b2_ax,   gradleg_quad_p5p1_b2_ax,   gradleg_quad_p5p2_b2_ax,   gradleg_quad_p5p3_b2_ax,   gradleg_quad_p5p4_b2_ax,   gradleg_quad_p5p5_b2_ax,   gradleg_quad_p5p6_b2_ax,   gradleg_quad_p5p7_b2_ax,   gradleg_quad_p5p8_b2_ax,   gradleg_quad_p5p9_b2_ax,   gradleg_quad_p5p10_b2_ax,   gradleg_quad_p6p0_b2_ax,   gradleg_quad_p6p1_b2_ax,   gradleg_quad_p6p2_b2_ax,   gradleg_quad_p6p3_b2_ax,   gradleg_quad_p6p4_b2_ax,   gradleg_quad_p6p5_b2_ax,   gradleg_quad_p6p6_b2_ax,   gradleg_quad_p6p7_b2_ax,   gradleg_quad_p6p8_b2_ax,   gradleg_quad_p6p9_b2_ax,   gradleg_quad_p6p10_b2_ax,   gradleg_quad_p7p0_b2_ax,   gradleg_quad_p7p1_b2_ax,   gradleg_quad_p7p2_b2_ax,   gradleg_quad_p7p3_b2_ax,   gradleg_quad_p7p4_b2_ax,   gradleg_quad_p7p5_b2_ax,   gradleg_quad_p7p6_b2_ax,   gradleg_quad_p7p7_b2_ax,   gradleg_quad_p7p8_b2_ax,   gradleg_quad_p7p9_b2_ax,   gradleg_quad_p7p10_b2_ax,   gradleg_quad_p8p0_b2_ax,   gradleg_quad_p8p1_b2_ax,   gradleg_quad_p8p2_b2_ax,   gradleg_quad_p8p3_b2_ax,   gradleg_quad_p8p4_b2_ax,   gradleg_quad_p8p5_b2_ax,   gradleg_quad_p8p6_b2_ax,   gradleg_quad_p8p7_b2_ax,   gradleg_quad_p8p8_b2_ax,   gradleg_quad_p8p9_b2_ax,   gradleg_quad_p8p10_b2_ax,   gradleg_quad_p9p0_b2_ax,   gradleg_quad_p9p1_b2_ax,   gradleg_quad_p9p2_b2_ax,   gradleg_quad_p9p3_b2_ax,   gradleg_quad_p9p4_b2_ax,   gradleg_quad_p9p5_b2_ax,   gradleg_quad_p9p6_b2_ax,   gradleg_quad_p9p7_b2_ax,   gradleg_quad_p9p8_b2_ax,   gradleg_quad_p9p9_b2_ax,   gradleg_quad_p9p10_b2_ax,   gradleg_quad_p10p0_b2_ax,   gradleg_quad_p10p1_b2_ax,   gradleg_quad_p10p2_b2_ax,   gradleg_quad_p10p3_b2_ax,   gradleg_quad_p10p4_b2_ax,   gradleg_quad_p10p5_b2_ax,   gradleg_quad_p10p6_b2_ax,   gradleg_quad_p10p7_b2_ax,   gradleg_quad_p10p8_b2_ax,   gradleg_quad_p10p9_b2_ax,   gradleg_quad_p10p10_b2_ax,   gradleg_quad_p11p0_b2_ax,   gradleg_quad_p11p1_b2_ax,   gradleg_quad_p11p2_b2_ax,   gradleg_quad_p11p3_b2_ax,   gradleg_quad_p11p4_b2_ax,   gradleg_quad_p11p5_b2_ax,   gradleg_quad_p11p6_b2_ax,   gradleg_quad_p11p7_b2_ax,   gradleg_quad_p11p8_b2_ax,   gradleg_quad_p11p9_b2_ax,   gradleg_quad_p11p10_b2_ax, };

    static Shapeset::shape_fn_t gradleg_quad_fn_bx[] =
    {
      gradleg_quad_p0_e1_bx, gradleg_quad_p0_e1_bx, gradleg_quad_p0_e2_bx_0, gradleg_quad_p0_e2_bx_1,  gradleg_quad_p0_e3_bx, gradleg_quad_p0_e3_bx, gradleg_quad_p0_e4_bx_0, gradleg_quad_p0_e4_bx_1,
      gradleg_quad_l2_l0_bx, gradleg_quad_l2_l0_bx, gradleg_quad_l1_l2_bx, gradleg_quad_l1_l2_bx, gradleg_quad_l2_l1_bx, gradleg_quad_l2_l1_bx, gradleg_quad_l0_l2_bx,  gradleg_quad_l0_l2_bx,
      gradleg_quad_l3_l0_bx_0, gradleg_quad_l3_l0_bx_1, gradleg_quad_l1_l3_bx_0, gradleg_quad_l1_l3_bx_1, gradleg_quad_l3_l1_bx_0, gradleg_quad_l3_l1_bx_1, gradleg_quad_l0_l3_bx_0, gradleg_quad_l0_l3_bx_1,
      gradleg_quad_l4_l0_bx, gradleg_quad_l4_l0_bx, gradleg_quad_l1_l4_bx, gradleg_quad_l1_l4_bx, gradleg_quad_l4_l1_bx, gradleg_quad_l4_l1_bx, gradleg_quad_l0_l4_bx,  gradleg_quad_l0_l4_bx,
      gradleg_quad_l5_l0_bx_0, gradleg_quad_l5_l0_bx_1, gradleg_quad_l1_l5_bx_0, gradleg_quad_l1_l5_bx_1, gradleg_quad_l5_l1_bx_0, gradleg_quad_l5_l1_bx_1, gradleg_quad_l0_l5_bx_0, gradleg_quad_l0_l5_bx_1,
      gradleg_quad_l6_l0_bx, gradleg_quad_l6_l0_bx, gradleg_quad_l1_l6_bx, gradleg_quad_l1_l6_bx, gradleg_quad_l6_l1_bx, gradleg_quad_l6_l1_bx, gradleg_quad_l0_l6_bx,  gradleg_quad_l0_l6_bx,
      gradleg_quad_l7_l0_bx_0, gradleg_quad_l7_l0_bx_1, gradleg_quad_l1_l7_bx_0, gradleg_quad_l1_l7_bx_1, gradleg_quad_l7_l1_bx_0, gradleg_quad_l7_l1_bx_1, gradleg_quad_l0_l7_bx_0, gradleg_quad_l0_l7_bx_1,
      gradleg_quad_l8_l0_bx, gradleg_quad_l8_l0_bx, gradleg_quad_l1_l8_bx, gradleg_quad_l1_l8_bx, gradleg_quad_l8_l1_bx, gradleg_quad_l8_l1_bx, gradleg_quad_l0_l8_bx,  gradleg_quad_l0_l8_bx,
      gradleg_quad_l9_l0_bx_0, gradleg_quad_l9_l0_bx_1, gradleg_quad_l1_l9_bx_0, gradleg_quad_l1_l9_bx_1, gradleg_quad_l9_l1_bx_0, gradleg_quad_l9_l1_bx_1, gradleg_quad_l0_l9_bx_0, gradleg_quad_l0_l9_bx_1,
      gradleg_quad_l10_l0_bx, gradleg_quad_l10_l0_bx, gradleg_quad_l1_l10_bx, gradleg_quad_l1_l10_bx, gradleg_quad_l10_l1_bx, gradleg_quad_l10_l1_bx, gradleg_quad_l0_l10_bx,  gradleg_quad_l0_l10_bx,
      gradleg_quad_l11_l0_bx_0, gradleg_quad_l11_l0_bx_1, gradleg_quad_l1_l11_bx_0, gradleg_quad_l1_l11_bx_1, gradleg_quad_l11_l1_bx_0, gradleg_quad_l11_l1_bx_1, gradleg_quad_l0_l11_bx_0, gradleg_quad_l0_l11_bx_1,

      gradleg_quad_p0p2_b1_bx,   gradleg_quad_p0p3_b1_bx,   gradleg_quad_p0p4_b1_bx,   gradleg_quad_p0p5_b1_bx,   gradleg_quad_p0p6_b1_bx,   gradleg_quad_p0p7_b1_bx,   gradleg_quad_p0p8_b1_bx,   gradleg_quad_p0p9_b1_bx,   gradleg_quad_p0p10_b1_bx,   gradleg_quad_p0p11_b1_bx,   gradleg_quad_p1p2_b1_bx,   gradleg_quad_p1p3_b1_bx,   gradleg_quad_p1p4_b1_bx,   gradleg_quad_p1p5_b1_bx,   gradleg_quad_p1p6_b1_bx,   gradleg_quad_p1p7_b1_bx,   gradleg_quad_p1p8_b1_bx,   gradleg_quad_p1p9_b1_bx,   gradleg_quad_p1p10_b1_bx,   gradleg_quad_p1p11_b1_bx,   gradleg_quad_p2p2_b1_bx,   gradleg_quad_p2p3_b1_bx,   gradleg_quad_p2p4_b1_bx,   gradleg_quad_p2p5_b1_bx,   gradleg_quad_p2p6_b1_bx,   gradleg_quad_p2p7_b1_bx,   gradleg_quad_p2p8_b1_bx,   gradleg_quad_p2p9_b1_bx,   gradleg_quad_p2p10_b1_bx,   gradleg_quad_p2p11_b1_bx,   gradleg_quad_p3p2_b1_bx,   gradleg_quad_p3p3_b1_bx,   gradleg_quad_p3p4_b1_bx,   gradleg_quad_p3p5_b1_bx,   gradleg_quad_p3p6_b1_bx,   gradleg_quad_p3p7_b1_bx,   gradleg_quad_p3p8_b1_bx,   gradleg_quad_p3p9_b1_bx,   gradleg_quad_p3p10_b1_bx,   gradleg_quad_p3p11_b1_bx,   gradleg_quad_p4p2_b1_bx,   gradleg_quad_p4p3_b1_bx,   gradleg_quad_p4p4_b1_bx,   gradleg_quad_p4p5_b1_bx,   gradleg_quad_p4p6_b1_bx,   gradleg_quad_p4p7_b1_bx,   gradleg_quad_p4p8_b1_bx,   gradleg_quad_p4p9_b1_bx,   gradleg_quad_p4p10_b1_bx,   gradleg_quad_p4p11_b1_bx,   gradleg_quad_p5p2_b1_bx,   gradleg_quad_p5p3_b1_bx,   gradleg_quad_p5p4_b1_bx,   gradleg_quad_p5p5_b1_bx,   gradleg_quad_p5p6_b1_bx,   gradleg_quad_p5p7_b1_bx,   gradleg_quad_p5p8_b1_bx,   gradleg_quad_p5p9_b1_bx,   gradleg_quad_p5p10_b1_bx,   gradleg_quad_p5p11_b1_bx,   gradleg_quad_p6p2_b1_bx,   gradleg_quad_p6p3_b1_bx,   gradleg_quad_p6p4_b1_bx,   gradleg_quad_p6p5_b1_bx,   gradleg_quad_p6p6_b1_bx,   gradleg_quad_p6p7_b1_bx,   gradleg_quad_p6p8_b1_bx,   gradleg_quad_p6p9_b1_bx,   gradleg_quad_p6p10_b1_bx,   gradleg_quad_p6p11_b1_bx,   gradleg_quad_p7p2_b1_bx,   gradleg_quad_p7p3_b1_bx,   gradleg_quad_p7p4_b1_bx,   gradleg_quad_p7p5_b1_bx,   gradleg_quad_p7p6_b1_bx,   gradleg_quad_p7p7_b1_bx,   gradleg_quad_p7p8_b1_bx,   gradleg_quad_p7p9_b1_bx,   gradleg_quad_p7p10_b1_bx,   gradleg_quad_p7p11_b1_bx,   gradleg_quad_p8p2_b1_bx,   gradleg_quad_p8p3_b1_bx,   gradleg_quad_p8p4_b1_bx,   gradleg_quad_p8p5_b1_bx,   gradleg_quad_p8p6_b1_bx,   gradleg_quad_p8p7_b1_bx,   gradleg_quad_p8p8_b1_bx,   gradleg_quad_p8p9_b1_bx,   gradleg_quad_p8p10_b1_bx,   gradleg_quad_p8p11_b1_bx,   gradleg_quad_p9p2_b1_bx,   gradleg_quad_p9p3_b1_bx,   gradleg_quad_p9p4_b1_bx,   gradleg_quad_p9p5_b1_bx,   gradleg_quad_p9p6_b1_bx,   gradleg_quad_p9p7_b1_bx,   gradleg_quad_p9p8_b1_bx,   gradleg_quad_p9p9_b1_bx,   gradleg_quad_p9p10_b1_bx,   gradleg_quad_p9p11_b1_bx,   gradleg_quad_p10p2_b1_bx,   gradleg_quad_p10p3_b1_bx,   gradleg_quad_p10p4_b1_bx,   gradleg_quad_p10p5_b1_bx,   gradleg_quad_p10p6_b1_bx,   gradleg_quad_p10p7_b1_bx,   gradleg_quad_p10p8_b1_bx,   gradleg_quad_p10p9_b1_bx,   gradleg_quad_p10p10_b1_bx,   gradleg_quad_p10p11_b1_bx,   gradleg_quad_p2p0_b2_bx,   gradleg_quad_p2p1_b2_bx,   gradleg_quad_p2p2_b2_bx,   gradleg_quad_p2p3_b2_bx,   gradleg_quad_p2p4_b2_bx,   gradleg_quad_p2p5_b2_bx,   gradleg_quad_p2p6_b2_bx,   gradleg_quad_p2p7_b2_bx,   gradleg_quad_p2p8_b2_bx,   gradleg_quad_p2p9_b2_bx,   gradleg_quad_p2p10_b2_bx,   gradleg_quad_p3p0_b2_bx,   gradleg_quad_p3p1_b2_bx,   gradleg_quad_p3p2_b2_bx,   gradleg_quad_p3p3_b2_bx,   gradleg_quad_p3p4_b2_bx,   gradleg_quad_p3p5_b2_bx,   gradleg_quad_p3p6_b2_bx,   gradleg_quad_p3p7_b2_bx,   gradleg_quad_p3p8_b2_bx,   gradleg_quad_p3p9_b2_bx,   gradleg_quad_p3p10_b2_bx,   gradleg_quad_p4p0_b2_bx,   gradleg_quad_p4p1_b2_bx,   gradleg_quad_p4p2_b2_bx,   gradleg_quad_p4p3_b2_bx,   gradleg_quad_p4p4_b2_bx,   gradleg_quad_p4p5_b2_bx,   gradleg_quad_p4p6_b2_bx,   gradleg_quad_p4p7_b2_bx,   gradleg_quad_p4p8_b2_bx,   gradleg_quad_p4p9_b2_bx,   gradleg_quad_p4p10_b2_bx,   gradleg_quad_p5p0_b2_bx,   gradleg_quad_p5p1_b2_bx,   gradleg_quad_p5p2_b2_bx,   gradleg_quad_p5p3_b2_bx,   gradleg_quad_p5p4_b2_bx,   gradleg_quad_p5p5_b2_bx,   gradleg_quad_p5p6_b2_bx,   gradleg_quad_p5p7_b2_bx,   gradleg_quad_p5p8_b2_bx,   gradleg_quad_p5p9_b2_bx,   gradleg_quad_p5p10_b2_bx,   gradleg_quad_p6p0_b2_bx,   gradleg_quad_p6p1_b2_bx,   gradleg_quad_p6p2_b2_bx,   gradleg_quad_p6p3_b2_bx,   gradleg_quad_p6p4_b2_bx,   gradleg_quad_p6p5_b2_bx,   gradleg_quad_p6p6_b2_bx,   gradleg_quad_p6p7_b2_bx,   gradleg_quad_p6p8_b2_bx,   gradleg_quad_p6p9_b2_bx,   gradleg_quad_p6p10_b2_bx,   gradleg_quad_p7p0_b2_bx,   gradleg_quad_p7p1_b2_bx,   gradleg_quad_p7p2_b2_bx,   gradleg_quad_p7p3_b2_bx,   gradleg_quad_p7p4_b2_bx,   gradleg_quad_p7p5_b2_bx,   gradleg_quad_p7p6_b2_bx,   gradleg_quad_p7p7_b2_bx,   gradleg_quad_p7p8_b2_bx,   gradleg_quad_p7p9_b2_bx,   gradleg_quad_p7p10_b2_bx,   gradleg_quad_p8p0_b2_bx,   gradleg_quad_p8p1_b2_bx,   gradleg_quad_p8p2_b2_bx,   gradleg_quad_p8p3_b2_bx,   gradleg_quad_p8p4_b2_bx,   gradleg_quad_p8p5_b2_bx,   gradleg_quad_p8p6_b2_bx,   gradleg_quad_p8p7_b2_bx,   gradleg_quad_p8p8_b2_bx,   gradleg_quad_p8p9_b2_bx,   gradleg_quad_p8p10_b2_bx,   gradleg_quad_p9p0_b2_bx,   gradleg_quad_p9p1_b2_bx,   gradleg_quad_p9p2_b2_bx,   gradleg_quad_p9p3_b2_bx,   gradleg_quad_p9p4_b2_bx,   gradleg_quad_p9p5_b2_bx,   gradleg_quad_p9p6_b2_bx,   gradleg_quad_p9p7_b2_bx,   gradleg_quad_p9p8_b2_bx,   gradleg_quad_p9p9_b2_bx,   gradleg_quad_p9p10_b2_bx,   gradleg_quad_p10p0_b2_bx,   gradleg_quad_p10p1_b2_bx,   gradleg_quad_p10p2_b2_bx,   gradleg_quad_p10p3_b2_bx,   gradleg_quad_p10p4_b2_bx,   gradleg_quad_p10p5_b2_bx,   gradleg_quad_p10p6_b2_bx,   gradleg_quad_p10p7_b2_bx,   gradleg_quad_p10p8_b2_bx,   gradleg_quad_p10p9_b2_bx,   gradleg_quad_p10p10_b2_bx,   gradleg_quad_p11p0_b2_bx,   gradleg_quad_p11p1_b2_bx,   gradleg_quad_p11p2_b2_bx,   gradleg_quad_p11p3_b2_bx,   gradleg_quad_p11p4_b2_bx,   gradleg_quad_p11p5_b2_bx,   gradleg_quad_p11p6_b2_bx,   gradleg_quad_p11p7_b2_bx,   gradleg_quad_p11p8_b2_bx,   gradleg_quad_p11p9_b2_bx,   gradleg_quad_p11p10_b2_bx, };

    static Shapeset::shape_fn_t gradleg_quad_fn_ay[] =
    {
      gradleg_quad_p0_e1_ay_0, gradleg_quad_p0_e1_ay_1, gradleg_quad_p0_e2_ay, gradleg_quad_p0_e2_ay, gradleg_quad_p0_e3_ay_0, gradleg_quad_p0_e3_ay_1, gradleg_quad_p0_e4_ay, gradleg_quad_p0_e4_ay,
      gradleg_quad_l2_l0_ay, gradleg_quad_l2_l0_ay, gradleg_quad_l1_l2_ay, gradleg_quad_l1_l2_ay, gradleg_quad_l2_l1_ay, gradleg_quad_l2_l1_ay, gradleg_quad_l0_l2_ay, gradleg_quad_l0_l2_ay,
      gradleg_quad_l3_l0_ay_0, gradleg_quad_l3_l0_ay_1, gradleg_quad_l1_l3_ay_0, gradleg_quad_l1_l3_ay_1,  gradleg_quad_l3_l1_ay_0, gradleg_quad_l3_l1_ay_1, gradleg_quad_l0_l3_ay_0, gradleg_quad_l0_l3_ay_1,
      gradleg_quad_l4_l0_ay, gradleg_quad_l4_l0_ay, gradleg_quad_l1_l4_ay, gradleg_quad_l1_l4_ay, gradleg_quad_l4_l1_ay, gradleg_quad_l4_l1_ay, gradleg_quad_l0_l4_ay, gradleg_quad_l0_l4_ay,
      gradleg_quad_l5_l0_ay_0, gradleg_quad_l5_l0_ay_1, gradleg_quad_l1_l5_ay_0, gradleg_quad_l1_l5_ay_1,  gradleg_quad_l5_l1_ay_0, gradleg_quad_l5_l1_ay_1, gradleg_quad_l0_l5_ay_0, gradleg_quad_l0_l5_ay_1,
      gradleg_quad_l6_l0_ay, gradleg_quad_l6_l0_ay, gradleg_quad_l1_l6_ay, gradleg_quad_l1_l6_ay, gradleg_quad_l6_l1_ay, gradleg_quad_l6_l1_ay, gradleg_quad_l0_l6_ay, gradleg_quad_l0_l6_ay,
      gradleg_quad_l7_l0_ay_0, gradleg_quad_l7_l0_ay_1, gradleg_quad_l1_l7_ay_0, gradleg_quad_l1_l7_ay_1,  gradleg_quad_l7_l1_ay_0, gradleg_quad_l7_l1_ay_1, gradleg_quad_l0_l7_ay_0, gradleg_quad_l0_l7_ay_1,
      gradleg_quad_l8_l0_ay, gradleg_quad_l8_l0_ay, gradleg_quad_l1_l8_ay, gradleg_quad_l1_l8_ay, gradleg_quad_l8_l1_ay, gradleg_quad_l8_l1_ay, gradleg_quad_l0_l8_ay, gradleg_quad_l0_l8_ay,
      gradleg_quad_l9_l0_ay_0, gradleg_quad_l9_l0_ay_1, gradleg_quad_l1_l9_ay_0, gradleg_quad_l1_l9_ay_1,  gradleg_quad_l9_l1_ay_0, gradleg_quad_l9_l1_ay_1, gradleg_quad_l0_l9_ay_0, gradleg_quad_l0_l9_ay_1,
      gradleg_quad_l10_l0_ay, gradleg_quad_l10_l0_ay, gradleg_quad_l1_l10_ay, gradleg_quad_l1_l10_ay, gradleg_quad_l10_l1_ay, gradleg_quad_l10_l1_ay, gradleg_quad_l0_l10_ay, gradleg_quad_l0_l10_ay,
      gradleg_quad_l11_l0_ay_0, gradleg_quad_l11_l0_ay_1, gradleg_quad_l1_l11_ay_0, gradleg_quad_l1_l11_ay_1,  gradleg_quad_l11_l1_ay_0, gradleg_quad_l11_l1_ay_1, gradleg_quad_l0_l11_ay_0, gradleg_quad_l0_l11_ay_1,

      gradleg_quad_p0p2_b1_ay,   gradleg_quad_p0p3_b1_ay,   gradleg_quad_p0p4_b1_ay,   gradleg_quad_p0p5_b1_ay,   gradleg_quad_p0p6_b1_ay,   gradleg_quad_p0p7_b1_ay,   gradleg_quad_p0p8_b1_ay,   gradleg_quad_p0p9_b1_ay,   gradleg_quad_p0p10_b1_ay,   gradleg_quad_p0p11_b1_ay,   gradleg_quad_p1p2_b1_ay,   gradleg_quad_p1p3_b1_ay,   gradleg_quad_p1p4_b1_ay,   gradleg_quad_p1p5_b1_ay,   gradleg_quad_p1p6_b1_ay,   gradleg_quad_p1p7_b1_ay,   gradleg_quad_p1p8_b1_ay,   gradleg_quad_p1p9_b1_ay,   gradleg_quad_p1p10_b1_ay,   gradleg_quad_p1p11_b1_ay,   gradleg_quad_p2p2_b1_ay,   gradleg_quad_p2p3_b1_ay,   gradleg_quad_p2p4_b1_ay,   gradleg_quad_p2p5_b1_ay,   gradleg_quad_p2p6_b1_ay,   gradleg_quad_p2p7_b1_ay,   gradleg_quad_p2p8_b1_ay,   gradleg_quad_p2p9_b1_ay,   gradleg_quad_p2p10_b1_ay,   gradleg_quad_p2p11_b1_ay,   gradleg_quad_p3p2_b1_ay,   gradleg_quad_p3p3_b1_ay,   gradleg_quad_p3p4_b1_ay,   gradleg_quad_p3p5_b1_ay,   gradleg_quad_p3p6_b1_ay,   gradleg_quad_p3p7_b1_ay,   gradleg_quad_p3p8_b1_ay,   gradleg_quad_p3p9_b1_ay,   gradleg_quad_p3p10_b1_ay,   gradleg_quad_p3p11_b1_ay,   gradleg_quad_p4p2_b1_ay,   gradleg_quad_p4p3_b1_ay,   gradleg_quad_p4p4_b1_ay,   gradleg_quad_p4p5_b1_ay,   gradleg_quad_p4p6_b1_ay,   gradleg_quad_p4p7_b1_ay,   gradleg_quad_p4p8_b1_ay,   gradleg_quad_p4p9_b1_ay,   gradleg_quad_p4p10_b1_ay,   gradleg_quad_p4p11_b1_ay,   gradleg_quad_p5p2_b1_ay,   gradleg_quad_p5p3_b1_ay,   gradleg_quad_p5p4_b1_ay,   gradleg_quad_p5p5_b1_ay,   gradleg_quad_p5p6_b1_ay,   gradleg_quad_p5p7_b1_ay,   gradleg_quad_p5p8_b1_ay,   gradleg_quad_p5p9_b1_ay,   gradleg_quad_p5p10_b1_ay,   gradleg_quad_p5p11_b1_ay,   gradleg_quad_p6p2_b1_ay,   gradleg_quad_p6p3_b1_ay,   gradleg_quad_p6p4_b1_ay,   gradleg_quad_p6p5_b1_ay,   gradleg_quad_p6p6_b1_ay,   gradleg_quad_p6p7_b1_ay,   gradleg_quad_p6p8_b1_ay,   gradleg_quad_p6p9_b1_ay,   gradleg_quad_p6p10_b1_ay,   gradleg_quad_p6p11_b1_ay,   gradleg_quad_p7p2_b1_ay,   gradleg_quad_p7p3_b1_ay,   gradleg_quad_p7p4_b1_ay,   gradleg_quad_p7p5_b1_ay,   gradleg_quad_p7p6_b1_ay,   gradleg_quad_p7p7_b1_ay,   gradleg_quad_p7p8_b1_ay,   gradleg_quad_p7p9_b1_ay,   gradleg_quad_p7p10_b1_ay,   gradleg_quad_p7p11_b1_ay,   gradleg_quad_p8p2_b1_ay,   gradleg_quad_p8p3_b1_ay,   gradleg_quad_p8p4_b1_ay,   gradleg_quad_p8p5_b1_ay,   gradleg_quad_p8p6_b1_ay,   gradleg_quad_p8p7_b1_ay,   gradleg_quad_p8p8_b1_ay,   gradleg_quad_p8p9_b1_ay,   gradleg_quad_p8p10_b1_ay,   gradleg_quad_p8p11_b1_ay,   gradleg_quad_p9p2_b1_ay,   gradleg_quad_p9p3_b1_ay,   gradleg_quad_p9p4_b1_ay,   gradleg_quad_p9p5_b1_ay,   gradleg_quad_p9p6_b1_ay,   gradleg_quad_p9p7_b1_ay,   gradleg_quad_p9p8_b1_ay,   gradleg_quad_p9p9_b1_ay,   gradleg_quad_p9p10_b1_ay,   gradleg_quad_p9p11_b1_ay,   gradleg_quad_p10p2_b1_ay,   gradleg_quad_p10p3_b1_ay,   gradleg_quad_p10p4_b1_ay,   gradleg_quad_p10p5_b1_ay,   gradleg_quad_p10p6_b1_ay,   gradleg_quad_p10p7_b1_ay,   gradleg_quad_p10p8_b1_ay,   gradleg_quad_p10p9_b1_ay,   gradleg_quad_p10p10_b1_ay,   gradleg_quad_p10p11_b1_ay,   gradleg_quad_p2p0_b2_ay,   gradleg_quad_p2p1_b2_ay,   gradleg_quad_p2p2_b2_ay,   gradleg_quad_p2p3_b2_ay,   gradleg_quad_p2p4_b2_ay,   gradleg_quad_p2p5_b2_ay,   gradleg_quad_p2p6_b2_ay,   gradleg_quad_p2p7_b2_ay,   gradleg_quad_p2p8_b2_ay,   gradleg_quad_p2p9_b2_ay,   gradleg_quad_p2p10_b2_ay,   gradleg_quad_p3p0_b2_ay,   gradleg_quad_p3p1_b2_ay,   gradleg_quad_p3p2_b2_ay,   gradleg_quad_p3p3_b2_ay,   gradleg_quad_p3p4_b2_ay,   gradleg_quad_p3p5_b2_ay,   gradleg_quad_p3p6_b2_ay,   gradleg_quad_p3p7_b2_ay,   gradleg_quad_p3p8_b2_ay,   gradleg_quad_p3p9_b2_ay,   gradleg_quad_p3p10_b2_ay,   gradleg_quad_p4p0_b2_ay,   gradleg_quad_p4p1_b2_ay,   gradleg_quad_p4p2_b2_ay,   gradleg_quad_p4p3_b2_ay,   gradleg_quad_p4p4_b2_ay,   gradleg_quad_p4p5_b2_ay,   gradleg_quad_p4p6_b2_ay,   gradleg_quad_p4p7_b2_ay,   gradleg_quad_p4p8_b2_ay,   gradleg_quad_p4p9_b2_ay,   gradleg_quad_p4p10_b2_ay,   gradleg_quad_p5p0_b2_ay,   gradleg_quad_p5p1_b2_ay,   gradleg_quad_p5p2_b2_ay,   gradleg_quad_p5p3_b2_ay,   gradleg_quad_p5p4_b2_ay,   gradleg_quad_p5p5_b2_ay,   gradleg_quad_p5p6_b2_ay,   gradleg_quad_p5p7_b2_ay,   gradleg_quad_p5p8_b2_ay,   gradleg_quad_p5p9_b2_ay,   gradleg_quad_p5p10_b2_ay,   gradleg_quad_p6p0_b2_ay,   gradleg_quad_p6p1_b2_ay,   gradleg_quad_p6p2_b2_ay,   gradleg_quad_p6p3_b2_ay,   gradleg_quad_p6p4_b2_ay,   gradleg_quad_p6p5_b2_ay,   gradleg_quad_p6p6_b2_ay,   gradleg_quad_p6p7_b2_ay,   gradleg_quad_p6p8_b2_ay,   gradleg_quad_p6p9_b2_ay,   gradleg_quad_p6p10_b2_ay,   gradleg_quad_p7p0_b2_ay,   gradleg_quad_p7p1_b2_ay,   gradleg_quad_p7p2_b2_ay,   gradleg_quad_p7p3_b2_ay,   gradleg_quad_p7p4_b2_ay,   gradleg_quad_p7p5_b2_ay,   gradleg_quad_p7p6_b2_ay,   gradleg_quad_p7p7_b2_ay,   gradleg_quad_p7p8_b2_ay,   gradleg_quad_p7p9_b2_ay,   gradleg_quad_p7p10_b2_ay,   gradleg_quad_p8p0_b2_ay,   gradleg_quad_p8p1_b2_ay,   gradleg_quad_p8p2_b2_ay,   gradleg_quad_p8p3_b2_ay,   gradleg_quad_p8p4_b2_ay,   gradleg_quad_p8p5_b2_ay,   gradleg_quad_p8p6_b2_ay,   gradleg_quad_p8p7_b2_ay,   gradleg_quad_p8p8_b2_ay,   gradleg_quad_p8p9_b2_ay,   gradleg_quad_p8p10_b2_ay,   gradleg_quad_p9p0_b2_ay,   gradleg_quad_p9p1_b2_ay,   gradleg_quad_p9p2_b2_ay,   gradleg_quad_p9p3_b2_ay,   gradleg_quad_p9p4_b2_ay,   gradleg_quad_p9p5_b2_ay,   gradleg_quad_p9p6_b2_ay,   gradleg_quad_p9p7_b2_ay,   gradleg_quad_p9p8_b2_ay,   gradleg_quad_p9p9_b2_ay,   gradleg_quad_p9p10_b2_ay,   gradleg_quad_p10p0_b2_ay,   gradleg_quad_p10p1_b2_ay,   gradleg_quad_p10p2_b2_ay,   gradleg_quad_p10p3_b2_ay,   gradleg_quad_p10p4_b2_ay,   gradleg_quad_p10p5_b2_ay,   gradleg_quad_p10p6_b2_ay,   gradleg_quad_p10p7_b2_ay,   gradleg_quad_p10p8_b2_ay,   gradleg_quad_p10p9_b2_ay,   gradleg_quad_p10p10_b2_ay,   gradleg_quad_p11p0_b2_ay,   gradleg_quad_p11p1_b2_ay,   gradleg_quad_p11p2_b2_ay,   gradleg_quad_p11p3_b2_ay,   gradleg_quad_p11p4_b2_ay,   gradleg_quad_p11p5_b2_ay,   gradleg_quad_p11p6_b2_ay,   gradleg_quad_p11p7_b2_ay,   gradleg_quad_p11p8_b2_ay,   gradleg_quad_p11p9_b2_ay,   gradleg_quad_p11p10_b2_ay, };

    static Shapeset::shape_fn_t gradleg_quad_fn_by[] =
    {
      gradleg_quad_p0_e1_by, gradleg_quad_p0_e1_by, gradleg_quad_p0_e2_by_0, gradleg_quad_p0_e2_by_1,  gradleg_quad_p0_e3_by, gradleg_quad_p0_e3_by, gradleg_quad_p0_e4_by_0, gradleg_quad_p0_e4_by_1,
      gradleg_quad_l2_l0_by, gradleg_quad_l2_l0_by, gradleg_quad_l1_l2_by, gradleg_quad_l1_l2_by, gradleg_quad_l2_l1_by, gradleg_quad_l2_l1_by, gradleg_quad_l0_l2_by,  gradleg_quad_l0_l2_by,
      gradleg_quad_l3_l0_by_0, gradleg_quad_l3_l0_by_1, gradleg_quad_l1_l3_by_0, gradleg_quad_l1_l3_by_1, gradleg_quad_l3_l1_by_0, gradleg_quad_l3_l1_by_1, gradleg_quad_l0_l3_by_0, gradleg_quad_l0_l3_by_1,
      gradleg_quad_l4_l0_by, gradleg_quad_l4_l0_by, gradleg_quad_l1_l4_by, gradleg_quad_l1_l4_by, gradleg_quad_l4_l1_by, gradleg_quad_l4_l1_by, gradleg_quad_l0_l4_by,  gradleg_quad_l0_l4_by,
      gradleg_quad_l5_l0_by_0, gradleg_quad_l5_l0_by_1, gradleg_quad_l1_l5_by_0, gradleg_quad_l1_l5_by_1, gradleg_quad_l5_l1_by_0, gradleg_quad_l5_l1_by_1, gradleg_quad_l0_l5_by_0, gradleg_quad_l0_l5_by_1,
      gradleg_quad_l6_l0_by, gradleg_quad_l6_l0_by, gradleg_quad_l1_l6_by, gradleg_quad_l1_l6_by, gradleg_quad_l6_l1_by, gradleg_quad_l6_l1_by, gradleg_quad_l0_l6_by,  gradleg_quad_l0_l6_by,
      gradleg_quad_l7_l0_by_0, gradleg_quad_l7_l0_by_1, gradleg_quad_l1_l7_by_0, gradleg_quad_l1_l7_by_1, gradleg_quad_l7_l1_by_0, gradleg_quad_l7_l1_by_1, gradleg_quad_l0_l7_by_0, gradleg_quad_l0_l7_by_1,
      gradleg_quad_l8_l0_by, gradleg_quad_l8_l0_by, gradleg_quad_l1_l8_by, gradleg_quad_l1_l8_by, gradleg_quad_l8_l1_by, gradleg_quad_l8_l1_by, gradleg_quad_l0_l8_by,  gradleg_quad_l0_l8_by,
      gradleg_quad_l9_l0_by_0, gradleg_quad_l9_l0_by_1, gradleg_quad_l1_l9_by_0, gradleg_quad_l1_l9_by_1, gradleg_quad_l9_l1_by_0, gradleg_quad_l9_l1_by_1, gradleg_quad_l0_l9_by_0, gradleg_quad_l0_l9_by_1,
      gradleg_quad_l10_l0_by, gradleg_quad_l10_l0_by, gradleg_quad_l1_l10_by, gradleg_quad_l1_l10_by, gradleg_quad_l10_l1_by, gradleg_quad_l10_l1_by, gradleg_quad_l0_l10_by,  gradleg_quad_l0_l10_by,
      gradleg_quad_l11_l0_by_0, gradleg_quad_l11_l0_by_1, gradleg_quad_l1_l11_by_0, gradleg_quad_l1_l11_by_1, gradleg_quad_l11_l1_by_0, gradleg_quad_l11_l1_by_1, gradleg_quad_l0_l11_by_0, gradleg_quad_l0_l11_by_1,

      gradleg_quad_p0p2_b1_by,   gradleg_quad_p0p3_b1_by,   gradleg_quad_p0p4_b1_by,   gradleg_quad_p0p5_b1_by,   gradleg_quad_p0p6_b1_by,   gradleg_quad_p0p7_b1_by,   gradleg_quad_p0p8_b1_by,   gradleg_quad_p0p9_b1_by,   gradleg_quad_p0p10_b1_by,   gradleg_quad_p0p11_b1_by,   gradleg_quad_p1p2_b1_by,   gradleg_quad_p1p3_b1_by,   gradleg_quad_p1p4_b1_by,   gradleg_quad_p1p5_b1_by,   gradleg_quad_p1p6_b1_by,   gradleg_quad_p1p7_b1_by,   gradleg_quad_p1p8_b1_by,   gradleg_quad_p1p9_b1_by,   gradleg_quad_p1p10_b1_by,   gradleg_quad_p1p11_b1_by,   gradleg_quad_p2p2_b1_by,   gradleg_quad_p2p3_b1_by,   gradleg_quad_p2p4_b1_by,   gradleg_quad_p2p5_b1_by,   gradleg_quad_p2p6_b1_by,   gradleg_quad_p2p7_b1_by,   gradleg_quad_p2p8_b1_by,   gradleg_quad_p2p9_b1_by,   gradleg_quad_p2p10_b1_by,   gradleg_quad_p2p11_b1_by,   gradleg_quad_p3p2_b1_by,   gradleg_quad_p3p3_b1_by,   gradleg_quad_p3p4_b1_by,   gradleg_quad_p3p5_b1_by,   gradleg_quad_p3p6_b1_by,   gradleg_quad_p3p7_b1_by,   gradleg_quad_p3p8_b1_by,   gradleg_quad_p3p9_b1_by,   gradleg_quad_p3p10_b1_by,   gradleg_quad_p3p11_b1_by,   gradleg_quad_p4p2_b1_by,   gradleg_quad_p4p3_b1_by,   gradleg_quad_p4p4_b1_by,   gradleg_quad_p4p5_b1_by,   gradleg_quad_p4p6_b1_by,   gradleg_quad_p4p7_b1_by,   gradleg_quad_p4p8_b1_by,   gradleg_quad_p4p9_b1_by,   gradleg_quad_p4p10_b1_by,   gradleg_quad_p4p11_b1_by,   gradleg_quad_p5p2_b1_by,   gradleg_quad_p5p3_b1_by,   gradleg_quad_p5p4_b1_by,   gradleg_quad_p5p5_b1_by,   gradleg_quad_p5p6_b1_by,   gradleg_quad_p5p7_b1_by,   gradleg_quad_p5p8_b1_by,   gradleg_quad_p5p9_b1_by,   gradleg_quad_p5p10_b1_by,   gradleg_quad_p5p11_b1_by,   gradleg_quad_p6p2_b1_by,   gradleg_quad_p6p3_b1_by,   gradleg_quad_p6p4_b1_by,   gradleg_quad_p6p5_b1_by,   gradleg_quad_p6p6_b1_by,   gradleg_quad_p6p7_b1_by,   gradleg_quad_p6p8_b1_by,   gradleg_quad_p6p9_b1_by,   gradleg_quad_p6p10_b1_by,   gradleg_quad_p6p11_b1_by,   gradleg_quad_p7p2_b1_by,   gradleg_quad_p7p3_b1_by,   gradleg_quad_p7p4_b1_by,   gradleg_quad_p7p5_b1_by,   gradleg_quad_p7p6_b1_by,   gradleg_quad_p7p7_b1_by,   gradleg_quad_p7p8_b1_by,   gradleg_quad_p7p9_b1_by,   gradleg_quad_p7p10_b1_by,   gradleg_quad_p7p11_b1_by,   gradleg_quad_p8p2_b1_by,   gradleg_quad_p8p3_b1_by,   gradleg_quad_p8p4_b1_by,   gradleg_quad_p8p5_b1_by,   gradleg_quad_p8p6_b1_by,   gradleg_quad_p8p7_b1_by,   gradleg_quad_p8p8_b1_by,   gradleg_quad_p8p9_b1_by,   gradleg_quad_p8p10_b1_by,   gradleg_quad_p8p11_b1_by,   gradleg_quad_p9p2_b1_by,   gradleg_quad_p9p3_b1_by,   gradleg_quad_p9p4_b1_by,   gradleg_quad_p9p5_b1_by,   gradleg_quad_p9p6_b1_by,   gradleg_quad_p9p7_b1_by,   gradleg_quad_p9p8_b1_by,   gradleg_quad_p9p9_b1_by,   gradleg_quad_p9p10_b1_by,   gradleg_quad_p9p11_b1_by,   gradleg_quad_p10p2_b1_by,   gradleg_quad_p10p3_b1_by,   gradleg_quad_p10p4_b1_by,   gradleg_quad_p10p5_b1_by,   gradleg_quad_p10p6_b1_by,   gradleg_quad_p10p7_b1_by,   gradleg_quad_p10p8_b1_by,   gradleg_quad_p10p9_b1_by,   gradleg_quad_p10p10_b1_by,   gradleg_quad_p10p11_b1_by,   gradleg_quad_p2p0_b2_by,   gradleg_quad_p2p1_b2_by,   gradleg_quad_p2p2_b2_by,   gradleg_quad_p2p3_b2_by,   gradleg_quad_p2p4_b2_by,   gradleg_quad_p2p5_b2_by,   gradleg_quad_p2p6_b2_by,   gradleg_quad_p2p7_b2_by,   gradleg_quad_p2p8_b2_by,   gradleg_quad_p2p9_b2_by,   gradleg_quad_p2p10_b2_by,   gradleg_quad_p3p0_b2_by,   gradleg_quad_p3p1_b2_by,   gradleg_quad_p3p2_b2_by,   gradleg_quad_p3p3_b2_by,   gradleg_quad_p3p4_b2_by,   gradleg_quad_p3p5_b2_by,   gradleg_quad_p3p6_b2_by,   gradleg_quad_p3p7_b2_by,   gradleg_quad_p3p8_b2_by,   gradleg_quad_p3p9_b2_by,   gradleg_quad_p3p10_b2_by,   gradleg_quad_p4p0_b2_by,   gradleg_quad_p4p1_b2_by,   gradleg_quad_p4p2_b2_by,   gradleg_quad_p4p3_b2_by,   gradleg_quad_p4p4_b2_by,   gradleg_quad_p4p5_b2_by,   gradleg_quad_p4p6_b2_by,   gradleg_quad_p4p7_b2_by,   gradleg_quad_p4p8_b2_by,   gradleg_quad_p4p9_b2_by,   gradleg_quad_p4p10_b2_by,   gradleg_quad_p5p0_b2_by,   gradleg_quad_p5p1_b2_by,   gradleg_quad_p5p2_b2_by,   gradleg_quad_p5p3_b2_by,   gradleg_quad_p5p4_b2_by,   gradleg_quad_p5p5_b2_by,   gradleg_quad_p5p6_b2_by,   gradleg_quad_p5p7_b2_by,   gradleg_quad_p5p8_b2_by,   gradleg_quad_p5p9_b2_by,   gradleg_quad_p5p10_b2_by,   gradleg_quad_p6p0_b2_by,   gradleg_quad_p6p1_b2_by,   gradleg_quad_p6p2_b2_by,   gradleg_quad_p6p3_b2_by,   gradleg_quad_p6p4_b2_by,   gradleg_quad_p6p5_b2_by,   gradleg_quad_p6p6_b2_by,   gradleg_quad_p6p7_b2_by,   gradleg_quad_p6p8_b2_by,   gradleg_quad_p6p9_b2_by,   gradleg_quad_p6p10_b2_by,   gradleg_quad_p7p0_b2_by,   gradleg_quad_p7p1_b2_by,   gradleg_quad_p7p2_b2_by,   gradleg_quad_p7p3_b2_by,   gradleg_quad_p7p4_b2_by,   gradleg_quad_p7p5_b2_by,   gradleg_quad_p7p6_b2_by,   gradleg_quad_p7p7_b2_by,   gradleg_quad_p7p8_b2_by,   gradleg_quad_p7p9_b2_by,   gradleg_quad_p7p10_b2_by,   gradleg_quad_p8p0_b2_by,   gradleg_quad_p8p1_b2_by,   gradleg_quad_p8p2_b2_by,   gradleg_quad_p8p3_b2_by,   gradleg_quad_p8p4_b2_by,   gradleg_quad_p8p5_b2_by,   gradleg_quad_p8p6_b2_by,   gradleg_quad_p8p7_b2_by,   gradleg_quad_p8p8_b2_by,   gradleg_quad_p8p9_b2_by,   gradleg_quad_p8p10_b2_by,   gradleg_quad_p9p0_b2_by,   gradleg_quad_p9p1_b2_by,   gradleg_quad_p9p2_b2_by,   gradleg_quad_p9p3_b2_by,   gradleg_quad_p9p4_b2_by,   gradleg_quad_p9p5_b2_by,   gradleg_quad_p9p6_b2_by,   gradleg_quad_p9p7_b2_by,   gradleg_quad_p9p8_b2_by,   gradleg_quad_p9p9_b2_by,   gradleg_quad_p9p10_b2_by,   gradleg_quad_p10p0_b2_by,   gradleg_quad_p10p1_b2_by,   gradleg_quad_p10p2_b2_by,   gradleg_quad_p10p3_b2_by,   gradleg_quad_p10p4_b2_by,   gradleg_quad_p10p5_b2_by,   gradleg_quad_p10p6_b2_by,   gradleg_quad_p10p7_b2_by,   gradleg_quad_p10p8_b2_by,   gradleg_quad_p10p9_b2_by,   gradleg_quad_p10p10_b2_by,   gradleg_quad_p11p0_b2_by,   gradleg_quad_p11p1_b2_by,   gradleg_quad_p11p2_b2_by,   gradleg_quad_p11p3_b2_by,   gradleg_quad_p11p4_b2_by,   gradleg_quad_p11p5_b2_by,   gradleg_quad_p11p6_b2_by,   gradleg_quad_p11p7_b2_by,   gradleg_quad_p11p8_b2_by,   gradleg_quad_p11p9_b2_by,   gradleg_quad_p11p10_b2_by, };

    static int qb_0_1[] = { 88, };
      static int qb_0_2[] = { 88, 89, };
      static int qb_0_3[] = { 88, 89, 90, };
      static int qb_0_4[] = { 88, 89, 90, 91, };
      static int qb_0_5[] = { 88, 89, 90, 91, 92, };
      static int qb_0_6[] = { 88, 89, 90, 91, 92, 93, };
      static int qb_0_7[] = { 88, 89, 90, 91, 92, 93, 94, };
      static int qb_0_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, };
      static int qb_0_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, };
      static int qb_0_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, };
      static int qb_1_0[] = { 198, };
      static int qb_1_1[] = { 88, 198, 98, 199, };
      static int qb_1_2[] = { 88, 89, 198, 98, 199, 99, 200, };
      static int qb_1_3[] = { 88, 89, 90, 198, 98, 199, 99, 200, 100, 201, };
      static int qb_1_4[] = { 88, 89, 90, 91, 198, 98, 199, 99, 200, 100, 201, 101, 202, };
      static int qb_1_5[] = { 88, 89, 90, 91, 92, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, };
      static int qb_1_6[] = { 88, 89, 90, 91, 92, 93, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, };
      static int qb_1_7[] = { 88, 89, 90, 91, 92, 93, 94, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, };
      static int qb_1_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, };
      static int qb_1_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, };
      static int qb_1_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 107, 208, };
      static int qb_2_0[] = { 198, 209, };
      static int qb_2_1[] = { 88, 198, 98, 199, 209, 108, 210, };
      static int qb_2_2[] = { 88, 89, 198, 98, 199, 99, 200, 209, 108, 210, 109, 211, };
      static int qb_2_3[] = { 88, 89, 90, 198, 98, 199, 99, 200, 100, 201, 209, 108, 210, 109, 211, 110, 212, };
      static int qb_2_4[] = { 88, 89, 90, 91, 198, 98, 199, 99, 200, 100, 201, 101, 202, 209, 108, 210, 109, 211, 110, 212, 111, 213, };
      static int qb_2_5[] = { 88, 89, 90, 91, 92, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, };
      static int qb_2_6[] = { 88, 89, 90, 91, 92, 93, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, };
      static int qb_2_7[] = { 88, 89, 90, 91, 92, 93, 94, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, };
      static int qb_2_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, };
      static int qb_2_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, };
      static int qb_2_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 107, 208, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 117, 219, };
      static int qb_3_0[] = { 198, 209, 220, };
      static int qb_3_1[] = { 88, 198, 98, 199, 209, 108, 210, 220, 118, 221, };
      static int qb_3_2[] = { 88, 89, 198, 98, 199, 99, 200, 209, 108, 210, 109, 211, 220, 118, 221, 119, 222, };
      static int qb_3_3[] = { 88, 89, 90, 198, 98, 199, 99, 200, 100, 201, 209, 108, 210, 109, 211, 110, 212, 220, 118, 221, 119, 222, 120, 223, };
      static int qb_3_4[] = { 88, 89, 90, 91, 198, 98, 199, 99, 200, 100, 201, 101, 202, 209, 108, 210, 109, 211, 110, 212, 111, 213, 220, 118, 221, 119, 222, 120, 223, 121, 224, };
      static int qb_3_5[] = { 88, 89, 90, 91, 92, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, };
      static int qb_3_6[] = { 88, 89, 90, 91, 92, 93, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, };
      static int qb_3_7[] = { 88, 89, 90, 91, 92, 93, 94, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, };
      static int qb_3_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, };
      static int qb_3_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, };
      static int qb_3_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 107, 208, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 117, 219, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 127, 230, };
      static int qb_4_0[] = { 198, 209, 220, 231, };
      static int qb_4_1[] = { 88, 198, 98, 199, 209, 108, 210, 220, 118, 221, 231, 128, 232, };
      static int qb_4_2[] = { 88, 89, 198, 98, 199, 99, 200, 209, 108, 210, 109, 211, 220, 118, 221, 119, 222, 231, 128, 232, 129, 233, };
      static int qb_4_3[] = { 88, 89, 90, 198, 98, 199, 99, 200, 100, 201, 209, 108, 210, 109, 211, 110, 212, 220, 118, 221, 119, 222, 120, 223, 231, 128, 232, 129, 233, 130, 234, };
      static int qb_4_4[] = { 88, 89, 90, 91, 198, 98, 199, 99, 200, 100, 201, 101, 202, 209, 108, 210, 109, 211, 110, 212, 111, 213, 220, 118, 221, 119, 222, 120, 223, 121, 224, 231, 128, 232, 129, 233, 130, 234, 131, 235, };
      static int qb_4_5[] = { 88, 89, 90, 91, 92, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, };
      static int qb_4_6[] = { 88, 89, 90, 91, 92, 93, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, };
      static int qb_4_7[] = { 88, 89, 90, 91, 92, 93, 94, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, };
      static int qb_4_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, };
      static int qb_4_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, };
      static int qb_4_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 107, 208, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 117, 219, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 127, 230, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 137, 241, };
      static int qb_5_0[] = { 198, 209, 220, 231, 242, };
      static int qb_5_1[] = { 88, 198, 98, 199, 209, 108, 210, 220, 118, 221, 231, 128, 232, 242, 138, 243, };
      static int qb_5_2[] = { 88, 89, 198, 98, 199, 99, 200, 209, 108, 210, 109, 211, 220, 118, 221, 119, 222, 231, 128, 232, 129, 233, 242, 138, 243, 139, 244, };
      static int qb_5_3[] = { 88, 89, 90, 198, 98, 199, 99, 200, 100, 201, 209, 108, 210, 109, 211, 110, 212, 220, 118, 221, 119, 222, 120, 223, 231, 128, 232, 129, 233, 130, 234, 242, 138, 243, 139, 244, 140, 245, };
      static int qb_5_4[] = { 88, 89, 90, 91, 198, 98, 199, 99, 200, 100, 201, 101, 202, 209, 108, 210, 109, 211, 110, 212, 111, 213, 220, 118, 221, 119, 222, 120, 223, 121, 224, 231, 128, 232, 129, 233, 130, 234, 131, 235, 242, 138, 243, 139, 244, 140, 245, 141, 246, };
      static int qb_5_5[] = { 88, 89, 90, 91, 92, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, };
      static int qb_5_6[] = { 88, 89, 90, 91, 92, 93, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, };
      static int qb_5_7[] = { 88, 89, 90, 91, 92, 93, 94, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, };
      static int qb_5_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, };
      static int qb_5_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, };
      static int qb_5_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 107, 208, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 117, 219, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 127, 230, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 137, 241, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 147, 252, };
      static int qb_6_0[] = { 198, 209, 220, 231, 242, 253, };
      static int qb_6_1[] = { 88, 198, 98, 199, 209, 108, 210, 220, 118, 221, 231, 128, 232, 242, 138, 243, 253, 148, 254, };
      static int qb_6_2[] = { 88, 89, 198, 98, 199, 99, 200, 209, 108, 210, 109, 211, 220, 118, 221, 119, 222, 231, 128, 232, 129, 233, 242, 138, 243, 139, 244, 253, 148, 254, 149, 255, };
      static int qb_6_3[] = { 88, 89, 90, 198, 98, 199, 99, 200, 100, 201, 209, 108, 210, 109, 211, 110, 212, 220, 118, 221, 119, 222, 120, 223, 231, 128, 232, 129, 233, 130, 234, 242, 138, 243, 139, 244, 140, 245, 253, 148, 254, 149, 255, 150, 256, };
      static int qb_6_4[] = { 88, 89, 90, 91, 198, 98, 199, 99, 200, 100, 201, 101, 202, 209, 108, 210, 109, 211, 110, 212, 111, 213, 220, 118, 221, 119, 222, 120, 223, 121, 224, 231, 128, 232, 129, 233, 130, 234, 131, 235, 242, 138, 243, 139, 244, 140, 245, 141, 246, 253, 148, 254, 149, 255, 150, 256, 151, 257, };
      static int qb_6_5[] = { 88, 89, 90, 91, 92, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, };
      static int qb_6_6[] = { 88, 89, 90, 91, 92, 93, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, };
      static int qb_6_7[] = { 88, 89, 90, 91, 92, 93, 94, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, };
      static int qb_6_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, };
      static int qb_6_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 156, 262, };
      static int qb_6_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 107, 208, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 117, 219, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 127, 230, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 137, 241, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 147, 252, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 156, 262, 157, 263, };
      static int qb_7_0[] = { 198, 209, 220, 231, 242, 253, 264, };
      static int qb_7_1[] = { 88, 198, 98, 199, 209, 108, 210, 220, 118, 221, 231, 128, 232, 242, 138, 243, 253, 148, 254, 264, 158, 265, };
      static int qb_7_2[] = { 88, 89, 198, 98, 199, 99, 200, 209, 108, 210, 109, 211, 220, 118, 221, 119, 222, 231, 128, 232, 129, 233, 242, 138, 243, 139, 244, 253, 148, 254, 149, 255, 264, 158, 265, 159, 266, };
      static int qb_7_3[] = { 88, 89, 90, 198, 98, 199, 99, 200, 100, 201, 209, 108, 210, 109, 211, 110, 212, 220, 118, 221, 119, 222, 120, 223, 231, 128, 232, 129, 233, 130, 234, 242, 138, 243, 139, 244, 140, 245, 253, 148, 254, 149, 255, 150, 256, 264, 158, 265, 159, 266, 160, 267, };
      static int qb_7_4[] = { 88, 89, 90, 91, 198, 98, 199, 99, 200, 100, 201, 101, 202, 209, 108, 210, 109, 211, 110, 212, 111, 213, 220, 118, 221, 119, 222, 120, 223, 121, 224, 231, 128, 232, 129, 233, 130, 234, 131, 235, 242, 138, 243, 139, 244, 140, 245, 141, 246, 253, 148, 254, 149, 255, 150, 256, 151, 257, 264, 158, 265, 159, 266, 160, 267, 161, 268, };
      static int qb_7_5[] = { 88, 89, 90, 91, 92, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, };
      static int qb_7_6[] = { 88, 89, 90, 91, 92, 93, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, };
      static int qb_7_7[] = { 88, 89, 90, 91, 92, 93, 94, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, };
      static int qb_7_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, };
      static int qb_7_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 156, 262, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 166, 273, };
      static int qb_7_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 107, 208, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 117, 219, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 127, 230, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 137, 241, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 147, 252, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 156, 262, 157, 263, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 166, 273, 167, 274, };
      static int qb_8_0[] = { 198, 209, 220, 231, 242, 253, 264, 275, };
      static int qb_8_1[] = { 88, 198, 98, 199, 209, 108, 210, 220, 118, 221, 231, 128, 232, 242, 138, 243, 253, 148, 254, 264, 158, 265, 275, 168, 276, };
      static int qb_8_2[] = { 88, 89, 198, 98, 199, 99, 200, 209, 108, 210, 109, 211, 220, 118, 221, 119, 222, 231, 128, 232, 129, 233, 242, 138, 243, 139, 244, 253, 148, 254, 149, 255, 264, 158, 265, 159, 266, 275, 168, 276, 169, 277, };
      static int qb_8_3[] = { 88, 89, 90, 198, 98, 199, 99, 200, 100, 201, 209, 108, 210, 109, 211, 110, 212, 220, 118, 221, 119, 222, 120, 223, 231, 128, 232, 129, 233, 130, 234, 242, 138, 243, 139, 244, 140, 245, 253, 148, 254, 149, 255, 150, 256, 264, 158, 265, 159, 266, 160, 267, 275, 168, 276, 169, 277, 170, 278, };
      static int qb_8_4[] = { 88, 89, 90, 91, 198, 98, 199, 99, 200, 100, 201, 101, 202, 209, 108, 210, 109, 211, 110, 212, 111, 213, 220, 118, 221, 119, 222, 120, 223, 121, 224, 231, 128, 232, 129, 233, 130, 234, 131, 235, 242, 138, 243, 139, 244, 140, 245, 141, 246, 253, 148, 254, 149, 255, 150, 256, 151, 257, 264, 158, 265, 159, 266, 160, 267, 161, 268, 275, 168, 276, 169, 277, 170, 278, 171, 279, };
      static int qb_8_5[] = { 88, 89, 90, 91, 92, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, };
      static int qb_8_6[] = { 88, 89, 90, 91, 92, 93, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, };
      static int qb_8_7[] = { 88, 89, 90, 91, 92, 93, 94, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, };
      static int qb_8_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 175, 283, };
      static int qb_8_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 156, 262, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 166, 273, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 175, 283, 176, 284, };
      static int qb_8_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 107, 208, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 117, 219, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 127, 230, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 137, 241, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 147, 252, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 156, 262, 157, 263, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 166, 273, 167, 274, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 175, 283, 176, 284, 177, 285, };
      static int qb_9_0[] = { 198, 209, 220, 231, 242, 253, 264, 275, 286, };
      static int qb_9_1[] = { 88, 198, 98, 199, 209, 108, 210, 220, 118, 221, 231, 128, 232, 242, 138, 243, 253, 148, 254, 264, 158, 265, 275, 168, 276, 286, 178, 287, };
      static int qb_9_2[] = { 88, 89, 198, 98, 199, 99, 200, 209, 108, 210, 109, 211, 220, 118, 221, 119, 222, 231, 128, 232, 129, 233, 242, 138, 243, 139, 244, 253, 148, 254, 149, 255, 264, 158, 265, 159, 266, 275, 168, 276, 169, 277, 286, 178, 287, 179, 288, };
      static int qb_9_3[] = { 88, 89, 90, 198, 98, 199, 99, 200, 100, 201, 209, 108, 210, 109, 211, 110, 212, 220, 118, 221, 119, 222, 120, 223, 231, 128, 232, 129, 233, 130, 234, 242, 138, 243, 139, 244, 140, 245, 253, 148, 254, 149, 255, 150, 256, 264, 158, 265, 159, 266, 160, 267, 275, 168, 276, 169, 277, 170, 278, 286, 178, 287, 179, 288, 180, 289, };
      static int qb_9_4[] = { 88, 89, 90, 91, 198, 98, 199, 99, 200, 100, 201, 101, 202, 209, 108, 210, 109, 211, 110, 212, 111, 213, 220, 118, 221, 119, 222, 120, 223, 121, 224, 231, 128, 232, 129, 233, 130, 234, 131, 235, 242, 138, 243, 139, 244, 140, 245, 141, 246, 253, 148, 254, 149, 255, 150, 256, 151, 257, 264, 158, 265, 159, 266, 160, 267, 161, 268, 275, 168, 276, 169, 277, 170, 278, 171, 279, 286, 178, 287, 179, 288, 180, 289, 181, 290, };
      static int qb_9_5[] = { 88, 89, 90, 91, 92, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, };
      static int qb_9_6[] = { 88, 89, 90, 91, 92, 93, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 183, 292, };
      static int qb_9_7[] = { 88, 89, 90, 91, 92, 93, 94, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 183, 292, 184, 293, };
      static int qb_9_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 175, 283, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 183, 292, 184, 293, 185, 294, };
      static int qb_9_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 156, 262, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 166, 273, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 175, 283, 176, 284, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 183, 292, 184, 293, 185, 294, 186, 295, };
      static int qb_9_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 107, 208, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 117, 219, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 127, 230, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 137, 241, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 147, 252, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 156, 262, 157, 263, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 166, 273, 167, 274, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 175, 283, 176, 284, 177, 285, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 183, 292, 184, 293, 185, 294, 186, 295, 187, 296, };
      static int qb_10_0[] = { 198, 209, 220, 231, 242, 253, 264, 275, 286, 297, };
      static int qb_10_1[] = { 88, 198, 98, 199, 209, 108, 210, 220, 118, 221, 231, 128, 232, 242, 138, 243, 253, 148, 254, 264, 158, 265, 275, 168, 276, 286, 178, 287, 297, 188, 298, };
      static int qb_10_2[] = { 88, 89, 198, 98, 199, 99, 200, 209, 108, 210, 109, 211, 220, 118, 221, 119, 222, 231, 128, 232, 129, 233, 242, 138, 243, 139, 244, 253, 148, 254, 149, 255, 264, 158, 265, 159, 266, 275, 168, 276, 169, 277, 286, 178, 287, 179, 288, 297, 188, 298, 189, 299, };
      static int qb_10_3[] = { 88, 89, 90, 198, 98, 199, 99, 200, 100, 201, 209, 108, 210, 109, 211, 110, 212, 220, 118, 221, 119, 222, 120, 223, 231, 128, 232, 129, 233, 130, 234, 242, 138, 243, 139, 244, 140, 245, 253, 148, 254, 149, 255, 150, 256, 264, 158, 265, 159, 266, 160, 267, 275, 168, 276, 169, 277, 170, 278, 286, 178, 287, 179, 288, 180, 289, 297, 188, 298, 189, 299, 190, 300, };
      static int qb_10_4[] = { 88, 89, 90, 91, 198, 98, 199, 99, 200, 100, 201, 101, 202, 209, 108, 210, 109, 211, 110, 212, 111, 213, 220, 118, 221, 119, 222, 120, 223, 121, 224, 231, 128, 232, 129, 233, 130, 234, 131, 235, 242, 138, 243, 139, 244, 140, 245, 141, 246, 253, 148, 254, 149, 255, 150, 256, 151, 257, 264, 158, 265, 159, 266, 160, 267, 161, 268, 275, 168, 276, 169, 277, 170, 278, 171, 279, 286, 178, 287, 179, 288, 180, 289, 181, 290, 297, 188, 298, 189, 299, 190, 300, 191, 301, };
      static int qb_10_5[] = { 88, 89, 90, 91, 92, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 297, 188, 298, 189, 299, 190, 300, 191, 301, 192, 302, };
      static int qb_10_6[] = { 88, 89, 90, 91, 92, 93, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 183, 292, 297, 188, 298, 189, 299, 190, 300, 191, 301, 192, 302, 193, 303, };
      static int qb_10_7[] = { 88, 89, 90, 91, 92, 93, 94, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 183, 292, 184, 293, 297, 188, 298, 189, 299, 190, 300, 191, 301, 192, 302, 193, 303, 194, 304, };
      static int qb_10_8[] = { 88, 89, 90, 91, 92, 93, 94, 95, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 175, 283, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 183, 292, 184, 293, 185, 294, 297, 188, 298, 189, 299, 190, 300, 191, 301, 192, 302, 193, 303, 194, 304, 195, 305, };
      static int qb_10_9[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 156, 262, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 166, 273, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 175, 283, 176, 284, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 183, 292, 184, 293, 185, 294, 186, 295, 297, 188, 298, 189, 299, 190, 300, 191, 301, 192, 302, 193, 303, 194, 304, 195, 305, 196, 306, };
      static int qb_10_10[] = { 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 198, 98, 199, 99, 200, 100, 201, 101, 202, 102, 203, 103, 204, 104, 205, 105, 206, 106, 207, 107, 208, 209, 108, 210, 109, 211, 110, 212, 111, 213, 112, 214, 113, 215, 114, 216, 115, 217, 116, 218, 117, 219, 220, 118, 221, 119, 222, 120, 223, 121, 224, 122, 225, 123, 226, 124, 227, 125, 228, 126, 229, 127, 230, 231, 128, 232, 129, 233, 130, 234, 131, 235, 132, 236, 133, 237, 134, 238, 135, 239, 136, 240, 137, 241, 242, 138, 243, 139, 244, 140, 245, 141, 246, 142, 247, 143, 248, 144, 249, 145, 250, 146, 251, 147, 252, 253, 148, 254, 149, 255, 150, 256, 151, 257, 152, 258, 153, 259, 154, 260, 155, 261, 156, 262, 157, 263, 264, 158, 265, 159, 266, 160, 267, 161, 268, 162, 269, 163, 270, 164, 271, 165, 272, 166, 273, 167, 274, 275, 168, 276, 169, 277, 170, 278, 171, 279, 172, 280, 173, 281, 174, 282, 175, 283, 176, 284, 177, 285, 286, 178, 287, 179, 288, 180, 289, 181, 290, 182, 291, 183, 292, 184, 293, 185, 294, 186, 295, 187, 296, 297, 188, 298, 189, 299, 190, 300, 191, 301, 192, 302, 193, 303, 194, 304, 195, 305, 196, 306, 197, 307, };

    #define NULL16 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,

    static int* gradleg_quad_bubble_indices[] =
    {
      NULL, qb_0_1, qb_0_2, qb_0_3, qb_0_4, qb_0_5, qb_0_6, qb_0_7, qb_0_8, qb_0_9, qb_0_10,  NULL, NULL, NULL, NULL, NULL, NULL16
      qb_1_0, qb_1_1, qb_1_2, qb_1_3, qb_1_4, qb_1_5, qb_1_6, qb_1_7, qb_1_8, qb_1_9, qb_1_10,  NULL, NULL, NULL, NULL, NULL, NULL16
      qb_2_0, qb_2_1, qb_2_2, qb_2_3, qb_2_4, qb_2_5, qb_2_6, qb_2_7, qb_2_8, qb_2_9, qb_2_10,  NULL, NULL, NULL, NULL, NULL, NULL16
      qb_3_0, qb_3_1, qb_3_2, qb_3_3, qb_3_4, qb_3_5, qb_3_6, qb_3_7, qb_3_8, qb_3_9, qb_3_10,  NULL, NULL, NULL, NULL, NULL, NULL16
      qb_4_0, qb_4_1, qb_4_2, qb_4_3, qb_4_4, qb_4_5, qb_4_6, qb_4_7, qb_4_8, qb_4_9, qb_4_10,  NULL, NULL, NULL, NULL, NULL, NULL16
      qb_5_0, qb_5_1, qb_5_2, qb_5_3, qb_5_4, qb_5_5, qb_5_6, qb_5_7, qb_5_8, qb_5_9, qb_5_10,  NULL, NULL, NULL, NULL, NULL, NULL16
      qb_6_0, qb_6_1, qb_6_2, qb_6_3, qb_6_4, qb_6_5, qb_6_6, qb_6_7, qb_6_8, qb_6_9, qb_6_10,  NULL, NULL, NULL, NULL, NULL, NULL16
      qb_7_0, qb_7_1, qb_7_2, qb_7_3, qb_7_4, qb_7_5, qb_7_6, qb_7_7, qb_7_8, qb_7_9, qb_7_10,  NULL, NULL, NULL, NULL, NULL, NULL16
      qb_8_0, qb_8_1, qb_8_2, qb_8_3, qb_8_4, qb_8_5, qb_8_6, qb_8_7, qb_8_8, qb_8_9, qb_8_10,  NULL, NULL, NULL, NULL, NULL, NULL16
      qb_9_0, qb_9_1, qb_9_2, qb_9_3, qb_9_4, qb_9_5, qb_9_6, qb_9_7, qb_9_8, qb_9_9, qb_9_10,  NULL, NULL, NULL, NULL, NULL, NULL16
      qb_10_0, qb_10_1, qb_10_2, qb_10_3, qb_10_4, qb_10_5, qb_10_6, qb_10_7, qb_10_8, qb_10_9, qb_10_10,  NULL, NULL, NULL, NULL, NULL, NULL16
    };

    #define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

    static int gradleg_quad_bubble_count[] =
    {
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  0,  0,  0,  0,  0, zero16
      1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 0,  0,  0,  0,  0, zero16
      2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 0,  0,  0,  0,  0, zero16
      3, 10, 17, 24, 31, 38, 45, 52, 59, 66, 73, 0,  0,  0,  0,  0, zero16
      4, 13, 22, 31, 40, 49, 58, 67, 76, 85, 94, 0,  0,  0,  0,  0, zero16
      5, 16, 27, 38, 49, 60, 71, 82, 93, 104, 115, 0,  0,  0,  0,  0, zero16
      6, 19, 32, 45, 58, 71, 84, 97, 110, 123, 136, 0,  0,  0,  0,  0, zero16
      7, 22, 37, 52, 67, 82, 97, 112, 127, 142, 157, 0,  0,  0,  0,  0, zero16
      8, 25, 42, 59, 76, 93, 110, 127, 144, 161, 178, 0,  0,  0,  0,  0, zero16
      9, 28, 47, 66, 85, 104, 123, 142, 161, 180, 199, 0,  0,  0,  0,  0, zero16
      10, 31, 52, 73, 94, 115, 136, 157, 178, 199, 220, 0,  0,  0,  0,  0, zero16
    };

    static int gradleg_quad_vertex_indices[4] ={-1, -1, -1, -1};

    static int gradleg_quad_edge_indices_0[] = { 0, 1, 8, 9, 16, 17, 24, 25, 32, 33, 40, 41, 48, 49, 56, 57, 64, 65, 72, 73, 80, 81, };

    static int gradleg_quad_edge_indices_1[] = { 2, 3, 10, 11, 18, 19, 26, 27, 34, 35, 42, 43, 50, 51, 58, 59, 66, 67, 74, 75, 82, 83, };

    static int gradleg_quad_edge_indices_2[] = { 4, 5, 12, 13, 20, 21, 28, 29, 36, 37, 44, 45, 52, 53, 60, 61, 68, 69, 76, 77, 84, 85, };

    static int gradleg_quad_edge_indices_3[] = { 6, 7, 14, 15, 22, 23, 30, 31, 38, 39, 46, 47, 54, 55, 62, 63, 70, 71, 78, 79, 86, 87, };

    static int* gradleg_quad_edge_indices[4] =
    {
      gradleg_quad_edge_indices_0,
      gradleg_quad_edge_indices_1,
      gradleg_quad_edge_indices_2,
      gradleg_quad_edge_indices_3,
    };

    #define oo H2D_MAKE_QUAD_ORDER

    static int gradleg_quad_index_to_order[] =
    {
      oo(0, 1), oo(0, 1), oo(1, 0), oo(1, 0), oo(0, 1), oo(0, 1), oo(1, 0), oo(1, 0),
      oo(2, 1), oo(2, 1), oo(1, 2), oo(1, 2), oo(2, 1), oo(2, 1), oo(1, 2), oo(1, 2),
      oo(3, 1), oo(3, 1), oo(1, 3), oo(1, 3), oo(3, 1), oo(3, 1), oo(1, 3), oo(1, 3),
      oo(4, 1), oo(4, 1), oo(1, 4), oo(1, 4), oo(4, 1), oo(4, 1), oo(1, 4), oo(1, 4),
      oo(5, 1), oo(5, 1), oo(1, 5), oo(1, 5), oo(5, 1), oo(5, 1), oo(1, 5), oo(1, 5),
      oo(6, 1), oo(6, 1), oo(1, 6), oo(1, 6), oo(6, 1), oo(6, 1), oo(1, 6), oo(1, 6),
      oo(7, 1), oo(7, 1), oo(1, 7), oo(1, 7), oo(7, 1), oo(7, 1), oo(1, 7), oo(1, 7),
      oo(8, 1), oo(8, 1), oo(1, 8), oo(1, 8), oo(8, 1), oo(8, 1), oo(1, 8), oo(1, 8),
      oo(9, 1), oo(9, 1), oo(1, 9), oo(1, 9), oo(9, 1), oo(9, 1), oo(1, 9), oo(1, 9),
      oo(10, 1), oo(10, 1), oo(1, 10), oo(1, 10), oo(10, 1), oo(10, 1), oo(1, 10), oo(1, 10),
      oo(11, 1), oo(11, 1), oo(1, 11), oo(1, 11), oo(11, 1), oo(11, 1), oo(1, 11), oo(1, 11),

      oo(0, 2),  oo(0, 3),  oo(0, 4),  oo(0, 5),  oo(0, 6),  oo(0, 7),  oo(0, 8),  oo(0, 9),  oo(0, 10),  oo(0, 11),  oo(1, 2),  oo(1, 3),  oo(1, 4),  oo(1, 5),  oo(1, 6),  oo(1, 7),  oo(1, 8),  oo(1, 9),  oo(1, 10),  oo(1, 11),  oo(2, 2),  oo(2, 3),  oo(2, 4),  oo(2, 5),  oo(2, 6),  oo(2, 7),  oo(2, 8),  oo(2, 9),  oo(2, 10),  oo(2, 11),  oo(3, 2),  oo(3, 3),  oo(3, 4),  oo(3, 5),  oo(3, 6),  oo(3, 7),  oo(3, 8),  oo(3, 9),  oo(3, 10),  oo(3, 11),  oo(4, 2),  oo(4, 3),  oo(4, 4),  oo(4, 5),  oo(4, 6),  oo(4, 7),  oo(4, 8),  oo(4, 9),  oo(4, 10),  oo(4, 11),  oo(5, 2),  oo(5, 3),  oo(5, 4),  oo(5, 5),  oo(5, 6),  oo(5, 7),  oo(5, 8),  oo(5, 9),  oo(5, 10),  oo(5, 11),  oo(6, 2),  oo(6, 3),  oo(6, 4),  oo(6, 5),  oo(6, 6),  oo(6, 7),  oo(6, 8),  oo(6, 9),  oo(6, 10),  oo(6, 11),  oo(7, 2),  oo(7, 3),  oo(7, 4),  oo(7, 5),  oo(7, 6),  oo(7, 7),  oo(7, 8),  oo(7, 9),  oo(7, 10),  oo(7, 11),  oo(8, 2),  oo(8, 3),  oo(8, 4),  oo(8, 5),  oo(8, 6),  oo(8, 7),  oo(8, 8),  oo(8, 9),  oo(8, 10),  oo(8, 11),  oo(9, 2),  oo(9, 3),  oo(9, 4),  oo(9, 5),  oo(9, 6),  oo(9, 7),  oo(9, 8),  oo(9, 9),  oo(9, 10),  oo(9, 11),  oo(10, 2),  oo(10, 3),  oo(10, 4),  oo(10, 5),  oo(10, 6),  oo(10, 7),  oo(10, 8),  oo(10, 9),  oo(10, 10),  oo(10, 11),  oo(2, 0),  oo(2, 1),  oo(2, 2),  oo(2, 3),  oo(2, 4),  oo(2, 5),  oo(2, 6),  oo(2, 7),  oo(2, 8),  oo(2, 9),  oo(2, 10),  oo(3, 0),  oo(3, 1),  oo(3, 2),  oo(3, 3),  oo(3, 4),  oo(3, 5),  oo(3, 6),  oo(3, 7),  oo(3, 8),  oo(3, 9),  oo(3, 10),  oo(4, 0),  oo(4, 1),  oo(4, 2),  oo(4, 3),  oo(4, 4),  oo(4, 5),  oo(4, 6),  oo(4, 7),  oo(4, 8),  oo(4, 9),  oo(4, 10),  oo(5, 0),  oo(5, 1),  oo(5, 2),  oo(5, 3),  oo(5, 4),  oo(5, 5),  oo(5, 6),  oo(5, 7),  oo(5, 8),  oo(5, 9),  oo(5, 10),  oo(6, 0),  oo(6, 1),  oo(6, 2),  oo(6, 3),  oo(6, 4),  oo(6, 5),  oo(6, 6),  oo(6, 7),  oo(6, 8),  oo(6, 9),  oo(6, 10),  oo(7, 0),  oo(7, 1),  oo(7, 2),  oo(7, 3),  oo(7, 4),  oo(7, 5),  oo(7, 6),  oo(7, 7),  oo(7, 8),  oo(7, 9),  oo(7, 10),  oo(8, 0),  oo(8, 1),  oo(8, 2),  oo(8, 3),  oo(8, 4),  oo(8, 5),  oo(8, 6),  oo(8, 7),  oo(8, 8),  oo(8, 9),  oo(8, 10),  oo(9, 0),  oo(9, 1),  oo(9, 2),  oo(9, 3),  oo(9, 4),  oo(9, 5),  oo(9, 6),  oo(9, 7),  oo(9, 8),  oo(9, 9),  oo(9, 10),  oo(10, 0),  oo(10, 1),  oo(10, 2),  oo(10, 3),  oo(10, 4),  oo(10, 5),  oo(10, 6),  oo(10, 7),  oo(10, 8),  oo(10, 9),  oo(10, 10),  oo(11, 0),  oo(11, 1),  oo(11, 2),  oo(11, 3),  oo(11, 4),  oo(11, 5),  oo(11, 6),  oo(11, 7),  oo(11, 8),  oo(11, 9),  oo(11, 10),
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Shapeset::shape_fn_t* gradleg_quad_shape_fn_table[2] =
    {
      gradleg_quad_fn_a,
      gradleg_quad_fn_b
    };

    static Shapeset::shape_fn_t* gradleg_quad_shape_fn_table_x[2] =
    {
      gradleg_quad_fn_ax,
      gradleg_quad_fn_bx
    };

    static Shapeset::shape_fn_t* gradleg_quad_shape_fn_table_y[2] =
    {
      gradleg_quad_fn_ay,
      gradleg_quad_fn_by
    };

    //// triangle tables and class constructor ///////////////////////////////////////////////

    static Shapeset::shape_fn_t** gradleg_shape_fn_table[2] =
    {
      gradleg_tri_shape_fn_table,
      gradleg_quad_shape_fn_table
    };

    static Shapeset::shape_fn_t** gradleg_shape_fn_table_x[2] =
    {
      gradleg_tri_shape_fn_table_x,
      gradleg_quad_shape_fn_table_x
    };

    static Shapeset::shape_fn_t** gradleg_shape_fn_table_y[2] =
    {
      gradleg_tri_shape_fn_table_y,
      gradleg_quad_shape_fn_table_y
    };

    static int* gradleg_vertex_indices[2] =
    {
      gradleg_tri_vertex_indices,
      gradleg_quad_vertex_indices
    };

    static int** gradleg_edge_indices[2] =
    {
      gradleg_tri_edge_indices,
      gradleg_quad_edge_indices
    };

    static int** gradleg_bubble_indices[2] =
    {
      gradleg_tri_bubble_indices,
      gradleg_quad_bubble_indices
    };

    static int* gradleg_bubble_count[2] =
    {
      gradleg_tri_bubble_count,
      gradleg_quad_bubble_count
    };

    static int* gradleg_index_to_order[2] =
    {
      gradleg_tri_index_to_order,
      gradleg_quad_index_to_order
    };

    void check_gradleg_tri(Shapeset* shapeset)
    {
      for (int i = 1; i <= 10; i++)
      {
        int nb = shapeset->get_num_bubbles(i, HERMES_MODE_TRIANGLE);
        if(nb != 3*(i-1) + (i-1)*(i-2))
          throw Hermes::Exceptions::Exception("Wrong bubble count");
      }

      int size_a  = sizeof(gradleg_tri_fn_a)  / sizeof(Shapeset::shape_fn_t);
      int size_b  = sizeof(gradleg_tri_fn_b)  / sizeof(Shapeset::shape_fn_t);
      int size_ax = sizeof(gradleg_tri_fn_ax) / sizeof(Shapeset::shape_fn_t);
      int size_bx = sizeof(gradleg_tri_fn_bx) / sizeof(Shapeset::shape_fn_t);
      int size_ay = sizeof(gradleg_tri_fn_ay) / sizeof(Shapeset::shape_fn_t);
      int size_by = sizeof(gradleg_tri_fn_by) / sizeof(Shapeset::shape_fn_t);

      if(size_a != size_b || size_a != size_ax || size_a != size_bx || size_a != size_ay || size_a != size_by)
        throw Hermes::Exceptions::Exception("Function tables dont have equal length.");

      if(size_a != gradleg_tri_bubble_indices[10][gradleg_tri_bubble_count[10]-1] + 1)
        throw Hermes::Exceptions::Exception("Bad index of last bubble");
    }

    HcurlShapesetGradLeg::HcurlShapesetGradLeg()
    {
      shape_table[0] = gradleg_shape_fn_table;
      shape_table[1] = gradleg_shape_fn_table_x;
      shape_table[2] = gradleg_shape_fn_table_y;
      shape_table[3] = NULL;
      shape_table[4] = NULL;
      shape_table[5] = NULL;

      vertex_indices = gradleg_vertex_indices;
      edge_indices = gradleg_edge_indices;
      bubble_indices = gradleg_bubble_indices;
      bubble_count = gradleg_bubble_count;
      index_to_order = gradleg_index_to_order;

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
      num_components = 2;

      max_index[0] = 164;
      max_index[1] = 307;

      ebias = 0;

      comb_table = NULL;

      check_gradleg_tri(this);
    }
  }
}