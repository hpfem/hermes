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
    // Shape functions for the curl operator for triangles, based on Legendre polynomials

    //Shape functions

    // ORDER 0

    // Edge functions, order 0 (Whitney functions)

    // number 1
    inline double leg_tri_f1_a0(double x, double y)
    {
      return (psi0e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f1_b0(double x, double y)
    {
      return (psi0e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f1_a1(double x, double y)
    {
      return -leg_tri_f1_a0(x, y);
    }

    inline double leg_tri_f1_b1(double x, double y)
    {
      return -leg_tri_f1_b0(x, y);
    }

    inline double leg_tri_f1_ax0(double x, double y)
    {
      return (psi0e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f1_ay0(double x, double y)
    {
      return (psi0e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f1_bx0(double x, double y)
    {
      return (psi0e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f1_by0(double x, double y)
    {
      return (psi0e1y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f1_ax1(double x, double y)
    {
     return -leg_tri_f1_ax0(x, y);
    }

    inline double leg_tri_f1_ay1(double x, double y)
    {
     return -leg_tri_f1_ay0(x, y);
    }

    inline double leg_tri_f1_bx1(double x, double y)
    {
     return -leg_tri_f1_bx0(x, y);
    }

    inline double leg_tri_f1_by1(double x, double y)
    {
     return -leg_tri_f1_by0(x, y);
    }

    // number 2
    inline double leg_tri_f2_a0(double x, double y)
    {
      return (psi0e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f2_b0(double x, double y)
    {
      return (psi0e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f2_a1(double x, double y)
    {
      return -leg_tri_f2_a0(x, y);
    }

    inline double leg_tri_f2_b1(double x, double y)
    {
      return -leg_tri_f2_b0(x, y);
    }

    inline double leg_tri_f2_ax0(double x, double y)
    {
      return (psi0e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f2_ay0(double x, double y)
    {
      return (psi0e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f2_bx0(double x, double y)
    {
      return (psi0e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f2_by0(double x, double y)
    {
      return (psi0e2y_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f2_ax1(double x, double y)
    {
     return -leg_tri_f2_ax0(x, y);
    }

    inline double leg_tri_f2_ay1(double x, double y)
    {
     return -leg_tri_f2_ay0(x, y);
    }

    inline double leg_tri_f2_bx1(double x, double y)
    {
     return -leg_tri_f2_bx0(x, y);
    }

    inline double leg_tri_f2_by1(double x, double y)
    {
     return -leg_tri_f2_by0(x, y);
    }

    // number 3
    inline double leg_tri_f3_a0(double x, double y)
    {
      return (psi0e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f3_b0(double x, double y)
    {
      return (psi0e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f3_a1(double x, double y)
    {
      return -leg_tri_f3_a0(x, y);
    }

    inline double leg_tri_f3_b1(double x, double y)
    {
      return -leg_tri_f3_b0(x, y);
    }

    inline double leg_tri_f3_ax0(double x, double y)
    {
      return (psi0e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f3_ay0(double x, double y)
    {
      return (psi0e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f3_bx0(double x, double y)
    {
      return (psi0e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f3_by0(double x, double y)
    {
      return (psi0e3y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f3_ax1(double x, double y)
    {
     return -leg_tri_f3_ax0(x, y);
    }

    inline double leg_tri_f3_ay1(double x, double y)
    {
     return -leg_tri_f3_ay0(x, y);
    }

    inline double leg_tri_f3_bx1(double x, double y)
    {
     return -leg_tri_f3_bx0(x, y);
    }

    inline double leg_tri_f3_by1(double x, double y)
    {
     return -leg_tri_f3_by0(x, y);
    }

    // ORDER 1

    // Edge functions, order 1

    // number 4
    inline double leg_tri_f4_a(double x, double y)
    {
      return (psi1e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f4_b(double x, double y)
    {
      return (psi1e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f4_ax(double x, double y)
    {
      return (psi1e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f4_ay(double x, double y)
    {
      return (psi1e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f4_bx(double x, double y)
    {
      return (psi1e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f4_by(double x, double y)
    {
      return (psi1e1y_2(x, y)) / 1.0;
    }

    // number 5
    inline double leg_tri_f5_a(double x, double y)
    {
      return (psi1e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f5_b(double x, double y)
    {
      return (psi1e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f5_ax(double x, double y)
    {
      return (psi1e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f5_ay(double x, double y)
    {
      return (psi1e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f5_bx(double x, double y)
    {
      return (psi1e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f5_by(double x, double y)
    {
      return (psi1e2y_2(x, y)) / 1.4142135623731;
    }

    // number 6
    inline double leg_tri_f6_a(double x, double y)
    {
      return (psi1e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f6_b(double x, double y)
    {
      return (psi1e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f6_ax(double x, double y)
    {
      return (psi1e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f6_ay(double x, double y)
    {
      return (psi1e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f6_bx(double x, double y)
    {
      return (psi1e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f6_by(double x, double y)
    {
      return (psi1e3y_2(x, y)) / 1.0;
    }

    // ORDER 2

    // Edge functions, order 2

    // number 7
    inline double leg_tri_f7_a0(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 3.0 / 2.0 * Legendre1(L3 - L2); q2 = 1.0 / 2.0 * Legendre0(L3 - L2);
      return (q1 * psi1e1_1(x, y) - q2 * psi0e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f7_b0(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 3.0 / 2.0 * Legendre1(L3 - L2); q2 = 1.0 / 2.0 * Legendre0(L3 - L2);
      return (q1 * psi1e1_2(x, y) - q2 * psi0e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f7_a1(double x, double y)
    {
      return -leg_tri_f7_a0(x, y);
    }

    inline double leg_tri_f7_b1(double x, double y)
    {
      return -leg_tri_f7_b0(x, y);
    }

    inline double leg_tri_f7_ax0(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 3.0 / 2.0 * Legendre1x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 3.0 / 2.0 * Legendre1(L3 - L2);
      cx = 1.0 / 2.0 * Legendre0x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 1.0 / 2.0 * Legendre0(L3 - L2);
      return (ax * psi1e1_1(x, y) + bx * psi1e1x_1(x, y) - cx * psi0e1_1(x, y) - dx * psi0e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f7_ay0(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 3.0 / 2.0 * Legendre1x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 3.0 / 2.0 * Legendre1(L3 - L2);
      cy = 1.0 / 2.0 * Legendre0x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 1.0 / 2.0 * Legendre0(L3 - L2);
      return (ay * psi1e1_1(x, y) + by * psi1e1y_1(x, y) - cy * psi0e1_1(x, y) - dy * psi0e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f7_bx0(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 3.0 / 2.0 * Legendre1x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 3.0 / 2.0 * Legendre1(L3 - L2);
      cx = 1.0 / 2.0 * Legendre0x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 1.0 / 2.0 * Legendre0(L3 - L2);
      return (ax * psi1e1_2(x, y) + bx * psi1e1x_2(x, y) - cx * psi0e1_2(x, y) - dx * psi0e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f7_by0(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 3.0 / 2.0 * Legendre1x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 3.0 / 2.0 * Legendre1(L3 - L2);
      cy = 1.0 / 2.0 * Legendre0x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 1.0 / 2.0 * Legendre0(L3 - L2);
      return (ay * psi1e1_2(x, y) + by * psi1e1y_2(x, y) - cy * psi0e1_2(x, y) - dy * psi0e1y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f7_ax1(double x, double y)
    {
      return -leg_tri_f7_ax0(x, y);
    }

    inline double leg_tri_f7_ay1(double x, double y)
    {
      return -leg_tri_f7_ay0(x, y);
    }

    inline double leg_tri_f7_bx1(double x, double y)
    {
      return -leg_tri_f7_bx0(x, y);
    }

    inline double leg_tri_f7_by1(double x, double y)
    {
      return -leg_tri_f7_by0(x, y);
    }

    // number 8
    inline double leg_tri_f8_a0(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 3.0 / 2.0 * Legendre1(L1 - L3); q2 = 1.0 / 2.0 * Legendre0(L1 - L3);
      return (q1 * psi1e2_1(x, y) - q2 * psi0e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f8_b0(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 3.0 / 2.0 * Legendre1(L1 - L3); q2 = 1.0 / 2.0 * Legendre0(L1 - L3);
      return (q1 * psi1e2_2(x, y) - q2 * psi0e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f8_a1(double x, double y)
    {
      return -leg_tri_f8_a0(x, y);
    }

    inline double leg_tri_f8_b1(double x, double y)
    {
      return -leg_tri_f8_b0(x, y);
    }

    inline double leg_tri_f8_ax0(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 3.0 / 2.0 * Legendre1x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 3.0 / 2.0 * Legendre1(L1 - L3);
      cx = 1.0 / 2.0 * Legendre0x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 1.0 / 2.0 * Legendre0(L1 - L3);
      return (ax * psi1e2_1(x, y) + bx * psi1e2x_1(x, y) - cx * psi0e2_1(x, y) - dx * psi0e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f8_ay0(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 3.0 / 2.0 * Legendre1x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 3.0 / 2.0 * Legendre1(L1 - L3);
      cy = 1.0 / 2.0 * Legendre0x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 1.0 / 2.0 * Legendre0(L1 - L3);
      return (ay * psi1e2_1(x, y) + by * psi1e2y_1(x, y) - cy * psi0e2_1(x, y) - dy * psi0e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f8_bx0(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 3.0 / 2.0 * Legendre1x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 3.0 / 2.0 * Legendre1(L1 - L3);
      cx = 1.0 / 2.0 * Legendre0x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 1.0 / 2.0 * Legendre0(L1 - L3);
      return (ax * psi1e2_2(x, y) + bx * psi1e2x_2(x, y) - cx * psi0e2_2(x, y) - dx * psi0e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f8_by0(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 3.0 / 2.0 * Legendre1x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 3.0 / 2.0 * Legendre1(L1 - L3);
      cy = 1.0 / 2.0 * Legendre0x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 1.0 / 2.0 * Legendre0(L1 - L3);
      return (ay * psi1e2_2(x, y) + by * psi1e2y_2(x, y) - cy * psi0e2_2(x, y) - dy * psi0e2y_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f8_ax1(double x, double y)
    {
      return -leg_tri_f8_ax0(x, y);
    }

    inline double leg_tri_f8_ay1(double x, double y)
    {
      return -leg_tri_f8_ay0(x, y);
    }

    inline double leg_tri_f8_bx1(double x, double y)
    {
      return -leg_tri_f8_bx0(x, y);
    }

    inline double leg_tri_f8_by1(double x, double y)
    {
      return -leg_tri_f8_by0(x, y);
    }

    // number 9
    inline double leg_tri_f9_a0(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 3.0 / 2.0 * Legendre1(L2 - L1); q2 = 1.0 / 2.0 * Legendre0(L2 - L1);
      return (q1 * psi1e3_1(x, y) - q2 * psi0e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f9_b0(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 3.0 / 2.0 * Legendre1(L2 - L1); q2 = 1.0 / 2.0 * Legendre0(L2 - L1);
      return (q1 * psi1e3_2(x, y) - q2 * psi0e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f9_a1(double x, double y)
    {
      return -leg_tri_f9_a0(x, y);
    }

    inline double leg_tri_f9_b1(double x, double y)
    {
      return -leg_tri_f9_b0(x, y);
    }

    inline double leg_tri_f9_ax0(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 3.0 / 2.0 * Legendre1x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 3.0 / 2.0 * Legendre1(L2 - L1);
      cx = 1.0 / 2.0 * Legendre0x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 1.0 / 2.0 * Legendre0(L2 - L1);
      return (ax * psi1e3_1(x, y) + bx * psi1e3x_1(x, y) - cx * psi0e3_1(x, y) - dx * psi0e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f9_ay0(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 3.0 / 2.0 * Legendre1x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 3.0 / 2.0 * Legendre1(L2 - L1);
      cy = 1.0 / 2.0 * Legendre0x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 1.0 / 2.0 * Legendre0(L2 - L1);
      return (ay * psi1e3_1(x, y) + by * psi1e3y_1(x, y) - cy * psi0e3_1(x, y) - dy * psi0e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f9_bx0(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 3.0 / 2.0 * Legendre1x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 3.0 / 2.0 * Legendre1(L2 - L1);
      cx = 1.0 / 2.0 * Legendre0x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 1.0 / 2.0 * Legendre0(L2 - L1);
      return (ax * psi1e3_2(x, y) + bx * psi1e3x_2(x, y) - cx * psi0e3_2(x, y) - dx * psi0e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f9_by0(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 3.0 / 2.0 * Legendre1x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 3.0 / 2.0 * Legendre1(L2 - L1);
      cy = 1.0 / 2.0 * Legendre0x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 1.0 / 2.0 * Legendre0(L2 - L1);
      return (ay * psi1e3_2(x, y) + by * psi1e3y_2(x, y) - cy * psi0e3_2(x, y) - dy * psi0e3y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f9_ax1(double x, double y)
    {
      return -leg_tri_f9_ax0(x, y);
    }

    inline double leg_tri_f9_ay1(double x, double y)
    {
      return -leg_tri_f9_ay0(x, y);
    }

    inline double leg_tri_f9_bx1(double x, double y)
    {
      return -leg_tri_f9_bx0(x, y);
    }

    inline double leg_tri_f9_by1(double x, double y)
    {
      return -leg_tri_f9_by0(x, y);
    }

    // Edge-base bubble functions (normal functions), order 2

    // number 10
    inline double leg_tri_f10_a(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre0(L3 - L2);
      return (k * n11) / 1.0;
    }

    inline double leg_tri_f10_b(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre0(L3 - L2);
      return (k * n12) / 1.0;
    }

    inline double leg_tri_f10_ax(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre0(L3 - L2);
      Legx = Legendre0x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n11) / 1.0;
    }

    inline double leg_tri_f10_ay(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy,  ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre0(L3 - L2);
      Legy = Legendre0x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n11) / 1.0;
    }

    inline double leg_tri_f10_bx(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre0(L3 - L2);
      Legx = Legendre0x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n12) / 1.0;
    }

    inline double leg_tri_f10_by(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy, ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre0(L3 - L2);
      Legy = Legendre0x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n12) / 1.0;
    }

    // number 11
    inline double leg_tri_f11_a(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre0(L1 - L3);
      return (k * n21) / 1.4142135623731 *2.0;
    }

    inline double leg_tri_f11_b(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre0(L1 - L3);
      return (k * n22) / 1.4142135623731 *2.0;
    }

    inline double leg_tri_f11_ax(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre0(L1 - L3);
      Legx = Legendre0x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n21) / 1.4142135623731 *2.0;
    }

    inline double leg_tri_f11_ay(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy,  ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre0(L1 - L3);
      Legy = Legendre0x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n21) / 1.4142135623731 *2.0;
    }

    inline double leg_tri_f11_bx(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre0(L1 - L3);
      Legx = Legendre0x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n22) / 1.4142135623731 *2.0;
    }

    inline double leg_tri_f11_by(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy, ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre0(L1 - L3);
      Legy = Legendre0x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n22) / 1.4142135623731 *2.0;
    }

    // number 12
    inline double leg_tri_f12_a(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre0(L2 - L1);
      return (k * n31) / 1.0;
    }

    inline double leg_tri_f12_b(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre0(L2 - L1);
      return (k * n32) / 1.0;
    }

    inline double leg_tri_f12_ax(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre0(L2 - L1);
      Legx = Legendre0x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n31) / 1.0;
    }

    inline double leg_tri_f12_ay(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy,  ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre0(L2 - L1);
      Legy = Legendre0x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n31) / 1.0;
    }

    inline double leg_tri_f12_bx(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre0(L2 - L1);
      Legx = Legendre0x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n32) / 1.0;
    }

    inline double leg_tri_f12_by(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy, ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre0(L2 - L1);
      Legy = Legendre0x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n32) / 1.0;
    }

    // ORDER 3

    // Edge functions, order 3

    // number 13
    inline double leg_tri_f13_a(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 5.0 / 3.0 * Legendre2(L3 - L2); q2 = 2.0 / 3.0 * Legendre1(L3 - L2);
      return (q1 * psi1e1_1(x, y) - q2 * psi0e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f13_b(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 5.0 / 3.0 * Legendre2(L3 - L2); q2 = 2.0 / 3.0 * Legendre1(L3 - L2);
      return (q1 * psi1e1_2(x, y) - q2 * psi0e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f13_ax(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 5.0 / 3.0 * Legendre2x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 5.0 / 3.0 * Legendre2(L3 - L2);
      cx = 2.0 / 3.0 * Legendre1x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 2.0 / 3.0 * Legendre1(L3 - L2);
      return (ax * psi1e1_1(x, y) + bx * psi1e1x_1(x, y) - cx * psi0e1_1(x, y) - dx * psi0e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f13_ay(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 5.0 / 3.0 * Legendre2x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 5.0 / 3.0 * Legendre2(L3 - L2);
      cy = 2.0 / 3.0 * Legendre1x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 2.0 / 3.0 * Legendre1(L3 - L2);
      return (ay * psi1e1_1(x, y) + by * psi1e1y_1(x, y) - cy * psi0e1_1(x, y) - dy * psi0e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f13_bx(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 5.0 / 3.0 * Legendre2x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 5.0 / 3.0 * Legendre2(L3 - L2);
      cx = 2.0 / 3.0 * Legendre1x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 2.0 / 3.0 * Legendre1(L3 - L2);
      return (ax * psi1e1_2(x, y) + bx * psi1e1x_2(x, y) - cx * psi0e1_2(x, y) - dx * psi0e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f13_by(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 5.0 / 3.0 * Legendre2x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 5.0 / 3.0 * Legendre2(L3 - L2);
      cy = 2.0 / 3.0 * Legendre1x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 2.0 / 3.0 * Legendre1(L3 - L2);
      return (ay * psi1e1_2(x, y) + by * psi1e1y_2(x, y) - cy * psi0e1_2(x, y) - dy * psi0e1y_2(x, y)) / 1.0;
    }

    // number 14
    inline double leg_tri_f14_a(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 5.0 / 3.0 * Legendre2(L1 - L3); q2 = 2.0 / 3.0 * Legendre1(L1 - L3);
      return (q1 * psi1e2_1(x, y) - q2 * psi0e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f14_b(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 5.0 / 3.0 * Legendre2(L1 - L3); q2 = 2.0 / 3.0 * Legendre1(L1 - L3);
      return (q1 * psi1e2_2(x, y) - q2 * psi0e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f14_ax(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 5.0 / 3.0 * Legendre2x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 5.0 / 3.0 * Legendre2(L1 - L3);
      cx = 2.0 / 3.0 * Legendre1x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 2.0 / 3.0 * Legendre1(L1 - L3);
      return (ax * psi1e2_1(x, y) + bx * psi1e2x_1(x, y) - cx * psi0e2_1(x, y) - dx * psi0e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f14_ay(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 5.0 / 3.0 * Legendre2x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 5.0 / 3.0 * Legendre2(L1 - L3);
      cy = 2.0 / 3.0 * Legendre1x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 2.0 / 3.0 * Legendre1(L1 - L3);
      return (ay * psi1e2_1(x, y) + by * psi1e2y_1(x, y) - cy * psi0e2_1(x, y) - dy * psi0e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f14_bx(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 5.0 / 3.0 * Legendre2x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 5.0 / 3.0 * Legendre2(L1 - L3);
      cx = 2.0 / 3.0 * Legendre1x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 2.0 / 3.0 * Legendre1(L1 - L3);
      return (ax * psi1e2_2(x, y) + bx * psi1e2x_2(x, y) - cx * psi0e2_2(x, y) - dx * psi0e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f14_by(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 5.0 / 3.0 * Legendre2x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 5.0 / 3.0 * Legendre2(L1 - L3);
      cy = 2.0 / 3.0 * Legendre1x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 2.0 / 3.0 * Legendre1(L1 - L3);
      return (ay * psi1e2_2(x, y) + by * psi1e2y_2(x, y) - cy * psi0e2_2(x, y) - dy * psi0e2y_2(x, y)) / 1.4142135623731;
    }

    // number 15
    inline double leg_tri_f15_a(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 5.0 / 3.0 * Legendre2(L2 - L1); q2 = 2.0 / 3.0 * Legendre1(L2 - L1);
      return (q1 * psi1e3_1(x, y) - q2 * psi0e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f15_b(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 5.0 / 3.0 * Legendre2(L2 - L1); q2 = 2.0 / 3.0 * Legendre1(L2 - L1);
      return (q1 * psi1e3_2(x, y) - q2 * psi0e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f15_ax(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 5.0 / 3.0 * Legendre2x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 5.0 / 3.0 * Legendre2(L2 - L1);
      cx = 2.0 / 3.0 * Legendre1x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 2.0 / 3.0 * Legendre1(L2 - L1);
      return (ax * psi1e3_1(x, y) + bx * psi1e3x_1(x, y) - cx * psi0e3_1(x, y) - dx * psi0e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f15_ay(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 5.0 / 3.0 * Legendre2x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 5.0 / 3.0 * Legendre2(L2 - L1);
      cy = 2.0 / 3.0 * Legendre1x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 2.0 / 3.0 * Legendre1(L2 - L1);
      return (ay * psi1e3_1(x, y) + by * psi1e3y_1(x, y) - cy * psi0e3_1(x, y) - dy * psi0e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f15_bx(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 5.0 / 3.0 * Legendre2x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 5.0 / 3.0 * Legendre2(L2 - L1);
      cx = 2.0 / 3.0 * Legendre1x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 2.0 / 3.0 * Legendre1(L2 - L1);
      return (ax * psi1e3_2(x, y) + bx * psi1e3x_2(x, y) - cx * psi0e3_2(x, y) - dx * psi0e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f15_by(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 5.0 / 3.0 * Legendre2x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 5.0 / 3.0 * Legendre2(L2 - L1);
      cy = 2.0 / 3.0 * Legendre1x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 2.0 / 3.0 * Legendre1(L2 - L1);
      return (ay * psi1e3_2(x, y) + by * psi1e3y_2(x, y) - cy * psi0e3_2(x, y) - dy * psi0e3y_2(x, y)) / 1.0;
    }

    // Edge-base bubble functions (normal functions), order 3

    // number 16
    inline double leg_tri_f16_a(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre1(L3 - L2);
      return (k * n11) / 1.0;
    }

    inline double leg_tri_f16_b(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre1(L3 - L2);
      return (k * n12) / 1.0;
    }

    inline double leg_tri_f16_ax(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre1(L3 - L2);
      Legx = Legendre1x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n11) / 1.0;
    }

    inline double leg_tri_f16_ay(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy,  ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre1(L3 - L2);
      Legy = Legendre1x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n11) / 1.0;
    }

    inline double leg_tri_f16_bx(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre1(L3 - L2);
      Legx = Legendre1x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n12) / 1.0;
    }

    inline double leg_tri_f16_by(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy, ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre1(L3 - L2);
      Legy = Legendre1x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n12) / 1.0;
    }

    // number 17
    inline double leg_tri_f17_a(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre1(L1 - L3);
      return (k * n21) / 1.4142135623731;
    }

    inline double leg_tri_f17_b(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre1(L1 - L3);
      return (k * n22) / 1.4142135623731;
    }

    inline double leg_tri_f17_ax(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre1(L1 - L3);
      Legx = Legendre1x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n21) / 1.4142135623731;
    }

    inline double leg_tri_f17_ay(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy,  ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre1(L1 - L3);
      Legy = Legendre1x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n21) / 1.4142135623731;
    }

    inline double leg_tri_f17_bx(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre1(L1 - L3);
      Legx = Legendre1x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n22) / 1.4142135623731;
    }

    inline double leg_tri_f17_by(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy, ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre1(L1 - L3);
      Legy = Legendre1x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n22) / 1.4142135623731;
    }

    // number 18
    inline double leg_tri_f18_a(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre1(L2 - L1);
      return (k * n31) / 1.0;
    }

    inline double leg_tri_f18_b(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre1(L2 - L1);
      return (k * n32) / 1.0;
    }

    inline double leg_tri_f18_ax(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre1(L2 - L1);
      Legx = Legendre1x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n31) / 1.0;
    }

    inline double leg_tri_f18_ay(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy,  ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre1(L2 - L1);
      Legy = Legendre1x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n31) / 1.0;
    }

    inline double leg_tri_f18_bx(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre1(L2 - L1);
      Legx = Legendre1x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n32) / 1.0;
    }

    inline double leg_tri_f18_by(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy, ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre1(L2 - L1);
      Legy = Legendre1x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n32) / 1.0;
    }

    // Genuine bubble functions, order 3

    // number 19
    inline double leg_tri_f19_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f19_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f19_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f19_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f19_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f19_by(double x, double y)
    {
      return 0.0;
    }

    // number 20
    inline double leg_tri_f20_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f20_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f20_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f20_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f20_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f20_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // ORDER 4

    // Edge functions, order 4

    // number 21
    inline double leg_tri_f21_a0(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 7.0 / 4.0 * Legendre3(L3 - L2); q2 = 3.0 / 4.0 * Legendre2(L3 - L2);
      return (q1 * psi1e1_1(x, y) - q2 * psi0e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f21_b0(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 7.0 / 4.0 * Legendre3(L3 - L2); q2 = 3.0 / 4.0 * Legendre2(L3 - L2);
      return (q1 * psi1e1_2(x, y) - q2 * psi0e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f21_a1(double x, double y)
    {
      return -leg_tri_f21_a0(x, y);
    }

    inline double leg_tri_f21_b1(double x, double y)
    {
      return -leg_tri_f21_b0(x, y);
    }

    inline double leg_tri_f21_ax0(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 7.0 / 4.0 * Legendre3x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 7.0 / 4.0 * Legendre3(L3 - L2);
      cx = 3.0 / 4.0 * Legendre2x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 3.0 / 4.0 * Legendre2(L3 - L2);
      return (ax * psi1e1_1(x, y) + bx * psi1e1x_1(x, y) - cx * psi0e1_1(x, y) - dx * psi0e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f21_ay0(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 7.0 / 4.0 * Legendre3x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 7.0 / 4.0 * Legendre3(L3 - L2);
      cy = 3.0 / 4.0 * Legendre2x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 3.0 / 4.0 * Legendre2(L3 - L2);
      return (ay * psi1e1_1(x, y) + by * psi1e1y_1(x, y) - cy * psi0e1_1(x, y) - dy * psi0e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f21_bx0(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 7.0 / 4.0 * Legendre3x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 7.0 / 4.0 * Legendre3(L3 - L2);
      cx = 3.0 / 4.0 * Legendre2x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 3.0 / 4.0 * Legendre2(L3 - L2);
      return (ax * psi1e1_2(x, y) + bx * psi1e1x_2(x, y) - cx * psi0e1_2(x, y) - dx * psi0e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f21_by0(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 7.0 / 4.0 * Legendre3x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 7.0 / 4.0 * Legendre3(L3 - L2);
      cy = 3.0 / 4.0 * Legendre2x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 3.0 / 4.0 * Legendre2(L3 - L2);
      return (ay * psi1e1_2(x, y) + by * psi1e1y_2(x, y) - cy * psi0e1_2(x, y) - dy * psi0e1y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f21_ax1(double x, double y)
    {
      return -leg_tri_f21_ax0(x, y);
    }

    inline double leg_tri_f21_ay1(double x, double y)
    {
      return -leg_tri_f21_ay0(x, y);
    }

    inline double leg_tri_f21_bx1(double x, double y)
    {
      return -leg_tri_f21_bx0(x, y);
    }

    inline double leg_tri_f21_by1(double x, double y)
    {
      return -leg_tri_f21_by0(x, y);
    }

    // number 22
    inline double leg_tri_f22_a0(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 7.0 / 4.0 * Legendre3(L1 - L3); q2 = 3.0 / 4.0 * Legendre2(L1 - L3);
      return (q1 * psi1e2_1(x, y) - q2 * psi0e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f22_b0(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 7.0 / 4.0 * Legendre3(L1 - L3); q2 = 3.0 / 4.0 * Legendre2(L1 - L3);
      return (q1 * psi1e2_2(x, y) - q2 * psi0e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f22_a1(double x, double y)
    {
      return -leg_tri_f22_a0(x, y);
    }

    inline double leg_tri_f22_b1(double x, double y)
    {
      return -leg_tri_f22_b0(x, y);
    }

    inline double leg_tri_f22_ax0(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 7.0 / 4.0 * Legendre3x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 7.0 / 4.0 * Legendre3(L1 - L3);
      cx = 3.0 / 4.0 * Legendre2x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 3.0 / 4.0 * Legendre2(L1 - L3);
      return (ax * psi1e2_1(x, y) + bx * psi1e2x_1(x, y) - cx * psi0e2_1(x, y) - dx * psi0e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f22_ay0(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 7.0 / 4.0 * Legendre3x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 7.0 / 4.0 * Legendre3(L1 - L3);
      cy = 3.0 / 4.0 * Legendre2x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 3.0 / 4.0 * Legendre2(L1 - L3);
      return (ay * psi1e2_1(x, y) + by * psi1e2y_1(x, y) - cy * psi0e2_1(x, y) - dy * psi0e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f22_bx0(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 7.0 / 4.0 * Legendre3x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 7.0 / 4.0 * Legendre3(L1 - L3);
      cx = 3.0 / 4.0 * Legendre2x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 3.0 / 4.0 * Legendre2(L1 - L3);
      return (ax * psi1e2_2(x, y) + bx * psi1e2x_2(x, y) - cx * psi0e2_2(x, y) - dx * psi0e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f22_by0(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 7.0 / 4.0 * Legendre3x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 7.0 / 4.0 * Legendre3(L1 - L3);
      cy = 3.0 / 4.0 * Legendre2x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 3.0 / 4.0 * Legendre2(L1 - L3);
      return (ay * psi1e2_2(x, y) + by * psi1e2y_2(x, y) - cy * psi0e2_2(x, y) - dy * psi0e2y_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f22_ax1(double x, double y)
    {
      return -leg_tri_f22_ax0(x, y);
    }

    inline double leg_tri_f22_ay1(double x, double y)
    {
      return -leg_tri_f22_ay0(x, y);
    }

    inline double leg_tri_f22_bx1(double x, double y)
    {
      return -leg_tri_f22_bx0(x, y);
    }

    inline double leg_tri_f22_by1(double x, double y)
    {
      return -leg_tri_f22_by0(x, y);
    }

    // number 23
    inline double leg_tri_f23_a0(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 7.0 / 4.0 * Legendre3(L2 - L1); q2 = 3.0 / 4.0 * Legendre2(L2 - L1);
      return (q1 * psi1e3_1(x, y) - q2 * psi0e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f23_b0(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 7.0 / 4.0 * Legendre3(L2 - L1); q2 = 3.0 / 4.0 * Legendre2(L2 - L1);
      return (q1 * psi1e3_2(x, y) - q2 * psi0e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f23_a1(double x, double y)
    {
      return -leg_tri_f23_a0(x, y);
    }

    inline double leg_tri_f23_b1(double x, double y)
    {
      return -leg_tri_f23_b0(x, y);
    }

    inline double leg_tri_f23_ax0(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 7.0 / 4.0 * Legendre3x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 7.0 / 4.0 * Legendre3(L2 - L1);
      cx = 3.0 / 4.0 * Legendre2x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 3.0 / 4.0 * Legendre2(L2 - L1);
      return (ax * psi1e3_1(x, y) + bx * psi1e3x_1(x, y) - cx * psi0e3_1(x, y) - dx * psi0e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f23_ay0(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 7.0 / 4.0 * Legendre3x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 7.0 / 4.0 * Legendre3(L2 - L1);
      cy = 3.0 / 4.0 * Legendre2x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 3.0 / 4.0 * Legendre2(L2 - L1);
      return (ay * psi1e3_1(x, y) + by * psi1e3y_1(x, y) - cy * psi0e3_1(x, y) - dy * psi0e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f23_bx0(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 7.0 / 4.0 * Legendre3x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 7.0 / 4.0 * Legendre3(L2 - L1);
      cx = 3.0 / 4.0 * Legendre2x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 3.0 / 4.0 * Legendre2(L2 - L1);
      return (ax * psi1e3_2(x, y) + bx * psi1e3x_2(x, y) - cx * psi0e3_2(x, y) - dx * psi0e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f23_by0(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 7.0 / 4.0 * Legendre3x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 7.0 / 4.0 * Legendre3(L2 - L1);
      cy = 3.0 / 4.0 * Legendre2x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 3.0 / 4.0 * Legendre2(L2 - L1);
      return (ay * psi1e3_2(x, y) + by * psi1e3y_2(x, y) - cy * psi0e3_2(x, y) - dy * psi0e3y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f23_ax1(double x, double y)
    {
      return -leg_tri_f23_ax0(x, y);
    }

    inline double leg_tri_f23_ay1(double x, double y)
    {
      return -leg_tri_f23_ay0(x, y);
    }

    inline double leg_tri_f23_bx1(double x, double y)
    {
      return -leg_tri_f23_bx0(x, y);
    }

    inline double leg_tri_f23_by1(double x, double y)
    {
      return -leg_tri_f23_by0(x, y);
    }

    // Edge-base bubble functions (normal functions), order 4

    // number 24
    inline double leg_tri_f24_a(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre2(L3 - L2);
      return (k * n11) / 1.0;
    }

    inline double leg_tri_f24_b(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre2(L3 - L2);
      return (k * n12) / 1.0;
    }

    inline double leg_tri_f24_ax(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre2(L3 - L2);
      Legx = Legendre2x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n11) / 1.0;
    }

    inline double leg_tri_f24_ay(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy,  ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre2(L3 - L2);
      Legy = Legendre2x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n11) / 1.0;
    }

    inline double leg_tri_f24_bx(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre2(L3 - L2);
      Legx = Legendre2x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n12) / 1.0;
    }

    inline double leg_tri_f24_by(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy, ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre2(L3 - L2);
      Legy = Legendre2x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n12) / 1.0;
    }

    // number 25
    inline double leg_tri_f25_a(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre2(L1 - L3);
      return (k * n21) / 1.4142135623731;
    }

    inline double leg_tri_f25_b(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre2(L1 - L3);
      return (k * n22) / 1.4142135623731;
    }

    inline double leg_tri_f25_ax(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre2(L1 - L3);
      Legx = Legendre2x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n21) / 1.4142135623731;
    }

    inline double leg_tri_f25_ay(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy,  ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre2(L1 - L3);
      Legy = Legendre2x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n21) / 1.4142135623731;
    }

    inline double leg_tri_f25_bx(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre2(L1 - L3);
      Legx = Legendre2x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n22) / 1.4142135623731;
    }

    inline double leg_tri_f25_by(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy, ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre2(L1 - L3);
      Legy = Legendre2x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n22) / 1.4142135623731;
    }

    // number 26
    inline double leg_tri_f26_a(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre2(L2 - L1);
      return (k * n31) / 1.0;
    }

    inline double leg_tri_f26_b(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre2(L2 - L1);
      return (k * n32) / 1.0;
    }

    inline double leg_tri_f26_ax(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre2(L2 - L1);
      Legx = Legendre2x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n31) / 1.0;
    }

    inline double leg_tri_f26_ay(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy,  ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre2(L2 - L1);
      Legy = Legendre2x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n31) / 1.0;
    }

    inline double leg_tri_f26_bx(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre2(L2 - L1);
      Legx = Legendre2x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n32) / 1.0;
    }

    inline double leg_tri_f26_by(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy, ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre2(L2 - L1);
      Legy = Legendre2x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n32) / 1.0;
    }

    // Genuine bubble functions, order 4

    // number 27
    inline double leg_tri_f27_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f27_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f27_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f27_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f27_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f27_by(double x, double y)
    {
      return 0.0;
    }

    // number 28
    inline double leg_tri_f28_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f28_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f28_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f28_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f28_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f28_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 29
    inline double leg_tri_f29_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f29_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f29_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f29_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f29_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f29_by(double x, double y)
    {
      return 0.0;
    }

    // number 30
    inline double leg_tri_f30_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f30_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f30_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f30_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f30_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f30_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // ORDER 5

    // Edge functions, order 5

    // number 31
    inline double leg_tri_f31_a(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 9.0 / 5.0 * Legendre4(L3 - L2); q2 = 4.0 / 5.0 * Legendre3(L3 - L2);
      return (q1 * psi1e1_1(x, y) - q2 * psi0e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f31_b(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 9.0 / 5.0 * Legendre4(L3 - L2); q2 = 4.0 / 5.0 * Legendre3(L3 - L2);
      return (q1 * psi1e1_2(x, y) - q2 * psi0e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f31_ax(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 9.0 / 5.0 * Legendre4x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 9.0 / 5.0 * Legendre4(L3 - L2);
      cx = 4.0 / 5.0 * Legendre3x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 4.0 / 5.0 * Legendre3(L3 - L2);
      return (ax * psi1e1_1(x, y) + bx * psi1e1x_1(x, y) - cx * psi0e1_1(x, y) - dx * psi0e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f31_ay(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 9.0 / 5.0 * Legendre4x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 9.0 / 5.0 * Legendre4(L3 - L2);
      cy = 4.0 / 5.0 * Legendre3x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 4.0 / 5.0 * Legendre3(L3 - L2);
      return (ay * psi1e1_1(x, y) + by * psi1e1y_1(x, y) - cy * psi0e1_1(x, y) - dy * psi0e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f31_bx(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 9.0 / 5.0 * Legendre4x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 9.0 / 5.0 * Legendre4(L3 - L2);
      cx = 4.0 / 5.0 * Legendre3x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 4.0 / 5.0 * Legendre3(L3 - L2);
      return (ax * psi1e1_2(x, y) + bx * psi1e1x_2(x, y) - cx * psi0e1_2(x, y) - dx * psi0e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f31_by(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 9.0 / 5.0 * Legendre4x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 9.0 / 5.0 * Legendre4(L3 - L2);
      cy = 4.0 / 5.0 * Legendre3x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 4.0 / 5.0 * Legendre3(L3 - L2);
      return (ay * psi1e1_2(x, y) + by * psi1e1y_2(x, y) - cy * psi0e1_2(x, y) - dy * psi0e1y_2(x, y)) / 1.0;
    }

    // number 32
    inline double leg_tri_f32_a(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 9.0 / 5.0 * Legendre4(L1 - L3); q2 = 4.0 / 5.0 * Legendre3(L1 - L3);
      return (q1 * psi1e2_1(x, y) - q2 * psi0e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f32_b(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 9.0 / 5.0 * Legendre4(L1 - L3); q2 = 4.0 / 5.0 * Legendre3(L1 - L3);
      return (q1 * psi1e2_2(x, y) - q2 * psi0e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f32_ax(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 9.0 / 5.0 * Legendre4x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 9.0 / 5.0 * Legendre4(L1 - L3);
      cx = 4.0 / 5.0 * Legendre3x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 4.0 / 5.0 * Legendre3(L1 - L3);
      return (ax * psi1e2_1(x, y) + bx * psi1e2x_1(x, y) - cx * psi0e2_1(x, y) - dx * psi0e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f32_ay(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 9.0 / 5.0 * Legendre4x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 9.0 / 5.0 * Legendre4(L1 - L3);
      cy = 4.0 / 5.0 * Legendre3x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 4.0 / 5.0 * Legendre3(L1 - L3);
      return (ay * psi1e2_1(x, y) + by * psi1e2y_1(x, y) - cy * psi0e2_1(x, y) - dy * psi0e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f32_bx(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 9.0 / 5.0 * Legendre4x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 9.0 / 5.0 * Legendre4(L1 - L3);
      cx = 4.0 / 5.0 * Legendre3x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 4.0 / 5.0 * Legendre3(L1 - L3);
      return (ax * psi1e2_2(x, y) + bx * psi1e2x_2(x, y) - cx * psi0e2_2(x, y) - dx * psi0e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f32_by(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 9.0 / 5.0 * Legendre4x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 9.0 / 5.0 * Legendre4(L1 - L3);
      cy = 4.0 / 5.0 * Legendre3x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 4.0 / 5.0 * Legendre3(L1 - L3);
      return (ay * psi1e2_2(x, y) + by * psi1e2y_2(x, y) - cy * psi0e2_2(x, y) - dy * psi0e2y_2(x, y)) / 1.4142135623731;
    }

    // number 33
    inline double leg_tri_f33_a(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 9.0 / 5.0 * Legendre4(L2 - L1); q2 = 4.0 / 5.0 * Legendre3(L2 - L1);
      return (q1 * psi1e3_1(x, y) - q2 * psi0e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f33_b(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 9.0 / 5.0 * Legendre4(L2 - L1); q2 = 4.0 / 5.0 * Legendre3(L2 - L1);
      return (q1 * psi1e3_2(x, y) - q2 * psi0e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f33_ax(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 9.0 / 5.0 * Legendre4x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 9.0 / 5.0 * Legendre4(L2 - L1);
      cx = 4.0 / 5.0 * Legendre3x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 4.0 / 5.0 * Legendre3(L2 - L1);
      return (ax * psi1e3_1(x, y) + bx * psi1e3x_1(x, y) - cx * psi0e3_1(x, y) - dx * psi0e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f33_ay(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 9.0 / 5.0 * Legendre4x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 9.0 / 5.0 * Legendre4(L2 - L1);
      cy = 4.0 / 5.0 * Legendre3x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 4.0 / 5.0 * Legendre3(L2 - L1);
      return (ay * psi1e3_1(x, y) + by * psi1e3y_1(x, y) - cy * psi0e3_1(x, y) - dy * psi0e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f33_bx(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 9.0 / 5.0 * Legendre4x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 9.0 / 5.0 * Legendre4(L2 - L1);
      cx = 4.0 / 5.0 * Legendre3x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 4.0 / 5.0 * Legendre3(L2 - L1);
      return (ax * psi1e3_2(x, y) + bx * psi1e3x_2(x, y) - cx * psi0e3_2(x, y) - dx * psi0e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f33_by(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 9.0 / 5.0 * Legendre4x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 9.0 / 5.0 * Legendre4(L2 - L1);
      cy = 4.0 / 5.0 * Legendre3x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 4.0 / 5.0 * Legendre3(L2 - L1);
      return (ay * psi1e3_2(x, y) + by * psi1e3y_2(x, y) - cy * psi0e3_2(x, y) - dy * psi0e3y_2(x, y)) / 1.0;
    }

    // Edge-base bubble functions (normal functions), order 5

    // number 34
    inline double leg_tri_f34_a(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre3(L3 - L2);
      return (k * n11) / 1.0;
    }

    inline double leg_tri_f34_b(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre3(L3 - L2);
      return (k * n12) / 1.0;
    }

    inline double leg_tri_f34_ax(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre3(L3 - L2);
      Legx = Legendre3x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n11) / 1.0;
    }

    inline double leg_tri_f34_ay(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy,  ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre3(L3 - L2);
      Legy = Legendre3x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n11) / 1.0;
    }

    inline double leg_tri_f34_bx(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre3(L3 - L2);
      Legx = Legendre3x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n12) / 1.0;
    }

    inline double leg_tri_f34_by(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy, ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre3(L3 - L2);
      Legy = Legendre3x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n12) / 1.0;
    }

    // number 35
    inline double leg_tri_f35_a(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre3(L1 - L3);
      return (k * n21) / 1.4142135623731;
    }

    inline double leg_tri_f35_b(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre3(L1 - L3);
      return (k * n22) / 1.4142135623731;
    }

    inline double leg_tri_f35_ax(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre3(L1 - L3);
      Legx = Legendre3x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n21) / 1.4142135623731;
    }

    inline double leg_tri_f35_ay(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy,  ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre3(L1 - L3);
      Legy = Legendre3x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n21) / 1.4142135623731;
    }

    inline double leg_tri_f35_bx(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre3(L1 - L3);
      Legx = Legendre3x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n22) / 1.4142135623731;
    }

    inline double leg_tri_f35_by(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy, ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre3(L1 - L3);
      Legy = Legendre3x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n22) / 1.4142135623731;
    }

    // number 36
    inline double leg_tri_f36_a(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre3(L2 - L1);
      return (k * n31) / 1.0;
    }

    inline double leg_tri_f36_b(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre3(L2 - L1);
      return (k * n32) / 1.0;
    }

    inline double leg_tri_f36_ax(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre3(L2 - L1);
      Legx = Legendre3x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n31) / 1.0;
    }

    inline double leg_tri_f36_ay(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy,  ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre3(L2 - L1);
      Legy = Legendre3x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n31) / 1.0;
    }

    inline double leg_tri_f36_bx(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre3(L2 - L1);
      Legx = Legendre3x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n32) / 1.0;
    }

    inline double leg_tri_f36_by(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy, ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre3(L2 - L1);
      Legy = Legendre3x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n32) / 1.0;
    }

    // Genuine bubble functions, order 5

    // number 37
    inline double leg_tri_f37_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f37_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f37_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f37_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f37_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f37_by(double x, double y)
    {
      return 0.0;
    }

    // number 38
    inline double leg_tri_f38_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f38_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f38_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f38_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f38_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f38_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 39
    inline double leg_tri_f39_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f39_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f39_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f39_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f39_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f39_by(double x, double y)
    {
      return 0.0;
    }

    // number 40
    inline double leg_tri_f40_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f40_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f40_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f40_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f40_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f40_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 41
    inline double leg_tri_f41_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f41_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f41_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f41_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f41_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f41_by(double x, double y)
    {
      return 0.0;
    }

    // number 42
    inline double leg_tri_f42_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f42_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f42_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f42_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f42_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f42_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // ORDER 6

    // Edge functions, order 6

    // number 43
    inline double leg_tri_f43_a0(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 11.0 / 6.0 * Legendre5(L3 - L2); q2 = 5.0 / 6.0 * Legendre4(L3 - L2);
      return (q1 * psi1e1_1(x, y) - q2 * psi0e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f43_b0(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 11.0 / 6.0 * Legendre5(L3 - L2); q2 = 5.0 / 6.0 * Legendre4(L3 - L2);
      return (q1 * psi1e1_2(x, y) - q2 * psi0e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f43_a1(double x, double y)
    {
      return -leg_tri_f43_a0(x, y);
    }

    inline double leg_tri_f43_b1(double x, double y)
    {
      return -leg_tri_f43_b0(x, y);
    }

    inline double leg_tri_f43_ax0(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 11.0 / 6.0 * Legendre5x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 11.0 / 6.0 * Legendre5(L3 - L2);
      cx = 5.0 / 6.0 * Legendre4x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 5.0 / 6.0 * Legendre4(L3 - L2);
      return (ax * psi1e1_1(x, y) + bx * psi1e1x_1(x, y) - cx * psi0e1_1(x, y) - dx * psi0e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f43_ay0(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 11.0 / 6.0 * Legendre5x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 11.0 / 6.0 * Legendre5(L3 - L2);
      cy = 5.0 / 6.0 * Legendre4x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 5.0 / 6.0 * Legendre4(L3 - L2);
      return (ay * psi1e1_1(x, y) + by * psi1e1y_1(x, y) - cy * psi0e1_1(x, y) - dy * psi0e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f43_bx0(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 11.0 / 6.0 * Legendre5x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 11.0 / 6.0 * Legendre5(L3 - L2);
      cx = 5.0 / 6.0 * Legendre4x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 5.0 / 6.0 * Legendre4(L3 - L2);
      return (ax * psi1e1_2(x, y) + bx * psi1e1x_2(x, y) - cx * psi0e1_2(x, y) - dx * psi0e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f43_by0(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 11.0 / 6.0 * Legendre5x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 11.0 / 6.0 * Legendre5(L3 - L2);
      cy = 5.0 / 6.0 * Legendre4x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 5.0 / 6.0 * Legendre4(L3 - L2);
      return (ay * psi1e1_2(x, y) + by * psi1e1y_2(x, y) - cy * psi0e1_2(x, y) - dy * psi0e1y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f43_ax1(double x, double y)
    {
      return -leg_tri_f43_ax0(x, y);
    }

    inline double leg_tri_f43_ay1(double x, double y)
    {
      return -leg_tri_f43_ay0(x, y);
    }

    inline double leg_tri_f43_bx1(double x, double y)
    {
      return -leg_tri_f43_bx0(x, y);
    }

    inline double leg_tri_f43_by1(double x, double y)
    {
      return -leg_tri_f43_by0(x, y);
    }

    // number 44
    inline double leg_tri_f44_a0(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 11.0 / 6.0 * Legendre5(L1 - L3); q2 = 5.0 / 6.0 * Legendre4(L1 - L3);
      return (q1 * psi1e2_1(x, y) - q2 * psi0e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f44_b0(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 11.0 / 6.0 * Legendre5(L1 - L3); q2 = 5.0 / 6.0 * Legendre4(L1 - L3);
      return (q1 * psi1e2_2(x, y) - q2 * psi0e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f44_a1(double x, double y)
    {
      return -leg_tri_f44_a0(x, y);
    }

    inline double leg_tri_f44_b1(double x, double y)
    {
      return -leg_tri_f44_b0(x, y);
    }

    inline double leg_tri_f44_ax0(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 11.0 / 6.0 * Legendre5x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 11.0 / 6.0 * Legendre5(L1 - L3);
      cx = 5.0 / 6.0 * Legendre4x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 5.0 / 6.0 * Legendre4(L1 - L3);
      return (ax * psi1e2_1(x, y) + bx * psi1e2x_1(x, y) - cx * psi0e2_1(x, y) - dx * psi0e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f44_ay0(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 11.0 / 6.0 * Legendre5x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 11.0 / 6.0 * Legendre5(L1 - L3);
      cy = 5.0 / 6.0 * Legendre4x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 5.0 / 6.0 * Legendre4(L1 - L3);
      return (ay * psi1e2_1(x, y) + by * psi1e2y_1(x, y) - cy * psi0e2_1(x, y) - dy * psi0e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f44_bx0(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 11.0 / 6.0 * Legendre5x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 11.0 / 6.0 * Legendre5(L1 - L3);
      cx = 5.0 / 6.0 * Legendre4x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 5.0 / 6.0 * Legendre4(L1 - L3);
      return (ax * psi1e2_2(x, y) + bx * psi1e2x_2(x, y) - cx * psi0e2_2(x, y) - dx * psi0e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f44_by0(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 11.0 / 6.0 * Legendre5x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 11.0 / 6.0 * Legendre5(L1 - L3);
      cy = 5.0 / 6.0 * Legendre4x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 5.0 / 6.0 * Legendre4(L1 - L3);
      return (ay * psi1e2_2(x, y) + by * psi1e2y_2(x, y) - cy * psi0e2_2(x, y) - dy * psi0e2y_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f44_ax1(double x, double y)
    {
      return -leg_tri_f44_ax0(x, y);
    }

    inline double leg_tri_f44_ay1(double x, double y)
    {
      return -leg_tri_f44_ay0(x, y);
    }

    inline double leg_tri_f44_bx1(double x, double y)
    {
      return -leg_tri_f44_bx0(x, y);
    }

    inline double leg_tri_f44_by1(double x, double y)
    {
      return -leg_tri_f44_by0(x, y);
    }

    // number 45
    inline double leg_tri_f45_a0(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 11.0 / 6.0 * Legendre5(L2 - L1); q2 = 5.0 / 6.0 * Legendre4(L2 - L1);
      return (q1 * psi1e3_1(x, y) - q2 * psi0e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f45_b0(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 11.0 / 6.0 * Legendre5(L2 - L1); q2 = 5.0 / 6.0 * Legendre4(L2 - L1);
      return (q1 * psi1e3_2(x, y) - q2 * psi0e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f45_a1(double x, double y)
    {
      return -leg_tri_f45_a0(x, y);
    }

    inline double leg_tri_f45_b1(double x, double y)
    {
      return -leg_tri_f45_b0(x, y);
    }

    inline double leg_tri_f45_ax0(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 11.0 / 6.0 * Legendre5x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 11.0 / 6.0 * Legendre5(L2 - L1);
      cx = 5.0 / 6.0 * Legendre4x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 5.0 / 6.0 * Legendre4(L2 - L1);
      return (ax * psi1e3_1(x, y) + bx * psi1e3x_1(x, y) - cx * psi0e3_1(x, y) - dx * psi0e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f45_ay0(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 11.0 / 6.0 * Legendre5x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 11.0 / 6.0 * Legendre5(L2 - L1);
      cy = 5.0 / 6.0 * Legendre4x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 5.0 / 6.0 * Legendre4(L2 - L1);
      return (ay * psi1e3_1(x, y) + by * psi1e3y_1(x, y) - cy * psi0e3_1(x, y) - dy * psi0e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f45_bx0(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 11.0 / 6.0 * Legendre5x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 11.0 / 6.0 * Legendre5(L2 - L1);
      cx = 5.0 / 6.0 * Legendre4x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 5.0 / 6.0 * Legendre4(L2 - L1);
      return (ax * psi1e3_2(x, y) + bx * psi1e3x_2(x, y) - cx * psi0e3_2(x, y) - dx * psi0e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f45_by0(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 11.0 / 6.0 * Legendre5x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 11.0 / 6.0 * Legendre5(L2 - L1);
      cy = 5.0 / 6.0 * Legendre4x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 5.0 / 6.0 * Legendre4(L2 - L1);
      return (ay * psi1e3_2(x, y) + by * psi1e3y_2(x, y) - cy * psi0e3_2(x, y) - dy * psi0e3y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f45_ax1(double x, double y)
    {
      return -leg_tri_f45_ax0(x, y);
    }

    inline double leg_tri_f45_ay1(double x, double y)
    {
      return -leg_tri_f45_ay0(x, y);
    }

    inline double leg_tri_f45_bx1(double x, double y)
    {
      return -leg_tri_f45_bx0(x, y);
    }

    inline double leg_tri_f45_by1(double x, double y)
    {
      return -leg_tri_f45_by0(x, y);
    }

    // Edge-base bubble functions (normal functions), order 6

    // number 46
    inline double leg_tri_f46_a(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre4(L3 - L2);
      return (k * n11) / 1.0;
    }

    inline double leg_tri_f46_b(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre4(L3 - L2);
      return (k * n12) / 1.0;
    }

    inline double leg_tri_f46_ax(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre4(L3 - L2);
      Legx = Legendre4x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n11) / 1.0;
    }

    inline double leg_tri_f46_ay(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy,  ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre4(L3 - L2);
      Legy = Legendre4x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n11) / 1.0;
    }

    inline double leg_tri_f46_bx(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre4(L3 - L2);
      Legx = Legendre4x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n12) / 1.0;
    }

    inline double leg_tri_f46_by(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy, ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre4(L3 - L2);
      Legy = Legendre4x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n12) / 1.0;
    }

    // number 47
    inline double leg_tri_f47_a(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre4(L1 - L3);
      return (k * n21) / 1.4142135623731;
    }

    inline double leg_tri_f47_b(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre4(L1 - L3);
      return (k * n22) / 1.4142135623731;
    }

    inline double leg_tri_f47_ax(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre4(L1 - L3);
      Legx = Legendre4x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n21) / 1.4142135623731;
    }

    inline double leg_tri_f47_ay(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy,  ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre4(L1 - L3);
      Legy = Legendre4x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n21) / 1.4142135623731;
    }

    inline double leg_tri_f47_bx(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre4(L1 - L3);
      Legx = Legendre4x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n22) / 1.4142135623731;
    }

    inline double leg_tri_f47_by(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy, ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre4(L1 - L3);
      Legy = Legendre4x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n22) / 1.4142135623731;
    }

    // number 48
    inline double leg_tri_f48_a(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre4(L2 - L1);
      return (k * n31) / 1.0;
    }

    inline double leg_tri_f48_b(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre4(L2 - L1);
      return (k * n32) / 1.0;
    }

    inline double leg_tri_f48_ax(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre4(L2 - L1);
      Legx = Legendre4x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n31) / 1.0;
    }

    inline double leg_tri_f48_ay(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy,  ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre4(L2 - L1);
      Legy = Legendre4x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n31) / 1.0;
    }

    inline double leg_tri_f48_bx(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre4(L2 - L1);
      Legx = Legendre4x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n32) / 1.0;
    }

    inline double leg_tri_f48_by(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy, ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre4(L2 - L1);
      Legy = Legendre4x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n32) / 1.0;
    }

    // Genuine bubble functions, order 6

    // number 49
    inline double leg_tri_f49_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre3(l2 - l1);
    }

    inline double leg_tri_f49_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f49_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre3x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f49_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre3x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f49_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f49_by(double x, double y)
    {
      return 0.0;
    }

    // number 50
    inline double leg_tri_f50_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f50_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre3(l2 - l1);
    }

    inline double leg_tri_f50_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f50_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f50_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre3x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f50_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre3x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 51
    inline double leg_tri_f51_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f51_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f51_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f51_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f51_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f51_by(double x, double y)
    {
      return 0.0;
    }

    // number 52
    inline double leg_tri_f52_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f52_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f52_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f52_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f52_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f52_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 53
    inline double leg_tri_f53_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f53_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f53_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f53_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f53_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f53_by(double x, double y)
    {
      return 0.0;
    }

    // number 54
    inline double leg_tri_f54_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f54_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f54_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f54_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f54_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f54_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 55
    inline double leg_tri_f55_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f55_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f55_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre3x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f55_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre3x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f55_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f55_by(double x, double y)
    {
      return 0.0;
    }

    // number 56
    inline double leg_tri_f56_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f56_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f56_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f56_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f56_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre3x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f56_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre3x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // ORDER 7

    // Edge functions, order 7

    // number 57
    inline double leg_tri_f57_a(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 13.0 / 7.0 * Legendre6(L3 - L2); q2 = 6.0 / 7.0 * Legendre5(L3 - L2);
      return (q1 * psi1e1_1(x, y) - q2 * psi0e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f57_b(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 13.0 / 7.0 * Legendre6(L3 - L2); q2 = 6.0 / 7.0 * Legendre5(L3 - L2);
      return (q1 * psi1e1_2(x, y) - q2 * psi0e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f57_ax(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 13.0 / 7.0 * Legendre6x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 13.0 / 7.0 * Legendre6(L3 - L2);
      cx = 6.0 / 7.0 * Legendre5x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 6.0 / 7.0 * Legendre5(L3 - L2);
      return (ax * psi1e1_1(x, y) + bx * psi1e1x_1(x, y) - cx * psi0e1_1(x, y) - dx * psi0e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f57_ay(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 13.0 / 7.0 * Legendre6x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 13.0 / 7.0 * Legendre6(L3 - L2);
      cy = 6.0 / 7.0 * Legendre5x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 6.0 / 7.0 * Legendre5(L3 - L2);
      return (ay * psi1e1_1(x, y) + by * psi1e1y_1(x, y) - cy * psi0e1_1(x, y) - dy * psi0e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f57_bx(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 13.0 / 7.0 * Legendre6x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 13.0 / 7.0 * Legendre6(L3 - L2);
      cx = 6.0 / 7.0 * Legendre5x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 6.0 / 7.0 * Legendre5(L3 - L2);
      return (ax * psi1e1_2(x, y) + bx * psi1e1x_2(x, y) - cx * psi0e1_2(x, y) - dx * psi0e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f57_by(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 13.0 / 7.0 * Legendre6x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 13.0 / 7.0 * Legendre6(L3 - L2);
      cy = 6.0 / 7.0 * Legendre5x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 6.0 / 7.0 * Legendre5(L3 - L2);
      return (ay * psi1e1_2(x, y) + by * psi1e1y_2(x, y) - cy * psi0e1_2(x, y) - dy * psi0e1y_2(x, y)) / 1.0;
    }

    // number 58
    inline double leg_tri_f58_a(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 13.0 / 7.0 * Legendre6(L1 - L3); q2 = 6.0 / 7.0 * Legendre5(L1 - L3);
      return (q1 * psi1e2_1(x, y) - q2 * psi0e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f58_b(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 13.0 / 7.0 * Legendre6(L1 - L3); q2 = 6.0 / 7.0 * Legendre5(L1 - L3);
      return (q1 * psi1e2_2(x, y) - q2 * psi0e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f58_ax(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 13.0 / 7.0 * Legendre6x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 13.0 / 7.0 * Legendre6(L1 - L3);
      cx = 6.0 / 7.0 * Legendre5x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 6.0 / 7.0 * Legendre5(L1 - L3);
      return (ax * psi1e2_1(x, y) + bx * psi1e2x_1(x, y) - cx * psi0e2_1(x, y) - dx * psi0e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f58_ay(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 13.0 / 7.0 * Legendre6x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 13.0 / 7.0 * Legendre6(L1 - L3);
      cy = 6.0 / 7.0 * Legendre5x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 6.0 / 7.0 * Legendre5(L1 - L3);
      return (ay * psi1e2_1(x, y) + by * psi1e2y_1(x, y) - cy * psi0e2_1(x, y) - dy * psi0e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f58_bx(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 13.0 / 7.0 * Legendre6x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 13.0 / 7.0 * Legendre6(L1 - L3);
      cx = 6.0 / 7.0 * Legendre5x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 6.0 / 7.0 * Legendre5(L1 - L3);
      return (ax * psi1e2_2(x, y) + bx * psi1e2x_2(x, y) - cx * psi0e2_2(x, y) - dx * psi0e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f58_by(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 13.0 / 7.0 * Legendre6x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 13.0 / 7.0 * Legendre6(L1 - L3);
      cy = 6.0 / 7.0 * Legendre5x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 6.0 / 7.0 * Legendre5(L1 - L3);
      return (ay * psi1e2_2(x, y) + by * psi1e2y_2(x, y) - cy * psi0e2_2(x, y) - dy * psi0e2y_2(x, y)) / 1.4142135623731;
    }

    // number 59
    inline double leg_tri_f59_a(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 13.0 / 7.0 * Legendre6(L2 - L1); q2 = 6.0 / 7.0 * Legendre5(L2 - L1);
      return (q1 * psi1e3_1(x, y) - q2 * psi0e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f59_b(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 13.0 / 7.0 * Legendre6(L2 - L1); q2 = 6.0 / 7.0 * Legendre5(L2 - L1);
      return (q1 * psi1e3_2(x, y) - q2 * psi0e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f59_ax(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 13.0 / 7.0 * Legendre6x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 13.0 / 7.0 * Legendre6(L2 - L1);
      cx = 6.0 / 7.0 * Legendre5x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 6.0 / 7.0 * Legendre5(L2 - L1);
      return (ax * psi1e3_1(x, y) + bx * psi1e3x_1(x, y) - cx * psi0e3_1(x, y) - dx * psi0e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f59_ay(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 13.0 / 7.0 * Legendre6x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 13.0 / 7.0 * Legendre6(L2 - L1);
      cy = 6.0 / 7.0 * Legendre5x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 6.0 / 7.0 * Legendre5(L2 - L1);
      return (ay * psi1e3_1(x, y) + by * psi1e3y_1(x, y) - cy * psi0e3_1(x, y) - dy * psi0e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f59_bx(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 13.0 / 7.0 * Legendre6x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 13.0 / 7.0 * Legendre6(L2 - L1);
      cx = 6.0 / 7.0 * Legendre5x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 6.0 / 7.0 * Legendre5(L2 - L1);
      return (ax * psi1e3_2(x, y) + bx * psi1e3x_2(x, y) - cx * psi0e3_2(x, y) - dx * psi0e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f59_by(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 13.0 / 7.0 * Legendre6x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 13.0 / 7.0 * Legendre6(L2 - L1);
      cy = 6.0 / 7.0 * Legendre5x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 6.0 / 7.0 * Legendre5(L2 - L1);
      return (ay * psi1e3_2(x, y) + by * psi1e3y_2(x, y) - cy * psi0e3_2(x, y) - dy * psi0e3y_2(x, y)) / 1.0;
    }

    // Edge-base bubble functions (normal functions), order 7

    // number 60
    inline double leg_tri_f60_a(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre5(L3 - L2);
      return (k * n11) / 1.0;
    }

    inline double leg_tri_f60_b(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre5(L3 - L2);
      return (k * n12) / 1.0;
    }

    inline double leg_tri_f60_ax(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre5(L3 - L2);
      Legx = Legendre5x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n11) / 1.0;
    }

    inline double leg_tri_f60_ay(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy,  ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre5(L3 - L2);
      Legy = Legendre5x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n11) / 1.0;
    }

    inline double leg_tri_f60_bx(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre5(L3 - L2);
      Legx = Legendre5x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n12) / 1.0;
    }

    inline double leg_tri_f60_by(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy, ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre5(L3 - L2);
      Legy = Legendre5x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n12) / 1.0;
    }

    // number 61
    inline double leg_tri_f61_a(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre5(L1 - L3);
      return (k * n21) / 1.4142135623731;
    }

    inline double leg_tri_f61_b(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre5(L1 - L3);
      return (k * n22) / 1.4142135623731;
    }

    inline double leg_tri_f61_ax(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre5(L1 - L3);
      Legx = Legendre5x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n21) / 1.4142135623731;
    }

    inline double leg_tri_f61_ay(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy,  ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre5(L1 - L3);
      Legy = Legendre5x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n21) / 1.4142135623731;
    }

    inline double leg_tri_f61_bx(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre5(L1 - L3);
      Legx = Legendre5x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n22) / 1.4142135623731;
    }

    inline double leg_tri_f61_by(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy, ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre5(L1 - L3);
      Legy = Legendre5x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n22) / 1.4142135623731;
    }

    // number 62
    inline double leg_tri_f62_a(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre5(L2 - L1);
      return (k * n31) / 1.0;
    }

    inline double leg_tri_f62_b(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre5(L2 - L1);
      return (k * n32) / 1.0;
    }

    inline double leg_tri_f62_ax(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre5(L2 - L1);
      Legx = Legendre5x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n31) / 1.0;
    }

    inline double leg_tri_f62_ay(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy,  ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre5(L2 - L1);
      Legy = Legendre5x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n31) / 1.0;
    }

    inline double leg_tri_f62_bx(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre5(L2 - L1);
      Legx = Legendre5x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n32) / 1.0;
    }

    inline double leg_tri_f62_by(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy, ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre5(L2 - L1);
      Legy = Legendre5x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n32) / 1.0;
    }

    // Genuine bubble functions, order 7

    // number 63
    inline double leg_tri_f63_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre4(l2 - l1);
    }

    inline double leg_tri_f63_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f63_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre4x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f63_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre4x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f63_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f63_by(double x, double y)
    {
      return 0.0;
    }

    // number 64
    inline double leg_tri_f64_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f64_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre4(l2 - l1);
    }

    inline double leg_tri_f64_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f64_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f64_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre4x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f64_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre4x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 65
    inline double leg_tri_f65_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre3(l2 - l1);
    }

    inline double leg_tri_f65_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f65_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre3x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f65_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre3x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f65_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f65_by(double x, double y)
    {
      return 0.0;
    }

    // number 66
    inline double leg_tri_f66_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f66_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre3(l2 - l1);
    }

    inline double leg_tri_f66_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f66_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f66_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre3x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f66_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre3x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 67
    inline double leg_tri_f67_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f67_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f67_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f67_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f67_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f67_by(double x, double y)
    {
      return 0.0;
    }

    // number 68
    inline double leg_tri_f68_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f68_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f68_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f68_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f68_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f68_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 69
    inline double leg_tri_f69_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f69_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f69_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre3x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f69_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre3x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f69_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f69_by(double x, double y)
    {
      return 0.0;
    }

    // number 70
    inline double leg_tri_f70_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f70_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f70_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f70_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f70_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre3x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f70_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre3x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 71
    inline double leg_tri_f71_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f71_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f71_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre4x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f71_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre4x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f71_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f71_by(double x, double y)
    {
      return 0.0;
    }

    // number 72
    inline double leg_tri_f72_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f72_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f72_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f72_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f72_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre4x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f72_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre4x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // ORDER 8

    // Edge functions, order 8

    // number 73
    inline double leg_tri_f73_a0(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 15.0 / 8.0 * Legendre7(L3 - L2); q2 = 7.0 / 8.0 * Legendre6(L3 - L2);
      return (q1 * psi1e1_1(x, y) - q2 * psi0e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f73_b0(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 15.0 / 8.0 * Legendre7(L3 - L2); q2 = 7.0 / 8.0 * Legendre6(L3 - L2);
      return (q1 * psi1e1_2(x, y) - q2 * psi0e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f73_a1(double x, double y)
    {
      return -leg_tri_f73_a0(x, y);
    }

    inline double leg_tri_f73_b1(double x, double y)
    {
      return -leg_tri_f73_b0(x, y);
    }

    inline double leg_tri_f73_ax0(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 15.0 / 8.0 * Legendre7x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 15.0 / 8.0 * Legendre7(L3 - L2);
      cx = 7.0 / 8.0 * Legendre6x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 7.0 / 8.0 * Legendre6(L3 - L2);
      return (ax * psi1e1_1(x, y) + bx * psi1e1x_1(x, y) - cx * psi0e1_1(x, y) - dx * psi0e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f73_ay0(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 15.0 / 8.0 * Legendre7x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 15.0 / 8.0 * Legendre7(L3 - L2);
      cy = 7.0 / 8.0 * Legendre6x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 7.0 / 8.0 * Legendre6(L3 - L2);
      return (ay * psi1e1_1(x, y) + by * psi1e1y_1(x, y) - cy * psi0e1_1(x, y) - dy * psi0e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f73_bx0(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 15.0 / 8.0 * Legendre7x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 15.0 / 8.0 * Legendre7(L3 - L2);
      cx = 7.0 / 8.0 * Legendre6x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 7.0 / 8.0 * Legendre6(L3 - L2);
      return (ax * psi1e1_2(x, y) + bx * psi1e1x_2(x, y) - cx * psi0e1_2(x, y) - dx * psi0e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f73_by0(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 15.0 / 8.0 * Legendre7x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 15.0 / 8.0 * Legendre7(L3 - L2);
      cy = 7.0 / 8.0 * Legendre6x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 7.0 / 8.0 * Legendre6(L3 - L2);
      return (ay * psi1e1_2(x, y) + by * psi1e1y_2(x, y) - cy * psi0e1_2(x, y) - dy * psi0e1y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f73_ax1(double x, double y)
    {
      return -leg_tri_f73_ax0(x, y);
    }

    inline double leg_tri_f73_ay1(double x, double y)
    {
      return -leg_tri_f73_ay0(x, y);
    }

    inline double leg_tri_f73_bx1(double x, double y)
    {
      return -leg_tri_f73_bx0(x, y);
    }

    inline double leg_tri_f73_by1(double x, double y)
    {
      return -leg_tri_f73_by0(x, y);
    }

    // number 74
    inline double leg_tri_f74_a0(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 15.0 / 8.0 * Legendre7(L1 - L3); q2 = 7.0 / 8.0 * Legendre6(L1 - L3);
      return (q1 * psi1e2_1(x, y) - q2 * psi0e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f74_b0(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 15.0 / 8.0 * Legendre7(L1 - L3); q2 = 7.0 / 8.0 * Legendre6(L1 - L3);
      return (q1 * psi1e2_2(x, y) - q2 * psi0e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f74_a1(double x, double y)
    {
      return -leg_tri_f74_a0(x, y);
    }

    inline double leg_tri_f74_b1(double x, double y)
    {
      return -leg_tri_f74_b0(x, y);
    }

    inline double leg_tri_f74_ax0(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 15.0 / 8.0 * Legendre7x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 15.0 / 8.0 * Legendre7(L1 - L3);
      cx = 7.0 / 8.0 * Legendre6x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 7.0 / 8.0 * Legendre6(L1 - L3);
      return (ax * psi1e2_1(x, y) + bx * psi1e2x_1(x, y) - cx * psi0e2_1(x, y) - dx * psi0e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f74_ay0(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 15.0 / 8.0 * Legendre7x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 15.0 / 8.0 * Legendre7(L1 - L3);
      cy = 7.0 / 8.0 * Legendre6x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 7.0 / 8.0 * Legendre6(L1 - L3);
      return (ay * psi1e2_1(x, y) + by * psi1e2y_1(x, y) - cy * psi0e2_1(x, y) - dy * psi0e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f74_bx0(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 15.0 / 8.0 * Legendre7x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 15.0 / 8.0 * Legendre7(L1 - L3);
      cx = 7.0 / 8.0 * Legendre6x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 7.0 / 8.0 * Legendre6(L1 - L3);
      return (ax * psi1e2_2(x, y) + bx * psi1e2x_2(x, y) - cx * psi0e2_2(x, y) - dx * psi0e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f74_by0(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 15.0 / 8.0 * Legendre7x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 15.0 / 8.0 * Legendre7(L1 - L3);
      cy = 7.0 / 8.0 * Legendre6x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 7.0 / 8.0 * Legendre6(L1 - L3);
      return (ay * psi1e2_2(x, y) + by * psi1e2y_2(x, y) - cy * psi0e2_2(x, y) - dy * psi0e2y_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f74_ax1(double x, double y)
    {
      return -leg_tri_f74_ax0(x, y);
    }

    inline double leg_tri_f74_ay1(double x, double y)
    {
      return -leg_tri_f74_ay0(x, y);
    }

    inline double leg_tri_f74_bx1(double x, double y)
    {
      return -leg_tri_f74_bx0(x, y);
    }

    inline double leg_tri_f74_by1(double x, double y)
    {
      return -leg_tri_f74_by0(x, y);
    }

    // number 75
    inline double leg_tri_f75_a0(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 15.0 / 8.0 * Legendre7(L2 - L1); q2 = 7.0 / 8.0 * Legendre6(L2 - L1);
      return (q1 * psi1e3_1(x, y) - q2 * psi0e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f75_b0(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 15.0 / 8.0 * Legendre7(L2 - L1); q2 = 7.0 / 8.0 * Legendre6(L2 - L1);
      return (q1 * psi1e3_2(x, y) - q2 * psi0e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f75_a1(double x, double y)
    {
      return -leg_tri_f75_a0(x, y);
    }

    inline double leg_tri_f75_b1(double x, double y)
    {
      return -leg_tri_f75_b0(x, y);
    }

    inline double leg_tri_f75_ax0(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 15.0 / 8.0 * Legendre7x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 15.0 / 8.0 * Legendre7(L2 - L1);
      cx = 7.0 / 8.0 * Legendre6x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 7.0 / 8.0 * Legendre6(L2 - L1);
      return (ax * psi1e3_1(x, y) + bx * psi1e3x_1(x, y) - cx * psi0e3_1(x, y) - dx * psi0e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f75_ay0(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 15.0 / 8.0 * Legendre7x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 15.0 / 8.0 * Legendre7(L2 - L1);
      cy = 7.0 / 8.0 * Legendre6x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 7.0 / 8.0 * Legendre6(L2 - L1);
      return (ay * psi1e3_1(x, y) + by * psi1e3y_1(x, y) - cy * psi0e3_1(x, y) - dy * psi0e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f75_bx0(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 15.0 / 8.0 * Legendre7x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 15.0 / 8.0 * Legendre7(L2 - L1);
      cx = 7.0 / 8.0 * Legendre6x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 7.0 / 8.0 * Legendre6(L2 - L1);
      return (ax * psi1e3_2(x, y) + bx * psi1e3x_2(x, y) - cx * psi0e3_2(x, y) - dx * psi0e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f75_by0(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 15.0 / 8.0 * Legendre7x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 15.0 / 8.0 * Legendre7(L2 - L1);
      cy = 7.0 / 8.0 * Legendre6x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 7.0 / 8.0 * Legendre6(L2 - L1);
      return (ay * psi1e3_2(x, y) + by * psi1e3y_2(x, y) - cy * psi0e3_2(x, y) - dy * psi0e3y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f75_ax1(double x, double y)
    {
      return -leg_tri_f75_ax0(x, y);
    }

    inline double leg_tri_f75_ay1(double x, double y)
    {
      return -leg_tri_f75_ay0(x, y);
    }

    inline double leg_tri_f75_bx1(double x, double y)
    {
      return -leg_tri_f75_bx0(x, y);
    }

    inline double leg_tri_f75_by1(double x, double y)
    {
      return -leg_tri_f75_by0(x, y);
    }

    // Edge-base bubble functions (normal functions), order 8

    // number 76
    inline double leg_tri_f76_a(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre6(L3 - L2);
      return (k * n11) / 1.0;
    }

    inline double leg_tri_f76_b(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre6(L3 - L2);
      return (k * n12) / 1.0;
    }

    inline double leg_tri_f76_ax(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre6(L3 - L2);
      Legx = Legendre6x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n11) / 1.0;
    }

    inline double leg_tri_f76_ay(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy,  ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre6(L3 - L2);
      Legy = Legendre6x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n11) / 1.0;
    }

    inline double leg_tri_f76_bx(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre6(L3 - L2);
      Legx = Legendre6x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n12) / 1.0;
    }

    inline double leg_tri_f76_by(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy, ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre6(L3 - L2);
      Legy = Legendre6x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n12) / 1.0;
    }

    // number 77
    inline double leg_tri_f77_a(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre6(L1 - L3);
      return (k * n21) / 1.4142135623731;
    }

    inline double leg_tri_f77_b(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre6(L1 - L3);
      return (k * n22) / 1.4142135623731;
    }

    inline double leg_tri_f77_ax(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre6(L1 - L3);
      Legx = Legendre6x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n21) / 1.4142135623731;
    }

    inline double leg_tri_f77_ay(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy,  ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre6(L1 - L3);
      Legy = Legendre6x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n21) / 1.4142135623731;
    }

    inline double leg_tri_f77_bx(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre6(L1 - L3);
      Legx = Legendre6x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n22) / 1.4142135623731;
    }

    inline double leg_tri_f77_by(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy, ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre6(L1 - L3);
      Legy = Legendre6x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n22) / 1.4142135623731;
    }

    // number 78
    inline double leg_tri_f78_a(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre6(L2 - L1);
      return (k * n31) / 1.0;
    }

    inline double leg_tri_f78_b(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre6(L2 - L1);
      return (k * n32) / 1.0;
    }

    inline double leg_tri_f78_ax(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre6(L2 - L1);
      Legx = Legendre6x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n31) / 1.0;
    }

    inline double leg_tri_f78_ay(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy,  ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre6(L2 - L1);
      Legy = Legendre6x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n31) / 1.0;
    }

    inline double leg_tri_f78_bx(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre6(L2 - L1);
      Legx = Legendre6x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n32) / 1.0;
    }

    inline double leg_tri_f78_by(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy, ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre6(L2 - L1);
      Legy = Legendre6x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n32) / 1.0;
    }

    // Genuine bubble functions, order 8

    // number 79
    inline double leg_tri_f79_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre5(l2 - l1);
    }

    inline double leg_tri_f79_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f79_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre5x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f79_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre5x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f79_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f79_by(double x, double y)
    {
      return 0.0;
    }

    // number 80
    inline double leg_tri_f80_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f80_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre5(l2 - l1);
    }

    inline double leg_tri_f80_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f80_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f80_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre5x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f80_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre5x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 81
    inline double leg_tri_f81_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre4(l2 - l1);
    }

    inline double leg_tri_f81_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f81_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre4x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f81_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre4x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f81_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f81_by(double x, double y)
    {
      return 0.0;
    }

    // number 82
    inline double leg_tri_f82_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f82_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre4(l2 - l1);
    }

    inline double leg_tri_f82_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f82_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f82_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre4x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f82_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre4x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 83
    inline double leg_tri_f83_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre3(l2 - l1);
    }

    inline double leg_tri_f83_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f83_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre3x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f83_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre3x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f83_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f83_by(double x, double y)
    {
      return 0.0;
    }

    // number 84
    inline double leg_tri_f84_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f84_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre3(l2 - l1);
    }

    inline double leg_tri_f84_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f84_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f84_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre3x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f84_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre3x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 85
    inline double leg_tri_f85_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f85_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f85_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre3x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f85_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre3x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f85_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f85_by(double x, double y)
    {
      return 0.0;
    }

    // number 86
    inline double leg_tri_f86_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f86_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f86_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f86_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f86_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre3x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f86_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre3x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 87
    inline double leg_tri_f87_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f87_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f87_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre4x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f87_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre4x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f87_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f87_by(double x, double y)
    {
      return 0.0;
    }

    // number 88
    inline double leg_tri_f88_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f88_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f88_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f88_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f88_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre4x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f88_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre4x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 89
    inline double leg_tri_f89_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f89_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f89_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre5x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f89_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre5x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f89_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f89_by(double x, double y)
    {
      return 0.0;
    }

    // number 90
    inline double leg_tri_f90_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f90_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f90_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f90_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f90_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre5x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f90_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre5x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // ORDER 9

    // Edge functions, order 9

    // number 91
    inline double leg_tri_f91_a(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 17.0 / 9.0 * Legendre8(L3 - L2); q2 = 8.0 / 9.0 * Legendre7(L3 - L2);
      return (q1 * psi1e1_1(x, y) - q2 * psi0e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f91_b(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 17.0 / 9.0 * Legendre8(L3 - L2); q2 = 8.0 / 9.0 * Legendre7(L3 - L2);
      return (q1 * psi1e1_2(x, y) - q2 * psi0e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f91_ax(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 17.0 / 9.0 * Legendre8x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 17.0 / 9.0 * Legendre8(L3 - L2);
      cx = 8.0 / 9.0 * Legendre7x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 8.0 / 9.0 * Legendre7(L3 - L2);
      return (ax * psi1e1_1(x, y) + bx * psi1e1x_1(x, y) - cx * psi0e1_1(x, y) - dx * psi0e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f91_ay(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 17.0 / 9.0 * Legendre8x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 17.0 / 9.0 * Legendre8(L3 - L2);
      cy = 8.0 / 9.0 * Legendre7x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 8.0 / 9.0 * Legendre7(L3 - L2);
      return (ay * psi1e1_1(x, y) + by * psi1e1y_1(x, y) - cy * psi0e1_1(x, y) - dy * psi0e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f91_bx(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 17.0 / 9.0 * Legendre8x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 17.0 / 9.0 * Legendre8(L3 - L2);
      cx = 8.0 / 9.0 * Legendre7x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 8.0 / 9.0 * Legendre7(L3 - L2);
      return (ax * psi1e1_2(x, y) + bx * psi1e1x_2(x, y) - cx * psi0e1_2(x, y) - dx * psi0e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f91_by(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 17.0 / 9.0 * Legendre8x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 17.0 / 9.0 * Legendre8(L3 - L2);
      cy = 8.0 / 9.0 * Legendre7x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 8.0 / 9.0 * Legendre7(L3 - L2);
      return (ay * psi1e1_2(x, y) + by * psi1e1y_2(x, y) - cy * psi0e1_2(x, y) - dy * psi0e1y_2(x, y)) / 1.0;
    }

    // number 92
    inline double leg_tri_f92_a(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 17.0 / 9.0 * Legendre8(L1 - L3); q2 = 8.0 / 9.0 * Legendre7(L1 - L3);
      return (q1 * psi1e2_1(x, y) - q2 * psi0e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f92_b(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 17.0 / 9.0 * Legendre8(L1 - L3); q2 = 8.0 / 9.0 * Legendre7(L1 - L3);
      return (q1 * psi1e2_2(x, y) - q2 * psi0e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f92_ax(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 17.0 / 9.0 * Legendre8x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 17.0 / 9.0 * Legendre8(L1 - L3);
      cx = 8.0 / 9.0 * Legendre7x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 8.0 / 9.0 * Legendre7(L1 - L3);
      return (ax * psi1e2_1(x, y) + bx * psi1e2x_1(x, y) - cx * psi0e2_1(x, y) - dx * psi0e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f92_ay(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 17.0 / 9.0 * Legendre8x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 17.0 / 9.0 * Legendre8(L1 - L3);
      cy = 8.0 / 9.0 * Legendre7x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 8.0 / 9.0 * Legendre7(L1 - L3);
      return (ay * psi1e2_1(x, y) + by * psi1e2y_1(x, y) - cy * psi0e2_1(x, y) - dy * psi0e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f92_bx(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 17.0 / 9.0 * Legendre8x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 17.0 / 9.0 * Legendre8(L1 - L3);
      cx = 8.0 / 9.0 * Legendre7x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 8.0 / 9.0 * Legendre7(L1 - L3);
      return (ax * psi1e2_2(x, y) + bx * psi1e2x_2(x, y) - cx * psi0e2_2(x, y) - dx * psi0e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f92_by(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 17.0 / 9.0 * Legendre8x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 17.0 / 9.0 * Legendre8(L1 - L3);
      cy = 8.0 / 9.0 * Legendre7x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 8.0 / 9.0 * Legendre7(L1 - L3);
      return (ay * psi1e2_2(x, y) + by * psi1e2y_2(x, y) - cy * psi0e2_2(x, y) - dy * psi0e2y_2(x, y)) / 1.4142135623731;
    }

    // number 93
    inline double leg_tri_f93_a(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 17.0 / 9.0 * Legendre8(L2 - L1); q2 = 8.0 / 9.0 * Legendre7(L2 - L1);
      return (q1 * psi1e3_1(x, y) - q2 * psi0e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f93_b(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 17.0 / 9.0 * Legendre8(L2 - L1); q2 = 8.0 / 9.0 * Legendre7(L2 - L1);
      return (q1 * psi1e3_2(x, y) - q2 * psi0e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f93_ax(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 17.0 / 9.0 * Legendre8x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 17.0 / 9.0 * Legendre8(L2 - L1);
      cx = 8.0 / 9.0 * Legendre7x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 8.0 / 9.0 * Legendre7(L2 - L1);
      return (ax * psi1e3_1(x, y) + bx * psi1e3x_1(x, y) - cx * psi0e3_1(x, y) - dx * psi0e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f93_ay(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 17.0 / 9.0 * Legendre8x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 17.0 / 9.0 * Legendre8(L2 - L1);
      cy = 8.0 / 9.0 * Legendre7x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 8.0 / 9.0 * Legendre7(L2 - L1);
      return (ay * psi1e3_1(x, y) + by * psi1e3y_1(x, y) - cy * psi0e3_1(x, y) - dy * psi0e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f93_bx(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 17.0 / 9.0 * Legendre8x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 17.0 / 9.0 * Legendre8(L2 - L1);
      cx = 8.0 / 9.0 * Legendre7x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 8.0 / 9.0 * Legendre7(L2 - L1);
      return (ax * psi1e3_2(x, y) + bx * psi1e3x_2(x, y) - cx * psi0e3_2(x, y) - dx * psi0e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f93_by(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 17.0 / 9.0 * Legendre8x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 17.0 / 9.0 * Legendre8(L2 - L1);
      cy = 8.0 / 9.0 * Legendre7x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 8.0 / 9.0 * Legendre7(L2 - L1);
      return (ay * psi1e3_2(x, y) + by * psi1e3y_2(x, y) - cy * psi0e3_2(x, y) - dy * psi0e3y_2(x, y)) / 1.0;
    }

    // Edge-base bubble functions (normal functions), order 9

    // number 94
    inline double leg_tri_f94_a(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre7(L3 - L2);
      return (k * n11) / 1.0;
    }

    inline double leg_tri_f94_b(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre7(L3 - L2);
      return (k * n12) / 1.0;
    }

    inline double leg_tri_f94_ax(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre7(L3 - L2);
      Legx = Legendre7x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n11) / 1.0;
    }

    inline double leg_tri_f94_ay(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy,  ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre7(L3 - L2);
      Legy = Legendre7x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n11) / 1.0;
    }

    inline double leg_tri_f94_bx(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre7(L3 - L2);
      Legx = Legendre7x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n12) / 1.0;
    }

    inline double leg_tri_f94_by(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy, ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre7(L3 - L2);
      Legy = Legendre7x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n12) / 1.0;
    }

    // number 95
    inline double leg_tri_f95_a(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre7(L1 - L3);
      return (k * n21) / 1.4142135623731;
    }

    inline double leg_tri_f95_b(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre7(L1 - L3);
      return (k * n22) / 1.4142135623731;
    }

    inline double leg_tri_f95_ax(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre7(L1 - L3);
      Legx = Legendre7x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n21) / 1.4142135623731;
    }

    inline double leg_tri_f95_ay(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy,  ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre7(L1 - L3);
      Legy = Legendre7x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n21) / 1.4142135623731;
    }

    inline double leg_tri_f95_bx(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre7(L1 - L3);
      Legx = Legendre7x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n22) / 1.4142135623731;
    }

    inline double leg_tri_f95_by(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy, ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre7(L1 - L3);
      Legy = Legendre7x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n22) / 1.4142135623731;
    }

    // number 96
    inline double leg_tri_f96_a(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre7(L2 - L1);
      return (k * n31) / 1.0;
    }

    inline double leg_tri_f96_b(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre7(L2 - L1);
      return (k * n32) / 1.0;
    }

    inline double leg_tri_f96_ax(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre7(L2 - L1);
      Legx = Legendre7x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n31) / 1.0;
    }

    inline double leg_tri_f96_ay(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy,  ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre7(L2 - L1);
      Legy = Legendre7x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n31) / 1.0;
    }

    inline double leg_tri_f96_bx(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre7(L2 - L1);
      Legx = Legendre7x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n32) / 1.0;
    }

    inline double leg_tri_f96_by(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy, ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre7(L2 - L1);
      Legy = Legendre7x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n32) / 1.0;
    }

    // Genuine bubble functions, order 9

    // number 97
    inline double leg_tri_f97_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre6(l2 - l1);
    }

    inline double leg_tri_f97_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f97_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre6(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre6x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f97_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre6(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre6x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f97_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f97_by(double x, double y)
    {
      return 0.0;
    }

    // number 98
    inline double leg_tri_f98_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f98_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre6(l2 - l1);
    }

    inline double leg_tri_f98_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f98_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f98_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre6(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre6x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f98_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre6(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre6x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 99
    inline double leg_tri_f99_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre5(l2 - l1);
    }

    inline double leg_tri_f99_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f99_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre5x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f99_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre5x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f99_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f99_by(double x, double y)
    {
      return 0.0;
    }

    // number 100
    inline double leg_tri_f100_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f100_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre5(l2 - l1);
    }

    inline double leg_tri_f100_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f100_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f100_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre5x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f100_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre5x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 101
    inline double leg_tri_f101_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre4(l2 - l1);
    }

    inline double leg_tri_f101_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f101_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre4x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f101_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre4x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f101_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f101_by(double x, double y)
    {
      return 0.0;
    }

    // number 102
    inline double leg_tri_f102_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f102_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre4(l2 - l1);
    }

    inline double leg_tri_f102_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f102_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f102_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre4x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f102_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre4x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 103
    inline double leg_tri_f103_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre3(l2 - l1);
    }

    inline double leg_tri_f103_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f103_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1x = Legendre3x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre3x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f103_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1y = Legendre3x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre3x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f103_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f103_by(double x, double y)
    {
      return 0.0;
    }

    // number 104
    inline double leg_tri_f104_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f104_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre3(l2 - l1);
    }

    inline double leg_tri_f104_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f104_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f104_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1x = Legendre3x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre3x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f104_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1y = Legendre3x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre3x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 105
    inline double leg_tri_f105_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f105_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f105_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre4x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f105_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre4x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f105_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f105_by(double x, double y)
    {
      return 0.0;
    }

    // number 106
    inline double leg_tri_f106_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f106_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f106_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f106_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f106_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre4x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f106_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre4x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 107
    inline double leg_tri_f107_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f107_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f107_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre5x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f107_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre5x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f107_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f107_by(double x, double y)
    {
      return 0.0;
    }

    // number 108
    inline double leg_tri_f108_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f108_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f108_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f108_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f108_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre5x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f108_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre5x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 109
    inline double leg_tri_f109_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre6(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f109_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f109_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre6(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre6x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f109_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre6(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre6x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f109_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f109_by(double x, double y)
    {
      return 0.0;
    }

    // number 110
    inline double leg_tri_f110_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f110_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre6(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f110_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f110_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f110_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre6(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre6x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f110_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre6(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre6x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // ORDER 10

    // Edge functions, order 10

    // number 111
    inline double leg_tri_f111_a0(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 19.0 / 10.0 * Legendre9(L3 - L2); q2 = 9.0 / 10.0 * Legendre8(L3 - L2);
      return (q1 * psi1e1_1(x, y) - q2 * psi0e1_1(x, y)) / 1.0;
    }

    inline double leg_tri_f111_b0(double x, double y)
    {
      double L3, L2, q1, q2;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      q1 = 19.0 / 10.0 * Legendre9(L3 - L2); q2 = 9.0 / 10.0 * Legendre8(L3 - L2);
      return (q1 * psi1e1_2(x, y) - q2 * psi0e1_2(x, y)) / 1.0;
    }

    inline double leg_tri_f111_a1(double x, double y)
    {
      return -leg_tri_f111_a0(x, y);
    }

    inline double leg_tri_f111_b1(double x, double y)
    {
      return -leg_tri_f111_b0(x, y);
    }

    inline double leg_tri_f111_ax0(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 19.0 / 10.0 * Legendre9x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 19.0 / 10.0 * Legendre9(L3 - L2);
      cx = 9.0 / 10.0 * Legendre8x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 9.0 / 10.0 * Legendre8(L3 - L2);
      return (ax * psi1e1_1(x, y) + bx * psi1e1x_1(x, y) - cx * psi0e1_1(x, y) - dx * psi0e1x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f111_ay0(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 19.0 / 10.0 * Legendre9x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 19.0 / 10.0 * Legendre9(L3 - L2);
      cy = 9.0 / 10.0 * Legendre8x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 9.0 / 10.0 * Legendre8(L3 - L2);
      return (ay * psi1e1_1(x, y) + by * psi1e1y_1(x, y) - cy * psi0e1_1(x, y) - dy * psi0e1y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f111_bx0(double x, double y)
    {
      double L3, L2, ax, bx, cx, dx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ax = 19.0 / 10.0 * Legendre9x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      bx = 19.0 / 10.0 * Legendre9(L3 - L2);
      cx = 9.0 / 10.0 * Legendre8x(L3 - L2) * (lambda3x(x, y) - lambda2x(x, y));
      dx = 9.0 / 10.0 * Legendre8(L3 - L2);
      return (ax * psi1e1_2(x, y) + bx * psi1e1x_2(x, y) - cx * psi0e1_2(x, y) - dx * psi0e1x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f111_by0(double x, double y)
    {
      double L3, L2, ay, by, cy, dy;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      ay = 19.0 / 10.0 * Legendre9x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      by = 19.0 / 10.0 * Legendre9(L3 - L2);
      cy = 9.0 / 10.0 * Legendre8x(L3 - L2) * (lambda3y(x, y) - lambda2y(x, y));
      dy = 9.0 / 10.0 * Legendre8(L3 - L2);
      return (ay * psi1e1_2(x, y) + by * psi1e1y_2(x, y) - cy * psi0e1_2(x, y) - dy * psi0e1y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f111_ax1(double x, double y)
    {
      return -leg_tri_f111_ax0(x, y);
    }

    inline double leg_tri_f111_ay1(double x, double y)
    {
      return -leg_tri_f111_ay0(x, y);
    }

    inline double leg_tri_f111_bx1(double x, double y)
    {
      return -leg_tri_f111_bx0(x, y);
    }

    inline double leg_tri_f111_by1(double x, double y)
    {
      return -leg_tri_f111_by0(x, y);
    }

    // number 112
    inline double leg_tri_f112_a0(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 19.0 / 10.0 * Legendre9(L1 - L3); q2 = 9.0 / 10.0 * Legendre8(L1 - L3);
      return (q1 * psi1e2_1(x, y) - q2 * psi0e2_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f112_b0(double x, double y)
    {
      double L1, L3, q1, q2;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      q1 = 19.0 / 10.0 * Legendre9(L1 - L3); q2 = 9.0 / 10.0 * Legendre8(L1 - L3);
      return (q1 * psi1e2_2(x, y) - q2 * psi0e2_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f112_a1(double x, double y)
    {
      return -leg_tri_f112_a0(x, y);
    }

    inline double leg_tri_f112_b1(double x, double y)
    {
      return -leg_tri_f112_b0(x, y);
    }

    inline double leg_tri_f112_ax0(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 19.0 / 10.0 * Legendre9x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 19.0 / 10.0 * Legendre9(L1 - L3);
      cx = 9.0 / 10.0 * Legendre8x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 9.0 / 10.0 * Legendre8(L1 - L3);
      return (ax * psi1e2_1(x, y) + bx * psi1e2x_1(x, y) - cx * psi0e2_1(x, y) - dx * psi0e2x_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f112_ay0(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 19.0 / 10.0 * Legendre9x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 19.0 / 10.0 * Legendre9(L1 - L3);
      cy = 9.0 / 10.0 * Legendre8x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 9.0 / 10.0 * Legendre8(L1 - L3);
      return (ay * psi1e2_1(x, y) + by * psi1e2y_1(x, y) - cy * psi0e2_1(x, y) - dy * psi0e2y_1(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f112_bx0(double x, double y)
    {
      double L1, L3, ax, bx, cx, dx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ax = 19.0 / 10.0 * Legendre9x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      bx = 19.0 / 10.0 * Legendre9(L1 - L3);
      cx = 9.0 / 10.0 * Legendre8x(L1 - L3) * (lambda1x(x, y) - lambda3x(x, y));
      dx = 9.0 / 10.0 * Legendre8(L1 - L3);
      return (ax * psi1e2_2(x, y) + bx * psi1e2x_2(x, y) - cx * psi0e2_2(x, y) - dx * psi0e2x_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f112_by0(double x, double y)
    {
      double L1, L3, ay, by, cy, dy;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      ay = 19.0 / 10.0 * Legendre9x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      by = 19.0 / 10.0 * Legendre9(L1 - L3);
      cy = 9.0 / 10.0 * Legendre8x(L1 - L3) * (lambda1y(x, y) - lambda3y(x, y));
      dy = 9.0 / 10.0 * Legendre8(L1 - L3);
      return (ay * psi1e2_2(x, y) + by * psi1e2y_2(x, y) - cy * psi0e2_2(x, y) - dy * psi0e2y_2(x, y)) / 1.4142135623731;
    }

    inline double leg_tri_f112_ax1(double x, double y)
    {
      return -leg_tri_f112_ax0(x, y);
    }

    inline double leg_tri_f112_ay1(double x, double y)
    {
      return -leg_tri_f112_ay0(x, y);
    }

    inline double leg_tri_f112_bx1(double x, double y)
    {
      return -leg_tri_f112_bx0(x, y);
    }

    inline double leg_tri_f112_by1(double x, double y)
    {
      return -leg_tri_f112_by0(x, y);
    }

    // number 113
    inline double leg_tri_f113_a0(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 19.0 / 10.0 * Legendre9(L2 - L1); q2 = 9.0 / 10.0 * Legendre8(L2 - L1);
      return (q1 * psi1e3_1(x, y) - q2 * psi0e3_1(x, y)) / 1.0;
    }

    inline double leg_tri_f113_b0(double x, double y)
    {
      double L2, L1, q1, q2;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      q1 = 19.0 / 10.0 * Legendre9(L2 - L1); q2 = 9.0 / 10.0 * Legendre8(L2 - L1);
      return (q1 * psi1e3_2(x, y) - q2 * psi0e3_2(x, y)) / 1.0;
    }

    inline double leg_tri_f113_a1(double x, double y)
    {
      return -leg_tri_f113_a0(x, y);
    }

    inline double leg_tri_f113_b1(double x, double y)
    {
      return -leg_tri_f113_b0(x, y);
    }

    inline double leg_tri_f113_ax0(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 19.0 / 10.0 * Legendre9x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 19.0 / 10.0 * Legendre9(L2 - L1);
      cx = 9.0 / 10.0 * Legendre8x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 9.0 / 10.0 * Legendre8(L2 - L1);
      return (ax * psi1e3_1(x, y) + bx * psi1e3x_1(x, y) - cx * psi0e3_1(x, y) - dx * psi0e3x_1(x, y)) / 1.0;
    }

    inline double leg_tri_f113_ay0(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 19.0 / 10.0 * Legendre9x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 19.0 / 10.0 * Legendre9(L2 - L1);
      cy = 9.0 / 10.0 * Legendre8x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 9.0 / 10.0 * Legendre8(L2 - L1);
      return (ay * psi1e3_1(x, y) + by * psi1e3y_1(x, y) - cy * psi0e3_1(x, y) - dy * psi0e3y_1(x, y)) / 1.0;
    }

    inline double leg_tri_f113_bx0(double x, double y)
    {
      double L2, L1, ax, bx, cx, dx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ax = 19.0 / 10.0 * Legendre9x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      bx = 19.0 / 10.0 * Legendre9(L2 - L1);
      cx = 9.0 / 10.0 * Legendre8x(L2 - L1) * (lambda2x(x, y) - lambda1x(x, y));
      dx = 9.0 / 10.0 * Legendre8(L2 - L1);
      return (ax * psi1e3_2(x, y) + bx * psi1e3x_2(x, y) - cx * psi0e3_2(x, y) - dx * psi0e3x_2(x, y)) / 1.0;
    }

    inline double leg_tri_f113_by0(double x, double y)
    {
      double L2, L1, ay, by, cy, dy;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      ay = 19.0 / 10.0 * Legendre9x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      by = 19.0 / 10.0 * Legendre9(L2 - L1);
      cy = 9.0 / 10.0 * Legendre8x(L2 - L1) * (lambda2y(x, y) - lambda1y(x, y));
      dy = 9.0 / 10.0 * Legendre8(L2 - L1);
      return (ay * psi1e3_2(x, y) + by * psi1e3y_2(x, y) - cy * psi0e3_2(x, y) - dy * psi0e3y_2(x, y)) / 1.0;
    }

    inline double leg_tri_f113_ax1(double x, double y)
    {
      return -leg_tri_f113_ax0(x, y);
    }

    inline double leg_tri_f113_ay1(double x, double y)
    {
      return -leg_tri_f113_ay0(x, y);
    }

    inline double leg_tri_f113_bx1(double x, double y)
    {
      return -leg_tri_f113_bx0(x, y);
    }

    inline double leg_tri_f113_by1(double x, double y)
    {
      return -leg_tri_f113_by0(x, y);
    }

    // Edge-base bubble functions (normal functions), order 10

    // number 114
    inline double leg_tri_f114_a(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre8(L3 - L2);
      return (k * n11) / 1.0;
    }

    inline double leg_tri_f114_b(double x, double y)
    {
      double L3, L2, k;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      k = L3 * L2 * Legendre8(L3 - L2);
      return (k * n12) / 1.0;
    }

    inline double leg_tri_f114_ax(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre8(L3 - L2);
      Legx = Legendre8x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n11) / 1.0;
    }

    inline double leg_tri_f114_ay(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy,  ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre8(L3 - L2);
      Legy = Legendre8x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n11) / 1.0;
    }

    inline double leg_tri_f114_bx(double x, double y)
    {
      double L3, L2, L3x, L2x,
      Leg, Legx, kx;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3x = lambda3x(x, y); L2x = lambda2x(x, y);
      Leg = Legendre8(L3 - L2);
      Legx = Legendre8x(L3 - L2) * (L3x - L2x);
      kx = L3x * L2 * Leg + L3 * L2x * Leg + L3 * L2 * Legx;
      return (kx * n12) / 1.0;
    }

    inline double leg_tri_f114_by(double x, double y)
    {
      double L3, L2,
      L3y, L2y, Leg, Legy, ky;
      L3 = lambda3(x, y); L2 = lambda2(x, y);
      L3y = lambda3y(x, y); L2y = lambda2y(x, y);
      Leg = Legendre8(L3 - L2);
      Legy = Legendre8x(L3 - L2) * (L3y - L2y);
      ky = L3y * L2 * Leg + L3 * L2y * Leg + L3 * L2 * Legy;
      return  (ky * n12) / 1.0;
    }

    // number 115
    inline double leg_tri_f115_a(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre8(L1 - L3);
      return (k * n21) / 1.4142135623731;
    }

    inline double leg_tri_f115_b(double x, double y)
    {
      double L1, L3, k;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      k = L1 * L3 * Legendre8(L1 - L3);
      return (k * n22) / 1.4142135623731;
    }

    inline double leg_tri_f115_ax(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre8(L1 - L3);
      Legx = Legendre8x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n21) / 1.4142135623731;
    }

    inline double leg_tri_f115_ay(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy,  ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre8(L1 - L3);
      Legy = Legendre8x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n21) / 1.4142135623731;
    }

    inline double leg_tri_f115_bx(double x, double y)
    {
      double L1, L3, L1x, L3x,
      Leg, Legx, kx;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1x = lambda1x(x, y); L3x = lambda3x(x, y);
      Leg = Legendre8(L1 - L3);
      Legx = Legendre8x(L1 - L3) * (L1x - L3x);
      kx = L1x * L3 * Leg + L1 * L3x * Leg + L1 * L3 * Legx;
      return (kx * n22) / 1.4142135623731;
    }

    inline double leg_tri_f115_by(double x, double y)
    {
      double L1, L3,
      L1y, L3y, Leg, Legy, ky;
      L1 = lambda1(x, y); L3 = lambda3(x, y);
      L1y = lambda1y(x, y); L3y = lambda3y(x, y);
      Leg = Legendre8(L1 - L3);
      Legy = Legendre8x(L1 - L3) * (L1y - L3y);
      ky = L1y * L3 * Leg + L1 * L3y * Leg + L1 * L3 * Legy;
      return  (ky * n22) / 1.4142135623731;
    }

    // number 116
    inline double leg_tri_f116_a(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre8(L2 - L1);
      return (k * n31) / 1.0;
    }

    inline double leg_tri_f116_b(double x, double y)
    {
      double L2, L1, k;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      k = L2 * L1 * Legendre8(L2 - L1);
      return (k * n32) / 1.0;
    }

    inline double leg_tri_f116_ax(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre8(L2 - L1);
      Legx = Legendre8x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n31) / 1.0;
    }

    inline double leg_tri_f116_ay(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy,  ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre8(L2 - L1);
      Legy = Legendre8x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n31) / 1.0;
    }

    inline double leg_tri_f116_bx(double x, double y)
    {
      double L2, L1, L2x, L1x,
      Leg, Legx, kx;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2x = lambda2x(x, y); L1x = lambda1x(x, y);
      Leg = Legendre8(L2 - L1);
      Legx = Legendre8x(L2 - L1) * (L2x - L1x);
      kx = L2x * L1 * Leg + L2 * L1x * Leg + L2 * L1 * Legx;
      return (kx * n32) / 1.0;
    }

    inline double leg_tri_f116_by(double x, double y)
    {
      double L2, L1,
      L2y, L1y, Leg, Legy, ky;
      L2 = lambda2(x, y); L1 = lambda1(x, y);
      L2y = lambda2y(x, y); L1y = lambda1y(x, y);
      Leg = Legendre8(L2 - L1);
      Legy = Legendre8x(L2 - L1) * (L2y - L1y);
      ky = L2y * L1 * Leg + L2 * L1y * Leg + L2 * L1 * Legy;
      return  (ky * n32) / 1.0;
    }

    // Genuine bubble functions, order 10

    // number 117
    inline double leg_tri_f117_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre7(l2 - l1);
    }

    inline double leg_tri_f117_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f117_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre7(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre7x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f117_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre7(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre7x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f117_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f117_by(double x, double y)
    {
      return 0.0;
    }

    // number 118
    inline double leg_tri_f118_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f118_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre0(l3 - l2) * Legendre7(l2 - l1);
    }

    inline double leg_tri_f118_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f118_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f118_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre7(l2 - l1);
      Leg1x = Legendre0x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre7x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f118_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre0(l3 - l2); Leg2 = Legendre7(l2 - l1);
      Leg1y = Legendre0x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre7x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 119
    inline double leg_tri_f119_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre6(l2 - l1);
    }

    inline double leg_tri_f119_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f119_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre6(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre6x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f119_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre6(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre6x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f119_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f119_by(double x, double y)
    {
      return 0.0;
    }

    // number 120
    inline double leg_tri_f120_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f120_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre1(l3 - l2) * Legendre6(l2 - l1);
    }

    inline double leg_tri_f120_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f120_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f120_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre6(l2 - l1);
      Leg1x = Legendre1x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre6x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f120_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre1(l3 - l2); Leg2 = Legendre6(l2 - l1);
      Leg1y = Legendre1x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre6x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 121
    inline double leg_tri_f121_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre5(l2 - l1);
    }

    inline double leg_tri_f121_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f121_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre5x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f121_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre5x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f121_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f121_by(double x, double y)
    {
      return 0.0;
    }

    // number 122
    inline double leg_tri_f122_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f122_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre2(l3 - l2) * Legendre5(l2 - l1);
    }

    inline double leg_tri_f122_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f122_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f122_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1x = Legendre2x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre5x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f122_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre2(l3 - l2); Leg2 = Legendre5(l2 - l1);
      Leg1y = Legendre2x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre5x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 123
    inline double leg_tri_f123_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre4(l2 - l1);
    }

    inline double leg_tri_f123_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f123_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1x = Legendre3x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre4x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f123_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1y = Legendre3x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre4x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f123_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f123_by(double x, double y)
    {
      return 0.0;
    }

    // number 124
    inline double leg_tri_f124_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f124_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre3(l3 - l2) * Legendre4(l2 - l1);
    }

    inline double leg_tri_f124_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f124_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f124_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1x = Legendre3x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre4x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f124_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre3(l3 - l2); Leg2 = Legendre4(l2 - l1);
      Leg1y = Legendre3x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre4x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 125
    inline double leg_tri_f125_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre3(l2 - l1);
    }

    inline double leg_tri_f125_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f125_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1x = Legendre4x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre3x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f125_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1y = Legendre4x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre3x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f125_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f125_by(double x, double y)
    {
      return 0.0;
    }

    // number 126
    inline double leg_tri_f126_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f126_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre4(l3 - l2) * Legendre3(l2 - l1);
    }

    inline double leg_tri_f126_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f126_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f126_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1x = Legendre4x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre3x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f126_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre4(l3 - l2); Leg2 = Legendre3(l2 - l1);
      Leg1y = Legendre4x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre3x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 127
    inline double leg_tri_f127_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f127_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f127_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre5x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f127_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre5x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f127_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f127_by(double x, double y)
    {
      return 0.0;
    }

    // number 128
    inline double leg_tri_f128_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f128_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre5(l3 - l2) * Legendre2(l2 - l1);
    }

    inline double leg_tri_f128_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f128_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f128_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1x = Legendre5x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre2x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f128_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre5(l3 - l2); Leg2 = Legendre2(l2 - l1);
      Leg1y = Legendre5x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre2x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 129
    inline double leg_tri_f129_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre6(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f129_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f129_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre6(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre6x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f129_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre6(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre6x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f129_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f129_by(double x, double y)
    {
      return 0.0;
    }

    // number 130
    inline double leg_tri_f130_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f130_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre6(l3 - l2) * Legendre1(l2 - l1);
    }

    inline double leg_tri_f130_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f130_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f130_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre6(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1x = Legendre6x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre1x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f130_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre6(l3 - l2); Leg2 = Legendre1(l2 - l1);
      Leg1y = Legendre6x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre1x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    // number 131
    inline double leg_tri_f131_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre7(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f131_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f131_ax(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre7(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre7x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f131_ay(double x, double y)
    {
      double l1, l2, l3, l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre7(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre7x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    inline double leg_tri_f131_bx(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f131_by(double x, double y)
    {
      return 0.0;
    }

    // number 132
    inline double leg_tri_f132_a(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return 0.0;
    }

    inline double leg_tri_f132_b(double x, double y)
    {
      double l1, l2, l3;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      return l1 * l2 * l3 * Legendre7(l3 - l2) * Legendre0(l2 - l1);
    }

    inline double leg_tri_f132_ax(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f132_ay(double x, double y)
    {
      return 0.0;
    }

    inline double leg_tri_f132_bx(double x, double y)
    {
      double l1, l2, l3, l1x, l2x, l3x, Leg1, Leg2, Leg1x, Leg2x;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1x = lambda1x(x, y); l2x = lambda2x(x, y); l3x = lambda3x(x, y);
      Leg1 = Legendre7(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1x = Legendre7x(l3 - l2) * (l3x - l2x);
      Leg2x = Legendre0x(l2 - l1) * (l2x - l1x);
      return l1x * l2 * l3 * Leg1 * Leg2 + l1 * l2x * l3 * Leg1 * Leg2 + l1 * l2 * l3x * Leg1 * Leg2 + l1 * l2 * l3 * Leg1x * Leg2 + l1 * l2 * l3 * Leg1 * Leg2x;
    }

    inline double leg_tri_f132_by(double x, double y)
    {
      double l1, l2, l3,  l1y, l2y, l3y, Leg1, Leg2, Leg1y, Leg2y;
      l1 = lambda1(x, y); l2 = lambda2(x, y); l3 = lambda3(x, y);
      l1y = lambda1y(x, y); l2y = lambda2y(x, y); l3y = lambda3y(x, y);
      Leg1 = Legendre7(l3 - l2); Leg2 = Legendre0(l2 - l1);
      Leg1y = Legendre7x(l3 - l2) * (l3y - l2y);
      Leg2y = Legendre0x(l2 - l1) * (l2y - l1y);
      return l1y * l2 * l3 * Leg1 * Leg2 + l1 * l2y * l3 * Leg1 * Leg2 + l1 * l2 * l3y * Leg1 * Leg2 + l1 * l2 * l3 * Leg1y * Leg2 + l1 * l2 * l3 * Leg1 * Leg2y;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////

    static Shapeset::shape_fn_t leg_tri_fn_a[] =
    {
      leg_tri_f1_a0,  leg_tri_f1_a1,  leg_tri_f2_a0,  leg_tri_f2_a1,  leg_tri_f3_a0,  leg_tri_f3_a1,
      leg_tri_f4_a,  leg_tri_f5_a,   leg_tri_f6_a,   leg_tri_f7_a0,  leg_tri_f7_a1,  leg_tri_f8_a0,
      leg_tri_f8_a1,  leg_tri_f9_a0, leg_tri_f9_a1,  leg_tri_f10_a,  leg_tri_f11_a,  leg_tri_f12_a,
      leg_tri_f13_a,  leg_tri_f14_a,  leg_tri_f15_a, leg_tri_f16_a,  leg_tri_f17_a,  leg_tri_f18_a,
      leg_tri_f19_a,  leg_tri_f20_a,  leg_tri_f21_a0, leg_tri_f21_a1, leg_tri_f22_a0, leg_tri_f22_a1,
      leg_tri_f23_a0, leg_tri_f23_a1, leg_tri_f24_a,  leg_tri_f25_a,  leg_tri_f26_a,  leg_tri_f27_a,
      leg_tri_f28_a,  leg_tri_f29_a,  leg_tri_f30_a,  leg_tri_f31_a,  leg_tri_f32_a,  leg_tri_f33_a,
      leg_tri_f34_a,  leg_tri_f35_a,  leg_tri_f36_a,  leg_tri_f37_a,  leg_tri_f38_a,  leg_tri_f39_a,
      leg_tri_f40_a,  leg_tri_f41_a,  leg_tri_f42_a,  leg_tri_f43_a0, leg_tri_f43_a1, leg_tri_f44_a0,
      leg_tri_f44_a1, leg_tri_f45_a0, leg_tri_f45_a1, leg_tri_f46_a,  leg_tri_f47_a,  leg_tri_f48_a,
      leg_tri_f49_a,  leg_tri_f50_a,  leg_tri_f51_a,  leg_tri_f52_a,  leg_tri_f53_a,  leg_tri_f54_a,
      leg_tri_f55_a,  leg_tri_f56_a,  leg_tri_f57_a,  leg_tri_f58_a, leg_tri_f59_a,  leg_tri_f60_a,
      leg_tri_f61_a,  leg_tri_f62_a,  leg_tri_f63_a,  leg_tri_f64_a,  leg_tri_f65_a,  leg_tri_f66_a,
      leg_tri_f67_a,  leg_tri_f68_a,  leg_tri_f69_a,  leg_tri_f70_a,  leg_tri_f71_a,  leg_tri_f72_a,
      leg_tri_f73_a0, leg_tri_f73_a1, leg_tri_f74_a0, leg_tri_f74_a1, leg_tri_f75_a0, leg_tri_f75_a1,
      leg_tri_f76_a,  leg_tri_f77_a,  leg_tri_f78_a,  leg_tri_f79_a,  leg_tri_f80_a,  leg_tri_f81_a,
      leg_tri_f82_a,  leg_tri_f83_a,  leg_tri_f84_a,  leg_tri_f85_a,  leg_tri_f86_a,  leg_tri_f87_a,
      leg_tri_f88_a,  leg_tri_f89_a,  leg_tri_f90_a,  leg_tri_f91_a,  leg_tri_f92_a,  leg_tri_f93_a,
      leg_tri_f94_a,  leg_tri_f95_a,  leg_tri_f96_a,  leg_tri_f97_a,  leg_tri_f98_a,  leg_tri_f99_a,
      leg_tri_f100_a,  leg_tri_f101_a,  leg_tri_f102_a,  leg_tri_f103_a,  leg_tri_f104_a,  leg_tri_f105_a,
      leg_tri_f106_a,  leg_tri_f107_a,  leg_tri_f108_a,  leg_tri_f109_a,  leg_tri_f110_a,  leg_tri_f111_a0,
      leg_tri_f111_a1, leg_tri_f112_a0, leg_tri_f112_a1, leg_tri_f113_a0, leg_tri_f113_a1, leg_tri_f114_a,
      leg_tri_f115_a,  leg_tri_f116_a,  leg_tri_f117_a,  leg_tri_f118_a,  leg_tri_f119_a,  leg_tri_f120_a,
      leg_tri_f121_a,  leg_tri_f122_a,  leg_tri_f123_a,  leg_tri_f124_a,  leg_tri_f125_a,  leg_tri_f126_a,
      leg_tri_f127_a,  leg_tri_f128_a,  leg_tri_f129_a,  leg_tri_f130_a,  leg_tri_f131_a,  leg_tri_f132_a
    };

    static Shapeset::shape_fn_t leg_tri_fn_b[] =
    {
      leg_tri_f1_b0,  leg_tri_f1_b1,  leg_tri_f2_b0,  leg_tri_f2_b1,  leg_tri_f3_b0,
      leg_tri_f3_b1,  leg_tri_f4_b,  leg_tri_f5_b,   leg_tri_f6_b,   leg_tri_f7_b0,
      leg_tri_f7_b1,  leg_tri_f8_b0,  leg_tri_f8_b1,  leg_tri_f9_b0, leg_tri_f9_b1,
      leg_tri_f10_b,  leg_tri_f11_b,  leg_tri_f12_b,  leg_tri_f13_b,  leg_tri_f14_b,
      leg_tri_f15_b, leg_tri_f16_b,  leg_tri_f17_b,  leg_tri_f18_b,  leg_tri_f19_b,
      leg_tri_f20_b,  leg_tri_f21_b0, leg_tri_f21_b1, leg_tri_f22_b0, leg_tri_f22_b1,
      leg_tri_f23_b0, leg_tri_f23_b1, leg_tri_f24_b,  leg_tri_f25_b,  leg_tri_f26_b,
      leg_tri_f27_b,  leg_tri_f28_b,  leg_tri_f29_b,  leg_tri_f30_b,  leg_tri_f31_b,
      leg_tri_f32_b,  leg_tri_f33_b, leg_tri_f34_b,  leg_tri_f35_b,  leg_tri_f36_b,
      leg_tri_f37_b,  leg_tri_f38_b,  leg_tri_f39_b,  leg_tri_f40_b,  leg_tri_f41_b,
      leg_tri_f42_b,  leg_tri_f43_b0, leg_tri_f43_b1, leg_tri_f44_b0, leg_tri_f44_b1,
      leg_tri_f45_b0, leg_tri_f45_b1, leg_tri_f46_b,  leg_tri_f47_b,  leg_tri_f48_b,
      leg_tri_f49_b,  leg_tri_f50_b,  leg_tri_f51_b,  leg_tri_f52_b,  leg_tri_f53_b,
      leg_tri_f54_b,  leg_tri_f55_b,  leg_tri_f56_b,  leg_tri_f57_b,  leg_tri_f58_b,
      leg_tri_f59_b,  leg_tri_f60_b,  leg_tri_f61_b,  leg_tri_f62_b,  leg_tri_f63_b,
      leg_tri_f64_b,  leg_tri_f65_b,  leg_tri_f66_b,  leg_tri_f67_b,  leg_tri_f68_b,
      leg_tri_f69_b,  leg_tri_f70_b,  leg_tri_f71_b,  leg_tri_f72_b,  leg_tri_f73_b0,
      leg_tri_f73_b1, leg_tri_f74_b0, leg_tri_f74_b1, leg_tri_f75_b0, leg_tri_f75_b1,
      leg_tri_f76_b,  leg_tri_f77_b,  leg_tri_f78_b,  leg_tri_f79_b,  leg_tri_f80_b,
      leg_tri_f81_b,  leg_tri_f82_b,  leg_tri_f83_b,  leg_tri_f84_b,  leg_tri_f85_b,
      leg_tri_f86_b,  leg_tri_f87_b,  leg_tri_f88_b,  leg_tri_f89_b,  leg_tri_f90_b,
      leg_tri_f91_b,  leg_tri_f92_b,  leg_tri_f93_b,  leg_tri_f94_b,  leg_tri_f95_b,
      leg_tri_f96_b,  leg_tri_f97_b,  leg_tri_f98_b,  leg_tri_f99_b,  leg_tri_f100_b,
      leg_tri_f101_b,  leg_tri_f102_b,  leg_tri_f103_b,  leg_tri_f104_b,  leg_tri_f105_b,
      leg_tri_f106_b,  leg_tri_f107_b,  leg_tri_f108_b,  leg_tri_f109_b,  leg_tri_f110_b,
      leg_tri_f111_b0, leg_tri_f111_b1, leg_tri_f112_b0, leg_tri_f112_b1, leg_tri_f113_b0,
      leg_tri_f113_b1, leg_tri_f114_b,  leg_tri_f115_b,  leg_tri_f116_b,  leg_tri_f117_b,
      leg_tri_f118_b,  leg_tri_f119_b,  leg_tri_f120_b,  leg_tri_f121_b,  leg_tri_f122_b,
      leg_tri_f123_b,  leg_tri_f124_b,  leg_tri_f125_b,  leg_tri_f126_b,  leg_tri_f127_b,
      leg_tri_f128_b,  leg_tri_f129_b,  leg_tri_f130_b,  leg_tri_f131_b,  leg_tri_f132_b
    };

    static Shapeset::shape_fn_t leg_tri_fn_ax[] =
    {
      leg_tri_f1_ax0,  leg_tri_f1_ax1,  leg_tri_f2_ax0,  leg_tri_f2_ax1,  leg_tri_f3_ax0,  leg_tri_f3_ax1,
      leg_tri_f4_ax,  leg_tri_f5_ax,   leg_tri_f6_ax,   leg_tri_f7_ax0,  leg_tri_f7_ax1,  leg_tri_f8_ax0,
      leg_tri_f8_ax1,  leg_tri_f9_ax0, leg_tri_f9_ax1,  leg_tri_f10_ax,  leg_tri_f11_ax,  leg_tri_f12_ax,
      leg_tri_f13_ax,  leg_tri_f14_ax,  leg_tri_f15_ax, leg_tri_f16_ax,  leg_tri_f17_ax,  leg_tri_f18_ax,
      leg_tri_f19_ax,  leg_tri_f20_ax,  leg_tri_f21_ax0, leg_tri_f21_ax1, leg_tri_f22_ax0, leg_tri_f22_ax1,
      leg_tri_f23_ax0, leg_tri_f23_ax1, leg_tri_f24_ax,  leg_tri_f25_ax,  leg_tri_f26_ax,  leg_tri_f27_ax,
      leg_tri_f28_ax,  leg_tri_f29_ax,  leg_tri_f30_ax,  leg_tri_f31_ax,  leg_tri_f32_ax,  leg_tri_f33_ax,
      leg_tri_f34_ax,  leg_tri_f35_ax,  leg_tri_f36_ax,  leg_tri_f37_ax,  leg_tri_f38_ax,  leg_tri_f39_ax,
      leg_tri_f40_ax,  leg_tri_f41_ax,  leg_tri_f42_ax,  leg_tri_f43_ax0, leg_tri_f43_ax1, leg_tri_f44_ax0,
      leg_tri_f44_ax1, leg_tri_f45_ax0, leg_tri_f45_ax1, leg_tri_f46_ax,  leg_tri_f47_ax,  leg_tri_f48_ax,
      leg_tri_f49_ax,  leg_tri_f50_ax,  leg_tri_f51_ax,  leg_tri_f52_ax,  leg_tri_f53_ax,  leg_tri_f54_ax,
      leg_tri_f55_ax,  leg_tri_f56_ax,  leg_tri_f57_ax,  leg_tri_f58_ax, leg_tri_f59_ax,  leg_tri_f60_ax,
      leg_tri_f61_ax,  leg_tri_f62_ax,  leg_tri_f63_ax,  leg_tri_f64_ax,  leg_tri_f65_ax,  leg_tri_f66_ax,
      leg_tri_f67_ax,  leg_tri_f68_ax,  leg_tri_f69_ax,  leg_tri_f70_ax,  leg_tri_f71_ax,  leg_tri_f72_ax,
      leg_tri_f73_ax0, leg_tri_f73_ax1, leg_tri_f74_ax0, leg_tri_f74_ax1, leg_tri_f75_ax0, leg_tri_f75_ax1,
      leg_tri_f76_ax,  leg_tri_f77_ax,  leg_tri_f78_ax,  leg_tri_f79_ax,  leg_tri_f80_ax,  leg_tri_f81_ax,
      leg_tri_f82_ax,  leg_tri_f83_ax,  leg_tri_f84_ax,  leg_tri_f85_ax,  leg_tri_f86_ax,  leg_tri_f87_ax,
      leg_tri_f88_ax,  leg_tri_f89_ax,  leg_tri_f90_ax,  leg_tri_f91_ax,  leg_tri_f92_ax,  leg_tri_f93_ax,
      leg_tri_f94_ax,  leg_tri_f95_ax,  leg_tri_f96_ax,  leg_tri_f97_ax,  leg_tri_f98_ax,  leg_tri_f99_ax,
      leg_tri_f100_ax,  leg_tri_f101_ax,  leg_tri_f102_ax,  leg_tri_f103_ax,  leg_tri_f104_ax,  leg_tri_f105_ax,
      leg_tri_f106_ax,  leg_tri_f107_ax,  leg_tri_f108_ax,  leg_tri_f109_ax,  leg_tri_f110_ax,  leg_tri_f111_ax0,
      leg_tri_f111_ax1, leg_tri_f112_ax0, leg_tri_f112_ax1, leg_tri_f113_ax0, leg_tri_f113_ax1, leg_tri_f114_ax,
      leg_tri_f115_ax,  leg_tri_f116_ax,  leg_tri_f117_ax,  leg_tri_f118_ax,  leg_tri_f119_ax,  leg_tri_f120_ax,
      leg_tri_f121_ax,  leg_tri_f122_ax,  leg_tri_f123_ax,  leg_tri_f124_ax,  leg_tri_f125_ax,  leg_tri_f126_ax,
      leg_tri_f127_ax,  leg_tri_f128_ax,  leg_tri_f129_ax,  leg_tri_f130_ax,  leg_tri_f131_ax,  leg_tri_f132_ax
    };

    static Shapeset::shape_fn_t leg_tri_fn_ay[] =
    {
      leg_tri_f1_ay0,  leg_tri_f1_ay1,  leg_tri_f2_ay0,  leg_tri_f2_ay1,  leg_tri_f3_ay0,  leg_tri_f3_ay1,
      leg_tri_f4_ay,  leg_tri_f5_ay,   leg_tri_f6_ay,   leg_tri_f7_ay0,  leg_tri_f7_ay1,  leg_tri_f8_ay0,
      leg_tri_f8_ay1,  leg_tri_f9_ay0, leg_tri_f9_ay1,  leg_tri_f10_ay,  leg_tri_f11_ay,  leg_tri_f12_ay,
      leg_tri_f13_ay,  leg_tri_f14_ay,  leg_tri_f15_ay, leg_tri_f16_ay,  leg_tri_f17_ay,  leg_tri_f18_ay,
      leg_tri_f19_ay,  leg_tri_f20_ay,  leg_tri_f21_ay0, leg_tri_f21_ay1, leg_tri_f22_ay0, leg_tri_f22_ay1,
      leg_tri_f23_ay0, leg_tri_f23_ay1, leg_tri_f24_ay,  leg_tri_f25_ay,  leg_tri_f26_ay,  leg_tri_f27_ay,
      leg_tri_f28_ay,  leg_tri_f29_ay,  leg_tri_f30_ay,  leg_tri_f31_ay,  leg_tri_f32_ay,  leg_tri_f33_ay,
      leg_tri_f34_ay,  leg_tri_f35_ay,  leg_tri_f36_ay,  leg_tri_f37_ay,  leg_tri_f38_ay,  leg_tri_f39_ay,
      leg_tri_f40_ay,  leg_tri_f41_ay,  leg_tri_f42_ay,  leg_tri_f43_ay0, leg_tri_f43_ay1, leg_tri_f44_ay0,
      leg_tri_f44_ay1, leg_tri_f45_ay0, leg_tri_f45_ay1, leg_tri_f46_ay,  leg_tri_f47_ay,  leg_tri_f48_ay,
      leg_tri_f49_ay,  leg_tri_f50_ay,  leg_tri_f51_ay,  leg_tri_f52_ay,  leg_tri_f53_ay,  leg_tri_f54_ay,
      leg_tri_f55_ay,  leg_tri_f56_ay,  leg_tri_f57_ay,  leg_tri_f58_ay, leg_tri_f59_ay,  leg_tri_f60_ay,
      leg_tri_f61_ay,  leg_tri_f62_ay,  leg_tri_f63_ay,  leg_tri_f64_ay,  leg_tri_f65_ay,  leg_tri_f66_ay,
      leg_tri_f67_ay,  leg_tri_f68_ay,  leg_tri_f69_ay,  leg_tri_f70_ay,  leg_tri_f71_ay,  leg_tri_f72_ay,
      leg_tri_f73_ay0, leg_tri_f73_ay1, leg_tri_f74_ay0, leg_tri_f74_ay1, leg_tri_f75_ay0, leg_tri_f75_ay1,
      leg_tri_f76_ay,  leg_tri_f77_ay,  leg_tri_f78_ay,  leg_tri_f79_ay,  leg_tri_f80_ay,  leg_tri_f81_ay,
      leg_tri_f82_ay,  leg_tri_f83_ay,  leg_tri_f84_ay,  leg_tri_f85_ay,  leg_tri_f86_ay,  leg_tri_f87_ay,
      leg_tri_f88_ay,  leg_tri_f89_ay,  leg_tri_f90_ay,  leg_tri_f91_ay,  leg_tri_f92_ay,  leg_tri_f93_ay,
      leg_tri_f94_ay,  leg_tri_f95_ay,  leg_tri_f96_ay,  leg_tri_f97_ay,  leg_tri_f98_ay,  leg_tri_f99_ay,
      leg_tri_f100_ay,  leg_tri_f101_ay,  leg_tri_f102_ay,  leg_tri_f103_ay,  leg_tri_f104_ay,  leg_tri_f105_ay,
      leg_tri_f106_ay,  leg_tri_f107_ay,  leg_tri_f108_ay,  leg_tri_f109_ay,  leg_tri_f110_ay,  leg_tri_f111_ay0,
      leg_tri_f111_ay1, leg_tri_f112_ay0, leg_tri_f112_ay1, leg_tri_f113_ay0, leg_tri_f113_ay1, leg_tri_f114_ay,
      leg_tri_f115_ay,  leg_tri_f116_ay,  leg_tri_f117_ay,  leg_tri_f118_ay,  leg_tri_f119_ay,  leg_tri_f120_ay,
      leg_tri_f121_ay,  leg_tri_f122_ay,  leg_tri_f123_ay,  leg_tri_f124_ay,  leg_tri_f125_ay,  leg_tri_f126_ay,
      leg_tri_f127_ay,  leg_tri_f128_ay,  leg_tri_f129_ay,  leg_tri_f130_ay,  leg_tri_f131_ay,  leg_tri_f132_ay
    };

    static Shapeset::shape_fn_t leg_tri_fn_bx[] =
    {
    leg_tri_f1_bx0,  leg_tri_f1_bx1,  leg_tri_f2_bx0,  leg_tri_f2_bx1,  leg_tri_f3_bx0,
      leg_tri_f3_bx1,  leg_tri_f4_bx,  leg_tri_f5_bx,   leg_tri_f6_bx,   leg_tri_f7_bx0,
      leg_tri_f7_bx1,  leg_tri_f8_bx0,  leg_tri_f8_bx1,  leg_tri_f9_bx0, leg_tri_f9_bx1,
      leg_tri_f10_bx,  leg_tri_f11_bx,  leg_tri_f12_bx,  leg_tri_f13_bx,  leg_tri_f14_bx,
      leg_tri_f15_bx, leg_tri_f16_bx,  leg_tri_f17_bx,  leg_tri_f18_bx,  leg_tri_f19_bx,
      leg_tri_f20_bx,  leg_tri_f21_bx0, leg_tri_f21_bx1, leg_tri_f22_bx0, leg_tri_f22_bx1,
     leg_tri_f23_bx0, leg_tri_f23_bx1, leg_tri_f24_bx,  leg_tri_f25_bx,  leg_tri_f26_bx,
      leg_tri_f27_bx,  leg_tri_f28_bx,  leg_tri_f29_bx,  leg_tri_f30_bx,  leg_tri_f31_bx,
      leg_tri_f32_bx,  leg_tri_f33_bx, leg_tri_f34_bx,  leg_tri_f35_bx,  leg_tri_f36_bx,
      leg_tri_f37_bx,  leg_tri_f38_bx,  leg_tri_f39_bx,  leg_tri_f40_bx,  leg_tri_f41_bx,
      leg_tri_f42_bx,  leg_tri_f43_bx0, leg_tri_f43_bx1, leg_tri_f44_bx0, leg_tri_f44_bx1,
     leg_tri_f45_bx0, leg_tri_f45_bx1, leg_tri_f46_bx,  leg_tri_f47_bx,  leg_tri_f48_bx,
      leg_tri_f49_bx,  leg_tri_f50_bx,  leg_tri_f51_bx,  leg_tri_f52_bx,  leg_tri_f53_bx,
      leg_tri_f54_bx,  leg_tri_f55_bx,  leg_tri_f56_bx,  leg_tri_f57_bx,  leg_tri_f58_bx,
     leg_tri_f59_bx,  leg_tri_f60_bx,  leg_tri_f61_bx,  leg_tri_f62_bx,  leg_tri_f63_bx,
      leg_tri_f64_bx,  leg_tri_f65_bx,  leg_tri_f66_bx,  leg_tri_f67_bx,  leg_tri_f68_bx,
      leg_tri_f69_bx,  leg_tri_f70_bx,  leg_tri_f71_bx,  leg_tri_f72_bx,  leg_tri_f73_bx0,
     leg_tri_f73_bx1, leg_tri_f74_bx0, leg_tri_f74_bx1, leg_tri_f75_bx0, leg_tri_f75_bx1,
     leg_tri_f76_bx,  leg_tri_f77_bx,  leg_tri_f78_bx,  leg_tri_f79_bx,  leg_tri_f80_bx,
      leg_tri_f81_bx,  leg_tri_f82_bx,  leg_tri_f83_bx,  leg_tri_f84_bx,  leg_tri_f85_bx,
      leg_tri_f86_bx,  leg_tri_f87_bx,  leg_tri_f88_bx,  leg_tri_f89_bx,  leg_tri_f90_bx,
      leg_tri_f91_bx,  leg_tri_f92_bx,  leg_tri_f93_bx,  leg_tri_f94_bx,  leg_tri_f95_bx,
      leg_tri_f96_bx,  leg_tri_f97_bx,  leg_tri_f98_bx,  leg_tri_f99_bx,  leg_tri_f100_bx,
      leg_tri_f101_bx,  leg_tri_f102_bx,  leg_tri_f103_bx,  leg_tri_f104_bx,  leg_tri_f105_bx,
      leg_tri_f106_bx,  leg_tri_f107_bx,  leg_tri_f108_bx,  leg_tri_f109_bx,  leg_tri_f110_bx,
      leg_tri_f111_bx0, leg_tri_f111_bx1, leg_tri_f112_bx0, leg_tri_f112_bx1, leg_tri_f113_bx0,
     leg_tri_f113_bx1, leg_tri_f114_bx,  leg_tri_f115_bx,  leg_tri_f116_bx,  leg_tri_f117_bx,
      leg_tri_f118_bx,  leg_tri_f119_bx,  leg_tri_f120_bx,  leg_tri_f121_bx,  leg_tri_f122_bx,
      leg_tri_f123_bx,  leg_tri_f124_bx,  leg_tri_f125_bx,  leg_tri_f126_bx,  leg_tri_f127_bx,
      leg_tri_f128_bx,  leg_tri_f129_bx,  leg_tri_f130_bx,  leg_tri_f131_bx,  leg_tri_f132_bx
    };

    static Shapeset::shape_fn_t leg_tri_fn_by[] =
    {
    leg_tri_f1_by0,  leg_tri_f1_by1,  leg_tri_f2_by0,  leg_tri_f2_by1,  leg_tri_f3_by0,
      leg_tri_f3_by1,  leg_tri_f4_by,  leg_tri_f5_by,   leg_tri_f6_by,   leg_tri_f7_by0,
      leg_tri_f7_by1,  leg_tri_f8_by0,  leg_tri_f8_by1,  leg_tri_f9_by0, leg_tri_f9_by1,
      leg_tri_f10_by,  leg_tri_f11_by,  leg_tri_f12_by,  leg_tri_f13_by,  leg_tri_f14_by,
      leg_tri_f15_by, leg_tri_f16_by,  leg_tri_f17_by,  leg_tri_f18_by,  leg_tri_f19_by,
      leg_tri_f20_by,  leg_tri_f21_by0, leg_tri_f21_by1, leg_tri_f22_by0, leg_tri_f22_by1,
     leg_tri_f23_by0, leg_tri_f23_by1, leg_tri_f24_by,  leg_tri_f25_by,  leg_tri_f26_by,
      leg_tri_f27_by,  leg_tri_f28_by,  leg_tri_f29_by,  leg_tri_f30_by,  leg_tri_f31_by,
      leg_tri_f32_by,  leg_tri_f33_by, leg_tri_f34_by,  leg_tri_f35_by,  leg_tri_f36_by,
      leg_tri_f37_by,  leg_tri_f38_by,  leg_tri_f39_by,  leg_tri_f40_by,  leg_tri_f41_by,
      leg_tri_f42_by,  leg_tri_f43_by0, leg_tri_f43_by1, leg_tri_f44_by0, leg_tri_f44_by1,
     leg_tri_f45_by0, leg_tri_f45_by1, leg_tri_f46_by,  leg_tri_f47_by,  leg_tri_f48_by,
      leg_tri_f49_by,  leg_tri_f50_by,  leg_tri_f51_by,  leg_tri_f52_by,  leg_tri_f53_by,
      leg_tri_f54_by,  leg_tri_f55_by,  leg_tri_f56_by,  leg_tri_f57_by,  leg_tri_f58_by,
     leg_tri_f59_by,  leg_tri_f60_by,  leg_tri_f61_by,  leg_tri_f62_by,  leg_tri_f63_by,
      leg_tri_f64_by,  leg_tri_f65_by,  leg_tri_f66_by,  leg_tri_f67_by,  leg_tri_f68_by,
      leg_tri_f69_by,  leg_tri_f70_by,  leg_tri_f71_by,  leg_tri_f72_by,  leg_tri_f73_by0,
     leg_tri_f73_by1, leg_tri_f74_by0, leg_tri_f74_by1, leg_tri_f75_by0, leg_tri_f75_by1,
     leg_tri_f76_by,  leg_tri_f77_by,  leg_tri_f78_by,  leg_tri_f79_by,  leg_tri_f80_by,
      leg_tri_f81_by,  leg_tri_f82_by,  leg_tri_f83_by,  leg_tri_f84_by,  leg_tri_f85_by,
      leg_tri_f86_by,  leg_tri_f87_by,  leg_tri_f88_by,  leg_tri_f89_by,  leg_tri_f90_by,
      leg_tri_f91_by,  leg_tri_f92_by,  leg_tri_f93_by,  leg_tri_f94_by,  leg_tri_f95_by,
      leg_tri_f96_by,  leg_tri_f97_by,  leg_tri_f98_by,  leg_tri_f99_by,  leg_tri_f100_by,
      leg_tri_f101_by,  leg_tri_f102_by,  leg_tri_f103_by,  leg_tri_f104_by,  leg_tri_f105_by,
      leg_tri_f106_by,  leg_tri_f107_by,  leg_tri_f108_by,  leg_tri_f109_by,  leg_tri_f110_by,
      leg_tri_f111_by0, leg_tri_f111_by1, leg_tri_f112_by0, leg_tri_f112_by1, leg_tri_f113_by0,
     leg_tri_f113_by1, leg_tri_f114_by,  leg_tri_f115_by,  leg_tri_f116_by,  leg_tri_f117_by,
      leg_tri_f118_by,  leg_tri_f119_by,  leg_tri_f120_by,  leg_tri_f121_by,  leg_tri_f122_by,
      leg_tri_f123_by,  leg_tri_f124_by,  leg_tri_f125_by,  leg_tri_f126_by,  leg_tri_f127_by,
      leg_tri_f128_by,  leg_tri_f129_by,  leg_tri_f130_by,  leg_tri_f131_by,  leg_tri_f132_by
    };

    static int leg_tri_bubble_indices_all_orders[]=
    {
      15, 16, 17,
      21, 22, 23, 24, 25,
      32, 33, 34, 35, 36, 37, 38 ,
      42, 43, 44, 45, 46, 47, 48, 49, 50,
      57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
      71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104,
      108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124 ,
      131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149
    };

    static int* leg_tri_bubble_indices[11] =
    {
      NULL, NULL,
      leg_tri_bubble_indices_all_orders,
      leg_tri_bubble_indices_all_orders,
      leg_tri_bubble_indices_all_orders,
      leg_tri_bubble_indices_all_orders,
      leg_tri_bubble_indices_all_orders,
      leg_tri_bubble_indices_all_orders,
      leg_tri_bubble_indices_all_orders,
      leg_tri_bubble_indices_all_orders,
      leg_tri_bubble_indices_all_orders
    };

    static int leg_tri_bubble_count[11] = {0, 0, 3, 8, 15, 24, 35, 48, 63, 80, 99};

    static int leg_tri_edge_indices_0[22] =  { 0, 1, 6, 6,  9, 10, 18, 18, 26, 27, 39, 39, 51, 52, 68, 68, 84, 85, 105, 105 , 125, 126};
    static int leg_tri_edge_indices_1[22] =  { 2, 3, 7, 7, 11, 12, 19, 19, 28, 29, 40, 40, 53, 54, 69, 69, 86, 87, 106, 106 , 127, 128};
    static int leg_tri_edge_indices_2[22] =  { 4, 5, 8, 8, 13, 14, 20, 20, 30, 31, 41, 41, 55, 56, 70, 70, 88, 89, 107, 107 , 129, 130};

    static int* leg_tri_edge_indices[3] =
    {
      leg_tri_edge_indices_0,
      leg_tri_edge_indices_1,
      leg_tri_edge_indices_2
    };

    static int leg_tri_vertex_indices[3] = { -1, -1, -1 };

    static int leg_tri_index_to_order[150] =
    {
      /*6*/  1, 1, 1, 1, 1, 1,
      /*3*/  1, 1, 1,
      /*9*/  2, 2, 2, 2, 2, 2, 2, 2, 2,
      /*8*/  3, 3, 3, 3, 3, 3, 3, 3,
      /*13*/ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      /*12*/ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
      /*17*/ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
      /*16*/ 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
      /*21*/ 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
      /*20*/ 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
      /*25*/ 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
    };

    static Shapeset::shape_fn_t* leg_tri_shape_fn_table[2] =
    {
      leg_tri_fn_a,
      leg_tri_fn_b
    };

    static Shapeset::shape_fn_t* leg_tri_shape_fn_table_x[2] =
    {
      leg_tri_fn_ax,
      leg_tri_fn_bx
    };

    static Shapeset::shape_fn_t* leg_tri_shape_fn_table_y[2] =
    {
      leg_tri_fn_ay,
      leg_tri_fn_by
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // QUADY

    // Shapeset for the curl operator for quads, based on Legendre polynomials

    /* EDGE FUNCTIONS */

    /* ORDER 0 */

    static double leg_quad_p0_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e1_b_1(double x, double y)
    {
      return l0(x) * Legendre0(y);
    }

    static double leg_quad_p0_e1_b_0(double x, double y)
    {
      return -(l0(x) * Legendre0(y));
    }

    static double leg_quad_p0_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e1_bx_1(double x, double y)
    {
      return dl0(x) * Legendre0(y);
    }

    static double leg_quad_p0_e1_bx_0(double x, double y)
    {
      return -(dl0(x) * Legendre0(y));
    }

    static double leg_quad_p0_e1_by_1(double x, double y)
    {
      return l0(x) * Legendre0x(y);
    }

    static double leg_quad_p0_e1_by_0(double x, double y)
    {
      return -(l0(x) * Legendre0x(y));
    }

    static double leg_quad_p0_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e2_b_0(double x, double y)
    {
      return l1(x) * Legendre0(y);
    }

    static double leg_quad_p0_e2_b_1(double x, double y)
    {
      return -(l1(x) * Legendre0(y));
    }

    static double leg_quad_p0_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e2_bx_0(double x, double y)
    {
      return dl1(x) * Legendre0(y);
    }

    static double leg_quad_p0_e2_bx_1(double x, double y)
    {
      return -(dl1(x) * Legendre0(y));
    }

    static double leg_quad_p0_e2_by_0(double x, double y)
    {
      return l1(x) * Legendre0x(y);
    }

    static double leg_quad_p0_e2_by_1(double x, double y)
    {
      return -(l1(x) * Legendre0x(y));
    }

    static double leg_quad_p0_e3_a_0(double x, double y)
    {
      return Legendre0(x) * l0(y);
    }

    static double leg_quad_p0_e3_a_1(double x, double y)
    {
      return -(Legendre0(x) * l0(y));
    }

    static double leg_quad_p0_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e3_ax_0(double x, double y)
    {
      return Legendre0x(x) * l0(y);
    }

    static double leg_quad_p0_e3_ax_1(double x, double y)
    {
      return -(Legendre0x(x) * l0(y));
    }

    static double leg_quad_p0_e3_ay_0(double x, double y)
    {
      return Legendre0(x) * dl0(y);
    }

    static double leg_quad_p0_e3_ay_1(double x, double y)
    {
      return -(Legendre0(x) * dl0(y));
    }

    static double leg_quad_p0_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e4_a_1(double x, double y)
    {
      return Legendre0(x) * l1(y);
    }

    static double leg_quad_p0_e4_a_0(double x, double y)
    {
      return -(Legendre0(x) * l1(y));
    }

    static double leg_quad_p0_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e4_ax_1(double x, double y)
    {
      return Legendre0x(x) * l1(y);
    }

    static double leg_quad_p0_e4_ax_0(double x, double y)
    {
      return -(Legendre0x(x) * l1(y));
    }

    static double leg_quad_p0_e4_ay_1(double x, double y)
    {
      return Legendre0(x) * dl1(y);
    }

    static double leg_quad_p0_e4_ay_0(double x, double y)
    {
      return -(Legendre0(x) * dl1(y));
    }

    static double leg_quad_p0_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 1 */

    static double leg_quad_p1_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e1_b(double x, double y)
    {
      return l0(x) * Legendre1(y);
    }

    static double leg_quad_p1_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e1_bx(double x, double y)
    {
      return dl0(x) * Legendre1(y);
    }

    static double leg_quad_p1_e1_by(double x, double y)
    {
      return l0(x) * Legendre1x(y);
    }

    static double leg_quad_p1_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e2_b(double x, double y)
    {
      return l1(x) * Legendre1(y);
    }

    static double leg_quad_p1_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e2_bx(double x, double y)
    {
      return dl1(x) * Legendre1(y);
    }

    static double leg_quad_p1_e2_by(double x, double y)
    {
      return l1(x) * Legendre1x(y);
    }

    static double leg_quad_p1_e3_a(double x, double y)
    {
      return Legendre1(x) * l0(y);
    }

    static double leg_quad_p1_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e3_ax(double x, double y)
    {
      return Legendre1x(x) * l0(y);
    }

    static double leg_quad_p1_e3_ay(double x, double y)
    {
      return Legendre1(x) * dl0(y);
    }

    static double leg_quad_p1_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e4_a(double x, double y)
    {
      return Legendre1(x) * l1(y);
    }

    static double leg_quad_p1_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e4_ax(double x, double y)
    {
      return Legendre1x(x) * l1(y);
    }

    static double leg_quad_p1_e4_ay(double x, double y)
    {
      return Legendre1(x) * dl1(y);
    }

    static double leg_quad_p1_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 2 */

    static double leg_quad_p2_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e1_b_1(double x, double y)
    {
      return l0(x) * Legendre2(y);
    }

    static double leg_quad_p2_e1_b_0(double x, double y)
    {
      return -(l0(x) * Legendre2(y));
    }

    static double leg_quad_p2_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e1_bx_1(double x, double y)
    {
      return dl0(x) * Legendre2(y);
    }

    static double leg_quad_p2_e1_bx_0(double x, double y)
    {
      return -(dl0(x) * Legendre2(y));
    }

    static double leg_quad_p2_e1_by_1(double x, double y)
    {
      return l0(x) * Legendre2x(y);
    }

    static double leg_quad_p2_e1_by_0(double x, double y)
    {
      return -(l0(x) * Legendre2x(y));
    }

    static double leg_quad_p2_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e2_b_0(double x, double y)
    {
      return l1(x) * Legendre2(y);
    }

    static double leg_quad_p2_e2_b_1(double x, double y)
    {
      return -(l1(x) * Legendre2(y));
    }

    static double leg_quad_p2_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e2_bx_0(double x, double y)
    {
      return dl1(x) * Legendre2(y);
    }

    static double leg_quad_p2_e2_bx_1(double x, double y)
    {
      return -(dl1(x) * Legendre2(y));
    }

    static double leg_quad_p2_e2_by_0(double x, double y)
    {
      return l1(x) * Legendre2x(y);
    }

    static double leg_quad_p2_e2_by_1(double x, double y)
    {
      return -(l1(x) * Legendre2x(y));
    }

    static double leg_quad_p2_e3_a_0(double x, double y)
    {
      return Legendre2(x) * l0(y);
    }

    static double leg_quad_p2_e3_a_1(double x, double y)
    {
      return -(Legendre2(x) * l0(y));
    }

    static double leg_quad_p2_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e3_ax_0(double x, double y)
    {
      return Legendre2x(x) * l0(y);
    }

    static double leg_quad_p2_e3_ax_1(double x, double y)
    {
      return -(Legendre2x(x) * l0(y));
    }

    static double leg_quad_p2_e3_ay_0(double x, double y)
    {
      return Legendre2(x) * dl0(y);
    }

    static double leg_quad_p2_e3_ay_1(double x, double y)
    {
      return -(Legendre2(x) * dl0(y));
    }

    static double leg_quad_p2_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e4_a_1(double x, double y)
    {
      return Legendre2(x) * l1(y);
    }

    static double leg_quad_p2_e4_a_0(double x, double y)
    {
      return -(Legendre2(x) * l1(y));
    }

    static double leg_quad_p2_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e4_ax_1(double x, double y)
    {
      return Legendre2x(x) * l1(y);
    }

    static double leg_quad_p2_e4_ax_0(double x, double y)
    {
      return -(Legendre2x(x) * l1(y));
    }

    static double leg_quad_p2_e4_ay_1(double x, double y)
    {
      return Legendre2(x) * dl1(y);
    }

    static double leg_quad_p2_e4_ay_0(double x, double y)
    {
      return -(Legendre2(x) * dl1(y));
    }

    static double leg_quad_p2_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 3 */

    static double leg_quad_p3_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e1_b(double x, double y)
    {
      return l0(x) * Legendre3(y);
    }

    static double leg_quad_p3_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e1_bx(double x, double y)
    {
      return dl0(x) * Legendre3(y);
    }

    static double leg_quad_p3_e1_by(double x, double y)
    {
      return l0(x) * Legendre3x(y);
    }

    static double leg_quad_p3_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e2_b(double x, double y)
    {
      return l1(x) * Legendre3(y);
    }

    static double leg_quad_p3_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e2_bx(double x, double y)
    {
      return dl1(x) * Legendre3(y);
    }

    static double leg_quad_p3_e2_by(double x, double y)
    {
      return l1(x) * Legendre3x(y);
    }

    static double leg_quad_p3_e3_a(double x, double y)
    {
      return Legendre3(x) * l0(y);
    }

    static double leg_quad_p3_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e3_ax(double x, double y)
    {
      return Legendre3x(x) * l0(y);
    }

    static double leg_quad_p3_e3_ay(double x, double y)
    {
      return Legendre3(x) * dl0(y);
    }

    static double leg_quad_p3_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e4_a(double x, double y)
    {
      return Legendre3(x) * l1(y);
    }

    static double leg_quad_p3_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e4_ax(double x, double y)
    {
      return Legendre3x(x) * l1(y);
    }

    static double leg_quad_p3_e4_ay(double x, double y)
    {
      return Legendre3(x) * dl1(y);
    }

    static double leg_quad_p3_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 4 */

    static double leg_quad_p4_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e1_b_1(double x, double y)
    {
      return l0(x) * Legendre4(y);
    }

    static double leg_quad_p4_e1_b_0(double x, double y)
    {
      return -(l0(x) * Legendre4(y));
    }

    static double leg_quad_p4_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e1_bx_1(double x, double y)
    {
      return dl0(x) * Legendre4(y);
    }

    static double leg_quad_p4_e1_bx_0(double x, double y)
    {
      return -(dl0(x) * Legendre4(y));
    }

    static double leg_quad_p4_e1_by_1(double x, double y)
    {
      return l0(x) * Legendre4x(y);
    }

    static double leg_quad_p4_e1_by_0(double x, double y)
    {
      return -(l0(x) * Legendre4x(y));
    }

    static double leg_quad_p4_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e2_b_0(double x, double y)
    {
      return l1(x) * Legendre4(y);
    }

    static double leg_quad_p4_e2_b_1(double x, double y)
    {
      return -(l1(x) * Legendre4(y));
    }

    static double leg_quad_p4_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e2_bx_0(double x, double y)
    {
      return dl1(x) * Legendre4(y);
    }

    static double leg_quad_p4_e2_bx_1(double x, double y)
    {
      return -(dl1(x) * Legendre4(y));
    }

    static double leg_quad_p4_e2_by_0(double x, double y)
    {
      return l1(x) * Legendre4x(y);
    }

    static double leg_quad_p4_e2_by_1(double x, double y)
    {
      return -(l1(x) * Legendre4x(y));
    }

    static double leg_quad_p4_e3_a_0(double x, double y)
    {
      return Legendre4(x) * l0(y);
    }

    static double leg_quad_p4_e3_a_1(double x, double y)
    {
      return -(Legendre4(x) * l0(y));
    }

    static double leg_quad_p4_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e3_ax_0(double x, double y)
    {
      return Legendre4x(x) * l0(y);
    }

    static double leg_quad_p4_e3_ax_1(double x, double y)
    {
      return -(Legendre4x(x) * l0(y));
    }

    static double leg_quad_p4_e3_ay_0(double x, double y)
    {
      return Legendre4(x) * dl0(y);
    }

    static double leg_quad_p4_e3_ay_1(double x, double y)
    {
      return -(Legendre4(x) * dl0(y));
    }

    static double leg_quad_p4_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e4_a_1(double x, double y)
    {
      return Legendre4(x) * l1(y);
    }

    static double leg_quad_p4_e4_a_0(double x, double y)
    {
      return -(Legendre4(x) * l1(y));
    }

    static double leg_quad_p4_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e4_ax_1(double x, double y)
    {
      return Legendre4x(x) * l1(y);
    }

    static double leg_quad_p4_e4_ax_0(double x, double y)
    {
      return -(Legendre4x(x) * l1(y));
    }

    static double leg_quad_p4_e4_ay_1(double x, double y)
    {
      return Legendre4(x) * dl1(y);
    }

    static double leg_quad_p4_e4_ay_0(double x, double y)
    {
      return -(Legendre4(x) * dl1(y));
    }

    static double leg_quad_p4_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 5 */

    static double leg_quad_p5_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e1_b(double x, double y)
    {
      return l0(x) * Legendre5(y);
    }

    static double leg_quad_p5_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e1_bx(double x, double y)
    {
      return dl0(x) * Legendre5(y);
    }

    static double leg_quad_p5_e1_by(double x, double y)
    {
      return l0(x) * Legendre5x(y);
    }

    static double leg_quad_p5_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e2_b(double x, double y)
    {
      return l1(x) * Legendre5(y);
    }

    static double leg_quad_p5_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e2_bx(double x, double y)
    {
      return dl1(x) * Legendre5(y);
    }

    static double leg_quad_p5_e2_by(double x, double y)
    {
      return l1(x) * Legendre5x(y);
    }

    static double leg_quad_p5_e3_a(double x, double y)
    {
      return Legendre5(x) * l0(y);
    }

    static double leg_quad_p5_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e3_ax(double x, double y)
    {
      return Legendre5x(x) * l0(y);
    }

    static double leg_quad_p5_e3_ay(double x, double y)
    {
      return Legendre5(x) * dl0(y);
    }

    static double leg_quad_p5_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e4_a(double x, double y)
    {
      return Legendre5(x) * l1(y);
    }

    static double leg_quad_p5_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e4_ax(double x, double y)
    {
      return Legendre5x(x) * l1(y);
    }

    static double leg_quad_p5_e4_ay(double x, double y)
    {
      return Legendre5(x) * dl1(y);
    }

    static double leg_quad_p5_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 6 */

    static double leg_quad_p6_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e1_b_1(double x, double y)
    {
      return l0(x) * Legendre6(y);
    }

    static double leg_quad_p6_e1_b_0(double x, double y)
    {
      return -(l0(x) * Legendre6(y));
    }

    static double leg_quad_p6_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e1_bx_1(double x, double y)
    {
      return dl0(x) * Legendre6(y);
    }

    static double leg_quad_p6_e1_bx_0(double x, double y)
    {
      return -(dl0(x) * Legendre6(y));
    }

    static double leg_quad_p6_e1_by_1(double x, double y)
    {
      return l0(x) * Legendre6x(y);
    }

    static double leg_quad_p6_e1_by_0(double x, double y)
    {
      return -(l0(x) * Legendre6x(y));
    }

    static double leg_quad_p6_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e2_b_0(double x, double y)
    {
      return l1(x) * Legendre6(y);
    }

    static double leg_quad_p6_e2_b_1(double x, double y)
    {
      return -(l1(x) * Legendre6(y));
    }

    static double leg_quad_p6_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e2_bx_0(double x, double y)
    {
      return dl1(x) * Legendre6(y);
    }

    static double leg_quad_p6_e2_bx_1(double x, double y)
    {
      return -(dl1(x) * Legendre6(y));
    }

    static double leg_quad_p6_e2_by_0(double x, double y)
    {
      return l1(x) * Legendre6x(y);
    }

    static double leg_quad_p6_e2_by_1(double x, double y)
    {
      return -(l1(x) * Legendre6x(y));
    }

    static double leg_quad_p6_e3_a_0(double x, double y)
    {
      return Legendre6(x) * l0(y);
    }

    static double leg_quad_p6_e3_a_1(double x, double y)
    {
      return -(Legendre6(x) * l0(y));
    }

    static double leg_quad_p6_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e3_ax_0(double x, double y)
    {
      return Legendre6x(x) * l0(y);
    }

    static double leg_quad_p6_e3_ax_1(double x, double y)
    {
      return -(Legendre6x(x) * l0(y));
    }

    static double leg_quad_p6_e3_ay_0(double x, double y)
    {
      return Legendre6(x) * dl0(y);
    }

    static double leg_quad_p6_e3_ay_1(double x, double y)
    {
      return -(Legendre6(x) * dl0(y));
    }

    static double leg_quad_p6_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e4_a_1(double x, double y)
    {
      return Legendre6(x) * l1(y);
    }

    static double leg_quad_p6_e4_a_0(double x, double y)
    {
      return -(Legendre6(x) * l1(y));
    }

    static double leg_quad_p6_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e4_ax_1(double x, double y)
    {
      return Legendre6x(x) * l1(y);
    }

    static double leg_quad_p6_e4_ax_0(double x, double y)
    {
      return -(Legendre6x(x) * l1(y));
    }

    static double leg_quad_p6_e4_ay_1(double x, double y)
    {
      return Legendre6(x) * dl1(y);
    }

    static double leg_quad_p6_e4_ay_0(double x, double y)
    {
      return -(Legendre6(x) * dl1(y));
    }

    static double leg_quad_p6_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 7 */

    static double leg_quad_p7_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e1_b(double x, double y)
    {
      return l0(x) * Legendre7(y);
    }

    static double leg_quad_p7_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e1_bx(double x, double y)
    {
      return dl0(x) * Legendre7(y);
    }

    static double leg_quad_p7_e1_by(double x, double y)
    {
      return l0(x) * Legendre7x(y);
    }

    static double leg_quad_p7_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e2_b(double x, double y)
    {
      return l1(x) * Legendre7(y);
    }

    static double leg_quad_p7_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e2_bx(double x, double y)
    {
      return dl1(x) * Legendre7(y);
    }

    static double leg_quad_p7_e2_by(double x, double y)
    {
      return l1(x) * Legendre7x(y);
    }

    static double leg_quad_p7_e3_a(double x, double y)
    {
      return Legendre7(x) * l0(y);
    }

    static double leg_quad_p7_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e3_ax(double x, double y)
    {
      return Legendre7x(x) * l0(y);
    }

    static double leg_quad_p7_e3_ay(double x, double y)
    {
      return Legendre7(x) * dl0(y);
    }

    static double leg_quad_p7_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e4_a(double x, double y)
    {
      return Legendre7(x) * l1(y);
    }

    static double leg_quad_p7_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e4_ax(double x, double y)
    {
      return Legendre7x(x) * l1(y);
    }

    static double leg_quad_p7_e4_ay(double x, double y)
    {
      return Legendre7(x) * dl1(y);
    }

    static double leg_quad_p7_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 8 */

    static double leg_quad_p8_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e1_b_1(double x, double y)
    {
      return l0(x) * Legendre8(y);
    }

    static double leg_quad_p8_e1_b_0(double x, double y)
    {
      return -(l0(x) * Legendre8(y));
    }

    static double leg_quad_p8_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e1_bx_1(double x, double y)
    {
      return dl0(x) * Legendre8(y);
    }

    static double leg_quad_p8_e1_bx_0(double x, double y)
    {
      return -(dl0(x) * Legendre8(y));
    }

    static double leg_quad_p8_e1_by_1(double x, double y)
    {
      return l0(x) * Legendre8x(y);
    }

    static double leg_quad_p8_e1_by_0(double x, double y)
    {
      return -(l0(x) * Legendre8x(y));
    }

    static double leg_quad_p8_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e2_b_0(double x, double y)
    {
      return l1(x) * Legendre8(y);
    }

    static double leg_quad_p8_e2_b_1(double x, double y)
    {
      return -(l1(x) * Legendre8(y));
    }

    static double leg_quad_p8_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e2_bx_0(double x, double y)
    {
      return dl1(x) * Legendre8(y);
    }

    static double leg_quad_p8_e2_bx_1(double x, double y)
    {
      return -(dl1(x) * Legendre8(y));
    }

    static double leg_quad_p8_e2_by_0(double x, double y)
    {
      return l1(x) * Legendre8x(y);
    }

    static double leg_quad_p8_e2_by_1(double x, double y)
    {
      return -(l1(x) * Legendre8x(y));
    }

    static double leg_quad_p8_e3_a_0(double x, double y)
    {
      return Legendre8(x) * l0(y);
    }

    static double leg_quad_p8_e3_a_1(double x, double y)
    {
      return -(Legendre8(x) * l0(y));
    }

    static double leg_quad_p8_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e3_ax_0(double x, double y)
    {
      return Legendre8x(x) * l0(y);
    }

    static double leg_quad_p8_e3_ax_1(double x, double y)
    {
      return -(Legendre8x(x) * l0(y));
    }

    static double leg_quad_p8_e3_ay_0(double x, double y)
    {
      return Legendre8(x) * dl0(y);
    }

    static double leg_quad_p8_e3_ay_1(double x, double y)
    {
      return -(Legendre8(x) * dl0(y));
    }

    static double leg_quad_p8_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e4_a_1(double x, double y)
    {
      return Legendre8(x) * l1(y);
    }

    static double leg_quad_p8_e4_a_0(double x, double y)
    {
      return -(Legendre8(x) * l1(y));
    }

    static double leg_quad_p8_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e4_ax_1(double x, double y)
    {
      return Legendre8x(x) * l1(y);
    }

    static double leg_quad_p8_e4_ax_0(double x, double y)
    {
      return -(Legendre8x(x) * l1(y));
    }

    static double leg_quad_p8_e4_ay_1(double x, double y)
    {
      return Legendre8(x) * dl1(y);
    }

    static double leg_quad_p8_e4_ay_0(double x, double y)
    {
      return -(Legendre8(x) * dl1(y));
    }

    static double leg_quad_p8_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 9 */

    static double leg_quad_p9_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e1_b(double x, double y)
    {
      return l0(x) * Legendre9(y);
    }

    static double leg_quad_p9_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e1_bx(double x, double y)
    {
      return dl0(x) * Legendre9(y);
    }

    static double leg_quad_p9_e1_by(double x, double y)
    {
      return l0(x) * Legendre9x(y);
    }

    static double leg_quad_p9_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e2_b(double x, double y)
    {
      return l1(x) * Legendre9(y);
    }

    static double leg_quad_p9_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e2_bx(double x, double y)
    {
      return dl1(x) * Legendre9(y);
    }

    static double leg_quad_p9_e2_by(double x, double y)
    {
      return l1(x) * Legendre9x(y);
    }

    static double leg_quad_p9_e3_a(double x, double y)
    {
      return Legendre9(x) * l0(y);
    }

    static double leg_quad_p9_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e3_ax(double x, double y)
    {
      return Legendre9x(x) * l0(y);
    }

    static double leg_quad_p9_e3_ay(double x, double y)
    {
      return Legendre9(x) * dl0(y);
    }

    static double leg_quad_p9_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e4_a(double x, double y)
    {
      return Legendre9(x) * l1(y);
    }

    static double leg_quad_p9_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e4_ax(double x, double y)
    {
      return Legendre9x(x) * l1(y);
    }

    static double leg_quad_p9_e4_ay(double x, double y)
    {
      return Legendre9(x) * dl1(y);
    }

    static double leg_quad_p9_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 10 */

    static double leg_quad_p10_e1_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e1_b_1(double x, double y)
    {
      return l0(x) * Legendre10(y);
    }

    static double leg_quad_p10_e1_b_0(double x, double y)
    {
      return -(l0(x) * Legendre10(y));
    }

    static double leg_quad_p10_e1_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e1_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e1_bx_1(double x, double y)
    {
      return dl0(x) * Legendre10(y);
    }

    static double leg_quad_p10_e1_bx_0(double x, double y)
    {
      return -(dl0(x) * Legendre10(y));
    }

    static double leg_quad_p10_e1_by_1(double x, double y)
    {
      return l0(x) * Legendre10x(y);
    }

    static double leg_quad_p10_e1_by_0(double x, double y)
    {
      return -(l0(x) * Legendre10x(y));
    }

    static double leg_quad_p10_e2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e2_b_0(double x, double y)
    {
      return l1(x) * Legendre10(y);
    }

    static double leg_quad_p10_e2_b_1(double x, double y)
    {
      return -(l1(x) * Legendre10(y));
    }

    static double leg_quad_p10_e2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e2_bx_0(double x, double y)
    {
      return dl1(x) * Legendre10(y);
    }

    static double leg_quad_p10_e2_bx_1(double x, double y)
    {
      return -(dl1(x) * Legendre10(y));
    }

    static double leg_quad_p10_e2_by_0(double x, double y)
    {
      return l1(x) * Legendre10x(y);
    }

    static double leg_quad_p10_e2_by_1(double x, double y)
    {
      return -(l1(x) * Legendre10x(y));
    }

    static double leg_quad_p10_e3_a_0(double x, double y)
    {
      return Legendre10(x) * l0(y);
    }

    static double leg_quad_p10_e3_a_1(double x, double y)
    {
      return -(Legendre10(x) * l0(y));
    }

    static double leg_quad_p10_e3_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e3_ax_0(double x, double y)
    {
      return Legendre10x(x) * l0(y);
    }

    static double leg_quad_p10_e3_ax_1(double x, double y)
    {
      return -(Legendre10x(x) * l0(y));
    }

    static double leg_quad_p10_e3_ay_0(double x, double y)
    {
      return Legendre10(x) * dl0(y);
    }

    static double leg_quad_p10_e3_ay_1(double x, double y)
    {
      return -(Legendre10(x) * dl0(y));
    }

    static double leg_quad_p10_e3_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e3_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e4_a_1(double x, double y)
    {
      return Legendre10(x) * l1(y);
    }

    static double leg_quad_p10_e4_a_0(double x, double y)
    {
      return -(Legendre10(x) * l1(y));
    }

    static double leg_quad_p10_e4_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e4_ax_1(double x, double y)
    {
      return Legendre10x(x) * l1(y);
    }

    static double leg_quad_p10_e4_ax_0(double x, double y)
    {
      return -(Legendre10x(x) * l1(y));
    }

    static double leg_quad_p10_e4_ay_1(double x, double y)
    {
      return Legendre10(x) * dl1(y);
    }

    static double leg_quad_p10_e4_ay_0(double x, double y)
    {
      return -(Legendre10(x) * dl1(y));
    }

    static double leg_quad_p10_e4_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10_e4_by(double x, double y)
    {
      return 0.0;
    }

    /* BUBBLE */

    /* BUBBLE ( 1 , 0 ) */

    static double leg_quad_p0p2_b1_a(double x, double y)
    {
      return Legendre0(x) * l2(y);
    }

    static double leg_quad_p0p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p2_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l2(y);
    }

    static double leg_quad_p0p2_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl2(y);
    }

    static double leg_quad_p0p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p3_b1_a(double x, double y)
    {
      return Legendre0(x) * l3(y);
    }

    static double leg_quad_p0p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p3_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l3(y);
    }

    static double leg_quad_p0p3_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl3(y);
    }

    static double leg_quad_p0p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p4_b1_a(double x, double y)
    {
      return Legendre0(x) * l4(y);
    }

    static double leg_quad_p0p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p4_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l4(y);
    }

    static double leg_quad_p0p4_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl4(y);
    }

    static double leg_quad_p0p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p5_b1_a(double x, double y)
    {
      return Legendre0(x) * l5(y);
    }

    static double leg_quad_p0p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p5_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l5(y);
    }

    static double leg_quad_p0p5_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl5(y);
    }

    static double leg_quad_p0p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p6_b1_a(double x, double y)
    {
      return Legendre0(x) * l6(y);
    }

    static double leg_quad_p0p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p6_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l6(y);
    }

    static double leg_quad_p0p6_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl6(y);
    }

    static double leg_quad_p0p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p7_b1_a(double x, double y)
    {
      return Legendre0(x) * l7(y);
    }

    static double leg_quad_p0p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p7_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l7(y);
    }

    static double leg_quad_p0p7_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl7(y);
    }

    static double leg_quad_p0p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p8_b1_a(double x, double y)
    {
      return Legendre0(x) * l8(y);
    }

    static double leg_quad_p0p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p8_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l8(y);
    }

    static double leg_quad_p0p8_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl8(y);
    }

    static double leg_quad_p0p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p9_b1_a(double x, double y)
    {
      return Legendre0(x) * l9(y);
    }

    static double leg_quad_p0p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p9_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l9(y);
    }

    static double leg_quad_p0p9_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl9(y);
    }

    static double leg_quad_p0p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p10_b1_a(double x, double y)
    {
      return Legendre0(x) * l10(y);
    }

    static double leg_quad_p0p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p10_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l10(y);
    }

    static double leg_quad_p0p10_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl10(y);
    }

    static double leg_quad_p0p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p11_b1_a(double x, double y)
    {
      return Legendre0(x) * l11(y);
    }

    static double leg_quad_p0p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p11_b1_ax(double x, double y)
    {
      return Legendre0x(x) * l11(y);
    }

    static double leg_quad_p0p11_b1_ay(double x, double y)
    {
      return Legendre0(x) * dl11(y);
    }

    static double leg_quad_p0p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p0p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p2_b1_a(double x, double y)
    {
      return Legendre1(x) * l2(y);
    }

    static double leg_quad_p1p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p2_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l2(y);
    }

    static double leg_quad_p1p2_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl2(y);
    }

    static double leg_quad_p1p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p3_b1_a(double x, double y)
    {
      return Legendre1(x) * l3(y);
    }

    static double leg_quad_p1p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p3_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l3(y);
    }

    static double leg_quad_p1p3_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl3(y);
    }

    static double leg_quad_p1p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p4_b1_a(double x, double y)
    {
      return Legendre1(x) * l4(y);
    }

    static double leg_quad_p1p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p4_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l4(y);
    }

    static double leg_quad_p1p4_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl4(y);
    }

    static double leg_quad_p1p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p5_b1_a(double x, double y)
    {
      return Legendre1(x) * l5(y);
    }

    static double leg_quad_p1p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p5_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l5(y);
    }

    static double leg_quad_p1p5_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl5(y);
    }

    static double leg_quad_p1p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p6_b1_a(double x, double y)
    {
      return Legendre1(x) * l6(y);
    }

    static double leg_quad_p1p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p6_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l6(y);
    }

    static double leg_quad_p1p6_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl6(y);
    }

    static double leg_quad_p1p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p7_b1_a(double x, double y)
    {
      return Legendre1(x) * l7(y);
    }

    static double leg_quad_p1p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p7_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l7(y);
    }

    static double leg_quad_p1p7_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl7(y);
    }

    static double leg_quad_p1p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p8_b1_a(double x, double y)
    {
      return Legendre1(x) * l8(y);
    }

    static double leg_quad_p1p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p8_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l8(y);
    }

    static double leg_quad_p1p8_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl8(y);
    }

    static double leg_quad_p1p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p9_b1_a(double x, double y)
    {
      return Legendre1(x) * l9(y);
    }

    static double leg_quad_p1p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p9_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l9(y);
    }

    static double leg_quad_p1p9_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl9(y);
    }

    static double leg_quad_p1p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p10_b1_a(double x, double y)
    {
      return Legendre1(x) * l10(y);
    }

    static double leg_quad_p1p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p10_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l10(y);
    }

    static double leg_quad_p1p10_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl10(y);
    }

    static double leg_quad_p1p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p11_b1_a(double x, double y)
    {
      return Legendre1(x) * l11(y);
    }

    static double leg_quad_p1p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p11_b1_ax(double x, double y)
    {
      return Legendre1x(x) * l11(y);
    }

    static double leg_quad_p1p11_b1_ay(double x, double y)
    {
      return Legendre1(x) * dl11(y);
    }

    static double leg_quad_p1p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p1p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p2_b1_a(double x, double y)
    {
      return Legendre2(x) * l2(y);
    }

    static double leg_quad_p2p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p2_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l2(y);
    }

    static double leg_quad_p2p2_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl2(y);
    }

    static double leg_quad_p2p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p3_b1_a(double x, double y)
    {
      return Legendre2(x) * l3(y);
    }

    static double leg_quad_p2p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p3_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l3(y);
    }

    static double leg_quad_p2p3_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl3(y);
    }

    static double leg_quad_p2p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p4_b1_a(double x, double y)
    {
      return Legendre2(x) * l4(y);
    }

    static double leg_quad_p2p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p4_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l4(y);
    }

    static double leg_quad_p2p4_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl4(y);
    }

    static double leg_quad_p2p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p5_b1_a(double x, double y)
    {
      return Legendre2(x) * l5(y);
    }

    static double leg_quad_p2p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p5_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l5(y);
    }

    static double leg_quad_p2p5_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl5(y);
    }

    static double leg_quad_p2p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p6_b1_a(double x, double y)
    {
      return Legendre2(x) * l6(y);
    }

    static double leg_quad_p2p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p6_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l6(y);
    }

    static double leg_quad_p2p6_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl6(y);
    }

    static double leg_quad_p2p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p7_b1_a(double x, double y)
    {
      return Legendre2(x) * l7(y);
    }

    static double leg_quad_p2p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p7_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l7(y);
    }

    static double leg_quad_p2p7_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl7(y);
    }

    static double leg_quad_p2p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p8_b1_a(double x, double y)
    {
      return Legendre2(x) * l8(y);
    }

    static double leg_quad_p2p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p8_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l8(y);
    }

    static double leg_quad_p2p8_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl8(y);
    }

    static double leg_quad_p2p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p9_b1_a(double x, double y)
    {
      return Legendre2(x) * l9(y);
    }

    static double leg_quad_p2p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p9_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l9(y);
    }

    static double leg_quad_p2p9_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl9(y);
    }

    static double leg_quad_p2p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p10_b1_a(double x, double y)
    {
      return Legendre2(x) * l10(y);
    }

    static double leg_quad_p2p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p10_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l10(y);
    }

    static double leg_quad_p2p10_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl10(y);
    }

    static double leg_quad_p2p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p11_b1_a(double x, double y)
    {
      return Legendre2(x) * l11(y);
    }

    static double leg_quad_p2p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p11_b1_ax(double x, double y)
    {
      return Legendre2x(x) * l11(y);
    }

    static double leg_quad_p2p11_b1_ay(double x, double y)
    {
      return Legendre2(x) * dl11(y);
    }

    static double leg_quad_p2p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p2_b1_a(double x, double y)
    {
      return Legendre3(x) * l2(y);
    }

    static double leg_quad_p3p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p2_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l2(y);
    }

    static double leg_quad_p3p2_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl2(y);
    }

    static double leg_quad_p3p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p3_b1_a(double x, double y)
    {
      return Legendre3(x) * l3(y);
    }

    static double leg_quad_p3p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p3_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l3(y);
    }

    static double leg_quad_p3p3_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl3(y);
    }

    static double leg_quad_p3p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p4_b1_a(double x, double y)
    {
      return Legendre3(x) * l4(y);
    }

    static double leg_quad_p3p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p4_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l4(y);
    }

    static double leg_quad_p3p4_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl4(y);
    }

    static double leg_quad_p3p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p5_b1_a(double x, double y)
    {
      return Legendre3(x) * l5(y);
    }

    static double leg_quad_p3p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p5_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l5(y);
    }

    static double leg_quad_p3p5_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl5(y);
    }

    static double leg_quad_p3p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p6_b1_a(double x, double y)
    {
      return Legendre3(x) * l6(y);
    }

    static double leg_quad_p3p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p6_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l6(y);
    }

    static double leg_quad_p3p6_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl6(y);
    }

    static double leg_quad_p3p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p7_b1_a(double x, double y)
    {
      return Legendre3(x) * l7(y);
    }

    static double leg_quad_p3p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p7_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l7(y);
    }

    static double leg_quad_p3p7_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl7(y);
    }

    static double leg_quad_p3p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p8_b1_a(double x, double y)
    {
      return Legendre3(x) * l8(y);
    }

    static double leg_quad_p3p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p8_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l8(y);
    }

    static double leg_quad_p3p8_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl8(y);
    }

    static double leg_quad_p3p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p9_b1_a(double x, double y)
    {
      return Legendre3(x) * l9(y);
    }

    static double leg_quad_p3p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p9_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l9(y);
    }

    static double leg_quad_p3p9_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl9(y);
    }

    static double leg_quad_p3p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p10_b1_a(double x, double y)
    {
      return Legendre3(x) * l10(y);
    }

    static double leg_quad_p3p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p10_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l10(y);
    }

    static double leg_quad_p3p10_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl10(y);
    }

    static double leg_quad_p3p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p11_b1_a(double x, double y)
    {
      return Legendre3(x) * l11(y);
    }

    static double leg_quad_p3p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p11_b1_ax(double x, double y)
    {
      return Legendre3x(x) * l11(y);
    }

    static double leg_quad_p3p11_b1_ay(double x, double y)
    {
      return Legendre3(x) * dl11(y);
    }

    static double leg_quad_p3p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p2_b1_a(double x, double y)
    {
      return Legendre4(x) * l2(y);
    }

    static double leg_quad_p4p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p2_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l2(y);
    }

    static double leg_quad_p4p2_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl2(y);
    }

    static double leg_quad_p4p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p3_b1_a(double x, double y)
    {
      return Legendre4(x) * l3(y);
    }

    static double leg_quad_p4p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p3_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l3(y);
    }

    static double leg_quad_p4p3_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl3(y);
    }

    static double leg_quad_p4p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p4_b1_a(double x, double y)
    {
      return Legendre4(x) * l4(y);
    }

    static double leg_quad_p4p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p4_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l4(y);
    }

    static double leg_quad_p4p4_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl4(y);
    }

    static double leg_quad_p4p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p5_b1_a(double x, double y)
    {
      return Legendre4(x) * l5(y);
    }

    static double leg_quad_p4p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p5_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l5(y);
    }

    static double leg_quad_p4p5_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl5(y);
    }

    static double leg_quad_p4p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p6_b1_a(double x, double y)
    {
      return Legendre4(x) * l6(y);
    }

    static double leg_quad_p4p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p6_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l6(y);
    }

    static double leg_quad_p4p6_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl6(y);
    }

    static double leg_quad_p4p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p7_b1_a(double x, double y)
    {
      return Legendre4(x) * l7(y);
    }

    static double leg_quad_p4p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p7_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l7(y);
    }

    static double leg_quad_p4p7_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl7(y);
    }

    static double leg_quad_p4p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p8_b1_a(double x, double y)
    {
      return Legendre4(x) * l8(y);
    }

    static double leg_quad_p4p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p8_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l8(y);
    }

    static double leg_quad_p4p8_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl8(y);
    }

    static double leg_quad_p4p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p9_b1_a(double x, double y)
    {
      return Legendre4(x) * l9(y);
    }

    static double leg_quad_p4p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p9_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l9(y);
    }

    static double leg_quad_p4p9_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl9(y);
    }

    static double leg_quad_p4p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p10_b1_a(double x, double y)
    {
      return Legendre4(x) * l10(y);
    }

    static double leg_quad_p4p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p10_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l10(y);
    }

    static double leg_quad_p4p10_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl10(y);
    }

    static double leg_quad_p4p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p11_b1_a(double x, double y)
    {
      return Legendre4(x) * l11(y);
    }

    static double leg_quad_p4p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p11_b1_ax(double x, double y)
    {
      return Legendre4x(x) * l11(y);
    }

    static double leg_quad_p4p11_b1_ay(double x, double y)
    {
      return Legendre4(x) * dl11(y);
    }

    static double leg_quad_p4p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p2_b1_a(double x, double y)
    {
      return Legendre5(x) * l2(y);
    }

    static double leg_quad_p5p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p2_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l2(y);
    }

    static double leg_quad_p5p2_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl2(y);
    }

    static double leg_quad_p5p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p3_b1_a(double x, double y)
    {
      return Legendre5(x) * l3(y);
    }

    static double leg_quad_p5p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p3_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l3(y);
    }

    static double leg_quad_p5p3_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl3(y);
    }

    static double leg_quad_p5p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p4_b1_a(double x, double y)
    {
      return Legendre5(x) * l4(y);
    }

    static double leg_quad_p5p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p4_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l4(y);
    }

    static double leg_quad_p5p4_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl4(y);
    }

    static double leg_quad_p5p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p5_b1_a(double x, double y)
    {
      return Legendre5(x) * l5(y);
    }

    static double leg_quad_p5p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p5_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l5(y);
    }

    static double leg_quad_p5p5_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl5(y);
    }

    static double leg_quad_p5p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p6_b1_a(double x, double y)
    {
      return Legendre5(x) * l6(y);
    }

    static double leg_quad_p5p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p6_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l6(y);
    }

    static double leg_quad_p5p6_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl6(y);
    }

    static double leg_quad_p5p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p7_b1_a(double x, double y)
    {
      return Legendre5(x) * l7(y);
    }

    static double leg_quad_p5p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p7_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l7(y);
    }

    static double leg_quad_p5p7_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl7(y);
    }

    static double leg_quad_p5p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p8_b1_a(double x, double y)
    {
      return Legendre5(x) * l8(y);
    }

    static double leg_quad_p5p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p8_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l8(y);
    }

    static double leg_quad_p5p8_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl8(y);
    }

    static double leg_quad_p5p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p9_b1_a(double x, double y)
    {
      return Legendre5(x) * l9(y);
    }

    static double leg_quad_p5p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p9_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l9(y);
    }

    static double leg_quad_p5p9_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl9(y);
    }

    static double leg_quad_p5p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p10_b1_a(double x, double y)
    {
      return Legendre5(x) * l10(y);
    }

    static double leg_quad_p5p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p10_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l10(y);
    }

    static double leg_quad_p5p10_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl10(y);
    }

    static double leg_quad_p5p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p11_b1_a(double x, double y)
    {
      return Legendre5(x) * l11(y);
    }

    static double leg_quad_p5p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p11_b1_ax(double x, double y)
    {
      return Legendre5x(x) * l11(y);
    }

    static double leg_quad_p5p11_b1_ay(double x, double y)
    {
      return Legendre5(x) * dl11(y);
    }

    static double leg_quad_p5p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p2_b1_a(double x, double y)
    {
      return Legendre6(x) * l2(y);
    }

    static double leg_quad_p6p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p2_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l2(y);
    }

    static double leg_quad_p6p2_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl2(y);
    }

    static double leg_quad_p6p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p3_b1_a(double x, double y)
    {
      return Legendre6(x) * l3(y);
    }

    static double leg_quad_p6p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p3_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l3(y);
    }

    static double leg_quad_p6p3_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl3(y);
    }

    static double leg_quad_p6p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p4_b1_a(double x, double y)
    {
      return Legendre6(x) * l4(y);
    }

    static double leg_quad_p6p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p4_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l4(y);
    }

    static double leg_quad_p6p4_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl4(y);
    }

    static double leg_quad_p6p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p5_b1_a(double x, double y)
    {
      return Legendre6(x) * l5(y);
    }

    static double leg_quad_p6p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p5_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l5(y);
    }

    static double leg_quad_p6p5_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl5(y);
    }

    static double leg_quad_p6p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p6_b1_a(double x, double y)
    {
      return Legendre6(x) * l6(y);
    }

    static double leg_quad_p6p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p6_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l6(y);
    }

    static double leg_quad_p6p6_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl6(y);
    }

    static double leg_quad_p6p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p7_b1_a(double x, double y)
    {
      return Legendre6(x) * l7(y);
    }

    static double leg_quad_p6p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p7_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l7(y);
    }

    static double leg_quad_p6p7_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl7(y);
    }

    static double leg_quad_p6p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p8_b1_a(double x, double y)
    {
      return Legendre6(x) * l8(y);
    }

    static double leg_quad_p6p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p8_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l8(y);
    }

    static double leg_quad_p6p8_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl8(y);
    }

    static double leg_quad_p6p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p9_b1_a(double x, double y)
    {
      return Legendre6(x) * l9(y);
    }

    static double leg_quad_p6p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p9_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l9(y);
    }

    static double leg_quad_p6p9_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl9(y);
    }

    static double leg_quad_p6p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p10_b1_a(double x, double y)
    {
      return Legendre6(x) * l10(y);
    }

    static double leg_quad_p6p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p10_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l10(y);
    }

    static double leg_quad_p6p10_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl10(y);
    }

    static double leg_quad_p6p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p11_b1_a(double x, double y)
    {
      return Legendre6(x) * l11(y);
    }

    static double leg_quad_p6p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p11_b1_ax(double x, double y)
    {
      return Legendre6x(x) * l11(y);
    }

    static double leg_quad_p6p11_b1_ay(double x, double y)
    {
      return Legendre6(x) * dl11(y);
    }

    static double leg_quad_p6p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p2_b1_a(double x, double y)
    {
      return Legendre7(x) * l2(y);
    }

    static double leg_quad_p7p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p2_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l2(y);
    }

    static double leg_quad_p7p2_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl2(y);
    }

    static double leg_quad_p7p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p3_b1_a(double x, double y)
    {
      return Legendre7(x) * l3(y);
    }

    static double leg_quad_p7p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p3_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l3(y);
    }

    static double leg_quad_p7p3_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl3(y);
    }

    static double leg_quad_p7p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p4_b1_a(double x, double y)
    {
      return Legendre7(x) * l4(y);
    }

    static double leg_quad_p7p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p4_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l4(y);
    }

    static double leg_quad_p7p4_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl4(y);
    }

    static double leg_quad_p7p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p5_b1_a(double x, double y)
    {
      return Legendre7(x) * l5(y);
    }

    static double leg_quad_p7p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p5_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l5(y);
    }

    static double leg_quad_p7p5_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl5(y);
    }

    static double leg_quad_p7p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p6_b1_a(double x, double y)
    {
      return Legendre7(x) * l6(y);
    }

    static double leg_quad_p7p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p6_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l6(y);
    }

    static double leg_quad_p7p6_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl6(y);
    }

    static double leg_quad_p7p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p7_b1_a(double x, double y)
    {
      return Legendre7(x) * l7(y);
    }

    static double leg_quad_p7p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p7_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l7(y);
    }

    static double leg_quad_p7p7_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl7(y);
    }

    static double leg_quad_p7p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p8_b1_a(double x, double y)
    {
      return Legendre7(x) * l8(y);
    }

    static double leg_quad_p7p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p8_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l8(y);
    }

    static double leg_quad_p7p8_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl8(y);
    }

    static double leg_quad_p7p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p9_b1_a(double x, double y)
    {
      return Legendre7(x) * l9(y);
    }

    static double leg_quad_p7p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p9_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l9(y);
    }

    static double leg_quad_p7p9_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl9(y);
    }

    static double leg_quad_p7p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p10_b1_a(double x, double y)
    {
      return Legendre7(x) * l10(y);
    }

    static double leg_quad_p7p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p10_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l10(y);
    }

    static double leg_quad_p7p10_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl10(y);
    }

    static double leg_quad_p7p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p11_b1_a(double x, double y)
    {
      return Legendre7(x) * l11(y);
    }

    static double leg_quad_p7p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p11_b1_ax(double x, double y)
    {
      return Legendre7x(x) * l11(y);
    }

    static double leg_quad_p7p11_b1_ay(double x, double y)
    {
      return Legendre7(x) * dl11(y);
    }

    static double leg_quad_p7p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p2_b1_a(double x, double y)
    {
      return Legendre8(x) * l2(y);
    }

    static double leg_quad_p8p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p2_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l2(y);
    }

    static double leg_quad_p8p2_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl2(y);
    }

    static double leg_quad_p8p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p3_b1_a(double x, double y)
    {
      return Legendre8(x) * l3(y);
    }

    static double leg_quad_p8p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p3_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l3(y);
    }

    static double leg_quad_p8p3_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl3(y);
    }

    static double leg_quad_p8p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p4_b1_a(double x, double y)
    {
      return Legendre8(x) * l4(y);
    }

    static double leg_quad_p8p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p4_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l4(y);
    }

    static double leg_quad_p8p4_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl4(y);
    }

    static double leg_quad_p8p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p5_b1_a(double x, double y)
    {
      return Legendre8(x) * l5(y);
    }

    static double leg_quad_p8p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p5_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l5(y);
    }

    static double leg_quad_p8p5_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl5(y);
    }

    static double leg_quad_p8p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p6_b1_a(double x, double y)
    {
      return Legendre8(x) * l6(y);
    }

    static double leg_quad_p8p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p6_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l6(y);
    }

    static double leg_quad_p8p6_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl6(y);
    }

    static double leg_quad_p8p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p7_b1_a(double x, double y)
    {
      return Legendre8(x) * l7(y);
    }

    static double leg_quad_p8p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p7_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l7(y);
    }

    static double leg_quad_p8p7_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl7(y);
    }

    static double leg_quad_p8p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p8_b1_a(double x, double y)
    {
      return Legendre8(x) * l8(y);
    }

    static double leg_quad_p8p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p8_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l8(y);
    }

    static double leg_quad_p8p8_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl8(y);
    }

    static double leg_quad_p8p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p9_b1_a(double x, double y)
    {
      return Legendre8(x) * l9(y);
    }

    static double leg_quad_p8p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p9_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l9(y);
    }

    static double leg_quad_p8p9_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl9(y);
    }

    static double leg_quad_p8p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p10_b1_a(double x, double y)
    {
      return Legendre8(x) * l10(y);
    }

    static double leg_quad_p8p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p10_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l10(y);
    }

    static double leg_quad_p8p10_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl10(y);
    }

    static double leg_quad_p8p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p11_b1_a(double x, double y)
    {
      return Legendre8(x) * l11(y);
    }

    static double leg_quad_p8p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p11_b1_ax(double x, double y)
    {
      return Legendre8x(x) * l11(y);
    }

    static double leg_quad_p8p11_b1_ay(double x, double y)
    {
      return Legendre8(x) * dl11(y);
    }

    static double leg_quad_p8p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p2_b1_a(double x, double y)
    {
      return Legendre9(x) * l2(y);
    }

    static double leg_quad_p9p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p2_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l2(y);
    }

    static double leg_quad_p9p2_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl2(y);
    }

    static double leg_quad_p9p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p3_b1_a(double x, double y)
    {
      return Legendre9(x) * l3(y);
    }

    static double leg_quad_p9p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p3_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l3(y);
    }

    static double leg_quad_p9p3_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl3(y);
    }

    static double leg_quad_p9p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p4_b1_a(double x, double y)
    {
      return Legendre9(x) * l4(y);
    }

    static double leg_quad_p9p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p4_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l4(y);
    }

    static double leg_quad_p9p4_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl4(y);
    }

    static double leg_quad_p9p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p5_b1_a(double x, double y)
    {
      return Legendre9(x) * l5(y);
    }

    static double leg_quad_p9p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p5_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l5(y);
    }

    static double leg_quad_p9p5_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl5(y);
    }

    static double leg_quad_p9p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p6_b1_a(double x, double y)
    {
      return Legendre9(x) * l6(y);
    }

    static double leg_quad_p9p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p6_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l6(y);
    }

    static double leg_quad_p9p6_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl6(y);
    }

    static double leg_quad_p9p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p7_b1_a(double x, double y)
    {
      return Legendre9(x) * l7(y);
    }

    static double leg_quad_p9p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p7_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l7(y);
    }

    static double leg_quad_p9p7_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl7(y);
    }

    static double leg_quad_p9p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p8_b1_a(double x, double y)
    {
      return Legendre9(x) * l8(y);
    }

    static double leg_quad_p9p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p8_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l8(y);
    }

    static double leg_quad_p9p8_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl8(y);
    }

    static double leg_quad_p9p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p9_b1_a(double x, double y)
    {
      return Legendre9(x) * l9(y);
    }

    static double leg_quad_p9p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p9_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l9(y);
    }

    static double leg_quad_p9p9_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl9(y);
    }

    static double leg_quad_p9p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p10_b1_a(double x, double y)
    {
      return Legendre9(x) * l10(y);
    }

    static double leg_quad_p9p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p10_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l10(y);
    }

    static double leg_quad_p9p10_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl10(y);
    }

    static double leg_quad_p9p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p11_b1_a(double x, double y)
    {
      return Legendre9(x) * l11(y);
    }

    static double leg_quad_p9p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p11_b1_ax(double x, double y)
    {
      return Legendre9x(x) * l11(y);
    }

    static double leg_quad_p9p11_b1_ay(double x, double y)
    {
      return Legendre9(x) * dl11(y);
    }

    static double leg_quad_p9p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p2_b1_a(double x, double y)
    {
      return Legendre10(x) * l2(y);
    }

    static double leg_quad_p10p2_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p2_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l2(y);
    }

    static double leg_quad_p10p2_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl2(y);
    }

    static double leg_quad_p10p2_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p2_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p3_b1_a(double x, double y)
    {
      return Legendre10(x) * l3(y);
    }

    static double leg_quad_p10p3_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p3_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l3(y);
    }

    static double leg_quad_p10p3_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl3(y);
    }

    static double leg_quad_p10p3_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p3_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p4_b1_a(double x, double y)
    {
      return Legendre10(x) * l4(y);
    }

    static double leg_quad_p10p4_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p4_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l4(y);
    }

    static double leg_quad_p10p4_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl4(y);
    }

    static double leg_quad_p10p4_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p4_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p5_b1_a(double x, double y)
    {
      return Legendre10(x) * l5(y);
    }

    static double leg_quad_p10p5_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p5_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l5(y);
    }

    static double leg_quad_p10p5_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl5(y);
    }

    static double leg_quad_p10p5_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p5_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p6_b1_a(double x, double y)
    {
      return Legendre10(x) * l6(y);
    }

    static double leg_quad_p10p6_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p6_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l6(y);
    }

    static double leg_quad_p10p6_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl6(y);
    }

    static double leg_quad_p10p6_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p6_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p7_b1_a(double x, double y)
    {
      return Legendre10(x) * l7(y);
    }

    static double leg_quad_p10p7_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p7_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l7(y);
    }

    static double leg_quad_p10p7_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl7(y);
    }

    static double leg_quad_p10p7_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p7_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p8_b1_a(double x, double y)
    {
      return Legendre10(x) * l8(y);
    }

    static double leg_quad_p10p8_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p8_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l8(y);
    }

    static double leg_quad_p10p8_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl8(y);
    }

    static double leg_quad_p10p8_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p8_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p9_b1_a(double x, double y)
    {
      return Legendre10(x) * l9(y);
    }

    static double leg_quad_p10p9_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p9_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l9(y);
    }

    static double leg_quad_p10p9_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl9(y);
    }

    static double leg_quad_p10p9_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p9_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p10_b1_a(double x, double y)
    {
      return Legendre10(x) * l10(y);
    }

    static double leg_quad_p10p10_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p10_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l10(y);
    }

    static double leg_quad_p10p10_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl10(y);
    }

    static double leg_quad_p10p10_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p10_b1_by(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p11_b1_a(double x, double y)
    {
      return Legendre10(x) * l11(y);
    }

    static double leg_quad_p10p11_b1_b(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p11_b1_ax(double x, double y)
    {
      return Legendre10x(x) * l11(y);
    }

    static double leg_quad_p10p11_b1_ay(double x, double y)
    {
      return Legendre10(x) * dl11(y);
    }

    static double leg_quad_p10p11_b1_bx(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p11_b1_by(double x, double y)
    {
      return 0.0;
    }

    /* BUBBLE ( 0 , 1 ) */

    static double leg_quad_p2p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p0_b2_b(double x, double y)
    {
      return l2(x) * Legendre0(y);
    }

    static double leg_quad_p2p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p0_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre0(y);
    }

    static double leg_quad_p2p0_b2_by(double x, double y)
    {
      return l2(x) * Legendre0x(y);
    }

    static double leg_quad_p2p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p1_b2_b(double x, double y)
    {
      return l2(x) * Legendre1(y);
    }

    static double leg_quad_p2p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p1_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre1(y);
    }

    static double leg_quad_p2p1_b2_by(double x, double y)
    {
      return l2(x) * Legendre1x(y);
    }

    static double leg_quad_p2p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p2_b2_b(double x, double y)
    {
      return l2(x) * Legendre2(y);
    }

    static double leg_quad_p2p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p2_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre2(y);
    }

    static double leg_quad_p2p2_b2_by(double x, double y)
    {
      return l2(x) * Legendre2x(y);
    }

    static double leg_quad_p2p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p3_b2_b(double x, double y)
    {
      return l2(x) * Legendre3(y);
    }

    static double leg_quad_p2p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p3_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre3(y);
    }

    static double leg_quad_p2p3_b2_by(double x, double y)
    {
      return l2(x) * Legendre3x(y);
    }

    static double leg_quad_p2p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p4_b2_b(double x, double y)
    {
      return l2(x) * Legendre4(y);
    }

    static double leg_quad_p2p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p4_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre4(y);
    }

    static double leg_quad_p2p4_b2_by(double x, double y)
    {
      return l2(x) * Legendre4x(y);
    }

    static double leg_quad_p2p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p5_b2_b(double x, double y)
    {
      return l2(x) * Legendre5(y);
    }

    static double leg_quad_p2p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p5_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre5(y);
    }

    static double leg_quad_p2p5_b2_by(double x, double y)
    {
      return l2(x) * Legendre5x(y);
    }

    static double leg_quad_p2p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p6_b2_b(double x, double y)
    {
      return l2(x) * Legendre6(y);
    }

    static double leg_quad_p2p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p6_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre6(y);
    }

    static double leg_quad_p2p6_b2_by(double x, double y)
    {
      return l2(x) * Legendre6x(y);
    }

    static double leg_quad_p2p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p7_b2_b(double x, double y)
    {
      return l2(x) * Legendre7(y);
    }

    static double leg_quad_p2p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p7_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre7(y);
    }

    static double leg_quad_p2p7_b2_by(double x, double y)
    {
      return l2(x) * Legendre7x(y);
    }

    static double leg_quad_p2p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p8_b2_b(double x, double y)
    {
      return l2(x) * Legendre8(y);
    }

    static double leg_quad_p2p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p8_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre8(y);
    }

    static double leg_quad_p2p8_b2_by(double x, double y)
    {
      return l2(x) * Legendre8x(y);
    }

    static double leg_quad_p2p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p9_b2_b(double x, double y)
    {
      return l2(x) * Legendre9(y);
    }

    static double leg_quad_p2p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p9_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre9(y);
    }

    static double leg_quad_p2p9_b2_by(double x, double y)
    {
      return l2(x) * Legendre9x(y);
    }

    static double leg_quad_p2p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p10_b2_b(double x, double y)
    {
      return l2(x) * Legendre10(y);
    }

    static double leg_quad_p2p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p2p10_b2_bx(double x, double y)
    {
      return dl2(x) * Legendre10(y);
    }

    static double leg_quad_p2p10_b2_by(double x, double y)
    {
      return l2(x) * Legendre10x(y);
    }

    static double leg_quad_p3p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p0_b2_b(double x, double y)
    {
      return l3(x) * Legendre0(y);
    }

    static double leg_quad_p3p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p0_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre0(y);
    }

    static double leg_quad_p3p0_b2_by(double x, double y)
    {
      return l3(x) * Legendre0x(y);
    }

    static double leg_quad_p3p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p1_b2_b(double x, double y)
    {
      return l3(x) * Legendre1(y);
    }

    static double leg_quad_p3p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p1_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre1(y);
    }

    static double leg_quad_p3p1_b2_by(double x, double y)
    {
      return l3(x) * Legendre1x(y);
    }

    static double leg_quad_p3p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p2_b2_b(double x, double y)
    {
      return l3(x) * Legendre2(y);
    }

    static double leg_quad_p3p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p2_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre2(y);
    }

    static double leg_quad_p3p2_b2_by(double x, double y)
    {
      return l3(x) * Legendre2x(y);
    }

    static double leg_quad_p3p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p3_b2_b(double x, double y)
    {
      return l3(x) * Legendre3(y);
    }

    static double leg_quad_p3p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p3_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre3(y);
    }

    static double leg_quad_p3p3_b2_by(double x, double y)
    {
      return l3(x) * Legendre3x(y);
    }

    static double leg_quad_p3p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p4_b2_b(double x, double y)
    {
      return l3(x) * Legendre4(y);
    }

    static double leg_quad_p3p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p4_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre4(y);
    }

    static double leg_quad_p3p4_b2_by(double x, double y)
    {
      return l3(x) * Legendre4x(y);
    }

    static double leg_quad_p3p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p5_b2_b(double x, double y)
    {
      return l3(x) * Legendre5(y);
    }

    static double leg_quad_p3p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p5_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre5(y);
    }

    static double leg_quad_p3p5_b2_by(double x, double y)
    {
      return l3(x) * Legendre5x(y);
    }

    static double leg_quad_p3p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p6_b2_b(double x, double y)
    {
      return l3(x) * Legendre6(y);
    }

    static double leg_quad_p3p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p6_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre6(y);
    }

    static double leg_quad_p3p6_b2_by(double x, double y)
    {
      return l3(x) * Legendre6x(y);
    }

    static double leg_quad_p3p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p7_b2_b(double x, double y)
    {
      return l3(x) * Legendre7(y);
    }

    static double leg_quad_p3p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p7_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre7(y);
    }

    static double leg_quad_p3p7_b2_by(double x, double y)
    {
      return l3(x) * Legendre7x(y);
    }

    static double leg_quad_p3p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p8_b2_b(double x, double y)
    {
      return l3(x) * Legendre8(y);
    }

    static double leg_quad_p3p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p8_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre8(y);
    }

    static double leg_quad_p3p8_b2_by(double x, double y)
    {
      return l3(x) * Legendre8x(y);
    }

    static double leg_quad_p3p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p9_b2_b(double x, double y)
    {
      return l3(x) * Legendre9(y);
    }

    static double leg_quad_p3p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p9_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre9(y);
    }

    static double leg_quad_p3p9_b2_by(double x, double y)
    {
      return l3(x) * Legendre9x(y);
    }

    static double leg_quad_p3p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p10_b2_b(double x, double y)
    {
      return l3(x) * Legendre10(y);
    }

    static double leg_quad_p3p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p3p10_b2_bx(double x, double y)
    {
      return dl3(x) * Legendre10(y);
    }

    static double leg_quad_p3p10_b2_by(double x, double y)
    {
      return l3(x) * Legendre10x(y);
    }

    static double leg_quad_p4p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p0_b2_b(double x, double y)
    {
      return l4(x) * Legendre0(y);
    }

    static double leg_quad_p4p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p0_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre0(y);
    }

    static double leg_quad_p4p0_b2_by(double x, double y)
    {
      return l4(x) * Legendre0x(y);
    }

    static double leg_quad_p4p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p1_b2_b(double x, double y)
    {
      return l4(x) * Legendre1(y);
    }

    static double leg_quad_p4p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p1_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre1(y);
    }

    static double leg_quad_p4p1_b2_by(double x, double y)
    {
      return l4(x) * Legendre1x(y);
    }

    static double leg_quad_p4p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p2_b2_b(double x, double y)
    {
      return l4(x) * Legendre2(y);
    }

    static double leg_quad_p4p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p2_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre2(y);
    }

    static double leg_quad_p4p2_b2_by(double x, double y)
    {
      return l4(x) * Legendre2x(y);
    }

    static double leg_quad_p4p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p3_b2_b(double x, double y)
    {
      return l4(x) * Legendre3(y);
    }

    static double leg_quad_p4p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p3_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre3(y);
    }

    static double leg_quad_p4p3_b2_by(double x, double y)
    {
      return l4(x) * Legendre3x(y);
    }

    static double leg_quad_p4p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p4_b2_b(double x, double y)
    {
      return l4(x) * Legendre4(y);
    }

    static double leg_quad_p4p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p4_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre4(y);
    }

    static double leg_quad_p4p4_b2_by(double x, double y)
    {
      return l4(x) * Legendre4x(y);
    }

    static double leg_quad_p4p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p5_b2_b(double x, double y)
    {
      return l4(x) * Legendre5(y);
    }

    static double leg_quad_p4p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p5_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre5(y);
    }

    static double leg_quad_p4p5_b2_by(double x, double y)
    {
      return l4(x) * Legendre5x(y);
    }

    static double leg_quad_p4p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p6_b2_b(double x, double y)
    {
      return l4(x) * Legendre6(y);
    }

    static double leg_quad_p4p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p6_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre6(y);
    }

    static double leg_quad_p4p6_b2_by(double x, double y)
    {
      return l4(x) * Legendre6x(y);
    }

    static double leg_quad_p4p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p7_b2_b(double x, double y)
    {
      return l4(x) * Legendre7(y);
    }

    static double leg_quad_p4p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p7_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre7(y);
    }

    static double leg_quad_p4p7_b2_by(double x, double y)
    {
      return l4(x) * Legendre7x(y);
    }

    static double leg_quad_p4p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p8_b2_b(double x, double y)
    {
      return l4(x) * Legendre8(y);
    }

    static double leg_quad_p4p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p8_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre8(y);
    }

    static double leg_quad_p4p8_b2_by(double x, double y)
    {
      return l4(x) * Legendre8x(y);
    }

    static double leg_quad_p4p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p9_b2_b(double x, double y)
    {
      return l4(x) * Legendre9(y);
    }

    static double leg_quad_p4p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p9_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre9(y);
    }

    static double leg_quad_p4p9_b2_by(double x, double y)
    {
      return l4(x) * Legendre9x(y);
    }

    static double leg_quad_p4p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p10_b2_b(double x, double y)
    {
      return l4(x) * Legendre10(y);
    }

    static double leg_quad_p4p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p4p10_b2_bx(double x, double y)
    {
      return dl4(x) * Legendre10(y);
    }

    static double leg_quad_p4p10_b2_by(double x, double y)
    {
      return l4(x) * Legendre10x(y);
    }

    static double leg_quad_p5p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p0_b2_b(double x, double y)
    {
      return l5(x) * Legendre0(y);
    }

    static double leg_quad_p5p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p0_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre0(y);
    }

    static double leg_quad_p5p0_b2_by(double x, double y)
    {
      return l5(x) * Legendre0x(y);
    }

    static double leg_quad_p5p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p1_b2_b(double x, double y)
    {
      return l5(x) * Legendre1(y);
    }

    static double leg_quad_p5p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p1_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre1(y);
    }

    static double leg_quad_p5p1_b2_by(double x, double y)
    {
      return l5(x) * Legendre1x(y);
    }

    static double leg_quad_p5p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p2_b2_b(double x, double y)
    {
      return l5(x) * Legendre2(y);
    }

    static double leg_quad_p5p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p2_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre2(y);
    }

    static double leg_quad_p5p2_b2_by(double x, double y)
    {
      return l5(x) * Legendre2x(y);
    }

    static double leg_quad_p5p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p3_b2_b(double x, double y)
    {
      return l5(x) * Legendre3(y);
    }

    static double leg_quad_p5p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p3_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre3(y);
    }

    static double leg_quad_p5p3_b2_by(double x, double y)
    {
      return l5(x) * Legendre3x(y);
    }

    static double leg_quad_p5p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p4_b2_b(double x, double y)
    {
      return l5(x) * Legendre4(y);
    }

    static double leg_quad_p5p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p4_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre4(y);
    }

    static double leg_quad_p5p4_b2_by(double x, double y)
    {
      return l5(x) * Legendre4x(y);
    }

    static double leg_quad_p5p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p5_b2_b(double x, double y)
    {
      return l5(x) * Legendre5(y);
    }

    static double leg_quad_p5p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p5_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre5(y);
    }

    static double leg_quad_p5p5_b2_by(double x, double y)
    {
      return l5(x) * Legendre5x(y);
    }

    static double leg_quad_p5p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p6_b2_b(double x, double y)
    {
      return l5(x) * Legendre6(y);
    }

    static double leg_quad_p5p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p6_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre6(y);
    }

    static double leg_quad_p5p6_b2_by(double x, double y)
    {
      return l5(x) * Legendre6x(y);
    }

    static double leg_quad_p5p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p7_b2_b(double x, double y)
    {
      return l5(x) * Legendre7(y);
    }

    static double leg_quad_p5p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p7_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre7(y);
    }

    static double leg_quad_p5p7_b2_by(double x, double y)
    {
      return l5(x) * Legendre7x(y);
    }

    static double leg_quad_p5p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p8_b2_b(double x, double y)
    {
      return l5(x) * Legendre8(y);
    }

    static double leg_quad_p5p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p8_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre8(y);
    }

    static double leg_quad_p5p8_b2_by(double x, double y)
    {
      return l5(x) * Legendre8x(y);
    }

    static double leg_quad_p5p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p9_b2_b(double x, double y)
    {
      return l5(x) * Legendre9(y);
    }

    static double leg_quad_p5p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p9_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre9(y);
    }

    static double leg_quad_p5p9_b2_by(double x, double y)
    {
      return l5(x) * Legendre9x(y);
    }

    static double leg_quad_p5p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p10_b2_b(double x, double y)
    {
      return l5(x) * Legendre10(y);
    }

    static double leg_quad_p5p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p5p10_b2_bx(double x, double y)
    {
      return dl5(x) * Legendre10(y);
    }

    static double leg_quad_p5p10_b2_by(double x, double y)
    {
      return l5(x) * Legendre10x(y);
    }

    static double leg_quad_p6p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p0_b2_b(double x, double y)
    {
      return l6(x) * Legendre0(y);
    }

    static double leg_quad_p6p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p0_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre0(y);
    }

    static double leg_quad_p6p0_b2_by(double x, double y)
    {
      return l6(x) * Legendre0x(y);
    }

    static double leg_quad_p6p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p1_b2_b(double x, double y)
    {
      return l6(x) * Legendre1(y);
    }

    static double leg_quad_p6p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p1_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre1(y);
    }

    static double leg_quad_p6p1_b2_by(double x, double y)
    {
      return l6(x) * Legendre1x(y);
    }

    static double leg_quad_p6p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p2_b2_b(double x, double y)
    {
      return l6(x) * Legendre2(y);
    }

    static double leg_quad_p6p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p2_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre2(y);
    }

    static double leg_quad_p6p2_b2_by(double x, double y)
    {
      return l6(x) * Legendre2x(y);
    }

    static double leg_quad_p6p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p3_b2_b(double x, double y)
    {
      return l6(x) * Legendre3(y);
    }

    static double leg_quad_p6p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p3_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre3(y);
    }

    static double leg_quad_p6p3_b2_by(double x, double y)
    {
      return l6(x) * Legendre3x(y);
    }

    static double leg_quad_p6p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p4_b2_b(double x, double y)
    {
      return l6(x) * Legendre4(y);
    }

    static double leg_quad_p6p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p4_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre4(y);
    }

    static double leg_quad_p6p4_b2_by(double x, double y)
    {
      return l6(x) * Legendre4x(y);
    }

    static double leg_quad_p6p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p5_b2_b(double x, double y)
    {
      return l6(x) * Legendre5(y);
    }

    static double leg_quad_p6p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p5_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre5(y);
    }

    static double leg_quad_p6p5_b2_by(double x, double y)
    {
      return l6(x) * Legendre5x(y);
    }

    static double leg_quad_p6p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p6_b2_b(double x, double y)
    {
      return l6(x) * Legendre6(y);
    }

    static double leg_quad_p6p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p6_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre6(y);
    }

    static double leg_quad_p6p6_b2_by(double x, double y)
    {
      return l6(x) * Legendre6x(y);
    }

    static double leg_quad_p6p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p7_b2_b(double x, double y)
    {
      return l6(x) * Legendre7(y);
    }

    static double leg_quad_p6p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p7_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre7(y);
    }

    static double leg_quad_p6p7_b2_by(double x, double y)
    {
      return l6(x) * Legendre7x(y);
    }

    static double leg_quad_p6p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p8_b2_b(double x, double y)
    {
      return l6(x) * Legendre8(y);
    }

    static double leg_quad_p6p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p8_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre8(y);
    }

    static double leg_quad_p6p8_b2_by(double x, double y)
    {
      return l6(x) * Legendre8x(y);
    }

    static double leg_quad_p6p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p9_b2_b(double x, double y)
    {
      return l6(x) * Legendre9(y);
    }

    static double leg_quad_p6p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p9_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre9(y);
    }

    static double leg_quad_p6p9_b2_by(double x, double y)
    {
      return l6(x) * Legendre9x(y);
    }

    static double leg_quad_p6p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p10_b2_b(double x, double y)
    {
      return l6(x) * Legendre10(y);
    }

    static double leg_quad_p6p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p6p10_b2_bx(double x, double y)
    {
      return dl6(x) * Legendre10(y);
    }

    static double leg_quad_p6p10_b2_by(double x, double y)
    {
      return l6(x) * Legendre10x(y);
    }

    static double leg_quad_p7p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p0_b2_b(double x, double y)
    {
      return l7(x) * Legendre0(y);
    }

    static double leg_quad_p7p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p0_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre0(y);
    }

    static double leg_quad_p7p0_b2_by(double x, double y)
    {
      return l7(x) * Legendre0x(y);
    }

    static double leg_quad_p7p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p1_b2_b(double x, double y)
    {
      return l7(x) * Legendre1(y);
    }

    static double leg_quad_p7p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p1_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre1(y);
    }

    static double leg_quad_p7p1_b2_by(double x, double y)
    {
      return l7(x) * Legendre1x(y);
    }

    static double leg_quad_p7p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p2_b2_b(double x, double y)
    {
      return l7(x) * Legendre2(y);
    }

    static double leg_quad_p7p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p2_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre2(y);
    }

    static double leg_quad_p7p2_b2_by(double x, double y)
    {
      return l7(x) * Legendre2x(y);
    }

    static double leg_quad_p7p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p3_b2_b(double x, double y)
    {
      return l7(x) * Legendre3(y);
    }

    static double leg_quad_p7p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p3_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre3(y);
    }

    static double leg_quad_p7p3_b2_by(double x, double y)
    {
      return l7(x) * Legendre3x(y);
    }

    static double leg_quad_p7p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p4_b2_b(double x, double y)
    {
      return l7(x) * Legendre4(y);
    }

    static double leg_quad_p7p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p4_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre4(y);
    }

    static double leg_quad_p7p4_b2_by(double x, double y)
    {
      return l7(x) * Legendre4x(y);
    }

    static double leg_quad_p7p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p5_b2_b(double x, double y)
    {
      return l7(x) * Legendre5(y);
    }

    static double leg_quad_p7p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p5_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre5(y);
    }

    static double leg_quad_p7p5_b2_by(double x, double y)
    {
      return l7(x) * Legendre5x(y);
    }

    static double leg_quad_p7p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p6_b2_b(double x, double y)
    {
      return l7(x) * Legendre6(y);
    }

    static double leg_quad_p7p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p6_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre6(y);
    }

    static double leg_quad_p7p6_b2_by(double x, double y)
    {
      return l7(x) * Legendre6x(y);
    }

    static double leg_quad_p7p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p7_b2_b(double x, double y)
    {
      return l7(x) * Legendre7(y);
    }

    static double leg_quad_p7p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p7_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre7(y);
    }

    static double leg_quad_p7p7_b2_by(double x, double y)
    {
      return l7(x) * Legendre7x(y);
    }

    static double leg_quad_p7p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p8_b2_b(double x, double y)
    {
      return l7(x) * Legendre8(y);
    }

    static double leg_quad_p7p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p8_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre8(y);
    }

    static double leg_quad_p7p8_b2_by(double x, double y)
    {
      return l7(x) * Legendre8x(y);
    }

    static double leg_quad_p7p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p9_b2_b(double x, double y)
    {
      return l7(x) * Legendre9(y);
    }

    static double leg_quad_p7p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p9_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre9(y);
    }

    static double leg_quad_p7p9_b2_by(double x, double y)
    {
      return l7(x) * Legendre9x(y);
    }

    static double leg_quad_p7p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p10_b2_b(double x, double y)
    {
      return l7(x) * Legendre10(y);
    }

    static double leg_quad_p7p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p7p10_b2_bx(double x, double y)
    {
      return dl7(x) * Legendre10(y);
    }

    static double leg_quad_p7p10_b2_by(double x, double y)
    {
      return l7(x) * Legendre10x(y);
    }

    static double leg_quad_p8p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p0_b2_b(double x, double y)
    {
      return l8(x) * Legendre0(y);
    }

    static double leg_quad_p8p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p0_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre0(y);
    }

    static double leg_quad_p8p0_b2_by(double x, double y)
    {
      return l8(x) * Legendre0x(y);
    }

    static double leg_quad_p8p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p1_b2_b(double x, double y)
    {
      return l8(x) * Legendre1(y);
    }

    static double leg_quad_p8p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p1_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre1(y);
    }

    static double leg_quad_p8p1_b2_by(double x, double y)
    {
      return l8(x) * Legendre1x(y);
    }

    static double leg_quad_p8p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p2_b2_b(double x, double y)
    {
      return l8(x) * Legendre2(y);
    }

    static double leg_quad_p8p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p2_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre2(y);
    }

    static double leg_quad_p8p2_b2_by(double x, double y)
    {
      return l8(x) * Legendre2x(y);
    }

    static double leg_quad_p8p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p3_b2_b(double x, double y)
    {
      return l8(x) * Legendre3(y);
    }

    static double leg_quad_p8p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p3_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre3(y);
    }

    static double leg_quad_p8p3_b2_by(double x, double y)
    {
      return l8(x) * Legendre3x(y);
    }

    static double leg_quad_p8p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p4_b2_b(double x, double y)
    {
      return l8(x) * Legendre4(y);
    }

    static double leg_quad_p8p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p4_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre4(y);
    }

    static double leg_quad_p8p4_b2_by(double x, double y)
    {
      return l8(x) * Legendre4x(y);
    }

    static double leg_quad_p8p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p5_b2_b(double x, double y)
    {
      return l8(x) * Legendre5(y);
    }

    static double leg_quad_p8p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p5_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre5(y);
    }

    static double leg_quad_p8p5_b2_by(double x, double y)
    {
      return l8(x) * Legendre5x(y);
    }

    static double leg_quad_p8p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p6_b2_b(double x, double y)
    {
      return l8(x) * Legendre6(y);
    }

    static double leg_quad_p8p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p6_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre6(y);
    }

    static double leg_quad_p8p6_b2_by(double x, double y)
    {
      return l8(x) * Legendre6x(y);
    }

    static double leg_quad_p8p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p7_b2_b(double x, double y)
    {
      return l8(x) * Legendre7(y);
    }

    static double leg_quad_p8p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p7_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre7(y);
    }

    static double leg_quad_p8p7_b2_by(double x, double y)
    {
      return l8(x) * Legendre7x(y);
    }

    static double leg_quad_p8p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p8_b2_b(double x, double y)
    {
      return l8(x) * Legendre8(y);
    }

    static double leg_quad_p8p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p8_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre8(y);
    }

    static double leg_quad_p8p8_b2_by(double x, double y)
    {
      return l8(x) * Legendre8x(y);
    }

    static double leg_quad_p8p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p9_b2_b(double x, double y)
    {
      return l8(x) * Legendre9(y);
    }

    static double leg_quad_p8p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p9_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre9(y);
    }

    static double leg_quad_p8p9_b2_by(double x, double y)
    {
      return l8(x) * Legendre9x(y);
    }

    static double leg_quad_p8p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p10_b2_b(double x, double y)
    {
      return l8(x) * Legendre10(y);
    }

    static double leg_quad_p8p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p8p10_b2_bx(double x, double y)
    {
      return dl8(x) * Legendre10(y);
    }

    static double leg_quad_p8p10_b2_by(double x, double y)
    {
      return l8(x) * Legendre10x(y);
    }

    static double leg_quad_p9p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p0_b2_b(double x, double y)
    {
      return l9(x) * Legendre0(y);
    }

    static double leg_quad_p9p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p0_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre0(y);
    }

    static double leg_quad_p9p0_b2_by(double x, double y)
    {
      return l9(x) * Legendre0x(y);
    }

    static double leg_quad_p9p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p1_b2_b(double x, double y)
    {
      return l9(x) * Legendre1(y);
    }

    static double leg_quad_p9p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p1_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre1(y);
    }

    static double leg_quad_p9p1_b2_by(double x, double y)
    {
      return l9(x) * Legendre1x(y);
    }

    static double leg_quad_p9p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p2_b2_b(double x, double y)
    {
      return l9(x) * Legendre2(y);
    }

    static double leg_quad_p9p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p2_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre2(y);
    }

    static double leg_quad_p9p2_b2_by(double x, double y)
    {
      return l9(x) * Legendre2x(y);
    }

    static double leg_quad_p9p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p3_b2_b(double x, double y)
    {
      return l9(x) * Legendre3(y);
    }

    static double leg_quad_p9p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p3_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre3(y);
    }

    static double leg_quad_p9p3_b2_by(double x, double y)
    {
      return l9(x) * Legendre3x(y);
    }

    static double leg_quad_p9p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p4_b2_b(double x, double y)
    {
      return l9(x) * Legendre4(y);
    }

    static double leg_quad_p9p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p4_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre4(y);
    }

    static double leg_quad_p9p4_b2_by(double x, double y)
    {
      return l9(x) * Legendre4x(y);
    }

    static double leg_quad_p9p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p5_b2_b(double x, double y)
    {
      return l9(x) * Legendre5(y);
    }

    static double leg_quad_p9p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p5_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre5(y);
    }

    static double leg_quad_p9p5_b2_by(double x, double y)
    {
      return l9(x) * Legendre5x(y);
    }

    static double leg_quad_p9p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p6_b2_b(double x, double y)
    {
      return l9(x) * Legendre6(y);
    }

    static double leg_quad_p9p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p6_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre6(y);
    }

    static double leg_quad_p9p6_b2_by(double x, double y)
    {
      return l9(x) * Legendre6x(y);
    }

    static double leg_quad_p9p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p7_b2_b(double x, double y)
    {
      return l9(x) * Legendre7(y);
    }

    static double leg_quad_p9p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p7_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre7(y);
    }

    static double leg_quad_p9p7_b2_by(double x, double y)
    {
      return l9(x) * Legendre7x(y);
    }

    static double leg_quad_p9p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p8_b2_b(double x, double y)
    {
      return l9(x) * Legendre8(y);
    }

    static double leg_quad_p9p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p8_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre8(y);
    }

    static double leg_quad_p9p8_b2_by(double x, double y)
    {
      return l9(x) * Legendre8x(y);
    }

    static double leg_quad_p9p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p9_b2_b(double x, double y)
    {
      return l9(x) * Legendre9(y);
    }

    static double leg_quad_p9p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p9_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre9(y);
    }

    static double leg_quad_p9p9_b2_by(double x, double y)
    {
      return l9(x) * Legendre9x(y);
    }

    static double leg_quad_p9p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p10_b2_b(double x, double y)
    {
      return l9(x) * Legendre10(y);
    }

    static double leg_quad_p9p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p9p10_b2_bx(double x, double y)
    {
      return dl9(x) * Legendre10(y);
    }

    static double leg_quad_p9p10_b2_by(double x, double y)
    {
      return l9(x) * Legendre10x(y);
    }

    static double leg_quad_p10p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p0_b2_b(double x, double y)
    {
      return l10(x) * Legendre0(y);
    }

    static double leg_quad_p10p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p0_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre0(y);
    }

    static double leg_quad_p10p0_b2_by(double x, double y)
    {
      return l10(x) * Legendre0x(y);
    }

    static double leg_quad_p10p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p1_b2_b(double x, double y)
    {
      return l10(x) * Legendre1(y);
    }

    static double leg_quad_p10p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p1_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre1(y);
    }

    static double leg_quad_p10p1_b2_by(double x, double y)
    {
      return l10(x) * Legendre1x(y);
    }

    static double leg_quad_p10p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p2_b2_b(double x, double y)
    {
      return l10(x) * Legendre2(y);
    }

    static double leg_quad_p10p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p2_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre2(y);
    }

    static double leg_quad_p10p2_b2_by(double x, double y)
    {
      return l10(x) * Legendre2x(y);
    }

    static double leg_quad_p10p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p3_b2_b(double x, double y)
    {
      return l10(x) * Legendre3(y);
    }

    static double leg_quad_p10p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p3_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre3(y);
    }

    static double leg_quad_p10p3_b2_by(double x, double y)
    {
      return l10(x) * Legendre3x(y);
    }

    static double leg_quad_p10p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p4_b2_b(double x, double y)
    {
      return l10(x) * Legendre4(y);
    }

    static double leg_quad_p10p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p4_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre4(y);
    }

    static double leg_quad_p10p4_b2_by(double x, double y)
    {
      return l10(x) * Legendre4x(y);
    }

    static double leg_quad_p10p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p5_b2_b(double x, double y)
    {
      return l10(x) * Legendre5(y);
    }

    static double leg_quad_p10p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p5_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre5(y);
    }

    static double leg_quad_p10p5_b2_by(double x, double y)
    {
      return l10(x) * Legendre5x(y);
    }

    static double leg_quad_p10p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p6_b2_b(double x, double y)
    {
      return l10(x) * Legendre6(y);
    }

    static double leg_quad_p10p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p6_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre6(y);
    }

    static double leg_quad_p10p6_b2_by(double x, double y)
    {
      return l10(x) * Legendre6x(y);
    }

    static double leg_quad_p10p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p7_b2_b(double x, double y)
    {
      return l10(x) * Legendre7(y);
    }

    static double leg_quad_p10p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p7_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre7(y);
    }

    static double leg_quad_p10p7_b2_by(double x, double y)
    {
      return l10(x) * Legendre7x(y);
    }

    static double leg_quad_p10p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p8_b2_b(double x, double y)
    {
      return l10(x) * Legendre8(y);
    }

    static double leg_quad_p10p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p8_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre8(y);
    }

    static double leg_quad_p10p8_b2_by(double x, double y)
    {
      return l10(x) * Legendre8x(y);
    }

    static double leg_quad_p10p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p9_b2_b(double x, double y)
    {
      return l10(x) * Legendre9(y);
    }

    static double leg_quad_p10p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p9_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre9(y);
    }

    static double leg_quad_p10p9_b2_by(double x, double y)
    {
      return l10(x) * Legendre9x(y);
    }

    static double leg_quad_p10p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p10_b2_b(double x, double y)
    {
      return l10(x) * Legendre10(y);
    }

    static double leg_quad_p10p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p10p10_b2_bx(double x, double y)
    {
      return dl10(x) * Legendre10(y);
    }

    static double leg_quad_p10p10_b2_by(double x, double y)
    {
      return l10(x) * Legendre10x(y);
    }

    static double leg_quad_p11p0_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p0_b2_b(double x, double y)
    {
      return l11(x) * Legendre0(y);
    }

    static double leg_quad_p11p0_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p0_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p0_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre0(y);
    }

    static double leg_quad_p11p0_b2_by(double x, double y)
    {
      return l11(x) * Legendre0x(y);
    }

    static double leg_quad_p11p1_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p1_b2_b(double x, double y)
    {
      return l11(x) * Legendre1(y);
    }

    static double leg_quad_p11p1_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p1_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p1_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre1(y);
    }

    static double leg_quad_p11p1_b2_by(double x, double y)
    {
      return l11(x) * Legendre1x(y);
    }

    static double leg_quad_p11p2_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p2_b2_b(double x, double y)
    {
      return l11(x) * Legendre2(y);
    }

    static double leg_quad_p11p2_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p2_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p2_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre2(y);
    }

    static double leg_quad_p11p2_b2_by(double x, double y)
    {
      return l11(x) * Legendre2x(y);
    }

    static double leg_quad_p11p3_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p3_b2_b(double x, double y)
    {
      return l11(x) * Legendre3(y);
    }

    static double leg_quad_p11p3_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p3_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p3_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre3(y);
    }

    static double leg_quad_p11p3_b2_by(double x, double y)
    {
      return l11(x) * Legendre3x(y);
    }

    static double leg_quad_p11p4_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p4_b2_b(double x, double y)
    {
      return l11(x) * Legendre4(y);
    }

    static double leg_quad_p11p4_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p4_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p4_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre4(y);
    }

    static double leg_quad_p11p4_b2_by(double x, double y)
    {
      return l11(x) * Legendre4x(y);
    }

    static double leg_quad_p11p5_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p5_b2_b(double x, double y)
    {
      return l11(x) * Legendre5(y);
    }

    static double leg_quad_p11p5_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p5_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p5_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre5(y);
    }

    static double leg_quad_p11p5_b2_by(double x, double y)
    {
      return l11(x) * Legendre5x(y);
    }

    static double leg_quad_p11p6_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p6_b2_b(double x, double y)
    {
      return l11(x) * Legendre6(y);
    }

    static double leg_quad_p11p6_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p6_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p6_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre6(y);
    }

    static double leg_quad_p11p6_b2_by(double x, double y)
    {
      return l11(x) * Legendre6x(y);
    }

    static double leg_quad_p11p7_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p7_b2_b(double x, double y)
    {
      return l11(x) * Legendre7(y);
    }

    static double leg_quad_p11p7_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p7_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p7_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre7(y);
    }

    static double leg_quad_p11p7_b2_by(double x, double y)
    {
      return l11(x) * Legendre7x(y);
    }

    static double leg_quad_p11p8_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p8_b2_b(double x, double y)
    {
      return l11(x) * Legendre8(y);
    }

    static double leg_quad_p11p8_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p8_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p8_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre8(y);
    }

    static double leg_quad_p11p8_b2_by(double x, double y)
    {
      return l11(x) * Legendre8x(y);
    }

    static double leg_quad_p11p9_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p9_b2_b(double x, double y)
    {
      return l11(x) * Legendre9(y);
    }

    static double leg_quad_p11p9_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p9_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p9_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre9(y);
    }

    static double leg_quad_p11p9_b2_by(double x, double y)
    {
      return l11(x) * Legendre9x(y);
    }

    static double leg_quad_p11p10_b2_a(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p10_b2_b(double x, double y)
    {
      return l11(x) * Legendre10(y);
    }

    static double leg_quad_p11p10_b2_ax(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p10_b2_ay(double x, double y)
    {
      return 0.0;
    }

    static double leg_quad_p11p10_b2_bx(double x, double y)
    {
      return dl11(x) * Legendre10(y);
    }

    static double leg_quad_p11p10_b2_by(double x, double y)
    {
      return l11(x) * Legendre10x(y);
    }

    static Shapeset::shape_fn_t leg_quad_fn_a[] =
    {
      leg_quad_p0_e1_a, leg_quad_p0_e1_a, leg_quad_p0_e2_a, leg_quad_p0_e2_a, leg_quad_p0_e3_a_0, leg_quad_p0_e3_a_1, leg_quad_p0_e4_a_0, leg_quad_p0_e4_a_1,
      leg_quad_p1_e1_a, leg_quad_p1_e1_a, leg_quad_p1_e2_a, leg_quad_p1_e2_a, leg_quad_p1_e3_a, leg_quad_p1_e3_a, leg_quad_p1_e4_a, leg_quad_p1_e4_a,
      leg_quad_p2_e1_a, leg_quad_p2_e1_a, leg_quad_p2_e2_a, leg_quad_p2_e2_a, leg_quad_p2_e3_a_0, leg_quad_p2_e3_a_1, leg_quad_p2_e4_a_0, leg_quad_p2_e4_a_1,
      leg_quad_p3_e1_a, leg_quad_p3_e1_a, leg_quad_p3_e2_a, leg_quad_p3_e2_a, leg_quad_p3_e3_a, leg_quad_p3_e3_a, leg_quad_p3_e4_a, leg_quad_p3_e4_a,
      leg_quad_p4_e1_a, leg_quad_p4_e1_a, leg_quad_p4_e2_a, leg_quad_p4_e2_a, leg_quad_p4_e3_a_0, leg_quad_p4_e3_a_1, leg_quad_p4_e4_a_0, leg_quad_p4_e4_a_1,
      leg_quad_p5_e1_a, leg_quad_p5_e1_a, leg_quad_p5_e2_a, leg_quad_p5_e2_a, leg_quad_p5_e3_a, leg_quad_p5_e3_a, leg_quad_p5_e4_a, leg_quad_p5_e4_a,
      leg_quad_p6_e1_a, leg_quad_p6_e1_a, leg_quad_p6_e2_a, leg_quad_p6_e2_a, leg_quad_p6_e3_a_0, leg_quad_p6_e3_a_1, leg_quad_p6_e4_a_0, leg_quad_p6_e4_a_1,
      leg_quad_p7_e1_a, leg_quad_p7_e1_a, leg_quad_p7_e2_a, leg_quad_p7_e2_a, leg_quad_p7_e3_a, leg_quad_p7_e3_a, leg_quad_p7_e4_a, leg_quad_p7_e4_a,
      leg_quad_p8_e1_a, leg_quad_p8_e1_a, leg_quad_p8_e2_a, leg_quad_p8_e2_a, leg_quad_p8_e3_a_0, leg_quad_p8_e3_a_1, leg_quad_p8_e4_a_0, leg_quad_p8_e4_a_1,
      leg_quad_p9_e1_a, leg_quad_p9_e1_a, leg_quad_p9_e2_a, leg_quad_p9_e2_a, leg_quad_p9_e3_a, leg_quad_p9_e3_a, leg_quad_p9_e4_a, leg_quad_p9_e4_a,
      leg_quad_p10_e1_a, leg_quad_p10_e1_a, leg_quad_p10_e2_a, leg_quad_p10_e2_a, leg_quad_p10_e3_a_0, leg_quad_p10_e3_a_1, leg_quad_p10_e4_a_0, leg_quad_p10_e4_a_1,

      leg_quad_p0p2_b1_a,   leg_quad_p0p3_b1_a,   leg_quad_p0p4_b1_a,   leg_quad_p0p5_b1_a,   leg_quad_p0p6_b1_a,   leg_quad_p0p7_b1_a,   leg_quad_p0p8_b1_a,   leg_quad_p0p9_b1_a,   leg_quad_p0p10_b1_a,   leg_quad_p0p11_b1_a,   leg_quad_p1p2_b1_a,   leg_quad_p1p3_b1_a,   leg_quad_p1p4_b1_a,   leg_quad_p1p5_b1_a,   leg_quad_p1p6_b1_a,   leg_quad_p1p7_b1_a,   leg_quad_p1p8_b1_a,   leg_quad_p1p9_b1_a,   leg_quad_p1p10_b1_a,   leg_quad_p1p11_b1_a,   leg_quad_p2p2_b1_a,   leg_quad_p2p3_b1_a,   leg_quad_p2p4_b1_a,   leg_quad_p2p5_b1_a,   leg_quad_p2p6_b1_a,   leg_quad_p2p7_b1_a,   leg_quad_p2p8_b1_a,   leg_quad_p2p9_b1_a,   leg_quad_p2p10_b1_a,   leg_quad_p2p11_b1_a,   leg_quad_p3p2_b1_a,   leg_quad_p3p3_b1_a,   leg_quad_p3p4_b1_a,   leg_quad_p3p5_b1_a,   leg_quad_p3p6_b1_a,   leg_quad_p3p7_b1_a,   leg_quad_p3p8_b1_a,   leg_quad_p3p9_b1_a,   leg_quad_p3p10_b1_a,   leg_quad_p3p11_b1_a,   leg_quad_p4p2_b1_a,   leg_quad_p4p3_b1_a,   leg_quad_p4p4_b1_a,   leg_quad_p4p5_b1_a,   leg_quad_p4p6_b1_a,   leg_quad_p4p7_b1_a,   leg_quad_p4p8_b1_a,   leg_quad_p4p9_b1_a,   leg_quad_p4p10_b1_a,   leg_quad_p4p11_b1_a,   leg_quad_p5p2_b1_a,   leg_quad_p5p3_b1_a,   leg_quad_p5p4_b1_a,   leg_quad_p5p5_b1_a,   leg_quad_p5p6_b1_a,   leg_quad_p5p7_b1_a,   leg_quad_p5p8_b1_a,   leg_quad_p5p9_b1_a,   leg_quad_p5p10_b1_a,   leg_quad_p5p11_b1_a,   leg_quad_p6p2_b1_a,   leg_quad_p6p3_b1_a,   leg_quad_p6p4_b1_a,   leg_quad_p6p5_b1_a,   leg_quad_p6p6_b1_a,   leg_quad_p6p7_b1_a,   leg_quad_p6p8_b1_a,   leg_quad_p6p9_b1_a,   leg_quad_p6p10_b1_a,   leg_quad_p6p11_b1_a,   leg_quad_p7p2_b1_a,   leg_quad_p7p3_b1_a,   leg_quad_p7p4_b1_a,   leg_quad_p7p5_b1_a,   leg_quad_p7p6_b1_a,   leg_quad_p7p7_b1_a,   leg_quad_p7p8_b1_a,   leg_quad_p7p9_b1_a,   leg_quad_p7p10_b1_a,   leg_quad_p7p11_b1_a,   leg_quad_p8p2_b1_a,   leg_quad_p8p3_b1_a,   leg_quad_p8p4_b1_a,   leg_quad_p8p5_b1_a,   leg_quad_p8p6_b1_a,   leg_quad_p8p7_b1_a,   leg_quad_p8p8_b1_a,   leg_quad_p8p9_b1_a,   leg_quad_p8p10_b1_a,   leg_quad_p8p11_b1_a,   leg_quad_p9p2_b1_a,   leg_quad_p9p3_b1_a,   leg_quad_p9p4_b1_a,   leg_quad_p9p5_b1_a,   leg_quad_p9p6_b1_a,   leg_quad_p9p7_b1_a,   leg_quad_p9p8_b1_a,   leg_quad_p9p9_b1_a,   leg_quad_p9p10_b1_a,   leg_quad_p9p11_b1_a,   leg_quad_p10p2_b1_a,   leg_quad_p10p3_b1_a,   leg_quad_p10p4_b1_a,   leg_quad_p10p5_b1_a,   leg_quad_p10p6_b1_a,   leg_quad_p10p7_b1_a,   leg_quad_p10p8_b1_a,   leg_quad_p10p9_b1_a,   leg_quad_p10p10_b1_a,   leg_quad_p10p11_b1_a,   leg_quad_p2p0_b2_a,   leg_quad_p2p1_b2_a,   leg_quad_p2p2_b2_a,   leg_quad_p2p3_b2_a,   leg_quad_p2p4_b2_a,   leg_quad_p2p5_b2_a,   leg_quad_p2p6_b2_a,   leg_quad_p2p7_b2_a,   leg_quad_p2p8_b2_a,   leg_quad_p2p9_b2_a,   leg_quad_p2p10_b2_a,   leg_quad_p3p0_b2_a,   leg_quad_p3p1_b2_a,   leg_quad_p3p2_b2_a,   leg_quad_p3p3_b2_a,   leg_quad_p3p4_b2_a,   leg_quad_p3p5_b2_a,   leg_quad_p3p6_b2_a,   leg_quad_p3p7_b2_a,   leg_quad_p3p8_b2_a,   leg_quad_p3p9_b2_a,   leg_quad_p3p10_b2_a,   leg_quad_p4p0_b2_a,   leg_quad_p4p1_b2_a,   leg_quad_p4p2_b2_a,   leg_quad_p4p3_b2_a,   leg_quad_p4p4_b2_a,   leg_quad_p4p5_b2_a,   leg_quad_p4p6_b2_a,   leg_quad_p4p7_b2_a,   leg_quad_p4p8_b2_a,   leg_quad_p4p9_b2_a,   leg_quad_p4p10_b2_a,   leg_quad_p5p0_b2_a,   leg_quad_p5p1_b2_a,   leg_quad_p5p2_b2_a,   leg_quad_p5p3_b2_a,   leg_quad_p5p4_b2_a,   leg_quad_p5p5_b2_a,   leg_quad_p5p6_b2_a,   leg_quad_p5p7_b2_a,   leg_quad_p5p8_b2_a,   leg_quad_p5p9_b2_a,   leg_quad_p5p10_b2_a,   leg_quad_p6p0_b2_a,   leg_quad_p6p1_b2_a,   leg_quad_p6p2_b2_a,   leg_quad_p6p3_b2_a,   leg_quad_p6p4_b2_a,   leg_quad_p6p5_b2_a,   leg_quad_p6p6_b2_a,   leg_quad_p6p7_b2_a,   leg_quad_p6p8_b2_a,   leg_quad_p6p9_b2_a,   leg_quad_p6p10_b2_a,   leg_quad_p7p0_b2_a,   leg_quad_p7p1_b2_a,   leg_quad_p7p2_b2_a,   leg_quad_p7p3_b2_a,   leg_quad_p7p4_b2_a,   leg_quad_p7p5_b2_a,   leg_quad_p7p6_b2_a,   leg_quad_p7p7_b2_a,   leg_quad_p7p8_b2_a,   leg_quad_p7p9_b2_a,   leg_quad_p7p10_b2_a,   leg_quad_p8p0_b2_a,   leg_quad_p8p1_b2_a,   leg_quad_p8p2_b2_a,   leg_quad_p8p3_b2_a,   leg_quad_p8p4_b2_a,   leg_quad_p8p5_b2_a,   leg_quad_p8p6_b2_a,   leg_quad_p8p7_b2_a,   leg_quad_p8p8_b2_a,   leg_quad_p8p9_b2_a,   leg_quad_p8p10_b2_a,   leg_quad_p9p0_b2_a,   leg_quad_p9p1_b2_a,   leg_quad_p9p2_b2_a,   leg_quad_p9p3_b2_a,   leg_quad_p9p4_b2_a,   leg_quad_p9p5_b2_a,   leg_quad_p9p6_b2_a,   leg_quad_p9p7_b2_a,   leg_quad_p9p8_b2_a,   leg_quad_p9p9_b2_a,   leg_quad_p9p10_b2_a,   leg_quad_p10p0_b2_a,   leg_quad_p10p1_b2_a,   leg_quad_p10p2_b2_a,   leg_quad_p10p3_b2_a,   leg_quad_p10p4_b2_a,   leg_quad_p10p5_b2_a,   leg_quad_p10p6_b2_a,   leg_quad_p10p7_b2_a,   leg_quad_p10p8_b2_a,   leg_quad_p10p9_b2_a,   leg_quad_p10p10_b2_a,   leg_quad_p11p0_b2_a,   leg_quad_p11p1_b2_a,   leg_quad_p11p2_b2_a,   leg_quad_p11p3_b2_a,   leg_quad_p11p4_b2_a,   leg_quad_p11p5_b2_a,   leg_quad_p11p6_b2_a,   leg_quad_p11p7_b2_a,   leg_quad_p11p8_b2_a,   leg_quad_p11p9_b2_a,   leg_quad_p11p10_b2_a, };

    static Shapeset::shape_fn_t leg_quad_fn_b[] =
    {
      leg_quad_p0_e1_b_0, leg_quad_p0_e1_b_1, leg_quad_p0_e2_b_0, leg_quad_p0_e2_b_1,  leg_quad_p0_e3_b, leg_quad_p0_e3_b, leg_quad_p0_e4_b, leg_quad_p0_e4_b,
      leg_quad_p1_e1_b, leg_quad_p1_e1_b, leg_quad_p1_e2_b, leg_quad_p1_e2_b, leg_quad_p1_e3_b, leg_quad_p1_e3_b, leg_quad_p1_e4_b, leg_quad_p1_e4_b,
      leg_quad_p2_e1_b_0, leg_quad_p2_e1_b_1, leg_quad_p2_e2_b_0, leg_quad_p2_e2_b_1, leg_quad_p2_e3_b, leg_quad_p2_e3_b, leg_quad_p2_e4_b, leg_quad_p2_e4_b,
      leg_quad_p3_e1_b, leg_quad_p3_e1_b, leg_quad_p3_e2_b, leg_quad_p3_e2_b, leg_quad_p3_e3_b, leg_quad_p3_e3_b, leg_quad_p3_e4_b, leg_quad_p3_e4_b,
      leg_quad_p4_e1_b_0, leg_quad_p4_e1_b_1, leg_quad_p4_e2_b_0, leg_quad_p4_e2_b_1, leg_quad_p4_e3_b, leg_quad_p4_e3_b, leg_quad_p4_e4_b, leg_quad_p4_e4_b,
      leg_quad_p5_e1_b, leg_quad_p5_e1_b, leg_quad_p5_e2_b, leg_quad_p5_e2_b, leg_quad_p5_e3_b, leg_quad_p5_e3_b, leg_quad_p5_e4_b, leg_quad_p5_e4_b,
      leg_quad_p6_e1_b_0, leg_quad_p6_e1_b_1, leg_quad_p6_e2_b_0, leg_quad_p6_e2_b_1, leg_quad_p6_e3_b, leg_quad_p6_e3_b, leg_quad_p6_e4_b, leg_quad_p6_e4_b,
      leg_quad_p7_e1_b, leg_quad_p7_e1_b, leg_quad_p7_e2_b, leg_quad_p7_e2_b, leg_quad_p7_e3_b, leg_quad_p7_e3_b, leg_quad_p7_e4_b, leg_quad_p7_e4_b,
      leg_quad_p8_e1_b_0, leg_quad_p8_e1_b_1, leg_quad_p8_e2_b_0, leg_quad_p8_e2_b_1, leg_quad_p8_e3_b, leg_quad_p8_e3_b, leg_quad_p8_e4_b, leg_quad_p8_e4_b,
      leg_quad_p9_e1_b, leg_quad_p9_e1_b, leg_quad_p9_e2_b, leg_quad_p9_e2_b, leg_quad_p9_e3_b, leg_quad_p9_e3_b, leg_quad_p9_e4_b, leg_quad_p9_e4_b,
      leg_quad_p10_e1_b_0, leg_quad_p10_e1_b_1, leg_quad_p10_e2_b_0, leg_quad_p10_e2_b_1, leg_quad_p10_e3_b, leg_quad_p10_e3_b, leg_quad_p10_e4_b, leg_quad_p10_e4_b,

      leg_quad_p0p2_b1_b,   leg_quad_p0p3_b1_b,   leg_quad_p0p4_b1_b,   leg_quad_p0p5_b1_b,   leg_quad_p0p6_b1_b,   leg_quad_p0p7_b1_b,   leg_quad_p0p8_b1_b,   leg_quad_p0p9_b1_b,   leg_quad_p0p10_b1_b,   leg_quad_p0p11_b1_b,   leg_quad_p1p2_b1_b,   leg_quad_p1p3_b1_b,   leg_quad_p1p4_b1_b,   leg_quad_p1p5_b1_b,   leg_quad_p1p6_b1_b,   leg_quad_p1p7_b1_b,   leg_quad_p1p8_b1_b,   leg_quad_p1p9_b1_b,   leg_quad_p1p10_b1_b,   leg_quad_p1p11_b1_b,   leg_quad_p2p2_b1_b,   leg_quad_p2p3_b1_b,   leg_quad_p2p4_b1_b,   leg_quad_p2p5_b1_b,   leg_quad_p2p6_b1_b,   leg_quad_p2p7_b1_b,   leg_quad_p2p8_b1_b,   leg_quad_p2p9_b1_b,   leg_quad_p2p10_b1_b,   leg_quad_p2p11_b1_b,   leg_quad_p3p2_b1_b,   leg_quad_p3p3_b1_b,   leg_quad_p3p4_b1_b,   leg_quad_p3p5_b1_b,   leg_quad_p3p6_b1_b,   leg_quad_p3p7_b1_b,   leg_quad_p3p8_b1_b,   leg_quad_p3p9_b1_b,   leg_quad_p3p10_b1_b,   leg_quad_p3p11_b1_b,   leg_quad_p4p2_b1_b,   leg_quad_p4p3_b1_b,   leg_quad_p4p4_b1_b,   leg_quad_p4p5_b1_b,   leg_quad_p4p6_b1_b,   leg_quad_p4p7_b1_b,   leg_quad_p4p8_b1_b,   leg_quad_p4p9_b1_b,   leg_quad_p4p10_b1_b,   leg_quad_p4p11_b1_b,   leg_quad_p5p2_b1_b,   leg_quad_p5p3_b1_b,   leg_quad_p5p4_b1_b,   leg_quad_p5p5_b1_b,   leg_quad_p5p6_b1_b,   leg_quad_p5p7_b1_b,   leg_quad_p5p8_b1_b,   leg_quad_p5p9_b1_b,   leg_quad_p5p10_b1_b,   leg_quad_p5p11_b1_b,   leg_quad_p6p2_b1_b,   leg_quad_p6p3_b1_b,   leg_quad_p6p4_b1_b,   leg_quad_p6p5_b1_b,   leg_quad_p6p6_b1_b,   leg_quad_p6p7_b1_b,   leg_quad_p6p8_b1_b,   leg_quad_p6p9_b1_b,   leg_quad_p6p10_b1_b,   leg_quad_p6p11_b1_b,   leg_quad_p7p2_b1_b,   leg_quad_p7p3_b1_b,   leg_quad_p7p4_b1_b,   leg_quad_p7p5_b1_b,   leg_quad_p7p6_b1_b,   leg_quad_p7p7_b1_b,   leg_quad_p7p8_b1_b,   leg_quad_p7p9_b1_b,   leg_quad_p7p10_b1_b,   leg_quad_p7p11_b1_b,   leg_quad_p8p2_b1_b,   leg_quad_p8p3_b1_b,   leg_quad_p8p4_b1_b,   leg_quad_p8p5_b1_b,   leg_quad_p8p6_b1_b,   leg_quad_p8p7_b1_b,   leg_quad_p8p8_b1_b,   leg_quad_p8p9_b1_b,   leg_quad_p8p10_b1_b,   leg_quad_p8p11_b1_b,   leg_quad_p9p2_b1_b,   leg_quad_p9p3_b1_b,   leg_quad_p9p4_b1_b,   leg_quad_p9p5_b1_b,   leg_quad_p9p6_b1_b,   leg_quad_p9p7_b1_b,   leg_quad_p9p8_b1_b,   leg_quad_p9p9_b1_b,   leg_quad_p9p10_b1_b,   leg_quad_p9p11_b1_b,   leg_quad_p10p2_b1_b,   leg_quad_p10p3_b1_b,   leg_quad_p10p4_b1_b,   leg_quad_p10p5_b1_b,   leg_quad_p10p6_b1_b,   leg_quad_p10p7_b1_b,   leg_quad_p10p8_b1_b,   leg_quad_p10p9_b1_b,   leg_quad_p10p10_b1_b,   leg_quad_p10p11_b1_b,   leg_quad_p2p0_b2_b,   leg_quad_p2p1_b2_b,   leg_quad_p2p2_b2_b,   leg_quad_p2p3_b2_b,   leg_quad_p2p4_b2_b,   leg_quad_p2p5_b2_b,   leg_quad_p2p6_b2_b,   leg_quad_p2p7_b2_b,   leg_quad_p2p8_b2_b,   leg_quad_p2p9_b2_b,   leg_quad_p2p10_b2_b,   leg_quad_p3p0_b2_b,   leg_quad_p3p1_b2_b,   leg_quad_p3p2_b2_b,   leg_quad_p3p3_b2_b,   leg_quad_p3p4_b2_b,   leg_quad_p3p5_b2_b,   leg_quad_p3p6_b2_b,   leg_quad_p3p7_b2_b,   leg_quad_p3p8_b2_b,   leg_quad_p3p9_b2_b,   leg_quad_p3p10_b2_b,   leg_quad_p4p0_b2_b,   leg_quad_p4p1_b2_b,   leg_quad_p4p2_b2_b,   leg_quad_p4p3_b2_b,   leg_quad_p4p4_b2_b,   leg_quad_p4p5_b2_b,   leg_quad_p4p6_b2_b,   leg_quad_p4p7_b2_b,   leg_quad_p4p8_b2_b,   leg_quad_p4p9_b2_b,   leg_quad_p4p10_b2_b,   leg_quad_p5p0_b2_b,   leg_quad_p5p1_b2_b,   leg_quad_p5p2_b2_b,   leg_quad_p5p3_b2_b,   leg_quad_p5p4_b2_b,   leg_quad_p5p5_b2_b,   leg_quad_p5p6_b2_b,   leg_quad_p5p7_b2_b,   leg_quad_p5p8_b2_b,   leg_quad_p5p9_b2_b,   leg_quad_p5p10_b2_b,   leg_quad_p6p0_b2_b,   leg_quad_p6p1_b2_b,   leg_quad_p6p2_b2_b,   leg_quad_p6p3_b2_b,   leg_quad_p6p4_b2_b,   leg_quad_p6p5_b2_b,   leg_quad_p6p6_b2_b,   leg_quad_p6p7_b2_b,   leg_quad_p6p8_b2_b,   leg_quad_p6p9_b2_b,   leg_quad_p6p10_b2_b,   leg_quad_p7p0_b2_b,   leg_quad_p7p1_b2_b,   leg_quad_p7p2_b2_b,   leg_quad_p7p3_b2_b,   leg_quad_p7p4_b2_b,   leg_quad_p7p5_b2_b,   leg_quad_p7p6_b2_b,   leg_quad_p7p7_b2_b,   leg_quad_p7p8_b2_b,   leg_quad_p7p9_b2_b,   leg_quad_p7p10_b2_b,   leg_quad_p8p0_b2_b,   leg_quad_p8p1_b2_b,   leg_quad_p8p2_b2_b,   leg_quad_p8p3_b2_b,   leg_quad_p8p4_b2_b,   leg_quad_p8p5_b2_b,   leg_quad_p8p6_b2_b,   leg_quad_p8p7_b2_b,   leg_quad_p8p8_b2_b,   leg_quad_p8p9_b2_b,   leg_quad_p8p10_b2_b,   leg_quad_p9p0_b2_b,   leg_quad_p9p1_b2_b,   leg_quad_p9p2_b2_b,   leg_quad_p9p3_b2_b,   leg_quad_p9p4_b2_b,   leg_quad_p9p5_b2_b,   leg_quad_p9p6_b2_b,   leg_quad_p9p7_b2_b,   leg_quad_p9p8_b2_b,   leg_quad_p9p9_b2_b,   leg_quad_p9p10_b2_b,   leg_quad_p10p0_b2_b,   leg_quad_p10p1_b2_b,   leg_quad_p10p2_b2_b,   leg_quad_p10p3_b2_b,   leg_quad_p10p4_b2_b,   leg_quad_p10p5_b2_b,   leg_quad_p10p6_b2_b,   leg_quad_p10p7_b2_b,   leg_quad_p10p8_b2_b,   leg_quad_p10p9_b2_b,   leg_quad_p10p10_b2_b,   leg_quad_p11p0_b2_b,   leg_quad_p11p1_b2_b,   leg_quad_p11p2_b2_b,   leg_quad_p11p3_b2_b,   leg_quad_p11p4_b2_b,   leg_quad_p11p5_b2_b,   leg_quad_p11p6_b2_b,   leg_quad_p11p7_b2_b,   leg_quad_p11p8_b2_b,   leg_quad_p11p9_b2_b,   leg_quad_p11p10_b2_b, };

    static Shapeset::shape_fn_t leg_quad_fn_ax[] =
    {
      leg_quad_p0_e1_ax, leg_quad_p0_e1_ax, leg_quad_p0_e2_ax, leg_quad_p0_e2_ax, leg_quad_p0_e3_ax_0, leg_quad_p0_e3_ax_1, leg_quad_p0_e4_ax_0, leg_quad_p0_e4_ax_1,
      leg_quad_p1_e1_ax, leg_quad_p1_e1_ax, leg_quad_p1_e2_ax, leg_quad_p1_e2_ax, leg_quad_p1_e3_ax, leg_quad_p1_e3_ax, leg_quad_p1_e4_ax, leg_quad_p1_e4_ax,
      leg_quad_p2_e1_ax, leg_quad_p2_e1_ax, leg_quad_p2_e2_ax, leg_quad_p2_e2_ax, leg_quad_p2_e3_ax_0, leg_quad_p2_e3_ax_1, leg_quad_p2_e4_ax_0, leg_quad_p2_e4_ax_1,
      leg_quad_p3_e1_ax, leg_quad_p3_e1_ax, leg_quad_p3_e2_ax, leg_quad_p3_e2_ax, leg_quad_p3_e3_ax, leg_quad_p3_e3_ax, leg_quad_p3_e4_ax, leg_quad_p3_e4_ax,
      leg_quad_p4_e1_ax, leg_quad_p4_e1_ax, leg_quad_p4_e2_ax, leg_quad_p4_e2_ax, leg_quad_p4_e3_ax_0, leg_quad_p4_e3_ax_1, leg_quad_p4_e4_ax_0, leg_quad_p4_e4_ax_1,
      leg_quad_p5_e1_ax, leg_quad_p5_e1_ax, leg_quad_p5_e2_ax, leg_quad_p5_e2_ax, leg_quad_p5_e3_ax, leg_quad_p5_e3_ax, leg_quad_p5_e4_ax, leg_quad_p5_e4_ax,
      leg_quad_p6_e1_ax, leg_quad_p6_e1_ax, leg_quad_p6_e2_ax, leg_quad_p6_e2_ax, leg_quad_p6_e3_ax_0, leg_quad_p6_e3_ax_1, leg_quad_p6_e4_ax_0, leg_quad_p6_e4_ax_1,
      leg_quad_p7_e1_ax, leg_quad_p7_e1_ax, leg_quad_p7_e2_ax, leg_quad_p7_e2_ax, leg_quad_p7_e3_ax, leg_quad_p7_e3_ax, leg_quad_p7_e4_ax, leg_quad_p7_e4_ax,
      leg_quad_p8_e1_ax, leg_quad_p8_e1_ax, leg_quad_p8_e2_ax, leg_quad_p8_e2_ax, leg_quad_p8_e3_ax_0, leg_quad_p8_e3_ax_1, leg_quad_p8_e4_ax_0, leg_quad_p8_e4_ax_1,
      leg_quad_p9_e1_ax, leg_quad_p9_e1_ax, leg_quad_p9_e2_ax, leg_quad_p9_e2_ax, leg_quad_p9_e3_ax, leg_quad_p9_e3_ax, leg_quad_p9_e4_ax, leg_quad_p9_e4_ax,
      leg_quad_p10_e1_ax, leg_quad_p10_e1_ax, leg_quad_p10_e2_ax, leg_quad_p10_e2_ax, leg_quad_p10_e3_ax_0, leg_quad_p10_e3_ax_1, leg_quad_p10_e4_ax_0, leg_quad_p10_e4_ax_1,

      leg_quad_p0p2_b1_ax,   leg_quad_p0p3_b1_ax,   leg_quad_p0p4_b1_ax,   leg_quad_p0p5_b1_ax,   leg_quad_p0p6_b1_ax,   leg_quad_p0p7_b1_ax,   leg_quad_p0p8_b1_ax,   leg_quad_p0p9_b1_ax,   leg_quad_p0p10_b1_ax,   leg_quad_p0p11_b1_ax,   leg_quad_p1p2_b1_ax,   leg_quad_p1p3_b1_ax,   leg_quad_p1p4_b1_ax,   leg_quad_p1p5_b1_ax,   leg_quad_p1p6_b1_ax,   leg_quad_p1p7_b1_ax,   leg_quad_p1p8_b1_ax,   leg_quad_p1p9_b1_ax,   leg_quad_p1p10_b1_ax,   leg_quad_p1p11_b1_ax,   leg_quad_p2p2_b1_ax,   leg_quad_p2p3_b1_ax,   leg_quad_p2p4_b1_ax,   leg_quad_p2p5_b1_ax,   leg_quad_p2p6_b1_ax,   leg_quad_p2p7_b1_ax,   leg_quad_p2p8_b1_ax,   leg_quad_p2p9_b1_ax,   leg_quad_p2p10_b1_ax,   leg_quad_p2p11_b1_ax,   leg_quad_p3p2_b1_ax,   leg_quad_p3p3_b1_ax,   leg_quad_p3p4_b1_ax,   leg_quad_p3p5_b1_ax,   leg_quad_p3p6_b1_ax,   leg_quad_p3p7_b1_ax,   leg_quad_p3p8_b1_ax,   leg_quad_p3p9_b1_ax,   leg_quad_p3p10_b1_ax,   leg_quad_p3p11_b1_ax,   leg_quad_p4p2_b1_ax,   leg_quad_p4p3_b1_ax,   leg_quad_p4p4_b1_ax,   leg_quad_p4p5_b1_ax,   leg_quad_p4p6_b1_ax,   leg_quad_p4p7_b1_ax,   leg_quad_p4p8_b1_ax,   leg_quad_p4p9_b1_ax,   leg_quad_p4p10_b1_ax,   leg_quad_p4p11_b1_ax,   leg_quad_p5p2_b1_ax,   leg_quad_p5p3_b1_ax,   leg_quad_p5p4_b1_ax,   leg_quad_p5p5_b1_ax,   leg_quad_p5p6_b1_ax,   leg_quad_p5p7_b1_ax,   leg_quad_p5p8_b1_ax,   leg_quad_p5p9_b1_ax,   leg_quad_p5p10_b1_ax,   leg_quad_p5p11_b1_ax,   leg_quad_p6p2_b1_ax,   leg_quad_p6p3_b1_ax,   leg_quad_p6p4_b1_ax,   leg_quad_p6p5_b1_ax,   leg_quad_p6p6_b1_ax,   leg_quad_p6p7_b1_ax,   leg_quad_p6p8_b1_ax,   leg_quad_p6p9_b1_ax,   leg_quad_p6p10_b1_ax,   leg_quad_p6p11_b1_ax,   leg_quad_p7p2_b1_ax,   leg_quad_p7p3_b1_ax,   leg_quad_p7p4_b1_ax,   leg_quad_p7p5_b1_ax,   leg_quad_p7p6_b1_ax,   leg_quad_p7p7_b1_ax,   leg_quad_p7p8_b1_ax,   leg_quad_p7p9_b1_ax,   leg_quad_p7p10_b1_ax,   leg_quad_p7p11_b1_ax,   leg_quad_p8p2_b1_ax,   leg_quad_p8p3_b1_ax,   leg_quad_p8p4_b1_ax,   leg_quad_p8p5_b1_ax,   leg_quad_p8p6_b1_ax,   leg_quad_p8p7_b1_ax,   leg_quad_p8p8_b1_ax,   leg_quad_p8p9_b1_ax,   leg_quad_p8p10_b1_ax,   leg_quad_p8p11_b1_ax,   leg_quad_p9p2_b1_ax,   leg_quad_p9p3_b1_ax,   leg_quad_p9p4_b1_ax,   leg_quad_p9p5_b1_ax,   leg_quad_p9p6_b1_ax,   leg_quad_p9p7_b1_ax,   leg_quad_p9p8_b1_ax,   leg_quad_p9p9_b1_ax,   leg_quad_p9p10_b1_ax,   leg_quad_p9p11_b1_ax,   leg_quad_p10p2_b1_ax,   leg_quad_p10p3_b1_ax,   leg_quad_p10p4_b1_ax,   leg_quad_p10p5_b1_ax,   leg_quad_p10p6_b1_ax,   leg_quad_p10p7_b1_ax,   leg_quad_p10p8_b1_ax,   leg_quad_p10p9_b1_ax,   leg_quad_p10p10_b1_ax,   leg_quad_p10p11_b1_ax,   leg_quad_p2p0_b2_ax,   leg_quad_p2p1_b2_ax,   leg_quad_p2p2_b2_ax,   leg_quad_p2p3_b2_ax,   leg_quad_p2p4_b2_ax,   leg_quad_p2p5_b2_ax,   leg_quad_p2p6_b2_ax,   leg_quad_p2p7_b2_ax,   leg_quad_p2p8_b2_ax,   leg_quad_p2p9_b2_ax,   leg_quad_p2p10_b2_ax,   leg_quad_p3p0_b2_ax,   leg_quad_p3p1_b2_ax,   leg_quad_p3p2_b2_ax,   leg_quad_p3p3_b2_ax,   leg_quad_p3p4_b2_ax,   leg_quad_p3p5_b2_ax,   leg_quad_p3p6_b2_ax,   leg_quad_p3p7_b2_ax,   leg_quad_p3p8_b2_ax,   leg_quad_p3p9_b2_ax,   leg_quad_p3p10_b2_ax,   leg_quad_p4p0_b2_ax,   leg_quad_p4p1_b2_ax,   leg_quad_p4p2_b2_ax,   leg_quad_p4p3_b2_ax,   leg_quad_p4p4_b2_ax,   leg_quad_p4p5_b2_ax,   leg_quad_p4p6_b2_ax,   leg_quad_p4p7_b2_ax,   leg_quad_p4p8_b2_ax,   leg_quad_p4p9_b2_ax,   leg_quad_p4p10_b2_ax,   leg_quad_p5p0_b2_ax,   leg_quad_p5p1_b2_ax,   leg_quad_p5p2_b2_ax,   leg_quad_p5p3_b2_ax,   leg_quad_p5p4_b2_ax,   leg_quad_p5p5_b2_ax,   leg_quad_p5p6_b2_ax,   leg_quad_p5p7_b2_ax,   leg_quad_p5p8_b2_ax,   leg_quad_p5p9_b2_ax,   leg_quad_p5p10_b2_ax,   leg_quad_p6p0_b2_ax,   leg_quad_p6p1_b2_ax,   leg_quad_p6p2_b2_ax,   leg_quad_p6p3_b2_ax,   leg_quad_p6p4_b2_ax,   leg_quad_p6p5_b2_ax,   leg_quad_p6p6_b2_ax,   leg_quad_p6p7_b2_ax,   leg_quad_p6p8_b2_ax,   leg_quad_p6p9_b2_ax,   leg_quad_p6p10_b2_ax,   leg_quad_p7p0_b2_ax,   leg_quad_p7p1_b2_ax,   leg_quad_p7p2_b2_ax,   leg_quad_p7p3_b2_ax,   leg_quad_p7p4_b2_ax,   leg_quad_p7p5_b2_ax,   leg_quad_p7p6_b2_ax,   leg_quad_p7p7_b2_ax,   leg_quad_p7p8_b2_ax,   leg_quad_p7p9_b2_ax,   leg_quad_p7p10_b2_ax,   leg_quad_p8p0_b2_ax,   leg_quad_p8p1_b2_ax,   leg_quad_p8p2_b2_ax,   leg_quad_p8p3_b2_ax,   leg_quad_p8p4_b2_ax,   leg_quad_p8p5_b2_ax,   leg_quad_p8p6_b2_ax,   leg_quad_p8p7_b2_ax,   leg_quad_p8p8_b2_ax,   leg_quad_p8p9_b2_ax,   leg_quad_p8p10_b2_ax,   leg_quad_p9p0_b2_ax,   leg_quad_p9p1_b2_ax,   leg_quad_p9p2_b2_ax,   leg_quad_p9p3_b2_ax,   leg_quad_p9p4_b2_ax,   leg_quad_p9p5_b2_ax,   leg_quad_p9p6_b2_ax,   leg_quad_p9p7_b2_ax,   leg_quad_p9p8_b2_ax,   leg_quad_p9p9_b2_ax,   leg_quad_p9p10_b2_ax,   leg_quad_p10p0_b2_ax,   leg_quad_p10p1_b2_ax,   leg_quad_p10p2_b2_ax,   leg_quad_p10p3_b2_ax,   leg_quad_p10p4_b2_ax,   leg_quad_p10p5_b2_ax,   leg_quad_p10p6_b2_ax,   leg_quad_p10p7_b2_ax,   leg_quad_p10p8_b2_ax,   leg_quad_p10p9_b2_ax,   leg_quad_p10p10_b2_ax,   leg_quad_p11p0_b2_ax,   leg_quad_p11p1_b2_ax,   leg_quad_p11p2_b2_ax,   leg_quad_p11p3_b2_ax,   leg_quad_p11p4_b2_ax,   leg_quad_p11p5_b2_ax,   leg_quad_p11p6_b2_ax,   leg_quad_p11p7_b2_ax,   leg_quad_p11p8_b2_ax,   leg_quad_p11p9_b2_ax,   leg_quad_p11p10_b2_ax, };

    static Shapeset::shape_fn_t leg_quad_fn_bx[] =
    {
      leg_quad_p0_e1_bx_0, leg_quad_p0_e1_bx_1, leg_quad_p0_e2_bx_0, leg_quad_p0_e2_bx_1,  leg_quad_p0_e3_bx, leg_quad_p0_e3_bx, leg_quad_p0_e4_bx, leg_quad_p0_e4_bx,
      leg_quad_p1_e1_bx, leg_quad_p1_e1_bx, leg_quad_p1_e2_bx, leg_quad_p1_e2_bx, leg_quad_p1_e3_bx, leg_quad_p1_e3_bx, leg_quad_p1_e4_bx, leg_quad_p1_e4_bx,
      leg_quad_p2_e1_bx_0, leg_quad_p2_e1_bx_1, leg_quad_p2_e2_bx_0, leg_quad_p2_e2_bx_1, leg_quad_p2_e3_bx, leg_quad_p2_e3_bx, leg_quad_p2_e4_bx, leg_quad_p2_e4_bx,
      leg_quad_p3_e1_bx, leg_quad_p3_e1_bx, leg_quad_p3_e2_bx, leg_quad_p3_e2_bx, leg_quad_p3_e3_bx, leg_quad_p3_e3_bx, leg_quad_p3_e4_bx, leg_quad_p3_e4_bx,
      leg_quad_p4_e1_bx_0, leg_quad_p4_e1_bx_1, leg_quad_p4_e2_bx_0, leg_quad_p4_e2_bx_1, leg_quad_p4_e3_bx, leg_quad_p4_e3_bx, leg_quad_p4_e4_bx, leg_quad_p4_e4_bx,
      leg_quad_p5_e1_bx, leg_quad_p5_e1_bx, leg_quad_p5_e2_bx, leg_quad_p5_e2_bx, leg_quad_p5_e3_bx, leg_quad_p5_e3_bx, leg_quad_p5_e4_bx, leg_quad_p5_e4_bx,
      leg_quad_p6_e1_bx_0, leg_quad_p6_e1_bx_1, leg_quad_p6_e2_bx_0, leg_quad_p6_e2_bx_1, leg_quad_p6_e3_bx, leg_quad_p6_e3_bx, leg_quad_p6_e4_bx, leg_quad_p6_e4_bx,
      leg_quad_p7_e1_bx, leg_quad_p7_e1_bx, leg_quad_p7_e2_bx, leg_quad_p7_e2_bx, leg_quad_p7_e3_bx, leg_quad_p7_e3_bx, leg_quad_p7_e4_bx, leg_quad_p7_e4_bx,
      leg_quad_p8_e1_bx_0, leg_quad_p8_e1_bx_1, leg_quad_p8_e2_bx_0, leg_quad_p8_e2_bx_1, leg_quad_p8_e3_bx, leg_quad_p8_e3_bx, leg_quad_p8_e4_bx, leg_quad_p8_e4_bx,
      leg_quad_p9_e1_bx, leg_quad_p9_e1_bx, leg_quad_p9_e2_bx, leg_quad_p9_e2_bx, leg_quad_p9_e3_bx, leg_quad_p9_e3_bx, leg_quad_p9_e4_bx, leg_quad_p9_e4_bx,
      leg_quad_p10_e1_bx_0, leg_quad_p10_e1_bx_1, leg_quad_p10_e2_bx_0, leg_quad_p10_e2_bx_1, leg_quad_p10_e3_bx, leg_quad_p10_e3_bx, leg_quad_p10_e4_bx, leg_quad_p10_e4_bx,

      leg_quad_p0p2_b1_bx,   leg_quad_p0p3_b1_bx,   leg_quad_p0p4_b1_bx,   leg_quad_p0p5_b1_bx,   leg_quad_p0p6_b1_bx,   leg_quad_p0p7_b1_bx,   leg_quad_p0p8_b1_bx,   leg_quad_p0p9_b1_bx,   leg_quad_p0p10_b1_bx,   leg_quad_p0p11_b1_bx,   leg_quad_p1p2_b1_bx,   leg_quad_p1p3_b1_bx,   leg_quad_p1p4_b1_bx,   leg_quad_p1p5_b1_bx,   leg_quad_p1p6_b1_bx,   leg_quad_p1p7_b1_bx,   leg_quad_p1p8_b1_bx,   leg_quad_p1p9_b1_bx,   leg_quad_p1p10_b1_bx,   leg_quad_p1p11_b1_bx,   leg_quad_p2p2_b1_bx,   leg_quad_p2p3_b1_bx,   leg_quad_p2p4_b1_bx,   leg_quad_p2p5_b1_bx,   leg_quad_p2p6_b1_bx,   leg_quad_p2p7_b1_bx,   leg_quad_p2p8_b1_bx,   leg_quad_p2p9_b1_bx,   leg_quad_p2p10_b1_bx,   leg_quad_p2p11_b1_bx,   leg_quad_p3p2_b1_bx,   leg_quad_p3p3_b1_bx,   leg_quad_p3p4_b1_bx,   leg_quad_p3p5_b1_bx,   leg_quad_p3p6_b1_bx,   leg_quad_p3p7_b1_bx,   leg_quad_p3p8_b1_bx,   leg_quad_p3p9_b1_bx,   leg_quad_p3p10_b1_bx,   leg_quad_p3p11_b1_bx,   leg_quad_p4p2_b1_bx,   leg_quad_p4p3_b1_bx,   leg_quad_p4p4_b1_bx,   leg_quad_p4p5_b1_bx,   leg_quad_p4p6_b1_bx,   leg_quad_p4p7_b1_bx,   leg_quad_p4p8_b1_bx,   leg_quad_p4p9_b1_bx,   leg_quad_p4p10_b1_bx,   leg_quad_p4p11_b1_bx,   leg_quad_p5p2_b1_bx,   leg_quad_p5p3_b1_bx,   leg_quad_p5p4_b1_bx,   leg_quad_p5p5_b1_bx,   leg_quad_p5p6_b1_bx,   leg_quad_p5p7_b1_bx,   leg_quad_p5p8_b1_bx,   leg_quad_p5p9_b1_bx,   leg_quad_p5p10_b1_bx,   leg_quad_p5p11_b1_bx,   leg_quad_p6p2_b1_bx,   leg_quad_p6p3_b1_bx,   leg_quad_p6p4_b1_bx,   leg_quad_p6p5_b1_bx,   leg_quad_p6p6_b1_bx,   leg_quad_p6p7_b1_bx,   leg_quad_p6p8_b1_bx,   leg_quad_p6p9_b1_bx,   leg_quad_p6p10_b1_bx,   leg_quad_p6p11_b1_bx,   leg_quad_p7p2_b1_bx,   leg_quad_p7p3_b1_bx,   leg_quad_p7p4_b1_bx,   leg_quad_p7p5_b1_bx,   leg_quad_p7p6_b1_bx,   leg_quad_p7p7_b1_bx,   leg_quad_p7p8_b1_bx,   leg_quad_p7p9_b1_bx,   leg_quad_p7p10_b1_bx,   leg_quad_p7p11_b1_bx,   leg_quad_p8p2_b1_bx,   leg_quad_p8p3_b1_bx,   leg_quad_p8p4_b1_bx,   leg_quad_p8p5_b1_bx,   leg_quad_p8p6_b1_bx,   leg_quad_p8p7_b1_bx,   leg_quad_p8p8_b1_bx,   leg_quad_p8p9_b1_bx,   leg_quad_p8p10_b1_bx,   leg_quad_p8p11_b1_bx,   leg_quad_p9p2_b1_bx,   leg_quad_p9p3_b1_bx,   leg_quad_p9p4_b1_bx,   leg_quad_p9p5_b1_bx,   leg_quad_p9p6_b1_bx,   leg_quad_p9p7_b1_bx,   leg_quad_p9p8_b1_bx,   leg_quad_p9p9_b1_bx,   leg_quad_p9p10_b1_bx,   leg_quad_p9p11_b1_bx,   leg_quad_p10p2_b1_bx,   leg_quad_p10p3_b1_bx,   leg_quad_p10p4_b1_bx,   leg_quad_p10p5_b1_bx,   leg_quad_p10p6_b1_bx,   leg_quad_p10p7_b1_bx,   leg_quad_p10p8_b1_bx,   leg_quad_p10p9_b1_bx,   leg_quad_p10p10_b1_bx,   leg_quad_p10p11_b1_bx,   leg_quad_p2p0_b2_bx,   leg_quad_p2p1_b2_bx,   leg_quad_p2p2_b2_bx,   leg_quad_p2p3_b2_bx,   leg_quad_p2p4_b2_bx,   leg_quad_p2p5_b2_bx,   leg_quad_p2p6_b2_bx,   leg_quad_p2p7_b2_bx,   leg_quad_p2p8_b2_bx,   leg_quad_p2p9_b2_bx,   leg_quad_p2p10_b2_bx,   leg_quad_p3p0_b2_bx,   leg_quad_p3p1_b2_bx,   leg_quad_p3p2_b2_bx,   leg_quad_p3p3_b2_bx,   leg_quad_p3p4_b2_bx,   leg_quad_p3p5_b2_bx,   leg_quad_p3p6_b2_bx,   leg_quad_p3p7_b2_bx,   leg_quad_p3p8_b2_bx,   leg_quad_p3p9_b2_bx,   leg_quad_p3p10_b2_bx,   leg_quad_p4p0_b2_bx,   leg_quad_p4p1_b2_bx,   leg_quad_p4p2_b2_bx,   leg_quad_p4p3_b2_bx,   leg_quad_p4p4_b2_bx,   leg_quad_p4p5_b2_bx,   leg_quad_p4p6_b2_bx,   leg_quad_p4p7_b2_bx,   leg_quad_p4p8_b2_bx,   leg_quad_p4p9_b2_bx,   leg_quad_p4p10_b2_bx,   leg_quad_p5p0_b2_bx,   leg_quad_p5p1_b2_bx,   leg_quad_p5p2_b2_bx,   leg_quad_p5p3_b2_bx,   leg_quad_p5p4_b2_bx,   leg_quad_p5p5_b2_bx,   leg_quad_p5p6_b2_bx,   leg_quad_p5p7_b2_bx,   leg_quad_p5p8_b2_bx,   leg_quad_p5p9_b2_bx,   leg_quad_p5p10_b2_bx,   leg_quad_p6p0_b2_bx,   leg_quad_p6p1_b2_bx,   leg_quad_p6p2_b2_bx,   leg_quad_p6p3_b2_bx,   leg_quad_p6p4_b2_bx,   leg_quad_p6p5_b2_bx,   leg_quad_p6p6_b2_bx,   leg_quad_p6p7_b2_bx,   leg_quad_p6p8_b2_bx,   leg_quad_p6p9_b2_bx,   leg_quad_p6p10_b2_bx,   leg_quad_p7p0_b2_bx,   leg_quad_p7p1_b2_bx,   leg_quad_p7p2_b2_bx,   leg_quad_p7p3_b2_bx,   leg_quad_p7p4_b2_bx,   leg_quad_p7p5_b2_bx,   leg_quad_p7p6_b2_bx,   leg_quad_p7p7_b2_bx,   leg_quad_p7p8_b2_bx,   leg_quad_p7p9_b2_bx,   leg_quad_p7p10_b2_bx,   leg_quad_p8p0_b2_bx,   leg_quad_p8p1_b2_bx,   leg_quad_p8p2_b2_bx,   leg_quad_p8p3_b2_bx,   leg_quad_p8p4_b2_bx,   leg_quad_p8p5_b2_bx,   leg_quad_p8p6_b2_bx,   leg_quad_p8p7_b2_bx,   leg_quad_p8p8_b2_bx,   leg_quad_p8p9_b2_bx,   leg_quad_p8p10_b2_bx,   leg_quad_p9p0_b2_bx,   leg_quad_p9p1_b2_bx,   leg_quad_p9p2_b2_bx,   leg_quad_p9p3_b2_bx,   leg_quad_p9p4_b2_bx,   leg_quad_p9p5_b2_bx,   leg_quad_p9p6_b2_bx,   leg_quad_p9p7_b2_bx,   leg_quad_p9p8_b2_bx,   leg_quad_p9p9_b2_bx,   leg_quad_p9p10_b2_bx,   leg_quad_p10p0_b2_bx,   leg_quad_p10p1_b2_bx,   leg_quad_p10p2_b2_bx,   leg_quad_p10p3_b2_bx,   leg_quad_p10p4_b2_bx,   leg_quad_p10p5_b2_bx,   leg_quad_p10p6_b2_bx,   leg_quad_p10p7_b2_bx,   leg_quad_p10p8_b2_bx,   leg_quad_p10p9_b2_bx,   leg_quad_p10p10_b2_bx,   leg_quad_p11p0_b2_bx,   leg_quad_p11p1_b2_bx,   leg_quad_p11p2_b2_bx,   leg_quad_p11p3_b2_bx,   leg_quad_p11p4_b2_bx,   leg_quad_p11p5_b2_bx,   leg_quad_p11p6_b2_bx,   leg_quad_p11p7_b2_bx,   leg_quad_p11p8_b2_bx,   leg_quad_p11p9_b2_bx,   leg_quad_p11p10_b2_bx, };

    static Shapeset::shape_fn_t leg_quad_fn_ay[] =
    {
      leg_quad_p0_e1_ay, leg_quad_p0_e1_ay, leg_quad_p0_e2_ay, leg_quad_p0_e2_ay, leg_quad_p0_e3_ay_0, leg_quad_p0_e3_ay_1, leg_quad_p0_e4_ay_0, leg_quad_p0_e4_ay_1,
      leg_quad_p1_e1_ay, leg_quad_p1_e1_ay, leg_quad_p1_e2_ay, leg_quad_p1_e2_ay, leg_quad_p1_e3_ay, leg_quad_p1_e3_ay, leg_quad_p1_e4_ay, leg_quad_p1_e4_ay,
      leg_quad_p2_e1_ay, leg_quad_p2_e1_ay, leg_quad_p2_e2_ay, leg_quad_p2_e2_ay, leg_quad_p2_e3_ay_0, leg_quad_p2_e3_ay_1, leg_quad_p2_e4_ay_0, leg_quad_p2_e4_ay_1,
      leg_quad_p3_e1_ay, leg_quad_p3_e1_ay, leg_quad_p3_e2_ay, leg_quad_p3_e2_ay, leg_quad_p3_e3_ay, leg_quad_p3_e3_ay, leg_quad_p3_e4_ay, leg_quad_p3_e4_ay,
      leg_quad_p4_e1_ay, leg_quad_p4_e1_ay, leg_quad_p4_e2_ay, leg_quad_p4_e2_ay, leg_quad_p4_e3_ay_0, leg_quad_p4_e3_ay_1, leg_quad_p4_e4_ay_0, leg_quad_p4_e4_ay_1,
      leg_quad_p5_e1_ay, leg_quad_p5_e1_ay, leg_quad_p5_e2_ay, leg_quad_p5_e2_ay, leg_quad_p5_e3_ay, leg_quad_p5_e3_ay, leg_quad_p5_e4_ay, leg_quad_p5_e4_ay,
      leg_quad_p6_e1_ay, leg_quad_p6_e1_ay, leg_quad_p6_e2_ay, leg_quad_p6_e2_ay, leg_quad_p6_e3_ay_0, leg_quad_p6_e3_ay_1, leg_quad_p6_e4_ay_0, leg_quad_p6_e4_ay_1,
      leg_quad_p7_e1_ay, leg_quad_p7_e1_ay, leg_quad_p7_e2_ay, leg_quad_p7_e2_ay, leg_quad_p7_e3_ay, leg_quad_p7_e3_ay, leg_quad_p7_e4_ay, leg_quad_p7_e4_ay,
      leg_quad_p8_e1_ay, leg_quad_p8_e1_ay, leg_quad_p8_e2_ay, leg_quad_p8_e2_ay, leg_quad_p8_e3_ay_0, leg_quad_p8_e3_ay_1, leg_quad_p8_e4_ay_0, leg_quad_p8_e4_ay_1,
      leg_quad_p9_e1_ay, leg_quad_p9_e1_ay, leg_quad_p9_e2_ay, leg_quad_p9_e2_ay, leg_quad_p9_e3_ay, leg_quad_p9_e3_ay, leg_quad_p9_e4_ay, leg_quad_p9_e4_ay,
      leg_quad_p10_e1_ay, leg_quad_p10_e1_ay, leg_quad_p10_e2_ay, leg_quad_p10_e2_ay, leg_quad_p10_e3_ay_0, leg_quad_p10_e3_ay_1, leg_quad_p10_e4_ay_0, leg_quad_p10_e4_ay_1,

      leg_quad_p0p2_b1_ay,   leg_quad_p0p3_b1_ay,   leg_quad_p0p4_b1_ay,   leg_quad_p0p5_b1_ay,   leg_quad_p0p6_b1_ay,   leg_quad_p0p7_b1_ay,   leg_quad_p0p8_b1_ay,   leg_quad_p0p9_b1_ay,   leg_quad_p0p10_b1_ay,   leg_quad_p0p11_b1_ay,   leg_quad_p1p2_b1_ay,   leg_quad_p1p3_b1_ay,   leg_quad_p1p4_b1_ay,   leg_quad_p1p5_b1_ay,   leg_quad_p1p6_b1_ay,   leg_quad_p1p7_b1_ay,   leg_quad_p1p8_b1_ay,   leg_quad_p1p9_b1_ay,   leg_quad_p1p10_b1_ay,   leg_quad_p1p11_b1_ay,   leg_quad_p2p2_b1_ay,   leg_quad_p2p3_b1_ay,   leg_quad_p2p4_b1_ay,   leg_quad_p2p5_b1_ay,   leg_quad_p2p6_b1_ay,   leg_quad_p2p7_b1_ay,   leg_quad_p2p8_b1_ay,   leg_quad_p2p9_b1_ay,   leg_quad_p2p10_b1_ay,   leg_quad_p2p11_b1_ay,   leg_quad_p3p2_b1_ay,   leg_quad_p3p3_b1_ay,   leg_quad_p3p4_b1_ay,   leg_quad_p3p5_b1_ay,   leg_quad_p3p6_b1_ay,   leg_quad_p3p7_b1_ay,   leg_quad_p3p8_b1_ay,   leg_quad_p3p9_b1_ay,   leg_quad_p3p10_b1_ay,   leg_quad_p3p11_b1_ay,   leg_quad_p4p2_b1_ay,   leg_quad_p4p3_b1_ay,   leg_quad_p4p4_b1_ay,   leg_quad_p4p5_b1_ay,   leg_quad_p4p6_b1_ay,   leg_quad_p4p7_b1_ay,   leg_quad_p4p8_b1_ay,   leg_quad_p4p9_b1_ay,   leg_quad_p4p10_b1_ay,   leg_quad_p4p11_b1_ay,   leg_quad_p5p2_b1_ay,   leg_quad_p5p3_b1_ay,   leg_quad_p5p4_b1_ay,   leg_quad_p5p5_b1_ay,   leg_quad_p5p6_b1_ay,   leg_quad_p5p7_b1_ay,   leg_quad_p5p8_b1_ay,   leg_quad_p5p9_b1_ay,   leg_quad_p5p10_b1_ay,   leg_quad_p5p11_b1_ay,   leg_quad_p6p2_b1_ay,   leg_quad_p6p3_b1_ay,   leg_quad_p6p4_b1_ay,   leg_quad_p6p5_b1_ay,   leg_quad_p6p6_b1_ay,   leg_quad_p6p7_b1_ay,   leg_quad_p6p8_b1_ay,   leg_quad_p6p9_b1_ay,   leg_quad_p6p10_b1_ay,   leg_quad_p6p11_b1_ay,   leg_quad_p7p2_b1_ay,   leg_quad_p7p3_b1_ay,   leg_quad_p7p4_b1_ay,   leg_quad_p7p5_b1_ay,   leg_quad_p7p6_b1_ay,   leg_quad_p7p7_b1_ay,   leg_quad_p7p8_b1_ay,   leg_quad_p7p9_b1_ay,   leg_quad_p7p10_b1_ay,   leg_quad_p7p11_b1_ay,   leg_quad_p8p2_b1_ay,   leg_quad_p8p3_b1_ay,   leg_quad_p8p4_b1_ay,   leg_quad_p8p5_b1_ay,   leg_quad_p8p6_b1_ay,   leg_quad_p8p7_b1_ay,   leg_quad_p8p8_b1_ay,   leg_quad_p8p9_b1_ay,   leg_quad_p8p10_b1_ay,   leg_quad_p8p11_b1_ay,   leg_quad_p9p2_b1_ay,   leg_quad_p9p3_b1_ay,   leg_quad_p9p4_b1_ay,   leg_quad_p9p5_b1_ay,   leg_quad_p9p6_b1_ay,   leg_quad_p9p7_b1_ay,   leg_quad_p9p8_b1_ay,   leg_quad_p9p9_b1_ay,   leg_quad_p9p10_b1_ay,   leg_quad_p9p11_b1_ay,   leg_quad_p10p2_b1_ay,   leg_quad_p10p3_b1_ay,   leg_quad_p10p4_b1_ay,   leg_quad_p10p5_b1_ay,   leg_quad_p10p6_b1_ay,   leg_quad_p10p7_b1_ay,   leg_quad_p10p8_b1_ay,   leg_quad_p10p9_b1_ay,   leg_quad_p10p10_b1_ay,   leg_quad_p10p11_b1_ay,   leg_quad_p2p0_b2_ay,   leg_quad_p2p1_b2_ay,   leg_quad_p2p2_b2_ay,   leg_quad_p2p3_b2_ay,   leg_quad_p2p4_b2_ay,   leg_quad_p2p5_b2_ay,   leg_quad_p2p6_b2_ay,   leg_quad_p2p7_b2_ay,   leg_quad_p2p8_b2_ay,   leg_quad_p2p9_b2_ay,   leg_quad_p2p10_b2_ay,   leg_quad_p3p0_b2_ay,   leg_quad_p3p1_b2_ay,   leg_quad_p3p2_b2_ay,   leg_quad_p3p3_b2_ay,   leg_quad_p3p4_b2_ay,   leg_quad_p3p5_b2_ay,   leg_quad_p3p6_b2_ay,   leg_quad_p3p7_b2_ay,   leg_quad_p3p8_b2_ay,   leg_quad_p3p9_b2_ay,   leg_quad_p3p10_b2_ay,   leg_quad_p4p0_b2_ay,   leg_quad_p4p1_b2_ay,   leg_quad_p4p2_b2_ay,   leg_quad_p4p3_b2_ay,   leg_quad_p4p4_b2_ay,   leg_quad_p4p5_b2_ay,   leg_quad_p4p6_b2_ay,   leg_quad_p4p7_b2_ay,   leg_quad_p4p8_b2_ay,   leg_quad_p4p9_b2_ay,   leg_quad_p4p10_b2_ay,   leg_quad_p5p0_b2_ay,   leg_quad_p5p1_b2_ay,   leg_quad_p5p2_b2_ay,   leg_quad_p5p3_b2_ay,   leg_quad_p5p4_b2_ay,   leg_quad_p5p5_b2_ay,   leg_quad_p5p6_b2_ay,   leg_quad_p5p7_b2_ay,   leg_quad_p5p8_b2_ay,   leg_quad_p5p9_b2_ay,   leg_quad_p5p10_b2_ay,   leg_quad_p6p0_b2_ay,   leg_quad_p6p1_b2_ay,   leg_quad_p6p2_b2_ay,   leg_quad_p6p3_b2_ay,   leg_quad_p6p4_b2_ay,   leg_quad_p6p5_b2_ay,   leg_quad_p6p6_b2_ay,   leg_quad_p6p7_b2_ay,   leg_quad_p6p8_b2_ay,   leg_quad_p6p9_b2_ay,   leg_quad_p6p10_b2_ay,   leg_quad_p7p0_b2_ay,   leg_quad_p7p1_b2_ay,   leg_quad_p7p2_b2_ay,   leg_quad_p7p3_b2_ay,   leg_quad_p7p4_b2_ay,   leg_quad_p7p5_b2_ay,   leg_quad_p7p6_b2_ay,   leg_quad_p7p7_b2_ay,   leg_quad_p7p8_b2_ay,   leg_quad_p7p9_b2_ay,   leg_quad_p7p10_b2_ay,   leg_quad_p8p0_b2_ay,   leg_quad_p8p1_b2_ay,   leg_quad_p8p2_b2_ay,   leg_quad_p8p3_b2_ay,   leg_quad_p8p4_b2_ay,   leg_quad_p8p5_b2_ay,   leg_quad_p8p6_b2_ay,   leg_quad_p8p7_b2_ay,   leg_quad_p8p8_b2_ay,   leg_quad_p8p9_b2_ay,   leg_quad_p8p10_b2_ay,   leg_quad_p9p0_b2_ay,   leg_quad_p9p1_b2_ay,   leg_quad_p9p2_b2_ay,   leg_quad_p9p3_b2_ay,   leg_quad_p9p4_b2_ay,   leg_quad_p9p5_b2_ay,   leg_quad_p9p6_b2_ay,   leg_quad_p9p7_b2_ay,   leg_quad_p9p8_b2_ay,   leg_quad_p9p9_b2_ay,   leg_quad_p9p10_b2_ay,   leg_quad_p10p0_b2_ay,   leg_quad_p10p1_b2_ay,   leg_quad_p10p2_b2_ay,   leg_quad_p10p3_b2_ay,   leg_quad_p10p4_b2_ay,   leg_quad_p10p5_b2_ay,   leg_quad_p10p6_b2_ay,   leg_quad_p10p7_b2_ay,   leg_quad_p10p8_b2_ay,   leg_quad_p10p9_b2_ay,   leg_quad_p10p10_b2_ay,   leg_quad_p11p0_b2_ay,   leg_quad_p11p1_b2_ay,   leg_quad_p11p2_b2_ay,   leg_quad_p11p3_b2_ay,   leg_quad_p11p4_b2_ay,   leg_quad_p11p5_b2_ay,   leg_quad_p11p6_b2_ay,   leg_quad_p11p7_b2_ay,   leg_quad_p11p8_b2_ay,   leg_quad_p11p9_b2_ay,   leg_quad_p11p10_b2_ay, };

    static Shapeset::shape_fn_t leg_quad_fn_by[] =
    {
      leg_quad_p0_e1_by_0, leg_quad_p0_e1_by_1, leg_quad_p0_e2_by_0, leg_quad_p0_e2_by_1,  leg_quad_p0_e3_by, leg_quad_p0_e3_by, leg_quad_p0_e4_by, leg_quad_p0_e4_by,
      leg_quad_p1_e1_by, leg_quad_p1_e1_by, leg_quad_p1_e2_by, leg_quad_p1_e2_by, leg_quad_p1_e3_by, leg_quad_p1_e3_by, leg_quad_p1_e4_by, leg_quad_p1_e4_b,
      leg_quad_p2_e1_by_0, leg_quad_p2_e1_by_1, leg_quad_p2_e2_by_0, leg_quad_p2_e2_by_1, leg_quad_p2_e3_by, leg_quad_p2_e3_by, leg_quad_p2_e4_by, leg_quad_p2_e4_by,
      leg_quad_p3_e1_by, leg_quad_p3_e1_by, leg_quad_p3_e2_by, leg_quad_p3_e2_by, leg_quad_p3_e3_by, leg_quad_p3_e3_by, leg_quad_p3_e4_by, leg_quad_p3_e4_b,
      leg_quad_p4_e1_by_0, leg_quad_p4_e1_by_1, leg_quad_p4_e2_by_0, leg_quad_p4_e2_by_1, leg_quad_p4_e3_by, leg_quad_p4_e3_by, leg_quad_p4_e4_by, leg_quad_p4_e4_by,
      leg_quad_p5_e1_by, leg_quad_p5_e1_by, leg_quad_p5_e2_by, leg_quad_p5_e2_by, leg_quad_p5_e3_by, leg_quad_p5_e3_by, leg_quad_p5_e4_by, leg_quad_p5_e4_b,
      leg_quad_p6_e1_by_0, leg_quad_p6_e1_by_1, leg_quad_p6_e2_by_0, leg_quad_p6_e2_by_1, leg_quad_p6_e3_by, leg_quad_p6_e3_by, leg_quad_p6_e4_by, leg_quad_p6_e4_by,
      leg_quad_p7_e1_by, leg_quad_p7_e1_by, leg_quad_p7_e2_by, leg_quad_p7_e2_by, leg_quad_p7_e3_by, leg_quad_p7_e3_by, leg_quad_p7_e4_by, leg_quad_p7_e4_b,
      leg_quad_p8_e1_by_0, leg_quad_p8_e1_by_1, leg_quad_p8_e2_by_0, leg_quad_p8_e2_by_1, leg_quad_p8_e3_by, leg_quad_p8_e3_by, leg_quad_p8_e4_by, leg_quad_p8_e4_by,
      leg_quad_p9_e1_by, leg_quad_p9_e1_by, leg_quad_p9_e2_by, leg_quad_p9_e2_by, leg_quad_p9_e3_by, leg_quad_p9_e3_by, leg_quad_p9_e4_by, leg_quad_p9_e4_b,
      leg_quad_p10_e1_by_0, leg_quad_p10_e1_by_1, leg_quad_p10_e2_by_0, leg_quad_p10_e2_by_1, leg_quad_p10_e3_by, leg_quad_p10_e3_by, leg_quad_p10_e4_by, leg_quad_p10_e4_by,

      leg_quad_p0p2_b1_by,   leg_quad_p0p3_b1_by,   leg_quad_p0p4_b1_by,   leg_quad_p0p5_b1_by,   leg_quad_p0p6_b1_by,   leg_quad_p0p7_b1_by,   leg_quad_p0p8_b1_by,   leg_quad_p0p9_b1_by,   leg_quad_p0p10_b1_by,   leg_quad_p0p11_b1_by,   leg_quad_p1p2_b1_by,   leg_quad_p1p3_b1_by,   leg_quad_p1p4_b1_by,   leg_quad_p1p5_b1_by,   leg_quad_p1p6_b1_by,   leg_quad_p1p7_b1_by,   leg_quad_p1p8_b1_by,   leg_quad_p1p9_b1_by,   leg_quad_p1p10_b1_by,   leg_quad_p1p11_b1_by,   leg_quad_p2p2_b1_by,   leg_quad_p2p3_b1_by,   leg_quad_p2p4_b1_by,   leg_quad_p2p5_b1_by,   leg_quad_p2p6_b1_by,   leg_quad_p2p7_b1_by,   leg_quad_p2p8_b1_by,   leg_quad_p2p9_b1_by,   leg_quad_p2p10_b1_by,   leg_quad_p2p11_b1_by,   leg_quad_p3p2_b1_by,   leg_quad_p3p3_b1_by,   leg_quad_p3p4_b1_by,   leg_quad_p3p5_b1_by,   leg_quad_p3p6_b1_by,   leg_quad_p3p7_b1_by,   leg_quad_p3p8_b1_by,   leg_quad_p3p9_b1_by,   leg_quad_p3p10_b1_by,   leg_quad_p3p11_b1_by,   leg_quad_p4p2_b1_by,   leg_quad_p4p3_b1_by,   leg_quad_p4p4_b1_by,   leg_quad_p4p5_b1_by,   leg_quad_p4p6_b1_by,   leg_quad_p4p7_b1_by,   leg_quad_p4p8_b1_by,   leg_quad_p4p9_b1_by,   leg_quad_p4p10_b1_by,   leg_quad_p4p11_b1_by,   leg_quad_p5p2_b1_by,   leg_quad_p5p3_b1_by,   leg_quad_p5p4_b1_by,   leg_quad_p5p5_b1_by,   leg_quad_p5p6_b1_by,   leg_quad_p5p7_b1_by,   leg_quad_p5p8_b1_by,   leg_quad_p5p9_b1_by,   leg_quad_p5p10_b1_by,   leg_quad_p5p11_b1_by,   leg_quad_p6p2_b1_by,   leg_quad_p6p3_b1_by,   leg_quad_p6p4_b1_by,   leg_quad_p6p5_b1_by,   leg_quad_p6p6_b1_by,   leg_quad_p6p7_b1_by,   leg_quad_p6p8_b1_by,   leg_quad_p6p9_b1_by,   leg_quad_p6p10_b1_by,   leg_quad_p6p11_b1_by,   leg_quad_p7p2_b1_by,   leg_quad_p7p3_b1_by,   leg_quad_p7p4_b1_by,   leg_quad_p7p5_b1_by,   leg_quad_p7p6_b1_by,   leg_quad_p7p7_b1_by,   leg_quad_p7p8_b1_by,   leg_quad_p7p9_b1_by,   leg_quad_p7p10_b1_by,   leg_quad_p7p11_b1_by,   leg_quad_p8p2_b1_by,   leg_quad_p8p3_b1_by,   leg_quad_p8p4_b1_by,   leg_quad_p8p5_b1_by,   leg_quad_p8p6_b1_by,   leg_quad_p8p7_b1_by,   leg_quad_p8p8_b1_by,   leg_quad_p8p9_b1_by,   leg_quad_p8p10_b1_by,   leg_quad_p8p11_b1_by,   leg_quad_p9p2_b1_by,   leg_quad_p9p3_b1_by,   leg_quad_p9p4_b1_by,   leg_quad_p9p5_b1_by,   leg_quad_p9p6_b1_by,   leg_quad_p9p7_b1_by,   leg_quad_p9p8_b1_by,   leg_quad_p9p9_b1_by,   leg_quad_p9p10_b1_by,   leg_quad_p9p11_b1_by,   leg_quad_p10p2_b1_by,   leg_quad_p10p3_b1_by,   leg_quad_p10p4_b1_by,   leg_quad_p10p5_b1_by,   leg_quad_p10p6_b1_by,   leg_quad_p10p7_b1_by,   leg_quad_p10p8_b1_by,   leg_quad_p10p9_b1_by,   leg_quad_p10p10_b1_by,   leg_quad_p10p11_b1_by,   leg_quad_p2p0_b2_by,   leg_quad_p2p1_b2_by,   leg_quad_p2p2_b2_by,   leg_quad_p2p3_b2_by,   leg_quad_p2p4_b2_by,   leg_quad_p2p5_b2_by,   leg_quad_p2p6_b2_by,   leg_quad_p2p7_b2_by,   leg_quad_p2p8_b2_by,   leg_quad_p2p9_b2_by,   leg_quad_p2p10_b2_by,   leg_quad_p3p0_b2_by,   leg_quad_p3p1_b2_by,   leg_quad_p3p2_b2_by,   leg_quad_p3p3_b2_by,   leg_quad_p3p4_b2_by,   leg_quad_p3p5_b2_by,   leg_quad_p3p6_b2_by,   leg_quad_p3p7_b2_by,   leg_quad_p3p8_b2_by,   leg_quad_p3p9_b2_by,   leg_quad_p3p10_b2_by,   leg_quad_p4p0_b2_by,   leg_quad_p4p1_b2_by,   leg_quad_p4p2_b2_by,   leg_quad_p4p3_b2_by,   leg_quad_p4p4_b2_by,   leg_quad_p4p5_b2_by,   leg_quad_p4p6_b2_by,   leg_quad_p4p7_b2_by,   leg_quad_p4p8_b2_by,   leg_quad_p4p9_b2_by,   leg_quad_p4p10_b2_by,   leg_quad_p5p0_b2_by,   leg_quad_p5p1_b2_by,   leg_quad_p5p2_b2_by,   leg_quad_p5p3_b2_by,   leg_quad_p5p4_b2_by,   leg_quad_p5p5_b2_by,   leg_quad_p5p6_b2_by,   leg_quad_p5p7_b2_by,   leg_quad_p5p8_b2_by,   leg_quad_p5p9_b2_by,   leg_quad_p5p10_b2_by,   leg_quad_p6p0_b2_by,   leg_quad_p6p1_b2_by,   leg_quad_p6p2_b2_by,   leg_quad_p6p3_b2_by,   leg_quad_p6p4_b2_by,   leg_quad_p6p5_b2_by,   leg_quad_p6p6_b2_by,   leg_quad_p6p7_b2_by,   leg_quad_p6p8_b2_by,   leg_quad_p6p9_b2_by,   leg_quad_p6p10_b2_by,   leg_quad_p7p0_b2_by,   leg_quad_p7p1_b2_by,   leg_quad_p7p2_b2_by,   leg_quad_p7p3_b2_by,   leg_quad_p7p4_b2_by,   leg_quad_p7p5_b2_by,   leg_quad_p7p6_b2_by,   leg_quad_p7p7_b2_by,   leg_quad_p7p8_b2_by,   leg_quad_p7p9_b2_by,   leg_quad_p7p10_b2_by,   leg_quad_p8p0_b2_by,   leg_quad_p8p1_b2_by,   leg_quad_p8p2_b2_by,   leg_quad_p8p3_b2_by,   leg_quad_p8p4_b2_by,   leg_quad_p8p5_b2_by,   leg_quad_p8p6_b2_by,   leg_quad_p8p7_b2_by,   leg_quad_p8p8_b2_by,   leg_quad_p8p9_b2_by,   leg_quad_p8p10_b2_by,   leg_quad_p9p0_b2_by,   leg_quad_p9p1_b2_by,   leg_quad_p9p2_b2_by,   leg_quad_p9p3_b2_by,   leg_quad_p9p4_b2_by,   leg_quad_p9p5_b2_by,   leg_quad_p9p6_b2_by,   leg_quad_p9p7_b2_by,   leg_quad_p9p8_b2_by,   leg_quad_p9p9_b2_by,   leg_quad_p9p10_b2_by,   leg_quad_p10p0_b2_by,   leg_quad_p10p1_b2_by,   leg_quad_p10p2_b2_by,   leg_quad_p10p3_b2_by,   leg_quad_p10p4_b2_by,   leg_quad_p10p5_b2_by,   leg_quad_p10p6_b2_by,   leg_quad_p10p7_b2_by,   leg_quad_p10p8_b2_by,   leg_quad_p10p9_b2_by,   leg_quad_p10p10_b2_by,   leg_quad_p11p0_b2_by,   leg_quad_p11p1_b2_by,   leg_quad_p11p2_b2_by,   leg_quad_p11p3_b2_by,   leg_quad_p11p4_b2_by,   leg_quad_p11p5_b2_by,   leg_quad_p11p6_b2_by,   leg_quad_p11p7_b2_by,   leg_quad_p11p8_b2_by,   leg_quad_p11p9_b2_by,   leg_quad_p11p10_b2_by, };

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

    static int* leg_quad_bubble_indices[] =
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

    static int leg_quad_bubble_count[] =
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

    static int leg_quad_vertex_indices[4] ={-1, -1, -1, -1};

    static int leg_quad_edge_indices_0[] = { 4, 5, 12, 13, 20, 21, 28, 29, 36, 37, 44, 45, 52, 53, 60, 61, 68, 69, 76, 77, 84, 85, };

    static int leg_quad_edge_indices_1[] = { 2, 3, 10, 11, 18, 19, 26, 27, 34, 35, 42, 43, 50, 51, 58, 59, 66, 67, 74, 75, 82, 83, };

    static int leg_quad_edge_indices_2[] = { 6, 7, 14, 15, 22, 23, 30, 31, 38, 39, 46, 47, 54, 55, 62, 63, 70, 71, 78, 79, 86, 87, };

    static int leg_quad_edge_indices_3[] = { 0, 1, 8, 9, 16, 17, 24, 25, 32, 33, 40, 41, 48, 49, 56, 57, 64, 65, 72, 73, 80, 81, };

    static int* leg_quad_edge_indices[4] =
    {
      leg_quad_edge_indices_0,
      leg_quad_edge_indices_1,
      leg_quad_edge_indices_2,
      leg_quad_edge_indices_3,
    };

    #define oo H2D_MAKE_QUAD_ORDER

    static int leg_quad_index_to_order[] =
    {
     oo(1, 0), oo(1, 0),  oo(1, 0), oo(1, 0), oo(0, 1), oo(0, 1), oo(0, 1), oo(0, 1),
     oo(1, 1), oo(1, 1),  oo(1, 1), oo(1, 1), oo(1, 1), oo(1, 1), oo(1, 1), oo(1, 1),
     oo(1, 2), oo(1, 2),  oo(1, 2), oo(1, 2), oo(2, 1), oo(2, 1), oo(2, 1), oo(2, 1),
     oo(1, 3), oo(1, 3),  oo(1, 3), oo(1, 3), oo(3, 1), oo(3, 1), oo(3, 1), oo(3, 1),
     oo(1, 4), oo(1, 4),  oo(1, 4), oo(1, 4), oo(4, 1), oo(4, 1), oo(4, 1), oo(4, 1),
     oo(1, 5), oo(1, 5),  oo(1, 5), oo(1, 5), oo(5, 1), oo(5, 1), oo(5, 1), oo(5, 1),
     oo(1, 6), oo(1, 6),  oo(1, 6), oo(1, 6), oo(6, 1), oo(6, 1), oo(6, 1), oo(6, 1),
     oo(1, 7), oo(1, 7),  oo(1, 7), oo(1, 7), oo(7, 1), oo(7, 1), oo(7, 1), oo(7, 1),
     oo(1, 8), oo(1, 8),  oo(1, 8), oo(1, 8), oo(8, 1), oo(8, 1), oo(8, 1), oo(8, 1),
     oo(1, 9), oo(1, 9),  oo(1, 9), oo(1, 9), oo(9, 1), oo(9, 1), oo(9, 1), oo(9, 1),
     oo(1, 10), oo(1, 10),  oo(1, 10), oo(1, 10), oo(10, 1), oo(10, 1), oo(10, 1), oo(10, 1),

      oo(0, 2),  oo(0, 3),  oo(0, 4),  oo(0, 5),  oo(0, 6),  oo(0, 7),  oo(0, 8),  oo(0, 9),  oo(0, 10),  oo(0, 11),  oo(1, 2),  oo(1, 3),  oo(1, 4),  oo(1, 5),  oo(1, 6),  oo(1, 7),  oo(1, 8),  oo(1, 9),  oo(1, 10),  oo(1, 11),  oo(2, 2),  oo(2, 3),  oo(2, 4),  oo(2, 5),  oo(2, 6),  oo(2, 7),  oo(2, 8),  oo(2, 9),  oo(2, 10),  oo(2, 11),  oo(3, 2),  oo(3, 3),  oo(3, 4),  oo(3, 5),  oo(3, 6),  oo(3, 7),  oo(3, 8),  oo(3, 9),  oo(3, 10),  oo(3, 11),  oo(4, 2),  oo(4, 3),  oo(4, 4),  oo(4, 5),  oo(4, 6),  oo(4, 7),  oo(4, 8),  oo(4, 9),  oo(4, 10),  oo(4, 11),  oo(5, 2),  oo(5, 3),  oo(5, 4),  oo(5, 5),  oo(5, 6),  oo(5, 7),  oo(5, 8),  oo(5, 9),  oo(5, 10),  oo(5, 11),  oo(6, 2),  oo(6, 3),  oo(6, 4),  oo(6, 5),  oo(6, 6),  oo(6, 7),  oo(6, 8),  oo(6, 9),  oo(6, 10),  oo(6, 11),  oo(7, 2),  oo(7, 3),  oo(7, 4),  oo(7, 5),  oo(7, 6),  oo(7, 7),  oo(7, 8),  oo(7, 9),  oo(7, 10),  oo(7, 11),  oo(8, 2),  oo(8, 3),  oo(8, 4),  oo(8, 5),  oo(8, 6),  oo(8, 7),  oo(8, 8),  oo(8, 9),  oo(8, 10),  oo(8, 11),  oo(9, 2),  oo(9, 3),  oo(9, 4),  oo(9, 5),  oo(9, 6),  oo(9, 7),  oo(9, 8),  oo(9, 9),  oo(9, 10),  oo(9, 11),  oo(10, 2),  oo(10, 3),  oo(10, 4),  oo(10, 5),  oo(10, 6),  oo(10, 7),  oo(10, 8),  oo(10, 9),  oo(10, 10),  oo(10, 11),  oo(2, 0),  oo(2, 1),  oo(2, 2),  oo(2, 3),  oo(2, 4),  oo(2, 5),  oo(2, 6),  oo(2, 7),  oo(2, 8),  oo(2, 9),  oo(2, 10),  oo(3, 0),  oo(3, 1),  oo(3, 2),  oo(3, 3),  oo(3, 4),  oo(3, 5),  oo(3, 6),  oo(3, 7),  oo(3, 8),  oo(3, 9),  oo(3, 10),  oo(4, 0),  oo(4, 1),  oo(4, 2),  oo(4, 3),  oo(4, 4),  oo(4, 5),  oo(4, 6),  oo(4, 7),  oo(4, 8),  oo(4, 9),  oo(4, 10),  oo(5, 0),  oo(5, 1),  oo(5, 2),  oo(5, 3),  oo(5, 4),  oo(5, 5),  oo(5, 6),  oo(5, 7),  oo(5, 8),  oo(5, 9),  oo(5, 10),  oo(6, 0),  oo(6, 1),  oo(6, 2),  oo(6, 3),  oo(6, 4),  oo(6, 5),  oo(6, 6),  oo(6, 7),  oo(6, 8),  oo(6, 9),  oo(6, 10),  oo(7, 0),  oo(7, 1),  oo(7, 2),  oo(7, 3),  oo(7, 4),  oo(7, 5),  oo(7, 6),  oo(7, 7),  oo(7, 8),  oo(7, 9),  oo(7, 10),  oo(8, 0),  oo(8, 1),  oo(8, 2),  oo(8, 3),  oo(8, 4),  oo(8, 5),  oo(8, 6),  oo(8, 7),  oo(8, 8),  oo(8, 9),  oo(8, 10),  oo(9, 0),  oo(9, 1),  oo(9, 2),  oo(9, 3),  oo(9, 4),  oo(9, 5),  oo(9, 6),  oo(9, 7),  oo(9, 8),  oo(9, 9),  oo(9, 10),  oo(10, 0),  oo(10, 1),  oo(10, 2),  oo(10, 3),  oo(10, 4),  oo(10, 5),  oo(10, 6),  oo(10, 7),  oo(10, 8),  oo(10, 9),  oo(10, 10),  oo(11, 0),  oo(11, 1),  oo(11, 2),  oo(11, 3),  oo(11, 4),  oo(11, 5),  oo(11, 6),  oo(11, 7),  oo(11, 8),  oo(11, 9),  oo(11, 10),
    };

    static Shapeset::shape_fn_t* leg_quad_shape_fn_table[2] =
    {
      leg_quad_fn_a,
      leg_quad_fn_b
    };

    static Shapeset::shape_fn_t* leg_quad_shape_fn_table_x[2] =
    {
      leg_quad_fn_ax,
      leg_quad_fn_bx
    };

    static Shapeset::shape_fn_t* leg_quad_shape_fn_table_y[2] =
    {
      leg_quad_fn_ay,
      leg_quad_fn_by
    };

    //// triangle tables and class constructor ///////////////////////////////////////////////

    static Shapeset::shape_fn_t** leg_shape_fn_table[2] =
    {
      leg_tri_shape_fn_table,
      leg_quad_shape_fn_table
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_x[2] =
    {
      leg_tri_shape_fn_table_x,
      leg_quad_shape_fn_table_x
    };

    static Shapeset::shape_fn_t** leg_shape_fn_table_y[2] =
    {
      leg_tri_shape_fn_table_y,
      leg_quad_shape_fn_table_y
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

    void check_leg_tri(Shapeset* shapeset)
    {
      for (int i = 1; i <= 10; i++)
      {
        int nb = shapeset->get_num_bubbles(i, HERMES_MODE_TRIANGLE);
        if(nb != 3*(i-1) + (i-1)*(i-2))
          throw Hermes::Exceptions::Exception("Wrong bubble count");
      }

      int size_a  = sizeof(leg_tri_fn_a)  / sizeof(Shapeset::shape_fn_t);
      int size_b  = sizeof(leg_tri_fn_b)  / sizeof(Shapeset::shape_fn_t);
      int size_ax = sizeof(leg_tri_fn_ax) / sizeof(Shapeset::shape_fn_t);
      int size_bx = sizeof(leg_tri_fn_bx) / sizeof(Shapeset::shape_fn_t);
      int size_ay = sizeof(leg_tri_fn_ay) / sizeof(Shapeset::shape_fn_t);
      int size_by = sizeof(leg_tri_fn_by) / sizeof(Shapeset::shape_fn_t);

      if(size_a != size_b || size_a != size_ax || size_a != size_bx || size_a != size_ay || size_a != size_by)
        throw Hermes::Exceptions::Exception("Function tables dont have equal length.");

      if(size_a != leg_tri_bubble_indices[10][leg_tri_bubble_count[10]-1] + 1)
        throw Hermes::Exceptions::Exception("Bad index of last bubble");
    }

    HcurlShapesetLegendre::HcurlShapesetLegendre()
    {
      shape_table[0] = leg_shape_fn_table;
      shape_table[1] = leg_shape_fn_table_x;
      shape_table[2] = leg_shape_fn_table_y;
      shape_table[3] = NULL;
      shape_table[4] = NULL;
      shape_table[5] = NULL;

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
      num_components = 2;

      max_index[0] = 149;
      max_index[1] = 307;

      ebias = 0;

      comb_table = NULL;

      check_leg_tri(this);
    }
  }
}