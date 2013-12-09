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
#include "shapeset_hd_all.h"
namespace Hermes
{
  namespace Hermes2D
  {
    /* EDGE FUNCTIONS */

    /* ORDER 0 */

    static double hdiv_leg_quad_p0_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e1_a_1(double x, double y)
    {
      return -l0(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p0_e1_a_0(double x, double y)
    {
      return (l0(x) * Legendre0(y));
    }

    static double hdiv_leg_quad_p0_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e1_ax_1(double x, double y)
    {
      return -dl0(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p0_e1_ax_0(double x, double y)
    {
      return (dl0(x) * Legendre0(y));
    }

    static double hdiv_leg_quad_p0_e1_ay_1(double x, double y)
    {
      return -l0(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p0_e1_ay_0(double x, double y)
    {
      return (l0(x) * Legendre0x(y));
    }

    static double hdiv_leg_quad_p0_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e2_a_0(double x, double y)
    {
      return -l1(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p0_e2_a_1(double x, double y)
    {
      return (l1(x) * Legendre0(y));
    }

    static double hdiv_leg_quad_p0_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e2_ax_0(double x, double y)
    {
      return -dl1(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p0_e2_ax_1(double x, double y)
    {
      return (dl1(x) * Legendre0(y));
    }

    static double hdiv_leg_quad_p0_e2_ay_0(double x, double y)
    {
      return -l1(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p0_e2_ay_1(double x, double y)
    {
      return (l1(x) * Legendre0x(y));
    }

    static double hdiv_leg_quad_p0_e3_b_0(double x, double y)
    {
      return Legendre0(x) * l0(y);
    }

    static double hdiv_leg_quad_p0_e3_b_1(double x, double y)
    {
      return -(Legendre0(x) * l0(y));
    }

    static double hdiv_leg_quad_p0_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e3_bx_0(double x, double y)
    {
      return Legendre0x(x) * l0(y);
    }

    static double hdiv_leg_quad_p0_e3_bx_1(double x, double y)
    {
      return -(Legendre0x(x) * l0(y));
    }

    static double hdiv_leg_quad_p0_e3_by_0(double x, double y)
    {
      return Legendre0(x) * dl0(y);
    }

    static double hdiv_leg_quad_p0_e3_by_1(double x, double y)
    {
      return -(Legendre0(x) * dl0(y));
    }

    static double hdiv_leg_quad_p0_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e4_b_1(double x, double y)
    {
      return Legendre0(x) * l1(y);
    }

    static double hdiv_leg_quad_p0_e4_b_0(double x, double y)
    {
      return -(Legendre0(x) * l1(y));
    }

    static double hdiv_leg_quad_p0_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e4_bx_1(double x, double y)
    {
      return Legendre0x(x) * l1(y);
    }

    static double hdiv_leg_quad_p0_e4_bx_0(double x, double y)
    {
      return -(Legendre0x(x) * l1(y));
    }

    static double hdiv_leg_quad_p0_e4_by_1(double x, double y)
    {
      return Legendre0(x) * dl1(y);
    }

    static double hdiv_leg_quad_p0_e4_by_0(double x, double y)
    {
      return -(Legendre0(x) * dl1(y));
    }

    static double hdiv_leg_quad_p0_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 1 */

    static double hdiv_leg_quad_p1_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e1_a(double x, double y)
    {
      return -l0(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p1_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e1_ax(double x, double y)
    {
      return -dl0(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p1_e1_ay(double x, double y)
    {
      return -l0(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p1_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e2_a(double x, double y)
    {
      return -l1(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p1_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e2_ax(double x, double y)
    {
      return -dl1(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p1_e2_ay(double x, double y)
    {
      return -l1(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p1_e3_b(double x, double y)
    {
      return Legendre1(x) * l0(y);
    }

    static double hdiv_leg_quad_p1_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e3_bx(double x, double y)
    {
      return Legendre1x(x) * l0(y);
    }

    static double hdiv_leg_quad_p1_e3_by(double x, double y)
    {
      return Legendre1(x) * dl0(y);
    }

    static double hdiv_leg_quad_p1_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e4_b(double x, double y)
    {
      return Legendre1(x) * l1(y);
    }

    static double hdiv_leg_quad_p1_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e4_bx(double x, double y)
    {
      return Legendre1x(x) * l1(y);
    }

    static double hdiv_leg_quad_p1_e4_by(double x, double y)
    {
      return Legendre1(x) * dl1(y);
    }

    static double hdiv_leg_quad_p1_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 2 */

    static double hdiv_leg_quad_p2_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e1_a_1(double x, double y)
    {
      return -l0(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p2_e1_a_0(double x, double y)
    {
      return (l0(x) * Legendre2(y));
    }

    static double hdiv_leg_quad_p2_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e1_ax_1(double x, double y)
    {
      return -dl0(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p2_e1_ax_0(double x, double y)
    {
      return (dl0(x) * Legendre2(y));
    }

    static double hdiv_leg_quad_p2_e1_ay_1(double x, double y)
    {
      return -l0(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p2_e1_ay_0(double x, double y)
    {
      return (l0(x) * Legendre2x(y));
    }

    static double hdiv_leg_quad_p2_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e2_a_0(double x, double y)
    {
      return -l1(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p2_e2_a_1(double x, double y)
    {
      return (l1(x) * Legendre2(y));
    }

    static double hdiv_leg_quad_p2_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e2_ax_0(double x, double y)
    {
      return -dl1(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p2_e2_ax_1(double x, double y)
    {
      return (dl1(x) * Legendre2(y));
    }

    static double hdiv_leg_quad_p2_e2_ay_0(double x, double y)
    {
      return -l1(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p2_e2_ay_1(double x, double y)
    {
      return (l1(x) * Legendre2x(y));
    }

    static double hdiv_leg_quad_p2_e3_b_0(double x, double y)
    {
      return Legendre2(x) * l0(y);
    }

    static double hdiv_leg_quad_p2_e3_b_1(double x, double y)
    {
      return -(Legendre2(x) * l0(y));
    }

    static double hdiv_leg_quad_p2_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e3_bx_0(double x, double y)
    {
      return Legendre2x(x) * l0(y);
    }

    static double hdiv_leg_quad_p2_e3_bx_1(double x, double y)
    {
      return -(Legendre2x(x) * l0(y));
    }

    static double hdiv_leg_quad_p2_e3_by_0(double x, double y)
    {
      return Legendre2(x) * dl0(y);
    }

    static double hdiv_leg_quad_p2_e3_by_1(double x, double y)
    {
      return -(Legendre2(x) * dl0(y));
    }

    static double hdiv_leg_quad_p2_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e4_b_1(double x, double y)
    {
      return Legendre2(x) * l1(y);
    }

    static double hdiv_leg_quad_p2_e4_b_0(double x, double y)
    {
      return -(Legendre2(x) * l1(y));
    }

    static double hdiv_leg_quad_p2_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e4_bx_1(double x, double y)
    {
      return Legendre2x(x) * l1(y);
    }

    static double hdiv_leg_quad_p2_e4_bx_0(double x, double y)
    {
      return -(Legendre2x(x) * l1(y));
    }

    static double hdiv_leg_quad_p2_e4_by_1(double x, double y)
    {
      return Legendre2(x) * dl1(y);
    }

    static double hdiv_leg_quad_p2_e4_by_0(double x, double y)
    {
      return -(Legendre2(x) * dl1(y));
    }

    static double hdiv_leg_quad_p2_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 3 */

    static double hdiv_leg_quad_p3_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e1_a(double x, double y)
    {
      return -l0(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p3_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e1_ax(double x, double y)
    {
      return -dl0(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p3_e1_ay(double x, double y)
    {
      return -l0(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p3_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e2_a(double x, double y)
    {
      return -l1(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p3_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e2_ax(double x, double y)
    {
      return -dl1(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p3_e2_ay(double x, double y)
    {
      return -l1(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p3_e3_b(double x, double y)
    {
      return Legendre3(x) * l0(y);
    }

    static double hdiv_leg_quad_p3_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e3_bx(double x, double y)
    {
      return Legendre3x(x) * l0(y);
    }

    static double hdiv_leg_quad_p3_e3_by(double x, double y)
    {
      return Legendre3(x) * dl0(y);
    }

    static double hdiv_leg_quad_p3_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e4_b(double x, double y)
    {
      return Legendre3(x) * l1(y);
    }

    static double hdiv_leg_quad_p3_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e4_bx(double x, double y)
    {
      return Legendre3x(x) * l1(y);
    }

    static double hdiv_leg_quad_p3_e4_by(double x, double y)
    {
      return Legendre3(x) * dl1(y);
    }

    static double hdiv_leg_quad_p3_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 4 */

    static double hdiv_leg_quad_p4_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e1_a_1(double x, double y)
    {
      return -l0(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p4_e1_a_0(double x, double y)
    {
      return (l0(x) * Legendre4(y));
    }

    static double hdiv_leg_quad_p4_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e1_ax_1(double x, double y)
    {
      return -dl0(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p4_e1_ax_0(double x, double y)
    {
      return (dl0(x) * Legendre4(y));
    }

    static double hdiv_leg_quad_p4_e1_ay_1(double x, double y)
    {
      return -l0(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p4_e1_ay_0(double x, double y)
    {
      return (l0(x) * Legendre4x(y));
    }

    static double hdiv_leg_quad_p4_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e2_a_0(double x, double y)
    {
      return -l1(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p4_e2_a_1(double x, double y)
    {
      return (l1(x) * Legendre4(y));
    }

    static double hdiv_leg_quad_p4_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e2_ax_0(double x, double y)
    {
      return -dl1(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p4_e2_ax_1(double x, double y)
    {
      return (dl1(x) * Legendre4(y));
    }

    static double hdiv_leg_quad_p4_e2_ay_0(double x, double y)
    {
      return -l1(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p4_e2_ay_1(double x, double y)
    {
      return (l1(x) * Legendre4x(y));
    }

    static double hdiv_leg_quad_p4_e3_b_0(double x, double y)
    {
      return Legendre4(x) * l0(y);
    }

    static double hdiv_leg_quad_p4_e3_b_1(double x, double y)
    {
      return -(Legendre4(x) * l0(y));
    }

    static double hdiv_leg_quad_p4_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e3_bx_0(double x, double y)
    {
      return Legendre4x(x) * l0(y);
    }

    static double hdiv_leg_quad_p4_e3_bx_1(double x, double y)
    {
      return -(Legendre4x(x) * l0(y));
    }

    static double hdiv_leg_quad_p4_e3_by_0(double x, double y)
    {
      return Legendre4(x) * dl0(y);
    }

    static double hdiv_leg_quad_p4_e3_by_1(double x, double y)
    {
      return -(Legendre4(x) * dl0(y));
    }

    static double hdiv_leg_quad_p4_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e4_b_1(double x, double y)
    {
      return Legendre4(x) * l1(y);
    }

    static double hdiv_leg_quad_p4_e4_b_0(double x, double y)
    {
      return -(Legendre4(x) * l1(y));
    }

    static double hdiv_leg_quad_p4_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e4_bx_1(double x, double y)
    {
      return Legendre4x(x) * l1(y);
    }

    static double hdiv_leg_quad_p4_e4_bx_0(double x, double y)
    {
      return -(Legendre4x(x) * l1(y));
    }

    static double hdiv_leg_quad_p4_e4_by_1(double x, double y)
    {
      return Legendre4(x) * dl1(y);
    }

    static double hdiv_leg_quad_p4_e4_by_0(double x, double y)
    {
      return -(Legendre4(x) * dl1(y));
    }

    static double hdiv_leg_quad_p4_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 5 */

    static double hdiv_leg_quad_p5_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e1_a(double x, double y)
    {
      return -l0(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p5_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e1_ax(double x, double y)
    {
      return -dl0(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p5_e1_ay(double x, double y)
    {
      return -l0(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p5_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e2_a(double x, double y)
    {
      return -l1(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p5_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e2_ax(double x, double y)
    {
      return -dl1(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p5_e2_ay(double x, double y)
    {
      return -l1(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p5_e3_b(double x, double y)
    {
      return Legendre5(x) * l0(y);
    }

    static double hdiv_leg_quad_p5_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e3_bx(double x, double y)
    {
      return Legendre5x(x) * l0(y);
    }

    static double hdiv_leg_quad_p5_e3_by(double x, double y)
    {
      return Legendre5(x) * dl0(y);
    }

    static double hdiv_leg_quad_p5_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e4_b(double x, double y)
    {
      return Legendre5(x) * l1(y);
    }

    static double hdiv_leg_quad_p5_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e4_bx(double x, double y)
    {
      return Legendre5x(x) * l1(y);
    }

    static double hdiv_leg_quad_p5_e4_by(double x, double y)
    {
      return Legendre5(x) * dl1(y);
    }

    static double hdiv_leg_quad_p5_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 6 */

    static double hdiv_leg_quad_p6_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e1_a_1(double x, double y)
    {
      return -l0(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p6_e1_a_0(double x, double y)
    {
      return (l0(x) * Legendre6(y));
    }

    static double hdiv_leg_quad_p6_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e1_ax_1(double x, double y)
    {
      return -dl0(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p6_e1_ax_0(double x, double y)
    {
      return (dl0(x) * Legendre6(y));
    }

    static double hdiv_leg_quad_p6_e1_ay_1(double x, double y)
    {
      return -l0(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p6_e1_ay_0(double x, double y)
    {
      return (l0(x) * Legendre6x(y));
    }

    static double hdiv_leg_quad_p6_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e2_a_0(double x, double y)
    {
      return -l1(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p6_e2_a_1(double x, double y)
    {
      return (l1(x) * Legendre6(y));
    }

    static double hdiv_leg_quad_p6_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e2_ax_0(double x, double y)
    {
      return -dl1(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p6_e2_ax_1(double x, double y)
    {
      return (dl1(x) * Legendre6(y));
    }

    static double hdiv_leg_quad_p6_e2_ay_0(double x, double y)
    {
      return -l1(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p6_e2_ay_1(double x, double y)
    {
      return (l1(x) * Legendre6x(y));
    }

    static double hdiv_leg_quad_p6_e3_b_0(double x, double y)
    {
      return Legendre6(x) * l0(y);
    }

    static double hdiv_leg_quad_p6_e3_b_1(double x, double y)
    {
      return -(Legendre6(x) * l0(y));
    }

    static double hdiv_leg_quad_p6_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e3_bx_0(double x, double y)
    {
      return Legendre6x(x) * l0(y);
    }

    static double hdiv_leg_quad_p6_e3_bx_1(double x, double y)
    {
      return -(Legendre6x(x) * l0(y));
    }

    static double hdiv_leg_quad_p6_e3_by_0(double x, double y)
    {
      return Legendre6(x) * dl0(y);
    }

    static double hdiv_leg_quad_p6_e3_by_1(double x, double y)
    {
      return -(Legendre6(x) * dl0(y));
    }

    static double hdiv_leg_quad_p6_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e4_b_1(double x, double y)
    {
      return Legendre6(x) * l1(y);
    }

    static double hdiv_leg_quad_p6_e4_b_0(double x, double y)
    {
      return -(Legendre6(x) * l1(y));
    }

    static double hdiv_leg_quad_p6_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e4_bx_1(double x, double y)
    {
      return Legendre6x(x) * l1(y);
    }

    static double hdiv_leg_quad_p6_e4_bx_0(double x, double y)
    {
      return -(Legendre6x(x) * l1(y));
    }

    static double hdiv_leg_quad_p6_e4_by_1(double x, double y)
    {
      return Legendre6(x) * dl1(y);
    }

    static double hdiv_leg_quad_p6_e4_by_0(double x, double y)
    {
      return -(Legendre6(x) * dl1(y));
    }

    static double hdiv_leg_quad_p6_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 7 */

    static double hdiv_leg_quad_p7_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e1_a(double x, double y)
    {
      return -l0(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p7_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e1_ax(double x, double y)
    {
      return -dl0(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p7_e1_ay(double x, double y)
    {
      return -l0(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p7_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e2_a(double x, double y)
    {
      return -l1(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p7_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e2_ax(double x, double y)
    {
      return -dl1(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p7_e2_ay(double x, double y)
    {
      return -l1(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p7_e3_b(double x, double y)
    {
      return Legendre7(x) * l0(y);
    }

    static double hdiv_leg_quad_p7_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e3_bx(double x, double y)
    {
      return Legendre7x(x) * l0(y);
    }

    static double hdiv_leg_quad_p7_e3_by(double x, double y)
    {
      return Legendre7(x) * dl0(y);
    }

    static double hdiv_leg_quad_p7_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e4_b(double x, double y)
    {
      return Legendre7(x) * l1(y);
    }

    static double hdiv_leg_quad_p7_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e4_bx(double x, double y)
    {
      return Legendre7x(x) * l1(y);
    }

    static double hdiv_leg_quad_p7_e4_by(double x, double y)
    {
      return Legendre7(x) * dl1(y);
    }

    static double hdiv_leg_quad_p7_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 8 */

    static double hdiv_leg_quad_p8_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e1_a_1(double x, double y)
    {
      return -l0(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p8_e1_a_0(double x, double y)
    {
      return (l0(x) * Legendre8(y));
    }

    static double hdiv_leg_quad_p8_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e1_ax_1(double x, double y)
    {
      return -dl0(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p8_e1_ax_0(double x, double y)
    {
      return (dl0(x) * Legendre8(y));
    }

    static double hdiv_leg_quad_p8_e1_ay_1(double x, double y)
    {
      return -l0(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p8_e1_ay_0(double x, double y)
    {
      return (l0(x) * Legendre8x(y));
    }

    static double hdiv_leg_quad_p8_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e2_a_0(double x, double y)
    {
      return -l1(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p8_e2_a_1(double x, double y)
    {
      return (l1(x) * Legendre8(y));
    }

    static double hdiv_leg_quad_p8_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e2_ax_0(double x, double y)
    {
      return -dl1(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p8_e2_ax_1(double x, double y)
    {
      return (dl1(x) * Legendre8(y));
    }

    static double hdiv_leg_quad_p8_e2_ay_0(double x, double y)
    {
      return -l1(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p8_e2_ay_1(double x, double y)
    {
      return (l1(x) * Legendre8x(y));
    }

    static double hdiv_leg_quad_p8_e3_b_0(double x, double y)
    {
      return Legendre8(x) * l0(y);
    }

    static double hdiv_leg_quad_p8_e3_b_1(double x, double y)
    {
      return -(Legendre8(x) * l0(y));
    }

    static double hdiv_leg_quad_p8_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e3_bx_0(double x, double y)
    {
      return Legendre8x(x) * l0(y);
    }

    static double hdiv_leg_quad_p8_e3_bx_1(double x, double y)
    {
      return -(Legendre8x(x) * l0(y));
    }

    static double hdiv_leg_quad_p8_e3_by_0(double x, double y)
    {
      return Legendre8(x) * dl0(y);
    }

    static double hdiv_leg_quad_p8_e3_by_1(double x, double y)
    {
      return -(Legendre8(x) * dl0(y));
    }

    static double hdiv_leg_quad_p8_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e4_b_1(double x, double y)
    {
      return Legendre8(x) * l1(y);
    }

    static double hdiv_leg_quad_p8_e4_b_0(double x, double y)
    {
      return -(Legendre8(x) * l1(y));
    }

    static double hdiv_leg_quad_p8_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e4_bx_1(double x, double y)
    {
      return Legendre8x(x) * l1(y);
    }

    static double hdiv_leg_quad_p8_e4_bx_0(double x, double y)
    {
      return -(Legendre8x(x) * l1(y));
    }

    static double hdiv_leg_quad_p8_e4_by_1(double x, double y)
    {
      return Legendre8(x) * dl1(y);
    }

    static double hdiv_leg_quad_p8_e4_by_0(double x, double y)
    {
      return -(Legendre8(x) * dl1(y));
    }

    static double hdiv_leg_quad_p8_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 9 */

    static double hdiv_leg_quad_p9_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e1_a(double x, double y)
    {
      return -l0(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p9_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e1_ax(double x, double y)
    {
      return -dl0(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p9_e1_ay(double x, double y)
    {
      return -l0(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p9_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e2_a(double x, double y)
    {
      return -l1(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p9_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e2_ax(double x, double y)
    {
      return -dl1(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p9_e2_ay(double x, double y)
    {
      return -l1(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p9_e3_b(double x, double y)
    {
      return Legendre9(x) * l0(y);
    }

    static double hdiv_leg_quad_p9_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e3_bx(double x, double y)
    {
      return Legendre9x(x) * l0(y);
    }

    static double hdiv_leg_quad_p9_e3_by(double x, double y)
    {
      return Legendre9(x) * dl0(y);
    }

    static double hdiv_leg_quad_p9_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e4_b(double x, double y)
    {
      return Legendre9(x) * l1(y);
    }

    static double hdiv_leg_quad_p9_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e4_bx(double x, double y)
    {
      return Legendre9x(x) * l1(y);
    }

    static double hdiv_leg_quad_p9_e4_by(double x, double y)
    {
      return Legendre9(x) * dl1(y);
    }

    static double hdiv_leg_quad_p9_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* ORDER 10 */

    static double hdiv_leg_quad_p10_e1_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e1_a_1(double x, double y)
    {
      return -l0(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p10_e1_a_0(double x, double y)
    {
      return (l0(x) * Legendre10(y));
    }

    static double hdiv_leg_quad_p10_e1_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e1_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e1_ax_1(double x, double y)
    {
      return -dl0(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p10_e1_ax_0(double x, double y)
    {
      return (dl0(x) * Legendre10(y));
    }

    static double hdiv_leg_quad_p10_e1_ay_1(double x, double y)
    {
      return -l0(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p10_e1_ay_0(double x, double y)
    {
      return (l0(x) * Legendre10x(y));
    }

    static double hdiv_leg_quad_p10_e2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e2_a_0(double x, double y)
    {
      return -l1(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p10_e2_a_1(double x, double y)
    {
      return (l1(x) * Legendre10(y));
    }

    static double hdiv_leg_quad_p10_e2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e2_ax_0(double x, double y)
    {
      return -dl1(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p10_e2_ax_1(double x, double y)
    {
      return (dl1(x) * Legendre10(y));
    }

    static double hdiv_leg_quad_p10_e2_ay_0(double x, double y)
    {
      return -l1(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p10_e2_ay_1(double x, double y)
    {
      return (l1(x) * Legendre10x(y));
    }

    static double hdiv_leg_quad_p10_e3_b_0(double x, double y)
    {
      return Legendre10(x) * l0(y);
    }

    static double hdiv_leg_quad_p10_e3_b_1(double x, double y)
    {
      return -(Legendre10(x) * l0(y));
    }

    static double hdiv_leg_quad_p10_e3_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e3_bx_0(double x, double y)
    {
      return Legendre10x(x) * l0(y);
    }

    static double hdiv_leg_quad_p10_e3_bx_1(double x, double y)
    {
      return -(Legendre10x(x) * l0(y));
    }

    static double hdiv_leg_quad_p10_e3_by_0(double x, double y)
    {
      return Legendre10(x) * dl0(y);
    }

    static double hdiv_leg_quad_p10_e3_by_1(double x, double y)
    {
      return -(Legendre10(x) * dl0(y));
    }

    static double hdiv_leg_quad_p10_e3_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e3_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e4_b_1(double x, double y)
    {
      return Legendre10(x) * l1(y);
    }

    static double hdiv_leg_quad_p10_e4_b_0(double x, double y)
    {
      return -(Legendre10(x) * l1(y));
    }

    static double hdiv_leg_quad_p10_e4_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e4_bx_1(double x, double y)
    {
      return Legendre10x(x) * l1(y);
    }

    static double hdiv_leg_quad_p10_e4_bx_0(double x, double y)
    {
      return -(Legendre10x(x) * l1(y));
    }

    static double hdiv_leg_quad_p10_e4_by_1(double x, double y)
    {
      return Legendre10(x) * dl1(y);
    }

    static double hdiv_leg_quad_p10_e4_by_0(double x, double y)
    {
      return -(Legendre10(x) * dl1(y));
    }

    static double hdiv_leg_quad_p10_e4_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10_e4_ay(double x, double y)
    {
      return 0.0;
    }

    /* BUBBLE */

    /* BUBBLE ( 1 , 0 ) */

    static double hdiv_leg_quad_p0p2_b1_b(double x, double y)
    {
      return Legendre0(x) * l2(y);
    }

    static double hdiv_leg_quad_p0p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p2_b1_bx(double x, double y)
    {
      return Legendre0x(x) * l2(y);
    }

    static double hdiv_leg_quad_p0p2_b1_by(double x, double y)
    {
      return Legendre0(x) * dl2(y);
    }

    static double hdiv_leg_quad_p0p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p3_b1_b(double x, double y)
    {
      return Legendre0(x) * l3(y);
    }

    static double hdiv_leg_quad_p0p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p3_b1_bx(double x, double y)
    {
      return Legendre0x(x) * l3(y);
    }

    static double hdiv_leg_quad_p0p3_b1_by(double x, double y)
    {
      return Legendre0(x) * dl3(y);
    }

    static double hdiv_leg_quad_p0p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p4_b1_b(double x, double y)
    {
      return Legendre0(x) * l4(y);
    }

    static double hdiv_leg_quad_p0p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p4_b1_bx(double x, double y)
    {
      return Legendre0x(x) * l4(y);
    }

    static double hdiv_leg_quad_p0p4_b1_by(double x, double y)
    {
      return Legendre0(x) * dl4(y);
    }

    static double hdiv_leg_quad_p0p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p5_b1_b(double x, double y)
    {
      return Legendre0(x) * l5(y);
    }

    static double hdiv_leg_quad_p0p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p5_b1_bx(double x, double y)
    {
      return Legendre0x(x) * l5(y);
    }

    static double hdiv_leg_quad_p0p5_b1_by(double x, double y)
    {
      return Legendre0(x) * dl5(y);
    }

    static double hdiv_leg_quad_p0p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p6_b1_b(double x, double y)
    {
      return Legendre0(x) * l6(y);
    }

    static double hdiv_leg_quad_p0p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p6_b1_bx(double x, double y)
    {
      return Legendre0x(x) * l6(y);
    }

    static double hdiv_leg_quad_p0p6_b1_by(double x, double y)
    {
      return Legendre0(x) * dl6(y);
    }

    static double hdiv_leg_quad_p0p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p7_b1_b(double x, double y)
    {
      return Legendre0(x) * l7(y);
    }

    static double hdiv_leg_quad_p0p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p7_b1_bx(double x, double y)
    {
      return Legendre0x(x) * l7(y);
    }

    static double hdiv_leg_quad_p0p7_b1_by(double x, double y)
    {
      return Legendre0(x) * dl7(y);
    }

    static double hdiv_leg_quad_p0p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p8_b1_b(double x, double y)
    {
      return Legendre0(x) * l8(y);
    }

    static double hdiv_leg_quad_p0p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p8_b1_bx(double x, double y)
    {
      return Legendre0x(x) * l8(y);
    }

    static double hdiv_leg_quad_p0p8_b1_by(double x, double y)
    {
      return Legendre0(x) * dl8(y);
    }

    static double hdiv_leg_quad_p0p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p9_b1_b(double x, double y)
    {
      return Legendre0(x) * l9(y);
    }

    static double hdiv_leg_quad_p0p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p9_b1_bx(double x, double y)
    {
      return Legendre0x(x) * l9(y);
    }

    static double hdiv_leg_quad_p0p9_b1_by(double x, double y)
    {
      return Legendre0(x) * dl9(y);
    }

    static double hdiv_leg_quad_p0p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p10_b1_b(double x, double y)
    {
      return Legendre0(x) * l10(y);
    }

    static double hdiv_leg_quad_p0p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p10_b1_bx(double x, double y)
    {
      return Legendre0x(x) * l10(y);
    }

    static double hdiv_leg_quad_p0p10_b1_by(double x, double y)
    {
      return Legendre0(x) * dl10(y);
    }

    static double hdiv_leg_quad_p0p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p11_b1_b(double x, double y)
    {
      return Legendre0(x) * l11(y);
    }

    static double hdiv_leg_quad_p0p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p11_b1_bx(double x, double y)
    {
      return Legendre0x(x) * l11(y);
    }

    static double hdiv_leg_quad_p0p11_b1_by(double x, double y)
    {
      return Legendre0(x) * dl11(y);
    }

    static double hdiv_leg_quad_p0p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p0p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p2_b1_b(double x, double y)
    {
      return Legendre1(x) * l2(y);
    }

    static double hdiv_leg_quad_p1p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p2_b1_bx(double x, double y)
    {
      return Legendre1x(x) * l2(y);
    }

    static double hdiv_leg_quad_p1p2_b1_by(double x, double y)
    {
      return Legendre1(x) * dl2(y);
    }

    static double hdiv_leg_quad_p1p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p3_b1_b(double x, double y)
    {
      return Legendre1(x) * l3(y);
    }

    static double hdiv_leg_quad_p1p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p3_b1_bx(double x, double y)
    {
      return Legendre1x(x) * l3(y);
    }

    static double hdiv_leg_quad_p1p3_b1_by(double x, double y)
    {
      return Legendre1(x) * dl3(y);
    }

    static double hdiv_leg_quad_p1p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p4_b1_b(double x, double y)
    {
      return Legendre1(x) * l4(y);
    }

    static double hdiv_leg_quad_p1p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p4_b1_bx(double x, double y)
    {
      return Legendre1x(x) * l4(y);
    }

    static double hdiv_leg_quad_p1p4_b1_by(double x, double y)
    {
      return Legendre1(x) * dl4(y);
    }

    static double hdiv_leg_quad_p1p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p5_b1_b(double x, double y)
    {
      return Legendre1(x) * l5(y);
    }

    static double hdiv_leg_quad_p1p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p5_b1_bx(double x, double y)
    {
      return Legendre1x(x) * l5(y);
    }

    static double hdiv_leg_quad_p1p5_b1_by(double x, double y)
    {
      return Legendre1(x) * dl5(y);
    }

    static double hdiv_leg_quad_p1p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p6_b1_b(double x, double y)
    {
      return Legendre1(x) * l6(y);
    }

    static double hdiv_leg_quad_p1p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p6_b1_bx(double x, double y)
    {
      return Legendre1x(x) * l6(y);
    }

    static double hdiv_leg_quad_p1p6_b1_by(double x, double y)
    {
      return Legendre1(x) * dl6(y);
    }

    static double hdiv_leg_quad_p1p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p7_b1_b(double x, double y)
    {
      return Legendre1(x) * l7(y);
    }

    static double hdiv_leg_quad_p1p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p7_b1_bx(double x, double y)
    {
      return Legendre1x(x) * l7(y);
    }

    static double hdiv_leg_quad_p1p7_b1_by(double x, double y)
    {
      return Legendre1(x) * dl7(y);
    }

    static double hdiv_leg_quad_p1p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p8_b1_b(double x, double y)
    {
      return Legendre1(x) * l8(y);
    }

    static double hdiv_leg_quad_p1p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p8_b1_bx(double x, double y)
    {
      return Legendre1x(x) * l8(y);
    }

    static double hdiv_leg_quad_p1p8_b1_by(double x, double y)
    {
      return Legendre1(x) * dl8(y);
    }

    static double hdiv_leg_quad_p1p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p9_b1_b(double x, double y)
    {
      return Legendre1(x) * l9(y);
    }

    static double hdiv_leg_quad_p1p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p9_b1_bx(double x, double y)
    {
      return Legendre1x(x) * l9(y);
    }

    static double hdiv_leg_quad_p1p9_b1_by(double x, double y)
    {
      return Legendre1(x) * dl9(y);
    }

    static double hdiv_leg_quad_p1p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p10_b1_b(double x, double y)
    {
      return Legendre1(x) * l10(y);
    }

    static double hdiv_leg_quad_p1p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p10_b1_bx(double x, double y)
    {
      return Legendre1x(x) * l10(y);
    }

    static double hdiv_leg_quad_p1p10_b1_by(double x, double y)
    {
      return Legendre1(x) * dl10(y);
    }

    static double hdiv_leg_quad_p1p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p11_b1_b(double x, double y)
    {
      return Legendre1(x) * l11(y);
    }

    static double hdiv_leg_quad_p1p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p11_b1_bx(double x, double y)
    {
      return Legendre1x(x) * l11(y);
    }

    static double hdiv_leg_quad_p1p11_b1_by(double x, double y)
    {
      return Legendre1(x) * dl11(y);
    }

    static double hdiv_leg_quad_p1p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p1p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p2_b1_b(double x, double y)
    {
      return Legendre2(x) * l2(y);
    }

    static double hdiv_leg_quad_p2p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p2_b1_bx(double x, double y)
    {
      return Legendre2x(x) * l2(y);
    }

    static double hdiv_leg_quad_p2p2_b1_by(double x, double y)
    {
      return Legendre2(x) * dl2(y);
    }

    static double hdiv_leg_quad_p2p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p3_b1_b(double x, double y)
    {
      return Legendre2(x) * l3(y);
    }

    static double hdiv_leg_quad_p2p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p3_b1_bx(double x, double y)
    {
      return Legendre2x(x) * l3(y);
    }

    static double hdiv_leg_quad_p2p3_b1_by(double x, double y)
    {
      return Legendre2(x) * dl3(y);
    }

    static double hdiv_leg_quad_p2p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p4_b1_b(double x, double y)
    {
      return Legendre2(x) * l4(y);
    }

    static double hdiv_leg_quad_p2p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p4_b1_bx(double x, double y)
    {
      return Legendre2x(x) * l4(y);
    }

    static double hdiv_leg_quad_p2p4_b1_by(double x, double y)
    {
      return Legendre2(x) * dl4(y);
    }

    static double hdiv_leg_quad_p2p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p5_b1_b(double x, double y)
    {
      return Legendre2(x) * l5(y);
    }

    static double hdiv_leg_quad_p2p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p5_b1_bx(double x, double y)
    {
      return Legendre2x(x) * l5(y);
    }

    static double hdiv_leg_quad_p2p5_b1_by(double x, double y)
    {
      return Legendre2(x) * dl5(y);
    }

    static double hdiv_leg_quad_p2p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p6_b1_b(double x, double y)
    {
      return Legendre2(x) * l6(y);
    }

    static double hdiv_leg_quad_p2p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p6_b1_bx(double x, double y)
    {
      return Legendre2x(x) * l6(y);
    }

    static double hdiv_leg_quad_p2p6_b1_by(double x, double y)
    {
      return Legendre2(x) * dl6(y);
    }

    static double hdiv_leg_quad_p2p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p7_b1_b(double x, double y)
    {
      return Legendre2(x) * l7(y);
    }

    static double hdiv_leg_quad_p2p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p7_b1_bx(double x, double y)
    {
      return Legendre2x(x) * l7(y);
    }

    static double hdiv_leg_quad_p2p7_b1_by(double x, double y)
    {
      return Legendre2(x) * dl7(y);
    }

    static double hdiv_leg_quad_p2p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p8_b1_b(double x, double y)
    {
      return Legendre2(x) * l8(y);
    }

    static double hdiv_leg_quad_p2p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p8_b1_bx(double x, double y)
    {
      return Legendre2x(x) * l8(y);
    }

    static double hdiv_leg_quad_p2p8_b1_by(double x, double y)
    {
      return Legendre2(x) * dl8(y);
    }

    static double hdiv_leg_quad_p2p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p9_b1_b(double x, double y)
    {
      return Legendre2(x) * l9(y);
    }

    static double hdiv_leg_quad_p2p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p9_b1_bx(double x, double y)
    {
      return Legendre2x(x) * l9(y);
    }

    static double hdiv_leg_quad_p2p9_b1_by(double x, double y)
    {
      return Legendre2(x) * dl9(y);
    }

    static double hdiv_leg_quad_p2p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p10_b1_b(double x, double y)
    {
      return Legendre2(x) * l10(y);
    }

    static double hdiv_leg_quad_p2p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p10_b1_bx(double x, double y)
    {
      return Legendre2x(x) * l10(y);
    }

    static double hdiv_leg_quad_p2p10_b1_by(double x, double y)
    {
      return Legendre2(x) * dl10(y);
    }

    static double hdiv_leg_quad_p2p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p11_b1_b(double x, double y)
    {
      return Legendre2(x) * l11(y);
    }

    static double hdiv_leg_quad_p2p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p11_b1_bx(double x, double y)
    {
      return Legendre2x(x) * l11(y);
    }

    static double hdiv_leg_quad_p2p11_b1_by(double x, double y)
    {
      return Legendre2(x) * dl11(y);
    }

    static double hdiv_leg_quad_p2p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p2_b1_b(double x, double y)
    {
      return Legendre3(x) * l2(y);
    }

    static double hdiv_leg_quad_p3p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p2_b1_bx(double x, double y)
    {
      return Legendre3x(x) * l2(y);
    }

    static double hdiv_leg_quad_p3p2_b1_by(double x, double y)
    {
      return Legendre3(x) * dl2(y);
    }

    static double hdiv_leg_quad_p3p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p3_b1_b(double x, double y)
    {
      return Legendre3(x) * l3(y);
    }

    static double hdiv_leg_quad_p3p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p3_b1_bx(double x, double y)
    {
      return Legendre3x(x) * l3(y);
    }

    static double hdiv_leg_quad_p3p3_b1_by(double x, double y)
    {
      return Legendre3(x) * dl3(y);
    }

    static double hdiv_leg_quad_p3p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p4_b1_b(double x, double y)
    {
      return Legendre3(x) * l4(y);
    }

    static double hdiv_leg_quad_p3p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p4_b1_bx(double x, double y)
    {
      return Legendre3x(x) * l4(y);
    }

    static double hdiv_leg_quad_p3p4_b1_by(double x, double y)
    {
      return Legendre3(x) * dl4(y);
    }

    static double hdiv_leg_quad_p3p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p5_b1_b(double x, double y)
    {
      return Legendre3(x) * l5(y);
    }

    static double hdiv_leg_quad_p3p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p5_b1_bx(double x, double y)
    {
      return Legendre3x(x) * l5(y);
    }

    static double hdiv_leg_quad_p3p5_b1_by(double x, double y)
    {
      return Legendre3(x) * dl5(y);
    }

    static double hdiv_leg_quad_p3p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p6_b1_b(double x, double y)
    {
      return Legendre3(x) * l6(y);
    }

    static double hdiv_leg_quad_p3p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p6_b1_bx(double x, double y)
    {
      return Legendre3x(x) * l6(y);
    }

    static double hdiv_leg_quad_p3p6_b1_by(double x, double y)
    {
      return Legendre3(x) * dl6(y);
    }

    static double hdiv_leg_quad_p3p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p7_b1_b(double x, double y)
    {
      return Legendre3(x) * l7(y);
    }

    static double hdiv_leg_quad_p3p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p7_b1_bx(double x, double y)
    {
      return Legendre3x(x) * l7(y);
    }

    static double hdiv_leg_quad_p3p7_b1_by(double x, double y)
    {
      return Legendre3(x) * dl7(y);
    }

    static double hdiv_leg_quad_p3p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p8_b1_b(double x, double y)
    {
      return Legendre3(x) * l8(y);
    }

    static double hdiv_leg_quad_p3p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p8_b1_bx(double x, double y)
    {
      return Legendre3x(x) * l8(y);
    }

    static double hdiv_leg_quad_p3p8_b1_by(double x, double y)
    {
      return Legendre3(x) * dl8(y);
    }

    static double hdiv_leg_quad_p3p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p9_b1_b(double x, double y)
    {
      return Legendre3(x) * l9(y);
    }

    static double hdiv_leg_quad_p3p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p9_b1_bx(double x, double y)
    {
      return Legendre3x(x) * l9(y);
    }

    static double hdiv_leg_quad_p3p9_b1_by(double x, double y)
    {
      return Legendre3(x) * dl9(y);
    }

    static double hdiv_leg_quad_p3p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p10_b1_b(double x, double y)
    {
      return Legendre3(x) * l10(y);
    }

    static double hdiv_leg_quad_p3p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p10_b1_bx(double x, double y)
    {
      return Legendre3x(x) * l10(y);
    }

    static double hdiv_leg_quad_p3p10_b1_by(double x, double y)
    {
      return Legendre3(x) * dl10(y);
    }

    static double hdiv_leg_quad_p3p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p11_b1_b(double x, double y)
    {
      return Legendre3(x) * l11(y);
    }

    static double hdiv_leg_quad_p3p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p11_b1_bx(double x, double y)
    {
      return Legendre3x(x) * l11(y);
    }

    static double hdiv_leg_quad_p3p11_b1_by(double x, double y)
    {
      return Legendre3(x) * dl11(y);
    }

    static double hdiv_leg_quad_p3p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p2_b1_b(double x, double y)
    {
      return Legendre4(x) * l2(y);
    }

    static double hdiv_leg_quad_p4p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p2_b1_bx(double x, double y)
    {
      return Legendre4x(x) * l2(y);
    }

    static double hdiv_leg_quad_p4p2_b1_by(double x, double y)
    {
      return Legendre4(x) * dl2(y);
    }

    static double hdiv_leg_quad_p4p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p3_b1_b(double x, double y)
    {
      return Legendre4(x) * l3(y);
    }

    static double hdiv_leg_quad_p4p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p3_b1_bx(double x, double y)
    {
      return Legendre4x(x) * l3(y);
    }

    static double hdiv_leg_quad_p4p3_b1_by(double x, double y)
    {
      return Legendre4(x) * dl3(y);
    }

    static double hdiv_leg_quad_p4p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p4_b1_b(double x, double y)
    {
      return Legendre4(x) * l4(y);
    }

    static double hdiv_leg_quad_p4p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p4_b1_bx(double x, double y)
    {
      return Legendre4x(x) * l4(y);
    }

    static double hdiv_leg_quad_p4p4_b1_by(double x, double y)
    {
      return Legendre4(x) * dl4(y);
    }

    static double hdiv_leg_quad_p4p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p5_b1_b(double x, double y)
    {
      return Legendre4(x) * l5(y);
    }

    static double hdiv_leg_quad_p4p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p5_b1_bx(double x, double y)
    {
      return Legendre4x(x) * l5(y);
    }

    static double hdiv_leg_quad_p4p5_b1_by(double x, double y)
    {
      return Legendre4(x) * dl5(y);
    }

    static double hdiv_leg_quad_p4p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p6_b1_b(double x, double y)
    {
      return Legendre4(x) * l6(y);
    }

    static double hdiv_leg_quad_p4p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p6_b1_bx(double x, double y)
    {
      return Legendre4x(x) * l6(y);
    }

    static double hdiv_leg_quad_p4p6_b1_by(double x, double y)
    {
      return Legendre4(x) * dl6(y);
    }

    static double hdiv_leg_quad_p4p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p7_b1_b(double x, double y)
    {
      return Legendre4(x) * l7(y);
    }

    static double hdiv_leg_quad_p4p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p7_b1_bx(double x, double y)
    {
      return Legendre4x(x) * l7(y);
    }

    static double hdiv_leg_quad_p4p7_b1_by(double x, double y)
    {
      return Legendre4(x) * dl7(y);
    }

    static double hdiv_leg_quad_p4p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p8_b1_b(double x, double y)
    {
      return Legendre4(x) * l8(y);
    }

    static double hdiv_leg_quad_p4p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p8_b1_bx(double x, double y)
    {
      return Legendre4x(x) * l8(y);
    }

    static double hdiv_leg_quad_p4p8_b1_by(double x, double y)
    {
      return Legendre4(x) * dl8(y);
    }

    static double hdiv_leg_quad_p4p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p9_b1_b(double x, double y)
    {
      return Legendre4(x) * l9(y);
    }

    static double hdiv_leg_quad_p4p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p9_b1_bx(double x, double y)
    {
      return Legendre4x(x) * l9(y);
    }

    static double hdiv_leg_quad_p4p9_b1_by(double x, double y)
    {
      return Legendre4(x) * dl9(y);
    }

    static double hdiv_leg_quad_p4p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p10_b1_b(double x, double y)
    {
      return Legendre4(x) * l10(y);
    }

    static double hdiv_leg_quad_p4p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p10_b1_bx(double x, double y)
    {
      return Legendre4x(x) * l10(y);
    }

    static double hdiv_leg_quad_p4p10_b1_by(double x, double y)
    {
      return Legendre4(x) * dl10(y);
    }

    static double hdiv_leg_quad_p4p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p11_b1_b(double x, double y)
    {
      return Legendre4(x) * l11(y);
    }

    static double hdiv_leg_quad_p4p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p11_b1_bx(double x, double y)
    {
      return Legendre4x(x) * l11(y);
    }

    static double hdiv_leg_quad_p4p11_b1_by(double x, double y)
    {
      return Legendre4(x) * dl11(y);
    }

    static double hdiv_leg_quad_p4p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p2_b1_b(double x, double y)
    {
      return Legendre5(x) * l2(y);
    }

    static double hdiv_leg_quad_p5p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p2_b1_bx(double x, double y)
    {
      return Legendre5x(x) * l2(y);
    }

    static double hdiv_leg_quad_p5p2_b1_by(double x, double y)
    {
      return Legendre5(x) * dl2(y);
    }

    static double hdiv_leg_quad_p5p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p3_b1_b(double x, double y)
    {
      return Legendre5(x) * l3(y);
    }

    static double hdiv_leg_quad_p5p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p3_b1_bx(double x, double y)
    {
      return Legendre5x(x) * l3(y);
    }

    static double hdiv_leg_quad_p5p3_b1_by(double x, double y)
    {
      return Legendre5(x) * dl3(y);
    }

    static double hdiv_leg_quad_p5p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p4_b1_b(double x, double y)
    {
      return Legendre5(x) * l4(y);
    }

    static double hdiv_leg_quad_p5p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p4_b1_bx(double x, double y)
    {
      return Legendre5x(x) * l4(y);
    }

    static double hdiv_leg_quad_p5p4_b1_by(double x, double y)
    {
      return Legendre5(x) * dl4(y);
    }

    static double hdiv_leg_quad_p5p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p5_b1_b(double x, double y)
    {
      return Legendre5(x) * l5(y);
    }

    static double hdiv_leg_quad_p5p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p5_b1_bx(double x, double y)
    {
      return Legendre5x(x) * l5(y);
    }

    static double hdiv_leg_quad_p5p5_b1_by(double x, double y)
    {
      return Legendre5(x) * dl5(y);
    }

    static double hdiv_leg_quad_p5p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p6_b1_b(double x, double y)
    {
      return Legendre5(x) * l6(y);
    }

    static double hdiv_leg_quad_p5p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p6_b1_bx(double x, double y)
    {
      return Legendre5x(x) * l6(y);
    }

    static double hdiv_leg_quad_p5p6_b1_by(double x, double y)
    {
      return Legendre5(x) * dl6(y);
    }

    static double hdiv_leg_quad_p5p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p7_b1_b(double x, double y)
    {
      return Legendre5(x) * l7(y);
    }

    static double hdiv_leg_quad_p5p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p7_b1_bx(double x, double y)
    {
      return Legendre5x(x) * l7(y);
    }

    static double hdiv_leg_quad_p5p7_b1_by(double x, double y)
    {
      return Legendre5(x) * dl7(y);
    }

    static double hdiv_leg_quad_p5p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p8_b1_b(double x, double y)
    {
      return Legendre5(x) * l8(y);
    }

    static double hdiv_leg_quad_p5p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p8_b1_bx(double x, double y)
    {
      return Legendre5x(x) * l8(y);
    }

    static double hdiv_leg_quad_p5p8_b1_by(double x, double y)
    {
      return Legendre5(x) * dl8(y);
    }

    static double hdiv_leg_quad_p5p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p9_b1_b(double x, double y)
    {
      return Legendre5(x) * l9(y);
    }

    static double hdiv_leg_quad_p5p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p9_b1_bx(double x, double y)
    {
      return Legendre5x(x) * l9(y);
    }

    static double hdiv_leg_quad_p5p9_b1_by(double x, double y)
    {
      return Legendre5(x) * dl9(y);
    }

    static double hdiv_leg_quad_p5p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p10_b1_b(double x, double y)
    {
      return Legendre5(x) * l10(y);
    }

    static double hdiv_leg_quad_p5p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p10_b1_bx(double x, double y)
    {
      return Legendre5x(x) * l10(y);
    }

    static double hdiv_leg_quad_p5p10_b1_by(double x, double y)
    {
      return Legendre5(x) * dl10(y);
    }

    static double hdiv_leg_quad_p5p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p11_b1_b(double x, double y)
    {
      return Legendre5(x) * l11(y);
    }

    static double hdiv_leg_quad_p5p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p11_b1_bx(double x, double y)
    {
      return Legendre5x(x) * l11(y);
    }

    static double hdiv_leg_quad_p5p11_b1_by(double x, double y)
    {
      return Legendre5(x) * dl11(y);
    }

    static double hdiv_leg_quad_p5p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p2_b1_b(double x, double y)
    {
      return Legendre6(x) * l2(y);
    }

    static double hdiv_leg_quad_p6p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p2_b1_bx(double x, double y)
    {
      return Legendre6x(x) * l2(y);
    }

    static double hdiv_leg_quad_p6p2_b1_by(double x, double y)
    {
      return Legendre6(x) * dl2(y);
    }

    static double hdiv_leg_quad_p6p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p3_b1_b(double x, double y)
    {
      return Legendre6(x) * l3(y);
    }

    static double hdiv_leg_quad_p6p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p3_b1_bx(double x, double y)
    {
      return Legendre6x(x) * l3(y);
    }

    static double hdiv_leg_quad_p6p3_b1_by(double x, double y)
    {
      return Legendre6(x) * dl3(y);
    }

    static double hdiv_leg_quad_p6p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p4_b1_b(double x, double y)
    {
      return Legendre6(x) * l4(y);
    }

    static double hdiv_leg_quad_p6p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p4_b1_bx(double x, double y)
    {
      return Legendre6x(x) * l4(y);
    }

    static double hdiv_leg_quad_p6p4_b1_by(double x, double y)
    {
      return Legendre6(x) * dl4(y);
    }

    static double hdiv_leg_quad_p6p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p5_b1_b(double x, double y)
    {
      return Legendre6(x) * l5(y);
    }

    static double hdiv_leg_quad_p6p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p5_b1_bx(double x, double y)
    {
      return Legendre6x(x) * l5(y);
    }

    static double hdiv_leg_quad_p6p5_b1_by(double x, double y)
    {
      return Legendre6(x) * dl5(y);
    }

    static double hdiv_leg_quad_p6p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p6_b1_b(double x, double y)
    {
      return Legendre6(x) * l6(y);
    }

    static double hdiv_leg_quad_p6p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p6_b1_bx(double x, double y)
    {
      return Legendre6x(x) * l6(y);
    }

    static double hdiv_leg_quad_p6p6_b1_by(double x, double y)
    {
      return Legendre6(x) * dl6(y);
    }

    static double hdiv_leg_quad_p6p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p7_b1_b(double x, double y)
    {
      return Legendre6(x) * l7(y);
    }

    static double hdiv_leg_quad_p6p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p7_b1_bx(double x, double y)
    {
      return Legendre6x(x) * l7(y);
    }

    static double hdiv_leg_quad_p6p7_b1_by(double x, double y)
    {
      return Legendre6(x) * dl7(y);
    }

    static double hdiv_leg_quad_p6p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p8_b1_b(double x, double y)
    {
      return Legendre6(x) * l8(y);
    }

    static double hdiv_leg_quad_p6p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p8_b1_bx(double x, double y)
    {
      return Legendre6x(x) * l8(y);
    }

    static double hdiv_leg_quad_p6p8_b1_by(double x, double y)
    {
      return Legendre6(x) * dl8(y);
    }

    static double hdiv_leg_quad_p6p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p9_b1_b(double x, double y)
    {
      return Legendre6(x) * l9(y);
    }

    static double hdiv_leg_quad_p6p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p9_b1_bx(double x, double y)
    {
      return Legendre6x(x) * l9(y);
    }

    static double hdiv_leg_quad_p6p9_b1_by(double x, double y)
    {
      return Legendre6(x) * dl9(y);
    }

    static double hdiv_leg_quad_p6p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p10_b1_b(double x, double y)
    {
      return Legendre6(x) * l10(y);
    }

    static double hdiv_leg_quad_p6p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p10_b1_bx(double x, double y)
    {
      return Legendre6x(x) * l10(y);
    }

    static double hdiv_leg_quad_p6p10_b1_by(double x, double y)
    {
      return Legendre6(x) * dl10(y);
    }

    static double hdiv_leg_quad_p6p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p11_b1_b(double x, double y)
    {
      return Legendre6(x) * l11(y);
    }

    static double hdiv_leg_quad_p6p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p11_b1_bx(double x, double y)
    {
      return Legendre6x(x) * l11(y);
    }

    static double hdiv_leg_quad_p6p11_b1_by(double x, double y)
    {
      return Legendre6(x) * dl11(y);
    }

    static double hdiv_leg_quad_p6p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p2_b1_b(double x, double y)
    {
      return Legendre7(x) * l2(y);
    }

    static double hdiv_leg_quad_p7p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p2_b1_bx(double x, double y)
    {
      return Legendre7x(x) * l2(y);
    }

    static double hdiv_leg_quad_p7p2_b1_by(double x, double y)
    {
      return Legendre7(x) * dl2(y);
    }

    static double hdiv_leg_quad_p7p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p3_b1_b(double x, double y)
    {
      return Legendre7(x) * l3(y);
    }

    static double hdiv_leg_quad_p7p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p3_b1_bx(double x, double y)
    {
      return Legendre7x(x) * l3(y);
    }

    static double hdiv_leg_quad_p7p3_b1_by(double x, double y)
    {
      return Legendre7(x) * dl3(y);
    }

    static double hdiv_leg_quad_p7p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p4_b1_b(double x, double y)
    {
      return Legendre7(x) * l4(y);
    }

    static double hdiv_leg_quad_p7p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p4_b1_bx(double x, double y)
    {
      return Legendre7x(x) * l4(y);
    }

    static double hdiv_leg_quad_p7p4_b1_by(double x, double y)
    {
      return Legendre7(x) * dl4(y);
    }

    static double hdiv_leg_quad_p7p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p5_b1_b(double x, double y)
    {
      return Legendre7(x) * l5(y);
    }

    static double hdiv_leg_quad_p7p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p5_b1_bx(double x, double y)
    {
      return Legendre7x(x) * l5(y);
    }

    static double hdiv_leg_quad_p7p5_b1_by(double x, double y)
    {
      return Legendre7(x) * dl5(y);
    }

    static double hdiv_leg_quad_p7p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p6_b1_b(double x, double y)
    {
      return Legendre7(x) * l6(y);
    }

    static double hdiv_leg_quad_p7p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p6_b1_bx(double x, double y)
    {
      return Legendre7x(x) * l6(y);
    }

    static double hdiv_leg_quad_p7p6_b1_by(double x, double y)
    {
      return Legendre7(x) * dl6(y);
    }

    static double hdiv_leg_quad_p7p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p7_b1_b(double x, double y)
    {
      return Legendre7(x) * l7(y);
    }

    static double hdiv_leg_quad_p7p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p7_b1_bx(double x, double y)
    {
      return Legendre7x(x) * l7(y);
    }

    static double hdiv_leg_quad_p7p7_b1_by(double x, double y)
    {
      return Legendre7(x) * dl7(y);
    }

    static double hdiv_leg_quad_p7p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p8_b1_b(double x, double y)
    {
      return Legendre7(x) * l8(y);
    }

    static double hdiv_leg_quad_p7p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p8_b1_bx(double x, double y)
    {
      return Legendre7x(x) * l8(y);
    }

    static double hdiv_leg_quad_p7p8_b1_by(double x, double y)
    {
      return Legendre7(x) * dl8(y);
    }

    static double hdiv_leg_quad_p7p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p9_b1_b(double x, double y)
    {
      return Legendre7(x) * l9(y);
    }

    static double hdiv_leg_quad_p7p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p9_b1_bx(double x, double y)
    {
      return Legendre7x(x) * l9(y);
    }

    static double hdiv_leg_quad_p7p9_b1_by(double x, double y)
    {
      return Legendre7(x) * dl9(y);
    }

    static double hdiv_leg_quad_p7p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p10_b1_b(double x, double y)
    {
      return Legendre7(x) * l10(y);
    }

    static double hdiv_leg_quad_p7p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p10_b1_bx(double x, double y)
    {
      return Legendre7x(x) * l10(y);
    }

    static double hdiv_leg_quad_p7p10_b1_by(double x, double y)
    {
      return Legendre7(x) * dl10(y);
    }

    static double hdiv_leg_quad_p7p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p11_b1_b(double x, double y)
    {
      return Legendre7(x) * l11(y);
    }

    static double hdiv_leg_quad_p7p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p11_b1_bx(double x, double y)
    {
      return Legendre7x(x) * l11(y);
    }

    static double hdiv_leg_quad_p7p11_b1_by(double x, double y)
    {
      return Legendre7(x) * dl11(y);
    }

    static double hdiv_leg_quad_p7p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p2_b1_b(double x, double y)
    {
      return Legendre8(x) * l2(y);
    }

    static double hdiv_leg_quad_p8p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p2_b1_bx(double x, double y)
    {
      return Legendre8x(x) * l2(y);
    }

    static double hdiv_leg_quad_p8p2_b1_by(double x, double y)
    {
      return Legendre8(x) * dl2(y);
    }

    static double hdiv_leg_quad_p8p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p3_b1_b(double x, double y)
    {
      return Legendre8(x) * l3(y);
    }

    static double hdiv_leg_quad_p8p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p3_b1_bx(double x, double y)
    {
      return Legendre8x(x) * l3(y);
    }

    static double hdiv_leg_quad_p8p3_b1_by(double x, double y)
    {
      return Legendre8(x) * dl3(y);
    }

    static double hdiv_leg_quad_p8p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p4_b1_b(double x, double y)
    {
      return Legendre8(x) * l4(y);
    }

    static double hdiv_leg_quad_p8p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p4_b1_bx(double x, double y)
    {
      return Legendre8x(x) * l4(y);
    }

    static double hdiv_leg_quad_p8p4_b1_by(double x, double y)
    {
      return Legendre8(x) * dl4(y);
    }

    static double hdiv_leg_quad_p8p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p5_b1_b(double x, double y)
    {
      return Legendre8(x) * l5(y);
    }

    static double hdiv_leg_quad_p8p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p5_b1_bx(double x, double y)
    {
      return Legendre8x(x) * l5(y);
    }

    static double hdiv_leg_quad_p8p5_b1_by(double x, double y)
    {
      return Legendre8(x) * dl5(y);
    }

    static double hdiv_leg_quad_p8p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p6_b1_b(double x, double y)
    {
      return Legendre8(x) * l6(y);
    }

    static double hdiv_leg_quad_p8p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p6_b1_bx(double x, double y)
    {
      return Legendre8x(x) * l6(y);
    }

    static double hdiv_leg_quad_p8p6_b1_by(double x, double y)
    {
      return Legendre8(x) * dl6(y);
    }

    static double hdiv_leg_quad_p8p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p7_b1_b(double x, double y)
    {
      return Legendre8(x) * l7(y);
    }

    static double hdiv_leg_quad_p8p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p7_b1_bx(double x, double y)
    {
      return Legendre8x(x) * l7(y);
    }

    static double hdiv_leg_quad_p8p7_b1_by(double x, double y)
    {
      return Legendre8(x) * dl7(y);
    }

    static double hdiv_leg_quad_p8p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p8_b1_b(double x, double y)
    {
      return Legendre8(x) * l8(y);
    }

    static double hdiv_leg_quad_p8p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p8_b1_bx(double x, double y)
    {
      return Legendre8x(x) * l8(y);
    }

    static double hdiv_leg_quad_p8p8_b1_by(double x, double y)
    {
      return Legendre8(x) * dl8(y);
    }

    static double hdiv_leg_quad_p8p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p9_b1_b(double x, double y)
    {
      return Legendre8(x) * l9(y);
    }

    static double hdiv_leg_quad_p8p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p9_b1_bx(double x, double y)
    {
      return Legendre8x(x) * l9(y);
    }

    static double hdiv_leg_quad_p8p9_b1_by(double x, double y)
    {
      return Legendre8(x) * dl9(y);
    }

    static double hdiv_leg_quad_p8p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p10_b1_b(double x, double y)
    {
      return Legendre8(x) * l10(y);
    }

    static double hdiv_leg_quad_p8p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p10_b1_bx(double x, double y)
    {
      return Legendre8x(x) * l10(y);
    }

    static double hdiv_leg_quad_p8p10_b1_by(double x, double y)
    {
      return Legendre8(x) * dl10(y);
    }

    static double hdiv_leg_quad_p8p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p11_b1_b(double x, double y)
    {
      return Legendre8(x) * l11(y);
    }

    static double hdiv_leg_quad_p8p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p11_b1_bx(double x, double y)
    {
      return Legendre8x(x) * l11(y);
    }

    static double hdiv_leg_quad_p8p11_b1_by(double x, double y)
    {
      return Legendre8(x) * dl11(y);
    }

    static double hdiv_leg_quad_p8p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p2_b1_b(double x, double y)
    {
      return Legendre9(x) * l2(y);
    }

    static double hdiv_leg_quad_p9p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p2_b1_bx(double x, double y)
    {
      return Legendre9x(x) * l2(y);
    }

    static double hdiv_leg_quad_p9p2_b1_by(double x, double y)
    {
      return Legendre9(x) * dl2(y);
    }

    static double hdiv_leg_quad_p9p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p3_b1_b(double x, double y)
    {
      return Legendre9(x) * l3(y);
    }

    static double hdiv_leg_quad_p9p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p3_b1_bx(double x, double y)
    {
      return Legendre9x(x) * l3(y);
    }

    static double hdiv_leg_quad_p9p3_b1_by(double x, double y)
    {
      return Legendre9(x) * dl3(y);
    }

    static double hdiv_leg_quad_p9p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p4_b1_b(double x, double y)
    {
      return Legendre9(x) * l4(y);
    }

    static double hdiv_leg_quad_p9p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p4_b1_bx(double x, double y)
    {
      return Legendre9x(x) * l4(y);
    }

    static double hdiv_leg_quad_p9p4_b1_by(double x, double y)
    {
      return Legendre9(x) * dl4(y);
    }

    static double hdiv_leg_quad_p9p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p5_b1_b(double x, double y)
    {
      return Legendre9(x) * l5(y);
    }

    static double hdiv_leg_quad_p9p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p5_b1_bx(double x, double y)
    {
      return Legendre9x(x) * l5(y);
    }

    static double hdiv_leg_quad_p9p5_b1_by(double x, double y)
    {
      return Legendre9(x) * dl5(y);
    }

    static double hdiv_leg_quad_p9p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p6_b1_b(double x, double y)
    {
      return Legendre9(x) * l6(y);
    }

    static double hdiv_leg_quad_p9p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p6_b1_bx(double x, double y)
    {
      return Legendre9x(x) * l6(y);
    }

    static double hdiv_leg_quad_p9p6_b1_by(double x, double y)
    {
      return Legendre9(x) * dl6(y);
    }

    static double hdiv_leg_quad_p9p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p7_b1_b(double x, double y)
    {
      return Legendre9(x) * l7(y);
    }

    static double hdiv_leg_quad_p9p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p7_b1_bx(double x, double y)
    {
      return Legendre9x(x) * l7(y);
    }

    static double hdiv_leg_quad_p9p7_b1_by(double x, double y)
    {
      return Legendre9(x) * dl7(y);
    }

    static double hdiv_leg_quad_p9p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p8_b1_b(double x, double y)
    {
      return Legendre9(x) * l8(y);
    }

    static double hdiv_leg_quad_p9p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p8_b1_bx(double x, double y)
    {
      return Legendre9x(x) * l8(y);
    }

    static double hdiv_leg_quad_p9p8_b1_by(double x, double y)
    {
      return Legendre9(x) * dl8(y);
    }

    static double hdiv_leg_quad_p9p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p9_b1_b(double x, double y)
    {
      return Legendre9(x) * l9(y);
    }

    static double hdiv_leg_quad_p9p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p9_b1_bx(double x, double y)
    {
      return Legendre9x(x) * l9(y);
    }

    static double hdiv_leg_quad_p9p9_b1_by(double x, double y)
    {
      return Legendre9(x) * dl9(y);
    }

    static double hdiv_leg_quad_p9p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p10_b1_b(double x, double y)
    {
      return Legendre9(x) * l10(y);
    }

    static double hdiv_leg_quad_p9p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p10_b1_bx(double x, double y)
    {
      return Legendre9x(x) * l10(y);
    }

    static double hdiv_leg_quad_p9p10_b1_by(double x, double y)
    {
      return Legendre9(x) * dl10(y);
    }

    static double hdiv_leg_quad_p9p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p11_b1_b(double x, double y)
    {
      return Legendre9(x) * l11(y);
    }

    static double hdiv_leg_quad_p9p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p11_b1_bx(double x, double y)
    {
      return Legendre9x(x) * l11(y);
    }

    static double hdiv_leg_quad_p9p11_b1_by(double x, double y)
    {
      return Legendre9(x) * dl11(y);
    }

    static double hdiv_leg_quad_p9p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p2_b1_b(double x, double y)
    {
      return Legendre10(x) * l2(y);
    }

    static double hdiv_leg_quad_p10p2_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p2_b1_bx(double x, double y)
    {
      return Legendre10x(x) * l2(y);
    }

    static double hdiv_leg_quad_p10p2_b1_by(double x, double y)
    {
      return Legendre10(x) * dl2(y);
    }

    static double hdiv_leg_quad_p10p2_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p2_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p3_b1_b(double x, double y)
    {
      return Legendre10(x) * l3(y);
    }

    static double hdiv_leg_quad_p10p3_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p3_b1_bx(double x, double y)
    {
      return Legendre10x(x) * l3(y);
    }

    static double hdiv_leg_quad_p10p3_b1_by(double x, double y)
    {
      return Legendre10(x) * dl3(y);
    }

    static double hdiv_leg_quad_p10p3_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p3_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p4_b1_b(double x, double y)
    {
      return Legendre10(x) * l4(y);
    }

    static double hdiv_leg_quad_p10p4_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p4_b1_bx(double x, double y)
    {
      return Legendre10x(x) * l4(y);
    }

    static double hdiv_leg_quad_p10p4_b1_by(double x, double y)
    {
      return Legendre10(x) * dl4(y);
    }

    static double hdiv_leg_quad_p10p4_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p4_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p5_b1_b(double x, double y)
    {
      return Legendre10(x) * l5(y);
    }

    static double hdiv_leg_quad_p10p5_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p5_b1_bx(double x, double y)
    {
      return Legendre10x(x) * l5(y);
    }

    static double hdiv_leg_quad_p10p5_b1_by(double x, double y)
    {
      return Legendre10(x) * dl5(y);
    }

    static double hdiv_leg_quad_p10p5_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p5_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p6_b1_b(double x, double y)
    {
      return Legendre10(x) * l6(y);
    }

    static double hdiv_leg_quad_p10p6_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p6_b1_bx(double x, double y)
    {
      return Legendre10x(x) * l6(y);
    }

    static double hdiv_leg_quad_p10p6_b1_by(double x, double y)
    {
      return Legendre10(x) * dl6(y);
    }

    static double hdiv_leg_quad_p10p6_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p6_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p7_b1_b(double x, double y)
    {
      return Legendre10(x) * l7(y);
    }

    static double hdiv_leg_quad_p10p7_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p7_b1_bx(double x, double y)
    {
      return Legendre10x(x) * l7(y);
    }

    static double hdiv_leg_quad_p10p7_b1_by(double x, double y)
    {
      return Legendre10(x) * dl7(y);
    }

    static double hdiv_leg_quad_p10p7_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p7_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p8_b1_b(double x, double y)
    {
      return Legendre10(x) * l8(y);
    }

    static double hdiv_leg_quad_p10p8_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p8_b1_bx(double x, double y)
    {
      return Legendre10x(x) * l8(y);
    }

    static double hdiv_leg_quad_p10p8_b1_by(double x, double y)
    {
      return Legendre10(x) * dl8(y);
    }

    static double hdiv_leg_quad_p10p8_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p8_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p9_b1_b(double x, double y)
    {
      return Legendre10(x) * l9(y);
    }

    static double hdiv_leg_quad_p10p9_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p9_b1_bx(double x, double y)
    {
      return Legendre10x(x) * l9(y);
    }

    static double hdiv_leg_quad_p10p9_b1_by(double x, double y)
    {
      return Legendre10(x) * dl9(y);
    }

    static double hdiv_leg_quad_p10p9_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p9_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p10_b1_b(double x, double y)
    {
      return Legendre10(x) * l10(y);
    }

    static double hdiv_leg_quad_p10p10_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p10_b1_bx(double x, double y)
    {
      return Legendre10x(x) * l10(y);
    }

    static double hdiv_leg_quad_p10p10_b1_by(double x, double y)
    {
      return Legendre10(x) * dl10(y);
    }

    static double hdiv_leg_quad_p10p10_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p10_b1_ay(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p11_b1_b(double x, double y)
    {
      return Legendre10(x) * l11(y);
    }

    static double hdiv_leg_quad_p10p11_b1_a(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p11_b1_bx(double x, double y)
    {
      return Legendre10x(x) * l11(y);
    }

    static double hdiv_leg_quad_p10p11_b1_by(double x, double y)
    {
      return Legendre10(x) * dl11(y);
    }

    static double hdiv_leg_quad_p10p11_b1_ax(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p11_b1_ay(double x, double y)
    {
      return 0.0;
    }

    /* BUBBLE ( 0 , 1 ) */

    static double hdiv_leg_quad_p2p0_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p0_b2_a(double x, double y)
    {
      return -l2(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p2p0_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p0_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p0_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p2p0_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p2p1_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p1_b2_a(double x, double y)
    {
      return -l2(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p2p1_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p1_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p1_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p2p1_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p2p2_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p2_b2_a(double x, double y)
    {
      return -l2(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p2p2_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p2_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p2_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p2p2_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p2p3_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p3_b2_a(double x, double y)
    {
      return -l2(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p2p3_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p3_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p3_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p2p3_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p2p4_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p4_b2_a(double x, double y)
    {
      return -l2(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p2p4_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p4_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p4_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p2p4_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p2p5_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p5_b2_a(double x, double y)
    {
      return -l2(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p2p5_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p5_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p5_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p2p5_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p2p6_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p6_b2_a(double x, double y)
    {
      return -l2(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p2p6_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p6_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p6_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p2p6_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p2p7_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p7_b2_a(double x, double y)
    {
      return -l2(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p2p7_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p7_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p7_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p2p7_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p2p8_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p8_b2_a(double x, double y)
    {
      return -l2(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p2p8_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p8_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p8_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p2p8_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p2p9_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p9_b2_a(double x, double y)
    {
      return -l2(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p2p9_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p9_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p9_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p2p9_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p2p10_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p10_b2_a(double x, double y)
    {
      return -l2(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p2p10_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p10_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p2p10_b2_ax(double x, double y)
    {
      return -dl2(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p2p10_b2_ay(double x, double y)
    {
      return -l2(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p3p0_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p0_b2_a(double x, double y)
    {
      return -l3(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p3p0_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p0_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p0_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p3p0_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p3p1_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p1_b2_a(double x, double y)
    {
      return -l3(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p3p1_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p1_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p1_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p3p1_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p3p2_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p2_b2_a(double x, double y)
    {
      return -l3(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p3p2_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p2_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p2_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p3p2_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p3p3_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p3_b2_a(double x, double y)
    {
      return -l3(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p3p3_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p3_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p3_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p3p3_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p3p4_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p4_b2_a(double x, double y)
    {
      return -l3(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p3p4_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p4_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p4_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p3p4_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p3p5_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p5_b2_a(double x, double y)
    {
      return -l3(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p3p5_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p5_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p5_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p3p5_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p3p6_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p6_b2_a(double x, double y)
    {
      return -l3(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p3p6_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p6_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p6_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p3p6_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p3p7_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p7_b2_a(double x, double y)
    {
      return -l3(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p3p7_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p7_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p7_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p3p7_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p3p8_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p8_b2_a(double x, double y)
    {
      return -l3(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p3p8_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p8_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p8_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p3p8_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p3p9_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p9_b2_a(double x, double y)
    {
      return -l3(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p3p9_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p9_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p9_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p3p9_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p3p10_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p10_b2_a(double x, double y)
    {
      return -l3(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p3p10_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p10_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p3p10_b2_ax(double x, double y)
    {
      return -dl3(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p3p10_b2_ay(double x, double y)
    {
      return -l3(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p4p0_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p0_b2_a(double x, double y)
    {
      return -l4(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p4p0_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p0_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p0_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p4p0_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p4p1_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p1_b2_a(double x, double y)
    {
      return -l4(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p4p1_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p1_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p1_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p4p1_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p4p2_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p2_b2_a(double x, double y)
    {
      return -l4(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p4p2_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p2_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p2_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p4p2_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p4p3_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p3_b2_a(double x, double y)
    {
      return -l4(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p4p3_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p3_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p3_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p4p3_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p4p4_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p4_b2_a(double x, double y)
    {
      return -l4(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p4p4_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p4_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p4_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p4p4_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p4p5_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p5_b2_a(double x, double y)
    {
      return -l4(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p4p5_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p5_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p5_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p4p5_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p4p6_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p6_b2_a(double x, double y)
    {
      return -l4(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p4p6_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p6_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p6_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p4p6_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p4p7_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p7_b2_a(double x, double y)
    {
      return -l4(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p4p7_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p7_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p7_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p4p7_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p4p8_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p8_b2_a(double x, double y)
    {
      return -l4(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p4p8_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p8_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p8_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p4p8_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p4p9_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p9_b2_a(double x, double y)
    {
      return -l4(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p4p9_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p9_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p9_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p4p9_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p4p10_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p10_b2_a(double x, double y)
    {
      return -l4(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p4p10_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p10_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p4p10_b2_ax(double x, double y)
    {
      return -dl4(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p4p10_b2_ay(double x, double y)
    {
      return -l4(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p5p0_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p0_b2_a(double x, double y)
    {
      return -l5(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p5p0_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p0_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p0_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p5p0_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p5p1_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p1_b2_a(double x, double y)
    {
      return -l5(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p5p1_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p1_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p1_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p5p1_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p5p2_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p2_b2_a(double x, double y)
    {
      return -l5(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p5p2_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p2_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p2_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p5p2_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p5p3_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p3_b2_a(double x, double y)
    {
      return -l5(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p5p3_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p3_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p3_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p5p3_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p5p4_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p4_b2_a(double x, double y)
    {
      return -l5(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p5p4_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p4_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p4_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p5p4_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p5p5_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p5_b2_a(double x, double y)
    {
      return -l5(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p5p5_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p5_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p5_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p5p5_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p5p6_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p6_b2_a(double x, double y)
    {
      return -l5(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p5p6_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p6_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p6_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p5p6_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p5p7_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p7_b2_a(double x, double y)
    {
      return -l5(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p5p7_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p7_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p7_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p5p7_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p5p8_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p8_b2_a(double x, double y)
    {
      return -l5(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p5p8_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p8_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p8_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p5p8_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p5p9_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p9_b2_a(double x, double y)
    {
      return -l5(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p5p9_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p9_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p9_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p5p9_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p5p10_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p10_b2_a(double x, double y)
    {
      return -l5(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p5p10_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p10_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p5p10_b2_ax(double x, double y)
    {
      return -dl5(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p5p10_b2_ay(double x, double y)
    {
      return -l5(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p6p0_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p0_b2_a(double x, double y)
    {
      return -l6(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p6p0_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p0_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p0_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p6p0_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p6p1_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p1_b2_a(double x, double y)
    {
      return -l6(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p6p1_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p1_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p1_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p6p1_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p6p2_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p2_b2_a(double x, double y)
    {
      return -l6(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p6p2_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p2_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p2_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p6p2_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p6p3_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p3_b2_a(double x, double y)
    {
      return -l6(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p6p3_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p3_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p3_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p6p3_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p6p4_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p4_b2_a(double x, double y)
    {
      return -l6(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p6p4_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p4_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p4_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p6p4_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p6p5_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p5_b2_a(double x, double y)
    {
      return -l6(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p6p5_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p5_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p5_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p6p5_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p6p6_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p6_b2_a(double x, double y)
    {
      return -l6(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p6p6_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p6_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p6_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p6p6_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p6p7_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p7_b2_a(double x, double y)
    {
      return -l6(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p6p7_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p7_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p7_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p6p7_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p6p8_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p8_b2_a(double x, double y)
    {
      return -l6(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p6p8_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p8_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p8_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p6p8_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p6p9_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p9_b2_a(double x, double y)
    {
      return -l6(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p6p9_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p9_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p9_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p6p9_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p6p10_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p10_b2_a(double x, double y)
    {
      return -l6(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p6p10_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p10_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p6p10_b2_ax(double x, double y)
    {
      return -dl6(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p6p10_b2_ay(double x, double y)
    {
      return -l6(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p7p0_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p0_b2_a(double x, double y)
    {
      return -l7(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p7p0_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p0_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p0_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p7p0_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p7p1_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p1_b2_a(double x, double y)
    {
      return -l7(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p7p1_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p1_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p1_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p7p1_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p7p2_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p2_b2_a(double x, double y)
    {
      return -l7(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p7p2_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p2_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p2_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p7p2_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p7p3_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p3_b2_a(double x, double y)
    {
      return -l7(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p7p3_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p3_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p3_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p7p3_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p7p4_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p4_b2_a(double x, double y)
    {
      return -l7(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p7p4_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p4_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p4_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p7p4_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p7p5_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p5_b2_a(double x, double y)
    {
      return -l7(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p7p5_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p5_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p5_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p7p5_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p7p6_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p6_b2_a(double x, double y)
    {
      return -l7(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p7p6_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p6_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p6_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p7p6_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p7p7_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p7_b2_a(double x, double y)
    {
      return -l7(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p7p7_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p7_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p7_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p7p7_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p7p8_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p8_b2_a(double x, double y)
    {
      return -l7(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p7p8_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p8_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p8_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p7p8_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p7p9_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p9_b2_a(double x, double y)
    {
      return -l7(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p7p9_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p9_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p9_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p7p9_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p7p10_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p10_b2_a(double x, double y)
    {
      return -l7(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p7p10_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p10_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p7p10_b2_ax(double x, double y)
    {
      return -dl7(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p7p10_b2_ay(double x, double y)
    {
      return -l7(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p8p0_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p0_b2_a(double x, double y)
    {
      return -l8(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p8p0_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p0_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p0_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p8p0_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p8p1_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p1_b2_a(double x, double y)
    {
      return -l8(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p8p1_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p1_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p1_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p8p1_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p8p2_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p2_b2_a(double x, double y)
    {
      return -l8(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p8p2_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p2_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p2_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p8p2_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p8p3_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p3_b2_a(double x, double y)
    {
      return -l8(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p8p3_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p3_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p3_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p8p3_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p8p4_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p4_b2_a(double x, double y)
    {
      return -l8(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p8p4_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p4_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p4_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p8p4_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p8p5_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p5_b2_a(double x, double y)
    {
      return -l8(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p8p5_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p5_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p5_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p8p5_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p8p6_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p6_b2_a(double x, double y)
    {
      return -l8(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p8p6_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p6_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p6_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p8p6_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p8p7_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p7_b2_a(double x, double y)
    {
      return -l8(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p8p7_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p7_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p7_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p8p7_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p8p8_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p8_b2_a(double x, double y)
    {
      return -l8(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p8p8_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p8_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p8_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p8p8_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p8p9_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p9_b2_a(double x, double y)
    {
      return -l8(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p8p9_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p9_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p9_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p8p9_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p8p10_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p10_b2_a(double x, double y)
    {
      return -l8(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p8p10_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p10_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p8p10_b2_ax(double x, double y)
    {
      return -dl8(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p8p10_b2_ay(double x, double y)
    {
      return -l8(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p9p0_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p0_b2_a(double x, double y)
    {
      return -l9(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p9p0_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p0_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p0_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p9p0_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p9p1_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p1_b2_a(double x, double y)
    {
      return -l9(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p9p1_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p1_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p1_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p9p1_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p9p2_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p2_b2_a(double x, double y)
    {
      return -l9(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p9p2_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p2_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p2_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p9p2_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p9p3_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p3_b2_a(double x, double y)
    {
      return -l9(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p9p3_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p3_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p3_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p9p3_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p9p4_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p4_b2_a(double x, double y)
    {
      return -l9(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p9p4_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p4_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p4_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p9p4_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p9p5_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p5_b2_a(double x, double y)
    {
      return -l9(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p9p5_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p5_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p5_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p9p5_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p9p6_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p6_b2_a(double x, double y)
    {
      return -l9(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p9p6_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p6_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p6_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p9p6_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p9p7_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p7_b2_a(double x, double y)
    {
      return -l9(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p9p7_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p7_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p7_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p9p7_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p9p8_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p8_b2_a(double x, double y)
    {
      return -l9(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p9p8_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p8_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p8_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p9p8_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p9p9_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p9_b2_a(double x, double y)
    {
      return -l9(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p9p9_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p9_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p9_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p9p9_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p9p10_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p10_b2_a(double x, double y)
    {
      return -l9(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p9p10_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p10_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p9p10_b2_ax(double x, double y)
    {
      return -dl9(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p9p10_b2_ay(double x, double y)
    {
      return -l9(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p10p0_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p0_b2_a(double x, double y)
    {
      return -l10(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p10p0_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p0_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p0_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p10p0_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p10p1_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p1_b2_a(double x, double y)
    {
      return -l10(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p10p1_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p1_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p1_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p10p1_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p10p2_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p2_b2_a(double x, double y)
    {
      return -l10(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p10p2_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p2_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p2_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p10p2_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p10p3_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p3_b2_a(double x, double y)
    {
      return -l10(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p10p3_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p3_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p3_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p10p3_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p10p4_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p4_b2_a(double x, double y)
    {
      return -l10(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p10p4_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p4_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p4_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p10p4_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p10p5_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p5_b2_a(double x, double y)
    {
      return -l10(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p10p5_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p5_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p5_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p10p5_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p10p6_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p6_b2_a(double x, double y)
    {
      return -l10(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p10p6_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p6_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p6_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p10p6_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p10p7_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p7_b2_a(double x, double y)
    {
      return -l10(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p10p7_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p7_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p7_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p10p7_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p10p8_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p8_b2_a(double x, double y)
    {
      return -l10(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p10p8_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p8_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p8_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p10p8_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p10p9_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p9_b2_a(double x, double y)
    {
      return -l10(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p10p9_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p9_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p9_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p10p9_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p10p10_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p10_b2_a(double x, double y)
    {
      return -l10(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p10p10_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p10_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p10p10_b2_ax(double x, double y)
    {
      return -dl10(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p10p10_b2_ay(double x, double y)
    {
      return -l10(x) * Legendre10x(y);
    }

    static double hdiv_leg_quad_p11p0_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p0_b2_a(double x, double y)
    {
      return -l11(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p11p0_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p0_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p0_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre0(y);
    }

    static double hdiv_leg_quad_p11p0_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre0x(y);
    }

    static double hdiv_leg_quad_p11p1_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p1_b2_a(double x, double y)
    {
      return -l11(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p11p1_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p1_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p1_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre1(y);
    }

    static double hdiv_leg_quad_p11p1_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre1x(y);
    }

    static double hdiv_leg_quad_p11p2_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p2_b2_a(double x, double y)
    {
      return -l11(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p11p2_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p2_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p2_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre2(y);
    }

    static double hdiv_leg_quad_p11p2_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre2x(y);
    }

    static double hdiv_leg_quad_p11p3_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p3_b2_a(double x, double y)
    {
      return -l11(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p11p3_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p3_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p3_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre3(y);
    }

    static double hdiv_leg_quad_p11p3_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre3x(y);
    }

    static double hdiv_leg_quad_p11p4_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p4_b2_a(double x, double y)
    {
      return -l11(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p11p4_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p4_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p4_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre4(y);
    }

    static double hdiv_leg_quad_p11p4_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre4x(y);
    }

    static double hdiv_leg_quad_p11p5_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p5_b2_a(double x, double y)
    {
      return -l11(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p11p5_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p5_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p5_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre5(y);
    }

    static double hdiv_leg_quad_p11p5_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre5x(y);
    }

    static double hdiv_leg_quad_p11p6_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p6_b2_a(double x, double y)
    {
      return -l11(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p11p6_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p6_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p6_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre6(y);
    }

    static double hdiv_leg_quad_p11p6_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre6x(y);
    }

    static double hdiv_leg_quad_p11p7_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p7_b2_a(double x, double y)
    {
      return -l11(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p11p7_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p7_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p7_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre7(y);
    }

    static double hdiv_leg_quad_p11p7_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre7x(y);
    }

    static double hdiv_leg_quad_p11p8_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p8_b2_a(double x, double y)
    {
      return -l11(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p11p8_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p8_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p8_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre8(y);
    }

    static double hdiv_leg_quad_p11p8_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre8x(y);
    }

    static double hdiv_leg_quad_p11p9_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p9_b2_a(double x, double y)
    {
      return -l11(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p11p9_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p9_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p9_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre9(y);
    }

    static double hdiv_leg_quad_p11p9_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre9x(y);
    }

    static double hdiv_leg_quad_p11p10_b2_b(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p10_b2_a(double x, double y)
    {
      return -l11(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p11p10_b2_bx(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p10_b2_by(double x, double y)
    {
      return 0.0;
    }

    static double hdiv_leg_quad_p11p10_b2_ax(double x, double y)
    {
      return -dl11(x) * Legendre10(y);
    }

    static double hdiv_leg_quad_p11p10_b2_ay(double x, double y)
    {
      return -l11(x) * Legendre10x(y);
    }

    static Shapeset::shape_fn_t hdiv_leg_quad_fn_a[] =
    {
      hdiv_leg_quad_p0_e1_a_0, hdiv_leg_quad_p0_e1_a_1, hdiv_leg_quad_p0_e2_a_0, hdiv_leg_quad_p0_e2_a_1,  hdiv_leg_quad_p0_e3_a, hdiv_leg_quad_p0_e3_a, hdiv_leg_quad_p0_e4_a, hdiv_leg_quad_p0_e4_a,
      hdiv_leg_quad_p1_e1_a, hdiv_leg_quad_p1_e1_a, hdiv_leg_quad_p1_e2_a, hdiv_leg_quad_p1_e2_a, hdiv_leg_quad_p1_e3_a, hdiv_leg_quad_p1_e3_a, hdiv_leg_quad_p1_e4_a, hdiv_leg_quad_p1_e4_a,
      hdiv_leg_quad_p2_e1_a_0, hdiv_leg_quad_p2_e1_a_1, hdiv_leg_quad_p2_e2_a_0, hdiv_leg_quad_p2_e2_a_1, hdiv_leg_quad_p2_e3_a, hdiv_leg_quad_p2_e3_a, hdiv_leg_quad_p2_e4_a, hdiv_leg_quad_p2_e4_a,
      hdiv_leg_quad_p3_e1_a, hdiv_leg_quad_p3_e1_a, hdiv_leg_quad_p3_e2_a, hdiv_leg_quad_p3_e2_a, hdiv_leg_quad_p3_e3_a, hdiv_leg_quad_p3_e3_a, hdiv_leg_quad_p3_e4_a, hdiv_leg_quad_p3_e4_a,
      hdiv_leg_quad_p4_e1_a_0, hdiv_leg_quad_p4_e1_a_1, hdiv_leg_quad_p4_e2_a_0, hdiv_leg_quad_p4_e2_a_1, hdiv_leg_quad_p4_e3_a, hdiv_leg_quad_p4_e3_a, hdiv_leg_quad_p4_e4_a, hdiv_leg_quad_p4_e4_a,
      hdiv_leg_quad_p5_e1_a, hdiv_leg_quad_p5_e1_a, hdiv_leg_quad_p5_e2_a, hdiv_leg_quad_p5_e2_a, hdiv_leg_quad_p5_e3_a, hdiv_leg_quad_p5_e3_a, hdiv_leg_quad_p5_e4_a, hdiv_leg_quad_p5_e4_a,
      hdiv_leg_quad_p6_e1_a_0, hdiv_leg_quad_p6_e1_a_1, hdiv_leg_quad_p6_e2_a_0, hdiv_leg_quad_p6_e2_a_1, hdiv_leg_quad_p6_e3_a, hdiv_leg_quad_p6_e3_a, hdiv_leg_quad_p6_e4_a, hdiv_leg_quad_p6_e4_a,
      hdiv_leg_quad_p7_e1_a, hdiv_leg_quad_p7_e1_a, hdiv_leg_quad_p7_e2_a, hdiv_leg_quad_p7_e2_a, hdiv_leg_quad_p7_e3_a, hdiv_leg_quad_p7_e3_a, hdiv_leg_quad_p7_e4_a, hdiv_leg_quad_p7_e4_a,
      hdiv_leg_quad_p8_e1_a_0, hdiv_leg_quad_p8_e1_a_1, hdiv_leg_quad_p8_e2_a_0, hdiv_leg_quad_p8_e2_a_1, hdiv_leg_quad_p8_e3_a, hdiv_leg_quad_p8_e3_a, hdiv_leg_quad_p8_e4_a, hdiv_leg_quad_p8_e4_a,
      hdiv_leg_quad_p9_e1_a, hdiv_leg_quad_p9_e1_a, hdiv_leg_quad_p9_e2_a, hdiv_leg_quad_p9_e2_a, hdiv_leg_quad_p9_e3_a, hdiv_leg_quad_p9_e3_a, hdiv_leg_quad_p9_e4_a, hdiv_leg_quad_p9_e4_a,
      hdiv_leg_quad_p10_e1_a_0, hdiv_leg_quad_p10_e1_a_1, hdiv_leg_quad_p10_e2_a_0, hdiv_leg_quad_p10_e2_a_1, hdiv_leg_quad_p10_e3_a, hdiv_leg_quad_p10_e3_a, hdiv_leg_quad_p10_e4_a, hdiv_leg_quad_p10_e4_a,

      hdiv_leg_quad_p0p2_b1_a,   hdiv_leg_quad_p0p3_b1_a,   hdiv_leg_quad_p0p4_b1_a,   hdiv_leg_quad_p0p5_b1_a,   hdiv_leg_quad_p0p6_b1_a,   hdiv_leg_quad_p0p7_b1_a,   hdiv_leg_quad_p0p8_b1_a,   hdiv_leg_quad_p0p9_b1_a,   hdiv_leg_quad_p0p10_b1_a,   hdiv_leg_quad_p0p11_b1_a,   hdiv_leg_quad_p1p2_b1_a,   hdiv_leg_quad_p1p3_b1_a,   hdiv_leg_quad_p1p4_b1_a,   hdiv_leg_quad_p1p5_b1_a,   hdiv_leg_quad_p1p6_b1_a,   hdiv_leg_quad_p1p7_b1_a,   hdiv_leg_quad_p1p8_b1_a,   hdiv_leg_quad_p1p9_b1_a,   hdiv_leg_quad_p1p10_b1_a,   hdiv_leg_quad_p1p11_b1_a,   hdiv_leg_quad_p2p2_b1_a,   hdiv_leg_quad_p2p3_b1_a,   hdiv_leg_quad_p2p4_b1_a,   hdiv_leg_quad_p2p5_b1_a,   hdiv_leg_quad_p2p6_b1_a,   hdiv_leg_quad_p2p7_b1_a,   hdiv_leg_quad_p2p8_b1_a,   hdiv_leg_quad_p2p9_b1_a,   hdiv_leg_quad_p2p10_b1_a,   hdiv_leg_quad_p2p11_b1_a,   hdiv_leg_quad_p3p2_b1_a,   hdiv_leg_quad_p3p3_b1_a,   hdiv_leg_quad_p3p4_b1_a,   hdiv_leg_quad_p3p5_b1_a,   hdiv_leg_quad_p3p6_b1_a,   hdiv_leg_quad_p3p7_b1_a,   hdiv_leg_quad_p3p8_b1_a,   hdiv_leg_quad_p3p9_b1_a,   hdiv_leg_quad_p3p10_b1_a,   hdiv_leg_quad_p3p11_b1_a,   hdiv_leg_quad_p4p2_b1_a,   hdiv_leg_quad_p4p3_b1_a,   hdiv_leg_quad_p4p4_b1_a,   hdiv_leg_quad_p4p5_b1_a,   hdiv_leg_quad_p4p6_b1_a,   hdiv_leg_quad_p4p7_b1_a,   hdiv_leg_quad_p4p8_b1_a,   hdiv_leg_quad_p4p9_b1_a,   hdiv_leg_quad_p4p10_b1_a,   hdiv_leg_quad_p4p11_b1_a,   hdiv_leg_quad_p5p2_b1_a,   hdiv_leg_quad_p5p3_b1_a,   hdiv_leg_quad_p5p4_b1_a,   hdiv_leg_quad_p5p5_b1_a,   hdiv_leg_quad_p5p6_b1_a,   hdiv_leg_quad_p5p7_b1_a,   hdiv_leg_quad_p5p8_b1_a,   hdiv_leg_quad_p5p9_b1_a,   hdiv_leg_quad_p5p10_b1_a,   hdiv_leg_quad_p5p11_b1_a,   hdiv_leg_quad_p6p2_b1_a,   hdiv_leg_quad_p6p3_b1_a,   hdiv_leg_quad_p6p4_b1_a,   hdiv_leg_quad_p6p5_b1_a,   hdiv_leg_quad_p6p6_b1_a,   hdiv_leg_quad_p6p7_b1_a,   hdiv_leg_quad_p6p8_b1_a,   hdiv_leg_quad_p6p9_b1_a,   hdiv_leg_quad_p6p10_b1_a,   hdiv_leg_quad_p6p11_b1_a,   hdiv_leg_quad_p7p2_b1_a,   hdiv_leg_quad_p7p3_b1_a,   hdiv_leg_quad_p7p4_b1_a,   hdiv_leg_quad_p7p5_b1_a,   hdiv_leg_quad_p7p6_b1_a,   hdiv_leg_quad_p7p7_b1_a,   hdiv_leg_quad_p7p8_b1_a,   hdiv_leg_quad_p7p9_b1_a,   hdiv_leg_quad_p7p10_b1_a,   hdiv_leg_quad_p7p11_b1_a,   hdiv_leg_quad_p8p2_b1_a,   hdiv_leg_quad_p8p3_b1_a,   hdiv_leg_quad_p8p4_b1_a,   hdiv_leg_quad_p8p5_b1_a,   hdiv_leg_quad_p8p6_b1_a,   hdiv_leg_quad_p8p7_b1_a,   hdiv_leg_quad_p8p8_b1_a,   hdiv_leg_quad_p8p9_b1_a,   hdiv_leg_quad_p8p10_b1_a,   hdiv_leg_quad_p8p11_b1_a,   hdiv_leg_quad_p9p2_b1_a,   hdiv_leg_quad_p9p3_b1_a,   hdiv_leg_quad_p9p4_b1_a,   hdiv_leg_quad_p9p5_b1_a,   hdiv_leg_quad_p9p6_b1_a,   hdiv_leg_quad_p9p7_b1_a,   hdiv_leg_quad_p9p8_b1_a,   hdiv_leg_quad_p9p9_b1_a,   hdiv_leg_quad_p9p10_b1_a,   hdiv_leg_quad_p9p11_b1_a,   hdiv_leg_quad_p10p2_b1_a,   hdiv_leg_quad_p10p3_b1_a,   hdiv_leg_quad_p10p4_b1_a,   hdiv_leg_quad_p10p5_b1_a,   hdiv_leg_quad_p10p6_b1_a,   hdiv_leg_quad_p10p7_b1_a,   hdiv_leg_quad_p10p8_b1_a,   hdiv_leg_quad_p10p9_b1_a,   hdiv_leg_quad_p10p10_b1_a,   hdiv_leg_quad_p10p11_b1_a,   hdiv_leg_quad_p2p0_b2_a,   hdiv_leg_quad_p2p1_b2_a,   hdiv_leg_quad_p2p2_b2_a,   hdiv_leg_quad_p2p3_b2_a,   hdiv_leg_quad_p2p4_b2_a,   hdiv_leg_quad_p2p5_b2_a,   hdiv_leg_quad_p2p6_b2_a,   hdiv_leg_quad_p2p7_b2_a,   hdiv_leg_quad_p2p8_b2_a,   hdiv_leg_quad_p2p9_b2_a,   hdiv_leg_quad_p2p10_b2_a,   hdiv_leg_quad_p3p0_b2_a,   hdiv_leg_quad_p3p1_b2_a,   hdiv_leg_quad_p3p2_b2_a,   hdiv_leg_quad_p3p3_b2_a,   hdiv_leg_quad_p3p4_b2_a,   hdiv_leg_quad_p3p5_b2_a,   hdiv_leg_quad_p3p6_b2_a,   hdiv_leg_quad_p3p7_b2_a,   hdiv_leg_quad_p3p8_b2_a,   hdiv_leg_quad_p3p9_b2_a,   hdiv_leg_quad_p3p10_b2_a,   hdiv_leg_quad_p4p0_b2_a,   hdiv_leg_quad_p4p1_b2_a,   hdiv_leg_quad_p4p2_b2_a,   hdiv_leg_quad_p4p3_b2_a,   hdiv_leg_quad_p4p4_b2_a,   hdiv_leg_quad_p4p5_b2_a,   hdiv_leg_quad_p4p6_b2_a,   hdiv_leg_quad_p4p7_b2_a,   hdiv_leg_quad_p4p8_b2_a,   hdiv_leg_quad_p4p9_b2_a,   hdiv_leg_quad_p4p10_b2_a,   hdiv_leg_quad_p5p0_b2_a,   hdiv_leg_quad_p5p1_b2_a,   hdiv_leg_quad_p5p2_b2_a,   hdiv_leg_quad_p5p3_b2_a,   hdiv_leg_quad_p5p4_b2_a,   hdiv_leg_quad_p5p5_b2_a,   hdiv_leg_quad_p5p6_b2_a,   hdiv_leg_quad_p5p7_b2_a,   hdiv_leg_quad_p5p8_b2_a,   hdiv_leg_quad_p5p9_b2_a,   hdiv_leg_quad_p5p10_b2_a,   hdiv_leg_quad_p6p0_b2_a,   hdiv_leg_quad_p6p1_b2_a,   hdiv_leg_quad_p6p2_b2_a,   hdiv_leg_quad_p6p3_b2_a,   hdiv_leg_quad_p6p4_b2_a,   hdiv_leg_quad_p6p5_b2_a,   hdiv_leg_quad_p6p6_b2_a,   hdiv_leg_quad_p6p7_b2_a,   hdiv_leg_quad_p6p8_b2_a,   hdiv_leg_quad_p6p9_b2_a,   hdiv_leg_quad_p6p10_b2_a,   hdiv_leg_quad_p7p0_b2_a,   hdiv_leg_quad_p7p1_b2_a,   hdiv_leg_quad_p7p2_b2_a,   hdiv_leg_quad_p7p3_b2_a,   hdiv_leg_quad_p7p4_b2_a,   hdiv_leg_quad_p7p5_b2_a,   hdiv_leg_quad_p7p6_b2_a,   hdiv_leg_quad_p7p7_b2_a,   hdiv_leg_quad_p7p8_b2_a,   hdiv_leg_quad_p7p9_b2_a,   hdiv_leg_quad_p7p10_b2_a,   hdiv_leg_quad_p8p0_b2_a,   hdiv_leg_quad_p8p1_b2_a,   hdiv_leg_quad_p8p2_b2_a,   hdiv_leg_quad_p8p3_b2_a,   hdiv_leg_quad_p8p4_b2_a,   hdiv_leg_quad_p8p5_b2_a,   hdiv_leg_quad_p8p6_b2_a,   hdiv_leg_quad_p8p7_b2_a,   hdiv_leg_quad_p8p8_b2_a,   hdiv_leg_quad_p8p9_b2_a,   hdiv_leg_quad_p8p10_b2_a,   hdiv_leg_quad_p9p0_b2_a,   hdiv_leg_quad_p9p1_b2_a,   hdiv_leg_quad_p9p2_b2_a,   hdiv_leg_quad_p9p3_b2_a,   hdiv_leg_quad_p9p4_b2_a,   hdiv_leg_quad_p9p5_b2_a,   hdiv_leg_quad_p9p6_b2_a,   hdiv_leg_quad_p9p7_b2_a,   hdiv_leg_quad_p9p8_b2_a,   hdiv_leg_quad_p9p9_b2_a,   hdiv_leg_quad_p9p10_b2_a,   hdiv_leg_quad_p10p0_b2_a,   hdiv_leg_quad_p10p1_b2_a,   hdiv_leg_quad_p10p2_b2_a,   hdiv_leg_quad_p10p3_b2_a,   hdiv_leg_quad_p10p4_b2_a,   hdiv_leg_quad_p10p5_b2_a,   hdiv_leg_quad_p10p6_b2_a,   hdiv_leg_quad_p10p7_b2_a,   hdiv_leg_quad_p10p8_b2_a,   hdiv_leg_quad_p10p9_b2_a,   hdiv_leg_quad_p10p10_b2_a,   hdiv_leg_quad_p11p0_b2_a,   hdiv_leg_quad_p11p1_b2_a,   hdiv_leg_quad_p11p2_b2_a,   hdiv_leg_quad_p11p3_b2_a,   hdiv_leg_quad_p11p4_b2_a,   hdiv_leg_quad_p11p5_b2_a,   hdiv_leg_quad_p11p6_b2_a,   hdiv_leg_quad_p11p7_b2_a,   hdiv_leg_quad_p11p8_b2_a,   hdiv_leg_quad_p11p9_b2_a,   hdiv_leg_quad_p11p10_b2_a, };

    static Shapeset::shape_fn_t hdiv_leg_quad_fn_b[] =
    {
      hdiv_leg_quad_p0_e1_b, hdiv_leg_quad_p0_e1_b, hdiv_leg_quad_p0_e2_b, hdiv_leg_quad_p0_e2_b, hdiv_leg_quad_p0_e3_b_0, hdiv_leg_quad_p0_e3_b_1, hdiv_leg_quad_p0_e4_b_0, hdiv_leg_quad_p0_e4_b_1,
      hdiv_leg_quad_p1_e1_b, hdiv_leg_quad_p1_e1_b, hdiv_leg_quad_p1_e2_b, hdiv_leg_quad_p1_e2_b, hdiv_leg_quad_p1_e3_b, hdiv_leg_quad_p1_e3_b, hdiv_leg_quad_p1_e4_b, hdiv_leg_quad_p1_e4_b,
      hdiv_leg_quad_p2_e1_b, hdiv_leg_quad_p2_e1_b, hdiv_leg_quad_p2_e2_b, hdiv_leg_quad_p2_e2_b, hdiv_leg_quad_p2_e3_b_0, hdiv_leg_quad_p2_e3_b_1, hdiv_leg_quad_p2_e4_b_0, hdiv_leg_quad_p2_e4_b_1,
      hdiv_leg_quad_p3_e1_b, hdiv_leg_quad_p3_e1_b, hdiv_leg_quad_p3_e2_b, hdiv_leg_quad_p3_e2_b, hdiv_leg_quad_p3_e3_b, hdiv_leg_quad_p3_e3_b, hdiv_leg_quad_p3_e4_b, hdiv_leg_quad_p3_e4_b,
      hdiv_leg_quad_p4_e1_b, hdiv_leg_quad_p4_e1_b, hdiv_leg_quad_p4_e2_b, hdiv_leg_quad_p4_e2_b, hdiv_leg_quad_p4_e3_b_0, hdiv_leg_quad_p4_e3_b_1, hdiv_leg_quad_p4_e4_b_0, hdiv_leg_quad_p4_e4_b_1,
      hdiv_leg_quad_p5_e1_b, hdiv_leg_quad_p5_e1_b, hdiv_leg_quad_p5_e2_b, hdiv_leg_quad_p5_e2_b, hdiv_leg_quad_p5_e3_b, hdiv_leg_quad_p5_e3_b, hdiv_leg_quad_p5_e4_b, hdiv_leg_quad_p5_e4_b,
      hdiv_leg_quad_p6_e1_b, hdiv_leg_quad_p6_e1_b, hdiv_leg_quad_p6_e2_b, hdiv_leg_quad_p6_e2_b, hdiv_leg_quad_p6_e3_b_0, hdiv_leg_quad_p6_e3_b_1, hdiv_leg_quad_p6_e4_b_0, hdiv_leg_quad_p6_e4_b_1,
      hdiv_leg_quad_p7_e1_b, hdiv_leg_quad_p7_e1_b, hdiv_leg_quad_p7_e2_b, hdiv_leg_quad_p7_e2_b, hdiv_leg_quad_p7_e3_b, hdiv_leg_quad_p7_e3_b, hdiv_leg_quad_p7_e4_b, hdiv_leg_quad_p7_e4_b,
      hdiv_leg_quad_p8_e1_b, hdiv_leg_quad_p8_e1_b, hdiv_leg_quad_p8_e2_b, hdiv_leg_quad_p8_e2_b, hdiv_leg_quad_p8_e3_b_0, hdiv_leg_quad_p8_e3_b_1, hdiv_leg_quad_p8_e4_b_0, hdiv_leg_quad_p8_e4_b_1,
      hdiv_leg_quad_p9_e1_b, hdiv_leg_quad_p9_e1_b, hdiv_leg_quad_p9_e2_b, hdiv_leg_quad_p9_e2_b, hdiv_leg_quad_p9_e3_b, hdiv_leg_quad_p9_e3_b, hdiv_leg_quad_p9_e4_b, hdiv_leg_quad_p9_e4_b,
      hdiv_leg_quad_p10_e1_b, hdiv_leg_quad_p10_e1_b, hdiv_leg_quad_p10_e2_b, hdiv_leg_quad_p10_e2_b, hdiv_leg_quad_p10_e3_b_0, hdiv_leg_quad_p10_e3_b_1, hdiv_leg_quad_p10_e4_b_0, hdiv_leg_quad_p10_e4_b_1,

      hdiv_leg_quad_p0p2_b1_b,   hdiv_leg_quad_p0p3_b1_b,   hdiv_leg_quad_p0p4_b1_b,   hdiv_leg_quad_p0p5_b1_b,   hdiv_leg_quad_p0p6_b1_b,   hdiv_leg_quad_p0p7_b1_b,   hdiv_leg_quad_p0p8_b1_b,   hdiv_leg_quad_p0p9_b1_b,   hdiv_leg_quad_p0p10_b1_b,   hdiv_leg_quad_p0p11_b1_b,   hdiv_leg_quad_p1p2_b1_b,   hdiv_leg_quad_p1p3_b1_b,   hdiv_leg_quad_p1p4_b1_b,   hdiv_leg_quad_p1p5_b1_b,   hdiv_leg_quad_p1p6_b1_b,   hdiv_leg_quad_p1p7_b1_b,   hdiv_leg_quad_p1p8_b1_b,   hdiv_leg_quad_p1p9_b1_b,   hdiv_leg_quad_p1p10_b1_b,   hdiv_leg_quad_p1p11_b1_b,   hdiv_leg_quad_p2p2_b1_b,   hdiv_leg_quad_p2p3_b1_b,   hdiv_leg_quad_p2p4_b1_b,   hdiv_leg_quad_p2p5_b1_b,   hdiv_leg_quad_p2p6_b1_b,   hdiv_leg_quad_p2p7_b1_b,   hdiv_leg_quad_p2p8_b1_b,   hdiv_leg_quad_p2p9_b1_b,   hdiv_leg_quad_p2p10_b1_b,   hdiv_leg_quad_p2p11_b1_b,   hdiv_leg_quad_p3p2_b1_b,   hdiv_leg_quad_p3p3_b1_b,   hdiv_leg_quad_p3p4_b1_b,   hdiv_leg_quad_p3p5_b1_b,   hdiv_leg_quad_p3p6_b1_b,   hdiv_leg_quad_p3p7_b1_b,   hdiv_leg_quad_p3p8_b1_b,   hdiv_leg_quad_p3p9_b1_b,   hdiv_leg_quad_p3p10_b1_b,   hdiv_leg_quad_p3p11_b1_b,   hdiv_leg_quad_p4p2_b1_b,   hdiv_leg_quad_p4p3_b1_b,   hdiv_leg_quad_p4p4_b1_b,   hdiv_leg_quad_p4p5_b1_b,   hdiv_leg_quad_p4p6_b1_b,   hdiv_leg_quad_p4p7_b1_b,   hdiv_leg_quad_p4p8_b1_b,   hdiv_leg_quad_p4p9_b1_b,   hdiv_leg_quad_p4p10_b1_b,   hdiv_leg_quad_p4p11_b1_b,   hdiv_leg_quad_p5p2_b1_b,   hdiv_leg_quad_p5p3_b1_b,   hdiv_leg_quad_p5p4_b1_b,   hdiv_leg_quad_p5p5_b1_b,   hdiv_leg_quad_p5p6_b1_b,   hdiv_leg_quad_p5p7_b1_b,   hdiv_leg_quad_p5p8_b1_b,   hdiv_leg_quad_p5p9_b1_b,   hdiv_leg_quad_p5p10_b1_b,   hdiv_leg_quad_p5p11_b1_b,   hdiv_leg_quad_p6p2_b1_b,   hdiv_leg_quad_p6p3_b1_b,   hdiv_leg_quad_p6p4_b1_b,   hdiv_leg_quad_p6p5_b1_b,   hdiv_leg_quad_p6p6_b1_b,   hdiv_leg_quad_p6p7_b1_b,   hdiv_leg_quad_p6p8_b1_b,   hdiv_leg_quad_p6p9_b1_b,   hdiv_leg_quad_p6p10_b1_b,   hdiv_leg_quad_p6p11_b1_b,   hdiv_leg_quad_p7p2_b1_b,   hdiv_leg_quad_p7p3_b1_b,   hdiv_leg_quad_p7p4_b1_b,   hdiv_leg_quad_p7p5_b1_b,   hdiv_leg_quad_p7p6_b1_b,   hdiv_leg_quad_p7p7_b1_b,   hdiv_leg_quad_p7p8_b1_b,   hdiv_leg_quad_p7p9_b1_b,   hdiv_leg_quad_p7p10_b1_b,   hdiv_leg_quad_p7p11_b1_b,   hdiv_leg_quad_p8p2_b1_b,   hdiv_leg_quad_p8p3_b1_b,   hdiv_leg_quad_p8p4_b1_b,   hdiv_leg_quad_p8p5_b1_b,   hdiv_leg_quad_p8p6_b1_b,   hdiv_leg_quad_p8p7_b1_b,   hdiv_leg_quad_p8p8_b1_b,   hdiv_leg_quad_p8p9_b1_b,   hdiv_leg_quad_p8p10_b1_b,   hdiv_leg_quad_p8p11_b1_b,   hdiv_leg_quad_p9p2_b1_b,   hdiv_leg_quad_p9p3_b1_b,   hdiv_leg_quad_p9p4_b1_b,   hdiv_leg_quad_p9p5_b1_b,   hdiv_leg_quad_p9p6_b1_b,   hdiv_leg_quad_p9p7_b1_b,   hdiv_leg_quad_p9p8_b1_b,   hdiv_leg_quad_p9p9_b1_b,   hdiv_leg_quad_p9p10_b1_b,   hdiv_leg_quad_p9p11_b1_b,   hdiv_leg_quad_p10p2_b1_b,   hdiv_leg_quad_p10p3_b1_b,   hdiv_leg_quad_p10p4_b1_b,   hdiv_leg_quad_p10p5_b1_b,   hdiv_leg_quad_p10p6_b1_b,   hdiv_leg_quad_p10p7_b1_b,   hdiv_leg_quad_p10p8_b1_b,   hdiv_leg_quad_p10p9_b1_b,   hdiv_leg_quad_p10p10_b1_b,   hdiv_leg_quad_p10p11_b1_b,   hdiv_leg_quad_p2p0_b2_b,   hdiv_leg_quad_p2p1_b2_b,   hdiv_leg_quad_p2p2_b2_b,   hdiv_leg_quad_p2p3_b2_b,   hdiv_leg_quad_p2p4_b2_b,   hdiv_leg_quad_p2p5_b2_b,   hdiv_leg_quad_p2p6_b2_b,   hdiv_leg_quad_p2p7_b2_b,   hdiv_leg_quad_p2p8_b2_b,   hdiv_leg_quad_p2p9_b2_b,   hdiv_leg_quad_p2p10_b2_b,   hdiv_leg_quad_p3p0_b2_b,   hdiv_leg_quad_p3p1_b2_b,   hdiv_leg_quad_p3p2_b2_b,   hdiv_leg_quad_p3p3_b2_b,   hdiv_leg_quad_p3p4_b2_b,   hdiv_leg_quad_p3p5_b2_b,   hdiv_leg_quad_p3p6_b2_b,   hdiv_leg_quad_p3p7_b2_b,   hdiv_leg_quad_p3p8_b2_b,   hdiv_leg_quad_p3p9_b2_b,   hdiv_leg_quad_p3p10_b2_b,   hdiv_leg_quad_p4p0_b2_b,   hdiv_leg_quad_p4p1_b2_b,   hdiv_leg_quad_p4p2_b2_b,   hdiv_leg_quad_p4p3_b2_b,   hdiv_leg_quad_p4p4_b2_b,   hdiv_leg_quad_p4p5_b2_b,   hdiv_leg_quad_p4p6_b2_b,   hdiv_leg_quad_p4p7_b2_b,   hdiv_leg_quad_p4p8_b2_b,   hdiv_leg_quad_p4p9_b2_b,   hdiv_leg_quad_p4p10_b2_b,   hdiv_leg_quad_p5p0_b2_b,   hdiv_leg_quad_p5p1_b2_b,   hdiv_leg_quad_p5p2_b2_b,   hdiv_leg_quad_p5p3_b2_b,   hdiv_leg_quad_p5p4_b2_b,   hdiv_leg_quad_p5p5_b2_b,   hdiv_leg_quad_p5p6_b2_b,   hdiv_leg_quad_p5p7_b2_b,   hdiv_leg_quad_p5p8_b2_b,   hdiv_leg_quad_p5p9_b2_b,   hdiv_leg_quad_p5p10_b2_b,   hdiv_leg_quad_p6p0_b2_b,   hdiv_leg_quad_p6p1_b2_b,   hdiv_leg_quad_p6p2_b2_b,   hdiv_leg_quad_p6p3_b2_b,   hdiv_leg_quad_p6p4_b2_b,   hdiv_leg_quad_p6p5_b2_b,   hdiv_leg_quad_p6p6_b2_b,   hdiv_leg_quad_p6p7_b2_b,   hdiv_leg_quad_p6p8_b2_b,   hdiv_leg_quad_p6p9_b2_b,   hdiv_leg_quad_p6p10_b2_b,   hdiv_leg_quad_p7p0_b2_b,   hdiv_leg_quad_p7p1_b2_b,   hdiv_leg_quad_p7p2_b2_b,   hdiv_leg_quad_p7p3_b2_b,   hdiv_leg_quad_p7p4_b2_b,   hdiv_leg_quad_p7p5_b2_b,   hdiv_leg_quad_p7p6_b2_b,   hdiv_leg_quad_p7p7_b2_b,   hdiv_leg_quad_p7p8_b2_b,   hdiv_leg_quad_p7p9_b2_b,   hdiv_leg_quad_p7p10_b2_b,   hdiv_leg_quad_p8p0_b2_b,   hdiv_leg_quad_p8p1_b2_b,   hdiv_leg_quad_p8p2_b2_b,   hdiv_leg_quad_p8p3_b2_b,   hdiv_leg_quad_p8p4_b2_b,   hdiv_leg_quad_p8p5_b2_b,   hdiv_leg_quad_p8p6_b2_b,   hdiv_leg_quad_p8p7_b2_b,   hdiv_leg_quad_p8p8_b2_b,   hdiv_leg_quad_p8p9_b2_b,   hdiv_leg_quad_p8p10_b2_b,   hdiv_leg_quad_p9p0_b2_b,   hdiv_leg_quad_p9p1_b2_b,   hdiv_leg_quad_p9p2_b2_b,   hdiv_leg_quad_p9p3_b2_b,   hdiv_leg_quad_p9p4_b2_b,   hdiv_leg_quad_p9p5_b2_b,   hdiv_leg_quad_p9p6_b2_b,   hdiv_leg_quad_p9p7_b2_b,   hdiv_leg_quad_p9p8_b2_b,   hdiv_leg_quad_p9p9_b2_b,   hdiv_leg_quad_p9p10_b2_b,   hdiv_leg_quad_p10p0_b2_b,   hdiv_leg_quad_p10p1_b2_b,   hdiv_leg_quad_p10p2_b2_b,   hdiv_leg_quad_p10p3_b2_b,   hdiv_leg_quad_p10p4_b2_b,   hdiv_leg_quad_p10p5_b2_b,   hdiv_leg_quad_p10p6_b2_b,   hdiv_leg_quad_p10p7_b2_b,   hdiv_leg_quad_p10p8_b2_b,   hdiv_leg_quad_p10p9_b2_b,   hdiv_leg_quad_p10p10_b2_b,   hdiv_leg_quad_p11p0_b2_b,   hdiv_leg_quad_p11p1_b2_b,   hdiv_leg_quad_p11p2_b2_b,   hdiv_leg_quad_p11p3_b2_b,   hdiv_leg_quad_p11p4_b2_b,   hdiv_leg_quad_p11p5_b2_b,   hdiv_leg_quad_p11p6_b2_b,   hdiv_leg_quad_p11p7_b2_b,   hdiv_leg_quad_p11p8_b2_b,   hdiv_leg_quad_p11p9_b2_b,   hdiv_leg_quad_p11p10_b2_b, };

    static Shapeset::shape_fn_t hdiv_leg_quad_fn_ax[] =
    {
      hdiv_leg_quad_p0_e1_ax_0, hdiv_leg_quad_p0_e1_ax_1, hdiv_leg_quad_p0_e2_ax_0, hdiv_leg_quad_p0_e2_ax_1,  hdiv_leg_quad_p0_e3_ax, hdiv_leg_quad_p0_e3_ax, hdiv_leg_quad_p0_e4_ax, hdiv_leg_quad_p0_e4_ax,
      hdiv_leg_quad_p1_e1_ax, hdiv_leg_quad_p1_e1_ax, hdiv_leg_quad_p1_e2_ax, hdiv_leg_quad_p1_e2_ax, hdiv_leg_quad_p1_e3_ax, hdiv_leg_quad_p1_e3_ax, hdiv_leg_quad_p1_e4_ax, hdiv_leg_quad_p1_e4_ax,
      hdiv_leg_quad_p2_e1_ax_0, hdiv_leg_quad_p2_e1_ax_1, hdiv_leg_quad_p2_e2_ax_0, hdiv_leg_quad_p2_e2_ax_1, hdiv_leg_quad_p2_e3_ax, hdiv_leg_quad_p2_e3_ax, hdiv_leg_quad_p2_e4_ax, hdiv_leg_quad_p2_e4_ax,
      hdiv_leg_quad_p3_e1_ax, hdiv_leg_quad_p3_e1_ax, hdiv_leg_quad_p3_e2_ax, hdiv_leg_quad_p3_e2_ax, hdiv_leg_quad_p3_e3_ax, hdiv_leg_quad_p3_e3_ax, hdiv_leg_quad_p3_e4_ax, hdiv_leg_quad_p3_e4_ax,
      hdiv_leg_quad_p4_e1_ax_0, hdiv_leg_quad_p4_e1_ax_1, hdiv_leg_quad_p4_e2_ax_0, hdiv_leg_quad_p4_e2_ax_1, hdiv_leg_quad_p4_e3_ax, hdiv_leg_quad_p4_e3_ax, hdiv_leg_quad_p4_e4_ax, hdiv_leg_quad_p4_e4_ax,
      hdiv_leg_quad_p5_e1_ax, hdiv_leg_quad_p5_e1_ax, hdiv_leg_quad_p5_e2_ax, hdiv_leg_quad_p5_e2_ax, hdiv_leg_quad_p5_e3_ax, hdiv_leg_quad_p5_e3_ax, hdiv_leg_quad_p5_e4_ax, hdiv_leg_quad_p5_e4_ax,
      hdiv_leg_quad_p6_e1_ax_0, hdiv_leg_quad_p6_e1_ax_1, hdiv_leg_quad_p6_e2_ax_0, hdiv_leg_quad_p6_e2_ax_1, hdiv_leg_quad_p6_e3_ax, hdiv_leg_quad_p6_e3_ax, hdiv_leg_quad_p6_e4_ax, hdiv_leg_quad_p6_e4_ax,
      hdiv_leg_quad_p7_e1_ax, hdiv_leg_quad_p7_e1_ax, hdiv_leg_quad_p7_e2_ax, hdiv_leg_quad_p7_e2_ax, hdiv_leg_quad_p7_e3_ax, hdiv_leg_quad_p7_e3_ax, hdiv_leg_quad_p7_e4_ax, hdiv_leg_quad_p7_e4_ax,
      hdiv_leg_quad_p8_e1_ax_0, hdiv_leg_quad_p8_e1_ax_1, hdiv_leg_quad_p8_e2_ax_0, hdiv_leg_quad_p8_e2_ax_1, hdiv_leg_quad_p8_e3_ax, hdiv_leg_quad_p8_e3_ax, hdiv_leg_quad_p8_e4_ax, hdiv_leg_quad_p8_e4_ax,
      hdiv_leg_quad_p9_e1_ax, hdiv_leg_quad_p9_e1_ax, hdiv_leg_quad_p9_e2_ax, hdiv_leg_quad_p9_e2_ax, hdiv_leg_quad_p9_e3_ax, hdiv_leg_quad_p9_e3_ax, hdiv_leg_quad_p9_e4_ax, hdiv_leg_quad_p9_e4_ax,
      hdiv_leg_quad_p10_e1_ax_0, hdiv_leg_quad_p10_e1_ax_1, hdiv_leg_quad_p10_e2_ax_0, hdiv_leg_quad_p10_e2_ax_1, hdiv_leg_quad_p10_e3_ax, hdiv_leg_quad_p10_e3_ax, hdiv_leg_quad_p10_e4_ax, hdiv_leg_quad_p10_e4_ax,

      hdiv_leg_quad_p0p2_b1_ax,   hdiv_leg_quad_p0p3_b1_ax,   hdiv_leg_quad_p0p4_b1_ax,   hdiv_leg_quad_p0p5_b1_ax,   hdiv_leg_quad_p0p6_b1_ax,   hdiv_leg_quad_p0p7_b1_ax,   hdiv_leg_quad_p0p8_b1_ax,   hdiv_leg_quad_p0p9_b1_ax,   hdiv_leg_quad_p0p10_b1_ax,   hdiv_leg_quad_p0p11_b1_ax,   hdiv_leg_quad_p1p2_b1_ax,   hdiv_leg_quad_p1p3_b1_ax,   hdiv_leg_quad_p1p4_b1_ax,   hdiv_leg_quad_p1p5_b1_ax,   hdiv_leg_quad_p1p6_b1_ax,   hdiv_leg_quad_p1p7_b1_ax,   hdiv_leg_quad_p1p8_b1_ax,   hdiv_leg_quad_p1p9_b1_ax,   hdiv_leg_quad_p1p10_b1_ax,   hdiv_leg_quad_p1p11_b1_ax,   hdiv_leg_quad_p2p2_b1_ax,   hdiv_leg_quad_p2p3_b1_ax,   hdiv_leg_quad_p2p4_b1_ax,   hdiv_leg_quad_p2p5_b1_ax,   hdiv_leg_quad_p2p6_b1_ax,   hdiv_leg_quad_p2p7_b1_ax,   hdiv_leg_quad_p2p8_b1_ax,   hdiv_leg_quad_p2p9_b1_ax,   hdiv_leg_quad_p2p10_b1_ax,   hdiv_leg_quad_p2p11_b1_ax,   hdiv_leg_quad_p3p2_b1_ax,   hdiv_leg_quad_p3p3_b1_ax,   hdiv_leg_quad_p3p4_b1_ax,   hdiv_leg_quad_p3p5_b1_ax,   hdiv_leg_quad_p3p6_b1_ax,   hdiv_leg_quad_p3p7_b1_ax,   hdiv_leg_quad_p3p8_b1_ax,   hdiv_leg_quad_p3p9_b1_ax,   hdiv_leg_quad_p3p10_b1_ax,   hdiv_leg_quad_p3p11_b1_ax,   hdiv_leg_quad_p4p2_b1_ax,   hdiv_leg_quad_p4p3_b1_ax,   hdiv_leg_quad_p4p4_b1_ax,   hdiv_leg_quad_p4p5_b1_ax,   hdiv_leg_quad_p4p6_b1_ax,   hdiv_leg_quad_p4p7_b1_ax,   hdiv_leg_quad_p4p8_b1_ax,   hdiv_leg_quad_p4p9_b1_ax,   hdiv_leg_quad_p4p10_b1_ax,   hdiv_leg_quad_p4p11_b1_ax,   hdiv_leg_quad_p5p2_b1_ax,   hdiv_leg_quad_p5p3_b1_ax,   hdiv_leg_quad_p5p4_b1_ax,   hdiv_leg_quad_p5p5_b1_ax,   hdiv_leg_quad_p5p6_b1_ax,   hdiv_leg_quad_p5p7_b1_ax,   hdiv_leg_quad_p5p8_b1_ax,   hdiv_leg_quad_p5p9_b1_ax,   hdiv_leg_quad_p5p10_b1_ax,   hdiv_leg_quad_p5p11_b1_ax,   hdiv_leg_quad_p6p2_b1_ax,   hdiv_leg_quad_p6p3_b1_ax,   hdiv_leg_quad_p6p4_b1_ax,   hdiv_leg_quad_p6p5_b1_ax,   hdiv_leg_quad_p6p6_b1_ax,   hdiv_leg_quad_p6p7_b1_ax,   hdiv_leg_quad_p6p8_b1_ax,   hdiv_leg_quad_p6p9_b1_ax,   hdiv_leg_quad_p6p10_b1_ax,   hdiv_leg_quad_p6p11_b1_ax,   hdiv_leg_quad_p7p2_b1_ax,   hdiv_leg_quad_p7p3_b1_ax,   hdiv_leg_quad_p7p4_b1_ax,   hdiv_leg_quad_p7p5_b1_ax,   hdiv_leg_quad_p7p6_b1_ax,   hdiv_leg_quad_p7p7_b1_ax,   hdiv_leg_quad_p7p8_b1_ax,   hdiv_leg_quad_p7p9_b1_ax,   hdiv_leg_quad_p7p10_b1_ax,   hdiv_leg_quad_p7p11_b1_ax,   hdiv_leg_quad_p8p2_b1_ax,   hdiv_leg_quad_p8p3_b1_ax,   hdiv_leg_quad_p8p4_b1_ax,   hdiv_leg_quad_p8p5_b1_ax,   hdiv_leg_quad_p8p6_b1_ax,   hdiv_leg_quad_p8p7_b1_ax,   hdiv_leg_quad_p8p8_b1_ax,   hdiv_leg_quad_p8p9_b1_ax,   hdiv_leg_quad_p8p10_b1_ax,   hdiv_leg_quad_p8p11_b1_ax,   hdiv_leg_quad_p9p2_b1_ax,   hdiv_leg_quad_p9p3_b1_ax,   hdiv_leg_quad_p9p4_b1_ax,   hdiv_leg_quad_p9p5_b1_ax,   hdiv_leg_quad_p9p6_b1_ax,   hdiv_leg_quad_p9p7_b1_ax,   hdiv_leg_quad_p9p8_b1_ax,   hdiv_leg_quad_p9p9_b1_ax,   hdiv_leg_quad_p9p10_b1_ax,   hdiv_leg_quad_p9p11_b1_ax,   hdiv_leg_quad_p10p2_b1_ax,   hdiv_leg_quad_p10p3_b1_ax,   hdiv_leg_quad_p10p4_b1_ax,   hdiv_leg_quad_p10p5_b1_ax,   hdiv_leg_quad_p10p6_b1_ax,   hdiv_leg_quad_p10p7_b1_ax,   hdiv_leg_quad_p10p8_b1_ax,   hdiv_leg_quad_p10p9_b1_ax,   hdiv_leg_quad_p10p10_b1_ax,   hdiv_leg_quad_p10p11_b1_ax,   hdiv_leg_quad_p2p0_b2_ax,   hdiv_leg_quad_p2p1_b2_ax,   hdiv_leg_quad_p2p2_b2_ax,   hdiv_leg_quad_p2p3_b2_ax,   hdiv_leg_quad_p2p4_b2_ax,   hdiv_leg_quad_p2p5_b2_ax,   hdiv_leg_quad_p2p6_b2_ax,   hdiv_leg_quad_p2p7_b2_ax,   hdiv_leg_quad_p2p8_b2_ax,   hdiv_leg_quad_p2p9_b2_ax,   hdiv_leg_quad_p2p10_b2_ax,   hdiv_leg_quad_p3p0_b2_ax,   hdiv_leg_quad_p3p1_b2_ax,   hdiv_leg_quad_p3p2_b2_ax,   hdiv_leg_quad_p3p3_b2_ax,   hdiv_leg_quad_p3p4_b2_ax,   hdiv_leg_quad_p3p5_b2_ax,   hdiv_leg_quad_p3p6_b2_ax,   hdiv_leg_quad_p3p7_b2_ax,   hdiv_leg_quad_p3p8_b2_ax,   hdiv_leg_quad_p3p9_b2_ax,   hdiv_leg_quad_p3p10_b2_ax,   hdiv_leg_quad_p4p0_b2_ax,   hdiv_leg_quad_p4p1_b2_ax,   hdiv_leg_quad_p4p2_b2_ax,   hdiv_leg_quad_p4p3_b2_ax,   hdiv_leg_quad_p4p4_b2_ax,   hdiv_leg_quad_p4p5_b2_ax,   hdiv_leg_quad_p4p6_b2_ax,   hdiv_leg_quad_p4p7_b2_ax,   hdiv_leg_quad_p4p8_b2_ax,   hdiv_leg_quad_p4p9_b2_ax,   hdiv_leg_quad_p4p10_b2_ax,   hdiv_leg_quad_p5p0_b2_ax,   hdiv_leg_quad_p5p1_b2_ax,   hdiv_leg_quad_p5p2_b2_ax,   hdiv_leg_quad_p5p3_b2_ax,   hdiv_leg_quad_p5p4_b2_ax,   hdiv_leg_quad_p5p5_b2_ax,   hdiv_leg_quad_p5p6_b2_ax,   hdiv_leg_quad_p5p7_b2_ax,   hdiv_leg_quad_p5p8_b2_ax,   hdiv_leg_quad_p5p9_b2_ax,   hdiv_leg_quad_p5p10_b2_ax,   hdiv_leg_quad_p6p0_b2_ax,   hdiv_leg_quad_p6p1_b2_ax,   hdiv_leg_quad_p6p2_b2_ax,   hdiv_leg_quad_p6p3_b2_ax,   hdiv_leg_quad_p6p4_b2_ax,   hdiv_leg_quad_p6p5_b2_ax,   hdiv_leg_quad_p6p6_b2_ax,   hdiv_leg_quad_p6p7_b2_ax,   hdiv_leg_quad_p6p8_b2_ax,   hdiv_leg_quad_p6p9_b2_ax,   hdiv_leg_quad_p6p10_b2_ax,   hdiv_leg_quad_p7p0_b2_ax,   hdiv_leg_quad_p7p1_b2_ax,   hdiv_leg_quad_p7p2_b2_ax,   hdiv_leg_quad_p7p3_b2_ax,   hdiv_leg_quad_p7p4_b2_ax,   hdiv_leg_quad_p7p5_b2_ax,   hdiv_leg_quad_p7p6_b2_ax,   hdiv_leg_quad_p7p7_b2_ax,   hdiv_leg_quad_p7p8_b2_ax,   hdiv_leg_quad_p7p9_b2_ax,   hdiv_leg_quad_p7p10_b2_ax,   hdiv_leg_quad_p8p0_b2_ax,   hdiv_leg_quad_p8p1_b2_ax,   hdiv_leg_quad_p8p2_b2_ax,   hdiv_leg_quad_p8p3_b2_ax,   hdiv_leg_quad_p8p4_b2_ax,   hdiv_leg_quad_p8p5_b2_ax,   hdiv_leg_quad_p8p6_b2_ax,   hdiv_leg_quad_p8p7_b2_ax,   hdiv_leg_quad_p8p8_b2_ax,   hdiv_leg_quad_p8p9_b2_ax,   hdiv_leg_quad_p8p10_b2_ax,   hdiv_leg_quad_p9p0_b2_ax,   hdiv_leg_quad_p9p1_b2_ax,   hdiv_leg_quad_p9p2_b2_ax,   hdiv_leg_quad_p9p3_b2_ax,   hdiv_leg_quad_p9p4_b2_ax,   hdiv_leg_quad_p9p5_b2_ax,   hdiv_leg_quad_p9p6_b2_ax,   hdiv_leg_quad_p9p7_b2_ax,   hdiv_leg_quad_p9p8_b2_ax,   hdiv_leg_quad_p9p9_b2_ax,   hdiv_leg_quad_p9p10_b2_ax,   hdiv_leg_quad_p10p0_b2_ax,   hdiv_leg_quad_p10p1_b2_ax,   hdiv_leg_quad_p10p2_b2_ax,   hdiv_leg_quad_p10p3_b2_ax,   hdiv_leg_quad_p10p4_b2_ax,   hdiv_leg_quad_p10p5_b2_ax,   hdiv_leg_quad_p10p6_b2_ax,   hdiv_leg_quad_p10p7_b2_ax,   hdiv_leg_quad_p10p8_b2_ax,   hdiv_leg_quad_p10p9_b2_ax,   hdiv_leg_quad_p10p10_b2_ax,   hdiv_leg_quad_p11p0_b2_ax,   hdiv_leg_quad_p11p1_b2_ax,   hdiv_leg_quad_p11p2_b2_ax,   hdiv_leg_quad_p11p3_b2_ax,   hdiv_leg_quad_p11p4_b2_ax,   hdiv_leg_quad_p11p5_b2_ax,   hdiv_leg_quad_p11p6_b2_ax,   hdiv_leg_quad_p11p7_b2_ax,   hdiv_leg_quad_p11p8_b2_ax,   hdiv_leg_quad_p11p9_b2_ax,   hdiv_leg_quad_p11p10_b2_ax, };

    static Shapeset::shape_fn_t hdiv_leg_quad_fn_bx[] =
    {
      hdiv_leg_quad_p0_e1_bx, hdiv_leg_quad_p0_e1_bx, hdiv_leg_quad_p0_e2_bx, hdiv_leg_quad_p0_e2_bx, hdiv_leg_quad_p0_e3_bx_0, hdiv_leg_quad_p0_e3_bx_1, hdiv_leg_quad_p0_e4_bx_0, hdiv_leg_quad_p0_e4_bx_1,
      hdiv_leg_quad_p1_e1_bx, hdiv_leg_quad_p1_e1_bx, hdiv_leg_quad_p1_e2_bx, hdiv_leg_quad_p1_e2_bx, hdiv_leg_quad_p1_e3_bx, hdiv_leg_quad_p1_e3_bx, hdiv_leg_quad_p1_e4_bx, hdiv_leg_quad_p1_e4_bx,
      hdiv_leg_quad_p2_e1_bx, hdiv_leg_quad_p2_e1_bx, hdiv_leg_quad_p2_e2_bx, hdiv_leg_quad_p2_e2_bx, hdiv_leg_quad_p2_e3_bx_0, hdiv_leg_quad_p2_e3_bx_1, hdiv_leg_quad_p2_e4_bx_0, hdiv_leg_quad_p2_e4_bx_1,
      hdiv_leg_quad_p3_e1_bx, hdiv_leg_quad_p3_e1_bx, hdiv_leg_quad_p3_e2_bx, hdiv_leg_quad_p3_e2_bx, hdiv_leg_quad_p3_e3_bx, hdiv_leg_quad_p3_e3_bx, hdiv_leg_quad_p3_e4_bx, hdiv_leg_quad_p3_e4_bx,
      hdiv_leg_quad_p4_e1_bx, hdiv_leg_quad_p4_e1_bx, hdiv_leg_quad_p4_e2_bx, hdiv_leg_quad_p4_e2_bx, hdiv_leg_quad_p4_e3_bx_0, hdiv_leg_quad_p4_e3_bx_1, hdiv_leg_quad_p4_e4_bx_0, hdiv_leg_quad_p4_e4_bx_1,
      hdiv_leg_quad_p5_e1_bx, hdiv_leg_quad_p5_e1_bx, hdiv_leg_quad_p5_e2_bx, hdiv_leg_quad_p5_e2_bx, hdiv_leg_quad_p5_e3_bx, hdiv_leg_quad_p5_e3_bx, hdiv_leg_quad_p5_e4_bx, hdiv_leg_quad_p5_e4_bx,
      hdiv_leg_quad_p6_e1_bx, hdiv_leg_quad_p6_e1_bx, hdiv_leg_quad_p6_e2_bx, hdiv_leg_quad_p6_e2_bx, hdiv_leg_quad_p6_e3_bx_0, hdiv_leg_quad_p6_e3_bx_1, hdiv_leg_quad_p6_e4_bx_0, hdiv_leg_quad_p6_e4_bx_1,
      hdiv_leg_quad_p7_e1_bx, hdiv_leg_quad_p7_e1_bx, hdiv_leg_quad_p7_e2_bx, hdiv_leg_quad_p7_e2_bx, hdiv_leg_quad_p7_e3_bx, hdiv_leg_quad_p7_e3_bx, hdiv_leg_quad_p7_e4_bx, hdiv_leg_quad_p7_e4_bx,
      hdiv_leg_quad_p8_e1_bx, hdiv_leg_quad_p8_e1_bx, hdiv_leg_quad_p8_e2_bx, hdiv_leg_quad_p8_e2_bx, hdiv_leg_quad_p8_e3_bx_0, hdiv_leg_quad_p8_e3_bx_1, hdiv_leg_quad_p8_e4_bx_0, hdiv_leg_quad_p8_e4_bx_1,
      hdiv_leg_quad_p9_e1_bx, hdiv_leg_quad_p9_e1_bx, hdiv_leg_quad_p9_e2_bx, hdiv_leg_quad_p9_e2_bx, hdiv_leg_quad_p9_e3_bx, hdiv_leg_quad_p9_e3_bx, hdiv_leg_quad_p9_e4_bx, hdiv_leg_quad_p9_e4_bx,
      hdiv_leg_quad_p10_e1_bx, hdiv_leg_quad_p10_e1_bx, hdiv_leg_quad_p10_e2_bx, hdiv_leg_quad_p10_e2_bx, hdiv_leg_quad_p10_e3_bx_0, hdiv_leg_quad_p10_e3_bx_1, hdiv_leg_quad_p10_e4_bx_0, hdiv_leg_quad_p10_e4_bx_1,

      hdiv_leg_quad_p0p2_b1_bx,   hdiv_leg_quad_p0p3_b1_bx,   hdiv_leg_quad_p0p4_b1_bx,   hdiv_leg_quad_p0p5_b1_bx,   hdiv_leg_quad_p0p6_b1_bx,   hdiv_leg_quad_p0p7_b1_bx,   hdiv_leg_quad_p0p8_b1_bx,   hdiv_leg_quad_p0p9_b1_bx,   hdiv_leg_quad_p0p10_b1_bx,   hdiv_leg_quad_p0p11_b1_bx,   hdiv_leg_quad_p1p2_b1_bx,   hdiv_leg_quad_p1p3_b1_bx,   hdiv_leg_quad_p1p4_b1_bx,   hdiv_leg_quad_p1p5_b1_bx,   hdiv_leg_quad_p1p6_b1_bx,   hdiv_leg_quad_p1p7_b1_bx,   hdiv_leg_quad_p1p8_b1_bx,   hdiv_leg_quad_p1p9_b1_bx,   hdiv_leg_quad_p1p10_b1_bx,   hdiv_leg_quad_p1p11_b1_bx,   hdiv_leg_quad_p2p2_b1_bx,   hdiv_leg_quad_p2p3_b1_bx,   hdiv_leg_quad_p2p4_b1_bx,   hdiv_leg_quad_p2p5_b1_bx,   hdiv_leg_quad_p2p6_b1_bx,   hdiv_leg_quad_p2p7_b1_bx,   hdiv_leg_quad_p2p8_b1_bx,   hdiv_leg_quad_p2p9_b1_bx,   hdiv_leg_quad_p2p10_b1_bx,   hdiv_leg_quad_p2p11_b1_bx,   hdiv_leg_quad_p3p2_b1_bx,   hdiv_leg_quad_p3p3_b1_bx,   hdiv_leg_quad_p3p4_b1_bx,   hdiv_leg_quad_p3p5_b1_bx,   hdiv_leg_quad_p3p6_b1_bx,   hdiv_leg_quad_p3p7_b1_bx,   hdiv_leg_quad_p3p8_b1_bx,   hdiv_leg_quad_p3p9_b1_bx,   hdiv_leg_quad_p3p10_b1_bx,   hdiv_leg_quad_p3p11_b1_bx,   hdiv_leg_quad_p4p2_b1_bx,   hdiv_leg_quad_p4p3_b1_bx,   hdiv_leg_quad_p4p4_b1_bx,   hdiv_leg_quad_p4p5_b1_bx,   hdiv_leg_quad_p4p6_b1_bx,   hdiv_leg_quad_p4p7_b1_bx,   hdiv_leg_quad_p4p8_b1_bx,   hdiv_leg_quad_p4p9_b1_bx,   hdiv_leg_quad_p4p10_b1_bx,   hdiv_leg_quad_p4p11_b1_bx,   hdiv_leg_quad_p5p2_b1_bx,   hdiv_leg_quad_p5p3_b1_bx,   hdiv_leg_quad_p5p4_b1_bx,   hdiv_leg_quad_p5p5_b1_bx,   hdiv_leg_quad_p5p6_b1_bx,   hdiv_leg_quad_p5p7_b1_bx,   hdiv_leg_quad_p5p8_b1_bx,   hdiv_leg_quad_p5p9_b1_bx,   hdiv_leg_quad_p5p10_b1_bx,   hdiv_leg_quad_p5p11_b1_bx,   hdiv_leg_quad_p6p2_b1_bx,   hdiv_leg_quad_p6p3_b1_bx,   hdiv_leg_quad_p6p4_b1_bx,   hdiv_leg_quad_p6p5_b1_bx,   hdiv_leg_quad_p6p6_b1_bx,   hdiv_leg_quad_p6p7_b1_bx,   hdiv_leg_quad_p6p8_b1_bx,   hdiv_leg_quad_p6p9_b1_bx,   hdiv_leg_quad_p6p10_b1_bx,   hdiv_leg_quad_p6p11_b1_bx,   hdiv_leg_quad_p7p2_b1_bx,   hdiv_leg_quad_p7p3_b1_bx,   hdiv_leg_quad_p7p4_b1_bx,   hdiv_leg_quad_p7p5_b1_bx,   hdiv_leg_quad_p7p6_b1_bx,   hdiv_leg_quad_p7p7_b1_bx,   hdiv_leg_quad_p7p8_b1_bx,   hdiv_leg_quad_p7p9_b1_bx,   hdiv_leg_quad_p7p10_b1_bx,   hdiv_leg_quad_p7p11_b1_bx,   hdiv_leg_quad_p8p2_b1_bx,   hdiv_leg_quad_p8p3_b1_bx,   hdiv_leg_quad_p8p4_b1_bx,   hdiv_leg_quad_p8p5_b1_bx,   hdiv_leg_quad_p8p6_b1_bx,   hdiv_leg_quad_p8p7_b1_bx,   hdiv_leg_quad_p8p8_b1_bx,   hdiv_leg_quad_p8p9_b1_bx,   hdiv_leg_quad_p8p10_b1_bx,   hdiv_leg_quad_p8p11_b1_bx,   hdiv_leg_quad_p9p2_b1_bx,   hdiv_leg_quad_p9p3_b1_bx,   hdiv_leg_quad_p9p4_b1_bx,   hdiv_leg_quad_p9p5_b1_bx,   hdiv_leg_quad_p9p6_b1_bx,   hdiv_leg_quad_p9p7_b1_bx,   hdiv_leg_quad_p9p8_b1_bx,   hdiv_leg_quad_p9p9_b1_bx,   hdiv_leg_quad_p9p10_b1_bx,   hdiv_leg_quad_p9p11_b1_bx,   hdiv_leg_quad_p10p2_b1_bx,   hdiv_leg_quad_p10p3_b1_bx,   hdiv_leg_quad_p10p4_b1_bx,   hdiv_leg_quad_p10p5_b1_bx,   hdiv_leg_quad_p10p6_b1_bx,   hdiv_leg_quad_p10p7_b1_bx,   hdiv_leg_quad_p10p8_b1_bx,   hdiv_leg_quad_p10p9_b1_bx,   hdiv_leg_quad_p10p10_b1_bx,   hdiv_leg_quad_p10p11_b1_bx,   hdiv_leg_quad_p2p0_b2_bx,   hdiv_leg_quad_p2p1_b2_bx,   hdiv_leg_quad_p2p2_b2_bx,   hdiv_leg_quad_p2p3_b2_bx,   hdiv_leg_quad_p2p4_b2_bx,   hdiv_leg_quad_p2p5_b2_bx,   hdiv_leg_quad_p2p6_b2_bx,   hdiv_leg_quad_p2p7_b2_bx,   hdiv_leg_quad_p2p8_b2_bx,   hdiv_leg_quad_p2p9_b2_bx,   hdiv_leg_quad_p2p10_b2_bx,   hdiv_leg_quad_p3p0_b2_bx,   hdiv_leg_quad_p3p1_b2_bx,   hdiv_leg_quad_p3p2_b2_bx,   hdiv_leg_quad_p3p3_b2_bx,   hdiv_leg_quad_p3p4_b2_bx,   hdiv_leg_quad_p3p5_b2_bx,   hdiv_leg_quad_p3p6_b2_bx,   hdiv_leg_quad_p3p7_b2_bx,   hdiv_leg_quad_p3p8_b2_bx,   hdiv_leg_quad_p3p9_b2_bx,   hdiv_leg_quad_p3p10_b2_bx,   hdiv_leg_quad_p4p0_b2_bx,   hdiv_leg_quad_p4p1_b2_bx,   hdiv_leg_quad_p4p2_b2_bx,   hdiv_leg_quad_p4p3_b2_bx,   hdiv_leg_quad_p4p4_b2_bx,   hdiv_leg_quad_p4p5_b2_bx,   hdiv_leg_quad_p4p6_b2_bx,   hdiv_leg_quad_p4p7_b2_bx,   hdiv_leg_quad_p4p8_b2_bx,   hdiv_leg_quad_p4p9_b2_bx,   hdiv_leg_quad_p4p10_b2_bx,   hdiv_leg_quad_p5p0_b2_bx,   hdiv_leg_quad_p5p1_b2_bx,   hdiv_leg_quad_p5p2_b2_bx,   hdiv_leg_quad_p5p3_b2_bx,   hdiv_leg_quad_p5p4_b2_bx,   hdiv_leg_quad_p5p5_b2_bx,   hdiv_leg_quad_p5p6_b2_bx,   hdiv_leg_quad_p5p7_b2_bx,   hdiv_leg_quad_p5p8_b2_bx,   hdiv_leg_quad_p5p9_b2_bx,   hdiv_leg_quad_p5p10_b2_bx,   hdiv_leg_quad_p6p0_b2_bx,   hdiv_leg_quad_p6p1_b2_bx,   hdiv_leg_quad_p6p2_b2_bx,   hdiv_leg_quad_p6p3_b2_bx,   hdiv_leg_quad_p6p4_b2_bx,   hdiv_leg_quad_p6p5_b2_bx,   hdiv_leg_quad_p6p6_b2_bx,   hdiv_leg_quad_p6p7_b2_bx,   hdiv_leg_quad_p6p8_b2_bx,   hdiv_leg_quad_p6p9_b2_bx,   hdiv_leg_quad_p6p10_b2_bx,   hdiv_leg_quad_p7p0_b2_bx,   hdiv_leg_quad_p7p1_b2_bx,   hdiv_leg_quad_p7p2_b2_bx,   hdiv_leg_quad_p7p3_b2_bx,   hdiv_leg_quad_p7p4_b2_bx,   hdiv_leg_quad_p7p5_b2_bx,   hdiv_leg_quad_p7p6_b2_bx,   hdiv_leg_quad_p7p7_b2_bx,   hdiv_leg_quad_p7p8_b2_bx,   hdiv_leg_quad_p7p9_b2_bx,   hdiv_leg_quad_p7p10_b2_bx,   hdiv_leg_quad_p8p0_b2_bx,   hdiv_leg_quad_p8p1_b2_bx,   hdiv_leg_quad_p8p2_b2_bx,   hdiv_leg_quad_p8p3_b2_bx,   hdiv_leg_quad_p8p4_b2_bx,   hdiv_leg_quad_p8p5_b2_bx,   hdiv_leg_quad_p8p6_b2_bx,   hdiv_leg_quad_p8p7_b2_bx,   hdiv_leg_quad_p8p8_b2_bx,   hdiv_leg_quad_p8p9_b2_bx,   hdiv_leg_quad_p8p10_b2_bx,   hdiv_leg_quad_p9p0_b2_bx,   hdiv_leg_quad_p9p1_b2_bx,   hdiv_leg_quad_p9p2_b2_bx,   hdiv_leg_quad_p9p3_b2_bx,   hdiv_leg_quad_p9p4_b2_bx,   hdiv_leg_quad_p9p5_b2_bx,   hdiv_leg_quad_p9p6_b2_bx,   hdiv_leg_quad_p9p7_b2_bx,   hdiv_leg_quad_p9p8_b2_bx,   hdiv_leg_quad_p9p9_b2_bx,   hdiv_leg_quad_p9p10_b2_bx,   hdiv_leg_quad_p10p0_b2_bx,   hdiv_leg_quad_p10p1_b2_bx,   hdiv_leg_quad_p10p2_b2_bx,   hdiv_leg_quad_p10p3_b2_bx,   hdiv_leg_quad_p10p4_b2_bx,   hdiv_leg_quad_p10p5_b2_bx,   hdiv_leg_quad_p10p6_b2_bx,   hdiv_leg_quad_p10p7_b2_bx,   hdiv_leg_quad_p10p8_b2_bx,   hdiv_leg_quad_p10p9_b2_bx,   hdiv_leg_quad_p10p10_b2_bx,   hdiv_leg_quad_p11p0_b2_bx,   hdiv_leg_quad_p11p1_b2_bx,   hdiv_leg_quad_p11p2_b2_bx,   hdiv_leg_quad_p11p3_b2_bx,   hdiv_leg_quad_p11p4_b2_bx,   hdiv_leg_quad_p11p5_b2_bx,   hdiv_leg_quad_p11p6_b2_bx,   hdiv_leg_quad_p11p7_b2_bx,   hdiv_leg_quad_p11p8_b2_bx,   hdiv_leg_quad_p11p9_b2_bx,   hdiv_leg_quad_p11p10_b2_bx, };

    static Shapeset::shape_fn_t hdiv_leg_quad_fn_ay[] =
    {
      hdiv_leg_quad_p0_e1_ay_0, hdiv_leg_quad_p0_e1_ay_1, hdiv_leg_quad_p0_e2_ay_0, hdiv_leg_quad_p0_e2_ay_1,  hdiv_leg_quad_p0_e3_ay, hdiv_leg_quad_p0_e3_ay, hdiv_leg_quad_p0_e4_ay, hdiv_leg_quad_p0_e4_ay,
      hdiv_leg_quad_p1_e1_ay, hdiv_leg_quad_p1_e1_ay, hdiv_leg_quad_p1_e2_ay, hdiv_leg_quad_p1_e2_ay, hdiv_leg_quad_p1_e3_ay, hdiv_leg_quad_p1_e3_ay, hdiv_leg_quad_p1_e4_ay, hdiv_leg_quad_p1_e4_ay,
      hdiv_leg_quad_p2_e1_ay_0, hdiv_leg_quad_p2_e1_ay_1, hdiv_leg_quad_p2_e2_ay_0, hdiv_leg_quad_p2_e2_ay_1, hdiv_leg_quad_p2_e3_ay, hdiv_leg_quad_p2_e3_ay, hdiv_leg_quad_p2_e4_ay, hdiv_leg_quad_p2_e4_ay,
      hdiv_leg_quad_p3_e1_ay, hdiv_leg_quad_p3_e1_ay, hdiv_leg_quad_p3_e2_ay, hdiv_leg_quad_p3_e2_ay, hdiv_leg_quad_p3_e3_ay, hdiv_leg_quad_p3_e3_ay, hdiv_leg_quad_p3_e4_ay, hdiv_leg_quad_p3_e4_ay,
      hdiv_leg_quad_p4_e1_ay_0, hdiv_leg_quad_p4_e1_ay_1, hdiv_leg_quad_p4_e2_ay_0, hdiv_leg_quad_p4_e2_ay_1, hdiv_leg_quad_p4_e3_ay, hdiv_leg_quad_p4_e3_ay, hdiv_leg_quad_p4_e4_ay, hdiv_leg_quad_p4_e4_ay,
      hdiv_leg_quad_p5_e1_ay, hdiv_leg_quad_p5_e1_ay, hdiv_leg_quad_p5_e2_ay, hdiv_leg_quad_p5_e2_ay, hdiv_leg_quad_p5_e3_ay, hdiv_leg_quad_p5_e3_ay, hdiv_leg_quad_p5_e4_ay, hdiv_leg_quad_p5_e4_ay,
      hdiv_leg_quad_p6_e1_ay_0, hdiv_leg_quad_p6_e1_ay_1, hdiv_leg_quad_p6_e2_ay_0, hdiv_leg_quad_p6_e2_ay_1, hdiv_leg_quad_p6_e3_ay, hdiv_leg_quad_p6_e3_ay, hdiv_leg_quad_p6_e4_ay, hdiv_leg_quad_p6_e4_ay,
      hdiv_leg_quad_p7_e1_ay, hdiv_leg_quad_p7_e1_ay, hdiv_leg_quad_p7_e2_ay, hdiv_leg_quad_p7_e2_ay, hdiv_leg_quad_p7_e3_ay, hdiv_leg_quad_p7_e3_ay, hdiv_leg_quad_p7_e4_ay, hdiv_leg_quad_p7_e4_ay,
      hdiv_leg_quad_p8_e1_ay_0, hdiv_leg_quad_p8_e1_ay_1, hdiv_leg_quad_p8_e2_ay_0, hdiv_leg_quad_p8_e2_ay_1, hdiv_leg_quad_p8_e3_ay, hdiv_leg_quad_p8_e3_ay, hdiv_leg_quad_p8_e4_ay, hdiv_leg_quad_p8_e4_ay,
      hdiv_leg_quad_p9_e1_ay, hdiv_leg_quad_p9_e1_ay, hdiv_leg_quad_p9_e2_ay, hdiv_leg_quad_p9_e2_ay, hdiv_leg_quad_p9_e3_ay, hdiv_leg_quad_p9_e3_ay, hdiv_leg_quad_p9_e4_ay, hdiv_leg_quad_p9_e4_ay,
      hdiv_leg_quad_p10_e1_ay_0, hdiv_leg_quad_p10_e1_ay_1, hdiv_leg_quad_p10_e2_ay_0, hdiv_leg_quad_p10_e2_ay_1, hdiv_leg_quad_p10_e3_ay, hdiv_leg_quad_p10_e3_ay, hdiv_leg_quad_p10_e4_ay, hdiv_leg_quad_p10_e4_ay,

      hdiv_leg_quad_p0p2_b1_ay,   hdiv_leg_quad_p0p3_b1_ay,   hdiv_leg_quad_p0p4_b1_ay,   hdiv_leg_quad_p0p5_b1_ay,   hdiv_leg_quad_p0p6_b1_ay,   hdiv_leg_quad_p0p7_b1_ay,   hdiv_leg_quad_p0p8_b1_ay,   hdiv_leg_quad_p0p9_b1_ay,   hdiv_leg_quad_p0p10_b1_ay,   hdiv_leg_quad_p0p11_b1_ay,   hdiv_leg_quad_p1p2_b1_ay,   hdiv_leg_quad_p1p3_b1_ay,   hdiv_leg_quad_p1p4_b1_ay,   hdiv_leg_quad_p1p5_b1_ay,   hdiv_leg_quad_p1p6_b1_ay,   hdiv_leg_quad_p1p7_b1_ay,   hdiv_leg_quad_p1p8_b1_ay,   hdiv_leg_quad_p1p9_b1_ay,   hdiv_leg_quad_p1p10_b1_ay,   hdiv_leg_quad_p1p11_b1_ay,   hdiv_leg_quad_p2p2_b1_ay,   hdiv_leg_quad_p2p3_b1_ay,   hdiv_leg_quad_p2p4_b1_ay,   hdiv_leg_quad_p2p5_b1_ay,   hdiv_leg_quad_p2p6_b1_ay,   hdiv_leg_quad_p2p7_b1_ay,   hdiv_leg_quad_p2p8_b1_ay,   hdiv_leg_quad_p2p9_b1_ay,   hdiv_leg_quad_p2p10_b1_ay,   hdiv_leg_quad_p2p11_b1_ay,   hdiv_leg_quad_p3p2_b1_ay,   hdiv_leg_quad_p3p3_b1_ay,   hdiv_leg_quad_p3p4_b1_ay,   hdiv_leg_quad_p3p5_b1_ay,   hdiv_leg_quad_p3p6_b1_ay,   hdiv_leg_quad_p3p7_b1_ay,   hdiv_leg_quad_p3p8_b1_ay,   hdiv_leg_quad_p3p9_b1_ay,   hdiv_leg_quad_p3p10_b1_ay,   hdiv_leg_quad_p3p11_b1_ay,   hdiv_leg_quad_p4p2_b1_ay,   hdiv_leg_quad_p4p3_b1_ay,   hdiv_leg_quad_p4p4_b1_ay,   hdiv_leg_quad_p4p5_b1_ay,   hdiv_leg_quad_p4p6_b1_ay,   hdiv_leg_quad_p4p7_b1_ay,   hdiv_leg_quad_p4p8_b1_ay,   hdiv_leg_quad_p4p9_b1_ay,   hdiv_leg_quad_p4p10_b1_ay,   hdiv_leg_quad_p4p11_b1_ay,   hdiv_leg_quad_p5p2_b1_ay,   hdiv_leg_quad_p5p3_b1_ay,   hdiv_leg_quad_p5p4_b1_ay,   hdiv_leg_quad_p5p5_b1_ay,   hdiv_leg_quad_p5p6_b1_ay,   hdiv_leg_quad_p5p7_b1_ay,   hdiv_leg_quad_p5p8_b1_ay,   hdiv_leg_quad_p5p9_b1_ay,   hdiv_leg_quad_p5p10_b1_ay,   hdiv_leg_quad_p5p11_b1_ay,   hdiv_leg_quad_p6p2_b1_ay,   hdiv_leg_quad_p6p3_b1_ay,   hdiv_leg_quad_p6p4_b1_ay,   hdiv_leg_quad_p6p5_b1_ay,   hdiv_leg_quad_p6p6_b1_ay,   hdiv_leg_quad_p6p7_b1_ay,   hdiv_leg_quad_p6p8_b1_ay,   hdiv_leg_quad_p6p9_b1_ay,   hdiv_leg_quad_p6p10_b1_ay,   hdiv_leg_quad_p6p11_b1_ay,   hdiv_leg_quad_p7p2_b1_ay,   hdiv_leg_quad_p7p3_b1_ay,   hdiv_leg_quad_p7p4_b1_ay,   hdiv_leg_quad_p7p5_b1_ay,   hdiv_leg_quad_p7p6_b1_ay,   hdiv_leg_quad_p7p7_b1_ay,   hdiv_leg_quad_p7p8_b1_ay,   hdiv_leg_quad_p7p9_b1_ay,   hdiv_leg_quad_p7p10_b1_ay,   hdiv_leg_quad_p7p11_b1_ay,   hdiv_leg_quad_p8p2_b1_ay,   hdiv_leg_quad_p8p3_b1_ay,   hdiv_leg_quad_p8p4_b1_ay,   hdiv_leg_quad_p8p5_b1_ay,   hdiv_leg_quad_p8p6_b1_ay,   hdiv_leg_quad_p8p7_b1_ay,   hdiv_leg_quad_p8p8_b1_ay,   hdiv_leg_quad_p8p9_b1_ay,   hdiv_leg_quad_p8p10_b1_ay,   hdiv_leg_quad_p8p11_b1_ay,   hdiv_leg_quad_p9p2_b1_ay,   hdiv_leg_quad_p9p3_b1_ay,   hdiv_leg_quad_p9p4_b1_ay,   hdiv_leg_quad_p9p5_b1_ay,   hdiv_leg_quad_p9p6_b1_ay,   hdiv_leg_quad_p9p7_b1_ay,   hdiv_leg_quad_p9p8_b1_ay,   hdiv_leg_quad_p9p9_b1_ay,   hdiv_leg_quad_p9p10_b1_ay,   hdiv_leg_quad_p9p11_b1_ay,   hdiv_leg_quad_p10p2_b1_ay,   hdiv_leg_quad_p10p3_b1_ay,   hdiv_leg_quad_p10p4_b1_ay,   hdiv_leg_quad_p10p5_b1_ay,   hdiv_leg_quad_p10p6_b1_ay,   hdiv_leg_quad_p10p7_b1_ay,   hdiv_leg_quad_p10p8_b1_ay,   hdiv_leg_quad_p10p9_b1_ay,   hdiv_leg_quad_p10p10_b1_ay,   hdiv_leg_quad_p10p11_b1_ay,   hdiv_leg_quad_p2p0_b2_ay,   hdiv_leg_quad_p2p1_b2_ay,   hdiv_leg_quad_p2p2_b2_ay,   hdiv_leg_quad_p2p3_b2_ay,   hdiv_leg_quad_p2p4_b2_ay,   hdiv_leg_quad_p2p5_b2_ay,   hdiv_leg_quad_p2p6_b2_ay,   hdiv_leg_quad_p2p7_b2_ay,   hdiv_leg_quad_p2p8_b2_ay,   hdiv_leg_quad_p2p9_b2_ay,   hdiv_leg_quad_p2p10_b2_ay,   hdiv_leg_quad_p3p0_b2_ay,   hdiv_leg_quad_p3p1_b2_ay,   hdiv_leg_quad_p3p2_b2_ay,   hdiv_leg_quad_p3p3_b2_ay,   hdiv_leg_quad_p3p4_b2_ay,   hdiv_leg_quad_p3p5_b2_ay,   hdiv_leg_quad_p3p6_b2_ay,   hdiv_leg_quad_p3p7_b2_ay,   hdiv_leg_quad_p3p8_b2_ay,   hdiv_leg_quad_p3p9_b2_ay,   hdiv_leg_quad_p3p10_b2_ay,   hdiv_leg_quad_p4p0_b2_ay,   hdiv_leg_quad_p4p1_b2_ay,   hdiv_leg_quad_p4p2_b2_ay,   hdiv_leg_quad_p4p3_b2_ay,   hdiv_leg_quad_p4p4_b2_ay,   hdiv_leg_quad_p4p5_b2_ay,   hdiv_leg_quad_p4p6_b2_ay,   hdiv_leg_quad_p4p7_b2_ay,   hdiv_leg_quad_p4p8_b2_ay,   hdiv_leg_quad_p4p9_b2_ay,   hdiv_leg_quad_p4p10_b2_ay,   hdiv_leg_quad_p5p0_b2_ay,   hdiv_leg_quad_p5p1_b2_ay,   hdiv_leg_quad_p5p2_b2_ay,   hdiv_leg_quad_p5p3_b2_ay,   hdiv_leg_quad_p5p4_b2_ay,   hdiv_leg_quad_p5p5_b2_ay,   hdiv_leg_quad_p5p6_b2_ay,   hdiv_leg_quad_p5p7_b2_ay,   hdiv_leg_quad_p5p8_b2_ay,   hdiv_leg_quad_p5p9_b2_ay,   hdiv_leg_quad_p5p10_b2_ay,   hdiv_leg_quad_p6p0_b2_ay,   hdiv_leg_quad_p6p1_b2_ay,   hdiv_leg_quad_p6p2_b2_ay,   hdiv_leg_quad_p6p3_b2_ay,   hdiv_leg_quad_p6p4_b2_ay,   hdiv_leg_quad_p6p5_b2_ay,   hdiv_leg_quad_p6p6_b2_ay,   hdiv_leg_quad_p6p7_b2_ay,   hdiv_leg_quad_p6p8_b2_ay,   hdiv_leg_quad_p6p9_b2_ay,   hdiv_leg_quad_p6p10_b2_ay,   hdiv_leg_quad_p7p0_b2_ay,   hdiv_leg_quad_p7p1_b2_ay,   hdiv_leg_quad_p7p2_b2_ay,   hdiv_leg_quad_p7p3_b2_ay,   hdiv_leg_quad_p7p4_b2_ay,   hdiv_leg_quad_p7p5_b2_ay,   hdiv_leg_quad_p7p6_b2_ay,   hdiv_leg_quad_p7p7_b2_ay,   hdiv_leg_quad_p7p8_b2_ay,   hdiv_leg_quad_p7p9_b2_ay,   hdiv_leg_quad_p7p10_b2_ay,   hdiv_leg_quad_p8p0_b2_ay,   hdiv_leg_quad_p8p1_b2_ay,   hdiv_leg_quad_p8p2_b2_ay,   hdiv_leg_quad_p8p3_b2_ay,   hdiv_leg_quad_p8p4_b2_ay,   hdiv_leg_quad_p8p5_b2_ay,   hdiv_leg_quad_p8p6_b2_ay,   hdiv_leg_quad_p8p7_b2_ay,   hdiv_leg_quad_p8p8_b2_ay,   hdiv_leg_quad_p8p9_b2_ay,   hdiv_leg_quad_p8p10_b2_ay,   hdiv_leg_quad_p9p0_b2_ay,   hdiv_leg_quad_p9p1_b2_ay,   hdiv_leg_quad_p9p2_b2_ay,   hdiv_leg_quad_p9p3_b2_ay,   hdiv_leg_quad_p9p4_b2_ay,   hdiv_leg_quad_p9p5_b2_ay,   hdiv_leg_quad_p9p6_b2_ay,   hdiv_leg_quad_p9p7_b2_ay,   hdiv_leg_quad_p9p8_b2_ay,   hdiv_leg_quad_p9p9_b2_ay,   hdiv_leg_quad_p9p10_b2_ay,   hdiv_leg_quad_p10p0_b2_ay,   hdiv_leg_quad_p10p1_b2_ay,   hdiv_leg_quad_p10p2_b2_ay,   hdiv_leg_quad_p10p3_b2_ay,   hdiv_leg_quad_p10p4_b2_ay,   hdiv_leg_quad_p10p5_b2_ay,   hdiv_leg_quad_p10p6_b2_ay,   hdiv_leg_quad_p10p7_b2_ay,   hdiv_leg_quad_p10p8_b2_ay,   hdiv_leg_quad_p10p9_b2_ay,   hdiv_leg_quad_p10p10_b2_ay,   hdiv_leg_quad_p11p0_b2_ay,   hdiv_leg_quad_p11p1_b2_ay,   hdiv_leg_quad_p11p2_b2_ay,   hdiv_leg_quad_p11p3_b2_ay,   hdiv_leg_quad_p11p4_b2_ay,   hdiv_leg_quad_p11p5_b2_ay,   hdiv_leg_quad_p11p6_b2_ay,   hdiv_leg_quad_p11p7_b2_ay,   hdiv_leg_quad_p11p8_b2_ay,   hdiv_leg_quad_p11p9_b2_ay,   hdiv_leg_quad_p11p10_b2_ay, };

    static Shapeset::shape_fn_t hdiv_leg_quad_fn_by[] =
    {
      hdiv_leg_quad_p0_e1_by, hdiv_leg_quad_p0_e1_by, hdiv_leg_quad_p0_e2_by, hdiv_leg_quad_p0_e2_by, hdiv_leg_quad_p0_e3_by_0, hdiv_leg_quad_p0_e3_by_1, hdiv_leg_quad_p0_e4_by_0, hdiv_leg_quad_p0_e4_by_1,
      hdiv_leg_quad_p1_e1_by, hdiv_leg_quad_p1_e1_by, hdiv_leg_quad_p1_e2_by, hdiv_leg_quad_p1_e2_by, hdiv_leg_quad_p1_e3_by, hdiv_leg_quad_p1_e3_by, hdiv_leg_quad_p1_e4_by, hdiv_leg_quad_p1_e4_by,
      hdiv_leg_quad_p2_e1_by, hdiv_leg_quad_p2_e1_by, hdiv_leg_quad_p2_e2_by, hdiv_leg_quad_p2_e2_by, hdiv_leg_quad_p2_e3_by_0, hdiv_leg_quad_p2_e3_by_1, hdiv_leg_quad_p2_e4_by_0, hdiv_leg_quad_p2_e4_by_1,
      hdiv_leg_quad_p3_e1_by, hdiv_leg_quad_p3_e1_by, hdiv_leg_quad_p3_e2_by, hdiv_leg_quad_p3_e2_by, hdiv_leg_quad_p3_e3_by, hdiv_leg_quad_p3_e3_by, hdiv_leg_quad_p3_e4_by, hdiv_leg_quad_p3_e4_by,
      hdiv_leg_quad_p4_e1_by, hdiv_leg_quad_p4_e1_by, hdiv_leg_quad_p4_e2_by, hdiv_leg_quad_p4_e2_by, hdiv_leg_quad_p4_e3_by_0, hdiv_leg_quad_p4_e3_by_1, hdiv_leg_quad_p4_e4_by_0, hdiv_leg_quad_p4_e4_by_1,
      hdiv_leg_quad_p5_e1_by, hdiv_leg_quad_p5_e1_by, hdiv_leg_quad_p5_e2_by, hdiv_leg_quad_p5_e2_by, hdiv_leg_quad_p5_e3_by, hdiv_leg_quad_p5_e3_by, hdiv_leg_quad_p5_e4_by, hdiv_leg_quad_p5_e4_by,
      hdiv_leg_quad_p6_e1_by, hdiv_leg_quad_p6_e1_by, hdiv_leg_quad_p6_e2_by, hdiv_leg_quad_p6_e2_by, hdiv_leg_quad_p6_e3_by_0, hdiv_leg_quad_p6_e3_by_1, hdiv_leg_quad_p6_e4_by_0, hdiv_leg_quad_p6_e4_by_1,
      hdiv_leg_quad_p7_e1_by, hdiv_leg_quad_p7_e1_by, hdiv_leg_quad_p7_e2_by, hdiv_leg_quad_p7_e2_by, hdiv_leg_quad_p7_e3_by, hdiv_leg_quad_p7_e3_by, hdiv_leg_quad_p7_e4_by, hdiv_leg_quad_p7_e4_by,
      hdiv_leg_quad_p8_e1_by, hdiv_leg_quad_p8_e1_by, hdiv_leg_quad_p8_e2_by, hdiv_leg_quad_p8_e2_by, hdiv_leg_quad_p8_e3_by_0, hdiv_leg_quad_p8_e3_by_1, hdiv_leg_quad_p8_e4_by_0, hdiv_leg_quad_p8_e4_by_1,
      hdiv_leg_quad_p9_e1_by, hdiv_leg_quad_p9_e1_by, hdiv_leg_quad_p9_e2_by, hdiv_leg_quad_p9_e2_by, hdiv_leg_quad_p9_e3_by, hdiv_leg_quad_p9_e3_by, hdiv_leg_quad_p9_e4_by, hdiv_leg_quad_p9_e4_by,
      hdiv_leg_quad_p10_e1_by, hdiv_leg_quad_p10_e1_by, hdiv_leg_quad_p10_e2_by, hdiv_leg_quad_p10_e2_by, hdiv_leg_quad_p10_e3_by_0, hdiv_leg_quad_p10_e3_by_1, hdiv_leg_quad_p10_e4_by_0, hdiv_leg_quad_p10_e4_by_1,

      hdiv_leg_quad_p0p2_b1_by,   hdiv_leg_quad_p0p3_b1_by,   hdiv_leg_quad_p0p4_b1_by,   hdiv_leg_quad_p0p5_b1_by,   hdiv_leg_quad_p0p6_b1_by,   hdiv_leg_quad_p0p7_b1_by,   hdiv_leg_quad_p0p8_b1_by,   hdiv_leg_quad_p0p9_b1_by,   hdiv_leg_quad_p0p10_b1_by,   hdiv_leg_quad_p0p11_b1_by,   hdiv_leg_quad_p1p2_b1_by,   hdiv_leg_quad_p1p3_b1_by,   hdiv_leg_quad_p1p4_b1_by,   hdiv_leg_quad_p1p5_b1_by,   hdiv_leg_quad_p1p6_b1_by,   hdiv_leg_quad_p1p7_b1_by,   hdiv_leg_quad_p1p8_b1_by,   hdiv_leg_quad_p1p9_b1_by,   hdiv_leg_quad_p1p10_b1_by,   hdiv_leg_quad_p1p11_b1_by,   hdiv_leg_quad_p2p2_b1_by,   hdiv_leg_quad_p2p3_b1_by,   hdiv_leg_quad_p2p4_b1_by,   hdiv_leg_quad_p2p5_b1_by,   hdiv_leg_quad_p2p6_b1_by,   hdiv_leg_quad_p2p7_b1_by,   hdiv_leg_quad_p2p8_b1_by,   hdiv_leg_quad_p2p9_b1_by,   hdiv_leg_quad_p2p10_b1_by,   hdiv_leg_quad_p2p11_b1_by,   hdiv_leg_quad_p3p2_b1_by,   hdiv_leg_quad_p3p3_b1_by,   hdiv_leg_quad_p3p4_b1_by,   hdiv_leg_quad_p3p5_b1_by,   hdiv_leg_quad_p3p6_b1_by,   hdiv_leg_quad_p3p7_b1_by,   hdiv_leg_quad_p3p8_b1_by,   hdiv_leg_quad_p3p9_b1_by,   hdiv_leg_quad_p3p10_b1_by,   hdiv_leg_quad_p3p11_b1_by,   hdiv_leg_quad_p4p2_b1_by,   hdiv_leg_quad_p4p3_b1_by,   hdiv_leg_quad_p4p4_b1_by,   hdiv_leg_quad_p4p5_b1_by,   hdiv_leg_quad_p4p6_b1_by,   hdiv_leg_quad_p4p7_b1_by,   hdiv_leg_quad_p4p8_b1_by,   hdiv_leg_quad_p4p9_b1_by,   hdiv_leg_quad_p4p10_b1_by,   hdiv_leg_quad_p4p11_b1_by,   hdiv_leg_quad_p5p2_b1_by,   hdiv_leg_quad_p5p3_b1_by,   hdiv_leg_quad_p5p4_b1_by,   hdiv_leg_quad_p5p5_b1_by,   hdiv_leg_quad_p5p6_b1_by,   hdiv_leg_quad_p5p7_b1_by,   hdiv_leg_quad_p5p8_b1_by,   hdiv_leg_quad_p5p9_b1_by,   hdiv_leg_quad_p5p10_b1_by,   hdiv_leg_quad_p5p11_b1_by,   hdiv_leg_quad_p6p2_b1_by,   hdiv_leg_quad_p6p3_b1_by,   hdiv_leg_quad_p6p4_b1_by,   hdiv_leg_quad_p6p5_b1_by,   hdiv_leg_quad_p6p6_b1_by,   hdiv_leg_quad_p6p7_b1_by,   hdiv_leg_quad_p6p8_b1_by,   hdiv_leg_quad_p6p9_b1_by,   hdiv_leg_quad_p6p10_b1_by,   hdiv_leg_quad_p6p11_b1_by,   hdiv_leg_quad_p7p2_b1_by,   hdiv_leg_quad_p7p3_b1_by,   hdiv_leg_quad_p7p4_b1_by,   hdiv_leg_quad_p7p5_b1_by,   hdiv_leg_quad_p7p6_b1_by,   hdiv_leg_quad_p7p7_b1_by,   hdiv_leg_quad_p7p8_b1_by,   hdiv_leg_quad_p7p9_b1_by,   hdiv_leg_quad_p7p10_b1_by,   hdiv_leg_quad_p7p11_b1_by,   hdiv_leg_quad_p8p2_b1_by,   hdiv_leg_quad_p8p3_b1_by,   hdiv_leg_quad_p8p4_b1_by,   hdiv_leg_quad_p8p5_b1_by,   hdiv_leg_quad_p8p6_b1_by,   hdiv_leg_quad_p8p7_b1_by,   hdiv_leg_quad_p8p8_b1_by,   hdiv_leg_quad_p8p9_b1_by,   hdiv_leg_quad_p8p10_b1_by,   hdiv_leg_quad_p8p11_b1_by,   hdiv_leg_quad_p9p2_b1_by,   hdiv_leg_quad_p9p3_b1_by,   hdiv_leg_quad_p9p4_b1_by,   hdiv_leg_quad_p9p5_b1_by,   hdiv_leg_quad_p9p6_b1_by,   hdiv_leg_quad_p9p7_b1_by,   hdiv_leg_quad_p9p8_b1_by,   hdiv_leg_quad_p9p9_b1_by,   hdiv_leg_quad_p9p10_b1_by,   hdiv_leg_quad_p9p11_b1_by,   hdiv_leg_quad_p10p2_b1_by,   hdiv_leg_quad_p10p3_b1_by,   hdiv_leg_quad_p10p4_b1_by,   hdiv_leg_quad_p10p5_b1_by,   hdiv_leg_quad_p10p6_b1_by,   hdiv_leg_quad_p10p7_b1_by,   hdiv_leg_quad_p10p8_b1_by,   hdiv_leg_quad_p10p9_b1_by,   hdiv_leg_quad_p10p10_b1_by,   hdiv_leg_quad_p10p11_b1_by,   hdiv_leg_quad_p2p0_b2_by,   hdiv_leg_quad_p2p1_b2_by,   hdiv_leg_quad_p2p2_b2_by,   hdiv_leg_quad_p2p3_b2_by,   hdiv_leg_quad_p2p4_b2_by,   hdiv_leg_quad_p2p5_b2_by,   hdiv_leg_quad_p2p6_b2_by,   hdiv_leg_quad_p2p7_b2_by,   hdiv_leg_quad_p2p8_b2_by,   hdiv_leg_quad_p2p9_b2_by,   hdiv_leg_quad_p2p10_b2_by,   hdiv_leg_quad_p3p0_b2_by,   hdiv_leg_quad_p3p1_b2_by,   hdiv_leg_quad_p3p2_b2_by,   hdiv_leg_quad_p3p3_b2_by,   hdiv_leg_quad_p3p4_b2_by,   hdiv_leg_quad_p3p5_b2_by,   hdiv_leg_quad_p3p6_b2_by,   hdiv_leg_quad_p3p7_b2_by,   hdiv_leg_quad_p3p8_b2_by,   hdiv_leg_quad_p3p9_b2_by,   hdiv_leg_quad_p3p10_b2_by,   hdiv_leg_quad_p4p0_b2_by,   hdiv_leg_quad_p4p1_b2_by,   hdiv_leg_quad_p4p2_b2_by,   hdiv_leg_quad_p4p3_b2_by,   hdiv_leg_quad_p4p4_b2_by,   hdiv_leg_quad_p4p5_b2_by,   hdiv_leg_quad_p4p6_b2_by,   hdiv_leg_quad_p4p7_b2_by,   hdiv_leg_quad_p4p8_b2_by,   hdiv_leg_quad_p4p9_b2_by,   hdiv_leg_quad_p4p10_b2_by,   hdiv_leg_quad_p5p0_b2_by,   hdiv_leg_quad_p5p1_b2_by,   hdiv_leg_quad_p5p2_b2_by,   hdiv_leg_quad_p5p3_b2_by,   hdiv_leg_quad_p5p4_b2_by,   hdiv_leg_quad_p5p5_b2_by,   hdiv_leg_quad_p5p6_b2_by,   hdiv_leg_quad_p5p7_b2_by,   hdiv_leg_quad_p5p8_b2_by,   hdiv_leg_quad_p5p9_b2_by,   hdiv_leg_quad_p5p10_b2_by,   hdiv_leg_quad_p6p0_b2_by,   hdiv_leg_quad_p6p1_b2_by,   hdiv_leg_quad_p6p2_b2_by,   hdiv_leg_quad_p6p3_b2_by,   hdiv_leg_quad_p6p4_b2_by,   hdiv_leg_quad_p6p5_b2_by,   hdiv_leg_quad_p6p6_b2_by,   hdiv_leg_quad_p6p7_b2_by,   hdiv_leg_quad_p6p8_b2_by,   hdiv_leg_quad_p6p9_b2_by,   hdiv_leg_quad_p6p10_b2_by,   hdiv_leg_quad_p7p0_b2_by,   hdiv_leg_quad_p7p1_b2_by,   hdiv_leg_quad_p7p2_b2_by,   hdiv_leg_quad_p7p3_b2_by,   hdiv_leg_quad_p7p4_b2_by,   hdiv_leg_quad_p7p5_b2_by,   hdiv_leg_quad_p7p6_b2_by,   hdiv_leg_quad_p7p7_b2_by,   hdiv_leg_quad_p7p8_b2_by,   hdiv_leg_quad_p7p9_b2_by,   hdiv_leg_quad_p7p10_b2_by,   hdiv_leg_quad_p8p0_b2_by,   hdiv_leg_quad_p8p1_b2_by,   hdiv_leg_quad_p8p2_b2_by,   hdiv_leg_quad_p8p3_b2_by,   hdiv_leg_quad_p8p4_b2_by,   hdiv_leg_quad_p8p5_b2_by,   hdiv_leg_quad_p8p6_b2_by,   hdiv_leg_quad_p8p7_b2_by,   hdiv_leg_quad_p8p8_b2_by,   hdiv_leg_quad_p8p9_b2_by,   hdiv_leg_quad_p8p10_b2_by,   hdiv_leg_quad_p9p0_b2_by,   hdiv_leg_quad_p9p1_b2_by,   hdiv_leg_quad_p9p2_b2_by,   hdiv_leg_quad_p9p3_b2_by,   hdiv_leg_quad_p9p4_b2_by,   hdiv_leg_quad_p9p5_b2_by,   hdiv_leg_quad_p9p6_b2_by,   hdiv_leg_quad_p9p7_b2_by,   hdiv_leg_quad_p9p8_b2_by,   hdiv_leg_quad_p9p9_b2_by,   hdiv_leg_quad_p9p10_b2_by,   hdiv_leg_quad_p10p0_b2_by,   hdiv_leg_quad_p10p1_b2_by,   hdiv_leg_quad_p10p2_b2_by,   hdiv_leg_quad_p10p3_b2_by,   hdiv_leg_quad_p10p4_b2_by,   hdiv_leg_quad_p10p5_b2_by,   hdiv_leg_quad_p10p6_b2_by,   hdiv_leg_quad_p10p7_b2_by,   hdiv_leg_quad_p10p8_b2_by,   hdiv_leg_quad_p10p9_b2_by,   hdiv_leg_quad_p10p10_b2_by,   hdiv_leg_quad_p11p0_b2_by,   hdiv_leg_quad_p11p1_b2_by,   hdiv_leg_quad_p11p2_b2_by,   hdiv_leg_quad_p11p3_b2_by,   hdiv_leg_quad_p11p4_b2_by,   hdiv_leg_quad_p11p5_b2_by,   hdiv_leg_quad_p11p6_b2_by,   hdiv_leg_quad_p11p7_b2_by,   hdiv_leg_quad_p11p8_b2_by,   hdiv_leg_quad_p11p9_b2_by,   hdiv_leg_quad_p11p10_b2_by, };

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

    static int* hdiv_leg_quad_bubble_indices[] =
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

    static int hdiv_leg_quad_bubble_count[] =
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

    static int hdiv_leg_quad_vertex_indices[4] ={-1, -1, -1, -1};

    static int hdiv_leg_quad_edge_indices_0[] = { 4, 5, 12, 13, 20, 21, 28, 29, 36, 37, 44, 45, 52, 53, 60, 61, 68, 69, 76, 77, 84, 85, };

    static int hdiv_leg_quad_edge_indices_1[] = { 2, 3, 10, 11, 18, 19, 26, 27, 34, 35, 42, 43, 50, 51, 58, 59, 66, 67, 74, 75, 82, 83, };

    static int hdiv_leg_quad_edge_indices_2[] = { 6, 7, 14, 15, 22, 23, 30, 31, 38, 39, 46, 47, 54, 55, 62, 63, 70, 71, 78, 79, 86, 87, };

    static int hdiv_leg_quad_edge_indices_3[] = { 0, 1, 8, 9, 16, 17, 24, 25, 32, 33, 40, 41, 48, 49, 56, 57, 64, 65, 72, 73, 80, 81, };

    static int* hdiv_leg_quad_edge_indices[4] =
    {
      hdiv_leg_quad_edge_indices_0,
      hdiv_leg_quad_edge_indices_1,
      hdiv_leg_quad_edge_indices_2,
      hdiv_leg_quad_edge_indices_3,
    };

    #define oo H2D_MAKE_QUAD_ORDER

    static int hdiv_leg_quad_index_to_order[] =
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

    static Shapeset::shape_fn_t* hdiv_leg_quad_shape_fn_table[2] =
    {
      hdiv_leg_quad_fn_a,
      hdiv_leg_quad_fn_b
    };

    static Shapeset::shape_fn_t* hdiv_leg_quad_shape_fn_table_x[2] =
    {
      hdiv_leg_quad_fn_ax,
      hdiv_leg_quad_fn_bx
    };

    static Shapeset::shape_fn_t* hdiv_leg_quad_shape_fn_table_y[2] =
    {
      hdiv_leg_quad_fn_ay,
      hdiv_leg_quad_fn_by
    };

    //// tables and class constructor ///////////////////////////////////////////////

    static Shapeset::shape_fn_t** hdiv_leg_shape_fn_table[2] =
    {
      NULL,
      hdiv_leg_quad_shape_fn_table
    };

    static Shapeset::shape_fn_t** hdiv_leg_shape_fn_table_x[2] =
    {
      NULL,
      hdiv_leg_quad_shape_fn_table_x
    };

    static Shapeset::shape_fn_t** hdiv_leg_shape_fn_table_y[2] =
    {
      NULL,
      hdiv_leg_quad_shape_fn_table_y
    };

    static int* hdiv_leg_vertex_indices[2] =
    {
      NULL,
      hdiv_leg_quad_vertex_indices
    };

    static int** hdiv_leg_edge_indices[2] =
    {
      NULL,
      hdiv_leg_quad_edge_indices
    };

    static int** hdiv_leg_bubble_indices[2] =
    {
      NULL,
      hdiv_leg_quad_bubble_indices
    };

    static int* hdiv_leg_bubble_count[2] =
    {
      NULL,
      hdiv_leg_quad_bubble_count
    };

    static int* hdiv_leg_index_to_order[2] =
    {
      NULL,
      hdiv_leg_quad_index_to_order
    };

    int HdivShapesetLegendre::get_max_index(ElementMode2D mode) { return max_index[mode]; }

    HdivShapesetLegendre::HdivShapesetLegendre()
    {
      shape_table[0] = hdiv_leg_shape_fn_table;
      shape_table[1] = hdiv_leg_shape_fn_table_x;
      shape_table[2] = hdiv_leg_shape_fn_table_y;
      shape_table[3] = NULL;
      shape_table[4] = NULL;
      shape_table[5] = NULL;

      vertex_indices = hdiv_leg_vertex_indices;
      edge_indices = hdiv_leg_edge_indices;
      bubble_indices = hdiv_leg_bubble_indices;
      bubble_count = hdiv_leg_bubble_count;
      index_to_order = hdiv_leg_index_to_order;

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
      min_order = 0;
      num_components = 2;

      ebias = 0;  // TODO

      comb_table = NULL;
    }

    const int HdivShapesetLegendre::max_index[2] = { 149, 307 };
  }
}