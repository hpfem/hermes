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

#ifndef __H2D_SHAPESET_COMMON_H
#define __H2D_SHAPESET_COMMON_H

// Common definitions used by the shapesets...


#define lambda1(x,y) (((y) + 1) / 2)
#define lambda2(x,y) (-((x) + (y)) / 2)
#define lambda3(x,y) (((x) + 1) / 2)

// x derivatives of affine coordinates
#define lambda1x(x,y) (0.0)
#define lambda2x(x,y) (-1.0 / 2.0)
#define lambda3x(x,y) (1.0 / 2.0)

// y derivatives of affine coordinates
#define lambda1y(x,y) (1.0 / 2.0)
#define lambda2y(x,y) (-1.0 / 2.0)
#define lambda3y(x,y) (0.0)

// kernel functions
#define phi0(x) (-2.0 * 1.22474487139158904909864203735)
#define phi1(x) (-2.0 * 1.58113883008418966599944677222 * (x))
#define phi2(x) (-1.0 / 2.0 * 1.87082869338697069279187436616 * (5 * (x) * (x) - 1))
#define phi3(x) (-1.0 / 2.0 * 2.12132034355964257320253308631 * (7 * (x) * (x) - 3) * (x))
#define phi4(x) (-1.0 / 4.0 * 2.34520787991171477728281505677 * (21 * (x) * (x) * (x) * (x) - 14 * (x) * (x) + 1))
#define phi5(x) (-1.0 / 4.0 * 2.54950975679639241501411205451 * ((33 * (x) * (x) - 30) * (x) * (x) + 5) * (x))
#define phi6(x) (-1.0 / 32.0 * 2.73861278752583056728484891400 * (((429 * (x) * (x) - 495) * (x) * (x) + 135) * (x) * (x) - 5))
#define phi7(x) (-1.0 / 32.0 * 2.91547594742265023543707643877 * (((715 * (x) * (x) - 1001) * (x) * (x) + 385) * (x) * (x) - 35) * (x))
#define phi8(x) (-1.0 / 64.0 * 3.08220700148448822512509619073 * ((((2431 * (x) * (x) - 4004) * (x) * (x) + 2002) * (x) * (x) - 308) * (x) * (x) + 7))
#define phi9(x) (-1.0 / 128.0 * 6.4807406984078603784382721642 * ((((4199 * (x) * (x) - 7956) * (x) * (x) + 4914) * (x) * (x) - 1092) * (x) * (x) + 63) * (x))

// derivatives of kernel functions
#define phi0x(x) (0)
#define phi1x(x) (-2.0 * 1.58113883008418966599944677222)
#define phi2x(x) (-1.0 / 2.0 * 1.87082869338697069279187436616 * (10 * (x)))
#define phi3x(x) (-1.0 / 2.0 * 2.12132034355964257320253308631 * (21.0*(x)*(x)-3.0))
#define phi4x(x) (-1.0 / 4.0 * 2.34520787991171477728281505677 * ((84.0*(x)*(x)-28.0)*(x)))
#define phi5x(x) (-1.0 / 4.0 * 2.54950975679639241501411205451 * ((165.0*(x)*(x)-90.0)*(x)*(x)+5.0))
#define phi6x(x) (-1.0 / 32.0 * 2.73861278752583056728484891400 * (((2574.0*(x)*(x)-1980.0)*(x)*(x)+270.0)*(x)))
#define phi7x(x) (-1.0 / 32.0 * 2.91547594742265023543707643877 * (((5005.0*(x)*(x)-5005.0)*(x)*(x)+1155.0)*(x)*(x)-35.0))
#define phi8x(x) (-1.0 / 64.0 * 3.08220700148448822512509619073 * ((((19448.0*(x)*(x)-24024.0)*(x)*(x)+8008.0)*(x)*(x)-616.0)*(x)))
#define phi9x(x) (-1.0 / 128.0 * 6.4807406984078603784382721642 * ((((37791.0*(x)*(x)-55692.0)*(x)*(x)+24570.0)*(x)*(x)-3276.0)*(x)*(x)-63.0))

// second derivatives of kernel functions
#define phi0xx(x) (0)
#define phi1xx(x) (0)
#define phi2xx(x) (-1.0 / 2.0 * 1.87082869338697069279187436616 * 10)
#define phi3xx(x) (-1.0 / 2.0 * 2.12132034355964257320253308631 * (42.0*(x)))
#define phi4xx(x) (-1.0 / 4.0 * 2.34520787991171477728281505677 * ((252.0*(x)*(x)-28.0)))
#define phi5xx(x) (-1.0 / 4.0 * 2.54950975679639241501411205451 * ((660.0*(x)*(x)*(x)-180.0*(x))))
#define phi6xx(x) (-1.0 / 32.0 * 2.73861278752583056728484891400 * (((12870.0*(x)*(x)-5940.0)*(x)*(x)+270.0)))
#define phi7xx(x) (-1.0 / 32.0 * 2.91547594742265023543707643877 * (((30030.0*(x)*(x)-20020.0)*(x)*(x)+2310.0)*(x)))
#define phi8xx(x) (-1.0 / 64.0 * 3.08220700148448822512509619073 * ((((136136.0*(x)*(x)-120120.0)*(x)*(x)+24024.0)*(x)*(x)-616.0)))
#define phi9xx(x) (-1.0 / 128.0 * 6.4807406984078603784382721642 * ((((302328.0*(x)*(x)-334152.0)*(x)*(x)+98280.0)*(x)*(x)-6552.0)*(x)*(x)))

// Legendre polynomials
#define Legendre0(x) (1.0)
#define Legendre1(x) (x)
#define Legendre2(x) (1.0 / 2.0 * (3 * (x) * (x) - 1))
#define Legendre3(x) (1.0 / 2.0 * (5 * (x) * (x) - 3) * (x))
#define Legendre4(x) (1.0 / 8.0 * ((35 * (x) * (x) - 30) * (x) * (x) + 3))
#define Legendre5(x) (1.0 / 8.0 * ((63 * (x) * (x) - 70) * (x) * (x) + 15) * (x))
#define Legendre6(x) (1.0 / 16.0 * (((231 * (x) * (x) - 315) * (x) * (x) + 105) * (x) * (x) - 5))
#define Legendre7(x) (1.0 / 16.0 * (((429 * (x) * (x) - 693) * (x) * (x) + 315) * (x) * (x) - 35) * (x))
#define Legendre8(x) (1.0 / 128.0 * ((((6435 * (x) * (x) - 12012) * (x) * (x) + 6930) * (x) * (x) - 1260) * (x) * (x) + 35))
#define Legendre9(x) (1.0 / 128.0 * ((((12155 * (x) * (x) - 25740) * (x) * (x) + 18018) * (x) * (x) - 4620) * (x) * (x) + 315) * (x))
#define Legendre10(x) (1.0 / 256.0 * (((((46189 * (x) * (x) - 109395) * (x) * (x) + 90090) * (x) * (x) - 30030) * (x) * (x) + 3465) * (x) * (x) - 63))

// derivatives of Legendre polynomials
#define Legendre0x(x) (0.0)
#define Legendre1x(x) (1.0)
#define Legendre2x(x) (3.0 * (x))
#define Legendre3x(x) (15.0 / 2.0 * (x) * (x) - 3.0 / 2.0)
#define Legendre4x(x) (5.0 / 2.0 * (x) * (7.0 * (x) * (x) - 3.0))
#define Legendre5x(x) ((315.0 / 8.0 * (x) * (x) - 105.0 / 4.0) * (x) * (x) + 15.0 / 8.0)
#define Legendre6x(x) (21.0 / 8.0 * (x) * ((33.0 * (x) * (x) - 30.0) * (x) * (x) + 5.0))
#define Legendre7x(x) (((3003.0 / 16.0 * (x) * (x) - 3465.0 / 16.0) * (x) * (x) + 945.0 / 16.0) * (x) * (x) - 35.0 / 16.0)
#define Legendre8x(x) (9.0 / 16.0 * (x) * (((715.0 * (x) * (x) - 1001.0) * (x) * (x) + 385.0) * (x) * (x) - 35.0))
#define Legendre9x(x) ((((109395.0 / 128.0 * (x) * (x) - 45045.0 / 32.0) * (x) * (x) + 45045.0 / 64.0) * (x) * (x) - 3465.0 / 32.0) * (x) * (x) + 315.0 / 128.0)
#define Legendre10x(x) (2.0 / 256.0 * (x) * ((((230945.0 * (x) * (x) - 437580.0) * (x) * (x) + 270270.0) * (x) * (x) - 60060.0) * (x) * (x) + 3465.0))

// second derivatives of Legendre polynomials
#define Legendre0xx(x) (0.0)
#define Legendre1xx(x) (0.0)
#define Legendre2xx(x) (3.0)
#define Legendre3xx(x) (15.0 * (x))
#define Legendre4xx(x) (105.0 / 2.0 * (x) * (x) - 15.0 / 2.0)
#define Legendre5xx(x) (105.0 / 2.0 * (x) * (3.0 * (x) * (x) - 1.0))
#define Legendre6xx(x) ((3465.0 / 8.0 * (x) * (x) - 1890.0 / 8.0) * (x) * (x) + 105.0 / 8.0)
#define Legendre7xx(x) (63.0 / 8.0 * (x) * ((143.0 * (x) * (x) - 110.0) * (x) * (x) + 15.0))
#define Legendre8xx(x) (((45045.0 / 16.0 * (x) * (x) - 45045.0 / 16.0) * (x) * (x) + 10395.0 / 16.0) * (x) * (x) - 315.0 / 16.0)
#define Legendre9xx(x) (45.0 / 16.0 * (x) * (((2431.0 * (x) * (x) - 3003.0) * (x) * (x) + 1001.0) * (x) * (x) - 77.0))
#define Legendre10xx(x) ((((2078505.0 / 128.0 * (x) * (x) - 765765.0 / 32.0) * (x) * (x) + 675675.0 / 64.0) * (x) * (x) - 45045.0 / 32.0) * (x) * (x) + 3465.0 / 128.0)

// first two Lobatto shape functions
#define l0(x) ((1.0 - (x)) * 0.5)
#define l1(x) ((1.0 + (x)) * 0.5)

#define l0l1(x) ((1.0 - (x)*(x)) * 0.25)

// other Lobatto shape functions
#define l2(x)  (phi0(x) * l0l1(x))
#define l3(x)  (phi1(x) * l0l1(x))
#define l4(x)  (phi2(x) * l0l1(x))
#define l5(x)  (phi3(x) * l0l1(x))
#define l6(x)  (phi4(x) * l0l1(x))
#define l7(x)  (phi5(x) * l0l1(x))
#define l8(x)  (phi6(x) * l0l1(x))
#define l9(x)  (phi7(x) * l0l1(x))
#define l10(x) (phi8(x) * l0l1(x))
#define l11(x) (phi9(x) * l0l1(x))

// derivatives of Lobatto functions
#define dl0(x)  (-0.5)
#define dl1(x)  (0.5)
#define dl2(x)  (sqrt(3.0/2.0) * Legendre1(x))
#define dl3(x)  (sqrt(5.0/2.0) * Legendre2(x))
#define dl4(x)  (sqrt(7.0/2.0) * Legendre3(x))
#define dl5(x)  (sqrt(9.0/2.0) * Legendre4(x))
#define dl6(x)  (sqrt(11.0/2.0) * Legendre5(x))
#define dl7(x)  (sqrt(13.0/2.0) * Legendre6(x))
#define dl8(x)  (sqrt(15.0/2.0) * Legendre7(x))
#define dl9(x)  (sqrt(17.0/2.0) * Legendre8(x))
#define dl10(x) (sqrt(19.0/2.0) * Legendre9(x))
#define dl11(x) (sqrt(21.0/2.0) * Legendre10(x))

// second derivatives of Lobatto functions
#define d2l0(x)  (0.0)
#define d2l1(x)  (0.0)
#define d2l2(x)  (sqrt(3.0/2.0) * Legendre1x(x))
#define d2l3(x)  (sqrt(5.0/2.0) * Legendre2x(x))
#define d2l4(x)  (sqrt(7.0/2.0) * Legendre3x(x))
#define d2l5(x)  (sqrt(9.0/2.0) * Legendre4x(x))
#define d2l6(x)  (sqrt(11.0/2.0) * Legendre5x(x))
#define d2l7(x)  (sqrt(13.0/2.0) * Legendre6x(x))
#define d2l8(x)  (sqrt(15.0/2.0) * Legendre7x(x))
#define d2l9(x)  (sqrt(17.0/2.0) * Legendre8x(x))
#define d2l10(x) (sqrt(19.0/2.0) * Legendre9x(x))
#define d2l11(x) (sqrt(21.0/2.0) * Legendre10x(x))


//// Hcurl /////////////////////////////////////////////////////////////////////////////////////////

// outward unit normal vectors to the reference triangle
#define n11  0.0
#define n12 -1.0
#define n21  0.707106781186547524400844362105
#define n22  0.707106781186547524400844362105
#define n31 -1.0
#define n32  0.0

// inner products of normals and tangents
#define n2t1  0.707106781186547524400844362105
#define n3t1 -1.0
#define n3t2  0.707106781186547524400844362105
#define n1t2 -0.707106781186547524400844362105
#define n1t3  1.0
#define n2t3 -0.707106781186547524400844362105

/*
// outward unit normal vectors to the reference triangle
#define n11  0.0
#define n12 -1.0
#define n21  1.0
#define n22  1.0
#define n31 -1.0
#define n32  0.0

// inner products of normals and tangents
#define n2t1  1.0
#define n3t1 -1.0
#define n3t2  1.0
#define n1t2 -1.0
#define n1t3  1.0
#define n2t3 -1.0
*/
/*
// Whitney functions
#define psi0e1_1(x,y) ((lambda3(x,y) * n21 / n2t1 + lambda2(x,y) * n31 / n3t1)/ 2.0)
#define psi0e1_2(x,y) ((lambda3(x,y) * n22 / n2t1 + lambda2(x,y) * n32 / n3t1)/ 2.0)
#define psi0e2_1(x,y) ((lambda1(x,y) * n31 / n3t2 + lambda3(x,y) * n11 / n1t2)/ 2.82842712474619)
#define psi0e2_2(x,y) ((lambda1(x,y) * n32 / n3t2 + lambda3(x,y) * n12 / n1t2)/ 2.82842712474619)
#define psi0e3_1(x,y) ((lambda2(x,y) * n11 / n1t3 + lambda1(x,y) * n21 / n2t3)/ 2.0)
#define psi0e3_2(x,y) ((lambda2(x,y) * n12 / n1t3 + lambda1(x,y) * n22 / n2t3)/ 2.0)

// x derivatives of Whitney functions
#define psi0e1x_1(x,y) ((lambda3x(x,y) * n21 / n2t1 + lambda2x(x,y) * n31 / n3t1)/2.0)
#define psi0e1x_2(x,y) ((lambda3x(x,y) * n22 / n2t1 + lambda2x(x,y) * n32 / n3t1)/2.0)
#define psi0e2x_1(x,y) ((lambda1x(x,y) * n31 / n3t2 + lambda3x(x,y) * n11 / n1t2)/ 2.82842712474619)
#define psi0e2x_2(x,y) ((lambda1x(x,y) * n32 / n3t2 + lambda3x(x,y) * n12 / n1t2)/ 2.82842712474619)
#define psi0e3x_1(x,y) ((lambda2x(x,y) * n11 / n1t3 + lambda1x(x,y) * n21 / n2t3)/2.0)
#define psi0e3x_2(x,y) ((lambda2x(x,y) * n12 / n1t3 + lambda1x(x,y) * n22 / n2t3)/2.0)

// y derivatives of Whitney functions
#define psi0e1y_1(x,y) ((lambda3y(x,y) * n21 / n2t1 + lambda2y(x,y) * n31 / n3t1)/2.0)
#define psi0e1y_2(x,y) ((lambda3y(x,y) * n22 / n2t1 + lambda2y(x,y) * n32 / n3t1)/2.0)
#define psi0e2y_1(x,y) ((lambda1y(x,y) * n31 / n3t2 + lambda3y(x,y) * n11 / n1t2)/ 2.82842712474619)
#define psi0e2y_2(x,y) ((lambda1y(x,y) * n32 / n3t2 + lambda3y(x,y) * n12 / n1t2)/ 2.82842712474619)
#define psi0e3y_1(x,y) ((lambda2y(x,y) * n11 / n1t3 + lambda1y(x,y) * n21 / n2t3)/2.0)
#define psi0e3y_2(x,y) ((lambda2y(x,y) * n12 / n1t3 + lambda1y(x,y) * n22 / n2t3)/2.0)

// linear edge functions
#define psi1e1_1(x,y) ((lambda3(x,y) * n21 / n2t1 - lambda2(x,y) * n31 / n3t1)/2.0)
#define psi1e1_2(x,y) ((lambda3(x,y) * n22 / n2t1 - lambda2(x,y) * n32 / n3t1)/2.0)
#define psi1e2_1(x,y) ((lambda1(x,y) * n31 / n3t2 - lambda3(x,y) * n11 / n1t2)/ 2.82842712474619)
#define psi1e2_2(x,y) ((lambda1(x,y) * n32 / n3t2 - lambda3(x,y) * n12 / n1t2)/ 2.82842712474619)
#define psi1e3_1(x,y) ((lambda2(x,y) * n11 / n1t3 - lambda1(x,y) * n21 / n2t3)/2.0)
#define psi1e3_2(x,y) ((lambda2(x,y) * n12 / n1t3 - lambda1(x,y) * n22 / n2t3)/2.0)

// x derivatives of linear edge functions
#define psi1e1x_1(x,y) ((lambda3x(x,y) * n21 / n2t1 - lambda2x(x,y) * n31 / n3t1)/2.0)
#define psi1e1x_2(x,y) ((lambda3x(x,y) * n22 / n2t1 - lambda2x(x,y) * n32 / n3t1)/2.0)
#define psi1e2x_1(x,y) ((lambda1x(x,y) * n31 / n3t2 - lambda3x(x,y) * n11 / n1t2)/ 2.82842712474619)
#define psi1e2x_2(x,y) ((lambda1x(x,y) * n32 / n3t2 - lambda3x(x,y) * n12 / n1t2)/ 2.82842712474619)
#define psi1e3x_1(x,y) ((lambda2x(x,y) * n11 / n1t3 - lambda1x(x,y) * n21 / n2t3)/2.0)
#define psi1e3x_2(x,y) ((lambda2x(x,y) * n12 / n1t3 - lambda1x(x,y) * n22 / n2t3)/2.0)

// y derivatives of linear edge functions
#define psi1e1y_1(x,y) ((lambda3y(x,y) * n21 / n2t1 - lambda2y(x,y) * n31 / n3t1)/2.0)
#define psi1e1y_2(x,y) ((lambda3y(x,y) * n22 / n2t1 - lambda2y(x,y) * n32 / n3t1)/2.0)
#define psi1e2y_1(x,y) ((lambda1y(x,y) * n31 / n3t2 - lambda3y(x,y) * n11 / n1t2)/ 2.82842712474619)
#define psi1e2y_2(x,y) ((lambda1y(x,y) * n32 / n3t2 - lambda3y(x,y) * n12 / n1t2)/ 2.82842712474619)
#define psi1e3y_1(x,y) ((lambda2y(x,y) * n11 / n1t3 - lambda1y(x,y) * n21 / n2t3)/2.0)
#define psi1e3y_2(x,y) ((lambda2y(x,y) * n12 / n1t3 - lambda1y(x,y) * n22 / n2t3)/2.0)
*/

// Whitney functions
#define psi0e1_1(x,y) ((lambda3(x,y) * n21 / n2t1 + lambda2(x,y) * n31 / n3t1))
#define psi0e1_2(x,y) ((lambda3(x,y) * n22 / n2t1 + lambda2(x,y) * n32 / n3t1))
#define psi0e2_1(x,y) ((lambda1(x,y) * n31 / n3t2 + lambda3(x,y) * n11 / n1t2))
#define psi0e2_2(x,y) ((lambda1(x,y) * n32 / n3t2 + lambda3(x,y) * n12 / n1t2))
#define psi0e3_1(x,y) ((lambda2(x,y) * n11 / n1t3 + lambda1(x,y) * n21 / n2t3))
#define psi0e3_2(x,y) ((lambda2(x,y) * n12 / n1t3 + lambda1(x,y) * n22 / n2t3))

// x derivatives of Whitney functions
#define psi0e1x_1(x,y) ((lambda3x(x,y) * n21 / n2t1 + lambda2x(x,y) * n31 / n3t1))
#define psi0e1x_2(x,y) ((lambda3x(x,y) * n22 / n2t1 + lambda2x(x,y) * n32 / n3t1))
#define psi0e2x_1(x,y) ((lambda1x(x,y) * n31 / n3t2 + lambda3x(x,y) * n11 / n1t2))
#define psi0e2x_2(x,y) ((lambda1x(x,y) * n32 / n3t2 + lambda3x(x,y) * n12 / n1t2))
#define psi0e3x_1(x,y) ((lambda2x(x,y) * n11 / n1t3 + lambda1x(x,y) * n21 / n2t3))
#define psi0e3x_2(x,y) ((lambda2x(x,y) * n12 / n1t3 + lambda1x(x,y) * n22 / n2t3))

// y derivatives of Whitney functions
#define psi0e1y_1(x,y) ((lambda3y(x,y) * n21 / n2t1 + lambda2y(x,y) * n31 / n3t1))
#define psi0e1y_2(x,y) ((lambda3y(x,y) * n22 / n2t1 + lambda2y(x,y) * n32 / n3t1))
#define psi0e2y_1(x,y) ((lambda1y(x,y) * n31 / n3t2 + lambda3y(x,y) * n11 / n1t2))
#define psi0e2y_2(x,y) ((lambda1y(x,y) * n32 / n3t2 + lambda3y(x,y) * n12 / n1t2))
#define psi0e3y_1(x,y) ((lambda2y(x,y) * n11 / n1t3 + lambda1y(x,y) * n21 / n2t3))
#define psi0e3y_2(x,y) ((lambda2y(x,y) * n12 / n1t3 + lambda1y(x,y) * n22 / n2t3))

// linear edge functions
#define psi1e1_1(x,y) ((lambda3(x,y) * n21 / n2t1 - lambda2(x,y) * n31 / n3t1))
#define psi1e1_2(x,y) ((lambda3(x,y) * n22 / n2t1 - lambda2(x,y) * n32 / n3t1))
#define psi1e2_1(x,y) ((lambda1(x,y) * n31 / n3t2 - lambda3(x,y) * n11 / n1t2))
#define psi1e2_2(x,y) ((lambda1(x,y) * n32 / n3t2 - lambda3(x,y) * n12 / n1t2))
#define psi1e3_1(x,y) ((lambda2(x,y) * n11 / n1t3 - lambda1(x,y) * n21 / n2t3))
#define psi1e3_2(x,y) ((lambda2(x,y) * n12 / n1t3 - lambda1(x,y) * n22 / n2t3))

// x derivatives of linear edge functions
#define psi1e1x_1(x,y) ((lambda3x(x,y) * n21 / n2t1 - lambda2x(x,y) * n31 / n3t1))
#define psi1e1x_2(x,y) ((lambda3x(x,y) * n22 / n2t1 - lambda2x(x,y) * n32 / n3t1))
#define psi1e2x_1(x,y) ((lambda1x(x,y) * n31 / n3t2 - lambda3x(x,y) * n11 / n1t2))
#define psi1e2x_2(x,y) ((lambda1x(x,y) * n32 / n3t2 - lambda3x(x,y) * n12 / n1t2))
#define psi1e3x_1(x,y) ((lambda2x(x,y) * n11 / n1t3 - lambda1x(x,y) * n21 / n2t3))
#define psi1e3x_2(x,y) ((lambda2x(x,y) * n12 / n1t3 - lambda1x(x,y) * n22 / n2t3))

// y derivatives of linear edge functions
#define psi1e1y_1(x,y) ((lambda3y(x,y) * n21 / n2t1 - lambda2y(x,y) * n31 / n3t1))
#define psi1e1y_2(x,y) ((lambda3y(x,y) * n22 / n2t1 - lambda2y(x,y) * n32 / n3t1))
#define psi1e2y_1(x,y) ((lambda1y(x,y) * n31 / n3t2 - lambda3y(x,y) * n11 / n1t2))
#define psi1e2y_2(x,y) ((lambda1y(x,y) * n32 / n3t2 - lambda3y(x,y) * n12 / n1t2))
#define psi1e3y_1(x,y) ((lambda2y(x,y) * n11 / n1t3 - lambda1y(x,y) * n21 / n2t3))
#define psi1e3y_2(x,y) ((lambda2y(x,y) * n12 / n1t3 - lambda1y(x,y) * n22 / n2t3))


#endif
