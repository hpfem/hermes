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

#ifndef __H2D_NURBS_H
#define __H2D_NURBS_H

#include "nurbs_matrix.h"
//#include "../../../hermes_common/matrix.h"

/*
  Structures needed for NURBS. 
  Things that are commented out are being developed. I want to 
  overload arithmetic operators and use them on my Point2D and HPoint2D objects.
*/
  struct Point2D {
    double x, y;
/*
    Point2D& operator + (Point2D&);
    Point2D& operator - (Point2D&);
    Point2D& operator* (double param);
    Point2D& operator/= (double param);
    Point2D& operator = (const Point2D&);
    Point2D& operator= (double param);
*/
  };
/*
  Point2D& Point2D::operator+ (Point2D& param) {
  x = x + param.x;
  y = y + param.y;
  return (*this);
  }
 
  Point2D& Point2D::operator- (Point2D& param) {
  x = x - param.x;
  y = y - param.y;
  return (*this);
  }

  Point2D& Point2D::operator* (double param) {
  x = x*param;
  y = y*param;
  return (*this);
  }

  Point2D& Point2D::operator/= (double param) {
  x = x/param;
  y = y/param;
  return (*this);
  }

  Point2D& Point2D::operator= (const Point2D& param) {
  x = param.x;
  y = param.y;
  return (*this);
  }

  Point2D& Point2D::operator= (double param) {
  x = param;
  y = param;
  return (*this);
  }
  
*/
  struct HPoint2D {
    double wx, wy, w;
/*
    HPoint2D& operator + (HPoint2D&);
    HPoint2D& operator - (HPoint2D&);
    HPoint2D& operator* (double param);
    HPoint2D& operator/ (double param);
    HPoint2D& operator = (const HPoint2D&);
    HPoint2D& operator= (double param);
*/
  };

/*  
  HPoint2D& HPoint2D::operator+ (HPoint2D& param) {
  wx = wx + param.wx;
  wy = wy + param.wy;
  w  = w + param.w;
  return (*this);
  }
 
  HPoint2D& HPoint2D::operator- (HPoint2D& param) {
  wx = wx - param.wx;
  wy = wy - param.wy;
  w  = w - param.w;
  return (*this);
  }

  HPoint2D& HPoint2D::operator* (double param) {
  wx = wx*param;
  wy = wy*param;
  w  = w*param;
  return (*this);
  }

  HPoint2D& HPoint2D::operator/ (double param) {
  wx = wx/param;
  wy = wy/param;
  w  = w/param;
  return (*this);
  }

  HPoint2D& HPoint2D::operator= (const HPoint2D& param) {
  wx = param.wx;
  wy = param.wy;
  w  = param.w;
  return (*this);
  }

  HPoint2D& HPoint2D::operator= (double param) {
  wx = param;
  wy = param;
  w  = param;
  return (*this);
  }

*/


/*
    A NURBS curve class
    
    This class is used to represent and manipulate NURBS curve. 
    The curves can have any degree and have any number of control points.
*/
  class NurbsCurve {

    public: 
    NurbsCurve();
    NurbsCurve(int deg, int CPtMaxIndx, double Px[], double Py[], double w[], int KntMaxIndx, double Uvals[]);
    ~NurbsCurve();
    Point2D CurvePoint(double u) const;
    Point2D CurvePoint(double u, int span, double N[]) const;
    HPoint2D CurvePointH(double u, int span) const;
    int FindSpan(double u) const;
    void BasisFuns(int span, double u, double N[]) const;
    void nurbsBasisFuns(int span, double u, int p, double U[], double N[]) const;
    void DersBasisFuns(int d, double u, int span, Matrix<double> &ders) const;
    double SingleBasisFun( double u, int i, int p) const;
    void AllBasisFuns(int span, double u, Matrix<double> &N) const;
    // double DersSingleBasisFun(double u, int i, int p) const;

    // Derivatives.
    void deriveAtH(double u, int d, struct HPoint2D ders[]) const;
    void deriveAtH(double u, int d, int span, struct HPoint2D ders[]) const;
    void DerivCpts(int d, int r1, int r2, Matrix<HPoint2D> &PK) const;
    void deriveAtH2(double u, int d, struct HPoint2D ders[]) const;  
    void deriveAt(double u, int d, struct Point2D ders[]) const;
    HPoint2D firstD(double u) const; 
    HPoint2D firstD(double u, int span) const;
    Point2D firstDn(double u) const;

    void FirstD_at_ends(struct Point2D Cd[]) const;
    Point2D normal_at(double u) const;
    Point2D normal_at_end(int which_end) const;


    private:
    int p; // Degree of the NURBS curve.
    int n; // The highest index of control points.
    HPoint2D *P; // Array of control points. Length n+1.
    int m; // The highest index of knot.
    double *U; // Knot vector. An array of length m+1.

  };

  // Misc. functions.
  void binomialCoef(int d, Matrix<double> &Bin);
  Point2D project(HPoint2D &Cw);

#endif

