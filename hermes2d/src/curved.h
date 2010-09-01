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

#ifndef __H2D_CURVED_H
#define __H2D_CURVED_H

#include "common.h"

struct Element;


/// \brief Represents one NURBS curve.
///
/// The structure Nurbs defines one curved edge, or, more precisely,
/// the control points and other data for one NURBS curve.
///
/// A good introduction to NURBS curves can be found <a href=
/// "http://devworld.apple.com/dev/techsupport/develop/issue25/schneider.html">
/// here</a>.
///
struct Nurbs
{
  Nurbs() { ref = 0; twin = false; };
  void unref();

  int degree;  ///< curve degree (2=quadratic, etc.)
  int np;      ///< number of control points
  double3* pt; ///< control points and their weights
  int nk;      ///< knot vector length
  double* kv;  ///< knot vector
  int ref;     ///< reference counter (the structure is deleted when this reaches zero)
  bool twin;   ///< true on internal curved edges for the second (artificial) Nurbs
  bool arc;     ///< true if this is in fact a circular arc
  double angle; ///< arc angle
};


/// CurvMap is a structure storing complete information on the curved edges of
/// an element. There are two variants of this structure. The first if for
/// top-level (master mesh) elements.
///
struct CurvMap
{
  CurvMap() { coefs = NULL; };
  CurvMap(CurvMap* cm);
  ~CurvMap();

  // this structure defines a curved mapping of an element; it has two
  // modes, depending on the value of 'toplevel'
  bool toplevel;
  union
  {
    // if toplevel=true, this structure belongs to a base mesh element
    // and the array 'nurbs' points to (up to four) NURBS curved edges
    Nurbs* nurbs[4];
    struct
    {
      // if toplevel=false, this structure belongs to a refined element
      // and 'parent' points to the base mesh element CurvMap structure;
      Element* parent;
      uint64_t part;
    };
  };

  // current polynomial degree of the refmap approximation
  int order;

  // finally here are the coefficients of the higher-order basis functions
  // that constitute the projected reference mapping:
  int nc; // number of coefficients (todo: mozna spis polyn. rad zobrazeni)
  double2* coefs; // array of the coefficients

  // this is called for every curvilinear element when it is created
  // or when it is necessary to re-calculate coefficients for another
  // order: 'e' is a pointer to the element to which this CurvMap
  // belongs to. First, old "coefs" are removed if they are not NULL,
  // then new coefficients are projected.
  void update_refmap_coefs(Element* e);

  void get_mid_edge_points(Element* e, double2* pt, int n);

};


void nurbs_edge(Element* e, Nurbs* nurbs, int edge, double t, double& x, double& y);


#endif
