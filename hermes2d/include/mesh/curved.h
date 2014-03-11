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

#include "../global.h"
#include "../shapeset/shapeset_common.h"
#include "../shapeset/precalc.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class Element;
    class RefMap;
    class Quad1DStd;
    class Quad2DStd;
    struct Trf;

    enum CurvType
    {
      ArcType,
      NurbsType
    };

    class HERMES_API Curve
    {
    public:
      Curve(CurvType type);
      virtual ~Curve();
      CurvType type;
    };

    class HERMES_API Arc : public Curve
    {
    public:
      Arc();
      Arc(double angle);
      Arc(const Arc* other);

      double angle; ///< arc angle

      /// Arc degree is 2.
      static const unsigned char degree = 2;

      /// Arc has 3 control points
      static const unsigned char np = 3;
      // there are 6 knots: {0, 0, 0, 1, 1, 1}
      static const unsigned char nk = 6;
      double kv[6];
      double3 pt[3];
    };

    /// \brief Represents one NURBS curve.
    ///
    /// The structure Nurbs defines one curved edge, or, more precisely,
    /// the control points and other data for one NURBS curve.
    ///
    /// A good introduction to NURBS curves can be found <a href=
    /// "http://devworld.apple.com/dev/techsupport/develop/issue25/schneider.html">
    /// here</a>.
    ///
    class HERMES_API Nurbs : public Curve
    {
    public:
      Nurbs();
      Nurbs(const Nurbs* other);
      ~Nurbs();

      unsigned char degree;  ///< curve degree (2=quadratic, etc.)
      unsigned char np;      ///< number of control points
      double3* pt; ///< control points and their weights
      unsigned char nk;      ///< knot vector length
      double* kv;  ///< knot vector
    };

    /// CurvMap is a structure storing complete information on the curved edges of
    /// an element. There are two variants of this structure. The first is for
    /// top-level (master mesh) elements.
    ///
    class HERMES_API CurvMap
    {
    public:
      CurvMap();
      CurvMap(const CurvMap* cm);
      ~CurvMap();
      void free();

      /// this structure defines a curved mapping of an element; it has two
      /// modes, depending on the value of 'toplevel'
      bool toplevel;
      union
      {
        // if toplevel=true, this structure belongs to a base mesh element
        // and the array 'nurbs' points to (up to four) NURBS curved edges
        Curve* curves[H2D_MAX_NUMBER_EDGES];
        struct
        {
          // if toplevel=false, this structure belongs to a refined element
          // and 'parent' points to the base mesh element CurvMap structure;
          Element* parent;
          uint64_t sub_idx;
        };
      };

      /// current polynomial degree of the refmap approximation
      int order;

    private:
      PrecalcShapesetAssembling ref_map_pss;

      /// Transformation (2x2) matrix.
      Trf* ctm;

      /// finally here are the coefficients of the higher-order basis functions
      /// that constitute the projected reference mapping:
      unsigned short nc; ///< number of coefficients
      double2* coeffs; ///< array of the coefficients

      /// this is called for every curvilinear element when it is created
      /// or when it is necessary to re-calculate coefficients for another
      /// order: 'e' is a pointer to the element to which this CurvMap
      /// belongs to. First, old "coeffs" are removed if they are not nullptr,
      /// then new_ coefficients are projected.
      void update_refmap_coeffs(Element* e);

      void get_mid_edge_points(Element* e, double2* pt, unsigned short n);

      /// Recursive calculation of the basis function N_i,k(int i, int k, double t, double* knot).
      static double nurbs_basis_fn(unsigned short i, unsigned short k, double t, double* knot);

      /// Nurbs curve: t goes from -1 to 1, function returns x, y coordinates in plane
      /// as well as the unit normal and unit tangential vectors. This is done using
      /// the Wikipedia page http://en.wikipedia.org/wiki/Non-uniform_rational_B-spline.
      static void nurbs_edge(Element* e, Curve* curve, int edge, double t, double& x,
        double& y);

      //// non-polynomial reference map //////////////////////////////////////////////////////////////////////////////////
      static const double2 ref_vert[2][H2D_MAX_NUMBER_VERTICES];

      /// Subtraction of straight edge and nurbs curve.
      static void nurbs_edge_0(Element* e, Curve* nurbs, unsigned short edge, double t, double& x, double& y, double& n_x, double& n_y, double& t_x, double& t_y);

      /// Calculation of nonpolynomial reference mapping on curved element
      static void calc_ref_map_tri(Element* e, Curve** nurbs, double xi_1, double xi_2, double& x, double& y);
      static void calc_ref_map_quad(Element* e, Curve** nurbs, double xi_1, double xi_2,
        double& x, double& y);

      static void calc_ref_map(Element* e, Curve** nurbs, double xi_1, double xi_2, double2& f);

      /// Edge part of projection based interpolation ///////////////////////////////////////////////////
      /// Compute point (x, y) in reference element, edge vector (v1, v2)
      void edge_coord(Element* e, unsigned short edge, double t, double2& x) const;
      void calc_edge_projection(Element* e, unsigned short edge, Curve** nurbs, unsigned short order, double2* proj) const;

      //// Bubble part of projection based interpolation /////////////////////////////////////////////////
      void old_projection(Element* e, unsigned short order, double2* proj, double* old[2]);
      void calc_bubble_projection(Element* e, Curve** nurbs, unsigned short order, double2* proj);

      static CurvMap* create_son_curv_map(Element* e, int son);

      template<typename T> friend class Space;
      template<typename T> friend class H1Space;
      template<typename T> friend class L2Space;
      template<typename T> friend class HcurlSpace;
      template<typename T> friend class HdivSpace;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class Adapt;
      template<typename T> friend class KellyTypeAdapt;
      friend class RefMap;
      friend class Mesh;
      friend class MeshReader;
      friend class MeshReaderH2D;
      friend class MeshReaderH2DXML;
      friend class MeshReaderH2DBSON;
    };

    class CurvMapStatic
    {
    public:
      CurvMapStatic();
      ~CurvMapStatic();

      /// Projection based interpolation
      /// Preparation of projection matrices, Cholesky factorization
      void precalculate_cholesky_projection_matrix_edge();
      /// Calculate the H1 seminorm products (\phi_i, \phi_j) for all 0 <= i, j < n, n is the number of bubble functions
      double** calculate_bubble_projection_matrix(unsigned short* indices, ElementMode2D mode);
      void precalculate_cholesky_projection_matrices_bubble();

      double** edge_proj_matrix;  ///< projection matrix for each edge is the same
      unsigned short edge_proj_matrix_size;
      double** bubble_proj_matrix_tri; ///< projection matrix for triangle bubbles
      double** bubble_proj_matrix_quad; ///< projection matrix for quad bubbles

      double* edge_p;  ///<  diagonal vector in cholesky factorization
      double* bubble_tri_p; ///<  diagonal vector in cholesky factorization
      unsigned short tri_bubble_np;
      double* bubble_quad_p; ///<  diagonal vector in cholesky factorization
      unsigned short quad_bubble_np;
    };

    /// Global instance used inside Hermes which is also accessible to users.
    extern CurvMapStatic curvMapStatic;
  }
}
#endif