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

#include "../hermes2d_common_defs.h"
#include "../shapeset/shapeset_common.h"
namespace Hermes
{
  namespace Hermes2D
  {
    class Element;
    class H1ShapesetJacobi;
    class PrecalcShapeset;
    class Quad1DStd;
    class Quad2DStd;
    struct Trf;

    /// \brief Represents one NURBS curve.
    ///
    /// The structure Nurbs defines one curved edge, or, more precisely,
    /// the control points and other data for one NURBS curve.
    ///
    /// A good introduction to NURBS curves can be found <a href=
    /// "http://devworld.apple.com/dev/techsupport/develop/issue25/schneider.html">
    /// here</a>.
    ///

    struct HERMES_API Nurbs
    {
      Nurbs()
      {
        ref = 0; twin = false;
      };
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
    /// an element. There are two variants of this structure. The first is for
    /// top-level (master mesh) elements.
    ///
    class HERMES_API CurvMap
    {
    public:
      CurvMap()
      {
        coeffs = NULL;};
        CurvMap(CurvMap* cm);
        ~CurvMap();

        /// this structure defines a curved mapping of an element; it has two
        /// modes, depending on the value of 'toplevel'
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

        /// current polynomial degree of the refmap approximation
        int order;

        /// finally here are the coefficients of the higher-order basis functions
        /// that constitute the projected reference mapping:
        int nc; // number of coefficients (todo: mozna spis polyn. rad zobrazeni)
        double2* coeffs; // array of the coefficients

        /// this is called for every curvilinear element when it is created
        /// or when it is necessary to re-calculate coefficients for another
        /// order: 'e' is a pointer to the element to which this CurvMap
        /// belongs to. First, old "coeffs" are removed if they are not NULL,
        /// then new coefficients are projected.
        void update_refmap_coeffs(Element* e);

        void get_mid_edge_points(Element* e, double2* pt, int n);

        static H1ShapesetJacobi ref_map_shapeset;
        static PrecalcShapeset ref_map_pss;

        static double** edge_proj_matrix;  ///< projection matrix for each edge is the same
        static double** bubble_proj_matrix_tri; ///< projection matrix for triangle bubbles
        static double** bubble_proj_matrix_quad; ///< projection matrix for quad bubbles

        static double* edge_p;  ///<  diagonal vector in cholesky factorization
        static double* bubble_tri_p; ///<  diagonal vector in cholesky factorization
        static double* bubble_quad_p; ///<  diagonal vector in cholesky factorization

        static Quad1DStd quad1d;
        static Quad2DStd quad2d; ///<  fixme: g_quad_2d_std

        static Trf ctm;

        /// Recursive calculation of the basis function N_i,k(int i, int k, double t, double* knot).
        static double nurbs_basis_fn(int i, int k, double t, double* knot);

        static void nurbs_edge(Element* e, Nurbs* nurbs, int edge, double t, double& x, 
          double& y, double& n_x, double& n_y, double& t_x, double& t_y);

        ///  Definition of vertex basis functions for triangles.
        static double lambda_0(double x, double y)
        {
          return -0.5 * (x + y);
        }
        static double lambda_1(double x, double y)
        {
          return  0.5 * (x + 1);
        }
        static double lambda_2(double x, double y)
        {
          return  0.5 * (y + 1);
        }

        ///  1D Lobatto functions.
        static double lob0(double x)
        {
          return l0(x);
        }
        ///  1D Lobatto functions.
        static double lob1(double x)
        {
          return l1(x);
        }
        ///  1D Lobatto functions.
        static double lob2(double x)
        {
          return l2(x);
        }
        ///  1D Lobatto functions.
        static double lob3(double x)
        {
          return l3(x);
        }
        ///  1D Lobatto functions.
        static double lob4(double x)
        {
          return l4(x);
        }
        ///  1D Lobatto functions.
        static double lob5(double x)
        {
          return l5(x);
        }
        ///  1D Lobatto functions.
        static double lob6(double x)
        {
          return l6(x);
        }
        ///  1D Lobatto functions.
        static double lob7(double x)
        {
          return l7(x);
        }
        ///  1D Lobatto functions.
        static double lob8(double x)
        {
          return l8(x);
        }
        ///  1D Lobatto functions.
        static double lob9(double x)
        {
          return l9(x);
        }
        ///  1D Lobatto functions.
        static double lob10(double x)
        {
          return l10(x);
        }
        ///  1D Lobatto functions.
        static double lob11(double x)
        {
          return l11(x);
        }

        static const double2 ref_vert[2][4];

        /// Subtraction of straight edge and nurbs curve.
        static void nurbs_edge_0(Element* e, Nurbs* nurbs, int edge, double t, double& x, double& y, double& n_x, double& n_y, double& t_x, double& t_y);
        static void calc_ref_map_tri(Element* e, Nurbs** nurbs, double xi_1, double xi_2, double& x, double& y);
        static void calc_ref_map_quad(Element* e, Nurbs** nurbs, double xi_1, double xi_2,
          double& x, double& y);

        static void calc_ref_map(Element* e, Nurbs** nurbs, double xi_1, double xi_2, double2& f);

        static void precalculate_cholesky_projection_matrix_edge();
        static double** calculate_bubble_projection_matrix(int nb, int* indices);
        static void precalculate_cholesky_projection_matrices_bubble();

        static void edge_coord(Element* e, int edge, double t, double2& x, double2& v);
        static void calc_edge_projection(Element* e, int edge, Nurbs** nurbs, int order, double2* proj);

        static void old_projection(Element* e, int order, double2* proj, double* old[2]);
        static void calc_bubble_projection(Element* e, Nurbs** nurbs, int order, double2* proj);

        static void ref_map_projection(Element* e, Nurbs** nurbs, int order, double2* proj);

        static bool warning_issued;
    };
  }
}
#endif
