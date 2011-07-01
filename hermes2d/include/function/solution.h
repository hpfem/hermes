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

#ifndef __H2D_SOLUTION_H
#define __H2D_SOLUTION_H

#include "mesh_function.h"

namespace Hermes
{
  class Ord;

  namespace Hermes2D
  {
    static double3* cheb_tab_tri[11];
    static double3* cheb_tab_quad[11];
    static int      cheb_np_tri[11];
    static int      cheb_np_quad[11];

    extern double3** cheb_tab[2];
    extern int*      cheb_np[2];

    /// Quad2DCheb is a special "quadrature" consisting of product Chebyshev
    /// points on the reference triangle and quad. It is used for expressing
    /// the solution on an element as a linear combination of monomials.
    ///
    class Quad2DCheb : public Quad2D
    {
    public:

      Quad2DCheb()
      {
        mode = HERMES_MODE_TRIANGLE;
        max_order[0]  = max_order[1]  = 10;
        num_tables[0] = num_tables[1] = 11;
        tables = cheb_tab;
        np = cheb_np;

        tables[0][0] = tables[1][0] = NULL;
        np[0][0] = np[1][0] = 0;

        int i, j, k, n, m;
        double3* pt;
        for (mode = 0; mode <= 1; mode++)
        {
          for (k = 0; k <= 10; k++)
          {
            np[mode][k] = n = mode ? sqr(k+1) : (k+1)*(k+2)/2;
            tables[mode][k] = pt = new double3[n];

            for (i = k, m = 0; i >= 0; i--)
              for (j = k; j >= (mode ? 0 : k-i); j--, m++) 
              {
                pt[m][0] = k ? cos(j * M_PI / k) : 1.0;
                pt[m][1] = k ? cos(i * M_PI / k) : 1.0;
                pt[m][2] = 1.0;
              }
          }
        }
      };

      ~Quad2DCheb()
      {
        for (int mode = 0; mode <= 1; mode++)
          for (int k = 0; k <= 10; k++){
            if (tables[mode][k])
            {
              delete[] tables[mode][k];
              tables[mode][k]=NULL;
            }
          }
      }

      virtual void dummy_fn() {}

    };

    enum SolutionType {
      HERMES_UNDEF = -1,
      HERMES_SLN = 0,
      HERMES_EXACT = 1,
      HERMES_CONST = 2
    };

    /// \brief Represents the solution of a PDE.
    ///
    /// The Solution class represents the solution of a PDE. Given a space and a solution vector,
    /// it calculates the appropriate linear combination of basis functions at the specified
    /// element and integration points.

    ///  The higher-order solution on elements is best calculated not as a linear  combination
    ///  of shape functions (the usual approach), but as a linear combination of monomials.
    ///  This has the advantage that no shape function table calculation and look-ups are
    ///  necessary (except for the conversion of the solution coefficients), which means that
    ///  visualization and multi-mesh calculations are much faster (all the push_transforms
    ///  and table searches take the most time when evaluating the solution).
    ///
    ///  The linear combination of monomials can be calculated using the Horner's scheme, which
    ///  requires the same number of multiplications as the calculation of the linear combination
    ///  of shape functions. However, sub-element transforms are trivial and cheap. Moreover,
    ///  after the solution on all elements is expressed as a combination of monomials, the
    ///  Space can be forgotten. This is comfortable for the user, since the Solution class acts
    ///  as a self-contained unit, internally containing just a copy of the mesh and a table of
    ///  monomial coefficients. It is also very straight-forward to obtain all derivatives of
    ///  a solution defined in this way. Finally, it is possible to store the solution on the
    ///  disk easily (no need to store the Space, which is difficult).
    ///
    ///  The following is an example of the set of monomials for a cubic quad and a cubic triangle.
    ///  (Note that these are actually the definitions of the polynomial spaces on these elements.)
    ///
    ///    [ x^3*y^3  x^2*y^3  x*y^3  y^3 ]       [                    y^3 ]
    ///    [ x^3*y^2  x^2*y^2  x*y^2  y^2 ]       [             x*y^2  y^2 ]
    ///    [ x^3*y    x^2*y    x*y    y   ]       [      x^2*y  x*y    y   ]
    ///    [ x^3      x^2      x      1   ]       [ x^3  x^2    x      1   ]
    ///
    ///  (The number of monomials is (n+1)^2 for quads and (n+1)*(n+2)/2 for triangles, where
    ///   'n' is the polynomial degree.)
    ///
    template<typename Scalar>
    class HERMES_API Solution : public MeshFunction<Scalar>
    {
    public:

      void init();

      Solution();

      Solution(Mesh *mesh);

      Solution(Mesh *mesh, Scalar init_const);

      Solution(Mesh *mesh, Scalar init_const_0, Scalar init_const_1);

      Solution (Space<Scalar>* s, Vector<Scalar>* coeff_vec);

      Solution (Space<Scalar>* s, Scalar* coeff_vec);

      virtual ~Solution();

      virtual void free();

      void assign(Solution<Scalar>* sln);

      Solution& operator = (Solution& sln) { assign(&sln); return *this; }

      void copy(const Solution<Scalar>* sln);

      int* get_element_orders() { return this->elem_orders;}

      void set_const(Mesh* mesh, Scalar c);

      void set_const(Mesh* mesh, Scalar c0, Scalar c1); ///< two-component (Hcurl) const

      void set_zero(Mesh* mesh);

      void set_zero_2(Mesh* mesh); ///< two-component (Hcurl) zero

      int get_edge_fn_order(int edge);

      int get_edge_fn_order(int edge, Space<Scalar>* space, Element* e = NULL);

      /// Enables or disables transformation of the solution derivatives (H1 case)
      /// or values (vector (Hcurl) case). This means H2D_FN_DX_0 and H2D_FN_DY_0 or
      /// H2D_FN_VAL_0 and H2D_FN_VAL_1 will or will not be returned premultiplied by the reference
      /// mapping matrix. The default is enabled (true).
      void enable_transform(bool enable = true);

      /// Saves the complete solution (i.e., including the internal copy of the mesh and
      /// element orders) to a binary file. On Linux, if `compress` is true, the file is
      /// compressed with gzip and a ".gz" suffix added to the file name.
      void save(const char* filename, bool compress = true);

      /// Loads the solution from a file previously created by Solution::save(). This completely
      /// restores the solution in the memory. The file name has to include the ".gz" suffix,
      /// in which case the file is piped through gzip to decompress the data (Linux only).
      void load(const char* filename);

      /// Returns solution value or derivatives at element e, in its reference domain point (xi1, xi2).
      /// 'item' controls the returned value: 0 = value, 1 = dx, 2 = dy, 3 = dxx, 4 = dyy, 5 = dxy.
      /// NOTE: This function should be used for postprocessing only, it is not effective
      /// enough for calculations.
      Scalar get_ref_value(Element* e, double xi1, double xi2, int component = 0, int item = 0);

      /// Returns solution value or derivatives (correctly transformed) at element e, in its reference
      /// domain point (xi1, xi2). 'item' controls the returned value: 0 = value, 1 = dx, 2 = dy,
      /// 3 = dxx, 4 = dyy, 5 = dxy.
      /// NOTE: This function should be used for postprocessing only, it is not effective
      /// enough for calculations.
      Scalar get_ref_value_transformed(Element* e, double xi1, double xi2, int a, int b);

      /// Returns solution value or derivatives at the physical domain point (x, y).
      /// 'item' controls the returned value: H2D_FN_VAL_0, H2D_FN_VAL_1, H2D_FN_DX_0, H2D_FN_DX_1, H2D_FN_DY_0,....
      /// NOTE: This function should be used for postprocessing only, it is not effective
      /// enough for calculations. Since it searches for an element sequentinally, it is extremelly
      /// slow. Prefer Solution::get_ref_value if possible.
      virtual Scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0);

      /// Returns the number of degrees of freedom of the solution.
      /// Returns -1 for exact or constant solutions.
      int get_num_dofs() const { return num_dofs; };

      /// Returns solution type.
      SolutionType get_type() const { return sln_type; };

      /// Returns space type.
      SpaceType get_space_type() const { return space_type; };

    public:
      /// Internal.
      virtual void set_active_element(Element* e);

      /// Passes solution components calculated from solution vector as Solutions.
      static void vector_to_solutions(Scalar* solution_vector, Hermes::vector<Space<Scalar>*> spaces,
        Hermes::vector<Solution<Scalar>*> solutions,
        Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>());

      static void vector_to_solution(Scalar* solution_vector, Space<Scalar>* space, Solution<Scalar>* solution,
        bool add_dir_lift = true);

      static void vector_to_solutions(Vector<Scalar>* vec, Hermes::vector<Space<Scalar>*> spaces,
        Hermes::vector<Solution<Scalar>*> solutions,
        Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>());

      static void vector_to_solution(Vector<Scalar>* vec, Space<Scalar>* space, Solution<Scalar>* solution,
        bool add_dir_lift = true);

      static void vector_to_solutions(Scalar* solution_vector, Hermes::vector<Space<Scalar>*> spaces,
        Hermes::vector<Solution<Scalar>*> solutions, Hermes::vector<PrecalcShapeset *> pss,
        Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>());

      static void vector_to_solution(Scalar* solution_vector, Space<Scalar>* space, Solution<Scalar>* solution,
        PrecalcShapeset* pss, bool add_dir_lift = true);

      /// If this is set to true, the mesh was created by this instance of this class.
      bool own_mesh;
    protected:

      static Quad2DCheb g_quad_2d_cheb;

      /// Converts a coefficient vector into a Solution.
      virtual void set_coeff_vector(Space<Scalar>* space, Vector<Scalar>* vec, bool add_dir_lift);

      virtual void set_coeff_vector(Space<Scalar>* space, PrecalcShapeset* pss, Scalar* coeffs, bool add_dir_lift);

      virtual void set_coeff_vector(Space<Scalar>* space, Scalar* coeffs, bool add_dir_lift);

      SolutionType sln_type;

      bool transform;

      /// Precalculated tables for last four used elements.
      /// There is a 2-layer structure of the precalculated tables.
      /// The first (the lowest) one is the layer where mapping of integral orders to
      /// Function::Node takes place. See function.h for details.
      /// The second one is the layer with mapping of sub-element transformation to
      /// a table from the lowest layer.
      /// The highest layer (in contrast to the PrecalcShapeset class) is represented
      /// here only by this array.
      std::map<uint64_t, LightArray<struct Function<Scalar>::Node*>*>* tables[4][4];

      Element* elems[4][4];

      int cur_elem, oldest[4];

      Scalar* mono_coefs;  ///< monomial coefficient array

      int* elem_coefs[2];  ///< array of pointers into mono_coefs

      int* elem_orders;    ///< stored element orders

      int num_coefs, num_elems;

      int num_dofs;

      SpaceType space_type;

      void transform_values(int order, struct Function<Scalar>::Node* node, int newmask, int oldmask, int np);

      Scalar   cnst[2];

      virtual void precalculate(int order, int mask);

      Scalar* dxdy_coefs[2][6];

      Scalar* dxdy_buffer;

      double** calc_mono_matrix(int o, int*& perm);

      void init_dxdy_buffer();

      void free_tables();

      Element* e_last; ///< last visited element when getting solution values at specific points

    };
  }
}
#endif
