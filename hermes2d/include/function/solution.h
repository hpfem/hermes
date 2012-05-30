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

#include "../function/mesh_function.h"
#include "../space/space.h"
#include "../mesh/refmap.h"
#include "exceptions.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Quad2DCheb is a special "quadrature" consisting of product Chebyshev
    /// points on the reference triangle and quad. It is used for expressing
    /// the solution on an element as a linear combination of monomials.
    ///
    class Quad2DCheb;

    /// \brief Represents the solution of a PDE.
    ///
    /// The Solution class represents the solution of a PDE. Given a space and a solution vector,
    /// it calculates the appropriate linear combination of basis functions at the specified
    /// element and integration points.
    ///
    /// The higher-order solution on elements is best calculated not as a linear combination
    /// of shape functions (the usual approach), but as a linear combination of monomials.
    /// This has the advantage that no shape function table calculation and look-ups are
    /// necessary (except for the conversion of the solution coefficients), which means that
    /// visualization and multi-mesh calculations are much faster (all the push_transforms
    /// and table searches take the most time when evaluating the solution).
    ///
    /// The linear combination of monomials can be calculated using the Horner's scheme, which
    /// requires the same number of multiplications as the calculation of the linear combination
    /// of shape functions. However, sub-element transforms are trivial and cheap. Moreover,
    /// after the solution on all elements is expressed as a combination of monomials, the
    /// Space can be forgotten. This is comfortable for the user, since the Solution class acts
    /// as a self-contained unit, internally containing just a copy of the mesh and a table of
    /// monomial coefficients. It is also very straight-forward to obtain all derivatives of
    /// a solution defined in this way. Finally, it is possible to store the solution on the
    /// disk easily (no need to store the Space, which is difficult).
    ///
    /// The following is an example of the set of monomials for a cubic quad and a cubic triangle.
    /// (Note that these are actually the definitions of the polynomial spaces on these elements.)
    ///
    ///   [ x^3*y^3  x^2*y^3  x*y^3  y^3 ]       [                    y^3 ]
    ///   [ x^3*y^2  x^2*y^2  x*y^2  y^2 ]       [             x*y^2  y^2 ]
    ///   [ x^3*y    x^2*y    x*y    y   ]       [      x^2*y  x*y    y   ]
    ///   [ x^3      x^2      x      1   ]       [ x^3  x^2    x      1   ]
    ///
    /// The number of monomials is (p + 1)^2 for quads and (p + 1)*(p + 2)/2 for triangles, where
    /// 'p' is the polynomial degree.
    ///

    enum SolutionType {
      HERMES_UNDEF = -1,
      HERMES_SLN = 0,
      HERMES_EXACT = 1
    };

    template<typename Scalar>
    class HERMES_API Solution : public MeshFunction<Scalar>
    {
    public:
      Solution();
      Solution(Mesh *mesh);
      Solution (Space<Scalar>* s, Vector<Scalar>* coeff_vec);
      Solution (Space<Scalar>* s, Scalar* coeff_vec);
      virtual ~Solution();

      void assign(Solution* sln);
      inline Solution& operator = (Solution& sln) { assign(&sln); return *this; }

      virtual void copy(const Solution<Scalar>* sln);

      /// Sets solution equal to Dirichlet lift only, solution vector = 0.
      void set_dirichlet_lift(const Space<Scalar>* space, PrecalcShapeset* pss = NULL);

      /// Saves the complete solution (i.e., including the internal copy of the mesh and
      /// element orders) to an XML file.
      void save(const char* filename) const;

      /// Loads the solution from a file previously created by Solution::save(). This completely
      /// restores the solution in the memory.
      void load(const char* filename, Mesh* mesh);

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
      /// 'item' controls the returned value: H2D_FN_VAL_0, H2D_FN_VAL_1, H2D_FN_DX_0, H2D_FN_DX_1, H2D_FN_DY_0, ....
      /// NOTE: This function should be used for postprocessing only, it is not effective
      /// enough for calculations. Since it searches for an element sequentinally, it is extremelly
      /// slow. Prefer Solution::get_ref_value if possible.
      virtual Scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0);

      /// Multiplies the function represented by this class by the given coefficient.
      void multiply(Scalar coef);

      /// Returns solution type.
      inline SolutionType get_type() const { return sln_type; };

      /// Returns space type.
      inline SpaceType get_space_type() const { return space_type; };

      /// In case there is a space this solution belongs to, this returns the seq number of the space.
      /// This is used to check if the pointer returned by get_space() points to the same space, or if the space has changed.
      int get_space_seq();

      /// In case this is valid it returns a pointer to the space this solution belongs to.
      /// Only use when get_space() == get_space_seq();
      const Space<Scalar>* get_space();

      /// In case this is valid it returns a vector of coefficient wrt. to the basis of the finite dimensional space this solution belongs to.
      /// Only use when get_space() == get_space_seq();
      Scalar* get_sln_vector();

      /// Passes solution components calculated from solution vector as Solutions.
      static void vector_to_solutions(const Scalar* solution_vector, Hermes::vector<const Space<Scalar> *> spaces,
          Hermes::vector<Solution<Scalar>*> solutions, 
          Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>(),
          Hermes::vector<int> start_indices = Hermes::vector<int>());

      static void vector_to_solution(const Scalar* solution_vector, const Space<Scalar>* space, Solution<Scalar>* solution,
          bool add_dir_lift = true, int start_index = 0);

      static void vector_to_solutions(const Vector<Scalar>* vec, Hermes::vector<const Space<Scalar> *> spaces,
          Hermes::vector<Solution<Scalar>*> solutions, 
          Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>(),
          Hermes::vector<int> start_indices = Hermes::vector<int>());

      static void vector_to_solution(const Vector<Scalar>* vec, const Space<Scalar>* space, Solution<Scalar>* solution,
          bool add_dir_lift = true, int start_index = 0);

      static void vector_to_solutions(const Scalar* solution_vector, Hermes::vector<const Space<Scalar> *> spaces,
          Hermes::vector<Solution<Scalar>*> solutions, Hermes::vector<PrecalcShapeset *> pss, 
          Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>(),
          Hermes::vector<int> start_indices = Hermes::vector<int>());

      static void vector_to_solution(const Scalar* solution_vector, const Space<Scalar>* space, Solution<Scalar>* solution,
          PrecalcShapeset* pss, bool add_dir_lift = true, int start_index = 0);

      /// If this is set to true, the mesh was created by this instance of this class.
      bool own_mesh;

      /// Internal.
      virtual void set_active_element(Element* e);

      virtual MeshFunction<Scalar>* clone();
    protected:
      virtual int get_edge_fn_order(int edge) { return MeshFunction<Scalar>::get_edge_fn_order(edge); }

      /// Enables or disables transformation of the solution derivatives (H1 case)
      /// or values (vector (Hcurl) case). This means H2D_FN_DX_0 and H2D_FN_DY_0 or
      /// H2D_FN_VAL_0 and H2D_FN_VAL_1 will or will not be returned premultiplied by the reference
      /// mapping matrix. The default is enabled (true).
      void enable_transform(bool enable = true);

      virtual void init();

      virtual void free();

      /// In case this is valid it is a vector of coefficient wrt. to the basis of the finite dimensional space this solution belongs to.
      Scalar* sln_vector;

      /// Converts a coefficient vector into a Solution.
      virtual void set_coeff_vector(const Space<Scalar>* space, const Vector<Scalar>* vec, bool add_dir_lift, int start_index);

      virtual void set_coeff_vector(const Space<Scalar>* space, PrecalcShapeset* pss, const Scalar* coeffs, bool add_dir_lift, int start_index);

      virtual void set_coeff_vector(const Space<Scalar>* space, const Scalar* coeffs, bool add_dir_lift, int start_index);

      SolutionType sln_type;
      SpaceType space_type;

      /// In case this is valid it contains a pointer to the space this solution belongs to.
      const Space<Scalar>* space;

      /// In case there is a space this solution belongs to, this is the seq number of the space.
      int space_seq;

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

      Scalar* mono_coeffs;  ///< monomial coefficient array
      int* elem_coeffs[2];  ///< array of pointers into mono_coeffs
      /// Stored element orders in the mathematical sense. The polynomial degree of the highest basis function + increments due to the element shape, etc.  .
      int* elem_orders;
      int num_coeffs, num_elems;
      int num_dofs;

      void transform_values(int order, struct Function<Scalar>::Node* node, int newmask, int oldmask, int np);

      virtual void precalculate(int order, int mask);

      Scalar* dxdy_coeffs[2][6];

      Scalar* dxdy_buffer;

      double** calc_mono_matrix(int o, int*& perm);

      void init_dxdy_buffer();

      void free_tables();

      Element* e_last; ///< last visited element when getting solution values at specific points

      friend class RefMap;
      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class CalculationContinuity;
      template<typename T> friend class OGProjection;
      template<typename T> friend class OGProjectionNOX;
      template<typename T> friend class Adapt;
      template<typename T> friend class Func;
      template<typename T> friend class DiscontinuousFunc;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
      template<typename T> friend class NeighborSearch;
      template<typename T> friend HERMES_API Func<T>* init_fn(Solution<T>*fu, const int order);
      template<typename T> friend class RefinementSelectors::ProjBasedSelector;
      template<typename T> friend class RefinementSelectors::H1ProjBasedSelector;
      template<typename T> friend class RefinementSelectors::L2ProjBasedSelector;
      friend class RefinementSelectors::HcurlProjBasedSelector;
    };
  }
}

#endif
