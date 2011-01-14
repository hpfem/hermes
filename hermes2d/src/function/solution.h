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

#include "function.h"
#include "../space/space.h"
#include "../mesh/refmap.h"
#include "../../../hermes_common/matrix.h"

class PrecalcShapeset;


/// \brief Represents a function defined on a mesh.
///
/// MeshFunction is a base class for all classes representing an arbitrary function
/// superimposed on a mesh (ie., domain). These include the Solution, ExactSolution
/// and Filter classes, which define the concrete behavior and the way the function
/// is (pre)calculated. Any such function can later be visualized.
///
/// (This is an abstract class and cannot be instantiated.)
///
class HERMES_API MeshFunction : public ScalarFunction
{
public:

  MeshFunction();
  MeshFunction(Mesh *mesh);
  virtual ~MeshFunction();

  virtual void init() {};
  virtual void reinit() {free(); init();};

  virtual void set_quad_2d(Quad2D* quad_2d);
  virtual void set_active_element(Element* e);
  
  virtual int get_edge_fn_order(int edge) { return ScalarFunction::get_edge_fn_order(edge); }

  Mesh*   get_mesh() const { return mesh; }
  RefMap* get_refmap() { update_refmap(); return refmap; }

  virtual scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0) = 0;

  /// Virtual function handling overflows. Has to be virtual, because
  /// the necessary iterators in the templated class do not work with GCC.
  virtual void handle_overflow_idx();

  /// See Transformable::push_transform.
	virtual void push_transform(int son);

protected:

  int mode;
  Mesh* mesh;
  RefMap* refmap;

public:

  /// For internal use only.
  void force_transform(MeshFunction* mf)
    { ScalarFunction::force_transform(mf->get_transform(), mf->get_ctm()); }
  void update_refmap()
    { refmap->force_transform(sub_idx, ctm); }
  void force_transform(uint64_t sub_idx, Trf* ctm)
  {
    this->sub_idx = sub_idx;
    this->ctm = ctm;
  }

};

/// \brief Represents the solution of a PDE.
///
/// The Solution class represents the solution of a PDE. Given a space and a solution vector,
/// it calculates the appropriate linear combination of basis functions at the specified
/// element and integration points.
///
/// TODO: write how to obtain solution values, maybe include inherited methods from Function as comments.
///
class HERMES_API Solution : public MeshFunction
{
public:

  void init();
  Solution();
  Solution(Mesh *mesh);
  Solution(Mesh *mesh, ExactFunction exactfn);
  Solution(Mesh *mesh, scalar init_const); 
  Solution (Space* s, Vector* coeff_vec);
  Solution (Space* s, scalar* coeff_vec);
  virtual ~Solution();
  virtual void free();

  void assign(Solution* sln);
  Solution& operator = (Solution& sln) { assign(&sln); return *this; }
  void copy(const Solution* sln);

  int* get_element_orders() { return this->elem_orders;}

  void set_exact(Mesh* mesh, ExactFunction exactfn);
  void set_exact(Mesh* mesh, ExactFunction2 exactfn);

  void set_const(Mesh* mesh, scalar c);
  void set_const(Mesh* mesh, scalar c0, scalar c1); // two-component (Hcurl) const

  void set_zero(Mesh* mesh);
  void set_zero_2(Mesh* mesh); // two-component (Hcurl) zero

  virtual int get_edge_fn_order(int edge) { return MeshFunction::get_edge_fn_order(edge); }
  int get_edge_fn_order(int edge, Space* space, Element* e = NULL);
  
  /// Sets solution equal to Dirichlet lift only, solution vector = 0
  void set_dirichlet_lift(Space* space, PrecalcShapeset* pss = NULL);
  
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
  scalar get_ref_value(Element* e, double xi1, double xi2, int component = 0, int item = 0);

  /// Returns solution value or derivatives (correctly transformed) at element e, in its reference
  /// domain point (xi1, xi2). 'item' controls the returned value: 0 = value, 1 = dx, 2 = dy,
  /// 3 = dxx, 4 = dyy, 5 = dxy.
  /// NOTE: This function should be used for postprocessing only, it is not effective
  /// enough for calculations.
  scalar get_ref_value_transformed(Element* e, double xi1, double xi2, int a, int b);

  /// Returns solution value or derivatives at the physical domain point (x, y).
  /// 'item' controls the returned value: H2D_FN_VAL_0, H2D_FN_VAL_1, H2D_FN_DX_0, H2D_FN_DX_1, H2D_FN_DY_0,....
  /// NOTE: This function should be used for postprocessing only, it is not effective
  /// enough for calculations. Since it searches for an element sequentinally, it is extremelly
  /// slow. Prefer Solution::get_ref_value if possible.
  virtual scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0);

  /// Returns the number of degrees of freedom of the solution.
  /// Returns -1 for exact or constant solutions.
  int get_num_dofs() const { return num_dofs; };

  /// Multiplies the function represented by this class by the given coefficient.
  void multiply(scalar coef);

  /// Returns solution type.
  ESolutionType get_type() const { return sln_type; };

  /// Returns space type.
  ESpaceType get_space_type() const { return space_type; };



public:
  /// Internal.
  virtual void set_active_element(Element* e);

  /// Passes solution components calculated from solution vector as Solutions.
  static void vector_to_solutions(scalar* solution_vector, Hermes::vector<Space *> spaces, 
                                  Hermes::vector<Solution *> solutions, 
                                  Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>());
  static void vector_to_solution(scalar* solution_vector, Space* space, Solution* solution, 
                                 bool add_dir_lift = true);
  static void vector_to_solutions(Vector* vec, Hermes::vector<Space *> spaces, 
                                  Hermes::vector<Solution*> solutions, 
                                  Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>());
  static void vector_to_solution(Vector* vec, Space* space, Solution* solution, 
                                 bool add_dir_lift = true);
  static void vector_to_solutions(scalar* solution_vector, Hermes::vector<Space *> spaces, 
                                  Hermes::vector<Solution *> solutions, Hermes::vector<PrecalcShapeset *> pss, 
                                  Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>());
  static void vector_to_solution(scalar* solution_vector, Space* space, Solution* solution, 
                                 PrecalcShapeset* pss, bool add_dir_lift = true);

  bool own_mesh;
protected:

  /// Converts a coefficient vector into a Solution.
  virtual void set_coeff_vector(Space* space, Vector* vec, bool add_dir_lift);
  virtual void set_coeff_vector(Space* space, PrecalcShapeset* pss, scalar* coeffs, bool add_dir_lift);
  virtual void set_coeff_vector(Space* space, scalar* coeffs, bool add_dir_lift);

  ESolutionType sln_type;

  bool transform;

  /// Precalculated tables for last four used elements.
  /// There is a 2-layer structure of the precalculated tables.
  /// The first (the lowest) one is the layer where mapping of integral orders to 
  /// Function::Node takes place. See function.h for details.
  /// The second one is the layer with mapping of sub-element transformation to
  /// a table from the lowest layer.
  /// The highest layer (in contrast to the PrecalcShapeset class) is represented
  /// here only by this array.
  std::map<uint64_t, std::map<unsigned int, Node*>*>* tables[4][4];   

  Element* elems[4][4];
  int cur_elem, oldest[4];

  scalar* mono_coefs;  ///< monomial coefficient array
  int* elem_coefs[2];  ///< array of pointers into mono_coefs
  int* elem_orders;    ///< stored element orders
  int num_coefs, num_elems;
  int num_dofs;

  ESpaceType space_type;
  void transform_values(int order, Node* node, int newmask, int oldmask, int np);

  ExactFunction exactfn1;
  ExactFunction2 exactfn2;
  scalar   cnst[2];
  scalar   exact_mult;

  virtual void precalculate(int order, int mask);

  scalar* dxdy_coefs[2][6];
  scalar* dxdy_buffer;

  double** calc_mono_matrix(int o, int*& perm);
  void init_dxdy_buffer();
  void free_tables();

  Element* e_last; ///< last visited element when getting solution values at specific points

};


/// \brief Represents an exact solution of a PDE.
///
/// ExactSolution represents an arbitrary user-specified function defined on a domain (mesh),
/// typically an exact solution to a PDE. This can be used to compare an approximate solution
/// with an exact solution (see DiffFilter).
///
/// Please note that the same functionality can be obtained by using Solution::set_exact().
/// This class is provided merely for convenience.
///
class HERMES_API ExactSolution : public Solution
{
public:

  ExactSolution(Mesh* mesh, ExactFunction exactfn)
    { set_exact(mesh, exactfn); }

  ExactSolution(Mesh* mesh, ExactFunction2 exactfn)
    { set_exact(mesh, exactfn); }

  int update(Mesh* mesh, ExactFunction exactfn)
    { set_exact(mesh, exactfn);  return 1; }

  int update(Mesh* mesh, ExactFunction2 exactfn)
    { set_exact(mesh, exactfn);  return 1; }

};



#endif
