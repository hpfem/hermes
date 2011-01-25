// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include "function.h"
#include "shapefn.h"
#include "space/space.h"
#include "asmlist.h"
#include "refmap.h"

/// @defgroup solutions Solutions
///
/// TODO: description

typedef scalar (*exact_fn_t)(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz);

typedef scalar3 &(*exact_vec_fn_t)(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz);

/// Represents a function defined on a mesh.
///
/// MeshFunction is a base class for all classes representing an arbitrary function
/// superimposed on a mesh (ie., domain). These include the Solution, ExactSolution
/// and Filter classes, which define the concrete behavior and the way the function
/// is (pre)calculated. Any such function can later be visualized.
///
/// (This is an abstract class and cannot be instantiated.)
///
/// @ingroup solutions
class HERMES_API MeshFunction : public ScalarFunction {
public:
  MeshFunction(Mesh *mesh);
  virtual ~MeshFunction();

  virtual void set_active_element(Element *e);

  ESpaceType get_sptype() { return sptype; }

  Mesh *get_mesh() const { return mesh; }
  RefMap *get_refmap() { update_refmap(); return refmap; }

  /// @return Order of the function on the active element
  virtual Ord3 get_order() = 0;

protected:
  Mesh *mesh;
  RefMap *refmap;

  int mode, seq;
  bool noinc;

  ESpaceType sptype;

public:
  /// For internal use only.
  void update_refmap() { refmap->force_transform(sub_idx, ctm); }

  friend class DiscreteProblem;
  friend class LinearProblem;
};


/// Represents the solution of a PDE.
///
/// Solution represents the solution of a PDE. Given a space and a solution vector,
/// it calculates the appripriate linear combination of basis function at the specified
/// element and integration points.
///
/// @ingroup solutions
class HERMES_API Solution : public MeshFunction {
public:
  Solution(Mesh *mesh);
  virtual ~Solution();
  virtual void free();

  void assign(Solution *sln);
  Solution &operator=(Solution &sln) { assign(&sln); return *this; }
  void copy(const Solution *sln);

  void set_exact(exact_fn_t exactfn);
  void set_exact(exact_vec_fn_t exactfn);

  void set_const(scalar c);
  void set_const(scalar c0, scalar c1, scalar c2);		// three-component const

  void set_zero();
  void set_zero_3(); // three-component zero

  virtual void set_active_element(Element *e);

  /// Enables or disables transformation of the solution derivatives (H1 case)
  /// or values (vector (Hcurl) case). This means FN_DX_0 and FN_DY_0 or
  /// FN_VAL_0 and FN_VAL_1 will or will not be returned premultiplied by the reference
  /// mapping matrix. The default is enabled (true).
  void enable_transform(bool enable);

  /// Calculates the value of the solution at a point that is given in terms 
  /// of its physical mesh coordinates. 
  virtual scalar get_pt_value(double x, double y, double z, int comp = 0) {
    QuadPt3D pt(x, y, z, 1.0);
    precalculate(1, &pt, FN_VAL);
    return get_fn_values(comp)[0];
  }

  /// Calculates the values or derivatives of the solution at 'np'
  /// points. The points are provided in the array 'pt'. What exactly is 
  /// calculated is determined by the mask. In the scalar (H1) case, 
  /// possible values of mask are FN_VAL (for function values) and 
  /// FN_DX, FN_DY, FN_DZ (first partial derivatives). In the vector-valued
  /// Hcurl case, function value masks are FN_VAL_0, FN_VAL_1, FN_VAL_2
  /// for function values and FN_DX_0, FN_DX_1, FN_DX_2, FN_DY_0, ...
  /// for partial derivatives.
  virtual void precalculate(const int np, const QuadPt3D *pt, int mask);

  virtual Ord3 get_order();

  static void vector_to_solutions(scalar* solution_vector, Hermes::vector<Space *> spaces, Hermes::vector<Solution *> solutions, Hermes::vector<double> dir = Hermes::vector<double>());
  static void vector_to_solution(scalar* solution_vector, Space* space, Solution* solution, double dir = 1.0);
  
protected:
	
  virtual void set_coeff_vector(Space *space, scalar *vec, double dir = 1.0);
    static const int NUM_ELEMENTS = 4;

    enum { HERMES_UNDEF = -1, HERMES_SLN, HERMES_EXACT, HERMES_CONST } type;


    bool transform;
    bool own_mesh;

    scalar *mono_coefs;				/// monomial coefficient array
    int *elem_coefs[3];				/// array of pointers into mono_coefs
    Ord3 *elem_orders;			        /// stored element orders (copied from space)
    int num_coefs, num_elems;
    int num_dofs;

    scalar   cnst[3];				/// constant solution

    union {
      exact_fn_t exact_fn;			/// exact function
      exact_vec_fn_t exact_vec_fn;		/// exact vector function
    };

    scalar *dxdydz_coefs[3][4];			/// 3 components, 4 types of values (FN, DX, DY, DZ)
    scalar *dxdydz_buffer;

    void init_dxdydz_buffer();

    void precalculate_fe(const int np, const QuadPt3D *pt, int mask);
    void precalculate_exact(const int np, const QuadPt3D *pt, int mask);
    void precalculate_const(const int np, const QuadPt3D *pt, int mask);
};


/// Represents an exact solution to a PDE.
///
/// ExactSolution represents an arbitrary user-specified function defined on
/// a domain (mesh), typically an exact solution to a PDE. This can be used to
/// compare an approximate solution with an exact solution (see DiffFilter).
///
/// @ingroup solutions
class HERMES_API ExactSolution : public Solution {
public:
  ExactSolution(Mesh *mesh, exact_fn_t fn) : Solution(mesh) { set_exact(fn); }
  ExactSolution(Mesh *mesh, exact_vec_fn_t fn) : Solution(mesh) { set_exact(fn); }
};

#endif
