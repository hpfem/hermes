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


#ifndef __H2D_WEAKFORM_H
#define __H2D_WEAKFORM_H

#include "function.h"
#include "solution.h"

class RefMap;
class DiscreteProblem;
class Space;
class MeshFunction;
struct SurfPos;
class Ord;

struct Element;
class Shapeset;
template<typename T> class Func;
template<typename T> class Geom;
template<typename T> class ExtData;

// Bilinear form symmetry flag, see WeakForm::add_matrix_form
enum SymFlag
{
  HERMES_ANTISYM = -1,
  HERMES_UNSYM = 0,
  HERMES_SYM = 1
};

/// \brief Represents the weak formulation of a problem.
///
/// The WeakForm class represents the weak formulation of a system of linear PDEs.
/// The number of equations ("neq") in the system is fixed and is passed to the constructor.
/// The weak formulation of the system A(U,V) = L(V) has a block structure. A(U,V) is
/// a (neq x neq) matrix of bilinear forms a_mn(u,v) and L(V) is a neq-component vector
/// of linear forms l(v). U and V are the vectors of basis and test functions.
///
///
///

class HERMES_API WeakForm
{
public:

  WeakForm(int neq = 1, bool mat_free = false);

  // general case
  typedef scalar (*matrix_form_val_t)(int n, double *wt, Func<scalar> *u[], Func<double> *vi, Func<double> *vj, Geom<double> *e, ExtData<scalar> *);
  typedef Ord (*matrix_form_ord_t)(int n, double *wt, Func<Ord> *u[], Func<Ord> *vi, Func<Ord> *vj, Geom<Ord> *e, ExtData<Ord> *);
  typedef scalar (*vector_form_val_t)(int n, double *wt, Func<scalar> *u[], Func<double> *vi, Geom<double> *e, ExtData<scalar> *);
  typedef Ord (*vector_form_ord_t)(int n, double *wt, Func<Ord> *u[], Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *);

  // general case
  void add_matrix_form(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, 
		   SymFlag sym = HERMES_UNSYM, int area = HERMES_ANY, Hermes::Tuple<MeshFunction*>ext = Hermes::Tuple<MeshFunction*>());
  void add_matrix_form(matrix_form_val_t fn, matrix_form_ord_t ord, 
		   SymFlag sym = HERMES_UNSYM, int area = HERMES_ANY, Hermes::Tuple<MeshFunction*>ext = Hermes::Tuple<MeshFunction*>()); // single equation case
  void add_matrix_form_surf(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, 
			int area = HERMES_ANY, Hermes::Tuple<MeshFunction*>ext = Hermes::Tuple<MeshFunction*>());
  void add_matrix_form_surf(matrix_form_val_t fn, matrix_form_ord_t ord, 
			int area = HERMES_ANY, Hermes::Tuple<MeshFunction*>ext = Hermes::Tuple<MeshFunction*>()); // single equation case
  void add_vector_form(int i, vector_form_val_t fn, vector_form_ord_t ord, 
		   int area = HERMES_ANY, Hermes::Tuple<MeshFunction*>ext = Hermes::Tuple<MeshFunction*>());
  void add_vector_form(vector_form_val_t fn, vector_form_ord_t ord, 
		   int area = HERMES_ANY, Hermes::Tuple<MeshFunction*>ext = Hermes::Tuple<MeshFunction*>()); // single equation case
  void add_vector_form_surf(int i, vector_form_val_t fn, vector_form_ord_t ord, 
			int area = HERMES_ANY, Hermes::Tuple<MeshFunction*>ext = Hermes::Tuple<MeshFunction*>());
  void add_vector_form_surf(vector_form_val_t fn, vector_form_ord_t ord, 
			int area = HERMES_ANY, Hermes::Tuple<MeshFunction*>ext = Hermes::Tuple<MeshFunction*>()); // single equation case

  void set_ext_fns(void* fn, Hermes::Tuple<MeshFunction*>ext = Hermes::Tuple<MeshFunction*>());

  /// Returns the number of equations
  int get_neq() { return neq; }

  /// Internal. Used by DiscreteProblem to detect changes in the weakform.
  int get_seq() const { return seq; }

  bool is_matrix_free() { return is_matfree; }

protected:
  int neq;
  int seq;
  bool is_matfree;

  struct Area  {  /*std::string name;*/  std::vector<int> markers;  };

  HERMES_API_USED_STL_VECTOR(Area);
  std::vector<Area> areas;
  HERMES_API_USED_STL_VECTOR(MeshFunction*);

  public:
    scalar evaluate_fn(int point_cnt, double *weights, Func<double> *values_v, Geom<double> *geometry, ExtData<scalar> *values_ext_fnc, Element* element, Shapeset* shape_set, int shape_inx); ///< Evaluate value of the user defined function.
    Ord evaluate_ord(int point_cnt, double *weights, Func<Ord> *values_v, Geom<Ord> *geometry, ExtData<Ord> *values_ext_fnc, Element* element, Shapeset* shape_set, int shape_inx); ///< Evaluate order of the user defined function.

  // general case
  struct MatrixFormVol  {  int i, j, sym, area;  matrix_form_val_t fn;  matrix_form_ord_t ord;  std::vector<MeshFunction *> ext; };
  struct MatrixFormSurf {  int i, j, area;       matrix_form_val_t fn;  matrix_form_ord_t ord;  std::vector<MeshFunction *> ext; };
  struct VectorFormVol  {  int i, area;          vector_form_val_t fn;  vector_form_ord_t ord;  std::vector<MeshFunction *> ext; };
  struct VectorFormSurf {  int i, area;          vector_form_val_t fn;  vector_form_ord_t ord;  std::vector<MeshFunction *> ext; };

  // general case
  std::vector<MatrixFormVol>  mfvol;
  std::vector<MatrixFormSurf> mfsurf;
  std::vector<VectorFormVol>  vfvol;
  std::vector<VectorFormSurf> vfsurf;

  struct Stage
  {
    std::vector<int> idx;
    std::vector<Mesh*> meshes;
    std::vector<Transformable*> fns;
    std::vector<MeshFunction*> ext;

    // general case
    std::vector<MatrixFormVol *>  mfvol;
    std::vector<MatrixFormSurf *> mfsurf;
    std::vector<VectorFormVol *>  vfvol;
    std::vector<VectorFormSurf *> vfsurf;

    std::set<int> idx_set;
    std::set<unsigned> seq_set;
    std::set<MeshFunction*> ext_set;
  };

  void get_stages(Hermes::Tuple< Space* > spaces, Hermes::Tuple< Solution* >& u_ext, 
                  std::vector< WeakForm::Stage >& stages, bool rhsonly);
  bool** get_blocks();

  bool is_in_area(int marker, int area) const
    { return area >= 0 ? area == marker : is_in_area_2(marker, area); }

  bool is_sym() const { return false; /* not impl. yet */ }

//  friend class DiscreteProblem;
//  friend class RefDiscreteProblem;
  friend class LinearProblem;
  friend class DiscreteProblem;
  friend class Precond;


private:

  Stage* find_stage(std::vector<WeakForm::Stage>& stages, int ii, int jj,
                    Mesh* m1, Mesh* m2, 
                    std::vector<MeshFunction*>& ext, Hermes::Tuple<Solution*>& u_ext);

  bool is_in_area_2(int marker, int area) const;
};

#endif
