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

#include "../function/function.h"
#include "../function/solution.h"
#include <string>

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
  HERMES_NONSYM = 0,
  HERMES_SYM = 1
};

/// \brief Represents the weak formulation of a PDE problem.
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

  // General case.
  typedef scalar (*matrix_form_val_t)(int n, double *wt, Func<scalar> *u[], Func<double> *vi, Func<double> *vj,
                  Geom<double> *e, ExtData<scalar> *);
  typedef Ord (*matrix_form_ord_t)(int n, double *wt, Func<Ord> *u[], Func<Ord> *vi, Func<Ord> *vj,
               Geom<Ord> *e, ExtData<Ord> *);
  typedef scalar (*vector_form_val_t)(int n, double *wt, Func<scalar> *u[], Func<double> *vi,
                  Geom<double> *e, ExtData<scalar> *);
  typedef Ord (*vector_form_ord_t)(int n, double *wt, Func<Ord> *u[], Func<Ord> *vi,
                                   Geom<Ord> *e, ExtData<Ord> *);

// Matrix forms for error calculation.
  typedef scalar (*error_matrix_form_val_t) (int n, double *wt, Func<scalar> *u_ext[],
                                             Func<scalar> *u, Func<scalar> *v, Geom<double> *e,
                                             ExtData<scalar> *); ///< Error bilinear form callback function.
  typedef Ord (*error_matrix_form_ord_t) (int n, double *wt, Func<Ord> *u_ext[],
                                          Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                                          ExtData<Ord> *); ///< Error bilinear form to estimate order of a function.

// Vector forms for error calculation.
  typedef scalar (*error_vector_form_val_t) (int n, double *wt, Func<scalar> *u_ext[],
                                             Func<scalar> *u, Geom<double> *e,
                                             ExtData<scalar> *); ///< Error linear form callback function.
  typedef Ord (*error_vector_form_ord_t) (int n, double *wt, Func<Ord> *u_ext[],
                                          Func<Ord> *u, Geom<Ord> *e,
                                          ExtData<Ord> *); ///< Error linear form to estimate order of a function.

  // General case.
  struct MatrixFormVol  {
    int i, j, sym, area;
    matrix_form_val_t fn;
    matrix_form_ord_t ord;
    Hermes::vector<MeshFunction *> ext;
    double scaling_factor;     // Form will be always multiplied (scaled) with this number.
    int u_ext_offset;          // External solutions for this form will start 
                               // with u_ext[u_ext_offset] where u_ext[] are external
                               // solutions coming to the assembling procedure via the 
                               // external coefficient vector.
  };
  struct MatrixFormSurf {
    int i, j, area;
    matrix_form_val_t fn;
    matrix_form_ord_t ord;
    Hermes::vector<MeshFunction *> ext;
    double scaling_factor;     // Form will be always multiplied (scaled) with this number.
    int u_ext_offset;          // External solutions for this form will start 
                               // with u_ext[u_ext_offset] where u_ext[] are external
                               // solutions coming to the assembling procedure via the 
                               // external coefficient vector.
  };
  struct VectorFormVol  {
    int i, area;
    vector_form_val_t fn;
    vector_form_ord_t ord;
    Hermes::vector<MeshFunction *> ext;
    double scaling_factor;     // Form will be always multiplied (scaled) with this number.
    int u_ext_offset;          // External solutions for this form will start 
                               // with u_ext[u_ext_offset] where u_ext[] are external
                               // solutions coming to the assembling procedure via the 
                               // external coefficient vector.
  };
  struct VectorFormSurf {
    int i, area;
    vector_form_val_t fn;
    vector_form_ord_t ord;
    Hermes::vector<MeshFunction *> ext;
    double scaling_factor;     // Form will be always multiplied (scaled) with this number.
    int u_ext_offset;          // External solutions for this form will start 
                               // with u_ext[u_ext_offset] where u_ext[] are external
                               // solutions coming to the assembling procedure via the 
                               // external coefficient vector.
  };

  // General case.
  void add_matrix_form(MatrixFormVol* mfv);
  void add_matrix_form(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord,
		       SymFlag sym = HERMES_NONSYM, int area = HERMES_ANY,
                       Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());
  void add_matrix_form(matrix_form_val_t fn, matrix_form_ord_t ord,
		       SymFlag sym = HERMES_NONSYM, int area = HERMES_ANY,
                       Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>()); // single equation case
  void add_matrix_form_surf(MatrixFormSurf* mfs);
  void add_matrix_form_surf(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord,
			    int area = HERMES_ANY,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());
  void add_matrix_form_surf(matrix_form_val_t fn, matrix_form_ord_t ord,
			    int area = HERMES_ANY,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>()); // single equation case
  void add_vector_form(VectorFormVol* vfv);
  void add_vector_form(int i, vector_form_val_t fn, vector_form_ord_t ord,
		       int area = HERMES_ANY,
                       Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());
  void add_vector_form(vector_form_val_t fn, vector_form_ord_t ord,
		       int area = HERMES_ANY,
                       Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>()); // single equation case
  void add_vector_form_surf(VectorFormSurf* vfs);
  void add_vector_form_surf(int i, vector_form_val_t fn, vector_form_ord_t ord,
			    int area = HERMES_ANY,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());
  void add_vector_form_surf(vector_form_val_t fn, vector_form_ord_t ord,
			    int area = HERMES_ANY,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>()); // single equation case

  // Wrapper functions utilizing the MarkersConversion class.
  void add_matrix_form(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord,
                       SymFlag sym, std::string area,
                       Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());
  void add_matrix_form(matrix_form_val_t fn, matrix_form_ord_t ord,
		       SymFlag sym, std::string area,
                       Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>()); // single equation case
  void add_matrix_form_surf(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord,
			    std::string area,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());
  void add_matrix_form_surf(matrix_form_val_t fn, matrix_form_ord_t ord,
                            std::string area,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>()); // single equation case
  void add_vector_form(int i, vector_form_val_t fn, vector_form_ord_t ord,
                       std::string area,
                       Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());
  void add_vector_form(vector_form_val_t fn, vector_form_ord_t ord,
                       std::string area,
                       Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>()); // single equation case
  void add_vector_form_surf(int i, vector_form_val_t fn, vector_form_ord_t ord,
                            std::string area,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());
  void add_vector_form_surf(vector_form_val_t fn, vector_form_ord_t ord,
			    std::string area,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>()); // single equation case

  void set_ext_fns(void* fn, Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());

  /// Returns the number of equations.
  int get_neq() { return neq; }

  /// Returns volumetric and surface weak forms.
  Hermes::vector<MatrixFormVol> get_mfvol() { return mfvol; }
  Hermes::vector<MatrixFormSurf> get_mfsurf() { return mfsurf; }
  Hermes::vector<VectorFormVol> get_vfvol() { return vfvol; }
  Hermes::vector<VectorFormSurf> get_vfsurf() { return vfsurf; }

  /// Sets volumetric and surface weak forms.
  void set_mfvol(Hermes::vector<MatrixFormVol> mfvol) { this->mfvol = mfvol; }
  void set_mfsurf(Hermes::vector<MatrixFormSurf> mfvol) { this->mfsurf = mfsurf; }
  void set_vfvol(Hermes::vector<VectorFormVol> vfvol) { this->vfvol = vfvol; }
  void set_vfsurf(Hermes::vector<VectorFormSurf> vfvol) { this->vfsurf = vfsurf; }

  /// Deletes all volumetric and surface forms.
  void delete_all()
  {
    mfvol.clear();
    mfsurf.clear();
    vfvol.clear();
    vfsurf.clear();
  };

  /// Internal. Used by DiscreteProblem to detect changes in the weakform.
  int get_seq() const { return seq; }

  bool is_matrix_free() { return is_matfree; }

protected:
  int neq;
  int seq;
  bool is_matfree;

  struct Area  {  /*std::string name;*/  Hermes::vector<int> markers;  };

  Hermes::vector<Area> areas;

public:
  scalar evaluate_fn(int point_cnt, double *weights, Func<double> *values_v,
                     Geom<double> *geometry, ExtData<scalar> *values_ext_fnc, Element* element,
                     Shapeset* shape_set, int shape_inx); ///< Evaluate value of the user defined function.
  Ord evaluate_ord(int point_cnt, double *weights, Func<Ord> *values_v,
                   Geom<Ord> *geometry, ExtData<Ord> *values_ext_fnc, Element* element,
                   Shapeset* shape_set, int shape_inx); ///< Evaluate order of the user defined function.

  // General case.
  Hermes::vector<MatrixFormVol>  mfvol;
  Hermes::vector<MatrixFormSurf> mfsurf;
  Hermes::vector<VectorFormVol>  vfvol;
  Hermes::vector<VectorFormSurf> vfsurf;

  // These members are used temporarily for storing markers defined by user-supplied strings.
  std::map<std::string, MatrixFormVol>  mfvol_string_temp;
  std::map<std::string, MatrixFormSurf> mfsurf_string_temp;
  std::map<std::string, VectorFormVol>  vfvol_string_temp;
  std::map<std::string, VectorFormSurf> vfsurf_string_temp;

  // Function which according to the conversion table provided, updates the above members.
  void update_markers_acc_to_conversion(Mesh::MarkersConversion* markers_conversion);

  struct Stage
  {
    Hermes::vector<int> idx;
    Hermes::vector<Mesh*> meshes;
    Hermes::vector<Transformable*> fns;
    Hermes::vector<MeshFunction*> ext;

    // general case
    Hermes::vector<MatrixFormVol *>  mfvol;
    Hermes::vector<MatrixFormSurf *> mfsurf;
    Hermes::vector<VectorFormVol *>  vfvol;
    Hermes::vector<VectorFormSurf *> vfsurf;

    std::set<int> idx_set;
    std::set<unsigned> seq_set;
    std::set<MeshFunction*> ext_set;
  };

  void get_stages(Hermes::vector< Space* > spaces, Hermes::vector< Solution* >& u_ext,
                  std::vector< WeakForm::Stage >& stages, bool rhsonly);
  bool** get_blocks(bool force_diagonal_blocks);

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
                    Hermes::vector<MeshFunction*>& ext, Hermes::vector<Solution*>& u_ext);

  bool is_in_area_2(int marker, int area) const;
};

#endif
