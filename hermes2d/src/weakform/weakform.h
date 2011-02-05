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

class BoundaryConditions;

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

  WeakForm(unsigned int neq = 1, bool mat_free = false);

  // General case.
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

  class HERMES_API Form
  {
  public:
    Form(int area = HERMES_ANY, Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);
<<<<<<< HEAD

=======
    
>>>>>>> Fixes & comments & old code deletion.
    inline void set_weakform(WeakForm* wf) { this->wf = wf; }

    int area;
    Hermes::vector<MeshFunction *> ext;
    // Form will be always multiplied (scaled) with this number.
    double scaling_factor;
    // External solutions for this form will start
    // with u_ext[u_ext_offset] where u_ext[] are external
    // solutions coming to the assembling procedure via the
    // external coefficient vector.
    int u_ext_offset;

  protected:
    WeakForm* wf;
  };

  class HERMES_API MatrixFormVol : public Form
  {
  public:
    MatrixFormVol(unsigned int i, unsigned int j, SymFlag sym = HERMES_NONSYM, int area = HERMES_ANY, Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

    unsigned int i, j;
    int sym;
<<<<<<< HEAD

=======
    
>>>>>>> Fixes & comments & old code deletion.
    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, ExtData<scalar> *ext);
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext);
  };

  class HERMES_API MatrixFormSurf : public Form
  {
  public:
    MatrixFormSurf(unsigned int i, unsigned int j, int area = HERMES_ANY, Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

    unsigned int i, j;
<<<<<<< HEAD

=======
    
>>>>>>> Fixes & comments & old code deletion.
    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext);
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext);
  };

  class HERMES_API VectorFormVol : public Form
  {
  public:
    VectorFormVol(unsigned int i, int area = HERMES_ANY, Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

    unsigned int i;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext);
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext);
  };

  class HERMES_API VectorFormSurf : public Form
  {
  public:
    VectorFormSurf(unsigned int i, int area = HERMES_ANY, Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

    unsigned int i;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext);
<<<<<<< HEAD
<<<<<<< HEAD
      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext);
=======
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext);
>>>>>>> Fixes & comments & old code deletion.
=======
      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext);
>>>>>>> Fixed forms and 07-general
  };

  // General case.
  void add_matrix_form(MatrixFormVol* mfv);
  void add_matrix_form_surf(MatrixFormSurf* mfs);
  void add_vector_form(VectorFormVol* vfv);
  void add_vector_form_surf(VectorFormSurf* vfs);

  void set_ext_fns(void* fn, Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());

  /// Returns the number of equations.
  unsigned int get_neq() { return neq; }

  /// Returns volumetric and surface weak forms.
  Hermes::vector<MatrixFormVol *> get_mfvol() { return mfvol; }
  Hermes::vector<MatrixFormSurf *> get_mfsurf() { return mfsurf; }
  Hermes::vector<VectorFormVol *> get_vfvol() { return vfvol; }
  Hermes::vector<VectorFormSurf *> get_vfsurf() { return vfsurf; }

  /// Sets volumetric and surface weak forms.
  void set_mfvol(Hermes::vector<MatrixFormVol *> mfvol) { this->mfvol = mfvol; }
  void set_mfsurf(Hermes::vector<MatrixFormSurf *> mfvol) { this->mfsurf = mfsurf; }
  void set_vfvol(Hermes::vector<VectorFormVol *> vfvol) { this->vfvol = vfvol; }
  void set_vfsurf(Hermes::vector<VectorFormSurf *> vfvol) { this->vfsurf = vfsurf; }

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
  unsigned int neq;
  int seq;
  bool is_matfree;
  BoundaryConditions* boundary_conditions;

  struct Area  {  /*std::string name;*/  Hermes::vector<int> markers;  };

  Hermes::vector<Area> areas;

<<<<<<< HEAD
public:
=======
public:  
>>>>>>> Added BoundaryConditions class to the Space class
  // General case.
  Hermes::vector<MatrixFormVol *> mfvol;
  Hermes::vector<MatrixFormSurf *> mfsurf;
  Hermes::vector<VectorFormVol *> vfvol;
  Hermes::vector<VectorFormSurf *> vfsurf;

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
    Hermes::vector<MatrixFormVol *> mfvol;
    Hermes::vector<MatrixFormSurf *> mfsurf;
    Hermes::vector<VectorFormVol *> vfvol;
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

  friend class DiscreteProblem;
  friend class Precond;

private:

  Stage* find_stage(std::vector<WeakForm::Stage>& stages, int ii, int jj,
                    Mesh* m1, Mesh* m2,
                    Hermes::vector<MeshFunction*>& ext, Hermes::vector<Solution*>& u_ext);

  bool is_in_area_2(int marker, int area) const;
};

#endif
