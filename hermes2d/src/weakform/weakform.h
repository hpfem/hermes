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
template<typename Scalar> class DiscreteProblem;
template<typename Scalar> class Space;
template<typename Scalar> class MeshFunction;
struct SurfPos;
class Ord;
template<typename Scalar> class Precond;

struct Element;
class Shapeset;
template<typename T> class Func;
template<typename T> class Geom;
template<typename T> class ExtData;

template<typename Scalar> class Stage;
template<typename Scalar> class Form;
template<typename Scalar> class MatrixFormVol;
template<typename Scalar> class VectorFormVol;
template<typename Scalar> class MatrixFormSurf;
template<typename Scalar> class VectorFormSurf;

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


template<typename Scalar>
class HERMES_API WeakForm
{
public:

  WeakForm(unsigned int neq = 1, bool mat_free = false);
  ~WeakForm();

  Mesh::ElementMarkersConversion* get_element_markers_conversion() { 
    return element_markers_conversion; 
  };
  Mesh::BoundaryMarkersConversion* get_boundary_markers_conversion() { 
    return boundary_markers_conversion; 
  };

  // General case.
  void add_matrix_form(MatrixFormVol<Scalar>* mfv);
  void add_matrix_form_surf(MatrixFormSurf<Scalar>* mfs);
  void add_vector_form(VectorFormVol<Scalar>* vfv);
  void add_vector_form_surf(VectorFormSurf<Scalar>* vfs);

  void set_ext_fns(void* fn, Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>());

  /// Returns the number of equations.
  unsigned int get_neq() { return neq; }

  /// Returns volumetric and surface weak forms.
  Hermes::vector<MatrixFormVol<Scalar> *> get_mfvol() { return mfvol; }
  Hermes::vector<MatrixFormSurf<Scalar> *> get_mfsurf() { return mfsurf; }
  Hermes::vector<VectorFormVol<Scalar> *> get_vfvol() { return vfvol; }
  Hermes::vector<VectorFormSurf<Scalar> *> get_vfsurf() { return vfsurf; }

  /// Sets volumetric and surface weak forms.
  void set_mfvol(Hermes::vector<MatrixFormVol<Scalar> *> mfvol) { this->mfvol = mfvol; }
  void set_mfsurf(Hermes::vector<MatrixFormSurf<Scalar> *> mfvol) { this->mfsurf = mfsurf; }
  void set_vfvol(Hermes::vector<VectorFormVol<Scalar> *> vfvol) { this->vfvol = vfvol; }
  void set_vfsurf(Hermes::vector<VectorFormSurf<Scalar> *> vfvol) { this->vfsurf = vfsurf; }

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

  /// For time-dependent right-hand side functions.
  void set_current_time(double time);
  double get_current_time();

protected:
  double current_time;
  unsigned int neq;
  int seq;
  bool is_matfree;

  struct Area  { Hermes::vector<std::string> markers;  };

  Hermes::vector<Area> areas;

public:
  // General case.
  Hermes::vector<MatrixFormVol<Scalar> *> mfvol;
  Hermes::vector<MatrixFormSurf<Scalar> *> mfsurf;
  Hermes::vector<VectorFormVol<Scalar> *> vfvol;
  Hermes::vector<VectorFormSurf<Scalar> *> vfsurf;

  typename void get_stages(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*>& u_ext,
                  std::vector<Stage<Scalar>>& stages, bool rhsonly);
  bool** get_blocks(bool force_diagonal_blocks);

  bool is_in_area(std::string marker, std::string area) const
  { return area == marker; }

  bool is_sym() const { return false; /* not impl. yet */ }

  friend class DiscreteProblem<Scalar>;
  friend class Precond<Scalar>;

  // To be called only by the constructor of DiscreteProblem.
  void set_markers_conversion(Mesh::ElementMarkersConversion* element_markers_conversion, 
                              Mesh::BoundaryMarkersConversion* boundary_markers_conversion)
  {
    this->element_markers_conversion = element_markers_conversion;
    this->boundary_markers_conversion = boundary_markers_conversion;
  }

protected:

  typename Stage<Scalar>* find_stage(std::vector<Stage<Scalar>>& stages, int ii, int jj,
                    Mesh* m1, Mesh* m2,
                    Hermes::vector<MeshFunction<Scalar>*>& ext, Hermes::vector<Solution<Scalar>*>& u_ext);

  Mesh::ElementMarkersConversion* element_markers_conversion;
  Mesh::BoundaryMarkersConversion* boundary_markers_conversion;
};

template<typename Scalar>
class HERMES_API Form
{
public:
  Form(std::string area = HERMES_ANY, Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

  inline void set_weakform(WeakForm<Scalar>* wf) { this->wf = wf; }

  std::string area;
  Hermes::vector<MeshFunction<Scalar>*> ext;
  Hermes::vector<Scalar> param;
  // Form will be always multiplied (scaled) with this number.
  double scaling_factor;
  // External solutions for this form will start
  // with u_ext[u_ext_offset] where u_ext[] are external
  // solutions coming to the assembling procedure via the
  // external coefficient vector.
  int u_ext_offset;
  // If true, the form will be evaluated using adaptive
  // numerical integration.
  bool adapt_eval;
  // To obtain reference value, the element is split into
  // four sons. In addition, the order is increased by this value.
  int adapt_order_increase;
  // Max. allowed relative error (stopping criterion for adaptive
  // numerical quadrature.
  double adapt_rel_error_tol;

protected:
  WeakForm<Scalar>* wf;
};

template<typename Scalar>
class HERMES_API MatrixFormVol : public Form<Scalar>
{
public:
  MatrixFormVol(unsigned int i, unsigned int j, SymFlag sym = HERMES_NONSYM, 
                std::string area = HERMES_ANY, 
                Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
                Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
                double scaling_factor = 1.0, int u_ext_offset = 0);

  virtual MatrixFormVol* clone();

  unsigned int i, j;
  int sym;

  virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
                        Geom<double> *e, ExtData<Scalar> *ext);
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext);
};

template<typename Scalar>
class HERMES_API MatrixFormSurf : public Form<Scalar>
{
public:
  MatrixFormSurf(unsigned int i, unsigned int j, std::string area = HERMES_ANY, 
                  Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
                  Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
                  double scaling_factor = 1.0, int u_ext_offset = 0);

  virtual MatrixFormSurf* clone();

  unsigned int i, j;

  virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
                        Geom<double> *e, ExtData<Scalar> *ext);
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext);
};

template<typename Scalar>
class HERMES_API VectorFormVol : public Form<Scalar>
{
public:
  VectorFormVol(unsigned int i, std::string area = HERMES_ANY, 
                Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
                Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
                double scaling_factor = 1.0, int u_ext_offset = 0);

  virtual VectorFormVol* clone();

  unsigned int i;

  virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, 
                        Geom<double> *e, ExtData<Scalar> *ext);
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                  ExtData<Ord> *ext);
};

template<typename Scalar>
class HERMES_API VectorFormSurf : public Form<Scalar>
{
public:
  VectorFormSurf(unsigned int i, std::string area = HERMES_ANY, 
                  Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
                  Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
                  double scaling_factor = 1.0, int u_ext_offset = 0);
    
  virtual VectorFormSurf* clone();

  unsigned int i;

  virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, 
                        Geom<double> *e, ExtData<Scalar> *ext);
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext);
};

template<typename Scalar>
class Stage
{
  Hermes::vector<int> idx;
  Hermes::vector<Mesh*> meshes;
  Hermes::vector<Transformable*> fns;
  Hermes::vector<MeshFunction<Scalar>*> ext;

  // general case
  Hermes::vector<MatrixFormVol<Scalar> *> mfvol;
  Hermes::vector<MatrixFormSurf<Scalar> *> mfsurf;
  Hermes::vector<VectorFormVol<Scalar> *> vfvol;
  Hermes::vector<VectorFormSurf<Scalar> *> vfsurf;

  std::set<int> idx_set;
  std::set<unsigned> seq_set;
  std::set<MeshFunction<Scalar>*> ext_set;
};


#endif
