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

class Element;
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

// Geometrical type of weak forms.
enum GeomType
{
  HERMES_PLANAR = 0,         // Planar problem.
  HERMES_AXISYM_X = 1,       // Axisymmetric problem where x-axis is the axis of symmetry.
  HERMES_AXISYM_Y = 2        // Axisymmetric problem where y-axis is the axis of symmetry.
};

/// \brief Represents the weak formulation of a PDE problem.
///
/// The WeakForm class represents the weak formulation of a system of linear PDEs.
/// The number of equations ("neq") in the system is fixed and is passed to the constructor.
/// The weak formulation of the system A(U,V) = L(V) has a block structure. A(U,V) is
/// a (neq x neq) matrix of bilinear forms a_mn(u,v) and L(V) is a neq-component vector
/// of linear forms l(v). U and V are the vectors of basis and test functions.
///


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

  class HERMES_API Form
  {
  public:
    Form(std::string area = HERMES_ANY, Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
         Hermes::vector<scalar> param = Hermes::vector<scalar>(),
         double scaling_factor = 1.0, int u_ext_offset = 0);

    inline void set_weakform(WeakForm* wf) { this->wf = wf; }

    std::string area;
    Hermes::vector<MeshFunction *> ext;
    Hermes::vector<scalar> param;
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

    /// For time-dependent right-hand side functions.
    /// E.g. for Runge-Kutta methods. Otherwise the one time for the whole WeakForm can be used.
    void set_current_stage_time(double time);
    double get_current_stage_time() const;

  protected:
    WeakForm* wf;
    double stage_time;
  };

  class HERMES_API MatrixFormVol : public Form
  {
  public:
    MatrixFormVol(unsigned int i, unsigned int j, SymFlag sym = HERMES_NONSYM, 
                  std::string area = HERMES_ANY, 
                  Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
                  Hermes::vector<scalar> param = Hermes::vector<scalar>(),
                  double scaling_factor = 1.0, int u_ext_offset = 0);

    virtual MatrixFormVol* clone();

    unsigned int i, j;
    int sym;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
  };

  class HERMES_API MatrixFormSurf : public Form
  {
  public:
    MatrixFormSurf(unsigned int i, unsigned int j, std::string area = HERMES_ANY, 
                   Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
                   Hermes::vector<scalar> param = Hermes::vector<scalar>(),
                   double scaling_factor = 1.0, int u_ext_offset = 0);

    virtual MatrixFormSurf* clone();

    unsigned int i, j;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
  };

  class HERMES_API VectorFormVol : public Form
  {
  public:
    VectorFormVol(unsigned int i, std::string area = HERMES_ANY, 
                  Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
                  Hermes::vector<scalar> param = Hermes::vector<scalar>(),
                  double scaling_factor = 1.0, int u_ext_offset = 0);

    virtual VectorFormVol* clone();

    unsigned int i;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;
  };

  class HERMES_API VectorFormSurf : public Form
  {
  public:
    VectorFormSurf(unsigned int i, std::string area = HERMES_ANY, 
                   Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
                   Hermes::vector<scalar> param = Hermes::vector<scalar>(),
                   double scaling_factor = 1.0, int u_ext_offset = 0);
    
    virtual VectorFormSurf* clone();

    unsigned int i;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                      Geom<Ord> *e, ExtData<Ord> *ext) const;
  };

  class HERMES_API MultiComponentMatrixFormVol : public Form
  {
  public:
    MultiComponentMatrixFormVol(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, SymFlag sym = HERMES_NONSYM, 
                  std::string area = HERMES_ANY, 
                  Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
                  Hermes::vector<scalar> param = Hermes::vector<scalar>(),
                  double scaling_factor = 1.0, int u_ext_offset = 0);

    virtual MultiComponentMatrixFormVol* clone();

    Hermes::vector<std::pair<unsigned int, unsigned int> > coordinates;
    int sym;

    virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const = 0;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const = 0;
  };

  class HERMES_API MultiComponentMatrixFormSurf : public Form
  {
  public:
    MultiComponentMatrixFormSurf(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, std::string area = HERMES_ANY, 
                   Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
                   Hermes::vector<scalar> param = Hermes::vector<scalar>(),
                   double scaling_factor = 1.0, int u_ext_offset = 0);

    virtual MultiComponentMatrixFormSurf* clone();

    Hermes::vector<std::pair<unsigned int, unsigned int> > coordinates;

    virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const = 0;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const = 0;
  };

  class HERMES_API MultiComponentVectorFormVol : public Form
  {
  public:
    MultiComponentVectorFormVol(Hermes::vector<unsigned int> coordinates, std::string area = HERMES_ANY, 
                  Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
                  Hermes::vector<scalar> param = Hermes::vector<scalar>(),
                  double scaling_factor = 1.0, int u_ext_offset = 0);

    virtual MultiComponentVectorFormVol* clone();

    Hermes::vector<unsigned int> coordinates;

    virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const = 0;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const = 0;
  };

  class HERMES_API MultiComponentVectorFormSurf : public Form
  {
  public:
    MultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates, std::string area = HERMES_ANY, 
                   Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction*>(),
                   Hermes::vector<scalar> param = Hermes::vector<scalar>(),
                   double scaling_factor = 1.0, int u_ext_offset = 0);
    
    virtual MultiComponentVectorFormSurf* clone();

    Hermes::vector<unsigned int> coordinates;

    virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const = 0;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                      Geom<Ord> *e, ExtData<Ord> *ext) const = 0;
  };

  // General case.
  void add_matrix_form(MatrixFormVol* mfv);
  void add_matrix_form_surf(MatrixFormSurf* mfs);
  void add_vector_form(VectorFormVol* vfv);
  void add_vector_form_surf(VectorFormSurf* vfs);

  void add_multicomponent_matrix_form(MultiComponentMatrixFormVol* mfv);
  void add_multicomponent_matrix_form_surf(MultiComponentMatrixFormSurf* mfs);
  void add_multicomponent_vector_form(MultiComponentVectorFormVol* vfv);
  void add_multicomponent_vector_form_surf(MultiComponentVectorFormSurf* vfs);

  void set_ext_fns(void* fn, Hermes::vector<MeshFunction*> ext = Hermes::vector<MeshFunction*>());

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

  /// For time-dependent right-hand side functions.
  void set_current_time(double time);
  virtual double get_current_time() const;

protected:
  double current_time;
  unsigned int neq;
  int seq;
  bool is_matfree;

  struct Area  { Hermes::vector<std::string> markers;  };

  Hermes::vector<Area> areas;

public:
  // General case.
  Hermes::vector<MatrixFormVol *> mfvol;
  Hermes::vector<MatrixFormSurf *> mfsurf;
  Hermes::vector<VectorFormVol *> vfvol;
  Hermes::vector<VectorFormSurf *> vfsurf;

  Hermes::vector<MultiComponentMatrixFormVol *> mfvol_mc;
  Hermes::vector<MultiComponentMatrixFormSurf *> mfsurf_mc;
  Hermes::vector<MultiComponentVectorFormVol *> vfvol_mc;
  Hermes::vector<MultiComponentVectorFormSurf *> vfsurf_mc;

  // Storage of forms according to user-supplied strings.
  std::map<std::string, MatrixFormVol>  mfvol_string_temp;
  std::map<std::string, MatrixFormSurf> mfsurf_string_temp;
  std::map<std::string, VectorFormVol>  vfvol_string_temp;
  std::map<std::string, VectorFormSurf> vfsurf_string_temp;


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

    Hermes::vector<MultiComponentMatrixFormVol *> mfvol_mc;
    Hermes::vector<MultiComponentMatrixFormSurf *> mfsurf_mc;
    Hermes::vector<MultiComponentVectorFormVol *> vfvol_mc;
    Hermes::vector<MultiComponentVectorFormSurf *> vfsurf_mc;

    std::set<int> idx_set;
    std::set<unsigned> seq_set;
    std::set<MeshFunction*> ext_set;
  };

  void get_stages(Hermes::vector< Space* > spaces, Hermes::vector< Solution* >& u_ext,
                  std::vector< WeakForm::Stage >& stages, bool want_matrix, bool want_vector);
  bool** get_blocks(bool force_diagonal_blocks);

  bool is_in_area(std::string marker, std::string area) const
  { return area == marker; }

  bool is_sym() const { return false; /* not impl. yet */ }

  friend class DiscreteProblem;
  friend class Precond;

  // To be called only by the constructor of DiscreteProblem.
  void set_markers_conversion(Mesh::ElementMarkersConversion* element_markers_conversion, 
                              Mesh::BoundaryMarkersConversion* boundary_markers_conversion)
  {
    this->element_markers_conversion = element_markers_conversion;
    this->boundary_markers_conversion = boundary_markers_conversion;
  }

protected:

  Stage* find_stage(std::vector<WeakForm::Stage>& stages, int ii, int jj,
                    Mesh* m1, Mesh* m2,
                    Hermes::vector<MeshFunction*>& ext, Hermes::vector<Solution*>& u_ext);

  Stage* find_stage(std::vector<WeakForm::Stage>& stages, Hermes::vector<std::pair<unsigned int, unsigned int> > coordinates,
                    Mesh* m1, Mesh* m2,
                    Hermes::vector<MeshFunction*>& ext, Hermes::vector<Solution*>& u_ext);

  Stage* find_stage(std::vector<WeakForm::Stage>& stages, Hermes::vector<unsigned int> coordinates,
                    Mesh* m1, Mesh* m2,
                    Hermes::vector<MeshFunction*>& ext, Hermes::vector<Solution*>& u_ext);

  Mesh::ElementMarkersConversion* element_markers_conversion;
  Mesh::BoundaryMarkersConversion* boundary_markers_conversion;
};

#endif
