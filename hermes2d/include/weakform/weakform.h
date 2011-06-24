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

namespace Hermes
{
  namespace Hermes2D
  {
    class RefMap;
    template<typename Scalar> class DiscreteProblem;
    template<typename Scalar> class Space;
    template<typename Scalar> class MeshFunction;
    struct SurfPos;
    class Hermes::Ord;

    class Element;
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
    template<typename Scalar> class MultiComponentMatrixFormVol;
    template<typename Scalar> class MultiComponentVectorFormVol;
    template<typename Scalar> class MultiComponentMatrixFormSurf;
    template<typename Scalar> class MultiComponentVectorFormSurf;


    /// Bilinear form symmetry flag, see WeakForm::add_matrix_form
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

      /// Adds volumetric matrix form.
      void add_matrix_form(MatrixFormVol<Scalar>* mfv);
      /// Adds surface matrix form.
      void add_matrix_form_surf(MatrixFormSurf<Scalar>* mfs);
      /// Adds volumetric vector form.
      void add_vector_form(VectorFormVol<Scalar>* vfv);
      /// Adds surface vector form.
      void add_vector_form_surf(VectorFormSurf<Scalar>* vfs);

      /// Adds multicomponent volumetric matrix form.
      void add_multicomponent_matrix_form(MultiComponentMatrixFormVol<Scalar>* mfv);
      /// Adds multicomponent surface matrix form.
      void add_multicomponent_matrix_form_surf(MultiComponentMatrixFormSurf<Scalar>* mfs);
      /// Adds multicomponent volumetric vector form.
      void add_multicomponent_vector_form(MultiComponentVectorFormVol<Scalar>* vfv);
      /// Adds multicomponent surface vector form.
      void add_multicomponent_vector_form_surf(MultiComponentVectorFormSurf<Scalar>* vfs);

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
      virtual double get_current_time() const;

    protected:
      double current_time;
      unsigned int neq;
      int seq;
      bool is_matfree;

      struct Area  { Hermes::vector<std::string> markers;  };

      Hermes::vector<Area> areas;

    public:
      /// Holds volumetric matrix forms.
      Hermes::vector<MatrixFormVol<Scalar> *> mfvol;
      /// Holds surface matrix forms.
      Hermes::vector<MatrixFormSurf<Scalar> *> mfsurf;
      /// Holds volumetric vector forms.
      Hermes::vector<VectorFormVol<Scalar> *> vfvol;
      /// Holds surface vector forms.
      Hermes::vector<VectorFormSurf<Scalar> *> vfsurf;

      /// Holds multicomponent volumetric matrix forms.
      Hermes::vector<MultiComponentMatrixFormVol<Scalar> *> mfvol_mc;
      /// Holds multicomponent surface matrix forms.
      Hermes::vector<MultiComponentMatrixFormSurf<Scalar> *> mfsurf_mc;
      /// Holds multicomponent volumetric vector forms.
      Hermes::vector<MultiComponentVectorFormVol<Scalar> *> vfvol_mc;
      /// Holds multicomponent surface vector forms.
      Hermes::vector<MultiComponentVectorFormSurf<Scalar> *> vfsurf_mc;


      void get_stages(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*>& u_ext,
        std::vector<Stage<Scalar> >& stages, bool want_matrix, bool want_vector);

      bool** get_blocks(bool force_diagonal_blocks);

      bool is_in_area(std::string marker, std::string area) const
      { return area == marker; }

      bool is_in_areas(std::string marker, Hermes::vector<std::string> areas) const
      { 
        bool found = false;
        for (unsigned int i = 0; i < areas.size(); i++) {
          if (areas[i] == marker) {
            found = true;
            break;
          }
        }
        return found;
      }

      bool is_sym() const { return false; /* not impl. yet */ }

      friend class DiscreteProblem<Scalar>;
      friend class Hermes::Preconditioners::Precond<Scalar>;

      /// To be called only by the constructor of DiscreteProblem.
      void set_markers_conversion(Mesh::ElementMarkersConversion* element_markers_conversion, 
        Mesh::BoundaryMarkersConversion* boundary_markers_conversion)
      {
        this->element_markers_conversion = element_markers_conversion;
        this->boundary_markers_conversion = boundary_markers_conversion;
      }

    protected:

      Stage<Scalar>* find_stage(std::vector<Stage<Scalar> >& stages, int ii, int jj, Mesh* m1, Mesh* m2,
        Hermes::vector<MeshFunction<Scalar>*>& ext, Hermes::vector<Solution<Scalar>*>& u_ext);

      Stage<Scalar>* find_stage(std::vector<Stage<Scalar> >& stages, Hermes::vector<std::pair<unsigned int, unsigned int> > coordinates,
        Mesh* m1, Mesh* m2,
        Hermes::vector<MeshFunction<Scalar>*>& ext, Hermes::vector<Solution<Scalar>*>& u_ext);
      Stage<Scalar>* find_stage(std::vector<Stage<Scalar> >& stages, Hermes::vector<unsigned int> coordinates,
        Mesh* m1, Mesh* m2,
        Hermes::vector<MeshFunction<Scalar>*>& ext, Hermes::vector<Solution<Scalar>*>& u_ext);

      Mesh::ElementMarkersConversion* element_markers_conversion;
      Mesh::BoundaryMarkersConversion* boundary_markers_conversion;
    };

    template<typename Scalar>
    class HERMES_API Form
    {
    public:
      /// One area constructor.
      Form(std::string area = HERMES_ANY, Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);
      /// Multiple areas constructor.
      Form(Hermes::vector<std::string> areas, Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      inline void set_weakform(WeakForm<Scalar>* wf) { this->wf = wf; }

      std::string area;
      Hermes::vector<std::string> areas;
      Hermes::vector<MeshFunction<Scalar>*> ext;
      Hermes::vector<Scalar> param;
      /// Form will be always multiplied (scaled) with this number.
      double scaling_factor;
      /// External solutions for this form will start
      /// with u_ext[u_ext_offset] where u_ext[] are external
      /// solutions coming to the assembling procedure via the
      /// external coefficient vector.
      int u_ext_offset;
      /// If true, the form will be evaluated using adaptive
      /// numerical integration.
      bool adapt_eval;
      /// To obtain reference value, the element is split into
      /// four sons. In addition, the order is increased by this value.
      int adapt_order_increase;
      /// Max. allowed relative error (stopping criterion for adaptive
      /// numerical quadrature.
      double adapt_rel_error_tol;
      /// For time-dependent right-hand side functions.
      /// E.g. for Runge-Kutta methods. Otherwise the one time for the whole WeakForm can be used.
      void set_current_stage_time(double time);
      double get_current_stage_time() const;

    protected:
      WeakForm<Scalar>* wf;
      double stage_time;
    };

    template<typename Scalar>
    class HERMES_API MatrixFormVol : public Form<Scalar>
    {
    public:
      /// One area constructor.
      MatrixFormVol(unsigned int i, unsigned int j, 
        std::string area = HERMES_ANY, SymFlag sym = HERMES_NONSYM, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);
      /// Multiple areas constructor..
      MatrixFormVol(unsigned int i, unsigned int j, 
        Hermes::vector<std::string> areas, SymFlag sym = HERMES_NONSYM, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual MatrixFormVol* clone();

      unsigned int i, j;
      int sym;

      virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const;
      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;
    };

    template<typename Scalar>
    class HERMES_API MatrixFormSurf : public Form<Scalar>
    {
    public:
      /// One area constructor.
      MatrixFormSurf(unsigned int i, unsigned int j, std::string area = HERMES_ANY, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);
      /// Multiple areas constructor..
      MatrixFormSurf(unsigned int i, unsigned int j, Hermes::vector<std::string> areas, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual MatrixFormSurf* clone();

      unsigned int i, j;

      virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const;
      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;
    };

    template<typename Scalar>
    class VectorFormVol : public Form<Scalar>
    {
    public:
      /// One area constructor.
      VectorFormVol(unsigned int i, std::string area = HERMES_ANY, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);
      /// Multiple areas constructor..
      VectorFormVol(unsigned int i, Hermes::vector<std::string> areas, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual VectorFormVol* clone();

      unsigned int i;

      virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<Scalar> *ext) const;
      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, 
        ExtData<Hermes::Ord> *ext) const;
    };

    template<typename Scalar>
    class VectorFormSurf : public Form<Scalar>
    {
    public:
      /// One area constructor.
      VectorFormSurf(unsigned int i, std::string area = HERMES_ANY, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);
      /// Multiple areas constructor..
      VectorFormSurf(unsigned int i, Hermes::vector<std::string> areas, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual VectorFormSurf* clone();

      unsigned int i;

      virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<Scalar> *ext) const;
      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v, 
        Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;
    };

    /// Multi-component forms.
    /// The principle of functioning of multicomponent forms is as follows.
    /// The form is registered on coordinates in the vector 'coordinates', e.g. {[0,0], [0,1], [2,1]}.
    /// The method 'value' then accepts the parameter 'result', which is a vector with resulting values for
    /// the form on one coordinate. The vectors coordinates and result must have 1-1 relationship, i.e. if the
    /// form is to be registered on three different coordinates, the method value must insert 3 values for these
    /// three coordinates in the same order.
    template<typename Scalar>
    class HERMES_API MultiComponentMatrixFormVol : public Form<Scalar>
    {
    public:
      /// One area constructor.
      MultiComponentMatrixFormVol(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates,  
        std::string area = HERMES_ANY, SymFlag sym = HERMES_NONSYM,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);
      /// Multiple areas constructor..
      MultiComponentMatrixFormVol(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates,  
        Hermes::vector<std::string> areas, SymFlag sym = HERMES_NONSYM,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual MultiComponentMatrixFormVol* clone();

      Hermes::vector<std::pair<unsigned int, unsigned int> > coordinates;
      int sym;
      virtual void value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const = 0;
      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const = 0;
    };

    template<typename Scalar>
    class HERMES_API MultiComponentMatrixFormSurf : public Form<Scalar>
    {
    public:
      /// One area constructor.
      MultiComponentMatrixFormSurf(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, 
        std::string area = HERMES_ANY, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);
      /// Multiple areas constructor..
      MultiComponentMatrixFormSurf(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, 
        Hermes::vector<std::string> areas,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual MultiComponentMatrixFormSurf* clone();

      Hermes::vector<std::pair<unsigned int, unsigned int> > coordinates;

      virtual void value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const = 0;
      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const = 0;
    };

    template<typename Scalar>
    class HERMES_API MultiComponentVectorFormVol : public Form<Scalar>
    {
    public:
      /// One area constructor.
      MultiComponentVectorFormVol(Hermes::vector<unsigned int> coordinates, 
        std::string area = HERMES_ANY, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);
      /// Multiple areas constructor..
      MultiComponentVectorFormVol(Hermes::vector<unsigned int> coordinates, 
        Hermes::vector<std::string> areas,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual MultiComponentVectorFormVol* clone();
      Hermes::vector<unsigned int> coordinates;

      virtual void value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const = 0;
      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, 
        ExtData<Hermes::Ord> *ext) const = 0;
    };

    template<typename Scalar>
    class HERMES_API MultiComponentVectorFormSurf : public Form<Scalar>
    {
    public:
      /// One area constructor.
      MultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates, 
        std::string area = HERMES_ANY, 
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);
      /// Multiple areas constructor..
      MultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates, 
        Hermes::vector<std::string> areas,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        Hermes::vector<Scalar> param = Hermes::vector<Scalar>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual MultiComponentVectorFormSurf* clone();

      Hermes::vector<unsigned int> coordinates;

      virtual void value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const = 0;
      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v, 
        Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const = 0;
    };

    template<typename Scalar>
    class HERMES_API Stage
    {
    public:
      Stage() {};
      Hermes::vector<int> idx;
      Hermes::vector<Mesh*> meshes;
      Hermes::vector<Transformable*> fns;
      Hermes::vector<MeshFunction<Scalar>*> ext;

      /// Holds volumetric matrix forms.
      Hermes::vector<MatrixFormVol<Scalar> *> mfvol;
      /// Holds surface matrix forms.
      Hermes::vector<MatrixFormSurf<Scalar> *> mfsurf;
      /// Holds volumetric vector forms.
      Hermes::vector<VectorFormVol<Scalar> *> vfvol;
      /// Holds surface vector forms.
      Hermes::vector<VectorFormSurf<Scalar> *> vfsurf;

      /// Holds multicomponent volumetric matrix forms.
      Hermes::vector<MultiComponentMatrixFormVol<Scalar> *> mfvol_mc;
      /// Holds multicomponent surface matrix forms.
      Hermes::vector<MultiComponentMatrixFormSurf<Scalar> *> mfsurf_mc;
      /// Holds multicomponent volumetric vector forms.
      Hermes::vector<MultiComponentVectorFormVol<Scalar> *> vfvol_mc;
      /// Holds multicomponent surface vector forms.
      Hermes::vector<MultiComponentVectorFormSurf<Scalar> *> vfsurf_mc;

      std::set<int> idx_set;
      std::set<unsigned> seq_set;
      std::set<MeshFunction<Scalar>*> ext_set;
    };
  }
}
#endif
