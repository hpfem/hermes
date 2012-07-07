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

#include "../function/solution.h"
#include <string>

namespace Hermes
{
  namespace Hermes2D
  {
    /// Geometrical type of weak forms.
    enum GeomType
    {
      HERMES_PLANAR = 0,         // Planar problem.
      HERMES_AXISYM_X = 1,       // Axisymmetric problem where x-axis is the axis of symmetry.
      HERMES_AXISYM_Y = 2        // Axisymmetric problem where y-axis is the axis of symmetry.
    };

    class RefMap;
    template<typename Scalar> class DiscreteProblem;
    template<typename Scalar> class DiscreteProblemLinear;
    template<typename Scalar> class RungeKutta;
    template<typename Scalar> class Space;
    template<typename Scalar> class MeshFunction;
    struct SurfPos;

    class Element;
    class Shapeset;
    template<typename T> class Func;
    template<typename T> class Geom;
    template<typename T> class ExtData;

    template<typename Scalar> class Form;
    template<typename Scalar> class MatrixFormVol;
    template<typename Scalar> class VectorFormVol;
    template<typename Scalar> class MatrixFormSurf;
    template<typename Scalar> class VectorFormSurf;

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
    /// The weak formulation of the system A(U, V) = L(V) has a block structure. A(U, V) is
    /// a (neq x neq) matrix of bilinear forms a_mn(u, v) and L(V) is a neq-component vector
    /// of linear forms l(v). U and V are the vectors of basis and test functions.
    ///
    template<typename Scalar>
    class HERMES_API WeakForm : public Hermes::Mixins::Loggable
    {
    public:

      /// Constructor.
      ///
      /// \param[in]  neq   Number of equations.
      ///
      /// \param[in]  mat_free    If this weak formulation does not include a matrix - e.g. JFNK method.
      WeakForm(unsigned int neq = 1, bool mat_free = false);

      /// Destructor.
      ~WeakForm();

      /// Adds volumetric matrix form.
      void add_matrix_form(MatrixFormVol<Scalar>* mfv);

      /// Adds surface matrix form.
      void add_matrix_form_surf(MatrixFormSurf<Scalar>* mfs);

      /// Adds volumetric vector form.
      void add_vector_form(VectorFormVol<Scalar>* vfv);

      /// Adds surface vector form.
      void add_vector_form_surf(VectorFormSurf<Scalar>* vfs);

      /// Returns the number of equations.
      unsigned int get_neq() const { return neq; }

      bool is_matrix_free() const { return is_matfree; }

      /// For time-dependent right-hand side functions.
      void set_current_time(double time);

      virtual double get_current_time() const;

      Hermes::vector<MatrixFormVol<Scalar> *> get_mfvol();
      Hermes::vector<MatrixFormSurf<Scalar> *> get_mfsurf();
      Hermes::vector<VectorFormVol<Scalar> *> get_vfvol();
      Hermes::vector<VectorFormSurf<Scalar> *> get_vfsurf();

      /// Deletes all volumetric and surface forms.
      void delete_all();

    protected:
      /// Internal. Used by DiscreteProblem to detect changes in the weakform.
      int get_seq() const { return seq; }

      bool** get_blocks(bool force_diagonal_blocks) const;

      double current_time;

      unsigned int neq;

      int seq;

      bool is_matfree;

      /// Holds volumetric matrix forms.
      Hermes::vector<MatrixFormVol<Scalar> *> mfvol;

      /// Holds surface matrix forms.
      Hermes::vector<MatrixFormSurf<Scalar> *> mfsurf;

      /// Holds volumetric vector forms.
      Hermes::vector<VectorFormVol<Scalar> *> vfvol;

      /// Holds surface vector forms.
      Hermes::vector<VectorFormSurf<Scalar> *> vfsurf;

      friend class DiscreteProblem<Scalar>;
      friend class DiscreteProblemLinear<Scalar>;
      friend class RungeKutta<Scalar>;
      friend class Hermes::Preconditioners::Precond<Scalar>;
    };

    template<typename Scalar>
    class HERMES_API Form
    {
    public:
      /// One area constructor.
      Form(std::string area = HERMES_ANY, Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual ~Form() {};

      /// Multiple areas constructor.
      Form(Hermes::vector<std::string> areas, Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      inline void set_weakform(WeakForm<Scalar>* wf) { this->wf = wf; }

      /// For time-dependent right-hand side functions.
      /// E.g. for Runge-Kutta methods. Otherwise the one time for the whole WeakForm can be used.
      void set_current_stage_time(double time);

      double get_current_stage_time() const;

      Hermes::vector<std::string> areas;

      Hermes::vector<MeshFunction<Scalar>*> ext;

      /// Form will be always multiplied (scaled) with this number.
      double scaling_factor;

      /// External solutions for this form will start
      /// with u_ext[u_ext_offset] where u_ext[] are external
      /// solutions coming to the assembling procedure via the
      /// external coefficient vector.
      int u_ext_offset;

      unsigned int i;

    protected:
      WeakForm<Scalar>* wf;
      double stage_time;
    };

    template<typename Scalar>
    class HERMES_API MatrixForm : public Form<Scalar>
    {
    public:
      /// One area constructor.
      MatrixForm(unsigned int i, unsigned int j,
        std::string area = HERMES_ANY, Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      /// Multiple areas constructor..
      MatrixForm(unsigned int i, unsigned int j,
        Hermes::vector<std::string> areas, Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual ~MatrixForm() {};

      unsigned int j;

      int sym;

      virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const;

      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;
    };

    template<typename Scalar>
    class HERMES_API MatrixFormVol : public MatrixForm<Scalar>
    {
    public:
      /// One area constructor.
      MatrixFormVol(unsigned int i, unsigned int j,
        std::string area = HERMES_ANY, SymFlag sym = HERMES_NONSYM,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      /// Multiple areas constructor..
      MatrixFormVol(unsigned int i, unsigned int j,
        Hermes::vector<std::string> areas, SymFlag sym = HERMES_NONSYM,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual ~MatrixFormVol() {};

      virtual MatrixFormVol* clone();
    };

    template<typename Scalar>
    class HERMES_API MatrixFormSurf : public MatrixForm<Scalar>
    {
    public:
      /// One area constructor.
      MatrixFormSurf(unsigned int i, unsigned int j, std::string area = HERMES_ANY,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      /// Multiple areas constructor..
      MatrixFormSurf(unsigned int i, unsigned int j, Hermes::vector<std::string> areas,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual ~MatrixFormSurf() {};

      virtual MatrixFormSurf* clone();
    };

    template<typename Scalar>
    class VectorForm : public Form<Scalar>
    {
    public:
      /// One area constructor.
      VectorForm(unsigned int i, std::string area = HERMES_ANY,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      /// Multiple areas constructor..
      VectorForm(unsigned int i, Hermes::vector<std::string> areas,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual ~VectorForm() {};

      virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const;

      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
        ExtData<Hermes::Ord> *ext) const;
    };

    template<typename Scalar>
    class VectorFormVol : public VectorForm<Scalar>
    {
    public:
      /// One area constructor.
      VectorFormVol(unsigned int i, std::string area = HERMES_ANY,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      /// Multiple areas constructor..
      VectorFormVol(unsigned int i, Hermes::vector<std::string> areas,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual ~VectorFormVol() {};

      virtual VectorFormVol* clone();
    };

    template<typename Scalar>
    class VectorFormSurf : public VectorForm<Scalar>
    {
    public:
      /// One area constructor.
      VectorFormSurf(unsigned int i, std::string area = HERMES_ANY,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      /// Multiple areas constructor..
      VectorFormSurf(unsigned int i, Hermes::vector<std::string> areas,
        Hermes::vector<MeshFunction<Scalar>*> ext = Hermes::vector<MeshFunction<Scalar>*>(),
        double scaling_factor = 1.0, int u_ext_offset = 0);

      virtual ~VectorFormSurf() {};

      virtual VectorFormSurf* clone();
    };
  }
}
#endif