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
// You should have received a copy of the GNU General Public Licenserix
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

    template<typename Scalar> class Form;
    template<typename Scalar> class OGProjection;
    template<typename Scalar> class MatrixFormVol;
    template<typename Scalar> class VectorFormVol;
    template<typename Scalar> class MatrixFormSurf;
    template<typename Scalar> class VectorFormSurf;
    template<typename Scalar> class MatrixFormDG;
    template<typename Scalar> class VectorFormDG;

    /// Bilinear form symmetry flag, see WeakForm::add_matrix_form
    enum SymFlag
    {
      HERMES_ANTISYM = -1,
      HERMES_NONSYM = 0,
      HERMES_SYM = 1
    };

    /// \brief Represents the weak formulation of a PDE problem.
    ///
    /// The WeakForm class represents the weak formulation of a system of linear PDEs.<br>
    /// The number of equations ("neq") in the system is fixed and is passed to the constructor.<br>
    /// The weak formulation of the system A(U, V) = L(V) has a block structure. A(U, V) is<br>
    /// a (neq x neq) matrix of bilinear forms a_mn(u, v) and L(V) is a neq-component vector<br>
    /// of linear forms l(v). U and V are the vectors of basis and test functions.<br>
    /// <br>
    /// There is a single tutorial just on implementing the weak formulation.<br>
    /// For some basic ideas, please see the examples provided.<br>
    template<typename Scalar>
    class HERMES_API WeakForm : public Hermes::Mixins::IntegrableWithGlobalOrder, public Hermes::Mixins::Loggable
    {
    public:
      /// Constructor.
      ///
      /// \param[in]  neq   Number of equations.
      ///
      /// \param[in]  mat_free    If this weak formulation does not include a matrix - e.g. JFNK method.
      WeakForm(unsigned int neq = 1, bool mat_free = false);

      /// Destructor.
      virtual ~WeakForm();

      /// Adds volumetric matrix form.
      void add_matrix_form(MatrixFormVol<Scalar>* mfv);

      /// Adds surface matrix form.
      void add_matrix_form_surf(MatrixFormSurf<Scalar>* mfs);

      /// Adds DG matrix form.
      void add_matrix_form_DG(MatrixFormDG<Scalar>* mfDG);

      /// Adds volumetric vector form.
      void add_vector_form(VectorFormVol<Scalar>* vfv);

      /// Adds surface vector form.
      void add_vector_form_surf(VectorFormSurf<Scalar>* vfs);

      /// Adds DG vector form.
      void add_vector_form_DG(VectorFormDG<Scalar>* vfDG);

      /// Provides possibility of setup element-wise parameters.
      /// For parameters that only depend on element and that do
      /// not have to be calculated for every form.
      /// This is rarely used and typically only for multi-physical tasks where there is a multitude of forms.
      virtual void set_active_state(Element** e);

      /// Provides possibility of setup edge-wise parameters.
      /// For parameters that only depend on element and edge and that do
      /// not have to be calculated for every form.
      /// This is rarely used and typically only for multi-physical tasks where there is a multitude of forms.
      virtual void set_active_edge_state(Element** e, int isurf);

      /// Provides possibility of setup edge-wise parameters.
      /// For parameters that only depend on element and inner edge and that do
      /// not have to be calculated for every form.
      /// This is rarely used and typically only for multi-physical tasks where there is a multitude of forms.
      virtual void set_active_DG_state(Element** e, int isurf);

      /// Returns the number of equations.
      unsigned int get_neq() const { return neq; }

      /// This weakform is matrix-free.
      bool is_matrix_free() const { return is_matfree; }

      /// For time-dependent right-hand side functions.
      /// Sets current time.
      void set_current_time(double time);

      /// For time-dependent right-hand side functions.
      /// Sets current time step.
      void set_current_time_step(double time_step);

      /// For time-dependent right-hand side functions.
      /// Gets current time.
      virtual double get_current_time() const;

      /// For time-dependent right-hand side functions.
      /// Gets current time step.
      virtual double get_current_time_step() const;

      /// External functions.
      /// Set one external function.
      void set_ext(MeshFunction<Scalar>* ext);

      /// External functions.
      /// Set more external functions.
      void set_ext(Hermes::vector<MeshFunction<Scalar>*> ext);

      /// External functions.
      /// Get external functions.
      Hermes::vector<MeshFunction<Scalar>*> get_ext() const;

      /// Cloning.
      virtual WeakForm* clone() const;

      /// Internal.
      Hermes::vector<Form<Scalar> *> get_forms() const;
      Hermes::vector<MatrixFormVol<Scalar> *> get_mfvol() const;
      Hermes::vector<MatrixFormSurf<Scalar> *> get_mfsurf() const;
      Hermes::vector<MatrixFormDG<Scalar> *> get_mfDG() const;
      Hermes::vector<VectorFormVol<Scalar> *> get_vfvol() const;
      Hermes::vector<VectorFormSurf<Scalar> *> get_vfsurf() const;
      Hermes::vector<VectorFormDG<Scalar> *> get_vfDG() const;

      /// Deletes all volumetric and surface forms.
      void delete_all();

    protected:
      /// External solutions.
      Hermes::vector<MeshFunction<Scalar>*> ext;

      double current_time;
      double current_time_step;

      unsigned int neq;

      bool is_matfree;

      /// Holds all forms.
      Hermes::vector<Form<Scalar> *> forms;

      /// Holds volumetric matrix forms.
      Hermes::vector<MatrixFormVol<Scalar> *> mfvol;

      /// Holds surface matrix forms.
      Hermes::vector<MatrixFormSurf<Scalar> *> mfsurf;

      /// Holds DG matrix forms.
      Hermes::vector<MatrixFormDG<Scalar> *> mfDG;

      /// Holds volumetric vector forms.
      Hermes::vector<VectorFormVol<Scalar> *> vfvol;

      /// Holds surface vector forms.
      Hermes::vector<VectorFormSurf<Scalar> *> vfsurf;

      /// Holds DG vector forms.
      Hermes::vector<VectorFormDG<Scalar> *> vfDG;

      bool** get_blocks(bool force_diagonal_blocks) const;

      friend class DiscreteProblem<Scalar>;
      friend class DiscreteProblemLinear<Scalar>;
      friend class RungeKutta<Scalar>;
      friend class OGProjection<Scalar>;
      friend class Hermes::Preconditioners::Precond<Scalar>;

      bool warned_nonOverride;

      // internal.
      virtual void cloneMembers(const WeakForm<Scalar>* otherWf);

    private:
      void free_ext();
    };

    /// \brief Abstract, base class for any form - i.e. integral in the weak formulation of (a system of) PDE<br>
    /// By default, the form is initialized with the following natural attributes:<br>
    /// - no external functions.<br>
    /// - assembled over any (parameter 'HERMES_ANY') element/boundary marker.<br>
    /// Internal.
    template<typename Scalar>
    class HERMES_API Form
    {
    public:
      /// Constructor with coordinates.
      Form();
      virtual ~Form();

      /// get-set methods
      /// areas
      void set_area(std::string area);
      void set_areas(Hermes::vector<std::string> areas);
      Hermes::vector<std::string> getAreas() const;

      /// external functions - dual functionality with the overall WeakForm.
      /// For Agros, this approach is better in some way, for e.g. Euler equations,
      /// the other one (for the whole WeakForm) is faster.
      void set_ext(MeshFunction<Scalar>* ext);
      void set_ext(Hermes::vector<MeshFunction<Scalar>*> ext);
      Hermes::vector<MeshFunction<Scalar>*> get_ext() const;

    protected:
      /// Set pointer to a WeakForm.
      inline void set_weakform(WeakForm<Scalar>* wf) { this->wf = wf; }

      /// Markers of the areas where this form will be assembled.
      Hermes::vector<std::string> areas;

      /// External solutions for this form will start
      /// with u_ext[u_ext_offset] where u_ext[] are external
      /// solutions coming to the assembling procedure via the
      /// external coefficient vector.
      int u_ext_offset;

      /// External solutions.
      Hermes::vector<MeshFunction<Scalar>*> ext;

      /// For time-dependent right-hand side functions.
      /// E.g. for Runge-Kutta methods. Otherwise the one time for the whole WeakForm can be used.
      void set_current_stage_time(double time);

      double get_current_stage_time() const;

      /// Form will be always multiplied (scaled) with this number.
      double scaling_factor;

      WeakForm<Scalar>* wf;
      double stage_time;
      void setScalingFactor(double scalingFactor);
      void set_uExtOffset(int u_ext_offset);
      friend class WeakForm<Scalar>;
      friend class RungeKutta<Scalar>;
      friend class DiscreteProblem<Scalar>;
      friend class DiscreteProblemLinear<Scalar>;
    };

    /// \brief Abstract, base class for matrix form - i.e. a single integral in the bilinear form on the left hand side of the variational formulation of a (system of) PDE.<br>
    /// By default, the matrix form is initialized with the following natural attribute:<br>
    /// - nonsymmetrical (if the user omits the HERMES_SYM / HERMES_ANTISYM parameters, nothing worse than a non-necessary calculations happen).
    template<typename Scalar>
    class HERMES_API MatrixForm : public Form<Scalar>
    {
    public:
      /// Constructor with coordinates.
      MatrixForm(unsigned int i, unsigned int j);

      virtual ~MatrixForm();

      unsigned int i;
      unsigned int j;
      unsigned int previous_iteration_space_index;

      SymFlag sym;

      typedef Scalar valueFunction(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const;

      virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const;

      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

    protected:
      friend class DiscreteProblem<Scalar>;
    };

    /// \brief Abstract, base class for matrix Volumetric form - i.e. MatrixForm, where the integration is with respect to 2D-Lebesgue measure (elements).
    template<typename Scalar>
    class HERMES_API MatrixFormVol : public MatrixForm<Scalar>
    {
    public:
      /// Constructor with coordinates.
      MatrixFormVol(unsigned int i, unsigned int j);

      void setSymFlag(SymFlag sym);
      SymFlag getSymFlag() const;

      virtual ~MatrixFormVol();

      virtual MatrixFormVol* clone() const;
    };

    /// \brief Abstract, base class for matrix Surface form - i.e. MatrixForm, where the integration is with respect to 1D-Lebesgue measure (element domain-boundary edges).
    template<typename Scalar>
    class HERMES_API MatrixFormSurf : public MatrixForm<Scalar>
    {
    public:
      /// Constructor with coordinates.
      MatrixFormSurf(unsigned int i, unsigned int j);

      virtual ~MatrixFormSurf();

      virtual MatrixFormSurf* clone() const;
    };

    /// \brief Abstract, base class for matrix DG form - i.e. MatrixForm, where the integration is with respect to 1D-Lebesgue measure (element inner-domain edges).
    template<typename Scalar>
    class HERMES_API MatrixFormDG : public MatrixForm<Scalar>
    {
    public:
      /// Constructor with coordinates.
      MatrixFormDG(unsigned int i, unsigned int j);

      virtual ~MatrixFormDG();

      virtual MatrixFormDG* clone() const;
    };

    /// \brief Abstract, base class for vector form - i.e. a single integral in the linear form on the right hand side of the variational formulation of a (system of) PDE.
    template<typename Scalar>
    class VectorForm : public Form<Scalar>
    {
    public:
      /// Constructor with coordinates.
      VectorForm(unsigned int i);

      virtual ~VectorForm();

      virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const;

      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
        Func<Ord> **ext) const;
      unsigned int i;

    protected:
      friend class DiscreteProblem<Scalar>;
    };

    /// \brief Abstract, base class for vector Volumetric form - i.e. VectorForm, where the integration is with respect to 2D-Lebesgue measure (elements).
    template<typename Scalar>
    class VectorFormVol : public VectorForm<Scalar>
    {
    public:
      /// Constructor with coordinates.
      VectorFormVol(unsigned int i);

      virtual ~VectorFormVol();

      virtual VectorFormVol* clone() const;
    };

    /// \brief Abstract, base class for vector Surface form - i.e. VectorForm, where the integration is with respect to 1D-Lebesgue measure (element domain-boundary edges).
    template<typename Scalar>
    class VectorFormSurf : public VectorForm<Scalar>
    {
    public:
      /// Constructor with coordinates.
      VectorFormSurf(unsigned int i);

      virtual ~VectorFormSurf();

      virtual VectorFormSurf* clone() const;
    };

    /// \brief Abstract, base class for vector DG form - i.e. VectorForm, where the integration is with respect to 1D-Lebesgue measure (element inner-domain edges).
    template<typename Scalar>
    class VectorFormDG : public VectorForm<Scalar>
    {
    public:
      /// Constructor with coordinates.
      VectorFormDG(unsigned int i);

      virtual ~VectorFormDG();

      virtual VectorFormDG* clone() const;
    };
  }
}
#endif
