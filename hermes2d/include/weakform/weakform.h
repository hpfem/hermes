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

#include "../function/exact_solution.h"

namespace Hermes
{
  namespace Hermes2D
  {
#pragma region forward-declarations
    class RefMap;
    template<typename Scalar> class DiscreteProblem;
    template<typename Scalar> class DiscreteProblemSelectiveAssembler;
    template<typename Scalar> class DiscreteProblemIntegrationOrderCalculator;
    template<typename Scalar> class RungeKutta;
    template<typename Scalar> class Space;
    template<typename Scalar> class MeshFunction;
    struct SurfPos;

    class Element;
    class Shapeset;
    template<typename T> class Func;
    template<typename T> class DiscontinuousFunc;
    template<typename T> class GeomVol;
    template<typename T> class GeomSurf;

    template<typename Scalar> class Form;
    template<typename Scalar> class OGProjection;
    template<typename Scalar> class MatrixFormVol;
    template<typename Scalar> class VectorFormVol;
    template<typename Scalar> class MatrixFormSurf;
    template<typename Scalar> class VectorFormSurf;
    template<typename Scalar> class MatrixFormDG;
    template<typename Scalar> class VectorFormDG;
#pragma endregion

    /// \brief Used to pass the instances of WeakForm around.
    template<typename Scalar>
    class HERMES_API WeakFormSharedPtr : public std::tr1::shared_ptr<Hermes::Hermes2D::WeakForm<Scalar> >
    {
    public:
      WeakFormSharedPtr(Hermes::Hermes2D::WeakForm<Scalar>* ptr = nullptr);

      WeakFormSharedPtr(const WeakFormSharedPtr<Scalar>& other);

      void operator=(const WeakFormSharedPtr<Scalar>& other);
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

      enum FormIntegrationDimension
      {
        FormVol = 0,
        FormSurf = 1
      };

      enum FormEquationSide
      {
        MatrixForm = 0,
        VectorForm = 1
      };

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
      virtual void set_active_edge_state(Element** e, unsigned char isurf);

      /// Provides possibility of setup edge-wise parameters.
      /// For parameters that only depend on element and inner edge and that do
      /// not have to be calculated for every form.
      /// This is rarely used and typically only for multi-physical tasks where there is a multitude of forms.
      virtual void set_active_DG_state(Element** e, unsigned char isurf);

      /// Returns the number of equations.
      inline unsigned char get_neq() const { return neq; }

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
      /// Set one function acting on the u_ext functions in assembling (for fast assembling of nonlinear problems).
      /// IMPORTANT: This function will appear at the beginning of the Func<Scalar>** ext array in the value(), and ord() methods of individual forms.
      void set_u_ext_fn(UExtFunctionSharedPtr<Scalar> ext);

      /// External functions.
      /// Set one external function.
      /// IMPORTANT: This function will appear at the END (after those functions coming via set_u_ext_fn) of the Func<Scalar>** ext array in the value(), and ord() methods of individual forms.
      void set_ext(MeshFunctionSharedPtr<Scalar> ext);

      /// External functions.
      /// Set functions acting on the u_ext functions in assembling (for fast assembling of nonlinear problems).
      /// IMPORTANT: These functions will appear at the beginning of the Func<Scalar>** ext array in the value(), and ord() methods of individual forms.
      void set_u_ext_fn(std::vector<UExtFunctionSharedPtr<Scalar> > ext);

      /// External functions.
      /// Set external functions.
      /// IMPORTANT: These functions will appear at the END (after those functions coming via set_u_ext_fn) of the Func<Scalar>** ext array in the value(), and ord() methods of individual forms.
      void set_ext(std::vector<MeshFunctionSharedPtr<Scalar> > ext);

      /// External functions.
      /// Get external functions.
      std::vector<MeshFunctionSharedPtr<Scalar> > get_ext() const;

      /// Cloning.
      virtual WeakForm* clone() const;

      // Checks presence of DG forms.
      bool is_DG() const;

      /// Internal.
      std::vector<Form<Scalar> *> get_forms() const;
      std::vector<MatrixFormVol<Scalar> *> get_mfvol() const;
      std::vector<MatrixFormSurf<Scalar> *> get_mfsurf() const;
      std::vector<MatrixFormDG<Scalar> *> get_mfDG() const;
      std::vector<VectorFormVol<Scalar> *> get_vfvol() const;
      std::vector<VectorFormSurf<Scalar> *> get_vfsurf() const;
      std::vector<VectorFormDG<Scalar> *> get_vfDG() const;

      /// Deletes all volumetric and surface forms.
      void delete_all();

    protected:
      /// External solutions.
      std::vector<MeshFunctionSharedPtr<Scalar> > ext;
      std::vector<UExtFunctionSharedPtr<Scalar> > u_ext_fn;

      double current_time;
      double current_time_step;

      /// Number of equations.
      unsigned int neq;

      /// Original number of equations in case this is a Runge-Kutta enlarged system.
      unsigned int original_neq;

      bool is_matfree;

      /// Holds all forms.
      std::vector<Form<Scalar> *> forms;

      /// Holds volumetric matrix forms.
      std::vector<MatrixFormVol<Scalar> *> mfvol;

      /// Holds surface matrix forms.
      std::vector<MatrixFormSurf<Scalar> *> mfsurf;

      /// Holds DG matrix forms.
      std::vector<MatrixFormDG<Scalar> *> mfDG;

      /// Holds volumetric vector forms.
      std::vector<VectorFormVol<Scalar> *> vfvol;

      /// Holds surface vector forms.
      std::vector<VectorFormSurf<Scalar> *> vfsurf;

      /// Holds DG vector forms.
      std::vector<VectorFormDG<Scalar> *> vfDG;

      bool** get_blocks(bool force_diagonal_blocks) const;

      friend class DiscreteProblem<Scalar>;
      friend class Form<Scalar>;
      friend class DiscreteProblemDGAssembler<Scalar>;
      friend class DiscreteProblemThreadAssembler<Scalar>;
      friend class DiscreteProblemIntegrationOrderCalculator<Scalar>;
      friend class DiscreteProblemSelectiveAssembler<Scalar>;
      friend class RungeKutta<Scalar>;
      friend class OGProjection<Scalar>;
      friend class Hermes::Preconditioners::Precond<Scalar>;

      // Internal.
      virtual void cloneMembers(const WeakFormSharedPtr<Scalar>& other_wf);
      // Internal.
      void cloneMemberExtFunctions(std::vector<MeshFunctionSharedPtr<Scalar> > source_ext, std::vector<MeshFunctionSharedPtr<Scalar> >& cloned_ext);

      // Internal - processes markers, translates from strings to ints.
      template<typename FormType>
      void processFormMarkers(const std::vector<SpaceSharedPtr<Scalar> > spaces, bool surface, std::vector<FormType> forms_to_process);
      void processFormMarkers(const std::vector<SpaceSharedPtr<Scalar> > spaces);

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
      Form(int i = 0);
      virtual ~Form();

      /// get-set methods
      /// areas
      void set_area(std::string area);
      void set_areas(std::vector<std::string> areas);
      std::vector<std::string> getAreas() const;

      /// external functions - dual functionality with the overall WeakForm.
      /// For Agros, this approach is better in some way, for e.g. Euler equations,
      /// the other one (for the whole WeakForm) is faster.
      /// External functions.
      /// Set one external function.
      void set_ext(MeshFunctionSharedPtr<Scalar> ext);
      void set_u_ext_fn(UExtFunctionSharedPtr<Scalar> ext);

      /// External functions.
      /// Set more external functions.
      void set_ext(std::vector<MeshFunctionSharedPtr<Scalar> > ext);
      void set_u_ext_fn(std::vector<UExtFunctionSharedPtr<Scalar> > ext);
      std::vector<MeshFunctionSharedPtr<Scalar> > get_ext() const;

      /// scaling factor
      void setScalingFactor(double scalingFactor);

      unsigned int i;

    protected:
      /// Set pointer to a WeakForm + handling of internal data.
      void set_weakform(WeakForm<Scalar>* wf);

      /// Markers of the areas where this form will be assembled.
      std::vector<std::string> areas;

      /// Internal - this structure is being filled anew with every assembling.
      std::vector<int> areas_internal;

      /// Internal - this structure is being filled anew with every assembling.
      /// True iff areas contain HERMES_ANY - meaning that this form represents an integral over the whole domain (whole boundary in case of surface forms).
      bool assembleEverywhere;

      /// External solutions for this form will start
      /// with u_ext[u_ext_offset] where u_ext[] are external
      /// solutions coming to the assembling procedure via the
      /// external coefficient vector.
      int u_ext_offset;

      /// When dealing with nonlinear problems of multiple equations, sometimes the nonlinearity has to access different quantity previous iterations. This selector
      /// servers for that purpose - according to this, the correct one will be selected - see e.g.DefaultJacobianDiffusion class.
      /// Defaults to 'i' for VectorForm and 'j' for MatrixForm. This can be changed in derived forms.
      unsigned int previous_iteration_space_index;

      /// External solutions.
      std::vector<MeshFunctionSharedPtr<Scalar> > ext;
      std::vector<UExtFunctionSharedPtr<Scalar> > u_ext_fn;

      double get_current_stage_time() const;

      WeakForm<Scalar>* wf;
    private:
      double stage_time;
      void set_uExtOffset(int u_ext_offset);
      /// Form will be always multiplied (scaled) with this number.
      double scaling_factor;
      /// For time-dependent right-hand side functions.
      /// E.g. for Runge-Kutta methods. Otherwise the one time for the whole WeakForm can be used.
      void set_current_stage_time(double time);
      /// Copy the basic data from other_form - used in cloning.
      void copy_base(Form<Scalar>* other_form);
      friend class WeakForm<Scalar>;
      friend class RungeKutta<Scalar>;
      friend class DiscreteProblem<Scalar>;
      friend class DiscreteProblemDGAssembler<Scalar>;
      friend class DiscreteProblemIntegrationOrderCalculator<Scalar>;
      friend class DiscreteProblemSelectiveAssembler<Scalar>;
      friend class DiscreteProblemThreadAssembler<Scalar>;
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

      unsigned int j;

      SymFlag sym;

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

      virtual Scalar value(int n, double *wt, Func<Scalar> **u_ext, Func<double> *u, Func<double> *v,
        GeomVol<double> *e, Func<Scalar> **ext) const;

      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> **u_ext, Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        GeomVol<Hermes::Ord> *e, Func<Ord> **ext) const;

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

      virtual Scalar value(int n, double *wt, Func<Scalar> **u_ext, Func<double> *u, Func<double> *v,
        GeomSurf<double> *e, Func<Scalar> **ext) const;

      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> **u_ext, Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        GeomSurf<Hermes::Ord> *e, Func<Ord> **ext) const;

      virtual MatrixFormSurf* clone() const;
    };

    /// \brief Abstract, base class for matrix DG form - i.e. bilinear form, where the integration is with respect to 1D-Lebesgue measure (element inner-domain edges).
    template<typename Scalar>
    class HERMES_API MatrixFormDG : public Form<Scalar>
    {
    public:
      /// Constructor with coordinates.
      MatrixFormDG(unsigned int i, unsigned int j);

      virtual ~MatrixFormDG();

      unsigned int j;

      virtual Scalar value(int n, double *wt, DiscontinuousFunc<Scalar> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v,
        GeomSurf<double> *e, DiscontinuousFunc<Scalar> **ext) const;

      virtual Hermes::Ord ord(int n, double *wt, DiscontinuousFunc<Hermes::Ord> **u_ext, DiscontinuousFunc<Hermes::Ord> *u, DiscontinuousFunc<Hermes::Ord> *v,
        GeomSurf<Hermes::Ord> *e, DiscontinuousFunc<Ord> **ext) const;

      virtual MatrixFormDG* clone() const;
    protected:
      friend class DiscreteProblem<Scalar>;
    };

    /// \brief Abstract, base class for vector form - i.e. a single integral in the linear form on the right hand side of the variational formulation of a (system of) PDE.
    template<typename Scalar>
    class VectorForm : public Form<Scalar>
    {
    public:
      /// Constructor with coordinates.
      VectorForm(unsigned int i);

      virtual ~VectorForm();

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

      virtual Scalar value(int n, double *wt, Func<Scalar> **u_ext, Func<double> *v,
        GeomVol<double> *e, Func<Scalar> **ext) const;

      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> **u_ext, Func<Hermes::Ord> *v, GeomVol<Hermes::Ord> *e,
        Func<Ord> **ext) const;

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

      virtual Scalar value(int n, double *wt, Func<Scalar> **u_ext, Func<double> *v,
        GeomSurf<double> *e, Func<Scalar> **ext) const;

      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> **u_ext, Func<Hermes::Ord> *v, GeomSurf<Hermes::Ord> *e,
        Func<Ord> **ext) const;

      virtual VectorFormSurf* clone() const;
    };

    /// \brief Abstract, base class for vector DG form - i.e. linear Form, where the integration is with respect to 1D-Lebesgue measure (element inner-domain edges).
    template<typename Scalar>
    class VectorFormDG : public Form<Scalar>
    {
    public:
      /// Constructor with coordinates.
      VectorFormDG(unsigned int i);

      virtual ~VectorFormDG();

      virtual Scalar value(int n, double *wt, DiscontinuousFunc<Scalar> **u_ext, Func<double> *v,
        GeomSurf<double> *e, DiscontinuousFunc<Scalar> **ext) const;

      virtual Hermes::Ord ord(int n, double *wt, DiscontinuousFunc<Hermes::Ord> **u_ext, Func<Hermes::Ord> *v, GeomSurf<Hermes::Ord> *e,
        DiscontinuousFunc<Ord> **ext) const;

      virtual VectorFormDG* clone() const;
    protected:
      friend class DiscreteProblem<Scalar>;
    };
  }
}
#endif
