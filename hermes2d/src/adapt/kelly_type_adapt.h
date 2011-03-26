#ifndef KELLY_TYPE_ADAPT_H
#define KELLY_TYPE_ADAPT_H

#include "adapt.h"
#include "../neighbor.h"
#include "../discrete_problem.h"

// #ifdef KELLY_TYPE_ADAPT_H_IS_REWORKED

/// Pre-defined function used for scaling interface error estimates (see the KellyTypeAdapt constructor).
inline double scale_by_element_diameter(double e_diam)
{
  return e_diam;
}

/// Type of a pointer to the interface estimator scaling function.
typedef double (*interface_estimator_scaling_fn_t)(double e_diam);

/// \class KellyTypeAdapt
/// \ingroup g_adapt
/// \brief A framework for explicit aposteriori error estimators.
///
/// Explicit error estimators estimate the error of approximate solution on an element by evaluating
/// element residuals and jumps of the solution across element edges ([2]). A typical example is
/// the Kelly error estimator ([1]) where a sum of the L2 norms of element residual and jumps of
/// solution gradients across the element boundaries defines the element error.
///
/// References:
///   [1] Kelly D. W., Gago O. C., Zienkiewicz O. C., Babuska I.:
///       A posteriori error analysis and adaptive processes in the finite element method: Part I—error analysis.
///       Int. J. Numer. Methods Engng. 1983;19:1593–619.
///   [2] Gratsch T., Bathe K. J.:
///       A posteriori error estimation techniques in practical finite element analysis.
///       Computers and Structures 83 (2005) 235–265.
///   [3] Zienkiewicz O. C., Taylor R. L., Zhu J. Z.:
///       The finite element method: its basis and fundamentals (Section 13.7.1).
///       6th ed. (2005), Elsevier.
///
class HERMES_API KellyTypeAdapt : public Adapt
{
public:
  /// Class representing the weak form of an error estimator.
  /// 
  /// A user must derive his own representation of the estimator from this class (an example is provided by the class
  /// \c BasicKellyAdapt below). The three attributes have the following meaning:
  ///
  ///   - i     ... with a multi-component solution, this defines for which component this estimate applies,
  ///   - area  ... defines in which geometric parts of the domain should the estimate be used - e.g. by defining 
  ///               area = H2D_DG_INNER_EDGE, errors at element interfaces will be tracked by the estimator,
  ///   - ext   ... vector with external functions possibly used within the estimator (e.g. previous time-level 
  ///               solutions appearing in the residual) - currently not used.
  ///
  /// Every estimator form must also implement the two methods \c ord and \c value, which are used for determining
  /// the integration order and for the actual evaluation of the form, respectively. During the evaluation, their 
  /// parameters will be interpreted as follows:
  ///
  ///   - int n,                 ... number of integration points in the currently processed element
  ///   - double *wt,            ... corresponding integration weights
  ///   - Func\<scalar\> *u[],   ... all solution components
  ///   - Func\<double\> *u,     ... currently processed solution component
  ///   - Geom\<double\> *e,     ... geometric data of the currently processed element
  ///   - ExtData\<scalar\> *ext ... external functions (currently unused).
  ///
  class HERMES_API ErrorEstimatorForm
  {
  public:
    int i;
    std::string area;
    Hermes::vector<MeshFunction *> ext;
    Hermes::vector<scalar> param;

    ErrorEstimatorForm(int i, std::string area = HERMES_ANY, 
                       Hermes::vector<MeshFunction *> ext = Hermes::vector<MeshFunction *>(), 
                       Hermes::vector<scalar> param = Hermes::vector<scalar>()) :
    i(i), area(area), ext(ext), param(param) {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                         Func<scalar> *u, Geom<double> *e,
                         ExtData<scalar> *ext) const
    {
      error("KellyTypeAdapt::ErrorEstimatorForm::value() must be overridden.");
      return 0.0;
    }
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                    Func<Ord> *u, Geom<Ord> *e,
                    ExtData<Ord> *ext) const
    {
      error("KellyTypeAdapt::ErrorEstimatorForm::ord() must be overridden.");
      return Ord();
    }

    // FIXME - temporary
    KellyTypeAdapt *adapt;
  };

  protected:  
    DiscreteProblem dp; // Only needed for gaining access to NeighborSearch methods.

    ///
    /// Functions used for evaluating the actual error estimator forms for an active element or edge segment.
    ///
    double eval_volumetric_estimator(KellyTypeAdapt::ErrorEstimatorForm* err_est_form,
                                    RefMap* rm);
    double eval_boundary_estimator(KellyTypeAdapt::ErrorEstimatorForm* err_est_form,
                                   RefMap* rm,
                                   SurfPos* surf_pos);
    double eval_interface_estimator(KellyTypeAdapt::ErrorEstimatorForm* err_est_form,
                                    RefMap *rm,
                                    SurfPos* surf_pos,
                                    LightArray<NeighborSearch*>& neighbor_searches,
                                    int neighbor_index);
    double eval_solution_norm(Adapt::MatrixFormVolError* form,
                              RefMap* rm,
                              MeshFunction* sln);

    ///
    /// Linear forms used to calculate the error estimator value for each component.
    ///
    Hermes::vector<KellyTypeAdapt::ErrorEstimatorForm *> error_estimators_vol;
    Hermes::vector<KellyTypeAdapt::ErrorEstimatorForm *> error_estimators_surf;

    Mesh::ElementMarkersConversion element_markers_conversion;
    Mesh::BoundaryMarkersConversion boundary_markers_conversion;

    /// Scaling of the interface error estimates. May be specified by the user during construction.
    ///
    Hermes::vector<interface_estimator_scaling_fn_t> interface_scaling_fns;
    bool use_aposteriori_interface_scaling; ///< Specifies whether the interface error estimators for each
                                            ///< component will be multiplied by \c interface_scaling_fns
                                            ///< after being evaluated.

    ///
    /// Constant scaling. Reserved for the derived classes, not to be used by the user explicitly.
    ///
    double interface_scaling_const;   ///< Constant scaling of the boundary error estimates.
    double volumetric_scaling_const;  ///< Constant scaling of the volumetric error estimates (like the residual norm).
    double boundary_scaling_const;    ///< Constant scaling of the boundary error estimates.

    /// Specifies whether the interface error estimator will be evaluated from each side of each interface
    /// (when <c>ignore_visited_segments == false</c> ), or only once for each interface
    /// (<c>ignore_visited_segments == true</c>).
    bool ignore_visited_segments;

    /// Calculates error estimates for each solution component, the total error estimate, and possibly also
    /// their normalizations. If called with a pair of solutions, the version from Adapt is used (this is e.g.
    /// done when comparing approximate solution to the exact one - in this case, we do not want to compute 
    /// the Kelly estimator value, but rather the ordinary difference between the solutions).
    virtual double calc_err_internal(Hermes::vector<Solution *> slns,
                                     Hermes::vector<double>* component_errors,
                                     unsigned int error_flags);

  public:

    /// Constructor.
    ///
    /// \param[in]  spaces_   Approximation space of each solution component.
    ///
    /// \param[in]  norms_    Norms used for making relative error estimates.
    ///                       If not specified, they are defined according to the spaces.
    ///
    /// \param[in]  ignore_visited_segments_ If true, error estimator for each inner edge will be evaluated only
    ///                                      once. It will be added to the total error estimate for both the active
    ///                                      element and its neighbors across that edge, after possible scaling by
    ///                                      \c interface_scaling_fns_ for the current component (with the diameter of
    ///                                      the appropriate element). This saves duplicate evaluations with same
    ///                                      results when the estimator is given e.g. by the jumps of the solution.
    ///
    ///                                      If false, error estimator for each surface of each element will be
    ///                                      evaluated, regardless of whether the neighbor side of the interface
    ///                                      has already been processed.
    ///
    ///                                      Note that if \c interface_scaling_fns_ is empty (or unspecified) then the
    ///                                      default scaling by element diameter will be always performed unless it is
    ///                                      switched off by a call to \c disable_aposteriori_interface_scaling.
    ///
    /// \param[in]  interface_scaling_fns_  Specifies functions used for scaling the interface error estimator for
    ///                                     each component. The scale is defined as a real function of the element
    ///                                     diameter and multiplies the result of the interface estimators.
    ///                                     It may thus be already present in the interface estimator forms
    ///                                     themselves, in which case call \c disable_aposteriori_interface_scaling.
    ///                                     In this case, it may also be required that \c ignore_visited_segments be
    ///                                     false in order to always ensure that the diameter belongs to the element
    ///                                     whose error is being calculated.
    ///
    KellyTypeAdapt(Hermes::vector<Space *> spaces_,
                   Hermes::vector<ProjNormType> norms_ = Hermes::vector<ProjNormType>(),
                   bool ignore_visited_segments_ = true,
                   Hermes::vector<interface_estimator_scaling_fn_t>
                        interface_scaling_fns_ = Hermes::vector<interface_estimator_scaling_fn_t>());

    KellyTypeAdapt(Space * space_,
                   ProjNormType norm_ = HERMES_UNSET_NORM,
                   bool ignore_visited_segments_ = true,
                   interface_estimator_scaling_fn_t interface_scaling_fn_ = NULL);
                   
    /// Destructor.
    virtual ~KellyTypeAdapt()
    {
      error_estimators_surf.clear();
      error_estimators_vol.clear();
    }

    Mesh::ElementMarkersConversion* get_element_markers_conversion() { return &this->element_markers_conversion; };
    Mesh::BoundaryMarkersConversion* get_boundary_markers_conversion() { return &this->boundary_markers_conversion; };

    /// Append volumetric error estimator form.
    ///
    /// For example, element residual norms may be represented by such a form.
    ///
    /// \param[in]  form ... object representing the form. A class derived from \c KellyTypeAdapt::ErrorEstimatorForm 
    ///                      defines its datatype.
    ///
    void add_error_estimator_vol(KellyTypeAdapt::ErrorEstimatorForm* form);

    /// Append boundary or interface error estimator form.
    ///
    /// Interface form is defined by <c> form::area == H2D_DG_INNER_EDGE </c>. The effective types for \c u_ext, \c u 
    /// and \c e (three of the obligatory parameters of form::value() and form::ord()) will then be, respectively 
    /// \c DiscontinuousFunc*[], \c DiscontinuousFunc* and \c InterfaceGeom*.
    ///
    void add_error_estimator_surf(KellyTypeAdapt::ErrorEstimatorForm* form); 
    
    ///
    /// The following two methods calculate the error of the given \c sln, using \code calc_err_internal \endcode.
    ///

    double calc_err_est(Solution *sln,
                        unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL)
    {
      if (num != 1) EXIT("Wrong number of solutions.");
      Hermes::vector<Solution *> slns;
      slns.push_back(sln);
      return calc_err_est(slns, NULL, error_flags);
    }


    double calc_err_est(Hermes::vector<Solution *> slns,
                        Hermes::vector<double>* component_errors = NULL,
                        unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL)
    {
      return calc_err_internal(slns, component_errors, error_flags);
    }

    /// Refines the elements selected by the \code RefinementSelectors::HOnlySelector \endcode according
    /// to the errors calculated by \code calc_err_est \endcode.
    ///
    bool adapt(double thr, int strat = 0, int regularize = -1, double to_be_processed = 0.0);

    void disable_aposteriori_interface_scaling() { use_aposteriori_interface_scaling = false; }
};

/// \class BasicKellyAdapt
/// \ingroup g_adapt
/// \brief Simple Kelly-estimator based adaptivity for elliptic problems.
///
/// Original error estimator that Kelly et. al. ([1]) derived for the Laplace equation with constant
/// coefficient, approximated on a quadrilateral mesh. The error of each element is estimated by the
/// L2 norm of jumps of gradients across element faces (the contribution of the residual norm is
/// relatively insignificant and is neglected, see [3]). Note that the estimator has been successfully
/// used also for other problems than that for which it had been originally derived.
///
/// TODO: Add handling of boundary conditions.
///       Currently, the forms for the Neumann and Newton boundary conditions must be specified by
///       the user, see the example \c poisson-kelly-adapt.
///

class HERMES_API BasicKellyAdapt : public KellyTypeAdapt
{
public:
  class HERMES_API ErrorEstimatorFormKelly : public KellyTypeAdapt::ErrorEstimatorForm
  {
  public:
    ErrorEstimatorFormKelly(int i) : ErrorEstimatorForm(i, H2D_DG_INNER_EDGE) {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                         Func<scalar> *u, Geom<double> *e,
                         ExtData<scalar> *ext) const
    {
      return original_kelly_interface_estimator<double, scalar>(n, wt, u_ext, u, e, ext);
    }
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                    Func<Ord> *u, Geom<Ord> *e,
                    ExtData<Ord> *ext) const
    {
      return original_kelly_interface_estimator<Ord, Ord>(n, wt, u_ext, u, e, ext);
    }

  private:
    template<typename Real, typename Scalar>
    static Scalar original_kelly_interface_estimator(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u,
                                                     Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0.;
      for (int i = 0; i < n; i++)
        result += wt[i] * sqr( e->nx[i] * (u->get_dx_central(i) - u->get_dx_neighbor(i)) +
                               e->ny[i] * (u->get_dy_central(i) - u->get_dy_neighbor(i))  );
      return result;
    }
  };
  
  /// Constructor.
  ///
  /// For the equation \f$ K \Delta u = f \f$, the argument \c const_by_laplacian is equal to \$ K \$.
  ///
  BasicKellyAdapt(Hermes::vector<Space *> spaces_,
                  Hermes::vector<ProjNormType> norms_ = Hermes::vector<ProjNormType>(),
                  double const_by_laplacian = 1.0) : KellyTypeAdapt(spaces_, norms_)
  {
    set_scaling_consts(const_by_laplacian);   
    for (int i = 0; i < num; i++)
      this->error_estimators_surf.push_back(new ErrorEstimatorFormKelly(i));
  }
  
  BasicKellyAdapt(Space* space_, ProjNormType norm_ = HERMES_UNSET_NORM, double const_by_laplacian = 1.0) : KellyTypeAdapt(space_, norm_) 
  {
    set_scaling_consts(const_by_laplacian);
    this->error_estimators_surf.push_back(new ErrorEstimatorFormKelly(0));
  }
  
private:
  void set_scaling_consts(double C)
  {
    interface_scaling_const = 1./(24.*C);
    volumetric_scaling_const = interface_scaling_const;
    boundary_scaling_const = interface_scaling_const;
  }
};

// #endif
#endif // KELLY_TYPE_ADAPT_H
