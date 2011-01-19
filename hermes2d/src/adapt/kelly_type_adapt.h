#ifndef KELLY_TYPE_ADAPT_H
#define KELLY_TYPE_ADAPT_H

#include "adapt.h"
#include "../neighbor.h"

// The following enumeration defines 
//    1/ how to accumulate values computed at individual quadrature points into an element integral result
//    2/ element-wise values into a full domain-integral value
// ACCUMULATE_BY_ADDITION defines the standard accumulation procedure used for calculating L^p norms with
//    1 <= p < oo.
// ACCUMULATE_MAX_VALUE may be used for obtaining the L^oo norm.
//
enum ElementAccumulationMethod
{
  ACCUMULATE_BY_ADDITION,
  ACCUMULATE_MAX_VALUE
};

/// Pre-defined function used for scaling interface error estimates (see the KellyTypeAdapt constructor).
inline double scale_by_element_diameter(double e_diam)
{
  return e_diam;
}

/// Type of a pointer to the interface estimator scaling function.
typedef double (*interface_estimator_scaling_fn_t)(double e_diam);

class HERMES_API KellyTypeAdapt : public Adapt
{
  protected:
    /// Accumulates \c current_contribution into \c accumulated according to the specified \c accum_type.
    void accumulate(double *accumulated, double current_contribution, ElementAccumulationMethod accum_type)
    {
      if (accum_type == ACCUMULATE_BY_ADDITION) *accumulated += current_contribution;
      else *accumulated = std::max(*accumulated, current_contribution);
    }
    
    ///
    /// Functions used for evaluating the actual error estimator forms for an active element or edge segment.
    ///
    double eval_volumetric_estimator(WeakForm::VectorFormVol* err_est_form, 
                                    RefMap* rm);
    double eval_boundary_estimator(WeakForm::VectorFormSurf* err_est_form, RefMap* rm, 
                                  SurfPos* surf_pos);
    double eval_interface_estimator(WeakForm::VectorFormSurf* err_est_form, 
                                    RefMap* rm, SurfPos* surf_pos, NeighborSearch* nbs);
    double eval_estimator_normalization(error_matrix_form_val_t val, error_matrix_form_ord_t ord, 
                                        RefMap* rm, MeshFunction* sln);

    ///
    /// Linear forms used to calculate the error estimator value for each component.
    ///
    Hermes::vector<WeakForm::VectorFormVol> error_estimators_vol;
    Hermes::vector<WeakForm::VectorFormSurf> error_estimators_surf;
    
    ElementAccumulationMethod estimator_normalization_accum_types[H2D_MAX_COMPONENTS];
    ElementAccumulationMethod total_norm_accum_type;

    interface_estimator_scaling_fn_t interface_scaling_fn;
    double interface_scaling_const;
    double volumetric_scaling_const;
    double boundary_scaling_const;
    
    bool ignore_visited_segments;

    virtual double calc_err_internal(Hermes::vector< Solution* > slns,
                                     Hermes::vector< double >* component_errors,
                                     unsigned int error_flags);
  public:
    KellyTypeAdapt(Hermes::vector<Space *> spaces_,
                   Hermes::vector<ProjNormType> norms_ = Hermes::vector<ProjNormType>(),
                   interface_estimator_scaling_fn_t interface_scaling_fn_ = scale_by_element_diameter,
                   bool ignore_visited_segments = true);

    virtual ~KellyTypeAdapt()
    {
      error_estimators_surf.clear();
      error_estimators_vol.clear();
    }

    void add_error_form_vol(int i,
                            WeakForm::vector_form_val_t vfv, WeakForm::vector_form_ord_t vfo,
                            int area = HERMES_ANY,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>(),
                            double scaling_factor = 0.0);
    void add_error_form_vol(WeakForm::vector_form_val_t vfv, WeakForm::vector_form_ord_t vfo,
                            int area = HERMES_ANY,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>(),
                            double scaling_factor = 0.0)
    {
      add_error_form_vol(0, vfv, vfo, area, ext, scaling_factor);
    }

    void add_error_form_surf(int i,
                             WeakForm::vector_form_val_t vfv, WeakForm::vector_form_ord_t vfo,
                             int area = H2D_DG_INNER_EDGE,
                             Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>(),
                             double scaling_factor = 0.0);
    void add_error_form_surf(WeakForm::vector_form_val_t vfv, WeakForm::vector_form_ord_t vfo,
                             int area = H2D_DG_INNER_EDGE,
                             Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>(),
                             double scaling_factor = 0.0)
    {
      add_error_form_surf(0, vfv, vfo, area, ext, scaling_factor);
    }

    void set_normalization_form(int i,
                                error_matrix_form_val_t bi_fn, error_matrix_form_ord_t bi_ord,
                                ElementAccumulationMethod accum_type);
    void set_normalization_form(error_matrix_form_val_t bi_fn, error_matrix_form_ord_t bi_ord,
                                ElementAccumulationMethod accum_type)
    {
      set_normalization_form(0, bi_fn, bi_ord, accum_type);
    }

    void set_total_norm_accumulation_type(ElementAccumulationMethod total_norm_accum_type_)
    {
      this->total_norm_accum_type = total_norm_accum_type_;
    }

    double calc_err_est(Solution *sln,
                        unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL)
    {
      if (num != 1) EXIT("Wrong number of solutions.");
      return calc_err_est(Hermes::vector<Solution *> (sln),
                          (Hermes::vector<double>*) NULL, error_flags);
    }

    double calc_err_est(Hermes::vector<Solution *> slns,
                        Hermes::vector<double>* component_errors = NULL,
                        unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL)
    {
      return calc_err_internal(slns, component_errors, error_flags);
    }

    bool adapt(double thr, int strat = 0, int regularize = -1, double to_be_processed = 0.0);
};

#endif // KELLY_TYPE_ADAPT_H
