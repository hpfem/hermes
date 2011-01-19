#ifndef KELLY_TYPE_ADAPT_H
#define KELLY_TYPE_ADAPT_H

#include "adapt.h"
#include "../neighbor.h"

enum ElementAccumulationMethod
{
  ACCUMULATE_BY_ADDITION,
  ACCUMULATE_MAX_VALUE
};

inline double original_kelly_scaling_factor(double e_diam)
{
  return e_diam/24.;
}
inline double scale_by_element_diameter(double e_diam)
{
  return e_diam;
}

template<typename Real, typename Scalar>
Scalar original_kelly_interface_estimator(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, 
                                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * sqr( e->nx[i] * (u->get_dx_central(i) - u->get_dx_neighbor(i)) +
                           e->ny[i] * (u->get_dy_central(i) - u->get_dy_neighbor(i))  );
  return result;
}

class HERMES_API KellyTypeAdapt : public Adapt
{
  void accumulate(double *accumulated, double element_contribution, ElementAccumulationMethod accum_type)
  {
    if (accum_type == ACCUMULATE_BY_ADDITION) *accumulated += element_contribution;
    else *accumulated = std::max(*accumulated, element_contribution);
  }
  
  typedef double (*scaling_factor_t)(double e_diam);
  struct ErrorEstimatorForm
  {  
    int i, area;          
    WeakForm::vector_form_val_t fn;  
    WeakForm::vector_form_ord_t ord;  
    Hermes::vector<MeshFunction *> ext;
  };
  
  double eval_volumetric_estimator(KellyTypeAdapt::ErrorEstimatorForm* err_est_form, RefMap* rm);
  double eval_boundary_estimator(ErrorEstimatorForm* err_est_form, RefMap* rm, SurfPos* surf_pos);
  double eval_interface_estimator(KellyTypeAdapt::ErrorEstimatorForm* err_est_form, RefMap* rm, SurfPos* surf_pos, NeighborSearch* nbs);
  double eval_estimator_normalization(matrix_form_val_t val, matrix_form_ord_t ord, RefMap* rm, MeshFunction* sln);
    
  protected:
    /// Linear forms to calculate the error estimator for each component.
    Hermes::vector<ErrorEstimatorForm> error_estimators_vol;     
    Hermes::vector<ErrorEstimatorForm> error_estimators_surf;
    
    ElementAccumulationMethod estimator_normalization_accum_types[H2D_MAX_COMPONENTS];    
    ElementAccumulationMethod total_norm_accum_type;
    
    scaling_factor_t interface_scaling_factor;
    bool ignore_visited_segments;
    
    virtual double calc_err_internal(Hermes::vector< Solution* > slns, 
                                     Hermes::vector< double >* component_errors,  
                                     unsigned int error_flags);
  public:
    KellyTypeAdapt(Hermes::vector<Space *> spaces_, 
                   Hermes::vector<ProjNormType> norms_ = Hermes::vector<ProjNormType>(),
                   scaling_factor_t interface_scaling_factor_ = original_kelly_scaling_factor,
                   bool ignore_visited_segments = true);
                   
    virtual ~KellyTypeAdapt()
    {
      error_estimators_surf.clear();
      error_estimators_vol.clear();
    }
    
    void add_error_form_vol(int i, 
                            WeakForm::vector_form_val_t vfv, WeakForm::vector_form_ord_t vfo, 
                            int area = HERMES_ANY,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());
    void add_error_form_vol(WeakForm::vector_form_val_t vfv, WeakForm::vector_form_ord_t vfo, 
                            int area = HERMES_ANY,
                            Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>())
    { 
      add_error_form_vol(0, vfv, vfo, area, ext); 
    }
    
    void add_error_form_surf(int i, 
                             WeakForm::vector_form_val_t vfv, WeakForm::vector_form_ord_t vfo,
                             int area = H2D_DG_INNER_EDGE,
                             Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>());
    void add_error_form_surf(WeakForm::vector_form_val_t vfv, WeakForm::vector_form_ord_t vfo,
                             int area = H2D_DG_INNER_EDGE,
                             Hermes::vector<MeshFunction*>ext = Hermes::vector<MeshFunction*>()) 
    { 
      add_error_form_surf(0, vfv, vfo, area, ext); 
    }
    
    void set_normalization_form(int i,
                                matrix_form_val_t bi_fn, matrix_form_ord_t bi_ord,
                                ElementAccumulationMethod accum_type);
    void set_normalization_form(matrix_form_val_t bi_fn, matrix_form_ord_t bi_ord,
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
