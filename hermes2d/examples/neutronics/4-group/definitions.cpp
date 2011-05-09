////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialPropertyMaps& matprop,
                                Hermes::vector<MeshFunction*>& iterates,
                                double init_keff, std::string bdy_vacuum )
  : DefaultWeakFormSourceIteration(matprop, iterates, init_keff, HERMES_AXISYM_Y)
{
  for (unsigned int g = 0; g < matprop.get_G(); g++)
  {
    add_matrix_form_surf(new VacuumBoundaryCondition::Jacobian(g, bdy_vacuum, HERMES_AXISYM_Y));
    add_vector_form_surf(new VacuumBoundaryCondition::Residual(g, bdy_vacuum, HERMES_AXISYM_Y));
  }
}

//////  Filters and functionals needed for the K_eff eigenvalue iteration. ////////////////////

void SourceFilter::filter_fn(int n, Hermes::vector<scalar*> values, scalar* result)
{
  for (int i = 0; i < n; i++) 
  {
    result[i] = 0;
    for (unsigned int j = 0; j < values.size(); j++)
      result[i] += nu[j] * Sigma_f[j] * values.at(j)[i];
  }
} 

// Integral over the active core.
double integrate(MeshFunction* sln, std::string area)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  
  double integral = 0.0;
  Element* e;
  Mesh* mesh = sln->get_mesh();
  int marker = mesh->get_element_markers_conversion().get_internal_marker(area);
  
  for_all_active_elements(e, mesh)
  {
    if (e->marker == marker)
    {
      update_limit_table(e->get_mode());
      sln->set_active_element(e);
      RefMap* ru = sln->get_refmap();
      int o = sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o);
      sln->set_quad_order(o, H2D_FN_VAL);
      scalar *uval = sln->get_fn_values();
      double* x = ru->get_phys_x(o);
      double result = 0.0;
      h1_integrate_expression(x[i] * uval[i]);
      integral += result;
    }
  }
  
  return 2.0 * M_PI * integral;
}