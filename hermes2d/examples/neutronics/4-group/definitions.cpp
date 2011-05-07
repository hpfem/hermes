////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "weakform_library/weakforms_neutronics.h"
#include "function/filter.h"

using namespace WeakFormsNeutronics::Multigroup::CompleteWeakForms::Diffusion; 

class CustomWeakForm : public DefaultWeakFormSourceIteration
{
  public:
    CustomWeakForm(const MaterialPropertyMaps& matprop,
                   Hermes::vector<MeshFunction*>& iterates,
                   double init_keff, 
                   std::string bdy_vacuum)
      : DefaultWeakFormSourceIteration(matprop, iterates, init_keff, HERMES_AXISYM_Y)
    {
      for (unsigned int g = 0; g < matprop.get_G(); g++)
      {
        add_matrix_form_surf(new VacuumBoundaryCondition::Jacobian(g, bdy_vacuum, HERMES_AXISYM_Y));
        add_vector_form_surf(new VacuumBoundaryCondition::Residual(g, bdy_vacuum, HERMES_AXISYM_Y));
      }
    }
};



//////  Filters and functionals needed for the K_eff eigenvalue iteration. ////////////////////

class SourceFilter : public SimpleFilter
{
  public: 
    SourceFilter(Hermes::vector<MeshFunction*> solutions, MaterialPropertyMaps& matprop)
      : SimpleFilter(solutions, Hermes::vector<int>())
    {
      nu = matprop.get_nu().at(core);
      Sigma_f = matprop.get_Sigma_f().at(core);
    }
    SourceFilter(Hermes::vector<Solution*> solutions, MaterialPropertyMaps& matprop)
    : SimpleFilter(solutions, Hermes::vector<int>())
    {
      nu = matprop.get_nu().at(core);
      Sigma_f = matprop.get_Sigma_f().at(core);
    }
    
  private:
    rank1 nu;
    rank1 Sigma_f;
    
    void filter_fn(int n, Hermes::vector<scalar*> values, scalar* result)
    {
      for (int i = 0; i < n; i++) {
        result[i] = 0;
        for (unsigned int j = 0; j < values.size(); j++)
          result[i] += nu[j] * Sigma_f[j] * values.at(j)[i];
      }
    } 
};

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

// Reference k_effective reactor eigenvalue for the material properties below
// and geometry from the file 'reactor.mesh'. For this example, it was obtained
// simplistically by a reference calculation on a 3x uniformly refined mesh
// with uniform distribution of polynomial degrees (=4), with convergence
// tolerance set to 5e-11.
const double REF_K_EFF = 1.1409144;



//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const MaterialPropertyMap1 D = material_property_map<rank1>
(
  reflector, grow(0.0164)(0.0085)(0.00832)(0.00821)
)(
  core, grow(0.0235)(0.0121)(0.0119)(0.0116)
);

const MaterialPropertyMap1 Sa = material_property_map<rank1>
(
  reflector, grow(0.00139)(0.000218)(0.00197)(0.0106)
)(
  core, grow(0.00977)(0.162)(0.156)(0.535)
);

const MaterialPropertyMap1 Sr = material_property_map<rank1>
(
  reflector, grow(1.77139)(0.533218)(3.31197)(0.0106)
)(
  core, grow(1.23977)(0.529)(2.436)(0.535)
);

const MaterialPropertyMap1 Sf = material_property_map<rank1>
(
  reflector, grow(0.0)(0.0)(0.0)(0.0)
)(
  core, grow(0.00395)(0.0262)(0.0718)(0.346)
);

const MaterialPropertyMap1 nu = material_property_map<rank1>
(
  reflector, grow(0.0)(0.0)(0.0)(0.0) 
)(
  core, grow(2.49)(2.43)(2.42)(2.42)
);

const MaterialPropertyMap1 chi = material_property_map<rank1>
(
  reflector, grow(0.0)(0.0)(0.0)(0.0)
)(
  core, grow(0.9675)(0.03250)(0.0)(0.0)
);

const MaterialPropertyMap2 Ss = material_property_map<rank2>
(
  reflector,
  gmat
  (
    grow(0.0)(0.0)(0.0)(0.0)
  )(
    grow(1.77)(0.0)(0.0)(0.0)
  )(
    grow(0.0)(0.533)(0.0)(0.0)
  )(
    grow(0.0)(0.0)(3.31)(0.0)
  )
)(
  core,
  gmat
  (
    grow(0.0)(0.0)(0.0)(0.0)
  )(
    grow(1.23)(0.0)(0.0)(0.0)
  )(
    grow(0.0)(0.367)(0.0)(0.0)
  )(
    grow(0.0)(0.0)(2.28)(0.0)
  )
);

const bool2 Ss_nnz = bool2
(
  bool_mat
  (
    bool_row(0)(0)(0)(0)
  )(
    bool_row(1)(0)(0)(0)
  )(
    bool_row(0)(1)(0)(0)
  )(
    bool_row(0)(0)(1)(0)
  )
);