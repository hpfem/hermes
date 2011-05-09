////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "hermes2d.h"

using namespace WeakFormsNeutronics::Multigroup::CompleteWeakForms::Diffusion; 

class CustomWeakForm : public DefaultWeakFormSourceIteration
{
  public:
    CustomWeakForm(const MaterialPropertyMaps& matprop,
                   Hermes::vector<MeshFunction*>& iterates,
                   double init_keff, std::string bdy_vacuum);
};



//////  Filters and functionals needed for the K_eff eigenvalue iteration. ////////////////////

class SourceFilter : public SimpleFilter
{
  public: 
    SourceFilter(Hermes::vector<MeshFunction*> solutions, MaterialPropertyMaps& matprop,
                 const std::string& source_area)
      : SimpleFilter(solutions, Hermes::vector<int>())
    {
      nu = matprop.get_nu().at(source_area);
      Sigma_f = matprop.get_Sigma_f().at(source_area);
    }
    SourceFilter(Hermes::vector<Solution*> solutions, MaterialPropertyMaps& matprop,
                 const std::string& source_area)
    : SimpleFilter(solutions, Hermes::vector<int>())
    {
      nu = matprop.get_nu().at(source_area);
      Sigma_f = matprop.get_Sigma_f().at(source_area);
    }
    
  private:
    rank1 nu;
    rank1 Sigma_f;
    
    void filter_fn(int n, Hermes::vector<scalar*> values, scalar* result);
};

// Integral over the active core.
double integrate(MeshFunction* sln, std::string area);