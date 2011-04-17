#include "weakform/weakform.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

class CustomWeakFormNeutronics : public WeakForm
{ 
public:
  CustomWeakFormNeutronics(Hermes::vector<std::string> layers, 
                                    Hermes::vector<double> D_map, 
                                    Hermes::vector<double> Sigma_a_map, 
                                    Hermes::vector<double> Sources_map)
  : WeakForm(1) {
    
    for (unsigned int i = 0; i < layers.size(); i++)
    {
      // Diffusion.
      add_matrix_form(new DefaultLinearDiffusion(0, 0, layers[i], D_map[i], HERMES_SYM));
    
      // Absorption.
      add_matrix_form(new DefaultLinearMass(0, 0, layers[i], Sigma_a_map[i], HERMES_SYM));
      
      // Sources.
      add_vector_form(new DefaultVectorFormConst(0, layers[i], Sources_map[i]));
    }
  }
};
