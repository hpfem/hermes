#ifndef __H2D_NEUTRONICS_WEAK_FORMS_H
#define __H2D_NEUTRONICS_WEAK_FORMS_H

#include "h1.h"

/* Default weak form for neutron diffusion equations
   with Dirichlet and/or zero Neumann BC.

   Nonzero Neumann or Newton boundary conditions can be enabled
   by creating a descendant and adding surface forms to it.
*/

namespace WeakFormsNeutronDiffusion 
{ 
  /* Simple monoenergetic neutron diffusion, with the following weak formulation within each
     homogeneous region:
  
        \int_{region} D \nabla\phi \cdot \nabla\psi d\bfx + \int_{region} \Sigma_a \phi\psi d\bfx
          = \int_{region} Q_{ext}\psi d\bfx
    
     where 
    
        D         ... diffusion coefficient, 
        \Sigma_a  ... absorption cross-section, 
        Q_{ext}   ... external neutron sources 
      
     are region-wise constant physical parameters of the problem. Each region has one entry in vector
     'regions', which is the marker used for all elements it is composed of (usually specified in the
     mesh file). A corresponding entry in the *_map arguments is the value of the particular physical 
     parameter for that marker.
  */
  class DefaultWeakFormSimpleMonoenergetic : public WeakForm
  {        
    public:
      DefaultWeakFormSimpleMonoenergetic( Hermes::vector<std::string> regions, 
                                          Hermes::vector<double> D_map, 
                                          Hermes::vector<double> Sigma_a_map, 
                                          Hermes::vector<double> Q_map ) : WeakForm(1) 
      {
        using namespace WeakFormsH1::VolumetricMatrixForms;
        using namespace WeakFormsH1::VolumetricVectorForms;
        
        for (unsigned int i = 0; i < regions.size(); i++)
        {
          // Diffusion.
          add_matrix_form(new DefaultLinearDiffusion(0, 0, regions[i], D_map[i], HERMES_SYM));
        
          // Absorption.
          add_matrix_form(new DefaultLinearMass(0, 0, regions[i], Sigma_a_map[i], HERMES_SYM));
          
          // Sources.
          add_vector_form(new DefaultVectorFormConst(0, regions[i], Q_map[i]));
        }
      }
  };
}

#endif