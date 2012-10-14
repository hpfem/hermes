#include "prev_solution.h"
      
      
				void PrevSolution::set_own_mesh(const Mesh* mesh){						
						if(this->mesh == mesh){
							Mesh* new_mesh = new Mesh;
							new_mesh->copy(mesh);
							this->mesh = new_mesh;
							own_mesh = true;						
						}else throw Hermes::Exceptions::Exception("Solution mesh unequal own_mesh.");
				}
    
           MeshFunction<double>* PrevSolution::clone(){
                if(this->get_type() == HERMES_SLN)
									return Solution<double>::clone();
								PrevSolution* sln = new PrevSolution(this->mesh);
								return sln;         
          
          };
    

  
    void PrevSolution::free()
    {
      if(mono_coeffs  != NULL) { delete [] mono_coeffs;   mono_coeffs = NULL;  }
      if(elem_orders != NULL) { delete [] elem_orders;  elem_orders = NULL; }
      if(dxdy_buffer != NULL) { delete [] dxdy_buffer;  dxdy_buffer = NULL; }

      for (int i = 0; i < this->num_components; i++)
        if(elem_coeffs[i] != NULL)
        { delete [] elem_coeffs[i];  elem_coeffs[i] = NULL; }
        
 				if((own_mesh == true)&&(mesh!=NULL)) delete mesh;
 				mesh = NULL;
 				
        e_last = NULL;

        free_tables();
    }




  
