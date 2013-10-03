#include "prev_solution.h"
      
      
				void PrevSolution::set_own_mesh(const MeshSharedPtr mesh){						
						if(this->mesh == mesh)
            {
							MeshSharedPtr new_mesh(new Mesh);
							new_mesh->copy(mesh);
							this->mesh = new_mesh;
							own_mesh = true;						
						}
            else
              throw Hermes::Exceptions::Exception("Solution mesh unequal own_mesh.");
				}
    
           MeshFunction<double>* PrevSolution::clone(){
                if(this->get_type() == HERMES_SLN)
									return Solution<double>::clone();
								PrevSolution* sln = new PrevSolution(this->mesh);
								return sln;         
          
          };
    

  
    void PrevSolution::free()
    {
      if(mono_coeffs  != nullptr) { delete [] mono_coeffs;   mono_coeffs = nullptr;  }
      if(elem_orders != nullptr) { delete [] elem_orders;  elem_orders = nullptr; }
      if(dxdy_buffer != nullptr) { delete [] dxdy_buffer;  dxdy_buffer = nullptr; }

      for (int i = 0; i < this->num_components; i++)
        if(elem_coeffs[i] != nullptr)
        { delete [] elem_coeffs[i];  elem_coeffs[i] = nullptr; }
        
 				
        e_last = nullptr;

        free_tables();
    }




  
