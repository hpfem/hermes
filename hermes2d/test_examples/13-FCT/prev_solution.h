#ifndef __PREV_SOLUTION_H
#define __PREV_SOLUTION_H
#include "hermes2d.h"
using namespace Hermes;
using namespace Hermes::Solvers;
using namespace Hermes::Algebra;
using namespace Hermes::Hermes2D;

    class PrevSolution : public Solution<double>
    {
    public:
      PrevSolution() : Solution<double>(){ };
       PrevSolution(MeshSharedPtr mesh): Solution<double>(mesh){
						own_mesh = false;     
       };
    
    ~PrevSolution()
    {
    }
    
    virtual MeshFunction<double>* clone();      
		void set_own_mesh(const MeshSharedPtr mesh);					

    
    
  protected:  
  bool own_mesh;  
    void free();




    };
    
    
    
#endif
