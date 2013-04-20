#ifndef __HP_ADAPT_H
#define __HP_ADAPT_H
#include "hermes2d.h"


using namespace Hermes;
using namespace Hermes::Hermes2D;


class  HPAdapt : public Adapt<double>
{
public:
	  HPAdapt(SpaceSharedPtr<double> space): Adapt<double>(space, new DefaultErrorCalculator<double, HERMES_L2_NORM>(RelativeErrorToGlobalNorm, 1))
    {	
	  }; 
	~HPAdapt(){};
 

bool adapt_smooth(int* smooth_dof, int max_p);


};







#endif
