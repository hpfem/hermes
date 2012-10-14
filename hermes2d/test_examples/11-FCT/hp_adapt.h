#ifndef __HP_ADAPT_H
#define __HP_ADAPT_H
#include "hermes2d.h"


using namespace Hermes;
using namespace Hermes::Hermes2D;


class  HPAdapt : public Adapt<double>
{
public:
	  HPAdapt(Space<double>* space, ProjNormType proj_norm): Adapt<double>(space, proj_norm){	
	} ; 
	~HPAdapt(){};
 

bool adapt_smooth(int* smooth_dof, int max_p);


};







#endif
