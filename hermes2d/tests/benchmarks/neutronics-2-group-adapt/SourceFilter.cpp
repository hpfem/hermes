#include "hermes2d.h"
#include "SourceFilter.h"

void SourceFilter::precalculate(int order, int mask)
{
  if (mask & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
    error("SourceFilter not defined for derivatives.");

  Quad2D* quad = quads[cur_quad];
  int np = quad->get_num_points(order);
  Node* node = new_node(H2D_FN_VAL, np);

  sln[0]->set_quad_order(order, H2D_FN_VAL);
  sln[1]->set_quad_order(order, H2D_FN_VAL);

  scalar *val1 = sln[0]->get_fn_values();
  scalar *val2 = sln[1]->get_fn_values();
  update_refmap();
  
  double *x = refmap->get_phys_x(order);
	double *y = refmap->get_phys_y(order);
	
  for (int i = 0; i < np; i++)
  {
  	int m = this->get_material(x[i],y[i]);
  	node->values[0][0][i] = Sf[m][0] * val1[i] + Sf[m][1] * val2[i];
 	}
   
  // remove the old node and attach the new one
  replace_cur_node(node);
}

