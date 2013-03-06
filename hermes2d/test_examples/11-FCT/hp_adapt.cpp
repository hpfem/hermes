#include "hp_adapt.h"

//int ref = no. of refinement steps
bool refine_elem(SpaceSharedPtr<double> space, Element* e, int ref)
{
  bool refined = true;
  int order = space->get_element_order(e->id);
  order = H2D_MAKE_QUAD_ORDER(1, 1);  //non_smooth elements -> polynomial degree 1
  if (e->active)
  {
    space->get_mesh()->refine_element_id(e->id);				
    for (int j = 0; j < 4; j++)
      space->set_element_order_internal(e->sons[j]->id, order);
  }
  if(ref>1) 
  {
    for (int j = 0; j < 4; j++)
      refined = refine_elem(space, e->sons[j], ref-1);
  }

  return refined;
}


bool HPAdapt::adapt_smooth(int* smooth_elem, int max_p)
{  
  if (this->num >1) 
    throw Hermes::Exceptions::Exception("adapt_smooth: Only for one space .");
  bool changed = false;
  SpaceSharedPtr<double> space = this->spaces[0];
  int order, v_ord, h_ord;
  int n_dof = space->get_num_dofs();
  Element* e;
  // Just set the verbose output to true so that the next calls to this->info() produce output.
  this->set_verbose_output(true);
  this->info("Elements are being refined in HPAdapt::adapt_smooth().");
    
  for_all_active_elements(e, space->get_mesh())
  {
    if(smooth_elem[e->id] == 1) //smooth => p increase	
    {						
      if(e->is_triangle()==true)
      {
        order = space->get_element_order(e->id); 
        order++;
        if(order>max_p) order = max_p;
        space->set_element_order_internal(e->id, order);
      }
      else if(e->is_quad()==true)
      {
        v_ord = H2D_GET_V_ORDER(space->get_element_order(e->id));
        h_ord = H2D_GET_H_ORDER(space->get_element_order(e->id));
        h_ord++; if(h_ord>max_p) h_ord =max_p;
        v_ord++; if(v_ord>max_p) v_ord =max_p;
        order = H2D_MAKE_QUAD_ORDER(h_ord, v_ord);
        space->set_element_order_internal(e->id, order);
        changed = true;
      }					
    }
    else if((smooth_elem[e->id] == 0)||(smooth_elem[e->id] == 2))	//non-smooth=>h refine					
      changed = refine_elem(space, e, 1);		
  }

  this->info("Assigning DOFs in HPAdapt::adapt_smooth().");
  space->assign_dofs();

  if(changed==false)
    return false;
  
  return true;
}

