#include "hermes2d_common_defs.h"
#include "range.h"
#include "element_to_refine.h"

namespace Hermes
{
  namespace Hermes2D
  {
    ElementToRefine::ElementToRefine() : id(-1), comp(-1) 
    {
    };

    ElementToRefine::ElementToRefine(int id, int comp) : id(id), comp(comp), split(H2D_REFINEMENT_H) 
    {
    };

    ElementToRefine::ElementToRefine(const ElementToRefine &orig) : id(orig.id), comp(orig.comp), split(orig.split) 
    {
      copy_orders(p, orig.p);
      copy_orders(q, orig.q);
    };

    int ElementToRefine::get_num_sons() const
    {
      return get_refin_sons(split); 
    };

    void ElementToRefine::copy_orders(int* dest, const int* src)
    {
      memcpy(dest, src, sizeof(int) * H2D_MAX_ELEMENT_SONS);
    }

    ElementToRefine& ElementToRefine::operator=(const ElementToRefine& orig)
    {
      id = orig.id;
      comp = orig.comp;
      split = orig.split;
      copy_orders(p, orig.p);
      copy_orders(q, orig.q);
      return *this;
    }

    HERMES_API std::ostream& operator<<(std::ostream& stream, const ElementToRefine& elem_ref)
    {
      stream << "id:" << elem_ref.id << ";comp:" << elem_ref.comp << "; split:" << get_refin_str(elem_ref.split) << "; orders:[";
      int num_sons = elem_ref.get_num_sons();
      for(int i = 0; i < num_sons; i++)
      {
        if (i > 0)
          stream << " ";
        stream << Global<double>::get_quad_order_str(elem_ref.p[i]);
      }
      stream << "]";
      return stream;
    }
  }
}