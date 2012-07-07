#include "global.h"
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
  }
}