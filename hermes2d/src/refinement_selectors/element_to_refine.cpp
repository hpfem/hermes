#include "element_to_refine.h"
#include "refinement_selectors/candidates.h"

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
      copy_orders(this->refinement_polynomial_order, orig.refinement_polynomial_order);
      for(int i = 0; i < 4; i++)
        copy_orders(this->best_refinement_polynomial_order_type[i], orig.best_refinement_polynomial_order_type[i]);
      copy_errors(errors, orig.errors);
    };

    int ElementToRefine::get_num_sons() const
    {
      return get_refin_sons(split);
    };

    void ElementToRefine::copy_orders(int* dest, const int* src)
    {
      memcpy(dest, src, sizeof(int) * H2D_MAX_ELEMENT_SONS);
    }

    void ElementToRefine::copy_errors(double* dest, const double* src)
    {
      memcpy(dest, src, sizeof(double) * H2D_MAX_ELEMENT_SONS);
    }

    ElementToRefine& ElementToRefine::operator=(const ElementToRefine& orig)
    {
      id = orig.id;
      comp = orig.comp;
      split = orig.split;
      copy_orders(this->refinement_polynomial_order, orig.refinement_polynomial_order);
      for(int i = 0; i < 5; i++)
        copy_orders(this->best_refinement_polynomial_order_type[i], orig.best_refinement_polynomial_order_type[i]);
      return *this;
    }
  }
}