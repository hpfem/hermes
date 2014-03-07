#include "element_to_refine.h"

namespace Hermes
{
  namespace Hermes2D
  {
    bool HERMES_API is_refin_aniso(const RefinementType refin_type)
    {
      if (refin_type == H2D_REFINEMENT_H_ANISO_H || refin_type == H2D_REFINEMENT_H_ANISO_V)
        return true;
      else
        return false;
    }

    int HERMES_API get_refin_sons(const RefinementType refin_type)
    {
      switch (refin_type)
      {
      case H2D_REFINEMENT_P: return 1; break;
      case H2D_REFINEMENT_H: return 4; break;
      case H2D_REFINEMENT_H_ANISO_H:
      case H2D_REFINEMENT_H_ANISO_V: return 2; break;
      default: throw Hermes::Exceptions::Exception("Invalid refinement type %d", (int)refin_type); return -1;
      }
    }

    const HERMES_API std::string get_refin_str(const RefinementType refin_type)
    {
      switch (refin_type)
      {
      case H2D_REFINEMENT_P: return "P"; break;
      case H2D_REFINEMENT_H: return "H"; break;
      case H2D_REFINEMENT_H_ANISO_H: return "AnisoH"; break;
      case H2D_REFINEMENT_H_ANISO_V: return "AnisoV"; break;
      default:
        std::stringstream str;
        str << "Unknown(" << refin_type << ")";
        return str.str();
      }
    }

    ElementToRefine::ElementToRefine() : valid(false)
    {
    };

    ElementToRefine::ElementToRefine(int id, unsigned short comp) : id(id), comp(comp), split(H2D_REFINEMENT_H), valid(false)
    {
    };

    unsigned short ElementToRefine::get_num_sons() const
    {
      return get_refin_sons(split);
    };

    void ElementToRefine::copy_orders(unsigned short* dest, const unsigned short* src)
    {
      memcpy(dest, src, sizeof(unsigned short)* H2D_MAX_ELEMENT_SONS);
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
      valid = orig.valid;
      copy_orders(this->refinement_polynomial_order, orig.refinement_polynomial_order);
      for(int i = 0; i < 5; i++)
        copy_orders(this->best_refinement_polynomial_order_type[i], orig.best_refinement_polynomial_order_type[i]);
      return *this;
    }
  }
}