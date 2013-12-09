#include "refinement_selectors/candidates.h"

namespace Hermes
{
  namespace Hermes2D
  {
    HERMES_API bool is_refin_aniso(const int refin_type)
    {
      if(refin_type == H2D_REFINEMENT_ANISO_H || refin_type == H2D_REFINEMENT_ANISO_V)
        return true;
      else
        return false;
    }

    HERMES_API int get_refin_sons(const int refin_type)
    {
      switch(refin_type)
      {
      case H2D_REFINEMENT_P: return 1; break;
      case H2D_REFINEMENT_H: return 4; break;
      case H2D_REFINEMENT_ANISO_H:
      case H2D_REFINEMENT_ANISO_V: return 2; break;
      default: throw Hermes::Exceptions::Exception("Invalid refinement type %d", (int)refin_type); return -1;
      }
    }

    HERMES_API const std::string get_refin_str(const int refin_type)
    {
      switch(refin_type)
      {
      case H2D_REFINEMENT_P: return "P"; break;
      case H2D_REFINEMENT_H: return "H"; break;
      case H2D_REFINEMENT_ANISO_H: return "AnisoH"; break;
      case H2D_REFINEMENT_ANISO_V: return "AnisoV"; break;
      default:
        std::stringstream str;
        str << "Unknown(" << refin_type << ")";
        return str.str();
      }
    }

    namespace RefinementSelectors
    {
      HERMES_API const char* get_cand_list_str(const CandList cand_list)
      {
        switch(cand_list)
        {
        case H2D_NONE: return "Custom";
        case H2D_P_ISO: return "P_ISO";
        case H2D_P_ANISO: return "P_ANISO";
        case H2D_H_ISO: return "H_ISO";
        case H2D_H_ANISO: return "H_ANISO";
        case H2D_HP_ISO: return "HP_ISO";
        case H2D_HP_ANISO_H: return "HP_ANISO_H";
        case H2D_HP_ANISO_P: return "HP_ANISO_P";
        case H2D_HP_ANISO: return "HP_ANISO";
        default:
          throw Hermes::Exceptions::Exception("Invalid adapt type %d.", cand_list);
          return NULL;
        }
      }

      HERMES_API bool is_hp(const CandList cand_list)
      {
        switch(cand_list)
        {
        case H2D_P_ISO:
        case H2D_P_ANISO:
        case H2D_H_ISO:
        case H2D_H_ANISO: return false; break;
        case H2D_NONE:
        case H2D_HP_ISO:
        case H2D_HP_ANISO_H:
        case H2D_HP_ANISO_P:
        case H2D_HP_ANISO: return true; break;
        default: throw Hermes::Exceptions::Exception("Invalid adapt type %d.", cand_list); return false;
        }
      }

       HERMES_API bool is_p(const CandList cand_list)
      {
        switch(cand_list)
        {
        case H2D_H_ISO:
        case H2D_H_ANISO: return false; break;
        case H2D_NONE:
        case H2D_P_ISO:
        case H2D_P_ANISO:
        case H2D_HP_ISO:
        case H2D_HP_ANISO_H:
        case H2D_HP_ANISO_P:
        case H2D_HP_ANISO: return true; break;
        default: throw Hermes::Exceptions::Exception("Invalid adapt type %d.", cand_list); return false;
        }
      }

      HERMES_API bool is_p_aniso(const CandList cand_list)
      {
        switch(cand_list)
        {
        case H2D_NONE: return false;
        case H2D_P_ISO: return false;
        case H2D_P_ANISO: return true;
        case H2D_H_ISO: return false;
        case H2D_H_ANISO: return false;
        case H2D_HP_ISO: return false;
        case H2D_HP_ANISO_H: return false;
        case H2D_HP_ANISO_P: return true;
        case H2D_HP_ANISO: return true;
        default: throw Hermes::Exceptions::Exception("Invalid adapt type %d.", cand_list); return false;
        }
      }

      Cand::Cand(const int split, const int order_elems[H2D_MAX_ELEMENT_SONS])
        : dofs(-1), split(split), score(0) {
          p[0] = order_elems[0];
          p[1] = order_elems[1];
          p[2] = order_elems[2];
          p[3] = order_elems[3];
      };

      Cand::Cand(const int split, const int order_elem0, const int order_elem1, const int order_elem2, const int order_elem3)
        : dofs(-1), split(split), score(0) {
          p[0] = order_elem0;
          p[1] = order_elem1;
          p[2] = order_elem2;
          p[3] = order_elem3;
      };

      int Cand::get_num_elems() const {
        switch (split) {
        case H2D_REFINEMENT_H: return 4;
        case H2D_REFINEMENT_P: return 1;
        case H2D_REFINEMENT_ANISO_H:
        case H2D_REFINEMENT_ANISO_V:
          return 2;
        default:
          throw Hermes::Exceptions::Exception("Invalid refinement type %d.", split);
          return -1;
          break;
        }
      }
    };
  }
}