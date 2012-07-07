#include "global.h"
#include "refinement_type.h"

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