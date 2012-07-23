#ifndef __H2D_MIXINS_H
#define __H2D_MIXINS_H
#include "space/space.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace Mixins
    {
      template<typename Scalar>
      class HERMES_API SettableSpaces
      {
      public:
        /// Sets new spaces for the instance.
        virtual void set_spaces(Hermes::vector<const Space<Scalar>*> spaces) = 0;
        virtual void set_space(const Space<Scalar>* space) = 0;
        /// Get all spaces as a Hermes::vector.
        virtual Hermes::vector<const Space<Scalar>*> get_spaces() const = 0;
        virtual const Space<Scalar>* get_space(int n) const;
      };
    }
  }
}
#endif