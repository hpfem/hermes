#ifndef __H2D_MIXINS_H
#define __H2D_MIXINS_H
#include "api2d.h"
#include "global.h"
namespace Hermes
{
  namespace Hermes2D
  {
    /** \defgroup g_mixins2d Mixins
    *  \brief Mixins are utility classes used for all kinds of other classes.
    *
    *  Mixin classes provide a single piece of functionality.
    *
    */

    namespace Mixins
    {

      /// \ingroup g_mixins2d
      /// Any XML parsing class should inherit from this mixin.
      /// It serves various purposes, first of which is disabling / re-enabling of validation
      /// against the schema referenced in a file being loaded.
      class HERMES_API XMLParsing
      {
      public:
        /// Constructor.
        XMLParsing();

        /// Set to validate / not to validate.
        void set_validation(bool to_set);

      protected:
        /// Internal.
        bool validate;
      };

      /// \brief Class utilizes parallel calculation
      class HERMES_API Parallel
      {
      protected:
        Parallel();
      protected:
        int num_threads_used;
        std::string exceptionMessageCaughtInParallelBlock;
      };
    }
  }
}
#endif
