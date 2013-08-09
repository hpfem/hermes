// This file is part of Hermes2D.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes; if not, see <http://www.gnu.prg/licenses/>.
#include "mixins2d.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Mixins
    {
      void StateQueryable::check() const
      {
        if(!this->isOkay())
        {
          std::stringstream ss;
          ss << "The instance of " << this->getClassName() << " is not OK.";
          throw Hermes::Exceptions::Exception(ss.str().c_str());
        }
      }

      XMLParsing::XMLParsing() : validate(false)
      {
      }

      void XMLParsing::set_validation(bool to_set)
      {
        this->validate = to_set;
      }


      Parallel::Parallel() : num_threads_used(HermesCommonApi.get_integral_param_value(numThreads))
      {
      }
    }
  }
}
