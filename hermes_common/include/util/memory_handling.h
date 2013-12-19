// This file is part of HermesCommon.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

/*! \file memory_handling.h
\brief File containing common definitions, and basic global enums etc. for HermesCommon.
*/
#ifndef __HERMES_COMMON_MEMORY_HANDLING_H
#define __HERMES_COMMON_MEMORY_HANDLING_H

#include "exceptions.h"

namespace Hermes
{

  template<typename Caller, typename ArrayItem>
  ArrayItem* calloc_with_check(int size, Caller* const caller)
  {
    ArrayItem* new_array = (ArrayItem*)calloc(size, sizeof(ArrayItem));
    if (new_array)
      return new_array;
    else
    {
      caller->free();
      throw Hermes::Exceptions::Exception("Hermes::realloc_with_check() failed to reallocate.", size * sizeof(ArrayItem));
      return nullptr;
    }
  }

  template<typename Caller, typename ArrayItem>
  ArrayItem* malloc_with_check(int size, Caller* const caller)
  {
    ArrayItem* new_array = (ArrayItem*)malloc(size * sizeof(ArrayItem));
    if (new_array)
      return new_array;
    else
    {
      caller->free();
      throw Hermes::Exceptions::Exception("Hermes::realloc_with_check() failed to reallocate.", size * sizeof(ArrayItem));
      return nullptr;
    }
  }

  template<typename Caller, typename ArrayItem>
  ArrayItem* realloc_with_check(ArrayItem* original_array, int new_size, Caller* const caller)
  {
    ArrayItem* new_array = (ArrayItem*)realloc(original_array, new_size * sizeof(ArrayItem));
    if (new_array)
      return original_array = new_array;
    else
    {
      caller->free();
      throw Hermes::Exceptions::Exception("Hermes::realloc_with_check() failed to reallocate.", new_size * sizeof(ArrayItem));
      return nullptr;
    }
  }
}
#endif
