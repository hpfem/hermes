// This file is part of Hermes2D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file api.h
\brief Main Hermes API
*/
#ifndef __HERMES_API_H_
#define __HERMES_API_H_

#include "compat.h"
#include <map>

namespace Hermes
{
  /// Enumeration of potential keys in the Api::parameters storage.
  enum HermesCommonApiParam
  {
    exceptionsPrintCallstack,
    numThreads,
    matrixSolverType,
    directMatrixSolverType,
    showInternalWarnings
  };

  /// API Class containing settings for the whole HermesCommon.
  class HERMES_API Api
  {
  public:
    Api();
    ~Api();

    /// Internal.
    /// Setter handler type.
    /// Serve for a custom reaction to some parameter settings.
    /// Such a handler must be registered in the map setter_handlers.
    typedef void (*SetterHandler)();

  protected:
    /// Parameter class, representing one parameter.
    /// Its identifier is a string identifier according to which, the instance is inserted into Api::parameters.
    class HERMES_API Parameter
    {
    public:
      /// Constructor.
      /// \param[in] default_val Default value, if the user does not specify his own.
      Parameter(int default_val);
      bool user_set;
      int user_val;
      int default_val;
    };

    /// The storage of parameters.
    /// This storage is not optimized for speed, but for comfort of users.
    /// There should not be any parameters, values of which are sought very often, because of the above reason.
    std::map<HermesCommonApiParam, Parameter*> parameters;

    /// Internal.
    /// Setter handlers.
    /// Purpose: when a parameter is set (such as the linear solver), some action might have to be taken to
    /// serve the event.
    std::map<std::pair<HermesCommonApiParam, int>, SetterHandler> setter_handlers;

    /// Internal.
    /// Change handlers.
    /// Used when the value of a particular parameter changes from the served one to another.
    /// Also used in destructor.
    std::map<std::pair<HermesCommonApiParam, int>, SetterHandler> change_handlers;

  public:
    int get_integral_param_value(HermesCommonApiParam);
    void set_integral_param_value(HermesCommonApiParam, int value);
  };

  /// Global instance used inside Hermes which is also accessible to users.
  HERMES_COMMON_API extern Hermes::Api HermesCommonApi;
}
#endif