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
    matrixSolverType
  };

  /// API Class containing settings for the whole HermesCommon.
  class HERMES_API Api
  {
  public:
    Api();
    ~Api();
  protected:
    /// Parameter class, representing one parameter.
    /// Its identifier is a string identifier according to which, the instance is inserted into Api::parameters.
    class HERMES_API Parameter
    {
    public:
      /// Constructor.
      /// \param[in] defaultVal Default value, if the user does not specify his own.
      Parameter(int defaultVal);
      bool userSet;
      int userVal;
      int defaultVal;
    };
    /// The storage of parameters.
    /// This storage is not optimized for speed, but for comfort of users.
    /// There should not be any parameters, values of which are sought very often, because of the above reason.
    std::map<HermesCommonApiParam, Parameter*> parameters;

  public:
    int getParamValue(HermesCommonApiParam);
    void setParamValue(HermesCommonApiParam, int value);
  };

  /// Global instance used inside Hermes which is also accessible to users.
  extern HERMES_API Hermes::Api HermesCommonApi;
}
#endif