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

#include "callstack.h"
#include "api.h"
#include "common.h"
#include "exceptions.h"
#include "api2d.h"
#include <xercesc/util/PlatformUtils.hpp>

using namespace xercesc;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename T>
    Api2D::Parameter<T>::Parameter(T default_val)
    {
      this->default_val = default_val;
      this->user_set = false;
    }

    Api2D::Api2D()
    {
      int* asf = (int*)malloc(1000);
      signal(SIGABRT, CallStack::dump);
      signal(SIGFPE, CallStack::dump);
      signal(SIGILL, CallStack::dump);
      signal(SIGSEGV, CallStack::dump);
      signal(SIGTERM, CallStack::dump);
      
      // Xerces initialization - for better performance.
      XMLPlatformUtils::Initialize();   

      this->integral_parameters.insert(std::pair<Hermes2DApiParam, Parameter<int>*> (Hermes::Hermes2D::numThreads,new Parameter<int>(NUM_THREADS)));
      this->text_parameters.insert(std::pair<Hermes2DApiParam, Parameter<std::string>*> (Hermes::Hermes2D::xmlSchemasDirPath,new Parameter<std::string>(*(new std::string(H2D_XML_SCHEMAS_DIRECTORY)))));
      std::stringstream ss;
      ss << H2D_PRECALCULATED_FORMS_DIRECTORY;
      if(ss.str().at(ss.str().length() - 1) == '\\' || ss.str().at(ss.str().length() - 1) == '/')
        this->text_parameters.insert(std::pair<Hermes2DApiParam, Parameter<std::string>*> (Hermes::Hermes2D::precalculatedFormsDirPath,new Parameter<std::string>(*(new std::string(H2D_PRECALCULATED_FORMS_DIRECTORY)))));
      else
      {
        ss << '/';
        this->text_parameters.insert(std::pair<Hermes2DApiParam, Parameter<std::string>*> (Hermes::Hermes2D::precalculatedFormsDirPath,new Parameter<std::string>(*(new std::string(ss.str())))));
      }

      XMLPlatformUtils::Terminate();
    }

    Api2D::~Api2D()
    {
      for(std::map<Hermes2DApiParam, Parameter<std::string>*>::const_iterator it = this->text_parameters.begin(); it != this->text_parameters.end(); ++it)
        delete it->second;

      for(std::map<Hermes2DApiParam, Parameter<int>*>::const_iterator it = this->integral_parameters.begin(); it != this->integral_parameters.end(); ++it)
        delete it->second;
    }

    int Api2D::get_integral_param_value(Hermes2DApiParam param)
    {
      if(this->integral_parameters.find(param) == integral_parameters.end())
        throw Hermes::Exceptions::Exception("Wrong Hermes::Api parameter name:%i", param);
      if(this->integral_parameters.find(param)->second->user_set)
        return this->integral_parameters.find(param)->second->user_val;
      else
        return this->integral_parameters.find(param)->second->default_val;
    }

    void Api2D::set_integral_param_value(Hermes2DApiParam param, int value)
    {
      if(this->integral_parameters.find(param) == integral_parameters.end())
        throw Hermes::Exceptions::Exception("Wrong Hermes::Api parameter name:%i", param);
      this->integral_parameters.find(param)->second->user_set = true;
      this->integral_parameters.find(param)->second->user_val = value;
    }

    std::string Api2D::get_text_param_value(Hermes2DApiParam param)
    {
      if(this->text_parameters.find(param) == text_parameters.end())
        throw Hermes::Exceptions::Exception("Wrong Hermes::Api parameter name:%i", param);
      if(this->text_parameters.find(param)->second->user_set)
        return this->text_parameters.find(param)->second->user_val;
      else
        return this->text_parameters.find(param)->second->default_val;
    }

    void Api2D::set_text_param_value(Hermes2DApiParam param, std::string value)
    {
      if(this->text_parameters.find(param) == text_parameters.end())
        throw Hermes::Exceptions::Exception("Wrong Hermes::Api parameter name:%i", param);
      this->text_parameters.find(param)->second->user_set = true;
      this->text_parameters.find(param)->second->user_val = value;
    }

    Hermes::Hermes2D::Api2D HERMES_API Hermes2DApi;
  }
}