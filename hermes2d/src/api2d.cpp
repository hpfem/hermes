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
namespace Hermes
{
  namespace Hermes2D
  {

    template<typename T>
    PointerCalculator<T>::PointerCalculator() : count(0)
    {
    }

    template<typename T>
    unsigned int PointerCalculator<T>::getNumber() const
    {
      return this->count;
    }

    template<typename T>
    void PointerCalculator<T>::operator+(unsigned int increaseBy)
    {
      this->count += increaseBy;
    }

    template<typename T>
    void PointerCalculator<T>::operator++()
    {
      this->count++;
    }

    template<typename T>
    void PointerCalculator<T>::operator-(unsigned int decreaseBy)
    {
			if(this->count < decreaseBy)
				Hermes::Mixins::Loggable::Static::warn("PointerCalculator: it was detected that with this decrement, the count would be < 0. Probably some increment was accidentally omitted.");
      this->count -= decreaseBy;
    }

    template<typename T>
    void PointerCalculator<T>::operator--()
    {
			if(this->count < 1)
				Hermes::Mixins::Loggable::Static::warn("PointerCalculator: it was detected that with this decrement, the count would be < 0. Probably some increment was accidentally omitted.");
      this->count--;
    }


    Api2D::Parameter::Parameter(int default_val)
    {
      this->default_val = default_val;
      this->user_set = false;
    }

    Api2D::Api2D()
    {
      signal(SIGABRT, CallStack::dump);
      signal(SIGFPE, CallStack::dump);
      signal(SIGILL, CallStack::dump);
      signal(SIGSEGV, CallStack::dump);
      signal(SIGTERM, CallStack::dump);

      this->parameters.insert(std::pair<Hermes2DApiParam, Parameter*> (Hermes::Hermes2D::numThreads,new Parameter(NUM_THREADS)));
      this->parameters.insert(std::pair<Hermes2DApiParam, Parameter*> (Hermes::Hermes2D::secondDerivatives,new Parameter(0)));
    }

    Api2D::~Api2D()
    {
      this->parameters.clear();
    }

    int Api2D::get_param_value(Hermes2DApiParam param)
    {
      if(this->parameters.find(param) == parameters.end())
        throw Hermes::Exceptions::Exception("Wrong Hermes::Api parameter name:%i", param);
      if(this->parameters.find(param)->second->user_set)
        return this->parameters.find(param)->second->user_val;
      else
        return this->parameters.find(param)->second->default_val;
    }

    void Api2D::set_param_value(Hermes2DApiParam param, int value)
    {
      if(this->parameters.find(param) == parameters.end())
        throw Hermes::Exceptions::Exception("Wrong Hermes::Api parameter name:%i", param);
      this->parameters.find(param)->second->user_set = true;
      this->parameters.find(param)->second->user_val = value;
    }

    unsigned int Api2D::getNumberMeshPointers() const
    {
      return this->meshPointerCalculator.getNumber();
    }

    unsigned int Api2D::getNumberSpacePointers() const
    {
      return this->realSpacePointerCalculator.getNumber() + this->complexSpacePointerCalculator.getNumber();
    }

    unsigned int Api2D::getNumberRealSpacePointers() const
    {
      return this->realSpacePointerCalculator.getNumber();
    }

    unsigned int Api2D::getNumberComplexSpacePointers() const
    {
      return this->complexSpacePointerCalculator.getNumber();
    }

    unsigned int Api2D::getNumberSolutionPointers() const
    {
      return this->realSolutionPointerCalculator.getNumber() + this->complexSolutionPointerCalculator.getNumber();
    }

    unsigned int Api2D::getNumberRealSolutionPointers() const
    {
      return this->realSolutionPointerCalculator.getNumber();
    }

    unsigned int Api2D::getNumberComplexSolutionPointers() const
    {
      return this->complexSolutionPointerCalculator.getNumber();
    }

		unsigned int Api2D::getNumberMeshData() const
		{
			return this->meshDataPointerCalculator.getNumber();
		}

		unsigned int Api2D::getNumberSpaceData() const
		{
			return this->realSpaceDataPointerCalculator.getNumber() + this->complexSpaceDataPointerCalculator.getNumber();
		}

		unsigned int Api2D::getNumberRealSpaceData() const
		{
			return this->realSpaceDataPointerCalculator.getNumber();
		}

		unsigned int Api2D::getNumberComplexSpaceData() const
		{
			return this->complexSpaceDataPointerCalculator.getNumber();
		}

		unsigned int Api2D::getNumberSolutionData() const
		{
			return this->realSolutionDataPointerCalculator.getNumber() + this->complexSolutionDataPointerCalculator.getNumber();
		}

		unsigned int Api2D::getNumberRealSolutionData() const
		{
			return this->realSolutionDataPointerCalculator.getNumber();
		}

		unsigned int Api2D::getNumberComplexSolutionData() const
		{
			return this->complexSolutionDataPointerCalculator.getNumber();
		}

    Hermes::Hermes2D::Api2D Hermes2DApi;
    template class HERMES_API PointerCalculator<Mesh>;
    template class HERMES_API PointerCalculator<Space<double> >;
    template class HERMES_API PointerCalculator<Space<std::complex<double> > >;
    template class HERMES_API PointerCalculator<Solution<double> >;
    template class HERMES_API PointerCalculator<Solution<std::complex<double> > >;
  }
}