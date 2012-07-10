// This file is part of HermesCommon
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
/*! \file dpinterface.cpp
\brief Interface for DiscreteProblem required by NoxProblemInterface.
*/
#include "dp_interface.h"

/// \todo If there are really no bodies, delete this file.
namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    DiscreteProblemInterface<Scalar>::DiscreteProblemInterface() : globalIntegrationOrderSet(false), globalIntegrationOrder(0)
    {}

    template<typename Scalar>
    void DiscreteProblemInterface<Scalar>::setGlobalIntegrationOrder(unsigned int order)
    {
      this->globalIntegrationOrder = order;
      this->globalIntegrationOrderSet = true;
    }

    template class HERMES_API DiscreteProblemInterface<double>;
    template class HERMES_API DiscreteProblemInterface<std::complex<double> >;
  }
}