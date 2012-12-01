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
#ifndef __HERMES_API_2D_H_
#define __HERMES_API_2D_H_
// This is here because of an operator-related error in gcc
#pragma GCC diagnostic warning "-fpermissive"

#include "compat.h"
#include "hermes_common.h"
#include "config.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class Mesh;
    template<typename Scalar> class Space;
    template<typename Scalar> class Solution;

    /// Enumeration of potential keys in the Api2D::parameters storage.
    enum Hermes2DApiParam
    {
      numThreads,
      secondDerivatives,
			xmlSchemasDirPath
    };


    /// Class for calculating pointers of instance T.
    template<typename T>
    class HERMES_API PointerCalculator
    {
    public:
      PointerCalculator();
      unsigned int getNumber() const;
    private:
      void operator+(unsigned int increaseBy);
      void operator++();
      void operator-(unsigned int decreaseBy);
      void operator--();
      unsigned int count;
			friend class Mesh;
      friend class MeshReaderH2DXML;
      template<typename T1> friend class Space;
      template<typename T1> friend class Solution;
    };

    /// API Class containing settings for the whole Hermes2D.
    class HERMES_API Api2D
    {
    public:
      Api2D();
      ~Api2D();
    protected:
      /// Parameter class, representing one parameter.
      /// Its identifier is a string identifier according to which, the instance is inserted into Api2D::parameters.
			template<typename T>
      class HERMES_API Parameter
      {
      public:
        /// Constructor.
        /// \param[in] default_val Default value, if the user does not specify his own.
				Parameter(T default_val);
        bool user_set;
        T user_val;
        T default_val;
      };

      /// The storage of parameters.
      /// This storage is not optimized for speed, but for comfort of users.
      /// There should not be any parameters, values of which are sought very often, because of the above reason.

			std::map<Hermes2DApiParam, Parameter<int>*> integral_parameters;
      std::map<Hermes2DApiParam, Parameter<std::string>*> text_parameters;
    public:
			int get_integral_param_value(Hermes2DApiParam);
      std::string get_text_param_value(Hermes2DApiParam);

			void set_integral_param_value(Hermes2DApiParam, int value);
			void set_text_param_value(Hermes2DApiParam, std::string value);

			/// Returns the number of Mesh pointers that exist.
			/// Defined as the difference between constructor calls and destructor calls.
      unsigned int getNumberMeshPointers() const;
      
			/// Returns the number of Space pointers that exist.
			/// Defined as the difference between constructor calls and destructor calls.
			unsigned int getNumberSpacePointers() const;

			/// Returns the number of Space<double> pointers that exist.
			/// Defined as the difference between constructor calls and destructor calls.
      unsigned int getNumberRealSpacePointers() const;

			/// Returns the number of Space<std::complex<double> > pointers that exist.
			/// Defined as the difference between constructor calls and destructor calls.
      unsigned int getNumberComplexSpacePointers() const;

			/// Returns the number of Solution pointers that exist.
			/// Defined as the difference between constructor calls and destructor calls.
      unsigned int getNumberSolutionPointers() const;

			/// Returns the number of Solution<double> pointers that exist.
			/// Defined as the difference between constructor calls and destructor calls.
      unsigned int getNumberRealSolutionPointers() const;

			/// Returns the number of Solution<std::complex<double> > pointers that exist.
			/// Defined as the difference between constructor calls and destructor calls.
      unsigned int getNumberComplexSolutionPointers() const;

			/// Returns the number of Mesh data that exist.
			/// Defined as the difference between init() calls and free() calls
			unsigned int getNumberMeshData() const;

			/// Returns the number of Space data that exist.
			/// Defined as the difference between init() calls and free() calls
			unsigned int getNumberSpaceData() const;

			/// Returns the number of Space<double> data that exist.
			/// Defined as the difference between init() calls and free() calls
			unsigned int getNumberRealSpaceData() const;

			/// Returns the number of Space<std::complex<double> > data that exist.
			/// Defined as the difference between init() calls and free() calls
			unsigned int getNumberComplexSpaceData() const;

			/// Returns the number of Solution data that exist.
			/// Defined as the difference between init() calls and free() calls
			unsigned int getNumberSolutionData() const;

			/// Returns the number of Solution<double> data that exist.
			/// Defined as the difference between init() calls and free() calls
			unsigned int getNumberRealSolutionData() const;

			/// Returns the number of Solution<std::complex<double> > data that exist.
			/// Defined as the difference between init() calls and free() calls
			unsigned int getNumberComplexSolutionData() const;

    private:
			// Storage for handling constructor / destructor calls.
			PointerCalculator<Mesh> meshPointerCalculator;
      PointerCalculator<Space<double> > realSpacePointerCalculator;
      PointerCalculator<Space<std::complex<double> > > complexSpacePointerCalculator;
      PointerCalculator<Solution<double> > realSolutionPointerCalculator;
      PointerCalculator<Solution<std::complex<double> > > complexSolutionPointerCalculator;

			// Storage for handling calls to init() / free().
			PointerCalculator<Mesh> meshDataPointerCalculator;
			PointerCalculator<Space<double> > realSpaceDataPointerCalculator;
			PointerCalculator<Space<std::complex<double> > > complexSpaceDataPointerCalculator;
			PointerCalculator<Solution<double> > realSolutionDataPointerCalculator;
			PointerCalculator<Solution<std::complex<double> > > complexSolutionDataPointerCalculator;

			friend class Mesh;
			friend class MeshReaderH2DXML;
      template<typename T1> friend class Space;
      template<typename T1> friend class Solution;
    };

    /// Global instance used inside Hermes which is also accessible to users.
    extern HERMES_API Hermes::Hermes2D::Api2D Hermes2DApi;
  }
}
#endif