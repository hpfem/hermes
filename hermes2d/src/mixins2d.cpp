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
      template<typename Scalar>
      const Space<Scalar>* SettableSpaces<Scalar>::get_space(int n) const
      {
        return this->get_spaces()[n];
      }


      template<typename Scalar>
      MatrixRhsOutput<Scalar>::MatrixRhsOutput() : outputMatrixOn(false), outputMatrixIterations(-1), matrixFilename("Matrix_"), 
        matrixVarname("A"), matrixFormat(Hermes::Algebra::DF_MATLAB_SPARSE), outputRhsOn(false), outputRhsIterations(-1),
        RhsFilename("Rhs_"), RhsVarname("b"), RhsFormat(Hermes::Algebra::DF_MATLAB_SPARSE)
      {
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::outputMatrix(int firstIterations)
      {
        outputMatrixOn = true;
        this->outputMatrixIterations = firstIterations;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::setMatrixFilename(std::string name)
      {
        this->matrixFilename = name;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::setMatrixVarname(std::string name)
      {
        this->matrixVarname = name;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::setMatrixEMatrixDumpFormat(EMatrixDumpFormat format)
      {
        this->matrixFormat = format;
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::outputRhs(int firstIterations)
      {
        this->outputRhsOn = true;
        this->outputRhsIterations = firstIterations;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::setRhsFilename(std::string name)
      {
        this->RhsFilename = name;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::setRhsVarname(std::string name)
      {
        this->RhsVarname = name;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::setRhsEMatrixDumpFormat(EMatrixDumpFormat format)
      {
        this->RhsFormat = format;
      }

      template HERMES_API class SettableSpaces<double>;
      template HERMES_API class SettableSpaces<std::complex<double> >;
      template HERMES_API class MatrixRhsOutput<double>;
      template HERMES_API class MatrixRhsOutput<std::complex<double> >;
    }
  }
}