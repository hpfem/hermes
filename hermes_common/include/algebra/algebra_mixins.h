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
/*! \file algebra_mixins.h
\brief Mixins classes for algebraic purposes.
*/
#ifndef __HERMES_COMMON_ALGEBRA_MIXINS_H
#define __HERMES_COMMON_ALGEBRA_MIXINS_H

#include "common.h"
#include "util/compat.h"
#include "algebra_utilities.h"
#include "exceptions.h"

namespace Hermes
{
  namespace Algebra
  {
    template<typename Scalar> class SparseMatrix;
    template<typename Scalar> class Vector;

    namespace Mixins
    {
      /// \ingroup g_mixins2d
      /// Mixin that interfaces linear algebra structures output on higher levels (for solvers).
      template<typename Scalar>
      class HERMES_API MatrixRhsOutput
      {
      public:
        /// Constructor.
        /// Sets defaults (see individual set methods for values of those).
        MatrixRhsOutput();

        /// Processes the matrix.
        void process_matrix_output(Hermes::Algebra::SparseMatrix<Scalar>* matrix, int iteration);
        void process_matrix_output(Hermes::Algebra::SparseMatrix<Scalar>* matrix);

        /// Processes the matrix.
        void process_vector_output(Hermes::Algebra::Vector<Scalar>* rhs, int iteration);
        void process_vector_output(Hermes::Algebra::Vector<Scalar>* rhs);

        /// Sets this instance to output the matrix in several first iterations.
        /// \param[in] only_last_iteration If true, only the last iteration is outputted, and the next parameter is ignored.
        /// \param[in] firstIterations Only during so many first iterations. Default: -1 meaning, that during all iterations, the matrix will be saved.
        void output_matrix(bool only_last_iteration = true, int firstIterations = -1);
        /// Sets this instance to output matrix entries even though they are zero or not.
        void set_print_zero_matrix_entries(bool to_set);
        /// Sets filename for the matrix
        /// Default: Matrix_'iteration number' with the ".m" extension in the case of matlab format.
        /// \param[in] name sets the main part of the name, i.e. replacement for "Matrix_" in the default name.
        void set_matrix_filename(std::string name);
        /// Sets varname for the matrix
        /// Default: "A".
        void set_matrix_varname(std::string name);
        /// Sets varname for the matrix
        /// Default: "DF_MATLAB_SPARSE - matlab file".
        void set_matrix_export_format(Hermes::Algebra::MatrixExportFormat format);
        /// Sets number format for the matrix output.
        /// Default: "%lf".
        void set_matrix_number_format(char* number_format);

        /// Sets this instance to output the rhs in several first iterations.
        /// \param[in] only_last_iteration If true, only the last iteration is outputted, and the next parameter is ignored.
        /// \param[in] firstIterations Only during so many first iterations. Default: -1 meaning, that during all iterations, the rhs will be saved.
        void output_rhs(bool only_last_iteration = true, int firstIterations = -1);
        /// Sets filename for the rhs
        /// Default: Rhs_'iteration number' with the ".m" extension in the case of matlab format.
        /// \param[in] name sets the main part of the name, i.e. replacement for "Rhs_" in the default name.
        void set_rhs_filename(std::string name);
        /// Sets varname for the rhs
        /// Default: "b".
        void set_rhs_varname(std::string name);
        /// Sets varname for the rhs
        /// Default: "DF_MATLAB_SPARSE - matlab file".
        void set_rhs_export_format(Hermes::Algebra::MatrixExportFormat format);
        /// Sets number format for the vector output.
        /// Default: "%lf".
        void set_rhs_number_format(char* number_format);

      protected:
        bool print_matrix_zero_values;
        bool output_matrixOn;
        bool only_lastMatrixIteration;
        int output_matrixIterations;
        std::string matrixFilename;
        std::string matrixVarname;
        Hermes::Algebra::MatrixExportFormat matrixFormat;
        char* matrix_number_format;

        bool output_rhsOn;
        bool only_lastRhsIteration;
        int output_rhsIterations;
        std::string RhsFilename;
        std::string RhsVarname;
        Hermes::Algebra::MatrixExportFormat RhsFormat;
        char* rhs_number_format;
      };

      /// \ingroup g_mixins2d
      /// Mixin that interfaces basic linear algebra structures output.
      template<typename Scalar>
      class HERMES_API MatrixRhsImportExport
      {
      public:
        /// writing matrix.
        /// @param[in] filename obvious
        /// @param[in] var_name name of variable (will be written to output file)
        /// @param[in] fmt output file format
        /// @param[in] number_format specifies the number format where possible
        virtual void export_to_file(const char* filename, const char* var_name, Algebra::MatrixExportFormat fmt, char* number_format = "%lf") = 0;
        void export_to_file(std::string filename, const char* var_name, Algebra::MatrixExportFormat fmt, char* number_format = "%lf");
        void export_to_file(std::string filename, std::string var_name, Algebra::MatrixExportFormat fmt, char* number_format = "%lf");

        /// reading matrix
        /// @param[in] filename obvious
        /// @param[in] var_name name of variable (will be searched for to output file)
        /// @param[in] fmt input file format
        virtual void import_from_file(const char* filename, const char* var_name, Algebra::MatrixExportFormat fmt) { throw Hermes::Exceptions::MethodNotOverridenException("MatrixRhsImportExport<Scalar>::import_from_file"); };
        void import_from_file(std::string filename, const char* var_name, Algebra::MatrixExportFormat fmt);
        void import_from_file(std::string filename, std::string var_name, Algebra::MatrixExportFormat fmt);
      };
    }
  }
}
#endif
