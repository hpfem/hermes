#ifndef __H2D_MIXINS_H
#define __H2D_MIXINS_H
#include "global.h"
namespace Hermes
{
  namespace Hermes2D
  {
    /** \defgroup g_mixins2d Mixins
    *  \brief Mixins are utility classes used for all kinds of other classes.
    *
    *  Mixin classes provide a single piece of functionality.
    *
    */

    /// \ingroup g_mixins2d
    /// \brief Namespace for mixin classes.
    /// These classes always serve one particular purpose that multiple classes of the Hermes2D library
    /// could use - setting of spaces, output of linear algebraic structures, ...
    namespace Mixins
    {
       /// \ingroup g_mixins2d
      /// Mixin that allows for asking about the instance state (ok / not ok).
      class HERMES_API StateQueryable
      {
      public:
        /// Ask if the instance is fine.
        virtual bool isOkay() const = 0;

        /// Get class name, for the purpose of messaging.
        virtual std::string getClassName() const = 0;

        /// Method to handle the state.
        void check() const;
      };

      /// \ingroup g_mixins2d
      /// Any XML parsing class should inherit from this mixin.
      /// It serves various purposes, first of which is disabling / re-enabling of validation
      /// against the schema referenced in a file being loaded.
      class HERMES_API XMLParsing
      {
      public:
        /// Constructor.
        XMLParsing();

        /// Set to validate / not to validate.
        void set_validation(bool to_set);

      protected:
        /// Internal.
        bool validate;
      };

      /// \ingroup g_mixins2d
      /// Mixin that interfaces linear algebra structures output.
      template<typename Scalar>
      class HERMES_API MatrixRhsOutput
      {
      public:
        /// Constructor.
        /// Sets defaults (see individual set methods for values of those).
        MatrixRhsOutput();

        /// Processes the matrix.
        void process_matrix_output(SparseMatrix<Scalar>* matrix, int iteration);
        
        /// Processes the matrix.
        void process_vector_output(Vector<Scalar>* rhs, int iteration);

        /// Sets this instance to output the matrix in several first iterations.
        /// \param[in] firstIterations Only during so many first iterations. Default: -1 meaning, that during all iterations, the matrix will be saved.
        void output_matrix(int firstIterations = -1);
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
        void set_matrix_E_matrix_dump_format(EMatrixDumpFormat format);
        /// Sets number format for the matrix output.
        /// Default: "%lf".
        void set_matrix_number_format(char* number_format);
        
        /// Sets this instance to output the rhs in several first iterations.
        /// \param[in] firstIterations Only during so many first iterations. Default: -1 meaning, that during all iterations, the rhs will be saved.
        void output_rhs(int firstIterations = -1);
        /// Sets filename for the rhs
        /// Default: Rhs_'iteration number' with the ".m" extension in the case of matlab format.
        /// \param[in] name sets the main part of the name, i.e. replacement for "Rhs_" in the default name.
        void set_rhs_filename(std::string name);
        /// Sets varname for the rhs
        /// Default: "b".
        void set_rhs_varname(std::string name);
        /// Sets varname for the rhs
        /// Default: "DF_MATLAB_SPARSE - matlab file".
        void set_rhs_E_matrix_dump_format(EMatrixDumpFormat format);
        /// Sets number format for the vector output.
        /// Default: "%lf".
        void set_rhs_number_format(char* number_format);
        
      protected:
        bool print_matrix_zero_values;
        bool output_matrixOn;
        int output_matrixIterations;
        std::string matrixFilename;
        std::string matrixVarname;
        EMatrixDumpFormat matrixFormat;
		char* matrix_number_format;

        bool output_rhsOn;
        int output_rhsIterations;
        std::string RhsFilename;
        std::string RhsVarname;
        EMatrixDumpFormat RhsFormat;
		char* rhs_number_format;
      };
    }
  }
}
#endif
