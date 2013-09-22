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
/*! \file vector.h
\brief Basic vector classes and operations.
*/
#ifndef __HERMES_COMMON_VECTOR_H
#define __HERMES_COMMON_VECTOR_H

#include "algebra_utilities.h"
#include "algebra_mixins.h"
#include "mixins.h"

namespace Hermes
{
  /// \brief Namespace containing classes for vector / matrix operations.
  namespace Algebra
  {
    /// \brief General (abstract) vector representation in Hermes.
    template<typename Scalar>
    class HERMES_API Vector : public Hermes::Mixins::Loggable, public Algebra::Mixins::MatrixRhsImportExport<Scalar>
    {
    public:
      /// Default constructor.
      Vector();
      /// Constructor of vector with specific size.
      /// @param[in] size size of vector
      Vector(unsigned int size);
      virtual ~Vector() { }

      /// allocate memory for storing ndofs elements
      ///
      /// @param[in] ndofs - number of elements of the vector
      virtual void alloc(unsigned int ndofs) = 0;
      /// free the memory
      virtual void free() = 0;
      /// finish the assembly of the vector
      virtual void finish() { }

      /// Get the value from a position
      /// @return the value form the specified index
      /// @param[in] idx - index which to obtain the value from
      virtual Scalar get(unsigned int idx) const = 0;

      /// Extract vector values into user-provided array.
      /// @param[out] v - array which will contain extracted values
      virtual void extract(Scalar *v) const = 0;

      /// Zero the vector
      virtual void zero() = 0;

      /// Multiply by minus one.
      virtual Vector<Scalar>* change_sign() = 0;

      /// set the entry on a specified position
      ///
      /// @param[in] idx - indices where to update
      /// @param[in] y   - value
      virtual void set(unsigned int idx, Scalar y) = 0;

      /// update element on the specified position
      ///
      /// @param[in] idx - indices where to update
      /// @param[in] y   - value
      virtual void add(unsigned int idx, Scalar y) = 0;

      /// Set values from a user-provided vector.
      virtual Vector<Scalar>* set_vector(Vector<Scalar>* vec);
      /// Set values from a user-provided array.
      virtual Vector<Scalar>* set_vector(Scalar* vec);

      /// Add a vector.
      virtual Vector<Scalar>* add_vector(Vector<Scalar>* vec);
      /// Add a vector.
      virtual Vector<Scalar>* add_vector(Scalar* vec);

      /// update subset of the elements
      ///
      /// @param[in] n   - number of positions to update
      /// @param[in] idx - indices where to update
      /// @param[in] y   - values
      virtual void add(unsigned int n, unsigned int *idx, Scalar *y) = 0;

      /// Get vector length.
      unsigned int get_size() const {return this->size;}
    protected:
      /// size of vector
      unsigned int size;
    };

    /** \brief Vector used with MUMPS solver */
    template <typename Scalar>
    class HERMES_API SimpleVector : public Vector<Scalar>
    {
    public:
      SimpleVector();
      SimpleVector(unsigned int size);
      virtual ~SimpleVector();

      virtual void alloc(unsigned int ndofs);
      virtual void free();
      virtual Scalar get(unsigned int idx) const;
      virtual void extract(Scalar *v) const;
      virtual void zero();
      virtual Vector<Scalar>* change_sign();
      virtual void set(unsigned int idx, Scalar y);
      virtual void add(unsigned int idx, Scalar y);
      virtual void add(unsigned int n, unsigned int *idx, Scalar *y);
      virtual Vector<Scalar>* add_vector(Vector<Scalar>* vec);
      virtual Vector<Scalar>* add_vector(Scalar* vec);
      virtual Vector<Scalar>* set_vector(Vector<Scalar>* vec);
      virtual Vector<Scalar>* set_vector(Scalar* vec);
      
      using Vector<Scalar>::export_to_file;
      using Vector<Scalar>::import_from_file;
      virtual void export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format = "%lf");
      virtual void import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt);

      /// Raw data.
      Scalar *v;
    };

    /// \brief Function returning a vector according to the users's choice.
    /// @return created vector
    template<typename Scalar> HERMES_API
      Vector<Scalar>* create_vector(bool use_direct_solver = false);
  }
}
#endif
