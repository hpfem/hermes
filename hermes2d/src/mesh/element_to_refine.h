// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_ELEMENT_TO_REFINE_H
#define __H2D_ELEMENT_TO_REFINE_H

#include "refinement_type.h"

/// A refinement record. \ingroup g_adapt
/** Except the attribute ElementToRefine::q, the class is able to dump its content to a stringstream
 *  in a human readable form, see operator<<(std::ostream& stream, const ElementToRefine& elem_ref). */
class HERMES_API ElementToRefine {
public:
  int id; ///< An ID of the element. -1 if invalid.
  int comp; ///< An index of the component. -1 if invalid.
  int split; ///< Proposed refinement. Possible values are defined in the enum ::RefinementType.
  int p[H2D_MAX_ELEMENT_SONS]; ///< Encoded orders of sons.
  int q[H2D_MAX_ELEMENT_SONS]; ///< Encoded orders of the best H-refinement. These orders are used in a case multiple components shares a single mesh.

public:
  /// Constructor. Creates an invalid refinement.
  ElementToRefine() : id(-1), comp(-1) {};

  /// Constructor.
  /** \param[in] id An ID of the element.
   *  \param[in] comp An index of a component. */
  ElementToRefine(int id, int comp) : id(id), comp(comp), split(H2D_REFINEMENT_H) {};

  /// Copy-contructor.
  ElementToRefine(const ElementToRefine &orig) : id(orig.id), comp(orig.comp), split(orig.split) {
    copy_orders(p, orig.p);
    copy_orders(q, orig.q);
  };

  /// Assignment operator.
  ElementToRefine& operator=(const ElementToRefine& orig);

  /// Returns a number of sons.
  /** \return A number of sons of a given refinement. */
  int get_num_sons() const { return get_refin_sons(split); };

  /// Copies array of orders.
  /** The length of the array is defubed by ::H2D_MAX_ELEMENT_SONS.
   *  \param[in] dest A destination array.
   *  \param[in] src A source arrapy. */
  static inline void copy_orders(int* dest, const int* src) {
    memcpy(dest, src, sizeof(int) * H2D_MAX_ELEMENT_SONS);
  }
};

/// Dumps contents of the structure to a stream in a text form. \ingroup g_adapt
/** Used for debugging purposes to dump the structure to a stringstream. */
extern HERMES_API std::ostream& operator<<(std::ostream& stream, const ElementToRefine& elem_ref);

/// Checks whether the input contains a given tag. Throws an std::runtime_error if the tag is not found.
/** \section s_use Suggested Use
 *  The class works with std::istream. The following code shows checking for tag 'TAG' in the input file stream:
 *  \code
    std::ifstream fin("my_file.txt");
    fin >> TagChecker("TAG");
 *  \endcode */
class HERMES_API TagChecker {
  const std::string& tag; ///< A tag.
public:
  /// Constructor.
  /** \param[in] tag A tag. */
  explicit TagChecker(const std::string& tag) : tag(tag) {};

  /// Returns tag.
  /** \return A tag. Returned reference becomes invalid after the instance is released */
  const std::string& get_tag() const { return tag; };
};

/// Peformes a check for a tag. Throws exception otherwise. \sa TagChecker.
extern HERMES_API std::istream& operator>>(std::istream& stream, const TagChecker& checker);

/// Refinement stream capable of input and output. \ingroup g_adapt
/** This class allows to store a vector of refinements to a file.
 *  It is intended for debugging purposes when a certain problem
 *  appearts at higher iteration. The stream is designed to work with
 *  the class Adapt.
 *
 *  \section s_use Suggested Use
 *  The stream is able to read, store and append. If appending, an existing file
 *  has to have a valid format. This is checked when the stream is opened.
 *
 *  Assuming that \c it=1 is the first iteration and the variable \c hp is an instance of
 *  the class Adapt, the code that writes the refinements is \code
    std::ios_base::openmode mode = std::ios_base::app;
    if (it == 2)
      mode = std::ios_base::out;
    ElementToRefineStream erstream("data.ers", mode | std::ios_base::binary);
    erstream << hp.get_last_refinements();
    erstream.close();
 *  \endcode
 *
 *  The code that reads and applies all refinements in a file is \code
    ElementToRefineStream erstream("data.ers", std::ios_base::in | std::ios_base::binary);
    if (erstream.is_open()) {
      while (erstream.good()) {
        std::vector<ElementToRefine> refinements;
        errstream >> refinements;
        hp.apply_refinements(refinements);
      }
      errstream.close();
    }
 *  \endcode
 *
 *  \section s_format File Format
 *  Currently, only binary format is supported. The file consist of a header and
 *  a sequence of refinements organized to blocks. All binary data are little-endian, always.
 *  All numbers are signed if not noted otherwise and all sizes are in bytes.
 *
 *  The file header is: \verbatim
    ------------------------------------------------------------------
    field         size          contents
    ------------------------------------------------------------------
    Start Tag     3             'ERS' (i.e. Element to Refine Stream)
    Delimiter     1             A space (ASCII: 32)
    Format Tag    3             'BIN' (i.e. Uncompressed Binary format)
    New Line      1             LF or '\n' \endverbatim
    ------------------------------------------------------------------\endverbatim
 *
 *  The header of a refinement block is: \verbatim
    ------------------------------------------------------------------
    name          size          contents
    ------------------------------------------------------------------
    Start Tag     3             'ERV' (i.e. Element to Refine Vector)
    SZ_NUM        1             A byte length of NUM
    SZ_ID_ROOT    1             A byte length of ID_ROOT
    SZ_ID_OFFS    1             A byte length of ID_OFFS
    SZ_COMP_ROOT  1             A byte length of COMP_ROOT
    SZ_COMP_OFFS  1             A byte length of COMP_OFFS
    SZ_ORDER      1             A byte length of an order
    NUM           SZ_NUM        A number of refinements in the block.
    ID_ROOT       SZ_ID_ROOT    A root element ID. All IDs
                                in the block are relative.
    COMP_ROOT     SZ_COMP_ROOT  A root component index.
                                All indices in the block are relative.
    ------------------------------------------------------------------\endverbatim
 *
 *  The format of a refinement is: \verbatim
    ------------------------------------------------------------------
    name          size          contents
    ------------------------------------------------------------------
    ID_OFFS       SZ_ID_OFFS    An offset of an element ID.
    COMP_OFFS     SZ_COMP_OFFS  An offset of a component index.
                                If SZ_COMP_OFFS is zero, COMP_OFFS is
                                zero and it is not present in the file.
    SPLIT         1             A refinement type, see RefinementTypes.
    ORDERS        2*SZ_ORDER*N  A sequence of pairs of orders in
                                the attribute ElementToRefine::p, where
                                N is a number of sons generated by the
                                refinement.
                                A pair consists of a horizontal and
                                a vertical orders. The vertical order
                                is zero in a case of a triangle.
    ------------------------------------------------------------------\endverbatim
 */
class HERMES_API ElementToRefineStream {
private:
  static const char* H2DER_START_TAG; ///< A tag which start file.
  static const char* H2DER_BIN_TAG; ///< A tag which defines uncompressed binary contents.
  static const char* H2DER_VECTOR_TAG; ///< A tag which starts a block of refinements.
  static const int H2DER_SIZE_BYTESIZE = 1; ///< A size of all size specifications.

  std::fstream stream; ///< An underlaying stream.
  bool little_endian; ///< True if the system is little endian.

  /// Returns a number of bytes necessary to store a given value in a signed form.
  /** \return A number of bytes necessary to store a given value in a signed form. */
  static uint8_t get_byte_size(int value);

  /// Writes an integer in a little-edinan form to a file.
  /** \todo Test on a BIG-endian machine.
   *  \param[in] integer_ptr_c A pointer to a first byte of an integer in a system-specific form.
   *  \param[in] num_bytes A number of bytes to write from the integers. The fist byte is the least-significant byte. */
  void write_bytes(const void* integer_ptr, int num_bytes);

  /// Writes an integer in a little-endian form to a file
  /** \overload void write_bytes(const int integer, int num_bytes) */
  void write_bytes(const int integer, int num_bytes);

  /// Reads bytes from a file an convert them to a system-specific signed integer.
  /** \param[in] num_bytes A number of bytes to read. Notice \a num_bytes < \c sizeof(int).
   *  \return A decoded and signed integer */
  int read_bytes(int num_bytes);

  /// Writes a header.
  /** \param[in] mode A stream open mode. */
  void write_header(std::ios_base::openmode mode);

  /// Reads a header.
  /** \param[in] mode A stream open mode.
   *  \return True if header was read. */
  bool read_header(std::ios_base::openmode mode);

public:
  /// Constructor. Opens a stream.
  /** \param[in] filename A file name.
   *  \param[in] mode An open mode. The file has to be opened in a binary mode (std::ios_base::binary) and it can be opened either for reading (std::ios_base::in), for writing (std::ios_base::out) or for appending (std::ios_base::append). */
  explicit ElementToRefineStream(const char* filename, std::ios_base::openmode mode);

  /// Opens a file stream.
  /** \param[in] filename A file name.
   *  \param[in] mode An open mode. The file has to be opened in a binary mode (std::ios_base::binary) and it can be opened either for reading (std::ios_base::in), for writing (std::ios_base::out) or for appending (std::ios_base::append). */
  void open(const char* filename, std::ios_base::openmode mode);

  /// Returns true if the stream is open. Used to test success of the method open().
  /** \return True if the stream is open. */
  bool is_open() const { return stream.is_open(); };

  /// Returns true if EOF has been reached.
  /** \return True if EOF has been reached. */
  bool eof() const { return stream.eof(); };

  /// Returns true if the stream did not failed and EOF has not been reached.
  /** \return True if the stream did not failed and EOF has not been reached, i.e., neither failbit nor badbit nor oefbit are set. */
  bool good() const { return stream.good(); };

  /// Returns true if the stream failed (to read header, ...).
  /** \return True if the stream failed, i.e., failbit or badbit is set. */
  bool fail() const { return stream.fail(); };

  /// Returns true if the an error occured.
  /** \return True if the stream fialed (to read a header, etc.) */
  bool operator!() const { return stream.fail(); };

  /// Closes the stream.
  void close() { stream.close(); };

  friend HERMES_API ElementToRefineStream& operator<<(ElementToRefineStream& stream, const std::vector<ElementToRefine>& elem_refs);
  friend HERMES_API ElementToRefineStream& operator>>(ElementToRefineStream& stream, std::vector<ElementToRefine>& elem_refs);
};

/// Operator. Stores a vector of refinements to a stream. \ingroup g_adapt
extern HERMES_API ElementToRefineStream& operator<<(ElementToRefineStream& stream, const std::vector<ElementToRefine>& elem_refs);

/// Operator. Reads a vector of refinemens. It erases the contents of the vector that it is reading to. \ingroup g_adapt
extern HERMES_API ElementToRefineStream& operator>>(ElementToRefineStream& stream, std::vector<ElementToRefine>& elem_refs);

#endif
