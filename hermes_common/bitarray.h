// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _BITARRAY_H_
#define _BITARRAY_H_

/// \file bitarray.h

#include <Judy.h>

#ifndef INVALID_IDX
#define INVALID_IDX					((unsigned int) -1)
#endif


/// \class BitArray
/// \brief Implementation of a hash table mapping an Index to a bit (0/1).
/// C++ encapsulation of Judy1 functions.
/// Provides functionality of a sparse bit vector.
class BitArray {
	void *judy;

public:
	BitArray() : judy(NULL) {}
	~BitArray() { free(); }

	/// Set <c>index</c>'s bit in the array.
	/// \param[in] index index to set.
	/// \retval false if the bit was already set.
	/// \retval true if the bit was previously unset and was succesfully set.
	bool set(Word_t index) {
		int rc;
		J1S(rc, judy, index);
		return rc == 1;
	}

	/// Unset <c>index</c>'s bit in the array.
	/// \param[in] index index to unset.
	/// \retval false if the bit was already unset.
	/// \retval true if the bit was previously set and was succesfully unset (removet from the array).
	bool unset(Word_t index) {
		int rc;
		J1U(rc, judy, index);
		return rc == 1;
	}

	/// Test whether <c>index</c>'s bit is set.
	/// \param[in] index Index to test.
	/// \retval false if the index is unset (absent in the array).
	/// \retval true if the index is set (present in the array).
	bool is_set(Word_t index) {
		int rc;
		J1T(rc, judy, index);
		return rc == 1;
	}

	/// Count the number of indices present in the array between index1 and index2 (inclusive).
	/// \param[in] index1 Starting index. Optional, default value \c 0.
	/// \param[in] index2 End index. Optional, default value <c>(Word_t) -1</c>.
	/// \return Number of indices present in the array between index1 and index2 (inclusive).
	Word_t count(Word_t index1 = 0, Word_t index2 = (Word_t) -1) const {
		Word_t count;
		J1C(count, judy, index1, index2);
		return count;
	}

	/// Locate the <c>nth</c> index that is present in the array.
	/// \param[in] nth N-th index we want to locate. \c nth equal to 1 returns the first index set.
	/// \return N-th set index in the array if it is found, \c INVALID_IDX otherwise.
	Word_t by_count(Word_t nth) const {
		int rc;
		Word_t index;
		J1BC(rc, judy, nth, index);
		return rc ? index : INVALID_IDX;
	}

	/// Get the number of bytes currently in use for the array.
	/// \return Number of bytes currently in use for the array.
	Word_t mem_used() const {
		Word_t memused;
		JLMU(memused, judy);
		return memused;
	}

	/// Get the first index that is set and is equal to or greater than the passed \c index.
	/// Typically used to begin an iteration over all set indices in the array.
	/// \param[in] index Optional, default value \c 0 (finds the first set index).
	/// \return \li First set index that is equal or greater than the passed \c index (if found), \li \c INVALID_IDX (if not found).
	Word_t first(Word_t index = 0) const {
		int rc;
		J1F(rc, judy, index);
		return rc ? index : INVALID_IDX;
	}

	/// Get the first set index that is greater than the passed \c index.
	/// Typically used to continue an iteration over all set indices in the array.
	/// \param[in] index Index whose succesor we want to find. Optional, default value \c 0.
	/// \return \li First set index that is greater than the passed \c index (if found), \li \c INVALID_IDX (if not found).
	Word_t next(Word_t index = 0) const {
		int rc;
		J1N(rc, judy, index);
		return rc ? index : INVALID_IDX;
	}

	/// Get the last index that is set and is equal to or less than the passed \c index.
	/// Typically used to begin a reverse iteration over all set indices in the array.
	/// \param[in] index Optional, default value <c>(Word_t) -1</c> (finds the last set index).
	/// \return \li Last set index that is equal or less than the passed \c index (if found), \li \c INVALID_IDX (if not found).
	Word_t last(Word_t index = (Word_t) -1) const {
		int rc;
		J1L(rc, judy, index);
		return rc ? index : INVALID_IDX;
	}

	/// Get the last set index that is less than the passed \c index.
	/// Typically used to continue a reverse iteration over all set indices in the array.
	/// \param[in] index Index whose predecessor we want to find. Optional, default value <c>(Word_t)-1</c>.
	/// \return \li Last set index that is less than the passed \c index (if found), \li \c INVALID_IDX (if not found).
	Word_t prev(Word_t index = (Word_t) -1) const {
		int rc;
		J1P(rc, judy, index);
		return rc ? index : INVALID_IDX;
	}

	/// Get the first unset index in the array that is equal to or greater than the passed \c index.
	/// \param[in] index Optional, default value \c 0 (finds the first unset index).
	/// \return \li First unset index that is equal or greater than the passed \c index (if found), \li \c INVALID_IDX (if not found).
	Word_t first_empty(Word_t index = 0) const {
		int rc;
		J1FE(rc, judy, index);
		return rc ? index : INVALID_IDX;
	}

	/// Get the first unset index that is greater than passed \c index.
	/// \param[in] index Index whose unset succesor we want to find. Optional, default value \c 0.
	/// \return \li First unset index that is greater than the passed \c index (if found), \li \c INVALID_IDX (if not found).
	Word_t next_empty(Word_t index = 0) const {
		int rc;
		J1NE(rc, judy, index);
		return rc ? index : INVALID_IDX;
	}

	/// Get the last unset index that is equal to or less than the passed \c index.
	/// \param[in] index Optional, default value <c>(Word_t)-1</c> (finds the last unset index).
	/// \return \li Last unset index that is equal or less than the passed \c index (if found), \li \c INVALID_IDX (if not found).
	Word_t last_empty(Word_t index = (Word_t) -1) const {
		int rc;
		J1LE(rc, judy, index);
		return rc ? index : INVALID_IDX;
	}

	/// Get the last unset index that is less than passed \c index.
	/// \param[in] index Index whose unset predecessor we want to find. Optional, default value <c>(Word_t)-1</c>.
	/// \return \li Last unset index that is less than the passed \c index (if found), \li \c INVALID_IDX (if not found).
	Word_t prev_empty(Word_t index = (Word_t) -1) const {
		int rc;
		J1PE(rc, judy, index);
		return rc ? index : INVALID_IDX;
	}

	/// Make a copy of array
	bool copy(BitArray *original) {
		for (Word_t ind = original->first(); ind != INVALID_IDX; ind = original->next(ind)) {
			if (!set(ind))
				return false;
		}
		return true;
	}

	/// unset all elements
	void free() {
		int val;
		J1FA(val, judy);
		judy = NULL;
	}
};

#endif
