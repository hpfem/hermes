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

#ifndef _ARRAY_H_
#define _ARRAY_H_

#include <assert.h>

/// \file array.h

#include <Judy.h>

#ifndef INVALID_IDX
#define INVALID_IDX       ((unsigned int) -1)
#endif

/// \class Array
/// \brief Implementation of a dynamic array of items.
/// C++ encapsulation of JudyL functions.
/// Provides functionality of a dynamic array.
/// Items are classes of TYPE
template<class TYPE>
class Array {
protected:
	void *judy;

public:
	Array();
	virtual ~Array();

	/// Insert an item at position index.
	/// \param[in] idx Index to insert.
	/// \return true, if ok, else false
	bool set(Word_t idx, TYPE item);

	/// Add an item to the end of the array
	///
	Word_t add(TYPE item);
	bool exists(Word_t idx) const;
	TYPE get(Word_t idx) const;
	/// overloaded operator helpers
	TYPE operator[](Word_t idx) const;
	TYPE &operator[](Word_t idx);

	/// Delete an item at index from the array.
	/// \param[in] idx Index to delete.
	void remove(Word_t idx);

	/// Count the number of indices present in the array between index1 and index2 (inclusive).
	/// \param[in] index1 Starting index. Optional, default value \c 0.
	/// \param[in] index2 End index. Optional, default value <c>(Word_t) -1</c>.
	/// \return Number of indices present in the array between index1 and index2 (inclusive).
	Word_t count(Word_t index1 = 0, Word_t index2 = (Word_t) -1) const;

	/// Locate the <c>nth</c> index that is present in the array.
	/// \param[in] nth Number of the index we want to locate. \c nth equal to 1 returns the first index present in the array.
	/// \return
	/// 	\li N-th index in the array if it is found,
	/// 	\li \c INVALID_IDX otherwise.
	Word_t by_count(Word_t nth) const;

	/// Get the number of bytes currently in use for the array.
	/// \return Number of bytes currently in use for the array.
	Word_t mem_used() const;

	/// Remove all items from the array
	/// Does NOT destruct items in the array
	void remove_all();

	// Iterators

	/// Get the first index that is present and is equal to or greater than the passed \c idx.
	/// Typically used to begin an iteration over all indices present in the array.
	/// \param[in] idx Optional, default value \c 0 (finds the first present index).
	/// \return
	/// 	\li First index present in the array that is equal or greater than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t first(Word_t idx = 0) const;

	/// Get the first index that is present and is greater than the passed \c idx.
	/// Typically used to continue an iteration over all indices present in the array.
	/// \param[in] idx Index whose succesor we want to find. Optional, default value \c 0.
	/// \return
	/// 	\li First idx present in the array that is greater than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t next(Word_t idx = 0) const;

	/// Get the last index present in the array that is equal to or less than the passed \c idx.
	/// Typically used to begin a reverse iteration over all indices present in the array.
	/// \param[in] idx Optional, default value <c>(Word_t) -1</c> (finds the last index present in the array).
	/// \return
	///		\li Last index present in the array that is equal or less than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t last(Word_t idx = (Word_t) -1) const;

	/// Get the last index present in the array that is less than the passed \c idx.
	/// Typically used to continue a reverse iteration over all indices present in the array.
	/// \param[in] idx Index whose predecessor we want to find. Optional, default value <c>(Word_t) -1</c>.
	/// \return
	/// 	\li Last index present in the array that is less than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t prev(Word_t idx = (Word_t) -1) const;

protected:
	void free_item(Word_t idx);
};

// implementation

template<class TYPE>
Array<TYPE>::Array() {
	judy = NULL;
};

template<class TYPE>
Array<TYPE>::~Array() {
	remove_all();
}

template<class TYPE>
bool Array<TYPE>::set(Word_t idx, TYPE item) {
	void *pval;
	JLG(pval, judy, idx);
	if (pval == NULL) {
		JLI(pval, judy, idx);
		if (pval == NULL)
			return false;
		else {
			*((TYPE **) pval) = new TYPE;
			**((TYPE **) pval) = item;
			return true;
		}
	}
	else {
		**((TYPE **) pval) = item;
		return true;
	}
}

template<class TYPE>
Word_t Array<TYPE>::add(TYPE item) {
	int rc;
	Word_t idx = last();
	if (idx == INVALID_IDX) idx = 0;
	JLFE(rc, judy, idx);
	if (rc) {
		set(idx, item);
		return idx;
	}
	else
		return INVALID_IDX;
}

template<class TYPE>
bool Array<TYPE>::exists(Word_t idx) const {
	void *pval;
	JLG(pval, judy, idx);
	return pval != NULL;
}

template<class TYPE>
TYPE Array<TYPE>::get(Word_t idx) const {
	void *pval;
	JLG(pval, judy, idx);
	assert(pval != NULL);
	return **((TYPE **) pval);
}

template<class TYPE>
TYPE Array<TYPE>::operator[](Word_t idx) const {
	return get(idx);
}

template<class TYPE>
TYPE &Array<TYPE>::operator[](Word_t idx) {
	void *pval;
	JLG(pval, judy, idx);
	if (pval == NULL) {
		JLI(pval, judy, idx);
		*((TYPE **) pval) = new TYPE;
	}

	return **((TYPE **) pval);
}

template<class TYPE>
void Array<TYPE>::free_item(Word_t idx) {
	void *pval;
	JLG(pval, judy, idx);
	if (pval != NULL) {
		TYPE *n = *((TYPE **) pval);
		delete n;
	}
}

template<class TYPE>
void Array<TYPE>::remove(Word_t idx) {
	free_item(idx);
	int rc;
	JLD(rc, judy, idx);
}

template<class TYPE>
Word_t Array<TYPE>::count(Word_t index1/* = 0*/, Word_t index2/* = (Word_t) -1*/) const {
	Word_t count;
	JLC(count, judy, index1, index2);
	return count;
}

template<class TYPE>
Word_t Array<TYPE>::by_count(Word_t nth) const {
	void *pval;
	Word_t index;
	JLBC(pval, judy, nth, index);
	return pval ? index : INVALID_IDX;
}

template<class TYPE>
Word_t Array<TYPE>::mem_used() const {
	Word_t memused;
	JLMU(memused, judy);
	return memused;
}

template<class TYPE>
void Array<TYPE>::remove_all() {
	for (Word_t i = first(); i != INVALID_IDX; i = next(i))
		free_item(i);

	int val;
	JLFA(val, judy);
}

// Iterators

template<class TYPE>
Word_t Array<TYPE>::first(Word_t idx/* = 0*/) const {
	void *pval;
	JLF(pval, judy, idx);
	return pval ? idx : INVALID_IDX;
}

template<class TYPE>
Word_t Array<TYPE>::next(Word_t idx/* = 0*/) const {
	void *pval;
	JLN(pval, judy, idx);
	return pval ? idx : INVALID_IDX;
}

template<class TYPE>
Word_t Array<TYPE>::last(Word_t idx/* = (Word_t) -1*/) const {
	void *pval;
	JLL(pval, judy, idx);
	return pval ? idx : INVALID_IDX;
}

template<class TYPE>
Word_t Array<TYPE>::prev(Word_t idx/* = (Word_t) -1*/) const {
	void *pval;
	JLP(pval, judy, idx);
	return pval ? idx : INVALID_IDX;
}

#endif
