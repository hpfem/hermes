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

#ifndef _ARRAYPTR_H_
#define _ARRAYPTR_H_

#include <assert.h>

/// \file arrayptr.h

#include <Judy.h>

#ifndef INVALID_IDX
#define INVALID_IDX					((unsigned int) -1)
#endif

/// \class JudyArrayPtr
/// \brief Implementation of a dynamic array of pointers to items.
/// C++ encapsulation of JudyL functions.
/// Provides functionality of a dynamic array.
/// Items are pointers to class TYPE
template<class TYPE>
class JudyArrayPtr {
protected:
	void *judy;

public:
	JudyArrayPtr();
	virtual ~JudyArrayPtr();

	/// Inserts an item at position index.
	/// \param[in] idx Index to insert.
	/// \return true, if ok, else false
	bool set(Word_t idx, TYPE *item);

	/// Adds an item to the end of the array
	///
	Word_t add(TYPE *item);
	bool exists(Word_t idx) const;
	TYPE *get(Word_t idx) const;
        /// overloaded operator helpers
	TYPE *operator[](Word_t idx) const;
	TYPE *&operator[](Word_t idx);

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

	/// Removes all items from the array
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
};

// implementation

template<class TYPE>
JudyArrayPtr<TYPE>::JudyArrayPtr() {
	judy = NULL;
};

template<class TYPE>
JudyArrayPtr<TYPE>::~JudyArrayPtr() {
	remove_all();
}

template<class TYPE>
bool JudyArrayPtr<TYPE>::set(Word_t idx, TYPE *item) {
	void *pval;
	JLG(pval, judy, idx);
	if (pval == NULL) {
		JLI(pval, judy, idx);
		if (pval == NULL)
			return false;
		else {
			*((void **) pval) = (void *) item;
			return true;
		}
	}
	else {
		*((void **) pval) = (void *) item;
		return true;
	}
}

template<class TYPE>
Word_t JudyArrayPtr<TYPE>::add(TYPE *item) {
	int rc;
	Word_t idx = last();
	JLFE(rc, judy, idx);
	if (rc) {
		set(idx, item);
		return idx;
	}
	else
		return INVALID_IDX;
}

template<class TYPE>
bool JudyArrayPtr<TYPE>::exists(Word_t idx) const {
	void *pval;
	JLG(pval, judy, idx);
	return pval != NULL;
}

template<class TYPE>
TYPE *JudyArrayPtr<TYPE>::get(Word_t idx) const {
	void *pval;
	JLG(pval, judy, idx);
	assert(pval != NULL);
	return *((TYPE **) pval);
}

template<class TYPE>
TYPE *JudyArrayPtr<TYPE>::operator[](Word_t idx) const {
	return get(idx);
}

template<class TYPE>
TYPE *&JudyArrayPtr<TYPE>::operator[](Word_t idx) {
	void *pval;
	JLG(pval, judy, idx);
	if (pval == NULL) {
		JLI(pval, judy, idx);
	}

	return *((TYPE **) pval);
}

template<class TYPE>
void JudyArrayPtr<TYPE>::remove(Word_t idx) {
	int rc;
	JLD(rc, judy, idx);
}

template<class TYPE>
Word_t JudyArrayPtr<TYPE>::count(Word_t index1/* = 0*/, Word_t index2/* = (Word_t) -1*/) const {
	Word_t count;
	JLC(count, judy, index1, index2);
	return count;
}

template<class TYPE>
Word_t JudyArrayPtr<TYPE>::by_count(Word_t nth) const {
	void *pval;
	Word_t index;
	JLBC(pval, judy, nth, index);
	return pval ? index : INVALID_IDX;
}

template<class TYPE>
Word_t JudyArrayPtr<TYPE>::mem_used() const {
	Word_t memused;
	JLMU(memused, judy);
	return memused;
}

template<class TYPE>
void JudyArrayPtr<TYPE>::remove_all() {
	int val;
	JLFA(val, judy);
}

// Iterators

template<class TYPE>
Word_t JudyArrayPtr<TYPE>::first(Word_t idx/* = 0*/) const {
	void *pval;
	JLF(pval, judy, idx);
	return pval ? idx : INVALID_IDX;
}

template<class TYPE>
Word_t JudyArrayPtr<TYPE>::next(Word_t idx/* = 0*/) const {
	void *pval;
	JLN(pval, judy, idx);
	return pval ? idx : INVALID_IDX;
}

template<class TYPE>
Word_t JudyArrayPtr<TYPE>::last(Word_t idx/* = (Word_t) -1*/) const {
	void *pval;
	JLL(pval, judy, idx);
	return pval ? idx : INVALID_IDX;
}

template<class TYPE>
Word_t JudyArrayPtr<TYPE>::prev(Word_t idx/* = (Word_t) -1*/) const {
	void *pval;
	JLP(pval, judy, idx);
	return pval ? idx : INVALID_IDX;
}

#endif
