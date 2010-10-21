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

#ifndef _FREQ_MAP_H_
#define _FREQ_MAP_H_

#include <assert.h>

#include <Judy.h>

#ifndef INVALID_IDX
#define INVALID_IDX					((unsigned int) -1)
#endif

template<class KEY, class TYPE>
class FreqMap {
protected:
	void *judy_hs;
	void *judy_l;

	struct Node {
		int freq;
		TYPE data;
	};

public:
	FreqMap();
	virtual ~FreqMap();

	// Get the number of values in the FreqMap
	// @return number of (key, item) pairs
	Word_t count() const;

	/// Test if the FreqMap is empty.
	/// @return true if the FreqMap has no items, otherwise false
	bool is_empty() const;

	bool lookup(const KEY &key, TYPE &item) const;

	Word_t get_idx(const KEY &key) const;

	/// Get a value as \c iter position
	/// \param[in] iter Iterator obtained by first(), next(), last() or prev().
	/// \return item at position \c iter
	TYPE get(Word_t iter) const;

	TYPE operator[](int idx) const;

	/// Add a new (key, item) pair into the FreqMap
	/// \param[in] key Pointer to the array-of-bytes.
	/// \param[in] length Length of \c key
	/// \param[in] item Item to insert
	bool set(const KEY &key, TYPE item);

	/// Delete an item with key \c key from the FreqMap.
	/// \param[in] key Pointer to the array-of-bytes.
	/// \param[in] length Length of \c key
	bool remove(const KEY &key);

	/// Remove all items from the array
	/// Does NOT destruct items in the array
	void remove_all();

	// Iterators

	// NOTE: during the iteration, keys of items are not available

	/// Get the first index that is present and is equal to or greater than the passed \c idx.
	/// Typically used to begin an iteration over all indices present in the array.
	/// \param[in] idx Optional, default value \c 0 (finds the first present index).
	/// \return
	/// 	\li First index present in the array that is equal or greater than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t first() const;

	/// Get the first index that is present and is greater than the passed \c idx.
	/// Typically used to continue an iteration over all indices present in the array.
	/// \param[in] idx Index whose succesor we want to find. Optional, default value \c 0.
	/// \return
	/// 	\li First idx present in the array that is greater than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t next(Word_t idx) const;

	/// Get the last index present in the array that is equal to or less than the passed \c idx.
	/// Typically used to begin a reverse iteration over all indices present in the array.
	/// \param[in] idx Optional, default value <c>(Word_t) -1</c> (finds the last index present in the array).
	/// \return
	///		\li Last index present in the array that is equal or less than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t last() const;

	/// Get the last index present in the array that is less than the passed \c idx.
	/// Typically used to continue a reverse iteration over all indices present in the array.
	/// \param[in] idx Index whose predecessor we want to find. Optional, default value <c>(Word_t) -1</c>.
	/// \return
	/// 	\li Last index present in the array that is less than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t prev(Word_t idx) const;

protected:
	void free_item(Word_t idx);
};

// implementation

template<class KEY, class TYPE>
FreqMap<KEY, TYPE>::FreqMap() {
	judy_hs = NULL;
	judy_l = NULL;
};

template<class KEY, class TYPE>
FreqMap<KEY, TYPE>::~FreqMap() {
	remove_all();
};

template<class KEY, class TYPE>
Word_t FreqMap<KEY, TYPE>::count() const {
	Word_t count;
	JLC(count, judy_l, 0, -1);
	return count;
}

template<class KEY, class TYPE>
bool FreqMap<KEY, TYPE>::is_empty() const {
	return count() == 0;
}

template<class KEY, class TYPE>
bool FreqMap<KEY, TYPE>::lookup(const KEY &key, TYPE &item) const {
	void *pval;
	// check if the key exists
	JHSG(pval, judy_hs, (char *) &key, sizeof(KEY));
	if (pval == NULL)
		return false;
	else {
		// get associated value
		Word_t idx = *(Word_t *) pval;
		JLG(pval, judy_l, idx);
		if (pval == NULL)
			return false;
		else {
			item = **((TYPE **) pval);
			return true;
		}
	}
}

template<class KEY, class TYPE>
Word_t FreqMap<KEY, TYPE>::get_idx(const KEY &key) const {
	void *pval;
	// check if the key exists
	JHSG(pval, judy_hs, (char *) &key, sizeof(KEY));
	if (pval == NULL)
		return INVALID_IDX;
	else
		// get associated value
		return *(Word_t *) pval;
}

template<class KEY, class TYPE>
TYPE FreqMap<KEY, TYPE>::get(Word_t iter) const {
	void *pval;
	JLG(pval, judy_l, iter);
	assert(pval != NULL);
	return **((TYPE **) pval);
}

template<class KEY, class TYPE>
TYPE FreqMap<KEY, TYPE>::operator[](int idx) const {
	return get(idx);
}

template<class KEY, class TYPE>
bool FreqMap<KEY, TYPE>::set(const KEY &key, TYPE item) {
	void *pval;
	int rc;

	// check if the key exists
	JHSG(pval, judy_hs, (char *) &key, sizeof(KEY));
	if (pval == NULL) {
		// add to array
		Word_t idx = 0;
		JLFE(rc, judy_l, idx);
		if (!rc)
			return false;
		JLI(pval, judy_l, idx);			// insert into the array
		if (pval == NULL)
			return false;

		*((TYPE **) pval) = new TYPE;
		**((TYPE **) pval) = item;

		// store the key -> value association
		JHSI(pval, judy_hs, (char *) &key, sizeof(KEY));
		*(Word_t *) pval = idx;

		return true;
	}
	else {
		// replace value in the array
		Word_t idx = *(Word_t *) pval;
		JLG(pval, judy_l, idx);			// insert into the array
		if (pval == NULL)
			return false;

		**((TYPE **) pval) = item;
		return true;
	}
}

template<class KEY, class TYPE>
bool FreqMap<KEY, TYPE>::remove(const KEY &key) {
	void *pval;
	// check if the key exists
	JHSG(pval, judy_hs, (char *) &key, sizeof(KEY));
	if (pval == NULL) {
		// removing non-existent item
		return false;
	}
	else {
		bool res = true;
		Word_t idx = *(Word_t *) pval;
		free_item(idx);
		int rc;
		JLD(rc, judy_l, idx);
		res &= rc == 1;
		JHSD(rc, judy_hs, (char *) &key, sizeof(KEY));
		res &= rc == 1;
		return res;
	}
}

template<class KEY, class TYPE>
void FreqMap<KEY, TYPE>::remove_all() {
	// free associated memory
	void *pval;
	Word_t idx = 0;
	JLF(pval, judy_l, idx);
	for (; idx != -1 && pval != NULL; ) {
		free_item(idx);
		JLN(pval, judy_l, idx);
	}

	int val;
	// clean array
	JLFA(val, judy_l);
	// clean index hash FreqMap
	JHSFA(val, judy_hs);
}

template<class KEY, class TYPE>
void FreqMap<KEY, TYPE>::free_item(Word_t idx) {
	void *pval;
	JLG(pval, judy_l, idx);
	if (pval != NULL) {
		TYPE *n = *((TYPE **) pval);
		delete n;
	}
}

template<class KEY, class TYPE>
Word_t FreqMap<KEY, TYPE>::first() const {
	void *pval;
	Word_t idx = 0;
	JLF(pval, judy_l, idx);
	return pval ? idx : INVALID_IDX;
}

template<class KEY, class TYPE>
Word_t FreqMap<KEY, TYPE>::next(Word_t idx) const {
	void *pval;
	JLN(pval, judy_l, idx);
	return pval ? idx : INVALID_IDX;
}

template<class KEY, class TYPE>
Word_t FreqMap<KEY, TYPE>::last() const {
	void *pval;
	Word_t idx = -1;
	JLL(pval, judy_l, idx);
	return pval ? idx : INVALID_IDX;
}

template<class KEY, class TYPE>
Word_t FreqMap<KEY, TYPE>::prev(Word_t idx) const {
	void *pval;
	JLP(pval, judy_l, idx);
	return pval ? idx : INVALID_IDX;
}

#endif
