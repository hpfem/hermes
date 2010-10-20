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

#ifndef _MAPORD_H_
#define _MAPORD_H_

#include <assert.h>

#include <Judy.h>
#include <climits>

#ifndef INVALID_IDX
#define INVALID_IDX					((unsigned int) UINT_MAX) // ((unsigned int) -1) was here
#endif

#include "map.h"
#include "maphs.h"

/// compare procedure
static int index_cmp(const void *a, const void *b) {
	const unsigned int *da = (unsigned int *) a;
	const unsigned int *db = (unsigned int *) b;
	return (*da > *db) - (*da < *db);
}

static void sort_key(unsigned int *key, int length) {
	if (length > 2) {
		qsort(key, length, sizeof(unsigned int), index_cmp);
	}
	else if (length == 2 && key[0] > key[1]) {
		unsigned int dummy = key[1];
		key[1] = key[0];
		key[0] = dummy;
	}
}

/// Implementation of a hash map (unsigned int key[]) => TYPE
template<class TYPE>
class MapOrd {
protected:
	void *judy_hs;
	void *judy_l;

public:
	MapOrd();
	virtual ~MapOrd();

	// Get the number of values in the Map
	// @return number of (key, item) pairs
	unsigned int count() const;

	/// Test if the Map is empty.
	/// @return true if the Map has no items, otherwise false
	bool is_empty() const;

	/// Lookup for item with key \c key.
	/// \param[in] key Indices that build up the key.
	/// \param[in] length Number of indices.
	/// \param[out] item Item with key \c key
	/// \return
	/// 	\li true if the key exists, item then contains the value,
	/// 	\li false otherwise
	bool lookup(unsigned int *key, int length, TYPE &item) const;

	/// Lookup for index of item with key \c key.
	/// \param[in] key Pointer to the array-of-bytes.
	/// \param[in] length Length of \c key
	/// \return
	/// 	\li true if the key exists, \c item then contains the item,
	/// 	\li false otherwise
	unsigned int get_idx(unsigned int *key, int length) const;

	/// Get a value as \c iter position
	/// \param[in] iter Iterator obtained by first(), next(), last() or prev().
	/// \return item at position \c iter
	TYPE get(unsigned int iter) const;

	TYPE operator[](unsigned int idx) const;
//	TYPE &operator[](int idx);

	/// Add a new (key, item) pair into the map
	/// @param[in] key Indices that build up the key.
	/// @param[in] length Number of indices.
	/// @param[in] item Item to insert
	/// @return Index where the (key, item) pair was inserted
	bool set(unsigned int *key, int length, TYPE item);

	/// Delete an item with key \c key from the map.
	/// @param[in] key Indices that build up the key.
	/// @param[in] length Number of indices.
	bool remove(unsigned int *key, int length);

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
	unsigned int first() const;

	/// Get the first index that is present and is greater than the passed \c idx.
	/// Typically used to continue an iteration over all indices present in the array.
	/// \param[in] idx Index whose succesor we want to find. Optional, default value \c 0.
	/// \return
	/// 	\li First idx present in the array that is greater than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	unsigned int next(Word_t idx) const;

	/// Get the last index present in the array that is equal to or less than the passed \c idx.
	/// Typically used to begin a reverse iteration over all indices present in the array.
	/// \param[in] idx Optional, default value <c>(unsigned int) -1</c> (finds the last index present in the array).
	/// \return
	///		\li Last index present in the array that is equal or less than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	unsigned int last() const;

	/// Get the last index present in the array that is less than the passed \c idx.
	/// Typically used to continue a reverse iteration over all indices present in the array.
	/// \param[in] idx Index whose predecessor we want to find. Optional, default value <c>(unsigned int) -1</c>.
	/// \return
	/// 	\li Last index present in the array that is less than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	unsigned int prev(unsigned int idx) const;

protected:
	void free_item(unsigned int idx);
};


// implementation

template<class TYPE>
MapOrd<TYPE>::MapOrd() {
	judy_hs = NULL;
	judy_l = NULL;
};

template<class TYPE>
MapOrd<TYPE>::~MapOrd() {
	remove_all();
};

template<class TYPE>
unsigned int MapOrd<TYPE>::count() const {
	unsigned int count;
	JLC(count, judy_l, 1, -1);
	return count;
}

template<class TYPE>
bool MapOrd<TYPE>::is_empty() const {
	return count() == 0;
}

template<class TYPE>
bool MapOrd<TYPE>::lookup(unsigned int *key, int length, TYPE &item) const {
	void *pval;

	sort_key(key, length);
	// check if the key exists
	JHSG(pval, judy_hs, key, length * sizeof(unsigned int));
	if (pval == NULL)
		return false;
	else {
		// get associated value
		unsigned int idx = *(unsigned int *) pval;
		JLG(pval, judy_l, idx);
		if (pval == NULL)
			return false;
		else {
			item = **((TYPE **) pval);
			return true;
		}
	}
}

template<class TYPE>
unsigned int MapOrd<TYPE>::get_idx(unsigned int *key, int length) const {
	void *pval;

	sort_key(key, length);
	// check if the key exists
	JHSG(pval, judy_hs, key, length * sizeof(unsigned int));
	if (pval == NULL)
		return INVALID_IDX;
	else
		// get associated value
		return *(unsigned int *) pval;
}

template<class TYPE>
TYPE MapOrd<TYPE>::get(unsigned int iter) const {
	void *pval;
	JLG(pval, judy_l, iter);
	assert(pval != NULL);
	return **((TYPE **) pval);
}

template<class TYPE>
TYPE MapOrd<TYPE>::operator[](unsigned int idx) const {
	return get(idx);
}

//template<class TYPE>
//TYPE &MapOrd<TYPE>::operator[](int idx) {
//	void *pval;
//	JLG(pval, judy_l, idx);
//	assert(pval != NULL);
//
//	return **((TYPE **) pval);
//}

template<class TYPE>
bool MapOrd<TYPE>::set(unsigned int *key, int length, TYPE item) {
	void *pval;
	int rc;

	sort_key(key, length);
	// check if the key exists
	JHSG(pval, judy_hs, key, length * sizeof(unsigned int));
	if (pval == NULL) {
		// add to array
		Word_t idx = 1;
		JLFE(rc, judy_l, idx);
		if (!rc)
			return false;
		JLI(pval, judy_l, idx);			// insert into the array
		if (pval == NULL)
			return false;

		*((TYPE **) pval) = new TYPE;
		**((TYPE **) pval) = item;

		// store the key -> value association
		JHSI(pval, judy_hs, key, length * sizeof(unsigned int));
		*(unsigned int *) pval = idx;

		return true;
	}
	else {
		// replace value in the array
		unsigned int idx = *(unsigned int *) pval;
		JLG(pval, judy_l, idx);			// insert into the array
		if (pval == NULL)
			return false;

		**((TYPE **) pval) = item;
		return true;
	}
}

template<class TYPE>
bool MapOrd<TYPE>::remove(unsigned int *key, int length) {
	void *pval;

	sort_key(key, length);
	// check if the key exists
	JHSG(pval, judy_hs, key, length * sizeof(unsigned int));
	if (pval == NULL) {
		// removing non-existent item
		return false;
	}
	else {
		bool res = true;
		unsigned int idx = *(unsigned int *) pval;
		free_item(idx);
		int rc;
		JLD(rc, judy_l, idx);
		res &= rc == 1;
		JHSD(rc, judy_hs, key, length * sizeof(unsigned int));
		res &= rc == 1;
		return res;
	}
}

template<class TYPE>
void MapOrd<TYPE>::remove_all() {
	// free associated memory
	void *pval;
	Word_t idx = 1;
	JLF(pval, judy_l, idx);
	for (; idx != INVALID_IDX && pval != NULL; ) {
		free_item(idx);
		JLN(pval, judy_l, idx);
	}

	int val;
	// clean array
	JLFA(val, judy_l);
	// clean index hash map
	JHSFA(val, judy_hs);
}

template<class TYPE>
void MapOrd<TYPE>::free_item(unsigned int idx) {
	void *pval;
	JLG(pval, judy_l, idx);
	if (pval != NULL) {
		TYPE *n = *((TYPE **) pval);
		delete n;
	}
}

template<class TYPE>
unsigned int MapOrd<TYPE>::first() const {
	void *pval;
	Word_t idx = 1;
	JLF(pval, judy_l, idx);
	return pval ? idx : INVALID_IDX;
}

template<class TYPE>
unsigned int MapOrd<TYPE>::next(Word_t idx) const {
	void *pval;
	JLN(pval, judy_l, idx);
	return pval ? idx : INVALID_IDX;
}

template<class TYPE>
unsigned int MapOrd<TYPE>::last() const {
	void *pval;
	unsigned int idx = -1;
	JLL(pval, judy_l, (Word_t)idx);
	return pval ? idx : INVALID_IDX;
}

template<class TYPE>
unsigned int MapOrd<TYPE>::prev(unsigned int idx) const {
	void *pval;
	JLP(pval, judy_l, (Word_t)idx);
	return pval ? idx : INVALID_IDX;
}

//
//
//

static unsigned int *index(unsigned int *key, int length) {
	unsigned int *p = new unsigned int[length];
	for (int k = 0; k < length; k++)
		p[k] = key[k];

	if (length > 2) {
		qsort(p, length, sizeof(unsigned int), index_cmp);
	}
	else if (length == 2 && p[0] > p[1]) {
		p[1] = key[0];
		p[0] = key[1];
	}

	return p;
}

/// Implementation of a hash map (Word_t key[]) => Word_t
class MapHSOrd : public MapHS {
public:
	MapHSOrd() : MapHS() {
	}

	virtual ~MapHSOrd() {
	}

	/// Lookup for item with key \c key.
	/// \param[in] key Indices that build up the key.
	/// \param[in] length Number of indices.
	/// \param[in] item Item to insert
	/// \return
	/// 	\li true if the key exists, item then contains the value,
	/// 	\li false otherwise
	bool lookup(unsigned int *key, int length, unsigned int &item) const {
		unsigned int *p = index(key, length);
		bool ret = MapHS::lookup((uint8_t *) p, length * sizeof(unsigned int), item);
		delete [] p;

		return ret;
	}

	/// Add a new (key, item) pair into the map
	/// \param[in] key Indices that build up the key.
	/// \param[in] length Number of indices.
	/// \param[in] item Item to insert
	bool set(unsigned int *key, int length, unsigned int item) {
		unsigned int *p = index(key, length);
		bool ret = MapHS::set((uint8_t *) p, length * sizeof(unsigned int), item);
		delete [] p;

		return ret;
	}

	/// Delete an item with key \c key from the map.
	/// \param[in] key Indices that build up the key.
	/// \param[in] length Number of indices.
	bool remove(unsigned int *key, int length) {
		unsigned int *p = index(key, length);
		bool ret = MapHS::remove((uint8_t *) p, length * sizeof(unsigned int));
		delete [] p;

		return ret;
	}
};

#endif
