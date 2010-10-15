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

#ifndef _MAPINT_H_
#define _MAPINT_H_

#include <assert.h>

/// \file mapint.h
/// \brief Class for a dynamic array, using an array-of-bytes of Length as an Index and a word as a Value.

#include <Judy.h>

#ifndef JUDY_INVALID
#define JUDY_INVALID					((Word_t) -1)
#endif

/// \class MapInt
/// \brief Implementation of a dynamic array, using an array-of-bytes of Length as an Index and a word as a Value.
/// C++ encapsulation of JudyHS functions.
template<class TYPE>
class MapInt {
protected:
	void *judy;

public:
	MapInt() {
		judy = NULL;
	}

	virtual ~MapInt() {
		remove_all();
	}

	/// Lookup for item with key \c key.
	/// \param[in] key Pointer to the array-of-bytes.
	/// \param[in] item Item to insert
	/// \return
	/// 	\li true if the key exists, item then contains the value,
	/// 	\li false otherwise
	bool lookup(const char *key, TYPE &item) {
		void *pval;
		// check if the key exists
		JSLG(pval, judy, (const uint8_t *) key);
		if (pval == NULL) {
			return false;
		}
		else {
			item = **(TYPE **) pval;
			return true;
		}
	}

	/// Add a new (key, item) pair into the map
	/// \param[in] key Pointer to the array-of-bytes.
	/// \param[in] item Item to insert
	bool set(const char *key, TYPE item) {
		void *pval;
		JSLG(pval, judy, (const uint8_t *) key);
		if (pval == NULL) {
			// store the key -> value association
			JSLI(pval, judy, (const uint8_t *) key);
			if (pval == NULL) return false;
			*(TYPE **) pval = new TYPE;
		}
		**(TYPE **) pval = item;
		return true;
	}

	/// Delete an item with key \c key from the map.
	/// \param[in] key Pointer to the array-of-bytes.
	bool remove(char *key) {
		void *pval;
		JSLG(pval, judy, (const uint8_t *) key);
		if (pval == NULL) {
			return false;
		}
		else {
			delete *(TYPE **) pval;

			int rc;
			JSLD(rc, judy, (const uint8_t *) key);
			return (rc == 1);
		}
	}

	/// Remove all items from the array
	void remove_all() {
		int val;
		JHSFA(val, judy);
	}

	// Iterators

#define MAXLINELEN 1000000

	/// Get the first index that is present and is equal to or greater than the passed \c idx.
	/// Typically used to begin an iteration over all indices present in the array.
	/// \param[in] idx Optional, default value \c 0 (finds the first present index).
	/// \return
	/// 	\li First index present in the array that is equal or greater than the passed \c idx (if found),
	/// 	\li \c JUDY_INVALID (if not found).
	char *first() const {
//		return NULL;
		void *pval = NULL;
		char index[MAXLINELEN];
		strcpy(index, "");
		JSLF(pval, judy, (uint8_t *) index);
		return pval ? index : NULL;
	}

	/// Get the first index that is present and is greater than the passed \c idx.
	/// Typically used to continue an iteration over all indices present in the array.
	/// \param[in] idx Index whose succesor we want to find. Optional, default value \c 0.
	/// \return
	/// 	\li First idx present in the array that is greater than the passed \c idx (if found),
	/// 	\li \c JUDY_INVALID (if not found).
	char *next(char *index) const {
		void *pval;
		JSLF(pval, judy, (uint8_t *) index);
		return pval ? index : NULL;
	}

	/// Get the last index present in the array that is equal to or less than the passed \c idx.
	/// Typically used to begin a reverse iteration over all indices present in the array.
	/// \param[in] idx Optional, default value <c>(Word_t) -1</c> (finds the last index present in the array).
	/// \return
	///		\li Last index present in the array that is equal or less than the passed \c idx (if found),
	/// 	\li \c JUDY_INVALID (if not found).
	Word_t last(Word_t idx = (Word_t) -1) const;

	/// Get the last index present in the array that is less than the passed \c idx.
	/// Typically used to continue a reverse iteration over all indices present in the array.
	/// \param[in] idx Index whose predecessor we want to find. Optional, default value <c>(Word_t) -1</c>.
	/// \return
	/// 	\li Last index present in the array that is less than the passed \c idx (if found),
	/// 	\li \c JUDY_INVALID (if not found).
	Word_t prev(Word_t idx = (Word_t) -1) const;
};

#endif
