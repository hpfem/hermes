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
/*! \file array.h
\brief File containing primarily the class Array<class TYPE> and LightArray<class TYPE>.
*/
#ifndef __HERMES_COMMON_ARRAY_H
#define __HERMES_COMMON_ARRAY_H

#include <vector>
#include <limits.h>
#include "tcmalloc.h"

#ifndef INVALID_IDX
#define INVALID_IDX      INT_MAX
#endif
namespace Hermes
{
  namespace Hermes2D
  {
    /// \brief A generic, inflatable array.
    ///
    /// This class is a generic dynamic array for storing nodes and elements of a mesh.
    /// All items contained in the array are assigned a unique id number. Internally,
    /// a list of unused items is maintained. Unused items (and their id numbers) are
    /// reused when new items are added to the array. The type 'TYPE' must contain the
    /// members 'id' and 'unused' in order to be usable by this class.
    /// \todo Is this dimension independent?
    template<class TYPE>
    class Array
    {
    protected:
      Hermes::vector<TYPE*> pages; ///< \todo standard array for maximum access speed
      Hermes::vector<int> unused;
      int  size, nitems;
      bool append_only;

      static const int HERMES_PAGE_BITS = 10;
      static const int HERMES_PAGE_SIZE = 1 << HERMES_PAGE_BITS;
      static const int HERMES_PAGE_MASK = HERMES_PAGE_SIZE-1;

    public:

      Array()
      {
        size = nitems = 0;
        append_only = false;
      }

      Array(Array& array) { copy(array); }

      ~Array() { free(); }

      /// Makes this array to hold a copy of another one.
      void copy(const Array& array)
      {
        free();

        pages = array.pages;
        unused = array.unused;
        size = array.size;
        nitems = array.nitems;
        append_only = array.append_only;

        for (unsigned i = 0; i < pages.size(); i++)
        {
          TYPE* new_page = new TYPE[HERMES_PAGE_SIZE];
          memcpy(new_page, pages[i], sizeof(TYPE) * HERMES_PAGE_SIZE);
          pages[i] = new_page;
        }
      }

      /// Removes all elements from the array.
      void free()
      {
        for (unsigned i = 0; i < pages.size(); i++) delete [] pages[i];
        pages.clear();
        unused.clear();
        size = nitems = 0;
      }

      /// Sets or resets the append-only mode. In append-only mode new
      /// elements are only added to the end of the array.
      /// This can be useful eg. when refining all elements of a mesh
      /// in which case the newly added elements must not be processed
      /// again by the for-loop. Normally this option should be off.
      void set_append_only(bool append_only)
      {
        this->append_only = append_only;
      }

      /// Wrapper function for Hermes::vector::add() for compatibility purposes.
      int add(TYPE item)
      {
        TYPE* ptr = this->add();
        *ptr = item;
        return ptr->id;
      }

      /// Adds a new item to the array: either it is appended at the
      /// end or an unused item is reused.
      /// \return A reference to the newly allocated item of the array.
      /// The item is assigned an id and its used flag is set to 1.
      TYPE* add()
      {
        TYPE* item;
        if (unused.empty() || append_only)
        {
          if (!(size & HERMES_PAGE_MASK))
          {
            TYPE* new_page = new TYPE[HERMES_PAGE_SIZE];
            pages.push_back(new_page);
          }
          item = pages[size >> HERMES_PAGE_BITS] + (size & HERMES_PAGE_MASK);
          item->id = size++;
          item->used = 1;
        }
        else
        {
          int id = unused.back();
          unused.pop_back();
          item = pages[id >> HERMES_PAGE_BITS] + (id & HERMES_PAGE_MASK);
          item->used = 1;
        }
        nitems++;
        return item;
      }

      /// Removes the given item from the array, ie., marks it as unused.
      /// Note that the array is never physically shrinked. This should not
      /// be a problem, since meshes tend to grow rather than become smaller.
      /// \param id [in] Item id number.
      void remove(int id)
      {
        assert(id >= 0 && id < size);
        TYPE* item = pages[id >> HERMES_PAGE_BITS] + (id & HERMES_PAGE_MASK);
        assert(item->used);
        item->used = 0;
        unused.push_back(id);
        nitems--;
      }

      // Iterators

      /// Get the first index that is present and is equal to or greater than the passed \c idx.
      /// Typically used to begin an iteration over all indices present in the array.
      /// \param[in] idx Optional, default value \c 0 (finds the first present index).
      /// \return
      /// 	\li First index present in the array that is equal or greater than the passed \c idx (if found),
      /// 	\li \c INVALID_IDX (if not found).
      int first(int idx = 0)
      {
        int index = idx;
        while (get(index).used == false)
        {
          index++;
          if (index >= nitems) return INVALID_IDX;
        }
        return index;
      }

      /// Get the first index that is present and is greater than the passed \c idx.
      /// Typically used to continue an iteration over all indices present in the array.
      /// \param[in] idx Index whose succesor we want to find. Optional, default value \c 0.
      /// \return
      /// 	\li First idx present in the array that is greater than the passed \c idx (if found),
      /// 	\li \c INVALID_IDX (if not found).
      int next(int idx = 0)
      {
        int index = idx + 1;
        while (get(index).used == false)
        {
          index++;
          if (index >= nitems) return INVALID_IDX;
        }
        return index;
      }

      /// Get the last index present in the array that is equal to or less than the passed \c idx.
      /// Typically used to begin a reverse iteration over all indices present in the array.
      /// \param[in] idx Optional, default value <c>(Word_t) -1</c> (finds the last index present in the array).
      /// \return
      ///		\li Last index present in the array that is equal or less than the passed \c idx (if found),
      /// 	\li \c INVALID_IDX (if not found).
      int last(int idx = INT_MAX)
      {
        int index = idx;
        if (index > nitems - 1) index = nitems - 1;
        while (get(index).used == false)
        {
          index--;
          if (index < 0) return INVALID_IDX;
        }
        return index;
      }

      /// Get the last index present in the array that is less than the passed \c idx.
      /// Typically used to continue a reverse iteration over all indices present in the array.
      /// \param[in] idx Index whose predecessor we want to find. Optional, default value <c>(Word_t) -1</c>.
      /// \return
      /// 	\li Last index present in the array that is less than the passed \c idx (if found),
      /// 	\li \c INVALID_IDX (if not found).
      int prev(int idx = INT_MAX)
      {
        int index = idx - 1;
        if (index > nitems - 1) index = nitems - 1;
        while (get(index).used == false)
        {
          index--;
          if (index < 0) return INVALID_IDX;
        }
        return index;
      }

      /// Checks whether an element exists at position idx.
      bool exists(int idx)
      {
        if (get(idx).used == true) return true;
        else return false;
      }

      /// Cleans the array and reserves space for up to 'size' items.
      /// This is a special-purpose function, used for loading the array
      /// from file.
      void force_size(int size)
      {
        free();
        while (size > 0)
        {
          TYPE* new_page = new TYPE[HERMES_PAGE_SIZE];
          memset(new_page, 0, sizeof(TYPE) * HERMES_PAGE_SIZE);
          pages.push_back(new_page);
          size -= HERMES_PAGE_SIZE;
        }
        this->size = pages.size() * HERMES_PAGE_SIZE;
      }

      /// Counts the items in the array and registers unused items.
      /// This is a special-purpose function, used after loading the array
      /// from file.
      void post_load_scan(int start = 0)
      {
        nitems = 0;
        for (int i = start; i < size; i++)
          if (get(i).used) nitems++;
          else unused.push_back(i);
      }

      /// Adds an unused item at the end of the array and skips its ID forever.
      /// This is a special-purpose function used to create empty element slots.
      void skip_slot()
      {
        if (!(size & HERMES_PAGE_MASK))
        {
          TYPE* new_page = new TYPE[HERMES_PAGE_SIZE];
          pages.push_back(new_page);
        }
        TYPE* item = pages[size >> HERMES_PAGE_BITS] + (size & HERMES_PAGE_MASK);
        item->id = size++;
        item->used = 0;
        nitems++;
      }

      int get_size() const { return size; }
      int get_num_items() const { return nitems; }

      TYPE& get(int id) const { return pages[id >> HERMES_PAGE_BITS][id & HERMES_PAGE_MASK]; }
      TYPE& operator[] (int id) const { return get(id); }

    };

    /// \brief A light version of the array.
    /// For ordinal types or pointers, does not provide memory handling.
    template<class TYPE>
    class LightArray
    {
    protected:
      TYPE** pages;
      int page_count;
      bool ** presence;
      unsigned int size;

      const unsigned int page_bits;
      const unsigned int page_size;
      const unsigned int page_mask;

    public:

      LightArray(unsigned int page_bits = 9, unsigned int default_page_count = 512) : page_bits(page_bits), page_size(1 << page_bits), page_mask((1 << page_bits) - 1), page_count(default_page_count)
      {
        size = 0;
        pages = (TYPE**)malloc(page_count * sizeof(TYPE*));
        presence = (bool**)malloc(page_count * sizeof(bool*));

        for(int i = 0; i < page_count; i++)
        {
          pages[i] = (TYPE*)malloc(page_size*sizeof(TYPE));
          presence[i] = (bool*)calloc(page_size, sizeof(bool));
        }
      }

      void clear()
      {
        for(unsigned int i = 0; i < page_count; i++)
        {
          memset(presence[i], 0, page_size * sizeof(bool));
        }
      }

      ~LightArray()
      {
        for(unsigned int i = 0; i < page_count; i++)
        {
          free(pages[i]);
          free(presence[i]);
        }
        free(pages);
        free(presence);
      }

      /// Adds a new item to the array.
      void add(TYPE item, unsigned int id)
      {
        TYPE* temp;

        if(id >= page_count * page_size)
        {
          int new_page_count = page_count + ((id - (page_count * page_size)) / page_size) + 2;

          pages = (TYPE**)realloc(pages, new_page_count * sizeof(TYPE*));
          presence = (bool**)realloc(presence, new_page_count * sizeof(bool*));

          for(int i = page_count; i < new_page_count; i++)
          {
            pages[i] = (TYPE*)malloc(page_size*sizeof(TYPE));
            presence[i] = (bool*)calloc(page_size, sizeof(bool));
          }

          page_count = new_page_count;
        }

        temp = pages[id >> page_bits] + (id & page_mask);
        *temp = item;

        presence[id >> page_bits][id & page_mask] = true;

        if(id >= size)
          size = id + 1;
        return;
      }

      unsigned int get_size() const
      {
        return size;
      }

      /// Checks the id position for presence.
      bool present(unsigned int id) const
      {
        if(id >= size)
          return false;
        return presence[id >> page_bits][id & page_mask];
      }

      /// After successful check for presence, the value can be retrieved.
      TYPE& get(unsigned int id) const
      {
        assert(id < size);
        return pages[id >> page_bits][id & page_mask];
      }

    };
  }
}
#endif
