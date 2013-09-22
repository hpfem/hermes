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

#include "common.h"

#ifndef INVALID_IDX
#define INVALID_IDX      INT_MAX
#endif
namespace Hermes
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
    TYPE** pages; ///< \todo standard array for maximum access speed
    int* unused;
    int  page_count, size, nitems, unused_size, nunused;
    bool append_only;

    static const int HERMES_PAGE_BITS = 8;
    static const int HERMES_PAGE_SIZE = 1 << HERMES_PAGE_BITS;
    static const int HERMES_PAGE_MASK = HERMES_PAGE_SIZE-1;
  public:

    Array(int initial_page_count = 0)
    {
      size = nitems = nunused = 0;
      page_count = initial_page_count;
      unused_size = initial_page_count;

      if(page_count)
      {
        this->pages = (TYPE**)malloc(sizeof(TYPE*) * page_count);
        for (unsigned i = 0; i < this->page_count; i++)
          this->pages[i] = new TYPE[HERMES_PAGE_SIZE];
        this->unused = (int*)malloc(sizeof(int));
      }
      else
      {
        this->pages = NULL;
        this->unused = NULL;
      }
      append_only = false;
    }

    Array(Array& array) { copy(array); }

    ~Array()
    {
      free();
      ::free(this->pages);
      ::free(this->unused);
    }

    /// Makes this array to hold a copy of another one.
    void copy(const Array& array)
    {
      free();

      this->pages = (TYPE**)realloc(this->pages, array.page_count * sizeof(TYPE*));
      this->unused = (int*)realloc(this->unused, array.unused_size * sizeof(int));
      memcpy(this->unused, array.unused, array.unused_size * sizeof(int));
      
      this->page_count = array.page_count;
      this->size = array.size;
      this->nitems = array.nitems;
      this->unused_size = array.unused_size;
      this->nunused = array.nunused;
      this->append_only = array.append_only;

      for (unsigned i = 0; i < this->page_count; i++)
      {
        TYPE* new_page = new TYPE[HERMES_PAGE_SIZE];
        memcpy(new_page, array.pages[i], sizeof(TYPE) * HERMES_PAGE_SIZE);
        this->pages[i] = new_page;
      }
    }

    /// Removes all elements from the array.
    void free()
    {
      for (unsigned i = 0; i < this->page_count; i++)
        delete [] pages[i];
      size = nitems = nunused = page_count = unused_size = 0;
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
      if (!nunused || append_only)
      {
        if (!(size & HERMES_PAGE_MASK))
        {
          this->pages = (TYPE**)realloc(this->pages, (this->page_count + 1) * sizeof(TYPE*));
          TYPE* new_page = new TYPE[HERMES_PAGE_SIZE];
          pages[this->page_count++] = new_page;
        }
        item = pages[size >> HERMES_PAGE_BITS] + (size & HERMES_PAGE_MASK);
        item->id = size++;
        item->used = 1;
      }
      else
      {
        int id = unused[nunused-- - 1];
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
      TYPE* item = pages[id >> HERMES_PAGE_BITS] + (id & HERMES_PAGE_MASK);
      item->used = 0;
      if(nunused >= unused_size)
      {
        this->unused = (int*)realloc(this->unused, ++this->unused_size * sizeof(int));
      }
      unused[nunused++] = id;
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

    /// Adds an unused item at the end of the array and skips its ID forever.
    /// This is a special-purpose function used to create empty element slots.
    void skip_slot()
    {
      if (!(size & HERMES_PAGE_MASK))
      {
        this->pages = (TYPE**)realloc(this->pages, (this->page_count + 1) * sizeof(TYPE*));
        TYPE* new_page = new TYPE[HERMES_PAGE_SIZE];
        pages[this->page_count++] = new_page;
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
      return pages[id >> page_bits][id & page_mask];
    }

  };
}
#endif
