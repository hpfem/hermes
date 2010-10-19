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

#ifndef __H2D_ARRAY_H
#define __H2D_ARRAY_H

#include "h2d_common.h"
#include <vector>

/// \brief A generic, inflatable array.
///
/// This class is a generic dynamic array for storing nodes and elements of a mesh.
/// All items contained in the array are assigned a unique id number. Internally,
/// a list of unused items is maintained. Unused items (and their id numbers) are
/// reused when new items are added to the array. The type 'T' must contain the
/// members 'id' and 'unused' in order to be usable by this class.
///
template<class T>
class Array
{
protected:
  HERMES_API_USED_STL_VECTOR(T*);
  std::vector<T*>  pages; // todo: standard array for maximum access speed
  HERMES_API_USED_STL_VECTOR(int);
  std::vector<int> unused;
  int  size, nitems;
  bool append_only;

  static const int H2D_PAGE_BITS = 10;
  static const int H2D_PAGE_SIZE = 1 << H2D_PAGE_BITS;
  static const int H2D_PAGE_MASK = H2D_PAGE_SIZE-1;

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
      T* new_page = new T[H2D_PAGE_SIZE];
      memcpy(new_page, pages[i], sizeof(T) * H2D_PAGE_SIZE);
      pages[i] = new_page;
    }
  }

  /// Removes all elements from the array.
  void free()
  {
    for (unsigned i = 0; i < pages.size(); i++)
      delete [] pages[i];
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

  /// Adds a new item to the array: either it is appended at the
  /// end or an unused item is reused.
  /// \return A reference to the newly allocated item of the array.
  /// The item is assigned an id and its used flag is set to 1.
  T* add()
  {
    T* item;
    if (unused.empty() || append_only)
    {
      if (!(size & H2D_PAGE_MASK))
      {
        T* new_page = new T[H2D_PAGE_SIZE];
        pages.push_back(new_page);
      }
      item = pages[size >> H2D_PAGE_BITS] + (size & H2D_PAGE_MASK);
      item->id = size++;
      item->used = 1;
    }
    else
    {
      int id = unused.back();
      unused.pop_back();
      item = pages[id >> H2D_PAGE_BITS] + (id & H2D_PAGE_MASK);
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
    T* item = pages[id >> H2D_PAGE_BITS] + (id & H2D_PAGE_MASK);
    assert(item->used);
    item->used = 0;
    unused.push_back(id);
    nitems--;
  }

  /// Cleans the array and reserves space for up to 'size' items.
  /// This is a special-purpose function, used for loading the array
  /// from file.
  void force_size(int size)
  {
    free();
    while (size > 0)
    {
      T* new_page = new T[H2D_PAGE_SIZE];
      memset(new_page, 0, sizeof(T) * H2D_PAGE_SIZE);
      pages.push_back(new_page);
      size -= H2D_PAGE_SIZE;
    }
    this->size = pages.size() * H2D_PAGE_SIZE;
  }

  /// Counts the items in the array and registers unused items.
  /// This is a special-purpose function, used after loading the array
  /// from file.
  void post_load_scan(int start = 0)
  {
    nitems = 0;
    for (int i = start; i < size; i++)
      if (get_item(i).used)
        nitems++;
      else
        unused.push_back(i);
  }

  /// Adds an unused item at the end of the array and skips its ID forever.
  /// This is a special-purpose function used to create empty element slots.
  void skip_slot()
  {
    if (!(size & H2D_PAGE_MASK))
    {
      T* new_page = new T[H2D_PAGE_SIZE];
      pages.push_back(new_page);
    }
    T* item = pages[size >> H2D_PAGE_BITS] + (size & H2D_PAGE_MASK);
    item->id = size++;
    item->used = 0;
    nitems++;
  }

  int get_size() const { return size; }
  int get_num_items() const { return nitems; }

  T& get_item(int id) const { return pages[id >> H2D_PAGE_BITS][id & H2D_PAGE_MASK]; }
  T& operator[] (int id) const { return get_item(id); }

};



#endif
