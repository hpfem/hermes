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

#ifndef __H2D_PRECALC_H
#define __H2D_PRECALC_H

#include "function.h"
#include "shapeset/shapeset.h"


/// \brief Caches precalculated shape function values.
///
/// PrecalcShapeset is a cache of precalculated shape function values.
///
///
class HERMES_API PrecalcShapeset : public RealFunction
{
public:

  /// \brief Constructs a standard (master) precalculated shapeset class.
  /// \param shapeset [in] Pointer to the shapeset to be precalculated.
  PrecalcShapeset(Shapeset* shapeset);

  /// \brief Constructs a slave precalculated shapeset class.
  /// \details The slave instance does not hold any precalculated tables.
  /// Instead, it refers to those contained in the master instance. However,
  /// the slave can have different shape function active, different transform
  /// selected, etc. Slave pss's are used for test functions when calling
  /// bilinear forms, inside Solution so as not to disrupt user's pss, etc.
  /// \param master_pss [in] Master precalculated shapeset pointer.
  PrecalcShapeset(PrecalcShapeset* master_pss);

  /// Destructor.
  virtual ~PrecalcShapeset();

  virtual void set_quad_2d(Quad2D* quad_2d);

  /// \brief Frees all precalculated tables.
  virtual void free();

  /// Ensures subsequent calls to get_active_element() will be returning 'e'.
  /// Switches the class to the appropriate mode (triangle, quad).
  virtual void set_active_element(Element* e);

  /// Virtual function handling overflows. Has to be virtual, because
  /// the necessary iterators in the templated class do not work with GCC.
  virtual void handle_overflow_idx();

  /// Activates a shape function given by its index. The values of the shape function
  /// can then be obtained by setting the required integration rule order by calling
  /// set_quad_order() and after that calling get_values(), get_dx_values(), etc.
  /// \param index [in] Shape index.
  void set_active_shape(int index);

  /// Returns the index of the active shape (can be negative if the shape is constrained).
  int get_active_shape() const { return index; };

  /// Returns a pointer to the shapeset which is being precalculated.
  Shapeset* get_shapeset() const { return shapeset; }

  /// Returns type of space
  ESpaceType get_space_type() const { return shapeset->get_space_type(); }

  /// Internal. Use set_active_element() instead.
  void set_mode(int mode);

  /// For internal use only.
  void set_master_transform()
  {
    assert(master_pss != NULL);
    sub_idx = master_pss->sub_idx;
    top = master_pss->top;
    stack[top] = *(master_pss->ctm);
    ctm = stack + top;
  }

  /// Returns the polynomial order of the active shape function on given edge. 
  virtual int get_edge_fn_order(int edge) { return h2d_make_edge_order(mode, edge, shapeset->get_order(index)); }
  
protected:

  Shapeset* shapeset;

  /// Main structure.
  /// There is a 3-layer structure of the precalculated tables.
  /// The first (the lowest) one is the layer where mapping of integral orders to 
  /// Function::Node takes place. See function.h for details.
  /// The second one is the layer with mapping of sub-element transformation to
  /// a table from the lowest layer.
  /// The highest and most complicated one maps a key formed by
  /// quadrature table selector (0-7), mode of the shape function (triangle/quad),
  /// and shape function index to a table from the middle layer.
  std::map<unsigned int, std::map<uint64_t, std::map<unsigned int, Node*>*>*> tables;

  int mode;
  int index;
  int max_index[2];

  PrecalcShapeset* master_pss;

  bool is_slave() const { return master_pss != NULL; }

  virtual void precalculate(int order, int mask);

  void update_max_index();

  /// Forces a transform without using push_transform() etc.
  /// Used by the Solution class. <b>For internal use only</b>.
  void force_transform(uint64_t sub_idx, Trf* ctm)
  {
    this->sub_idx = sub_idx;
    this->ctm = ctm;
  }

  friend class Solution;
  friend class RefMap;
  
public:
    
    /// Key for caching precalculated shapeset values on transformed elements.
    struct Key
    {
      int index;
      int order;
#ifdef _MSC_VER
      UINT64 sub_idx;
#else
      unsigned int sub_idx;
#endif
      int shapeset_type;
#ifdef _MSC_VER
      Key(int index, int order, UINT64 sub_idx, int shapeset_type)
      {
        this->index = index;
        this->order = order;
        this->sub_idx = sub_idx;
        this->shapeset_type = shapeset_type;
      }
#else
      Key(int index, int order, unsigned int sub_idx, int shapeset_type)
      {
        this->index = index;
        this->order = order;
        this->sub_idx = sub_idx;
        this->shapeset_type = shapeset_type;
      }
#endif
    };
    
    /// Functor that compares two PrecalcShapeset keys (needed e.g. to create a std::map indexed by these keys);
    struct Compare
    {
      bool operator()(Key a, Key b) const
      {
        if (a.index < b.index) return true;
        else if (a.index > b.index) return false;
        else
        {
          if (a.order < b.order) return true;
          else if (a.order > b.order) return false;
          else
          {
            if (a.sub_idx < b.sub_idx) return true;
            else if (a.sub_idx > b.sub_idx) return false;
            else
            {
              if (a.shapeset_type < b.shapeset_type) return true;
              else return false;
            }
          }
        }
      }
    };

};



#endif
