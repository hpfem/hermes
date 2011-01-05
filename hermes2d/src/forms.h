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


#ifndef __H2D_FORMS_H
#define __H2D_FORMS_H

#include "h2d_common.h"
#include "quad.h"
#include "function.h"
#include "solution.h"
#include "refmap.h"

#define callback(a)	a<double, scalar>, a<Ord, Ord>

/// Base type for orders of functions.
///
/// We defined a special arithmetics with this type to be able to analyze forms
/// and determine the necessary integration order.  This works for forms, but it also
/// works for user-defined functions.
class Ord
{
public:

  Ord(): order(0) {}
  explicit Ord(int o): order(o) {}
  Ord(double d): order(0) {}

  int get_order() const { return order; }
  int get_max_order() const {return 30;}

  Ord operator+(const Ord &o) { return Ord(std::max(this->order, o.order)); }
  Ord operator+(double d) { return *this; }
  Ord operator-(const Ord &o) { return Ord(std::max(this->order, o.order)); }
  Ord operator-(double d) { return *this; }
  Ord operator*(const Ord &o) { return Ord(this->order + o.order); }
  Ord operator*(double d) { return *this; }
  Ord operator/(const Ord &o) { return Ord(this->get_max_order()); }
  Ord operator/(double d) { return *this; }

  Ord operator+=(const Ord &o) { this->order = std::max(this->order, o.order); return *this; }

  Ord operator+=(const double &d) { return *this; }
  Ord operator-=(const double &d) { return *this; }
  Ord operator*=(const double &d) { return *this; }
  Ord operator/=(const double &d) { return *this; }

  bool operator<(double d) { return true; }
  bool operator>(double d) { return false; }

protected:
  int order;

};

inline Ord operator/(const scalar &a, const Ord &b) { return Ord(b.get_max_order()); }
inline Ord operator*(const scalar &a, const Ord &b) { return b; }
inline Ord operator+(const scalar &a, const Ord &b) { return b; }
inline Ord operator-(const scalar &a, const Ord &b) { return b; }
inline Ord operator-(const Ord &a) { return a; }

inline Ord pow(const Ord &a, const double &b) { return Ord((int) ceil(fabs(b)) * a.get_order()); }
inline Ord sqrt(const Ord &a) { return a; }
inline Ord sqr(const Ord &a) { return Ord(2 * a.get_order()); }
inline Ord conj(const Ord &a) { return a; }
inline Ord abs(const Ord &a) { return a; }

inline Ord atan2(const Ord &a, const Ord &b) { return Ord(a.get_max_order()); }
inline Ord atan(const Ord &a) { return Ord(a.get_max_order()); }
inline Ord sin(const Ord &a) { return Ord(a.get_max_order()); }
inline Ord cos(const Ord &a) { return Ord(a.get_max_order()); }
inline Ord log(const Ord &a) { return Ord(a.get_max_order()); }
inline Ord exp(const Ord &a) { return Ord(3 * a.get_order()); }

static const char* ERR_UNDEFINED_NEIGHBORING_ELEMENTS = 
  "Neighboring elements are not defined and so are not function traces on their interface. "
  "Did you forget setting H2D_ANY_INNER_EDGE in add_matrix/vector_form?";
  
// Function
template<typename T>
class Func
{
public:
  const int num_gip; ///< Number of integration points used by this intance.
  const int nc;      ///< Number of components. Currently accepted values are 1 (H1, L2 space) and 2 (Hcurl, Hdiv space).
  T *val;            ///< Function values. If T == Ord and orders vary with direction, this returns max(h_order, v_order).
  T *dx, *dy;        ///< First-order partial derivatives.
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  T *laplace;        ///< Sum of second-order partial derivatives. Enabled by defining H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h.
#endif
  T *val0, *val1;    ///< Components of a vector field.
  T *dx0, *dx1;      ///< Components of the gradient of a vector field.
  T *dy0, *dy1;

  T *curl;           ///< Components of the curl of a vector field.
  T *div;            ///< Components of the div of a vector field.

  /// Constructor.
  /** \param[in] num_gip A number of integration points.
   *  \param[in] num_comps A number of components. */
  explicit Func(int num_gip, int num_comps) : num_gip(num_gip), nc(num_comps) {
    val = val0 = val1 = NULL;
    dx = dx0 = dx1 = NULL;
    dy = dy0 = dy1 = NULL;
    curl = NULL;
    div = NULL;
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    laplace = NULL;
#endif
  };

/// Subtract arrays stored in a given attribute from the same array in provided function.
#define H2D_SUBTRACT_IF_NOT_NULL(__ATTRIB, __OTHER_FUNC) { if (__ATTRIB != NULL) { \
  assert_msg(__OTHER_FUNC.__ATTRIB != NULL, "Unable to subtract a function expansion " #__ATTRIB " is NULL in the other function."); \
  for(int i = 0; i < num_gip; i++) __ATTRIB[i] -= __OTHER_FUNC.__ATTRIB[i]; } }

  /// Calculate this -= func for each function expations and each integration point.
  /** \param[in] func A function which is subtracted from *this. A number of integratio points and a number of component has to match. */
  //FIXME : It should be 'virtual', but then it doesn't compile.
  void subtract(const Func<T>& func) {
    assert_msg(num_gip == func.num_gip, "Unable to subtract a function due to a different number of integration points (this: %d, other: %d)", num_gip, func.num_gip);
    assert_msg(nc == func.nc, "Unable to subtract a function due to a different number of components (this: %d, other: %d)", nc, func.nc);
    H2D_SUBTRACT_IF_NOT_NULL(val, func)
    H2D_SUBTRACT_IF_NOT_NULL(dx, func)
    H2D_SUBTRACT_IF_NOT_NULL(dy, func)
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    H2D_SUBTRACT_IF_NOT_NULL(laplace, func)
#endif
    if (nc > 1) {
      H2D_SUBTRACT_IF_NOT_NULL(val0, func)
      H2D_SUBTRACT_IF_NOT_NULL(val1, func)
      H2D_SUBTRACT_IF_NOT_NULL(dx0, func)
      H2D_SUBTRACT_IF_NOT_NULL(dx1, func)
      H2D_SUBTRACT_IF_NOT_NULL(dy0, func)
      H2D_SUBTRACT_IF_NOT_NULL(dy1, func)
      H2D_SUBTRACT_IF_NOT_NULL(curl, func)
      H2D_SUBTRACT_IF_NOT_NULL(div, func)
    }
  };
#undef H2D_SUBTRACT_IF_NOT_NULL

  virtual void free_ord() {
    delete val;
    val = val0 = val1 = NULL;
    dx = dx0 = dx1 = NULL;
    dy = dy0 = dy1 = NULL;
    curl = NULL;
    div = NULL;
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    laplace = NULL;
#endif
  }
  virtual void free_fn()
  {
    delete [] val; val = NULL;
    delete [] dx; dx = NULL;
    delete [] dy; dy = NULL;
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    delete [] laplace; laplace = NULL;
#endif

    delete [] val0; delete [] val1; val0 = val1 = NULL;
    delete [] dx0;  delete [] dx1; dx0 = dx1 = NULL;
    delete [] dy0;  delete [] dy1; dy0 = dy1 = NULL;
    delete [] curl; curl = NULL;
    delete [] div; div = NULL;
  }
  
  virtual ~Func() { }; // All deallocation done via free_fn / free_ord. 
                       // This is to allow proper destruction of DiscontinuousFunc by applying delete on a Func pointer.
  
  // NOTE: An error is raised if the user tries to use a Func object for a discontinuous function.
  // Alternatively, both Func::get_*_central and Func::get_*_neighbor could return the central values as
  // expected from a continuous function. 
  virtual T& get_val_central(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; } 
  virtual T& get_val_neighbor(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }
  virtual T& get_dx_central(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }
  virtual T& get_dx_neighbor(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }
  virtual T& get_dy_central(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS);  return * new T; }
  virtual T& get_dy_neighbor(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }

};

#define GET_CENT_FN(__ATTRIB) return (fn_central != NULL) ? fn_central->__ATTRIB[k] : zero;
#define GET_NEIB_FN(__ATTRIB) return (fn_neighbor != NULL) ? fn_neighbor->__ATTRIB[ reverse_neighbor_side ? fn_neighbor->num_gip-k-1 : k ] : zero;

/** \class DiscontinuousFunc forms.h "src/forms.h"
 *  \brief This class represents a function with jump discontinuity on an interface of two elements.
 *
 *  We will refer to one of the elements sharing the interface of discontinuity as to the \em central element,
 *  while to the other one as to the \em neighbor element.
 *
 *  Instance of the class may be constructed either with two \c Func objects, which represent the continuous 
 *  components on the central and the neighbor element, respectively, or with only one \c Func object and
 *  information about its support (where it attains non-zero value). The discontinuous function is in the latter
 *  case constructed by extending the supplied function by zero to the other element. Values and derivatives from 
 *  both elements may then be obtained by querying the corresponding \c Func object, using methods
 *  \c get_val_central, \c get_val_neighbor, etc.
 **/
template<typename T>
class DiscontinuousFunc : public Func<T>
{
private:
  bool reverse_neighbor_side; ///< True if values from the neighbor have to be retrieved in reverse order
                              ///< (when retrieving values on an edge that is oriented differently in both elements).
  static T zero;              ///< Zero value used for the zero-extension.
  
public: 
  Func<T> *fn_central;        ///< Central element's component.
  Func<T> *fn_neighbor;       ///< Neighbor element's component.
  
  /// One-component constructor.
  ///
  /// \param[in]  fn                  Function defined either on the central or the neighbor element.
  /// \param[in]  support_on_neighbor True if \c fn is defined on the neighbor element, false if on the central element.
  /// \param[in]  reverse             Same meaning as \c reverse_neighbor_side.
  ///
  DiscontinuousFunc(Func<T>* fn, bool support_on_neighbor = false, bool reverse = false) :
    Func<T>(fn->num_gip, fn->nc),
    reverse_neighbor_side(reverse),
    fn_central(NULL),
    fn_neighbor(NULL)
  { 
    assert_msg(fn != NULL, "Invalid arguments to DiscontinuousFunc constructor.");
    if (support_on_neighbor) fn_neighbor = fn; else fn_central = fn;
  }
  
  /// Two-component constructor.
  ///
  /// \param[in]  fn_c                Function defined on the central element.
  /// \param[in]  fn_n                Function defined on the neighbor element.
  /// \param[in]  reverse             Same meaning as \c reverse_neighbor_side.
  ///
  DiscontinuousFunc(Func<T>* fn_c, Func<T>* fn_n, bool reverse = false) : 
    Func<T>(fn_c->num_gip, fn_c->nc),
    reverse_neighbor_side(reverse),
    fn_central(fn_c),
    fn_neighbor(fn_n)
  { 
    assert_msg(fn_c != NULL && fn_n != NULL, "Invalid arguments to DiscontinuousFunc constructor.");
    assert_msg(fn_c->num_gip == fn_n->num_gip && fn_c->nc == fn_n->nc,
                "DiscontinuousFunc must be formed by two Func's with same number of integration points and components.");
  }
  
  // Get values, derivatives, etc. in both elements adjacent to the discontinuity.
  
  virtual T& get_val_central(int k) const { GET_CENT_FN(val); }
  virtual T& get_val_neighbor(int k) const { GET_NEIB_FN(val); }
  virtual T& get_dx_central(int k) const { GET_CENT_FN(dx); }
  virtual T& get_dx_neighbor(int k) const { GET_NEIB_FN(dx); }
  virtual T& get_dy_central(int k) const { GET_CENT_FN(dy); }
  virtual T& get_dy_neighbor(int k) const { GET_NEIB_FN(dy); }
  
  #ifdef H2D_SECOND_DERIVATIVES_ENABLED
    virtual T& get_laplace_central(int k) { GET_CENT_FN(laplace); }
    virtual T& get_laplace_neighbor(int k) { GET_NEIB_FN(laplace); }
  #endif
  
  void subtract(const DiscontinuousFunc<T>& func) 
  { 
    // TODO: Add sanity checks, revise for adaptivity.
    if (fn_central != NULL && func.fn_central != NULL)
      fn_central->subtract(func.fn_central);
    if (fn_neighbor != NULL && func.fn_neighbor != NULL)
      fn_neighbor->subtract(func.fn_neighbor);
  }
  
  // Default destructor may be used. Deallocation is done using the following functions.
  // FIXME: This is not safe since it allows calling free_ord in a Func<scalar> object. Template-specialized
  //  destructors should be used instead (also in Func).
  virtual void free_fn() { 
    if (fn_central != NULL) {
      fn_central->free_fn(); 
      delete fn_central;
      fn_central = NULL;
    }
    if (fn_neighbor != NULL) {
      fn_neighbor->free_fn();
      delete fn_neighbor;
      fn_neighbor = NULL;
    }
  }
  virtual void free_ord() {
    if (fn_central != NULL) {
      fn_central->free_ord(); 
      delete fn_central;
      fn_central = NULL;
    }
    if (fn_neighbor != NULL) {
      fn_neighbor->free_ord();
      delete fn_neighbor;
      fn_neighbor = NULL;
    }
  }
};


/// Geometry (coordinates, normals, tangents) of either an element or an edge.
template<typename T>
class Geom
{
public:
  int elem_marker;       ///< Element marker (for both volumetric and surface forms).
  int edge_marker;       ///< Edge marker (for surface forms only).
  int id;           ///< ID number of the element (undefined for edge).
  T diam;           ///< Element diameter (for edge, diameter of the parent element).
  //Element *element;   // Active element. NOTE: We used this for some time but
                        // decided against it because (a) it disables automatic
                        // order parsing and (b) if the form is called with T
                        // == Ord, element is not initialized, so the user has
                        // to be aware of this and test it in his weak form.

  T *x, *y;         ///< Coordinates [in physical domain].
  T *nx, *ny;       ///< Normals [in physical domain] (locally oriented
                    ///< to point outside the element). Only for edge 
                    ///< (undefined for element).
  T *tx, *ty;       ///< Tangents [in physical domain]. Only for edge.
  int orientation;  ///< 0 .... if (nx, ny) is equal to the global normal,
                    ///< otherwise 1 (each edge has a unique global normal).
                    ///< Only for edge.

  Geom()
  {
    elem_marker = -1;
    edge_marker = -1;
    id = 0;
    x = y = NULL;
    nx = ny = NULL;
    tx = ty = NULL;
    diam = 0;
  }
    
  virtual ~Geom() { };

  void free()
  {
    delete [] tx;    delete [] ty;
    delete [] nx;    delete [] ny;
  }
  
  virtual int get_neighbor_marker() const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return -1; }
  virtual int get_neighbor_id()     const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return -1; }
  virtual T   get_neighbor_diam()   const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return  T(); }
};


/// Small class which contains information about the element on the other side of an interface.
///
/// It just appends three new parameters to an instance of Geom. During destruction, the wrapped
/// instance is also automatically destroyed.
///
template<typename T>
class InterfaceGeom : public Geom<T>
{
private:  
  const Geom<T>* central_geom;  // The wrapped instance of Geom (representing geometric data for the
                                // central element).
  
public:
  int neighb_marker;
  int neighb_id;
  T   neighb_diam;
  
  InterfaceGeom(const Geom<T>* geom, int n_marker, int n_id, T n_diam) : 
      Geom<T>(), central_geom(geom), neighb_marker(n_marker), neighb_id(n_id), neighb_diam(n_diam)
  {
    // Let this class expose the standard Geom interface.
    this->edge_marker = geom->edge_marker;
    this->elem_marker = geom->elem_marker;
    this->id = geom->id;
    this->diam = geom->diam;
    this->x = geom->x;
    this->y = geom->y;
    this->tx = geom->tx; 
    this->ty = geom->ty; 
    this->nx = geom->nx; 
    this->ny = geom->ny; 
    this->orientation = geom->orientation;
  }
  
  ~InterfaceGeom() { delete central_geom; }
  
  virtual int get_neighbor_marker() const { return neighb_marker; }
  virtual int get_neighbor_id()     const { return neighb_id; }
  virtual T   get_neighbor_diam()   const { return neighb_diam; }
};

/// Init element geometry for calculating the integration order.
Geom<Ord>* init_geom_ord();
/// Init element geometry for volumetric integrals.
Geom<double>* init_geom_vol(RefMap *rm, const int order);
/// Init element geometry for surface integrals.
Geom<double>* init_geom_surf(RefMap *rm, SurfPos* surf_pos, const int order);


/// Init the function for calculation the integration order.
Func<Ord>* init_fn_ord(const int order);
/// Init the shape function for the evaluation of the volumetric/surface integral (transformation of values).
Func<double>* init_fn(PrecalcShapeset *fu, RefMap *rm, const int order);
/// Init the mesh-function for the evaluation of the volumetric/surface integral.
Func<scalar>* init_fn(MeshFunction *fu, RefMap *rm, const int order);



/// User defined data that can go to the bilinear and linear forms.
/// It also holds arbitraty number of functions, that user can use.
/// Typically, these functions are solutions from the previous time/iteration levels.
template<typename T>
class ExtData 
{
public:
  int nf;         ///< Number of functions in 'fn' array.
  Func<T>** fn;   ///< Array of pointers to functions.

  ExtData() 
  {
    nf = 0;
    fn = NULL;
  }

  void free()
  {
    for (int i = 0; i < nf; i++)
    {
      fn[i]->free_fn();
      delete fn[i];
    }
    delete [] fn;
  }

  void free_ord()
  {
    for (int i = 0; i < nf; i++)
    {
      fn[i]->free_ord();
      delete fn[i];
    }
    delete [] fn;
  }
};

#endif
