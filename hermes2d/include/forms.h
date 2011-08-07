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

/// \file This file contains definition of classes for form evaluation (hence the name) Func, DiscontinuousFunc, Geom, and InterfaceGeom.

#ifndef __H2D_FORMS_H
#define __H2D_FORMS_H

#include "hermes2d_common_defs.h"
#include "quadrature/quad.h"
#include "function/function.h"
#include "function/solution.h"
#include "mesh/refmap.h"
#include "mesh/traverse.h"
#include <complex>

namespace Hermes
{
  namespace Hermes2D
  {
    static const char* ERR_UNDEFINED_NEIGHBORING_ELEMENTS =
      "Neighboring elements are not defined and so are not function traces on their interface. "
      "Did you forget setting H2D_ANY_INNER_EDGE in add_matrix/vector_form?";

    /// Calculated function values (from the class Function) on an element for assembling.
    template<typename T>
    class Func
    {
    public:
      const int num_gip; ///< Number of integration points used by this intance.
      const int nc;      ///< Number of components. Currently accepted values are 1 (H1, L2 space) and 2 (Hcurl, Hdiv space).
      T *val;            ///< Function values. If T == Hermes::Ord and orders vary with direction, this returns max(h_order, v_order).
      T *dx, *dy;        ///< First-order partial derivatives.
      T *laplace;        ///< Sum of second-order partial derivatives. Enabled by defining H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h.
      T *val0, *val1;    ///< Components of a vector field.
      T *dx0, *dx1;      ///< Components of the gradient of a vector field.
      T *dy0, *dy1;      ///< Components of the gradient of a vector field.
      T *curl;           ///< Components of the curl of a vector field.
      T *div;            ///< Components of the div of a vector field.

      /// Constructor.
      /** \param[in] num_gip A number of integration points.
      *  \param[in] num_comps A number of components. */
      Func(int num_gip, int num_comps);

      /// Calculate this -= func for each function expations and each integration point.
      /** \param[in] func A function which is subtracted from *this. A number of integratioN points and a number of component has to match. */
      void subtract(const Func<T>& func);
      void subtract(T* attribute, T* other_attribute);

      /// Calculate this += func for each function expations and each integration point.
      /** \param[in] func A function which is added to *this. A number of integratioN points and a number of component has to match. */
      void add(const Func<T>& func);
      void add(T* attribute, T* other_attribute);

      /// Dellocates an instance of Func<Ord>
      virtual void free_ord();

      /// Dellocates an instance of Func<double> / Func<std::complex<double> >
      virtual void free_fn();

      /// All deallocation done via free_fn / free_ord.
      /// This is to allow proper destruction of DiscontinuousFunc by applying delete on a Func pointer.
      /// NOTE: An error is raised if the user tries to use a Func object for a discontinuous function.
      /// Alternatively, both Func::get_*_central and Func::get_*_neighbor could return the central values as
      /// expected from a continuous function.
      virtual ~Func() { }; 

      /// Methods designed for discontinuous functions, return errors here.
      virtual T& get_val_central(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual T& get_val_neighbor(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual T& get_dx_central(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual T& get_dx_neighbor(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual T& get_dy_central(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS);  return * new T; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual T& get_dy_neighbor(int k) const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual T& get_laplace_central(int k) { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual T& get_laplace_neighbor(int k) { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return * new T; }
    };

    /** \class DiscontinuousFunc forms.h "src/form/forms.h"
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
    class HERMES_API DiscontinuousFunc : public Func<T>
    {
    public:
      Func<T> *fn_central;        ///< Central element's component.
      Func<T> *fn_neighbor;       ///< Neighbor element's component.

      /// One-component constructor.
      ///
      /// \param[in]  fn                  Function defined either on the central or the neighbor element.
      /// \param[in]  support_on_neighbor True if \c fn is defined on the neighbor element, false if on the central element.
      /// \param[in]  reverse             Same meaning as \c reverse_neighbor_side.
      ///
      DiscontinuousFunc(Func<T>* fn, bool support_on_neighbor = false, bool reverse = false);

      /// Two-component constructor.
      ///
      /// \param[in]  fn_c                Function defined on the central element.
      /// \param[in]  fn_n                Function defined on the neighbor element.
      /// \param[in]  reverse             Same meaning as \c reverse_neighbor_side.
      ///
      DiscontinuousFunc(Func<T>* fn_c, Func<T>* fn_n, bool reverse = false);

      /// Get values on the central element.
      virtual T& get_val_central(int k) const;
      /// Get values on the neighboring element.
      virtual T& get_val_neighbor(int k) const;
      /// Get derivatives wrt. x on the central element.
      virtual T& get_dx_central(int k) const;
      /// Get derivatives wrt. x on the neighboring element.
      virtual T& get_dx_neighbor(int k) const;
      /// Get derivatives wrt. y on the central element.
      virtual T& get_dy_central(int k) const;
      /// Get derivatives wrt. y on the neighboring element.
      virtual T& get_dy_neighbor(int k) const;
      /// Get second derivatives (laplacian) on the central element.
      virtual T& get_laplace_central(int k);
      /// Get second derivatives (laplacian) on the neighboring element.
      virtual T& get_laplace_neighbor(int k);

      void subtract(const DiscontinuousFunc<T>& func);

      /// Default destructor may be used. Deallocation is done using the following functions.
      /// FIXME: This is not safe since it allows calling free_ord in a Func<Scalar> object. Template-specialized
      ///  destructors should be used instead (also in Func).
      virtual void free_fn();

      virtual void free_ord();

    private:
      bool reverse_neighbor_side; ///< True if values from the neighbor have to be retrieved in reverse order
      ///< (when retrieving values on an edge that is oriented differently in both elements).
      static T zero;              ///< Zero value used for the zero-extension.

    };

    /// Geometry (coordinates, normals, tangents) of either an element or an edge.
    template<typename T>
    class HERMES_API Geom
    {
    public:
      int elem_marker;       ///< Element marker (for both volumetric and surface forms).
      int edge_marker;       ///< Edge marker (for surface forms only).
      int id;           ///< ID number of the element (undefined for edge).
      T diam;           ///< Element diameter (for edge, diameter of the parent element).
      T *x, *y;         ///< Coordinates [in physical domain].
      T *nx, *ny;       ///< Normals [in physical domain] (locally oriented
      ///< to point outside the element). Only for edge
      ///< (undefined for element).
      T *tx, *ty;       ///< Tangents [in physical domain]. Only for edge.
      int orientation;  ///< 0 .... if (nx, ny) is equal to the global normal,
      ///< otherwise 1 (each edge has a unique global normal).
      ///< Only for edge.

      /// Constructor.
      Geom();
      
      /// Virtual destructor allowing deallocation of inherited classes (InterfaceGeom) in polymorphic cases.
      virtual ~Geom() {};

      /// Deallocation.
      virtual void free();
      virtual void free_ord() {};

      /// Methods designed for discontinuous functions, return errors here.
      virtual int get_neighbor_marker() const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return -1; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual int get_neighbor_id()     const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return -1; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual T   get_neighbor_diam()   const { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return  T(); }
    };


    /// Small class which contains information about the element on the other side of an interface.
    ///
    /// It just appends three new parameters to an instance of Geom. During destruction, the wrapped
    /// instance is not touched - it must be destroyed separately. You may call the overriden methods 
    /// \c free or \c free_ord in order to do this via the instance of InterfaceGeom.
    ///
    template<typename T>
    class HERMES_API InterfaceGeom : public Geom<T>
    {
    public:
      /// Constructor.
      InterfaceGeom(Geom<T>* geom, int n_marker, int n_id, T n_diam);
      
      int get_neighbor_marker() const;
      int get_neighbor_id()  const;
      T get_neighbor_diam() const;
      
      virtual void free() { wrapped_geom->free(); delete wrapped_geom; }
      virtual void free_ord() { delete wrapped_geom; }

    private:
      int neighb_marker;
      int neighb_id;
      T   neighb_diam;
      Geom<T>* wrapped_geom;
    };

    /// Init element geometry for calculating the integration order.
    HERMES_API Geom<Hermes::Ord>* init_geom_ord();
    /// Init element geometry for volumetric integrals.
    HERMES_API Geom<double>* init_geom_vol(RefMap *rm, const int order);
    /// Init element geometry for surface integrals.
    HERMES_API Geom<double>* init_geom_surf(RefMap *rm, SurfPos* surf_pos, const int order);

    /// Init the function for calculation the integration order.
    HERMES_API Func<Hermes::Ord>* init_fn_ord(const int order);
    /// Init the shape function for the evaluation of the volumetric/surface integral (transformation of values).
    HERMES_API Func<double>* init_fn(PrecalcShapeset *fu, RefMap *rm, const int order);
    /// Init the mesh-function for the evaluation of the volumetric/surface integral.
    template<typename Scalar>
    HERMES_API Func<Scalar>* init_fn(MeshFunction<Scalar>*fu, const int order);
    /// Init the solution for the evaluation of the volumetric/surface integral.
    template<typename Scalar>
    HERMES_API Func<Scalar>* init_fn(Solution<Scalar>*fu, const int order);

    /// User defined data that can go to the bilinear and linear forms.
    /// It also holds arbitraty number of functions, that user can use.
    /// Typically, these functions are solutions from the previous time/iteration levels.
    template<typename T>
    class HERMES_API ExtData
    {
    public:
      int nf;           ///< Number of functions in 'fn' array.
      Func<T>** fn;     ///< Array of pointers to functions.

      /// Constructor.
      ExtData();

      /// Deallocation for numerical types.
      void free();

      /// Deallocation for integration order.
      void free_ord();
    };
  }
}
#endif
