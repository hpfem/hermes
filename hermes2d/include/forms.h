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

#include "global.h"
#include "quadrature/quad.h"
#include "function/function.h"
#include "function/exact_solution.h"
#include "mesh/refmap.h"
#include "mesh/traverse.h"

namespace Hermes
{
  namespace Hermes2D
  {
    static const char* ERR_UNDEFINED_NEIGHBORING_ELEMENTS =
      "Neighboring elements are not defined and so are not function traces on their interface. "
      "Did you forget setting H2D_ANY_INNER_EDGE in add_matrix/vector_form?";



#pragma region Geometry
    /// Geometry (coordinates, normals, tangents) of either an element or an edge.
    /// @ingroup inner
    template<typename T>
    class HERMES_API Geom
    {
    public:
      /// Constructor.
      Geom();

      T *x, *y;         ///< Coordinates[in physical domain].
      T *nx, *ny;       ///< Normals[in physical domain] (locally oriented
      ///< to point outside the element). Only for edge
      ///< (undefined for element).
      T *tx, *ty;       ///< Tangents[in physical domain]. Only for edge.
      int id;           ///< ID number of the element (undefined for edge).
      int isurf;        ///< Order number of an edge of the element.

      /// Methods designed for discontinuous functions, return errors here.
      virtual int get_neighbor_marker() const { throw Hermes::Exceptions::Exception(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return -1; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual int get_neighbor_id()     const { throw Hermes::Exceptions::Exception(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return -1; }
      /// Methods designed for discontinuous functions, return errors here.
      virtual T   get_neighbor_diam()   const { throw Hermes::Exceptions::Exception(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); return  T(); }

      /// Element diameter (for edge, diameter of the parent element).
      double get_diam_approximation(int n);
      /// Element area (for edge, area of the parent element).
      double get_area(int n, double* wt);

      /// Virtual destructor allowing deallocation of inherited classes (InterfaceGeom) in polymorphic cases.
      virtual ~Geom() {};

      /// Deallocation.
      virtual void free();
      virtual void free_ord() {};
      int elem_marker;       ///< Element marker (for both volumetric and surface forms).
      int edge_marker;       ///< Edge marker (for surface forms only).

      int orientation;  ///< 0 .... if(nx, ny) is equal to the global normal,
      ///< otherwise 1 (each edge has a unique global normal).
      ///< Only for edge.
    };

    template<>
    class HERMES_API Geom<Hermes::Ord>
    {
    public:
      Geom()
      {
        x[0] = y[0] = tx[0] = ty[0] = nx[0] = ny[0] = diam = area = Hermes::Ord(1);
      }
      Hermes::Ord x[1];
      Hermes::Ord y[1];
      Hermes::Ord tx[1];
      Hermes::Ord ty[1];
      Hermes::Ord nx[1];
      Hermes::Ord ny[1];

      /// Element diameter (for edge, diameter of the parent element).
      Hermes::Ord get_diam_approximation(int n) {
        return this->diam;
      }
      /// Element area (for edge, area of the parent element).
      Hermes::Ord get_area(int n, double* wt) {
        return this->area;
      }

      Hermes::Ord diam;           ///< Element diameter (for edge, diameter of the parent element).
      Hermes::Ord area;           ///< Element area (for edge, area of the parent element).
      int id;           ///< ID number of the element (undefined for edge).
      int isurf;        ///< Order number of an edge of the element.

      int elem_marker;       ///< Element marker (for both volumetric and surface forms).
      int edge_marker;       ///< Edge marker (for surface forms only).

      int orientation;  ///< 0 .... if(nx, ny) is equal to the global normal,
    };

    /// Small class which contains information about the element on the other side of an interface.
    ///
    /// It just appends three new_ parameters to an instance of Geom. During destruction, the wrapped
    /// instance is not touched - it must be destroyed separately. You may call the overriden methods
    /// \c free or \c free_ord in order to do this via the instance of InterfaceGeom.
    ///
    /// @ingroup inner
    template<typename T>
    class HERMES_API InterfaceGeom : public Geom<T>
    {
    public:
      int neighb_id;
      T   neighb_diam;
      int get_neighbor_marker() const;
      int get_neighbor_id()  const;
      T get_neighbor_diam() const;

      /// Constructor.
      InterfaceGeom(Geom<T>* geom, int n_marker, int n_id, T n_diam);

      void free();
      void free_ord();

    private:
      Geom<T>* wrapped_geom;
      int neighb_marker;
      template<typename Scalar> friend class KellyTypeAdapt;
    };

    /// Init element geometry for volumetric integrals.
    HERMES_API Geom<double>* init_geom_vol(RefMap *rm, const int order);
    /// Init element geometry for surface integrals.
    HERMES_API Geom<double>* init_geom_surf(RefMap *rm, int isurf, int marker, const int order, double3*& tan);
#pragma endregion

#pragma region Func
    /// Calculated function values (from the class Function) on an element for assembling.
    /// Internal.
    /// @ingroup inner
    template<typename Scalar>
    class HERMES_API Func
    {
    };

    /// Calculated function values (from the class Function) on an element for assembling.
    /// @ingroup inner
    template<>
    class HERMES_API Func<double>
    {
    public:
      /// Constructor.
      Func();
      /// Constructor.
      /** \param[in] num_gip A number of integration points.
      *  \param[in] num_comps A number of components. */
      Func(int np, int nc);
      union
      {
        double val[H2D_MAX_INTEGRATION_POINTS_COUNT];
        double val0[H2D_MAX_INTEGRATION_POINTS_COUNT];
      };

      union
      {
        double dx[H2D_MAX_INTEGRATION_POINTS_COUNT];
        double val1[H2D_MAX_INTEGRATION_POINTS_COUNT];
      };

      union
      {
        double dy[H2D_MAX_INTEGRATION_POINTS_COUNT];
        double curl[H2D_MAX_INTEGRATION_POINTS_COUNT];
      };

      union
      {
        double laplace[H2D_MAX_INTEGRATION_POINTS_COUNT];
        double div[H2D_MAX_INTEGRATION_POINTS_COUNT];
      };

      /// Number of integration points used by this intance.
      int np;
      /// Number of components. Currently accepted values are 1 (H1, L2 space) and 2 (Hcurl, Hdiv space).
      int nc;
      /// Calculate this -= func for each function expations and each integration point.
      /** \param[in] func A function which is added to *this. A number of integratioN points and a number of component has to match. */
      void subtract(Func<double>* func);
      /// Subtract version specifying just one attribute.
      void subtract(double* attribute, double* other_attribute);
      /// Calculate this += func for each function expations and each integration point.
      /** \param[in] func A function which is added to *this. A number of integratioN points and a number of component has to match. */
      void add(Func<double>* func);
      /// Add version specifying just one attribute.
      void add(double* attribute, double* other_attribute);
    };

    /// Calculated function values (from the class Function) on an element for assembling.
    /// @ingroup inner
    template<>
    class HERMES_API Func<std::complex<double> >
    {
    public:
      /// Constructor.
      Func();
      /// Constructor.
      /** \param[in] num_gip A number of integration points.
      *  \param[in] num_comps A number of components. */
      Func(int np, int nc);

      std::complex<double> val[H2D_MAX_INTEGRATION_POINTS_COUNT];
      std::complex<double> val0[H2D_MAX_INTEGRATION_POINTS_COUNT];

      std::complex<double> dx[H2D_MAX_INTEGRATION_POINTS_COUNT];
      std::complex<double> val1[H2D_MAX_INTEGRATION_POINTS_COUNT];

      std::complex<double> dy[H2D_MAX_INTEGRATION_POINTS_COUNT];
      std::complex<double> curl[H2D_MAX_INTEGRATION_POINTS_COUNT];

      std::complex<double> laplace[H2D_MAX_INTEGRATION_POINTS_COUNT];
      std::complex<double> div[H2D_MAX_INTEGRATION_POINTS_COUNT];

      /// Number of integration points used by this intance.
      int np;
      /// Number of components. Currently accepted values are 1 (H1, L2 space) and 2 (Hcurl, Hdiv space).
      int nc;
      /// Calculate this -= func for each function expations and each integration point.
      /** \param[in] func A function which is added to *this. A number of integratioN points and a number of component has to match. */
      void subtract(Func<std::complex<double> >* func);
      /// Subtract version specifying just one attribute.
      void subtract(std::complex<double> * attribute, std::complex<double> * other_attribute);
      /// Calculate this += func for each function expations and each integration point.
      /** \param[in] func A function which is added to *this. A number of integratioN points and a number of component has to match. */
      void add(Func<std::complex<double> >* func);
      /// Add version specifying just one attribute.
      void add(std::complex<double> * attribute, std::complex<double> * other_attribute);
    };

    template<>
    class HERMES_API Func<Ord>
    {
    public:
      Ord val;
      Ord dx;
      Ord dy;
      Ord laplace;

      Ord val0;
      Ord val1;
      Ord curl;
      Ord div;
      void subtract(Func<Ord>* func);
    };


    /// @ingroup inner
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

      T *val;              ///< Function values. If T == Hermes::Ord and orders vary with direction, this returns max(h_order, v_order).
      T *dx, *dy; ///< First-order partial derivatives.
      T *val_neighbor;              ///< Function values. If T == Hermes::Ord and orders vary with direction, this returns max(h_order, v_order).
      T *dx_neighbor, *dy_neighbor; ///< First-order partial derivatives.

      /// One-component constructor.
      ///
      /// \param[in]  fn                  Function defined either on the central or the neighbor element.
      /// \param[in]  support_on_neighbor True if \c fn is defined on the neighbor element, false if on the central element.
      /// \param[in]  reverse             Same meaning as \c reverse_neighbor_side.
      ///
      DiscontinuousFunc(Func<T>* fn, bool support_on_neighbor, bool reverse = false);

      /// Two-component constructor.
      ///
      /// \param[in]  fn_c                Function defined on the central element.
      /// \param[in]  fn_n                Function defined on the neighbor element.
      /// \param[in]  reverse             Same meaning as \c reverse_neighbor_side.
      ///
      DiscontinuousFunc(Func<T>* fn_c, Func<T>* fn_n, bool reverse = false);

      virtual ~DiscontinuousFunc();
      void free();

      using Func<T>::subtract;
      void subtract(const DiscontinuousFunc<T>& func);

      bool reverse_neighbor_side; ///< True if values from the neighbor have to be retrieved in reverse order
      ///< (when retrieving values on an edge that is oriented differently in both elements).
      static T zero;              ///< Zero value used for the zero-extension.
    };

    template<>
    class HERMES_API DiscontinuousFunc<Ord> : public Func<Ord>
    {
    public:
      Func<Ord> *fn_central;        ///< Central element's component.
      Func<Ord> *fn_neighbor;       ///< Neighbor element's component.

      Ord val;
      Ord dx, dy;
      Ord val_neighbor;
      Ord dx_neighbor, dy_neighbor;

      /// One-component constructor.
      ///
      /// \param[in]  fn                  Function defined either on the central or the neighbor element.
      /// \param[in]  support_on_neighbor True if \c fn is defined on the neighbor element, false if on the central element.
      /// \param[in]  reverse             Same meaning as \c reverse_neighbor_side.
      ///
      DiscontinuousFunc(Func<Ord>* fn, bool support_on_neighbor, bool reverse = false);

      /// Two-component constructor.
      ///
      /// \param[in]  fn_c                Function defined on the central element.
      /// \param[in]  fn_n                Function defined on the neighbor element.
      /// \param[in]  reverse             Same meaning as \c reverse_neighbor_side.
      ///
      DiscontinuousFunc(Func<Ord>* fn_c, Func<Ord>* fn_n, bool reverse = false);

      using Func<Ord>::subtract;
      void subtract(const DiscontinuousFunc<Ord>& func);

      bool reverse_neighbor_side; ///< True if values from the neighbor have to be retrieved in reverse order
      ///< (when retrieving values on an edge that is oriented differently in both elements).
      static Ord zero;              ///< Zero value used for the zero-extension.
    };


    /// Init the function for calculation the integration order.
    HERMES_API Func<Hermes::Ord>* init_fn_ord(const int order);

    /// Init the shape function for the evaluation of the volumetric/surface integral (transformation of values).
    HERMES_API Func<double>* init_fn(PrecalcShapeset *fu, RefMap *rm, const int order);
    /// Init the mesh-function for the evaluation of the volumetric/surface integral.
    template<typename Scalar>
    HERMES_API Func<Scalar>* init_fn(MeshFunction<Scalar>* fu, const int order);


    /// Preallocate the Func (all we need is np & nc).
    template<typename Scalar>
    HERMES_API Func<Scalar>* preallocate_fn(pj_pool_t* memoryPool = nullptr);

    /// Init the shape function for the evaluation of the volumetric/surface integral (transformation of values) - preallocated version.
    HERMES_API void init_fn_preallocated(Func<double>* u, PrecalcShapeset *fu, RefMap *rm, const int order);
    /// Init the mesh-function for the evaluation of the volumetric/surface integral - preallocated version.
    template<typename Scalar>
    HERMES_API void init_fn_preallocated(Func<Scalar>* u, MeshFunction<Scalar>* fu, const int order);
    /// Init UExt function - preallocated version.
    template<typename Scalar>
    HERMES_API void init_fn_preallocated(Func<Scalar>* u, UExtFunction<Scalar>* fu, Func<Scalar>** ext, Func<Scalar>** u_ext, const int order, Geom<double>* geometry, ElementMode2D mode);


    /// Utilities follow
    /// Init zero function
    template<typename Scalar>
    HERMES_API Func<Scalar>* init_zero_fn(ElementMode2D mode, int order, Quad2D* quad_2d = nullptr, int nc = 1);

    /// Init UExt function
    template<typename Scalar>
    HERMES_API Func<Scalar>* init_fn(UExtFunction<Scalar>* fu, Func<Scalar>** ext, Func<Scalar>** u_ext, const int order, Geom<double>* geometry, ElementMode2D mode);
#pragma endregion
  }
}
#endif
