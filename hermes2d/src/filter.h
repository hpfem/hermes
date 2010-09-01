// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D. If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_FILTER_H
#define __H2D_FILTER_H

#include "solution.h"

struct UniData;


/// Filter is a general postprocessing class, intended for visualization.
/// The output of Filter is an arbitrary combination of up to three input functions,
/// which usually are Solutions to PDEs, but can also be other Filters.
///
/// (This class cannot be instantiated.)
///
class H2D_API Filter : public MeshFunction
{
public:

  Filter() {};

  Filter(Tuple<MeshFunction*> solutions);
  virtual ~Filter();
	
	void init(Tuple<MeshFunction*> solutions);
  
	virtual void set_quad_2d(Quad2D* quad_2d);
  virtual void set_active_element(Element* e);
  virtual void free();
  virtual void reinit();

  virtual void push_transform(int son);
  virtual void pop_transform();

  virtual void init();

protected:

  int num;
  MeshFunction* sln[10];
  uint64_t sln_sub[10];
  void* tables[10];

  bool unimesh;
  UniData** unidata;

  void copy_base(Filter* flt);

};


/// SimpleFilter is a base class for predefined simple filters (MagFilter, DiffFilter...).
/// The 'simplicity' lies in the fact that only one value per input function can be
/// combined (e.g., not a value and a derivative). If this is not sufficient, a full-fledged
/// filter must be derived from the Filter class (see VonMisesFilter). SimpleFilter is also
/// intended for the user to be able to easily create custom filters only by supplying the
/// combining function.
///
/// The user specifies the combining function, the arguments ('sln1', 'sln2', 'sln3'), and
/// optionally the 'item' for each argument, which can be any of H2D_FN_VAL_0, H2D_FN_DX_0, H2D_FN_DY_0
/// etc.
///
/// SimpleFilter is vector-valued, if at least one input function is vector-valued and
/// both components are specified in 'item', e.g., item1 = H2D_FN_DX (which is H2D_FN_DX_0 | H2D_FN_DX_1).
/// Otherwise it is scalar-valued.
///
class H2D_API SimpleFilter : public Filter
{
public:

  SimpleFilter() {};

  SimpleFilter(void (*filter_fn)(int n, Tuple<scalar*> values, scalar* result),
               Tuple<MeshFunction*> solutions, Tuple<int> items = *(new Tuple<int>));


  virtual scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0);

protected:

  int item[3];

  void (*filter_fn)(int n, Tuple<scalar*>, scalar*);

  void init_components();
  virtual void precalculate(int order, int mask);

};


/// DXDYFilter is a more advanced version of SimpleFilter. It allows combining derivatives
/// of the inputs and also, unlike SimpleFilter, it defines derivatives of the filtered
/// result. The user-supplied combining function has a different format: it takes and must
/// return also the DX and DY values.
///
class H2D_API DXDYFilter : public Filter
{
public:

  // one result (rslt), all inputs and result including derivatives
  typedef void (*filter_fn_)(int n, Tuple<scalar *> values, Tuple<scalar *> dx, Tuple<scalar *> dy, scalar* rslt, scalar* rslt_dx, scalar* rslt_dy);

	DXDYFilter() {};
  DXDYFilter(filter_fn_ filter_fn, Tuple<MeshFunction*> solutions);
	
	void init(filter_fn_ filter_fn, Tuple<MeshFunction*> solutions);
  
		virtual scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0)
  { error("Not implemented yet"); return 0; }

protected:

  filter_fn_ filter_fn;

  void init_components();
  virtual void precalculate(int order, int mask);

};


/// MagFilter takes two functions representing the components of a vector function and
/// calculates the vector magnitude, sqrt(x^2 + y^2).
/// \brief Calculates the magnitude of a vector function.
class H2D_API MagFilter : public SimpleFilter
{
  public: 
		MagFilter() {};
		MagFilter(Tuple<MeshFunction*> solutions, Tuple<int> items = *(new Tuple<int>));
		MagFilter(MeshFunction* sln1, int item1 = H2D_FN_VAL); // for vector-valued sln1
};


/// Calculates the difference of two functions.
class H2D_API DiffFilter : public SimpleFilter
{
  public: DiffFilter(Tuple<MeshFunction*> solutions, Tuple<int> items = *(new Tuple<int>));
};


/// Calculates the sum of two functions.
class H2D_API SumFilter : public SimpleFilter
{
  public: SumFilter(Tuple<MeshFunction*> solutions, Tuple<int> items = *(new Tuple<int>));
};


/// Calculates the square of a function.
class H2D_API SquareFilter : public SimpleFilter
{
  public: SquareFilter(Tuple<MeshFunction*> solutions, Tuple<int> items = *(new Tuple<int>));
};


/// Removes the imaginary part from a function.
class H2D_API RealFilter : public SimpleFilter
{
  public: RealFilter(Tuple<MeshFunction*> solutions, Tuple<int> items = *(new Tuple<int>));
};


/// ImagFilter puts the imaginary part of the input function to the real part of the
/// output, allowing it to be visualized.
class H2D_API ImagFilter : public SimpleFilter
{
  public: ImagFilter(Tuple<MeshFunction*> solutions, Tuple<int> items = *(new Tuple<int>));
};


/// Computes the absolute value of a complex solution.
class H2D_API AbsFilter : public SimpleFilter
{
  public: AbsFilter(Tuple<MeshFunction*> solutions, Tuple<int> items = *(new Tuple<int>));
};

/// Computes the angle of a complex solution.
class H2D_API AngleFilter : public SimpleFilter
{
  public: AngleFilter(Tuple<MeshFunction*> solutions, Tuple<int> items = *(new Tuple<int>));
};


/// VonMisesFilter is a postprocessing filter for visualizing elastic stresses in a body.
/// It calculates the stress tensor and applies the Von Mises equivalent stress formula
/// to obtain the resulting stress measure.
/// \brief Calculates the Von Mises stress.
class H2D_API VonMisesFilter : public Filter
{
public: // TODO: cylindrical coordinates

  VonMisesFilter(Tuple<MeshFunction*> solutions, double lambda, double mu,
                 int cyl = 0, int item1 = H2D_FN_VAL, int item2 = H2D_FN_VAL);

  virtual scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0)
  { error("Not implemented yet"); return 0; }

protected:

  double lambda, mu;
  int cyl, item1, item2;

  virtual void precalculate(int order, int mask);

};


/// Linearization filter for use in nonlinear problems. From one or two previous
/// solution values it extrapolates an estimate of the new one.
/// With adaptive time step: tau_frac = tau_new / tau_old
class H2D_API LinearFilter : public Filter
{
  public: LinearFilter(MeshFunction* old);
          LinearFilter(MeshFunction* older, MeshFunction* old, double tau_frac = 1);

  virtual scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0)
  { error("Not implemented yet"); return 0; }

  protected:

    double tau_frac;

    virtual void precalculate(int order, int mask);
    void init_components();
    virtual void set_active_element(Element* e);


};




/// todo: divergence and curl (vorticity) filtr


#endif

