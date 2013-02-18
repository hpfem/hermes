Hermes C++ object model - deriving your own specialized classes
---------------------------------------------------------------

There are several classes that represent some piece of the whole FEM discretization and calculation process that hold the custom informaion for a specific problem.

Be sure to have checked the section "Typical example" that will give you an idea of the basics. This section is more technically focused.

These are

 * Forms (= integrals) of the weak formulation
 * Essential (= Dirichlet) boundary conditions
 * Mesh functions (= initial conditions, exact solutions, etc.)
 * Mathematical functions (= nonlinear relations, etc.)
 * Filters (= functions of solutions, for post-processing)
 
Each of these is described in what follows.

Notes on templating
~~~~~~~~~~~~~~~~~~~
In the section of advanced Object-Oriented aspects of C++ and especially its template metaprogramming capabilities, you should find enough information on the subject to start with, should you need.
The templates in Hermes are very simple. Most classes (and all that are discussed in this section) come only in two forms - there exist only two *template instantiations* - for real and for complex numbers.
Also, since we do not see much use in using single-precision calculations, only the two following instantiations exist for most classes:

  * double
  * std::complex<double>
  
So whenever there is a class header that looks like this::
  
  template<typename Scalar>
  class AClass
  {
    ...
  };

It is exactly this case, you can use AClass<double> for real problems, and AClass<std::complex<double> > for complex problems.

The explicit instantiations, should you need to add such, need to be at the end of the **source** file (.cpp), **not** the header file.

Forms
~~~~~~~~
Weak formulation is the key part of the FEM discretization.
The whole weak formulation with all the integrals that it contains and their specification (over what part of the domain is the integration, on which coordinates of the system of PDEs it is located, etc.) is represented by the class template

  WeakForm<Scalar> (file: hermes2d/include/weakform/weakform.h)

In the implementation of your programs, you need to subclass (or derive from) this class as was seen in the "Weak formulation" part of the "Typical example" section:
  
  class MyWeakForm : public WeakForm<double>
  
The important methods and attributes to note in this class are:

  - the constructor WeakForm(unsigned int neq = 1, bool mat_free = false)
  
    - neq stands for Number Of Equations
    - in the body, you need to use the following methods to add single forms (integrals) to the weak formulation
    
      - void add_matrix_form(MatrixFormVol<Scalar>* mfv);
      - void add_matrix_form_surf(MatrixFormSurf<Scalar>* mfs);
      - ...
    - other methods (set_ext(...), set_current_time(...)) are described in the Doxygen documentation.
    - **clone()** - the cloning method, for parallelization
    
      - the purpose of this method is to create an exact copy of the instance at hand (as simple as that - return a copy, but with no shared data)
      - the reason is that in the openMP parallel assembling, each thread will receive a copy of the instance for its own processing (and this will eliminate the potential parallelization issues)
      - this is there just for the case that you have any data that are not thread-safe and you do not want to take care about this on your own (using openMP directives like **pragma omp critical** etc.)
    
Right after the WeakForm come the most difficutlt and problem causing classes - **Forms** - here is the list of form types:

  * MatrixFormVol - a bilinear form that represent integration over a volume
  * MatrixFormSurf - a bilinear form that represent integration over a surface
  * MatrixFormDG - a bilinear form that represent integration over an inner edge - for DG
  * VectorFormVol - a linear form that represent integration over a volume
  * VectorFormSurf - a linear form that represent integration over a surface
  * VectorFormDG - a linear form that represent integration over an inner edge - for DG
  
How to subclass these should be obvious from the "Weak formulation" part of the "Typical example" section. 
  
Essential boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There is three classes to note for handling of Dirichlet boundary conditions (all in hermes2d/include/boundary_conditions/essential_boundary_conditions.h)

  * DefaultEssentialBCConst - a class you may use directly, without any subclassing
  
    * only thing is to specify the constant value in the constructor::
      
      DefaultEssentialBCConst(std::string marker, Scalar **value_const**);
  * DefaultEssentialBCNonConst - a class you again may use directly, **or** subclass, should you wish to use a non-Constant Dirichlet BC

    * if you wish to use the class directly, you need to specify an instance of the following in the constructor::
    
      ExactSolutionScalar<Scalar>
    * if you wish to use a specific value function, you may just subclass this class and override the value() method
    * instructions for subclassing can be found in the file hermes2d/include/boundary_conditions/essential_boundary_conditions.h just above the class
  
  * EssentialBCs - a class you will never subclass, it is a container that you pass to a constructor of a Space (see the "Space" part of the "Typical example" section).
  
An example usage of the non-constant boundary condition with subclassing is in the test examples::

    hermes2d\test_examples\03-navier-stokes\definitions.cpp at the very bottom

Mesh functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
An example of this is the following code from the test example 06::

  hermes2d\test_examples\06-system-adapt\definitions.cpp (.h)
  
The point here are the two classes ExactSolutionFitzHughNagumo1, ExactSolutionFitzHughNagumo2. This is the definition of the class **ExactSolutionFitzHughNagumo1** and the declaration of its methods in definitions.h::

  class ExactSolutionFitzHughNagumo1 : public ExactSolutionScalar<double>
  {
  public:
    ExactSolutionFitzHughNagumo1(const Mesh* mesh);

    virtual double value(double x, double y) const;

    virtual void derivatives(double x, double y, double& dx, double& dy) const;

    virtual Ord ord(Ord x, Ord y) const;

    ~ExactSolutionFitzHughNagumo1();
    
    virtual MeshFunction<double>* clone() const;

    CustomExactFunction1* cef1;
  };
  
Note the subclassing line::

  // This should be obvious for any C++ user
  class ExactSolutionFitzHughNagumo1 : public ExactSolutionScalar<double>
  
Then we can see the important methods are overriden in the source file definitions.cpp:
  * value(double x, double y) const 

    ::
  
      double ExactSolutionFitzHughNagumo1::value(double x, double y) const 
      {
        return cef1->val(x)*cef1->val(y);
      }
  
  * derivatives(double x, double y, double& dx, double& dy) const 

    ::
    
      void ExactSolutionFitzHughNagumo1::derivatives(double x, double y, 
        double& dx, double& dy) const 
      {
        dx = cef1->dx(x)*cef1->val(y);
        dy = cef1->val(x)*cef1->dx(y);
      }
    
  * ord(Ord x, Ord y) const 

    ::
    
      Ord ExactSolutionFitzHughNagumo1::ord(Ord x, Ord y) const 
      {
        return Ord(10);
      }
    
  * clone() const

    ::
    
      MeshFunction<double>* ExactSolutionFitzHughNagumo1::clone() const
      {
        return new ExactSolutionFitzHughNagumo1(this->mesh);
      }
      
      
Mathematical functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once again we shall use the example 06::

  hermes2d\test_examples\06-system-adapt\definitions.cpp (.h)
  
In the header file (definitions.h) we can see the following class definition::

  class CustomRightHandSide1: public Hermes2DFunction<double>
  {
  public:
    CustomRightHandSide1(double K, double d_u, double sigma);

    virtual double value(double x, double y) const;

    virtual Ord value(Ord x, Ord y) const;

    ~CustomRightHandSide1();

    CustomExactFunction1* cef1;
    CustomExactFunction2* cef2;
    double d_u, sigma;
  };
  
The important methods here are (definitions - method bodies from definitions.cpp):
  * value(double x, double y) const 

    ::
  
      double CustomRightHandSide1::value(double x, double y) const 
      {
        double Laplace_u = cef1->ddxx(x) * cef1->val(y)
                           + cef1->val(x) * cef1->ddxx(y);
        double u = cef1->val(x) * cef1->val(y);
        double v = cef2->val(x) * cef2->val(y);
        return -d_u * d_u * Laplace_u - u + sigma * v;
      }
  
  * value(Ord x, Ord y) const
  
    ::
      
      // Note here that we are saying that this function is ok to be integrated with a 
      // quadrature rule precise for polynomials of order 10.
      Ord CustomRightHandSide1::value(Ord x, Ord y) const 
      {
        return Ord(10);
      }
      
  * Note that there is no **clone** method here. That is because these classes - mathematical functions - are used in OpenMP paralell blocks only inside methods of already cloned class instances - like Form::value() etc.
  

Filters
~~~~~~~
There is a number of pre-defined Filters for you in::

  hermes2d\include\function\filter.h
  
These include
  * AngleFilter
  * VonMisesFilter
  * LinearFilter
  * ValFilter
  * MagFilter
  * ...
  
They all come from the base class template Filter<Scalar>.

Underneath there is a distinction between the filters that come from the classes SimpleFilter or DXDYFilter (real function of real solutions, or complex one of complex) and those coming from ComplexFilter (real function of complex solutions).

The difference is obvious, the Solution<Scalar> template instances it operates with differ: for SimpleFilter / DXDYFilter successors, the type (real vs. complex) depends on the type of the filter, and in the case of ComplexFilter successors it is always::

  MeshFunction<std::complex<double> >* solution
  
Also note that whereas SimpleFilter / DXDYFilter are class **templates** - as explained in the previous paragraph, the ComplexFilter is just a class, and it inherits from Filter<double>.

SimpleFilter serves for functions of the solutions(s) values, DXDYFilter for functions of the solution(s) derivatives.

The common method all filters must override is::

  virtual MeshFunction<Scalar>* clone() const

Then there is always the method **filter_fn(...)** that comes in the following versions::

  // SimpleFilter - values here represent the solution values, n is the number of points.
  virtual void filter_fn(int n, Hermes::vector<Scalar*> values, Scalar* result) = 0;
  
  // DXDYFilter - contains values, dx - derivatives w.r.t. x, dy - derivatives w.r.t. y,
  // and also the resulting derivatives, should those be necessary.
  virtual void filter_fn (int n, Hermes::vector<Scalar *> values, Hermes::vector<Scalar *> dx, Hermes::vector<Scalar *> dy, Scalar* rslt, Scalar* rslt_dx, Scalar* rslt_dy) = 0;
  
  // ComplexFilter - values here represent the solution values, n is the number of points,
  // note that here, the values are complex.
  virtual void filter_fn(int n, std::complex<double>* values, double* result) = 0;
  
  
  