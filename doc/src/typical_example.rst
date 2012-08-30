================================
Hermes typical example structure
================================
- a nice thing about Hermes is that it follows the math tightly, so the steps taken in solving the example with Hermes correspond to those taken in theory.

A beginning of each example can look like this::
    
    include "hermes2d.h"
    // Two basic namespaces.
    using namespace Hermes;
    using namespace Hermes::Hermes2D;
    // For adaptivity.
    using namespace Hermes::Hermes2D::RefinementSelectors;
    // For visualization.
    using namespace Hermes::Hermes2D::Views;

Mesh
-----------------------------------------
First part one needs to handle is the computational mesh, typically the following would be used::

    // Native Hermes mesh format.
    MeshReaderH2D mloader;
    mloader.load("domain.mesh", &mesh);
    
    // XML Hermes mesh format.
    MeshReaderH2DXML mloader;  
    mloader.load("domain.xml", &mesh);
    
More about meshes can be found in the 'hermes-tutorial' documentation, section 'A-linear', chapter '01-mesh'.

Space
-----
Secondly, the Finite Element space must be set on the computational mesh. One of the following is typically used (including setting of Dirichlet boundary conditions)::
    
    // H1 Space.
    // Polynomial order.
    int POLYNOMIAL_ORDER = 3;
    // Initialize boundary conditions.
    CustomDirichletCondition bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"),
      BDY_A_PARAM, BDY_B_PARAM, BDY_C_PARAM);
    EssentialBCs<double> bcs(&bc_essential);
    // Create an H1 space.
    H1Space<double> space(&mesh, &bcs, POLYNOMIAL_ORDER);
    
    // HCurl Space.
    // Polynomial order.
    int POLYNOMIAL_ORDER = 5;
    // Initialize boundary conditions.
    Hermes::Hermes2D::DefaultEssentialBCConst<std::complex<double> > bc_essential(Hermes::vector<std::string>("Corner_horizontal",
      "Corner_vertical"), 0);
    EssentialBCs<std::complex<double> > bcs(&bc_essential);
    // Create an Hcurl space.
    HcurlSpace<std::complex<double> > space(&mesh, &bcs, POLYNOMIAL_ORDER);
    
    // HDiv Space. This example does not use any Dirichlet boundary conditions.
    int POLYNOMIAL_ORDER = 2;
    HdivSpace<double> space(&mesh, POLYNOMIAL_ORDER);
    
    // L2 Space. This Space does not take any boundary conditions which corresponds to the fact that the FE space is a space of discontinuous functions.
    // If we for example use polynomial order = 0, we use just piecewise constant basis functions.
    L2Space<double> space(&mesh, 0);
    
More about spaces can be found in the 'hermes-tutorial' documentation, section 'A-linear', chapter '02-space'.

More about Dirichlet boundary conditions can be found in the 'hermes-tutorial' documentation, section 'A-linear', chapter '04-bc-dirichlet', and for defining a non-constant custom boundary condition, see the chapter '07-general'.

Weak formulation
----------------
When we already have a mesh and a space, we have to know what equations we will be solving on those. And that is where the weak formulation comes to the light.
Of course, there is a vast mathematical background of differential equations, their weak solutions, Sobolev spaces, etc., but we assume of those, our users already have a good knowledge. Right here we are concerned with the implementation. A typical creation of a weak formulation for the use with Hermes might look like this::
 
    // Initialize the weak formulation.
    // This is a weak formulation for linear elasticity, with custom parameters of the constructor.
    // There is a lot of documentation for using some predefined weak forms, as well as creating your own. See the info below.
    CustomWeakFormLinearElasticity wf(E, nu, rho*g1, "Top", f0, f1);
    
More about a typical basic weak form can be found in the 'hermes-tutorial' documentation, section 'A-linear', chapter '03-poisson'.

More about creating a custom weak form can be found in other tutorial examples. One always needs to subclass the Hermes::Hermes2D::WeakForm<Scalar> class template.
For defining custom forms (integrals), one needs to subclass templates Hermes::Hermes2D::MatrixFormVol<Scalar>, MatrixFormSurf<Scalar>, VectorFormVol<Scalar>, VectorFormSurf<Scalar>.
A typical constructor of a derived class::

    template<> DefaultMatrixFormVol<std::complex<double> >::DefaultMatrixFormVol
      (int i, int j, std::string area, Hermes2DFunction<std::complex<double> >* coeff, SymFlag sym, GeomType gt)
      : MatrixFormVol<std::complex<double> >(i, j, area, sym), coeff(coeff), gt(gt)
    {
      // If coeff is HERMES_ONE, initialize it to be constant 1.0.
      if(coeff == HERMES_ONE)
        this->coeff = new Hermes2DFunction<std::complex<double> >(std::complex<double>(1.0, 1.0));
    }
    
In this constructor::

    template<> DefaultMatrixFormVol<std::complex<double> > - means that this is an explicit instantiation of a template for complex numbers (DefaultMatrixFormVol, the derived class, is actually also a template, as the prent class is).
    i, j - coordinates in the system of equations, first is the row (basis functions), second the column (test functions).
    area (typically optional) - either a std::string for the marker on which this form will be evaluated, of HERMES_ANY constant for 'any', i.e. all markers (this is the default in the parent class constructor).
    coeff (typically optional) - custom function having its meaning specified in the calculating methods (see further). The constant HERMES_ONE, that really represents the number 1.0, is the default in the parent class constructor.
    sym (typically optional) - symmetry flag - see the 'hermes-tutorial' documentation, section 'A-linear', chapter '03-poisson'.
    gt - type of geometry: HERMES_PLANAR, HERMES_AXISYM_X, HERMES_AXISYM_Z, to distinguish between the normal 2D settings (HERMES_PLANAR), or an axisymmetric one. See the 'hermes-tutorial' documentation, section 'A-linear', chapter '09-axisym' for more details.
    
In those, the main methods to override are value(...), and ord(...), calculating the value and integration order respectively. It is a good idea to refer to the default forms (located in the library repository, with headers in hermes2d/include/weakform_library/*.h and the sources in hermes2d/src/weakform_library/*.cpp).
The header is pretty self-explanatory::

    // MatrixForm.
    virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const;
        
    // A typical implementation.
    template<typename Scalar> Scalar DefaultMatrixFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, ExtData<Scalar> *ext) const
    {
      Scalar result = 0;
      
      for (int i = 0; i < n; i++)
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
    }
        
    // VectorForm.
    // Identical to MatrixForm, only the basis function is missing for obvious reasons.
    virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const;
        
In these::
    
    n - number of integration points.
    wt - integration weights (an array containing 'n' values).
    u_ext - values from previous Newton iterations, as many as there are spaces (equations) in the system.
    u - the basis function, represented by the class Func. For more info about the class, see the developers documentation (in doxygen). How to get that, see the previous page.
    v - the test function, represented by the class Func. For more info about the class, see the developers documentation (in doxygen). How to get that, see the previous page.
    e - geometry attributes: coordinates, element size, normal directions (for surface forms), you name it. For more info about the class, see the developers documentation (in doxygen). How to get that, see the previous page.
    ext - external functions, as many as you like (provided you set it up in constructor of your weak formulation derived from the class WeakForm. For more info about the class, see the developers documentation (in doxygen). How to get that, see the previous page.
    
Now we have a space and a weak formulation, we are ready to calculate!

Calculation
-----------
We are going to get a solution vector from what we already have in the most general setup. This means for a time-dependent, adaptive example. 
This is to illustrate the various classes and methods, and the best thing about them, they are used pretty much the same way. 

For details about time-dependent examples, and various aspects of that, see the 'hermes-tutorial' documentation, section 'C-transient'.
For details about adaptive examples, and various aspects of that, see the 'hermes-tutorial' documentation, section 'D-adaptivity'.
Right here we focus on the calculation::

    double current_time = 'something';
    double current_time_step = 'also something';
    Time-loop
    {
      Adaptive-loop // not necessarily on each time step.
      {
        // create reference space(s), see the adaptivity section of hermes-tutorial documentation for this.
        // e.g.
        Space<double>* ref_space = construct_refined_space(&space);
      
        // WE ARE NOW HERE.
        The calculation
        // WE ARE NOW HERE
        
        // do the adaptivity thing, see the adaptivity section of hermes-tutorial documentation for this.
        // This would change the 'coarse' Space instance: 'Space<double> space'.
        
        // Do some cleaning.
      }
      // adjust time, and time step any way you want (stability conditions, time-adaptivity, ...).
    }

In the following, some parameters are supposed to be passed as 'const', and everywhere, where one can pass an instance of Space<Scalar>*, one can pass 
an instance of Hermes::vector<Space<Scalar>*> (Hermes::vector<const Space<Scalar>*>), if the problem is a system of equations.

We shall start from the simplest case. 

1 - linear example
~~~~~~~~~~~~~~~~~~

Once we created the Space(s), and the WeakForm (always one!), we create (outside of the loops!) an instance of LinearSolver<Scalar>::

    // Let us say that for real numbers, but for complex, it would be analogic.
    LinearSolver<double> solver(&space, &wf).

Then the following code inside the loops would do the trick::

    solver.setTime(current_time);
    solver.setTimeStep(current_time_step);
    // Yes! You are right, these can be used outside of the adaptivity loop!
    
    // Set the new Space.
    solver.set_space(ref_space);
    
    // This is usually the place where something can go wrong, so we use the try-catch block. Note that the exceptions we use in Hermes are std::exception descendants (so only one catch block is enough).
    try
    {
      // Do the magic (assemble the matrix, right-hand side, solve).
      solver.solve();
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }
    
    // Get the solution (The Solution class is described in the developers documentation, The method vector_to_solution(s) too)
    Solution<double> ref_sln;
    Solution<double>::vector_to_solution(solver.get_sln_vector(), ref_sln, &ref_sln);
    
And that is it, we have the solution of the problem in that adaptivity step on that time level. What to do with it (visualize, do some calculations, projections, limiting, whatever) and how to do it is described in various points in the tutorial.
But let us say that we would like to see it, the following will make us happy::

    ScalarView<double> view("My Solution");
    view.show(&ref_sln);

2 - nonlinear example
~~~~~~~~~~~~~~~~~~~~~

In this surprisingly short section, we will learn how to use NewtonSolver, and PicardSolver.

It is literaly the same as in the previous section, just take out 'LinearSolver', and pass 'NewtonSolver', or 'PicardSolver'.
    
There is one more thing, if you want your NewtonSolver not to start from a zero initial guess, the following helps::

    // Initialize the vector for initial guess. Real case.
    // Also do not forget to 'delete []' this after you do not need it.
    double* coeff_vec = new double[Space<double>::get_num_dofs(&ref_space)];
    
    // For example let us project the previous time level solution and use it as initial guess.
    OGProjection<double> ogProjection;
    ogProjection.project_global(ref_space, &previous_time_level_sln, coeff_vec);
    
    // And now use it in the NewtonSolver<Scalar>::solve (solver is now NewtonSolver<double>) method.
    solver.solve(coeff_vec);
    
One can also use the NOX solver from the Trilinos package (with analogic, but not exactly same methods). One needs Trilinos for that. And documentation for that is coming.

3 - RungeKutta solver.
~~~~~~~~~~~~~~~~~~~~~~

Again, pretty much the same as in the LinearSolver case, but the solve() method will now take the previous time level Solution(s) and return the new Solution(s), so there is no need for using the vector_to_solution(s) method::

    // Initialize the solution(it can be outside of the loops, the solution would always be rewritten when it is natural)
    Solution<double> ref_sln;
    
    // "solver" is now an instance of RungeKutta<double>.
    solver.setTime(current_time);
    solver.setTimeStep(current_time_step);
    // Yes! You are right, these can be used outside of the adaptivity loop!
    
    // Set the new Space.
    solver.set_space(ref_space);
    
    // This is usually the place where something can go wrong, so we use the try-catch block. Note that the exceptions we use in Hermes are std::exception descendants (so only one catch block is enough).
    try
    {
      // Do the usual magic, plus put the result in the ref_sln instance.
      solver.solve(&previous_time_level_sln, &ref_sln);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }

4 - use DiscreteProblem class directly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For special purposes, like DG or FVM (Finite Volume Method), where one needs to access the matrix or right-hand side, or needs to have the solution in hand before projection (to do limiting etc.), one can also directly use this class.

It shares some methods with the above 'calculation' classes, but of course does not do any calculation. The usage would look like this::

    // We assume we have an instance DiscreteProblem<double> dp(&wf, &space);
    
    // These can be outside the loop, the memory would get properly freed / reallocated every time without worrying about it.
    SparseMatrix<double>* matrix = create_matrix<double>();
    Vector<double>* rhs = create_vector<double>();
    LinearMatrixSolver<double>* linear_matrix_solver = create_linear_solver<double>(matrix, rhs);
    
    dp.setTime(current_time);
    dp.setTimeStep(current_time_step);
    
    // Set the new Space.
    dp.set_space(ref_space);
    
    // This is usually the place where something can go wrong, so we use the try-catch block. Note that the exceptions we use in Hermes are std::exception descendants (so only one catch block is enough).
    try
    {
      dp.assemble(matrix, rhs);
      
      // NOW WE HAVE THE MATRIX and RHS ASSEMBLED and we can do whatever we want with it.
      linear_matrix_solver.solve();
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }
    
    // Get the solution (The Solution class is described in the developers documentation, The method vector_to_solution(s) too)
    Solution<double> ref_sln;
    Solution<double>::vector_to_solution(linear_matrix_solver.get_sln_vector(), ref_sln, &ref_sln);
    
And that is it. There is not much more to it. See the 'transient', and 'adaptivity' sections of the hermes-tutorial documentation and all will fall into place.

Of course every problem is different, such as in the case of DG, one needs to do some limiting, shock capturing etc...
One can also save / load various entities (Spaces, Solutions, Meshes, time steps, ...) during calculation.

And especially, one needs to be careful not to forget deallocating stuff. How to do that, see the hermes-tutorial, and hermes-examples repositories. The examples there should be done properly.