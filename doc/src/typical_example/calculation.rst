Calculation
~~~~~~~~~~~~
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