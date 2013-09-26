PARALUTION
-------

.. _PARALUTION home page: http://www.paralution.com

1) Download PARALUTION from http://www.paralution.com/
2) Compile, build PARALUTION, and copy headers, and libraries so that Hermes's CMake system can find them (as with other dependencies)
2.1) On Linux, installing PARALUTION to default install directories is sufficient, on Windows some paths have to be set
3) In your CMake.vars file (or directly in CMakeLists.txt) in the root of Hermes (see step 0) add "set(WITH_PARALUTION YES)"
3.1) It is on by default, so by default one has to include PARALUTION to build Hermes
4) That is it, build Hermes, it will automatically link to PARALUTION, include headers and make it usable.

How to use PARALUTION
--------------
5) Read the doxygen manual of the classes
5.1) Hermes::Algebra::ParalutionMatrix
5.2) Hermes::Algebra::ParalutionVector
5.3) Hermes::Preconditioners::ParalutionPrecond
5.4) Hermes::Solvers::IterativeParalutionLinearMatrixSolver
5.5) Hermes::Solvers::AMGParalutionLinearMatrixSolver
5.6) and all classes that these inherit from
6) If you want to see Hermes & PARALUTION readily work together, take any test example in the /hermes2d folder in the Hermes root and add one of these lines at the beginning of your main()
6.1) HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_ITERATIVE); // to use iterative solver
6.2) HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_AMG); // to use AMG solver
7) Solver classes of Hermes (NewtonSolver, PicardSolver, LinearSolver, ...) will then take this API setting into account and use PARALUTION as the matrix solver.