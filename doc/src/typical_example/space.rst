Space
~~~~~~~~
Secondly, the Finite Element space must be set on the computational mesh. One of the following is typically used (including setting of Dirichlet boundary conditions)::
    
    // H1 Space.
    // Polynomial order.
    int POLYNOMIAL_ORDER = 3;
    
    // Initialize boundary conditions.
    // This is a custom (derived) boundary condition. More about this in the section 'Object model - deriving your own specialized classes'.
    CustomDirichletCondition bc_essential(
      Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"),
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
