Weak formulation
~~~~~~~~~~~~~~~~
When we already have a mesh and a space, we have to know what equations we will be solving on those. And that is where the weak formulation comes to the light.
Of course, there is a vast mathematical background of differential equations, their weak solutions, Sobolev spaces, etc., but we assume of those, our users already have a good knowledge. Right here we are concerned with the implementation. A typical creation of a weak formulation for the use with Hermes might look like this::
 
    // Initialize the weak formulation.
    // This is a weak formulation for linear elasticity, with custom
    // parameters of the constructor.
    // There is a lot of documentation for using some predefined
    // weak forms, as well as creating your own. See the info below.
    CustomWeakFormLinearElasticity wf(E, nu, rho*g1, "Top", f0, f1);
    
More about a typical basic weak form can be found in the 'hermes-tutorial' documentation, section 'A-linear', chapter '03-poisson'.

More about creating a custom weak form can be found in other tutorial examples. One always needs to subclass the Hermes::Hermes2D::WeakForm<Scalar> class template.
For defining custom forms (integrals), one needs to subclass templates Hermes::Hermes2D::MatrixFormVol<Scalar>, MatrixFormSurf<Scalar>, VectorFormVol<Scalar>, VectorFormSurf<Scalar>.
A typical constructor of a derived class::

    template<> DefaultMatrixFormVol<std::complex<double> >::DefaultMatrixFormVol
      (int i, int j, std::string area, Hermes2DFunction<std::complex<double> >* coeff,
      SymFlag sym, GeomType gt)
      : MatrixFormVol<std::complex<double> >(i, j, area, sym), coeff(coeff), gt(gt)
    {
      // If coeff is HERMES_ONE, initialize it to be constant 1.0.
      if(coeff == HERMES_ONE)
        this->coeff = new Hermes2DFunction<std::complex<double> >
        (std::complex<double>(1.0, 1.0));
    }
    
In this constructor::

    template<> DefaultMatrixFormVol<std::complex<double> >
    // means that this is an explicit instantiation of a template for
    // complex numbers (DefaultMatrixFormVol, the derived class, is actually also
    // a template, as the prent class is).
    
    int i, j
    // coordinates in the system of equations, first is the row (basis functions),
    // second the column (test functions).
    
    std::string area //(typically optional)
    // either a std::string for the marker on which this form will be evaluated,
    // or HERMES_ANY constant for 'any', i.e. all markers 
    // (this is the default in the parent class constructor).
    
    Hermes2DFunction<std::complex<double> >*coeff //(typically optional)
    // custom function having its meaning specified in the 
    // calculating methods (see further). The constant HERMES_ONE,
    // that really represents the number 1.0, is the default 
    // in the parent class constructor.
   
    SymFlag sym //(typically optional)
    // symmetry flag 
    // see the 'hermes-tutorial' documentation, section 'A-linear', chapter '03-poisson'.
    
    GeomType gt //(typically optional)
    // type of geometry: HERMES_PLANAR, HERMES_AXISYM_X, HERMES_AXISYM_Z,
    // to distinguish between the normal 2D settings (HERMES_PLANAR),
    // or an axisymmetric one. See the 'hermes-tutorial' documentation, 
    // section 'A-linear', chapter '09-axisym' for more details.
    
In those, the main methods to override are value(...), and ord(...), calculating the value and integration order respectively. It is a good idea to refer to the default forms (located in the library repository, with headers in hermes2d/include/weakform_library/*.h and the sources in hermes2d/src/weakform_library/*.cpp).
The header is pretty self-explanatory::

    // MatrixForm.
    virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const;
        
    // A typical implementation.
    template<typename Scalar> Scalar DefaultMatrixFormVol<Scalar>::value(int n, 
      double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, Func<Scalar> **ext) const
    {
      Scalar result = 0;
      
      for (int i = 0; i < n; i++)
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
    }
        
    // VectorForm.
    // Identical to MatrixForm, only the basis function is missing for obvious reasons.
    virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const;
        
In these::
    
    int n
    // number of integration points.
    
    double *wt
    // integration weights (an array containing 'n' values).
    
    Func<Scalar> *u_ext
    // values from previous Newton iterations, as many as there are spaces 
    // (equations) in the system.
    
    Func<double> *u
    // the basis function, represented by the class Func. 
    // For more info about the class, see the developers documentation (in doxygen).
    // How to get that, see the documentation section.
    
    Func<double> *v
    // the test function, represented by the class Func.
    // For more info about the class, see the developers documentation (in doxygen).
    
    Geom<double> *e
    // geometry attributes: coordinates, element size,
    // normal directions (for surface forms), you name it.
    // For more info about the class, see the developers documentation (in doxygen).
    
    Func<Scalar> **ext
    // external functions, as many as you like 
    // (provided you set it up in constructor of your weak formulation 
    // derived from the class WeakForm). 
    // For more info about the class, see the developers documentation (in doxygen).
    
Now we have a space and a weak formulation, we are ready to calculate!