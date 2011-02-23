Time-Harmonic Maxwell's Equations (05-hcurl-adapt)
--------------------------------------

**Git reference:** Tutorial example `05-hcurl-adapt <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P04-linear-adapt/05-hcurl-adapt>`_. 

Model problem
~~~~~~~~~~~~~

This example solves the time-harmonic Maxwell's equations in an L-shaped domain and it 
describes the diffraction of an electromagnetic wave from a re-entrant corner. It comes with an 
exact solution that contains a strong singularity.

Equation solved: Time-harmonic Maxwell's equations

.. math::
    :label: example-05-hcurl-adapt

    \frac{1}{\mu_r} \nabla \times \nabla \times E - \kappa^2 \epsilon_r E = \Phi.

Domain of interest is the square $(-10, 10)^2$ missing the quarter lying in the 
fourth quadrant. It is filled with air:

.. image:: maxwell-adapt/domain.png
   :align: center
   :width: 490
   :height: 490
   :alt: Computational domain.

Boundary conditions: Combined essential and natural, see the 
`main.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/tutorial/P04-linear-adapt/05-hcurl-adapt/main.cpp>`_ file.

Exact solution
~~~~~~~~~~~~~~

.. math::
    :label: example-05-hcurl-adapt-exact

    E(x, y) = \nabla \times J_{\alpha} (r) \cos(\alpha \theta)

where $J_{\alpha}$ is the Bessel function of the first kind, 
$(r, \theta)$ the polar coordinates and $\alpha = 2/3$. In 
computer code, this reads:

::

    void exact_sol(double x, double y, scalar& e0, scalar& e1)
    {
      double t1 = x*x;
      double t2 = y*y;
      double t4 = sqrt(t1+t2);
      double t5 = jv(-1.0/3.0,t4);
      double t6 = 1/t4;
      double t7 = jv(2.0/3.0,t4);
      double t11 = (t5-2.0/3.0*t6*t7)*t6;
      double t12 = atan2(y,x);
      if (t12 < 0) t12 += 2.0*M_PI;
      double t13 = 2.0/3.0*t12;
      double t14 = cos(t13);
      double t17 = sin(t13);
      double t18 = t7*t17;
      double t20 = 1/t1;
      double t23 = 1/(1.0+t2*t20);
      e0 = t11*y*t14-2.0/3.0*t18/x*t23;
      e1 = -t11*x*t14-2.0/3.0*t18*y*t20*t23;
    }  

Here jv() is the Bessel function $\bfJ_{\alpha}$. For its source code see the 
`forms.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/tutorial/P04-linear-adapt/05-hcurl-adapt/forms.cpp>`_ file.

Creating an Hcurl space
~~~~~~~~~~~~~~~~~~~~~~~

New in this example is the fact that we solve in the Hcurl space::

    // Create an Hcurl space with default shapeset.
    HcurlSpace space(&mesh, bc_types, essential_bc_values, P_INIT);

Choosing refinement selector for the Hcurl space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Therefore we need to use a refinement selector for the Hcurl space::

    // Initialize refinement selector.
    HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

Projecting in the Hcurl norm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The H2D_HCURL_NORM needs to be used in the global projection::

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver, HERMES_HCURL_NORM); 

Initializing the Adapt class with the Hcurl norm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The H2D_HCURL_NORM is also used in the initialization of the Adapt class::

    // Calculate element errors and total error estimate.
    info("Calculating error estimate and exact error."); 
    Adapt* adaptivity = new Adapt(&space, HERMES_HCURL_NORM);
    bool solutions_for_adapt = true;
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

Exact error calculation and the 'solutions_for_adapt' flag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Above, solutions_for_adapt=true means that 'sln' and 'ref_sln' will be used to calculate element 
errors to guide adaptivity. In the following call to calc_err_exact() we set solutions_for_adapt=false 
so that just the total error is calculated::

    // Calculate exact error,
    solutions_for_adapt = false;
    double err_exact_rel = adaptivity->calc_err_exact(&sln, &sln_exact, solutions_for_adapt, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

Weak forms
~~~~~~~~~~

::

    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
    return 1.0/mu_r * int_curl_e_curl_f<Real, Scalar>(n, wt, u, v) -
           sqr(kappa) * int_e_f<Real, Scalar>(n, wt, u, v);
    }
   
    template<typename Real, typename Scalar>
    Scalar bilinear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      cplx ii = cplx(0.0, 1.0);
      return ii * (-kappa) * int_e_tau_f_tau<Real, Scalar>(n, wt, u, v, e);
    }
   
    scalar linear_form_surf(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      scalar result = 0;
      for (int i = 0; i < n; i++)
      {
        double r = sqrt(e->x[i] * e->x[i] + e->y[i] * e->y[i]);
        double theta = atan2(e->y[i], e->x[i]);
        if (theta < 0) theta += 2.0*M_PI;
        double j13    = jv(-1.0/3.0, r),    j23    = jv(+2.0/3.0, r);
        double cost   = cos(theta),         sint   = sin(theta);
        double cos23t = cos(2.0/3.0*theta), sin23t = sin(2.0/3.0*theta);
   
        double Etau = e->tx[i] * (cos23t*sint*j13 - 2.0/(3.0*r)*j23*(cos23t*sint + sin23t*cost)) +
                      e->ty[i] * (-cos23t*cost*j13 + 2.0/(3.0*r)*j23*(cos23t*cost - sin23t*sint));
  
        result += wt[i] * cplx(cos23t*j23, -Etau) * ((v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
      }
      return result;
    }

    // Maximal polynomial order to integrate surface linear form.
    Ord linear_form_surf_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {  return Ord(v->val[0].get_max_order());  }


Sample results
~~~~~~~~~~~~~~

Solution:

.. image:: maxwell-adapt/solution.png
   :align: center
   :width: 500
   :height: 420
   :alt: Solution.

Final mesh (h-FEM with linear elements):

.. image:: maxwell-adapt/mesh-h1.png
   :align: center
   :width: 460
   :height: 390
   :alt: Final mesh (h-FEM with linear elements).

Note that the polynomial order indicated corresponds to the tangential components 
of approximation on element interfaces, not to polynomial degrees inside the elements
(those are one higher).

Final mesh (h-FEM with quadratic elements):

.. image:: maxwell-adapt/mesh-h2.png
   :align: center
   :width: 460
   :height: 390
   :alt: Final mesh (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: maxwell-adapt/mesh-hp.png
   :align: center
   :width: 460
   :height: 390
   :alt: Final mesh (hp-FEM).

DOF convergence graphs:

.. image:: maxwell-adapt/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: maxwell-adapt/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

