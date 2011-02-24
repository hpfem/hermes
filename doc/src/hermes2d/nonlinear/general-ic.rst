General Initial Condition (03-newton-2)
------------------------------

**Git reference:** Tutorial example `03-newton-2 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P02-nonlinear/03-newton-2>`_.

Model problem
~~~~~~~~~~~~~

We will still stay with the nonlinear model problem from the previous sections,

.. math::

    -\nabla \cdot (\lambda(u)\nabla u) - f(x,y) = 0 \ \ \ \mbox{in } \Omega = (-10,10)^2

but now we will change the boundary conditions to nonhomogeneous Dirichlet,

.. math::

    u(x, y) = (x+10)(y+10)/100 \ \ \ \mbox{on } \partial \Omega

and with a general initial guess init_guess(x,y).

Defining nonhomogeneous Dirichlet boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    // This function is used to define Dirichlet boundary conditions.
    double dir_lift(double x, double y, double& dx, double& dy) {
      dx = (y+10)/10.;
      dy = (x+10)/10.;
      return (x+10)*(y+10)/100.;
    }

    // Boundary condition types.
    BCType bc_types(int marker)
    {
      return BC_ESSENTIAL;
    }

    // Essential (Dirichlet) boundary condition values.
    scalar essential_bc_values(int ess_bdy_marker, double x, double y)
    {
      double dx, dy;
      return dir_lift(x, y, dx, dy);
    }

Setting a nonconstant initial condition for the Newton's method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The initial condition has the form::

    // Initial condition. It will be projected on the FE mesh 
    // to obtain initial coefficient vector for the Newton's method.
    scalar init_cond(double x, double y, double& dx, double& dy)
    {
      // Using the Dirichlet lift elevated by two
      double val = dir_lift(x, y, dx, dy) + 2;
      return val;
    }

As in the previous example, the initial condition is projected on the finite element space 
to obtain an initial coefficient vector $\bfY_0$ for the Newton's iteration::

    // Project the initial condition on the FE space to obtain initial 
    // coefficient vector for the Newton's method.
    info("Projecting to obtain initial vector for the Newton's method.");
    scalar* coeff_vec = new scalar[Space::get_num_dofs(&space)] ;
    Solution* init_sln = new Solution(&mesh, init_cond);
    OGProjection::project_global(&space, init_sln, coeff_vec, matrix_solver); 
    delete init_sln;

Sample results
~~~~~~~~~~~~~~

The following figure shows the $H^1$-projection of the initial condition init_cond():

.. image:: general-ic/proj-h1.png
   :align: center
   :width: 600
   :height: 350
   :alt: H1 projection

The converged solution looks as follows:

.. image:: general-ic/solution.png
   :align: center
   :width: 600
   :height: 350
   :alt: approximate solution

