General Initial Condition (16)
------------------------------

**Git reference:** Tutorial example `16-newton-elliptic-2 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/16-newton-elliptic-2>`_.

We will solve the nonlinear model problem from the previous section again,

.. math::

    -\nabla \cdot (\lambda(u)\nabla u) - f(x,y) = 0 \ \ \ \mbox{in } \Omega = (-10,10)^2

but now with nonhomogeneous Dirichlet boundary conditions 

.. math::

    u(x, y) = (x+10)(y+10)/100 \ \ \ \mbox{on } \partial \Omega

and with a general initial guess init_guess(x,y).

The treatment of the Dirichlet boundary conditions in the code looks as follows::

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

.. image:: 16/proj-h1.png
   :align: center
   :width: 600
   :height: 350
   :alt: H1 projection

The converged solution looks as follows:

.. image:: 16/solution.png
   :align: center
   :width: 600
   :height: 350
   :alt: approximate solution
