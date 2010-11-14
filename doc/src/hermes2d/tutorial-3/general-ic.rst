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

The treatment of the Dirichlet boundary conditions in the code looks as follows:

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

The initial condition has the form::

    // Initial condition. It will be projected on the FE mesh 
    // to obtain initial coefficient vector for the Newton's method.
    scalar init_cond(double x, double y, double& dx, double& dy)
    {
      // Using the Dirichlet lift elevated by two
      double val = dir_lift(x, y, dx, dy) + 2;
      return val;
    }

The initial condition must be projected on the finite element space 
in order to obtain the initial coefficient vector $\bfY_0$ for the Newton's
iteration::

    // Project the function init_cond() on the FE space
    // to obtain initial coefficient vector for the Newton's method.
    info("Projecting initial condition to obtain initial vector for the Newton'w method.");
    nls.project_global(init_cond, &u_prev);

Recall that the vector $\bfY_0$ can be retrieved from the NonLinSystem
class using the method get_solution_vector(). 

The following figure shows the $H^1$-projection of the initial condition init_cond():

.. image:: 16/proj-h1.png
   :align: center
   :width: 600
   :height: 350
   :alt: H1 projection

The Newton's iteration is again performed using

::

  // Perform Newton's iteration.
  info("Performing Newton's iteration.");
  bool verbose = true; // Default is false.
  if (!nls.solve_newton(&u_prev, NEWTON_TOL, NEWTON_MAX_ITER, verbose)) 
    error("Newton's method did not converge.");

The converged solution looks as follows:

.. image:: 16/solution.png
   :align: center
   :width: 600
   :height: 350
   :alt: approximate solution
