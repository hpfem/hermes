Picard's Method (15)
--------------------

**Git reference:** Tutorial example `15-picard 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/15-picard>`_.

Picard's Method
~~~~~~~~~~~~~~~

The Picard's method is a simple approach to the solution of nonlinear problems
where nonlinear products are linearized by moving the source of nonlinearity 
to the previous iteration level. For example, a nonlinear product of the form 
$g(u)u$ would be linearized as $g(u^n) u^{n+1}$. Let us illustrate this on a 
simple model problem.

Model problem
~~~~~~~~~~~~~

Let us solve a nonlinear equation

.. math::

    -\nabla \cdot (\lambda(u)\nabla u) = f(x,y), \ \ \ u = 0 \ \mbox{on}\ \partial \Omega.

One possible interpretation of this equation is stationary heat transfer where the thermal
conductivity $\lambda$ depends on the temperature $u$, and $f(x,y)$ are heat sources/losses.
Our domain is a square $\Omega = (-10,10)^2$, $f(x,y) = 1$, and the nonlinearity $\lambda$ has the form 

.. math::

    \lambda(u) = 1 + u^\alpha.

Recall that $\lambda$ must be entirely positive or entirely negative for the problem to be solvable, so it is safe 
to restrict $\alpha$ to be an even nonnegative integer. The linearized equation has the form 

.. math::

    -\nabla \cdot (\lambda(u^n)\nabla u^{n+1}) = f(x,y), \ \ \ u = 0 \ \mbox{on}\ \partial \Omega.

The Picard's iteration begins from some initial guess $u^0$, in our case a constant 
function, and runs until a convergence criterion is satisfied. Most widely used 
convergence criteria are the relative error between two consecutive iterations, or 
residual of the equation. In this example we will use the former.

Defining initial condition
~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses a constant initial guess::

    // Initialize previous iteration solution for the Picard method.
    sln_prev = new Solution();
    sln_prev->set_const(&mesh, INIT_COND_CONST);


Weak forms
~~~~~~~~~~

The weak forms are as the user expects::

    // Stiffness matrix.
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      Func<Scalar>* u_prev = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (lam(u_prev->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));

      return result;
    }

    // Right-hand side.
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      Func<Scalar>* u_prev = u_ext[0];
      for (int i = 0; i < n; i++) result += wt[i] * rhs(e->x[i], e->y[i]) * v->val[i];
      return result;
    }

Notice that the solution $u$ is accessed through

::

    Func<Scalar>* u_prev = ext->fn[0];

**Do not** use the variable u_ext[] which is not initialized for linear problems.

Registering weak forms
~~~~~~~~~~~~~~~~~~~~~~

The weak forms are registered with the previous iteration level solution 
being an external function::

    // Initialize the weak formulation.
    WeakForm wf;
    wf.add_matrix_form(callback(jac), HERMES_UNSYM, HERMES_ANY, sln_prev);
    wf.add_vector_form(callback(res), HERMES_ANY);

Initializing a linear DiscreteProblem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Picard's method deals with a linear problem, therefore the DiscreteProblem 
class is initialized with is_linear = true::

    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem dp(&wf, &space, is_linear);

Picard's iteration loop
~~~~~~~~~~~~~~~~~~~~~~~

The iteration loop is straightforward::

    // Picard's iteration.
    double rel_err;
    int iter_count = 0;
    while (true) {
      // Assemble the stiffness matrix and right-hand side.
      info("Assembling the stiffness matrix and right-hand side vector.");
      dp.assemble(matrix, rhs);

      // Solve the linear system and if successful, obtain the solution.
      info("Solving the matrix problem.");
      if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
      else error ("Matrix solver failed.\n");

      double rel_error = calc_abs_error(sln_prev, &sln, HERMES_H1_NORM) / calc_norm(&sln, HERMES_H1_NORM) * 100;
      info("Relative error: %g%%", rel_error);

      // Stopping criterion.
      if (rel_error < PICARD_TOL || iter_count >= MAX_PICARD_ITER_NUM) break;

      // Saving solution for the next iteration;
      sln_prev->copy(&sln);
   
      iter_count++;
    }

As a last step, we clean up as usual::

    // Cleanup.
    delete [] coeff_vec;
    delete matrix;
    delete rhs;
    delete solver;

Sample results
~~~~~~~~~~~~~~

Approximate solution $u$ for $\alpha = 2$: 

.. image:: 16/newton-ellipt-1-2.png
   :align: center
   :width: 600
   :height: 400
   :alt: result for alpha = 2

Approximate solution $u$ for $\alpha = 4$: 

.. image:: 16/newton-ellipt-1-4.png
   :align: center
   :width: 600
   :height: 400
   :alt: result for alpha = 4
