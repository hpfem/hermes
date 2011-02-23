Constant Initial Condition (02-newton-1)
-------------------------------

**Git reference:** Tutorial example `02-newton-1 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P02-nonlinear/02-newton-1>`_.

Model problem
~~~~~~~~~~~~~

Let us solve the nonlinear model problem from the previous section,

.. math::

    -\nabla \cdot (\lambda(u)\nabla u) - f(x,y) = 0, \ \ \ u = 0 \ \mbox{on}\ \partial \Omega.

Recall that the domain is a square $\Omega = (-10,10)^2$, $f(x,y) = 1$, and the nonlinearity $\lambda$ 
has the form 

.. math::

    \lambda(u) = 1 + u^\alpha

where $\alpha$ is an even nonnegative integer. Also recall from the previous section that 

.. math::

    F_i(\bfY) =  \int_{\Omega} \lambda(u)\nabla u \cdot \nabla v_i 
    - f v_i \, \mbox{d}x\mbox{d}y.

and

.. math::

    J_{ij}(\bfY) =
    \int_{\Omega} \left[ \frac{\partial \lambda}{\partial u}(u) v_j 
    \nabla u + \lambda(u)\nabla v_j \right] \cdot \nabla v_i \, \mbox{d}x\mbox{d}y.

Defining Jacobian matrix and residual vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the `code <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/tutorial/16-newton-1/forms.cpp>`_, 
the above formulas become::

    // Heat sources (can be a general function of 'x' and 'y').
    template<typename Real>
    Real heat_src(Real x, Real y)
    {
      return 1.0;
    }

    // Jacobian matrix
    template<typename Real, typename Scalar>
    Scalar jac(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      Func<Scalar>* u_prev = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (dlam_du(u_prev->val[i]) * u->val[i] * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
                           + lam(u_prev->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));

      return result;
    }

    // Fesidual vector
    template<typename Real, typename Scalar>
    Scalar res(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      Func<Scalar>* u_prev = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (lam(u_prev->val[i]) * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
	    	           - heat_src(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

Notice that the solution $u$ is accessed through

::

    Func<Scalar>* u_prev = u_ext[0];

In principle, one could register a previous iteration level solution as in the Picard
iteration, but this is simpler and more elegant. 

Also notice that the basis function $v_j$ and the test function 
$v_i$ are entering the weak forms via the parameters u and v, respectively (same as for linear 
problems). The user does not have to 
take care about their indices $i$ and $j$, this is handled by Hermes outside the weak forms. 

Using the ExtData and Geom structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The code snippet above also shows how values and derivatives of the solution $u$ can be accessed via 
the ExtData structure, and the coordinates of the integration points via the Geom structure. 
The contents of ExtData is user-defined and the Geom structure contains geometrical information 
including the unit normal and tangential vectors to the boundary at the integration points 
(also for curved boundaries). See the file 
`src/forms.h <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/src/forms.h>`_ for more details. 

Defining nonlinearities and their derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function $\lambda(u)$ and its derivative are defined as follows::

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u) 
    { 
      return 1 + pow(u, 4); 
    }

    // Derivative of the thermal conductivity with respect to 'u'.
    template<typename Real>
    Real dlam_du(Real u) 
    { 
      return 4*pow(u, 3); 
    }

Setting a constant initial condition for the Newton's method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Newton's method always starts from an initial coefficient vector $\bfY_0$.
If we want to start from a zero initial function, we can just define this vector 
to be zero. However, more often we want to start from some nonzero function
(such as, for example, the previous time-level solution in time-dependent problems). 
In such a case, the initial coefficient vector is obtained by projecting the 
initial solution on the finite element space. In the present example, the initial 
solution is a constant function::

    // Project the initial condition on the FE space to obtain initial
    // coefficient vector for the Newton's method.
    info("Projecting to obtain initial vector for the Newton's method.");
    scalar* coeff_vec = new scalar[Space::get_num_dofs(&space)];
    Solution* init_sln = new Solution();
    init_sln->set_const(&mesh, INIT_COND_CONST);
    OGProjection::project_global(&space, init_sln, coeff_vec, matrix_solver);
    delete init_sln;

The method project_global() has an optional parameter which is the projection 
norm. Its default value is HERMES_H1_NORM but other norms such as HERMES_HCURL_NORM,
HERMES_HDIV_NORM, and HERMES_L2_NORM are also possible. This will be explained 
later and we'll also show how to handle projections for systems of equations.

The user is at liberty to use for the (always symmetric positive definite) 
projection matrix a different matrix solver
than for the solution of the matrix problems arising in the Newton's iteration. 

Registering weak forms
~~~~~~~~~~~~~~~~~~~~~~

The weak forms are registered as usual::

    // Initialize the weak formulation.
    WeakForm wf;
    wf.add_matrix_form(callback(jac), HERMES_NONSYM, HERMES_ANY);
    wf.add_vector_form(callback(res), HERMES_ANY);

Recall that by HERMES_NONSYM we declare that the Jacobian bilinear form is not symmetric,
and by HERMES_ANY that the form should be used for elements with any material marker.

Initializing a nonlinear DiscreteProblem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As opposed to linear problems, the DiscreteProblem class is now initialized with 
the boolean flag is_linear=false::

    // Initialize the FE problem.
    bool is_linear = false;
    DiscreteProblem dp(&wf, &space, is_linear);

The Newton's iteration loop
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Newton's iteration loop is very similar in all examples, therefore we 
provide a simple function solve_newton() that is called as follows::

    // Perform Newton's iteration.
    bool verbose = true;
    if (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
        NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

The same written in full would be::

    // Perform Newton's iteration.
    bool verbose = true;
    int it = 1;
    while (1)
    {
      // Obtain the number of degrees of freedom.
      int ndof = Space::get_num_dofs(&space);

      // Assemble the Jacobian matrix and residual vector.
      dp.assemble(coeff_vec, matrix, rhs, false);

      // Multiply the residual vector with -1 since the matrix 
      // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
      for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));
    
      // Calculate the l2-norm of residual vector.
      double res_l2_norm = get_l2_norm(rhs);

      // Info for the user.
      if (verbose) info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, 
                   Space::get_num_dofs(&space), res_l2_norm);

      // If l2 norm of the residual vector is within tolerance, or the maximum number 
      // of iteration has been reached, then quit.
      if (res_l2_norm < NEWTON_TOL || it > NEWTON_MAX_ITER) break;

      // Solve the linear system.
      if(!solver->solve()) error ("Matrix solver failed.\n");

      // Add \deltaY^{n+1} to Y^n.
      for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];
    
      if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");

      it++;
    }

Note that the Newton's loop always handles a coefficient vector, not 
solutions. 

Translating the resulting vector into a Solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the Newton's loop is finished, the resulting coefficient vector 
is translated into a Solution as follows::

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solution(coeff_vec, &space, &sln);

As a last step, we clean up as usual::

    // Cleanup.
    delete [] coeff_vec;
    delete matrix;
    delete rhs;
    delete solver;

Sample results
~~~~~~~~~~~~~~

The results are exactly the same as in the Picard's example 01-picard. 
Notice that the Newton's method uses very few iterations compared
to Picard.
