Nonlinear Parabolic Problem (22)
--------------------------------

**Git reference:** Tutorial example `22-newton-timedep-heat-adapt 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/22-newton-timedep-heat-adapt>`_.

Model problem
~~~~~~~~~~~~~

We consider a nonlinear time-dependent equation of the form 

.. math::
    \frac{\partial u}{\partial t} - \mbox{div}(\lambda(u)\nabla u) = f

that is equipped with Dirichlet boundary 
conditions and solved in a square $\Omega = (-10, 10)^2$. The nonlinearity $\lambda(u)$
has the form 

.. math::
    \lambda(u) = 1 + u^4.

Weak forms for the implicit Euler method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We discretize the time derivative using the implicit Euler method, and 
obtain the weak formulation

.. math::
    \int_{\Omega} \frac{u^{n+1}}{\tau}v \, \mbox{d}\bfx - 
    \int_{\Omega} \frac{u^{n}}{\tau}v \, \mbox{d}\bfx +
    \int_{\Omega} \lambda(u^{n+1})\nabla u^{n+1}\cdot \nabla v \, \mbox{d}\bfx -
    \int_{\Omega} f^{n+1}v \, \mbox{d}\bfx = 0.

In this example the heat sources $f$ do not depend on time, thus we can 
drop the index and just use $f \equiv f^{n+1}$. Expressing 

.. math::
    u(\bfY) = \sum_{k=1}^N y_k v_k

as usual in the Newton's method, where $\bfY = (y_1, y_2, \ldots, y_N)^T$,
we can write the residual function $\bfF(\bfY)$ in the form 

.. math::
    F_i(\bfY) = \int_{\Omega} \frac{u(\bfY)}{\tau}v_i \, \mbox{d}\bfx - 
    \int_{\Omega} \frac{u^{n}}{\tau}v_i \, \mbox{d}\bfx +
    \int_{\Omega} \lambda(u(\bfY))\nabla u(\bfY)\cdot \nabla v_i \, \mbox{d}\bfx -
    \int_{\Omega} fv_i \, \mbox{d}\bfx = 0.

The Jacobian matrix $\bfJ(\bfY)$ has the form 

.. math::
    J_{ij}(\bfY) = \frac{\partial F_i}{\partial y_j}(\bfY) = \int_{\Omega} \frac{v_j}{\tau}v_i \, \mbox{d}\bfx +
    \int_{\Omega} \frac{\partial \lambda}{\partial u}(u(\bfY))v_j \nabla u(\bfY)\cdot \nabla v_i \, 
    \mbox{d}\bfx +
    \int_{\Omega} \lambda(u(\bfY))\nabla v_j\cdot \nabla v_i \, \mbox{d}\bfx.

In the code, this looks as follows::

    // Residual for the implicit Euler time discretization
    template<typename Real, typename Scalar>
    Scalar F_euler(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      Func<Scalar>* u_prev_newton = u_ext[0];
      Func<Scalar>* u_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((u_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / TAU +
                          lam(u_prev_newton->val[i]) * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
		          - heat_src(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    // Jacobian for the implicit Euler time discretization
    template<typename Real, typename Scalar>
    Scalar J_euler(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      Func<Scalar>* u_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] / TAU + dlam_du(u_prev_newton->val[i]) * u->val[i] *
                           (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                           + lam(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      return result;
    }

Notice in the above code that the previous-level Newton iteration solution is obtained via::

    Func<Scalar>* u_prev_newton = u_ext[0];

Previous time level solution is accessed via the ExtData structure::

    Func<Scalar>* u_prev_time = ext->fn[0];

Physical coordinates of the i-th integration point are extracted from the Geom 
structure::

    x_phys = e->x[i];
    y_phys = e->y[i];

Weak forms for the Crank-Nicolson method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Crank-Nicolson method is an average of the implicit and explicit Euler methods,
thus of 

.. math::
    \int_{\Omega} \frac{u^{n+1}}{\tau}v \, \mbox{d}\bfx - 
    \int_{\Omega} \frac{u^{n}}{\tau}v \, \mbox{d}\bfx +
    \int_{\Omega} \lambda(u^{n+1})\nabla u^{n+1}\cdot \nabla v \, \mbox{d}\bfx -
    \int_{\Omega} fv \, \mbox{d}\bfx = 0.

and

.. math::
    \int_{\Omega} \frac{u^{n+1}}{\tau}v \, \mbox{d}\bfx - 
    \int_{\Omega} \frac{u^{n}}{\tau}v \, \mbox{d}\bfx +
    \int_{\Omega} \lambda(u^{n})\nabla u^{n}\cdot \nabla v \, \mbox{d}\bfx -
    \int_{\Omega} f v \, \mbox{d}\bfx = 0.

Weighting the last two relations with 1/2 and adding up, we get

.. math::
    \int_{\Omega} \frac{u^{n+1}}{\tau}v \, \mbox{d}\bfx - 
    \int_{\Omega} \frac{u^{n}}{\tau}v \, \mbox{d}\bfx +
    \int_{\Omega} \frac{1}{2}\lambda(u^{n+1})\nabla u^{n+1}\cdot \nabla v \, \mbox{d}\bfx +
    \int_{\Omega} \frac{1}{2}\lambda(u^{n})\nabla u^{n}\cdot \nabla v \, \mbox{d}\bfx -
    \int_{\Omega} fv \, \mbox{d}\bfx = 0.

The residual function $\bfF(\bfY)$ has the form 

.. math::
    F_i(\bfY) = \int_{\Omega} \frac{u(\bfY)}{\tau}v_i \, \mbox{d}\bfx - 
    \int_{\Omega} \frac{u^{n}}{\tau}v_i \, \mbox{d}\bfx +
    \int_{\Omega} \frac{1}{2}\lambda(u(\bfY))\nabla u(\bfY)\cdot \nabla v_i \, \mbox{d}\bfx +
    \int_{\Omega} \frac{1}{2}\lambda(u^n)\nabla u^n\cdot \nabla v_i \, \mbox{d}\bfx -
    \int_{\Omega} fv_i \, \mbox{d}\bfx = 0.

The Jacobian matrix $\bfJ(\bfY)$ has the form 

.. math::
    J_{ij}(\bfY) = \frac{\partial F_i}{\partial y_j}(\bfY) = \int_{\Omega} \frac{v_j}{\tau}v_i \, \mbox{d}\bfx +
    \int_{\Omega} \frac{1}{2}\frac{\partial \lambda}{\partial u}(u(\bfY))v_j \nabla u(\bfY)\cdot \nabla v_i \, 
    \mbox{d}\bfx +
    \int_{\Omega} \frac{1}{2}\lambda(u(\bfY))\nabla v_j\cdot \nabla v_i \, \mbox{d}\bfx.

In the code, this looks as follows::

    // Residual for the Crank-Nicolson time discretization
    template<typename Real, typename Scalar>
    Scalar F_cranic(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      Func<Scalar>* u_prev_newton = u_ext[0];
      Func<Scalar>* u_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((u_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / TAU
                           + 0.5 * lam(u_prev_newton->val[i]) * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                           + 0.5 * lam(u_prev_time->val[i]) * (u_prev_time->dx[i] * v->dx[i] + u_prev_time->dy[i] * v->dy[i])
                           - heat_src(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    // Jacobian for the Crank-Nicolson time discretization
    template<typename Real, typename Scalar>
    Scalar J_cranic(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      Func<Scalar>* u_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] / TAU +
                           0.5 * dlam_du(u_prev_newton->val[i]) * u->val[i] * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                           + 0.5 * lam(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      return result;
    }

Preparation for adaptivity
~~~~~~~~~~~~~~~~~~~~~~~~~~

At the beginning we convert the initial condition into a Solution::

    // Convert initial condition into a Solution.
    Solution sln_prev_time;
    sln_prev_time.set_exact(&mesh, init_cond);

In order to obtain an initial vector for the Newton's method, we have to project the 
initial condition on the FE space::

    // Project the initial condition on the FE space to obtain initial
    // coefficient vector for the Newton's method.
    info("Projecting initial condition to obtain initial vector for the Newton's method.");
    scalar* coeff_vec_coarse = new scalar[ndof];
    OGProjection::project_global(&space, &sln_prev_time, coeff_vec_coarse, matrix_solver);

Next we perform the Newton's method on the coarse mesh and translate the resulting 
coefficient vector into a Solution::

    // Newton's loop on the coarse mesh.
    info("Solving on coarse mesh:");
    bool verbose = true;
    if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, 
        NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solution(coeff_vec_coarse, &space, &sln);

Time stepping and periodic mesh derefinement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The time stepping loop begins with a periodic global mesh derefinement,
after which the last reference solution is projected on the globally
derefined mesh. The derefinement frequency is set by the user via the 
parameter UNREF_FREQ::

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);

      // Project on globally derefined mesh.
      info("Projecting previous fine mesh solution on derefined mesh.");
      OGProjection::project_global(&space, &sln_prev_time, &sln);
    }

The code above resets the actual mesh to the basemesh. Alternatively,
one could just remove a few layers of refinement (this is not so clean 
from the mathematical point of view but faster in practice). Speed 
optimization is not the main goal of the present example.

Adaptivity loop
~~~~~~~~~~~~~~~

The adaptivity loop begins by constructing a globally refined reference 
mesh::

    // Construct globally refined reference mesh
    // and setup reference space.
    Space* ref_space = construct_refined_space(&space);

In the first adaptivity step, a projection of the coarse mesh solution is used as 
an initial guess for the Newton's method on the reference mesh. Starting with the 
second adaptivity step, the previous reference mesh solution is projected instead::

    // Calculate initial coefficient vector for Newton on the fine mesh.
    if (as == 1) {
      info("Projecting coarse mesh solution to obtain coefficient vector on new fine mesh.");
      OGProjection::project_global(ref_space, &sln, coeff_vec, matrix_solver);
    }
    else {
      info("Projecting previous fine mesh solution to obtain coefficient vector on new fine mesh.");
      OGProjection::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver);
    }

Next we perform the Newton's loop on the reference mesh::

    // Newton's loop on the fine mesh.
    info("Solving on fine mesh:");
    if (!solve_newton(coeff_vec, dp, solver, matrix, rhs, 
		      NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

    // Store the result in ref_sln.
    Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);

The reference solution is projected on the coarse mesh for error calculation::

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

With the coarse and reference mesh approximations in hand, the coarse mesh is adapted 
as usual. At the end of each time step, the reference mesh solution is saved for the 
next time step::

    // Copy last reference solution into sln_prev_time.
    sln_prev_time.copy(&ref_sln);

Sample results
~~~~~~~~~~~~~~

TO BE CONTINUED.

