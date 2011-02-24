Gross-Pitaevski Equation (08-gross-pitaevski)
-----------------------------

**Git reference:** Tutorial example `08-gross-pitaevski 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P03-timedep/08-gross-pitaevski>`_.

In this example we use the Newton's method to solve the nonlinear complex-valued 
time-dependent Gross-Pitaevski equation. This equation describes the ground state of 
a quantum system of identical bosons using the Hartree–Fock approximation and the 
pseudopotential interaction model. For time-discretization one can use either
the first-order implicit Euler method or the second-order Crank-Nicolson
method. 

Model problem
~~~~~~~~~~~~~

The computational domain is the square $(-1,1)^2$ and boundary conditions are zero Dirichlet. The equation has the form 

.. math::

    i\hbar \frac{\partial \psi}{\partial t} = -\frac{\hbar^2}{2m} \Delta \psi + g \psi |\psi|^2 + \frac{m}{2} \omega^2 (x^2 + y^2) \psi

where $\psi(x,y)$ is the unknown solution (wave function), $i$ the complex unit, 
$\hbar$ the Planck constant, $m$ the mass of the boson, 
$g$ the coupling constant (proportional to the scattering length of two interacting bosons) and 
$\omega$ the frequency.

Complex-valued Jacobian matrix and residual vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We start by defining complex-valued  weak forms for the Newton's method::

    // Residuum for the implicit Euler time discretization
    template<typename Real, typename Scalar>
    Scalar F_euler(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

      Scalar result = 0;
      Func<Scalar>* psi_prev_newton = u_ext[0];
      Func<Scalar>* psi_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (ii * H * (psi_prev_newton->val[i] - psi_prev_time->val[i]) * v->val[i] / TAU
                - H*H/(2*M) * (psi_prev_newton->dx[i] * v->dx[i] + psi_prev_newton->dy[i] * v->dy[i])
                - G * psi_prev_newton->val[i] *  psi_prev_newton->val[i] * conj(psi_prev_newton->val[i]) * v->val[i]
                - .5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * psi_prev_newton->val[i] * v->val[i]);

      return result;
    }

    // Jacobian for the implicit Euler time discretization
    template<typename Real, typename Scalar>
    Scalar J_euler(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

      Scalar result = 0;
      Func<Scalar>* psi_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (ii * H * u->val[i] * v->val[i] / TAU
                         - H*H/(2*M) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                         - 2* G * u->val[i] *  psi_prev_newton->val[i] * conj(psi_prev_newton->val[i]) * v->val[i]
                         - G * psi_prev_newton->val[i] * psi_prev_newton->val[i] * u->val[i] * v->val[i]
                         - .5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * u->val[i] * v->val[i]);
      return result;
    }

    // Residuum for the Crank-Nicolson method
    template<typename Real, typename Scalar>
    Scalar F_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

      Scalar result = 0;
      Func<Scalar>* psi_prev_newton = u_ext[0];
      Func<Scalar>* psi_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (ii * H * (psi_prev_newton->val[i] - psi_prev_time->val[i]) * v->val[i] / TAU
                - 0.5*H*H/(2*M) * (psi_prev_newton->dx[i] * v->dx[i] + psi_prev_newton->dy[i] * v->dy[i])
                - 0.5*H*H/(2*M) * (psi_prev_time->dx[i] * v->dx[i] + psi_prev_time->dy[i] * v->dy[i])
                - 0.5*G * psi_prev_newton->val[i] *  psi_prev_newton->val[i] * conj(psi_prev_newton->val[i]) * v->val[i]
                - 0.5*G * psi_prev_time->val[i] *  psi_prev_time->val[i] * conj(psi_prev_time->val[i]) * v->val[i]
                - 0.5*0.5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * (psi_prev_newton->val[i] + psi_prev_time->val[i]) * v->val[i]);

      return result;
    }

    // Jacobian for the Crank-Nicolson method
    template<typename Real, typename Scalar>
    Scalar J_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

      Scalar result = 0;
      Func<Scalar>* psi_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (ii * H * u->val[i] * v->val[i] / TAU
                         - 0.5*H*H/(2*M) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                         - 0.5*2* G * u->val[i] *  psi_prev_newton->val[i] * conj(psi_prev_newton->val[i]) * v->val[i]
                         - 0.5*G * psi_prev_newton->val[i] *  psi_prev_newton->val[i] * u->val[i] * v->val[i]
                         - 0.5*.5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * u->val[i] * v->val[i]);
      return result;
    }

Registration of weak forms
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    // Initialize the weak formulation.
    WeakForm wf;
    if(TIME_DISCR == 1) {
      wf.add_matrix_form(callback(J_euler), HERMES_NONSYM, HERMES_ANY);
      wf.add_vector_form(callback(F_euler), HERMES_ANY, &psi_prev_time);
    }
    else {
      wf.add_matrix_form(callback(J_cranic), HERMES_NONSYM, HERMES_ANY);
      wf.add_vector_form(callback(F_cranic), HERMES_ANY, &psi_prev_time);
    }

Time stepping loop
~~~~~~~~~~~~~~~~~~

::

    // Time stepping loop:
    int nstep = (int)(T_FINAL/TAU + 0.5);
    for(int ts = 1; ts <= nstep; ts++)
    {
      info("Time step %d:", ts);

      // Perform Newton's iteration.
      info("Solving nonlinear problem:");
      bool verbose = true;
      if (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
          NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      // Update previous time level solution.
      Solution::vector_to_solution(coeff_vec, &space, &psi_prev_time);
    }

Sample results
~~~~~~~~~~~~~~

Snapshot 1:

.. image:: gross-pitaevski/sol_1.png
   :align: center
   :width: 600
   :alt: solution

Snapshot 2:

.. image:: gross-pitaevski/sol_2.png
   :align: center
   :width: 600
   :alt: solution

Snapshot 3:

.. image:: gross-pitaevski/sol_3.png
   :align: center
   :width: 600
   :alt: solution

