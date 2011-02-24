Time-Integration with Runge-Kutta Methods (02-cathedral-rk)
-------------------------------------------------

**Git reference:** Tutorial example `02-cathedral-rk <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P03-timedep/02-cathedral-rk>`_. 

This example solves the same model problem as example `01-cathedral-ie <http://hpfem.org/hermes/doc/src/hermes2d/timedep/cathedral-ie.html>`_ but it shows how various Runge-Kutta methods can be used for time stepping. Let us begin with a brief introduction 
to the Runge-Kutta methods and Butcher's tables before we explain implementation details.

Runge-Kutta methods and Butcher's tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Runge-Kutta methods are very nicely described on their `Wikipedia page <http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods>`_
and we recommend the reader to have a brief look there before reading further. In particular, read about the representation 
of Runge-Kutta methods via Butcher's tables, and learn how the table looks like if the method is explicit, diagonally-implicit,
and fully implicit. You should understand basic properties of Runge-Kutta methods such as that explicit methods always need
very small time step or they will blow up, that implicit methods can use much larger time steps, and that among implicit methods, 
diagonally-implicit ones are especially desirable because of relatively low computatonal cost.

Butcher's tables currently available in Hermes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a list of predefined Runge-Kutta methods that can be found 
in the file `hermes_common/tables.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes_common/tables.cpp>`_.
The names of the tables are self-explanatory. The last number is the order of the 
method (a pair of orders for embedded ones). The second-to-last, if provided, is the number of stages.

Explicit methods::

* Explicit_RK_1, 
* Explicit_RK_2, 
* Explicit_RK_3, 
* Explicit_RK_4. 

Implicit methods::

* Implicit_RK_1, 
* Implicit_Crank_Nicolson_2_2, 
* Implicit_SIRK_2_2, 
* Implicit_ESIRK_2_2, 
* Implicit_SDIRK_2_2, 
* Implicit_Lobatto_IIIA_2_2, 
* Implicit_Lobatto_IIIB_2_2, 
* Implicit_Lobatto_IIIC_2_2, 
* Implicit_Lobatto_IIIA_3_4, 
* Implicit_Lobatto_IIIB_3_4, 
* Implicit_Lobatto_IIIC_3_4, 
* Implicit_Radau_IIA_3_5, 
* Implicit_SDIRK_4_5.

Embedded explicit methods::

* Explicit_HEUN_EULER_2_12_embedded, 
* Explicit_BOGACKI_SHAMPINE_4_23_embedded, 
* Explicit_FEHLBERG_6_45_embedded,
* Explicit_CASH_KARP_6_45_embedded, 
* Explicit_DORMAND_PRINCE_7_45_embedded.

Embedded implicit methods::

* Implicit_SDIRK_CASH_3_23_embedded,
* Implicit_SDIRK_BILLINGTON_3_23_embedded,
* Implicit_ESDIRK_TRBDF2_3_23_embedded, 
* Implicit_ESDIRK_TRX2_3_23_embedded,
* Implicit_SDIRK_CASH_5_24_embedded,
* Implicit_SDIRK_CASH_5_34_embedded,
* Implicit_DIRK_FUDZIAH_7_45_embedded. 

Plus, the user is free to define any Butcher's table of his own.

Model problem
~~~~~~~~~~~~~

The model problem is the same as in example
`01-cathedral-ie <http://hpfem.org/hermes/doc/src/hermes2d/timedep/cathedral-ie.html>`_ 
However, for the purpose of using Runge-Kutta methods, the equation has to be 
formulated in such a way that the time derivative stands solo on the left-hand side and 
everything else is on the right

.. math::
    :label: eqvit1

       \frac{\partial T}{\partial t} = \frac{\lambda}{c \varrho} \Delta T.

Weak formulation
~~~~~~~~~~~~~~~~

The weak formulation is only needed for the right-hand side

.. math::

     F(T) = - \int_{\Omega} \frac{\lambda}{c \varrho} \nabla T\cdot \nabla v
            + \int_{\Gamma_{air}} \frac{\alpha \lambda}{c \varrho} (T_{ext}(t) - T)v.

This is different from example `01-cathedral-ie <http://hpfem.org/hermes/doc/src/hermes2d/timedep/cathedral-ie.html>`_ 
where the discretization of the time derivative term was part of the weak formulation. The approach presented
here makes the method more modular.

The function $F$ above is the stationary residual of the equation (i.e., the weak form of the right-hand side).
Since the Runge-Kutta equations are solved using the Newton's method, the reader may want to have a brief 
look ahead to the `Newton's method section <http://hpfem.org/hermes/doc/src/hermes2d/nonlinear/newton.html>`_.
Then it will be easy to see that the weak form of the Jacobian matrix of the stationary residual is

.. math::

     \frac{\partial F_i(Y)}{\partial y_j} = - \int_{\Omega} \frac{\lambda}{c \varrho} \nabla v_j\cdot \nabla v_i 
                  - \int_{\Gamma_{air}} \frac{\alpha \lambda}{c \varrho} v_j v_i.

Defining weak forms
~~~~~~~~~~~~~~~~~~~

Bilinear and linear forms are defined as usual::

    template<typename Real, typename Scalar>
    Scalar stac_jacobian_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
			     Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
	result += -wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      }

      return result * LAMBDA / HEATCAP / RHO;
    }

    template<typename Real, typename Scalar>
    Scalar stac_residual_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
		             Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Func<Scalar>* u_prev = u_ext[0];

      Scalar result = 0;
      for (int i = 0; i < n; i++) {
	result += -wt[i] * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i]);
	result += wt[i] * heat_src(e->x[i], e->y[i]) * v->val[i];	       
      }

      return result * LAMBDA / HEATCAP / RHO;
    }

    template<typename Real, typename Scalar>
    Scalar stac_jacobian_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
			      Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return - LAMBDA / HEATCAP / RHO * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    template<typename Real, typename Scalar>
    Scalar stac_residual_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
			      Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Func<Scalar>* u_prev = u_ext[0];

      // This is a temporary workaround. The stage time t_n + h * c_i
      // can be accessed via u_stage_time->val[0];
      Func<Scalar>* u_stage_time = ext->fn[0]; 

      Scalar stage_time = u_stage_time->val[0];
      Real stage_ext_temp = temp_ext<Real>(stage_time);

      Scalar result = 0;
      for (int i = 0; i < n; i++) {
	result += wt[i] * (stage_ext_temp - u_prev->val[i]) * v->val[i];		       
      }

      return LAMBDA / HEATCAP / RHO * ALPHA * result;
    }
 
The previous-level solution is accessed via::

    Func<Scalar>* u_prev = u_ext[0];

and the stage time as::

  Scalar stage_time = u_stage_time->val[0];

The latter is a temporary solution and it will be replaced in due course
by passing a real number as it ought to be.

Selecting a Butcher's table
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unless the user wants to define a Butcher's table on his/her own, he/she can select 
a predefined one - for example a second-order diagonally implicit SDIRK-22
method::

    ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

This is followed in main.cpp by creating an instance of the table::

    ButcherTable bt(butcher_table_type);

Registering weak forms
~~~~~~~~~~~~~~~~~~~~~~

The weak forms are registered as follows::

    // Initialize weak formulation.
    WeakForm wf;
    wf.add_matrix_form(callback(stac_jacobian_vol));
    wf.add_vector_form(callback(stac_residual_vol));
    wf.add_matrix_form_surf(callback(stac_jacobian_surf), BDY_AIR);
    wf.add_vector_form_surf(callback(stac_residual_surf), BDY_AIR);

Setting initial condition
~~~~~~~~~~~~~~~~~~~~~~~~~ 

Before time stepping, one needs to obtain the coefficient vector of the initial
condition::

    // Project the initial condition on the FE space to obtain initial solution coefficient vector.
    info("Projecting initial condition to translate initial condition into a vector.");
    scalar* coeff_vec = new scalar[ndof];
    OGProjection::project_global(&space, &u_prev_time, coeff_vec, matrix_solver);

Initializing the discrete problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The discrete problem is initialized with is_linear = false (the default value), 
disregarding whether it is linear or not::

    // Initialize the FE problem.
    bool is_linear = false;
    DiscreteProblem dp(&wf, &space, is_linear);

Time-stepping loop
~~~~~~~~~~~~~~~~~~

Finally, the time-stepping loop takes the form::

    // Time stepping loop:
    double current_time = 0.0; int ts = 1;
    do 
    {
      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      info("Runge-Kutta time step (t = %g, tau = %g, stages: %d).", 
           current_time, time_step, bt.get_size());
      bool verbose = true;
      bool is_linear = true;
      if (!rk_time_step(current_time, time_step, &bt, coeff_vec, &dp, matrix_solver,
	  	        verbose, is_linear)) {
        error("Runge-Kutta time step failed, try to decrease time step size.");
      }

      // Convert coeff_vec into a new time level solution.
      Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);

      // Update time.
      current_time += time_step;

      // Show the new time level solution.
      char title[100];
      sprintf(title, "Time %3.2f, exterior temperature %3.5f", current_time, temp_ext(current_time));
      Tview.set_title(title);
      Tview.show(&u_prev_time);

      // Increase counter of time steps.
      ts++;
    } 
    while (current_time < T_FINAL);


