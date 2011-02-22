Systems of Equations (08-system)
-------------------------

**Git reference:** Tutorial example `08-system <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P01-linear/08-system>`_. 

So far we have just solved single linear PDE problems with a weak formulation
of the form $a(u,v) = l(v)$, where $u$ was a continuous approximation in the
$H^1$ space. Hermes can also solve equations whose solutions lie in the spaces
$Hcurl$, $Hdiv$ or $L^2$, and one can combine these spaces for PDE systems.

General scheme
~~~~~~~~~~~~~~

First let us understand how Hermes handles systems of linear PDE whose weak formulation 
is written as

.. math::
    :label: weaksystem

      a_{11}(u_1,v_1)\,+ a_{12}(u_2,v_1)\,+ \cdots\,+ a_{1n}(u_n,v_1) = l_1(v_1),

      a_{21}(u_1,v_2)\,+ a_{22}(u_2,v_2)\,+ \cdots\,+ a_{2n}(u_n,v_2) = l_2(v_2),

                                                          \vdots

      a_{n1}(u_1,v_n) + a_{n2}(u_2,v_n) + \cdots + a_{nn}(u_n,v_n) = l_n(v_n).

The solution $u = (u_1, u_2, \dots, u_n)$ and test functions $v =
(v_1, v_2, \dots, v_n)$ belong to the space $W = V_1 \times V_2 \times \dots
\times V_n$, where each $V_i$ is one of the available function spaces $H^1$, 
$H(curl)$, $H(div)$ or $L^2$. The resulting discrete matrix problem will have 
an $n \times n$ block structure.

Model problem of linear elasticity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us illustrate this by solving a simple problem of linear elasticity. Consider a
two-dimensional elastic body shown in the following figure (the bottom edge is
axis of planar symmetry):

.. image:: 08/elastsample.png
   :align: center
   :width: 500
   :height: 300
   :alt: Geometry and boundary conditions.

In the plane-strain model of linear elasticity the goal is to determine the
deformation of the body subject to the forces $f$. The deformation is sought
as a vector function $u(x) = (u_1, u_2)^T$, describing the displacement of each point
$x \in \Omega$ after the load $f = (f_1, f_2)^T$ is applied.

Boundary conditions
~~~~~~~~~~~~~~~~~~~

The boundary conditions are

.. math::
    :nowrap:

    \begin{eqnarray*}
    \frac{\partial u_1}{\partial n} &=& f_1 \ \text{on $\Gamma_3$,} \\
    \frac{\partial u_1}{\partial n} &=& 0 \ \text{on $\Gamma_2$, $\Gamma_4$, $\Gamma_5$,} \\
    \frac{\partial u_2}{\partial n} &=& f_2 \ \text{on $\Gamma_3$,} \\
    \frac{\partial u_2}{\partial n} &=& 0 \ \text{on $\Gamma_2$, $\Gamma_4$, $\Gamma_5$,} \\
    u_1 &=& u_2 = 0 \ \mbox{on} \ \Gamma_1. 
    \end{eqnarray*}

They are implemented as follows::

    // Boundary condition types.
    BCType bc_types(int marker)
      { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL;; }

    // Essential (Dirichlet) boundary condition values.
    scalar essential_bc_values(int ess_bdy_marker, double x, double y)
      { return 0; }

As usual, the Neumann boundary conditions will be incorporated into the 
weak formulation.

Displacement spaces
~~~~~~~~~~~~~~~~~~~

Next let us define function spaces for the two solution
components, $u_1$ and $u_2$ (the $x$ and $y$ displacement)::

    // Create x- and y- displacement spaces using default H1 shapesets.
    H1Space xdisp(&mesh, bc_types, essential_bc_values, P_INIT);
    H1Space ydisp(&mesh, bc_types, essential_bc_values, P_INIT);

Weak formulation
~~~~~~~~~~~~~~~~

Applying the standard procedure to the elastostatic equilibrium equations, we 
arrive at the following weak formulation:

.. math::
    :nowrap:

    \begin{eqnarray*}   \int_\Omega     (2\mu\!+\!\lambda)\dd{u_1}{x_1}\dd{v_1}{x_1} + \mu\dd{u_1}{x_2}\dd{v_1}{x_2} +     \mu\dd{u_2}{x_1}\dd{v_1}{x_2} + \lambda\dd{u_2}{x_2}\dd{v_1}{x_1}     \,\mbox{d}\bfx \!\!&=&\!\!\!     \int_{\Gamma_3} \!\!f_1 v_1 \,\mbox{d}S, \\ \smallskip   \int_\Omega     \mu\dd{u_1}{x_2}\dd{v_2}{x_1} + \lambda\dd{u_1}{x_1}\dd{v_2}{x_2} +     (2\mu\!+\!\lambda)\dd{u_2}{x_2}\dd{v_2}{x_2} + \mu\dd{u_2}{x_1}\dd{v_2}{x_1}     \,\mbox{d}\bfx \!\!&=&\!\!\!     \int_{\Gamma_3} \!\!f_2 v_2 \,\mbox{d}S. \end{eqnarray*}


We see that the weak formulation can indeed be written in the form :eq:`weaksystem`:

.. math::
    :nowrap:

    \begin{eqnarray*}
      a_{11}(u_1, v_1) \!&=&\! \int_\Omega (2\mu+\lambda)\dd{u_1}{x_1}\dd{v_1}{x_1} + \mu\dd{u_1}{x_2}\dd{v_1}{x_2} \,\mbox{d}\bfx,  \\
      a_{12}(u_2, v_1) \!&=&\! \int_\Omega \mu\dd{u_2}{x_1}\dd{v_1}{x_2} + \lambda\dd{u_2}{x_2}\dd{v_1}{x_1} \,\mbox{d}\bfx,\\
      a_{21}(u_1, v_2) \!&=&\! \int_\Omega \mu\dd{u_1}{x_2}\dd{v_2}{x_1} + \lambda\dd{u_1}{x_1}\dd{v_2}{x_2} \,\mbox{d}\bfx,\\
      a_{22}(u_2, v_2) \!&=&\! \int_\Omega (2\mu+\lambda)\dd{u_2}{x_2}\dd{v_2}{x_2} + \mu\dd{u_2}{x_1}\dd{v_2}{x_1} \,\mbox{d}\bfx,  \\
      l_{1}(v_1) \!&=&\!
      \int_{\Gamma_3} \!\!f_1 v_1 \,\mbox{d}S, \\
      l_{2}(v_2) \!&=&\!
      \int_{\Gamma_3} \!\!f_2 v_2 \,\mbox{d}S.
    \end{eqnarray*}

Here, $\mu$ and $\lambda$ are material constants (Lame coefficients) defined as

.. math::

    \mu = \frac{E}{2(1+\nu)}, \ \ \ \ \  \lambda = \frac{E\nu}{(1+\nu)(1-2\nu)},

where $E$ is the Young modulus and $\nu$ the Poisson ratio of the material. For
steel, we have $E = 200$ GPa and $\nu = 0.3$. The load is $f = (0, 10^4)^T$ N.

Block-wise WeakForm
~~~~~~~~~~~~~~~~~~~

The WeakForm instance is initialized for a system of two equations::

    // Initialize the weak formulation.
    WeakForm wf(2);
    wf.add_matrix_form(0, 0, callback(bilinear_form_0_0), HERMES_SYM);  // Note that only one symmetric part is
    wf.add_matrix_form(0, 1, callback(bilinear_form_0_1), HERMES_SYM);  // added in the case of symmetric bilinear
    wf.add_matrix_form(1, 1, callback(bilinear_form_1_1), HERMES_SYM);  // forms.
    wf.add_vector_form_surf(0, callback(linear_form_surf_0), GAMMA_3_BDY);
    wf.add_vector_form_surf(1, callback(linear_form_surf_1), GAMMA_3_BDY);

In the registration of matrix and vector forms,  
the block index 0, 0 means that bilinear_form_0_0() takes basis functions from 
space 0 (x-displacement space) and test functions from space 0. The block index 
0, 1 means that bilinear_form_0_1 takes basis functions from space 0 and test functions 
from space 1 (y-displacement space), etc. This yields a 2x2 block structure in the 
resulting matrix system.

Flags HERMES_SYM, HERMES_NONSYM, HERMES_ANTISYM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

Also explanation of the extra parameter HERMES_SYM in add_matrix_form() is in order.
Since the two diagonal forms $a_{11}$ and $a_{22}$ are symmetric, i.e.,
$a_{ii}(u,v) = a_{ii}(v,u)$, Hermes can be told to only evaluate them once for the
two cases $a_{ii}(u,v)$ and $a_{ii}(v,u)$ to speed up assembly. In fact, we should have
used the HERMES_SYM flag already in the previous sections, since the form
$a(u,v) = \nabla u \cdot \nabla v$ was symmetric. Of course this is not the case
for all forms and so the default value of the fourth parameter of add_matrix_form() 
is HERMES_NONSYM.

The off-diagonal forms $a_{12}(u_2, v_1)$ and $a_{21}(u_1, v_2)$ are not
(and cannot) be symmetric, since their arguments come from different spaces in general.
However, we can see that $a_{12}(u, v) = a_{21}(v, u)$, i.e., the corresponding blocks
of the local stiffness matrix are transposes of each other. Here, the HERMES_SYM flag
has a different effect: it tells Hermes to take the block of the local stiffness
matrix corresponding to the form $a_{12}$, transpose it and copy it where a block
corresponding to $a_{21}$ would belong, without evaluating $a_{21}$ at all (this is why
we don't add bilinear_form_1_0). This again speeds up the matrix assembly.
You can also use the flag HERMES_ANTISYM, which moreover inverts the sign of the block.
This makes sense in the case where $a_{ij}(u, v) = -a_{ji}(v, u)$.

It is recommended that you start with the default (and safe) HERMES_NONSYM flag for all
forms when developing your project, and only optimize the evaluation of the forms when
the code works well.

Assembling and solving the discrete problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the spaces and weak forms are ready, one can initialize the 
discrete problem::

    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem dp(&wf, Tuple<Space *>(&u_space, &v_space), is_linear);

Next we initialize the matrix solver::

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

And assemble and solve the matrix problem::

    // Assemble the stiffness matrix and right-hand side vector.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the linear system and if successful, obtain the solutions.
    info("Solving the matrix problem.");
    if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), Tuple<Space *>(&u_space, &v_space), 
                                                      Tuple<Solution *>(&u_sln, &v_sln));
    else error ("Matrix solver failed.\n");

Visualizing Von Mises stress
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Von Mises stress can be visualized via the VonMises filter as follows::

    // Visualize the solution.
    WinGeom* sln_win_geom = new WinGeom(0, 0, 800, 400);
    ScalarView view("Von Mises stress [Pa]", sln_win_geom);
    VonMisesFilter stress(Tuple<MeshFunction*>(&u_sln, &v_sln), lambda, mu);
    view.show_mesh(false);
    view.show(&stress, HERMES_EPS_HIGH, HERMES_FN_VAL_0, &u_sln, &v_sln, 1.5e5);

More about visualization and Filters will be said in the following section,
where we will also show sample results for the present model problem.
