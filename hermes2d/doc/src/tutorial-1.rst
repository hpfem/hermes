=================================
Tutorial Part I (Linear Problems)
=================================

This tutorial should give you a good idea of how Hermes2D works. After reading it, you will
be able to create your own applications and/or adjust existing Hermes examples for your 
purposes. At the beginning of every section we give a reference to the corresponding example in the 
Hermes git repository -- there you will always find the corresponding main.cpp file, mesh file(s) etc.

This document is under continuous development. If you find bugs, typos, dead links or such,
let us know through the `mailing list <http://groups.google.com/group/hermes2d/>`_.

Finite Element Mesh (01)
--------------------------------

**Git reference:** Tutorial example `01-mesh <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/01-mesh>`_. 

Every finite element computation starts with partitioning the domain
into a finite element mesh. Hermes uses triangles and quadrilaterals, and 
can combine both element types in one mesh. While complicated meshes need 
to be constructed using specialized mesh generation software, in many cases 
we only need a simple initial mesh that can be created by hand. In Hermes, all you 
need to do is partition the domain very coarsely into several large elements.
Then you can use a number of elementary mesh refinement functions and/or
let the automatic adaptivity algorithm take care of the rest. 

.. image:: img/tutorial-01/simplemesh.png
   :align: center
   :width: 400
   :height: 400
   :alt: Sample finite element mesh.

The domain in this example is defined via four macroelements -- two
quadrilaterals and two curvilinear triangles. The elements are enumerated from 0 to 3. 
One also needs to enumerate all mesh vertices and assign markers to all boundary edges. 
Boundary markers are used to link boundary conditions with the boundary edges. 

Mesh File Format
~~~~~~~~~~~~~~~~

Hermes can read meshes in its own generic format as well as in the
`ExodusII <http://sourceforge.net/projects/exodusii/>`_ format
(this is, for example, the output of `Cubit <http://cubit.sandia.gov/>`_).
First let us discuss the generic Hermes mesh data format. Reading
of ExodusII mesh files is very simple as we will see in example 
`iron-water <http://hpfem.org/hermes2d/doc/src/examples.html#iron-water-neutronics>`_. 

Generic Hermes mesh file consists of variable assignments. Each variable can hold a real number, 
list of real numbers, or list of lists. The following are all valid definitions in 
the Hermes mesh file format::

    # comments start with a hash
    var = 5.0 + cos(pi)  # number
    list = { 1, 2, 3, 4, var }  # list
    pairs = { {1, 2}, {1, var}, {0, list} }  # list of lists

Every mesh file must contain at least the variables ``vertices``, ``elements``
and ``boundaries``. The variable ``vertices`` defines the coordinates
of all mesh vertices (in any order). In our case it looks like this::

    a = 1.0  # size of the mesh
    b = sqrt(2)/2

    vertices =
    {
      { 0, -a },    # vertex 0
      { a, -a },    # vertex 1
      { -a, 0 },    # vertex 2
      { 0, 0 },     # vertex 3
      { a, 0 },     # vertex 4
      { -a, a },    # vertex 5
      { 0, a },     # vertex 6
      { a*b, a*b }  # vertex 7
    }

The variable ``elements`` defines all elements in the mesh via zero-based indices of their vertices in counter-clockwise order, plus an extra number denoting the element (material) marker. Element markers allow you to use different material parameters in areas with different material parameters. Moreover, Hermes allows you to assign different weak formulations to those areas, which can be very useful for some types of multiphysics problems. If the domain is composed of only one material, as in our case, all elements may be assigned a zero marker:
::

    elements =
    {
      { 0, 1, 4, 3, 0 },  # quad 0
      { 3, 4, 7, 0 },     # tri 1
      { 3, 7, 6, 0 },     # tri 2
      { 2, 3, 6, 5, 0 }   # quad 3
    }

The last mandatory variable, ``boundaries``, defines boundary markers for all
boundary edges. By default, all edges have zero markers. Only those with
positive markers are considered to be part of the domain boundary and can be
assigned a boundary condition, as we will see later. An edge is identified by
two vertex indices. In our case, we have
::

    boundaries =
    {
      { 0, 1, 1 },
      { 1, 4, 2 },
      { 3, 0, 4 },
      { 4, 7, 2 },
      { 7, 6, 2 },
      { 2, 3, 4 },
      { 6, 5, 2 },
      { 5, 2, 3 }
    }

Finally, the file can also include the variable ``curves``, which lists all
curved edges.  Each curved edge is described by one NURBS curve, defined by its
degree, control points and knot vector. Simplified syntax is available for
circular arcs.

NURBS Curves
~~~~~~~~~~~~

Every NURBS curve is defined by its degree, control points with weights and the
knot vector. The degree $d$ is a positive integer, usually 1, 2, 3 or 5. Lines
and polylines are of degree 1, circles have degree 2 and free-form curves are
of degree 3 or 5. The control points $p_i$, $i = 0 \ldots n$, are the main tool for changing the
shape of the curve. A curve of degree $d$ must have at least $d+1$ control
points. In Hermes, the endpoints of the edge are always assumed to be the
first and last control points and therefore only the inner control points are
listed in the mesh file. There is a weight $w_i \geq 0$ for every control point,
that influences the shape of the curve in its vicinity. If $w_i = 0$ then 
$p_i$ has no effect on the shape.  As $w_i$ increases, the curve is pulled 
towards $p_i$.

The knot vector is a sequence of $m+1$ values that determines how much and
where the control points influence the shape. The relation $m = n+d+1$ must
hold. The sequence is nondecreasing, $t_i \leq t_{i+1}$, and divides the whole
interval $[0,1]$ into smaller intervals which determine the area of influence
of the control points. Since the curve has to start and end at the edge
vertices, the knot vector in Hermes always starts with $d+1$ zeros and ends
with $d+1$ ones. Only the inner knots are listed in the above definition of the
variable ``curves``, where $knots$ is a simple list of real values. For the 
above example, we have
::

    curves =
    {
      { 4, 7, 45 },  # +45 degree circular arcs
      { 7, 6, 45 }
    }
    # EOF


Loading Mesh
~~~~~~~~~~~~

As a ''Hello world'' example, let us load the mesh we have just created, and display it in a window. 
Every main.cpp file in the git repository contains lots of comments and instructions. Skipping those, 
the `main.cpp <http://git.hpfem.org/hermes2d.git/blob/HEAD:/tutorial/01-mesh/main.cpp>`_ 
file begins with creating an instance of the class Mesh. In order to load
the mesh file, you have to create a mesh loader class (in our case that is ``H2DReader``) and
call the method ``load()``:
::

    #include "hermes2d.h"

    int main(int argc, char* argv[])
    {
      // load the mesh file
      Mesh mesh;
      H2DReader mloader;
      mloader.load("domain.mesh", &mesh);

Note: To load the exodus-II mesh file, one has to use ``ExodusIIReader`` class instead.

The following portion of code illustrates various types of initial mesh refinements.
It does not matter if the mesh becomes irregular, in fact, arbitrarily irregular
meshes are at the heart of Hermes: 
::

      // perform some sample initial refinements
      mesh.refine_all_elements();          // refines all elements
      mesh.refine_towards_vertex(3, 4);    // refines mesh towards
                                           // vertex #3 (4x)
      mesh.refine_towards_boundary(2, 4);  // refines all elements
                                           // along boundary 2 (4x)
      mesh.refine_element(86, 0);          // refines element #86
                                           // isotropically
      mesh.refine_element(112, 0);         // refines element #112
                                           // isotropically
      mesh.refine_element(84, 2);          // refines element #84
                                           // anisotropically
      mesh.refine_element(114, 1);         // refines element #114
                                           // anisotropically

Other ways of modifying meshes on the fly include
::

    Mesh::refine_element(int id, int refinement = 0);
    Mesh::convert_quads_to_triangles();
    Mesh::convert_triangles_to_quads();
    Mesh::refine_by_criterion(int (*criterion)(Element* e), int depth);
    Mesh::refine_towards_vertex(int vertex_id, int depth);
    Mesh::regularize(int n);
    Mesh::unrefine_element(int id);
    Mesh::unrefine_all_elements();

See the file `src/mesh.cpp <http://git.hpfem.org/hermes2d.git/blob/HEAD:/src/mesh.cpp>`_ for more details. 
The following code illustrates how to visualize the mesh using the class MeshView:
::

    // display the mesh
    // (100, 100) is the upper left corner position
    // 500 x 500 is the window size
    MeshView mview("Hello world!", 100, 100, 500, 500);
    mview.show(&mesh);

You can initialize it by supplying the title of the window and its initial position and size (all of these
parameters are optional). The class MeshView provides the method show() that displays a window showing the mesh:

.. image:: img/tutorial-01/meshview2.png
   :align: center
   :width: 400
   :height: 400
   :alt: Image of the mesh created via the MeshView class.

Every main.cpp file is finished with 
::

    // wait for keyboard or mouse input
    View::wait();
    return 0;
  }

so that you have a chance to see the graphical output.



Setting Up Finite Element Space (02)
------------------------------------

**Git reference:** Tutorial example `02-space <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/02-space>`_. 

Hermes follows the mathematical concept of FEM closely -- after creating a mesh,
in the next step one needs to construct a finite element space on it.
The following predefined spaces are currently available:

* H1Space - the most common space of continuous, piecewise-polynomial functions belonging to $H^1(\Omega) = \{ v \in L^2(\Omega); \nabla u \in [L^2(\Omega)]^2 \}$,
* HcurlSpace - the space of vector-valued functions discontinuous along mesh edges, with continuous tangential component on the edges $H(\mbox{curl},\Omega) = \{ E \in [L^2(\Omega)]^2; \nabla \times E \in L^2(\Omega)\}$,
* HdivSpace - the space of vector-valued functions discontinuous along mesh edges, with continuous normal component on the edges $H(\mbox{div},\Omega) = \{ v \in [L^2(\Omega)^2; \nabla \cdot v \in L^2(\Omega)\}$,
* L2Space -  the space of functions discontinuous along mesh edges, belonging to the space $L^2(\Omega)$.

All these spaces allow for higher-order elements and meshes with arbitrary-level hanging nodes.
If you are not familiar with higher-order FEM, let us just say that the spaces can contain
quadratic, cubic, etc., *edge functions* that generate higher-degree
polynomials along mesh edges, and *bubble functions* that complete the higher-order
approximation in element interiors. Edge functions are associated with mesh edges,
and bubble functions with element interiors. The next figure shows a patch consisting of two triangular elements. An edge function is shown on the left, and a bubble function on one of the triangles on the right:

.. image:: img/tutorial-02/basisfn.jpg
   :align: center
   :width: 600
   :height: 200
   :alt: Fourth-order edge function  (left) and one of the fifth-order bubble functions (right).

There are many possible ways of defining the
higher-order basis functions. A particular set of polynomials is called
*shapeset*. Using good shapeset is crucial for the
performance of the *hp*-FEM. No shapeset can be optimal for all possible operators.
Therefore, Hermes offers several shapesets from which
you need to choose one when creating a FE space. The ones which perform best
in most computations (according to our experience) are simply called
H1Shapeset, HcurlShapeset, HdivShapeset and L2Shapeset.
Others can be found in the files `src/shapeset* <http://git.hpfem.org/hermes2d.git/tree/HEAD:/src>`_ in the git repo.
Any shapeset can be used for more than one space.

We are now ready for an example. The following is (up to some comments) the complete
`main.cpp <http://git.hpfem.org/hermes2d.git/blob/HEAD:/tutorial/02-space/main.cpp>`_ file
of the example 02-space::

    #include "hermes2d.h"
    int P_INIT = 3;
    int main(int argc, char* argv[])
    {
      // Load the mesh.
      Mesh mesh;
      H2DReader mloader;
      mloader.load("domain.mesh", &mesh);

      // Create an H1 space with default shapeset and natural BC.
      H1Space space(&mesh, NULL, NULL, P_INIT);

      // View FE basis functions.
      BaseView bview("FE Space", 0, 0, 600, 600);
      bview.show(&space);

      // Wait for the view to be closed.
      View::wait();
      return 0;
    }

An instance of H1Space is initialized with four arguments: 

* pointer to a mesh, 
* function providing the type of boundary condition for various boundary markers 
  (NULL means natural boundary conditions on the entire boundary),
* function providing values of essential boundary conditions (not relevant for natural BC),
* uniform initial polynomial degree of all mesh elements.

If only linear elements are used, then the initialization of the $H^1$ space is even simpler::

    // Create an H1 space with default shapeset,
    // natural BC, and linear elements.
    H1Space space(&mesh);

The polynomial degree of elements can also be set individually by calling 
the method Space::set_element_order() or for all elements at once using
Space::set_uniform_order(). Note that element degrees
are stored in Space, not in Mesh. The reason is that in Hermes one can
have multiple spaces with different element degrees and even types 
over the same mesh. In Hermes, Mesh only stores geometrical information.
A space created in this way is ready for use. 

As a debugging/learning feature, Hermes can visualize the basis of each Space.
Similarly to MeshView, one can create a BaseView object and use it 
to display the entire basis (VectorBaseView has to be used for vector-valued 
approximations in spaces Hcurl and Hdiv - this will be discussed later). 
One can cycle through all basis functions in the window using the arrow keys. 
If you press the left mouse button at the beginning, you will see the Dirichlet 
lift (a function that represents Dirichlet boundary conditions).

This is how the last figure above was obtained (press the '3' key for 3D mode).
We suggest that you spend some time experimenting with element refinements and 
hanging nodes to see how basis functions on irregular meshes look like.

Solving Poisson Equation (03)
-----------------------------

**Git reference:** Tutorial example `03-poisson <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/03-poisson>`_. 

Let us solve the Poisson equation

.. math::
    :label: poisson1

       -\Delta u = CONST_F

on the L-shaped domain $\Omega$ from the previous example,
equipped with homogeneous (zero) Dirichlet boundary conditions

.. math::
    :label: poisson2

       u = 0\ \ \  \mbox{on}\  \partial \Omega,

where $CONST_F$ is a real number. The weak formulation 
is derived in the standard way, first by multiplying equation :eq:`poisson1` with a test
function $v$, then integrating over the domain $\Omega$, and then applying the Green's
theorem (integration by parts) to the second derivatives.
Because of the homogeneous Dirichlet condition :eq:`poisson2`,
the proper space for the solution is $V = H^1_0(\Omega)$. The weak formulation reads:
Find $u \in V$ such that

.. math::
    :label: poissonweak

         \int_\Omega \nabla u \cdot \nabla v \;\mbox{d\bfx} = CONST_F \int_\Omega v \;\mbox{d\bfx} \ \ \ \mbox{for all}\ v \in V.

Equation :eq:`poissonweak` has the standard form $a(u,v) = l(v)$ and thus in Hermes
we need a way to define the bilinear form $a(u,v)$ and the linear form $l(v)$.
This is done by implementing the following two functions:
::

    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext);

    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext);

These functions are called for each element during the stiffness matrix
assembly and must return the values of the bilinear and linear forms for the given arguments.
RealFunction represents one of the basis functions restricted to the
current element and RefMap represents the reference mapping of the current element.
There are methods for extracting the values of the basis functions at integration points,
which allows you to evaluate the integrals by yourself, but this is normally not needed,
since many common weak forms have already been implemented.
In this case, we can simply use the predefined functions
int_grad_u_grad_v and int_v:
::

    // Return the value \int \nabla u . \nabla v dx.
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }
   
    // Return the value \int v dx.
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return CONST_F * int_v<Real, Scalar>(n, wt, v);
    }

Later we will learn how to compose arbitrary integrals using function values and derivatives, and integration points and weights. The weak forms are registered as follows::

    // Initialize the weak formulation.
    WeakForm wf();
    wf.add_matrix_form(callback(bilinear_form));
    wf.add_vector_form(callback(linear_form));

Later we will learn how to register Jacobian and residual forms for nonlinear problems. If the PDE is more complicated, we can add multiple matrix and vector forms.

With the space and weak formulation in hand, the problem is solved via::

    // Solve the linear problem.
    Solution sln;
    solve_linear(&space, &wf, &sln, SOLVER_UMFPACK);

The parameter SOLVER_UMFPACK indicates that we are using the direct sparse matrix solver UMFpack. Other options include SOLVER_PETSC, SOLVER_MUMPS, a variety of SciPy matrix solvers and others - the choice of matrix solvers will be discussed in more detail later. 

The solution can be visualized via the ScalarView class::

    // Visualize the solution.
    WinGeom* sln_win_geom = new WinGeom(0, 0, 440, 350);
    ScalarView view("Solution", sln_win_geom);
    view.show(&sln);

The following figure shows the output of this example (again, press '3' for 3D view).

.. image:: img/tutorial-03/poisson.png
   :align: center
   :width: 400
   :height: 350
   :alt: Solution of the Poisson equation.

Short and Long Versions of Examples
-----------------------------------

Most tutorial examples come in two versions: A short one that is intended for effortless basic use, and a long one that is more explicit and thus more convenient for development. The first example with a long version is 03-poisson.

**Git reference:** Tutorial example `03-poisson-long <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/03-poisson-long>`_. 

The long version does not employ the function solve_linear(). Instead, after initializing the weak formulation, one initializes the LinearProblem class::

      // Initialize the linear problem.
      LinearProblem lp(&wf, &space);

This class is a descendant of a more general DiscreteProblem class that handles nonlinear problems. Next we initialize the matrix solver and the corresponding matrix and vector structures::

      // Select matrix solver.
      Matrix* mat; Vector* rhs; CommonSolver* solver;
      init_matrix_solver(SOLVER_UMFPACK, ndof, mat, rhs, solver);

Again, other matrix solvers besides SOLVER_UMFPACK can be used. The variable ndof stands for the number of degrees of greedom (unknowns in the discrete problem) that can be calculated after initializing a Space::

      int ndof = get_num_dofs(&space);

Assembling is done into the user-provided data structures::

      // Assemble stiffness matrix and rhs.
      lp.assemble(mat, rhs);

After this, the matrix problem is solved::

      // Solve the matrix problem.
      if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

And finally, the solution vector is translated into a Solution::

      // Convert coefficient vector into a Solution.
      Solution sln;
      sln.set_fe_solution(&space, rhs);

Visualization and the rest of the main() function are the same as in the short version.

Boundary Conditions (04, 05, 06)
--------------------------------

Hermes recognizes two basic types of boundary conditions: *essential* and *natural*.
Essential boundary conditions (prescribed values on the boundary) influence the finite element 
space while natural conditions do not - they are incorporated into boundary integrals in the weak formulation.
In the context of elliptic problems, Dirichlet conditions are essential and Neumann/Newton
conditions are natural.

Examples 04, 05 and 06 also come in long versions but we will not discuss them explicitly since they are analogous to the long version of example 03.

Dirichlet BC
~~~~~~~~~~~~

**Git reference:** Tutorial example `04-bc-dirichlet <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/04-bc-dirichlet>`_. Long version: `04-bc-dirichlet-long <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/04-bc-dirichlet-long>`_. 

Since essential boundary conditions eliminate degrees of freedom (DOF) from the FE space, 
they need to be incorporated while the space is set up.
The user has to provide the following two callback functions::

    BCType bc_types(int marker);
    scalar essential_bc_values(int ess_bdy_marker, double x, double y);

The first one takes as argument a boundary marker number, and it determines the type of BC 
for the corresponding portion of the domain boundary, by returning one of the predefined constants 
BC_ESSENTIAL, BC_NATURAL. The second callback needs to return the boundary value for a given marker
and position on the boundary (only needed for essential boundary condition markers - for natural
boundary conditions this value is ignored). The space initialization then consists of the following 
line::

    H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

Here P_INIT is the initial polynomial degree of all elements in the mesh as before. 
Suppose that we would like to modify the boundary conditions for the previous Poisson 
model problem as follows:

.. math::
         u(x,y) = -\frac{CONST_F}{4}(x^2 + y^2)\,\ \mbox{on}\,\ \partial \Omega.

This is done by defining

::

    BCType bc_types(int marker)
    {
      return BC_ESSENTIAL;
    }

and setting the essential BC values callback to return the value of the Dirichlet BC::

    scalar essential_bc_values(int ess_bdy_marker, double x, double y)
    {
      return (-CONST_F/4)*(x*x + y*y);
    }

It is easy to see that the solution to this problem is the function

.. math::
         u(x,y) = -\frac{CONST_F}{4}(x^2 + y^2). 

For the value $CONST_F = -4$, the output is shown below:

.. image:: img/tutorial-04/dirichlet.png
   :align: center
   :width: 400
   :height: 350
   :alt: Solution of the Dirichlet problem.

Neumann BC
~~~~~~~~~~

**Git reference:** Tutorial example `05-bc-neumann <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/05-bc-neumann>`_. Long version: `05-bc-neumann-long <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/05-bc-neumann-long>`_.

Next, let us consider Neumann boundary conditions. The new model problem
will have the form

.. math::
    :nowrap:

    \begin{eqnarray*}   -\Delta u = CONST_F,\ \ \ \ \ &&u = 0\,\ \mbox{on}\,\ \Gamma_4,\\                            &&\dd{u}{n} = C_1\,\ \mbox{on}\,\ \Gamma_1,\\                            &&\dd{u}{n} = C_2\,\ \mbox{on}\,\ \Gamma_2,\\                            &&\dd{u}{n} = C_3\,\ \mbox{on}\,\ \Gamma_3. \end{eqnarray*}

where $\Gamma_1 \dots \Gamma_4$ correspond to the edges marked $1 \dots 4$. Now, the weak formulation contains some surface integrals:

.. math::

    \int_\Omega \nabla u \cdot \nabla v \;\mbox{d\bfx} =   CONST_F\int_\Omega v \;\mbox{d\bfx}   + C_1\int_{\Gamma_1} \!v \;\mbox{d}l   + C_2\int_{\Gamma_2} \!v \;\mbox{d}l   + C_3\int_{\Gamma_3} \!v \;\mbox{d}l


In Hermes, all forms in the standard weak formulation $a(u,v) = l(v)$
are in fact defined as a sum of contributions from volume integrals and from
surface integrals. In the case of the linear form $l(v)$, this means

.. math::

    l(v) = \sum_m l_m^{\,\rm vol}(v) + \sum_n l_n^{\,\rm surf}(v).

We have already seen volumetric linear forms in example 
`03-poisson <http://hpfem.org/hermes2d/doc/src/tutorial-1.html#solving-poisson-equation-03>`_. 
Surface linear forms are implemented similarly. Our new right-hand side is
represented by two functions with the following prototypes::

    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    
    template<typename Real, typename Scalar>
    Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext);

and registered as follows::

    // Initialize the weak formulation
    WeakForm wf();
    wf.add_matrix_form(callback(bilinear_form));
    wf.add_vector_form(callback(linear_form));
    wf.add_vector_form_surf(callback(linear_form_surf));

The surface linear form is defined as::

    template<typename Real, typename Scalar>
    Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return CONST_GAMMA[e->marker - 1] * int_v<Real, Scalar>(n, wt, v);
    }

Here, we have used the predefined surface integral int_v (see the
file `src/integrals_h1.h <http://git.hpfem.org/hermes2d.git/blob/HEAD:/src/integrals_h1.h>`_). 
If the boundary conditions were more complicated, we could also
have used int_F_v, where F stands for an arbitrary user-supplied
function returning the value $\partial u/\partial n$.

Note that in this example, the mesh is a-priori refined towards the re-entrant corner 
to capture the singular gradient::

    mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);  // '3' is the vertex index from the mesh file.

The gradient magnitude can be visualized via a MagFilter::

    // Compute and show gradient magnitude
    // (note that the infinite gradient at the re-entrant
    // corner will be truncated for visualization purposes)
    ScalarView gradview("Gradient", grad_win_geom);
    MagFilter grad(Tuple<MeshFunction>(&sln, &sln), Tuple<int>(H2D_FN_DX, H2D_FN_DY));
    gradview.show(&grad);

Here we first meet Tuple - a construction designed to avoid variable argument 
lists. The first Tuple is used to pass a pair of pointers to the same MeshFunction,
and the next Tuple says that the vector components for the magnitude calculation 
are the x- and y- partial derivatives. The class Solution that represents a piecewise-polynomial
finite element function on a Mesh, is descendant of a more general class MeshFunction
that can represent constants, general functions given via an analytic formula, 
finite element solutions, etc. 

The approximate solution for the values $C_1 = -1/2$, $C_2 = 1$, $C_3 = -1/2$,
along with the singularity of gradient at the re-entrant corner are
shown in the following figures:

.. image:: img/tutorial-05/neumann2.png
   :align: left
   :width: 530
   :height: 400
   :alt: Solution of the Neumann problem.

.. image:: img/tutorial-05/neumann3.png
   :align: right
   :width: 400
   :height: 400
   :alt: Detail of gradient singularity at the re-entrant corner.

.. raw:: html

   <hr style="clear: both; visibility: hidden;">

Newton BC
~~~~~~~~~

**Git reference:** Tutorial example `06-bc-newton <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/06-bc-newton>`_. Long version: `06-bc-newton-long <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/06-bc-newton-long>`_.

Another common natural boundary condition is the Newton (sometimes called Robin) condition
of the form

.. math::

    \dd{u}{n} + c_1 u = c_2, \ \ \ \ c_1 \ne 0.

Analogously to Neumann conditions, also Newton conditions yield surface integrals. However,
this time they are both in the bilinear form and in the linear form,
The bilinear form is
a sum of volume and surface forms that can be added to the weak formulation using the methods
add_matrix_form() and add_matrix_form_surf(). 
The surface bilinear form must have the following prototype:
::

    template<typename Real, typename Scalar>
    Scalar bilinear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext);

Inside this function you can use predefined
forms such as int_u_v, int_F_u_v (see the
file `src/integrals_h1.h <http://git.hpfem.org/hermes2d.git/blob/HEAD:/src/integrals_h1.h>`_) or your custom forms.

The following code snippet contains the linear and bilinear forms:
::

    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    template<typename Real, typename Scalar>
    Scalar bilinear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return H * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    template<typename Real, typename Scalar>
    Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return T0 * H * int_v<Real, Scalar>(n, wt, v);
    }

  

Here, $T_0$ is the exterior temperature, and $H$ is the heat flux.
The above forms are registered using::

    // Initialize the weak formulation.
    WeakForm wf;
    wf.add_matrix_form(callback(bilinear_form));
    wf.add_matrix_form_surf(callback(bilinear_form_surf), NEWTON_BDY);
    wf.add_vector_form_surf(callback(linear_form_surf), NEWTON_BDY);

Here NEWTON_BDY is the boundary marker for the Newton boundary. The following figures 
show the solution and singularity of gradient at the re-entrant corner:

.. image:: img/tutorial-06/newton1.png
   :align: left
   :width: 530
   :height: 400
   :alt: Solution of the Newton problem.

.. image:: img/tutorial-06/newton2.png
   :align: right
   :width: 400
   :height: 400
   :alt: Detail of gradient singularity at the re-entrant corner.

.. raw:: html

   <hr style="clear: both; visibility: hidden;">

Determination of Quadrature Orders in Weak Forms
------------------------------------------------

You may wonder why templates are used in the definition of weak forms. As a matter of fact, 
they do not have to be, as we will see in a moment. However, if the weak form only contains 
algebraic operations (without if-then statements and such), templates help to determine
numerical integration orders automatically. In higher-order FEM, basis and test functions may 
have very different polynomial degrees, ranging from one and some maximum polynomial 
degree (currently 10 in Hermes). The basis and test functions can be combined inside the 
weak forms in many different ways. As a result, the minimum quadrature order which is needed 
to evaluate a weak form accurately may vary between zero (product of gradients of 
two linear functions) to infinity (whenever a nonpolynomial expression is present). 
Numerical quadrature is one of the trickiest issues in higher-order FEM.

A brute-force solution to this problem would be to integrate everything using 
a maximum order, but this would lead to tremendous computing times. Therefore Hermes offers 
two options: the polynomial degree of the integrated expressions can be detected 
automatically (via templates), or the user can define for each weak form the 
quadrature order explicitly. If the weak form only contains polynomial expressions, 
the former approach works very well. If the form is more complicated, it is recommended 
to handle the integration orders explicitly. 

Automatic determination of quadrature order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In example 03-poisson, the bilinear and linear forms were defined using templates,

::

    // return the value \int \nabla u . \nabla v dx
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    // return the value \int v dx
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return CONST_F * int_v<Real, Scalar>(n, wt, v);
    }

and registered using the callback() macro,

::

    // initialize the weak formulation
    WeakForm wf();
    wf.add_matrix_form(callback(bilinear_form));
    wf.add_vector_form(callback(linear_form));
   
The callback() macro, defined in `src/forms.h 
<http://git.hpfem.org/hermes2d.git/blob/HEAD:/src/forms.h>`_ by

::

    #define callback(a)     a<double, scalar>, a<Ord, Ord>

expands the above add_matrix_form() and add_vector_form() functions into

::

    // initialize the weak formulation
    WeakForm wf();
    wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>);
    wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>);

For those who are not familiar with templates, they make it possible to 
call the same function with different parameter types. In particular, 
using bilinear_form<double, scalar> and bilinear_form<Ord, Ord> for
the bilinear form defined above gives 

::

    scalar bilinear_form(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return int_grad_u_grad_v<double, scalar>(n, wt, u, v);
    }

    Ord bilinear_form(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
    }

The <double, scalar> copy is used to obtain the result of the numerical integration,
the <Ord, Ord> copy for automatic evaluation of the quadrature order. 
The parser (see `src/forms.h 
<http://git.hpfem.org/hermes2d.git/blob/HEAD:/src/forms.h>`_) 
works well for algebraic expressions. If the weak form bilinear_form() is complicated, 
one can create and register a simpler weak form bilinear_form_order() for the parser,
that provides an arbitrary expression with the same polynomial degree as 
the integrand in bilinear_form(). Then the two functions would be registered as 

::

    wf.add_matrix_form(bilinear_form, bilinear_form_order);

Of course the same holds for linear forms.
If the bilinear form contains things like the if-then statement, it cannot 
be parsed. Whenever the weak form contains non-polynomial expressions or
is otherwise very complicated, it is recommended to handle the quadrature 
orders manually.

Manual determination of quadrature order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The polynomial degree of basis and test functions inside a bilinear or linear form 
can be handled manually as follows

::

    Ord bilinear_form_order(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      int uo = u->val[0].get_order();
      int vo = v->val[0].get_order();
      return Ord(uo + vo);            // this would correspond to integral of u times v
    }

It is also possible to return a constant order (for example 5) by using 

::

    Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return Ord(5);
    }

Currently, one cannot make the integration order dependent on spatial coordinates and such. However,
one can assign different weak forms to elements with different material markers. This is
described in examples `iron-water <http://git.hpfem.org/hermes2d.git/tree/HEAD:/examples/iron-water>`_,
`saphir <http://git.hpfem.org/hermes2d.git/blob/HEAD:/examples/saphir/main.cpp>`_ and others.

The following example handles quadrature orders manually. 

General 2nd-Order Linear Equation (07)
--------------------------------------

**Git reference:** Tutorial example `07-general <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/07-general>`_. Long version: `07-general-long <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/07-general-long>`_.

This example deals with a linear second-order equation of the form 

.. math::

         -\frac{\partial}{\partial x}\left(a_{11}(x,y)\frac{\partial u}{\partial x}\right) - \frac{\partial}{\partial x}\left(a_{12}(x,y)\frac{\partial u}{\partial y}\right) - \frac{\partial}{\partial y}\left(a_{21}(x,y)\frac{\partial u}{\partial x}\right) - \frac{\partial}{\partial y}\left(a_{22}(x,y)\frac{\partial u}{\partial y}\right) + a_1(x,y)\frac{\partial u}{\partial x} + a_{21}(x,y)\frac{\partial u}{\partial y} + a_0(x,y)u = rhs(x,y),

equipped with Dirichlet and/or Neumann boundary conditions. Its goal is to show how to 
use space-dependent coefficients and how to define quadrature orders explicitly. 

First we define the (generally) non-constant equation coefficients:
::

    double a_11(double x, double y) {
      if (y > 0) return 1 + x*x + y*y;
      else return 1;
    }

and so on. Then we define boundary conditions as usual. The weak formulation contains
both volumetric and surface integrals. 

The Ord class in Hermes (see the file `src/forms.h 
<http://git.hpfem.org/hermes2d.git/blob/HEAD:/src/forms.h>`_) provides
an automatic parser of weak forms that is able to determine the integration orders for 
algebraic expressions. So, in order to define an integration order explicitly, one can 
provide on top the weak form another function that defines a simple algebraic expression 
that leads the parser to the desired polynomial degree. The values defined in this  
additional function are not used for computation. 

::

    // (Volumetric) bilinear form
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      for (int i=0; i < n; i++) {
        double x = e->x[i];
        double y = e->y[i];
        result += (a_11(x, y)*u->dx[i]*v->dx[i] + 
                   a_12(x, y)*u->dy[i]*v->dx[i] +
                   a_21(x, y)*u->dx[i]*v->dy[i] +
                   a_22(x, y)*u->dy[i]*v->dy[i] +
                   a_1(x, y)*u->dx[i]*v->val[i] +
                   a_2(x, y)*u->dy[i]*v->val[i] +
                   a_0(x, y)*u->val[i]*v->val[i]) * wt[i];
      }
      return result;
    }

    // Integration order for the bilinear form
    Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return u->val[0] * v->val[0] * e->x[0] * e->x[0]; // returning the sum of the degrees of the basis 
                                                        // and test function plus two
    }

    // Surface linear form (natural boundary conditions)
    template<typename Real, typename Scalar>
    Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_F_v<Real, Scalar>(n, wt, g_N, v, e);
    }
  
    // Integration order for surface linear form
    Ord linear_form_surf_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
    }
  
    // Volumetric linear form (right-hand side)
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
    }
  
    // Integration order for the volumetric linear form
    Ord linear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
    }

Note the sign of the surface linear form - when using the LinearProblem class, all linear forms have to be on the right-hand side and all bilinear forms on the left. 

The output of this example is shown below:

.. image:: img/tutorial-07/general.png
   :align: center
   :width: 500
   :height: 400
   :alt: Output of example 07-general.

Systems of Equations (08)
-------------------------

**Git reference:** Tutorial example `08-system <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/08-system>`_. Long version `08-system-long <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/08-system-long>`_.

So far we have just solved single linear PDE problems with a weak formulation
of the form $a(u,v) = l(v)$, where $u, v$ were continuous approximations in the
$H^1$ space. One can also solve equations whose solutions lie in the spaces
$Hcurl$, $Hdiv$ or $L^2$, and one can combine these spaces for PDE systems.

Here we show how to handle systems of linear PDE whose weak formulation is written as

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

Let us illustrate this by solving a simple problem of linear elasticity. Consider a
two-dimensional elastic body shown in the following figure (the bottom edge is
axis of planar symmetry):

.. image:: img/tutorial-08/elastsample.png
   :align: center
   :width: 500
   :height: 300
   :alt: Geometry and boundary conditions.

In the plane-strain model of linear elasticity the goal is to determine the
deformation of the body subject to the forces $f$. The deformation is sought
as a vector function $u(x) = (u_1, u_2)^T$, describing the displacement of each point
$x \in \Omega$ after the load $f = (f_1, f_2)^T$ is applied.


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

Applying the standard procedure to the elastostatic equilibrium equations, we arrive at the following weak formulation:

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

We begin with defining the function spaces for the two solution
components, $u_1$ and $u_2$ (the $x$ and $y$ displacement). The boundary
conditions can be implemented as follows::

    // Boundary condition types.
    BCType bc_types(int marker)
      { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL;; }

    // Essential (Dirichlet) boundary condition values.
    scalar essential_bc_values(int ess_bdy_marker, double x, double y)
      { return 0; }

Next we create two displacement spaces::

    // Create x- and y- displacement spaces using default H1 shapesets.
    H1Space xdisp(&mesh, bc_types, essential_bc_values, P_INIT);
    H1Space ydisp(&mesh, bc_types, essential_bc_values, P_INIT);

The WeakForm instance is initialized for a system of two equations::

    // initialize the weak formulation
    WeakForm wf(2);
    wf.add_matrix_form(0, 0, callback(bilinear_form_0_0), H2D_SYM);  // Note that only one symmetric part is
    wf.add_matrix_form(0, 1, callback(bilinear_form_0_1), H2D_SYM);  // added in the case of symmetric bilinear
    wf.add_matrix_form(1, 1, callback(bilinear_form_1_1), H2D_SYM);  // forms.
    wf.add_vector_form_surf(0, callback(linear_form_surf_0), GAMMA_3_BDY);
    wf.add_vector_form_surf(1, callback(linear_form_surf_1), GAMMA_3_BDY);

In the registration of matrix and vector forms,  
the block index 0, 0 means that bilinear_form_0_0() takes basis functions from 
space 0 (x-displacement space) and test functions from space 0. The block index 
0, 1 means that bilinear_form_0_1 takes basis functions from space 0 and test functions 
from space 1 (y-displacement space), etc. This yields a 2x2 block structure in the 
resulting matrix system.

Also explanation of the extra parameter H2D_SYM in add_matrix_form() is in order.
Since the two diagonal forms $a_{11}$ and $a_{22}$ are symmetric, i.e.,
$a_{ii}(u,v) = a_{ii}(v,u)$, Hermes can be told to only evaluate them once for the
two cases $a_{ii}(u,v)$ and $a_{ii}(v,u)$ to speed up assembly. In fact, we should have
used the H2D_SYM flag already in the previous sections, since the form
$a(u,v) = \nabla u \cdot \nabla v$ was symmetric. Of course this is not the case
for all forms and so the default value of the fourth parameter of add_matrix_form() 
is H2D_UNSYM.

The off-diagonal forms $a_{12}(u_2, v_1)$ and $a_{21}(u_1, v_2)$ are not
(and cannot) be symmetric, since their arguments come from different spaces in general.
However, we can see that $a_{12}(u, v) = a_{21}(v, u)$, i.e., the corresponding blocks
of the local stiffness matrix are transposes of each other. Here, the H2D_SYM flag
has a different effect: it tells Hermes to take the block of the local stiffness
matrix corresponding to the form $a_{12}$, transpose it and copy it where a block
corresponding to $a_{21}$ would belong, without evaluating $a_{21}$ at all (this is why
we don't add bilinear_form_1_0). This again speeds up the matrix assembly.
You can also use the flag H2D_ANTISYM, which moreover inverts the sign of the block.
This makes sense in the case where $a_{ij}(u, v) = -a_{ji}(v, u)$.

It is recommended that you start with the default (and safe) H2D_UNSYM flag for all
forms when developing your project, and only optimize the evaluation of the forms when
the code works well.

When the spaces and weak forms are ready, one can use the function solve_linear() to
assemble and solve the discrete problem::

    // Solve the linear problem.
    Solution u_sln, v_sln;
    solve_linear(Tuple<Space *>(&u_space, &v_space), &wf, 
                 Tuple<Solution*>(&u_sln, &v_sln), matrix_solver);

Von Mises stress can be visualized via the VonMises filter as follows::

    // Visualize the solution.
    WinGeom* sln_win_geom = new WinGeom(0, 0, 800, 400);
    ScalarView view("Von Mises stress [Pa]", sln_win_geom);
    VonMisesFilter stress(Tuple<MeshFunction*>(&u_sln, &v_sln), lambda, mu);
    view.show_mesh(false);
    view.show(&stress, H2D_EPS_HIGH, H2D_FN_VAL_0, &u_sln, &v_sln, 1.5e5);

We will say more about visualization and Filters in a moment, after showing the long version of this example.

Long Version of Example 08
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Git reference:** Tutorial example `08-system-long <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/08-system-long>`_.

As in example 03, the long version of this example does not employ the function solve_linear(). Instead, after initializing the weak formulation, one initializes the LinearProblem class, selects a matrix solver, assembles the matrix problem, solves it, and translates the resulting coefficient vector into Solutions::

    // Initialize the linear problem.
    LinearProblem lp(&wf, Tuple<Space *>(&u_space, &v_space));

    // Select matrix solver.
    Matrix* mat; Vector* rhs; CommonSolver* solver;
    init_matrix_solver(matrix_solver, ndof, mat, rhs, solver);

    // Assemble stiffness matrix and rhs.
    lp.assemble(mat, rhs);

    // Solve the matrix problem.
    if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

    // Convert coefficient vector into a Solution.
    Solution u_sln, v_sln;
    u_sln.set_fe_solution(&u_space, rhs);
    v_sln.set_fe_solution(&v_space, rhs);

Visualization and Filters
-------------------------

In elasticity problems one often wants to see the material
stress, which is obtained by a formula that combines the derivatives 
of the two displacement components.
Hermes implements postprocessing through Filters. Filter is a special class
which takes up to three Solutions, performs some computation and in the end acts
as another Solution (which can be visualized, passed into another Filter,
passed into a weak form, etc.). More advanced usage of Filters will be discussed 
later. In elasticity examples we typically use the predefined VonMisesFilter::

    VonMisesFilter stress(&u_sln, &v_sln, lambda, mu);
    view.show_mesh(false);
    view.show(&stress, H2D_EPS_HIGH);

The second line tells Hermes not to display mesh edges.
The second parameter of show() is the visualization accuracy. It can have the 
values H2D_EPS_LOW, H2D_EPS_NORMAL (default) and H2D_EPS_HIGH. This parameter 
influences the number of linear triangles that Hermes uses to approximate 
higher-order polynomial solutions within finite elements. Using linear 
triangles is required by OpenGL, so Hermes at least performs automatic 
adaptivity to reduce their number to a minimum. The above parameters
are used to set the accuracy of this piecewise-linear approximation. 

The method show() has an optional third parameter to indicate whether 
function values or partial derivatives should be displayed. For example,
H2D_FN_VAL_0 stands for the function value of solution component 0
(first solution component which in this case is the VonMises stress).
H2D_FN_VAL_1 would mean the function value of the second solution component
(relevant for vector-valued $Hcurl$ or $Hdiv$ elements only), 
H2D_FN_DX_0 means the x-derivative of the first solution component, etc.

Finally, in elasticity problems it may be desirable to deform the computational
domain according to the calculated displacements. The method View::show() has
additional three optional parameters for this::

    VonMisesFilter stress(Tuple<MeshFunction*>(&u_sln, &v_sln), mu, lambda);
    view.show(&stress, H2D_EPS_HIGH, H2D_FN_VAL_0, &u_sln, &v_sln, 1.5e5);

Here the fourth and fifth parameters are the displacement components used to 
distort the domain geometry, and the sixth parameter is a scaling factor to multiply the 
displacements. Of course, the color map still shows the Von Mises stress as before. 

.. image:: img/tutorial-08/mises.png
   :align: center
   :width: 550
   :height: 300
   :alt: Elastic stress plotted on deformed domain.

Time-Dependent Problems (09)
----------------------------

**Git reference:** Tutorial example `09-timedep <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/09-timedep>`_. 

This section describes the implementation of a simple time-dependent
heat transfer model that describes, in a naive approximation, how the St. Vitus cathedral
in Prague responds to changes in the surrounding air temperature
during one 24-hour cycle. The geometry is shown below:

.. image:: img/tutorial-09/vitus1.png
   :align: center
   :width: 400
   :height: 500
   :alt: Model geometry and temperature distribution after 24 hours.

We will solve the standard heat transfer equation

.. math::
    :label: eqvit1

       c \varrho\frac{\partial T}{\partial t} - \lambda \Delta T = 0

equipped with a Dirichlet condition

.. math::

     T = T_{init}

on the bottom edge $\Gamma_{ground}$ and a Newton condition

.. math::

     \frac{\partial T}{\partial \nu} = \alpha(T_{ext}(t) - T)

on the rest of the boundary $\Gamma_{air}$. Here, $c$ is the heat capacity of the material,
$\varrho$ the material density, $\lambda$ the thermal conductivity,
$T_{init}$ the fixed temperature on the
ground (same as the initial temperature of the building), and $\alpha$
the heat transfer coefficient 
between the building and the surrounding air. The surrounding air temperature
$T_{ext}$ is time-dependent of the form

.. math::

     T_{ext}(t) = T_{init} + 10\sin(2\pi t/T_{final}),

where $T_{final}$ is 24 hours (translated into seconds).

Equation :eq:`eqvit1` is also equipped with an initial condition of the
form

.. math::

     T(x,y,0) = T_{init}(x,y) \ \ \ \mbox{in} \ \Omega.



For simplicity we will use the implicit Euler method with a constant
time step $\tau$, which transforms equation :eq:`eqvit1` into


.. math::

     c \varrho\frac{T^{n+1} - T^n}{\tau} - \lambda \Delta T^{n+1} = 0.

The corresponding weak formulation is

.. math::

     \int_{\Omega} c \varrho\frac{T^{n+1}}{\tau} + \int_{\Omega} \lambda \nabla T^{n+1}\cdot \nabla v + \int_{\Gamma_{air}} \alpha \lambda T^{n+1}v = \int_{\Omega} c \varrho\frac{T^{n}}{\tau} + \int_{\Gamma_{air}} \alpha \lambda T_{ext}(t^{n+1})v.

The implementation starts by defining the
boundary condition types::

    BCType bc_types(int marker)
    {
      if (marker == marker_ground) return BC_ESSENTIAL;
      else return BC_NATURAL;
    }

and values::

    scalar essential_bc_values(int ess_bdy_marker, double x, double y)
    {
      if (ess_bdy_marker == marker_ground) return T_INIT;
    }

Then the space for the temperature $T$ is set up::

    // Initialize an H1 space with default shepeset.
    H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
    info("ndof = %d", get_num_dofs(&space));

Bilinear and linear forms are defined as follows::

    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, u, v) / TAU +
             LAMBDA * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }
  
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, ext->fn[0], v) / TAU;
    }
  
    template<typename Real, typename Scalar>
    Scalar bilinear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return LAMBDA * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
    }
  
    template<typename Real, typename Scalar>
    Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return LAMBDA * ALPHA * temp_ext(TIME) * int_v<Real, Scalar>(n, wt, v);
    }

Next we need to initialize the previous solution tsln with the initial condition $T_{init}$.
Besides holding the finite element solution, the Solution class
can be forced to return zero, to return a constant, or to return an arbitrary function
using the methods set_zero(), set_const() and set_exact(), respectively.
Here we simply call set_const() and supply the initial temperature::

    // Set constant initial condition.
    Solution tsln;
    tsln.set_const(&mesh, T_INIT);

The weak forms are registered as follows::

    // Initialize weak formulation.
    WeakForm wf();
    wf.add_matrix_form(bilinear_form<double, double>, bilinear_form<Ord, Ord>);
    wf.add_matrix_form_surf(bilinear_form_surf<double, double>, bilinear_form_surf<Ord, Ord>, marker_air);
    wf.add_vector_form(linear_form<double, double>, linear_form<Ord, Ord>, H2D_ANY, 1, &tsln);
    wf.add_vector_form_surf(linear_form_surf<double, double>, linear_form_surf<Ord, Ord>, marker_air);

Next, the LinearProblem class and the matrix solver structures are initialized::

    // Initialize the linear problem.
    LinearProblem lp(&wf, &space);

    // Initialize matrix solver.
    Matrix* mat; Vector* rhs; CommonSolver* solver;  
    init_matrix_solver(matrix_solver, ndof, mat, rhs, solver);

We are now ready to start the iterative process. Since the stiffness matrix does
not depend on the solution, it only needs to be assembled once in the first time
step. For all remaining time steps it will be the same, and we just need to
re-construct the load vector. This is done via the Boolean variable rhsonly
which is set to false before the time stepping begins. For completeness, we show 
the entire time stepping loop below::

    // Time stepping:
    int nsteps = (int)(FINAL_TIME/TAU + 0.5);
    bool rhsonly = false;
    for(int ts = 1; ts <= nsteps; ts++)
    {
      info("---- Time step %d, time %3.5f, ext_temp %g", ts, TIME, temp_ext(TIME));

      // Assemble stiffness matrix and rhs.
      lp.assemble(mat, rhs, rhsonly);
      rhsonly = true;

      // Solve the matrix problem.
      if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

      // Update tsln.
      tsln.set_fe_solution(&space, rhs);

      // Update the time variable.
      TIME += TAU;

      // Visualize the solution.
      sprintf(title, "Time %3.2f, exterior temperature %3.5f", TIME, temp_ext(TIME));
      Tview.set_title(title);
      Tview.show(&tsln);
    }



