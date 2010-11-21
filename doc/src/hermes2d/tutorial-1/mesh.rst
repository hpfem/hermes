Finite Element Mesh (01)
------------------------

**Git reference:** Tutorial example `01-mesh <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/01-mesh>`_. 

Every finite element computation starts with partitioning the domain
into a finite element mesh. Hermes uses (possibly curvilinear) triangles and 
quadrilaterals that can be combined together in one mesh. While non-adaptive
or low-order finite element codes need fine initial meshes constructed using 
specialized mesh generation software, in Hermes it usually suffices to
create a simple initial mesh and use a variety of built-in functions for 
a-priori mesh refinement. In most cases, automatic adaptivity will take 
care of the rest. 

.. image:: 01/simplemesh.png
   :align: center
   :width: 400
   :height: 400
   :alt: Sample finite element mesh.

The domain in this example is defined via four macroelements -- two
quadrilaterals and two curvilinear triangles. The elements are enumerated from 0 to 3. 
One also needs to enumerate all mesh vertices and assign markers to all boundary edges. 
Boundary markers are used to link boundary conditions with the boundary edges. 

Hermes2D mesh file format
~~~~~~~~~~~~~~~~~~~~~~~~~

Generic Hermes mesh file consists of variable assignments. Each variable can hold a real number, 
list of real numbers, or list of lists. The following are all valid definitions in 
the Hermes mesh file format::

    # comments start with a hash
    var = 5.0 + cos(pi)  # number
    list = { 1, 2, 3, 4, var }  # list
    pairs = { {1, 2}, {1, var}, {0, list} }  # list of lists

Every mesh file must contain at least the variables ``vertices``, ``elements``
and ``boundaries``. The variable ``vertices`` defines the coordinates
of all mesh vertices (in any order). For the above geometry it looks like this::

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

The variable ``elements`` defines all elements in the mesh via zero-based indices of their vertices in counter-clockwise order, plus an extra number denoting the element (material) marker. Element markers allow you to use different material parameters in areas with different material parameters. Moreover, Hermes allows you to assign different weak formulations to those areas, which can be very useful for some types of multiphysics problems. If the domain is composed of only one material, as the above geometry, all elements may be assigned a zero marker:
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
three numbers: two vertex indices and its marker. For the above geometry, we have
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

Finally, the mesh file can also include the variable ``curves`` which lists all
curved edges.  Each curved edge is described by one NURBS curve, defined by its
degree, control points and knot vector. Simplified syntax is available for
circular arcs.

NURBS curves
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
variable ``curves``, where $knots$ is a simple list of real values. 

Circular arcs
~~~~~~~~~~~~~

Circular arcs are very easy to define. For the above example, we have::

    curves =
    {
      { 4, 7, 45 },  # circular arcs with central angle 45 degrees
      { 7, 6, 45 }   # circular arcs with central angle 45 degrees
    }
    # EOF


Loading meshes in Hermes2D format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As a ''Hello world'' example, let us load the mesh we have just created, and display it in a window. 
Every main.cpp file in the git repository contains lots of comments and instructions. Skipping those, 
the `main.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/tutorial/01-mesh/main.cpp>`_ 
file begins with creating an instance of the class Mesh. In order to load
the mesh file, you have to create a mesh loader class (in our case that is ``H2DReader``) and
call the method ``load()``::

    #include "hermes2d.h"

    int main(int argc, char* argv[])
    {
      // Load the mesh file.
      Mesh mesh;
      H2DReader mloader;
      mloader.load("domain.mesh", &mesh);

Loading meshes in ExodusII format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hermes can read meshes in the `ExodusII <http://sourceforge.net/projects/exodusii/>`_ format.
This is a widely used format that can be generated, for example, 
with `Cubit <http://cubit.sandia.gov/>`_. To load an ExodusII mesh file, 
one has to use the ``ExodusIIReader`` class instead of the ``H2DReader`` class above.
We will encounter meshes in the ExodusII format in example 
`iron-water <http://hpfem.org/hermes/doc/src/hermes2d/examples/iron-water.html>`_. 

Manual mesh refinements
~~~~~~~~~~~~~~~~~~~~~~~

Below are examples of manual mesh refinements the user can do after loading the mesh.
All of them work for triangular, quadrilateral, and curvilinear elements. 

To begin with, here is how to refine element with index 'id'. If the element
is a quad, 0 means refine in both directions, 1 means refine
horizontally (with respect to the reference domain), 2 means refine vertically::

    void Mesh::refine_element(int id, int refinement = 0);

The mesh can be refined uniformly (multiple times if needed). The parameter 
'refinement' has the same meaning as in refine_element() above::

    void Mesh::refine_all_elements(int refinement = 0);

The mesh can be refined 'depth' times towards a vertex with index 'vertex_id'. In this
way a graded mesh towards the vertex is created::

    void Mesh::refine_towards_vertex(int vertex_id, int depth);

The following function performs repeated refinements of elements touching 
the boundary with boundary marker 'marker'. Elements touching with an 
edge or with a vertex are refined. 'aniso' allows or disables anisotropic
splitting of quads::

    void refine_towards_boundary(int marker, int depth, bool aniso = true);

The following will convert all quadrilateral elements in a triangular or 
triangular-quadrilateral mesh into triangles::

    void Mesh::convert_quads_to_triangles();

This will convert all triangular elements into quadrilaterals::

    void Mesh::convert_triangles_to_quads();

The following function selects elements to refine according to a given criterion and
performs 'depth' levels of refinements. The criterion function
receives a pointer to an element to be considered.
It must return -1 if the element is not to be refined, 0 if it
should be refined uniformly, 1 if it is a quad and should be split
horizontally or 2 if it is a quad and should be split vertically::

    void Mesh::refine_by_criterion(int (*criterion)(Element* e), int depth);

Meshes in Hermes can be arbitrarily irregular. The following function 
regularizes the mesh by refining elements with hanging nodes of
degree more than 'n'. As a result, n-irregular mesh is obtained.
If n = 0, completely regular mesh is created. In this case, however,
due to incompatible refinements, the element refinement hierarchy
is removed and all elements become top-level elements. Also, total
regularization does not work on curved elements. Returns an array of 
new element parents which can be passed to
Space::distribute_orders()::

    int* Mesh::regularize(int n);

The following function recursively removes all son elements 
of the given element and makes it active:: 

    Mesh::unrefine_element(int id);

All elements in the mesh can be unrefined using::

    Mesh::unrefine_all_elements();

See the file `src/mesh.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/src/mesh.cpp>`_ for more details. 

Visualizing the mesh
~~~~~~~~~~~~~~~~~~~~

The following code illustrates how to visualize the mesh using the MeshView class::

    // Display the mesh.
    // (0, 0) is the upper left corner position
    // 350 x 350 is the window size
    MeshView mview("Hello world!", new WinGeom(0, 0, 350, 350));
    mview.show(&mesh);

The class MeshView provides the method show() that displays a window showing the mesh:

.. image:: 01/meshview2.png
   :align: center
   :width: 400
   :height: 400
   :alt: Image of the mesh created via the MeshView class.

To see the graphical output, the main.cpp file should be finished with::

    // Wait for the view to be closed.
    View::wait();
    return 0;
  }
