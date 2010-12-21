Using NURBS (37)
----------------

**Git reference:** Tutorial example `37-nurbs <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/37-nurbs>`_. 

This example illustrates how to use full-featured NURBS
to define curved boundary edges. Recall that simplified 
format is available for circular arcs, as was shown 
in example 03-poisson. 

Solved is a simple Poisson equation with constant right-hand
side and homogeneous Dirichlet boundary conditions.

The domain is a rectangle (0,2) x (0, 1) whose upper
edge is a NURBS curve. There are three mesh files
in this example: domain-1.mesh (one control point),
domain-2.mesh (two control points), and domain-3.mesh
(three control points). One of these files needs to be 
selected on line 15 in main.cpp::

    const char* mesh_file = "domain-1.mesh";          // One control point.
    //const char* mesh_file = "domain-2.mesh";          // Two control points.
    //const char* mesh_file = "domain-3.mesh";          // Three control points.

**Example of a NURBS with one inner control point**

Snippet from the file domain-1.mesh::

    degree = 2                # Degree should be equal to the 
    num_inner_points = 1      # number of inner points plus one.
                              
    inner_points =
    {
      { 1.5, 2.0, 2.0 }       # x, y, weight
    } 
    knots = 
    {
      0, 0, 0, 1, 1, 1        
    }
    curves =
    {
      {2, 3, degree, inner_points, knots} 
    }

Result:

.. image:: 37/1.png
   :align: center
   :width: 500
   :alt: NURBS with one control point.

**Example of a NURBS with two inner control points**

Snippet from the file domain-3.mesh::

    degree = 3
    num_inner_points = 2
    inner_points =
    {
      { 1.5, 1.5, 1.0 },
      { 0.5, 0.5, 1.0 }
    } 
    knots = 
    {
      0, 0, 0, 1, 1, 1
    }
    curves =
    {
      {2, 3, degree, inner_points, knots} 
    }

Result:

.. image:: 37/2.png
   :align: center
   :width: 500
   :alt: NURBS with two control points.


**Example of a NURBS with three inner control points**

Snippet from the file domain-2.mesh::

    degree = 4
    num_inner_points = 3
    inner_points =
    {
      { 1.5, 1.5, 1.0 },
      { 1.0, -1.0, 1.0 },
      { 0.5, 1.5, 1.0 }
    } 
    knots = 
    {
      0, 0, 0, 1, 1, 1
    }
    curves =
    {
      {2, 3, degree, inner_points, knots} 
    }

Result:

.. image:: 37/3.png
   :align: center
   :width: 500
   :alt: NURBS with three control points.




