Finite Element Space (02-space)
-------------------------

**Git reference:** Tutorial example `02-space <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P01-linear/02-space>`_. 

Hermes follows the mathematical concept of FEM closely -- after creating a mesh,
in the next step one needs to construct a finite element space on it.

Spaces available in Hermes
~~~~~~~~~~~~~~~~~~~~~~~~~~

The following predefined spaces are currently available:

* H1Space - the most common space of continuous, piecewise-polynomial functions belonging to $H^1(\Omega) = \{ v \in L^2(\Omega); \nabla u \in [L^2(\Omega)]^2 \}$,
* HcurlSpace - space of vector-valued functions discontinuous along mesh edges, with continuous tangential component on the edges $H(\mbox{curl},\Omega) = \{ E \in [L^2(\Omega)]^2; \nabla \times E \in L^2(\Omega)\}$,
* HdivSpace - space of vector-valued functions discontinuous along mesh edges, with continuous normal component on the edges $H(\mbox{div},\Omega) = \{ v \in [L^2(\Omega)^2; \nabla \cdot v \in L^2(\Omega)\}$,
* L2Space - space of functions discontinuous along mesh edges, belonging to the space $L^2(\Omega)$.

All these spaces allow for higher-order elements, meshes with arbitrary-level hanging nodes,
and automatic *hp*-adaptivity. 

Structure of higher-order basis functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are not familiar with higher-order FEM, let us just say that the spaces can contain
quadratic, cubic, etc., *edge functions* that generate higher-degree
polynomials along mesh edges, and *bubble functions* that complete the higher-order
approximation in element interiors. Edge functions are associated with mesh edges,
and bubble functions with element interiors. The next figure shows a higher-order  
edge function (left) and a higher-order bubble function (right). 

.. image:: 02/basisfn.jpg
   :align: center
   :width: 600
   :height: 200
   :alt: Fourth-order edge function  (left) and one of the fifth-order bubble functions (right).

Higher-order basis functions can be defined in many different ways. 
A particular set of polynomials is called *shapeset*. Using a good shapeset is crucial for the
performance of the *hp*-FEM. No shapeset can be optimal for all possible operators.
Therefore, Hermes offers several shapesets from which
you need to choose one when creating a FE space. The ones which perform best
in most computations (according to our experience) are simply called
H1Shapeset, HcurlShapeset, HdivShapeset and L2Shapeset.
Others can be found in the directory `src/shapeset/ <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/src/shapeset>`_. 

Example
~~~~~~~

Up to some comments, the following is the complete
`main.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/tutorial/P01-linear/02-space/main.cpp>`_ file
of the example 02-space::

    #include "hermes2d.h"
    int P_INIT = 3;
    int main(int argc, char* argv[])
    {
      // Load the mesh.
      Mesh mesh;
      H2DReader mloader;
      mloader.load("domain.mesh", &mesh);

      // Enter boundary markers.
      // (If no markers are entered, default is a natural BC).
      BCTypes bc_types;

      // Enter Dirichlet boundary values (default is zero).
      BCValues bc_values;

      // Create an H1 space with default shapeset and natural BC.
      H1Space space(&mesh, &bc_types, &bc_values, P_INIT);

      // View FE basis functions.
      BaseView bview("FE Space", new WinGeom(0, 0, 440, 350));
      bview.show(&space);

      // Wait for the view to be closed.
      View::wait();
      return 0;
    }

Initializing a Space
~~~~~~~~~~~~~~~~~~~~

An instance of H1Space is initialized with four arguments: 

* Pointer to a mesh. 
* Pointer to class BCTypes that links boundary markers to boundary condition types.
  If a marker is not linked to any type, then the default is Neumann. This does not 
  really matter in this example since no PDE is solved. Usage of the BCTypes class 
  will be explained in more detail in the tutorial examples 04, 05 and 06. 
* Pointer to class BCValues that provides values of Dirichlet boundary conditions. 
  Again this does not matter in this example since no PDE is solved, and details will
  be given in the tutorial examples 04, 05 and 06.
* Uniform initial polynomial degree of all mesh elements.

Setting element orders individually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Polynomial degrees of elements can also be set individually by calling 
the method Space::set_element_order() or for all elements at once using
Space::set_uniform_order(). Note that element degrees
are stored in Space, not in Mesh. The reason is that in Hermes one can
have multiple spaces with different element degrees and even types 
over the same mesh. In Hermes, Mesh only stores geometrical information.
A space created in this way is ready for use. 

Visualizing basis functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

As a debugging/learning feature, Hermes can visualize the basis of each Space.
Similarly to MeshView, one can create a BaseView object and use it 
to display the entire basis (VectorBaseView has to be used for vector-valued 
approximations in spaces Hcurl and Hdiv - this will be discussed later). 
One can cycle through all basis functions in the window using the arrow keys. 
If you press the left mouse button at the beginning, you will see the Dirichlet 
lift (a function that represents Dirichlet boundary conditions).

3D view
~~~~~~~

This is how the last figure above was obtained (press the '3' key for 3D mode).
We suggest that you spend some time experimenting with element refinements and 
hanging nodes to see how basis functions on irregular meshes look like.
