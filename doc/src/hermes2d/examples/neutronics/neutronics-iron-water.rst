Iron-Water
----------

**Git reference:** Example `neutronics-iron-water <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/examples/neutronics/neutronics-iron-water>`_.

This example is very similar to the example "saphir", the main difference being that 
it reads a mesh file in the exodusii format (created by Cubit). This example only builds 
if you have the `ExodusII <http://sourceforge.net/projects/exodusii/>`_ and 
`NetCDF <http://www.unidata.ucar.edu/software/netcdf/>`_ libraries installed on your system and 
the variables WITH_EXODUSII, EXODUSII_ROOT and NETCDF_ROOT defined properly. 
The latter can be done, for example, in the CMake.vars file as follows:

::

    SET(WITH_EXODUSII YES)
    SET(EXODUSII_ROOT /opt/packages/exodusii)
    SET(NETCDF_ROOT   /opt/packages/netcdf)

The mesh is now loaded using the ExodusIIReader (see 
the `mesh_loader.h <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/src/mesh/mesh_loader.h>`_ file):

::

    // Load the mesh
    Mesh mesh;
    ExodusIIReader mloader;
    if (!mloader.load("iron-water.e", &mesh)) error("ExodusII mesh load failed.");

The model describes an external-force-driven configuration without fissile materials present.
We will solve the one-group neutron diffusion equation

.. math::
    :label: iron-water

       -\nabla \cdot (D(x,y) \nabla \Phi) + \Sigma_a(x,y) \Phi = Q_{ext}(x,y).

The domain of interest is a 30 x 30 cm square consisting of four regions.
A uniform volumetric source is placed in water in the lower-left corner 
of the domain, surrounded with a layer of water, a layer of iron, and finally
another layer of water:

.. image:: example-iron-water/iron-water.png
   :align: center
   :width: 400
   :height: 400
   :alt: Schematic picture for the iron-water example.

The unknown is the neutron flux $\Phi(x, y)$. The values of the diffusion coefficient 
$D(x, y)$, absorption cross-section $\Sigma_a(x, y)$ and the source term $Q_{ext}(x,y)$
are constant in the subdomains. The source $Q_{ext} = 1$ in area 1 and zero 
elsewhere. The boundary conditions for this problem are zero Dirichlet (right and top edges)
and zero Neumann (bottom and left edges). Sample results of this computation are shown below.

Solution:

.. image:: example-iron-water/iron-water-sol.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution to the iron-water example.


Final mesh (h-FEM with linear elements):

.. image:: example-iron-water/iron-water-mesh-h1.png
   :align: center
   :width: 440
   :height: 400
   :alt: Final finite element mesh for the iron-water example (h-FEM with linear elements).

Final mesh (h-FEM with quadratic elements):

.. image:: example-iron-water/iron-water-mesh-h2.png
   :align: center
   :width: 440
   :height: 400
   :alt: Final finite element mesh for the iron-water example (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: example-iron-water/iron-water-mesh-hp.png
   :align: center
   :width: 440
   :height: 400
   :alt: Final finite element mesh for the iron-water example (hp-FEM).

DOF convergence graphs:

.. image:: example-iron-water/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph for example iron-water.

CPU time convergence graphs:

.. image:: example-iron-water/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph for example iron-water.

