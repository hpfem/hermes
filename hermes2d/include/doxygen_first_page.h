// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D. If not, see <http://www.gnu.org/licenses/>.

/**

\mainpage

This manual documents the source code of hermes2d. It is intended for the developers of
the library. If you are only interested in using hermes2d, please refer to the User's Manual.

This page provides an overview of the structure of the project, in a bottom-up fashion.
You can also visit the <a href="annotated.html">Class List</a> and the
<a href="hierarchy.html">Class Hierarchy</a> pages.

One is also advised to refer to
<h1><a href="modules.html">Modules</a></h1>


<br>
<h2>Mesh Classes</h2>

A finite element mesh consists of elements and nodes, defined by the structures 
<a href="classHermes_1_1Hermes2D_1_1Element.html">Element</a> and <a href="structHermes_1_1Hermes2D_1_1Node.html">Node</a>.
Elements can be triangles or quadrilaterals (quads). Each element contains pointers to three or
four vertex nodes and to three or four edge nodes. Vertex nodes store the physical coordinates
of mesh vertices and edge nodes provide element connectivity across mesh edges. Our meshes can
be irregular, i.e., we allow hanging vertex nodes (sometimes called T-junctions). Elements can be
active or inactive. Inactive elements are those which have been refined and are no longer part
of the actual mesh. Their purpose is to provide pointers to the refined (son) elements.

There are two main classes that define the finite element mesh as a whole, 
<a href="classHermes_1_1Hermes2D_1_1HashTable.html">HashTable</a> and <a href="classHermes_1_1Hermes2D_1_1Mesh.html">Mesh</a>.
The first serves as a container for nodes, plus provides node search functions via hashing, the
latter contains an Array of elements and various methods.

<img src="classHermes_1_1Hermes2D_1_1Mesh.png">

Except for top-level vertex nodes, all nodes can be reached by providing the id numbers of their
"parents". For example, an edge node exists halfway between two vertex nodes, and is created by
calling HashTable::get_edge_node() and passing the id's of the two vertex nodes. If the edge node
does not exist, it is created first. This greatly simplifies mesh initialization and hanging nodes
in general, since one does not have to care whether the edge node already exists and from which
element to reach it. The class HashTable achieves this by hashing the id numbers of the parents.
The variable "ref" is maintained in each node, specifying the number of elements using the node.
This determines when a node should be deleted (ref = 0), but can also be used to tell the type
of the node. Standard edge nodes have ref = 2, hanging or boundary edge nodes have ref = 1. Similarly
for vertex nodes, where ref can be 2 or 3 (hanging), 4 or 6 (regular) or a large number (top-level
vertex node). A detailed description of these mechanisms can be found in Chapter 5 of the thesis of Jakub Cerveny.

The code supports elements with curved edges defined by NURBS curves. If an element is
curvilinear, its member "cm" points to the structure CurvMap, which in turn defines the curvature
of the edges using the structures Nurbs. Since NURBS curves are non-polynomial, the curvilinear
reference mapping has to be projected to a polynomial one (see <a href="curved_8h_source.html">curved.h</a>). Once this is done, the
curvilinear mapping is defined by adding edge and bubble shape functions to the standard affine
reference mapping.

Throughout the code, the macros for_all_active_elements(), for_all_vertex_nodes(), etc. (see <a href="mesh_8h_source.html">mesh.h</a>),
are used to iterate through elements and nodes of the mesh.

Relevant files are located in the <b>mesh</b> directory (see <a href="dirs.html">directory list</a>)


<br>
<h2>Shapeset Classes</h2>

In higher-order FEM, shapeset is a collection of the so-called shape functions, which are polynomials
that can be used to construct higher-order basis functions. There are three types of shape functions:
vertex, edge and bubble functions. Vertex functions are standard linear functions for the construction
of the well-known pyramid basis functions. By gluing together two edge shape functions, one obtains
edge basis functions, which generate higher-order solutions along edges. Finally, bubble functions
complete the higher-order approximation on individual elements. For more information, refer to
Chapter 3 of Jakub's thesis. Hermes2D uses hierarchic functions
exclusively.

<img src="classHermes_1_1Hermes2D_1_1Shapeset.png">

Shapeset is a base class for all shapesets, both Scalar (H1) and vector-valued (H(curl)). It defines
the interface (basically it is just a bunch of tables) and the functionality for constrained shape
functions.

The actual shapesets are H1ShapesetOrtho, H1ShapesetJacobi, HcurlShapesetLegendre and
HcurlShapesetGradLeg. The other ("eigen") shapesets are experimental and not recommended for normal
use.

Relevant files are located in the <b>shapeset</b> directory (see <a href="dirs.html">directory list</a>)

<br>
<h2>Space Classes</h2>

The Space class corresponds to the mathematical notion of a finite element space over a 2D domain.
Given a mesh and a shapeset, the Space class constructs a specific higher-order FE space. Currently,
four space types are defined: H1Space, HcurlSpace, HdivSpace and L2Space. These represent their respective
mathematical counterparts.

<img src="classHermes_1_1Hermes2D_1_1Space.png">

The main purpose of these classes is (1) to hold polynomial degrees of mesh elements, (2) to construct
and enumerate the basis functions of the space and (3) to return assembly lists (see AsmList<Scalar>),
which are used in the stiffness matrix assembly.

Multiple spaces can share the same mesh. This is used when solving systems of PDEs. In multi-mesh
computations, each space usually has its own mesh, but all meshes must descend from a single
master mesh.

All of the three spaces are able to construct higher-order basis functions on irregular meshes.
For details, see the documentation for the Space class.

Relevant files are located in the <b>space</b> directory (see <a href="dirs.html">directory list</a>)

<br>
<h2>Quadrature Classes</h2>

Classes holding quadrature formulas.
<center>
  <table cellpadding=10 border=0>
    <tr>
      <td><img src="classHermes_1_1Hermes2D_1_1Quad1D.png"></td>
      <td><img src="classHermes_1_1Hermes2D_1_1Quad2D.png"></td>
    </tr>
  </table>
</center>
Relevant files are located in the <b>quadrature</b> directory (see <a href="dirs.html">directory list</a>)

<br>
<h2>Functions in Hermes</h2>

Types of functions representations in Hermes.

<h3>MeshFunction class and its descendants</h3>

<img src="classHermes_1_1Hermes2D_1_1Function.png" width=100%>

Relevant files are located in the <b>function</b> directory, except for hermes_function.h (see <a href="dirs.html">directory list</a>)

<h3>HermesFunction</h3>

<a href="classHermes_1_1Hermes2D_1_1HermesFunction.html">HermesFunction</a>

Relevant file is hermes_function.h in <b>function</b> directory (see <a href="dirs.html">directory list</a>)

<br>
<h2>RefMap Class</h2>

Describes reference mapping of an actual physical Element onto the reference one.

<img src="classHermes_1_1Hermes2D_1_1RefMap.png">

Relevant file is refmap.h in <b>mesh</b> directory (see <a href="dirs.html">directory list</a>)

<br>
<h2>Multi-mesh Assembling</h2>

Class facilitates monolithic multi-mesh traversing.

<a href="classHermes_1_1Hermes2D_1_1Traverse">Traverse</a>

Relevant file is traverse.h in <b>mesh</b> directory (see <a href="dirs.html">directory list</a>)

<br>
<h2>Visualization Classes</h2>

<h3>View classes</h3>

<img src="classHermes_1_1Hermes2D_1_1Views_1_1View.png" width=90%>

Relevant files are located in the <b>views</b> directory (see <a href="dirs.html">directory list</a>)

<h3>Linearizer classes</h3>

<img src="classHermes_1_1Hermes2D_1_1Views_1_1LinearizerBase.png">

Relevant files are located in the <b>views</b> directory (see <a href="dirs.html">directory list</a>)
*/
