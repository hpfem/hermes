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

// $Id$

#ifndef __HERMES_2D_H
#define __HERMES_2D_H

// hermes_common solvers
#include "../hermes_common/solver/amesos.h"
#include "../hermes_common/solver/aztecoo.h"
#include "../hermes_common/solver/epetra.h"
#include "../hermes_common/solver/mumps.h"
#include "../hermes_common/solver/nox.h"
#include "../hermes_common/solver/pardiso.h"
#include "../hermes_common/solver/petsc.h"
#include "../hermes_common/solver/umfpack_solver.h"
#include "../hermes_common/solver/superlu.h"

// preconditioners
#include "../hermes_common/solver/precond.h"
#include "../hermes_common/solver/precond_ifpack.h"
#include "../hermes_common/solver/precond_ml.h"

// boundary conditions
#include "../hermes_common/bctypes.h"

#include "h2d_common.h"
#include "hermes_logging.h"

#include "range.h"
#include "limit_order.h"

#include "mesh.h"
#include "mesh_loader.h"
#include "h2d_reader.h"
#include "exodusii.h"

#include "space/space_h1.h"
#include "space/space_hcurl.h"
#include "space/space_l2.h"
#include "space/space_hdiv.h"

#include "quad_all.h"
#include "shapeset/shapeset_h1_all.h"
#include "shapeset/shapeset_hc_all.h"
#include "shapeset/shapeset_hd_all.h"
#include "shapeset/shapeset_l2_all.h"

#include "refmap.h"
#include "traverse.h"
#include "trans.h"

#include "weakform.h"
#include "discrete_problem.h"
#include "forms.h"

#include "integrals_h1.h"
#include "integrals_hcurl.h"
#include "integrals_hdiv.h"

#include "solution.h"
#include "filter.h"

#include "norm.h"
#include "graph.h"

#include "views/view.h"
#include "views/base_view.h"
#include "views/mesh_view.h"
#include "views/order_view.h"
#include "views/scalar_view.h"
#include "views/stream_view.h"
#include "views/vector_base_view.h"
#include "views/vector_view.h"

#include "refinement_type.h"
#include "element_to_refine.h"
#include "ref_selectors/selector.h"
#include "ref_selectors/order_permutator.h"
#include "ref_selectors/optimum_selector.h"
#include "ref_selectors/proj_based_selector.h"
#include "ref_selectors/l2_proj_based_selector.h"
#include "ref_selectors/h1_proj_based_selector.h"
#include "ref_selectors/hcurl_proj_based_selector.h"

#include "adapt/adapt.h"
#include "neighbor.h"
#include "ogprojection.h"

#include "numerical_flux.h"
#include "runge_kutta.h"
#include "tables.h"
/**

\mainpage

This manual documents the source code of hermes2d. It is intended for the developers of
the library. If you are only interested in using hermes2d, please refer to the User's Manual.

This page provides an overview of the structure of the project, in a bottom-up fashion.
You can also visit the <a href="annotated.html">Class List</a> and the
<a href="hierarchy.html">Class Hierarchy</a> pages.


<br>
<h2>Mesh Classes</h2>

A finite element mesh consists of elements and nodes, defined by the structures Element and Node.
Elements can be triangles or quadrilaterals (quads). Each element contains pointers to three or
four vertex nodes and to three or four edge nodes. Vertex nodes store the physical coordinates
of mesh vertices and edge nodes provide element connectivity across mesh edges. Our meshes can
be irregular, i.e., we allow hanging vertex nodes (sometimes called T-junctions). Elements can be
active or inactive. Inactive elements are those which have been refined and are no longer part
of the actual mesh. Their purpose is to provide pointers to the refined (son) elements.

There are two main classes that define the finite element mesh as a whole, HashTable and Mesh.
The first serves as a container for nodes, plus provides node search functions via hashing, the
latter contains an Array of elements and functions for loading the mesh from a file and functions
for refining elements.

<img src="classMesh.png">

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
vertex node). A detailed description of these mechanisms can be found in <a href="../data/chap5.pdf">
Chapter 5</a> of the thesis of Jakub Cerveny.

The code supports elements with curved edges defined by NURBS curves. If an element is
curvilinear, its member "cm" points to the structure CurvMap, which in turn defines the curvature
of the edges using the structures Nurbs. Since NURBS curves are non-polynomial, the curvilinear
reference mapping has to be projected to a polynomial one (see curved.cpp). Once this is done, the
curvilinear mapping is defined by adding edge and bubble shape functions to the standard affine
reference mapping.

Throughout the code, the macros for_all_active_elements(), for_all_vertex_nodes(), etc. (see mesh.h),
are used to iterate through elements and nodes of the mesh.

A commented example .mesh file can be found <a href="../data/example.mesh">here</a>.

Relevant files: array.h, hash.h, mesh.h, curved.h


<br>
<h2>Shapeset Classes</h2>

In higher-order FEM, shapeset is a collection of the so-called shape functions, which are polynomials
that can be used to construct higher-order basis functions. There are three types of shape functions:
vertex, edge and bubble functions. Vertex functions are standard linear functions for the construction
of the well-known pyramid basis functions. By gluing together two edge shape functions, one obtains
edge basis functions, which generate higher-order solutions along edges. Finally, bubble functions
complete the higher-order approximation on individual elements. For more information, refer to
<a href="../data/chap3.pdf">Chapter 3</a> of Jakub's thesis. Hermes2D uses hierarchic functions
exclusively.

<img src="classShapeset.png">

Shapeset is a base class for all shapesets, both scalar (H1) and vector-valued (H(curl)). It defines
the interface (basically it is just a bunch of tables) and the functionality for constrained shape
functions.

The actual shapesets are H1ShapesetOrtho, H1ShapesetJacobi, HcurlShapesetLegendre and
HcurlShapesetGradLeg. The other ("eigen") shapesets are experimental and not recommended for normal
use.

Relevant files: shapeset.h, shapeset_h1_all.h, shapeset_hc_all.h



<br>
<h2>Space Classes</h2>

The Space class corresponds to the mathematical notion of a finite element space over a 2D domain.
Given a mesh and a shapeset, the Space class constructs a specific higher-order FE space. Currently,
three space types are defined: H1Space, HcurlSpace and L2Space. These represent their respective
mathematical counterparts.

<img src="classSpace.png">

The main purpose of these classes is (1) to hold polynomial degrees of mesh elements, (2) to construct
and enumerate the basis functions of the space and (3) to return assembly lists (see AsmList),
which are used in the stiffness matrix assembly.

Multiple spaces can share the same mesh. This is used when solving systems of PDEs. In multi-mesh
computations, each space usually has its own mesh, but all meshes must descend from a single
master mesh.

All of the three spaces are able to construct higher-order basis functions on irregular meshes.
For details, see the documentation for the Space class.

Relevant files: space.h, space_h1.h, space_hcurl.h, space_l2.h, asmlist.h



<br>
<h2>Quadrature Classes</h2>

Quad1D, Quad2D, edge tables (24)

Relevant files: quad.h, quad_all.h



<br>
<h2>PrecalcShapeset</h2>

TODO

<img src="classPrecalcShapeset.png">

Relevant files: transform.h, function.h, precalc.h



<br>
<h2>RefMap Class</h2>

TODO

<img src="classRefMap.png">

Relevant files: refmap.h





<br>
<h2>Multi-mesh Assembling</h2>

TODO

Relevant file: traverse.h



<br>
<h2>Solution and Filter Classes</h2>

TODO

<img src="classMeshFunction.png">

Relevant files: solution.h, filter.h



<br>
<h2>Linearizer Classes</h2>

TODO

<img src="classLinearizer.png">

Relevant files: linear.h



<br>
<h2>View Classes</h2>

TODO

<img src="classView.png">

Relevant files: view.h




*/




#endif
