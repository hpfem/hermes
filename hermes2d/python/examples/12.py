#! /usr/bin/env python

#  This example uses automatic adaptivity to solve a general second-order linear
#  equation with non-constant coefficients.
#
#  PDE: -d/dx(a_11(x,y)du/dx) - d/dx(a_12(x,y)du/dy) - d/dy(a_21(x,y)du/dx) - d/dy(a_22(x,y)du/dy)
#       + a_1(x,y)du/dx + a_21(x,y)du/dy + a_0(x,y)u = rhs(x,y)
#
#  Domain: arbitrary
#
#  BC:  Dirichlet for boundary marker 1: u = g_D(x,y)
#       Natural for any other boundary marker:   (a_11(x,y)*nu_1 + a_21(x,y)*nu_2) * dudx
#                                              + (a_12(x,y)*nu_1 + s_22(x,y)*nu_2) * dudy = g_N(x,y)

# Import modules
from hermes2d import Mesh, MeshView, VectorView, OrderView, H1Shapeset, PrecalcShapeset, H1Space, \
WeakForm, Solution, ScalarView, LinSystem, DummySolver, RefSystem, \
H1Adapt, H1ProjBasedSelector, CandList

from hermes2d.examples.c12 import set_bc, set_forms
from hermes2d.examples import get_12_mesh

#  The following parameters can be changed:
SOLVE_ON_COARSE_MESH = True   # if true, coarse mesh FE problem is solved in every adaptivity step
P_INIT = 2              # Initial polynomial degree of all mesh elements.
THRESHOLD = 0.6         # This is a quantitative parameter of the adapt(...) function and
                        # it has different meanings for various adaptive strategies (see below).
STRATEGY = 0            # Adaptive strategy:
                            # STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                            #   error is processed. If more elements have similar errors, refine
                            #   all to keep the mesh symmetric.
                            # STRATEGY = 1 ... refine all elements whose error is larger
                            #   than THRESHOLD times maximum element error.
                            # STRATEGY = 2 ... refine all elements whose error is larger
                            #   than THRESHOLD.
                            # More adaptive strategies can be created in adapt_ortho_h1.cpp.
CAND_LIST = CandList.H2D_HP_ANISO  # Predefined list of element refinement candidates.
                        # Possible values are are attributes of the class CandList:
                        # P_ISO, P_ANISO, H_ISO, H_ANISO, HP_ISO, HP_ANISO_H, HP_ANISO_P, HP_ANISO
                        # See the Sphinx tutorial (http://hpfem.org/hermes2d/doc/src/tutorial-2.html#adaptive-h-fem-and-hp-fem) for details.
MESH_REGULARITY = -1    # Maximum allowed level of hanging nodes:
                            # MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                            # MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                            # MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                            # Note that regular meshes are not supported, this is due to
                            # their notoriously bad performance.
CONV_EXP = 1.0
ERR_STOP = 0.1         # Stopping criterion for adaptivity (rel. error tolerance between the
                            # fine mesh and coarse mesh solution in percent).
NDOF_STOP = 60000       # Adaptivity process stops when the number of degrees of freedom grows
                            # over this limit. This is to prevent h-adaptivity to go on forever.

H2DRS_DEFAULT_ORDER = -1 # A default order. Used to indicate an unkonwn order or a maximum support order

# Boundary markers
BDY_DIRICHLET = 1
BDY_NEUMANN = 2

# Load the mesh
mesh = Mesh()
mesh.load(get_12_mesh())

# Perform initial mesh refinements
mesh.refine_all_elements()

# Create an H1 space with default shapeset
space = H1Space(mesh, P_INIT)
set_bc(space)

# Initialize the weak formulation
wf = WeakForm()
set_forms(wf)

# Initialize views
sview = ScalarView("Coarse solution", 0, 0, 600, 1000)
oview = OrderView("Polynomial orders", 1220, 0, 600, 1000)

# Initialize refinement selector
selector = H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER)

# Initialize the linear system.
ls = LinSystem(wf)
ls.set_spaces(space)

# Adaptivity loop
it = 0
done = False
sln_coarse = Solution()
sln_fine = Solution()

while (not done):
    print("\n---- Adaptivity step %d ---------------------------------------------\n" % (it+1))
    it += 1

    # Assemble and solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(sln_fine)
    
    # Either solve on coarse mesh or project the fine mesh solution 
    # on the coarse mesh.   
    if SOLVE_ON_COARSE_MESH:
        ls.assemble()
        ls.solve_system(sln_coarse)
    else:
        ls.project_global(sln_fine, sln_coarse)
        
    # View the solution and mesh
    sview.show(sln_coarse);
    mesh.plot(space=space)

    # Calculate error estimate wrt. fine mesh solution
    hp = H1Adapt(ls)
    hp.set_solutions([sln_coarse], [sln_fine])
    err_est = hp.calc_error() * 100
    print("Error estimate: %d" % err_est)

    # If err_est too large, adapt the mesh
    if (err_est < ERR_STOP):
        done = True
    else:
        done = hp.adapt(selector, THRESHOLD, STRATEGY, MESH_REGULARITY)
        if (ls.get_num_dofs() >= NDOF_STOP):
            done = True

sview.show(sln_fine)
