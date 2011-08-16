#! /usr/bin/env python

# This example explains how to use the multimesh adaptive hp-FEM,
# where different physical fields (or solution components) can be
# approximated using different meshes and equipped with mutually
# independent adaptivity mechanisms. Here we consider linear elasticity
# and will approximate each displacement components using an individual
# mesh.
#
# PDE: Lame equations of linear elasticity, treated as a coupled system
#      of two PDEs
#
# BC: u_1 = u_2 = 0 on Gamma_1
#     du_2/dn = f on Gamma_2
#     du_1/dn = du_2/dn = 0 elsewhere

# Import modules
from hermes2d import Mesh, MeshView, VectorView, OrderView, H1Shapeset, PrecalcShapeset, H1Space, \
        WeakForm, Solution, ScalarView, LinSystem, DummySolver, RefSystem, \
        H1Adapt, H1ProjBasedSelector, CandList, \
        VonMisesFilter

from hermes2d.examples.c11 import set_bc, set_wf_forms, set_hp_forms
from hermes2d.examples import get_bracket_mesh

# The following parameters can be changed: In particular, compare hp- and
# h-adaptivity via the CAND_LIST option, and compare the multi-mesh vs. single-mesh
# using the MULTI parameter.

SOLVE_ON_COARSE_MESH = True  # If true, coarse mesh FE problem is solved in every adaptivity step.
P_INIT_U = 2             # Initial polynomial degree for u
P_INIT_V = 2             # Initial polynomial degree for v
INIT_REF_BDY = 3         # Number of initial boundary refinements
MULTI = True             # MULTI = true  ... use multi-mesh,
                            # MULTI = false ... use single-mesh.
                            # Note: In the single mesh option, the meshes are
                            # forced to be geometrically the same but the
                            # polynomial degrees can still vary.
THRESHOLD = 0.3          # This is a quantitative parameter of the adapt(...) function and
                                 # it has different meanings for various adaptive strategies (see below).
STRATEGY = 1             # Adaptive strategy:
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

MESH_REGULARITY = -1     # Maximum allowed level of hanging nodes:
                            # MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                            # MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                            # MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                            # Note that regular meshes are not supported, this is due to
                            # their notoriously bad performance.
CONV_EXP = 1             # Default value is 1.0. This parameter influences the selection of
                            # cancidates in hp-adaptivity. See get_optimal_refinement() for details.
MAX_ORDER = 10           # Maximum allowed element degree
ERR_STOP = 0.5           # Stopping criterion for adaptivity (rel. error tolerance between the
                            # fine mesh and coarse mesh solution in percent).
NDOF_STOP = 60000        # Adaptivity process stops when the number of degrees of freedom grows over
                            # this limit. This is mainly to prevent h-adaptivity to go on forever.

H2DRS_DEFAULT_ORDER = -1 # A default order. Used to indicate an unkonwn order or a maximum support order

# Load the mesh
umesh = Mesh()
vmesh = Mesh()
umesh.load(get_bracket_mesh())
if MULTI == False:
    umesh.refine_towards_boundary(1, INIT_REF_BDY)
    
# Create initial mesh (master mesh).
vmesh.copy(umesh)

# Initial mesh refinements in the vmesh towards the boundary
if MULTI == True:
    vmesh.refine_towards_boundary(1, INIT_REF_BDY)

# Create the x displacement space
uspace = H1Space(umesh, P_INIT_U)
vspace = H1Space(vmesh, P_INIT_V)

# Initialize the weak formulation
wf = WeakForm(2)
set_wf_forms(wf)

# Initialize views
uoview = OrderView("Coarse mesh for u", 0, 0, 360, 300)
voview = OrderView("Coarse mesh for v", 370, 0, 360, 300)
uview = ScalarView("Coarse mesh solution u", 740, 0, 400, 300)
vview = ScalarView("Coarse mesh solution v", 1150, 0, 400, 300)

# Initialize refinement selector
selector = H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER)

# Initialize the coarse mesh problem
ls = LinSystem(wf)
ls.set_spaces(uspace, vspace)

# adaptivity loop
it = 1
done = False
u_sln_coarse = Solution()
v_sln_coarse = Solution()
u_sln_fine = Solution()
v_sln_fine = Solution()

while(not done):

    print ("\n---- Adaptivity step %d ---------------------------------------------\n" % it)
    it += 1
    
    # Assemble and Solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(u_sln_fine, v_sln_fine, lib="scipy")

    # Either solve on coarse mesh or project the fine mesh solution 
    # on the coarse mesh.
    if SOLVE_ON_COARSE_MESH:
        ls.assemble()
        ls.solve_system(u_sln_coarse, v_sln_coarse, lib="scipy")
    else:
        ls.project_global()

    # View the solution and meshes
    uview.show(u_sln_coarse)
    vview.show(v_sln_coarse)
    umesh.plot(space=uspace)
    vmesh.plot(space=vspace)

    # Calculate element errors and total error estimate
    hp = H1Adapt(ls)
    hp.set_solutions([u_sln_coarse, v_sln_coarse], [u_sln_fine, v_sln_fine]);
    set_hp_forms(hp)
    err_est = hp.calc_error() * 100

    print("Error estimate: %s" % err_est)

    # If err_est too large, adapt the mesh
    if err_est < ERR_STOP:
        done = True
    else:
        MULTI = False if MULTI == True else True
        hp.adapt(selector, THRESHOLD, STRATEGY, MESH_REGULARITY, MULTI) 
        if ls.get_num_dofs() >= NDOF_STOP:
            done = True
    

# Show the fine solution - this is the final result
uview.show(u_sln_fine)
vview.show(v_sln_fine)
