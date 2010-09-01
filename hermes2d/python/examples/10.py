#! /usr/bin/env python

# This example shows how to run adaptive hp-FEM, h-FEM and p-FEM with
# basic control parameters. The underlying problem is a planar model
# of an electrostatic micromotor (MEMS). You may want to first run the model
# with hp-FEM (CAND_LIST = CandList.HP_ANISO), then with h-FEM (CAND_LIST = CandList.H_ANISO),
# and last with p-FEM (CAND_LIST = CandList.P_ANISO). The last option
# is there solely for completeness, since adaptive p-FEM is not really
# useful in practice. Uniform initial polynomial degree of mesh elements
# can be set using the variable P_INIT. 
#
# The function adapt(...) takes the parameters selector, THRESHOLD,
# STRATEGY and MESH_REGULARITY whose meaning is explained below.
# Only the first two parameters selector and THRESHOLD are mandatory, all
# others can be left out and in that case default values are used.
#
# Additional control parameters are possible, these are demonstrated
# in the next tutorial example. Two types of convergence graphs are created
# -- error estimate wrt. the number of degrees of freedom (DOF), and error
# estimate wrt. CPU time. Later you will learn that also the error wrt.
# exact solution can be created for problems where exact solution is
# available. In this example you also can see how
# to define different material parameters in various parts of the
# computational domain.
#
# PDE: -div[eps_r(x,y) grad phi] = 0
#      eps_r = EPS1 in Omega_1 (surrounding air)
#      eps_r = EPS2 in Omega_2 (moving part of the motor)
#
# BC: phi = 0 V on Gamma_1 (left edge and also the rest of the outer boundary
#     phi = VOLTAGE on Gamma_2 (boundary of stator)

# Import modules
from hermes2d import Mesh, MeshView, VectorView, OrderView, H1Shapeset, PrecalcShapeset, H1Space, \
        WeakForm, Solution, ScalarView, LinSystem, DummySolver, RefSystem, \
    H1Adapt, H1ProjBasedSelector, CandList, \
        H2D_EPS_HIGH, H2D_FN_DX, H2D_FN_DY

from hermes2d.examples.c10 import set_bc, set_forms
from hermes2d.examples import get_motor_mesh

# The following parameters can be changed:
SOLVE_ON_COARSE_MESH = True   # If true, coarse mesh FE problem is solved in every adaptivity step
P_INIT = 2              # Initial polynomial degree of all mesh elements.
THRESHOLD = 0.2         # This is a quantitative parameter of the adapt(...) function and
                        # it has different meanings for various adaptive strategies (see below).

STRATEGY = 1            # Adaptive strategy:
                        # STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                        #   error is processed. If more elements have similar errors, refine
                        #   all to keep the mesh symmetric.
                        # STRATEGY = 1 ... refine all elements whose error is larger
                        #   than THRESHOLD times maximum element error.
                        # STRATEGY = 2 ... refine all elements whose error is larger
                        #   than THRESHOLD.
                        # More adaptive strategies can be created in adapt_ortho_h1.cpp.

CAND_LIST = CandList.H2D_HP_ANISO_H  # Predefined list of element refinement candidates.
                        # Possible values are are attributes of the class CandList:
                        # H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO
                        # See User Documentation for details.                      

MESH_REGULARITY = -1    # Maximum allowed level of hanging nodes:
                        # MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                        # MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                        # MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                        # Note that regular meshes are not supported, this is due to
                        # their notoriously bad performance.

ERR_STOP = 1.0          # Stopping criterion for adaptivity (rel. error tolerance between the
                        # fine mesh and coarse mesh solution in percent).
CONV_EXP = 1.0;         # Default value is 1.0. This parameter influences the selection of
                        # cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                        # fine mesh and coarse mesh solution in percent).
NDOF_STOP = 60000       # Adaptivity process stops when the number of degrees of freedom grows
                        # over this limit. This is to prevent h-adaptivity to go on forever.

H2DRS_DEFAULT_ORDER = -1 # A default order. Used to indicate an unkonwn order or a maximum support order

# Load the mesh
mesh = Mesh()
mesh.load(get_motor_mesh())

# Create an H1 space with default shapeset
space = H1Space(mesh, P_INIT)
set_bc(space)

# Initialize the discrete problem
wf = WeakForm()
set_forms(wf)

# Initialize views
sview = ScalarView("Coarse solution", 0, 0, 600, 1000)
#gview = VectorView("Gradient", 610, 0, 600, 1000)
oview = OrderView("Polynomial orders", 1220, 0, 600, 1000)

# Initialize refinement selector.
selector = H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER)

# Initialize the linear system.
ls = LinSystem(wf)
ls.set_spaces(space)

# Adaptivity loop:
sln_coarse = Solution()
sln_fine = Solution()
done = False
it = 1

while(not done):

    print("\n---- Adaptivity step %d ---------------------------------------------\n" % it)
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
        
    # View the coarse mesh solution
    sview.show(sln_coarse)

    # View the mesh
    mesh.plot()

    # Calculate element errors and total error estimate
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

#print ("Total running time: %d sec" % cpu)

# Show the fine solution - this is the final result
sview.show(sln_fine)
#gview.show(sln_fine, sln_fine, H2D_EPS_HIGH, H2D_FN_DX_0, H2D_FN_DY_0)

