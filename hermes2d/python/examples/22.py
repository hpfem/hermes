#! /usr/bin/env python

#  This is another example that allows you to compare h- and hp-adaptivity from the point of view
#  of both CPU time requirements and discrete problem size, look at the quality of the a-posteriori
#  error estimator used by Hermes (exact error is provided), etc. You can also change
#  the parameter MESH_REGULARITY to see the influence of hanging nodes on the adaptive process.
#  The problem is made harder for adaptive algorithms by increasing the parameter SLOPE.
#
#  PDE: -Laplace u = f
#
#  Known exact solution, see functions fn() and fndd()
#
#  Domain: unit square (0, 1)x(0, 1), see the file square.mesh
#
#  BC:  Dirichlet, given by exact solution

# Import modules
from hermes2d import (Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space,
        WeakForm, Solution, DummySolver, LinSystem, ScalarView, RefSystem,
        H1Adapt, H1ProjBasedSelector, CandList, \
    set_verbose)
from hermes2d.examples.c22 import set_bc, set_forms

#  The following parameters can be changed:
SOLVE_ON_COARSE_MESH = True   # if true, coarse mesh FE problem is solved in every adaptivity step
INIT_REF_NUM = 1               # Number of initial uniform mesh refinements
P_INIT = 2              # Initial polynomial degree of all mesh elements.
THRESHOLD = 0.3         # This is a quantitative parameter of the adapt(...) function and
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
CONV_EXP = 0.5
ERR_STOP = 0.1         # Stopping criterion for adaptivity (rel. error tolerance between the
                            # fine mesh and coarse mesh solution in percent).
NDOF_STOP = 60000       # Adaptivity process stops when the number of degrees of freedom grows
                            # over this limit. This is to prevent h-adaptivity to go on forever.

H2DRS_DEFAULT_ORDER = -1 # A default order. Used to indicate an unkonwn order or a maximum support order

# Problem parameters.
SLOPE = 60           # Slope of the layer.

# Load the mesh
mesh = Mesh()
mesh.create([
        [0, 0],
        [1, 0],
        [1, 1],
        [0, 1],
    ], [
        [2, 3, 0, 1, 0],
    ], [
        [0, 1, 1],
        [1, 2, 1],
        [2, 3, 1],
        [3, 0, 1],
    ], [])

# Perform initial mesh refinements
mesh.refine_all_elements()

# Create an H1 space with default shapeset
space = H1Space(mesh, P_INIT)
set_bc(space)

# Initialize the weak formulation
wf = WeakForm()
set_forms(wf)

# Initialize views
sview = ScalarView("Solution")
mview = MeshView("Mesh")
graph = []

# Initialize refinement selector
selector = H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER)

# Initialize the coarse mesh problem
ls = LinSystem(wf)
ls.set_spaces(space)

# Adaptivity loop
iter = 0
done =  False
print "Calculating..."
sln_coarse = Solution()
sln_fine = Solution()

while (not done):
    
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

# Show the fine solution - the final result
sview.show(sln_fine)
