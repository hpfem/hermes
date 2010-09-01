from hermes2d import Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, WeakForm, Solution, ScalarView, VonMisesFilter, \
        set_verbose, set_warn_integration, DummySolver, H2D_EPS_HIGH, H2D_FN_DX, \
        H2D_FN_DY, \
    H1Adapt, H1ProjBasedSelector, CandList, \
    RefSystem
from hermes2d.forms import set_forms
from hermes2d.examples import get_example_mesh, get_sample_mesh, \
        get_cylinder_mesh, get_07_mesh, get_cathedral_mesh, get_bracket_mesh

domain_mesh = get_example_mesh()
sample_mesh = get_sample_mesh()
cylinder_mesh = get_cylinder_mesh()

def test_example_01():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    mesh.refine_all_elements()

def test_example_02():
    set_verbose(False)    
    P_INIT = 3

    # Load the mesh file
    domain_mesh = get_example_mesh()    # Original L-shape domain
    mesh = Mesh()
    mesh.load(domain_mesh)

    # Refine all elements (optional)
    mesh.refine_all_elements()

    # Create a shapeset and an H1 space
    space = H1Space(mesh)

    # Assign element orders and initialize the space
    space.set_uniform_order(P_INIT)    # Set uniform polynomial order
                                       # P_INIT to all mesh elements

def test_example_03():
    from hermes2d.examples.c03 import set_bc

    set_verbose(False)

    P_INIT = 5                # Uniform polynomial degree of mesh elements.

    # Problem parameters.
    CONST_F = 2.0

    # Load the mesh file
    mesh = Mesh()
    mesh.load(get_example_mesh())

    # Sample "manual" mesh refinement
    mesh.refine_all_elements()

    # Create an H1 space with default shapeset
    space = H1Space(mesh, P_INIT)
    set_bc(space)

    # Initialize the weak formulation
    wf = WeakForm(1)
    set_forms(wf)

    # Initialize the linear system
    ls = LinSystem(wf)
    ls.set_spaces(space)

    # Assemble and solve the matrix problem.
    sln = Solution()
    ls.assemble()
    ls.solve_system(sln)

    #assert abs(sln.l2_norm() - 0.25493) < 1e-4
    #assert abs(sln.h1_norm() - 0.89534) < 1e-4

def test_example_04():
    from hermes2d.examples.c04 import set_bc

    set_verbose(False)

    # Below you can play with the parameters CONST_F, P_INIT, and UNIFORM_REF_LEVEL.
    INIT_REF_NUM = 2         # number of initial uniform mesh refinements
    P_INIT = 2               # initial polynomial degree in all elements

    # Load the mesh file
    mesh = Mesh()
    mesh.load(get_example_mesh())

    # Perform initial mesh refinements
    for i in range(INIT_REF_NUM):
        mesh.refine_all_elements()

    # Create an H1 space with default shapeset
    space = H1Space(mesh, P_INIT)
    set_bc(space)

    # Initialize the weak formulation
    wf = WeakForm()
    set_forms(wf)

    # Initialize the linear system
    ls = LinSystem(wf)
    ls.set_spaces(space)

    # Assemble and solve the matrix problem
    sln = Solution()
    ls.assemble()
    ls.solve_system(sln)
    #assert abs(sln.l2_norm() - 1.22729) < 1e-4
    #assert abs(sln.h1_norm() - 2.90006) < 1e-4

def test_example_05():
    from hermes2d.examples.c05 import set_bc
    from hermes2d.examples.c05 import set_forms as set_forms_surf

    set_verbose(False)

    P_INIT = 4                           # initial polynomial degree in all elements
    CORNER_REF_LEVEL = 12                # number of mesh refinements towards the re-entrant corner

    # Load the mesh file
    mesh = Mesh()
    mesh.load(get_example_mesh())

    # Perform initial mesh refinements.
    mesh.refine_towards_vertex(3, CORNER_REF_LEVEL)

    # Create an H1 space with default shapeset
    space = H1Space(mesh, P_INIT)
    set_bc(space)

    # Initialize the weak formulation
    wf = WeakForm()
    set_forms(wf)

    # Initialize the linear system.
    ls = LinSystem(wf)
    ls.set_spaces(space)

    # Assemble and solve the matrix problem
    sln = Solution()
    ls.assemble()
    ls.solve_system(sln)
    #assert abs(sln.l2_norm() - 0.535833) < 1e-4
    #assert abs(sln.h1_norm() - 1.332908) < 1e-4

def test_example_06():
    from hermes2d.examples.c06 import set_bc, set_forms

    set_verbose(False)

    # The following parameters can be changed:

    UNIFORM_REF_LEVEL = 2;   # Number of initial uniform mesh refinements.
    CORNER_REF_LEVEL = 12;   # Number of mesh refinements towards the re-entrant corner.
    P_INIT = 6;              # Uniform polynomial degree of all mesh elements.

    # Boundary markers
    NEWTON_BDY = 1

    # Load the mesh file
    mesh = Mesh()
    mesh.load(get_example_mesh())

    # Perform initial mesh refinements.
    for i in range(UNIFORM_REF_LEVEL):
        mesh.refine_all_elements()
    mesh.refine_towards_vertex(3, CORNER_REF_LEVEL)

    # Create an H1 space with default shapeset
    space = H1Space(mesh, P_INIT)
    set_bc(space)

    # Initialize the weak formulation
    wf = WeakForm()
    set_forms(wf)

    # Initialize the linear system.
    ls = LinSystem(wf)
    ls.set_spaces(space)

    # Assemble and solve the matrix problem
    sln = Solution()
    ls.assemble()
    ls.solve_system(sln)
    #assert abs(sln.l2_norm() - 121.78788) < 1e-4
    #assert abs(sln.h1_norm() - 126.96528) < 1e-4

def test_example_07():
    from hermes2d.examples.c07 import set_bc, set_forms

    set_verbose(False)

    # The following parameters can be changed:
    P_INIT = 2             # Initial polynomial degree of all mesh elements.
    INIT_REF_NUM = 4       # Number of initial uniform refinements

    # Load the mesh
    mesh = Mesh()
    mesh.load(get_07_mesh())

    # Perform initial mesh refinements.
    for i in range(INIT_REF_NUM):
        mesh.refine_all_elements()

    # Create an H1 space with default shapeset
    space = H1Space(mesh, P_INIT)
    set_bc(space)

    # Initialize the weak formulation
    wf = WeakForm()
    set_forms(wf)

    # Initialize the linear system.
    ls = LinSystem(wf)
    ls.set_spaces(space)

    # Assemble and solve the matrix problem
    sln = Solution()
    ls.assemble()
    ls.solve_system(sln)

def test_example_08():
    from hermes2d.examples.c08 import set_bc, set_forms

    set_verbose(False)

    # The following parameter can be changed:
    P_INIT = 4

    # Load the mesh file
    mesh = Mesh()
    mesh.load(get_sample_mesh())

    # Perform uniform mesh refinement
    mesh.refine_all_elements()

    # Create the x- and y- displacement space using the default H1 shapeset
    xdisp = H1Space(mesh, P_INIT)
    ydisp = H1Space(mesh, P_INIT)
    set_bc(xdisp, ydisp)

    # Initialize the weak formulation
    wf = WeakForm(2)
    set_forms(wf)

    # Initialize the linear system.
    ls = LinSystem(wf)
    ls.set_spaces(xdisp, ydisp)

    # Assemble and solve the matrix problem
    xsln = Solution()
    ysln = Solution()
    ls.assemble()
    ls.solve_system(xsln, ysln, lib="scipy")

def test_example_09():
    from hermes2d.examples.c09 import set_bc, temp_ext, set_forms

    # The following parameters can be changed:
    INIT_REF_NUM = 4      # number of initial uniform mesh refinements
    INIT_REF_NUM_BDY = 1  # number of initial uniform mesh refinements towards the boundary
    P_INIT = 4            # polynomial degree of all mesh elements
    TAU = 300.0           # time step in seconds

    # Problem constants
    T_INIT = 10           # temperature of the ground (also initial temperature)
    FINAL_TIME = 86400    # length of time interval (24 hours) in seconds

    # Global variable
    TIME = 0;

    # Boundary markers.
    bdy_ground = 1
    bdy_air = 2

    # Load the mesh
    mesh = Mesh()
    mesh.load(get_cathedral_mesh())

    # Perform initial mesh refinements
    for i in range(INIT_REF_NUM):
        mesh.refine_all_elements()
    mesh.refine_towards_boundary(bdy_air, INIT_REF_NUM_BDY)

    # Create an H1 space with default shapeset
    space = H1Space(mesh, P_INIT)
    set_bc(space)

    # Set initial condition
    tsln = Solution()
    tsln.set_const(mesh, T_INIT)

    # Initialize the weak formulation
    wf = WeakForm()
    set_forms(wf)

    # Initialize the linear system.
    ls = LinSystem(wf)
    ls.set_spaces(space)

    # Time stepping
    nsteps = int(FINAL_TIME/TAU + 0.5)
    rhsonly = False;

    # Assemble and solve
    ls.assemble()
    rhsonly = True
    ls.solve_system(tsln, lib="scipy")

def test_example_10():
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

    # Initialize refinement selector.
    selector = H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER)

    # Initialize the linear system.
    ls = LinSystem(wf)
    ls.set_spaces(space)

    sln_coarse = Solution()
    sln_fine = Solution()

    # Assemble and solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(sln_fine)
    
    # Either solve on coarse mesh or project the fine mesh solution 
    # on the coarse mesh.
    if SOLVE_ON_COARSE_MESH:
        ls.assemble()
        ls.solve_system(sln_coarse)

    # Calculate element errors and total error estimate
    hp = H1Adapt(ls)
    hp.set_solutions([sln_coarse], [sln_fine])
    err_est = hp.calc_error() * 100
            
def test_example_11():
    from hermes2d.examples.c11 import set_bc, set_wf_forms, set_hp_forms

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

    # Initialize refinement selector
    selector = H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER)

    # Initialize the coarse mesh problem
    ls = LinSystem(wf)
    ls.set_spaces(uspace, vspace)

    u_sln_coarse = Solution()
    v_sln_coarse = Solution()
    u_sln_fine = Solution()
    v_sln_fine = Solution()
    
    # Assemble and Solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(u_sln_fine, v_sln_fine, lib="scipy")

    # Either solve on coarse mesh or project the fine mesh solution 
    # on the coarse mesh.
    if SOLVE_ON_COARSE_MESH:
        ls.assemble()
        ls.solve_system(u_sln_coarse, v_sln_coarse, lib="scipy")

    # Calculate element errors and total error estimate
    hp = H1Adapt(ls)
    hp.set_solutions([u_sln_coarse, v_sln_coarse], [u_sln_fine, v_sln_fine]);
    set_hp_forms(hp)
    err_est = hp.calc_error() * 100

def test_example_12():
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

    # Initialize refinement selector
    selector = H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER)

    # Initialize the linear system.
    ls = LinSystem(wf)
    ls.set_spaces(space)

    sln_coarse = Solution()
    sln_fine = Solution()

    # Assemble and solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(sln_fine)
    
    # Either solve on coarse mesh or project the fine mesh solution 
    # on the coarse mesh.   
    if SOLVE_ON_COARSE_MESH:
        ls.assemble()
        ls.solve_system(sln_coarse)

    # Calculate error estimate wrt. fine mesh solution
    hp = H1Adapt(ls)
    hp.set_solutions([sln_coarse], [sln_fine])
    err_est = hp.calc_error() * 100

def test_example_22():
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

    # Initialize refinement selector
    selector = H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER)

    # Initialize the coarse mesh problem
    ls = LinSystem(wf)
    ls.set_spaces(space)

    # Adaptivity loop
    iter = 0
    done =  False
    sln_coarse = Solution()
    sln_fine = Solution()
      
    # Assemble and solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(sln_fine)
    
    # Either solve on coarse mesh or project the fine mesh solution 
    # on the coarse mesh.   
    if SOLVE_ON_COARSE_MESH:
        ls.assemble()
        ls.solve_system(sln_coarse)

    # Calculate error estimate wrt. fine mesh solution
    hp = H1Adapt(ls)
    hp.set_solutions([sln_coarse], [sln_fine])
    err_est = hp.calc_error() * 100
