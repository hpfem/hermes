from numpy import array

from hermes2d import Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        WeakForm, Solution, ScalarView, set_verbose, LinSystem, DummySolver, \
        set_warn_integration
from hermes2d.forms import set_forms
from hermes2d.examples import get_example_mesh

domain_mesh = get_example_mesh()
set_warn_integration(False)

def equal_arrays(a, b, prec=1e-10):
    if len(a) == len(b):
        d = a - b
        s = sum(d**2)
        if s < prec:
            return True
    return False

def test_matrix():
    set_verbose(False)

    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_element_id(0)

    # create an H1 space with default shapeset
    space = H1Space(mesh, 1)

    # initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf)

    sys = LinSystem(wf)
    sys.set_spaces(space)

    # assemble the stiffness matrix and solve the system
    sln = Solution()
    sys.assemble()
    A = sys.get_matrix()
    #dp.solve_system(sln)

def test_fe_solutions():
    mesh = Mesh()
    mesh.load(domain_mesh)

    space = H1Space(mesh, 1)
    space.set_uniform_order(2)
    space.assign_dofs()

    a = array([1, 2, 3, 8, 0.1])

    sln = Solution()

    # the sln.get_fe_solution() is not yet implemented in hermes2d
    #b = sln.get_fe_solution()
    #assert equal_arrays(a, b)
