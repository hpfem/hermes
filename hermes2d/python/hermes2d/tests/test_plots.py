disabled=True

import matplotlib
matplotlib.use('Agg')

from hermes2d import (Mesh, H1Shapeset, PrecalcShapeset, H1Space, WeakForm,
        Solution, ScalarView, LinSystem, DummySolver, raises, MeshView,
        set_verbose, plot_mesh_mpl)
from hermes2d.forms import set_forms
from hermes2d.examples import get_example_mesh
set_verbose(False)

domain_mesh = get_example_mesh()

def test_ScalarView_mpl_default():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_element_id(0)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(5)
    space.assign_dofs()

    # initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)

    # assemble the stiffness matrix and solve the system
    sys.assemble()
    A = sys.get_matrix()
    b = sys.get_rhs()
    from scipy.sparse.linalg import cg
    x, res = cg(A, b)
    sln = Solution()
    sln.set_fe_solution(space, pss, x)

    view = ScalarView("Solution")
    view.show(sln, show=False)

def test_ScalarView_mpl_default():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_element_id(0)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(5)
    space.assign_dofs()

    # initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)

    # assemble the stiffness matrix and solve the system
    sys.assemble()
    A = sys.get_matrix()
    b = sys.get_rhs()
    from scipy.sparse.linalg import cg
    x, res = cg(A, b)
    sln = Solution()
    sln.set_fe_solution(space, pss, x)

    view = ScalarView("Solution")
    view.show(sln, show=False, method="contour")

def test_ScalarView_mpl_unknown():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_element_id(0)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(5)
    space.assign_dofs()

    # initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)

    # assemble the stiffness matrix and solve the system
    sys.assemble()
    A = sys.get_matrix()
    b = sys.get_rhs()
    from scipy.sparse.linalg import cg
    x, res = cg(A, b)
    sln = Solution()
    sln.set_fe_solution(space, pss, x)

    view = ScalarView("Solution")
    #This fails currently:
    #assert raises(ValueError, 'view.show(sln, show=False, method="something_unknown_123")')

def test_plot_mesh1a():
    mesh = Mesh()
    mesh.load(domain_mesh)

    view = MeshView("Solution")
    view.show(mesh, lib="mpl", method="simple", show=False)
    plot_mesh_mpl(mesh.nodes, mesh.elements)
    plot_mesh_mpl(mesh.nodes_dict, mesh.elements)
    plot_mesh_mpl(mesh.nodes, mesh.elements, plot_nodes=False)
    plot_mesh_mpl(mesh.nodes_dict, mesh.elements, plot_nodes=False)

def test_plot_mesh1b():
    mesh = Mesh()
    mesh.load(domain_mesh)

    view = MeshView("Solution")
    view.show(mesh, lib="mpl", method="orders", show=False)
    plot_mesh_mpl(mesh.nodes, mesh.elements)
    plot_mesh_mpl(mesh.nodes_dict, mesh.elements)

def test_plot_mesh1c():
    mesh = Mesh()
    mesh.load(domain_mesh)

    view = MeshView("Solution")
    assert raises(ValueError, 'view.show(mesh, lib="mpl", method="something_unknown_123")')

def test_plot_mesh2():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_element_id(0)

    view = MeshView("Solution")
    view.show(mesh, lib="mpl", method="simple", show=False)
    plot_mesh_mpl(mesh.nodes_dict, mesh.elements)
    plot_mesh_mpl(mesh.nodes_dict, mesh.elements, plot_nodes=False)
    view.show(mesh, lib="mpl", method="orders", show=False)
    plot_mesh_mpl(mesh.nodes_dict, mesh.elements)

def test_plot_mesh3a():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_all_elements()
    mesh.refine_all_elements()

    view = MeshView("Solution")
    view.show(mesh, lib="mpl", show=False)

def test_plot_mesh3b():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_all_elements()
    mesh.refine_all_elements()

    plot_mesh_mpl(mesh.nodes_dict, mesh.elements)

def test_plot_mesh3c():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_all_elements()
    mesh.refine_all_elements()

    plot_mesh_mpl(mesh.nodes_dict, mesh.elements, plot_nodes=False)

def test_plot_mesh3d():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_all_elements()
    mesh.refine_all_elements()

    view = MeshView("Solution")
    view.show(mesh, lib="mpl", method="orders", show=False)

def test_plot_mesh3e():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_all_elements()
    mesh.refine_all_elements()

    plot_mesh_mpl(mesh.nodes_dict, mesh.elements)

def test_plot_mesh4():
    mesh = Mesh()
    mesh.load(domain_mesh)

    view = MeshView("Solution")
    view.show(mesh, lib="mpl", show=False, method="orders")
