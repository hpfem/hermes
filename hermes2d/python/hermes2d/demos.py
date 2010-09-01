"""
This file contains a few demo functions for easily showing some nice examples
in hermes and also showing that things work.
"""

def demo_layer(lib="mayavi"):
    """
    Shows off the example 22-layer.

    It adaptively refines things and shows the final solution and a convergence
    graph.
    """
    from hermes2d import (Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space,
           WeakForm, Solution, DummySolver, LinSystem, ScalarView, RefSystem,
           H1Adapt, H1ProjBasedSelector, CandList,
	   set_verbose)
    from hermes2d.examples.c22 import set_bc, set_forms

    set_verbose(False)

    def calc(threshold=0.3, strategy=0, h_only=False, error_tol=1,
            interactive_plotting=False, show_mesh=False, show_graph=True):
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

       mesh.refine_all_elements()

       shapeset = H1Shapeset()
       pss = PrecalcShapeset(shapeset)

       space = H1Space(mesh, shapeset)
       set_bc(space)
       space.set_uniform_order(1)

       wf = WeakForm(1)
       set_forms(wf)

       sln = Solution()
       rsln = Solution()
       solver = DummySolver()

       selector = H1ProjBasedSelector(CandList.HP_ANISO, 1.0, -1, shapeset)

       view = ScalarView("Solution")
       iter = 0
       graph = []
       while 1:
           space.assign_dofs()

           sys = LinSystem(wf, solver)
           sys.set_spaces(space)
           sys.set_pss(pss)
           sys.assemble()
           sys.solve_system(sln)
           dofs = sys.get_matrix().shape[0]
           if interactive_plotting:
               view.show(sln, lib=lib, notebook=True,
                       filename="a%02d.png" % iter)

           rsys = RefSystem(sys)
           rsys.assemble()

           rsys.solve_system(rsln)

           hp = H1Adapt([space])
	   hp.set_solutions([sln], [rsln])
	   err_est = hp.calc_error() * 100

           err_est =  hp.calc_error(sln, rsln)*100
           print "iter=%02d, err_est=%5.2f%%, DOFS=%d" % (iter, err_est, dofs)
           graph.append([dofs, err_est])
           if err_est < error_tol:
               break
           hp.adapt(selector, threshold, strategy)
           iter += 1


       if not interactive_plotting:
           view.show(sln, lib=lib, notebook=True)

       if show_mesh:
           mview = MeshView("Mesh")
           mview.show(mesh, lib="mpl", notebook=True, filename="b.png")

       if show_graph:
           from numpy import array
           graph = array(graph)
           import pylab
           pylab.clf()
           pylab.plot(graph[:, 0], graph[:, 1], "ko", label="error estimate")
           pylab.plot(graph[:, 0], graph[:, 1], "k-")
           pylab.title("Error Convergence for the Inner Layer Problem")
           pylab.legend()
           pylab.xlabel("Degrees of Freedom")
           pylab.ylabel("Error [%]")
           pylab.yscale("log")
           pylab.grid()
           pylab.savefig("graph.png")

    calc()
