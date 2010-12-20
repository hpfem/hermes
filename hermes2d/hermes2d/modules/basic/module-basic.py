import sys
sys.path.append("../../../")

from hermes2d.modules.basic import ModuleBasic
from hermes2d.hermes2d import Linearizer
from hermes2d.plot import sln2png

def main():
    e = ModuleBasic()
    e.set_mesh_str("\na = 1.0  # size of the mesh\nb = sqrt(2)/2\n\nvertices =\n{\n  { 0, -a },    # vertex 0\n  { a, -a },    # vertex 1\n  { -a, 0 },    # vertex 2\n  { 0, 0 },     # vertex 3\n  { a, 0 },     # vertex 4\n  { -a, a },    # vertex 5\n  { 0, a },     # vertex 6\n  { a*b, a*b }  # vertex 7\n}\n\nelements =\n{\n  { 0, 1, 4, 3, 0 },  # quad 0\n  { 3, 4, 7, 0 },     # tri 1\n  { 3, 7, 6, 0 },     # tri 2\n  { 2, 3, 6, 5, 0 }   # quad 3\n}\n\nboundaries =\n{\n  { 0, 1, 1 },\n  { 1, 4, 2 },\n  { 3, 0, 4 },\n  { 4, 7, 2 },\n  { 7, 6, 2 },\n  { 2, 3, 4 },\n  { 6, 5, 2 },\n  { 5, 2, 3 }\n}\n\ncurves =\n{\n  { 4, 7, 45 },  # +45 degree circular arcs\n  { 7, 6, 45 }\n}\n");
    e.set_initial_mesh_refinement(2)
    e.set_initial_poly_degree(4)
    e.set_matrix_solver("umfpack")
    e.set_material_markers([0])
    e.set_c1_array([1])
    e.set_c2_array([0])
    e.set_c3_array([0])
    e.set_c4_array([0])
    e.set_c5_array([1])
    e.set_dirichlet_markers([4])
    e.set_dirichlet_values([4], [0])
    e.set_neumann_markers([1, 3])
    e.set_neumann_values([0, 0])
    e.set_newton_markers([2])
    e.set_newton_values( [(1, 1)] )

    success = e.calculate()
    assert success is True
    sln = e.get_solution()
    print "Assembly time:", e.get_assembly_time()
    print "Solver time:", e.get_solver_time()
    print "Saving solution to 'solution.png'"
    sln2png(sln, "solution.png")

if __name__ == "__main__":
    main()
