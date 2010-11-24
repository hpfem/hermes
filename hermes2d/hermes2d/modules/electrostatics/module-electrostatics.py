import sys
sys.path.append("../../../")

from hermes2d.modules.electrostatics import Electrostatics
from hermes2d.hermes2d import Linearizer
from hermes2d.plot import sln2png

def main():
    e = Electrostatics()
    e.set_mesh_filename("domain.mesh")
    e.set_initial_mesh_refinement(2)
    e.set_initial_poly_degree(4)
    e.set_material_markers([8, 2])
    e.set_permittivity_array([4, 3.1, 5])
    e.set_charge_density_array([4, 3.1, 5])
    e.set_boundary_markers_value([1, 3])
    e.set_boundary_values([1, 5])
    e.set_boundary_markers_derivative([2, 4])
    e.set_boundary_derivatives([1, 5])
    r, sln = e.calculate()
    assert r is True
    print "Saving solution to 'solution.png'"
    sln2png(sln, "solution.png")

if __name__ == "__main__":
    main()
