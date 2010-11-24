from hermes2d.modules.electrostatics import Electrostatics

def test_basic():
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
    phi = e.calculate()
    print phi
