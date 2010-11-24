from hermes2d.modules.electrostatics import Electrostatics

def test_basic():
    e = Electrostatics()
    e.set_mesh_filename("domain.mesh")
    e.set_initial_poly_degree(4)
