import agros2d

# problem
problem = agros2d.problem(clear = True)
problem.coordinate_type = "planar"
problem.mesh_type = "triangle"
problem.matrix_solver = "umfpack"
problem.time_step_method = "fixed"
problem.time_method_order = 2
problem.time_method_tolerance = 1
problem.time_total = 0.001
problem.time_steps = 250

# disable view
agros2d.view.mesh.initial_mesh = False
agros2d.view.mesh.solution_mesh = False
agros2d.view.mesh.order = False
agros2d.view.post2d.scalar = False
agros2d.view.post2d.contours = False
agros2d.view.post2d.vectors = False

# fields
# acoustic
acoustic = agros2d.field("acoustic")
acoustic.analysis_type = "transient"
acoustic.initial_condition = 0
acoustic.number_of_refinements = 0
acoustic.polynomial_order = 2
acoustic.linearity_type = "linear"
acoustic.adaptivity_type = "disabled"

# boundaries
acoustic.add_boundary("Matched bundary", "acoustic_impedance", {"acoustic_impedance" : 345*1.25 })
acoustic.add_boundary("Source", "acoustic_pressure", {"acoustic_pressure_real" : { "expression" : "sin(2*pi*(time/(1.0/1000)))" }, "acoustic_pressure_time_derivative" : { "expression" : "2*pi*(1.0/(1.0/1000))*cos(2*pi*(time/(1.0/1000)))" }})
acoustic.add_boundary("Hard wall", "acoustic_normal_acceleration", {"acoustic_normal_acceleration_real" : 0})
acoustic.add_boundary("Soft wall", "acoustic_pressure", {"acoustic_pressure_real" : 0, "acoustic_pressure_time_derivative" : 0})

# materials
acoustic.add_material("Air", {"acoustic_density" : 1.25, "acoustic_speed" : 343})

# geometry
geometry = agros2d.geometry
geometry.add_edge(-0.4, 0.05, 0.1, 0.2, boundaries = {"acoustic" : "Matched bundary"})
geometry.add_edge(0.1, -0.2, -0.4, -0.05, boundaries = {"acoustic" : "Matched bundary"})
geometry.add_edge(-0.4, 0.05, -0.4, -0.05, boundaries = {"acoustic" : "Soft wall"})
geometry.add_edge(-0.18, -0.06, -0.17, -0.05, angle = 90, boundaries = {"acoustic" : "Source"})
geometry.add_edge(-0.17, -0.05, -0.18, -0.04, angle = 90, boundaries = {"acoustic" : "Source"})
geometry.add_edge(-0.18, -0.04, -0.19, -0.05, angle = 90, boundaries = {"acoustic" : "Source"})
geometry.add_edge(-0.19, -0.05, -0.18, -0.06, angle = 90, boundaries = {"acoustic" : "Source"})
geometry.add_edge(0.1, -0.2, 0.1, 0.2, angle = 90, boundaries = {"acoustic" : "Matched bundary"})
geometry.add_edge(0.03, 0.1, -0.04, -0.05, angle = 90, boundaries = {"acoustic" : "Hard wall"})
geometry.add_edge(-0.04, -0.05, 0.08, -0.04, boundaries = {"acoustic" : "Hard wall"})
geometry.add_edge(0.08, -0.04, 0.03, 0.1, boundaries = {"acoustic" : "Hard wall"})

geometry.add_label(-0.0814934, 0.0707097, area = 10e-05, materials = {"acoustic" : "Air"})
geometry.add_label(-0.181474, -0.0504768, materials = {"acoustic" : "none"})
geometry.add_label(0.0314514, 0.0411749, materials = {"acoustic" : "none"})

agros2d.view.zoom_best_fit()

# solve problem
problem.solve()

# point 
point = acoustic.local_values(0.042132, -0.072959)
testp = agros2d.test("Acoustic pressure", point["pr"], 0.200436)
# testSPL = agros2d.test("Acoustic sound level", point["SPL"], 77.055706)

print("Test: Acoustic - transient - planar: " + str(testp))