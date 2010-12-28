# elasticity
add_subdirectory(bracket)
add_subdirectory(crack)

# Thermoelasticity
add_subdirectory(thermoelasticity)

# Maxwell
add_subdirectory(waveguide)

# Fluid mechanics
add_subdirectory(linear-advection-diffusion)
add_subdirectory(linear-advection-dg-adapt)
add_subdirectory(euler-gamm-explicit)
add_subdirectory(euler-gamm-implicit)
add_subdirectory(euler-forward-step)
add_subdirectory(euler-gamm-adapt)

# Neutronics
add_subdirectory(iron-water)
add_subdirectory(saphir)
add_subdirectory(neutronics-4-group)
add_subdirectory(neutronics-4-group-adapt)

# Nernst-Planck
add_subdirectory(newton-np-timedep-adapt-system)

# Navier-Stokes
add_subdirectory(ns-bearing)
add_subdirectory(ns-timedep-adapt)
add_subdirectory(ns-two-phase)

# Schroedinger
add_subdirectory(quantum-billiard)

# Richards
add_subdirectory(richards)
add_subdirectory(richards-tracy-adapt)
add_subdirectory(richards-seepage-adapt)

# Eigenproblems
add_subdirectory(schroedinger-eigen)

# Miscellaneous
add_subdirectory(singular-perturbation)
add_subdirectory(optimal-meshes)

