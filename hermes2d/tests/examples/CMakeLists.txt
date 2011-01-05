# elasticity
add_subdirectory(elasticity-bracket)
add_subdirectory(elasticity-crack)

# thermoelasticity
add_subdirectory(thermoelasticity)

# maxwell
add_subdirectory(maxwell-waveguide)
add_subdirectory(linear-advection-diffusion)
add_subdirectory(singular-perturbation)
add_subdirectory(optimal-meshes)

# neutronics
add_subdirectory(neutronics-iron-water)
add_subdirectory(neutronics-saphir)
add_subdirectory(neutronics-4-group)
add_subdirectory(neutronics-4-group-adapt)

# nernst-planck
add_subdirectory(nernst-planck-timedep-adapt)

# navier-stokes
add_subdirectory(navier-stokes-bearing)
add_subdirectory(navier-stokes-timedep-adapt)
add_subdirectory(navier-stokes-two-phase)

# schroedinger
add_subdirectory(schroedinger-quantum-billiard)

# richards
add_subdirectory(richards)
add_subdirectory(richards-tracy-adapt)
add_subdirectory(richards-seepage-adapt)
