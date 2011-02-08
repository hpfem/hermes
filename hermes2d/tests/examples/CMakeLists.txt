# Higher-order Runge-Kutta methods (via Butcher's tables)
add_subdirectory(butcher)
add_subdirectory(butcher-adapt)
add_subdirectory(butcher-comparison)

# Elasticity
add_subdirectory(elasticity-bracket)
add_subdirectory(elasticity-crack)

# Thermoelasticity
add_subdirectory(thermoelasticity)

# Maxwell
add_subdirectory(maxwell-waveguide)

# Fluid mechanics
add_subdirectory(linear-advection-diffusion)
add_subdirectory(euler-gamm-explicit)

# Miscellaneous
add_subdirectory(singular-perturbation)
add_subdirectory(optimal-meshes)

# Neutronics
add_subdirectory(neutronics-iron-water)
add_subdirectory(neutronics-saphir)
add_subdirectory(neutronics-4-group)
add_subdirectory(neutronics-4-group-adapt)

# Nernst-Planck
add_subdirectory(nernst-planck-timedep-adapt)
add_subdirectory(nernst-planck-ipmc)

# Navier-Stokes
add_subdirectory(navier-stokes-bearing)
add_subdirectory(navier-stokes-timedep-adapt)
add_subdirectory(navier-stokes-two-phase)

# Schroedinger
add_subdirectory(schroedinger-quantum-billiard)

# Richards
add_subdirectory(richards)
add_subdirectory(richards-tracy-adapt)
add_subdirectory(richards-seepage-adapt)
#add_subdirectory(richards-capillary-barrier-adapt)
#add_subdirectory(richards-capillary-barrier-rk)

# Eigenproblems
add_subdirectory(eigen-schroedinger)
add_subdirectory(eigen-schroedinger-iter)

# Helmholtz
add_subdirectory(helmholtz-waveguide)

