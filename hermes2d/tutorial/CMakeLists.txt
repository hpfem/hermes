# Linear problems
add_subdirectory(01-mesh)
add_subdirectory(02-space)
add_subdirectory(03-poisson)
add_subdirectory(04-bc-dirichlet)
add_subdirectory(05-bc-neumann)
add_subdirectory(06-bc-newton)
add_subdirectory(07-general)
add_subdirectory(08-system)
add_subdirectory(09-timedep-basic)
add_subdirectory(09-timedep-rk)

# Adaptivity
add_subdirectory(10-adapt)
add_subdirectory(11-system-adapt)
add_subdirectory(12-general-adapt)
add_subdirectory(13-complex-adapt)
add_subdirectory(14-hcurl-adapt)

# Nonlinear problems
add_subdirectory(15-picard)
add_subdirectory(16-newton-1)
add_subdirectory(17-newton-2)
add_subdirectory(18-newton-adapt)
add_subdirectory(19-newton-timedep-heat-basic)
add_subdirectory(19-newton-timedep-heat-rk)
add_subdirectory(20-newton-timedep-flame)
add_subdirectory(21-newton-timedep-ns)
add_subdirectory(22-newton-timedep-gp)
add_subdirectory(23-newton-timedep-heat-adapt-basic)
add_subdirectory(23-newton-timedep-heat-adapt-rk)
add_subdirectory(24-heat-and-moisture-adapt)
add_subdirectory(25-newton-timedep-gp-adapt)

# Miscellaneous
add_subdirectory(30-space-hcurl)
add_subdirectory(31-space-hdiv)
add_subdirectory(32-space-l2)
add_subdirectory(33-exact-adapt)
if(H2D_WITH_GLUT)
  add_subdirectory(34-remote-computing)
endif(H2D_WITH_GLUT)
add_subdirectory(35-markers-typical)
add_subdirectory(36-markers-nontypical)
add_subdirectory(37-nurbs)
add_subdirectory(38-matrix-solvers)

# Trilinos
if(WITH_TRILINOS)
  add_subdirectory(40-trilinos-linear)
  add_subdirectory(41-trilinos-nonlinear)
  add_subdirectory(42-trilinos-timedep)
  add_subdirectory(43-trilinos-adapt)
  add_subdirectory(44-trilinos-coupled)
endif(WITH_TRILINOS)

# Eigenvalue problems
add_subdirectory(50-eigenvalue)
add_subdirectory(51-eigenvalue-adapt)
add_subdirectory(52-eigenvalue-h2p)

# DG
add_subdirectory(60-linear-advection-dg)

