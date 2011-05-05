# Trilinos
if(WITH_TRILINOS)
  add_subdirectory(01-trilinos-linear)
  add_subdirectory(02-trilinos-nonlinear)
  add_subdirectory(03-trilinos-timedep)
  add_subdirectory(04-trilinos-adapt)
  #add_subdirectory(05-trilinos-coupled) TODO: Convert to new forms.
endif(WITH_TRILINOS)

