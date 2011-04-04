# Benchmarks.
add_subdirectory(smooth-iso)
add_subdirectory(smooth-aniso-x)
add_subdirectory(smooth-aniso-y)
add_subdirectory(lshape)
add_subdirectory(layer-interior)
add_subdirectory(layer-boundary)
#add_subdirectory(neutronics-heat-conduction) TODO: Confert to new forms.
#add_subdirectory(neutronics-heat-conduction-adapt) TODO: Confert to new forms.
#add_subdirectory(neutronics-2-group-adapt) TODO: Confert to new forms.
#add_subdirectory(stabilized-advection-diffusion) TODO: Confert to new forms.
#add_subdirectory(stabilized-advection-reaction) TODO: Confert to new forms.
#add_subdirectory(sdirk-22) TODO: Confert to new forms.
add_subdirectory(nonsym-check)

#if(NOT WITH_TRILINOS)
  add_subdirectory(screen)
#endif(NOT WITH_TRILINOS)
