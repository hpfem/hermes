add_subdirectory(adapt-quad-1)
add_subdirectory(adapt-quad-2)
#add_subdirectory(eigen-adapt-iter) TODO: Convert to new forms.
add_subdirectory(optimal-meshes)
#add_subdirectory(singular-perturbation) TODO: Convert to new forms.
if(H2D_WITH_GLUT)
  add_subdirectory(view)
endif(H2D_WITH_GLUT)
