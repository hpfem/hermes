# Miscellaneous techniques
#add_subdirectory(05-space-hcurl) TODO: Convert to new forms.
#add_subdirectory(10-space-hdiv) TODO: Convert to new forms.
add_subdirectory(15-space-l2)
if(H2D_WITH_GLUT)
  add_subdirectory(20-remote-computing)
endif(H2D_WITH_GLUT)
#add_subdirectory(25-markers) ALMOST READY BUT THROWS ERRORS
add_subdirectory(30-nurbs)
add_subdirectory(35-matrix-solvers)
