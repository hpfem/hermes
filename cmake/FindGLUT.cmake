#
# GLUT
#

if(WIN64)
  FIND_LIBRARY(GLUT_LIBRARY NAMES freeglut glut PATHS ${GLUT_ROOT}/lib/x64 ${GLUT_ROOT}/lib)
else(WIN64)  
  FIND_LIBRARY(GLUT_LIBRARY NAMES freeglut glut PATHS ${GLUT_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
endif(WIN64)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GLUT DEFAULT_MSG GLUT_LIBRARY)