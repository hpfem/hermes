#
# GLUT
#

if(WIN64)
  FIND_LIBRARY(GLUT_LIBRARY NAMES freeglut glut PATHS ${GLUT_ROOT}/lib/x64 ${GLUT_ROOT}/lib)
else(WIN64)  
  FIND_LIBRARY(GLUT_LIBRARY NAMES freeglut glut PATHS ${GLUT_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
  if(WIN32)
  else(WIN32)
    FIND_LIBRARY(GL_LIBRARY NAMES GL PATHS ${GLUT_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
  endif(WIN32)
endif(WIN64)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GLUT DEFAULT_MSG GLUT_LIBRARY)
if(WIN32)
else(WIN32)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(GLUT DEFAULT_MSG GL_LIBRARY)
endif(WIN32)