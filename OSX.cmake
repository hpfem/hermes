#

## ADDED BY COLMAN

## MODIFIED BY MILAN 

## BLAS
FIND_LIBRARY(BLAS_LIBRARY NAMES blas  PATHS /usr/lib)

IF (H2D_WITH_GLUT)
   
    ## OSMesa
    FIND_LIBRARY(MESA_LIBRARY NAMES OSMesa PATHS /usr/X11/lib /opt/local/lib)

    ## GLU 
    FIND_LIBRARY(GLU_LIBRARY NAMES GLU PATHS /usr/X11/lib /opt/local/lib)

    SET (ADDITIONAL_LIBS ${MESA_LIBRARY} ${GLU_LIBRARY})

    ## FIND_LIBRARY does not work here
    FIND_FILE(GLUT_FILE "libglut.dylib" /usr/lib /opt/local/lib)

    IF (${GLUT_FILE} MATCHES ".*libglut.dylib")
        SET (GLUT_LIBRARY ${GLUT_FILE} CACHE STRING "libglut shows strange behaviour" FORCE)
        EXECUTE_PROCESS(COMMAND file ${GLUT_LIBRARY} OUTPUT_VARIABLE RESULT)

        MESSAGE(STATUS "Using ${GLUT_LIBRARY} as ${RESULT}")
        
        ## The output is seen as one string
        EXECUTE_PROCESS(COMMAND "nm" ${GLUT_LIBRARY} OUTPUT_VARIABLE RESULT)
        STRING(REGEX MATCH ".*T([ ])+_glutSetOption.*" SYMBOL ${RESULT})
        IF (${SYMBOL} STREQUAL ${RESULT})
            MESSAGE(STATUS "SYMBOL _glutSetOption found and defined")
        ELSE()
            MESSAGE(FATAL_ERROR "SYMBOL _glutSetOption not found or not defined")
        ENDIF()
    ELSE ()
        MESSAGE(FATAL_ERROR "GLUT_LIBRARY was not found! ${GLUT_LIBRARY}")
    ENDIF()

    ## GLEW 
    ## FIND_LIBRARY does not work here
    FIND_FILE(GLEW_FILE "libGLEW.dylib" /opt/local/lib)

    IF (${GLEW_FILE} MATCHES ".*libGLEW.dylib")
        SET (GLEW_LIBRARY ${GLEW_FILE} CACHE STRING "libGLEW shows strange behaviour" FORCE)
    ELSE ()
        MESSAGE(FATAL_ERROR "GLEW_LIBRARY was not found!")
    ENDIF()
    
ENDIF(H2D_WITH_GLUT)
