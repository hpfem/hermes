#

## ADDED BY COLMAN 

## BLAS
FIND_LIBRARY(BLAS_LIBRARY NAMES blas  PATHS /usr/lib)

IF (${BLAS_LIBRARY} MATCHES ".*libblas.dylib")
    EXECUTE_PROCESS(COMMAND file ${BLAS_LIBRARY} OUTPUT_VARIABLE RESULT)
            
    IF (${RESULT} MATCHES ".*i386.*")
        MESSAGE(STATUS "Using ${BLAS_LIBRARY} with architecture i386")
    ELSE ()
        MESSAGE(FATAL_ERROR "BLAS_LIBRARY ${BLAS_LIBRARY} is NOT i386!")
    ENDIF ()
ELSE ()
    MESSAGE(FATAL_ERROR "BLAS_LIBRARY was not found!")
ENDIF()

IF (H2D_WITH_GLUT)
   
    ## OSMesa
    FIND_LIBRARY(MESA_LIBRARY NAMES OSMesa PATHS /usr/X11/lib)
    IF (${MESA_LIBRARY} MATCHES ".*libOSMesa.dylib")
        EXECUTE_PROCESS(COMMAND file ${MESA_LIBRARY} OUTPUT_VARIABLE RESULT)
            
        IF (${RESULT} MATCHES ".*i386.*")
            MESSAGE(STATUS "Using ${MESA_LIBRARY} with architecture i386")
        ELSE ()
            MESSAGE(FATAL_ERROR "MESA_LIBRARY ${MESA_LIBRARY} is NOT i386!")
        ENDIF ()
    ELSE ()
        MESSAGE(FATAL_ERROR "MESA_LIBRARY was not found!")
    ENDIF()

    ## GLU 
    FIND_LIBRARY(GLU_LIBRARY NAMES GLU PATHS /usr/X11/lib)
    IF (${GLU_LIBRARY} MATCHES ".*libGLU.dylib")
        EXECUTE_PROCESS(COMMAND file ${GLU_LIBRARY} OUTPUT_VARIABLE RESULT)
            
        IF (${RESULT} MATCHES ".*i386.*")
            MESSAGE(STATUS "Using ${GLU_LIBRARY} with architecture i386")
        ELSE ()
            MESSAGE(FATAL_ERROR "GLU_LIBRARY ${MESA_LIBRARY} is NOT i386!")
        ENDIF ()
    ELSE ()
        MESSAGE(FATAL_ERROR "GLU_LIBRARY was not found!")
    ENDIF()

    SET (ADDITIONAL_LIBS ${MESA_LIBRARY} ${GLU_LIBRARY})

    ## FIND_LIBRARY does not work here
    FIND_FILE(GLUT_FILE "libglut.dylib" /usr/lib)

    IF (${GLUT_FILE} MATCHES ".*libglut.dylib")
        SET (GLUT_LIBRARY ${GLUT_FILE} CACHE STRING "libglut shows strange behaviour" FORCE)
        EXECUTE_PROCESS(COMMAND file ${GLUT_LIBRARY} OUTPUT_VARIABLE RESULT)

        IF (${RESULT} MATCHES ".*i386.*")
            MESSAGE(STATUS "Using ${GLUT_LIBRARY} with architecture i386")
            
            ## The output is seen as one string
            EXECUTE_PROCESS(COMMAND "nm" ${GLUT_LIBRARY} OUTPUT_VARIABLE RESULT)
            STRING(REGEX MATCH ".*T([ ])+_glutSetOption.*" SYMBOL ${RESULT})
            IF (${SYMBOL} STREQUAL ${RESULT})
                MESSAGE(STATUS "SYMBOL _glutSetOption found and defined")
            ELSE()
                MESSAGE(FATAL_ERROR "SYMBOL _glutSetOption not found or not defined")
            ENDIF()
        ELSE ()
            MESSAGE(FATAL_ERROR "GLUT_LIBRARY ${BLAS_LIBRARY} is NOT i386!")
        ENDIF ()
    ELSE ()
        MESSAGE(FATAL_ERROR "GLUT_LIBRARY was not found! ${GLUT_LIBRARY}")
    ENDIF()

    ## GLEW 
    ## FIND_LIBRARY does not work here
    FIND_FILE(GLEW_FILE "libGLEW.dylib" /opt/local/lib)

    IF (${GLEW_FILE} MATCHES ".*libGLEW.dylib")
        SET (GLEW_LIBRARY ${GLEW_FILE} CACHE STRING "libGLEW shows strange behaviour" FORCE)
        EXECUTE_PROCESS(COMMAND file ${GLEW_LIBRARY} OUTPUT_VARIABLE RESULT)
            
        IF (${RESULT} MATCHES ".*i386.*")
            MESSAGE(STATUS "Using ${GLEW_LIBRARY} with architecture i386")
        ELSE ()
            MESSAGE(FATAL_ERROR "GLEW_LIBRARY ${GLEW_LIBRARY} is NOT i386!")
        ENDIF ()
    ELSE ()
        MESSAGE(FATAL_ERROR "GLEW_LIBRARY was not found!")
    ENDIF()
    
ENDIF(H2D_WITH_GLUT)
