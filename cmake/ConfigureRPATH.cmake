# Required to find the required .so files in correct folders both during 
# local building and after installing to a custom location.
#
# (Adapted from http://www.vtk.org/Wiki/CMake_RPATH_handling)
#

macro (ConfigureRPATH)

    if(${ARGC} EQUAL 1)
        set(ADDITIONAL_RPATH ${ARGV0})
    endif(${ARGC} EQUAL 1) 
  
    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH yet
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # the RPATH to be used when installing, but only if it's not a system directory
    if (ADDITIONAL_PATH)
        LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${ADDITIONAL_RPATH}" isSystemDir)
        IF("${isSystemDir}" STREQUAL "-1")
           SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} ${ADDITIONAL_RPATH})
        ENDIF("${isSystemDir}" STREQUAL "-1")
    endif (ADDITIONAL_PATH)
        
endmacro (ConfigureRPATH)
