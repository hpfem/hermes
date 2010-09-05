# This allows to link Cython files
# Examples:
# 1) to compile assembly.pyx to assembly.so:
#   CYTHON_ADD_MODULE(assembly)
# 2) to compile assembly.pyx and something.cpp to assembly.so:
#   CYTHON_ADD_MODULE(assembly something.cpp)

if(NOT CYTHON_INCLUDE_DIRECTORIES)
    set(CYTHON_INCLUDE_DIRECTORIES .)
endif(NOT CYTHON_INCLUDE_DIRECTORIES)

macro(CYTHON_ADD_MODULE_COMPILE name)
    # When linking Python extension modules, a special care must be taken about
    # the link flags, which is platform dependent:
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # on Mac, we need to use the "-bundle" gcc flag, which is what MODULE
        # does:
        add_library(${name} MODULE ${name}.cpp ${ARGN})
        # and "-flat_namespace -undefined suppress" link flags, that we need
        # to add by hand:
        set_target_properties(${name} PROPERTIES
            LINK_FLAGS "-flat_namespace -undefined suppress")
    ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # on Linux, we need to use the "-shared" gcc flag, which is what SHARED
        # does:
        add_library(${name} SHARED ${name}.cpp ${ARGN})
    ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${name} PROPERTIES PREFIX "")
    IF(${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
            target_link_libraries(${name} ${PYTHON_LIBRARY})
    ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endmacro(CYTHON_ADD_MODULE_COMPILE)

macro(CYTHON_ADD_MODULE name)
    add_custom_command(
        OUTPUT ${name}.cpp
        COMMAND cython
        ARGS -I ${CYTHON_INCLUDE_DIRECTORIES} -o ${name}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${name}.pyx
        DEPENDS ${name}.pyx ${name}.pxd
        COMMENT "Cythonizing ${name}.pyx")
    CYTHON_ADD_MODULE_COMPILE(${name} ${ARGN})
endmacro(CYTHON_ADD_MODULE)
