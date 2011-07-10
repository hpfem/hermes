execute_process(
	COMMAND python -c "from distutils.sysconfig import get_python_inc; print get_python_inc()"
	OUTPUT_VARIABLE PYTHON_SYS_PATH
	)
string(STRIP ${PYTHON_SYS_PATH} PYTHON_SYS_PATH)
FIND_PATH(PYTHON_INCLUDE_PATH Python.h
    PATHS ${PYTHON_SYS_PATH}
    NO_DEFAULT_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    )
	
execute_process(
	COMMAND python -c "from distutils.sysconfig import get_config_var; print get_config_var('LIBDIR')"
	OUTPUT_VARIABLE PYTHON_LIB_PATH
	)
string(STRIP ${PYTHON_LIB_PATH} PYTHON_LIB_PATH)
IF(APPLE)
    FIND_LIBRARY(PYTHON_LIBRARY NAMES python)
ELSE(APPLE)
    FIND_LIBRARY(PYTHON_LIBRARY NAMES python2.7 python2.6
        PATHS ${PYTHON_LIB_PATH}
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
    )
ENDIF(APPLE)

# Python packages directory. 
#
# Find the system directory into which Python puts all the modules by default.
execute_process(
    COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" 
    OUTPUT_VARIABLE DEFAULT_PYTHON_INSTALL_PATH 
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

IF("${PYTHON_INSTALL_PATH}" STREQUAL "USE_SYSTEM_PYTHON_DIRECTORY")
    set(PYTHON_INSTALL_PATH ${DEFAULT_PYTHON_INSTALL_PATH})
ELSEIF(NOT PYTHON_INSTALL_PATH)
    # Use a local subdirectory of CMAKE_INSTALL_PREFIX.
    set(PYTHON_INSTALL_PATH lib/python)
ENDIF("${PYTHON_INSTALL_PATH}" STREQUAL "USE_SYSTEM_PYTHON_DIRECTORY") 

# To make hermes a Python module (individual parts of the library would then be
# impported as 'from hermes.hermes2d import' instead of 'from hermes2d import')
#   set(PYTHON_INSTALL_PATH "${PYTHON_INSTALL_PATH}/hermes")

# For MSVC (on Win, the function get_config_var does not accept the parameter 'LIBDIR').
IF(MSVC)	
	IF(PYTHON_ROOT)
		FIND_LIBRARY(PYTHON_LIBRARY NAMES python26_d python27_d PATHS ${PYTHON_ROOT}/lib)
	ENDIF(PYTHON_ROOT)
ENDIF(MSVC)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PYTHON DEFAULT_MSG PYTHON_LIBRARY PYTHON_INCLUDE_PATH)


# Links a Python extension module.
#
# The exact link flags are platform dependent and this macro makes it possible
# to write platform independent cmakefiles. All you have to do is to change
# this:
#
# add_library(simple_wrapper SHARED ${SRC})  # Linux only
# set_target_properties(simple_wrapper PROPERTIES PREFIX "")
#
# to this:
#
# add_python_library(simple_wrapper ${SRC})  # Platform independent
#
# Full example:
#
# set(SRC
#     iso_c_utilities.f90
#     pde_pointers.f90
#     example1.f90
#     example2.f90
#     example_eigen.f90
#     simple.f90
#     simple_wrapper.c
# )
# add_python_library(simple_wrapper ${SRC})

macro(ADD_PYTHON_LIBRARY name)
    # When linking Python extension modules, a special care must be taken about
    # the link flags, which are platform dependent:
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # on Mac, we need to use the "-bundle" gcc flag, which is what MODULE
        # does:
        add_library(${name} MODULE ${ARGN})
        # and "-flat_namespace -undefined suppress" link flags, that we need
        # to add by hand:
        set_target_properties(${name} PROPERTIES
            LINK_FLAGS "-flat_namespace -undefined suppress")
    ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # on Linux, we need to use the "-shared" gcc flag, which is what SHARED
        # does:
        add_library(${name} SHARED ${ARGN})
    ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${name} PROPERTIES PREFIX "")
    IF(${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
            target_link_libraries(${name} ${PYTHON_LIBRARY})
    ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
endmacro(ADD_PYTHON_LIBRARY)
