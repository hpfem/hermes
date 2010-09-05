project(hermes_common)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# optional functionality
set(COMMON_WITH_SCIPY YES)
set(COMMON_WITH_UMFPACK NO)
set(COMMON_WITH_SUPERLU NO)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

include(UseCython)
# if the user specified PYTHON_INCLUDE_PATH, let's use that, otherwise let's
# try to guess using cmake/FindPythonLibs.cmake
if (NOT PYTHON_INCLUDE_PATH)
    find_package(PythonLibs REQUIRED)
endif (NOT PYTHON_INCLUDE_PATH)
find_package(NumPy REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})

enable_testing()

#PYTHONPATH=${hermes_common_SOURCE_DIR}
set(HERMES_COMMON hermes_common)

set(SRC
    _hermes_common_api_new.cpp
    matrix.cpp
    matrixio.cpp
    solvers.cpp
    python_solvers.cpp
    python_api.cpp
    umfpack_solver.cpp
    superlu_solver.cpp
    sparselib_solver.cpp
    common_time_period.cpp
    )

if(MSVC)
	add_library(${HERMES_COMMON} STATIC ${SRC})
else(MSVC)
	add_library(${HERMES_COMMON} SHARED ${SRC})
endif(MSVC)

set_property(TARGET ${HERMES_COMMON} PROPERTY COMPILE_FLAGS "-g")


set(SPARSELIB_LIBRARIES mv sparse spblas)
set(SPARSELIB_INCLUDE_DIRS sparselib sparselib/mv sparselib/iml)
set(SPARSELIB_DEFINITIONS -DCOMPLEX=std::complex<double>)

add_subdirectory(sparselib)
add_definitions(${SPARSELIB_DEFINITIONS})
include_directories(${SPARSELIB_INCLUDE_DIRS})

target_link_libraries(${HERMES_COMMON} ${PYTHON_LIBRARIES} ${SPARSELIB_LIBRARIES})
set_source_files_properties(matrix.cpp PROPERTIES
    OBJECT_DEPENDS ${hermes_common_SOURCE_DIR}/_hermes_common_api_new.h
    )
#add_custom_command(
#    OUTPUT _hermes_common_api.h
#    COMMAND cython _hermes_common.pyx
#    )

cython_add_module_compile(_hermes_common stdcython.cpp)
target_link_libraries(_hermes_common ${HERMES_COMMON})
set_property(TARGET _hermes_common PROPERTY COMPILE_FLAGS "-g")


if(MSVC)
    IF(NOT PYTHON_VERSION)
        execute_process(COMMAND python -c "import sys;print '%d.%d' % sys.version_info[:2]" OUTPUT_VARIABLE PYTHON_VERSION)
        string(STRIP ${PYTHON_VERSION} PYTHON_VERSION)
    ENDIF(NOT PYTHON_VERSION)

    get_filename_component(THIS_DIR . ABSOLUTE) 

    STRING(REPLACE "/" "\\" THIS_DIR ${THIS_DIR})   #function get_filename_component returns the path in Unix format, so we need to adjust it

    add_custom_command(
	TARGET _hermes_common
	POST_BUILD

	#needed so that the function PyImport_Import can find the _hermes_common library
	COMMAND if exist Debug\\_hermes_common.dll copy /Y Debug\\_hermes_common.dll Debug\\_hermes_common.pyd 
	COMMAND if exist Release\\_hermes_common.dll copy /Y Release\\_hermes_common.dll Release\\_hermes_common.pyd

	#this is a very tricky thing on windows, we set the varible after build, but the variable's value is not updated in the current process
	#and for the same reason it is not possible to set it completely using cmake
        COMMAND setx PYTHONPATH \"%PYTHONPATH%\;${THIS_DIR}\\Debug\;${THIS_DIR}\\Release\" /M
    )
    #if the variable exists, set it empty
    execute_process(COMMAND if defined PYTHONPATH set PYTHONPATH=)
    #fill the variable from windows registry
    execute_process(COMMAND setx PYTHONPATH /K HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\${PYTHON_VERSION}\\PythonPath /M)
endif(MSVC)

if(COMMON_WITH_SCIPY)
    add_definitions(-DCOMMON_WITH_SCIPY)
endif(COMMON_WITH_SCIPY)

if(COMMON_WITH_UMFPACK)
    add_definitions(-DCOMMON_WITH_UMFPACK)
    find_package(UMFPACK REQUIRED)
    include_directories(${UMFPACK_INCLUDE_DIR})
    include_directories(${AMD_INCLUDE_DIR})
    target_link_libraries(${HERMES_COMMON} ${UMFPACK_LIBRARY})
    target_link_libraries(${HERMES_COMMON} ${AMD_LIBRARY})
endif(COMMON_WITH_UMFPACK)

if(COMMON_WITH_SPOOLES)
    add_definitions(-DCOMMON_WITH_SPOOLES)
    find_package(SPOOLES REQUIRED)
    include_directories(${SPOOLES_INCLUDE_DIR})
    target_link_libraries(${HERMES_COMMON} ${SPOOLES_LIBRARY})
endif(COMMON_WITH_SPOOLES)

if(COMMON_WITH_SUPERLU)
    add_definitions(-DCOMMON_WITH_SUPERLU)
    find_package(SUPERLU REQUIRED)
    include_directories(${SUPERLU_INCLUDE_DIR})
    target_link_libraries(${HERMES_COMMON} ${SUPERLU_LIBRARY})
endif(COMMON_WITH_SUPERLU)

if(NOT MSVC)
  IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries(${HERMES_COMMON} rt dl util)
  ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(NOT MSVC)

add_subdirectory(tests)
