project(hermes)

cmake_minimum_required(VERSION 2.6)

#---------------------------------------------------------
# begin of configuration section
#
# The values below are sane defaults to produce a build on all platforms,
# if you want to override them, set them in CMake.vars

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/hermes_common/cmake)
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/hermes_common)

# Optional parts of the library:
set(WITH_H1D           YES)
set(WITH_H2D           YES)
set(WITH_H3D           YES)

# Allow to override the default values in CMake.vars:
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
  include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

if(WITH_H1D)
    add_subdirectory(hermes1d)
endif(WITH_H1D)

if(WITH_H2D)
    add_subdirectory(hermes2d)
endif(WITH_H2D)

if(WITH_H3D)
    add_subdirectory(hermes3d)
endif(WITH_H3D)

message("\nConfiguration results")
message("---------------------")
message("Build Hermes1D: ${WITH_H1D}")
message("Build Hermes2D: ${WITH_H2D}")
message("Build Hermes3D: ${WITH_H3D}")

message("\n")
