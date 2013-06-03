#
# HDF5
#

FIND_PATH(HDF5_INCLUDE_DIR hdf5.h ${HDF5_ROOT}/include /usr/include/ /usr/local/include/hdf5)

if(WIN64)
  FIND_LIBRARY(HDF5_LIBRARY hdf5 ${HDF5_ROOT}/lib/x64 ${HDF5_ROOT}/lib/x64)
else(WIN64)
  FIND_LIBRARY(HDF5_LIBRARY hdf5 ${HDF5_ROOT}/lib /usr/lib /usr/lib/hdf5 /usr/local/lib /usr/local/lib/hfd5)
endif(WIN64)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( HDF5 HDF5_LIBRARY HDF5_INCLUDE_DIR )
