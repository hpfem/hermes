#
# HDF5
#

#FIND_PATH(HDF5_INCLUDE_DIR hdf5.h /usr/include/ /usr/local/include/hdf5 ${HDF5_ROOT})
FIND_PATH(HDF5_INCLUDE_DIR hdf5.h ${HDF5_ROOT}/include)

#FIND_LIBRARY(HDF5_LIBRARY hdf5 /usr/lib/ /usr/local/lib/hdf5 ${HDF5_ROOT}) 

if(64_BIT)
  FIND_LIBRARY(HDF5_LIBRARY hdf5 ${HDF5_ROOT}/lib/x64)
else(64_BIT)  
  FIND_LIBRARY(HDF5_LIBRARY hdf5 ${HDF5_ROOT}/lib)
endif(64_BIT)

#FIND_FILE(HDF5_LIBRARY libhdf5.a /usr/lib /usr/local/lib/hfd5 ${HDF5_ROOT}) 

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( HDF5 HDF5_LIBRARY HDF5_INCLUDE_DIR )
