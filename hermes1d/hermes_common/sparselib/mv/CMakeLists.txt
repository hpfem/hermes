include_directories(..)

add_library(mv STATIC
	mvvd.cc
	mvmd.cc
	mvvf.cc
	mvmf.cc
	mvvdio.cc
	mvblasd.cc
	mvblasf.cc
	mvblasi.cc
	mvvi.cc
	mvmi.cc
	mvvc.cc
	mvmc.cc
	mvvcio.cc
	mvblasc.cc
)

set_property(TARGET mv PROPERTY COMPILE_FLAGS "-g -fPIC")
