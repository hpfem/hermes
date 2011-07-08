add_subdirectory(linear-solvers)
# unistd.h needed for this test does not have to be present.
if (NOT MSVC)
	add_subdirectory(timer)
endif( NOT MSVC)
