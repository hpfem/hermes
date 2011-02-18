#add_subdirectory(euler-forward-step)
#add_subdirectory(euler-gamm-adapt)
add_subdirectory(euler-gamm-explicit)

if (WITH_TRILINOS)
#  add_subdirectory(euler-gamm-implicit)
endif (WITH_TRILINOS)
