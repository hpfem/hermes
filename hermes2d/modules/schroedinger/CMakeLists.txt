add_subdirectory(oscillator)

# Create file with preprocessor definitions exposing the build settings to the source code.
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

project(module-schroedinger)

set(SRC
    schroedinger.cpp
    )

include_directories(${hermes2d_SOURCE_DIR}/src)

add_library(hermes_schroedinger SHARED ${SRC})
target_link_libraries(hermes_schroedinger ${HERMES2D_REAL})

install(TARGETS hermes_schroedinger
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
