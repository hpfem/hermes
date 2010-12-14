project(module-schroedinger)

set(SRC
    schroedinger.cpp
    forms.cpp
    )

include_directories(${hermes2d_SOURCE_DIR}/src)

add_library(hermes_schroedinger SHARED ${SRC})
target_link_libraries(hermes_schroedinger ${HERMES_REAL_BIN})

install(TARGETS hermes_schroedinger
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
