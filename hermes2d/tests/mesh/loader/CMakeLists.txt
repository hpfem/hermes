project(mesh-loader)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

add_test(mesh-loader-h2d-1 "${BIN}" h2d domain.mesh domain.dump)

add_test(mesh-loader-h2d-old-1 "${BIN}" h2d-old old-fmt.mesh old-fmt.dump)
add_test(mesh-loader-h2d-str-1 "${BIN}" h2d-str old-fmt.mesh old-fmt.dump)

if (WITH_EXODUSII)
  add_test(mesh-loader-exodus-ii-1 "${BIN}" exII square.e square.dump)
  add_test(mesh-loader-exodus-ii-2 "${BIN}" exII square-tri.e square-tri.dump)
  add_test(mesh-loader-exodus-ii-3 "${BIN}" exII iron-water.e iron-water.dump)
endif(WITH_EXODUSII)
