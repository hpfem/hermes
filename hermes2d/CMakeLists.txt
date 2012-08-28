project(hermes2d)
  # Generating classes for XML mesh parsing according to your version of XSD.
  message("Hermes2D: Generating classes for XML mesh parsing according to your version of XSD.")
  if(WIN32)
    execute_process(COMMAND "${PROJECT_SOURCE_DIR}/xml_schemas/xml_schemas_generate_win.bat" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/xml_schemas/" OUTPUT_QUIET)
  else(WIN32)
    execute_process(COMMAND "chmod 700 ${PROJECT_SOURCE_DIR}/xml_schemas/xml_schemas_generate_linux.sh" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/xml_schemas/")
    execute_process(COMMAND "${PROJECT_SOURCE_DIR}/xml_schemas/xml_schemas_generate_linux.sh" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/xml_schemas/" OUTPUT_QUIET)
  endif(WIN32)

  SET(H2D_XML_SCHEMAS_DIRECTORY "\"${PROJECT_SOURCE_DIR}/xml_schemas\"")

  configure_file(
    config.h.cmake
    include/config.h
  )

  # Determine which version (debug x release) to build.
  if (NOT H2D_DEBUG AND NOT H2D_RELEASE)
    message(FATAL_ERROR "WITH_H2D is set to YES, but no versions are specified.
                         Re-run with H2D_DEBUG and/or H2D_RELEASE set to YES.")
  endif (NOT H2D_DEBUG AND NOT H2D_RELEASE)

  if(MSVC)
    set(HERMES2D "hermes2d")              # Base name of the library.
    set(HERMES2D_RELEASE ${HERMES2D})     # Name of the release version.
    set(HERMES2D_DEBUG ${HERMES2D}-debug) # Name of the debug version.
      
  else(MSVC)
    if(H2D_DEBUG)
      set(HERMES2D_DEBUG "hermes2d-debug")
      set(HERMES2D ${HERMES2D_DEBUG})
      # set(CMAKE_BUILD_TYPE Debug)
    # This does not work with multi-version configurations,
    # because it is linked to the project (i.e. we would have
    # to have separate projects for Debug and Release)
    endif(H2D_DEBUG)

    if(H2D_RELEASE)
      set(HERMES2D_RELEASE "hermes2d")
      set(HERMES2D ${HERMES2D_RELEASE})
      # set(CMAKE_BUILD_TYPE Release)
    endif(H2D_RELEASE)
  endif(MSVC)
  
  # Linker settings.
  include(ConfigureRPATH)
  ConfigureRPATH()

  #
  # Optional libraries and settings for H2D.
  #
  if(H2D_WITH_GLUT)
    if(NOT APPLE)   # See hermes/OSX.cmake for APPLE-specific search procedures.
      find_package(GLUT REQUIRED)
      find_package(GLEW REQUIRED)
    endif(NOT APPLE)
  else(H2D_WITH_GLUT)
    add_definitions(-DNOGLUT)
  endif(H2D_WITH_GLUT)

  if(H2D_WITH_VIEWER_GUI)
    find_package(ANTTWEAKBAR REQUIRED)
    include_directories(${ANTTWEAKBAR_INCLUDE_DIR})
  endif(H2D_WITH_VIEWER_GUI)
  
  # Mesh format.
  if(WITH_HDF5)
    find_package(HDF5 REQUIRED)
    include_directories(${HDF5_INCLUDE_DIR})
  endif(WITH_HDF5)

  if(WITH_EXODUSII)
    find_package(EXODUSII REQUIRED)
    include_directories(${EXODUSII_INCLUDE_DIR})
  endif(WITH_EXODUSII)
  include_directories(${XSD_INCLUDE_DIR})
  include_directories(${XERCES_INCLUDE_DIR})
  
  #
  # Header files location for the Hermes2d library.
  #
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes_common/include)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes_common/include/solvers)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/form)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/space)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/mesh)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/boundary_conditions)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/shapeset)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/function)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/integrals)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/linearizer)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/adapt)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/quadrature)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/refinement_selectors)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/views)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/weakform)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/weakform_library)
  include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/projections)

  #
  # Source files for the Hermes2d library.
  #
  set(SRC
    src/forms.cpp
    src/asmlist.cpp
    src/newton_solver.cpp
    src/picard_solver.cpp
    src/linear_solver.cpp
    
    src/calculation_continuity.cpp

    src/adapt/adapt.cpp
    src/adapt/kelly_type_adapt.cpp

    src/boundary_conditions/essential_boundary_conditions.cpp

    src/function/transformable.cpp
    src/function/function.cpp
    src/function/exact_solution.cpp
    src/function/solution.cpp
    src/function/filter.cpp
    src/function/mesh_function.cpp
    src/function/solution_h2d_xml.cpp

    src/api2d.cpp
	src/mixins2d.cpp
    
    src/mesh/refmap.cpp
    src/mesh/curved.cpp
    src/mesh/refinement_type.cpp
    src/mesh/element_to_refine.cpp
    src/mesh/mesh_reader_exodusii.cpp
    src/mesh/hash.cpp
    src/mesh/mesh_reader_h2d.cpp
    src/mesh/mesh_reader_h2d_xml.cpp
    src/mesh/mesh_reader_h1d_xml.cpp
    src/mesh/mesh_h2d_xml.cpp
    src/mesh/mesh_h1d_xml.cpp
    src/mesh/subdomains_h2d_xml.cpp
    src/mesh/mesh.cpp
    src/mesh/traverse.cpp
    src/mesh/mesh_data.cpp

    src/quadrature/limit_order.cpp
    src/quadrature/quad_std.cpp

    src/refinement_selectors/selector.cpp
    src/refinement_selectors/order_permutator.cpp
    src/refinement_selectors/optimum_selector.cpp
    src/refinement_selectors/proj_based_selector.cpp
    src/refinement_selectors/l2_proj_based_selector.cpp
    src/refinement_selectors/h1_proj_based_selector.cpp
    src/refinement_selectors/hcurl_proj_based_selector.cpp

    src/shapeset/shapeset.cpp
    src/shapeset/shapeset_h1_ortho.cpp
    src/shapeset/shapeset_h1_jacobi.cpp
    src/shapeset/shapeset_h1_quad.cpp
    src/shapeset/shapeset_hc_legendre.cpp
    src/shapeset/shapeset_hc_gradleg.cpp
    src/shapeset/shapeset_hd_legendre.cpp
    src/shapeset/shapeset_l2_legendre.cpp
    src/shapeset/precalc.cpp

    src/space/space.cpp
    src/space/space_h1.cpp
    src/space/space_hcurl.cpp
    src/space/space_l2.cpp
    src/space/space_hdiv.cpp
    src/space/space_h2d_xml.cpp

    src/views/base_view.cpp
    src/views/mesh_view.cpp
    src/views/order_view.cpp
    src/views/scalar_view.cpp
    src/views/stream_view.cpp
    src/views/vector_base_view.cpp
    src/views/vector_view.cpp
    src/views/view.cpp
    src/views/view_data.cpp
    src/views/view_support.cpp
    src/views/linearizer.cpp
    src/views/linearizer_base.cpp
    src/views/orderizer.cpp
    src/views/vectorizer.cpp

    src/weakform/weakform.cpp

    src/neighbor.cpp
    src/graph.cpp
    src/global.cpp
    src/discrete_problem.cpp
    src/discrete_problem_linear.cpp
    src/runge_kutta.cpp
    src/spline.cpp

    src/projections/ogprojection.cpp
    src/projections/ogprojection_nox.cpp
    src/projections/localprojection.cpp

    src/weakform_library/weakforms_elasticity.cpp
    src/weakform_library/weakforms_h1.cpp
    src/weakform_library/weakforms_hcurl.cpp
    src/weakform_library/weakforms_maxwell.cpp
    src/weakform_library/weakforms_neutronics.cpp
  )
  
  set(HEADERS
    include/forms.h
    include/asmlist.h
    include/newton_solver.h
    include/picard_solver.h
    include/linear_solver.h

    include/calculation_continuity.h

    include/adapt/adapt.h
    include/adapt/kelly_type_adapt.h

    include/boundary_conditions/essential_boundary_conditions.h

    include/function/transformable.h
    include/function/function.h
    include/function/exact_solution.h
    include/function/solution.h
    include/function/filter.h
    include/function/mesh_function.h
    include/function/solution_h2d_xml.h

    include/integrals/l2.h
    include/integrals/h1.h
    include/integrals/hdiv.h
    include/integrals/hcurl.h
    
    include/api2d.h
	include/mixins2d.h
    
    include/mesh/refmap.h
    include/mesh/curved.h
    include/mesh/refinement_type.h
    include/mesh/element_to_refine.h
    include/mesh/mesh_reader_exodusii.h
    include/mesh/hash.h
    include/mesh/mesh_reader_h2d.h
    include/mesh/mesh_reader_h2d_xml.h
    include/mesh/mesh_reader_h1d_xml.h
    include/mesh/mesh_h2d_xml.h
    include/mesh/mesh_h1d_xml.h
    include/mesh/subdomains_h2d_xml.h
    include/mesh/mesh.h
    include/mesh/traverse.h
    include/mesh/mesh_data.h

    include/quadrature/limit_order.h
    include/quadrature/quad.h
    include/quadrature/quad_all.h

    include/refinement_selectors/selector.h
    include/refinement_selectors/order_permutator.h
    include/refinement_selectors/optimum_selector.h
    include/refinement_selectors/proj_based_selector.h
    include/refinement_selectors/l2_proj_based_selector.h
    include/refinement_selectors/h1_proj_based_selector.h
    include/refinement_selectors/hcurl_proj_based_selector.h

    include/shapeset/shapeset.h
    include/shapeset/shapeset_h1_all.h
    include/shapeset/shapeset_h1_quad.h
    include/shapeset/shapeset_hc_all.h
    include/shapeset/shapeset_hd_all.h
    include/shapeset/shapeset_l2_all.h
    include/shapeset/precalc.h

    include/space/space.h
    include/space/space_h1.h
    include/space/space_hcurl.h
    include/space/space_l2.h
    include/space/space_hdiv.h
    include/space/space_h2d_xml.h

    include/views/base_view.h
    include/views/mesh_view.h
    include/views/order_view.h
    include/views/scalar_view.h
    include/views/stream_view.h
    include/views/vector_base_view.h
    include/views/vector_view.h
    include/views/view.h
    include/views/view_support.h
    include/views/linearizer.h
    include/views/linearizer_base.h
    include/views/orderizer.h
    include/views/vectorizer.h

    include/weakform/weakform.h

    include/neighbor.h
    include/graph.h
    include/global.h
    include/discrete_problem.h
    include/discrete_problem_linear.h
    include/runge_kutta.h
    include/spline.h

    include/projections/ogprojection.h
    include/projections/ogprojection_nox.h
    include/projections/localprojection.h

    include/weakform_library/weakforms_elasticity.h
    include/weakform_library/weakforms_h1.h
    include/weakform_library/weakforms_hcurl.h
    include/weakform_library/weakforms_maxwell.h
    include/weakform_library/weakforms_neutronics.h
  )
    
  #
  # Build and install the library.
  #
  include(BuildAndInstallScripts)

  # Add a target for a specified version of the Hermes library,
  # and set its build properties (compile and link flags, installation directories). Due to
  # different optional libraries for H1D, H2D and H3D, this macro is dimension-specific.
  macro(BUILD_2D_LIB HERMES_LIB HERMES_COMMON_LIB BUILD_FLAGS)
    add_library(${HERMES_LIB} SHARED ${SRC} ${HEADERS})

    set(_FLAGS ${BUILD_FLAGS})

    if(NOT REPORT_WITH_LOGO)
      set(_FLAGS "${_FLAGS} -DH2D_NO_LOGO")
    endif(NOT REPORT_WITH_LOGO)
    if(NOT REPORT_TO_FILE)
      set(_FLAGS "${_FLAGS} -DH2D_REPORT_NO_FILE")
    endif(NOT REPORT_TO_FILE)
    if(REPORT_RUNTIME_CONTROL)
      set(_FLAGS "${_FLAGS} -DH2D_REPORT_RUNTIME_CONTROL")
    endif(REPORT_RUNTIME_CONTROL)
    if(REPORT_ALL)
      set(_FLAGS "${_FLAGS} -DHERMES_REPORT_ALL")
    else(REPORT_ALL)
      if(REPORT_WARNING)
        set(_FLAGS "${_FLAGS} -DHERMES_REPORT_WARNING")
      endif(REPORT_WARNING)
      if(REPORT_INTR_WARNING)
        set(_FLAGS "${_FLAGS} -DHERMES_REPORT_INTR_WARNING")
      endif(REPORT_INTR_WARNING)
      if(REPORT_INFO)
        set(_FLAGS "${_FLAGS} -DHERMES_REPORT_INFO")
      endif(REPORT_INFO)
      if(REPORT_VERBOSE)
        set(_FLAGS "${_FLAGS} -DHERMES_REPORT_VERBOSE")
      endif(REPORT_VERBOSE)
      if(REPORT_TRACE)
        set(_FLAGS "${_FLAGS} -DHERMES_REPORT_TRACE")
      endif(REPORT_TRACE)
      if(REPORT_TIME)
        set(_FLAGS "${_FLAGS} -DHERMES_REPORT_TIME")
      endif(REPORT_TIME)
    endif(REPORT_ALL)

    if(WITH_VIEWER_GUI)
      set(_FLAGS "${_FLAGS} -DENABLE_VIEWER_GUI")
    endif(WITH_VIEWER_GUI)

    set_target_properties(${HERMES_LIB} PROPERTIES COMPILE_FLAGS ${_FLAGS})

    target_link_libraries(  ${HERMES_LIB}
      ${HERMES_COMMON_LIB}
      ${GLUT_LIBRARY} ${GLEW_LIBRARY}
      ${ANTTWEAKBAR_LIBRARY}
      ${XSD_LIBRARY}
      ${XERCES_LIBRARY}
      ${LAPACK_LIBRARY}
      ${CLAPACK_LIBRARY} ${BLAS_LIBRARY}
    )
  endmacro(BUILD_2D_LIB)

  # Build and install:

  if(MSVC)
    if(H2D_DEBUG)
      set(FLAGS "${DEBUG_FLAGS}")
    elseif(H2D_DEBUG)
      set(FLAGS "${RELEASE_FLAGS}")
    endif(H2D_DEBUG)

    BUILD_2D_LIB(${HERMES2D} ${HERMES_COMMON_LIB} ${FLAGS})
    ADD_MSVC_BUILD_FLAGS(${HERMES2D} ${HERMES2D_DEBUG} ${HERMES2D_RELEASE})
    INSTALL_LIB(${HERMES2D})
  else(MSVC)
    if(H2D_DEBUG)
      BUILD_2D_LIB(${HERMES2D_DEBUG} ${HERMES_COMMON_LIB_DEBUG} ${DEBUG_FLAGS})
      INSTALL_LIB(${HERMES2D_DEBUG})
    endif(H2D_DEBUG)
    if(H2D_RELEASE)
      BUILD_2D_LIB(${HERMES2D_RELEASE} ${HERMES_COMMON_LIB_RELEASE} ${RELEASE_FLAGS})
      INSTALL_LIB(${HERMES2D_RELEASE})
    endif(H2D_RELEASE)
  endif(MSVC)

  # Install header files.

  file(GLOB INC            "include/*.h")
  file(GLOB INC_FORM       "include/form/*.h")
  file(GLOB INC_SPACE      "include/space/*.h")
  file(GLOB INC_MESH       "include/mesh/*.h")
  file(GLOB INC_XSD        "xml_schemas/*.xsd")
  file(GLOB INC_BC         "include/boundary_conditions/*.h")
  file(GLOB INC_SHAPESET   "include/shapeset/*.h")
  file(GLOB INC_FUNCTION   "include/function/*.h")
  file(GLOB INC_INTEGRALS  "include/integrals/*.h")
  file(GLOB INC_LINEARIZER "include/linearizer/*.h")
  file(GLOB INC_ADAPT      "include/adapt/*.h")
  file(GLOB INC_QUAD       "include/quadrature/*.h")
  file(GLOB INC_REF_SEL    "include/refinement_selectors/*.h")
  file(GLOB INC_VIEWS      "include/views/*.h")
  file(GLOB INC_WEAKFORM   "include/weakform/*.h")
  file(GLOB INC_WEAKFORM_LIBRARY "include/weakform_library/*.h")
  file(GLOB INC_PROJECTIONS "include/projections/*.h")

  install(FILES ${INC}           DESTINATION ${TARGET_ROOT}/include/hermes2d)
  install(FILES ${INC_SPACE}     DESTINATION ${TARGET_ROOT}/include/hermes2d/space)
  install(FILES ${INC_MESH}      DESTINATION ${TARGET_ROOT}/include/hermes2d/mesh)
  install(FILES ${XSD}           DESTINATION ${TARGET_ROOT}/include/hermes2d/xml_schemas)
  install(FILES ${INC_BC}        DESTINATION ${TARGET_ROOT}/include/hermes2d/boundary_conditions)
  install(FILES ${INC_SHAPESET}  DESTINATION ${TARGET_ROOT}/include/hermes2d/shapeset)
  install(FILES ${INC_FUNCTION}  DESTINATION ${TARGET_ROOT}/include/hermes2d/function)
  install(FILES ${INC_INTEGRALS} DESTINATION ${TARGET_ROOT}/include/hermes2d/integrals)
  install(FILES ${INC_LINEARIZER}DESTINATION ${TARGET_ROOT}/include/hermes2d/linearizer)
  install(FILES ${INC_ADAPT}     DESTINATION ${TARGET_ROOT}/include/hermes2d/adapt)
  install(FILES ${INC_QUAD}      DESTINATION ${TARGET_ROOT}/include/hermes2d/quadrature)
  install(FILES ${INC_REF_SEL}   DESTINATION ${TARGET_ROOT}/include/hermes2d/refinement_selectors)
  install(FILES ${INC_VIEWS}     DESTINATION ${TARGET_ROOT}/include/hermes2d/views)
  install(FILES ${INC_WEAKFORM}  DESTINATION ${TARGET_ROOT}/include/hermes2d/weakform)
  install(FILES ${INC_WEAKFORM_LIBRARY}  DESTINATION ${TARGET_ROOT}/include/hermes2d/weakform_library)
  install(FILES ${INC_PROJECTIONS}  DESTINATION ${TARGET_ROOT}/include/hermes2d/projections)

IF(EXISTS test_examples)
  if(H2D_WITH_TEST_EXAMPLES)
    enable_testing()
    add_subdirectory(test_examples)
  endif(H2D_WITH_TEST_EXAMPLES)
ENDIF(EXISTS test_examples)
IF(EXISTS tests)
  if(H2D_WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
  endif(H2D_WITH_TESTS)
ENDIF(EXISTS tests)