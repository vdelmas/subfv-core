cmake_minimum_required(VERSION 3.20)

project(subfvcore Fortran)

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library(subfvcore
  src/utils/subfv_precision_module.F90

  src/mpi/subfv_mpi_module.F90

  src/mesh/subfv_mesh_module.F90
  src/mesh/subfv_mesh_reading_module.F90
  src/mesh/subfv_mesh_geometry_module.F90
  src/mesh/subfv_mesh_connectivity_module.F90

  src/linear_algebra/subfv_sparse_linear_module.F90
  src/linear_algebra/subfv_sparse_csr_linear_module.F90
  src/linear_algebra/subfv_sparse_bcsr_linear_module.F90
  src/linear_algebra/subfv_linear_solver_module.F90

  src/io/subfv_io_module.F90
)

set_target_properties(subfvcore PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
)

target_link_libraries(subfvcore
  PUBLIC
  MPI::MPI_Fortran
  LAPACK::LAPACK
  BLAS::BLAS
)

target_include_directories(subfvcore
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/mod>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

include(GNUInstallDirs)

install(TARGETS subfvcore
  EXPORT subfvcore-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Installer les fichiers modules Fortran (.mod)
install(DIRECTORY ${CMAKE_BINARY_DIR}/mod/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT subfvcore-targets
  FILE subfvcoreTargets.cmake
  NAMESPACE subfvcore::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/subfvcore
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/subfvcoreConfigVersion.cmake"
  VERSION 1.0.0
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/subfvcoreConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/subfvcoreConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/subfvcore
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/subfvcoreConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/subfvcoreConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/subfvcore
)

# =========================
# Examples
# =========================

option(BUILD_EXAMPLES "Build examples" ON)

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# Installer les exemples
install(DIRECTORY examples/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/subfvcore/examples
)
