cmake_minimum_required(VERSION 3.20)

project(subfvcore Fortran)

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library(subfvcore
  src/utils/precision_module.F90

  src/mpi/mpi_module.F90

  src/mesh/mesh_module.F90
  src/mesh/mesh_reading_module.F90
  src/mesh/mesh_geometry_module.F90
  src/mesh/mesh_connectivity_module.F90

  src/linear_algebra/sparse_linear_module.F90
  src/linear_algebra/sparse_csr_linear_module.F90
  src/linear_algebra/sparse_bcsr_linear_module.F90
  src/linear_algebra/linear_solver_module.F90
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
