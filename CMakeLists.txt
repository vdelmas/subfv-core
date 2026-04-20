cmake_minimum_required(VERSION 3.20)

project(subfv-core Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_FLAGS "-O3 -fimplicit-none")

find_package(MPI REQUIRED)

set(CORE_SRC
  src/precision_module.F90
  src/mpi_module.F90
  src/mesh_module.F90
  src/mesh_reading_module.F90
  src/mesh_geometry_module.F90
  src/mesh_connectivity_module.F90
)

add_library(subfvcore ${CORE_SRC})

set_target_properties(subfvcore PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
)

target_include_directories(subfvcore PUBLIC
  ${CMAKE_BINARY_DIR}/mod
)

target_link_libraries(subfvcore PUBLIC MPI::MPI_Fortran)
