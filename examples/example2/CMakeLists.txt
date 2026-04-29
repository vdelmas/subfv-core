cmake_minimum_required(VERSION 3.20)

project(example2 Fortran)

# Si la lib n'est pas déjà dans le projet → on la cherche
if(NOT TARGET subfvcore AND NOT TARGET subfvcore::subfvcore)
  find_package(subfvcore REQUIRED)
endif()

add_executable(example2 example2_module.F90 example2.F90)

# Linking compatible dev + install
if(TARGET subfvcore::subfvcore)
  target_link_libraries(example2 subfvcore::subfvcore)
else()
  target_link_libraries(example2 subfvcore)
endif()
