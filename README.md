# subfv-core

Subface based Finite Volume core library written in Fortran built with CMake.

## Requirements

- Fortran compiler gfortran ifort etc
- CMake greater or equal 3.20
- MPI implementation OpenMPI or MPICH for parallel features
- Gmsh optional
- Paraview optional

## Configure Build

mkdir build
cd build
cmake ..

If MPI is not automatically detected

cmake .. -DMPI_Fortran_COMPILER=mpif90

You can also set the compiler explicitly

cmake .. -DCMAKE_Fortran_COMPILER=gfortran

CMake will generate
- Makefiles
- Fortran module directory
- Build configuration for library and examples

## Compile

make -j

This builds
- libsubfvcore.a
- Fortran module files mod files
- example executables

## Install

Recommended local installation no sudo required

cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local
make -j
make install

This installs
- libsubfvcore.a
- mod files
- subfvcoreConfig.cmake
- subfvcoreConfigVersion.cmake

## Environment setup

Add this to your shell configuration file bashrc or zshrc

export PATH=$HOME/.local/bin:$PATH
export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH
export CMAKE_PREFIX_PATH=$HOME/.local:$CMAKE_PREFIX_PATH

## Using the library in your project

In CMake

find_package subfvcore REQUIRED

If not found automatically

set(subfvcore_DIR $HOME/.local/lib/cmake/subfvcore)
    find_package subfvcore REQUIRED

## Running Examples

    Example 1 serial

    cd examples/example1/
    gmsh -3 example1.geo
    cd build/
    cmake ..
    make -j
    ./example1

    Example 2 Serial

    cd examples/example2/
    gmsh -3 example2.geo
    cd build/
    cmake ..
    make -j
    ./example2

    Example 2 MPI

    cd examples/example2/
    gmsh -3 example2.geo
    ../../tools/scripts/partition_mesh.sh example2.msh 4
    cd build/
    cmake ..
    make -j
    mpirun -np 4 ./example2

## Project Structure

    src/ core library source code
    src/io/ input output utilities
    src/mesh mesh handling geometry connectivity reading
    src/mpi MPI abstraction layer
    src/linear algebra sparse solvers and algebra
    src/utils precision and helper tools

    examples/ examples of usage
    example1 serial case to check the mesh reading
    example2 simple Rusanov Euler code for Mach 20 half cylinder

    tools utilities 
    tools/scripts scripts 

## Module Overview

    Mesh modules
    subfv mesh module
    subfv mesh geometry module
    subfv mesh connectivity module
    subfv mesh reading module

    Linear algebra modules
    subfv sparse linear module
    subfv sparse csr linear module
    subfv sparse bcsr linear module
    subfv linear solver module

    MPI
    subfv mpi module

    Utilities
    subfv io module
    subfv precision module

    ---

## Cleaning Build

    rm -rf build
    mkdir build
    cd build
    cmake ..
    make

    Remove artifacts

    find . name .o delete
    find . name .mod delete

    ---

## Visualization

Open results with Paraview, one vtu file is generated for each MPI subdomain, always open the header pvtu file. One may need to use the D3 filter in Paraview (Filters -> Search) to make sure the iso-surfaces and contours are plotted correctly across the subdomains boundaries.

## Author

Developed by
- Vincent Delmas, University of Bordeaux, Institut de Mathematiques de Bordeaux (IMB)

