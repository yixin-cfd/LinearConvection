cmake_minimum_required(VERSION 3.8)

enable_language(Fortran)

project(LinearConvection Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fallow-argument-mismatch -ffree-line-length-none")

file(GLOB src_path "*.f90" )

add_executable(LC ${src_path})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
