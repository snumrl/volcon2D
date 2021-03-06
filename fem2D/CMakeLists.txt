cmake_minimum_required(VERSION 2.8.6)
project(fem2D)
add_compile_options(-std=gnu++11)

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

file(GLOB srcs "*.h" "*.cpp")
file(GLOB constraint_srcs "Constraint/*.h" "Constraint/*.cpp")
file(GLOB mesh_srcs "Mesh/*.h" "Mesh/*.cpp")

add_library(fem2D ${srcs} ${constraint_srcs} ${mesh_srcs})
						   
