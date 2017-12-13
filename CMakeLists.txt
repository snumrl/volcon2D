cmake_minimum_required(VERSION 2.8.6)
project(volumetric)
add_compile_options(-std=gnu++11)
#add_compile_options(-mavx)
#add_compile_options(-mfma)
add_compile_options(-Wdeprecated-declarations)

SET(CMAKE_BUILD_TYPE Release CACHE STRING
	"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
#	FORCE
	)
SET( ENABLE_VMCON2D ON CACHE BOOL "Enable Vmcon2D" )
SET( ENABLE_OCTOCON2D ON CACHE BOOL "Enable Octocon2D" )
SET( ENABLE_FEM_TEST ON CACHE BOOL "Enable fem test" )


add_subdirectory( fem2D )
add_subdirectory( GUI )
include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLUT_INCLUDE_DIR})


include_directories(${CMAKE_SOURCE_DIR})

link_directories(./fem2D/)
link_directories(./GUI/)

file(GLOB vmcon2d_srcs "vmcon2D/*.h" "vmcon2D/*.cpp" "vmcon2D/iLQR/*.h" "vmcon2D/iLQR/*.cpp")
file(GLOB octocon2d_srcs "octocon2D/*.h" "octocon2D/*.cpp")
file(GLOB fem_test_srcs "fem_test/*.h" "fem_test/*.cpp")

if (ENABLE_VMCON2D)
	project(vmcon2D)
	find_package(DART REQUIRED COMPONENTS utils-urdf gui optimizer-ipopt CONFIG)
	add_compile_options(-DHAVE_CSTDDEF)
	include_directories(${DART_INCLUDE_DIRS})

	add_executable(vmcon2D ${vmcon2d_srcs})
	target_link_libraries(vmcon2D ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${DART_LIBRARIES} fem2D GUI)
endif()

if (ENABLE_OCTOCON2D)
	project(octocon2D)

	    find_library(TINYXML_LIBRARIES
      NAMES tinyxml
    )

	add_executable(octocon2D ${octocon2d_srcs})
	target_link_libraries(octocon2D ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${TINYXML_LIBRARIES} fem2D GUI)
endif()


if (ENABLE_FEM_TEST)
	project(fem_test)

	find_library(TINYXML_LIBRARIES
      NAMES tinyxml
    )

	add_executable(fem_test ${fem_test_srcs})
	target_link_libraries(fem_test ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${TINYXML_LIBRARIES} fem2D GUI)
endif()


