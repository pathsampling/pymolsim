cmake_minimum_required(VERSION 2.8.12)

project(pymolsim)
#------------------------------------------------------

# Set source directory
#------------------------------------------------------
# main directory for steinhardt cpp modules
set(SOURCE_DIR "src/pymolsim")


#Headers
#-------------------------------------------------------
# Tell CMake that headers are also in SOURCE_DIR
include_directories(${SOURCE_DIR})
set(SYSSOURCES "${SOURCE_DIR}/analyse.cpp"
			   "${SOURCE_DIR}/forces.cpp"
			   "${SOURCE_DIR}/init.cpp"
			   "${SOURCE_DIR}/md.cpp"
			   "${SOURCE_DIR}/ran3.cpp"
			   "${SOURCE_DIR}/output.cpp"
			   )

#Now add subdir for the modules
#--------------------------------------------------------
add_subdirectory(lib/pybind11)
#message(STATUS "include dir: ${CMAKE_BUILD_TYPE}")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O2 -DNDEBUG")
#set(CXX_STANDARD 11)

#Finally create the package
#--------------------------------------------------------
pybind11_add_module(cmolsim ${SYSSOURCES} "${SOURCE_DIR}/bindings.cpp" )
