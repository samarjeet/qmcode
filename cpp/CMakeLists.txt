cmake_minimum_required (VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)
project(psi4_crawford)


set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ") 
include_directories("${CMAKE_SOURCE_DIR}/include")
include_directories("${CMAKE_SOURCE_DIR}/include/Eigen")
include_directories("${CMAKE_SOURCE_DIR}/src")

FILE(GLOB SRC "${CMAKE_SOURCE_DIR}/src/*.cpp")
set(PROJ1_SRC "${CMAKE_SOURCE_DIR}/test/proj1.cpp")
#add_executable(proj1 ${SRC} ${PROJ1_SRC})
set(PROJ2_SRC "${CMAKE_SOURCE_DIR}/test/proj2.cpp")
#add_executable(proj2 ${SRC} ${PROJ2_SRC})
set(PROJ3_SRC "${CMAKE_SOURCE_DIR}/test/proj3.cpp")
add_executable(proj3 ${SRC} ${PROJ3_SRC})


