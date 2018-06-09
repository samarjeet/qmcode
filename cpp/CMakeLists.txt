cmake_minimum_required (VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)
project(psi4_crawford)


set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ") 
include_directories("${CMAKE_SOURCE_DIR}/include")
include_directories("${CMAKE_SOURCE_DIR}/src")

FILE(GLOB SRC "${CMAKE_SOURCE_DIR}/src/*.cpp")
add_executable(psitest ${SRC})


