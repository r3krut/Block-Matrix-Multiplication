cmake_minimum_required(VERSION 2.8)
project(MultiplicationBlockMatrix)

include_directories(headers)
file(GLOB SOURCES "src/*.cpp")

set (CMAKE_CXX_STANDARD 11)
#clang compiler
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_C_COMPILER clang)

set(CMAKE_BUILD_TYPE Release)
#set(CXX_FLAGS -O3)
set(CMAKE_CXX_FLAGS_MYREL "-O3")

find_package(OpenMP)
if (OPENMP_FOUND)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(${PROJECT_NAME} ${SOURCES})
