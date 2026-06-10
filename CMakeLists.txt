cmake_minimum_required(VERSION 3.9.0)

# set(CMAKE_CXX_COMPILER "/public1/soft/gcc/9.3.0/bin/g++")
# set(CMAKE_C_COMPILER "/public1/soft/gcc/9.3.0/bin/gcc")
set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
# set(CMAKE_CXX_COMPILER "/usr/bin/g++")
# set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-O3")

project(sphsol VERSION 1.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Eigen path on MacOS
include_directories(/opt/homebrew/opt/eigen/include)

find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

aux_source_directory(${PROJECT_SOURCE_DIR}/src DIR_SRCS)
add_executable(${PROJECT_NAME} ${DIR_SRCS})

add_definitions(-w)
