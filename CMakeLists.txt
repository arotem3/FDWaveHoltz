cmake_minimum_required(VERSION 3.16)
project(fdwh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
set(COMMON_FLAGS "-Wall -fopenmp")
set(RELEASE_FLAGS "-O3")
set(DEBUG_FLAGS "-g -DWH_DEBUG")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${COMMON_FLAGS} ${DEBUG_FLAGS}")
else()
    set(CMAKE_CXX_FLAGS "${COMMON_FLAGS} ${RELEASE_FLAGS}")
endif()

# Directories
set(SRC_DIR ${CMAKE_SOURCE_DIR}/source)
set(BUILD_DIR ${CMAKE_BINARY_DIR}/build)
set(DRIVERS_DIR ${CMAKE_SOURCE_DIR}/drivers)
set(PAPER_DIR ${CMAKE_SOURCE_DIR}/paper)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

include_directories(${INCLUDE_DIR})
file(GLOB SRC_FILES ${SRC_DIR}/*.cpp)

# Find LAPACK, BLAS, and OpenMP
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)

# add linsolve
add_subdirectory(${CMAKE_SOURCE_DIR}/LinSolve)

# Library
add_library(fdwh STATIC ${SRC_FILES})
target_link_libraries(fdwh PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(fdwh PUBLIC linsolve)

# Drivers
file(GLOB DRIVERS_SRC ${DRIVERS_DIR}/*.cpp)
foreach(driver ${DRIVERS_SRC})
    get_filename_component(driver_name ${driver} NAME_WE)
    add_executable(${driver_name} ${driver})
    target_link_libraries(${driver_name} fdwh)
endforeach()

# Paper executables
file(GLOB PAPER_SRC ${PAPER_DIR}/*.cpp)
foreach(paper ${PAPER_SRC})
    get_filename_component(paper_name ${paper} NAME_WE)
    add_executable(${paper_name} ${paper})
    target_link_libraries(${paper_name} fdwh LAPACK::LAPACK BLAS::BLAS)
endforeach()

# solution subdirectory
# Solution subdirectory
set(SOLUTION_DIR ${BUILD_DIR}/solution)
file(MAKE_DIRECTORY ${SOLUTION_DIR})