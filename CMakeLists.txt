cmake_minimum_required(VERSION 2.8)

if(NOT CMAKE_BUILD_TYPE)
    message("** CMAKE_BUILD_TYPE not passed. Assuming Debug")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

#add the customised package searches to the module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

#find_package (Eigen3 3.0.0 REQUIRED)
#include_directories(${EIGEN3_INCLUDE_DIR})

enable_testing()
add_subdirectory(tests)
