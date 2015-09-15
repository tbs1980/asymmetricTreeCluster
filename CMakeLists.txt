cmake_minimum_required(VERSION 2.8)

if(NOT CMAKE_BUILD_TYPE)
    message("** CMAKE_BUILD_TYPE not passed. Assuming Debug")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

enable_testing()
add_subdirectory(tests)
