## CMakeLists.txt - Main CMake file that defines how cmake should process and generate the necessary build files
 # written by C. D. Degawa, last modified in 2023 - see License.htm for legal notices
 #
 # The copyright in this software is being made available under the exhale Copyright License
 # and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 # party rights, including patent rights. No such rights are granted under this License.
 #
 # Copyright (c) 2018-2024 Christian R. Helmrich, project ecodis. All rights reserved.
 ##

cmake_minimum_required(VERSION 3.5) # Default version of cmake on Ubuntu 16.04

if("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    message(FATAL_ERROR "Building in the source tree is not supported.\n"
                        "Please re-run cmake from a build folder.")
endif()


project(exhale VERSION 1.2.1 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release
        CACHE
        STRING "Build type: Debug, Release, RelWithDebInfo or MinSizeRel"
        FORCE)
endif()

set(CMAKE_THREAD_PREFER_PTHREAD true)
find_package(Threads) # For Threads::Threads
include(GNUInstallDirs)

set(CMAKE_POSITION_INDEPENDENT_CODE true)
set(CXX_STANDARD 11)
if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "/WX /D_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")
else()
    set(CMAKE_CXX_FLAGS "-Wall -Wshadow ${CMAKE_CXX_FLAGS}")
    if(WIN32 AND MINGW)
        set(CMAKE_CXX_FLAGS "-municode ${CMAKE_CXX_FLAGS}")
    endif()
    set(CMAKE_CXX_FLAGS_RELEASE "-Wuninitialized ${CMAKE_CXX_FLAGS_RELEASE}")
    set(CMAKE_CXX_FLAGS_DEBUG "-Werror -D_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    add_definitions(-D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64)
endif()


add_subdirectory(src/lib)
add_subdirectory(src/app)
