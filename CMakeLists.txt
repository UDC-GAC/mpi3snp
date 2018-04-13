
# This file is part of MPI3SNP.
# Copyright (C) 2018 by Christian Ponte
#
# MPI3SNP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MPI3SNP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MPI3SNP. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.1.3)
project(MPI3SNP VERSION 1.0)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

find_package(Threads REQUIRED)
find_package(MPI REQUIRED)
find_package(CUDA)

if (NOT DEFINED TARGET_ARCH)
    if (CUDA_FOUND)
        set(TARGET_ARCH "GPU" CACHE STRING "Select target architecture")
    else ()
        set(TARGET_ARCH "CPU" CACHE STRING "Select target architecture")
    endif ()
    set_property(CACHE TARGET_ARCH PROPERTY STRINGS "CPU" "GPU")
endif ()

set(DEFAULT_BUILD_TYPE "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to ${DEFAULT_BUILD_TYPE} as none was specified.")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
            STRING "Choose the type of build" FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif ()

# The definitions of stoi, stol, etc., as well as the to_string functions, are guarded by the condition
#    if ((__cplusplus >= 201103L) && defined(_GLIBCXX_USE_C99) && !defined(_GLIBCXX_HAVE_BROKEN_VSWPRINTF))
# To prevent compatibility issues, the macro _GLIBCXX_USE_C99 will be defined manually
set(CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_C99")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")
if (${TARGET_ARCH} STREQUAL "GPU")
    set(CUDA_NVCC_FLAGS "-arch=sm_30 -std=c++11 -D_GLIBCXX_USE_C99")
    set(CUDA_NVCC_FLAGS_DEBUG "-g -O0")
    set(CUDA_NVCC_FLAGS_RELEASE "-DNDEBUG -O3")
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
endif ()

# Third party libraries (common to GPU and CPU builds)
add_subdirectory(thirdparty)
add_subdirectory(src)
