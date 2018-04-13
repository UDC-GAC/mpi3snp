
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

list(APPEND COMMON_SOURCE_FILES
        main.cpp Arg_parser.cpp Arg_parser.h IOMpi.cpp IOMpi.h Node_information.cpp Node_information.h
        Cpu_node_information.cpp Cpu_node_information.h Distributor.cpp Distributor.h Individual.h SNP.h
        Dataset.cpp Dataset.h Search.cpp Search.h Statistics.cpp Statistics.h)

configure_file(${PROJECT_SOURCE_DIR}/src/Definitions.h.in ${PROJECT_BINARY_DIR}/src/Definitions.h)

if (${TARGET_ARCH} STREQUAL "CPU")
    include_directories(${MPI_CXX_INCLUDE_PATH} ${PROJECT_BINARY_DIR}/src/)
    add_executable(MPI3SNP ${COMMON_SOURCE_FILES} cpu/CPUEngine.cpp cpu/CPUEngine.h cpu/EntropySearch.cpp cpu/EntropySearch.h)
    target_link_libraries(MPI3SNP ${MPI_CXX_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} args)
else ()
    include_directories(${MPI_CXX_INCLUDE_PATH} ${PROJECT_BINARY_DIR}/src/)
    add_library(obj OBJECT ${COMMON_SOURCE_FILES})
    target_include_directories(obj PRIVATE $<TARGET_PROPERTY:args,INTERFACE_INCLUDE_DIRECTORIES>)
    cuda_add_executable(MPI3SNP $<TARGET_OBJECTS:obj> gpu/Gpu_node_information.cpp gpu/Gpu_node_information.h gpu/GPUEngine.cpp gpu/GPUEngine.h gpu/EntropySearch.cu gpu/EntropySearch.h)
    target_link_libraries(MPI3SNP ${MPI_CXX_LIBRARIES})
endif ()
