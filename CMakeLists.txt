cmake_minimum_required(VERSION 3.0)
project(MPI3SNP VERSION 1.0)
set(CMAKE_CXX_STANDARD 14)
set(PROJECT_NAME MPI3SNP)
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
#set(CMAKE_CXX_FLAGS "-O2")
add_subdirectory(thirdparty/args)

find_package(MPI REQUIRED)
find_package(CUDA REQUIRED)

configure_file(${PROJECT_SOURCE_DIR}/src/Definitions.h.in ${PROJECT_BINARY_DIR}/src/Definitions.h)
set(MPI3SNP_c++_files src/main.cpp src/GPUSearchMI.cpp src/LineReader.cpp src/Utils.cpp
        src/IOMpi.cpp src/Statistics.cpp)
set(MPI3SNP_cuda_files src/EntropySearch.cu src/GPUEngine.cpp src/GPUMacros.cpp src/GPUSNPDistributor.cpp)
#src/GPUSNPDistributorStatic.cpp

include_directories(${MPI_CXX_INCLUDE_PATH} ${ARGS_SOURCE_DIR} ${PROJECT_BINARY_DIR}/src/)
add_library(cpp_objs OBJECT ${MPI3SNP_c++_files})
cuda_add_executable(MPI3SNP $<TARGET_OBJECTS:cpp_objs> ${MPI3SNP_cuda_files})
target_link_libraries(MPI3SNP ${MPI_CXX_LIBRARIES})