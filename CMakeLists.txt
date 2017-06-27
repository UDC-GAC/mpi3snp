cmake_minimum_required(VERSION 3.0)
project(MPI3SNP VERSION 1.0)
set(CMAKE_CXX_STANDARD 14)
set(PROJECT_NAME MPI3SNP)
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
set(CMAKE_CXX_FLAGS "-g")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -g")
add_subdirectory(thirdparty/args)


# GPU ------------------
find_package(MPI REQUIRED)
find_package(CUDA REQUIRED)

configure_file(${PROJECT_SOURCE_DIR}/src/Definitions.h.in ${PROJECT_BINARY_DIR}/src/Definitions.h)
set(MPI3SNP_c++_files src/main.cpp src/Dataset.cpp src/Search.cpp src/IOMpi.cpp src/Statistics.cpp)
set(MPI3SNP_cuda_files src/Distributor.cpp src/gpu/EntropySearch.cu src/gpu/GPUEngine.cpp)

include_directories(${MPI_CXX_INCLUDE_PATH} ${ARGS_SOURCE_DIR} ${PROJECT_BINARY_DIR}/src/)
add_library(cpp_objs OBJECT ${MPI3SNP_c++_files})
cuda_add_executable(MPI3SNP $<TARGET_OBJECTS:cpp_objs> ${MPI3SNP_cuda_files})
target_link_libraries(MPI3SNP ${MPI_CXX_LIBRARIES})