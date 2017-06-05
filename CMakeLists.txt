cmake_minimum_required(VERSION 3.0)
project(MPI3SNP)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_FLAGS "-O2")

find_package(MPI REQUIRED)
find_package(CUDA REQUIRED)

# Falta cambiar archivos
#set(MPI3SNP_c++_files src/main.cpp .cpp)
#set(MPI3SNP_cuda_files .cu)

#include_directories(${MPI_CXX_INCLUDE_PATH})
#add_library(cpp_objs OBJECT ${MPI3SNP_c++_files})
#cuda_add_executable(MPI3SNP $<TARGET_OBJECTS:cpp_objs> ${MPI3SNP_cuda_files})
#target_link_libraries(MPI3SNP ${MPI_CXX_LIBRARIES})

file(GLOB SOURCE_FILES src/*.cpp src/*.cu)
cuda_add_executable(MPI3SNP ${SOURCE_FILES})