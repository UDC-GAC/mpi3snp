cmake_minimum_required(VERSION 3.0)
project(MPI3SNP)

file(GLOB SOURCE_FILES src/*.cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER, "mpicxx")
set(CMAKE_CXX_FLAGS "-O2")

find_package(Threads REQUIRED)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(THREADS_PREFER_PTHREAD_FLAG ON)

add_executable(MPI3SNP ${SOURCE_FILES})
target_link_libraries(MPI3SNP pthread)
target_link_libraries(MPI3SNP ${MPI_LIBRARIES})
