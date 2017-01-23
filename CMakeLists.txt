cmake_minimum_required(VERSION 3.6)
project(MPI3SNP)

file(GLOB SOURCE_FILES src/*.cpp)

# set(CMAKE_C_COMPILER "")
# set(GCC_COVERAGE_COMPILE_FLAGS "-O3")
# set(GCC_COVERAGE_LINK_FLAGS    "-lgcov")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER, "mpic++")
set(CMAKE_CXX_FLAGS "-O2")

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(THREADS_PREFER_PTHREAD_FLAG ON)
if(OPENMP_FOUND)
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(MPI3SNP ${SOURCE_FILES})
target_link_libraries(MPI3SNP Threads::Threads)
target_link_libraries(MPI3SNP ${MPI_LIBRARIES})
