cmake_minimum_required(VERSION 3.0)
project(MPI3SNP VERSION 1.0)
set(CMAKE_CXX_STANDARD 14)
set(PROJECT_NAME MPI3SNP)
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)

find_package(Threads REQUIRED)
find_package(MPI REQUIRED)
find_package(CUDA)

# CXX compiler flags (common to GPU and CPU builds)
# The definitions of stoi, stol, etc., as well as the to_string functions, are guarded by the condition
#    if ((__cplusplus >= 201103L) && defined(_GLIBCXX_USE_C99) && !defined(_GLIBCXX_HAVE_BROKEN_VSWPRINTF))
# To prevent compatibility issues, the macro _GLIBCXX_USE_C99 will be defined manually
list(APPEND CMAKE_CXX_FLAGS "-std=c++14 -D_GLIBCXX_USE_C99 -O3")
# Third party libraries (common to GPU and CPU builds)
add_subdirectory(thirdparty/args)

if (CUDA_FOUND AND NOT IGNORE_GPU)
    # GPU build
    list(APPEND CUDA_NVCC_FLAGS "-arch=sm_30 -std=c++11 -D_GLIBCXX_USE_C99 -O3")
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)

    set(MPI3SNP_USE_GPU ON)
    set(MPI3SNP_cuda_files src/gpu/Gpu_node_information.cpp src/gpu/EntropySearch.cu src/gpu/GPUEngine.cpp)
    set(MPI3SNP_c++_files src/main.cpp src/Arg_parser.cpp src/IOMpi.cpp src/Node_information.cpp src/Cpu_node_information.cpp
            src/Distributor.cpp src/Individual.h src/SNP.h src/Dataset.cpp src/Search.cpp src/Statistics.cpp)

    configure_file(${PROJECT_SOURCE_DIR}/src/Definitions.h.in ${PROJECT_BINARY_DIR}/src/Definitions.h)
    include_directories(${MPI_CXX_INCLUDE_PATH} ${ARGS_SOURCE_DIR} ${PROJECT_BINARY_DIR}/src/)
    add_library(cpp_objs OBJECT ${MPI3SNP_c++_files})
    cuda_add_executable(MPI3SNP $<TARGET_OBJECTS:cpp_objs> ${MPI3SNP_cuda_files})
    target_link_libraries(MPI3SNP ${MPI_CXX_LIBRARIES})
else (CUDA_FOUND AND NOT IGNORE_GPU)
    # CPU build
    set(SOURCE_FILES src/main.cpp src/Arg_parser.cpp src/IOMpi.cpp src/Node_information.cpp src/Cpu_node_information.cpp
            src/Statistics.cpp src/Distributor.cpp src/Dataset.cpp src/Search.cpp src/cpu/CPUEngine.cpp src/cpu/EntropySearch.cpp)

    configure_file(${PROJECT_SOURCE_DIR}/src/Definitions.h.in ${PROJECT_BINARY_DIR}/src/Definitions.h)
    include_directories(${MPI_CXX_INCLUDE_PATH} ${ARGS_SOURCE_DIR} ${PROJECT_BINARY_DIR}/src/)
    add_executable(MPI3SNP ${SOURCE_FILES})
    target_link_libraries(MPI3SNP ${MPI_CXX_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif (CUDA_FOUND AND NOT IGNORE_GPU)
