cmake_minimum_required(VERSION 3.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
add_compile_options(-Wall -Wno-unused-result -O3 -g)
#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)
project(pivot)

aux_source_directory(. SRC)

add_executable(pivot ${SRC})

target_link_libraries(pivot m)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
target_link_libraries(pivot MPI::MPI_CXX)
target_link_libraries(pivot OpenMP::OpenMP_CXX)

