cmake_minimum_required(VERSION 3.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
add_compile_options(-Wall -Wno-unused-result -O3 -g -mavx2)
#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)
project(pivot)


add_executable(pivot pivot.cc)

target_link_libraries(pivot m)

find_package(OpenMP REQUIRED)
target_link_libraries(pivot OpenMP::OpenMP_CXX)

