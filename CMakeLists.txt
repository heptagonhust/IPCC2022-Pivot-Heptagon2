cmake_minimum_required(VERSION 3.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_compile_options(-Wall -Werror -Wno-unused-result -O3 -g)
# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)
project(pivot)

aux_source_directory(. SRC)

add_executable(pivot ${SRC})
find_package(OpenMP REQUIRED)

target_link_libraries(pivot m)
target_link_libraries(pivot OpenMP::OpenMP_CXX)
