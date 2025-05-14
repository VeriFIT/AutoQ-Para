cmake_minimum_required(VERSION 3.30)

set(CMAKE_PREFIX_PATH
    ${CMAKE_PREFIX_PATH}
    /usr/lib64/cmake/Catch2/)

find_package(Catch2 REQUIRED)

project(
    AlgebraicComplexNumbers
    VERSION 0.1
    LANGUAGES C CXX
)

add_library(acn STATIC arith.cpp)
add_executable(runner runner.cpp)

target_link_libraries(runner PRIVATE acn)

add_executable(tests tests.cpp)
