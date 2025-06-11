cmake_minimum_required(VERSION 3.30)

set(CMAKE_PREFIX_PATH
    ${CMAKE_PREFIX_PATH}
    /usr/lib64/cmake/Catch2/)

find_package(Catch2 REQUIRED)
find_library(GMP_LIBRARY gmp)

project(
    SWTAs
    VERSION 0.1
    LANGUAGES C CXX
)

add_library(acn STATIC arith.cpp)
add_executable(runner runner.cpp)

target_link_libraries(runner PRIVATE acn gmp)

add_executable(tests tests.cpp arith.cpp weighted_automata.cpp swta.cpp predefined_wtts.cpp)
target_link_libraries(tests PRIVATE acn gmp)
