cmake_minimum_required(VERSION 3.30)

set(CMAKE_PREFIX_PATH
    ${CMAKE_PREFIX_PATH}
    /usr/lib64/cmake/Catch2/)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Catch2)
find_library(GMP_LIBRARY gmp)

project(
    SWTAs
    VERSION 0.1
    LANGUAGES C CXX
)

add_library(acn STATIC arith.cpp)
add_executable(experiments experiments.cpp arith.cpp weighted_automata.cpp swta.cpp predefined_automata.cpp bit_set.cpp nfa.cpp swta_builders.cpp quantum_program.cpp)

target_link_libraries(experiments PRIVATE acn gmp)

add_executable(tests tests.cpp arith.cpp weighted_automata.cpp swta.cpp predefined_automata.cpp bit_set.cpp nfa.cpp swta_builders.cpp quantum_program.cpp)
target_link_libraries(tests PRIVATE acn gmp)
