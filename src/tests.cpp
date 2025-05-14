#include "arith.hpp"
#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

TEST_CASE( "Simple multiplication", "[Algebraic complex numbers]" ) {
    Algebraic_Complex_Number acn = {.a = 0, .b = 1, .c = 0, .d = 0, .k = 0};
    auto result = acn * acn;

    REQUIRE(result.a == 0);
    REQUIRE(result.b == 0);
    REQUIRE(result.c == 1);
    REQUIRE(result.d == 0);
    REQUIRE(result.k == 0);
}

TEST_CASE( "Multiplication comutativity", "[Algebraic complex numbers]" ) {
    Algebraic_Complex_Number left  = {.a = 1, .b = 2, .c = 3, .d = 4, .k = 0};
    Algebraic_Complex_Number right = {.a = 0, .b = 1, .c = 2, .d = 3, .k = 1};

    auto result = left*right - right*left;
    assert(result.is_zero());
}
