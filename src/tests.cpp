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

TEST_CASE( "Scaling during addition", "[Algebraic complex numbers]" ) {
    {
        // Represents: 1 + 2/sqrt(2) + 2*i/sqrt(2)
        Algebraic_Complex_Number left  = {.a = 1, .b = 2, .c = 0, .d = 0, .k = 0};

        // Represents: -10 -8/sqrt(2) - 10*i
        Algebraic_Complex_Number right = {.a = -2, .b = -5, .c = 0, .d = 0, .k = -2};

        Algebraic_Complex_Number result = left + right;
        Direct_ACN direct_result = convert_acn_into_direct_repr(result);

        // Result should be  -9 - 6/sqrt(2) - 10i + 2i/sqrt(2)
        Direct_ACN expected = {-3, -8, 0, -8, 0};
        REQUIRE(direct_result == expected);
    }

    {
        // Represents: 1 + 2/sqrt(2) + 2*i/sqrt(2)
        Algebraic_Complex_Number left  = {.a = 1, .b = 2, .c = 0, .d = 0, .k = 0};

        Algebraic_Complex_Number right = {.a = -2, .b = -5, .c = 0, .d = 0, .k = -3};

        Algebraic_Complex_Number result = left + right;

        REQUIRE(result.a == 11);
        REQUIRE(result.b == -2);
        REQUIRE(result.c == -10);
        REQUIRE(result.d == 4);
        REQUIRE(result.k == 0);
    }

    {
        // Represents: 1 + 2/sqrt(2) + 2*i/sqrt(2) + 1*i
        Algebraic_Complex_Number left  = {.a = 1, .b = 2, .c = 1, .d = 0, .k = 0};

        // Represents: -7 -4/sqrt(2) -7i/sqrt(2)
        Algebraic_Complex_Number right = {.a = -2, .b = -5, .c = 0, .d = 2, .k = -1};

        Algebraic_Complex_Number result = left + right;

        REQUIRE(result.a == 4);
        REQUIRE(result.b == 0);
        REQUIRE(result.c == -2);
        REQUIRE(result.d == 2);
        REQUIRE(result.k == 0);
    }
}

TEST_CASE( "Conversion into direct representation", "[Algebraic complex numbers]") {

    // Represents: 1 + -2/sqrt(2) + 3i + 6i/sqrt(2)
    {
        Algebraic_Complex_Number number = {.a = 1, .b = 2, .c = 3, .d = 4, .k = 0};
        Direct_ACN direct_repr = convert_acn_into_direct_repr(number);

        Direct_ACN expected_result = {1, -2, 3, 6, 0};
        REQUIRE(expected_result == direct_repr);
    }

    // Represents: -1 + 0 + 2i + 2i/sqrt(2)
    {
        Algebraic_Complex_Number number = {.a = 0, .b = 1, .c = 2, .d = 3, .k = 1};
        Direct_ACN direct_repr = convert_acn_into_direct_repr(number);

        Direct_ACN expected_result = {-1, 0, 2, 2, 0};
        REQUIRE(expected_result == direct_repr);
    }

    // Represents: (-1 - 2/sqrt(2) - 4i - 2i/sqrt(2)) * 4
    {
        Algebraic_Complex_Number number = {.a = -2, .b = -5, .c = 2, .d = -3, .k = -3};
        Direct_ACN direct_repr = convert_acn_into_direct_repr(number);

        Direct_ACN expected_result = {-1, -2, -4, 2, -2};
        REQUIRE(expected_result == direct_repr);
    }
}

TEST_CASE( "Add row to a row-echelon-form matrix", "[ACN Matrix]") {
    {
        ACN_Matrix matrix = square_acn_matrix_from_ints({
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
        });
        ACN_Matrix row = row_from_ints({1, 0, 0});

        s64 row_slot_idx        = add_row_to_row_echelon_matrix(matrix, row);
        s64 subsequent_slot_idx = add_row_to_row_echelon_matrix(matrix, row);

        REQUIRE(row_slot_idx == 0);
        REQUIRE(subsequent_slot_idx == -1);
    }

    {
        ACN_Matrix matrix = square_acn_matrix_from_ints({
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
        });
        ACN_Matrix row = row_from_ints({1, 2, 3});

        s64 row_slot_idx = add_row_to_row_echelon_matrix(matrix, row);

        REQUIRE(row_slot_idx == -1);
    }

    {
        ACN_Matrix matrix = square_acn_matrix_from_ints({
            1, 0, 0,
            0, 3, 0,
            0, 0, 8,
        });
        ACN_Matrix row = row_from_ints({1, 2, 0});

        s64 row_slot_idx = add_row_to_row_echelon_matrix(matrix, row);

        REQUIRE(row_slot_idx == -1);
    }
}
