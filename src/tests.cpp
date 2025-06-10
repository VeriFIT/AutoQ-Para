#include "arith.hpp"
#include "bit_set.hpp"
#include "weighted_automata.hpp"
#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

TEST_CASE( "Simple multiplication", "[Algebraic complex numbers]" ) {
    Algebraic_Complex_Number acn(0, 1, 0, 0, 0);
    auto result = acn * acn;
    auto fp_result = result.into_fixed_precision();

    REQUIRE(fp_result.a == 0);
    REQUIRE(fp_result.b == 0);
    REQUIRE(fp_result.c == 1);
    REQUIRE(fp_result.d == 0);
    REQUIRE(fp_result.k == 0);
}

TEST_CASE( "Multiplication comutativity", "[Algebraic complex numbers]" ) {
    Algebraic_Complex_Number left  (1, 2, 3, 4, 0);
    Algebraic_Complex_Number right (0, 1, 2, 3, 1);
    auto left_imm  = left*right;
    auto right_imm = right*left;
    auto result = left_imm - right_imm;

    REQUIRE(result.is_zero());
}

TEST_CASE( "Scaling during addition", "[Algebraic complex numbers]" ) {
    {
        Algebraic_Complex_Number left (1, 2, 0, 0, 0);
        Algebraic_Complex_Number right (-2, -5, 0, 0, -2);

        Algebraic_Complex_Number result = left + right;
        Fixed_Precision_ACN direct_result = result.into_fixed_precision();

        // Result should be  -9 - 6/sqrt(2) - 10i + 2i/sqrt(2)
        Fixed_Precision_ACN expected (-3, -8, 0, 0, 0);
        REQUIRE(direct_result == expected);
    }

    {
        Algebraic_Complex_Number left (1, 2, 0, 0, 0);
        Algebraic_Complex_Number right (-2, -5, 0, 0, -3);

        Algebraic_Complex_Number result = left + right;
        auto fp_result = result.into_fixed_precision();

        Fixed_Precision_ACN expected_result (11, -2, -10, 4, 0);
        REQUIRE(expected_result == fp_result);
    }

    {
        // Represents: 1 + 2/sqrt(2) + 2*i/sqrt(2) + 1*i
        Algebraic_Complex_Number left(1, 2, 1, 0, 0);

        // Represents: -7 -4/sqrt(2) -7i/sqrt(2)
        Algebraic_Complex_Number right(-2, -5, 0, 2, -1);

        auto result = left + right;
        auto fp_result = result.into_fixed_precision();

        Fixed_Precision_ACN expected_result (4, 0, -2, 2, 0);

        REQUIRE(expected_result == fp_result);
    }
}

TEST_CASE( "Conversion into direct representation", "[Algebraic complex numbers]") {

    // Represents: 1 + -2/sqrt(2) + 3i + 6i/sqrt(2)
    {
        Algebraic_Complex_Number number (1, 2, 3, 4, 0);
        Direct_ACN direct_repr = convert_acn_into_direct_repr(number);

        Direct_ACN expected_result = {1, -2, 3, 6, 0};
        REQUIRE(expected_result == direct_repr);
    }

    // Represents: -1 + 0 + 2i + 2i/sqrt(2)
    {
        Algebraic_Complex_Number number (0, 1, 2, 3, 1);
        Direct_ACN direct_repr = convert_acn_into_direct_repr(number);

        Direct_ACN expected_result = {-1, 0, 2, 2, 0};
        REQUIRE(expected_result == direct_repr);
    }

    // Represents: (-1 - 2/sqrt(2) - 4i - 2i/sqrt(2)) * 4
    {
        Algebraic_Complex_Number number (-2, -5, 2, -3, -3);
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

TEST_CASE( "MMUL", "[ACN Matrix]") {
    {
       ACN_Matrix row = row_from_ints({1, -2});
       ACN_Matrix mat = square_acn_matrix_from_ints({
           1, 2,
           3, 4
       });
       auto result = row * mat;

       auto expected = row_from_ints({-5, -6});
       REQUIRE(result == expected);
    }
    {
       ACN_Matrix row = row_from_ints({0, -2, 1, 0});
       ACN_Matrix mat = square_acn_matrix_from_ints({
           0, -2, 1, 0,
           0,  1, 0, 1,
           0,  0, 1, 2,
           0,  0, 0, 0,
       });
       auto result = row * mat;

       auto expected = row_from_ints({0, -2, 1, 0});
       REQUIRE(result == expected);
    }
}

TEST_CASE( "Intersection emptiness", "[Bit sets]") {
    {
        Bit_Set bit_set_a(10);
        bit_set_a.set_bit(0, true);
        bit_set_a.set_bit(1, true);

        Bit_Set bit_set_b(8);
        bit_set_b.set_bit(0, false);
        bit_set_b.set_bit(1, false);
        bit_set_b.set_bit(2, true);

        bool result = bit_set_a.is_intersection_empty(bit_set_b);
        REQUIRE(result);
    }

    {
        Bit_Set bit_set_a(127);
        bit_set_a.set_bit(84, true);
        bit_set_a.set_bit(85, true);

        Bit_Set bit_set_b(133);
        bit_set_b.set_bit(67, false);
        bit_set_b.set_bit(85, true);

        bool result = bit_set_a.is_intersection_empty(bit_set_b);
        REQUIRE(!result);
    }
}

TEST_CASE( "Superset checking", "[Bit sets]") {
    {
        Bit_Set bit_set_a(10);
        bit_set_a.set_bit(0, true);
        bit_set_a.set_bit(1, true);

        Bit_Set bit_set_b(65);
        bit_set_b.set_bit(0, false);
        bit_set_b.set_bit(1, false);
        bit_set_b.set_bit(2, true);

        bool result_a = bit_set_a.is_superset(bit_set_b);
        REQUIRE(!result_a);

        bool result_b = bit_set_b.is_superset(bit_set_a);
        REQUIRE(!result_b);
    }

    {
        Bit_Set bit_set_a(99);
        bit_set_a.set_bit(0,  true);
        bit_set_a.set_bit(1,  true);
        bit_set_a.set_bit(65, true);
        bit_set_a.set_bit(66, true);

        Bit_Set bit_set_b(66);
        bit_set_b.set_bit(0,  true);
        bit_set_b.set_bit(1,  true);
        bit_set_b.set_bit(65, true);

        bool result = bit_set_a.is_superset(bit_set_b);
        REQUIRE(result);
    }
}

TEST_CASE( "Zero Tests", "[Weighted automata]") {
    {
        u64 state_cnt = 3;
        ACN_Matrix initial_vec = row_from_ints({1, 0, 0});
        ACN_Matrix final_vec   = column_from_ints({0, 0, 1});
        ACN_Matrix c0_transitions = square_acn_matrix_from_ints({
            0, 1, 0,
            0, 0, 1,
            0, 0, 1,
        });

        Bit_Set c0_support (state_cnt);
        c0_support.set_all(true);

        Weighted_Automaton wa (
            initial_vec,
            final_vec,
            {c0_transitions},
            {c0_support}
        );

        bool result = is_weighted_automaton_zero(wa);

        // We can reach non-zero by, e.g., "C0, C0"
        REQUIRE(!result);
    }
    {
        u64 state_cnt = 4;
        ACN_Matrix initial_vec = row_from_ints({1, 0, 0, 0});
        ACN_Matrix final_vec   = column_from_ints({0, 0, 0, 1});
        ACN_Matrix c0_transitions = square_acn_matrix_from_ints({
            0, -2, 1, 0,
            0,  1, 0, 1,
            0,  0, 1, 2,
            0,  0, 0, 0,
        });

        Bit_Set c0_support (state_cnt);
        c0_support.set_all(true);

        Weighted_Automaton wa (
            initial_vec,
            final_vec,
            {c0_transitions},
            {c0_support}
        );

        bool result = is_weighted_automaton_zero(wa);

        // We can reach non-zero by, e.g., "C0, C0"
        REQUIRE(result);
    }
}
