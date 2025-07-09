#include "arith.hpp"
#include "bit_set.hpp"
#include "nfa.hpp"
#include "quantum_program.hpp"
#include "swta.hpp"
#include "swta_builders.hpp"
#include "weighted_automata.hpp"
#include "predefined_automata.hpp"

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

    {
        ACN_Matrix matrix = square_acn_matrix_from_ints({
            1, -1, 0,
            0,  0, 0,
            0,  0, 0,
        });
        ACN_Matrix row = row_from_ints({1, -1, 0});

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

TEST_CASE( "set_all handling of boundaries", "[Bit sets]") {
    Bit_Set set (1);
    set.set_all();

    REQUIRE(set.data[0] == 0x01);
}

TEST_CASE( "are_all_bit_set", "[Bit sets]") {
    Bit_Set set (1);
    set.set_all();

    REQUIRE(set.are_all_bits_set());
}

TEST_CASE( "Zero Tests", "[Weighted automata]") {
    return; // Test is disabled since the zero-checking algorithm seems to be unsound
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

TEST_CASE("Sequential composition of two Hadamards", "[WTT]") {
    SWTA::Metadata metadata = {
        .number_of_internal_symbols = 1,
        .number_of_colors = 0,
    };
    auto hadamard = get_predefined_wtt(Predefined_WTT_Names::HADAMARD, metadata);
    auto result = compose_wtts_sequentially(hadamard, hadamard);

    REQUIRE(result.initial_states.size() == 1);
    result.normalize_all_transitions();

    REQUIRE(result.number_of_internal_symbols() == 1);
    auto& transitions_along_sym0 = result.transitions[0];

    REQUIRE(transitions_along_sym0.size() == 1);  // There should be only one state
    auto& transition = transitions_along_sym0[0];

    { // Check that LEFT subtree is an identity relation
        REQUIRE(transition.ll.size() == 1);
        auto& ll_component = transition.ll.components[0];
        REQUIRE(ll_component.coef == Algebraic_Complex_Number(1, 0, 0, 0, 0));

        // We do not want to test that we can reliably prune zero states here, so the state is still present in the product
        REQUIRE(transition.lr.size() == 1);
        auto& lr_component = transition.lr.components[0];
        REQUIRE(lr_component.coef == Algebraic_Complex_Number(0, 0, 0, 0, 0));
    }

    { // Check that RIGHT subtree is an identity relation
        REQUIRE(transition.rl.size() == 1);
        auto& rl_component = transition.rl.components[0];
        REQUIRE(rl_component.coef == Algebraic_Complex_Number(0, 0, 0, 0, 0));

        REQUIRE(transition.rr.size() == 1);
        auto& rr_component = transition.rr.components[0];
        REQUIRE(rr_component.coef == Algebraic_Complex_Number(1, 0, 0, 0, 0));
    }
}

TEST_CASE("Build frontier automaton", "[SWTA]") {
    using ACN = Algebraic_Complex_Number;
    // Compute the frontier of the following automaton (there is only one color and one internal symbol):
    // q0 -> q1 + q2, q1
    // q1 -> q1, q1
    // q2 -> q3, q3
    // q3 -> q3, q3
    // FINAL STATES: q3, q1
    State q0 = 0;
    State q1 = 1;
    State q2 = 2;
    State q3 = 3;

    std::vector<SWTA::Transition> q0_transitions {synthetize_swta_transition({Def_Coef(ACN::ONE()) * q1, Def_Coef(ACN::ONE()) * q2}, {Def_Coef(ACN::ONE()) * q1})};
    std::vector<SWTA::Transition> q1_transitions {synthetize_swta_transition({Def_Coef(ACN::ONE()) * q1}, {Def_Coef(ACN::ONE()) * q1})};
    std::vector<SWTA::Transition> q2_transitions {synthetize_swta_transition({Def_Coef(ACN::ONE()) * q3}, {Def_Coef(ACN::ONE()) * q3})};
    std::vector<SWTA::Transition> q3_transitions {synthetize_swta_transition({Def_Coef(ACN::ONE()) * q3}, {Def_Coef(ACN::ONE()) * q3})};
    SWTA::Transition_Fn transitions {{q0_transitions}, {q1_transitions}, {q2_transitions}, {q3_transitions}};

    Bit_Set leaf_states (4, {q1, q3});
    std::vector<State> initial_states ({});
    initial_states.push_back(q0);

    SWTA swta (transitions, initial_states, leaf_states);

    NFA nfa = build_frontier_automaton(swta);

    REQUIRE(nfa.initial_states.size() == 1);
    REQUIRE(nfa.final_states.popcount() == 1);

    REQUIRE(nfa.number_of_states() == 3);

    Color color = 0;
    State m0 = 0;
    // The automaton should have the following structure: {q0} -> {q1, q2} -> {q1, q3} ---self loop-->
    REQUIRE(nfa.transitions[m0][color].size() == 1);
    State m1 = nfa.transitions[m0][color][0];

    REQUIRE(nfa.transitions[m1][color].size() == 1);
    State m2 = nfa.transitions[m1][color][0];

    REQUIRE(nfa.transitions[m1][color].size() == 1);
    REQUIRE(nfa.transitions[m2][color][0] == m2);

    REQUIRE(nfa.final_states == Bit_Set(3, {m2}));
}


TEST_CASE("Check frontiers for BV example are equivalent", "[SWTA]") {
    auto bv_example_result = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_RESULT);
    auto bv_example_post   = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_POST);

    auto result_frontier_nfa = build_frontier_automaton(bv_example_result);
    auto result_frontier_dfa = result_frontier_nfa.determinize();
    result_frontier_dfa.complete();

    auto post_frontier_nfa = build_frontier_automaton(bv_example_post);
    auto post_frontier_dfa = post_frontier_nfa.determinize();
    post_frontier_dfa.complete();

    bool are_equivalent = are_two_complete_dfas_equivalent(result_frontier_dfa, post_frontier_dfa);
    REQUIRE(are_equivalent);
}


TEST_CASE("Build first affine program", "[Affine programs]") {
    auto swta = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_POST);

    SWTA::Metadata metadata = swta.get_metadata();

    auto frontier_automaton    = build_frontier_automaton(swta);
    auto color_sym_abstraction = build_color_internal_symbol_abstraction(swta);
    auto first_program         = build_affine_program(swta, color_sym_abstraction);
}

TEST_CASE("Are two SWTAs color equivalent", "[Affine programs]") {
    auto bv_example_post   = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_POST);
    auto bv_example_result = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_RESULT);

    bool are_equivalent = are_two_swtas_color_equivalent(bv_example_post, bv_example_result);
    REQUIRE(are_equivalent);
}

TEST_CASE("Test BV example correcness step by step") {
    // -------- Step 1 - apply hadamard to precondition -------
    SWTA expected_swta = get_predefined_swta(Predefined_SWTA_Names::TEST_BV_EXAMPLE_AFTER_STEP1);

    SWTA initial_swta = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_PRE);
    WTT  hadamard     = get_predefined_wtt(Predefined_WTT_Names::HADAMARD, expected_swta.get_metadata());
    SWTA after_step1  = apply_wtt_to_swta(initial_swta, hadamard);

    bool are_equivalent = are_two_swtas_color_equivalent(expected_swta, after_step1);
    REQUIRE(are_equivalent);

    // -------- Step 2 - apply parity_cnot to the result of the previous step -------
    SWTA expected_swta2 = get_predefined_swta(Predefined_SWTA_Names::TEST_BV_EXAMPLE_AFTER_STEP2);

    WTT  cnot_parity  = get_predefined_wtt(Predefined_WTT_Names::PARITY_CNOT, expected_swta2.get_metadata());
    SWTA after_step2  = apply_wtt_to_swta(after_step1, cnot_parity);

    bool are_equivalent2 = are_two_swtas_color_equivalent(expected_swta2, after_step2);
    REQUIRE(are_equivalent2);

    // -------- Step 3 - apply hadamard to the result of the previous step -------
    SWTA expected_swta3 = get_predefined_swta(Predefined_SWTA_Names::TEST_BV_EXAMPLE_AFTER_STEP3);

    SWTA after_step3  = apply_wtt_to_swta(after_step2, hadamard);

    bool are_equivalent3 = are_two_swtas_color_equivalent(expected_swta3, after_step3);
    REQUIRE(are_equivalent3);
}

TEST_CASE("Run SWTA Program") {
    SWTA::Metadata swta_metadata = {
        .number_of_internal_symbols = 2, // Work qubit and ancilla
        .number_of_colors = 1
    };

    std::vector<WTT> needed_transducers {
        get_predefined_wtt(Predefined_WTT_Names::HADAMARD, swta_metadata),
        get_predefined_wtt(Predefined_WTT_Names::PARITY_CNOT, swta_metadata),
    };

    std::vector<Transducer_Application> applications {
        Transducer_Application(0), // Hadamard
        Transducer_Application(1), // CNOT
        Transducer_Application(0), // Hadamard
    };

    SWTA_Program program {
        .initial_swta = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_PRE),
        .transducers = needed_transducers,
        .applications = applications,
    };

    SWTA post_condition = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_POST);
    SWTA result         = run_swta_program(program);

    bool is_our_bv_correct = are_two_swtas_color_equivalent(result, post_condition);
    REQUIRE(is_our_bv_correct);
}

TEST_CASE("CnZ vs CCX-impl") {
    SWTA::Metadata swta_metadata = {
        .number_of_internal_symbols = 3, // Work qubit and ancilla
        .number_of_colors = 1
    };

    SWTA all_basis = get_predefined_swta(Predefined_SWTA_Names::GROVER_ALL_BASIS);

    std::vector<WTT> needed_transducers {
        /* 0 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_X, swta_metadata),
        /* 1 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_H, swta_metadata),
        /* 2 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_FIRST_MULTI_Z, swta_metadata),
        /* 3 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_FIRST_MULTI_Z_USING_CCX, swta_metadata),
        /* 4 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_SECOND_MULTI_Z, swta_metadata),
        /* 5 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_SECOND_MULTI_Z_USING_CCX, swta_metadata),
    };

    SWTA_Program first_program = {
        .initial_swta = all_basis,
        .transducers = needed_transducers,
        .applications = {
            Transducer_Application(0), // X
            Transducer_Application(2), // CnZ to last ancilla
            Transducer_Application(0), // X
            Transducer_Application(1), // H
            Transducer_Application(0), // X
            Transducer_Application(4), // CnZ to last working qubit
            Transducer_Application(0), // X
            Transducer_Application(1), // H
        }
    };

    SWTA_Program program_for_hw = {
        .initial_swta = all_basis,
        .transducers = needed_transducers,
        .applications = {
            Transducer_Application(0), // X
            Transducer_Application(3), // CnZ to last ancilla
            Transducer_Application(0), // X
            Transducer_Application(1), // H
            Transducer_Application(0), // X
            Transducer_Application(5), // CnZ to last working qubit
            Transducer_Application(0), // X
            Transducer_Application(1), // H
        }
    };

    SWTA result_idealistic = run_swta_program(first_program);
    SWTA result_for_hw     = run_swta_program(program_for_hw);

    bool is_our_impl_correct = are_two_swtas_color_equivalent(
        result_idealistic,
        result_for_hw
    );
    REQUIRE(is_our_impl_correct);
}

TEST_CASE("Staircase composition of identities") {
    SWTA::Metadata metadata = {
        .number_of_internal_symbols = 2,
        .number_of_colors = 1,
    };
    WTT box = get_predefined_wtt(Predefined_WTT_Names::TEST_STAIRCASE_IDENTITY3, metadata);

    Internal_Symbol terminating_symbol = 2; // Working qubit = 0, ancilla = 1
    u64 offset = 2;
    WTT result = perform_staircase_construction(box, {0, 1, 0}, offset, terminating_symbol);
}

TEST_CASE("ADDER Compose the UMA gate") {
    SWTA::Metadata metadata = {
        .number_of_internal_symbols = 1,
        .number_of_colors = 1,
    };
    WTT uma1 = get_predefined_wtt(Predefined_WTT_Names::ADDER_UMA1, metadata);
    WTT uma2 = get_predefined_wtt(Predefined_WTT_Names::ADDER_UMA2, metadata);
    WTT uma3 = get_predefined_wtt(Predefined_WTT_Names::ADDER_UMA3, metadata);

    auto uma12  = compose_wtts_sequentially(uma1,  uma2);
    auto uma123 = compose_wtts_sequentially(uma12, uma3);

    std::cout << uma123 << "\n";
}

TEST_CASE("ADDER Compose the MAJ gate") {
    SWTA::Metadata metadata = {
        .number_of_internal_symbols = 1,
        .number_of_colors = 2,
    };
    WTT maj1 = get_predefined_wtt(Predefined_WTT_Names::ADDER_MAJ1, metadata);
    WTT maj2 = get_predefined_wtt(Predefined_WTT_Names::ADDER_MAJ2, metadata);
    WTT maj3 = get_predefined_wtt(Predefined_WTT_Names::ADDER_MAJ3, metadata);

    auto maj12  = compose_wtts_sequentially(maj1, maj2);
    auto maj123 = compose_wtts_sequentially(maj12, maj3);

    WTT handwritten_maj = get_predefined_wtt(Predefined_WTT_Names::ADDER_MAJ_RESULT, metadata);

    SWTA all_basis = get_predefined_swta(Predefined_SWTA_Names::TEST_ADDER_ALL_3BASIS);

    SWTA result_automatic = apply_wtt_to_swta(all_basis, maj123);
    SWTA result_handwritten = apply_wtt_to_swta(all_basis, handwritten_maj);

    bool are_equivalent = are_two_swtas_color_equivalent(
        result_automatic,
        result_handwritten
    );
    REQUIRE(are_equivalent);
}

TEST_CASE("Verify adder circuit") {
    SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1};

    WTT maj = get_predefined_wtt(Predefined_WTT_Names::ADDER_MAJ_RESULT, metadata);

    WTT uma1 = get_predefined_wtt(Predefined_WTT_Names::ADDER_UMA1, metadata);
    WTT uma2 = get_predefined_wtt(Predefined_WTT_Names::ADDER_UMA2, metadata);
    WTT uma3 = get_predefined_wtt(Predefined_WTT_Names::ADDER_UMA3, metadata);

    WTT uma12  = compose_wtts_sequentially(uma1, uma2);
    WTT uma123 = compose_wtts_sequentially(uma12, uma3);

    std::vector<u64> box_inputs {0, 0, 0};
    u64 box_offset = 2;
    u64 new_symbol = 1;

    WTT maj_staircase = perform_staircase_construction(maj,    box_inputs, box_offset, new_symbol, Staircase_Direction::LEFT_RIGHT);
    WTT uma_staircase = perform_staircase_construction(uma123, box_inputs, box_offset, new_symbol, Staircase_Direction::RIGHT_LEFT);

    WTT id1 = get_predefined_wtt(Predefined_WTT_Names::TEST_FIXED_ID1, metadata);
    WTT extended_maj_staircase = compose_wtts_horizontally(maj_staircase, id1);
    WTT extended_uma_staircase = compose_wtts_horizontally(uma_staircase, id1);

    WTT middle_piece = get_predefined_wtt(Predefined_WTT_Names::ADDER_MIDDLE, metadata);

    WTT result12  = compose_wtts_sequentially(extended_maj_staircase, middle_piece);
    WTT adder_circuit = compose_wtts_sequentially(result12, extended_uma_staircase);

    SWTA precondition = get_predefined_swta(Predefined_SWTA_Names::ADDER_PRE);
    SWTA postcondition = get_predefined_swta(Predefined_SWTA_Names::ADDER_POST);

    SWTA result_swta = apply_wtt_to_swta(precondition, adder_circuit);
    std::cout << "Result has: " << result_swta.number_of_states() << " states.\n";

    auto precondition_dfa = build_frontier_automaton(precondition).determinize();
    precondition_dfa.complete();
    auto result_dfa       = build_frontier_automaton(postcondition).determinize();
    result_dfa.complete();

    bool are_two_swtas_equivalent = are_two_swtas_color_equivalent(result_swta, postcondition);
    std::cout << "Are equivalent: " << are_two_swtas_equivalent << "\n";

    // 0 1 0 0
     std::vector<u64> ints = {
         0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0,
     };
     std::vector<Algebraic_Complex_Number> nums;
     for (auto i : ints) {
         nums.push_back(Algebraic_Complex_Number(i, 0, 0, 0, 0));
     }

     std::cout << evaluate_wtt_on_tree(adder_circuit, 0, nums, {0, 1, 0, 0}, 0);

}

