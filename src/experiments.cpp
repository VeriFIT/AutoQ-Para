#include "basics.hpp"
#include "predefined_automata.hpp"
#include "quantum_program.hpp"
#include "swta.hpp"
#include "swta_builders.hpp"

#include <chrono>

#define RUN_EXPERIMENT(name, experiment_function) \
{ \
    std::cout << "Running experiment \"" << name << "\"..."; \
    auto begin = std::chrono::steady_clock::now(); \
    auto status = experiment_function(); \
    std::cout << " Done. Summary:\n"; \
    auto end = std::chrono::steady_clock::now(); \
    std::cout << " > Status : " << get_word_for_status(status) << "\n"; \
    std::cout << " > Runtime: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << "\n"; \
}

enum Verification_Status : u32 {
    FAILURE = 0,
    SUCCESS = 1,
};

const char* get_word_for_status(Verification_Status status) {
    return status == Verification_Status::FAILURE ? "FAIL" : "SUCCESS";
}

Verification_Status verify_bv() {
    // -------- Step 1 - apply hadamard to precondition -------
    SWTA precondition = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_PRE);
    auto metadata = precondition.get_metadata();

    WTT hadamard  = get_predefined_wtt(Predefined_WTT_Names::HADAMARD, metadata);
    WTT parity_cnot = get_predefined_wtt(Predefined_WTT_Names::PARITY_CNOT, metadata);
    std::vector<WTT> required_transducers = {
        hadamard,
        parity_cnot
    };

    SWTA_Program program = {
        .initial_swta = precondition,
        .transducers = required_transducers,
        .applications = {
            Transducer_Application(0),
            Transducer_Application(1),
            Transducer_Application(0),
        }
    };

    SWTA result = run_swta_program(program);
    SWTA postcondition = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_POST);

    bool are_equivalent = are_two_swtas_color_equivalent(result, postcondition);

    return are_equivalent ? Verification_Status::SUCCESS : Verification_Status::FAILURE;
}


Verification_Status verify_grover() {
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

    bool are_equivalent = are_two_swtas_color_equivalent(
        result_idealistic,
        result_for_hw
    );

    return are_equivalent ? Verification_Status::SUCCESS : Verification_Status::FAILURE;
}

Verification_Status verify_adder() {
    
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

    auto precondition_dfa = build_frontier_automaton(precondition).determinize();
    precondition_dfa.complete();
    auto result_dfa       = build_frontier_automaton(postcondition).determinize();
    result_dfa.complete();

    bool are_equivalent = are_two_swtas_color_equivalent(result_swta, postcondition);
    return are_equivalent ? Verification_Status::SUCCESS : Verification_Status::FAILURE;
}

Verification_Status verify_ecc() {
    auto ecc_pre  = get_predefined_swta(Predefined_SWTA_Names::ECC_PRE);
    auto ecc_post  = get_predefined_swta(Predefined_SWTA_Names::ECC_POST);

    SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 4 };

    auto box_part1 = get_predefined_wtt(Predefined_WTT_Names::ECC_BOX1, metadata);
    auto box_part2 = get_predefined_wtt(Predefined_WTT_Names::ECC_BOX2, metadata);
    auto box = compose_wtts_sequentially(box_part1, box_part2);

    u64 new_symbol = 1;
    u64 offset = 2;
    auto staircase = perform_staircase_construction(box, {0, 0, 0, 0}, offset, new_symbol);

    auto ecc_result = apply_wtt_to_swta(ecc_pre, staircase);

    bool are_equivalent = are_two_swtas_color_equivalent(ecc_result, ecc_post);
    return are_equivalent ? Verification_Status::SUCCESS : Verification_Status::FAILURE;
}

Verification_Status verify_hamiltonian_simulation() {
    using ACN = Algebraic_Complex_Number;

    SWTA::Metadata metadata;

    WTT rzz_box   = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_RZZ, metadata);
    WTT rzz_stage = perform_staircase_construction(rzz_box, {0, 0}, 1, 1);

    WTT rxx_box   = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_RXX, metadata);
    WTT rxx_stage = perform_staircase_construction(rxx_box, {0, 0}, 1, 1);

    WTT ryy_box   = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_RYY, metadata);
    WTT ryy_stage = perform_staircase_construction(ryy_box, {0, 0}, 1, 1);

    WTT uzz_box   = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_UZZ, metadata);
    WTT uzz_stage = perform_staircase_construction(uzz_box, {0, 0}, 1, 1);

    WTT sqrt_x_stage = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_SQRT_X_STAGE, metadata);
    WTT s_stage = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_S_STAGE, metadata);
    WTT h_stage = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_H_STAGE, metadata);
    WTT last_x_stage = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_LAST_X_STAGE, metadata);

    std::vector<WTT> circuit_stages = {
        h_stage,       // 0
        rzz_stage,     // 1
        sqrt_x_stage,  // 2
        uzz_stage,     // 3
        last_x_stage,  // 4
        s_stage,       // 5
        rxx_stage,     // 6
        ryy_stage,     // 7
    };

    SWTA precondition = get_predefined_swta(Predefined_SWTA_Names::HAMILTONIAN_ALL_BASIS);

    SWTA_Program naive_program = {
        .initial_swta = precondition,
        .transducers = circuit_stages,
        .applications = {
            Transducer_Application(0),
            Transducer_Application(1),
            Transducer_Application(2),
            Transducer_Application(3),
            Transducer_Application(4),
            Transducer_Application(0),
            Transducer_Application(5),
            Transducer_Application(1),
        },
    };

    SWTA_Program optimized_program = {
        .initial_swta = precondition,
        .transducers = circuit_stages,
        .applications = {
            Transducer_Application(6),
            Transducer_Application(7),
            Transducer_Application(1),
        },
    };

    auto naive_program_result     = run_swta_program(naive_program);
    auto optimized_program_result = run_swta_program(optimized_program);

    bool are_equivalent = are_two_swtas_color_equivalent(naive_program_result, optimized_program_result);
    return are_equivalent ? Verification_Status::SUCCESS : Verification_Status::FAILURE;
}


int main() {

    RUN_EXPERIMENT("Bernstein-Vazirani", verify_bv);
    RUN_EXPERIMENT("Grover (amplitude amplification)", verify_grover);
    RUN_EXPERIMENT("Adder", verify_adder);
    RUN_EXPERIMENT("Error correction code", verify_ecc);
    RUN_EXPERIMENT("Hamiltonian simulation", verify_hamiltonian_simulation);

    return 0;
}
