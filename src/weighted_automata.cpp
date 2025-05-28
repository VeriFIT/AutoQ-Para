#include "weighted_automata.hpp"
#include "arith.hpp"


typedef std::vector<u64> Word;


bool is_weighted_automaton_zero(const Weighted_Automaton& automaton) {
    u64 max_color = automaton.transitions.size();

    const u64 state_cnt = automaton.number_of_states();
    ACN_Matrix forward_space_basis (state_cnt, state_cnt);

    std::vector<std::tuple<Word, ACN_Matrix>> worklist;

    { // Seed the worklist
        Word       empty_word;
        ACN_Matrix state_distribution = automaton.initial_vector;

        worklist.push_back({empty_word, state_distribution});
    }

    Bit_Set active_states(state_cnt);

    while (!worklist.empty()) {
        auto worklist_tuple = worklist.at(worklist.size() - 1);
        worklist.pop_back();

        auto word         = std::get<0>(worklist_tuple);
        auto distribution = std::get<1>(worklist_tuple);

        s64 row_destination = add_row_to_row_echelon_matrix(forward_space_basis, distribution);
        if (row_destination < 0) continue;

        if (row_destination == static_cast<s64>(state_cnt)) break; // The matrix has already a full rank, there is no need to continue

        active_states.clear();
        for (u64 state = 0; state < state_cnt; state++) {
            if (!distribution.at(0, state).is_zero()) active_states.set_bit(state, true);
        }

        for (u64 color = 0; color < max_color; color++) {
            bool can_take_transition = automaton.can_make_transition_from_active_states(active_states, color);
            if (!can_take_transition) continue;

            Word extended_word (word);
            extended_word.push_back(color);

            ACN_Matrix extended_distribution = distribution * automaton.transitions[color];
            worklist.push_back({extended_word, extended_distribution});
        }
    }

    auto product = forward_space_basis * automaton.final_vector;

    return product.contains_only_zeros();
}
