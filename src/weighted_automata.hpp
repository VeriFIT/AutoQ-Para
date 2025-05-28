#pragma once

#include "arith.hpp"
#include "bit_set.hpp"
#include <vector>


struct Weighted_Automaton {
    ACN_Matrix initial_vector;
    ACN_Matrix final_vector;    // @Todo(m): What about the 'e' situation?

    /**
     * Transition matrices for colors.
     *
     * Indexed by color value (int).
     */
    std::vector<ACN_Matrix> transitions;

    /**
     * States that allow transition for a given color.
     *
     * Indexed by color value (int).
     */
    std::vector<Bit_Set> supported_colors;


    Weighted_Automaton(const ACN_Matrix& initial, const ACN_Matrix& final, const std::vector<ACN_Matrix>& transitions, const std::vector<Bit_Set>& support) : initial_vector(initial), final_vector(final), transitions(transitions), supported_colors(support) {

    };


    u64 number_of_states() const {
        return initial_vector.width;
    }

    bool can_make_transition_from_active_states(const Bit_Set& active_states, u64 color) const {
        auto& states_allowing_transition = supported_colors.at(color);
        return states_allowing_transition.is_superset(active_states);
    }
};

bool is_weighted_automaton_zero(const Weighted_Automaton& automaton);
