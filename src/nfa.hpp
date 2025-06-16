#pragma once

#include "arith.hpp"
#include "bit_set.hpp"


struct NFA {
    using State = u64;
    using Transitions_From_State = std::vector<std::vector<State>>;

    std::vector<State> initial_states;
    Bit_Set final_states;
    std::vector<Transitions_From_State> transitions;

    NFA() : initial_states({}), final_states(0), transitions({}) {}
    NFA(const std::vector<State>& init_states, const std::vector<State>& fin_states, const std::vector<Transitions_From_State>& aut_transitions) :
        initial_states(init_states),
        final_states(aut_transitions.size()),
        transitions(aut_transitions)
    {
        for (auto state: fin_states) {
            final_states.set_bit(state, true);
        }
    }

    NFA(const std::vector<State>& init_states, const Bit_Set& fin_states, const std::vector<Transitions_From_State>& aut_transitions) :
        initial_states(init_states),
        final_states(fin_states),
        transitions(aut_transitions) {}

    u64 number_of_states() const {
        return this->transitions.size();
    }
};
