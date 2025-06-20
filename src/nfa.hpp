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

    u64 alphabet_size() const {
        return this->transitions[0].size();
    }

    NFA determinize() const;
    bool is_every_state_accepting() const;
};


struct Macrostate {
    using State = u64;

    std::vector<State> state_names; // Accelerates computation of Post
    Bit_Set            state_set;   // To test what states are present in the macrostate
    State              handle;      // Assigned after initialization

    Macrostate(u64 state_cnt) : state_names({}), state_set(state_cnt) {}
    Macrostate(u64 size, const std::vector<State>& content) : state_names(content), state_set(size, content) {}

    Macrostate(Macrostate&& other) : state_names(other.state_names), state_set(other.state_set), handle(other.handle) {}

    Macrostate(const Macrostate& other) : state_names(other.state_names), state_set(other.state_set), handle(other.handle) {}

    bool operator<(const Macrostate& other) const {
        return state_set < other.state_set;
    }

    Macrostate& operator=(Macrostate&& other) {
        this->handle      = other.handle;
        this->state_set   = std::move(other.state_set);
        this->state_names = std::move(other.state_names);

        return *this;
    }

    Macrostate& operator=(const Macrostate& other) {
        this->handle      = other.handle;
        this->state_set   = other.state_set;
        this->state_names = other.state_names;

        return *this;
    }

    bool empty() const {
        return state_names.empty();
    }

    void populate_state_names_from_set() {
        this->state_names.clear();
        for (State state = 0; state < this->state_set.size; state++) {
            if (this->state_set.get_bit_value(state)) {
                this->state_names.push_back(state);
            }
        }
    }
};

namespace std {
    template <typename T>
    ostream& operator<<(ostream& os, std::vector<T> vec) {
        os << "[";
        for (u64 idx = 0; idx < vec.size(); idx++) {
            auto& elem = vec[idx];
            std::cout << elem;
            if (idx + 1 != vec.size()) {
                os << ", ";
            }
        }
        os << "]";
        return os;
    }
}

std::ostream& operator<<(std::ostream& os, const Macrostate& macrostate);

