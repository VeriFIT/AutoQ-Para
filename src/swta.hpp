#pragma once

#include "arith.hpp"
#include "bit_set.hpp"
#include "nfa.hpp"

#include <vector>

typedef u64 Color;
typedef u64 State;
typedef u64 Internal_Symbol;


struct Linear_Form {
    struct Component {
        Algebraic_Complex_Number coef;
        State state;

        Component() : coef(), state(0) {}
        Component(const Algebraic_Complex_Number& number, State state) : coef(number), state(state) {}

        void swap(Component& other) noexcept {
            this->coef.swap(other.coef);
            std::swap(this->state, other.state);
        }
    };

    std::vector<Component> components;  // Kept sorted by component.state

    Linear_Form() {}
    Linear_Form(const std::vector<Component>& components) : components(components) {}

    size_t size() const {
        return components.size();
    }

    bool empty() const {
        return components.empty();
    }

    void normalize() {
        for (auto& component : this->components) {
            component.coef.normalize();
        }
    }
};

/**
 * Synchronized Weighted Tree Automaton
 */
struct SWTA {
    /**
     * Transition from a state along a color.
     *
     * If any of the linear forms are empty, then the state does not have a transition for the given color.
     */
    struct Transition {
        Linear_Form left;
        Linear_Form right;

        Transition(Linear_Form& l, Linear_Form& r) : left(l), right(r) {}

        bool is_present() const {
            bool is_anything_missing = left.empty() || right.empty();
            return !is_anything_missing;
        }
    };

    using Transitions_From_State = std::vector<Transition>;

    std::vector<Transitions_From_State> transitions;  // transitions[state] are the transitions from a state
    std::vector<State> initial_states;
    Bit_Set states_with_leaf_transitions;

    SWTA(const std::vector<Transitions_From_State>& aut_transitions, std::vector<State>& init_states, const Bit_Set& final_state_set) :
        transitions(aut_transitions),
        initial_states(init_states),
        states_with_leaf_transitions(final_state_set) {}

    SWTA(const std::vector<Transitions_From_State>& aut_transitions, std::vector<State>& init_states, const std::vector<State>& final_state_set) :
        transitions(aut_transitions),
        initial_states(init_states),
        states_with_leaf_transitions(aut_transitions.size(), final_state_set) {}

    u64 number_of_states() const {
        return transitions.size();
    }

    u64 number_of_colors() const {
        if (this->transitions.empty()) { // There are no states, so the automaton is weird
            return 0;
        }
        return transitions[0].size();
    }
};

std::ostream& operator<<(std::ostream& os, const Linear_Form& form);
std::ostream& operator<<(std::ostream& os, const SWTA::Transition& swta_transition);
std::ostream& operator<<(std::ostream& os, const SWTA& swta);


/**
 * Weighted Tree Transducer
 */
struct WTT {
    struct Transition {
        Linear_Form ll;  // Left successor, its left  subtree inputs
        Linear_Form lr;  // Left successor, its right subtree inputs
        Linear_Form rl;  // Right successor, its left subtree inputs
        Linear_Form rr;  // Right successor, its right subtree inputs

        Transition(Linear_Form ll, Linear_Form lr, Linear_Form rl, Linear_Form rr) : ll(ll), lr(lr), rl(rl), rr(rr) {}

        bool is_present() const {
            bool is_left_present  = !ll.empty() || !lr.empty();
            bool is_right_present = !rl.empty() || !rr.empty();

            return is_left_present && is_right_present;
        }
    };

    using Transitions = std::vector<std::vector<Transition>>;

    /**
     *  Indexed by state: transitions[internal_symbol][state] are transitions from the state `state` along the internal symbol `internal_symbol`.
     */
    Transitions transitions;

    std::vector<State> initial_states;

    /**
     *  Tells which states can transform leaf.
     *  If bit_set[state] = 1, then the given state has a transition of the form state(leaf) -> leaf
     */
    Bit_Set states_with_leaf_transitions;

    WTT(const Transitions& transitions, const std::vector<State>& states_with_leaf_transitions, const std::vector<State>& initial_states) : transitions(transitions), states_with_leaf_transitions(0, nullptr), initial_states(initial_states)
    {
        Bit_Set state_set (transitions.size());
        for (State state : states_with_leaf_transitions) {
            state_set.set_bit(state, true);
        }
        this->states_with_leaf_transitions = state_set;
    };

    size_t number_of_states() const {
        return transitions[0].size();
    }

    size_t number_of_colors() const {
        return 1;  // @Todo: Figure out what are colors and what are internal symbols
    }

    void normalize_all_transitions() {
        for (Internal_Symbol internal_symbol = 0; internal_symbol < this->transitions.size(); internal_symbol++) {
            std::vector<Transition>& transitions_for_symbol = this->transitions[internal_symbol];

            for (auto& transition : transitions_for_symbol) {
                transition.ll.normalize();
                transition.lr.normalize();
                transition.rl.normalize();
                transition.rr.normalize();
            }
        }
    }

    u64 get_number_of_internal_symbols() const {
        return this->transitions.size();
    }

    bool does_state_accept_trees_for_any_colored_sequence(State state) const;
    void remove_zeros_from_transitions();
};

std::ostream& operator<<(std::ostream& os, const WTT::Transition& wtt_transition);
std::ostream& operator<<(std::ostream& os, const WTT& wtt);

/**
 * Compute a transducer that is equivalent to first appyling the `first` transducer and then the `second`.
 */
WTT compose_wtts_sequentially(WTT& first, WTT& second);

/**
 * Compute frontier NFA for a given SWTA. If root > 0, start with the given root state.
 */
template <typename Tree_Transition_System>
NFA build_frontier_automaton(const Tree_Transition_System& tts, s64 root = -1);
