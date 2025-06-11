#pragma once

#include "arith.hpp"
#include "bit_set.hpp"

#include <vector>
#include <map>

typedef u64 Color;
typedef u64 State;


struct Linear_Form {
    struct Component {
        Algebraic_Complex_Number coef;
        State state;
        Component(const Algebraic_Complex_Number& number, State state) : coef(number), state(state) {}
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
};

/**
 * Synchronized Weighted Tree Automaton
 */
struct SWTA {
    struct Transition {
        Linear_Form left;
        Linear_Form right;
    };

    std::vector< std::map<Color, Transition> > transitions;  // transitions[state] are the transitions from a state
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
    };

    /**
     *  Indexed by state: transitions[state] are transitions from the state `state`.
     */
    std::vector<Transition> transitions;

    std::vector<State> initial_states;

    /**
     *  Tells which states can transform leaf.
     *  If bit_set[state] = 1, then the given state has a transition of the form state(leaf) -> leaf
     */
    Bit_Set states_with_leaf_transitions;

    WTT(const std::vector<Transition>& transitions, const std::vector<State>& states_with_leaf_transitions, const std::vector<State>& initial_states) : transitions(transitions), states_with_leaf_transitions(0, nullptr), initial_states(initial_states)
    {
        Bit_Set state_set (transitions.size());
        for (State state : states_with_leaf_transitions) {
            state_set.set_bit(state, true);
        }
        this->states_with_leaf_transitions = state_set;
    };

    size_t number_of_states() const {
        return transitions.size();
    }
};

std::ostream& operator<<(std::ostream& os, const WTT::Transition& wtt_transition);
std::ostream& operator<<(std::ostream& os, const WTT& wtt);

/**
 * Compute a transducer that is equivalent to first appyling the `first` transducer and then the `second`.
 */
WTT compose_wtts_sequentially(WTT& first, WTT& second);
