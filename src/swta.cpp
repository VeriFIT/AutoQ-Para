#include "swta.hpp"
#include "arith.hpp"
#include "bit_set.hpp"
#include "nfa.hpp"

#include <vector>
#include <map>

std::ostream& operator<<(std::ostream& os, const Linear_Form& form) {
    os << "Linear_Form{ ";
    for (u64 i = 0; i < form.components.size(); i++) {
        auto& component = form.components[i];
        os << "(" << component.coef << ")" << "*q" << component.state;
        if (i != form.components.size() - 1) {
            os << " ";
        }
    }
    os << "}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const SWTA::Transition& swta_transition) {
    os << "LEFT: " << swta_transition.left << ", "
       << "RIGHT: " << swta_transition.right;
    return os;
}

std::ostream& operator<<(std::ostream& os, const SWTA& swta) {
    for (u64 state = 0; state < swta.transitions.size(); state++) {
        auto& transitions = swta.transitions[state];
        s64 i = -1;
        os << "State=" << state << "{ ";
        for (auto& [color, transition] : transitions) {
            i++;

            os << "Color@" << color << "[ " << transition << "]";

            if (i != transitions.size()) {
                os << ", ";
            }
        }
        os << "}"; // End of State=x { ... transitions ...
    }
    return os;
}


std::ostream& operator<<(std::ostream& os, const Linear_Form::Component& component) {
     os << component.coef << " * q" << component.state;
     return os;
}


void write_norm_with_subtree_info(std::ostream& target, const Linear_Form& form, const char* subtree_info, bool needs_leading_plus) {
    if (needs_leading_plus) {
        target << " + ";
    }

    for (u64 i = 0; i < form.size(); i++) {
        target << form.components[i] << subtree_info;
        if (i < form.size() - 1) target << " + ";
    }
}

std::ostream& operator<<(std::ostream& os, const WTT::Transition& wtt_transition) {
    os << "LEFT SUBTREE: ";
    bool needs_plus = false; // True, if anything has been written to the output stream
    if (!wtt_transition.ll.empty()) {
        write_norm_with_subtree_info(os, wtt_transition.ll, "(L)", needs_plus);
        needs_plus = true;
    }

    if (!wtt_transition.lr.empty()) {
        write_norm_with_subtree_info(os, wtt_transition.lr, "(R)", needs_plus);
        needs_plus = true;
    }

    if (wtt_transition.ll.empty() && wtt_transition.lr.empty()) {
        os << "0";
    }

    os << "; RIGHT SUBTREE: ";
    needs_plus = false; // Reset

    if (!wtt_transition.rl.empty()) {
        write_norm_with_subtree_info(os, wtt_transition.rl, "(L)", needs_plus);
        needs_plus = true;
    }

    if (!wtt_transition.rr.empty()) {
        write_norm_with_subtree_info(os, wtt_transition.rr, "(R)", needs_plus);
        needs_plus = true;
    }

    if (wtt_transition.rl.empty() && wtt_transition.rr.empty()) {
        os << "0";
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const WTT& wtt) {
    os << "WTT {\n";
    os << "  initial states: " << wtt.initial_states << "\n";

    for (Internal_Symbol sym = 0; sym < wtt.transitions.size(); sym++) {
        const std::vector<WTT::Transition>& transition_for_sym = wtt.transitions[sym];
        for (u64 state = 0; state < wtt.number_of_states(); state++) {
            os << "  " << state << "--(sym=" << sym << ")-->: " << transition_for_sym[state] << "\n";
        }
    }

    os << "}";
    return os;
}


struct State_Pair {
    State first, second;

    /**
     * Lexigraphical ordering.
     */
    bool operator<(const State_Pair& other) const {
        if (this->first < other.first)
            return true;
        else if (this->first > other.first)
            return false;
        return this->second < other.second;
    }
};


void extend_form_with_product_and_node_discoveries(Linear_Form& destination, Linear_Form& first, Linear_Form& second, std::map<State_Pair, u64>& handles, std::vector<State_Pair>& worklist) {
    for (auto& outer_comp : second.components) {
        for (auto& inner_comp : first.components) {
            State_Pair state {inner_comp.state, outer_comp.state};
            auto [insert_pos, was_inserted] = handles.emplace(state, handles.size());
            State handle = insert_pos->second;

            if (was_inserted) {
                worklist.push_back(state);
            }

            Algebraic_Complex_Number coef = inner_comp.coef * outer_comp.coef;

            bool already_present = false;
            for (auto& component : destination.components) {
                if (component.state == handle) {
                    component.coef += coef;
                    already_present = true;
                    break;
                }
            }

            if (already_present) {
                continue;  // We are done here
            }

            destination.components.push_back({coef, handle});
        }
    }
}


WTT compose_wtts_sequentially(WTT& first, WTT& second) {
    std::vector<State_Pair>   worklist;

    std::map<State_Pair, u64>        state_handles;
    std::vector<State>               initial_states;
    std::vector<std::map<State, WTT::Transition>> transitions;  // indexed by internal symbols

    const u64 num_of_internal_symbols = first.get_number_of_internal_symbols();
    assert(num_of_internal_symbols == second.get_number_of_internal_symbols());
    transitions.resize(num_of_internal_symbols);

    u64 init_state_cnt = first.initial_states.size()*second.initial_states.size();
    initial_states.reserve(init_state_cnt);
    for (State state = 0; state < init_state_cnt; state++) {
        initial_states.push_back(state);
    }

    for (auto first_state : first.initial_states) {
        for (auto second_state : second.initial_states) {
            worklist.push_back({first_state, second_state});

            State_Pair pair = {first_state, second_state};
            state_handles.emplace(pair, state_handles.size());
        }
    }

    while (!worklist.empty()) {
        State_Pair state_pair = worklist.back();
        worklist.pop_back();

        State handle = state_handles.at(state_pair);

        for (Internal_Symbol internal_symbol = 0; internal_symbol < num_of_internal_symbols; internal_symbol++) {
            WTT::Transition first_transition   = first.transitions[internal_symbol][state_pair.first];
            WTT::Transition second_transitions = second.transitions[internal_symbol][state_pair.second];

            Linear_Form ll, lr, rl, rr;
            { // ll
                extend_form_with_product_and_node_discoveries(ll, first_transition.ll, second_transitions.ll, state_handles, worklist);
                extend_form_with_product_and_node_discoveries(ll, first_transition.rl, second_transitions.lr, state_handles, worklist);
            }

            { // lr
                extend_form_with_product_and_node_discoveries(lr, first_transition.ll, second_transitions.lr, state_handles, worklist);
                extend_form_with_product_and_node_discoveries(lr, first_transition.lr, second_transitions.rr, state_handles, worklist);
            }

            { // rl
                extend_form_with_product_and_node_discoveries(rl, first_transition.ll, second_transitions.rl, state_handles, worklist);
                extend_form_with_product_and_node_discoveries(rl, first_transition.rl, second_transitions.rr, state_handles, worklist);
            }

            { // rr
                extend_form_with_product_and_node_discoveries(rr, first_transition.lr, second_transitions.rl, state_handles, worklist);
                extend_form_with_product_and_node_discoveries(rr, first_transition.rr, second_transitions.rr, state_handles, worklist);
            }

            WTT::Transition resulting_transition (ll, lr, rl, rr);
            transitions[internal_symbol].emplace(handle, resulting_transition);
        }
    }

    WTT::Transitions resulting_transitions;
    resulting_transitions.reserve(num_of_internal_symbols);

    for (Internal_Symbol sym = 0; sym < num_of_internal_symbols; sym++) {
        std::vector<WTT::Transition> transitions_for_sym;

        for (auto& [state, transition] : transitions[sym]) {
            transitions_for_sym.push_back(transition);
        }

        resulting_transitions.push_back(transitions_for_sym);
    }

    WTT result (resulting_transitions, {}, initial_states);
    return result;
}



Macrostate compute_post(const Macrostate* macrostate, const SWTA& swta, Color color) {
    Macrostate post(swta.number_of_states());

    for (State state : macrostate->state_names) {
        auto& transitions_from_state = swta.transitions[state];
        auto& transitions_for_color  = transitions_from_state[color];

        if (!transitions_for_color.is_present()) {
            post.state_set.clear();
            break;
        }

        for (auto& component : transitions_for_color.left.components) {
            post.state_set.set_bit(component.state);
        }

        for (auto& component : transitions_for_color.right.components) {
            post.state_set.set_bit(component.state);
        }
    }

    for (State state = 0; state < swta.number_of_states(); state++) {
        if (!post.state_set.get_bit_value(state)) {
            continue;
        }

        post.state_names.push_back(state);
    }

    return post;
}

Macrostate compute_post(const Macrostate* macrostate, const WTT& wtt, Color color) {
    Macrostate post(wtt.number_of_states());

    for (State state : macrostate->state_names) {
        auto& transitions_from_state = wtt.transitions[state];
        auto& transitions_for_color  = transitions_from_state[color];

        if (!transitions_for_color.is_present()) {
            post.state_set.clear();
            break;
        }

        for (auto& component : transitions_for_color.ll.components) {
            post.state_set.set_bit(component.state);
        }

        for (auto& component : transitions_for_color.lr.components) {
            post.state_set.set_bit(component.state);
        }

        for (auto& component : transitions_for_color.rl.components) {
            post.state_set.set_bit(component.state);
        }

        for (auto& component : transitions_for_color.rr.components) {
            post.state_set.set_bit(component.state);
        }
    }

    for (State state = 0; state < wtt.number_of_states(); state++) {
        if (!post.state_set.get_bit_value(state)) {
            continue;
        }

        post.state_names.push_back(state);
    }

    return post;
}

void dump_discovered_transitions(const std::map<State, std::vector<std::vector<State>>>& transitions) {
    std::cout << "Known states: ";
    for (const auto& [state, transitions_from_state]: transitions) {
        std::cout << state << ", ";
    }
    std::cout << "\n";
}

void initialize_frontier_with_initial_states(std::vector<const Macrostate*>& worklist, std::map<Macrostate, NFA::State>& handles, const std::vector<State>& initial_states, u64 total_number_of_states, s64 root) {
    if (root < 0) {
        Macrostate initial_macrostate (total_number_of_states, initial_states);
        initial_macrostate.handle = 0;
        auto [insert_pos, was_inserted] = handles.emplace(initial_macrostate, 0); // We do not have any handles, so we know that the first will have value 0

        worklist.push_back(&insert_pos->first);
    } else { // Start the construction from the provided root
        std::vector<State> root_states ({static_cast<State>(root)});
        Macrostate initial_states (total_number_of_states, root_states);
        initial_states.handle = 0;

        auto [insert_pos, was_inserted] = handles.emplace(initial_states, 0); // We do not have any handles, so we know that the first will have value 0
        worklist.push_back(&insert_pos->first);
    }
}

template <typename Tree_Transition_System>
NFA build_frontier_automaton(const Tree_Transition_System& tts, s64 root) {
    std::map<Macrostate, NFA::State> handles;

    // @Note: Use pointers to avoid copying Bit_Sets into the worklist -- only one copy present in handles should be sufficient
    std::vector<const Macrostate*> worklist;

    initialize_frontier_with_initial_states(worklist, handles, tts.initial_states, tts.number_of_states(), root);

    u64 color_cnt = tts.number_of_colors();

    std::map<State, std::vector<std::vector<State>>> resulting_transitions; // Use an associative container here, because we do now know the final number of states
    Bit_Set final_macrostates (0);

    while (!worklist.empty()) {
        auto macrostate = worklist.back();
        worklist.pop_back();

        if (tts.states_with_leaf_transitions.is_superset(macrostate->state_set)) { // All of the states in macrostate can make a leaf transition
            final_macrostates.grow_and_set_bit(macrostate->handle);
        }

        for (Color color = 0; color < color_cnt; color++) {
            Macrostate post = compute_post(macrostate, tts, color);

            if (post.empty()) {
                continue;
            }

            post.handle = handles.size(); // Set this speculatively, before we store it in the handles map

            const auto& [insert_pos, was_inserted] = handles.emplace(post, handles.size());
            if (was_inserted) { // the macrostate is new, we need to explore it
                worklist.push_back(&insert_pos->first);
            } else { // There already is such a macrostate, so our speculation with the handle was incorrect
                post.handle = insert_pos->second;
            }

            auto& transitions_from_this_macrostate = resulting_transitions[macrostate->handle];
            if (transitions_from_this_macrostate.size() < color_cnt) {
                transitions_from_this_macrostate.resize(color_cnt);
            }

            transitions_from_this_macrostate[color].push_back(post.handle);
        }
    }

    std::vector<NFA::Transitions_From_State> ordered_resulting_transitions;
    ordered_resulting_transitions.resize(resulting_transitions.size());

    for (auto& [state, transitions_from_state] : resulting_transitions) {
        ordered_resulting_transitions[state] = transitions_from_state;
    }

    return NFA({0}, final_macrostates, ordered_resulting_transitions);
}

template
NFA build_frontier_automaton<SWTA>(const SWTA& tts, s64 root = -1);

template
NFA build_frontier_automaton<WTT>(const WTT& tts, s64 root = -1);

bool WTT::does_state_accept_trees_for_any_colored_sequence(State state) const {
    NFA accepted_colored_sequences_abstraction = build_frontier_automaton(*this);
    NFA determinized_abstraction = accepted_colored_sequences_abstraction.determinize();

    return determinized_abstraction.is_every_state_accepting();
}


enum class State_Universality_Status : u8 {
    UNKNOWN = 0,
    UNIVERSAL = 1,
    NONUNIVERSAL = 2,
};

bool can_component_be_removed(Linear_Form::Component& component, const WTT& wtt, std::vector<State_Universality_Status>& cache) {

    if (cache[component.state] != State_Universality_Status::UNKNOWN) {
        bool is_coef_zero = component.coef.is_zero();

        if (is_coef_zero && cache[component.state] == State_Universality_Status::UNIVERSAL) {
            return true;
        }
        return false;
    }

    bool is_coef_zero = component.coef.is_zero();
    bool is_state_universal = wtt.does_state_accept_trees_for_any_colored_sequence(component.state);

    cache[component.state] = is_state_universal ? State_Universality_Status::UNIVERSAL : State_Universality_Status::NONUNIVERSAL;

    return is_coef_zero && is_state_universal;
}

void remove_zeros_from_form(const WTT& wtt, Linear_Form& form, std::vector<State_Universality_Status>& cache) {
    s64 nonzero_idx = form.size() - 1;
    s64 zero_idx    = 0;

    if (nonzero_idx == zero_idx) { // there is only one element
        auto& component = form.components[zero_idx];

        if (can_component_be_removed(component, wtt, cache)) {
            form.components.clear();
        }
        return;
    }

    while (zero_idx < nonzero_idx) {
        // Search for the next zero slot that needs to be filled
        for (; zero_idx < form.size(); zero_idx++) {
            if (can_component_be_removed(form.components[zero_idx], wtt, cache)) {
                break;
            }
        }

        // Search for the next coef to fill the slot from the back
        for (; nonzero_idx >= 0; nonzero_idx--) {
            if (!can_component_be_removed(form.components[zero_idx], wtt, cache)) {
                break;
            }
        }
        if (zero_idx < nonzero_idx) break;

        form.components[zero_idx].swap(form.components[nonzero_idx]);
    }

    if (zero_idx >= form.size()) {
        form.components.resize(nonzero_idx);
    }
}

void WTT::remove_zeros_from_transitions() {
    std::vector<State_Universality_Status> cache;
    cache.resize(this->number_of_states());

    for (Internal_Symbol internal_symbol = 0; internal_symbol < this->transitions.size(); internal_symbol++) {
        std::vector<Transition>& transitions_for_symbol = this->transitions[internal_symbol];

        for (auto& transition : transitions_for_symbol) {
            remove_zeros_from_form(*this, transition.ll, cache);
            remove_zeros_from_form(*this, transition.lr, cache);
            remove_zeros_from_form(*this, transition.rl, cache);
            remove_zeros_from_form(*this, transition.rr, cache);
        }
    }
}
