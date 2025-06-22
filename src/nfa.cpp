#include "nfa.hpp"

#include <map>


std::ostream& operator<<(std::ostream& os, const Macrostate& macrostate) {
    os << "Macrostate " << macrostate.state_names << " (handle=" << macrostate.handle << ")";
    return os;
}


Macrostate compute_post(const Macrostate* macrostate, const NFA& nfa, u64 symbol) {
    Macrostate post (nfa.number_of_states());

    for (NFA::State state : macrostate->state_names) {
        auto& transitions_from_state = nfa.transitions[state];
        auto& transitions_for_symbol = transitions_from_state[symbol];

        for (auto post_state : transitions_for_symbol) {
            post.state_set.set_bit(post_state);
        }
    }

    for (NFA::State state = 0; state < nfa.number_of_states(); state++) {
        if (post.state_set.get_bit_value(state)) {
            post.state_names.push_back(state);
        }
    }

    return post;
}

bool NFA::complete() {
    u64 state_cnt = this->number_of_states();
    State sink_state = state_cnt;

    bool was_sink_state_needed = false;
    for (State state = 0; state < state_cnt; state++) {
        for (u64 color = 0; color < this->alphabet_size(); color++) {
            if (this->transitions[state][color].empty()) {
                this->transitions[state][color].push_back(sink_state);
                was_sink_state_needed = true;
            }
        }
    }

    if (was_sink_state_needed) {
        std::vector<std::vector<State>> sink_state_transitions;
        sink_state_transitions.resize(this->alphabet_size());

        for (u64 color = 0; color < this->alphabet_size(); color++) {
            sink_state_transitions[color].push_back(sink_state);
        }

        this->transitions.push_back(sink_state_transitions);
    }

    return was_sink_state_needed;
}

NFA NFA::determinize() const {
    std::map<Macrostate, NFA::State> handles;

    std::vector<const Macrostate*> worklist;

    { // Initialize frontier
        Macrostate initial_macrostate (this->number_of_states(), initial_states);
        initial_macrostate.handle = 0;
        auto [insert_pos, was_inserted] = handles.emplace(initial_macrostate, 0);
        worklist.push_back(&insert_pos->first);
    }

    std::map<State, std::vector<std::vector<State>>> resulting_transitions; // Use an associative container here, because we do now know the final number of states
    Bit_Set final_macrostates (0);

    while (!worklist.empty()) {
        auto macrostate = worklist.back();
        worklist.pop_back();

        if (this->final_states.is_superset(macrostate->state_set)) { // All of the states in macrostate can make a leaf transition
            final_macrostates.grow_and_set_bit(macrostate->handle);
        }

        for (u64 symbol = 0; symbol < this->alphabet_size(); symbol++) {
            Macrostate post = compute_post(macrostate, *this, symbol);

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
            if (transitions_from_this_macrostate.size() < this->alphabet_size()) {
                transitions_from_this_macrostate.resize(this->alphabet_size());
            }

            transitions_from_this_macrostate[symbol].push_back(post.handle);
        }
    }

    std::vector<NFA::Transitions_From_State> ordered_resulting_transitions;
    ordered_resulting_transitions.resize(resulting_transitions.size());

    for (auto& [state, transitions_from_state] : resulting_transitions) {
        ordered_resulting_transitions[state] = transitions_from_state;
    }

    return NFA({0}, final_macrostates, ordered_resulting_transitions);
}

bool NFA::is_every_state_accepting() const {
    return this->final_states.are_all_bits_set();
}

void NFA::write_dot(std::ostream& stream) const {
    stream << "digraph NFA {\n";

    for (State initial_state : this->initial_states) {
        stream << "  qInit" << initial_state << " [shape=none, label=\"\"]\n";
    }
    for (State state = 0; state < this->number_of_states(); state++) {
        stream << "  q" << state;
        if (this->debug_data != nullptr && this->debug_data->state_names.contains(state)) {
            stream << " [label=\"" << this->debug_data->state_names.at(state) << "\"]";
        }
        stream << "\n";
    }

    for (State initial_state : this->initial_states) {
        stream << "  qInit" << initial_state << " -> q" << initial_state << "\n";
    }

    for (State state = 0; state < number_of_states(); state++) {
        for (u64 symbol = 0; symbol < transitions[state].size(); symbol++) {
            for (State destination : transitions[state][symbol]) {
                stream << "  q" << state << " -> q" << destination << " [label=\"" << symbol << "\"]\n";
            }
        }
    }

    stream << "}";
}

