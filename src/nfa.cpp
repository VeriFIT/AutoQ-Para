#include "nfa.hpp"
#include "basics.hpp"

#include <map>
#include <algorithm>
#include <sstream>


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

        if (!this->final_states.is_intersection_empty(macrostate->state_set)) {
            // There is at least one state that can make a leaf transition
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
    ordered_resulting_transitions.resize(handles.size());

    for (auto& [state, transitions_from_state] : resulting_transitions) {
        ordered_resulting_transitions[state] = transitions_from_state;
    }

    for (State state = 0; state < handles.size(); state++) {
        if (ordered_resulting_transitions[state].empty()) {
            ordered_resulting_transitions[state].resize(this->alphabet_size());
        }
    }

    NFA result ({0}, final_macrostates, ordered_resulting_transitions);

    do_on_debug({
        result.debug_data = new NFA::Debug_Data;

        for (auto& [macrostate, handle] : handles) {
            std::stringstream string_stream;
            string_stream << "{";

            u64 written_states_cnt = 0;
            for (State state : macrostate.state_names) {
                std::string state_name = (this->debug_data->state_names.contains(state)) ?  this->debug_data->state_names.at(state) : std::to_string(state);
                string_stream << state_name;

                if (written_states_cnt + 1 != macrostate.state_names.size()) {
                    string_stream << ", ";
                }

                written_states_cnt += 1;
            }

            string_stream << "}";

            std::string macrostate_name = string_stream.str();
            result.debug_data->state_names[handle] = macrostate_name;
        }
    });

    return result;
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
        const char* color = this->final_states.get_bit_value(state) ? "green" : "black";
        if (this->debug_data != nullptr && this->debug_data->state_names.contains(state)) {
            stream << " ["
                   << "label=\"" << this->debug_data->state_names.at(state) << "\", "
                   << "color=\"" << color << "\""
                   << "]";
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

bool are_two_complete_dfas_equivalent(const NFA& first_nfa, NFA& second_nfa) {
    u64 alphabet_size = first_nfa.alphabet_size();
    assert(alphabet_size == second_nfa.alphabet_size());

    Worklist_Construction_Context<State_Pair> context;

    {
        State_Pair initial_pair = {.first = first_nfa.initial_states.at(0), .second = second_nfa.initial_states.at(0), .handle = 0};
        context.mark_discovery(initial_pair);
    }

    while (!context.worklist.empty()) {
        auto current_pair = context.extract();

        bool is_final_in_first  = first_nfa.final_states.get_bit_value(current_pair->first);
        bool is_final_in_second = second_nfa.final_states.get_bit_value(current_pair->second);

        if (is_final_in_first != is_final_in_second) {
            return false;
        }

        for (u64 color = 0; color < alphabet_size; color++) {
            auto first_post  = first_nfa.transitions[current_pair->first][color][0];
            auto second_post = second_nfa.transitions[current_pair->second][color][0];

            State_Pair discovery = {.first = first_post, .second = second_post};
            context.mark_discovery(discovery);
        }
    }

    return true;
}

std::ostream& operator<<(std::ostream& os, const State_Pair& state) {
    os << "(" << state.first << ", " << state.second << ", handle=" << state.handle << ")";
    return os;
}


NFA NFA_Builder::build(s64 state_cnt) {
    std::vector<NFA::Transitions_From_State> result_transitions;

    if (state_cnt < 0) {
        for (auto& [state, outgoing_transitions] : this->transitions) {
            for (auto& post_vector : outgoing_transitions) {
                NFA::State post_vector_max = *std::max_element(post_vector.begin(), post_vector.end());
                state_cnt = std::max(state_cnt, static_cast<s64>(post_vector_max));
            }
        }
    }

    result_transitions.resize(state_cnt);

    for (auto& [state, transitions_from_state] : this->transitions) {
        result_transitions[state].resize(this->alphabet_size);

        for (u64 symbol = 0; symbol < this->alphabet_size; symbol++) {
            result_transitions[state][symbol] = std::vector<NFA::State>(transitions_from_state[symbol].begin(), transitions_from_state[symbol].end());
        }
    }

    for (NFA::State state = 0; state < state_cnt; state++) {
        if (result_transitions[state].empty()) {
            result_transitions[state].resize(this->alphabet_size);
        }
    }

    final_states.grow(state_cnt);

    NFA result ({0}, final_states, result_transitions);
    return result;
}

void NFA_Builder::add_transition(NFA::State source_state, u64 symbol, NFA::State destination) {
    auto& transitions_from_state = this->transitions[source_state];
    if (transitions_from_state.empty()) {
        transitions_from_state.resize(this->alphabet_size);
    }

    transitions_from_state[symbol].insert(destination);
}
