#include "swta.hpp"
#include "arith.hpp"

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

    for (u64 state = 0; state < wtt.number_of_states(); state++) {
        os << "  " << state << ": " << wtt.transitions[state] << "\n";
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
    std::map<State, WTT::Transition> transitions;

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

        WTT::Transition first_transition   = first.transitions[state_pair.first];
        WTT::Transition second_transitions = second.transitions[state_pair.second];

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
        transitions.emplace(handle, resulting_transition);
    }

    std::vector<WTT::Transition> ordered_transitions;
    ordered_transitions.reserve(transitions.size());
    for (auto& [state, transition] : transitions) {
        ordered_transitions.push_back(transition);
    }

    WTT result (ordered_transitions, {}, initial_states);
    return result;
}
