#include "swta.hpp"
#include "arith.hpp"
#include "bit_set.hpp"
#include "nfa.hpp"
#include "basics.hpp"

#include <vector>
#include <map>
#include <set>

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
    os << "SWTA {\n";
    for (u64 state = 0; state < swta.transitions.size(); state++) {
        auto& transitions_from_state = swta.transitions[state];

        for (Internal_Symbol sym = 0; sym < swta.number_of_internal_symbols(); sym++) {
            const auto& transitions_along_symbol = transitions_from_state[sym];
                for (Color color = 0; color < swta.number_of_colors(); color++) {
                auto& transition = transitions_along_symbol[color];

                if (!transition.is_present()) continue;

                os << "  q" << state << " --symbol=" << sym << ", color=" << color << "--> [ " << transition << "]\n";
            }
        }
    }
    os << "}";
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

void extend_form_with_product_and_node_discoveries(Linear_Form& destination, const Linear_Form& first, const Linear_Form& second, std::map<State_Pair, u64>& handles, std::vector<State_Pair>& worklist) {
    for (auto& outer_comp : second.components) {
        for (auto& inner_comp : first.components) {

            State_Pair state = {.first = inner_comp.state, .second = outer_comp.state, .handle = handles.size()};

            auto [insert_pos, was_inserted] = handles.emplace(state, state.handle);

            if (was_inserted) {
                worklist.push_back(state);
            } else {
                state.handle = insert_pos->second;
            }

            Algebraic_Complex_Number coef = inner_comp.coef * outer_comp.coef;

            bool already_present = false;
            for (auto& component : destination.components) {
                if (component.state == state.handle) {
                    component.coef += coef;
                    already_present = true;
                    break;
                }
            }

            if (already_present) {
                continue;  // We are done here
            }

            destination.components.push_back({coef, state.handle});
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

        for (Internal_Symbol sym = 0; sym < swta.number_of_internal_symbols(); sym++) {

            auto& transitions_for_sym = transitions_from_state[sym];
            auto& transition          = transitions_for_sym[color];

            if (!transition.is_present()) {
                post.state_set.clear();
                break;
            }

            for (auto& component : transition.left.components) {
                post.state_set.set_bit(component.state);
            }

            for (auto& component : transition.right.components) {
                post.state_set.set_bit(component.state);
            }
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

    NFA result ({0}, final_macrostates, ordered_resulting_transitions);

    do_on_debug({
        result.debug_data = new NFA::Debug_Data;

        for (auto& [macrostate, handle] : handles) {
            result.debug_data->state_names.emplace(handle, macrostate.to_string());
        }
    });

    return result;
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

SWTA::Transition compose_swta_transition_with_wtt(const SWTA::Transition& swta_transition, const WTT::Transition& wtt_transition, std::map<State_Pair, State>& handles, std::vector<State_Pair>& worklist) {
    SWTA::Transition result;

    extend_form_with_product_and_node_discoveries(result.left, wtt_transition.ll, swta_transition.left, handles, worklist);
    extend_form_with_product_and_node_discoveries(result.left, wtt_transition.lr, swta_transition.right, handles, worklist);

    extend_form_with_product_and_node_discoveries(result.right, wtt_transition.rr, swta_transition.right, handles, worklist);
    extend_form_with_product_and_node_discoveries(result.right, wtt_transition.rl, swta_transition.left, handles, worklist);

    return result;
}

SWTA apply_wtt_to_swta(const SWTA& swta, const WTT& wtt) {
    std::vector<State_Pair> worklist; // (SWTA state, WTT state)

    std::map<State_Pair, State> handles;

    Bit_Set leaf_states (0);
    std::vector<State> initial_states;
    SWTA::Transition_Builder transition_builder ( swta.get_metadata() );


    for (auto swta_state : swta.initial_states) {
        for (auto wtt_state : wtt.initial_states) {
            State_Pair product_state = { .first = swta_state, .second = wtt_state, .handle = handles.size() };

            initial_states.push_back(product_state.handle);

            worklist.push_back(product_state);
            handles.emplace(product_state, product_state.handle);
        }
    }

    while (!worklist.empty()) {
        auto& product_state = worklist.back();
        worklist.pop_back();

        bool is_swta_state_leaf = swta.states_with_leaf_transitions.get_bit_value(product_state.first);
        bool is_wtt_state_leaf  = wtt.states_with_leaf_transitions.get_bit_value(product_state.second);

        if (is_swta_state_leaf && is_wtt_state_leaf) {
            leaf_states.grow_and_set_bit(product_state.handle);
        }

        for (Internal_Symbol sym = 0; sym < swta.number_of_internal_symbols(); sym++) {
            const auto& wtt_transition = wtt.transitions[product_state.second][sym];

            for (Color color = 0; color < swta.number_of_colors(); color++) {
                const auto& swta_transition = swta.transitions[product_state.first][sym][color];

                if (!swta_transition.is_present()) continue;
                if (!wtt_transition.is_present())  continue;

                auto result_form = compose_swta_transition_with_wtt(swta_transition, wtt_transition, handles, worklist);
                transition_builder.add_transition(product_state.handle, color, sym, result_form);
            }
        }

    }

    auto transition_fn = transition_builder.build();
    SWTA result (transition_fn, initial_states, leaf_states);
    return result;
}

void put_form_into_matrix(ACN_Matrix& matrix, State source_state, const Linear_Form& form) {
    for (const auto& component : form.components) {
        matrix.set(source_state, component.state, component.coef);
    }
}

std::pair<Affine_Program<Branch_Selector>::Symbol_Handles, Affine_Program<Branch_Selector>::Symbol_Store>
extract_transition_matrices_from_swta(const SWTA& swta) {
    Affine_Program<Branch_Selector>::Symbol_Store symbol_store;
    std::map<Branch_Selector, u64> symbol_handles;

    for (Internal_Symbol internal_sym = 0; internal_sym < swta.number_of_internal_symbols(); internal_sym++) {
        for (Color color = 0; color < swta.number_of_colors(); color++) {
            Branch_Selector left_symbol  (internal_sym, color, Subtree_Tag::LEFT);
            Branch_Selector right_symbol (internal_sym, color, Subtree_Tag::RIGHT);

            ACN_Matrix left_matrix  (swta.number_of_states(), swta.number_of_states());
            ACN_Matrix right_matrix (swta.number_of_states(), swta.number_of_states());

            for (State state = 0; state < swta.number_of_states(); state++) {
                auto& transition = swta.get_transition(state, internal_sym, color);

                put_form_into_matrix(left_matrix,  state, transition.left);
                put_form_into_matrix(right_matrix, state, transition.right);
            }

            const auto [left_insert_pos, left_was_inserted]   = symbol_handles.emplace(left_symbol,  symbol_handles.size());
            const auto [right_insert_pos, right_was_inserted] = symbol_handles.emplace(right_symbol, symbol_handles.size());

            symbol_store.push_back( {left_symbol, left_matrix} );
            symbol_store.push_back( {right_symbol, right_matrix} );
        }
    }

    return {symbol_handles, symbol_store};
}


struct AP_Build_Context {
    const SWTA& swta;
    Affine_Program_Builder<Branch_Selector> builder;
    std::map<Macrostate, u64> macrostate_handles;
    std::vector<Macrostate> worklist;

    AP_Build_Context (const SWTA& swta, const Affine_Program<Branch_Selector>::Symbol_Handles& handles, const Affine_Program<Branch_Selector>::Symbol_Store& store) :
        swta(swta),
        builder(handles, store)
    {}
};

void compute_post_and_store_result_in_ap(const Macrostate& macrostate, const Linear_Form& form, Branch_Selector& symbol_info, AP_Build_Context& context) {
    Macrostate post (context.swta.number_of_states());

    for (const auto& component : form.components) {
        post.state_set.set_bit(component.state);
    }

    post.handle = context.macrostate_handles.size();
    auto [insert_pos, was_inserted] = context.macrostate_handles.emplace(post, post.handle);
    if (was_inserted) { // Brand new macrostate, we need to explore it
        Macrostate& macrostate_in_map = const_cast<Macrostate&>(insert_pos->first);
        macrostate_in_map.populate_state_names_from_set();
        context.worklist.push_back(macrostate_in_map);
    } else {
        post.handle = insert_pos->second;
    }

    u64 symbol_handle = context.builder.symbol_handles.at(symbol_info);
    context.builder.add_transition(macrostate.handle, symbol_handle, post.handle);
}

void extend_affine_program_with_post(const Macrostate& macrostate, Branch_Selector& symbol_info, AP_Build_Context& context) {
    Macrostate right_post (context.swta.number_of_states());

    for (State state : macrostate.state_names) {
        auto& transition = context.swta.get_transition(state, symbol_info.symbol, symbol_info.color);

        if (!transition.is_present()) continue;

        symbol_info.tag = Subtree_Tag::LEFT;
        compute_post_and_store_result_in_ap(macrostate, transition.left, symbol_info, context);

        symbol_info.tag = Subtree_Tag::RIGHT;
        compute_post_and_store_result_in_ap(macrostate, transition.right, symbol_info, context);
    }
}

void rename_linear_form_states(Linear_Form& form, State state_offset) {
    for (auto& component : form.components) {
        component.state += state_offset;
    }
}


SWTA build_difference_swta(const SWTA& first, const SWTA& second) {
    SWTA::Transition_Fn resulting_transitions;
    u64 resulting_state_cnt = 1 + first.number_of_states() + second.number_of_states();
    resulting_transitions.resize(resulting_state_cnt);

    for (State first_state = 0; first_state < first.number_of_states(); first_state++) {
        resulting_transitions[first_state] = first.transitions[first_state];
    }

    u64 offset = first.number_of_states();
    for (State second_state = 0; second_state < second.number_of_states(); second_state++) {
        State new_second_state = second_state + offset;
        resulting_transitions[new_second_state] = second.transitions[second_state];

        for (auto& transition_along_symbol : resulting_transitions[new_second_state]) {
            for (auto& transition_along_color : transition_along_symbol) {
                rename_linear_form_states(transition_along_color.left,  offset);
                rename_linear_form_states(transition_along_color.right, offset);
            }
        }
    }

    State new_initial_state = first.number_of_states() + second.number_of_states();
    { // Add transition (different) from the new initial state
        resulting_transitions[new_initial_state].resize(first.number_of_internal_symbols());

        for (Internal_Symbol symbol = 0; symbol < first.number_of_internal_symbols(); symbol++) {
            resulting_transitions[new_initial_state][symbol].resize(first.number_of_colors());
        }

        Linear_Form::Component first_init_component (Algebraic_Complex_Number::ONE(), first.initial_states[0]);
        Linear_Form::Component second_init_component (-Algebraic_Complex_Number::ONE(), second.initial_states[0] + offset);

        Linear_Form form ({first_init_component, second_init_component});
        resulting_transitions[new_initial_state][0][0] = SWTA::Transition(form, form);
    }

    Bit_Set new_final_states (resulting_state_cnt);
    for (u64 bucket_idx = 0; bucket_idx < first.states_with_leaf_transitions.get_bucket_count(); bucket_idx++) {
        new_final_states.data[bucket_idx] = first.states_with_leaf_transitions.data[bucket_idx];
    }

    for (State second_state = 0; second_state < second.number_of_states(); second_state++) {
        if (second.states_with_leaf_transitions.get_bit_value(second_state)) {
            new_final_states.set_bit(second_state + offset);
        }
    }

    return SWTA(resulting_transitions, {new_initial_state}, new_final_states);
}

Affine_Program<Branch_Selector> build_affine_program(const SWTA& swta) {
    auto [symbol_handles, symbol_store] = extract_transition_matrices_from_swta(swta);
    AP_Build_Context build_context (swta, std::move(symbol_handles), std::move(symbol_store));

    {
        Macrostate init_macrostate (swta.number_of_states(), swta.initial_states);
        init_macrostate.handle = 0;

        build_context.macrostate_handles.emplace(init_macrostate, init_macrostate.handle);
        build_context.worklist.push_back(init_macrostate);

        build_context.builder.mark_state_initial(init_macrostate.handle);
    }

    while (!build_context.worklist.empty()) {
        Macrostate macrostate = build_context.worklist.back();
        build_context.worklist.pop_back();

        if (build_context.swta.states_with_leaf_transitions.is_superset(macrostate.state_set)) {
            build_context.builder.mark_state_final(macrostate.handle);
        }

        for (Internal_Symbol symbol = 0; symbol < swta.number_of_internal_symbols(); symbol++) {
            for (Color color = 0; color < swta.number_of_colors(); color++) {
                Branch_Selector symbol_info(symbol, color, Subtree_Tag::NONE);

                extend_affine_program_with_post(macrostate, symbol_info, build_context);
            }
        }
    }

    Affine_Program program = build_context.builder.build(build_context.macrostate_handles.size());

    do_on_debug({
        program.debug_data = new Affine_Program<Branch_Selector>::Debug_Data;
        for (auto& [macrostate, handle] : build_context.macrostate_handles) {
            program.debug_data->state_names.emplace(handle, macrostate.to_string());
        }
    });

    return program;
}

template <typename T>
void write_affine_program_into_dot(std::ostream& stream, const Affine_Program<T>& program) {
    stream << "digraph Affine_Program {\n";

    stream << "  qInit [shape=none, label=\"\"]\n";
    for (State state = 0; state < program.number_of_states(); state++) {
        stream << "  q" << state;
        if (program.debug_data != nullptr && program.debug_data->state_names.contains(state)) {
            const char* color = program.final_states.get_bit_value(state) ? "green" : "black";
            stream << " [label=\"" << program.debug_data->state_names.at(state) << "\", color=" << color << "]";
        }
        stream << "\n";
    }

    stream << "  qInit -> q" << program.initial_state << "\n";

    for (State state = 0; state < program.number_of_states(); state++) {
        for (u64 symbol = 0; symbol < program.transition_fn[state].size(); symbol++) {
            typename Affine_Program<T>::Transition_Symbol symbol_desc = program.symbol_store.at(symbol);
            for (State destination : program.transition_fn[state][symbol]) {
                stream << "  q" << state << " -> q" << destination << " [label=\"" << symbol_desc.info << "\"]\n";
            }
        }
    }

    stream << "}\n";
}


std::ostream& operator<<(std::ostream& os, const Branch_Selector& info) {
    const char* tag_name = info.tag == Subtree_Tag::LEFT ? "L" : "R";
    // os << "(symbol=" << info.symbol << ", color=" << info.color << ", " << tag_name << ")";
    os << info.symbol << "-" << info.color << "-" << tag_name;
    return os;
}

NFA build_color_language_abstraction(const SWTA& swta) {
    auto abstraction_nfa = build_frontier_automaton(swta);
    auto abstraction_dfa = abstraction_nfa.determinize();
    abstraction_dfa.complete();
    return std::move(abstraction_dfa);
}

bool are_two_swtas_color_equivalent(const SWTA& first, const SWTA& second) {
    NFA first_swta_abstraction  = build_color_language_abstraction(first);
    NFA second_swta_abstraction = build_color_language_abstraction(second);
    bool are_colored_languages_equivalent = are_two_complete_dfas_equivalent(first_swta_abstraction, second_swta_abstraction);

    if (!are_colored_languages_equivalent) return false;

    auto first_program  = build_affine_program(first);
    auto second_program = build_affine_program(second);

    auto product = build_colored_product_of_affine_programs(first_program, second_program, first_swta_abstraction);

    Underlying_SWTA_Info first_swta_info (first);
    Underlying_SWTA_Info second_swta_info (second);
    Underlying_SWTA_Info_Pair swta_info_pair (first_swta_info, second_swta_info);

    bool are_equivalent = !does_affine_program_reach_nonzero_final_states(product, swta_info_pair);

    return are_equivalent;
}

struct State_Delta_Info {
    u32 new_vector_count = 0;  // How many new vectors did we discover for the state
    s32 new_vector_start = 0;  // Where in the matrix is the first newly discovered vector
};

// For debugging
struct Propagation_Info {
    State source;
    State target;
    ACN_Matrix source_row;
    ACN_Matrix propagated_row;
    ACN_Matrix propagated_row_after_insert;
    ACN_Matrix symbol_matrix;
    Branch_Product_Sym symbol_info;
};

struct Affine_Program_Propagation_Context {
    const Affine_Program<Branch_Product_Sym>& program;
    std::vector<State_Delta_Info> state_deltas;
    std::vector<State>            worklist;
    std::vector<ACN_Matrix>       state_vector_spaces;
    u64                           state_space_dimension; 
    std::vector<Propagation_Info> propagation_log;
    s64                           final_state_with_nonzero = -1;
};

struct Propagation_Stats {
    u32 propagation_cnt = 0;
};

std::vector<Propagation_Info*> filter_propagations_for_those_that_affect_state(std::vector<Propagation_Info>& propagations, State state_of_interest) {
    std::set<State> interesting_states;
    interesting_states.insert(state_of_interest);

    std::vector<Propagation_Info*> filtered_propagations;
    
    for (auto it = propagations.rbegin(); it != propagations.rend(); ++it) {
        auto& prop_info = *it;

        if (interesting_states.contains(prop_info.target)) {
            interesting_states.insert(prop_info.source);
            filtered_propagations.push_back(&prop_info);
        }
    }

    return filtered_propagations;
}

void write_propagation_info(std::map<State, std::string> state_names, Propagation_Info& info) {
    std::cout << "Propagated from "
              << info.source << " aka " << state_names.at(info.source)
              << " to " << info.target << " aka " << state_names.at(info.target)
              <<" along " << info.symbol_info << "\n";
    std::cout << "  row entering pipe : " << info.source_row << "\n";
    std::cout << "  row exiting pipe(0):" << info.propagated_row << "\n";
    std::cout << "  row exiting pipe(1):" << info.propagated_row_after_insert << "\n";
    // std::cout << "  pipe matrix:" << propagation_ptr->symbol_matrix << "\n";
}

void dump_propagations(std::vector<Propagation_Info*> propagations, Affine_Program_Propagation_Context& context) {
    std::map<State, std::string>& state_names = context.program.debug_data->state_names;
    for (auto propagation_ptr : propagations) {
        write_propagation_info(state_names, *propagation_ptr);
    }
}

void dump_propagations(std::vector<Propagation_Info> propagations, Affine_Program_Propagation_Context& context) {
    std::map<State, std::string>& state_names = context.program.debug_data->state_names;
    for (auto propagation_ptr : propagations) {
        write_propagation_info(state_names, propagation_ptr);
    }
}

void propagate_vector_from_state_matrix_to_successors(Affine_Program_Propagation_Context& context, State_Delta_Info& old_state_info, State current_state, u64 symbol) {
    auto& successors_along_symbol = context.program.transition_fn[current_state][symbol];

    u64 last_row_idx_exc = old_state_info.new_vector_start + old_state_info.new_vector_count;

    auto& current_state_matrix  = context.state_vector_spaces[current_state];

    const auto& transition_matrix = context.program.symbol_store[symbol].matrix;

    for (State successor : successors_along_symbol) {
        auto& successor_delta = context.state_deltas[successor];
        if (successor_delta.new_vector_start >= context.state_space_dimension) {
            continue;  // The matrix for the successor is already saturated 
        }

        for (u64 row_idx = old_state_info.new_vector_start; row_idx < last_row_idx_exc; row_idx++) {
            auto& target_matrix = context.state_vector_spaces[successor];

            auto row = current_state_matrix.extract_row(row_idx); // @Optimize: maybe we do not need to do a copy here since we already will be allocating during matrix multiply
            auto propagated_row = row * transition_matrix;

            if (true) {
                Propagation_Info propagation_info = {
                    .source = current_state,
                    .target = successor,
                    .source_row = row,
                    .propagated_row = row * transition_matrix,
                    .propagated_row_after_insert = propagated_row,
                    .symbol_matrix = transition_matrix,
                    .symbol_info = context.program.symbol_store[symbol].info
                };
                write_propagation_info(context.program.debug_data->state_names, propagation_info);
                context.propagation_log.push_back(propagation_info);
            };

            if (propagated_row.contains_only_zeros()) continue;

            s64 new_row_position = add_row_to_row_echelon_matrix_no_copy(target_matrix, propagated_row);

            if (new_row_position < 0) {
                continue;
            }

            if (new_row_position == context.state_space_dimension) {
                successor_delta.new_vector_start = context.program.number_of_states();
                break;  // The matrix became satured, no need to try propagate anything
            }

            if (successor_delta.new_vector_count == 0) {
                context.worklist.push_back(successor);
            }

            successor_delta.new_vector_count += 1;
        }
    }
}

bool does_state_vector_space_contain_nonzeros(Affine_Program_Propagation_Context& context, const ACN_Matrix& final_vector, State state) {
    auto& current_state_matrix  = context.state_vector_spaces[state];
    auto pontential_leaf_values = current_state_matrix * final_vector; // [n x 1] vector

    if (!pontential_leaf_values.contains_only_zeros()) { // A final state is reachable with non-zero value
        context.final_state_with_nonzero = state;
        return true;
    };

    return false;
}

bool does_affine_program_reach_nonzero_final_states(const Affine_Program<Branch_Product_Sym>& program, const Underlying_SWTA_Info_Pair& swta_pair_info) {
    u64 swta_state_cnt = swta_pair_info.first_swta_info.state_cnt + swta_pair_info.second_swta_info.state_cnt;

    ACN_Matrix final_vector (swta_state_cnt, 1);
    for (State final_state : swta_pair_info.first_swta_info.final_states) {
        Algebraic_Complex_Number acn_one = Algebraic_Complex_Number::ONE();
        final_vector.set(final_state, 0, acn_one);
    }
    for (State final_state : swta_pair_info.second_swta_info.final_states) {
        Algebraic_Complex_Number acn_one = Algebraic_Complex_Number::ONE();
        final_vector.set(final_state + swta_pair_info.first_swta_info.state_cnt, 0, acn_one);
    }
    
    Propagation_Stats stats;
    Affine_Program_Propagation_Context context(program);
    {  // Initialize the context
        context.state_space_dimension = swta_state_cnt;
        context.state_vector_spaces.reserve(program.number_of_states());

        for (State state = 0; state < program.number_of_states(); state++) {
            context.state_vector_spaces.push_back(ACN_Matrix(context.state_space_dimension, context.state_space_dimension));
        }
        context.state_deltas.resize(program.number_of_states());
    }

    { // Push the initial state with the initial vector into the worklist
        context.worklist.push_back(program.initial_state);
        context.state_deltas[program.initial_state].new_vector_count = 1;

        auto& initial_state_matrix = context.state_vector_spaces[program.initial_state];
        Algebraic_Complex_Number init_vector_value = Algebraic_Complex_Number::ONE();
        initial_state_matrix.set(0, swta_pair_info.first_swta_info.initial_states[0], init_vector_value);
        initial_state_matrix.set(0, swta_pair_info.second_swta_info.initial_states[0], -init_vector_value);
    }

    while (!context.worklist.empty()) {
        stats.propagation_cnt += 1;
        do_on_debug({
            std::cout << "Processing to " <<  context.worklist.back() << " state matrix:\n " <<  context.state_vector_spaces[context.worklist.back()] << "\n";
        });
        
        State current_state = context.worklist.back();
        State_Delta_Info old_delta_info = context.state_deltas[current_state];
        context.worklist.pop_back();

        // Mark in state deltas that we have already processe the new vectors
        context.state_deltas[current_state].new_vector_start += old_delta_info.new_vector_count;
        context.state_deltas[current_state].new_vector_count = 0;

        if (program.final_states.get_bit_value(current_state)) {
            bool contains_nonzeros = does_state_vector_space_contain_nonzeros(context, final_vector, current_state);
            if (contains_nonzeros) {
                goto terminate;
            }
        }

        for (u64 symbol = 0; symbol < program.number_of_symbols(); symbol++) {
            propagate_vector_from_state_matrix_to_successors(context, old_delta_info, current_state, symbol);
        }
    }
    
    for (State state = 0; state < program.number_of_states(); state++) {
        if (program.final_states.get_bit_value(state)) {
            bool contains_nonzeros = does_state_vector_space_contain_nonzeros(context, final_vector, state);
            if (contains_nonzeros) {
                goto terminate;
            }
        }
    }

    terminate:
        bool reaches_nonzero = context.final_state_with_nonzero >= 0;

        do_on_debug({std::cout << "Performed " << context.propagation_log.size() << " propagations.\n"; });

        dump_propagations(context.propagation_log, context);
        if (reaches_nonzero) {
            auto filtered_propagations = filter_propagations_for_those_that_affect_state(context.propagation_log, context.final_state_with_nonzero);
            dump_propagations(filtered_propagations, context);
        }

        return reaches_nonzero;
}

struct Colored_Product_State {
    State in_first;
    State in_second;
    State in_frontier;

    s64 handle = -1;

    bool operator<(const Colored_Product_State& other) const {
        if (this->in_first < other.in_first) return true;
        if (this->in_first > other.in_first) return false;

        if (this->in_second < other.in_second) return true;
        if (this->in_second > other.in_second) return false;

        return this->in_frontier < other.in_frontier;
    }
};

ACN_Matrix compose_transition_matrices(const ACN_Matrix& first_matrix, const ACN_Matrix& second_matrix) {
    u64 width  = first_matrix.width + second_matrix.width;
    u64 height = first_matrix.height + second_matrix.height;
    ACN_Matrix result(height, width);

    for (u64 first_matrix_row_idx = 0; first_matrix_row_idx < first_matrix.height; first_matrix_row_idx++) {
        for (u64 first_matrix_col_idx = 0; first_matrix_col_idx < first_matrix.width; first_matrix_col_idx++) {
            auto& elem = first_matrix.at(first_matrix_row_idx, first_matrix_col_idx);
            result.set(first_matrix_row_idx, first_matrix_col_idx, elem);
        }
    }

    for (u64 second_matrix_row_idx = 0; second_matrix_row_idx < second_matrix.height; second_matrix_row_idx++) {
        for (u64 second_matrix_col_idx = 0; second_matrix_col_idx < second_matrix.width; second_matrix_col_idx++) {
            auto& elem = second_matrix.at(second_matrix_row_idx, second_matrix_col_idx);
            result.set(second_matrix_row_idx + first_matrix.height, second_matrix_col_idx + second_matrix.width, elem);
        }
    }

    return result;
}

std::vector<std::vector<const Affine_Program<Branch_Selector>::Transition_Symbol*>>
group_transition_symbols_by_color(const Affine_Program<Branch_Selector>::Symbol_Store& symbol_store, u64 color_cnt) {
    std::vector<std::vector<const Affine_Program<Branch_Selector>::Transition_Symbol*>> result;
    result.resize(color_cnt);

    for (auto& symbol : symbol_store) {
        result[symbol.info.color].push_back(&symbol);
    }

    return result;
}

struct Product_Alphabet_Prep {
    std::vector<std::vector<const Affine_Program<Branch_Selector>::Transition_Symbol*>> first_transition_syms_by_color;
    std::vector<std::vector<const Affine_Program<Branch_Selector>::Transition_Symbol*>> second_transition_syms_by_color;
    Affine_Program<Branch_Product_Sym>::Symbol_Handles symbol_handles;
    Affine_Program<Branch_Product_Sym>::Symbol_Store   symbol_store;
};

Product_Alphabet_Prep prepare_transition_symbols_for_product_automaton(const Affine_Program<Branch_Selector>& first_ap, const Affine_Program<Branch_Selector>& second_ap, u64 color_cnt) {
    auto selectors_in_first_by_color  = group_transition_symbols_by_color(first_ap.symbol_store, color_cnt);
    auto selectors_in_second_by_color = group_transition_symbols_by_color(second_ap.symbol_store, color_cnt);

    Affine_Program<Branch_Product_Sym>::Symbol_Handles handles;
    Affine_Program<Branch_Product_Sym>::Symbol_Store   store;

    for (Color color = 0; color < color_cnt; color++) {
        for (auto& first_ap_symbol : selectors_in_first_by_color[color]) {
            for (auto& second_ap_symbol : selectors_in_second_by_color[color]) {
                Branch_Product_Sym sym = {
                    .color = color,
                    .first_sym = first_ap_symbol->info.symbol, .first_tag = first_ap_symbol->info.tag,
                    .second_sym = second_ap_symbol->info.symbol, .second_tag = second_ap_symbol->info.tag
                };

                handles.emplace(sym, handles.size());

                ACN_Matrix matrix = compose_transition_matrices(first_ap_symbol->matrix, second_ap_symbol->matrix);
                Affine_Program<Branch_Product_Sym>::Transition_Symbol symbol_to_store = {.info = sym, .matrix = matrix};
                store.push_back(symbol_to_store);
            }
        }
    }

    Product_Alphabet_Prep result(selectors_in_first_by_color, selectors_in_second_by_color, handles, store);
    return result;
}


struct Colored_Product_Build_Context {
    Worklist_Construction_Context<Colored_Product_State> worklist_info;
    const Affine_Program<Branch_Selector>& first_ap;
    const Affine_Program<Branch_Selector>& second_ap;
    const NFA& frontier;
    Affine_Program_Builder<Branch_Product_Sym> builder;

    Colored_Product_Build_Context(
        const Affine_Program<Branch_Selector>& first_ap,
        const Affine_Program<Branch_Selector>& second_ap,
        const NFA& frontier,
        const Affine_Program<Branch_Product_Sym>::Symbol_Handles& handles,
        const Affine_Program<Branch_Product_Sym>::Symbol_Store& store) :
        first_ap(first_ap),
        second_ap(second_ap),
        frontier(frontier),
        builder(handles, store) {}
};

void product_build_step_along_color(Colored_Product_Build_Context& context, Product_Alphabet_Prep& alphabet_prep, const Colored_Product_State* current_state, Color color, State frontier_nfa_state) {
    for (auto first_branch_selector_ptr : alphabet_prep.first_transition_syms_by_color[color]) {
        for (auto second_branch_selector_ptr : alphabet_prep.second_transition_syms_by_color[color]) {
            Branch_Product_Sym product_symbol = {
                .color = color,
                .first_sym = first_branch_selector_ptr->info.symbol, .first_tag = first_branch_selector_ptr->info.tag,
                .second_sym = second_branch_selector_ptr->info.symbol, .second_tag = second_branch_selector_ptr->info.tag
            };
            u64 product_symbol_handle = context.builder.symbol_handles.at(product_symbol);

            u64 first_branch_selector_handle  = context.first_ap.symbol_handles.at(first_branch_selector_ptr->info);
            u64 second_branch_selector_handle = context.second_ap.symbol_handles.at(second_branch_selector_ptr->info);

            auto& first_ap_successors  = context.first_ap.transition_fn[current_state->handle][first_branch_selector_handle];
            auto& second_ap_successors = context.second_ap.transition_fn[current_state->handle][second_branch_selector_handle];

            for (auto first_successor : first_ap_successors) {
                for (auto second_successor : second_ap_successors) {
                    Colored_Product_State imm_product_state = {.in_first = first_successor, .in_second = second_successor, .in_frontier = frontier_nfa_state};
                    auto dest_state = context.worklist_info.mark_discovery(imm_product_state);

                    context.builder.add_transition(current_state->handle, product_symbol_handle, dest_state->handle);
                }
            }
        }
    }
}

Affine_Program<Branch_Product_Sym>
build_colored_product_of_affine_programs(const Affine_Program<Branch_Selector>& first_ap, const Affine_Program<Branch_Selector>& second_ap, const NFA& frontier) {
    u64 color_cnt = frontier.alphabet_size();
    auto alphabet_prep = prepare_transition_symbols_for_product_automaton(first_ap, second_ap, color_cnt);

    Colored_Product_Build_Context context (first_ap, second_ap, frontier, alphabet_prep.symbol_handles, alphabet_prep.symbol_store);
    {
        Colored_Product_State imm_initial_state {.in_first = first_ap.initial_state, .in_second = second_ap.initial_state, .in_frontier = frontier.initial_states[0]};
        auto initial_state = context.worklist_info.mark_discovery(imm_initial_state);        
        context.builder.mark_state_initial(initial_state->handle);
    }

    while (!context.worklist_info.worklist.empty()) {
        auto current_state = context.worklist_info.extract();

        bool first_accepts    = context.first_ap.final_states.get_bit_value(current_state->in_first);
        bool second_accepts   = context.second_ap.final_states.get_bit_value(current_state->in_second);
        bool frontier_accepts = context.frontier.final_states.get_bit_value(current_state->in_frontier);

        if (first_accepts && second_accepts && frontier_accepts) {
            context.builder.mark_state_final(current_state->handle);
        };

        for (Color color = 0; color < color_cnt; color++) {
            State nfa_post_state = context.frontier.transitions[current_state->in_frontier][color].at(0);
            product_build_step_along_color(context, alphabet_prep, current_state, color, nfa_post_state);
        }
    }

    auto result = context.builder.build();
    return result;
}


std::ostream& operator<<(std::ostream& stream, const Branch_Product_Sym& sym) {
    const char* first_sym_str  = sym.first_tag == Subtree_Tag::LEFT ? "L" : "R";
    const char* second_sym_str = sym.second_tag == Subtree_Tag::LEFT ? "L" : "R";

    stream << "<c=" << sym.color << "-" << sym.first_sym << first_sym_str << "-" << sym.second_sym << second_sym_str << ">";

    return stream;    
}
