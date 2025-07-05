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
    std::cout << "  initial states: " << swta.initial_states << "\n";
    std::cout << "  final  states:  " << swta.states_with_leaf_transitions.into_vector() << "\n";
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

    for (u64 state = 0; state < wtt.number_of_states(); state++) {
        const std::vector<WTT::Transition>& transitions_from_state = wtt.transitions[state];
        for (Internal_Symbol sym = 0; sym < wtt.number_of_internal_symbols(); sym++) {
            auto& transition = transitions_from_state[sym];
            if (!transition.is_present()) continue;
            os << "  " << state << "--(sym=" << sym << ")-->: " << transitions_from_state[sym] << "\n";
        }
    }

    os << "}";
    return os;
}

void extend_form_with_product_and_node_discoveries(Linear_Form& destination, const Linear_Form& first, const Linear_Form& second, Worklist_Construction_Context<State_Pair>& worklist_state) {
    for (auto& outer_comp : second.components) {
        for (auto& inner_comp : first.components) {

            const State_Pair* imm_state = nullptr;
            {
                State_Pair state = { .first = inner_comp.state, .second = outer_comp.state };
                imm_state = worklist_state.mark_discovery(state);   
            }

            Algebraic_Complex_Number coef = inner_comp.coef * outer_comp.coef;

            bool already_present = false;
            for (auto& component : destination.components) {
                if (component.state == imm_state->handle) {
                    component.coef += coef;
                    already_present = true;
                    break;
                }
            }

            if (already_present) {
                continue;  // We are done here
            }

            Linear_Form::Component new_component(coef, imm_state->handle);
            destination.components.push_back(new_component);
        }
    }
}


WTT compose_wtts_sequentially(WTT& first, WTT& second) {

    const u64 num_of_internal_symbols = first.number_of_internal_symbols();
    assert(num_of_internal_symbols == second.number_of_internal_symbols());

    SWTA::Metadata metadata = {.number_of_internal_symbols = num_of_internal_symbols, .number_of_colors = 1};
    WTT_Builder builder(metadata);

    Worklist_Construction_Context<State_Pair> worklist_state;

    State_Pair initial_state = { .first = first.initial_states[0], .second = second.initial_states[0] };
    auto imm_initial_state = worklist_state.mark_discovery(initial_state);
    builder.mark_state_initial(imm_initial_state->handle);

    while (worklist_state.has_more_to_explore()) {
        auto state_pair = worklist_state.extract();

        auto is_first_leaf  = first.states_with_leaf_transitions.get_bit_value(state_pair->first);
        auto is_second_leaf = second.states_with_leaf_transitions.get_bit_value(state_pair->second);

        if (is_first_leaf && is_second_leaf) {
            builder.mark_state_final(state_pair->handle);
        }

        for (Internal_Symbol internal_symbol = 0; internal_symbol < num_of_internal_symbols; internal_symbol++) {
            WTT::Transition& first_transition   = first.transitions[state_pair->first][internal_symbol];
            WTT::Transition& second_transitions = second.transitions[state_pair->second][internal_symbol];

            Linear_Form ll, lr, rl, rr;
            { // ll
                extend_form_with_product_and_node_discoveries(ll, first_transition.ll, second_transitions.ll, worklist_state);
                extend_form_with_product_and_node_discoveries(ll, first_transition.rl, second_transitions.lr, worklist_state);
            }

            { // lr
                extend_form_with_product_and_node_discoveries(lr, first_transition.lr, second_transitions.ll, worklist_state);
                extend_form_with_product_and_node_discoveries(lr, first_transition.rr, second_transitions.lr, worklist_state);
            }

            { // rl
                extend_form_with_product_and_node_discoveries(rl, first_transition.ll, second_transitions.rl, worklist_state);
                extend_form_with_product_and_node_discoveries(rl, first_transition.rl, second_transitions.rr, worklist_state);
            }

            { // rr
                extend_form_with_product_and_node_discoveries(rr, first_transition.lr, second_transitions.rl, worklist_state);
                extend_form_with_product_and_node_discoveries(rr, first_transition.rr, second_transitions.rr, worklist_state);
            }

            WTT::Transition resulting_transition (ll, lr, rl, rr);

            builder.add_transition(state_pair->handle, internal_symbol, resulting_transition);
        }
    }

    WTT result = builder.build(worklist_state.handles.size());
    return result;
}



Macrostate compute_post(const Macrostate* macrostate, const SWTA& swta, Color color, Internal_Symbol symbol) {
    Macrostate post(swta.number_of_states());

    for (State state : macrostate->state_names) {
        auto& transitions_from_state = swta.transitions[state];
        auto& transitions_for_sym    = transitions_from_state[symbol];
        auto& transition             = transitions_for_sym[color];

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

    post.populate_state_names_from_set();

    return post;
}

Macrostate compute_post(const Macrostate* macrostate, const WTT& wtt, Color color, Internal_Symbol symbol) {
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

    post.populate_state_names_from_set();

    return post;
}

void dump_discovered_transitions(const std::map<State, std::vector<std::vector<State>>>& transitions) {
    std::cout << "Known states: ";
    for (const auto& [state, transitions_from_state]: transitions) {
        std::cout << state << ", ";
    }
    std::cout << "\n";
}

void initialize_frontier_with_initial_states(Worklist_Construction_Context<Macrostate>& worklist_state, NFA_Builder& builder, const std::vector<State>& initial_states, u64 total_number_of_states, s64 root) {
    if (root < 0) {
        Macrostate initial_macrostate (total_number_of_states, initial_states);
        worklist_state.mark_discovery(initial_macrostate);
    } else { // Start the construction from the provided root
        std::vector<State> root_states ({static_cast<State>(root)});
        Macrostate initial_macrostate (total_number_of_states, root_states);
        worklist_state.mark_discovery(initial_macrostate);
    }

    builder.mark_state_initial(0);
}

template <typename Tree_Transition_System>
NFA build_frontier_automaton(const Tree_Transition_System& tts, s64 root) {

    Worklist_Construction_Context<Macrostate> worklist_state;
    NFA_Builder builder(tts.number_of_colors());

    initialize_frontier_with_initial_states(worklist_state, builder, tts.initial_states, tts.number_of_states(), root);

    u64 color_cnt = tts.number_of_colors();

    while (worklist_state.has_more_to_explore()) {
        auto macrostate = worklist_state.extract();

        if (tts.states_with_leaf_transitions.is_superset(macrostate->state_set)) { // All of the states in macrostate can make a leaf transition
            builder.mark_state_final(macrostate->handle);
        }

        for (Internal_Symbol internal_symbol = 0; internal_symbol < tts.number_of_internal_symbols(); internal_symbol++) {
            for (Color color = 0; color < color_cnt; color++) {
                Macrostate imm_post = compute_post(macrostate, tts, color, internal_symbol);

                if (imm_post.empty()) {
                    continue;
                }

                auto post = worklist_state.mark_discovery(imm_post);
                builder.add_transition(macrostate->handle, color, post->handle);
            }
        }
    }

    NFA result = builder.build(worklist_state.handles.size());

    do_on_debug({
        result.debug_data = new NFA::Debug_Data;

        for (auto& [macrostate, handle] : worklist_state.handles) {
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

SWTA::Transition compose_swta_transition_with_wtt(const SWTA::Transition& swta_transition, const WTT::Transition& wtt_transition, Worklist_Construction_Context<State_Pair>& worklist_state) {
    SWTA::Transition result;

    extend_form_with_product_and_node_discoveries(result.left, swta_transition.left, wtt_transition.ll, worklist_state);
    extend_form_with_product_and_node_discoveries(result.left, swta_transition.right, wtt_transition.lr, worklist_state);

    extend_form_with_product_and_node_discoveries(result.right, swta_transition.right, wtt_transition.rr, worklist_state);
    extend_form_with_product_and_node_discoveries(result.right, swta_transition.left, wtt_transition.rl, worklist_state);

    return result;
}

std::ostream& operator<<(std::ostream& out, const SWTA::Transition_Builder& builder) {
    std::cout << "Builder (stored transitions): ";
    for (auto& [handle, transitions] : builder.transitions) {
        out << transitions << "\n";
    }

    return out;
}

SWTA apply_wtt_to_swta(const SWTA& swta, const WTT& wtt) {
    Worklist_Construction_Context<State_Pair> worklist_state;

    Bit_Set leaf_states (0);
    std::vector<State> initial_states;
    SWTA::Transition_Builder transition_builder ( swta.get_metadata() );


    State_Pair initial_state = { .first = swta.initial_states[0], .second = wtt.initial_states[0] };
    auto imm_initial_state = worklist_state.mark_discovery(initial_state);
    initial_states.push_back(imm_initial_state->handle);
   
    while (worklist_state.has_more_to_explore()) {
        auto product_state = worklist_state.extract();

        bool is_swta_state_leaf = swta.states_with_leaf_transitions.get_bit_value(product_state->first);
        bool is_wtt_state_leaf  = wtt.states_with_leaf_transitions.get_bit_value(product_state->second);

        if (is_swta_state_leaf && is_wtt_state_leaf) {
            leaf_states.grow_and_set_bit(product_state->handle);
        }

        for (Internal_Symbol sym = 0; sym < swta.number_of_internal_symbols(); sym++) {
            const auto& wtt_transition = wtt.transitions[product_state->second][sym];

            for (Color color = 0; color < swta.number_of_colors(); color++) {
                const auto& swta_transition = swta.transitions[product_state->first][sym][color];

                if (!swta_transition.is_present()) continue;
                if (!wtt_transition.is_present())  continue;

                auto result_form = compose_swta_transition_with_wtt(swta_transition, wtt_transition, worklist_state);
                transition_builder.add_transition(product_state->handle, sym, color, result_form);
            }
        }

    }

    auto transition_fn = transition_builder.build(worklist_state.handles.size());
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


struct AP_State_Info {
    Macrostate macrostate; // What states are active in a branch
    State color_sym_abstraction_state;

    AP_State_Info(u64 swta_state_cnt, State color_sym_abstr_state) :
        macrostate(swta_state_cnt),
        color_sym_abstraction_state(color_sym_abstr_state) {}

    AP_State_Info(const Macrostate& macrostate, State color_sym_abstr_state) :
        macrostate(macrostate),
        color_sym_abstraction_state(color_sym_abstr_state) {}

    bool operator<(const AP_State_Info& other) const {
        INSERT_LEX_LT_CODE(color_sym_abstraction_state, other.color_sym_abstraction_state);
        return this->macrostate < other.macrostate;
    }
};

struct AP_Build_Context {
    const SWTA& swta;
    const Color_Symbol_Abstraction               color_sym_abstraction;
    Affine_Program_Builder<Branch_Selector>      builder;
    Worklist_Construction_Context<AP_State_Info> worklist_state;

    AP_Build_Context (const SWTA& swta, const Color_Symbol_Abstraction& abstraction, const Affine_Program<Branch_Selector>::Symbol_Handles& handles, const Affine_Program<Branch_Selector>::Symbol_Store& store) :
        swta(swta),
        color_sym_abstraction(abstraction),
        builder(handles, store)
    {}
};

void add_transition_and_add_to_worklist_if_new(const AP_State_Info& source, AP_State_Info& post, Branch_Selector& symbol_info, AP_Build_Context& context) {

    post.macrostate.handle = context.worklist_state.handles.size();
    auto [insert_pos, was_inserted] = context.worklist_state.handles.emplace(post, post.macrostate.handle);

    if (was_inserted) {
        context.worklist_state.worklist.push_back(&insert_pos->first);
        const_cast<AP_State_Info*>(&insert_pos->first)->macrostate.populate_state_names_from_set();
    } else {
        post.macrostate.handle = insert_pos->second;
    }

    u64 symbol_handle = context.builder.symbol_handles.at(symbol_info);
    context.builder.add_transition(source.macrostate.handle, symbol_handle, post.macrostate.handle);
}

void extend_affine_program_with_post(const AP_State_Info& source_ap_state, Branch_Selector& symbol_info, AP_Build_Context& context) {
    for (State state : source_ap_state.macrostate.state_names) {
        auto& transition = context.swta.get_transition(state, symbol_info.symbol, symbol_info.color);

        if (!transition.is_present()) return; // There is no post to compute, all of the states have to have a transition along the symbol
    }

    Color_Symbol color_symbol = {.color = symbol_info.color, .symbol = symbol_info.symbol};
    u64 abstraction_color_sym_handle = context.color_sym_abstraction.symbol_handles.at(color_symbol);
    auto& abstraction_successors = context.color_sym_abstraction.abstraction.transitions[source_ap_state.color_sym_abstraction_state][abstraction_color_sym_handle];
    if (abstraction_successors.empty()) return;
    State abstraction_post = abstraction_successors[0];

    AP_State_Info left_post (context.swta.number_of_states(), abstraction_post);
    AP_State_Info right_post (context.swta.number_of_states(), abstraction_post);

    for (State state : source_ap_state.macrostate.state_names) {
        auto& transition = context.swta.get_transition(state, symbol_info.symbol, symbol_info.color);

        for (auto& component : transition.left.components) {
            left_post.macrostate.state_set.set_bit(component.state);
        }

        for (auto& component : transition.right.components) {
            right_post.macrostate.state_set.set_bit(component.state);
        }
    }

    symbol_info.tag = Subtree_Tag::LEFT;
    add_transition_and_add_to_worklist_if_new(source_ap_state, left_post, symbol_info, context);

    symbol_info.tag = Subtree_Tag::RIGHT;
    add_transition_and_add_to_worklist_if_new(source_ap_state, right_post, symbol_info, context);
}

Affine_Program<Branch_Selector> build_affine_program(const SWTA& swta, const Color_Symbol_Abstraction& color_sym_abstraction) {
    auto [symbol_handles, symbol_store] = extract_transition_matrices_from_swta(swta);
    AP_Build_Context build_context (swta, color_sym_abstraction, std::move(symbol_handles), std::move(symbol_store));

    {
        Macrostate init_macrostate (swta.number_of_states(), swta.initial_states);
        init_macrostate.handle = 0;
        AP_State_Info initial_ap_state (init_macrostate, build_context.color_sym_abstraction.abstraction.initial_states[0]);
        auto [insert_pos, was_inserted] = build_context.worklist_state.handles.emplace(initial_ap_state, 0);

        build_context.worklist_state.worklist.push_back(&insert_pos->first);
        build_context.builder.mark_state_initial(0);
    }

    while (build_context.worklist_state.has_more_to_explore()) {
        auto current_ap_state = build_context.worklist_state.extract();

        bool swta_branch_contains_only_leaves = build_context.swta.states_with_leaf_transitions.is_superset(current_ap_state->macrostate.state_set);
        bool swta_accepts_in_remaining_branches = build_context.color_sym_abstraction.abstraction.final_states.get_bit_value(current_ap_state->color_sym_abstraction_state);

        if (swta_accepts_in_remaining_branches && swta_accepts_in_remaining_branches) {
            build_context.builder.mark_state_final(current_ap_state->macrostate.handle);
        }

        for (Internal_Symbol symbol = 0; symbol < swta.number_of_internal_symbols(); symbol++) {
            for (Color color = 0; color < swta.number_of_colors(); color++) {
                Branch_Selector symbol_info(symbol, color, Subtree_Tag::NONE);

                extend_affine_program_with_post(*current_ap_state, symbol_info, build_context);
            }
        }
    }

    Affine_Program program = build_context.builder.build(build_context.worklist_state.handles.size());

    do_on_debug({
        program.debug_data = new Affine_Program<Branch_Selector>::Debug_Data;
        for (auto& [ap_state, handle] : build_context.worklist_state.handles) {
            program.debug_data->state_names.emplace(handle, ap_state.macrostate.to_string());
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
    // abstraction_nfa.write_dot(std::cout);
    auto abstraction_dfa = abstraction_nfa.determinize();
    abstraction_dfa.complete();
    return std::move(abstraction_dfa);
}

bool are_two_swtas_color_equivalent(const SWTA& first, const SWTA& second) {
    NFA first_swta_abstraction  = build_color_language_abstraction(first);

    NFA second_swta_abstraction = build_color_language_abstraction(second);
    bool are_colored_languages_equivalent = are_two_complete_dfas_equivalent(first_swta_abstraction, second_swta_abstraction);

    if (!are_colored_languages_equivalent) {
        do_on_debug({
            std::cout << "Provided SWTAs do not have equal colored languages:" << "\n";
            first_swta_abstraction.write_dot(std::cout);
            std::cout << "\n------------------------" << "\n";
            second_swta_abstraction.write_dot(std::cout);
        });
        return false;
    }

    auto first_color_sym_abstraction = build_color_internal_symbol_abstraction(first);
    auto first_program  = build_affine_program(first, first_color_sym_abstraction);

    auto second_color_sym_abstraction = build_color_internal_symbol_abstraction(second);
    auto second_program = build_affine_program(second, second_color_sym_abstraction);

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
    ACN_Matrix                    final_vector;
    std::vector<State_Delta_Info> state_deltas;
    std::vector<State>            worklist;
    std::vector<ACN_Matrix>       state_vector_spaces;
    u64                           state_space_dimension;
    std::vector<Propagation_Info> propagation_log;
    s64                           final_state_with_nonzero = -1;
    bool                          check_final_state_early = true;
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

void write_propagation_info(std::map<State, std::string> state_names, Propagation_Info& info, bool write_state_names = false) {
    std::cout << "Propagated from "
              << info.source;
    if (write_state_names) {
        std::cout << " aka " << state_names.at(info.source);
    }
    std::cout << " to " << info.target;
    if (write_state_names) {
        std::cout << " aka " << state_names.at(info.target);
    }
    std::cout <<" along " << info.symbol_info << "\n";
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

bool does_state_vector_space_contain_nonzeros(Affine_Program_Propagation_Context& context, const ACN_Matrix& final_vector, State state) {
    auto& current_state_matrix  = context.state_vector_spaces[state];
    auto pontential_leaf_values = current_state_matrix * final_vector; // [n x 1] vector

    if (!pontential_leaf_values.contains_only_zeros()) { // A final state is reachable with non-zero value
        context.final_state_with_nonzero = state;
        return true;
    };

    return false;
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

            if (propagated_row.contains_only_zeros()) continue;

            s64 new_row_position = add_row_to_row_echelon_matrix_no_copy(target_matrix, propagated_row);

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

                if (context.check_final_state_early) {
                     if (does_state_vector_space_contain_nonzeros(context, context.final_vector, current_state)) {
                         return;
                     }
                }
                // write_propagation_info(context.program.debug_data->state_names, propagation_info);
                context.propagation_log.push_back(propagation_info);
            };

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
    Affine_Program_Propagation_Context context(program, final_vector);
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

        State second_initial_state_with_offset = swta_pair_info.second_swta_info.initial_states[0] + swta_pair_info.first_swta_info.state_cnt;
        initial_state_matrix.set(0, second_initial_state_with_offset, -init_vector_value);
    }

    while (!context.worklist.empty()) {
        stats.propagation_cnt += 1;
        // do_on_debug({
            // std::cout << "Processing to " <<  context.worklist.back() << " state matrix:\n " <<  context.state_vector_spaces[context.worklist.back()] << "\n";
        // });

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
            if (context.final_state_with_nonzero > -1) goto terminate;
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
            result.set(second_matrix_row_idx + first_matrix.height, second_matrix_col_idx + first_matrix.width, elem);
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
                if (first_ap_symbol->info.tag != second_ap_symbol->info.tag) continue;

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
            if (first_branch_selector_ptr->info.tag != second_branch_selector_ptr->info.tag) {
                // We have to go along the same branches in the tree, otherwise we would be comparing leaves from completely different paths
                continue;
            }

            Branch_Product_Sym product_symbol = {
                .color = color,
                .first_sym = first_branch_selector_ptr->info.symbol, .first_tag = first_branch_selector_ptr->info.tag,
                .second_sym = second_branch_selector_ptr->info.symbol, .second_tag = second_branch_selector_ptr->info.tag
            };
            u64 product_symbol_handle = context.builder.symbol_handles.at(product_symbol);

            u64 first_branch_selector_handle  = context.first_ap.symbol_handles.at(first_branch_selector_ptr->info);
            u64 second_branch_selector_handle = context.second_ap.symbol_handles.at(second_branch_selector_ptr->info);

            auto& first_ap_successors  = context.first_ap.transition_fn[current_state->in_first][first_branch_selector_handle];
            auto& second_ap_successors = context.second_ap.transition_fn[current_state->in_second][second_branch_selector_handle];

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

    auto result = context.builder.build(context.worklist_info.handles.size());

    do_on_debug({
        result.debug_data = new Affine_Program<Branch_Product_Sym>::Debug_Data;

        auto& state_names = result.debug_data->state_names;
        for (auto& [product_state, handle] : context.worklist_info.handles) {
            std::string& first_ap_state_name  = first_ap.debug_data->state_names[product_state.in_first];
            std::string& second_ap_state_name = second_ap.debug_data->state_names[product_state.in_second];
            std::string& frontier_state_name  = frontier.debug_data->state_names[product_state.in_frontier];
            state_names[product_state.handle] = "(" + first_ap_state_name + ", " + second_ap_state_name + ", " + frontier_state_name + ")";
        }
    });

    return result;
}


std::ostream& operator<<(std::ostream& stream, const Branch_Product_Sym& sym) {
    const char* first_sym_str  = sym.first_tag == Subtree_Tag::LEFT ? "L" : "R";
    const char* second_sym_str = sym.second_tag == Subtree_Tag::LEFT ? "L" : "R";

    stream << "<c=" << sym.color << "-" << sym.first_sym << first_sym_str << "-" << sym.second_sym << second_sym_str << ">";

    return stream;
}

struct Color_Symbol_Build_Context {
    const SWTA& swta;
    Worklist_Construction_Context<Macrostate> worklist_state;
    std::map<Color_Symbol, u64> symbol_handles;
    std::map<State, std::vector<std::vector<State>> > result_transitions;
    Bit_Set final_states;
    u64 number_of_symbols;

    Color_Symbol_Build_Context(const SWTA& swta_ref) :
        swta(swta_ref),
        final_states(0),
        number_of_symbols(swta.number_of_internal_symbols() * swta.number_of_colors()) {}
};
void compute_post_in_color_sym_abstraction(Color_Symbol_Build_Context& context, const Macrostate& source, Color_Symbol& color_symbol) {
    auto [insert_pos, was_inserted] = context.symbol_handles.emplace(color_symbol, context.symbol_handles.size());
    u64 symbol_handle = insert_pos->second;

    for (State state : source.state_names) {
        auto& transition = context.swta.transitions[state][color_symbol.symbol][color_symbol.color];
        if (!transition.is_present()) return;  // All of the members need to define a transition
    }

    Macrostate imm_post (context.swta.number_of_states());

    for (State state : source.state_names) {
        auto& transition = context.swta.transitions[state][color_symbol.symbol][color_symbol.color];

        for (auto& component : transition.left.components) {
            imm_post.state_set.set_bit(component.state);
        }

        for (auto& component : transition.right.components) {
            imm_post.state_set.set_bit(component.state);
        }
    }

    auto post = context.worklist_state.mark_discovery(imm_post);
    {
        auto mut_post_ptr = const_cast<Macrostate*>(post);
        if (mut_post_ptr->state_names.empty()) {
            mut_post_ptr->populate_state_names_from_set();
        }
    }

    auto& transitions_from_state = context.result_transitions[source.handle];
    if (transitions_from_state.empty()) {
        transitions_from_state.resize(context.number_of_symbols);
    }

    transitions_from_state [symbol_handle].push_back(post->handle);
}

Color_Symbol_Abstraction build_color_internal_symbol_abstraction(const SWTA& swta) {
    Color_Symbol_Build_Context context(swta);

    {
        Macrostate initial_macrostate (context.swta.number_of_states(), swta.initial_states);
        context.worklist_state.mark_discovery(initial_macrostate);
    }

    while (context.worklist_state.has_more_to_explore()) {
        auto current_macrostate = context.worklist_state.extract();

        if (context.swta.states_with_leaf_transitions.is_superset(current_macrostate->state_set)) {
            context.final_states.grow_and_set_bit(current_macrostate->handle);
        }

        for (Internal_Symbol internal_sym = 0; internal_sym < context.swta.number_of_internal_symbols(); internal_sym++) {
            for (Color color = 0; color < context.swta.number_of_colors(); color++) {
                Color_Symbol color_symbol = {.color = color, .symbol = internal_sym};
                compute_post_in_color_sym_abstraction(context, *current_macrostate, color_symbol);
            }
        }
    }

    context.final_states.grow(context.worklist_state.handles.size());

    NFA::Transition_Fn transitions;
    transitions.resize(context.worklist_state.handles.size());

    for (State state = 0; state < context.worklist_state.handles.size(); state++) {
        auto& created_transitions = context.result_transitions[state];
        if (!created_transitions.empty()) {
            transitions[state] = std::move(created_transitions);
        } else {
            transitions[state].resize(context.number_of_symbols);
        }
    }

    NFA resulting_abstraction({0}, context.final_states, transitions);

    do_on_debug({
        resulting_abstraction.debug_data = new NFA::Debug_Data;

        for (auto& [macrostate, handle] : context.worklist_state.handles) {
            resulting_abstraction.debug_data->state_names[handle] = macrostate.to_string();
        }
    });

    Color_Symbol_Abstraction result = {
        .abstraction = resulting_abstraction,
        .symbol_handles = context.symbol_handles,
    };

    return result;
}
