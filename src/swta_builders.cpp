#include "swta_builders.hpp"
#include "basics.hpp"
#include "nfa.hpp"
#include "swta.hpp"


struct Staircase_State {
    u16 first_nth_qubit  = 0;
    u16 second_nth_qubit = 0;

    u64 handle          =  0;
    s64 first           = -1;
    s64 second          = -1;

    bool operator<(const Staircase_State& other) const {
        INSERT_LEX_LT_CODE(this->first, other.first);
        return this->second < other.second;
    }

    bool operator==(const Staircase_State& other) const {
        return (first == other.first) && (second == other.second);
    }
};

std::ostream& operator<<(std::ostream& stream, const Staircase_State& state) {
    stream << "(";

    if (state.first >= 0) {
        stream << state.first;
    } else {
        stream << "ID";
    }

    stream << ", ";

    if (state.second >= 0) {
        stream << state.second;
    } else {
        stream << "ID";
    }

    stream << ")";

    return stream;
}

struct Staircase_Construction_Ctx {
    Worklist_Construction_Context<Staircase_State> worklist_state;
    WTT& box;
    const std::vector<Internal_Symbol>& box_inputs;
    u64 box_offset;
    u64 terminating_symbol;
    std::vector<const Staircase_State*> states_by_handle;

    Staircase_Construction_Ctx(WTT& box, const std::vector<Internal_Symbol>& inputs, u64 box_offset, u64 term_symbol) : box(box), box_inputs(inputs), box_offset(box_offset), terminating_symbol(term_symbol) {}

    const Staircase_State* mark_discovery(Staircase_State& state) {
        u64 size_before_insert = this->worklist_state.handles.size();
        auto result = this->worklist_state.mark_discovery(state);
        u64 size_after_insert = this->worklist_state.handles.size();

        if (size_before_insert != size_after_insert) {
            this->states_by_handle.push_back(result);
        }

        return result;
    }

};


u64 count_different_symbols_in_box_inputs(const std::vector<Internal_Symbol>& inputs, u64 terminating_symbol) {
    Bit_Set bit_set(inputs.size() + 1);

    bit_set.grow_and_set_bit(terminating_symbol);

    for (Internal_Symbol internal_symbol : inputs) {
        bit_set.grow_and_set_bit(internal_symbol);
    }

    return bit_set.popcount();
}

void extend_form_with_product_and_note_discoveries(const Staircase_State& source_state, Linear_Form& destination, const Linear_Form& first, const Linear_Form& second, Staircase_Construction_Ctx& ctx) {
    for (auto& outer_comp : second.components) {
        for (auto& inner_comp : first.components) {

            Staircase_State state = {
                .first_nth_qubit  = source_state.first_nth_qubit + 1, // The product state will read the (n+1)-th qubit
                .second_nth_qubit = source_state.second_nth_qubit + 1,
                .first  = static_cast<s64>(inner_comp.state),
                .second = static_cast<s64>(outer_comp.state),
            };

            if (state.first_nth_qubit > ctx.box_inputs.size()) { // Since qubit number if 0-based, if the equiality hold, we have already exhausted all qubits read by the box
                assert (ctx.box.states_with_leaf_transitions.get_bit_value(state.first));
                state.first  = state.second;
                state.second = -1;
            }

            auto imm_state = ctx.mark_discovery(state);

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

            destination.components.push_back({coef, imm_state->handle});
        }
    }
}

WTT::Transition calculate_product_of_two_forms(Staircase_Construction_Ctx& ctx, const Staircase_State& source_state, const WTT::Transition& first_transition, const WTT::Transition& second_transition) {
    Linear_Form ll, lr, rl, rr;
    { // ll
        extend_form_with_product_and_note_discoveries(source_state, ll, first_transition.ll, second_transition.ll, ctx);
        extend_form_with_product_and_note_discoveries(source_state, ll, first_transition.rl, second_transition.lr, ctx);
    }

    { // lr
        extend_form_with_product_and_note_discoveries(source_state, lr, first_transition.lr, second_transition.ll, ctx);
        extend_form_with_product_and_note_discoveries(source_state, lr, first_transition.rr, second_transition.lr, ctx);
    }

    { // rl
        extend_form_with_product_and_note_discoveries(source_state, rl, first_transition.ll, second_transition.rl, ctx);
        extend_form_with_product_and_note_discoveries(source_state, rl, first_transition.rl, second_transition.rr, ctx);
    }

    { // rr
        extend_form_with_product_and_note_discoveries(source_state, rr, first_transition.lr, second_transition.rl, ctx);
        extend_form_with_product_and_note_discoveries(source_state, rr, first_transition.rr, second_transition.rr, ctx);
    }
    WTT::Transition result (ll, lr, rl, rr);
    return result;
}

template <typename Padder>
void pad_linear_form_with_ids(Staircase_Construction_Ctx& ctx, Padder padder, Linear_Form& destination, const Linear_Form& source) {
    for (auto& component : source.components) {
        Staircase_State padded_state = padder(component.state);

        auto discovered_state = ctx.mark_discovery(padded_state);

        // We don't have to check whether there already exists such a padded state within the destination
        // since the source form should contain only unique states, so by padding with identity we get
        // again unique states
        Linear_Form::Component new_component(component.coef, discovered_state->handle);
        destination.components.push_back(new_component);
    }
}

template <typename Padder>
WTT::Transition pad_transition_using_padder(Staircase_Construction_Ctx& ctx, WTT::Transition& transition_to_pad, const Staircase_State* source_state, Padder padder) {
    Linear_Form ll, lr, rl, rr;
    pad_linear_form_with_ids(ctx, padder, ll, transition_to_pad.ll);
    pad_linear_form_with_ids(ctx, padder, lr, transition_to_pad.lr);
    pad_linear_form_with_ids(ctx, padder, rl, transition_to_pad.rl);
    pad_linear_form_with_ids(ctx, padder, rr, transition_to_pad.rr);
    return WTT::Transition(ll, lr, rl, rr);
}

WTT::Transition pad_transition_with_ids_right(Staircase_Construction_Ctx& ctx, WTT::Transition& transition_to_pad, const Staircase_State* source_state) {
    auto right_padder = [&ctx, &source_state](State present_state) {
        u64 nth_qubit_result_will_read = source_state->first_nth_qubit + 1;

        s64 pad_value = -1;
        if (nth_qubit_result_will_read == ctx.box_offset) { // We can start the next box, so instead of padding with -1, we put there the initial state of the box
            pad_value = static_cast<s64>(ctx.box.initial_states[0]);
        }

        Staircase_State padded_state = {
            .first_nth_qubit = nth_qubit_result_will_read,
            .second_nth_qubit = 0,
            .first = static_cast<s64>(present_state),
            .second = pad_value,
        };

        if (padded_state.second >= 0) {
            std::cout << "From state: " << *source_state << " on "
                      << nth_qubit_result_will_read << "-th qubit we active the second component, producing a state: "
                      << padded_state << "\n";
        }

        return padded_state;
    };
    auto result = pad_transition_using_padder(ctx, transition_to_pad, source_state, right_padder);
    return result;
}

void unrestart_linear_form(Staircase_Construction_Ctx& ctx, Linear_Form& dest, const Linear_Form& source) {
    for (auto& component : source.components) {
        auto staircase_state = ctx.states_by_handle[component.state];
        Staircase_State unrestared_state = {
            .first_nth_qubit = staircase_state->first_nth_qubit,
            .second_nth_qubit = staircase_state->second_nth_qubit,
            .first = staircase_state->first,
            .second = -1,
        };
        auto imm_state = ctx.mark_discovery(unrestared_state);

        Linear_Form::Component new_component (component.coef, imm_state->handle);
        dest.components.push_back(new_component);
    }
}

WTT::Transition craft_unrestarted_transition(Staircase_Construction_Ctx& ctx, WTT::Transition& transition_to_unrestart) {
    Linear_Form ll, lr, rl, rr;

    unrestart_linear_form(ctx, ll, transition_to_unrestart.ll);
    unrestart_linear_form(ctx, lr, transition_to_unrestart.lr);
    unrestart_linear_form(ctx, rl, transition_to_unrestart.rl);
    unrestart_linear_form(ctx, rr, transition_to_unrestart.rr);

    WTT::Transition result (ll, lr, rl, rr);
    std::cout << "UNRestarting transition " << transition_to_unrestart << " yielded " << result << "\n";
    return result;
}


WTT perform_staircase_construction(WTT& box, const std::vector<Internal_Symbol>& box_inputs, u64 box_offset, u64 terminating_symbol) {
    // Check whether the overlap is sufficiently small for products to be 2 tuples
    assert (box_inputs.size() <= 2*box_offset);

    u64 number_of_internal_symbols = count_different_symbols_in_box_inputs(box_inputs, terminating_symbol);
    SWTA::Metadata metadata = { .number_of_internal_symbols = number_of_internal_symbols, .number_of_colors = 0 };
    WTT_Builder builder (metadata);

    Staircase_Construction_Ctx ctx(box, box_inputs, box_offset, terminating_symbol);

    Staircase_State initial_state = { .first_nth_qubit = 0, .first = static_cast<s64>(box.initial_states[0]), .second = -1, };
    auto imm_initial_state = ctx.mark_discovery(initial_state);
    builder.mark_state_initial(imm_initial_state->handle);

    while (ctx.worklist_state.has_more_to_explore()) {
        auto imm_state = ctx.worklist_state.extract();

        std::cout << "Exploring " << *imm_state << " :: "<< imm_state->handle << "\n";

        if (imm_state->second < 0) {

            if (imm_state->first < 0 || ctx.box.states_with_leaf_transitions.get_bit_value(imm_state->first)) {
                builder.mark_state_final(imm_state->handle);
                continue;;
            }

            Internal_Symbol transition_symbol = box_inputs[imm_state->first_nth_qubit];
            auto& first_transition = box.transitions[imm_state->first][transition_symbol];
            WTT::Transition padded_transition = pad_transition_with_ids_right(ctx, first_transition, imm_state);

            bool did_restart_happen =  (imm_state->first_nth_qubit + 1 == ctx.box_offset);
            if (did_restart_happen) {
                auto unrestarted_transition = craft_unrestarted_transition(ctx, padded_transition);
                builder.add_transition(imm_state->handle, ctx.terminating_symbol, unrestarted_transition);
            }

            builder.add_transition(imm_state->handle, transition_symbol, padded_transition);
            continue;
        }

        // We have both of the states being active
        Internal_Symbol transition_symbol = box_inputs[imm_state->second_nth_qubit]; // The first might be in leaf state, which would lead to out-of-bounds access

        auto& first_transition  = box.transitions[imm_state->first][transition_symbol];
        auto& second_transition = box.transitions[imm_state->second][transition_symbol];

        if (!first_transition.is_present() || !second_transition.is_present()) continue;

        WTT::Transition transition = calculate_product_of_two_forms(ctx, *imm_state, first_transition, second_transition);
        std::cout << transition << "\n";
        builder.add_transition(imm_state->handle, transition_symbol, transition);
    }

    WTT result = builder.build(ctx.states_by_handle.size());

    do_on_debug({
       for (auto& [state, handle] : ctx.worklist_state.handles) {
           std::cout << state << " = " << handle << "\n";
       }
    });

    return result;
}
