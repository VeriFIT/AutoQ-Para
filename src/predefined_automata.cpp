#include "predefined_automata.hpp"
#include "arith.hpp"
#include "swta.hpp"

#include <stdexcept>
#include <string>


Def_Linear_Form Def_Coef::operator*(const Def_State& other) {
    return Def_Linear_Form(*this, other);
}


Linear_Form collect_tagged_subtrees_into_form(const std::vector<Def_Linear_Form>& components, Subtree_Tag tag) {
    Linear_Form form;
    for (auto& component : components) {
        if (component.tag == tag) {
            form.components.push_back({component.coef, component.state});
        }
    }
    return form;
}

WTT::Transition synthetize_wtt_transition(const std::vector<Def_Linear_Form>& left_subtree, const std::vector<Def_Linear_Form>& right_subtree) {
    Linear_Form ll = collect_tagged_subtrees_into_form(left_subtree, Subtree_Tag::LEFT);
    Linear_Form lr = collect_tagged_subtrees_into_form(left_subtree, Subtree_Tag::RIGHT);
    Linear_Form rl = collect_tagged_subtrees_into_form(right_subtree, Subtree_Tag::LEFT);
    Linear_Form rr = collect_tagged_subtrees_into_form(right_subtree, Subtree_Tag::RIGHT);

    return WTT::Transition(ll, lr, rl, rr);
}

SWTA::Transition synthetize_swta_transition(const std::vector<Def_Linear_Form>& left_subtree, const std::vector<Def_Linear_Form>& right_subtree) {
    Linear_Form ll = collect_tagged_subtrees_into_form(left_subtree,  Subtree_Tag::NONE);
    Linear_Form rr = collect_tagged_subtrees_into_form(right_subtree, Subtree_Tag::NONE);

    return SWTA::Transition(ll, rr);
}

WTT get_predefined_wtt(Predefined_WTT_Names name) {
    if (name == Predefined_WTT_Names::HADAMARD) {
        // q0 -> LEFT{ 1/sqrt(2)(q0, L) + 1/sqrt(2)(q0, R) }, RIGHT{ 1/sqrt(2)(q0, L) - 1/sqrt(2)(q0, R) }
        // q0(left) -> (left)

        Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE_OVER_SQRT2(), 0);
        Linear_Form ll ({ll_component});

        Linear_Form::Component lr_component (Algebraic_Complex_Number::ONE_OVER_SQRT2(), 0);
        Linear_Form lr ({lr_component});

        Linear_Form::Component rl_component (Algebraic_Complex_Number::ONE_OVER_SQRT2(), 0);
        Linear_Form rl ({rl_component});

        Linear_Form::Component rr_component (-Algebraic_Complex_Number::ONE_OVER_SQRT2(), 0);
        Linear_Form rr ({rr_component});

        WTT::Transition transition (ll, lr, rl, rr);

        WTT transducer ({{transition}}, {0}, {0});
        return transducer;
    }

    if (name == Predefined_WTT_Names::PARITY_CNOT) {
        // Computes parity as in BV circuit with secred (10)*
        // Qubits/states are labeled with A, B
        State a0 = 0; // Odd qubit, even parity
        State b0 = 1; // Even qubit, even parity
        State a1 = 2; // Odd qubit, odd parity
        State b1 = 3; // Even qubit, odd parity

        std::vector<WTT::Transition> work_qubit_transitions;  // @Care: The order of pushing is important, indexes are sourcestates
        std::vector<WTT::Transition> ancilla_transitions;

        // MOVE: a0 -w-> b0(L), b1(r)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), b0);
            Linear_Form ll ({ll_component});

            Linear_Form::Component rr_component (Algebraic_Complex_Number::ONE(), b1);
            Linear_Form rr ({rr_component});

            WTT::Transition transition (ll, {}, {}, rr);
            work_qubit_transitions.push_back(transition);
        }

        // MOVE: b0 -w-> a0(L), a0(r)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), a0);
            Linear_Form ll ({ll_component});

            WTT::Transition transition (ll, {}, {}, ll);
            work_qubit_transitions.push_back(transition);
        }

        // MOVE: a1 -w-> b1(L), b0(r)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), b1);
            Linear_Form ll ({ll_component});

            Linear_Form::Component rr_component (Algebraic_Complex_Number::ONE(), b0);
            Linear_Form rr ({rr_component});

            WTT::Transition transition (ll, {}, {}, rr);
            work_qubit_transitions.push_back(transition);
        }

        // MOVE: b1 -w-> a1(L), a1(r)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), a1);
            Linear_Form ll ({ll_component});

            WTT::Transition transition (ll, {}, {}, ll);
            work_qubit_transitions.push_back(transition);
        }

        // MOVE: a0 -a-> b0(L), b0(r), b0 -a-> b0(L), b0(r)
        {
            std::vector<Def_Linear_Form> left_subtree  {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::LEFT};
            std::vector<Def_Linear_Form> right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::RIGHT};
            WTT::Transition transition = synthetize_wtt_transition(left_subtree, right_subtree);
            ancilla_transitions.push_back(transition);
            ancilla_transitions.push_back(transition);
        }

        // MOVE: a1 -a-> b0(R), b0(L), b1 -a-> b0(R), b0(L)
        {
            std::vector<Def_Linear_Form> left_subtree  {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::RIGHT};
            std::vector<Def_Linear_Form> right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::LEFT};
            WTT::Transition transition = synthetize_wtt_transition(left_subtree, right_subtree);
            ancilla_transitions.push_back(transition);
            ancilla_transitions.push_back(transition);
        }

        WTT transducer ({work_qubit_transitions, ancilla_transitions}, {b0}, {a0});
        return transducer;
    }

    throw std::runtime_error("Unknown WTT. " + std::to_string(static_cast<u64>(name)));
}


SWTA get_predefined_swta(Predefined_SWTA_Names name) {
    using DLF = std::vector<Def_Linear_Form>;
    using ACN = Algebraic_Complex_Number;

    if (name == Predefined_SWTA_Names::BV_EXAMPLE_10STAR_PRE) {
        SWTA::Metadata metadata {.number_of_internal_symbols = 2, .number_of_colors = 1};
        SWTA::Transition_Builder builder (metadata);

        State q0    = 0;
        State q1    = 1;
        State q_bot = 2;  // Accepts zero trees of any height

        Internal_Symbol sym_w = 0;
        Internal_Symbol sym_a = 1;

        Color c = 0;

        { // Transition q0 ----> w(q0, q_bot)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ONE()) * q0};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q0, sym_w, c, transition);
        }

        { // Transition q0 ----> a(q_bot, q1)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * q1};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q0, sym_a, c, transition);
        }

        { // Transition q_bot ----> a(q_bot, q_bot)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_bot, sym_w, c, transition);
             builder.add_transition(q_bot, sym_a, c, transition);
        }

        Bit_Set leaf_states (3, {q1, q_bot});
        std::vector<State> initial_states {q0};
        auto transition_fn = builder.build(3);

        SWTA result (transition_fn, initial_states, leaf_states);
    }

    if (name == Predefined_SWTA_Names::BV_EXAMPLE_10STAR_POST) {
        State q_g   = 0;
        State q_h   = 1;
        State q_c   = 2;
        State q_bot = 3;

        Internal_Symbol sym_w = 0;
        Internal_Symbol sym_a = 1;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        Color c = 0;

        { // g ----> w(0 q_bot, h)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * q_h};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_g, sym_w, c, transition);
        }

        { // g ----> a(0 q_bot, c)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * q_c};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_g, sym_a, c, transition);
        }

        { // h ----> w(g, 0 q_bot)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ONE()) * q_g};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_h, sym_w, c, transition);
        }

        { // h ----> a(0 q_bot, c)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * q_c};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_g, sym_a, c, transition);
        }

        builder.add_bot_state_transitions(q_bot);

        std::vector<State> initial_states ({q_g});
        Bit_Set leaf_states (4, {q_bot, q_c});
        auto transition_fn = builder.build(4);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::BV_EXAMPLE_10STAR_RESULT) {
        State q_gamma   = 0;
        State q_delta   = 1;
        State q_epsilon = 2;
        State q_mu      = 3;
        State q_sigma   = 4;
        State q_bot     = 5;

        u64 state_cnt = 6;

        Internal_Symbol sym_w = 0;
        Internal_Symbol sym_a = 1;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        Color c = 0;

        ACN one_half (1, 0, 0, 0, 2);

        { // gamma ----> w(1/2 delta + 1/2 epsilon, 1/2 delta - 1/2 epsilon)
             DLF left_subtree  {Def_Coef(one_half) * q_delta, Def_Coef(one_half) * q_epsilon};
             DLF right_subtree {Def_Coef(one_half) * q_delta, Def_Coef(-one_half) * q_epsilon};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_gamma, sym_w, c, transition);
        }

        { // gamma ----> a(0 q_bot, q_sigma)
             DLF left_subtree  {Def_Coef(ACN::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(ACN::ONE()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_gamma, sym_a, c, transition);
        }

        { // delta ----> w(q_gamma, 0 q_bot)
             DLF left_subtree  {Def_Coef(ACN::ONE()) * q_gamma};
             DLF right_subtree {Def_Coef(ACN::ZERO()) * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_delta, sym_w, c, transition);
        }

        { // delta ----> a(0 q_bot, q_sigma)
             DLF left_subtree  {Def_Coef(ACN::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(ACN::ONE()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_delta, sym_a, c, transition);
        }

        { // epsilon ----> w(q_mu, 0 q_bot)
             DLF left_subtree  {Def_Coef(ACN::ONE()) * q_mu};
             DLF right_subtree {Def_Coef(ACN::ZERO()) * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_epsilon, sym_w, c, transition);
        }

        { // epsilon ----> a(0 q_bot, - q_sigma)
             DLF left_subtree  {Def_Coef(ACN::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(-ACN::ONE()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_epsilon, sym_a, c, transition);
        }

        { // mu ----> w(1/2 epsilon + 1/2 delta, 1/2 delta - 1/2 epsilon)
             DLF left_subtree  {Def_Coef(one_half) * q_epsilon, Def_Coef(one_half) * q_delta};
             DLF right_subtree {Def_Coef(one_half) * q_epsilon, Def_Coef(-one_half) * q_delta};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_mu, sym_w, c, transition);
        }

        { // mu ----> a(0 q_bot, q_sigma)
             DLF left_subtree  {Def_Coef(ACN::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(-ACN::ONE()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_mu, sym_a, c, transition);
        }

        builder.add_bot_state_transitions(q_bot);

        std::vector<State> initial_states ({q_gamma});
        Bit_Set leaf_states (state_cnt, {q_bot, q_sigma});
        auto transition_fn = builder.build(state_cnt);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::TRIVIAL_BOT) {
        State q_bot   = 0;

        Internal_Symbol sym_w = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        Color c = 0;

        builder.add_bot_state_transitions(q_bot);

        std::vector<State> initial_states ({q_bot});
        Bit_Set leaf_states (1, {q_bot});
        auto transition_fn = builder.build(1);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::TRIVIAL_ONES) {
        State q_one   = 0;

        Internal_Symbol sym_w = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        Color c = 0;

        { // q_one ----> w(q_one, q_one)
             DLF left_subtree  {Def_Coef(ACN::ONE()) * q_one};
             DLF right_subtree {Def_Coef(ACN::ONE()) * q_one};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_one, sym_w, c, transition);
        }

        std::vector<State> initial_states ({q_one});
        Bit_Set leaf_states (1, {q_one});
        auto transition_fn = builder.build(1);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }


    throw std::runtime_error("No definition for the predefined SWTA: " + std::to_string(static_cast<u64>(name)));
}
