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

WTT get_predefined_wtt(Predefined_WTT_Names name, const SWTA::Metadata& swta_metadata) {
    using DLF = std::vector<Def_Linear_Form>;
    using ACN = Algebraic_Complex_Number;

    if (name == Predefined_WTT_Names::HADAMARD) {
        // q0 -> LEFT{ 1/sqrt(2)(q0, L) + 1/sqrt(2)(q0, R) }, RIGHT{ 1/sqrt(2)(q0, L) - 1/sqrt(2)(q0, R) }
        // q0(left) -> (left)

        State q0 = 0;

        std::vector<std::vector<WTT::Transition>> transitions_by_symbol;
        transitions_by_symbol.resize(1);

        Def_Coef one_over_sqrt2       ( Algebraic_Complex_Number::ONE_OVER_SQRT2());
        Def_Coef minus_one_over_sqrt2 (-Algebraic_Complex_Number::ONE_OVER_SQRT2());

        for (Internal_Symbol symbol = 0; symbol < swta_metadata.number_of_internal_symbols; symbol++) {
            DLF left_subtree  {one_over_sqrt2 * q0 * Subtree_Tag::LEFT, one_over_sqrt2 * q0 * Subtree_Tag::RIGHT};
            DLF right_subtree {one_over_sqrt2 * q0 * Subtree_Tag::LEFT, minus_one_over_sqrt2 * q0 * Subtree_Tag::RIGHT};
            WTT::Transition transition = synthetize_wtt_transition(left_subtree, right_subtree);

            transitions_by_symbol[0].push_back(transition);
        }

        WTT transducer (transitions_by_symbol, {0}, {0});
        return transducer;
    }

    if (name == Predefined_WTT_Names::PARITY_CNOT) {
        // Computes parity as in BV circuit with secred (10)*
        // Qubits/states are labeled with A, B
        State a0 = 0; // Odd qubit, even parity
        State b0 = 1; // Even qubit, even parity
        State a1 = 2; // Odd qubit, odd parity
        State b1 = 3; // Even qubit, odd parity

        size_t state_cnt = 4;
        size_t internal_symbol_cnt = 2;

        Internal_Symbol work_qubit = 0, ancilla = 1;

        std::vector<std::vector<WTT::Transition>> transitions;
        transitions.resize(4);
        for (State state = 0; state < state_cnt; state++) {
            transitions[state].resize(internal_symbol_cnt);
        }

        // MOVE: a0 -w-> b0(L), b1(r)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), b0);
            Linear_Form ll ({ll_component});

            Linear_Form::Component rr_component (Algebraic_Complex_Number::ONE(), b1);
            Linear_Form rr ({rr_component});

            WTT::Transition transition (ll, {}, {}, rr);

            transitions[a0][work_qubit] = transition;
        }

        // MOVE: b0 -w-> a0(L), a0(r)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), a0);
            Linear_Form ll ({ll_component});

            WTT::Transition transition (ll, {}, {}, ll);

            transitions[b0][work_qubit] = transition;
        }

        // MOVE: a1 -w-> b1(L), b0(r)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), b1);
            Linear_Form ll ({ll_component});

            Linear_Form::Component rr_component (Algebraic_Complex_Number::ONE(), b0);
            Linear_Form rr ({rr_component});

            WTT::Transition transition (ll, {}, {}, rr);

            transitions[a1][work_qubit] = transition;
        }

        // MOVE: b1 -w-> a1(L), a1(r)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), a1);
            Linear_Form ll ({ll_component});

            WTT::Transition transition (ll, {}, {}, ll);

            transitions[b1][work_qubit] = transition;
        }

        // MOVE: a0 -a-> b0(L), b0(r)
        // MOVE: b0 -a-> b0(L), b0(r)
        {
            std::vector<Def_Linear_Form> left_subtree  {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::LEFT};
            std::vector<Def_Linear_Form> right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::RIGHT};
            WTT::Transition transition = synthetize_wtt_transition(left_subtree, right_subtree);

            transitions[a0][ancilla] = transition;
            transitions[b0][ancilla] = transition;
        }

        // MOVE: a1 -a-> b0(R), b0(L)
        // MOVE: b1 -a-> b0(R), b0(L)
        {
            std::vector<Def_Linear_Form> left_subtree  {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::RIGHT};
            std::vector<Def_Linear_Form> right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::LEFT};
            WTT::Transition transition = synthetize_wtt_transition(left_subtree, right_subtree);

            transitions[a1][ancilla] = transition;
            transitions[b1][ancilla] = transition;
        }

        WTT transducer (transitions, {b0}, {a0});
        return transducer;
    }

    if (name == Predefined_WTT_Names::GROVER_FIRST_MULTI_Z) {
        State q_init         = 0;
        State q_w            = 1;
        State q_a            = 2;
        State q_last_ancilla = 3;
        State q_id           = 4;
        State q_leaf         = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);
        builder.mark_state_final(q_id);

        Def_Coef one (ACN::ONE());
        { // q_init --w--> ( Id(L), q_w(R) )
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_w  * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_w --w--> ( Id(L), q_a(R) )
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_a  * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_w, working_qubit, transition);
        }
        { // q_w --W'--> ( Id(L), q_last_ancilla(R) )
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_last_ancilla  * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_w, last_working_qubit, transition);
        }

        { // q_a --a--> ( q_w(L), q_w(R) )
             DLF left_subtree  {one * q_w * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_w * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_a, ancilla, transition);
        }

        // q_id --w-->  ( q_id(L), q_id(R) )
        // q_id --W'--> ( q_id(L), q_id(R) )
        // q_id --a-->  ( q_id(L), q_id(R) )
        {
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_id * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);

             builder.add_transition(q_id, working_qubit, transition);
             builder.add_transition(q_id, last_working_qubit, transition);
             builder.add_transition(q_id, ancilla, transition);
        }

        // q_last_ancilla --a-->  (q_leaf(L), -1 q_leaf(R) )
        {
             DLF left_subtree  {one * q_leaf * Subtree_Tag::LEFT};
             DLF right_subtree {Def_Coef(-ACN::ONE()) * q_leaf * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);

             builder.add_transition(q_last_ancilla, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::GROVER_SECOND_MULTI_Z) {
        State q_init         = 0;
        State q_w            = 1;
        State q_a            = 2;
        State q_last_ancilla = 3;
        State q_id           = 4;
        State q_leaf         = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);
        builder.mark_state_final(q_id);

        Def_Coef one (ACN::ONE());
        { // q_init --w--> ( Id(L), q_w(R) )
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_w  * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_w --w--> ( Id(L), q_a(R) )
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_a  * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_w, working_qubit, transition);
        }
        { // q_w --W'--> ( 1 * q_last_ancilla(L), -1 q_last_ancilla(R) )
             DLF left_subtree  {                  one * q_last_ancilla * Subtree_Tag::LEFT};
             DLF right_subtree {Def_Coef(-ACN::ONE()) * q_last_ancilla * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_w, last_working_qubit, transition);
        }

        { // q_a --a--> ( q_w(L), q_w(R) )
             DLF left_subtree  {one * q_w * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_w * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_a, ancilla, transition);
        }

        // q_id --w-->  ( q_id(L), q_id(R) )
        // q_id --W'--> ( q_id(L), q_id(R) )
        // q_id --a-->  ( q_id(L), q_id(R) )
        {
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_id * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);

             builder.add_transition(q_id, working_qubit, transition);
             builder.add_transition(q_id, last_working_qubit, transition);
             builder.add_transition(q_id, ancilla, transition);
        }

        // q_last_ancilla --a-->  (q_leaf(L), q_leaf(R) )
        {
             DLF left_subtree  {one * q_leaf * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_leaf * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);

             builder.add_transition(q_last_ancilla, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::GROVER_X) {
        State q_init         = 0;
        State q_x            = 1;
        State q_a            = 2;
        State q_last_ancilla = 3;
        State q_leaf         = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);

        Def_Coef one (ACN::ONE());
        { // q_init --w--> ( q_x(R), q_x(L) )
             DLF left_subtree  {one * q_x * Subtree_Tag::RIGHT };
             DLF right_subtree {one * q_x * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_x --w--> ( q_a(R), q_a(L) )
             DLF left_subtree  {one * q_a * Subtree_Tag::RIGHT };
             DLF right_subtree {one * q_a * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_x, working_qubit, transition);
        }

        { // q_x --W'--> ( q_last_ancilla(R), q_last_ancilla(L) )
             DLF left_subtree  {one * q_last_ancilla * Subtree_Tag::RIGHT };
             DLF right_subtree {one * q_last_ancilla * Subtree_Tag::LEFT  };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_x, last_working_qubit, transition);
        }

        { // q_a --a--> ( q_x(L), q_x(R) )
             DLF left_subtree  {one * q_x * Subtree_Tag::LEFT  };
             DLF right_subtree {one * q_x * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_a, ancilla, transition);
        }

        { // q_last_ancilla --a--> ( q_leaf(L), q_leafw(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_last_ancilla, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::GROVER_H) {
        State q_init         = 0;
        State q_h            = 1;
        State q_a            = 2;
        State q_last_ancilla = 3;
        State q_leaf         = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);

        Def_Coef omega       ( ACN::ONE_OVER_SQRT2());
        Def_Coef minus_omega (-ACN::ONE_OVER_SQRT2());
        Def_Coef one         (ACN::ONE());

        { // q_init --w--> ( omega*q_h(L) + omega*q_h(R), omega*q_h(L) - omega*q_h(R) )
             DLF left_subtree  { omega * q_h * Subtree_Tag::LEFT, omega * q_h * Subtree_Tag::RIGHT};
             DLF right_subtree { omega * q_h * Subtree_Tag::LEFT, minus_omega * q_h * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_h --w--> ( omega*q_a(L) + omega*q_a(R), omega*q_a(L) - omega*q_a(R) )
             DLF left_subtree  { omega * q_a * Subtree_Tag::LEFT, omega * q_a * Subtree_Tag::RIGHT};
             DLF right_subtree { omega * q_a * Subtree_Tag::LEFT, minus_omega * q_a * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_h, working_qubit, transition);
        }

        { // q_h --W'--> ( omega*q_last_ancilla(L) + omega*q_last_ancilla(R), omega*q_last_ancilla(L) - omega*q_last_ancilla(R) )
             DLF left_subtree  { omega * q_last_ancilla * Subtree_Tag::LEFT, omega * q_last_ancilla * Subtree_Tag::RIGHT};
             DLF right_subtree { omega * q_last_ancilla * Subtree_Tag::LEFT, minus_omega * q_last_ancilla * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_h, last_working_qubit, transition);
        }

        { // q_a --a--> ( q_h(L), q_h(R) )
             DLF left_subtree  { one * q_h * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_h * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_a, ancilla, transition);
        }

        { // q_last_ancilla --a--> ( q_leaf(L), q_leaf(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_last_ancilla, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::GROVER_FIRST_MULTI_Z_USING_CCX) {
        State q_init         = 0;
        State q_no_swap      = 1;
        State q_maybe_swap   = 2;
        State q_do_swap      = 3;
        State q_no_apply     = 4;
        State q_apply        = 5;
        State q_leaf         = 6;

        u64 number_of_states = 7;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);

        Def_Coef one       ( ACN::ONE());
        Def_Coef minus_one (-ACN::ONE());

        { // q_init --w--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  {one * q_no_swap    * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_maybe_swap * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_no_swap --w--> ( q_no_swap(L), q_no_swap(R) )
             DLF left_subtree  {one * q_no_swap * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_no_swap * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, working_qubit, transition);
        }

        { // q_no_swap --a--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  { one * q_no_swap    * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_maybe_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, ancilla, transition);
        }

        { // q_maybe_swap --w--> ( q_no_swap(L), q_swap(R) )
             DLF left_subtree  { one * q_no_swap * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_do_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap, working_qubit, transition);
        }

        { // q_swap --a--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  { one * q_no_swap    * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_maybe_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_do_swap, ancilla, transition);
        }

        { // q_no_swap --W'--> ( q_no_apply(L), q_no_apply(R) )
             DLF left_subtree  { one * q_no_apply * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_no_apply * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, last_working_qubit, transition);
        }

        { // q_maybe_swap --W'--> ( q_no_apply(L), q_apply(R) )
             DLF left_subtree  { one * q_no_apply * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_apply * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap, last_working_qubit, transition);
        }

        { // q_no_apply --a--> ( leaf(L), leaf(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_apply, ancilla, transition);
        }

        { // q_apply --a--> ( leaf(L), -leaf(R) )
             DLF left_subtree  {       one * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { minus_one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_apply, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::GROVER_SECOND_MULTI_Z_USING_CCX) {
        State q_init         = 0;
        State q_no_swap      = 1;
        State q_maybe_swap   = 2;
        State q_do_swap      = 3;
        State q_last_anc     = 4;
        State q_leaf         = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);

        Def_Coef one       ( ACN::ONE());
        Def_Coef minus_one (-ACN::ONE());

        { // q_init --w--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  {one * q_no_swap    * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_maybe_swap * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_no_swap --w--> ( q_no_swap(L), q_no_swap(R) )
             DLF left_subtree  {one * q_no_swap * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_no_swap * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, working_qubit, transition);
        }

        { // q_no_swap --a--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  { one * q_no_swap    * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_maybe_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, ancilla, transition);
        }

        { // q_maybe_swap --w--> ( q_no_swap(L), q_swap(R) )
             DLF left_subtree  { one * q_no_swap * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_do_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap, working_qubit, transition);
        }

        { // q_swap --a--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  { one * q_no_swap    * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_maybe_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_do_swap, ancilla, transition);
        }

        { // q_no_swap --W'--> ( q_last_anc(L), q_last_anc(R) )
             DLF left_subtree  { one * q_last_anc * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_last_anc * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, last_working_qubit, transition);
        }

        { // q_maybe_swap --W'--> ( q_last_anc(L), -q_last_anc(R) )
             DLF left_subtree  {       one * q_last_anc * Subtree_Tag::LEFT };
             DLF right_subtree { minus_one * q_last_anc * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap, last_working_qubit, transition);
        }

        { // q_last_anc --a--> ( leaf(L), leaf(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_last_anc, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::TEST_STAIRCASE_IDENTITY3) {
        State q3     = 0;
        State q2     = 1;
        State q1     = 2;
        State q_leaf = 3;

        u64 number_of_states = 4;

        Internal_Symbol working_qubit = 0, ancilla = 1;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q3);
        builder.mark_state_final(q_leaf);

        Def_Coef one (  ACN::ONE());
        Def_Coef two (ACN(2, 0, 0, 0, 0));

        { // q3 --w--> ( q2(L), q2(R) )
             DLF left_subtree  {two * q2 * Subtree_Tag::LEFT};
             DLF right_subtree {two * q2 * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q3, working_qubit, transition);
        }

        { // q2 --a--> ( q1(L), q1(R) )
             DLF left_subtree  {two * q1 * Subtree_Tag::LEFT};
             DLF right_subtree {two * q1 * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q2, ancilla, transition);
        }

        { // q1 --w--> ( q_leaf(L), q_leaf(R) )
             DLF left_subtree  {one * q_leaf * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_leaf * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q1, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_UMA1 || name == Predefined_WTT_Names::ADDER_MAJ3) {
        State q_maybe_swap0 = 0;
        State q_maybe_swap1 = 1;
        State q_no_swap1    = 2;
        State q_no_swap2    = 3;
        State q_swap        = 4;
        State q_leaf        = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_maybe_swap0);
        builder.mark_state_final(q_leaf);

        Def_Coef one (  ACN::ONE());

        { // q_maybe_swap0 --w--> ( q_no_swap1(L), q_maybe_swap1(R) )
             DLF left_subtree  {one * q_no_swap1    * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_maybe_swap1 * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap0, working_qubit, transition);
        }

        { // q_no_swap1 --w--> ( q_no_swap2(L), q_no_swap2(R) )
             DLF left_subtree  {one * q_no_swap2 * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_no_swap2 * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap1, working_qubit, transition);
        }

        { // q_no_swap2 --w--> ( q_leaf(L), q_leaf(R) )
             DLF left_subtree  {one * q_leaf * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_leaf * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap2, working_qubit, transition);
        }

        { // q_maybe_swap1 --w--> ( q_no_swap2(L), q_swap(R) )
             DLF left_subtree  {one * q_no_swap2 * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_swap     * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap1, working_qubit, transition);
        }

        { // q_swap --w--> ( q_leaf(R), q_leaf(L) )
             DLF left_subtree  {one * q_leaf * Subtree_Tag::RIGHT};
             DLF right_subtree {one * q_leaf * Subtree_Tag::LEFT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_swap, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_UMA2 || name == Predefined_WTT_Names::ADDER_MAJ2) {
        State root = 0;
        State p_star1  = 1; // Project away all branches *1
        State p_1      = 2;
        State e_star1  = 3; // Extract (keep values of) the branches *1
        State e_1      = 4;
        State q_leaf   = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(root);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // root --w--> ( p_star1(L) + e_star1(R), p_star1(R) + e_star1(L) )
             DLF left_subtree  { one * p_star1 * Subtree_Tag::LEFT, one * e_star1 * Subtree_Tag::RIGHT };
             DLF right_subtree { one * p_star1 * Subtree_Tag::RIGHT, one * e_star1 * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(root, working_qubit, transition);
        }

        { // p_star1 --w--> ( p_1(L), p_1(R))
             DLF left_subtree  { one * p_1 * Subtree_Tag::LEFT };
             DLF right_subtree { one * p_1 * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(p_star1, working_qubit, transition);
        }

        { // p_1 --w--> ( q_leaf(L), 0 q_leaf(R) )
             DLF left_subtree  { one  * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { zero * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(p_1, working_qubit, transition);
        }

        { // e_star1 --w--> ( e_1(L), e_1(R))
             DLF left_subtree  { one * e_1 * Subtree_Tag::LEFT };
             DLF right_subtree { one * e_1 * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(e_star1, working_qubit, transition);
        }

        { // e_1 --w--> ( 0 q_leaf(L), q_leaf(R))
             DLF left_subtree  { zero * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { one  * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(e_1, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_UMA3) {
        State root      = 0;
        State q_no_swap = 1;
        State q_swap    = 2;
        State q_id_1    = 3;
        State q_leaf    = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(root);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());

        { // root --w--> ( p_no_swap(L), q_swap(R) )
             DLF left_subtree  { one * q_no_swap * Subtree_Tag::LEFT};
             DLF right_subtree { one * q_swap    * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(root, working_qubit, transition);
        }

        { // q_no_swap --w--> ( q_id_1(L), q_id_1(R) )
             DLF left_subtree  { one * q_id_1 * Subtree_Tag::LEFT};
             DLF right_subtree { one * q_id_1 * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, working_qubit, transition);
        }

        { // q_swap --w--> ( q_id_1(R), q_id_1(L) )
             DLF left_subtree  { one * q_id_1 * Subtree_Tag::RIGHT};
             DLF right_subtree { one * q_id_1 * Subtree_Tag::LEFT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_swap, working_qubit, transition);
        }

        { // q_id_1 --w--> ( q_leaf(L), q_leaf(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::RIGHT};
             DLF right_subtree { one * q_leaf * Subtree_Tag::LEFT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_id_1, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_MAJ1) {
        State root   = 0;
        State q_id   = 1;
        State q_p1   = 2;  // Project away branch "1"
        State q_e1   = 3;  // Extract branch "1"
        State q_leaf = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(root);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // root --w--> ( q_id(L), q_id(R) )
             DLF left_subtree  { one * q_id * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_id * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(root, working_qubit, transition);
        }

        { // q_id --w--> ( q_p1(L) + q_e1(R), q_p1(R) + q_e1(L) )
             DLF left_subtree  { one * q_p1 * Subtree_Tag::LEFT, one * q_e1 * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q_p1 * Subtree_Tag::RIGHT, one * q_e1 * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_id, working_qubit, transition);
        }

        { // q_p1 --w--> ( q_leaf, 0 q_leaf )
             DLF left_subtree  { one  * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { zero * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_p1, working_qubit, transition);
        }

        { // q_e1 --w--> ( 0 q_leaf, q_leaf )
             DLF left_subtree  { zero * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { one  * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_e1, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_MAJ_RESULT) {
        State q_a = 0;
        State q_b = 1;
        State q_c = 2;
        State q_d = 3;
        State q_e = 4;
        State q_f = 5;
        State q_g = 6;
        State q_h = 7;
        State q_i = 8;
        State q_j = 9;

        u64 number_of_states = 10;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_a);
        builder.mark_state_final(q_j);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // q_a --w--> ( q_b(L) + q_c(R), q_e(L) + q_d(R) )
             DLF left_subtree  { one * q_b * Subtree_Tag::LEFT, one * q_c * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q_e * Subtree_Tag::LEFT, one * q_d * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_a, working_qubit, transition);
        }

        { // q_b --w--> ( q_f(L), q_f(R) )
             DLF left_subtree  { one * q_f * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_f * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_b, working_qubit, transition);
        }

        { // q_c --w--> ( q_h(R), q_h(L) )
             DLF left_subtree  { one * q_h * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q_h * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_c, working_qubit, transition);
        }

        { // q_d --w--> ( q_f(L), q_i(R) )
             DLF left_subtree  { one * q_f * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_i * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_d, working_qubit, transition);
        }

        { // q_e --w--> ( q_h(R), q_g(L) )
             DLF left_subtree  { one * q_h * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q_g * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_e, working_qubit, transition);
        }

        { // q_f --w--> ( q_j(L), 0 q_j(R) )
             DLF left_subtree  {  one * q_j * Subtree_Tag::LEFT };
             DLF right_subtree { zero * q_j * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_f, working_qubit, transition);
        }

        { // q_g --w--> ( q_j(L), 0 q_j(R) )
             DLF left_subtree  {  one * q_j * Subtree_Tag::RIGHT };
             DLF right_subtree { zero * q_j * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_g, working_qubit, transition);
        }

        { // q_h --w--> ( 0 q_j(L), q_j(R) )
             DLF left_subtree  { zero * q_j * Subtree_Tag::LEFT };
             DLF right_subtree { one  * q_j * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_h, working_qubit, transition);
        }

        { // q_i --w--> ( 0 q_j(R), q_j(L) )
             DLF left_subtree  { zero * q_j * Subtree_Tag::RIGHT };
             DLF right_subtree { one  * q_j * Subtree_Tag::LEFT  };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_i, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);

        do_on_debug({
             result.debug_data = new WTT::Debug_Data();
             result.debug_data->state_names[q_a] = "A";
             result.debug_data->state_names[q_b] = "B";
             result.debug_data->state_names[q_c] = "C";
             result.debug_data->state_names[q_d] = "D";
             result.debug_data->state_names[q_e] = "E";
             result.debug_data->state_names[q_f] = "F";
             result.debug_data->state_names[q_g] = "G";
             result.debug_data->state_names[q_h] = "H";
             result.debug_data->state_names[q_i] = "I";
             result.debug_data->state_names[q_j] = "J";
        });

        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_MAJ_RESULT_12) {
        State q0_a0  = 0;
        State q1_p0  = 1;
        State q1_p1  = 2;
        State p0_p0  = 3;
        State p1_p1  = 4;
        State q_leaf = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q0_a0);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // q0_a0 --w--> ( q1_p0(L) + q1_p1(R), q1_p0(R) + q1_p1(L) )
             DLF left_subtree  { one * q1_p0 * Subtree_Tag::LEFT, one * q1_p1 * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q1_p0 * Subtree_Tag::RIGHT, one * q1_p1 * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q0_a0, working_qubit, transition);
        }

        { // q1_p0 --w--> ( p0_p0(L), p0_p0(R) )
             DLF left_subtree  { one * p0_p0 * Subtree_Tag::LEFT  };
             DLF right_subtree { one * p0_p0 * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q1_p0, working_qubit, transition);
        }

        { // q1_p1 --w--> ( p0_p0(L), p0_p0(R) )
             DLF left_subtree  { one * p1_p1 * Subtree_Tag::RIGHT };
             DLF right_subtree { one * p1_p1 * Subtree_Tag::LEFT  };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q1_p1, working_qubit, transition);
        }

        { // p0_p0 --w--> ( q_leaf(L), 0 q_leaf(R) )
             DLF left_subtree  {  one * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { zero * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(p0_p0, working_qubit, transition);
        }

        { // p1_p1 --w--> ( q_leaf(L), 0 q_leaf(R) )
             DLF left_subtree  { zero * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { one  * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(p1_p1, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
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
        return result;
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
             builder.add_transition(q_h, sym_a, c, transition);
        }

        builder.add_bot_state_transitions(q_bot);

        std::vector<State> initial_states ({q_g});
        Bit_Set leaf_states (4, {q_bot, q_c});
        auto transition_fn = builder.build(4);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::BV_EXAMPLE_10STAR_RESULT || name == Predefined_SWTA_Names::TEST_BV_EXAMPLE_AFTER_STEP3) {
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

    if (name == Predefined_SWTA_Names::TEST_BV_EXAMPLE_AFTER_STEP1) {
        State q_alpha  = 0;
        State q_beta   = 1;

        Internal_Symbol working_qubit = 0, ancilla = 1;
        Color c = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        { // q_alpha --w--> (1/sqrt(2) q_alpha, 1/sqrt(2) q_alpha)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_alpha};
             DLF right_subtree {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_alpha};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_alpha, working_qubit, c, transition);
        }

        { // q_alpha --a--> (1/sqrt(2) q_beta, 1/sqrt(2) q_beta)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_beta};
             DLF right_subtree {Def_Coef(-ACN::ONE_OVER_SQRT2()) * q_beta};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_alpha, ancilla, c, transition);
        }

        std::vector<State> initial_states ({q_alpha});
        Bit_Set leaf_states (2, {q_beta});
        auto transition_fn = builder.build(2);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::TEST_BV_EXAMPLE_AFTER_STEP2) {
        State q_gamma   = 0;
        State q_delta   = 1;
        State q_epsilon = 2;
        State q_mu      = 3;
        State q_sigma   = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0, ancilla = 1;
        Color c = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        { // q_gamma --w--> (1/sqrt(2) q_delta, 1/sqrt(2) q_epsilon)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_delta};
             DLF right_subtree {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_epsilon};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_gamma, working_qubit, c, transition);
        }

        { // q_delta --w--> (1/sqrt(2) q_beta, 1/sqrt(2) q_beta)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_gamma};
             DLF right_subtree {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_gamma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_delta, working_qubit, c, transition);
        }

        { // q_epsilon --w--> (1/sqrt(2) q_mu, 1/sqrt(2) q_mu)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_mu};
             DLF right_subtree {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_mu};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_epsilon, working_qubit, c, transition);
        }

        { // q_mu --w--> (1/sqrt(2) q_epsilon, 1/sqrt(2) q_delta)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_epsilon};
             DLF right_subtree {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_delta};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_mu, working_qubit, c, transition);
        }

        { // q_gamma --a--> (1/sqrt(2) q_sigma, -1/sqrt(2) q_sigma)
             DLF left_subtree  {Def_Coef( ACN::ONE_OVER_SQRT2()) * q_sigma};
             DLF right_subtree {Def_Coef(-ACN::ONE_OVER_SQRT2()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_gamma, ancilla, c, transition);
        }

        { // q_delta --a--> (1/sqrt(2) q_sigma, -1/sqrt(2) q_sigma)
             DLF left_subtree  {Def_Coef( ACN::ONE_OVER_SQRT2()) * q_sigma};
             DLF right_subtree {Def_Coef(-ACN::ONE_OVER_SQRT2()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_delta, ancilla, c, transition);
        }

        { // q_epsilon --a--> (-1/sqrt(2) q_sigma, 1/sqrt(2) q_sigma)
             DLF left_subtree  {Def_Coef(-ACN::ONE_OVER_SQRT2()) * q_sigma};
             DLF right_subtree {Def_Coef( ACN::ONE_OVER_SQRT2()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_epsilon, ancilla, c, transition);
        }

        { // q_mu --a--> (-1/sqrt(2) q_sigma, 1/sqrt(2) q_sigma)
             DLF left_subtree  {Def_Coef(-ACN::ONE_OVER_SQRT2()) * q_sigma};
             DLF right_subtree {Def_Coef( ACN::ONE_OVER_SQRT2()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_mu, ancilla, c, transition);
        }

        std::vector<State> initial_states ({q_gamma});
        Bit_Set leaf_states (number_of_states, {q_sigma});
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::GROVER_ALL_BASIS) {
        State q_no_one_seen = 0;
        State q_all_zeros   = 1;
        State q_bot         = 2;
        State q_anc         = 3;
        State q_last_anc    = 4;
        State q_leaf        = 5;

        Internal_Symbol sym_w = 0;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // q_no_one_seen --w--> (q_no_one_seen, q_all_zeros)
             DLF left_subtree  {one * q_no_one_seen};
             DLF right_subtree {one * q_all_zeros};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_one_seen, working_qubit, color, transition);
        }

        { // q_no_one_seen --W--> (q_last_anc, 0 q_bot)
             DLF left_subtree  {one  * q_last_anc};
             DLF right_subtree {zero * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_one_seen, last_working_qubit, color, transition);
        }

        { // q_all_zeros --w--> (q_anc, 0 q_bot)
             DLF left_subtree  {one  * q_anc};
             DLF right_subtree {zero * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_all_zeros, working_qubit, color, transition);
        }

        { // q_anc --a--> (q_all_zeros, q_all_zeros)
             DLF left_subtree  {one * q_all_zeros};
             DLF right_subtree {one * q_all_zeros};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_anc, ancilla, color, transition);
        }

        { // q_all_zeros --W'--> (q_last_anc, 0 q_bot)
             DLF left_subtree  {one  * q_last_anc};
             DLF right_subtree {zero * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_all_zeros, last_working_qubit, color, transition);
        }

        { // q_last_anc --a--> (q_leaf, q_leaf)
             DLF left_subtree  {one * q_leaf};
             DLF right_subtree {one * q_leaf};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_last_anc, ancilla, color, transition);
        }

        builder.add_bot_state_transitions(q_bot);

        std::vector<State> initial_states ({q_no_one_seen});
        Bit_Set leaf_states (number_of_states, {q_bot, q_leaf});
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::TEST_ADDER_ALL_3BASIS) {
        State q3 = 0;
        State q2 = 1;
        State q1 = 2;
        State q0 = 3;

        u64 number_of_states = 4;

        Internal_Symbol working_qubit = 0;
        Color one_left  = 0;
        Color one_right = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 2 };
        SWTA::Transition_Builder builder (metadata);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // q3 --one_left--> (q2, 0q2)
             DLF left_subtree  {one  * q2};
             DLF right_subtree {zero * q2};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q3, working_qubit, one_left, transition);
        }

        { // q3 --one_right--> (0q2, q2)
             DLF left_subtree  {zero * q2};
             DLF right_subtree {one  * q2};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q3, working_qubit, one_right, transition);
        }

        { // q2 --one_left--> (q1, 0q1)
             DLF left_subtree  {one  * q1};
             DLF right_subtree {zero * q1};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q2, working_qubit, one_left, transition);
        }

        { // q2 --one_right--> (0q1, q1)
             DLF left_subtree  {zero * q1};
             DLF right_subtree {one  * q1};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q2, working_qubit, one_right, transition);
        }

        { // q1 --one_left--> (q0, 0 q0)
             DLF left_subtree  {one  * q0};
             DLF right_subtree {zero * q0};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q1, working_qubit, one_left, transition);
        }

        { // q1 --one_right--> (0 q0, q0)
             DLF left_subtree  {zero * q0};
             DLF right_subtree {one  * q0};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q1, working_qubit, one_right, transition);
        }

        std::vector<State> initial_states ({q3});
        Bit_Set leaf_states (number_of_states, { q0 });
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }


    throw std::runtime_error("No definition for the predefined SWTA: " + std::to_string(static_cast<u64>(name)));
}
