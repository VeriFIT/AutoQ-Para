#include "predefined_wtts.hpp"
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

WTT get_predefined_wtt(Predefined_WTT_Name name) {
    if (name == Predefined_WTT_Name::HADAMARD) {
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

    if (name == Predefined_WTT_Name::PARITY_CNOT) {
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
