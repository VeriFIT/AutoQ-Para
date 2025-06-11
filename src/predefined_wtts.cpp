#include "predefined_wtts.hpp"
#include "arith.hpp"
#include <stdexcept>
#include <string>

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

        WTT transducer ({transition}, {0}, {0});
        return transducer;
    }

    throw std::runtime_error("Unknown WTT. " + std::to_string(static_cast<u64>(name)));
}
