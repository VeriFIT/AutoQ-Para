#pragma once

#include "arith.hpp"
#include "swta.hpp"


enum class Predefined_WTT_Name : u64 {
    HADAMARD = 0,
    PARITY_CNOT = 1,
};

enum class Subtree_Tag : u64 {
    NONE = 0, // The transition DSL definition is not finish - its missing its subtree information
    LEFT = 1,
    RIGHT = 2,
};

struct Def_State;
struct Def_Coef;
struct Def_Linear_Form;

/**
 * Wrapper around WTT/SWTA state that allows easier definition
 * of linear forms/transitions.
 */
struct Def_State {
    State state;

    Def_State(State state) : state(state) {}
};

struct Def_Coef {
    using ACN = Algebraic_Complex_Number;

    ACN number;

    Def_Coef(const ACN& number) : number(number) {}

    Def_Linear_Form operator*(const Def_State& other);
};

struct Def_Linear_Form {
    using ACN = Algebraic_Complex_Number;

    ACN         coef;
    State       state;
    Subtree_Tag tag;

    Def_Linear_Form(const Def_Coef& def_coef, const Def_State& def_state) : coef(def_coef.number), state(def_state.state), tag(Subtree_Tag::NONE) {};

    Def_Linear_Form operator*(const Subtree_Tag& tag) {
        this->tag = tag;
        return *this;
    }
};



WTT get_predefined_wtt(Predefined_WTT_Name name);
