#pragma once

#include "arith.hpp"
#include "swta.hpp"


enum class Predefined_WTT_Names : u64 {
    HADAMARD    = 0,
    PARITY_CNOT = 1,

    GROVER_FIRST_MULTI_Z  =  8,
    GROVER_SECOND_MULTI_Z =  9,
    GROVER_X              = 10,
    GROVER_H              = 11,
    GROVER_FIRST_MULTI_Z_USING_CCX  = 12,
    GROVER_SECOND_MULTI_Z_USING_CCX = 13,
};

enum class Predefined_SWTA_Names : u64 {
    BV_EXAMPLE_10STAR_PRE = 0,
    BV_EXAMPLE_10STAR_POST = 1,
    BV_EXAMPLE_10STAR_RESULT = 2,

    TRIVIAL_BOT = 3,   // Used for tests
    TRIVIAL_ONES = 4,  // Used for tests
    TEST_BV_EXAMPLE_AFTER_STEP1 = 5,
    TEST_BV_EXAMPLE_AFTER_STEP2 = 6,
    TEST_BV_EXAMPLE_AFTER_STEP3 = 7,

    GROVER_ALL_BASIS = 8,
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



WTT get_predefined_wtt(Predefined_WTT_Names name, const SWTA::Metadata& metadata);

SWTA get_predefined_swta(Predefined_SWTA_Names name);

SWTA::Transition synthetize_swta_transition(const std::vector<Def_Linear_Form>& left_subtree, const std::vector<Def_Linear_Form>& right_subtree);

WTT::Transition synthetize_wtt_transition(const std::vector<Def_Linear_Form>& left_subtree, const std::vector<Def_Linear_Form>& right_subtree);
