#pragma once

#include "basics.hpp"
#include "swta.hpp"


enum Staircase_Direction : u64 {
    LEFT_RIGHT = 0,
    RIGHT_LEFT = 1,
};

WTT perform_staircase_construction(WTT& box, const std::vector<Internal_Symbol>& box_inputs, u64 box_offset, u64 terminating_symbol, Staircase_Direction direction = Staircase_Direction::LEFT_RIGHT);
