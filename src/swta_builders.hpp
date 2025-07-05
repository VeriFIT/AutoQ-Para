#pragma once

#include "basics.hpp"
#include "swta.hpp"


WTT perform_staircase_construction(WTT& box, const std::vector<Internal_Symbol>& box_inputs, u64 box_offset, u64 terminating_symbol);
