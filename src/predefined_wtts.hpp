#pragma once

#include "swta.hpp"


enum class Predefined_WTT_Name : u64 {
    HADAMARD = 0,
};


WTT get_predefined_wtt(Predefined_WTT_Name name);
