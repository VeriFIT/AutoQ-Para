#pragma once

#include "swta.hpp"

WTT makeProjectionWTT(const int qubitCount, const std::vector<int>& partialBasis, const std::vector<int>& mask);

WTT makeRotationWTT(
    const int qubitCount,
    const std::vector<int>& sourceBasis,
    const std::vector<int>& targetBasis,
    const std::vector<int>& mask);
