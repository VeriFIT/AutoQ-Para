#pragma once
#include "basics.hpp"
#include "swta.hpp"


struct Transducer_Application {
    /**
     * The id (index into the vector of synthetized transducers) of the transducer to apply.
     */
    u16 transducer_handle;
};

struct SWTA_Program {
    SWTA initial_swta;
    std::vector<WTT>          transducers;
    std::vector<Transducer_Application> applications;
};

SWTA run_swta_program(const SWTA_Program& program);
