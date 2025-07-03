#include "quantum_program.hpp"
#include "swta.hpp"


SWTA run_swta_program(const SWTA_Program &program) {
    SWTA imm = program.initial_swta;

    for (auto& application : program.applications) {
        const auto& transducer = program.transducers[application.transducer_handle];
        imm = apply_wtt_to_swta(imm, transducer);
    }

    return imm;
}
