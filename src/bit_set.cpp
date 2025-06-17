#include "bit_set.hpp"

std::vector<u64> Bit_Set::into_vector() const {
    std::vector<u64> bit_set_contents;

    for (size_t bit_idx = 0; bit_idx < this->size; bit_idx++) {
        if (this->get_bit_value(bit_idx)) {
            bit_set_contents.push_back(bit_idx);
        }
    }

    return bit_set_contents;
}
