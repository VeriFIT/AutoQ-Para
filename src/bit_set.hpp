#pragma once

#include "arith.hpp"
#include <cstring>

struct Bit_Set {
    u64  size; // Number of *bits* stored
    u64* data;

    Bit_Set(u64 size, u64* data_ptr = nullptr) : size(size), data(data_ptr) {
        if (this->data == nullptr) {
            u64 bucket_count = size / sizeof(u64);
            bucket_count += (size % sizeof(u64)) > 0;
            this->data = new u64[bucket_count];

            memset(this->data, 0, sizeof(u64)*bucket_count);
        }
    }

    Bit_Set(const Bit_Set& other) {
        this->size = other.size;
        this->data = new u64[other.get_bucket_count()];
        std::memcpy(this->data, other.data, sizeof(u64) * other.get_bucket_count());
    }

    bool is_intersection_empty(const Bit_Set& other) const {
        u64 min_size = std::min(this->size, other.size);

        u64 chunks_to_compare = min_size / sizeof(u64);

        for (u64 chunk_idx = 0; chunk_idx < chunks_to_compare; chunk_idx++) {
            if (this->data[chunk_idx] & other.data[chunk_idx]) return false;
        }


        u64 spilled_size = min_size % sizeof(u64);
        if (spilled_size == 0) return true;


        u64 spilled_contents = this->data[chunks_to_compare] & other.data[chunks_to_compare];
        spilled_contents = spilled_contents & ~((~0) << spilled_size); // Extract only the relevant bits

        return spilled_contents == 0; // The result of and is 0, meaning that the intersection is empty
    }

    bool is_superset(const Bit_Set& other) const {
        u64 min_size = std::min(this->size, other.size);

        u64 chunks_to_compare = min_size / sizeof(u64);

        for (u64 chunk_idx = 0; chunk_idx < chunks_to_compare; chunk_idx++) {
            u64 interleaved = this->data[chunk_idx] & other.data[chunk_idx];

            // There is a bit set in other chunk that is not set in this
            if (interleaved != this->data[chunk_idx]) return false;
        }

        u64 spilled_size = min_size % sizeof(u64);
        if (spilled_size == 0) return true;

        u64 mask =~((~0) << spilled_size); // Extract only the relevant bits
        u64 spilled_interleaved = this->data[chunks_to_compare] | other.data[chunks_to_compare];

        spilled_interleaved = spilled_interleaved & mask;

        return spilled_interleaved == (this->data[chunks_to_compare] & mask);
    }

    u64 calc_target_bucket(u64 bit_idx) {
        u64 bucket = bit_idx / sizeof(u64);
        return bucket;
    }

    void set_bit(u64 bit_idx, bool value) {
        assert (bit_idx < this->size);

        u64 target_bucket_idx = this->calc_target_bucket(bit_idx);
        u64 bucket_offset     = bit_idx % sizeof(u64);

        if (value) {
            this->data[target_bucket_idx] = this->data[target_bucket_idx] | (1u << bucket_offset);
        } else {
            this->data[target_bucket_idx] = this->data[target_bucket_idx] & ~(1u << bucket_offset);
        }
    }

    void set_all(bool value) {
        u64 bucket_value = value ? (~0) : 0;

        u64 bucket_cnt = this->get_bucket_count();
        for (u64 bucket_idx = 0; bucket_idx < bucket_cnt; bucket_idx++) {
            this->data[bucket_idx] = bucket_value;
        }
    }

    void clear() {
        this->set_all(false);
    }

    u64 get_bucket_count() const {
        u64 bucket_count = size / sizeof(u64);
        bucket_count += (size % sizeof(u64)) > 0;

        return bucket_count;
    }

    ~Bit_Set() {
        if (this->data != nullptr) delete[] this->data;
        this->data = nullptr;
    }
};
