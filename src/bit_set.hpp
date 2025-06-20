#pragma once

#include "arith.hpp"

#include <cstring>

#include <bit>  // For popcount

struct Bit_Set {
    u64  size; // Number of *bits* stored
    u64* data;

    Bit_Set(u64 size, const std::vector<u64>& content) : size(size), data(nullptr) {
        u64 bucket_count = this->get_bucket_count();
        this->data = new u64[bucket_count];

        memset(this->data, 0, sizeof(u64)*bucket_count);

        for (auto& bit_to_set : content) {
            this->set_bit(bit_to_set, true);
        }
    }

    Bit_Set(u64 size, u64* data_ptr = nullptr) : size(size), data(data_ptr) {
        if (this->data == nullptr && this->size > 0) {
            u64 bucket_count = size / bits_in_bucket();
            bucket_count += (size % bits_in_bucket()) > 0;
            this->data = new u64[bucket_count];

            memset(this->data, 0, sizeof(u64)*bucket_count);
        }

        if (this->size == 0) this->data = nullptr;
    }

    Bit_Set(const Bit_Set& other) {
        this->size = other.size;
        this->data = new u64[other.get_bucket_count()];
        std::memcpy(this->data, other.data, sizeof(u64) * other.get_bucket_count());
    }

    Bit_Set(Bit_Set&& other) {
        this->size = other.size;
        this->data = other.data;

        other.data = nullptr;
        other.size = 0;
    }

    Bit_Set& operator=(const Bit_Set& other) {
        this->size = other.size;
        this->data = new u64[other.get_bucket_count()];
        std::memcpy(this->data, other.data, sizeof(u64) * other.get_bucket_count());

        return *this;
    }

    Bit_Set& operator=(Bit_Set&& other) {
        this->size = other.size;
        other.size = 0;

        this->data = other.data;
        other.data = nullptr;

        return *this;
    }

    u64 bits_in_bucket() const {
        return sizeof(u64) * 8;
    }

    bool operator==(const Bit_Set& other) const {
        if (this->size != other.size) return false;

        u64 bucket_cnt = this->get_bucket_count();

        for (u64 bucket_idx = 0; bucket_idx < bucket_cnt; bucket_idx++) {
            if (this->data[bucket_idx] != other.data[bucket_idx]) return false;
        }

        return true;
    }

    u64 popcount() const {
        u64 popcount = 0;
        for (u64 bucket_idx = 0; bucket_idx < this->get_bucket_count(); bucket_idx++) {
            u64 bucket_val = this->data[bucket_idx];
            popcount += std::popcount(bucket_val);
        }
        return popcount;
    }

    bool are_all_bits_set() const {
        if (this->get_bucket_count() == 0) return true;

        for (u64 bucket_idx = 0; bucket_idx < this->get_bucket_count() - 1; bucket_idx++) {
            if (this->data[bucket_idx] != ~static_cast<u64>(0)) return false;
        }

        u64 last_bucket = this->data[this->get_bucket_count() - 1];

        u64 last_valid_bit_offset = (this->size - 1) % bits_in_bucket();
        u64 bits_to_zero_out = bits_in_bucket() - (last_valid_bit_offset + 1);

        u64 all_ones = ~static_cast<u64>(0);
        u64 mask = (all_ones >> bits_to_zero_out);

        return last_bucket == mask;
    }

    std::vector<u64> into_vector() const;

    /**
     * Lexicographical ordering on the underlying bit representation
     */
    bool operator<(const Bit_Set& other) const {
        assert(this->size == other.size);

        for (u64 bucket_idx = 0; bucket_idx < this->get_bucket_count(); bucket_idx++) {
            if (this->data[bucket_idx] < other.data[bucket_idx]) return true;
            if (this->data[bucket_idx] > other.data[bucket_idx]) return false;
        }
        return false; // All buckets have equal value
    }

    bool is_intersection_empty(const Bit_Set& other) const {
        u64 min_size = std::min(this->size, other.size);

        u64 chunks_to_compare = min_size / bits_in_bucket();

        for (u64 chunk_idx = 0; chunk_idx < chunks_to_compare; chunk_idx++) {
            if (this->data[chunk_idx] & other.data[chunk_idx]) return false;
        }

        u64 spilled_size = min_size % bits_in_bucket();
        if (spilled_size == 0) return true;


        u64 spilled_contents = this->data[chunks_to_compare] & other.data[chunks_to_compare];
        spilled_contents = spilled_contents & ~((~0) << spilled_size); // Extract only the relevant bits

        return spilled_contents == 0; // The result of and is 0, meaning that the intersection is empty
    }

    bool is_superset(const Bit_Set& other) const {
        u64 min_size = std::min(this->size, other.size);

        u64 chunks_to_compare = min_size / bits_in_bucket();

        for (u64 chunk_idx = 0; chunk_idx < chunks_to_compare; chunk_idx++) {
            u64 interleaved = this->data[chunk_idx] & other.data[chunk_idx];

            // There is a bit set in other chunk that is not set in this
            if (interleaved != this->data[chunk_idx]) return false;
        }

        u64 spilled_size = min_size % bits_in_bucket();
        if (spilled_size == 0) return true;

        u64 mask =~((~0) << spilled_size); // Extract only the relevant bits
        u64 spilled_interleaved = this->data[chunks_to_compare] | other.data[chunks_to_compare];

        spilled_interleaved = spilled_interleaved & mask;

        return spilled_interleaved == (this->data[chunks_to_compare] & mask);
    }

    u64 calc_target_bucket(u64 bit_idx) const {
        u64 bucket = bit_idx / bits_in_bucket();
        return bucket;
    }

    void set_bit(u64 bit_idx, bool value = true) {
        assert (bit_idx < this->size);

        u64 target_bucket_idx = this->calc_target_bucket(bit_idx);
        u64 bucket_offset     = bit_idx % bits_in_bucket();

        if (value) {
            this->data[target_bucket_idx] = this->data[target_bucket_idx] | (1u << bucket_offset);
        } else {
            this->data[target_bucket_idx] = this->data[target_bucket_idx] & ~(1u << bucket_offset);
        }
    }

   void grow(u64 new_size) {
        if (new_size <= size) return;

        u64 target_bucket_idx = this->calc_target_bucket(new_size - 1);

        if (target_bucket_idx < this->get_bucket_count()) { // We have enough buckets, no need to alloc
            this->size = new_size; // Note that we can accomodate more bits now
            return;
        }

        // There is not enough space; we need to grow the underlying storage
        u64 old_bucket_cnt = this->get_bucket_count();
        this->size = new_size;
        u64 new_bucket_cnt = this->get_bucket_count();

        u64* old_data = this->data;
        u64* new_data = new u64[new_bucket_cnt];

        for (u64 i = 0; i < old_bucket_cnt; i++) {
            new_data[i] = this->data[i];
        }

        for (u64 i = old_bucket_cnt; i < new_bucket_cnt; i++) {
            new_data[i] = 0;
        }

        this->data = new_data;
        delete[] old_data;
    }

    void grow_and_set_bit(u64 bit_idx, bool bit_value = true) {
        u64 required_size = bit_idx + 1;
        this->grow(required_size);
        this->set_bit(bit_idx, bit_value);
    }

    bool get_bit_value(u64 bit_idx) const {
        u64 target_bucket_idx = this->calc_target_bucket(bit_idx);
        u64 bucket_offset     = bit_idx % bits_in_bucket();

        u64 bucket_value = this->data[target_bucket_idx];
        u64 bit_value = (bucket_value & (1u << bucket_offset)) > 0;

        return bit_value;
    }

    void set_all(bool value = true) {
        u64 bucket_value = value ? (~0) : 0;

        u64 bucket_cnt = this->get_bucket_count();
        for (u64 bucket_idx = 0; bucket_idx < bucket_cnt; bucket_idx++) {
            this->data[bucket_idx] = bucket_value;
        }

        u64 last_bucket_idx = bucket_cnt - 1;
        u64 last_valid_bit_offset = (this->size - 1) % bits_in_bucket();
        u64 bits_to_zero_out = bits_in_bucket() - (last_valid_bit_offset + 1);

        this->data[last_bucket_idx] = (this->data[last_bucket_idx] >> bits_to_zero_out);
    }

    void clear() {
        this->set_all(false);
    }

    u64 get_bucket_count() const {
        u64 bucket_count = size / bits_in_bucket();
        bucket_count += (size %  bits_in_bucket()) > 0;

        return bucket_count;
    }

    ~Bit_Set() {
        if (this->data != nullptr) delete[] this->data;
        this->data = nullptr;
    }
};
