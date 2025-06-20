#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <gmp-x86_64.h>
#include <iostream>
#include <ostream>
#include <vector>

#include "gmp.h"

typedef int64_t  s64;
typedef uint8_t  u8;
typedef uint64_t u64;

struct Complex_Number {
    float real, im;
};

struct Fixed_Precision_ACN {
    s64 a, b, c, d, k;

    Fixed_Precision_ACN(s64 a, s64 b, s64 c, s64 d, s64 k) :
        a(a), b(b), c(c), d(d), k(k) {}

    bool operator==(const Fixed_Precision_ACN& other) const {
        return (a == other.a) && (b == other.b) && (c == other.c) && (d == other.d) && (k == other.k);
    }
};

std::ostream& operator<<(std::ostream& os, const Complex_Number& number);
std::ostream& operator<<(std::ostream& os, const Fixed_Precision_ACN& number);

/**
 * Algebraic representation of a complex number with arbitrary precision. Number C is represented
 * as a 5-tuple (a, b, c, d, k) such that: C = (1/sqrt(2))^k (a + bO + cO^2 + dO^3) where
 * O (Omega) is e^(PI*i/4).
 */
struct Algebraic_Complex_Number {
    mpz_t a, b, c, d, k;

    Algebraic_Complex_Number() {
        mpz_inits(a, b, c, d, 0);
        mpz_init_set_si(k, 0);
    }

    Algebraic_Complex_Number(s64 a, s64 b, s64 c, s64 d, s64 k = 1) {
        mpz_init_set_si(this->a, a);
        mpz_init_set_si(this->b, b);
        mpz_init_set_si(this->c, c);
        mpz_init_set_si(this->d, d);
        mpz_init_set_si(this->k, k);
    }

    Algebraic_Complex_Number(const Algebraic_Complex_Number& other) {
        mpz_init_set(this->a, other.a);
        mpz_init_set(this->b, other.b);
        mpz_init_set(this->c, other.c);
        mpz_init_set(this->d, other.d);
        mpz_init_set(this->k, other.k);
    }

    Algebraic_Complex_Number(Algebraic_Complex_Number&& other) {
        mpz_init_set(this->a, other.a);
        mpz_init_set(this->b, other.b);
        mpz_init_set(this->c, other.c);
        mpz_init_set(this->d, other.d);
        mpz_init_set(this->k, other.k);
    };

    ~Algebraic_Complex_Number() {
        mpz_clears(a, b, c, d, k, 0);
    }

    // Useful constants:
    static Algebraic_Complex_Number ONE_OVER_SQRT2() {
        return Algebraic_Complex_Number(1, 0, 0, 0, 1);
    }

    static Algebraic_Complex_Number ONE() {
        return Algebraic_Complex_Number(1, 0, 0, 0, 0);
    }

    static Algebraic_Complex_Number ZERO() {
        return Algebraic_Complex_Number(0, 0, 0, 0, 0);
    }

    Algebraic_Complex_Number operator*(const Algebraic_Complex_Number& other) const {

        Algebraic_Complex_Number result;

        mpz_t immediate;
        mpz_init(immediate);

        mpz_add(result.k, this->k, other.k);

        { // (this->a*other.a) - (this->b*other.d) - (this->c*other.c) - (this->d*other.b)
            mpz_mul(immediate, this->a, other.a);
            mpz_add(result.a, result.a, immediate);

            mpz_mul(immediate, this->b, other.d);
            mpz_sub(result.a, result.a, immediate);

            mpz_mul(immediate, this->c, other.c);
            mpz_sub(result.a, result.a, immediate);

            mpz_mul(immediate, this->d, other.b);
            mpz_sub(result.a, result.a, immediate);
        }

        { // (this->a*other.b) + (this->b*other.a) - (this->c*other.d) - (this->d*other.c)
            mpz_mul(immediate, this->a, other.b);
            mpz_add(result.b, result.b, immediate);

            mpz_mul(immediate, this->b, other.a);
            mpz_add(result.b, result.b, immediate);

            mpz_mul(immediate, this->c, other.d);
            mpz_sub(result.b, result.b, immediate);

            mpz_mul(immediate, this->d, other.c);
            mpz_sub(result.b, result.b, immediate);
        }

        { // (this->a*other.c) + (this->b*other.b) + (this->c*other.a) - (this->d*other.d)
            mpz_mul(immediate, this->a, other.c);
            mpz_add(result.c, result.c, immediate);

            mpz_mul(immediate, this->b, other.b);
            mpz_add(result.c, result.c, immediate);

            mpz_mul(immediate, this->c, other.a);
            mpz_add(result.c, result.c, immediate);

            mpz_mul(immediate, this->d, other.d);
            mpz_sub(result.c, result.c, immediate);
        }

        { // result_d = (this->a*other.d) + (this->b*other.c) + (this->c*other.b) + (this->d*other.a);
            mpz_mul(immediate, this->a, other.d);
            mpz_add(result.d, result.d, immediate);

            mpz_mul(immediate, this->b, other.c);
            mpz_add(result.d, result.d, immediate);

            mpz_mul(immediate, this->c, other.b);
            mpz_add(result.d, result.d, immediate);

            mpz_mul(immediate, this->d, other.a);
            mpz_add(result.d, result.d, immediate);
        }

        mpz_clear(immediate);

        return result;
    }

    Algebraic_Complex_Number rescale(const mpz_t larger_k) const {
        assert(mpz_cmp(larger_k, this->k) >= 0);

        s64 scale_difference_int = mpz_get_si(larger_k) - mpz_get_si(this->k);
        s64 half_scale_diff = scale_difference_int / 2;

        Algebraic_Complex_Number rescaled;

        mpz_mul_2exp(rescaled.a, this->a, half_scale_diff);
        mpz_mul_2exp(rescaled.b, this->b, half_scale_diff);
        mpz_mul_2exp(rescaled.c, this->c, half_scale_diff);
        mpz_mul_2exp(rescaled.d, this->d, half_scale_diff);

        if (scale_difference_int % 2) { // Multiply by sqrt(2) if needed
            mpz_t imm;
            mpz_init(imm);

            Algebraic_Complex_Number multiplied_by_sqrt2;

            // Use scale_difference as immediate value
            // COMPUTE: s64 new_a = -rescaled.b - rescaled.d;
            mpz_set_ui(imm, 0);
            mpz_sub(imm, imm, rescaled.b);
            mpz_sub(multiplied_by_sqrt2.a, imm, rescaled.d);

            // COMPUTE: s64 new_b = rescaled.a + rescaled.c;
            mpz_add(multiplied_by_sqrt2.b, rescaled.a, rescaled.c);

            // COMPUTE: s64 new_c = rescaled.b + rescaled.d;
            mpz_add(multiplied_by_sqrt2.c, rescaled.b, rescaled.d);

            // COMPUTE: s64 new_d = rescaled.c - rescaled.a;
            mpz_sub(multiplied_by_sqrt2.d, rescaled.c, rescaled.a);

            mpz_clear(imm);

            return multiplied_by_sqrt2;
        }

        return rescaled;
    }

    Algebraic_Complex_Number operator-() const {
        Algebraic_Complex_Number result;

        mpz_neg(result.a, this->a);
        mpz_neg(result.b, this->b);
        mpz_neg(result.c, this->c);
        mpz_neg(result.d, this->d);
        mpz_set(result.k, this->k);

        return result;
    }

    Algebraic_Complex_Number operator+(const Algebraic_Complex_Number& other) const {
        const Algebraic_Complex_Number* smaller = this;
        const Algebraic_Complex_Number* larger  = &other;
        if (mpz_cmp(larger->k, smaller->k) >= 1) {
            std::swap(smaller, larger);
        }

        Algebraic_Complex_Number larger_rescaled = larger->rescale(smaller->k);

        Algebraic_Complex_Number result;

        mpz_add(result.a, smaller->a, larger_rescaled.a);
        mpz_add(result.b, smaller->b, larger_rescaled.b);
        mpz_add(result.c, smaller->c, larger_rescaled.c);
        mpz_add(result.d, smaller->d, larger_rescaled.d);
        mpz_set(result.k, smaller->k);

        return result;
    }

    Algebraic_Complex_Number operator-(const Algebraic_Complex_Number& other) const {
        auto other_negated = -other;
        return *this + other_negated;
    }

    void operator=(const Algebraic_Complex_Number& other) {
        mpz_set(this->a, other.a);
        mpz_set(this->b, other.b);
        mpz_set(this->c, other.c);
        mpz_set(this->d, other.d);
        mpz_set(this->k, other.k);
    }

    void operator+=(const Algebraic_Complex_Number& other) {
        auto addition_result = *this + other;

        mpz_set(this->a, addition_result.a);
        mpz_set(this->b, addition_result.b);
        mpz_set(this->c, addition_result.c);
        mpz_set(this->d, addition_result.d);
        mpz_set(this->k, addition_result.k);
    }

    bool operator==(const Algebraic_Complex_Number& other) const {
        return (mpz_cmp(this->a, other.a) == 0) &&
               (mpz_cmp(this->b, other.b) == 0) &&
               (mpz_cmp(this->c, other.c) == 0) &&
               (mpz_cmp(this->d, other.d) == 0) &&
               (mpz_cmp(this->k, other.k) == 0);
    }

    bool operator!=(const Algebraic_Complex_Number& other) const {
        return !(*this == other);
    }

    bool is_zero() const {
        return (mpz_cmp_ui(this->a, 0) == 0) &&
               (mpz_cmp_ui(this->b, 0) == 0) &&
               (mpz_cmp_ui(this->c, 0) == 0) &&
               (mpz_cmp_ui(this->d, 0) == 0);
    }

    void swap(Algebraic_Complex_Number& other) noexcept {
        mpz_swap(this->a, other.a);
        mpz_swap(this->b, other.b);
        mpz_swap(this->c, other.c);
        mpz_swap(this->d, other.d);
        mpz_swap(this->k, other.k);
    }

    Complex_Number into_approx() const {
        float one_over_sqrt2 = 1.0f / std::sqrt(2);

        float real = static_cast<float>(mpz_get_si(this->a)) + one_over_sqrt2*(static_cast<float>(mpz_get_si(b)) - static_cast<float>(mpz_get_si(this->d)));
        float im   = static_cast<float>(mpz_get_si(this->c)) + one_over_sqrt2*(static_cast<float>(mpz_get_si(b)) + static_cast<float>(mpz_get_si(this->d)));

        return {.real = real, .im = im};
    }

    Fixed_Precision_ACN into_fixed_precision() const {
        return Fixed_Precision_ACN(
            mpz_get_si(this->a),
            mpz_get_si(this->b),
            mpz_get_si(this->c),
            mpz_get_si(this->d),
            mpz_get_si(this->k)
        );
    }

    void normalize() {
        if (this->is_zero()) {
            mpz_set_ui(this->k, 0);
        }

        u64 a_2pow = mpz_scan1(this->a, 0);
        u64 b_2pow = mpz_scan1(this->b, 0);
        u64 c_2pow = mpz_scan1(this->c, 0);
        u64 d_2pow = mpz_scan1(this->d, 0);

        u64 spare_2pow;
        {
            u64 spare_2pow_a = std::min(a_2pow, b_2pow);
            u64 spare_2pow_b = std::min(c_2pow, d_2pow);
            spare_2pow = std::min(spare_2pow_a, spare_2pow_b);
        }

        if (spare_2pow == 0) return; // There is no common power of 2 dividing all of the components (multiplying omega)

        s64 available_k = mpz_get_ui(this->k) / 2;
        if (available_k < 0) return;

        u64 normalization_2pow = std::min(spare_2pow, static_cast<u64>(available_k));

        // Divide (a, b, c, d) by the agreed power of 2
        mpz_div_2exp(this->a, this->a, normalization_2pow);
        mpz_div_2exp(this->b, this->b, normalization_2pow);
        mpz_div_2exp(this->c, this->c, normalization_2pow);
        mpz_div_2exp(this->d, this->d, normalization_2pow);

        // Subtract from k what has been used in in the process
        mpz_sub_ui(this->k, this->k, 2*normalization_2pow);
    }
};

std::ostream& operator<<(std::ostream& os, const Algebraic_Complex_Number& number);

Algebraic_Complex_Number acn_zero();

/**
 * Represents a complex number as:
 *    (1/2)^k * (a + b*(1/sqrt(2)) + i*c + i*d*(1/sqrt(2)))
 * Note: The scaling factor k is necessary, without it we cannot represent arbitrary small/large numbers
 */
struct Direct_ACN {
    s64 a, b, c, d, k;

    bool operator==(const Direct_ACN& other) const {
        return (this->a == other.a) && (this->b == other.b) && (this->c == other.c) && (this->d == other.d) && (this->k == other.k);
    }
};

std::ostream& operator<<(std::ostream& os, const Direct_ACN& number);
Direct_ACN convert_acn_into_direct_repr(const Algebraic_Complex_Number& number);


struct ACN_Matrix {
    u64 height, width;
    Algebraic_Complex_Number* data = nullptr;

    ACN_Matrix(u64 height, u64 width, Algebraic_Complex_Number* data_ptr = nullptr) :
        height(height), width(width), data(data_ptr)
    {
        if (this->data == nullptr) {
            this->data = new Algebraic_Complex_Number[this->width*this->height]; // Zero initialized
        }
    }

    ACN_Matrix(const ACN_Matrix& other) {
        this->width  = other.width;
        this->height = other.height;

        this->data   = new Algebraic_Complex_Number[this->width*this->height];

        for (u64 idx = 0; idx < other.width*other.height; idx++) {
            this->data[idx] = other.data[idx];
        }
    }

    ACN_Matrix(ACN_Matrix&& other) {
        this->width  = other.width;
        this->height = other.height;

        this->data = other.data;
        other.data = nullptr;
    }

    ACN_Matrix operator*(const ACN_Matrix& other) const {
        assert(this->width == other.height);

        u64 result_height = this->height;
        u64 result_width  = other.width;

        Algebraic_Complex_Number* result_data = new Algebraic_Complex_Number[result_height*result_width];
        u64 target_cell_idx = 0;

        for (u64 row_idx = 0; row_idx < this->height; row_idx++) {
            for (u64 col_idx = 0; col_idx < other.width; col_idx++) {

                Algebraic_Complex_Number dot_product;
                for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
                    Algebraic_Complex_Number fragment = this->at(row_idx, elem_idx) * other.at(elem_idx, col_idx);
                    dot_product += fragment;
                }

                result_data[target_cell_idx] = dot_product;
                target_cell_idx += 1;
            }
        }

        return ACN_Matrix(result_height, result_width, result_data);
    }

    Algebraic_Complex_Number& at(u64 row_idx, u64 col_idx) const {
        return this->data[row_idx*this->width + col_idx];
    }

    void set(u64 row_idx, u64 col_idx, const Algebraic_Complex_Number& value) {
        this->data[row_idx*this->width + col_idx] = value;
    }

    u64 find_nonzero_elem_in_row(u64 row_idx) const {
        assert (row_idx < this->height);

        for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
            Algebraic_Complex_Number& elem = this->at(row_idx, elem_idx);
            if (!elem.is_zero()) {
                return elem_idx;
            }
        }

        return this->width;
    }

    void insert_row_at(const ACN_Matrix& row, u64 row_idx, bool skip_shifting_subsequent_rows = false) {
        assert(row.width == this->width);

        if (!skip_shifting_subsequent_rows) {
            s64 last_elem_idx  = (this->height - 1) * this->width - 1;
            s64 first_elem_idx = row_idx * this->width;

            for (s64 elem_idx = last_elem_idx; elem_idx >= first_elem_idx; elem_idx--) {
                this->data[elem_idx + this->width] = this->data[elem_idx];
            }
        }

        for (u64 elem_idx = 0; elem_idx < row.width; elem_idx++) {
            this->data[row_idx*this->width + elem_idx] = row.at(0, elem_idx);
        }
    }

    void subtract_from_ith_row(u64 row_idx, Algebraic_Complex_Number& row_coef, ACN_Matrix& rows_to_subtract, u64 row_to_subtract_idx, Algebraic_Complex_Number& row_to_subtract_coef) const {
        for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
            auto subtractee_weighted = this->at(row_idx, elem_idx) * row_coef;
            auto subtractor_weighted = rows_to_subtract.at(row_to_subtract_idx, elem_idx) * row_to_subtract_coef;
            auto result_elem = subtractee_weighted - subtractor_weighted;

            this->data[row_idx*this->width + elem_idx] = result_elem;
        }
    }

    bool contains_only_zeros() const {
        for (u64 idx = 0; idx < this->width*this->height; idx++) {
            if (!this->data[idx].is_zero()) return false;
        }
        return true;
    }

    bool operator==(const ACN_Matrix& other) const {
        if (width != other.width || height != other.height) return false;

        for (u64 idx = 0; idx < height*width; idx++) {
            if (data[idx] != other.data[idx]) return false;
        }

        return true;
    }

    ~ACN_Matrix() {
        if (this->data != nullptr) delete[] this->data;
    }
};

std::ostream& operator<<(std::ostream& os, const ACN_Matrix& matrix);

s64 add_row_to_row_echelon_matrix(ACN_Matrix& matrix, ACN_Matrix& row);
ACN_Matrix square_acn_matrix_from_ints(const std::vector<s64>& ints);
ACN_Matrix row_from_ints(const std::vector<s64>& row_data);
ACN_Matrix column_from_ints(const std::vector<s64>& column_data);
