#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <ostream>
#include <vector>

typedef int64_t  s64;
typedef uint64_t u64;

struct Complex_Number {
    float real, im;
};

std::ostream& operator<<(std::ostream& os, const Complex_Number& number);

struct Algebraic_Complex_Number {
    s64 a, b, c, d, k;

    Algebraic_Complex_Number operator*(const Algebraic_Complex_Number& other) const {
        s64 result_k = this->k + other.k;
        s64 result_a = (this->a*other.a) - (this->b*other.d) - (this->c*other.c) - (this->d*other.b);
        s64 result_b = (this->a*other.b) + (this->b*other.a) - (this->c*other.d) - (this->d*other.c);
        s64 result_c = (this->a*other.c) + (this->b*other.b) + (this->c*other.a) - (this->d*other.d);
        s64 result_d = (this->a*other.d) + (this->b*other.c) + (this->c*other.b) + (this->d*other.a);

        return {
            .a = result_a,
            .b = result_b,
            .c = result_c,
            .d = result_d,
            .k = result_k
        };
    }

    Algebraic_Complex_Number rescale(s64 larger_k) const {
        assert(larger_k >= this->k);

        s64 scale_difference = larger_k - this->k;
        s64 half_scale = scale_difference / 2;

        Algebraic_Complex_Number rescaled = { // First multiply by (sqrt(2))^2k
            .a = this->a << half_scale,
            .b = this->b << half_scale,
            .c = this->c << half_scale,
            .d = this->d << half_scale,
            .k = larger_k
        };

        if (scale_difference % 2) { // Multiply by sqrt(2) if needed
            s64 new_a = -rescaled.b - rescaled.d;
            s64 new_b = rescaled.a + rescaled.c;
            s64 new_c = rescaled.b + rescaled.d;
            s64 new_d = rescaled.c - rescaled.a;

            rescaled.a = new_a;
            rescaled.b = new_b;
            rescaled.c = new_c;
            rescaled.d = new_d;
        }

        return rescaled;
    }

    Algebraic_Complex_Number negate() const {
        Algebraic_Complex_Number result = {
            .a = -this->a,
            .b = -this->b,
            .c = -this->c,
            .d = -this->d,
            .k = this->k
        };
        return result;
    }

    Algebraic_Complex_Number operator+(const Algebraic_Complex_Number& other) const {
        const Algebraic_Complex_Number* smaller = this;
        const Algebraic_Complex_Number* larger  = &other;
        if (larger->k > smaller->k) {
            std::swap(smaller, larger);
        }

        Algebraic_Complex_Number larger_rescaled = larger->rescale(smaller->k);

        Algebraic_Complex_Number result = {
            .a = smaller->a + larger_rescaled.a,
            .b = smaller->b + larger_rescaled.b,
            .c = smaller->c + larger_rescaled.c,
            .d = smaller->d + larger_rescaled.d,
            .k = smaller->k
        };

        return result;
    }

    Algebraic_Complex_Number operator-(const Algebraic_Complex_Number& other) const {
        auto other_negated = other.negate();
        return *this + other.negate();
    }

    void operator+=(const Algebraic_Complex_Number& other) {
        auto addition_result = *this + other;

        this->a = addition_result.a;
        this->b = addition_result.b;
        this->c = addition_result.c;
        this->d = addition_result.d;
        this->k = addition_result.k;
    }

    bool is_zero() const {
        return (this->a == 0) && (this->b == 0) && (this->c == 0) && (this->d == 0);
    }

    Complex_Number into_approx() {
        float one_over_sqrt2 = 1.0f / std::sqrt(2);

        float real = static_cast<float>(this->a) + one_over_sqrt2*(static_cast<float>(b) - static_cast<float>(this->d));
        float im   = static_cast<float>(this->c) + one_over_sqrt2*(static_cast<float>(b) + static_cast<float>(this->d));

        return {.real = real, .im = im};
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

    ACN_Matrix(u64 height, u64 width, Algebraic_Complex_Number* data_ptr = nullptr) : height(height), width(width), data(data_ptr) {}

    ACN_Matrix(const ACN_Matrix& other) {
        this->width  = other.width;
        this->height = other.height;

        this->data   = new Algebraic_Complex_Number[this->width*this->height];
        std::memcpy(this->data, other.data, this->width*this->height*sizeof(Algebraic_Complex_Number));
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

                Algebraic_Complex_Number dot_product = acn_zero();
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

    ~ACN_Matrix() {
        if (this->data != nullptr) delete[] this->data;
    }
};

std::ostream& operator<<(std::ostream& os, const ACN_Matrix& matrix);

s64 add_row_to_row_echelon_matrix(ACN_Matrix& matrix, ACN_Matrix& row);
ACN_Matrix square_acn_matrix_from_ints(const std::vector<s64>& ints);
ACN_Matrix row_from_ints(const std::vector<s64>& row_data);
