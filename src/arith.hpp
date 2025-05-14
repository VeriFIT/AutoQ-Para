#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <ostream>

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
            .a = this->a * (1u << half_scale),
            .b = this->b * (1u << half_scale),
            .c = this->c * (1u << half_scale),
            .d = this->d * (1u << half_scale),
            .k = larger_k
        };

        if (scale_difference % 2) { // Multiply by sqrt(2) if needed
            rescaled.a = rescaled.b - rescaled.d;
            rescaled.b = rescaled.a + rescaled.c;
            rescaled.c = rescaled.b + rescaled.d;
            rescaled.d = rescaled.c - rescaled.a;
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
