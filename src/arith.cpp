#include <cassert>
#include <cmath>
#include <cstdlib>
#include <new>
#include <ostream>
#include <iostream>
#include <vector>

#include "arith.hpp"


std::ostream& operator<<(std::ostream& os, const Complex_Number& number) {
    const char* sign = number.im > 0 ? "+" : "-";

    os << number.real << " "
       << sign
       << std::abs(number.im) << " i";

    return os;
}

Algebraic_Complex_Number acn_zero() {
    return Algebraic_Complex_Number {
        .a = 0, .b = 0, .c = 0, .d = 0, .k = 0
    };
}

std::ostream& operator<<(std::ostream& os, const Algebraic_Complex_Number& number) {
    os << "("
       << number.a << ", "
       << number.b << ", "
       << number.c << ", "
       << number.d << ", "
       << number.k << ")";
    return os;
}


s64 add_row_to_row_echelon_matrix(ACN_Matrix& matrix, ACN_Matrix& row) {
    assert(matrix.width == row.width);
    assert(row.height == 1);

    u64 inserted_row_pivot_idx = row.find_nonzero_elem_in_row(0);

    s64 row_inserted_at = -1; // Not inserted

    for (u64 row_idx = 0; row_idx < matrix.height; row_idx++) {
        u64 row_pivot_idx = matrix.find_nonzero_elem_in_row(row_idx);

        if (row_pivot_idx < inserted_row_pivot_idx) {
            continue;
        }

        if (row_pivot_idx > inserted_row_pivot_idx) {
            matrix.insert_row_at(row, row_idx);
            row_inserted_at = row_idx;
            return row_inserted_at;
        }

        // Subtract the current matrix row and continue;
        auto& matrix_pivot = matrix.at(row_idx, row_pivot_idx);
        auto& row_pivot    = row.at(0, inserted_row_pivot_idx);
        row.subtract_from_ith_row(0, matrix_pivot, matrix, row_idx, row_pivot);
        inserted_row_pivot_idx = row.find_nonzero_elem_in_row(0);

        if (inserted_row_pivot_idx >= row.width) {
            return -1; // Not inserted, Gaussian elimination reduced the row to 0
        }
    }

    assert(false);  // Should be unreachable - the row is indepentent and we insert it, or our matrix has full rank and hence we do not insert it
}

u64 compute_square_matrix_dim_from_1d_repr(const std::vector<s64>& repr) {
    u64 dimension = static_cast<s64>(std::sqrt(repr.size()));

    assert(dimension*dimension == repr.size());

    return dimension;
}

ACN_Matrix square_acn_matrix_from_ints(const std::vector<s64>& ints) {
    u64 dim = compute_square_matrix_dim_from_1d_repr(ints);

    ACN_Matrix result = {.height = dim, .width = dim, .data = nullptr};

    return result;
}


Direct_ACN convert_acn_into_direct_repr(const Algebraic_Complex_Number& number) {
    if (number.is_zero()) {
        return {}; // All zeros
    }

    s64 a = number.a;
    s64 b = number.b - number.d;
    s64 c = number.c;
    s64 d = number.b + number.d;

    s64 k = number.k / 2;

    if (number.k % 2) {
        k += (number.k > 0);  // We want to avoid dividing integers, so in case the number is small, we make it scaling factor even smaller so we do not divide by 1/2

        // Multiply everything by 2/sqrt(2)
        s64 new_a = b;
        s64 new_b = 2*a;
        s64 new_c = d;
        s64 new_d = 2*c;

        a = new_a;
        b = new_b;
        c = new_c;
        d = new_d;
    }

    // Count trailing zeros to normalize K
    s64 product = a | b | c | d;
    u64 trailing_zeros = 0;
    while (!(product % 2)) { // Until the last bit is 1
        trailing_zeros += 1;
        product = product >> 1;
    }

    a = a >> trailing_zeros;
    b = b >> trailing_zeros;
    c = c >> trailing_zeros;
    d = d >> trailing_zeros;
    k = k - trailing_zeros;

    return {.a = a, .b = b, .c = c, .d = d, .k = k};
}

std::ostream& operator<<(std::ostream& os, const Direct_ACN& number) {
    const char* b_sign = number.b < 0 ? " - " : " + ";
    const char* c_sign = number.c < 0 ? " - " : " + ";
    const char* d_sign = number.d < 0 ? " - " : " + ";

    os << "("
       << number.a
       << b_sign << std::abs(number.b) << "*/sqrt(2)"
       << c_sign << std::abs(number.c) << "*i"
       << d_sign << std::abs(number.d) << "*i/sqrt(2)"
       << ") * (1/2)^(" << number.k << ")";

    return os;
}
