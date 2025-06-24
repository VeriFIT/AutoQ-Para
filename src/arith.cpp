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

std::ostream& operator<<(std::ostream& os, const Algebraic_Complex_Number& number) {
    os << "("
       << mpz_get_si(number.a) << ", "
       << mpz_get_si(number.b) << ", "
       << mpz_get_si(number.c) << ", "
       << mpz_get_si(number.d) << ", "
       << mpz_get_si(number.k) << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Fixed_Precision_ACN& number) {
    os << "("
       << number.a << ", "
       << number.b << ", "
       << number.c << ", "
       << number.d << ", "
       << number.k << ")";
    return os;
}

s64 add_row_to_row_echelon_matrix_no_copy(ACN_Matrix& matrix, ACN_Matrix& row) {
    assert(matrix.width == row.width);

    if (row.contains_only_zeros()) return -1;

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

        row.subtract_from_ith_row(
            0,             // Which row of the `row` matrix to subtract from
            matrix_pivot,  // What coefficient should multiply the row before multiplication
            matrix,        // Matrix containing rows that we can subtract from
            row_idx,       // Which of the matrix rows are we subtracting
            row_pivot);    // What coeffcient should multiply the row

        inserted_row_pivot_idx = row.find_nonzero_elem_in_row(0);

        if (inserted_row_pivot_idx >= row.width) {
            return -1; // Not inserted, Gaussian elimination reduced the row to 0
        }
    }

    assert(false);  // Should be unreachable - the row is indepentent and we insert it, or our matrix has full rank and hence we do not insert it
}

s64 add_row_to_row_echelon_matrix(ACN_Matrix& matrix, const ACN_Matrix& row) {
    ACN_Matrix row_copy (1, row.width); // Make a local copy, since we will be modifying it
    for (auto row_elem_idx = 0; row_elem_idx < row.width; row_elem_idx++) {
        auto& elem_value = row.at(0, row_elem_idx);
        row_copy.set(0, row_elem_idx, elem_value);
    }
    return add_row_to_row_echelon_matrix_no_copy(matrix, row_copy);
}

u64 compute_square_matrix_dim_from_1d_repr(const std::vector<s64>& repr) {
    u64 dimension = static_cast<s64>(std::sqrt(repr.size()));

    assert(dimension*dimension == repr.size());

    return dimension;
}

ACN_Matrix square_acn_matrix_from_ints(const std::vector<s64>& ints) {
    u64 dim = compute_square_matrix_dim_from_1d_repr(ints);

    ACN_Matrix result(dim, dim);
    Algebraic_Complex_Number* matrix_slots = new Algebraic_Complex_Number[dim*dim];

    for (u64 elem_idx = 0; elem_idx < ints.size(); elem_idx++) {
        s64 elem_int_value = ints[elem_idx];
        mpz_set_si(matrix_slots[elem_idx].a, elem_int_value);
    }

    result.data = matrix_slots;
    return result;
}

ACN_Matrix row_from_ints(const std::vector<s64>& row_data) {
    Algebraic_Complex_Number* matrix_slots = new Algebraic_Complex_Number[row_data.size()];

    for (u64 elem_idx = 0; elem_idx < row_data.size(); elem_idx++) {
        s64 elem_int_value = row_data[elem_idx];
        mpz_set_si(matrix_slots[elem_idx].a, elem_int_value);
    }

    ACN_Matrix result(1, row_data.size(), matrix_slots);
    return result;
}

ACN_Matrix column_from_ints(const std::vector<s64>& column_data) {
    Algebraic_Complex_Number* matrix_slots = new Algebraic_Complex_Number[column_data.size()];

    for (u64 elem_idx = 0; elem_idx < column_data.size(); elem_idx++) {
        s64 elem_int_value = column_data[elem_idx];
        mpz_set_si(matrix_slots[elem_idx].a, elem_int_value);
    }

    ACN_Matrix result(column_data.size(), 1, matrix_slots);
    return result;
}


Direct_ACN convert_acn_into_direct_repr(const Algebraic_Complex_Number& number) {
    if (number.is_zero()) {
        return {}; // All zeros
    }

    s64 a = mpz_get_si(number.a);
    s64 b = mpz_get_si(number.b) - mpz_get_si(number.d);
    s64 c = mpz_get_si(number.c);
    s64 d = mpz_get_si(number.b) + mpz_get_si(number.d);

    s64 input_k = mpz_get_si(number.k);
    s64 k = input_k / 2;

    if (input_k % 2) {
        k += (input_k > 0);  // We want to avoid dividing integers, so in case the number is small, we make it scaling factor even smaller so we do not divide by 1/2

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


std::ostream& operator<<(std::ostream& os, const ACN_Matrix& matrix) {
    os << "[";
    if (matrix.height > 1) os << "\n";
    for (u64 row_idx = 0; row_idx < matrix.height; row_idx++) {
        for (u64 col_idx = 0; col_idx < matrix.width; col_idx++) {
            os << matrix.at(row_idx, col_idx) << ", ";
        }

        if (matrix.height > 1) os << "\n";
    }
    os << "]";

    return os;
}

