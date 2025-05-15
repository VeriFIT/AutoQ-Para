#include <cassert>
#include <cmath>
#include <cstdlib>
#include <new>
#include <ostream>
#include <iostream>

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

struct ACN_Matrix {
    u64 height, width;
    Algebraic_Complex_Number* data = nullptr;

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

        return ACN_Matrix {
            .height = result_height,
            .width  = result_width,
            .data   = result_data
        };
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

        return this->width + 1;
    }

    void insert_row_at(const ACN_Matrix& row, u64 row_idx, bool skip_shifting_subsequent_rows = false) {
        assert(row.width == this->width);

        if (!skip_shifting_subsequent_rows) {
            u64 last_elem_idx  = (this->height - 1) * this->width - 1;
            u64 first_elem_idx = row_idx * this->width;

            for (u64 elem_idx = last_elem_idx; elem_idx >= first_elem_idx; elem_idx--) {
                this->data[elem_idx + this->width] = this->data[elem_idx];
            }
        }

        for (u64 elem_idx = 0; elem_idx < row.width; elem_idx++) {
            this->data[row_idx*this->width] = row.at(0, elem_idx);
        }
    }

    void subtract_from_ith_row(u64 row_idx, Algebraic_Complex_Number& row_coef, ACN_Matrix& rows_to_subtract, u64 row_to_subtract_idx, Algebraic_Complex_Number& row_to_subtract_coef) const {
        for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
            auto subtractee_weighted = this->at(row_idx, elem_idx) * row_coef;
            auto subtractor_weighted = rows_to_subtract.at(row_idx, elem_idx) * row_to_subtract_coef;
            auto result_elem = subtractee_weighted - subtractor_weighted;

            this->data[row_idx*this->width + elem_idx] = result_elem;
        }
    }

    ~ACN_Matrix() {
        if (this->data != nullptr) delete[] this->data;
    }
};


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
       << b_sign << number.b << "*(1/sqrt(2))"
       << c_sign << number.c << "*i"
       << d_sign << number.d << "*(1/sqrt(2))*i"
       << ") * (1/2)^(" << number.k << ")";

    return os;
}
