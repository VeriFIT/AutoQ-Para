from swta.extended_ring_arith import (
    Extended_Ring_Elem,
    Extended_Ring_Matrix,
    reduce_matrix_into_echelon_form,
    add_vector_to_echelon_matrix
)
from fractions import Fraction

import pytest


@pytest.mark.skip()
def test_add_vector_to_echelon_matrix():
    a11 = Extended_Ring_Elem(Fraction(1), Fraction(0))
    a12 = Extended_Ring_Elem(Fraction(2), Fraction(0))
    a13 = Extended_Ring_Elem(Fraction(3), Fraction(0))

    a21 = Extended_Ring_Elem(Fraction(2), Fraction(0))
    a22 = Extended_Ring_Elem(Fraction(4), Fraction(0))
    a23 = Extended_Ring_Elem(Fraction(6), Fraction(0))

    test_matrix = [
        [a11, a12, a13],
    ]
    test_row = [a21, a22, a23]
    add_vector_to_echelon_matrix(test_matrix, test_row, add_empty=False)
    assert test_matrix == [[a11, a12, a13]]


def test_matrix_left_multiply_by_row():
    matrix = Extended_Ring_Matrix.square_from_ints(
        0, 2, 0,
        0, 1, 2,
        0, 0, 0,
    )
    raw_vector = [1, 0, 0]
    vector = [Extended_Ring_Elem(Fraction(value)) for value in raw_vector]

    result = matrix.left_multiply_by_row(vector)
    print(result)

    assert result == [Extended_Ring_Elem(Fraction(0)), Extended_Ring_Elem(Fraction(2)), Extended_Ring_Elem(Fraction(0))]


def test_add_row_to_echelon_form_matrix():
    matrix = Extended_Ring_Matrix.square_from_ints(
        1, 2, 0,
        0, 1, 2,
        0, 0, 0,
    )
    raw_vector = [0, 1, 2]
    vector = [Extended_Ring_Elem(Fraction(value)) for value in raw_vector]

    was_added = add_vector_to_echelon_matrix(matrix, vector)
    assert not was_added
