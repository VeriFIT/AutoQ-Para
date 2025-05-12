from __future__ import annotations
from dataclasses import dataclass, field
from fractions import Fraction
import math
from typing import Any, Sequence, cast



@dataclass
class Extended_Ring_Elem:
    """
    Element of the ring Rationals[Sqrt(2)]

    Represents the number `rational_part + (1/sqrt(2))*irrational_part`

    """
    rational:   Fraction
    irrational: Fraction = field(default_factory=lambda: Fraction(0))

    def __mul__(self, other: Extended_Ring_Elem) -> Extended_Ring_Elem:
        rational   = self.rational*other.rational + Fraction(1, 2)*self.irrational*other.irrational
        irrational = self.rational*other.irrational + self.irrational*other.rational
        return Extended_Ring_Elem(rational=rational, irrational=irrational)

    def __add__(self, other: Extended_Ring_Elem) -> Extended_Ring_Elem:
        rational   = self.rational + other.rational
        irrational = self.irrational + other.irrational
        return Extended_Ring_Elem(rational=rational, irrational=irrational)

    def __sub__(self, other: Extended_Ring_Elem) -> Extended_Ring_Elem:
        rational   = self.rational - other.rational
        irrational = self.irrational - other.irrational
        return Extended_Ring_Elem(rational=rational, irrational=irrational)

    def __neg__(self) -> Extended_Ring_Elem:
        return Extended_Ring_Elem(rational=-self.rational, irrational=-self.irrational)

    def __repr__(self) -> str:
        def fmt_fraction(fraction: Fraction, take_abs=False):
            numerator = abs(fraction.numerator) if take_abs else fraction.numerator

            if fraction.denominator == 1:
                return '{0}'.format(numerator)
            return '{0}/{1}'.format(numerator, fraction.denominator)

        if not self:
            return '0'

        sign = '+' if self.irrational.numerator > 0 else '-'

        if not self.rational:
            return f'{fmt_fraction(self.irrational, take_abs=True)}*(1/√2)'
        if not self.irrational:
            return f'{fmt_fraction(self.rational)}'

        return f'{fmt_fraction(self.rational)} {sign} {fmt_fraction(self.irrational, take_abs=True)}*(1/√2)'

    def __bool__(self) -> bool:
        return bool(self.rational.numerator or self.irrational.numerator)

    def into_float(self) -> float:
        return fraction_into_float(self.rational) + fraction_into_float(self.irrational)/math.sqrt(2)

    @staticmethod
    def zero() -> Extended_Ring_Elem:
        return Extended_Ring_Elem(rational=Fraction(0), irrational=Fraction(0))

    @staticmethod
    def one() -> Extended_Ring_Elem:
        return Extended_Ring_Elem(rational=Fraction(1), irrational=Fraction(0))


Row    = list[Extended_Ring_Elem]

@dataclass
class Extended_Ring_Matrix:
    width: int
    height: int
    data: list[Extended_Ring_Elem]

    nonzero_rows: int = field(default=-1, compare=False)
    ''' Used only when dealing with row-echelon matrices. '''

    def at(self, row_idx: int, column_idx: int) -> Extended_Ring_Elem:
        idx = row_idx*self.width + column_idx
        return self.data[idx]

    def left_multiply_by_row(self, row: Row) -> Row:
        assert len(row) == self.height, f'Attempting to multiply {len(row)} row vector with {self.height} x {self.width} matrix'

        result: Row = []
        for target_column_idx in range(self.width):
            # Compute dot product with `row_elem_idx`-th column
            dot_product = Extended_Ring_Elem.zero()
            for row_elem_idx, row_elem in enumerate(row):
                column_elem = self.at(row_elem_idx, target_column_idx)
                dot_product += row_elem * column_elem

            result.append(dot_product)

        return result

    def multiply_by_matrix(self, other: Extended_Ring_Matrix) -> Extended_Ring_Matrix:
        assert self.width == other.height, f'Cannot multply {self.height} x {self.width} matrix with {other.height} x {other.width}'
        result_height = self.height
        result_width  = other.width

        result_data: list[Extended_Ring_Elem | None] = cast(list[Extended_Ring_Elem | None], [None] * result_height * result_width)
        target_elem_idx = 0
        for row_idx in range(self.height):
            for column_idx in range(other.width):
                # Compute dot product of self row and other's column
                dot_product = Extended_Ring_Elem.zero()

                for elem_idx in range(self.width):
                    self_elem  = self.at(row_idx, elem_idx)
                    other_elem = other.at(elem_idx, column_idx)
                    dot_product += self_elem * other_elem

                result_data[target_elem_idx] = dot_product

                target_elem_idx += 1

        result = Extended_Ring_Matrix(width=result_width, height=result_height, data=cast(list[Extended_Ring_Elem], result_data))
        return result

    def multiply_by_e_column_vector(self, vector: E_Column_Vector) -> E_Column_Vector:
        assert self.width == vector.dimension, f'Cannot multiply {self.height} x {self.width} matrix by a column vector of size {vector.dimension}'

        result_data: list[None | Extended_Ring_Elem] = [None] * vector.dimension
        for row_idx in range(self.height):
            dot_product: Extended_Ring_Elem | None = Extended_Ring_Elem.zero()

            for elem_idx in range(vector.dimension):
                vector_elem = vector.data[elem_idx]

                if vector_elem is None:
                    continue

                matrix_elem = self.at(row_idx, elem_idx)
                dot_product += vector_elem * matrix_elem

            result_data[row_idx] = dot_product

        return E_Column_Vector(data=result_data)


    def __str__(self) -> str:
        lines = ['[']
        for row_idx in range(self.height):
            row_str = ', '.join(repr(self.at(row_idx, col_idx)) for col_idx in range(self.width))
            lines.append(row_str)
        lines.append(']')
        lines.append('')

        return '\n'.join(lines)


    def insert_row(self, row_idx: int, row: Row, skip_shifting_subseq_rows: bool = False):
        ''' Insert a given row on position, shifting all subsequent rows by one position '''
        if not skip_shifting_subseq_rows:
            last_elem  = (self.height - 1) * self.width - 1
            first_elem = row_idx * self.width

            for elem_idx in range(last_elem, first_elem - 1, -1):
                self.data[elem_idx + self.width] = self.data[elem_idx]

        for inserted_elem_idx, inserted_elem in enumerate(row):
            self.data[row_idx * self.width + inserted_elem_idx] = inserted_elem

    def find_nonzero_elem_pos_in_row(self, row_idx: int) -> int:
        start_idx = row_idx*self.width
        for elem_idx in range(start_idx, start_idx + self.width):
            elem = self.data[elem_idx]

            if elem:
                return elem_idx - start_idx
        return self.width + 1

    def subtract_ith_row_from_row_vector(self,
                                         row_idx: int,
                                         row_coef: Extended_Ring_Elem,
                                         row_to_subtract_from: Row,
                                         row_to_subtract_from_coef: Extended_Ring_Elem) -> Row:
        new_row = []
        for elem_idx, elem_to_subtract_from in enumerate(row_to_subtract_from):
            subtractee = elem_to_subtract_from * row_to_subtract_from_coef
            subtractor = self.data[row_idx*self.width + elem_idx] * row_coef
            result_elem = subtractee - subtractor
            new_row.append(result_elem)
        return new_row

    @staticmethod
    def new_filled_with_zeros(width: int, height: int) -> Extended_Ring_Matrix:
        data = [Extended_Ring_Elem.zero() for _ in range(width*height)]
        return Extended_Ring_Matrix(width=width, height=height, data=data, nonzero_rows=0)

    @staticmethod
    def _calculcate_square_dimension_from_array(arr: Sequence[Any]) -> int:
        float_dim = math.sqrt(len(arr))
        int_dim = math.floor(float_dim)

        assert math.ceil(float_dim) == int_dim
        return int_dim

    @staticmethod
    def square_from_fractions(*fractions: Fraction) -> Extended_Ring_Matrix:
        int_dim = Extended_Ring_Matrix._calculcate_square_dimension_from_array(fractions)
        data = [Extended_Ring_Elem(rational=value, irrational=Fraction(0)) for value in fractions]
        return Extended_Ring_Matrix(width=int_dim, height=int_dim, data=data)

    @staticmethod
    def square_from_ints(*ints: int) -> Extended_Ring_Matrix:
        int_dim = Extended_Ring_Matrix._calculcate_square_dimension_from_array(ints)
        data = [Extended_Ring_Elem(rational=Fraction(value), irrational=Fraction(0)) for value in ints]
        return Extended_Ring_Matrix(width=int_dim, height=int_dim, data=data)

    @staticmethod
    def square_from_tuples(*extended_ring_elem_tuples: tuple[int, int]) -> Extended_Ring_Matrix:
        int_dim = Extended_Ring_Matrix._calculcate_square_dimension_from_array(extended_ring_elem_tuples)
        data = [Extended_Ring_Elem(rational=Fraction(rat), irrational=Fraction(irr)) for rat, irr in extended_ring_elem_tuples]
        return Extended_Ring_Matrix(width=int_dim, height=int_dim, data=data)

    @staticmethod
    def column_vector(data: list[Extended_Ring_Elem]) -> Extended_Ring_Matrix:
        return Extended_Ring_Matrix(width=1, height=len(data), data=data)

    @staticmethod
    def row_vector(data: list[Extended_Ring_Elem]) -> Extended_Ring_Matrix:
        return Extended_Ring_Matrix(width=len(data), height=1, data=data)

    @staticmethod
    def column_vector_from_ints(data: list[int]) -> Extended_Ring_Matrix:
        return Extended_Ring_Matrix(width=1, height=len(data), data=[Extended_Ring_Elem(Fraction(value)) for value in data])


def fraction_into_float(fraction: Fraction) -> float:
    return fraction.numerator / fraction.denominator


def subtract_rows(row_to_subtract_from: Row, row_to_subtract_from_coef: Extended_Ring_Elem,
                  row_to_subtract: Row, row_to_subtract_coef: Extended_Ring_Elem) -> Row:

    new_row = []
    for subtract_from_elem, to_subtract_elem in zip(row_to_subtract_from, row_to_subtract):
        subtractee = subtract_from_elem * row_to_subtract_from_coef
        subtractor = to_subtract_elem * row_to_subtract_coef
        result_elem = subtractee - subtractor
        new_row.append(result_elem)
    return new_row


def find_first_nonzero_elem_idx(row: Row) -> int:
    pos = len(row)
    for idx, elem in enumerate(row):
        if elem:
            pos = idx
            break

    if pos < len(row):
        return pos
    return pos


def reduce_matrix_into_echelon_form(matrix: Matrix):

    matrix = sorted(matrix, key=find_first_nonzero_elem_idx)

    for row_idx, row in enumerate(matrix):
        pivot_idx, _ = find_first_nonzero_elem_idx(row)

        if pivot_idx == len(row):
            break

        pivot = row[pivot_idx]

        for row_below_idx, row_below in enumerate(matrix[row_idx+1:]):
            row_below_idx += row_idx + 1

            row_below_pivot_idx, _ = find_first_nonzero_elem_idx(row_below)
            row_below_pivot = row_below[row_below_pivot_idx]

            if row_below_pivot_idx != pivot_idx:
                break  # All subsequent rows contain 0 in the column

            subtraction_result = subtract_rows(row_below, pivot, row, row_below_pivot)
            matrix[row_below_idx] = subtraction_result

    return matrix


def add_vector_to_echelon_matrix(matrix: Extended_Ring_Matrix, row: Row) -> int | None:
    '''
    Add a row to a matrix in the row-echelon form.

    The row is modified by subtracting non-empty rows of the matrix while a correct
    place for the row is identified.

    Returns:
        True, if the row has been inserted.
    '''
    row_pivot_idx = find_first_nonzero_elem_idx(row)

    row_inserted_at: int | None = None
    for matrix_row_idx in range(matrix.width):
        matrix_row_pivot_idx = matrix.find_nonzero_elem_pos_in_row(matrix_row_idx)

        if matrix_row_pivot_idx < row_pivot_idx:
            continue

        if matrix_row_pivot_idx > row_pivot_idx:
            row_inserted_at = matrix_row_idx
            matrix.insert_row(row_inserted_at, row)
            break

        # Subtract the given row and continue
        matrix_pivot = matrix.at(matrix_row_idx, matrix_row_pivot_idx)
        row_pivot    = row[row_pivot_idx]
        row = matrix.subtract_ith_row_from_row_vector(matrix_row_idx, row_pivot, row, matrix_pivot)
        row_pivot_idx = find_first_nonzero_elem_idx(row)

        if row_pivot_idx == len(row):
            return None

    if row_inserted_at is None:  # The has not been inserted yet
        is_row_empty = bool(row[-1])
        if is_row_empty:
            return None

        insert_row_at = matrix.nonzero_rows
        matrix.insert_row(insert_row_at, row, skip_shifting_subseq_rows=True)
        matrix.nonzero_rows += 1
        return insert_row_at

    return row_inserted_at


@dataclass
class E_Column_Vector:
    data: list[Extended_Ring_Elem | None]

    @property
    def dimension(self) -> int:
        return len(self.data)

    @staticmethod
    def from_ints(ints: list[int | None]):
        def into_ring_elem(maybe_value: int | None) -> Extended_Ring_Elem | None:
            if maybe_value is None:
                return maybe_value
            return Extended_Ring_Elem(Fraction(maybe_value))

        data = [into_ring_elem(it) for it in ints]
        return E_Column_Vector(data=data)

    def is_zero_vector_mod_e(self) -> bool:
        ''' Returns True if all of the non-e elements are zero '''
        contains_non_zero = any(self.data)
        return not contains_non_zero


Column = Extended_Ring_Matrix
Vector = Extended_Ring_Matrix
