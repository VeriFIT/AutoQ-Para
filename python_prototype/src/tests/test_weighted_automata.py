from fractions import Fraction
from swta.extended_ring_arith import E_Column_Vector, Extended_Ring_Elem, Extended_Ring_Matrix
from swta.weighted_automata import (
    construct_difference_of_weighted_automata,
    Weighted_Automaton,
    is_weighted_automaton_identic_to_zero,
)

import pytest


@pytest.fixture
def left_simple_automaton():
    left = Weighted_Automaton(
        initial_vector=[Extended_Ring_Elem.one(), Extended_Ring_Elem.zero(), Extended_Ring_Elem.zero()],
        transitions={
            0: Extended_Ring_Matrix.square_from_fractions(
                Fraction(0), Fraction(2), Fraction(0),
                Fraction(0), Fraction(1), Fraction(1, 2),
                Fraction(0), Fraction(0), Fraction(0),
            )
        },
        final_vector=E_Column_Vector.from_ints([0, 0, 1])
    )
    return left


@pytest.fixture
def right_simple_automaton():
    right = Weighted_Automaton(
        initial_vector=[Extended_Ring_Elem.one(), Extended_Ring_Elem.zero(), Extended_Ring_Elem.zero()],
        transitions={
            0: Extended_Ring_Matrix.square_from_fractions(
                Fraction(0), Fraction(3), Fraction(0),
                Fraction(0), Fraction(1), Fraction(1, 3),
                Fraction(0), Fraction(0), Fraction(0),
            )
        },
        final_vector=E_Column_Vector.from_ints([0, 0, 1])
    )
    return right


def test_difference_construction(left_simple_automaton: Weighted_Automaton, right_simple_automaton: Weighted_Automaton):
    '''
    Both automata are over single-letter alphabet.

    Both automata accept runs of the form  `a * 1 * .. * 1 * 1/a`
    '''
    left, right = left_simple_automaton, right_simple_automaton

    result = construct_difference_of_weighted_automata(left, right)

    assert result.initial_vector == [
        Extended_Ring_Elem(Fraction(1)), Extended_Ring_Elem(Fraction(0)), Extended_Ring_Elem(Fraction(0)),
        Extended_Ring_Elem(Fraction(-1)), Extended_Ring_Elem(Fraction(0)), Extended_Ring_Elem(Fraction(0))
    ]

    raw_expected_transition_matrix = [
        Fraction(0), Fraction(2), Fraction(0),    Fraction(0), Fraction(0), Fraction(0),
        Fraction(0), Fraction(1), Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0),
        Fraction(0), Fraction(0), Fraction(0),    Fraction(0), Fraction(0), Fraction(0),
        Fraction(0), Fraction(0), Fraction(0),    Fraction(0), Fraction(3), Fraction(0),
        Fraction(0), Fraction(0), Fraction(0),    Fraction(0), Fraction(1), Fraction(1, 3),
        Fraction(0), Fraction(0), Fraction(0),    Fraction(0), Fraction(0), Fraction(0),
    ]
    expected_transition_matrix = Extended_Ring_Matrix.square_from_fractions(*raw_expected_transition_matrix)

    assert len(result.transitions) == 1
    assert result.transitions[0] == expected_transition_matrix


def test_zero_checks():
    raw_transition_matrix = [
        Fraction(0), Fraction(2), Fraction(0),    Fraction(0), Fraction(0), Fraction(0),
        Fraction(0), Fraction(1), Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0),
        Fraction(0), Fraction(0), Fraction(0),    Fraction(0), Fraction(0), Fraction(0),
        Fraction(0), Fraction(0), Fraction(0),    Fraction(0), Fraction(3), Fraction(0),
        Fraction(0), Fraction(0), Fraction(0),    Fraction(0), Fraction(1), Fraction(1, 3),
        Fraction(0), Fraction(0), Fraction(0),    Fraction(0), Fraction(0), Fraction(0),
    ]
    transition_matrix = Extended_Ring_Matrix.square_from_fractions(*raw_transition_matrix)

    raw_initial_vector = [1, 0, 0, -1, 0, 0]
    initial_vector = [Extended_Ring_Elem(Fraction(value)) for value in raw_initial_vector]

    zero_automaton = Weighted_Automaton(
        initial_vector=initial_vector,
        transitions={
            0: transition_matrix
        },
        final_vector=E_Column_Vector.from_ints([0, 0, 1, 0, 0, 1])
    )

    assert is_weighted_automaton_identic_to_zero(zero_automaton)
