from enum import IntEnum
from fractions import Fraction
from swta.extended_ring_arith import E_Column_Vector, Extended_Ring_Elem, Extended_Ring_Matrix
from swta.weighted_automata import Weighted_Automaton, construct_difference_of_weighted_automata, is_weighted_automaton_identic_to_zero

import pytest


class Colors(IntEnum):
    wL = 0
    wR = 1
    aL = 2
    aR = 3


@pytest.fixture
def bv_post_condition() -> Weighted_Automaton:
    wL_raw = [
        0, 1, 0, 0,
        0, 1, 0, 0,
        1, 0, 0, 0,
        0, 0, 0, 0,
    ]

    wR_raw = [
        0, 0, 1, 0,
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 0, 0,
    ]

    aL_raw = [
        0, 1, 0, 0,
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 0, 0,
    ]

    aR_raw = [
        0, 0, 0, 1,
        1, 0, 0, 0,
        0, 0, 0, 1,
        0, 0, 0, 0,
    ]

    final_vector = E_Column_Vector.from_ints([None, 0, None, 1])
    raw_initial_vector = [1, 0, 0, 0]
    initial_vector = [Extended_Ring_Elem(Fraction(it)) for it in raw_initial_vector]

    automaton = Weighted_Automaton(
        initial_vector=initial_vector,
        transitions={
            Colors.wL: Extended_Ring_Matrix.square_from_ints(*wL_raw),
            Colors.wR: Extended_Ring_Matrix.square_from_ints(*wR_raw),
            Colors.aL: Extended_Ring_Matrix.square_from_ints(*aL_raw),
            Colors.aR: Extended_Ring_Matrix.square_from_ints(*aR_raw),
        },
        final_vector=final_vector
    )

    return automaton


@pytest.fixture
def bv_constructed_simplified() -> Weighted_Automaton:
    wL_raw = [
        0, 0, 0,
        1, 0, 0,
        0, 0, 0,
    ]

    wR_raw = [
        0, 1, 0,
        0, 0, 0,
        0, 0, 0,
    ]

    aL_raw = [
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
    ]

    aR_raw = [
        0, 0, 1,
        0, 0, 1,
        0, 0, 0,
    ]

    final_vector = E_Column_Vector.from_ints([None, None, 1])
    raw_initial_vector = [1, 0, 0]
    initial_vector = [Extended_Ring_Elem(Fraction(it)) for it in raw_initial_vector]

    automaton = Weighted_Automaton(
        initial_vector=initial_vector,
        transitions={
            Colors.wL: Extended_Ring_Matrix.square_from_ints(*wL_raw),
            Colors.wR: Extended_Ring_Matrix.square_from_ints(*wR_raw),
            Colors.aL: Extended_Ring_Matrix.square_from_ints(*aL_raw),
            Colors.aR: Extended_Ring_Matrix.square_from_ints(*aR_raw),
        },
        final_vector=final_vector
    )

    return automaton


class BV_Complex_State_Names(IntEnum):
    ''' Documents how are states in the complex automaton numbered. '''
    GAMMA   = 0
    DELTA   = 1
    EPSILON = 2
    MU      = 3
    OMEGA   = 4


@pytest.fixture
def bv_constructed_complex() -> Weighted_Automaton:
    wL_raw = [
        (0, 0), (0, 1), (0, 1), (0, 0), (0, 0),
        (1, 0), (0, 0), (0, 0), (0, 0), (0, 0),
        (0, 0), (0, 0), (0, 0), (1, 0), (0, 0),
        (0, 0), (0, 1), (0, 1), (0, 0), (0, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),
    ]

    wR_raw = [
        (0, 0), (0, 1), (0,-1), (0, 0), (0, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),
        (0, 0), (0,-1), (0, 1), (0, 0), (0, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),
    ]

    aL_raw = [
        (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),
    ]

    aR_raw = [
        (0, 0), (0, 0), (0, 0), (0, 0), ( 1, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), ( 1, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), (-1, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), (-1, 0),
        (0, 0), (0, 0), (0, 0), (0, 0), ( 0, 0),
    ]


    final_vector = E_Column_Vector.from_ints([None, None, None, None, 1])
    raw_initial_vector = [1, 0, 0, 0, 0]
    initial_vector = [Extended_Ring_Elem(Fraction(it)) for it in raw_initial_vector]


    leaf_states = {
        BV_Complex_State_Names.OMEGA.value,
    }

    unsuppored_colors = {
        Colors.wL.value: leaf_states,
        Colors.wR.value: leaf_states,
        Colors.aL.value: leaf_states,
        Colors.aR.value: leaf_states,
    }

    automaton = Weighted_Automaton(
        initial_vector=initial_vector,
        transitions={
            Colors.wL: Extended_Ring_Matrix.square_from_tuples(*wL_raw),
            Colors.wR: Extended_Ring_Matrix.square_from_tuples(*wR_raw),
            Colors.aL: Extended_Ring_Matrix.square_from_tuples(*aL_raw),
            Colors.aR: Extended_Ring_Matrix.square_from_tuples(*aR_raw),
        },
        final_vector=final_vector,
        unsupported_colors=unsuppored_colors
    )

    return automaton


def test_equivalence_on_bv_automata(bv_post_condition: Weighted_Automaton, bv_constructed_simplified: Weighted_Automaton):
    difference = construct_difference_of_weighted_automata(bv_post_condition, bv_constructed_simplified)
    assert difference

    res = is_weighted_automaton_identic_to_zero(difference)
    assert res


def test_equivalence_on_bv_automata_with_complex(bv_post_condition: Weighted_Automaton, bv_constructed_complex: Weighted_Automaton):
    difference = construct_difference_of_weighted_automata(bv_post_condition, bv_constructed_complex)
    assert difference

    res = is_weighted_automaton_identic_to_zero(difference)
    assert res

