from _pytest.logging import DEFAULT_LOG_FORMAT
from swta.extended_ring_arith import Column, E_Column_Vector, Extended_Ring_Matrix, Row, add_vector_to_echelon_matrix
from dataclasses import dataclass, field


Word =  tuple[int, ...]


@dataclass
class Weighted_Automaton:
    initial_vector: Row
    transitions: dict[int, Extended_Ring_Matrix]
    final_vector: E_Column_Vector

    unsupported_colors: dict[int, set[int]] = field(default_factory=dict)

    @property
    def number_of_states(self) -> int:
        return len(self.initial_vector)


    def can_continue_with_color(self, active_states: set[int], color: int) -> bool:
        if self.unsupported_colors is None:
            return True

        states_that_lack_transition_for_color: set[int] = self.unsupported_colors.get(color, set())
        has_state_with_lacking_transition = bool(active_states.intersection(states_that_lack_transition_for_color))

        return has_state_with_lacking_transition



def construct_difference_of_weighted_automata(left: Weighted_Automaton, right: Weighted_Automaton) -> Weighted_Automaton:
    assert len(left.transitions) == len(right.transitions), '[WA Subtraction] Automata are not using the same alphabets'

    initial_vector = left.initial_vector + [-elem for elem in right.initial_vector]
    number_of_states = left.number_of_states + right.number_of_states
    new_transitions: dict[int, Extended_Ring_Matrix] = {}

    for symbol in left.transitions:
        matrix_for_symbol = Extended_Ring_Matrix.new_filled_with_zeros(width=number_of_states, height=number_of_states)

        left_matrix_for_symbol  = left.transitions[symbol]
        for left_row_idx in range(left.number_of_states):
            for left_column_idx in range(left.number_of_states):
                matrix_for_symbol.data[left_row_idx * number_of_states + left_column_idx] = left_matrix_for_symbol.at(left_row_idx,  left_column_idx)

        right_matrix_for_symbol = right.transitions[symbol]

        for right_row_idx in range(right.number_of_states):
            for right_column_idx in range(right.number_of_states):

                target_row_idx    = right_row_idx + left.number_of_states
                target_column_idx = right_column_idx + left.number_of_states

                matrix_for_symbol.data[target_row_idx * number_of_states + target_column_idx] = right_matrix_for_symbol.at(right_row_idx,  right_column_idx)

        new_transitions[symbol] = matrix_for_symbol

    raw_final_vector = left.final_vector.data + right.final_vector.data
    final_vector = E_Column_Vector(raw_final_vector)

    unsupported_colors = {}
    for color in left.transitions:
        left_unsupported  = left.unsupported_colors.get(color, set())
        right_unsupported = right.unsupported_colors.get(color, set())

        result_unsuppored = left_unsupported.union(right_state + left.number_of_states for right_state in right_unsupported)
        if result_unsuppored:
            unsupported_colors[color] = result_unsuppored

    return Weighted_Automaton(initial_vector=initial_vector, transitions=new_transitions, final_vector=final_vector)


def is_weighted_automaton_identic_to_zero(automaton: Weighted_Automaton) -> bool:
    alphabet: set[int] = set(automaton.transitions.keys())

    words_corresponding_to_basis: list[tuple[Word, Row]] = []
    basis: Extended_Ring_Matrix = Extended_Ring_Matrix.new_filled_with_zeros(width=automaton.number_of_states,
                                                                             height=automaton.number_of_states)  # @Todo(m): Is this correct

    def make_vectors_for_single_letters(letter: int) -> tuple[Word, Row]:
        word = (letter, )
        letter_matrix = automaton.transitions[letter]
        vector = letter_matrix.left_multiply_by_row(automaton.initial_vector)
        return (word, vector)

    worklist: list[tuple[Word, Row]]= [make_vectors_for_single_letters(letter) for letter in alphabet]

    while worklist:
        word, word_vector = worklist.pop(-1)

        dest_row_idx: int | None = add_vector_to_echelon_matrix(basis, word_vector)
        was_vector_added = dest_row_idx is not None

        if was_vector_added:
            words_corresponding_to_basis.insert(dest_row_idx, (word, word_vector))

            active_states = set(state for state, state_value in enumerate(word_vector) if state_value)

            for letter in alphabet:
                if not automaton.can_continue_with_color(active_states, letter):
                    continue

                extended_word = word + (letter, )

                letter_matrix = automaton.transitions[letter]
                extended_word_vector = letter_matrix.left_multiply_by_row(word_vector)

                if any(extended_word_vector):
                    worklist.append((extended_word, extended_word_vector))

        if basis.nonzero_rows == automaton.number_of_states:
            break

    resulting_dot_product = basis.multiply_by_e_column_vector(automaton.final_vector)

    is_resulting_dot_product_zero = resulting_dot_product.is_zero_vector_mod_e()
    return is_resulting_dot_product_zero
