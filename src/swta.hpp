#pragma once

#include "arith.hpp"
#include "bit_set.hpp"
#include "nfa.hpp"
#include "basics.hpp"

#include <map>
#include <vector>

typedef u64 Color;
typedef u64 State;
typedef u64 Internal_Symbol;

struct Linear_Form {
    struct Component {
        Algebraic_Complex_Number coef;
        State state;

        Component() : coef(), state(0) {}
        Component(const Algebraic_Complex_Number& number, State state) : coef(number), state(state) {}

        void swap(Component& other) noexcept {
            this->coef.swap(other.coef);
            std::swap(this->state, other.state);
        }

        void operator=(const Component& other) {
            this->coef = other.coef;
            this->state = other.state;
        }
    };

    std::vector<Component> components;  // Kept sorted by component.state

    Linear_Form() {}
    Linear_Form(const std::vector<Component>& components) : components(components) {}

    void operator=(const Linear_Form& other) {
        this->components = other.components;
    }

    size_t size() const {
        return components.size();
    }

    bool empty() const {
        return components.empty();
    }

    void normalize() {
        for (auto& component : this->components) {
            component.coef.normalize();
        }
    }
};

/**
 * Synchronized Weighted Tree Automaton
 */
struct SWTA {
    /**
     * Metadata about transition relation.
     */
    struct Metadata {
        u64 number_of_internal_symbols;
        u64 number_of_colors;
    };

    /**
     * Transition from a state along a color.
     *
     * If any of the linear forms are empty, then the state does not have a transition for the given color.
     */
    struct Transition {
        Linear_Form left;
        Linear_Form right;

        Transition() {}
        Transition(Linear_Form& l, Linear_Form& r) : left(l), right(r) {}

        void operator=(const Transition& other) {
            this->left  = other.left;
            this->right = other.right;
        }

        bool is_present() const {
            bool is_anything_missing = left.empty() || right.empty();
            return !is_anything_missing;
        }
    };

    using Transitions_Along_Color           = std::vector<Transition>;
    using Transitions_Along_Internal_Symbol = std::vector<Transitions_Along_Color>;
    using Transition_Fn                     = std::vector<Transitions_Along_Internal_Symbol>;

    /**
     * Abstraction that allows to build the transition in worklist-like algorithms.
     */
    struct Transition_Builder {
        Metadata metadata;
        std::map<State, Transitions_Along_Internal_Symbol> transitions;

        Transition_Builder(const Metadata& metadata) : metadata(metadata) {}

        void add_transition(State source_state, Internal_Symbol symbol, Color color, const Transition& transition) {
            if (!this->transitions.contains(source_state)) {
                auto& transitions_from_state = this->transitions[source_state];
                transitions_from_state.resize(this->metadata.number_of_internal_symbols);

                for (auto& transitions_along_symbol : transitions_from_state) {
                    transitions_along_symbol.resize(this->metadata.number_of_colors);
                }
            }

            auto& transitions_from_state = this->transitions[source_state];
            transitions_from_state[symbol][color] = transition;
        }

        Transition_Fn build(s64 state_cnt = -1) const {
            Transition_Fn resulting_transitions;
            state_cnt = (state_cnt < 0) ? this->transitions.rbegin()->first + 1 : state_cnt;
            resulting_transitions.resize(state_cnt);

            u64 next_state = 0;
            for (const auto& [source_state, transitions_from_state] : transitions) {
                if (source_state > next_state) {
                    for (State hole_state = next_state; hole_state < source_state; hole_state++) {
                        init_empty_transitions_for_state(resulting_transitions, hole_state);
                    }
                }

                resulting_transitions[source_state] = transitions_from_state;
                next_state = source_state + 1;
            }

            for (u64 hole_state = next_state; hole_state < static_cast<u64>(state_cnt); hole_state++) {
                init_empty_transitions_for_state(resulting_transitions, hole_state);
            }

            return resulting_transitions;
        }

        void init_empty_transitions_for_state(Transition_Fn& resulting_transitions, State state) const {
            resulting_transitions[state].resize(this->metadata.number_of_internal_symbols);
            for (Internal_Symbol sym = 0; sym < this->metadata.number_of_internal_symbols; sym++) {
                resulting_transitions[state][sym].resize(this->metadata.number_of_colors);
            }
        }

        void add_bot_state_transitions(State bot_state) {
            if (this->transitions.contains(bot_state)) return; // There already is a bot state ?!

            auto& bot_state_transitions = this->transitions[bot_state];
            bot_state_transitions.resize(this->metadata.number_of_internal_symbols);

            Linear_Form left_right ({ Linear_Form::Component(Algebraic_Complex_Number::ZERO(), bot_state) }); // Left and right subtrees have the same form
            Transition transition (left_right, left_right);

            for (Internal_Symbol symbol = 0; symbol < this->metadata.number_of_internal_symbols; symbol++) {
                bot_state_transitions[symbol].resize(this->metadata.number_of_colors);
                for (Color color = 0; color < this->metadata.number_of_colors; color++) {
                    bot_state_transitions[symbol][color] = transition;
                }
            }
        }
    };

    /**
     * Transitions of the SWTA indexed first by state, then by color, then by internal symbols.
     *
     * Details:
     *   transitions[state] are transitions from the state `state`.
     *   transitions[state][internal_symbol] are transitions from state along the symbol.
     *   transitions[state][internal_symbol][color] Is the transition ... along the color.
     *
     * @Todo: Check whether indexing is truly done in this way.
     */
    std::vector<Transitions_Along_Internal_Symbol> transitions;

    std::vector<State> initial_states;
    Bit_Set            states_with_leaf_transitions;

    SWTA(const Transition_Fn& aut_transitions, const std::vector<State>& init_states, const Bit_Set& final_state_set) :
        transitions(aut_transitions),
        initial_states(init_states),
        states_with_leaf_transitions(final_state_set) {}

    SWTA(const Transition_Fn& aut_transitions, const std::vector<State>& init_states, const std::vector<State>& final_state_set) :
        transitions(aut_transitions),
        initial_states(init_states),
        states_with_leaf_transitions(aut_transitions.size(), final_state_set) {}

    u64 number_of_states() const {
        return transitions.size();
    }

    u64 number_of_colors() const {
        if (this->transitions.empty()) { // There are no states, so the automaton is weird
            return 0;
        }
        return transitions[0][0].size();
    }

    u64 number_of_internal_symbols() const {
        return this->transitions[0].size();
    }

    Metadata get_metadata() const {
        return {.number_of_internal_symbols = number_of_internal_symbols(), .number_of_colors = number_of_colors() };
    }

    const Transition& get_transition(State source, Internal_Symbol symbol, Color color) const {
        return this->transitions[source][symbol][color];
    }
};

std::ostream& operator<<(std::ostream& os, const Linear_Form& form);
std::ostream& operator<<(std::ostream& os, const SWTA::Transition& swta_transition);
std::ostream& operator<<(std::ostream& os, const SWTA& swta);

/**
 * Weighted Tree Transducer
 */
struct WTT {
    struct Transition {
        Linear_Form ll;  // Left successor, its left  subtree inputs
        Linear_Form lr;  // Left successor, its right subtree inputs
        Linear_Form rl;  // Right successor, its left subtree inputs
        Linear_Form rr;  // Right successor, its right subtree inputs

        Transition() {}
        Transition(Linear_Form ll, Linear_Form lr, Linear_Form rl, Linear_Form rr) : ll(ll), lr(lr), rl(rl), rr(rr) {}

        Transition(Transition&& other) {
            ll = std::move(other.ll);
            lr = std::move(other.lr);
            rl = std::move(other.rl);
            rr = std::move(other.rr);
        }

        Transition(const Transition& other) : ll(other.ll), lr(other.lr), rl(other.rl), rr(other.rr) {}

        Transition& operator=(Transition&& other) {
            this->ll = other.ll;
            this->lr = other.lr;
            this->rl = other.rl;
            this->rr = other.rr;
            return *this;
        }

        Transition& operator=(const Transition& other) {
            this->ll = other.ll;
            this->lr = other.lr;
            this->rl = other.rl;
            this->rr = other.rr;
            return *this;
        }

        bool is_present() const {
            bool is_left_present  = !ll.empty() || !lr.empty();
            bool is_right_present = !rl.empty() || !rr.empty();

            return is_left_present && is_right_present;
        }
    };

    using Transitions = std::vector<std::vector<Transition>>;

    /**
     * Transitions of the WTT indexed first by state, then by internal symbol.
     *
     * Details:
     *   transitions[state] are transitions from the state `state`
     *   transitions[state][internal_state] are transitions from the state along the internal symbol.
     *
     * @Todo: Fix this - current implementation performs the indexing the wrong way (first internal sym, then state).
     */
    Transitions transitions;

    std::vector<State> initial_states;

    /**
     *  Tells which states can transform leaf.
     *  If bit_set[state] = 1, then the given state has a transition of the form state(leaf) -> leaf
     */
    Bit_Set states_with_leaf_transitions;

    WTT(const Transitions& transitions, const std::vector<State>& states_with_leaf_transitions, const std::vector<State>& initial_states) : transitions(transitions), states_with_leaf_transitions(0, nullptr), initial_states(initial_states)
    {
        Bit_Set state_set (transitions.size());
        for (State state : states_with_leaf_transitions) {
            state_set.set_bit(state, true);
        }
        this->states_with_leaf_transitions = state_set;
    };

    WTT(const Transitions& transitions, Bit_Set& leaf_states, const std::vector<State>& initial_states) : transitions(transitions), states_with_leaf_transitions(leaf_states), initial_states(initial_states) {};

    size_t number_of_states() const {
        return transitions.size();
    }

    size_t number_of_colors() const {
        return 0;
    }

    void normalize_all_transitions() {
        for (Internal_Symbol internal_symbol = 0; internal_symbol < this->transitions.size(); internal_symbol++) {
            std::vector<Transition>& transitions_for_symbol = this->transitions[internal_symbol];

            for (auto& transition : transitions_for_symbol) {
                transition.ll.normalize();
                transition.lr.normalize();
                transition.rl.normalize();
                transition.rr.normalize();
            }
        }
    }

    u64 number_of_internal_symbols() const {
        return this->transitions.size();
    }

    bool does_state_accept_trees_for_any_colored_sequence(State state) const;
    void remove_zeros_from_transitions();
};

std::ostream& operator<<(std::ostream& os, const WTT::Transition& wtt_transition);
std::ostream& operator<<(std::ostream& os, const WTT& wtt);


struct WTT_Builder {
    s64 initial_state = -1;
    Bit_Set leaf_states;
    std::map<State, std::vector<WTT::Transition>> transitions;
    u64 internal_symbol_cnt;

    WTT_Builder(SWTA::Metadata& metadata) : leaf_states(0), internal_symbol_cnt(metadata.number_of_internal_symbols) {};

    void mark_state_initial(State state) {
        this->initial_state = state;
    }

    void mark_state_final(State state) {
        this->leaf_states.grow_and_set_bit(state);
    }

    void add_transition(State source, Internal_Symbol symbol, const WTT::Transition& transition) {
        auto& transitions_from_source = this->transitions[source];
        if (transitions_from_source.empty()) transitions_from_source.resize(this->internal_symbol_cnt);

        transitions_from_source[symbol] = transition;
    }

    State find_largest_referenced_state_in_transitions() const {
        s64 max_state = -1;
        for (auto& [source, transitions]: this->transitions) {
            for (auto& transition : transitions) {
                for (auto& component : transition.ll.components) {
                    max_state = std::max(max_state, static_cast<s64>(component.state));
                }
                for (auto& component : transition.lr.components) {
                    max_state = std::max(max_state, static_cast<s64>(component.state));
                }
                for (auto& component : transition.rl.components) {
                    max_state = std::max(max_state, static_cast<s64>(component.state));
                }
                for (auto& component : transition.rr.components) {
                    max_state = std::max(max_state, static_cast<s64>(component.state));
                }
            }
        }
        assert (max_state >= 0);
        return max_state;
    }

    WTT build(s64 state_cnt = -1) {
        if (state_cnt < 0) {
            state_cnt = this->find_largest_referenced_state_in_transitions();
        }

        assert(this->initial_state >= 0);

        std::vector<std::vector<WTT::Transition>> result_transitions;
        result_transitions.resize(state_cnt);

        for (auto& [source, transitions_from_state] : this->transitions) {
            result_transitions[source] = std::move(transitions_from_state);
        }

        for (auto& transitions_from_state : result_transitions) {
            if (transitions_from_state.empty()) {
                transitions_from_state.resize(this->internal_symbol_cnt);
            }
        }

        this->leaf_states.grow(state_cnt);

        WTT transducer (result_transitions, this->leaf_states, {static_cast<u64>(this->initial_state)});
        return transducer;
    }
};



/**
 * Compute a transducer that is equivalent to first appyling the `first` transducer and then the `second`.
 */
WTT compose_wtts_sequentially(WTT& first, WTT& second);

/**
 * Compute frontier NFA for a given SWTA. If root > 0, start with the given root state.
 */
template <typename Tree_Transition_System>
NFA build_frontier_automaton(const Tree_Transition_System& tts, s64 root = -1);


enum class Subtree_Tag : u64 {
    NONE = 0, // The transition DSL definition is not finish - its missing its subtree information
    LEFT = 1,
    RIGHT = 2,
};

struct Branch_Selector {
    Internal_Symbol symbol;
    Color           color;
    Subtree_Tag     tag;

    Branch_Selector(Internal_Symbol sym, Color col, Subtree_Tag subtree_tag) : symbol(sym), color(col), tag(subtree_tag) {}

    bool operator<(const Branch_Selector& other) const {
        if (symbol < other.symbol) return true;
        if (symbol > other.symbol) return false;

        if (color < other.color) return true;
        if (color > other.color) return false;

        return static_cast<u64>(this->tag) < static_cast<u64>(other.tag);
    }
};

template <typename Transition_Tag>
struct Affine_Program {
    struct Debug_Data {
         std::map<State, std::string> state_names;
    };

    struct Transition_Symbol {
        Transition_Tag info;
        ACN_Matrix     matrix;
    };

    using Symbol_Handles = std::map<Transition_Tag, u64>;
    using Symbol_Store   = std::vector<Transition_Symbol>;

    /**
     * Given a transition symbol handle, we have a vector of possible target states - the program might
     * be non-deterministic.
     */
    using Transitions_From_State = std::vector<std::vector<State>>;
    using Transition_Fn          = std::vector<Transitions_From_State>;  // Indexed by states

    Transition_Fn  transition_fn;
    Bit_Set        final_states;
    State          initial_state;
    Symbol_Handles symbol_handles;
    Symbol_Store   symbol_store;
    Debug_Data*    debug_data = nullptr;

    Affine_Program(Transition_Fn& transitions, Bit_Set& final_states, State initial_state, Symbol_Handles& handles, Symbol_Store& store) : transition_fn(transitions), symbol_handles(handles), symbol_store(store), final_states(final_states), initial_state(initial_state) {}

    Affine_Program(const Affine_Program& other) :
        transition_fn(other.transition_fn),
        final_states(other.final_states),
        initial_state(other.initial_state),
        symbol_handles(other.symbol_handles),
        symbol_store(other.symbol_store)
    {
        if (other.debug_data != nullptr) {
            this->debug_data = new Debug_Data(*other.debug_data);
        }
    }

    u64 number_of_states() const {
        return this->transition_fn.size();
    }

    u64 number_of_symbols() const {
        return this->symbol_store.size();
    }

    ~Affine_Program() {
        if (this->debug_data != nullptr) delete this->debug_data;
        this->debug_data = nullptr;
    }
};

template <typename Transition_Tag>
struct Affine_Program_Builder {
    Affine_Program<Transition_Tag>::Symbol_Handles symbol_handles;
    Affine_Program<Transition_Tag>::Symbol_Store   symbol_store;
    Bit_Set                        final_states;
    s64                            initial_state = -1;

    std::map<State, typename Affine_Program<Transition_Tag>::Transitions_From_State> pending_transitions;

    Affine_Program_Builder (const Affine_Program<Transition_Tag>::Symbol_Handles& handles, const Affine_Program<Transition_Tag>::Symbol_Store& store) :
        symbol_handles(handles),
        symbol_store(store),
        final_states(0) {}

    void add_transition(State source_state, u64 transition_symbol, State destination_state) {
        auto& transitions_from_state = this->pending_transitions[source_state];

        if (transitions_from_state.size() < this->symbol_store.size()) {
            transitions_from_state.resize(this->symbol_store.size());
        }

        transitions_from_state[transition_symbol].push_back(destination_state);
    }

    void mark_state_final(State state) {
        this->final_states.grow_and_set_bit(state);
    }

    void mark_state_initial(State state) {
        this->initial_state = state;
    }

    void check_all_fields_filled_out() {
        if (this->initial_state < 0) throw std::runtime_error("Initial state is not set when attempting to build an affine program!");
    }

    Affine_Program<Transition_Tag> build(s64 state_cnt = -1) {
        check_all_fields_filled_out();

        // Some states might not have transitions, so we have to be careful if making decisions around the number of states based solely on discovered transitions.
        state_cnt = state_cnt > 0 ? state_cnt : this->pending_transitions.rbegin()->first + 1;

        typename Affine_Program<Transition_Tag>::Transition_Fn transitions;
        transitions.resize(state_cnt);

        State next_state = 0;
        for (auto& [state, state_transitions] : this->pending_transitions) {
            if (state > next_state) { // So we must have skipped over some states which do not have any transitions.
                for (State hole_state = next_state; hole_state < state; hole_state++) {
                    transitions[hole_state].resize(symbol_store.size());
                }
            }

            transitions[state] = std::move(state_transitions);
        }

        State last_state = this->pending_transitions.rbegin()->first;
        for (State hole_state = last_state + 1; static_cast<s64>(hole_state) < state_cnt; hole_state++) {
            transitions[hole_state].resize(symbol_store.size());
        }

        final_states.grow(state_cnt);

        return Affine_Program(transitions, final_states, initial_state, symbol_handles, symbol_store);
    }
};

struct Branch_Product_Sym {
    Color color;
    Internal_Symbol first_sym;
    Subtree_Tag     first_tag;
    Internal_Symbol second_sym;
    Subtree_Tag     second_tag;

    bool operator<(const Branch_Product_Sym& other) const {
        INSERT_LEX_LT_CODE(this->color, other.color);
        INSERT_LEX_LT_CODE(this->first_sym, other.first_sym);
        INSERT_LEX_LT_CODE(this->first_tag, other.first_tag);
        INSERT_LEX_LT_CODE(this->second_sym, other.second_sym);
        return this->second_tag < other.second_tag;
    }
};

std::ostream& operator<<(std::ostream& stream, const Branch_Product_Sym& sym);


std::ostream& operator<<(std::ostream& os, const Branch_Selector& info);

template <typename T>
void write_affine_program_into_dot(std::ostream& stream, const Affine_Program<T>& program);

bool are_two_swtas_color_equivalent(const SWTA& first, const SWTA& second);
Affine_Program<Branch_Product_Sym>
build_colored_product_of_affine_programs(const Affine_Program<Branch_Selector>& first_ap, const Affine_Program<Branch_Selector>& second_ap, const NFA& frontier);


struct Underlying_SWTA_Info {
    std::vector<State> initial_states;
    std::vector<State> final_states;
    u64 state_cnt;

    Underlying_SWTA_Info(const SWTA& swta) : initial_states(swta.initial_states), final_states(swta.states_with_leaf_transitions.into_vector()), state_cnt(swta.number_of_states()) {}
};

struct Underlying_SWTA_Info_Pair {
    Underlying_SWTA_Info first_swta_info;
    Underlying_SWTA_Info second_swta_info;
};

bool does_affine_program_reach_nonzero_final_states(const Affine_Program<Branch_Product_Sym>& program, const Underlying_SWTA_Info_Pair& swta_pair_info);

struct Color_Symbol {
    Color color;
    Internal_Symbol symbol;

    bool operator<(const Color_Symbol& other) const {
        INSERT_LEX_LT_CODE(color, other.color);
        return symbol < other.symbol;
    }
};

struct Color_Symbol_Abstraction {
    NFA abstraction;
    std::map<Color_Symbol, u64> symbol_handles;
};

Affine_Program<Branch_Selector> build_affine_program(const SWTA& swta, const Color_Symbol_Abstraction& color_sym_abstraction);

Color_Symbol_Abstraction build_color_internal_symbol_abstraction(const SWTA& swta);

SWTA apply_wtt_to_swta(const SWTA& swta, const WTT& wtt);
