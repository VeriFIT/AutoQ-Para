#pragma once

#include "arith.hpp"
#include "bit_set.hpp"

#include <map>
#include <set>
#include <string>
#include <sstream>

namespace std {
    template <typename T>
    ostream& operator<<(ostream& os, std::vector<T> vec) {
        os << "[";
        for (u64 idx = 0; idx < vec.size(); idx++) {
            auto& elem = vec[idx];
            os << elem;
            if (idx + 1 != vec.size()) {
                os << ", ";
            }
        }
        os << "]";
        return os;
    }
}

struct NFA {
    using State = u64;
    using Transitions_From_State = std::vector<std::vector<State>>;
    using Transition_Fn = std::vector<Transitions_From_State>;

    struct Debug_Data {
        std::map<State, std::string> state_names;
    };

    std::vector<State> initial_states;
    Bit_Set            final_states;
    Transition_Fn      transitions;
    Debug_Data*        debug_data = nullptr;

    NFA() : initial_states({}), final_states(0), transitions({}) {}
    NFA(const std::vector<State>& init_states, const std::vector<State>& fin_states, const std::vector<Transitions_From_State>& aut_transitions) :
        initial_states(init_states),
        final_states(aut_transitions.size()),
        transitions(aut_transitions)
    {
        for (auto state: fin_states) {
            final_states.set_bit(state, true);
        }
    }

    NFA(const NFA& other) :
        initial_states(other.initial_states),
        final_states(other.final_states),
        transitions(other.transitions)
    {
        if (other.debug_data != nullptr) {
            this->debug_data = new Debug_Data(*other.debug_data);
        }
    };

    NFA(const std::vector<State>& init_states, const Bit_Set& fin_states, const std::vector<Transitions_From_State>& aut_transitions) :
        initial_states(init_states),
        final_states(fin_states),
        transitions(aut_transitions) {}

    ~NFA() {
        if (this->debug_data != nullptr) delete this->debug_data;
        this->debug_data = nullptr;
    }

    u64 number_of_states() const {
        return this->transitions.size();
    }

    u64 alphabet_size() const {
        return this->transitions[0].size();
    }

    NFA determinize() const;
    bool complete();
    bool is_every_state_accepting() const;
    void write_dot(std::ostream& stream) const;
};


struct NFA_Builder {
    std::map<NFA::State, std::vector<std::set<NFA::State>>> transitions;
    std::vector<NFA::State> initial_states;
    Bit_Set final_states;
    u64 alphabet_size;

    NFA_Builder(u64 alphabet_size) : alphabet_size(alphabet_size), final_states(0) {};

    void mark_state_initial(NFA::State state) {
        this->initial_states.push_back(state);
    }

    void mark_state_final(NFA::State state) {
        this->final_states.grow_and_set_bit(state);
    }

    void add_transition(NFA::State source_state, u64 symbol, NFA::State destination);
    NFA build(s64 state_cnt = -1);
};


struct Macrostate {
    using State = u64;

    std::vector<State> state_names; // Accelerates computation of Post
    Bit_Set            state_set;   // To test what states are present in the macrostate
    State              handle;      // Assigned after initialization

    Macrostate(u64 state_cnt) : state_names({}), state_set(state_cnt) {}
    Macrostate(u64 size, const std::vector<State>& content) : state_names(content), state_set(size, content) {}

    Macrostate(Macrostate&& other) : state_names(other.state_names), state_set(other.state_set), handle(other.handle) {}

    Macrostate(const Macrostate& other) : state_names(other.state_names), state_set(other.state_set), handle(other.handle) {}

    bool operator<(const Macrostate& other) const {
        return state_set < other.state_set;
    }

    Macrostate& operator=(Macrostate&& other) {
        this->handle      = other.handle;
        this->state_set   = std::move(other.state_set);
        this->state_names = std::move(other.state_names);

        return *this;
    }

    Macrostate& operator=(const Macrostate& other) {
        this->handle      = other.handle;
        this->state_set   = other.state_set;
        this->state_names = other.state_names;

        return *this;
    }

    bool empty() const {
        return state_names.empty();
    }

    void populate_state_names_from_set() {
        this->state_names.clear();
        for (State state = 0; state < this->state_set.size; state++) {
            if (this->state_set.get_bit_value(state)) {
                this->state_names.push_back(state);
            }
        }
    }

    std::string to_string() const {
        std::stringstream ss;
        ss << "Macrostate " << state_names << " (handle=" << handle << ")";
        return ss.str();
    }
};

template <typename T>
struct Worklist_Construction_Context {
    std::map<T, u64>      handles;
    std::vector<const T*> worklist;

    const T* extract() {
        auto elem = worklist.back();
        worklist.pop_back();
        return elem;
    }

    const T* mark_discovery(T& discovery) {
        discovery.handle = this->handles.size();
        auto [insert_pos, was_inserted] = this->handles.emplace(discovery, discovery.handle);

        const T* result_ptr = &insert_pos->first;

        if (was_inserted) {
            this->worklist.push_back(result_ptr);
        }

        return result_ptr;
    }

    bool has_more_to_explore() const {
        return !this->worklist.empty();
    }
};

struct State_Pair {
    u64 first, second;
    u64 handle = -1;

    /**
     * Lexigraphical ordering.
     */
    bool operator<(const State_Pair& other) const {
        if (this->first < other.first)
            return true;
        else if (this->first > other.first)
            return false;
        return this->second < other.second;
    }

    bool operator==(const State_Pair& other) const {
        return first == other.first && second == other.second;
    }
};

std::ostream& operator<<(std::ostream& os, const State_Pair& state);
std::ostream& operator<<(std::ostream& os, const Macrostate& macrostate);

bool are_two_complete_dfas_equivalent(const NFA& first_nfa, NFA& other_nfa);
