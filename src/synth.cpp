#include <functional>

#include "swta.hpp"

WTT makeWTTApplyingOpQubits(
    const int qubitCount,
    std::function<void(std::vector<std::vector<WTT::Transition>>, int, int)> transitionMaker,
    const std::vector<int>& mask)
{
    size_t stateCnt = qubitCount;

    std::vector<std::vector<WTT::Transition>> wttTransitions;

    State state = 0;

    int interestingQubitsRead = 0;
    int nextInterestingQubit  = mask[interestingQubitsRead];

    for (int qubitIdx = 0; qubitIdx < qubitCount; qubitIdx++) {
        Linear_Form anihilator({Linear_Form::Component(Algebraic_Complex_Number::ZERO(), state + 1)});
        Linear_Form unit({Linear_Form::Component(Algebraic_Complex_Number::ONE(), state + 1)});

        bool isQubitInteresting = (qubitIdx == nextInterestingQubit);
        bool areAllInterestingQubitsSeen = interestingQubitsRead >= mask.size();

        if (areAllInterestingQubitsSeen || !isQubitInteresting) {
            // Just add an identity state (preserve both subtrees as they are)
            WTT::Transition transition (unit, anihilator, anihilator, unit);
            wttTransitions.push_back({transition});
            continue;
        }

        // We see the next interesting qubit, let the caller make the transition since it is interesting to him
        transitionMaker(wttTransitions, state, interestingQubitsRead);

        interestingQubitsRead += 1;

        if (interestingQubitsRead < mask.size()) {
            nextInterestingQubit = mask[interestingQubitsRead];
        }
        state += 1;
    }
    Bit_Set leaves (stateCnt, {state});
    std::vector<State> initialStates ({0});
    WTT ret (wttTransitions, leaves, initialStates);

    return ret;
}

/**
 * Make a projection with qubit values from mask being equal to basis.
 */
WTT makeProjectionWTT(const int qubitCount, const std::vector<int>& partialBasis, const std::vector<int>& mask) {
    assert(partialBasis.size() == mask.size());

    auto transitionMaker = [&partialBasis](std::vector<std::vector<WTT::Transition>> transitions, int currentState, int interestingQubitIdx) {
        Linear_Form anihilator({Linear_Form::Component(Algebraic_Complex_Number::ZERO(), currentState + 1)});
        Linear_Form unit({Linear_Form::Component(Algebraic_Complex_Number::ONE(), currentState + 1)});

        int desiredQubitValue = partialBasis[interestingQubitIdx];
        if (desiredQubitValue) {
            // We want to see 1 in the projection, so anihilate the left subtree
            WTT::Transition transition (anihilator, anihilator, anihilator, unit);
            transitions.push_back({transition});
        } else {
            // We want to see 0 in the projection, so anihilate the right subtree
            WTT::Transition transition (unit, anihilator, anihilator, anihilator);
            transitions.push_back({transition});
        }
    };

    WTT result = makeWTTApplyingOpQubits(qubitCount, transitionMaker, mask);
    return result;
}

WTT makeRotationWTT(
    const int qubitCount,
    const std::vector<int>& sourceBasis,
    const std::vector<int>& targetBasis,
    const std::vector<int>& mask)
{
    assert(sourceBasis.size() == targetBasis.size());
    assert(sourceBasis.size() == mask.size());

    auto transitionMaker = [&sourceBasis, &targetBasis](std::vector<std::vector<WTT::Transition>> transitions, int currentState, int interestingQubitIdx) {
        Linear_Form emptyForm;
        Linear_Form unit({Linear_Form::Component(Algebraic_Complex_Number::ONE(), currentState + 1)});

        int sourceBit = sourceBasis.at(interestingQubitIdx);
        int targetBit = targetBasis.at(interestingQubitIdx);

        if (sourceBit == targetBit) { // No swapping
            WTT::Transition idTransition (unit, emptyForm, emptyForm, unit);
            transitions.push_back({idTransition});
        } else {
            WTT::Transition swapTransition (emptyForm, unit, unit, emptyForm);
            transitions.push_back({swapTransition});
        }
    };

    WTT result = makeWTTApplyingOpQubits(qubitCount, transitionMaker, mask);
    return result;
}
