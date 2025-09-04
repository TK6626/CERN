#pragma once
#include "evaluate_operator.h"

template<typename T>
struct ScalarConstraint {
    Operator op;
    T value;
};

template <typename T>
struct VectorConstraint {
    std::function<bool(const T&)> vector_condition;
    Operator op;     // Now uses your existing enum
    int element_count = 0;   // Used in evaluate()
};

// template<typename T>
// struct VectorConstraint {
//     T match_value;
//     Operator op;
//     int count;
// };
