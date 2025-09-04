#pragma once

// Enum for comparison
enum class Operator {
    equal,
    not_equal,
    greater,
    less,
    greater_equal,
    less_equal,
    all,
	abs_less,
	abs_greater,
	abs_equal,
	abs_less_equal,
	abs_greater_equal
};

// Template declaration
template<typename T>
bool evaluate(const T& actual, Operator op, const T& expected);

// Include implementation
#include "evaluate_operator.tpp"
