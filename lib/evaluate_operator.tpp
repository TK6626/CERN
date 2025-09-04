
// EVALUATE_OPERATOR.CPP // 
#pragma once
#include "TMath.h"

template<typename T>
bool evaluate(const T& actual, Operator op, const T& nominal) {
    switch (op) {
        case Operator::equal:           return actual == nominal;
        case Operator::not_equal:       return actual != nominal;
        case Operator::greater:         return actual > nominal;
        case Operator::less:            return actual < nominal;
        case Operator::greater_equal:   return actual >= nominal;
        case Operator::less_equal:      return actual <= nominal;

		case Operator::abs_less:          return TMath::Abs(actual) <  TMath::Abs(nominal);
        case Operator::abs_equal:         return TMath::Abs(actual) == TMath::Abs(nominal);
        case Operator::abs_greater:       return TMath::Abs(actual) >  TMath::Abs(nominal);
        case Operator::abs_less_equal:    return TMath::Abs(actual) <= TMath::Abs(nominal);
        case Operator::abs_greater_equal: return TMath::Abs(actual) >= TMath::Abs(nominal);
		default:                      return false;
    }
}
