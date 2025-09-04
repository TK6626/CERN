#pragma once
#include "custom_definitions.h"
#include "constraints.h"

template<typename DF, typename T>
RN cut_branch(
    DF df,
    const std::string& branch_name,
    const ScalarConstraint<T>& constraint,
    bool verbose)
{
    auto lambda = [constraint](const T& value) {
        return evaluate(value, constraint.op, constraint.value);
    };

    auto filtered_df = df.Filter(lambda, {branch_name});
    if (verbose)
        std::cout << "Scalar cut on '" << branch_name
                  << "' passed " << *filtered_df.Count() << " events.\n";

    return filtered_df;
}


// example - 
// auto constraint = VectorConstraint<float>{
//     .vector_condition = [](float x) { return x > 25; }, // is the element greater than 25 
//     .op = Operator::greater_equal,
//     .element_count = 3               /at least 3 entries satisfy the top
// };

template<typename DF, typename T>
RN cut_branch(
    DF df,
    const std::string& branch_name,
    const VectorConstraint<T>& constraint,
    bool verbose)
{
    auto lambda = [constraint](const ROOT::RVec<T>& vec) {
        int count = std::count_if(vec.begin(), vec.end(), constraint.vector_condition);
        if (constraint.op == Operator::all)
            return count == static_cast<int>(vec.size());
        return evaluate(count, constraint.op, constraint.element_count);
    };

    auto filtered_df = df.Filter(lambda, {branch_name});
    if (verbose)
        std::cout << "Vector cut on '" << branch_name
                  << "' passed " << *filtered_df.Count() << " events.\n";

    return filtered_df;
}
