#pragma once

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <iostream>
#include "custom_definitions.h"
#include "constraints.h"

// Template declarations
template<typename DF, typename T>
RN cut_branch(
    DF df,
    const std::string& branch_name,
    const ScalarConstraint<T>& constraint,
    bool verbose = true);

template<typename DF, typename T>
RN cut_branch(
    DF df,
    const std::string& branch_name,
    const VectorConstraint<T>& constraint,
    bool verbose = true);

// Include template implementations
#include "cut_branch.tpp"
