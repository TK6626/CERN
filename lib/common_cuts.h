
#pragma once
// common_cuts.h

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <iostream>
#include "custom_definitions.h"
#include "constraints.h"

// Template declarations

template<typename DF>
RN DefineVectorBranch(DF df,
                         const std::string& value_branch,
                         const std::string& error_branch,
                         const std::string& count_branch,
                         const std::string& new_branch_name = "");

template<typename DF>
RN DefineScalarBranch(DF df,
                          const std::string& value_branch,
                          const std::string& error_branch,
                          const std::string& new_branch_name = "");

// Include template implementations
#include "common_cuts.tpp"

