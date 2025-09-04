#pragma once
#include "custom_definitions.h"
#include <string>

template<typename DF>
RN DefineVectorBranch(DF df,
                         const std::string& value_branch,
                         const std::string& error_branch,
                         const std::string& count_branch,
                         const std::string& new_branch_name) {
    std::string out_branch = new_branch_name.empty()
        ? value_branch + "_SNR"
        : new_branch_name;

    return df.Define(out_branch, [=](const Int_t n_trk,
                                     const ROOT::RVec<float>& val,
                                     const ROOT::RVec<float>& err) {
        ROOT::RVec<float> result(n_trk);
        for (size_t i = 0; i < static_cast<size_t>(n_trk); ++i) {
            result[i] = (err[i] != 0.0f) ? val[i] / err[i] : 0.0f;
        }
        return result;
    }, {count_branch, value_branch, error_branch});
}

template<typename DF>
RN DefineScalarBranch(DF df,
                          const std::string& value_branch,
                          const std::string& error_branch,
                          const std::string& new_branch_name) {
    std::string out_branch = new_branch_name.empty()
        ? value_branch + "_SNR"
        : new_branch_name;

    return df.Define(out_branch, [=](float val, float err) {
        return (err != 0.0f) ? val / err : 0.0f;
    }, {value_branch, error_branch});
}
