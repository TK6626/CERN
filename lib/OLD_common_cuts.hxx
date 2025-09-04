#include "common_cuts.h"
#include <string>
#include "ROOT/RDataFrame.hxx"
using ROOT::RDF::RNode;
/**
 * @brief Filter a DataFrame on a scalar branch using a comparison operation.
 *
 * This function filters entries in a ROOT DataFrame based on a scalar branch value
 * and a comparison operation (e.g., "==", ">=", "<=", ">", "<").
 *
 * @tparam DF         The type of the input DataFrame (ROOT::RDataFrame or a filtered node).
 * @tparam T          The type of the branch and threshold (e.g., int, float).
 * @param data_frame  The input DataFrame.
 * @param threshold   The value to compare against. MUST BE THE SAME DATA TYPE AS THE COLUMN!
 * @param branch_name The name of the scalar branch to cut on.
 * @param operation   The comparison operation as a string: "==", ">=", "<=", ">", "<".
 * @param verbose     If true, print the number of events passing the cut.
 * @return            A new filtered DataFrame node.
 *
 * @example
 *     auto df_cut = cut_branch_single(df, 1.0, "my_branch", ">=", true);
 */
template<typename DF, typename T>
RNode  cut_branch_single(
    DF data_frame,
    T threshold,
    const char* branch_name,
    const std::string& operation,
    Bool_t verbose)
{
    auto track_cut = [threshold, operation](T x) {
        if (operation == "==") return x == threshold;
        if (operation == ">=") return x >= threshold;
        if (operation == "<=") return x <= threshold;
        if (operation == ">")  return x > threshold;
        if (operation == "<")  return x < threshold;
        return false;
    };
    auto df_cut = data_frame.Filter(track_cut, {branch_name});
    if (verbose) {
        std::cout <<"events such that " << branch_name << " has value " << operation << " " << threshold << " = "<<  *df_cut.Count() << std::endl;
    }
    return df_cut;
}
template RNode cut_branch_single<ROOT::RDataFrame, int>(
    ROOT::RDataFrame, int, const char*, const std::string&, bool);


/**
 * @brief Filter a DataFrame on the sum of an RVec branch using a comparison operation.
 *
 * This function filters entries in a ROOT DataFrame by summing the values of an RVec branch
 * and applying a comparison operation (e.g., "==", ">=", "<=", ">", "<") to the sum.
 *
 * @tparam DF         The type of the input DataFrame (ROOT::RDataFrame or a filtered node).
 * @tparam T          The numeric type of the RVec elements and threshold (e.g., float, int).
 * @param data_frame  The input DataFrame.
 * @param threshold   The value to compare the sum against. MUST BE THE SAME DATA TYPE AS THE COLUMN!
 * @param branch_name The name of the RVec branch to sum and cut on.
 * @param operation   The comparison operation as a string: "==", ">=", "<=", ">", "<".
 * @param verbose     If true, print the number of events passing the cut.
 * @return            A new filtered DataFrame node.
 *
 * @example
 *     auto df_cut = cut_branch_sum<float>(df, 10.0, "my_RVec_branch", ">", true);
 */
template<typename DF, typename T>
RNode cut_branch_sum(
    DF data_frame,
    T threshold,
    const char* branch_name,
    const std::string& operation,
    Bool_t verbose)
{
    auto sum_cut = [threshold, operation](const ROOT::VecOps::RVec<T>& vec) {
        T sum = Sum(vec);
        if (operation == "==") return sum == threshold;
        if (operation == ">=") return sum >= threshold;
        if (operation == "<=") return sum <= threshold;
        if (operation == ">")  return sum > threshold;
        if (operation == "<")  return sum < threshold;
        return false;
    };
    auto df_cut = data_frame.Filter(sum_cut, {branch_name});
    if (verbose) {
        std::cout <<"number of events such that sum(" << branch_name << ") " << operation << " " << threshold << " = "<<  *df_cut.Count() << std::endl;
    }
    return df_cut;
}

