#include "custom_definitions.h"
#include "computations.h"
#include "cut_branch.h"



/**
 * @class CutApplier
 * @brief Applies a series of physics-related selection cuts to an input dataframe (RN).
 *
 * The CutApplier class provides a clean, modular Int_terface to apply individual
 * filtering operations (cuts) on a ROOT-style dataframe used in physics analysis.
 * Each cut corresponds to a method, and users may call only the cuts they require.
 * The class supports method chaining for fluent-style usage.
 *
 * Example usage:
 * @code
 * RN final_df = CutApplier(my_df)
 *     .apply_ntrk_cut()          // Require 4 tracks
 *     .apply_pt_cut(120000.0f)   // pT threshold of 120 GeV
 *     .apply_charge_cut()        // 2 positive and 2 negative tracks
 *     .result();                 // Get the final filtered dataframe
 * @endcode
 *
 * To apply cuts step-by-step:
 * @code
 * CutApplier cutter(my_df);
 * cutter.apply_ntrk_cut();
 * cutter.apply_dxy_cut();
 * RN output = cutter.result();
 * @endcode
 *
 * @note All methods return a reference to the CutApplier object itself,
 *       enabling method chaining.
 */
class CutApplier {
public:
    explicit CutApplier(RN input_df) : df(std::move(input_df)) {}

	// use chaining method



	// apply cut on track number
    CutApplier& apply_ntrk_cut(Int_t track_num = 4) {
        VectorConstraint<Int_t> constraint{
            .vector_condition = [track_num](Int_t ntrk) { return ntrk == track_num; },
            .op = Operator::all
        };
        df = cut_branch(df, "ntrk", constraint, true);
        return *this;
    }
	
	// apply cut on preceisely 2+ve and 2-ve tracks
    CutApplier& apply_charge_cut(Int_t num_pos = 2, Int_t num_neg = 2) {
        VectorConstraint<Int_t> pos_constraint{
            .vector_condition = [](Int_t q) { return q == +1; }, // condition each individual entry in tthe RVec that must be satisfied
            .op = Operator::equal, // the condition that the number of etries within a given RVec must be EQUAL to the NUM_POS
            .element_count = num_pos
        };
        VectorConstraint<Int_t> neg_constraint{
            .vector_condition = [](Int_t q) { return q == -1; },
            .op = Operator::equal,
            .element_count = num_neg
        };
        df = cut_branch(df, "trk_q", pos_constraint, true);
        df = cut_branch(df, "trk_q", neg_constraint, true);
        return *this;
    }
	
	// apply cut on transverse (momentum)
    CutApplier& apply_pt_cut(Float_t low = 0, Float_t high = 5) {
        VectorConstraint<Float_t> constraint{
            .vector_condition = [low, high](Float_t pt) { return low < pt && pt < high; },
            .op = Operator::all
        };
        df = cut_branch(df, "trk_pt", constraint, true);
        return *this;
    }

	// apply cut on total(momentum)
    CutApplier& apply_p_cut(Float_t low = 0, Float_t high = 5) {
        VectorConstraint<Float_t> constraint{
            .vector_condition = [low, high](Float_t p) { return low < p && p < high; },
            .op = Operator::all
        };
        df = cut_branch(df, "trk_p", constraint, true);
        return *this;
    }

	// apply cut on pseudo-rapidity (eta)
	 CutApplier& apply_eta_cut(Float_t low = 0, Float_t high = 2.5) {
        VectorConstraint<Float_t> constraint{
            .vector_condition = [low, high](Float_t eta) { return low < eta && eta < high; },
            .op = Operator::all
        };
        df = cut_branch(df, "trk_eta", constraint, true);
        return *this;
    }	
	
	
	// note for these characteristic variable plots they 
	// are structured as mu ± n * sigma

	// Apply cut on z momentum (3 sigma)
    CutApplier& apply_zPV_cut(Float_t low = (0 - 3 * 4.4), Float_t high = (0 +  3 * 4.4)) {
        VectorConstraint<Float_t> constraint{
            .vector_condition = [low, high](Float_t zpv) { return low < zpv && zpv < high; },
            .op = Operator::all
        };
        df = cut_branch(df, "zPV", constraint, true);
        return *this;
    }
	
	// Apply cut on masses outside of a given mas bound using RLorVectors
    CutApplier& apply_mass_bound_cut(
			const char* p4_branch, Float_t m_0, Float_t bound = 0.1){
        VectorConstraint<TLorentzVector> constraint{
            .vector_condition = [m_0, bound](TLorentzVector p) { return ( (p.M()-m_0)* (p.M()-m_0) < bound * bound); },
            .op = Operator::all,
			.element_count = 0
        };
        df = cut_branch(df, p4_branch, constraint, true);
        return *this;
    }


	// Apply cut on reconstructed z Int_teraction poInt_t (3 sigma)
    CutApplier& apply_dz_cut(Float_t low = (-0.5), Float_t high = (0.5)){
        VectorConstraint<Float_t> constraint{
            .vector_condition = [low, high](Float_t dz) { return low < dz && dz < high; },
            .op = Operator::all,
			.element_count = 0
        };
        df = cut_branch(df, "trk_dz", constraint, true);
        return *this;
    }

	// Apply cut on reconstructed z_snr Int_teraction poInt_t (3 sigma)
    CutApplier& apply_dz_snr_cut(Float_t low = (-4.2e-3 - 3 * 0.78), Float_t high = (-4.2e-3 + 3 * 0.78)){
        VectorConstraint<Float_t> constraint{
            .vector_condition = [low, high](Float_t dz) { return low < dz && dz < high; },
            .op = Operator::all,
			.element_count = 0
        };
        df = cut_branch(df, "trk_dz_snr", constraint, true);
        return *this;
    }


	// Apply cut on reconstructed xy Int_teraction poInt_t (3 sigma)
    CutApplier& apply_dxy_cut(Float_t low = -0.3, Float_t high = 0.3) {
        VectorConstraint<Float_t> constraint{
            .vector_condition = [low, high](Float_t dxy) { return low < dxy && dxy < high; },
            .op = Operator::all,
			.element_count = 0
        };
        df = cut_branch(df, "trk_dxy", constraint, true);
        return *this;
    }

	// apply cut on fractional momentum loss

	// Apply cut on reconstructed xy_snr Int_teraction poInt_t (3 sigma)
    CutApplier& apply_dxy_snr_cut(Float_t low = (3.9e-3 - 3 * 0.94), Float_t high = (3.9e-3 + 3 * 0.94)) {
        VectorConstraint<Float_t> constraint{
            .vector_condition = [low, high](Float_t dxy) { return low < dxy && dxy < high; },
            .op = Operator::all,
			.element_count = 0
        };
        df = cut_branch(df, "trk_dxy", constraint, true);
        return *this;
    }

	// apply cut on fractional momentum loss





	// create result to allow 
    RN result() const {
        return df;
    }



private:
    RN df;
};

