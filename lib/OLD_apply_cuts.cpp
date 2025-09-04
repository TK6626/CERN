#include "computations.h"
#include "cut_branch.h"

struct CutOptions {
    bool apply_ntrk_cut = true;
    bool apply_pt_cut = true;
    bool apply_charge_cut = true;
    bool apply_zpv_cut = false;
    bool apply_dz_cut = false;
    bool apply_dxy_cut = false;
};



/**
 * Takes in a data frame and applies the desired cuts 
 *
 */
RN apply_cuts(RN df, const CutOptions& options) {
    // Cut on number of tracks
    if (options.apply_ntrk_cut) {
        Int_t track_num = 4; // 4 track_pions
        VectorConstraint<Int_t> track_constraint{
            .vector_condition = [track_num](Int_t ntrk) { return ntrk == track_num; },
            .op = Operator::all
        };
        df = cut_branch(df, "ntrk", track_constraint, true);
    }



    // Cut on transverse momentum
    if (options.apply_pt_cut) {
        float low_pt_threshold = 100000;
        VectorConstraint<float> momentum_constraint{
            .vector_condition = [low_pt_threshold](float pt) { return pt < low_pt_threshold; },
            .op = Operator::all
        };
        df = cut_branch(df, "alltrk_pt", momentum_constraint, true);
    }



    // Cut on charge
    if (options.apply_charge_cut) {
        int num_pos = 2, num_neg = 2;
        VectorConstraint<int> charge_constraint_pos{
            .vector_condition = [](int q) { return q == +1; },
            .op = Operator::equal,
            .element_count = num_pos
        };
        VectorConstraint<int> charge_constraint_neg{
            .vector_condition = [](int q) { return q == -1; },
            .op = Operator::equal,
            .element_count = num_neg
        };
        df = cut_branch(df, "trk_q", charge_constraint_pos, true);
        df = cut_branch(df, "trk_q", charge_constraint_neg, true);
    }

    

	// Cut on zPV
    if (options.apply_zpv_cut) {
        float ZPV_up = 6, ZPV_down = -6;
        VectorConstraint<float> zpv_cut{
            .vector_condition = [ZPV_down, ZPV_up](float z_pv) {
                return ZPV_down < z_pv && z_pv < ZPV_up;
            },
            .op = Operator::all
        };
        df = cut_branch(df, "zPV", zpv_cut, true);
    }

    

	// Cut on dz
    if (options.apply_dz_cut) {
        float dz_up = 2.5, dz_down = -2.5;
        VectorConstraint<float> dz_cut{
            .vector_condition = [dz_down, dz_up](float dz_val) {
                return dz_down < dz_val && dz_val < dz_up;
            },
            .op = Operator::all
        };
        df = cut_branch(df, "trk_dz", dz_cut, true);
    }

    

	// Cut on dxy currently 3 sigma
    if (options.apply_dxy_cut) {
        float dxy_up = 2.4, dxy_down = -2.4;
        VectorConstraint<float> dxy_cut{
            .vector_condition = [dxy_down, dxy_up](float dxy_val) {
                return dxy_down < dxy_val && dxy_val < dxy_up;
            },
            .op = Operator::all
        };
        df = cut_branch(df, "trk_dxy", dxy_cut, true);
    }

    return df;
}
