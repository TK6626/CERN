#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/plotting_params.h"
#include "../../../../lib/admin_utils.h"

// !./bashing/run_file.sh reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/generic_data



/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */

void AddTrkSNR(RN &df);
VectorConstraint<Float_t> CreateChiConstraint(Float_t bound, Operator op = Operator::all, size_t element_count = 0);

int main() {

	// relevent constnats
	const Float_t pion_mass = 0.135;
	const Int_t track_num = 4;



	// Load in aggregated uncut data
	
	ROOT::EnableImplicitMT();
	const char* file_20 = "data/glueball_mass_reconstruction/data_20_uncut.root";
	const char* file_40 = "data/glueball_mass_reconstruction/data_40_uncut.root";
	RDF df_20_df("tree", file_20);
	RDF df_40_df("tree", file_40);

	std::cout << df_20_df.Describe() << std::endl;
	
	// applly standard cuts
	RN df_20 = CutApplier(df_20_df)
		.apply_ntrk_cut(4) 			//ensure number of tracks is 4 
    	.apply_pt_cut(0.8f) // cut on transverse mometum
		.apply_eta_cut(2.5) // cut on pseudorapidity
    	.apply_charge_cut()	// cut ensureing precisesly 2 +ve, 2 -ve charges
		.apply_dxy_cut(-2.8, 2.8)
		.apply_dz_cut(-2.38, 2.38)
		.apply_zPV_cut(-13.1, 13.1)
		.result();
	RN df_40 = CutApplier(df_40_df)
		.apply_ntrk_cut(4)
		.apply_eta_cut(0.8)
    	.apply_pt_cut(00000.8f)
    	.apply_charge_cut()
		.apply_dxy_cut(-2.8, 2.8)
		.apply_dz_cut(-2.38, 2.38)
		.apply_zPV_cut(-13.1, 13.1)
		.result();

    df_20.Define("four_momentum",
        [track_num, pion_mass](const RVecI &trk_q, const RVecF &trk_pt, const RVecF &trk_eta, const RVecF &trk_phi) {
            std::array<Int_t, 2> pion_neg;
            std::array<Int_t, 2> pion_pos;
            Int_t pos_count = 0, neg_count = 0;
            ROOT::VecOps::RVec<TLorentzVector> pion_4p(track_num);

            for (Int_t i = 0; i < track_num; ++i) {
                if (trk_q[i] == +1 && pos_count < 2) {
                    pion_pos[pos_count++] = i;
                } else if (trk_q[i] == -1 && neg_count < 2) {
                    pion_neg[neg_count++] = i;
                }
            }

            if (pos_count < 2 || neg_count < 2) {
                return ROOT::VecOps::RVec<TLorentzVector>{};
            }

            for (Int_t i = 0; i < track_num; ++i) {
                TLorentzVector p_pion;
                p_pion.SetPtEtaPhiM(trk_pt[i], trk_eta[i], trk_phi[i], pion_mass);
                pion_4p[i] = p_pion;
            }

            TLorentzVector P_00 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[0]];
            TLorentzVector P_11 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[1]];
            TLorentzVector P_01 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[1]];
            TLorentzVector P_10 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[0]];

            return ROOT::VecOps::RVec<TLorentzVector>{P_00, P_01, P_10, P_11};
        },
        {"trk_q", "trk_pt", "trk_eta", "trk_phi"}
    ).Snapshot("tree", "data/glueball_mass_reconstruction/data_20_SNR.root");

    df_40.Define("four_momentum",
        [track_num, pion_mass](const RVecI &trk_q, const RVecF &trk_pt, const RVecF &trk_eta, const RVecF &trk_phi) {
            std::array<Int_t, 2> pion_neg;
            std::array<Int_t, 2> pion_pos;
            Int_t pos_count = 0, neg_count = 0;
            ROOT::VecOps::RVec<TLorentzVector> pion_4p(track_num);

            for (Int_t i = 0; i < track_num; ++i) {
                if (trk_q[i] == +1 && pos_count < 2) {
                    pion_pos[pos_count++] = i;
                } else if (trk_q[i] == -1 && neg_count < 2) {
                    pion_neg[neg_count++] = i;
                }
            }

            if (pos_count < 2 || neg_count < 2) {
                return ROOT::VecOps::RVec<TLorentzVector>{};
            }

            for (Int_t i = 0; i < track_num; ++i) {
                TLorentzVector p_pion;
                p_pion.SetPtEtaPhiM(trk_pt[i], trk_eta[i], trk_phi[i], pion_mass);
                pion_4p[i] = p_pion;
            }

            TLorentzVector P_00 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[0]];
            TLorentzVector P_11 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[1]];
            TLorentzVector P_01 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[1]];
            TLorentzVector P_10 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[0]];

            return ROOT::VecOps::RVec<TLorentzVector>{P_00, P_01, P_10, P_11};
        },
        {"trk_q", "trk_pt", "trk_eta", "trk_phi"}
    ).Snapshot("tree", "data/glueball_mass_reconstruction/data_40_SNR.root");
	

	return 0;
}


VectorConstraint<Float_t> CreateChiConstraint(Float_t bound, Operator op, size_t element_count) {
    VectorConstraint<Float_t> constraint;
    constraint.vector_condition = [bound](const Float_t& val) {
        return std::abs(val) < bound;
    };
    constraint.op = op;
    constraint.element_count = element_count;
    return constraint;
}


void AddTrkSNR(RN &df) {
    df = df
        .Define("trk_dz_snr", [] (Int_t ntrk, RVecF z, RVecF err)
        {
            RVecF result(ntrk);
            for (Int_t i = 0; i < ntrk; ++i) {
                result[i] = (err[i] != 0.0f) ? z[i] / err[i] : 0.0f;
            }
            return result;
        }, {"ntrk", "trk_dz", "trk_dzerr"})
        .Define("trk_dxy_snr", [] (Int_t ntrk, RVecF z, RVecF err)
        {
            RVecF result(ntrk);
            for (Int_t i = 0; i < ntrk; ++i) {
                result[i] = (err[i] != 0.0f) ? z[i] / err[i] : 0.0f;
            }
            return result;
        }, {"ntrk", "trk_dxy", "trk_dxyerr"});
}
