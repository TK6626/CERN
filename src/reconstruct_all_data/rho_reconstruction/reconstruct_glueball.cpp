#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

#include "../../../lib/cut_branch.h"
#include "../../../lib/computations.h"
#include "../../../lib/apply_cuts.cpp"
#include "../../../lib/plotting_params.h"


//  !g++ src/reconstruct_all_data/rho_reconstruction/reconstruct_glueball.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/rho_reconstruction/reconstruct_glueball; ./bin/reconstruct_all_data/rho_reconstruction/reconstruct_glueball



/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */
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
    	.apply_pt_cut(00000.8f)
		.apply_eta_cut(2.5)
    	.apply_charge_cut()
		.result();
	RN df_40 = CutApplier(df_40_df)
		.apply_ntrk_cut(4)
		.apply_eta_cut(2.5)
    	.apply_pt_cut(00000.8f)
    	.apply_charge_cut()
		.result();

df_20 = 	df_20.Define("trk_dz_snr", [] (Int_t ntrk, RVecF z, RVecF err)
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
	


df_40 = 	df_40.Define("trk_dz_snr", [] (Int_t ntrk, RVecF z, RVecF err)
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


	// apply seperate cuts on each of the z, xy and zPV to see if any have significant impact
	
	Float_t trk_z_bound = 2.38;
VectorConstraint<Float_t> trk_z_constraint { 
    .vector_condition = [trk_z_bound](const Float_t& val) {
        return std::abs(val) < trk_z_bound;
    },                                          
    .op = Operator::all,
    .element_count = 0
};

Float_t trk_xy_bound = 2.8;
VectorConstraint<Float_t> trk_xy_constraint { 
    .vector_condition = [trk_xy_bound](const Float_t& val) {
        return std::abs(val) < trk_xy_bound;
    },                                          
    .op = Operator::all,
    .element_count = 0
};

Float_t Z_bound = 13.1;
VectorConstraint<Float_t> Z_constraint {  
    .vector_condition = [Z_bound](const Float_t& val) {
        return std::abs(val) < Z_bound;
    },                     
    .op = Operator::all,
    .element_count = 0
};


df_20 = cut_branch(df_20, "trk_dz_snr", trk_z_constraint, true);
df_40 = cut_branch(df_40, "trk_dz_snr", trk_z_constraint, true);

df_20 = cut_branch(df_20, "trk_dxy_snr", trk_xy_constraint, true);
df_40 = cut_branch(df_40, "trk_dxy_snr", trk_xy_constraint, true);

df_20 = cut_branch(df_20, "trk_dz_snr", Z_constraint, true);
df_40 = cut_branch(df_40, "trk_dz_snr", Z_constraint, true);






    df_20 = df_20.Define("four_momentum",
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
    );

    df_40 = df_40.Define("four_momentum",
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
    );




	Float_t mass_bound = 30; // MeV/c^2


	// Plot the invariant mass of the reconstructed particles
	TH1F* histo_20 = new TH1F("histo_20", "Invariant Mass;M (MeV/c^{2});Counts", 70, 1600, 2600);

	df_20.Foreach(
    	[&histo_20, mass_bound](const ROOT::VecOps::RVec<TLorentzVector>& p4) {
        	if (TMath::Abs(p4[0].M()*1e3 - m_p) < mass_bound && TMath::Abs(p4[1].M()*1e3 - m_p) < mass_bound) {
            	histo_20->Fill((p4[0] + p4[1]).M() * 1e3);
       	 	}
       	 	if (TMath::Abs(p4[2].M()*1e3 - m_p) < mass_bound && TMath::Abs(p4[3].M()*1e3 - m_p) < mass_bound) {
            	histo_20->Fill((p4[2] + p4[3]).M() * 1e3);
        	}
    	},
    	{"four_momentum"}
	);

	TH1F* histo_40 = new TH1F("histo_40", "Invariant Mass;M (MeV/c^{2});Counts", 70, 1600, 2600);

	df_40.Foreach(
    	[&histo_40, mass_bound](const ROOT::VecOps::RVec<TLorentzVector>& p4) {
        	if (TMath::Abs(p4[0].M()*1e3 - m_p) < mass_bound && TMath::Abs(p4[1].M()*1e3 - m_p) < mass_bound) {
            	histo_40->Fill((p4[0] + p4[1]).M() * 1e3);
        	}
        	if (TMath::Abs(p4[2].M()*1e3 - m_p) < mass_bound && TMath::Abs(p4[3].M()*1e3 - m_p) < mass_bound) {
            	histo_40->Fill((p4[2] + p4[3]).M() * 1e3);
        	}
    	},
    	{"four_momentum"}
	);
    
	TCanvas* can1 = new TCanvas("can1");
	histo_20->Draw("E1");
	can1->Update();

    // Save histograms as ROOT files
    TFile* f20 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/glueball_20.root", "RECREATE");
    can1->Write();
    f20->Close();
	
	can1->Clear();
	histo_40->Draw("E1");
	can1->Update();

    TFile* f40 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/glueball_40.root", "RECREATE");
    can1->Write();
    f40->Close();






	return 0;
}
