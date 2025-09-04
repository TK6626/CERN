#include <stdio.h>
#include <stdlib.h>

#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TLorentzVector.h"

#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/apply_cuts.cpp"
#include "../../../../lib/plotting_params.h"


//  !g++ src/reconstruct_all_data/rho_reconstruction/elastic_background/proton_momentum_y.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/rho_reconstruction/elastic_background/proton_momentum_y; ./bin/reconstruct_all_data/rho_reconstruction/elastic_background/proton_momentum_y



/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */
int main() {
	
	SetPlotStyle();
	// relevent constnats
	const Int_t track_num = 4;
	const Float_t pion_mass = m_pi_char * 1e-3; // GeV


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
    //	.apply_pt_cut(0.8f)
	//	.apply_eta_cut(2.5)
    	.apply_charge_cut()
		.result();
	RN df_40 = CutApplier(df_40_df)
		.apply_ntrk_cut(4)
	//	.apply_eta_cut(2.5)
    //	.apply_pt_cut(0.8f)
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


//df_20 = cut_branch(df_20, "trk_dz_snr", trk_z_constraint, true);
//df_40 = cut_branch(df_40, "trk_dz_snr", trk_z_constraint, true);

//df_20 = cut_branch(df_20, "trk_dxy_snr", trk_xy_constraint, true);
//df_40 = cut_branch(df_40, "trk_dxy_snr", trk_xy_constraint, true);

//df_20 = cut_branch(df_20, "trk_dz_snr", Z_constraint, true);
//df_40 = cut_branch(df_40, "trk_dz_snr", Z_constraint, true);




	//Decide wheter to create the 4ps 
	Bool_t create_4p_20 = false;
	Bool_t create_4p_40 = false;
	

	if(create_4p_20){
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
	};

	if(create_4p_40){
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
	};

// Look at the distribution of the proton momentum in the a and b bins
// suggestive that if the difference in momentum of the protons is 'too' in a given direction
// then no interaction occured...
//
// Simmilarly should the difference in momentum be 'too large' then the interaction may not be
// elastic but inelastic



// look at distribution of y momentum 
// 1D difference histogram and 2D visual distribution
	TH1F* y_20_1d = new TH1F("y_20_1d", "Difference in p_{y};p_{y} (MeV/c^{2});Events", 500, -2, 2);
	TH2F* y_20_2d = new TH2F("y_20_2d", "Comparison of p_{y};Pot_a p_{y} (MeV/c^{2});Pot_b p_{y} (MeV/c^{2})", 500, -1, 1, 500, -1, 1);
	TCanvas* c_20 = new TCanvas("c_20", "Title");
	df_20.Foreach(
			[&y_20_1d, &y_20_2d] (Double_t py_a, Double_t py_b) {
				y_20_1d->Fill(py_a - py_b);
				y_20_2d->Fill(py_a, py_b);	
			},{"pr_py_a", "pr_py_b"});
	y_20_1d->GetYaxis()->SetTitle(Form("Events [%.2f MeV/c^{2}]", y_20_1d->GetBinWidth(1)));
	y_20_2d->GetYaxis()->SetTitle(Form("Pot_a p_{y} MeV/c^{2} [%.2g]", y_20_2d->GetYaxis()->GetBinWidth(1)));
	y_20_2d->GetXaxis()->SetTitle(Form("Pot_b p_{y} MeV/c^{2} [%.2f]", y_20_2d->GetXaxis()->GetBinWidth(1)));

	// fit the double gtlp-stat -saussian
	TF1* y_20_1d_fit = new TF1("y_20_1d_fit", fit_double_gaussian, -2, 2, 8);
	y_20_1d_fit->SetParameters(-0.26, 0.1, 5600, 0, 0.2, 0.1, 4800, 0);
	y_20_1d_fit->SetParNames("μ_1","σ_1", "A_1", "C_1");
	y_20_1d_fit->SetParLimits(0, -0.3, -0.25);
	y_20_1d_fit->SetParLimits(1, 0.15, 0.3);
	y_20_1d_fit->SetParLimits(4, 0.15, 0.25);
	y_20_1d_fit->SetParLimits(5, 0.15, 0.3);
	y_20_1d_fit->SetLineColor(kRed);
	y_20_1d_fit->SetLineWidth(2);
	y_20_1d_fit->SetLineStyle(2);
	TFitResultPtr r_y_20_1d = y_20_1d->Fit("y_20_1d_fit", "S0");

	y_20_1d->SetMarkerStyle(20);
	y_20_1d->SetMarkerColor(red_tran);
	y_20_1d->SetLineColor(red_tran);


	y_20_1d->Draw("E1");
	y_20_1d_fit->Draw("SAME");
	c_20->Update();
	TFile* f_py_20_1d = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/proton_momentum_y_20_difference.root", "RECREATE");
	c_20->Write();
	f_py_20_1d->Close();
	c_20->Clear();

	y_20_2d->Draw("COLZ");
	c_20->Update();
	TFile* f_py_20_2d = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/proton_momentum_y_20_comparison.root", "RECREATE");
	c_20->Write();
	f_py_20_1d->Close();
	c_20->Clear();

// 40 Topology
// 1D difference histogram and 2D visual distribution
	TH1F* y_40_1d = new TH1F("y_40_1d", "Difference in p_{y};p_{y} (MeV/c^{2});Events", 500, -0.75, 0.75);
	TH2F* y_40_2d = new TH2F("y_40_2d", "Comparison of p_{y};Pot_a p_{y} (MeV/c^{2});Pot_b p_{y} (MeV/c^{2})", 500, -1, 1 !g++ src/reconstruct_all_data/rho_reconstruction/elastic_background/proton_momentum_y.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/rho_reconstruction/elastic_background/    proton_momentum_y; ./bin/reconstruct_all_data/rho_reconstruction/elastic_background/proton_momentum_y
, 500, -1, 1);
	TCanvas* c_40 = new TCanvas("c_40", "Title");
	df_40.Foreach(
			[&y_40_1d, &y_40_2d] (Double_t py_a, Double_t py_b) {
				y_40_1d->Fill(py_a - py_b);
				y_40_2d->Fill(py_a, py_b);	
			},{"pr_py_a", "pr_py_b"});
	y_40_1d->GetYaxis()->SetTitle(Form("Events [%.2f MeV/c^{2}]", y_40_1d->GetBinWidth(1)));
	y_40_2d->GetYaxis()->SetTitle(Form("Pot_a p_{y} MeV/c^{2} [%.2g]", y_40_2d->GetYaxis()->GetBinWidth(1)));
	y_40_2d->GetXaxis()->SetTitle(Form("Pot_b p_{y} MeV/c^{2} [%.2f]", y_40_2d->GetXaxis()->GetBinWidth(1)));

	// fit the double gtlp-stat -saussian
	TF1* y_40_1d_fit = new TF1("y_40_1d_fit", fit_double_gaussian, -2, 2, 8);
	y_40_1d_fit->SetParameters(-0.26, 0.1, 5600, 0, 0.2, 0.1, 4800, 0);
	y_40_1d_fit->SetParNames("μ_1","σ_1", "A_1", "C_1");
	y_40_1d_fit->SetParLimits(0, -0.3, -0.25);
	y_40_1d_fit->SetParLimits(1, 0.15, 0.3);


	y_40_1d_fit->SetParLimits(4, 0.15, 0.25);
	y_40_1d_fit->SetParLimits(5, 0.15, 0.3);
	y_40_1d_fit->SetLineColor(kRed);
	y_40_1d_fit->SetLineWidth(2);
	y_40_1d_fit->SetLineStyle(2);
	TFitResultPtr r_y_40_1d = y_40_1d->Fit("y_40_1d_fit", "S0");

	y_40_1d->SetMarkerStyle(40);
	y_40_1d->SetMarkerColor(red_tran);
	y_40_1d->SetLineColor(red_tran);


	y_40_1d->Draw("E1");
	y_40_1d_fit->Draw("SAME");
	c_40->Update();
	TFile* f_py_40_1d = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/proton_momentum_y_40_difference.root", "RECREATE");
	c_40->Write();
	f_py_40_1d->Close();
	c_40->Clear();

	y_40_2d->Draw("COLZ");
	c_40->Update();
	TFile* f_py_40_2d = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/elastic_background/proton_momentum_y_40_comparison.root", "RECREATE");
	c_40->Write();
	f_py_40_1d->Close();
	c_40->Clear();


	
	return 0;
}
