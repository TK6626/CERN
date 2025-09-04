#include <stdio.h>
#include <stdlib.h>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

#include "../../../../lib/cut_branch.h"
#include "../../../../lib/computations.h"
#include "../../../../lib/common_cuts.h"
#include "../../../../lib/apply_cuts.cpp"
//  !g++ src/reconstruct_all_data/rho_reconstruction/kk_reconstruction/kk_distribution.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/rho_reconstruction/kk_reconstruction/kk_distribution; ./bin/reconstruct_all_data/rho_reconstruction/kk_reconstruction/kk_distribution


/**
 * Plots and evaluates trk_dz, trk_dxy and zPV of 
 * Kaons around the dominant prodcution (reconstructed mass around 497MeV)
 * to check of presence of KK-bar pair
 **/




int main() {

	// relevent constnats
	const Int_t track_num = 4;
	const Float_t m_pi_char = m_pi_char * 1e-3; // in MeV

	// Load in aggregated uncut data
	
	ROOT::EnableImplicitMT();
	const char* file_20 = "data/glueball_mass_reconstruction/data_20_uncut.root";
	const char* file_40 = "data/glueball_mass_reconstruction/data_40_uncut.root";
	RDF df_20_df("tree", file_20);
	RDF df_40_df("tree", file_40);

	//std::cout << df_20_df.Describe() << std::endl;
	
	// applly standard cuts
	RN df_20 = CutApplier(df_20_df)
		.apply_ntrk_cut(4) 			//ensure number of tracks is 4 
    //	.apply_pt_cut(0.8)
	//	.apply_eta_cut(2.5)
    	.apply_charge_cut()
		.result();
	RN df_40 = CutApplier(df_40_df)
		.apply_ntrk_cut(4)
	//	.apply_eta_cut(0.8)
    //	.apply_pt_cut(2.5)
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






    df_20 = df_20.Define("four_momentum",
        [track_num, m_pi_char](const RVecI &trk_q, const RVecF &trk_pt, const RVecF &trk_eta, const RVecF &trk_phi) {
            std::array<Int_t, 2> pion_neg;
            std::array<Int_t, 2> pion_pos;
            Int_t pos_count = 0, neg_count = 0;
            RVecLor pion_4p(track_num);

            for (Int_t i = 0; i < track_num; ++i) {
                if (trk_q[i] == +1 && pos_count < 2) {
                    pion_pos[pos_count++] = i;
                } else if (trk_q[i] == -1 && neg_count < 2) {
                    pion_neg[neg_count++] = i;
                }
            }

            if (pos_count < 2 || neg_count < 2) {
                return RVecLor{};
            }

            for (Int_t i = 0; i < track_num; ++i) {
                TLorentzVector p_pion;
                p_pion.SetPtEtaPhiM(trk_pt[i], trk_eta[i], trk_phi[i], m_pi_char);
                pion_4p[i] = p_pion;
            }

            TLorentzVector P_00 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[0]];
            TLorentzVector P_11 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[1]];
            TLorentzVector P_01 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[1]];
            TLorentzVector P_10 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[0]];

            return RVecLor{P_00, P_01, P_10, P_11};
        },
        {"trk_q", "trk_pt", "trk_eta", "trk_phi"}
    );

    df_40 = df_40.Define("four_momentum",
        [track_num, m_pi_char](const RVecI &trk_q, const RVecF &trk_pt, const RVecF &trk_eta, const RVecF &trk_phi) {
            std::array<Int_t, 2> pion_neg;
            std::array<Int_t, 2> pion_pos;
            Int_t pos_count = 0, neg_count = 0;
            RVecLor pion_4p(track_num);

            for (Int_t i = 0; i < track_num; ++i) {
                if (trk_q[i] == +1 && pos_count < 2) {
                    pion_pos[pos_count++] = i;
                } else if (trk_q[i] == -1 && neg_count < 2) {
                    pion_neg[neg_count++] = i;
                }
            }

            if (pos_count < 2 || neg_count < 2) {
                return RVecLor{};
            }

            for (Int_t i = 0; i < track_num; ++i) {
                TLorentzVector p_pion;
                p_pion.SetPtEtaPhiM(trk_pt[i], trk_eta[i], trk_phi[i], m_pi_char);
                pion_4p[i] = p_pion;
            }

            TLorentzVector P_00 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[0]];
            TLorentzVector P_11 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[1]];
            TLorentzVector P_01 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[1]];
            TLorentzVector P_10 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[0]];

            return RVecLor{P_00, P_01, P_10, P_11};
        },
        {"trk_q", "trk_pt", "trk_eta", "trk_phi"}
    );

	TH2F* histo_20 = new TH2F("histo_20", "Topology 2", 100, 0, 1300, 100, 0, 1300);
	TH2F* histo_40 = new TH2F("histo_40", "Topology 4", 100, 0, 1300, 100, 0, 1300);



// once cuts are all done create plots for all the data and create cuts on them
// want to see if this correlates to the production of KK- kaons

	TH1F* histo_20_dz = new TH1F("histo_20_dz", "title", 100, -10, 10);
	TH1F* histo_20_dxy = new TH1F("histo_20_dxy", "title", 100, -10, 10);
	TH1F* histo_20_zpv = new TH1F("histo_20_zpv", "title", 100, -10, 10);
	
	TH1F* histo_40_dz = new TH1F("histo_40_dz", "title", 100, -10, 10);
	TH1F* histo_40_dxy = new TH1F("histo_40_dxy", "title", 100, -10, 10);
	TH1F* histo_40_zpv = new TH1F("histo_40_zpv", "title", 100, -10, 10);


	// for the kaons look for up to 1.5 sigma around the centre - heavily covered in background

	Float_t mu = 497;
	Float_t sig = 17.1;
	Float_t mass_bound = mu + sig * 1.5;

	df_20.Filter(
			[&mass_bound, &mu]
			(RVecLor p)
			{
				Bool_t pair1 = TMath::Abs(p[0].M() - mu) < mass_bound && TMath::Abs(p[1].M() - mu) < mass_bound;
				Bool_t pair2 = TMath::Abs(p[2].M() - mu) < mass_bound && TMath::Abs(p[3].M() - mu) < mass_bound;
				return pair1 || pair2;
			},
			{"four_momentum"});
		


	df_20.Foreach([&track_num, &histo_20_dxy](const RVecF& var, const RVecF& err){for (size_t i = 0; i < track_num; ++i) {if (err[i] != 0) {histo_20_dxy->Fill(var[i] / err[i]);}}}, {"trk_dxy", "trk_dxyerr"});
	df_20.Foreach([&track_num, &histo_20_dz](const RVecF& var, const RVecF& err){for (size_t i = 0; i < track_num; ++i) {if (err[i] != 0) {histo_20_dz->Fill(var[i] / err[i]);}}}, {"trk_dz", "trk_dzerr"});
	df_20.Foreach([&track_num, &histo_20_zpv](const Float_t &var){histo_20_zpv->Fill(var);}, {"zPV"});
	
	df_40.Foreach([&track_num, &histo_40_dxy](const RVecF& var, const RVecF& err){for (size_t i = 0; i < track_num; ++i) {if (err[i] != 0) {histo_40_dxy->Fill(var[i] / err[i]);}}}, {"trk_dxy", "trk_dxyerr"});
	df_40.Foreach([&track_num, &histo_40_dz](const RVecF& var, const RVecF& err){for (size_t i = 0; i < track_num; ++i) {if (err[i] != 0) {histo_40_dz->Fill(var[i] / err[i]);}}}, {"trk_dz", "trk_dzerr"});
	df_40.Foreach([&track_num, &histo_40_zpv](const Float_t &var){histo_40_zpv->Fill(var);}, {"zPV"});

	
	TCanvas* c = new TCanvas("c", "title");

	histo_20_dxy->Draw();
	c->Update();
	TFile* f1 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/K_characteristics/K_20_dxy.root", "RECREATE");
	c->Write();
	f1->Close();
	c->Clear();

	histo_20_dz->Draw();
	c->Update();
	TFile* f2 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/K_characteristics/K_20_dz.root", "RECREATE");
	c->Write();
	f2->Close();
	c->Clear();

	histo_20_zpv->Draw();
	c->Update();
	TFile* f3 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/K_characteristics/K_20_zpv.root", "RECREATE");
	c->Write();
	f3->Close();
	c->Clear();
	

	histo_40_dxy->Draw();
	c->Update();
	TFile* f4 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/K_characteristics/K_40_dxy.root", "RECREATE");
	c->Write();
	f4->Close();
	c->Clear();

	histo_40_dz->Draw();
	c->Update();
	TFile* f5 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/K_characteristics/K_40_dz.root", "RECREATE");
	c->Write();
	f5->Close();
	c->Clear();

	histo_40_zpv->Draw();
	c->Update();
	TFile* f6 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/K_characteristics/K_40_zpv.root", "RECREATE");
	c->Write();
	f5->Close();
	c->Clear();


	//cut now on data such that only pi pi pairs close to the centre are analysed
	//df_20.Foreach(
	//[&histo_20_1, &histo_20_2]
	//(RVecLor p) {
	//histo_20_1->Fill(p[0].M() *1e3, p[1].M() *1e3);
	//histo_20_2->Fill(p[2].M() *1e3, p[3].M() *1e3);
	
//	}, {"four_momentum"}
//	);

	



	//histo_20->Draw("COLZ");
	//c->Update();

	//TFile* f10 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/K_20.root", "RECREATE");
	//c->Write();
	//f10->Close();
	//c->Clear();
	
	//histo_40->Draw("COLZ");
	//c->Update();
	//TFile* f20 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/K_40.root", "RECREATE");
	//c->Write();
	//f20->Close();
	//c->Clear();

	return 0;
}
