#include <array>
#include <stdio.h>
#include <cmath>
#include <vector>

#include "TMath.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TLegend.h"

// include headers containing useful functions and definitions
#include "../lib/computations.h"
#include "../lib/common_cuts.h"
using RD = ROOT::RDataFrameh;



/*
 * 	Construct the mass of the 2-pair mesons from the 4-meson tracks 
 * 
 *	cut1:
 *	cut2:
 *	cut3:
 *
 *
 * */

int main(){
	
	// constants
	const Float_t pion_mass = 0.139; // GeV/c^2

	const Float_t x_min = 260;
	const Float_t x_max = 1000;
	const Int_t x_bins = 300;
	const Float_t y_min = x_min;
	const Float_t y_max = x_max;
	const Int_t y_bins = x_bins;

	TH2* graph = new  TH2F("graph", "Reconstructed  #rho mesons;#rho_{1} (MeV / c^{2}) ;#rho_{1} (MeV / c^{2})",
								x_bins, x_min, x_max, y_bins,y_min, y_max);
	TCanvas* canvas = new TCanvas("canvas");


	// load in all 40 type events ,get column names, enable multithreading
    const char* file_name = "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM40.root";
	ROOT::EnableImplicitMT();

	std::vector<std::string> files = 
	{"/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM40.root",
	"/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM41.root",
	"/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM42.root",
	"/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM43.root"};
	
	RD df("tree", files);
	std::cout << df.Describe() << std::endl;
	
	//get total number of origional entries
	std::cout <<*df.Count() << std::endl;
	
	//cut based on track number (4)
	Int_t track_num = 4;
	auto df_cut1 = cut_branch_single(df, track_num, "ntrk", "==", true);

	// Cut data based on low momentum threshold
	Float_t low_pt_threshold = 100000;
	auto df_cut2 = cut_branch_single(df_cut1, low_pt_threshold, "alltrk_pt", "<", true);

	// enforce charge neutrality
	Int_t trk_charge = 0;
	auto df_cut3 = cut_branch_sum(df_cut2, trk_charge, "trk_q", "==", true);

	
	auto df_cut4 = df_cut3.Define("trk_rap", [pion_mass](const ROOT::RVec<float>& pt, 
                             const ROOT::RVec<float>& eta) {
    ROOT::RVec<float> rap(pt.size());
    for (size_t i = 0; i < pt.size(); ++i) {
        float mt = sqrt(pion_mass * pion_mass + pt[i] * pt[i]);
        float num = sqrt(pion_mass * pion_mass + pt[i] * pt[i] * pow(cosh(eta[i]), 2)) + pt[i] * sinh(eta[i]);
        rap[i] = log(num / mt);
    }
    return rap;
	}, {"trk_pt", "trk_eta"});


	//find 4-momentum of pions and keep track of which ones are are +ve and -ve
	Int_t len3 = *df_cut3.Count();
	std::array<Int_t, 2> pion_neg;
	std::array<Int_t, 2> pion_pos;

	
	df_cut4.Foreach
	(
		[track_num, pion_mass, &pion_neg, &pion_pos, &graph]
		(const RVecI &trk_q, const RVecF &trk_pt, const RVecF &trk_eta, const RVecF &trk_phi)
		
		{
		Int_t pos_count = 0;
        Int_t neg_count = 0;
		std::array<TLorentzVector,4> pion_4p;
		for (Int_t i = 0; i < track_num; ++i)
			{
			switch(trk_q[i]){
			
			case +1:
				if (pos_count < 2){
					pion_pos[pos_count] = i;
					++pos_count;
					}
				break;

			case -1:
				if (neg_count < 2){
					pion_neg[neg_count] = i;
					++neg_count;
					}
				break;

			default:
				std::cout << " this is wrong" << std::endl;
			}
		if (pos_count == 2 && neg_count == 2)
		{for (Int_t i = 0; i < track_num; ++i){

		//	Float_t px = trk_pt[i] * TMath::Cos(trk_phi[i]);
		//	Float_t	py = trk_pt[i] * TMath::Sin(trk_phi[i]);
	    //	Float_t pz = trk_pt[i] * std::sinh(trk_eta[i]);
		// 	Float_t E = TMath::Sqrt(pion_mass*pion_mass + trk_pt[i]*trk_pt[i] + pz*pz);
		//	TLorentzVector p_pion(px, py, pz, E); 
			
			TLorentzVector p_pion;
            p_pion.SetPtEtaPhiM(trk_pt[i], trk_eta[i], trk_phi[i], pion_mass);
            pion_4p[i] = p_pion;

		
		}
		
		// The two possible combinations of pions that are charge neutral	
		TLorentzVector P_00 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[0]];
		TLorentzVector P_11 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[1]];
		
		TLorentzVector P_01 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[1]];
		TLorentzVector P_10 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[0]];	
		
		// Convert MeV for easier reading and comparison
	
		graph->Fill(P_00.M()*1e3, P_11.M()*1e3);	
		graph->Fill(P_01.M()*1e3, P_10.M()*1e3);	
	}
	}
			
		},
		{"trk_q", "trk_pt", "trk_eta", "trk_phi"}
	);


	//Y Projection
	Double_t x_1 = 484, x_2 = 504; // MeV/c^2
	Int_t x_bin_1 = graph->GetXaxis()->FindFixBin(x_1);
	Int_t x_bin_2 = graph->GetXaxis()->FindFixBin(x_2);
	TH1D *y_projection = graph->ProjectionY("Y_projection", x_bin_1, x_bin_2);
	y_projection->GetXaxis()->SetTitle("#rho_{1} (MeV / c^{2})");
	y_projection->GetYaxis()->SetTitle(Form("Events / %.2g MeV)",y_projection->GetXaxis()->GetBinWidth(1)));
	
	// Use breitwivner fit to find invarient mass of the resonances	
		
	Double_t fx_1 = 480, fx_2 = 515;
	TF1* fit_y = new TF1("fit_y", fit_gaussian, fx_1, fx_2, 4); // min, max num of params

    // fit to the histogram
	fit_y->SetParameters(495, 30, 100, 50); //create initial guesses for the params
	fit_y->SetParNames("Mass","Gamma", "normalisation", "constant");


	fit_y->SetParLimits(0, 492, 498);
	fit_y->SetParLimits(1, 5, 15);
	fit_y->SetParLimits(2, 160, 500);
	fit_y->SetParLimits(3, 50, 100);
	TFitResultPtr r = y_projection->Fit("fit_y", "S0");
	r->Print("V");
	
	auto legend_y = new TLegend(0.6,0.7,0.48,0.9);
	legend_y->SetHeader("Legend");
	legend_y->AddEntry("fit_y", "Gaussian fit", "l");
	legend_y->AddEntry("y_projection", "data", "lep");
	legend_y->Draw();


	y_projection->Draw("E*");
	fit_y->Draw("SAME");
	canvas->Update();
	canvas->SaveAs("media/photos/Younes_NTuples_TOTEM40/glueball_mass_reconstruction/kaons_from_rho_y_proj_all_40.pdf");
	canvas->Clear();
	
	TFile* outFile_y = new TFile("media/root_files/Younes_NTuples_TOTEM40/glueball_mass_reconstruction/kaons_from_rho_y_proj_all_40.root", "RECREATE");
	y_projection->Draw("E*");
	fit_y->Draw("SAME");
	legend_y->Draw();
    canvas->Update();
	canvas->Write();
	outFile_y->Close();
	canvas->Clear();
	

	//X projections
	Double_t y_1 = 484, y_2 = 504; // MeV/c^2
	Int_t y_bin_1 = graph->GetYaxis()->FindFixBin(y_1);
	Int_t y_bin_2 = graph->GetYaxis()->FindFixBin(y_2);
	TH1D *x_projection = graph->ProjectionX("X_projection", y_bin_1, y_bin_2);
	x_projection->GetXaxis()->SetTitle("#rho_{1} (MeV / c^{2})");
	x_projection->GetYaxis()->SetTitle(Form("Events (.%1f)",x_projection->GetXaxis()->GetBinWidth(1)));
	

	Double_t fy_1 = 480, fy_2 = 515;
	TF1* fit_x = new TF1("fit_x", fit_gaussian, fy_1, fy_2, 4); // min, max num of params

    // fit to the histogram
	fit_x->SetParameters(495, 30, 100, 50); //create initial guesses for the params
	fit_x->SetParNames("Mass","Gamma", "normalisation", "constant");
	fit_x->SetParLimits(0, 492, 498);
	fit_x->SetParLimits(1, 5, 15);
	fit_x->SetParLimits(2, 160, 500);
	fit_x->SetParLimits(3, 50, 100);
	TFitResultPtr r_x =x_projection->Fit("fit_x", "S0");
	r_x ->Print("V");
	
	
	x_projection->Draw("E*");
	fit_x->Draw("SAME");
	canvas->Update();
	canvas->SaveAs("media/photos/Younes_NTuples_TOTEM40/glueball_mass_reconstruction/kaons_from_rho_x_proj_all_40.pdf");
	canvas->Clear();
	
	TFile* outFile_x = new TFile("media/root_files/Younes_NTuples_TOTEM40/glueball_mass_reconstruction/kaons_from_rho_x_proj_all_40.root", "RECREATE");
	x_projection->Draw("E*");
	fit_x->Draw("SAME");
    canvas->Update();
	canvas->Write();
	outFile_x->Close();
    canvas->Clear();
	return 0;
}
