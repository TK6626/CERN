#include <array>
#include <stdio.h>
#include <cmath>

#include "TMath.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH2.h"


using array4 = std::array<Float_t, 4>; 
using RVecF = ROOT::VecOps::RVec<Float_t>;
using RVecD = ROOT::VecOps::RVec<Double_t>;
using RVecI = ROOT::VecOps::RVec<Int_t>;



/*
 * 	Construct the mass of the 2-pair mesons from the 4-meson tracks 
 * 
 *	cut1:
 *	cut2:
 *	cut3:
 *
 *
 * */

void reconstruct_mass(){
	
	// constants
	Float_t pion_mass = 0.139; // GeV/c^2

	// Create graph and canvas for later
	
	TH2* pair_meson_plot = new  TH2F("pair_meson_plot", "Reconstructed #rho mesons;m_{#rho_{1}} (MeV / c^{2} ;m_{#rho_{2}} (MeV / c^{2})",
								100, 240, 2000, 100, 240, 2000);
	TCanvas* canvas = new TCanvas("canvas");



	// load in tree,get column names, enable multithreading
    const char* file_name = "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM20.root";
    // const char* file_name = "/afs/cern.ch/user/j/jloder/public/simple_cutted_data/TOTEM20simplecut.root";
	const char* tree_name = "tree";
	ROOT::EnableImplicitMT();
	ROOT::RDataFrame df(tree_name, file_name);
	std::cout << df.Describe() << std::endl;
	
	
	//create cut on data to ensure 4-tracks	
	Int_t track_num = 4;
	auto track_cut = [track_num](Int_t trk_n){return trk_n == track_num;};
	auto df_cut1 = df.Filter(track_cut, {"ntrk"});
	std::cout << "number of entries with" << track_num << " tracks is " << *df_cut1.Count() << endl; 	
	
	// create cut on transverse momentum (low)
	Float_t pt_cutoff = 1.5;
	auto pt_cut = [pt_cutoff](const RVecF &trk_pt){return Sum(trk_pt) < pt_cutoff;};
	auto df_cut2 = df_cut1.Filter(pt_cut, {"trk_pt"});
	std::cout << "number of entries with pt below" << pt_cutoff << " is " << *df_cut2.Count() << endl; 	
	
	// create cut on charge neutrality 
	auto charge_cut = [](const RVecI &trk_q){return Sum(trk_q) == 0;};
	auto df_cut3 = df_cut2.Filter(charge_cut, {"trk_q"});
	std::cout << "number of charge neutral events is " << *df_cut3.Count() << endl; 	
	
	// introduce cut on trk_dxy / trk_dxyerr
	// introduce cut on trk_dz / trk_z err
	// cut on zPV
	// trk_phi / trk_dxy
	// trk_eta / trk_d


	
	//find 4-momentum of pions and keep track of which ones are are +ve and -ve
	Int_t len3 = *df_cut3.Count();
	std::array<Int_t, 2> pion_neg;
	std::array<Int_t, 2> pion_pos;

	
	df_cut3.Foreach
	(
		[track_num, pion_mass, &pion_neg, &pion_pos, &pair_meson_plot]
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

			Float_t px = trk_pt[i] * TMath::Cos(trk_phi[i]);
			Float_t	py = trk_pt[i] * TMath::Sin(trk_phi[i]);
			Float_t pz = trk_pt[i] * std::sinh(trk_eta[i]);
			Float_t E = TMath::Sqrt(pion_mass*pion_mass + trk_pt[i]*trk_pt[i] + pz*pz);
			TLorentzVector p_pion(px, py, pz, E); 
			pion_4p[i] = p_pion;	
			}		
		
		// The two possible combinations of pions that are charge neutral	
		TLorentzVector P_00 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[0]];
		TLorentzVector P_11 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[1]];
		
		TLorentzVector P_01 = pion_4p[pion_pos[0]] + pion_4p[pion_neg[1]];
		TLorentzVector P_10 = pion_4p[pion_pos[1]] + pion_4p[pion_neg[0]];	
		
		// Convert MeV for easier reading and comparison
	
		pair_meson_plot->Fill(P_00.M()*1e3, P_11.M()*1e3);	
		pair_meson_plot->Fill(P_10.M()*1e3, P_10.M()*1e3);	
			
		},
		{"trk_q", "trk_pt", "trk_eta", "trk_phi"}
	);
		 
	pair_meson_plot->Draw("COLZ");
	canvas->Update();
	canvas->SaveAs("media/photos/Younes_NTuples_TOTEM20/reconstucted_mass_pair_mesons.pdf");


}
