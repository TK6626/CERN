#include <typeinfo>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <TMath.h>
#include <array>

#include "TLine.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include "TTree.h"
#include "ROOT/RVec.hxx"
#include "root/TLorentzVector.h"


// include headers containing useful functions and definitions
// #include "../lib/computations"
#include "../lib/common_cuts.h"

/*
	Plots distribution of 4-momentum components of the initial state
	
	Cut1: cut data based on n (4) tracks
	Cut2: cut based on the the value of the transverse momentum
	Cut3: cut to ensure charge neutrality (4 meson decay)
	
*/
const float pion_mass =0.1396; //GeV/c^2


int main()
{   
	const char* file_name = "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM20.root";
	const char* tree_name = "tree";

    // Load in the Tree

	ROOT::EnableImplicitMT();
	ROOT::RDataFrame df(tree_name, file_name);

	for (const auto& name : df.GetColumnNames())
    	std::cout << name << std::endl;

	//cut based on track number (4)
	Int_t track_num = 4;
	auto df_cut1 = cut_branch_single(df, track_num, "ntrk", "==", true);

	// Cut data based on low momentum threshold
	Float_t low_pt_threshold = 0.25;
	auto df_cut2 = cut_branch_single(df_cut1, low_pt_threshold, "alltrk_pt", "<", true);

	// enforce charge neutrality
	Int_t trk_charge = 0;
	auto df_cut3 = cut_branch_sum(df_cut2, trk_charge, "trk_q", "==", true);


// Create canvases and graphs for graphs
	TCanvas *c1 = new TCanvas("c1");
TH1 *px_plot = new TH1F("px_plot", "p_{x}", 300, 7000, 0);
	px_plot->GetXaxis()->SetTitle("Momentum (MeV / c)");
	px_plot->GetYaxis()->SetTitle("Counts");
TH1 *py_plot = new TH1F("py_plot", "p_{y}", 300, 7000, 0);
	py_plot->GetXaxis()->SetTitle("Momentum (MeV / c)");
	py_plot->GetYaxis()->SetTitle("Counts");
TH1 *pz_plot = new TH1F("pz_plot", "p_{z}", 300, 7000, 0);
	pz_plot->GetXaxis()->SetTitle("Momentum (MeV / c)");
	pz_plot->GetYaxis()->SetTitle("Counts");
TH1 *E_plot = new TH1F("E_plot", "E", 300, 7000, 0);
	E_plot->GetXaxis()->SetTitle("Energy (MeV)");
	E_plot->GetYaxis()->SetTitle("Counts");
	
	//create a Rvector that contains the 4 momentum of each of the pions
	
	int n_plus = 0;
	int n_min = 0;
	df_cut3.Foreach(
	[track_num, &c1, &px_plot, &py_plot, &pz_plot, &E_plot]
	(const RVecF &trk_eta, const RVecF &trk_phi, const RVecF &trk_pt)
	{

	//instantiate the momentum for each row
	//this contains the 4 momenta for each track
	TLorentzVector total_p(0, 0, 0, 0);
	//loop through each track to calculate kinematical variables 
    for (int i = 0; i < track_num; ++i) {
		Float_t px = trk_pt[i] * TMath::Cos(trk_phi[i]);
		Float_t	py = trk_pt[i] * TMath::Sin(trk_phi[i]);
		Float_t pz = trk_pt[i] * std::sinh(trk_eta[i]);
		Float_t E = TMath::Sqrt(pion_mass*pion_mass + trk_pt[i]*trk_pt[i] + pz*pz);
		total_p += TLorentzVector(px, py, pz, E);
		}	

	// fill graph but CONVERT TO GeV ->MeV
	total_p *= 1e3;
		
	// fill histograms of 4-momentumj	 
	px_plot->Fill(total_p[0]); 
	py_plot->Fill(total_p[1]); 
	pz_plot->Fill(total_p[2]); 
	E_plot->Fill(total_p[3]);
}, {"trk_eta", "trk_phi", "trk_pt"}
);

E_plot->Draw("HIST");
c1->Update();
c1->SaveAs("media/photos/Younes_NTuples_TOTEM20/com_energy_distribution.pdf");
c1->Clear();

px_plot->Draw("HIST");
c1->Update();
c1->SaveAs("media/photos/Younes_NTuples_TOTEM20/com_px_distribution.pdf");
c1->Clear();

py_plot->Draw("HIST");
c1->Update();
c1->SaveAs("media/photos/Younes_NTuples_TOTEM20/com_py_distribution.pdf");
c1->Clear();

pz_plot->Draw("HIST");
c1->Update();
c1->SaveAs("media/photos/Younes_NTuples_TOTEM20/com_pz_distribution.pdf");			
c1->Clear();

return 0;

}

