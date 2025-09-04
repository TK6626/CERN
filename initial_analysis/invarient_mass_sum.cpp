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

using RVecF = ROOT::VecOps::RVec<float>;
using RVecI = ROOT::VecOps::RVec<int>;
using array_4 = std::array<float, 4>;
using array_5 = std::array<float, 5>;
using RVecF_4 = ROOT::VecOps::RVec<array_4>;
using RVecF_5 = ROOT::VecOps::RVec<array_5>;

float calc_inv_mass(const array_4 &p);
array_4 add_4_vec(const array_4 &p1, const array_4 &p2);

/*
	calculates the specific energy loss against the rigidity 

	Cut1: cut data based on n (4) tracks
	Cut2: cut based on the the value of the transverse momentum
	
	Plots dE/dx against p 
*/
const float pion_mass =0.1396; //GeV/c^2

array_4 add_4_vec(const array_4 &p1, const array_4 &p2);
float calc_inv_mass(const array_4 &p);


void invarient_mass_sum()
{  
 
	const char* file_name = "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM20.root";
	const char* tree_name = "tree";


	// check the identifier type of the entries to know 
	// how to use RVec
	// TFile *f = TFile::Open(file_name);
	// TTree *t = (TTree*)f->Get(tree_name);
	// t->Print();

    // Load in the Tree

	ROOT::EnableImplicitMT();
	ROOT::RDataFrame df(tree_name, file_name);

	for (const auto& name : df.GetColumnNames())
    	std::cout << name << std::endl;

	//cut based on track number (4)
	const int track_num = 4;
	auto track_cut = [track_num](Int_t x) { return x == track_num; };
	auto df_cut1 = df.Filter(track_cut, {"ntrk"});
	std::cout <<"number of events with 4 tracks =" <<  *df_cut1.Count() << std::endl;

	// Cut data based on low momentum threshold
	Float_t low_pt_threshold = 0.25; 
	auto pt_cut = [low_pt_threshold](Float_t x){ return x < low_pt_threshold; };
	auto df_cut2 = df_cut1.Filter(pt_cut, {"alltrk_pt"});
	std::cout << "number of tracks with low momentum = " << *df_cut2.Count() << std::endl;


	// want to check that the event has equal amoubts of
	// both positive and negative charge
	auto neutral_event_cut = [] (const RVecI &trk_q) {return Sum(trk_q) == 0; };
	auto df_cut3 = df_cut2.Filter(neutral_event_cut, {"trk_q"});
	int num_neutral_events = *df_cut3.Count();
	std::cout << "number of 4 track events with low transverse momentum and neutral =" << num_neutral_events << std::endl;



	TCanvas *c1 = new TCanvas("c1");
	TH1 *mass_plot = new  TH1F("mass_plot", "Invarient Mass 4-track event", 500, 0, 4000);
	mass_plot->GetXaxis()->SetTitle("Mass (MeV / c^{2})");
	mass_plot->GetYaxis()->SetTitle("Counts");

	//create a Rvector that contains the 4 momentum of each of the pions
	
	int n_plus = 0;
	int n_min = 0;
	auto df_cut4 = df_cut3.Define("invariant_mass",
	[track_num]( 
    const RVecF &trk_eta, const RVecF &trk_phi, const RVecF &trk_pt) {
	
	//instantiate the momentum for each row
	//this contains the 4 momenta for each track
    std::array<array_4, track_num> p;
	
	//loop through each track to calculate kinematical variables 
    for (int i = 0; i < track_num; ++i) {
        float pt = trk_pt[i];
        float phi = trk_phi[i];
        float eta = trk_eta[i];
	
	//calculate the 4 momenta for each track
        float px = pt * std::cos(phi);
        float py = pt * std::sin(phi);
        float pz = pt * std::sinh(eta);
        float E  = std::sqrt(pt * pt + pz * pz + pion_mass * pion_mass);
	
		// fill the mometum array of 4-momenta
        p[i] = {E, px, py, pz};
    }
	
	//create an array that sums the total 4 momenta of all the tracks
    	array_4 total_p = {0, 0, 0, 0}; // [E, px, py, pz]
    	for (int i = 0; i < track_num; ++i) {
        	total_p[0] += p[i][0];  // E
        	total_p[1] += p[i][1];  // px
        	total_p[2] += p[i][2];  // py
        	total_p[3] += p[i][3];  // pz
    	}
	// calculat the invariant mass of the origional particle
	float invariant_mass = calc_inv_mass(total_p); 

    return invariant_mass;
}, {"trk_eta", "trk_phi", "trk_pt"});


// fill graph but CONVERT TO MeV/c^2
df_cut4.Foreach(
	[&mass_plot]
	(float &mass)
	{
	mass *= 1e3;
	mass_plot->Fill(mass);
	},
	{"invariant_mass"}
);	
std::cout << "Entries in mass_plot: " << mass_plot->GetEntries() << std::endl;


float y_min = mass_plot->GetMinimum();
float y_max = mass_plot->GetMaximum();
std::cout << "max counts =" << y_max << "\nmin counts =" << y_min << std::endl;

TLine *theoretical_mass = new TLine(2220, y_min, 2220, y_max);
theoretical_mass->SetLineColor(kRed);
theoretical_mass->SetLineStyle(2);
theoretical_mass->SetLineWidth(2);

TLatex *label = new TLatex(2225, y_max * 0.9, "Theoretical Mass");
label->SetTextColor(kRed);
label->SetTextSize(0.03);
label->Draw();


	mass_plot->Draw("HIST");	
	theoretical_mass->Draw("same");
	c1->Update();
	c1->SaveAs("media/photos/Younes_NTuples_TOTEM20/inavr_mass.pdf");
				
}

float calc_inv_mass(const array_4 &p) {
	return sqrt(p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3]);
}

array_4 add_4_vec(const array_4 &p1, const array_4 &p2){
	array_4 p;
	for (int i=0;i<4;++i){
		p[i] = p1[i]+p2[i];
		};
	return p;
	}



