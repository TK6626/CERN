#include <typeinfo>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <TMath.h>
#include <array>

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
const float pion_mass = 130;//MeV/c^2

array_4 add_4_vec(const array_4 &p1, const array_4 &p2);
float calc_inv_mass(const array_4 &p);


void invarient_mass()
{  
 
  	//const char* file_name= "data/TOTEM20.root";
    //const char* tree_name = "tree;1";
	
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
	int track_num = 4;
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
	auto neutral_event_cut = [] (const RVecI &trk_q) {
	return Sum(trk_q) == 0; };

	auto df_cut3 = df_cut2.Filter(neutral_event_cut, {"trk_q"});
	int num_neutral_events = *df_cut3.Count();
	std::cout << "number of 4 track events with low transverse momentum and neutral =" << num_neutral_events << std::endl;



	TCanvas *c1 = new TCanvas("c1");
	TH2 *m = new TH2F("m", "Invarient Mass; ", 400, 259.98, 260.11, 400, 259.98, 260.11);
	TH1 *px_plot = new TH1F("px_plot", "X Momentum Distribution;#p_x (MeV/c);Number", 1000, -5, 5);
	TH1 *py_plot = new TH1F("py_plot", "X Momentum Distribution;#p_y (MeV/c);Number", 100, -1e3, 1e3);
	TH1 *pz_plot = new TH1F("pz_plot", "X Momentum Distribution;#p_z (MeV/c);Number", 100, -1e3, 1e3);
	TH1 *E_plot = new TH1F("E_plot", "X Momentum Distribution;#E (MeV);Number", 100, -1e3, 1e3);

	//create a Rvector that contains the 4 momentum of each of the pions
	
	int n_plus = 0;
	int n_min = 0;
	RVecF_4 pion_plus;
	RVecF_4 pion_min;

	
	df_cut3.Foreach(
		[&m, &px_plot](const RVecI &trk_q, const RVecF &trk_eta, const RVecF &trk_phi, const RVecF &trk_pt){
			RVecF_4 pion_plus, pion_min;
			int n_plus = 0, n_min = 0;
			for (int i = 0; i < trk_q.size(); ++i) {
				float px = trk_pt[i]*cos(trk_phi[i]);
				float py = trk_pt[i]*sin(trk_phi[i]);
				float pz = trk_pt[i]*sinh(trk_eta[i]);
				float E = sqrt(pion_mass*pion_mass + trk_pt[i]*trk_pt[i] + pz*pz);	
				if (trk_q[i] == 1) {
					pion_plus.push_back({E, px, py, pz});
					n_plus++;
				} else if (trk_q[i] == -1) {
					pion_min.push_back({E, px, py, pz});
					n_min++;
				}
			}
			// std::cout << " number of +ve muons = " << n_plus << std::endl;
			// std::cout << " number of -ve muons = " << n_min << std::endl;
	
			array_4 p11 = add_4_vec(pion_plus[0], pion_min[0]);
			array_4 p12 = add_4_vec(pion_plus[1], pion_min[0]);
			array_4 p21 = add_4_vec(pion_plus[0], pion_min[1]);
			array_4 p22 = add_4_vec(pion_plus[1], pion_min[1]);

			float p11_mass = calc_inv_mass(p11);
			float p12_mass = calc_inv_mass(p12);
			float p21_mass = calc_inv_mass(p21);
			float p22_mass = calc_inv_mass(p22);

			m->Fill(p11_mass, p22_mass);
			m->Fill(p12_mass, p21_mass);
			px_plot->Fill(p11[1] + p22[1] + p12[1] + p22[1]);

		}, {"trk_q", "trk_eta", "trk_phi", "trk_pt"});
		m->Draw("COLZ");
		c1->Update();
		c1->SaveAs("media/photos/invarient_mass.pdf");
				
		c1->Clear();
		px_plot->Draw("AL");
		c1->Update();
		c1->SaveAs("media/photos/total_x_momentum.pdf");	
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



