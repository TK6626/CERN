#include <typeinfo>
#include <stdio.h>
#include <iostream>
#include <vector>
#include "TMath.h"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include "TTree.h"
#include "ROOT/RVec.hxx"

using RVecF = ROOT::VecOps::RVec<float>;

/*
	calculates the specific energy loss against the rigidity 

	Cut1: cut data based on n (4) tracks
	Cut2: cut based on the the value of the transverse momentum
	
	Plots dE/dx against p 
*/
void dedx_p()
{  

 
  	const char* file_name= "data/TOTEM20.root";
    const char* tree_name = "tree;1";
	
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
	auto track_cut = [](Int_t x) { return x == 4; };
	auto df_cut1 = df.Filter(track_cut, {"ntrk"});
	std::cout <<"number of events with 4 tracks =" <<  *df_cut1.Count() << std::endl;

	// Cut data based on low momentum threshold
	Float_t low_pt_threshold = 0.25;
	auto pt_cut = [low_pt_threshold](Float_t x){ return x < low_pt_threshold; };
	auto df_cut2 = df_cut1.Filter(pt_cut, {"alltrk_pt"});
	std::cout << "number of tracks with low momentum = " << *df_cut2.Count() << std::endl;


	//Want to plot  momemetum against de/dx

    TCanvas *c1 = new TCanvas("c1", "dE/dx vs p");

    // Create the 2D histogram
    TH2F *plot_dedx = new TH2F("plot_dedx", "dE/dx vs p;Momentum (MeV);dE/dx (MeV/cm)",
                               1000, 0, 10,  // x-axis bins and range
                               1000, 0, 10); // y-axis bins and range

    // Fill the histogram from RVecs
    df_cut2.Foreach(
        [plot_dedx] (const RVecF &p, const RVecF &dedx) 
        {int n = p.size();
		for (int i = 0; i < n; ++i)
    	plot_dedx->Fill(p[i], dedx[i]);
		},
	{"trk_p", "trk_dedx"}
);

plot_dedx->Draw("COLZ");
c1->Update();
c1->SaveAs("media/photos/momentum_against_dedx.pdf");
		 
}
