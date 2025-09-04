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

const float pi = 3.14159;
/*

calculate the specific energy loss against phi

*/
void dedx_phi()
{  

 
  	const char* file_name= "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM20.root";
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


	//Want to plot de/dx agains trk_phi

    TCanvas *c1 = new TCanvas("c1", "dE/dx vs #phi");

    // Create the 2D histogram
    TH2F *plot_dedx = new TH2F("plot_dedx", "dE/dx vs #phi;dE/dx (MeV/cm); #phi",
                               1000, - pi, pi,  // x-axis bins and range
                               1000, -pi, pi); // y-axis bins and range


	df_cut2.Foreach(
	[plot_dedx, track_num] (const RVecF &E, const RVecF &phi)
	{
		for (int i = 0; i < track_num; ++i)
			plot_dedx->Fill(E[i], phi[i]);	
	},
	{"trk_dxy", "trk_phi"} 
	);

// plot the graph
plot_dedx->Draw("COLZ");
c1->Update();
c1->SaveAs("media/photos/dxy_phi.pdf");
		 
}
