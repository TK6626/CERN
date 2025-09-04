#include <stdio.h>
#include <cmath>
#include <array>

#include "TCanvas.h"
#include "TF1.h"
#include "TH2.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TStyle.h"

#include "../lib/common_cuts.h"
#include "../lib/computations.h"


/**
 * Plot the track azimuthal angle against the distance to the interaction
 * vertex in the xy plane
 *
 *
 * Reults:
 * The resulting graph is approximately a gaussia with resulting 
 * Chi^2 = n
 */

int main()
{
	Float_t pi = TMath::Pi();

    // Number of bins, range of x values for the histogram
    Int_t eta_bins = 200;
    Int_t z_bins = 200;
    Float_t eta_min = -pi; // radians
    Float_t eta_max = pi; // radians
    Float_t z_min = -0.15; // 
    Float_t z_max = 0.15; // 

    // Create histogram
    TH2F* graph = new TH2F("graph", "#eta against Z", z_bins, z_min, z_max, eta_bins, eta_min, eta_max);
    graph->GetXaxis()->SetTitle("Z");
    graph->GetYaxis()->SetTitle("#eta");
 	TCanvas* canvas = new TCanvas("canvas"); 
  
	// load in tree,get column names, enable multithreading
    const char* file_name = "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM20.root";
    // const char* file_name = "/afs/cern.ch/user/j/jloder/public/simple_cutted_data/TOTEM20simplecut.root";
    const char* tree_name = "tree";
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(tree_name, file_name);
    std::cout << df.Describe() << std::endl;
  
  
	df.Foreach(
	[graph]
	(const Int_t track_num, const RVecF &z, const RVecF &err, const RVecF &eta)
	{
		for (Int_t i = 0; i < track_num; ++i){
			graph->Fill(z[i], eta[i]);
		};
	}, {"ntrk", "trk_dz", "trk_dzerr", "trk_eta"}
	);

	// Use gaussian fitting function from computations header	
	//TF1* fit = new TF1("fit", fit_gaussian, -10, 10, 4); // min, max, number of parameters
	
	// fit to the histogram
	//fit->SetParameters(0, xy_distance_err->GetMean(),xy_distance_err->GetRMS());
	//fit->SetParNames("Constant","scaling factor", "#mu", "#sigma");
	//TFitResultPtr r = xy_distance_err->Fit("fit", "S");
	//r->Print("V");	

	graph->Draw("COLZ");
	canvas->Update();
    canvas->SaveAs("media/photos/Younes_NTuples_TOTEM20/z_eta.pdf");
	canvas->Clear();

	// Save canvas as a root File
	TFile* outFile = new TFile("media/root_files/Younes_NTuples_TOTEM20/z_eta.root", "RECREATE");
	gStyle->SetPalette(1);
   	graph->Draw("surf1");
    canvas->Update();
	canvas->Write();
	outFile->Close();

	return 0;

}
