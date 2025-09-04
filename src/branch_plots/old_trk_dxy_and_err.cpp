#include <stdio.h>
#include <cmath>
#include <array>

#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../lib/common_cuts.h"
#include "../lib/computations.h"


/**
 * d_xy is the "Closest approach of reconstructed trajectory 
 * to primary vertex, with the x,y distance metric" - Younes E
 *
 * Plot this agains divided by the error in the distance 
 *
 *
 * Reults:
 * The resulting graph is approximately a gaussia with resulting 
 * Chi^2 = n
 */

int main()
{
    // constants
    Float_t pion_mass = 0.139; // GeV/c^2

    // Number of bins, range of x values for the histogram
    Int_t n_bins = 400;
    Float_t x_min = -10.0; // Unitless (ratio)
    Float_t x_max = 10.0; // Unitless (ratio)

    // Create histogram
    TH1F* graph = new TH1F("graph", "XY Distance / Error", n_bins, x_min, x_max);
    graph->GetXaxis()->SetTitle("#frac{xy}{  xy_{err}}");
    graph->GetYaxis()->SetTitle(Form("Events / %.1f (unitless)", (x_max - x_min) / n_bins));  
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
	(const Int_t track_num, const RVecF &dxy, const RVecF &err)
	{
		for (Int_t i = 0; i < track_num; ++i){
			graph->Fill(dxy[i] / err[i]);
		};
	}, {"ntrk", "trk_dxy", "trk_dxyerr"}
	);
	
	// Use gaussian fitting function from computations header	
	TF1* fit = new TF1("fit", fit_gaussian, -10, 10, 4); // min, max, number of parameters
	
	// fit to the histogram
	fit->SetParameters(0, 1e5, graph->GetMean(), graph->GetRMS());
	fit->SetParNames("Constant","scaling factor", "#mu", "#sigma");
	TFitResultPtr r = graph->Fit("fit", "S");
	r->Print("V");	

	// Update Canvas and Save
	graph->Draw("HISTE");
	fit->Draw("SAME");
    canvas->Update();
    canvas->SaveAs("media/photos/Younes_NTuples_TOTEM20/dxy_over_err.pdf");
	canvas->Clear();

	TFile* outFile = new TFile("media/root_files/Younes_NTuples_TOTEM20/dxy_over_err.root", "RECREATE");
	graph->Draw("HISTE");
	fit->Draw("SAME");
    canvas->Update();
	canvas->Write();
	outFile->Close();

	return 0;
}
