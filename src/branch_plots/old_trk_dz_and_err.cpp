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
 * Closest approach of the reconstucted tracjectory to the primary vertex,
 * with the z distance metric. Plot this diveded by its error.
 * We can use this to cut on tracks that are close to IP5
 *	
 *
 *	Results
 * Note that the resulting graph is a very peaked gaussian about (≈) 0. 
 * This may require further investigation. 
 */
int main()
{
    // constants
    Float_t pion_mass = 0.139; // GeV/c^2

    // Number of bins, range of values for the histogram
    Int_t n_bins = 400;
    Float_t z_min = -10.0; // unitless (ratio)
    Float_t z_max = 10.0; //  unitless (ratio)

    // Create histogram
    TH1F* graph = new TH1F("graph", "Z Distance / Error", n_bins, z_min, z_max);
    graph->GetXaxis()->SetTitle("#frac{z}{  z_{err}}");
    graph->GetYaxis()->SetTitle(Form("Events / %.1f (unitless)", (z_max - z_min) / n_bins));  
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
	(const Int_t track_num, const RVecF &dz, const RVecF &err)
	{
		for (Int_t i = 0; i < track_num; ++i){
			graph->Fill(dz[i] / err[i]);
		}
	}, {"ntrk", "trk_dz", "trk_dzerr"}
	);
	
	// Use more general gaussian fit from header file (includes scaling factor) to fit
	TF1* fit = new TF1("fit", fit_gaussian, -2, 2, 4); // min, max, number of parameters
	
	fit->SetParameters(0, 1e5, graph->GetMean(), graph->GetRMS());
	fit->SetParNames("Constant","scaling factor", "#mu", "#sigma");
	TFitResultPtr r = graph->Fit("fit", "S");
	r->Print("V");	

	// Update Canvas and Save
	graph->Draw("HISTE");
	fit->Draw("SAME");
    canvas->Update();
    canvas->SaveAs("media/photos/Younes_NTuples_TOTEM20/dz_over_err.pdf");
	canvas->Clear();

	TFile* outFile = new TFile("media/root_files/Younes_NTuples_TOTEM20/dz_over_err.root", "RECREATE");
	graph->Draw("HISTE");
	fit->Draw("SAME");
    canvas->Update();
	canvas->Write();
	outFile->Close();

	return 0;

}
