#include <stdio.h>
#include <cmath>
#include <array>

#include "TCanvas.h"
#include "TH1.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TStyle.h"
#include "TF1.h"
#include "TFitResult.h"
#include "../lib/common_cuts.h"
#include "../lib/computations.h"


/**
 * zPV is the z distance to the primary vertex
 *
 *
 * Reults:
 */
int main()
{

    // Number of bins, range of x values for the histogram
    Int_t n_bins = 400;
    Float_t z_min = -10.0; // cm (? is distance)
    Float_t z_max = 10.0; // cm (? is a distance)

    // Create histogram
    TH1F* graph = new TH1F("graph", "Z Distance / Error", n_bins, z_min, z_max);
    graph->GetXaxis()->SetTitle("Z distance to Primary Vertex Distribution");
    graph->GetYaxis()->SetTitle(Form("Events / %.1f (cm?)", (z_max - z_min) / n_bins));  
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
	(const Float_t zPV)
	{
			graph->Fill(zPV);
	}, {"zPV"}
	);

	// Use gaussian fitting function from computations header	
	TF1* fit = new TF1("fit", fit_gaussian, z_min, z_max, 4); // min, max, number of parameters
	
	fit->SetParameters(0, 1e5, graph->GetMean(), graph->GetRMS());
	fit->SetParNames("Constant","scaling factor", "#mu", "#sigma");
	TFitResultPtr r = graph->Fit("fit", "S");
	r->Print("V");	

	// Update Canvas and Save
	graph->Draw("HISTE");
	fit->Draw("SAME");
    canvas->Update();
    canvas->SaveAs("media/photos/Younes_NTuples_TOTEM20/z_vertex_point.pdf");
	canvas->Clear();

	TFile* outFile = new TFile("media/root_files/Younes_NTuples_TOTEM20/z_vertex_point.root", "RECREATE");
	graph->Draw("HISTE");
	fit->Draw("SAME");
    canvas->Update();
	canvas->Write();
	outFile->Close();

	return 0;
}
