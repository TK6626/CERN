#include <stdio.h>
#include <cmath>
#include <array>
#include <cmath>

#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TFitResult.h"

#include "../lib/common_cuts.h"
#include "../lib/computations.h"


/**
 * Look at distribution of energy deposited
 *
 *
 * Reults:
 */
int main()
{

    // Number of bins, range of x values for the histogram
    Int_t n_bins = 100;
    Float_t x_min = 0; // GeV / distance (?)
    Float_t x_max = 4; // GeV / (distance ?)

  
	// load in tree,get column names, enable multithreading
    //const char* file_name = "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM20.root";
    const char* file_name = "/afs/cern.ch/user/j/jloder/public/simple_cutted_data/TOTEM20simplecut.root";
    const char* tree_name = "tree";
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(tree_name, file_name);
    std::cout << df.Describe() << std::endl;
  	
    // Create histogram
 	TCanvas* canvas = new TCanvas("canvas"); 
    TH1F* graph = new TH1F("graph", "Specific Energy Loss ", n_bins, x_min, x_max);
    graph->GetXaxis()->SetTitle("#frac{dE}{dx}");
    graph->GetYaxis()->SetTitle(Form("Events / %.1f (GeV/cm)", (x_max - x_min) / n_bins));  
	

    std::mutex mtx;
	// initialise bins and their error
	RVecF bin_total(n_bins, 0);
	RVecF bin_err(n_bins, 0);	
    const Double_t bin_width = (x_max - x_min) / n_bins;
	
	df.Foreach(
	[&]
	(const RVecF& dedx, const RVecF& dedx_err)
	{
		for (size_t i = 0; i < dedx.size(); ++i) {
        float x = dedx[i];
        float err = dedx_err[i];

        Int_t bin = static_cast<Int_t>((x - x_min) / bin_width);
        if (bin < 0 || bin >= n_bins) continue;
		
		Float_t err1 = err/5000;
        std::lock_guard<std::mutex> lock(mtx);
        bin_total[bin] += 1.0;
        bin_err[bin]     +=TMath::Abs(err1 * err1) ;
        }
    }, {"trk_dedx", "trk_dedxerr"}); 


	for (int i = 0; i < n_bins; ++i) {
        graph->SetBinContent(i, bin_total[i]);
        graph->SetBinError(i, TMath::Sqrt(bin_err[i]));
    	
		std::cout << graph->GetBinContent(i) <<std::endl;
		std::cout << graph->GetBinError(i) <<std::endl;
	}
	for (int i = 1; i <= graph->GetNbinsX(); ++i) {
    double err = graph->GetBinError(i);
    
    // Check for NaN or inf
    if (!std::isfinite(err) || err < 0) {
        graph->SetBinError(i, 0);  // Set to 0 or another default safe value
    }
}
	float total_entries = Sum(bin_total);
	std::cout << "Total entries across bins: " << total_entries << std::endl;
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
    canvas->SaveAs("media/photos/initial_TOTEM20/dedx_distribution.pdf");
	canvas->Clear();

	TFile* outFile = new TFile("media/root_files/initial_TOTEM20/dedx_distribution.root", "RECREATE");
	graph->Draw("HISTE");
	fit->Draw("SAME");
    canvas->Update();
	canvas->Write();
	outFile->Close();

    return 0;
}
