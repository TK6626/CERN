#include <stdio.h>
#include <cmath>
#include <array>

#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"


using VecF = ROOT::VecOps::RVec<Float_t>;

void trk_dxy_and_err()
{
    // constants
    Float_t pion_mass = 0.139; // GeV/c^2

    // Number of bins, range of x values for the histogram
    int n_bins = 400;
    float x_min = -10.0;
    float x_max = 10.0;

    // Create histogram
    TH1F* xy_distance_err = new TH1F("xy_distance_err", "XY Distance Error Distribution", n_bins, x_min, x_max);
    xy_distance_err->GetXaxis()->SetTitle("#frac{xy}{  xy_{err}}");
    xy_distance_err->GetYaxis()->SetTitle("Counts");  
 	TCanvas* canvas = new TCanvas("canvas"); 
  
	// load in tree,get column names, enable multithreading
    const char* file_name = "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM20.root";
    // const char* file_name = "/afs/cern.ch/user/j/jloder/public/simple_cutted_data/TOTEM20simplecut.root";
    const char* tree_name = "tree";
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(tree_name, file_name);
    std::cout << df.Describe() << std::endl;
  
  
	df.Foreach(
	[xy_distance_err]
	(const Int_t track_num, const VecF &dxy, const VecF &err)
	{
		for (Int_t i = 0; i < track_num; ++i){
			xy_distance_err->Fill(dxy[i] / err[i]);
		};
	}, {"ntrk", "trk_dxy", "trk_dxyerr"}
	);
	

	xy_distance_err->Draw("COLZ");
    canvas->Update();
    canvas->SaveAs("media/photos/Younes_NTuples_TOTEM20/dxy_over_err.pdf");

}
