#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "ROOT/RDataFrame.hxx"
#include "TLatex.h"

#include<stdio.h>
#include <iomanip>

void plot_2D_Histogram(const char* axis_name,
  const std::vector<float>& theta_left,
  const std::vector<float>& theta_right,
  Int_t nbins = 100,
  Double_t min_x = -2e-4,
  Double_t max_x = 2e-4);

void right_left_correlation()
{

  /*
  Generate RDataFrame
  */ 

  const char* file_name = "data/TOTEM20.root";
  const char* tree_name = "tree;1";
  
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df(tree_name, file_name);
  auto theta_left_x = df.Take<Float_t>("ThxL");     
  auto theta_right_x = df.Take<Float_t>("ThxR");     
  auto theta_left_y = df.Take<Float_t>("ThyL");      
  auto theta_right_y = df.Take<Float_t>("ThyR");      
  

  plot_2D_Histogram("x", *theta_left_x, *theta_right_x, 1000);






}

void plot_2D_Histogram(const char* axis_name,
    const std::vector<float>& theta_left,
    const std::vector<float>& theta_right,
    Int_t nbins = 2000,
    Double_t min_x = -2e-4,
    Double_t max_x = 2e-4)
{
// Print min/max for diagnostics
auto [min, max] = std::minmax_element(theta_left.begin(), theta_left.end());
std::cout << std::setprecision(20) << "min = " << *min << ", max = " << *max << std::endl;

TCanvas* c1 = new TCanvas("c1");

// Create histogram before filling
TH2F* left_right_collinearity_histogram = new TH2F(
  "left_right_collinearity_histogram",
  "Right Left Collinearity",
  nbins, min_x, max_x, nbins, min_x, max_x);

// Fill histogram
size_t n = std::min(theta_left.size(), theta_right.size());
for (size_t i = 0; i < n; ++i) {
  left_right_collinearity_histogram->Fill(theta_left[i], theta_right[i]);
}


// Draw and save

left_right_collinearity_histogram->GetXaxis()->SetTitle("#theta_{#scale[0.9]{L}}*");
left_right_collinearity_histogram->GetYaxis()->SetTitle("#theta_{#scale[0.9]{R}}*");
left_right_collinearity_histogram->Draw("COLZ");
c1->Update();
c1->SaveAs("media/photos/Left_right_correlation_xaxis.pdf");
}




// zPV
// EventNum
// LumiSection
// Run
// ThxL
// ThxR
// ThyL
// ThyR
// alltrk_mass
// alltrk_pt
// ntrk
// trk_dedx
// trk_dedxerr
// trk_dxy
// trk_dz
// trk_eta
// trk_isK
// trk_isP
// trk_isPi
// trk_nMeasure
// trk_nMeasureLayer
// trk_nSaturMeasure
// trk_p
// trk_phi
// trk_pt
// trk_q