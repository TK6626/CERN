#include <array>
#include <stdio.h>
#include <cmath>
#include <vector>

#include "TMath.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TColor.h"

#include "../../lib/plotting_style.h"
#include "../../lib/computations.h"

/** get all data and plot the y projections of the
 * reconstructed 4 track pion data
 */

// !g++ src/reconstruct_all_data/y_projection_plot.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/y_projection_plot; ./bin/reconstruct_all_data/y_projection_plot


int main() {

// set standard plotting styles
SetPlotStyle();

TFile* file20 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/histo_20.root", "READ");
TFile* file40 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/histo_40.root", "READ");


TH2F* histo_20 = (TH2F*)file20->Get("histo_20");
TH2F* histo_40 = (TH2F*)file40->Get("histo_40");

// set up transparant colours for use later
Int_t red_tran = TColor::GetFreeColorIndex();
new TColor(red_tran, 1.0, 0.0, 0.0, "", 0.1); 
Int_t blue_tran = TColor::GetFreeColorIndex();
new TColor(blue_tran, 0.0, 0.0, 1.0, "", 0.1); 


	// Set up canvas
	TCanvas* canvas = new TCanvas("canvas", "Fitted Mass Peaks", 800, 600);
	canvas->SetGrid();
	canvas->SetMargin(0.12, 0.05, 0.12, 0.05); // left, right, bottom, top


	//Y Projection (20 and 40)
	Double_t x_1_20 = 475, x_2_20 = 520; // MeV/c^2
	Int_t x_bin_1_20 = histo_20->GetXaxis()->FindFixBin(x_1_20);
	Int_t x_bin_2_20 = histo_20->GetXaxis()->FindFixBin(x_2_20);
	TH1D *y_projection_20 = histo_20->ProjectionY("y_projection_20", x_bin_1_20, x_bin_2_20);
	
	Double_t x_1_40 = 475, x_2_40 = 520; // MeV/c^2
	Int_t x_bin_1_40 = histo_40->GetXaxis()->FindFixBin(x_1_40);
	Int_t x_bin_2_40 = histo_40->GetXaxis()->FindFixBin(x_2_40);
	TH1D *y_projection_40 = histo_40->ProjectionY("y_projection_40", x_bin_1_40, x_bin_2_40);

  	y_projection_20->GetXaxis()->SetTitle("m (MeV/c^{2})");
  	y_projection_20->GetYaxis()->SetTitle(Form("Events (%.2f MeV/c^{2})", y_projection_20->GetBinWidth(1)));


	// Fit range
	Double_t fx_1 = 482, fx_2 = 510;

	TF1* fit_y_20 = new TF1("fit_y_20", fit_gaussian, fx_1, fx_2, 4);
	fit_y_20->SetParameters(495, 10, 600, 300);
	fit_y_20->SetParNames("Mass","Gamma", "normalisation", "constant");
	fit_y_20->SetParLimits(0, 497, 500);
	fit_y_20->SetParLimits(1, 5, 15);
	//fit_y_20->SetParLimits(2, 160, 600);
	//fit_y_20->SetParLimits(3, 290, 370);
	fit_y_20->SetLineColor(kRed);
	fit_y_20->SetLineWidth(2);
	fit_y_20->SetLineStyle(2);
	TFitResultPtr r_20_y = y_projection_20->Fit("fit_y_20", "S0");
	
	// Set styles for y_projection_20
	y_projection_20->SetMarkerStyle(20);
	y_projection_20->SetMarkerColor(red_tran);
	y_projection_20->SetLineColor(red_tran);



	TF1* fit_y_40 = new TF1("fit_y_40", fit_gaussian, fx_1, fx_2, 4);
	fit_y_40->SetParameters(495, 12, 500, 300);
	fit_y_40->SetParNames("Mass","Gamma", "normalisation", "constant");
	fit_y_40->SetParLimits(0, 497, 500);
	fit_y_40->SetParLimits(1, 5, 15);
	//fit_y_40->SetParLimits(2, 290, 600);
	//fit_y_40->SetParLimits(3, 290, 370);
	fit_y_40->SetLineColor(kBlue);
	fit_y_40->SetLineWidth(2);
	fit_y_40->SetLineStyle(2);
	TFitResultPtr r_40_y = y_projection_40->Fit("fit_y_40", "S0");


	// Set styles for y_projection_40
	y_projection_40->SetMarkerStyle(21);
	y_projection_40->SetMarkerColor(blue_tran);
	y_projection_40->SetLineColor(blue_tran);

	// Add entries to the legend setting styles
	auto legend_y = new TLegend(0.60, 0.70, 0.88, 0.88);
	legend_y->SetHeader("Invariant Mass Fits", "C");
	legend_y->AddEntry(y_projection_20, "20 Geometry", "lep");
	legend_y->AddEntry(fit_y_20, "Gaussian Fit (20)", "l");
	legend_y->AddEntry(y_projection_40, "40 Geometry", "lep");
	legend_y->AddEntry(fit_y_40, "Gaussian Fit (40)", "l");

	// Draw projections and fits
	y_projection_20->Draw("E1");
	fit_y_20->Draw("SAME");
	y_projection_40->Draw("E1 SAME");
	fit_y_40->Draw("SAME");

	// Draw the legend
	legend_y->SetTextSize(0.035);
	legend_y->Draw();

	canvas->Update();
	TFile* outFile_y = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/kaons_from_rho_y_proj_all_20_and_40.root", "RECREATE");
	canvas->Write();
	outFile_y->Close();

	canvas->SaveAs("media/photos/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/kaons_from_rho_y_proj_all_20_and_40.pdf");
	canvas->Clear();

return 0;
}
