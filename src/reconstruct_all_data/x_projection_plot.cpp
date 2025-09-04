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

#include "../../lib/computations.h"
#include "../../lib/plotting_style.h"
/** Loads in all geometries and plots the x projection 
 * the projections of the reconstructed 4 pion track
 *
 */


// !g++ src/reconstruct_all_data/x_projection_plot.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/x_projection_plot; ./bin/reconstruct_all_data/x_projection_plot
int main(){


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



	//X Projection (20 and 40)
	Double_t y_1_20 = 480, y_2_20 = 512; // MeV/c^2
	Int_t y_bin_1_20 = histo_20->GetXaxis()->FindFixBin(y_1_20);
	Int_t y_bin_2_20 = histo_20->GetXaxis()->FindFixBin(y_2_20);
	TH1D *x_projection_20 = histo_20->ProjectionX("x_projection_20", y_bin_1_20, y_bin_2_20);
	
	Double_t y_1_40 = 480, y_2_40 = 512; // MeV/c^2
	Int_t y_bin_1_40 = histo_40->GetXaxis()->FindFixBin(y_1_40);
	Int_t y_bin_2_40 = histo_40->GetXaxis()->FindFixBin(y_2_40);
	TH1D *x_projection_40 = histo_40->ProjectionX("x_projection_40", y_bin_1_40, y_bin_2_40);

  	x_projection_20->GetXaxis()->SetTitle("m (MeV/c^{2})");
  	x_projection_20->GetYaxis()->SetTitle(Form("Events (%.2f MeV/c^{2})", x_projection_20->GetBinWidth(1)));
	

	// Set up canvas
	TCanvas* canvas_x = new TCanvas("canvas_x", "Fitted Mass Peaks", 800, 600);
	canvas_x->SetGrid();
	canvas_x->SetMargin(0.12, 0.05, 0.12, 0.05); // left, right, bottom, top

	// Fit range
	Double_t fy_1 = 480, fy_2 = 515;

	TF1* fit_x_20 = new TF1("fit_x_20", fit_gaussian, fy_1, fy_2, 4);
	//fit_x_20->SetParameters(495, 10, 500, 300);
	fit_x_20->SetParNames("Mass","Gamma", "normalisation", "constant");
	fit_x_20->SetParLimits(0, 492, 498);
	fit_x_20->SetParLimits(1, 10, 15);
	//fit_x_20->SetParLimits(2, 160, 500);
	//fit_x_20->SetParLimits(3, 50, 200);
	fit_x_20->SetLineColor(kRed);
	fit_x_20->SetLineWidth(2);
	fit_x_20->SetLineStyle(2);
	TFitResultPtr r_20_x = x_projection_20->Fit("fit_x_20", "S0");
		
	// Set stxles for x_projection_20
	x_projection_20->SetMarkerStyle(20);
	x_projection_20->SetMarkerColor(red_tran);
	x_projection_20->SetLineColor(red_tran);

	

	TF1* fit_x_40 = new TF1("fit_x_40", fit_gaussian, fy_1, fy_2, 4);
	//fit_x_40->SetParameters(495, 10, 500, 300);
	fit_x_40->SetParNames("Mass","Gamma", "normalisation", "constant");
	fit_x_40->SetParLimits(0, 492, 498);
	fit_x_40->SetParLimits(1, 10, 15);
	//fit_x_40->SetParLimits(2, 160, 500);
	//fit_x_40->SetParLimits(3, 50, 200);
	fit_x_40->SetLineColor(kBlue);
	fit_x_40->SetLineWidth(2);
	fit_x_40->SetLineStyle(2);
	TFitResultPtr r_40_x = x_projection_40->Fit("fit_x_40", "S0");

	// Set stxles for x_projection_40
	x_projection_40->SetMarkerStyle(21);
	x_projection_40->SetMarkerColor(blue_tran);
	x_projection_40->SetLineColor(blue_tran);



	// Add entries to the legend AFTER setting styles
	auto legend_x = new TLegend(0.60, 0.70, 0.88, 0.88);
	// legend_x->SetHeader("Invariant Mass Fits", "C");
	legend_x->AddEntry(x_projection_20, "20 Geometry", "lep");
	legend_x->AddEntry(fit_x_20, "Gaussian Fit (20)", "l");
	legend_x->AddEntry(x_projection_40, "40 Geometry", "lep");
	legend_x->AddEntry(fit_x_40, "Gaussian Fit (40)", "l");

	// Draw projections and fits
	x_projection_20->Draw("E1P");
	fit_x_20->Draw("SAME");
	x_projection_40->Draw("E1P SAME");
	fit_x_40->Draw("SAME");

	// Draw the legend
	legend_x->SetTextSize(0.035);
	legend_x->Draw();

	canvas_x->Update();
	TFile* outFile_x = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/kaons_from_rho_x_proj_all_20_and_40.root", "RECREATE");
	canvas_x->Write();
	outFile_x->Close();

	canvas_x->SaveAs("media/photos/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/kaons_from_rho_x_proj_all_20_and_40.pdf");
	canvas_x->Clear();


	return 0;
}


