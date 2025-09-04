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

#include "../../../lib/computations.h"
#include "../../../lib/plotting_params.h"

/** Loads in all geometries and plots the y projection of the rho reconstrution.
 * Kaons have been (mostly) removed with the help of the cuts on interaction verticies
 */


// !g++ src/reconstruct_all_data/rho_reconstruction/y_projection_breit_wigner_bg.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/rho_reconstruction/y_projection_breit_wigner_bg; ./bin/reconstruct_all_data/rho_reconstruction/y_projection_breit_wigner_bg

int main() {


	// set standard plotting styles
	SetPlotStyle();

	TFile* file20 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/histo_20.root", "READ");
	TFile* file40 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/histo_40.root", "READ");


	TH2F* histo_20 = (TH2F*)file20->Get("histo_20");
	TH2F* histo_40 = (TH2F*)file40->Get("histo_40");
	histo_20->SetTitle("Invariant #rho_{2} mass");

	// set up transparant colours for use later
	Int_t red_tran = TColor::GetFreeColorIndex();
	new TColor(red_tran, 1.0, 0.0, 0.0, "", 0.1); 
	Int_t blue_tran = TColor::GetFreeColorIndex();
	new TColor(blue_tran, 0.0, 0.0, 1.0, "", 0.1); 

	
	// make cut of a given number of sigma around m_p
	Float_t mass_constraint = 125; // MeV/c^2
	Float_t lower_bound = m_p - mass_constraint; // MeV/c^2
	Float_t upper_bound = m_p + mass_constraint; // MeV/c^2



	//Y Projection (20 and 40)
	Double_t x_1_20 = lower_bound, x_2_20 = upper_bound; // MeV/c^2
	Int_t y_bin_1_20 = histo_20->GetYaxis()->FindFixBin(x_1_20);
	Int_t y_bin_2_20 = histo_20->GetYaxis()->FindFixBin(x_2_20);
	TH1D *y_projection_20 = histo_20->ProjectionY("y_projection_20", y_bin_1_20, y_bin_2_20);
	
	Double_t x_1_40 = lower_bound, x_2_40 = upper_bound; // MeV/c^2
	Int_t y_bin_1_40 = histo_40->GetYaxis()->FindFixBin(x_1_40);
	Int_t y_bin_2_40 = histo_40->GetYaxis()->FindFixBin(x_2_40);
	TH1D *y_projection_40 = histo_40->ProjectionY("y_projection_40", y_bin_1_40, y_bin_2_40);

  	y_projection_20->GetYaxis()->SetTitle("m_{#pi^{+}#pi^{-}} (MeV/c^{2})");
  	y_projection_20->GetYaxis()->SetTitle(Form("Events (%.2f MeV/c^{2})", y_projection_20->GetBinWidth(1)));
	

	// Set up canvas
	TCanvas* canvas_y = new TCanvas("canvas_y", "Fitted Mass Peaks", 800, 600);
	canvas_y->SetGrid();
	canvas_y->SetMargin(0.12, 0.05, 0.12, 0.05); // left, right, bottom, top


// relativistic breit-wigner formula

	// Fit range
	Double_t fy_1 = 400, fy_2 = 1150; // what part of the projection is actully fitted to 

	TF1* fit_y_20 = new TF1("fit_y_20", fit_breit_wigner_bg, fy_1, fy_2, 8);
	fit_y_20->SetParameters(766, 150, 1.5e14, 100, 0.046, 296,  1.68, -0.0036);
	//fit_y_20->SetParNames("m","Γ", "A_0", "C");
	fit_y_20->SetParLimits(0, 745, 780);
	fit_y_20->SetParLimits(1, 150, 180);
	//fit_y_20->SetParLimits(2, 6.5e14, 10e15);
	//fit_y_20->SetParLimits(3, 50, 200);
	fit_y_20->SetLineColor(kRed);
	fit_y_20->SetLineWidth(2);
	fit_y_20->SetLineStyle(2);
	TFitResultPtr r_20_y = y_projection_20->Fit("fit_y_20", "S0");
		
	// Set styles for y_projection_20
	y_projection_20->SetMarkerStyle(20);
	y_projection_20->SetMarkerColor(red_tran);
	y_projection_20->SetLineColor(red_tran);

	

	TF1* fit_y_40 = new TF1("fit_y_40", fit_breit_wigner_bg, fy_1, fy_2, 8);
	fit_y_40->SetParameters(766, 150, 1.2e14, 100, 0.046, 296,  1.68, -0.0036);
	//fit_y_40->SetParNames("m","Γ", "A_0", "C");
	fit_y_40->SetParLimits(0, 745, 770);
	//fit_y_40->SetParLimits(1, 100, 180);
	//fit_y_40->SetParLimits(2, 160, 500);
	//fit_y_40->SetParLimits(3, 50, 200);
	fit_y_40->SetLineColor(kBlue);
	fit_y_40->SetLineWidth(2);
	fit_y_40->SetLineStyle(2);
	TFitResultPtr r_40_y = y_projection_40->Fit("fit_y_40", "S0");

	// Set styles for y_projection_40
	y_projection_40->SetMarkerStyle(21);
	y_projection_40->SetMarkerColor(blue_tran);
	y_projection_40->SetLineColor(blue_tran);


	// Add entries to the legend AFTER setting styles
	auto legend_y = new TLegend(0.70, 0.50, 0.95, 0.65); // bottom left, top right
	// legend_y->SetHeader("Invariant Mass Fits", "C");
	legend_y->AddEntry(y_projection_20, "20 Topology", "lep");
	legend_y->AddEntry(fit_y_20, "Gaussian Fit (20)", "l");
	legend_y->AddEntry(y_projection_40, "40 Topology", "lep");
	legend_y->AddEntry(fit_y_40, "Gaussian Fit (40)", "l");

	// Draw projections and fits
	y_projection_20->Draw("E1P");
	fit_y_20->Draw("SAME");
	y_projection_40->Draw("E1P SAME");
	fit_y_40->Draw("SAME");

	// Draw the legend
	legend_y->SetTextSize(0.035);
	legend_y->Draw();


	// manually set ayis range if not drawn properly
	y_projection_20->SetAxisRange(50, 700, "Y");

	canvas_y->Update();
	TFile* outFile_y = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/kaons_from_rho_y_proj_all_20_and_40_bg.root", "RECREATE");
	canvas_y->Write();
	outFile_y->Close();

	canvas_y->SaveAs("media/photos/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/kaons_from_rho_y_proj_all_20_and_40.pdf");
	canvas_y->Clear();


	return 0;
}



