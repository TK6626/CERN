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

/** Loads in all geometries and plots the x projection of the rho reconstrution.
 * Kaons have been (mostly) removed with the help of the cuts on interaction verticies
 */


// !g++ src/reconstruct_all_data/rho_reconstruction/x_projection_breit_wigner_bg.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/rho_reconstruction/x_projection_breit_wigner_bg; ./bin/reconstruct_all_data/rho_reconstruction/x_projection_breit_wigner_bg

int main() {


	// set standard plotting styles
	SetPlotStyle();

	TFile* file20 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/histo_20.root", "READ");


	TH2F* histo_20 = (TH2F*)file20->Get("histo_20");
	

	histo_20->SetTitle("Invariant #rho_{1} mass");

	// set up transparant colours for use later
	Int_t red_tran = TColor::GetFreeColorIndex();
	new TColor(red_tran, 1.0, 0.0, 0.0, "", 0.1); 
	Int_t blue_tran = TColor::GetFreeColorIndex();
	new TColor(blue_tran, 0.0, 0.0, 1.0, "", 0.1); 

	
	// make cut of a given number of sigma around m_p
	Float_t mass_constraint = 125; // MeV/c^2
	Float_t lower_bound = m_p - mass_constraint; // MeV/c^2
	Float_t upper_bound = m_p + mass_constraint; // MeV/c^2



	//X Projection (20 and 40)
	Double_t y_1_20 = lower_bound, y_2_20 = upper_bound; // MeV/c^2
	Int_t y_bin_1_20 = histo_20->GetXaxis()->FindFixBin(y_1_20);
	Int_t y_bin_2_20 = histo_20->GetXaxis()->FindFixBin(y_2_20);
	TH1D *x_projection_20 = histo_20->ProjectionX("x_projection_20", y_bin_1_20, y_bin_2_20);
	
	Double_t x_1_40 = lower_bound, x_2_40 = upper_bound; // MeV/c^2
	Int_t x_bin_1_40 = histo_20->GetXaxis()->FindFixBin(x_1_40);
	Int_t x_bin_2_40 = histo_20->GetXaxis()->FindFixBin(x_2_40);
	TH1D *x_projection_20 = histo_20->ProjectionX("x_projection_20", x_bin_1_40, x_bin_2_40);

  	x_projection_20->GetXaxis()->SetTitle("m_{#pi^{+}#pi^{-}} (MeV/c^{2})");
  	x_projection_20->GetYaxis()->SetTitle(Form("Events (%.2f MeV/c^{2})", x_projection_20->GetBinWidth(1)));
	

	// Set up canvas
	TCanvas* canvas_x = new TCanvas("canvas_x", "Fitted Mass Peaks", 800, 600);
	canvas_x->SetGrid();
	canvas_x->SetMargin(0.12, 0.05, 0.12, 0.05); // left, right, bottom, top


// relativistic breit-wigner formula

	// Fit range
	Double_t fy_1 = 400, fy_2 = 1150; // what part of the projection is actully fitted to 

	TF1* fit_x_20 = new TF1("fit_x_20", fit_breit_wigner_bg, fy_1, fy_2, 8);
	fit_x_20->SetParameters(750, 150, 1e14, 100, 0.046, 296,  1.68, -0.0036);
	//fit_x_20->SetParNames("m","Γ", "A_0", "C");
	fit_x_20->SetParLimits(0, 745, 780);
	fit_x_20->SetParLimits(1, 150, 180);
	//fit_x_20->SetParLimits(2, 6.5e14, 10e15);
	//fit_x_20->SetParLimits(3, 50, 200);
	fit_x_20->SetLineColor(kRed);
	fit_x_20->SetLineWidth(2);
	fit_x_20->SetLineStyle(2);
	TFitResultPtr r_20_x = x_projection_20->Fit("fit_x_20", "S0");
		
	// Set stxles for x_projection_20
	x_projection_20->SetMarkerStyle(20);
	x_projection_20->SetMarkerColor(red_tran);
	x_projection_20->SetLineColor(red_tran);

	

	TF1* fit_x_40 = new TF1("fit_x_40", fit_breit_wigner_bg, fy_1, fy_2, 8);
	fit_x_40->SetParameters(750, 150, 1e14, 100, 0.046, 296,  1.68, -0.0036);
	//fit_x_40->SetParNames("m","Γ", "A_0", "C");
	fit_x_40->SetParLimits(0, 745, 770);
	//fit_x_40->SetParLimits(1, 100, 180);
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
	auto legend_x = new TLegend(0.70, 0.50, 0.95, 0.65); // bottom left, top right
	// legend_x->SetHeader("Invariant Mass Fits", "C");
	legend_x->AddEntry(x_projection_20, "20 Topology", "lep");
	legend_x->AddEntry(fit_x_20, "Gaussian Fit (20)", "l");
	legend_x->AddEntry(x_projection_40, "40 Topology", "lep");
	legend_x->AddEntry(fit_x_40, "Gaussian Fit (40)", "l");

	// Draw projections and fits
	x_projection_20->Draw("E1P");
	fit_x_20->Draw("SAME");
	x_projection_40->Draw("E1P SAME");
	fit_x_40->Draw("SAME");

	// Draw the legend
	legend_x->SetTextSize(0.035);
	legend_x->Draw();


	// manually set axis range if not drawn properly
	x_projection_20->SetAxisRange(50, 700, "Y");

	canvas_x->Update();
	TFile* outFile_x = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/kaons_from_rho_x_proj_all_20_and_40_bg.root", "RECREATE");
	canvas_x->Write();
	outFile_x->Close();

	canvas_x->SaveAs("media/photos/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/rho_reconstruction/kaons_from_rho_x_proj_all_20_and_40.pdf");
	canvas_x->Clear();


	return 0;
}



