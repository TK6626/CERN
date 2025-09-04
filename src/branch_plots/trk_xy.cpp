#include <stdio.h>
#include <cmath>
#include <vector>


#include "TPad.h"
#include "TMath.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TColor.h"
#include "TLegend.h"
#include "TFile.h"   
#include "ROOT/RVec.hxx"
#include "TLine.h"

#include "../../lib/computations.h"
#include "../../lib/plotting_params.h"



//	Command to run file
//	!g++ src/branch_plots/trk_xy.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/branch_plots/trk_xy; ./bin/branch_plots/trk_xy

/**
 * Look at cuts on the xy interaction point?
 * as secondary particle detectoins such as 
 * kaon reconstructions will be created further fron 
 * the interaction point

 */


void xy_xyerr_plot(RN df, TH1F* h);

int main(){

	SetPlotStyle();
	Int_t red_tran = TColor::GetFreeColorIndex();
	new TColor(red_tran, 1.0, 0.0, 0.0, "", 0.1); 
	Int_t blue_tran = TColor::GetFreeColorIndex();
	new TColor(blue_tran, 0.0, 0.0, 1.0, "", 0.1); 
	// load in 20 and 40 geometries
	
	const char* file_20 = "data/glueball_mass_reconstruction/data_20.root";
	const char* file_40 = "data/glueball_mass_reconstruction/data_40.root";



	ROOT::EnableImplicitMT();
	RDF df_20("tree1", file_20);
	RDF df_40("tree1", file_40);


	TCanvas* c = new TCanvas("c");
	Float_t c_x1 = -10;
	Float_t c_x2 = 10;
	TH1F* h_20 = new TH1F("h_20", "dxy / dxy_{err} with 3#sigma bounds", 1000, c_x1, c_x2);
	
	h_20->SetMarkerStyle(20);
	h_20->SetMarkerColor(red_tran);
	h_20->SetLineColor(red_tran);
	TH1F* h_40 = (TH1F*)h_20->Clone("h_40");
	h_40->SetMarkerColor(blue_tran);
	h_40->SetLineColor(blue_tran);
	
	xy_xyerr_plot(df_20, h_20);
	xy_xyerr_plot(df_40, h_40);

	h_20->GetXaxis()->SetTitle("#frac{xy}{xy_{err}}");
	h_20->GetYaxis()->SetTitle(Form("Events / bin width = %.2f", (double)h_20->GetXaxis()->GetBinWidth(1)));
	h_20->Draw();
	c->Update();

	Float_t x1_fit = -5; 	// no units 
	Float_t x2_fit = 5;		// no units

	// Create a TF1 for h_20
	TF1* fit_20 = new TF1("fit_20", fit_exp_quartic, x1_fit, x2_fit, 6);
	fit_20->SetParameters(0, 0, 0.003, 0,  1e5, 400);
	fit_20->SetParNames("a4", "a3", "a2", "a1", "A", "C");
	fit_20->SetParLimits(1, -1, 1);
	fit_20->SetLineColor(kRed);
	fit_20->SetLineWidth(2);
	fit_20->SetLineStyle(1);

	// Create a separate TF1 for h_40
	TF1* fit_40 = new TF1("fit_40", fit_exp_quartic, x1_fit, x2_fit, 6);
	fit_40->SetParameters(0, 0, 0.003, 0,  1e5, 400);
	fit_40->SetParNames("a4", "a3", "a2", "a1", "A", "C");
	fit_40->SetParLimits(1, -1, 1);
	fit_40->SetLineColor(kBlue + 1);
	fit_40->SetLineWidth(2);
	fit_40->SetLineStyle(1);

	// Fit histograms independently
	TFitResultPtr r_20 = h_20->Fit(fit_20, "S0");
	TFitResultPtr r_40 = h_40->Fit(fit_40, "S0");


	// Draw histograms and fits
	h_20->Draw("E1P");
	fit_20->Draw("SAME");
	h_40->Draw("E1P SAME");
	fit_40->Draw("SAME");


	// Get Lines for Illustrating where cuts are taken place
	Float_t sig = 3;
	Float_t x1 = r_20->Parameter(0) - sig * r_20->Parameter(1);
	Float_t x2 = r_20->Parameter(0) + sig * r_20->Parameter(1);

	Float_t y1 = h_20->GetMinimum();
	Float_t y2 = h_20->GetMaximum()* 1.05;

	TLine *upper_line = new TLine(x1, y1, x1, y2);
	//upper_line->SetNDC(true);
	upper_line->SetLineColor(kBlack);
  	upper_line->SetLineWidth(4);
	upper_line->SetLineStyle(2);
	TLine *lower_line = (TLine*)upper_line->Clone();


	lower_line->SetX1(x2); lower_line->SetX2(x2);
	upper_line->Draw("SAME");
	lower_line->Draw("SAME");
	
	auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
	legend->AddEntry(h_20, "Geometry 20", "l");
	legend->AddEntry(h_40, "Geometry 40", "l");
	legend->Draw("SAME");
	h_20->SetAxisRange(0, 35e3, "Y");
	c->Update();

	TFile* f = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/trk_xy_cut.root", "RECREATE");
	c->Write();
	f->Close();
	return 0;
}


void xy_xyerr_plot(RN df, TH1F* h){
	
	df.Foreach(
		[&h]
		(RVecF xy, RVecF xy_err, Int_t trk_num)
		{
			for (Int_t i = 0; i < trk_num; ++i){
				h->Fill(xy[i] / xy_err[i]);
			};
		},
		{"trk_dxy", "trk_dxyerr", "ntrk"}
	);

}

