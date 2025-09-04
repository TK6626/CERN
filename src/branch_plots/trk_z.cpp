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
#include "TLegend.h"
#include "TFile.h"   
#include "ROOT/RVec.hxx"
#include "TLine.h"
#include "../../lib/computations.h"
using RN = ROOT::RDF::RNode;

//	Command to run file
//	!g++ src/reconstruct_all_data/trk_z.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/trk_z; ./bin/reconstruct_all_data/trk_z

/**
 * Look at cuts on the z interaction point?
 * as secondary particle detectoins such as 
 * kaon reconstructions will be created further fron 
 * the interaction point
 */


void z_zerr_plot(RN df, TH1F* h);

int main(){


	// load in 20 and 40 geometries
	
	const char* file_20 = "data/glueball_mass_reconstruction/data_20.root";
	const char* file_40 = "data/glueball_mass_reconstruction/data_40.root";



	ROOT::EnableImplicitMT();
	RDF df_20("tree1", file_20);
	RDF df_40("tree1", file_40);


	TCanvas* c = new TCanvas("c");
	Float_t c_x1 = -10;
	Float_t c_x2 = 10;
	TH1F* h_20 = new TH1F("h_20", "dz / dz_{err} with 3#sigma bounds", 400, c_x1, c_x2);
	
	h_20->SetLineColor(kRed);
	h_20->SetLineWidth(2);
	TH1F* h_40 = (TH1F*)h_20->Clone("h_40");
	h_40->SetLineColor(kBlue);
	
	z_zerr_plot(df_20, h_20);
	z_zerr_plot(df_40, h_40);

	h_20->GetXaxis()->SetTitle("#frac{z}{z_{err}}");
	h_20->GetYaxis()->SetTitle(Form("Events / bin width = %.2f", (double)h_20->GetXaxis()->GetBinWidth(1)));
	h_20->Draw();
	c->Update();

	Float_t x1_fit = -3; 	// no units 
	Float_t x2_fit = 3;		// no units

	// Create a TF1 for h_20
	TF1* fit_20 = new TF1("fit_20", fit_gaussian, x1_fit, x2_fit, 4);
	fit_20->SetParameters(0, 1, 65e3, 250);
	fit_20->SetParNames("Mean", "Sigma", "Normal", "C");
	fit_20->SetParLimits(1, -1, 1);
	fit_20->SetLineColor(kRed);
	fit_20->SetLineWidth(2);
	fit_20->SetLineStyle(1);

	// Create a separate TF1 for h_40
	TF1* fit_40 = new TF1("fit_40", fit_gaussian, x1_fit, x2_fit, 4);
	fit_40->SetParameters(0, 1, 65e3, 250);
	fit_40->SetParNames("Mean", "Sigma", "Normal", "C");
	fit_40->SetParLimits(1, -1, 1);
	fit_40->SetLineColor(kBlue + 1);
	fit_40->SetLineWidth(2);
	fit_40->SetLineStyle(1);

	// Fit histograms independently
	TFitResultPtr r_20 = h_20->Fit(fit_20, "S0");
	TFitResultPtr r_40 = h_40->Fit(fit_40, "S0");

	// Draw histograms and fits
	h_20->Draw();
	fit_20->Draw("SAME");
	h_40->Draw("SAME");
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
	c->Update();

	TFile* f = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/glueball_mass_reconstruction/trk_z_cut.root", "RECREATE");
	c->Write();
	f->Close();
	return 0;
}


void z_zerr_plot(RN df, TH1F* h){
	
	df.Foreach(
		[&h]
		(RVecF z, RVecF z_err, Int_t trk_num)
		{
			for (Int_t i = 0; i < trk_num; ++i){
				h->Fill(z[i] / z_err[i]);
			};
		},
		{"trk_dz", "trk_dzerr", "ntrk"}
	);	
}

