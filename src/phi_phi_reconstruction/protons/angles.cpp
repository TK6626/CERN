#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"

#include "../../../lib/computations.h"
#include "../../../lib/admin_utils.h"
#include "../../../lib/plotting_params.h"
#include "../../../lib/roofithelper.h"


// !./bashing/run_file.sh phi_phi_reconstruction/protons/angles

int main() {

	SetPlotStyle();
	// set up transparant colours for use later
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi / 1e3;
	RVecF mass = {50000};
	mass /= 1e3; 
	Int_t bins = 400;

	RVecStr topology = {"20", "40"};
	for (Float_t mass_bound : mass) {
	ROOT::EnableImplicitMT();
	
		
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		TCanvas* cx = new TCanvas("cx", "cx", 600, 700);
		TCanvas* cy = new TCanvas("cy", "cy", 600, 700);
		TH1F* hist_x = new TH1F("hist_x", "", bins, -0.35, 0.35);
		TH1F* hist_y = new TH1F("hist_y", "", bins, -0.3, 0.3);
		
		cx->SetTopMargin(0.02);
		cy->SetTopMargin(0.02);
		cx->SetRightMargin(0.05);
		cy->SetRightMargin(0.05);
		hist_x->GetXaxis()->SetTitleOffset(1.2);
		hist_y->GetXaxis()->SetTitleOffset(1.2);
		
		gStyle->SetOptStat(111110); 

		RN df = df_df.Filter([phi_mass, mass_bound] (const RVecLorCyl p) {
			return 	( 
					((TMath::Abs(p[0].M() - phi_mass) < mass_bound) 
					|| (TMath::Abs(p[1].M() - phi_mass) < mass_bound))
				||
					((TMath::Abs(p[2].M() - phi_mass) < mass_bound) 
					|| (TMath::Abs(p[3].M() - phi_mass) < mass_bound))
				);
			},{"phi_four_momentum"});
		
		df.Foreach(
			[&hist_x, &hist_y] (Float_t x1, Float_t y1, Float_t x2, Float_t y2) {
				
				Float_t dx = x2 - x1; 
				Float_t dy = y2 - y1;
		
				hist_x->Fill(dx*1e3);
				hist_y->Fill(dy*1e3);
			}, {"ThxL", "ThyL", "ThxR", "ThyR"}
		);


		RooRealVar x("x","x observable", hist_x->GetXaxis()->GetXmin(), hist_x->GetXaxis()->GetXmax());
		RooRealVar y("y","y observable", hist_y->GetYaxis()->GetXmin(), hist_y->GetYaxis()->GetXmax());

		RooRealVar mux("mux", "mean of x", 0, -5e-3, 5e-3);
		RooRealVar sigx("sigx", "sigma of x", 0.07, 1e-3, 0.5);
		
		RooRealVar muy("muy", "mean of y", 0, -5e-3, 5e-3);
		RooRealVar sigy("sigy", "sigma of y", 0.07, 1e-3, 0.5);
		
		RooGaussian* gauss_x = new RooGaussian("gauss_x", "gauss_x", x, mux, sigx);
		RooGaussian* gauss_y = new RooGaussian("gauss_y", "gauss_y", y, muy, sigy);

		Double_t integral_x = hist_x->Integral("width");
		Double_t integral_y = hist_y->Integral("width");
		
		RooDataHist dataHist_x("dataHist x","gauss x", RooArgList(x), hist_x);
		gauss_x->fitTo(dataHist_x);
		RooDataHist dataHist_y("dataHist y","gauss y", RooArgList(y), hist_y);
		gauss_y->fitTo(dataHist_y);


		TH1* modelHist_x = (TH1*) hist_x->Clone("modelHist_y");
		modelHist_x->Reset();
		gauss_x->fillHistogram(modelHist_x, RooArgList(x));
		modelHist_x->SetLineColor(kRed);  modelHist_x->Scale(integral_x, "width");
		calc_chi2(x, dataHist_x, gauss_x, bins);

		TH1* modelHist_y = (TH1*) hist_y->Clone("modelHist_y");
		modelHist_y->Reset();
		gauss_y->fillHistogram(modelHist_y, RooArgList(y));
		modelHist_y->SetLineColor(kRed);  modelHist_y->Scale(integral_y, "width");
		calc_chi2(y, dataHist_y, gauss_y, bins);


		hist_x->SetTitle(TString::Format("; #Delta#theta^{}_{x} (mrad); Events [%.3g]", hist_x->GetXaxis()->GetBinWidth(1)));
		hist_y->SetTitle(TString::Format("; #Delta#theta^{}_{y} (mrad); Events [%.3g]", hist_y->GetXaxis()->GetBinWidth(1)));

		for (int i=1; i<=modelHist_x->GetNbinsX(); i++) {
        modelHist_x->SetBinError(i, 0);
    	}
		for (int i=1; i<=modelHist_y->GetNbinsX(); i++) {
        modelHist_y->SetBinError(i, 0);
    	}

		RVecDraw dat1 = {
			{hist_x, "EP"},
			{modelHist_x, "C"}
			};
		RVecDraw dat2 = {
			{hist_y, "EP"},
			{modelHist_x, "C"}
			};
		
		SaveCanvas(cx, dat1, TString::Format("media/root_files/phi_phi_reconstruction/protons/angles/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "RECREATE");
		SaveCanvas(cy, dat2, TString::Format("media/root_files/phi_phi_reconstruction/protons/angles/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "UPDATE");
		
		delete cx;
		delete hist_x;
		delete cy;
		delete hist_y;
		delete gauss_x;
		delete modelHist_x;
		delete modelHist_y;
		delete gauss_y;
	}
	}
	return 0;
}


