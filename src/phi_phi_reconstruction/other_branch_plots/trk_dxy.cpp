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


// !./bashing/run_file.sh phi_phi_reconstruction/other_branch_plots/trk_dxy

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
		TH1F* hist_x = new TH1F("hist_x", "", bins, -50, 50);
		
		cx->SetTopMargin(0.02);
		cx->SetRightMargin(0.05);
		hist_x->GetXaxis()->SetTitleOffset(1.2);
		
		gStyle->SetOptStat(111110); 

		df_df.Foreach(
				[&hist_x] (const RVecF dxy) {
				for (const Float_t d : dxy) {
				hist_x->Fill(d*1e2);
				}
			}, {"trk_dxy"});
		
		RooRealVar x("x","x observable", hist_x->GetXaxis()->GetXmin(), hist_x->GetXaxis()->GetXmax());
		RooRealVar mux("mux", "mean of x", 0, 5, 5);
		RooRealVar sigx("sigx", "sigma of x", 2, 0.5, 10);
		RooGaussian* gauss_x = new RooGaussian("gauss_x", "gauss_x", x, mux, sigx);

		Double_t integral_x = hist_x->Integral("width");
		
		RooDataHist dataHist_x("dataHist x","gauss x", RooArgList(x), hist_x);
		gauss_x->fitTo(dataHist_x);


		TH1* modelHist_x = (TH1*) hist_x->Clone("modelHist_y");
		modelHist_x->Reset();
		gauss_x->fillHistogram(modelHist_x, RooArgList(x));
		modelHist_x->SetLineColor(kRed);  modelHist_x->Scale(integral_x, "width");
		calc_chi2(x, dataHist_x, gauss_x, bins);

		hist_x->SetTitle(TString::Format("; #Delta#theta^{}_{x} (mrad); Events [%.3g]", hist_x->GetXaxis()->GetBinWidth(1)));

		for (int i=1; i<=modelHist_x->GetNbinsX(); i++) {
        modelHist_x->SetBinError(i, 0);
    	}
		RVecDraw dat1 = {
			{hist_x, "EP"},
			{modelHist_x, "C"}
			};
		SaveCanvas(cx, dat1, TString::Format("media/root_files/phi_phi_reconstruction/other_branch_plots/dxy/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "RECREATE");
		
		delete cx;
		delete hist_x;
		delete gauss_x;
		delete modelHist_x;
	}
	}
	return 0;
}



