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


// !./bashing/run_file.sh phi_phi_reconstruction/protons/pt_track_proton_comparison

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
		TH2F* hist_x = new TH2F("hist_x", "", bins, -1.6, 1.6, bins, -1.6, 1.6);
		TH2F* hist_y = new TH2F("hist_y", "", bins, -0.6, 0.6, bins, -0.6, 0.6);
		
		cx->SetTopMargin(0.02);
		cy->SetTopMargin(0.02);
		hist_x->GetXaxis()->SetTitleOffset(1.2);
		hist_y->GetXaxis()->SetTitleOffset(1.2);
		
		gStyle->SetOptStat(101110); // include overflow bin counter to verify info not lost
		
		if (topo == "40") {
			delete hist_y;
			TH2F* hist_y = new TH2F("hist_y", "", bins, -1.5, 1.5, bins, -1.5, 1.5);
		}

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
			[&hist_x, &hist_y] (const Double_t px_a, const Double_t py_a, const Double_t px_b, const Double_t py_b, const RVecF p, const RVecF phi) {
				Double_t px = px_a + px_b; 	// proton momentua 
				Double_t py = py_a + py_b;
				
				Float_t pk_x = 0;			// 4-track momenta
				Float_t pk_y = 0;
				for (Int_t i = 0; i < 4; ++i) {
					pk_x += p[i]*TMath::Cos(phi[i]); 
					pk_y += p[i]*TMath::Sin(phi[i]);
				}
				
				hist_x->Fill(px, pk_x);
				hist_y->Fill(py, pk_y);

			}, {"pr_px_a", "pr_py_a", "pr_px_b", "pr_py_b", "trk_pt", "trk_phi"}
		);


		RooRealVar x("x","x observable", hist_x->GetXaxis()->GetXmin(), hist_x->GetXaxis()->GetXmax());
		RooRealVar y("y","y observable", hist_x->GetYaxis()->GetXmin(), hist_x->GetYaxis()->GetXmax());

		RooRealVar mux("mux", "mean of x", 0, -0.1, 0.1);
		RooRealVar sigx("sigx", "sigma of x", 0.4, 0.1, 0.7);
		RooRealVar muy("muy", "mean of y", 0, -0.1, 0.1);
		RooRealVar sigy("sigy", "sigma of y", 0.4, 0.1, 0.7);
		RooRealVar rho("rho", "correlation", 0, -0.99, 0.99);
		RooGenericPdf* biv = BivariateGaussianPdf("biv", "Bivariate Gaussian", x, mux, sigx, y, muy, sigy, rho);

		Double_t integral = hist_x->Integral("width");
		RooDataHist dataHist("dataHist","dataset from TH2", RooArgList(x,y), hist_x);
		biv->fitTo(dataHist);


		TH2* modelHist = dynamic_cast<TH2*>(hist_x->Clone("modelHist"));
		modelHist->Reset();
		biv->fillHistogram(modelHist, RooArgList(x,y));
		modelHist->SetLineColor(kRed);  modelHist->Scale(integral, "width");
		calc_chi2(x, dataHist, biv, bins * bins);


		hist_x->SetTitle(TString::Format(";p_{x}^{4-track} (GeV/c) Events [%.3g];p_{x}^{proton} (GeV/c) Events [%.3g]", hist_x->GetXaxis()->GetBinWidth(1), hist_x->GetYaxis()->GetBinWidth(1)));
		hist_y->SetTitle(TString::Format(";p_{y}^{4-track} (GeV/c) Events [%.3g];p_{y}^{proton} (GeV/c) Events [%.3g]", hist_y->GetXaxis()->GetBinWidth(1), hist_y->GetYaxis()->GetBinWidth(1)));

		RVecDraw dat1 = {{hist_x, "COLZ"}
		//	,{modelHist, "CONT3"}
			};
		RVecDraw dat2 = {{hist_y, "COLZ"}};
		
		SaveCanvas(cx, dat1, TString::Format("media/root_files/phi_phi_reconstruction/protons/pt_track_proton_comparison/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "RECREATE");
		SaveCanvas(cy, dat2, TString::Format("media/root_files/phi_phi_reconstruction/protons/pt_track_proton_comparison/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "UPDATE");
		
		delete cx;
		delete hist_x;
		delete cy;
		delete hist_y;
	}
	}
	return 0;
}


