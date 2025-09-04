#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFitResult.h"


#include "../../../lib/computations.h"
#include "../../../lib/admin_utils.h"
#include "../../../lib/plotting_params.h"


// !./bashing/run_file.sh phi_phi_reconstruction/phi_phi_plot/projections


int main() {

	SetPlotStyle();
	
	// set up transparant colours for use later
	
	Float_t kaon_mass = m_kaon_char;
	Float_t phi_mass = m_phi;
	Float_t mass_bound = 10; 
	Float_t lower_bound = phi_mass - mass_bound; 
	Float_t upper_bound = phi_mass + mass_bound;


	RVecStr projection_axis = {"x", "y"};
	RVecStr topology = {"40", "20"};
	ROOT::EnableImplicitMT();
	
		
	for (TString topo : topology) {
	
		TFile* file = new TFile(TString::Format("media/root_files/phi_phi_reconstruction/phi_phi_plot/dxy/" + topo + "/0.10_cut.root"), "READ");
		TH2F* phi_phi_h = (TH2F*)file->Get("h");
		for (TString ax : projection_axis) {
			
			Int_t red_tran  = TColor::GetColorTransparent(kRed, 0.1);
			Int_t blue_tran = TColor::GetColorTransparent(kBlue, 0.1);

			TCanvas* c = new TCanvas("c", "c", 600, 800);
			
			//look for a window only 20 MeV wide	
			TH1D *h = nullptr;
			if (ax == "y"){
				Int_t lower_x_bin = phi_phi_h->GetXaxis()->FindFixBin(lower_bound);
				Int_t upper_x_bin = phi_phi_h->GetXaxis()->FindFixBin(upper_bound);
				//h = phi_phi_h->ProjectionY("h_proj", lower_x_bin, upper_x_bin);
				h = phi_phi_h->ProjectionY("h_proj", lower_x_bin, upper_x_bin);
			}
			if (ax == "x") {
				Int_t lower_y_bin = phi_phi_h->GetYaxis()->FindFixBin(lower_bound);
				Int_t upper_y_bin = phi_phi_h->GetYaxis()->FindFixBin(upper_bound);
				h = phi_phi_h->ProjectionX("h_pro", lower_y_bin, upper_y_bin);
			}
			
			h->SetTitle(TString::Format(ax + " Projection of Invariant K^{+}K^{-} Mass;M (MeV/c^{2});Events [%.2g MeV/c^{2}]", h->GetBinWidth(1)));
		
			Double_t fy_1 = 1008, fy_2 = 1031; // what part of the projection is actully fitted to 

			TF1* fit = new TF1("fit", fit_gaussian, fy_1, fy_2, 4);
			fit->SetParameters(phi_mass, 3, 80, 100);
			fit->SetParNames("m","Γ", "A_0", "C");
			fit->SetParLimits(0, 1017, 1021);
			fit->SetParLimits(1, 4, 8);
			fit->SetParLimits(2, 55, 150);
			fit->SetParLimits(3, 45, 60);
			fit->SetLineColor(kRed);
			fit->SetLineWidth(4);
			fit->SetLineStyle(2);
			TFitResultPtr r = h->Fit("fit", "S0");
		
			h->SetMarkerStyle(20);
			h->SetMarkerColor(17);
			h->SetLineColor(17);

			// Add entries to the legend AFTER setting styles
			TLegend* l = new TLegend(0.65, 0.55, 0.90, 0.70); // bottom left, top right
			// l->SetHeader("Invariant Mass Fits", "C");
			l->AddEntry(h, "Data", "lep");
			l->AddEntry(fit, "Gaussian Fit", "l");
			l->SetTextSize(0.035);
			
			h->Draw("E1P");
			l->Draw("SAME");
			fit->Draw("SAME");

			RVecDraw dat = {{h, "E1P"}, {l, ""}, {fit, ""}};
			SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/"+ax+"_projection/"+ topo +"/mass_bound=%.4gMeV.root", mass_bound), "RECREATE");
		
			delete c;
			delete fit;
			delete l;
			delete h;
		}
		delete file;
	}

	return 0;
}
