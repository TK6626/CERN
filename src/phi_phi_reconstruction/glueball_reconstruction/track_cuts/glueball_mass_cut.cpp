#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TColor.h"
#include "TFitResult.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TLegend.h"
#include "TLine.h"


#include "../../../lib/cut_branch.h"
#include "../../../lib/computations.h"
#include "../../../lib/apply_cuts.cpp"
#include "../../../lib/plotting_params.h"
#include "../../../lib/admin_utils.h"

// !./bashing/run_file.sh phi_phi_reconstruction/glueball_reconstruction/glueball_mass_cut


/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */
int main() {

	
	// constants 
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi /1e3;
	RVecF mass= {25, 20, 15, 10, 7, 5};
	mass /= 1e3;
	
	for (const Float_t mass_bound : mass) {

	RVecStr topology = {"20", "40"};
	ROOT::EnableImplicitMT();
	for (TString topo : topology) {
		
		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
	
		TCanvas* c = new TCanvas("c");
		TH1F* h = new TH1F("h", "h", 300, 2000, 3000);
		TLegend* leg = new TLegend(0.6,0.6,  0.9,0.9, "Legend");

		RN df = df_df;


		df.Foreach(
			[&h, mass_bound, phi_mass](const RVecLorCyl p) {
			
			if (
				(TMath::Abs((p[0].M() - phi_mass)) < mass_bound)  
				&& 
				(TMath::Abs((p[1].M() - phi_mass)) < mass_bound)
			) {h->Fill((p[0] + p[1]).M() * 1e3);}
			
			if (
				(TMath::Abs((p[2].M() - phi_mass)) < mass_bound)  
				&& 
				(TMath::Abs((p[3].M() - phi_mass)) < mass_bound)
			) {h->Fill((p[2] + p[3]).M() * 1e3);}
			}, {"phi_four_momentum"}
		);
		h->SetTitle(TString::Format("Invariant X Mass;M (MeV/c^{2});Events [%.2g MeV]", h->GetXaxis()->GetBinWidth(1)));
	

	TF1* fit = new TF1("fit", fit_gaussian, 2203, 2238, 4);
			fit->SetParameters(2220, 1, 10, 10);
			fit->SetParNames("m","Γ", "A_0", "C");
			fit->SetParLimits(0, 2218, 2222);
			fit->SetParLimits(1, 0, 10);
			fit->SetParLimits(2, 0, 100);
			fit->SetParLimits(3, 0, 100);
			fit->SetLineColor(kRed);
			fit->SetLineWidth(4);
			fit->SetLineStyle(2);
			TFitResultPtr r = h->Fit("fit", "S0");
		
			h->SetMarkerStyle(20);
			h->SetMarkerColor(kBlack);
			h->SetLineColor(kBlack);

			// Add entries to the legend AFTER setting styles
			TLegend* l = new TLegend(0.65, 0.55, 0.90, 0.70); // bottom left, top right
			// l->SetHeader("Invariant Mass Fits", "C");
			l->AddEntry(h, "Data", "lep");
			l->AddEntry(fit, "Gaussian Fit", "l");
			l->SetTextSize(0.035);
			
		TLine* l1 = DrawLine(2220, 0, h, kVertical); 
		RVecDraw dat = {{h, "HIST"}, {fit, ""}, {l1, ""}, {l, ""}};
		SaveCanvas(c, dat, TString::Format("media/root_files/phi_phi_reconstruction/glueball_reconstruction/mass_glueball/"+ topo +"/mass_bound=%.3g_MeV.root", mass_bound*1e3), "RECREATE");
		
		delete h;
		delete c;
	}
	}
	return 0;
}





