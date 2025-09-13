#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"


#include "../../../lib/computations.h"
#include "../../../lib/admin_utils.h"
#include "../../../lib/plotting_params.h"


// !./bashing/run_file.sh phi_phi_reconstruction/protons/elastic_pt_sum

int main() {

	SetPlotStyle();
	// set up transparant colours for use later
	Float_t kaon_mass = m_kaon_char / 1e3;
	Float_t phi_mass = m_phi / 1e3;
	RVecF mass = {20, 50, 50000};
	mass /= 1e3; 

	RVecStr topology = {"20", "40"};
	for (Float_t mass_bound : mass) {
	ROOT::EnableImplicitMT();
	
		
	for (TString topo : topology) {

		TString file = TString::Format("data/phi_phi_reconstruction/uncut_SNR_"	+ topo + ".root");
		RDF df_df("tree", file);
		TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
		TCanvas* c2 = new TCanvas("c2", "c2", 600, 800);
		TCanvas* c3 = new TCanvas("c3", "c3", 600, 800);
		TH1F* hist_both = new TH1F("hist_both", "", 200, -0.3, 0.3);
		TH1F* hist_m = new TH1F("hist_m", "", 200, 2, 5);
		TH2F* hist_comp = new TH2F("hist_comp", "", 400, -0.02, 2, 400, -0.02,  2);
		gStyle->SetOptStat(101111); // include overflow bin counter to verify info not lost

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
			[&hist_both, &hist_comp] (const Double_t px_a, const Double_t py_a, const Double_t px_b, const Double_t py_b, const RVecLorCyl p) {	
				Double_t px = px_a + px_b;
				Double_t py = py_a + py_b;
				Double_t pt = TMath::Sqrt(px*px + py*py);
			
				LorCyl p_k = p[0] + p[1] + p[2] + p[3];
				hist_both->Fill(pt - p_k.Pt());
				hist_comp->Fill(p_k.Pt(), pt);
			}, {"pr_px_a", "pr_py_a", "pr_px_b", "pr_py_b", "kaon_four_momentum"}
		);


		df.Foreach(
			[&hist_m] (const Double_t px_a, const Double_t py_a, const Double_t px_b, const Double_t py_b, const RVecLorCyl pk, const RVecLorCyl pp) {	
				Double_t px = px_a + px_b;
				Double_t py = py_a + py_b;
				Double_t pt2 = px*px + py*py;
				Float_t pk_pt = (pk[0] + pk[1] + pk[2] + pk[3]).Pt();
		
				if ((pk_pt*pk_pt /(0.15 * 0.15)) + (pt2 / 0.1) < 1){
					hist_m->Fill((pp[0] + pp[1]).M());
					hist_m->Fill((pp[2] + pp[3]).M());
					}
			}, {"pr_px_a", "pr_py_a", "pr_px_b", "pr_py_b", "kaon_four_momentum", "phi_four_momentum"}
		);



		hist_both->SetTitle(TString::Format("Proton - Kaon p_{t}; p_{t} (GeV/c); Events [%.3g]", hist_both->GetXaxis()->GetBinWidth(1)));
		hist_comp->SetTitle(TString::Format("Proton Kaon p_{t}Comparison; Kaon p_{t} (GeV/c);Proton p_{t} (GeV/c)"));
		hist_both->SetTitle(TString::Format("X mass; X (GeV/c); Events [%.3g]", hist_both->GetXaxis()->GetBinWidth(1)));

		RVecDraw dat1 = {{hist_both, "HIST"}};
		RVecDraw dat2 = {{hist_comp, "HIST"}};
		RVecDraw dat3 = {{hist_m, "HIST"}};
		SaveCanvas(c1, dat1, TString::Format("media/root_files/phi_phi_reconstruction/protons/elastic_pt_sum/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "RECREATE");
		SaveCanvas(c2, dat2, TString::Format("media/root_files/phi_phi_reconstruction/protons/elastic_pt_sum/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "UPDATE");
		SaveCanvas(c3, dat3, TString::Format("media/root_files/phi_phi_reconstruction/protons/elastic_pt_sum/"+topo+"/mass_cut=%.4gMeV.root", mass_bound*1e3), "UPDATE");
		
		delete c1;
		delete hist_both;
		delete c2;
		delete hist_comp;
		delete c3;
		delete hist_m;
	}
	}
	return 0;
}

